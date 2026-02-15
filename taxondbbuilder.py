import io
import logging
import os
import re
import time
import csv
import hashlib
from dataclasses import dataclass
from datetime import date, datetime
from enum import Enum
from pathlib import Path
from queue import Queue
from string import Formatter
from threading import Event, Lock, Thread
from typing import Dict, Iterable, List, Optional, Tuple
from urllib.error import HTTPError

import typer
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from rich.console import Console
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    import tomli as tomllib

app = typer.Typer(
    add_completion=False,
    help="TaxonDBBuilder - build a generic NCBI FASTA database by taxon and marker.",
)
console = Console()

DEFAULT_FEATURE_TYPES = ["rRNA", "gene", "CDS"]
DEFAULT_FEATURE_FIELDS = ["gene", "product", "note", "standard_name"]
DEFAULT_HEADER_FORMAT = "{acc_id}|{organism}|{marker}|{label}|{type}|{loc}|{strand}"
FORMATTER = Formatter()
IUPAC_DNA_VALUES = {k.upper(): v.upper() for k, v in ambiguous_dna_values.items()}
IUPAC_DNA_VALUES["U"] = "T"


class PostPrepStep(str, Enum):
    PRIMER_TRIM = "primer_trim"
    LENGTH_FILTER = "length_filter"
    DUPLICATE_REPORT = "duplicate_report"


POST_PREP_STEP_ORDER = [
    PostPrepStep.PRIMER_TRIM.value,
    PostPrepStep.LENGTH_FILTER.value,
    PostPrepStep.DUPLICATE_REPORT.value,
]


def print_header() -> None:
    console.print(
        Panel(
            "Generic NCBI GenBank fetch + feature extraction",
            title="TaxonDBBuilder",
            subtitle="DB FASTA builder",
            expand=False,
        )
    )


def render_run_table(
    config: Path,
    taxids: List[str],
    markers: List[str],
    out_path: Path,
    filters_cfg: Dict,
    dump_gb: Optional[Path] = None,
    from_gb: Optional[Path] = None,
    resume: bool = False,
) -> None:
    table = Table(title="Run Summary", show_header=True, header_style="bold")
    table.add_column("Item")
    table.add_column("Value", overflow="fold")
    table.add_row("Config", str(config))
    table.add_row("Taxids", ", ".join(taxids))
    table.add_row("Markers", ", ".join(markers))
    table.add_row("Output", str(out_path))
    if filters_cfg:
        table.add_row("Filters", ", ".join(sorted(filters_cfg.keys())))
    else:
        table.add_row("Filters", "none")
    if from_gb:
        table.add_row("From GB", str(from_gb))
    if dump_gb:
        table.add_row("Dump GB", str(dump_gb))
    if resume:
        table.add_row("Resume", "true")
    console.print(table)


def render_result_table(
    total_records: int,
    matched_records: int,
    matched_features: int,
    kept_records: int,
    skipped_same: int,
    duplicated_diff: int,
    out_path: Path,
    log_path: Path,
) -> None:
    table = Table(title="Result Summary", show_header=True, header_style="bold")
    table.add_column("Metric")
    table.add_column("Value")
    table.add_row("Total records", str(total_records))
    table.add_row("Matched records", str(matched_records))
    table.add_row("Matched features", str(matched_features))
    table.add_row("Kept records", str(kept_records))
    table.add_row("Skipped duplicates (same acc+seq)", str(skipped_same))
    table.add_row("Kept duplicates (same acc, diff seq)", str(duplicated_diff))
    table.add_row("Output", str(out_path))
    table.add_row("Log", str(log_path))
    console.print(table)


def process_genbank_chunk(
    chunk: str,
    marker_rules: List[Dict],
    acc_to_seqs: Dict[str, set],
    out_f,
    counters: Dict[str, int],
    dup_accessions: Dict[str, int],
    lock: Lock,
    progress: Progress,
    task_id: int,
    taxid: str,
    dump_gb_dir: Optional[Path],
    emitted_records: List[Dict[str, str]],
) -> None:
    handle = io.StringIO(chunk)
    for record in SeqIO.parse(handle, "genbank"):
        with lock:
            counters["total_records"] += 1
        progress.update(task_id, advance=1)

        if not record.seq:
            continue

        acc = record.id
        organism = record.annotations.get("organism", "unknown")
        organism_safe = sanitize_header(str(organism))
        record_matched = False

        if dump_gb_dir:
            gb_root = dump_gb_dir / f"taxid{taxid}"
            gb_root.mkdir(parents=True, exist_ok=True)
            acc_safe = sanitize_header(acc)
            gb_path = gb_root / f"{acc_safe}.gb"
            with lock:
                if not gb_path.exists():
                    with gb_path.open("w", encoding="utf-8") as f:
                        SeqIO.write(record, f, "genbank")

        for feature in record.features:
            matched_label = None
            matched_marker = None
            matched_type = feature.type

            header_format = None
            for rule in marker_rules:
                if rule["feature_types"] and feature.type not in rule["feature_types"]:
                    continue
                matched_label = match_feature(
                    feature,
                    rule["patterns"],
                    rule["feature_fields"],
                )
                if matched_label:
                    matched_marker = rule["key"]
                    header_format = rule["header_format"]
                    break

            if not matched_label or not feature.location:
                continue

            try:
                seq = str(feature.extract(record.seq)).upper()
            except Exception:
                continue

            if not seq:
                continue

            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            strand = feature.location.strand or 0

            label_safe = sanitize_header(str(matched_label))
            marker_safe = sanitize_header(str(matched_marker or "marker"))
            type_safe = sanitize_header(str(matched_type))
            loc = f"{start}-{end}"

            with lock:
                counters["matched_features"] += 1
                record_matched = True

                seqs = acc_to_seqs.setdefault(acc, set())
                if seq in seqs:
                    counters["skipped_same"] += 1
                    continue

                dup_index = None
                if seqs:
                    dup_index = len(seqs) + 1
                    counters["duplicated_diff"] += 1
                    dup_accessions[acc] = dup_accessions.get(acc, 1) + 1

                seqs.add(seq)

                acc_id = f"{acc}_dup{dup_index}" if dup_index else acc
                dup_tag = f"dup{dup_index}" if dup_index else ""
                header_values = {
                    "acc": acc,
                    "acc_id": acc_id,
                    "organism": organism_safe,
                    "organism_raw": str(organism),
                    "marker": marker_safe,
                    "marker_raw": str(matched_marker or ""),
                    "label": label_safe,
                    "label_raw": str(matched_label),
                    "type": type_safe,
                    "type_raw": str(matched_type),
                    "start": str(start),
                    "end": str(end),
                    "loc": loc,
                    "strand": str(strand),
                    "dup": dup_tag,
                }
                header = build_header(header_format or DEFAULT_HEADER_FORMAT, header_values).strip()
                if not header:
                    header = acc_id

                out_f.write(f">{header}\n")
                out_f.write(f"{seq}\n")
                counters["kept_records"] += 1
                emitted_records.append(
                    {
                        "acc_id": acc_id,
                        "accession": acc,
                        "organism_name": str(organism),
                        "header": header,
                    }
                )
        if record_matched:
            with lock:
                counters["matched_records"] += 1



def load_config(path: Path) -> Dict:
    if not path.exists():
        raise typer.BadParameter(f"Config file not found: {path}")
    with path.open("rb") as f:
        data = tomllib.load(f)

    if "markers_file" in data:
        raise typer.BadParameter("markers_file must be defined under [markers].file (top-level is not supported).")

    markers_file = None
    markers_section = data.get("markers")
    if isinstance(markers_section, dict):
        if "markers_file" in markers_section:
            raise typer.BadParameter("Use [markers].file instead of [markers].markers_file.")
        if "file" in markers_section:
            markers_file = markers_section.get("file")
            if not isinstance(markers_file, str):
                raise typer.BadParameter("markers.file must be a string path.")

    inline_markers: Dict[str, Dict] = {}
    if markers_section is None:
        markers_section = {}
    if not isinstance(markers_section, dict):
        raise typer.BadParameter("[markers] must be a table (dict).")
    for key, value in markers_section.items():
        if key in ("file", "markers_file"):
            continue
        if not isinstance(value, dict):
            raise typer.BadParameter(f"markers.{key} must be a table (dict).")
        inline_markers[key] = value

    markers_from_file: Dict[str, Dict] = {}
    if markers_file:
        if not isinstance(markers_file, str):
            raise typer.BadParameter("markers.file must be a string path.")
        markers_path = Path(os.path.expandvars(os.path.expanduser(markers_file)))
        candidates: List[Path] = []
        if markers_path.is_absolute():
            candidates.append(markers_path)
        else:
            candidates.append(path.parent / markers_path)
            candidates.append(Path.cwd() / markers_path)
            candidates.append(Path(__file__).resolve().parent / markers_path)

        markers_path = next((p for p in candidates if p.exists()), None)
        if not markers_path:
            tried = ", ".join(str(p) for p in candidates)
            raise typer.BadParameter(f"Markers file not found. Tried: {tried}")
        with markers_path.open("rb") as f:
            markers_data = tomllib.load(f)
        markers_from_file = markers_data.get("markers")
        if not isinstance(markers_from_file, dict) or not markers_from_file:
            raise typer.BadParameter("Markers file must define a non-empty [markers] section.")

    merged = {}
    merged.update(markers_from_file)
    merged.update(inline_markers)
    if merged:
        data["markers"] = merged

    if "ncbi" not in data:
        raise typer.BadParameter("Missing [ncbi] section in config.")
    if "markers" not in data or not data["markers"]:
        raise typer.BadParameter("Missing [markers] section in config.")

    post_prep = data.get("post_prep")
    if post_prep is not None:
        if not isinstance(post_prep, dict):
            raise typer.BadParameter("[post_prep] must be a table (dict).")

        min_len = post_prep.get("sequence_length_min")
        max_len = post_prep.get("sequence_length_max")
        if min_len is not None:
            try:
                min_len = int(min_len)
                post_prep["sequence_length_min"] = min_len
            except (TypeError, ValueError) as exc:
                raise typer.BadParameter("post_prep.sequence_length_min must be an integer.") from exc
        if max_len is not None:
            try:
                max_len = int(max_len)
                post_prep["sequence_length_max"] = max_len
            except (TypeError, ValueError) as exc:
                raise typer.BadParameter("post_prep.sequence_length_max must be an integer.") from exc
        if min_len is not None and max_len is not None and min_len > max_len:
            raise typer.BadParameter("post_prep.sequence_length_min must be <= post_prep.sequence_length_max.")

        primer_file = post_prep.get("primer_file")
        primer_set_raw = post_prep.get("primer_set")
        primer_set_list: Optional[List[str]] = None
        if primer_file is not None:
            if not isinstance(primer_file, str) or not primer_file.strip():
                raise typer.BadParameter("post_prep.primer_file must be a non-empty string path.")
            primer_file = primer_file.strip()
            post_prep["primer_file"] = primer_file
        if primer_set_raw is not None:
            if isinstance(primer_set_raw, str):
                primer_set_list = [primer_set_raw]
            elif isinstance(primer_set_raw, list) and all(isinstance(v, str) for v in primer_set_raw):
                primer_set_list = list(primer_set_raw)
            else:
                raise typer.BadParameter("post_prep.primer_set must be a string or list of strings.")

            normalized_sets: List[str] = []
            for value in primer_set_list:
                name = value.strip()
                if not name:
                    raise typer.BadParameter("post_prep.primer_set cannot contain empty values.")
                if name not in normalized_sets:
                    normalized_sets.append(name)
            primer_set_list = normalized_sets
            post_prep["primer_set"] = primer_set_list

        if primer_set_list is not None and primer_file is None:
            raise typer.BadParameter("post_prep.primer_set requires post_prep.primer_file.")

        if primer_file and primer_set_list:
            primer_path = resolve_support_file_path(primer_file, path, "Primer file")
            primer_sets_data = load_primer_sets_from_file(primer_path)
            forward, reverse = combine_primer_set_sequences(primer_sets_data, primer_set_list)

            post_prep["_primer_forward"] = forward
            post_prep["_primer_reverse"] = reverse
            post_prep["_primer_file_resolved"] = str(primer_path)
            post_prep["_primer_set_names"] = primer_set_list
            post_prep["_primer_set_candidates"] = sorted(primer_sets_data.keys())
        elif primer_file:
            primer_path = resolve_support_file_path(primer_file, path, "Primer file")
            primer_sets_data = load_primer_sets_from_file(primer_path)
            post_prep["_primer_file_resolved"] = str(primer_path)
            post_prep["_primer_set_candidates"] = sorted(primer_sets_data.keys())

    return data


def setup_entrez(ncbi_cfg: Dict) -> None:
    email = ncbi_cfg.get("email") or os.environ.get("NCBI_EMAIL")
    api_key = ncbi_cfg.get("api_key") or os.environ.get("NCBI_API_KEY")

    if not email:
        console.print("[yellow]WARNING:[/yellow] NCBI email is not set. Set ncbi.email or NCBI_EMAIL.")
    Entrez.email = email or ""

    if api_key:
        Entrez.api_key = api_key


def normalize_marker_map(markers_cfg: Dict) -> Dict[str, Dict]:
    marker_map: Dict[str, Dict] = {}
    for key, cfg in markers_cfg.items():
        phrases = cfg.get("phrases") or []
        terms = cfg.get("terms") or []
        region_patterns = cfg.get("region_patterns") or []
        header_format = cfg.get("header_format")
        feature_types = cfg.get("feature_types")
        feature_fields = cfg.get("feature_fields")
        if not isinstance(phrases, list):
            raise typer.BadParameter(f"markers.{key}.phrases must be a list of strings.")
        if not isinstance(terms, list):
            raise typer.BadParameter(f"markers.{key}.terms must be a list of strings.")
        if not isinstance(region_patterns, list):
            raise typer.BadParameter(f"markers.{key}.region_patterns must be a list of strings.")
        if feature_types is not None and not isinstance(feature_types, list):
            raise typer.BadParameter(f"markers.{key}.feature_types must be a list of strings.")
        if feature_fields is not None and not isinstance(feature_fields, list):
            raise typer.BadParameter(f"markers.{key}.feature_fields must be a list of strings.")
        if not phrases and not terms:
            raise typer.BadParameter(f"markers.{key} must define phrases or terms.")
        aliases = cfg.get("aliases") or []
        marker_map[key] = {
            "phrases": phrases,
            "terms": terms,
            "aliases": aliases,
            "region_patterns": region_patterns,
            "header_format": header_format,
            "feature_types": feature_types,
            "feature_fields": feature_fields,
        }
    return marker_map


def resolve_marker_key(value: str, marker_map: Dict[str, Dict]) -> str:
    value_l = value.lower()
    exact_matches = []
    prefix_matches = []

    for key, cfg in marker_map.items():
        key_l = key.lower()
        aliases = [a.lower() for a in cfg.get("aliases", [])]
        if value_l == key_l or value_l in aliases:
            exact_matches.append(key)
            continue
        if key_l.startswith(value_l) or any(a.startswith(value_l) for a in aliases):
            prefix_matches.append(key)

    if exact_matches:
        if len(exact_matches) > 1:
            raise typer.BadParameter(f"Marker '{value}' matches multiple entries: {exact_matches}")
        return exact_matches[0]
    if prefix_matches:
        if len(prefix_matches) > 1:
            raise typer.BadParameter(f"Marker '{value}' matches multiple entries: {prefix_matches}")
        return prefix_matches[0]

    raise typer.BadParameter(f"Marker '{value}' not found in config.")


def is_raw_term(value: str) -> bool:
    return "[" in value and "]" in value


def build_marker_query(marker_keys: List[str], marker_map: Dict[str, Dict]) -> str:
    phrase_terms: List[str] = []
    for key in marker_keys:
        for term in marker_map[key].get("terms", []):
            phrase_terms.append(term)
        for phrase in marker_map[key].get("phrases", []):
            if is_raw_term(phrase):
                phrase_terms.append(phrase)
            else:
                phrase_escaped = phrase.replace('"', '\\"')
                phrase_terms.append(f'"{phrase_escaped}"[All Fields]')
    if not phrase_terms:
        raise typer.BadParameter("No marker phrases resolved.")
    if len(phrase_terms) == 1:
        return phrase_terms[0]
    return "(" + " OR ".join(f"({t})" for t in phrase_terms) + ")"


def strip_field_spec(value: str) -> str:
    stripped = re.sub(r"\s*\[[^\]]+\]", "", value)
    return stripped.replace('"', "").strip()


def build_region_patterns(marker_cfg: Dict) -> List[str]:
    patterns = marker_cfg.get("region_patterns") or []
    if patterns:
        return patterns

    raw = []
    raw.extend(marker_cfg.get("terms", []))
    raw.extend(marker_cfg.get("phrases", []))
    cleaned = []
    for item in raw:
        stripped = strip_field_spec(item)
        if stripped:
            cleaned.append(re.escape(stripped))
    return cleaned


def compile_patterns(patterns: List[str]) -> List[re.Pattern]:
    compiled = []
    for pat in patterns:
        if not pat:
            continue
        compiled.append(re.compile(pat, re.IGNORECASE))
    return compiled


def sanitize_header(text: str) -> str:
    safe = re.sub(r"\s+", "_", text.strip())
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", safe)
    return safe


def resolve_header_format(cfg: Dict, output_cfg: Dict) -> str:
    header_formats = output_cfg.get("header_formats") or {}
    if not isinstance(header_formats, dict):
        raise typer.BadParameter("[output].header_formats must be a table (dict).")
    default_format = output_cfg.get("default_header_format", DEFAULT_HEADER_FORMAT)
    if not isinstance(default_format, str):
        raise typer.BadParameter("[output].default_header_format must be a string.")
    header_key = cfg.get("header_format")
    if not header_key:
        return default_format
    if not isinstance(header_key, str):
        raise typer.BadParameter("markers.<id>.header_format must be a string.")
    if header_key in header_formats:
        value = header_formats[header_key]
        if not isinstance(value, str):
            raise typer.BadParameter(f"header_formats.{header_key} must be a string.")
        return value
    return header_key


class SafeFormatDict(dict):
    def __missing__(self, key: str) -> str:
        return ""


def build_header(template: str, values: Dict[str, str]) -> str:
    return template.format_map(SafeFormatDict(values))


def template_has_field(template: str, field_name: str) -> bool:
    for _, parsed_field, _, _ in FORMATTER.parse(template):
        if parsed_field == field_name:
            return True
    return False


@dataclass(frozen=True)
class HeaderExtractor:
    pattern: re.Pattern
    captures_organism: bool


def compile_header_extractors(header_formats: Iterable[str]) -> Tuple[List[HeaderExtractor], bool, bool]:
    extractors: List[HeaderExtractor] = []
    has_acc_id_template = False
    has_organism_template = False
    for template in sorted(set(header_formats)):
        has_acc_id = template_has_field(template, "acc_id")
        organism_field = None
        if template_has_field(template, "organism_raw"):
            organism_field = "organism_raw"
        elif template_has_field(template, "organism"):
            organism_field = "organism"

        has_acc_id_template = has_acc_id_template or has_acc_id
        has_organism_template = has_organism_template or organism_field is not None
        if not has_acc_id:
            continue

        parts: List[str] = []
        seen_acc_id = False
        seen_organism = False
        for literal, field, _, _ in FORMATTER.parse(template):
            parts.append(re.escape(literal))
            if field is None:
                continue
            if field == "acc_id":
                if not seen_acc_id:
                    parts.append(r"(?P<acc_id>[A-Za-z0-9._-]+)")
                    seen_acc_id = True
                else:
                    parts.append(r"(?P=acc_id)")
            elif organism_field and field == organism_field:
                if not seen_organism:
                    parts.append(r"(?P<organism_name>.*?)")
                    seen_organism = True
                else:
                    parts.append(r"(?P=organism_name)")
            else:
                parts.append(r".*?")
        if seen_acc_id:
            extractors.append(
                HeaderExtractor(
                    pattern=re.compile("^" + "".join(parts) + "$"),
                    captures_organism=seen_organism,
                )
            )
    return extractors, has_acc_id_template, has_organism_template


def extract_header_fields_from_header(header: str, extractors: List[HeaderExtractor]) -> Tuple[Optional[str], Optional[str]]:
    for extractor in extractors:
        match = extractor.pattern.match(header)
        if not match:
            continue
        acc_id = match.groupdict().get("acc_id")
        if not acc_id:
            continue
        organism_name = match.groupdict().get("organism_name")
        if organism_name is not None:
            organism_name = organism_name.strip()
        return acc_id, organism_name or None
    return None, None


def apply_post_prep_length_filter(
    fasta_path: Path,
    min_len: Optional[int],
    max_len: Optional[int],
) -> Dict[str, int]:
    before_count = 0
    after_count = 0
    tmp_path = fasta_path.with_suffix(fasta_path.suffix + ".postprep.tmp")
    try:
        with fasta_path.open("r", encoding="utf-8") as in_f, tmp_path.open("w", encoding="utf-8") as out_f:
            for record in SeqIO.parse(in_f, "fasta"):
                before_count += 1
                seq_len = len(record.seq)
                if min_len is not None and seq_len < min_len:
                    continue
                if max_len is not None and seq_len > max_len:
                    continue
                SeqIO.write(record, out_f, "fasta")
                after_count += 1
        tmp_path.replace(fasta_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()
    return {
        "before": before_count,
        "after": after_count,
        "removed": before_count - after_count,
    }


def primer_to_regex(primer: str) -> str:
    parts: List[str] = []
    for ch in primer.upper():
        values = IUPAC_DNA_VALUES.get(ch)
        if values is None:
            raise ValueError(f"Unsupported primer base: {ch}")
        uniq = "".join(sorted(set(values)))
        if len(uniq) == 1:
            parts.append(uniq)
        else:
            parts.append(f"[{uniq}]")
    return "".join(parts)


def compile_primer_end_patterns(
    primers: List[str],
    anchor: str,
) -> List[Tuple[int, re.Pattern]]:
    compiled: List[Tuple[int, re.Pattern]] = []
    for primer in primers:
        regex = primer_to_regex(primer)
        if anchor == "start":
            pat = re.compile("^" + regex, re.IGNORECASE)
        elif anchor == "end":
            pat = re.compile(regex + "$", re.IGNORECASE)
        else:
            raise ValueError(f"Unknown primer anchor: {anchor}")
        compiled.append((len(primer), pat))
    return compiled


def matched_prefix_length(seq: str, patterns: List[Tuple[int, re.Pattern]]) -> int:
    max_len = 0
    for primer_len, pat in patterns:
        if pat.match(seq):
            max_len = max(max_len, primer_len)
    return max_len


def matched_suffix_length(seq: str, patterns: List[Tuple[int, re.Pattern]]) -> int:
    max_len = 0
    for primer_len, pat in patterns:
        if pat.search(seq):
            max_len = max(max_len, primer_len)
    return max_len


def choose_trim_lengths(
    seq: str,
    left_patterns: List[Tuple[int, re.Pattern]],
    right_patterns: List[Tuple[int, re.Pattern]],
) -> Tuple[int, int]:
    left = matched_prefix_length(seq, left_patterns)
    right = matched_suffix_length(seq, right_patterns)
    return left, right


def apply_post_prep_primer_trim(
    fasta_path: Path,
    forward_primers: List[str],
    reverse_primers: List[str],
) -> Dict[str, int]:
    forward_rc = [str(Seq(p).reverse_complement()).upper() for p in forward_primers]
    reverse_rc = [str(Seq(p).reverse_complement()).upper() for p in reverse_primers]

    canonical_left = compile_primer_end_patterns(forward_primers, anchor="start")
    canonical_right = compile_primer_end_patterns(reverse_rc, anchor="end")
    reverse_left = compile_primer_end_patterns(reverse_primers, anchor="start")
    reverse_right = compile_primer_end_patterns(forward_rc, anchor="end")

    before_count = 0
    after_count = 0
    trimmed_both = 0
    trimmed_left_only = 0
    trimmed_right_only = 0
    untrimmed = 0
    dropped_empty = 0
    canonical_orientation = 0
    reverse_orientation = 0

    tmp_path = fasta_path.with_suffix(fasta_path.suffix + ".postprep.primer.tmp")
    try:
        with fasta_path.open("r", encoding="utf-8") as in_f, tmp_path.open("w", encoding="utf-8") as out_f:
            for record in SeqIO.parse(in_f, "fasta"):
                before_count += 1
                seq = str(record.seq).upper().replace("U", "T")
                seq_len = len(seq)

                can_left, can_right = choose_trim_lengths(seq, canonical_left, canonical_right)
                rev_left, rev_right = choose_trim_lengths(seq, reverse_left, reverse_right)

                can_ends = (1 if can_left else 0) + (1 if can_right else 0)
                rev_ends = (1 if rev_left else 0) + (1 if rev_right else 0)
                can_trim = can_left + can_right
                rev_trim = rev_left + rev_right

                if (can_ends, can_trim) >= (rev_ends, rev_trim):
                    left_len, right_len = can_left, can_right
                    chosen_orientation = "canonical"
                else:
                    left_len, right_len = rev_left, rev_right
                    chosen_orientation = "reverse"

                ends_trimmed = (1 if left_len else 0) + (1 if right_len else 0)
                if ends_trimmed == 0:
                    untrimmed += 1
                    trimmed_seq = seq
                else:
                    if ends_trimmed == 2:
                        trimmed_both += 1
                    elif left_len:
                        trimmed_left_only += 1
                    else:
                        trimmed_right_only += 1

                    if chosen_orientation == "canonical":
                        canonical_orientation += 1
                    else:
                        reverse_orientation += 1

                    right_index = seq_len - right_len if right_len else seq_len
                    trimmed_seq = seq[left_len:right_index]

                if not trimmed_seq:
                    dropped_empty += 1
                    continue

                out_f.write(f">{record.description}\n")
                out_f.write(f"{trimmed_seq}\n")
                after_count += 1

        tmp_path.replace(fasta_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()

    return {
        "before": before_count,
        "after": after_count,
        "removed": before_count - after_count,
        "trimmed_both": trimmed_both,
        "trimmed_left_only": trimmed_left_only,
        "trimmed_right_only": trimmed_right_only,
        "untrimmed": untrimmed,
        "dropped_empty": dropped_empty,
        "canonical_orientation": canonical_orientation,
        "reverse_orientation": reverse_orientation,
    }


def write_duplicate_acc_reports_csv(
    fasta_path: Path,
    header_formats: Iterable[str],
) -> Tuple[Optional[Path], Optional[Path], Dict[str, int], Optional[str]]:
    extractors, has_acc_id_template, has_organism_template = compile_header_extractors(header_formats)
    if not has_acc_id_template:
        return None, None, {}, "selected header format does not include {acc_id}"
    if not has_organism_template:
        return None, None, {}, "selected header format does not include {organism_raw} or {organism}"
    if not extractors:
        return None, None, {}, "could not compile header extractor"
    if not any(extractor.captures_organism for extractor in extractors):
        return None, None, {}, "selected header format does not include {organism_raw} or {organism} alongside {acc_id}"

    seq_groups: Dict[str, List[Dict[str, str]]] = {}
    total_records = 0
    parsed_records = 0
    unparsed_records = 0
    with fasta_path.open("r", encoding="utf-8") as fasta_f:
        for record in SeqIO.parse(fasta_f, "fasta"):
            total_records += 1
            header = str(record.description).strip()
            acc_id, organism_name = extract_header_fields_from_header(header, extractors)
            if not acc_id or not organism_name:
                unparsed_records += 1
                continue
            parsed_records += 1
            accession = re.sub(r"_dup\d+$", "", acc_id)
            seq = str(record.seq).upper()
            seq_groups.setdefault(seq, []).append(
                {
                    "acc_id": acc_id,
                    "accession": accession,
                    "organism_name": organism_name,
                    "header": header,
                }
            )

    records_path = fasta_path.with_suffix(fasta_path.suffix + ".duplicate_acc.records.csv")
    groups_path = fasta_path.with_suffix(fasta_path.suffix + ".duplicate_acc.groups.csv")
    duplicate_groups = 0
    duplicate_records = 0
    cross_organism_groups = 0
    with records_path.open("w", newline="", encoding="utf-8") as records_f, groups_path.open(
        "w", newline="", encoding="utf-8"
    ) as groups_f:
        records_writer = csv.DictWriter(
            records_f,
            fieldnames=[
                "group_id",
                "sequence_hash",
                "sequence_length",
                "records_in_group",
                "unique_accessions",
                "unique_organisms",
                "cross_organism_duplicate",
                "acc_id",
                "accession",
                "organism_name",
                "header",
            ],
        )
        records_writer.writeheader()
        groups_writer = csv.DictWriter(
            groups_f,
            fieldnames=[
                "group_id",
                "sequence_hash",
                "sequence_length",
                "records_in_group",
                "unique_accessions",
                "unique_organisms",
                "cross_organism_duplicate",
                "accessions",
                "organism_names",
            ],
        )
        groups_writer.writeheader()

        group_id = 0
        for seq, entries in sorted(seq_groups.items(), key=lambda item: len(item[1]), reverse=True):
            if len(entries) < 2:
                continue
            unique_accessions = sorted({item["accession"] for item in entries})
            if len(unique_accessions) < 2:
                continue
            unique_organisms = sorted({item["organism_name"] for item in entries})
            is_cross_organism = len(unique_organisms) > 1
            group_id += 1
            duplicate_groups += 1
            duplicate_records += len(entries)
            if is_cross_organism:
                cross_organism_groups += 1
            seq_hash = hashlib.sha1(seq.encode("utf-8")).hexdigest()
            groups_writer.writerow(
                {
                    "group_id": group_id,
                    "sequence_hash": seq_hash,
                    "sequence_length": len(seq),
                    "records_in_group": len(entries),
                    "unique_accessions": len(unique_accessions),
                    "unique_organisms": len(unique_organisms),
                    "cross_organism_duplicate": str(is_cross_organism).lower(),
                    "accessions": ";".join(unique_accessions),
                    "organism_names": ";".join(unique_organisms),
                }
            )
            for item in entries:
                records_writer.writerow(
                    {
                        "group_id": group_id,
                        "sequence_hash": seq_hash,
                        "sequence_length": len(seq),
                        "records_in_group": len(entries),
                        "unique_accessions": len(unique_accessions),
                        "unique_organisms": len(unique_organisms),
                        "cross_organism_duplicate": str(is_cross_organism).lower(),
                        "acc_id": item["acc_id"],
                        "accession": item["accession"],
                        "organism_name": item["organism_name"],
                        "header": item["header"],
                    }
                )

    stats = {
        "total_records": total_records,
        "parsed_records": parsed_records,
        "unparsed_records": unparsed_records,
        "duplicate_groups": duplicate_groups,
        "duplicate_records": duplicate_records,
        "cross_organism_groups": cross_organism_groups,
    }
    return records_path, groups_path, stats, None


def write_acc_organism_mapping_csv(
    fasta_path: Path,
    emitted_records: List[Dict[str, str]],
) -> Tuple[Path, Dict[str, int]]:
    mapping_path = fasta_path.with_suffix(fasta_path.suffix + ".acc_organism.csv")
    records_by_header: Dict[str, List[Dict[str, str]]] = {}
    for row in emitted_records:
        header = row.get("header", "")
        records_by_header.setdefault(header, []).append(row)

    offsets: Dict[str, int] = {}
    total_records = 0
    mapped_records = 0
    unmapped_records = 0
    unique_accessions = set()
    unique_organisms = set()
    mapped_rows: List[Dict[str, str]] = []
    with fasta_path.open("r", encoding="utf-8") as fasta_f:
        for record in SeqIO.parse(fasta_f, "fasta"):
            total_records += 1
            header = str(record.description).strip()
            rows = records_by_header.get(header)
            idx = offsets.get(header, 0)
            if not rows or idx >= len(rows):
                unmapped_records += 1
                continue
            row = rows[idx]
            offsets[header] = idx + 1
            mapped_records += 1
            unique_accessions.add(row["accession"])
            unique_organisms.add(row["organism_name"])
            mapped_rows.append(
                {
                    "acc_id": row["acc_id"],
                    "accession": row["accession"],
                    "organism_name": row["organism_name"],
                    "header": header,
                }
            )

    mapped_rows.sort(key=lambda row: (row["acc_id"], row["accession"], row["organism_name"], row["header"]))

    with mapping_path.open("w", newline="", encoding="utf-8") as out_f:
        writer = csv.DictWriter(
            out_f,
            fieldnames=["acc_id", "accession", "organism_name", "header"],
        )
        writer.writeheader()
        writer.writerows(mapped_rows)

    unused_records = 0
    for header, rows in records_by_header.items():
        unused_records += max(0, len(rows) - offsets.get(header, 0))

    stats = {
        "total_records": total_records,
        "mapped_records": mapped_records,
        "unmapped_records": unmapped_records,
        "unused_records": unused_records,
        "unique_accessions": len(unique_accessions),
        "unique_organisms": len(unique_organisms),
    }
    return mapping_path, stats


def feature_texts(feature, fields: List[str]) -> List[str]:
    texts: List[str] = []
    for field in fields:
        val = feature.qualifiers.get(field)
        if not val:
            continue
        if isinstance(val, list):
            texts.extend([str(v) for v in val])
        else:
            texts.append(str(val))
    return texts


def match_feature(feature, compiled_patterns: List[re.Pattern], fields: List[str]) -> Optional[str]:
    texts = feature_texts(feature, fields)
    for text in texts:
        for pat in compiled_patterns:
            if pat.search(text):
                return text
    return None


def build_filter_terms(filters: Dict) -> List[str]:
    def as_str_list(value: object, name: str) -> List[str]:
        if value is None:
            return []
        if isinstance(value, str):
            return [value]
        if isinstance(value, list) and all(isinstance(v, str) for v in value):
            return value
        raise typer.BadParameter(f"filters.{name} must be a string or list of strings.")

    def as_date_str(value: object, name: str) -> Optional[str]:
        if value is None:
            return None
        if isinstance(value, str):
            return value
        if isinstance(value, (date, datetime)):
            return value.strftime("%Y/%m/%d")
        raise typer.BadParameter(f"filters.{name} must be a string or date.")

    terms: List[str] = []
    if filters is None:
        return terms
    if not isinstance(filters, dict):
        raise typer.BadParameter("[filters] must be a table (dict).")
    if not filters:
        return terms

    filters_list = as_str_list(filters.get("filter"), "filter")
    for value in filters_list:
        terms.append(f"{value}[filter]")

    properties_list = as_str_list(filters.get("properties"), "properties")
    for value in properties_list:
        terms.append(f"{value}[prop]")

    length_min = filters.get("sequence_length_min")
    length_max = filters.get("sequence_length_max")
    if length_min is not None or length_max is not None:
        try:
            lmin = int(length_min) if length_min is not None else 0
            lmax = int(length_max) if length_max is not None else 1000000000
        except (TypeError, ValueError) as exc:
            raise typer.BadParameter("filters.sequence_length_min/max must be integers.") from exc
        terms.append(f"{lmin}[SLEN] : {lmax}[SLEN]")

    date_from = as_date_str(filters.get("publication_date_from"), "publication_date_from")
    date_to = as_date_str(filters.get("publication_date_to"), "publication_date_to")
    if date_from or date_to:
        dfrom = date_from or "1800/01/01"
        dto = date_to or "3000/12/31"
        terms.append(f"{dfrom}[PDAT] : {dto}[PDAT]")

    date_from = as_date_str(filters.get("modification_date_from"), "modification_date_from")
    date_to = as_date_str(filters.get("modification_date_to"), "modification_date_to")
    if date_from or date_to:
        dfrom = date_from or "1800/01/01"
        dto = date_to or "3000/12/31"
        terms.append(f"{dfrom}[MDAT] : {dto}[MDAT]")

    include_keywords = as_str_list(filters.get("all_fields_include"), "all_fields_include")
    if include_keywords:
        inc = [f'"{k}"[All Fields]' for k in include_keywords]
        terms.append("(" + " OR ".join(inc) + ")")

    exclude_keywords = as_str_list(filters.get("all_fields_exclude"), "all_fields_exclude")
    if exclude_keywords:
        exc = [f'"{k}"[All Fields]' for k in exclude_keywords]
        terms.append("NOT (" + " OR ".join(exc) + ")")

    raw = filters.get("raw")
    if raw is not None:
        raw_terms = as_str_list(raw, "raw")
        for term in raw_terms:
            terms.append(str(term))

    return terms


def build_query(taxid: str, marker_query: str, filters: Dict, taxon_noexp: bool) -> str:
    tax_term = f"txid{taxid}[Organism:noexp]" if taxon_noexp else f"txid{taxid}[Organism]"
    parts = [tax_term, marker_query]
    parts.extend(build_filter_terms(filters))
    return " AND ".join(f"({p})" for p in parts)


def resolve_taxid(taxon: str) -> Tuple[str, Optional[str]]:
    if re.fullmatch(r"\d+", taxon):
        return taxon, None

    term = f'"{taxon}"[Scientific Name]'
    handle = Entrez.esearch(db="taxonomy", term=term, retmax=5)
    record = Entrez.read(handle)
    ids = record.get("IdList", [])
    if not ids:
        raise typer.BadParameter(f"Taxon not found in NCBI Taxonomy: {taxon}")

    if len(ids) > 1:
        return ids[0], f"Multiple taxids found for '{taxon}'. Using {ids[0]}."

    return ids[0], None


def default_delay(ncbi_cfg: Dict) -> float:
    delay = ncbi_cfg.get("delay_sec")
    if delay is not None:
        return float(delay)
    if getattr(Entrez, "api_key", None):
        return 0.11
    return 0.34


def fetch_genbank(
    query: str,
    ncbi_cfg: Dict,
    delay_sec: float,
    start_at: int = 0,
    dump_dir: Optional[Path] = None,
    resume: bool = False,
    taxid: Optional[str] = None,
) -> Tuple[int, Iterable[Tuple[int, str]]]:
    db = ncbi_cfg.get("db", "nucleotide")
    rettype = ncbi_cfg.get("rettype", "fasta")
    retmode = ncbi_cfg.get("retmode", "text")
    per_query = int(ncbi_cfg.get("per_query", 100))
    use_history = bool(ncbi_cfg.get("use_history", True))

    if rettype not in {"gb", "gbwithparts"}:
        raise typer.BadParameter("ncbi.rettype must be 'gb' or 'gbwithparts' for region extraction.")

    handle = Entrez.esearch(db=db, term=query, retmax=0, usehistory="y" if use_history else "n")
    record = Entrez.read(handle)
    count = int(record.get("Count", 0))
    webenv = record.get("WebEnv")
    query_key = record.get("QueryKey")

    if count == 0:
        return 0, []

    cache_root = None
    if dump_dir:
        cache_root = dump_dir / ".cache"
        cache_root.mkdir(parents=True, exist_ok=True)

    def dump_path_for(start: int) -> Optional[Path]:
        if not cache_root:
            return None
        return cache_root / f"start{start:09d}_count{per_query:04d}.cache"

    def load_cached(start: int) -> Optional[str]:
        path = dump_path_for(start)
        if path and path.exists():
            return path.read_text(encoding="utf-8", errors="ignore")
        return None

    def save_cached(start: int, data: str) -> None:
        path = dump_path_for(start)
        if path:
            path.write_text(data, encoding="utf-8")

    def fetch_chunk_by_ids(start: int) -> Optional[str]:
        search_handle = Entrez.esearch(
            db=db,
            term=query,
            retstart=start,
            retmax=per_query,
            usehistory="n",
        )
        search_record = Entrez.read(search_handle)
        ids = search_record.get("IdList", [])
        if not ids:
            return None
        fetch_handle = Entrez.efetch(
            db=db,
            rettype=rettype,
            retmode=retmode,
            id=",".join(ids),
        )
        return fetch_handle.read()

    def gen() -> Iterable[Tuple[int, str]]:
        if use_history and webenv and query_key:
            for start in range(start_at, count, per_query):
                cached = load_cached(start) if resume else None
                if cached is not None:
                    yield start, cached
                    continue
                try:
                    fetch_handle = Entrez.efetch(
                        db=db,
                        rettype=rettype,
                        retmode=retmode,
                        retstart=start,
                        retmax=per_query,
                        webenv=webenv,
                        query_key=query_key,
                    )
                    data = fetch_handle.read()
                except HTTPError as exc:
                    if exc.code != HTTPStatus.BAD_REQUEST:
                        raise
                    fallback = fetch_chunk_by_ids(start)
                    if fallback is None:
                        continue
                    data = fallback
                if cache_root:
                    save_cached(start, data)
                yield start, data
                time.sleep(delay_sec)
            return

        for start in range(start_at, count, per_query):
            cached = load_cached(start) if resume else None
            if cached is not None:
                yield start, cached
                continue
            data = fetch_chunk_by_ids(start)
            if data is None:
                continue
            if cache_root:
                save_cached(start, data)
            yield start, data
            time.sleep(delay_sec)

    return count, gen()


def iter_genbank_files(from_dir: Path, taxid: Optional[str] = None) -> Iterable[Tuple[int, str]]:
    gb_root = from_dir
    if taxid:
        candidate = from_dir / f"taxid{taxid}"
        if candidate.exists():
            gb_root = candidate
    files = sorted(gb_root.glob("*.gb"))
    for idx, path in enumerate(files):
        yield idx, path.read_text(encoding="utf-8", errors="ignore")


def build_output_path(
    out: Optional[Path],
    taxids: List[str],
    markers: List[str],
    output_prefix: str = "",
) -> Path:
    run_date = datetime.now().strftime("%Y%m%d")
    if out:
        if out.suffix in {".fa", ".fasta", ".fas"}:
            out.parent.mkdir(parents=True, exist_ok=True)
            return out
        out_dir = out
    else:
        out_dir = Path("Results") / "db" / run_date

    out_dir.mkdir(parents=True, exist_ok=True)
    taxon_label = f"taxid{'+'.join(taxids)}" if len(taxids) == 1 else "multi_taxon"
    marker_label = "+".join(markers) if len(markers) == 1 else "multi_marker"
    prefix = output_prefix
    if prefix and not prefix.endswith("_"):
        prefix = prefix + "_"
    return out_dir / f"{prefix}{taxon_label}__{marker_label}.fasta"


def setup_run_logger(log_path: Path) -> logging.Logger:
    logger_name = f"taxondbbuilder.run.{os.getpid()}.{int(time.time() * 1_000_000)}"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(handler)
    return logger


def close_run_logger(logger: logging.Logger) -> None:
    for handler in list(logger.handlers):
        handler.flush()
        handler.close()
        logger.removeHandler(handler)


def resolve_support_file_path(raw_path: str, config_path: Path, label: str) -> Path:
    path = Path(os.path.expandvars(os.path.expanduser(raw_path)))
    candidates: List[Path] = []
    if path.is_absolute():
        candidates.append(path)
    else:
        candidates.append(config_path.parent / path)
        candidates.append(Path.cwd() / path)
        candidates.append(Path(__file__).resolve().parent / path)

    resolved = next((p for p in candidates if p.exists()), None)
    if not resolved:
        tried = ", ".join(str(p) for p in candidates)
        raise typer.BadParameter(f"{label} not found. Tried: {tried}")
    return resolved


def normalize_primer_values(values: List[str], primer_set_name: str, field_name: str) -> List[str]:
    normalized: List[str] = []
    for value in values:
        primer = value.strip().upper().replace("U", "T")
        if not primer:
            raise typer.BadParameter(f"primer_sets.{primer_set_name}.{field_name} cannot contain empty primers.")
        invalid = sorted({ch for ch in primer if ch not in IUPAC_DNA_VALUES})
        if invalid:
            chars = "".join(invalid)
            raise typer.BadParameter(
                f"primer_sets.{primer_set_name}.{field_name} contains unsupported IUPAC chars: {chars}"
            )
        normalized.append(primer)
    return normalized


def load_primer_sets_from_file(primer_path: Path) -> Dict[str, Dict[str, List[str]]]:
    with primer_path.open("rb") as f:
        primer_data = tomllib.load(f)
    primer_sets = primer_data.get("primer_sets")
    if not isinstance(primer_sets, dict) or not primer_sets:
        raise typer.BadParameter("Primer file must define a non-empty [primer_sets] section.")

    normalized_sets: Dict[str, Dict[str, List[str]]] = {}
    for set_name, set_cfg in primer_sets.items():
        if not isinstance(set_cfg, dict):
            raise typer.BadParameter(f"primer_sets.{set_name} must be a table (dict).")
        forward = set_cfg.get("forward")
        reverse = set_cfg.get("reverse")
        if not isinstance(forward, list) or not forward or not all(isinstance(v, str) for v in forward):
            raise typer.BadParameter(f"primer_sets.{set_name}.forward must be a non-empty list of strings.")
        if not isinstance(reverse, list) or not reverse or not all(isinstance(v, str) for v in reverse):
            raise typer.BadParameter(f"primer_sets.{set_name}.reverse must be a non-empty list of strings.")
        normalized_sets[set_name] = {
            "forward": normalize_primer_values(forward, set_name, "forward"),
            "reverse": normalize_primer_values(reverse, set_name, "reverse"),
        }
    return normalized_sets


def combine_primer_set_sequences(
    primer_sets_data: Dict[str, Dict[str, List[str]]],
    selected_sets: List[str],
) -> Tuple[List[str], List[str]]:
    forward: List[str] = []
    reverse: List[str] = []
    for set_name in selected_sets:
        set_data = primer_sets_data.get(set_name)
        if not set_data:
            available = ", ".join(sorted(primer_sets_data.keys()))
            raise typer.BadParameter(f"Primer set '{set_name}' was not found. Available: {available}")
        for p in set_data["forward"]:
            if p not in forward:
                forward.append(p)
        for p in set_data["reverse"]:
            if p not in reverse:
                reverse.append(p)
    return forward, reverse


@app.command("list-markers")
def list_markers(
    config: Path = typer.Option(..., "--config", "-c", help="Path to TOML config file."),
) -> None:
    """
    List marker IDs and aliases from the config (including markers.file).
    """
    cfg = load_config(config)
    markers = cfg.get("markers", {})
    table = Table(title="Markers", show_header=True, header_style="bold")
    table.add_column("Marker ID")
    table.add_column("Aliases")
    table.add_column("Header Format")
    for key in sorted(markers.keys()):
        entry = markers[key] or {}
        aliases = entry.get("aliases") or []
        header_format = entry.get("header_format") or ""
        aliases_text = ", ".join(aliases) if aliases else "-"
        table.add_row(str(key), aliases_text, str(header_format) or "-")
    console.print(table)


@app.command("list-primer-sets")
def list_primer_sets(
    config: Path = typer.Option(..., "--config", "-c", help="Path to TOML config file."),
) -> None:
    """
    List primer set IDs from [post_prep].primer_file.
    """
    if not config.exists():
        raise typer.BadParameter(f"Config file not found: {config}")
    with config.open("rb") as f:
        cfg = tomllib.load(f)

    post_prep_cfg = cfg.get("post_prep") or {}
    if not isinstance(post_prep_cfg, dict):
        raise typer.BadParameter("[post_prep] must be a table (dict).")
    primer_file = post_prep_cfg.get("primer_file")
    if not primer_file:
        raise typer.BadParameter(
            "[post_prep].primer_file is not set in config. Set it to use list-primer-sets."
        )
    primer_path = resolve_support_file_path(str(primer_file), config, "Primer file")
    primer_sets = load_primer_sets_from_file(primer_path)

    table = Table(title=f"Primer Sets ({primer_path})", show_header=True, header_style="bold")
    table.add_column("Primer Set")
    table.add_column("Forward")
    table.add_column("Reverse")
    for key in sorted(primer_sets.keys()):
        entry = primer_sets[key]
        table.add_row(key, str(len(entry["forward"])), str(len(entry["reverse"])))
    console.print(table)


@app.command()
def build(
    config: Path = typer.Option(..., "--config", "-c", help="Path to TOML config file."),
    taxon: List[str] = typer.Option(..., "--taxon", "-t", help="Taxon (taxid or scientific name)."),
    marker: List[str] = typer.Option(..., "--marker", "-m", help="Marker key or prefix."),
    out: Optional[Path] = typer.Option(None, "--out", "-o", help="Output file or directory."),
    dump_gb: Optional[Path] = typer.Option(
        None,
        "--dump-gb",
        help="Directory to store GenBank chunks for caching.",
    ),
    from_gb: Optional[Path] = typer.Option(
        None,
        "--from-gb",
        help="Directory of GenBank chunks to extract without downloading.",
    ),
    resume: bool = typer.Option(False, "--resume", help="Resume using cached GenBank chunks."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print query and exit."),
    workers: int = typer.Option(2, "--workers", "-w", help="Number of extraction workers."),
    output_prefix: str = typer.Option(
        "taxondbbuilder_",
        "--output-prefix",
        help="Prefix added to output FASTA filename.",
    ),
    post_prep: bool = typer.Option(
        False,
        "--post-prep",
        help="Apply [post_prep] FASTA processing after extraction.",
    ),
    post_prep_step: Optional[List[PostPrepStep]] = typer.Option(
        None,
        "--post-prep-step",
        help=(
            "Post-prep step(s) to run. Repeat to select multiple. "
            "Choices: primer_trim, length_filter, duplicate_report."
        ),
    ),
    post_prep_primer_set: Optional[List[str]] = typer.Option(
        None,
        "--post-prep-primer-set",
        help="Primer set name(s) for primer_trim. Repeat to select multiple (overrides config).",
    ),
):
    """
    Build a FASTA database by downloading GenBank records and extracting features.

    Examples (Teleostomi):
      taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s
      taxondbbuilder.py build -c configs/db.toml -t "Salmo salar" -m mitogenome
      taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s --workers 2
      taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s --output-prefix "mifish"
      taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s --dump-gb Results/gb --resume
      taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s --from-gb Results/gb
    """
    cfg = load_config(config)
    ncbi_cfg = cfg.get("ncbi", {})
    filters_cfg = cfg.get("filters", {})
    output_cfg = cfg.get("output", {})
    post_prep_cfg = cfg.get("post_prep") or {}
    taxon_noexp = bool(cfg.get("taxon", {}).get("noexp", False))

    setup_entrez(ncbi_cfg)
    marker_map = normalize_marker_map(cfg.get("markers", {}))

    marker_keys = [resolve_marker_key(m, marker_map) for m in marker]
    marker_query = build_marker_query(marker_keys, marker_map)
    output_prefix = output_prefix.strip()
    selected_header_formats: List[str] = []
    marker_rules = []
    for key in marker_keys:
        cfg_m = marker_map[key]
        region_patterns = build_region_patterns(cfg_m)
        if not region_patterns:
            raise typer.BadParameter(f"markers.{key} has no patterns for region extraction.")
        compiled = compile_patterns(region_patterns)
        if not compiled:
            raise typer.BadParameter(f"markers.{key} patterns did not compile.")

        feature_types = cfg_m.get("feature_types")
        if feature_types is None:
            feature_types = DEFAULT_FEATURE_TYPES
        elif not feature_types:
            feature_types = None

        feature_fields = cfg_m.get("feature_fields")
        if feature_fields is None:
            feature_fields = DEFAULT_FEATURE_FIELDS
        elif not feature_fields:
            raise typer.BadParameter(f"markers.{key}.feature_fields cannot be empty.")

        marker_rules.append(
            {
                "key": key,
                "patterns": compiled,
                "feature_types": feature_types,
                "feature_fields": feature_fields,
                "header_format": resolve_header_format(cfg_m, output_cfg),
            }
        )
        selected_header_formats.append(marker_rules[-1]["header_format"])

    requested_post_prep_steps = [step.value for step in (post_prep_step or [])]
    requested_primer_sets: List[str] = []
    for value in (post_prep_primer_set or []):
        name = value.strip()
        if not name:
            raise typer.BadParameter("--post-prep-primer-set cannot include empty values.")
        if name not in requested_primer_sets:
            requested_primer_sets.append(name)

    if post_prep:
        post_min = post_prep_cfg.get("sequence_length_min")
        post_max = post_prep_cfg.get("sequence_length_max")
        has_length_filter = post_min is not None or post_max is not None
        post_primer_forward = list(post_prep_cfg.get("_primer_forward") or [])
        post_primer_reverse = list(post_prep_cfg.get("_primer_reverse") or [])
        post_primer_file = post_prep_cfg.get("_primer_file_resolved") or post_prep_cfg.get("primer_file")
        post_primer_set_names = list(post_prep_cfg.get("_primer_set_names") or [])
        if not post_primer_set_names:
            primer_set_cfg = post_prep_cfg.get("primer_set")
            if isinstance(primer_set_cfg, str) and primer_set_cfg.strip():
                post_primer_set_names = [primer_set_cfg.strip()]
            elif isinstance(primer_set_cfg, list):
                post_primer_set_names = [str(v).strip() for v in primer_set_cfg if str(v).strip()]

        if requested_primer_sets:
            if not post_primer_file:
                raise typer.BadParameter(
                    "--post-prep-primer-set requires [post_prep].primer_file in config."
                )
            primer_path = resolve_support_file_path(str(post_primer_file), config, "Primer file")
            primer_sets_data = load_primer_sets_from_file(primer_path)
            post_primer_forward, post_primer_reverse = combine_primer_set_sequences(
                primer_sets_data, requested_primer_sets
            )
            post_primer_file = str(primer_path)
            post_primer_set_names = requested_primer_sets

        has_primer_trim = bool(post_primer_forward and post_primer_reverse)

        if PostPrepStep.PRIMER_TRIM.value in requested_post_prep_steps and not has_primer_trim:
            raise typer.BadParameter(
                "post-prep step 'primer_trim' requires post_prep.primer_file and post_prep.primer_set."
            )
        if PostPrepStep.LENGTH_FILTER.value in requested_post_prep_steps and not has_length_filter:
            raise typer.BadParameter(
                "post-prep step 'length_filter' requires post_prep.sequence_length_min or post_prep.sequence_length_max."
            )

        if requested_post_prep_steps:
            post_prep_steps_run = [
                step for step in POST_PREP_STEP_ORDER if step in requested_post_prep_steps
            ]
        else:
            post_prep_steps_run = []
            if has_primer_trim:
                post_prep_steps_run.append(PostPrepStep.PRIMER_TRIM.value)
            if has_length_filter:
                post_prep_steps_run.append(PostPrepStep.LENGTH_FILTER.value)
            post_prep_steps_run.append(PostPrepStep.DUPLICATE_REPORT.value)

        post_min = int(post_min) if post_min is not None else None
        post_max = int(post_max) if post_max is not None else None
    else:
        if requested_post_prep_steps or requested_primer_sets:
            raise typer.BadParameter("--post-prep-step/--post-prep-primer-set requires --post-prep.")
        has_length_filter = False
        has_primer_trim = False
        post_prep_steps_run = []
        post_min = None
        post_max = None
        post_primer_forward = []
        post_primer_reverse = []
        post_primer_file = None
        post_primer_set_names = []

    taxids: List[str] = []
    warnings: List[str] = []
    for t in taxon:
        taxid, warn = resolve_taxid(t)
        taxids.append(taxid)
        if warn:
            warnings.append(warn)

    out_path = build_output_path(out, taxids, marker_keys, output_prefix=output_prefix)
    log_path = out_path.with_suffix(out_path.suffix + ".log")
    if resume and not dump_gb and not from_gb:
        raise typer.BadParameter("--resume requires --dump-gb or --from-gb.")
    if from_gb and not from_gb.exists():
        raise typer.BadParameter(f"--from-gb not found: {from_gb}")
    if dump_gb:
        dump_gb.mkdir(parents=True, exist_ok=True)

    print_header()
    render_run_table(
        config,
        taxids,
        marker_keys,
        out_path,
        filters_cfg,
        dump_gb=dump_gb,
        from_gb=from_gb,
        resume=resume,
    )
    for w in warnings:
        console.print(f"[yellow]WARNING:[/yellow] {w}")

    if dry_run:
        for taxid in taxids:
            query = build_query(taxid, marker_query, filters_cfg, taxon_noexp)
            console.print(query)
        return

    acc_to_seqs: Dict[str, set] = {}
    counters = {
        "total_records": 0,
        "matched_records": 0,
        "matched_features": 0,
        "kept_records": 0,
        "skipped_same": 0,
        "duplicated_diff": 0,
    }
    dup_accessions: Dict[str, int] = {}
    emitted_records: List[Dict[str, str]] = []

    progress = Progress(
        SpinnerColumn(),
        TextColumn("[bold]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
        disable=not console.is_terminal,
    )

    run_logger = setup_run_logger(log_path)
    try:
        run_logger.info(f"# started: {datetime.now().isoformat()}")
        run_logger.info(f"# config: {config}")
        run_logger.info(f"# taxon input: {taxon}")
        run_logger.info(f"# taxids: {taxids}")
        run_logger.info(f"# markers: {marker_keys}")
        run_logger.info(f"# output_prefix: {output_prefix}")
        run_logger.info(f"# dump_gb: {dump_gb}" if dump_gb else "# dump_gb: none")
        run_logger.info(f"# from_gb: {from_gb}" if from_gb else "# from_gb: none")
        run_logger.info(f"# resume: {resume}")
        run_logger.info(f"# post_prep: {post_prep}")
        if post_prep:
            steps_text = ", ".join(post_prep_steps_run) if post_prep_steps_run else "none"
            run_logger.info(f"# post_prep.steps: {steps_text}")
        if post_prep and has_length_filter:
            if post_min is not None:
                run_logger.info(f"# post_prep.sequence_length_min: {post_min}")
            if post_max is not None:
                run_logger.info(f"# post_prep.sequence_length_max: {post_max}")
        if post_prep and has_primer_trim:
            run_logger.info(f"# post_prep.primer_file: {post_primer_file}")
            run_logger.info(f"# post_prep.primer_set: {', '.join(post_primer_set_names)}")
            run_logger.info(f"# post_prep.primer_forward_count: {len(post_primer_forward)}")
            run_logger.info(f"# post_prep.primer_reverse_count: {len(post_primer_reverse)}")
        if warnings:
            run_logger.info("# warnings:")
            for warning in warnings:
                run_logger.info(f"# - {warning}")

        lock = Lock()
        run_logger.info(f"# workers: {workers}")

        with out_path.open("w", encoding="utf-8") as out_f, progress:
            for taxid in taxids:
                query = build_query(taxid, marker_query, filters_cfg, taxon_noexp)
                run_logger.info(f"# query taxid={taxid}: {query}")
                delay_sec = default_delay(ncbi_cfg)
                if from_gb:
                    data_iter = iter_genbank_files(from_gb, taxid)
                    count = None
                    run_logger.info(f"# query count taxid={taxid}: from-gb")
                else:
                    count, data_iter = fetch_genbank(
                        query,
                        ncbi_cfg,
                        delay_sec,
                        dump_dir=dump_gb,
                        resume=resume,
                        taxid=taxid,
                    )
                    run_logger.info(f"# query count taxid={taxid}: {count}")
                    if count == 0:
                        console.print(f"[yellow]taxid {taxid}: 0 records[/yellow]")
                        continue
                    run_logger.info(f"# fetch progress taxid={taxid}: 0/{count}")

                task_id = progress.add_task(f"taxid {taxid}", total=count)
                if workers < 1:
                    raise typer.BadParameter("--workers must be >= 1.")

                q: Queue = Queue(maxsize=max(1, workers * 2))
                stop_event = Event()
                errors: List[Exception] = []

                def worker() -> None:
                    while True:
                        item = q.get()
                        if item is None:
                            q.task_done()
                            break
                        try:
                            start, chunk = item
                            process_genbank_chunk(
                                chunk,
                                marker_rules,
                                acc_to_seqs,
                                out_f,
                                counters,
                                dup_accessions,
                                lock,
                                progress,
                                task_id,
                                taxid,
                                dump_gb,
                                emitted_records,
                            )
                        except Exception as exc:
                            errors.append(exc)
                            stop_event.set()
                        finally:
                            q.task_done()

                threads = [Thread(target=worker, daemon=True) for _ in range(workers)]
                for t in threads:
                    t.start()

                for start, chunk in data_iter:
                    if stop_event.is_set():
                        break
                    if not chunk:
                        continue
                    if count is not None and count > 0:
                        per_query = int(ncbi_cfg.get("per_query", 100))
                        fetched = min(start + per_query, count)
                        run_logger.info(f"# fetch progress taxid={taxid}: {fetched}/{count}")
                    q.put((start, chunk))

                for _ in threads:
                    q.put(None)
                q.join()
                for t in threads:
                    t.join()
                if errors:
                    raise errors[0]

        duplicate_records_report_path: Optional[Path] = None
        duplicate_groups_report_path: Optional[Path] = None
        if post_prep:
            before_post_prep = counters["kept_records"]
            run_logger.info(f"# kept records before post_prep: {before_post_prep}")

            if PostPrepStep.PRIMER_TRIM.value in post_prep_steps_run:
                primer_stats = apply_post_prep_primer_trim(
                    out_path,
                    post_primer_forward,
                    post_primer_reverse,
                )
                counters["kept_records"] = primer_stats["after"]
                run_logger.info(
                    "# post_prep primer trim:"
                    f" before={primer_stats['before']} after={primer_stats['after']}"
                    f" removed={primer_stats['removed']} trimmed_both={primer_stats['trimmed_both']}"
                    f" trimmed_left_only={primer_stats['trimmed_left_only']}"
                    f" trimmed_right_only={primer_stats['trimmed_right_only']}"
                    f" untrimmed={primer_stats['untrimmed']}"
                    f" dropped_empty={primer_stats['dropped_empty']}"
                    f" canonical_orientation={primer_stats['canonical_orientation']}"
                    f" reverse_orientation={primer_stats['reverse_orientation']}"
                )

            if PostPrepStep.LENGTH_FILTER.value in post_prep_steps_run:
                length_stats = apply_post_prep_length_filter(out_path, post_min, post_max)
                counters["kept_records"] = length_stats["after"]
                run_logger.info(
                    "# post_prep length filter:"
                    f" before={length_stats['before']} after={length_stats['after']} removed={length_stats['removed']}"
                )

            if PostPrepStep.DUPLICATE_REPORT.value in post_prep_steps_run:
                (
                    duplicate_records_report_path,
                    duplicate_groups_report_path,
                    dup_stats,
                    dup_reason,
                ) = write_duplicate_acc_reports_csv(out_path, selected_header_formats)
                if dup_reason:
                    run_logger.info(f"# post_prep duplicate_acc_report: skipped ({dup_reason})")
                    console.print(f"[yellow]post_prep:[/yellow] duplicate ACC report skipped ({dup_reason}).")
                else:
                    run_logger.info(
                        "# post_prep duplicate_acc_report:"
                        f" total={dup_stats['total_records']} parsed={dup_stats['parsed_records']}"
                        f" unparsed={dup_stats['unparsed_records']} groups={dup_stats['duplicate_groups']}"
                        f" records={dup_stats['duplicate_records']}"
                        f" cross_organism_groups={dup_stats['cross_organism_groups']}"
                    )
                    if duplicate_records_report_path:
                        run_logger.info(f"# post_prep duplicate_acc_records_csv: {duplicate_records_report_path}")
                        console.print(f"post_prep duplicate ACC records CSV: {duplicate_records_report_path}")
                    if duplicate_groups_report_path:
                        run_logger.info(f"# post_prep duplicate_acc_groups_csv: {duplicate_groups_report_path}")
                        console.print(f"post_prep duplicate ACC groups CSV: {duplicate_groups_report_path}")
            else:
                run_logger.info("# post_prep duplicate_acc_report: skipped (step disabled)")

        acc_species_map_path, acc_species_stats = write_acc_organism_mapping_csv(out_path, emitted_records)
        run_logger.info(
            "# acc_organism_map:"
            f" total={acc_species_stats['total_records']} mapped={acc_species_stats['mapped_records']}"
            f" unmapped={acc_species_stats['unmapped_records']}"
            f" unused_source_records={acc_species_stats['unused_records']}"
            f" unique_accessions={acc_species_stats['unique_accessions']}"
            f" unique_organisms={acc_species_stats['unique_organisms']}"
        )
        run_logger.info(f"# acc_organism_map_csv: {acc_species_map_path}")

        run_logger.info(f"# total records: {counters['total_records']}")
        run_logger.info(f"# matched records: {counters['matched_records']}")
        run_logger.info(f"# matched features: {counters['matched_features']}")
        run_logger.info(f"# kept records: {counters['kept_records']}")
        run_logger.info(f"# skipped duplicates (same accession+sequence): {counters['skipped_same']}")
        run_logger.info(f"# kept duplicates (same accession, different sequence): {counters['duplicated_diff']}")
        if dup_accessions:
            run_logger.info("# duplicate accessions with different sequences:")
            for acc, count in sorted(dup_accessions.items()):
                run_logger.info(f"# - {acc}: {count} sequences")
        run_logger.info(f"# output: {out_path}")
        run_logger.info(f"# finished: {datetime.now().isoformat()}")
    finally:
        close_run_logger(run_logger)

    render_result_table(
        counters["total_records"],
        counters["matched_records"],
        counters["matched_features"],
        counters["kept_records"],
        counters["skipped_same"],
        counters["duplicated_diff"],
        out_path,
        log_path,
    )
    if dup_accessions:
        console.print(
            "[yellow]WARNING:[/yellow] duplicate accessions with different sequences were kept. See log for details."
        )
    console.print(f"ACC-organism mapping CSV: {acc_species_map_path}")
    if acc_species_stats["unmapped_records"] > 0:
        console.print(
            "[yellow]WARNING:[/yellow] Some final FASTA records could not be mapped to source ACC/organism. "
            "See .log for details."
        )


if __name__ == "__main__":
    app()

import io
import json
import logging
import os
import re
import subprocess
import tempfile
import time
import csv
import hashlib
import shutil
import sys
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from dataclasses import dataclass
from datetime import date, datetime
from enum import Enum
from http import HTTPStatus
from http.client import HTTPException, RemoteDisconnected
from pathlib import Path
from queue import Queue
from string import Formatter
from threading import Event, Lock, Thread
from typing import Any, Dict, Iterable, List, Optional, Tuple
from urllib.error import HTTPError, URLError

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
from taxondb_bold import BoldApiError, fetch_bold_records_for_taxon, parse_accession_tokens

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
DEFAULT_BOLD_HEADER_FORMAT = "bold|{acc_id}|{organism}"
FORMATTER = Formatter()
IUPAC_DNA_VALUES = {k.upper(): v.upper() for k, v in ambiguous_dna_values.items()}
IUPAC_DNA_VALUES["U"] = "T"


class PostPrepStep(str, Enum):
    PRIMER_TRIM = "primer_trim"
    LENGTH_FILTER = "length_filter"
    DUPLICATE_REPORT = "duplicate_report"


class BuildSource(str, Enum):
    NCBI = "ncbi"
    BOLD = "bold"
    BOTH = "both"


@dataclass(frozen=True)
class ResolvedTaxon:
    input_value: str
    taxid: str
    scientific_name: str
    warning: Optional[str] = None


@dataclass
class CanonicalRecord:
    source: str
    source_record_id: str
    accession: Optional[str]
    processid: Optional[str]
    sampleid: Optional[str]
    taxon_name: Optional[str]
    marker_key: str
    marker_label: Optional[str]
    sequence: str
    header_values: Dict[str, str]
    metadata: Dict[str, str]
    linked_to_ncbi: bool = False
    emitted_to_fasta: bool = True
    skip_reason: Optional[str] = None


POST_PREP_STEP_ORDER = [
    PostPrepStep.PRIMER_TRIM.value,
    PostPrepStep.LENGTH_FILTER.value,
    PostPrepStep.DUPLICATE_REPORT.value,
]

PRIMER_TRIM_MODE_BOTH_REQUIRED = "both_required"
PRIMER_TRIM_MODE_ONE_OR_BOTH = "one_or_both"
PRIMER_TRIM_MODE_ONE_END_MARK_ONLY = "one_end_mark_only"
PRIMER_TRIM_MODES = {
    PRIMER_TRIM_MODE_BOTH_REQUIRED,
    PRIMER_TRIM_MODE_ONE_OR_BOTH,
    PRIMER_TRIM_MODE_ONE_END_MARK_ONLY,
}


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
    source: BuildSource,
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
    table.add_row("Source", source.value)
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


def extract_ncbi_records_from_genbank_chunk(
    chunk: str,
    marker_rules: List[Dict],
    acc_to_seqs: Dict[str, set],
    counters: Dict[str, int],
    dup_accessions: Dict[str, int],
    lock: Lock,
    progress: Progress,
    task_id: int,
    taxid: str,
    dump_gb_dir: Optional[Path],
    source: BuildSource = BuildSource.NCBI,
) -> List[CanonicalRecord]:
    extracted_records: List[CanonicalRecord] = []
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
                    "source": source.value,
                    "source_id": acc,
                }
                extracted_records.append(
                    CanonicalRecord(
                        source=source.value,
                        source_record_id=acc,
                        accession=acc,
                        processid=None,
                        sampleid=None,
                        taxon_name=str(organism),
                        marker_key=str(matched_marker or ""),
                        marker_label=str(matched_label),
                        sequence=seq,
                        header_values=header_values,
                        metadata={
                            "organism_name": str(organism),
                            "matched_type": str(matched_type),
                            "header_format": header_format or DEFAULT_HEADER_FORMAT,
                        },
                    )
                )
        if record_matched:
            with lock:
                counters["matched_records"] += 1
    return extracted_records


def emit_records_to_fasta(
    records: List[CanonicalRecord],
    out_f,
    counters: Dict[str, int],
    emitted_records: List[Dict[str, str]],
    lock: Lock,
    source_merge_rows: Optional[List[Dict[str, str]]] = None,
) -> None:
    for record in records:
        if not record.emitted_to_fasta:
            continue

        header_template = record.metadata.get("header_format") or DEFAULT_HEADER_FORMAT
        header = build_header(header_template, record.header_values).strip()
        if not header:
            header = record.header_values.get("acc_id") or record.source_record_id

        with lock:
            out_f.write(f">{header}\n")
            out_f.write(f"{record.sequence}\n")
            counters["kept_records"] += 1
            emitted_row = build_source_merge_row(record, header=header)
            emitted_records.append(dict(emitted_row))
            if source_merge_rows is not None:
                source_merge_rows.append(dict(emitted_row))


def build_bold_canonical_record(
    normalized_row: Dict[str, Any],
    marker_map: Dict[str, Dict],
    output_cfg: Dict,
) -> CanonicalRecord:
    marker_key = str(normalized_row["marker_key"])
    marker_cfg = marker_map[marker_key]
    header_key = marker_cfg.get("header_format")
    if not header_key:
        header_format = DEFAULT_BOLD_HEADER_FORMAT
    elif header_key == "mifish_pipeline":
        header_format = DEFAULT_BOLD_HEADER_FORMAT
    else:
        header_format = resolve_header_format(marker_cfg, output_cfg)
    processid = normalized_row.get("processid")
    sampleid = normalized_row.get("sampleid")
    accession = normalized_row.get("accession")
    taxon_name = normalized_row.get("taxon_name") or "unknown"
    source_record_id = str(normalized_row["source_record_id"])
    acc_id = f"BOLD_{sanitize_header(processid or source_record_id)}"
    marker_label = str(normalized_row.get("marker_label") or "")
    marker_safe = sanitize_header(marker_key)
    marker_label_safe = sanitize_header(marker_label or marker_key)
    organism_safe = sanitize_header(str(taxon_name))

    header_values = {
        "acc": str(accession or ""),
        "acc_id": acc_id,
        "organism": organism_safe,
        "organism_raw": str(taxon_name),
        "marker": marker_safe,
        "marker_raw": marker_key,
        "label": marker_label_safe,
        "label_raw": marker_label,
        "type": "barcode",
        "type_raw": "barcode",
        "start": "",
        "end": "",
        "loc": "",
        "strand": "",
        "dup": "",
        "source": BuildSource.BOLD.value,
        "source_id": source_record_id,
    }
    return CanonicalRecord(
        source=BuildSource.BOLD.value,
        source_record_id=source_record_id,
        accession=str(accession or "") or None,
        processid=str(processid or "") or None,
        sampleid=str(sampleid or "") or None,
        taxon_name=str(taxon_name),
        marker_key=marker_key,
        marker_label=marker_label,
        sequence=str(normalized_row["sequence"]),
        header_values=header_values,
        metadata={
            "header_format": header_format,
            "raw_row_json": json.dumps(normalized_row.get("raw_row") or {}, ensure_ascii=False),
        },
    )


def build_source_merge_row(record: CanonicalRecord, header: str = "") -> Dict[str, str]:
    return {
        "source": record.source,
        "source_record_id": record.source_record_id,
        "acc_id": record.header_values.get("acc_id", ""),
        "accession": record.accession or "",
        "processid": record.processid or "",
        "sampleid": record.sampleid or "",
        "organism_name": record.taxon_name or "",
        "marker_key": record.marker_key,
        "marker_label": record.marker_label or "",
        "linked_to_ncbi": str(record.linked_to_ncbi).lower(),
        "emitted_to_fasta": str(record.emitted_to_fasta).lower(),
        "skip_reason": record.skip_reason or "",
        "header": header,
    }


def canonical_record_sort_key(record: CanonicalRecord) -> Tuple[str, str, str]:
    return (
        record.taxon_name or "",
        record.marker_key,
        record.source_record_id,
    )


def canonical_record_to_dict(record: CanonicalRecord) -> Dict[str, Any]:
    return {
        "source": record.source,
        "source_record_id": record.source_record_id,
        "accession": record.accession,
        "processid": record.processid,
        "sampleid": record.sampleid,
        "taxon_name": record.taxon_name,
        "marker_key": record.marker_key,
        "marker_label": record.marker_label,
        "sequence": record.sequence,
        "header_values": record.header_values,
        "metadata": record.metadata,
        "linked_to_ncbi": record.linked_to_ncbi,
        "emitted_to_fasta": record.emitted_to_fasta,
        "skip_reason": record.skip_reason,
    }


def canonical_record_from_dict(data: Dict[str, Any]) -> CanonicalRecord:
    return CanonicalRecord(
        source=str(data.get("source", "")),
        source_record_id=str(data.get("source_record_id", "")),
        accession=data.get("accession"),
        processid=data.get("processid"),
        sampleid=data.get("sampleid"),
        taxon_name=data.get("taxon_name"),
        marker_key=str(data.get("marker_key", "")),
        marker_label=data.get("marker_label"),
        sequence=str(data.get("sequence", "")),
        header_values=dict(data.get("header_values") or {}),
        metadata=dict(data.get("metadata") or {}),
        linked_to_ncbi=bool(data.get("linked_to_ncbi", False)),
        emitted_to_fasta=bool(data.get("emitted_to_fasta", True)),
        skip_reason=data.get("skip_reason"),
    )


def append_records_to_spool(records: List[CanonicalRecord], spool_f, lock: Lock) -> None:
    if not records:
        return
    with lock:
        for record in records:
            spool_f.write(json.dumps(canonical_record_to_dict(record), ensure_ascii=False) + "\n")


def load_records_from_spool(spool_path: Path) -> List[CanonicalRecord]:
    records: List[CanonicalRecord] = []
    if not spool_path.exists():
        return records
    with spool_path.open("r", encoding="utf-8") as in_f:
        for line in in_f:
            raw = line.strip()
            if not raw:
                continue
            records.append(canonical_record_from_dict(json.loads(raw)))
    return records


def write_source_merge_csv(fasta_path: Path, rows: List[Dict[str, str]]) -> Path:
    merge_path = fasta_path.with_suffix(fasta_path.suffix + ".source_merge.csv")
    sorted_rows = sorted(
        rows,
        key=lambda row: (
            row.get("source", ""),
            row.get("organism_name", ""),
            row.get("marker_key", ""),
            row.get("source_record_id", ""),
        ),
    )
    with merge_path.open("w", newline="", encoding="utf-8") as out_f:
        writer = csv.DictWriter(
            out_f,
            fieldnames=[
                "source",
                "source_record_id",
                "acc_id",
                "accession",
                "processid",
                "sampleid",
                "organism_name",
                "marker_key",
                "marker_label",
                "linked_to_ncbi",
                "emitted_to_fasta",
                "skip_reason",
                "header",
            ],
        )
        writer.writeheader()
        writer.writerows(sorted_rows)
    return merge_path



def load_config(path: Path, source: BuildSource = BuildSource.NCBI) -> Dict:
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

    requires_ncbi = source in {BuildSource.NCBI, BuildSource.BOTH}
    if requires_ncbi and "ncbi" not in data:
        raise typer.BadParameter("Missing [ncbi] section in config.")
    if "markers" not in data or not data["markers"]:
        raise typer.BadParameter("Missing [markers] section in config.")
    bold_cfg = data.get("bold")
    if bold_cfg is not None and not isinstance(bold_cfg, dict):
        raise typer.BadParameter("[bold] must be a table (dict).")
    if source == BuildSource.BOTH and bold_cfg is None:
        data["bold"] = {}

    post_prep = data.get("post_prep")
    if post_prep is not None:
        if not isinstance(post_prep, dict):
            raise typer.BadParameter("[post_prep] must be a table (dict).")

        def parse_int_option(name: str, raw: Any, min_value: Optional[int] = None) -> int:
            try:
                value = int(raw)
            except (TypeError, ValueError) as exc:
                raise typer.BadParameter(f"{name} must be an integer.") from exc
            if min_value is not None and value < min_value:
                raise typer.BadParameter(f"{name} must be >= {min_value}.")
            return value

        def parse_float_option(
            name: str,
            raw: Any,
            min_value: Optional[float] = None,
            max_value: Optional[float] = None,
        ) -> float:
            try:
                value = float(raw)
            except (TypeError, ValueError) as exc:
                raise typer.BadParameter(f"{name} must be a number.") from exc
            if min_value is not None and value < min_value:
                raise typer.BadParameter(f"{name} must be >= {min_value}.")
            if max_value is not None and value > max_value:
                raise typer.BadParameter(f"{name} must be <= {max_value}.")
            return value

        def parse_bool_option(name: str, raw: Any) -> bool:
            if isinstance(raw, bool):
                return raw
            if isinstance(raw, str):
                value = raw.strip().lower()
                if value in {"1", "true", "yes", "on"}:
                    return True
                if value in {"0", "false", "no", "off"}:
                    return False
            raise typer.BadParameter(f"{name} must be a boolean.")

        min_len = post_prep.get("sequence_length_min")
        max_len = post_prep.get("sequence_length_max")
        if min_len is not None:
            min_len = parse_int_option("post_prep.sequence_length_min", min_len)
            post_prep["sequence_length_min"] = min_len
        if max_len is not None:
            max_len = parse_int_option("post_prep.sequence_length_max", max_len)
            post_prep["sequence_length_max"] = max_len
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

        primer_max_mismatch_raw = post_prep.get("primer_max_mismatch", 0)
        primer_max_error_rate_raw = post_prep.get("primer_max_error_rate", 0.0)
        primer_min_overlap_bp_raw = post_prep.get("primer_min_overlap_bp")
        primer_min_overlap_ratio_raw = post_prep.get("primer_min_overlap_ratio", 1.0)
        primer_end_max_offset_raw = post_prep.get("primer_end_max_offset", 0)
        primer_trim_mode_raw = post_prep.get("primer_trim_mode", PRIMER_TRIM_MODE_ONE_OR_BOTH)
        primer_keep_retained_raw = post_prep.get("primer_keep_retained_fasta", True)
        primer_iter_enable_raw = post_prep.get("primer_iter_enable", False)
        primer_iter_max_rounds_raw = post_prep.get("primer_iter_max_rounds", 3)
        primer_iter_stop_delta_raw = post_prep.get("primer_iter_stop_delta", 0.002)
        primer_iter_target_conf_raw = post_prep.get("primer_iter_target_conf", 0.98)
        primer_recheck_tool_raw = post_prep.get("primer_recheck_tool", "off")
        primer_recheck_min_identity_raw = post_prep.get("primer_recheck_min_identity", 0.85)
        primer_recheck_min_query_cov_raw = post_prep.get("primer_recheck_min_query_cov", 0.7)
        primer_phylo_check_raw = post_prep.get("primer_phylo_check", "off")
        primer_phylo_target_raw = post_prep.get("primer_phylo_target_confidence", "medium")
        primer_sidecar_format_raw = post_prep.get("primer_sidecar_format", "tsv")

        primer_max_mismatch = parse_int_option("post_prep.primer_max_mismatch", primer_max_mismatch_raw, 0)
        primer_max_error_rate = parse_float_option(
            "post_prep.primer_max_error_rate", primer_max_error_rate_raw, 0.0, 1.0
        )
        primer_min_overlap_bp = None
        if primer_min_overlap_bp_raw is not None:
            primer_min_overlap_bp = parse_int_option("post_prep.primer_min_overlap_bp", primer_min_overlap_bp_raw, 1)
        primer_min_overlap_ratio = parse_float_option(
            "post_prep.primer_min_overlap_ratio", primer_min_overlap_ratio_raw, 0.0, 1.0
        )
        primer_end_max_offset = parse_int_option("post_prep.primer_end_max_offset", primer_end_max_offset_raw, 0)

        if not isinstance(primer_trim_mode_raw, str):
            raise typer.BadParameter("post_prep.primer_trim_mode must be a string.")
        primer_trim_mode = primer_trim_mode_raw.strip().lower()
        if primer_trim_mode not in PRIMER_TRIM_MODES:
            modes = ", ".join(sorted(PRIMER_TRIM_MODES))
            raise typer.BadParameter(f"post_prep.primer_trim_mode must be one of: {modes}")

        primer_keep_retained_fasta = parse_bool_option("post_prep.primer_keep_retained_fasta", primer_keep_retained_raw)
        primer_iter_enable = parse_bool_option("post_prep.primer_iter_enable", primer_iter_enable_raw)
        primer_iter_max_rounds = parse_int_option("post_prep.primer_iter_max_rounds", primer_iter_max_rounds_raw, 1)
        primer_iter_stop_delta = parse_float_option("post_prep.primer_iter_stop_delta", primer_iter_stop_delta_raw, 0.0)
        primer_iter_target_conf = parse_float_option(
            "post_prep.primer_iter_target_conf", primer_iter_target_conf_raw, 0.0, 1.0
        )

        if not isinstance(primer_recheck_tool_raw, str):
            raise typer.BadParameter("post_prep.primer_recheck_tool must be a string.")
        primer_recheck_tool = primer_recheck_tool_raw.strip().lower()
        if primer_recheck_tool not in {"off", "vsearch", "blast"}:
            raise typer.BadParameter("post_prep.primer_recheck_tool must be one of: off, vsearch, blast")
        primer_recheck_min_identity = parse_float_option(
            "post_prep.primer_recheck_min_identity", primer_recheck_min_identity_raw, 0.0, 1.0
        )
        primer_recheck_min_query_cov = parse_float_option(
            "post_prep.primer_recheck_min_query_cov", primer_recheck_min_query_cov_raw, 0.0, 1.0
        )

        if not isinstance(primer_phylo_check_raw, str):
            raise typer.BadParameter("post_prep.primer_phylo_check must be a string.")
        primer_phylo_check = primer_phylo_check_raw.strip().lower()
        if primer_phylo_check not in {"off", "flag_only"}:
            raise typer.BadParameter("post_prep.primer_phylo_check must be one of: off, flag_only")

        if not isinstance(primer_phylo_target_raw, str):
            raise typer.BadParameter("post_prep.primer_phylo_target_confidence must be a string.")
        primer_phylo_target = primer_phylo_target_raw.strip().lower()
        if primer_phylo_target not in {"low", "medium"}:
            raise typer.BadParameter("post_prep.primer_phylo_target_confidence must be one of: low, medium")

        if not isinstance(primer_sidecar_format_raw, str):
            raise typer.BadParameter("post_prep.primer_sidecar_format must be a string.")
        primer_sidecar_format = primer_sidecar_format_raw.strip().lower()
        if primer_sidecar_format not in {"tsv", "jsonl"}:
            raise typer.BadParameter("post_prep.primer_sidecar_format must be one of: tsv, jsonl")

        post_prep["primer_max_mismatch"] = primer_max_mismatch
        post_prep["primer_max_error_rate"] = primer_max_error_rate
        post_prep["primer_min_overlap_bp"] = primer_min_overlap_bp
        post_prep["primer_min_overlap_ratio"] = primer_min_overlap_ratio
        post_prep["primer_end_max_offset"] = primer_end_max_offset
        post_prep["primer_trim_mode"] = primer_trim_mode
        post_prep["primer_keep_retained_fasta"] = primer_keep_retained_fasta
        post_prep["primer_iter_enable"] = primer_iter_enable
        post_prep["primer_iter_max_rounds"] = primer_iter_max_rounds
        post_prep["primer_iter_stop_delta"] = primer_iter_stop_delta
        post_prep["primer_iter_target_conf"] = primer_iter_target_conf
        post_prep["primer_recheck_tool"] = primer_recheck_tool
        post_prep["primer_recheck_min_identity"] = primer_recheck_min_identity
        post_prep["primer_recheck_min_query_cov"] = primer_recheck_min_query_cov
        post_prep["primer_phylo_check"] = primer_phylo_check
        post_prep["primer_phylo_target_confidence"] = primer_phylo_target
        post_prep["primer_sidecar_format"] = primer_sidecar_format

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


def setup_entrez(ncbi_cfg: Dict, warn_if_missing: bool = True) -> None:
    email = ncbi_cfg.get("email") or os.environ.get("NCBI_EMAIL")
    api_key = ncbi_cfg.get("api_key") or os.environ.get("NCBI_API_KEY")

    if not email and warn_if_missing:
        console.print("[yellow]WARNING:[/yellow] NCBI email is not set. Set ncbi.email or NCBI_EMAIL.")
    Entrez.email = email or ""

    if api_key:
        Entrez.api_key = api_key


def normalize_marker_map(markers_cfg: Dict, source: BuildSource = BuildSource.NCBI) -> Dict[str, Dict]:
    marker_map: Dict[str, Dict] = {}
    uses_ncbi = source in {BuildSource.NCBI, BuildSource.BOTH}
    uses_bold = source in {BuildSource.BOLD, BuildSource.BOTH}
    for key, cfg in markers_cfg.items():
        phrases = cfg.get("phrases") or []
        terms = cfg.get("terms") or []
        region_patterns = cfg.get("region_patterns") or []
        header_format = cfg.get("header_format")
        feature_types = cfg.get("feature_types")
        feature_fields = cfg.get("feature_fields")
        aliases = cfg.get("aliases") or []
        bold_cfg = cfg.get("bold") or {}
        if not isinstance(phrases, list):
            raise typer.BadParameter(f"markers.{key}.phrases must be a list of strings.")
        if any(not isinstance(item, str) for item in phrases):
            raise typer.BadParameter(f"markers.{key}.phrases must contain only strings.")
        if not isinstance(terms, list):
            raise typer.BadParameter(f"markers.{key}.terms must be a list of strings.")
        if any(not isinstance(item, str) for item in terms):
            raise typer.BadParameter(f"markers.{key}.terms must contain only strings.")
        if not isinstance(region_patterns, list):
            raise typer.BadParameter(f"markers.{key}.region_patterns must be a list of strings.")
        if any(not isinstance(item, str) for item in region_patterns):
            raise typer.BadParameter(f"markers.{key}.region_patterns must contain only strings.")
        if not isinstance(aliases, list):
            raise typer.BadParameter(f"markers.{key}.aliases must be a list of strings.")
        if any(not isinstance(item, str) for item in aliases):
            raise typer.BadParameter(f"markers.{key}.aliases must contain only strings.")
        if feature_types is not None and not isinstance(feature_types, list):
            raise typer.BadParameter(f"markers.{key}.feature_types must be a list of strings.")
        if feature_types is not None and any(not isinstance(item, str) for item in feature_types):
            raise typer.BadParameter(f"markers.{key}.feature_types must contain only strings.")
        if feature_fields is not None and not isinstance(feature_fields, list):
            raise typer.BadParameter(f"markers.{key}.feature_fields must be a list of strings.")
        if feature_fields is not None and any(not isinstance(item, str) for item in feature_fields):
            raise typer.BadParameter(f"markers.{key}.feature_fields must contain only strings.")
        if bold_cfg and not isinstance(bold_cfg, dict):
            raise typer.BadParameter(f"markers.{key}.bold must be a table (dict).")
        marker_codes = bold_cfg.get("marker_codes") if isinstance(bold_cfg, dict) else None
        if marker_codes is None:
            marker_codes = []
        if not isinstance(marker_codes, list):
            raise typer.BadParameter(f"markers.{key}.bold.marker_codes must be a list of strings.")
        if any(not isinstance(item, str) for item in marker_codes):
            raise typer.BadParameter(f"markers.{key}.bold.marker_codes must contain only strings.")
        if uses_ncbi and not phrases and not terms:
            raise typer.BadParameter(f"markers.{key} must define phrases or terms when source uses NCBI.")
        if uses_bold and not str(key).strip():
            raise typer.BadParameter(
                f"markers.{key} must define bold.marker_codes, aliases, phrases, terms, or marker key fallback when source uses BOLD."
            )
        marker_map[key] = {
            "phrases": phrases,
            "terms": terms,
            "aliases": aliases,
            "region_patterns": region_patterns,
            "header_format": header_format,
            "feature_types": feature_types,
            "feature_fields": feature_fields,
            "bold": {"marker_codes": marker_codes},
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


@dataclass
class PrimerEndpointMatch:
    primer: str
    overlap_bp: int
    mismatches: int
    score: int
    full_len_match: bool
    trim_bp: int


@dataclass
class OrientationScore:
    name: str
    left: Optional[PrimerEndpointMatch]
    right: Optional[PrimerEndpointMatch]

    @property
    def matched_ends(self) -> int:
        return (1 if self.left else 0) + (1 if self.right else 0)

    @property
    def score_total(self) -> int:
        left_score = self.left.score if self.left else 0
        right_score = self.right.score if self.right else 0
        return left_score + right_score

    @property
    def mismatch_total(self) -> int:
        left_mm = self.left.mismatches if self.left else 0
        right_mm = self.right.mismatches if self.right else 0
        return left_mm + right_mm

    def rank(self) -> Tuple[int, int, int]:
        return (self.matched_ends, self.score_total, -self.mismatch_total)


def primer_base_matches(seq_base: str, primer_base: str) -> bool:
    values = IUPAC_DNA_VALUES.get(primer_base.upper())
    if values is None:
        return False
    return seq_base.upper() in values


def required_overlap_bp(primer_len: int, min_overlap_bp: Optional[int], min_overlap_ratio: float) -> int:
    ratio_required = int(round(primer_len * min_overlap_ratio))
    ratio_required = max(1, min(primer_len, ratio_required))
    if min_overlap_bp is None:
        return ratio_required
    return max(1, min(primer_len, max(min_overlap_bp, ratio_required)))


def count_mismatches(seq_segment: str, primer_segment: str) -> int:
    return sum(1 for s, p in zip(seq_segment.upper(), primer_segment.upper()) if not primer_base_matches(s, p))


def find_best_prefix_match(
    seq: str,
    primers: List[str],
    max_mismatch: int,
    max_error_rate: float,
    min_overlap_bp: Optional[int],
    min_overlap_ratio: float,
    max_offset: int = 0,
) -> Optional[PrimerEndpointMatch]:
    best: Optional[PrimerEndpointMatch] = None
    for primer in primers:
        max_overlap = min(len(primer), len(seq))
        required = required_overlap_bp(len(primer), min_overlap_bp, min_overlap_ratio)
        if max_overlap < required:
            continue
        max_shift = min(max_offset, max(0, len(seq) - required))
        for offset in range(max_shift + 1):
            available = len(seq) - offset
            if available < required:
                continue
            for overlap in range(min(max_overlap, available), required - 1, -1):
                primer_seg = primer[-overlap:]
                seq_seg = seq[offset : offset + overlap]
                mismatches = count_mismatches(seq_seg, primer_seg)
                if mismatches > max_mismatch:
                    continue
                error_rate = mismatches / overlap if overlap else 1.0
                if error_rate > max_error_rate:
                    continue
                candidate = PrimerEndpointMatch(
                    primer=primer,
                    overlap_bp=overlap,
                    mismatches=mismatches,
                    score=overlap - mismatches,
                    full_len_match=(overlap == len(primer)),
                    trim_bp=offset + overlap,
                )
                if best is None:
                    best = candidate
                    continue
                current_key = (
                    candidate.overlap_bp,
                    candidate.score,
                    -candidate.mismatches,
                    int(candidate.full_len_match),
                    -offset,
                )
                best_key = (
                    best.overlap_bp,
                    best.score,
                    -best.mismatches,
                    int(best.full_len_match),
                    -(best.trim_bp - best.overlap_bp),
                )
                if current_key > best_key:
                    best = candidate
    return best


def find_best_suffix_match(
    seq: str,
    primers: List[str],
    max_mismatch: int,
    max_error_rate: float,
    min_overlap_bp: Optional[int],
    min_overlap_ratio: float,
    max_offset: int = 0,
) -> Optional[PrimerEndpointMatch]:
    best: Optional[PrimerEndpointMatch] = None
    for primer in primers:
        max_overlap = min(len(primer), len(seq))
        required = required_overlap_bp(len(primer), min_overlap_bp, min_overlap_ratio)
        if max_overlap < required:
            continue
        max_shift = min(max_offset, max(0, len(seq) - required))
        for shift in range(max_shift + 1):
            end = len(seq) - shift
            if end < required:
                continue
            for overlap in range(min(max_overlap, end), required - 1, -1):
                start = end - overlap
                primer_seg = primer[:overlap]
                seq_seg = seq[start:end]
                mismatches = count_mismatches(seq_seg, primer_seg)
                if mismatches > max_mismatch:
                    continue
                error_rate = mismatches / overlap if overlap else 1.0
                if error_rate > max_error_rate:
                    continue
                candidate = PrimerEndpointMatch(
                    primer=primer,
                    overlap_bp=overlap,
                    mismatches=mismatches,
                    score=overlap - mismatches,
                    full_len_match=(overlap == len(primer)),
                    trim_bp=len(seq) - start,
                )
                if best is None:
                    best = candidate
                    continue
                current_key = (
                    candidate.overlap_bp,
                    candidate.score,
                    -candidate.mismatches,
                    int(candidate.full_len_match),
                    -shift,
                )
                best_shift = max(0, best.trim_bp - best.overlap_bp)
                best_key = (
                    best.overlap_bp,
                    best.score,
                    -best.mismatches,
                    int(best.full_len_match),
                    -best_shift,
                )
                if current_key > best_key:
                    best = candidate
    return best


def resolve_orientation(seq: str, canonical: OrientationScore, reverse: OrientationScore) -> Tuple[OrientationScore, bool]:
    can_rank = canonical.rank()
    rev_rank = reverse.rank()
    if can_rank > rev_rank:
        return canonical, False
    if rev_rank > can_rank:
        return reverse, False
    return canonical, True


def confidence_label(matched_ends: int, mismatch_total: int, ambiguous_orientation: bool) -> str:
    if matched_ends == 2 and mismatch_total <= 1 and not ambiguous_orientation:
        return "high"
    if matched_ends >= 1:
        return "medium"
    return "low"


def compute_trim_lengths_from_row(row: Dict[str, Any], trim_mode: str) -> Tuple[int, int]:
    left_hit = int(row.get("left_hit", 0)) > 0
    right_hit = int(row.get("right_hit", 0)) > 0
    left_trim_bp = int(row.get("left_trim_bp", 0) or 0)
    right_trim_bp = int(row.get("right_trim_bp", 0) or 0)
    left_overlap_bp = int(row.get("left_overlap_bp", 0) or 0)
    right_overlap_bp = int(row.get("right_overlap_bp", 0) or 0)
    left_overlap = left_trim_bp if left_trim_bp > 0 else left_overlap_bp
    right_overlap = right_trim_bp if right_trim_bp > 0 else right_overlap_bp
    if trim_mode == PRIMER_TRIM_MODE_ONE_OR_BOTH:
        return (left_overlap if left_hit else 0, right_overlap if right_hit else 0)
    if trim_mode == PRIMER_TRIM_MODE_BOTH_REQUIRED:
        if left_hit and right_hit:
            return left_overlap, right_overlap
        return 0, 0
    return 0, 0


def summarize_rows_to_records(
    rows: List[Dict[str, Any]],
    seq_by_header: Dict[str, str],
    trim_mode: str,
) -> Tuple[List[Tuple[str, str]], Dict[str, Any]]:
    trimmed_records: List[Tuple[str, str]] = []
    trimmed_both = 0
    trimmed_left_only = 0
    trimmed_right_only = 0
    untrimmed = 0
    dropped_empty = 0
    canonical_orientation = 0
    reverse_orientation = 0
    confidence_high = 0
    confidence_medium = 0
    confidence_low = 0

    for row in rows:
        header = str(row.get("record_id", ""))
        seq = seq_by_header.get(header, "")
        if not seq:
            continue

        confidence = str(row.get("confidence", "low"))
        if confidence == "high":
            confidence_high += 1
        elif confidence == "medium":
            confidence_medium += 1
        else:
            confidence_low += 1

        matched_ends = int(row.get("left_hit", 0)) + int(row.get("right_hit", 0))
        if matched_ends > 0:
            if str(row.get("orientation_chosen")) == "canonical":
                canonical_orientation += 1
            else:
                reverse_orientation += 1

        left_len, right_len = compute_trim_lengths_from_row(row, trim_mode)
        ends_trimmed = (1 if left_len else 0) + (1 if right_len else 0)
        if ends_trimmed == 0:
            untrimmed += 1
        elif ends_trimmed == 2:
            trimmed_both += 1
        elif left_len:
            trimmed_left_only += 1
        else:
            trimmed_right_only += 1

        right_index = len(seq) - right_len if right_len else len(seq)
        trimmed_seq = seq[left_len:right_index]
        row["trim_start"] = left_len
        row["trim_end"] = right_index
        if not trimmed_seq:
            dropped_empty += 1
            row["dropped_empty"] = 1
            continue
        row["dropped_empty"] = 0
        trimmed_records.append((header, trimmed_seq))

    before = len(rows)
    after = len(trimmed_records)
    summary = {
        "before": before,
        "after": after,
        "removed": before - after,
        "trimmed_both": trimmed_both,
        "trimmed_left_only": trimmed_left_only,
        "trimmed_right_only": trimmed_right_only,
        "untrimmed": untrimmed,
        "dropped_empty": dropped_empty,
        "canonical_orientation": canonical_orientation,
        "reverse_orientation": reverse_orientation,
        "confidence_high": confidence_high,
        "confidence_medium": confidence_medium,
        "confidence_low": confidence_low,
        "high_conf_rate": (confidence_high / before) if before else 0.0,
    }
    return trimmed_records, summary


def run_vsearch_endpoint_recheck(
    rows: List[Dict[str, Any]],
    seq_by_header: Dict[str, str],
    forward_primers: List[str],
    reverse_primers: List[str],
    forward_rc: List[str],
    reverse_rc: List[str],
    min_identity: float,
    min_query_cov: float,
) -> Tuple[int, int, Optional[str]]:
    vsearch_bin = shutil.which("vsearch")
    if not vsearch_bin:
        return 0, 0, "vsearch_not_found"

    candidates: List[Dict[str, Any]] = []
    for idx, row in enumerate(rows):
        if str(row.get("confidence", "")) not in {"low", "medium"}:
            row["recheck_status"] = "not_target_confidence"
            continue
        if int(row.get("left_hit", 0)) and int(row.get("right_hit", 0)):
            row["recheck_status"] = "already_both_ends"
            continue
        header = str(row.get("record_id", ""))
        seq = seq_by_header.get(header, "")
        if not seq:
            row["recheck_status"] = "sequence_not_found"
            continue
        candidates.append({"idx": idx, "header": header, "seq": seq, "row": row})

    if not candidates:
        return 0, 0, None

    primer_db_entries: List[Tuple[str, str]] = []
    for i, primer in enumerate(forward_primers):
        primer_db_entries.append((f"CL_{i}", primer))
    for i, primer in enumerate(reverse_rc):
        primer_db_entries.append((f"CR_{i}", primer))
    for i, primer in enumerate(reverse_primers):
        primer_db_entries.append((f"RL_{i}", primer))
    for i, primer in enumerate(forward_rc):
        primer_db_entries.append((f"RR_{i}", primer))

    max_primer_len = max((len(p) for _, p in primer_db_entries), default=30)
    window_len = max_primer_len + 8
    if window_len < 20:
        window_len = 20

    rescued = 0
    attempted = 0
    with tempfile.TemporaryDirectory(prefix="taxondb-vsearch-") as tmpdir:
        tmp = Path(tmpdir)
        db_fa = tmp / "primers.fa"
        q_fa = tmp / "queries.fa"
        out_path = tmp / "hits.tsv"

        with db_fa.open("w", encoding="utf-8") as f:
            for pid, seq in primer_db_entries:
                f.write(f">{pid}\n{seq}\n")

        query_rows: List[Tuple[int, str, str]] = []
        with q_fa.open("w", encoding="utf-8") as f:
            for c in candidates:
                row = c["row"]
                seq = c["seq"]
                if not int(row.get("left_hit", 0)):
                    qid = f"{c['idx']}|L"
                    segment = seq[:window_len]
                    f.write(f">{qid}\n{segment}\n")
                    query_rows.append((c["idx"], "L", qid))
                    attempted += 1
                if not int(row.get("right_hit", 0)):
                    qid = f"{c['idx']}|R"
                    segment = seq[-window_len:]
                    f.write(f">{qid}\n{segment}\n")
                    query_rows.append((c["idx"], "R", qid))
                    attempted += 1

        cmd = [
            vsearch_bin,
            "--usearch_global",
            str(q_fa),
            "--db",
            str(db_fa),
            "--id",
            f"{min_identity:.4f}",
            "--strand",
            "plus",
            "--blast6out",
            str(out_path),
            "--maxaccepts",
            "16",
            "--maxrejects",
            "64",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            return attempted, rescued, "vsearch_failed"
        if not out_path.exists():
            return attempted, rescued, None

        primer_len_map = {pid: len(seq) for pid, seq in primer_db_entries}
        best_hits: Dict[str, Dict[str, Any]] = {}
        with out_path.open("r", encoding="utf-8") as f:
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) < 12:
                    continue
                qid, sid = cols[0], cols[1]
                try:
                    pident = float(cols[2]) / 100.0
                    aln_len = int(cols[3])
                    mismatch = int(cols[4])
                except ValueError:
                    continue
                primer_len = primer_len_map.get(sid, 0)
                if primer_len <= 0:
                    continue
                cov = aln_len / primer_len
                if pident < min_identity or cov < min_query_cov:
                    continue
                hit = {
                    "sid": sid,
                    "pident": pident,
                    "aln_len": aln_len,
                    "mismatch": mismatch,
                    "cov": cov,
                }
                prev = best_hits.get(qid)
                if prev is None:
                    best_hits[qid] = hit
                    continue
                current_key = (hit["pident"], hit["cov"], hit["aln_len"], -hit["mismatch"])
                prev_key = (prev["pident"], prev["cov"], prev["aln_len"], -prev["mismatch"])
                if current_key > prev_key:
                    best_hits[qid] = hit

        hit_applied: Dict[Tuple[int, str], bool] = {}
        for idx, side, qid in query_rows:
            hit = best_hits.get(qid)
            row = rows[idx]
            if hit is None:
                continue
            sid = hit["sid"]
            is_left_subject = sid.startswith("CL_") or sid.startswith("RL_")
            is_right_subject = sid.startswith("CR_") or sid.startswith("RR_")
            if side == "L" and not is_left_subject:
                continue
            if side == "R" and not is_right_subject:
                continue

            row["recheck_status"] = "rescued_by_vsearch"
            if side == "L":
                row["left_hit"] = 1
                row["left_overlap_bp"] = int(hit["aln_len"])
                row["left_trim_bp"] = int(hit["aln_len"])
                row["left_mismatch"] = int(hit["mismatch"])
            else:
                row["right_hit"] = 1
                row["right_overlap_bp"] = int(hit["aln_len"])
                row["right_trim_bp"] = int(hit["aln_len"])
                row["right_mismatch"] = int(hit["mismatch"])
            mism_total = int(row.get("left_mismatch", 0)) + int(row.get("right_mismatch", 0))
            matched_ends = int(row.get("left_hit", 0)) + int(row.get("right_hit", 0))
            row["confidence"] = confidence_label(matched_ends, mism_total, bool(int(row.get("orientation_ambiguous", 0))))
            rescued += 1
            hit_applied[(idx, side)] = True

        for idx, side, _ in query_rows:
            if hit_applied.get((idx, side)):
                continue
            row = rows[idx]
            if str(row.get("recheck_status", "")) != "rescued_by_vsearch":
                row["recheck_status"] = "attempted_no_hit"

    return attempted, rescued, None


def apply_post_prep_primer_trim(
    fasta_path: Path,
    forward_primers: List[str],
    reverse_primers: List[str],
    options: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    opts = options or {}
    trim_mode = str(opts.get("trim_mode", PRIMER_TRIM_MODE_ONE_OR_BOTH)).strip().lower()
    max_mismatch = int(opts.get("max_mismatch", 0))
    max_error_rate = float(opts.get("max_error_rate", 0.0))
    min_overlap_bp = opts.get("min_overlap_bp")
    if min_overlap_bp is not None:
        min_overlap_bp = int(min_overlap_bp)
    min_overlap_ratio = float(opts.get("min_overlap_ratio", 1.0))
    end_max_offset = int(opts.get("end_max_offset", 0))
    keep_retained_fasta = bool(opts.get("keep_retained_fasta", True))
    iter_enable = bool(opts.get("iter_enable", False))
    iter_max_rounds = int(opts.get("iter_max_rounds", 3))
    iter_stop_delta = float(opts.get("iter_stop_delta", 0.002))
    iter_target_conf = float(opts.get("iter_target_conf", 0.98))
    sidecar_format = str(opts.get("sidecar_format", "tsv")).strip().lower()
    recheck_tool = str(opts.get("recheck_tool", "off")).strip().lower()
    recheck_min_identity = float(opts.get("recheck_min_identity", 0.85))
    recheck_min_query_cov = float(opts.get("recheck_min_query_cov", 0.7))
    phylo_target_confidence = str(opts.get("phylo_target_confidence", "medium")).strip().lower()

    forward_rc = [str(Seq(p).reverse_complement()).upper() for p in forward_primers]
    reverse_rc = [str(Seq(p).reverse_complement()).upper() for p in reverse_primers]

    retained_path = fasta_path.with_suffix(fasta_path.suffix + ".postprep.primer.retained.fasta")
    if keep_retained_fasta:
        shutil.copyfile(fasta_path, retained_path)

    records: List[Tuple[str, str]] = []
    with fasta_path.open("r", encoding="utf-8") as in_f:
        for record in SeqIO.parse(in_f, "fasta"):
            records.append((record.description, str(record.seq).upper().replace("U", "T")))

    before_count = len(records)
    if before_count == 0:
        return {
            "before": 0,
            "after": 0,
            "removed": 0,
            "trimmed_both": 0,
            "trimmed_left_only": 0,
            "trimmed_right_only": 0,
            "untrimmed": 0,
            "dropped_empty": 0,
            "canonical_orientation": 0,
            "reverse_orientation": 0,
            "confidence_high": 0,
            "confidence_medium": 0,
            "confidence_low": 0,
            "rounds_run": 0,
            "best_round": 0,
            "high_conf_rate": 0.0,
            "sidecar_path": None,
            "retained_path": str(retained_path) if keep_retained_fasta else None,
            "recheck_tool": recheck_tool,
            "phylo_target_confidence": phylo_target_confidence,
        }

    round_limit = iter_max_rounds if iter_enable else 1
    round_limit = max(1, round_limit)
    all_round_rows: List[Dict[str, Any]] = []
    round_summaries: List[Dict[str, Any]] = []
    prev_high_rate: Optional[float] = None

    for round_idx in range(1, round_limit + 1):
        # Phase 2: deterministic threshold relaxation to make iterative rescue measurable.
        relax = round_idx - 1
        round_max_mismatch = max_mismatch + relax
        round_overlap_ratio = max(0.5, min_overlap_ratio - (0.05 * relax))
        round_overlap_bp = min_overlap_bp
        if round_overlap_bp is not None:
            round_overlap_bp = max(8, round_overlap_bp - relax)

        trimmed_records: List[Tuple[str, str]] = []
        trimmed_both = 0
        trimmed_left_only = 0
        trimmed_right_only = 0
        untrimmed = 0
        dropped_empty = 0
        canonical_orientation = 0
        reverse_orientation = 0
        confidence_high = 0
        confidence_medium = 0
        confidence_low = 0

        for header, seq in records:
            can = OrientationScore(
                name="canonical",
                left=find_best_prefix_match(
                    seq,
                    forward_primers,
                    round_max_mismatch,
                    max_error_rate,
                    round_overlap_bp,
                    round_overlap_ratio,
                    max_offset=end_max_offset,
                ),
                right=find_best_suffix_match(
                    seq,
                    reverse_rc,
                    round_max_mismatch,
                    max_error_rate,
                    round_overlap_bp,
                    round_overlap_ratio,
                    max_offset=end_max_offset,
                ),
            )
            rev = OrientationScore(
                name="reverse",
                left=find_best_prefix_match(
                    seq,
                    reverse_primers,
                    round_max_mismatch,
                    max_error_rate,
                    round_overlap_bp,
                    round_overlap_ratio,
                    max_offset=end_max_offset,
                ),
                right=find_best_suffix_match(
                    seq,
                    forward_rc,
                    round_max_mismatch,
                    max_error_rate,
                    round_overlap_bp,
                    round_overlap_ratio,
                    max_offset=end_max_offset,
                ),
            )
            chosen, ambiguous_orientation = resolve_orientation(seq, can, rev)
            confidence = confidence_label(chosen.matched_ends, chosen.mismatch_total, ambiguous_orientation)
            if confidence == "high":
                confidence_high += 1
            elif confidence == "medium":
                confidence_medium += 1
            else:
                confidence_low += 1

            if chosen.matched_ends > 0:
                if chosen.name == "canonical":
                    canonical_orientation += 1
                else:
                    reverse_orientation += 1

            do_left_trim = False
            do_right_trim = False
            if trim_mode == PRIMER_TRIM_MODE_ONE_OR_BOTH:
                do_left_trim = chosen.left is not None
                do_right_trim = chosen.right is not None
            elif trim_mode == PRIMER_TRIM_MODE_BOTH_REQUIRED:
                if chosen.left is not None and chosen.right is not None:
                    do_left_trim = True
                    do_right_trim = True
            elif trim_mode == PRIMER_TRIM_MODE_ONE_END_MARK_ONLY:
                do_left_trim = False
                do_right_trim = False

            left_len = chosen.left.trim_bp if (do_left_trim and chosen.left) else 0
            right_len = chosen.right.trim_bp if (do_right_trim and chosen.right) else 0
            ends_trimmed = (1 if left_len else 0) + (1 if right_len else 0)
            if ends_trimmed == 0:
                untrimmed += 1
            elif ends_trimmed == 2:
                trimmed_both += 1
            elif left_len:
                trimmed_left_only += 1
            else:
                trimmed_right_only += 1

            right_index = len(seq) - right_len if right_len else len(seq)
            trimmed_seq = seq[left_len:right_index]
            if not trimmed_seq:
                dropped_empty += 1
                all_round_rows.append(
                    {
                        "record_id": header,
                        "round": round_idx,
                        "orientation_chosen": chosen.name,
                        "orientation_ambiguous": int(ambiguous_orientation),
                        "left_hit": int(chosen.left is not None),
                        "right_hit": int(chosen.right is not None),
                        "left_overlap_bp": chosen.left.overlap_bp if chosen.left else 0,
                        "right_overlap_bp": chosen.right.overlap_bp if chosen.right else 0,
                        "left_trim_bp": left_len,
                        "right_trim_bp": right_len,
                        "left_mismatch": chosen.left.mismatches if chosen.left else 0,
                        "right_mismatch": chosen.right.mismatches if chosen.right else 0,
                        "left_primer_name": chosen.left.primer if chosen.left else "",
                        "right_primer_name": chosen.right.primer if chosen.right else "",
                        "trim_start": left_len,
                        "trim_end": right_index,
                        "trim_mode": trim_mode,
                        "confidence": confidence,
                        "dropped_empty": 1,
                        "recheck_tool": recheck_tool,
                        "recheck_status": "not_run_phase2",
                        "phylo_mismatch_flag": 0,
                    }
                )
                continue

            trimmed_records.append((header, trimmed_seq))
            all_round_rows.append(
                {
                    "record_id": header,
                    "round": round_idx,
                    "orientation_chosen": chosen.name,
                    "orientation_ambiguous": int(ambiguous_orientation),
                    "left_hit": int(chosen.left is not None),
                    "right_hit": int(chosen.right is not None),
                    "left_overlap_bp": chosen.left.overlap_bp if chosen.left else 0,
                    "right_overlap_bp": chosen.right.overlap_bp if chosen.right else 0,
                    "left_trim_bp": left_len,
                    "right_trim_bp": right_len,
                    "left_mismatch": chosen.left.mismatches if chosen.left else 0,
                    "right_mismatch": chosen.right.mismatches if chosen.right else 0,
                    "left_primer_name": chosen.left.primer if chosen.left else "",
                    "right_primer_name": chosen.right.primer if chosen.right else "",
                    "trim_start": left_len,
                    "trim_end": right_index,
                    "trim_mode": trim_mode,
                    "confidence": confidence,
                    "dropped_empty": 0,
                    "recheck_tool": recheck_tool,
                    "recheck_status": "not_run_phase2",
                    "phylo_mismatch_flag": 0,
                }
            )

        after_count = len(trimmed_records)
        high_rate = confidence_high / before_count if before_count else 0.0
        round_summaries.append(
            {
                "round": round_idx,
                "records": trimmed_records,
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
                "confidence_high": confidence_high,
                "confidence_medium": confidence_medium,
                "confidence_low": confidence_low,
                "high_conf_rate": high_rate,
                "round_max_mismatch": round_max_mismatch,
                "round_min_overlap_ratio": round_overlap_ratio,
                "round_min_overlap_bp": round_overlap_bp or 0,
            }
        )

        if not iter_enable:
            break
        if high_rate >= iter_target_conf:
            break
        if prev_high_rate is not None and (high_rate - prev_high_rate) < iter_stop_delta:
            break
        prev_high_rate = high_rate

    best_summary = max(
        round_summaries,
        key=lambda row: (
            row["high_conf_rate"],
            row["after"],
            -row["dropped_empty"],
            row["round"],
        ),
    )
    best_round = int(best_summary["round"])
    seq_by_header = {header: seq for header, seq in records}
    sidecar_rows = [dict(row) for row in all_round_rows if int(row["round"]) == best_round]

    recheck_attempted = 0
    recheck_rescued = 0
    recheck_error: Optional[str] = None
    if recheck_tool == "vsearch":
        recheck_attempted, recheck_rescued, recheck_error = run_vsearch_endpoint_recheck(
            sidecar_rows,
            seq_by_header,
            forward_primers,
            reverse_primers,
            forward_rc,
            reverse_rc,
            min_identity=recheck_min_identity,
            min_query_cov=recheck_min_query_cov,
        )
        rebuilt_records, rebuilt_summary = summarize_rows_to_records(sidecar_rows, seq_by_header, trim_mode)
        best_summary.update(rebuilt_summary)
        best_summary["records"] = rebuilt_records
    elif recheck_tool != "off":
        recheck_error = f"unsupported_tool:{recheck_tool}"

    tmp_path = fasta_path.with_suffix(fasta_path.suffix + ".postprep.primer.tmp")
    try:
        with tmp_path.open("w", encoding="utf-8") as out_f:
            for header, seq in best_summary["records"]:
                out_f.write(f">{header}\n")
                out_f.write(f"{seq}\n")
        tmp_path.replace(fasta_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()

    sidecar_path = None
    if sidecar_format == "jsonl":
        sidecar_path = fasta_path.with_suffix(fasta_path.suffix + ".postprep.primer.jsonl")
        with sidecar_path.open("w", encoding="utf-8") as f:
            for row in sidecar_rows:
                f.write(json.dumps(row, ensure_ascii=False) + "\n")
    else:
        sidecar_path = fasta_path.with_suffix(fasta_path.suffix + ".postprep.primer.tsv")
        fields = [
            "record_id",
            "round",
            "orientation_chosen",
            "orientation_ambiguous",
            "left_hit",
            "right_hit",
            "left_overlap_bp",
            "right_overlap_bp",
            "left_trim_bp",
            "right_trim_bp",
            "left_mismatch",
            "right_mismatch",
            "left_primer_name",
            "right_primer_name",
            "trim_start",
            "trim_end",
            "trim_mode",
            "confidence",
            "dropped_empty",
            "recheck_tool",
            "recheck_status",
            "phylo_mismatch_flag",
        ]
        with sidecar_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(sidecar_rows)

    return {
        "before": int(best_summary["before"]),
        "after": int(best_summary["after"]),
        "removed": int(best_summary["removed"]),
        "trimmed_both": int(best_summary["trimmed_both"]),
        "trimmed_left_only": int(best_summary["trimmed_left_only"]),
        "trimmed_right_only": int(best_summary["trimmed_right_only"]),
        "untrimmed": int(best_summary["untrimmed"]),
        "dropped_empty": int(best_summary["dropped_empty"]),
        "canonical_orientation": int(best_summary["canonical_orientation"]),
        "reverse_orientation": int(best_summary["reverse_orientation"]),
        "confidence_high": int(best_summary["confidence_high"]),
        "confidence_medium": int(best_summary["confidence_medium"]),
        "confidence_low": int(best_summary["confidence_low"]),
        "high_conf_rate": float(best_summary["high_conf_rate"]),
        "rounds_run": len(round_summaries),
        "best_round": best_round,
        "sidecar_path": str(sidecar_path) if sidecar_path else None,
        "retained_path": str(retained_path) if keep_retained_fasta else None,
        "recheck_tool": recheck_tool,
        "recheck_attempted": recheck_attempted,
        "recheck_rescued": recheck_rescued,
        "recheck_error": recheck_error,
        "phylo_target_confidence": phylo_target_confidence,
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
                    "source": row.get("source", ""),
                    "source_record_id": row.get("source_record_id", ""),
                    "processid": row.get("processid", ""),
                    "sampleid": row.get("sampleid", ""),
                    "marker_key": row.get("marker_key", ""),
                    "linked_to_ncbi": row.get("linked_to_ncbi", ""),
                    "emitted_to_fasta": row.get("emitted_to_fasta", ""),
                    "skip_reason": row.get("skip_reason", ""),
                }
            )

    mapped_rows.sort(key=lambda row: (row["acc_id"], row["accession"], row["organism_name"], row["header"]))

    with mapping_path.open("w", newline="", encoding="utf-8") as out_f:
        writer = csv.DictWriter(
            out_f,
            fieldnames=[
                "acc_id",
                "accession",
                "organism_name",
                "header",
                "source",
                "source_record_id",
                "processid",
                "sampleid",
                "marker_key",
                "linked_to_ncbi",
                "emitted_to_fasta",
                "skip_reason",
            ],
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


def fetch_taxonomy_scientific_name(taxid: str) -> str:
    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    record = Entrez.read(handle)

    if isinstance(record, list) and record:
        scientific_name = record[0].get("ScientificName")
        if scientific_name:
            return str(scientific_name)
    if isinstance(record, dict):
        scientific_name = record.get("ScientificName")
        if scientific_name:
            return str(scientific_name)

    raise typer.BadParameter(f"Could not resolve scientific name for taxid {taxid}.")


def resolve_taxon(taxon: str, require_scientific_name: bool = False) -> ResolvedTaxon:
    warning = None
    if re.fullmatch(r"\d+", taxon):
        taxid = taxon
        scientific_name = fetch_taxonomy_scientific_name(taxid) if require_scientific_name else taxon
    else:
        term = f'"{taxon}"[Scientific Name]'
        handle = Entrez.esearch(db="taxonomy", term=term, retmax=5)
        record = Entrez.read(handle)
        ids = record.get("IdList", [])
        if not ids:
            raise typer.BadParameter(f"Taxon not found in NCBI Taxonomy: {taxon}")

        taxid = ids[0]
        if len(ids) > 1:
            warning = f"Multiple taxids found for '{taxon}'. Using {taxid}."
        scientific_name = fetch_taxonomy_scientific_name(taxid) if require_scientific_name else taxon
    return ResolvedTaxon(
        input_value=taxon,
        taxid=taxid,
        scientific_name=scientific_name,
        warning=warning,
    )


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
    fetch_retries = max(1, int(ncbi_cfg.get("fetch_retries", 4)))
    retry_delay_sec = float(ncbi_cfg.get("retry_delay_sec", max(delay_sec, 0.5)))

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

    def read_efetch_with_retries(start: int, **fetch_kwargs: Any) -> str:
        retryable_http_codes = {
            HTTPStatus.REQUEST_TIMEOUT,
            HTTPStatus.TOO_MANY_REQUESTS,
            HTTPStatus.BAD_GATEWAY,
            HTTPStatus.SERVICE_UNAVAILABLE,
            HTTPStatus.GATEWAY_TIMEOUT,
        }
        last_exc: Optional[Exception] = None
        for attempt in range(1, fetch_retries + 1):
            try:
                fetch_handle = Entrez.efetch(**fetch_kwargs)
                try:
                    return fetch_handle.read()
                finally:
                    fetch_handle.close()
            except HTTPError as exc:
                last_exc = exc
                if exc.code == HTTPStatus.BAD_REQUEST:
                    raise
                if exc.code not in retryable_http_codes or attempt >= fetch_retries:
                    raise
            except (HTTPException, RemoteDisconnected, URLError, TimeoutError, OSError) as exc:
                last_exc = exc
                if attempt >= fetch_retries:
                    raise

            wait_sec = retry_delay_sec * attempt
            print(
                f"# efetch retry start={start} attempt={attempt}/{fetch_retries}"
                f" wait={wait_sec:.2f}s reason={last_exc}",
                file=sys.stderr,
                flush=True,
            )
            time.sleep(wait_sec)

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
        return read_efetch_with_retries(
            start,
            db=db,
            rettype=rettype,
            retmode=retmode,
            id=",".join(ids),
        )

    def gen() -> Iterable[Tuple[int, str]]:
        if use_history and webenv and query_key:
            for start in range(start_at, count, per_query):
                cached = load_cached(start) if resume else None
                if cached is not None:
                    yield start, cached
                    continue
                try:
                    data = read_efetch_with_retries(
                        start,
                        db=db,
                        rettype=rettype,
                        retmode=retmode,
                        retstart=start,
                        retmax=per_query,
                        webenv=webenv,
                        query_key=query_key,
                    )
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


class TeeStream:
    def __init__(self, *streams) -> None:
        self._streams = streams

    def write(self, data: str) -> int:
        for stream in self._streams:
            stream.write(data)
        return len(data)

    def flush(self) -> None:
        for stream in self._streams:
            stream.flush()

    def isatty(self) -> bool:
        return bool(getattr(self._streams[0], "isatty", lambda: False)())

    @property
    def encoding(self) -> str:
        return getattr(self._streams[0], "encoding", "utf-8")


@contextmanager
def tee_console_output(log_path: Path):
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    with log_path.open("a", encoding="utf-8") as log_file:
        tee_stdout = TeeStream(original_stdout, log_file)
        tee_stderr = TeeStream(original_stderr, log_file)
        with redirect_stdout(tee_stdout), redirect_stderr(tee_stderr):
            yield


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
    source: BuildSource = typer.Option(
        BuildSource.NCBI,
        "--source",
        help="Sequence source: ncbi, bold, or both.",
    ),
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
      taxondbbuilder.py build -c configs/db.toml -t 117570 -m coi --source ncbi
    """
    cfg = load_config(config, source=source)
    ncbi_cfg = cfg.get("ncbi", {})
    filters_cfg = cfg.get("filters", {})
    output_cfg = cfg.get("output", {})
    post_prep_cfg = cfg.get("post_prep") or {}
    taxon_noexp = bool(cfg.get("taxon", {}).get("noexp", False))
    uses_ncbi = source in {BuildSource.NCBI, BuildSource.BOTH}
    uses_bold = source in {BuildSource.BOLD, BuildSource.BOTH}

    needs_entrez = uses_ncbi or uses_bold
    if needs_entrez:
        setup_entrez(ncbi_cfg if isinstance(ncbi_cfg, dict) else {}, warn_if_missing=uses_ncbi)
    marker_map = normalize_marker_map(cfg.get("markers", {}), source=source)

    marker_keys = [resolve_marker_key(m, marker_map) for m in marker]
    marker_query = build_marker_query(marker_keys, marker_map) if uses_ncbi else ""
    output_prefix = output_prefix.strip()
    selected_header_formats: List[str] = []
    marker_rules = []
    for key in marker_keys:
        cfg_m = marker_map[key]
        header_format = resolve_header_format(cfg_m, output_cfg)
        selected_header_formats.append(header_format)
        if not uses_ncbi:
            continue

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
                "header_format": header_format,
            }
        )

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
        post_primer_trim_options: Dict[str, Any] = {
            "trim_mode": post_prep_cfg.get("primer_trim_mode", PRIMER_TRIM_MODE_ONE_OR_BOTH),
            "max_mismatch": int(post_prep_cfg.get("primer_max_mismatch", 0)),
            "max_error_rate": float(post_prep_cfg.get("primer_max_error_rate", 0.0)),
            "min_overlap_bp": post_prep_cfg.get("primer_min_overlap_bp"),
            "min_overlap_ratio": float(post_prep_cfg.get("primer_min_overlap_ratio", 1.0)),
            "end_max_offset": int(post_prep_cfg.get("primer_end_max_offset", 0)),
            "keep_retained_fasta": bool(post_prep_cfg.get("primer_keep_retained_fasta", True)),
            "iter_enable": bool(post_prep_cfg.get("primer_iter_enable", False)),
            "iter_max_rounds": int(post_prep_cfg.get("primer_iter_max_rounds", 3)),
            "iter_stop_delta": float(post_prep_cfg.get("primer_iter_stop_delta", 0.002)),
            "iter_target_conf": float(post_prep_cfg.get("primer_iter_target_conf", 0.98)),
            "sidecar_format": post_prep_cfg.get("primer_sidecar_format", "tsv"),
            "recheck_tool": post_prep_cfg.get("primer_recheck_tool", "off"),
            "recheck_min_identity": float(post_prep_cfg.get("primer_recheck_min_identity", 0.85)),
            "recheck_min_query_cov": float(post_prep_cfg.get("primer_recheck_min_query_cov", 0.7)),
            "phylo_target_confidence": post_prep_cfg.get("primer_phylo_target_confidence", "medium"),
        }

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
            if source != BuildSource.BOTH:
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
        post_primer_trim_options = {
            "trim_mode": PRIMER_TRIM_MODE_ONE_OR_BOTH,
            "max_mismatch": 0,
            "max_error_rate": 0.0,
            "min_overlap_bp": None,
            "min_overlap_ratio": 1.0,
            "end_max_offset": 0,
            "keep_retained_fasta": True,
            "iter_enable": False,
            "iter_max_rounds": 1,
            "iter_stop_delta": 0.002,
            "iter_target_conf": 0.98,
            "sidecar_format": "tsv",
            "recheck_tool": "off",
            "recheck_min_identity": 0.85,
            "recheck_min_query_cov": 0.7,
            "phylo_target_confidence": "medium",
        }

    resolved_taxa: List[ResolvedTaxon] = []
    warnings: List[str] = []
    for t in taxon:
        resolved = resolve_taxon(t, require_scientific_name=uses_bold)
        resolved_taxa.append(resolved)
        if resolved.warning:
            warnings.append(resolved.warning)

    taxids = [item.taxid for item in resolved_taxa]

    out_path = build_output_path(out, taxids, marker_keys, output_prefix=output_prefix)
    log_path = out_path.with_suffix(out_path.suffix + ".log")
    if source == BuildSource.BOLD and (dump_gb or from_gb or resume):
        raise typer.BadParameter("--dump-gb, --from-gb, and --resume are not supported with --source bold.")
    if resume and not dump_gb and not from_gb:
        raise typer.BadParameter("--resume requires --dump-gb or --from-gb.")
    if from_gb and not from_gb.exists():
        raise typer.BadParameter(f"--from-gb not found: {from_gb}")
    if dump_gb:
        dump_gb.mkdir(parents=True, exist_ok=True)

    run_logger = setup_run_logger(log_path)
    try:
        with tee_console_output(log_path):
            print_header()
            render_run_table(
                config,
                source,
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
                if uses_ncbi:
                    for taxid in taxids:
                        query = build_query(taxid, marker_query, filters_cfg, taxon_noexp)
                        console.print(query)
                if uses_bold:
                    for item in resolved_taxa:
                        console.print(f"BOLD query taxon: {item.scientific_name}")
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
            source_merge_rows: List[Dict[str, str]] = []

            progress = Progress(
                SpinnerColumn(),
                TextColumn("[bold]{task.description}"),
                BarColumn(),
                MofNCompleteColumn(),
                TimeElapsedColumn(),
                console=console,
                disable=not console.is_terminal,
            )

            run_logger.info(f"# started: {datetime.now().isoformat()}")
            run_logger.info(f"# config: {config}")
            run_logger.info(f"# source: {source.value}")
            run_logger.info(f"# taxon input: {taxon}")
            run_logger.info(f"# taxids: {taxids}")
            run_logger.info(f"# scientific_names: {[item.scientific_name for item in resolved_taxa]}")
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
                run_logger.info(f"# post_prep.primer_trim_mode: {post_primer_trim_options['trim_mode']}")
                run_logger.info(f"# post_prep.primer_max_mismatch: {post_primer_trim_options['max_mismatch']}")
                run_logger.info(f"# post_prep.primer_max_error_rate: {post_primer_trim_options['max_error_rate']}")
                run_logger.info(f"# post_prep.primer_min_overlap_bp: {post_primer_trim_options['min_overlap_bp']}")
                run_logger.info(
                    f"# post_prep.primer_min_overlap_ratio: {post_primer_trim_options['min_overlap_ratio']}"
                )
                run_logger.info(f"# post_prep.primer_end_max_offset: {post_primer_trim_options['end_max_offset']}")
                run_logger.info(f"# post_prep.primer_iter_enable: {post_primer_trim_options['iter_enable']}")
                run_logger.info(f"# post_prep.primer_iter_max_rounds: {post_primer_trim_options['iter_max_rounds']}")
                run_logger.info(f"# post_prep.primer_sidecar_format: {post_primer_trim_options['sidecar_format']}")
                run_logger.info(f"# post_prep.primer_recheck_tool: {post_primer_trim_options['recheck_tool']}")
                run_logger.info(
                    f"# post_prep.primer_recheck_min_identity: {post_primer_trim_options['recheck_min_identity']}"
                )
                run_logger.info(
                    f"# post_prep.primer_recheck_min_query_cov: {post_primer_trim_options['recheck_min_query_cov']}"
                )
            if warnings:
                run_logger.info("# warnings:")
                for warning in warnings:
                    run_logger.info(f"# - {warning}")

            lock = Lock()
            run_logger.info(f"# workers: {workers}")

            with tempfile.TemporaryDirectory(prefix="taxondbbuilder_spool_") as spool_dir_name:
                spool_dir = Path(spool_dir_name)
                ncbi_spool_path = spool_dir / "ncbi_records.jsonl"
                bold_spool_path = spool_dir / "bold_records.jsonl"

                with ncbi_spool_path.open("w", encoding="utf-8") as ncbi_spool_f, bold_spool_path.open(
                    "w", encoding="utf-8"
                ) as bold_spool_f, progress:
                    if uses_ncbi:
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
                                        records = extract_ncbi_records_from_genbank_chunk(
                                            chunk,
                                            marker_rules,
                                            acc_to_seqs,
                                            counters,
                                            dup_accessions,
                                            lock,
                                            progress,
                                            task_id,
                                            taxid,
                                            dump_gb,
                                        )
                                        append_records_to_spool(records, ncbi_spool_f, lock)
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

                    if uses_bold:
                        ncbi_accessions = set(acc_to_seqs.keys())
                        for resolved_taxon in resolved_taxa:
                            try:
                                bold_rows, bold_stats = fetch_bold_records_for_taxon(
                                    resolved_taxon.scientific_name,
                                    marker_keys,
                                    marker_map,
                                    cfg.get("bold"),
                                )
                            except BoldApiError as exc:
                                raise typer.BadParameter(str(exc)) from exc
                            counters["total_records"] += int(bold_stats.get("downloaded_rows") or 0)
                            counters["matched_records"] += int(bold_stats.get("matched_rows") or 0)
                            counters["matched_features"] += int(bold_stats.get("matched_rows") or 0)
                            run_logger.info(
                                "# bold query:"
                                f" taxon={resolved_taxon.scientific_name}"
                                f" normalized={bold_stats.get('normalized_query')}"
                                f" specimens={bold_stats.get('specimen_count')}"
                                f" downloaded={bold_stats.get('downloaded_rows')}"
                                f" matched={bold_stats.get('matched_rows')}"
                            )
                            if bold_stats.get("specimen_count") == 0:
                                console.print(
                                    f"[yellow]BOLD {resolved_taxon.scientific_name}: 0 records "
                                    f"(NCBI scientific name mismatch may be involved)[/yellow]"
                                )
                                continue

                            bold_records = [
                                build_bold_canonical_record(row, marker_map, output_cfg) for row in bold_rows
                            ]
                            kept_bold_records: List[CanonicalRecord] = []
                            for record in bold_records:
                                accession_tokens = parse_accession_tokens(record.accession)
                                if source == BuildSource.BOTH and accession_tokens and any(
                                    token in ncbi_accessions for token in accession_tokens
                                ):
                                    record.linked_to_ncbi = True
                                    record.emitted_to_fasta = False
                                    record.skip_reason = "linked_by_insdcacs"
                                    source_merge_rows.append(build_source_merge_row(record))
                                    continue
                                kept_bold_records.append(record)
                            append_records_to_spool(kept_bold_records, bold_spool_f, lock)

                ncbi_records = load_records_from_spool(ncbi_spool_path)
                bold_records = load_records_from_spool(bold_spool_path)
                ncbi_records.sort(key=canonical_record_sort_key)
                bold_records.sort(key=canonical_record_sort_key)

                with out_path.open("w", encoding="utf-8") as out_f:
                    if ncbi_records:
                        emit_records_to_fasta(
                            ncbi_records,
                            out_f,
                            counters,
                            emitted_records,
                            lock,
                            source_merge_rows=source_merge_rows,
                        )
                    if bold_records:
                        emit_records_to_fasta(
                            bold_records,
                            out_f,
                            counters,
                            emitted_records,
                            lock,
                            source_merge_rows=source_merge_rows,
                        )

            duplicate_records_report_path: Optional[Path] = None
            duplicate_groups_report_path: Optional[Path] = None
            source_merge_path = write_source_merge_csv(out_path, source_merge_rows)
            run_logger.info(f"# source_merge_csv: {source_merge_path}")
            if post_prep:
                before_post_prep = counters["kept_records"]
                run_logger.info(f"# kept records before post_prep: {before_post_prep}")

                if PostPrepStep.PRIMER_TRIM.value in post_prep_steps_run:
                    primer_stats = apply_post_prep_primer_trim(
                        out_path,
                        post_primer_forward,
                        post_primer_reverse,
                        options=post_primer_trim_options,
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
                        f" confidence_high={primer_stats['confidence_high']}"
                        f" confidence_medium={primer_stats['confidence_medium']}"
                        f" confidence_low={primer_stats['confidence_low']}"
                        f" rounds_run={primer_stats['rounds_run']}"
                        f" best_round={primer_stats['best_round']}"
                        f" high_conf_rate={primer_stats['high_conf_rate']:.4f}"
                    )
                    if primer_stats.get("sidecar_path"):
                        run_logger.info(f"# post_prep primer sidecar: {primer_stats['sidecar_path']}")
                    if primer_stats.get("retained_path"):
                        run_logger.info(f"# post_prep primer retained_fasta: {primer_stats['retained_path']}")
                    run_logger.info(
                        "# post_prep primer recheck:"
                        f" tool={primer_stats.get('recheck_tool', 'off')}"
                        f" attempted={primer_stats.get('recheck_attempted', 0)}"
                        f" rescued={primer_stats.get('recheck_rescued', 0)}"
                        f" error={primer_stats.get('recheck_error') or 'none'}"
                    )

                if PostPrepStep.LENGTH_FILTER.value in post_prep_steps_run:
                    length_stats = apply_post_prep_length_filter(out_path, post_min, post_max)
                    counters["kept_records"] = length_stats["after"]
                    run_logger.info(
                        "# post_prep length filter:"
                        f" before={length_stats['before']} after={length_stats['after']}"
                        f" removed={length_stats['removed']}"
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
            console.print(f"Source merge CSV: {source_merge_path}")
            if acc_species_stats["unmapped_records"] > 0:
                console.print(
                    "[yellow]WARNING:[/yellow] Some final FASTA records could not be mapped to source ACC/organism. "
                    "See .log for details."
                )
    finally:
        close_run_logger(run_logger)


if __name__ == "__main__":
    app()

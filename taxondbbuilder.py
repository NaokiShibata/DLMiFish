import io
import os
import re
import time
from datetime import datetime
from pathlib import Path
from queue import Queue
from threading import Event, Lock, Thread
from typing import Dict, Iterable, List, Optional, Tuple

import typer
from Bio import Entrez, SeqIO
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

        if record_matched:
            with lock:
                counters["matched_records"] += 1


def load_config(path: Path) -> Dict:
    if not path.exists():
        raise typer.BadParameter(f"Config file not found: {path}")
    with path.open("rb") as f:
        data = tomllib.load(f)

    markers_file = data.get("markers_file")
    if markers_file:
        if not isinstance(markers_file, str):
            raise typer.BadParameter("markers_file must be a string path.")
        markers_path = Path(markers_file)
        if not markers_path.is_absolute():
            markers_path = path.parent / markers_path
        if not markers_path.exists():
            raise typer.BadParameter(f"Markers file not found: {markers_path}")
        with markers_path.open("rb") as f:
            markers_data = tomllib.load(f)
        markers_from_file = markers_data.get("markers")
        if not isinstance(markers_from_file, dict) or not markers_from_file:
            raise typer.BadParameter("Markers file must define a non-empty [markers] section.")
        data.setdefault("markers", {})
        if not isinstance(data["markers"], dict):
            raise typer.BadParameter("[markers] must be a table (dict).")
        merged = dict(markers_from_file)
        merged.update(data["markers"])
        data["markers"] = merged

    if "ncbi" not in data:
        raise typer.BadParameter("Missing [ncbi] section in config.")
    if "markers" not in data or not data["markers"]:
        raise typer.BadParameter("Missing [markers] section in config.")

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
    terms: List[str] = []
    if not filters:
        return terms

    organelle = filters.get("organelle")
    if organelle:
        terms.append(f"{organelle}[filter]")

    source = filters.get("source")
    if source:
        terms.append(f"{source}[filter]")

    length_min = filters.get("length_min")
    length_max = filters.get("length_max")
    if length_min is not None or length_max is not None:
        lmin = int(length_min) if length_min is not None else 0
        lmax = int(length_max) if length_max is not None else 1000000000
        terms.append(f"{lmin}[SLEN] : {lmax}[SLEN]")

    date_from = filters.get("date_from")
    date_to = filters.get("date_to")
    if date_from or date_to:
        dfrom = date_from or "1800/01/01"
        dto = date_to or "3000/12/31"
        terms.append(f"{dfrom}[PDAT] : {dto}[PDAT]")

    include_keywords = filters.get("include_keywords")
    if include_keywords:
        inc = [f'"{k}"[All Fields]' for k in include_keywords]
        terms.append("(" + " OR ".join(inc) + ")")

    exclude_keywords = filters.get("exclude_keywords")
    if exclude_keywords:
        exc = [f'"{k}"[All Fields]' for k in exclude_keywords]
        terms.append("NOT (" + " OR ".join(exc) + ")")

    extra = filters.get("extra")
    if extra:
        terms.append(str(extra))

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
) -> Iterable[str]:
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

    yield f"__COUNT__={count}"
    if count == 0:
        return

    if use_history and webenv and query_key:
        for start in range(0, count, per_query):
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
            yield data
            time.sleep(delay_sec)
        return

    for start in range(0, count, per_query):
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
            continue
        fetch_handle = Entrez.efetch(
            db=db,
            rettype=rettype,
            retmode=retmode,
            id=",".join(ids),
        )
        data = fetch_handle.read()
        yield data
        time.sleep(delay_sec)


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


def write_log(log_path: Path, lines: List[str]) -> None:
    with log_path.open("w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


@app.command()
def build(
    config: Path = typer.Option(..., "--config", "-c", help="Path to TOML config file."),
    taxon: List[str] = typer.Option(..., "--taxon", "-t", help="Taxon (taxid or scientific name)."),
    marker: List[str] = typer.Option(..., "--marker", "-m", help="Marker key or prefix."),
    out: Optional[Path] = typer.Option(None, "--out", "-o", help="Output file or directory."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print query and exit."),
    workers: int = typer.Option(2, "--workers", "-w", help="Number of extraction workers."),
    output_prefix: str = typer.Option(
        "taxondbbuilder_",
        "--output-prefix",
        help="Prefix added to output FASTA filename.",
    ),
):
    """
    Build a FASTA database by downloading GenBank records and extracting features.

    Examples (Teleostomi):
      taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s
      taxondbbuilder.py build -c configs/db.toml -t "Salmo salar" -m mitogenome
      taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s --workers 2
      taxondbbuilder.py build -c configs/db.toml -t 117570 -m 12s --output-prefix "mifish"
    """
    cfg = load_config(config)
    ncbi_cfg = cfg.get("ncbi", {})
    filters_cfg = cfg.get("filters", {})
    output_cfg = cfg.get("output", {})
    taxon_noexp = bool(cfg.get("taxon", {}).get("noexp", False))

    setup_entrez(ncbi_cfg)
    marker_map = normalize_marker_map(cfg.get("markers", {}))

    marker_keys = [resolve_marker_key(m, marker_map) for m in marker]
    marker_query = build_marker_query(marker_keys, marker_map)
    output_prefix = output_prefix.strip()
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

    taxids: List[str] = []
    warnings: List[str] = []
    for t in taxon:
        taxid, warn = resolve_taxid(t)
        taxids.append(taxid)
        if warn:
            warnings.append(warn)

    out_path = build_output_path(out, taxids, marker_keys, output_prefix=output_prefix)
    log_path = out_path.with_suffix(out_path.suffix + ".log")

    log_lines = []
    log_lines.append(f"# started: {datetime.now().isoformat()}")
    log_lines.append(f"# config: {config}")
    log_lines.append(f"# taxon input: {taxon}")
    log_lines.append(f"# taxids: {taxids}")
    log_lines.append(f"# markers: {marker_keys}")
    if warnings:
        log_lines.append("# warnings:")
        log_lines.extend([f"# - {w}" for w in warnings])

    print_header()
    render_run_table(config, taxids, marker_keys, out_path, filters_cfg)
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

    progress = Progress(
        SpinnerColumn(),
        TextColumn("[bold]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
        disable=not console.is_terminal,
    )

    lock = Lock()
    log_lines.append(f"# workers: {workers}")

    with out_path.open("w", encoding="utf-8") as out_f, progress:
        for taxid in taxids:
            query = build_query(taxid, marker_query, filters_cfg, taxon_noexp)
            log_lines.append(f"# query taxid={taxid}: {query}")
            delay_sec = default_delay(ncbi_cfg)
            data_iter = fetch_genbank(query, ncbi_cfg, delay_sec)
            count_line = next(data_iter)
            count = int(count_line.split("=", 1)[1])
            log_lines.append(f"# query count taxid={taxid}: {count}")
            if count == 0:
                console.print(f"[yellow]taxid {taxid}: 0 records[/yellow]")
                continue

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
                        process_genbank_chunk(
                            item,
                            marker_rules,
                            acc_to_seqs,
                            out_f,
                            counters,
                            dup_accessions,
                            lock,
                            progress,
                            task_id,
                        )
                    except Exception as exc:
                        errors.append(exc)
                        stop_event.set()
                    finally:
                        q.task_done()

            threads = [Thread(target=worker, daemon=True) for _ in range(workers)]
            for t in threads:
                t.start()

            for chunk in data_iter:
                if stop_event.is_set():
                    break
                if not chunk:
                    continue
                q.put(chunk)

            for _ in threads:
                q.put(None)
            q.join()
            for t in threads:
                t.join()
            if errors:
                raise errors[0]

    log_lines.append(f"# total records: {counters['total_records']}")
    log_lines.append(f"# matched records: {counters['matched_records']}")
    log_lines.append(f"# matched features: {counters['matched_features']}")
    log_lines.append(f"# kept records: {counters['kept_records']}")
    log_lines.append(f"# skipped duplicates (same accession+sequence): {counters['skipped_same']}")
    log_lines.append(f"# kept duplicates (same accession, different sequence): {counters['duplicated_diff']}")
    if dup_accessions:
        log_lines.append("# duplicate accessions with different sequences:")
        for acc, count in sorted(dup_accessions.items()):
            log_lines.append(f"# - {acc}: {count} sequences")
    log_lines.append(f"# output: {out_path}")
    log_lines.append(f"# finished: {datetime.now().isoformat()}")

    write_log(log_path, log_lines)
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


if __name__ == "__main__":
    app()

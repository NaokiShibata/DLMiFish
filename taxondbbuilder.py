# -*- coding: utf-8 -*-

import io
import os
import re
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import typer
from Bio import Entrez
from Bio import SeqIO

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    import tomli as tomllib

app = typer.Typer(
    add_completion=False,
    help="TaxonDBBuilder - build a generic NCBI FASTA database by taxon and marker.",
)


def load_config(path: Path) -> Dict:
    if not path.exists():
        raise typer.BadParameter(f"Config file not found: {path}")
    with path.open("rb") as f:
        data = tomllib.load(f)

    if "ncbi" not in data:
        raise typer.BadParameter("Missing [ncbi] section in config.")
    if "markers" not in data or not data["markers"]:
        raise typer.BadParameter("Missing [markers] section in config.")

    return data


def setup_entrez(ncbi_cfg: Dict) -> None:
    email = ncbi_cfg.get("email") or os.environ.get("NCBI_EMAIL")
    api_key = ncbi_cfg.get("api_key") or os.environ.get("NCBI_API_KEY")

    if not email:
        typer.echo("WARNING: NCBI email is not set. Set ncbi.email or NCBI_EMAIL.")
    Entrez.email = email or ""

    if api_key:
        Entrez.api_key = api_key


def normalize_marker_map(markers_cfg: Dict) -> Dict[str, Dict]:
    marker_map: Dict[str, Dict] = {}
    for key, cfg in markers_cfg.items():
        phrases = cfg.get("phrases")
        if not phrases or not isinstance(phrases, list):
            raise typer.BadParameter(f"markers.{key}.phrases must be a list of strings.")
        aliases = cfg.get("aliases") or []
        marker_map[key] = {"phrases": phrases, "aliases": aliases}
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


def build_marker_query(marker_keys: List[str], marker_map: Dict[str, Dict]) -> str:
    phrase_terms: List[str] = []
    for key in marker_keys:
        for phrase in marker_map[key]["phrases"]:
            phrase_escaped = phrase.replace('"', '\\"')
            phrase_terms.append(f'"{phrase_escaped}"[All Fields]')
    if not phrase_terms:
        raise typer.BadParameter("No marker phrases resolved.")
    if len(phrase_terms) == 1:
        return phrase_terms[0]
    return "(" + " OR ".join(phrase_terms) + ")"


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


def fetch_fasta(
    query: str,
    ncbi_cfg: Dict,
    delay_sec: float,
) -> Iterable[str]:
    db = ncbi_cfg.get("db", "nucleotide")
    rettype = ncbi_cfg.get("rettype", "fasta")
    retmode = ncbi_cfg.get("retmode", "text")
    per_query = int(ncbi_cfg.get("per_query", 100))
    use_history = bool(ncbi_cfg.get("use_history", True))

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


def build_output_path(out: Optional[Path], taxids: List[str], markers: List[str]) -> Path:
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
    return out_dir / f"{taxon_label}__{marker_label}.fasta"


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
):
    cfg = load_config(config)
    ncbi_cfg = cfg.get("ncbi", {})
    filters_cfg = cfg.get("filters", {})
    taxon_noexp = bool(cfg.get("taxon", {}).get("noexp", False))

    setup_entrez(ncbi_cfg)
    marker_map = normalize_marker_map(cfg.get("markers", {}))

    marker_keys = [resolve_marker_key(m, marker_map) for m in marker]
    marker_query = build_marker_query(marker_keys, marker_map)

    taxids: List[str] = []
    warnings: List[str] = []
    for t in taxon:
        taxid, warn = resolve_taxid(t)
        taxids.append(taxid)
        if warn:
            warnings.append(warn)

    out_path = build_output_path(out, taxids, marker_keys)
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

    if dry_run:
        for taxid in taxids:
            query = build_query(taxid, marker_query, filters_cfg, taxon_noexp)
            typer.echo(query)
        return

    acc_to_seqs: Dict[str, set] = {}
    total_records = 0
    kept_records = 0
    skipped_same = 0
    duplicated_diff = 0
    dup_accessions: Dict[str, int] = {}

    with out_path.open("w", encoding="utf-8") as out_f:
        for taxid in taxids:
            query = build_query(taxid, marker_query, filters_cfg, taxon_noexp)
            log_lines.append(f"# query taxid={taxid}: {query}")
            delay_sec = default_delay(ncbi_cfg)
            data_iter = fetch_fasta(query, ncbi_cfg, delay_sec)
            count_line = next(data_iter)
            count = int(count_line.split("=", 1)[1])
            log_lines.append(f"# query count taxid={taxid}: {count}")
            typer.echo(f"taxid {taxid}: {count} records")

            for chunk in data_iter:
                if not chunk:
                    continue
                handle = io.StringIO(chunk)
                for record in SeqIO.parse(handle, "fasta"):
                    total_records += 1
                    acc = record.id
                    seq = str(record.seq).upper()
                    seqs = acc_to_seqs.setdefault(acc, set())
                    if seq in seqs:
                        skipped_same += 1
                        continue

                    dup_index = None
                    if seqs:
                        dup_index = len(seqs) + 1
                        duplicated_diff += 1
                        dup_accessions[acc] = dup_accessions.get(acc, 1) + 1

                    seqs.add(seq)

                    desc_parts = record.description.split(None, 1)
                    rest = desc_parts[1] if len(desc_parts) > 1 else ""
                    if dup_index:
                        new_id = f"{acc}_dup{dup_index}"
                        header = f"{new_id} {rest} dup{dup_index}".strip()
                    else:
                        header = f"{acc} {rest}".strip()

                    out_f.write(f">{header}\n")
                    out_f.write(f"{seq}\n")
                    kept_records += 1

    log_lines.append(f"# total records: {total_records}")
    log_lines.append(f"# kept records: {kept_records}")
    log_lines.append(f"# skipped duplicates (same accession+sequence): {skipped_same}")
    log_lines.append(f"# kept duplicates (same accession, different sequence): {duplicated_diff}")
    if dup_accessions:
        log_lines.append("# duplicate accessions with different sequences:")
        for acc, count in sorted(dup_accessions.items()):
            log_lines.append(f"# - {acc}: {count} sequences")
    log_lines.append(f"# output: {out_path}")
    log_lines.append(f"# finished: {datetime.now().isoformat()}")

    write_log(log_path, log_lines)
    typer.echo(f"Output FASTA: {out_path}")
    typer.echo(f"Log: {log_path}")
    if dup_accessions:
        typer.echo("WARNING: duplicate accessions with different sequences were kept. See log for details.")


if __name__ == "__main__":
    app()

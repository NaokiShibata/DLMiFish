#!/usr/bin/env python3
"""
Build taxonomy lookup artifacts for Tauri GUI from NCBI taxdump names.dmp.

Outputs:
  - taxid_scientific_name.csv
  - taxonomy.db (SQLite, table: taxonomy)
"""

from __future__ import annotations

import argparse
import csv
import sqlite3
from pathlib import Path
from typing import Iterable, Iterator, Tuple


def iter_scientific_names(names_dmp: Path) -> Iterator[Tuple[int, str]]:
    """
    Yield (tax_id, scientific_name) from names.dmp where name_class == scientific name.
    """
    with names_dmp.open("r", encoding="utf-8") as fin:
        for raw in fin:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            # taxdump lines end with '\t|'
            if line.endswith("\t|"):
                line = line[:-2]
            parts = line.split("\t|\t")
            if len(parts) < 4:
                continue

            tax_id_s = parts[0].strip()
            name_txt = parts[1].strip()
            name_class = parts[3].strip()
            if name_class != "scientific name":
                continue
            if not tax_id_s.isdigit() or not name_txt:
                continue
            yield int(tax_id_s), name_txt


def extract_scientific_names_csv(names_dmp: Path, output_csv: Path) -> int:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with output_csv.open("w", newline="", encoding="utf-8") as fout:
        writer = csv.writer(fout)
        writer.writerow(["tax_id", "scientific_name"])
        for tax_id, sci_name in iter_scientific_names(names_dmp):
            writer.writerow([tax_id, sci_name])
            count += 1
    return count


def iter_csv_rows(input_csv: Path) -> Iterator[Tuple[int, str]]:
    with input_csv.open("r", newline="", encoding="utf-8") as fin:
        reader = csv.DictReader(fin)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {input_csv}")

        # Support both canonical and reversed column order.
        headers = {h.strip(): h for h in reader.fieldnames if h}
        tax_col = headers.get("tax_id")
        sci_col = headers.get("scientific_name")
        if not tax_col or not sci_col:
            raise ValueError("CSV must contain headers: tax_id, scientific_name")

        for row in reader:
            tax_s = str(row.get(tax_col, "")).strip()
            sci_name = str(row.get(sci_col, "")).strip()
            if not tax_s.isdigit() or not sci_name:
                continue
            yield int(tax_s), sci_name


def _chunked(rows: Iterable[Tuple[int, str]], size: int = 10000) -> Iterator[list[Tuple[int, str]]]:
    buf: list[Tuple[int, str]] = []
    for row in rows:
        buf.append(row)
        if len(buf) >= size:
            yield buf
            buf = []
    if buf:
        yield buf


def build_sqlite_from_csv(input_csv: Path, output_db: Path) -> int:
    output_db.parent.mkdir(parents=True, exist_ok=True)
    if output_db.exists():
        output_db.unlink()

    conn = sqlite3.connect(str(output_db))
    inserted = 0
    try:
        cur = conn.cursor()
        cur.execute(
            """
            CREATE TABLE taxonomy (
                tax_id INTEGER PRIMARY KEY,
                scientific_name TEXT NOT NULL
            )
            """
        )

        for batch in _chunked(iter_csv_rows(input_csv), size=20000):
            cur.executemany(
                "INSERT INTO taxonomy (tax_id, scientific_name) VALUES (?, ?)",
                batch,
            )
            inserted += len(batch)

        cur.execute(
            "CREATE INDEX idx_taxonomy_scientific_name_nocase "
            "ON taxonomy(scientific_name COLLATE NOCASE)"
        )
        conn.commit()
    finally:
        conn.close()
    return inserted


def build_from_names_dmp(names_dmp: Path, output_csv: Path, output_db: Path) -> tuple[int, int]:
    extracted = extract_scientific_names_csv(names_dmp, output_csv)
    inserted = build_sqlite_from_csv(output_csv, output_db)
    return extracted, inserted


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build tax_id <-> scientific_name artifacts from NCBI taxdump names.dmp"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_extract = sub.add_parser("extract-csv", help="Extract scientific names from names.dmp to CSV")
    p_extract.add_argument("--names-dmp", type=Path, required=True, help="Path to names.dmp")
    p_extract.add_argument("--out-csv", type=Path, required=True, help="Output CSV path")

    p_sqlite = sub.add_parser("build-sqlite", help="Build SQLite taxonomy.db from CSV")
    p_sqlite.add_argument("--in-csv", type=Path, required=True, help="Input CSV path")
    p_sqlite.add_argument("--out-db", type=Path, required=True, help="Output SQLite DB path")

    p_all = sub.add_parser("all", help="Run extract-csv and build-sqlite in one step")
    p_all.add_argument("--names-dmp", type=Path, required=True, help="Path to names.dmp")
    p_all.add_argument("--out-csv", type=Path, required=True, help="Output CSV path")
    p_all.add_argument("--out-db", type=Path, required=True, help="Output SQLite DB path")

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.command == "extract-csv":
        count = extract_scientific_names_csv(args.names_dmp, args.out_csv)
        print(f"extracted {count} rows -> {args.out_csv}")
        return 0

    if args.command == "build-sqlite":
        count = build_sqlite_from_csv(args.in_csv, args.out_db)
        print(f"inserted {count} rows -> {args.out_db}")
        return 0

    if args.command == "all":
        extracted, inserted = build_from_names_dmp(args.names_dmp, args.out_csv, args.out_db)
        print(f"extracted {extracted} rows -> {args.out_csv}")
        print(f"inserted {inserted} rows -> {args.out_db}")
        return 0

    raise RuntimeError(f"unknown command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())

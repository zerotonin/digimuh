#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  Gouna respiration data loader                               ║
# ║  « append new data, skip existing rows »                     ║
# ╚══════════════════════════════════════════════════════════════╝
"""Ingest additional gouna respiration CSVs into the existing cow.db,
skipping any rows that already exist (by animal_id + timestamp).

Expected CSV format (one file per animal)::

    timestamp,respirationFrequency
    2022-11-03 00:00:22.701000+00:00,31.0
    ...

Animal ID is extracted from the filename (first numeric block matching
the 15-digit EU ear tag pattern, e.g.
``276001263589826_gouna_2021-04-01_2024-09-30.csv``).

Usage::

    python ingest_gouna.py --db cow.db --input /path/to/new_gouna_csvs/

    # Or a single file:
    python ingest_gouna.py --db cow.db --input 276001263589826_gouna_*.csv

The script:
1. Reads each CSV
2. Extracts animal_id from the filename
3. Registers the animal in the ``animals`` dimension table if new
4. Registers a source file entry in ``source_files`` if new
5. Inserts rows into ``gouna`` using INSERT OR IGNORE to skip duplicates
   (requires a UNIQUE index on animal_id + timestamp — the script
   creates one if it doesn't exist)
6. Reports how many rows were inserted vs skipped per file
"""

from __future__ import annotations

import argparse
import logging
import re
import sqlite3
from pathlib import Path

import pandas as pd

log = logging.getLogger("digimuh.ingest_gouna")

# Pattern to extract 15-digit EU ear tag from filename
EAR_TAG_RE = re.compile(r"(276\d{12})")


def ensure_unique_index(con: sqlite3.Connection) -> None:
    """Create a UNIQUE index on gouna(animal_id, timestamp) if missing.

    This enables INSERT OR IGNORE to skip duplicates efficiently.
    If a non-unique index already exists, it is dropped and replaced.
    """
    cur = con.cursor()

    # Check existing indexes on gouna
    cur.execute("PRAGMA index_list('gouna')")
    indexes = cur.fetchall()

    for idx in indexes:
        idx_name = idx[1]
        is_unique = idx[2]
        cur.execute(f"PRAGMA index_info('{idx_name}')")
        cols = [row[2] for row in cur.fetchall()]
        if cols == ["animal_id", "timestamp"]:
            if is_unique:
                log.info("UNIQUE index '%s' already exists on gouna(animal_id, timestamp)", idx_name)
                return
            else:
                log.info("Dropping non-unique index '%s' to replace with UNIQUE", idx_name)
                cur.execute(f"DROP INDEX IF EXISTS {idx_name}")
                break

    log.info("Creating UNIQUE index on gouna(animal_id, timestamp) …")
    con.execute("""
        CREATE UNIQUE INDEX IF NOT EXISTS idx_gouna_animal_ts_unique
        ON gouna(animal_id, "timestamp")
    """)
    con.commit()
    log.info("  Done.")


def extract_animal_id(filepath: Path) -> int | None:
    """Extract the 15-digit EU ear tag from a filename.

    Args:
        filepath: Path to the CSV file.

    Returns:
        Animal ID as integer, or None if not found.
    """
    match = EAR_TAG_RE.search(filepath.stem)
    if match:
        return int(match.group(1))
    return None


def ensure_animal(con: sqlite3.Connection, animal_id: int) -> None:
    """Register an animal in the dimension table if not present."""
    existing = con.execute(
        "SELECT animal_id FROM animals WHERE animal_id = ?", (animal_id,)
    ).fetchone()
    if existing is None:
        con.execute("INSERT INTO animals (animal_id) VALUES (?)", (animal_id,))
        log.info("  Registered new animal: %d", animal_id)


def ensure_source_file(con: sqlite3.Connection, filepath: Path) -> int:
    """Register a source file and return its ID."""
    fname = filepath.name
    existing = con.execute(
        "SELECT file_id FROM source_files WHERE file_name = ?", (fname,)
    ).fetchone()
    if existing:
        return existing[0]
    cur = con.execute(
        "INSERT INTO source_files (file_name, file_path, system) VALUES (?, ?, ?)",
        (fname, str(filepath), "gouna"),
    )
    return cur.lastrowid


def ingest_file(con: sqlite3.Connection, filepath: Path, dry_run: bool = False) -> dict:
    """Ingest a single gouna CSV, skipping existing rows.

    Args:
        con: Database connection.
        filepath: Path to the CSV file.
        dry_run: If True, report what would happen without inserting.

    Returns:
        Dict with animal_id, total_rows, inserted, skipped.
    """
    animal_id = extract_animal_id(filepath)
    if animal_id is None:
        log.warning("  Could not extract animal ID from: %s", filepath.name)
        return {"file": filepath.name, "error": "no animal ID in filename"}

    df = pd.read_csv(filepath)
    if df.empty:
        return {"file": filepath.name, "animal_id": animal_id, "total_rows": 0,
                "inserted": 0, "skipped": 0}

    # Normalise column names
    col_map = {}
    for col in df.columns:
        if "timestamp" in col.lower():
            col_map[col] = "timestamp"
        elif "respiration" in col.lower() or "freq" in col.lower():
            col_map[col] = "respirationfrequency"
    df = df.rename(columns=col_map)

    if "timestamp" not in df.columns or "respirationfrequency" not in df.columns:
        log.warning("  Missing expected columns in %s. Found: %s",
                    filepath.name, list(df.columns))
        return {"file": filepath.name, "error": f"missing columns: {list(df.columns)}"}

    total = len(df)

    if dry_run:
        # Count how many already exist
        existing = con.execute(
            "SELECT COUNT(*) FROM gouna WHERE animal_id = ?", (animal_id,)
        ).fetchone()[0]
        return {"file": filepath.name, "animal_id": animal_id,
                "total_rows": total, "existing_in_db": existing,
                "dry_run": True}

    # Ensure animal and source file exist
    ensure_animal(con, animal_id)
    file_id = ensure_source_file(con, filepath)

    # Prepare data
    df["animal_id"] = animal_id
    df["source_file_id"] = file_id

    inserted_before = con.execute("SELECT COUNT(*) FROM gouna WHERE animal_id = ?",
                                   (animal_id,)).fetchone()[0]

    # Use INSERT OR IGNORE — duplicates on UNIQUE(animal_id, timestamp)
    # are silently skipped.  executemany is orders of magnitude faster
    # than row-by-row inserts.
    rows = [
        (animal_id, row["timestamp"], row["respirationfrequency"], file_id)
        for _, row in df.iterrows()
    ]
    con.executemany(
        'INSERT OR IGNORE INTO gouna (animal_id, "timestamp", respirationfrequency, source_file_id) '
        'VALUES (?, ?, ?, ?)',
        rows,
    )
    con.commit()

    inserted_after = con.execute("SELECT COUNT(*) FROM gouna WHERE animal_id = ?",
                                  (animal_id,)).fetchone()[0]
    inserted = inserted_after - inserted_before
    skipped = total - inserted

    return {
        "file": filepath.name,
        "animal_id": animal_id,
        "total_rows": total,
        "inserted": inserted,
        "skipped": skipped,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ingest additional gouna respiration CSVs into cow.db "
                    "(duplicates are skipped)",
    )
    parser.add_argument("--db", type=Path, required=True,
                        help="Path to existing cow.db")
    parser.add_argument("--input", type=Path, required=True, nargs="+",
                        help="CSV files or directory containing gouna CSVs")
    parser.add_argument("--dry-run", action="store_true",
                        help="Report what would be inserted without writing")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(name)-20s │ %(message)s",
        datefmt="%H:%M:%S",
    )

    # Collect CSV files
    csv_files = []
    for p in args.input:
        if p.is_dir():
            csv_files.extend(sorted(p.glob("*.csv")))
        elif p.is_file() and p.suffix == ".csv":
            csv_files.append(p)
        else:
            log.warning("Skipping: %s", p)

    if not csv_files:
        log.error("No CSV files found.")
        return

    log.info("Found %d gouna CSV file(s)", len(csv_files))

    con = sqlite3.connect(str(args.db))
    con.execute("PRAGMA journal_mode = WAL")
    con.execute("PRAGMA cache_size = -256000")  # 256 MB

    if not args.dry_run:
        ensure_unique_index(con)

    total_inserted = 0
    total_skipped = 0
    results = []

    for filepath in csv_files:
        log.info("Processing: %s", filepath.name)
        result = ingest_file(con, filepath, dry_run=args.dry_run)
        results.append(result)

        if "error" in result:
            log.warning("  Error: %s", result["error"])
        elif args.dry_run:
            log.info("  Animal %s: %d rows in CSV, %d already in DB",
                     result.get("animal_id"), result.get("total_rows", 0),
                     result.get("existing_in_db", 0))
        else:
            ins = result.get("inserted", 0)
            skip = result.get("skipped", 0)
            total_inserted += ins
            total_skipped += skip
            log.info("  Animal %s: %d inserted, %d skipped (duplicates)",
                     result.get("animal_id"), ins, skip)

    if not args.dry_run:
        log.info("─" * 50)
        log.info("Total: %d inserted, %d skipped across %d files",
                 total_inserted, total_skipped, len(csv_files))
    else:
        log.info("─" * 50)
        log.info("Dry run complete. No rows were written.")

    con.close()


if __name__ == "__main__":
    main()

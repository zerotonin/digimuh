#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════╗
║  DigiMuh — WAL recovery + index builder                      ║
║  « verbose edition with stall detection »                    ║
╚══════════════════════════════════════════════════════════════╝

Recovers from the OOM-killed ingestion by:
1. Checkpointing the WAL back into the main DB
2. Switching from WAL to DELETE journal mode
3. Building each index individually with per-index timing

Usage:
    python recover_and_index.py /path/to/cow.db
"""

import os
import sqlite3
import sys
import threading
import time
from pathlib import Path


def fmt_size(nbytes: int) -> str:
    """Human-readable file size."""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(nbytes) < 1024:
            return f"{nbytes:.1f} {unit}"
        nbytes /= 1024
    return f"{nbytes:.1f} PB"


def file_monitor(db_path: str, stop_event: threading.Event, label: str) -> None:
    """Background thread that prints DB + WAL size every 30 seconds.

    Args:
        db_path: Path to the .db file.
        stop_event: Set this to stop the monitor.
        label: Current operation label to display.
    """
    wal_path = db_path + "-wal"
    last_db_size = 0
    last_wal_size = 0
    stall_count = 0

    while not stop_event.is_set():
        try:
            db_size = os.path.getsize(db_path) if os.path.exists(db_path) else 0
            wal_size = os.path.getsize(wal_path) if os.path.exists(wal_path) else 0

            db_delta = db_size - last_db_size
            wal_delta = wal_size - last_wal_size

            # Stall detection
            if db_delta == 0 and wal_delta == 0:
                stall_count += 1
            else:
                stall_count = 0

            stall_warn = ""
            if stall_count >= 4:  # 2 minutes no change
                stall_warn = "  ⚠ no size change for 2+ min (may be I/O bound)"

            print(
                f"  [{label}]  "
                f"DB: {fmt_size(db_size)} (Δ {fmt_size(db_delta)})  "
                f"WAL: {fmt_size(wal_size)} (Δ {fmt_size(wal_delta)})"
                f"{stall_warn}",
                flush=True,
            )

            last_db_size = db_size
            last_wal_size = wal_size

        except OSError:
            pass

        stop_event.wait(30)


def start_monitor(db_path: str, label: str) -> tuple[threading.Event, threading.Thread]:
    """Start the background file-size monitor.

    Returns:
        (stop_event, thread) — call stop_event.set() to stop.
    """
    stop = threading.Event()
    t = threading.Thread(target=file_monitor, args=(db_path, stop, label), daemon=True)
    t.start()
    return stop, t


def main() -> None:
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} /path/to/cow.db")
        sys.exit(1)

    db_path = sys.argv[1]
    wal_path = db_path + "-wal"

    if not os.path.exists(db_path):
        print(f"ERROR: {db_path} not found")
        sys.exit(1)

    db_size = os.path.getsize(db_path)
    wal_size = os.path.getsize(wal_path) if os.path.exists(wal_path) else 0

    print("=" * 60)
    print("  DigiMuh WAL RECOVERY + INDEX BUILDER")
    print(f"  DB:   {db_path}  ({fmt_size(db_size)})")
    print(f"  WAL:  {wal_path}  ({fmt_size(wal_size)})")
    print("=" * 60)
    print()

    con = sqlite3.connect(db_path)
    con.execute("PRAGMA cache_size = -128000")  # 128 MB — gentle on RAM
    con.execute("PRAGMA temp_store = FILE")      # temp sort to disk, not RAM

    # ── Step 1: verify data is intact ────────────────────────
    print("Step 1: Checking table row counts …", flush=True)
    cur = con.cursor()
    tables = [
        "animals", "sensors", "barns", "source_files",
        "allocations", "bcs", "diseases", "dwd_weather",
        "gouna", "herdeplus", "hobo_weather", "lorawan",
        "smaxtec_barns", "smaxtec_derived", "smaxtec_events",
        "smaxtec_water_intake",
    ]
    for table in tables:
        try:
            cur.execute(f'SELECT COUNT(*) FROM "{table}"')
            count = cur.fetchone()[0]
            print(f"  {table:30s} {count:>15,}", flush=True)
        except sqlite3.OperationalError as e:
            print(f"  {table:30s} ERROR: {e}", flush=True)
    print()

    # ── Step 2: checkpoint WAL ───────────────────────────────
    if wal_size > 0:
        print("Step 2: Checkpointing WAL → main DB …", flush=True)
        print(f"  WAL is {fmt_size(wal_size)} — this may take 20-60 min on USB", flush=True)

        stop_evt, monitor_thread = start_monitor(db_path, "checkpoint")

        t0 = time.time()
        con.execute("PRAGMA wal_checkpoint(TRUNCATE)")
        con.commit()
        elapsed = time.time() - t0

        stop_evt.set()
        monitor_thread.join(timeout=2)

        print(f"  Checkpoint complete in {elapsed:.0f}s ({elapsed/60:.1f} min)", flush=True)
    else:
        print("Step 2: No WAL file — skipping checkpoint.", flush=True)
    print()

    # ── Step 3: switch to DELETE journal mode ────────────────
    print("Step 3: Switching journal mode WAL → DELETE …", flush=True)
    con.execute("PRAGMA journal_mode = DELETE")
    con.commit()
    print("  Done.", flush=True)
    print()

    # ── Step 4: check existing indexes ───────────────────────
    print("Step 4: Checking existing indexes …", flush=True)
    cur.execute("SELECT name FROM sqlite_master WHERE type='index' AND name LIKE 'idx_%'")
    existing = {row[0] for row in cur.fetchall()}
    if existing:
        print(f"  Found {len(existing)} existing indexes:", flush=True)
        for name in sorted(existing):
            print(f"    {name}", flush=True)
    else:
        print("  No indexes found — will create all.", flush=True)
    print()

    # ── Step 5: build indexes one by one ─────────────────────
    indexes = [
        # Biggest first — get the pain out of the way
        ("idx_smaxtec_derived_animal_ts", "smaxtec_derived",
         '"animal_id", "timestamp"'),
        ("idx_gouna_animal_ts", "gouna",
         '"animal_id", "timestamp"'),
        ("idx_lorawan_sensor_ts", "lorawan",
         '"sensor_id", "timestamp"'),
        ("idx_smaxtec_barns_barn_ts", "smaxtec_barns",
         '"barn_id", "timestamp"'),
        ("idx_herdeplus_animal_ts", "herdeplus",
         '"animal_id", "timestamp"'),
        ("idx_bcs_animal_ts", "bcs",
         '"animal_id", "timestamp"'),
        ("idx_smaxtec_events_animal_ts", "smaxtec_events",
         '"animal_id", "timestamp"'),
        ("idx_smaxtec_water_animal_ts", "smaxtec_water_intake",
         '"animal_id", "timestamp"'),
        ("idx_alloc_animal", "allocations",
         '"animal_id"'),
        ("idx_alloc_enter", "allocations",
         '"datetime_enter"'),
        ("idx_disease_animal", "diseases",
         '"animal_id"'),
        ("idx_disease_start", "diseases",
         '"disease_first_day"'),
        ("idx_dwd_dt", "dwd_weather",
         '"dt"'),
        ("idx_hobo_ts", "hobo_weather",
         '"datetime"'),
        ("idx_source_folder", "source_files",
         '"folder"'),
    ]

    print(f"Step 5: Building {len(indexes)} indexes …", flush=True)
    t_total = time.time()

    for i, (name, table, cols) in enumerate(indexes, 1):
        if name in existing:
            print(f"  [{i:2d}/{len(indexes)}] {name:45s} SKIP (already exists)",
                  flush=True)
            continue

        print(f"  [{i:2d}/{len(indexes)}] {name:45s} building …", flush=True)

        stop_evt, monitor_thread = start_monitor(db_path, name)
        t1 = time.time()

        try:
            con.execute(
                f'CREATE INDEX IF NOT EXISTS "{name}" ON "{table}" ({cols})',
            )
            con.commit()
            elapsed = time.time() - t1
            print(
                f"  [{i:2d}/{len(indexes)}] {name:45s} DONE  {elapsed:.1f}s "
                f"({elapsed/60:.1f} min)",
                flush=True,
            )
        except Exception as e:
            elapsed = time.time() - t1
            print(
                f"  [{i:2d}/{len(indexes)}] {name:45s} FAILED after {elapsed:.1f}s: {e}",
                flush=True,
            )
        finally:
            stop_evt.set()
            monitor_thread.join(timeout=2)

    total_elapsed = time.time() - t_total
    print()
    print("=" * 60)
    print(f"  ALL INDEXES COMPLETE")
    print(f"  Total index time: {total_elapsed:.0f}s ({total_elapsed/60:.1f} min)")
    print(f"  Final DB size:    {fmt_size(os.path.getsize(db_path))}")

    wal_exists = os.path.exists(wal_path)
    if wal_exists:
        print(f"  WAL still exists: {fmt_size(os.path.getsize(wal_path))}")
    else:
        print("  WAL: removed (clean)")
    print("=" * 60)

    con.close()


if __name__ == "__main__":
    main()
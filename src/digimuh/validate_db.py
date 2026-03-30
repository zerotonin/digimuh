#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  DigiMuh post-ingestion sanity check                         ║
# ║  « trust but verify »                                        ║
# ╚══════════════════════════════════════════════════════════════╝
"""Quick validation of the ingested cow database.

Checks row counts, null rates, value ranges, referential integrity,
and temporal coverage.  Run immediately after ingestion to catch
problems before launching analysis.

Usage::

    python -m digimuh.validate_db --db cow.db
"""

from __future__ import annotations

import argparse
import logging
import sqlite3
import sys
from pathlib import Path

log = logging.getLogger("digimuh.validate")

# ─────────────────────────────────────────────────────────────
#  « expected ranges for plausibility checks »
# ─────────────────────────────────────────────────────────────

RANGE_CHECKS = {
    "smaxtec_derived": {
        "temp": (30.0, 45.0, "°C — rumen temperature"),
        "ph": (3.0, 8.5, "pH units"),
        "act": (0.0, 500.0, "activity index"),
    },
    "herdeplus": {
        "herdeplus_milked_mkg": (0.0, 80.0, "kg — single milking yield"),
        "herdeplus_mlp_fat_percent": (1.0, 9.0, "% — milk fat"),
        "herdeplus_mlp_protein_percent": (1.5, 6.0, "% — milk protein"),
        "herdeplus_mlp_f_e": (0.5, 3.0, "— fat:protein ratio"),
        "herdeplus_mlp_cell_count": (0.0, 10000.0, "×1000 cells/mL — SCC"),
    },
    "bcs": {
        "bcs_wert": (1.0, 5.0, "BCS score"),
    },
    "gouna": {
        "respirationfrequency": (5.0, 120.0, "breaths/min"),
    },
    "smaxtec_barns": {
        "temp": (-20.0, 55.0, "°C — barn temperature"),
        "hum": (0.0, 100.0, "% — relative humidity"),
        "temp_hum_index": (20.0, 100.0, "THI index"),
    },
    "smaxtec_water_intake": {
        "water_intake_liter": (0.0, 300.0, "litres — daily water"),
    },
    "dwd_weather": {
        "thi_max": (20.0, 100.0, "THI max"),
        "enthalpy_max": (-10.0, 120.0, "kJ/kg — enthalpy"),
    },
}


def check_table_counts(con: sqlite3.Connection) -> list[str]:
    """Verify all expected tables exist and report row counts.

    Returns:
        List of warning/error messages (empty = all good).
    """
    issues: list[str] = []
    cur = con.cursor()

    expected_tables = [
        "animals", "sensors", "barns", "source_files",
        "allocations", "bcs", "diseases", "dwd_weather",
        "gouna", "herdeplus", "hobo_weather", "lorawan",
        "smaxtec_barns", "smaxtec_derived", "smaxtec_events",
        "smaxtec_water_intake",
    ]

    log.info("─" * 50)
    log.info("TABLE ROW COUNTS")
    log.info("─" * 50)

    for table in expected_tables:
        try:
            cur.execute(f'SELECT COUNT(*) FROM "{table}"')
            count = cur.fetchone()[0]
            status = "OK" if count > 0 else "EMPTY"
            log.info("  %-30s %12s  [%s]", table, f"{count:,}", status)
            if count == 0:
                issues.append(f"WARN: Table '{table}' is empty")
        except sqlite3.OperationalError:
            log.info("  %-30s %12s  [MISSING]", table, "—")
            issues.append(f"ERROR: Table '{table}' does not exist")

    # Dimension table cross-check
    cur.execute("SELECT COUNT(*) FROM animals")
    n_animals = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM source_files")
    n_files = cur.fetchone()[0]
    log.info("")
    log.info("  Animals: %d,  Source files: %d", n_animals, n_files)

    return issues


def check_null_rates(con: sqlite3.Connection) -> list[str]:
    """Report null rates for key columns.

    Returns:
        List of warning messages for suspiciously high null rates.
    """
    issues: list[str] = []
    cur = con.cursor()

    checks = [
        ("smaxtec_derived", "temp"),
        ("smaxtec_derived", "ph"),
        ("smaxtec_derived", "act"),
        ("smaxtec_derived", "rum_index"),
        ("smaxtec_derived", "mot_period"),
        ("herdeplus", "herdeplus_milked_mkg"),
        ("gouna", "respirationfrequency"),
        ("smaxtec_water_intake", "water_intake_liter"),
    ]

    log.info("")
    log.info("─" * 50)
    log.info("NULL RATES (key columns)")
    log.info("─" * 50)

    for table, col in checks:
        try:
            cur.execute(
                f'SELECT COUNT(*) AS total, '
                f'SUM(CASE WHEN "{col}" IS NULL THEN 1 ELSE 0 END) AS nulls '
                f'FROM "{table}"',
            )
            row = cur.fetchone()
            total, nulls = row[0], row[1]
            if total > 0:
                pct = 100.0 * nulls / total
                status = "OK" if pct < 80 else "HIGH"
                log.info(
                    "  %s.%-25s %6.1f%% null  (%s / %s)  [%s]",
                    table, col, pct,
                    f"{nulls:,}", f"{total:,}", status,
                )
                if pct > 95:
                    issues.append(
                        f"WARN: {table}.{col} is {pct:.1f}% null",
                    )
        except sqlite3.OperationalError:
            pass

    return issues


def check_value_ranges(con: sqlite3.Connection) -> list[str]:
    """Check that numeric values fall within plausible ranges.

    Returns:
        List of warnings for out-of-range values.
    """
    issues: list[str] = []
    cur = con.cursor()

    log.info("")
    log.info("─" * 50)
    log.info("VALUE RANGE CHECKS")
    log.info("─" * 50)

    for table, columns in RANGE_CHECKS.items():
        for col, (lo, hi, desc) in columns.items():
            try:
                cur.execute(
                    f'SELECT MIN(CAST("{col}" AS REAL)), '
                    f'MAX(CAST("{col}" AS REAL)), '
                    f'AVG(CAST("{col}" AS REAL)), '
                    f'COUNT(*) '
                    f'FROM "{table}" WHERE "{col}" IS NOT NULL',
                )
                row = cur.fetchone()
                vmin, vmax, vmean, count = row
                if count == 0:
                    continue

                out_lo = vmin is not None and vmin < lo
                out_hi = vmax is not None and vmax > hi
                status = "OK"
                if out_lo or out_hi:
                    status = "RANGE"
                    issues.append(
                        f"WARN: {table}.{col} range [{vmin:.2f}, {vmax:.2f}] "
                        f"exceeds expected [{lo}, {hi}] ({desc})",
                    )

                log.info(
                    "  %s.%-25s min=%8.2f  max=%8.2f  mean=%8.2f  [%s]",
                    table, col, vmin or 0, vmax or 0, vmean or 0, status,
                )
            except sqlite3.OperationalError:
                pass

    return issues


def check_temporal_coverage(con: sqlite3.Connection) -> list[str]:
    """Report the date range of each timestamped table.

    Returns:
        List of warnings for unexpected gaps or ranges.
    """
    issues: list[str] = []
    cur = con.cursor()

    tables_ts = [
        ("smaxtec_derived", "timestamp"),
        ("herdeplus", "timestamp"),
        ("gouna", "timestamp"),
        ("bcs", "timestamp"),
        ("smaxtec_barns", "timestamp"),
        ("smaxtec_events", "timestamp"),
        ("smaxtec_water_intake", "timestamp"),
        ("lorawan", "timestamp"),
        ("hobo_weather", "datetime"),
        ("dwd_weather", "dt"),
    ]

    log.info("")
    log.info("─" * 50)
    log.info("TEMPORAL COVERAGE")
    log.info("─" * 50)

    for table, ts_col in tables_ts:
        try:
            cur.execute(
                f'SELECT MIN("{ts_col}"), MAX("{ts_col}"), COUNT(*) '
                f'FROM "{table}" WHERE "{ts_col}" IS NOT NULL',
            )
            row = cur.fetchone()
            ts_min, ts_max, count = row
            if count > 0:
                log.info(
                    "  %-30s %s  →  %s  (%s rows)",
                    table,
                    str(ts_min)[:10] if ts_min else "?",
                    str(ts_max)[:10] if ts_max else "?",
                    f"{count:,}",
                )
        except sqlite3.OperationalError:
            pass

    return issues


def check_referential_integrity(con: sqlite3.Connection) -> list[str]:
    """Check foreign key relationships between fact and dimension tables.

    Returns:
        List of warnings for orphaned references.
    """
    issues: list[str] = []
    cur = con.cursor()

    fk_checks = [
        ("bcs", "animal_id", "animals", "animal_id"),
        ("gouna", "animal_id", "animals", "animal_id"),
        ("herdeplus", "animal_id", "animals", "animal_id"),
        ("smaxtec_derived", "animal_id", "animals", "animal_id"),
        ("smaxtec_events", "animal_id", "animals", "animal_id"),
        ("smaxtec_water_intake", "animal_id", "animals", "animal_id"),
        ("lorawan", "sensor_id", "sensors", "sensor_id"),
        ("smaxtec_barns", "barn_id", "barns", "barn_id"),
    ]

    log.info("")
    log.info("─" * 50)
    log.info("REFERENTIAL INTEGRITY")
    log.info("─" * 50)

    for fact_table, fk_col, dim_table, pk_col in fk_checks:
        try:
            cur.execute(
                f'SELECT COUNT(*) FROM "{fact_table}" f '
                f'LEFT JOIN "{dim_table}" d ON f."{fk_col}" = d."{pk_col}" '
                f'WHERE d."{pk_col}" IS NULL AND f."{fk_col}" IS NOT NULL',
            )
            orphans = cur.fetchone()[0]
            status = "OK" if orphans == 0 else f"{orphans:,} ORPHANS"
            log.info(
                "  %s.%s → %s.%s  [%s]",
                fact_table, fk_col, dim_table, pk_col, status,
            )
            if orphans > 0:
                issues.append(
                    f"WARN: {fact_table}.{fk_col} has {orphans:,} "
                    f"orphaned references to {dim_table}",
                )
        except sqlite3.OperationalError:
            pass

    return issues


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────

def main() -> None:
    """Run all validation checks and print a summary."""
    parser = argparse.ArgumentParser(
        description="Validate the ingested DigiMuh database.",
    )
    parser.add_argument("--db", type=Path, required=True)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(message)s",
        datefmt="%H:%M:%S",
    )

    if not args.db.is_file():
        log.error("Database not found: %s", args.db)
        sys.exit(1)

    con = sqlite3.connect(str(args.db))
    con.row_factory = sqlite3.Row

    log.info("═" * 50)
    log.info("  DigiMuh DATABASE VALIDATION")
    log.info("  DB: %s", args.db)
    db_size = args.db.stat().st_size / (1024 ** 3)
    log.info("  Size: %.2f GB", db_size)
    log.info("═" * 50)

    all_issues: list[str] = []
    all_issues.extend(check_table_counts(con))
    all_issues.extend(check_null_rates(con))
    all_issues.extend(check_value_ranges(con))
    all_issues.extend(check_temporal_coverage(con))
    all_issues.extend(check_referential_integrity(con))

    log.info("")
    log.info("═" * 50)
    log.info("  SUMMARY")
    log.info("═" * 50)

    errors = [i for i in all_issues if i.startswith("ERROR")]
    warnings = [i for i in all_issues if i.startswith("WARN")]

    if errors:
        log.error("  %d ERRORS:", len(errors))
        for e in errors:
            log.error("    %s", e)

    if warnings:
        log.warning("  %d WARNINGS:", len(warnings))
        for w in warnings:
            log.warning("    %s", w)

    if not all_issues:
        log.info("  All checks passed ✓")

    con.close()
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()

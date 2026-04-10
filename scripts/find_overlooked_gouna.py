#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  Find overlooked gouna animals                                ║
# ║  « animals with respiration data in Neubau but not in TA »   ║
# ╚══════════════════════════════════════════════════════════════╝
"""For each summer (Jun-Sep), find animals that have gouna respiration
data AND were allocated to the Neubau (groups 1005/1006) but are NOT
in the Tierauswahl.

Usage::

    python scripts/find_overlooked_gouna.py --db ~/cow.db --tierauswahl Tierauswahl.xlsx
    python scripts/find_overlooked_gouna.py --db ~/cow.db --tierauswahl Tierauswahl.xlsx --year 2023
"""

from __future__ import annotations

import argparse
import logging
import sqlite3
from pathlib import Path

import pandas as pd

log = logging.getLogger("digimuh.find_gouna")


def find_overlooked(
    con: sqlite3.Connection, tierauswahl: pd.DataFrame, year: int,
) -> pd.DataFrame:
    """Find gouna-equipped Neubau animals not in Tierauswahl for a year.

    Args:
        con: Database connection.
        tierauswahl: Cleaned Tierauswahl DataFrame.
        year: Summer year to check.

    Returns:
        DataFrame with animal_id, n_gouna, first/last timestamp,
        group, Neubau enter/exit.
    """
    summer_start = f"{year}-06-01"
    summer_end = f"{year}-09-30"

    # All animals with gouna data in this summer
    gouna = pd.read_sql_query(
        """SELECT animal_id, COUNT(*) AS n_readings,
                  MIN(timestamp) AS first_ts, MAX(timestamp) AS last_ts
           FROM gouna
           WHERE timestamp >= ? AND timestamp <= ?
             AND respirationfrequency IS NOT NULL
             AND respirationfrequency BETWEEN 5 AND 150
           GROUP BY animal_id
           ORDER BY n_readings DESC""",
        con, params=(summer_start, summer_end),
    )

    if gouna.empty:
        return pd.DataFrame()

    # Animals in Neubau (group 1005/1006) during this summer
    placeholders = ",".join(["?"] * len(gouna))
    alloc = pd.read_sql_query(
        f"""SELECT animal_id, "group",
                   MIN(datetime_enter) AS neubau_enter,
                   MAX(datetime_exit) AS neubau_exit
            FROM allocations
            WHERE animal_id IN ({placeholders})
              AND "group" IN (1005, 1006)
              AND datetime_enter <= ?
              AND datetime_exit >= ?
            GROUP BY animal_id, "group"
            ORDER BY animal_id""",
        con,
        params=[*gouna["animal_id"].tolist(), summer_end, summer_start],
    )

    if alloc.empty:
        return pd.DataFrame()

    # Merge: keep gouna animals that are also in Neubau
    # Collapse multiple group allocations to widest window per animal
    alloc_merged = alloc.groupby("animal_id").agg(
        neubau_enter=("neubau_enter", "min"),
        neubau_exit=("neubau_exit", "max"),
        groups=("group", lambda x: "/".join(str(int(g)) for g in sorted(x.unique()))),
    ).reset_index()

    merged = gouna.merge(alloc_merged, on="animal_id", how="inner")

    # Exclude animals already in Tierauswahl for this year
    ta_animals = set(
        tierauswahl[tierauswahl["year"] == year]["animal_id"].astype(int)
    )
    merged["in_tierauswahl"] = merged["animal_id"].isin(ta_animals)
    overlooked = merged[~merged["in_tierauswahl"]].drop(columns=["in_tierauswahl"])

    return overlooked.sort_values("n_readings", ascending=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Find gouna animals in Neubau but not in Tierauswahl")
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--tierauswahl", type=Path, required=True)
    parser.add_argument("--year", type=int, default=None,
                        help="Check a specific year (default: all years 2021-2024)")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(name)-20s │ %(message)s",
        datefmt="%H:%M:%S",
    )

    con = sqlite3.connect(str(args.db))

    # Load Tierauswahl
    xls = pd.ExcelFile(args.tierauswahl)
    for sheet in xls.sheet_names:
        ta = pd.read_excel(xls, sheet_name=sheet)
        if "animal_id" in ta.columns:
            break
    ta = ta[ta["Auswahl"] == "Ja"].copy()
    ta["animal_id"] = pd.to_numeric(ta["animal_id"], errors="coerce").astype(int)
    ta["year"] = (
        pd.to_numeric(ta["Versuchsjahr"], errors="coerce")
        if "Versuchsjahr" in ta.columns
        else pd.to_datetime(ta["datetime_enter"]).dt.year
    )

    years = [args.year] if args.year else [2021, 2022, 2023, 2024]

    for year in years:
        n_ta = (ta["year"] == year).sum()
        log.info("=" * 60)
        log.info("Year %d — Tierauswahl has %d animals", year, n_ta)

        overlooked = find_overlooked(con, ta, year)
        if overlooked.empty:
            log.info("  No overlooked animals found.")
            continue

        log.info("  Found %d animals with gouna + Neubau but NOT in Tierauswahl:",
                 len(overlooked))
        for _, r in overlooked.iterrows():
            log.info(
                "    %d  n=%6d  groups=%-9s  Neubau %s to %s",
                int(r["animal_id"]), int(r["n_readings"]),
                r["groups"], r["neubau_enter"], r["neubau_exit"],
            )

        # Summary by data quality
        good = overlooked[overlooked["n_readings"] >= 25000]
        marginal = overlooked[
            (overlooked["n_readings"] >= 1000) & (overlooked["n_readings"] < 25000)
        ]
        sparse = overlooked[overlooked["n_readings"] < 1000]
        log.info("  Summary: %d good (>25k), %d marginal (1k-25k), %d sparse (<1k)",
                 len(good), len(marginal), len(sparse))

    con.close()


if __name__ == "__main__":
    main()

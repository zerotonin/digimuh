#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — run_annual_yield                                      ║
# ║  « orchestrator for the year-round milk-yield analysis »         ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Reads the large-cohort extract CSVs from a results tree         ║
# ║  (default results/broken_stick), runs the three annual-yield     ║
# ║  analyses, and writes CSVs + figures into 07_annual_yield.       ║
# ╚══════════════════════════════════════════════════════════════════╝
"""CLI: year-round milk-yield analysis (z-scores, surfaces, duration)."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from digimuh.console import kv, section, setup_logging
from digimuh.paths import resolve_input, resolve_output
from digimuh.stats_annual_yield import (
    aggregate_per_cow,
    build_median_surface,
    compute_annual_zscores,
    compute_heatstress_duration_deltas,
    summarise_duration_deltas,
    test_duration_medians,
)
from digimuh.stats_lactation_curve import load_calvings
from digimuh.viz_annual_yield import plot_heatstress_duration, plot_yield_heatmap

log = logging.getLogger("digimuh.annual_yield")

_SURFACES = (
    ("lactation_nr", "z_yield",        "z",  "lactation", "heatmap_yield_z_lactation"),
    ("lactation_nr", "daily_yield_kg", "kg", "lactation", "heatmap_yield_kg_lactation"),
    ("dim",          "z_yield",        "z",  "dim",       "heatmap_yield_z_dim"),
    ("dim",          "daily_yield_kg", "kg", "dim",       "heatmap_yield_kg_dim"),
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Year-round milk-yield analysis (large Neubau cohort).")
    parser.add_argument("--data", type=str, default="results/broken_stick",
                        help="Results tree holding the extract CSVs.")
    parser.add_argument("--db", type=str, default=None,
                        help="DigiMuh SQLite DB (only to cache calvings.csv "
                             "if it is missing).")
    args = parser.parse_args()
    setup_logging()
    data_dir = Path(args.data)
    db_path = Path(args.db) if args.db else None

    # ── load ──────────────────────────────────────────────
    section("Annual milk-yield analysis", str(data_dir))
    daily_full = pd.read_csv(resolve_input(data_dir, "daily_milk_yield_full.csv"),
                             parse_dates=["date"])
    calvings = load_calvings(data_dir, db_path=db_path)
    kv("daily cow-days (full history)", len(daily_full))
    kv("cohort animals", daily_full["animal_id"].nunique())

    # ── analysis 1: z-scores ──────────────────────────────
    section("1 · per-cow-year z-scores")
    z_df = compute_annual_zscores(daily_full, calvings)
    z_df.to_csv(resolve_output(data_dir, "annual_zscores.csv"), index=False)
    kv("z-scored cow-days", len(z_df))
    kv("cow-years", z_df.groupby(["animal_id", "year"]).ngroups)

    # ── analysis 2: median surfaces ───────────────────────
    section("2 · median-cow yield surfaces")
    for y_col, value_col, kind, y_kind, name in _SURFACES:
        grid = build_median_surface(z_df, y_col=y_col, value_col=value_col)
        plot_yield_heatmap(grid, value_kind=kind, y_kind=y_kind,
                           name=name, out_dir=data_dir)
        kv(name, f"{grid['median_value'].notna().sum()} cells")

    # ── analysis 3: heat-stress duration response ─────────
    section("3 · heat-stress duration response (summer)")
    yield_wood = pd.read_csv(
        resolve_input(data_dir, "daily_milk_yield_wood.csv"),
        parse_dates=["date"])
    breakpoints = pd.read_csv(
        resolve_input(data_dir, "broken_stick_results.csv"))
    climate_daily = pd.read_csv(
        resolve_input(data_dir, "climate_daily.csv"))

    deltas = compute_heatstress_duration_deltas(
        yield_wood, breakpoints, climate_daily)
    if deltas.empty:
        log.warning("no heat-stress runs — skipping duration figures")
        return
    deltas.to_csv(
        resolve_output(data_dir, "heatstress_duration_deltas.csv"), index=False)

    # Pooled median (every run-day weighted equally) — reference.
    pooled = summarise_duration_deltas(deltas)
    pooled.to_csv(
        resolve_output(data_dir, "heatstress_duration_pooled.csv"), index=False)

    # Per-cow → inter-cow median (every cow weighted equally) — headline.
    per_cow = aggregate_per_cow(deltas)
    per_cow.to_csv(
        resolve_output(data_dir, "heatstress_duration_per_cow.csv"), index=False)
    summary = summarise_duration_deltas(per_cow)

    # Wilcoxon signed-rank (median ≠ 0) per streak day, BH-FDR corrected.
    tests = test_duration_medians(per_cow)
    tests.to_csv(
        resolve_output(data_dir, "heatstress_duration_tests.csv"), index=False)

    for tag, name in (("residual", "heatstress_duration_residual"),
                      ("kg", "heatstress_duration_kg")):
        plot_heatstress_duration(
            summary, per_cow, value_tag=tag, name=name, out_dir=data_dir,
            tests=tests[tests["metric"] == tag])

    kv("run-day deltas", len(deltas))
    kv("cow-year runs",
       deltas.groupby(["animal_id", "year", "run_id"]).ngroups)
    kv("cows in duration analysis", per_cow["animal_id"].nunique())
    for _, r in tests[tests["metric"] == "residual"].iterrows():
        kv(f"day {int(r['streak_day'])} (residual)",
           f"med={r['median']:+.2f}  p_fdr={r['p_fdr']:.1e} {r['stars']}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Three-cohort ICC sensitivity sweep for §3.13 of the methods doc.

Runs the breakpoint repeatability analysis on three nested cohorts so
the paper can show that the negative-repeatability finding does not
depend on which cohort filter we adopt:

  1. ``strict``    — Tierauswahl.xlsx (Gundi's tight selection,
                      220 cow-years / 181 animals).
  2. ``extended``  — Tierauswahl_extended.xlsx (current default,
                      255 cow-years / 200 animals).
  3. ``neubau``    — every (animal_id, summer) pair from the DB
                      ``allocations`` table where the animal spent
                      ≥ 30 days in the Neubau experimental groups
                      (1005 / 1006) during Jun-Sep — ~ 771 cow-years
                      / 480 animals before convergence filtering.

For each cohort the script

  * extracts ``rumen_barn`` for that cohort,
  * fits broken-sticks (Frontiers configuration),
  * computes ``breakpoint_icc_<cohort>.csv``,
  * renders ``breakpoint_icc_<cohort>_forest.{svg,png}``.

Outputs go to ``results/broken_stick/`` under the existing subject
subfolders (02_breakpoints, 06_longitudinal); cohort-specific files
carry a ``_<cohort>`` suffix.
"""

from __future__ import annotations

import logging
import os
import sqlite3
from pathlib import Path

# Headless backend before any matplotlib import
os.environ.setdefault("MPLBACKEND", "Agg")

import pandas as pd

from digimuh.config import load_config
from digimuh.extract import (
    extract_production,
    extract_rumen_barn,
    load_tierauswahl,
)
from digimuh.paths import resolve_output
from digimuh.stats_core import run_broken_stick_fits
from digimuh.stats_lactation_curve import attach_dim, load_calvings
from digimuh.stats_longitudinal import compute_breakpoint_icc
from digimuh.viz_base import setup_figure
from digimuh.viz_longitudinal import plot_breakpoint_icc

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(name)s | %(message)s")
log = logging.getLogger("digimuh.icc_sweep")


# ── Configuration ────────────────────────────────────────────────
DB_PATH        = Path(load_config().database)
DATA_DIR       = Path("results/broken_stick")
TIERAUSWAHL    = Path("Tierauswahl.xlsx")
TIERAUSWAHL_X  = Path("Tierauswahl_extended.xlsx")
NEUBAU_GROUPS  = (1005, 1006)
SUMMER_MONTHS  = (6, 7, 8, 9)
MIN_DAYS_IN_GROUP = 30


# ── Helper: build a synthetic Tierauswahl frame from allocations ──

def neubau_pen_cohort(con: sqlite3.Connection) -> pd.DataFrame:
    """Build a Tierauswahl-shaped frame from `allocations`.

    For every animal × summer where the animal spent ≥ 30 days in
    groups 1005 or 1006 during June-September, emit one row with
    ``datetime_enter`` = first day in pen that summer and
    ``datetime_exit`` = last day in pen that summer.  The rumen
    extraction will further restrict via SQL_RUMEN's timestamp
    bounds.
    """
    qmarks = ",".join("?" * len(NEUBAU_GROUPS))
    alloc = pd.read_sql(
        f"""SELECT animal_id, "group" AS grp, datetime_enter, datetime_exit
             FROM allocations
             WHERE "group" IN ({qmarks})
               AND datetime_enter <= '2024-09-30'
               AND datetime_exit  >= '2021-04-01'""",
        con, params=list(NEUBAU_GROUPS),
    )
    alloc["datetime_enter"] = pd.to_datetime(alloc["datetime_enter"])
    alloc["datetime_exit"]  = pd.to_datetime(alloc["datetime_exit"])

    rows: list[dict] = []
    for year in (2021, 2022, 2023, 2024):
        s_start = pd.Timestamp(f"{year}-06-01")
        s_end   = pd.Timestamp(f"{year}-09-30")
        for aid, grp in alloc.groupby("animal_id"):
            covers = []
            for _, r in grp.iterrows():
                lo = max(r["datetime_enter"], s_start)
                hi = min(r["datetime_exit"],  s_end)
                if hi > lo:
                    covers.append((lo, hi))
            if not covers:
                continue
            days = sum((hi - lo).days for lo, hi in covers)
            if days < MIN_DAYS_IN_GROUP:
                continue
            enter = min(lo for lo, _ in covers)
            exit_ = max(hi for _, hi in covers)
            rows.append({
                "animal_id": int(aid),
                "datetime_enter": enter,
                "datetime_exit":  exit_,
                "year": int(year),
                # "Auswahl"/"group" left out — not used downstream.
            })
    df = pd.DataFrame(rows)
    log.info("Neubau-pen cohort: %d cow-years / %d animals (≥%dd in groups %s, summer)",
             len(df), df.animal_id.nunique(), MIN_DAYS_IN_GROUP, NEUBAU_GROUPS)
    return df


# ── DIM + parity lookup (independent of Wood fits) ───────────────

def compute_dim_parity_frame(rumen_barn: pd.DataFrame,
                              calvings: pd.DataFrame) -> pd.DataFrame:
    """Per-cow-day ``(animal_id, year, dim, lactation_nr)`` frame.

    Only requires the rumen_barn extraction and the calvings table
    (no milk-yield data, no Wood fits).  The §3.13 residual-mode
    ICC needs only this — Wood is for milk-yield stratification
    in §3.15, not here.
    """
    if rumen_barn.empty or calvings.empty:
        return pd.DataFrame(columns=["animal_id", "year", "dim", "lactation_nr"])
    ts = pd.to_datetime(rumen_barn["timestamp"], utc=True).dt.tz_convert(None)
    cow_days = (rumen_barn
                .assign(date=ts.dt.normalize())
                [["animal_id", "year", "date"]]
                .drop_duplicates()
                .reset_index(drop=True))
    # ``attach_dim`` expects a daily_yield-shaped frame; the yield
    # column is unused for our purposes but the API requires it.
    cow_days["daily_yield_kg"] = 0.0
    return attach_dim(cow_days, calvings)


# ── Per-cohort runner ────────────────────────────────────────────

def run_cohort(con: sqlite3.Connection, name: str, ta: pd.DataFrame,
               calvings: pd.DataFrame, data_dir: Path) -> None:
    log.info("=" * 70)
    log.info("Cohort '%s': %d cow-years / %d animals",
             name, len(ta), ta.animal_id.nunique())

    # 1. Extract rumen + barn for this cohort.
    log.info("Extracting rumen + barn climate ...")
    rumen = extract_rumen_barn(con, ta)
    log.info("  rumen_barn rows: %d", len(rumen))
    if rumen.empty:
        log.warning("Cohort '%s' produced no rumen rows — skipping", name)
        return

    # 2. Fit broken-sticks (Frontiers configuration: BS only).
    log.info("Fitting broken-sticks (Frontiers mode, no respiration)...")
    bs = run_broken_stick_fits(rumen, pd.DataFrame(), frontiers_only=True)

    # 2b. Attach DIM + parity per cow-year (calvings-only, no Wood).
    dim_frame = compute_dim_parity_frame(rumen, calvings)
    if not dim_frame.empty:
        per_cow_year = (dim_frame
                        .groupby(["animal_id", "year"], as_index=False)
                        .agg(lactation_nr=("lactation_nr",
                                           lambda s: int(s.mode().iloc[0]))))
        bs = bs.merge(per_cow_year, on=["animal_id", "year"], how="left")
        log.info("  attached lactation_nr to %d / %d cow-years",
                 bs["lactation_nr"].notna().sum(), len(bs))

    # 2c. Extract milk yield (mean kg/d during summer window).
    log.info("Extracting milk yield per cow-year ...")
    prod = extract_production(con, ta)
    if not prod.empty:
        prod = prod[["animal_id", "year", "mean_milk_yield_kg"]].copy()
        log.info("  yield records: %d (%d with yield)",
                 len(prod), prod["mean_milk_yield_kg"].notna().sum())
    bs_out = resolve_output(data_dir, f"broken_stick_results_{name}.csv")
    bs.to_csv(bs_out, index=False)
    log.info("  broken-stick rows: %d  →  %s", len(bs), bs_out)

    # 3. ICC — raw + measurement-corrected (SE-cohort) +
    #         parity+DIM residual + parity+DIM+yield residual +
    #         each residual stacked with measurement correction.
    log.info("Computing breakpoint ICC ...")
    icc_df = compute_breakpoint_icc(bs, wood_yield=dim_frame, production=prod)
    icc_out = resolve_output(data_dir, f"breakpoint_icc_{name}.csv")
    icc_df.to_csv(icc_out, index=False)
    log.info("  ICC rows: %d  →  %s", len(icc_df), icc_out)
    log.info("\n%s", icc_df[["predictor", "mode", "n_animals", "n_obs",
                              "icc", "ci_lower", "ci_upper", "p"]]
                       .round(3).to_string(index=False))

    # 4. Forest plot.
    log.info("Rendering forest ...")
    title_suffix_map = {
        "strict":   "Tierauswahl strict (Gundi's selection)",
        "extended": "Tierauswahl extended",
        "neubau":   f"all animals in Neubau pens {NEUBAU_GROUPS}",
    }
    plot_breakpoint_icc(
        data_dir,
        csv_stem=f"breakpoint_icc_{name}",
        figure_stem=f"breakpoint_icc_{name}_forest",
        title_suffix=title_suffix_map.get(name, name),
    )
    log.info("  forest saved to 06_longitudinal/breakpoint_icc_%s_forest.{svg,png}", name)


# ── Main ─────────────────────────────────────────────────────────

def main() -> None:
    setup_figure()
    if not DB_PATH.exists():
        raise SystemExit(
            f"DB not found at {DB_PATH} — set DIGIMUH_DB or the `database` "
            "key in your config to point at cow.db")
    con = sqlite3.connect(DB_PATH)

    # Calvings: needed for DIM + lactation_nr (residual ICC mode).
    # Pulls every calving_confirmation event in the DB once and
    # reuses across all three cohorts.
    calvings = load_calvings(DATA_DIR / "01_extract", db_path=DB_PATH)
    log.info("Calving events loaded: %d", len(calvings))

    # Cohort 1: Gundi's strict Tierauswahl.
    ta_strict = load_tierauswahl(TIERAUSWAHL)
    run_cohort(con, "strict", ta_strict, calvings, DATA_DIR)

    # Cohort 2: Extended Tierauswahl (current default).
    ta_extended = load_tierauswahl(TIERAUSWAHL_X)
    run_cohort(con, "extended", ta_extended, calvings, DATA_DIR)

    # Cohort 3: Synthetic full Neubau-pen cohort.
    ta_neubau = neubau_pen_cohort(con)
    run_cohort(con, "neubau", ta_neubau, calvings, DATA_DIR)

    con.close()
    log.info("All three cohort sweeps complete.")


if __name__ == "__main__":
    main()

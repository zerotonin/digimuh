#!/usr/bin/env python3
"""Table 1 — yearly climate × cohort × RRT × milk-yield summary.

For each summer (1 June – 30 September) of 2021–2024 plus an "Overall"
row, compute mean ± SD of:

* Barn temperature (°C) — Neubau sensors only (barn_id ∈ {1, 2})
* Relative humidity (%) — Neubau sensors only
* Temperature–Humidity Index — Neubau sensors only
* Number of cows included (Neubau-pen cohort, ≥ 30 d in groups 1005/1006
  during the summer window — same definition as the manuscript cohort)
* Reticulorumen temperature (RRT, °C) — `temp_without_drink_cycles` from
  smaxtec_derived, restricted to the Neubau cohort animals × that summer,
  filtered to 30–43 °C, milking-hours excluded (04–07, 16–19)
* Daily milk yield (kg/d) — Neubau-cohort animals × that summer

Writes Markdown to ``results/broken_stick/01_extract/table1_yearly_climate.md``.
"""

from __future__ import annotations

import logging
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.config import load_config

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s | %(message)s")
log = logging.getLogger("digimuh.table1")


DB_PATH       = Path(load_config().database)
DATA_DIR      = Path("results/broken_stick")
OUT_PATH      = DATA_DIR / "01_extract" / "table1_yearly_climate.md"
NEUBAU_GROUPS = (1005, 1006)
NEUBAU_BARNS  = (1, 2)
SUMMER_START  = "06-01"
SUMMER_END    = "09-30"
YEARS         = (2021, 2022, 2023, 2024)


# ── Cohort builder (re-uses the manuscript cohort) ───────────────

def neubau_pen_cohort(con: sqlite3.Connection,
                       min_days: int = 30) -> pd.DataFrame:
    """Per-(animal, year) cohort: ≥ ``min_days`` d in groups 1005/1006
    during 1 June – 30 September of each study year."""
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
    for year in YEARS:
        s_start = pd.Timestamp(f"{year}-{SUMMER_START}")
        s_end   = pd.Timestamp(f"{year}-{SUMMER_END}")
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
            if days < min_days:
                continue
            enter = min(lo for lo, _ in covers)
            exit_ = max(hi for _, hi in covers)
            rows.append({
                "animal_id": int(aid), "year": int(year),
                "datetime_enter": enter, "datetime_exit": exit_,
            })
    return pd.DataFrame(rows)


# ── Climate per summer (Neubau sensors) ──────────────────────────

def climate_summary(con: sqlite3.Connection) -> dict[int, dict]:
    """Mean ± SD barn temp, RH, THI per summer.  Reads at native
    sensor resolution from smaxtec_barns, filtered to Neubau."""
    qmarks = ",".join("?" * len(NEUBAU_BARNS))
    out: dict[int, dict] = {}
    all_T, all_H, all_TH = [], [], []
    for year in YEARS:
        s_start = f"{year}-{SUMMER_START}"
        s_end   = f"{year}-{SUMMER_END} 23:59:59"
        df = pd.read_sql(
            f"""SELECT temp, hum, temp_hum_index
                 FROM smaxtec_barns
                 WHERE barn_id IN ({qmarks})
                   AND "timestamp" >= ? AND "timestamp" <= ?
                   AND temp IS NOT NULL
                   AND hum  IS NOT NULL
                   AND temp_hum_index IS NOT NULL""",
            con, params=[*NEUBAU_BARNS, s_start, s_end],
        )
        out[year] = {
            "barn_T_mean":  df["temp"].mean(),
            "barn_T_sd":    df["temp"].std(),
            "barn_H_mean":  df["hum"].mean(),
            "barn_H_sd":    df["hum"].std(),
            "barn_THI_mean": df["temp_hum_index"].mean(),
            "barn_THI_sd":   df["temp_hum_index"].std(),
            "n_climate":    len(df),
        }
        all_T.append(df["temp"]); all_H.append(df["hum"]); all_TH.append(df["temp_hum_index"])
    overall = {
        "barn_T_mean":  pd.concat(all_T).mean(),
        "barn_T_sd":    pd.concat(all_T).std(),
        "barn_H_mean":  pd.concat(all_H).mean(),
        "barn_H_sd":    pd.concat(all_H).std(),
        "barn_THI_mean": pd.concat(all_TH).mean(),
        "barn_THI_sd":   pd.concat(all_TH).std(),
        "n_climate":    sum(len(s) for s in all_T),
    }
    out["overall"] = overall
    return out


# ── RRT per summer (Neubau cohort × milking-hours-excluded) ──────

def rrt_summary(con: sqlite3.Connection,
                cohort: pd.DataFrame) -> dict[int, dict]:
    """Mean ± SD reticulorumen temperature per summer, restricted to
    the Neubau cohort animals × each summer window.  Filtered to
    30–43 °C; milking hours excluded (04–07, 16–19) to match the
    manuscript pipeline (§3.1 of analysis_00_methods.md).  Reads in
    chunks to keep memory bounded."""
    out: dict[int, dict] = {}
    grand_n = 0
    grand_sum = 0.0
    grand_sumsq = 0.0
    for year in YEARS:
        ids = cohort.loc[cohort["year"] == year, "animal_id"].astype(int).unique()
        if len(ids) == 0:
            out[year] = {"rrt_mean": np.nan, "rrt_sd": np.nan, "n_rrt": 0}
            continue
        s_start = f"{year}-{SUMMER_START}"
        s_end   = f"{year}-{SUMMER_END} 23:59:59"
        qmarks = ",".join("?" * len(ids))
        log.info(f"  RRT query for {year} ({len(ids)} animals) ...")
        df = pd.read_sql(
            f"""SELECT "timestamp",
                       CAST(temp_without_drink_cycles AS REAL) AS rrt
                 FROM smaxtec_derived
                 WHERE animal_id IN ({qmarks})
                   AND "timestamp" >= ? AND "timestamp" <= ?
                   AND temp_without_drink_cycles IS NOT NULL
                   AND CAST(temp_without_drink_cycles AS REAL)
                       BETWEEN 30 AND 43""",
            con, params=[*[int(x) for x in ids], s_start, s_end],
        )
        # Milking-hour exclusion (UTC hours 4-7 and 16-19, matching
        # extract_rumen_barn).
        ts = pd.to_datetime(df["timestamp"], utc=True)
        hr = ts.dt.hour
        df = df[~hr.between(4, 7) & ~hr.between(16, 19)]
        rrt = df["rrt"].dropna().to_numpy()
        out[year] = {
            "rrt_mean": float(np.mean(rrt)) if rrt.size else np.nan,
            "rrt_sd":   float(np.std(rrt, ddof=1)) if rrt.size > 1 else np.nan,
            "n_rrt":    int(rrt.size),
        }
        # Online accumulator for overall mean/sd (avoid concatenating
        # the year-arrays of millions of rows).
        grand_n += rrt.size
        grand_sum += float(np.sum(rrt))
        grand_sumsq += float(np.sum(rrt ** 2))
    if grand_n > 1:
        mean_o = grand_sum / grand_n
        var_o  = (grand_sumsq - grand_n * mean_o ** 2) / (grand_n - 1)
        out["overall"] = {"rrt_mean": mean_o,
                          "rrt_sd":   float(np.sqrt(max(var_o, 0.0))),
                          "n_rrt":    grand_n}
    else:
        out["overall"] = {"rrt_mean": np.nan, "rrt_sd": np.nan, "n_rrt": 0}
    return out


# ── Daily milk yield per summer (Neubau cohort) ──────────────────

def yield_summary(con: sqlite3.Connection,
                   cohort: pd.DataFrame) -> dict[int, dict]:
    """Mean ± SD of per-cow-day daily milk yield per summer.  One
    record per (animal, date) — daily totals from herdeplus."""
    out: dict[int, dict] = {}
    all_y: list[pd.Series] = []
    for year in YEARS:
        ids = cohort.loc[cohort["year"] == year, "animal_id"].astype(int).unique()
        if len(ids) == 0:
            out[year] = {"y_mean": np.nan, "y_sd": np.nan, "n_y": 0}
            continue
        s_start = f"{year}-{SUMMER_START}"
        s_end   = f"{year}-{SUMMER_END} 23:59:59"
        qmarks = ",".join("?" * len(ids))
        df = pd.read_sql(
            f"""SELECT animal_id,
                       DATE("timestamp") AS date,
                       SUM(herdeplus_milked_mkg) AS daily_yield_kg
                 FROM herdeplus
                 WHERE animal_id IN ({qmarks})
                   AND "timestamp" >= ? AND "timestamp" <= ?
                   AND herdeplus_milked_mkg IS NOT NULL
                   AND herdeplus_milked_mkg > 0
                 GROUP BY animal_id, DATE("timestamp")""",
            con, params=[*[int(x) for x in ids], s_start, s_end],
        )
        y = df["daily_yield_kg"].dropna()
        out[year] = {
            "y_mean": float(y.mean()) if not y.empty else np.nan,
            "y_sd":   float(y.std(ddof=1)) if len(y) > 1 else np.nan,
            "n_y":    int(len(y)),
        }
        all_y.append(y)
    big = pd.concat(all_y) if all_y else pd.Series(dtype=float)
    out["overall"] = {
        "y_mean": float(big.mean()) if not big.empty else np.nan,
        "y_sd":   float(big.std(ddof=1)) if len(big) > 1 else np.nan,
        "n_y":    int(len(big)),
    }
    return out


# ── Markdown writer ──────────────────────────────────────────────

def fmt(mean, sd, decimals=1):
    if not np.isfinite(mean) or not np.isfinite(sd):
        return "—"
    return f"{mean:.{decimals}f} ± {sd:.{decimals}f}"


def render_markdown(climate, rrt, yield_, cohort) -> str:
    n_cows = {y: int(cohort.loc[cohort.year == y, "animal_id"].nunique())
              for y in YEARS}
    n_cows["overall"] = int(cohort["animal_id"].nunique())

    rows = []
    for key in [*YEARS, "overall"]:
        c = climate[key]
        r = rrt[key]
        ye = yield_[key]
        label = "Overall" if key == "overall" else str(key)
        rows.append({
            "Year":               label,
            "Barn T (°C)":        fmt(c["barn_T_mean"], c["barn_T_sd"]),
            "RH (%)":             fmt(c["barn_H_mean"], c["barn_H_sd"]),
            "THI":                fmt(c["barn_THI_mean"], c["barn_THI_sd"]),
            "n cows":             str(n_cows[key]),
            "RRT (°C)":           fmt(r["rrt_mean"], r["rrt_sd"], 2),
            "Daily yield (kg)":   fmt(ye["y_mean"], ye["y_sd"]),
        })
    df = pd.DataFrame(rows)

    md = []
    md.append("# Table 1 — Yearly climate, cohort size, RRT, and milk yield")
    md.append("")
    md.append("Summer period 1 June – 30 September each year.  Cohort = Neubau pen cohort")
    md.append("(animals with ≥ 30 days in groups 1005 or 1006 during the summer window).")
    md.append("Climate from smaXtec barn sensors restricted to the Neubau barns")
    md.append("(`barn_id ∈ {1, 2}`).  RRT from smaxtec_derived `temp_without_drink_cycles`")
    md.append("filtered to 30–43 °C with milking hours (04–07, 16–19 UTC) excluded.")
    md.append("Daily milk yield is the per-(animal, date) sum of `herdeplus_milked_mkg`.")
    md.append("Mean ± SD shown.")
    md.append("")
    md.append("| Year | Barn temperature (°C) | Relative air humidity (%) | "
              "Temperature-humidity index | Number of cows included | "
              "RRT (°C) | Daily milk yield (kg) |")
    md.append("|---|---|---|---|---:|---|---|")
    for _, r in df.iterrows():
        md.append(f"| {r['Year']} | {r['Barn T (°C)']} | {r['RH (%)']} | "
                  f"{r['THI']} | {r['n cows']} | {r['RRT (°C)']} | "
                  f"{r['Daily yield (kg)']} |")
    md.append("")
    clim_n = ", ".join(f"{y}: {climate[y]['n_climate']:,}" for y in YEARS)
    rrt_n  = ", ".join(f"{y}: {rrt[y]['n_rrt']:,}" for y in YEARS)
    yld_n  = ", ".join(f"{y}: {yield_[y]['n_y']:,}" for y in YEARS)
    md.append(
        f"*Underlying counts: climate readings n = {clim_n} "
        f"(overall {climate['overall']['n_climate']:,}); "
        f"RRT readings n = {rrt_n} "
        f"(overall {rrt['overall']['n_rrt']:,}); "
        f"daily-yield records n = {yld_n} "
        f"(overall {yield_['overall']['n_y']:,}).*")
    md.append("")
    return "\n".join(md)


def main() -> None:
    if not DB_PATH.exists():
        log.error("Database not found at %s — set DIGIMUH_DB or the "
                  "`database` key in your config to point at cow.db", DB_PATH)
        sys.exit(1)
    con = sqlite3.connect(DB_PATH)

    log.info("Building Neubau cohort ...")
    cohort = neubau_pen_cohort(con)
    log.info("  cohort: %d cow-years across %d animals",
             len(cohort), cohort.animal_id.nunique())

    log.info("Aggregating barn climate ...")
    clim = climate_summary(con)

    log.info("Aggregating reticulorumen temperature ...")
    rrt = rrt_summary(con, cohort)

    log.info("Aggregating daily milk yield ...")
    ye = yield_summary(con, cohort)

    md = render_markdown(clim, rrt, ye, cohort)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(md)
    log.info("Wrote %s", OUT_PATH)
    print()
    print(md)


if __name__ == "__main__":
    main()

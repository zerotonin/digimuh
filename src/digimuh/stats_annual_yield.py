#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_annual_yield                                    ║
# ║  « year-round yield: z-scores, surfaces, heat-stress runs »      ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Full-calendar-year milk-yield analysis on the large Neubau      ║
# ║  cohort (480 cows, 2021-2024).  Three products:                  ║
# ║                                                                  ║
# ║    1. per-cow-year z-scored daily yield                          ║
# ║    2. median-cow yield surfaces over (day-of-year x              ║
# ║       lactation) and (day-of-year x DIM)                         ║
# ║    3. yield delta as a function of consecutive heat-stress       ║
# ║       days, on the Frontiers summer breakpoints                  ║
# ║                                                                  ║
# ║  Reuses the Wood/DIM machinery in stats_lactation_curve.         ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Year-round milk-yield statistics for the large Neubau cohort."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from digimuh.stats_lactation_curve import DIM_MAX, DIM_MIN, attach_dim

log = logging.getLogger("digimuh.annual_yield")


# ─────────────────────────────────────────────────────────────
#  « Defaults »
# ─────────────────────────────────────────────────────────────

MIN_DAYS_PER_COW_YEAR: int = 30
"""Fewer milking days than this in a calendar year and the per-cow-year
z-score (mean/SD) is too unstable to trust — that cow-year is dropped."""

DOY_BIN_DAYS: int = 7
"""Day-of-year bin width for the median surfaces (weekly)."""

DIM_BIN_DAYS: int = 15
"""DIM bin width for the (day-of-year × DIM) surface."""

LACTATION_CAP: int = 8
"""Lactation numbers ≥ this are pooled into one top bin (sparse tail)."""

MIN_COWS_PER_CELL: int = 3
"""Surface cells backed by fewer cows than this are left empty (NaN)."""

BASELINE_LOOKBACK_DAYS: int = 3
"""How far back to search for a valid pre-run baseline yield when the
calendar day immediately before a heat-stress run has no milking record."""


# ─────────────────────────────────────────────────────────────
#  « Analysis 1 — per-cow-year z-scores »
# ─────────────────────────────────────────────────────────────

def compute_annual_zscores(
    daily_full: pd.DataFrame,
    calvings: pd.DataFrame,
    dim_min: int = DIM_MIN,
    dim_max: int = DIM_MAX,
    min_days: int = MIN_DAYS_PER_COW_YEAR,
) -> pd.DataFrame:
    """Z-score each cow's daily yield within each calendar year.

    For every (animal, calendar-year) the daily yield is standardised
    to that cow-year's own mean and SD, so cows of different production
    level and lactation stage become comparable before pooling.  DIM and
    lactation number are attached for the downstream surfaces.

    Args:
        daily_full: Full-history daily yield (``animal_id``, ``date``,
            ``daily_yield_kg``).
        calvings: From :func:`stats_lactation_curve.load_calvings`.
        dim_min, dim_max: DIM clipping (defaults 5–305).
        min_days: Minimum milking days in a cow-year to keep it.

    Returns:
        DataFrame with ``animal_id``, ``year``, ``date``, ``doy``,
        ``dim``, ``lactation_nr``, ``daily_yield_kg``, ``z_yield``.
    """
    df = attach_dim(daily_full, calvings, dim_min=dim_min, dim_max=dim_max)
    df["year"] = df["date"].dt.year.astype(int)
    df["doy"] = df["date"].dt.dayofyear.astype(int)

    grp = df.groupby(["animal_id", "year"])["daily_yield_kg"]
    mu = grp.transform("mean")
    sd = grp.transform("std")
    n = grp.transform("size")
    df["z_yield"] = (df["daily_yield_kg"] - mu) / sd.replace(0.0, np.nan)
    df = df[n >= min_days].copy()

    cols = ["animal_id", "year", "date", "doy", "dim",
            "lactation_nr", "daily_yield_kg", "z_yield"]
    log.info("z-scored %d cow-days across %d cow-years",
             len(df), df.groupby(["animal_id", "year"]).ngroups)
    return df[cols].reset_index(drop=True)


# ─────────────────────────────────────────────────────────────
#  « Analysis 2 — median-cow yield surfaces »
# ─────────────────────────────────────────────────────────────

def build_median_surface(
    z_df: pd.DataFrame,
    y_col: str,
    value_col: str,
    doy_bin: int = DOY_BIN_DAYS,
    dim_bin: int = DIM_BIN_DAYS,
    min_cows: int = MIN_COWS_PER_CELL,
    n_boot: int = 2000,
    seed: int = 12345,
) -> pd.DataFrame:
    """Median-across-cows yield per (day-of-year bin × y bin) cell.

    Each cow contributes once per cell (her mean ``value_col`` in that
    cell), so cows with many days do not dominate; the cell value is the
    median over those per-cow means — "the median cow's yield".  A
    bootstrap 95% CI of that median (resampling over cows) is attached so
    the same surface can be drawn as a heatmap or as CI-banded line plots.

    Args:
        z_df: Output of :func:`compute_annual_zscores`.
        y_col: ``"lactation_nr"`` or ``"dim"``.
        value_col: ``"z_yield"`` or ``"daily_yield_kg"``.
        doy_bin: Day-of-year bin width.
        dim_bin: DIM bin width (only used when ``y_col == "dim"``).
        min_cows: Cells with fewer backing cows are set NaN.
        n_boot: Bootstrap resamples for the per-cell median CI.
        seed: RNG seed (reproducible).

    Returns:
        Long-form DataFrame: ``doy_bin`` (left edge), ``y_bin``,
        ``median_value``, ``ci_lo``, ``ci_hi``, ``n_cows``, ``n_obs``.
    """
    d = z_df.dropna(subset=[value_col, y_col]).copy()
    d["doy_bin"] = ((d["doy"] - 1) // doy_bin) * doy_bin + 1

    if y_col == "lactation_nr":
        d["y_bin"] = d["lactation_nr"].clip(upper=LACTATION_CAP).astype(int)
    elif y_col == "dim":
        d["y_bin"] = (d["dim"] // dim_bin) * dim_bin
    else:
        raise ValueError(f"unsupported y_col {y_col!r}")

    # 1. per-cow mean within each cell, 2. median across cows (+ boot CI).
    per_cow = (d.groupby(["doy_bin", "y_bin", "animal_id"])[value_col]
                 .mean().reset_index())
    rng = np.random.default_rng(seed)
    records: list[dict] = []
    for (doy_b, y_b), g in per_cow.groupby(["doy_bin", "y_bin"]):
        vals = g[value_col].to_numpy()
        lo, hi = _median_ci(vals, n_boot, rng)
        records.append({"doy_bin": int(doy_b), "y_bin": int(y_b),
                        "median_value": float(np.median(vals)),
                        "ci_lo": lo, "ci_hi": hi,
                        "n_cows": int(g["animal_id"].nunique())})
    grid = pd.DataFrame.from_records(records)
    n_obs = (d.groupby(["doy_bin", "y_bin"]).size()
               .rename("n_obs").reset_index())
    grid = grid.merge(n_obs, on=["doy_bin", "y_bin"], how="left")
    grid.loc[grid["n_cows"] < min_cows,
             ["median_value", "ci_lo", "ci_hi"]] = np.nan
    return grid.sort_values(["y_bin", "doy_bin"]).reset_index(drop=True)


def _median_ci(
    values: np.ndarray,
    n_boot: int,
    rng: np.random.Generator,
    alpha: float = 0.05,
) -> tuple[float, float]:
    """Bootstrap ``(1-alpha)`` CI of the median over ``values``."""
    x = np.asarray(values, dtype=float)
    x = x[~np.isnan(x)]
    if x.size == 0:
        return (np.nan, np.nan)
    if x.size == 1:
        return (float(x[0]), float(x[0]))
    boot = np.median(x[rng.integers(0, x.size, size=(n_boot, x.size))], axis=1)
    return (float(np.percentile(boot, 100 * alpha / 2)),
            float(np.percentile(boot, 100 * (1 - alpha / 2))))


# ─────────────────────────────────────────────────────────────
#  « Analysis 3 — yield vs consecutive heat-stress days »
# ─────────────────────────────────────────────────────────────

def compute_heatstress_duration_deltas(
    yield_wood: pd.DataFrame,
    breakpoints: pd.DataFrame,
    climate_daily: pd.DataFrame,
    lookback: int = BASELINE_LOOKBACK_DAYS,
) -> pd.DataFrame:
    """Yield change across runs of consecutive heat-stress days.

    A day is *heat-stress* for a cow when the barn's daily-maximum THI
    reaches or exceeds that cow's individual (converged) THI breakpoint.
    Within each cow-year we find maximal runs of consecutive calendar
    days that are all heat-stress, take the last non-stress day before
    the run as the baseline, and record the yield change on the 1st,
    2nd, 3rd … day of the run.

    The primary response is the Wood residual (DIM-adjusted), so the
    lactation decline over the spell is not mistaken for a heat effect;
    raw kg is carried alongside.

    Args:
        yield_wood: Summer Wood frame (``animal_id``, ``date``, ``year``,
            ``dim``, ``daily_yield_kg``, ``yield_residual``).
        breakpoints: ``broken_stick_results`` with ``animal_id``,
            ``year``, ``thi_breakpoint``, ``thi_converged``.
        climate_daily: ``day``, ``barn_thi_max`` (herd-level barn climate).
        lookback: Max days to search back for a valid baseline yield.

    Returns:
        Long DataFrame, one row per (cow-year run-day with valid yield):
        ``animal_id``, ``year``, ``run_id``, ``streak_day``, ``date``,
        ``dim``, ``daily_yield_kg``, ``yield_residual``,
        ``baseline_residual``, ``baseline_kg``, ``delta_residual``,
        ``delta_kg``, ``run_length``.
    """
    bp = breakpoints.loc[
        breakpoints["thi_converged"].astype(bool),
        ["animal_id", "year", "thi_breakpoint"]].copy()

    clim = climate_daily[["day", "barn_thi_max"]].copy()
    clim["day"] = pd.to_datetime(clim["day"]).dt.normalize()

    yw = yield_wood[["animal_id", "year", "date", "dim",
                     "daily_yield_kg", "yield_residual"]].copy()
    yw["date"] = pd.to_datetime(yw["date"]).dt.normalize()

    records: list[dict] = []
    for (aid, yr), thr in (
        bp.set_index(["animal_id", "year"])["thi_breakpoint"].items()
    ):
        cow = yw[(yw["animal_id"] == aid) & (yw["year"] == yr)]
        if cow.empty:
            continue
        records.extend(_cow_year_runs(aid, yr, float(thr), cow, clim, lookback))

    out = pd.DataFrame.from_records(records)
    if out.empty:
        log.warning("no heat-stress runs found")
    else:
        log.info("%d run-day deltas across %d cow-year runs",
                 len(out), out.groupby(["animal_id", "year", "run_id"]).ngroups)
    return out


def _cow_year_runs(aid, yr, threshold, cow, clim, lookback) -> list[dict]:
    """Detect heat-stress runs for one cow-year and emit run-day deltas."""
    # Daily calendar spanning this cow-year's observed window.
    start, end = cow["date"].min(), cow["date"].max()
    cal = pd.DataFrame({"date": pd.date_range(start, end, freq="D")})
    cal = cal.merge(clim, left_on="date", right_on="day", how="left")
    cal = cal.merge(
        cow[["date", "dim", "daily_yield_kg", "yield_residual"]],
        on="date", how="left")
    cal["heat"] = cal["barn_thi_max"] >= threshold
    cal = cal.reset_index(drop=True)

    # Run id increments at every False→True / True→False boundary.
    run_break = cal["heat"].ne(cal["heat"].shift()).cumsum()
    out: list[dict] = []
    run_counter = 0
    for _, seg in cal.groupby(run_break):
        if not bool(seg["heat"].iloc[0]):
            continue
        run_counter += 1
        run_len = len(seg)
        first_idx = seg.index[0]
        baseline = _baseline_before(cal, first_idx, lookback)
        if baseline is None:
            continue
        base_res, base_kg = baseline
        for pos, (_, day) in enumerate(seg.iterrows(), start=1):
            if pd.isna(day["yield_residual"]):
                continue
            out.append(dict(
                animal_id=int(aid), year=int(yr), run_id=run_counter,
                streak_day=pos, date=day["date"], dim=day["dim"],
                daily_yield_kg=day["daily_yield_kg"],
                yield_residual=day["yield_residual"],
                baseline_residual=base_res, baseline_kg=base_kg,
                delta_residual=day["yield_residual"] - base_res,
                delta_kg=day["daily_yield_kg"] - base_kg,
                run_length=run_len,
            ))
    return out


def _baseline_before(cal, first_idx, lookback) -> tuple[float, float] | None:
    """Last non-heat day with a valid yield within ``lookback`` of a run."""
    for back in range(1, lookback + 1):
        j = first_idx - back
        if j < 0:
            break
        row = cal.iloc[j]
        if not bool(row["heat"]) and pd.notna(row["yield_residual"]):
            return float(row["yield_residual"]), float(row["daily_yield_kg"])
    return None


# ─────────────────────────────────────────────────────────────
#  « Duration-response summary + bootstrap CI »
# ─────────────────────────────────────────────────────────────

def aggregate_per_cow(
    deltas: pd.DataFrame,
    max_streak: int = 7,
) -> pd.DataFrame:
    """Collapse run-day deltas to one median per cow per streak-day.

    Stage one of the per-cow → inter-cow estimator: each cow contributes
    her own median delta at each streak position, so a cow with many
    heat-stress runs does not dominate the herd-level summary.  Feed the
    result to :func:`summarise_duration_deltas` to obtain the inter-cow
    median (median across cows) with a cow-level bootstrap CI.

    Args:
        deltas: Output of :func:`compute_heatstress_duration_deltas`.
        max_streak: Keep streak days 1…max_streak.

    Returns:
        DataFrame: ``animal_id``, ``streak_day``, ``delta_residual``,
        ``delta_kg`` (each the cow's median at that streak day), ``n_obs``.
    """
    d = deltas[deltas["streak_day"] <= max_streak]
    return (d.groupby(["animal_id", "streak_day"])
             .agg(delta_residual=("delta_residual", "median"),
                  delta_kg=("delta_kg", "median"),
                  n_obs=("delta_residual", "size"))
             .reset_index())


def summarise_duration_deltas(
    deltas: pd.DataFrame,
    max_streak: int = 7,
    n_boot: int = 20000,
    seed: int = 12345,
) -> pd.DataFrame:
    """Median yield delta per streak-day with bootstrap 95% CI.

    Operates on whatever rows it is given: pass the raw run-day deltas
    for the pooled median, or the :func:`aggregate_per_cow` table for the
    inter-cow median (median of per-cow medians, bootstrap over cows).

    Args:
        deltas: Run-day deltas, or per-cow medians from
            :func:`aggregate_per_cow`.
        max_streak: Report streak days 1…max_streak (deeper runs pooled
            out — the sample thins quickly).
        n_boot: Bootstrap resamples for the median CI.
        seed: RNG seed (reproducible).

    Returns:
        DataFrame: ``streak_day``, ``n``, ``n_cows``, ``median_delta_residual``,
        ``ci_lo_residual``, ``ci_hi_residual``, ``median_delta_kg``,
        ``ci_lo_kg``, ``ci_hi_kg``.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for k in range(1, max_streak + 1):
        sub = deltas[deltas["streak_day"] == k]
        if sub.empty:
            continue
        row = {"streak_day": k, "n": len(sub),
               "n_cows": sub["animal_id"].nunique()}
        for tag, col in (("residual", "delta_residual"), ("kg", "delta_kg")):
            x = sub[col].dropna().to_numpy()
            med = float(np.median(x))
            boot = np.median(x[rng.integers(0, len(x), size=(n_boot, len(x)))],
                             axis=1)
            row[f"median_delta_{tag}"] = med
            row[f"ci_lo_{tag}"] = float(np.percentile(boot, 2.5))
            row[f"ci_hi_{tag}"] = float(np.percentile(boot, 97.5))
        rows.append(row)
    return pd.DataFrame(rows)


def _pvalue_stars(p: float) -> str:
    """Standard significance stars (n.s. / * / ** / ***)."""
    if not np.isfinite(p):
        return "n.s."
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def test_duration_medians(
    per_cow: pd.DataFrame,
    max_streak: int = 7,
) -> pd.DataFrame:
    """Wilcoxon signed-rank test that each streak-day median ≠ 0, BH-FDR.

    One test per streak day, on the per-cow median deltas (one value per
    cow — independent observations), asking whether the centre of that
    per-cow distribution differs from zero.  Raw p-values are
    Benjamini-Hochberg FDR corrected *within each response metric*
    (the seven streak-day tests form the family).

    Args:
        per_cow: Output of :func:`aggregate_per_cow` (per-cow medians).
        max_streak: Number of streak days tested (the test family size).

    Returns:
        DataFrame: ``metric`` ("residual"/"kg"), ``streak_day``,
        ``n_cows``, ``median``, ``wilcoxon_stat``, ``p_raw``, ``p_fdr``,
        ``stars``.
    """
    from rerandomstats import correct_pvalues_array
    from scipy.stats import wilcoxon

    rows: list[dict] = []
    for tag, col in (("residual", "delta_residual"), ("kg", "delta_kg")):
        block: list[dict] = []
        for k in range(1, max_streak + 1):
            x = per_cow.loc[per_cow["streak_day"] == k, col].dropna().to_numpy()
            if x.size == 0:
                continue
            try:
                stat, p = wilcoxon(x)
            except ValueError:
                # All values zero → no evidence against H0.
                stat, p = np.nan, 1.0
            block.append(dict(metric=tag, streak_day=k, n_cows=int(x.size),
                              median=float(np.median(x)),
                              wilcoxon_stat=float(stat), p_raw=float(p)))
        adj = correct_pvalues_array(
            np.array([b["p_raw"] for b in block]), method="fdr_bh")
        for b, pa in zip(block, adj):
            b["p_fdr"] = float(pa)
            b["stars"] = _pvalue_stars(float(pa))
            rows.append(b)
    return pd.DataFrame(rows)

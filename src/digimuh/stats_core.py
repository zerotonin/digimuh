#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_core                                           ║
# ║  « FDR correction, broken-stick fits, correlations, tests »    ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Core statistical pipeline: run broken-stick fits across all    ║
# ║  animals, compute Spearman correlations, below/above bp means, ║
# ║  and Wilcoxon tests with BH-FDR correction.                    ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Core statistical functions for the broken-stick analysis pipeline."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from digimuh.constants import THI_RANGE, TEMP_RANGE, MIN_READINGS

log = logging.getLogger("digimuh.stats")

def benjamini_hochberg(p_values: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """Benjamini-Hochberg FDR correction.

    Args:
        p_values: Array of raw p-values.
        alpha: Target FDR level (for reference; returns adjusted p).

    Returns:
        Array of BH-adjusted p-values.
    """
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return p

    # Handle NaNs
    notnan = ~np.isnan(p)
    p_valid = p[notnan]
    m = len(p_valid)
    if m == 0:
        return p

    order = np.argsort(p_valid)
    ranked = np.empty(m)
    ranked[order] = np.arange(1, m + 1)

    adjusted = p_valid * m / ranked
    # Enforce monotonicity (step-down)
    adjusted[order] = np.minimum.accumulate(adjusted[order[::-1]])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    result = np.full(n, np.nan)
    result[notnan] = adjusted
    return result


def p_to_stars(p: float) -> str:
    """Convert p-value to significance stars.

    Returns:
        ``'***'`` if p < 0.001, ``'**'`` if p < 0.01,
        ``'*'`` if p < 0.05, ``'n.s.'`` otherwise.
    """
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


# ─────────────────────────────────────────────────────────────
#  « broken-stick fits for all animals »
# ─────────────────────────────────────────────────────────────

def run_broken_stick_fits(
    rumen: pd.DataFrame, resp: pd.DataFrame,
    frontiers_only: bool = False,
) -> pd.DataFrame:
    """Fit broken-stick models per animal-year.

    Args:
        rumen: rumen_barn.csv DataFrame.
        resp: respiration_barn.csv DataFrame.
        frontiers_only: If True, skip Davies/pscore/Hill fits
            (Frontiers paper: broken-stick only).  These methods
            are reserved for the COMPAG companion papers.

    Returns:
        One row per animal-year with breakpoint results.
    """
    results = []
    groups = rumen.groupby(["animal_id", "year", "date_enter"])
    total = len(groups)

    try:
        from tqdm import tqdm
        iterator = tqdm(groups, desc="  Fitting animals", total=total,
                        bar_format="  {l_bar}{bar:30}{r_bar}")
    except ImportError:
        iterator = groups
        log.info("  Fitting %d animal-years …", total)

    for (aid, year, enter), grp in iterator:

        n_readings = len(grp)
        if n_readings < 50:
            results.append({
                "animal_id": aid, "year": year, "date_enter": enter,
                "n_readings": n_readings,
                "thi_breakpoint": np.nan, "thi_converged": False,
                "temp_breakpoint": np.nan, "temp_converged": False,
                "resp_thi_breakpoint": np.nan, "resp_thi_converged": False,
                "resp_temp_breakpoint": np.nan, "resp_temp_converged": False,
                "comment": f"insufficient data ({n_readings})",
            })
            continue

        # Body temp fits
        thi_fit = broken_stick_fit(
            grp["barn_thi"].values, grp["body_temp"].values, x_range=(45, 80))
        temp_fit = broken_stick_fit(
            grp["barn_temp"].values, grp["body_temp"].values, x_range=(5, 35))

        # Breakpoint existence tests (body temp) — COMPAG companion papers
        if not frontiers_only:
            thi_davies = davies_test(
                grp["barn_thi"].values, grp["body_temp"].values, x_range=(45, 80))
            thi_pscore = pscore_test(
                grp["barn_thi"].values, grp["body_temp"].values, x_range=(45, 80))
            temp_davies = davies_test(
                grp["barn_temp"].values, grp["body_temp"].values, x_range=(5, 35))
            temp_pscore = pscore_test(
                grp["barn_temp"].values, grp["body_temp"].values, x_range=(5, 35))
            thi_hill = hill_fit(
                grp["barn_thi"].values, grp["body_temp"].values, x_range=(45, 80))
            temp_hill = hill_fit(
                grp["barn_temp"].values, grp["body_temp"].values, x_range=(5, 35))
        else:
            _empty_test = {"p_value": np.nan}
            thi_davies = thi_pscore = temp_davies = temp_pscore = _empty_test
            _empty_hill = {"ec50": np.nan, "hill_n": np.nan,
                           "lower_bend": np.nan, "r_squared": np.nan,
                           "aic": np.nan, "converged": False}
            thi_hill = temp_hill = _empty_hill

        # Respiration fits
        resp_grp = resp[
            (resp["animal_id"] == aid) & (resp["year"] == year)
        ] if not resp.empty else pd.DataFrame()
        has_resp = len(resp_grp) >= 50

        if has_resp:
            resp_thi_fit = broken_stick_fit(
                resp_grp["barn_thi"].values, resp_grp["resp_rate"].values,
                x_range=(45, 80))
            resp_temp_fit = broken_stick_fit(
                resp_grp["barn_temp"].values, resp_grp["resp_rate"].values,
                x_range=(5, 35))
            # Breakpoint existence tests (respiration) — COMPAG
            if not frontiers_only:
                resp_thi_davies = davies_test(
                    resp_grp["barn_thi"].values, resp_grp["resp_rate"].values,
                    x_range=(45, 80))
                resp_thi_pscore = pscore_test(
                    resp_grp["barn_thi"].values, resp_grp["resp_rate"].values,
                    x_range=(45, 80))
                resp_temp_davies = davies_test(
                    resp_grp["barn_temp"].values, resp_grp["resp_rate"].values,
                    x_range=(5, 35))
                resp_temp_pscore = pscore_test(
                    resp_grp["barn_temp"].values, resp_grp["resp_rate"].values,
                    x_range=(5, 35))
                resp_thi_hill = hill_fit(
                    resp_grp["barn_thi"].values, resp_grp["resp_rate"].values,
                    x_range=(45, 80))
                resp_temp_hill = hill_fit(
                    resp_grp["barn_temp"].values, resp_grp["resp_rate"].values,
                    x_range=(5, 35))
            else:
                _empty_test = {"p_value": np.nan}
                resp_thi_davies = resp_thi_pscore = _empty_test
                resp_temp_davies = resp_temp_pscore = _empty_test
                _empty_hill = {"ec50": np.nan, "hill_n": np.nan,
                               "lower_bend": np.nan, "r_squared": np.nan,
                               "aic": np.nan, "converged": False}
                resp_thi_hill = resp_temp_hill = _empty_hill
        else:
            resp_thi_fit = {"breakpoint": np.nan, "converged": False, "n": 0}
            resp_temp_fit = {"breakpoint": np.nan, "converged": False, "n": 0}
            _empty_test = {"p_value": np.nan}
            resp_thi_davies = resp_thi_pscore = _empty_test
            resp_temp_davies = resp_temp_pscore = _empty_test
            _empty_hill = {"ec50": np.nan, "hill_n": np.nan,
                           "lower_bend": np.nan, "r_squared": np.nan,
                           "aic": np.nan, "converged": False}
            resp_thi_hill = resp_temp_hill = _empty_hill

        rec = {
            "animal_id": aid, "year": int(year), "date_enter": enter,
            "n_readings": n_readings,
            "n_resp_readings": len(resp_grp),
            # Body temp vs THI
            "thi_breakpoint": thi_fit["breakpoint"],
            "thi_slope_below": thi_fit.get("slope_below"),
            "thi_slope_above": thi_fit.get("slope_above"),
            "thi_r_squared": thi_fit.get("r_squared"),
            "thi_converged": thi_fit["converged"],
            "thi_davies_p": thi_davies["p_value"],
            "thi_pscore_p": thi_pscore["p_value"],
            "thi_hill_ec50": thi_hill.get("ec50"),
            "thi_hill_n": thi_hill.get("hill_n"),
            "thi_hill_bend": thi_hill.get("lower_bend"),
            "thi_hill_r2": thi_hill.get("r_squared"),
            "thi_hill_converged": thi_hill.get("converged", False),
            # Body temp vs barn temp
            "temp_breakpoint": temp_fit["breakpoint"],
            "temp_slope_below": temp_fit.get("slope_below"),
            "temp_slope_above": temp_fit.get("slope_above"),
            "temp_r_squared": temp_fit.get("r_squared"),
            "temp_converged": temp_fit["converged"],
            "temp_davies_p": temp_davies["p_value"],
            "temp_pscore_p": temp_pscore["p_value"],
            "temp_hill_ec50": temp_hill.get("ec50"),
            "temp_hill_n": temp_hill.get("hill_n"),
            "temp_hill_bend": temp_hill.get("lower_bend"),
            "temp_hill_r2": temp_hill.get("r_squared"),
            "temp_hill_converged": temp_hill.get("converged", False),
            # Resp vs THI
            "resp_thi_breakpoint": resp_thi_fit["breakpoint"],
            "resp_thi_slope_below": resp_thi_fit.get("slope_below"),
            "resp_thi_slope_above": resp_thi_fit.get("slope_above"),
            "resp_thi_r_squared": resp_thi_fit.get("r_squared"),
            "resp_thi_converged": resp_thi_fit.get("converged", False),
            "resp_thi_davies_p": resp_thi_davies["p_value"],
            "resp_thi_pscore_p": resp_thi_pscore["p_value"],
            "resp_thi_hill_ec50": resp_thi_hill.get("ec50"),
            "resp_thi_hill_n": resp_thi_hill.get("hill_n"),
            "resp_thi_hill_bend": resp_thi_hill.get("lower_bend"),
            "resp_thi_hill_r2": resp_thi_hill.get("r_squared"),
            "resp_thi_hill_converged": resp_thi_hill.get("converged", False),
            # Resp vs barn temp
            "resp_temp_breakpoint": resp_temp_fit["breakpoint"],
            "resp_temp_r_squared": resp_temp_fit.get("r_squared"),
            "resp_temp_converged": resp_temp_fit.get("converged", False),
            "resp_temp_davies_p": resp_temp_davies["p_value"],
            "resp_temp_pscore_p": resp_temp_pscore["p_value"],
            "resp_temp_hill_ec50": resp_temp_hill.get("ec50"),
            "resp_temp_hill_n": resp_temp_hill.get("hill_n"),
            "resp_temp_hill_bend": resp_temp_hill.get("lower_bend"),
            "resp_temp_hill_r2": resp_temp_hill.get("r_squared"),
            "resp_temp_hill_converged": resp_temp_hill.get("converged", False),
        }
        results.append(rec)

    return pd.DataFrame(results)


# ─────────────────────────────────────────────────────────────
#  « Spearman correlations »
# ─────────────────────────────────────────────────────────────

def compute_spearman(rumen: pd.DataFrame, resp: pd.DataFrame) -> pd.DataFrame:
    """Per-animal Spearman correlations."""
    records = []
    for (aid, year), grp in rumen.groupby(["animal_id", "year"]):
        rec = {"animal_id": aid, "year": int(year), "n": len(grp)}
        for x_col, prefix in [("barn_thi", "thi"), ("barn_temp", "temp")]:
            sub = grp.dropna(subset=[x_col, "body_temp"])
            if len(sub) > 20:
                rs, p = spearmanr(sub[x_col], sub["body_temp"])
                rec[f"{prefix}_rs"] = rs
                rec[f"{prefix}_p"] = p
            else:
                rec[f"{prefix}_rs"] = np.nan
                rec[f"{prefix}_p"] = np.nan
        # Respiration Spearman
        resp_grp = resp[(resp["animal_id"] == aid) & (resp["year"] == year)] if not resp.empty else pd.DataFrame()
        for x_col, prefix in [("barn_thi", "resp_thi"), ("barn_temp", "resp_temp")]:
            sub = resp_grp.dropna(subset=[x_col, "resp_rate"]) if not resp_grp.empty else pd.DataFrame()
            if len(sub) > 20:
                rs, p = spearmanr(sub[x_col], sub["resp_rate"])
                rec[f"{prefix}_rs"] = rs
                rec[f"{prefix}_p"] = p
            else:
                rec[f"{prefix}_rs"] = np.nan
                rec[f"{prefix}_p"] = np.nan
        records.append(rec)
    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────
#  « below/above breakpoint comparisons »
# ─────────────────────────────────────────────────────────────

def compute_below_above(
    rumen: pd.DataFrame, resp: pd.DataFrame,
    bs_results: pd.DataFrame,
) -> pd.DataFrame:
    """Per-animal means below/above their individual THI breakpoint."""
    converged = bs_results[bs_results["thi_converged"] == True]
    records = []

    for _, row in converged.iterrows():
        aid = int(row["animal_id"])
        year = int(row["year"])
        bp = row["thi_breakpoint"]

        grp = rumen[(rumen["animal_id"] == aid) & (rumen["year"] == year)]
        if len(grp) < 30:
            continue

        below = grp[grp["barn_thi"] <= bp]
        above = grp[grp["barn_thi"] > bp]
        if len(below) < 10 or len(above) < 10:
            continue

        rec = {
            "animal_id": aid, "year": year,
            "thi_breakpoint": bp,
            "body_temp_below": below["body_temp"].mean(),
            "body_temp_above": above["body_temp"].mean(),
            "n_below": len(below),
            "n_above": len(above),
        }

        # Respiration
        resp_grp = resp[(resp["animal_id"] == aid) & (resp["year"] == year)] if not resp.empty else pd.DataFrame()
        if len(resp_grp) >= 20:
            rb = resp_grp[resp_grp["barn_thi"] <= bp]
            ra = resp_grp[resp_grp["barn_thi"] > bp]
            if len(rb) >= 5 and len(ra) >= 5:
                rec["resp_below"] = rb["resp_rate"].mean()
                rec["resp_above"] = ra["resp_rate"].mean()
            else:
                rec["resp_below"] = rec["resp_above"] = np.nan
        else:
            rec["resp_below"] = rec["resp_above"] = np.nan

        records.append(rec)

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────
#  « within-year Fisher resampling tests with BH-FDR »
# ─────────────────────────────────────────────────────────────

def run_statistical_tests(beh: pd.DataFrame) -> pd.DataFrame:
    """Run Fisher resampling tests within each year, BH-FDR corrected.

    Uses reRandomStats FisherResamplingTest with medianDiff as the test
    statistic and 20,000 permutations.

    Tests per year:
    - Body temp below vs above breakpoint

    Args:
        beh: behavioural_response DataFrame.

    Returns:
        DataFrame: one row per test, with raw p, adjusted p, stars.
    """
    from rerandomstats import FisherResamplingTest

    test_rows = []
    years = sorted(beh["year"].dropna().unique().astype(int))

    for year in years:
        yr = beh[beh["year"] == year]
        year_tests = []

        # Body temperature below vs above breakpoint
        paired = yr.dropna(subset=["body_temp_below", "body_temp_above"])
        if len(paired) >= 10:
            p = FisherResamplingTest(
                data_a=paired["body_temp_below"].tolist(),
                data_b=paired["body_temp_above"].tolist(),
                func="medianDiff",
                combination_n=20_000,
            ).main()
            year_tests.append({
                "year": year,
                "test": "body_temp below vs above (Fisher medianDiff)",
                "n": len(paired),
                "p_raw": p,
                "median_below": paired["body_temp_below"].median(),
                "median_above": paired["body_temp_above"].median(),
                "median_diff": (paired["body_temp_above"] - paired["body_temp_below"]).median(),
            })

        # Respiration (if columns exist and have data)
        if "resp_below" in yr.columns:
            paired_r = yr.dropna(subset=["resp_below", "resp_above"])
            if len(paired_r) >= 10:
                p = FisherResamplingTest(
                    data_a=paired_r["resp_below"].tolist(),
                    data_b=paired_r["resp_above"].tolist(),
                    func="medianDiff",
                    combination_n=20_000,
                ).main()
                year_tests.append({
                    "year": year,
                    "test": "respiration below vs above (Fisher medianDiff)",
                    "n": len(paired_r),
                    "p_raw": p,
                    "median_below": paired_r["resp_below"].median(),
                    "median_above": paired_r["resp_above"].median(),
                    "median_diff": (paired_r["resp_above"] - paired_r["resp_below"]).median(),
                })

        # BH-FDR within this year
        if year_tests:
            raw_ps = np.array([t["p_raw"] for t in year_tests])
            adj_ps = benjamini_hochberg(raw_ps)
            for t, adj_p in zip(year_tests, adj_ps):
                t["p_adj"] = adj_p
                t["stars"] = p_to_stars(adj_p)
            test_rows.extend(year_tests)

    return pd.DataFrame(test_rows)


# ─────────────────────────────────────────────────────────────
#  « cross-correlation / cross-covariance below/above bp »

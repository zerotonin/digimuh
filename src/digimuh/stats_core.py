#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_core                                            ║
# ║  « broken-stick fits, correlations, below/above tests »          ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Core statistical orchestration: per-animal broken-stick /       ║
# ║  Davies / pscore / Hill fits, Spearman correlations, below /     ║
# ║  above-breakpoint means, and Fisher resampling tests with        ║
# ║  BH-FDR correction.  The four fitters and the FDR routine        ║
# ║  themselves live in reRandomStats (>= 0.2.0); this module        ║
# ║  contains only the DigiMuh-specific wiring.                      ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Core statistical functions for the broken-stick analysis pipeline."""

from __future__ import annotations

import logging
from collections.abc import Sequence

import numpy as np
import pandas as pd
from rerandomstats import (
    broken_stick_fit,
    correct_pvalues_array,
    davies_test,
    hill_fit,
    pscore_test,
)
from scipy.stats import spearmanr

from digimuh.constants import BREAKPOINT_SEARCH_WINDOW, RESAMPLING_SEED
from digimuh.stats_breakpoint_summary import check_breakpoint_continuity

log = logging.getLogger("digimuh.stats")


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


def apply_fdr_within(
    df: pd.DataFrame,
    family_cols: Sequence[str] = (),
    *,
    p_col: str = "p",
    out_col: str = "p_fdr",
) -> pd.DataFrame:
    """Add Benjamini-Hochberg adjusted p-values, corrected within families.

    A screen that reports dozens of correlations and stars each one on its
    raw p-value will manufacture significance; this is the correction that
    stops it.  ``p_col`` is deliberately left untouched so its meaning never
    changes underneath an existing reader — the adjusted value lands in
    ``out_col``, and the family each test belonged to is recorded in
    ``family`` so the correction can be audited from the CSV alone.

    Args:
        df:          Long results table, one row per test.
        family_cols: Columns whose combination defines a family.  Empty means
                     the whole table is a single family.
        p_col:       Column holding the raw p-values.
        out_col:     Column to write the adjusted p-values into.

    Returns:
        A copy with ``out_col`` and ``family`` added.  NaN p-values are
        preserved in place and excluded from the family size.
    """
    if df.empty or p_col not in df.columns:
        return df
    out = df.copy()
    keys = [c for c in family_cols if c in out.columns]
    out["family"] = (out[keys].astype(str).agg(" | ".join, axis=1)
                     if keys else "all")
    out[out_col] = np.nan
    for _, idx in out.groupby("family").groups.items():
        raw = out.loc[idx, p_col].to_numpy(dtype=float)
        out.loc[idx, out_col] = correct_pvalues_array(raw, method="fdr_bh")
    return out


# ─────────────────────────────────────────────────────────────
#  « broken-stick fits for all animals »
# ─────────────────────────────────────────────────────────────

def flag_unreliable_breakpoints(bs: pd.DataFrame) -> pd.DataFrame:
    """Add ``<predictor>_bp_reliable`` flags marking unidentified fits.

    A converged fit is *unreliable* when its breakpoint is not
    identifiable — the CI is wider than :data:`BP_CI_WIDTH_FRAC` of the
    window, the knee sits within :data:`BP_EDGE_FRAC` of the lower edge, or
    the sub-threshold slope is more negative than
    :data:`BP_SLOPE_BELOW_MIN`.  Non-converged fits are reliable = False.
    Nothing is dropped; consumers filter on the flag.

    A CI *truncated at the upper boundary* is deliberately not a flag: it
    marks a high, well-estimated threshold (a heat-tolerant cow whose knee
    sits near the top of the range), not an unidentified fit — the CI-width
    rule already catches genuinely wide, non-identifiable CIs in either
    direction.
    """
    from digimuh.constants import (
        BP_CI_WIDTH_FRAC,
        BP_EDGE_FRAC,
        BP_SLOPE_BELOW_MIN,
        BREAKPOINT_SEARCH_WINDOW,
    )

    out = bs.copy()
    for pred, (lo, hi) in BREAKPOINT_SEARCH_WINDOW.items():
        bp_col = f"{pred}_breakpoint"
        if bp_col not in out.columns:
            continue
        w = hi - lo
        idx = out.index
        conv = out.get(f"{pred}_converged", pd.Series(False, index=idx)).fillna(False) == True  # noqa: E712
        bp = out[bp_col]
        ci_lo = out.get(f"{pred}_breakpoint_ci_lo")
        ci_hi = out.get(f"{pred}_breakpoint_ci_hi")
        wide = ((ci_hi - ci_lo) > BP_CI_WIDTH_FRAC * w
                if ci_lo is not None and ci_hi is not None
                else pd.Series(False, index=idx))
        lower_edge = bp <= lo + BP_EDGE_FRAC * w
        neg_slope = out.get(f"{pred}_slope_below",
                            pd.Series(np.nan, index=idx)) < BP_SLOPE_BELOW_MIN
        unreliable = (wide.fillna(False) | lower_edge.fillna(False)
                      | neg_slope.fillna(False))
        out[f"{pred}_bp_reliable"] = conv & ~unreliable
    return out


def _linear_r2(x: np.ndarray, y: np.ndarray) -> float:
    """R² of a straight-line fit — the smooth reaction-norm baseline.

    Recorded per animal-year so the segmented model can be compared by
    AIC against a no-threshold linear response without re-fitting.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3 or np.ptp(x[mask]) == 0 or np.ptp(y[mask]) == 0:
        return np.nan
    r = np.corrcoef(x[mask], y[mask])[0, 1]
    return float(r * r)


def _constraint_bound(free_fit: dict) -> bool | float:
    """Did the slope constraint bind?  True when the free optimum violates it."""
    if not free_fit.get("converged"):
        return np.nan
    return not (free_fit["slope_above"] > free_fit["slope_below"]
                and free_fit["slope_above"] > 0)


def run_broken_stick_fits(
    rumen: pd.DataFrame, resp: pd.DataFrame,
) -> pd.DataFrame:
    """Fit broken-stick models per animal-year.

    Every cow-summer gets the constrained broken-stick fit, the
    unconstrained sensitivity fit, the Davies and pseudo-score tests for
    the existence of a breakpoint, and the Hill fit — unconditionally,
    so the decision of what to report is made downstream, not here.  The
    returned frame is checked for continuity at every breakpoint.

    Args:
        rumen: rumen_barn.csv DataFrame.
        resp: respiration_barn.csv DataFrame.

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
                "thi_linear_r2": np.nan, "temp_linear_r2": np.nan,
                "temp_breakpoint": np.nan, "temp_converged": False,
                "resp_thi_breakpoint": np.nan, "resp_thi_converged": False,
                "resp_temp_breakpoint": np.nan, "resp_temp_converged": False,
                "comment": f"insufficient data ({n_readings})",
            })
            continue

        # Body temp fits: constrained, unconstrained sensitivity, existence
        # tests and the Hill alternative — all of them, every cow-summer
        thi_x, temp_x, y = (grp["barn_thi"].values, grp["barn_temp"].values,
                            grp["body_temp"].values)
        thi_win, temp_win = (BREAKPOINT_SEARCH_WINDOW["thi"],
                             BREAKPOINT_SEARCH_WINDOW["temp"])
        thi_fit = broken_stick_fit(thi_x, y, x_range=thi_win)
        temp_fit = broken_stick_fit(temp_x, y, x_range=temp_win)
        thi_free = broken_stick_fit(thi_x, y, x_range=thi_win, constrain=False)
        temp_free = broken_stick_fit(temp_x, y, x_range=temp_win, constrain=False)
        thi_davies = davies_test(thi_x, y, x_range=thi_win)
        thi_pscore = pscore_test(thi_x, y, x_range=thi_win)
        temp_davies = davies_test(temp_x, y, x_range=temp_win)
        temp_pscore = pscore_test(temp_x, y, x_range=temp_win)
        thi_hill = hill_fit(thi_x, y, x_range=thi_win)
        temp_hill = hill_fit(temp_x, y, x_range=temp_win)

        # Respiration fits
        resp_grp = resp[
            (resp["animal_id"] == aid) & (resp["year"] == year)
        ] if not resp.empty else pd.DataFrame()
        has_resp = len(resp_grp) >= 50

        if has_resp:
            rthi_x, rtemp_x, ry = (resp_grp["barn_thi"].values,
                                   resp_grp["barn_temp"].values,
                                   resp_grp["resp_rate"].values)
            resp_thi_fit = broken_stick_fit(rthi_x, ry, x_range=thi_win)
            resp_temp_fit = broken_stick_fit(rtemp_x, ry, x_range=temp_win)
            resp_thi_davies = davies_test(rthi_x, ry, x_range=thi_win)
            resp_thi_pscore = pscore_test(rthi_x, ry, x_range=thi_win)
            resp_temp_davies = davies_test(rtemp_x, ry, x_range=temp_win)
            resp_temp_pscore = pscore_test(rtemp_x, ry, x_range=temp_win)
            resp_thi_hill = hill_fit(rthi_x, ry, x_range=thi_win)
            resp_temp_hill = hill_fit(rtemp_x, ry, x_range=temp_win)
        else:
            resp_thi_fit = {"breakpoint": np.nan, "converged": False, "n": 0}
            resp_temp_fit = {"breakpoint": np.nan, "converged": False, "n": 0}
            _empty_test = {"pvalue": np.nan}
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
            "thi_breakpoint_ci_lo": thi_fit.get("breakpoint_ci_lo"),
            "thi_breakpoint_ci_hi": thi_fit.get("breakpoint_ci_hi"),
            "thi_breakpoint_se": thi_fit.get("breakpoint_se"),
            "thi_breakpoint_ci_truncated": thi_fit.get("breakpoint_ci_truncated", False),
            "thi_slope_below": thi_fit.get("slope_below"),
            "thi_slope_above": thi_fit.get("slope_above"),
            "thi_intercept_below": thi_fit.get("intercept_below"),
            "thi_intercept_above": thi_fit.get("intercept_above"),
            "thi_r_squared": thi_fit.get("r_squared"),
            "thi_linear_r2": _linear_r2(grp["barn_thi"].values,
                                        grp["body_temp"].values),
            "thi_converged": thi_fit["converged"],
            "thi_davies_p": thi_davies["pvalue"],
            "thi_pscore_p": thi_pscore["pvalue"],
            "thi_hill_ec50": thi_hill.get("ec50"),
            "thi_hill_n": thi_hill.get("hill_n"),
            "thi_hill_bend": thi_hill.get("lower_bend"),
            "thi_hill_r2": thi_hill.get("r_squared"),
            "thi_hill_converged": thi_hill.get("converged", False),
            "thi_breakpoint_unconstrained": thi_free["breakpoint"],
            "thi_unconstrained_converged": thi_free["converged"],
            "thi_unconstrained_slope_below": thi_free.get("slope_below"),
            "thi_unconstrained_slope_above": thi_free.get("slope_above"),
            "thi_constraint_bound": _constraint_bound(thi_free),
            # Body temp vs barn temp
            "temp_breakpoint": temp_fit["breakpoint"],
            "temp_breakpoint_ci_lo": temp_fit.get("breakpoint_ci_lo"),
            "temp_breakpoint_ci_hi": temp_fit.get("breakpoint_ci_hi"),
            "temp_breakpoint_se": temp_fit.get("breakpoint_se"),
            "temp_breakpoint_ci_truncated": temp_fit.get("breakpoint_ci_truncated", False),
            "temp_slope_below": temp_fit.get("slope_below"),
            "temp_slope_above": temp_fit.get("slope_above"),
            "temp_intercept_below": temp_fit.get("intercept_below"),
            "temp_intercept_above": temp_fit.get("intercept_above"),
            "temp_r_squared": temp_fit.get("r_squared"),
            "temp_linear_r2": _linear_r2(grp["barn_temp"].values,
                                         grp["body_temp"].values),
            "temp_converged": temp_fit["converged"],
            "temp_davies_p": temp_davies["pvalue"],
            "temp_pscore_p": temp_pscore["pvalue"],
            "temp_hill_ec50": temp_hill.get("ec50"),
            "temp_hill_n": temp_hill.get("hill_n"),
            "temp_hill_bend": temp_hill.get("lower_bend"),
            "temp_hill_r2": temp_hill.get("r_squared"),
            "temp_hill_converged": temp_hill.get("converged", False),
            "temp_breakpoint_unconstrained": temp_free["breakpoint"],
            "temp_unconstrained_converged": temp_free["converged"],
            "temp_unconstrained_slope_below": temp_free.get("slope_below"),
            "temp_unconstrained_slope_above": temp_free.get("slope_above"),
            "temp_constraint_bound": _constraint_bound(temp_free),
            # Resp vs THI
            "resp_thi_breakpoint": resp_thi_fit["breakpoint"],
            "resp_thi_slope_below": resp_thi_fit.get("slope_below"),
            "resp_thi_slope_above": resp_thi_fit.get("slope_above"),
            "resp_thi_r_squared": resp_thi_fit.get("r_squared"),
            "resp_thi_converged": resp_thi_fit.get("converged", False),
            "resp_thi_davies_p": resp_thi_davies["pvalue"],
            "resp_thi_pscore_p": resp_thi_pscore["pvalue"],
            "resp_thi_hill_ec50": resp_thi_hill.get("ec50"),
            "resp_thi_hill_n": resp_thi_hill.get("hill_n"),
            "resp_thi_hill_bend": resp_thi_hill.get("lower_bend"),
            "resp_thi_hill_r2": resp_thi_hill.get("r_squared"),
            "resp_thi_hill_converged": resp_thi_hill.get("converged", False),
            # Resp vs barn temp
            "resp_temp_breakpoint": resp_temp_fit["breakpoint"],
            "resp_temp_r_squared": resp_temp_fit.get("r_squared"),
            "resp_temp_converged": resp_temp_fit.get("converged", False),
            "resp_temp_davies_p": resp_temp_davies["pvalue"],
            "resp_temp_pscore_p": resp_temp_pscore["pvalue"],
            "resp_temp_hill_ec50": resp_temp_hill.get("ec50"),
            "resp_temp_hill_n": resp_temp_hill.get("hill_n"),
            "resp_temp_hill_bend": resp_temp_hill.get("lower_bend"),
            "resp_temp_hill_r2": resp_temp_hill.get("r_squared"),
            "resp_temp_hill_converged": resp_temp_hill.get("converged", False),
        }
        results.append(rec)

    out = flag_unreliable_breakpoints(pd.DataFrame(results))
    worst = check_breakpoint_continuity(out)
    log.info("  continuity at the breakpoint — max |jump|: %s",
             ", ".join(f"{k} {v:.1e}" for k, v in worst.items()))
    return out


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
        resp_grp = (
            resp[(resp["animal_id"] == aid) & (resp["year"] == year)]
            if not resp.empty else pd.DataFrame()
        )
        for x_col, prefix in [("barn_thi", "resp_thi"), ("barn_temp", "resp_temp")]:
            sub = (
                resp_grp.dropna(subset=[x_col, "resp_rate"])
                if not resp_grp.empty else pd.DataFrame()
            )
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
    keep_col = ("thi_bp_reliable" if "thi_bp_reliable" in bs_results.columns
                else "thi_converged")
    converged = bs_results[bs_results[keep_col] == True]
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
        resp_grp = (
            resp[(resp["animal_id"] == aid) & (resp["year"] == year)]
            if not resp.empty else pd.DataFrame()
        )
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
                combination_n=20_000, seed=RESAMPLING_SEED).main()
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
                    combination_n=20_000, seed=RESAMPLING_SEED).main()
                year_tests.append({
                    "year": year,
                    "test": "respiration below vs above (Fisher medianDiff)",
                    "n": len(paired_r),
                    "p_raw": p,
                    "median_below": paired_r["resp_below"].median(),
                    "median_above": paired_r["resp_above"].median(),
                    "median_diff": (paired_r["resp_above"] - paired_r["resp_below"]).median(),
                })

        # BH-FDR within this year (delegated to rerandomstats — shared core)
        if year_tests:
            raw_ps = np.array([t["p_raw"] for t in year_tests])
            adj_ps = correct_pvalues_array(raw_ps, method="fdr_bh")
            for t, adj_p in zip(year_tests, adj_ps):
                t["p_adj"] = adj_p
                t["stars"] = p_to_stars(adj_p)
            test_rows.extend(year_tests)

    return pd.DataFrame(test_rows)


# ─────────────────────────────────────────────────────────────
#  « cross-correlation / cross-covariance below/above bp »

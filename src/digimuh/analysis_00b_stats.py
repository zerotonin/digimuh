#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  Analysis 0b — Statistical analysis ¹                        ║
# ║  « broken-stick fits, correlations, BH-FDR tests »          ║
# ╚══════════════════════════════════════════════════════════════╝
"""Reads CSVs from ``analysis_00a_extract`` and runs all statistics.

Outputs (in ``--data`` directory)::

    broken_stick_results.csv   Per-animal breakpoints (rumen + resp)
    spearman_correlations.csv  Per-animal Spearman rₛ
    behavioural_response.csv   Below/above breakpoint means per animal
    statistical_tests.csv      All Wilcoxon tests with BH-FDR correction
    breakpoint_stability.csv   Repeat-animal year pairs + ICC
    summary_table.csv          Table 1 for manuscript

Usage::

    digimuh-stats --data results/broken_stick

¹ Analysis led by Dr. med. vet. Gundula Hoffmann, ATB Potsdam.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr, wilcoxon

log = logging.getLogger("digimuh.stats")


# ─────────────────────────────────────────────────────────────
#  « broken-stick regression model »
# ─────────────────────────────────────────────────────────────

def broken_stick_fit(
    x: np.ndarray, y: np.ndarray,
    x_range: tuple[float, float] | None = None,
    n_grid: int = 200,
) -> dict:
    """Fit a two-segment piecewise linear model with one breakpoint.

    Model: flat-ish below, steeper above.  Biological constraint:
    ``slope_above > slope_below`` and ``slope_above > 0``.

    Args:
        x: Predictor (THI or barn temp).
        y: Response (body temp or resp rate).
        x_range: Search bounds for breakpoint.
        n_grid: Grid resolution for initial search.

    Returns:
        Dict with breakpoint, slopes, R², convergence status.
    """
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)
    if n < 30:
        return {"breakpoint": np.nan, "converged": False, "n": n}

    if x_range is None:
        x_range = (np.percentile(x, 10), np.percentile(x, 90))
    if x_range[0] >= x_range[1]:
        return {"breakpoint": np.nan, "converged": False, "n": n}

    def rss_at_bp(bp):
        left, right = x <= bp, x > bp
        if left.sum() < 10 or right.sum() < 10:
            return np.inf
        try:
            Xl = np.column_stack([np.ones(left.sum()), x[left]])
            cl, *_ = np.linalg.lstsq(Xl, y[left], rcond=None)
            Xr = np.column_stack([np.ones(right.sum()), x[right]])
            cr, *_ = np.linalg.lstsq(Xr, y[right], rcond=None)
        except np.linalg.LinAlgError:
            return np.inf
        return np.sum((y[left] - Xl @ cl) ** 2) + np.sum((y[right] - Xr @ cr) ** 2)

    grid = np.linspace(x_range[0], x_range[1], n_grid)
    rss_vals = np.array([rss_at_bp(bp) for bp in grid])
    best_idx = np.argmin(rss_vals)
    if not np.isfinite(rss_vals[best_idx]):
        return {"breakpoint": np.nan, "converged": False, "n": n}

    lo = grid[max(0, best_idx - 2)]
    hi = grid[min(len(grid) - 1, best_idx + 2)]
    result = minimize_scalar(rss_at_bp, bounds=(lo, hi), method="bounded")
    bp = result.x
    final_rss = result.fun

    left, right = x <= bp, x > bp
    Xl = np.column_stack([np.ones(left.sum()), x[left]])
    cl, *_ = np.linalg.lstsq(Xl, y[left], rcond=None)
    Xr = np.column_stack([np.ones(right.sum()), x[right]])
    cr, *_ = np.linalg.lstsq(Xr, y[right], rcond=None)

    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_sq = 1 - final_rss / ss_tot if ss_tot > 0 else np.nan

    slope_below, slope_above = cl[1], cr[1]
    valid = slope_above > slope_below and slope_above > 0

    return {
        "breakpoint": bp if valid else np.nan,
        "slope_below": slope_below,
        "intercept_below": cl[0],
        "slope_above": slope_above,
        "intercept_above": cr[0],
        "r_squared": r_sq,
        "n": n,
        "n_below": int(left.sum()),
        "n_above": int(right.sum()),
        "converged": valid,
        "rejected_reason": (
            None if valid
            else f"slope_above ({slope_above:.4f}) <= slope_below ({slope_below:.4f})"
        ),
    }


# ─────────────────────────────────────────────────────────────
#  « Benjamini-Hochberg FDR »
# ─────────────────────────────────────────────────────────────

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
) -> pd.DataFrame:
    """Fit broken-stick models per animal-year.

    Args:
        rumen: rumen_barn.csv DataFrame.
        resp: respiration_barn.csv DataFrame.

    Returns:
        One row per animal-year with breakpoint results.
    """
    results = []
    groups = rumen.groupby(["animal_id", "year", "date_enter"])
    total = len(groups)

    for i, ((aid, year, enter), grp) in enumerate(groups):
        if (i + 1) % 20 == 0 or i == 0:
            log.info("  [%d/%d] Fitting animal %d (%s)", i + 1, total, aid, year)

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
        else:
            resp_thi_fit = {"breakpoint": np.nan, "converged": False, "n": 0}
            resp_temp_fit = {"breakpoint": np.nan, "converged": False, "n": 0}

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
            # Body temp vs barn temp
            "temp_breakpoint": temp_fit["breakpoint"],
            "temp_slope_below": temp_fit.get("slope_below"),
            "temp_slope_above": temp_fit.get("slope_above"),
            "temp_r_squared": temp_fit.get("r_squared"),
            "temp_converged": temp_fit["converged"],
            # Resp vs THI
            "resp_thi_breakpoint": resp_thi_fit["breakpoint"],
            "resp_thi_slope_below": resp_thi_fit.get("slope_below"),
            "resp_thi_slope_above": resp_thi_fit.get("slope_above"),
            "resp_thi_r_squared": resp_thi_fit.get("r_squared"),
            "resp_thi_converged": resp_thi_fit.get("converged", False),
            # Resp vs barn temp
            "resp_temp_breakpoint": resp_temp_fit["breakpoint"],
            "resp_temp_r_squared": resp_temp_fit.get("r_squared"),
            "resp_temp_converged": resp_temp_fit.get("converged", False),
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
#  « within-year Wilcoxon tests with BH-FDR »
# ─────────────────────────────────────────────────────────────

def run_statistical_tests(beh: pd.DataFrame) -> pd.DataFrame:
    """Run Wilcoxon signed-rank tests within each year, BH-FDR corrected.

    Tests per year:
    - Body temp below vs above breakpoint
    - Respiration below vs above breakpoint (if available)

    Args:
        beh: behavioural_response DataFrame.

    Returns:
        DataFrame: one row per test, with raw p, adjusted p, stars.
    """
    test_rows = []
    years = sorted(beh["year"].dropna().unique().astype(int))

    for year in years:
        yr = beh[beh["year"] == year]
        year_tests = []

        # Body temperature
        paired = yr.dropna(subset=["body_temp_below", "body_temp_above"])
        if len(paired) >= 10:
            stat, p = wilcoxon(paired["body_temp_below"], paired["body_temp_above"])
            year_tests.append({
                "year": year, "test": "body_temp below vs above",
                "n": len(paired), "statistic": stat, "p_raw": p,
                "median_below": paired["body_temp_below"].median(),
                "median_above": paired["body_temp_above"].median(),
                "median_diff": (paired["body_temp_above"] - paired["body_temp_below"]).median(),
            })

        # Respiration
        paired_r = yr.dropna(subset=["resp_below", "resp_above"])
        if len(paired_r) >= 10:
            stat, p = wilcoxon(paired_r["resp_below"], paired_r["resp_above"])
            year_tests.append({
                "year": year, "test": "respiration below vs above",
                "n": len(paired_r), "statistic": stat, "p_raw": p,
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
#  « breakpoint stability (ICC) »
# ─────────────────────────────────────────────────────────────

def compute_stability(bs_results: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    """ICC and paired data for repeat animals.

    Returns:
        (pairs DataFrame, ICC value).
    """
    conv = bs_results[bs_results["thi_converged"] == True]
    counts = conv.groupby("animal_id").size()
    repeats = counts[counts >= 2].index

    if len(repeats) < 3:
        return pd.DataFrame(), np.nan

    repeat_data = conv[conv["animal_id"].isin(repeats)].sort_values(["animal_id", "year"])
    groups = [g["thi_breakpoint"].values for _, g in repeat_data.groupby("animal_id")]
    all_vals = repeat_data["thi_breakpoint"].dropna()
    grand_mean = all_vals.mean()
    n_animals = len(groups)
    k_mean = len(all_vals) / n_animals

    ss_b = sum(len(g) * (g.mean() - grand_mean) ** 2 for g in groups)
    ss_w = sum(np.sum((g - g.mean()) ** 2) for g in groups)
    ms_b = ss_b / max(n_animals - 1, 1)
    ms_w = ss_w / max(len(all_vals) - n_animals, 1)
    denom = ms_b + (k_mean - 1) * ms_w
    icc = (ms_b - ms_w) / denom if denom > 0 else np.nan

    pairs = []
    for aid, grp in repeat_data.groupby("animal_id"):
        sg = grp.sort_values("year")
        for i in range(len(sg) - 1):
            pairs.append({
                "animal_id": aid,
                "year_1": int(sg.iloc[i]["year"]),
                "year_2": int(sg.iloc[i + 1]["year"]),
                "thi_bp_1": sg.iloc[i]["thi_breakpoint"],
                "thi_bp_2": sg.iloc[i + 1]["thi_breakpoint"],
                "temp_bp_1": sg.iloc[i]["temp_breakpoint"],
                "temp_bp_2": sg.iloc[i + 1]["temp_breakpoint"],
            })

    return pd.DataFrame(pairs), icc


# ─────────────────────────────────────────────────────────────
#  « summary table »
# ─────────────────────────────────────────────────────────────

def make_summary_table(bs: pd.DataFrame) -> pd.DataFrame:
    """Generate Table 1 for the manuscript."""
    years = sorted(bs["year"].dropna().unique().astype(int))
    rows = []
    for year in [*years, "All"]:
        yr = bs if year == "All" else bs[bs["year"] == year]
        for prefix, label in [
            ("thi", "THI→body_temp"),
            ("temp", "barn_temp→body_temp"),
            ("resp_thi", "THI→resp"),
            ("resp_temp", "barn_temp→resp"),
        ]:
            conv = yr[yr[f"{prefix}_converged"] == True]
            bp_col = f"{prefix}_breakpoint"
            rows.append({
                "year": year,
                "model": label,
                "n_total": len(yr),
                "n_converged": len(conv),
                "bp_median": conv[bp_col].median() if len(conv) else np.nan,
                "bp_q25": conv[bp_col].quantile(0.25) if len(conv) else np.nan,
                "bp_q75": conv[bp_col].quantile(0.75) if len(conv) else np.nan,
                "bp_min": conv[bp_col].min() if len(conv) else np.nan,
                "bp_max": conv[bp_col].max() if len(conv) else np.nan,
            })
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Statistical analysis for broken-stick")
    parser.add_argument("--data", type=Path, required=True,
                        help="Directory with CSVs from extraction step")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(name)-20s │ %(message)s",
        datefmt="%H:%M:%S",
    )
    d = args.data

    log.info("Loading extracted CSVs from %s …", d)
    rumen = pd.read_csv(d / "rumen_barn.csv")
    resp_path = d / "respiration_barn.csv"
    resp = pd.read_csv(resp_path) if resp_path.exists() else pd.DataFrame()
    prod = pd.read_csv(d / "production.csv") if (d / "production.csv").exists() else pd.DataFrame()

    # 1. Broken-stick fits
    log.info("═" * 50)
    log.info("1. Broken-stick fits …")
    bs = run_broken_stick_fits(rumen, resp)
    # Merge production data
    if not prod.empty:
        bs = bs.merge(
            prod[["animal_id", "year", "mean_milk_yield_kg", "lactation_nr"]],
            on=["animal_id", "year"], how="left",
        )
    bs.to_csv(d / "broken_stick_results.csv", index=False)
    for prefix, label in [("thi", "THI→body"), ("temp", "Temp→body"),
                          ("resp_thi", "THI→resp"), ("resp_temp", "Temp→resp")]:
        n_conv = (bs[f"{prefix}_converged"] == True).sum()
        log.info("  %s: %d/%d converged", label, n_conv, len(bs))

    # 2. Spearman correlations
    log.info("═" * 50)
    log.info("2. Spearman correlations …")
    spearman = compute_spearman(rumen, resp)
    spearman.to_csv(d / "spearman_correlations.csv", index=False)
    log.info("  Body temp vs THI: median rₛ = %.3f", spearman["thi_rs"].median())
    log.info("  Body temp vs barn temp: median rₛ = %.3f", spearman["temp_rs"].median())

    # 3. Below/above breakpoint
    log.info("═" * 50)
    log.info("3. Below/above breakpoint comparisons …")
    beh = compute_below_above(rumen, resp, bs)
    beh.to_csv(d / "behavioural_response.csv", index=False)
    log.info("  %d animals with paired data", len(beh))

    # 4. Statistical tests with BH-FDR
    log.info("═" * 50)
    log.info("4. Wilcoxon tests (within-year, BH-FDR corrected) …")
    tests = run_statistical_tests(beh)
    tests.to_csv(d / "statistical_tests.csv", index=False)
    for _, t in tests.iterrows():
        log.info(
            "  %d | %-30s | n=%3d | p_raw=%.4f | p_adj=%.4f | %s",
            t["year"], t["test"], t["n"], t["p_raw"], t["p_adj"], t["stars"],
        )

    # 5. Breakpoint stability
    log.info("═" * 50)
    log.info("5. Breakpoint stability (repeat animals) …")
    pairs, icc = compute_stability(bs)
    if not pairs.empty:
        pairs.to_csv(d / "breakpoint_stability.csv", index=False)
    log.info("  ICC = %.3f, %d pairs from %d animals",
             icc, len(pairs), pairs["animal_id"].nunique() if not pairs.empty else 0)

    # 6. Summary table
    log.info("═" * 50)
    log.info("6. Summary table …")
    summary = make_summary_table(bs)
    summary.to_csv(d / "summary_table.csv", index=False)

    log.info("═" * 50)
    log.info("All statistics complete.")


if __name__ == "__main__":
    main()

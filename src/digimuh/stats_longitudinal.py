#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_longitudinal                                   ║
# ║  « breakpoint stability, ICC, and year-over-year tests »       ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Repeat-animal stability analysis: ICC, Fisher resampling of    ║
# ║  year pairs, and manuscript summary table.                      ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Longitudinal breakpoint stability analysis."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger("digimuh.stats")

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

def _run_longitudinal_tests(bs: pd.DataFrame, d: Path) -> None:
    """Fisher resampling tests on longitudinal breakpoints."""
    from digimuh.console import result_table, kv, stars_styled
    from rerandomstats import FisherResamplingTest

    for bp_col, conv_col, label in [
        ("thi_breakpoint", "thi_converged", "THI"),
        ("temp_breakpoint", "temp_converged", "Barn temp"),
    ]:
        conv = bs[bs[conv_col] == True].dropna(subset=[bp_col])
        if conv.empty:
            continue

        # Animals in 2+ years
        year_counts = conv.groupby("animal_id")["year"].nunique()
        repeat_ids = year_counts[year_counts >= 2].index
        if len(repeat_ids) < 5:
            kv(f"{label} repeat animals", f"{len(repeat_ids)} (too few)")
            continue

        repeat = conv[conv["animal_id"].isin(repeat_ids)].copy()
        repeat = repeat.sort_values(["animal_id", "year"])
        first_bp = repeat.groupby("animal_id").first()[bp_col]
        repeat["bp_change"] = repeat[bp_col] - repeat["animal_id"].map(first_bp)
        years = sorted(repeat["year"].unique().astype(int))

        # Pairwise year comparisons (absolute breakpoints)
        test_rows = []
        raw_ps = []
        for i, y1 in enumerate(years):
            for y2 in years[i + 1:]:
                # Paired: same animals in both years
                ids_both = set(repeat[repeat["year"] == y1]["animal_id"]) & \
                           set(repeat[repeat["year"] == y2]["animal_id"])
                if len(ids_both) < 5:
                    continue
                d1 = repeat[(repeat["year"] == y1) & (repeat["animal_id"].isin(ids_both))][bp_col].tolist()
                d2 = repeat[(repeat["year"] == y2) & (repeat["animal_id"].isin(ids_both))][bp_col].tolist()
                p = FisherResamplingTest(
                    data_a=d1, data_b=d2,
                    func="medianDiff", combination_n=20_000,
                ).main()
                test_rows.append([f"{y1} vs {y2}", len(ids_both),
                                  np.median(d1), np.median(d2),
                                  np.median(d2) - np.median(d1), p, ""])
                raw_ps.append(p)

        if test_rows:
            adj_ps = benjamini_hochberg(np.array(raw_ps))
            for r, adj_p in zip(test_rows, adj_ps):
                r[5] = adj_p
                r[6] = stars_styled(p_to_stars(adj_p))
            result_table(
                f"{label} breakpoints: pairwise year comparisons (BH-FDR)",
                ["Comparison", "n paired", "Median Y1", "Median Y2",
                 "Diff", "p adj", "Sig."],
                test_rows,
            )

        # Relative change from first year: test if different from 0
        for year in years[1:]:  # skip first year (change = 0)
            changes = repeat[repeat["year"] == year]["bp_change"].dropna()
            if len(changes) >= 5:
                # Test median change vs zero
                p = FisherResamplingTest(
                    data_a=changes.tolist(),
                    data_b=[0.0] * len(changes),
                    func="medianDiff", combination_n=20_000,
                ).main()
                kv(f"{label} {year} change from baseline",
                   f"median = {changes.median():.1f}, p = {p:.4f} {p_to_stars(p)}")

    # Save longitudinal test results alongside stability
    log.info("  Longitudinal tests complete.")


if __name__ == "__main__":
    main()

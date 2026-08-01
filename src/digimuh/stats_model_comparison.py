#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_model_comparison                                ║
# ║  « threshold vs smooth reaction norm, chosen by AIC »            ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Body temperature is a homeostatically defended variable: it is  ║
# ║  held near a set-point until thermoregulation is overwhelmed,    ║
# ║  then rises.  That is a threshold, not a smooth reaction norm.    ║
# ║  This module makes the claim testable — it compares the          ║
# ║  segmented (broken-stick) model against a linear reaction norm   ║
# ║  and a Hill curve by AIC, per animal-summer and as population    ║
# ║  mixed models — answering the reviewer request for a formal      ║
# ║  model comparison against a mixed random regression model.       ║
# ╚══════════════════════════════════════════════════════════════════╝
"""AIC model comparison: broken-stick threshold vs smooth reaction norm."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.constants import RESAMPLING_SEED

log = logging.getLogger("digimuh.stats")

# Mean-parameter counts (a shared +1 variance parameter cancels in ΔAIC).
_PARAMS = {"linear": 2, "broken-stick": 4, "hill": 4}


# ─────────────────────────────────────────────────────────────
#  « per animal-summer AIC from R² »
#
#  For nested Gaussian models on the same response, the AIC gap
#  is  ΔAIC = n·ln((1−R²_a)/(1−R²_b)) + 2(k_a−k_b), so the whole
#  comparison follows from the three R² already stored per fit
#  plus the reading count — no re-fitting is needed.
# ─────────────────────────────────────────────────────────────

def _aic_score(n: float, r2: float, model: str) -> float:
    """Relative Gaussian AIC (dropping the shared TSS/variance constant)."""
    if not np.isfinite(r2) or r2 >= 1.0 or n < 5:
        return np.nan
    return float(n) * np.log(1.0 - r2) + 2.0 * _PARAMS[model]


def compute_model_comparison(
    bs: pd.DataFrame, predictor: str = "thi",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per animal-summer AIC comparison of the three response models.

    Args:
        bs:        ``broken_stick_results`` frame; must carry
                   ``<predictor>_r_squared`` (broken-stick),
                   ``<predictor>_hill_r2`` (Hill) and
                   ``<predictor>_linear_r2`` (linear reaction norm).
        predictor: ``"thi"`` or ``"temp"``.

    Returns:
        ``(per_animal, summary)``.  ``per_animal`` has one row per
        converged fit with the three AIC scores, the winning model and
        the sub/supra-threshold slopes; ``summary`` has one row per model
        with its lowest-AIC share, how often the broken-stick beats it,
        and the median ΔAIC in favour of the broken-stick.
    """
    conv = bs[(bs[f"{predictor}_converged"] == True)  # noqa: E712
              & bs[f"{predictor}_r_squared"].notna()].copy()
    if conv.empty:
        return pd.DataFrame(), pd.DataFrame()

    n = conv["n_readings"].to_numpy(dtype=float)
    conv["aic_linear"] = [_aic_score(ni, r2, "linear")
                          for ni, r2 in zip(n, conv[f"{predictor}_linear_r2"])]
    conv["aic_broken_stick"] = [_aic_score(ni, r2, "broken-stick")
                                for ni, r2 in zip(n, conv[f"{predictor}_r_squared"])]
    conv["aic_hill"] = [_aic_score(ni, r2, "hill")
                        for ni, r2 in zip(n, conv[f"{predictor}_hill_r2"])]

    aic_cols = {"linear": "aic_linear", "broken-stick": "aic_broken_stick",
                "hill": "aic_hill"}

    def _winner(row) -> str | float:
        scores = {m: row[c] for m, c in aic_cols.items() if np.isfinite(row[c])}
        return min(scores, key=scores.get) if scores else np.nan

    conv["winner"] = conv.apply(_winner, axis=1)
    keep = ["animal_id", "year", "n_readings",
            f"{predictor}_slope_below", f"{predictor}_slope_above",
            f"{predictor}_breakpoint", f"{predictor}_linear_r2",
            f"{predictor}_r_squared", f"{predictor}_hill_r2",
            "aic_linear", "aic_broken_stick", "aic_hill", "winner"]
    per_animal = conv[[c for c in keep if c in conv.columns]].reset_index(drop=True)

    n_fits = len(per_animal)
    win_share = per_animal["winner"].value_counts()
    rows = []
    for model, col in aic_cols.items():
        delta = per_animal[col] - per_animal["aic_broken_stick"]
        delta = delta.replace([np.inf, -np.inf], np.nan).dropna()
        rows.append({
            "predictor": predictor,
            "model": model,
            "params": _PARAMS[model],
            "n_fits": n_fits,
            "lowest_aic_share_pct": round(100.0 * win_share.get(model, 0) / n_fits, 1),
            "pct_beaten_by_broken_stick": (
                round(float((delta > 0).mean()) * 100.0, 1)
                if model != "broken-stick" and len(delta) else 0.0),
            "median_delta_aic_vs_broken_stick": (
                round(float(delta.median()), 1) if len(delta) else 0.0),
        })
    summary = pd.DataFrame(rows)
    return per_animal, summary


# ─────────────────────────────────────────────────────────────
#  « population mixed models — smooth RRM vs threshold hinge »
# ─────────────────────────────────────────────────────────────

def compute_population_rrm(
    rumen: pd.DataFrame, predictor: str = "thi",
    knot_grid: tuple[int, int] = (58, 82), max_per_group: int = 120,
) -> pd.DataFrame:
    """Compare a smooth linear random-regression against a threshold hinge.

    Both are mixed models with an identical random cow intercept, so only
    the fixed functional form differs; the knot of the threshold form is
    chosen by AIC over ``knot_grid``.  Fitted on a balanced subsample of up
    to ``max_per_group`` readings per cow-summer to bound computation.
    """
    import statsmodels.formula.api as smf

    x_col = "barn_thi" if predictor == "thi" else "barn_temp"
    d = rumen.dropna(subset=[x_col, "body_temp"]).copy()
    d["cowyr"] = d["animal_id"].astype(str) + "_" + d["year"].astype(str)
    idx = (d.groupby("cowyr")
           .sample(n=max_per_group, replace=True, random_state=RESAMPLING_SEED)
           .index.unique())
    s = d.loc[idx].reset_index(drop=True)
    x = s[x_col].to_numpy(dtype=float)
    s = s.rename(columns={x_col: "x"})

    m_lin = smf.mixedlm("body_temp ~ x", s, groups=s["animal_id"]).fit(reml=False)
    best = None
    for knot in range(knot_grid[0], knot_grid[1] + 1):
        s["_hinge"] = np.maximum(x - knot, 0.0)
        m = smf.mixedlm("body_temp ~ x + _hinge", s,
                        groups=s["animal_id"]).fit(reml=False)
        if best is None or m.aic < best[1]:
            best = (knot, float(m.aic))
    knot, aic_thr = best
    return pd.DataFrame([
        {"predictor": predictor, "model": "smooth linear random regression",
         "formula": "body_temp ~ THI + (1|cow)", "knot": np.nan,
         "aic": round(float(m_lin.aic), 1),
         "delta_aic_vs_threshold": round(float(m_lin.aic) - aic_thr, 1),
         "n_readings": len(s), "n_cow_summers": int(s["cowyr"].nunique())},
        {"predictor": predictor, "model": "threshold hinge",
         "formula": "body_temp ~ THI + max(THI-knot,0) + (1|cow)", "knot": knot,
         "aic": round(aic_thr, 1), "delta_aic_vs_threshold": 0.0,
         "n_readings": len(s), "n_cow_summers": int(s["cowyr"].nunique())},
    ])


# ─────────────────────────────────────────────────────────────
#  « CLI »
# ─────────────────────────────────────────────────────────────

def main() -> None:
    """Write the AIC model-comparison tables next to the longitudinal figures."""
    from digimuh.paths import resolve_input, resolve_output

    parser = argparse.ArgumentParser(
        description="Threshold vs smooth reaction-norm AIC comparison.")
    parser.add_argument("--data", type=Path, default=Path("results/broken_stick"),
                        help="Analysis output directory (holds 02_breakpoints/).")
    parser.add_argument("--predictor", choices=["thi", "temp"], default="thi")
    parser.add_argument("--population", action="store_true",
                        help="Also fit the population mixed models (slower; "
                             "needs rumen_barn.csv).")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    bs_path = resolve_input(args.data, "broken_stick_results.csv")
    bs = pd.read_csv(bs_path)
    if f"{args.predictor}_linear_r2" not in bs.columns:
        raise SystemExit(
            f"{bs_path} lacks the '{args.predictor}_linear_r2' column — re-run "
            "`digimuh-stats` so the broken-stick fits record the linear R².")

    per_animal, summary = compute_model_comparison(bs, args.predictor)
    per_out = resolve_output(args.data, "breakpoint_model_aic_peranimal.csv")
    sum_out = resolve_output(args.data, "breakpoint_model_aic_comparison.csv")
    per_animal.to_csv(per_out, index=False)
    summary.to_csv(sum_out, index=False)
    log.info("Per-animal AIC   → %s  (%d fits)", per_out, len(per_animal))
    log.info("AIC summary      → %s", sum_out)
    log.info("\n%s", summary.to_string(index=False))

    if args.population:
        rumen = pd.read_csv(resolve_input(args.data, "rumen_barn.csv"))
        pop = compute_population_rrm(rumen, args.predictor)
        pop_out = resolve_output(args.data, "breakpoint_model_aic_population.csv")
        pop.to_csv(pop_out, index=False)
        log.info("Population RRM   → %s", pop_out)
        log.info("\n%s", pop.to_string(index=False))


if __name__ == "__main__":
    main()

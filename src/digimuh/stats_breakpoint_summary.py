#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_breakpoint_summary                              ║
# ║  « descriptive summaries of the breakpoint distribution »        ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Continuity guard, share of readings below each cow's threshold, ║
# ║  and herd percentiles.  Kept out of stats_core so the fitting    ║
# ║  module stays about fitting.                                     ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Descriptive summaries of fitted breakpoints: continuity, imbalance, percentiles."""

from __future__ import annotations

import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────
#  « continuity at the breakpoint »
#
#  The segmented model of Equation 2 is continuous: both segments
#  share the value a + b₁·ψ at the breakpoint.  An earlier fitter
#  fitted two independent lines and jumped there.  This guard runs
#  on every results frame so that can never silently return.
# ─────────────────────────────────────────────────────────────

_FIT_COLUMNS = ("breakpoint", "intercept_below", "intercept_above",
                "slope_below", "slope_above", "converged")


def breakpoint_jump(bs: pd.DataFrame, predictor: str = "thi") -> pd.Series:
    """Discontinuity of each converged fit at its breakpoint (NaN otherwise)."""
    psi = bs[f"{predictor}_breakpoint"]
    left = bs[f"{predictor}_intercept_below"] + bs[f"{predictor}_slope_below"] * psi
    right = bs[f"{predictor}_intercept_above"] + bs[f"{predictor}_slope_above"] * psi
    return (right - left).where(bs[f"{predictor}_converged"] == True)


def check_breakpoint_continuity(bs: pd.DataFrame,
                                predictors: tuple[str, ...] = ("thi", "temp"),
                                tol: float = 1e-6) -> dict[str, float]:
    """Raise if any converged fit is discontinuous; return max |jump| per predictor."""
    worst: dict[str, float] = {}
    for pred in predictors:
        if not {f"{pred}_{c}" for c in _FIT_COLUMNS} <= set(bs.columns):
            continue
        jump = breakpoint_jump(bs, pred).abs()
        worst[pred] = float(jump.max()) if jump.notna().any() else 0.0
        if worst[pred] > tol:
            n_bad = int((jump > tol).sum())
            raise RuntimeError(
                f"{pred}: {n_bad} broken-stick fit(s) are discontinuous at the "
                f"breakpoint (max |jump| = {worst[pred]:.3g}); the fitter is not "
                f"the continuous segmented model of Equation 2")
    return worst


# ─────────────────────────────────────────────────────────────
#  « imbalance: how much of the record sits below the threshold »
# ─────────────────────────────────────────────────────────────

def _reliable_column(bs: pd.DataFrame, predictor: str) -> str:
    flag = f"{predictor}_bp_reliable"
    return flag if flag in bs.columns else f"{predictor}_converged"


def compute_fraction_below_breakpoint(rumen: pd.DataFrame, bs: pd.DataFrame,
                                      predictor: str = "thi") -> pd.DataFrame:
    """Share of paired readings at or below each cow's breakpoint.

    Reliable fits only.  Most readings sit below the threshold, where the
    fitted slope is near zero by construction, so a weak overall
    correlation between climate and rumen temperature is expected even
    when the response above the threshold is strong; this puts a number
    on that imbalance.
    """
    env = "barn_thi" if predictor == "thi" else "barn_temp"
    bp = f"{predictor}_breakpoint"
    fits = bs.loc[bs[_reliable_column(bs, predictor)] == True,
                  ["animal_id", "year", bp]].dropna()
    merged = rumen[["animal_id", "year", env]].merge(fits, on=["animal_id", "year"])
    merged["below"] = merged[env] <= merged[bp]
    out = (merged.groupby(["animal_id", "year"])
                 .agg(breakpoint=(bp, "first"), n_readings=(env, "size"),
                      n_below=("below", "sum"))
                 .reset_index())
    out["fraction_below"] = out["n_below"] / out["n_readings"]
    out.insert(2, "predictor", predictor)
    return out


def summarise_fraction_below(frac: pd.DataFrame) -> dict[str, float]:
    """Median / IQR of the per-cow fraction plus the pooled share of all readings."""
    if frac.empty:
        return {}
    f = frac["fraction_below"]
    return {"n_cow_summers": int(len(frac)),
            "median": float(f.median()),
            "q25": float(f.quantile(0.25)), "q75": float(f.quantile(0.75)),
            "pooled": float(frac["n_below"].sum() / frac["n_readings"].sum())}


# ─────────────────────────────────────────────────────────────
#  « herd percentiles: the median is a description, not a trigger »
# ─────────────────────────────────────────────────────────────

def compute_breakpoint_percentiles(bs: pd.DataFrame,
                                   predictors: tuple[str, ...] = ("thi", "temp"),
                                   percentiles: tuple[int, ...] = (10, 25, 50, 75, 90),
                                   ) -> pd.DataFrame:
    """Herd breakpoint percentiles over reliable fits, pooled and per summer.

    The herd median marks conditions at which half the herd is above its
    own threshold; a management trigger would sit at a lower percentile.
    """
    rows = []
    for pred in predictors:
        col = f"{pred}_breakpoint"
        if col not in bs.columns:
            continue
        conv = bs.loc[bs[_reliable_column(bs, pred)] == True, ["year", col]].dropna()
        if conv.empty:
            continue
        groups = [("all", conv[col])] + [(int(y), g[col])
                                          for y, g in conv.groupby("year")]
        for year, vals in groups:
            row: dict = {"predictor": pred, "year": year, "n": int(len(vals))}
            row.update({f"p{p}": float(np.percentile(vals, p)) for p in percentiles})
            rows.append(row)
    return pd.DataFrame(rows)

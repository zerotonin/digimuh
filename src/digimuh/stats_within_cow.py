#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_within_cow                                      ║
# ║  « pairwise summer contrasts that stay inside the cow »          ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Stratified (Freedman-Lane) within-cow permutation: for two       ║
# ║  summers the cow effect is removed by centring, and permuting    ║
# ║  summer labels within cow is a paired sign-flip on the           ║
# ║  within-cow differences.  Repeated-measures valid where the      ║
# ║  pooled Fisher resampling of cow-summers was not.                ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Within-cow (paired sign-flip) permutation contrasts between summers."""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from rerandomstats import correct_pvalues_array

from digimuh.constants import RESAMPLING_SEED
from digimuh.stats_core import p_to_stars

N_PERMUTATIONS = 20_000
MIN_PAIRED = 3


def within_cow_signflip(values_a: np.ndarray, values_b: np.ndarray, *,
                        n_perm: int = N_PERMUTATIONS,
                        seed: int = RESAMPLING_SEED) -> tuple[float, float]:
    """Paired sign-flip permutation on within-cow differences ``a − b``.

    Returns the observed mean difference and the two-sided p-value
    ``(k + 1) / (n_perm + 1)``, where ``k`` counts permuted mean
    differences at least as extreme as the observed one.
    """
    d = np.asarray(values_a, dtype=float) - np.asarray(values_b, dtype=float)
    d = d[np.isfinite(d)]
    if d.size < MIN_PAIRED:
        return np.nan, np.nan
    observed = float(d.mean())
    rng = np.random.default_rng(seed)
    flips = rng.choice(np.array([-1.0, 1.0]), size=(n_perm, d.size))
    permuted = (flips * d).mean(axis=1)
    k = int(np.sum(np.abs(permuted) >= abs(observed)))
    return observed, (k + 1) / (n_perm + 1)


def compute_within_cow_posthoc(long: pd.DataFrame, *, metric: str,
                               predictor: str,
                               n_perm: int = N_PERMUTATIONS,
                               seed: int = RESAMPLING_SEED) -> pd.DataFrame:
    """All pairwise summer contrasts of one per-cow-summer outcome.

    Args:
        long:      One row per cow-summer with ``animal_id``, ``year``,
                   ``value``.
        metric:    Outcome label written into the table.
        predictor: ``"thi"`` or ``"temp"``.

    Returns:
        One row per year pair: group sizes, ``n_paired`` (cows measured in
        both summers — the only ones that identify the contrast), the mean
        within-cow difference ``A − B``, the raw p-value, and the BH-FDR
        corrected p-value over the pairs of this metric × predictor.
    """
    long = long.dropna(subset=["value"])
    years = sorted(int(y) for y in long["year"].unique())
    wide = long.pivot_table(index="animal_id", columns="year",
                            values="value", aggfunc="first")
    rows = []
    for a, b in combinations(years, 2):
        both = wide[[a, b]].dropna() if a in wide and b in wide else pd.DataFrame()
        diff, p = (within_cow_signflip(both[a].to_numpy(), both[b].to_numpy(),
                                       n_perm=n_perm, seed=seed)
                   if len(both) >= MIN_PAIRED else (np.nan, np.nan))
        rows.append({"metric": metric, "predictor": predictor,
                     "groupA": a, "groupB": b,
                     "nA": int((long["year"] == a).sum()),
                     "nB": int((long["year"] == b).sum()),
                     "n_paired": int(len(both)),
                     "mean_within_cow_diff": diff, "p_value": p})
    out = pd.DataFrame(rows)
    out["p_value_fdr"] = np.nan
    valid = out["p_value"].notna()
    if valid.any():
        out.loc[valid, "p_value_fdr"] = correct_pvalues_array(
            out.loc[valid, "p_value"].to_numpy(), method="fdr_bh")
    out["sig_level"] = [p_to_stars(p) if np.isfinite(p) else "n.a."
                        for p in out["p_value_fdr"]]
    return out


# ─────────────────────────────────────────────────────────────
#  « the three raincloud outcomes »
# ─────────────────────────────────────────────────────────────

def breakpoint_value_long(bs: pd.DataFrame, predictor: str) -> pd.DataFrame:
    """Reliable breakpoint per cow-summer as a long frame."""
    keep = (f"{predictor}_bp_reliable" if f"{predictor}_bp_reliable" in bs.columns
            else f"{predictor}_converged")
    conv = bs[bs[keep] == True]
    return (conv[["animal_id", "year", f"{predictor}_breakpoint"]]
            .rename(columns={f"{predictor}_breakpoint": "value"}))


def crossing_count_long(crossing_times: pd.DataFrame, predictor: str) -> pd.DataFrame:
    """Upward crossings per cow-summer (cow-summers with no crossing are absent)."""
    sub = crossing_times[crossing_times["predictor"] == predictor]
    return (sub.groupby(["animal_id", "year"]).size()
            .rename("value").reset_index())


def minutes_over_long(minutes_over: pd.DataFrame, predictor: str) -> pd.DataFrame:
    """Minutes above the breakpoint per cow-summer."""
    sub = minutes_over[minutes_over["predictor"] == predictor]
    return sub[["animal_id", "year", "minutes_over"]].rename(
        columns={"minutes_over": "value"})


def compute_all_within_cow_posthoc(bs: pd.DataFrame, crossing_times: pd.DataFrame,
                                   minutes_over: pd.DataFrame) -> pd.DataFrame:
    """Within-cow pairwise contrasts for all three outcomes × both predictors."""
    frames = []
    for predictor in ("thi", "temp"):
        for metric, long in (
            ("breakpoint_value", breakpoint_value_long(bs, predictor)),
            ("crossing_count", crossing_count_long(crossing_times, predictor)),
            ("minutes_over", minutes_over_long(minutes_over, predictor)),
        ):
            if long.empty or long["year"].nunique() < 2:
                continue
            frames.append(compute_within_cow_posthoc(
                long, metric=metric, predictor=predictor))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

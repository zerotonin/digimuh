#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_retest                                          ║
# ║  « does last summer's threshold predict this summer's? »         ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Test-retest statistics on breakpoints of consecutive summers:   ║
# ║  correlation, regression slope against 0 and against 1, the      ║
# ║  within-cow shift with its limits of agreement, and how much of  ║
# ║  that shift belongs to the summer rather than to the cow.        ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Test-retest statistics for breakpoints of consecutive summers."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from digimuh.stats_breakpoint_summary import reliable_fits

MIN_PAIRS = 5
"""Fewest consecutive-summer pairs for which the statistics are computed."""

LOA_Z = 1.96
"""Multiplier of the SD of the shift for the 95 % limits of agreement."""


# ─────────────────────────────────────────────────────────────
#  « pairs of consecutive summers »
# ─────────────────────────────────────────────────────────────

def build_consecutive_pairs(bs: pd.DataFrame,
                            predictor: str = "thi") -> pd.DataFrame:
    """One row per cow and pair of consecutive summers.

    Only identified fits enter, and only summers that follow each other
    directly: a cow measured in 2021 and 2023 but not 2022 gives no pair.
    A cow measured in three consecutive summers gives two pairs.

    Returns:
        ``animal_id``, ``from_year``, ``to_year``, ``bp_last``,
        ``bp_this``, and ``se_last`` / ``se_this`` when the results carry
        the breakpoint standard error.
    """
    bp, se = f"{predictor}_breakpoint", f"{predictor}_breakpoint_se"
    has_se = se in bs.columns
    fits = reliable_fits(bs, predictor).dropna(subset=[bp])
    fits = fits[["animal_id", "year", bp] + ([se] if has_se else [])].copy()
    fits["year"] = fits["year"].astype(int)

    earlier = fits.rename(columns={"year": "from_year", bp: "bp_last",
                                   se: "se_last"})
    earlier["to_year"] = earlier["from_year"] + 1
    later = fits.rename(columns={"year": "to_year", bp: "bp_this",
                                 se: "se_this"})
    pairs = earlier.merge(later, on=["animal_id", "to_year"])
    columns = ["animal_id", "from_year", "to_year", "bp_last", "bp_this"]
    if has_se:
        columns += ["se_last", "se_this"]
    return (pairs[columns].sort_values(["animal_id", "from_year"])
            .reset_index(drop=True))


# ─────────────────────────────────────────────────────────────
#  « the statistics »
#
#  The slope is read against two references.  Slope 1 is a
#  breakpoint that belongs to the cow; slope 0 is a breakpoint
#  that last summer says nothing about.  The test of r = 0 and
#  the test of slope = 0 are the same test.
# ─────────────────────────────────────────────────────────────

def compute_retest_statistics(pairs: pd.DataFrame) -> dict | None:
    """Agreement between the breakpoints of consecutive summers.

    Pairs of the same cow are not independent of each other, so the
    p-values are somewhat optimistic; ``n_animals`` is reported beside
    ``n_pairs`` for that reason.

    Args:
        pairs: Output of :func:`build_consecutive_pairs`.

    Returns:
        ``None`` for fewer than :data:`MIN_PAIRS` pairs or a constant
        variable.  Otherwise a dict with

        * ``pearson_r`` and its 95 % CI and p, ``spearman_rs`` and p;
        * ``slope`` of this summer on last summer with its 95 % CI,
          ``slope_p_vs_0`` and ``slope_p_vs_1``, ``intercept``,
          ``r_squared``;
        * ``shift_mean`` (later minus earlier) with its 95 % CI,
          ``shift_sd`` and ``loa``, the half-width of the 95 % limits of
          agreement (Bland & Altman 1986);
        * ``anova_*`` and ``between_year_share``: one-way analysis of
          variance of the shift with the pair of summers as the factor,
          and the share of the shift's variance that lies between pairs
          of summers (herd level) rather than within them (individual);
        * ``reliability`` and ``slope_disattenuated``: the share of the
          variance of last summer's breakpoints that is not fit noise,
          and the slope corrected for it.  The standard errors ignore the
          serial correlation of the readings, so the correction is a
          lower bound.
    """
    if len(pairs) < MIN_PAIRS:
        return None
    x = pairs["bp_last"].to_numpy(dtype=float)
    y = pairs["bp_this"].to_numpy(dtype=float)
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    n = len(pairs)

    pearson = stats.pearsonr(x, y)
    pearson_ci = pearson.confidence_interval(0.95)
    spearman = stats.spearmanr(x, y)
    fit = stats.linregress(x, y)
    t_slope = stats.t.ppf(0.975, n - 2)
    p_vs_one = 2.0 * stats.t.sf(abs((fit.slope - 1.0) / fit.stderr), n - 2)

    shift = y - x
    shift_sd = float(shift.std(ddof=1))
    half_width = stats.t.ppf(0.975, n - 1) * shift_sd / np.sqrt(n)

    result = {
        "n_pairs": int(n),
        "n_animals": int(pairs["animal_id"].nunique()),
        "pearson_r": float(pearson.statistic),
        "pearson_ci_lo": float(pearson_ci.low),
        "pearson_ci_hi": float(pearson_ci.high),
        "pearson_p": float(pearson.pvalue),
        "spearman_rs": float(spearman.statistic),
        "spearman_p": float(spearman.pvalue),
        "slope": float(fit.slope),
        "slope_ci_lo": float(fit.slope - t_slope * fit.stderr),
        "slope_ci_hi": float(fit.slope + t_slope * fit.stderr),
        "slope_p_vs_0": float(fit.pvalue),
        "slope_p_vs_1": float(p_vs_one),
        "intercept": float(fit.intercept),
        "r_squared": float(fit.rvalue ** 2),
        "shift_mean": float(shift.mean()),
        "shift_ci_lo": float(shift.mean() - half_width),
        "shift_ci_hi": float(shift.mean() + half_width),
        "shift_sd": shift_sd,
        "loa": LOA_Z * shift_sd,
    }
    result.update(_partition_shift(shift, pairs["to_year"].to_numpy()))
    result.update(_disattenuate(pairs, x, fit.slope))
    return result


def _partition_shift(shift: np.ndarray, to_year: np.ndarray) -> dict:
    """One-way ANOVA of the shift by pair of summers."""
    empty = {"anova_f": np.nan, "anova_df1": np.nan, "anova_df2": np.nan,
             "anova_p": np.nan, "between_year_share": np.nan}
    groups = [shift[to_year == year] for year in np.unique(to_year)]
    groups = [g for g in groups if g.size >= 2]
    if len(groups) < 2:
        return empty
    pooled = np.concatenate(groups)
    ss_total = float(np.sum((pooled - pooled.mean()) ** 2))
    if ss_total == 0:
        return empty
    ss_between = float(sum(g.size * (g.mean() - pooled.mean()) ** 2
                           for g in groups))
    anova = stats.f_oneway(*groups)
    return {"anova_f": float(anova.statistic),
            "anova_df1": len(groups) - 1,
            "anova_df2": int(pooled.size - len(groups)),
            "anova_p": float(anova.pvalue),
            "between_year_share": ss_between / ss_total}


def _disattenuate(pairs: pd.DataFrame, x: np.ndarray, slope: float) -> dict:
    """Slope corrected for the fit noise in last summer's breakpoint."""
    empty = {"reliability": np.nan, "slope_disattenuated": np.nan}
    if "se_last" not in pairs.columns:
        return empty
    se = pairs["se_last"].to_numpy(dtype=float)
    se = se[np.isfinite(se)]
    if se.size == 0:
        return empty
    reliability = 1.0 - float(np.mean(se ** 2)) / float(x.var(ddof=1))
    if reliability <= 0:
        return {"reliability": reliability, "slope_disattenuated": np.nan}
    return {"reliability": reliability,
            "slope_disattenuated": float(slope) / reliability}


def compute_retest_summary(bs: pd.DataFrame,
                           predictors: tuple[str, ...] = ("thi", "temp"),
                           ) -> pd.DataFrame:
    """Test-retest statistics per predictor, one row each."""
    rows = []
    for predictor in predictors:
        if f"{predictor}_breakpoint" not in bs.columns:
            continue
        result = compute_retest_statistics(build_consecutive_pairs(bs, predictor))
        if result is not None:
            rows.append({"predictor": predictor, **result})
    return pd.DataFrame(rows)

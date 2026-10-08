"""Tests for the test-retest statistics on consecutive-summer breakpoints."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from digimuh.stats_retest import (
    build_consecutive_pairs,
    compute_retest_statistics,
    compute_retest_summary,
)


def _results(rows: list[tuple]) -> pd.DataFrame:
    """Breakpoint results from (animal_id, year, breakpoint, reliable)."""
    return pd.DataFrame([
        {"animal_id": a, "year": y, "thi_breakpoint": bp,
         "thi_breakpoint_se": 0.3, "thi_converged": True,
         "thi_bp_reliable": ok}
        for a, y, bp, ok in rows])


def _pairs(slope: float, n: int = 200, noise: float = 1.0, shift: float = 0.0,
           seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    last = rng.normal(76.0, 3.0, n)
    this = 76.0 + slope * (last - 76.0) + shift + rng.normal(0.0, noise, n)
    return pd.DataFrame({"animal_id": np.arange(n), "from_year": 2021,
                         "to_year": 2022, "bp_last": last, "bp_this": this,
                         "se_last": 0.3, "se_this": 0.3})


# ─────────────────────────────────────────────────────────────────
#  pairs
# ─────────────────────────────────────────────────────────────────


def test_pairs_need_summers_that_follow_each_other() -> None:
    bs = _results([
        (1, 2021, 75.0, True), (1, 2022, 76.0, True), (1, 2023, 77.0, True),
        (2, 2021, 74.0, True), (2, 2023, 78.0, True),      # gap: no pair
        (3, 2022, 73.0, True), (3, 2023, 72.0, False),     # unidentified
        (4, 2024, 80.0, True),                             # single summer
    ])
    pairs = build_consecutive_pairs(bs, "thi")
    assert list(pairs["animal_id"]) == [1, 1]
    assert list(pairs["from_year"]) == [2021, 2022]
    assert list(pairs["bp_last"]) == [75.0, 76.0]
    assert list(pairs["bp_this"]) == [76.0, 77.0]
    assert {"se_last", "se_this"} <= set(pairs.columns)


# ─────────────────────────────────────────────────────────────────
#  statistics
# ─────────────────────────────────────────────────────────────────


def test_weak_carry_over_is_told_apart_from_both_references() -> None:
    res = compute_retest_statistics(_pairs(slope=0.2, noise=2.0))
    assert res["slope_ci_lo"] < 0.2 < res["slope_ci_hi"]
    assert res["slope_p_vs_0"] < 0.05          # not flat
    assert res["slope_p_vs_1"] < 1e-10         # nowhere near stable
    assert res["n_pairs"] == 200 and res["n_animals"] == 200


def test_stable_trait_has_slope_one() -> None:
    res = compute_retest_statistics(_pairs(slope=1.0, noise=0.3))
    assert res["slope"] == pytest.approx(1.0, abs=0.05)
    assert res["slope_p_vs_1"] > 0.05
    assert res["pearson_r"] > 0.95


def test_correlation_and_slope_share_one_test() -> None:
    res = compute_retest_statistics(_pairs(slope=0.3, noise=2.0, seed=2))
    assert res["pearson_p"] == pytest.approx(res["slope_p_vs_0"], rel=1e-9)
    assert res["r_squared"] == pytest.approx(res["pearson_r"] ** 2)
    assert res["pearson_ci_lo"] < res["pearson_r"] < res["pearson_ci_hi"]
    assert -1.0 <= res["spearman_rs"] <= 1.0


def test_shift_and_limits_of_agreement() -> None:
    pairs = _pairs(slope=1.0, noise=1.0, shift=1.5, seed=3)
    res = compute_retest_statistics(pairs)
    shift = pairs["bp_this"] - pairs["bp_last"]
    assert res["shift_mean"] == pytest.approx(shift.mean())
    assert res["shift_ci_lo"] < 1.5 < res["shift_ci_hi"]
    assert res["loa"] == pytest.approx(1.96 * shift.std(ddof=1))


def test_shift_variance_is_attributed_to_the_summer() -> None:
    """Every cow moves with her herd: the shift is a herd-level matter."""
    rng = np.random.default_rng(4)
    frames = []
    for to_year, herd_shift in ((2022, -2.0), (2023, 0.0), (2024, 3.0)):
        last = rng.normal(76.0, 3.0, 60)
        frames.append(pd.DataFrame({
            "animal_id": [f"{to_year}-{i}" for i in range(60)],
            "from_year": to_year - 1, "to_year": to_year, "bp_last": last,
            "bp_this": last + herd_shift + rng.normal(0.0, 0.3, 60)}))
    res = compute_retest_statistics(pd.concat(frames, ignore_index=True))
    assert res["between_year_share"] > 0.9
    assert res["anova_p"] < 1e-6
    assert (res["anova_df1"], res["anova_df2"]) == (2, 177)


def test_one_pair_of_summers_has_no_partition() -> None:
    res = compute_retest_statistics(_pairs(slope=0.5))
    assert np.isnan(res["between_year_share"]) and np.isnan(res["anova_p"])


def test_fit_noise_correction() -> None:
    pairs = _pairs(slope=0.5, noise=1.0, seed=5)
    res = compute_retest_statistics(pairs)
    expected = 1.0 - 0.3 ** 2 / pairs["bp_last"].var(ddof=1)
    assert res["reliability"] == pytest.approx(expected)
    assert res["slope_disattenuated"] == pytest.approx(res["slope"] / expected)
    without_se = compute_retest_statistics(pairs.drop(columns=["se_last", "se_this"]))
    assert np.isnan(without_se["reliability"])


def test_too_few_pairs_or_constant_values_give_nothing() -> None:
    assert compute_retest_statistics(_pairs(slope=0.5).head(4)) is None
    flat = _pairs(slope=0.5).assign(bp_last=75.0)
    assert compute_retest_statistics(flat) is None


def test_summary_has_one_row_per_predictor() -> None:
    rng = np.random.default_rng(6)
    rows = []
    for cow in range(40):
        level = rng.normal(76.0, 3.0)
        for year in (2021, 2022):
            rows.append((cow, year, level + rng.normal(0.0, 1.0), True))
    summary = compute_retest_summary(_results(rows))
    assert list(summary["predictor"]) == ["thi"]      # no temp columns given
    assert summary.loc[0, "n_pairs"] == 40

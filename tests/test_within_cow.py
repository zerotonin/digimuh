"""Tests for the within-cow (paired sign-flip) permutation contrasts."""

from __future__ import annotations

import numpy as np
import pandas as pd

from digimuh.stats_within_cow import (
    compute_all_within_cow_posthoc,
    compute_within_cow_posthoc,
    within_cow_signflip,
)


def test_signflip_detects_a_consistent_within_cow_shift() -> None:
    rng = np.random.default_rng(1)
    a = rng.normal(75.0, 3.0, 40)
    b = a + 1.0 + rng.normal(0.0, 0.3, 40)        # every cow ~1 unit higher in b
    diff, p = within_cow_signflip(a, b, n_perm=5_000)
    assert abs(diff + 1.0) < 0.2
    assert p < 0.001


def test_signflip_is_null_under_symmetric_noise() -> None:
    rng = np.random.default_rng(2)
    a = rng.normal(75.0, 3.0, 40)
    b = a + rng.normal(0.0, 1.0, 40)
    _, p = within_cow_signflip(a, b, n_perm=5_000)
    assert p > 0.05


def test_signflip_is_seeded_and_bounded() -> None:
    a = np.arange(10.0)
    b = a + np.array([0.5, -0.4, 0.6, -0.2, 0.3, -0.5, 0.1, 0.2, -0.3, 0.4])
    p1 = within_cow_signflip(a, b, n_perm=2_000)[1]
    p2 = within_cow_signflip(a, b, n_perm=2_000)[1]
    assert p1 == p2
    assert 1.0 / 2_001 <= p1 <= 1.0
    assert np.isnan(within_cow_signflip(a[:2], b[:2])[1])   # too few pairs


def _long() -> pd.DataFrame:
    rng = np.random.default_rng(3)
    rows = []
    for cow in range(30):
        base = rng.normal(75.0, 2.0)
        rows.append({"animal_id": cow, "year": 2021, "value": base})
        rows.append({"animal_id": cow, "year": 2022, "value": base + 2.0})
        if cow < 10:                                  # only ten cows reach 2023
            rows.append({"animal_id": cow, "year": 2023, "value": base + 2.0})
    return pd.DataFrame(rows)


def test_pairwise_table_counts_pairs_and_corrects() -> None:
    out = compute_within_cow_posthoc(_long(), metric="x", predictor="thi",
                                     n_perm=2_000)
    assert list(out["groupA"]) == [2021, 2021, 2022]
    assert list(out["n_paired"]) == [30, 10, 10]
    assert list(out["nA"]) == [30, 30, 30] and list(out["nB"]) == [30, 10, 10]
    r = out.set_index(["groupA", "groupB"])
    assert r.loc[(2021, 2022), "mean_within_cow_diff"] == -2.0
    assert r.loc[(2021, 2022), "p_value_fdr"] < 0.01
    assert r.loc[(2022, 2023), "p_value"] == 1.0          # identical values
    assert (out["p_value_fdr"] >= out["p_value"]).all()
    assert set(out["sig_level"]) <= {"n.s.", "*", "**", "***"}


def test_all_outcomes_cover_metrics_and_predictors() -> None:
    bs = pd.DataFrame({
        "animal_id": [1, 1, 2, 2, 3, 3], "year": [2021, 2022] * 3,
        "thi_breakpoint": [75, 76, 74, 77, 73, 75], "thi_converged": True,
        "thi_bp_reliable": True,
        "temp_breakpoint": [28, 29, 27, 30, 26, 28], "temp_converged": True,
        "temp_bp_reliable": True,
    })
    crossings = pd.DataFrame({
        "animal_id": [1, 1, 1, 2, 2, 3, 3, 3, 3],
        "year": [2021, 2021, 2022, 2021, 2022, 2021, 2022, 2022, 2022],
        "predictor": "thi",
    })
    minutes = pd.DataFrame({
        "animal_id": [1, 1, 2, 2, 3, 3], "year": [2021, 2022] * 3,
        "predictor": "thi", "minutes_over": [100, 200, 150, 250, 120, 300],
    })
    out = compute_all_within_cow_posthoc(bs, crossings, minutes)
    got = set(zip(out["metric"], out["predictor"]))
    assert got == {("breakpoint_value", "thi"), ("crossing_count", "thi"),
                   ("minutes_over", "thi"), ("breakpoint_value", "temp")}
    cc = out[(out["metric"] == "crossing_count")].iloc[0]
    assert cc["n_paired"] == 3
    # per-cow 2021 − 2022 crossing differences: 2−1, 1−1, 1−3 → mean −1/3
    assert abs(cc["mean_within_cow_diff"] + 1.0 / 3.0) < 1e-12

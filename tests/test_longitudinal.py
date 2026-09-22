"""Tests for stats_longitudinal: weighted across-summer test, convergence
audit, and the breakpoint-CI overlap statistic.  Synthetic data only."""

from __future__ import annotations

import numpy as np
import pandas as pd

from digimuh.stats_longitudinal import (
    across_summer_test,
    compute_ci_overlap,
    compute_convergence_audit,
    compute_convergence_rates,
    format_across_summer,
)


def _panel(seed: int = 0) -> pd.DataFrame:
    """30 cows × 3 summers with a drifting outcome and a per-row SE."""
    rng = np.random.default_rng(seed)
    rows = []
    for cow in range(30):
        base = rng.normal(75.0, 2.0)
        for i, year in enumerate((2021, 2022, 2023)):
            rows.append({"animal_id": cow, "year": year,
                         "y": base + 0.8 * i + rng.normal(0.0, 1.0),
                         "se": rng.uniform(0.1, 1.0)})
    return pd.DataFrame(rows)


def _bs() -> pd.DataFrame:
    """12 cows × 2 summers; every third cow never converges."""
    rows = []
    for cow in range(12):
        conv = cow % 3 != 0
        for year in (2021, 2022):
            bp = 75.0 + (cow % 4)
            rows.append({
                "animal_id": cow, "year": year,
                "n_readings": 500 if conv else 80,
                "lactation_nr": 2, "mean_milk_yield_kg": 30.0,
                "thi_converged": conv, "thi_bp_reliable": conv,
                "thi_breakpoint": bp if conv else np.nan,
                "thi_breakpoint_ci_lo": bp - 0.5 if conv else np.nan,
                "thi_breakpoint_ci_hi": bp + 0.5 if conv else np.nan,
            })
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────
#  across_summer_test with weights
# ─────────────────────────────────────────────────────────────────


def test_across_summer_weighted_adds_gee_line() -> None:
    df = _panel()
    res = across_summer_test(df["y"], df["year"], df["animal_id"],
                             kind="continuous", weights=1.0 / df["se"] ** 2)
    assert res["w_model"].startswith("Weighted GEE")
    assert 0.0 <= res["w_p"] <= 1.0
    assert res["n_weighted"] == len(df)
    assert "Weighted GEE" in format_across_summer(res)


def test_across_summer_weights_drop_nonpositive_rows_only() -> None:
    df = _panel()
    w = 1.0 / df["se"] ** 2
    w.iloc[0] = np.nan
    w.iloc[1] = 0.0
    res = across_summer_test(df["y"], df["year"], df["animal_id"], weights=w)
    assert res["n_obs"] == len(df)          # unweighted fit keeps every row
    assert res["n_weighted"] == len(df) - 2


def test_across_summer_without_weights_has_no_w_keys() -> None:
    df = _panel()
    res = across_summer_test(df["y"], df["year"], df["animal_id"])
    assert "w_model" not in res


# ─────────────────────────────────────────────────────────────────
#  convergence rates + audit
# ─────────────────────────────────────────────────────────────────


def test_convergence_rates_per_year_and_pooled() -> None:
    rates = compute_convergence_rates(_bs())
    thi = rates[rates["predictor"] == "thi"]
    assert set(thi["year"]) == {2021, 2022, "all"}
    assert thi.loc[thi["year"] == "all", "n_converged"].item() == 16
    assert thi.loc[thi["year"] == 2021, "rate"].item() == 8 / 12


def test_convergence_audit_detects_short_records() -> None:
    bs = _bs()
    ts = pd.date_range("2021-06-01", periods=48, freq="h")
    rumen = pd.concat([
        pd.DataFrame({"animal_id": c, "year": y, "timestamp": ts,
                      "barn_thi": 60.0 + (c % 5)})
        for c in range(12) for y in (2021, 2022)
    ], ignore_index=True)
    audit = compute_convergence_audit(bs, rumen, "thi")
    assert {"n_days", "max_predictor"} <= set(audit["covariate"])
    row = audit[audit["covariate"] == "n_readings"].iloc[0]
    assert row["median_converged"] > row["median_not_converged"]
    assert row["p_adj"] < 0.05
    assert row["n_converged"] == 16 and row["n_not_converged"] == 8


# ─────────────────────────────────────────────────────────────────
#  CI overlap
# ─────────────────────────────────────────────────────────────────


def test_ci_overlap_flags_disjoint_consecutive_pair() -> None:
    bs = _bs()
    moved = (bs["animal_id"] == 1) & (bs["year"] == 2022)
    bs.loc[moved, ["thi_breakpoint", "thi_breakpoint_ci_lo",
                   "thi_breakpoint_ci_hi"]] = [79.0, 78.5, 79.5]
    pairs, proportion = compute_ci_overlap(bs, "thi")
    assert len(pairs) == 8                    # 8 converged cows, one pair each
    by_cow = pairs.set_index("animal_id")
    assert not by_cow.loc[1, "overlap"]
    assert by_cow.loc[1, "gap"] > 0
    assert by_cow.loc[2, "overlap"]
    assert abs(proportion - 7 / 8) < 1e-12


def test_ci_overlap_empty_when_no_consecutive_years() -> None:
    bs = _bs()
    bs.loc[bs["year"] == 2022, "year"] = 2024
    pairs, proportion = compute_ci_overlap(bs, "thi")
    assert pairs.empty and np.isnan(proportion)

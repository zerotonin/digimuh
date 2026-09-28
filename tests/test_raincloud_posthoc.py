"""Tests for the herd-level post-hoc behind the raincloud letters."""

from __future__ import annotations

import numpy as np
import pandas as pd

from digimuh.viz_longitudinal import _compact_letters, _posthoc_group_medians


def _three_summers(shifted_year: int = 2022, n_repeat: int = 40,
                   n_single: int = 20) -> pd.DataFrame:
    """A herd over three summers: repeat cows plus single-summer cows.

    Repeat cows carry their own level into every summer, so their values
    are correlated between years; one summer sits higher for everyone.
    """
    rng = np.random.default_rng(4)
    level = rng.normal(8.4, 0.5, n_repeat)
    rows = []
    for year in (2021, 2022, 2023):
        shift = 0.8 if year == shifted_year else 0.0
        for cow in range(n_repeat):
            rows.append({"animal_id": f"repeat-{cow}", "year": year,
                         "value": float(np.exp(level[cow] + shift
                                               + rng.normal(0, 0.2)))})
        for cow in range(n_single):
            rows.append({"animal_id": f"{year}-{cow}", "year": year,
                         "value": float(rng.lognormal(8.4 + shift, 0.5))})
    return pd.DataFrame(rows)


def test_herd_posthoc_separates_the_shifted_summer() -> None:
    long = _three_summers()
    res = _posthoc_group_medians(long, [2021, 2022, 2023], combination_n=2_000)
    assert {"groupA", "groupB", "n_paired", "n_only_A", "n_only_B",
            "p value", "p value corrected", "h",
            "median_A", "median_B", "median_diff"} <= set(res.columns)
    assert len(res) == 3
    assert (res["n_paired"] == 40).all()
    assert (res["n_only_A"] == 20).all() and (res["n_only_B"] == 20).all()
    assert (res["p value corrected"] >= res["p value"] - 1e-12).all()
    letters = _compact_letters([2021, 2022, 2023], res)
    # letters run alphabetically from the first year down
    assert letters == {"2021": "a", "2022": "b", "2023": "a"}


def test_herd_posthoc_reports_the_median_difference() -> None:
    long = _three_summers()
    res = _posthoc_group_medians(long, [2021, 2022, 2023], combination_n=500)
    row = res[(res["groupA"] == "2021") & (res["groupB"] == "2022")].iloc[0]
    med = long.groupby("year")["value"].median()
    assert row["median_diff"] == med[2021] - med[2022]
    assert row["median_A"] == med[2021] and row["median_B"] == med[2022]
    assert row["median_diff"] < 0


def test_herd_posthoc_needs_two_populated_summers() -> None:
    long = _three_summers()
    only_one = long[long["year"] == 2021]
    assert _posthoc_group_medians(only_one, [2021], combination_n=500) is None
    sparse = pd.concat([only_one, long[long["year"] == 2022].head(2)])
    assert _posthoc_group_medians(sparse, [2021, 2022], combination_n=500) is None

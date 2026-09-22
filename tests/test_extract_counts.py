"""Tests for the extract-stage filter accounting (no database needed)."""

from __future__ import annotations

import pandas as pd

from digimuh.extract import FILTER_STAGES, summarise_filter_counts


def test_summarise_filter_counts_totals_and_shares() -> None:
    counts = pd.DataFrame([
        {"animal_id": 1, "year": 2021, "n_raw": 1000, "n_in_range": 990,
         "n_after_drink": 900, "n_matched": 880, "n_kept": 600},
        {"animal_id": 2, "year": 2021, "n_raw": 1000, "n_in_range": 1000,
         "n_after_drink": 950, "n_matched": 900, "n_kept": 580},
    ])
    out = summarise_filter_counts(counts)
    assert list(out["stage"]) == list(FILTER_STAGES)
    assert list(out["n"]) == [2000, 1990, 1850, 1780, 1180]
    assert out["pct_of_raw"].iloc[0] == 100.0
    assert out["pct_of_raw"].iloc[-1] == 59.0
    assert list(out["n_removed"]) == [0, 10, 140, 70, 600]

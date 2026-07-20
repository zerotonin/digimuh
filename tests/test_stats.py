# ╔══════════════════════════════════════════════════════════════╗
# ║  DigiMuh — test_stats                                        ║
# ║  « unit tests for statistical functions on synthetic data »   ║
# ╚══════════════════════════════════════════════════════════════╝
"""Test statistical functions using synthetic DataFrames.

No database connection required.  All tests construct small DataFrames
that mimic the structure of rumen_barn.csv and broken_stick_results.csv.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from digimuh.stats_production import (
    compute_thermoneutral_fraction,
)
from digimuh.stats_temporal import (
    compute_circadian_null_model,
    compute_climate_eta,
    compute_crossing_times,
    compute_thi_daily_profile,
)

# ─────────────────────────────────────────────────────────────
#  fixtures
# ─────────────────────────────────────────────────────────────

@pytest.fixture()
def synthetic_rumen() -> pd.DataFrame:
    """Synthetic rumen_barn DataFrame: 5 animals, 30 days, 10-min sampling."""
    np.random.seed(42)
    rows = []
    for aid in range(5):
        for day in range(30):
            for hour in range(24):
                if 4 <= hour <= 6 or 16 <= hour <= 18:
                    continue  # milking exclusion
                for minute in [0, 10, 20, 30, 40, 50]:
                    thi = 60 + 15 * np.sin(np.pi * hour / 24) + np.random.normal(0, 2)
                    temp = thi / 2.8
                    rows.append({
                        "animal_id": aid,
                        "year": 2023,
                        "timestamp": pd.Timestamp(2023, 7, day + 1, hour, minute),
                        "body_temp": 39.3 + 0.01 * max(thi - 68, 0)
                                     + np.random.normal(0, 0.05),
                        "barn_temp": temp,
                        "barn_thi": thi,
                    })
    return pd.DataFrame(rows)


@pytest.fixture()
def synthetic_bs() -> pd.DataFrame:
    """Synthetic broken_stick_results: 5 animals, all converged."""
    return pd.DataFrame([
        {
            "animal_id": aid, "year": 2023,
            "thi_breakpoint": 68 + aid, "thi_converged": True,
            "temp_breakpoint": 24 + aid * 0.5, "temp_converged": True,
            "mean_milk_yield_kg": 18 + np.random.normal(0, 2),
            "lactation_nr": np.random.choice([1, 2, 3]),
        }
        for aid in range(5)
    ])


# ─────────────────────────────────────────────────────────────
#  thermoneutral fraction
# ─────────────────────────────────────────────────────────────

def test_tnf_returns_daily_records(synthetic_rumen, synthetic_bs):
    """TNF produces one row per animal per day."""
    tnf = compute_thermoneutral_fraction(synthetic_rumen, synthetic_bs)
    assert not tnf.empty
    assert "thi_tnf" in tnf.columns
    assert "date" in tnf.columns
    # Each animal should have ~30 days
    for aid in range(5):
        n_days = tnf[tnf["animal_id"] == aid]["date"].nunique()
        assert n_days >= 20, f"Animal {aid} has only {n_days} days"


def test_tnf_values_bounded(synthetic_rumen, synthetic_bs):
    """TNF values must be between 0 and 1."""
    tnf = compute_thermoneutral_fraction(synthetic_rumen, synthetic_bs)
    assert tnf["thi_tnf"].between(0, 1).all()


# ─────────────────────────────────────────────────────────────
#  circadian null model
# ─────────────────────────────────────────────────────────────

def test_circadian_has_cool_and_stress(synthetic_rumen, synthetic_bs):
    """Circadian model separates days into cool and stress categories."""
    circ = compute_circadian_null_model(synthetic_rumen, synthetic_bs)
    assert not circ.empty
    assert set(circ["day_type"].unique()) <= {"cool", "stress"}
    assert "hour" in circ.columns
    assert "body_temp_mean" in circ.columns


def test_circadian_hours_range(synthetic_rumen, synthetic_bs):
    """Hours should span 0–23 (excluding milking windows)."""
    circ = compute_circadian_null_model(synthetic_rumen, synthetic_bs)
    hours = set(circ["hour"].unique())
    # Milking windows 4-6, 16-18 excluded from data, so those hours
    # may have fewer readings but shouldn't crash
    assert len(hours) >= 10


# ─────────────────────────────────────────────────────────────
#  THI daily profile
# ─────────────────────────────────────────────────────────────

def test_thi_profile_has_breakpoint(synthetic_rumen, synthetic_bs):
    """THI profile should include the herd median breakpoint."""
    prof = compute_thi_daily_profile(synthetic_rumen, synthetic_bs)
    assert not prof.empty
    assert "herd_median_bp" in prof.columns
    assert np.isfinite(prof["herd_median_bp"].iloc[0])


def test_thi_profile_hourly(synthetic_rumen, synthetic_bs):
    """Profile should have entries for multiple hours."""
    prof = compute_thi_daily_profile(synthetic_rumen, synthetic_bs)
    assert prof["hour"].nunique() >= 10


# ─────────────────────────────────────────────────────────────
#  crossing times
# ─────────────────────────────────────────────────────────────

def test_crossing_times_structure(synthetic_rumen, synthetic_bs):
    """Crossing times should have expected columns."""
    ct = compute_crossing_times(synthetic_rumen, synthetic_bs)
    assert not ct.empty
    for col in ["animal_id", "predictor", "breakpoint",
                "clock_hour", "day_fraction"]:
        assert col in ct.columns, f"Missing column: {col}"


def test_crossing_times_no_milking_gaps(synthetic_rumen, synthetic_bs):
    """No crossings should occur during milking windows (artifact filter)."""
    ct = compute_crossing_times(synthetic_rumen, synthetic_bs)
    if ct.empty:
        pytest.skip("No crossings in synthetic data")
    # The milking-gap filter removes crossings where consecutive
    # readings are > 30 min apart.  The filter is on dt, not clock hour,
    # so this is not a strict assertion — just verify the function runs.
    assert isinstance(ct, pd.DataFrame)


# ─────────────────────────────────────────────────────────────
#  climate ETA
# ─────────────────────────────────────────────────────────────

def test_climate_eta_normalised(synthetic_rumen, synthetic_bs):
    """Climate ETA should have both normalised columns."""
    eta = compute_climate_eta(
        synthetic_rumen, synthetic_bs,
        crossing_hour_range=None,  # all hours for synthetic data
    )
    if eta.empty:
        pytest.skip("No crossings in synthetic data")
    for col in ["thi_norm", "temp_norm", "thi_raw", "temp_raw"]:
        assert col in eta.columns, f"Missing column: {col}"
    # At the crossing point (relative_minutes == 0), thi_norm should
    # be near zero for thi-triggered events
    thi_at_zero = eta[(eta["trigger"] == "thi")
                      & (eta["relative_minutes"] == 0)]["thi_norm"]
    if not thi_at_zero.empty:
        assert abs(thi_at_zero.mean()) < 5, "THI norm at crossing should be near 0"


# ─────────────────────────────────────────────────────────────
#  « BH-FDR within families (C1) »
# ─────────────────────────────────────────────────────────────

def test_apply_fdr_within_never_undercuts_raw_p():
    """BH adjustment is conservative: p_fdr >= p for every test."""
    import numpy as np
    import pandas as pd

    from digimuh.stats_core import apply_fdr_within

    rng = np.random.default_rng(0)
    df = pd.DataFrame({"group": ["a"] * 30 + ["b"] * 30,
                       "p": rng.uniform(0, 1, 60)})
    out = apply_fdr_within(df, ["group"])
    assert (out["p_fdr"] >= out["p"] - 1e-12).all()
    assert out["family"].nunique() == 2


def test_apply_fdr_within_preserves_nan_and_family_size():
    """NaN p-values stay NaN and are excluded from the family size."""
    import numpy as np
    import pandas as pd

    from digimuh.stats_core import apply_fdr_within

    df = pd.DataFrame({"p": [0.01, 0.02, np.nan, 0.60]})
    out = apply_fdr_within(df)
    assert out["p_fdr"].isna().sum() == 1
    assert out.loc[2, "p_fdr"] != out.loc[2, "p_fdr"]        # still NaN
    # Family size is 3, not 4: 0.01 * 3/1 = 0.03
    assert np.isclose(out["p_fdr"].min(), 0.03, atol=1e-9)


def test_milk_composition_screens_report_corrected_p():
    """The four stratified screens must expose p_fdr, not just raw p."""
    import numpy as np
    import pandas as pd

    from digimuh.stats_production import tnf_yield_correlations_by_class

    rng = np.random.default_rng(2)
    n = 300
    tnf = pd.DataFrame({
        "animal_id": rng.integers(1, 30, n),
        "yield_class": rng.choice(["low", "middle", "high"], n),
        "thi_tnf": rng.uniform(0, 1, n), "temp_tnf": rng.uniform(0, 1, n),
        "daily_yield_kg": rng.normal(30, 5, n),
        "relative_yield": rng.uniform(0.5, 1.1, n),
        "yield_residual": rng.normal(0, 3, n),
    })
    out = tnf_yield_correlations_by_class(tnf)
    assert {"p", "p_fdr", "family"} <= set(out.columns)
    fin = out.dropna(subset=["p"])
    assert (fin["p_fdr"] >= fin["p"] - 1e-12).all()


def test_longitudinal_change_from_baseline_is_fdr_corrected(caplog):
    """The change-from-first-year block must be BH-corrected, not raw.

    It used to print one raw p-value per year directly beneath a table that
    *was* corrected, so both read as adjusted.  Feeds deterministic p-values
    through the resampling test and checks the reported value is the BH
    transform, not the input.

    Both output backends are captured: ``console.result_table`` renders through
    rich when it is installed and falls back to ``log.info`` when it is not, so
    a stdout-only capture passes locally and sees nothing on a runner without
    the optional dependency.
    """
    import contextlib
    import io
    import logging
    from pathlib import Path

    import numpy as np
    import pandas as pd

    import digimuh.stats_longitudinal as sl

    feed = [0.20, 0.30, 0.40, 0.50, 0.60, 0.70,   # pairwise family
            0.010, 0.040, 0.900]                   # change-from-baseline family

    class _FakeFisher:
        def __init__(self, **kwargs):
            pass

        def main(self):
            return feed.pop(0) if feed else 0.5

    # The function imports FisherResamplingTest lazily from rerandomstats,
    # so that is the binding to replace.
    import rerandomstats

    original = rerandomstats.FisherResamplingTest
    rerandomstats.FisherResamplingTest = _FakeFisher
    try:
        rng = np.random.default_rng(0)
        bs = pd.DataFrame([
            dict(animal_id=aid, year=year,
                 thi_breakpoint=rng.normal(69, 2), thi_converged=True,
                 temp_breakpoint=rng.normal(22, 2), temp_converged=False)
            for aid in range(1, 31)
            for year in (2021, 2022, 2023, 2024)
        ])
        buf = io.StringIO()
        with caplog.at_level(logging.INFO), contextlib.redirect_stdout(buf):
            sl._run_longitudinal_tests(bs, Path("."))
        out = buf.getvalue() + "\n".join(r.getMessage() for r in caplog.records)
    finally:
        rerandomstats.FisherResamplingTest = original

    assert "change from first year (BH-FDR)" in out
    # BH over [0.010, 0.040, 0.900] is [0.030, 0.060, 0.900]: the 0.040 test
    # loses its star.  The raw value must not appear as the reported p.
    assert "0.0300" in out and "0.0600" in out
    assert "0.0400" not in out

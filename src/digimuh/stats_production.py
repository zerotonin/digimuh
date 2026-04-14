#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_production                                     ║
# ║  « thermoneutral fraction and milk yield coupling »             ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Compute daily thermoneutral fraction (TNF) and correlate with  ║
# ║  P95-normalised daily milk yield.                               ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Production impact analysis for the broken-stick pipeline."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

log = logging.getLogger("digimuh.stats")

def compute_thermoneutral_fraction(
    rumen: pd.DataFrame, bs_results: pd.DataFrame,
) -> pd.DataFrame:
    """Compute daily thermoneutral fraction (TNF) per animal.

    For each animal with a converged THI breakpoint, count the fraction
    of 10-min readings per calendar day where barn THI is at or below
    the individual breakpoint.  This is the fraction of the day the
    animal spends within its thermoneutral zone.

    Also computes analogous fraction for barn temperature breakpoints.

    Args:
        rumen: rumen_barn.csv DataFrame (needs timestamp, barn_thi, barn_temp).
        bs_results: broken_stick_results.csv DataFrame.

    Returns:
        DataFrame with columns: animal_id, year, date, thi_tnf, temp_tnf,
        n_readings, mean_thi, mean_body_temp.
    """
    records = []
    converged = bs_results[bs_results["thi_converged"] == True]

    for _, row in converged.iterrows():
        aid = int(row["animal_id"])
        year = int(row["year"])
        thi_bp = row["thi_breakpoint"]

        # Barn temp breakpoint (may not converge)
        temp_bp = row.get("temp_breakpoint", np.nan)
        temp_conv = row.get("temp_converged", False)

        grp = rumen[(rumen["animal_id"] == aid) & (rumen["year"] == year)].copy()
        if len(grp) < 50:
            continue

        grp["date"] = pd.to_datetime(grp["timestamp"]).dt.date

        for date, day_data in grp.groupby("date"):
            n = len(day_data)
            if n < 6:  # at least 1 hour of data
                continue

            thi_tnf = (day_data["barn_thi"] <= thi_bp).mean()
            temp_tnf = (day_data["barn_temp"] <= temp_bp).mean() if temp_conv else np.nan

            records.append({
                "animal_id": aid,
                "year": year,
                "date": date,
                "thi_tnf": thi_tnf,
                "temp_tnf": temp_tnf,
                "n_readings": n,
                "mean_thi": day_data["barn_thi"].mean(),
                "mean_body_temp": day_data["body_temp"].mean(),
            })

    return pd.DataFrame(records)


def compute_tnf_yield_analysis(
    tnf: pd.DataFrame, daily_yield: pd.DataFrame,
) -> pd.DataFrame:
    """Merge daily TNF with daily milk yield for per-cow, per-day pairs.

    For each cow, each day produces two floats:
    - thi_tnf: fraction of the day below the cow's THI breakpoint (0-1)
    - relative_yield: that day's milk yield / cow-specific P95 (0-1ish)

    The P95 is the 95th percentile of each individual cow's daily yields
    across her entire dataset (all years).  This is a robust, cow-specific
    reference maximum that avoids outlier sensitivity of the absolute max.

    Args:
        tnf: Daily TNF DataFrame from compute_thermoneutral_fraction.
        daily_yield: Daily milk yield DataFrame (animal_id, date,
            daily_yield_kg, year).

    Returns:
        DataFrame with one row per cow-day: animal_id, year, date,
        thi_tnf, temp_tnf, daily_yield_kg, yield_p95, relative_yield.
    """
    if tnf.empty or daily_yield.empty:
        return pd.DataFrame()

    # Ensure date columns are comparable
    tnf = tnf.copy()
    daily_yield = daily_yield.copy()
    tnf["date"] = pd.to_datetime(tnf["date"]).dt.date
    daily_yield["date"] = pd.to_datetime(daily_yield["date"]).dt.date

    # Merge TNF and yield on (animal_id, date)
    merged = tnf.merge(
        daily_yield[["animal_id", "date", "daily_yield_kg", "n_milkings"]],
        on=["animal_id", "date"],
        how="inner",
    )

    if merged.empty:
        return pd.DataFrame()

    # Cow-specific P95: robust reference maximum from ALL that cow's daily yields
    cow_p95 = daily_yield.groupby("animal_id")["daily_yield_kg"].agg(
        yield_p95=lambda x: np.nanpercentile(x.dropna(), 95),
        yield_median=lambda x: x.median(),
        n_yield_days=lambda x: x.notna().sum(),
    ).reset_index()

    merged = merged.merge(cow_p95, on="animal_id", how="left")

    # Relative yield: this day's yield / this cow's P95
    merged["relative_yield"] = merged["daily_yield_kg"] / merged["yield_p95"]

    return merged


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────

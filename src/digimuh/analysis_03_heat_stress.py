#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  Analysis 3 — Heat stress multi-sensor fusion                ║
# ║  « rumen temp × THI × water intake × respiration »           ║
# ╚══════════════════════════════════════════════════════════════╝
"""Per-animal heat stress response modelling.

Uses the ``v_analysis_heat_stress`` view which combines daily
smaXtec rumen data, DWD weather, gouna respiration, and milking
production.

The analysis:

1.  Builds **per-animal Z-scored rumen temperature** following
    the approach of the NZ smaXtec study (JDS Communications,
    2024): each cow's temperature distribution is scaled to a
    common mean and SD before thresholding.
2.  Fits per-animal **thermoregulatory dose-response curves**:
    rumen_temp_z = f(THI) using sigmoidal regression.  The
    inflection point and slope characterise each animal's
    heat tolerance.
3.  Computes a daily **heat load index** fusing rumen temp,
    respiration rate, activity suppression, and water intake.
4.  Quantifies the **production impact**: milk yield loss per
    unit of heat load.

Usage::

    python -m digimuh.analysis_03_heat_stress --db cow.db --out results/heat

References:
    Identifying and predicting heat stress events for grazing
    dairy cows using rumen temperature boluses. JDS Comm. 2024.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from digimuh.analysis_utils import connect_db, query_df, setup_plotting, save_fig

log = logging.getLogger("digimuh.heat_stress")


# ─────────────────────────────────────────────────────────────
#  « data loading »
# ─────────────────────────────────────────────────────────────

SQL_HEAT = """
SELECT *
FROM v_analysis_heat_stress
WHERE rumen_temp_mean IS NOT NULL
  AND dwd_thi_max IS NOT NULL
  AND smaxtec_readings >= 20
ORDER BY animal_id, day
"""


def load_heat_data(con) -> pd.DataFrame:
    """Load heat stress view and add per-animal Z-scored temperature.

    Args:
        con: Database connection with views.

    Returns:
        DataFrame with per-animal Z-scored rumen temperature and
        a composite heat load index.
    """
    log.info("Loading heat stress data …")
    df = query_df(con, SQL_HEAT)
    log.info("  %d animal-days with rumen temp + THI", len(df))

    if df.empty:
        return df

    # ── per-animal Z-score (the NZ approach) ─────────────────
    stats = df.groupby("animal_id")["rumen_temp_clean_mean"].agg(
        ["mean", "std"],
    )
    stats.columns = ["cow_temp_mean", "cow_temp_std"]
    df = df.merge(stats, on="animal_id", how="left")
    df["rumen_temp_z"] = (
        (df["rumen_temp_clean_mean"] - df["cow_temp_mean"])
        / df["cow_temp_std"].replace(0, np.nan)
    )

    # ── heat stress indicator (Z > 1.5 = stressed) ──────────
    df["heat_stressed"] = (df["rumen_temp_z"] > 1.5).astype(int)

    # ── composite heat load index (multi-sensor) ─────────────
    # Each component Z-scored, positive = more stress
    components = {}

    for col, invert in [
        ("rumen_temp_z", False),      # high temp = stress
        ("resp_mean", False),         # high resp = stress
        ("act_index_mean", True),     # low activity = stress
        ("water_liter", False),       # high water = compensation
        ("rum_index_mean", True),     # low rumination = stress
    ]:
        if col in df.columns:
            vals = df[col].dropna()
            if len(vals) > 10:
                mu, sd = vals.mean(), vals.std()
                if sd > 0:
                    z = (df[col] - mu) / sd
                    if invert:
                        z = -z
                    components[col] = z

    if components:
        comp_df = pd.DataFrame(components)
        df["heat_load_index"] = comp_df.mean(axis=1)
    else:
        df["heat_load_index"] = np.nan

    return df


# ─────────────────────────────────────────────────────────────
#  « thermoregulatory dose-response »
# ─────────────────────────────────────────────────────────────

def sigmoid(x: np.ndarray, L: float, k: float, x0: float, b: float) -> np.ndarray:
    """Four-parameter sigmoid: L / (1 + exp(-k*(x - x0))) + b."""
    return L / (1.0 + np.exp(-k * (x - x0))) + b


def fit_dose_response(df: pd.DataFrame, min_days: int = 30) -> pd.DataFrame:
    """Fit per-animal sigmoid dose-response: rumen_temp_z = f(THI).

    The inflection point (x0) represents the THI at which the
    cow's temperature begins to rise sharply — its personal heat
    tolerance threshold.

    Args:
        df: DataFrame from :func:`load_heat_data`.
        min_days: Minimum number of observation days per animal.

    Returns:
        DataFrame with one row per animal, columns: ``animal_id``,
        ``thi_threshold`` (x0), ``slope`` (k), ``n_days``.
    """
    results = []

    for aid, grp in df.groupby("animal_id"):
        sub = grp.dropna(subset=["dwd_thi_max", "rumen_temp_z"])
        if len(sub) < min_days:
            continue

        x = sub["dwd_thi_max"].values
        y = sub["rumen_temp_z"].values

        try:
            popt, _ = curve_fit(
                sigmoid, x, y,
                p0=[2.0, 0.1, 68.0, -1.0],
                maxfev=5000,
                bounds=([0, 0.001, 40, -5], [10, 2, 90, 5]),
            )
            results.append({
                "animal_id": aid,
                "sigmoid_L": popt[0],
                "slope_k": popt[1],
                "thi_threshold": popt[2],
                "baseline_b": popt[3],
                "n_days": len(sub),
            })
        except (RuntimeError, ValueError):
            continue

    result_df = pd.DataFrame(results)
    log.info(
        "Fitted dose-response for %d / %d animals",
        len(result_df), df["animal_id"].nunique(),
    )
    return result_df


# ─────────────────────────────────────────────────────────────
#  « production impact »
# ─────────────────────────────────────────────────────────────

def compute_production_impact(df: pd.DataFrame) -> pd.DataFrame:
    """Estimate milk yield loss attributable to heat load.

    Bins the heat_load_index into quartiles and computes mean
    milk yield per bin.

    Args:
        df: DataFrame from :func:`load_heat_data`.

    Returns:
        Summary DataFrame with heat load bins and mean production.
    """
    sub = df.dropna(subset=["heat_load_index", "milk_yield_kg"])
    if sub.empty:
        return pd.DataFrame()

    sub = sub.copy()
    sub["heat_quartile"] = pd.qcut(
        sub["heat_load_index"], q=4, labels=["Q1 (cool)", "Q2", "Q3", "Q4 (hot)"],
    )
    summary = sub.groupby("heat_quartile", observed=True).agg(
        milk_mean=("milk_yield_kg", "mean"),
        milk_std=("milk_yield_kg", "std"),
        temp_z_mean=("rumen_temp_z", "mean"),
        resp_mean=("resp_mean", "mean"),
        water_mean=("water_liter", "mean"),
        n_days=("milk_yield_kg", "count"),
    ).reset_index()

    return summary


# ─────────────────────────────────────────────────────────────
#  « plots »
# ─────────────────────────────────────────────────────────────

def plot_heat_overview(df: pd.DataFrame, dose: pd.DataFrame,
                       out_dir: Path) -> None:
    """Generate heat stress analysis plots.

    Args:
        df: Full analysis DataFrame.
        dose: Per-animal dose-response fit results.
        out_dir: Output directory.
    """
    import matplotlib.pyplot as plt
    setup_plotting()

    # ── THI vs rumen temperature Z-score (scatter + sigmoid) ─
    sub = df.dropna(subset=["dwd_thi_max", "rumen_temp_z"])
    if len(sub) > 100:
        fig, ax = plt.subplots(figsize=(9, 6))
        ax.scatter(
            sub["dwd_thi_max"], sub["rumen_temp_z"],
            s=3, alpha=0.1, c="#888780",
        )
        ax.axhline(1.5, color="#A32D2D", linestyle="--",
                    label="Heat stress threshold (Z=1.5)")
        ax.axvline(68, color="#BA7517", linestyle="--",
                    label="THI=68 (mild stress onset)")
        ax.set_xlabel("DWD THI max")
        ax.set_ylabel("Per-animal Z-scored rumen temperature")
        ax.legend()
        ax.set_title("Rumen temperature response to ambient heat load")
        save_fig(fig, "heat_thi_vs_rumen_temp", out_dir)

    # ── THI threshold distribution ───────────────────────────
    if not dose.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(dose["thi_threshold"], bins=25, color="#1D9E75",
                edgecolor="white")
        ax.axvline(68, color="#A32D2D", linestyle="--",
                    label="Classical THI=68 threshold")
        ax.set_xlabel("Individual THI threshold (sigmoid inflection)")
        ax.set_ylabel("Number of animals")
        ax.legend()
        ax.set_title("Per-animal heat tolerance thresholds")
        save_fig(fig, "heat_individual_thresholds", out_dir)

    # ── production impact ────────────────────────────────────
    impact = compute_production_impact(df)
    if not impact.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.bar(
            impact["heat_quartile"], impact["milk_mean"],
            yerr=impact["milk_std"] / np.sqrt(impact["n_days"]),
            color=["#85B7EB", "#FAC775", "#F0997B", "#E24B4A"],
            edgecolor="white",
        )
        ax.set_ylabel("Mean daily milk yield (kg)")
        ax.set_xlabel("Heat load index quartile")
        ax.set_title("Milk production by heat load quartile")
        save_fig(fig, "heat_production_impact", out_dir)


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────

def main() -> None:
    """Entry point for heat stress analysis."""
    parser = argparse.ArgumentParser(description="Heat stress analysis")
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("results/heat"))
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(name)-20s │ %(message)s",
        datefmt="%H:%M:%S",
    )

    con = connect_db(args.db)
    df = load_heat_data(con)

    if df.empty:
        log.warning("No data — check ingestion.")
        return

    args.out.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out / "heat_stress_data.csv", index=False)

    dose = fit_dose_response(df)
    dose.to_csv(args.out / "dose_response_fits.csv", index=False)

    plot_heat_overview(df, dose, args.out)
    con.close()


if __name__ == "__main__":
    main()

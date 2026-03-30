#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  Analysis 6 — Digestive efficiency composite                 ║
# ║  « motility × pH × milk composition (time-lagged) »         ║
# ╚══════════════════════════════════════════════════════════════╝
"""Rumen mechanical–chemical–production coupling analysis.

Uses ``v_analysis_digestive`` which joins daily smaXtec motility/pH
with HerdePlus MLP milk composition.

The key insight: reticulorumen contraction patterns (motility) drive
mixing, mixing drives fermentation rate, fermentation determines the
volatile fatty acid profile, and VFA ratios directly shape milk fat
and protein.  This pipeline has a multi-day lag.

The analysis:

1.  Computes **time-lagged cross-correlations** between daily
    motility metrics and the next available MLP test-day values.
2.  Builds a **digestive efficiency score** from the motility–pH
    coupling: animals where motility and pH co-vary tightly
    have well-functioning rumens.
3.  Tests whether digestive efficiency predicts milk composition
    at the next MLP test day.

Usage::

    python -m digimuh.analysis_06_digestive --db cow.db --out results/digestive
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.analysis_utils import connect_db, query_df, setup_plotting, save_fig

log = logging.getLogger("digimuh.digestive")


# ─────────────────────────────────────────────────────────────
#  « data loading »
# ─────────────────────────────────────────────────────────────

SQL_DIGESTIVE = """
SELECT *
FROM v_analysis_digestive
WHERE smaxtec_readings >= 20
ORDER BY animal_id, day
"""


def load_digestive_data(con) -> pd.DataFrame:
    """Load digestive analysis view.

    Args:
        con: Database connection with views.

    Returns:
        DataFrame with daily motility/pH and sparse MLP composition.
    """
    log.info("Loading digestive efficiency data …")
    df = query_df(con, SQL_DIGESTIVE)
    log.info("  %d animal-days", len(df))
    return df


# ─────────────────────────────────────────────────────────────
#  « time-lagged cross-correlation »
# ─────────────────────────────────────────────────────────────

def compute_lagged_correlations(
    df: pd.DataFrame,
    predictor_cols: list[str],
    target_cols: list[str],
    max_lag_days: int = 14,
) -> pd.DataFrame:
    """Compute cross-correlations between daily rumen metrics and
    MLP test-day values at various time lags.

    For each animal, MLP test days are identified (non-null target),
    and the mean of each predictor over the preceding N days is
    correlated with the target value.

    Args:
        df: Full digestive DataFrame, sorted by animal_id + day.
        predictor_cols: Daily rumen metrics to use as predictors.
        target_cols: MLP test-day columns (sparse).
        max_lag_days: Maximum look-back window in days.

    Returns:
        DataFrame: lag × predictor × target → correlation.
    """
    results = []

    # Identify MLP test days
    available_targets = [c for c in target_cols if c in df.columns]
    available_preds = [c for c in predictor_cols if c in df.columns]

    if not available_targets or not available_preds:
        log.warning("Insufficient columns for lagged correlation")
        return pd.DataFrame()

    df = df.copy()
    df["day"] = pd.to_datetime(df["day"])

    for target in available_targets:
        # Rows where MLP value exists
        mlp_rows = df[df[target].notna()].copy()
        if len(mlp_rows) < 20:
            continue

        for lag in range(1, max_lag_days + 1):
            for pred in available_preds:
                # For each MLP test day, get mean of predictor
                # over [day-lag, day-1]
                vals = []
                for _, row in mlp_rows.iterrows():
                    window = df[
                        (df["animal_id"] == row["animal_id"])
                        & (df["day"] >= row["day"] - pd.Timedelta(days=lag))
                        & (df["day"] < row["day"])
                    ]
                    pred_mean = window[pred].mean()
                    if pd.notna(pred_mean) and pd.notna(row[target]):
                        vals.append((pred_mean, row[target]))

                if len(vals) >= 15:
                    arr = np.array(vals)
                    r = np.corrcoef(arr[:, 0], arr[:, 1])[0, 1]
                    results.append({
                        "predictor": pred,
                        "target": target,
                        "lag_days": lag,
                        "correlation": r,
                        "n_pairs": len(vals),
                    })

    return pd.DataFrame(results)


# ─────────────────────────────────────────────────────────────
#  « digestive efficiency score »
# ─────────────────────────────────────────────────────────────

def compute_digestive_efficiency(df: pd.DataFrame, window: int = 7) -> pd.DataFrame:
    """Compute a per-animal rolling digestive efficiency score.

    Efficiency is defined as the strength of coupling between
    motility (contraction interval) and rumen pH over a rolling
    window.  In a well-functioning rumen, shorter contraction
    intervals (faster mixing) correspond to lower pH (more
    active fermentation) — a negative correlation.

    Args:
        df: Full digestive DataFrame.
        window: Rolling window size in days.

    Returns:
        DataFrame with ``animal_id``, ``day``, ``digest_eff``,
        and ``digest_eff_rank`` (percentile within herd-day).
    """
    required = ["mot_period_mean", "ph_mean"]
    if not all(c in df.columns for c in required):
        log.warning("Missing motility/pH columns for efficiency score")
        return pd.DataFrame()

    df = df.copy()
    df["day"] = pd.to_datetime(df["day"])
    df = df.sort_values(["animal_id", "day"])

    records = []
    for aid, grp in df.groupby("animal_id"):
        sub = grp.dropna(subset=required).copy()
        if len(sub) < window:
            continue
        # Rolling correlation between mot_period and pH
        sub["digest_eff"] = (
            sub["mot_period_mean"]
            .rolling(window, min_periods=window // 2)
            .corr(sub["ph_mean"])
        )
        records.append(sub[["animal_id", "day", "digest_eff"]])

    if not records:
        return pd.DataFrame()

    result = pd.concat(records, ignore_index=True)

    # Herd-level percentile rank per day
    result["digest_eff_rank"] = result.groupby("day")["digest_eff"].rank(
        pct=True,
    )

    log.info(
        "Digestive efficiency computed for %d animals, %d days",
        result["animal_id"].nunique(),
        result["day"].nunique(),
    )
    return result


# ─────────────────────────────────────────────────────────────
#  « plots »
# ─────────────────────────────────────────────────────────────

def plot_digestive_results(
    lagged: pd.DataFrame, efficiency: pd.DataFrame, out_dir: Path,
) -> None:
    """Generate digestive analysis plots.

    Args:
        lagged: Lagged correlation results.
        efficiency: Digestive efficiency scores.
        out_dir: Output directory.
    """
    import matplotlib.pyplot as plt
    setup_plotting()

    # ── lag-correlation heatmap ───────────────────────────────
    if not lagged.empty:
        for target in lagged["target"].unique():
            sub = lagged[lagged["target"] == target]
            pivot = sub.pivot_table(
                index="predictor", columns="lag_days",
                values="correlation",
            )
            if pivot.empty:
                continue

            fig, ax = plt.subplots(figsize=(10, 5))
            im = ax.imshow(
                pivot.values, aspect="auto", cmap="RdBu_r",
                vmin=-0.5, vmax=0.5,
            )
            ax.set_xticks(range(pivot.shape[1]))
            ax.set_xticklabels(pivot.columns)
            ax.set_yticks(range(pivot.shape[0]))
            ax.set_yticklabels(pivot.index, fontsize=9)
            ax.set_xlabel("Lag (days before MLP test)")
            ax.set_title(f"Lagged correlations → {target}")
            fig.colorbar(im, ax=ax, label="Pearson r")
            save_fig(fig, f"digestive_lag_{target}", out_dir)

    # ── efficiency score distribution ────────────────────────
    if not efficiency.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(
            efficiency["digest_eff"].dropna(), bins=50,
            color="#0F6E56", edgecolor="white",
        )
        ax.axvline(0, color="#A32D2D", linestyle="--",
                    label="Zero coupling")
        ax.set_xlabel("Digestive efficiency (mot-pH correlation)")
        ax.set_ylabel("Count (animal-days)")
        ax.legend()
        ax.set_title("Distribution of motility–pH coupling strength")
        save_fig(fig, "digestive_efficiency_dist", out_dir)


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────

def main() -> None:
    """Entry point for digestive efficiency analysis."""
    parser = argparse.ArgumentParser(description="Digestive efficiency analysis")
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("results/digestive"))
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(name)-20s │ %(message)s",
        datefmt="%H:%M:%S",
    )

    con = connect_db(args.db)
    df = load_digestive_data(con)

    if df.empty:
        log.warning("No data.")
        return

    args.out.mkdir(parents=True, exist_ok=True)

    # Lagged cross-correlation
    lagged = compute_lagged_correlations(
        df,
        predictor_cols=[
            "mot_period_mean", "mot_pw_mean", "ph_mean",
            "ph_under_58_total", "rum_index_mean", "act_index_mean",
        ],
        target_cols=["mlp_fat_pct", "mlp_protein_pct", "mlp_fpr", "mlp_ecm"],
    )
    if not lagged.empty:
        lagged.to_csv(args.out / "lagged_correlations.csv", index=False)

    # Digestive efficiency
    efficiency = compute_digestive_efficiency(df)
    if not efficiency.empty:
        efficiency.to_csv(args.out / "digestive_efficiency.csv", index=False)

    plot_digestive_results(lagged, efficiency, args.out)
    con.close()


if __name__ == "__main__":
    main()

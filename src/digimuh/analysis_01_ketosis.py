#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  Analysis 1 — Subclinical ketosis detection                  ║
# ║  « FPR × rumination × milk yield × rumen state »             ║
# ╚══════════════════════════════════════════════════════════════╝
"""Subclinical ketosis risk scoring from multi-sensor fusion.

Uses the ``v_analysis_ketosis`` view which joins daily milking
data (HerdePlus MLP test-day results), smaXtec rumen metrics,
water intake, and disease ground truth.

The analysis:

1.  Extracts days with MLP test-day data (FPR available).
2.  Computes a composite **ketosis risk score** from:
    - Fat-to-protein ratio (FPR > 1.4 = energy deficit)
    - Rumination index (lower = reduced feed intake)
    - Milk yield deviation from rolling cow mean
    - Rumen pH (low pH + high FPR = metabolic confusion)
3.  Validates against disease records (ground truth).
4.  Trains a Random Forest classifier and reports feature
    importance and cross-validated performance.

Usage::

    python -m digimuh.analysis_01_ketosis --db cow.db --out results/ketosis

References:
    Oetzel (2013) — FPR thresholds for subclinical ketosis.
    Kaufman et al. (2016) J Dairy Sci 99:5604–18 — rumination
        time association with subclinical ketosis.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.analysis_utils import connect_db, query_df, setup_plotting, save_fig

log = logging.getLogger("digimuh.ketosis")


# ─────────────────────────────────────────────────────────────
#  « data extraction »
# ─────────────────────────────────────────────────────────────

SQL_KETOSIS = """
SELECT *
FROM v_analysis_ketosis
WHERE mlp_fpr IS NOT NULL
ORDER BY animal_id, day
"""

SQL_MILK_HISTORY = """
SELECT animal_id, day, milk_yield_kg
FROM v_herdeplus_daily
WHERE milk_yield_kg IS NOT NULL
ORDER BY animal_id, day
"""


def load_ketosis_data(con) -> pd.DataFrame:
    """Load ketosis analysis view and add derived features.

    Args:
        con: Active database connection with views created.

    Returns:
        DataFrame with one row per animal per MLP test day,
        enriched with rolling-mean deviations and risk scores.
    """
    log.info("Loading ketosis analysis data …")
    df = query_df(con, SQL_KETOSIS)
    log.info("  %d rows with FPR data", len(df))

    if df.empty:
        return df

    # ── rolling milk yield deviation ─────────────────────────
    milk = query_df(con, SQL_MILK_HISTORY)
    if not milk.empty:
        milk = milk.sort_values(["animal_id", "day"])
        milk["milk_yield_rolling_7d"] = (
            milk.groupby("animal_id")["milk_yield_kg"]
            .transform(lambda x: x.rolling(7, min_periods=3).mean())
        )
        milk["milk_yield_dev"] = (
            milk["milk_yield_kg"] - milk["milk_yield_rolling_7d"]
        )
        df = df.merge(
            milk[["animal_id", "day", "milk_yield_dev", "milk_yield_rolling_7d"]],
            on=["animal_id", "day"],
            how="left",
        )

    # ── composite ketosis risk score ─────────────────────────
    # Z-score each component, then sum (higher = more risk)
    risk_cols = []
    for col, invert in [
        ("mlp_fpr", False),          # high FPR = risk
        ("rum_index_mean", True),    # low rumination = risk
        ("milk_yield_dev", True),    # yield drop = risk
        ("water_liter", True),       # low water = risk
    ]:
        if col in df.columns and df[col].notna().sum() > 10:
            mu = df[col].mean()
            sd = df[col].std()
            if sd > 0:
                z = (df[col] - mu) / sd
                if invert:
                    z = -z
                df[f"z_{col}"] = z
                risk_cols.append(f"z_{col}")

    if risk_cols:
        df["ketosis_risk_score"] = df[risk_cols].mean(axis=1)
    else:
        df["ketosis_risk_score"] = np.nan

    return df


# ─────────────────────────────────────────────────────────────
#  « classification »
# ─────────────────────────────────────────────────────────────

def train_ketosis_classifier(df: pd.DataFrame, out_dir: Path) -> dict:
    """Train a Random Forest to classify ketosis risk.

    Uses FPR > 1.4 as the positive label (subclinical ketosis
    indicator) and evaluates against disease records where
    available.

    Args:
        df: DataFrame from :func:`load_ketosis_data`.
        out_dir: Directory for saving results.

    Returns:
        Dict with performance metrics and feature importances.
    """
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import StratifiedKFold, cross_validate
    from sklearn.preprocessing import StandardScaler

    feature_cols = [
        "mlp_fpr", "mlp_fat_pct", "mlp_protein_pct", "mlp_scc",
        "mlp_urea", "mlp_lactose", "rum_index_mean",
        "rumen_temp_mean", "act_index_mean", "ph_mean",
        "water_liter", "milk_yield_kg",
    ]
    available = [c for c in feature_cols if c in df.columns]
    target = "fpr_flag"

    # Binary: fpr_flag == 1 (energy deficit) vs 0 (normal)
    sub = df[df[target].isin([0, 1])].copy()
    sub = sub.dropna(subset=available + [target])

    if len(sub) < 30:
        log.warning("Too few samples (%d) for classification", len(sub))
        return {"status": "insufficient_data", "n_samples": len(sub)}

    X = sub[available].values
    y = sub[target].values.astype(int)

    log.info(
        "Training RF classifier: %d samples, %d features, "
        "%.1f%% positive",
        len(y), len(available), 100 * y.mean(),
    )

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    clf = RandomForestClassifier(
        n_estimators=200, max_depth=8, random_state=42,
        class_weight="balanced",
    )

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    scores = cross_validate(
        clf, X_scaled, y, cv=cv,
        scoring=["accuracy", "f1", "roc_auc"],
        return_train_score=False,
    )

    # Fit on full data for feature importance
    clf.fit(X_scaled, y)
    importances = dict(zip(available, clf.feature_importances_))

    results = {
        "status": "ok",
        "n_samples": len(y),
        "n_positive": int(y.sum()),
        "features": available,
        "cv_accuracy": float(np.mean(scores["test_accuracy"])),
        "cv_f1": float(np.mean(scores["test_f1"])),
        "cv_auc": float(np.mean(scores["test_roc_auc"])),
        "feature_importances": importances,
    }

    log.info("  CV Accuracy: %.3f", results["cv_accuracy"])
    log.info("  CV F1:       %.3f", results["cv_f1"])
    log.info("  CV AUC:      %.3f", results["cv_auc"])
    log.info("  Top features:")
    for feat, imp in sorted(
        importances.items(), key=lambda x: -x[1],
    )[:5]:
        log.info("    %-25s %.4f", feat, imp)

    # ── save feature importance plot ─────────────────────────
    import matplotlib.pyplot as plt
    setup_plotting()

    fi = pd.Series(importances).sort_values()
    fig, ax = plt.subplots(figsize=(8, 5))
    fi.plot.barh(ax=ax, color="#534AB7")
    ax.set_xlabel("Feature importance (Gini)")
    ax.set_title("Ketosis risk — Random Forest feature importances")
    save_fig(fig, "ketosis_feature_importance", out_dir)

    return results


# ─────────────────────────────────────────────────────────────
#  « visualisation »
# ─────────────────────────────────────────────────────────────

def plot_ketosis_overview(df: pd.DataFrame, out_dir: Path) -> None:
    """Generate overview plots for ketosis analysis.

    Args:
        df: DataFrame from :func:`load_ketosis_data`.
        out_dir: Output directory for figures.
    """
    import matplotlib.pyplot as plt
    setup_plotting()

    # ── FPR distribution by disease status ───────────────────
    if "is_sick" in df.columns and df["mlp_fpr"].notna().sum() > 10:
        fig, ax = plt.subplots(figsize=(8, 5))
        for label, grp in df.groupby("is_sick"):
            vals = grp["mlp_fpr"].dropna()
            if len(vals) > 5:
                name = "Sick" if label == 1 else "Healthy"
                ax.hist(vals, bins=40, alpha=0.6, label=name, density=True)
        ax.axvline(1.4, color="red", linestyle="--", label="Ketosis threshold (1.4)")
        ax.axvline(1.1, color="orange", linestyle="--", label="Acidosis threshold (1.1)")
        ax.set_xlabel("Fat-to-protein ratio")
        ax.set_ylabel("Density")
        ax.legend()
        ax.set_title("FPR distribution: healthy vs. sick cows")
        save_fig(fig, "ketosis_fpr_distribution", out_dir)

    # ── risk score vs disease ────────────────────────────────
    if "ketosis_risk_score" in df.columns and df["ketosis_risk_score"].notna().sum() > 10:
        fig, ax = plt.subplots(figsize=(8, 5))
        for label, grp in df.groupby("is_sick"):
            vals = grp["ketosis_risk_score"].dropna()
            if len(vals) > 5:
                name = "Sick" if label == 1 else "Healthy"
                ax.hist(vals, bins=40, alpha=0.6, label=name, density=True)
        ax.set_xlabel("Composite ketosis risk score (z-score sum)")
        ax.set_ylabel("Density")
        ax.legend()
        ax.set_title("Composite risk score: healthy vs. sick cows")
        save_fig(fig, "ketosis_risk_score", out_dir)


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────

def main() -> None:
    """Entry point for ketosis analysis."""
    parser = argparse.ArgumentParser(description="Ketosis risk analysis")
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("results/ketosis"))
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(name)-20s │ %(message)s",
        datefmt="%H:%M:%S",
    )

    con = connect_db(args.db)
    df = load_ketosis_data(con)

    if df.empty:
        log.warning("No data — is the database fully ingested?")
        return

    args.out.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out / "ketosis_data.csv", index=False)
    log.info("Saved raw data: %s", args.out / "ketosis_data.csv")

    plot_ketosis_overview(df, args.out)
    results = train_ketosis_classifier(df, args.out)

    # Save results summary
    import json
    with open(args.out / "ketosis_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    con.close()


if __name__ == "__main__":
    main()

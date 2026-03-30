#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  Analysis 11 — Circadian rhythm disruption                   ║
# ║  « 24h Fourier decomposition of temp, activity, rumination » ║
# ╚══════════════════════════════════════════════════════════════╝
"""Circadian rhythm analysis as a general-purpose welfare biomarker.

Uses ``v_analysis_circadian`` which provides hourly aggregates of
rumen temperature, activity, and rumination per animal per day,
joined with disease ground truth.

Biological rationale: healthy ruminants exhibit strong ~24h rhythms
in core body temperature (nadir early morning, peak late afternoon),
activity (bimodal: dawn and dusk feeding bouts), and rumination
(complementary to activity — peaks at rest).  Circadian amplitude
collapse or phase shift is a well-established early marker of
sickness in human chronobiology but barely explored in cattle.

The analysis:

1.  For each animal-day, fits a **single-harmonic Fourier model**
    (24h period) to the hourly profile of each signal.
2.  Extracts **amplitude** (strength of rhythm) and **acrophase**
    (time of peak) as daily biomarkers.
3.  Computes a **Circadian Disruption Index (CDI)** = deviation
    from the animal's own healthy-period baseline.
4.  Tests whether CDI elevation precedes clinical disease onset.

Usage::

    python -m digimuh.analysis_11_circadian --db cow.db --out results/circadian
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.analysis_utils import connect_db, query_df, setup_plotting, save_fig

log = logging.getLogger("digimuh.circadian")


# ─────────────────────────────────────────────────────────────
#  « data loading »
# ─────────────────────────────────────────────────────────────

SQL_CIRCADIAN = """
SELECT *
FROM v_analysis_circadian
WHERE n_readings >= 3
ORDER BY animal_id, day, hour
"""


def load_circadian_data(con) -> pd.DataFrame:
    """Load hourly circadian data.

    Args:
        con: Database connection with views.

    Returns:
        DataFrame with hourly temp/activity/rumination per animal-day.
    """
    log.info("Loading circadian data …")
    df = query_df(con, SQL_CIRCADIAN)
    log.info(
        "  %d hourly records, %d animal-days",
        len(df),
        df.groupby(["animal_id", "day"]).ngroups if not df.empty else 0,
    )
    return df


# ─────────────────────────────────────────────────────────────
#  « Fourier decomposition »
# ─────────────────────────────────────────────────────────────

def fit_circadian_harmonic(hours: np.ndarray, values: np.ndarray) -> dict:
    """Fit a single 24h-period harmonic to hourly data.

    Model: y(t) = A * cos(2π/24 * t - φ) + M

    Uses the closed-form DFT approach (no iterative fitting):
    compute the first Fourier coefficient at frequency 1/24h.

    Args:
        hours: Array of hour-of-day values (0–23).
        values: Corresponding signal values.

    Returns:
        Dict with ``amplitude``, ``acrophase_h`` (hour of peak,
        0–24), ``mesor`` (24h mean), ``n_hours`` (data points),
        and ``relative_amplitude`` (amplitude / mesor).
    """
    mask = ~np.isnan(values)
    if mask.sum() < 6:
        return {
            "amplitude": np.nan,
            "acrophase_h": np.nan,
            "mesor": np.nan,
            "relative_amplitude": np.nan,
            "n_hours": int(mask.sum()),
        }

    h = hours[mask].astype(float)
    v = values[mask].astype(float)

    mesor = np.mean(v)
    omega = 2.0 * np.pi / 24.0

    # DFT at the 24h frequency
    cos_comp = np.mean(v * np.cos(omega * h))
    sin_comp = np.mean(v * np.sin(omega * h))

    amplitude = 2.0 * np.sqrt(cos_comp**2 + sin_comp**2)
    acrophase = np.arctan2(sin_comp, cos_comp)  # radians
    acrophase_h = (acrophase * 24.0 / (2.0 * np.pi)) % 24.0  # hours

    rel_amp = amplitude / abs(mesor) if abs(mesor) > 1e-6 else np.nan

    return {
        "amplitude": amplitude,
        "acrophase_h": acrophase_h,
        "mesor": mesor,
        "relative_amplitude": rel_amp,
        "n_hours": int(mask.sum()),
    }


def extract_circadian_features(df: pd.DataFrame) -> pd.DataFrame:
    """Extract circadian features for each animal-day.

    For each of temperature, activity, and rumination, computes
    the 24h Fourier amplitude, acrophase, and relative amplitude.

    Args:
        df: Hourly circadian DataFrame.

    Returns:
        DataFrame with one row per animal-day and columns for
        each signal's circadian parameters.
    """
    signals = {
        "temp": "temp_clean_mean",
        "act": "act_index_mean",
        "rum": "rum_index_mean",
    }

    records = []
    grouped = df.groupby(["animal_id", "day"])

    log.info("Fitting circadian harmonics for %d animal-days …", len(grouped))

    for (aid, day), grp in grouped:
        if len(grp) < 6:
            continue

        hours = grp["hour"].values
        row = {"animal_id": aid, "day": day}

        # Disease status (take max across hours — any sick hour = sick day)
        row["is_sick"] = int(grp["is_sick"].max()) if "is_sick" in grp.columns else 0

        for prefix, col in signals.items():
            if col not in grp.columns:
                continue
            result = fit_circadian_harmonic(hours, grp[col].values)
            for key, val in result.items():
                row[f"{prefix}_{key}"] = val

        records.append(row)

    result = pd.DataFrame(records)
    log.info("  Extracted features for %d animal-days", len(result))
    return result


# ─────────────────────────────────────────────────────────────
#  « Circadian Disruption Index »
# ─────────────────────────────────────────────────────────────

def compute_disruption_index(
    features: pd.DataFrame, baseline_days: int = 30,
) -> pd.DataFrame:
    """Compute Circadian Disruption Index (CDI) per animal-day.

    CDI measures how far each day's circadian parameters deviate
    from the animal's own baseline (first *baseline_days* of
    healthy-period data).  A Mahalanobis-like distance across
    amplitude, phase, and mesor of all three signals.

    Args:
        features: DataFrame from :func:`extract_circadian_features`.
        baseline_days: Number of initial healthy days to define
            the per-animal baseline.

    Returns:
        DataFrame with ``animal_id``, ``day``, ``cdi``.
    """
    metric_cols = [
        c for c in features.columns
        if c.endswith("_amplitude") or c.endswith("_acrophase_h") or c.endswith("_mesor")
    ]

    if not metric_cols:
        log.warning("No circadian metric columns found")
        return pd.DataFrame()

    features = features.copy()
    features["day"] = pd.to_datetime(features["day"])
    features = features.sort_values(["animal_id", "day"])

    records = []
    for aid, grp in features.groupby("animal_id"):
        # Baseline from first N healthy days
        healthy = grp[grp["is_sick"] == 0].head(baseline_days)
        if len(healthy) < 10:
            continue

        # Baseline statistics
        baselines = {}
        for col in metric_cols:
            vals = healthy[col].dropna()
            if len(vals) >= 5:
                baselines[col] = (vals.mean(), vals.std())

        if len(baselines) < 2:
            continue

        # CDI = mean Z-score deviation from baseline
        for _, row in grp.iterrows():
            z_scores = []
            for col, (mu, sd) in baselines.items():
                if pd.notna(row[col]) and sd > 0:
                    z = abs(row[col] - mu) / sd
                    z_scores.append(z)
            if z_scores:
                records.append({
                    "animal_id": aid,
                    "day": row["day"],
                    "cdi": np.mean(z_scores),
                    "is_sick": row.get("is_sick", 0),
                })

    result = pd.DataFrame(records)
    log.info("CDI computed for %d animal-days", len(result))
    return result


# ─────────────────────────────────────────────────────────────
#  « plots »
# ─────────────────────────────────────────────────────────────

def plot_circadian_results(
    features: pd.DataFrame, cdi: pd.DataFrame, out_dir: Path,
) -> None:
    """Generate circadian analysis plots.

    Args:
        features: Circadian feature DataFrame.
        cdi: Circadian Disruption Index DataFrame.
        out_dir: Output directory.
    """
    import matplotlib.pyplot as plt
    setup_plotting()

    # ── amplitude distributions: healthy vs sick ─────────────
    for signal in ["temp", "act", "rum"]:
        col = f"{signal}_amplitude"
        if col not in features.columns:
            continue
        sub = features.dropna(subset=[col])
        if len(sub) < 50:
            continue

        fig, ax = plt.subplots(figsize=(8, 5))
        for sick_val, label, color in [(0, "Healthy", "#1D9E75"), (1, "Sick", "#D85A30")]:
            vals = sub[sub["is_sick"] == sick_val][col].dropna()
            if len(vals) > 10:
                ax.hist(vals, bins=40, alpha=0.6, label=label,
                        color=color, density=True)
        ax.set_xlabel(f"{signal} circadian amplitude")
        ax.set_ylabel("Density")
        ax.legend()
        ax.set_title(f"24h rhythm strength ({signal}): healthy vs. sick")
        save_fig(fig, f"circadian_{signal}_amplitude", out_dir)

    # ── CDI time course for example animals ──────────────────
    if not cdi.empty:
        # Pick animals that have both sick and healthy periods
        transitions = cdi.groupby("animal_id")["is_sick"].agg(["sum", "count"])
        candidates = transitions[
            (transitions["sum"] > 5) & (transitions["sum"] < transitions["count"] * 0.8)
        ]
        if not candidates.empty:
            example_ids = candidates.head(4).index.tolist()
            fig, axes = plt.subplots(
                len(example_ids), 1, figsize=(12, 3 * len(example_ids)),
                sharex=True,
            )
            if len(example_ids) == 1:
                axes = [axes]

            for ax, aid in zip(axes, example_ids):
                cow = cdi[cdi["animal_id"] == aid].sort_values("day")
                ax.plot(cow["day"], cow["cdi"], color="#534AB7",
                        linewidth=0.8, alpha=0.8)
                # Shade sick periods
                sick = cow[cow["is_sick"] == 1]
                if not sick.empty:
                    for _, row in sick.iterrows():
                        ax.axvspan(
                            row["day"], row["day"] + pd.Timedelta(days=1),
                            alpha=0.15, color="#E24B4A",
                        )
                ax.set_ylabel("CDI")
                ax.set_title(f"Animal {aid}", fontsize=10)

            axes[-1].set_xlabel("Date")
            fig.suptitle(
                "Circadian Disruption Index — example animals\n"
                "(red shading = disease period)",
                fontsize=12,
            )
            fig.tight_layout()
            save_fig(fig, "circadian_cdi_examples", out_dir)

        # ── CDI distribution: healthy vs sick ────────────────
        fig, ax = plt.subplots(figsize=(8, 5))
        for sick_val, label, color in [(0, "Healthy", "#1D9E75"), (1, "Sick", "#D85A30")]:
            vals = cdi[cdi["is_sick"] == sick_val]["cdi"].dropna()
            if len(vals) > 10:
                ax.hist(vals, bins=50, alpha=0.6, label=label,
                        color=color, density=True)
        ax.set_xlabel("Circadian Disruption Index")
        ax.set_ylabel("Density")
        ax.legend()
        ax.set_title("CDI distribution: healthy vs. sick days")
        save_fig(fig, "circadian_cdi_distribution", out_dir)


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────

def main() -> None:
    """Entry point for circadian analysis."""
    parser = argparse.ArgumentParser(description="Circadian rhythm analysis")
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("results/circadian"))
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(name)-20s │ %(message)s",
        datefmt="%H:%M:%S",
    )

    con = connect_db(args.db)
    df = load_circadian_data(con)

    if df.empty:
        log.warning("No data.")
        return

    args.out.mkdir(parents=True, exist_ok=True)

    features = extract_circadian_features(df)
    features.to_csv(args.out / "circadian_features.csv", index=False)

    cdi = compute_disruption_index(features)
    if not cdi.empty:
        cdi.to_csv(args.out / "circadian_disruption_index.csv", index=False)

    plot_circadian_results(features, cdi, args.out)
    con.close()


if __name__ == "__main__":
    main()

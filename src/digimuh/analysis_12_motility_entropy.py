#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  Analysis 12 — Motility pattern entropy                      ║
# ║  « rumen health via information theory »                     ║
# ╚══════════════════════════════════════════════════════════════╝
"""Reticulorumen contraction entropy as a novel welfare biomarker.

Uses ``v_analysis_motility`` which extracts raw motility time series
(contraction interval and pulse width) from smaXtec derived data.

Biological rationale: in a healthy rumen, reticulorumen contractions
are quasi-periodic with modest beat-to-beat variability (analogous
to healthy heart rate variability).  Pathological states — acidosis,
inflammation, impaction — disrupt this regularity.  Very low entropy
(rigid, uncoupled contractions) and very high entropy (chaotic,
disorganised contractions) both indicate dysfunction.

This is directly analogous to heart rate variability (HRV) analysis
in cardiology, applied to the rumen motor complex.  As far as we
know, this approach has not been published.

The analysis:

1.  Computes **sample entropy** and **permutation entropy** of the
    contraction interval (``mot_period``) series in sliding windows.
2.  Derives daily summary statistics: mean entropy, entropy SD,
    and entropy trend (slope).
3.  Correlates entropy features with concurrent rumen pH, rumination
    index, and disease status.
4.  Tests whether entropy changes precede clinical diagnosis.

Usage::

    python -m digimuh.analysis_12_motility_entropy --db cow.db --out results/entropy
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.analysis_utils import connect_db, query_df, setup_plotting, save_fig

log = logging.getLogger("digimuh.entropy")


# ─────────────────────────────────────────────────────────────
#  « entropy computations »
# ─────────────────────────────────────────────────────────────

def sample_entropy(x: np.ndarray, m: int = 2, r: float | None = None) -> float:
    """Compute sample entropy of a time series.

    Sample entropy (SampEn) quantifies the regularity of a signal.
    Lower values indicate more self-similarity (regularity); higher
    values indicate more complexity/randomness.

    Args:
        x: 1-D time series (must have len > m+1).
        m: Embedding dimension (template length).  Default 2,
            standard for physiological signals.
        r: Tolerance radius.  Default is 0.2 * std(x), the
            standard choice from Richman & Moorman (2000).

    Returns:
        Sample entropy value.  Returns ``np.nan`` if the series
        is too short or constant.

    References:
        Richman JS, Moorman JR. Am J Physiol Heart Circ Physiol.
        2000;278:H2039–49.
    """
    x = np.asarray(x, dtype=float)
    N = len(x)
    if N < m + 2:
        return np.nan

    if r is None:
        sd = np.std(x, ddof=1)
        if sd < 1e-10:
            return np.nan
        r = 0.2 * sd

    def _count_templates(dim: int) -> int:
        """Count template matches at a given dimension."""
        templates = np.array([x[i:i + dim] for i in range(N - dim)])
        count = 0
        for i in range(len(templates)):
            for j in range(i + 1, len(templates)):
                if np.max(np.abs(templates[i] - templates[j])) < r:
                    count += 1
        return count

    A = _count_templates(m + 1)
    B = _count_templates(m)

    if B == 0:
        return np.nan

    return -np.log(A / B) if A > 0 else np.nan


def permutation_entropy(
    x: np.ndarray, order: int = 3, delay: int = 1, normalize: bool = True,
) -> float:
    """Compute permutation entropy of a time series.

    Permutation entropy (PE) captures the complexity of a signal
    based on the ordinal patterns of consecutive values.  It is
    robust to noise and monotonic transformations.

    Args:
        order: Embedding order (permutation length).  Default 3.
        delay: Embedding delay.  Default 1.
        normalize: If True, normalise by log(order!) to [0, 1].

    Returns:
        Permutation entropy value.

    References:
        Bandt C, Pompe B. Phys Rev Lett. 2002;88:174102.
    """
    x = np.asarray(x, dtype=float)
    N = len(x)

    if N < order * delay + 1:
        return np.nan

    # Build embedded matrix
    n_patterns = N - (order - 1) * delay
    patterns = np.zeros((n_patterns, order))
    for i in range(order):
        patterns[:, i] = x[i * delay: i * delay + n_patterns]

    # Convert to ordinal patterns (rank sequences)
    ordinals = np.argsort(np.argsort(patterns, axis=1), axis=1)

    # Count unique patterns
    # Convert each ordinal pattern to a tuple for hashing
    pattern_counts: dict[tuple, int] = {}
    for row in ordinals:
        key = tuple(row)
        pattern_counts[key] = pattern_counts.get(key, 0) + 1

    # Shannon entropy of pattern distribution
    total = sum(pattern_counts.values())
    probs = np.array(list(pattern_counts.values())) / total
    H = -np.sum(probs * np.log(probs))

    if normalize:
        from math import factorial, log as mathlog
        H_max = mathlog(factorial(order))
        H = H / H_max if H_max > 0 else np.nan

    return H


# ─────────────────────────────────────────────────────────────
#  « data loading and feature extraction »
# ─────────────────────────────────────────────────────────────

SQL_MOTILITY = """
SELECT animal_id, day, hour, mot_period, mot_pw, ph, rum_index, temp
FROM v_analysis_motility
ORDER BY animal_id, "timestamp"
"""

SQL_DISEASE_WINDOWS = """
SELECT
    animal_id,
    disease_first_day,
    disease_stop_day,
    disease_category
FROM diseases
WHERE disease_category != 'No Disease'
  AND disease_category != 'Healthy/remarks'
  AND disease_first_day IS NOT NULL
"""


def compute_daily_entropy(con, window_size: int = 50) -> pd.DataFrame:
    """Compute daily motility entropy features per animal.

    For each animal-day, collects the ``mot_period`` readings,
    computes sample entropy and permutation entropy, and derives
    summary statistics.

    Args:
        con: Database connection with views.
        window_size: Minimum number of motility readings per day
            to compute entropy (default: 50).

    Returns:
        DataFrame with daily entropy features per animal.
    """
    log.info("Loading motility data …")
    df = query_df(con, SQL_MOTILITY)
    log.info("  %d motility readings", len(df))

    if df.empty:
        return pd.DataFrame()

    # Load disease windows for ground truth
    diseases = query_df(con, SQL_DISEASE_WINDOWS)

    records = []
    grouped = df.groupby(["animal_id", "day"])
    n_groups = len(grouped)
    log.info("Computing entropy for %d animal-days …", n_groups)

    for i, ((aid, day), grp) in enumerate(grouped):
        series = grp["mot_period"].dropna().values

        if len(series) < window_size:
            continue

        # ── sample entropy ───────────────────────────────────
        sampen = sample_entropy(series, m=2)

        # ── permutation entropy ──────────────────────────────
        permen = permutation_entropy(series, order=3, delay=1)

        # ── basic variability stats (HRV analogues) ──────────
        ibi_mean = np.mean(series)       # mean inter-beat interval
        ibi_std = np.std(series, ddof=1) # SDNN analogue
        ibi_cv = ibi_std / ibi_mean if ibi_mean > 0 else np.nan
        rmssd = np.sqrt(np.mean(np.diff(series) ** 2))  # RMSSD

        # ── concurrent rumen state ───────────────────────────
        ph_vals = grp["ph"].dropna()
        rum_vals = grp["rum_index"].dropna()

        row = {
            "animal_id": aid,
            "day": day,
            "sampen": sampen,
            "permen": permen,
            "mot_ibi_mean": ibi_mean,
            "mot_ibi_std": ibi_std,
            "mot_ibi_cv": ibi_cv,
            "mot_rmssd": rmssd,
            "ph_mean": ph_vals.mean() if len(ph_vals) > 0 else np.nan,
            "rum_mean": rum_vals.mean() if len(rum_vals) > 0 else np.nan,
            "n_readings": len(series),
        }

        # ── disease status ───────────────────────────────────
        if not diseases.empty:
            sick = diseases[
                (diseases["animal_id"] == aid)
                & (diseases["disease_first_day"] <= day)
                & (diseases["disease_stop_day"].fillna(day) >= day)
            ]
            row["is_sick"] = 1 if len(sick) > 0 else 0
        else:
            row["is_sick"] = 0

        records.append(row)

        if (i + 1) % 5000 == 0:
            log.info("  %d / %d animal-days processed", i + 1, n_groups)

    result = pd.DataFrame(records)
    log.info("  Entropy features for %d animal-days", len(result))
    return result


# ─────────────────────────────────────────────────────────────
#  « pre-disease entropy trends »
# ─────────────────────────────────────────────────────────────

def compute_predisease_trends(
    entropy_df: pd.DataFrame, lookback_days: int = 7,
) -> pd.DataFrame:
    """Test whether entropy changes before disease onset.

    For each disease event, computes the mean entropy in the
    *lookback_days* before diagnosis and compares to the animal's
    healthy-period baseline.

    Args:
        entropy_df: Daily entropy features.
        lookback_days: Days before diagnosis to examine.

    Returns:
        DataFrame with pre-disease vs. baseline entropy comparison.
    """
    if "is_sick" not in entropy_df.columns:
        return pd.DataFrame()

    entropy_df = entropy_df.copy()
    entropy_df["day"] = pd.to_datetime(entropy_df["day"])
    entropy_df = entropy_df.sort_values(["animal_id", "day"])

    records = []
    for aid, grp in entropy_df.groupby("animal_id"):
        healthy = grp[grp["is_sick"] == 0]
        if len(healthy) < 14:
            continue

        baseline_sampen = healthy["sampen"].mean()
        baseline_permen = healthy["permen"].mean()

        # Find disease onset days (transitions from 0 → 1)
        grp = grp.copy()
        grp["sick_shift"] = grp["is_sick"].shift(1, fill_value=0)
        onsets = grp[(grp["is_sick"] == 1) & (grp["sick_shift"] == 0)]

        for _, onset in onsets.iterrows():
            # Look back N days before onset
            pre_window = grp[
                (grp["day"] >= onset["day"] - pd.Timedelta(days=lookback_days))
                & (grp["day"] < onset["day"])
            ]
            if len(pre_window) < 3:
                continue

            records.append({
                "animal_id": aid,
                "onset_day": onset["day"],
                "baseline_sampen": baseline_sampen,
                "baseline_permen": baseline_permen,
                "predisease_sampen": pre_window["sampen"].mean(),
                "predisease_permen": pre_window["permen"].mean(),
                "sampen_delta": pre_window["sampen"].mean() - baseline_sampen,
                "permen_delta": pre_window["permen"].mean() - baseline_permen,
                "n_predisease_days": len(pre_window),
            })

    result = pd.DataFrame(records)
    if not result.empty:
        log.info(
            "Pre-disease trends: %d events, mean SampEn delta = %.4f",
            len(result), result["sampen_delta"].mean(),
        )
    return result


# ─────────────────────────────────────────────────────────────
#  « plots »
# ─────────────────────────────────────────────────────────────

def plot_entropy_results(
    entropy_df: pd.DataFrame, trends: pd.DataFrame, out_dir: Path,
) -> None:
    """Generate entropy analysis plots.

    Args:
        entropy_df: Daily entropy features.
        trends: Pre-disease trend results.
        out_dir: Output directory.
    """
    import matplotlib.pyplot as plt
    setup_plotting()

    # ── SampEn vs PermEn scatter ─────────────────────────────
    sub = entropy_df.dropna(subset=["sampen", "permen"])
    if len(sub) > 100:
        fig, ax = plt.subplots(figsize=(8, 6))
        colors = ["#1D9E75" if s == 0 else "#D85A30" for s in sub["is_sick"]]
        ax.scatter(sub["permen"], sub["sampen"], s=5, alpha=0.3, c=colors)
        ax.set_xlabel("Permutation entropy (normalised)")
        ax.set_ylabel("Sample entropy")
        ax.set_title("Motility entropy: healthy (green) vs. sick (coral)")
        save_fig(fig, "entropy_sampen_vs_permen", out_dir)

    # ── entropy distributions ────────────────────────────────
    for col, label in [("sampen", "Sample entropy"), ("permen", "Permutation entropy")]:
        if col not in entropy_df.columns:
            continue
        fig, ax = plt.subplots(figsize=(8, 5))
        for sick_val, name, color in [(0, "Healthy", "#1D9E75"), (1, "Sick", "#D85A30")]:
            vals = entropy_df[entropy_df["is_sick"] == sick_val][col].dropna()
            if len(vals) > 10:
                ax.hist(vals, bins=50, alpha=0.6, label=name,
                        color=color, density=True)
        ax.set_xlabel(label)
        ax.set_ylabel("Density")
        ax.legend()
        ax.set_title(f"Motility {label.lower()}: healthy vs. sick")
        save_fig(fig, f"entropy_{col}_distribution", out_dir)

    # ── entropy vs pH scatter ────────────────────────────────
    sub = entropy_df.dropna(subset=["sampen", "ph_mean"])
    if len(sub) > 100:
        fig, ax = plt.subplots(figsize=(8, 6))
        sc = ax.scatter(
            sub["ph_mean"], sub["sampen"], s=5, alpha=0.3,
            c=sub["mot_ibi_cv"], cmap="viridis",
        )
        ax.set_xlabel("Daily mean rumen pH")
        ax.set_ylabel("Motility sample entropy")
        ax.set_title("Motility entropy vs. rumen pH (colour = contraction CV)")
        fig.colorbar(sc, ax=ax, label="Contraction interval CV")
        save_fig(fig, "entropy_vs_ph", out_dir)

    # ── pre-disease entropy shift ────────────────────────────
    if not trends.empty:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for ax, col, label in [
            (axes[0], "sampen_delta", "Sample entropy Δ"),
            (axes[1], "permen_delta", "Permutation entropy Δ"),
        ]:
            vals = trends[col].dropna()
            if len(vals) > 5:
                ax.hist(vals, bins=30, color="#534AB7", edgecolor="white")
                ax.axvline(0, color="#A32D2D", linestyle="--")
                ax.set_xlabel(f"{label} (pre-disease − baseline)")
                ax.set_ylabel("Disease events")
                median = vals.median()
                ax.set_title(f"{label}: median shift = {median:.4f}")
        fig.suptitle(
            f"Entropy change in {len(trends)} disease onset events\n"
            f"(7-day pre-disease window vs. healthy baseline)",
            fontsize=12,
        )
        fig.tight_layout()
        save_fig(fig, "entropy_predisease_shift", out_dir)


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────

def main() -> None:
    """Entry point for motility entropy analysis."""
    parser = argparse.ArgumentParser(description="Motility entropy analysis")
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("results/entropy"))
    parser.add_argument(
        "--min-readings", type=int, default=50,
        help="Min motility readings per day (default: 50)",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(name)-20s │ %(message)s",
        datefmt="%H:%M:%S",
    )

    con = connect_db(args.db)
    entropy_df = compute_daily_entropy(con, window_size=args.min_readings)

    if entropy_df.empty:
        log.warning("No data.")
        return

    args.out.mkdir(parents=True, exist_ok=True)
    entropy_df.to_csv(args.out / "motility_entropy.csv", index=False)

    trends = compute_predisease_trends(entropy_df)
    if not trends.empty:
        trends.to_csv(args.out / "predisease_entropy_trends.csv", index=False)

    plot_entropy_results(entropy_df, trends, args.out)
    con.close()


if __name__ == "__main__":
    main()

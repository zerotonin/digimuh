#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  Analysis 0 — Broken-stick regression ¹                      ║
# ║  « individual thermoregulatory breakpoints »                 ║
# ╚══════════════════════════════════════════════════════════════╝
"""Piecewise linear (broken-stick) regression of rumen temperature
against barn Temperature-Humidity Index (THI) and barn air temperature.

Implements the workflow described by Hoffmann et al. (Frontiers
manuscript, 2026) for the DigiMuh dataset.  For each animal in the
collaborator-provided selection list (``Tierauswahl.xlsx``), the script:

1.  Extracts rumen temperature (``temp_without_drink_cycles``) from the
    smaXtec derived table.
2.  Joins with barn climate data (``smaxtec_barns``) to obtain
    concurrent barn air temperature and THI.
3.  Excludes milking hours (04:00–07:59 and 16:00–19:59) when cows
    are not in the barn.
4.  Fits a two-segment piecewise linear regression per animal to find
    the **breakpoint** — the environmental value at which rumen
    temperature begins to rise, indicating the onset of heat stress.
5.  Reports per-animal breakpoints and generates boxplot summaries
    across years.

The breakpoint represents the individual upper critical threshold:
below it, rumen temperature is constant (thermoneutral zone); above it,
rumen temperature increases linearly with environmental load.

Usage::

    python -m digimuh.analysis_00_broken_stick \\
        --db cow.db \\
        --tierauswahl Tierauswahl.xlsx \\
        --out results/broken_stick

¹ Analysis led by Dr. med. vet. Gundula Hoffmann, Head of working group
  "Digital monitoring of animal welfare", Leibniz Institute for
  Agricultural Engineering and Bioeconomy (ATB), Potsdam, Germany.
  https://www.atb-potsdam.de/en/

References:
    NRC (1971) — THI formula: THI = (1.8×AT+32) − (0.55−0.0055×RH) × (1.8×AT−26)
    Neira et al. (2026) PLOS Climate — thermal stress duration/load projections
    Hoffmann et al. (2020) — animal-based heat stress indicators
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from digimuh.analysis_utils import connect_db, query_df, setup_plotting, save_fig

log = logging.getLogger("digimuh.broken_stick")


# ─────────────────────────────────────────────────────────────
#  « broken-stick regression model »
# ─────────────────────────────────────────────────────────────

def broken_stick_fit(
    x: np.ndarray, y: np.ndarray,
    x_range: tuple[float, float] | None = None,
    n_grid: int = 200,
) -> dict:
    """Fit a two-segment piecewise linear model with one breakpoint.

    Model::

        y = { a1 + b1*x          if x <= bp
            { a1 + b1*bp + b2*(x - bp)  if x > bp

    The breakpoint *bp* is found by grid search + refinement: for each
    candidate bp, the two segments are fitted by OLS and the total
    residual sum of squares (RSS) is minimised.

    Args:
        x: Predictor values (e.g. THI or barn temperature).
        y: Response values (e.g. rumen temperature).
        x_range: (min, max) search range for the breakpoint.
            Defaults to the 10th–90th percentile of *x*.
        n_grid: Number of grid points for initial search.

    Returns:
        Dict with keys: ``breakpoint``, ``slope_below``, ``slope_above``,
        ``intercept``, ``r_squared``, ``rss``, ``n``, ``converged``.
    """
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    n = len(x)

    if n < 30:
        return {"breakpoint": np.nan, "converged": False, "n": n}

    if x_range is None:
        x_range = (np.percentile(x, 10), np.percentile(x, 90))

    if x_range[0] >= x_range[1]:
        return {"breakpoint": np.nan, "converged": False, "n": n}

    def rss_at_bp(bp: float) -> float:
        """Compute RSS for a given breakpoint."""
        left = x <= bp
        right = ~left

        if left.sum() < 10 or right.sum() < 10:
            return np.inf

        # Left segment: y = a1 + b1 * x
        X_left = np.column_stack([np.ones(left.sum()), x[left]])
        try:
            coef_left, _, _, _ = np.linalg.lstsq(X_left, y[left], rcond=None)
        except np.linalg.LinAlgError:
            return np.inf

        # Right segment: y = a2 + b2 * x
        X_right = np.column_stack([np.ones(right.sum()), x[right]])
        try:
            coef_right, _, _, _ = np.linalg.lstsq(X_right, y[right], rcond=None)
        except np.linalg.LinAlgError:
            return np.inf

        pred_left = X_left @ coef_left
        pred_right = X_right @ coef_right

        return np.sum((y[left] - pred_left) ** 2) + np.sum((y[right] - pred_right) ** 2)

    # ── grid search ──────────────────────────────────────────
    grid = np.linspace(x_range[0], x_range[1], n_grid)
    rss_vals = np.array([rss_at_bp(bp) for bp in grid])
    best_idx = np.argmin(rss_vals)

    if not np.isfinite(rss_vals[best_idx]):
        return {"breakpoint": np.nan, "converged": False, "n": n}

    # ── refine with bounded scalar minimisation ──────────────
    lo = grid[max(0, best_idx - 2)]
    hi = grid[min(len(grid) - 1, best_idx + 2)]
    result = minimize_scalar(rss_at_bp, bounds=(lo, hi), method="bounded")
    bp = result.x
    final_rss = result.fun

    # ── extract final coefficients ───────────────────────────
    left = x <= bp
    right = ~left

    X_left = np.column_stack([np.ones(left.sum()), x[left]])
    coef_left, _, _, _ = np.linalg.lstsq(X_left, y[left], rcond=None)

    X_right = np.column_stack([np.ones(right.sum()), x[right]])
    coef_right, _, _, _ = np.linalg.lstsq(X_right, y[right], rcond=None)

    # R² against global mean
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_sq = 1 - final_rss / ss_tot if ss_tot > 0 else np.nan

    # ── biological plausibility check ───────────────────────
    # Heat stress means body temp rises FASTER above the
    # breakpoint.  If the above-slope is negative or flatter
    # than the below-slope, the fit is spurious (e.g. peaked
    # shape from sparse data at the high end).
    slope_below = coef_left[1]
    slope_above = coef_right[1]
    biologically_valid = slope_above > slope_below and slope_above > 0

    return {
        "breakpoint": bp if biologically_valid else np.nan,
        "slope_below": slope_below,
        "intercept_below": coef_left[0],
        "slope_above": slope_above,
        "intercept_above": coef_right[0],
        "r_squared": r_sq,
        "rss": final_rss,
        "n": n,
        "n_below": int(left.sum()),
        "n_above": int(right.sum()),
        "converged": biologically_valid,
        "rejected_reason": (
            None if biologically_valid
            else f"slope_above ({slope_above:.4f}) <= slope_below ({slope_below:.4f})"
        ),
    }


# ─────────────────────────────────────────────────────────────
#  « data loading »
# ─────────────────────────────────────────────────────────────

def load_tierauswahl(path: Path) -> pd.DataFrame:
    """Load the collaborator-provided animal selection list.

    Args:
        path: Path to ``Tierauswahl.xlsx``.

    Returns:
        DataFrame with columns: ``animal_id``, ``datetime_enter``,
        ``datetime_exit``, ``group``, ``year``.
    """
    # Try different possible sheet names
    xls = pd.ExcelFile(path)
    for sheet_name in xls.sheet_names:
        df = pd.read_excel(xls, sheet_name=sheet_name)
        if "animal_id" in df.columns:
            break

    df = df[df["Auswahl"] == "Ja"].copy()
    df["animal_id"] = pd.to_numeric(df["animal_id"], errors="coerce")
    df = df.dropna(subset=["animal_id"])
    df["animal_id"] = df["animal_id"].astype(int)
    df["datetime_enter"] = pd.to_datetime(df["datetime_enter"])
    df["datetime_exit"] = pd.to_datetime(df["datetime_exit"])

    if "Versuchsjahr" in df.columns:
        df["year"] = pd.to_numeric(df["Versuchsjahr"], errors="coerce")
    else:
        df["year"] = df["datetime_enter"].dt.year

    log.info(
        "Tierauswahl: %d animal-year entries across years %s",
        len(df),
        sorted(df["year"].dropna().unique().astype(int)),
    )
    return df


SQL_RUMEN = """
SELECT
    animal_id,
    "timestamp",
    CAST("temp_without_drink_cycles" AS REAL) AS body_temp,
    CAST("drink_cycles_v2" AS REAL) AS drink_cycles
FROM smaxtec_derived
WHERE animal_id = ?
  AND "timestamp" >= ?
  AND "timestamp" <= ?
  AND "temp_without_drink_cycles" IS NOT NULL
  AND CAST("temp_without_drink_cycles" AS REAL) > 30
  AND CAST("temp_without_drink_cycles" AS REAL) < 43
"""

SQL_BARN = """
SELECT
    "timestamp",
    AVG(temp) AS barn_temp,
    AVG(temp_hum_index) AS barn_thi
FROM smaxtec_barns
WHERE "timestamp" >= ?
  AND "timestamp" <= ?
  AND temp IS NOT NULL
  AND temp_hum_index IS NOT NULL
GROUP BY "timestamp"
"""

DRINK_PAD = pd.Timedelta(minutes=15)


def _exclude_drinking_windows(df: pd.DataFrame) -> pd.DataFrame:
    """Remove rows during and 15 minutes after detected drinking events.

    Drinking causes a cold-water artifact in rumen temperature.
    smaXtec's ``temp_without_drink_cycles`` corrects the value, but
    the post-drinking recovery period still shows depressed temperatures.
    We exclude the drinking event timestamp plus a 15-minute padding
    to eliminate this residual effect.

    Args:
        df: DataFrame with ``timestamp`` and ``drink_cycles`` columns.

    Returns:
        Filtered DataFrame with drinking windows removed.
    """
    if "drink_cycles" not in df.columns:
        return df

    # Identify drinking event timestamps
    drink_times = df.loc[
        df["drink_cycles"].fillna(0) > 0, "timestamp"
    ].values

    if len(drink_times) == 0:
        return df.drop(columns=["drink_cycles"])

    # Build exclusion mask: True = keep this row
    keep = np.ones(len(df), dtype=bool)
    ts = df["timestamp"].values

    for dt in drink_times:
        # Exclude [drink_time, drink_time + 15 min]
        keep &= ~((ts >= dt) & (ts <= dt + DRINK_PAD))

    n_excluded = (~keep).sum()
    if n_excluded > 0:
        log.debug(
            "  Excluded %d readings in %d drinking windows (%.1f%%)",
            n_excluded, len(drink_times),
            100 * n_excluded / len(df),
        )

    return df.loc[keep].drop(columns=["drink_cycles"])


def load_animal_data(
    con, animal_id: int, date_enter: str, date_exit: str,
    barn_cache: dict | None = None,
) -> pd.DataFrame:
    """Load rumen + barn data for one animal, filtered for analysis.

    Queries smaxtec_derived and smaxtec_barns separately (both use
    indexes efficiently) and joins in pandas by hour-rounded timestamp.

    Filtering steps:
    1. Remove drinking episodes + 15-min post-drinking recovery window.
    2. Exclude milking hours (04:00–07:59 and 16:00–19:59).
    3. Join with barn climate by hour.

    Args:
        con: Database connection.
        animal_id: EU ear tag integer.
        date_enter: Start of observation window (ISO date).
        date_exit: End of observation window (ISO date).
        barn_cache: Optional dict mapping (enter, exit) → barn DataFrame
            to avoid re-querying barn data for the same date range.

    Returns:
        DataFrame with body_temp, barn_temp, barn_thi columns,
        drinking and milking periods excluded.
    """
    # ── rumen data (uses idx_smaxtec_derived_animal_ts) ──────
    rumen = query_df(con, SQL_RUMEN, (animal_id, date_enter, date_exit))
    if rumen.empty:
        return rumen

    rumen["timestamp"] = pd.to_datetime(rumen["timestamp"])

    # ── exclude drinking episodes + 15 min padding ───────────
    rumen = _exclude_drinking_windows(rumen)
    if rumen.empty:
        return rumen

    rumen["hour_key"] = rumen["timestamp"].dt.floor("h")

    # ── barn climate (cached across animals in same date range) ──
    cache_key = (date_enter, date_exit)
    if barn_cache is not None and cache_key in barn_cache:
        barn = barn_cache[cache_key]
    else:
        barn = query_df(con, SQL_BARN, (date_enter, date_exit))
        if not barn.empty:
            barn["timestamp"] = pd.to_datetime(barn["timestamp"])
            barn["hour_key"] = barn["timestamp"].dt.floor("h")
            # Average barn readings per hour (in case of sub-hourly data)
            barn = barn.groupby("hour_key").agg(
                barn_temp=("barn_temp", "mean"),
                barn_thi=("barn_thi", "mean"),
            ).reset_index()
        if barn_cache is not None:
            barn_cache[cache_key] = barn

    if barn.empty:
        return pd.DataFrame()

    # ── merge on hour ────────────────────────────────────────
    df = rumen.merge(barn, on="hour_key", how="inner")

    if df.empty:
        return df

    # ── exclude milking hours (04:00–07:59 and 16:00–19:59) ──
    hour = df["timestamp"].dt.hour
    df = df[~hour.between(4, 7) & ~hour.between(16, 19)]

    return df


# ─────────────────────────────────────────────────────────────
#  « main analysis loop »
# ─────────────────────────────────────────────────────────────

def run_broken_stick_analysis(
    con, tierauswahl: pd.DataFrame,
) -> pd.DataFrame:
    """Run broken-stick regression for all animals in the selection.

    Args:
        con: Database connection.
        tierauswahl: Animal selection DataFrame.

    Returns:
        DataFrame with one row per animal-year, columns for THI breakpoint,
        barn temp breakpoint, slopes, R², etc.
    """
    results = []
    total = len(tierauswahl)
    barn_cache: dict = {}  # cache barn data across animals in same date range

    for i, (_, row) in enumerate(tierauswahl.iterrows()):
        aid = int(row["animal_id"])
        enter = str(row["datetime_enter"])[:10]
        exit_ = str(row["datetime_exit"])[:10]
        year = int(row["year"]) if pd.notna(row.get("year")) else None

        if (i + 1) % 20 == 0 or i == 0:
            log.info("  [%d/%d] Processing animal %d (%s) …", i + 1, total, aid, year)

        df = load_animal_data(con, aid, enter, exit_, barn_cache=barn_cache)

        if len(df) < 50:
            log.debug("  Animal %d: only %d readings, skipping", aid, len(df))
            results.append({
                "animal_id": aid,
                "year": year,
                "n_readings": len(df),
                "thi_breakpoint": np.nan,
                "temp_breakpoint": np.nan,
                "comment": f"insufficient data ({len(df)} rows)",
            })
            continue

        # ── Fit 1: body temp vs barn THI ─────────────────────
        thi_fit = broken_stick_fit(
            df["barn_thi"].values,
            df["body_temp"].values,
            x_range=(45, 80),
        )

        # ── Fit 2: body temp vs barn temperature ─────────────
        temp_fit = broken_stick_fit(
            df["barn_temp"].values,
            df["body_temp"].values,
            x_range=(5, 35),
        )

        result = {
            "animal_id": aid,
            "year": year,
            "group": row.get("group"),
            "n_readings": len(df),
            "date_enter": enter,
            "date_exit": exit_,
            # THI breakpoint
            "thi_breakpoint": thi_fit.get("breakpoint"),
            "thi_slope_below": thi_fit.get("slope_below"),
            "thi_slope_above": thi_fit.get("slope_above"),
            "thi_r_squared": thi_fit.get("r_squared"),
            "thi_converged": thi_fit.get("converged", False),
            "thi_n_below": thi_fit.get("n_below"),
            "thi_n_above": thi_fit.get("n_above"),
            # Barn temp breakpoint
            "temp_breakpoint": temp_fit.get("breakpoint"),
            "temp_slope_below": temp_fit.get("slope_below"),
            "temp_slope_above": temp_fit.get("slope_above"),
            "temp_r_squared": temp_fit.get("r_squared"),
            "temp_converged": temp_fit.get("converged", False),
            "temp_n_below": temp_fit.get("n_below"),
            "temp_n_above": temp_fit.get("n_above"),
            "comment": "",
        }

        if not thi_fit.get("converged", False):
            reason = thi_fit.get("rejected_reason", "no breakpoint found")
            result["comment"] = f"THI: {reason}"
        if not temp_fit.get("converged", False):
            reason = temp_fit.get("rejected_reason", "no breakpoint found")
            prefix = "; " if result["comment"] else ""
            result["comment"] += f"{prefix}Temp: {reason}"

        results.append(result)

    return pd.DataFrame(results)


# ─────────────────────────────────────────────────────────────
#  « plots »
# ─────────────────────────────────────────────────────────────

def plot_breakpoint_boxplots(results: pd.DataFrame, out_dir: Path) -> None:
    """Generate boxplots of per-animal breakpoints across years.

    Args:
        results: DataFrame from :func:`run_broken_stick_analysis`.
        out_dir: Output directory.
    """
    import matplotlib.pyplot as plt
    setup_plotting()

    converged = results[results["thi_converged"] == True]
    if converged.empty:
        log.warning("No converged fits — skipping plots")
        return

    years = sorted(converged["year"].dropna().unique().astype(int))

    # ── THI breakpoint boxplot ───────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 5))
    data_thi = [
        converged[converged["year"] == y]["thi_breakpoint"].dropna().values
        for y in years
    ]
    bp = ax.boxplot(
        data_thi, labels=[str(y) for y in years],
        patch_artist=True,
        boxprops=dict(facecolor="#B7D4E8", edgecolor="#2B5E8C"),
        medianprops=dict(color="#A32D2D", linewidth=2),
    )
    ax.axhline(68.8, color="#BA7517", linestyle="--", linewidth=1,
               label="THI 68.8 (mild stress onset; Neira et al. 2026)")
    ax.set_xlabel("Year")
    ax.set_ylabel("Individual THI breakpoint")
    ax.set_title(
        "Per-animal heat stress thresholds — rumen temperature vs. barn THI\n"
        "(broken-stick regression breakpoints)",
    )
    ax.legend(fontsize=9)
    # Annotate sample sizes
    for i, y in enumerate(years):
        n = len(data_thi[i])
        ax.text(i + 1, ax.get_ylim()[0] + 0.5, f"n={n}",
                ha="center", fontsize=9, color="#555555")
    save_fig(fig, "broken_stick_thi_boxplot", out_dir)

    # ── Barn temperature breakpoint boxplot ──────────────────
    converged_temp = results[results["temp_converged"] == True]
    if not converged_temp.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        data_temp = [
            converged_temp[converged_temp["year"] == y]["temp_breakpoint"].dropna().values
            for y in years
        ]
        bp = ax.boxplot(
            data_temp, labels=[str(y) for y in years],
            patch_artist=True,
            boxprops=dict(facecolor="#D4E8B7", edgecolor="#3D6B2B"),
            medianprops=dict(color="#A32D2D", linewidth=2),
        )
        ax.set_xlabel("Year")
        ax.set_ylabel("Individual barn temperature breakpoint (°C)")
        ax.set_title(
            "Per-animal heat stress thresholds — rumen temperature vs. barn temperature\n"
            "(broken-stick regression breakpoints)",
        )
        for i, y in enumerate(years):
            n = len(data_temp[i])
            ax.text(i + 1, ax.get_ylim()[0] + 0.3, f"n={n}",
                    ha="center", fontsize=9, color="#555555")
        save_fig(fig, "broken_stick_temp_boxplot", out_dir)

    # ── Summary scatter: THI bp vs Temp bp ───────────────────
    both = results[results["thi_converged"] & results["temp_converged"]]
    if len(both) > 10:
        fig, ax = plt.subplots(figsize=(7, 6))
        colors = {2021: "#2B5E8C", 2022: "#1D9E75", 2023: "#BA7517", 2024: "#A32D2D"}
        for y in years:
            sub = both[both["year"] == y]
            if not sub.empty:
                ax.scatter(
                    sub["temp_breakpoint"], sub["thi_breakpoint"],
                    s=30, alpha=0.7, label=str(y),
                    color=colors.get(y, "#888888"),
                )
        ax.set_xlabel("Barn temperature breakpoint (°C)")
        ax.set_ylabel("THI breakpoint")
        ax.set_title("Individual heat stress thresholds: THI vs. barn temperature")
        ax.legend()
        save_fig(fig, "broken_stick_thi_vs_temp", out_dir)

    # ── Example animal: fitted curve ─────────────────────────
    # Pick the animal with the best THI R²
    best_row = converged.loc[converged["thi_r_squared"].idxmax()]
    log.info(
        "Example animal: %d (year %s, THI bp=%.1f, R²=%.3f)",
        best_row["animal_id"], best_row["year"],
        best_row["thi_breakpoint"], best_row["thi_r_squared"],
    )


def plot_example_animal(
    con, animal_id: int, date_enter: str, date_exit: str,
    thi_fit: dict, temp_fit: dict, out_dir: Path, year: int,
) -> None:
    """Plot the broken-stick fit for one example animal.

    Args:
        con: Database connection.
        animal_id: EU ear tag.
        date_enter, date_exit: Observation window.
        thi_fit, temp_fit: Fit result dicts.
        out_dir: Output directory.
        year: Observation year.
    """
    import matplotlib.pyplot as plt
    setup_plotting()

    df = load_animal_data(con, animal_id, date_enter, date_exit)
    if df.empty:
        return

    for x_col, fit, xlabel, name in [
        ("barn_thi", thi_fit, "Barn THI", "thi"),
        ("barn_temp", temp_fit, "Barn temperature (°C)", "temp"),
    ]:
        if not fit.get("converged", False):
            continue

        bp = fit["breakpoint"]
        fig, ax = plt.subplots(figsize=(9, 6))
        ax.scatter(df[x_col], df["body_temp"], s=2, alpha=0.15, c="#888780")

        # Draw fitted lines
        x_range = np.linspace(df[x_col].min(), df[x_col].max(), 200)
        y_pred = np.where(
            x_range <= bp,
            fit["intercept_below"] + fit["slope_below"] * x_range,
            fit["intercept_above"] + fit["slope_above"] * x_range,
        )
        ax.plot(x_range, y_pred, color="#A32D2D", linewidth=2, label="Broken-stick fit")
        ax.axvline(bp, color="#2B5E8C", linestyle="--", linewidth=1.5,
                   label=f"Breakpoint = {bp:.1f}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Rumen temperature (°C)")
        ax.set_title(
            f"Animal {animal_id} ({year}) — breakpoint at {xlabel} = {bp:.1f}\n"
            f"R² = {fit['r_squared']:.3f}, n = {fit['n']:,}",
        )
        ax.legend()
        save_fig(fig, f"broken_stick_example_{name}_{animal_id}_{year}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────

def main() -> None:
    """Entry point for broken-stick analysis."""
    parser = argparse.ArgumentParser(
        description="Broken-stick regression: individual heat stress thresholds",
    )
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--tierauswahl", type=Path, required=True,
                        help="Path to Tierauswahl.xlsx")
    parser.add_argument("--out", type=Path, default=Path("results/broken_stick"))
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(name)-20s │ %(message)s",
        datefmt="%H:%M:%S",
    )

    con = connect_db(args.db)
    tierauswahl = load_tierauswahl(args.tierauswahl)

    if tierauswahl.empty:
        log.warning("No animals in selection list.")
        return

    args.out.mkdir(parents=True, exist_ok=True)

    log.info("Running broken-stick regression for %d animal-year entries …", len(tierauswahl))
    results = run_broken_stick_analysis(con, tierauswahl)

    # Save results
    results.to_csv(args.out / "broken_stick_results.csv", index=False)
    log.info("Saved: %s", args.out / "broken_stick_results.csv")

    # Summary statistics
    converged_thi = results[results["thi_converged"] == True]
    converged_temp = results[results["temp_converged"] == True]
    log.info(
        "Converged: %d/%d THI fits, %d/%d barn temp fits",
        len(converged_thi), len(results),
        len(converged_temp), len(results),
    )
    if not converged_thi.empty:
        log.info(
            "THI breakpoints: median=%.1f, IQR=[%.1f, %.1f]",
            converged_thi["thi_breakpoint"].median(),
            converged_thi["thi_breakpoint"].quantile(0.25),
            converged_thi["thi_breakpoint"].quantile(0.75),
        )
    if not converged_temp.empty:
        log.info(
            "Barn temp breakpoints: median=%.1f °C, IQR=[%.1f, %.1f]",
            converged_temp["temp_breakpoint"].median(),
            converged_temp["temp_breakpoint"].quantile(0.25),
            converged_temp["temp_breakpoint"].quantile(0.75),
        )

    # Plots
    plot_breakpoint_boxplots(results, args.out)

    # Example animal plot (best THI fit)
    if not converged_thi.empty:
        best = converged_thi.loc[converged_thi["thi_r_squared"].idxmax()]
        df_ex = load_animal_data(
            con, int(best["animal_id"]),
            best["date_enter"], best["date_exit"],
        )
        if not df_ex.empty:
            thi_refit = broken_stick_fit(
                df_ex["barn_thi"].values, df_ex["body_temp"].values,
                x_range=(45, 80),
            )
            temp_refit = broken_stick_fit(
                df_ex["barn_temp"].values, df_ex["body_temp"].values,
                x_range=(5, 35),
            )
            plot_example_animal(
                con, int(best["animal_id"]),
                best["date_enter"], best["date_exit"],
                thi_refit, temp_refit, args.out,
                int(best["year"]),
            )

    con.close()
    log.info("Done.")


if __name__ == "__main__":
    main()

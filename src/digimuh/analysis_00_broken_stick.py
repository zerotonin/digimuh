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
#  « colour scheme — Wong (2011), colourblind-safe »
# ─────────────────────────────────────────────────────────────
# Reference: Wong B. (2011) Nature Methods 8:441
# All colours verified with Coblis, Sim Daltonism, and
# Color Oracle for deuteranopia, protanopia, tritanopia.

COLOURS = {
    # ── per-year palette (4 years, maximally distinct) ───────
    "year": {
        2021: "#0072B2",  # blue
        2022: "#E69F00",  # orange/amber
        2023: "#009E73",  # bluish green
        2024: "#CC79A7",  # reddish purple
    },
    # ── binary comparisons ───────────────────────────────────
    "below_bp": "#0072B2",    # blue — thermoneutral / below breakpoint
    "above_bp": "#D55E00",    # vermilion — heat stress / above breakpoint
    # ── general purpose ──────────────────────────────────────
    "healthy": "#009E73",     # bluish green
    "sick": "#D55E00",        # vermilion
    "fit_line": "#D55E00",    # vermilion — regression / fit lines
    "reference": "#E69F00",   # amber — reference thresholds (THI 68.8 etc.)
    "identity": "#999999",    # grey — identity / zero lines
    "scatter": "#56B4E9",     # sky blue — generic scatter points
    "scatter_alt": "#009E73", # bluish green — second scatter variable
    "hist_thi": "#0072B2",    # blue — THI-related histograms
    "hist_temp": "#009E73",   # green — barn temp histograms
    "box_thi": "#56B4E9",     # sky blue — THI boxplot fill
    "box_temp": "#009E73",    # green — barn temp boxplot fill
    "median": "#D55E00",      # vermilion — median lines in boxplots
    "paired_line": "#999999", # grey — paired animal lines
}


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
  AND barn_id IN (1, 2)
GROUP BY "timestamp"
"""

# Neubau barn sensor IDs — groups 1005 and 1006 are housed here.
# barn_id 1 = "NewBridge", barn_id 2 = "NewPillar".
# barn_id 3 ("OldBridge") and 4 ("OldRepro") are excluded because
# the Tierauswahl animals are not in those buildings.
NEUBAU_BARN_IDS = (1, 2)

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
        boxprops=dict(facecolor=COLOURS["box_thi"], edgecolor="#333333"),
        medianprops=dict(color=COLOURS["median"], linewidth=2),
    )
    ax.axhline(68.8, color=COLOURS["reference"], linestyle="--", linewidth=1,
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
            boxprops=dict(facecolor=COLOURS["box_temp"], edgecolor="#333333"),
            medianprops=dict(color=COLOURS["median"], linewidth=2),
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
        colors = COLOURS["year"]
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
        ax.scatter(df[x_col], df["body_temp"], s=2, alpha=0.15, c=COLOURS["identity"])

        # Draw fitted lines
        x_range = np.linspace(df[x_col].min(), df[x_col].max(), 200)
        y_pred = np.where(
            x_range <= bp,
            fit["intercept_below"] + fit["slope_below"] * x_range,
            fit["intercept_above"] + fit["slope_above"] * x_range,
        )
        ax.plot(x_range, y_pred, color=COLOURS["fit_line"], linewidth=2, label="Broken-stick fit")
        ax.axvline(bp, color=COLOURS["below_bp"], linestyle="--", linewidth=1.5,
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
#  « Spearman correlations (manuscript requirement) »
# ─────────────────────────────────────────────────────────────

def compute_spearman_correlations(
    con, tierauswahl: pd.DataFrame, barn_cache: dict | None = None,
) -> pd.DataFrame:
    """Compute per-animal Spearman correlations: body temp vs THI and barn temp.

    Args:
        con: Database connection.
        tierauswahl: Animal selection DataFrame.
        barn_cache: Shared barn data cache.

    Returns:
        DataFrame with one row per animal-year: rₛ for THI, rₛ for barn
        temp, p-values, and sample sizes.
    """
    from scipy.stats import spearmanr

    if barn_cache is None:
        barn_cache = {}

    records = []
    for _, row in tierauswahl.iterrows():
        aid = int(row["animal_id"])
        enter = str(row["datetime_enter"])[:10]
        exit_ = str(row["datetime_exit"])[:10]
        year = int(row["year"]) if pd.notna(row.get("year")) else None

        df = load_animal_data(con, aid, enter, exit_, barn_cache=barn_cache)
        if len(df) < 30:
            continue

        rec = {"animal_id": aid, "year": year, "n": len(df)}

        for x_col, prefix in [("barn_thi", "thi"), ("barn_temp", "temp")]:
            sub = df.dropna(subset=[x_col, "body_temp"])
            if len(sub) > 20:
                rs, p = spearmanr(sub[x_col], sub["body_temp"])
                rec[f"{prefix}_spearman_rs"] = rs
                rec[f"{prefix}_spearman_p"] = p
            else:
                rec[f"{prefix}_spearman_rs"] = np.nan
                rec[f"{prefix}_spearman_p"] = np.nan

        records.append(rec)

    return pd.DataFrame(records)


def plot_spearman_results(spearman: pd.DataFrame, out_dir: Path) -> None:
    """Plot Spearman correlation distributions.

    Args:
        spearman: DataFrame from :func:`compute_spearman_correlations`.
        out_dir: Output directory.
    """
    import matplotlib.pyplot as plt
    setup_plotting()

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, col, xlabel, color in [
        (axes[0], "thi_spearman_rs", "Spearman rₛ (body temp vs. barn THI)", COLOURS["hist_thi"]),
        (axes[1], "temp_spearman_rs", "Spearman rₛ (body temp vs. barn temp)", COLOURS["hist_temp"]),
    ]:
        vals = spearman[col].dropna()
        if len(vals) < 5:
            continue
        ax.hist(vals, bins=30, color=color, edgecolor="white", alpha=0.8)
        ax.axvline(vals.median(), color=COLOURS["fit_line"], linestyle="--",
                   label=f"Median = {vals.median():.3f}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Number of animals")
        ax.legend(fontsize=9)

    fig.suptitle("Per-animal Spearman correlations between rumen temperature and barn climate")
    fig.tight_layout()
    save_fig(fig, "spearman_correlations", out_dir)


# ─────────────────────────────────────────────────────────────
#  « descriptive climate summary »
# ─────────────────────────────────────────────────────────────

SQL_CLIMATE_SUMMARY = """
SELECT
    date("timestamp") AS day,
    AVG(temp) AS barn_temp_mean,
    MIN(temp) AS barn_temp_min,
    MAX(temp) AS barn_temp_max,
    AVG(hum) AS barn_rh_mean,
    AVG(temp_hum_index) AS barn_thi_mean,
    MIN(temp_hum_index) AS barn_thi_min,
    MAX(temp_hum_index) AS barn_thi_max
FROM smaxtec_barns
WHERE "timestamp" >= ? AND "timestamp" <= ?
  AND temp IS NOT NULL AND temp_hum_index IS NOT NULL
  AND barn_id IN (1, 2)
GROUP BY date("timestamp")
ORDER BY day
"""


def compute_climate_summary(con, tierauswahl: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    """Generate descriptive statistics of barn climate per summer.

    Args:
        con: Database connection.
        tierauswahl: Used to determine the date range per year.
        out_dir: Output directory for plots.

    Returns:
        DataFrame of daily barn climate summaries.
    """
    import matplotlib.pyplot as plt
    setup_plotting()

    years = sorted(tierauswahl["year"].dropna().unique().astype(int))
    all_climate = []

    for year in years:
        start = f"{year}-06-01"
        end = f"{year}-09-30"
        df = query_df(con, SQL_CLIMATE_SUMMARY, (start, end))
        if df.empty:
            continue
        df["year"] = year
        df["day"] = pd.to_datetime(df["day"])
        all_climate.append(df)

    if not all_climate:
        return pd.DataFrame()

    climate = pd.concat(all_climate, ignore_index=True)

    # ── summary table ────────────────────────────────────────
    summary = climate.groupby("year").agg(
        temp_mean=("barn_temp_mean", "mean"),
        temp_sd=("barn_temp_mean", "std"),
        temp_max_overall=("barn_temp_max", "max"),
        thi_mean=("barn_thi_mean", "mean"),
        thi_sd=("barn_thi_mean", "std"),
        thi_max_overall=("barn_thi_max", "max"),
        rh_mean=("barn_rh_mean", "mean"),
        n_days=("barn_temp_mean", "count"),
    ).reset_index()
    log.info("Climate summary per year:")
    for _, r in summary.iterrows():
        log.info(
            "  %d: Temp %.1f ± %.1f °C (max %.1f), "
            "THI %.1f ± %.1f (max %.1f), RH %.1f%%, %d days",
            r["year"], r["temp_mean"], r["temp_sd"], r["temp_max_overall"],
            r["thi_mean"], r["thi_sd"], r["thi_max_overall"],
            r["rh_mean"], r["n_days"],
        )

    # ── daily THI time series per year ───────────────────────
    fig, ax = plt.subplots(figsize=(12, 5))
    colors = COLOURS["year"]
    for year in years:
        sub = climate[climate["year"] == year]
        doy = sub["day"].dt.dayofyear
        ax.plot(doy, sub["barn_thi_mean"], linewidth=0.8, alpha=0.8,
                color=colors.get(year, "#888"), label=str(year))
        ax.fill_between(doy, sub["barn_thi_min"], sub["barn_thi_max"],
                        alpha=0.1, color=colors.get(year, "#888"))
    ax.axhline(68.8, color=COLOURS["reference"], linestyle="--", linewidth=1,
               label="THI 68.8 (mild stress)")
    ax.set_xlabel("Day of year (Jun–Sep)")
    ax.set_ylabel("Barn THI")
    ax.set_title("Daily barn THI during summer months (mean ± daily range)")
    ax.legend(fontsize=9)
    save_fig(fig, "climate_thi_timeseries", out_dir)

    # ── daily barn temperature ───────────────────────────────
    fig, ax = plt.subplots(figsize=(12, 5))
    for year in years:
        sub = climate[climate["year"] == year]
        doy = sub["day"].dt.dayofyear
        ax.plot(doy, sub["barn_temp_mean"], linewidth=0.8, alpha=0.8,
                color=colors.get(year, "#888"), label=str(year))
        ax.fill_between(doy, sub["barn_temp_min"], sub["barn_temp_max"],
                        alpha=0.1, color=colors.get(year, "#888"))
    ax.set_xlabel("Day of year (Jun–Sep)")
    ax.set_ylabel("Barn temperature (°C)")
    ax.set_title("Daily barn temperature during summer months (mean ± daily range)")
    ax.legend(fontsize=9)
    save_fig(fig, "climate_temp_timeseries", out_dir)

    climate.to_csv(out_dir / "climate_daily.csv", index=False)
    summary.to_csv(out_dir / "climate_summary.csv", index=False)
    return summary


# ─────────────────────────────────────────────────────────────
#  « breakpoint stability across years »
# ─────────────────────────────────────────────────────────────

def compute_breakpoint_stability(results: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    """Assess within-animal breakpoint stability for repeat animals.

    For animals appearing in multiple years, computes ICC (intraclass
    correlation) and paired statistics.

    Args:
        results: Full results DataFrame (all years).
        out_dir: Output directory.

    Returns:
        DataFrame of repeat-animal breakpoint pairs.
    """
    import matplotlib.pyplot as plt
    setup_plotting()

    converged = results[results["thi_converged"] == True].copy()
    counts = converged.groupby("animal_id").size()
    repeats = counts[counts >= 2].index

    if len(repeats) < 3:
        log.info("Only %d repeat animals — skipping stability analysis", len(repeats))
        return pd.DataFrame()

    log.info("Breakpoint stability: %d animals with 2+ years of data", len(repeats))

    repeat_data = converged[converged["animal_id"].isin(repeats)].copy()
    repeat_data = repeat_data.sort_values(["animal_id", "year"])

    # ── ICC for THI breakpoints ──────────────────────────────
    # Simple one-way random ICC
    groups = [
        grp["thi_breakpoint"].values
        for _, grp in repeat_data.groupby("animal_id")
    ]
    # Between-animal variance vs within-animal variance
    all_vals = repeat_data["thi_breakpoint"].dropna()
    grand_mean = all_vals.mean()
    n_animals = len(groups)

    ss_between = sum(len(g) * (g.mean() - grand_mean) ** 2 for g in groups)
    ss_within = sum(np.sum((g - g.mean()) ** 2) for g in groups)
    k_mean = len(all_vals) / n_animals  # mean group size

    ms_between = ss_between / max(n_animals - 1, 1)
    ms_within = ss_within / max(len(all_vals) - n_animals, 1)
    icc = (ms_between - ms_within) / (ms_between + (k_mean - 1) * ms_within) if ms_between + (k_mean - 1) * ms_within > 0 else np.nan

    log.info("  THI breakpoint ICC = %.3f (n=%d animals, %d observations)",
             icc, n_animals, len(all_vals))

    # ── paired year-to-year plot ─────────────────────────────
    # For animals with exactly 2 years: scatter year1 vs year2
    pairs = []
    for aid, grp in repeat_data.groupby("animal_id"):
        if len(grp) >= 2:
            sorted_grp = grp.sort_values("year")
            for i in range(len(sorted_grp) - 1):
                pairs.append({
                    "animal_id": aid,
                    "year_1": int(sorted_grp.iloc[i]["year"]),
                    "year_2": int(sorted_grp.iloc[i + 1]["year"]),
                    "thi_bp_1": sorted_grp.iloc[i]["thi_breakpoint"],
                    "thi_bp_2": sorted_grp.iloc[i + 1]["thi_breakpoint"],
                    "temp_bp_1": sorted_grp.iloc[i]["temp_breakpoint"],
                    "temp_bp_2": sorted_grp.iloc[i + 1]["temp_breakpoint"],
                })

    pairs_df = pd.DataFrame(pairs)
    if not pairs_df.empty:
        pairs_df.to_csv(out_dir / "breakpoint_stability_pairs.csv", index=False)

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(pairs_df["thi_bp_1"], pairs_df["thi_bp_2"],
                   s=30, alpha=0.6, color=COLOURS["scatter"])
        lims = [
            min(pairs_df["thi_bp_1"].min(), pairs_df["thi_bp_2"].min()) - 2,
            max(pairs_df["thi_bp_1"].max(), pairs_df["thi_bp_2"].max()) + 2,
        ]
        ax.plot(lims, lims, "--", color=COLOURS["identity"], linewidth=1, label="Identity line")
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel("THI breakpoint (earlier year)")
        ax.set_ylabel("THI breakpoint (later year)")
        ax.set_title(
            f"Within-animal breakpoint stability across years\n"
            f"ICC = {icc:.3f}, n = {len(pairs_df)} pairs from {len(repeats)} animals",
        )
        ax.legend()
        save_fig(fig, "breakpoint_stability_scatter", out_dir)

    return pairs_df


# ─────────────────────────────────────────────────────────────
#  « breakpoint predictors: lactation number, milk yield »
# ─────────────────────────────────────────────────────────────

SQL_MILK_YIELD = """
SELECT
    animal_id,
    AVG(herdeplus_milked_mkg) AS mean_milk_yield_kg,
    COUNT(*) AS n_milkings
FROM herdeplus
WHERE animal_id = ?
  AND "timestamp" >= ?
  AND "timestamp" <= ?
  AND herdeplus_milked_mkg IS NOT NULL
  AND herdeplus_milked_mkg > 0
GROUP BY animal_id
"""

SQL_LACTATION_NR = """
SELECT MAX(CAST(herdeplus_calving_lactation AS INTEGER)) AS lactation_nr
FROM herdeplus
WHERE animal_id = ?
  AND herdeplus_calving_lactation IS NOT NULL
  AND herdeplus_calving_lactation > 0
"""


def add_production_data(
    con, results: pd.DataFrame, tierauswahl: pd.DataFrame,
) -> pd.DataFrame:
    """Join breakpoint results with milk yield and lactation number.

    Milk yield is averaged over the observation window.  Lactation number
    is taken as the maximum recorded value across the animal's full
    history (since it is only populated on calving events, not every
    milking record).

    Args:
        con: Database connection.
        results: Breakpoint results DataFrame.
        tierauswahl: Animal selection with date windows.

    Returns:
        Results enriched with ``mean_milk_yield_kg`` and ``lactation_nr``.
    """
    records = []
    for _, row in tierauswahl.iterrows():
        aid = int(row["animal_id"])
        enter = str(row["datetime_enter"])[:10]
        exit_ = str(row["datetime_exit"])[:10]

        rec = {"animal_id": aid, "date_enter": enter}

        # Milk yield from observation window
        milk = query_df(con, SQL_MILK_YIELD, (aid, enter, exit_))
        if not milk.empty:
            rec["mean_milk_yield_kg"] = milk.iloc[0]["mean_milk_yield_kg"]
        else:
            rec["mean_milk_yield_kg"] = np.nan

        # Lactation number from full animal history
        lac = query_df(con, SQL_LACTATION_NR, (aid,))
        if not lac.empty and pd.notna(lac.iloc[0]["lactation_nr"]):
            rec["lactation_nr"] = int(lac.iloc[0]["lactation_nr"])
        else:
            rec["lactation_nr"] = np.nan

        records.append(rec)

    if not records:
        return results

    prod = pd.DataFrame(records)
    merged = results.merge(prod, on=["animal_id", "date_enter"], how="left")
    n_lac = merged["lactation_nr"].notna().sum()
    n_milk = merged["mean_milk_yield_kg"].notna().sum()
    log.info("  Production data joined: %d with milk yield, %d with lactation nr", n_milk, n_lac)
    return merged


def analyse_breakpoint_predictors(results: pd.DataFrame, out_dir: Path) -> None:
    """Test whether lactation number or milk yield predict breakpoints.

    Fits OLS: breakpoint ~ lactation_nr + mean_milk_yield_kg + year.
    Generates scatter plots with regression lines.

    Args:
        results: Results with production data merged.
        out_dir: Output directory.
    """
    import matplotlib.pyplot as plt
    from scipy.stats import pearsonr
    setup_plotting()

    converged = results[results["thi_converged"] == True].copy()

    for bp_col, bp_label in [
        ("thi_breakpoint", "THI breakpoint"),
        ("temp_breakpoint", "Barn temperature breakpoint (°C)"),
    ]:
        conv = converged.dropna(subset=[bp_col])

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # ── milk yield vs breakpoint (breakpoint on x) ──────
        ax = axes[0]
        sub = conv.dropna(subset=["mean_milk_yield_kg", bp_col])
        if len(sub) > 10:
            ax.scatter(sub[bp_col], sub["mean_milk_yield_kg"],
                       s=20, alpha=0.5, color=COLOURS["scatter"])
            r, p = pearsonr(sub[bp_col], sub["mean_milk_yield_kg"])
            z = np.polyfit(sub[bp_col], sub["mean_milk_yield_kg"], 1)
            x_line = np.linspace(sub[bp_col].min(), sub[bp_col].max(), 50)
            ax.plot(x_line, np.polyval(z, x_line), "--", color=COLOURS["fit_line"])
            ax.set_title(f"r = {r:.3f}, p = {p:.3f}, n = {len(sub)}", fontsize=10)
        ax.set_xlabel(bp_label)
        ax.set_ylabel("Mean milk yield (kg/milking)")

        # ── lactation number vs breakpoint ───────────────────
        ax = axes[1]
        sub = conv.dropna(subset=["lactation_nr", bp_col])
        sub = sub[sub["lactation_nr"] > 0]
        if len(sub) > 10:
            jitter = np.random.uniform(-0.2, 0.2, len(sub))
            ax.scatter(sub[bp_col], sub["lactation_nr"] + jitter,
                       s=20, alpha=0.5, color=COLOURS["scatter_alt"])
            r, p = pearsonr(sub[bp_col], sub["lactation_nr"])
            z = np.polyfit(sub[bp_col], sub["lactation_nr"], 1)
            x_line = np.linspace(sub[bp_col].min(), sub[bp_col].max(), 50)
            ax.plot(x_line, np.polyval(z, x_line), "--", color=COLOURS["fit_line"])
            ax.set_title(f"r = {r:.3f}, p = {p:.3f}, n = {len(sub)}", fontsize=10)
        else:
            ax.text(0.5, 0.5, f"Insufficient data\n(n={len(sub)} with lactation nr)",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=11, color=COLOURS["identity"])
        ax.set_xlabel(bp_label)
        ax.set_ylabel("Lactation number")

        name = "thi" if "thi" in bp_col else "temp"
        fig.suptitle(f"{bp_label}: association with production parameters")
        fig.tight_layout()
        save_fig(fig, f"predictors_{name}_breakpoint", out_dir)


# ─────────────────────────────────────────────────────────────
#  « behavioural response curves »
# ─────────────────────────────────────────────────────────────

SQL_BEHAVIOUR = """
SELECT
    "timestamp",
    CAST("act_index" AS REAL) AS act_index,
    CAST("rum_index" AS REAL) AS rum_index,
    CAST("temp_without_drink_cycles" AS REAL) AS body_temp
FROM smaxtec_derived
WHERE animal_id = ?
  AND "timestamp" >= ?
  AND "timestamp" <= ?
  AND "temp_without_drink_cycles" IS NOT NULL
  AND CAST("temp_without_drink_cycles" AS REAL) > 30
  AND CAST("temp_without_drink_cycles" AS REAL) < 43
"""


def compute_behavioural_response(
    con, results: pd.DataFrame, tierauswahl: pd.DataFrame,
    out_dir: Path,
) -> None:
    """Analyse activity and rumination relative to individual breakpoints.

    For each animal with a converged THI breakpoint, splits observations
    into below/above breakpoint and compares mean activity and rumination.

    Args:
        con: Database connection.
        results: Breakpoint results.
        tierauswahl: Animal selection.
        out_dir: Output directory.
    """
    import matplotlib.pyplot as plt
    setup_plotting()

    converged = results[results["thi_converged"] == True].copy()
    if converged.empty:
        return

    barn_cache: dict = {}
    records = []

    for _, row in converged.iterrows():
        aid = int(row["animal_id"])
        enter = row["date_enter"]
        exit_ = row["date_exit"]
        bp = row["thi_breakpoint"]
        year = int(row["year"])

        # Get behaviour data
        beh = query_df(con, SQL_BEHAVIOUR, (aid, enter, exit_))
        if len(beh) < 50:
            continue

        beh["timestamp"] = pd.to_datetime(beh["timestamp"])

        # Get barn THI
        barn_df = load_animal_data(con, aid, enter, exit_, barn_cache=barn_cache)
        if barn_df.empty:
            continue

        # Merge behaviour with barn THI
        beh["hour_key"] = beh["timestamp"].dt.floor("h")
        barn_hourly = barn_df[["hour_key", "barn_thi"]].drop_duplicates("hour_key")
        beh = beh.merge(barn_hourly, on="hour_key", how="inner")

        below = beh[beh["barn_thi"] <= bp]
        above = beh[beh["barn_thi"] > bp]

        if len(below) < 10 or len(above) < 10:
            continue

        records.append({
            "animal_id": aid,
            "year": year,
            "thi_breakpoint": bp,
            "act_below": below["act_index"].mean(),
            "act_above": above["act_index"].mean(),
            "rum_below": below["rum_index"].mean(),
            "rum_above": above["rum_index"].mean(),
            "body_temp_below": below["body_temp"].mean(),
            "body_temp_above": above["body_temp"].mean(),
            "n_below": len(below),
            "n_above": len(above),
        })

    if not records:
        log.warning("No behavioural data collected")
        return

    beh_df = pd.DataFrame(records)
    beh_df.to_csv(out_dir / "behavioural_response.csv", index=False)
    log.info("Behavioural response: %d animals", len(beh_df))

    # ── paired comparison plots ──────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    for ax, below_col, above_col, ylabel, title in [
        (axes[0], "body_temp_below", "body_temp_above",
         "Mean rumen temperature (°C)", "Rumen temperature"),
        (axes[1], "act_below", "act_above",
         "Mean activity index", "Activity"),
        (axes[2], "rum_below", "rum_above",
         "Mean rumination index", "Rumination"),
    ]:
        below_vals = beh_df[below_col].dropna()
        above_vals = beh_df[above_col].dropna()

        positions = [0, 1]
        bp_plot = ax.boxplot(
            [below_vals, above_vals],
            positions=positions,
            labels=["Below\nbreakpoint", "Above\nbreakpoint"],
            patch_artist=True,
            widths=0.5,
        )
        bp_plot["boxes"][0].set_facecolor(COLOURS["below_bp"])
        bp_plot["boxes"][0].set_alpha(0.5)
        bp_plot["boxes"][1].set_facecolor(COLOURS["above_bp"])
        bp_plot["boxes"][1].set_alpha(0.5)
        for median in bp_plot["medians"]:
            median.set_color("#333333")
            median.set_linewidth(2)

        # Paired lines
        for i in range(len(beh_df)):
            b = beh_df.iloc[i][below_col]
            a = beh_df.iloc[i][above_col]
            if pd.notna(b) and pd.notna(a):
                ax.plot([0, 1], [b, a], "-", color=COLOURS["paired_line"],
                        alpha=0.2, linewidth=0.5)

        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=11)

        # Wilcoxon test
        from scipy.stats import wilcoxon
        paired = beh_df[[below_col, above_col]].dropna()
        if len(paired) >= 10:
            try:
                stat, p = wilcoxon(paired[below_col], paired[above_col])
                ax.text(0.5, 0.02, f"Wilcoxon p = {p:.4f}\nn = {len(paired)}",
                        transform=ax.transAxes, ha="center", fontsize=9,
                        color="#555555")
            except Exception:
                pass

    fig.suptitle(
        "Behavioural and physiological parameters below vs. above\n"
        "individual THI breakpoint (per-animal means)",
    )
    fig.tight_layout()
    save_fig(fig, "behavioural_below_vs_above", out_dir)


# ─────────────────────────────────────────────────────────────
#  « summary table for manuscript »
# ─────────────────────────────────────────────────────────────

def generate_summary_table(results: pd.DataFrame, out_dir: Path) -> None:
    """Generate Table 1 for the manuscript: summary statistics per year.

    Args:
        results: Full results DataFrame.
        out_dir: Output directory.
    """
    years = sorted(results["year"].dropna().unique().astype(int))
    rows = []

    for year in years:
        yr = results[results["year"] == year]
        thi_conv = yr[yr["thi_converged"] == True]
        temp_conv = yr[yr["temp_converged"] == True]

        rows.append({
            "Year": year,
            "n_animals": len(yr),
            "n_THI_converged": len(thi_conv),
            "n_temp_converged": len(temp_conv),
            "THI_bp_median": thi_conv["thi_breakpoint"].median() if len(thi_conv) > 0 else np.nan,
            "THI_bp_IQR": (
                f"[{thi_conv['thi_breakpoint'].quantile(0.25):.1f}–"
                f"{thi_conv['thi_breakpoint'].quantile(0.75):.1f}]"
            ) if len(thi_conv) > 0 else "",
            "THI_bp_range": (
                f"[{thi_conv['thi_breakpoint'].min():.1f}–"
                f"{thi_conv['thi_breakpoint'].max():.1f}]"
            ) if len(thi_conv) > 0 else "",
            "temp_bp_median": temp_conv["temp_breakpoint"].median() if len(temp_conv) > 0 else np.nan,
            "temp_bp_IQR": (
                f"[{temp_conv['temp_breakpoint'].quantile(0.25):.1f}–"
                f"{temp_conv['temp_breakpoint'].quantile(0.75):.1f}]"
            ) if len(temp_conv) > 0 else "",
            "mean_n_readings": yr["n_readings"].mean(),
        })

    # Add overall row
    thi_all = results[results["thi_converged"] == True]
    temp_all = results[results["temp_converged"] == True]
    rows.append({
        "Year": "All",
        "n_animals": len(results),
        "n_THI_converged": len(thi_all),
        "n_temp_converged": len(temp_all),
        "THI_bp_median": thi_all["thi_breakpoint"].median() if len(thi_all) > 0 else np.nan,
        "THI_bp_IQR": (
            f"[{thi_all['thi_breakpoint'].quantile(0.25):.1f}–"
            f"{thi_all['thi_breakpoint'].quantile(0.75):.1f}]"
        ) if len(thi_all) > 0 else "",
        "THI_bp_range": (
            f"[{thi_all['thi_breakpoint'].min():.1f}–"
            f"{thi_all['thi_breakpoint'].max():.1f}]"
        ) if len(thi_all) > 0 else "",
        "temp_bp_median": temp_all["temp_breakpoint"].median() if len(temp_all) > 0 else np.nan,
        "temp_bp_IQR": (
            f"[{temp_all['temp_breakpoint'].quantile(0.25):.1f}–"
            f"{temp_all['temp_breakpoint'].quantile(0.75):.1f}]"
        ) if len(temp_all) > 0 else "",
        "mean_n_readings": results["n_readings"].mean(),
    })

    table = pd.DataFrame(rows)
    table.to_csv(out_dir / "summary_table.csv", index=False)
    log.info("Summary table saved.")

    # Log it nicely
    log.info("─" * 60)
    log.info("TABLE 1: Summary of broken-stick regression results")
    log.info("─" * 60)
    for _, r in table.iterrows():
        log.info(
            "  %s: n=%s, THI converged=%s, THI median=%.1f %s, "
            "Temp converged=%s, Temp median=%.1f %s",
            r["Year"], r["n_animals"],
            r["n_THI_converged"],
            r["THI_bp_median"] if pd.notna(r["THI_bp_median"]) else 0,
            r["THI_bp_IQR"],
            r["n_temp_converged"],
            r["temp_bp_median"] if pd.notna(r["temp_bp_median"]) else 0,
            r["temp_bp_IQR"],
        )



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

    # ── 1. Broken-stick regression (core analysis) ───────────
    log.info("=" * 60)
    log.info("1. Running broken-stick regression for %d entries …",
             len(tierauswahl))
    results = run_broken_stick_analysis(con, tierauswahl)
    results.to_csv(args.out / "broken_stick_results.csv", index=False)

    converged_thi = results[results["thi_converged"] == True]
    converged_temp = results[results["temp_converged"] == True]
    log.info(
        "  Converged: %d/%d THI fits, %d/%d barn temp fits",
        len(converged_thi), len(results),
        len(converged_temp), len(results),
    )

    # ── 2. Spearman correlations ─────────────────────────────
    log.info("=" * 60)
    log.info("2. Computing Spearman correlations …")
    spearman = compute_spearman_correlations(con, tierauswahl)
    if not spearman.empty:
        spearman.to_csv(args.out / "spearman_correlations.csv", index=False)
        thi_rs = spearman["thi_spearman_rs"].dropna()
        temp_rs = spearman["temp_spearman_rs"].dropna()
        log.info(
            "  THI rₛ: median=%.3f, range=[%.3f, %.3f], n=%d",
            thi_rs.median(), thi_rs.min(), thi_rs.max(), len(thi_rs),
        )
        log.info(
            "  Temp rₛ: median=%.3f, range=[%.3f, %.3f], n=%d",
            temp_rs.median(), temp_rs.min(), temp_rs.max(), len(temp_rs),
        )
        plot_spearman_results(spearman, args.out)

    # ── 3. Descriptive climate summary ───────────────────────
    log.info("=" * 60)
    log.info("3. Computing barn climate summary …")
    compute_climate_summary(con, tierauswahl, args.out)

    # ── 4. Breakpoint predictors (milk yield, lactation) ─────
    log.info("=" * 60)
    log.info("4. Joining production data and testing breakpoint predictors …")
    results = add_production_data(con, results, tierauswahl)
    results.to_csv(args.out / "broken_stick_results.csv", index=False)
    analyse_breakpoint_predictors(results, args.out)

    # ── 5. Breakpoint stability across years ─────────────────
    log.info("=" * 60)
    log.info("5. Assessing breakpoint stability for repeat animals …")
    compute_breakpoint_stability(results, args.out)

    # ── 6. Behavioural response below/above breakpoint ───────
    log.info("=" * 60)
    log.info("6. Computing behavioural response curves …")
    compute_behavioural_response(con, results, tierauswahl, args.out)

    # ── 7. Summary table and plots ───────────────────────────
    log.info("=" * 60)
    log.info("7. Generating summary table and breakpoint plots …")
    generate_summary_table(results, args.out)
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
    log.info("=" * 60)
    log.info("All analyses complete. Outputs in: %s", args.out)
    log.info("=" * 60)


if __name__ == "__main__":
    main()

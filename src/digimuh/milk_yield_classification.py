#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — milk_yield_classification                            ║
# ║  « daily MY distribution vs literature tercile boundaries »    ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Produces a pooled and a per-year histogram of daily milk       ║
# ║  yield for the extended 255-animal-year Tierauswahl, overlays   ║
# ║  the Müschner-Siemens (2020) and Yan (2021) low/middle/high     ║
# ║  boundaries plus our own 33rd / 67th percentile terciles, and   ║
# ║  prints the class counts each scheme would produce.  Used to    ║
# ║  pick a defensible MY stratification for the Frontiers paper.   ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Milk yield classification analysis for the Frontiers manuscript.

Reads ``daily_milk_yield.csv`` from the extract stage and compares
its pooled distribution against two published tercile schemes and
the herd's own 33rd / 67th percentiles.  Writes two figures
(pooled and per-year) plus a console summary.

Usage::

    digimuh-milk-yield --data results/broken_stick

References:
    Müschner-Siemens et al. (2020) *J Thermal Biol* 88:102484.
    Yan et al. (2021) *J Thermal Biol* 100:103041.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.constants import (
    WONG_BLUE, WONG_GREEN, WONG_ORANGE, WONG_VERMILLION, WONG_GREY,
    WONG_SKY, WONG_PINK,
    COLOURS,
)
from digimuh.viz_base import setup_figure, save_figure
from digimuh.paths import resolve_input, resolve_output

log = logging.getLogger("digimuh.milk_yield")


# ─────────────────────────────────────────────────────────────
#  « classification schemes »
# ─────────────────────────────────────────────────────────────
#  Each scheme: (low_upper_bound, middle_upper_bound, linestyle, colour)
#  A cow-day with y ≤ low_upper goes to "low"; ≤ middle_upper → "middle";
#  otherwise "high".

MUESCHNER_SIEMENS: tuple[float, float] = (28.8, 38.4)
"""Müschner-Siemens et al. (2020): terciles of their German Holstein
cohort at ATB Potsdam.  Same region / breed as our herd."""

YAN: tuple[float, float] = (26.0, 39.0)
"""Yan et al. (2021): 826-cow Chinese Holstein cohort."""


def classify(values: pd.Series, low_bd: float, middle_bd: float) -> pd.Series:
    """Bucket a yield series into low / middle / high.

    Args:
        values: Daily milk yield (kg/d) as a pandas Series.
        low_bd: Upper bound of the "low" bucket (inclusive).
        middle_bd: Upper bound of the "middle" bucket (inclusive).

    Returns:
        Categorical series of labels ``"low"``, ``"middle"``, ``"high"``.
    """
    out = pd.Series(
        np.where(values <= low_bd, "low",
                 np.where(values <= middle_bd, "middle", "high")),
        index=values.index, name="class",
    )
    return pd.Categorical(out, categories=["low", "middle", "high"], ordered=True)


# ─────────────────────────────────────────────────────────────
#  « data loading »
# ─────────────────────────────────────────────────────────────

def load_daily_yields(data_dir: Path) -> pd.DataFrame:
    """Load and clean ``daily_milk_yield.csv``.

    Drops rows with non-positive or NaN ``daily_yield_kg``.

    Args:
        data_dir: Directory containing the extract-stage CSVs.

    Returns:
        DataFrame with at least ``animal_id``, ``year``, ``date``,
        ``daily_yield_kg``.
    """
    path = resolve_input(data_dir, "daily_milk_yield.csv")
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found — run `digimuh-extract` first.")
    df = pd.read_csv(path)
    df = df[df["daily_yield_kg"].notna() & (df["daily_yield_kg"] > 0)].copy()
    df["year"] = df["year"].astype(int)
    df["date"] = pd.to_datetime(df["date"]).dt.date
    return df


# ─────────────────────────────────────────────────────────────
#  « summary table »
# ─────────────────────────────────────────────────────────────

def print_summary(df: pd.DataFrame, our_terciles: tuple[float, float]) -> None:
    """Print per-year, pooled, and class-count tables to the console.

    Uses the rich console helpers if available, otherwise plain text.
    """
    from digimuh.console import section, result_table, kv, banner

    banner("Milk yield classification analysis")

    # ── Per-year stats ────────────────────────────────────────
    section("Per-year daily milk yield (kg/d)",
            "n_days, n_animals, mean, median, Q33, Q67, min, max")
    year_rows = []
    for year in sorted(df["year"].unique()):
        sub = df[df["year"] == year]["daily_yield_kg"]
        year_rows.append([
            int(year),
            len(sub),
            sub.index.to_series().count() if False else
                df[df["year"] == year]["animal_id"].nunique(),
            f"{sub.mean():.2f}",
            f"{sub.median():.2f}",
            f"{sub.quantile(0.33):.2f}",
            f"{sub.quantile(0.67):.2f}",
            f"{sub.min():.2f}",
            f"{sub.max():.2f}",
        ])
    pooled = df["daily_yield_kg"]
    year_rows.append([
        "pooled",
        len(pooled),
        df["animal_id"].nunique(),
        f"{pooled.mean():.2f}",
        f"{pooled.median():.2f}",
        f"{pooled.quantile(0.33):.2f}",
        f"{pooled.quantile(0.67):.2f}",
        f"{pooled.min():.2f}",
        f"{pooled.max():.2f}",
    ])
    result_table(
        "Yield distribution",
        ["Year", "n_days", "n_animals", "Mean", "Median",
         "Q33", "Q67", "Min", "Max"],
        year_rows,
    )

    # ── Class counts per scheme ──────────────────────────────
    section("Class counts by scheme",
            "cow-days in low / middle / high buckets")
    schemes = [
        ("Müschner-Siemens 2020", MUESCHNER_SIEMENS),
        ("Yan 2021", YAN),
        (f"Our terciles ({our_terciles[0]:.1f} / {our_terciles[1]:.1f})",
         our_terciles),
    ]
    class_rows = []
    for name, (low_bd, mid_bd) in schemes:
        classes = classify(df["daily_yield_kg"], low_bd, mid_bd)
        counts = pd.Series(classes).value_counts()
        n_total = len(classes)
        class_rows.append([
            name,
            f"{low_bd:.1f}",
            f"{mid_bd:.1f}",
            f"{counts.get('low', 0):,} ({100*counts.get('low', 0)/n_total:.1f}%)",
            f"{counts.get('middle', 0):,} ({100*counts.get('middle', 0)/n_total:.1f}%)",
            f"{counts.get('high', 0):,} ({100*counts.get('high', 0)/n_total:.1f}%)",
        ])
    result_table(
        "Scheme comparison",
        ["Scheme", "Low ≤", "Middle ≤", "Low", "Middle", "High"],
        class_rows,
    )

    # ── Unique animals per class (for Frontiers sample-size check) ──
    section("Unique animals per class",
            "n_animals with ≥1 cow-day in each bucket — "
            "informs whether broken-stick by group is feasible")
    animal_rows = []
    for name, (low_bd, mid_bd) in schemes:
        dd = df.copy()
        dd["class"] = classify(dd["daily_yield_kg"], low_bd, mid_bd)
        n_low    = dd.loc[dd["class"] == "low",    "animal_id"].nunique()
        n_middle = dd.loc[dd["class"] == "middle", "animal_id"].nunique()
        n_high   = dd.loc[dd["class"] == "high",   "animal_id"].nunique()
        animal_rows.append([name, n_low, n_middle, n_high])
    result_table(
        "Animals per class",
        ["Scheme", "Low", "Middle", "High"],
        animal_rows,
    )


# ─────────────────────────────────────────────────────────────
#  « plotting »
# ─────────────────────────────────────────────────────────────

_SCHEME_COLOURS = {
    "muesch":    WONG_ORANGE,
    "yan":       WONG_GREEN,
    "terciles":  WONG_VERMILLION,
}
_SCHEME_LINESTYLES = {
    "muesch":    "--",
    "yan":       ":",
    "terciles":  "-",
}


def _add_boundary_lines(ax, our_terciles: tuple[float, float],
                        show_legend: bool = True) -> None:
    """Overlay the three scheme boundaries on an axes."""
    specs = [
        ("muesch",   MUESCHNER_SIEMENS,
         f"Müschner-Siemens 2020 ({MUESCHNER_SIEMENS[0]} / {MUESCHNER_SIEMENS[1]})"),
        ("yan",      YAN,
         f"Yan 2021 ({YAN[0]} / {YAN[1]})"),
        ("terciles", our_terciles,
         f"Our terciles ({our_terciles[0]:.1f} / {our_terciles[1]:.1f})"),
    ]
    for key, (low_bd, mid_bd), label in specs:
        colour = _SCHEME_COLOURS[key]
        ls = _SCHEME_LINESTYLES[key]
        ax.axvline(low_bd, color=colour, linestyle=ls, linewidth=1.8, alpha=0.9,
                   label=label)
        ax.axvline(mid_bd, color=colour, linestyle=ls, linewidth=1.8, alpha=0.9)


def plot_pooled_histogram(
    df: pd.DataFrame, our_terciles: tuple[float, float], out_dir: Path,
) -> None:
    """Single-panel histogram of all animal-days with all three schemes."""
    import matplotlib.pyplot as plt
    setup_figure()

    values = df["daily_yield_kg"].values
    log.info("  Plotting pooled histogram (%d cow-days) …", len(values))

    fig, ax = plt.subplots(figsize=(11, 6))

    max_y = float(np.ceil(values.max()))
    bins = np.arange(0, max_y + 1, 1.0)
    ax.hist(values, bins=bins, color=WONG_BLUE, edgecolor="white",
            alpha=0.78, linewidth=0.4)

    _add_boundary_lines(ax, our_terciles)

    n = len(values)
    mean_ = values.mean()
    median_ = np.median(values)
    anno = (f"n = {n:,} cow-days\n"
            f"mean = {mean_:.2f} kg/d\n"
            f"median = {median_:.2f} kg/d\n"
            f"animals = {df['animal_id'].nunique()}")
    ax.text(0.99, 0.97, anno, transform=ax.transAxes, va="top", ha="right",
            fontsize=10, color="#333",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      alpha=0.85, edgecolor=WONG_GREY))

    ax.set_xlabel("Daily milk yield (kg/d)")
    ax.set_ylabel("Cow-days")
    ax.set_title("Daily milk yield distribution — all animal-years pooled\n"
                 "(2021–2024, extended Tierauswahl)")
    ax.set_xlim(0, max_y + 1)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(False)
    fig.tight_layout()
    save_figure(fig, "milk_yield_histogram_pooled", out_dir)


def plot_per_year_histogram(
    df: pd.DataFrame, our_terciles: tuple[float, float], out_dir: Path,
) -> None:
    """Per-year histograms stacked vertically with a shared x-axis."""
    import matplotlib.pyplot as plt
    setup_figure()

    years = sorted(df["year"].unique())
    log.info("  Plotting per-year histogram (%d years) …", len(years))

    max_y = float(np.ceil(df["daily_yield_kg"].max()))
    bins = np.arange(0, max_y + 1, 1.0)

    fig, axes = plt.subplots(len(years), 1,
                             figsize=(11, 2.6 * len(years)),
                             sharex=True)
    if len(years) == 1:
        axes = [axes]

    for ax, year in zip(axes, years):
        sub = df[df["year"] == year]["daily_yield_kg"].values
        ax.hist(sub, bins=bins, color=WONG_BLUE, edgecolor="white",
                alpha=0.78, linewidth=0.4)
        _add_boundary_lines(ax, our_terciles, show_legend=(year == years[0]))

        n = len(sub)
        n_animals = df[df["year"] == year]["animal_id"].nunique()
        anno = (f"{year}\n"
                f"n = {n:,} cow-days\n"
                f"animals = {n_animals}\n"
                f"median = {np.median(sub):.1f} kg/d")
        ax.text(0.99, 0.95, anno, transform=ax.transAxes, va="top", ha="right",
                fontsize=9, color="#333",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          alpha=0.85, edgecolor=WONG_GREY))

        ax.set_ylabel("Cow-days")
        ax.grid(False)

    axes[0].legend(fontsize=8, loc="upper left")
    axes[-1].set_xlabel("Daily milk yield (kg/d)")
    axes[-1].set_xlim(0, max_y + 1)
    fig.suptitle("Daily milk yield distribution by year",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "milk_yield_histogram_by_year", out_dir)


# ─────────────────────────────────────────────────────────────
#  « Wood-residual (DIM-adjusted) mode »
# ─────────────────────────────────────────────────────────────

def print_wood_summary(wood: pd.DataFrame, fits: pd.DataFrame) -> None:
    """Pretty-print Wood-curve coverage, fit stats, and residual terciles."""
    from digimuh.console import section, result_table, kv, stars_styled

    section("Wood lactation-curve coverage",
            "per-lactation fits with ≥30 points use per_lactation, "
            "else fall back to the parity pool")
    n_total   = len(fits)
    n_indiv   = (fits["method"] == "per_lactation").sum()
    n_pool    = (fits["method"] == "parity_fallback").sum()
    n_animals = fits["animal_id"].nunique()
    kv("Lactations fitted", f"{n_total} from {n_animals} animals")
    kv("  individual fits", f"{n_indiv} ({100 * n_indiv / n_total:.1f}%)")
    kv("  parity fallback", f"{n_pool} ({100 * n_pool / n_total:.1f}%)")

    indiv = fits[fits["method"] == "per_lactation"]
    if not indiv.empty:
        kv("Per-lactation R² median [IQR]",
           f"{indiv['r_squared'].median():.3f} "
           f"[{indiv['r_squared'].quantile(0.25):.3f}"
           f"–{indiv['r_squared'].quantile(0.75):.3f}]")
        kv("Peak DIM median",
           f"{indiv['peak_dim'].median():.0f} d "
           f"(expected ~40–70 d for Holsteins)")
        kv("Peak yield median",
           f"{indiv['peak_yield'].median():.1f} kg/d")

    section("Residual distribution",
            "yield − Wood(DIM) — should be approximately symmetric")
    r = wood["yield_residual"].dropna()
    from scipy import stats as _st
    kv("n",       f"{len(r):,}")
    kv("mean",    f"{r.mean():.3f} kg/d")
    kv("median",  f"{r.median():.3f} kg/d")
    kv("SD",      f"{r.std():.3f} kg/d")
    kv("skewness",        f"{_st.skew(r):+.3f} (raw yield: +0.41)")
    kv("excess kurtosis", f"{_st.kurtosis(r, fisher=True):+.3f}")

    q33, q67 = float(r.quantile(0.33)), float(r.quantile(0.67))
    kv("Residual Q33",    f"{q33:+.2f} kg/d")
    kv("Residual Q67",    f"{q67:+.2f} kg/d")

    section("Class counts under Wood-residual terciles",
            "primary Frontiers stratifier (DIM & parity removed)")
    classes = classify(r, q33, q67)
    counts = pd.Series(classes).value_counts()
    n = len(r)
    class_rows = [[
        "Wood-residual terciles",
        f"{q33:+.2f}", f"{q67:+.2f}",
        f"{counts.get('low', 0):,} ({100 * counts.get('low', 0) / n:.1f}%)",
        f"{counts.get('middle', 0):,} ({100 * counts.get('middle', 0) / n:.1f}%)",
        f"{counts.get('high', 0):,} ({100 * counts.get('high', 0) / n:.1f}%)",
    ]]
    result_table(
        "Residual scheme",
        ["Scheme", "Low ≤", "Middle ≤", "Low", "Middle", "High"],
        class_rows,
    )

    # Animals per class (sample-size sanity check)
    wdf = wood.dropna(subset=["yield_residual"]).copy()
    wdf["class"] = classify(wdf["yield_residual"], q33, q67)
    section("Unique animals per residual class",
            "is there enough per bucket for stratified broken-stick?")
    result_table(
        "Animals per class",
        ["Scheme", "Low", "Middle", "High"],
        [[
            "Wood-residual terciles",
            wdf.loc[wdf["class"] == "low",    "animal_id"].nunique(),
            wdf.loc[wdf["class"] == "middle", "animal_id"].nunique(),
            wdf.loc[wdf["class"] == "high",   "animal_id"].nunique(),
        ]],
    )


def plot_residual_histogram(
    wood: pd.DataFrame, out_dir: Path,
) -> tuple[float, float]:
    """Histogram of Wood residuals with tercile and zero lines.

    Returns the (Q33, Q67) residual tercile boundaries so the caller
    can persist them or compare against the raw-yield terciles.
    """
    import matplotlib.pyplot as plt
    setup_figure()

    r = wood["yield_residual"].dropna().values
    q33 = float(np.quantile(r, 0.33))
    q67 = float(np.quantile(r, 0.67))
    log.info("  Plotting residual histogram (%d cow-days) …", len(r))

    lo = float(np.floor(np.quantile(r, 0.005)))
    hi = float(np.ceil(np.quantile(r, 0.995)))
    bins = np.arange(lo, hi + 1, 0.5)

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.hist(r, bins=bins, color=WONG_SKY, edgecolor="white",
            alpha=0.78, linewidth=0.4)

    ax.axvline(0, color=WONG_GREY, linestyle="--", linewidth=1.3,
               label="Expected (Wood curve)")
    ax.axvline(q33, color=WONG_VERMILLION, linewidth=1.8,
               label=f"Residual tercile Q33 = {q33:+.2f} kg/d")
    ax.axvline(q67, color=WONG_VERMILLION, linewidth=1.8)

    anno = (f"n = {len(r):,} cow-days\n"
            f"mean = {r.mean():+.2f} kg/d\n"
            f"median = {np.median(r):+.2f} kg/d\n"
            f"SD = {r.std():.2f} kg/d")
    ax.text(0.99, 0.97, anno, transform=ax.transAxes, va="top", ha="right",
            fontsize=10, color="#333",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      alpha=0.85, edgecolor=WONG_GREY))

    ax.set_xlabel("DIM-adjusted daily milk yield residual (kg/d)")
    ax.set_ylabel("Cow-days")
    ax.set_title("Wood-residual distribution — daily yield minus "
                 "fitted lactation curve\n(per-lactation Wood 1967, "
                 "parity-pooled fallback)")
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(False)
    fig.tight_layout()
    save_figure(fig, "milk_yield_wood_residual_histogram", out_dir)
    return q33, q67


def plot_wood_example_fits(
    wood: pd.DataFrame, fits: pd.DataFrame, out_dir: Path,
    n_examples: int = 9,
    fit_points: pd.DataFrame | None = None,
) -> None:
    """Grid of example per-lactation Wood fits with raw data overlay.

    Args:
        wood: ``daily_milk_yield_wood.csv`` (summer-window residual frame).
        fits: ``wood_curve_fits.csv`` (per-lactation parameters).
        out_dir: Figure output directory.
        n_examples: Number of panels to draw.
        fit_points: Optional full-history cow-days with ``dim``
            attached (e.g. the output of ``attach_dim`` on
            ``daily_milk_yield_full.csv``).  When supplied, the
            scatter uses these rows so the Wood curve sits on top of
            **all** the cow-days that actually informed the fit
            rather than the much smaller summer-window slice.
    """
    import matplotlib.pyplot as plt
    from digimuh.stats_lactation_curve import predict_wood
    setup_figure()

    indiv = fits[(fits["method"] == "per_lactation")
                 & fits["r_squared"].notna()].copy()
    if indiv.empty:
        return

    # Pick a spread: best, median, worst (by R²) plus random in between.
    rng = np.random.default_rng(0)
    indiv = indiv.sort_values("r_squared", ascending=False)
    anchors = indiv.head(3)
    if len(indiv) > 6:
        mid = indiv.iloc[len(indiv) // 2 - 1 : len(indiv) // 2 + 2]
    else:
        mid = pd.DataFrame()
    tail = indiv.tail(3)
    extras = indiv.drop(anchors.index).drop(mid.index, errors="ignore")
    extras = extras.drop(tail.index, errors="ignore")
    random_pick = (
        extras.sample(min(max(n_examples - len(anchors) - len(mid) - len(tail), 0),
                          len(extras)),
                      random_state=0)
        if not extras.empty else pd.DataFrame()
    )
    picks = pd.concat([anchors, mid, random_pick, tail]).head(n_examples)

    log.info("  Plotting %d example Wood fits …", len(picks))

    ncols = min(3, len(picks))
    nrows = int(np.ceil(len(picks) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4.5 * ncols, 3.5 * nrows),
                             squeeze=False)
    axes_flat = axes.flatten()

    scatter_source = fit_points if fit_points is not None else wood
    scatter_label = (
        "full history" if fit_points is not None else "summer window")
    scatter_points = scatter_source.copy()
    if "calving_date" in scatter_points.columns:
        scatter_points["calving_date"] = pd.to_datetime(
            scatter_points["calving_date"])

    for ax, (_, row) in zip(axes_flat, picks.iterrows()):
        aid, cdate = int(row["animal_id"]), pd.to_datetime(row["calving_date"])
        sub = scatter_points[
            (scatter_points["animal_id"] == aid)
            & (scatter_points["calving_date"] == cdate)
        ]
        ax.scatter(sub["dim"], sub["daily_yield_kg"],
                   s=6, alpha=0.35, color=WONG_SKY,
                   edgecolors="none",
                   label=f"cow-days ({scatter_label}, n={len(sub):,})")

        dim_grid = np.linspace(5, 305, 200)
        y_hat = predict_wood(dim_grid, row["a"], row["b"], row["c"])
        ax.plot(dim_grid, y_hat, color=WONG_VERMILLION, linewidth=2,
                label=f"a={row['a']:.1f} b={row['b']:.2f} c={row['c']:.4f}")

        ax.axvline(row["peak_dim"], color=WONG_GREY, linestyle="--",
                   linewidth=1, alpha=0.7)
        ax.set_xlabel("DIM (d)")
        ax.set_ylabel("Daily yield (kg/d)")
        ax.set_title(
            f"Animal {aid}, parity {row['parity']}\n"
            f"R²={row['r_squared']:.3f}, peak {row['peak_dim']:.0f} d "
            f"/ {row['peak_yield']:.1f} kg/d",
            fontsize=9,
        )
        ax.legend(fontsize=7, loc="lower right")
        ax.grid(False)

    for ax in axes_flat[len(picks):]:
        ax.set_axis_off()

    fig.suptitle("Example per-lactation Wood (1967) fits",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "milk_yield_wood_example_fits", out_dir)


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Daily milk yield distribution & classification analysis")
    parser.add_argument("--data", type=Path, default=Path("results/broken_stick"),
                        help="Directory with daily_milk_yield.csv "
                             "(default: results/broken_stick)")
    parser.add_argument("--wood", action="store_true",
                        help="Also fit per-lactation Wood (1967) curves, "
                             "write daily_milk_yield_wood.csv, and produce "
                             "DIM-adjusted residual figures + terciles.")
    parser.add_argument("--db", type=Path, default=None,
                        help="Path to cow.db (only used with --wood if "
                             "calvings.csv is missing from --data).")
    parser.add_argument("--min-points", type=int, default=30,
                        help="Per-lactation point threshold for an "
                             "individual Wood fit (default: 30).  "
                             "Lactations with fewer points fall back to "
                             "the parity-pooled curve.")
    args = parser.parse_args()

    from digimuh.console import setup_logging
    setup_logging()

    df = load_daily_yields(args.data)
    log.info("Loaded %d cow-days from %d animals across %d years",
             len(df), df["animal_id"].nunique(),
             df["year"].nunique())

    our_terciles = (
        float(df["daily_yield_kg"].quantile(0.33)),
        float(df["daily_yield_kg"].quantile(0.67)),
    )

    print_summary(df, our_terciles)

    plot_pooled_histogram(df, our_terciles, args.data)
    plot_per_year_histogram(df, our_terciles, args.data)

    if args.wood:
        from digimuh.stats_lactation_curve import (
            load_calvings, compute_wood_residuals,
            load_daily_yields_for_fitting, attach_dim,
        )
        calvings = load_calvings(args.data, db_path=args.db)
        fit_frame = load_daily_yields_for_fitting(args.data)
        if fit_frame is not None:
            log.info("Wood fit source: daily_milk_yield_full.csv "
                     "(%d cow-days across full HerdePlus history)",
                     len(fit_frame))
        else:
            log.info("Wood fit source: daily_milk_yield.csv "
                     "(summer-window only — re-run digimuh-extract to "
                     "produce daily_milk_yield_full.csv for better "
                     "lactation-curve coverage)")
        wood, fits = compute_wood_residuals(
            df, calvings, min_points=args.min_points,
            fit_yields=fit_frame,
        )
        wood.to_csv(resolve_output(args.data, "daily_milk_yield_wood.csv"), index=False)
        fits.to_csv(resolve_output(args.data, "wood_curve_fits.csv"), index=False)
        log.info("  → wrote daily_milk_yield_wood.csv (%d rows) "
                 "and wood_curve_fits.csv (%d fits)",
                 len(wood), len(fits))
        print_wood_summary(wood, fits)
        plot_residual_histogram(wood, args.data)

        # Build the fit-frame scatter source so the example-fit grid
        # shows ALL the cow-days that actually informed each Wood fit
        # (not just the summer-window slice used for residuals).
        fit_points = (
            attach_dim(fit_frame, calvings)
            if fit_frame is not None else None
        )
        plot_wood_example_fits(wood, fits, args.data,
                               fit_points=fit_points)

        # ── Stratified TNF × yield (needs tnf_yield.csv from stats_runner) ──
        tnf_path = resolve_input(args.data, "tnf_yield.csv")
        if tnf_path.exists():
            _run_stratified_tnf(wood, tnf_path, args.data)
        else:
            log.info("  tnf_yield.csv not found at %s — run "
                     "`digimuh-stats` first to enable the stratified "
                     "TNF × yield analysis.", tnf_path)

    log.info("Figures written to %s", args.data)


def _run_stratified_tnf(
    wood: pd.DataFrame, tnf_path: Path, out_dir: Path,
) -> None:
    """Per-yield-class Spearman correlations of TNF vs yield.

    Reports three response variants per (class × predictor):
    daily_yield_kg, relative_yield (daily / cow P95), and
    yield_residual (DIM-adjusted Wood residual).  The residual is
    the clean heat-stress signal because the within-cow lactation
    decline has been removed; the raw variants are retained for
    comparability with literature.
    """
    from digimuh.console import section, result_table, kv, stars_styled
    from digimuh.stats_production import (
        classify_cow_years_by_wood_residual,
        compute_tnf_yield_by_class,
        tnf_yield_correlations_by_class,
    )
    from digimuh.viz_production import plot_tnf_yield_by_class

    tnf_yield = pd.read_csv(tnf_path)

    class_table, (q33, q67) = classify_cow_years_by_wood_residual(wood)
    if class_table.empty:
        log.info("  No cow-years survived classification, skipping "
                 "stratified TNF × yield")
        return

    tnf_by_class = compute_tnf_yield_by_class(tnf_yield, class_table, wood=wood)
    corr = tnf_yield_correlations_by_class(tnf_by_class)

    class_table.to_csv(resolve_output(out_dir, "yield_class_per_cow_year.csv"), index=False)
    tnf_by_class.to_csv(resolve_output(out_dir, "tnf_yield_by_class.csv"), index=False)
    corr.to_csv(resolve_output(out_dir, "tnf_yield_correlations_by_class.csv"), index=False)

    section("Cow-year yield-class assignment (Wood-residual means)",
            f"tercile boundaries on cow-year means: "
            f"Q33 = {q33:+.2f} kg/d, Q67 = {q67:+.2f} kg/d")
    kv("Cow-years classified", len(class_table))
    for cls in ("low", "middle", "high"):
        sub = class_table[class_table["yield_class"] == cls]
        kv(f"  {cls:6s}", f"{len(sub)} cow-years "
           f"({sub['animal_id'].nunique()} unique animals)")

    response_labels = {
        "daily_yield_kg": ("daily yield (kg/d)", "kg/d per 1.0 TNF"),
        "relative_yield": ("relative yield (/P95)",
                           "(y/P95) per 1.0 TNF"),
        "yield_residual": ("DIM-adj. residual (kg/d)",
                           "kg/d per 1.0 TNF"),
    }

    from digimuh.stats_core import p_to_stars
    for resp, (resp_label, slope_unit) in response_labels.items():
        section(f"TNF × {resp_label} by class — Spearman",
                "low / middle / high producers; pooled = all cow-days"
                + (" · DIM-adjusted (heat-stress-only)"
                   if resp == "yield_residual" else ""))
        rows = []
        for cls in ("low", "middle", "high", "pooled"):
            for pred in ("thi_tnf", "temp_tnf"):
                stat = corr[(corr["yield_class"] == cls)
                            & (corr["predictor"] == pred)
                            & (corr["response"] == resp)]
                if stat.empty:
                    continue
                s = stat.iloc[0]
                rows.append([
                    cls, pred,
                    int(s["n"]), int(s["n_animals"]),
                    f"{s['rs']:+.3f}" if np.isfinite(s["rs"]) else "—",
                    f"{s['p']:.2e}" if np.isfinite(s["p"]) else "—",
                    (stars_styled(p_to_stars(s["p"]))
                     if np.isfinite(s["p"]) else ""),
                    f"{s['slope']:+.2f}" if np.isfinite(s["slope"]) else "—",
                ])
        result_table(
            f"TNF × {resp_label} by class",
            ["Class", "Predictor", "n days", "n cows",
             "rs", "p", "Sig.", f"Slope ({slope_unit})"],
            rows,
        )

    plot_tnf_yield_by_class(tnf_by_class, corr, out_dir,
                            response="daily_yield_kg")
    plot_tnf_yield_by_class(tnf_by_class, corr, out_dir,
                            response="relative_yield")
    plot_tnf_yield_by_class(tnf_by_class, corr, out_dir,
                            response="yield_residual")

    # ── Crossing-day raincloud + daily-mean climate fit ─────
    _run_crossing_and_climate_analyses(tnf_by_class, out_dir)

    # ── MLP composition × climate (thin-milk hypothesis) ────
    _run_milk_composition_analysis(
        wood=wood, tnf_yield=tnf_yield,
        class_table=class_table, out_dir=out_dir,
    )


def _run_crossing_and_climate_analyses(
    tnf_by_class: pd.DataFrame, out_dir: Path,
) -> None:
    """Binary crossed-vs-not raincloud + one daily-mean scatter fit.

    Extends the stratified TNF analysis with two questions the
    reviewers may prefer over a continuous TNF:

    1. *Does a cow produce differently on days where her
       breakpoint was actually crossed vs days where it was not?*
       — raincloud per climate predictor (THI, barn-temp) with
       a Mann-Whitney median-difference test.
    2. *Does the daily mean climate itself predict yield?* —
       single scatter + OLS fit per climate predictor.
    """
    from digimuh.console import section, result_table, kv, stars_styled
    from digimuh.stats_core import p_to_stars
    from digimuh.stats_production import (
        compute_daily_crossing_flags,
        attach_daily_climate_means,
        attach_crossing_flags,
        crossing_day_comparison,
        daily_climate_vs_yield_correlations,
    )
    from digimuh.viz_production import (
        plot_crossing_day_raincloud,
        plot_daily_climate_vs_yield,
    )

    crossing_path = resolve_input(out_dir, "crossing_times.csv")
    if not crossing_path.exists():
        log.info("  crossing_times.csv not found — skipping crossing-day "
                 "raincloud")
        return

    crossing_times = pd.read_csv(crossing_path)
    flags = compute_daily_crossing_flags(crossing_times)
    df = attach_crossing_flags(tnf_by_class, flags)

    # mean_barn_temp is a newer column — fall back to rumen_barn.csv.
    rumen_path = resolve_input(out_dir, "rumen_barn.csv")
    if "mean_barn_temp" not in df.columns and rumen_path.exists():
        log.info("  mean_barn_temp absent from tnf_yield.csv — "
                 "re-deriving from rumen_barn.csv …")
        rumen = pd.read_csv(rumen_path)
        df = attach_daily_climate_means(df, rumen)

    df.to_csv(resolve_output(out_dir, "tnf_yield_by_class.csv"), index=False)

    # ── Raincloud: crossed vs not crossed ───────────────────
    comparison = crossing_day_comparison(df, response="yield_residual")
    comparison.to_csv(resolve_output(out_dir, "crossing_day_comparison.csv"), index=False)

    section("Crossed-day vs not-crossed-day yield (DIM-adjusted residual)",
            "Mann-Whitney median-difference test per climate predictor")
    rows = []
    for _, r in comparison.iterrows():
        if not np.isfinite(r["p"]):
            continue
        rows.append([
            r["group"], r["predictor"],
            int(r["n_yes"]), int(r["n_no"]),
            f"{r['median_yes']:+.2f}",
            f"{r['median_no']:+.2f}",
            f"{r['median_diff']:+.2f}",
            f"{r['p']:.2e}",
            stars_styled(p_to_stars(r["p"])),
        ])
    if rows:
        result_table(
            "Crossing-day comparison",
            ["Group", "Predictor", "n crossed", "n not",
             "Median (crossed)", "Median (not)",
             "Δ median (kg/d)", "p", "Sig."],
            rows,
        )

    for group in ("pooled", "low", "middle", "high"):
        plot_crossing_day_raincloud(
            df, comparison, out_dir,
            response="yield_residual", group=group,
        )

    # ── Daily mean climate vs yield — one fit ───────────────
    clim_corr = daily_climate_vs_yield_correlations(df,
                                                    response="yield_residual")
    clim_corr.to_csv(resolve_output(out_dir, "daily_climate_vs_yield.csv"), index=False)

    section("Daily mean climate vs DIM-adjusted yield (pooled)",
            "single OLS + Spearman per climate predictor — no class split")
    if not clim_corr.empty:
        crows = []
        for _, r in clim_corr.iterrows():
            if not np.isfinite(r.get("rs", np.nan)):
                continue
            crows.append([
                r["predictor"],
                int(r["n"]), int(r["n_animals"]),
                f"{r['rs']:+.3f}", f"{r['p']:.2e}",
                stars_styled(p_to_stars(r["p"])),
                f"{r['slope']:+.3f}",
            ])
        if crows:
            result_table(
                "Daily mean × residual",
                ["Predictor", "n days", "n cows",
                 "rs", "p", "Sig.", "Slope (kg/d per unit)"],
                crows,
            )

    plot_daily_climate_vs_yield(df, clim_corr, out_dir,
                                response="yield_residual")


def _run_milk_composition_analysis(
    wood: pd.DataFrame, tnf_yield: pd.DataFrame,
    class_table: pd.DataFrame, out_dir: Path,
) -> None:
    """MLP composition × climate — tests the thin-milk hypothesis.

    Joins monthly MLP test-day rows (fat %, protein %, SCC, ECM, …)
    to the cow-day climate / residual / yield-class frames, then
    runs Spearman correlations between each composition metric and
    each climate predictor (pooled and per yield class).  Also
    renders a dilution-focused 2×2 figure and a signed-rs heatmap.
    """
    from digimuh.console import section, result_table, kv, stars_styled
    from digimuh.stats_core import p_to_stars
    from digimuh.stats_milk_composition import (
        load_mlp_test_days, merge_mlp_with_cowday,
        mlp_climate_correlations, thin_milk_verdict,
        MLP_RESPONSES, CLIMATE_PREDICTORS,
    )
    from digimuh.viz_milk_composition import (
        plot_thin_milk_hypothesis, plot_composition_heatmap,
    )

    try:
        mlp = load_mlp_test_days(out_dir)
    except FileNotFoundError:
        log.info("  mlp_test_days.csv not found — skipping milk-composition "
                 "analysis")
        return

    # Back-fill climate (mean_barn_temp specifically) from rumen_barn.csv
    # when legacy tnf_yield.csv lacks the column.
    rumen_path = resolve_input(out_dir, "rumen_barn.csv")
    rumen = pd.read_csv(rumen_path) if rumen_path.exists() else None

    merged = merge_mlp_with_cowday(mlp, wood, tnf_yield,
                                   class_table=class_table,
                                   rumen=rumen)
    if merged.empty:
        log.info("  No MLP test-days match the cow-day climate frame — "
                 "skipping milk-composition analysis")
        return

    merged.to_csv(resolve_output(out_dir, "mlp_composition_by_cowday.csv"), index=False)

    corr = mlp_climate_correlations(merged)
    corr.to_csv(resolve_output(out_dir, "mlp_composition_correlations.csv"), index=False)

    verdict = thin_milk_verdict(corr, predictor="mean_thi")

    section("MLP milk composition × climate (thin-milk hypothesis)",
            f"{len(merged):,} MLP test-days matched to cow-day climate "
            f"({merged['animal_id'].nunique()} animals)")
    kv("Hypothesis",
       "hot day → more milk volume, lower fat / protein % (dilution)")
    kv("Pooled rs (milk volume vs mean THI)",
       f"{verdict['rs_milk_volume']:+.3f}"
       if np.isfinite(verdict.get("rs_milk_volume", np.nan)) else "n/a")
    kv("Pooled rs (fat % vs mean THI)",
       f"{verdict['rs_fat_percent']:+.3f}"
       if np.isfinite(verdict.get("rs_fat_percent", np.nan)) else "n/a")
    kv("Pooled rs (protein % vs mean THI)",
       f"{verdict['rs_protein_percent']:+.3f}"
       if np.isfinite(verdict.get("rs_protein_percent", np.nan)) else "n/a")
    kv("Pooled rs (F/E ratio vs mean THI)",
       f"{verdict['rs_fat_protein_ratio']:+.3f}"
       if np.isfinite(verdict.get("rs_fat_protein_ratio", np.nan)) else "n/a")
    kv("Verdict (hypothesis test)", verdict["verdict"])

    # Full per-(group × predictor × response) table — one big rich table.
    section("MLP × climate correlations (full)",
            "rs, p, slope per (yield class × climate predictor × MLP metric)")
    for predictor, pred_label in CLIMATE_PREDICTORS:
        rows = []
        for group in ("low", "middle", "high", "pooled"):
            for resp, resp_label, unit in MLP_RESPONSES:
                stat = corr[(corr["group"] == group)
                            & (corr["predictor"] == predictor)
                            & (corr["response"] == resp)]
                if stat.empty:
                    continue
                s = stat.iloc[0]
                if not np.isfinite(s["rs"]):
                    continue
                rows.append([
                    group, resp_label, unit,
                    int(s["n"]), int(s["n_animals"]),
                    f"{s['rs']:+.3f}",
                    f"{s['p']:.2e}",
                    stars_styled(p_to_stars(s["p"])),
                    f"{s['slope']:+.4f}",
                ])
        if rows:
            result_table(
                f"Predictor: {pred_label}",
                ["Class", "MLP metric", "Unit", "n days", "n cows",
                 "rs", "p", "Sig.", "Slope (metric per unit)"],
                rows,
            )

    plot_thin_milk_hypothesis(merged, corr, out_dir, predictor="mean_thi")
    plot_composition_heatmap(corr, out_dir, predictor="mean_thi")

    # ── Dilution partition: pure water vs rumen suppression ─
    from digimuh.stats_milk_composition import (
        compute_dilution_partition, dilution_partition_summary,
    )
    from digimuh.viz_milk_composition import plot_dilution_partition

    partitioned = compute_dilution_partition(merged)
    if "fat_percent_diluted" not in partitioned.columns:
        return
    partitioned.to_csv(resolve_output(out_dir, "mlp_dilution_partition.csv"), index=False)

    dil_summary = dilution_partition_summary(partitioned, predictor="mean_thi")
    dil_summary.to_csv(resolve_output(out_dir, "mlp_dilution_partition_summary.csv"),
                       index=False)

    section("Dilution partition — how much is just added water?",
            "observed composition slope vs the slope you'd get from "
            "pure-volume dilution alone (per-cow reference fat/protein kg)")
    rows = []
    nutrients_in_summary = sorted(dil_summary["nutrient"].unique().tolist(),
                                  key=["fat", "protein", "lactose"].index)
    for nutrient in nutrients_in_summary:
        for comp in ("observed", "dilution_predicted", "rumen_residual"):
            s = dil_summary[(dil_summary["nutrient"] == nutrient)
                            & (dil_summary["component"] == comp)]
            if s.empty:
                continue
            s = s.iloc[0]
            rows.append([
                nutrient.capitalize(),
                comp.replace("_", " "),
                int(s["n"]),
                f"{s['rs']:+.3f}",
                f"{s['p']:.2e}",
                stars_styled(p_to_stars(s["p"])),
                f"{s['slope']:+.4f}",
            ])
    if rows:
        result_table(
            "Dilution partition (vs mean barn THI)",
            ["Nutrient", "Component", "n", "rs", "p", "Sig.",
             "Slope (%/unit THI)"],
            rows,
        )

    # Short interpretive sentence
    def _slope(nutrient, comp):
        r = dil_summary[(dil_summary["nutrient"] == nutrient)
                        & (dil_summary["component"] == comp)]
        return float(r["slope"].iloc[0]) if not r.empty else np.nan

    for nutrient in nutrients_in_summary:
        obs = _slope(nutrient, "observed")
        dil = _slope(nutrient, "dilution_predicted")
        rum = _slope(nutrient, "rumen_residual")
        if not all(np.isfinite(v) for v in (obs, dil, rum)):
            continue
        # If |observed| is smaller than either component in absolute
        # terms, the two components are essentially cancelling each
        # other — reporting "dilution explains XYZ%" is misleading.
        # Print the component-cancellation form instead.
        if abs(obs) < 0.5 * max(abs(dil), abs(rum)):
            kv(f"  {nutrient.capitalize()}: near-zero net slope",
               f"dilution {dil:+.4f} offset by rumen {rum:+.4f} "
               f"→ observed {obs:+.4f}")
            continue
        explained = 100 * dil / obs
        kv(f"  {nutrient.capitalize()}: dilution explains",
           f"{explained:+.0f}% of the observed slope "
           f"(obs {obs:+.4f} = dil {dil:+.4f} + rumen {rum:+.4f})")

    plot_dilution_partition(partitioned, dil_summary, out_dir,
                            predictor="mean_thi")


if __name__ == "__main__":
    main()

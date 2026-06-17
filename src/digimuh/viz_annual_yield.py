#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — viz_annual_yield                                      ║
# ║  « yield-vs-month lines + heat-stress duration response »        ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Median-cow daily yield as line plots over the calendar year     ║
# ║  (day-of-year x value, one line per lactation / DIM band, with   ║
# ║  bootstrap 95% CI bands) plus the heat-stress duration response  ║
# ║  figure (yield delta vs consecutive stress days).                ║
# ║                                                                  ║
# ║  SVG + PNG via viz_base.save_figure; CSV companions alongside.   ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Year-round yield line plots and the heat-stress duration-response figure."""

from __future__ import annotations

import calendar
import logging
from pathlib import Path

import pandas as pd

from digimuh.constants import (
    WONG_BLUE,
    WONG_GREEN,
    WONG_GREY,
    WONG_ORANGE,
    WONG_PINK,
    WONG_SKY,
    WONG_VERMILLION,
    WONG_YELLOW,
)
from digimuh.stats_annual_yield import LACTATION_CAP
from digimuh.viz_base import save_figure, setup_figure

log = logging.getLogger("digimuh.viz.annual_yield")

_MONTH_STARTS = [pd.Timestamp(2023, m, 1).dayofyear for m in range(1, 13)]
_MONTH_LABELS = [calendar.month_abbr[m] for m in range(1, 13)]

# Discrete line colours for the lactation bands (Wong, ordered L1…L8+).
_LACTATION_COLOURS = (
    WONG_BLUE, WONG_ORANGE, WONG_GREEN, WONG_VERMILLION,
    WONG_SKY, WONG_PINK, WONG_YELLOW, WONG_GREY,
)


# ─────────────────────────────────────────────────────────────
#  « Median-yield line plots (CI bands) »
# ─────────────────────────────────────────────────────────────

def plot_yield_lines(
    grid: pd.DataFrame,
    value_kind: str,
    y_kind: str,
    name: str,
    out_dir: Path,
) -> None:
    """Median-cow yield over the year as CI-banded lines per band.

    One line per y-band (median yield vs day-of-year) with a bootstrap
    95% CI ribbon.  Lactation bands get the discrete Wong palette and a
    legend; DIM bands get a viridis gradient and a colourbar (too many
    bands for a legend).

    Args:
        grid: Long-form grid from
            :func:`stats_annual_yield.build_median_surface`
            (``doy_bin``, ``y_bin``, ``median_value``, ``ci_lo``,
            ``ci_hi``, ``n_cows``).
        value_kind: ``"z"`` (z-scored, zero reference line) or
            ``"kg"`` (raw kg).
        y_kind: ``"lactation"`` or ``"dim"`` — sets band colouring.
        name: Output stem (also the CSV companion stem).
        out_dir: Results directory (routed to ``07_annual_yield``).
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    setup_figure()
    fig, ax = plt.subplots(figsize=(10, 5.5))
    if value_kind == "z":
        ax.axhline(0.0, color=WONG_GREY, lw=1.0, ls="--", zorder=1)

    y_bins = sorted(int(v) for v in grid["y_bin"].dropna().unique())

    if y_kind == "lactation":
        for i, yb in enumerate(y_bins):
            sub = _band(grid, yb)
            if sub.empty:
                continue
            colour = _LACTATION_COLOURS[i % len(_LACTATION_COLOURS)]
            label = f"L{yb}+" if yb >= LACTATION_CAP else f"L{yb}"
            _draw_band(ax, sub, colour, label, ci_alpha=0.18, lw=1.8)
        ax.legend(title="lactation", frameon=False, ncol=2, fontsize=8,
                  loc="best")
    else:
        cmap = plt.get_cmap("viridis")
        norm = mpl.colors.Normalize(vmin=min(y_bins), vmax=max(y_bins))
        for yb in y_bins:
            sub = _band(grid, yb)
            if sub.empty:
                continue
            _draw_band(ax, sub, cmap(norm(yb)), None, ci_alpha=0.10, lw=1.4)
        sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, pad=0.02).set_label("days in milk (DIM)")

    ax.set_xticks(_MONTH_STARTS)
    ax.set_xticklabels(_MONTH_LABELS)
    ax.set_xlim(1, 366)
    ax.set_xlabel("month of year")
    ax.set_ylabel("median z-scored daily yield" if value_kind == "z"
                  else "median daily yield (kg)")
    ax.set_title("Median-cow daily milk yield over the year")

    fig.tight_layout()
    save_figure(fig, name, out_dir)
    _write_csv(grid, name, out_dir)
    log.info("wrote yield-line figure %s", name)


def _band(grid: pd.DataFrame, y_bin: int) -> pd.DataFrame:
    """One y-band's cells, ordered by day-of-year, valid medians only."""
    return (grid[grid["y_bin"] == y_bin]
            .dropna(subset=["median_value"])
            .sort_values("doy_bin"))


def _draw_band(ax, sub: pd.DataFrame, colour, label: str | None,
               ci_alpha: float, lw: float) -> None:
    """Plot one band's median line plus its 95% CI ribbon."""
    has_ci = sub[["ci_lo", "ci_hi"]].notna().all(axis=1)
    if has_ci.any():
        s = sub[has_ci]
        ax.fill_between(s["doy_bin"], s["ci_lo"], s["ci_hi"],
                        color=colour, alpha=ci_alpha, lw=0, zorder=2)
    ax.plot(sub["doy_bin"], sub["median_value"], color=colour, lw=lw,
            label=label, zorder=4)


# ─────────────────────────────────────────────────────────────
#  « Heat-stress duration response »
# ─────────────────────────────────────────────────────────────

def plot_heatstress_duration(
    summary: pd.DataFrame,
    deltas: pd.DataFrame,
    value_tag: str,
    name: str,
    out_dir: Path,
    tests: pd.DataFrame | None = None,
) -> None:
    """Yield delta vs consecutive heat-stress days: median+CI over rainclouds.

    Args:
        summary: From :func:`stats_annual_yield.summarise_duration_deltas`.
        deltas: Per-cow medians from
            :func:`stats_annual_yield.aggregate_per_cow`.
        value_tag: ``"residual"`` (DIM-adjusted) or ``"kg"``.
        name: Output stem.
        out_dir: Results directory.
        tests: Optional per-streak-day significance table from
            :func:`stats_annual_yield.test_duration_medians` (already
            filtered to this metric); annotates BH-FDR stars on the median.
    """
    import matplotlib.pyplot as plt

    setup_figure()
    col = f"delta_{value_tag}"
    med = f"median_delta_{value_tag}"
    lo, hi = f"ci_lo_{value_tag}", f"ci_hi_{value_tag}"
    ks = summary["streak_day"].to_numpy()
    unit = "Wood residual (kg/d)" if value_tag == "residual" else "raw kg/d"

    fig, (ax_t, ax_b) = plt.subplots(
        2, 1, figsize=(8, 7), sharex=True,
        gridspec_kw={"height_ratios": [2, 1.4], "hspace": 0.08})

    # ── top: median Δ + 95% CI (the headline, zoomed) ──────
    ax_t.axhline(0.0, color=WONG_GREY, lw=1.0, ls="--", zorder=1)
    ax_t.fill_between(summary["streak_day"], summary[lo], summary[hi],
                      color=WONG_VERMILLION, alpha=0.20, zorder=2, label="95% CI")
    ax_t.plot(summary["streak_day"], summary[med], "-o", color=WONG_VERMILLION,
              lw=2, zorder=4, label="median Δ")
    span = float(summary[hi].max() - summary[lo].min())
    pad = 0.25 * span if span > 0 else 1.0
    ax_t.set_ylim(summary[lo].min() - pad, summary[hi].max() + pad)
    for _, r in summary.iterrows():
        ax_t.annotate(f"n={int(r['n'])}", (r["streak_day"], r[hi]),
                      textcoords="offset points", xytext=(0, 6),
                      ha="center", fontsize=8, color=WONG_GREY)
    if tests is not None:
        star_map = dict(zip(tests["streak_day"], tests["stars"]))
        for _, r in summary.iterrows():
            stars = star_map.get(r["streak_day"], "")
            if not stars or stars == "n.s.":
                continue
            ax_t.annotate(stars, (r["streak_day"], r[med]),
                          textcoords="offset points", xytext=(0, -14),
                          ha="center", fontsize=12, fontweight="bold",
                          color="black")
    ax_t.set_ylabel(f"median Δ yield — {unit}")
    ax_t.set_title("Milk-yield response to heat-stress duration")
    ax_t.legend(frameon=False, loc="lower left")

    # ── bottom: per-day distribution (box, 5–95% whiskers) ─
    data = [deltas.loc[deltas["streak_day"] == k, col].dropna().to_numpy()
            for k in ks]
    ax_b.axhline(0.0, color=WONG_GREY, lw=1.0, ls="--", zorder=1)
    bp = ax_b.boxplot(data, positions=ks, widths=0.6, whis=(5, 95),
                      showfliers=False, patch_artist=True, zorder=2)
    for patch in bp["boxes"]:
        patch.set(facecolor=WONG_SKY, alpha=0.55, edgecolor=WONG_BLUE)
    for med_ln in bp["medians"]:
        med_ln.set(color=WONG_BLUE, lw=1.5)
    ax_b.set_xlabel("consecutive heat-stress days")
    ax_b.set_ylabel(f"per-cow median Δ — {unit}")
    ax_b.set_xticks(ks)

    save_figure(fig, name, out_dir)
    _write_csv(summary, name, out_dir)
    log.info("wrote duration-response %s", name)


# ─────────────────────────────────────────────────────────────
#  « helpers »
# ─────────────────────────────────────────────────────────────

def _write_csv(df: pd.DataFrame, name: str, out_dir: Path) -> None:
    """Write the figure's data companion next to the SVG/PNG."""
    from digimuh.paths import resolve_output
    dest = resolve_output(out_dir, f"{name}.csv")
    df.to_csv(dest, index=False)

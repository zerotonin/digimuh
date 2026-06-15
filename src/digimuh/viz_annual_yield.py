#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — viz_annual_yield                                      ║
# ║  « heatmaps of median yield + heat-stress duration response »    ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  2-D heatmaps (day-of-year x lactation / DIM, colour = median    ║
# ║  yield, z-scored and raw kg) plus the heat-stress duration       ║
# ║  response figure (yield delta vs consecutive stress days).       ║
# ║                                                                  ║
# ║  SVG + PNG via viz_base.save_figure; CSV companions alongside.   ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Heatmaps and the heat-stress duration-response figure."""

from __future__ import annotations

import calendar
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.viz_base import save_figure, setup_figure

log = logging.getLogger("digimuh.viz.annual_yield")

_MONTH_STARTS = [pd.Timestamp(2023, m, 1).dayofyear for m in range(1, 13)]
_MONTH_LABELS = [calendar.month_abbr[m] for m in range(1, 13)]


# ─────────────────────────────────────────────────────────────
#  « Median-yield heatmaps »
# ─────────────────────────────────────────────────────────────

def plot_yield_heatmap(
    grid: pd.DataFrame,
    value_kind: str,
    y_kind: str,
    name: str,
    out_dir: Path,
) -> None:
    """Render one median-yield heatmap (day-of-year × lactation/DIM).

    Args:
        grid: Long-form grid from
            :func:`stats_annual_yield.build_median_surface`
            (``doy_bin``, ``y_bin``, ``median_value``, ``n_cows``).
        value_kind: ``"z"`` (diverging, centred on 0) or
            ``"kg"`` (sequential).
        y_kind: ``"lactation"`` or ``"dim"`` — sets the y-axis label.
        name: Output stem (also the CSV companion stem).
        out_dir: Results directory (routed to ``07_annual_yield``).
    """
    import matplotlib.pyplot as plt

    setup_figure()
    mat = grid.pivot(index="y_bin", columns="doy_bin", values="median_value")
    mat = mat.sort_index()
    x = mat.columns.to_numpy(dtype=float)
    y = mat.index.to_numpy(dtype=float)

    if value_kind == "z":
        vmax = float(np.nanpercentile(np.abs(mat.to_numpy()), 98))
        cmap, vmin, vctr = "RdBu_r", -vmax, 0.0
        cbar_label = "median z-scored daily yield"
    else:
        vmin = float(np.nanpercentile(mat.to_numpy(), 2))
        vmax = float(np.nanpercentile(mat.to_numpy(), 98))
        cmap, vctr = "viridis", None
        cbar_label = "median daily yield (kg)"

    fig, ax = plt.subplots(figsize=(10, 5))
    norm = (plt.matplotlib.colors.TwoSlopeNorm(vmin=vmin, vcenter=vctr, vmax=vmax)
            if vctr is not None else
            plt.matplotlib.colors.Normalize(vmin=vmin, vmax=vmax))
    mesh = ax.pcolormesh(x, y, mat.to_numpy(), cmap=cmap, norm=norm,
                         shading="nearest")
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(cbar_label)

    ax.set_xticks(_MONTH_STARTS)
    ax.set_xticklabels(_MONTH_LABELS)
    ax.set_xlabel("month of year")
    if y_kind == "lactation":
        ax.set_ylabel("lactation number")
        ax.set_yticks(sorted(int(v) for v in y))
    else:
        ax.set_ylabel("days in milk (DIM)")
    ax.set_title("Median-cow daily milk yield")

    fig.tight_layout()
    save_figure(fig, name, out_dir)
    _write_csv(grid, name, out_dir)
    log.info("wrote heatmap %s", name)


# ─────────────────────────────────────────────────────────────
#  « Heat-stress duration response »
# ─────────────────────────────────────────────────────────────

def plot_heatstress_duration(
    summary: pd.DataFrame,
    deltas: pd.DataFrame,
    value_tag: str,
    name: str,
    out_dir: Path,
) -> None:
    """Yield delta vs consecutive heat-stress days: median+CI over rainclouds.

    Args:
        summary: From :func:`stats_annual_yield.summarise_duration_deltas`.
        deltas: From
            :func:`stats_annual_yield.compute_heatstress_duration_deltas`.
        value_tag: ``"residual"`` (DIM-adjusted) or ``"kg"``.
        name: Output stem.
        out_dir: Results directory.
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
    ax_t.axhline(0.0, color="#888888", lw=1.0, ls="--", zorder=1)
    ax_t.fill_between(summary["streak_day"], summary[lo], summary[hi],
                      color="#D55E00", alpha=0.20, zorder=2, label="95% CI")
    ax_t.plot(summary["streak_day"], summary[med], "-o", color="#D55E00",
              lw=2, zorder=4, label="median Δ")
    span = float(summary[hi].max() - summary[lo].min())
    pad = 0.25 * span if span > 0 else 1.0
    ax_t.set_ylim(summary[lo].min() - pad, summary[hi].max() + pad)
    for _, r in summary.iterrows():
        ax_t.annotate(f"n={int(r['n'])}", (r["streak_day"], r[hi]),
                      textcoords="offset points", xytext=(0, 6),
                      ha="center", fontsize=8, color="#444444")
    ax_t.set_ylabel(f"median Δ yield — {unit}")
    ax_t.set_title("Milk-yield response to heat-stress duration")
    ax_t.legend(frameon=False, loc="lower left")

    # ── bottom: per-day distribution (box, 5–95% whiskers) ─
    data = [deltas.loc[deltas["streak_day"] == k, col].dropna().to_numpy()
            for k in ks]
    ax_b.axhline(0.0, color="#888888", lw=1.0, ls="--", zorder=1)
    bp = ax_b.boxplot(data, positions=ks, widths=0.6, whis=(5, 95),
                      showfliers=False, patch_artist=True, zorder=2)
    for patch in bp["boxes"]:
        patch.set(facecolor="#56B4E9", alpha=0.55, edgecolor="#33647d")
    for med_ln in bp["medians"]:
        med_ln.set(color="#33647d", lw=1.5)
    ax_b.set_xlabel("consecutive heat-stress days")
    ax_b.set_ylabel(f"Δ yield — {unit}")
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

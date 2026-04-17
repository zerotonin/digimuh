#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — viz_temporal                                         ║
# ║  « circadian, daily exceedance, and crossing activation plots »║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Circadian null model, THI daily exceedance profile, and        ║
# ║  breakpoint crossing raster + raincloud.                        ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Temporal analysis figures: circadian and crossing activation."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.constants import COLOURS
from digimuh.viz_base import setup_figure, save_figure
from digimuh.paths import resolve_input

log = logging.getLogger("digimuh.viz")

# ─────────────────────────────────────────────────────────────
#  « rumen circadian null model »
# ─────────────────────────────────────────────────────────────

_MILKING_WINDOWS = [(4, 7), (16, 19)]
_KDE_COLOUR = "#009E73"


def _hourly_by_day_type(df: pd.DataFrame, var: str) -> dict[str, pd.DataFrame]:
    """Aggregate a variable by clock hour within each day_type.

    Returns a dict keyed by ``"cool"`` / ``"stress"`` with columns
    ``mean``, ``std``, ``count``, ``sem`` indexed by hour.
    """
    out: dict[str, pd.DataFrame] = {}
    if var not in df.columns:
        return out
    for day_type in ("cool", "stress"):
        sub = df[df["day_type"] == day_type]
        if sub.empty:
            continue
        h = sub.groupby("hour")[var].agg(["mean", "std", "count"])
        h["sem"] = h["std"] / np.sqrt(h["count"])
        out[day_type] = h
    return out


def _draw_cool_stress_panel(ax, hourly: dict[str, pd.DataFrame]) -> None:
    """Draw cool-day (blue) and stress-day (vermillion) mean ± SEM ribbons."""
    for day_type, colour, label in [
        ("cool",   COLOURS["below_bp"], "Cool days (THI < bp all day)"),
        ("stress", COLOURS["above_bp"], "Heat stress days (THI exceeded bp)"),
    ]:
        h = hourly.get(day_type)
        if h is None:
            continue
        ax.fill_between(h.index, h["mean"] - h["sem"], h["mean"] + h["sem"],
                        alpha=0.2, color=colour)
        ax.plot(h.index, h["mean"], color=colour, linewidth=2,
                marker="o", markersize=4, label=label)
    for s, e in _MILKING_WINDOWS:
        ax.axvspan(s, e, alpha=0.08, color="#999", zorder=0)


def _overlay_crossing_density(ax, kde: tuple | None, legend_label: str):
    """Add crossing-density KDE ribbon on a twinx; returns the twin axis."""
    if kde is None:
        return None
    x_k, y_k, n_k = kde
    ax2 = ax.twinx()
    ax2.fill_between(x_k, 0, y_k, alpha=0.12, color=_KDE_COLOUR, zorder=0)
    ax2.plot(x_k, y_k, color=_KDE_COLOUR, linewidth=1.5, alpha=0.7,
             label=f"{legend_label} (n={n_k})")
    ax2.set_ylabel("Crossing density", color=_KDE_COLOUR)
    ax2.tick_params(axis="y", labelcolor=_KDE_COLOUR)
    ax2.set_ylim(0, y_k.max() * 2.5)
    ax2.grid(False)
    return ax2


def _merge_legends(primary, secondary, **kw):
    lines1, labels1 = primary.get_legend_handles_labels()
    if secondary is not None:
        lines2, labels2 = secondary.get_legend_handles_labels()
        primary.legend(lines1 + lines2, labels1 + labels2, **kw)
    else:
        primary.legend(**kw)


def plot_circadian_null_model(out_dir: Path) -> None:
    """Plot the rumen-temperature circadian null model (2×2 grid).

    Four panels, all plotted over clock hour on cool vs stress days
    (stratified by each animal's own THI breakpoint):

    * **A** — Rumen temperature, with THI-crossing density overlay.
    * **B** — Stress − cool Δ rumen temperature (heat-stress dose),
      with THI-crossing density overlay.
    * **C** — Barn THI on cool vs stress days + THI-crossing density.
    * **D** — Barn temperature on cool vs stress days +
      barn-temp-crossing density.

    Shaded regions:

    * Coloured ribbon around each mean line = ±1 SEM across animals
      (cool=blue, stress=vermillion, Δ=pink).
    * Grey vertical bands at 04:00–07:00 and 16:00–19:00 = excluded
      milking windows.
    * Green ribbon on the right y-axis = kernel density estimate of
      breakpoint-crossing clock times (smoothed histogram of when,
      during the 24 h cycle, each cow's individual breakpoint is
      crossed upward).

    Grid lines are disabled to keep the overlays readable.  A
    companion figure ``circadian_null_model_stacked`` is written with
    the same variables on rows sharing a single hour-of-day x-axis.
    """
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde
    setup_figure()

    path = resolve_input(out_dir, "circadian_null_model.csv")
    if not path.exists():
        log.info("  circadian_null_model.csv not found, skipping")
        return

    df = pd.read_csv(path)
    if df.empty:
        return

    log.info("  Plotting circadian null model (2×2) …")

    hourly = {
        "body_temp": _hourly_by_day_type(df, "body_temp_mean"),
        "thi":       _hourly_by_day_type(df, "thi_mean"),
        "barn_temp": _hourly_by_day_type(df, "barn_temp_mean"),
    }

    # ── Crossing density KDE per predictor ─────────────────
    kdes: dict[str, tuple] = {}
    crossing_path = resolve_input(out_dir, "crossing_times.csv")
    if crossing_path.exists():
        ct = pd.read_csv(crossing_path)
        for pred in ("thi", "temp"):
            cpred = ct[ct["predictor"] == pred] if "predictor" in ct.columns else ct
            vals = cpred["day_fraction"].dropna().values
            if len(vals) > 20:
                kde = gaussian_kde(vals, bw_method=0.3)
                x_k = np.linspace(0, 24, 200)
                kdes[pred] = (x_k, kde(x_k), len(vals))

    # ── 2×2 figure ─────────────────────────────────────────
    fig, ((ax_a, ax_c), (ax_b, ax_d)) = plt.subplots(
        2, 2, figsize=(16, 10))

    # (A) Rumen temperature
    _draw_cool_stress_panel(ax_a, hourly["body_temp"])
    ax_a2 = _overlay_crossing_density(ax_a, kdes.get("thi"),
                                      "THI crossing density")
    ax_a.set_ylabel("Rumen temperature (°C)")
    ax_a.set_title("(A) Rumen temperature circadian profile")
    _merge_legends(ax_a, ax_a2, fontsize=8, loc="upper left")

    # (B) Stress − cool Δ rumen
    cool_h = hourly["body_temp"].get("cool")
    stress_h = hourly["body_temp"].get("stress")
    if cool_h is not None and stress_h is not None:
        shared = cool_h.index.intersection(stress_h.index)
        if len(shared) > 5:
            diff = stress_h.loc[shared, "mean"] - cool_h.loc[shared, "mean"]
            dsem = np.sqrt(cool_h.loc[shared, "sem"] ** 2
                           + stress_h.loc[shared, "sem"] ** 2)
            ax_b.fill_between(shared, diff - dsem, diff + dsem,
                              alpha=0.2, color="#CC79A7")
            ax_b.plot(shared, diff, color="#CC79A7", linewidth=2.5,
                      marker="s", markersize=5,
                      label="Stress − cool (Δ rumen temp)")
            ax_b.axhline(0, color="#999", linewidth=0.8, linestyle="--")
            for s, e in _MILKING_WINDOWS:
                ax_b.axvspan(s, e, alpha=0.08, color="#999", zorder=0)
    ax_b2 = _overlay_crossing_density(ax_b, kdes.get("thi"),
                                      "THI crossing density")
    ax_b.set_ylabel("Δ rumen temperature (°C)")
    ax_b.set_title("(B) Heat stress effect (stress − cool)")
    _merge_legends(ax_b, ax_b2, fontsize=8, loc="upper left")

    # (C) Barn THI
    _draw_cool_stress_panel(ax_c, hourly["thi"])
    ax_c2 = _overlay_crossing_density(ax_c, kdes.get("thi"),
                                      "THI crossing density")
    ax_c.set_ylabel("Barn THI")
    ax_c.set_title("(C) Barn THI circadian profile")
    _merge_legends(ax_c, ax_c2, fontsize=8, loc="upper left")

    # (D) Barn temperature
    if hourly["barn_temp"]:
        _draw_cool_stress_panel(ax_d, hourly["barn_temp"])
        ax_d2 = _overlay_crossing_density(ax_d, kdes.get("temp"),
                                          "Barn-temp crossing density")
        ax_d.set_ylabel("Barn temperature (°C)")
        ax_d.set_title("(D) Barn temperature circadian profile")
        _merge_legends(ax_d, ax_d2, fontsize=8, loc="upper left")
    else:
        ax_d.text(0.5, 0.5,
                  "barn_temp_mean absent from CSV —\n"
                  "re-run digimuh-stats to populate",
                  transform=ax_d.transAxes, ha="center", va="center",
                  fontsize=10, color="#666")
        ax_d.set_axis_off()

    for ax in (ax_a, ax_b, ax_c, ax_d):
        ax.set_xlabel("Hour of day")
        ax.set_xticks(range(0, 24, 2))
        ax.set_xlim(-0.5, 23.5)
        ax.grid(False)

    fig.suptitle("Circadian null model — cool vs heat-stress days",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "circadian_null_model", out_dir)

    # ── Companion stacked figure: single shared x-axis ─────
    _plot_circadian_stacked(out_dir, hourly, kdes)


def _plot_circadian_stacked(
    out_dir: Path,
    hourly: dict[str, dict[str, pd.DataFrame]],
    kdes: dict[str, tuple],
) -> None:
    """Stacked companion plot: one shared hour-of-day x-axis, each
    variable on its own row, crossing densities in the bottom row."""
    import matplotlib.pyplot as plt

    rows: list[tuple[str, str]] = []
    if hourly["body_temp"]:
        rows.append(("body_temp", "Rumen temp (°C)"))
    if hourly["thi"]:
        rows.append(("thi", "Barn THI"))
    if hourly["barn_temp"]:
        rows.append(("barn_temp", "Barn temp (°C)"))
    if not rows:
        return

    include_kde_row = bool(kdes)
    n_rows = len(rows) + (1 if include_kde_row else 0)

    fig, axes = plt.subplots(n_rows, 1,
                             figsize=(10, 2.2 * n_rows),
                             sharex=True)
    if n_rows == 1:
        axes = [axes]

    for ax, (var_key, ylabel) in zip(axes, rows):
        _draw_cool_stress_panel(ax, hourly[var_key])
        ax.set_ylabel(ylabel)
        ax.grid(False)
        ax.legend(fontsize=7, loc="upper left")

    if include_kde_row:
        ax_k = axes[-1]
        for pred, colour in [("thi", _KDE_COLOUR), ("temp", "#E69F00")]:
            if pred not in kdes:
                continue
            x_k, y_k, n_k = kdes[pred]
            ax_k.fill_between(x_k, 0, y_k, alpha=0.22, color=colour)
            ax_k.plot(x_k, y_k, color=colour, linewidth=1.5,
                      label=f"{pred.upper()} crossings (n={n_k})")
        for s, e in _MILKING_WINDOWS:
            ax_k.axvspan(s, e, alpha=0.08, color="#999", zorder=0)
        ax_k.set_ylabel("Crossing density")
        ax_k.grid(False)
        ax_k.legend(fontsize=7, loc="upper left")

    axes[-1].set_xlabel("Hour of day")
    axes[-1].set_xticks(range(0, 24, 2))
    axes[-1].set_xlim(-0.5, 23.5)
    fig.suptitle("Circadian null model — stacked view",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "circadian_null_model_stacked", out_dir)


# ─────────────────────────────────────────────────────────────
#  « THI daily exceedance profile »
# ─────────────────────────────────────────────────────────────

def plot_thi_daily_profile(out_dir: Path) -> None:
    """Plot barn THI across 24h by month, with herd breakpoint line."""
    import matplotlib.pyplot as plt
    setup_figure()

    path = resolve_input(out_dir, "thi_daily_profile.csv")
    if not path.exists():
        log.info("  thi_daily_profile.csv not found, skipping")
        return

    df = pd.read_csv(path)
    if df.empty:
        return

    log.info("  Plotting THI daily profile …")

    herd_bp = df["herd_median_bp"].iloc[0] if "herd_median_bp" in df.columns else np.nan

    # One plot per year
    years = sorted(df["year"].unique().astype(int))

    for year in years:
        ydf = df[df["year"] == year]
        if ydf.empty:
            continue

        fig, ax = plt.subplots(figsize=(10, 6))

        month_colours = {6: "#009E73", 7: "#E69F00", 8: "#D55E00", 9: "#CC79A7"}
        month_names = {6: "June", 7: "July", 8: "August", 9: "September"}

        for month in sorted(ydf["month"].unique().astype(int)):
            msub = ydf[ydf["month"] == month]
            if msub.empty:
                continue

            colour = month_colours.get(month, "#888")
            ax.fill_between(msub["hour"], msub["thi_q25"], msub["thi_q75"],
                           alpha=0.1, color=colour)
            ax.plot(msub["hour"], msub["thi_mean"], color=colour, linewidth=2,
                    marker="o", markersize=3,
                    label=f"{month_names.get(month, str(month))}")

        # Herd median breakpoint
        if not np.isnan(herd_bp):
            ax.axhline(herd_bp, color="#333", linewidth=1.5, linestyle="--",
                       label=f"Herd median THI breakpoint ({herd_bp:.1f})")

        # Mark milking windows
        for start, end in [(4, 7), (16, 19)]:
            ax.axvspan(start, end, alpha=0.08, color="#999", zorder=0)

        ax.set_xlabel("Hour of day")
        ax.set_ylabel("Barn THI")
        ax.set_title(f"Barn THI daily profile ({year})\n"
                     f"(lines = mean, shading = IQR)")
        ax.set_xticks(range(0, 24, 2))
        ax.set_xlim(-0.5, 23.5)
        ax.legend(fontsize=9)
        fig.tight_layout()
        save_figure(fig, f"thi_daily_profile_{year}", out_dir)

    # All years combined
    fig, ax = plt.subplots(figsize=(12, 6))
    labels_seen = set()
    for ml in sorted(df["month_label"].unique()):
        msub = df[df["month_label"] == ml]
        month = int(msub["month"].iloc[0])
        year = int(msub["year"].iloc[0])
        colour = month_colours.get(month, "#888")
        # Vary line style by year
        ls = ["-", "--", ":", "-."][years.index(year) % 4]
        ax.plot(msub["hour"], msub["thi_mean"], color=colour, linewidth=1.5,
                linestyle=ls, alpha=0.7, label=ml)

    if not np.isnan(herd_bp):
        ax.axhline(herd_bp, color="#333", linewidth=1.5, linestyle="--",
                   label=f"Herd breakpoint ({herd_bp:.1f})")

    for start, end in [(4, 7), (16, 19)]:
        ax.axvspan(start, end, alpha=0.08, color="#999", zorder=0)

    ax.set_xlabel("Hour of day")
    ax.set_ylabel("Barn THI")
    ax.set_title("Barn THI daily profile — all years\n"
                 "(when does heat stress occur?)")
    ax.set_xticks(range(0, 24, 2))
    ax.set_xlim(-0.5, 23.5)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    save_figure(fig, "thi_daily_profile_all", out_dir)


# ─────────────────────────────────────────────────────────────
#  « breakpoint crossing raster + raincloud »
# ─────────────────────────────────────────────────────────────

def plot_crossing_raster(out_dir: Path) -> None:
    """Raster plot of breakpoint crossing events across the 24h cycle.

    Left panel: activation raster — each row is one cow (sorted by
    breakpoint value), each dot is a crossing event at that clock time.
    Dot colour = breakpoint value.

    Right panel: raincloud — half-violin (KDE) + jittered scatter +
    boxplot of crossing clock times.
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import PathCollection
    setup_figure()

    path = resolve_input(out_dir, "crossing_times.csv")
    if not path.exists():
        log.info("  crossing_times.csv not found, skipping raster plot")
        return

    df = pd.read_csv(path)
    if df.empty:
        return

    log.info("  Plotting crossing raster (%d events) …", len(df))

    for pred, pred_label, bp_label in [
        ("thi", "THI breakpoint crossing", "THI breakpoint"),
        ("temp", "Barn temp breakpoint crossing", "Barn temp breakpoint (°C)"),
    ]:
        sub = df[df["predictor"] == pred]
        if sub.empty or sub["animal_id"].nunique() < 5:
            continue

        # Sort animals by their breakpoint value
        animal_bps = sub.groupby("animal_id")["breakpoint"].first().sort_values()
        animal_order = {aid: i for i, aid in enumerate(animal_bps.index)}
        sub = sub.copy()
        sub["y_pos"] = sub["animal_id"].map(animal_order)
        n_animals = len(animal_order)

        fig, axes = plt.subplots(1, 2, figsize=(16, max(6, n_animals * 0.06)),
                                 gridspec_kw={"width_ratios": [3, 1]})

        # ── Left: Raster ─────────────────────────────────────
        ax = axes[0]
        sc = ax.scatter(
            sub["day_fraction"], sub["y_pos"],
            c=sub["breakpoint"], cmap="coolwarm",
            s=3, alpha=0.3, edgecolors="none",
            vmin=sub["breakpoint"].quantile(0.05),
            vmax=sub["breakpoint"].quantile(0.95),
        )
        cbar = fig.colorbar(sc, ax=ax, pad=0.01, shrink=0.7)
        cbar.set_label(bp_label, fontsize=9)

        # Milking windows
        for start, end in [(4, 7), (16, 19)]:
            ax.axvspan(start, end, alpha=0.08, color="#999", zorder=0)

        ax.set_xlabel("Hour of day")
        ax.set_ylabel(f"Animals (sorted by {bp_label}, n={n_animals})")
        ax.set_xlim(-0.5, 24)
        ax.set_xticks(range(0, 25, 2))
        ax.set_ylim(-1, n_animals)
        ax.set_yticks([])
        ax.set_title(f"Crossing event raster\n({len(sub)} events)")

        # ── Right: Raincloud ─────────────────────────────────
        ax = axes[1]

        # Half-violin (KDE)
        from scipy.stats import gaussian_kde
        vals = sub["day_fraction"].dropna().values
        if len(vals) > 10:
            kde = gaussian_kde(vals, bw_method=0.3)
            x_kde = np.linspace(0, 24, 200)
            y_kde = kde(x_kde)
            # Scale KDE to fit nicely
            y_kde_scaled = y_kde / y_kde.max() * 0.35
            ax.fill_betweenx(x_kde, 0, y_kde_scaled,
                            alpha=0.3, color=COLOURS["below_bp"])
            ax.plot(y_kde_scaled, x_kde, color=COLOURS["below_bp"], linewidth=1.5)

        # Jittered scatter
        jitter = np.random.uniform(0.4, 0.7, len(vals))
        ax.scatter(jitter, vals, s=2, alpha=0.15,
                   color=COLOURS["below_bp"], edgecolors="none")

        # Boxplot
        bp_plot = ax.boxplot(
            vals, positions=[0.85], widths=0.12,
            vert=True, patch_artist=True,
            boxprops=dict(facecolor=COLOURS["below_bp"], alpha=0.5, edgecolor="#333"),
            medianprops=dict(color="#333", linewidth=2),
            whiskerprops=dict(color="#333"),
            capprops=dict(color="#333"),
            flierprops=dict(marker="o", markersize=2, alpha=0.3),
        )

        ax.set_ylim(-0.5, 24)
        ax.set_yticks(range(0, 25, 2))
        ax.set_ylabel("Hour of day")
        ax.set_xlim(-0.05, 1.1)
        ax.set_xticks([])
        ax.set_title("Distribution")
        ax.invert_yaxis()  # 0h at top to match intuition

        # Milking windows on raincloud too
        for start, end in [(4, 7), (16, 19)]:
            ax.axhspan(start, end, alpha=0.08, color="#999", zorder=0)

        fig.suptitle(f"Breakpoint crossing activation: {pred_label}",
                     fontsize=13, fontweight="bold")
        fig.tight_layout()
        save_figure(fig, f"crossing_raster_{pred}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « cross-correlation below/above breakpoint »
# ─────────────────────────────────────────────────────────────


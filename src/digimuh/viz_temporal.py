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

log = logging.getLogger("digimuh.viz")

# ─────────────────────────────────────────────────────────────
#  « rumen circadian null model »
# ─────────────────────────────────────────────────────────────

def plot_circadian_null_model(out_dir: Path) -> None:
    """Plot 24h rumen temperature profile: cool days vs stress days,
    with breakpoint crossing density and stress-cool difference curve.

    Two panels:
    A) Cool-day and stress-day profiles + crossing density KDE
    B) Difference curve (stress minus cool) = additional rumen temperature
       attributable to heat stress, with crossing density for reference
    """
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde
    setup_figure()

    path = out_dir / "circadian_null_model.csv"
    if not path.exists():
        log.info("  circadian_null_model.csv not found, skipping")
        return

    df = pd.read_csv(path)
    if df.empty:
        return

    log.info("  Plotting circadian null model …")

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(10, 10),
                                          gridspec_kw={"height_ratios": [1.2, 1]})

    # ── Compute hourly profiles ──────────────────────────────
    hourly_data = {}
    for day_type in ["cool", "stress"]:
        sub = df[df["day_type"] == day_type]
        if sub.empty:
            continue
        hourly = sub.groupby("hour")["body_temp_mean"].agg(["mean", "std", "count"])
        hourly["sem"] = hourly["std"] / np.sqrt(hourly["count"])
        hourly_data[day_type] = hourly

    # Load crossing density
    crossing_path = out_dir / "crossing_times.csv"
    has_crossings = False
    if crossing_path.exists():
        ct = pd.read_csv(crossing_path)
        ct_thi = ct[ct["predictor"] == "thi"] if "predictor" in ct.columns else ct
        if len(ct_thi) > 20:
            has_crossings = True
            vals = ct_thi["day_fraction"].dropna().values
            kde = gaussian_kde(vals, bw_method=0.3)
            x_kde = np.linspace(0, 24, 200)
            y_kde = kde(x_kde)

    # ── Panel A: Cool vs stress profiles ─────────────────────
    for day_type, colour, label in [
        ("cool", COLOURS["below_bp"], "Cool days (THI < breakpoint all day)"),
        ("stress", COLOURS["above_bp"], "Heat stress days (THI exceeded breakpoint)"),
    ]:
        if day_type not in hourly_data:
            continue
        hourly = hourly_data[day_type]
        ax_top.fill_between(hourly.index, hourly["mean"] - hourly["sem"],
                           hourly["mean"] + hourly["sem"],
                           alpha=0.2, color=colour)
        ax_top.plot(hourly.index, hourly["mean"], color=colour, linewidth=2,
                    marker="o", markersize=4, label=label)

    for start, end in [(4, 7), (16, 19)]:
        ax_top.axvspan(start, end, alpha=0.08, color="#999", zorder=0)

    ax_top.set_ylabel("Rumen temperature (°C)")
    ax_top.set_xticks(range(0, 24, 2))
    ax_top.set_xlim(-0.5, 23.5)

    if has_crossings:
        ax_top2 = ax_top.twinx()
        ax_top2.fill_between(x_kde, 0, y_kde, alpha=0.12,
                            color="#009E73", zorder=0)
        ax_top2.plot(x_kde, y_kde, color="#009E73", linewidth=1.5,
                     alpha=0.7, label=f"Crossing density (n={len(vals)})")
        ax_top2.set_ylabel("Crossing density", color="#009E73")
        ax_top2.tick_params(axis="y", labelcolor="#009E73")
        ax_top2.set_ylim(0, y_kde.max() * 2.5)

        lines1, labels1 = ax_top.get_legend_handles_labels()
        lines2, labels2 = ax_top2.get_legend_handles_labels()
        ax_top.legend(lines1 + lines2, labels1 + labels2,
                      fontsize=9, loc="upper left")
    else:
        ax_top.legend(fontsize=9, loc="upper left")

    ax_top.set_title("(A) Rumen temperature circadian profile")

    # ── Panel B: Difference curve (stress minus cool) ────────
    if "cool" in hourly_data and "stress" in hourly_data:
        cool_h = hourly_data["cool"]
        stress_h = hourly_data["stress"]
        shared_hours = cool_h.index.intersection(stress_h.index)

        if len(shared_hours) > 5:
            diff_mean = stress_h.loc[shared_hours, "mean"] - cool_h.loc[shared_hours, "mean"]
            # Propagate SEM: sqrt(sem_cool² + sem_stress²)
            diff_sem = np.sqrt(
                cool_h.loc[shared_hours, "sem"] ** 2 +
                stress_h.loc[shared_hours, "sem"] ** 2
            )

            ax_bot.fill_between(shared_hours, diff_mean - diff_sem,
                               diff_mean + diff_sem,
                               alpha=0.2, color="#CC79A7")
            ax_bot.plot(shared_hours, diff_mean, color="#CC79A7", linewidth=2.5,
                        marker="s", markersize=5,
                        label="Stress − cool (additional rumen temp)")
            ax_bot.axhline(0, color="#999", linewidth=0.8, linestyle="--")

            for start, end in [(4, 7), (16, 19)]:
                ax_bot.axvspan(start, end, alpha=0.08, color="#999", zorder=0)

            if has_crossings:
                ax_bot2 = ax_bot.twinx()
                ax_bot2.fill_between(x_kde, 0, y_kde, alpha=0.10,
                                    color="#009E73", zorder=0)
                ax_bot2.plot(x_kde, y_kde, color="#009E73", linewidth=1.2,
                             alpha=0.5, label="Crossing density")
                ax_bot2.set_ylabel("Crossing density", color="#009E73")
                ax_bot2.tick_params(axis="y", labelcolor="#009E73")
                ax_bot2.set_ylim(0, y_kde.max() * 2.5)

                lines1, labels1 = ax_bot.get_legend_handles_labels()
                lines2, labels2 = ax_bot2.get_legend_handles_labels()
                ax_bot.legend(lines1 + lines2, labels1 + labels2,
                              fontsize=9, loc="upper left")
            else:
                ax_bot.legend(fontsize=9, loc="upper left")
    else:
        ax_bot.text(0.5, 0.5, "Need both cool and stress day data",
                    transform=ax_bot.transAxes, ha="center")

    ax_bot.set_xlabel("Hour of day")
    ax_bot.set_ylabel("Additional rumen temperature (°C)")
    ax_bot.set_xticks(range(0, 24, 2))
    ax_bot.set_xlim(-0.5, 23.5)
    ax_bot.set_title("(B) Heat stress effect (stress − cool day profile)")

    fig.suptitle("Rumen temperature circadian null model",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "circadian_null_model", out_dir)


# ─────────────────────────────────────────────────────────────
#  « THI daily exceedance profile »
# ─────────────────────────────────────────────────────────────

def plot_thi_daily_profile(out_dir: Path) -> None:
    """Plot barn THI across 24h by month, with herd breakpoint line."""
    import matplotlib.pyplot as plt
    setup_figure()

    path = out_dir / "thi_daily_profile.csv"
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

    path = out_dir / "crossing_times.csv"
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


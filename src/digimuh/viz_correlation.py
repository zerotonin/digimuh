#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — viz_correlation                                      ║
# ║  « cross-correlation, derivative CCF, and event-triggered avg »║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Climate–rumen coupling figures: raw and detrended cross-       ║
# ║  correlation, derivative CCF, event-triggered average (peri-    ║
# ║  event), and climate ETA at breakpoint crossings.               ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Correlation and event-triggered figures for the Frontiers manuscript."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.constants import COLOURS
from digimuh.viz_base import setup_figure, save_figure

log = logging.getLogger("digimuh.viz")

# ─────────────────────────────────────────────────────────────
#  « cross-correlation below/above breakpoint »
# ─────────────────────────────────────────────────────────────

def plot_cross_correlation(out_dir: Path) -> None:
    """Plot cross-correlation AND cross-covariance curves below vs above bp.

    Error bands are standard error of the mean (SEM).
    Peak lag annotated with vertical coloured line and horizontal arrow
    from lag=0 to peak, labelled 'rumen temperature N min later'.
    """
    import matplotlib.pyplot as plt
    setup_figure()

    xcorr_path = out_dir / "cross_correlation.csv"
    if not xcorr_path.exists():
        log.info("  cross_correlation.csv not found, skipping xcorr plots")
        return

    xcorr = pd.read_csv(xcorr_path)
    if xcorr.empty:
        log.info("  cross_correlation.csv is empty, skipping xcorr plots")
        return

    log.info("  Plotting cross-correlations (%d rows) …", len(xcorr))

    def _annotate_peak(ax, agg, colour, y_offset_frac=0.08):
        """Add vertical peak line + horizontal arrow from lag=0."""
        peak_lag = agg["mean"].idxmax()
        peak_val = agg["mean"].max()
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        arrow_y = peak_val + y_range * y_offset_frac

        if abs(peak_lag) > 5:  # only annotate if peak is away from zero
            # Vertical dashed line at peak
            ax.axvline(peak_lag, color=colour, linewidth=1.2,
                       linestyle="--", alpha=0.7)

            # Horizontal arrow from lag=0 to peak
            ax.annotate(
                "", xy=(peak_lag, arrow_y), xytext=(0, arrow_y),
                arrowprops=dict(arrowstyle="<->", color=colour, lw=1.5),
            )
            # Label above the arrow
            mid_x = peak_lag / 2
            ax.text(mid_x, arrow_y + y_range * 0.02,
                    f"rumen temperature\n{abs(peak_lag):.0f} min later",
                    ha="center", va="bottom", fontsize=8,
                    color=colour, fontstyle="italic")

    # Detect whether we have the variant column (raw vs detrended)
    has_variants = "variant" in xcorr.columns
    variants = [("raw", ""), ("detrended", " (diurnal removed)")] if has_variants else [("raw", "")]

    for variant, variant_subtitle in variants:
        vdata = xcorr[xcorr["variant"] == variant] if has_variants else xcorr
        if vdata.empty:
            continue

        variant_suffix = f"_{variant}" if variant != "raw" else ""

        # Plot both cross-correlation and cross-covariance
        for metric, metric_label, ylab, fname_suffix in [
            ("xcorr", "Cross-correlation", "Cross-correlation (r)", "xcorr"),
            ("xcov", "Cross-covariance", "Cross-covariance", "xcov"),
        ]:
            # ── Separate panels per region ────────────────────
            for pred, pred_label in [("thi", "Barn THI"), ("temp", "Barn temperature")]:
                sub = vdata[vdata["predictor"] == pred]
                if sub.empty:
                    continue

                fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

                for ax, (region, colour, label) in zip(axes, [
                    ("below", COLOURS["below_bp"], "Below breakpoint"),
                    ("above", COLOURS["above_bp"], "Above breakpoint"),
                ]):
                    rsub = sub[sub["region"] == region]
                    if rsub.empty:
                        continue

                    agg = rsub.groupby("lag_minutes")[metric].agg(["mean", "std", "count"])
                    agg["sem"] = agg["std"] / np.sqrt(agg["count"])

                    ax.fill_between(agg.index, agg["mean"] - agg["sem"],
                                   agg["mean"] + agg["sem"],
                                   alpha=0.2, color=colour, label="SEM")
                    ax.plot(agg.index, agg["mean"], color=colour, linewidth=2)
                    ax.axhline(0, color="#999", linewidth=0.5, linestyle="--")
                    ax.axvline(0, color="#999", linewidth=1, linestyle="--",
                               label="Lag = 0")

                    _annotate_peak(ax, agg, colour)

                    ax.set_xlabel("Lag (minutes)")
                    ax.set_ylabel(ylab)
                    ax.set_title(f"{label}\n(n={rsub['animal_id'].nunique()} animals)")

                fig.suptitle(
                    f"{metric_label}: {pred_label} vs rumen temperature{variant_subtitle}\n"
                    f"(shaded: ± 1 SEM across animals)",
                    fontsize=13, fontweight="bold")
                fig.tight_layout()
                save_figure(fig, f"{fname_suffix}_{pred}{variant_suffix}", out_dir)

            # ── Overlay: below vs above on same axes ─────────
            for pred, pred_label in [("thi", "Barn THI"), ("temp", "Barn temperature")]:
                sub = vdata[vdata["predictor"] == pred]
                if sub.empty:
                    continue

                fig, ax = plt.subplots(figsize=(10, 6))
                y_offsets = {"below": 0.05, "above": 0.12}

                for region, colour, label in [
                    ("below", COLOURS["below_bp"], "Below breakpoint"),
                    ("above", COLOURS["above_bp"], "Above breakpoint"),
                ]:
                    rsub = sub[sub["region"] == region]
                    if rsub.empty:
                        continue
                    agg = rsub.groupby("lag_minutes")[metric].agg(["mean", "std", "count"])
                    agg["sem"] = agg["std"] / np.sqrt(agg["count"])

                    ax.fill_between(agg.index, agg["mean"] - agg["sem"],
                                   agg["mean"] + agg["sem"],
                                   alpha=0.15, color=colour)
                    ax.plot(agg.index, agg["mean"], color=colour, linewidth=2,
                            label=label)

                    _annotate_peak(ax, agg, colour, y_offset_frac=y_offsets[region])

                ax.axhline(0, color="#999", linewidth=0.5, linestyle="--")
                ax.axvline(0, color="#999", linewidth=1, linestyle="--",
                           label="Lag = 0")
                ax.set_xlabel("Lag (minutes)")
                ax.set_ylabel(ylab)
                ax.set_title(
                    f"{metric_label}: {pred_label} vs rumen temperature{variant_subtitle}\n"
                    f"(shaded: ± 1 SEM)")
                ax.legend()
                fig.tight_layout()
                save_figure(fig, f"{fname_suffix}_{pred}_overlay{variant_suffix}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « derivative cross-correlation plots »
# ─────────────────────────────────────────────────────────────

def plot_derivative_ccf(out_dir: Path) -> None:
    """Plot derivative CCF: d(climate)/dt vs d(body_temp)/dt."""
    import matplotlib.pyplot as plt
    setup_figure()

    dccf_path = out_dir / "derivative_ccf.csv"
    if not dccf_path.exists():
        log.info("  derivative_ccf.csv not found, skipping")
        return

    dccf = pd.read_csv(dccf_path)
    if dccf.empty:
        return

    log.info("  Plotting derivative CCF (%d rows) …", len(dccf))

    def _annotate_peak(ax, agg, colour, y_offset_frac=0.08):
        """Add vertical peak line + horizontal arrow from lag=0."""
        peak_lag = agg["mean"].idxmax()
        peak_val = agg["mean"].max()
        ylo, yhi = ax.get_ylim()
        y_range = yhi - ylo
        arrow_y = peak_val + y_range * y_offset_frac

        if abs(peak_lag) > 5:
            ax.axvline(peak_lag, color=colour, linewidth=1.2,
                       linestyle="--", alpha=0.7)
            ax.annotate(
                "", xy=(peak_lag, arrow_y), xytext=(0, arrow_y),
                arrowprops=dict(arrowstyle="<->", color=colour, lw=1.5),
            )
            mid_x = peak_lag / 2
            ax.text(mid_x, arrow_y + y_range * 0.02,
                    f"rumen temperature\n{abs(peak_lag):.0f} min later",
                    ha="center", va="bottom", fontsize=8,
                    color=colour, fontstyle="italic")

    for pred, pred_label in [("thi", "Barn THI"), ("temp", "Barn temperature")]:
        sub = dccf[dccf["predictor"] == pred]
        if sub.empty:
            continue

        fig, ax = plt.subplots(figsize=(10, 6))
        y_offsets = {"below": 0.05, "above": 0.12}

        for region, colour, label in [
            ("below", COLOURS["below_bp"], "Below breakpoint"),
            ("above", COLOURS["above_bp"], "Above breakpoint"),
        ]:
            rsub = sub[sub["region"] == region]
            if rsub.empty:
                continue
            agg = rsub.groupby("lag_minutes")["dxcorr"].agg(["mean", "std", "count"])
            agg["sem"] = agg["std"] / np.sqrt(agg["count"])

            ax.fill_between(agg.index, agg["mean"] - agg["sem"],
                           agg["mean"] + agg["sem"],
                           alpha=0.15, color=colour)
            ax.plot(agg.index, agg["mean"], color=colour, linewidth=2,
                    label=label)

            _annotate_peak(ax, agg, colour, y_offset_frac=y_offsets[region])

        ax.axhline(0, color="#999", linewidth=0.5, linestyle="--")
        ax.axvline(0, color="#999", linewidth=1, linestyle="--",
                   label="Lag = 0")
        ax.set_xlabel("Lag (minutes)")
        ax.set_ylabel("Cross-correlation of derivatives (r)")
        ax.set_title(
            f"Derivative CCF: d({pred_label})/dt vs d(Rumen temp)/dt\n"
            f"(shaded: ± 1 SEM)")
        ax.legend()
        fig.tight_layout()
        save_figure(fig, f"dccf_{pred}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « event-triggered average plots »
# ─────────────────────────────────────────────────────────────

def plot_event_triggered_average(
    out_dir: Path,
    traces_file: str = "event_triggered_traces.csv",
    suffix: str = "",
    title_extra: str = "",
) -> None:
    """Plot peri-event average of rumen temp around breakpoint crossings.

    Three panels per predictor:
    A) Climate signal (THI or barn temp) aligned to crossing
    B) Raw rumen temperature
    C) Rumen temperature baseline-subtracted (acute onset)

    Args:
        out_dir: Output directory.
        traces_file: Name of the traces CSV file.
        suffix: Appended to output filename (e.g. '_filtered').
        title_extra: Appended to figure title (e.g. ' (8-11h crossings)').
    """
    import matplotlib.pyplot as plt
    setup_figure()

    traces_path = out_dir / traces_file
    if not traces_path.exists():
        log.info("  %s not found, skipping", traces_file)
        return

    traces = pd.read_csv(traces_path)
    if traces.empty:
        return

    log.info("  Plotting event-triggered averages (%d trace points) …", len(traces))

    for pred, pred_label, climate_label in [
        ("thi", "THI breakpoint crossing", "Barn THI"),
        ("temp", "Barn temp breakpoint crossing", "Barn temperature (°C)"),
    ]:
        sub = traces[traces["predictor"] == pred]
        if sub.empty:
            continue

        n_events = sub.groupby(["animal_id", "year", "event_id"]).ngroups
        n_animals = sub["animal_id"].nunique()

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # ── Panel A: Climate signal ──────────────────────────
        ax = axes[0, 0]
        agg = sub.groupby("relative_minutes")["climate_val"].agg(["mean", "std", "count"])
        agg["sem"] = agg["std"] / np.sqrt(agg["count"])

        ax.fill_between(agg.index, agg["mean"] - agg["sem"],
                       agg["mean"] + agg["sem"],
                       alpha=0.2, color=COLOURS["above_bp"])
        ax.plot(agg.index, agg["mean"], color=COLOURS["above_bp"], linewidth=2)
        ax.axvline(0, color="#333", linewidth=1.5, linestyle="--",
                   label="Breakpoint crossing")
        ax.set_xlabel("Time relative to crossing (minutes)")
        ax.set_ylabel(climate_label)
        ax.set_title(f"(A) {climate_label}\n(n={n_events} events, {n_animals} animals)")
        ax.legend(fontsize=8)

        # ── Panel B: Raw rumen temperature ──────────────────
        ax = axes[0, 1]
        agg_raw = sub.groupby("relative_minutes")["body_temp"].agg(
            ["mean", "std", "count"])
        agg_raw["sem"] = agg_raw["std"] / np.sqrt(agg_raw["count"])

        ax.fill_between(agg_raw.index, agg_raw["mean"] - agg_raw["sem"],
                       agg_raw["mean"] + agg_raw["sem"],
                       alpha=0.2, color=COLOURS["below_bp"])
        ax.plot(agg_raw.index, agg_raw["mean"], color=COLOURS["below_bp"], linewidth=2)
        ax.axvline(0, color="#333", linewidth=1.5, linestyle="--",
                   label="Breakpoint crossing")
        ax.set_xlabel("Time relative to crossing (minutes)")
        ax.set_ylabel("Rumen temperature (°C)")
        ax.set_title("(B) Rumen temperature")
        ax.legend(fontsize=8)

        # ── Panel C: Rumen temp baseline-subtracted ──────────
        ax = axes[1, 0]
        agg_bl = sub.groupby("relative_minutes")["body_temp_baseline"].agg(
            ["mean", "std", "count"])
        agg_bl["sem"] = agg_bl["std"] / np.sqrt(agg_bl["count"])

        ax.fill_between(agg_bl.index, agg_bl["mean"] - agg_bl["sem"],
                       agg_bl["mean"] + agg_bl["sem"],
                       alpha=0.2, color=COLOURS["below_bp"])
        ax.plot(agg_bl.index, agg_bl["mean"], color=COLOURS["below_bp"], linewidth=2)
        ax.axvline(0, color="#333", linewidth=1.5, linestyle="--",
                   label="Breakpoint crossing")
        ax.axhline(0, color="#999", linewidth=0.5, linestyle="--")

        # Find when body temp first exceeds 2 SEM above baseline
        post = agg_bl.loc[agg_bl.index >= 0]
        onset = post[post["mean"] > 2 * post["sem"].iloc[0]].index
        if len(onset) > 0:
            onset_min = onset[0]
            onset_val = agg_bl.loc[onset_min, "mean"]
            ax.axvline(onset_min, color=COLOURS["below_bp"], linewidth=1,
                       linestyle=":", alpha=0.7)
            ax.annotate(
                f"Response onset\n{onset_min:.0f} min",
                xy=(onset_min, onset_val),
                xytext=(onset_min + 40, onset_val + 0.02),
                fontsize=9, color=COLOURS["below_bp"],
                arrowprops=dict(arrowstyle="->", color=COLOURS["below_bp"], lw=1),
            )

        ax.set_xlabel("Time relative to crossing (minutes)")
        ax.set_ylabel("Rumen temp change from pre-event (°C)")
        ax.set_title("(C) Rumen temperature\n(pre-event baseline subtracted)")
        ax.legend(fontsize=8)

        # ── Panel D: Additional rumen temperature ────────────
        #    raw body_temp minus cool-day circadian profile at that clock hour
        ax = axes[1, 1]
        circadian_path = out_dir / "circadian_null_model.csv"
        has_clock = "clock_hour" in sub.columns
        has_circ = circadian_path.exists()

        if has_circ and has_clock:
            circ = pd.read_csv(circadian_path)
            cool = circ[circ["day_type"] == "cool"]

            if not cool.empty:
                # Per-animal hourly cool-day lookup
                cool_lookup = {}
                for _, crow in cool.iterrows():
                    cool_lookup[(int(crow["animal_id"]), int(crow["year"]),
                                 int(crow["hour"]))] = crow["body_temp_mean"]
                # Grand mean fallback per hour
                grand_cool = cool.groupby("hour")["body_temp_mean"].mean().to_dict()

                sub_d = sub.copy()

                def _get_cool(row):
                    key = (int(row["animal_id"]), int(row["year"]),
                           int(row["clock_hour"]))
                    v = cool_lookup.get(key)
                    if v is None:
                        v = grand_cool.get(int(row["clock_hour"]))
                    return v

                cool_vals = sub_d.apply(_get_cool, axis=1)
                sub_d["additional_temp"] = sub_d["body_temp"] - cool_vals
                valid = sub_d.dropna(subset=["additional_temp"])

                if len(valid) > 100:
                    agg_add = valid.groupby("relative_minutes")["additional_temp"].agg(
                        ["mean", "std", "count"])
                    agg_add["sem"] = agg_add["std"] / np.sqrt(agg_add["count"])

                    ax.fill_between(agg_add.index,
                                   agg_add["mean"] - agg_add["sem"],
                                   agg_add["mean"] + agg_add["sem"],
                                   alpha=0.2, color=COLOURS["above_bp"])
                    ax.plot(agg_add.index, agg_add["mean"],
                            color=COLOURS["above_bp"], linewidth=2)
                    ax.axvline(0, color="#333", linewidth=1.5, linestyle="--",
                               label="Breakpoint crossing")
                    ax.axhline(0, color="#999", linewidth=0.5, linestyle="--",
                               label="Cool-day circadian level")

                    # Onset detection
                    post = agg_add.loc[agg_add.index >= 0]
                    if len(post) > 2:
                        onset = post[post["mean"] > 2 * post["sem"].iloc[0]].index
                        if len(onset) > 0:
                            onset_min = onset[0]
                            onset_val = agg_add.loc[onset_min, "mean"]
                            ax.axvline(onset_min, color=COLOURS["above_bp"],
                                       linewidth=1, linestyle=":", alpha=0.7)
                            ax.annotate(
                                f"Response onset\n{onset_min:.0f} min",
                                xy=(onset_min, onset_val),
                                xytext=(onset_min + 40, onset_val + 0.01),
                                fontsize=9, color=COLOURS["above_bp"],
                                arrowprops=dict(arrowstyle="->",
                                               color=COLOURS["above_bp"], lw=1),
                            )
                else:
                    ax.text(0.5, 0.5, "Insufficient cool-day data",
                            transform=ax.transAxes, ha="center")
            else:
                ax.text(0.5, 0.5, "No cool-day data available",
                        transform=ax.transAxes, ha="center")
        else:
            missing = []
            if not has_circ:
                missing.append("circadian_null_model.csv")
            if not has_clock:
                missing.append("clock_hour column")
            ax.text(0.5, 0.5, f"Missing: {', '.join(missing)}\nRe-run stats.",
                    transform=ax.transAxes, ha="center", fontsize=9)

        ax.set_xlabel("Time relative to crossing (minutes)")
        ax.set_ylabel("Additional rumen temperature (°C)")
        ax.set_title("(D) Additional rumen temperature\n"
                     "(surplus above cool-day circadian profile)")
        ax.legend(fontsize=8)

        fig.suptitle(f"Event-triggered average: {pred_label}{title_extra}",
                     fontsize=13, fontweight="bold")
        fig.tight_layout()
        save_figure(fig, f"eta_{pred}{suffix}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « climate ETA: normalised THI + barn temp around crossings »
# ─────────────────────────────────────────────────────────────

def plot_climate_eta(out_dir: Path) -> None:
    """Climate signal around breakpoint crossings, normalised to breakpoint.

    Two figures:
    1. THI-triggered crossings: left y = THI − THI_breakpoint,
       right y = barn temperature (raw °C)
    2. Barn temp-triggered crossings: left y = barn_temp − temp_breakpoint,
       right y = THI (raw)

    In both cases y=0 on the left axis is the breakpoint threshold.
    """
    import matplotlib.pyplot as plt
    setup_figure()

    path = out_dir / "climate_eta.csv"
    if not path.exists():
        log.info("  climate_eta.csv not found, skipping")
        return

    df = pd.read_csv(path)
    if df.empty:
        return

    log.info("  Plotting climate ETA (%d trace points) …", len(df))

    configs = [
        {
            "trigger": "thi",
            "norm_col": "thi_norm",
            "companion_col": "temp_norm",
            "norm_label": "ΔTHI (THI − breakpoint)",
            "companion_label": "ΔBarn temperature (°C from breakpoint)",
            "title": "Climate change around THI breakpoint crossing",
            "fname": "climate_eta_thi",
            "norm_colour": COLOURS["above_bp"],
            "comp_colour": "#009E73",
        },
        {
            "trigger": "temp",
            "norm_col": "temp_norm",
            "companion_col": "thi_norm",
            "norm_label": "ΔBarn temperature (°C from breakpoint)",
            "companion_label": "ΔTHI (THI − breakpoint)",
            "title": "Climate change around barn temp breakpoint crossing",
            "fname": "climate_eta_temp",
            "norm_colour": COLOURS["above_bp"],
            "comp_colour": "#009E73",
        },
    ]

    for cfg in configs:
        sub = df[df["trigger"] == cfg["trigger"]]
        if sub.empty:
            continue

        n_events = sub.groupby(["animal_id", "year", "event_id"]).ngroups
        n_animals = sub["animal_id"].nunique()

        fig, ax1 = plt.subplots(figsize=(10, 6))

        # Left axis: normalised trigger predictor
        agg_norm = sub.groupby("relative_minutes")[cfg["norm_col"]].agg(
            ["mean", "std", "count"])
        agg_norm["sem"] = agg_norm["std"] / np.sqrt(agg_norm["count"])

        ax1.fill_between(agg_norm.index,
                        agg_norm["mean"] - agg_norm["sem"],
                        agg_norm["mean"] + agg_norm["sem"],
                        alpha=0.2, color=cfg["norm_colour"])
        ax1.plot(agg_norm.index, agg_norm["mean"],
                 color=cfg["norm_colour"], linewidth=2,
                 label=cfg["norm_label"])
        ax1.axhline(0, color=cfg["norm_colour"], linewidth=1, linestyle="--",
                     alpha=0.5, label="Breakpoint (threshold)")
        ax1.set_xlabel("Time relative to crossing (minutes)")
        ax1.set_ylabel(cfg["norm_label"], color=cfg["norm_colour"])
        ax1.tick_params(axis="y", labelcolor=cfg["norm_colour"])

        # Crossing marker
        ax1.axvline(0, color="#333", linewidth=1.5, linestyle="--")

        # Right axis: companion predictor (normalised to its breakpoint)
        ax2 = ax1.twinx()
        agg_comp = sub.groupby("relative_minutes")[cfg["companion_col"]].agg(
            ["mean", "std", "count"])
        agg_comp["sem"] = agg_comp["std"] / np.sqrt(agg_comp["count"])

        ax2.fill_between(agg_comp.index,
                        agg_comp["mean"] - agg_comp["sem"],
                        agg_comp["mean"] + agg_comp["sem"],
                        alpha=0.12, color=cfg["comp_colour"])
        ax2.plot(agg_comp.index, agg_comp["mean"],
                 color=cfg["comp_colour"], linewidth=2,
                 linestyle="-", label=cfg["companion_label"])
        ax2.axhline(0, color=cfg["comp_colour"], linewidth=1, linestyle="--",
                    alpha=0.4)
        ax2.set_ylabel(cfg["companion_label"], color=cfg["comp_colour"])
        ax2.tick_params(axis="y", labelcolor=cfg["comp_colour"])

        # Combined legend
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2,
                   fontsize=9, loc="upper left")

        ax1.set_title(f"{cfg['title']}\n"
                      f"(n={n_events} events, {n_animals} animals, "
                      f"y=0 = individual breakpoint)")
        fig.tight_layout()
        save_figure(fig, cfg["fname"], out_dir)


# ─────────────────────────────────────────────────────────────
#  « thermoneutral fraction vs milk yield »
# ─────────────────────────────────────────────────────────────


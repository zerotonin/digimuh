#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — viz_production                                       ║
# ║  « thermoneutral fraction vs milk yield figures »               ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Pooled and per-cow Spearman correlation between daily TNF      ║
# ║  and P95-normalised milk yield.                                 ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Production impact figures for the Frontiers manuscript."""

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
#  « thermoneutral fraction vs milk yield »
# ─────────────────────────────────────────────────────────────

def plot_tnf_yield(out_dir: Path) -> None:
    """Scatter plots: daily thermoneutral fraction vs daily milk yield.

    Each dot is one cow-day.  Panel A shows absolute yield, Panel B
    shows relative yield (daily yield / cow-specific P95).
    """
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr
    setup_figure()

    tnf_path = resolve_input(out_dir, "tnf_yield.csv")
    if not tnf_path.exists():
        log.info("  tnf_yield.csv not found, skipping TNF plots")
        return

    df = pd.read_csv(tnf_path)
    if df.empty:
        return

    valid = df.dropna(subset=["thi_tnf", "daily_yield_kg", "relative_yield"])
    if len(valid) < 20:
        log.info("  Too few cow-days for TNF vs yield plot (%d)", len(valid))
        return

    n_cows = valid["animal_id"].nunique()
    log.info("  Plotting TNF vs yield (%d cow-days, %d animals) …",
             len(valid), n_cows)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # ── Panel A: TNF vs absolute daily yield ──────────────────
    ax = axes[0]
    for year in sorted(valid["year"].unique().astype(int)):
        sub = valid[valid["year"] == year]
        ax.scatter(sub["thi_tnf"], sub["daily_yield_kg"],
                   s=4, alpha=0.15, label=str(year),
                   color=COLOURS["year"].get(year, "#888"))

    rs, p = spearmanr(valid["thi_tnf"], valid["daily_yield_kg"])
    z = np.polyfit(valid["thi_tnf"], valid["daily_yield_kg"], 1)
    xl = np.linspace(0, 1, 50)
    ax.plot(xl, np.polyval(z, xl), "--", color=COLOURS["fit_line"], linewidth=2)

    ax.set_xlabel("Daily thermoneutral fraction (TNF)")
    ax.set_ylabel("Daily milk yield (kg)")
    ax.set_title(f"TNF vs daily yield\n"
                 f"rs = {rs:.3f}, p = {p:.2e}, n = {len(valid):,} cow-days")
    ax.legend(fontsize=8, title="Year", markerscale=3)

    # ── Panel B: TNF vs relative yield (cow-specific P95) ─────
    ax = axes[1]
    for year in sorted(valid["year"].unique().astype(int)):
        sub = valid[valid["year"] == year]
        ax.scatter(sub["thi_tnf"], sub["relative_yield"],
                   s=4, alpha=0.15, label=str(year),
                   color=COLOURS["year"].get(year, "#888"))

    rs_rel, p_rel = spearmanr(valid["thi_tnf"], valid["relative_yield"])
    z = np.polyfit(valid["thi_tnf"], valid["relative_yield"], 1)
    ax.plot(xl, np.polyval(z, xl), "--", color=COLOURS["fit_line"], linewidth=2)
    ax.axhline(1.0, color="#999", linestyle=":", linewidth=0.8, label="P95 reference")

    ax.set_xlabel("Daily thermoneutral fraction (TNF)")
    ax.set_ylabel("Relative daily yield (yield / cow P95)")
    ax.set_title(f"TNF vs relative yield\n"
                 f"rs = {rs_rel:.3f}, p = {p_rel:.2e}, n = {len(valid):,} cow-days")
    ax.legend(fontsize=8, title="Year", markerscale=3)

    fig.suptitle(f"Thermoneutral fraction vs milk yield "
                 f"({n_cows} animals)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "tnf_yield", out_dir)


# ─────────────────────────────────────────────────────────────
#  « TNF vs yield, stratified by Wood-residual yield class »
# ─────────────────────────────────────────────────────────────

_CLASS_ORDER = ("low", "middle", "high")
_CLASS_COLOURS = {
    "low":    "#0072B2",   # Wong blue
    "middle": "#999999",   # Wong grey
    "high":   "#D55E00",   # Wong vermillion
}
_PREDICTOR_LABELS = {
    "thi_tnf":  "Daily time below THI breakpoint",
    "temp_tnf": "Daily time below barn-temp breakpoint",
}


def plot_tnf_yield_by_class(
    tnf_by_class: pd.DataFrame,
    correlations: pd.DataFrame,
    out_dir: Path,
    response: str = "daily_yield_kg",
) -> None:
    """Scatter grid — 3 yield classes × 2 TNF predictors.

    Each panel shows cow-day points coloured by class, an OLS
    reference line, and a Spearman annotation taken directly from
    the ``correlations`` table so the figure and the console
    summary never disagree.

    Args:
        tnf_by_class: Output of
            :func:`stats_production.compute_tnf_yield_by_class`.
        correlations: Output of
            :func:`stats_production.tnf_yield_correlations_by_class`.
        out_dir: Figure output directory.
        response: Column name for the y-axis — usually
            ``"daily_yield_kg"`` (default) or ``"relative_yield"``.
    """
    import matplotlib.pyplot as plt
    setup_figure()

    if tnf_by_class.empty:
        log.info("  No stratified TNF × yield data, skipping")
        return

    predictors = [p for p in ("thi_tnf", "temp_tnf")
                  if p in tnf_by_class.columns]
    if not predictors:
        return
    if response not in tnf_by_class.columns:
        log.info("  Response '%s' missing from tnf_by_class, skipping",
                 response)
        return

    y_label = {
        "daily_yield_kg":  "Daily milk yield (kg/d)",
        "relative_yield":  "Relative daily yield (yield / cow P95)",
        "yield_residual":  "DIM-adjusted yield residual (kg/d)",
    }.get(response, response)

    fig, axes = plt.subplots(
        len(_CLASS_ORDER), len(predictors),
        figsize=(6.0 * len(predictors), 3.8 * len(_CLASS_ORDER)),
        sharex="col", sharey="row",
        squeeze=False,
    )

    for row, cls in enumerate(_CLASS_ORDER):
        cls_df = tnf_by_class[tnf_by_class["yield_class"] == cls]
        colour = _CLASS_COLOURS[cls]
        n_animals = cls_df["animal_id"].nunique()

        for col, pred in enumerate(predictors):
            ax = axes[row, col]
            valid = cls_df.dropna(subset=[pred, response])

            if len(valid) < 5:
                ax.text(0.5, 0.5,
                        f"no data ({cls}, {pred})",
                        transform=ax.transAxes, ha="center", va="center",
                        color="#666")
                ax.set_axis_off()
                continue

            ax.scatter(valid[pred], valid[response],
                       s=4, alpha=0.22, color=colour, edgecolors="none")

            # OLS reference line (consistent with correlations table)
            stat = correlations[
                (correlations["yield_class"] == cls)
                & (correlations["predictor"] == pred)
                & (correlations["response"] == response)
            ]
            if not stat.empty and np.isfinite(stat["slope"].iloc[0]):
                slope = float(stat["slope"].iloc[0])
                intercept = float(stat["intercept"].iloc[0])
                xl = np.linspace(0, 1, 50)
                ax.plot(xl, intercept + slope * xl,
                        "--", color="#333", linewidth=1.6)

                rs = float(stat["rs"].iloc[0])
                p  = float(stat["p"].iloc[0])
                anno = (f"rs = {rs:+.3f}\n"
                        f"p = {p:.2e}\n"
                        f"n = {int(stat['n'].iloc[0]):,} cow-days\n"
                        f"animals = {n_animals}")
                ax.text(0.02, 0.97, anno, transform=ax.transAxes,
                        va="top", ha="left", fontsize=9, color="#222",
                        bbox=dict(boxstyle="round,pad=0.35",
                                  facecolor="white", alpha=0.85,
                                  edgecolor=_CLASS_COLOURS[cls]))

            if col == 0:
                ax.set_ylabel(f"{cls.capitalize()} yield class\n{y_label}",
                              fontsize=10)
            if row == len(_CLASS_ORDER) - 1:
                ax.set_xlabel(_PREDICTOR_LABELS.get(pred, pred))
            if row == 0:
                ax.set_title(_PREDICTOR_LABELS.get(pred, pred),
                             fontsize=11)

            ax.set_xlim(-0.02, 1.02)
            ax.grid(False)
            if response == "yield_residual":
                ax.axhline(0, color="#999", linewidth=0.8,
                           linestyle=":", zorder=0)

    suffix_map = {
        "daily_yield_kg":  "",
        "relative_yield":  "_relative",
        "yield_residual":  "_residual",
    }
    suffix = suffix_map.get(response, f"_{response}")
    title_map = {
        "daily_yield_kg":  "Thermoneutral fraction vs daily milk yield",
        "relative_yield":  "Thermoneutral fraction vs relative yield "
                           "(yield / cow P95)",
        "yield_residual":  "Thermoneutral fraction vs DIM-adjusted yield "
                           "residual (Wood 1967)",
    }
    fig.suptitle(
        f"{title_map.get(response, response)} — by Wood-residual "
        "yield class",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout()
    save_figure(fig, f"tnf_yield_by_class{suffix}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « crossing-day raincloud (crossed vs not-crossed) »
# ─────────────────────────────────────────────────────────────

def _raincloud_one(
    ax, y_yes, y_no, colour_yes, colour_no, labels=("Crossed", "Not crossed"),
):
    """Draw a two-group horizontal raincloud at positions (1, 2)."""
    from scipy.stats import gaussian_kde

    data = {labels[0]: np.asarray(y_yes, dtype=float),
            labels[1]: np.asarray(y_no, dtype=float)}
    colours = {labels[0]: colour_yes, labels[1]: colour_no}

    for i, (name, vals) in enumerate(data.items(), start=1):
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            continue
        colour = colours[name]
        # Half-violin to the LEFT of the centre line.
        if len(vals) > 10:
            try:
                kde = gaussian_kde(vals, bw_method=0.3)
                y_grid = np.linspace(vals.min(), vals.max(), 200)
                dens = kde(y_grid)
                dens = dens / dens.max() * 0.35
                ax.fill_betweenx(y_grid, i - dens, i,
                                 alpha=0.35, color=colour)
                ax.plot(i - dens, y_grid, color=colour, linewidth=1.2)
            except (np.linalg.LinAlgError, ValueError):
                pass
        # Jittered scatter to the RIGHT.
        rng = np.random.default_rng(i)
        jitter = rng.uniform(0.08, 0.32, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals,
                   s=2.5, alpha=0.15, color=colour, edgecolors="none")
        # Boxplot overlay on the centre.
        ax.boxplot(
            vals, positions=[i + 0.45], widths=0.12, vert=True,
            patch_artist=True,
            boxprops=dict(facecolor=colour, alpha=0.5, edgecolor="#333"),
            medianprops=dict(color="#333", linewidth=1.8),
            whiskerprops=dict(color="#333"),
            capprops=dict(color="#333"),
            flierprops=dict(marker="o", markersize=1.5, alpha=0.25),
        )


def plot_crossing_day_raincloud(
    df: pd.DataFrame,
    comparison: pd.DataFrame,
    out_dir: Path,
    response: str = "yield_residual",
    group: str = "pooled",
) -> None:
    """Raincloud pair per climate predictor — crossed-day vs not.

    One figure, two panels: THI crossings (left), barn-temp
    crossings (right).  For each, two rainclouds compare the
    response on cow-days with ≥1 breakpoint crossing of that
    predictor vs cow-days with none.

    Args:
        df: Cow-day DataFrame with ``thi_crossed``, ``temp_crossed``,
            the response column, and (if filtering) ``yield_class``.
        comparison: Output of
            :func:`stats_production.crossing_day_comparison` with
            the same ``response``.  Used to pull the Mann-Whitney
            annotation.
        out_dir: Figure output directory.
        response: Response column, default ``yield_residual``.
        group: Which class row to plot; ``"pooled"`` or one of
            ``"low"``, ``"middle"``, ``"high"``.
    """
    import matplotlib.pyplot as plt
    setup_figure()

    if df.empty or response not in df.columns:
        return
    sub = df if group == "pooled" else df[df["yield_class"] == group]
    if sub.empty:
        return

    y_label = {
        "daily_yield_kg":  "Daily milk yield (kg/d)",
        "relative_yield":  "Relative daily yield (yield / cow P95)",
        "yield_residual":  "DIM-adjusted yield residual (kg/d)",
    }.get(response, response)

    from digimuh.constants import WONG_VERMILLION, WONG_BLUE, WONG_GREY

    fig, (ax_thi, ax_temp) = plt.subplots(1, 2, figsize=(12, 6),
                                           sharey=True)

    for ax, flag, title_pred, pred_short in [
        (ax_thi,  "thi_crossed",  "THI breakpoint",          "thi"),
        (ax_temp, "temp_crossed", "Barn-temp breakpoint",    "temp"),
    ]:
        if flag not in sub.columns:
            ax.text(0.5, 0.5, f"{flag} missing", transform=ax.transAxes,
                    ha="center", va="center")
            ax.set_axis_off()
            continue
        valid = sub.dropna(subset=[flag, response])
        y_yes = valid.loc[valid[flag] == True, response].values
        y_no  = valid.loc[valid[flag] == False, response].values

        _raincloud_one(ax, y_yes, y_no,
                       colour_yes=WONG_VERMILLION, colour_no=WONG_BLUE,
                       labels=("Crossed", "Not crossed"))

        if response == "yield_residual":
            ax.axhline(0, color=WONG_GREY, linestyle=":", linewidth=0.8,
                       zorder=0)

        # Annotation from the comparison table
        stat = comparison[(comparison["group"] == group)
                          & (comparison["predictor"] == pred_short)]
        if not stat.empty and np.isfinite(stat["p"].iloc[0]):
            s = stat.iloc[0]
            from digimuh.stats_core import p_to_stars
            anno = (f"Crossed  n = {int(s['n_yes']):,}  "
                    f"median = {s['median_yes']:+.2f}\n"
                    f"Not crossed  n = {int(s['n_no']):,}  "
                    f"median = {s['median_no']:+.2f}\n"
                    f"Δ median = {s['median_diff']:+.2f} kg/d  "
                    f"(Mann-Whitney p = {s['p']:.2e} "
                    f"{p_to_stars(s['p'])})")
            ax.text(0.02, 0.97, anno, transform=ax.transAxes,
                    va="top", ha="left", fontsize=9, color="#222",
                    bbox=dict(boxstyle="round,pad=0.35",
                              facecolor="white", alpha=0.85,
                              edgecolor=WONG_GREY))

        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Crossed", "Not crossed"])
        ax.set_xlim(0.3, 2.8)
        ax.set_xlabel(f"{title_pred} on that cow-day")
        ax.set_title(f"{title_pred} crossings")
        ax.grid(False)

    ax_thi.set_ylabel(y_label)

    title_suffix = "all cow-days" if group == "pooled" else f"{group} yield class"
    fig.suptitle(
        f"Crossed-day vs non-crossed-day yield — {title_suffix}",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout()
    fname_suffix = "" if group == "pooled" else f"_{group}"
    save_figure(fig, f"crossing_day_raincloud{fname_suffix}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « daily-mean climate vs yield scatter fit »
# ─────────────────────────────────────────────────────────────

def plot_daily_climate_vs_yield(
    df: pd.DataFrame,
    correlations: pd.DataFrame,
    out_dir: Path,
    response: str = "yield_residual",
) -> None:
    """Scatter of daily-mean THI / barn-temp vs a yield response.

    One figure, two panels — mean THI left, mean barn-temp right.
    A single OLS fit across ALL cow-days per panel; the Spearman
    rs, p and slope come from ``correlations`` so figure and
    console agree.
    """
    import matplotlib.pyplot as plt
    from digimuh.constants import WONG_SKY, WONG_VERMILLION, WONG_GREY
    setup_figure()

    if df.empty or response not in df.columns:
        return
    if correlations.empty:
        return

    y_label = {
        "daily_yield_kg":  "Daily milk yield (kg/d)",
        "relative_yield":  "Relative daily yield (yield / cow P95)",
        "yield_residual":  "DIM-adjusted yield residual (kg/d)",
    }.get(response, response)

    fig, (ax_thi, ax_temp) = plt.subplots(1, 2, figsize=(13, 5.5))

    for ax, pred, x_label in [
        (ax_thi,  "mean_thi",       "Daily mean barn THI"),
        (ax_temp, "mean_barn_temp", "Daily mean barn temperature (°C)"),
    ]:
        if pred not in df.columns:
            ax.text(0.5, 0.5, f"{pred} missing",
                    transform=ax.transAxes, ha="center", va="center")
            ax.set_axis_off()
            continue
        valid = df.dropna(subset=[pred, response])
        if valid.empty:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                    ha="center", va="center")
            ax.set_axis_off()
            continue

        ax.scatter(valid[pred], valid[response],
                   s=3, alpha=0.15, color=WONG_SKY, edgecolors="none")

        stat = correlations[correlations["predictor"] == pred]
        if not stat.empty and np.isfinite(stat["slope"].iloc[0]):
            s = stat.iloc[0]
            xl = np.linspace(valid[pred].min(), valid[pred].max(), 100)
            ax.plot(xl, s["intercept"] + s["slope"] * xl,
                    "--", color=WONG_VERMILLION, linewidth=2,
                    label=f"OLS slope = {s['slope']:+.3f}")
            from digimuh.stats_core import p_to_stars
            anno = (f"rs = {s['rs']:+.3f}\n"
                    f"p = {s['p']:.2e} {p_to_stars(s['p'])}\n"
                    f"n = {int(s['n']):,} cow-days\n"
                    f"animals = {int(s['n_animals'])}")
            ax.text(0.02, 0.97, anno, transform=ax.transAxes,
                    va="top", ha="left", fontsize=9, color="#222",
                    bbox=dict(boxstyle="round,pad=0.35",
                              facecolor="white", alpha=0.85,
                              edgecolor=WONG_GREY))

        if response == "yield_residual":
            ax.axhline(0, color=WONG_GREY, linestyle=":",
                       linewidth=0.8, zorder=0)

        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.grid(False)
        ax.legend(fontsize=8, loc="lower left")

    fig.suptitle("Daily-mean climate vs DIM-adjusted yield"
                 if response == "yield_residual"
                 else f"Daily-mean climate vs {y_label}",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    suffix_map = {
        "daily_yield_kg":  "",
        "relative_yield":  "_relative",
        "yield_residual":  "_residual",
    }
    suffix = suffix_map.get(response, f"_{response}")
    save_figure(fig, f"daily_climate_vs_yield{suffix}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « longitudinal breakpoint tracking (repeat animals) »
# ─────────────────────────────────────────────────────────────


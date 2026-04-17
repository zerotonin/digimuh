#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — viz_milk_composition                                 ║
# ║  « MLP × climate figures (thin-milk hypothesis) »              ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Figures for the MLP composition × climate analysis.

Two figures:

* ``mlp_thin_milk_hypothesis.{svg,png}`` — a compact 2×2 scatter
  directly testing the hypothesis: milk volume and ECM on the top
  row (expected to rise with climate if dilution is real), fat %
  and protein % on the bottom row (expected to fall).  Each panel
  has an OLS reference line and a Spearman annotation.

* ``mlp_composition_heatmap.{svg,png}`` — a signed-rs heatmap
  across every MLP response × every climate predictor × each
  Wood-residual yield class, so patterns in SCC / urea / lactose
  and class-specific differences come out at a glance.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.constants import (
    WONG_BLUE, WONG_VERMILLION, WONG_SKY, WONG_ORANGE, WONG_GREY,
)
from digimuh.viz_base import setup_figure, save_figure

log = logging.getLogger("digimuh.viz")


_HYPO_PANELS: tuple[tuple[str, str, str, str], ...] = (
    # (response_col, response_label, expected_sign_hint, colour)
    ("herdeplus_mlp_mkg",             "Test-day milk (kg/d)",
     "↑ with heat if dilution",       WONG_BLUE),
    ("herdeplus_mlp_ecm",             "ECM (kg/d)",
     "↓ with heat if total solids drop", WONG_SKY),
    ("herdeplus_mlp_fat_percent",     "Fat (%)",
     "↓ with heat if dilution",       WONG_VERMILLION),
    ("herdeplus_mlp_protein_percent", "Protein (%)",
     "↓ with heat if dilution",       WONG_ORANGE),
)


def _spearman_annotation(ax, corr_row, xl, fit_colour):
    """Draw an OLS reference line + a rs/p/n annotation box."""
    if corr_row.empty:
        return
    row = corr_row.iloc[0]
    if not np.isfinite(row.get("slope", np.nan)):
        return
    slope = float(row["slope"])
    intercept = float(row["intercept"])
    ax.plot(xl, intercept + slope * xl, "--", color=fit_colour, linewidth=1.8,
            label=f"OLS slope = {slope:+.3f}")
    from digimuh.stats_core import p_to_stars
    rs = float(row["rs"])
    p  = float(row["p"])
    anno = (f"rs = {rs:+.3f}\n"
            f"p = {p:.2e} {p_to_stars(p)}\n"
            f"n = {int(row['n']):,} test-days\n"
            f"animals = {int(row['n_animals'])}")
    ax.text(0.02, 0.97, anno, transform=ax.transAxes,
            va="top", ha="left", fontsize=9, color="#222",
            bbox=dict(boxstyle="round,pad=0.35",
                      facecolor="white", alpha=0.85,
                      edgecolor=WONG_GREY))


def plot_thin_milk_hypothesis(
    merged: pd.DataFrame,
    correlations: pd.DataFrame,
    out_dir: Path,
    predictor: str = "mean_thi",
    predictor_label: str = "Daily mean barn THI",
) -> None:
    """2×2 scatter panel directly testing the dilution / thin-milk story."""
    import matplotlib.pyplot as plt
    setup_figure()

    if merged.empty or predictor not in merged.columns:
        log.info("  Thin-milk plot skipped: no MLP × climate merged data")
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    axes_flat = axes.flatten()

    xl = np.linspace(
        np.nanpercentile(merged[predictor].dropna(), 1),
        np.nanpercentile(merged[predictor].dropna(), 99),
        80,
    )

    for ax, (resp, label, hint, colour) in zip(axes_flat, _HYPO_PANELS):
        if resp not in merged.columns:
            ax.set_axis_off()
            continue
        valid = merged.dropna(subset=[predictor, resp])
        ax.scatter(valid[predictor], valid[resp],
                   s=10, alpha=0.35, color=colour, edgecolors="none")

        crow = correlations[(correlations["group"] == "pooled")
                            & (correlations["predictor"] == predictor)
                            & (correlations["response"] == resp)]
        _spearman_annotation(ax, crow, xl, fit_colour=WONG_VERMILLION)

        ax.set_ylabel(label)
        ax.set_title(f"{label}\n{hint}", fontsize=10)
        ax.grid(False)

    for ax in axes[-1]:
        ax.set_xlabel(predictor_label)

    fig.suptitle("Thin-milk hypothesis: does composition dilute as climate warms?",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "mlp_thin_milk_hypothesis", out_dir)


# ─────────────────────────────────────────────────────────────
#  « dilution partition: observed vs pure-water prediction »
# ─────────────────────────────────────────────────────────────

def plot_dilution_partition(
    df: pd.DataFrame,
    summary: pd.DataFrame,
    out_dir: Path,
    predictor: str = "mean_thi",
    predictor_label: str = "Daily mean barn THI",
) -> None:
    """Overlay observed vs dilution-predicted fat % / protein %.

    Two panels (fat % top, protein % bottom) showing three OLS
    reference lines per panel:
        * **observed**            (vermillion)  — the actual drop
        * **dilution predicted**  (blue, dashed) — the drop the
                                  volume increase alone can explain
                                  if the cow's reference fat/protein
                                  kg were held constant
        * **rumen residual**      (grey, dotted) — observed − dilution;
                                  how much extra drop comes from
                                  reduced absolute fat/protein output

    A zero-rs rumen residual means pure dilution; a negative
    residual means rumen suppression on top of dilution.
    """
    import matplotlib.pyplot as plt
    setup_figure()

    if df.empty or summary.empty or predictor not in df.columns:
        log.info("  Dilution partition plot skipped: no data")
        return

    if "fat_percent_diluted" not in df.columns:
        log.info("  Dilution partition plot skipped: compute_dilution_partition "
                 "was not run on this frame")
        return

    specs = [
        ("Fat %", "herdeplus_mlp_fat_percent",
         "fat_percent_diluted", "fat_percent_rumen", "fat"),
        ("Protein %", "herdeplus_mlp_protein_percent",
         "protein_percent_diluted", "protein_percent_rumen", "protein"),
        ("Lactose %", "herdeplus_mlp_lactose",
         "lactose_percent_diluted", "lactose_percent_rumen", "lactose"),
    ]

    fig, axes = plt.subplots(len(specs), 1, figsize=(11, 11), sharex=True)
    if len(specs) == 1:
        axes = [axes]

    x_valid = df[predictor].dropna()
    xl = np.linspace(
        float(np.percentile(x_valid, 1)),
        float(np.percentile(x_valid, 99)),
        80,
    )

    for ax, (label, obs_col, dil_col, rum_col, nutrient) in zip(axes, specs):
        valid = df.dropna(subset=[predictor, obs_col, dil_col, rum_col])
        ax.scatter(valid[predictor], valid[obs_col],
                   s=8, alpha=0.22, color=WONG_VERMILLION,
                   edgecolors="none", label="observed")

        # OLS lines from the summary table (consistent with the numbers
        # reported in the console table).
        from digimuh.stats_core import p_to_stars
        for comp_label, colour, linestyle in [
            ("observed",           WONG_VERMILLION, "-"),
            ("dilution_predicted", WONG_BLUE,       "--"),
            ("rumen_residual",     WONG_GREY,       ":"),
        ]:
            row = summary[(summary["nutrient"] == nutrient)
                          & (summary["component"] == comp_label)]
            if row.empty or not np.isfinite(row["slope"].iloc[0]):
                continue
            s = row.iloc[0]
            ax.plot(xl, s["intercept"] + s["slope"] * xl,
                    linestyle=linestyle, linewidth=2, color=colour,
                    label=(f"{comp_label} — "
                           f"rs={s['rs']:+.3f}{p_to_stars(s['p'])}, "
                           f"slope={s['slope']:+.4f}"))

        # Rumen residual curve as a lighter scatter so the reader can see
        # where the residual (observed − predicted) sits per point.
        ax.axhline(0, color=WONG_GREY, linewidth=0.6, linestyle=":",
                   alpha=0.4, zorder=0)

        ax.set_ylabel(label)
        ax.legend(fontsize=9, loc="upper right")
        ax.grid(False)

    axes[-1].set_xlabel(predictor_label)
    fig.suptitle("Dilution partition — is the composition drop pure added water?",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "mlp_dilution_partition", out_dir)


# ─────────────────────────────────────────────────────────────
#  « full heatmap of composition × climate »
# ─────────────────────────────────────────────────────────────

def plot_composition_heatmap(
    correlations: pd.DataFrame,
    out_dir: Path,
    predictor: str = "mean_thi",
) -> None:
    """Signed-rs heatmap: rows = MLP metrics, columns = yield class."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    setup_figure()

    if correlations.empty:
        return
    sub = correlations[correlations["predictor"] == predictor].copy()
    if sub.empty:
        return

    class_order = ["low", "middle", "high", "pooled"]
    present_classes = [c for c in class_order if c in sub["group"].unique()]
    if not present_classes:
        return

    pivot_rs = sub.pivot_table(
        index="response", columns="group", values="rs",
    ).reindex(columns=present_classes)
    pivot_p = sub.pivot_table(
        index="response", columns="group", values="p",
    ).reindex(columns=present_classes)
    pivot_label = sub.drop_duplicates(
        subset="response").set_index("response")["response_label"]

    row_order = [c for c, _, _ in
                 __import__("digimuh.stats_milk_composition",
                            fromlist=["MLP_RESPONSES"]).MLP_RESPONSES
                 if c in pivot_rs.index]
    pivot_rs = pivot_rs.loc[row_order]
    pivot_p  = pivot_p.loc[row_order]

    fig_h = max(4.0, 0.55 * len(row_order) + 1.2)
    fig, ax = plt.subplots(figsize=(1.2 + 1.4 * len(present_classes), fig_h))

    # Symmetric diverging norm, capped at ±0.4 so tiny rs don't wash out.
    max_abs = float(np.nanmax(np.abs(pivot_rs.values))) if pivot_rs.size else 0.1
    vmax = max(0.1, min(0.4, max_abs))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    im = ax.imshow(pivot_rs.values, cmap="RdBu_r", norm=norm, aspect="auto")

    # Annotate each cell with rs (stars inline).
    from digimuh.stats_core import p_to_stars
    for i, row_col in enumerate(row_order):
        for j, grp in enumerate(present_classes):
            rs = pivot_rs.iloc[i, j]
            p = pivot_p.iloc[i, j]
            if not np.isfinite(rs):
                continue
            stars = p_to_stars(p) if np.isfinite(p) else ""
            ax.text(j, i, f"{rs:+.2f}{stars}",
                    ha="center", va="center",
                    color="#111" if abs(rs) < 0.25 else "white",
                    fontsize=9)

    ax.set_xticks(range(len(present_classes)))
    ax.set_xticklabels([c.capitalize() for c in present_classes])
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels([pivot_label.get(c, c) for c in row_order])
    ax.set_xlabel("Yield class")
    ax.set_title(f"MLP composition × {predictor} — Spearman rs")
    ax.grid(False)

    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.04)
    cbar.set_label("Spearman rs")

    fig.tight_layout()
    save_figure(fig, f"mlp_composition_heatmap_{predictor}", out_dir)

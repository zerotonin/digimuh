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

    tnf_path = out_dir / "tnf_yield.csv"
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
#  « longitudinal breakpoint tracking (repeat animals) »
# ─────────────────────────────────────────────────────────────


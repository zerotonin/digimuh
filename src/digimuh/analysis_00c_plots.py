#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  Analysis 0c — Plotting ¹                                    ║
# ║  « publication figures with significance brackets »          ║
# ╚══════════════════════════════════════════════════════════════╝
"""Reads CSVs from ``analysis_00a_extract`` and ``analysis_00b_stats``,
generates all figures for the Frontiers manuscript.

Usage::

    digimuh-plots --data results/broken_stick

¹ Analysis led by Dr. med. vet. Gundula Hoffmann, ATB Potsdam.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger("digimuh.plots")


# ─────────────────────────────────────────────────────────────
#  « colour scheme — Wong (2011), colourblind-safe »
# ─────────────────────────────────────────────────────────────

COLOURS = {
    "year": {
        2021: "#0072B2",
        2022: "#E69F00",
        2023: "#009E73",
        2024: "#CC79A7",
    },
    "below_bp": "#0072B2",
    "above_bp": "#D55E00",
    "fit_line": "#D55E00",
    "reference": "#E69F00",
    "identity": "#999999",
    "scatter": "#56B4E9",
    "scatter_alt": "#009E73",
    "hist_thi": "#0072B2",
    "hist_temp": "#009E73",
    "box_thi": "#56B4E9",
    "box_temp": "#009E73",
    "median": "#D55E00",
    "paired_line": "#999999",
}


# ─────────────────────────────────────────────────────────────
#  « plotting setup and helpers »
# ─────────────────────────────────────────────────────────────

def _setup():
    """Configure matplotlib for publication-quality figures."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")
    mpl.rcParams.update({
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })


def _save(fig, name: str, out_dir: Path):
    """Save figure as SVG + PNG and close."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(out_dir / f"{name}.{ext}")
    import matplotlib.pyplot as plt
    plt.close(fig)


def add_significance_bracket(
    ax, x1: float, x2: float, y: float, stars: str,
    h: float = 0.02, lw: float = 1.2,
) -> None:
    """Draw a significance bracket with stars between two x positions.

    Args:
        ax: Matplotlib axes.
        x1, x2: Left and right x positions.
        y: Y position of the bracket bottom (data coords).
        stars: Text to display (e.g. '***', 'n.s.').
        h: Bracket height as fraction of y-range.
        lw: Line width.
    """
    if not stars or stars == "":
        return
    yrange = ax.get_ylim()[1] - ax.get_ylim()[0]
    dh = h * yrange
    ax.plot([x1, x1, x2, x2], [y, y + dh, y + dh, y], color="#333333",
            linewidth=lw, clip_on=False)
    ax.text((x1 + x2) / 2, y + dh, stars, ha="center", va="bottom",
            fontsize=11, fontweight="bold", color="#333333")


# ─────────────────────────────────────────────────────────────
#  « breakpoint boxplots (grouped: rumen vs respiration) »
# ─────────────────────────────────────────────────────────────

def plot_grouped_boxplots(bs: pd.DataFrame, out_dir: Path) -> None:
    """Side-by-side boxplots: rumen vs respiration breakpoints per year."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    _setup()

    years = sorted(bs["year"].dropna().unique().astype(int))

    for env_col, env_label, fname in [
        ("thi", "THI", "grouped_thi"),
        ("temp", "Barn temperature (°C)", "grouped_temp"),
    ]:
        rumen_col = f"{env_col}_breakpoint"
        resp_col = f"resp_{env_col}_breakpoint"
        rumen_conv = bs[bs[f"{env_col}_converged"] == True]
        resp_conv = bs[bs[f"resp_{env_col}_converged"] == True]

        fig, ax = plt.subplots(figsize=(10, 6))
        width = 0.35
        x_pos = np.arange(len(years))

        data_r = [rumen_conv[rumen_conv["year"] == y][rumen_col].dropna().values for y in years]
        data_resp = [resp_conv[resp_conv["year"] == y][resp_col].dropna().values for y in years]

        bp1 = ax.boxplot(
            data_r, positions=x_pos - width / 2, widths=width,
            patch_artist=True,
            boxprops=dict(facecolor=COLOURS["below_bp"], alpha=0.6, edgecolor="#333"),
            medianprops=dict(color="#333", linewidth=2),
            manage_ticks=False,
        )
        has_resp = any(len(d) > 0 for d in data_resp)
        if has_resp:
            bp2 = ax.boxplot(
                data_resp, positions=x_pos + width / 2, widths=width,
                patch_artist=True,
                boxprops=dict(facecolor=COLOURS["above_bp"], alpha=0.6, edgecolor="#333"),
                medianprops=dict(color="#333", linewidth=2),
                manage_ticks=False,
            )

        if "thi" in env_col:
            ax.axhline(68.8, color=COLOURS["reference"], linestyle="--",
                       linewidth=1, label="THI 68.8 (mild stress)")

        ax.set_xticks(x_pos)
        ax.set_xticklabels([str(y) for y in years])
        ax.set_xlabel("Year")
        ax.set_ylabel(f"Individual {env_label} breakpoint")
        ax.set_title(f"Per-animal {env_label} breakpoints: rumen temperature vs. respiration rate")
        ax.legend(handles=[
            Patch(facecolor=COLOURS["below_bp"], alpha=0.6, label="Rumen temperature"),
            Patch(facecolor=COLOURS["above_bp"], alpha=0.6, label="Respiration rate"),
        ], fontsize=10)

        for i, y in enumerate(years):
            ymin = ax.get_ylim()[0]
            ax.text(i - width / 2, ymin + 0.5, f"n={len(data_r[i])}",
                    ha="center", fontsize=8, color="#555")
            if has_resp and len(data_resp[i]) > 0:
                ax.text(i + width / 2, ymin + 0.5, f"n={len(data_resp[i])}",
                        ha="center", fontsize=8, color="#555")

        _save(fig, f"boxplot_{fname}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « paired rumen vs resp with significance brackets »
# ─────────────────────────────────────────────────────────────

def plot_paired_below_above(
    beh: pd.DataFrame, tests: pd.DataFrame, out_dir: Path,
) -> None:
    """Paired boxplots: below vs above breakpoint with significance brackets."""
    import matplotlib.pyplot as plt
    from scipy.stats import wilcoxon
    _setup()

    years = sorted(beh["year"].dropna().unique().astype(int))

    for year in years:
        yr = beh[beh["year"] == year]

        panels = [
            ("body_temp_below", "body_temp_above", "Rumen temperature (°C)",
             "body_temp below vs above"),
        ]
        # Add respiration if data exists
        if yr["resp_below"].notna().sum() >= 5:
            panels.append(
                ("resp_below", "resp_above", "Respiration rate (bpm)",
                 "respiration below vs above"),
            )

        fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 6))
        if len(panels) == 1:
            axes = [axes]

        for ax, (below_col, above_col, ylabel, test_name) in zip(axes, panels):
            paired = yr.dropna(subset=[below_col, above_col])
            if len(paired) < 5:
                ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes,
                        ha="center", color=COLOURS["identity"])
                continue

            below_vals = paired[below_col].values
            above_vals = paired[above_col].values

            bp_plot = ax.boxplot(
                [below_vals, above_vals],
                positions=[0, 1],
                labels=["Below\nbreakpoint", "Above\nbreakpoint"],
                patch_artist=True, widths=0.5,
            )
            bp_plot["boxes"][0].set_facecolor(COLOURS["below_bp"])
            bp_plot["boxes"][0].set_alpha(0.5)
            bp_plot["boxes"][1].set_facecolor(COLOURS["above_bp"])
            bp_plot["boxes"][1].set_alpha(0.5)
            for med in bp_plot["medians"]:
                med.set_color("#333")
                med.set_linewidth(2)

            # Paired lines
            for b, a in zip(below_vals, above_vals):
                ax.plot([0, 1], [b, a], "-", color=COLOURS["paired_line"],
                        alpha=0.3, linewidth=0.7)

            # Significance bracket from tests CSV
            test_row = tests[
                (tests["year"] == year) & (tests["test"].str.contains(test_name.split()[0]))
            ]
            if not test_row.empty:
                stars = test_row.iloc[0]["stars"]
                p_adj = test_row.iloc[0]["p_adj"]
            else:
                # Compute on the fly if not in tests
                try:
                    _, p_raw = wilcoxon(below_vals, above_vals)
                    stars = "***" if p_raw < 0.001 else "**" if p_raw < 0.01 else "*" if p_raw < 0.05 else "n.s."
                    p_adj = p_raw
                except Exception:
                    stars = ""
                    p_adj = np.nan

            if stars:
                ymax = max(below_vals.max(), above_vals.max())
                add_significance_bracket(ax, 0, 1, ymax * 1.02, stars)

            ax.set_ylabel(ylabel)
            ax.set_title(f"n = {len(paired)}", fontsize=10)

        fig.suptitle(
            f"Below vs. above individual THI breakpoint ({year})\n"
            f"(paired per-animal means, Wilcoxon signed-rank, BH-FDR corrected)",
            fontsize=12,
        )
        fig.tight_layout()
        _save(fig, f"paired_below_above_{year}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « paired rumen vs respiration breakpoints »
# ─────────────────────────────────────────────────────────────

def plot_paired_rumen_vs_resp(bs: pd.DataFrame, out_dir: Path) -> None:
    """Paired boxplot: rumen temp breakpoint vs respiration breakpoint."""
    import matplotlib.pyplot as plt
    from scipy.stats import wilcoxon
    _setup()

    both = bs[bs["thi_converged"] & bs["resp_thi_converged"]].dropna(
        subset=["thi_breakpoint", "resp_thi_breakpoint"])

    if len(both) < 5:
        log.info("Too few paired rumen+resp animals (%d) for paired plot", len(both))
        return

    fig, ax = plt.subplots(figsize=(7, 6))
    rumen_vals = both["thi_breakpoint"].values
    resp_vals = both["resp_thi_breakpoint"].values

    bp_plot = ax.boxplot(
        [rumen_vals, resp_vals], positions=[0, 1],
        labels=["Rumen\ntemperature", "Respiration\nrate"],
        patch_artist=True, widths=0.5,
    )
    bp_plot["boxes"][0].set_facecolor(COLOURS["below_bp"])
    bp_plot["boxes"][0].set_alpha(0.5)
    bp_plot["boxes"][1].set_facecolor(COLOURS["above_bp"])
    bp_plot["boxes"][1].set_alpha(0.5)
    for med in bp_plot["medians"]:
        med.set_color("#333")
        med.set_linewidth(2)

    for r, rp in zip(rumen_vals, resp_vals):
        ax.plot([0, 1], [r, rp], "-", color=COLOURS["paired_line"],
                alpha=0.3, linewidth=0.7)

    try:
        stat, p = wilcoxon(rumen_vals, resp_vals)
        diff = np.median(resp_vals - rumen_vals)
        stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
        ymax = max(rumen_vals.max(), resp_vals.max())
        add_significance_bracket(ax, 0, 1, ymax * 1.02, stars)
        ax.text(0.5, 0.02,
                f"n = {len(both)}, median Δ = {diff:+.1f}",
                transform=ax.transAxes, ha="center", fontsize=9, color="#555")
    except Exception:
        pass

    ax.set_ylabel("Individual THI breakpoint")
    ax.set_title("Paired comparison: rumen temperature vs. respiration rate\nTHI breakpoints")
    _save(fig, "paired_rumen_vs_resp", out_dir)


# ─────────────────────────────────────────────────────────────
#  « Spearman histograms »
# ─────────────────────────────────────────────────────────────

def plot_spearman(spearman: pd.DataFrame, out_dir: Path) -> None:
    """Spearman correlation distributions."""
    import matplotlib.pyplot as plt
    _setup()

    panels = [
        ("thi_rs", "body temp vs. barn THI", COLOURS["hist_thi"]),
        ("temp_rs", "body temp vs. barn temp", COLOURS["hist_temp"]),
    ]
    # Add respiration if available
    if "resp_thi_rs" in spearman.columns and spearman["resp_thi_rs"].notna().sum() > 5:
        panels.append(("resp_thi_rs", "resp rate vs. barn THI", COLOURS["above_bp"]))

    fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 5))
    if len(panels) == 1:
        axes = [axes]

    for ax, (col, label, color) in zip(axes, panels):
        vals = spearman[col].dropna()
        if len(vals) < 5:
            continue
        ax.hist(vals, bins=30, color=color, edgecolor="white", alpha=0.8)
        ax.axvline(vals.median(), color=COLOURS["fit_line"], linestyle="--",
                   label=f"Median = {vals.median():.3f}")
        ax.set_xlabel(f"Spearman rₛ ({label})")
        ax.set_ylabel("Number of animals")
        ax.legend(fontsize=9)

    fig.suptitle("Per-animal Spearman correlations")
    fig.tight_layout()
    _save(fig, "spearman_correlations", out_dir)


# ─────────────────────────────────────────────────────────────
#  « climate time series »
# ─────────────────────────────────────────────────────────────

def plot_climate(climate: pd.DataFrame, out_dir: Path) -> None:
    """Daily barn THI and temperature time series per summer."""
    import matplotlib.pyplot as plt
    _setup()

    years = sorted(climate["year"].unique().astype(int))
    colors = COLOURS["year"]

    for y_col, ylabel, fname, ref_line in [
        ("barn_thi_mean", "Barn THI", "climate_thi", 68.8),
        ("barn_temp_mean", "Barn temperature (°C)", "climate_temp", None),
    ]:
        min_col = y_col.replace("mean", "min")
        max_col = y_col.replace("mean", "max")

        fig, ax = plt.subplots(figsize=(12, 5))
        for year in years:
            sub = climate[climate["year"] == year]
            doy = pd.to_datetime(sub["day"]).dt.dayofyear
            ax.plot(doy, sub[y_col], linewidth=0.8, alpha=0.8,
                    color=colors.get(year, "#888"), label=str(year))
            if min_col in sub.columns:
                ax.fill_between(doy, sub[min_col], sub[max_col],
                                alpha=0.1, color=colors.get(year, "#888"))
        if ref_line is not None:
            ax.axhline(ref_line, color=COLOURS["reference"], linestyle="--",
                       linewidth=1, label=f"THI {ref_line}")
        ax.set_xlabel("Day of year (Jun–Sep)")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=9)
        _save(fig, fname, out_dir)


# ─────────────────────────────────────────────────────────────
#  « predictors: breakpoint vs milk yield / lactation »
# ─────────────────────────────────────────────────────────────

def plot_predictors(bs: pd.DataFrame, out_dir: Path) -> None:
    """Scatter plots: breakpoint vs production parameters."""
    import matplotlib.pyplot as plt
    from scipy.stats import pearsonr
    _setup()

    for bp_col, bp_label, fname in [
        ("thi_breakpoint", "THI breakpoint", "predictors_thi"),
        ("temp_breakpoint", "Barn temp breakpoint (°C)", "predictors_temp"),
    ]:
        conv = bs[bs[f"{bp_col.split('_')[0]}_converged"] == True].dropna(subset=[bp_col])
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Milk yield
        ax = axes[0]
        sub = conv.dropna(subset=["mean_milk_yield_kg"])
        if len(sub) > 10:
            ax.scatter(sub[bp_col], sub["mean_milk_yield_kg"],
                       s=20, alpha=0.5, color=COLOURS["scatter"])
            r, p = pearsonr(sub[bp_col], sub["mean_milk_yield_kg"])
            z = np.polyfit(sub[bp_col], sub["mean_milk_yield_kg"], 1)
            xl = np.linspace(sub[bp_col].min(), sub[bp_col].max(), 50)
            ax.plot(xl, np.polyval(z, xl), "--", color=COLOURS["fit_line"])
            ax.set_title(f"r = {r:.3f}, p = {p:.3f}, n = {len(sub)}", fontsize=10)
        ax.set_xlabel(bp_label)
        ax.set_ylabel("Mean milk yield (kg/milking)")

        # Lactation
        ax = axes[1]
        sub = conv.dropna(subset=["lactation_nr"])
        sub = sub[sub["lactation_nr"] > 0]
        if len(sub) > 10:
            jitter = np.random.uniform(-0.2, 0.2, len(sub))
            ax.scatter(sub[bp_col], sub["lactation_nr"] + jitter,
                       s=20, alpha=0.5, color=COLOURS["scatter_alt"])
            r, p = pearsonr(sub[bp_col], sub["lactation_nr"])
            z = np.polyfit(sub[bp_col], sub["lactation_nr"], 1)
            xl = np.linspace(sub[bp_col].min(), sub[bp_col].max(), 50)
            ax.plot(xl, np.polyval(z, xl), "--", color=COLOURS["fit_line"])
            ax.set_title(f"r = {r:.3f}, p = {p:.3f}, n = {len(sub)}", fontsize=10)
        else:
            ax.text(0.5, 0.5, f"Insufficient data (n={len(sub)})",
                    transform=ax.transAxes, ha="center", color=COLOURS["identity"])
        ax.set_xlabel(bp_label)
        ax.set_ylabel("Lactation number")

        fig.suptitle(f"{bp_label}: association with production parameters")
        fig.tight_layout()
        _save(fig, fname, out_dir)


# ─────────────────────────────────────────────────────────────
#  « breakpoint stability scatter »
# ─────────────────────────────────────────────────────────────

def plot_stability(pairs: pd.DataFrame, icc: float, out_dir: Path) -> None:
    """Year-to-year breakpoint stability scatter."""
    import matplotlib.pyplot as plt
    _setup()

    if pairs.empty:
        return

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(pairs["thi_bp_1"], pairs["thi_bp_2"],
               s=30, alpha=0.6, color=COLOURS["scatter"])
    lims = [
        min(pairs["thi_bp_1"].min(), pairs["thi_bp_2"].min()) - 2,
        max(pairs["thi_bp_1"].max(), pairs["thi_bp_2"].max()) + 2,
    ]
    ax.plot(lims, lims, "--", color=COLOURS["identity"], linewidth=1, label="Identity")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("THI breakpoint (earlier year)")
    ax.set_ylabel("THI breakpoint (later year)")
    ax.set_title(
        f"Within-animal breakpoint stability\n"
        f"ICC = {icc:.3f}, n = {len(pairs)} pairs",
    )
    ax.legend()
    _save(fig, "breakpoint_stability", out_dir)


# ─────────────────────────────────────────────────────────────
#  « example broken-stick fits (rumen + respiration) »
# ─────────────────────────────────────────────────────────────

def plot_examples(
    rumen: pd.DataFrame, resp: pd.DataFrame, bs: pd.DataFrame,
    out_dir: Path,
) -> None:
    """Plot broken-stick fit for the best example animal (rumen + resp)."""
    import matplotlib.pyplot as plt
    # Import the fit function from stats module
    from digimuh.analysis_00b_stats import broken_stick_fit
    _setup()

    for signal, data, response_col, ylabel, prefix, conv_col in [
        ("rumen", rumen, "body_temp", "Rumen temperature (°C)", "thi", "thi_converged"),
        ("rumen", rumen, "body_temp", "Rumen temperature (°C)", "temp", "temp_converged"),
        ("resp", resp, "resp_rate", "Respiration rate (bpm)", "resp_thi", "resp_thi_converged"),
        ("resp", resp, "resp_rate", "Respiration rate (bpm)", "resp_temp", "resp_temp_converged"),
    ]:
        if data.empty:
            continue

        conv = bs[bs[conv_col] == True]
        r2_col = f"{prefix}_r_squared"
        if r2_col not in conv.columns or conv.empty:
            continue

        best = conv.loc[conv[r2_col].idxmax()]
        aid = int(best["animal_id"])
        year = int(best["year"])

        grp = data[(data["animal_id"] == aid) & (data["year"] == year)]
        if len(grp) < 30:
            continue

        env_col = "barn_thi" if "thi" in prefix else "barn_temp"
        env_label = "Barn THI" if "thi" in prefix else "Barn temperature (°C)"
        x_range = (45, 80) if "thi" in prefix else (5, 35)

        fit = broken_stick_fit(grp[env_col].values, grp[response_col].values,
                               x_range=x_range)
        if not fit.get("converged"):
            continue

        bp = fit["breakpoint"]
        fig, ax = plt.subplots(figsize=(9, 6))
        ax.scatter(grp[env_col], grp[response_col], s=2, alpha=0.15,
                   c=COLOURS["identity"])

        xr = np.linspace(grp[env_col].min(), grp[env_col].max(), 200)
        yp = np.where(
            xr <= bp,
            fit["intercept_below"] + fit["slope_below"] * xr,
            fit["intercept_above"] + fit["slope_above"] * xr,
        )
        ax.plot(xr, yp, color=COLOURS["fit_line"], linewidth=2,
                label="Broken-stick fit")
        ax.axvline(bp, color=COLOURS["below_bp"], linestyle="--",
                   linewidth=1.5, label=f"Breakpoint = {bp:.1f}")
        ax.set_xlabel(env_label)
        ax.set_ylabel(ylabel)
        ax.set_title(
            f"Animal {aid} ({year}) — {ylabel.split('(')[0].strip()} "
            f"breakpoint at {env_label} = {bp:.1f}\n"
            f"R² = {fit['r_squared']:.3f}, n = {fit['n']:,}",
        )
        ax.legend()
        _save(fig, f"example_{prefix}_{aid}_{year}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « body temp vs resp breakpoint scatter »
# ─────────────────────────────────────────────────────────────

def plot_bodytemp_vs_resp_scatter(bs: pd.DataFrame, out_dir: Path) -> None:
    """Scatter: rumen temp THI breakpoint vs resp THI breakpoint."""
    import matplotlib.pyplot as plt
    _setup()

    both = bs[bs["thi_converged"] & bs["resp_thi_converged"]].dropna(
        subset=["thi_breakpoint", "resp_thi_breakpoint"])
    if len(both) < 5:
        return

    years = sorted(both["year"].unique().astype(int))
    fig, ax = plt.subplots(figsize=(7, 6))
    for y in years:
        sub = both[both["year"] == y]
        if not sub.empty:
            ax.scatter(sub["thi_breakpoint"], sub["resp_thi_breakpoint"],
                       s=30, alpha=0.7, label=str(y),
                       color=COLOURS["year"].get(y, "#888"))
    lims = [
        min(both["thi_breakpoint"].min(), both["resp_thi_breakpoint"].min()) - 2,
        max(both["thi_breakpoint"].max(), both["resp_thi_breakpoint"].max()) + 2,
    ]
    ax.plot(lims, lims, "--", color=COLOURS["identity"], linewidth=1, label="Identity")
    ax.set_xlabel("Rumen temp THI breakpoint")
    ax.set_ylabel("Respiration rate THI breakpoint")
    ax.set_title("Individual THI breakpoints: rumen temperature vs. respiration")
    ax.legend(fontsize=9)
    _save(fig, "scatter_bodytemp_vs_resp_bp", out_dir)


# ─────────────────────────────────────────────────────────────
#  « THI vs barn temp breakpoint scatter »
# ─────────────────────────────────────────────────────────────

def plot_thi_vs_temp_scatter(bs: pd.DataFrame, out_dir: Path) -> None:
    """Scatter: THI breakpoint vs barn temp breakpoint."""
    import matplotlib.pyplot as plt
    _setup()

    both = bs[bs["thi_converged"] & bs["temp_converged"]].dropna(
        subset=["thi_breakpoint", "temp_breakpoint"])
    if len(both) < 10:
        return

    years = sorted(both["year"].unique().astype(int))
    fig, ax = plt.subplots(figsize=(7, 6))
    for y in years:
        sub = both[both["year"] == y]
        if not sub.empty:
            ax.scatter(sub["temp_breakpoint"], sub["thi_breakpoint"],
                       s=30, alpha=0.7, label=str(y),
                       color=COLOURS["year"].get(y, "#888"))
    ax.set_xlabel("Barn temperature breakpoint (°C)")
    ax.set_ylabel("THI breakpoint")
    ax.set_title("Individual heat stress thresholds: THI vs. barn temperature")
    ax.legend()
    _save(fig, "scatter_thi_vs_temp", out_dir)


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Generate broken-stick analysis figures")
    parser.add_argument("--data", type=Path, required=True,
                        help="Directory with CSVs from extract + stats steps")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(name)-20s │ %(message)s",
        datefmt="%H:%M:%S",
    )
    d = args.data

    log.info("Loading CSVs from %s …", d)
    bs = pd.read_csv(d / "broken_stick_results.csv")
    beh = pd.read_csv(d / "behavioural_response.csv")
    tests = pd.read_csv(d / "statistical_tests.csv")
    spearman = pd.read_csv(d / "spearman_correlations.csv")
    climate = pd.read_csv(d / "climate_daily.csv")
    rumen = pd.read_csv(d / "rumen_barn.csv")

    resp_path = d / "respiration_barn.csv"
    resp = pd.read_csv(resp_path) if resp_path.exists() else pd.DataFrame()

    stability_path = d / "breakpoint_stability.csv"
    pairs = pd.read_csv(stability_path) if stability_path.exists() else pd.DataFrame()

    # Read ICC from stats output
    icc = np.nan
    if not pairs.empty:
        # Recompute ICC from pairs (quick)
        from digimuh.analysis_00b_stats import compute_stability
        _, icc = compute_stability(bs)

    log.info("Generating figures …")

    plot_grouped_boxplots(bs, d)
    plot_paired_below_above(beh, tests, d)
    plot_paired_rumen_vs_resp(bs, d)
    plot_spearman(spearman, d)
    plot_climate(climate, d)
    plot_predictors(bs, d)
    plot_stability(pairs, icc, d)
    plot_thi_vs_temp_scatter(bs, d)
    plot_bodytemp_vs_resp_scatter(bs, d)
    plot_examples(rumen, resp, bs, d)

    log.info("All figures saved to %s", d)


if __name__ == "__main__":
    main()

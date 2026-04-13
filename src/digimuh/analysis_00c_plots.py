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
    """Side-by-side boxplots: rumen vs respiration breakpoints per year,
    with within-year Wilcoxon signed-rank tests (BH-FDR corrected)."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from scipy.stats import wilcoxon
    from digimuh.analysis_00b_stats import benjamini_hochberg, p_to_stars
    _setup()

    years = sorted(bs["year"].dropna().unique().astype(int))

    for env_col, env_label, fname in [
        ("thi", "THI", "grouped_thi"),
        ("temp", "Barn temperature (°C)", "grouped_temp"),
    ]:
        rumen_col = f"{env_col}_breakpoint"
        resp_col = f"resp_{env_col}_breakpoint"
        rumen_conv_col = f"{env_col}_converged"
        resp_conv_col = f"resp_{env_col}_converged"
        rumen_conv = bs[bs[rumen_conv_col] == True]
        resp_conv = bs[bs[resp_conv_col] == True]

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

        # Within-year paired Wilcoxon: rumen vs resp breakpoints
        # (only for animals with both converged in the same year)
        raw_ps = []
        year_indices = []
        paired_ns = []
        for i, y in enumerate(years):
            paired = bs[
                (bs["year"] == y)
                & (bs[rumen_conv_col] == True)
                & (bs[resp_conv_col] == True)
            ].dropna(subset=[rumen_col, resp_col])
            if len(paired) >= 5:
                try:
                    _, p = wilcoxon(paired[rumen_col], paired[resp_col])
                    raw_ps.append(p)
                    year_indices.append(i)
                    paired_ns.append(len(paired))
                except ValueError:
                    pass

        # BH-FDR correction across years, draw brackets
        if raw_ps:
            adj_ps = benjamini_hochberg(np.array(raw_ps))
            yrange = ax.get_ylim()[1] - ax.get_ylim()[0]
            bracket_base = ax.get_ylim()[1] + yrange * 0.02

            for j, (yi, adj_p, n_paired) in enumerate(
                zip(year_indices, adj_ps, paired_ns)
            ):
                stars = p_to_stars(adj_p)
                x_left = yi - width / 2
                x_right = yi + width / 2
                bracket_y = bracket_base + j * yrange * 0.001
                add_significance_bracket(
                    ax, x_left, x_right, bracket_y, stars, h=0.02)
                # Annotate paired n below bracket
                ax.text((x_left + x_right) / 2, bracket_y - yrange * 0.015,
                        f"n={n_paired}", ha="center", fontsize=7,
                        color="#777777")

            # Expand y-axis for brackets + labels
            ax.set_ylim(ax.get_ylim()[0],
                        bracket_base + yrange * 0.12)

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
    """Diagnostic example panels for each signal/predictor combination.

    For each of the four models (body temp vs THI, body temp vs barn temp,
    resp vs THI, resp vs barn temp), plots three scenarios:

    A) BS converged + Hill converged (both methods agree)
    B) BS failed + Hill converged (Hill rescues the threshold)
    C) Davies n.s. (no threshold exists, models should not be trusted)

    Each panel shows the raw data, the broken-stick fit (if converged),
    the Hill fit with lower bend point, and the Davies p-value.
    """
    import matplotlib.pyplot as plt
    from digimuh.analysis_00b_stats import broken_stick_fit, hill_fit
    _setup()

    configs = [
        ("rumen", rumen, "body_temp", "Rumen temperature (°C)",
         "thi", "thi", "barn_thi", "Barn THI", (45, 80)),
        ("rumen", rumen, "body_temp", "Rumen temperature (°C)",
         "temp", "temp", "barn_temp", "Barn temperature (°C)", (5, 35)),
        ("resp", resp, "resp_rate", "Respiration rate (bpm)",
         "resp_thi", "resp_thi", "barn_thi", "Barn THI", (45, 80)),
        ("resp", resp, "resp_rate", "Respiration rate (bpm)",
         "resp_temp", "resp_temp", "barn_temp", "Barn temperature (°C)", (5, 35)),
    ]

    for signal, data, response_col, ylabel, prefix, bs_prefix, env_col, env_label, x_range in configs:
        if data.empty:
            continue

        conv_col = f"{bs_prefix}_converged"
        davies_col = f"{bs_prefix}_davies_p"
        hill_conv_col = f"{bs_prefix}_hill_converged"
        hill_bend_col = f"{bs_prefix}_hill_bend"
        r2_col = f"{bs_prefix}_r_squared"

        if conv_col not in bs.columns or davies_col not in bs.columns:
            continue

        # ── Select exemplar animals for three scenarios ──────

        # Scenario A: BS converged + Hill converged (best R2)
        mask_a = (bs[conv_col] == True)
        if hill_conv_col in bs.columns:
            mask_a = mask_a & (bs[hill_conv_col] == True)
        candidates_a = bs[mask_a]
        animal_a = None
        if not candidates_a.empty and r2_col in candidates_a.columns:
            valid_a = candidates_a[candidates_a[r2_col].notna()]
            if not valid_a.empty:
                animal_a = valid_a.loc[valid_a[r2_col].idxmax()]

        # Scenario B: BS failed + Hill converged (Hill rescues)
        mask_b = (bs[conv_col] != True) & (bs[davies_col].fillna(1) < 0.05)
        if hill_conv_col in bs.columns:
            mask_b = mask_b & (bs[hill_conv_col] == True) & bs[hill_bend_col].notna()
        candidates_b = bs[mask_b]
        animal_b = None
        if not candidates_b.empty:
            hill_r2_col = f"{bs_prefix}_hill_r2"
            if hill_r2_col in candidates_b.columns:
                valid_b = candidates_b[candidates_b[hill_r2_col].notna()]
                if not valid_b.empty:
                    animal_b = valid_b.loc[valid_b[hill_r2_col].idxmax()]

        # Scenario C: Davies n.s. (no threshold)
        mask_c = bs[davies_col].fillna(1) >= 0.05
        candidates_c = bs[mask_c]
        animal_c = None
        if not candidates_c.empty:
            # Pick one with most data
            n_col = "n_readings" if "resp" not in prefix else "n_resp_readings"
            if n_col in candidates_c.columns:
                valid_c = candidates_c[candidates_c[n_col] > 50]
                if not valid_c.empty:
                    animal_c = valid_c.loc[valid_c[n_col].idxmax()]

        scenarios = [
            ("A", animal_a, "BS + Hill converged"),
            ("B", animal_b, "BS failed, Hill rescues"),
            ("C", animal_c, "Davies n.s. (no threshold)"),
        ]
        valid_scenarios = [(s, a, t) for s, a, t in scenarios if a is not None]

        if not valid_scenarios:
            continue

        fig, axes = plt.subplots(1, len(valid_scenarios),
                                 figsize=(7 * len(valid_scenarios), 6))
        if len(valid_scenarios) == 1:
            axes = [axes]

        for ax, (scenario, animal_row, scenario_title) in zip(axes, valid_scenarios):
            aid = int(animal_row["animal_id"])
            year = int(animal_row["year"])

            grp = data[(data["animal_id"] == aid) & (data["year"] == year)]
            if len(grp) < 30:
                ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes, ha="center")
                continue

            x_vals = grp[env_col].values
            y_vals = grp[response_col].values

            # Scatter
            ax.scatter(x_vals, y_vals, s=2, alpha=0.12, c=COLOURS["identity"])

            xr = np.linspace(np.min(x_vals), np.max(x_vals), 300)

            # Broken-stick fit
            bs_fit = broken_stick_fit(x_vals, y_vals, x_range=x_range)
            if bs_fit.get("converged"):
                bp = bs_fit["breakpoint"]
                yp = np.where(
                    xr <= bp,
                    bs_fit["intercept_below"] + bs_fit["slope_below"] * xr,
                    bs_fit["intercept_above"] + bs_fit["slope_above"] * xr,
                )
                ax.plot(xr, yp, color=COLOURS["below_bp"], linewidth=2,
                        label=f"BS bp={bp:.1f}")
                ax.axvline(bp, color=COLOURS["below_bp"], linestyle="--",
                           linewidth=1, alpha=0.6)

            # Hill fit
            h = hill_fit(x_vals, y_vals, x_range=x_range)
            if h.get("converged"):
                y_min, y_max_fit = h["y_min"], h["y_max"]
                ec50, hill_n = h["ec50"], h["hill_n"]
                ratio = np.clip(ec50 / np.maximum(xr, 1e-10), 1e-10, 1e10)
                yh = y_min + (y_max_fit - y_min) / (1.0 + np.power(ratio, hill_n))
                ax.plot(xr, yh, color=COLOURS["above_bp"], linewidth=2,
                        linestyle="-", label=f"Hill EC50={ec50:.1f}")

                # Lower bend point
                bend = h.get("lower_bend")
                if bend is not None and not np.isnan(bend):
                    ax.axvline(bend, color=COLOURS["above_bp"], linestyle=":",
                               linewidth=1.5, alpha=0.8,
                               label=f"Hill onset={bend:.1f}")

            # Davies p annotation
            davies_p = animal_row.get(davies_col)
            pscore_p = animal_row.get(f"{bs_prefix}_pscore_p")
            anno = []
            if pd.notna(davies_p):
                anno.append(f"Davies p={davies_p:.4f}")
            if pd.notna(pscore_p):
                anno.append(f"Pscore p={pscore_p:.4f}")
            if anno:
                ax.text(0.02, 0.98, "\n".join(anno),
                        transform=ax.transAxes, va="top", fontsize=8,
                        color="#555", fontstyle="italic",
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                                  alpha=0.8, edgecolor="#ccc"))

            ax.set_xlabel(env_label)
            ax.set_ylabel(ylabel)
            ax.set_title(f"({scenario}) {scenario_title}\n"
                         f"Animal {aid} ({year}), n={len(grp):,}",
                         fontsize=10)
            ax.legend(fontsize=8, loc="lower right")

        fig.suptitle(
            f"Diagnostic examples: {ylabel.split('(')[0].strip()} vs {env_label}",
            fontsize=13, fontweight="bold",
        )
        fig.tight_layout()
        _save(fig, f"diagnostic_{prefix}", out_dir)


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
#  « rumen circadian null model »
# ─────────────────────────────────────────────────────────────

def plot_circadian_null_model(out_dir: Path) -> None:
    """Plot 24h rumen temperature profile: cool days vs stress days,
    with breakpoint crossing probability density on secondary y-axis."""
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde
    _setup()

    path = out_dir / "circadian_null_model.csv"
    if not path.exists():
        log.info("  circadian_null_model.csv not found, skipping")
        return

    df = pd.read_csv(path)
    if df.empty:
        return

    log.info("  Plotting circadian null model …")

    fig, ax = plt.subplots(figsize=(10, 6))

    for day_type, colour, label in [
        ("cool", COLOURS["below_bp"], "Cool days (THI < breakpoint all day)"),
        ("stress", COLOURS["above_bp"], "Heat stress days (THI exceeded breakpoint)"),
    ]:
        sub = df[df["day_type"] == day_type]
        if sub.empty:
            continue

        # Grand mean ± SEM across animals at each hour
        hourly = sub.groupby("hour")["body_temp_mean"].agg(["mean", "std", "count"])
        hourly["sem"] = hourly["std"] / np.sqrt(hourly["count"])

        ax.fill_between(hourly.index, hourly["mean"] - hourly["sem"],
                       hourly["mean"] + hourly["sem"],
                       alpha=0.2, color=colour)
        ax.plot(hourly.index, hourly["mean"], color=colour, linewidth=2,
                marker="o", markersize=4, label=label)

    # Mark milking exclusion windows
    for start, end in [(4, 7), (16, 19)]:
        ax.axvspan(start, end, alpha=0.08, color="#999", zorder=0)
        if start == 4:
            ax.text(5.5, ax.get_ylim()[0] + 0.001, "Milking",
                    ha="center", fontsize=8, color="#999", fontstyle="italic")

    ax.set_xlabel("Hour of day")
    ax.set_ylabel("Rumen temperature (°C)")
    ax.set_xticks(range(0, 24, 2))
    ax.set_xlim(-0.5, 23.5)

    # ── Secondary y-axis: crossing probability density ────────
    crossing_path = out_dir / "crossing_times.csv"
    if crossing_path.exists():
        ct = pd.read_csv(crossing_path)
        ct_thi = ct[ct["predictor"] == "thi"] if "predictor" in ct.columns else ct
        if len(ct_thi) > 20:
            ax2 = ax.twinx()

            vals = ct_thi["day_fraction"].dropna().values
            kde = gaussian_kde(vals, bw_method=0.3)
            x_kde = np.linspace(0, 24, 200)
            y_kde = kde(x_kde)

            ax2.fill_between(x_kde, 0, y_kde, alpha=0.12,
                            color="#009E73", zorder=0)
            ax2.plot(x_kde, y_kde, color="#009E73", linewidth=1.5,
                     linestyle="-", alpha=0.7,
                     label=f"Crossing density (n={len(vals)})")
            ax2.set_ylabel("Breakpoint crossing density", color="#009E73")
            ax2.tick_params(axis="y", labelcolor="#009E73")
            ax2.set_ylim(0, y_kde.max() * 2.5)  # leave room above

            # Combine legends
            lines1, labels1 = ax.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax.legend(lines1 + lines2, labels1 + labels2,
                      fontsize=9, loc="upper left")
        else:
            ax.legend(fontsize=9)
    else:
        ax.legend(fontsize=9)

    ax.set_title("Rumen temperature circadian profile\n"
                 "(cool days = null model, stress days = heat effect, "
                 "green = crossing density)")
    fig.tight_layout()
    _save(fig, "circadian_null_model", out_dir)


# ─────────────────────────────────────────────────────────────
#  « THI daily exceedance profile »
# ─────────────────────────────────────────────────────────────

def plot_thi_daily_profile(out_dir: Path) -> None:
    """Plot barn THI across 24h by month, with herd breakpoint line."""
    import matplotlib.pyplot as plt
    _setup()

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
        _save(fig, f"thi_daily_profile_{year}", out_dir)

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
    _save(fig, "thi_daily_profile_all", out_dir)


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
    _setup()

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
        _save(fig, f"crossing_raster_{pred}", out_dir)


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
    _setup()

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
                _save(fig, f"{fname_suffix}_{pred}{variant_suffix}", out_dir)

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
                _save(fig, f"{fname_suffix}_{pred}_overlay{variant_suffix}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « derivative cross-correlation plots »
# ─────────────────────────────────────────────────────────────

def plot_derivative_ccf(out_dir: Path) -> None:
    """Plot derivative CCF: d(climate)/dt vs d(body_temp)/dt."""
    import matplotlib.pyplot as plt
    _setup()

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
        _save(fig, f"dccf_{pred}", out_dir)


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
    _setup()

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

        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # ── Panel A: Climate signal ──────────────────────────
        ax = axes[0]
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
        ax = axes[1]
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
        ax = axes[2]
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

        fig.suptitle(f"Event-triggered average: {pred_label}",
                     fontsize=13, fontweight="bold")
        fig.tight_layout()
        _save(fig, f"eta_{pred}", out_dir)


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
    _setup()

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
    _save(fig, "tnf_yield", out_dir)


# ─────────────────────────────────────────────────────────────
#  « longitudinal breakpoint tracking (repeat animals) »
# ─────────────────────────────────────────────────────────────

def plot_longitudinal_breakpoints(bs: pd.DataFrame, out_dir: Path) -> None:
    """Track how individual breakpoints change across years.

    Only animals present in 2+ years are included.  Shows both absolute
    breakpoints and change relative to the animal's first year.
    """
    import matplotlib.pyplot as plt
    _setup()

    for bp_col, conv_col, env_label, fname in [
        ("thi_breakpoint", "thi_converged", "THI breakpoint", "longitudinal_thi"),
        ("temp_breakpoint", "temp_converged", "Barn temp breakpoint (°C)", "longitudinal_temp"),
    ]:
        conv = bs[bs[conv_col] == True].dropna(subset=[bp_col])
        if conv.empty:
            log.info("  %s: no converged animals, skipping", fname)
            continue

        # Find animals present in 2+ years
        year_counts = conv.groupby("animal_id")["year"].nunique()
        repeat_ids = year_counts[year_counts >= 2].index
        log.info("  %s: %d animals with 2+ years (of %d converged)",
                 fname, len(repeat_ids), len(conv["animal_id"].unique()))
        if len(repeat_ids) < 3:
            log.info("  %s: fewer than 3 repeat animals, skipping", fname)
            continue

        repeat = conv[conv["animal_id"].isin(repeat_ids)].copy()
        repeat = repeat.sort_values(["animal_id", "year"])

        # Compute relative change from first year
        first_bp = repeat.groupby("animal_id").first()[bp_col]
        repeat["first_bp"] = repeat["animal_id"].map(first_bp)
        repeat["bp_change"] = repeat[bp_col] - repeat["first_bp"]

        years = sorted(repeat["year"].unique().astype(int))

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # ── Panel A: Absolute breakpoints ────────────────────
        ax = axes[0]
        # Individual trajectories
        for aid in repeat_ids:
            animal = repeat[repeat["animal_id"] == aid].sort_values("year")
            ax.plot(animal["year"], animal[bp_col],
                    color="#aaaaaa", linewidth=0.5, alpha=0.4, zorder=1)
            ax.scatter(animal["year"], animal[bp_col],
                      s=8, color="#888888", alpha=0.3, zorder=2)

        # Boxplot overlay
        data_by_year = [repeat[repeat["year"] == y][bp_col].values for y in years]
        bp_plot = ax.boxplot(
            data_by_year, positions=years, widths=0.5,
            patch_artist=True,
            boxprops=dict(facecolor=COLOURS["below_bp"], alpha=0.5, edgecolor="#333"),
            medianprops=dict(color="#333", linewidth=2),
            manage_ticks=False, zorder=3,
        )
        ax.set_xlabel("Year")
        ax.set_ylabel(env_label)
        ax.set_title(f"Absolute breakpoints\n(n={len(repeat_ids)} animals, 2+ years)")
        ax.set_xticks(years)

        # Significance brackets: consecutive years (Fisher resampling, BH-FDR)
        try:
            from rerandomstats import FisherResamplingTest
            from digimuh.analysis_00b_stats import benjamini_hochberg, p_to_stars
            raw_ps_abs = []
            year_pairs_abs = []
            for j in range(len(years) - 1):
                y1, y2 = years[j], years[j + 1]
                ids_both = set(repeat[repeat["year"] == y1]["animal_id"]) & \
                           set(repeat[repeat["year"] == y2]["animal_id"])
                if len(ids_both) >= 5:
                    d1 = repeat[(repeat["year"] == y1) & (repeat["animal_id"].isin(ids_both))][bp_col].tolist()
                    d2 = repeat[(repeat["year"] == y2) & (repeat["animal_id"].isin(ids_both))][bp_col].tolist()
                    p = FisherResamplingTest(data_a=d1, data_b=d2,
                                            func="medianDiff", combination_n=20_000).main()
                    raw_ps_abs.append(p)
                    year_pairs_abs.append((y1, y2))
            if raw_ps_abs:
                adj_ps_abs = benjamini_hochberg(np.array(raw_ps_abs))
                ymax = ax.get_ylim()[1]
                step = (ymax - ax.get_ylim()[0]) * 0.05
                for k, ((y1, y2), adj_p) in enumerate(zip(year_pairs_abs, adj_ps_abs)):
                    stars = p_to_stars(adj_p)
                    add_significance_bracket(ax, y1, y2, ymax + step * k, stars)
                ax.set_ylim(ax.get_ylim()[0], ymax + step * (len(raw_ps_abs) + 1))
        except (ImportError, Exception) as e:
            log.debug("  Skipping longitudinal brackets: %s", e)

        # ── Panel B: Relative change from first year ─────────
        ax = axes[1]
        for aid in repeat_ids:
            animal = repeat[repeat["animal_id"] == aid].sort_values("year")
            ax.plot(animal["year"], animal["bp_change"],
                    color="#aaaaaa", linewidth=0.5, alpha=0.4, zorder=1)
            ax.scatter(animal["year"], animal["bp_change"],
                      s=8, color="#888888", alpha=0.3, zorder=2)

        data_by_year_rel = [repeat[repeat["year"] == y]["bp_change"].values for y in years]
        bp_plot = ax.boxplot(
            data_by_year_rel, positions=years, widths=0.5,
            patch_artist=True,
            boxprops=dict(facecolor=COLOURS["above_bp"], alpha=0.5, edgecolor="#333"),
            medianprops=dict(color="#333", linewidth=2),
            manage_ticks=False, zorder=3,
        )
        ax.axhline(0, color="#999", linestyle="--", linewidth=1)
        ax.set_xlabel("Year")
        ax.set_ylabel("Change from first year")

        # Significance: test if each year's change differs from 0
        try:
            from rerandomstats import FisherResamplingTest
            from digimuh.analysis_00b_stats import p_to_stars
            anno_lines = []
            for y in years[1:]:  # first year is always 0
                changes = repeat[repeat["year"] == y]["bp_change"].dropna()
                if len(changes) >= 5:
                    p = FisherResamplingTest(
                        data_a=changes.tolist(),
                        data_b=[0.0] * len(changes),
                        func="medianDiff", combination_n=20_000).main()
                    stars = p_to_stars(p)
                    if stars != "n.s.":
                        anno_lines.append(f"{y}: {stars}")
                    ax.text(y, changes.median(), stars,
                            ha="center", va="bottom", fontsize=8,
                            fontweight="bold", color="#333", zorder=5)
        except (ImportError, Exception):
            pass

        ax.set_title("Relative change\n(baseline = animal's first summer)")
        ax.set_xticks(years)

        fig.suptitle(f"Longitudinal {env_label} tracking",
                     fontsize=13, fontweight="bold")
        fig.tight_layout()
        _save(fig, fname, out_dir)


# ─────────────────────────────────────────────────────────────
#  « threshold detection pipeline Sankey (plotly) »
# ─────────────────────────────────────────────────────────────

def plot_threshold_sankey(bs: pd.DataFrame, out_dir: Path) -> None:
    """Plotly Sankey diagrams showing how animals flow through the
    threshold detection pipeline: Davies/pscore → broken-stick → Hill.

    Produces one figure for rumen temperature vs THI and one for
    respiration rate vs THI.

    Args:
        bs: broken_stick_results.csv DataFrame.
        out_dir: Output directory.
    """
    import plotly.graph_objects as go

    out_dir.mkdir(parents=True, exist_ok=True)

    configs = [
        ("thi", "rumen", "Rumen temperature vs. barn THI",
         lambda df: df),
        ("resp_thi", "respiration", "Respiration rate vs. barn THI",
         lambda df: df[df["n_resp_readings"] > 0]),
    ]

    for prefix, signal_name, title, filter_fn in configs:
        davies_col = f"{prefix}_davies_p"
        pscore_col = f"{prefix}_pscore_p"
        conv_col = f"{prefix}_converged"
        hill_conv_col = f"{prefix}_hill_converged"
        hill_bend_col = f"{prefix}_hill_bend"

        if davies_col not in bs.columns:
            continue

        subset = filter_fn(bs)
        if subset.empty or len(subset) < 5:
            continue

        n_total = len(subset)

        # Stage 1: Davies
        n_davies_sig = (subset[davies_col].fillna(1) < 0.05).sum()
        n_davies_ns = n_total - n_davies_sig

        # Stage 1b: Pscore (parallel)
        n_pscore_sig = (subset[pscore_col].fillna(1) < 0.05).sum()

        # Stage 2: Broken-stick (from Davies sig)
        davies_mask = subset[davies_col].fillna(1) < 0.05
        n_bs_conv = (davies_mask & (subset[conv_col] == True)).sum()
        n_bs_fail = n_davies_sig - n_bs_conv

        # Stage 3: Hill (from BS-failed)
        bs_fail_mask = davies_mask & (subset[conv_col] != True)
        if hill_conv_col in subset.columns:
            n_hill_ok = (
                bs_fail_mask
                & (subset[hill_conv_col] == True)
                & subset[hill_bend_col].notna()
            ).sum()
        else:
            n_hill_ok = 0
        n_hill_fail = n_bs_fail - n_hill_ok

        # Final tallies
        n_threshold = n_bs_conv + n_hill_ok
        n_no_threshold = n_davies_ns + n_hill_fail

        # ── Node definitions ─────────────────────────────────
        # 0: Total
        # 1: Davies significant (nonlinear)
        # 2: Davies n.s. (linear)
        # 3: BS converged
        # 4: BS failed
        # 5: Hill onset found
        # 6: Hill / no threshold
        # 7: Threshold identified (final)
        # 8: No threshold (final)

        node_labels = [
            f"All animals<br>n = {n_total}",
            f"Nonlinear (Davies p<0.05)<br>n = {n_davies_sig}"
            f"<br><i>Pscore: {n_pscore_sig}/{n_total}</i>",
            f"Linear (no breakpoint)<br>n = {n_davies_ns}",
            f"Broken-stick converged<br>n = {n_bs_conv}",
            f"Broken-stick failed<br>n = {n_bs_fail}",
            f"Hill onset found<br>n = {n_hill_ok}",
            f"Hill failed<br>n = {n_hill_fail}",
            f"Threshold identified<br>n = {n_threshold}",
            f"No threshold<br>n = {n_no_threshold}",
        ]

        # Wong palette as rgba
        node_colors = [
            "rgba(86, 180, 233, 0.9)",   # 0 sky blue (total)
            "rgba(0, 114, 178, 0.9)",     # 1 blue (Davies sig)
            "rgba(153, 153, 153, 0.9)",   # 2 grey (linear)
            "rgba(0, 114, 178, 0.9)",     # 3 blue (BS ok)
            "rgba(230, 159, 0, 0.9)",     # 4 amber (BS fail)
            "rgba(213, 94, 0, 0.9)",      # 5 vermilion (Hill ok)
            "rgba(204, 121, 167, 0.9)",   # 6 purple (Hill fail)
            "rgba(0, 158, 115, 0.9)",     # 7 green (threshold)
            "rgba(153, 153, 153, 0.9)",   # 8 grey (no threshold)
        ]

        # ── Link definitions (source → target, value) ────────
        links_source = []
        links_target = []
        links_value = []
        links_color = []

        def add_link(src, tgt, val, color):
            if val > 0:
                links_source.append(src)
                links_target.append(tgt)
                links_value.append(val)
                links_color.append(color)

        # Total → Davies sig / n.s.
        add_link(0, 1, n_davies_sig, "rgba(0, 114, 178, 0.25)")
        add_link(0, 2, n_davies_ns,  "rgba(153, 153, 153, 0.25)")

        # Davies sig → BS converged / failed
        add_link(1, 3, n_bs_conv, "rgba(0, 114, 178, 0.25)")
        add_link(1, 4, n_bs_fail, "rgba(230, 159, 0, 0.25)")

        # BS failed → Hill ok / fail
        add_link(4, 5, n_hill_ok,   "rgba(213, 94, 0, 0.25)")
        add_link(4, 6, n_hill_fail, "rgba(204, 121, 167, 0.25)")

        # → Final: threshold identified
        add_link(3, 7, n_bs_conv, "rgba(0, 158, 115, 0.25)")
        add_link(5, 7, n_hill_ok, "rgba(0, 158, 115, 0.25)")

        # → Final: no threshold
        add_link(2, 8, n_davies_ns,  "rgba(153, 153, 153, 0.20)")
        add_link(6, 8, n_hill_fail,  "rgba(153, 153, 153, 0.20)")

        fig = go.Figure(data=[go.Sankey(
            arrangement="snap",
            node=dict(
                pad=25,
                thickness=30,
                line=dict(color="#333333", width=1),
                label=node_labels,
                color=node_colors,
            ),
            link=dict(
                source=links_source,
                target=links_target,
                value=links_value,
                color=links_color,
            ),
        )])

        fig.update_layout(
            title=dict(
                text=title,
                font=dict(size=16, family="Arial"),
                x=0.5,
            ),
            font=dict(size=12, family="Arial"),
            width=1000,
            height=500,
            margin=dict(l=20, r=20, t=60, b=20),
        )

        # Save as interactive HTML (skip SVG/PNG: kaleido+chromium hangs)
        fig.write_html(str(out_dir / f"sankey_{prefix}.html"),
                        include_plotlyjs="cdn")
        log.info("  Saved sankey_%s.html", prefix)


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

    _plot_calls = [
        ("Grouped boxplots", lambda: plot_grouped_boxplots(bs, d)),
        ("Paired below/above", lambda: plot_paired_below_above(beh, tests, d)),
        ("Paired rumen vs resp", lambda: plot_paired_rumen_vs_resp(bs, d)),
        ("Spearman histograms", lambda: plot_spearman(spearman, d)),
        ("Climate time series", lambda: plot_climate(climate, d)),
        ("Predictors", lambda: plot_predictors(bs, d)),
        ("Stability scatter", lambda: plot_stability(pairs, icc, d)),
        ("THI vs temp scatter", lambda: plot_thi_vs_temp_scatter(bs, d)),
        ("Body temp vs resp scatter", lambda: plot_bodytemp_vs_resp_scatter(bs, d)),
        ("Diagnostic examples", lambda: plot_examples(rumen, resp, bs, d)),
        ("Cross-correlation", lambda: plot_cross_correlation(d)),
        ("Circadian null model", lambda: plot_circadian_null_model(d)),
        ("THI daily profile", lambda: plot_thi_daily_profile(d)),
        ("Crossing raster", lambda: plot_crossing_raster(d)),
        ("Derivative CCF", lambda: plot_derivative_ccf(d)),
        ("Event-triggered average", lambda: plot_event_triggered_average(d)),
        ("Event-triggered average (8-11h)", lambda: plot_event_triggered_average(
            d, traces_file="event_triggered_traces_filtered.csv",
            suffix="_8to11h", title_extra=" (crossings 8:00–11:00 only)")),
        ("TNF vs yield", lambda: plot_tnf_yield(d)),
        ("Longitudinal breakpoints", lambda: plot_longitudinal_breakpoints(bs, d)),
        ("Sankey diagrams", lambda: plot_threshold_sankey(bs, d)),
    ]

    for name, fn in _plot_calls:
        try:
            fn()
        except Exception as e:
            log.warning("  %s FAILED: %s", name, e)

    log.info("All figures saved to %s", d)


if __name__ == "__main__":
    main()

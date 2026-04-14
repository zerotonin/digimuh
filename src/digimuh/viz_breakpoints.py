#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — viz_breakpoints                                      ║
# ║  « breakpoint distribution, comparison, and diagnostic plots » ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Grouped boxplots, paired below/above, rumen vs respiration,    ║
# ║  Spearman histograms, climate time series, predictor panels,    ║
# ║  stability scatter, diagnostic examples, and scatter plots.     ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Breakpoint analysis figures for the Frontiers manuscript."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.constants import COLOURS
from digimuh.viz_base import setup_figure, save_figure, add_significance_bracket

log = logging.getLogger("digimuh.viz")

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
    setup_figure()

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

        save_figure(fig, f"boxplot_{fname}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « paired rumen vs resp with significance brackets »
# ─────────────────────────────────────────────────────────────

def plot_paired_below_above(
    beh: pd.DataFrame, tests: pd.DataFrame, out_dir: Path,
) -> None:
    """Paired boxplots: below vs above breakpoint with significance brackets."""
    import matplotlib.pyplot as plt
    from scipy.stats import wilcoxon
    setup_figure()

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
        save_figure(fig, f"paired_below_above_{year}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « paired rumen vs respiration breakpoints »
# ─────────────────────────────────────────────────────────────

def plot_paired_rumen_vs_resp(bs: pd.DataFrame, out_dir: Path) -> None:
    """Paired boxplot: rumen temp breakpoint vs respiration breakpoint."""
    import matplotlib.pyplot as plt
    from scipy.stats import wilcoxon
    setup_figure()

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
    save_figure(fig, "paired_rumen_vs_resp", out_dir)


# ─────────────────────────────────────────────────────────────
#  « Spearman histograms »
# ─────────────────────────────────────────────────────────────

def plot_spearman(spearman: pd.DataFrame, out_dir: Path) -> None:
    """Spearman correlation distributions."""
    import matplotlib.pyplot as plt
    setup_figure()

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
        ax.set_xlabel(f"Spearman $r_s$ ({label})")
        ax.set_ylabel("Number of animals")
        ax.legend(fontsize=9)

    fig.suptitle("Per-animal Spearman correlations")
    fig.tight_layout()
    save_figure(fig, "spearman_correlations", out_dir)


# ─────────────────────────────────────────────────────────────
#  « climate time series »
# ─────────────────────────────────────────────────────────────

def plot_climate(climate: pd.DataFrame, out_dir: Path) -> None:
    """Daily barn THI and temperature time series per summer."""
    import matplotlib.pyplot as plt
    setup_figure()

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
        save_figure(fig, fname, out_dir)


# ─────────────────────────────────────────────────────────────
#  « predictors: breakpoint vs milk yield / lactation »
# ─────────────────────────────────────────────────────────────

def plot_predictors(bs: pd.DataFrame, out_dir: Path) -> None:
    """Scatter plots: breakpoint vs production parameters."""
    import matplotlib.pyplot as plt
    from scipy.stats import pearsonr
    setup_figure()

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
        save_figure(fig, fname, out_dir)


# ─────────────────────────────────────────────────────────────
#  « breakpoint stability scatter »
# ─────────────────────────────────────────────────────────────

def plot_stability(pairs: pd.DataFrame, icc: float, out_dir: Path) -> None:
    """Year-to-year breakpoint stability scatter."""
    import matplotlib.pyplot as plt
    setup_figure()

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
    save_figure(fig, "breakpoint_stability", out_dir)


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
    setup_figure()

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
        save_figure(fig, f"diagnostic_{prefix}", out_dir)


# ─────────────────────────────────────────────────────────────
#  « body temp vs resp breakpoint scatter »
# ─────────────────────────────────────────────────────────────

def plot_bodytemp_vs_resp_scatter(bs: pd.DataFrame, out_dir: Path) -> None:
    """Scatter: rumen temp THI breakpoint vs resp THI breakpoint."""
    import matplotlib.pyplot as plt
    setup_figure()

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
    save_figure(fig, "scatter_bodytemp_vs_resp_bp", out_dir)


# ─────────────────────────────────────────────────────────────
#  « THI vs barn temp breakpoint scatter »
# ─────────────────────────────────────────────────────────────

def plot_thi_vs_temp_scatter(bs: pd.DataFrame, out_dir: Path) -> None:
    """Scatter: THI breakpoint vs barn temp breakpoint."""
    import matplotlib.pyplot as plt
    setup_figure()

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
    save_figure(fig, "scatter_thi_vs_temp", out_dir)


# ─────────────────────────────────────────────────────────────
#  « rumen circadian null model »
# ─────────────────────────────────────────────────────────────


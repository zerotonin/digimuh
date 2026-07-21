#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — viz_longitudinal                                     ║
# ║  « breakpoint trajectories, Sankey adaptation, crossing count »║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Longitudinal breakpoint tracking: individual trajectories,     ║
# ║  crossing count rainclouds, alluvial Sankey for adaptation,     ║
# ║  and threshold detection pipeline Sankey.                       ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Longitudinal and pipeline-summary figures."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.constants import COLOURS, RESAMPLING_SEED
from digimuh.paths import resolve_input, resolve_output
from digimuh.viz_base import add_significance_bracket, save_figure, setup_figure

log = logging.getLogger("digimuh.viz")

# ─────────────────────────────────────────────────────────────
#  « longitudinal breakpoint tracking (repeat animals) »
# ─────────────────────────────────────────────────────────────

def plot_longitudinal_breakpoints(bs: pd.DataFrame, out_dir: Path) -> None:
    """Track how individual breakpoints change across years.

    Only animals present in 2+ years are included.  Shows both absolute
    breakpoints and change relative to the animal's first year.
    """
    import matplotlib.pyplot as plt
    setup_figure()

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
        ax.boxplot(
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
            from rerandomstats import FisherResamplingTest, correct_pvalues_array

            from digimuh.stats_core import p_to_stars
            raw_ps_abs = []
            year_pairs_abs = []
            for j in range(len(years) - 1):
                y1, y2 = years[j], years[j + 1]
                ids_both = set(repeat[repeat["year"] == y1]["animal_id"]) & \
                           set(repeat[repeat["year"] == y2]["animal_id"])
                if len(ids_both) >= 5:
                    d1 = repeat[
                        (repeat["year"] == y1) & (repeat["animal_id"].isin(ids_both))
                    ][bp_col].tolist()
                    d2 = repeat[
                        (repeat["year"] == y2) & (repeat["animal_id"].isin(ids_both))
                    ][bp_col].tolist()
                    p = FisherResamplingTest(
                        data_a=d1, data_b=d2, func="medianDiff",
                        combination_n=20_000, seed=RESAMPLING_SEED).main()
                    raw_ps_abs.append(p)
                    year_pairs_abs.append((y1, y2))
            if raw_ps_abs:
                adj_ps_abs = correct_pvalues_array(np.array(raw_ps_abs), method="fdr_bh")
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
        ax.boxplot(
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

            from digimuh.stats_core import p_to_stars
            anno_lines = []
            for y in years[1:]:  # first year is always 0
                changes = repeat[repeat["year"] == y]["bp_change"].dropna()
                if len(changes) >= 5:
                    p = FisherResamplingTest(
                        data_a=changes.tolist(),
                        data_b=[0.0] * len(changes),
                        func="medianDiff", combination_n=20_000, seed=RESAMPLING_SEED).main()
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
        save_figure(fig, fname, out_dir)


# ─────────────────────────────────────────────────────────────
#  « post-hoc helpers for the per-year rainclouds »
# ─────────────────────────────────────────────────────────────

def _posthoc_pairwise(year_data: list, years: list,
                      *, test: str = "Fisher:medianDiff",
                      combination_n: int = 20_000):
    """All pairwise year comparisons with BH-FDR via MultiGroupTest.

    The omnibus Kruskal-Wallis only says *some* year differs; this is the
    follow-up that says *which* pairs.  Returns the reRandomStats result
    frame (groupA/groupB, n, raw + corrected p, ``h``, ``sig. level``),
    or ``None`` when fewer than two years carry enough data.
    """
    from rerandomstats import MultiGroupTest
    flat: list[float] = []
    labels: list[str] = []
    for y, vals in zip(years, year_data):
        v = np.asarray(vals, dtype=float)
        v = v[np.isfinite(v)]
        if v.size >= 3:
            flat.extend(v.tolist())
            labels.extend([str(int(y))] * v.size)
    if len(set(labels)) < 2:
        return None
    return MultiGroupTest(data=flat, group=labels, test=test,
                          combination_n=combination_n,
                          correction_type="fdr_bh", seed=RESAMPLING_SEED).main()


def _compact_letters(years: list, posthoc) -> dict[str, str]:
    """Compact letter display from a MultiGroupTest result frame.

    Groups sharing a letter are NOT significantly different after BH-FDR.
    Built as a maximal-clique cover of the 'not different' graph (cheap
    brute force; the number of years is small).
    """
    import itertools
    import string

    labels = [str(int(y)) for y in years]
    notdiff = {frozenset((a, b)) for a, b in itertools.combinations(labels, 2)}
    for _, row in posthoc.iterrows():
        if bool(row["h"]):
            notdiff.discard(frozenset((str(row["groupA"]), str(row["groupB"]))))

    def is_clique(members) -> bool:
        return all(frozenset((a, b)) in notdiff
                   for a, b in itertools.combinations(members, 2))

    cliques: list[set] = []
    for size in range(len(labels), 0, -1):
        for combo in itertools.combinations(labels, size):
            members = set(combo)
            if is_clique(members) and not any(members <= c for c in cliques):
                cliques.append(members)

    letters = {lab: "" for lab in labels}
    for i, clique in enumerate(cliques):
        for lab in clique:
            letters[lab] += string.ascii_lowercase[i]
    return {lab: "".join(sorted(letters[lab])) for lab in labels}


def _annotate_posthoc(ax, year_data: list, years: list, fname: str,
                      out_dir: Path) -> None:
    """Run pairwise post-hoc, draw a compact letter display, write CSV.

    Shared by both per-year rainclouds.  Years sharing a letter are not
    significantly different (Fisher resampling, BH-FDR).
    """
    posthoc = _posthoc_pairwise(year_data, years)
    if posthoc is None:
        return
    cld = _compact_letters(years, posthoc)
    for i, y in enumerate(years):
        ax.text(1.01, i, cld.get(str(int(y)), ""),
                transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=11, fontweight="bold", color="#333")
    ax.text(0.01, 0.02,
            "letters: post-hoc groups (Fisher resampling, BH-FDR; "
            "shared letter = n.s.)",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=7, color="#666")
    posthoc.to_csv(resolve_output(out_dir, f"{fname}_posthoc.csv"), index=False)


def _write_year_summary(year_data: list, years: list, fname: str,
                        out_dir: Path) -> None:
    """Write per-year and overall (all years pooled) n / median / IQR.

    Lands as ``<fname>_summary.csv`` next to the figure — the pooled
    ``overall`` row carries the cross-year median for the record (it is
    deliberately not drawn on the figure).
    """
    rows = []
    pooled = []
    for y, vals in zip(years, year_data):
        v = np.asarray(vals, dtype=float)
        v = v[np.isfinite(v)]
        pooled.append(v)
        if v.size:
            rows.append({"year": int(y), "n": int(v.size),
                         "median": float(np.median(v)),
                         "q25": float(np.percentile(v, 25)),
                         "q75": float(np.percentile(v, 75)),
                         "min": float(v.min()),
                         "max": float(v.max())})
    allv = np.concatenate(pooled) if pooled else np.array([])
    if allv.size:
        rows.append({"year": "overall", "n": int(allv.size),
                     "median": float(np.median(allv)),
                     "q25": float(np.percentile(allv, 25)),
                     "q75": float(np.percentile(allv, 75)),
                     "min": float(allv.min()),
                     "max": float(allv.max())})
    pd.DataFrame(rows).to_csv(
        resolve_output(out_dir, f"{fname}_summary.csv"), index=False)


# ─────────────────────────────────────────────────────────────
#  « breakpoint crossing count raincloud per year »
# ─────────────────────────────────────────────────────────────

def plot_breakpoint_raincloud(out_dir: Path) -> None:
    """Raincloud plots of annual breakpoint crossing counts per animal.

    For each animal-year, counts how many times the barn THI (or barn
    temp) crossed that animal's individual breakpoint upward.  Reads
    crossing events from ``crossing_times.csv`` (produced by stats).

    Each year gets a horizontal raincloud: half-violin (KDE), jittered
    scatter of integer counts, and boxplot.  Shows whether cows
    experience more threshold exceedances over time (climate signal).
    """
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde, kruskal
    setup_figure()

    path = resolve_input(out_dir, "crossing_times.csv")
    if not path.exists():
        log.info("  crossing_times.csv not found, skipping crossing count raincloud")
        return

    ct = pd.read_csv(path)
    if ct.empty:
        return

    for pred, pred_label, fname in [
        ("thi", "THI breakpoint", "raincloud_crossing_count_thi"),
        ("temp", "Barn temp breakpoint", "raincloud_crossing_count_temp"),
    ]:
        sub = ct[ct["predictor"] == pred]
        if sub.empty:
            continue

        # Count crossings per animal per year
        counts = (sub.groupby(["animal_id", "year"])
                     .size()
                     .reset_index(name="n_crossings"))

        years = sorted(counts["year"].unique().astype(int))
        if len(years) < 2:
            continue

        year_data = [counts[counts["year"] == y]["n_crossings"].values
                     for y in years]

        fig, ax = plt.subplots(figsize=(10, 1.5 + 1.4 * len(years)))

        for i, (y, vals) in enumerate(zip(years, year_data)):
            y_pos = i
            colour = COLOURS["year"].get(y, COLOURS["below_bp"])

            # Half-violin (KDE) — above row centre
            if len(vals) > 5:
                kde = gaussian_kde(vals, bw_method=0.3)
                x_kde = np.linspace(
                    max(0, vals.min() - 5), vals.max() + 5, 200)
                density = kde(x_kde)
                density_scaled = density / density.max() * 0.38
                ax.fill_between(x_kde, y_pos, y_pos + density_scaled,
                                alpha=0.3, color=colour)
                ax.plot(x_kde, y_pos + density_scaled,
                        color=colour, linewidth=1.2)

            # Jittered scatter — below row centre
            jitter = np.random.uniform(-0.15, -0.35, len(vals))
            ax.scatter(vals, y_pos + jitter, s=10, alpha=0.5,
                       color=colour, edgecolors="none", zorder=3)

            # Boxplot — at row centre
            ax.boxplot(
                vals, positions=[y_pos - 0.02], widths=0.12,
                vert=False, patch_artist=True,
                boxprops=dict(facecolor=colour, alpha=0.5,
                              edgecolor="#333"),
                medianprops=dict(color="#333", linewidth=2),
                whiskerprops=dict(color="#333"),
                capprops=dict(color="#333"),
                flierprops=dict(marker="o", markersize=2, alpha=0.3),
                manage_ticks=False,
            )

            # n animals + median label
            n_animals = len(vals)
            med = int(np.median(vals))
            ax.text(ax.get_xlim()[0] if ax.get_xlim()[0] != 0 else
                    max(0, vals.min() - 8),
                    y_pos - 0.38,
                    f"n={n_animals}, median={med}",
                    fontsize=8, color="#666", va="center")

        # Kruskal-Wallis across years
        valid_groups = [v for v in year_data if len(v) >= 3]
        if len(valid_groups) >= 2:
            try:
                stat, p = kruskal(*valid_groups)
                ax.text(0.99, 0.02,
                        f"Kruskal-Wallis H={stat:.1f}, p={p:.3g}",
                        transform=ax.transAxes, ha="right",
                        va="bottom", fontsize=9, color="#333",
                        bbox=dict(boxstyle="round,pad=0.3",
                                  facecolor="white", alpha=0.8))
            except ValueError:
                pass

        # Post-hoc: pairwise years (Fisher resampling, BH-FDR) + CLD
        _annotate_posthoc(ax, year_data, years,
                          f"raincloud_crossing_count_{pred}", out_dir)
        _write_year_summary(year_data, years,
                            f"raincloud_crossing_count_{pred}", out_dir)

        ax.set_yticks(range(len(years)))
        ax.set_yticklabels([str(y) for y in years], fontsize=11)
        ax.set_xlabel(f"Number of {pred_label} crossings per animal")
        ax.set_title(
            f"Annual {pred_label} crossing count per animal",
            fontsize=13, fontweight="bold")
        ax.set_xlim(left=0)
        ax.invert_yaxis()
        fig.tight_layout()
        save_figure(fig, fname, out_dir)


# ─────────────────────────────────────────────────────────────
#  « breakpoint value raincloud per year »
# ─────────────────────────────────────────────────────────────

def plot_breakpoint_value_raincloud(bs: pd.DataFrame, out_dir: Path) -> None:
    """Raincloud plots of the individual breakpoint *value* per year.

    Companion to :func:`plot_breakpoint_raincloud`, which shows crossing
    *counts*.  Here each converged animal contributes its fitted THI (or
    barn-temp) breakpoint for the year, so the rows show how the herd's
    threshold distribution shifts across 2021-2024.  Same per-year colour
    code and Kruskal-Wallis test as the crossing-count raincloud.
    """
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde, kruskal
    setup_figure()

    if bs.empty:
        log.info("  broken_stick results empty, skipping breakpoint value raincloud")
        return

    for bp_col, conv_col, pred_label, fname in [
        ("thi_breakpoint", "thi_converged", "THI breakpoint",
         "raincloud_breakpoint_value_thi"),
        ("temp_breakpoint", "temp_converged", "Barn temp breakpoint (°C)",
         "raincloud_breakpoint_value_temp"),
    ]:
        conv = bs[bs[conv_col] == True].dropna(subset=[bp_col])
        if conv.empty:
            continue

        years = sorted(conv["year"].unique().astype(int))
        if len(years) < 2:
            continue

        year_data = [conv[conv["year"] == y][bp_col].to_numpy() for y in years]

        all_vals = np.concatenate(year_data)
        span = all_vals.max() - all_vals.min()
        pad = 0.05 * span if span > 0 else 1.0
        x_lo, x_hi = all_vals.min() - pad, all_vals.max() + pad

        fig, ax = plt.subplots(figsize=(10, 1.5 + 1.4 * len(years)))

        for i, (y, vals) in enumerate(zip(years, year_data)):
            y_pos = i
            colour = COLOURS["year"].get(y, COLOURS["below_bp"])

            # Half-violin (KDE) — above row centre
            if len(vals) > 5:
                kde = gaussian_kde(vals)
                x_kde = np.linspace(x_lo, x_hi, 200)
                density = kde(x_kde)
                density_scaled = density / density.max() * 0.38
                ax.fill_between(x_kde, y_pos, y_pos + density_scaled,
                                alpha=0.3, color=colour)
                ax.plot(x_kde, y_pos + density_scaled,
                        color=colour, linewidth=1.2)

            # Jittered scatter — below row centre
            jitter = np.random.uniform(-0.15, -0.35, len(vals))
            ax.scatter(vals, y_pos + jitter, s=10, alpha=0.5,
                       color=colour, edgecolors="none", zorder=3)

            # Boxplot — at row centre
            ax.boxplot(
                vals, positions=[y_pos - 0.02], widths=0.12,
                vert=False, patch_artist=True,
                boxprops=dict(facecolor=colour, alpha=0.5, edgecolor="#333"),
                medianprops=dict(color="#333", linewidth=2),
                whiskerprops=dict(color="#333"),
                capprops=dict(color="#333"),
                flierprops=dict(marker="o", markersize=2, alpha=0.3),
                manage_ticks=False,
            )

            # n animals + median label
            ax.text(x_lo, y_pos - 0.38,
                    f"n={len(vals)}, median={np.median(vals):.1f}",
                    fontsize=8, color="#666", va="center")

        # Kruskal-Wallis across years (same test as the crossing raincloud)
        valid_groups = [v for v in year_data if len(v) >= 3]
        if len(valid_groups) >= 2:
            try:
                stat, p = kruskal(*valid_groups)
                ax.text(0.99, 0.02,
                        f"Kruskal-Wallis H={stat:.1f}, p={p:.3g}",
                        transform=ax.transAxes, ha="right", va="bottom",
                        fontsize=9, color="#333",
                        bbox=dict(boxstyle="round,pad=0.3",
                                  facecolor="white", alpha=0.8))
            except ValueError:
                pass

        # Post-hoc: pairwise years (Fisher resampling, BH-FDR) + CLD
        _annotate_posthoc(ax, year_data, years, fname, out_dir)
        _write_year_summary(year_data, years, fname, out_dir)

        ax.set_yticks(range(len(years)))
        ax.set_yticklabels([str(y) for y in years], fontsize=11)
        ax.set_xlabel(f"{pred_label} per animal")
        ax.set_title(f"Annual {pred_label} per animal", fontsize=13,
                     fontweight="bold")
        ax.set_xlim(x_lo, x_hi)
        ax.invert_yaxis()
        fig.tight_layout()
        save_figure(fig, fname, out_dir)

        # CSV companion (one row per animal-year) for reproducibility
        conv[["animal_id", "year", bp_col]].to_csv(
            resolve_output(out_dir, f"{fname}.csv"), index=False)


# ─────────────────────────────────────────────────────────────
#  « year-to-year predictability (test-retest scatter) »
# ─────────────────────────────────────────────────────────────

def plot_breakpoint_retest(bs: pd.DataFrame, out_dir: Path) -> None:
    """This-summer vs last-summer breakpoint scatter (predictability).

    Third companion to the per-year value and crossing-count rainclouds:
    each point is one cow's consecutive-summer pair (x = last summer,
    y = this summer), coloured by the later year.  The identity line is
    perfect stability; the flat fitted line and near-zero r show the
    breakpoint does not carry over between years.
    """
    import matplotlib.pyplot as plt
    from scipy import stats
    setup_figure()

    for bp_col, conv_col, label, unit, fname in [
        ("thi_breakpoint", "thi_converged", "THI breakpoint", "",
         "retest_breakpoint_thi"),
        ("temp_breakpoint", "temp_converged", "Barn-temp breakpoint", " °C",
         "retest_breakpoint_temp"),
    ]:
        d = bs[bs[conv_col] == True].dropna(subset=[bp_col])[
            ["animal_id", "year", bp_col]].copy()
        d["year"] = d["year"].astype(int)

        rows = []
        for aid, g in d.groupby("animal_id"):
            m = dict(zip(g["year"], g[bp_col]))
            for y in m:
                if y + 1 in m:
                    rows.append((aid, y, y + 1, m[y], m[y + 1]))
        if len(rows) < 5:
            continue
        pairs = pd.DataFrame(rows, columns=["animal_id", "from_year",
                                            "to_year", "bp_last", "bp_this"])
        x = pairs["bp_last"].to_numpy()
        y = pairs["bp_this"].to_numpy()
        r, _ = stats.pearsonr(x, y)
        slope, intercept, *_ = stats.linregress(x, y)
        diff = y - x
        loa = 1.96 * diff.std(ddof=1)

        lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
        pad = 0.04 * (hi - lo)
        lo, hi = lo - pad, hi + pad

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.plot([lo, hi], [lo, hi], ls="--", color=COLOURS["identity"],
                linewidth=1.2, zorder=1, label="identity (perfect stability)")
        for yr in sorted(pairs["to_year"].unique()):
            mask = pairs["to_year"].to_numpy() == yr
            ax.scatter(x[mask], y[mask], s=24, alpha=0.6, zorder=3,
                       color=COLOURS["year"].get(int(yr), COLOURS["below_bp"]),
                       edgecolors="white", linewidths=0.4, label=str(int(yr)))
        xs = np.array([lo, hi])
        ax.plot(xs, intercept + slope * xs, color="#333", linewidth=2,
                zorder=4, label=f"fit (slope {slope:.2f})")

        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_xlabel(f"{label} — last summer{unit}")
        ax.set_ylabel(f"{label} — this summer{unit}")
        ax.set_title(f"Year-to-year predictability of {label}",
                     fontsize=12, fontweight="bold")
        ax.text(0.03, 0.97,
                f"r = {r:+.2f}  (r² = {r * r:.3f})\n"
                f"slope = {slope:.2f}\n"
                f"95% limits of agreement: ±{loa:.1f}{unit}",
                transform=ax.transAxes, va="top", ha="left", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          alpha=0.85, edgecolor="#ccc"))
        ax.legend(fontsize=8, loc="lower right", framealpha=0.9)
        fig.tight_layout()
        save_figure(fig, fname, out_dir)
        pairs.to_csv(resolve_output(out_dir, f"{fname}.csv"), index=False)


# ─────────────────────────────────────────────────────────────
#  « breakpoint ICC forest plot »
# ─────────────────────────────────────────────────────────────

def plot_breakpoint_icc(out_dir: Path,
                         csv_stem: str = "breakpoint_icc",
                         figure_stem: str = "breakpoint_icc_forest",
                         title_suffix: str = "",
                         ) -> None:
    """Forest plot of one-way ICC(1,1) for multi-year breakpoints.

    Reads the cohort-specific ICC CSV (default ``breakpoint_icc.csv``)
    and renders a forest with point estimates, 95% CIs, the n / k
    label next to each row, and three reference lines: ICC = 0 (no
    repeatability — within = between), ICC = 0.5 ("moderate"
    repeatability, Koo & Li 2016), and ICC = 0.75 ("good").  Pass
    ``csv_stem='breakpoint_icc_strict'`` and
    ``figure_stem='breakpoint_icc_strict_forest'`` to render a
    cohort-suffixed variant.
    """
    import matplotlib.pyplot as plt

    icc_path = resolve_input(out_dir, f"{csv_stem}.csv")
    if not icc_path.exists():
        log.info("  %s.csv not found — skipping forest plot", csv_stem)
        return
    df = pd.read_csv(icc_path)
    if df.empty:
        return

    setup_figure()

    df = df.copy()
    df["row_label"] = df.apply(
        lambda r: f"{r['predictor']} — {r['mode']}", axis=1)
    # Plot rows top-to-bottom in the order they appear in the CSV;
    # invert the y-axis so the first row sits at the top.
    df = df.reset_index(drop=True)
    y = np.arange(len(df))[::-1]

    colours = {
        "raw":                                  "#0072B2",   # Wong blue
        "raw (SE-cohort)":                      "#56B4E9",   # Wong sky blue
        "measurement-corrected":                "#009E73",   # Wong green
        "residual":                             "#E69F00",   # Wong orange
        "parity + DIM residual":                "#E69F00",   # Wong orange
        "parity + DIM + meas-corr":             "#CC79A7",   # Wong pink
        "parity + DIM + yield residual":        "#D55E00",   # Wong vermillion
        "parity + DIM + yield + meas-corr":     "#000000",   # black (fully corrected)
    }
    point_colours = [colours.get(m, "#999999") for m in df["mode"]]

    fig, ax = plt.subplots(figsize=(8.0, 0.55 * len(df) + 2.2))

    # Reference bands / lines
    ax.axvspan(0.50, 0.75, color="#dddddd", alpha=0.45,
               label="moderate (0.50–0.75)")
    ax.axvspan(0.75, 1.00, color="#bbbbbb", alpha=0.55,
               label="good (≥ 0.75)")
    ax.axvline(0.0,  color="#444444", lw=1.2, ls="--", zorder=1)
    ax.axvline(0.50, color="#666666", lw=0.8, ls=":",  zorder=1)
    ax.axvline(0.75, color="#666666", lw=0.8, ls=":",  zorder=1)

    # CI bars + points
    for i, (_, r) in zip(y, df.iterrows()):
        ax.hlines(i, r["ci_lower"], r["ci_upper"],
                  color=colours.get(r["mode"], "#999999"),
                  linewidth=2.4, zorder=3)
        ax.plot([r["ci_lower"], r["ci_upper"]],
                [i, i], marker="|", markersize=10,
                color=colours.get(r["mode"], "#999999"),
                linestyle="none", zorder=4)
    ax.scatter(df["icc"], y, s=80, c=point_colours,
               edgecolors="white", linewidths=1.2, zorder=5)

    # Right-side annotation: n animals / n obs / p
    for i, (_, r) in zip(y, df.iterrows()):
        sig = ""
        if pd.notna(r["p"]):
            sig = (" ***" if r["p"] < .001 else
                   " **" if r["p"] < .01 else
                   " *" if r["p"] < .05 else "")
        ax.text(1.05, i,
                f"  n={int(r['n_animals'])}  obs={int(r['n_obs'])}  "
                f"p={r['p']:.2f}{sig}",
                va="center", ha="left", fontsize=9, color="#333333",
                transform=ax.get_yaxis_transform())

    ax.set_yticks(y)
    ax.set_yticklabels(df["row_label"], fontsize=10)
    ax.set_xlabel("ICC(1,1)  —  point estimate ± 95% CI")
    ax.set_xlim(-1.05, 1.05)
    title_main = "Repeatability of individual cow breakpoints across summers"
    if title_suffix:
        title_main = f"{title_main}  —  {title_suffix}"
    ax.set_title(
        f"{title_main}\n"
        "(Shrout & Fleiss 1979 / McGraw & Wong 1996; multi-year converged subset)",
        fontsize=11, loc="left", pad=10,
    )
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # Mode legend (raw vs residual) only — reference bands are
    # already labelled by the axvspan calls below the curves.
    from matplotlib.lines import Line2D
    # Build legend dynamically from the modes that actually appear,
    # keyed in a sensible visual order.
    visual_order = [
        "raw", "raw (SE-cohort)", "measurement-corrected",
        "parity + DIM residual", "parity + DIM + meas-corr",
        "parity + DIM + yield residual", "parity + DIM + yield + meas-corr",
        "residual",  # legacy alias
    ]
    seen = list(df["mode"].unique())
    handles = [
        Line2D([0], [0], marker="o", color=colours.get(m, "#999999"),
               lw=0, markersize=8, label=m)
        for m in visual_order if m in seen
    ]
    handles.append(Line2D([0], [0], color="#444444", lw=1.2, ls="--",
                          label="ICC = 0"))
    ax.legend(handles=handles, loc="upper left", frameon=False,
              fontsize=8, ncol=3, bbox_to_anchor=(0.0, -0.10))
    fig.tight_layout()
    save_figure(fig, figure_stem, out_dir)


# ─────────────────────────────────────────────────────────────
#  « longitudinal breakpoint stability Sankey »
# ─────────────────────────────────────────────────────────────

def plot_longitudinal_sankey(bs: pd.DataFrame, out_dir: Path) -> None:
    """Alluvial plot tracking breakpoint adaptation across years.

    Closed cohort: only animals with a converged breakpoint in every
    study year are included.  Column 0 is the baseline year (all
    animals start as "stable").  Subsequent columns show the year-
    to-year Δ category (strongly decreased / decreased / stable /
    increased / strongly increased).  Bezier bands track individual
    animals between columns.

    Duplicate entries per animal-year (multiple date_enter) are
    resolved by averaging breakpoints before computing deltas.
    """
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as MplPath
    setup_figure()

    CATS = [
        ("strongly decreased", -np.inf, -3),
        ("decreased",          -3,      -1),
        ("stable",             -1,       1),
        ("increased",           1,       3),
        ("strongly increased",  3,       np.inf),
    ]
    CAT_COLOURS = {
        "strongly decreased": "#0072B2",
        "decreased":          "#56B4E9",
        "stable":             "#999999",
        "increased":          "#E69F00",
        "strongly increased": "#D55E00",
    }
    cat_labels = [c[0] for c in CATS]

    def _classify(delta: float) -> str:
        for label, lo, hi in CATS:
            if hi == np.inf and delta > lo:
                return label
            if lo == -np.inf and delta <= hi:
                return label
            if lo < delta <= hi:
                return label
        return "stable"

    def _bezier_band(ax, x0, y0_bot, y0_top, x1, y1_bot, y1_top,
                     colour, alpha=0.35):
        """Draw a curved band between two vertical extents."""
        xm = (x0 + x1) / 2
        verts = [
            (x0, y0_bot), (xm, y0_bot), (xm, y1_bot), (x1, y1_bot),
            (x1, y1_top), (xm, y1_top), (xm, y0_top), (x0, y0_top),
            (x0, y0_bot),
        ]
        codes = [
            MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
            MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
            MplPath.CLOSEPOLY,
        ]
        ax.add_patch(mpatches.PathPatch(
            MplPath(verts, codes), facecolor=colour, edgecolor="none",
            alpha=alpha))

    for bp_col, conv_col, bp_label, fname in [
        ("thi_breakpoint", "thi_converged", "THI breakpoint",
         "sankey_longitudinal_thi"),
        ("temp_breakpoint", "temp_converged", "Barn temp breakpoint",
         "sankey_longitudinal_temp"),
    ]:
        conv = bs[bs[conv_col] == True].dropna(subset=[bp_col])
        if conv.empty:
            continue

        # ── Closed cohort: all years required ───────────────────
        mean_bp = (conv.groupby(["animal_id", "year"])[bp_col]
                       .mean().reset_index())
        years = sorted(mean_bp["year"].unique().astype(int))
        if len(years) < 2:
            continue

        year_counts = mean_bp.groupby("animal_id")["year"].nunique()
        cohort_ids = year_counts[year_counts == len(years)].index
        if len(cohort_ids) < 10:
            log.info("  %s: too few animals in all %d years (%d < 10), "
                     "skipping", fname, len(years), len(cohort_ids))
            continue

        cohort = (mean_bp[mean_bp["animal_id"].isin(cohort_ids)]
                  .pivot(index="animal_id", columns="year",
                         values=bp_col))
        n_animals = len(cohort)
        log.info("  %s: %d animals present in all %d years",
                 fname, n_animals, len(years))

        # Column 0 = baseline year (all "stable"), then year-to-year Δ
        col_labels = [str(years[0])] + [
            f"{years[i]}\u2192{years[i+1]}" for i in range(len(years) - 1)]
        n_cols = len(col_labels)

        animal_cats = {}
        for aid in cohort.index:
            cats = ["stable"]
            for i in range(len(years) - 1):
                cats.append(_classify(
                    cohort.loc[aid, years[i + 1]]
                    - cohort.loc[aid, years[i]]))
            animal_cats[aid] = cats

        # ── Count per column per category ───────────────────────
        col_counts = []
        for ci in range(n_cols):
            counts = {cat: 0 for cat in cat_labels}
            for cats in animal_cats.values():
                counts[cats[ci]] += 1
            col_counts.append(counts)

        # ── Flow counts between consecutive columns ─────────────
        flows = []
        for ci in range(n_cols - 1):
            flow = {}
            for cats in animal_cats.values():
                key = (cats[ci], cats[ci + 1])
                flow[key] = flow.get(key, 0) + 1
            flows.append(flow)

        # ── Draw ────────────────────────────────────────────────
        fig_width = 5 + 3.5 * n_cols
        fig, ax = plt.subplots(figsize=(fig_width, 7))

        col_width = 0.3
        gap = 0.02
        x_positions = list(range(n_cols))
        total_height = 1.0

        col_y_positions = []
        for ci in range(n_cols):
            total_n = sum(col_counts[ci].values())
            if total_n == 0:
                col_y_positions.append({})
                continue
            usable_height = total_height - gap * (len(cat_labels) - 1)
            positions = {}
            y_cursor = 0
            for cat in cat_labels:
                h = (usable_height * col_counts[ci][cat] / total_n
                     if total_n > 0 else 0)
                positions[cat] = (y_cursor, y_cursor + h)
                y_cursor += h + gap
            col_y_positions.append(positions)

        # Rectangles + labels
        for ci in range(n_cols):
            x = x_positions[ci]
            for cat in cat_labels:
                count = col_counts[ci][cat]
                if count == 0:
                    continue
                y_bot, y_top = col_y_positions[ci][cat]
                rect = mpatches.FancyBboxPatch(
                    (x - col_width / 2, y_bot), col_width,
                    y_top - y_bot, boxstyle="round,pad=0.01",
                    facecolor=CAT_COLOURS[cat], edgecolor="#333",
                    linewidth=0.8, alpha=0.85)
                ax.add_patch(rect)

                mid_y = (y_bot + y_top) / 2
                if y_top - y_bot > 0.04:
                    ax.text(x, mid_y, f"{count}",
                            ha="center", va="center", fontsize=10,
                            fontweight="bold", color="white")

            ax.text(x, -0.08, col_labels[ci],
                    ha="center", va="top", fontsize=11,
                    fontweight="bold")
            if ci == 0:
                ax.text(x, -0.13, "(baseline)",
                        ha="center", va="top", fontsize=9,
                        color="#666")

        # Bezier flow bands
        for ci in range(n_cols - 1):
            flow = flows[ci]
            if not flow:
                continue

            src_cursor = {cat: col_y_positions[ci][cat][0]
                         for cat in cat_labels
                         if cat in col_y_positions[ci]}
            tgt_cursor = {cat: col_y_positions[ci + 1][cat][0]
                         for cat in cat_labels
                         if cat in col_y_positions[ci + 1]}

            x0 = x_positions[ci] + col_width / 2
            x1 = x_positions[ci + 1] - col_width / 2

            for cat_from in cat_labels:
                for cat_to in cat_labels:
                    count = flow.get((cat_from, cat_to), 0)
                    if count == 0:
                        continue

                    src_total = col_y_positions[ci].get(
                        cat_from, (0, 0))
                    src_h = ((src_total[1] - src_total[0]) * count
                             / max(col_counts[ci][cat_from], 1))
                    y0_bot = src_cursor[cat_from]
                    y0_top = y0_bot + src_h
                    src_cursor[cat_from] = y0_top

                    tgt_total = col_y_positions[ci + 1].get(
                        cat_to, (0, 0))
                    tgt_h = ((tgt_total[1] - tgt_total[0]) * count
                             / max(col_counts[ci + 1][cat_to], 1))
                    y1_bot = tgt_cursor[cat_to]
                    y1_top = y1_bot + tgt_h
                    tgt_cursor[cat_to] = y1_top

                    _bezier_band(ax, x0, y0_bot, y0_top, x1, y1_bot,
                                 y1_top, CAT_COLOURS[cat_from])

                    if count >= 2:
                        mid_x = (x0 + x1) / 2
                        mid_y = ((y0_bot + y0_top + y1_bot + y1_top)
                                 / 4)
                        ax.text(mid_x, mid_y, str(count),
                                ha="center", va="center", fontsize=7,
                                color="#333", alpha=0.7)

        # Legend
        legend_handles = [mpatches.Patch(color=CAT_COLOURS[cat],
                                         label=cat)
                          for cat in cat_labels]
        ax.legend(handles=legend_handles, loc="upper right",
                  fontsize=8, title="\u0394 breakpoint category",
                  title_fontsize=9)

        ax.set_xlim(-0.6, n_cols - 0.4)
        ax.set_ylim(-0.2,
                    total_height + gap * len(cat_labels) + 0.05)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines[:].set_visible(False)

        ax.set_title(
            f"Breakpoint adaptation: {bp_label}\n"
            f"({n_animals} animals in all {len(years)} years, "
            f"\u0394 thresholds: \u00b11 stable, \u00b13 moderate, >3 strong)",
            fontsize=13, fontweight="bold")
        fig.tight_layout()
        save_figure(fig, fname, out_dir)


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
        fig.write_html(str(resolve_output(out_dir, f"sankey_{prefix}.html")),
                        include_plotlyjs="cdn")
        log.info("  Saved sankey_%s.html", prefix)


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────


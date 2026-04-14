#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — viz_runner                                           ║
# ║  « CLI entry point for figure generation »                     ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Reads CSVs from the stats stage and produces all Frontiers     ║
# ║  manuscript figures.  Imports plot functions from the viz_*     ║
# ║  family modules.                                                ║
# ╚══════════════════════════════════════════════════════════════════╝
"""CLI orchestration for figure generation.

Usage::

    digimuh-plots --data results/broken_stick
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger("digimuh.viz")

from digimuh.viz_breakpoints import (
    plot_grouped_boxplots,
    plot_paired_below_above,
    plot_paired_rumen_vs_resp,
    plot_spearman,
    plot_climate,
    plot_predictors,
    plot_stability,
    plot_thi_vs_temp_scatter,
    plot_bodytemp_vs_resp_scatter,
    plot_examples,
)
from digimuh.viz_temporal import (
    plot_circadian_null_model,
    plot_thi_daily_profile,
    plot_crossing_raster,
)
from digimuh.viz_correlation import (
    plot_cross_correlation,
    plot_derivative_ccf,
    plot_event_triggered_average,
    plot_climate_eta,
)
from digimuh.viz_production import plot_tnf_yield
from digimuh.viz_longitudinal import (
    plot_longitudinal_breakpoints,
    plot_breakpoint_raincloud,
    plot_longitudinal_sankey,
    plot_threshold_sankey,
)

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
        ("Event-triggered average (8-11h)", lambda: plot_event_triggered_average(
            d, traces_file="event_triggered_traces_filtered.csv",
            suffix="_8to11h", title_extra=" (crossings 8:00–11:00 only)")),
        ("Climate ETA", lambda: plot_climate_eta(d)),
        ("TNF vs yield", lambda: plot_tnf_yield(d)),
        ("Longitudinal breakpoints", lambda: plot_longitudinal_breakpoints(bs, d)),
        ("Breakpoint raincloud", lambda: plot_breakpoint_raincloud(d)),
        ("Longitudinal Sankey", lambda: plot_longitudinal_sankey(bs, d)),
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


if __name__ == "__main__":
    main()

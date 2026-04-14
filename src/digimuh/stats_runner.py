#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_runner                                         ║
# ║  « CLI entry point for the statistical analysis pipeline »     ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Reads CSVs from extract stage, runs all statistical analyses, ║
# ║  and writes result CSVs.  Imports functions from stats_core,    ║
# ║  stats_temporal, stats_production, and stats_longitudinal.      ║
# ╚══════════════════════════════════════════════════════════════════╝
"""CLI orchestration for the broken-stick statistical pipeline.

Usage::

    digimuh-stats --data results/broken_stick
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

log = logging.getLogger("digimuh.stats")


# ─────────────────────────────────────────────────────────────────
#  Imports from library modules (lazy — keeps startup fast)
# ─────────────────────────────────────────────────────────────────

from digimuh.stats_core import (
    benjamini_hochberg,
    p_to_stars,
    run_broken_stick_fits,
    compute_spearman,
    compute_below_above,
    run_statistical_tests,
)
from digimuh.stats_temporal import (
    compute_cross_correlation,
    compute_circadian_null_model,
    compute_thi_daily_profile,
    compute_derivative_ccf,
    compute_event_triggered_average,
    compute_crossing_times,
    compute_climate_eta,
)
from digimuh.stats_production import (
    compute_thermoneutral_fraction,
    compute_tnf_yield_analysis,
)
from digimuh.stats_longitudinal import (
    compute_stability,
    make_summary_table,
    _run_longitudinal_tests,
)

def main() -> None:
    from digimuh.console import (
        setup_logging, banner, section, result_table, kv, stars_styled,
        progress, done, reset_steps,
    )

    parser = argparse.ArgumentParser(description="Statistical analysis for broken-stick")
    parser.add_argument("--data", type=Path, required=True,
                        help="Directory with CSVs from extraction step")
    parser.add_argument("--no-resp", action="store_true",
                        help="Skip respiration analysis (Frontiers paper: "
                             "rumen temperature only)")
    parser.add_argument("--frontiers", action="store_true",
                        help="Frontiers paper mode: skip Davies/pscore/Hill "
                             "(reserved for COMPAG companion papers)")
    args = parser.parse_args()

    setup_logging()
    reset_steps()
    d = args.data

    if args.frontiers:
        banner("Broken-stick analysis (Frontiers mode)")
        log.info("Davies/pscore/Hill SKIPPED (--frontiers, reserved for COMPAG)")
    else:
        banner("Broken-stick statistical analysis")

    log.info("Loading CSVs from %s", d)
    rumen = pd.read_csv(d / "rumen_barn.csv")
    resp_path = d / "respiration_barn.csv"
    if args.no_resp:
        resp = pd.DataFrame()
        log.info("Respiration analysis SKIPPED (--no-resp)")
    else:
        resp = pd.read_csv(resp_path) if resp_path.exists() else pd.DataFrame()
    prod = pd.read_csv(d / "production.csv") if (d / "production.csv").exists() else pd.DataFrame()

    # ── 1. Model fitting ─────────────────────────────────────
    if args.frontiers:
        section("Model fitting", "Broken-stick regression only")
    else:
        section("Model fitting", "Broken-stick, Davies/pscore, Hill (4PL)")
    bs = run_broken_stick_fits(rumen, resp, frontiers_only=args.frontiers)
    if not prod.empty:
        bs = bs.merge(
            prod[["animal_id", "year", "mean_milk_yield_kg", "lactation_nr"]],
            on=["animal_id", "year"], how="left",
        )
    bs.to_csv(d / "broken_stick_results.csv", index=False)

    prefixes = [("thi", "THI > body temp"), ("temp", "Barn temp > body temp")]
    if not args.no_resp:
        prefixes += [("resp_thi", "THI > resp"), ("resp_temp", "Barn temp > resp")]

    if args.frontiers:
        # Simplified table: broken-stick only
        conv_rows = []
        for prefix, label in prefixes:
            n_conv = (bs[f"{prefix}_converged"] == True).sum()
            conv_rows.append([label, len(bs), f"{n_conv} ({100*n_conv/len(bs):.0f}%)"])
        result_table("Broken-stick convergence",
                     ["Model", "n", "Converged"], conv_rows)
    else:
        conv_rows = []
        for prefix, label in prefixes:
            n_conv = (bs[f"{prefix}_converged"] == True).sum()
            davies_col = f"{prefix}_davies_p"
            pscore_col = f"{prefix}_pscore_p"
            hill_col = f"{prefix}_hill_converged"
            hill_bend = f"{prefix}_hill_bend"

            n_davies = (bs[davies_col] < 0.05).sum() if davies_col in bs.columns else 0
            n_pscore = (bs[pscore_col] < 0.05).sum() if pscore_col in bs.columns else 0
            n_hill = (bs[hill_col] == True).sum() if hill_col in bs.columns else 0
            n_bend = bs.loc[bs.get(hill_col, pd.Series(dtype=bool)) == True, hill_bend].notna().sum() if hill_col in bs.columns else 0

            conv_rows.append([
                label, len(bs),
                f"{n_davies} ({100*n_davies/len(bs):.0f}%)",
                f"{n_pscore} ({100*n_pscore/len(bs):.0f}%)",
                f"{n_conv} ({100*n_conv/len(bs):.0f}%)",
                f"{n_hill} ({100*n_hill/len(bs):.0f}%)",
                n_bend,
            ])
        result_table(
            "Model convergence",
            ["Model", "n", "Davies p<.05", "Pscore p<.05",
             "BS converged", "Hill converged", "Bend valid"],
            conv_rows,
        )
        for prefix, label in prefixes:
            hill_col = f"{prefix}_hill_converged"
            hill_bend_col = f"{prefix}_hill_bend"
            if hill_col not in bs.columns:
                continue
            bends = bs.loc[bs[hill_col] == True, hill_bend_col].dropna()
            if len(bends) > 0:
                kv(f"{label} Hill bend median [IQR]",
                   f"{bends.median():.1f} [{bends.quantile(0.25):.1f} - {bends.quantile(0.75):.1f}]")

    # ── 1b. Test breakpoints against literature THI threshold ─
    section("Breakpoint vs literature threshold",
            "Fisher resampling: individual breakpoints vs THI 68.8")
    thi_ref = 68.8  # Hoffmann et al. (2020) mild stress threshold
    thi_conv = bs[bs["thi_converged"] == True]["thi_breakpoint"].dropna()
    if len(thi_conv) >= 10:
        from rerandomstats import FisherResamplingTest
        # Test: are the individual breakpoints centred on 68.8?
        # Create a "reference" sample of the same size, all = 68.8
        ref_sample = [thi_ref] * len(thi_conv)
        p_vs_ref = FisherResamplingTest(
            data_a=thi_conv.tolist(),
            data_b=ref_sample,
            func="medianDiff",
            combination_n=20_000,
        ).main()
        median_bp = thi_conv.median()
        diff = median_bp - thi_ref
        result_table(
            f"Individual THI breakpoints vs reference THI {thi_ref}",
            ["n", "Median BP", "Reference", "Diff", "p", "Sig."],
            [[len(thi_conv), median_bp, thi_ref, diff, p_vs_ref,
              stars_styled(p_to_stars(p_vs_ref))]],
        )
        # Per-year breakdown
        year_rows = []
        for year in sorted(bs["year"].unique().astype(int)):
            yr_bps = bs[(bs["year"] == year) & (bs["thi_converged"] == True)]["thi_breakpoint"].dropna()
            if len(yr_bps) >= 5:
                p_yr = FisherResamplingTest(
                    data_a=yr_bps.tolist(),
                    data_b=[thi_ref] * len(yr_bps),
                    func="medianDiff",
                    combination_n=20_000,
                ).main()
                year_rows.append([year, len(yr_bps), yr_bps.median(),
                                  yr_bps.median() - thi_ref, p_yr,
                                  stars_styled(p_to_stars(p_yr))])
        if year_rows:
            # BH-FDR across years
            raw_ps = np.array([r[4] for r in year_rows])
            adj_ps = benjamini_hochberg(raw_ps)
            for r, adj_p in zip(year_rows, adj_ps):
                r[4] = adj_p
                r[5] = stars_styled(p_to_stars(adj_p))
            result_table(
                "Per-year (BH-FDR corrected)",
                ["Year", "n", "Median", "Diff", "p adj", "Sig."],
                year_rows,
            )

    # ── 2. Spearman ──────────────────────────────────────────
    section("Spearman correlations")
    spearman = compute_spearman(rumen, resp)
    spearman.to_csv(d / "spearman_correlations.csv", index=False)
    kv("Body temp vs THI median rs", f"{spearman['thi_rs'].median():.3f}")
    kv("Body temp vs barn temp median rs", f"{spearman['temp_rs'].median():.3f}")

    # ── 3. Below/above ───────────────────────────────────────
    section("Below/above breakpoint",
            "Per-animal means split at individual THI breakpoint")
    beh = compute_below_above(rumen, resp, bs)
    beh.to_csv(d / "behavioural_response.csv", index=False)
    kv("Animals with paired data", len(beh))

    # ── 4. Fisher resampling tests ───────────────────────────
    section("Statistical tests",
            "Fisher resampling (reRandomStats), BH-FDR corrected")
    tests = run_statistical_tests(beh)
    tests.to_csv(d / "statistical_tests.csv", index=False)

    if not tests.empty:
        test_rows = []
        for _, t in tests.iterrows():
            test_rows.append([
                int(t["year"]),
                t["test"].replace(" (Fisher medianDiff)", ""),
                int(t["n"]),
                t["p_raw"],
                t["p_adj"],
                stars_styled(t["stars"]),
                t["median_diff"],
            ])
        result_table(
            "Within-year tests (BH-FDR corrected)",
            ["Year", "Test", "n", "p raw", "p adj", "Sig.", "Median diff"],
            test_rows,
            highlight_col=5,
        )

    # ── 5. Rumen temperature circadian null model ──────────
    section("Rumen circadian null model",
            "Hourly rumen temp profile on cool days (no heat stress) vs stress days")
    circadian = compute_circadian_null_model(rumen, bs)
    circadian.to_csv(d / "circadian_null_model.csv", index=False)

    if not circadian.empty:
        for day_type in ["cool", "stress"]:
            sub = circadian[circadian["day_type"] == day_type]
            if not sub.empty:
                # Grand mean across animals per hour
                hourly = sub.groupby("hour")["body_temp_mean"].agg(["mean", "std", "count"])
                peak_h = hourly["mean"].idxmax()
                trough_h = hourly["mean"].idxmin()
                amplitude = hourly["mean"].max() - hourly["mean"].min()
                kv(f"{day_type.capitalize()} days: n animals", sub["animal_id"].nunique())
                kv(f"  Amplitude (peak-trough)", f"{amplitude:.3f} °C")
                kv(f"  Peak hour", f"{peak_h}:00")
                kv(f"  Trough hour", f"{trough_h}:00")

    # ── 6. THI daily exceedance profile ──────────────────────
    section("THI daily exceedance",
            "Barn THI by clock hour and month — when does heat stress occur?")
    thi_profile = compute_thi_daily_profile(rumen, bs)
    thi_profile.to_csv(d / "thi_daily_profile.csv", index=False)

    if not thi_profile.empty:
        herd_bp = thi_profile["herd_median_bp"].iloc[0]
        kv("Herd median THI breakpoint", f"{herd_bp:.1f}")
        # For each month, find the hour range where mean THI exceeds breakpoint
        for ml in sorted(thi_profile["month_label"].unique()):
            msub = thi_profile[thi_profile["month_label"] == ml]
            exceeds = msub[msub["thi_mean"] > herd_bp]
            if not exceeds.empty:
                kv(f"  {ml}", f"THI > bp from {exceeds['hour'].min()}:00 "
                   f"to {exceeds['hour'].max()}:00 "
                   f"({len(exceeds)} hours)")
            else:
                kv(f"  {ml}", "THI stays below breakpoint all day")

    # ── 7. Breakpoint crossing raster ────────────────────────
    section("Breakpoint crossing times",
            "When during the 24h cycle does each cow's THI cross her breakpoint?")
    crossing_times = compute_crossing_times(rumen, bs)
    crossing_times.to_csv(d / "crossing_times.csv", index=False)

    if not crossing_times.empty:
        for pred in ["thi", "temp"]:
            psub = crossing_times[crossing_times["predictor"] == pred]
            if not psub.empty:
                kv(f"{pred.upper()} crossings",
                   f"{len(psub)} events from {psub['animal_id'].nunique()} animals")
                kv(f"  Median clock time",
                   f"{psub['day_fraction'].median():.1f}h "
                   f"[IQR {psub['day_fraction'].quantile(0.25):.1f}–"
                   f"{psub['day_fraction'].quantile(0.75):.1f}h]")

    # ── 8. Cross-correlation ─────────────────────────────────
    section("Cross-correlation",
            "Climate vs rumen temp, below/above breakpoint\n"
            "Raw = includes 24h diurnal cycle; Detrended = per-day mean removed")
    xcorr = compute_cross_correlation(rumen, bs)
    xcorr.to_csv(d / "cross_correlation.csv", index=False)

    for variant, variant_label in [("raw", "Raw"), ("detrended", "Detrended (diurnal removed)")]:
        xcorr_rows = []
        vsub = xcorr[xcorr["variant"] == variant] if "variant" in xcorr.columns else xcorr
        for pred in ["thi", "temp"]:
            for region in ["below", "above"]:
                sub = vsub[(vsub["predictor"] == pred) & (vsub["region"] == region)
                           & (vsub["lag"] == 0)]
                if not sub.empty:
                    xcorr_rows.append([
                        pred.upper(), region,
                        sub["animal_id"].nunique(),
                        sub["xcorr"].median(),
                        sub["xcorr"].quantile(0.25),
                        sub["xcorr"].quantile(0.75),
                    ])
        if xcorr_rows:
            result_table(
                f"Cross-correlation at lag 0 — {variant_label}",
                ["Predictor", "Region", "n animals",
                 "Median r", "Q25", "Q75"],
                xcorr_rows,
            )

    # ── 6. Derivative cross-correlation ──────────────────────
    section("Derivative cross-correlation",
            "d(climate)/dt vs d(body_temp)/dt — temporal coupling of changes")
    dccf = compute_derivative_ccf(rumen, bs)
    dccf.to_csv(d / "derivative_ccf.csv", index=False)

    dccf_rows = []
    for pred in ["thi", "temp"]:
        for region in ["below", "above"]:
            sub = dccf[(dccf["predictor"] == pred) & (dccf["region"] == region)]
            if sub.empty:
                continue
            # Find peak lag across animals
            agg = sub.groupby("lag_minutes")["dxcorr"].mean()
            peak_lag = agg.idxmax()
            peak_r = agg.max()
            lag0_r = agg.get(0, np.nan)
            dccf_rows.append([
                pred.upper(), region,
                sub["animal_id"].nunique(),
                lag0_r, peak_r, int(peak_lag),
            ])
    if dccf_rows:
        result_table(
            "Derivative CCF summary",
            ["Predictor", "Region", "n animals",
             "r at lag=0", "Peak r", "Peak lag (min)"],
            dccf_rows,
        )

    # ── 7. Event-triggered average ───────────────────────────
    section("Event-triggered average",
            "Rumen temp aligned to THI breakpoint upward crossings")

    # Full ETA (all crossing times)
    eta_traces, eta_summary = compute_event_triggered_average(rumen, bs)
    eta_traces.to_csv(d / "event_triggered_traces.csv", index=False)
    eta_summary.to_csv(d / "event_triggered_summary.csv", index=False)

    if not eta_summary.empty:
        for pred in ["thi", "temp"]:
            psub = eta_summary[eta_summary["predictor"] == pred]
            if not psub.empty:
                kv(f"{pred.upper()} crossing events (all hours)",
                   f"{psub['n_events'].sum()} events from "
                   f"{psub['animal_id'].nunique()} animals "
                   f"(median {psub['n_events'].median():.0f}/animal)")

    # Hour-filtered ETA (morning crossings only, configurable)
    eta_hour_start = 8
    eta_hour_end = 11
    kv("", "")
    kv(f"Hour-filtered ETA", f"{eta_hour_start}:00–{eta_hour_end}:00 only")

    eta_filt_traces, eta_filt_summary = compute_event_triggered_average(
        rumen, bs, crossing_hour_range=(eta_hour_start, eta_hour_end))
    eta_filt_traces.to_csv(d / "event_triggered_traces_filtered.csv", index=False)
    eta_filt_summary.to_csv(d / "event_triggered_summary_filtered.csv", index=False)

    if not eta_filt_summary.empty:
        for pred in ["thi", "temp"]:
            psub = eta_filt_summary[eta_filt_summary["predictor"] == pred]
            if not psub.empty:
                kv(f"{pred.upper()} crossings ({eta_hour_start}–{eta_hour_end}h)",
                   f"{psub['n_events'].sum()} events from "
                   f"{psub['animal_id'].nunique()} animals "
                   f"(median {psub['n_events'].median():.0f}/animal)")

    # ── Climate ETA: THI + barn temp around crossings ────────
    section("Climate ETA",
            "THI and barn temp around breakpoint crossings, normalised to breakpoint")
    climate_eta = compute_climate_eta(
        rumen, bs,
        crossing_hour_range=(eta_hour_start, eta_hour_end))
    climate_eta.to_csv(d / "climate_eta.csv", index=False)

    if not climate_eta.empty:
        for trigger in ["thi", "temp"]:
            tsub = climate_eta[climate_eta["trigger"] == trigger]
            if not tsub.empty:
                n_ev = tsub.groupby(["animal_id", "year", "event_id"]).ngroups
                kv(f"{trigger.upper()} trigger ({eta_hour_start}–{eta_hour_end}h)",
                   f"{n_ev} events from {tsub['animal_id'].nunique()} animals")

    # ── 8. Thermoneutral fraction vs daily milk yield ────────
    section("Thermoneutral fraction (TNF)",
            "Daily fraction below breakpoint vs daily milk yield (cow-specific P95)")
    tnf = compute_thermoneutral_fraction(rumen, bs)
    tnf.to_csv(d / "thermoneutral_fraction.csv", index=False)
    kv("Daily TNF records", len(tnf))
    kv("Animals with TNF", tnf["animal_id"].nunique() if not tnf.empty else 0)

    # Load daily milk yield
    daily_yield_path = d / "daily_milk_yield.csv"
    if daily_yield_path.exists() and not tnf.empty:
        daily_yield = pd.read_csv(daily_yield_path)
        kv("Daily yield records", len(daily_yield))

        tnf_yield = compute_tnf_yield_analysis(tnf, daily_yield)
        tnf_yield.to_csv(d / "tnf_yield.csv", index=False)
        kv("Matched cow-day pairs", len(tnf_yield))
        kv("Animals with pairs", tnf_yield["animal_id"].nunique() if not tnf_yield.empty else 0)

        if not tnf_yield.empty:
            valid = tnf_yield.dropna(subset=["thi_tnf", "relative_yield"])
            if len(valid) >= 20:
                # Spearman on all daily pairs
                rs, p = spearmanr(valid["thi_tnf"], valid["daily_yield_kg"])
                kv("TNF vs daily yield: Spearman rs", f"{rs:.3f} (p={p:.2e})")

                rs_rel, p_rel = spearmanr(valid["thi_tnf"], valid["relative_yield"])
                kv("TNF vs relative yield: Spearman rs", f"{rs_rel:.3f} (p={p_rel:.2e})")

                # Per-cow Spearman (within-animal correlation)
                per_cow_rs = []
                for aid, cow in valid.groupby("animal_id"):
                    if len(cow) >= 10:
                        r, _ = spearmanr(cow["thi_tnf"], cow["relative_yield"])
                        if np.isfinite(r):
                            per_cow_rs.append(r)
                if per_cow_rs:
                    kv("Per-cow Spearman median (within-animal)",
                       f"{np.median(per_cow_rs):.3f} [IQR {np.percentile(per_cow_rs, 25):.3f} - "
                       f"{np.percentile(per_cow_rs, 75):.3f}], n={len(per_cow_rs)} cows")

                # Quartile table
                valid = valid.copy()
                try:
                    valid["tnf_quartile"] = pd.qcut(
                        valid["thi_tnf"], 4, duplicates="drop",
                    )
                    n_bins = valid["tnf_quartile"].nunique()
                    # Rename bins to readable labels
                    bin_labels = {cat: f"Q{i+1}" for i, cat in
                                  enumerate(sorted(valid["tnf_quartile"].dropna().unique()))}
                    if n_bins >= 2:
                        first = sorted(bin_labels.keys())[0]
                        last = sorted(bin_labels.keys())[-1]
                        bin_labels[first] = f"Q1 (low TNF)"
                        bin_labels[last] = f"Q{n_bins} (high TNF)"
                    valid["tnf_quartile"] = valid["tnf_quartile"].map(bin_labels)
                except ValueError:
                    # Fallback: median split
                    med = valid["thi_tnf"].median()
                    valid["tnf_quartile"] = np.where(
                        valid["thi_tnf"] <= med, "Below median TNF", "Above median TNF")

                q_rows = []
                for q in sorted(valid["tnf_quartile"].dropna().unique()):
                    qdata = valid[valid["tnf_quartile"] == q]
                    if len(qdata) > 0:
                        q_rows.append([
                            q, len(qdata),
                            f"{qdata['thi_tnf'].median():.2f}",
                            f"{qdata['daily_yield_kg'].median():.1f}",
                            f"{qdata['relative_yield'].median():.3f}",
                        ])
                if q_rows:
                    result_table(
                        "Daily yield by thermoneutral fraction quartile",
                        ["Quartile", "n days", "Median TNF",
                         "Median yield (kg)", "Median rel. yield"],
                        q_rows,
                    )

                # ── Heat stress days only ─────────────────────────
                # Restrict to days where the cow actually experienced
                # some time above her breakpoint (TNF < threshold)
                kv("", "")  # spacer
                kv("Heat stress days only", "(TNF < threshold)")

                hs_rows = []
                for threshold in [0.95, 0.90, 0.80, 0.70]:
                    stressed = valid[valid["thi_tnf"] < threshold]
                    n_days = len(stressed)
                    n_cows = stressed["animal_id"].nunique()
                    if n_days < 20:
                        hs_rows.append([
                            f"< {threshold:.2f}", n_days, n_cows,
                            "", "", "", "",
                        ])
                        continue

                    rs_hs, p_hs = spearmanr(stressed["thi_tnf"],
                                            stressed["relative_yield"])
                    # Per-cow within-animal (stressed days only)
                    cow_rs_hs = []
                    for aid, cow in stressed.groupby("animal_id"):
                        if len(cow) >= 5:
                            r, _ = spearmanr(cow["thi_tnf"], cow["relative_yield"])
                            if np.isfinite(r):
                                cow_rs_hs.append(r)

                    hs_rows.append([
                        f"< {threshold:.2f}", n_days, n_cows,
                        f"{rs_hs:.3f}", f"{p_hs:.2e}",
                        f"{np.median(cow_rs_hs):.3f}" if cow_rs_hs else "",
                        len(cow_rs_hs),
                    ])

                result_table(
                    "TNF vs relative yield — heat stress days only",
                    ["TNF thresh", "n days", "n cows",
                     "Pooled rs", "p", "Per-cow rs", "n cows (rs)"],
                    hs_rows,
                )
    else:
        if not daily_yield_path.exists():
            log.info("  daily_milk_yield.csv not found — re-run digimuh-extract")

    # ── 8. Stability ─────────────────────────────────────────
    section("Breakpoint stability", "Repeat animals across years")
    pairs, icc = compute_stability(bs)
    if not pairs.empty:
        pairs.to_csv(d / "breakpoint_stability.csv", index=False)
    kv("ICC", f"{icc:.3f}")
    kv("Pairs", len(pairs))
    kv("Unique animals", pairs["animal_id"].nunique() if not pairs.empty else 0)

    # ── 7. Longitudinal statistical tests ────────────────────
    section("Longitudinal breakpoint tests",
            "Fisher resampling: absolute + relative change across years")
    _run_longitudinal_tests(bs, d)

    # ── 8. Summary ───────────────────────────────────────────
    section("Summary table")
    summary = make_summary_table(bs)
    summary.to_csv(d / "summary_table.csv", index=False)

    done(f"All statistics complete. Output in: {d}")



if __name__ == "__main__":
    main()

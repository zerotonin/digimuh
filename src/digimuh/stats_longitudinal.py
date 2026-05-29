#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_longitudinal                                   ║
# ║  « breakpoint stability, ICC, and year-over-year tests »       ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Repeat-animal stability analysis: ICC, Fisher resampling of    ║
# ║  year pairs, and manuscript summary table.                      ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Longitudinal breakpoint stability analysis."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from rerandomstats import correct_pvalues_array
from scipy.stats import f as f_dist

from digimuh.stats_core import p_to_stars

log = logging.getLogger("digimuh.stats")


# ─────────────────────────────────────────────────────────────
#  « ICC(1,1) — Shrout & Fleiss 1979, McGraw & Wong 1996 »
#
#  Single-rater absolute-agreement ICC for a one-way random-
#  effects design with unbalanced k.  The "rater" here is the
#  year a cow's breakpoint was estimated in; we treat year as
#  exchangeable measurement noise and ask how much of the total
#  variance in breakpoints sits *between* cows rather than
#  *within* a single cow across years.
# ─────────────────────────────────────────────────────────────

def _icc_one_way(values: np.ndarray, groups: np.ndarray,
                 alpha: float = 0.05) -> dict | None:
    """ICC(1,1) on a 1-d array of values + group labels.

    Implements Shrout & Fleiss (1979) Type 1 ICC with the
    unbalanced-design correction k0 from McGraw & Wong (1996).
    The 95% CI is derived from the F-distribution of MS_b / MS_w.

    Returns ``None`` when the data are too sparse (K<3 groups or
    df_w <= 0) or when MS_within is exactly zero (degenerate).
    """
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups)
    mask = np.isfinite(values)
    values, groups = values[mask], groups[mask]
    K = int(pd.Series(groups).nunique())
    N = len(values)
    if K < 3 or N - K < 1:
        return None

    s = pd.Series(values).groupby(pd.Series(groups))
    n_per = s.size().to_numpy(dtype=float)
    means = s.mean().to_numpy()
    grand = values.mean()

    ss_b = float(np.sum(n_per * (means - grand) ** 2))
    ss_w = float(np.sum((s.transform(lambda g: g - g.mean()).to_numpy()) ** 2))
    df_b, df_w = K - 1, N - K
    ms_b = ss_b / df_b
    ms_w = ss_w / df_w
    if ms_w <= 0:
        return None

    # Harmonic-mean k for unbalanced designs (McGraw & Wong 1996).
    k0 = (N - float(np.sum(n_per ** 2)) / N) / (K - 1)
    icc = (ms_b - ms_w) / (ms_b + (k0 - 1.0) * ms_w)

    F0 = ms_b / ms_w
    p = float(f_dist.sf(F0, df_b, df_w))
    f_lo = F0 / float(f_dist.ppf(1.0 - alpha / 2.0, df_b, df_w))
    f_hi = F0 * float(f_dist.ppf(1.0 - alpha / 2.0, df_w, df_b))
    icc_lo = (f_lo - 1.0) / (f_lo + k0 - 1.0)
    icc_hi = (f_hi - 1.0) / (f_hi + k0 - 1.0)

    return {
        "icc": float(icc),
        "ci_lower": float(icc_lo),
        "ci_upper": float(icc_hi),
        "f": float(F0),
        "df1": int(df_b),
        "df2": int(df_w),
        "p": p,
        "n_animals": K,
        "n_obs": N,
        "k_mean": float(k0),
    }


def _residualise_breakpoint(values: np.ndarray, parity: np.ndarray,
                            mean_dim: np.ndarray,
                            mean_yield: np.ndarray | None = None,
                            ) -> np.ndarray | None:
    """Return OLS residuals of breakpoint ~ parity (1/2/3+) + mean_DIM
    (+ optional mean_yield).

    Parity is bucketed into 1, 2, and 3+ to match the Wood-curve
    pooling scheme.  Adding ``mean_yield`` (the cow-year mean
    daily milk yield in kg) controls for the well-known confound
    that high-yielding cows tend to show lower THI breakpoints
    (Bernabucci et al. 2014; Tao et al. 2018), so the resulting
    residuals isolate cow-individual heat tolerance *beyond* what
    parity, lactation stage and production level predict.

    Rows with missing parity, mean_DIM, or mean_yield (when
    requested) are dropped; the returned array aligns with the
    input order, with NaN where a row was dropped.
    """
    import statsmodels.api as sm

    df = pd.DataFrame({
        "y": values,
        "parity_raw": parity,
        "dim": mean_dim,
    })
    if mean_yield is not None:
        df["yield"] = mean_yield
    df["parity_bucket"] = (
        df["parity_raw"].fillna(0).astype(int).clip(upper=3).astype(str)
    )
    drop_cols = ["y", "dim"] + (["yield"] if mean_yield is not None else [])
    fit_df = df.dropna(subset=drop_cols).query("parity_bucket != '0'")
    if len(fit_df) < 6:
        return None
    X = pd.get_dummies(fit_df["parity_bucket"], drop_first=True, dtype=float)
    X["dim"] = fit_df["dim"].to_numpy()
    if mean_yield is not None:
        X["yield"] = fit_df["yield"].to_numpy()
    X = sm.add_constant(X, has_constant="add")
    model = sm.OLS(fit_df["y"].to_numpy(), X.to_numpy()).fit()
    out = np.full(len(df), np.nan)
    out[fit_df.index] = model.resid
    return out


def _mean_dim_per_cow_year(wood_yield: pd.DataFrame) -> pd.DataFrame:
    """Mean DIM per (animal_id, year) from a daily-yield-with-DIM frame.

    ``daily_milk_yield_wood.csv`` already carries ``dim``,
    ``year`` and ``lactation_nr`` per cow-day inside the summer
    analysis window, so we just collapse to the cow-year mean.
    """
    if wood_yield is None or wood_yield.empty:
        return pd.DataFrame(columns=["animal_id", "year", "mean_dim"])
    cols = {"animal_id", "year", "dim"}
    if not cols.issubset(wood_yield.columns):
        return pd.DataFrame(columns=["animal_id", "year", "mean_dim"])
    g = (wood_yield.dropna(subset=["dim"])
                   .groupby(["animal_id", "year"], as_index=False)["dim"]
                   .mean()
                   .rename(columns={"dim": "mean_dim"}))
    g["year"] = g["year"].astype(int)
    return g


def _icc_measurement_corrected(values: np.ndarray, groups: np.ndarray,
                                se_values: np.ndarray,
                                alpha: float = 0.05,
                                n_boot: int = 2000,
                                rng_seed: int = 17,
                                ) -> dict | None:
    """Measurement-error-corrected ICC with parametric-bootstrap CI.

    Decomposes the observed within-cow variance into a fit-noise
    component (mean of squared per-fit SEs) and a residual drift
    component, then recomputes ICC = σ²_between / (σ²_between +
    σ²_drift).  σ²_meas is treated as the iid-residual SE on the
    breakpoint as returned by ``broken_stick_fit`` — a *lower
    bound* on real measurement uncertainty because it ignores
    temporal autocorrelation in the smaXtec readings.

    The 95% CI is obtained by parametric bootstrap on (MS_between,
    MS_within − σ²_meas) treating both as scaled-χ² draws under
    the one-way random-effects model.  ``n_boot=2000`` gives a
    monte-carlo standard error of ~0.01 on the CI quantiles.
    """
    raw = _icc_one_way(values, groups, alpha=alpha)
    if raw is None:
        return None
    # σ²_meas = mean of squared per-fit SEs across the multi-year cohort.
    se_arr = np.asarray(se_values, dtype=float)
    se_arr = se_arr[np.isfinite(se_arr)]
    if se_arr.size == 0:
        return None
    sigma2_meas = float(np.mean(se_arr ** 2))

    K, N = raw["n_animals"], raw["n_obs"]
    df_b, df_w = K - 1, N - K
    k0 = raw["k_mean"]

    # Method-of-moments variance components from the *observed* MS.
    ms_b_obs = raw["f"] * (raw["icc"] >= -1)  # placeholder; recompute from inputs
    # Recover MS values from the raw-ICC F + degrees of freedom.
    # raw['f'] = MS_b / MS_w; absolute scale comes from re-running ANOVA.
    grand = float(np.mean(values))
    s = pd.Series(values).groupby(pd.Series(groups))
    n_per = s.size().to_numpy(dtype=float)
    means = s.mean().to_numpy()
    ss_b = float(np.sum(n_per * (means - grand) ** 2))
    ss_w = float(np.sum(s.transform(lambda g: g - g.mean()).to_numpy() ** 2))
    ms_b_obs = ss_b / df_b
    ms_w_obs = ss_w / df_w

    sigma2_drift = max(ms_w_obs - sigma2_meas, 0.0)
    sigma2_between = max((ms_b_obs - ms_w_obs) / k0, 0.0)
    denom = sigma2_between + sigma2_drift
    icc_c = sigma2_between / denom if denom > 0 else np.nan

    # Parametric bootstrap on (MS_b, MS_w) under one-way random model.
    rng = np.random.default_rng(rng_seed)
    icc_boot = np.empty(n_boot)
    if df_b > 0 and df_w > 0:
        # MS_b ~ (sigma_w² + k0 * sigma_a²) * χ²_{df_b} / df_b under H_random
        # MS_w ~ sigma_w² * χ²_{df_w} / df_w
        sigma_a2_hat = sigma2_between
        sigma_w2_hat = ms_w_obs  # observed within-cow scatter (drift + meas)
        for b in range(n_boot):
            ms_b_b = (sigma_w2_hat + k0 * sigma_a2_hat) * \
                     rng.chisquare(df_b) / df_b
            ms_w_b = sigma_w2_hat * rng.chisquare(df_w) / df_w
            sa2_b = max((ms_b_b - ms_w_b) / k0, 0.0)
            sd2_b = max(ms_w_b - sigma2_meas, 0.0)
            den_b = sa2_b + sd2_b
            icc_boot[b] = sa2_b / den_b if den_b > 0 else 0.0
        ci_lo, ci_hi = np.quantile(icc_boot, [alpha / 2, 1 - alpha / 2])
    else:
        ci_lo = ci_hi = np.nan

    return {
        "icc": float(icc_c),
        "ci_lower": float(ci_lo),
        "ci_upper": float(ci_hi),
        "f": raw["f"],
        "df1": df_b,
        "df2": df_w,
        "p": raw["p"],  # Note: p still tests raw H0 (between > 0)
        "n_animals": K,
        "n_obs": N,
        "k_mean": k0,
        "sigma2_meas": sigma2_meas,
        "ms_within_obs": float(ms_w_obs),
        "sigma2_drift_est": float(sigma2_drift),
        "sigma2_between_est": float(sigma2_between),
    }


def compute_breakpoint_icc(bs_results: pd.DataFrame,
                           wood_yield: pd.DataFrame | None = None,
                           production: pd.DataFrame | None = None,
                           ) -> pd.DataFrame:
    """One-way random-effects ICC(1,1) on multi-year breakpoints.

    For each predictor in {THI, barn temperature} we run two
    flavours of the ICC:

      * ``raw`` — directly on the converged breakpoint values.
      * ``residual`` — on OLS residuals from regressing the
        breakpoint on parity bucket (1/2/3+) plus mean DIM
        across that cow-year's summer window.  Requires
        ``daily_milk_yield_wood.csv``-style frame; silently
        skipped when not available.

    The cohort for each row is the subset of animals with at
    least two converged breakpoints for that predictor.

    Returns one row per (predictor, mode) with point estimate,
    95% CI, F, df, p, n_animals, n_obs, mean k (the harmonic-
    mean k0 used by Shrout & Fleiss for unbalanced designs).
    """
    out_rows: list[dict] = []
    dim_lookup = _mean_dim_per_cow_year(wood_yield)

    # Optional milk-yield lookup (animal_id, year, mean_milk_yield_kg).
    yield_lookup = pd.DataFrame(columns=["animal_id", "year", "mean_milk_yield_kg"])
    if production is not None and not production.empty \
            and {"animal_id", "year", "mean_milk_yield_kg"}.issubset(production.columns):
        yield_lookup = (production[["animal_id", "year", "mean_milk_yield_kg"]]
                        .dropna(subset=["mean_milk_yield_kg"])
                        .copy())
        yield_lookup["year"] = yield_lookup["year"].astype(int)

    predictors = [
        ("thi",  "THI breakpoint",  "thi_breakpoint",  "thi_converged",
         "thi_breakpoint_se"),
        ("temp", "Barn temp breakpoint", "temp_breakpoint", "temp_converged",
         "temp_breakpoint_se"),
    ]

    def _pack(label, mode, result):
        if result is None:
            return None
        return {"predictor": label, "mode": mode, **result}

    for key, label, bp_col, conv_col, se_col in predictors:
        if bp_col not in bs_results.columns or conv_col not in bs_results.columns:
            continue
        conv = bs_results[bs_results[conv_col] == True].dropna(subset=[bp_col]).copy()
        if conv.empty:
            continue
        # Attach DIM, parity (lactation_nr) and mean yield once at the
        # cohort level so every downstream mode draws from the same
        # frame.  Each mode then chooses which subset / which
        # covariates to use; the cohort filter (≥ 2 converged years
        # per animal) is re-imposed separately for each mode after
        # any NaN-row drops.
        cohort = conv
        if not dim_lookup.empty:
            cohort = cohort.merge(dim_lookup, on=["animal_id", "year"], how="left")
        if not yield_lookup.empty:
            cohort = cohort.merge(yield_lookup, on=["animal_id", "year"], how="left")

        counts = conv.groupby("animal_id").size()
        keep = counts[counts >= 2].index
        sub = conv[conv["animal_id"].isin(keep)].copy()
        if sub.empty:
            out_rows.append({
                "predictor": label, "mode": "raw",
                "icc": np.nan, "ci_lower": np.nan, "ci_upper": np.nan,
                "f": np.nan, "df1": 0, "df2": 0, "p": np.nan,
                "n_animals": 0, "n_obs": 0, "k_mean": np.nan,
            })
            continue

        # ── Mode 1: raw ICC on the full multi-year cohort ───
        out_rows.append(_pack(label, "raw",
                              _icc_one_way(sub[bp_col].to_numpy(),
                                           sub["animal_id"].to_numpy())))

        # ── Mode 2 + 3: SE-cohort raw + measurement-corrected ─
        sub_se = pd.DataFrame()
        if se_col in sub.columns:
            sub_se = sub.dropna(subset=[se_col])
            counts_se = sub_se.groupby("animal_id").size()
            sub_se = sub_se[sub_se["animal_id"].isin(
                counts_se[counts_se >= 2].index)]
        if len(sub_se) >= 6:
            out_rows.append(_pack(label, "raw (SE-cohort)",
                                  _icc_one_way(sub_se[bp_col].to_numpy(),
                                               sub_se["animal_id"].to_numpy())))
            out_rows.append(_pack(label, "measurement-corrected",
                                  _icc_measurement_corrected(
                                      sub_se[bp_col].to_numpy(),
                                      sub_se["animal_id"].to_numpy(),
                                      sub_se[se_col].to_numpy())))

        # Helper: residualise on parity+DIM (+optional yield) using
        # the full converged cohort as the regression sample, then
        # ICC the residuals of the multi-year subset (optionally
        # SE-restricted, optionally measurement-corrected).
        def _residual_modes(use_yield: bool) -> None:
            tag = "parity + DIM"
            if use_yield:
                tag += " + yield"
            if "lactation_nr" not in cohort.columns or "mean_dim" not in cohort.columns:
                return
            if use_yield and "mean_milk_yield_kg" not in cohort.columns:
                return
            yield_arr = (cohort["mean_milk_yield_kg"].to_numpy()
                         if use_yield else None)
            cohort_local = cohort.copy()
            cohort_local["resid"] = _residualise_breakpoint(
                cohort_local[bp_col].to_numpy(),
                cohort_local["lactation_nr"].to_numpy(),
                cohort_local["mean_dim"].to_numpy(),
                mean_yield=yield_arr,
            )
            if cohort_local["resid"].isna().all():
                return

            # Mode N: residual ICC on multi-year subset (no SE filter)
            sub2 = cohort_local[cohort_local["animal_id"].isin(keep)] \
                       .dropna(subset=["resid"])
            c2 = sub2.groupby("animal_id").size()
            sub2 = sub2[sub2["animal_id"].isin(c2[c2 >= 2].index)]
            if not sub2.empty:
                out_rows.append(_pack(label, f"{tag} residual",
                                      _icc_one_way(sub2["resid"].to_numpy(),
                                                   sub2["animal_id"].to_numpy())))

            # Mode N+1: residualised AND measurement-corrected on the
            # SE-cohort.  This is the "fully corrected" estimate —
            # systematic parity / lactation-stage / yield trend AND
            # iid fit noise both subtracted.
            if se_col not in cohort_local.columns:
                return
            sub3 = cohort_local[cohort_local["animal_id"].isin(keep)] \
                       .dropna(subset=["resid", se_col])
            c3 = sub3.groupby("animal_id").size()
            sub3 = sub3[sub3["animal_id"].isin(c3[c3 >= 2].index)]
            if len(sub3) >= 6:
                out_rows.append(_pack(label, f"{tag} + meas-corr",
                                      _icc_measurement_corrected(
                                          sub3["resid"].to_numpy(),
                                          sub3["animal_id"].to_numpy(),
                                          sub3[se_col].to_numpy())))

        _residual_modes(use_yield=False)
        if not yield_lookup.empty:
            _residual_modes(use_yield=True)

    return pd.DataFrame(out_rows)


def compute_stability(bs_results: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    """ICC and paired data for repeat animals.

    Returns:
        (pairs DataFrame, ICC value).
    """
    conv = bs_results[bs_results["thi_converged"] == True]
    counts = conv.groupby("animal_id").size()
    repeats = counts[counts >= 2].index

    if len(repeats) < 3:
        return pd.DataFrame(), np.nan

    repeat_data = conv[conv["animal_id"].isin(repeats)].sort_values(["animal_id", "year"])
    groups = [g["thi_breakpoint"].values for _, g in repeat_data.groupby("animal_id")]
    all_vals = repeat_data["thi_breakpoint"].dropna()
    grand_mean = all_vals.mean()
    n_animals = len(groups)
    k_mean = len(all_vals) / n_animals

    ss_b = sum(len(g) * (g.mean() - grand_mean) ** 2 for g in groups)
    ss_w = sum(np.sum((g - g.mean()) ** 2) for g in groups)
    ms_b = ss_b / max(n_animals - 1, 1)
    ms_w = ss_w / max(len(all_vals) - n_animals, 1)
    denom = ms_b + (k_mean - 1) * ms_w
    icc = (ms_b - ms_w) / denom if denom > 0 else np.nan

    pairs = []
    for aid, grp in repeat_data.groupby("animal_id"):
        sg = grp.sort_values("year")
        for i in range(len(sg) - 1):
            pairs.append({
                "animal_id": aid,
                "year_1": int(sg.iloc[i]["year"]),
                "year_2": int(sg.iloc[i + 1]["year"]),
                "thi_bp_1": sg.iloc[i]["thi_breakpoint"],
                "thi_bp_2": sg.iloc[i + 1]["thi_breakpoint"],
                "temp_bp_1": sg.iloc[i]["temp_breakpoint"],
                "temp_bp_2": sg.iloc[i + 1]["temp_breakpoint"],
            })

    return pd.DataFrame(pairs), icc


# ─────────────────────────────────────────────────────────────
#  « summary table »
# ─────────────────────────────────────────────────────────────

def make_summary_table(bs: pd.DataFrame) -> pd.DataFrame:
    """Generate Table 1 for the manuscript."""
    years = sorted(bs["year"].dropna().unique().astype(int))
    rows = []
    for year in [*years, "All"]:
        yr = bs if year == "All" else bs[bs["year"] == year]
        for prefix, label in [
            ("thi", "THI→body_temp"),
            ("temp", "barn_temp→body_temp"),
            ("resp_thi", "THI→resp"),
            ("resp_temp", "barn_temp→resp"),
        ]:
            conv = yr[yr[f"{prefix}_converged"] == True]
            bp_col = f"{prefix}_breakpoint"
            rows.append({
                "year": year,
                "model": label,
                "n_total": len(yr),
                "n_converged": len(conv),
                "bp_median": conv[bp_col].median() if len(conv) else np.nan,
                "bp_q25": conv[bp_col].quantile(0.25) if len(conv) else np.nan,
                "bp_q75": conv[bp_col].quantile(0.75) if len(conv) else np.nan,
                "bp_min": conv[bp_col].min() if len(conv) else np.nan,
                "bp_max": conv[bp_col].max() if len(conv) else np.nan,
            })
    return pd.DataFrame(rows)

def _run_longitudinal_tests(bs: pd.DataFrame, d: Path) -> None:
    """Fisher resampling tests on longitudinal breakpoints."""
    from rerandomstats import FisherResamplingTest

    from digimuh.console import kv, result_table, stars_styled

    for bp_col, conv_col, label in [
        ("thi_breakpoint", "thi_converged", "THI"),
        ("temp_breakpoint", "temp_converged", "Barn temp"),
    ]:
        conv = bs[bs[conv_col] == True].dropna(subset=[bp_col])
        if conv.empty:
            continue

        # Animals in 2+ years
        year_counts = conv.groupby("animal_id")["year"].nunique()
        repeat_ids = year_counts[year_counts >= 2].index
        if len(repeat_ids) < 5:
            kv(f"{label} repeat animals", f"{len(repeat_ids)} (too few)")
            continue

        repeat = conv[conv["animal_id"].isin(repeat_ids)].copy()
        repeat = repeat.sort_values(["animal_id", "year"])
        first_bp = repeat.groupby("animal_id").first()[bp_col]
        repeat["bp_change"] = repeat[bp_col] - repeat["animal_id"].map(first_bp)
        years = sorted(repeat["year"].unique().astype(int))

        # Pairwise year comparisons (absolute breakpoints)
        test_rows = []
        raw_ps = []
        for i, y1 in enumerate(years):
            for y2 in years[i + 1:]:
                # Paired: same animals in both years
                ids_both = set(repeat[repeat["year"] == y1]["animal_id"]) & \
                           set(repeat[repeat["year"] == y2]["animal_id"])
                if len(ids_both) < 5:
                    continue
                d1 = repeat[
                    (repeat["year"] == y1) & (repeat["animal_id"].isin(ids_both))
                ][bp_col].tolist()
                d2 = repeat[
                    (repeat["year"] == y2) & (repeat["animal_id"].isin(ids_both))
                ][bp_col].tolist()
                p = FisherResamplingTest(
                    data_a=d1, data_b=d2,
                    func="medianDiff", combination_n=20_000,
                ).main()
                test_rows.append([f"{y1} vs {y2}", len(ids_both),
                                  np.median(d1), np.median(d2),
                                  np.median(d2) - np.median(d1), p, ""])
                raw_ps.append(p)

        if test_rows:
            adj_ps = correct_pvalues_array(np.array(raw_ps), method="fdr_bh")
            for r, adj_p in zip(test_rows, adj_ps):
                r[5] = adj_p
                r[6] = stars_styled(p_to_stars(adj_p))
            result_table(
                f"{label} breakpoints: pairwise year comparisons (BH-FDR)",
                ["Comparison", "n paired", "Median Y1", "Median Y2",
                 "Diff", "p adj", "Sig."],
                test_rows,
            )

        # Relative change from first year: test if different from 0
        for year in years[1:]:  # skip first year (change = 0)
            changes = repeat[repeat["year"] == year]["bp_change"].dropna()
            if len(changes) >= 5:
                # Test median change vs zero
                p = FisherResamplingTest(
                    data_a=changes.tolist(),
                    data_b=[0.0] * len(changes),
                    func="medianDiff", combination_n=20_000,
                ).main()
                kv(f"{label} {year} change from baseline",
                   f"median = {changes.median():.1f}, p = {p:.4f} {p_to_stars(p)}")

    # Save longitudinal test results alongside stability
    log.info("  Longitudinal tests complete.")

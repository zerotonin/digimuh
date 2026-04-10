#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  Analysis 0b — Statistical analysis ¹                        ║
# ║  « broken-stick fits, correlations, BH-FDR tests »          ║
# ╚══════════════════════════════════════════════════════════════╝
"""Reads CSVs from ``analysis_00a_extract`` and runs all statistics.

Outputs (in ``--data`` directory)::

    broken_stick_results.csv   Per-animal breakpoints (rumen + resp)
    spearman_correlations.csv  Per-animal Spearman rₛ
    behavioural_response.csv   Below/above breakpoint means per animal
    statistical_tests.csv      All Wilcoxon tests with BH-FDR correction
    breakpoint_stability.csv   Repeat-animal year pairs + ICC
    summary_table.csv          Table 1 for manuscript

Usage::

    digimuh-stats --data results/broken_stick

¹ Analysis led by Dr. med. vet. Gundula Hoffmann, ATB Potsdam.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr

log = logging.getLogger("digimuh.stats")


# ─────────────────────────────────────────────────────────────
#  « broken-stick regression model »
# ─────────────────────────────────────────────────────────────

def broken_stick_fit(
    x: np.ndarray, y: np.ndarray,
    x_range: tuple[float, float] | None = None,
    n_grid: int = 200,
) -> dict:
    """Fit a two-segment piecewise linear model with one breakpoint.

    Model: flat-ish below, steeper above.  Biological constraint:
    ``slope_above > slope_below`` and ``slope_above > 0``.

    Args:
        x: Predictor (THI or barn temp).
        y: Response (body temp or resp rate).
        x_range: Search bounds for breakpoint.
        n_grid: Grid resolution for initial search.

    Returns:
        Dict with breakpoint, slopes, R², convergence status.
    """
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)
    if n < 30:
        return {"breakpoint": np.nan, "converged": False, "n": n}

    if x_range is None:
        x_range = (np.percentile(x, 10), np.percentile(x, 90))
    if x_range[0] >= x_range[1]:
        return {"breakpoint": np.nan, "converged": False, "n": n}

    def rss_at_bp(bp):
        left, right = x <= bp, x > bp
        if left.sum() < 10 or right.sum() < 10:
            return np.inf
        try:
            Xl = np.column_stack([np.ones(left.sum()), x[left]])
            cl, *_ = np.linalg.lstsq(Xl, y[left], rcond=None)
            Xr = np.column_stack([np.ones(right.sum()), x[right]])
            cr, *_ = np.linalg.lstsq(Xr, y[right], rcond=None)
        except np.linalg.LinAlgError:
            return np.inf
        return np.sum((y[left] - Xl @ cl) ** 2) + np.sum((y[right] - Xr @ cr) ** 2)

    grid = np.linspace(x_range[0], x_range[1], n_grid)
    rss_vals = np.array([rss_at_bp(bp) for bp in grid])
    best_idx = np.argmin(rss_vals)
    if not np.isfinite(rss_vals[best_idx]):
        return {"breakpoint": np.nan, "converged": False, "n": n}

    lo = grid[max(0, best_idx - 2)]
    hi = grid[min(len(grid) - 1, best_idx + 2)]
    result = minimize_scalar(rss_at_bp, bounds=(lo, hi), method="bounded")
    bp = result.x
    final_rss = result.fun

    left, right = x <= bp, x > bp
    Xl = np.column_stack([np.ones(left.sum()), x[left]])
    cl, *_ = np.linalg.lstsq(Xl, y[left], rcond=None)
    Xr = np.column_stack([np.ones(right.sum()), x[right]])
    cr, *_ = np.linalg.lstsq(Xr, y[right], rcond=None)

    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_sq = 1 - final_rss / ss_tot if ss_tot > 0 else np.nan

    slope_below, slope_above = cl[1], cr[1]
    valid = slope_above > slope_below and slope_above > 0

    return {
        "breakpoint": bp if valid else np.nan,
        "slope_below": slope_below,
        "intercept_below": cl[0],
        "slope_above": slope_above,
        "intercept_above": cr[0],
        "r_squared": r_sq,
        "n": n,
        "n_below": int(left.sum()),
        "n_above": int(right.sum()),
        "converged": valid,
        "rejected_reason": (
            None if valid
            else f"slope_above ({slope_above:.4f}) <= slope_below ({slope_below:.4f})"
        ),
    }


# ─────────────────────────────────────────────────────────────
#  « Davies test for breakpoint existence »
# ─────────────────────────────────────────────────────────────

def davies_test(
    x: np.ndarray, y: np.ndarray,
    k: int = 10,
    x_range: tuple[float, float] | None = None,
) -> dict:
    """Davies (1987/2002) test for existence of a breakpoint.

    Tests H0: the relationship is a single straight line (no breakpoint)
    vs H1: there is an unknown breakpoint where the slope changes.

    The procedure evaluates k candidate breakpoints, computes a naive
    t-statistic for the slope-difference at each, takes the maximum,
    and corrects the p-value for the search using Davies' upper bound.

    Args:
        x: Predictor variable (THI or barn temp).
        y: Response variable (body temp or resp rate).
        k: Number of evaluation points (default 10, as in R segmented).
        x_range: Search range for candidate breakpoints.

    Returns:
        Dict with p_value, best_at (candidate with strongest signal),
        max_statistic, and n_sign_changes.

    References:
        Davies RB (1987) Biometrika 74:33-43.
        Davies RB (2002) Biometrika 89:484-489.
    """
    from scipy.stats import t as t_dist, norm

    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)
    if n < 30:
        return {"p_value": np.nan, "best_at": np.nan, "max_statistic": np.nan, "n": n}

    if x_range is None:
        x_range = (np.percentile(x, 5), np.percentile(x, 95))

    # k evaluation points spanning the predictor range
    psi_candidates = np.linspace(x_range[0], x_range[1], k)

    t_stats = np.full(k, np.nan)
    for j, psi in enumerate(psi_candidates):
        # Segmented variable: (x - psi)_+ = max(x - psi, 0)
        seg_var = np.maximum(x - psi, 0.0)

        # Design matrix: [1, x, (x-psi)+]
        X_design = np.column_stack([np.ones(n), x, seg_var])

        # OLS fit
        try:
            beta, residuals, rank, sv = np.linalg.lstsq(X_design, y, rcond=None)
        except np.linalg.LinAlgError:
            continue

        if rank < 3:
            continue

        # Residuals and variance
        y_hat = X_design @ beta
        resid = y - y_hat
        df = n - 3
        if df <= 0:
            continue
        sigma2 = np.sum(resid ** 2) / df

        # Standard error of the slope-difference coefficient (beta[2])
        try:
            XtX_inv = np.linalg.inv(X_design.T @ X_design)
        except np.linalg.LinAlgError:
            continue

        se_delta = np.sqrt(sigma2 * XtX_inv[2, 2])
        if se_delta <= 0:
            continue

        t_stats[j] = beta[2] / se_delta

    # Remove NaNs
    valid_t = t_stats[np.isfinite(t_stats)]
    if len(valid_t) < 2:
        return {"p_value": np.nan, "best_at": np.nan, "max_statistic": np.nan, "n": n}

    # Maximum absolute t-statistic
    abs_t = np.abs(valid_t)
    M = np.max(abs_t)
    best_idx = np.nanargmax(np.abs(t_stats))
    best_psi = psi_candidates[best_idx]

    # Count sign changes in the t-statistic sequence
    signs = np.sign(valid_t)
    sign_changes = np.sum(np.abs(np.diff(signs)) > 0)

    # Davies (1987) approximate upper bound for p-value
    # p <= Pr(|T| > M) + V * phi(M)
    # where V = number of sign changes, phi = standard normal density
    df = n - 3
    if df > 300:
        # Large sample: use normal approximation (exact gamma overflows)
        p_naive = 2.0 * (1.0 - norm.cdf(M))
        phi_M = norm.pdf(M)
    else:
        p_naive = 2.0 * (1.0 - t_dist.cdf(M, df))
        phi_M = t_dist.pdf(M, df)

    p_davies = p_naive + sign_changes * phi_M
    p_davies = min(p_davies, 1.0)

    return {
        "p_value": p_davies,
        "best_at": best_psi,
        "max_statistic": M,
        "n_sign_changes": int(sign_changes),
        "n": n,
    }


# ─────────────────────────────────────────────────────────────
#  « Pseudo-Score test (Muggeo 2016) »
# ─────────────────────────────────────────────────────────────

def pscore_test(
    x: np.ndarray, y: np.ndarray,
    k: int = 10,
    x_range: tuple[float, float] | None = None,
) -> dict:
    """Pseudo-Score test for breakpoint existence (Muggeo 2016).

    More powerful than Davies test for detecting a single changepoint.
    Averages the segmented variable over k candidate breakpoints and
    tests whether this averaged term is significant in the augmented
    linear model.

    The key insight: under H0 (no breakpoint), the averaged segmented
    variable phi_bar has zero coefficient.  Under H1 (one breakpoint),
    phi_bar picks up the slope change regardless of the true breakpoint
    location.

    Args:
        x: Predictor variable (THI or barn temp).
        y: Response variable (body temp or resp rate).
        k: Number of evaluation points.
        x_range: Search range for candidate breakpoints.

    Returns:
        Dict with p_value, t_statistic, best_at, and df.

    References:
        Muggeo VMR (2016) J Stat Comput Simul 86:3059-3067.
    """
    from scipy.stats import t as t_dist

    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)
    if n < 30:
        return {"p_value": np.nan, "t_statistic": np.nan, "n": n}

    if x_range is None:
        x_range = (np.percentile(x, 5), np.percentile(x, 95))

    # k evaluation points
    psi_candidates = np.linspace(x_range[0], x_range[1], k)

    # Averaged segmented variable:
    # phi_bar_i = (1/k) * sum_j max(x_i - psi_j, 0)
    phi_bar = np.zeros(n)
    for psi in psi_candidates:
        phi_bar += np.maximum(x - psi, 0.0)
    phi_bar /= k

    # Augmented model: y = alpha + beta*x + gamma*phi_bar
    X_design = np.column_stack([np.ones(n), x, phi_bar])

    try:
        beta, residuals, rank, sv = np.linalg.lstsq(X_design, y, rcond=None)
    except np.linalg.LinAlgError:
        return {"p_value": np.nan, "t_statistic": np.nan, "n": n}

    if rank < 3:
        return {"p_value": np.nan, "t_statistic": np.nan, "n": n}

    y_hat = X_design @ beta
    resid = y - y_hat
    df = n - 3
    if df <= 0:
        return {"p_value": np.nan, "t_statistic": np.nan, "n": n}

    sigma2 = np.sum(resid ** 2) / df

    try:
        XtX_inv = np.linalg.inv(X_design.T @ X_design)
    except np.linalg.LinAlgError:
        return {"p_value": np.nan, "t_statistic": np.nan, "n": n}

    se_gamma = np.sqrt(sigma2 * XtX_inv[2, 2])
    if se_gamma <= 0:
        return {"p_value": np.nan, "t_statistic": np.nan, "n": n}

    t_stat = beta[2] / se_gamma
    p_value = 2.0 * (1.0 - t_dist.cdf(abs(t_stat), df))

    # Also find which candidate psi gave the strongest individual signal
    # (for reporting, not for the test itself)
    best_t = 0.0
    best_psi = np.nan
    for psi in psi_candidates:
        seg = np.maximum(x - psi, 0.0)
        Xd = np.column_stack([np.ones(n), x, seg])
        try:
            b, *_ = np.linalg.lstsq(Xd, y, rcond=None)
            yh = Xd @ b
            r = y - yh
            s2 = np.sum(r ** 2) / (n - 3)
            inv_xx = np.linalg.inv(Xd.T @ Xd)
            se = np.sqrt(s2 * inv_xx[2, 2])
            if se > 0 and abs(b[2] / se) > abs(best_t):
                best_t = b[2] / se
                best_psi = psi
        except (np.linalg.LinAlgError, ValueError):
            continue

    return {
        "p_value": p_value,
        "t_statistic": t_stat,
        "best_at": best_psi,
        "df": df,
        "n": n,
    }


# ─────────────────────────────────────────────────────────────
#  « four-parameter logistic (Hill) fit »
# ─────────────────────────────────────────────────────────────

def hill_fit(
    x: np.ndarray, y: np.ndarray,
    x_range: tuple[float, float] | None = None,
) -> dict:
    """Fit a four-parameter logistic (Hill) dose-response model.

    Model::

        y = y_min + (y_max - y_min) / (1 + (EC50 / x)^n)

    The EC50 is the predictor value at half-maximal response.  The Hill
    coefficient *n* captures steepness: large *n* = sharp switch (like
    rumen temperature), small *n* = gradual transition (like respiration).

    In addition to EC50, we compute the **lower bend point** following
    Sebaugh & McCray (2003, Pharmaceutical Statistics 2:167-174).  This
    is the x-value where the second derivative of the Hill curve equals
    zero on the lower side, marking the transition from the baseline
    plateau into the rising phase::

        x_bend_lower = EC50 * ((n - 1) / (n + 1))^(1/n)    for n > 1

    This is the onset of the nonlinear rise and is directly comparable
    to the broken-stick breakpoint.

    Args:
        x: Predictor (THI or barn temperature).
        y: Response (body temp or respiration rate).
        x_range: Plausible range for EC50 initial guess.

    Returns:
        Dict with ec50, hill_n, y_min, y_max, lower_bend, r_squared,
        aic, converged.

    References:
        Sebaugh JL, McCray PD (2003) Defining the linear portion of a
        sigmoid-shaped curve: bend points. Pharm Stat 2:167-174.
    """
    from scipy.optimize import curve_fit

    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)
    if n < 30:
        return {"ec50": np.nan, "converged": False, "n": n}

    if x_range is None:
        x_range = (np.percentile(x, 10), np.percentile(x, 90))

    def hill_func(xv, y_min, y_max, ec50, hill_n):
        """Four-parameter logistic (Hill equation)."""
        ratio = np.clip(ec50 / np.maximum(xv, 1e-10), 1e-10, 1e10)
        return y_min + (y_max - y_min) / (1.0 + np.power(ratio, hill_n))

    # Initial guesses from data quantiles
    y_lo = np.percentile(y, 10)
    y_hi = np.percentile(y, 90)
    ec50_init = np.mean(x_range)

    # Bounds to keep parameters physiologically plausible
    bounds_lo = [np.min(y) - 5, np.median(y), x_range[0] * 0.5, 0.5]
    bounds_hi = [np.median(y), np.max(y) + 5, x_range[1] * 1.5, 30.0]

    try:
        popt, pcov = curve_fit(
            hill_func, x, y,
            p0=[y_lo, y_hi, ec50_init, 2.0],
            bounds=(bounds_lo, bounds_hi),
            maxfev=10000,
        )
    except (RuntimeError, ValueError):
        return {"ec50": np.nan, "converged": False, "n": n}

    y_min_fit, y_max_fit, ec50, hill_n = popt

    # Goodness of fit
    y_pred = hill_func(x, *popt)
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    # AIC (4 parameters + sigma)
    k_params = 5
    aic = n * np.log(ss_res / n) + 2 * k_params if n > k_params else np.nan

    # Lower bend point (Sebaugh & McCray 2003)
    # x_bend = EC50 * ((n-1)/(n+1))^(1/n)  for n > 1
    if hill_n > 1.0:
        lower_bend = ec50 * ((hill_n - 1.0) / (hill_n + 1.0)) ** (1.0 / hill_n)
    else:
        # For n <= 1 the curve has no inflection in the lower portion;
        # fall back to EC10 as a conventional onset marker
        lower_bend = ec50 * (0.10 / 0.90) ** (1.0 / hill_n)

    # Sanity: bend point should be within or near the data range
    x_lo, x_hi = np.min(x), np.max(x)
    plausible = (x_lo - 10) <= lower_bend <= x_hi

    return {
        "ec50": ec50,
        "hill_n": hill_n,
        "y_min": y_min_fit,
        "y_max": y_max_fit,
        "lower_bend": lower_bend if plausible else np.nan,
        "r_squared": r_sq,
        "aic": aic,
        "n": n,
        "converged": True,
        "bend_plausible": plausible,
    }


# ─────────────────────────────────────────────────────────────
#  « Benjamini-Hochberg FDR »
# ─────────────────────────────────────────────────────────────

def benjamini_hochberg(p_values: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """Benjamini-Hochberg FDR correction.

    Args:
        p_values: Array of raw p-values.
        alpha: Target FDR level (for reference; returns adjusted p).

    Returns:
        Array of BH-adjusted p-values.
    """
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return p

    # Handle NaNs
    notnan = ~np.isnan(p)
    p_valid = p[notnan]
    m = len(p_valid)
    if m == 0:
        return p

    order = np.argsort(p_valid)
    ranked = np.empty(m)
    ranked[order] = np.arange(1, m + 1)

    adjusted = p_valid * m / ranked
    # Enforce monotonicity (step-down)
    adjusted[order] = np.minimum.accumulate(adjusted[order[::-1]])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    result = np.full(n, np.nan)
    result[notnan] = adjusted
    return result


def p_to_stars(p: float) -> str:
    """Convert p-value to significance stars.

    Returns:
        ``'***'`` if p < 0.001, ``'**'`` if p < 0.01,
        ``'*'`` if p < 0.05, ``'n.s.'`` otherwise.
    """
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


# ─────────────────────────────────────────────────────────────
#  « broken-stick fits for all animals »
# ─────────────────────────────────────────────────────────────

def run_broken_stick_fits(
    rumen: pd.DataFrame, resp: pd.DataFrame,
) -> pd.DataFrame:
    """Fit broken-stick models per animal-year.

    Args:
        rumen: rumen_barn.csv DataFrame.
        resp: respiration_barn.csv DataFrame.

    Returns:
        One row per animal-year with breakpoint results.
    """
    results = []
    groups = rumen.groupby(["animal_id", "year", "date_enter"])
    total = len(groups)

    for i, ((aid, year, enter), grp) in enumerate(groups):
        if (i + 1) % 20 == 0 or i == 0:
            log.info("  [%d/%d] Fitting animal %d (%s)", i + 1, total, aid, year)

        n_readings = len(grp)
        if n_readings < 50:
            results.append({
                "animal_id": aid, "year": year, "date_enter": enter,
                "n_readings": n_readings,
                "thi_breakpoint": np.nan, "thi_converged": False,
                "temp_breakpoint": np.nan, "temp_converged": False,
                "resp_thi_breakpoint": np.nan, "resp_thi_converged": False,
                "resp_temp_breakpoint": np.nan, "resp_temp_converged": False,
                "comment": f"insufficient data ({n_readings})",
            })
            continue

        # Body temp fits
        thi_fit = broken_stick_fit(
            grp["barn_thi"].values, grp["body_temp"].values, x_range=(45, 80))
        temp_fit = broken_stick_fit(
            grp["barn_temp"].values, grp["body_temp"].values, x_range=(5, 35))

        # Breakpoint existence tests (body temp)
        thi_davies = davies_test(
            grp["barn_thi"].values, grp["body_temp"].values, x_range=(45, 80))
        thi_pscore = pscore_test(
            grp["barn_thi"].values, grp["body_temp"].values, x_range=(45, 80))
        temp_davies = davies_test(
            grp["barn_temp"].values, grp["body_temp"].values, x_range=(5, 35))
        temp_pscore = pscore_test(
            grp["barn_temp"].values, grp["body_temp"].values, x_range=(5, 35))

        # Hill / logistic fits (body temp)
        thi_hill = hill_fit(
            grp["barn_thi"].values, grp["body_temp"].values, x_range=(45, 80))
        temp_hill = hill_fit(
            grp["barn_temp"].values, grp["body_temp"].values, x_range=(5, 35))

        # Respiration fits
        resp_grp = resp[
            (resp["animal_id"] == aid) & (resp["year"] == year)
        ] if not resp.empty else pd.DataFrame()
        has_resp = len(resp_grp) >= 50

        if has_resp:
            resp_thi_fit = broken_stick_fit(
                resp_grp["barn_thi"].values, resp_grp["resp_rate"].values,
                x_range=(45, 80))
            resp_temp_fit = broken_stick_fit(
                resp_grp["barn_temp"].values, resp_grp["resp_rate"].values,
                x_range=(5, 35))
            # Breakpoint existence tests (respiration)
            resp_thi_davies = davies_test(
                resp_grp["barn_thi"].values, resp_grp["resp_rate"].values,
                x_range=(45, 80))
            resp_thi_pscore = pscore_test(
                resp_grp["barn_thi"].values, resp_grp["resp_rate"].values,
                x_range=(45, 80))
            resp_temp_davies = davies_test(
                resp_grp["barn_temp"].values, resp_grp["resp_rate"].values,
                x_range=(5, 35))
            resp_temp_pscore = pscore_test(
                resp_grp["barn_temp"].values, resp_grp["resp_rate"].values,
                x_range=(5, 35))
            # Hill / logistic fits (respiration)
            resp_thi_hill = hill_fit(
                resp_grp["barn_thi"].values, resp_grp["resp_rate"].values,
                x_range=(45, 80))
            resp_temp_hill = hill_fit(
                resp_grp["barn_temp"].values, resp_grp["resp_rate"].values,
                x_range=(5, 35))
        else:
            resp_thi_fit = {"breakpoint": np.nan, "converged": False, "n": 0}
            resp_temp_fit = {"breakpoint": np.nan, "converged": False, "n": 0}
            _empty_test = {"p_value": np.nan}
            resp_thi_davies = resp_thi_pscore = _empty_test
            resp_temp_davies = resp_temp_pscore = _empty_test
            _empty_hill = {"ec50": np.nan, "hill_n": np.nan,
                           "lower_bend": np.nan, "r_squared": np.nan,
                           "aic": np.nan, "converged": False}
            resp_thi_hill = resp_temp_hill = _empty_hill

        rec = {
            "animal_id": aid, "year": int(year), "date_enter": enter,
            "n_readings": n_readings,
            "n_resp_readings": len(resp_grp),
            # Body temp vs THI
            "thi_breakpoint": thi_fit["breakpoint"],
            "thi_slope_below": thi_fit.get("slope_below"),
            "thi_slope_above": thi_fit.get("slope_above"),
            "thi_r_squared": thi_fit.get("r_squared"),
            "thi_converged": thi_fit["converged"],
            "thi_davies_p": thi_davies["p_value"],
            "thi_pscore_p": thi_pscore["p_value"],
            "thi_hill_ec50": thi_hill.get("ec50"),
            "thi_hill_n": thi_hill.get("hill_n"),
            "thi_hill_bend": thi_hill.get("lower_bend"),
            "thi_hill_r2": thi_hill.get("r_squared"),
            "thi_hill_converged": thi_hill.get("converged", False),
            # Body temp vs barn temp
            "temp_breakpoint": temp_fit["breakpoint"],
            "temp_slope_below": temp_fit.get("slope_below"),
            "temp_slope_above": temp_fit.get("slope_above"),
            "temp_r_squared": temp_fit.get("r_squared"),
            "temp_converged": temp_fit["converged"],
            "temp_davies_p": temp_davies["p_value"],
            "temp_pscore_p": temp_pscore["p_value"],
            "temp_hill_ec50": temp_hill.get("ec50"),
            "temp_hill_n": temp_hill.get("hill_n"),
            "temp_hill_bend": temp_hill.get("lower_bend"),
            "temp_hill_r2": temp_hill.get("r_squared"),
            "temp_hill_converged": temp_hill.get("converged", False),
            # Resp vs THI
            "resp_thi_breakpoint": resp_thi_fit["breakpoint"],
            "resp_thi_slope_below": resp_thi_fit.get("slope_below"),
            "resp_thi_slope_above": resp_thi_fit.get("slope_above"),
            "resp_thi_r_squared": resp_thi_fit.get("r_squared"),
            "resp_thi_converged": resp_thi_fit.get("converged", False),
            "resp_thi_davies_p": resp_thi_davies["p_value"],
            "resp_thi_pscore_p": resp_thi_pscore["p_value"],
            "resp_thi_hill_ec50": resp_thi_hill.get("ec50"),
            "resp_thi_hill_n": resp_thi_hill.get("hill_n"),
            "resp_thi_hill_bend": resp_thi_hill.get("lower_bend"),
            "resp_thi_hill_r2": resp_thi_hill.get("r_squared"),
            "resp_thi_hill_converged": resp_thi_hill.get("converged", False),
            # Resp vs barn temp
            "resp_temp_breakpoint": resp_temp_fit["breakpoint"],
            "resp_temp_r_squared": resp_temp_fit.get("r_squared"),
            "resp_temp_converged": resp_temp_fit.get("converged", False),
            "resp_temp_davies_p": resp_temp_davies["p_value"],
            "resp_temp_pscore_p": resp_temp_pscore["p_value"],
            "resp_temp_hill_ec50": resp_temp_hill.get("ec50"),
            "resp_temp_hill_n": resp_temp_hill.get("hill_n"),
            "resp_temp_hill_bend": resp_temp_hill.get("lower_bend"),
            "resp_temp_hill_r2": resp_temp_hill.get("r_squared"),
            "resp_temp_hill_converged": resp_temp_hill.get("converged", False),
        }
        results.append(rec)

    return pd.DataFrame(results)


# ─────────────────────────────────────────────────────────────
#  « Spearman correlations »
# ─────────────────────────────────────────────────────────────

def compute_spearman(rumen: pd.DataFrame, resp: pd.DataFrame) -> pd.DataFrame:
    """Per-animal Spearman correlations."""
    records = []
    for (aid, year), grp in rumen.groupby(["animal_id", "year"]):
        rec = {"animal_id": aid, "year": int(year), "n": len(grp)}
        for x_col, prefix in [("barn_thi", "thi"), ("barn_temp", "temp")]:
            sub = grp.dropna(subset=[x_col, "body_temp"])
            if len(sub) > 20:
                rs, p = spearmanr(sub[x_col], sub["body_temp"])
                rec[f"{prefix}_rs"] = rs
                rec[f"{prefix}_p"] = p
            else:
                rec[f"{prefix}_rs"] = np.nan
                rec[f"{prefix}_p"] = np.nan
        # Respiration Spearman
        resp_grp = resp[(resp["animal_id"] == aid) & (resp["year"] == year)] if not resp.empty else pd.DataFrame()
        for x_col, prefix in [("barn_thi", "resp_thi"), ("barn_temp", "resp_temp")]:
            sub = resp_grp.dropna(subset=[x_col, "resp_rate"]) if not resp_grp.empty else pd.DataFrame()
            if len(sub) > 20:
                rs, p = spearmanr(sub[x_col], sub["resp_rate"])
                rec[f"{prefix}_rs"] = rs
                rec[f"{prefix}_p"] = p
            else:
                rec[f"{prefix}_rs"] = np.nan
                rec[f"{prefix}_p"] = np.nan
        records.append(rec)
    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────
#  « below/above breakpoint comparisons »
# ─────────────────────────────────────────────────────────────

def compute_below_above(
    rumen: pd.DataFrame, resp: pd.DataFrame,
    bs_results: pd.DataFrame,
) -> pd.DataFrame:
    """Per-animal means below/above their individual THI breakpoint."""
    converged = bs_results[bs_results["thi_converged"] == True]
    records = []

    for _, row in converged.iterrows():
        aid = int(row["animal_id"])
        year = int(row["year"])
        bp = row["thi_breakpoint"]

        grp = rumen[(rumen["animal_id"] == aid) & (rumen["year"] == year)]
        if len(grp) < 30:
            continue

        below = grp[grp["barn_thi"] <= bp]
        above = grp[grp["barn_thi"] > bp]
        if len(below) < 10 or len(above) < 10:
            continue

        rec = {
            "animal_id": aid, "year": year,
            "thi_breakpoint": bp,
            "body_temp_below": below["body_temp"].mean(),
            "body_temp_above": above["body_temp"].mean(),
            "n_below": len(below),
            "n_above": len(above),
        }

        # Respiration
        resp_grp = resp[(resp["animal_id"] == aid) & (resp["year"] == year)] if not resp.empty else pd.DataFrame()
        if len(resp_grp) >= 20:
            rb = resp_grp[resp_grp["barn_thi"] <= bp]
            ra = resp_grp[resp_grp["barn_thi"] > bp]
            if len(rb) >= 5 and len(ra) >= 5:
                rec["resp_below"] = rb["resp_rate"].mean()
                rec["resp_above"] = ra["resp_rate"].mean()
            else:
                rec["resp_below"] = rec["resp_above"] = np.nan
        else:
            rec["resp_below"] = rec["resp_above"] = np.nan

        records.append(rec)

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────
#  « within-year Fisher resampling tests with BH-FDR »
# ─────────────────────────────────────────────────────────────

def run_statistical_tests(beh: pd.DataFrame) -> pd.DataFrame:
    """Run Fisher resampling tests within each year, BH-FDR corrected.

    Uses reRandomStats FisherResamplingTest with medianDiff as the test
    statistic and 20,000 permutations.

    Tests per year:
    - Body temp below vs above breakpoint

    Args:
        beh: behavioural_response DataFrame.

    Returns:
        DataFrame: one row per test, with raw p, adjusted p, stars.
    """
    from rerandomstats import FisherResamplingTest

    test_rows = []
    years = sorted(beh["year"].dropna().unique().astype(int))

    for year in years:
        yr = beh[beh["year"] == year]
        year_tests = []

        # Body temperature below vs above breakpoint
        paired = yr.dropna(subset=["body_temp_below", "body_temp_above"])
        if len(paired) >= 10:
            p = FisherResamplingTest(
                data_a=paired["body_temp_below"].tolist(),
                data_b=paired["body_temp_above"].tolist(),
                func="medianDiff",
                combination_n=20_000,
            ).main()
            year_tests.append({
                "year": year,
                "test": "body_temp below vs above (Fisher medianDiff)",
                "n": len(paired),
                "p_raw": p,
                "median_below": paired["body_temp_below"].median(),
                "median_above": paired["body_temp_above"].median(),
                "median_diff": (paired["body_temp_above"] - paired["body_temp_below"]).median(),
            })

        # Respiration (if columns exist and have data)
        if "resp_below" in yr.columns:
            paired_r = yr.dropna(subset=["resp_below", "resp_above"])
            if len(paired_r) >= 10:
                p = FisherResamplingTest(
                    data_a=paired_r["resp_below"].tolist(),
                    data_b=paired_r["resp_above"].tolist(),
                    func="medianDiff",
                    combination_n=20_000,
                ).main()
                year_tests.append({
                    "year": year,
                    "test": "respiration below vs above (Fisher medianDiff)",
                    "n": len(paired_r),
                    "p_raw": p,
                    "median_below": paired_r["resp_below"].median(),
                    "median_above": paired_r["resp_above"].median(),
                    "median_diff": (paired_r["resp_above"] - paired_r["resp_below"]).median(),
                })

        # BH-FDR within this year
        if year_tests:
            raw_ps = np.array([t["p_raw"] for t in year_tests])
            adj_ps = benjamini_hochberg(raw_ps)
            for t, adj_p in zip(year_tests, adj_ps):
                t["p_adj"] = adj_p
                t["stars"] = p_to_stars(adj_p)
            test_rows.extend(year_tests)

    return pd.DataFrame(test_rows)


# ─────────────────────────────────────────────────────────────
#  « cross-correlation / cross-covariance below/above bp »
# ─────────────────────────────────────────────────────────────

def compute_cross_correlation(
    rumen: pd.DataFrame, bs_results: pd.DataFrame,
    max_lag: int = 12,
) -> pd.DataFrame:
    """Cross-correlation and cross-covariance of climate vs rumen temp,
    computed separately below and above each animal's breakpoint.

    For each animal with a converged breakpoint (THI and barn temp), the
    rumen temperature time series is split at the breakpoint.  The
    normalised cross-correlation function (CCF) and raw cross-covariance
    are computed for lags -max_lag to +max_lag (in units of the 10-min
    sampling interval, so max_lag=12 = 2 hours).

    Positive lags mean the climate signal leads the rumen temperature
    response.

    Args:
        rumen: rumen_barn.csv DataFrame (must have timestamp column).
        bs_results: broken_stick_results.csv DataFrame.
        max_lag: Maximum lag in samples (default 12 = 2 hours).

    Returns:
        DataFrame with columns: animal_id, year, predictor (thi/temp),
        region (below/above), lag, xcorr, xcov.
    """
    records = []

    for prefix, env_col, conv_col, bp_col in [
        ("thi", "barn_thi", "thi_converged", "thi_breakpoint"),
        ("temp", "barn_temp", "temp_converged", "temp_breakpoint"),
    ]:
        converged = bs_results[bs_results[conv_col] == True]

        for _, row in converged.iterrows():
            aid = int(row["animal_id"])
            year = int(row["year"])
            bp = row[bp_col]

            grp = rumen[(rumen["animal_id"] == aid) & (rumen["year"] == year)].copy()
            if len(grp) < 50:
                continue

            # Sort by timestamp for temporal lag analysis
            grp = grp.sort_values("timestamp").reset_index(drop=True)

            x_full = grp[env_col].values
            y_full = grp["body_temp"].values

            for region, mask_fn in [
                ("below", lambda x: x <= bp),
                ("above", lambda x: x > bp),
            ]:
                mask = mask_fn(x_full)
                if mask.sum() < max_lag * 3:
                    continue

                x = x_full[mask]
                y = y_full[mask]

                # Demean
                xc = x - np.mean(x)
                yc = y - np.mean(y)

                sx = np.std(x, ddof=1)
                sy = np.std(y, ddof=1)
                n = len(x)

                if sx < 1e-10 or sy < 1e-10:
                    continue

                for lag in range(-max_lag, max_lag + 1):
                    if lag >= 0:
                        x_seg = xc[:n - lag] if lag > 0 else xc
                        y_seg = yc[lag:] if lag > 0 else yc
                    else:
                        x_seg = xc[-lag:]
                        y_seg = yc[:n + lag]

                    m = len(x_seg)
                    if m < 10:
                        continue

                    # Cross-covariance (raw)
                    xcov = np.sum(x_seg * y_seg) / (m - 1)
                    # Normalised cross-correlation
                    xcorr = xcov / (sx * sy)

                    records.append({
                        "animal_id": aid, "year": year,
                        "predictor": prefix,
                        "region": region,
                        "lag": lag,
                        "lag_minutes": lag * 10,
                        "xcorr": xcorr,
                        "xcov": xcov,
                        "n": m,
                    })

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────
#  « breakpoint stability (ICC) »
# ─────────────────────────────────────────────────────────────

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


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Statistical analysis for broken-stick")
    parser.add_argument("--data", type=Path, required=True,
                        help="Directory with CSVs from extraction step")
    parser.add_argument("--no-resp", action="store_true",
                        help="Skip respiration analysis (Frontiers paper: "
                             "rumen temperature only)")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(name)-20s │ %(message)s",
        datefmt="%H:%M:%S",
    )
    d = args.data

    log.info("Loading extracted CSVs from %s …", d)
    rumen = pd.read_csv(d / "rumen_barn.csv")
    resp_path = d / "respiration_barn.csv"
    if args.no_resp:
        resp = pd.DataFrame()
        log.info("  Respiration analysis SKIPPED (--no-resp)")
    else:
        resp = pd.read_csv(resp_path) if resp_path.exists() else pd.DataFrame()
    prod = pd.read_csv(d / "production.csv") if (d / "production.csv").exists() else pd.DataFrame()

    # 1. Broken-stick + Davies/pscore + Hill fits
    log.info("═" * 50)
    log.info("1. Broken-stick + Hill fits …")
    bs = run_broken_stick_fits(rumen, resp)
    # Merge production data
    if not prod.empty:
        bs = bs.merge(
            prod[["animal_id", "year", "mean_milk_yield_kg", "lactation_nr"]],
            on=["animal_id", "year"], how="left",
        )
    bs.to_csv(d / "broken_stick_results.csv", index=False)

    prefixes = [("thi", "THI→body"), ("temp", "Temp→body")]
    if not args.no_resp:
        prefixes += [("resp_thi", "THI→resp"), ("resp_temp", "Temp→resp")]

    for prefix, label in prefixes:
        n_conv = (bs[f"{prefix}_converged"] == True).sum()
        log.info("  %s: %d/%d converged", label, n_conv, len(bs))

    # Breakpoint existence test summary
    log.info("─" * 50)
    log.info("  Breakpoint existence tests (p < 0.05):")
    for prefix, label in prefixes:
        davies_col = f"{prefix}_davies_p"
        pscore_col = f"{prefix}_pscore_p"
        if davies_col not in bs.columns:
            continue
        n_valid = bs[davies_col].notna().sum()
        if n_valid == 0:
            continue
        n_davies_sig = (bs[davies_col] < 0.05).sum()
        n_pscore_sig = (bs[pscore_col] < 0.05).sum()
        log.info(
            "  %s: Davies %d/%d (%.0f%%), pscore %d/%d (%.0f%%)",
            label, n_davies_sig, n_valid, 100 * n_davies_sig / n_valid,
            n_pscore_sig, n_valid, 100 * n_pscore_sig / n_valid,
        )

    # Hill fit summary
    log.info("─" * 50)
    log.info("  Hill (4PL) fits — lower bend point (Sebaugh & McCray 2003):")
    for prefix, label in prefixes:
        hill_conv_col = f"{prefix}_hill_converged"
        hill_bend_col = f"{prefix}_hill_bend"
        if hill_conv_col not in bs.columns:
            continue
        n_hill_conv = (bs[hill_conv_col] == True).sum()
        bends = bs.loc[bs[hill_conv_col] == True, hill_bend_col].dropna()
        n_bend_valid = len(bends)
        if n_bend_valid > 0:
            log.info(
                "  %s: converged %d/%d (%.0f%%), bend valid %d, "
                "median bend=%.1f [IQR %.1f–%.1f]",
                label, n_hill_conv, len(bs), 100 * n_hill_conv / len(bs),
                n_bend_valid, bends.median(), bends.quantile(0.25),
                bends.quantile(0.75),
            )
        else:
            log.info("  %s: converged %d/%d (%.0f%%), no valid bend points",
                     label, n_hill_conv, len(bs), 100 * n_hill_conv / len(bs))

    # 2. Spearman correlations
    log.info("═" * 50)
    log.info("2. Spearman correlations …")
    spearman = compute_spearman(rumen, resp)
    spearman.to_csv(d / "spearman_correlations.csv", index=False)
    log.info("  Body temp vs THI: median rₛ = %.3f", spearman["thi_rs"].median())
    log.info("  Body temp vs barn temp: median rₛ = %.3f", spearman["temp_rs"].median())

    # 3. Below/above breakpoint
    log.info("═" * 50)
    log.info("3. Below/above breakpoint comparisons …")
    beh = compute_below_above(rumen, resp, bs)
    beh.to_csv(d / "behavioural_response.csv", index=False)
    log.info("  %d animals with paired data", len(beh))

    # 4. Fisher resampling tests with BH-FDR
    log.info("═" * 50)
    log.info("4. Fisher resampling tests (reRandomStats, BH-FDR corrected) …")
    tests = run_statistical_tests(beh)
    tests.to_csv(d / "statistical_tests.csv", index=False)
    for _, t in tests.iterrows():
        log.info(
            "  %d | %-45s | n=%3d | p_raw=%.4f | p_adj=%.4f | %s",
            t["year"], t["test"], t["n"], t["p_raw"], t["p_adj"], t["stars"],
        )

    # 5. Cross-correlation / cross-covariance below/above breakpoint
    log.info("═" * 50)
    log.info("5. Cross-correlation below/above breakpoint …")
    xcorr = compute_cross_correlation(rumen, bs)
    xcorr.to_csv(d / "cross_correlation.csv", index=False)
    for pred in ["thi", "temp"]:
        for region in ["below", "above"]:
            sub = xcorr[(xcorr["predictor"] == pred) & (xcorr["region"] == region)
                        & (xcorr["lag"] == 0)]
            if not sub.empty:
                log.info("  %s %s bp: median r(lag=0) = %.3f, n=%d animals",
                         pred.upper(), region, sub["xcorr"].median(),
                         sub["animal_id"].nunique())

    # 6. Breakpoint stability
    log.info("═" * 50)
    log.info("6. Breakpoint stability (repeat animals) …")
    pairs, icc = compute_stability(bs)
    if not pairs.empty:
        pairs.to_csv(d / "breakpoint_stability.csv", index=False)
    log.info("  ICC = %.3f, %d pairs from %d animals",
             icc, len(pairs), pairs["animal_id"].nunique() if not pairs.empty else 0)

    # 7. Summary table
    log.info("═" * 50)
    log.info("7. Summary table …")
    summary = make_summary_table(bs)
    summary.to_csv(d / "summary_table.csv", index=False)

    log.info("═" * 50)
    log.info("All statistics complete.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — fitting                                              ║
# ║  « threshold detection models for climate–physiology coupling » ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Broken-stick regression, Davies test, pseudo-Score test,       ║
# ║  and Hill/4-parameter logistic fit.  Each function takes raw    ║
# ║  x/y arrays and returns a result dict — no I/O, no side        ║
# ║  effects.                                                       ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Statistical model-fitting functions for breakpoint detection."""

from __future__ import annotations

import logging

import numpy as np
from scipy.optimize import minimize_scalar

from digimuh.constants import GRID_STEPS

log = logging.getLogger("digimuh.fitting")

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

        se_delta = np.sqrt(max(sigma2 * XtX_inv[2, 2], 0.0))
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

    se_gamma = np.sqrt(max(sigma2 * XtX_inv[2, 2], 0.0))
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
            se = np.sqrt(max(s2 * inv_xx[2, 2], 0.0))
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

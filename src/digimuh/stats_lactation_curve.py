#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_lactation_curve                                ║
# ║  « Wood (1967) γ fits + DIM-adjusted milk-yield residuals »   ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Utilities for expressing each cow-day's daily milk yield as    ║
# ║  a deviation from the animal's expected yield given her days    ║
# ║  in milk (DIM).  The residuals decorrelate yield from the       ║
# ║  lactation curve — the dominant nuisance variance in raw        ║
# ║  daily yield — so downstream analyses (milk-yield classifica-  ║
# ║  tion, broken-stick stratification) can compare cows fairly    ║
# ║  across lactation stages.                                       ║
# ║                                                                 ║
# ║  Reusable by the Frontiers-2026 pipeline.  Standalone: reads    ║
# ║  daily_milk_yield.csv and calvings.csv (produced by extract),   ║
# ║  writes daily_milk_yield_wood.csv alongside.                    ║
# ╚══════════════════════════════════════════════════════════════════╝
"""DIM-adjusted milk yield via per-lactation Wood (1967) curves.

Wood model::

    y(t) = a · t^b · exp(-c · t)

fitted in log space by OLS (``ln y = ln a + b ln t − c t``).  Each
cow-lactation gets its own (a, b, c); lactations with too few points
fall back to a parity-pooled Wood curve.  The resulting residuals
``y − ŷ`` are approximately symmetric, freeing downstream
stratifications (e.g. terciles) from DIM × parity confounding.

References:
    Wood, P.D.P. (1967) Algebraic model of the lactation curve in
        cattle. *Nature* 216: 164–165. doi:10.1038/216164a0
    Wilmink, J.B.M. (1987) Adjustment of test-day milk, fat and
        protein yield for age, season and stage of lactation.
        *Livest Prod Sci* 16: 335–348.
    Macciotta, N.P.P. et al. (2011) Mathematical description of
        lactation curves in dairy cattle. *Ital J Anim Sci* 10: e51.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.paths import resolve_input, resolve_output

log = logging.getLogger("digimuh.lactation")


# ─────────────────────────────────────────────────────────────
#  « Defaults »
# ─────────────────────────────────────────────────────────────

DIM_MIN: int = 5
"""Minimum DIM to include (fresh-cow noise dominates 0–4 d)."""

DIM_MAX: int = 305
"""Maximum DIM to include (standard 305-day lactation)."""

MIN_POINTS_PER_LACTATION: int = 30
"""Fewer than this and we fall back to the parity-pooled curve."""


# ─────────────────────────────────────────────────────────────
#  « Calving data loading »
# ─────────────────────────────────────────────────────────────

def load_calvings(
    data_dir: Path,
    db_path: Path | None = None,
) -> pd.DataFrame:
    """Load one calving-confirmation event per row.

    Resolution order:

    1. ``<data_dir>/calvings.csv`` if present.
    2. Query ``smaxtec_events`` on ``db_path`` and cache to
       ``<data_dir>/calvings.csv``.

    Args:
        data_dir: Directory containing extract-stage CSVs.
        db_path: SQLite DigiMuh database (only used if the CSV
            is missing).

    Returns:
        DataFrame with ``animal_id``, ``calving_date``
        (pandas ``datetime64``).
    """
    csv_path = resolve_input(data_dir, "calvings.csv")
    if csv_path.exists():
        df = pd.read_csv(csv_path, parse_dates=["calving_date"])
        return df

    if db_path is None:
        raise FileNotFoundError(
            f"{csv_path} not found and no db_path given — "
            "re-run `digimuh-extract` or pass db_path explicitly.")

    import sqlite3
    log.info("Querying calving events from %s …", db_path)
    con = sqlite3.connect(str(db_path))
    try:
        df = pd.read_sql(
            "SELECT animal_id, timestamp AS calving_date "
            "FROM smaxtec_events "
            "WHERE event_type = 'calving_confirmation' "
            "ORDER BY animal_id, timestamp",
            con,
            parse_dates=["calving_date"],
        )
    finally:
        con.close()

    # Drop tz; calving day resolution is sufficient.
    df["calving_date"] = df["calving_date"].dt.tz_localize(None).dt.normalize()
    out_path = resolve_output(data_dir, "calvings.csv")
    df.to_csv(out_path, index=False)
    log.info("  wrote %d calving events to %s", len(df), out_path)
    return df


# ─────────────────────────────────────────────────────────────
#  « DIM & lactation-number attachment »
# ─────────────────────────────────────────────────────────────

def attach_dim(
    daily_yield: pd.DataFrame,
    calvings: pd.DataFrame,
    dim_min: int = DIM_MIN,
    dim_max: int = DIM_MAX,
) -> pd.DataFrame:
    """Attach ``dim`` and ``lactation_nr`` to each cow-day.

    For each cow-day, DIM is the number of days between the most
    recent prior calving and the observation date.  Cow-days with
    no prior calving on record, or whose DIM is outside
    ``[dim_min, dim_max]``, are dropped.

    Args:
        daily_yield: DataFrame with ``animal_id``, ``date``,
            ``daily_yield_kg``.
        calvings: Output of :func:`load_calvings`.
        dim_min: Lower DIM bound (inclusive).  Default 5.
        dim_max: Upper DIM bound (inclusive).  Default 305.

    Returns:
        Filtered DataFrame with ``dim`` (int) and ``lactation_nr``
        (1, 2, 3, …) added.  Cow-days are dropped if no prior
        calving exists.
    """
    dy = daily_yield.copy()
    dy["date"] = pd.to_datetime(dy["date"]).astype("datetime64[ns]")

    cal = calvings.copy()
    cal["calving_date"] = (
        pd.to_datetime(cal["calving_date"]).astype("datetime64[ns]"))
    # Assign lactation number first (needs per-animal order), then
    # re-sort globally by the merge key as pd.merge_asof requires.
    cal = cal.sort_values(["animal_id", "calving_date"]).reset_index(drop=True)
    cal["lactation_nr"] = cal.groupby("animal_id").cumcount() + 1
    cal = cal.sort_values("calving_date").reset_index(drop=True)
    dy = dy.sort_values("date").reset_index(drop=True)

    merged = pd.merge_asof(
        dy, cal,
        left_on="date", right_on="calving_date",
        by="animal_id", direction="backward",
    )
    merged = merged.dropna(subset=["calving_date"]).copy()
    merged["dim"] = (merged["date"] - merged["calving_date"]).dt.days.astype(int)
    merged = merged[(merged["dim"] >= dim_min) & (merged["dim"] <= dim_max)]
    merged["lactation_nr"] = merged["lactation_nr"].astype(int)
    return merged.reset_index(drop=True)


# ─────────────────────────────────────────────────────────────
#  « Wood fit »
# ─────────────────────────────────────────────────────────────

def fit_wood(dim: np.ndarray, yld: np.ndarray) -> dict:
    """Fit Wood's incomplete-gamma lactation model via log-linear OLS.

    Model::

        y = a · t^b · exp(-c · t)
        ln y = ln a + b · ln t − c · t

    Args:
        dim: Days in milk (must be > 0).
        yld: Daily yield (kg/d, must be > 0).

    Returns:
        Dict with ``a``, ``b``, ``c``, ``r_squared``, ``n``,
        ``peak_dim``, ``peak_yield`` and ``converged`` (True when
        the fit is biologically plausible: b>0, c>0, peak within
        [5, 120] days).
    """
    dim = np.asarray(dim, dtype=float)
    yld = np.asarray(yld, dtype=float)
    mask = np.isfinite(dim) & np.isfinite(yld) & (dim > 0) & (yld > 0)
    dim, yld = dim[mask], yld[mask]
    n = int(len(dim))

    empty = dict(a=np.nan, b=np.nan, c=np.nan, r_squared=np.nan, n=n,
                 peak_dim=np.nan, peak_yield=np.nan, converged=False)
    if n < 10:
        return empty

    X = np.column_stack([np.ones(n), np.log(dim), dim])
    try:
        coef, *_ = np.linalg.lstsq(X, np.log(yld), rcond=None)
    except np.linalg.LinAlgError:
        return empty

    ln_a, b, neg_c = coef
    a, c = float(np.exp(ln_a)), float(-neg_c)

    # Predicted yield and residuals in original space.
    y_hat_log = X @ coef
    y_hat = np.exp(y_hat_log)
    ss_res = float(np.sum((yld - y_hat) ** 2))
    ss_tot = float(np.sum((yld - yld.mean()) ** 2))
    r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    peak_dim = b / c if c > 0 else np.nan
    peak_yield = (a * peak_dim ** b * np.exp(-c * peak_dim)
                  if np.isfinite(peak_dim) and peak_dim > 0 else np.nan)
    # Biological plausibility AND minimum explanatory power.  R² < 0.1
    # almost always means the cow-year slice missed the curve's
    # informative early rise and the parity pool is a safer prior.
    converged = bool(
        b > 0 and c > 0
        and np.isfinite(peak_dim) and 15 <= peak_dim <= 100
        and np.isfinite(peak_yield) and peak_yield > 0
        and np.isfinite(r_sq) and r_sq >= 0.10
    )
    return dict(a=a, b=float(b), c=c, r_squared=float(r_sq), n=n,
                peak_dim=float(peak_dim), peak_yield=float(peak_yield),
                converged=converged)


def predict_wood(dim, a, b, c) -> np.ndarray:
    """Evaluate the Wood curve at the given DIMs.

    All four arguments broadcast together — scalars or arrays of
    matching length are accepted.  The result is NaN anywhere
    ``dim <= 0`` or any of the parameters is NaN, evaluated
    per-element rather than per-array, so a single NaN in the
    middle of ``a`` / ``b`` / ``c`` does not contaminate the rest
    of the output.
    """
    dim = np.asarray(dim, dtype=float)
    a_arr = np.broadcast_to(np.asarray(a, dtype=float), dim.shape)
    b_arr = np.broadcast_to(np.asarray(b, dtype=float), dim.shape)
    c_arr = np.broadcast_to(np.asarray(c, dtype=float), dim.shape)
    ok = (np.isfinite(dim) & (dim > 0)
          & np.isfinite(a_arr) & np.isfinite(b_arr) & np.isfinite(c_arr))
    out = np.full(dim.shape, np.nan, dtype=float)
    if ok.any():
        safe_dim = np.where(ok, dim, 1.0)
        safe_a   = np.where(ok, a_arr, 0.0)
        safe_b   = np.where(ok, b_arr, 0.0)
        safe_c   = np.where(ok, c_arr, 0.0)
        pred = safe_a * np.power(safe_dim, safe_b) * np.exp(-safe_c * safe_dim)
        out = np.where(ok, pred, np.nan)
    return out


# ─────────────────────────────────────────────────────────────
#  « Per-lactation fits + parity-pooled fallback »
# ─────────────────────────────────────────────────────────────

def _parity_bucket(lactation_nr: int) -> str:
    """Map a 1-based lactation number to ``"1"``, ``"2"``, or ``"3+"``."""
    if lactation_nr <= 1:
        return "1"
    if lactation_nr == 2:
        return "2"
    return "3+"


def fit_wood_per_lactation(
    df_with_dim: pd.DataFrame,
    min_points: int = MIN_POINTS_PER_LACTATION,
) -> pd.DataFrame:
    """Fit one Wood curve per (animal_id, calving_date).

    Lactations with fewer than ``min_points`` valid observations are
    marked ``method = "parity_fallback"`` and will use the
    parity-pooled curve at prediction time.

    Args:
        df_with_dim: Output of :func:`attach_dim`.
        min_points: Per-lactation point threshold for an
            individual fit.  Default 30.

    Returns:
        DataFrame of fit results, one row per (animal_id,
        calving_date), with columns ``a``, ``b``, ``c``,
        ``r_squared``, ``n``, ``peak_dim``, ``peak_yield``,
        ``converged``, ``method``, ``parity``.
    """
    df = df_with_dim.copy()
    df["parity"] = df["lactation_nr"].apply(_parity_bucket)

    # ── Parity-pooled reference curves ─────────────────────
    pool_fits: dict[str, dict] = {}
    for parity, sub in df.groupby("parity"):
        pool_fits[parity] = fit_wood(sub["dim"].values,
                                     sub["daily_yield_kg"].values)

    records = []
    for (aid, cdate), grp in df.groupby(["animal_id", "calving_date"]):
        parity = str(grp["parity"].iloc[0])
        if len(grp) >= min_points:
            fit = fit_wood(grp["dim"].values, grp["daily_yield_kg"].values)
            if fit["converged"]:
                method = "per_lactation"
                rec = {**fit, "method": method}
            else:
                rec = {**pool_fits.get(parity, {}), "method": "parity_fallback"}
        else:
            rec = {**pool_fits.get(parity, {}), "method": "parity_fallback"}
        rec.update({
            "animal_id": int(aid),
            "calving_date": pd.to_datetime(cdate),
            "lactation_nr": int(grp["lactation_nr"].iloc[0]),
            "parity": parity,
            "n_points": int(len(grp)),
        })
        records.append(rec)

    fits = pd.DataFrame(records)
    cols = ["animal_id", "calving_date", "lactation_nr", "parity",
            "method", "n_points", "a", "b", "c",
            "peak_dim", "peak_yield", "r_squared", "converged"]
    return fits.reindex(columns=cols)


# ─────────────────────────────────────────────────────────────
#  « Residual attachment (public entry point) »
# ─────────────────────────────────────────────────────────────

def load_daily_yields_for_fitting(data_dir: Path) -> pd.DataFrame | None:
    """Load the full-history daily yield CSV if present.

    ``daily_milk_yield_full.csv`` is the per-animal full HerdePlus
    history (no date filter) written by
    :func:`digimuh.extract.extract_daily_milk_yield_full`.  It is
    the recommended fit frame for :func:`compute_wood_residuals`:
    per-lactation Wood fits converge far more often when they see
    the full pre-peak / post-peak shape instead of a summer slice.

    Returns:
        DataFrame or ``None`` when the file does not exist.
    """
    path = resolve_input(data_dir, "daily_milk_yield_full.csv")
    if not path.exists():
        return None
    df = pd.read_csv(path)
    df = df[df["daily_yield_kg"].notna() & (df["daily_yield_kg"] > 0)].copy()
    df["date"] = pd.to_datetime(df["date"])
    # We don't need a year column for fitting — attach_dim only uses
    # (animal_id, date) — but some downstream code expects it, so
    # derive it from date when it's absent.
    if "year" not in df.columns:
        df["year"] = df["date"].dt.year.astype(int)
    return df


def compute_wood_residuals(
    daily_yield: pd.DataFrame,
    calvings: pd.DataFrame,
    min_points: int = MIN_POINTS_PER_LACTATION,
    dim_min: int = DIM_MIN,
    dim_max: int = DIM_MAX,
    fit_yields: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """End-to-end: attach DIM, fit per-lactation Wood curves, add residuals.

    The default behaviour (``fit_yields=None``) fits Wood curves on
    the same yield frame used to compute residuals — the original
    single-frame mode.  Passing ``fit_yields`` decouples the two:
    Wood parameters are estimated on ``fit_yields`` (typically the
    full HerdePlus history so each lactation has its pre-peak and
    post-peak tail visible), then residuals are evaluated on
    ``daily_yield`` (typically the analysis window — for the
    Frontiers broken-stick pipeline, the Tierauswahl summer slice).

    Args:
        daily_yield: The cow-days the caller wants residuals for.
            Must have ``animal_id``, ``date``, ``daily_yield_kg``,
            ``year``.
        calvings: From :func:`load_calvings`.
        min_points: Per-lactation point threshold for an individual
            Wood fit (lactations with fewer points use the parity
            pool).
        dim_min, dim_max: DIM clipping range applied to both frames.
        fit_yields: Optional frame to estimate Wood parameters on.
            When given, the returned ``fits`` table reports the
            full-history lactation coverage (typically tens of
            thousands more points per animal), and ``yield_expected``
            on each ``daily_yield`` cow-day comes from the curve
            fitted on that richer data.  When ``None``, falls back
            to the legacy single-frame behaviour.

    Returns:
        ``(yields, fits)`` where ``yields`` is the ``daily_yield``
        cow-day DataFrame enriched with ``dim``, ``lactation_nr``,
        ``parity``, ``yield_expected``, ``yield_residual``,
        ``yield_residual_rel`` (residual / expected), and
        ``method`` (``"per_lactation"`` / ``"parity_fallback"``).
        ``fits`` is the per-lactation fit table.
    """
    # 1. Estimate Wood parameters.
    fit_frame = fit_yields if fit_yields is not None else daily_yield
    fit_dim = attach_dim(fit_frame, calvings,
                         dim_min=dim_min, dim_max=dim_max)
    fits = fit_wood_per_lactation(fit_dim, min_points=min_points)

    # 2. Evaluate residuals on the analysis frame.
    dy = attach_dim(daily_yield, calvings,
                    dim_min=dim_min, dim_max=dim_max)
    dy["parity"] = dy["lactation_nr"].apply(_parity_bucket)
    merged = dy.merge(
        fits[["animal_id", "calving_date", "a", "b", "c", "method"]],
        on=["animal_id", "calving_date"], how="left",
    )
    merged["yield_expected"] = predict_wood(
        merged["dim"].values,
        merged["a"].values, merged["b"].values, merged["c"].values,
    )
    merged["yield_residual"] = merged["daily_yield_kg"] - merged["yield_expected"]
    merged["yield_residual_rel"] = (
        merged["yield_residual"] / merged["yield_expected"].replace(0, np.nan))
    # Keep the stable public columns first; drop the temporary params.
    cols = [
        "animal_id", "date", "year", "daily_yield_kg",
        "dim", "lactation_nr", "parity", "calving_date",
        "yield_expected", "yield_residual", "yield_residual_rel",
        "method",
    ]
    return merged.reindex(columns=cols), fits


# ─────────────────────────────────────────────────────────────
#  « Convenience: residual terciles »
# ─────────────────────────────────────────────────────────────

def residual_terciles(residuals: pd.Series) -> tuple[float, float]:
    """Return the (Q33, Q67) boundaries of a residual series."""
    s = residuals.dropna()
    return float(s.quantile(0.33)), float(s.quantile(0.67))

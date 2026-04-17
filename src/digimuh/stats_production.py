#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_production                                     ║
# ║  « thermoneutral fraction and milk yield coupling »             ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Compute daily thermoneutral fraction (TNF) and correlate with  ║
# ║  P95-normalised daily milk yield.                               ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Production impact analysis for the broken-stick pipeline."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

log = logging.getLogger("digimuh.stats")

def compute_thermoneutral_fraction(
    rumen: pd.DataFrame, bs_results: pd.DataFrame,
) -> pd.DataFrame:
    """Compute daily thermoneutral fraction (TNF) per animal.

    For each animal with a converged THI breakpoint, count the fraction
    of 10-min readings per calendar day where barn THI is at or below
    the individual breakpoint.  This is the fraction of the day the
    animal spends within its thermoneutral zone.

    Also computes analogous fraction for barn temperature breakpoints.

    Args:
        rumen: rumen_barn.csv DataFrame (needs timestamp, barn_thi, barn_temp).
        bs_results: broken_stick_results.csv DataFrame.

    Returns:
        DataFrame with columns: animal_id, year, date, thi_tnf, temp_tnf,
        n_readings, mean_thi, mean_body_temp.
    """
    records = []
    converged = bs_results[bs_results["thi_converged"] == True]

    for _, row in converged.iterrows():
        aid = int(row["animal_id"])
        year = int(row["year"])
        thi_bp = row["thi_breakpoint"]

        # Barn temp breakpoint (may not converge)
        temp_bp = row.get("temp_breakpoint", np.nan)
        temp_conv = row.get("temp_converged", False)

        grp = rumen[(rumen["animal_id"] == aid) & (rumen["year"] == year)].copy()
        if len(grp) < 50:
            continue

        grp["date"] = pd.to_datetime(grp["timestamp"]).dt.date

        for date, day_data in grp.groupby("date"):
            n = len(day_data)
            if n < 6:  # at least 1 hour of data
                continue

            thi_tnf = (day_data["barn_thi"] <= thi_bp).mean()
            temp_tnf = (day_data["barn_temp"] <= temp_bp).mean() if temp_conv else np.nan

            records.append({
                "animal_id": aid,
                "year": year,
                "date": date,
                "thi_tnf": thi_tnf,
                "temp_tnf": temp_tnf,
                "n_readings": n,
                "mean_thi": day_data["barn_thi"].mean(),
                "mean_barn_temp": day_data["barn_temp"].mean(),
                "mean_body_temp": day_data["body_temp"].mean(),
            })

    return pd.DataFrame(records)


def compute_tnf_yield_analysis(
    tnf: pd.DataFrame, daily_yield: pd.DataFrame,
) -> pd.DataFrame:
    """Merge daily TNF with daily milk yield for per-cow, per-day pairs.

    For each cow, each day produces two floats:
    - thi_tnf: fraction of the day below the cow's THI breakpoint (0-1)
    - relative_yield: that day's milk yield / cow-specific P95 (0-1ish)

    The P95 is the 95th percentile of each individual cow's daily yields
    across her entire dataset (all years).  This is a robust, cow-specific
    reference maximum that avoids outlier sensitivity of the absolute max.

    Args:
        tnf: Daily TNF DataFrame from compute_thermoneutral_fraction.
        daily_yield: Daily milk yield DataFrame (animal_id, date,
            daily_yield_kg, year).

    Returns:
        DataFrame with one row per cow-day: animal_id, year, date,
        thi_tnf, temp_tnf, daily_yield_kg, yield_p95, relative_yield.
    """
    if tnf.empty or daily_yield.empty:
        return pd.DataFrame()

    # Ensure date columns are comparable
    tnf = tnf.copy()
    daily_yield = daily_yield.copy()
    tnf["date"] = pd.to_datetime(tnf["date"]).dt.date
    daily_yield["date"] = pd.to_datetime(daily_yield["date"]).dt.date

    # Merge TNF and yield on (animal_id, date)
    merged = tnf.merge(
        daily_yield[["animal_id", "date", "daily_yield_kg", "n_milkings"]],
        on=["animal_id", "date"],
        how="inner",
    )

    if merged.empty:
        return pd.DataFrame()

    # Cow-specific P95: robust reference maximum from ALL that cow's daily yields
    cow_p95 = daily_yield.groupby("animal_id")["daily_yield_kg"].agg(
        yield_p95=lambda x: np.nanpercentile(x.dropna(), 95),
        yield_median=lambda x: x.median(),
        n_yield_days=lambda x: x.notna().sum(),
    ).reset_index()

    merged = merged.merge(cow_p95, on="animal_id", how="left")

    # Relative yield: this day's yield / this cow's P95
    merged["relative_yield"] = merged["daily_yield_kg"] / merged["yield_p95"]

    return merged


# ─────────────────────────────────────────────────────────────
#  « stratified TNF × yield analysis by Wood-residual class »
# ─────────────────────────────────────────────────────────────

def classify_cow_years_by_wood_residual(
    wood_residuals: pd.DataFrame,
    terciles: tuple[float, float] | None = None,
    min_days: int = 10,
) -> tuple[pd.DataFrame, tuple[float, float]]:
    """Assign each (animal_id, year) a low / middle / high class.

    The class is based on the cow-year's mean Wood residual so that
    a single label follows each animal-lactation across all her
    cow-days.  Tercile boundaries default to Q33/Q67 of the
    cow-year means themselves (not the cow-day level) so the three
    classes each hold roughly a third of cow-years.

    Args:
        wood_residuals: Output of ``compute_wood_residuals`` — one
            row per cow-day with ``animal_id``, ``year``,
            ``yield_residual``.
        terciles: Optional (Q33, Q67) boundary override.  If
            ``None`` the boundaries are computed from the cow-year
            means.
        min_days: Require this many cow-days per (animal_id, year)
            for the mean to be trusted; cow-years below this are
            dropped.

    Returns:
        ``(class_table, (q33, q67))`` — class_table has columns
        ``animal_id``, ``year``, ``mean_residual``, ``n_days``,
        ``yield_class`` (categorical, ordered low < middle < high).
    """
    wr = wood_residuals.dropna(subset=["yield_residual"])
    if wr.empty:
        return pd.DataFrame(), (np.nan, np.nan)

    agg = wr.groupby(["animal_id", "year"]).agg(
        mean_residual=("yield_residual", "mean"),
        n_days=("yield_residual", "size"),
    ).reset_index()
    agg = agg[agg["n_days"] >= min_days].copy()
    if agg.empty:
        return agg.assign(yield_class=[]), (np.nan, np.nan)

    if terciles is None:
        q33 = float(agg["mean_residual"].quantile(0.33))
        q67 = float(agg["mean_residual"].quantile(0.67))
    else:
        q33, q67 = terciles

    labels = np.where(
        agg["mean_residual"] <= q33, "low",
        np.where(agg["mean_residual"] <= q67, "middle", "high"),
    )
    agg["yield_class"] = pd.Categorical(
        labels, categories=["low", "middle", "high"], ordered=True)
    return agg, (q33, q67)


def compute_tnf_yield_by_class(
    tnf_yield: pd.DataFrame,
    class_table: pd.DataFrame,
    wood: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Attach ``yield_class`` (and optional DIM-adjusted residual) to a
    TNF × yield table.

    Args:
        tnf_yield: ``tnf_yield.csv`` DataFrame (one row per cow-day
            with ``thi_tnf``, ``temp_tnf``, ``daily_yield_kg``,
            ``relative_yield``, ``year``).
        class_table: Output of
            :func:`classify_cow_years_by_wood_residual`.
        wood: Optional ``daily_milk_yield_wood.csv`` DataFrame.  When
            given, ``yield_residual``, ``yield_expected`` and ``dim``
            are merged onto each cow-day on ``(animal_id, date)`` so
            downstream correlations can use the DIM-adjusted residual
            as the response and remove the within-cow-year lactation
            decline from the heat-stress signal.

    Returns:
        Filtered cow-day DataFrame restricted to cow-years that
        have a class; gains ``yield_class``, ``mean_residual`` and
        (when ``wood`` is supplied) ``yield_residual``,
        ``yield_expected``, ``dim``.
    """
    if tnf_yield.empty or class_table.empty:
        return pd.DataFrame()
    merged = tnf_yield.merge(
        class_table[["animal_id", "year", "yield_class", "mean_residual"]],
        on=["animal_id", "year"], how="inner",
    )
    if wood is not None and not wood.empty:
        w = wood[["animal_id", "date", "yield_residual",
                  "yield_expected", "dim"]].copy()
        w["date"] = pd.to_datetime(w["date"]).dt.date
        merged["date"] = pd.to_datetime(merged["date"]).dt.date
        merged = merged.merge(w, on=["animal_id", "date"], how="left")
    return merged


def tnf_yield_correlations_by_class(
    tnf_by_class: pd.DataFrame,
    predictors: tuple[str, ...] = ("thi_tnf", "temp_tnf"),
    responses: tuple[str, ...] = (
        "daily_yield_kg", "relative_yield", "yield_residual"),
) -> pd.DataFrame:
    """Per-class Spearman correlations of TNF vs yield.

    Args:
        tnf_by_class: Output of :func:`compute_tnf_yield_by_class`.
        predictors: TNF columns to test.
        responses: Yield columns to test.

    Returns:
        DataFrame of ``yield_class`` × ``predictor`` × ``response``
        with ``n``, ``n_animals``, ``rs``, ``p``, ``slope``,
        ``intercept``, ``median_y``.  The class ``"pooled"`` is
        included as a reference row.
    """
    rows: list[dict] = []
    classes = ["low", "middle", "high"]
    groups = [(c, tnf_by_class[tnf_by_class["yield_class"] == c])
              for c in classes]
    groups.append(("pooled", tnf_by_class))

    for cls, sub in groups:
        for pred in predictors:
            for resp in responses:
                if pred not in sub.columns or resp not in sub.columns:
                    continue
                valid = sub.dropna(subset=[pred, resp])
                n = len(valid)
                if n < 10:
                    rows.append(dict(
                        yield_class=cls, predictor=pred, response=resp,
                        n=n, n_animals=valid["animal_id"].nunique(),
                        rs=np.nan, p=np.nan,
                        slope=np.nan, intercept=np.nan, median_y=np.nan,
                    ))
                    continue
                rs, p = spearmanr(valid[pred], valid[resp])
                slope, intercept = np.polyfit(
                    valid[pred].astype(float),
                    valid[resp].astype(float), 1)
                rows.append(dict(
                    yield_class=cls, predictor=pred, response=resp,
                    n=n, n_animals=valid["animal_id"].nunique(),
                    rs=float(rs), p=float(p),
                    slope=float(slope), intercept=float(intercept),
                    median_y=float(valid[resp].median()),
                ))
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────
#  « crossing-day flag + daily climate means »
# ─────────────────────────────────────────────────────────────

def compute_daily_crossing_flags(crossing_times: pd.DataFrame) -> pd.DataFrame:
    """Per (animal_id, year, date): did any THI / barn-temp crossing occur?

    Consumes ``crossing_times.csv`` — one row per upward crossing of
    the cow's individual breakpoint — and collapses to a cow-day
    table with two boolean flags so downstream plots can split
    cow-days into "days with a crossing" vs "days without".

    Args:
        crossing_times: ``crossing_times.csv`` DataFrame with
            ``animal_id``, ``year``, ``date``, ``predictor``
            (``"thi"`` / ``"temp"``).

    Returns:
        DataFrame with ``animal_id``, ``year``, ``date``,
        ``thi_crossed``, ``temp_crossed`` (booleans),
        ``n_thi_crossings``, ``n_temp_crossings`` (ints).
    """
    if crossing_times.empty:
        return pd.DataFrame(columns=["animal_id", "year", "date",
                                     "thi_crossed", "temp_crossed",
                                     "n_thi_crossings", "n_temp_crossings"])
    ct = crossing_times.copy()
    ct["date"] = pd.to_datetime(ct["date"]).dt.date
    counts = (
        ct.groupby(["animal_id", "year", "date", "predictor"])
          .size().unstack(fill_value=0).reset_index()
    )
    counts.columns.name = None
    counts = counts.rename(columns={
        "thi": "n_thi_crossings",
        "temp": "n_temp_crossings",
    })
    for col in ("n_thi_crossings", "n_temp_crossings"):
        if col not in counts.columns:
            counts[col] = 0
    counts["thi_crossed"]  = counts["n_thi_crossings"]  > 0
    counts["temp_crossed"] = counts["n_temp_crossings"] > 0
    return counts[["animal_id", "year", "date",
                   "thi_crossed", "temp_crossed",
                   "n_thi_crossings", "n_temp_crossings"]]


def attach_daily_climate_means(
    tnf_by_class: pd.DataFrame,
    rumen: pd.DataFrame,
) -> pd.DataFrame:
    """Ensure cow-day rows carry ``mean_thi`` and ``mean_barn_temp``.

    ``mean_thi`` is usually already present from
    ``compute_thermoneutral_fraction``; ``mean_barn_temp`` was
    added later so legacy CSVs may lack it.  Missing columns are
    re-derived from ``rumen_barn.csv`` via a per-cow-day groupby.
    """
    if tnf_by_class.empty:
        return tnf_by_class
    if "mean_thi" in tnf_by_class.columns and "mean_barn_temp" in tnf_by_class.columns:
        return tnf_by_class
    if rumen.empty:
        return tnf_by_class
    r = rumen.copy()
    r["timestamp"] = pd.to_datetime(r["timestamp"])
    r["date"] = r["timestamp"].dt.date
    means = r.groupby(["animal_id", "year", "date"]).agg(
        mean_thi=("barn_thi", "mean"),
        mean_barn_temp=("barn_temp", "mean"),
    ).reset_index()
    merged = tnf_by_class.copy()
    merged["date"] = pd.to_datetime(merged["date"]).dt.date
    # Drop then re-merge so new values win over legacy ones.
    drop_cols = [c for c in ("mean_thi", "mean_barn_temp")
                 if c in merged.columns]
    if drop_cols:
        merged = merged.drop(columns=drop_cols)
    return merged.merge(means, on=["animal_id", "year", "date"], how="left")


def attach_crossing_flags(
    tnf_by_class: pd.DataFrame,
    crossing_flags: pd.DataFrame,
) -> pd.DataFrame:
    """Merge the cow-day crossing flag table onto a TNF × class table."""
    if tnf_by_class.empty or crossing_flags.empty:
        return tnf_by_class
    merged = tnf_by_class.copy()
    merged["date"] = pd.to_datetime(merged["date"]).dt.date
    cf = crossing_flags.copy()
    cf["date"] = pd.to_datetime(cf["date"]).dt.date
    return merged.merge(
        cf, on=["animal_id", "year", "date"], how="left",
    ).assign(
        thi_crossed=lambda d: d["thi_crossed"].fillna(False),
        temp_crossed=lambda d: d["temp_crossed"].fillna(False),
        n_thi_crossings=lambda d: d["n_thi_crossings"].fillna(0).astype(int),
        n_temp_crossings=lambda d: d["n_temp_crossings"].fillna(0).astype(int),
    )


def crossing_day_comparison(
    df: pd.DataFrame,
    response: str = "yield_residual",
) -> pd.DataFrame:
    """Per-predictor effect of "crossed that day" on a yield response.

    For each climate predictor (THI, barn-temp) and for the pooled
    cow-days plus each yield_class (if present), reports:

    * group sizes,
    * medians,
    * Mann-Whitney U statistic and p,
    * median-difference effect size (crossed − not-crossed).

    Args:
        df: Cow-day DataFrame that has ``thi_crossed``,
            ``temp_crossed``, the response column, and optionally
            ``yield_class``.
        response: Response column name.

    Returns:
        Long DataFrame, one row per (group × predictor).
    """
    from scipy.stats import mannwhitneyu
    rows = []
    groups = [("pooled", df)]
    if "yield_class" in df.columns:
        for cls in ("low", "middle", "high"):
            groups.append((cls, df[df["yield_class"] == cls]))

    for grp_name, sub in groups:
        for pred, flag in [("thi", "thi_crossed"),
                           ("temp", "temp_crossed")]:
            if flag not in sub.columns or response not in sub.columns:
                continue
            valid = sub.dropna(subset=[response, flag])
            y_yes = valid.loc[valid[flag] == True, response].values
            y_no  = valid.loc[valid[flag] == False, response].values
            if len(y_yes) < 5 or len(y_no) < 5:
                rows.append(dict(
                    group=grp_name, predictor=pred,
                    n_yes=len(y_yes), n_no=len(y_no),
                    median_yes=np.nan, median_no=np.nan,
                    median_diff=np.nan, U=np.nan, p=np.nan,
                ))
                continue
            U, p = mannwhitneyu(y_yes, y_no, alternative="two-sided")
            rows.append(dict(
                group=grp_name, predictor=pred,
                n_yes=int(len(y_yes)), n_no=int(len(y_no)),
                median_yes=float(np.median(y_yes)),
                median_no=float(np.median(y_no)),
                median_diff=float(np.median(y_yes) - np.median(y_no)),
                U=float(U), p=float(p),
            ))
    return pd.DataFrame(rows)


def daily_climate_vs_yield_correlations(
    df: pd.DataFrame,
    response: str = "yield_residual",
) -> pd.DataFrame:
    """Spearman rs of daily mean THI / barn-temp vs a yield response."""
    from scipy.stats import spearmanr
    rows = []
    for pred in ("mean_thi", "mean_barn_temp"):
        if pred not in df.columns or response not in df.columns:
            continue
        valid = df.dropna(subset=[pred, response])
        n = len(valid)
        if n < 10:
            rows.append(dict(predictor=pred, n=n, n_animals=0,
                             rs=np.nan, p=np.nan,
                             slope=np.nan, intercept=np.nan))
            continue
        rs, p = spearmanr(valid[pred], valid[response])
        slope, intercept = np.polyfit(
            valid[pred].astype(float),
            valid[response].astype(float), 1,
        )
        rows.append(dict(
            predictor=pred, n=n,
            n_animals=int(valid["animal_id"].nunique()),
            rs=float(rs), p=float(p),
            slope=float(slope), intercept=float(intercept),
        ))
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────

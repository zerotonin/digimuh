#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_milk_composition                               ║
# ║  « MLP composition × climate (thin-milk hypothesis) »          ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Joins MLP (Milchleistungsprüfung) monthly test-day rows with   ║
# ║  cow-day climate and yield-class info.  Tests the "thin milk"   ║
# ║  hypothesis: on hot days cows produce *more* milk (volume up)   ║
# ║  but with *lower* solids (fat %, protein % down) because rumen  ║
# ║  activity is suppressed while water output is increased for     ║
# ║  calf hydration / evaporative cooling.                          ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Milk composition vs climate analysis on MLP test-days.

The MLP channel (§3.20) supplies once-a-month per-cow test-day rows
with fat %, protein %, lactose, fat-kg, ECM, F/E ratio, SCC and
urea.  Matching these to the cow-day climate from ``tnf_yield.csv``
lets us ask whether heat exposure shifts milk *composition* in the
direction the "thin milk" hypothesis predicts:

=====================  =======================  =====================
Hypothesis prediction  Expected sign vs hot day  Column name
=====================  =======================  =====================
more milk volume       positive                 ``herdeplus_mlp_mkg``
less fat per kg        negative                 ``herdeplus_mlp_fat_percent``
less protein per kg    negative                 ``herdeplus_mlp_protein_percent``
total fat roughly flat ≈ zero                   ``herdeplus_mlp_fkg``
F/E ratio drops        negative                 ``herdeplus_mlp_f_e``
ECM drops              negative                 ``herdeplus_mlp_ecm``
=====================  =======================  =====================

This module exposes reusable helpers so the milk-yield classifier
can orchestrate the join, correlation and plot without
re-implementing the logic.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from digimuh.paths import resolve_input

log = logging.getLogger("digimuh.milk_comp")


# Column subset we care about in the MLP table (ordered for display).
# Notes on the two misleadingly-named HerdePlus columns:
#   * ``herdeplus_mlp_ekg_percent`` is actually *protein yield* in kg/d
#     (Eiweißkilogramm) despite the ``_percent`` suffix — confirmed by
#     ekg_percent/mkg × 100 ≈ protein_percent in the data.
#   * ``herdeplus_mlp_lkg`` is *lactose* yield in kg/d (Laktose-kg)
#     — not protein-kg as one might guess from the letter.
MLP_RESPONSES: tuple[tuple[str, str, str], ...] = (
    # (column, short-label, unit)
    ("herdeplus_mlp_mkg",             "Test-day milk",     "kg/d"),
    ("herdeplus_mlp_fat_percent",     "Fat %",             "%"),
    ("herdeplus_mlp_protein_percent", "Protein %",         "%"),
    ("herdeplus_mlp_lactose",         "Lactose %",         "%"),
    ("herdeplus_mlp_fkg",             "Fat",               "kg/d"),
    ("herdeplus_mlp_ekg_percent",     "Protein (Eiweiß)",  "kg/d"),
    ("herdeplus_mlp_lkg",             "Lactose",           "kg/d"),
    ("herdeplus_mlp_f_e",             "F/E ratio",         "—"),
    ("herdeplus_mlp_ecm",             "ECM (energy-corrected milk)", "kg/d"),
    ("herdeplus_mlp_cell_count",      "Somatic cell count", "10³/mL"),
    ("herdeplus_mlp_urea",            "Urea",              "mg/dL"),
)

CLIMATE_PREDICTORS: tuple[tuple[str, str], ...] = (
    # (column, label)
    ("mean_thi",       "Daily mean barn THI"),
    ("mean_barn_temp", "Daily mean barn temperature (°C)"),
    ("thi_tnf",        "Daily time below THI breakpoint (TNF)"),
)


# ─────────────────────────────────────────────────────────────
#  « data loading / joining »
# ─────────────────────────────────────────────────────────────

def load_mlp_test_days(data_dir: Path) -> pd.DataFrame:
    """Load ``mlp_test_days.csv`` and add a day-resolution ``date`` column."""
    path = resolve_input(data_dir, "mlp_test_days.csv")
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found — re-run digimuh-extract to produce it.")
    df = pd.read_csv(path)
    df["timestamp"] = pd.to_datetime(df["timestamp"])
    df["date"] = df["timestamp"].dt.date
    return df


def merge_mlp_with_cowday(
    mlp: pd.DataFrame,
    wood: pd.DataFrame,
    tnf_yield: pd.DataFrame,
    class_table: pd.DataFrame | None = None,
    rumen: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Join an MLP test-day table to the cow-day climate + residual frames.

    The ``(animal_id, date)`` join is exact: MLP timestamps land on
    the test-day itself.  Rows with no matching cow-day (e.g. MLP
    test-days outside the Tierauswahl summer window) are dropped —
    there is no climate to correlate against for those.

    Args:
        mlp: ``mlp_test_days.csv`` DataFrame.
        wood: ``daily_milk_yield_wood.csv`` DataFrame.  Provides
            ``yield_residual``, ``dim``, ``year``.
        tnf_yield: ``tnf_yield.csv`` DataFrame.  Provides
            ``thi_tnf``, ``temp_tnf``, ``mean_thi`` (and
            ``mean_barn_temp`` when available).
        class_table: Optional ``yield_class_per_cow_year.csv`` frame
            (from :func:`stats_production.classify_cow_years_by_wood_residual`);
            when supplied, the output carries a ``yield_class``
            column.

    Returns:
        DataFrame with one row per matched MLP test-day and all
        MLP, climate, residual and (optionally) class columns.
    """
    if mlp.empty or wood.empty or tnf_yield.empty:
        return pd.DataFrame()

    mlp = mlp.copy()
    mlp["date"] = pd.to_datetime(mlp["date"]).dt.date

    w = wood[["animal_id", "date", "year", "yield_residual",
              "dim", "lactation_nr"]].copy()
    w["date"] = pd.to_datetime(w["date"]).dt.date

    t_cols = [c for c in (
        "thi_tnf", "temp_tnf", "mean_thi", "mean_barn_temp", "mean_body_temp"
    ) if c in tnf_yield.columns]
    t = tnf_yield[["animal_id", "date"] + t_cols].copy()
    t["date"] = pd.to_datetime(t["date"]).dt.date

    merged = mlp.merge(w, on=["animal_id", "date"], how="inner")
    merged = merged.merge(t, on=["animal_id", "date"], how="left")

    # Back-fill mean_barn_temp when the loaded tnf_yield.csv pre-dates
    # the column (we added it in §3.17 but older CSVs are missing it).
    if "mean_barn_temp" not in merged.columns and rumen is not None \
            and not rumen.empty:
        from digimuh.stats_production import attach_daily_climate_means
        merged = attach_daily_climate_means(merged, rumen)

    if class_table is not None and not class_table.empty:
        merged = merged.merge(
            class_table[["animal_id", "year", "yield_class", "mean_residual"]],
            on=["animal_id", "year"], how="left",
        )
    return merged


# ─────────────────────────────────────────────────────────────
#  « correlations »
# ─────────────────────────────────────────────────────────────

def mlp_climate_correlations(
    df: pd.DataFrame,
    responses: tuple[tuple[str, str, str], ...] = MLP_RESPONSES,
    predictors: tuple[tuple[str, str], ...] = CLIMATE_PREDICTORS,
) -> pd.DataFrame:
    """Spearman rs + p + OLS slope per (group × predictor × response).

    ``group`` is ``"pooled"`` by default; when the input has a
    ``yield_class`` column each of its categories is added as a
    separate group row.

    Returns:
        Long DataFrame with columns ``group``, ``predictor``,
        ``response``, ``response_label``, ``unit``, ``n``,
        ``n_animals``, ``rs``, ``p``, ``slope``, ``intercept``.
    """
    if df.empty:
        return pd.DataFrame()

    groups = [("pooled", df)]
    if "yield_class" in df.columns:
        for cls in ("low", "middle", "high"):
            sub = df[df["yield_class"] == cls]
            if not sub.empty:
                groups.append((cls, sub))

    rows: list[dict] = []
    for grp_name, sub in groups:
        for pred, _pred_label in predictors:
            if pred not in sub.columns:
                continue
            for resp, resp_label, unit in responses:
                if resp not in sub.columns:
                    continue
                valid = sub.dropna(subset=[pred, resp])
                n = len(valid)
                if n < 10:
                    rows.append(dict(
                        group=grp_name, predictor=pred, response=resp,
                        response_label=resp_label, unit=unit,
                        n=n, n_animals=valid["animal_id"].nunique(),
                        rs=np.nan, p=np.nan,
                        slope=np.nan, intercept=np.nan,
                    ))
                    continue
                rs, p = spearmanr(valid[pred], valid[resp])
                slope, intercept = np.polyfit(
                    valid[pred].astype(float),
                    valid[resp].astype(float), 1,
                )
                rows.append(dict(
                    group=grp_name, predictor=pred, response=resp,
                    response_label=resp_label, unit=unit,
                    n=n, n_animals=int(valid["animal_id"].nunique()),
                    rs=float(rs), p=float(p),
                    slope=float(slope), intercept=float(intercept),
                ))
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────
#  « dilution partition (pure water addition vs rumen suppression) »
# ─────────────────────────────────────────────────────────────

def compute_dilution_partition(merged: pd.DataFrame) -> pd.DataFrame:
    """Partition each composition response into dilution + rumen effects.

    For every cow-day row we compute what the fat % and protein %
    *would be* if the cow produced her personal reference amount of
    fat and protein (in kg) but diluted them into the actually-
    observed volume.  Deviations from that dilution prediction tell
    us how much of the observed composition drop on hot days is
    explained by added water alone versus how much extra
    suppression comes from reduced rumen output.

    Per-cow references are the animal's own mean test-day
    ``fat_kg`` / ``protein_kg`` / ``volume`` across all of her MLP
    rows in the merged frame — so "dilution-predicted" is
    self-referential to the cow and does not borrow composition
    across animals.

    Added columns:
      - ``volume_ratio``         = observed_volume / ref_volume
      - ``fat_percent_diluted``  = fat_kg_ref / observed_volume × 100
      - ``protein_percent_diluted`` = protein_kg_ref / observed_volume × 100
      - ``fat_percent_rumen``    = observed − dilution-predicted fat %
      - ``protein_percent_rumen``= observed − dilution-predicted protein %
      - ``fat_kg_ref``, ``protein_kg_ref``, ``volume_ref``
        (the per-cow references, repeated on every row)
    """
    if merged.empty:
        return merged
    df = merged.copy()
    needed = {"herdeplus_mlp_mkg", "herdeplus_mlp_fkg",
              "herdeplus_mlp_ekg_percent",
              "herdeplus_mlp_lkg",
              "herdeplus_mlp_fat_percent",
              "herdeplus_mlp_protein_percent",
              "herdeplus_mlp_lactose"}
    missing = needed - set(df.columns)
    if missing:
        log.info("  dilution partition skipped — missing columns: %s",
                 missing)
        return df

    # Per-cow reference: mean over all her MLP test-days in the merged frame.
    # Lactose reference comes from the lkg column (lactose-kg/day) despite
    # its terse name; it's the only "sugar" channel in the MLP table and
    # accounts for ~99% of milk carbohydrate.
    ref = (
        df.groupby("animal_id")
          .agg(volume_ref=("herdeplus_mlp_mkg", "mean"),
               fat_kg_ref=("herdeplus_mlp_fkg", "mean"),
               protein_kg_ref=("herdeplus_mlp_ekg_percent", "mean"),
               lactose_kg_ref=("herdeplus_mlp_lkg", "mean"))
          .reset_index()
    )
    df = df.merge(ref, on="animal_id", how="left")

    v_obs = df["herdeplus_mlp_mkg"].astype(float)
    df["volume_ratio"] = v_obs / df["volume_ref"]

    # Pure-dilution predictions: keep the cow's reference absolute
    # output (fat_kg, protein_kg, lactose_kg) and divide by the
    # observed volume.
    safe_v = v_obs.replace(0, np.nan)
    df["fat_percent_diluted"] = df["fat_kg_ref"] / safe_v * 100.0
    df["protein_percent_diluted"] = df["protein_kg_ref"] / safe_v * 100.0
    df["lactose_percent_diluted"] = df["lactose_kg_ref"] / safe_v * 100.0

    # Rumen suppression component = what the observed composition
    # undershoots (or overshoots) the dilution prediction by.  For
    # lactose a *positive* residual means the cow is actively
    # synthesising more lactose than her baseline — the osmotic
    # reason volume rose in the first place.
    df["fat_percent_rumen"] = (
        df["herdeplus_mlp_fat_percent"] - df["fat_percent_diluted"])
    df["protein_percent_rumen"] = (
        df["herdeplus_mlp_protein_percent"] - df["protein_percent_diluted"])
    df["lactose_percent_rumen"] = (
        df["herdeplus_mlp_lactose"] - df["lactose_percent_diluted"])
    return df


def dilution_partition_summary(df: pd.DataFrame,
                               predictor: str = "mean_thi") -> pd.DataFrame:
    """Compare observed vs dilution-predicted composition slopes.

    For each of fat % and protein %, computes three correlations
    (rs + OLS slope) against the climate predictor:

      1. ``observed``           — the actual composition drop
      2. ``dilution-predicted`` — the drop implied by the cow's own
                                  reference fat/protein kg diluted
                                  into the observed volume
      3. ``rumen residual``     — observed − dilution-predicted;
                                  positive means rumen over-produced
                                  relative to the reference, negative
                                  means rumen was suppressed below
                                  the reference

    A cleanly *dilutive* explanation has a near-zero rumen-residual
    slope and a dilution-predicted slope that matches the observed
    one.  A *suppression-dominant* story has an observed slope
    noticeably steeper than the dilution-only one.

    Returns a long DataFrame with columns ``nutrient`` (``"fat"`` or
    ``"protein"``), ``component`` (``"observed"`` /
    ``"dilution_predicted"`` / ``"rumen_residual"``), ``rs``, ``p``,
    ``slope``, ``n``.
    """
    from scipy.stats import spearmanr
    rows: list[dict] = []
    if df.empty or predictor not in df.columns:
        return pd.DataFrame(rows)

    specs = [
        ("fat", "herdeplus_mlp_fat_percent",
         "fat_percent_diluted", "fat_percent_rumen"),
        ("protein", "herdeplus_mlp_protein_percent",
         "protein_percent_diluted", "protein_percent_rumen"),
        ("lactose", "herdeplus_mlp_lactose",
         "lactose_percent_diluted", "lactose_percent_rumen"),
    ]
    for nutrient, obs_col, dil_col, rum_col in specs:
        for comp_label, col in [
            ("observed",            obs_col),
            ("dilution_predicted",  dil_col),
            ("rumen_residual",      rum_col),
        ]:
            if col not in df.columns:
                continue
            valid = df.dropna(subset=[predictor, col])
            n = len(valid)
            if n < 10:
                continue
            rs, p = spearmanr(valid[predictor], valid[col])
            slope, intercept = np.polyfit(
                valid[predictor].astype(float),
                valid[col].astype(float), 1)
            rows.append(dict(
                nutrient=nutrient, component=comp_label,
                n=int(n), rs=float(rs), p=float(p),
                slope=float(slope), intercept=float(intercept),
            ))
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────
#  « thin-milk hypothesis summary »
# ─────────────────────────────────────────────────────────────

def thin_milk_verdict(correlations: pd.DataFrame,
                      predictor: str = "mean_thi") -> dict:
    """Distil the MLP × climate table into a one-line verdict.

    The "thin milk" hypothesis predicts, under a hotter cow-day:

    * milk volume (``mkg``) rises with climate  — positive rs
    * fat % and protein % drop                 — negative rs
    * F/E ratio drops                          — negative rs

    Returns a dict with the pooled rs for each component plus a
    single ``"verdict"`` string (``"supported"`` / ``"partial"`` /
    ``"refuted"``).
    """
    if correlations.empty:
        return {"verdict": "no-data"}

    pool = correlations[(correlations["group"] == "pooled")
                        & (correlations["predictor"] == predictor)]

    def rs_of(col):
        row = pool[pool["response"] == col]
        return float(row["rs"].iloc[0]) if not row.empty else np.nan

    rs_mkg = rs_of("herdeplus_mlp_mkg")
    rs_fat = rs_of("herdeplus_mlp_fat_percent")
    rs_prt = rs_of("herdeplus_mlp_protein_percent")
    rs_fe  = rs_of("herdeplus_mlp_f_e")

    supports_volume = np.isfinite(rs_mkg) and rs_mkg > 0
    supports_solids = (np.isfinite(rs_fat) and rs_fat < 0
                       and np.isfinite(rs_prt) and rs_prt < 0)
    n_support = int(supports_volume) + int(supports_solids)
    verdict = ("supported" if n_support == 2
               else "partial" if n_support == 1
               else "refuted")
    return dict(
        predictor=predictor,
        rs_milk_volume=rs_mkg,
        rs_fat_percent=rs_fat,
        rs_protein_percent=rs_prt,
        rs_fat_protein_ratio=rs_fe,
        volume_component=supports_volume,
        solids_component=supports_solids,
        verdict=verdict,
    )

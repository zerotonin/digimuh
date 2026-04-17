#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  Analysis 0a — Data extraction ¹                             ║
# ║  « DB → CSV for broken-stick analysis »                      ║
# ╚══════════════════════════════════════════════════════════════╝
"""Extract rumen temperature, respiration, barn climate, and production
data from the DigiMuh database into analysis-ready CSVs.

Run this once.  Downstream scripts (``stats_runner`` and
``viz_runner``) read the CSVs and never touch the database.

Outputs (in ``--out`` directory)::

    tierauswahl.csv          Animal selection (cleaned from xlsx)
    rumen_barn.csv           Rumen temp + barn THI/temp per 10-min tick
    respiration_barn.csv     Respiration + barn THI/temp per reading
    production.csv           Mean milk yield + lactation nr per animal
    climate_daily.csv        Daily barn climate summary (Jun–Sep)

Usage::

    digimuh-extract --db cow.db --tierauswahl Tierauswahl.xlsx \\
        --out results/broken_stick

¹ Analysis led by Dr. med. vet. Gundula Hoffmann, ATB Potsdam.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from digimuh.analysis_utils import connect_db, query_df
from digimuh.paths import resolve_output

log = logging.getLogger("digimuh.extract")

# ─────────────────────────────────────────────────────────────
#  « constants »
# ─────────────────────────────────────────────────────────────

NEUBAU_BARN_IDS = (1, 2)  # NewBridge, NewPillar
DRINK_PAD = pd.Timedelta(minutes=15)

# ─────────────────────────────────────────────────────────────
#  « SQL queries »
# ─────────────────────────────────────────────────────────────

SQL_RUMEN = """
SELECT
    animal_id,
    "timestamp",
    CAST("temp_without_drink_cycles" AS REAL) AS body_temp,
    CAST("drink_cycles_v2" AS REAL)           AS drink_cycles
FROM smaxtec_derived
WHERE animal_id = ?
  AND "timestamp" >= ? AND "timestamp" <= ?
  AND "temp_without_drink_cycles" IS NOT NULL
  AND CAST("temp_without_drink_cycles" AS REAL) BETWEEN 30 AND 43
"""

SQL_BARN = """
SELECT "timestamp",
       AVG(temp)           AS barn_temp,
       AVG(temp_hum_index) AS barn_thi
FROM smaxtec_barns
WHERE "timestamp" >= ? AND "timestamp" <= ?
  AND temp IS NOT NULL AND temp_hum_index IS NOT NULL
  AND barn_id IN (1, 2)
GROUP BY "timestamp"
"""

SQL_RESPIRATION = """
SELECT animal_id, "timestamp",
       respirationfrequency AS resp_rate
FROM gouna
WHERE animal_id = ?
  AND "timestamp" >= ? AND "timestamp" <= ?
  AND respirationfrequency IS NOT NULL
  AND respirationfrequency BETWEEN 5 AND 150
"""

SQL_MILK_YIELD = """
SELECT animal_id,
       AVG(herdeplus_milked_mkg) AS mean_milk_yield_kg,
       COUNT(*)                  AS n_milkings
FROM herdeplus
WHERE animal_id = ?
  AND "timestamp" >= ? AND "timestamp" <= ?
  AND herdeplus_milked_mkg IS NOT NULL
  AND herdeplus_milked_mkg > 0
GROUP BY animal_id
"""

SQL_DAILY_MILK_YIELD = """
SELECT animal_id,
       DATE("timestamp") AS date,
       SUM(herdeplus_milked_mkg) AS daily_yield_kg,
       COUNT(*)                   AS n_milkings
FROM herdeplus
WHERE animal_id = ?
  AND "timestamp" >= ? AND "timestamp" <= ?
  AND herdeplus_milked_mkg IS NOT NULL
  AND herdeplus_milked_mkg > 0
GROUP BY animal_id, DATE("timestamp")
"""

SQL_LACTATION_NR = """
SELECT MAX(CAST(herdeplus_calving_lactation AS INTEGER)) AS lactation_nr
FROM herdeplus
WHERE animal_id = ?
  AND herdeplus_calving_lactation IS NOT NULL
  AND herdeplus_calving_lactation > 0
"""

SQL_CALVINGS = """
SELECT animal_id,
       date("timestamp") AS calving_date
FROM smaxtec_events
WHERE event_type = 'calving_confirmation'
ORDER BY animal_id, calving_date
"""

SQL_DAILY_MILK_YIELD_FULL = """
SELECT animal_id,
       DATE("timestamp") AS date,
       SUM(herdeplus_milked_mkg) AS daily_yield_kg,
       COUNT(*)                   AS n_milkings
FROM herdeplus
WHERE animal_id = ?
  AND herdeplus_milked_mkg IS NOT NULL
  AND herdeplus_milked_mkg > 0
GROUP BY animal_id, DATE("timestamp")
ORDER BY DATE("timestamp")
"""

SQL_MLP_TEST_DAYS = """
SELECT animal_id,
       "timestamp",
       herdeplus_mlp_mkg,
       herdeplus_mlp_fat_percent,
       herdeplus_mlp_fkg,
       herdeplus_mlp_protein_percent,
       herdeplus_mlp_ekg_percent,
       herdeplus_mlp_lactose,
       herdeplus_mlp_cell_count,
       herdeplus_mlp_urea,
       herdeplus_mlp_f_e,
       herdeplus_mlp_lkg,
       herdeplus_mlp_ecm,
       herdeplus_calving_lactation
FROM herdeplus
WHERE animal_id = ?
  AND herdeplus_mlp_fat_percent IS NOT NULL
ORDER BY "timestamp"
"""

SQL_CLIMATE_DAILY = """
SELECT date("timestamp") AS day,
       AVG(temp)           AS barn_temp_mean,
       MIN(temp)           AS barn_temp_min,
       MAX(temp)           AS barn_temp_max,
       AVG(hum)            AS barn_rh_mean,
       AVG(temp_hum_index) AS barn_thi_mean,
       MIN(temp_hum_index) AS barn_thi_min,
       MAX(temp_hum_index) AS barn_thi_max
FROM smaxtec_barns
WHERE "timestamp" >= ? AND "timestamp" <= ?
  AND temp IS NOT NULL AND temp_hum_index IS NOT NULL
  AND barn_id IN (1, 2)
GROUP BY date("timestamp")
ORDER BY day
"""

# ─────────────────────────────────────────────────────────────
#  « helpers »
# ─────────────────────────────────────────────────────────────

def load_tierauswahl(path: Path) -> pd.DataFrame:
    """Load and clean the collaborator-provided animal selection list."""
    xls = pd.ExcelFile(path)
    for sheet_name in xls.sheet_names:
        df = pd.read_excel(xls, sheet_name=sheet_name)
        if "animal_id" in df.columns:
            break
    df = df[df["Auswahl"] == "Ja"].copy()
    df["animal_id"] = pd.to_numeric(df["animal_id"], errors="coerce")
    df = df.dropna(subset=["animal_id"])
    df["animal_id"] = df["animal_id"].astype(int)
    df["datetime_enter"] = pd.to_datetime(df["datetime_enter"])
    df["datetime_exit"] = pd.to_datetime(df["datetime_exit"])
    df["year"] = (
        pd.to_numeric(df["Versuchsjahr"], errors="coerce")
        if "Versuchsjahr" in df.columns
        else df["datetime_enter"].dt.year
    )
    return df


def _exclude_drinking_windows(df: pd.DataFrame) -> pd.DataFrame:
    """Remove rows during and 15 min after detected drinking events."""
    if "drink_cycles" not in df.columns:
        return df
    drink_times = df.loc[df["drink_cycles"].fillna(0) > 0, "timestamp"].values
    if len(drink_times) == 0:
        return df.drop(columns=["drink_cycles"])
    keep = np.ones(len(df), dtype=bool)
    ts = df["timestamp"].values
    for dt in drink_times:
        keep &= ~((ts >= dt) & (ts <= dt + DRINK_PAD))
    return df.loc[keep].drop(columns=["drink_cycles"])


def _get_barn(con, date_enter, date_exit, barn_cache):
    """Fetch and cache barn climate at native resolution.

    Returns barn data with a timestamp column, ready for nearest-time
    merge with rumen/respiration data via pd.merge_asof.
    """
    key = (date_enter, date_exit)
    if key in barn_cache:
        return barn_cache[key]
    barn = query_df(con, SQL_BARN, (date_enter, date_exit))
    if not barn.empty:
        barn["timestamp"] = pd.to_datetime(barn["timestamp"])
        barn = barn.sort_values("timestamp").reset_index(drop=True)
    barn_cache[key] = barn
    return barn


# ─────────────────────────────────────────────────────────────
#  « extraction routines »
# ─────────────────────────────────────────────────────────────

def extract_rumen_barn(
    con, tierauswahl: pd.DataFrame, exclude_drinking: bool = True,
) -> pd.DataFrame:
    """Extract rumen temp + barn climate for all selected animals.

    Args:
        con: Database connection.
        tierauswahl: Animal selection DataFrame.
        exclude_drinking: If True, exclude drinking events + 15 min
            padding.  If False, rely solely on smaXtec's built-in
            ``temp_without_drink_cycles`` correction.
    """
    barn_cache: dict = {}
    frames = []
    total = len(tierauswahl)

    for i, (_, row) in enumerate(tierauswahl.iterrows()):
        aid = int(row["animal_id"])
        enter = str(row["datetime_enter"])[:10]
        exit_ = str(row["datetime_exit"])[:10]
        year = int(row["year"]) if pd.notna(row.get("year")) else None

        if (i + 1) % 20 == 0 or i == 0:
            log.info("  [%d/%d] Rumen data: animal %d (%s)", i + 1, total, aid, year)

        rumen = query_df(con, SQL_RUMEN, (aid, enter, exit_))
        if rumen.empty:
            continue
        rumen["timestamp"] = pd.to_datetime(rumen["timestamp"])
        if exclude_drinking:
            rumen = _exclude_drinking_windows(rumen)
        elif "drink_cycles" in rumen.columns:
            rumen = rumen.drop(columns=["drink_cycles"])
        if rumen.empty:
            continue
        rumen = rumen.sort_values("timestamp").reset_index(drop=True)

        barn = _get_barn(con, enter, exit_, barn_cache)
        if barn.empty:
            continue

        # Nearest-time merge: each rumen reading gets the closest barn
        # reading within 30 min (works at any barn sensor resolution)
        df = pd.merge_asof(
            rumen, barn, on="timestamp",
            tolerance=pd.Timedelta("30min"),
            direction="nearest",
        )
        df = df.dropna(subset=["barn_temp", "barn_thi"])
        if df.empty:
            continue

        hour = df["timestamp"].dt.hour
        df = df[~hour.between(4, 7) & ~hour.between(16, 19)]
        df["year"] = year
        df["date_enter"] = enter
        df["date_exit"] = exit_
        frames.append(df[["animal_id", "timestamp", "year", "date_enter",
                          "date_exit", "body_temp", "barn_temp", "barn_thi"]])

    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def extract_respiration_barn(con, tierauswahl: pd.DataFrame) -> pd.DataFrame:
    """Extract respiration + barn climate for all selected animals."""
    barn_cache: dict = {}
    frames = []
    total = len(tierauswahl)

    for i, (_, row) in enumerate(tierauswahl.iterrows()):
        aid = int(row["animal_id"])
        enter = str(row["datetime_enter"])[:10]
        exit_ = str(row["datetime_exit"])[:10]
        year = int(row["year"]) if pd.notna(row.get("year")) else None

        resp = query_df(con, SQL_RESPIRATION, (aid, enter, exit_))
        if resp.empty:
            continue
        if (i + 1) % 50 == 0:
            log.info("  [%d/%d] Resp data: animal %d (%s)", i + 1, total, aid, year)

        resp["timestamp"] = pd.to_datetime(resp["timestamp"])
        resp = resp.sort_values("timestamp").reset_index(drop=True)

        barn = _get_barn(con, enter, exit_, barn_cache)
        if barn.empty:
            continue

        df = pd.merge_asof(
            resp, barn, on="timestamp",
            tolerance=pd.Timedelta("30min"),
            direction="nearest",
        )
        df = df.dropna(subset=["barn_temp", "barn_thi"])
        if df.empty:
            continue

        hour = df["timestamp"].dt.hour
        df = df[~hour.between(4, 7) & ~hour.between(16, 19)]
        df["year"] = year
        df["date_enter"] = enter
        df["date_exit"] = exit_
        frames.append(df[["animal_id", "timestamp", "year", "date_enter",
                          "date_exit", "resp_rate", "barn_temp", "barn_thi"]])

    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def extract_production(con, tierauswahl: pd.DataFrame) -> pd.DataFrame:
    """Extract milk yield and lactation number per animal."""
    records = []
    for _, row in tierauswahl.iterrows():
        aid = int(row["animal_id"])
        enter = str(row["datetime_enter"])[:10]
        exit_ = str(row["datetime_exit"])[:10]
        year = int(row["year"]) if pd.notna(row.get("year")) else None
        rec = {"animal_id": aid, "year": year, "date_enter": enter}

        milk = query_df(con, SQL_MILK_YIELD, (aid, enter, exit_))
        rec["mean_milk_yield_kg"] = (
            milk.iloc[0]["mean_milk_yield_kg"] if not milk.empty else np.nan
        )
        lac = query_df(con, SQL_LACTATION_NR, (aid,))
        rec["lactation_nr"] = (
            int(lac.iloc[0]["lactation_nr"])
            if not lac.empty and pd.notna(lac.iloc[0]["lactation_nr"])
            else np.nan
        )
        records.append(rec)
    return pd.DataFrame(records)


def extract_daily_milk_yield(con, tierauswahl: pd.DataFrame) -> pd.DataFrame:
    """Extract daily milk yield (sum of milkings per day) per animal.

    Each row is one animal-day with the total yield across all milkings
    that day.  Typically 2-3 milkings per day.

    Args:
        con: Database connection.
        tierauswahl: Cleaned Tierauswahl DataFrame.

    Returns:
        DataFrame with animal_id, year, date, daily_yield_kg, n_milkings.
    """
    frames = []
    for _, row in tierauswahl.iterrows():
        aid = int(row["animal_id"])
        enter = str(row["datetime_enter"])[:10]
        exit_ = str(row["datetime_exit"])[:10]
        year = int(row["year"]) if pd.notna(row.get("year")) else None

        df = query_df(con, SQL_DAILY_MILK_YIELD, (aid, enter, exit_))
        if not df.empty:
            df["year"] = year
            frames.append(df)

    if not frames:
        return pd.DataFrame()
    result = pd.concat(frames, ignore_index=True)
    result["date"] = pd.to_datetime(result["date"]).dt.date
    return result


def extract_daily_milk_yield_full(
    con, tierauswahl: pd.DataFrame,
) -> pd.DataFrame:
    """Full-history daily yield per Tierauswahl animal.

    Unlike :func:`extract_daily_milk_yield` this is not restricted
    to each animal's ``datetime_enter → datetime_exit`` observation
    window — we pull every day on which the animal has a valid
    ``herdeplus_milked_mkg`` entry, across her entire HerdePlus
    history in the DB.

    This is the data source for the full-history Wood (1967) fit
    used by :func:`digimuh.stats_lactation_curve.compute_wood_residuals`:
    a per-lactation curve converges far more often when it sees the
    full pre-peak / post-peak shape instead of a summer slice.

    Args:
        con: DB connection.
        tierauswahl: Cleaned Tierauswahl DataFrame (provides the
            ``animal_id`` set to query).

    Returns:
        DataFrame with ``animal_id``, ``date`` (datetime64[ns]),
        ``daily_yield_kg``, ``n_milkings``.
    """
    animal_ids = sorted({int(a) for a in tierauswahl["animal_id"].unique()})
    frames = []
    for aid in animal_ids:
        df = query_df(con, SQL_DAILY_MILK_YIELD_FULL, (aid,))
        if not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame(
            columns=["animal_id", "date", "daily_yield_kg", "n_milkings"])
    result = pd.concat(frames, ignore_index=True)
    result["date"] = pd.to_datetime(result["date"])
    return result


def extract_mlp_test_days(
    con, tierauswahl: pd.DataFrame,
) -> pd.DataFrame:
    """Per-(animal, test-day) MLP composition records for Tierauswahl animals.

    HerdePlus stores two interleaved channels in the same table:

    * high-frequency **per-milking** events
      (``herdeplus_milked_*`` populated, ``herdeplus_mlp_*`` null)
    * monthly **MLP (Milchleistungsprüfung) test-day** analytics
      (``herdeplus_mlp_*`` populated, ``herdeplus_milked_*`` typically null)

    This helper extracts the MLP channel: fat %, protein %, fat-kg,
    lactose, somatic cell count, urea, fat-to-protein ratio,
    lactose-kg, and energy-corrected milk (ECM) at one-month
    resolution.  These are the standard dairy health / composition
    analytics that would otherwise be missing from the analysis
    pipeline.

    Args:
        con: DB connection.
        tierauswahl: Cleaned Tierauswahl DataFrame.

    Returns:
        DataFrame with one row per (animal_id, test-day timestamp)
        and the eleven MLP columns plus ``herdeplus_calving_lactation``.
    """
    animal_ids = sorted({int(a) for a in tierauswahl["animal_id"].unique()})
    frames = []
    for aid in animal_ids:
        df = query_df(con, SQL_MLP_TEST_DAYS, (aid,))
        if not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    result = pd.concat(frames, ignore_index=True)
    result["timestamp"] = pd.to_datetime(result["timestamp"])
    return result


def extract_calvings(con) -> pd.DataFrame:
    """Extract all calving-confirmation events herd-wide.

    Pulls one row per calving from ``smaxtec_events`` across all
    animals — not just the Tierauswahl — so that downstream Wood
    (1967) lactation-curve fits can determine DIM for any cow-day
    including lactations that started before the observation window.

    Returns:
        DataFrame with ``animal_id`` and ``calving_date``
        (pandas datetime64[ns], day resolution).
    """
    df = query_df(con, SQL_CALVINGS, ())
    if df.empty:
        return df
    df["calving_date"] = pd.to_datetime(df["calving_date"])
    return df


def extract_climate(con, tierauswahl: pd.DataFrame) -> pd.DataFrame:
    """Extract daily barn climate for each summer in the dataset."""
    years = sorted(tierauswahl["year"].dropna().unique().astype(int))
    frames = []
    for year in years:
        df = query_df(con, SQL_CLIMATE_DAILY, (f"{year}-06-01", f"{year}-09-30"))
        if not df.empty:
            df["year"] = year
            df["day"] = pd.to_datetime(df["day"])
            frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# ─────────────────────────────────────────────────────────────
#  « main »
# ─────────────────────────────────────────────────────────────

def main() -> None:
    from digimuh.config import load_config, print_config

    parser = argparse.ArgumentParser(description="Extract data for broken-stick analysis")
    parser.add_argument("--db", type=Path, default=None,
                        help="Path to cow.db (default: from config)")
    parser.add_argument("--tierauswahl", type=Path, default=None,
                        help="Path to Tierauswahl.xlsx (default: from config)")
    parser.add_argument("--out", type=Path, default=None,
                        help="Output directory (default: from config)")
    parser.add_argument(
        "--smaxtec-drink-correction", action="store_true", default=None,
        help="Use only smaXtec's built-in temp_without_drink_cycles "
             "correction for drinking events.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(name)-20s │ %(message)s",
        datefmt="%H:%M:%S",
    )

    cfg = load_config(args)
    print_config(cfg)

    if cfg.database is None:
        parser.error("--db is required (or set 'database' in config)")
    if cfg.tierauswahl is None:
        parser.error("--tierauswahl is required (or set 'tierauswahl' in config)")

    cfg.output.mkdir(parents=True, exist_ok=True)

    exclude_drinking = not cfg.smaxtec_drink_correction
    log.info("Drink correction: %s",
             "smaXtec only (temp_without_drink_cycles)"
             if not exclude_drinking
             else "smaXtec + 15-min post-drinking exclusion window")

    con = connect_db(cfg.database)
    ta = load_tierauswahl(cfg.tierauswahl)
    ta.to_csv(resolve_output(cfg.output, "tierauswahl.csv"), index=False)
    log.info("Tierauswahl: %d entries, years %s",
             len(ta), sorted(ta["year"].dropna().unique().astype(int)))

    log.info("Extracting rumen + barn data …")
    rumen = extract_rumen_barn(con, ta, exclude_drinking=exclude_drinking)
    rumen.to_csv(resolve_output(cfg.output, "rumen_barn.csv"), index=False)
    log.info("  → %d rows, %d animals", len(rumen), rumen["animal_id"].nunique())

    log.info("Extracting respiration + barn data …")
    resp = extract_respiration_barn(con, ta)
    resp.to_csv(resolve_output(cfg.output, "respiration_barn.csv"), index=False)
    log.info("  → %d rows, %d animals", len(resp), resp["animal_id"].nunique())

    log.info("Extracting production data …")
    prod = extract_production(con, ta)
    prod.to_csv(resolve_output(cfg.output, "production.csv"), index=False)
    n_milk = prod["mean_milk_yield_kg"].notna().sum()
    n_lac = prod["lactation_nr"].notna().sum()
    log.info("  → %d with milk yield, %d with lactation nr", n_milk, n_lac)

    log.info("Extracting daily milk yield …")
    daily_yield = extract_daily_milk_yield(con, ta)
    daily_yield.to_csv(resolve_output(cfg.output, "daily_milk_yield.csv"),
                       index=False)
    log.info("  → %d daily records, %d animals",
             len(daily_yield), daily_yield["animal_id"].nunique() if not daily_yield.empty else 0)

    log.info("Extracting full-history daily milk yield "
             "(for Wood lactation-curve fitting) …")
    daily_yield_full = extract_daily_milk_yield_full(con, ta)
    daily_yield_full.to_csv(
        resolve_output(cfg.output, "daily_milk_yield_full.csv"), index=False)
    log.info(
        "  → %d daily records across full DB history, %d animals",
        len(daily_yield_full),
        daily_yield_full["animal_id"].nunique()
        if not daily_yield_full.empty else 0,
    )

    log.info("Extracting MLP test-day composition (fat%%, protein%%, "
             "SCC, urea, ECM, …) …")
    mlp = extract_mlp_test_days(con, ta)
    mlp.to_csv(resolve_output(cfg.output, "mlp_test_days.csv"), index=False)
    log.info(
        "  → %d MLP test-day rows, %d animals",
        len(mlp),
        mlp["animal_id"].nunique() if not mlp.empty else 0,
    )

    log.info("Extracting barn climate …")
    climate = extract_climate(con, ta)
    climate.to_csv(resolve_output(cfg.output, "climate_daily.csv"),
                   index=False)
    log.info("  → %d daily records", len(climate))

    log.info("Extracting calving events …")
    calvings = extract_calvings(con)
    calvings.to_csv(resolve_output(cfg.output, "calvings.csv"), index=False)
    log.info("  → %d calvings from %d animals",
             len(calvings),
             calvings["animal_id"].nunique() if not calvings.empty else 0)

    con.close()
    log.info("Extraction complete. CSVs in: %s", cfg.output)


if __name__ == "__main__":
    main()

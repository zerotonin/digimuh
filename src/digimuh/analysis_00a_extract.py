#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  Analysis 0a — Data extraction ¹                             ║
# ║  « DB → CSV for broken-stick analysis »                      ║
# ╚══════════════════════════════════════════════════════════════╝
"""Extract rumen temperature, respiration, barn climate, and production
data from the DigiMuh database into analysis-ready CSVs.

Run this once.  Downstream scripts (``analysis_00b_stats`` and
``analysis_00c_plots``) read the CSVs and never touch the database.

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

SQL_LACTATION_NR = """
SELECT MAX(CAST(herdeplus_calving_lactation AS INTEGER)) AS lactation_nr
FROM herdeplus
WHERE animal_id = ?
  AND herdeplus_calving_lactation IS NOT NULL
  AND herdeplus_calving_lactation > 0
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
    """Fetch and cache hourly barn climate."""
    key = (date_enter, date_exit)
    if key in barn_cache:
        return barn_cache[key]
    barn = query_df(con, SQL_BARN, (date_enter, date_exit))
    if not barn.empty:
        barn["timestamp"] = pd.to_datetime(barn["timestamp"])
        barn["hour_key"] = barn["timestamp"].dt.floor("h")
        barn = barn.groupby("hour_key").agg(
            barn_temp=("barn_temp", "mean"),
            barn_thi=("barn_thi", "mean"),
        ).reset_index()
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
        rumen["hour_key"] = rumen["timestamp"].dt.floor("h")

        barn = _get_barn(con, enter, exit_, barn_cache)
        if barn.empty:
            continue

        df = rumen.merge(barn, on="hour_key", how="inner")
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
        resp["hour_key"] = resp["timestamp"].dt.floor("h")

        barn = _get_barn(con, enter, exit_, barn_cache)
        if barn.empty:
            continue

        df = resp.merge(barn, on="hour_key", how="inner")
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
    parser = argparse.ArgumentParser(description="Extract data for broken-stick analysis")
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--tierauswahl", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("results/broken_stick"))
    parser.add_argument(
        "--smaxtec-drink-correction", action="store_true",
        help="Use only smaXtec's built-in temp_without_drink_cycles "
             "correction for drinking events.  Default behaviour applies "
             "an additional 15-min post-drinking exclusion window on top "
             "of smaXtec's correction to remove residual recovery artifacts.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(name)-20s │ %(message)s",
        datefmt="%H:%M:%S",
    )
    args.out.mkdir(parents=True, exist_ok=True)

    exclude_drinking = not args.smaxtec_drink_correction
    log.info("Drink correction: %s",
             "smaXtec only (temp_without_drink_cycles)"
             if not exclude_drinking
             else "smaXtec + 15-min post-drinking exclusion window")

    con = connect_db(args.db)
    ta = load_tierauswahl(args.tierauswahl)
    ta.to_csv(args.out / "tierauswahl.csv", index=False)
    log.info("Tierauswahl: %d entries, years %s",
             len(ta), sorted(ta["year"].dropna().unique().astype(int)))

    log.info("Extracting rumen + barn data …")
    rumen = extract_rumen_barn(con, ta, exclude_drinking=exclude_drinking)
    rumen.to_csv(args.out / "rumen_barn.csv", index=False)
    log.info("  → %d rows, %d animals", len(rumen), rumen["animal_id"].nunique())

    log.info("Extracting respiration + barn data …")
    resp = extract_respiration_barn(con, ta)
    resp.to_csv(args.out / "respiration_barn.csv", index=False)
    log.info("  → %d rows, %d animals", len(resp), resp["animal_id"].nunique())

    log.info("Extracting production data …")
    prod = extract_production(con, ta)
    prod.to_csv(args.out / "production.csv", index=False)
    n_milk = prod["mean_milk_yield_kg"].notna().sum()
    n_lac = prod["lactation_nr"].notna().sum()
    log.info("  → %d with milk yield, %d with lactation nr", n_milk, n_lac)

    log.info("Extracting barn climate …")
    climate = extract_climate(con, ta)
    climate.to_csv(args.out / "climate_daily.csv", index=False)
    log.info("  → %d daily records", len(climate))

    con.close()
    log.info("Extraction complete. CSVs in: %s", args.out)


if __name__ == "__main__":
    main()

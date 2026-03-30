# Database validation report

**Date:** 31 March 2026
**Database:** `cow.db` (114.62 GB)
**Ingestion script:** `digimuh.ingest_cow_db` v0.1.0
**Validated by:** `digimuh.validate_db`


## Table row counts

All 16 tables present and populated.

| Table | Rows | Status |
|---|---|---|
| animals | 5,572 | OK |
| sensors | 22 | OK |
| barns | 4 | OK |
| source_files | 4,312 | OK |
| smaxtec_derived | 762,042,093 | OK |
| gouna | 14,222,031 | OK |
| lorawan | 4,395,637 | OK |
| smaxtec_barns | 1,204,242 | OK |
| herdeplus | 983,384 | OK |
| smaxtec_water_intake | 503,087 | OK |
| bcs | 350,023 | OK |
| hobo_weather | 328,282 | OK |
| diseases | 100,821 | OK |
| allocations | 92,152 | OK |
| smaxtec_events | 62,970 | OK |
| dwd_weather | 1,323 | OK |

**Total:** ~799 million rows across 12 fact tables and 4 dimension tables.


## Null rates

| Table | Column | Null % | Interpretation |
|---|---|---|---|
| smaxtec_derived | `temp` | 90.5% | **Expected.** The smaxtec_derived table contains one row per measurement event.  Most rows correspond to motility readings (~10 s intervals) that do not include a temperature sample (~10 min intervals).  ~72 M temperature values are present. |
| smaxtec_derived | `ph` | 99.9% | **Expected.** Only ~6–10% of the herd received pH boluses.  ~1.1 M pH readings are present — sufficient for SARA analysis on the pH-bolus subset. |
| smaxtec_derived | `act` | 90.5% | **Expected.** Same sampling structure as temperature — activity is reported at ~10 min intervals, not every motility tick.  ~72 M activity values present. |
| smaxtec_derived | `rum_index` | 90.6% | **Expected.** Rumination index follows the same reporting cadence.  ~71.5 M values present. |
| smaxtec_derived | `mot_period` | 27.7% | **OK.** Motility is the highest-frequency signal.  ~551 M contraction interval readings present. |
| herdeplus | `herdeplus_milked_mkg` | 1.9% | **OK.** Small number of milking records without yield (sensor timeout or dry-off events). |
| gouna | `respirationfrequency` | 0.0% | **OK.** Complete. |
| smaxtec_water_intake | `water_intake_liter` | 0.0% | **OK.** Complete. |

**Conclusion:** No unexpected nulls.  The high null rates in smaxtec_derived
reflect the multi-rate sampling architecture of the bolus (motility at ~10 s,
temperature/activity/rumination at ~10 min, pH only on equipped animals), not
data loss.


## Value range checks

| Table | Column | Min | Max | Mean | Expected range | Assessment |
|---|---|---|---|---|---|---|
| smaxtec_derived | `temp` | 0.00 | 98.98 | 39.01 | 30–45 °C | **Sensor artefacts.** Mean of 39.01 °C is physiologically correct.  Values of 0.00 (sensor initialisation) and 98.98 (likely encoding error) are outliers to be filtered in analysis.  Recommend: exclude temp < 30 or temp > 43. |
| smaxtec_derived | `ph` | 3.13 | 26.64 | 6.15 | 3.0–8.5 pH | **Sensor artefact.** Mean of 6.15 is correct.  Max of 26.64 is an encoding error.  Recommend: exclude ph > 9.0. |
| smaxtec_derived | `act` | 0.00 | 54.71 | 3.10 | 0–500 | **OK.** Within range. |
| herdeplus | `milked_mkg` | 0.01 | 65.45 | 16.99 | 0–80 kg | **OK.** 65 kg is high but plausible for a high-yielding Holstein at peak lactation. |
| herdeplus | `mlp_fat_percent` | 0.00 | 10.98 | 3.51 | 1.0–9.0% | **Mostly OK.** Mean of 3.51% is normal.  Zero values are non-test days stored as 0 instead of NULL.  Max of 10.98% is rare but physiologically possible (colostrum, very late lactation).  Recommend: exclude fat_percent = 0 when analysing MLP data. |
| herdeplus | `mlp_protein_percent` | 0.00 | 7.26 | 3.07 | 1.5–6.0% | **Same pattern.** Mean of 3.07% is normal.  Zeros are non-test days.  7.26% is high but possible in colostrum. |
| herdeplus | `mlp_f_e` | 0.00 | 2.78 | 0.97 | 0.5–3.0 | **Same pattern.** Mean of 0.97 is suspicious (should be ~1.15 for healthy cows) — the zeros from non-test days drag the mean down.  After filtering zeros, expect mean ~1.15–1.25. |
| herdeplus | `mlp_cell_count` | 0.00 | 32,673 | 277.3 | 0–10,000 | **Extreme outlier.** 32,673 × 1000 cells/mL indicates severe clinical mastitis.  Mean of 277 is normal.  Value is physiologically possible but should be flagged in analysis. |
| bcs | `bcs_wert` | 1.10 | 5.00 | 3.23 | 1.0–5.0 | **OK.** Full range of BCS scores present. |
| gouna | `respirationfrequency` | 0.00 | 32,796 | 29.65 | 5–120 bpm | **Sensor artefact.** Values > 255 are encoding errors (likely 16-bit overflow, not just 0xFF).  Mean of 29.65 is correct after implicit filtering.  Recommend: exclude resp > 150. |
| smaxtec_barns | `temp` | −5.94 | 36.80 | 15.05 | −20–55 °C | **OK.** Winter lows and summer highs for a German barn. |
| smaxtec_barns | `hum` | 13.53 | 99.09 | 74.48 | 0–100% | **OK.** |
| smaxtec_barns | `temp_hum_index` | 23.92 | 84.74 | 58.52 | 20–100 | **OK.** Peak THI of 84.7 indicates severe heat stress days occurred. |
| smaxtec_water_intake | `water_intake_liter` | 0.00 | 232.00 | 93.57 | 0–300 L | **OK.** 232 L/day is high but plausible for a high-yielding cow in summer heat. |
| dwd_weather | `thi_max` | 25.90 | 81.30 | 59.37 | 20–100 | **OK.** |
| dwd_weather | `enthalpy_max` | −1.70 | 72.00 | 33.94 | −10–120 kJ/kg | **OK.** Negative values in winter are physically correct. |


## Temporal coverage

| Table | First record | Last record | Span |
|---|---|---|---|
| smaxtec_derived | 2021-04-28 | 2024-09-30 | 3.4 years |
| herdeplus | 2021-04-01 | 2024-09-30 | 3.5 years |
| gouna | 2022-08-24 | 2024-09-30 | 2.1 years |
| bcs | 2021-10-12 | 2024-09-30 | 3.0 years |
| smaxtec_barns | 2021-04-28 | 2024-09-30 | 3.4 years |
| smaxtec_events | 2019-06-21 | 2024-10-13 | 5.3 years |
| smaxtec_water_intake | 2021-04-28 | 2024-09-30 | 3.4 years |
| lorawan | 2021-07-06 | 2024-09-30 | 3.2 years |
| hobo_weather | 2021-08-23 | 2024-09-30 | 3.1 years |
| dwd_weather | 2021-03-01 | 2024-10-13 | 3.6 years |

**Notes:**

- smaXtec events extend back to June 2019 (historical calving/insemination records pre-dating the main monitoring period).
- Gouna respiration monitoring started ~16 months after smaXtec — overlap period is August 2022 to September 2024.
- All core systems (smaXtec, HerdePlus, barn sensors, weather) have near-complete overlap from April 2021 onward.
- Three full summers (2022, 2023, 2024) are available for heat stress analysis.


## Referential integrity

All foreign key relationships intact.  No orphaned references.

| Fact table | FK column | Dimension table | Status |
|---|---|---|---|
| bcs | animal_id | animals | OK |
| gouna | animal_id | animals | OK |
| herdeplus | animal_id | animals | OK |
| smaxtec_derived | animal_id | animals | OK |
| smaxtec_events | animal_id | animals | OK |
| smaxtec_water_intake | animal_id | animals | OK |
| lorawan | sensor_id | sensors | OK |
| smaxtec_barns | barn_id | barns | OK |


## Recommended pre-analysis filters

Based on the artefacts identified above, all analysis scripts should apply
these filters before computing statistics:

| Column | Filter | Reason |
|---|---|---|
| `temp` | 30.0 < temp < 43.0 | Remove sensor initialisation (0.0) and encoding errors (98.98) |
| `ph` | 3.0 < ph < 9.0 | Remove encoding errors (max 26.64) |
| `respirationfrequency` | 0 < resp < 150 | Remove overflow artefacts (values up to 32,796) |
| `herdeplus_mlp_fat_percent` | > 0.0 when analysing MLP | Zeros are non-test days, not real measurements |
| `herdeplus_mlp_protein_percent` | > 0.0 when analysing MLP | Same as fat |
| `herdeplus_mlp_f_e` | > 0.0 when analysing MLP | Same — zeros deflate the mean |
| `herdeplus_mlp_cell_count` | Consider flagging > 10,000 | Extreme clinical mastitis; valid but may dominate statistics |


## Indexes

All 15 indexes confirmed present:

`idx_alloc_animal`, `idx_alloc_enter`, `idx_bcs_animal_ts`,
`idx_disease_animal`, `idx_disease_start`, `idx_dwd_dt`,
`idx_gouna_animal_ts`, `idx_herdeplus_animal_ts`, `idx_hobo_ts`,
`idx_lorawan_sensor_ts`, `idx_smaxtec_barns_barn_ts`,
`idx_smaxtec_derived_animal_ts`, `idx_smaxtec_events_animal_ts`,
`idx_smaxtec_water_animal_ts`, `idx_source_folder`.


## Database size note

The database file is 114.62 GB.  This is larger than expected for ~799 M rows
due to WAL journal overhead during ingestion (the 82 GB WAL was checkpointed
back into the main file, leaving internal free pages).  Running `VACUUM` would
compact the file to an estimated 40–60 GB but requires ~120 GB of temporary
disk space and several hours.  This is recommended but not blocking for analysis.

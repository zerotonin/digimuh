# Database structure

This document describes the SQLite database produced by `digimuh-ingest`.
The schema follows a **star schema** design: small dimension tables hold entity
metadata; large fact tables hold time-series measurements with foreign keys
pointing back to the dimensions.


## Design decisions

### Why SQLite?

The dataset is large (~8.9 GB CSV, ~842 million rows) but single-user and
read-heavy.  SQLite handles this well: no server process, a single portable
`.db` file, built-in date/time functions, and indexed queries that are orders
of magnitude faster than scanning thousands of CSV files.

### Why a star schema?

The raw data arrives as thousands of per-animal CSV files with no relational
structure.  A star schema gives us:

- **Normalised entity storage** — animal IDs, sensor names, and barn names
  each appear once in their dimension table rather than repeated in every row.
- **Provenance** — every row carries a `file_id` foreign key to the
  `source_files` dimension, so any datum can be traced back to its original CSV.
- **Extensibility** — the `animals` dimension can later be enriched with breed,
  birth date, dam/sire, etc. without touching the fact tables.

### Why timestamps stay inline (not normalised)?

Timestamps are stored directly in each fact table as ISO-8601 TEXT strings.
A timestamp dimension table would itself contain hundreds of millions of rows
(the gouna table alone has sub-second precision) and every query would need a
join through it.  SQLite's built-in `date()`, `time()`, and `strftime()`
functions operate natively on ISO-8601 text, so inline storage is both simpler
and faster.

### Why the EU ear tag as INTEGER PRIMARY KEY?

EU ear tag numbers are 15-digit integers that are permanent, unique, and never
change.  SQLite's `INTEGER PRIMARY KEY` is special: it becomes an alias for the
internal `rowid`, giving the fastest possible lookup.  A surrogate key would add
a layer of indirection (an extra join) for zero space savings — the 8-byte
integer is stored in fact tables either way.


## Dimension tables

### animals

| Column | Type | Description |
|---|---|---|
| `animal_id` | INTEGER PRIMARY KEY | EU ear tag number (15-digit); IS the rowid |

Currently contains only the ID.  Future columns: breed, birth date, sex,
dam/sire lineage.

### sensors

| Column | Type | Description |
|---|---|---|
| `sensor_id` | INTEGER PRIMARY KEY | Auto-incrementing surrogate key |
| `sensor_name` | TEXT NOT NULL UNIQUE | LoRaWAN sensor name (e.g. "CU-1") |

### barns

| Column | Type | Description |
|---|---|---|
| `barn_id` | INTEGER PRIMARY KEY | Auto-incrementing surrogate key |
| `barn_name` | TEXT NOT NULL UNIQUE | smaXtec barn sensor name (e.g. "NewBridge") |

### source_files

| Column | Type | Description |
|---|---|---|
| `file_id` | INTEGER PRIMARY KEY | Auto-incrementing surrogate key |
| `filename` | TEXT NOT NULL | Basename of the original CSV file |
| `folder` | TEXT NOT NULL | Containing folder name, or "(standalone)" |

Unique constraint on `(filename, folder)`.


## Fact tables

Every fact table includes a `file_id INTEGER NOT NULL REFERENCES source_files(file_id)` column for provenance tracing.

### allocations

Records which barn group/pen an animal was assigned to and when.

| Column | Type | FK → |
|---|---|---|
| `animal_id` | INTEGER | animals |
| `datetime_enter` | TEXT | — |
| `datetime_exit` | TEXT | — |
| `group` | INTEGER | — |
| `file_id` | INTEGER | source_files |

**Source:** `output_allocations/allocations.csv`

### diseases

Health events and diagnoses from HerdePlus.

| Column | Type | FK → |
|---|---|---|
| `animal_id` | INTEGER | animals |
| `disease_first_day` | TEXT | — |
| `disease_stop_day` | TEXT | — |
| `disease_description` | TEXT | — |
| `herdeplus_tier_id` | INTEGER | — |
| `herdeplus_gesund_id` | INTEGER | — |
| `disease_category` | TEXT | — |
| `file_id` | INTEGER | source_files |

**Source:** `herdeplus_diseases.csv` (standalone)

### herdeplus

Milking events, MLP (Milchleistungsprüfung) test-day results, and calving
records from the HerdePlus herd management system.

| Column | Type | FK → |
|---|---|---|
| `animal_id` | INTEGER | animals |
| `timestamp` | TEXT | — |
| `herdeplus_milked_duration_sec` | REAL | — |
| `herdeplus_milked_milk_flow` | REAL | — |
| `herdeplus_milked_mkg` | REAL | — |
| `herdeplus_mlp_mkg` | REAL | — |
| `herdeplus_mlp_fat_percent` | REAL | — |
| `herdeplus_mlp_fkg` | REAL | — |
| `herdeplus_mlp_protein_percent` | REAL | — |
| `herdeplus_mlp_ekg_percent` | REAL | — |
| `herdeplus_mlp_lactose` | REAL | — |
| `herdeplus_mlp_cell_count` | REAL | — |
| `herdeplus_mlp_urea` | REAL | — |
| `herdeplus_mlp_f_e` | REAL | — |
| `herdeplus_mlp_lkg` | REAL | — |
| `herdeplus_mlp_ecm` | REAL | — |
| `herdeplus_calving_lactation` | TEXT | — |
| `file_id` | INTEGER | source_files |

**Source:** `outputs_herdeplus_mlp_gemelk_kalbung/{animal_id}_herdeplus_*.csv`

### bcs

Body Condition Score assessments.

| Column | Type | FK → |
|---|---|---|
| `animal_id` | INTEGER | animals |
| `timestamp` | TEXT | — |
| `bcs_wert` | REAL | — |
| `file_id` | INTEGER | source_files |

**Source:** `outputs_bcs/{animal_id}_bcs_*.csv`

### gouna

Respiration frequency from Gouna sensors.

| Column | Type | FK → |
|---|---|---|
| `animal_id` | INTEGER | animals |
| `timestamp` | TEXT | — |
| `respirationfrequency` | REAL | — |
| `file_id` | INTEGER | source_files |

**Source:** `outputs_gouna/{animal_id}_gouna_*.csv`

### smaxtec_derived

Derived metrics from the smaXtec rumen bolus.  This is the largest table
(~824 million rows at full ingestion) with 30+ columns of computed indices
covering activity, estrus detection, calving prediction, rumination, rumen pH,
temperature, and motility.

| Column | Type | FK → |
|---|---|---|
| `animal_id` | INTEGER | animals |
| `timestamp` | TEXT | — |
| `act` | REAL | — |
| `act_decrease_index` | REAL | — |
| `act_estrus_index` | REAL | — |
| `act_estrus_preprocess` | REAL | — |
| `act_group_heat_index` | REAL | — |
| `act_group_ratio` | REAL | — |
| `act_index` | REAL | — |
| `act_pasture_index` | REAL | — |
| `calving_index` | REAL | — |
| `drink_cycles_v2` | REAL | — |
| `heat_index` | REAL | — |
| `in_reticulum` | REAL | — |
| `mot_period` | REAL | — |
| `mot_pulse_width` | REAL | — |
| `mot_pulse_width_median` | REAL | — |
| `ph` | REAL | — |
| `ph_under_58` | REAL | — |
| `rum_classification` | TEXT | — |
| `rum_dec_index` | REAL | — |
| `rum_index` | REAL | — |
| `temp` | REAL | — |
| `temp_dec_index` | REAL | — |
| `temp_group_ratio_svm_inc_index` | REAL | — |
| `temp_height_index` | REAL | — |
| `temp_inc_index` | REAL | — |
| `temp_limit_crossing` | REAL | — |
| `temp_normal_index` | REAL | — |
| `temp_svm_inc_index` | REAL | — |
| `temp_without_drink_cycles` | REAL | — |
| `mot_period_rum_6h_we` | REAL | — |
| `mot_period_not_rum_6h_we` | REAL | — |
| `file_id` | INTEGER | source_files |

**Source:** `outputs_smaxtec_derived/{animal_id}_smaxtec_derived_*.csv`

**Note:** Due to sparse data in early rows, some columns may be typed as TEXT
by the auto-inference engine.  SQLite's dynamic typing means numeric values are
still stored and compared correctly regardless of the declared affinity.

### smaxtec_events

Discrete reproductive and health events from smaXtec.

| Column | Type | FK → |
|---|---|---|
| `animal_id` | INTEGER | animals |
| `timestamp` | TEXT | — |
| `cow` | INTEGER | — |
| `event_type` | TEXT | — |
| `value` | REAL | — |
| `file_id` | INTEGER | source_files |

**Source:** `outputs_smaxtec_events/{animal_id}_events.csv`

**Note:** The `cow` column contains the animal_id redundantly (the canonical
`animal_id` is extracted from the filename).  Both are kept for completeness.

### smaxtec_water_intake

Daily water intake estimates derived from rumen temperature drops.

| Column | Type | FK → |
|---|---|---|
| `animal_id` | INTEGER | animals |
| `timestamp` | TEXT | — |
| `water_intake_liter` | REAL | — |
| `file_id` | INTEGER | source_files |

**Source:** `outputs_smaxtec_water_intake/{animal_id}_smaxtec_derived_*.csv`

### smaxtec_barns

Barn climate measurements from smaXtec barn sensors.

| Column | Type | FK → |
|---|---|---|
| `barn_id` | INTEGER | barns |
| `timestamp` | TEXT | — |
| `rawtemp` | REAL | — |
| `rawhum` | REAL | — |
| `temp` | REAL | — |
| `hum` | REAL | — |
| `temp_hum_index` | REAL | — |
| `file_id` | INTEGER | source_files |

**Source:** `outputs_smaxtec_barns/{barn_name}_smaxtec_raw_*.csv`

### lorawan

LoRaWAN environmental sensor readings (battery and current).

| Column | Type | FK → |
|---|---|---|
| `sensor_id` | INTEGER | sensors |
| `timestamp` | TEXT | — |
| `battery_level` | REAL | — |
| `current_ampere` | REAL | — |
| `file_id` | INTEGER | source_files |

**Source:** `outputs_lorawan/{sensor_name}_LoRaWAN_raw_*.csv`

### dwd_weather

Daily weather summaries from the Deutscher Wetterdienst (DWD).

| Column | Type | FK → |
|---|---|---|
| `dt` | TEXT | — |
| `thi_max` | REAL | — |
| `qb_thi` | TEXT | — |
| `num_values_thi` | INTEGER | — |
| `enthalpy_max` | REAL | — |
| `qb_enthalpy` | TEXT | — |
| `num_values_enthalpy` | INTEGER | — |
| `file_id` | INTEGER | source_files |

**Source:** `outputs_dwd.csv` (standalone)

### hobo_weather

Weather station readings from HOBO loggers.  Column names include sensor
serial numbers from the HOBO export format.

| Column | Type | FK → |
|---|---|---|
| `datetime` | TEXT | — |
| `21136553_b_battery_v` | REAL | — |
| `21141733_1_temperature` | REAL | — |
| `21141733_2_rh` | REAL | — |
| `21141733_3_dew_point` | REAL | — |
| `21141733_b_battery_level` | REAL | — |
| `21141735_1_solar_radiation` | REAL | — |
| `21141735_b_battery_level` | REAL | — |
| `21141737_1_wetness` | REAL | — |
| `21141737_b_battery_level` | REAL | — |
| `21141734_1_wind_speed` | REAL | — |
| `21141734_2_gust_speed` | REAL | — |
| `21141734_3_wind_direction` | REAL | — |
| `21141734_b_battery_level` | REAL | — |
| `file_id` | INTEGER | source_files |

**Source:** `outputs_hobo/hobo_exports_*.csv`


## Indexes

All indexes are created after bulk insertion to avoid write overhead.

| Index | Table | Columns |
|---|---|---|
| `idx_bcs_animal_ts` | bcs | `animal_id, timestamp` |
| `idx_gouna_animal_ts` | gouna | `animal_id, timestamp` |
| `idx_herdeplus_animal_ts` | herdeplus | `animal_id, timestamp` |
| `idx_smaxtec_derived_animal_ts` | smaxtec_derived | `animal_id, timestamp` |
| `idx_smaxtec_events_animal_ts` | smaxtec_events | `animal_id, timestamp` |
| `idx_smaxtec_water_animal_ts` | smaxtec_water_intake | `animal_id, timestamp` |
| `idx_alloc_animal` | allocations | `animal_id` |
| `idx_alloc_enter` | allocations | `datetime_enter` |
| `idx_disease_animal` | diseases | `animal_id` |
| `idx_disease_start` | diseases | `disease_first_day` |
| `idx_lorawan_sensor_ts` | lorawan | `sensor_id, timestamp` |
| `idx_smaxtec_barns_barn_ts` | smaxtec_barns | `barn_id, timestamp` |
| `idx_dwd_dt` | dwd_weather | `dt` |
| `idx_hobo_ts` | hobo_weather | `datetime` |
| `idx_source_folder` | source_files | `folder` |

The composite `(entity_id, timestamp)` indexes accelerate the primary query
pattern: "give me all measurements for entity X between dates A and B".


## Planned additions

- **Lactation view** combining herdeplus milking + BCS + disease windows per lactation cycle
- **Materialised summary tables** for daily/weekly aggregates (if view performance is insufficient)


## Analysis views

The file `src/digimuh/create_views.sql` defines a three-layer view hierarchy.
Views are created automatically when any analysis script connects to the database.

### Layer 0 — hourly aggregates

| View | Purpose | Grouped by |
|---|---|---|
| `v_smaxtec_hourly` | Hourly means of temp, activity, rumination, motility, pH | animal × day × hour |

### Layer 1 — daily summaries

These views aggregate each data source into one row per animal (or sensor/barn) per day.

| View | Source table(s) | Key columns |
|---|---|---|
| `v_smaxtec_daily` | smaxtec_derived | temp (mean/min/max/range), activity, rumination, motility, pH, drinking, estrus/calving indices |
| `v_herdeplus_daily` | herdeplus | total milk yield, mean duration/flow, MLP test-day values (fat, protein, FPR, SCC, urea, lactose, ECM) |
| `v_gouna_daily` | gouna | mean/min/max respiration frequency |
| `v_water_daily` | smaxtec_water_intake | total water intake (litres) |
| `v_bcs_latest` | bcs | BCS value and assessment date |
| `v_barn_daily` | smaxtec_barns | barn temp/humidity/THI (mean and max) |

### Layer 2 — analysis-specific joins

Each analysis view joins the daily summaries needed for a specific research question.

| View | Joins | Purpose |
|---|---|---|
| `v_analysis_ketosis` | herdeplus + smaxtec + water + diseases | FPR, rumination, milk yield, rumen pH, with disease ground truth.  Includes a `fpr_flag` column (+1 = energy deficit, −1 = acidosis, 0 = normal) and `is_sick` binary from disease records. |
| `v_analysis_heat_stress` | smaxtec + water + gouna + herdeplus + DWD weather | Rumen temp, respiration, activity, water intake, milk yield, all aligned with ambient THI and enthalpy. |
| `v_analysis_digestive` | smaxtec + herdeplus + water | Motility/pH daily profiles alongside MLP test-day composition for time-lagged cross-correlation analysis. |
| `v_analysis_circadian` | smaxtec_hourly + diseases | Hourly temp/activity/rumination profiles per animal-day with disease status for circadian rhythm analysis. |
| `v_analysis_motility` | smaxtec_derived (minimal aggregation) | Raw motility series (contraction intervals, pulse widths) with concurrent pH, rumination, and temperature for entropy computation. |

### Example queries using views

```sql
-- Daily milk yield vs. heat load for a specific cow
SELECT day, milk_yield_kg, rumen_temp_mean, dwd_thi_max, water_liter
FROM v_analysis_heat_stress
WHERE animal_id = 276001260919234
ORDER BY day;

-- MLP test days with FPR above ketosis threshold
SELECT animal_id, day, mlp_fpr, rum_index_mean, milk_yield_kg, disease_category
FROM v_analysis_ketosis
WHERE fpr_flag = 1
ORDER BY day;

-- Hourly temperature profile for circadian analysis
SELECT hour, AVG(temp_clean_mean) AS temp, AVG(act_index_mean) AS activity
FROM v_analysis_circadian
WHERE animal_id = 276001260919234
  AND day = '2023-07-15'
GROUP BY hour
ORDER BY hour;
```

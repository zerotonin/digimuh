-- ╔══════════════════════════════════════════════════════════════╗
-- ║  DigiMuh analysis views                                     ║
-- ║  « layered views for welfare & production analytics »       ║
-- ╠══════════════════════════════════════════════════════════════╣
-- ║  Layer 0 — hourly summaries (for circadian analysis)        ║
-- ║  Layer 1 — daily summaries per animal per data source       ║
-- ║  Layer 2 — cross-table joins for each analysis question     ║
-- ╚══════════════════════════════════════════════════════════════╝
--
--  Run once after ingestion:
--    sqlite3 cow.db < create_views.sql
--
--  Or from Python:
--    con.executescript(open("create_views.sql").read())


-- ═════════════════════════════════════════════════════════════
--  LAYER 0 — hourly aggregates (needed for circadian analysis)
-- ═════════════════════════════════════════════════════════════

DROP VIEW IF EXISTS v_smaxtec_hourly;
CREATE VIEW v_smaxtec_hourly AS
SELECT
    animal_id,
    date("timestamp")                              AS day,
    CAST(strftime('%H', "timestamp") AS INTEGER)   AS hour,
    AVG(CAST("temp" AS REAL))                      AS temp_mean,
    AVG(CAST("temp_without_drink_cycles" AS REAL)) AS temp_clean_mean,
    AVG(CAST("act" AS REAL))                       AS act_mean,
    AVG(CAST("act_index" AS REAL))                 AS act_index_mean,
    AVG(CAST("rum_index" AS REAL))                 AS rum_index_mean,
    AVG(CAST("mot_period" AS REAL))                AS mot_period_mean,
    AVG(CAST("mot_pulse_width" AS REAL))           AS mot_pw_mean,
    AVG(CAST("ph" AS REAL))                        AS ph_mean,
    COUNT(*)                                       AS n_readings
FROM smaxtec_derived
WHERE "timestamp" IS NOT NULL
GROUP BY animal_id, day, hour;


-- ═════════════════════════════════════════════════════════════
--  LAYER 1 — daily animal summaries
-- ═════════════════════════════════════════════════════════════

-- ── smaXtec derived: daily per animal ───────────────────────

DROP VIEW IF EXISTS v_smaxtec_daily;
CREATE VIEW v_smaxtec_daily AS
SELECT
    animal_id,
    date("timestamp")                              AS day,

    -- temperature
    AVG(CAST("temp" AS REAL))                      AS temp_mean,
    MIN(CAST("temp" AS REAL))                      AS temp_min,
    MAX(CAST("temp" AS REAL))                      AS temp_max,
    MAX(CAST("temp" AS REAL))
      - MIN(CAST("temp" AS REAL))                  AS temp_range,
    AVG(CAST("temp_without_drink_cycles" AS REAL)) AS temp_clean_mean,

    -- activity
    AVG(CAST("act" AS REAL))                       AS act_mean,
    AVG(CAST("act_index" AS REAL))                 AS act_index_mean,

    -- rumination
    AVG(CAST("rum_index" AS REAL))                 AS rum_index_mean,

    -- motility
    AVG(CAST("mot_period" AS REAL))                AS mot_period_mean,
    AVG(CAST("mot_pulse_width" AS REAL))           AS mot_pw_mean,
    AVG(CAST("mot_pulse_width_median" AS REAL))    AS mot_pw_median_mean,

    -- rumen pH
    AVG(CAST("ph" AS REAL))                        AS ph_mean,
    MIN(CAST("ph" AS REAL))                        AS ph_min,
    SUM(CAST("ph_under_58" AS REAL))               AS ph_under_58_total,

    -- drinking
    SUM(CAST("drink_cycles_v2" AS REAL))           AS drink_cycles_total,

    -- estrus / calving indices
    MAX(CAST("heat_index" AS REAL))                AS heat_index_max,
    MAX(CAST("calving_index" AS REAL))             AS calving_index_max,

    COUNT(*)                                       AS n_readings

FROM smaxtec_derived
WHERE "timestamp" IS NOT NULL
GROUP BY animal_id, date("timestamp");


-- ── herdeplus: daily milking summary ────────────────────────

DROP VIEW IF EXISTS v_herdeplus_daily;
CREATE VIEW v_herdeplus_daily AS
SELECT
    animal_id,
    date("timestamp")                                AS day,
    SUM("herdeplus_milked_mkg")                      AS milk_yield_kg,
    AVG("herdeplus_milked_duration_sec")              AS milk_duration_mean_sec,
    AVG("herdeplus_milked_milk_flow")                 AS milk_flow_mean,
    -- MLP test-day values (sparse — will be NULL most days)
    MAX("herdeplus_mlp_fat_percent")                  AS mlp_fat_pct,
    MAX("herdeplus_mlp_protein_percent")              AS mlp_protein_pct,
    MAX("herdeplus_mlp_f_e")                          AS mlp_fpr,
    MAX("herdeplus_mlp_cell_count")                   AS mlp_scc,
    MAX("herdeplus_mlp_urea")                         AS mlp_urea,
    MAX("herdeplus_mlp_lactose")                      AS mlp_lactose,
    MAX("herdeplus_mlp_ecm")                          AS mlp_ecm,
    MAX("herdeplus_calving_lactation")                AS lactation_nr,
    COUNT(*)                                          AS n_milkings
FROM herdeplus
WHERE "timestamp" IS NOT NULL
GROUP BY animal_id, date("timestamp");


-- ── gouna: daily respiration ────────────────────────────────

DROP VIEW IF EXISTS v_gouna_daily;
CREATE VIEW v_gouna_daily AS
SELECT
    animal_id,
    date("timestamp")      AS day,
    AVG(respirationfrequency) AS resp_mean,
    MIN(respirationfrequency) AS resp_min,
    MAX(respirationfrequency) AS resp_max,
    COUNT(*)                  AS n_readings
FROM gouna
WHERE "timestamp" IS NOT NULL
  AND respirationfrequency IS NOT NULL
GROUP BY animal_id, date("timestamp");


-- ── water intake: already daily ─────────────────────────────

DROP VIEW IF EXISTS v_water_daily;
CREATE VIEW v_water_daily AS
SELECT
    animal_id,
    date("timestamp")     AS day,
    SUM(water_intake_liter)  AS water_liter
FROM smaxtec_water_intake
WHERE "timestamp" IS NOT NULL
GROUP BY animal_id, date("timestamp");


-- ── BCS: forward-fill via last observation carried forward ──

DROP VIEW IF EXISTS v_bcs_latest;
CREATE VIEW v_bcs_latest AS
SELECT
    b.animal_id,
    b."timestamp"         AS bcs_date,
    b.bcs_wert            AS bcs_value
FROM bcs b;


-- ── barn climate: daily per barn ────────────────────────────

DROP VIEW IF EXISTS v_barn_daily;
CREATE VIEW v_barn_daily AS
SELECT
    barn_id,
    date("timestamp")          AS day,
    AVG("temp")                AS barn_temp_mean,
    MAX("temp")                AS barn_temp_max,
    AVG("hum")                 AS barn_hum_mean,
    AVG("temp_hum_index")      AS barn_thi_mean,
    MAX("temp_hum_index")      AS barn_thi_max,
    COUNT(*)                   AS n_readings
FROM smaxtec_barns
WHERE "timestamp" IS NOT NULL
GROUP BY barn_id, date("timestamp");


-- ═════════════════════════════════════════════════════════════
--  LAYER 2 — analysis-specific views
-- ═════════════════════════════════════════════════════════════

-- ┌──────────────────────────────────────────────────────────┐
-- │  ANALYSIS 1: Subclinical ketosis detection               │
-- │  FPR + rumination + milk yield + temperature + SCC       │
-- │  joined with disease ground truth                        │
-- └──────────────────────────────────────────────────────────┘

DROP VIEW IF EXISTS v_analysis_ketosis;
CREATE VIEW v_analysis_ketosis AS
SELECT
    h.animal_id,
    h.day,
    -- milking
    h.milk_yield_kg,
    h.n_milkings,
    h.mlp_fat_pct,
    h.mlp_protein_pct,
    h.mlp_fpr,
    h.mlp_scc,
    h.mlp_urea,
    h.mlp_lactose,
    h.mlp_ecm,
    h.lactation_nr,
    -- rumen state (same day)
    s.rum_index_mean,
    s.temp_mean                AS rumen_temp_mean,
    s.act_index_mean,
    s.ph_mean,
    s.ph_under_58_total,
    s.drink_cycles_total,
    -- water intake
    w.water_liter,
    -- disease ground truth: active disease window
    d.disease_description,
    d.disease_category,
    CASE
        WHEN d.disease_category IS NOT NULL
             AND d.disease_category != 'No Disease'
             AND d.disease_category != 'Healthy/remarks'
        THEN 1 ELSE 0
    END                        AS is_sick,
    -- ketosis flag from FPR
    CASE
        WHEN h.mlp_fpr > 1.4 THEN 1
        WHEN h.mlp_fpr < 1.1 THEN -1
        ELSE 0
    END                        AS fpr_flag
FROM v_herdeplus_daily h
LEFT JOIN v_smaxtec_daily s
    ON h.animal_id = s.animal_id AND h.day = s.day
LEFT JOIN v_water_daily w
    ON h.animal_id = w.animal_id AND h.day = w.day
LEFT JOIN diseases d
    ON h.animal_id = d.animal_id
    AND h.day >= d.disease_first_day
    AND h.day <= COALESCE(d.disease_stop_day, h.day);


-- ┌──────────────────────────────────────────────────────────┐
-- │  ANALYSIS 3: Heat stress — multi-sensor fusion           │
-- │  Rumen temp × THI × water intake × respiration           │
-- └──────────────────────────────────────────────────────────┘

DROP VIEW IF EXISTS v_analysis_heat_stress;
CREATE VIEW v_analysis_heat_stress AS
SELECT
    s.animal_id,
    s.day,
    -- rumen signals
    s.temp_mean              AS rumen_temp_mean,
    s.temp_clean_mean        AS rumen_temp_clean_mean,
    s.temp_max               AS rumen_temp_max,
    s.temp_range             AS rumen_temp_range,
    s.act_index_mean,
    s.rum_index_mean,
    s.drink_cycles_total,
    -- water
    w.water_liter,
    -- respiration
    g.resp_mean,
    g.resp_max,
    -- weather (DWD)
    dwd.thi_max              AS dwd_thi_max,
    dwd.enthalpy_max         AS dwd_enthalpy_max,
    -- milking (production impact)
    h.milk_yield_kg,
    h.mlp_fat_pct,
    h.mlp_protein_pct,
    -- readings count for quality filter
    s.n_readings             AS smaxtec_readings,
    g.n_readings             AS gouna_readings
FROM v_smaxtec_daily s
LEFT JOIN v_water_daily w
    ON s.animal_id = w.animal_id AND s.day = w.day
LEFT JOIN v_gouna_daily g
    ON s.animal_id = g.animal_id AND s.day = g.day
LEFT JOIN v_herdeplus_daily h
    ON s.animal_id = h.animal_id AND s.day = h.day
LEFT JOIN dwd_weather dwd
    ON s.day = dwd.dt;


-- ┌──────────────────────────────────────────────────────────┐
-- │  ANALYSIS 6: Digestive efficiency composite              │
-- │  Motility × pH × milk composition (time-lagged)          │
-- └──────────────────────────────────────────────────────────┘

DROP VIEW IF EXISTS v_analysis_digestive;
CREATE VIEW v_analysis_digestive AS
SELECT
    s.animal_id,
    s.day,
    -- motility profile
    s.mot_period_mean,
    s.mot_pw_mean,
    s.mot_pw_median_mean,
    -- rumen chemistry
    s.ph_mean,
    s.ph_min,
    s.ph_under_58_total,
    -- rumination
    s.rum_index_mean,
    -- activity (feeding proxy)
    s.act_index_mean,
    -- drinking
    s.drink_cycles_total,
    w.water_liter,
    -- MLP test-day composition (sparse, carried from herdeplus)
    h.mlp_fat_pct,
    h.mlp_protein_pct,
    h.mlp_fpr,
    h.mlp_lactose,
    h.mlp_scc,
    h.mlp_urea,
    h.mlp_ecm,
    h.milk_yield_kg,
    h.lactation_nr,
    -- quality
    s.n_readings         AS smaxtec_readings
FROM v_smaxtec_daily s
LEFT JOIN v_herdeplus_daily h
    ON s.animal_id = h.animal_id AND s.day = h.day
LEFT JOIN v_water_daily w
    ON s.animal_id = w.animal_id AND s.day = w.day;


-- ┌──────────────────────────────────────────────────────────┐
-- │  ANALYSIS 11: Circadian rhythm disruption                │
-- │  Hourly temp/activity/rumination profiles per day        │
-- └──────────────────────────────────────────────────────────┘
--
--  v_smaxtec_hourly (Layer 0) IS the view for this analysis.
--  The Python script computes Fourier amplitude and phase
--  from the 24 hourly bins per animal per day.
--
--  Additional: join disease status for ground truth.

DROP VIEW IF EXISTS v_analysis_circadian;
CREATE VIEW v_analysis_circadian AS
SELECT
    sh.animal_id,
    sh.day,
    sh.hour,
    sh.temp_mean,
    sh.temp_clean_mean,
    sh.act_mean,
    sh.act_index_mean,
    sh.rum_index_mean,
    sh.mot_period_mean,
    sh.n_readings,
    -- disease status on that day
    CASE
        WHEN d.disease_category IS NOT NULL
             AND d.disease_category != 'No Disease'
             AND d.disease_category != 'Healthy/remarks'
        THEN 1 ELSE 0
    END                      AS is_sick,
    d.disease_category
FROM v_smaxtec_hourly sh
LEFT JOIN diseases d
    ON sh.animal_id = d.animal_id
    AND sh.day >= d.disease_first_day
    AND sh.day <= COALESCE(d.disease_stop_day, sh.day);


-- ┌──────────────────────────────────────────────────────────┐
-- │  ANALYSIS 12: Motility pattern entropy                   │
-- │  Raw motility series — minimal aggregation               │
-- │  Python does the entropy computation                     │
-- └──────────────────────────────────────────────────────────┘

DROP VIEW IF EXISTS v_analysis_motility;
CREATE VIEW v_analysis_motility AS
SELECT
    animal_id,
    "timestamp",
    date("timestamp")                            AS day,
    CAST(strftime('%H', "timestamp") AS INTEGER) AS hour,
    CAST("mot_period" AS REAL)                   AS mot_period,
    CAST("mot_pulse_width" AS REAL)              AS mot_pw,
    CAST("mot_pulse_width_median" AS REAL)       AS mot_pw_median,
    CAST("ph" AS REAL)                           AS ph,
    CAST("rum_index" AS REAL)                    AS rum_index,
    CAST("temp" AS REAL)                         AS temp
FROM smaxtec_derived
WHERE "mot_period" IS NOT NULL
  AND CAST("mot_period" AS REAL) > 0;

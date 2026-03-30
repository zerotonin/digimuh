# Column dictionary

Descriptions of every column across all tables in the DigiMuh database.
Where information comes from the smaXtec parameter documentation provided
by collaborators, it is marked **(confirmed)**.  Columns still pending
clarification are marked **(?)**.


## allocations

| Column | Unit | Description |
|---|---|---|
| `animal_id` | integer ID | EU ear tag number (15-digit unique identifier) |
| `datetime_enter` | ISO-8601 datetime | Date/time the animal entered the group/pen |
| `datetime_exit` | ISO-8601 datetime | Date/time the animal left the group/pen |
| `group` | integer code | Group or pen identifier (Gruppennummer) |


## diseases

| Column | Unit | Description |
|---|---|---|
| `animal_id` | integer ID | EU ear tag number |
| `disease_first_day` | date | Date the disease or health event was first recorded |
| `disease_stop_day` | date | Date the disease or health event resolved |
| `disease_description` | free text | Textual description of the diagnosis (e.g. BVD AG neg, mastitis, ketosis) |
| `herdeplus_tier_id` | integer | HerdePlus internal animal ID (software-specific identifier) |
| `herdeplus_gesund_id` | integer | HerdePlus internal health-record ID **(?)** |
| `disease_category` | categorical text | Broad classification (e.g. No Disease, Healthy/remarks, Other/Illness) |


## herdeplus

Milking events, MLP (Milchleistungsprüfung = official milk recording) test-day
results, and calving records from the HerdePlus herd management system.

| Column | Unit | Description |
|---|---|---|
| `timestamp` | ISO-8601 datetime | Date/time of the milking event or MLP sample |
| `herdeplus_milked_duration_sec` | seconds | Duration of the milking event |
| `herdeplus_milked_milk_flow` | kg/min | Average milk flow rate during milking |
| `herdeplus_milked_mkg` | kg | Milk yield per milking event (Milchkilogramm) |
| `herdeplus_mlp_mkg` | kg | Milk yield at MLP test day |
| `herdeplus_mlp_fat_percent` | % | Milk fat content from MLP test-day sample |
| `herdeplus_mlp_fkg` | kg | Fat yield from MLP test day (Fettkilogramm) |
| `herdeplus_mlp_protein_percent` | % | Milk protein content from MLP test-day sample |
| `herdeplus_mlp_ekg_percent` | % | Protein yield percentage from MLP (Eiweisskilogramm as fraction of milk) **(?)** |
| `herdeplus_mlp_lactose` | % | Lactose content from MLP test-day sample |
| `herdeplus_mlp_cell_count` | cells/mL (x1000) | Somatic cell count from MLP (indicator of udder health / mastitis) |
| `herdeplus_mlp_urea` | mg/dL or mg/L | Milk urea nitrogen from MLP (indicator of dietary protein-energy balance) |
| `herdeplus_mlp_f_e` | dimensionless ratio | Fat-to-protein ratio (Fett-Eiweiss-Quotient); values > 1.4 suggest energy deficit / ketosis risk, < 1.1 suggests SARA |
| `herdeplus_mlp_lkg` | --- | *Unknown -- possibly Liter-Kilogramm or lactose-kg; needs confirmation* **(?)** |
| `herdeplus_mlp_ecm` | kg | Energy-corrected milk yield, standardised to 4.0% fat and 3.4% protein |
| `herdeplus_calving_lactation` | integer | Lactation number (Laktationsnummer): which calving the current lactation follows (1 = primiparous) |


## bcs

| Column | Unit | Description |
|---|---|---|
| `timestamp` | ISO-8601 datetime | Date of BCS assessment |
| `bcs_wert` | score (1-5) | Body Condition Score (Wert = value), 1-5 scale in 0.25 increments |


## gouna

| Column | Unit | Description |
|---|---|---|
| `timestamp` | ISO-8601 datetime | Date/time of respiration measurement |
| `respirationfrequency` | breaths/min | Respiration rate (Atemfrequenz) **(confirmed)** |


## lorawan

| Column | Unit | Description |
|---|---|---|
| `timestamp` | ISO-8601 datetime | Date/time of LoRaWAN sensor reading |
| `battery_level` | % | Battery charge level of the LoRaWAN sensor node |
| `current_ampere` | Ampere | Electrical current drawn by the sensor |


## smaxtec_barns

Barn climate measurements from smaXtec barn sensors.  THI is calculated
according to the NRC (1971) formula **(confirmed)**.

| Column | Unit | Description |
|---|---|---|
| `timestamp` | ISO-8601 datetime | Date/time of barn climate measurement |
| `rawtemp` | deg C | Raw (uncalibrated) temperature reading from barn sensor |
| `rawhum` | % | Raw (uncalibrated) relative humidity from barn sensor |
| `temp` | deg C | Processed/calibrated barn temperature (Lufttemperatur) |
| `hum` | % | Processed/calibrated barn relative humidity (rel. Luftfeuchtigkeit) |
| `temp_hum_index` | dimensionless index | Temperature-Humidity Index (THI) per NRC (1971); values > 68 indicate mild heat stress **(confirmed)** |


## smaxtec_events

| Column | Unit | Description |
|---|---|---|
| `timestamp` | ISO-8601 datetime | Date/time of the recorded event |
| `cow` | integer ID | EU ear tag number (redundant with animal_id from filename) |
| `event_type` | categorical text | Type of event (e.g. calving_confirmation, insemination, pregnancy_result) |
| `value` | text or numeric | Associated value or outcome for the event (may be empty) |


## smaxtec_derived

Derived metrics computed by smaXtec software from raw rumen bolus
measurements.  The bolus sits in the reticulum and measures temperature,
pH, and three-axis acceleration.  All derived values are computed from
these raw signals.

Columns marked **(direct)** are directly measured or minimally processed.
Columns marked **(derived)** are computed by smaXtec's algorithms from
the direct measurements.

### Activity

| Column | Unit | Description |
|---|---|---|
| `act` | arbitrary units | Activity level from rumen bolus accelerometer **(direct)** |
| `act_decrease_index` | index | Activity decrease relative to individual baseline **(derived)** |
| `act_estrus_index` | index | Activity-based estrus (heat) detection index **(derived)** |
| `act_estrus_preprocess` | arbitrary units | Pre-processed/filtered activity signal used as input for estrus algorithm **(derived)** |
| `act_group_heat_index` | index | Group-level heat detection; compares animal's activity to herd average **(direct, herd-referenced)** |
| `act_group_ratio` | dimensionless ratio | Ratio of this animal's activity to group mean **(direct, herd-referenced)** |
| `act_index` | index | Normalised general activity index **(derived)**; example value: 8.96 |
| `act_pasture_index` | index | Activity index adjusted for pasture vs. barn context **(derived)** |

### Reproduction and calving

| Column | Unit | Description |
|---|---|---|
| `calving_index` | index | Imminent calving prediction index (based on temperature drop and activity changes) **(derived)** |
| `heat_index` | index | Combined heat/estrus detection score (fuses activity and temperature signals) **(derived)** |

### Drinking

| Column | Unit | Description |
|---|---|---|
| `drink_cycles_v2` | count | Number of detected drinking bouts, v2 detection algorithm **(derived)** |

### Rumen state

| Column | Unit | Description |
|---|---|---|
| `in_reticulum` | 0 or 1 | Whether the bolus is confirmed to be in the reticulum **(direct)**; example: 1.0 |
| `ph` | pH units | Rumen pH **(direct)**; only available in animals with smaXtec pH bolus (~6-10% of herd) |
| `ph_under_58` | minutes or count | Duration/count of pH readings below 5.8 -- the subacute ruminal acidosis (SARA) threshold **(derived)** |
| `rum_classification` | categorical | Rumination state classification (ruminating / not ruminating) **(direct)** |
| `rum_dec_index` | index | Rumination decrease index -- flags drops relative to individual baseline **(derived)** |
| `rum_index` | contraction count | Cumulative reticulorumen contraction count **(direct)**; example value: 34,852 -- this is a raw count, not a normalised 0-1 index |

### Motility

| Column | Unit | Description |
|---|---|---|
| `mot_period` | seconds | Reticulorumen contraction interval (motility period) **(direct)** |
| `mot_pulse_width` | seconds | Duration of a single reticulorumen contraction pulse **(direct)** |
| `mot_pulse_width_median` | seconds | Median contraction pulse width over a rolling window **(derived)** |
| `mot_period_rum_6h_we` | seconds | Mean contraction interval during rumination, 6-hour window **(derived, ?)** |
| `mot_period_not_rum_6h_we` | seconds | Mean contraction interval outside rumination, 6-hour window **(derived, ?)** |

### Temperature

| Column | Unit | Description |
|---|---|---|
| `temp` | deg C | Core body (rumen/reticulum) temperature **(direct)**; example: 39.14 deg C |
| `temp_without_drink_cycles` | deg C | Rumen temperature with drinking-event cold-water artifacts filtered out **(derived)**; example: 39.11 deg C |
| `temp_normal_index` | deg C | Animal's individual baseline/reference temperature **(derived)**; example: 39.48 deg C -- this is a temperature value, not a normalised score |
| `temp_height_index` | deg C | Deviation of current temperature from the animal's baseline (temp_normal_index) **(derived)**; example: -0.37 deg C (negative = below baseline) |
| `temp_dec_index` | index | Temperature decrease index -- flags abnormal temperature drops **(derived)** |
| `temp_inc_index` | index | Temperature increase index -- flags abnormal temperature rises (e.g. fever) **(derived)** |
| `temp_limit_crossing` | count | Number of times temperature crossed a predefined threshold **(derived)** |
| `temp_group_ratio_svm_inc_index` | index | SVM-based classifier for temperature increase relative to group **(derived)** |
| `temp_svm_inc_index` | index | SVM-based classifier for individual temperature increase **(derived)** |
| `temp_duration_index` | index | Duration of temperature deviation episodes **(derived)** -- present in smaXtec documentation but may not appear in all export files |
| `temp_repetitiveness_index` | index | Repetitiveness of temperature deviation patterns **(derived)** -- present in smaXtec documentation but may not appear in all export files |


## smaxtec_water_intake

Water intake is not measured directly but inferred from rumen temperature
drops caused by cold water entering the reticulum during drinking events
**(confirmed)**.

| Column | Unit | Description |
|---|---|---|
| `timestamp` | ISO-8601 datetime | Date/time of the daily water intake estimate |
| `water_intake_liter` | litres/day | Estimated daily water intake (und Wasser, l/Tag) **(confirmed)** |


## dwd_weather

Daily weather summaries from the Deutscher Wetterdienst (DWD).

| Column | Unit | Description |
|---|---|---|
| `dt` | date (YYYY-MM-DD) | Date of the weather observation |
| `thi_max` | dimensionless index | Maximum Temperature-Humidity Index for the day |
| `qb_thi` | categorical code | DWD quality flag (Qualitaetsbyte) **(?)** |
| `num_values_thi` | count | Number of hourly measurements used for daily THI |
| `enthalpy_max` | kJ/kg | Maximum air enthalpy for the day (heat load metric) |
| `qb_enthalpy` | categorical code | DWD quality flag for enthalpy **(?)** |
| `num_values_enthalpy` | count | Number of hourly measurements used for daily enthalpy |


## hobo_weather

Weather station readings from HOBO loggers.  Part of the barn climate
monitoring system (Stallklima) alongside LoRaWAN and smaXtec barn sensors
**(confirmed)**.  Measurements include air temperature at various positions
with barn mean (Lufttemperatur: versch. Positionen und Mittelwert von
Stall-Neubau-Sensoren), relative humidity, solar radiation (Solarstrahlung),
wind speed (Windgeschwindigkeit), and wind direction **(confirmed)**.

Column names include HOBO logger serial numbers from the export format.

| Column | Unit | Description |
|---|---|---|
| `datetime` | ISO-8601 datetime | Date/time of weather station reading |
| `21136553_b_battery_v` | Volt | Battery voltage of HOBO logger 21136553 |
| `21141733_1_temperature` | deg C | Ambient air temperature (Lufttemperatur) |
| `21141733_2_rh` | % | Relative humidity (rel. Luftfeuchtigkeit) |
| `21141733_3_dew_point` | deg C | Dew point temperature |
| `21141733_b_battery_level` | % | Battery level of logger 21141733 |
| `21141735_1_solar_radiation` | W/m2 | Incoming solar radiation (Solarstrahlung) |
| `21141735_b_battery_level` | % | Battery level of logger 21141735 |
| `21141737_1_wetness` | arbitrary units | Leaf/surface wetness sensor reading |
| `21141737_b_battery_level` | % | Battery level of logger 21141737 |
| `21141734_1_wind_speed` | m/s | Wind speed (Windgeschwindigkeit) |
| `21141734_2_gust_speed` | m/s | Peak wind gust speed |
| `21141734_3_wind_direction` | degrees (0-360) | Wind direction (Windrichtung) |
| `21141734_b_battery_level` | % | Battery level of logger 21141734 |


## Notes from collaborators

The following additional parameters were listed in the smaXtec parameter
documentation as desired for analysis but are not present as dedicated
columns in the current export:

- **Laktationstag** (days in milk) -- can be derived from calving dates in `smaxtec_events` (event_type = `calving_confirmation`) and `herdeplus_calving_lactation`
- **Trächtigkeitstag** (pregnancy day, if pregnant) -- can be derived from insemination dates in `smaxtec_events` (event_type = `insemination`) combined with `pregnancy_result` events
- **Lux** (illuminance) -- listed in barn climate requirements but not present in HOBO export; may require additional sensor

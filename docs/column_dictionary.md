# Column dictionary

Descriptions of every column across all tables in the DigiMuh database.
Columns marked with **(?)** are best guesses pending confirmation from
the smaXtec or HerdePlus documentation.


## allocations

| Column | Unit | Description |
|---|---|---|
| `animal_id` | integer ID | EU ear tag number (15-digit unique identifier) |
| `datetime_enter` | ISO-8601 datetime | Date/time the animal entered the group/pen |
| `datetime_exit` | ISO-8601 datetime | Date/time the animal left the group/pen |
| `group` | integer code | Group or pen identifier the animal was allocated to |


## diseases

| Column | Unit | Description |
|---|---|---|
| `animal_id` | integer ID | EU ear tag number (same as animal_id in other tables) |
| `disease_first_day` | date | Date the disease or health event was first recorded |
| `disease_stop_day` | date | Date the disease or health event resolved |
| `disease_description` | free text | Textual description of the diagnosis (e.g. BVD AG neg, mastitis, ketosis) |
| `herdeplus_tier_id` | integer | HerdePlus internal animal ID (software-specific identifier) |
| `herdeplus_gesund_id` | integer | HerdePlus internal health-record ID **(?)** |
| `disease_category` | categorical text | Broad classification (e.g. No Disease, Healthy/remarks, Other/Illness) |


## herdeplus

| Column | Unit | Description |
|---|---|---|
| `timestamp` | ISO-8601 datetime | Date/time of the milking event or MLP sample |
| `herdeplus_milked_duration_sec` | seconds | Duration of the milking event |
| `herdeplus_milked_milk_flow` | kg/min | Average milk flow rate during milking |
| `herdeplus_milked_mkg` | kg | Milk yield per milking event (Milchkilogramm) |
| `herdeplus_mlp_mkg` | kg | Milk yield at MLP test day (Milchleistungsprüfung = official milk recording) |
| `herdeplus_mlp_fat_percent` | % | Milk fat content from MLP test-day sample |
| `herdeplus_mlp_fkg` | kg | Fat yield from MLP test day (Fettkilogramm) |
| `herdeplus_mlp_protein_percent` | % | Milk protein content from MLP test-day sample |
| `herdeplus_mlp_ekg_percent` | % | Protein yield percentage from MLP (Eiweisskilogramm as fraction of milk) **(?)** |
| `herdeplus_mlp_lactose` | % | Lactose content from MLP test-day sample |
| `herdeplus_mlp_cell_count` | cells/mL (×1000) | Somatic cell count from MLP (indicator of udder health / mastitis) |
| `herdeplus_mlp_urea` | mg/dL or mg/L | Milk urea nitrogen from MLP (indicator of dietary protein–energy balance) |
| `herdeplus_mlp_f_e` | dimensionless ratio | Fat-to-protein ratio (Fett-Eiweiss-Quotient); values > 1.5 suggest energy deficit / ketosis risk |
| `herdeplus_mlp_lkg` | — | *Unknown — possibly Liter-Kilogramm or lactose-kg; needs confirmation* |
| `herdeplus_mlp_ecm` | kg | Energy-corrected milk yield, standardised to 4.0% fat and 3.4% protein |
| `herdeplus_calving_lactation` | integer | Lactation number (which calving the current lactation follows; 1 = primiparous) |


## bcs

| Column | Unit | Description |
|---|---|---|
| `timestamp` | ISO-8601 datetime | Date of BCS assessment |
| `bcs_wert` | score (1–5) | Body Condition Score (Wert = value), typically 1–5 scale in 0.25 increments |


## gouna

| Column | Unit | Description |
|---|---|---|
| `timestamp` | ISO-8601 datetime | Date/time of respiration measurement |
| `respirationfrequency` | breaths/min | Respiration rate |


## lorawan

| Column | Unit | Description |
|---|---|---|
| `timestamp` | ISO-8601 datetime | Date/time of LoRaWAN sensor reading |
| `battery_level` | % | Battery charge level of the LoRaWAN sensor node |
| `current_ampere` | Ampere | Electrical current drawn by the sensor |


## smaxtec_barns

| Column | Unit | Description |
|---|---|---|
| `timestamp` | ISO-8601 datetime | Date/time of barn climate measurement |
| `rawtemp` | °C | Raw (uncalibrated) temperature reading from barn sensor |
| `rawhum` | % | Raw (uncalibrated) relative humidity from barn sensor |
| `temp` | °C | Processed/calibrated barn temperature |
| `hum` | % | Processed/calibrated barn relative humidity |
| `temp_hum_index` | dimensionless index | Temperature-Humidity Index (THI); values > 68 indicate mild heat stress |


## smaxtec_events

| Column | Unit | Description |
|---|---|---|
| `timestamp` | ISO-8601 datetime | Date/time of the recorded event |
| `cow` | integer ID | EU ear tag number (redundant with animal_id from filename) |
| `event_type` | categorical text | Type of event (e.g. calving_confirmation, insemination, pregnancy_result) |
| `value` | text or numeric | Associated value or outcome for the event (may be empty) |


## smaxtec_derived

Derived metrics computed by smaXtec software from raw rumen bolus measurements.
The bolus sits in the reticulum and measures temperature, pH, and
three-axis acceleration; the derived values are computed from these raw signals.

### Activity

| Column | Unit | Description |
|---|---|---|
| `act` | arbitrary units | Activity level from rumen bolus accelerometer |
| `act_decrease_index` | index | Quantifies activity decrease relative to baseline |
| `act_estrus_index` | index | Activity-based estrus (heat) detection index |
| `act_estrus_preprocess` | arbitrary units | Pre-processed/filtered activity signal for estrus algorithm **(?)** |
| `act_group_heat_index` | index | Group-level heat detection; compares animal's activity to herd **(?)** |
| `act_group_ratio` | dimensionless ratio | Ratio of this animal's activity to group mean |
| `act_index` | index | Normalised general activity index |
| `act_pasture_index` | index | Activity index adjusted for pasture context **(?)** |

### Reproduction and calving

| Column | Unit | Description |
|---|---|---|
| `calving_index` | index | Predicts imminent calving (based on temperature drop and activity changes) |
| `heat_index` | index | Combined heat/estrus detection score (may fuse activity and temperature) |

### Drinking

| Column | Unit | Description |
|---|---|---|
| `drink_cycles_v2` | count | Number of detected drinking bouts (v2 detection algorithm) |

### Rumen state

| Column | Unit | Description |
|---|---|---|
| `in_reticulum` | boolean or 0/1 | Whether the bolus is confirmed in the reticulum **(?)** |
| `ph` | pH units | Rumen pH |
| `ph_under_58` | minutes or count | Duration/count of pH readings below 5.8 (SARA threshold) |
| `rum_classification` | categorical | Rumination state classification (ruminating / not ruminating) |
| `rum_dec_index` | index | Rumination decrease index — flags drops relative to baseline |
| `rum_index` | index | Overall rumination activity index |

### Motility

| Column | Unit | Description |
|---|---|---|
| `mot_period` | seconds | Reticulorumen contraction interval (motility period) |
| `mot_pulse_width` | seconds | Duration of a single reticulorumen contraction pulse |
| `mot_pulse_width_median` | seconds | Median contraction pulse width over a rolling window |
| `mot_period_rum_6h_we` | seconds | Mean contraction interval during rumination, 6-hour window **(?)** |
| `mot_period_not_rum_6h_we` | seconds | Mean contraction interval outside rumination, 6-hour window **(?)** |

### Temperature

| Column | Unit | Description |
|---|---|---|
| `temp` | °C | Core body (rumen) temperature |
| `temp_dec_index` | index | Temperature decrease index — flags abnormal drops |
| `temp_inc_index` | index | Temperature increase index — flags abnormal rises (e.g. fever) |
| `temp_limit_crossing` | count or boolean | Temperature crossing a predefined threshold **(?)** |
| `temp_normal_index` | index | How close temperature is to expected normal range **(?)** |
| `temp_without_drink_cycles` | °C | Rumen temperature with drinking-event artifacts filtered out |
| `temp_group_ratio_svm_inc_index` | — | *Unknown — SVM-based classifier output involving temperature and group ratios* |
| `temp_height_index` | — | *Unknown — possibly temperature amplitude or peak height* |
| `temp_svm_inc_index` | — | *Unknown — SVM-based temperature increase classifier* |


## smaxtec_water_intake

| Column | Unit | Description |
|---|---|---|
| `timestamp` | ISO-8601 datetime | Date/time of the daily water intake estimate |
| `water_intake_liter` | litres | Estimated daily water intake (derived from rumen temperature drops) |


## dwd_weather

| Column | Unit | Description |
|---|---|---|
| `dt` | date (YYYY-MM-DD) | Date of the weather observation |
| `thi_max` | dimensionless index | Maximum Temperature-Humidity Index for the day |
| `qb_thi` | categorical code | DWD quality flag for THI (Qualitätsbyte) **(?)** |
| `num_values_thi` | count | Number of hourly measurements used for daily THI |
| `enthalpy_max` | kJ/kg | Maximum air enthalpy for the day (heat load metric) |
| `qb_enthalpy` | categorical code | DWD quality flag for enthalpy **(?)** |
| `num_values_enthalpy` | count | Number of hourly measurements used for daily enthalpy |


## hobo_weather

Column names include HOBO logger serial numbers from the export format.

| Column | Unit | Description |
|---|---|---|
| `datetime` | ISO-8601 datetime | Date/time of weather station reading |
| `21136553_b_battery_v` | Volt | Battery voltage of HOBO logger 21136553 |
| `21141733_1_temperature` | °C | Ambient air temperature |
| `21141733_2_rh` | % | Relative humidity |
| `21141733_3_dew_point` | °C | Dew point temperature |
| `21141733_b_battery_level` | % | Battery level of logger 21141733 |
| `21141735_1_solar_radiation` | W/m² | Incoming solar radiation |
| `21141735_b_battery_level` | % | Battery level of logger 21141735 |
| `21141737_1_wetness` | arbitrary units | Leaf/surface wetness sensor reading |
| `21141737_b_battery_level` | % | Battery level of logger 21141737 |
| `21141734_1_wind_speed` | m/s | Wind speed |
| `21141734_2_gust_speed` | m/s | Peak wind gust speed |
| `21141734_3_wind_direction` | degrees (0–360) | Wind direction |
| `21141734_b_battery_level` | % | Battery level of logger 21141734 |

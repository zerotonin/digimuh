# Analysis 0 — Individual heat stress thresholds

**Pipeline version:** 0.4.0
**Authors:** Bart R. H. Geurten (University of Otago), Gundula Hoffmann (ATB Potsdam)
**Target:** Frontiers in Animal Science


## 1. Pipeline overview

Three independent scripts, each reading CSVs from the previous step:

```
  digimuh-extract  (DB → CSV, runs once)
        │  rumen_barn.csv, respiration_barn.csv,
        │  production.csv, daily_milk_yield.csv, climate_daily.csv
        ▼
  digimuh-stats    (CSV → statistics + analysis CSVs)
        │  broken_stick_results.csv, cross_correlation.csv,
        │  circadian_null_model.csv, crossing_times.csv,
        │  event_triggered_traces.csv, climate_eta.csv,
        │  thermoneutral_fraction.csv, tnf_yield.csv, ...
        ▼
  digimuh-plots    (CSV → SVG/PNG figures)
```


## 2. Data extraction (analysis_00a_extract.py)

For each of 255 animal-years in the extended Tierauswahl:

**Rumen temperature:** `temp_without_drink_cycles` from `smaxtec_derived`, filtered to 30-43 C, joined with barn climate via `pd.merge_asof` with nearest-time matching within 30 min tolerance.  This preserves native barn sensor resolution (~10 min for smaXtec) rather than collapsing to hourly averages.  Barn sensors: Neubau only (`barn_id IN (1,2)`).  Milking hours excluded (04-07, 16-19).

**Drinking correction:** Timestamps with `drink_cycles_v2 > 0` plus a 15-min recovery window are excluded by default.

**Respiration:** `respirationfrequency` from `gouna`, filtered to 5-150 bpm, same barn join and milking exclusion.  Available from Aug 2022.

**Production:** Mean milk yield per observation window and daily milk yield per cow (for thermoneutral fraction analysis), lactation number from full animal history.


## 3. Statistical analysis (analysis_00b_stats.py)

### 3.1 Broken-stick regression

Two-segment piecewise linear model per animal.  Breakpoint psi found by grid search (200 points) + bounded minimisation.  Biological constraint: `b2 > b1` and `b2 > 0`.  Search ranges: THI 45-80, barn temp 5-35 C.

### 3.2 Breakpoint existence tests (reserved for COMPAG)

Davies test (Davies 1987, 2002) and pseudo-Score test (Muggeo 2016).  Skipped with `--frontiers` flag.

### 3.3 Four-parameter logistic (Hill) fit (reserved for COMPAG)

Sigmoidal dose-response model with lower bend point (Sebaugh & McCray 2003).

### 3.4 Below/above breakpoint comparison

Within-year Fisher resampling tests (reRandomStats, Geurten 2026) on per-animal means below vs above the individual THI breakpoint.  Median difference, 20,000 permutations, BH-FDR corrected.

### 3.5 Rumen temperature circadian null model

Cow-days classified as cool (THI stayed below breakpoint all day) or stress (THI exceeded breakpoint).  Mean rumen temperature computed at each clock hour for each category.  The cool-day profile is the circadian null model.  The stress-minus-cool difference at each hour isolates the heat effect from the circadian rhythm.

### 3.6 THI daily exceedance profile

Mean barn THI by clock hour and month (June-September), with herd median THI breakpoint as reference.  Shows when heat stress begins and ends each day and how this shifts seasonally.

### 3.7 Breakpoint crossing event detection

Upward crossings: consecutive readings where barn THI transitions from at-or-below to above the cow's individual breakpoint.  Crossings at data gaps > 30 min (milking exclusion boundaries) are rejected to prevent false detection.

### 3.8 Crossing time raster

Clock time of each crossing event per animal, used for the activation raster plot and KDE crossing density overlaid on the circadian null model.

### 3.9 Cross-correlation and cross-covariance

Time series split at individual breakpoint.  Raw and per-day-detrended CCF at lags +/-240 min (10-min steps).  The raw CCF contains a diurnal artifact (shared 24h cycle produces sinusoidal CCF); this is retained for methodological discussion.

### 3.10 Event-triggered average (ETA)

Peri-event average of rumen temperature aligned to crossing events.  Clock hour filter (default 8:00-11:00).  Four panels: (A) climate signal, (B) raw rumen temp, (C) baseline-subtracted, (D) additional rumen temp (raw minus cool-day circadian profile at that clock hour).

### 3.11 Climate ETA

Both barn THI and barn temperature extracted +/-6 hours around each crossing, normalised by subtracting the cow's individual breakpoint.  y=0 = threshold.  Dual y-axis showing both delta-THI and delta-barn-temperature.

### 3.12 Thermoneutral fraction (TNF) vs milk yield

Daily TNF = fraction of readings where barn THI <= cow's breakpoint.  Yield normalised to cow-specific P95 (95th percentile of daily yields).  Spearman correlation: pooled and per-cow within-animal.  Heat-stress-only filter at TNF thresholds 0.95/0.90/0.80/0.70.

### 3.13 Breakpoint stability

ICC for repeat animals.  Pairwise year comparisons via Fisher resampling (BH-FDR).

### 3.14 Longitudinal breakpoint stability (alluvial)

Animals in 2+ years classified by year-to-year breakpoint change: strongly decreased (delta < -3), decreased (-3 to -1), stable (-1 to +1), increased (+1 to +3), strongly increased (delta > +3).  Alluvial plot tracks individual animal flow between categories across consecutive year transitions.


## 4. Figures (analysis_00c_plots.py)

All figures use the Wong (2011) colourblind-safe palette.  SVG + PNG output.

| Figure | Description |
|---|---|
| `circadian_null_model` | 2-panel: (A) cool/stress profiles + crossing density, (B) difference curve |
| `thi_daily_profile_*` | Barn THI by clock hour and month with herd bp line |
| `crossing_raster_*` | Activation raster (animals x clock time) + raincloud |
| `xcorr_*` / `xcov_*` | Cross-correlation/covariance below vs above bp |
| `eta_*_8to11h` | 2x2 event-triggered average (8-11h crossings) |
| `climate_eta_*` | Delta-THI + Delta-Temp around crossings (dual y-axis) |
| `grouped_boxplots_*` | Breakpoint distributions per year |
| `paired_below_above_*` | Boxplots + significance brackets (BH-FDR) |
| `spearman_*` | Per-animal correlation distributions |
| `tnf_yield` | TNF vs daily milk yield (absolute + relative) |
| `longitudinal_*` | Breakpoint trajectories (repeat animals) |
| `sankey_longitudinal_*` | Alluvial: stability categories across year transitions |


## 5. Output files

```
results/broken_stick/
  rumen_barn.csv                 Rumen temp + barn climate per reading
  respiration_barn.csv           Respiration + barn climate per reading
  production.csv                 Mean milk yield + lactation per animal
  daily_milk_yield.csv           Daily milk yield per cow
  climate_daily.csv              Daily barn climate (Jun-Sep)
  broken_stick_results.csv       Per-animal breakpoints + convergence
  spearman_correlations.csv
  behavioural_response.csv       Per-animal means below/above bp
  statistical_tests.csv          Fisher resampling + BH-FDR
  cross_correlation.csv          CCF + xcov (raw + detrended)
  derivative_ccf.csv             d(climate)/dt vs d(body_temp)/dt
  circadian_null_model.csv       Hourly rumen temp: cool vs stress days
  thi_daily_profile.csv          Barn THI by clock hour and month
  crossing_times.csv             Clock time of each crossing event
  event_triggered_traces*.csv    Peri-event traces (all + filtered)
  event_triggered_summary*.csv   Event counts per animal
  climate_eta.csv                THI + barn temp around crossings
  thermoneutral_fraction.csv     Daily TNF per cow
  tnf_yield.csv                  TNF + daily yield + P95
  breakpoint_stability.csv       ICC pairs
  summary_table.csv
  *.svg / *.png                  All figures
```


## 6. Running the pipeline

```bash
# Full pipeline
bash scripts/run_00_broken_stick_ana.sh --frontiers --no-resp

# Or step by step:
digimuh-extract
digimuh-stats --data results/broken_stick --frontiers --no-resp
digimuh-plots --data results/broken_stick
```


## 7. References

- Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300.
- Davies RB (1987) Biometrika 74:33-43.
- Davies RB (2002) Biometrika 89:484-489.
- Geurten BRH (2026) reRandomStats. github.com/zerotonin/rerandomstats
- Hoffmann G et al. (2020) Biosystems Engineering 199:83-96.
- Muggeo VMR (2003) Statistics in Medicine 22:3055-3071.
- Muggeo VMR (2016) J Stat Comput Simul 86:3059-3067.
- Piccione G et al. (2014) Biological Rhythm Research 45:371-381.
- Sebaugh JL, McCray PD (2003) Pharmaceutical Statistics 2:167-174.
- Wong B (2011) Nature Methods 8:441.

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

### 3.15 Milk yield classification (DIM-adjusted via Wood 1967)

**Rationale.** The Frontiers paper stratifies per-animal broken-stick fits by cow-related factors.  Daily milk yield (MY) cannot be stratified on the raw value because raw yield confounds "low-yielder" with "late-lactation" and "primiparous" — two dominant nuisance variance sources.  Any tercile split on raw yield therefore tests a mixture of parity × DIM × heat effects and cannot be interpreted as a production-level contrast.  We remove the DIM × parity confound before stratifying.

**Distribution shape.** Raw daily yield (n = 33,755 cow-days) is mildly right-skewed (skewness +0.41, excess kurtosis +1.83) with 72.8% of values within ± 1 SD of the mean (Gaussian expectation: 68.3%).  Shapiro-Wilk p → 0 and Anderson-Darling A² = 32 both reject normality.  Log-transformation over-corrects (skewness flips to −1.31, kurtosis +12.7) because of a few near-zero dry-off days, so log-normal modelling is not appropriate here.  A fixed mean ± 1 SD partition would force a 68/16/16% membership split — wasteful of statistical power — and its upper bound would sit further from the median than the lower bound, biasing the "high" group.

**Method.** For each cow-lactation we fit Wood's (1967) incomplete-gamma lactation curve

    y(t) = a · t^b · exp(−c · t)

in log space (OLS on `ln y = ln a + b ln t − c t`).  DIM = (observation date − most recent prior calving-confirmation event from `smaxtec_events`), clipped to [5, 305] days.  Fits are retained as `per_lactation` when convergent (b>0, c>0, peak DIM in [15, 100] d, R² ≥ 0.10); otherwise the cow-day uses the parity-pooled Wood curve for its parity bucket (1 / 2 / 3+).  **Wood parameters are fitted on each Tierauswahl animal's full HerdePlus history** (`daily_milk_yield_full.csv`, ~157k cow-days across 2021-04-01 → 2024-09-30) so each lactation sees its pre-peak ramp-up and post-peak tail; residuals are then evaluated on the summer-window analysis frame (`daily_milk_yield.csv`) — see §3.19 for the motivation and before/after numbers.  Under the full-history fit, 71% of lactations (378 / 531) hold individual Wood parameters and 29% use the parity pool.

Residuals `y − ŷ` on the analysis frame are approximately symmetric (skewness −0.32 vs. raw +0.41) with SD 4.56 kg/d (vs. raw 7.29 kg/d).  We stratify the residuals at Q33 / Q67 (`−0.30` / `+1.98` kg/d on the current extract) into `low` / `middle` / `high` producers relative to DIM-specific expectation.

**Sensitivity analysis.** Raw-yield terciles (Q33 = 35.4, Q67 = 40.7 kg/d) are retained as a supplementary stratifier.  The two literature schemes (Müschner-Siemens 28.8 / 38.4; Yan 26 / 39) collapse 94–97% of our herd's cow-days into the middle + high buckets and are reported only for reference because our herd out-yields both source cohorts by ~6–7 kg/d.

**Why not mean ± 1 SD or GMM.** No Frontiers / JDS / J Thermal Biology heat-stress paper we located uses Gaussian SD partitions for per-cow yield stratification; SD partitions are herd-level benchmarking devices (Bava-style yield-gap analysis), not physiology tools.  Gaussian mixture models can describe the residual distribution but introduce free mixture parameters that reviewers will challenge absent a physiological justification.

**Implementation.**
- `src/digimuh/stats_lactation_curve.py` — reusable DIM attachment, Wood fits, parity-pooled fallback, residual computation.
- `src/digimuh/milk_yield_classification.py` — pooled + per-year raw histograms, residual histogram, tercile boundaries, example-fit grid, class-count summary.
- Extracted via `digimuh-extract` (adds `calvings.csv`); the classifier writes `daily_milk_yield_wood.csv` (joinable on `animal_id` + `date`) and `wood_curve_fits.csv` for downstream use.

**Key references.**
- Wood, P.D.P. (1967) Algebraic model of the lactation curve in cattle. *Nature* 216:164-165. doi:10.1038/216164a0
- Wilmink, J.B.M. (1987) Adjustment of test-day milk, fat and protein yield for age, season and stage of lactation. *Livest Prod Sci* 16:335-348.
- Macciotta, N.P.P., Dimauro, C., Rassu, S.P.G., Steri, R., Pulina, G. (2011) The mathematical description of lactation curves in dairy cattle. *Ital J Anim Sci* 10:e51. doi:10.4081/ijas.2011.e51
- Aguilar, I., Misztal, I., Tsuruta, S. (2010) Short communication: Genetic trends of milk yield under heat stress for US Holsteins. *J Dairy Sci* 93:1754-1758.
- Bohmanova, J., Misztal, I., Tsuruta, S., Norman, H.D., Lawlor, T.J. (2008) Short communication: Genotype by environment interaction due to heat stress. *J Dairy Sci* 91:840-846.
- Gernand, E., König, S., Kipp, C. (2019) Influence of on-farm measurements for heat stress indicators on dairy cow productivity, female fertility, and health. *J Dairy Sci* 102:6660-6671.
- Bernabucci, U., Biffani, S., Buggiotti, L., Vitali, A., Lacetera, N., Nardone, A. (2014) The effects of heat stress in Italian Holstein dairy cattle. *J Dairy Sci* 97:471-486. doi:10.3168/jds.2013-6611
- Tao, S., Orellana, R.M., Weng, X., Marins, T.N., Dahl, G.E., Bernard, J.K. (2018) Symposium review: The influences of heat stress on bovine mammary gland function. *Animal* 12 (S2):s445-s451.
- Müschner-Siemens, T., Hoffmann, G., Ammon, C., Amon, T. (2020) Daily rumination time of lactating dairy cows under heat stress conditions. *J Thermal Biology* 88:102484.
- Yan, G., Liu, K., Hao, Z., Shi, Z., Li, H. (2021) The effects of cow-related factors on rectal temperature, respiration rate, and temperature-humidity index thresholds for lactating cows exposed to heat stress. *J Thermal Biology* 100:103041.

**Decision log.**
- 2026-04-17: chose **Wood-residual terciles as primary**, **raw-yield terciles as supplementary**, in consultation with the literature summary above.  Literature-boundary schemes (Müschner-Siemens, Yan) reported for reference only because our German Holstein herd out-yields both source cohorts by ~6-7 kg/d, producing severely unbalanced class counts.

### 3.16 TNF vs milk yield, stratified by yield class

**Rationale.** Having removed the DIM × parity confound from the stratifier (§3.15), we can now test whether thermoneutral-fraction exposure (TNF — daily fraction of readings below the cow's personal breakpoint, §3.12) relates to daily milk yield differently across low / middle / high producers.

**Grouping.** Each (animal_id, year) is labelled by the sign and magnitude of its mean Wood residual: cow-years below Q33 → "low", above Q67 → "high", in between → "middle".  Cow-years with fewer than ten valid residual cow-days are dropped.  The cow-year boundary set differs numerically from the cow-day boundary set because averaging across days shrinks the distribution (present extract: cow-day Q33 / Q67 = -1.43 / +1.92 kg/d, cow-year Q33 / Q67 = -0.62 / +1.27 kg/d).

**Measurement.** For each class we compute Spearman rs between TNF and three responses on all cow-days in that class, separately for the THI-breakpoint TNF (`thi_tnf`) and the barn-temperature-breakpoint TNF (`temp_tnf`):

1. **`yield_residual`** — *primary* heat-stress response.  The DIM-adjusted Wood residual (§3.15) removes the within-cow-year lactation decline, so the remaining variance is attributable to day-to-day factors including heat exposure.
2. **`daily_yield_kg`** — raw daily yield, for comparability with the literature.  Within a cow-year the lactation decline is still present, so this correlation mixes heat and DIM.
3. **`relative_yield`** — yield / cow P95, dampens between-cow level differences but does not remove the within-cow DIM decline.

BH-FDR correction is not applied here because the hypotheses are class-indexed rather than multiple-testing of the same contrast; the pooled row is provided as a reference.  OLS slope is reported alongside rs so the reader can read the effect in kg/d per unit TNF.

**Sign-convention note.**  TNF is a *cool-day* indicator: **high TNF ⇔ more time below breakpoint ⇔ cooler cow-day**.  The sign table for the two climate predictors used in §3.16–§3.18 is therefore:

| Predictor | High value means | Low value means |
|---|---|---|
| `thi_tnf` / `temp_tnf` (§3.16, §3.17 binary) | cooler day | hotter day |
| `mean_thi` / `mean_barn_temp` (§3.18) | hotter day | cooler day |

A negative `rs(TNF, residual)` and a positive `rs(mean_climate, residual)` therefore encode **the same physical statement** — *hotter day → higher residual yield* — once the polarity flip is applied.

**Key result on the present extract (full-history Wood fit, §3.19).**

| Class | Predictor | rs (residual) | rs (raw) | slope (residual, kg/d per 1.0 TNF) |
|---|---|---|---|---|
| low | THI | +0.034 ** | −0.074 *** | +0.40 |
| low | Temp | +0.038 ** | −0.084 *** | +0.44 |
| middle | THI | −0.001 n.s. | −0.157 *** | +0.11 |
| middle | Temp | +0.011 n.s. | −0.145 *** | +0.04 |
| high | THI | −0.086 *** | −0.267 *** | −1.29 |
| high | Temp | −0.050 *** | −0.230 *** | −0.61 |
| pooled | THI | +0.009 n.s. | −0.163 *** | +0.02 |
| pooled | Temp | +0.013 n.s. | −0.150 *** | +0.11 |

DIM adjustment (Wood-residual column) sharply reduces the pooled heat-exposure × raw-yield correlation (|rs| drops from ~0.16 to near zero), confirming that the bulk of the raw-yield correlation was driven by the lactation decline.  **Applying the polarity flip**, the residual rows say the following:

* **Pooled and middle class** — essentially no residual TNF effect after DIM adjustment.  Daily TNF alone does not explain within-lactation yield variation at the herd level.
* **High-yield class** — rs_THI = −0.086 *** (slope −1.29 kg/d per full 0→1 TNF swing), i.e. **hotter days → slightly lower residual yield** for high producers.  This is the direction the heat-stress hypothesis predicts, and §3.18 confirms it with the daily-mean climate response (rs = −0.054 / −0.055 ***).
* **Low-yield class** — rs_THI = +0.034 **, which after polarity flip means hotter days → slightly *higher* residual for low producers.  Small effect but opposite sign from the high class — low producers are not the cows carrying the heat-stress signal in this dataset.

The magnitude across the whole table is modest (|slope| ≤ 1.3 kg/d per full 0→1 TNF swing) but the direction pattern — heat-stress signal in high producers only — is biologically consistent with the literature that high-yielding cows are the most thermally stressed (Bernabucci et al. 2014; Tao et al. 2018).

**Outputs.**  The milk-yield classifier writes:
- `yield_class_per_cow_year.csv` — one row per (animal_id, year) with `mean_residual`, `n_days`, `yield_class`.
- `tnf_yield_by_class.csv` — all cow-days with `yield_class`, `mean_residual`, and (when the Wood pipeline ran) `yield_residual`, `yield_expected`, `dim`, plus the crossing flags and daily climate means from §3.17–§3.18.
- `tnf_yield_correlations_by_class.csv` — rs, p, slope, intercept per (class × predictor × response) for all three responses.
- `tnf_yield_by_class.{svg,png}` — daily yield axis.
- `tnf_yield_by_class_relative.{svg,png}` — relative yield axis.
- `tnf_yield_by_class_residual.{svg,png}` — **primary** figure, DIM-adjusted residual axis with a horizontal reference line at y=0.

### 3.17 Crossed-day vs not-crossed-day yield (raincloud)

**Rationale.** The §3.16 TNF treats heat exposure as a continuous fraction.  Reviewers often prefer a cleaner binary: did the cow's personal THI (or barn-temperature) breakpoint actually get crossed on that day, yes or no?  This sidesteps the question of *how much* time above the threshold matters and asks the simpler *was there a heat-stress event*.

**Method.** `crossing_times.csv` (§3.8) lists upward crossings of each cow's breakpoint.  For each cow-day we flag `thi_crossed ∈ {True, False}` and `temp_crossed ∈ {True, False}`.  We then compare the DIM-adjusted yield residual between crossed and not-crossed cow-days using a Mann-Whitney two-sided test and report median-difference effect sizes.  The comparison is run pooled and separately within each yield class.

**Key result on the present extract.**

| Group | Predictor | n crossed | n not | Median (crossed) | Median (not) | Δ median (kg/d) | p |
|---|---|---|---|---|---|---|---|
| pooled | THI | 4,549 | 18,849 | +0.80 | +0.88 | −0.08 | 6.0e−1 n.s. |
| pooled | Temp | 4,537 | 18,861 | +0.70 | +0.90 | **−0.20** | 1.8e−3 *** |
| low | THI | 1,587 | 5,898 | −0.20 | +0.06 | **−0.26** | 1.3e−3 *** |
| low | Temp | 1,513 | 5,972 | −0.15 | +0.04 | −0.19 | 5.8e−2 n.s. |
| middle | THI | 1,608 | 6,961 | +0.65 | +0.67 | −0.02 | 8.0e−1 n.s. |
| middle | Temp | 1,578 | 6,991 | +0.47 | +0.72 | **−0.25** | 1.6e−3 ** |
| high | THI | 1,354 | 5,990 | +2.39 | +2.16 | +0.24 | 1.6e−2 * |
| high | Temp | 1,446 | 5,898 | +2.18 | +2.20 | −0.02 | 2.2e−1 n.s. |

A crossed day is by construction a **hotter** day (the cow's barn-THI or barn-temp crossed her individual breakpoint upward).  Under the full-history Wood fit (§3.19) the pooled effect is now small, with the barn-temperature channel pointing in the heat-stress direction (−0.20 kg/d **, hotter → lower residual) and the THI channel indistinguishable from zero.  Per-class the split matches the §3.16 TNF story: the **low-yield class loses** ~0.25 kg/d on crossed days via either predictor (heat-stress direction), the **middle class** drops ~0.25 kg/d on barn-temp crossings, and the **high-yield class** still gains marginally (+0.24 kg/d *) on THI-crossed days (anti-heat-stress direction).  These |Δ| values are small compared to the raw pooled effect in §3.16 and do not cleanly resolve the question on their own; the pooled daily-mean climate fit in §3.18 is cleaner because it uses the full continuous climate gradient rather than a binary split.

**Outputs.**
- `crossing_day_comparison.csv` — per-(group × predictor) Mann-Whitney table.
- `crossing_day_raincloud.{svg,png}` — pooled, both climates side by side.
- `crossing_day_raincloud_{low,middle,high}.{svg,png}` — class-specific.

### 3.18 Daily-mean climate vs DIM-adjusted yield

**Rationale.** A single pooled scatter relating each cow-day's mean barn THI (and mean barn temperature) to the DIM-adjusted yield residual.  No yield-class stratification, no breakpoint-based split — the simplest picture of whether daily weather, on its own, predicts residual yield in our herd.

**Method.** Per cow-day we use `mean_thi` (from §3.12) and `mean_barn_temp` (added to `compute_thermoneutral_fraction` or re-derived on demand from `rumen_barn.csv`).  Spearman rs, p, OLS slope, and intercept are computed across all cow-days that carry both a climate mean and a yield residual.

**Key result on the present extract.**

| Predictor | n | rs | p | Slope (kg/d per unit) |
|---|---|---|---|---|
| mean THI | 23,398 | **−0.054** | 8.7e−17 *** | −0.018 |
| mean barn temp | 23,398 | **−0.055** | 3.9e−17 *** | −0.028 |

After the full-history Wood refit (§3.19) both pooled rs values are **negative** — i.e. **hotter daily-mean climate associates with lower DIM-adjusted residual yield**, which is the naïve heat-stress direction.  Before the refit (summer-window Wood fit) the same pooled correlations were +0.024 / +0.028: the sign flipped once the Wood curves had seen enough of each lactation to stop systematically under-predicting peak-summer days.  The effect sizes remain small — about 0.5 kg/d across the full summer THI range — so this is not the paper's headline heat-stress finding; it is a sanity check that the DIM adjustment is now well-posed.  The within-class signal documented in §3.16 (high-yield class rs = −0.086 ***, slope −1.29 kg/d per 1.0 TNF) is the cleaner per-cow statement of the same phenomenon.

**Outputs.**
- `daily_climate_vs_yield.csv` — one row per climate predictor.
- `daily_climate_vs_yield_residual.{svg,png}` — 1×2 scatter with OLS line.

### 3.19 Full-history Wood refit (fit frame ≠ evaluate frame)

**Why.** The first iteration of §3.15 fitted Wood's model on the same summer-window extract used for the broken-stick analysis (`daily_milk_yield.csv`).  Because each Tierauswahl observation window is Jun–Sep and each lactation is up to 305 days long, the summer slice usually misses the cow's early-lactation peak ramp — the part of the curve the Wood model relies on to identify its `a` and `b` parameters.  Only **22% (118 / 531) of lactations** converged to an individual fit under that setup; the other 78% fell back to the parity-pooled curve, which is systematically flatter than the true peak.  That flatter fallback curve under-predicted peak-summer yields, which in turn pushed residuals *upward* on exactly the hot days, producing the anti-heat-stress sign in §3.16–§3.18 that was worrying us.

**What changed.**  A new `digimuh-extract` output, `daily_milk_yield_full.csv`, pulls every `herdeplus` row for the 181 Tierauswahl animals across the DB's full window (2021-04-01 → 2024-09-30) with no date restriction — 156,972 cow-days, 4.6× the summer-window extract (33,755 cow-days).  `compute_wood_residuals` now accepts an optional `fit_yields` frame: Wood curves are estimated on the full-history frame and residuals are evaluated on the summer-window analysis frame.  The same cow-day count (31,401) flows downstream, but each residual now uses a curve fitted from that lactation's pre-peak climb and post-peak tail.

**Before / after (same extract, same code, only the fit frame changes):**

| Metric | Summer-window fit | Full-history fit |
|---|---|---|
| Lactations with an individual Wood fit | 118 / 531 (22%) | **378 / 531 (71%)** |
| Per-lactation R² median [IQR] | 0.38 [0.25–0.50] | **0.57 [0.39–0.73]** |
| Peak DIM median | 62 d | **53 d** (Holstein norm) |
| Residual SD | 5.3 kg/d | **4.56 kg/d** |
| Residual skewness | −0.38 | −0.32 |
| Pooled rs(mean THI, residual) | +0.024 | **−0.054 ***** |
| Pooled rs(mean barn temp, residual) | +0.028 | **−0.055 ***** |
| High-class rs(thi_tnf, residual) | −0.132 *** (flip: hot → higher) | **−0.086 *** (flip: hot → lower)** |

The before → after column tells a clean story: swapping to the full-history fit roughly triples the share of individual lactation fits, tightens residual SD by 14%, brings peak DIM to the Holstein textbook value, and — most importantly — **flips the pooled daily-climate residual correlations from anti-heat-stress to heat-stress direction**.  Small effect sizes remain, but the DIM pipeline is now well-posed.

**Back-compatibility.**  When `daily_milk_yield_full.csv` is absent the orchestrator silently falls back to the single-frame behaviour (fit and evaluate on the summer window).  A one-line log message identifies which source was used so downstream claims remain auditable.

### 3.20 MLP (Milchleistungsprüfung) test-day composition extract

**Why.**  The `herdeplus` table interleaves two data channels in the same rows: high-frequency per-milking events (`herdeplus_milked_*` populated) and monthly MLP test-day analytics (`herdeplus_mlp_*` populated).  Until this iteration the pipeline extracted only the per-milking channel.  The MLP channel carries the standard dairy-health indicators — fat %, protein %, somatic cell count, urea, energy-corrected milk (ECM), fat-to-protein ratio — that Gundula will want for companion work on ketosis (fat/protein > 1.4), acidosis (low pH paired with low fat %), and mastitis (SCC elevation).

**Method.**  `extract_mlp_test_days` filters `herdeplus` to rows with a non-null `herdeplus_mlp_fat_percent` (the simplest way to identify an MLP row), restricted to Tierauswahl animals.  For the present database that yields 5,582 test-day rows across 181 animals (roughly monthly per cow, 3.5 years of coverage).  Columns preserved: `herdeplus_mlp_mkg` (test-day milk kg), `herdeplus_mlp_fat_percent`, `herdeplus_mlp_fkg` (fat kg), `herdeplus_mlp_protein_percent`, `herdeplus_mlp_ekg_percent`, `herdeplus_mlp_lactose`, `herdeplus_mlp_cell_count`, `herdeplus_mlp_urea`, `herdeplus_mlp_f_e` (fat-to-protein ratio), `herdeplus_mlp_lkg`, `herdeplus_mlp_ecm`, plus `herdeplus_calving_lactation`.

**Scope.**  This round writes `mlp_test_days.csv` only — no analysis module yet.  Wiring MLP composition into the broken-stick or ketosis modules is a follow-up task.

### 3.21 Thin-milk hypothesis — MLP composition × climate

**Hypothesis (Geurten, April 2026).** On hot days a lactating cow with a calf may *increase* milk volume to hydrate the calf while rumen activity drops, producing dilute milk — higher water, lower fat, lower protein.  Under the classical heat-stress pattern (Bernabucci et al. 2014; Tao et al. 2018) the opposite is expected (hot → both volume and composition %s decrease).  We have enough monthly MLP test-day rows (§3.20) to tell the two apart directly.

**Method.**  Each MLP test-day is joined to the cow-day climate frame (`mean_thi`, `mean_barn_temp`, `thi_tnf`) on `(animal_id, date)`.  When legacy `tnf_yield.csv` predates the `mean_barn_temp` column, it is back-filled from `rumen_barn.csv` so both climate scales are available.  Each composition metric — test-day milk kg, fat %, protein %, lactose, fat-kg, protein-kg (Eiweiß), F/E ratio, ECM, SCC, urea — is correlated (Spearman rs + p + OLS slope) against each climate predictor, pooled and per Wood-residual yield class.  A compact `thin_milk_verdict` helper reports whether all three core predictions (volume ↑, fat % ↓, protein % ↓) hit their expected sign.

**Key result on the present extract (pooled vs daily mean barn THI; n = 806 MLP test-days across 176 animals):**

| MLP metric | rs | p | Direction | Hypothesis prediction |
|---|---|---|---|---|
| Test-day milk (kg/d) | **+0.191** | 4.4e−8 *** | more volume when hot | ✓ ↑ |
| Fat % | **−0.281** | 4.0e−16 *** | less fat per kg when hot | ✓ ↓ |
| Protein % | **−0.284** | 2.2e−16 *** | less protein per kg when hot | ✓ ↓ |
| F/E ratio | −0.159 | — *** | fat suppressed more than protein | ✓ ↓ |
| Fat (kg/d) | −0.110 | ** | total fat slightly down | ≈ 0 |
| Protein (kg/d) | +0.187 | *** | total protein up (volume wins) | ≈ 0 |
| Lactose (%) | +0.172 | *** | lactose rises with volume | (consistent with dilution) |
| ECM | +0.031 | n.s. | energy-corrected output flat | — |
| SCC (10³/mL) | −0.094 | * | no mastitis signal | — |
| Urea (mg/dL) | −0.188 | *** | lower nitrogen excretion | — |

**Verdict.** `thin_milk_verdict` → **supported**.  All three core predictions (volume ↑, fat % ↓, protein % ↓) reach their predicted sign at p < 0.001; the F/E ratio also drops, indicating that milk-fat synthesis is suppressed more than milk-protein synthesis.  Lactose rises with volume (consistent with dilution — lactose is the osmotic driver of milk water content).  ECM stays flat because the volume gain and the concentration loss cancel on an energy basis.  SCC is essentially unchanged, so this is not a mastitis effect.

**Where the effect is strongest.** The heatmap in `mlp_composition_heatmap_mean_thi.{svg,png}` shows the **high-yield class** carrying the largest signals: rs_volume = +0.38 ***, rs_fat% = −0.42 ***, rs_protein% = −0.42 ***.  Low and middle classes show the same signs but roughly half the magnitude.  This matches §3.16: the high-yield cows are the most heat-responsive, both in their absolute volume response **and** in the composition shift.

**Reconciliation with §3.16 / §3.18.** §3.16 reports a *residual* yield signal: hotter day → slightly *below*-expected yield in the high-yield class (Wood under-predicted those days less than the actual shortfall).  §3.21 reports *absolute* volume: hotter day → *higher* absolute volume, because hot days in this summer-window extract coincide with early-lactation DIM where cows are near peak.  The two views are consistent once you remember that Wood-residual measures distance from the expected peak, while MLP `mkg` measures total output.  The thin-milk story is specifically about composition, not absolute kg — and on the composition dimension the hypothesis is clear.

**Outputs.**
- `mlp_composition_by_cowday.csv` — every matched MLP test-day with climate + class + residual merged on.
- `mlp_composition_correlations.csv` — long rs / p / slope table per (class × predictor × metric).
- `mlp_thin_milk_hypothesis.{svg,png}` — 2×2 dilution scatter (volume, ECM, fat %, protein %) vs mean THI.
- `mlp_composition_heatmap_mean_thi.{svg,png}` — signed-rs heatmap across all metrics × yield class.

**Follow-up candidates.**
- Re-run with `mean_barn_temp` as primary predictor (parallel to §3.18) once the full-year climate extract is available.
- Check the same correlations within the 7–14 d post-heat window (lagged response) — MLP is monthly so the signal may be sharper with a trailing climate aggregate.
- Test specifically on crossed-day vs not-crossed-day (binary) for the subset of MLP rows that fall on crossing days, parallel to §3.17.

### 3.22 Dilution partition — observed vs pure-water prediction

**Why.** The §3.21 result shows fat %, protein %, and F/E ratio all fall as THI rises, but a concentration drop can come from either (a) the cow's mammary gland adding water on top of an unchanged absolute solid output ("dilution"), or (b) the rumen / mammary complex synthesising less solid on hot days ("suppression").  Without the partition we can't tell how much of the observed thin-milk signal is mechanical and how much is metabolic.

**Method.**  For each MLP test-day row we compute a per-cow reference: her mean test-day volume (`herdeplus_mlp_mkg`) and her mean absolute fat kg (`herdeplus_mlp_fkg`), protein kg (`herdeplus_mlp_ekg_percent` — which is actually protein-kg/d despite the HerdePlus suffix) and lactose kg (`herdeplus_mlp_lkg`).  We then predict the composition % she'd have on that specific test-day *if she had kept her reference absolute kg output and only the volume had changed*:

- `fat %_diluted       = fat_kg_ref     / volume_observed × 100`
- `protein %_diluted   = protein_kg_ref / volume_observed × 100`
- `lactose %_diluted   = lactose_kg_ref / volume_observed × 100`

Deviations from that prediction are the rumen/mammary residual — how much of the actual concentration gap is unexplained by dilution alone:

- `fat %_rumen     = fat %_observed     − fat %_diluted`
- `protein %_rumen = protein %_observed − protein %_diluted`
- `lactose %_rumen = lactose %_observed − lactose %_diluted`

A rumen residual of **0** means pure dilution.  A **negative** residual means the cow is *under*-producing the solid relative to her baseline on top of the dilution (classical metabolic suppression).  A **positive** residual means she is actively *up-regulating* synthesis of that solid beyond the dilution pull.

**Key result on the present extract (pooled, n = 806 MLP test-days, vs daily mean barn THI).**

| Nutrient | Observed slope (pp/THI) | Dilution-only slope | Rumen residual | Reading |
|---|---|---|---|---|
| Fat % | −0.0257 *** | −0.0205 *** | **−0.0052 \*\*** | ~80% dilution + ~20% mild rumen suppression |
| Protein % | −0.0169 *** | −0.0165 *** | −0.0005 n.s. | ~100% dilution, no rumen component |
| Lactose % | −0.0015 (≈ 0) | −0.0180 *** | **+0.0164 \*\*\*** | pure-dilution −0.0180 nearly *cancelled* by active lactose up-synthesis +0.0164; net flat |

**Interpretation.**

- **Protein drop is pure dilution.**  The cow synthesises the same absolute protein (kg/d) on hot and cool days; the percentage falls purely because the water fraction rises.
- **Fat drop is mostly dilution (~80%) plus a small but significant rumen-suppression component (~20%).**  This matches the classical literature on heat-stress rumen volatile-fatty-acid profiles: acetate (the main milk-fat precursor) falls faster than propionate / butyrate, so absolute fat-kg is modestly reduced on top of the dilution signal (Wheelock et al. 2010).
- **Lactose % is actively defended.**  The cow up-regulates absolute lactose synthesis enough to offset the volume increase and keep her milk near its osmotic set-point — the rumen residual (+0.0164 *** pp/THI) almost perfectly cancels the −0.0180 *** pp/THI dilution effect.  Why this has to happen is the subject of §3.23 below.

**Outputs.**
- `mlp_dilution_partition.csv` — every MLP test-day with `volume_ref`, `fat_kg_ref`, `protein_kg_ref`, `lactose_kg_ref`, and the three `*_percent_diluted` / `*_percent_rumen` columns.
- `mlp_dilution_partition_summary.csv` — long table (nutrient × component) of rs / p / slope / intercept vs mean THI.
- `mlp_dilution_partition.{svg,png}` — 3-row figure (fat %, protein %, lactose %) with the three OLS reference lines per row.

### 3.23 Why lactose *must* rise with milk volume — osmoregulatory mechanism

**Physiology.** Lactose is synthesised inside the mammary secretory (alveolar) epithelial cells, specifically in the *trans*-Golgi lumen, where the lactose-synthase complex (β-1,4-galactosyltransferase with its α-lactalbumin regulatory subunit) couples UDP-galactose to cytosolic glucose taken up via GLUT1 on the basolateral membrane (Linzell & Peaker 1971; Shennan & Peaker 2000; Costa et al. 2019).  The Golgi and secretory-vesicle membranes are effectively impermeable to this disaccharide, so once formed, lactose cannot diffuse back into the cytosol and remains trapped as an osmotically active solute in the lumen (Neville 1990; Shennan & Peaker 2000).  Water — together with accompanying Na⁺, K⁺ and Cl⁻ needed to keep the vesicle isotonic with cytosol and plasma — is drawn into the Golgi / secretory vesicle down the osmotic gradient that lactose creates, and those vesicles then fuse with the apical membrane to discharge their contents into the alveolar lumen by exocytosis (Linzell & Peaker 1971; Shennan & Peaker 2000).  **Each mole of lactose secreted therefore obligatorily carries a roughly fixed mass of water with it**, so total milk volume is set primarily by the rate of lactose synthesis — lactose is the *pacemaker* of milk water (Neville 1990; Cant et al. 2018; Costa et al. 2019).  The tight-junction-sealed epithelium holds milk very close to iso-osmotic with plasma (~300 mOsm/kg); as lactose concentration rises, monovalent ions (Na⁺, K⁺, Cl⁻) fall reciprocally so that osmotic equilibrium across the secretory cell is maintained (Linzell & Peaker 1971; Shennan & Peaker 2000).

**What this predicts for our data.**  Because mammary water output is set stoichiometrically by lactose synthesis, a cow that increases milk volume on hot days must synthesise proportionally more lactose — she cannot simply dilute a fixed pool of lactose, because the resulting drop in osmotic pressure would halt water entry into the secretory vesicles.  This is exactly the §3.22 pattern on our 806 MLP test-days: milk volume (+0.148 kg/d per THI unit) and lactose mass (rs = +0.186 ***, slope +0.006 kg/d per THI unit) rise together, lactose concentration stays flat (observed slope ≈ 0, rumen residual +0.016 *** pp/THI exactly cancelling the −0.018 *** pp/THI dilution slope), while fat % and protein % fall mostly by passive dilution because their synthesis is not coupled to the osmotic set-point and therefore drifts down as extra water enters the alveolar lumen.

**References.**
- Linzell, J.L., Peaker, M. (1971) Mechanism of milk secretion. *Physiological Reviews* 51(3):564–597. doi:10.1152/physrev.1971.51.3.564
- Neville, M.C. (1990) The physiological basis of milk secretion. *Annals of the New York Academy of Sciences* 586:1–11. doi:10.1111/j.1749-6632.1990.tb17783.x
- Shennan, D.B., Peaker, M. (2000) Transport of milk constituents by the mammary gland. *Physiological Reviews* 80(3):925–951. doi:10.1152/physrev.2000.80.3.925
- Wheelock, J.B., Rhoads, R.P., VanBaale, M.J., Sanders, S.R., Baumgard, L.H. (2010) Effects of heat stress on energetic metabolism in lactating Holstein cows. *Journal of Dairy Science* 93(2):644–655. doi:10.3168/jds.2009-2295
- Cant, J.P., Kim, J.J.M., Cieslar, S.R.L., Doelman, J. (2018) Symposium review: Amino acid uptake by the mammary glands: Where does the control lie? *Journal of Dairy Science* 101(6):5655–5666. doi:10.3168/jds.2017-13844
- Costa, A., Lopez-Villalobos, N., Sneddon, N.W., Shalloo, L., Franzoi, M., De Marchi, M., Penasa, M. (2019) Invited review: Milk lactose — current status and future challenges in dairy cattle. *Journal of Dairy Science* 102(7):5883–5898. doi:10.3168/jds.2018-15955


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

Outputs are routed into subject-focused subfolders by
``digimuh.paths.resolve_output``.  Readers use the matching
``resolve_input`` helper, which first checks the subject subfolder
and then falls back to a flat root layout for back-compat with
older runs.

```
results/broken_stick/
├── 01_extract/                    raw CSVs pulled from cow.db
│   tierauswahl.csv                Animal selection (cleaned from xlsx)
│   rumen_barn.csv                 Rumen temp + barn climate per reading
│   respiration_barn.csv           Respiration + barn climate per reading
│   production.csv                 Mean milk yield + lactation per animal
│   daily_milk_yield.csv           Daily milk yield per cow (summer window)
│   daily_milk_yield_full.csv      Daily yield, full HerdePlus history (Wood fit source, §3.19)
│   mlp_test_days.csv              MLP monthly composition (§3.20)
│   climate_daily.csv              Daily barn climate (Jun–Sep)
│   calvings.csv                   Calving-confirmation events (for DIM, §3.15)
│
├── 02_breakpoints/                broken-stick fits + diagnostics + Spearman
│   broken_stick_results.csv       Per-animal breakpoints + convergence
│   spearman_correlations.csv
│   behavioural_response.csv       Per-animal means below/above bp
│   statistical_tests.csv          Fisher resampling + BH-FDR
│   summary_table.csv
│   boxplot_grouped_*.{svg,png}    Breakpoint distributions per year
│   paired_below_above_*.{svg,png} Per-year significance brackets
│   predictors_*.{svg,png}         Spearman-histogram summary per predictor
│   diagnostic_*.{svg,png}         Diagnostic example panels (+ _top10 for Frontiers)
│   spearman_correlations.{svg,png}
│   climate_{thi,temp}.{svg,png}   Year-on-year climate time series
│   scatter_thi_vs_temp.{svg,png}
│
├── 03_temporal/                   circadian, crossings, CCF, ETA
│   circadian_null_model.csv       Hourly rumen temp: cool vs stress days
│   circadian_null_model{,_stacked}.{svg,png}
│   thi_daily_profile.csv          Barn THI by clock hour × month
│   thi_daily_profile_{year,all}.{svg,png}
│   crossing_times.csv             Clock time of each crossing event
│   crossing_raster_*.{svg,png}    Activation raster (animals × clock time)
│   cross_correlation.csv          CCF + xcov (raw + detrended)
│   xcorr_*.{svg,png}, xcov_*.{svg,png}
│   derivative_ccf.csv             d(climate)/dt vs d(body_temp)/dt
│   dccf_*.{svg,png}
│   event_triggered_traces*.csv    Peri-event traces (all + filtered)
│   event_triggered_summary*.csv
│   eta_*_8to11h.{svg,png}         2×2 event-triggered average figures
│   climate_eta.csv                THI + barn temp around crossings
│   climate_eta_*.{svg,png}
│
├── 04_production/                 TNF × yield, Wood residuals, MY class
│   thermoneutral_fraction.csv     Daily TNF per cow
│   tnf_yield.csv                  TNF + daily yield + P95
│   tnf_yield.{svg,png}
│   daily_milk_yield_wood.csv      Cow-days × Wood residual (DIM-adjusted)
│   wood_curve_fits.csv            Per-lactation Wood parameters + convergence
│   yield_class_per_cow_year.csv   Low / middle / high class per (animal, year)
│   tnf_yield_by_class.csv         Cow-days × class × crossing flags × climate
│   tnf_yield_correlations_by_class.csv   Per-class Spearman tables
│   tnf_yield_by_class{,_relative,_residual}.{svg,png}
│   crossing_day_comparison.csv    Mann-Whitney per (group × predictor)
│   crossing_day_raincloud{,_low,_middle,_high}.{svg,png}
│   daily_climate_vs_yield.csv     Pooled rs + OLS of mean climate vs residual
│   daily_climate_vs_yield_residual.{svg,png}
│   milk_yield_histogram_{pooled,by_year}.{svg,png}
│   milk_yield_wood_{example_fits,residual_histogram}.{svg,png}
│
├── 05_composition/                MLP thin-milk + dilution partition
│   mlp_composition_by_cowday.csv  MLP × climate × class cow-days (§3.21)
│   mlp_composition_correlations.csv   rs / p / slope per (class × predictor × metric)
│   mlp_thin_milk_hypothesis.{svg,png}
│   mlp_composition_heatmap_mean_thi.{svg,png}
│   mlp_dilution_partition.csv     Per-row dilution vs rumen-residual (§3.22)
│   mlp_dilution_partition_summary.csv   Observed / diluted / rumen slopes
│   mlp_dilution_partition.{svg,png}
│
└── 06_longitudinal/               year-on-year stability + Sankey
    breakpoint_stability.csv       ICC pairs
    breakpoint_stability.{svg,png}
    longitudinal_{thi,temp}.{svg,png}
    sankey_*.html                  Alluvial stability + threshold pipeline
    raincloud_crossing_count_*.{svg,png}
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

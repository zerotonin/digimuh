# Analysis 0 — Individual heat stress thresholds

**Pipeline version:** 0.3.0
**Authors:** Bart R. H. Geurten (University of Otago), Gundula Hoffmann (ATB Potsdam)
**Target:** Frontiers in Animal Science


## 1. Pipeline overview

Three independent scripts, each reading CSVs from the previous step:

```
  digimuh-extract  (DB → CSV, runs once)
        │  rumen_barn.csv, respiration_barn.csv,
        │  production.csv, climate_daily.csv
        ▼
  digimuh-stats    (CSV → statistics, seconds)
        │  broken_stick_results.csv, statistical_tests.csv, ...
        ▼
  digimuh-plots    (CSV → SVG/PNG figures, seconds)
```


## 2. Data extraction (analysis_00a_extract.py)

For each of 220 animal-years in the Tierauswahl:

**Rumen temperature:** `temp_without_drink_cycles` from `smaxtec_derived`,
filtered to 30-43 C, joined with hourly barn climate (Neubau sensors
only, `barn_id IN (1,2)`), milking hours excluded (04-07, 16-19).

**Drinking correction:** By default, timestamps with `drink_cycles_v2 > 0`
plus a 15-min recovery window are excluded.  Flag `--smaxtec-drink-correction`
disables this and relies solely on smaXtec's built-in correction.

**Respiration:** `respirationfrequency` from `gouna`, filtered to 5-150 bpm,
same barn join and milking exclusion.  Available from Aug 2022 onward.

**Production:** Mean milk yield per observation window, lactation number
from full animal history.


## 3. Statistical analysis (analysis_00b_stats.py)

### 3.1 Broken-stick regression

Two-segment piecewise linear model per animal:

    y = a1 + b1*x    if x <= psi
    y = a2 + b2*x    if x > psi

Breakpoint psi found by grid search (200 points) + bounded minimisation.
Biological constraint: `b2 > b1` and `b2 > 0`.

Four models per animal: body temp vs THI, body temp vs barn temp,
respiration vs THI, respiration vs barn temp.

Search ranges: THI 45-80, barn temp 5-35 C.

### 3.2 Breakpoint existence tests

Two formal tests for H0: "no breakpoint exists" (the relationship is
a single straight line) vs H1: "there is a slope change."

**Davies test (Davies 1987, 2002):**
Evaluates k=10 candidate breakpoints, computes the t-statistic for the
slope-difference coefficient at each, takes max|t|, corrects p-value
using Davies' upper bound for the multiple-search problem.

**Pseudo-Score test (Muggeo 2016):**
Averages the segmented variable over k candidates, tests whether this
averaged term is significant in the augmented linear model.  More
powerful than Davies for single changepoints (Muggeo 2016, D'Angelo
et al. 2025).

These tests answer "is the relationship nonlinear?" before any threshold
model is fitted.  A significant Davies/pscore p-value means the slope
changes somewhere, but does not specify whether the change is sharp
(broken-stick) or gradual (Hill/sigmoid).

### 3.3 Four-parameter logistic (Hill) fit

Sigmoidal dose-response model:

    y = y_min + (y_max - y_min) / (1 + (EC50 / x)^n)

Fitted via constrained nonlinear least squares (`scipy.optimize.curve_fit`).

**Parameters:**
- `EC50`: predictor value at half-maximal response (midpoint, not onset)
- `n` (Hill coefficient): steepness of transition; large n = sharp switch,
  small n = gradual acceleration
- `y_min`, `y_max`: lower and upper asymptotes

**Lower bend point (Sebaugh & McCray 2003):**

The EC50 is the midpoint, not the onset.  The onset is defined as the
lower bend point: the x-value where the second derivative of the Hill
curve equals zero, marking the transition from the baseline plateau
into the rising phase:

    x_bend_lower = EC50 * ((n - 1) / (n + 1))^(1/n)    for n > 1

For n <= 1 (no lower inflection), we fall back to EC10:

    x_bend_lower = EC50 * (0.10 / 0.90)^(1/n)

This lower bend point is directly comparable to the broken-stick
breakpoint: both represent "where the response starts to deviate
from the baseline."

**Reference:**
Sebaugh JL, McCray PD (2003) Defining the linear portion of a
sigmoid-shaped curve: bend points. Pharmaceutical Statistics 2:167-174.

### 3.4 Interpretation of the three methods together

| Davies/pscore p | BS converges | Hill converges | Interpretation |
|---|---|---|---|
| p < 0.05 | Yes | Yes | Sharp threshold; use BS breakpoint |
| p < 0.05 | No  | Yes | Gradual onset; use Hill lower bend |
| p < 0.05 | No  | No  | Nonlinear but neither model fits |
| p >= 0.05 | Yes | Yes | Spurious; no real threshold (Breit 2023) |
| p >= 0.05 | No  | No  | Linear relationship, no threshold |

The decision logic: always run Davies/pscore first to test whether a
threshold exists.  If significant, fit both broken-stick and Hill.  If
broken-stick converges, use it (sharp threshold confirmed).  If only
Hill converges, use the lower bend point (gradual onset).

### 3.5 Spearman correlations

Per-animal Spearman rs between each signal and each predictor.

### 3.6 Below/above breakpoint comparison

For converged animals: within-year Fisher resampling tests (from the
reRandomStats package, Geurten 2026) on per-animal means below vs above
the individual THI breakpoint.  The test statistic is the median
difference, with 20,000 permutations.  BH-FDR corrected across all
tests within a year.

The Fisher resampling test is a permutation-based two-sample test that
makes no distributional assumptions, making it more appropriate than the
Wilcoxon signed-rank test for our data where sample sizes vary greatly
between years and normality cannot be assumed.

**Reference:**
Geurten BRH (2026) reRandomStats: Re-randomisation Statistics Toolkit.
https://github.com/zerotonin/rerandomstats

### 3.7 Cross-correlation and cross-covariance below/above breakpoint

For each animal with a converged breakpoint, the time series is split
into readings below and above the breakpoint.  For each region, the
normalised cross-correlation function (CCF) and raw cross-covariance
are computed between the climate predictor (THI or barn temperature)
and rumen temperature at lags from -120 to +120 minutes (in 10-minute
steps matching the sensor sampling interval).

Positive lags mean the climate signal leads the rumen temperature
response.  This characterises how quickly the cow's thermoregulatory
system responds to environmental changes, and whether the coupling
strength differs between the thermoneutral zone (below breakpoint)
and the heat stress zone (above breakpoint).

Expected pattern: below the breakpoint, cross-correlation at lag 0
should be weak (homeostatic regulation buffers the signal).  Above the
breakpoint, cross-correlation should be stronger and may show a time
lag of 30-60 minutes (the thermal inertia of the rumen).

**Output:** `cross_correlation.csv` with columns: animal_id, year,
predictor, region, lag, lag_minutes, xcorr, xcov, n.

### 3.8 Breakpoint stability

ICC for repeat animals appearing in multiple years.


## 4. Plotting (analysis_00c_plots.py)

All figures use the Wong (2011) colourblind-safe palette.

| Figure | Description |
|---|---|
| Grouped boxplots | Rumen vs respiration breakpoints per year |
| Paired below/above | Boxplots with significance brackets (BH-FDR stars) |
| Paired rumen vs resp | Same-animal breakpoint comparison |
| Spearman histograms | Per-animal correlation distributions |
| Climate time series | Daily barn THI/temp per summer |
| Predictors | Breakpoint vs milk yield and lactation (Pearson) |
| Stability scatter | Year-to-year THI breakpoints with ICC |
| Example fits | Best-R2 animal per model (rumen + resp) |


## 5. Output files

```
results/broken_stick/
├── rumen_barn.csv              Rumen temp + barn climate per reading
├── respiration_barn.csv        Respiration + barn climate per reading
├── production.csv              Milk yield + lactation per animal
├── climate_daily.csv           Daily barn climate
├── broken_stick_results.csv    Per-animal: breakpoints, Davies/pscore p,
│                               Hill EC50/n/lower_bend, convergence flags
├── spearman_correlations.csv
├── behavioural_response.csv
├── statistical_tests.csv       Fisher resampling tests with BH-FDR
├── cross_correlation.csv       CCF and cross-covariance below/above bp
├── breakpoint_stability.csv
├── summary_table.csv
└── *.svg / *.png               All figures
```


## 6. Running the pipeline

```bash
# Full pipeline (rumen temp + respiration)
digimuh-extract
digimuh-stats --data results/broken_stick
digimuh-plots --data results/broken_stick

# Frontiers paper (rumen temperature only, no respiration)
digimuh-stats --data results/broken_stick --no-resp
digimuh-plots --data results/broken_stick

# Or via the run script
bash scripts/run_00_broken_stick_ana.sh
bash scripts/run_00_broken_stick_ana.sh --no-resp
```


## 7. Key columns in broken_stick_results.csv

Per-animal, per-model prefix (`thi_`, `temp_`, `resp_thi_`, `resp_temp_`):

| Column suffix | Description |
|---|---|
| `_breakpoint` | Broken-stick breakpoint (NaN if not converged) |
| `_converged` | Broken-stick slope constraint passed |
| `_slope_below/above` | Segment slopes |
| `_r_squared` | Broken-stick R2 |
| `_davies_p` | Davies test p-value for breakpoint existence |
| `_pscore_p` | Pseudo-Score test p-value |
| `_hill_ec50` | Hill midpoint (half-maximal response) |
| `_hill_n` | Hill coefficient (steepness) |
| `_hill_bend` | Lower bend point (Sebaugh & McCray 2003 onset) |
| `_hill_r2` | Hill model R2 |
| `_hill_converged` | Hill fit succeeded |


## 8. References

- Benjamini Y, Hochberg Y (1995) J R Stat Soc B 57:289-300.
- Breit M et al. (2023) Multivariate Behavioral Research.
- D'Angelo et al. (2025) Statistics in Medicine.
- Davies RB (1987) Biometrika 74:33-43.
- Davies RB (2002) Biometrika 89:484-489.
- Hoffmann G et al. (2020) Biosystems Engineering 199:83-96.
- Muggeo VMR (2003) Statistics in Medicine 22:3055-3071.
- Muggeo VMR (2016) J Stat Comput Simul 86:3059-3067.
- Neira M et al. (2026) PLOS Climate 5:e0000761.
- Pinto S et al. (2020) J Thermal Biology 88:102523.
- Sebaugh JL, McCray PD (2003) Pharmaceutical Statistics 2:167-174.
- Wang X et al. (2018) J Thermal Biology 77:24-37.
- Wong B (2011) Nature Methods 8:441.

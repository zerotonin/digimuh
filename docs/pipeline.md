# DigiMuh analysis pipeline

**Version:** 0.1.0
**Authors:** Bart R. H. Geurten, Claude (Anthropic)
**Date:** March 2026


## Overview

This document describes the data processing pipeline for the DigiMuh
dairy-cow sensor dataset.  The pipeline transforms ~8.9 GB of heterogeneous
CSV exports from six on-farm monitoring systems into a queryable SQLite
database and runs five analysis modules targeting animal welfare and
production metrics.

```
  CSV exports (8.9 GB, ~4 300 files, 6 sensor systems)
       │
       ▼
  ┌─────────────┐
  │  Ingestion   │──→  SQLite database (star schema)
  └──────┬──────┘     762 M rows in smaxtec_derived alone
         ▼
  ┌─────────────┐
  │  Validation  │──→  Row counts, null rates, value ranges,
  └──────┬──────┘     temporal coverage, referential integrity
         ▼
  ┌─────────────┐
  │  SQL Views   │──→  3-layer view hierarchy:
  └──────┬──────┘     hourly → daily → analysis-specific joins
         ▼
  ┌─────────────┐
  │  Analysis    │──→  5 modules (CSV features + SVG/PNG figures)
  └─────────────┘
```


## Data sources

| System | Location | Measurements | Resolution | Animals/Sensors |
|---|---|---|---|---|
| smaXtec Classic bolus | Reticulum | Temperature, activity, rumination, motility, estrus/calving indices | ~10 min | 837 animals |
| smaXtec pH bolus | Reticulum | Rumen pH, SARA duration | ~10 min | ~6–10% of herd |
| smaXtec water intake | Reticulum | Daily water intake (from temperature drops) | Daily | 837 animals |
| smaXtec barn sensors | Barn walls | Temperature, humidity, THI (NRC 1971) | ~10 min | 4 barns |
| smaXtec events | — | Calving, insemination, pregnancy results | Event-based | 837 animals |
| HerdePlus milking | Milking parlour | Yield, flow, duration, MLP composition | Per milking / monthly MLP | 965 animals |
| HerdePlus diseases | — | Health events, diagnoses, categories | Event-based | Herd-wide |
| Gouna | On-animal | Respiration frequency | ~1 min | 91 animals |
| BCS | Visual assessment | Body Condition Score (1–5) | ~biweekly | 715 animals |
| LoRaWAN | Environmental | Battery level, current | ~2 min | 22 sensors |
| HOBO weather station | Farm | Temp, RH, dew point, solar, wind, wetness | 5 min | 1 station |
| DWD (German Weather Service) | Regional | Daily max THI, max enthalpy | Daily | 1 station |

**Coverage:** April 2021 – September 2024 (3.5 years), ~5 500 animals.


## Step 1 — Ingestion

**Script:** `digimuh-ingest` (or `python -m digimuh.ingest_cow_db`)

**Input:** Directory tree of CSV files, one folder per data source.

**Output:** Single SQLite database file.

**Schema design:** Star schema with four dimension tables and twelve
fact tables.

**Dimension tables:**

| Table | Key | Content |
|---|---|---|
| `animals` | EU ear tag (INTEGER PRIMARY KEY) | 5 572 unique animals |
| `sensors` | Auto-increment ID | 22 LoRaWAN sensor names |
| `barns` | Auto-increment ID | 4 barn locations |
| `source_files` | Auto-increment ID | 4 312 CSV file provenance records |

**Fact tables:** One per data source.  Every row carries a `file_id`
foreign key for full provenance tracing.  Composite indexes on
`(entity_id, timestamp)` enable fast time-range queries per animal.

**Key numbers from full ingestion:**

| Table | Rows |
|---|---|
| `smaxtec_derived` | 762 042 093 |
| `gouna` | 14 222 031 |
| `lorawan` | 4 395 637 |
| `smaxtec_barns` | 1 204 242 |
| `herdeplus` | 983 384 |
| `smaxtec_water_intake` | 503 087 |
| `bcs` | 350 023 |
| `hobo_weather` | 328 282 |
| `allocations` | 92 152 |
| `diseases` | 100 821 |
| `smaxtec_events` | 62 970 |
| `dwd_weather` | 1 323 |

**Runtime:** ~2.7 hours for data insertion (dominated by smaxtec_derived),
plus ~30–60 minutes for index construction.  Performance is I/O-bound;
an internal NVMe SSD is strongly recommended over USB storage.


## Step 2 — Validation

**Script:** `digimuh-validate --db cow.db`

Five automated checks run immediately after ingestion:

1. **Table row counts** — verifies all 16 tables exist and are non-empty.
2. **Null rates** — reports missing-value percentages for key measurement
   columns.  Expected: `ph` is ~90% null (only pH-bolus animals have it);
   `temp` should be < 5% null.
3. **Value ranges** — plausibility checks against physiological bounds
   (rumen temp 30–45 °C, BCS 1–5, respiration 5–120 bpm, etc.).
   Known artefact: gouna reports 255 bpm values (0xFF sensor saturation).
4. **Temporal coverage** — date range per table, confirming April 2021 to
   September 2024.
5. **Referential integrity** — checks for orphaned foreign key references
   between fact and dimension tables.


## Step 3 — SQL views

**Script:** `create_views.sql` (auto-executed on first analysis run)

A three-layer view hierarchy pre-computes the aggregations and joins
needed by the analysis modules.

**Layer 0 — Hourly:** `v_smaxtec_hourly` aggregates the 762 M-row
smaxtec_derived table into hourly means per animal (for circadian
analysis).

**Layer 1 — Daily summaries:**

| View | Aggregates | Key columns |
|---|---|---|
| `v_smaxtec_daily` | smaxtec_derived → daily | temp (mean/min/max/range), activity, rumination, motility, pH, drinking |
| `v_herdeplus_daily` | herdeplus → daily | total milk yield, MLP test-day values (fat, protein, FPR, SCC, urea, ECM) |
| `v_gouna_daily` | gouna → daily | respiration rate (mean/min/max) |
| `v_water_daily` | water_intake → daily | total litres |
| `v_barn_daily` | smaxtec_barns → daily | barn temp, humidity, THI |

**Layer 2 — Analysis joins:** Five views, each joining the daily
summaries needed for a specific research question with disease ground
truth where applicable.


## Step 4 — Analysis modules

Six analysis scripts, each available as a CLI command.  Each produces
CSV feature tables and publication-ready SVG/PNG figures.


### Analysis 0 — Individual heat stress thresholds (broken-stick regression) ¹

**Command:** `digimuh-broken-stick --db cow.db --tierauswahl Tierauswahl.xlsx --out results/broken_stick`

**Rationale:** The Temperature-Humidity Index (THI) threshold at which
rumen temperature begins to rise varies between individuals.  Identifying
per-animal breakpoints enables precision management and provides
ground-truth data for validating population-level THI thresholds (e.g.
THI 68.8 for mild stress onset; Neira et al. 2026).

**Method:**

- For each animal in the collaborator-provided selection list
  (``Tierauswahl.xlsx``, 220 animal-year entries across 2021–2024), joins
  rumen temperature (``temp_without_drink_cycles``) with concurrent barn
  climate sensor readings (barn temperature, barn THI per NRC 1971).
- Excludes milking hours (04:00–07:59 and 16:00–19:59) when cows are
  not in the barn.
- Fits a two-segment piecewise linear (broken-stick) regression per
  animal via grid search + bounded refinement.
- Reports two breakpoints per animal-year:
  (a) rumen temperature vs. barn THI, and
  (b) rumen temperature vs. barn air temperature.
- Below the breakpoint, rumen temperature is stable (thermoneutral zone);
  above it, rumen temperature increases linearly with environmental load.

**Outputs:** Per-animal breakpoint table (CSV), boxplots of breakpoints
across years, THI vs. barn temperature breakpoint scatter, example
animal fit plots.

¹ Analysis led by Dr. med. vet. Gundula Hoffmann, Head of working group
  "Digital monitoring of animal welfare", Leibniz Institute for
  Agricultural Engineering and Bioeconomy (ATB), Potsdam.
  https://www.atb-potsdam.de/en/


### Analysis 1 — Subclinical ketosis detection

**Command:** `digimuh-ketosis --db cow.db --out results/ketosis`

**Rationale:** Subclinical ketosis (negative energy balance) is the most
common metabolic disorder in early-lactation dairy cows.  The milk
fat-to-protein ratio (FPR) is an established indirect indicator: FPR > 1.4
suggests energy deficit, FPR < 1.1 suggests subacute ruminal acidosis.

**Method:**

- Extracts MLP test days with FPR, milk yield, SCC, rumination, rumen
  temperature, pH, and water intake.
- Computes rolling 7-day milk yield deviation per animal.
- Builds a composite ketosis risk score from Z-scored FPR, rumination,
  milk yield deviation, and water intake.
- Trains a Random Forest classifier (200 trees, 5-fold stratified CV)
  to predict FPR > 1.4, validated against disease records.
- Reports accuracy, F1, AUC, and feature importances.

**Outputs:** Feature importance ranking, FPR distribution (healthy vs.
sick), composite risk score distribution.


### Analysis 3 — Heat stress multi-sensor fusion

**Command:** `digimuh-heat --db cow.db --out results/heat`

**Rationale:** Heat stress reduces milk production, fertility, and welfare.
Fixed THI thresholds (e.g. THI > 68) do not account for individual
variation in heat tolerance.

**Method:**

- Per-animal Z-scored rumen temperature (following the NZ smaXtec study,
  JDS Communications 2024): each cow's temperature distribution is scaled
  to a common mean and SD before thresholding (Z > 1.5 = heat stressed).
- Per-animal sigmoid dose-response curve: rumen_temp_z = f(THI).  The
  inflection point represents the animal's personal heat tolerance threshold.
- Composite heat load index fusing rumen temp Z, respiration rate, activity
  suppression, water intake, and rumination.
- Production impact: milk yield loss per heat load quartile.

**Outputs:** THI vs. rumen temperature scatter, per-animal heat tolerance
threshold distribution, milk production by heat load quartile.


### Analysis 6 — Digestive efficiency composite (novel)

**Command:** `digimuh-digestive --db cow.db --out results/digestive`

**Rationale:** Reticulorumen contractions drive mixing, mixing drives
fermentation rate, fermentation determines volatile fatty acid profiles,
and VFA ratios directly shape milk fat and protein.  This causal chain
has a multi-day time lag.

**Method:**

- Time-lagged cross-correlations (1–14 days) between daily motility/pH
  metrics and the next available MLP test-day composition (fat%, protein%,
  FPR, ECM).
- Rolling 7-day motility–pH correlation as a "digestive efficiency score":
  strong negative coupling (faster mixing → lower pH) indicates a
  well-functioning rumen.
- Per-animal herd-percentile ranking of efficiency.

**Outputs:** Lag-correlation heatmaps (predictor × lag → r), digestive
efficiency score distribution.


### Analysis 11 — Circadian rhythm disruption (novel)

**Command:** `digimuh-circadian --db cow.db --out results/circadian`

**Rationale:** Healthy ruminants show strong ~24 h rhythms in core body
temperature (nadir early morning, peak late afternoon), activity (bimodal
dawn/dusk), and rumination (complementary to activity).  Circadian
amplitude collapse or phase shift is a well-established early marker of
sickness in human chronobiology but has barely been explored in precision
dairy farming.

**Method:**

- Single-harmonic Fourier fit (24 h period) per animal-day for temperature,
  activity, and rumination.  Extracts amplitude, acrophase (hour of peak),
  and mesor (24 h mean).
- Circadian Disruption Index (CDI): mean absolute Z-score deviation from the
  animal's own healthy-period baseline across all circadian parameters.
- CDI validated against disease onset from HerdePlus records.

**Outputs:** Circadian amplitude distributions (healthy vs. sick), CDI time
courses for example animals with disease-period shading, CDI distribution
comparison.


### Analysis 12 — Motility pattern entropy (novel)

**Command:** `digimuh-entropy --db cow.db --out results/entropy`

**Rationale:** In a healthy rumen, reticulorumen contractions are
quasi-periodic with modest beat-to-beat variability.  Very low entropy
(rigid, uncoupled contractions) and very high entropy (chaotic,
disorganised contractions) both indicate dysfunction.  This is directly
analogous to heart rate variability (HRV) analysis in cardiology, applied
to the rumen motor complex.  To our knowledge, this approach has not been
published.

**Method:**

- Sample entropy (SampEn, m=2, r=0.2×SD) per Richman & Moorman (2000):
  quantifies self-similarity of the contraction interval series.
- Permutation entropy (PermEn, order=3) per Bandt & Pompe (2002):
  captures ordinal pattern complexity, robust to noise.
- HRV-analogue statistics: mean inter-beat interval, SDNN, coefficient
  of variation, RMSSD.
- Pre-disease entropy trend analysis: compares the 7-day pre-onset window
  against each animal's healthy-period baseline.

**Outputs:** SampEn vs. PermEn scatter (healthy vs. sick), entropy
distributions, entropy vs. rumen pH, pre-disease entropy shift histograms.


## Software and dependencies

| Component | Version |
|---|---|
| Python | >= 3.10 |
| SQLite | (standard library) |
| pandas | data manipulation |
| numpy | numerical computation |
| scipy | sigmoid fitting, statistical tests |
| scikit-learn | Random Forest, cross-validation |
| matplotlib | figure generation |
| tqdm | progress bars (optional) |

Install: `pip install -e ".[analysis]"` or `conda env create -f environment.yml`


## Outputs per analysis

Each analysis script writes to its output directory:

| File type | Content |
|---|---|
| `*_data.csv` / `*_features.csv` | Full feature matrix for further analysis |
| `*.svg` / `*.png` | Publication-ready figures |
| `*_results.json` | Performance metrics and key statistics (where applicable) |


## References

- Hoffmann et al. (2026, in preparation) — Frontiers manuscript: individual heat stress assessment via broken-stick regression of rumen temperature
- Hoffmann et al. (2020) — animal-based heat stress indicators in dairy cattle
- Oetzel (2013) — FPR thresholds for subclinical ketosis monitoring
- Kaufman et al. (2016) J Dairy Sci 99:5604–18 — rumination time and subclinical ketosis
- JDS Communications (2024) — per-animal Z-scored rumen temperature for heat stress
- NRC (1971) — Temperature-Humidity Index formula
- Richman & Moorman (2000) Am J Physiol 278:H2039–49 — sample entropy
- Bandt & Pompe (2002) Phys Rev Lett 88:174102 — permutation entropy

# DigiMuh

> **🚧 Repository under construction** — core ingestion pipeline is functional;
> analysis modules, views, and query helpers are coming next.

**DigiMuh** consolidates ~8.9 GB of heterogeneous dairy-cow CSV sensor data into
a single normalised SQLite database.  The data spans 3.5 years (April 2021 –
September 2024) of continuous monitoring from multiple on-farm systems:

| System | What it measures |
|---|---|
| **smaXtec** bolus | Rumen temperature, pH, activity, motility, rumination, water intake, estrus/calving indices |
| **smaXtec** barn sensors | Barn temperature, humidity, THI |
| **HerdePlus** | Milking events, MLP test-day results, calving/lactation records |
| **HerdePlus** diseases | Health events and diagnoses |
| **Gouna** | Respiration frequency |
| **BCS** | Body condition scores |
| **LoRaWAN** | Environmental sensor battery/current |
| **HOBO** | Weather station (temperature, humidity, solar radiation, wind, wetness) |
| **DWD** | German Weather Service THI and enthalpy |

The database uses a **star schema**: four dimension tables (`animals`, `sensors`,
`barns`, `source_files`) and twelve fact tables, connected by integer foreign
keys.  Every row carries a `file_id` for full provenance tracing back to the
original CSV.

See [`docs/database_structure.md`](docs/database_structure.md) for the full
schema and [`docs/column_dictionary.md`](docs/column_dictionary.md) for a
description of every column.


## Installation

```bash
# Clone the repository
git clone https://github.com/zerotonin/digimuh.git
cd digimuh

# Option A: conda (recommended)
conda env create -f environment.yml
conda activate digimuh

# Option B: pip
pip install -e ".[dev]"
```


## Quick start

```bash
# 1. Smoke test with 5 files per folder (~1 min)
digimuh-ingest /path/to/DigiMuh-Export --db cow_test.db --test-n 5

# 2. Full ingestion (~2–3 hours)
rm cow_test.db
digimuh-ingest /path/to/DigiMuh-Export --db cow.db

# 3. Query the database
python -c "
import sqlite3
con = sqlite3.connect('cow.db')
cur = con.execute('SELECT COUNT(*) FROM smaxtec_derived')
print(f'smaxtec_derived rows: {cur.fetchone()[0]:,}')
"
```


## Expected input layout

The ingestion script expects the DigiMuh CSV export directory to have this
structure:

```
DigiMuh-Export_2021-04-01_2024-09-30/
├── output_allocations/
│   └── allocations.csv
├── outputs_bcs/
│   └── {animal_id}_bcs_{date_range}.csv  ×715
├── outputs_gouna/
│   └── {animal_id}_gouna_{date_range}.csv  ×91
├── outputs_herdeplus_mlp_gemelk_kalbung/
│   └── {animal_id}_herdeplus_{date_range}.csv  ×965
├── outputs_hobo/
│   └── hobo_exports_{date_range}.csv
├── outputs_lorawan/
│   └── {sensor_name}_LoRaWAN_raw_{date_range}.csv  ×22
├── outputs_smaxtec_barns/
│   └── {barn_name}_smaxtec_raw_{date_range}.csv  ×4
├── outputs_smaxtec_derived/
│   └── {animal_id}_smaxtec_derived_{date_range}.csv  ×837
├── outputs_smaxtec_events/
│   └── {animal_id}_events.csv  ×837
├── outputs_smaxtec_water_intake/
│   └── {animal_id}_smaxtec_derived_{date_range}.csv  ×837
├── herdeplus_diseases.csv
└── outputs_dwd.csv
```

Animal IDs are 15-digit EU ear tag numbers.  The entity identifier is always
the first underscore-delimited segment of each filename.


## CLI reference

```
digimuh-ingest [-h] [--db DB] [--chunk-size N] [--verbose] [--test-n N] root_dir
```

| Argument | Description |
|---|---|
| `root_dir` | Root directory containing all CSV folders |
| `--db` | Output SQLite path (default: `cow.db`) |
| `--chunk-size` | Rows per INSERT batch (default: 50 000) |
| `--test-n N` | Only ingest first N files per folder |
| `--verbose`, `-v` | Print CREATE TABLE SQL and debug info |


## Running tests

```bash
python -m pytest
```


## Analysis pipeline

After ingestion, five analysis scripts are available as CLI commands.  Each
creates analysis views on first run, queries the database, and writes results
(CSV data + figures) to an output directory.

```bash
# Install with analysis dependencies
pip install -e ".[analysis]"

# 1. Subclinical ketosis risk — FPR × rumination × milk yield
digimuh-ketosis --db cow.db --out results/ketosis

# 2. Heat stress — rumen temp × THI × water × respiration
digimuh-heat --db cow.db --out results/heat

# 3. Digestive efficiency — motility × pH → milk composition (time-lagged)
digimuh-digestive --db cow.db --out results/digestive

# 4. Circadian disruption — 24h Fourier decomposition as welfare marker
digimuh-circadian --db cow.db --out results/circadian

# 5. Motility entropy — rumen HRV analogue via information theory
digimuh-entropy --db cow.db --out results/entropy
```

Each script writes:
- A CSV of the extracted features (for further analysis in R, Python, etc.)
- Publication-ready SVG + PNG figures
- A JSON summary of key results (where applicable)

See [`docs/database_structure.md`](docs/database_structure.md) for the SQL view
definitions that power these analyses.


## Roadmap

- [x] CSV → SQLite ingestion with star schema
- [x] SQL views for analysis (daily summaries + cross-table joins)
- [x] Analysis: subclinical ketosis detection (FPR + RF classifier)
- [x] Analysis: heat stress multi-sensor fusion
- [x] Analysis: digestive efficiency (motility–pH coupling)
- [x] Analysis: circadian rhythm disruption index
- [x] Analysis: motility pattern entropy (novel)
- [ ] Data validation and quality-check reports
- [ ] Parallelised entropy computation for full dataset
- [ ] Sphinx documentation on GitHub Pages


## Authors

**Bart R. H. Geurten** — Department of Zoology, University of Otago

## License

MIT

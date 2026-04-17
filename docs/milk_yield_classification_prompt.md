## Task

Analyse the daily milk yield distribution in the DigiMuh dataset and produce
a publication-quality histogram that shows how well existing literature
classification boundaries fit our herd.

## Context

Read `milk_yield_classification.md` in this directory for full background.
The repo is at `~/PyProjects/digimuh`.  Follow `CODING_STYLE_GUIDE.md` in
the repo root.

## Data sources

All CSVs are in `results/broken_stick/`:

1. `daily_milk_yield.csv` — one row per cow per day, columns include
   `animal_id`, `year`, `date`, `daily_yield_kg` (sum of all milkings that day).
   Check the actual column names with `head -1`.

2. `production.csv` — one row per animal-year, includes `mean_milk_yield_kg`.
   Use this as a secondary check.

## Steps

1. **Load and inspect** `daily_milk_yield.csv`.  Print shape, column names,
   and basic stats (`describe()`).  Drop rows where `daily_yield_kg` (or
   whatever the yield column is called) is NaN or ≤0.

2. **Histogram — all years pooled**.  Horizontal axis: daily milk yield (kg/d).
   Vertical axis: count (or density).  Bin width ~1 kg.  Use Wong palette
   from `constants.py`.  Overlay:
   - **Müschner-Siemens boundaries**: vertical dashed lines at 28.8 and 38.4
     (label "Müschner-Siemens 2020")
   - **Yan boundaries**: vertical dotted lines at 26 and 39
     (label "Yan 2021")
   - **Our tercile boundaries**: solid vertical lines at the 33rd and 67th
     percentiles of the pooled data (label with actual values)
   - Annotate n, median, mean on the plot.

3. **Per-year facet**.  Same histogram but faceted by year (2021–2024),
   stacked vertically sharing x-axis.  Same boundary lines.  Per-panel n
   labels.

4. **Summary table**.  Print to console:
   - Per year: n_days, n_animals, mean, median, Q33, Q67, min, max
   - Pooled: same stats
   - Müschner-Siemens class counts: how many cow-days fall in L/M/H
   - Yan class counts: same
   - Our tercile class counts: same

5. **Save figures** as SVG+PNG to `results/broken_stick/` using
   `save_figure()` from `digimuh.viz_base`.  File names:
   `milk_yield_histogram_pooled` and `milk_yield_histogram_by_year`.

## Style

- Use `setup_figure()` from `digimuh.viz_base` for matplotlib config
- Use `save_figure()` for output
- Wong (2011) colourblind-safe palette from `digimuh.constants`
- No bar charts — histogram is fine here since it's a distribution
- British English in labels
- SVG-first with `svg.fonttype = "none"`

## Output

After running, tell me:
- Where do our terciles fall relative to the literature boundaries?
- Which classification scheme best fits our data?
- Are there enough cows in each group for meaningful broken-stick analysis?

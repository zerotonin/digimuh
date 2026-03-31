# Step 1: extract (slow, hits DB — run once) | --no-drink-exclusion means that we rely on the smaxtech drink exclusion
digimuh-extract --db ~/cow.db --tierauswahl Tierauswahl.xlsx --out results/broken_stick --no-drink-exclusion

# Step 2: stats (fast, reads CSVs — re-run when changing model)
digimuh-stats --data results/broken_stick

# Step 3: plots (fast, reads CSVs — re-run when tweaking figures)
digimuh-plots --data results/broken_stick
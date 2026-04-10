#!/usr/bin/env bash
# ╔══════════════════════════════════════════════════════════════╗
# ║  Run Analysis 0 — Broken-stick pipeline                      ║
# ║  Uses config from ~/.config/digimuh/config.yaml              ║
# ║  CLI args override config. Run digimuh-config to set up.     ║
# ║                                                              ║
# ║  Usage:                                                      ║
# ║    bash scripts/run_00_broken_stick_ana.sh                   ║
# ║    bash scripts/run_00_broken_stick_ana.sh --tierauswahl X.xlsx  ║
# ╚══════════════════════════════════════════════════════════════╝

set -e

# Step 1: extract (slow, hits DB)
# All paths come from config; pass any extra CLI args through
digimuh-extract "$@"

# Step 2: stats (fast, reads CSVs)
digimuh-stats --data results/broken_stick

# Step 3: plots (fast, reads CSVs)
digimuh-plots --data results/broken_stick
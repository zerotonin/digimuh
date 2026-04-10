#!/usr/bin/env bash
# ╔══════════════════════════════════════════════════════════════╗
# ║  Run Analysis 0 — Broken-stick pipeline                      ║
# ║  Uses config from ~/.config/digimuh/config.yaml              ║
# ║  CLI args override config. Run digimuh-config to set up.     ║
# ║                                                              ║
# ║  Usage:                                                      ║
# ║    bash scripts/run_00_broken_stick_ana.sh                   ║
# ║    bash scripts/run_00_broken_stick_ana.sh --no-resp          ║
# ║    bash scripts/run_00_broken_stick_ana.sh --tierauswahl X.xlsx  ║
# ╚══════════════════════════════════════════════════════════════╝

set -e

# Step 1: extract (slow, hits DB)
digimuh-extract "$@"

# Step 2: stats (fast, reads CSVs)
# Pass --no-resp if present in args
NO_RESP=""
for arg in "$@"; do
    if [ "$arg" = "--no-resp" ]; then
        NO_RESP="--no-resp"
    fi
done
digimuh-stats --data results/broken_stick $NO_RESP

# Step 3: plots (fast, reads CSVs)
digimuh-plots --data results/broken_stick
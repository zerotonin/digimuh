#!/usr/bin/env bash
# ╔══════════════════════════════════════════════════════════════╗
# ║  DEPRECATED — use run_frontiers_2026.sh or digimuh-frontiers ║
# ╚══════════════════════════════════════════════════════════════╝
echo "⚠  DEPRECATED: use 'digimuh-frontiers' or 'scripts/run_frontiers_2026.sh'"
echo "   Forwarding arguments …"
echo ""

set -e

# Separate stats-only flags from extract args
NO_RESP=""
FRONTIERS=""
EXTRACT_ARGS=()
for arg in "$@"; do
    if [ "$arg" = "--no-resp" ]; then
        NO_RESP="--no-resp"
    elif [ "$arg" = "--frontiers" ]; then
        FRONTIERS="--frontiers"
    else
        EXTRACT_ARGS+=("$arg")
    fi
done

digimuh-extract "${EXTRACT_ARGS[@]}"
digimuh-stats --data results/broken_stick $NO_RESP $FRONTIERS
digimuh-plots --data results/broken_stick

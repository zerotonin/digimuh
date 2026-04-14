#!/usr/bin/env bash
# ╔══════════════════════════════════════════════════════════════╗
# ║  Run Frontiers 2026 analysis pipeline                        ║
# ║  « thin wrapper around digimuh-frontiers »                   ║
# ╚══════════════════════════════════════════════════════════════╝
set -e
digimuh-frontiers "$@"

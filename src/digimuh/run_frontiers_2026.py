#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — run_frontiers_2026                                   ║
# ║  « end-to-end pipeline for the Frontiers paper »               ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Orchestrates extract → stats → plots for the 2026 Frontiers   ║
# ║  in Animal Science submission (Hoffmann & Geurten).             ║
# ║                                                                 ║
# ║  Usage:                                                         ║
# ║    digimuh-frontiers                                            ║
# ║    digimuh-frontiers --no-resp --frontiers                      ║
# ╚══════════════════════════════════════════════════════════════════╝
"""End-to-end runner for the Frontiers in Animal Science paper.

Calls extract, stats, and plots in sequence.  Equivalent to the
old ``scripts/run_00_broken_stick_ana.sh`` but as a Python entry
point with proper argument forwarding.
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys

log = logging.getLogger("digimuh.runner")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the full Frontiers 2026 analysis pipeline")
    parser.add_argument("--no-resp", action="store_true",
                        help="Skip respiration analysis (only 65 animals)")
    parser.add_argument("--frontiers", action="store_true",
                        help="Skip Davies/pscore/Hill (reserved for COMPAG)")
    parser.add_argument("--data", type=str, default="results/broken_stick",
                        help="Output directory for CSVs and figures")
    parser.add_argument("--skip-extract", action="store_true",
                        help="Skip DB extraction (reuse existing CSVs)")
    # Forward remaining args to extract
    parser.add_argument("extract_args", nargs="*",
                        help="Additional arguments for digimuh-extract")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s │ %(name)-20s │ %(message)s",
        datefmt="%H:%M:%S",
    )

    # Step 1: Extract (slow, hits DB)
    if not args.skip_extract:
        log.info("═══ Step 1/3: Extract ═══")
        cmd = [sys.executable, "-m", "digimuh.extract"] + args.extract_args
        subprocess.run(cmd, check=True)

    # Step 2: Stats (fast, reads CSVs)
    log.info("═══ Step 2/3: Stats ═══")
    cmd = [sys.executable, "-m", "digimuh.stats_runner",
           "--data", args.data]
    if args.no_resp:
        cmd.append("--no-resp")
    if args.frontiers:
        cmd.append("--frontiers")
    subprocess.run(cmd, check=True)

    # Step 3: Plots (fast, reads CSVs)
    log.info("═══ Step 3/3: Plots ═══")
    cmd = [sys.executable, "-m", "digimuh.viz_runner",
           "--data", args.data]
    subprocess.run(cmd, check=True)

    log.info("═══ Pipeline complete ═══")


if __name__ == "__main__":
    main()

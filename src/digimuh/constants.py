#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — constants                                            ║
# ║  « one source of truth for colours, thresholds, figure rules »  ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Wong (2011) colourblind-safe palette with semantic mappings    ║
# ║  to breakpoint analysis categories.  Import this module         ║
# ║  instead of hardcoding hex values or magic numbers.             ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Central constants for all DigiMuh analyses and visualisations."""

from __future__ import annotations

# ─────────────────────────────────────────────────────────────────
#  Wong (2011) base palette
# ─────────────────────────────────────────────────────────────────

WONG_BLUE       = "#0072B2"
WONG_ORANGE     = "#E69F00"
WONG_GREEN      = "#009E73"
WONG_VERMILLION = "#D55E00"
WONG_SKY        = "#56B4E9"
WONG_PINK       = "#CC79A7"
WONG_YELLOW     = "#F0E442"
WONG_GREY       = "#999999"


# ─────────────────────────────────────────────────────────────────
#  Semantic colour mappings
# ─────────────────────────────────────────────────────────────────

COLOURS = {
    "year": {
        2021: WONG_BLUE,
        2022: WONG_ORANGE,
        2023: WONG_GREEN,
        2024: WONG_PINK,
    },
    "below_bp":    WONG_BLUE,
    "above_bp":    WONG_VERMILLION,
    "fit_line":    WONG_VERMILLION,
    "reference":   WONG_ORANGE,
    "identity":    WONG_GREY,
    "scatter":     WONG_SKY,
    "scatter_alt": WONG_GREEN,
    "hist_thi":    WONG_BLUE,
    "hist_temp":   WONG_GREEN,
    "box_thi":     WONG_SKY,
    "box_temp":    WONG_GREEN,
    "median":      WONG_VERMILLION,
    "paired_line": WONG_GREY,
}

SANKEY_COLOURS = {
    "strongly decreased": WONG_BLUE,
    "decreased":          WONG_SKY,
    "stable":             WONG_GREY,
    "increased":          WONG_ORANGE,
    "strongly increased": WONG_VERMILLION,
}


# ─────────────────────────────────────────────────────────────────
#  Breakpoint classification thresholds
# ─────────────────────────────────────────────────────────────────

THI_REFERENCE = 68.8
"""THI mild-stress threshold from Hoffmann et al. (2020)."""

DELTA_STABLE   = 1.0
"""Max |Δ breakpoint| classified as 'stable'."""

DELTA_MODERATE = 3.0
"""Max |Δ breakpoint| classified as 'decreased' / 'increased'."""

MIN_COHORT_SIZE = 10
"""Minimum animals for longitudinal Sankey plots."""

RESAMPLING_SEED = 20260721
"""Seed for every Fisher re-randomisation test in the pipeline.

Resampling draws a finite sample of the permutation space, so an unseeded
run returns a slightly different p-value each time — measured drift is
about 0.001-0.007, enough to move a borderline call across alpha.  Fixing
one seed makes every reported p-value reproducible; report it in the
methods section alongside the resampling count.
"""


# ─────────────────────────────────────────────────────────────────
#  Figure defaults
# ─────────────────────────────────────────────────────────────────

RCPARAMS = {
    # Editable text in vector exports — text survives the SVG/PDF
    # round-trip into Adobe Illustrator as real characters rather
    # than vector paths.  This is the strongest constraint the
    # Frontiers production team requires, and we honour it ahead of
    # font-family fidelity per author instruction.
    "svg.fonttype":   "none",   # SVG <text> elements (editable)
    "pdf.fonttype":    42,       # TrueType in PDF (editable in Illustrator)
    "ps.fonttype":     42,       # TrueType in EPS
    # Font family: Arial first per Frontiers production guidelines,
    # falling back to Liberation Sans (metrically identical to Arial)
    # on rendering machines without Arial installed.  The recorded
    # family name in the SVG is the first one matplotlib resolves;
    # the SVG opens with Arial in Illustrator on any Adobe-licensed
    # workstation regardless of which font matplotlib used at render.
    "font.family":     ["Arial", "Liberation Sans", "DejaVu Sans",
                        "sans-serif"],
    # Font sizes per Frontiers guideline (≥ 7-8 pt at final print
    # size; single-column max width 8.5 cm, double-column 17.8 cm,
    # absolute max 20 cm).
    "font.size":         8,    # default text on figures
    "axes.titlesize":    9,    # subplot titles
    "axes.labelsize":    8,    # axis labels
    "xtick.labelsize":   7,
    "ytick.labelsize":   7,
    "legend.fontsize":   7,
    "figure.titlesize":  10,   # suptitle
    # Resolution: Frontiers raster minimum is 300 DPI.  Even though
    # we export SVG/PNG (the user finalises layout in Illustrator),
    # the PNG previews must already meet the production threshold.
    "figure.dpi":      150,    # screen preview
    "savefig.dpi":     300,    # production-grade rasters
    "savefig.bbox":    "tight",
    "savefig.facecolor": "white",
    "savefig.edgecolor": "white",
    "savefig.transparent": False,
}

MPL_STYLE = "seaborn-v0_8-whitegrid"


# ─────────────────────────────────────────────────────────────────
#  Broken-stick fit defaults
# ─────────────────────────────────────────────────────────────────

THI_RANGE  = (45, 80)
"""Default x_range for THI → body temp fits."""

TEMP_RANGE = (5, 35)
"""Default x_range for barn temp → body temp fits."""

MIN_READINGS = 50
"""Minimum readings per animal-year to attempt a fit."""

GRID_STEPS = 200
"""Number of breakpoint candidates in grid search."""

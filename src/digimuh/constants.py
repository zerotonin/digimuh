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


# ─────────────────────────────────────────────────────────────────
#  Figure defaults
# ─────────────────────────────────────────────────────────────────

RCPARAMS = {
    "svg.fonttype":   "none",
    "font.family":    "sans-serif",
    "font.size":      11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "figure.dpi":     150,
    "savefig.dpi":    300,
    "savefig.bbox":   "tight",
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

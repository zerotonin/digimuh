#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — paths                                                ║
# ║  « subject-based subfolder routing for the results tree »      ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Every stats/viz module funnels its outputs through             ║
# ║  ``resolve_output(data_dir, filename)`` so the results          ║
# ║  directory stays organised by analysis subject instead of       ║
# ║  piling 150+ files into one flat folder.  Readers use the       ║
# ║  same helper to find files wherever they live.                  ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Central routing of result artefacts into subject-specific subfolders.

The six subjects correspond to the §3.x sections in
``docs/analysis_00_methods.md``:

    01_extract        — raw per-animal / per-cow-day CSVs pulled from the DB
    02_breakpoints    — broken-stick fits, diagnostics, Spearman
    03_temporal       — circadian, crossings, cross-correlation, ETA
    04_production     — TNF × yield, Wood residuals, milk-yield class
    05_composition    — MLP thin-milk + dilution partition
    06_longitudinal   — year-on-year stability, Sankey, raincloud

Call sites use ``resolve_output(data_dir, filename)`` on both write
and read paths.  When a file's subject cannot be determined the
helper falls back to the flat root so unrecognised outputs still
land somewhere predictable.
"""

from __future__ import annotations

from pathlib import Path


SUBDIRS: dict[str, str] = {
    "extract":      "01_extract",
    "breakpoints":  "02_breakpoints",
    "temporal":     "03_temporal",
    "production":   "04_production",
    "composition":  "05_composition",
    "longitudinal": "06_longitudinal",
}


# Exact stem → subject.  `Path(filename).stem` strips .csv / .png /
# .svg / .html, so all three extension variants of the same figure
# land in the same subfolder.
_STEM_TO_SUBJECT: dict[str, str] = {
    # ── 01 extract ────────────────────────────────────────
    "tierauswahl":              "extract",
    "rumen_barn":               "extract",
    "respiration_barn":         "extract",
    "production":               "extract",
    "daily_milk_yield":         "extract",
    "daily_milk_yield_full":    "extract",
    "climate_daily":            "extract",
    "calvings":                 "extract",
    "mlp_test_days":            "extract",

    # ── 02 breakpoints ────────────────────────────────────
    "broken_stick_results":     "breakpoints",
    "spearman_correlations":    "breakpoints",
    "behavioural_response":     "breakpoints",
    "statistical_tests":        "breakpoints",
    "summary_table":            "breakpoints",
    "scatter_thi_vs_temp":      "breakpoints",
    "scatter_bodytemp_vs_resp_bp": "breakpoints",

    # ── 03 temporal ───────────────────────────────────────
    "circadian_null_model":         "temporal",
    "circadian_null_model_stacked": "temporal",
    "thi_daily_profile":            "temporal",
    "crossing_times":               "temporal",
    "crossing_raster_thi":          "temporal",
    "crossing_raster_temp":         "temporal",
    "cross_correlation":            "temporal",
    "derivative_ccf":               "temporal",
    "event_triggered_traces":       "temporal",
    "event_triggered_traces_filtered":  "temporal",
    "event_triggered_summary":      "temporal",
    "event_triggered_summary_filtered": "temporal",
    "climate_eta":                  "temporal",
    "climate_eta_thi":              "temporal",
    "climate_eta_temp":             "temporal",
    "climate_thi":                  "temporal",
    "climate_temp":                 "temporal",

    # ── 04 production ─────────────────────────────────────
    "tnf_yield":                        "production",
    "thermoneutral_fraction":           "production",
    "daily_milk_yield_wood":            "production",
    "wood_curve_fits":                  "production",
    "yield_class_per_cow_year":         "production",
    "tnf_yield_by_class":               "production",
    "tnf_yield_by_class_relative":      "production",
    "tnf_yield_by_class_residual":      "production",
    "tnf_yield_correlations_by_class":  "production",
    "crossing_day_comparison":          "production",
    "crossing_day_raincloud":           "production",
    "crossing_day_raincloud_low":       "production",
    "crossing_day_raincloud_middle":    "production",
    "crossing_day_raincloud_high":      "production",
    "daily_climate_vs_yield":           "production",
    "daily_climate_vs_yield_residual":  "production",
    "daily_climate_vs_yield_relative":  "production",
    "milk_yield_histogram_pooled":      "production",
    "milk_yield_histogram_by_year":     "production",
    "milk_yield_wood_example_fits":     "production",
    "milk_yield_wood_residual_histogram": "production",

    # ── 05 composition ────────────────────────────────────
    "mlp_composition_by_cowday":        "composition",
    "mlp_composition_correlations":     "composition",
    "mlp_thin_milk_hypothesis":         "composition",
    "mlp_dilution_partition":           "composition",
    "mlp_dilution_partition_summary":   "composition",

    # ── 06 longitudinal ───────────────────────────────────
    "breakpoint_stability":             "longitudinal",
}


# Ordered prefix rules for parameterised names (years, per-predictor,
# per-class, detrended/overlay variants, …).  First match wins, so put
# the most specific prefix first.
_PREFIX_RULES: tuple[tuple[str, str], ...] = (
    # breakpoints
    ("diagnostic_",              "breakpoints"),
    ("boxplot_grouped_",         "breakpoints"),
    ("paired_below_above_",      "breakpoints"),
    ("predictors_",              "breakpoints"),

    # temporal
    ("thi_daily_profile_",       "temporal"),
    ("xcorr_",                   "temporal"),
    ("xcov_",                    "temporal"),
    ("dccf_",                    "temporal"),
    ("eta_thi_",                 "temporal"),
    ("eta_temp_",                "temporal"),
    ("climate_eta_",             "temporal"),

    # production
    ("tnf_yield_",               "production"),
    ("crossing_day_raincloud",   "production"),
    ("daily_climate_vs_yield",   "production"),
    ("milk_yield_",              "production"),

    # composition
    ("mlp_composition_heatmap",  "composition"),
    ("mlp_",                     "composition"),

    # longitudinal
    ("longitudinal_",            "longitudinal"),
    ("sankey_",                  "longitudinal"),
    ("raincloud_crossing_count_","longitudinal"),
)


def _subject_for_stem(stem: str) -> str | None:
    """Return the subject key for a file stem, or ``None`` if unknown."""
    if stem in _STEM_TO_SUBJECT:
        return _STEM_TO_SUBJECT[stem]
    for prefix, subject in _PREFIX_RULES:
        if stem.startswith(prefix):
            return subject
    return None


def resolve_output(data_dir: Path, filename: str, *,
                   create: bool = True) -> Path:
    """Return the canonical path for an output file inside ``data_dir``.

    When the filename's subject is in the routing table, the file
    lands in (and is read from) the corresponding subject subfolder.
    When it isn't, the file stays at the flat root.  The target
    directory is created on demand unless ``create=False``.

    Args:
        data_dir: The analysis output directory (e.g.
            ``results/broken_stick``).
        filename: A basename (no path components) like
            ``"rumen_barn.csv"`` or ``"diagnostic_thi_top10.svg"``.
        create: When True (default), ``mkdir(parents=True,
            exist_ok=True)`` is called on the resolved folder before
            it is returned, so callers never have to.

    Returns:
        ``data_dir / SUBDIRS[subject] / filename`` when a subject is
        known, otherwise ``data_dir / filename``.
    """
    stem = Path(filename).stem
    subject = _subject_for_stem(stem)
    folder = data_dir if subject is None else data_dir / SUBDIRS[subject]
    if create:
        folder.mkdir(parents=True, exist_ok=True)
    return folder / filename


def resolve_input(data_dir: Path, filename: str) -> Path:
    """Return the path to read ``filename`` from.

    Preference order:
      1. Subject subfolder (if known from the routing table).
      2. Flat root (legacy layout, for back-compat with older runs).

    Does not create directories.  Returns whichever path exists; if
    neither exists, returns the subject-subfolder path so callers
    get a consistent "canonical" error via ``FileNotFoundError``.
    """
    subject_path = resolve_output(data_dir, filename, create=False)
    if subject_path.exists():
        return subject_path
    legacy_path = data_dir / filename
    if legacy_path.exists():
        return legacy_path
    return subject_path

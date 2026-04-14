#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — test_imports                                         ║
# ║  « smoke tests for all public modules »                        ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Smoke tests: every public module must import without errors."""

from __future__ import annotations

import importlib
import pytest


@pytest.mark.parametrize("module_name", [
    "digimuh.config",
    "digimuh.console",
    "digimuh.constants",
    "digimuh.extract",
    "digimuh.fitting",
    "digimuh.ingest",
    "digimuh.run_frontiers_2026",
    "digimuh.stats_core",
    "digimuh.stats_longitudinal",
    "digimuh.stats_production",
    "digimuh.stats_runner",
    "digimuh.stats_temporal",
    "digimuh.validate_db",
    "digimuh.viz_base",
    "digimuh.viz_breakpoints",
    "digimuh.viz_correlation",
    "digimuh.viz_longitudinal",
    "digimuh.viz_production",
    "digimuh.viz_runner",
    "digimuh.viz_temporal",
    "digimuh.analysis_utils",
    "digimuh.analysis_01_ketosis",
    "digimuh.analysis_03_heat_stress",
    "digimuh.analysis_06_digestive",
    "digimuh.analysis_11_circadian",
    "digimuh.analysis_12_motility_entropy",
])
def test_import(module_name: str) -> None:
    """Module imports without raising."""
    importlib.import_module(module_name)


def test_version() -> None:
    import digimuh
    parts = digimuh.__version__.split(".")
    assert len(parts) >= 2

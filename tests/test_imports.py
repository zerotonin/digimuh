# ╔══════════════════════════════════════════════════════════════╗
# ║  DigiMuh — test_imports                                      ║
# ║  « verify every module loads without a database »             ║
# ╚══════════════════════════════════════════════════════════════╝
"""Smoke tests: every public module must import without errors."""
from __future__ import annotations

import importlib

import pytest


MODULES = [
    "digimuh",
    "digimuh.config",
    "digimuh.console",
    "digimuh.analysis_utils",
    "digimuh.analysis_00a_extract",
    "digimuh.analysis_00b_stats",
    "digimuh.analysis_00c_plots",
    "digimuh.analysis_00_broken_stick",
    "digimuh.analysis_01_ketosis",
    "digimuh.analysis_03_heat_stress",
    "digimuh.analysis_06_digestive",
    "digimuh.analysis_11_circadian",
    "digimuh.analysis_12_motility_entropy",
    "digimuh.ingest_cow_db",
    "digimuh.validate_db",
]


@pytest.mark.parametrize("module_name", MODULES)
def test_import(module_name: str) -> None:
    """Module imports without raising."""
    importlib.import_module(module_name)


def test_version_string() -> None:
    """Package exposes a semantic version."""
    import digimuh
    assert hasattr(digimuh, "__version__")
    parts = digimuh.__version__.split(".")
    assert len(parts) >= 2, "Version must be at least major.minor"

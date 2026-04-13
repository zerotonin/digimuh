# ╔══════════════════════════════════════════════════════════════╗
# ║  DigiMuh — Sphinx documentation configuration               ║
# ║  « Napoleon + autodoc + Furo »                               ║
# ╚══════════════════════════════════════════════════════════════╝
from __future__ import annotations

import sys
from pathlib import Path

# -- Path setup --------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# -- Project information -----------------------------------------------
project = "DigiMuh"
author = "Bart R. H. Geurten"
copyright = "2024–2026, Bart R. H. Geurten"  # noqa: A001

from digimuh import __version__  # noqa: E402
release = __version__
version = __version__

# -- General configuration ---------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "myst_parser",
]

# Napoleon settings (Google-style docstrings)
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_use_param = True
napoleon_use_rtype = True

# Autodoc settings
autodoc_member_order = "bysource"
autodoc_typehints = "description"

# Mock imports that are heavy or unavailable on CI runners
autodoc_mock_imports = [
    "kaleido",
    "plotly",
    "rich",
    "rerandomstats",
    "openpyxl",
    "yaml",
    "pyyaml",
]

# MyST for .md files
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Intersphinx
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# -- Options for HTML output -------------------------------------------
html_theme = "furo"
html_title = f"DigiMuh {version}"
html_static_path = []
html_theme_options = {
    "source_repository": "https://github.com/zerotonin/digimuh",
    "source_branch": "main",
    "source_directory": "docs/",
}

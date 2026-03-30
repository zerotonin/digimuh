#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  DigiMuh analysis utilities                                  ║
# ║  « shared plumbing for all analysis scripts »                ║
# ╚══════════════════════════════════════════════════════════════╝
"""Shared database connection, view initialisation, and plotting
defaults used across all DigiMuh analysis modules.

Usage::

    from digimuh.analysis_utils import connect_db, setup_plotting
    con = connect_db(Path("cow.db"))
    setup_plotting()
"""

from __future__ import annotations

import logging
import sqlite3
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

log = logging.getLogger("digimuh.analysis")

# ─────────────────────────────────────────────────────────────
#  « database helpers »
# ─────────────────────────────────────────────────────────────

VIEWS_SQL = Path(__file__).parent / "create_views.sql"


def connect_db(db_path: Path, create_views: bool = True) -> sqlite3.Connection:
    """Open the cow database and optionally create analysis views.

    Args:
        db_path: Path to the SQLite database file.
        create_views: If ``True``, execute ``create_views.sql`` to
            ensure all analysis views exist.

    Returns:
        An open ``sqlite3.Connection`` with ``Row`` factory enabled.
    """
    if not db_path.is_file():
        raise FileNotFoundError(f"Database not found: {db_path}")

    con = sqlite3.connect(str(db_path))
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA journal_mode = WAL")
    con.execute("PRAGMA cache_size = -256000")  # 256 MB

    if create_views and VIEWS_SQL.is_file():
        log.info("Creating/refreshing analysis views …")
        con.executescript(VIEWS_SQL.read_text(encoding="utf-8"))
        log.info("Views ready.")

    return con


def query_df(con: sqlite3.Connection, sql: str, params: tuple = ()) -> "pd.DataFrame":
    """Execute a query and return a pandas DataFrame.

    Args:
        con: Active database connection.
        sql: SQL query string.
        params: Bind parameters for the query.

    Returns:
        A ``pandas.DataFrame`` with the query results.
    """
    import pandas as pd
    return pd.read_sql_query(sql, con, params=params)


# ─────────────────────────────────────────────────────────────
#  « plotting defaults »
# ─────────────────────────────────────────────────────────────

def setup_plotting() -> None:
    """Configure matplotlib for publication-quality figures.

    Sets SVG text as editable (``svg.fonttype = "none"``),
    uses a clean style, and configures reasonable defaults for
    axis labels and tick sizes.
    """
    plt.style.use("seaborn-v0_8-whitegrid")
    mpl.rcParams.update({
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })


def save_fig(fig: plt.Figure, name: str, out_dir: Path) -> None:
    """Save a figure as SVG, PNG, and close it.

    Args:
        fig: The matplotlib figure to save.
        name: Base filename (without extension).
        out_dir: Output directory.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "png"):
        path = out_dir / f"{name}.{ext}"
        fig.savefig(path)
        log.info("Saved: %s", path)
    plt.close(fig)

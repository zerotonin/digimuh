#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — viz_base                                             ║
# ║  « shared plotting setup, save, and helper functions »          ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  All viz_*.py modules import setup_figure() and save_figure()   ║
# ║  from here.  Never call plt.style.use or rcParams directly in   ║
# ║  a plotting module — use these helpers instead.                 ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Shared matplotlib configuration and figure I/O for all plots."""

from __future__ import annotations

import logging
from pathlib import Path

from digimuh.constants import RCPARAMS, MPL_STYLE

log = logging.getLogger("digimuh.viz")


def setup_figure() -> None:
    """Configure matplotlib for publication-quality figures."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    plt.style.use(MPL_STYLE)
    mpl.rcParams.update(RCPARAMS)


def save_figure(fig, name: str, out_dir: Path) -> None:
    """Save figure as SVG + PNG and close."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(out_dir / f"{name}.{ext}")
    import matplotlib.pyplot as plt
    plt.close(fig)


def add_significance_bracket(
    ax, x1: float, x2: float, y: float, stars: str,
    h: float = 0.02, lw: float = 1.2,
) -> None:
    """Draw a significance bracket with stars between two x positions.

    Args:
        ax: Matplotlib axes.
        x1, x2: Left and right x positions.
        y: Y position of the bracket bottom (data coords).
        stars: Text to display (e.g. '***', 'n.s.').
        h: Bracket height as fraction of y-range.
        lw: Line width.
    """
    if not stars or stars == "":
        return
    yrange = ax.get_ylim()[1] - ax.get_ylim()[0]
    dh = h * yrange
    ax.plot([x1, x1, x2, x2], [y, y + dh, y + dh, y], color="#333333",
            linewidth=lw, clip_on=False)
    ax.text((x1 + x2) / 2, y + dh, stars, ha="center", va="bottom",
            fontsize=11, fontweight="bold", color="#333333")

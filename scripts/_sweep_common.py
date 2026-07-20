# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — _sweep_common                                         ║
# ║  « one grid, so two sweeps answer one question »                 ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Each memory-timescale sweep used to carry its own tau grid.     ║
# ║  Three shared a 19-point log grid; one swept 72 linear hours,    ║
# ║  so a peak from that script was not comparable point-for-point   ║
# ║  with a peak from the others.  Both grids are named here, and    ║
# ║  each script states which one it uses.                           ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Memory-timescale grids and progress reporting for the sweeps."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import TypeVar

T = TypeVar("T")

# Wide-range log grid: 0.5 h to 3 days, roughly even on a log axis, so the
# R2(tau) and convergence(tau) curves stay readable and every cross-index
# comparison is made on identical timescales.
TAU_GRID_LOG: list[float] = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0,
                             10.0, 12.0, 16.0, 20.0, 24.0, 30.0, 36.0,
                             48.0, 60.0, 72.0]

# Dense integer grid for localising a single peak.  Finer than TAU_GRID_LOG
# above ~6 h and coarser below 1 h, so an optimum found on this grid is NOT
# interchangeable with one found on the log grid — always report which was
# used before quoting a timescale.
TAU_GRID_DENSE: list[float] = [float(t) for t in range(1, 73)]

TAU_GRIDS: dict[str, list[float]] = {
    "log": TAU_GRID_LOG,
    "dense": TAU_GRID_DENSE,
}


def progress(items: Iterable[T], total: int | None = None,
             desc: str = "sweeping") -> Iterator[T]:
    """Wrap an iterable in a tqdm bar, degrading to a plain pass-through.

    Feed this a *lazy* iterator.  ``multiprocessing.Pool.map`` blocks until
    every task is finished, so wrapping it shows 0 % and then 100 %; use
    ``imap`` instead, which yields as results arrive.  (``ProcessPoolExecutor``
    ``.map`` is already lazy and needs no change.)
    """
    try:
        from tqdm import tqdm
    except ImportError:
        return iter(items)
    return tqdm(items, total=total, desc=desc,
                bar_format="  {l_bar}{bar:30}{r_bar}")

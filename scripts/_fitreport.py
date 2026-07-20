# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — _fitreport                                            ║
# ║  « who failed, and why — no silent survivorship »                ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  A fit that raises used to vanish on a bare `continue`, so a     ║
# ║  code error looked just like genuine non-convergence.  These     ║
# ║  helpers tally each swallowed exception by type, so the sweep    ║
# ║  CSVs carry an audit trail and the log says whether any fit      ║
# ║  actually broke.                                                 ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Failure accounting shared by the per-animal-year fitting sweeps."""

from __future__ import annotations

import logging
from collections import Counter

import pandas as pd

FailureTally = Counter


def tally(failures: FailureTally, exc: BaseException) -> None:
    """Record one raised fit, keyed by exception type."""
    failures[type(exc).__name__] += 1


def format_failures(failures: FailureTally) -> str:
    """Render a tally as ``Name:count`` pairs for a single CSV cell."""
    return ";".join(f"{name}:{n}" for name, n in sorted(failures.items()))


def parse_failures(cell: object) -> FailureTally:
    """Inverse of :func:`format_failures`; tolerates NaN and empty cells."""
    out: FailureTally = Counter()
    if not isinstance(cell, str) or not cell:
        return out
    for part in cell.split(";"):
        name, _, n = part.partition(":")
        if name and n.isdigit():
            out[name] += int(n)
    return out


def report(df: pd.DataFrame, names: list[str], log: logging.Logger) -> int:
    """Log whether any fit raised, and return the total.

    Silence here is the point: it tells the reader that the convergence
    rates in the same table reflect genuine non-convergence rather than
    code that broke.  Raised fits are *not* excluded from the denominator
    — they still count as non-converged — so a non-zero total means the
    reported rates are pessimistic by that many animal-years.
    """
    total = int(sum(int(df[f"{n}_fail"].sum()) for n in names))
    if total == 0:
        log.info("fit failures: none — convergence reflects the data, not bugs")
        return 0

    kinds: FailureTally = Counter()
    for name in names:
        for cell in df[f"{name}_fail_kinds"]:
            kinds.update(parse_failures(cell))
    log.warning(
        "fit failures: %d raised across the sweep (%s); they are counted as "
        "non-converged, so the rates below are pessimistic — see the *_fail "
        "columns", total,
        ", ".join(f"{k}:{v}" for k, v in kinds.most_common()),
    )
    return total

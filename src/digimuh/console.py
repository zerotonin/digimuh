#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════╗
# ║  DigiMuh rich console output                                  ║
# ║  « coloured panels, tables, progress bars »                  ║
# ╚══════════════════════════════════════════════════════════════╝
"""Pretty-printed console output for DigiMuh analysis pipelines.

Uses ``rich`` for coloured panels, tables, progress bars, and tree
views.  Falls back to plain logging if rich is not installed.

Usage::

    from digimuh.console import console, section, result_table, progress

    console.print("[bold blue]Starting analysis …[/]")
    with progress("Fitting animals") as pb:
        task = pb.add_task("Broken-stick", total=220)
        for animal in animals:
            fit(animal)
            pb.advance(task)
    section("Results", "Breakpoint analysis complete")
    result_table("Convergence", headers, rows)
"""

from __future__ import annotations

import logging
from contextlib import contextmanager
from typing import Any

log = logging.getLogger("digimuh")

try:
    from rich.console import Console
    from rich.logging import RichHandler
    from rich.panel import Panel
    from rich.progress import (
        Progress, SpinnerColumn, BarColumn, TextColumn,
        TimeElapsedColumn, TimeRemainingColumn, MofNCompleteColumn,
    )
    from rich.table import Table
    from rich.theme import Theme
    from rich.tree import Tree
    HAS_RICH = True
except ImportError:
    HAS_RICH = False


# ─────────────────────────────────────────────────────────────
#  « theme and console »
# ─────────────────────────────────────────────────────────────

if HAS_RICH:
    _theme = Theme({
        "info": "cyan",
        "warn": "yellow",
        "error": "bold red",
        "success": "bold green",
        "section": "bold blue",
        "value": "bold white",
        "dim": "dim",
        "stars.3": "bold green",
        "stars.2": "bold yellow",
        "stars.1": "yellow",
        "stars.ns": "dim",
    })
    console = Console(theme=_theme, highlight=False)
else:
    console = None


def setup_logging(level: int = logging.INFO) -> None:
    """Configure logging with rich handler if available."""
    if HAS_RICH:
        logging.basicConfig(
            level=level,
            format="%(message)s",
            datefmt="%H:%M:%S",
            handlers=[RichHandler(
                console=console,
                show_path=False,
                show_time=True,
                markup=True,
            )],
        )
    else:
        logging.basicConfig(
            level=level,
            format="%(asctime)s │ %(name)-20s │ %(message)s",
            datefmt="%H:%M:%S",
        )


# ─────────────────────────────────────────────────────────────
#  « section headers »
# ─────────────────────────────────────────────────────────────

_step_counter = 0


def section(title: str, subtitle: str = "") -> None:
    """Print a section header."""
    global _step_counter
    _step_counter += 1

    if HAS_RICH:
        label = f"[bold]{_step_counter}. {title}[/]"
        if subtitle:
            label += f"\n[dim]{subtitle}[/]"
        console.print()
        console.print(Panel(label, border_style="blue", width=60))
    else:
        log.info("═" * 50)
        log.info("%d. %s", _step_counter, title)
        if subtitle:
            log.info("   %s", subtitle)


def reset_steps() -> None:
    """Reset the step counter (for testing)."""
    global _step_counter
    _step_counter = 0


# ─────────────────────────────────────────────────────────────
#  « result tables »
# ─────────────────────────────────────────────────────────────

def result_table(
    title: str,
    headers: list[str],
    rows: list[list[Any]],
    highlight_col: int | None = None,
) -> None:
    """Print a formatted results table."""
    if HAS_RICH:
        table = Table(title=title, show_lines=False, pad_edge=True)
        for i, h in enumerate(headers):
            style = "bold" if i == highlight_col else ""
            table.add_column(h, style=style, justify="right" if i > 0 else "left")
        for row in rows:
            str_row = []
            for i, cell in enumerate(row):
                s = _format_cell(cell)
                if i == highlight_col and isinstance(cell, str) and cell.startswith("*"):
                    s = f"[bold green]{s}[/]"
                str_row.append(s)
            table.add_row(*str_row)
        console.print(table)
    else:
        log.info("  %s", title)
        hdr = " | ".join(f"{h:>12s}" for h in headers)
        log.info("  %s", hdr)
        log.info("  %s", "-" * len(hdr))
        for row in rows:
            log.info("  %s", " | ".join(f"{_format_cell(c):>12s}" for c in row))


def _format_cell(val: Any) -> str:
    """Format a single cell value for display."""
    if val is None:
        return ""
    if isinstance(val, float):
        if abs(val) < 0.0001 and val != 0:
            return f"{val:.2e}"
        return f"{val:.3f}"
    if isinstance(val, int):
        if abs(val) < 10000:
            return str(val)
        return f"{val:,}"
    return str(val)


def stars_styled(stars: str) -> str:
    """Return rich-styled significance stars."""
    if not HAS_RICH:
        return stars
    if stars == "***":
        return "[stars.3]***[/]"
    if stars == "**":
        return "[stars.2]**[/]"
    if stars == "*":
        return "[stars.1]*[/]"
    return "[stars.ns]n.s.[/]"


# ─────────────────────────────────────────────────────────────
#  « key-value displays »
# ─────────────────────────────────────────────────────────────

def kv(key: str, value: Any, indent: int = 2) -> None:
    """Print a key-value pair."""
    pad = " " * indent
    if HAS_RICH:
        console.print(f"{pad}[dim]{key}:[/] [value]{_format_cell(value)}[/]")
    else:
        log.info("%s%s: %s", pad, key, _format_cell(value))


def kv_pair(key: str, val1: Any, val2: Any, sep: str = " / ") -> None:
    """Print a key with two values (e.g. converged / total)."""
    if HAS_RICH:
        console.print(
            f"  [dim]{key}:[/] [value]{_format_cell(val1)}[/]"
            f"[dim]{sep}[/][value]{_format_cell(val2)}[/]"
        )
    else:
        log.info("  %s: %s%s%s", key, _format_cell(val1), sep, _format_cell(val2))


# ─────────────────────────────────────────────────────────────
#  « progress bars »
# ─────────────────────────────────────────────────────────────

@contextmanager
def progress(description: str = "Processing"):
    """Context manager for a rich progress bar."""
    if HAS_RICH:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(bar_width=30),
            MofNCompleteColumn(),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=console,
            transient=False,
        ) as pb:
            yield pb
    else:
        # Minimal fallback: a dummy that logs progress
        class _DummyProgress:
            def add_task(self, desc, total=None):
                self._total = total or 0
                self._done = 0
                self._desc = desc
                return 0
            def advance(self, task_id, advance=1):
                self._done += advance
                if self._done % 20 == 0 or self._done == self._total:
                    log.info("  [%d/%d] %s", self._done, self._total, self._desc)
        yield _DummyProgress()


# ─────────────────────────────────────────────────────────────
#  « banner »
# ─────────────────────────────────────────────────────────────

def banner(title: str, version: str = "0.1.0") -> None:
    """Print the DigiMuh startup banner."""
    if HAS_RICH:
        console.print(Panel(
            f"[bold blue]{title}[/]\n[dim]v{version}[/]",
            border_style="blue",
            width=60,
            title="[bold]DigiMuh[/]",
            subtitle="[dim]Precision Livestock Farming[/]",
        ))
    else:
        log.info("=" * 60)
        log.info("  DigiMuh — %s (v%s)", title, version)
        log.info("=" * 60)


def done(message: str = "All analyses complete.") -> None:
    """Print completion message."""
    if HAS_RICH:
        console.print()
        console.print(Panel(
            f"[bold green]{message}[/]",
            border_style="green",
            width=60,
        ))
    else:
        log.info("═" * 50)
        log.info(message)

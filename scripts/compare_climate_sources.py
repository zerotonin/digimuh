#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — compare_climate_sources                              ║
# ║  « what changes when the predictor is HOBO, not smaXtec »       ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Reads broken_stick_results.csv from two analysis output dirs   ║
# ║  (the smaXtec-barn run and the HOBO run) and writes a Markdown  ║
# ║  report comparing coverage, broken-stick convergence, threshold ║
# ║  existence (Davies), and breakpoint location per predictor and  ║
# ║  per year.                                                      ║
# ║                                                                  ║
# ║  Usage:                                                          ║
# ║    python scripts/compare_climate_sources.py \                  ║
# ║        --smaxtec results/broken_stick \                         ║
# ║        --hobo    results/broken_stick_hobo                      ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Compare smaXtec-barn vs HOBO broken-stick results into a report."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

PREDICTORS = (("thi", "THI breakpoint"), ("temp", "Barn-temp breakpoint (°C)"))


def _load_bs(data_dir: Path) -> pd.DataFrame:
    """Load broken_stick_results.csv from a run's 02_breakpoints folder."""
    for candidate in (data_dir / "02_breakpoints" / "broken_stick_results.csv",
                      data_dir / "broken_stick_results.csv"):
        if candidate.exists():
            return pd.read_csv(candidate)
    raise FileNotFoundError(f"broken_stick_results.csv not found under {data_dir}")


def _summarise(bs: pd.DataFrame, prefix: str) -> dict:
    """Coverage / convergence / location stats for one predictor prefix."""
    conv_col = f"{prefix}_converged"
    bp_col = f"{prefix}_breakpoint"
    davies_col = f"{prefix}_davies_p"

    n_rows = len(bs)
    conv = bs[bs[conv_col] == True] if conv_col in bs.columns else bs.iloc[0:0]
    vals = conv[bp_col].dropna() if bp_col in conv.columns else pd.Series(dtype=float)
    davies_sig = None
    if davies_col in bs.columns and bs[davies_col].notna().any():
        davies_sig = int((bs[davies_col] < 0.05).sum())
    return {
        "n_animal_years": n_rows,
        "n_animals": int(bs["animal_id"].nunique()),
        "n_converged": int(len(conv)),
        "pct_converged": 100.0 * len(conv) / n_rows if n_rows else float("nan"),
        "median": float(vals.median()) if len(vals) else float("nan"),
        "q25": float(vals.quantile(0.25)) if len(vals) else float("nan"),
        "q75": float(vals.quantile(0.75)) if len(vals) else float("nan"),
        "min": float(vals.min()) if len(vals) else float("nan"),
        "max": float(vals.max()) if len(vals) else float("nan"),
        "davies_sig": davies_sig,
        "vals": vals,
    }


def _per_year(bs: pd.DataFrame, prefix: str) -> pd.Series:
    """Median converged breakpoint per year."""
    conv_col, bp_col = f"{prefix}_converged", f"{prefix}_breakpoint"
    if conv_col not in bs.columns or bp_col not in bs.columns:
        return pd.Series(dtype=float)
    conv = bs[bs[conv_col] == True].dropna(subset=[bp_col])
    return conv.groupby(conv["year"].astype(int))[bp_col].median()


def _fmt(x: float, nd: int = 2) -> str:
    return "—" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{nd}f}"


def build_report(smaxtec_dir: Path, hobo_dir: Path) -> str:
    bs_s = _load_bs(smaxtec_dir)
    bs_h = _load_bs(hobo_dir)
    lines: list[str] = []
    add = lines.append

    add("# DigiMuh — climate-source comparison (smaXtec vs HOBO)\n")
    add(f"- **smaXtec run:** `{smaxtec_dir}`")
    add(f"- **HOBO run:** `{hobo_dir}`\n")
    add("HOBO has no native THI; it is derived with the NRC (1971) dairy "
        "formula from HOBO air temperature and relative humidity, whereas the "
        "smaXtec `temp_hum_index` is proprietary.  **THI values are therefore "
        "not directly comparable between sources** — read the barn-temperature "
        "rows as the clean comparison, and the THI rows for structure only.\n")

    for prefix, label in PREDICTORS:
        s = _summarise(bs_s, prefix)
        h = _summarise(bs_h, prefix)
        add(f"## {label}\n")
        add("| Metric | smaXtec | HOBO | Δ (HOBO−smaXtec) |")
        add("|---|---:|---:|---:|")
        add(f"| Animal-years (total) | {s['n_animal_years']} | "
            f"{h['n_animal_years']} | {h['n_animal_years']-s['n_animal_years']:+d} |")
        add(f"| Converged (n) | {s['n_converged']} | {h['n_converged']} | "
            f"{h['n_converged']-s['n_converged']:+d} |")
        add(f"| Converged (%) | {_fmt(s['pct_converged'],1)} | "
            f"{_fmt(h['pct_converged'],1)} | "
            f"{_fmt(h['pct_converged']-s['pct_converged'],1)} |")
        if s["davies_sig"] is not None and h["davies_sig"] is not None:
            add(f"| Davies p<0.05 (n) | {s['davies_sig']} | {h['davies_sig']} | "
                f"{h['davies_sig']-s['davies_sig']:+d} |")
        add(f"| Breakpoint median | {_fmt(s['median'])} | {_fmt(h['median'])} | "
            f"{_fmt(h['median']-s['median'])} |")
        add(f"| Breakpoint IQR | {_fmt(s['q25'])}–{_fmt(s['q75'])} | "
            f"{_fmt(h['q25'])}–{_fmt(h['q75'])} | |")
        add(f"| Breakpoint range | {_fmt(s['min'])}–{_fmt(s['max'])} | "
            f"{_fmt(h['min'])}–{_fmt(h['max'])} | |\n")

        py_s, py_h = _per_year(bs_s, prefix), _per_year(bs_h, prefix)
        years = sorted(set(py_s.index) | set(py_h.index))
        if years:
            add("Per-year median breakpoint:\n")
            add("| Year | smaXtec | HOBO | Δ |")
            add("|---|---:|---:|---:|")
            for y in years:
                vs = py_s.get(y, np.nan)
                vh = py_h.get(y, np.nan)
                dd = vh - vs if not (np.isnan(vs) or np.isnan(vh)) else np.nan
                add(f"| {y} | {_fmt(vs)} | {_fmt(vh)} | {_fmt(dd)} |")
            add("")

    add("## Reading the result\n")
    add("- **Coverage / convergence**: a drop in converged animal-years under "
        "HOBO points to thinner climate coverage or a weaker temperature "
        "signal at the loggers.")
    add("- **Barn-temp breakpoint shift**: HOBO is an outdoor / Stallklima "
        "logger set, so its temperature axis differs from the in-barn smaXtec "
        "sensors; a systematic offset in the median breakpoint is expected and "
        "is the headline number.")
    add("- **Davies count**: how many animals still show a detectable threshold "
        "at all — the test of whether the breakpoint phenomenon survives the "
        "change of predictor.")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smaxtec", type=Path, default=Path("results/broken_stick"))
    parser.add_argument("--hobo", type=Path, default=Path("results/broken_stick_hobo"))
    parser.add_argument("--out", type=Path, default=None,
                        help="Report path (default: <hobo>/COMPARISON_REPORT.md)")
    args = parser.parse_args()

    report = build_report(args.smaxtec, args.hobo)
    out = args.out or (args.hobo / "COMPARISON_REPORT.md")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(report, encoding="utf-8")
    print(report)
    print(f"\nReport written to: {out}")


if __name__ == "__main__":
    main()

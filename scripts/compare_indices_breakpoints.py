# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — compare_indices_breakpoints                           ║
# ║  « which thermal index best marks the rumen-temp threshold »     ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Joins the full HOBO index set onto the extracted rumen readings ║
# ║  and fits one broken-stick (+ Davies existence test) per         ║
# ║  animal-year for EVERY index, using the same reRandomStats       ║
# ║  primitives as the main pipeline.  The output is Part I Fig 2:   ║
# ║  one fixed method, many indices, ranked by how cleanly they      ║
# ║  produce a reticulorumen-temperature breakpoint.                 ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Compare every thermal index by the rumen-temperature breakpoint it yields.

Reads the HOBO index table (``compute_hobo_indices.py``) and the extracted
rumen readings, nearest-time-merges them, then fits a broken-stick per
animal-year for each index and summarises convergence, explained variance,
and threshold existence (Davies).  Run:

    PYTHONPATH=<path to heataxis>/src \\
        python scripts/compare_indices_breakpoints.py --jobs 20
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rerandomstats import broken_stick_fit, davies_test

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _sweep_common import progress  # noqa: E402

log = logging.getLogger("index_compare")

MIN_READINGS = 50           # per animal-year, matches the main pipeline
X_LO, X_HI = 5.0, 95.0      # percentile bounds for each index's breakpoint search

# Predictors: the barn air temperature (reference) + every heataxis index.
REFERENCE = "barn_temp"


def merged_readings(indices_csv: Path, rumen_csv: Path) -> tuple[pd.DataFrame, list[str]]:
    """Nearest-time-merge the HOBO indices onto the rumen readings."""
    idx = pd.read_csv(indices_csv)
    idx["datetime"] = pd.to_datetime(idx["datetime"], utc=True)
    reserved = {"datetime", "ta", "rh", "tdp", "sr", "u", "Tw", "Tbg", "year", "day"}
    index_cols = [c for c in idx.columns if c not in reserved]

    rumen = pd.read_csv(rumen_csv)
    rumen["timestamp"] = pd.to_datetime(rumen["timestamp"], utc=True)
    rumen = rumen.dropna(subset=["body_temp"]).sort_values("timestamp")

    idx = idx.sort_values("datetime")
    merged = pd.merge_asof(
        rumen, idx[["datetime", *index_cols]],
        left_on="timestamp", right_on="datetime",
        direction="nearest", tolerance=pd.Timedelta(minutes=10))
    merged = merged.dropna(subset=index_cols, how="all")
    predictors = [REFERENCE, *index_cols]
    log.info("merged %d rumen readings; %d predictors", len(merged), len(predictors))
    return merged, predictors


def x_ranges(df: pd.DataFrame, predictors: list[str]) -> dict[str, tuple[float, float]]:
    """Per-predictor breakpoint search window from robust percentiles."""
    out = {}
    for col in predictors:
        s = df[col].to_numpy()
        s = s[np.isfinite(s)]
        out[col] = (float(np.percentile(s, X_LO)), float(np.percentile(s, X_HI)))
    return out


# ─────────────────────────────────────────────────────────────
#  Per-animal-year worker (module level -> picklable for pools)
# ─────────────────────────────────────────────────────────────

def fit_animal_year(payload: dict) -> list[dict]:
    """Fit a broken-stick + Davies for every predictor on one animal-year."""
    grp = payload["grp"]
    ranges = payload["ranges"]
    y = grp["body_temp"].to_numpy()
    rows = []
    for col, xr in ranges.items():
        x = grp[col].to_numpy()
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < MIN_READINGS:
            rows.append({**payload["key"], "index": col, "n": int(ok.sum()),
                         "converged": False, "breakpoint": np.nan,
                         "r_squared": np.nan, "davies_p": np.nan})
            continue
        fit = broken_stick_fit(x[ok], y[ok], x_range=xr)
        try:
            dp, davies_error = davies_test(x[ok], y[ok], x_range=xr)["pvalue"], ""
        except Exception as exc:
            dp, davies_error = np.nan, type(exc).__name__
        rows.append({
            **payload["key"], "index": col, "n": int(ok.sum()),
            "converged": bool(fit.get("converged", False)),
            "breakpoint": fit.get("breakpoint", np.nan),
            "r_squared": fit.get("r_squared", np.nan),
            "slope_below": fit.get("slope_below", np.nan),
            "slope_above": fit.get("slope_above", np.nan),
            "davies_p": dp,
            "davies_error": davies_error,
        })
    return rows


def run_fits(merged: pd.DataFrame, predictors: list[str], jobs: int,
             limit: int | None) -> pd.DataFrame:
    ranges = x_ranges(merged, predictors)
    groups = list(merged.groupby(["animal_id", "year", "date_enter"]))
    if limit:
        groups = groups[:limit]
    payloads = [
        {"key": {"animal_id": aid, "year": year, "date_enter": enter},
         "grp": grp, "ranges": ranges}
        for (aid, year, enter), grp in groups
    ]
    log.info("fitting %d animal-years x %d predictors (jobs=%d) ...",
             len(payloads), len(predictors), jobs)

    results: list[dict] = []
    if jobs and jobs > 1:
        # ProcessPoolExecutor.map returns a lazy iterator, so the bar advances
        # as results arrive (unlike multiprocessing.Pool.map, which blocks).
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            for rows in progress(pool.map(fit_animal_year, payloads),
                                 total=len(payloads), desc="animal-years"):
                results.extend(rows)
    else:
        for payload in progress(payloads, total=len(payloads),
                                desc="animal-years"):
            results.extend(fit_animal_year(payload))
    return pd.DataFrame(results)


# ─────────────────────────────────────────────────────────────
#  Summary and figure
# ─────────────────────────────────────────────────────────────

def summarise(fits: pd.DataFrame) -> pd.DataFrame:
    """Per-index convergence, explained variance, and threshold existence."""
    rows = []
    for col, g in fits.groupby("index"):
        conv = g[g["converged"]]
        # A Davies test that raised leaves davies_p NaN, which compares False
        # and so deflates pct_davies_sig.  Report the error rate alongside it
        # rather than letting a code failure read as "no threshold".
        n_err = int((g.get("davies_error", "").astype(str) != "").sum())
        rows.append({
            "index": col,
            "n_animal_years": len(g),
            "pct_converged": 100.0 * g["converged"].mean(),
            "median_r2": float(conv["r_squared"].median()) if len(conv) else np.nan,
            "median_breakpoint": float(conv["breakpoint"].median()) if len(conv) else np.nan,
            "pct_davies_sig": 100.0 * (g["davies_p"] < 0.05).mean(),
            "n_davies_error": n_err,
        })
    out = pd.DataFrame(rows).sort_values("median_r2", ascending=False)
    total_err = int(out["n_davies_error"].sum())
    if total_err:
        log.warning("Davies test raised on %d index x animal-year fits; those "
                    "count as non-significant in pct_davies_sig", total_err)
    else:
        log.info("Davies test: no failures")
    return out.reset_index(drop=True)


def plot_summary(summ: pd.DataFrame, out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from heataxis.constants import WONG
    from heataxis.viz import save_figure, setup_figure
    setup_figure()

    s = summ.sort_values("median_r2")
    y = np.arange(len(s))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.5, 0.42 * len(s) + 1.6),
                                   sharey=True)

    def _colour(name: str) -> str:
        if name == REFERENCE:
            return WONG["vermilion"]            # air temperature (reference)
        if name.startswith("leakyTHI"):
            return WONG["bluish_green"]         # history — THI memory
        if name.startswith("leakyITSC"):
            return WONG["reddish_purple"]       # history — ITSC memory
        return WONG["blue"]                     # instantaneous indices

    colours = [_colour(n) for n in s["index"]]
    ax1.hlines(y, 0, s["median_r2"], color=colours, lw=1.2, alpha=0.5)
    ax1.scatter(s["median_r2"], y, color=colours, s=42, zorder=3)
    ax1.set_yticks(y, s["index"], fontsize=8)
    ax1.set_xlabel("median R² (broken-stick, converged)")
    ax1.set_title("explained variance", fontsize=10, loc="left")

    ax2.hlines(y, 0, s["pct_converged"], color=colours, lw=1.2, alpha=0.5)
    ax2.scatter(s["pct_converged"], y, color=colours, s=42, zorder=3)
    ax2.set_xlabel("% animal-years converged")
    ax2.set_title("breakpoint detectability", fontsize=10, loc="left")

    from matplotlib.lines import Line2D
    legend = [
        Line2D([0], [0], marker="o", color="w", label="instantaneous index",
               markerfacecolor=WONG["blue"], markersize=7),
        Line2D([0], [0], marker="o", color="w", label="history — leaky THI",
               markerfacecolor=WONG["bluish_green"], markersize=7),
        Line2D([0], [0], marker="o", color="w", label="history — leaky ITSC",
               markerfacecolor=WONG["reddish_purple"], markersize=7),
        Line2D([0], [0], marker="o", color="w", label="air temperature (ref.)",
               markerfacecolor=WONG["vermilion"], markersize=7),
    ]
    ax1.legend(handles=legend, loc="lower right", fontsize=7.5, frameon=False)
    fig.suptitle("Does the x-axis need a memory? Index vs rumen-temperature "
                 "threshold (broken-stick, HOBO 2021–2024)", fontsize=11)
    save_figure(fig, "index_breakpoint_comparison", out_dir, dpi=300)
    plt.close(fig)


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                        datefmt="%H:%M:%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indices", type=Path,
                        default=Path("results/hobo_indices/hobo_indices_5min.csv.gz"))
    parser.add_argument("--rumen", type=Path,
                        default=Path("results/broken_stick_hobo/01_extract/rumen_barn.csv"))
    parser.add_argument("--out", type=Path, default=Path("results/index_comparison"))
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--limit", type=int, default=None,
                        help="cap animal-years (quick test)")
    parser.add_argument("--months", type=str, default="5,6,7,8,9",
                        help="warm-season months to keep (fair window for all "
                             "predictors); empty string = whole year")
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    merged, predictors = merged_readings(args.indices, args.rumen)
    if args.months.strip():
        months = {int(m) for m in args.months.split(",")}
        before = len(merged)
        merged = merged[merged["timestamp"].dt.month.isin(months)]
        log.info("summer restriction %s: kept %d of %d readings",
                 sorted(months), len(merged), before)
    fits = run_fits(merged, predictors, args.jobs, args.limit)
    fits.to_csv(args.out / "per_animal_year_fits.csv", index=False, encoding="utf-8")
    summ = summarise(fits)
    summ.to_csv(args.out / "index_ranking.csv", index=False, encoding="utf-8")
    plot_summary(summ, args.out)

    log.info("done -> %s", args.out)
    print("\nIndex ranking (by median R² on the rumen-temperature breakpoint):")
    print(summ.to_string(index=False, float_format=lambda v: f"{v:.2f}"))


if __name__ == "__main__":
    main()

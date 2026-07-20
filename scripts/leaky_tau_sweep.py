# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — leaky_tau_sweep                                       ║
# ║  « how long a memory best marks the rumen-temp threshold? »      ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Sweeps the leaky-integrator memory constant tau = 1..72 h for   ║
# ║  THI and ITSC, fitting a broken-stick per animal-year at each    ║
# ║  tau, and plots explained variance and convergence vs tau.       ║
# ║  The peak is the biological memory timescale.                    ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Sweep the leaky-integrator memory tau for THI and ITSC.

Leaky-integrates each base index on the continuous HOBO station series (real
per-step gaps), maps every rumen reading to its station sample once, then fits
a broken-stick per animal-year at each tau.  Output: median R² and convergence
as a function of tau, for THI and ITSC.  Run:

    PYTHONPATH=<path to heataxis>/src \\
        python scripts/leaky_tau_sweep.py --jobs 20

``--tau-grid`` selects the timescales.  The default ``dense`` grid (72 linear
hours) localises a single peak finely; ``log`` uses the 19-point grid the
other sweeps share.  An optimum found on one grid is **not** interchangeable
with one found on the other, so always report which was used.
"""

from __future__ import annotations

import argparse
import logging
import multiprocessing as mp
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from heataxis import history
from heataxis.constants import THI_STRESS_ONSET
from rerandomstats import broken_stick_fit, hill_fit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fitreport import format_failures, report, tally  # noqa: E402
from _sweep_common import TAU_GRIDS, progress  # noqa: E402

LENSES = {"bs": broken_stick_fit, "hill": hill_fit}
SERIES_NAMES = ["THI", "ITSC", "TSD", "TSL"]

log = logging.getLogger("tau_sweep")

MIN_READINGS = 50
SUMMER = {5, 6, 7, 8, 9}
X_LO, X_HI = 5.0, 95.0

# Read-only globals shared with worker processes (via fork, no pickling).
_THI = _ITSC = _DT = None
_ITSC_BASE = 0.0
_ANIMAL_YEARS: list[tuple[np.ndarray, np.ndarray]] = []
_ALL_POS = None
_T_HOURS = _CUM_LOAD = _CUM_OVER = None   # for windowed TSD / TSL


def _init_globals(thi, itsc, dt, itsc_base, animal_years, all_pos,
                  t_hours, cum_load, cum_over) -> None:
    global _THI, _ITSC, _DT, _ITSC_BASE, _ANIMAL_YEARS, _ALL_POS
    global _T_HOURS, _CUM_LOAD, _CUM_OVER
    _THI, _ITSC, _DT = thi, itsc, dt
    _ITSC_BASE = itsc_base
    _ANIMAL_YEARS = animal_years
    _ALL_POS = all_pos
    _T_HOURS, _CUM_LOAD, _CUM_OVER = t_hours, cum_load, cum_over


def _trailing(cum: np.ndarray, w: float) -> np.ndarray:
    """Trailing-``w``-hour window sum of a series with cumulative array ``cum``.

    ``cum`` is ``concat([[0], cumsum(v)])``; the window is the interval
    ``(t-w, t]``, so a gap larger than ``w`` correctly drops the pre-gap data.
    """
    left = np.searchsorted(_T_HOURS, _T_HOURS - w, side="right")
    idx = np.arange(_T_HOURS.size)
    return cum[idx + 1] - cum[left]


def _fit_series(series: np.ndarray) -> dict:
    """Median R², % converged, and raised-fit tally per lens for one series."""
    xr = (float(np.percentile(series[_ALL_POS], X_LO)),
          float(np.percentile(series[_ALL_POS], X_HI)))
    acc = {lens: {"r2": [], "conv": 0, "fail": Counter()} for lens in LENSES}
    for pos, y in _ANIMAL_YEARS:
        x = series[pos]
        for lens, fit_fn in LENSES.items():
            try:
                fit = fit_fn(x, y, x_range=xr)
            except Exception as exc:
                tally(acc[lens]["fail"], exc)
                continue
            if fit.get("converged", False):
                acc[lens]["conv"] += 1
                acc[lens]["r2"].append(fit.get("r_squared", np.nan))
    n = len(_ANIMAL_YEARS)
    out = {}
    for lens, a in acc.items():
        out[f"{lens}_r2"] = float(np.nanmedian(a["r2"])) if a["r2"] else np.nan
        out[f"{lens}_conv"] = 100.0 * a["conv"] / n
        out[f"{lens}_fail"] = sum(a["fail"].values())
        out[f"{lens}_fail_kinds"] = format_failures(a["fail"])
    return out


def sweep_tau(tau: float) -> dict:
    """At timescale ``tau``: leaky THI/ITSC and windowed TSD/TSL, every lens."""
    series = {
        "THI": history.leaky_integrate(_THI, _DT, tau, baseline=THI_STRESS_ONSET),
        "ITSC": history.leaky_integrate(_ITSC, _DT, tau, baseline=_ITSC_BASE),
        "TSD": _trailing(_CUM_OVER, tau),    # hours above threshold in trailing tau
        "TSL": _trailing(_CUM_LOAD, tau),    # degree-hours above threshold in tau
    }
    row = {"tau_h": tau}
    for name, s in series.items():
        for key, val in _fit_series(s).items():
            row[f"{name}_{key}"] = val       # e.g. THI_bs_r2, THI_hill_conv
    return row


# ─────────────────────────────────────────────────────────────
#  Data preparation
# ─────────────────────────────────────────────────────────────

def prepare(indices_csv: Path, rumen_csv: Path):
    idx = pd.read_csv(indices_csv, usecols=["datetime", "THI_NRC1971", "ITSC"])
    idx["datetime"] = pd.to_datetime(idx["datetime"], utc=True)
    idx = idx.sort_values("datetime").reset_index(drop=True)
    idx["hobo_pos"] = np.arange(len(idx))
    thi = idx["THI_NRC1971"].to_numpy()
    itsc = idx["ITSC"].to_numpy()
    dt = idx["datetime"].diff().dt.total_seconds().to_numpy() / 3600.0
    dt[0] = np.nanmedian(dt[1:])
    frac = float((thi > THI_STRESS_ONSET).mean())
    itsc_base = float(np.nanquantile(itsc, 1.0 - frac))
    log.info("ITSC baseline (matched to THI>%.0f): %.1f", THI_STRESS_ONSET, itsc_base)

    # Windowed TSD/TSL support: cumulative time-in-stress and load, plus the
    # absolute hour of each station sample (for trailing-window sums).
    t_hours = ((idx["datetime"] - idx["datetime"].iloc[0])
               .dt.total_seconds().to_numpy() / 3600.0)
    load = np.clip(thi - THI_STRESS_ONSET, 0.0, None) * dt       # degree-hours/step
    over = (thi > THI_STRESS_ONSET).astype(float) * dt           # hours/step
    cum_load = np.concatenate([[0.0], np.cumsum(load)])
    cum_over = np.concatenate([[0.0], np.cumsum(over)])

    rumen = pd.read_csv(rumen_csv)
    rumen["timestamp"] = pd.to_datetime(rumen["timestamp"], utc=True)
    rumen = rumen.dropna(subset=["body_temp"])
    rumen = rumen[rumen["timestamp"].dt.month.isin(SUMMER)].sort_values("timestamp")
    merged = pd.merge_asof(
        rumen, idx[["datetime", "hobo_pos"]],
        left_on="timestamp", right_on="datetime",
        direction="nearest", tolerance=pd.Timedelta(minutes=10))
    merged = merged.dropna(subset=["hobo_pos"])
    merged["hobo_pos"] = merged["hobo_pos"].astype(int)

    # Aggregate to hourly per animal-year: 5-min rumen temp is heavily
    # autocorrelated, and hourly keeps the Hill fits tractable. Each hour keeps
    # its mean body temp and a representative station position (memory is smooth
    # within the hour).
    merged["hour"] = merged["timestamp"].dt.floor("h")
    hourly = (merged.groupby(["animal_id", "year", "date_enter", "hour"])
              .agg(body_temp=("body_temp", "mean"),
                   hobo_pos=("hobo_pos", "median"))
              .reset_index())
    hourly["hobo_pos"] = hourly["hobo_pos"].round().astype(int)

    animal_years = []
    for _, grp in hourly.groupby(["animal_id", "year", "date_enter"]):
        if len(grp) >= MIN_READINGS:
            animal_years.append((grp["hobo_pos"].to_numpy(),
                                 grp["body_temp"].to_numpy()))
    all_pos = hourly["hobo_pos"].to_numpy()
    log.info("summer: %d 5-min -> %d hourly points; %d animal-years >= %d hours",
             len(merged), len(hourly), len(animal_years), MIN_READINGS)
    return (thi, itsc, dt, itsc_base, animal_years, all_pos,
            t_hours, cum_load, cum_over)


# ─────────────────────────────────────────────────────────────
#  Figure
# ─────────────────────────────────────────────────────────────

def plot_sweep(df: pd.DataFrame, out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from heataxis.constants import WONG
    from heataxis.viz import save_figure, setup_figure
    setup_figure()

    series = [("THI", WONG["bluish_green"], "-", "leaky THI"),
              ("ITSC", WONG["reddish_purple"], "-", "leaky ITSC"),
              ("TSL", WONG["blue"], "--", "TSL (window)"),
              ("TSD", WONG["orange"], "--", "TSD (window)")]
    lenses = [("bs", "broken-stick lens"), ("hill", "Hill / 4PL lens")]
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0), sharey="row", sharex=True)
    for j, (lens, ltitle) in enumerate(lenses):
        for col, colour, ls, label in series:
            axes[0, j].plot(df["tau_h"], df[f"{col}_{lens}_r2"], color=colour,
                            lw=1.8, ls=ls, label=label)
            axes[1, j].plot(df["tau_h"], df[f"{col}_{lens}_conv"], color=colour,
                            lw=1.8, ls=ls, label=label)
        axes[0, j].set_title(f"R²  ·  {ltitle}", fontsize=10, loc="left")
        axes[1, j].set_title(f"convergence  ·  {ltitle}", fontsize=10, loc="left")
        axes[1, j].set_xlabel("memory / window length  τ  (hours)")
    axes[0, 0].set_ylabel("median R²  (converged)")
    axes[1, 0].set_ylabel("% animal-years converged")
    axes[0, 0].legend(fontsize=8, frameon=False, ncol=2)
    fig.suptitle("Memory systems under two lenses: does the Hill model rescue "
                 "the accumulated-load fit? (HOBO summers)", fontsize=11)
    save_figure(fig, "leaky_tau_sweep_lenses", out_dir, dpi=300)
    plt.close(fig)


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                        datefmt="%H:%M:%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indices", type=Path,
                        default=Path("results/hobo_indices/hobo_indices_5min.csv.gz"))
    parser.add_argument("--rumen", type=Path,
                        default=Path("results/broken_stick_hobo/01_extract/rumen_barn.csv"))
    parser.add_argument("--out", type=Path, default=Path("results/leaky_tau_sweep"))
    parser.add_argument("--jobs", type=int, default=20)
    parser.add_argument("--tau-grid", choices=sorted(TAU_GRIDS), default="dense",
                        help="'dense' = 72 linear hours (fine peak localisation, "
                             "the historical default); 'log' = the 19-point grid "
                             "the other sweeps use, for a like-for-like comparison")
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    data = prepare(args.indices, args.rumen)
    _init_globals(*data)
    taus = TAU_GRIDS[args.tau_grid]
    log.info("sweeping tau on the '%s' grid (%d points, %.1f-%.0f h; jobs=%d) ...",
             args.tau_grid, len(taus), taus[0], taus[-1], args.jobs)

    if args.jobs > 1:
        # The initialiser form works under spawn (macOS, Windows) as well as
        # fork; under fork the initargs are inherited, so Linux pays no copy.
        with mp.Pool(processes=args.jobs, initializer=_init_globals,
                     initargs=data) as pool:
            rows = list(progress(pool.imap(sweep_tau, taus),
                                 total=len(taus), desc="leaky tau"))
    else:
        rows = [sweep_tau(t) for t in progress(taus, desc="leaky tau")]

    df = pd.DataFrame(rows)
    report(df, [f"{s}_{lens}" for s in SERIES_NAMES for lens in LENSES], log)
    df.to_csv(args.out / "leaky_tau_sweep.csv", index=False, encoding="utf-8")
    plot_sweep(df, args.out)

    log.info("done -> %s", args.out)
    for lens, ltitle in [("bs", "broken-stick"), ("hill", "Hill / 4PL")]:
        print(f"\nPeak median R² under the {ltitle} lens:")
        for col, label in [("THI", "leaky THI "), ("ITSC", "leaky ITSC"),
                           ("TSL", "TSL window"), ("TSD", "TSD window")]:
            pk = df.loc[df[f"{col}_{lens}_r2"].idxmax()]
            print(f"  {label}: R²={pk[f'{col}_{lens}_r2']:.3f} at "
                  f"tau={pk['tau_h']:.0f} h "
                  f"({pk[f'{col}_{lens}_conv']:.0f}% converged)")


if __name__ == "__main__":
    main()

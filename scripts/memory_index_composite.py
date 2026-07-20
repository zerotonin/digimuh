# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — memory_index_composite                               ║
# ║  « memory systems (low-pass THI/ITSC, TSD, TSL) vs instantaneous »║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  All on the same summer animal-years, hourly, broken-stick.      ║
# ║  Left: memory systems as curves vs their timescale tau —         ║
# ║  low-pass THI, low-pass ITSC (linear Newton-cooling trackers),   ║
# ║  and TSD / TSL (trailing-window duration / load).  Right: the    ║
# ║  instantaneous indices as categorical points.  Rows share a y    ║
# ║  axis (R² top, convergence bottom).                              ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Memory systems vs instantaneous indices, broken-stick, rumen temperature."""

from __future__ import annotations

import argparse
import logging
import multiprocessing as mp
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from heataxis.constants import THI_STRESS_ONSET
from heataxis.history import low_pass
from rerandomstats import broken_stick_fit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fitreport import format_failures, report, tally  # noqa: E402
from _sweep_common import TAU_GRID_LOG, progress  # noqa: E402

log = logging.getLogger("composite")

SUMMER = {5, 6, 7, 8, 9}
MIN_HOURS = 50
X_LO, X_HI = 5.0, 95.0
TAUS = TAU_GRID_LOG
INSTANT = ["THI_NRC1971", "BGHI", "HLI", "ETI", "THI_adj", "CCI", "ETIC",
           "ITSC", "ESI", "air_temp"]

_THI = _ITSC = _DT = _T_HOURS = _CUM_LOAD = _CUM_OVER = None
_ANIMAL_YEARS: list[tuple[np.ndarray, np.ndarray]] = []
_ALL_POS = None


def _trailing(cum, w):
    left = np.searchsorted(_T_HOURS, _T_HOURS - w, side="right")
    return cum[np.arange(_T_HOURS.size) + 1] - cum[left]


def _fit(series):
    """Median R², % converged, and a tally of fits that raised."""
    xr = (float(np.percentile(series[_ALL_POS], X_LO)),
          float(np.percentile(series[_ALL_POS], X_HI)))
    r2, conv, failures = [], 0, Counter()
    for pos, y in _ANIMAL_YEARS:
        try:
            fit = broken_stick_fit(series[pos], y, x_range=xr)
        except Exception as exc:
            tally(failures, exc)
            continue
        if fit.get("converged", False):
            conv += 1
            r2.append(fit.get("r_squared", np.nan))
    return (float(np.nanmedian(r2)) if r2 else np.nan,
            100.0 * conv / len(_ANIMAL_YEARS),
            failures)


MEMORY_SYSTEMS = ["lowpass_THI", "lowpass_ITSC", "TSD", "TSL"]


def sweep_tau(tau):
    systems = {
        "lowpass_THI": low_pass(_THI, _DT, tau),
        "lowpass_ITSC": low_pass(_ITSC, _DT, tau),
        "TSD": _trailing(_CUM_OVER, tau),
        "TSL": _trailing(_CUM_LOAD, tau),
    }
    row = {"tau_h": tau}
    for name, s in systems.items():
        r2, conv, failures = _fit(s)
        row[f"{name}_r2"], row[f"{name}_conv"] = r2, conv
        row[f"{name}_fail"] = sum(failures.values())
        row[f"{name}_fail_kinds"] = format_failures(failures)
    return row


def _init(thi, itsc, dt, t_hours, cum_load, cum_over, animal_years,
          all_pos) -> None:
    """Publish the read-only sweep state into this process's globals."""
    global _THI, _ITSC, _DT, _T_HOURS, _CUM_LOAD, _CUM_OVER, _ANIMAL_YEARS, _ALL_POS
    _THI, _ITSC, _DT = thi, itsc, dt
    _T_HOURS, _CUM_LOAD, _CUM_OVER = t_hours, cum_load, cum_over
    _ANIMAL_YEARS, _ALL_POS = animal_years, all_pos


def prepare(indices_csv, rumen_csv):
    """Build the sweep state.  Returns ``(state, instantaneous series)``."""
    cols = ["datetime", "ta", *INSTANT[:-1]]
    idx = pd.read_csv(indices_csv, usecols=lambda c: c in cols or c == "ta")
    idx = idx.rename(columns={"ta": "air_temp"})
    idx["datetime"] = pd.to_datetime(idx["datetime"], utc=True)
    idx = idx.sort_values("datetime").reset_index(drop=True)
    idx["pos"] = np.arange(len(idx))
    dt = idx["datetime"].diff().dt.total_seconds().to_numpy() / 3600.0
    dt[0] = np.nanmedian(dt[1:])
    thi = idx["THI_NRC1971"].to_numpy()
    itsc = idx["ITSC"].to_numpy()
    t_hours = ((idx["datetime"] - idx["datetime"].iloc[0])
               .dt.total_seconds().to_numpy() / 3600.0)
    load = np.clip(thi - THI_STRESS_ONSET, 0.0, None) * dt
    over = (thi > THI_STRESS_ONSET).astype(float) * dt
    cum_load = np.concatenate([[0.0], np.cumsum(load)])
    cum_over = np.concatenate([[0.0], np.cumsum(over)])

    rumen = pd.read_csv(rumen_csv)
    rumen["timestamp"] = pd.to_datetime(rumen["timestamp"], utc=True)
    rumen = rumen.dropna(subset=["body_temp"])
    rumen = rumen[rumen["timestamp"].dt.month.isin(SUMMER)].sort_values("timestamp")
    keep = ["datetime", "pos", *INSTANT]
    merged = pd.merge_asof(rumen, idx[keep], left_on="timestamp",
                           right_on="datetime", direction="nearest",
                           tolerance=pd.Timedelta(minutes=10)).dropna(subset=["pos"])
    merged["pos"] = merged["pos"].astype(int)
    merged["hour"] = merged["timestamp"].dt.floor("h")
    agg = {"body_temp": ("body_temp", "mean"), "pos": ("pos", "median")}
    agg.update({c: (c, "mean") for c in INSTANT})
    hourly = (merged.groupby(["animal_id", "year", "date_enter", "hour"])
              .agg(**agg).reset_index())
    hourly["pos"] = hourly["pos"].round().astype(int)
    animal_years = [(g["pos"].to_numpy(), g["body_temp"].to_numpy())
                    for _, g in hourly.groupby(["animal_id", "year", "date_enter"])
                    if len(g) >= MIN_HOURS]
    all_pos = hourly["pos"].to_numpy()
    log.info("%d hourly points; %d animal-years", len(hourly), len(animal_years))
    # instantaneous index series aligned to positions, for the categorical panel
    inst_series = {c: idx[c].to_numpy() for c in INSTANT}
    state = (thi, itsc, dt, t_hours, cum_load, cum_over, animal_years, all_pos)
    return state, inst_series


def instantaneous_ranking(inst_series):
    rows = []
    # Single-threaded and one full fit pass per index — worth a bar.
    for name, s in progress(inst_series.items(), total=len(inst_series),
                            desc="instantaneous"):
        r2, conv, failures = _fit(s)
        rows.append({"index": name, "r2": r2, "conv": conv,
                     "fail": sum(failures.values()),
                     "fail_kinds": format_failures(failures)})
    out = pd.DataFrame(rows).sort_values("r2", ascending=False)
    if out["fail"].sum():
        log.warning("instantaneous ranking: %d fits raised — see fail_kinds",
                    int(out["fail"].sum()))
    return out.reset_index(drop=True)


def plot(sweep, inst, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from heataxis.constants import WONG
    from heataxis.viz import save_figure, setup_figure
    setup_figure()

    mem = [("lowpass_THI", WONG["bluish_green"], "-", "low-pass THI"),
           ("lowpass_ITSC", WONG["reddish_purple"], "-", "low-pass ITSC"),
           ("TSD", WONG["orange"], "--", "TSD (window)"),
           ("TSL", WONG["blue"], "--", "TSL (window)")]
    best_mem = max(sweep[f"{m}_r2"].max() for m, *_ in mem)
    x = range(len(inst))
    icol = [WONG["vermilion"] if n == "air_temp" else WONG["blue"] for n in inst["index"]]

    fig, ax = plt.subplots(2, 2, figsize=(10.0, 7.0), sharey="row",
                           gridspec_kw={"width_ratios": [1.35, 1.0],
                                        "hspace": 0.32, "wspace": 0.05})
    for m, colour, ls, label in mem:
        ax[0, 0].plot(sweep["tau_h"], sweep[f"{m}_r2"], color=colour, lw=1.8,
                      ls=ls, marker="o", ms=2.5, label=label)
        ax[1, 0].plot(sweep["tau_h"], sweep[f"{m}_conv"], color=colour, lw=1.8,
                      ls=ls, marker="o", ms=2.5, label=label)
        for panel, metric in ((ax[0, 0], "r2"), (ax[1, 0], "conv")):
            pk = sweep.loc[sweep[f"{m}_{metric}"].idxmax()]
            panel.scatter([pk["tau_h"]], [pk[f"{m}_{metric}"]], color=colour,
                          s=95, marker="*", zorder=5, edgecolor="white", lw=0.6)
            panel.annotate(f"{pk['tau_h']:.0f} h",
                           (pk["tau_h"], pk[f"{m}_{metric}"]),
                           textcoords="offset points", xytext=(4, 5), fontsize=7,
                           color=colour, fontweight="bold")
    for a in (ax[0, 0], ax[1, 0]):
        a.set_xscale("log")
        a.set_xlabel("memory / window τ (hours)")
    ax[0, 0].set_ylabel("median R²  (broken-stick)")
    ax[0, 0].set_title("memory systems  ·  R² vs τ", fontsize=10, loc="left")
    ax[0, 0].axhline(best_mem, color=WONG["grey"], lw=0.8, ls=":")
    ax[0, 0].legend(fontsize=7.5, frameon=False, ncol=2)
    ax[1, 0].set_ylabel("% animal-years converged")
    ax[1, 0].set_title("memory systems  ·  convergence vs τ", fontsize=10, loc="left")

    ax[0, 1].axhline(best_mem, color=WONG["grey"], lw=0.8, ls=":",
                     label="best memory")
    ax[0, 1].vlines(x, 0, inst["r2"], color=icol, lw=1.2, alpha=0.5)
    ax[0, 1].scatter(x, inst["r2"], color=icol, s=34, zorder=3)
    ax[0, 1].set_title("instantaneous indices", fontsize=10, loc="left")
    ax[0, 1].set_xticks(list(x))
    ax[0, 1].set_xticklabels([])
    ax[0, 1].legend(fontsize=7.5, frameon=False, loc="upper right")
    ax[1, 1].vlines(x, 0, inst["conv"], color=icol, lw=1.2, alpha=0.5)
    ax[1, 1].scatter(x, inst["conv"], color=icol, s=34, zorder=3)
    ax[1, 1].set_xticks(list(x))
    ax[1, 1].set_xticklabels(inst["index"], rotation=90, fontsize=7.5)

    fig.suptitle("Memory systems (low-pass THI/ITSC, TSD, TSL) vs instantaneous "
                 "indices — rumen temperature, HOBO summers", fontsize=10.5)
    save_figure(fig, "memory_index_composite", out_dir, dpi=300)
    plt.close(fig)


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                        datefmt="%H:%M:%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indices", type=Path,
                        default=Path("results/hobo_indices/hobo_indices_5min.csv.gz"))
    parser.add_argument("--rumen", type=Path,
                        default=Path("results/broken_stick_hobo/01_extract/rumen_barn.csv"))
    parser.add_argument("--out", type=Path, default=Path("results/memory_index_composite"))
    parser.add_argument("--jobs", type=int, default=20)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    state, inst_series = prepare(args.indices, args.rumen)
    _init(*state)
    inst = instantaneous_ranking(inst_series)
    inst.to_csv(args.out / "instantaneous_ranking.csv", index=False, encoding="utf-8")

    log.info("sweeping memory tau (jobs=%d) ...", args.jobs)
    # The initialiser form works under spawn (macOS, Windows) as well as fork;
    # under fork the initargs are inherited, so Linux pays no extra copy.
    with mp.Pool(args.jobs, initializer=_init, initargs=state) as pool:
        sweep = pd.DataFrame(progress(pool.imap(sweep_tau, TAUS),
                                      total=len(TAUS), desc="memory tau"))
    report(sweep, MEMORY_SYSTEMS, log)
    sweep.to_csv(args.out / "memory_sweep.csv", index=False, encoding="utf-8")
    plot(sweep, inst, args.out)

    print("\nmemory-system peaks (R²):")
    for m in ("lowpass_THI", "lowpass_ITSC", "TSD", "TSL"):
        pk = sweep.loc[sweep[f"{m}_r2"].idxmax()]
        print(f"  {m:13}: R²={pk[f'{m}_r2']:.3f} at tau={pk['tau_h']:.1f} h "
              f"({pk[f'{m}_conv']:.0f}% conv)")
    print("\ninstantaneous top:")
    print(inst.head(4).to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    log.info("done -> %s", args.out)


if __name__ == "__main__":
    main()

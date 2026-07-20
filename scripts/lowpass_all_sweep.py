# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — lowpass_all_sweep                                     ║
# ║  « low-pass every index: which memory, which index wins? »       ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Applies the linear low-pass (EMA / Newton cooling) to every     ║
# ║  instantaneous index across tau, fits broken-stick on rumen      ║
# ║  temperature, and reports R² and convergence vs tau.             ║
# ║  Left: sweep curves with each index's R²/convergence optimum     ║
# ║  starred.  Right: the best values ranked as bars (shared y).     ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Low-pass sweep across all indices, with per-index optima and best-value bars."""

from __future__ import annotations

import argparse
import logging
import multiprocessing as mp
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from heataxis.history import low_pass
from rerandomstats import broken_stick_fit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _fitreport import format_failures, report, tally  # noqa: E402
from _sweep_common import TAU_GRID_LOG, progress  # noqa: E402

log = logging.getLogger("lowpass_all")

SUMMER = {5, 6, 7, 8, 9}
MIN_HOURS = 50
X_LO, X_HI = 5.0, 95.0
TAUS = TAU_GRID_LOG
INDICES = ["THI_NRC1971", "BGHI", "HLI", "ETI", "THI_adj", "CCI", "ETIC",
           "ITSC", "ESI", "air_temp"]

_SERIES: dict[str, np.ndarray] = {}
_DT = None
_ANIMAL_YEARS: list[tuple[np.ndarray, np.ndarray]] = []
_ALL_POS = None


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


def sweep_tau(tau):
    row = {"tau_h": tau}
    for name, s in _SERIES.items():
        r2, conv, failures = _fit(low_pass(s, _DT, tau))
        row[f"{name}_r2"], row[f"{name}_conv"] = r2, conv
        row[f"{name}_fail"] = sum(failures.values())
        row[f"{name}_fail_kinds"] = format_failures(failures)
    return row


def _init(series, dt, animal_years, all_pos) -> None:
    """Publish the read-only sweep state into this process's globals."""
    global _SERIES, _DT, _ANIMAL_YEARS, _ALL_POS
    _SERIES, _DT, _ANIMAL_YEARS, _ALL_POS = series, dt, animal_years, all_pos


def prepare(indices_csv, rumen_csv):
    want = {"datetime", "ta", *INDICES}
    idx = pd.read_csv(indices_csv, usecols=lambda c: c in want)
    idx = idx.rename(columns={"ta": "air_temp"})
    idx["datetime"] = pd.to_datetime(idx["datetime"], utc=True)
    idx = idx.sort_values("datetime").reset_index(drop=True)
    idx["pos"] = np.arange(len(idx))
    dt = idx["datetime"].diff().dt.total_seconds().to_numpy() / 3600.0
    dt[0] = np.nanmedian(dt[1:])
    series = {c: idx[c].to_numpy() for c in INDICES}

    rumen = pd.read_csv(rumen_csv)
    rumen["timestamp"] = pd.to_datetime(rumen["timestamp"], utc=True)
    rumen = rumen.dropna(subset=["body_temp"])
    rumen = rumen[rumen["timestamp"].dt.month.isin(SUMMER)].sort_values("timestamp")
    merged = pd.merge_asof(rumen, idx[["datetime", "pos"]], left_on="timestamp",
                           right_on="datetime", direction="nearest",
                           tolerance=pd.Timedelta(minutes=10)).dropna(subset=["pos"])
    merged["pos"] = merged["pos"].astype(int)
    merged["hour"] = merged["timestamp"].dt.floor("h")
    hourly = (merged.groupby(["animal_id", "year", "date_enter", "hour"])
              .agg(body_temp=("body_temp", "mean"), pos=("pos", "median"))
              .reset_index())
    hourly["pos"] = hourly["pos"].round().astype(int)
    animal_years = [(g["pos"].to_numpy(), g["body_temp"].to_numpy())
                    for _, g in hourly.groupby(["animal_id", "year", "date_enter"])
                    if len(g) >= MIN_HOURS]
    all_pos = hourly["pos"].to_numpy()
    log.info("%d hourly points; %d animal-years", len(hourly), len(animal_years))
    return series, dt, animal_years, all_pos


def best_table(sweep: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name in INDICES:
        r2i = sweep[f"{name}_r2"].idxmax()
        ci = sweep[f"{name}_conv"].idxmax()
        rows.append({
            "index": name,
            "best_r2": sweep.loc[r2i, f"{name}_r2"],
            "best_r2_tau": sweep.loc[r2i, "tau_h"],
            "best_conv": sweep.loc[ci, f"{name}_conv"],
            "best_conv_tau": sweep.loc[ci, "tau_h"],
        })
    return pd.DataFrame(rows).sort_values("best_r2", ascending=False).reset_index(drop=True)


def plot(sweep, best, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from heataxis.viz import save_figure, setup_figure
    setup_figure()
    cmap = plt.cm.tab10
    colour = {name: cmap(i % 10) for i, name in enumerate(INDICES)}

    fig, ax = plt.subplots(2, 2, figsize=(11.0, 7.2), sharey="row",
                           gridspec_kw={"width_ratios": [1.35, 1.0],
                                        "hspace": 0.3, "wspace": 0.06})
    for name in INDICES:
        c = colour[name]
        ax[0, 0].plot(sweep["tau_h"], sweep[f"{name}_r2"], color=c, lw=1.3,
                      alpha=0.85, label=name)
        ax[1, 0].plot(sweep["tau_h"], sweep[f"{name}_conv"], color=c, lw=1.3,
                      alpha=0.85)
        for panel, metric in ((ax[0, 0], "r2"), (ax[1, 0], "conv")):
            pk = sweep.loc[sweep[f"{name}_{metric}"].idxmax()]
            panel.scatter([pk["tau_h"]], [pk[f"{name}_{metric}"]], color=c,
                          s=85, marker="*", zorder=5, edgecolor="white", lw=0.5)
    for a in (ax[0, 0], ax[1, 0]):
        a.set_xscale("log")
        a.set_xlabel("low-pass memory τ (hours)")
    ax[0, 0].set_ylabel("median R²  (broken-stick)")
    ax[0, 0].set_title("low-pass sweep  ·  R² vs τ  (★ = optimum)", fontsize=10,
                       loc="left")
    ax[0, 0].legend(fontsize=6.5, frameon=False, ncol=2, loc="upper left")
    ax[1, 0].set_ylabel("% animal-years converged")
    ax[1, 0].set_title("low-pass sweep  ·  convergence vs τ", fontsize=10,
                       loc="left")

    x = np.arange(len(best))
    rc = [colour[n] for n in best["index"]]
    ax[0, 1].bar(x, best["best_r2"], color=rc, width=0.7)
    for xi, (_, r) in zip(x, best.iterrows(), strict=False):
        ax[0, 1].text(xi, r["best_r2"], f"{r['best_r2_tau']:.0f}h", ha="center",
                      va="bottom", fontsize=6.5)
    ax[0, 1].set_title("best R² (label = τ)", fontsize=10, loc="left")
    ax[0, 1].set_xticks(x)
    ax[0, 1].set_xticklabels([])

    order = best["index"].tolist()
    conv_vals = best.set_index("index").loc[order]
    ax[1, 1].bar(x, conv_vals["best_conv"], color=rc, width=0.7)
    for xi, n in zip(x, order, strict=False):
        ax[1, 1].text(xi, conv_vals.loc[n, "best_conv"],
                      f"{conv_vals.loc[n, 'best_conv_tau']:.0f}h", ha="center",
                      va="bottom", fontsize=6.5)
    ax[1, 1].set_title("best convergence (label = τ)", fontsize=10, loc="left")
    ax[1, 1].set_xticks(x)
    ax[1, 1].set_xticklabels(order, rotation=90, fontsize=7.5)

    fig.suptitle("Low-pass version of every index: rumen-temperature breakpoint "
                 "(HOBO summers, hourly)", fontsize=11)
    save_figure(fig, "lowpass_all_sweep", out_dir, dpi=300)
    plt.close(fig)


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                        datefmt="%H:%M:%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indices", type=Path,
                        default=Path("results/hobo_indices/hobo_indices_5min.csv.gz"))
    parser.add_argument("--rumen", type=Path,
                        default=Path("results/broken_stick_hobo/01_extract/rumen_barn.csv"))
    parser.add_argument("--out", type=Path, default=Path("results/lowpass_all_sweep"))
    parser.add_argument("--jobs", type=int, default=20)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    state = prepare(args.indices, args.rumen)
    _init(*state)
    log.info("sweeping low-pass tau over %d indices (jobs=%d) ...",
             len(INDICES), args.jobs)
    # The initialiser form works under spawn (macOS, Windows) as well as fork;
    # under fork the initargs are inherited, so Linux pays no extra copy.
    with mp.Pool(args.jobs, initializer=_init, initargs=state) as pool:
        sweep = pd.DataFrame(progress(pool.imap(sweep_tau, TAUS),
                                      total=len(TAUS), desc="low-pass tau"))
    report(sweep, INDICES, log)
    sweep.to_csv(args.out / "lowpass_all_sweep.csv", index=False, encoding="utf-8")
    best = best_table(sweep)
    best.to_csv(args.out / "lowpass_all_best.csv", index=False, encoding="utf-8")
    plot(sweep, best, args.out)
    print("\nbest low-pass per index (ranked by R²):")
    print(best.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    log.info("done -> %s", args.out)


if __name__ == "__main__":
    main()

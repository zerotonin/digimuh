# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — heat_index_tau_sweep                                  ║
# ║  « a linear low-pass temperature index: how short is the memory? »║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  The body as a linear heat capacitor: T = low-pass_tau(driver),  ║
# ║  dT/dt = (driver - T)/tau (Newton cooling — constant tau,        ║
# ║  constant loss, NO thresholding, NO effectors, monotone          ║
# ║  carryover).  Sweeps tau for three drivers (heat-balance T_eq,   ║
# ║  air temperature, THI) and reports broken-stick convergence and  ║
# ║  R² of rumen temperature vs the low-passed index.                ║
# ╚══════════════════════════════════════════════════════════════════╝
"""How short must the memory of a low-pass temperature index be?"""

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
from heatbalance_model import BalanceParams, env_equilibrium  # noqa: E402

log = logging.getLogger("tau_index")

SUMMER = {5, 6, 7, 8, 9}
MIN_HOURS = 50
X_LO, X_HI = 5.0, 95.0
TAUS = TAU_GRID_LOG

_DRIVERS: dict[str, np.ndarray] = {}
_DT = None
_ANIMAL_YEARS: list[tuple[np.ndarray, np.ndarray]] = []
_ALL_POS = None


def _init(drivers, dt, animal_years, all_pos) -> None:
    global _DRIVERS, _DT, _ANIMAL_YEARS, _ALL_POS
    _DRIVERS, _DT, _ANIMAL_YEARS, _ALL_POS = drivers, dt, animal_years, all_pos


def _fit(series: np.ndarray) -> tuple[float, float, Counter]:
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
    pct = 100.0 * conv / len(_ANIMAL_YEARS)
    return (float(np.nanmedian(r2)) if r2 else np.nan, pct, failures)


def sweep_tau(tau: float) -> dict:
    row = {"tau_h": tau}
    for name, driver in _DRIVERS.items():
        r2, pct, failures = _fit(low_pass(driver, _DT, tau))
        row[f"{name}_r2"], row[f"{name}_conv"] = r2, pct
        row[f"{name}_fail"] = sum(failures.values())
        row[f"{name}_fail_kinds"] = format_failures(failures)
    return row


def prepare(indices_csv: Path, rumen_csv: Path):
    idx = pd.read_csv(indices_csv,
                      usecols=["datetime", "ta", "rh", "u", "sr", "THI_NRC1971"])
    idx["datetime"] = pd.to_datetime(idx["datetime"], utc=True)
    idx = idx.sort_values("datetime").reset_index(drop=True)
    dt = idx["datetime"].diff().dt.total_seconds().to_numpy() / 3600.0
    dt[0] = np.nanmedian(dt[1:])
    teq = env_equilibrium(idx["ta"], idx["rh"], idx["u"], idx["sr"], BalanceParams())
    drivers = {"T_eq": np.asarray(teq), "air": idx["ta"].to_numpy(),
               "THI": idx["THI_NRC1971"].to_numpy()}
    idx["pos"] = np.arange(len(idx))

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
    log.info("%d hourly points; %d animal-years", len(hourly), len(animal_years))
    return drivers, dt, animal_years, hourly["pos"].to_numpy()


def plot(df: pd.DataFrame, out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from heataxis.constants import WONG
    from heataxis.viz import save_figure, setup_figure
    setup_figure()

    styles = [("T_eq", WONG["vermilion"], "heat-balance T_eq"),
              ("air", WONG["sky_blue"], "air temperature"),
              ("THI", WONG["black"], "THI")]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.9))
    for col, colour, label in styles:
        ax1.plot(df["tau_h"], df[f"{col}_r2"], color=colour, lw=1.8, marker="o",
                 ms=3, label=label)
        pk = df.loc[df[f"{col}_r2"].idxmax()]
        ax1.annotate(f"{pk['tau_h']:.0f} h", (pk["tau_h"], pk[f"{col}_r2"]),
                     textcoords="offset points", xytext=(3, 4), fontsize=8,
                     color=colour)
        ax2.plot(df["tau_h"], df[f"{col}_conv"], color=colour, lw=1.8, marker="o",
                 ms=3, label=label)
    for ax, ylab, title in ((ax1, "median R² (broken-stick)", "explained variance"),
                            (ax2, "% animal-years converged", "convergence")):
        ax.set_xlabel("low-pass memory τ (hours)")
        ax.set_ylabel(ylab)
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_xscale("log")
    ax1.legend(fontsize=8, frameon=False)
    fig.suptitle("A linear low-pass temperature index: rumen temperature vs "
                 "low-pass_τ(driver), broken-stick (HOBO summers)", fontsize=10.5)
    save_figure(fig, "heat_index_tau_sweep", out_dir, dpi=300)
    plt.close(fig)


def main(argv=None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                        datefmt="%H:%M:%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indices", type=Path,
                        default=Path("results/hobo_indices/hobo_indices_5min.csv.gz"))
    parser.add_argument("--rumen", type=Path,
                        default=Path("results/broken_stick_hobo/01_extract/rumen_barn.csv"))
    parser.add_argument("--out", type=Path, default=Path("results/heat_index_tau_sweep"))
    parser.add_argument("--jobs", type=int, default=20)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    state = prepare(args.indices, args.rumen)
    _init(*state)
    log.info("sweeping tau over %s h (jobs=%d) ...", TAUS, args.jobs)
    if args.jobs > 1:
        # Pass the shared state through an initialiser rather than relying on
        # fork's copy-on-write: this is the only form that also works under
        # spawn (macOS, Windows).  Under fork the initargs are still inherited,
        # so there is no extra copy on Linux.
        with mp.Pool(args.jobs, initializer=_init, initargs=state) as pool:
            rows = list(progress(pool.imap(sweep_tau, TAUS),
                                 total=len(TAUS), desc="tau sweep"))
    else:
        rows = [sweep_tau(t) for t in progress(TAUS, desc="tau sweep")]
    df = pd.DataFrame(rows)
    report(df, list(_DRIVERS), log)
    df.to_csv(args.out / "heat_index_tau_sweep.csv", index=False, encoding="utf-8")
    plot(df, args.out)
    for col in ("T_eq", "air", "THI"):
        pk = df.loc[df[f"{col}_r2"].idxmax()]
        print(f"{col:5}: peak R²={pk[f'{col}_r2']:.3f} at tau={pk['tau_h']:.1f} h "
              f"({pk[f'{col}_conv']:.0f}% converged)")
    log.info("done -> %s", args.out)


if __name__ == "__main__":
    main()

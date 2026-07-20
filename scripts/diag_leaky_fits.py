# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — diag_leaky_fits                                       ║
# ║  « why does the x-axis change fit convergence? see the fits »    ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Fits broken-stick AND Hill per animal-year for rumen temp vs a  ║
# ║  leaky-integrated THI, classifies by which converges, records    ║
# ║  WHY broken-stick fails, and plots example fits superimposed for ║
# ║  "both converge" vs "only Hill converges".                       ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Diagnostic: broken-stick vs Hill fits of rumen temperature on leaky THI.

The response shape is a y-axis property (rumen temperature is buffered ->
broken-stick).  This asks whether swapping to a leaky x-axis changes broken-
stick convergence for a *distributional* reason (zero-inflation of the leaky
index) rather than a physiological one, and shows the actual fits.
"""

from __future__ import annotations

import argparse
import logging
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from heataxis import history
from heataxis.constants import THI_STRESS_ONSET
from rerandomstats import broken_stick_fit, hill_fit

log = logging.getLogger("diag")

SUMMER = {5, 6, 7, 8, 9}
MIN_HOURS = 50
X_LO, X_HI = 5.0, 95.0
ZERO_EPS = 1e-6


def prepare(indices_csv: Path, rumen_csv: Path, tau: float):
    idx = pd.read_csv(indices_csv, usecols=["datetime", "THI_NRC1971"])
    idx["datetime"] = pd.to_datetime(idx["datetime"], utc=True)
    idx = idx.sort_values("datetime").reset_index(drop=True)
    thi = idx["THI_NRC1971"].to_numpy()
    dt = idx["datetime"].diff().dt.total_seconds().to_numpy() / 3600.0
    dt[0] = np.nanmedian(dt[1:])
    leaky = history.leaky_integrate(thi, dt, tau, baseline=THI_STRESS_ONSET)
    idx["leaky"] = leaky
    idx["pos"] = np.arange(len(idx))

    rumen = pd.read_csv(rumen_csv)
    rumen["timestamp"] = pd.to_datetime(rumen["timestamp"], utc=True)
    rumen = rumen.dropna(subset=["body_temp"])
    rumen = rumen[rumen["timestamp"].dt.month.isin(SUMMER)].sort_values("timestamp")
    merged = pd.merge_asof(
        rumen, idx[["datetime", "leaky", "pos"]],
        left_on="timestamp", right_on="datetime",
        direction="nearest", tolerance=pd.Timedelta(minutes=10))
    merged = merged.dropna(subset=["leaky"])
    merged["hour"] = merged["timestamp"].dt.floor("h")
    hourly = (merged.groupby(["animal_id", "year", "date_enter", "hour"])
              .agg(body_temp=("body_temp", "mean"), leaky=("leaky", "mean"))
              .reset_index())
    xr = (float(np.percentile(hourly["leaky"], X_LO)),
          float(np.percentile(hourly["leaky"], X_HI)))
    log.info("tau=%.0f h; %d hourly points; leaky x_range=%.2f..%.2f",
             tau, len(hourly), *xr)
    return hourly, xr


def fit_all(hourly: pd.DataFrame, xr: tuple[float, float]) -> list[dict]:
    rows = []
    for (aid, year, enter), grp in hourly.groupby(["animal_id", "year", "date_enter"]):
        if len(grp) < MIN_HOURS:
            continue
        x = grp["leaky"].to_numpy()
        y = grp["body_temp"].to_numpy()
        # Name the exception rather than tagging every failure "exception":
        # this script exists to explain *why* fits fail, so a code error must
        # stay distinguishable from a genuine non-convergence.
        try:
            bs = broken_stick_fit(x, y, x_range=xr)
        except Exception as exc:
            bs = {"converged": False, "rejected_reason": type(exc).__name__}
        try:
            hl = hill_fit(x, y, x_range=xr)
        except Exception as exc:
            hl = {"converged": False, "rejected_reason": type(exc).__name__}
        rows.append({
            "animal_id": aid, "year": year, "date_enter": enter,
            "n": len(grp), "zero_frac": float(np.mean(x <= ZERO_EPS)),
            "bs_conv": bool(bs.get("converged", False)),
            "bs_reason": bs.get("rejected_reason", None),
            "bs_r2": bs.get("r_squared", np.nan),
            "hl_conv": bool(hl.get("converged", False)),
            "hl_reason": hl.get("rejected_reason", None),
            "hl_r2": hl.get("r_squared", np.nan),
            "x": x, "y": y, "bs": bs, "hl": hl,
        })
    log.info("fit_all: %d animal-years; bs converged %d, hill converged %d",
             len(rows), sum(r["bs_conv"] for r in rows),
             sum(r["hl_conv"] for r in rows))
    return rows


def _bs_curve(bs: dict, xs: np.ndarray) -> np.ndarray:
    bp = bs["breakpoint"]
    below = bs["intercept_below"] + bs["slope_below"] * xs
    above = bs["intercept_above"] + bs["slope_above"] * xs
    return np.where(xs <= bp, below, above)


def _hill_curve(hl: dict, xs: np.ndarray) -> np.ndarray:
    xg = np.clip(xs, ZERO_EPS, None)
    return hl["y_min"] + (hl["y_max"] - hl["y_min"]) / (
        1.0 + (hl["ec50"] / xg) ** hl["hill_n"])


def plot_examples(rows: list[dict], out_dir: Path, tau: float) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from heataxis.constants import WONG
    from heataxis.viz import save_figure, setup_figure
    setup_figure()

    both = [r for r in rows if r["bs_conv"] and r["hl_conv"]]
    only_hill = [r for r in rows if r["hl_conv"] and not r["bs_conv"]]
    both.sort(key=lambda r: r["n"], reverse=True)
    only_hill.sort(key=lambda r: r["n"], reverse=True)

    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.6), sharex=True)
    groups = [("both converge", both, axes[0]),
              ("only Hill converges", only_hill, axes[1])]
    for title, group, axrow in groups:
        for ax, r in zip(axrow, group[:3], strict=False):
            ax.scatter(r["x"], r["y"], s=4, color="#bbbbbb", alpha=0.4, zorder=1)
            xs = np.linspace(r["x"].min(), r["x"].max(), 300)
            if r["hl_conv"]:
                ax.plot(xs, _hill_curve(r["hl"], xs), color=WONG["vermilion"],
                        lw=2.0, label=f"Hill (R²={r['hl_r2']:.2f})", zorder=3)
            if r["bs_conv"]:
                ax.plot(xs, _bs_curve(r["bs"], xs), color=WONG["blue"], lw=2.0,
                        label=f"broken-stick (R²={r['bs_r2']:.2f})", zorder=2)
            else:
                ax.plot([], [], color=WONG["blue"], lw=2.0,
                        label=f"broken-stick FAILED\n({r['bs_reason']})")
            ax.set_title(f"{title}\ncow {str(r['animal_id'])[-4:]} {r['year']}, "
                         f"{100 * r['zero_frac']:.0f}% at zero load",
                         fontsize=8, loc="left")
            ax.legend(fontsize=7, frameon=False, loc="upper left")
            ax.set_ylabel("rumen temp (°C)", fontsize=8)
    for ax in axes[1]:
        ax.set_xlabel(f"leaky THI (τ={tau:.0f} h)", fontsize=9)
    fig.suptitle("Rumen temperature vs leaky THI: broken-stick and Hill "
                 "superimposed", fontsize=11)
    save_figure(fig, "diag_leaky_fits", out_dir, dpi=300)
    plt.close(fig)


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                        datefmt="%H:%M:%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indices", type=Path,
                        default=Path("results/hobo_indices/hobo_indices_5min.csv.gz"))
    parser.add_argument("--rumen", type=Path,
                        default=Path("results/broken_stick_hobo/01_extract/rumen_barn.csv"))
    parser.add_argument("--out", type=Path, default=Path("results/diag_leaky_fits"))
    parser.add_argument("--tau", type=float, default=24.0)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    hourly, xr = prepare(args.indices, args.rumen, args.tau)
    rows = fit_all(hourly, xr)

    both = [r for r in rows if r["bs_conv"] and r["hl_conv"]]
    only_hill = [r for r in rows if r["hl_conv"] and not r["bs_conv"]]
    only_bs = [r for r in rows if r["bs_conv"] and not r["hl_conv"]]
    neither = [r for r in rows if not r["bs_conv"] and not r["hl_conv"]]
    n = len(rows)
    print(f"\n{n} animal-years (tau={args.tau:.0f} h):")
    print(f"  both converge      : {len(both):3d} ({100*len(both)/n:.0f}%)  "
          f"median zero-load frac {np.median([r['zero_frac'] for r in both]):.2f}")
    print(f"  only Hill converges: {len(only_hill):3d} ({100*len(only_hill)/n:.0f}%)  "
          f"median zero-load frac {np.median([r['zero_frac'] for r in only_hill]):.2f}")
    print(f"  only broken-stick  : {len(only_bs):3d}")
    print(f"  neither            : {len(neither):3d}")
    print("\nWhy broken-stick failed (only-Hill group):")
    for reason, c in Counter(r["bs_reason"] for r in only_hill).most_common():
        print(f"  {reason}: {c}")

    plot_examples(rows, args.out, args.tau)
    log.info("done -> %s", args.out)


if __name__ == "__main__":
    main()

# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — heatbalance_model                                     ║
# ║  « continuous heat balance: metabolic (circadian) + environment »║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  A reduced, always-integrating Thompson-style heat balance.      ║
# ║  y = rumen temp minus each cow's cool-day circadian rhythm (the  ║
# ║      metabolic baseline) -> the heat-stress excess.              ║
# ║  x = leaky-integrated environmental equilibrium temperature      ║
# ║      T_eq = Ta + (metabolic + solar - evaporative)/sensible.     ║
# ║  No THI threshold, no rectification -> no zero-inflation.        ║
# ║  Compares broken-stick and Hill fits of y on x.                  ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Continuous heat-balance x-axis vs circadian-detrended rumen temperature.

Fixes the rectified-leaky model: the body always integrates a signed heat
balance (gain from metabolism + sun + warm air, loss by wind-driven sensible
and humidity-capped evaporation), and the metabolic term is each cow's own
cool-day circadian rhythm.  Tests whether this removes the zero-inflation and
restores the buffered-then-breakaway shape that broken-stick expects.
"""

from __future__ import annotations

import argparse
import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from heataxis import history
from rerandomstats import broken_stick_fit, hill_fit

log = logging.getLogger("heatbalance")

SUMMER = {5, 6, 7, 8, 9}
MIN_HOURS = 50
X_LO, X_HI = 5.0, 95.0


@dataclass(frozen=True)
class BalanceParams:
    """Reduced heat-balance coefficients (first-pass, Thompson-scaled)."""
    metabolic: float = 100.0   # W/m² endogenous heat (absorbed as a constant)
    alpha_sr: float = 0.15     # solar absorptivity x geometry
    g0: float = 8.0            # sensible+radiative conductance (W/m²/K)
    ku: float = 0.6            # wind enhancement of sensible loss (per √(m/s))
    e0: float = 40.0           # evaporative capacity per kPa VPD (W/m²/kPa)


def env_equilibrium(ta, rh, u, sr, p: BalanceParams) -> np.ndarray:
    """Equilibrium body temperature the environment drives (°C).

    Gain = metabolic + solar; loss = wind-enhanced sensible + humidity-capped
    evaporation (∝ vapour-pressure deficit).  Hot/humid/still -> high T_eq;
    cool/dry/windy -> low T_eq.  Always defined, never rectified.
    """
    ta = np.asarray(ta, float)
    es = 0.6108 * np.exp(17.27 * ta / (ta + 237.3))     # kPa (Tetens)
    vpd = es * (1.0 - np.asarray(rh, float) / 100.0)     # evaporative potential
    g_sens = p.g0 * (1.0 + p.ku * np.sqrt(np.clip(u, 0.0, None)))
    e_evap = p.e0 * vpd
    return ta + (p.metabolic + p.alpha_sr * np.asarray(sr, float) - e_evap) / g_sens


def load_circadian_baseline(path: Path) -> dict:
    """Cool-day rumen rhythm per (animal_id, year, hour) — the metabolic base."""
    c = pd.read_csv(path)
    cool = c[c["day_type"] == "cool"]
    return {(int(r.animal_id), int(r.year), int(r.hour)): float(r.body_temp_mean)
            for r in cool.itertuples()}


def prepare(indices_csv, rumen_csv, circadian_csv, tau, params):
    idx = pd.read_csv(indices_csv, usecols=["datetime", "ta", "rh", "u", "sr"])
    idx["datetime"] = pd.to_datetime(idx["datetime"], utc=True)
    idx = idx.sort_values("datetime").reset_index(drop=True)
    teq = env_equilibrium(idx["ta"], idx["rh"], idx["u"], idx["sr"], params)
    dt = idx["datetime"].diff().dt.total_seconds().to_numpy() / 3600.0
    dt[0] = np.nanmedian(dt[1:])
    teq_neutral = float(np.nanpercentile(teq, 10))       # thermoneutral reference
    load = history.leaky_integrate(teq, dt, tau, baseline=teq_neutral, rectify=False)
    idx["load"] = load
    idx["pos"] = np.arange(len(idx))
    log.info("T_eq range %.1f..%.1f °C (neutral %.1f); load range %.1f..%.1f",
             float(teq.min()), float(teq.max()), teq_neutral,
             float(load.min()), float(load.max()))

    base = load_circadian_baseline(circadian_csv)
    rumen = pd.read_csv(rumen_csv)
    rumen["timestamp"] = pd.to_datetime(rumen["timestamp"], utc=True)
    rumen = rumen.dropna(subset=["body_temp"])
    rumen = rumen[rumen["timestamp"].dt.month.isin(SUMMER)].sort_values("timestamp")
    merged = pd.merge_asof(
        rumen, idx[["datetime", "load", "pos"]], left_on="timestamp",
        right_on="datetime", direction="nearest", tolerance=pd.Timedelta(minutes=10))
    merged = merged.dropna(subset=["load"])
    merged["hour"] = merged["timestamp"].dt.floor("h")
    merged["hod"] = merged["timestamp"].dt.hour
    hourly = (merged.groupby(["animal_id", "year", "date_enter", "hour"])
              .agg(body_temp=("body_temp", "mean"), load=("load", "mean"),
                   hod=("hod", "first"))
              .reset_index())
    # Subtract each cow's cool-day circadian baseline -> heat-stress excess.
    hourly["base"] = [base.get((int(a), int(y), int(h)), np.nan)
                      for a, y, h in zip(hourly["animal_id"], hourly["year"],
                                         hourly["hod"], strict=False)]
    hourly = hourly.dropna(subset=["base"])
    hourly["excess"] = hourly["body_temp"] - hourly["base"]
    xr = (float(np.percentile(hourly["load"], X_LO)),
          float(np.percentile(hourly["load"], X_HI)))
    log.info("%d hourly points with circadian baseline; load x_range %.1f..%.1f",
             len(hourly), *xr)
    return hourly, xr


def fit_all(hourly, xr):
    rows = []
    lo = xr[0] + 0.02 * (xr[1] - xr[0])
    for (aid, year, enter), grp in hourly.groupby(["animal_id", "year", "date_enter"]):
        if len(grp) < MIN_HOURS:
            continue
        x = grp["load"].to_numpy()
        y = grp["excess"].to_numpy()
        # Both lenses record why they failed, and name the exception when one
        # raised — otherwise a bug is indistinguishable from a genuine
        # non-convergence in the summary tables downstream.
        try:
            bs = broken_stick_fit(x, y, x_range=xr)
        except Exception as exc:
            bs = {"converged": False, "rejected_reason": type(exc).__name__}
        try:
            hl = hill_fit(x, y, x_range=xr)
        except Exception as exc:
            hl = {"converged": False, "rejected_reason": type(exc).__name__}
        rows.append({
            "animal_id": aid, "year": year, "n": len(grp),
            "low_frac": float(np.mean(x <= lo)),
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


def _bs_curve(bs, xs):
    return np.where(xs <= bs["breakpoint"],
                    bs["intercept_below"] + bs["slope_below"] * xs,
                    bs["intercept_above"] + bs["slope_above"] * xs)


def _hill_curve(hl, xs):
    xg = np.clip(xs - xs.min() + 1e-6, 1e-6, None)
    return hl["y_min"] + (hl["y_max"] - hl["y_min"]) / (
        1.0 + (hl["ec50"] / xg) ** hl["hill_n"])


def plot_examples(rows, out_dir, tau):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from heataxis.constants import WONG
    from heataxis.viz import save_figure, setup_figure
    setup_figure()

    both = sorted([r for r in rows if r["bs_conv"] and r["hl_conv"]],
                  key=lambda r: r["n"], reverse=True)
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.6))
    for ax, r in zip(axes.ravel(), both[:6], strict=False):
        ax.scatter(r["x"], r["y"], s=4, color="#bbbbbb", alpha=0.4, zorder=1)
        xs = np.linspace(r["x"].min(), r["x"].max(), 300)
        if r["bs_conv"]:
            ax.plot(xs, _bs_curve(r["bs"], xs), color=WONG["blue"], lw=2.0,
                    label=f"broken-stick (R²={r['bs_r2']:.2f})", zorder=3)
        if r["hl_conv"]:
            ax.plot(xs, _hill_curve(r["hl"], xs), color=WONG["vermilion"], lw=2.0,
                    label=f"Hill (R²={r['hl_r2']:.2f})", zorder=2)
        ax.axhline(0, color="#dddddd", lw=0.8, zorder=0)
        ax.set_title(f"cow {str(r['animal_id'])[-4:]} {r['year']}", fontsize=8,
                     loc="left")
        ax.legend(fontsize=7, frameon=False, loc="upper left")
        ax.set_xlabel(f"integrated env. load (τ={tau:.0f} h)", fontsize=8)
        ax.set_ylabel("rumen temp − circadian (°C)", fontsize=8)
    fig.suptitle("Heat-stress excess vs continuous environmental load "
                 "(circadian-detrended, no rectification)", fontsize=11)
    save_figure(fig, "heatbalance_fits", out_dir, dpi=300)
    plt.close(fig)


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                        datefmt="%H:%M:%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indices", type=Path,
                        default=Path("results/hobo_indices/hobo_indices_5min.csv.gz"))
    parser.add_argument("--rumen", type=Path,
                        default=Path("results/broken_stick_hobo/01_extract/rumen_barn.csv"))
    parser.add_argument("--circadian", type=Path,
                        default=Path("results/broken_stick_hobo/03_temporal/circadian_null_model.csv"))
    parser.add_argument("--out", type=Path, default=Path("results/heatbalance_model"))
    parser.add_argument("--tau", type=float, default=24.0)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    hourly, xr = prepare(args.indices, args.rumen, args.circadian, args.tau,
                         BalanceParams())
    rows = fit_all(hourly, xr)
    n = len(rows)
    both = [r for r in rows if r["bs_conv"] and r["hl_conv"]]
    only_hill = [r for r in rows if r["hl_conv"] and not r["bs_conv"]]
    only_bs = [r for r in rows if r["bs_conv"] and not r["hl_conv"]]
    print(f"\n{n} animal-years (heat-balance load, tau={args.tau:.0f} h):")
    print(f"  both converge      : {len(both):3d} ({100*len(both)/n:.0f}%)")
    print(f"  only Hill          : {len(only_hill):3d} ({100*len(only_hill)/n:.0f}%)")
    print(f"  only broken-stick  : {len(only_bs):3d}")
    if both:
        print(f"  broken-stick median R² (both): {np.median([r['bs_r2'] for r in both]):.3f}")
    if only_hill:
        print("  only-Hill broken-stick failure reasons:")
        for reason, c in Counter(r["bs_reason"] for r in only_hill).most_common(4):
            print(f"    {reason}: {c}")
    plot_examples(rows, args.out, args.tau)
    log.info("done -> %s", args.out)


if __name__ == "__main__":
    main()

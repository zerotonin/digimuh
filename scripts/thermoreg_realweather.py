# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — thermoreg_realweather                                 ║
# ║  « drive the validated ODE with real HOBO weather »             ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Takes the synthetic-validated one-state heat-balance ODE        ║
# ║  (thermoreg_sim), drives it with the real HOBO weather over a    ║
# ║  hot window for several cow-summers, and overlays the simulated  ║
# ║  rumen temperature on the observed readings.  Each cow uses its  ║
# ║  own Frontiers cool-day circadian.                              ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Overlay ODE-simulated rumen temperature on observed, driven by real weather."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from thermoreg_sim import (  # noqa: E402
    Cow,
    load_cold_day_setpoints,
    simulate,
)

log = logging.getLogger("realweather")

WINDOW_DAYS = 14
GRID_MIN = 10.0        # regular grid for integration (min)


def load_weather(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, usecols=["datetime", "ta", "rh", "u", "sr",
                                    "THI_NRC1971"])
    df["datetime"] = pd.to_datetime(df["datetime"], utc=True)
    return df.sort_values("datetime").reset_index(drop=True)


def window_weather(weather: pd.DataFrame, start, end) -> dict:
    """Regular-grid HOBO weather over [start, end] for the integrator."""
    grid = pd.date_range(start, end, freq=f"{int(GRID_MIN)}min", tz="UTC")
    sub = weather[(weather["datetime"] >= start) & (weather["datetime"] <= end)]
    out = {}
    for col in ("ta", "rh", "u", "sr"):
        out[col] = np.interp(grid.astype("int64"),
                             sub["datetime"].astype("int64"), sub[col].to_numpy())
    t_h = (grid - grid[0]).total_seconds().to_numpy() / 3600.0
    return {"t_h": t_h, "hod": (grid.hour + grid.minute / 60.0).to_numpy(),
            "day": t_h / 24.0, "ta": out["ta"], "rh": out["rh"],
            "u": out["u"], "sr": out["sr"], "grid": grid}


def select_cow_summers(rumen: pd.DataFrame, weather: pd.DataFrame,
                       setpoints: dict, n: int) -> list[dict]:
    """Cow-summers with good coverage; window = ±7 d around their hottest day."""
    thi_daily = (weather.set_index("datetime")["THI_NRC1971"]
                 .resample("1D").max())
    counts = rumen.groupby(["animal_id", "year"]).size().sort_values(ascending=False)
    picks = []
    for (aid, year), _ in counts.items():
        if (int(aid), int(year)) not in setpoints:
            continue
        grp = rumen[(rumen["animal_id"] == aid) & (rumen["year"] == year)]
        t0, t1 = grp["timestamp"].min(), grp["timestamp"].max()
        yr_thi = thi_daily[(thi_daily.index >= t0) & (thi_daily.index <= t1)]
        if yr_thi.empty or not np.isfinite(yr_thi.max()):
            continue
        hot_day = yr_thi.idxmax()
        start = max(t0, hot_day - pd.Timedelta(days=WINDOW_DAYS // 2))
        end = min(t1, start + pd.Timedelta(days=WINDOW_DAYS))
        picks.append({"aid": int(aid), "year": int(year), "grp": grp,
                      "start": start, "end": end})
        if len(picks) >= n:
            break
    return picks


def plot(runs: list[dict], out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from heataxis.constants import WONG
    from heataxis.viz import save_figure, setup_figure
    setup_figure()

    fig, axes = plt.subplots(2, 3, figsize=(12.5, 6.8), sharey=True)
    for ax, r in zip(axes.ravel(), runs, strict=False):
        ax.scatter(r["obs_day"], r["obs_t"], s=4, color="#bbbbbb", alpha=0.5,
                   zorder=1, label="observed")
        ax.plot(r["day"], r["sim"], color=WONG["vermilion"], lw=1.3, zorder=3,
                label="simulated (real weather)")
        ax2 = ax.twinx()
        ax2.plot(r["day"], r["ta"], color="#cfe3f2", lw=0.7, zorder=0)
        ax2.set_ylabel("air (°C)", color="#a8c6df", fontsize=7)
        ax2.tick_params(labelsize=6)
        ax.set_zorder(ax2.get_zorder() + 1)
        ax.patch.set_visible(False)
        ax.set_title(f"cow {str(r['aid'])[-4:]} {r['year']}  "
                     f"(r={r['corr']:.2f})", fontsize=9, loc="left")
        ax.set_xlabel("day in window", fontsize=8)
    axes[0, 0].set_ylabel("rumen T (°C)")
    axes[1, 0].set_ylabel("rumen T (°C)")
    axes[0, 0].legend(fontsize=7, frameon=False, loc="upper left")
    fig.suptitle("ODE driven by real HOBO weather vs observed rumen temperature "
                 "(per-animal Frontiers circadian)", fontsize=11)
    save_figure(fig, "thermoreg_realweather", out_dir, dpi=300)
    plt.close(fig)


def main(argv=None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                        datefmt="%H:%M:%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indices", type=Path,
                        default=Path("results/hobo_indices/hobo_indices_5min.csv.gz"))
    parser.add_argument("--rumen", type=Path,
                        default=Path("results/broken_stick_hobo/01_extract/rumen_barn.csv"))
    parser.add_argument("--circadian", type=Path,
                        default=Path("results/broken_stick_hobo/03_temporal/circadian_null_model.csv"))
    parser.add_argument("--out", type=Path, default=Path("results/thermoreg_realweather"))
    parser.add_argument("--n", type=int, default=6)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    weather = load_weather(args.indices)
    _, per_animal = load_cold_day_setpoints(args.circadian)
    rumen = pd.read_csv(args.rumen)
    rumen["timestamp"] = pd.to_datetime(rumen["timestamp"], utc=True)
    rumen = rumen.dropna(subset=["body_temp"])

    cow = Cow()
    runs = []
    for pick in select_cow_summers(rumen, weather, per_animal, args.n):
        w = window_weather(weather, pick["start"], pick["end"])
        setpoint = per_animal[(pick["aid"], pick["year"])]
        sim, _ = simulate(w, cow, setpoint)

        obs = pick["grp"]
        obs = obs[(obs["timestamp"] >= pick["start"]) & (obs["timestamp"] <= pick["end"])]
        obs_hourly = (obs.set_index("timestamp")["body_temp"]
                      .resample("1h").mean().dropna())
        obs_day = (obs_hourly.index - w["grid"][0]).total_seconds().to_numpy() / 86400.0
        # correlation of simulated vs observed at the observed hours
        sim_at_obs = np.interp(obs_day, w["day"], sim)
        corr = float(np.corrcoef(sim_at_obs, obs_hourly.to_numpy())[0, 1]) \
            if len(obs_hourly) > 5 else np.nan
        runs.append({**pick, "day": w["day"], "sim": sim, "ta": w["ta"],
                     "obs_day": obs_day, "obs_t": obs_hourly.to_numpy(),
                     "corr": corr})
        log.info("cow %d %d: window %s..%s, r=%.2f", pick["aid"], pick["year"],
                 pick["start"].date(), pick["end"].date(), corr)

    plot(runs, args.out)
    med = np.nanmedian([r["corr"] for r in runs])
    log.info("done -> %s  (median r = %.2f)", args.out, med)


if __name__ == "__main__":
    main()

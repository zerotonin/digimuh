# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — thermoreg_sim                                         ║
# ║  « a one-state heat-capacitor cow: survive winter and summer »   ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  A simple, physically-sound Thompson-style ODE for core body     ║
# ║  temperature.  The environment is ALWAYS coupled through the     ║
# ║  heat-exchange gradient — never thresholded.  Regulated          ║
# ║  conductance (vasomotor), evaporation, and shivering defend the  ║
# ║  setpoint; heat capacity C buffers.  Validated on SYNTHETIC      ║
# ║  winter and summer before it ever sees real data.               ║
# ╚══════════════════════════════════════════════════════════════════╝
"""One-state body-heat ODE, tested on synthetic winter and summer.

    C dT_b/dt = HE(t) + a*SR - k(T_b,u)(T_b - T_a) - E(T_b,RH,u) + shiver

The question: does the integrator stay stable and hold T_b near setpoint with
its circadian rhythm across a cold winter and a hot summer (breaking away only
in an extreme humid heat wave)?
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger("thermoreg")


@dataclass(frozen=True)
class CircadianSetpoint:
    """Regulated thermoneutral rumen rhythm — a single 24 h harmonic.

    The cool-day rumen temperature *is* the metabolic (fermentation) baseline the
    animal defends; a sinusoid smooths its gappy hourly means.  One per herd (the
    mean) now; one per animal-year later — same interface either way.
    """
    mesor: float
    amplitude: float
    acrophase_h: float

    def at(self, hod: np.ndarray | float) -> np.ndarray | float:
        return self.mesor + self.amplitude * np.cos(
            2 * np.pi * (np.asarray(hod) - self.acrophase_h) / 24)


def fit_setpoint(hours, values) -> CircadianSetpoint:
    """Closed-form single-harmonic (DFT) fit of a 24 h rhythm."""
    h = np.asarray(hours, dtype=float)
    v = np.asarray(values, dtype=float)
    m = np.isfinite(v)
    h, v = h[m], v[m]
    w = 2 * np.pi / 24
    cos_c, sin_c = np.mean(v * np.cos(w * h)), np.mean(v * np.sin(w * h))
    amp = 2 * np.hypot(cos_c, sin_c)
    acro = (np.arctan2(sin_c, cos_c) * 24 / (2 * np.pi)) % 24
    return CircadianSetpoint(float(v.mean()), float(amp), float(acro))


def load_cold_day_setpoints(path: Path):
    """(mean, per-animal) cold-day rumen setpoints from the circadian null model."""
    cool = pd.read_csv(path).query("day_type == 'cool'")
    hourly_mean = cool.groupby("hour")["body_temp_mean"].mean()
    mean = fit_setpoint(hourly_mean.index.to_numpy(), hourly_mean.to_numpy())
    per_animal = {
        (int(aid), int(yr)): fit_setpoint(g["hour"].to_numpy(),
                                          g["body_temp_mean"].to_numpy())
        for (aid, yr), g in cool.groupby(["animal_id", "year"])
        if g["hour"].nunique() >= 6
    }
    return mean, per_animal


@dataclass(frozen=True)
class Cow:
    """Whole-body thermal parameters (SI: W, J/K, °C).

    Approach B: the circadian is the metabolic (fermentation) HEAT rhythm — it
    peaks in the evening and actively drives the rumen peak.  The effectors
    defend a constant setpoint (the circadian mesor); vasomotor control holds the
    core across the thermoneutral zone, shivering engages below the lower critical
    temp, and humidity-capped evaporation saturates in a humid heat wave ->
    breakaway.
    """
    C: float = 9.0e5             # body heat capacity (J/K)  -> tau = C/k
    he_mean: float = 1350.0      # mean metabolic heat (W)
    he_gain: float = 2200.0      # metabolic-heat rhythm per °C of circadian amp (W/°C)
    k_min: float = 12.0          # vasoconstricted conductance (W/K, cold)
    k_max: float = 130.0         # vasodilated conductance (W/K, hot)
    k_w: float = 0.30            # vasomotor gain width (°C)
    wind_k: float = 0.15         # wind enhancement of sensible loss (per √(m/s))
    e_max: float = 2000.0        # max evaporative heat loss (W)
    e_w: float = 0.25            # evaporative activation width (°C)
    alpha_sr: float = 0.30       # solar gain per unit SR (W per W/m²)
    lct_offset: float = 0.5      # °C below setpoint where shivering engages
    shiver: float = 1200.0       # cold thermogenesis gain (W/K)


def _sigmoid(z: np.ndarray | float) -> np.ndarray | float:
    return 1.0 / (1.0 + np.exp(-np.clip(z, -50, 50)))


def conductance(t_b: float, t_set: float, u: float, cow: Cow) -> float:
    """Vasomotor conductance about the setpoint: constrict cool, dilate warm."""
    vaso = cow.k_min + (cow.k_max - cow.k_min) * _sigmoid((t_b - t_set) / cow.k_w)
    return vaso * (1.0 + cow.wind_k * np.sqrt(max(u, 0.0)))


def evaporation(t_b: float, t_set: float, rh: float, u: float, cow: Cow) -> float:
    """Evaporative loss: activates above setpoint, humidity-capped, wind-aided."""
    activation = _sigmoid((t_b - t_set) / cow.e_w)
    humidity = max(1.0 - rh / 100.0, 0.0)
    return cow.e_max * activation * humidity * (1.0 + cow.wind_k * np.sqrt(max(u, 0.0)))


def d_tb(t_b: float, t_set: float, he: float, ta: float, rh: float, u: float,
         sr: float, cow: Cow) -> float:
    """Heat balance. ``he`` is the metabolic heat (already carries the rhythm);
    ``t_set`` is the constant setpoint the effectors defend."""
    he_total = he + cow.shiver * max((t_set - cow.lct_offset) - t_b, 0.0)
    q_sens = conductance(t_b, t_set, u, cow) * (t_b - ta)
    q_evap = evaporation(t_b, t_set, rh, u, cow)
    q_sol = cow.alpha_sr * sr
    return (he_total + q_sol - q_sens - q_evap) / cow.C


def simulate(weather: dict, cow: Cow,
             setpoint: CircadianSetpoint) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the ODE. Circadian drives the metabolic HEAT rhythm (approach B);
    effectors defend the constant setpoint mesor. Returns (T_b, cool-day target)."""
    t_h, ta, rh, u, sr = (weather[k] for k in ("t_h", "ta", "rh", "u", "sr"))
    hod = weather["hod"]
    dt_s = (t_h[1] - t_h[0]) * 3600.0
    # Metabolic heat with the fermentation circadian (phase & amplitude from data).
    he_series = cow.he_mean + cow.he_gain * setpoint.amplitude * np.cos(
        2 * np.pi * (hod - setpoint.acrophase_h) / 24)
    t_set = setpoint.mesor                       # constant effector setpoint
    t_target = np.asarray(setpoint.at(hod), dtype=float)   # cool-day rhythm (overlay)
    t_b = np.empty_like(ta)
    t_b[0] = t_set
    for i in range(1, len(ta)):
        t_b[i] = t_b[i - 1] + dt_s * d_tb(
            t_b[i - 1], t_set, he_series[i - 1], ta[i - 1], rh[i - 1], u[i - 1],
            sr[i - 1], cow)
    return t_b, t_target


# ─────────────────────────────────────────────────────────────
#  Synthetic weather
# ─────────────────────────────────────────────────────────────

def synth_weather(season: str, days: int = 24, dt_min: float = 5.0) -> dict:
    t_h = np.arange(0.0, days * 24.0, dt_min / 60.0)
    hod = t_h % 24.0
    day = t_h / 24.0
    diurnal = np.cos(2 * np.pi * (hod - 15.0) / 24.0)   # peak ~15:00
    daylight = np.clip(np.sin(np.pi * (hod - 6.0) / 12.0), 0.0, None)
    if season == "winter":
        ta = -2.0 + 5.0 * diurnal                        # -7 .. +3 °C
        rh = np.clip(85.0 - 8.0 * diurnal, 40, 98)
        u = np.full_like(t_h, 1.5)
        sr = 250.0 * daylight                             # short, weak sun
    else:  # summer with a hot, HUMID heat wave in the middle (evaporation fails)
        wave = np.where((day >= 8) & (day <= 15),
                        0.5 * (1 - np.cos(2 * np.pi * (day - 8) / 7)), 0.0)
        ta = 20.0 + 12.0 * wave + (7.0 - 3.0 * wave) * diurnal
        rh = np.clip(58.0 + 20.0 * wave - 8.0 * diurnal, 30, 95)
        u = np.full_like(t_h, 0.6)
        sr = 800.0 * daylight * (0.7 + 0.3 * wave)
    return {"t_h": t_h, "hod": hod, "day": day, "ta": ta, "rh": rh, "u": u, "sr": sr}


# ─────────────────────────────────────────────────────────────
#  Figure
# ─────────────────────────────────────────────────────────────

def plot(results: dict, setpoint: CircadianSetpoint, out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from heataxis.constants import WONG
    from heataxis.viz import save_figure, setup_figure
    setup_figure()

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 6.6),
                             gridspec_kw={"width_ratios": [2.2, 1.0]})
    for row, season in enumerate(("winter", "summer")):
        w, tb, tset = (results[season][k] for k in ("w", "tb", "tset"))
        colour = WONG["blue"] if season == "winter" else WONG["vermilion"]
        axenv, axzoom = axes[row]
        ax2 = axenv.twinx()
        ax2.plot(w["day"], w["ta"], color="#cccccc", lw=0.8)
        ax2.set_ylabel("air temp (°C)", color="#999999", fontsize=8)
        axenv.plot(w["day"], tb, color=colour, lw=1.2, label="simulated rumen T")
        axenv.plot(w["day"], tset, color=WONG["black"], lw=0.7, ls="--",
                   label="circadian setpoint")
        axenv.set_ylabel("rumen T (°C)")
        axenv.set_title(f"{season}: rumen temperature vs air "
                        f"(T range {tb.min():.2f}–{tb.max():.2f} °C)",
                        fontsize=10, loc="left")
        axenv.legend(fontsize=7, frameon=False, loc="upper left")
        axenv.set_zorder(ax2.get_zorder() + 1)
        axenv.patch.set_visible(False)
        mid = (w["day"] >= 5) & (w["day"] <= 8)
        axzoom.plot(w["day"][mid], tb[mid], color=colour, lw=1.2)
        axzoom.plot(w["day"][mid], tset[mid], color=WONG["black"], lw=0.7, ls="--")
        axzoom.set_title("3-day zoom (circadian)", fontsize=9, loc="left")
        axzoom.set_ylabel("rumen T (°C)", fontsize=8)
    for ax in axes[1]:
        ax.set_xlabel("day")
    fig.suptitle("Circadian setpoint wired in: rumen T tracks the cool-day rhythm "
                 f"(mesor {setpoint.mesor:.2f} °C, amp {setpoint.amplitude:.2f}, "
                 f"peak {setpoint.acrophase_h:.1f} h) and breaks away in heat",
                 fontsize=10.5)
    save_figure(fig, "thermoreg_sim", out_dir, dpi=300)
    plt.close(fig)


def _season_stats(tb: np.ndarray, w: dict, setpoint: CircadianSetpoint) -> str:
    settled = tb[w["day"] >= 2]
    return (f"T {settled.min():.2f}–{settled.max():.2f} °C "
            f"(mean {settled.mean():.2f}, setpoint mesor {setpoint.mesor:.2f})")


def main(argv=None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                        datefmt="%H:%M:%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("results/thermoreg_sim"))
    parser.add_argument("--circadian", type=Path,
                        default=Path("results/broken_stick_hobo/03_temporal/circadian_null_model.csv"))
    parser.add_argument("--days", type=int, default=24)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    setpoint, per_animal = load_cold_day_setpoints(args.circadian)
    log.info("mean cold-day setpoint: mesor %.2f °C, amp %.3f, peak %.1f h "
             "(%d per-animal curves available)", setpoint.mesor,
             setpoint.amplitude, setpoint.acrophase_h, len(per_animal))

    cow = Cow()
    results = {}
    for season in ("winter", "summer"):
        w = synth_weather(season, days=args.days)
        tb, tset = simulate(w, cow, setpoint)
        results[season] = {"w": w, "tb": tb, "tset": tset}
        log.info("%-7s: %s", season, _season_stats(tb, w, setpoint))
    plot(results, setpoint, args.out)
    log.info("done -> %s", args.out)


if __name__ == "__main__":
    main()

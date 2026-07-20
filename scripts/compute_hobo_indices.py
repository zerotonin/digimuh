# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — compute_hobo_indices                                  ║
# ║  « every thermal index, on the HOBO wind+solar weather station » ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Pulls hobo_weather (Ta, RH, dew point, solar, wind) from cow.db, ║
# ║  cleans it, derives Tw / Tbg, and computes the full heataxis      ║
# ║  index set — the ones that need wind/solar included, which the    ║
# ║  smaXtec climate could not support.                              ║
# ║                                                                  ║
# ║  Outputs a tidy 5-min table, a daily aggregate, a per-index       ║
# ║  summary, and a sample-week figure into results/hobo_indices/.    ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Compute every heataxis thermal index on the HOBO weather station.

The HOBO station carries wind speed and solar radiation, so it supports the
wind- and radiation-aware indices (ETI, HLI, CCI, ETIC, THI_adj, BGHI, ITSC,
ESI) that the temperature/humidity-only smaXtec climate cannot.  Run:

    PYTHONPATH=<path to heataxis>/src \\
        python scripts/compute_hobo_indices.py
    ... --out results/hobo_indices  --db /path/to/cow.db
"""

from __future__ import annotations

import argparse
import logging
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from heataxis import history
from heataxis import indices as idx
from heataxis.constants import THI_STRESS_ONSET

from digimuh.config import load_config

log = logging.getLogger("hobo_indices")

# HOBO channel -> tidy name (station node ids are opaque in the raw schema).
HOBO_COLUMNS = {
    "datetime": "datetime",
    "21141733_1_temperature": "ta",
    "21141733_2_rh": "rh",
    "21141733_3_dew_point": "tdp",
    "21141735_1_solar_radiation": "sr",
    "21141734_1_wind_speed": "u",
}
SOLAR_CAP = 1300.0        # W/m²; a couple of sensor spikes exceed this


# ─────────────────────────────────────────────────────────────
#  Load and clean
# ─────────────────────────────────────────────────────────────

def load_hobo(db_path: Path) -> pd.DataFrame:
    """Read hobo_weather into a tidy, cleaned frame (Ta, RH, Tdp, SR, u)."""
    cols = ", ".join(f'"{c}"' for c in HOBO_COLUMNS)
    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql(f"select {cols} from hobo_weather", conn)
    df = df.rename(columns=HOBO_COLUMNS)
    df["datetime"] = pd.to_datetime(df["datetime"], utc=True)

    before = len(df)
    df = df.dropna(subset=["ta", "rh"])
    df["rh"] = df["rh"].clip(1.0, 100.0)
    df["sr"] = df["sr"].clip(lower=0.0, upper=SOLAR_CAP)
    df["u"] = df["u"].clip(lower=0.0)
    # Fill the few missing wind/solar with the column median (indoor still air).
    for col in ("sr", "u"):
        df[col] = df[col].fillna(df[col].median())
    df = df.sort_values("datetime").reset_index(drop=True)
    log.info("loaded %d rows (%d dropped for missing Ta/RH); %s to %s",
             len(df), before - len(df),
             df["datetime"].min().date(), df["datetime"].max().date())
    return df


# ─────────────────────────────────────────────────────────────
#  Compute every index
# ─────────────────────────────────────────────────────────────

def compute_indices(df: pd.DataFrame) -> pd.DataFrame:
    """Derive Tw/Tbg and compute the full heataxis index set on the frame."""
    ta, rh, u, sr = df["ta"].to_numpy(), df["rh"].to_numpy(), \
        df["u"].to_numpy(), df["sr"].to_numpy()
    tdp = df["tdp"].to_numpy()
    tw = idx.wet_bulb(ta, rh)
    tbg = idx.tbg_estimate(ta, sr, u)

    out = df.copy()
    out["Tw"] = tw
    out["Tbg"] = tbg

    # Temperature-humidity family (incl. the coefficient variants).
    out["THI_NRC1971"] = idx.thi_nrc1971(ta, rh)
    out["THI_wetbulb"] = idx.thi_wetbulb(ta, tw=tw)
    out["THI_dewpoint"] = idx.thi_dewpoint(ta, tdp=tdp)
    for name, vals in idx.thi_variants(ta, rh).items():
        out[f"THIvar_{name}"] = vals

    # Radiation-, wind-, and full-microclimate indices.
    out["BGHI"] = idx.bghi(tbg, tdp=tdp)
    out["HLI"] = idx.hli_gaughan(tbg, rh, u)
    out["ETI"] = idx.eti_baeta(ta, rh, u)
    out["THI_adj"] = idx.thi_adj_mader(ta, rh, u, sr)
    out["CCI"] = idx.cci(ta, rh, u, sr)
    out["ETIC"] = idx.etic(ta, rh, u, sr)
    out["ITSC"] = idx.itsc(ta, u, rh=rh, sr=sr, t_rm=tbg)
    out["ESI"] = idx.esi(ta, rh, sr)

    # History (exposure-memory) axis: leaky-integrated load above the stress
    # onset, on the continuous station series with the real per-step gap
    # (irregular sampling -> a large gap correctly discharges the integrator).
    # ITSC is integrated too, above a baseline set to the SAME exposure fraction
    # as THI > 72 (its matching upper quantile), so leakyTHI and leakyITSC are
    # comparable.
    dt_h = out["datetime"].diff().dt.total_seconds().to_numpy() / 3600.0
    dt_h[0] = np.nanmedian(dt_h[1:])
    thi = out["THI_NRC1971"].to_numpy()
    itsc = out["ITSC"].to_numpy()
    frac = float((thi > THI_STRESS_ONSET).mean())
    itsc_base = float(np.nanquantile(itsc, 1.0 - frac))
    log.info("ITSC stress baseline (matched to THI>%.0f, %.1f%% of time): %.1f",
             THI_STRESS_ONSET, 100 * frac, itsc_base)
    for tau in (6.0, 12.0, 48.0):
        out[f"leakyTHI_{int(tau):02d}h"] = history.leaky_integrate(
            thi, dt_h, tau, baseline=THI_STRESS_ONSET)
        out[f"leakyITSC_{int(tau):02d}h"] = history.leaky_integrate(
            itsc, dt_h, tau, baseline=itsc_base)
    return out


def index_columns(df: pd.DataFrame) -> list[str]:
    reserved = {"datetime", "ta", "rh", "tdp", "sr", "u", "Tw", "Tbg", "year", "day"}
    return [c for c in df.columns if c not in reserved]


# ─────────────────────────────────────────────────────────────
#  Summaries and outputs
# ─────────────────────────────────────────────────────────────

def summarise(df: pd.DataFrame) -> pd.DataFrame:
    """Per-index coverage and distribution over the whole record."""
    rows = []
    for col in index_columns(df):
        s = df[col].to_numpy()
        s = s[np.isfinite(s)]
        rows.append({
            "index": col, "n": s.size,
            "min": np.min(s), "p50": np.percentile(s, 50),
            "p90": np.percentile(s, 90), "max": np.max(s),
        })
    return pd.DataFrame(rows)


def daily_aggregate(df: pd.DataFrame) -> pd.DataFrame:
    """Daily mean and max of each index (the cadence for breakpoint work)."""
    day = df["datetime"].dt.floor("D")
    cols = index_columns(df)
    agg = df.groupby(day)[cols].agg(["mean", "max"])
    agg.columns = [f"{c}_{stat}" for c, stat in agg.columns]
    agg.index.name = "day"
    return agg.reset_index()


def plot_sample_week(df: pd.DataFrame, out_dir: Path) -> None:
    """Normalised overlay of representative indices over the hottest week."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from heataxis.constants import WONG
    from heataxis.viz import save_figure, setup_figure
    setup_figure()

    peak = df.loc[df["THI_NRC1971"].idxmax(), "datetime"]
    lo, hi = peak - pd.Timedelta(days=3), peak + pd.Timedelta(days=3)
    wk = df[(df["datetime"] >= lo) & (df["datetime"] <= hi)]
    t = (wk["datetime"] - wk["datetime"].iloc[0]).dt.total_seconds() / 86400.0

    show = {"THI_NRC1971": WONG["black"], "ETI": WONG["sky_blue"],
            "HLI": WONG["yellow"], "CCI": WONG["blue"],
            "ETIC": WONG["vermilion"], "ESI": WONG["bluish_green"]}
    fig, ax = plt.subplots(figsize=(7.09, 3.6))
    for col, colour in show.items():
        z = (wk[col] - df[col].mean()) / df[col].std()
        ax.plot(t, z, color=colour, lw=1.4, label=col)
    ax.set_xlabel("time (days)")
    ax.set_ylabel("index (z over full record)")
    ax.set_title(f"HOBO indices over the hottest week (around {peak.date()})",
                 fontsize=10, loc="left")
    ax.legend(loc="upper right", ncol=3, fontsize=7.5, frameon=False)
    save_figure(fig, "hobo_indices_sample_week", out_dir, dpi=300)
    plt.close(fig)


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                        datefmt="%H:%M:%S")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=None,
                        help="path to cow.db (default: from digimuh config)")
    parser.add_argument("--out", type=Path, default=Path("results/hobo_indices"))
    args = parser.parse_args(argv)

    db_path = args.db or load_config().database
    out_dir = args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    log.info("database: %s", db_path)

    df = load_hobo(Path(db_path))
    log.info("computing %d-row index table ...", len(df))
    full = compute_indices(df)

    n_idx = len(index_columns(full))
    full.to_csv(out_dir / "hobo_indices_5min.csv.gz", index=False,
                encoding="utf-8", compression="gzip")
    daily_aggregate(full).to_csv(out_dir / "hobo_indices_daily.csv",
                                 index=False, encoding="utf-8")
    summ = summarise(full)
    summ.to_csv(out_dir / "index_summary.csv", index=False, encoding="utf-8")
    plot_sample_week(full, out_dir)

    in_stress = float((full["THI_NRC1971"] > THI_STRESS_ONSET).mean() * 100.0)
    log.info("done: %d indices over %d rows -> %s", n_idx, len(full), out_dir)
    log.info("time above THI %.0f: %.1f%%", THI_STRESS_ONSET, in_stress)
    print("\nPer-index summary:")
    print(summ.to_string(index=False,
                         float_format=lambda v: f"{v:.2f}"))


if __name__ == "__main__":
    main()

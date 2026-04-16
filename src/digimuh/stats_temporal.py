#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════╗
# ║  DigiMuh — stats_temporal                                       ║
# ║  « circadian, cross-correlation, ETA, and crossing analyses »  ║
# ╠══════════════════════════════════════════════════════════════════╣
# ║  Temporal coupling between barn climate and rumen temperature:  ║
# ║  cross-correlation, circadian null model, derivative CCF,       ║
# ║  event-triggered averages, and climate ETA.                     ║
# ╚══════════════════════════════════════════════════════════════════╝
"""Temporal analysis functions for the broken-stick pipeline."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

log = logging.getLogger("digimuh.stats")

def compute_cross_correlation(
    rumen: pd.DataFrame, bs_results: pd.DataFrame,
    max_lag: int = 24,
) -> pd.DataFrame:
    """Cross-correlation and cross-covariance of climate vs rumen temp,
    computed separately below and above each animal's breakpoint.

    For each animal with a converged breakpoint (THI and barn temp), the
    rumen temperature time series is split at the breakpoint.  The
    normalised cross-correlation function (CCF) and raw cross-covariance
    are computed for lags -max_lag to +max_lag (in units of the 10-min
    sampling interval, so max_lag=24 = 4 hours).

    Positive lags mean the climate signal leads the rumen temperature
    response.

    Args:
        rumen: rumen_barn.csv DataFrame (must have timestamp column).
        bs_results: broken_stick_results.csv DataFrame.
        max_lag: Maximum lag in samples (default 24 = 4 hours).

    Returns:
        DataFrame with columns: animal_id, year, predictor (thi/temp),
        region (below/above), lag, xcorr, xcov.
    """
    records = []

    for prefix, env_col, conv_col, bp_col in [
        ("thi", "barn_thi", "thi_converged", "thi_breakpoint"),
        ("temp", "barn_temp", "temp_converged", "temp_breakpoint"),
    ]:
        converged = bs_results[bs_results[conv_col] == True]

        for _, row in converged.iterrows():
            aid = int(row["animal_id"])
            year = int(row["year"])
            bp = row[bp_col]

            grp = rumen[(rumen["animal_id"] == aid) & (rumen["year"] == year)].copy()
            if len(grp) < 50:
                continue

            # Sort by timestamp for temporal lag analysis
            grp = grp.sort_values("timestamp").reset_index(drop=True)
            grp["_date"] = pd.to_datetime(grp["timestamp"]).dt.date

            x_full = grp[env_col].values
            y_full = grp["body_temp"].values

            # Detrended versions: subtract per-day mean to remove
            # 24h diurnal cycle (shared periodicity → sinusoidal CCF)
            x_detrended = np.empty_like(x_full)
            y_detrended = np.empty_like(y_full)
            for date_val, day_idx in grp.groupby("_date").groups.items():
                idx = day_idx.values
                x_detrended[idx] = x_full[idx] - np.mean(x_full[idx])
                y_detrended[idx] = y_full[idx] - np.mean(y_full[idx])

            for region, mask_fn in [
                ("below", lambda x: x <= bp),
                ("above", lambda x: x > bp),
            ]:
                mask = mask_fn(x_full)
                if mask.sum() < max_lag * 3:
                    continue

                # Run both raw and detrended
                for variant, xv, yv in [
                    ("raw", x_full[mask], y_full[mask]),
                    ("detrended", x_detrended[mask], y_detrended[mask]),
                ]:
                    x = xv
                    y = yv

                    # Demean (global, on top of per-day for detrended)
                    xc = x - np.mean(x)
                    yc = y - np.mean(y)

                    sx = np.std(x, ddof=1)
                    sy = np.std(y, ddof=1)
                    n = len(x)

                    if sx < 1e-10 or sy < 1e-10:
                        continue

                    for lag in range(-max_lag, max_lag + 1):
                        if lag >= 0:
                            x_seg = xc[:n - lag] if lag > 0 else xc
                            y_seg = yc[lag:] if lag > 0 else yc
                        else:
                            x_seg = xc[-lag:]
                            y_seg = yc[:n + lag]

                        m = len(x_seg)
                        if m < 10:
                            continue

                        # Cross-covariance (raw)
                        xcov = np.sum(x_seg * y_seg) / (m - 1)
                        # Normalised cross-correlation
                        xcorr = xcov / (sx * sy)

                        records.append({
                            "animal_id": aid, "year": year,
                            "predictor": prefix,
                            "region": region,
                            "variant": variant,
                            "lag": lag,
                            "lag_minutes": lag * 10,
                            "xcorr": xcorr,
                            "xcov": xcov,
                            "n": m,
                        })

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────
#  « rumen temperature circadian null model »
# ─────────────────────────────────────────────────────────────

def compute_circadian_null_model(
    rumen: pd.DataFrame, bs_results: pd.DataFrame,
) -> pd.DataFrame:
    """Rumen temperature circadian profile on non-stress days.

    For each animal with a converged THI breakpoint, identifies days
    where barn THI stayed below the breakpoint for the entire day
    (no heat stress).  Computes the mean rumen temperature at each
    clock hour across these cool days.

    This gives the circadian null model: what rumen temperature looks
    like when the cow is entirely in the thermoneutral zone.

    Also computes the profile on stress days (THI exceeded breakpoint
    at some point during the day) for comparison.

    Args:
        rumen: rumen_barn.csv DataFrame.
        bs_results: broken_stick_results.csv DataFrame.

    Returns:
        DataFrame with columns: animal_id, year, hour, day_type
        (cool/stress), body_temp_mean, body_temp_std, n_readings.
    """
    rumen = rumen.copy()
    rumen["timestamp"] = pd.to_datetime(rumen["timestamp"])
    rumen["date"] = rumen["timestamp"].dt.date
    rumen["hour"] = rumen["timestamp"].dt.hour

    converged = bs_results[bs_results["thi_converged"] == True]
    records = []

    for _, row in converged.iterrows():
        aid = int(row["animal_id"])
        year = int(row["year"])
        bp = row["thi_breakpoint"]

        grp = rumen[(rumen["animal_id"] == aid) & (rumen["year"] == year)]
        if len(grp) < 100:
            continue

        # Classify each day: did THI exceed the breakpoint at any point?
        daily_max_thi = grp.groupby("date")["barn_thi"].max()
        cool_days = set(daily_max_thi[daily_max_thi <= bp].index)
        stress_days = set(daily_max_thi[daily_max_thi > bp].index)

        for day_type, day_set in [("cool", cool_days), ("stress", stress_days)]:
            if not day_set:
                continue
            subset = grp[grp["date"].isin(day_set)]
            for hour, hdata in subset.groupby("hour"):
                if len(hdata) < 5:
                    continue
                records.append({
                    "animal_id": aid, "year": year,
                    "hour": int(hour),
                    "day_type": day_type,
                    "body_temp_mean": hdata["body_temp"].mean(),
                    "body_temp_std": hdata["body_temp"].std(),
                    "thi_mean": hdata["barn_thi"].mean(),
                    "barn_temp_mean": hdata["barn_temp"].mean(),
                    "n_readings": len(hdata),
                    "n_days": hdata["date"].nunique(),
                })

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────
#  « THI daily exceedance profile »
# ─────────────────────────────────────────────────────────────

def compute_thi_daily_profile(
    rumen: pd.DataFrame, bs_results: pd.DataFrame,
) -> pd.DataFrame:
    """Barn THI profile across 24h, by month.

    For each month of the observation period, computes the mean barn
    THI at each clock hour.  Overlaid with the herd median breakpoint,
    this shows when heat stress typically begins and ends each month.

    Args:
        rumen: rumen_barn.csv DataFrame.
        bs_results: broken_stick_results.csv DataFrame.

    Returns:
        DataFrame with columns: year, month, month_label, hour,
        thi_mean, thi_std, thi_q25, thi_q75, n_readings.
        Also includes herd_median_bp.
    """
    rumen = rumen.copy()
    rumen["timestamp"] = pd.to_datetime(rumen["timestamp"])
    rumen["hour"] = rumen["timestamp"].dt.hour
    rumen["month"] = rumen["timestamp"].dt.month
    rumen["year"] = rumen["year"].astype(int)

    month_names = {6: "June", 7: "July", 8: "August", 9: "September"}

    # Herd median THI breakpoint
    thi_conv = bs_results[bs_results["thi_converged"] == True]
    herd_median_bp = thi_conv["thi_breakpoint"].median() if not thi_conv.empty else np.nan

    records = []
    for (year, month), grp in rumen.groupby(["year", "month"]):
        if int(month) not in month_names:
            continue
        for hour, hdata in grp.groupby("hour"):
            records.append({
                "year": int(year),
                "month": int(month),
                "month_label": f"{month_names[int(month)]} {int(year)}",
                "hour": int(hour),
                "thi_mean": hdata["barn_thi"].mean(),
                "thi_std": hdata["barn_thi"].std(),
                "thi_q25": hdata["barn_thi"].quantile(0.25),
                "thi_q75": hdata["barn_thi"].quantile(0.75),
                "n_readings": len(hdata),
                "herd_median_bp": herd_median_bp,
            })

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────
#  « derivative cross-correlation »
# ─────────────────────────────────────────────────────────────

def compute_derivative_ccf(
    rumen: pd.DataFrame, bs_results: pd.DataFrame,
    max_lag: int = 24,
) -> pd.DataFrame:
    """Cross-correlation of rate-of-change signals: dTHI/dt vs dTbody/dt.

    Instead of correlating the raw levels (contaminated by shared
    diurnal cycle), we correlate the temporal derivatives.  This asks:
    "when the barn heats up, how long until the cow heats up?"

    The derivative removes DC offsets and slow trends while preserving
    the temporal coupling of changes.  The rumen's thermal inertia
    (~100 L) means the derivative response is delayed by 30-90 min.

    Args:
        rumen: rumen_barn.csv DataFrame.
        bs_results: broken_stick_results.csv DataFrame.
        max_lag: Maximum lag in 10-min samples (default 24 = 4 hours).

    Returns:
        DataFrame with columns: animal_id, year, predictor, region,
        lag, lag_minutes, dxcorr (normalised CCF of derivatives).
    """
    records = []

    for prefix, env_col, conv_col, bp_col in [
        ("thi", "barn_thi", "thi_converged", "thi_breakpoint"),
        ("temp", "barn_temp", "temp_converged", "temp_breakpoint"),
    ]:
        converged = bs_results[bs_results[conv_col] == True]

        for _, row in converged.iterrows():
            aid = int(row["animal_id"])
            year = int(row["year"])
            bp = row[bp_col]

            grp = rumen[(rumen["animal_id"] == aid) & (rumen["year"] == year)].copy()
            if len(grp) < 100:
                continue

            grp = grp.sort_values("timestamp").reset_index(drop=True)

            # Compute derivatives (forward difference)
            x_full = grp[env_col].values
            y_full = grp["body_temp"].values
            dx = np.diff(x_full)
            dy = np.diff(y_full)

            # Use midpoint climate values for region assignment
            x_mid = (x_full[:-1] + x_full[1:]) / 2.0

            for region, mask_fn in [
                ("below", lambda x: x <= bp),
                ("above", lambda x: x > bp),
            ]:
                mask = mask_fn(x_mid)
                if mask.sum() < max_lag * 3:
                    continue

                dxr = dx[mask]
                dyr = dy[mask]

                # Demean
                dxc = dxr - np.mean(dxr)
                dyc = dyr - np.mean(dyr)

                sx = np.std(dxr, ddof=1)
                sy = np.std(dyr, ddof=1)
                n = len(dxr)

                if sx < 1e-12 or sy < 1e-12:
                    continue

                for lag in range(-max_lag, max_lag + 1):
                    if lag >= 0:
                        x_seg = dxc[:n - lag] if lag > 0 else dxc
                        y_seg = dyc[lag:] if lag > 0 else dyc
                    else:
                        x_seg = dxc[-lag:]
                        y_seg = dyc[:n + lag]

                    m = len(x_seg)
                    if m < 10:
                        continue

                    xcov = np.sum(x_seg * y_seg) / (m - 1)
                    dxcorr = xcov / (sx * sy)

                    records.append({
                        "animal_id": aid, "year": year,
                        "predictor": prefix,
                        "region": region,
                        "lag": lag,
                        "lag_minutes": lag * 10,
                        "dxcorr": dxcorr,
                        "n": m,
                    })

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────
#  « event-triggered average at breakpoint crossing »
# ─────────────────────────────────────────────────────────────

def compute_event_triggered_average(
    rumen: pd.DataFrame, bs_results: pd.DataFrame,
    window: int = 36,
    min_gap: int = 6,
    crossing_hour_range: tuple[int, int] | None = None,
) -> pd.DataFrame:
    """Peri-event average of rumen temperature around THI breakpoint crossings.

    Finds moments when barn THI crosses the animal's breakpoint upward
    (heat stress onset).  Extracts a window of rumen temperature
    centred on each crossing event and averages across events.

    Args:
        rumen: rumen_barn.csv DataFrame.
        bs_results: broken_stick_results.csv DataFrame.
        window: Half-window in 10-min samples (default 36 = 6 hours).
        min_gap: Minimum samples between events to avoid overlap
            (default 6 = 1 hour).
        crossing_hour_range: If set, only include crossing events whose
            clock hour falls within [start, end).  E.g. (8, 11) keeps
            crossings at 8:00-10:59.  None = all hours.

    Returns:
        DataFrame with columns: animal_id, year, predictor, event_id,
        relative_lag (samples from crossing), relative_minutes,
        body_temp, climate_val.  Also a summary DataFrame.
    """
    traces = []
    summaries = []

    for prefix, env_col, conv_col, bp_col in [
        ("thi", "barn_thi", "thi_converged", "thi_breakpoint"),
        ("temp", "barn_temp", "temp_converged", "temp_breakpoint"),
    ]:
        converged = bs_results[bs_results[conv_col] == True]

        for _, row in converged.iterrows():
            aid = int(row["animal_id"])
            year = int(row["year"])
            bp = row[bp_col]

            grp = rumen[(rumen["animal_id"] == aid) & (rumen["year"] == year)].copy()
            if len(grp) < window * 3:
                continue

            grp = grp.sort_values("timestamp").reset_index(drop=True)
            x = grp[env_col].values
            y = grp["body_temp"].values
            ts = pd.to_datetime(grp["timestamp"]).values
            n = len(x)
            # Animal's overall mean body temp (for centering across animals)
            animal_mean_bt = np.mean(y)

            # Find upward crossings: x[i] <= bp and x[i+1] > bp
            # Only count if consecutive readings are < 30 min apart
            # (rejects milking gap artifacts where 3h of data is missing)
            below = x[:-1] <= bp
            above = x[1:] > bp
            dt = np.diff(ts).astype("timedelta64[m]").astype(float)
            consecutive = dt < 30
            crossings = np.where(below & above & consecutive)[0]

            if len(crossings) == 0:
                continue

            # Filter by clock hour if requested
            if crossing_hour_range is not None:
                h_start, h_end = crossing_hour_range
                hours = pd.DatetimeIndex(ts[crossings]).hour
                crossings = crossings[(hours >= h_start) & (hours < h_end)]
                if len(crossings) == 0:
                    continue

            # Filter: enforce minimum gap between events
            filtered = [crossings[0]]
            for c in crossings[1:]:
                if c - filtered[-1] >= min_gap:
                    filtered.append(c)
            crossings = np.array(filtered)

            # Extract windows around each crossing
            event_count = 0
            for c in crossings:
                start = c - window
                end = c + window + 1
                if start < 0 or end > n:
                    continue

                event_count += 1
                # Baseline-subtract body temp (mean of pre-event period)
                pre_mean = np.mean(y[start:c])
                for j, idx in enumerate(range(start, end)):
                    t_point = pd.Timestamp(ts[idx])
                    traces.append({
                        "animal_id": aid, "year": year,
                        "predictor": prefix,
                        "event_id": event_count,
                        "relative_lag": j - window,
                        "relative_minutes": (j - window) * 10,
                        "body_temp": y[idx],
                        "body_temp_centered": y[idx] - animal_mean_bt,
                        "body_temp_baseline": y[idx] - pre_mean,
                        "climate_val": x[idx],
                        "clock_hour": t_point.hour,
                    })

            if event_count > 0:
                summaries.append({
                    "animal_id": aid, "year": year,
                    "predictor": prefix,
                    "n_events": event_count,
                })

    traces_df = pd.DataFrame(traces)
    summaries_df = pd.DataFrame(summaries)
    return traces_df, summaries_df


# ─────────────────────────────────────────────────────────────
#  « breakpoint crossing clock-time raster »
# ─────────────────────────────────────────────────────────────

def compute_crossing_times(
    rumen: pd.DataFrame, bs_results: pd.DataFrame,
    min_gap: int = 6,
) -> pd.DataFrame:
    """Extract clock times of all breakpoint crossing events.

    For each animal, finds moments when barn THI (or barn temp) crosses
    the animal's individual breakpoint upward, and records the clock
    time.  Used for the activation raster plot.

    Args:
        rumen: rumen_barn.csv DataFrame.
        bs_results: broken_stick_results.csv DataFrame.
        min_gap: Minimum samples between events (default 6 = 1 hour).

    Returns:
        DataFrame with columns: animal_id, year, predictor, breakpoint,
        crossing_timestamp, clock_hour, clock_minute, day_fraction.
    """
    rumen = rumen.copy()
    rumen["timestamp"] = pd.to_datetime(rumen["timestamp"])

    records = []

    for prefix, env_col, conv_col, bp_col in [
        ("thi", "barn_thi", "thi_converged", "thi_breakpoint"),
        ("temp", "barn_temp", "temp_converged", "temp_breakpoint"),
    ]:
        converged = bs_results[bs_results[conv_col] == True]

        for _, row in converged.iterrows():
            aid = int(row["animal_id"])
            year = int(row["year"])
            bp = row[bp_col]

            grp = rumen[(rumen["animal_id"] == aid) & (rumen["year"] == year)].copy()
            if len(grp) < 50:
                continue

            grp = grp.sort_values("timestamp").reset_index(drop=True)
            x = grp[env_col].values
            ts = pd.to_datetime(grp["timestamp"]).values
            n = len(x)

            # Upward crossings — only if consecutive readings < 30 min apart
            below = x[:-1] <= bp
            above = x[1:] > bp
            dt = np.diff(ts).astype("timedelta64[m]").astype(float)
            consecutive = dt < 30
            crossings = np.where(below & above & consecutive)[0]

            if len(crossings) == 0:
                continue

            # Filter minimum gap
            filtered = [crossings[0]]
            for c in crossings[1:]:
                if c - filtered[-1] >= min_gap:
                    filtered.append(c)

            for c in filtered:
                t = pd.Timestamp(ts[c])
                records.append({
                    "animal_id": aid,
                    "year": year,
                    "predictor": prefix,
                    "breakpoint": bp,
                    "crossing_timestamp": t,
                    "date": t.date(),
                    "clock_hour": t.hour,
                    "clock_minute": t.minute,
                    "day_fraction": t.hour + t.minute / 60.0,
                })

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────
#  « climate ETA: THI + barn temp around breakpoint crossings »
# ─────────────────────────────────────────────────────────────

def compute_climate_eta(
    rumen: pd.DataFrame, bs_results: pd.DataFrame,
    window: int = 36,
    min_gap: int = 6,
    crossing_hour_range: tuple[int, int] | None = (8, 11),
) -> pd.DataFrame:
    """Climate signal around breakpoint crossings, normalised to the breakpoint.

    For each animal and each crossing event (upward), extracts both
    barn THI and barn temperature in a ±window around the crossing.
    The trigger predictor is normalised by subtracting the animal's
    breakpoint value, so y=0 corresponds to the threshold.

    Two modes:
    - 'thi' crossings: THI normalised, barn temp as companion
    - 'temp' crossings: barn temp normalised, THI as companion

    Args:
        rumen: rumen_barn.csv DataFrame.
        bs_results: broken_stick_results.csv DataFrame.
        window: Half-window in 10-min samples (default 36 = 6 hours).
        min_gap: Minimum samples between events (default 6 = 1 hour).
        crossing_hour_range: Restrict to crossings in this clock-hour
            range.  Default (8, 11) = 8:00–10:59.

    Returns:
        DataFrame with columns: animal_id, year, trigger (thi/temp),
        event_id, relative_minutes, thi_raw, temp_raw, thi_norm,
        temp_norm, breakpoint_thi, breakpoint_temp.
    """
    rumen = rumen.copy()
    rumen["timestamp"] = pd.to_datetime(rumen["timestamp"])
    records = []

    for trigger, env_col, conv_col, bp_col in [
        ("thi", "barn_thi", "thi_converged", "thi_breakpoint"),
        ("temp", "barn_temp", "temp_converged", "temp_breakpoint"),
    ]:
        converged = bs_results[bs_results[conv_col] == True]

        for _, row in converged.iterrows():
            aid = int(row["animal_id"])
            year = int(row["year"])
            bp = row[bp_col]

            # Also get the other predictor's breakpoint if available
            bp_thi = row.get("thi_breakpoint", np.nan)
            bp_temp = row.get("temp_breakpoint", np.nan)

            grp = rumen[(rumen["animal_id"] == aid) & (rumen["year"] == year)].copy()
            if len(grp) < window * 3:
                continue

            grp = grp.sort_values("timestamp").reset_index(drop=True)
            x = grp[env_col].values
            thi_vals = grp["barn_thi"].values
            temp_vals = grp["barn_temp"].values
            ts = grp["timestamp"].values
            n = len(x)

            # Upward crossings with gap filter
            below = x[:-1] <= bp
            above = x[1:] > bp
            dt = np.diff(ts).astype("timedelta64[m]").astype(float)
            consecutive = dt < 30
            crossings = np.where(below & above & consecutive)[0]

            if len(crossings) == 0:
                continue

            # Filter by clock hour
            if crossing_hour_range is not None:
                h_start, h_end = crossing_hour_range
                hours = pd.DatetimeIndex(ts[crossings]).hour
                crossings = crossings[(hours >= h_start) & (hours < h_end)]
                if len(crossings) == 0:
                    continue

            # Enforce minimum gap
            filtered = [crossings[0]]
            for c in crossings[1:]:
                if c - filtered[-1] >= min_gap:
                    filtered.append(c)
            crossings = np.array(filtered)

            event_count = 0
            for c in crossings:
                start = c - window
                end = c + window + 1
                if start < 0 or end > n:
                    continue

                event_count += 1
                for j, idx in enumerate(range(start, end)):
                    records.append({
                        "animal_id": aid,
                        "year": year,
                        "trigger": trigger,
                        "event_id": event_count,
                        "relative_minutes": (j - window) * 10,
                        "thi_raw": thi_vals[idx],
                        "temp_raw": temp_vals[idx],
                        "thi_norm": thi_vals[idx] - bp_thi,
                        "temp_norm": temp_vals[idx] - bp_temp,
                        "breakpoint_thi": bp_thi,
                        "breakpoint_temp": bp_temp,
                    })

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────
#  « breakpoint stability (ICC) »
# ─────────────────────────────────────────────────────────────

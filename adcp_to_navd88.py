#!/usr/bin/env python3
"""
ADCP Water Level On NAVD88
=============================
Converts the Signature 1000 ADCP water level (sig1000_waves_ALL.csv,
deployed off Marconi Dec 2024 - Mar 2025) into a NAVD88 water-level file
that extract_elevation_contours.py can read.

WHY A CONVERSION. The ADCP's `water_level` is water depth minus a
constant (depth - water_level = 21.02 m on every row), i.e. relative to
the MEAN level over the deployment, not to any survey datum. Beach
elevations must be NAVD88 to match the camera geometry and every other
product.

HOW. Mean water level over the deployment is taken from the NOAA
Chatham gauge (8447435, reported directly in NAVD88) over exactly the
same hours:

    level_navd88 = water_level + (mean Chatham NAVD88 - mean ADCP water_level)

The ADCP and the gauge share the same winter weather, so this carries the
seasonal sea-level anomaly that a fixed MSL-to-NAVD88 constant would
miss. The constant used elsewhere (+0.09 m, tide model to NAVD88, from
the 2026 GNSS-R/OPUS comparison) is printed alongside as a cross-check;
the two should agree to within about 0.1 m.

The same comparison also fits ADCP ~ a * gauge(t - lag) + b. It is a
CLOCK CHECK: Marconi's tide leads Chatham's by roughly 40-50 min (from
the GNSS-R fit). A lag of hours instead means the ADCP timestamps are
not UTC, and every water level would be assigned to the wrong image.

OUTPUT columns: time (UTC, ISO), water_level_navd88, water_level_adcp.

Usage:
    python3 adcp_to_navd88.py --adcp sig1000_waves_ALL.csv \\
        --output /mnt/I2Rgus_Data/Chelsea_calibration/adcp_water_level_navd88.csv
    (offline, or to skip the gauge:  --offset 0.09)
"""

import sys
import argparse
from datetime import timedelta

import numpy as np
import pandas as pd

MSL_TO_NAVD88 = 0.09      # tide model / mean sea level -> NAVD88, see process_historical.py
EXPECTED_LAG_MIN = (-90, 0)   # Marconi leads Chatham; GNSS-R fits gave -42 to -48 min


def to_epoch(times):
    """Seconds since 1970, independent of the pandas time resolution."""
    t = pd.to_datetime(pd.Series(times), utc=True)
    return ((t - pd.Timestamp("1970-01-01", tz="UTC")) // pd.Timedelta(seconds=1)).to_numpy(float)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--adcp", default="sig1000_waves_ALL.csv")
    ap.add_argument("--output", required=True)
    ap.add_argument("--station", default="8447435", help="NOAA gauge (default Chatham)")
    ap.add_argument("--offset", type=float, default=None,
                    help="Use this ADCP -> NAVD88 offset (m) instead of deriving it from the gauge.")
    args = ap.parse_args()

    adcp = pd.read_csv(args.adcp)
    adcp["time"] = pd.to_datetime(adcp["time"], utc=True)      # naive timestamps are UTC
    adcp = adcp.dropna(subset=["water_level"]).sort_values("time")
    a_ep = to_epoch(adcp["time"])
    a_wl = adcp["water_level"].to_numpy(float)
    print(f"ADCP: {len(adcp)} hourly readings, {adcp.time.min():%Y-%m-%d %H:%M} to "
          f"{adcp.time.max():%Y-%m-%d %H:%M} UTC; depth - water_level = "
          f"{(adcp.water_depth - adcp.water_level).mean():.3f} m (constant: level is relative "
          f"to the deployment mean)")

    offset = args.offset
    if offset is None:
        from compare_gnssr_to_gauge import fetch_gauge
        start = adcp.time.min().to_pydatetime() - timedelta(days=1)
        end = adcp.time.max().to_pydatetime() + timedelta(days=1)
        print(f"Downloading gauge {args.station} (NAVD88) for the same period ...")
        g = fetch_gauge(args.station, start, end)
        g_ep, g_lv = to_epoch(g["time"]), g["level"].to_numpy(float)

        def gauge_at(ep):
            out = np.interp(ep, g_ep, g_lv, left=np.nan, right=np.nan)
            i = np.clip(np.searchsorted(g_ep, ep), 1, len(g_ep) - 1)
            out[(g_ep[i] - g_ep[i - 1]) > 720] = np.nan
            return out

        # Clock check: best lag of ADCP against the gauge.
        best = None
        for lag_min in range(-360, 361, 6):
            x = gauge_at(a_ep - lag_min * 60)
            m = np.isfinite(x)
            if m.sum() < 200:
                continue
            A = np.c_[x[m], np.ones(m.sum())]
            coef, *_ = np.linalg.lstsq(A, a_wl[m], rcond=None)
            rms = float(np.sqrt(np.mean((a_wl[m] - A @ coef) ** 2)))
            if best is None or rms < best[3]:
                best = (lag_min, coef[0], coef[1], rms, int(m.sum()))
        if best is None:
            sys.exit("Too little overlap with the gauge; rerun with --offset 0.09")
        lag_min, a, b, rms, n = best
        print(f"Fit ADCP = {a:.2f} x gauge(t {lag_min:+d} min) {b:+.3f}: RMS {rms:.3f} m over {n} hours")
        if not (EXPECTED_LAG_MIN[0] <= lag_min <= EXPECTED_LAG_MIN[1]):
            print(f"  WARNING: lag {lag_min:+d} min is outside the expected "
                  f"{EXPECTED_LAG_MIN[0]}..{EXPECTED_LAG_MIN[1]} min. The ADCP clock may not be "
                  f"UTC; check before using these water levels.")
        else:
            print("  clock check OK (Marconi leads Chatham by the expected amount)")

        # Datum: same hours, mean level at the gauge (NAVD88) vs mean ADCP level.
        x = gauge_at(a_ep)
        m = np.isfinite(x)
        offset = float(np.mean(x[m]) - np.mean(a_wl[m]))
        print(f"Mean level over the deployment: Chatham {np.mean(x[m]):+.3f} m NAVD88, "
              f"ADCP {np.mean(a_wl[m]):+.3f} m (own zero)")
        print(f"ADCP -> NAVD88 offset: {offset:+.3f} m   (fixed MSL constant would give "
              f"{MSL_TO_NAVD88 - np.mean(a_wl):+.3f} m; difference "
              f"{offset - (MSL_TO_NAVD88 - np.mean(a_wl)):+.3f} m)")
    else:
        print(f"ADCP -> NAVD88 offset: {offset:+.3f} m (given)")

    out = pd.DataFrame({
        "time": adcp["time"].dt.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "water_level_navd88": np.round(a_wl + offset, 4),
        "water_level_adcp": np.round(a_wl, 4),
    })
    out.to_csv(args.output, index=False)
    print(f"Wrote {len(out)} rows to {args.output} "
          f"(range {out.water_level_navd88.min():+.2f} to {out.water_level_navd88.max():+.2f} m NAVD88)")


if __name__ == "__main__":
    main()

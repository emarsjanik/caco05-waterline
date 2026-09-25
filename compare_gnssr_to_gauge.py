#!/usr/bin/env python3
"""
Check The GNSS-R Vertical Datum Against A NOAA Tide Gauge
============================================================
Compares the station's GNSS-R water level (gnssrefl spline, treated as
NAVD88 throughout this project) with NOAA CO-OPS water levels requested
directly in NAVD88, and estimates the constant offset between them.

WHY. Every elevation in the DEM is a GNSS-R water level, so its datum
is the DEM's datum. Agreement with a tide model's SHAPE and TIMING
confirms the tide; only a NAVD88-referenced record confirms the constant
LEVEL. This is an independent, repeatable check of that level -- worth
re-running after any change to the GNSS station (antenna, height,
gnssrefl settings, geoid model) and before publishing elevations.

GAUGE. 8447435 Chatham, Lydia Cove, MA -- inside Chatham Harbor, ~30 km
south of Marconi. Its tide is not the open-coast tide at Marconi: the
amplitude and timing differ, and storm setup in the harbour differs.
So the script does not trust a point-by-point difference. It reports:

  1. TIDALLY AVERAGED OFFSET (the answer). Mean GNSS-R minus mean gauge
     over each complete day, with the gauge sampled at the GNSS-R times.
     Daily means average the tide out, so amplitude and timing
     differences largely cancel; what is left is datum offset plus a
     few cm of genuine mean-level difference between the two sites.
  2. REGRESSION. GNSS-R = a * gauge(t - lag) + b, with the lag that fits
     best. a and lag describe the tide difference between the sites;
     they are reported so an implausible fit is visible.
  3. MEAN-SEA-LEVEL SANITY CHECK. The GNSS-R mean over the whole record
     against NOAA's published MSL for the gauge (+0.018 m NAVD88,
     2015-2017 series; sea level has risen a few cm since).

READING THE RESULT.
  |offset| < 0.05 m   consistent: GNSS-R is on NAVD88 to within what two
                      sites 30 km apart can show.
  0.05 - 0.15 m       small offset: worth checking how the GNSS-R
                      orthometric height was computed.
  > 0.15 m            datum problem likely (e.g. EGM96 instead of
                      GEOID18, or a wrong antenna height). Every DEM
                      elevation is off by about this amount.

Usage:
    python3 compare_gnssr_to_gauge.py /path/to/usgs_spline_out.txt
    python3 compare_gnssr_to_gauge.py SPLINE --start 2026-08-01 --end 2026-09-23 \\
        --plot gnssr_vs_chatham.png
    python3 compare_gnssr_to_gauge.py SPLINE --gauge-csv chatham.csv   (offline)
"""

import sys
import argparse
import urllib.request
from io import StringIO
from datetime import timedelta

import numpy as np
import pandas as pd

from extract_elevation_contours import load_gnssr_spline

API = ("https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=water_level"
       "&station={station}&begin_date={begin}&end_date={end}&datum=NAVD&units=metric"
       "&time_zone=gmt&format=csv&application=caco05_waterline")
GAUGE_MSL_NAVD88 = 0.778 - 0.760   # MSL and NAVD88 above MLLW, NOAA datums page, 8447435
MIN_DAY_COVERAGE = 0.9             # a day counts only if GNSS-R covers 90% of it


def to_epoch(times):
    """Seconds since 1970 for a tz-aware datetime Series, independent of
    the pandas time resolution (ns/us/s differ between pandas versions)."""
    return ((times - pd.Timestamp("1970-01-01", tz="UTC")) // pd.Timedelta(seconds=1)).to_numpy()


def fetch_gauge(station, start, end):
    """6-minute water level in m NAVD88, UTC, in <=31-day requests (the API limit)."""
    frames = []
    t = start
    while t <= end:
        stop = min(t + timedelta(days=30), end)
        url = API.format(station=station, begin=t.strftime("%Y%m%d"), end=stop.strftime("%Y%m%d"))
        text = urllib.request.urlopen(url, timeout=60).read().decode()
        if text.lstrip().lower().startswith("error") or "Date Time" not in text:
            sys.exit(f"NOAA API returned no data for {t:%Y-%m-%d}..{stop:%Y-%m-%d}:\n{text[:300]}")
        frames.append(pd.read_csv(StringIO(text)))
        t = stop + timedelta(days=1)
    return tidy_gauge(pd.concat(frames, ignore_index=True))


def tidy_gauge(df):
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns={"Date Time": "time", "Water Level": "level"})
    df["time"] = pd.to_datetime(df["time"], utc=True)
    df["level"] = pd.to_numeric(df["level"], errors="coerce")
    df = df.dropna(subset=["level"]).drop_duplicates("time").sort_values("time")
    return df[["time", "level"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("gnssr", help="gnssrefl subdaily spline file (usgs_spline_out.txt)")
    ap.add_argument("--station", default="8447435", help="NOAA CO-OPS station (default Chatham)")
    ap.add_argument("--start", help="YYYY-MM-DD (default: start of GNSS-R record)")
    ap.add_argument("--end", help="YYYY-MM-DD (default: end of GNSS-R record)")
    ap.add_argument("--gauge-csv", help="Use a saved CO-OPS CSV (datum NAVD, metric, GMT) "
                                        "instead of downloading.")
    ap.add_argument("--save-gauge", help="Save the downloaded gauge data to this CSV.")
    ap.add_argument("--max-lag-minutes", type=int, default=180)
    ap.add_argument("--plot", help="Write a comparison figure (PNG).")
    args = ap.parse_args()

    g_ep, g_lv, _, _ = load_gnssr_spline(args.gnssr)
    gn = pd.DataFrame({"time": pd.to_datetime(g_ep, unit="s", utc=True), "level": g_lv})
    if args.start:
        gn = gn[gn["time"] >= pd.Timestamp(args.start, tz="UTC")]
    if args.end:
        gn = gn[gn["time"] < pd.Timestamp(args.end, tz="UTC") + pd.Timedelta(days=1)]
    if len(gn) < 50:
        sys.exit("Too few GNSS-R readings in the requested window.")
    start = gn["time"].min().to_pydatetime()
    end = gn["time"].max().to_pydatetime()

    if args.gauge_csv:
        gauge = tidy_gauge(pd.read_csv(args.gauge_csv))
    else:
        print(f"Downloading NOAA {args.station} water level (NAVD88) "
              f"{start:%Y-%m-%d} to {end:%Y-%m-%d} ...")
        gauge = fetch_gauge(args.station, start, end)
        if args.save_gauge:
            gauge.assign(time=gauge["time"].dt.strftime("%Y-%m-%d %H:%M")) \
                 .rename(columns={"time": "Date Time", "level": "Water Level"}) \
                 .to_csv(args.save_gauge, index=False)

    gs = to_epoch(gauge["time"])
    gl = gauge["level"].to_numpy(float)
    ts = to_epoch(gn["time"])
    lv = gn["level"].to_numpy(float)

    def gauge_at(epochs):
        """Gauge level at `epochs`, NaN where the gauge has a gap > 12 min."""
        out = np.interp(epochs, gs, gl, left=np.nan, right=np.nan)
        i = np.clip(np.searchsorted(gs, epochs), 1, len(gs) - 1)
        out[(gs[i] - gs[i - 1]) > 720] = np.nan
        return out

    at = gauge_at(ts)
    ok = np.isfinite(at)
    print()
    print("=" * 72)
    print(f"GNSS-R vs NOAA {args.station} (both m NAVD88)")
    print("=" * 72)
    print(f"GNSS-R readings compared : {int(ok.sum())} of {len(ts)}, "
          f"{start:%Y-%m-%d} to {end:%Y-%m-%d}")
    if ok.sum() < 50:
        sys.exit("Too little overlap with the gauge record.")

    # 1. Tidally averaged offset: complete days only.
    df = pd.DataFrame({"time": gn["time"].to_numpy()[ok], "gnssr": lv[ok], "gauge": at[ok]})
    df["day"] = df["time"].dt.strftime("%Y-%m-%d")
    per_day = df.groupby("day").agg(n=("gnssr", "size"), gnssr=("gnssr", "mean"),
                                    gauge=("gauge", "mean"),
                                    span=("time", lambda t: (t.max() - t.min()).total_seconds()))
    full = per_day[per_day["span"] >= MIN_DAY_COVERAGE * 86400 * (1 - 1 / max(per_day["n"].median(), 2))]
    daily = (full["gnssr"] - full["gauge"]).to_numpy()

    # 2. Regression with lag.
    best = None
    for lag_min in range(-args.max_lag_minutes, args.max_lag_minutes + 1, 6):
        x = gauge_at(ts - lag_min * 60)
        m = np.isfinite(x)
        if m.sum() < 50:
            continue
        A = np.c_[x[m], np.ones(m.sum())]
        coef, *_ = np.linalg.lstsq(A, lv[m], rcond=None)
        resid = lv[m] - A @ coef
        r2 = 1 - np.sum(resid ** 2) / np.sum((lv[m] - lv[m].mean()) ** 2)
        if best is None or r2 > best[0]:
            best = (r2, lag_min, coef[0], coef[1], float(np.sqrt(np.mean(resid ** 2))))

    raw = lv[ok] - at[ok]
    print()
    print("Point-by-point GNSS-R minus gauge (NOT the answer -- different tides):")
    print(f"   mean {raw.mean():+.3f} m, std {raw.std():.3f} m")
    print()
    if best:
        r2, lag, a, b, rms = best
        print("Regression  GNSS-R = a * gauge(t - lag) + b:")
        print(f"   lag {lag:+d} min, a = {a:.3f}, b = {b:+.3f} m, R2 {r2:.3f}, residual RMS {rms:.3f} m")
        print("   (a and lag describe how Marconi's tide differs from Chatham Harbor's;")
        print("    b mixes datum offset with a*mean-level, so use the daily figure below)")
    print()
    print(f"MEAN-SEA-LEVEL CHECK: GNSS-R mean over the record {lv.mean():+.3f} m NAVD88;")
    print(f"   NOAA MSL at the gauge {GAUGE_MSL_NAVD88:+.3f} m NAVD88 (2015-2017), a few cm")
    print(f"   higher today. Gauge mean over the same times {at[ok].mean():+.3f} m.")
    print()
    if len(daily) == 0:
        print("No complete days in common -- cannot give a tidally averaged offset.")
        return
    offset = float(np.median(daily))
    sem = float(np.std(daily, ddof=1) / np.sqrt(len(daily))) if len(daily) > 1 else float("nan")
    print(f"TIDALLY AVERAGED OFFSET (GNSS-R minus gauge, daily means, {len(daily)} complete day(s)):")
    print(f"   median {offset:+.3f} m   (spread of days {np.std(daily):.3f} m, "
          f"standard error ~{sem:.3f} m)")
    print()
    if abs(offset) < 0.05:
        verdict = ("CONSISTENT: GNSS-R agrees with NAVD88 at the gauge to within what two "
                   "sites 30 km apart can show. No datum correction indicated.")
    elif abs(offset) < 0.15:
        verdict = ("SMALL OFFSET: worth checking how the GNSS-R orthometric height was "
                   "computed (geoid model, antenna height) before correcting anything.")
    else:
        verdict = (f"DATUM PROBLEM LIKELY: GNSS-R reads {abs(offset):.2f} m "
                   f"{'LOW' if offset < 0 else 'HIGH'} against NAVD88. Every DEM elevation "
                   "carries about this error. Check the geoid model (EGM96 vs GEOID18) and "
                   "antenna height in the gnssrefl station settings.")
    print("VERDICT: " + verdict)

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7))
        ax1.plot(gauge["time"], gauge["level"], lw=0.6, color="0.5", label=f"NOAA {args.station}")
        ax1.plot(gn["time"], gn["level"], ".", ms=2, color="C0", label="GNSS-R")
        ax1.set_ylabel("m NAVD88"); ax1.legend(loc="upper right")
        ax1.set_xlim(gn["time"].min(), gn["time"].max())
        ax2.bar(pd.to_datetime(full.index), daily, color="C1", width=0.8)
        ax2.axhline(offset, color="k", lw=1, label=f"median {offset:+.3f} m")
        ax2.axhline(0, color="0.6", lw=0.8)
        ax2.set_ylabel("daily mean GNSS-R - gauge (m)"); ax2.legend(loc="upper right")
        fig.tight_layout(); fig.savefig(args.plot, dpi=110)
        print(f"Figure: {args.plot}")


if __name__ == "__main__":
    main()

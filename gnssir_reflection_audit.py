#!/usr/bin/env python3
"""
GNSS-IR Reflection Audit: Where Did Each Arc Reflect, And Was It The Sea?
==========================================================================
At Marconi the antenna is ~18.7 m above NAVD88 and the waterline only
~55-90 m seaward of it (more at low tide). A satellite arc at elevation e
reflects about RH / tan(e) metres out along its azimuth:

    5 deg -> ~210 m (sea)    10 deg -> ~105 m (sea)    15 deg -> ~70 m (beach at low tide)

and arcs near the edges of the 353-173 deg azimuth window run ALONG the
shore, not out to sea. A reflection from wet sand does not fail -- it
returns a plausible "water level" equal to the beach elevation. In surf the
sea's coherent reflection fades first at higher elevations (by ~0.4 m wave
height at 10 deg), leaving the smooth beach: storm readings 3-4 m above the
tide may be the upper beach, not wave setup.

elevation_quality.py compares arcs with the spline consensus; when most
arcs are contaminated, the consensus is too. This compares every arc with
an INDEPENDENT reference -- the Marconi still-water level from the Chatham
gauge (owg_live.py's fit) -- and puts each arc's reflection point on the
map at that moment's tide:

  1. residual (arc minus reference) by elevation band, azimuth band, and
     reflection point at sea vs on the beach, in calm and rough seas;
  2. candidate masks (upper elevation x azimuth window): how many arcs each
     keeps and how well they match the reference, calm and rough;
  3. the recommended gnssrefl settings for station.json.

Nothing is changed. Revert any mask you adopt by restoring station.json.

Usage (on the station):
    python3 gnssir_reflection_audit.py \\
        --subdaily /home/argus_user/GNSS/v4.1/products/refl_code/Files/usgs/usgs_2026_subdaily_edit.txt \\
        --station-json /home/argus_user/GNSS/v4.1/station/resources/station.json
"""

import sys
import json
import math
import argparse
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
# beach geometry from the Jan 2025 lidar (sea_patch.py): the 0 m NAVD88 contour
SHORE_ORIGIN = (420150.0, 4638400.0)
SHORE_BEARING = 351.0          # deg; seaward normal is 81 deg
FORESHORE_SLOPE = 0.08


def latlon_to_utm19(lat, lon):
    """WGS84 -> UTM zone 19N (m), standard series (mm-level, no pyproj needed)."""
    a, f, k0 = 6378137.0, 1 / 298.257223563, 0.9996
    e2 = f * (2 - f); ep2 = e2 / (1 - e2)
    lat, lon = math.radians(lat), math.radians(lon)
    lon0 = math.radians(-69.0)
    N = a / math.sqrt(1 - e2 * math.sin(lat) ** 2)
    T, C = math.tan(lat) ** 2, ep2 * math.cos(lat) ** 2
    A = math.cos(lat) * (lon - lon0)
    M = a * ((1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256) * lat
             - (3 * e2 / 8 + 3 * e2 ** 2 / 32 + 45 * e2 ** 3 / 1024) * math.sin(2 * lat)
             + (15 * e2 ** 2 / 256 + 45 * e2 ** 3 / 1024) * math.sin(4 * lat)
             - (35 * e2 ** 3 / 3072) * math.sin(6 * lat))
    E = k0 * N * (A + (1 - T + C) * A ** 3 / 6 + (5 - 18 * T + T * T + 72 * C - 58 * ep2) * A ** 5 / 120) + 500000
    Nn = k0 * (M + N * math.tan(lat) * (A * A / 2 + (5 - T + 9 * C + 4 * C * C) * A ** 4 / 24
                                        + (61 - 58 * T + T * T + 600 * C - 330 * ep2) * A ** 6 / 720))
    return E, Nn


def load_subdaily(path, hortho):
    """gnssrefl subdaily file -> per-arc table (same columns elevation_quality.py reads)."""
    rows = []
    for line in open(path, errors="replace"):
        if line.startswith("%") or not line.strip():
            continue
        c = line.split()
        if len(c) < 22:
            continue
        try:
            t = datetime(int(float(c[0])), int(float(c[17])), int(float(c[18])), int(float(c[19])),
                         int(float(c[20])), int(float(c[21])), tzinfo=timezone.utc)
            rows.append({"epoch": t.timestamp(), "rh": float(c[2]), "az": float(c[5]),
                         "emin": float(c[7]), "emax": float(c[8]), "freq": int(float(c[10])),
                         "pkn": float(c[13])})
        except (ValueError, IndexError):
            continue
    d = pd.DataFrame(rows)
    if len(d):
        d["wl"] = hortho - d["rh"]
        d["elev"] = (d["emin"] + d["emax"]) / 2
    return d


def marconi_reference(gauge_csv, epochs, fit=(1.24, -48.0, -0.10)):
    g = pd.read_csv(gauge_csv).dropna(subset=["level_navd88"]).sort_values("epoch")
    a, lag, b = fit
    t = np.asarray(epochs, float) - lag * 60
    ep, lv = g["epoch"].to_numpy(float), g["level_navd88"].to_numpy(float)
    out = a * np.interp(t, ep, lv) + b
    i = np.clip(np.searchsorted(ep, t), 1, len(ep) - 1)
    out[(t < ep[0]) | (t > ep[-1]) | (ep[i] - ep[i - 1] > 3 * 3600)] = np.nan
    return out


def wave_height(waves_csv, epochs):
    if not waves_csv or not Path(waves_csv).exists():
        return np.full(len(epochs), np.nan)
    w = pd.read_csv(waves_csv)
    col = "hs_best" if "hs_best" in w else "wvht_m"
    w = w[pd.to_numeric(w[col], errors="coerce").notna()]
    return np.interp(epochs, w["epoch"].astype(float), w[col].astype(float), left=np.nan, right=np.nan)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--subdaily", required=True, help="gnssrefl subdaily file (per-arc reflector heights)")
    ap.add_argument("--station-json", default="/home/argus_user/GNSS/v4.1/station/resources/station.json")
    ap.add_argument("--hortho", type=float, default=None, help="antenna height, m NAVD88 (default from station.json)")
    ap.add_argument("--antenna-en", nargs=2, type=float, default=None, metavar=("E", "N"),
                    help="antenna UTM 19N (default from station.json latitude/longitude)")
    ap.add_argument("--gauge-csv", default=str(HERE / "archive" / "gauge_8447435.csv"))
    ap.add_argument("--waves-csv", default=None, help="default archive/waves_marconi.csv, else waves_44008.csv")
    ap.add_argument("--rough", type=float, default=1.5, help="Hs (m) above which a sea counts as rough")
    ap.add_argument("--min-per-day", type=float, default=24.0,
                    help="arcs per day a mask must keep for the subdaily spline (default 24)")
    ap.add_argument("--margin", type=float, default=10.0,
                    help="m: a reflection point within this distance seaward of the waterline counts as 'beach'")
    args = ap.parse_args()

    cfg = json.loads(Path(args.station_json).read_text()) if Path(args.station_json).exists() else {}
    hortho = args.hortho if args.hortho is not None else cfg.get("gnssrefl_orthometric_height")
    if hortho is None:
        sys.exit("need --hortho (antenna m NAVD88) or a station.json with gnssrefl_orthometric_height")
    if args.antenna_en:
        E0, N0 = args.antenna_en
    elif cfg.get("latitude") is not None:
        E0, N0 = latlon_to_utm19(cfg["latitude"], cfg["longitude"])
    else:
        sys.exit("need --antenna-en E N or a station.json with latitude/longitude")

    d = load_subdaily(args.subdaily, float(hortho))
    if not len(d):
        sys.exit(f"no arcs read from {args.subdaily}")
    d["ref"] = marconi_reference(args.gauge_csv, d["epoch"])
    d = d[d["ref"].notna()].copy()
    if not len(d):
        sys.exit("no arcs inside the gauge record -- run fetch_tide_gauge.py --days 120 first")
    waves = args.waves_csv or (HERE / "archive" / "waves_marconi.csv"
                               if (HERE / "archive" / "waves_marconi.csv").exists()
                               else HERE / "archive" / "waves_44008.csv")
    d["hs"] = wave_height(waves, d["epoch"])
    d["resid"] = d["wl"] - d["ref"]

    # where did it reflect? horizontal distance RH/tan(e) along the azimuth
    rh_true = float(hortho) - d["ref"]
    dist = rh_true / np.tan(np.radians(d["elev"]))
    E = E0 + dist * np.sin(np.radians(d["az"]))
    N = N0 + dist * np.cos(np.radians(d["az"]))
    sea = np.radians(SHORE_BEARING + 90)
    seaward = (E - SHORE_ORIGIN[0]) * np.sin(sea) + (N - SHORE_ORIGIN[1]) * np.cos(sea)
    waterline = -d["ref"] / FORESHORE_SLOPE                 # 0 m contour moves seaward as the tide falls
    d["past_waterline_m"] = seaward - waterline
    d["on"] = np.where(d["past_waterline_m"] > args.margin, "sea", "beach")
    ant_seaward = (E0 - SHORE_ORIGIN[0]) * np.sin(sea) + (N0 - SHORE_ORIGIN[1]) * np.cos(sea)
    d["rough"] = d["hs"] > args.rough

    r = lambda x: float(np.sqrt(np.mean(np.square(x)))) if len(x) else np.nan
    print("=" * 72)
    print("  GNSS-IR REFLECTION AUDIT -- arcs against the Chatham-derived tide")
    print("=" * 72)
    print(f"  antenna {E0:.1f} E, {N0:.1f} N; {hortho:.3f} m NAVD88; {-ant_seaward:.0f} m landward of "
          f"the Jan 2025 0 m contour")
    print(f"  {len(d)} arcs, {pd.to_datetime(d.epoch.min(), unit='s'):%Y-%m-%d} to "
          f"{pd.to_datetime(d.epoch.max(), unit='s'):%Y-%m-%d}; rough (Hs > {args.rough} m): "
          f"{int(d.rough.sum())} arcs; reference itself is good to ~0.13 m")
    print()

    def table(title, key, bins, labels=None):
        print(f"  {title}")
        print(f"    {'':>12} {'calm: n':>8} {'bias':>7} {'RMS':>6}   {'rough: n':>8} {'bias':>7} {'RMS':>6}  {'on beach':>8}")
        groups = pd.cut(d[key], bins, labels=labels) if bins is not None else d[key]
        for g, x in d.groupby(groups, observed=True):
            c, s = x[~x.rough], x[x.rough]
            print(f"    {str(g):>12} {len(c):>8d} {c.resid.mean() if len(c) else np.nan:>+7.2f} {r(c.resid):>6.2f}"
                  f"   {len(s):>8d} {s.resid.mean() if len(s) else np.nan:>+7.2f} {r(s.resid):>6.2f}"
                  f"  {(x.on == 'beach').mean():>7.0%}")
        print()

    table("By arc mean elevation (deg)", "elev", [0, 7, 9, 11, 13, 15, 20, 30])
    table("By azimuth (deg; shore-normal is 81)", "az", [0, 25, 50, 75, 100, 125, 150, 175, 330, 360])
    table("By where the reflection point fell", "on", None)

    # candidate masks: keep arcs with emax <= e2 and az within the window
    days = max((d.epoch.max() - d.epoch.min()) / 86400, 1.0)
    rb = d[d.rough]
    if len(rb):
        hb = rb[rb.on == "beach"]
        print(f"  ROUGH SEAS: {len(hb) / len(rb):.0%} of rough-sea arcs reflected on the beach; they read "
              f"{hb.resid.mean() if len(hb) else np.nan:+.2f} m against the tide, the sea arcs "
              f"{rb[rb.on == 'sea'].resid.mean() if (rb.on == 'sea').any() else np.nan:+.2f} m.")
        print()
    print("  CANDIDATE MASKS (arcs kept / RMS against the tide, calm and rough)")
    print(f"    {'upper elev':>10} {'azimuths':>11} {'kept':>6} {'calm RMS':>9} {'rough RMS':>10} {'rough bias':>11}")
    best = None
    for e2 in (10, 12, 13, 15, 25):
        for lo, hi in ((353, 173), (15, 150), (25, 135), (35, 125)):
            az_ok = (d.az >= lo) | (d.az <= hi) if lo > hi else (d.az >= lo) & (d.az <= hi)
            k = d[(d.emax <= e2 + 0.01) & az_ok]
            if len(k) < 50:
                continue
            c, s = k[~k.rough], k[k.rough]
            print(f"    {e2:>10} {f'{lo}-{hi}':>11} {len(k):>6} {r(c.resid):>9.3f} {r(s.resid):>10.3f} "
                  f"{s.resid.mean() if len(s) else np.nan:>+11.2f}")
            score = r(c.resid) + (0.5 * r(s.resid) if len(s) >= 20 else 0)
            # the subdaily spline needs arcs all day: at least --min-per-day on average
            if len(k) / days >= args.min_per_day and (best is None or score < best[0]):
                best = (score, e2, lo, hi, len(k))
    print()
    if best:
        _, e2, lo, hi, n = best
        print(f"  RECOMMENDED: elevation 5-{e2} deg, azimuths {lo}-{hi} ({n} of {len(d)} arcs, "
              f"{n / days:.0f} a day).")
        print("  In station.json (keep a copy of the old file to revert):")
        print(f'      "gnssrefl_elevation_max": {e2},')
        az = [lo, 360, 0, hi] if lo > hi else [lo, hi]
        print(f'      "gnssrefl_azimuth_regions": {az},')
        print("  then rerun gnssir for the days to reprocess, and subdaily.")
    print()
    print("  Rough-sea arcs that read high AND reflected on the beach are beach elevations,")
    print("  not wave setup: they should not be exported as observed total water level.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())

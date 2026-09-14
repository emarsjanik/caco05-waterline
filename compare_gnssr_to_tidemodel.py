#!/usr/bin/env python3
"""
Compare GNSS-R Measured Water Level to the Tide Model
--------------------------------------------------------
Reads gnssrefl subdaily spline output (usgs_spline_out.txt) and a tide
model file, interpolates the model onto the GNSS-R timestamps, and
reports the difference.

WHY THIS MATTERS: the GNSS-R "quasi-sea-level" column is computed as
    water elevation = station orthometric height - reflector height
and orthometric height in the US is NAVD88 -- the same datum the
camera EO z-coordinates use. The global ocean tide models
(EOT20/GOT/FES), by contrast, are referenced to mean sea level or a
geoid. The mean difference between them over a long overlap is
therefore the MSL->NAVD88 offset at this site, measured rather than
looked up, and is exactly the value
    extract_elevation_contours.py --vertical-datum-offset
needs.

It is worth being precise about what that difference contains. It is
NOT purely a datum shift. It also absorbs:
  * any bias in the global tide model at this specific location
    (these are gridded ocean models, not a local gauge);
  * non-tidal water level -- storm surge, wind setup, seasonal steric
    effects -- which the astronomical model does not predict at all
    and the GNSS-R measurement does include;
  * any error in the station's assumed orthometric height, which
    propagates one-for-one into every GNSS-R elevation.
The MEAN over many tidal cycles is dominated by the datum term, since
surge is roughly zero-mean over weeks. The STANDARD DEVIATION reported
below is the more honest description of how well the model represents
this site, and the trend check flags whether the offset is stable or
drifting.

Usage:
    python3 compare_gnssr_to_tidemodel.py <gnssr_spline_file> <tide_model_file>

Example:
    python3 compare_gnssr_to_tidemodel.py \\
        /home/argus_user/GNSS/v4.1/products/refl_code/Files/usgs/usgs_spline_out.txt \\
        marconi_tides_sherwood.xlsx
"""

import sys
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import pandas as pd


GAP_SENTINEL = 999.0

HEIGHT_COL_HINTS = ("_heightm", "_height_m", "heightm")
TIME_COL_HINTS = ("time", "datetime", "date")


def load_gnssr_spline(path):
    """
    Parses gnssrefl subdaily spline output. Columns (1-indexed):
        1 MJD, 2 RH(m), 3 YYYY, 4 MM, 5 DD, 6 HH, 7 MM, 8 SS,
        9 quasi-sea-level(m)  [= station orthometric height - RH]
    Lines beginning with '%' are comments. Rows where column 9 equals
    999 are gap sentinels and are dropped (the file header documents
    this convention).

    Also extracts the station orthometric height from the header when
    present, so the datum assumption is surfaced rather than implied.
    """
    times, levels, rh = [], [], []
    hortho = None

    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("%"):
                if "Hortho" in s or "orthometric height" in s.lower():
                    for tok in s.replace(",", " ").split():
                        try:
                            val = float(tok)
                            if 0.0 < val < 1000.0:
                                hortho = val
                        except ValueError:
                            continue
                continue
            parts = s.split()
            if len(parts) < 9:
                continue
            try:
                yyyy, mm, dd = int(parts[2]), int(parts[3]), int(parts[4])
                hh, mi, ss = int(parts[5]), int(parts[6]), int(parts[7])
                level = float(parts[8])
                reflector = float(parts[1])
            except ValueError:
                continue
            if abs(level - GAP_SENTINEL) < 1e-6:
                continue
            times.append(datetime(yyyy, mm, dd, hh, mi, ss, tzinfo=timezone.utc))
            levels.append(level)
            rh.append(reflector)

    if not times:
        raise ValueError(f"No usable GNSS-R rows parsed from {path}")

    epoch = np.array([t.timestamp() for t in times])
    order = np.argsort(epoch)
    return epoch[order], np.array(levels)[order], np.array(rh)[order], hortho


def load_tide_model(path):
    """Loads the tide model, averaging all '<model>_heightm' columns."""
    suffix = Path(path).suffix.lower()
    df = pd.read_excel(path) if suffix in (".xlsx", ".xls") else pd.read_csv(path)

    time_col = None
    for c in df.columns:
        if any(h in str(c).lower() for h in TIME_COL_HINTS):
            time_col = c
            break
    if time_col is None:
        raise ValueError(f"No time column found in {path}. Columns: {list(df.columns)}")

    height_cols = [c for c in df.columns
                   if any(h in str(c).lower() for h in HEIGHT_COL_HINTS)]
    if not height_cols:
        raise ValueError(f"No '*_heightm' columns found in {path}. "
                         f"Columns: {list(df.columns)}")

    ts = pd.to_datetime(df[time_col], errors="coerce")
    valid = ts.notna()
    for c in height_cols:
        valid &= df[c].notna()
    df, ts = df[valid], ts[valid]

    if ts.dt.tz is None:
        ts = ts.dt.tz_localize("UTC")
    ts = ts.dt.tz_convert("UTC")

    epoch = np.array([t.timestamp() for t in ts])
    levels = df[height_cols].to_numpy(dtype=float).mean(axis=1)
    order = np.argsort(epoch)
    return epoch[order], levels[order], height_cols


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 compare_gnssr_to_tidemodel.py "
              "<gnssr_spline_file> <tide_model_file>")
        sys.exit(1)

    g_epoch, g_level, g_rh, hortho = load_gnssr_spline(sys.argv[1])
    t_epoch, t_level, height_cols = load_tide_model(sys.argv[2])

    def fmt(e):
        return datetime.fromtimestamp(e, tz=timezone.utc).strftime("%Y-%m-%d %H:%M")

    print("=" * 92)
    print("GNSS-R MEASURED WATER LEVEL  vs  TIDE MODEL")
    print("=" * 92)
    print(f"GNSS-R      : {len(g_epoch):5d} points   {fmt(g_epoch[0])} -> {fmt(g_epoch[-1])}")
    if hortho is not None:
        print(f"              station orthometric height (NAVD88) = {hortho:.3f} m")
        print(f"              water elevation = {hortho:.3f} - reflector height")
    print(f"Tide model  : {len(t_epoch):5d} points   {fmt(t_epoch[0])} -> {fmt(t_epoch[-1])}")
    print(f"              averaging {len(height_cols)} model(s): {height_cols}")
    print()

    lo = max(g_epoch[0], t_epoch[0])
    hi = min(g_epoch[-1], t_epoch[-1])
    if hi <= lo:
        print("NO TEMPORAL OVERLAP between the two datasets -- cannot compare.")
        sys.exit(1)

    mask = (g_epoch >= lo) & (g_epoch <= hi)
    ge, gl = g_epoch[mask], g_level[mask]
    tl = np.interp(ge, t_epoch, t_level)

    diff = gl - tl
    days = (hi - lo) / 86400.0

    print(f"Overlap     : {fmt(lo)} -> {fmt(hi)}  ({days:.1f} days, {len(ge)} matched points)")
    print()
    print("=" * 92)
    print("DIFFERENCE  (GNSS-R measured  minus  tide model)")
    print("=" * 92)
    print(f"  mean            : {diff.mean():+8.4f} m   <-- candidate --vertical-datum-offset")
    print(f"  median          : {np.median(diff):+8.4f} m")
    print(f"  std dev         : {diff.std():8.4f} m   <-- how well the model represents this site")
    print(f"  RMS             : {np.sqrt((diff**2).mean()):8.4f} m")
    print(f"  min / max       : {diff.min():+8.4f} / {diff.max():+8.4f} m")
    print(f"  5th / 95th pct  : {np.percentile(diff,5):+8.4f} / {np.percentile(diff,95):+8.4f} m")
    print()

    # Is the offset stable, or drifting? A drift would mean a single
    # constant offset is the wrong model.
    t_days = (ge - ge[0]) / 86400.0
    slope, intercept = np.polyfit(t_days, diff, 1)
    print(f"  linear trend    : {slope*365:+.4f} m/year "
          f"({slope*days:+.4f} m across this {days:.0f}-day window)")
    if abs(slope * days) > 0.05:
        print("    -> NOTE: the offset DRIFTS materially over the window. A single constant")
        print("       offset is then an approximation; investigate before relying on it.")
    else:
        print("    -> offset is stable across the window; a constant is a reasonable model.")
    print()

    # Month-by-month, to expose seasonal/steric structure.
    print("  Monthly breakdown:")
    months = np.array([datetime.fromtimestamp(e, tz=timezone.utc).strftime("%Y-%m") for e in ge])
    for m in sorted(set(months)):
        sel = months == m
        print(f"    {m}: n={sel.sum():5d}  mean {diff[sel].mean():+7.4f} m  "
              f"std {diff[sel].std():6.4f} m")
    print()

    print("=" * 92)
    print("HOW TO USE THIS")
    print("=" * 92)
    print("  The tide model is referenced to MSL/geoid; GNSS-R is NAVD88. To put tide-model")
    print("  elevations onto NAVD88, ADD the mean difference:")
    print()
    print(f"    python3 extract_elevation_contours.py <tide_file> "
          f"--vertical-datum-offset {diff.mean():.4f}")
    print()
    print("  Caveat: that mean also absorbs model bias at this site and any error in the")
    print("  assumed station orthometric height. Where GNSS-R covers the imagery directly,")
    print("  prefer it over the model -- it needs no offset at all.")
    print("=" * 92)


if __name__ == "__main__":
    main()


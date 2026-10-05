#!/usr/bin/env python3
"""
Re-Label Image Manifests From The ADCP, At Each Image's Mid-Time
==================================================================
The manifests were built by taking the NEAREST hourly ADCP record, and on
the 30-minute ties (every :30 image) the EARLIER one -- so half of all
training and holdback labels were the previous hour's wave height (a
17:30 frame carried 17:00's 1.71 m, while 18:00 read 1.87 m). This
re-labels every row by linear interpolation of the ADCP wave height at
the middle of the image's 10-minute exposure (filename epoch + 300 s).

Wave period stays the NEAREST record: the ADCP peak period is a spectral
bin, and interpolating between bins invents periods it never reported.

Works on both manifest layouts:
  filename,...,wave_height_m,wave_period_s     (pair_images_waves.py)
  id,H,T                                       (training lists)
and adds/updates a label_method column. Rows whose time is not bracketed
by two ADCP records at most --max-bracket hours apart are dropped (the
last frame of the record, 15:00 on 2025-03-10, is the only one).

Usage:
    python3 relabel_manifests.py manifest_c2*.csv owg_models/marconi-c2-clean-*.csv --in-place
    python3 relabel_manifests.py manifest_c2_bright_clean2.csv --output relabelled.csv
Undo: git checkout -- <file>
"""

import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def load_adcp(path, hcol, tcol):
    w = pd.read_csv(path)
    t = (pd.to_datetime(w["time"], utc=True) - pd.Timestamp("1970-01-01", tz="UTC")) // pd.Timedelta("1s")
    w = pd.DataFrame({"epoch": t, "h": pd.to_numeric(w[hcol], errors="coerce"),
                      "t": pd.to_numeric(w[tcol], errors="coerce")}).dropna(subset=["h"])
    return w.sort_values("epoch").reset_index(drop=True)


def relabel(df, adcp, offset_s, max_bracket_h):
    name = "filename" if "filename" in df else "id"
    hcol, tcol = ("wave_height_m", "wave_period_s") if "wave_height_m" in df else ("H", "T")
    ep = df[name].astype(str).str.split(".").str[0].astype(np.int64).to_numpy() + offset_s
    ae, ah = adcp["epoch"].to_numpy(), adcp["h"].to_numpy()
    i = np.clip(np.searchsorted(ae, ep), 1, len(ae) - 1)
    ok = (ep >= ae[0]) & (ep <= ae[-1]) & (ae[i] - ae[i - 1] <= max_bracket_h * 3600)
    new_h = np.where(ok, np.interp(ep, ae, ah), np.nan)
    j = np.clip(np.searchsorted(ae, ep), 0, len(ae) - 1)
    jm = np.clip(j - 1, 0, len(ae) - 1)
    near = np.where(np.abs(ae[jm] - ep) <= np.abs(ae[j] - ep), jm, j)
    old = df[hcol].to_numpy(float)
    df = df.copy()
    df[hcol] = new_h
    if tcol in df:
        df[tcol] = np.where(ok, adcp["t"].to_numpy()[near], np.nan)
    df["label_method"] = f"adcp linear at epoch+{offset_s:.0f}s"
    d = new_h - old
    return df, d[np.isfinite(d)], int((~ok).sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("manifests", nargs="+")
    ap.add_argument("--adcp", default=str(Path(__file__).resolve().parent / "sig1000_waves_ALL.csv"))
    ap.add_argument("--height-col", default="wh_4061")
    ap.add_argument("--period-col", default="wp_peak")
    ap.add_argument("--offset", type=float, default=300.0,
                    help="seconds after the filename epoch to label at (default 300: middle of 10 min)")
    ap.add_argument("--max-bracket", type=float, default=2.0, help="hours between bracketing records")
    ap.add_argument("--in-place", action="store_true")
    ap.add_argument("--output", default=None, help="with one manifest: where to write")
    args = ap.parse_args()
    if not args.in_place and not (args.output and len(args.manifests) == 1):
        sys.exit("give --in-place, or --output with a single manifest")
    adcp = load_adcp(args.adcp, args.height_col, args.period_col)
    for m in args.manifests:
        df = pd.read_csv(m)
        new, d, n_out = relabel(df, adcp, args.offset, args.max_bracket)
        if not len(d):
            print(f"{m}: no row inside the ADCP record -- left unchanged")
            continue
        out = m if args.in_place else args.output
        new = new[np.isfinite(new["wave_height_m" if "wave_height_m" in new else "H"])]   # trainers expect labels
        new.to_csv(out, index=False)
        print(f"{m}: {len(df)} rows, label change RMS {np.sqrt(np.mean(d ** 2)):.3f} m, "
              f"max {np.abs(d).max():.2f} m, {n_out} outside the ADCP record (dropped) -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

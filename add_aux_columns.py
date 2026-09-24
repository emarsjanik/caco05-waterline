#!/usr/bin/env python3
"""
Add Auxiliary Columns To OWG Split Files
===========================================
Joins extra measurements from the ADCP record -- water level in
particular -- onto the train/validation CSVs produced by
prepare_owg_data.py, matching on the capture time encoded in each
image id.

WHY WATER LEVEL: the crop previews showed that the most conspicuous
difference between frames is how far up the beach the water reaches
and how much wet sand is exposed, and that is governed by tide as much
as by waves. Without knowing the tide, the network has to infer it
from the scene before it can isolate any wave signal, spending
capacity on a quantity that is already measured. Supplying it directly
frees the network to use the image for what only the image can say.

Tide and wave height are close to uncorrelated, so this is not leakage
of the target by another route: the script reports the correlation so
that can be checked rather than assumed.

Usage:
    python3 add_aux_columns.py --splits ~/owg_marconi/marconi-c2-clean \\
        --waves sig1000_waves_ALL.csv --columns water_level
"""

import re
import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def epoch_from_id(i):
    m = re.match(r"^(\d{9,11})\.", str(i))
    return int(m.group(1)) if m else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", required=True,
                    help="Stem of the split files, e.g. '~/owg_marconi/marconi-c2-clean'. "
                         "Reads <stem>-train.csv and <stem>-val.csv and rewrites them.")
    ap.add_argument("--waves", required=True, help="The ADCP CSV.")
    ap.add_argument("--columns", nargs="+", default=["water_level"])
    ap.add_argument("--tz", default="UTC")
    ap.add_argument("--max-gap", type=float, default=60.0,
                    help="Minutes; rows further than this from a record are dropped.")
    args = ap.parse_args()

    w = pd.read_csv(args.waves)
    w["time"] = pd.to_datetime(w["time"])
    w["time"] = (w["time"].dt.tz_localize(args.tz) if w["time"].dt.tz is None
                 else w["time"].dt.tz_convert(args.tz))
    origin = pd.Timestamp("1970-01-01", tz="UTC")
    w["epoch"] = (w["time"].dt.tz_convert("UTC") - origin).dt.total_seconds().astype("int64")
    w = w.sort_values("epoch")
    missing = [c for c in args.columns if c not in w.columns]
    if missing:
        sys.exit(f"wave file has no column(s): {missing}")

    wep = w["epoch"].to_numpy()
    stem = Path(args.splits).expanduser()
    print("=" * 66)
    print("ADD AUXILIARY COLUMNS")
    print("=" * 66)

    for kind in ("train", "val"):
        path = Path(str(stem) + f"-{kind}.csv")
        if not path.exists():
            print(f"  {path} not found -- skipped")
            continue
        df = pd.read_csv(path)
        df.columns = [c.strip() for c in df.columns]
        df["_epoch"] = df["id"].map(epoch_from_id)
        if df["_epoch"].isna().any():
            sys.exit(f"{path}: could not read a timestamp from every id")

        idx = np.searchsorted(wep, df["_epoch"].to_numpy())
        idx = np.clip(idx, 1, len(wep) - 1)
        left, right = wep[idx - 1], wep[idx]
        pick = np.where(np.abs(df["_epoch"] - left) <= np.abs(df["_epoch"] - right),
                        idx - 1, idx)
        gap = np.abs(wep[pick] - df["_epoch"]) / 60.0
        for c in args.columns:
            df[c] = w[c].to_numpy()[pick]
        far = gap > args.max_gap
        if far.any():
            print(f"  {kind}: {int(far.sum())} row(s) beyond {args.max_gap:g} min -- dropped")
            df = df[~far]
        df = df.drop(columns=["_epoch"])
        df.to_csv(path, index=False)
        desc = ", ".join(f"{c} {df[c].min():.2f}..{df[c].max():.2f}" for c in args.columns)
        print(f"  {kind}: {len(df)} rows, added {desc}")

        # Tide must not simply stand in for the target. If it did, the
        # network could score well without using the image at all.
        if "H" in df.columns:
            for c in args.columns:
                r = df[c].corr(df["H"])
                note = "  -- FINE" if abs(r) < 0.3 else "  -- CHECK: this may substitute for the image"
                print(f"        corr({c}, H) = {r:+.2f}{note}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Which Days Was The Camera Aimed Differently? (Sea-Horizon Check)
===================================================================
For every day of imagery, finds the sea horizon in the photos and
compares it with where the 2025-02-19 calibration puts it. The horizon
is a fixed target (only the camera can move it), visible in almost
every daytime frame, and needs no snow-free beach or matching features.

Per day and camera it reports how far off the horizon is (median over
up to --per-day frames nearest 17:00 UTC) and the tilt and roll change
that best explains it. A camera that was knocked or re-aimed shows as a
run of days with the same large offset; those days need their own
calibration, which estimate_eo_rotation.py solves -- the exact command
is printed.

LIMITS. The horizon fixes tilt and roll. A pure left-right turn
(azimuth) barely moves it, so a small pan can pass unnoticed here;
estimate_eo_rotation.py solves all three angles. Haze or fog hides the
horizon: such days show as "unclear".

Usage:
    python3 horizon_check.py /mnt/I2Rgus_Data/Chelsea_calibration/work/original
"""

import re
import sys
import argparse
from pathlib import Path
from datetime import datetime, timezone
from collections import defaultdict

import numpy as np

from georectify import load_extrinsics, load_intrinsics
from estimate_eo_rotation import horizon_rows, observed_horizon

HERE = Path(__file__).resolve().parent
CAL = HERE / "calibration"
MOVED_PX = 12.0          # median horizon offset above this = different pointing
FIT_OK_PX = 8.0          # a tilt+roll fit must explain the horizon this well to trust it


def fit_tilt_roll(io, eo, cols, obs, iterations=15):
    """Tilt and roll change (deg) that best moves the predicted horizon onto obs.
    Gauss-Newton with robust (soft-L1) weights, numpy only: the station's scipy is
    too old for its numpy, so scipy.optimize cannot be imported there."""
    ok = np.isfinite(obs)
    p = np.zeros(2)

    def res(q):
        e = eo.copy()
        e[4] += np.deg2rad(q[0]); e[5] += np.deg2rad(q[1])
        r = horizon_rows(io, e, cols[ok]) - obs[ok]
        return np.where(np.isfinite(r), r, 200.0)

    for _ in range(iterations):
        r = res(p)
        J = np.column_stack([(res(p + d) - r) / 0.01 for d in (np.array([0.01, 0]), np.array([0, 0.01]))])
        w = 1.0 / np.sqrt(1.0 + (r / 5.0) ** 2)          # soft-L1, scale 5 px
        A, b = J * w[:, None], -r * w
        step = np.linalg.lstsq(A, b, rcond=None)[0]
        p = p + np.clip(step, -3, 3)
        if np.abs(step).max() < 1e-4:
            break
    return p, np.abs(res(p))


def main():
    import cv2
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("image_dir", help="folder of original *.timex.jpg frames")
    ap.add_argument("--eo-date", default="20250219", help="calibration to test against")
    ap.add_argument("--per-day", type=int, default=3)
    args = ap.parse_args()

    frames = defaultdict(list)
    for p in sorted(Path(args.image_dir).glob("*.timex.jpg")):
        m = re.match(r"^(\d{9,11})\.", p.name)
        cam = "c1" if ".c1." in p.name else "c2" if ".c2." in p.name else None
        if m and cam:
            t = datetime.fromtimestamp(int(m.group(1)), tz=timezone.utc)
            frames[(cam, t.strftime("%Y-%m-%d"))].append((abs(t.hour + t.minute / 60 - 17), p))
    if not frames:
        sys.exit(f"No *.timex.jpg in {args.image_dir}")

    moved = defaultdict(list)
    for cam in ("c1", "c2"):
        days = sorted(d for c, d in frames if c == cam)
        if not days:
            continue
        io = load_intrinsics(CAL / f"CACO05_{cam}_20240801_IO.yaml")
        eo = load_extrinsics(CAL / f"CACO05_{cam}_{args.eo_date}_EO.yaml")
        cols = np.arange(150, int(io[0]) - 150, 100)
        pred = horizon_rows(io, eo, cols)
        print(f"\n{cam}: sea horizon against the {args.eo_date} calibration")
        print("  date         frames  offset px   tilt change  roll change  fit px   verdict")
        for day in days:
            offs, fits = [], []
            for _, p in sorted(frames[(cam, day)])[:args.per_day]:
                g = cv2.imread(str(p), cv2.IMREAD_GRAYSCALE)
                if g is None:
                    continue
                obs = observed_horizon(g, cols, pred, window=250)
                (dt, dr), res = fit_tilt_roll(io, eo, cols, obs)
                offs.append(np.nanmedian(np.abs(obs - pred)))
                fits.append((dt, dr, np.median(res)))
            if not offs:
                continue
            off = float(np.median(offs))
            dt, dr, fres = (float(np.median([f[i] for f in fits])) for i in range(3))
            if fres > FIT_OK_PX:
                verdict = "unclear (haze? horizon not found cleanly)"
            elif off > MOVED_PX:
                verdict = "DIFFERENT POINTING"
                moved[cam].append(day)
            else:
                verdict = "matches"
            print(f"  {day}   {len(offs):4d}   {off:8.1f}   {dt:+9.2f} deg  {dr:+9.2f} deg  "
                  f"{fres:6.1f}   {verdict}")

    if not any(moved.values()):
        print("\nEvery day with a clear horizon matches the calibration.")
        return
    print("\nDays with a different pointing need their own calibration. For each run of such "
          "days, solve it from the photos (reference = days that match):")
    for cam, days in moved.items():
        runs, start, prev = [], days[0], days[0]
        for d in days[1:]:
            gap = (datetime.fromisoformat(d) - datetime.fromisoformat(prev)).days
            if gap > 1:
                runs.append((start, prev)); start = d
            prev = d
        runs.append((start, prev))
        for a, b in runs:
            print(f"  python3 {HERE / 'estimate_eo_rotation.py'} --camera {cam} "
                  f"--images-dir {args.image_dir} --first {a} --last {b} --add-to-setups")
    print("Check each run's dates against the table first: an 'unclear' day inside or next to a "
          "run may belong to it.")


if __name__ == "__main__":
    main()

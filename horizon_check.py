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
calibration. --write-setups writes one per run of days with the same
pointing (tilt and roll from the horizon fit, azimuth and position kept)
and lists them in calibration/chelsea_setups.csv, which
process_chelsea.py reads. --against-setups then checks every day against
its own calibration: all should match. estimate_eo_rotation.py solves
all three angles instead, when photos with a known pointing exist.

LIMITS. The horizon fixes tilt and roll. A pure left-right turn
(azimuth) barely moves it, so a small pan can pass unnoticed here;
estimate_eo_rotation.py solves all three angles. Haze or fog hides the
horizon: such days show as "unclear".

Usage:
    python3 horizon_check.py /mnt/I2Rgus_Data/Chelsea_calibration/work/original
    python3 horizon_check.py <same folder> --write-setups
    python3 horizon_check.py <same folder> --against-setups
"""

import re
import sys
import argparse
from pathlib import Path
from datetime import datetime, timezone
from collections import defaultdict

import numpy as np

from georectify import load_extrinsics, load_intrinsics
from estimate_eo_rotation import horizon_rows, observed_horizon, write_eo

HERE = Path(__file__).resolve().parent
CAL = HERE / "calibration"
MOVED_PX = 12.0          # median horizon offset above this = different pointing
FIT_OK_PX = 8.0          # a tilt+roll fit must explain the horizon this well to trust it
RUN_TOL_DEG = 0.3        # days whose tilt and roll agree this well share one pointing


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
    ap.add_argument("--write-setups", action="store_true",
                    help="write a calibration for each run of days with a different pointing "
                         "(tilt and roll from the horizon, azimuth kept) and list them in "
                         "calibration/chelsea_setups.csv for process_chelsea.py")
    ap.add_argument("--against-setups", action="store_true",
                    help="compare each day with ITS calibration from calibration/chelsea_setups.csv "
                         "(a check after --write-setups: every day should then match)")
    args = ap.parse_args()
    setups = []
    if args.against_setups:
        import csv
        sc = CAL / "chelsea_setups.csv"
        if sc.exists():
            with open(sc, newline="") as f:
                setups = list(csv.DictReader(f))

    frames = defaultdict(list)
    for p in sorted(Path(args.image_dir).glob("*.timex.jpg")):
        m = re.match(r"^(\d{9,11})\.", p.name)
        cam = "c1" if ".c1." in p.name else "c2" if ".c2." in p.name else None
        if m and cam:
            t = datetime.fromtimestamp(int(m.group(1)), tz=timezone.utc)
            frames[(cam, t.strftime("%Y-%m-%d"))].append((abs(t.hour + t.minute / 60 - 17), p))
    if not frames:
        sys.exit(f"No *.timex.jpg in {args.image_dir}")

    moved = defaultdict(list)          # cam -> [(day, dtilt, droll)]
    for cam in ("c1", "c2"):
        days = sorted(d for c, d in frames if c == cam)
        if not days:
            continue
        io = load_intrinsics(CAL / f"CACO05_{cam}_20240801_IO.yaml")
        eo_default = load_extrinsics(CAL / f"CACO05_{cam}_{args.eo_date}_EO.yaml")
        cols = np.arange(150, int(io[0]) - 150, 100)
        against = "each day's own calibration (chelsea_setups.csv)" if args.against_setups \
            else f"the {args.eo_date} calibration"
        print(f"\n{cam}: sea horizon against {against}")
        print("  date         frames  offset px   tilt change  roll change  fit px   verdict")
        for day in days:
            eo = eo_default
            for r in setups:
                if r["camera"] == cam and r["first_date"] <= day <= r["last_date"]:
                    eo = load_extrinsics(CAL / r["eo_file"])
            pred = horizon_rows(io, eo, cols)
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
                moved[cam].append((day, dt, dr))
            else:
                verdict = "matches"
            print(f"  {day}   {len(offs):4d}   {off:8.1f}   {dt:+9.2f} deg  {dr:+9.2f} deg  "
                  f"{fres:6.1f}   {verdict}")

    if not any(moved.values()):
        print("\nEvery day with a clear horizon matches the calibration.")
        return

    # Runs of consecutive days with the SAME pointing (a new run starts at a
    # gap of more than a day or a change of more than RUN_TOL_DEG).
    runs = []
    for cam, days in moved.items():
        cur = [days[0]]
        for d in days[1:]:
            gap = (datetime.fromisoformat(d[0]) - datetime.fromisoformat(cur[-1][0])).days
            same = abs(d[1] - cur[-1][1]) < RUN_TOL_DEG and abs(d[2] - cur[-1][2]) < RUN_TOL_DEG
            if gap > 1 or not same:
                runs.append((cam, cur)); cur = [d]
            else:
                cur.append(d)
        runs.append((cam, cur))
    print("\nPointings found (each run: median tilt and roll change against the calibration):")
    rows = []
    for cam, run in runs:
        a, b = run[0][0], run[-1][0]
        dt = float(np.median([d[1] for d in run])); dr = float(np.median([d[2] for d in run]))
        name = f"CACO05_{cam}_{a}_to_{b}_EO.yaml"
        print(f"  {cam}  {a} to {b}  ({len(run)} day(s))  tilt {dt:+.2f} deg, roll {dr:+.2f} deg  -> {name}")
        rows.append((cam, a, b, name, dt, dr))
    print("  The horizon fixes tilt and roll; the azimuth (left-right pan) is kept from the "
          f"{args.eo_date} calibration. A pan error moves points mostly ALONG the beach for these "
          "cameras, little across it.")
    if not args.write_setups:
        print("\nTo write these calibrations and use them in process_chelsea.py, rerun with "
              "--write-setups.")
        return
    import csv
    sc = CAL / "chelsea_setups.csv"
    old = []
    if sc.exists():
        with open(sc, newline="") as f:
            old = [r for r in csv.DictReader(f) if r["camera"] not in {r_[0] for r_ in rows}]
    for cam, a, b, name, dt, dr in rows:
        eo = load_extrinsics(CAL / f"CACO05_{cam}_{args.eo_date}_EO.yaml")
        ang = np.degrees(eo[3:6]) + [0.0, dt, dr]
        write_eo(CAL / name, eo, ang,
                 f"horizon_check.py: CACO05_{cam}_{args.eo_date}_EO.yaml with tilt {dt:+.3f} and "
                 f"roll {dr:+.3f} deg fitted to the sea horizon, {a} to {b}; azimuth and position kept.")
    with open(sc, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["camera", "first_date", "last_date", "eo_file"])
        w.writeheader()
        w.writerows(old)
        w.writerows({"camera": c, "first_date": a, "last_date": b, "eo_file": n} for c, a, b, n, _, _ in rows)
    print(f"\nWrote {len(rows)} calibration(s) and {sc}. Check: rerun horizon_check.py -- each "
          f"day should then show a small offset against ITS setup (--against-setups) -- and "
          f"then run process_chelsea.py.")

if __name__ == "__main__":
    main()

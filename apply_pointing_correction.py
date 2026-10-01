#!/usr/bin/env python3
"""
Carry A Survey-Fitted Pointing Correction To Every Historical Period
=======================================================================
fit_eo_to_survey.py solved one period's pointing against a lidar survey
(Jan 18-23 2025 against the 23 Jan lidar). The error it found -- ~22 deg
of pan on both cameras -- sits in the 2025-02-19 calibration that every
historical period starts from: horizon_check.py only changes tilt and
roll on top of it, period by period. So the same correction applies to
all periods: each period's pointing = its horizon-fitted pointing + the
(azimuth, tilt, roll) change the survey fit made to the period it fitted.

Writes CACO05_<cam>_<first>_to_<last>_corr_EO.yaml for every row of
calibration/chelsea_setups.csv for the camera (except the fitted period
itself, already correct) and points the rows at them. The default
pointing for dates not listed (2025-02-19) gets a corrected copy too,
CACO05_<cam>_20250219_corr_EO.yaml, and a row for each gap is added.

CHECK IT. The correction is measured on one period and assumed for the
others. Any period with a survey of its own should be checked -- e.g.
Feb 17-Mar 10 against the 6 Mar 2025 lidar with detect_original_view.py
and compare_dem_survey.py.

A correction can be limited to some periods with --dates, e.g. when another
survey gave a better fit for them: each run starts again from the periods'
original horizon-fitted pointings, so corrections never stack.

Usage:
    python3 apply_pointing_correction.py --camera c2 \\
        --fitted CACO05_c2_2025-01-18_to_2025-01-23_lidar_EO.yaml \\
        --fitted-from CACO05_c2_2025-01-08_to_2025-01-23_EO.yaml
"""

import sys
import csv
import argparse
from pathlib import Path
from datetime import date, timedelta

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from georectify import load_extrinsics                  # noqa: E402
from estimate_eo_rotation import write_eo               # noqa: E402

CAL = HERE / "calibration"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--camera", required=True, choices=["c1", "c2"])
    ap.add_argument("--fitted", required=True, help="EO written by fit_eo_to_survey.py")
    ap.add_argument("--fitted-from", required=True,
                    help="the EO that fit started from (the period's horizon-fitted pointing)")
    ap.add_argument("--base", default=None, help="default pointing (default CACO05_<cam>_20250219_EO.yaml)")
    ap.add_argument("--dates", nargs=2, default=None, metavar=("FIRST", "LAST"),
                    help="only correct periods overlapping these dates (e.g. a fit from the Mar "
                         "lidar for 2025-01-24 2025-03-31); periods outside keep what they have")
    ap.add_argument("--span", nargs=2, default=["2024-10-01", "2025-03-31"],
                    help="dates the historical imagery covers; gaps between listed periods in "
                         "this span get rows with the corrected default pointing")
    args = ap.parse_args()

    p = lambda n: Path(n) if Path(n).exists() else CAL / n
    fitted, start = load_extrinsics(p(args.fitted)), load_extrinsics(p(args.fitted_from))
    delta = np.degrees(fitted[3:6] - start[3:6])
    print(f"correction from {Path(args.fitted).name}: azimuth {delta[0]:+.2f}, tilt {delta[1]:+.2f}, "
          f"roll {delta[2]:+.2f} deg")

    sc = CAL / "chelsea_setups.csv"
    rows = list(csv.DictReader(open(sc, newline=""))) if sc.exists() else []
    fitted_name = Path(args.fitted).name

    def corrected(src_name, out_name, note):
        eo = load_extrinsics(p(src_name))
        write_eo(CAL / out_name, eo, np.degrees(eo[3:6]) + delta,
                 f"apply_pointing_correction.py: {src_name} + ({delta[0]:+.3f}, {delta[1]:+.3f}, "
                 f"{delta[2]:+.3f}) deg from {fitted_name}; {note}")
        return out_name

    lo, hi = args.dates if args.dates else ("0000", "9999")

    def original(name):
        """The horizon-fitted pointing a corrected file was made from."""
        return name.replace("_corr_EO", "_EO")

    out = []
    for r in rows:
        overlaps = r["first_date"] <= hi and r["last_date"] >= lo
        if r["camera"] != args.camera or not overlaps or r["eo_file"] == fitted_name \
                or "_lidar_EO" in r["eo_file"]:
            out.append(r)
            continue
        src = original(r["eo_file"])
        if src.startswith(f"CACO05_{args.camera}_20250219"):
            out.append(r)              # a gap row: refilled below
            continue
        name = f"CACO05_{args.camera}_{r['first_date']}_to_{r['last_date']}_corr_EO.yaml"
        corrected(src, name, f"{r['first_date']} to {r['last_date']}")
        print(f"  {r['first_date']} .. {r['last_date']}: {src} -> {name}")
        out.append(dict(r, eo_file=name))

    # gaps: dates in the span not covered by any row for this camera use the default pointing
    base = args.base or f"CACO05_{args.camera}_20250219_EO.yaml"
    tag = f"_{lo}_to_{hi}" if args.dates else ""
    base_corr = corrected(base, f"CACO05_{args.camera}_20250219{tag}_corr_EO.yaml", "default pointing")
    is_gap = lambda r: (r["camera"] == args.camera and original(r["eo_file"]).startswith(
        f"CACO05_{args.camera}_20250219") and r["first_date"] <= hi and r["last_date"] >= lo)
    out = [r for r in out if not is_gap(r)]
    cover = sorted((r["first_date"], r["last_date"]) for r in out if r["camera"] == args.camera)
    d = max(date.fromisoformat(args.span[0]), date.fromisoformat(lo) if args.dates else date.min)
    end = min(date.fromisoformat(args.span[1]), date.fromisoformat(hi) if args.dates else date.max)
    gaps, g0 = [], None
    while d <= end:
        iso = d.isoformat()
        inside = any(a <= iso <= b for a, b in cover)
        if not inside and g0 is None:
            g0 = iso
        if inside and g0 is not None:
            gaps.append((g0, (d - timedelta(days=1)).isoformat())); g0 = None
        d += timedelta(days=1)
    if g0 is not None:
        gaps.append((g0, end.isoformat()))
    for a, b in gaps:
        out.append({"camera": args.camera, "first_date": a, "last_date": b, "eo_file": base_corr})
        print(f"  {a} .. {b}: (default) -> {base_corr}")

    out.sort(key=lambda r: (r["camera"], r["first_date"]))
    with open(sc, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["camera", "first_date", "last_date", "eo_file"])
        w.writeheader()
        w.writerows(out)
    print(f"\nwrote {sc}: every {args.camera} date in {args.span[0]}..{args.span[1]} now has a "
          f"corrected pointing. Check a period with its own survey before relying on it.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Carry A Survey-Fitted Pointing Correction To Other Periods Of The Same Aim
=============================================================================
fit_eo_to_survey.py solved one period's pointing against a lidar survey
(Jan 18-23 2025 against the 23 Jan lidar): ~22 deg of pan on both
cameras. That is not an error in the 2025-02-19 calibration. The cameras
were re-aimed on 24 Jan 2025, the day after the lidar flew (station ID
CACO03 -> CACO04 that evening; horizon_check.py sees new tilt and roll
from that day), and the calibration was made on 19 Feb with the NEW aim.
Every period before 24 Jan inherited a calibration of a different aim,
and horizon_check.py, which only changes tilt and roll, cannot see a pan.
Evidence: the 23 Jan lidar drawn on 22 Jan and 25 Jan photos
(project_survey.py) fits the fitted pointing before the re-aim and the
calibration after it; c2's view gains the dune in its lower left from
25 Jan; the Mar 2025 lidar matches the calibration's pan to a few
degrees, and the Nov 2025 target calibration agrees with it to ~3 deg.

So the January correction belongs to the old aim only: each period's
pointing = its horizon-fitted pointing + the (azimuth, tilt, roll)
change the survey fit made. Without --dates it is applied from the start
of --span to 2025-01-23, the last day of the old aim; later periods are
left as they are. Periods of the new aim get their own fit (the Mar 2025
lidar) with --dates 2025-01-24 2025-03-31.

Writes CACO05_<cam>_<first>_to_<last>_corr_EO.yaml for every row of
calibration/chelsea_setups.csv for the camera in the dates (except the
fitted period itself, already correct) and points the rows at them.
Dates in the window that no row lists get a corrected copy of the default
pointing (2025-02-19) and a row of their own.

CHECK IT. The correction is measured on one period and assumed for the
others of the same aim. Check a period against a survey of its own, or
draw a survey on its photos with project_survey.py.

Each run starts again from the periods' original horizon-fitted
pointings, so corrections never stack.

Usage:
    python3 apply_pointing_correction.py --camera c2 \\
        --fitted CACO05_c2_2025-01-18_to_2025-01-23_lidar_EO.yaml \\
        --fitted-from CACO05_c2_2025-01-08_to_2025-01-23_EO.yaml
    python3 apply_pointing_correction.py --camera c2 \\
        --fitted CACO05_c2_2025-03-01_to_2025-03-10_lidar_EO.yaml \\
        --fitted-from CACO05_c2_2025-02-17_to_2025-03-10_EO.yaml \\
        --dates 2025-01-24 2025-03-31
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
# Last day of the aim the Jan 2025 lidar fit measured: the cameras were re-aimed
# on 24 Jan 2025 (CACO03 -> CACO04), so a correction without --dates stops here.
OLD_AIM_LAST_DAY = "2025-01-23"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--camera", required=True, choices=["c1", "c2"])
    ap.add_argument("--fitted", required=True, help="EO written by fit_eo_to_survey.py")
    ap.add_argument("--fitted-from", required=True,
                    help="the EO that fit started from (the period's horizon-fitted pointing)")
    ap.add_argument("--base", default=None, help="default pointing (default CACO05_<cam>_20250219_EO.yaml)")
    ap.add_argument("--dates", nargs=2, default=None, metavar=("FIRST", "LAST"),
                    help="only correct periods overlapping these dates (e.g. a fit from the Mar "
                         "lidar for 2025-01-24 2025-03-31); periods outside keep what they have. "
                         f"Default: the start of --span to {OLD_AIM_LAST_DAY}, the last day of the "
                         "aim the Jan 2025 lidar fit measured")
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

    dates_given = args.dates is not None
    lo, hi = args.dates if dates_given else (args.span[0], OLD_AIM_LAST_DAY)
    if not dates_given:
        print(f"no --dates: correcting {lo} .. {hi} only (the cameras were re-aimed on 24 Jan 2025)")

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
        if r["first_date"] < lo or r["last_date"] > hi:
            print(f"  WARNING: {r['first_date']} .. {r['last_date']} runs past {lo} .. {hi}; "
                  f"the whole period is corrected")
        name = f"CACO05_{args.camera}_{r['first_date']}_to_{r['last_date']}_corr_EO.yaml"
        corrected(src, name, f"{r['first_date']} to {r['last_date']}")
        print(f"  {r['first_date']} .. {r['last_date']}: {src} -> {name}")
        out.append(dict(r, eo_file=name))

    # gaps: dates in the span not covered by any row for this camera use the default pointing
    base = args.base or f"CACO05_{args.camera}_20250219_EO.yaml"
    tag = f"_{lo}_to_{hi}" if dates_given else ""
    base_corr = corrected(base, f"CACO05_{args.camera}_20250219{tag}_corr_EO.yaml", "default pointing")
    # the refill window: the dates this run corrects that the span covers
    w0 = max(date.fromisoformat(args.span[0]), date.fromisoformat(lo))
    w1 = min(date.fromisoformat(args.span[1]), date.fromisoformat(hi))
    ws, we = w0.isoformat(), w1.isoformat()
    is_gap = lambda r: (r["camera"] == args.camera and original(r["eo_file"]).startswith(
        f"CACO05_{args.camera}_20250219") and r["first_date"] <= we and r["last_date"] >= ws)
    if w0 <= w1:
        kept = []
        for r in out:
            if not is_gap(r):
                kept.append(r)
                continue
            # only the part inside the window is refilled; the rest keeps its pointing
            if r["first_date"] < ws:
                kept.append(dict(r, last_date=(w0 - timedelta(days=1)).isoformat()))
            if r["last_date"] > we:
                kept.append(dict(r, first_date=(w1 + timedelta(days=1)).isoformat()))
        out = kept
    cover = sorted((r["first_date"], r["last_date"]) for r in out if r["camera"] == args.camera)
    d, end = w0, w1
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
    if w0 <= w1:
        print(f"\nwrote {sc}: every {args.camera} date in {ws}..{we} now has a corrected pointing. "
              f"Check a period with its own survey before relying on it.")
    else:
        print(f"\nwrote {sc}: listed {args.camera} periods overlapping {lo}..{hi} corrected; "
              f"--span does not reach those dates, so no rows were added for unlisted days.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Camera Pointing From A Lidar Survey
======================================
Solves a camera's pan, tilt and roll for a period of historical imagery
by making its waterlines land on a survey (lidar DSM) where the sand is
exactly at each line's water level.

WHY. compare_dem_survey.py against the 23 Jan 2025 lidar (Oct 2026) found
the Jan 18-23 waterlines displaced horizontally -- c1 by ~15-18 m
landward, c2 onto the berm -- not a vertical datum shift, although the
overlays show the lines on the water's edge. Those frames were redrawn
into today's view with pointings from the sea horizon (horizon_check.py),
which fixes tilt and roll but says nothing about pan; c2's lines sit
rotated ~15-20 deg about the camera. A survey made during the period
constrains pan too: every waterline point at water level z must lie on
the survey's z contour.

HOW. process_chelsea.py detected the lines in frames resampled into the
2025-11-13 view, assuming the period's pointing (calibration/
chelsea_setups.csv). For each point that resampling is undone: the
current-view pixel, at its water level, goes back to the ORIGINAL pixel
under the assumed pointing. Then, for a trial pointing, the original
pixel's ray meets the plane z at a ground position, and the survey
elevation there should equal z. The cost is the mean |survey - z| over
all points, capped at 1 m (and 1 m off the survey) so lines that end up
on water or outside the survey neither help nor hurt. Pan, tilt and roll
changes are found by a coarse-to-fine grid search, pan over +/-25 deg and
tilt and roll within +/-1.5 deg of the horizon fit; position and lens are
kept.

Each (camera, pointing) group is solved separately. The result is
written as a new EO file; with --write-setups the group's row in
chelsea_setups.csv is pointed at it, so process_chelsea.py uses it.

CHECK THE ANSWER. The report shows the cost before and after, the
per-frame offsets, and how sharply the cost rises when each angle moves
0.25 deg from the optimum: a flat direction (often pan) is not
determined by these data and should not be trusted.

Usage:
    python3 fit_eo_to_survey.py /mnt/I2Rgus_Data/Chelsea_calibration/contour_points_ground.csv \\
        --survey /mnt/I2Rgus_Data/Chelsea_calibration/2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif \\
        --start-date 2025-01-18 --end-date 2025-01-23 [--camera c1] [--write-setups]
"""

import sys
import csv
import argparse
from pathlib import Path
from datetime import datetime, timezone
from collections import defaultdict

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from georectify import load_extrinsics, load_intrinsics, pixel_to_ground   # noqa: E402
from view_reproject import dst_to_src_points                              # noqa: E402
from compare_dem_survey import read_survey, sample                        # noqa: E402
from estimate_eo_rotation import write_eo                                 # noqa: E402

CAL = HERE / "calibration"
EO_NEW = "20251113"
CAP = 1.0


def cost_fn(io, eo0, u0, v0, z, survey):
    grid, x0, y0, cell = survey

    def cost(d_deg, per_point=False):
        eo = eo0.copy()
        eo[3:6] += np.deg2rad(d_deg)
        E, N = pixel_to_ground(u0, v0, z, io, eo)
        s = sample(grid, x0, y0, cell, E, N)
        r = s - z
        a = np.where(np.isfinite(r), np.minimum(np.abs(r), CAP), CAP)
        return (a, r) if per_point else float(a.mean())
    return cost


# Search ranges (deg): pan is wide because nothing else constrains it -- the
# sea horizon fixes tilt and roll to a few tenths of a degree but says nothing
# about pan, and a re-aimed camera can be turned 10-20 deg. Tilt and roll stay
# within MAX_PLAUSIBLE_DEG of the horizon fit.
SPAN = (25.0, 1.5, 1.5)
MAX_PLAUSIBLE_DEG = 1.5


def search(cost, span=SPAN):
    """Coarse-to-fine grid search over (d_azimuth, d_tilt, d_roll) in degrees."""
    best = np.zeros(3)
    lo, hi = -np.array(span), np.array(span)
    stages = ((np.array([1.0, 0.25, 0.25]), np.array(span)),
              (np.array([0.2, 0.05, 0.1]), np.array([1.0, 0.25, 0.25])),
              (np.array([0.04, 0.01, 0.02]), np.array([0.2, 0.05, 0.1])))
    for step, half in stages:
        axes = [np.arange(max(lo[k], best[k] - half[k]), min(hi[k], best[k] + half[k]) + step[k] / 2,
                          step[k]) for k in range(3)]
        best_c = np.inf
        for a in axes[0]:
            for t in axes[1]:
                for r in axes[2]:
                    c = cost(np.array([a, t, r]))
                    if c < best_c:
                        best_c, cand = c, np.array([a, t, r])
        best = cand
    return best, best_c


def main():
    from process_chelsea import load_setups, setup_of
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("contours", help="georectified contour points from process_chelsea.py")
    ap.add_argument("--survey", required=True, help="lidar DSM (.tif or .asc), NAVD88")
    ap.add_argument("--start-date"); ap.add_argument("--end-date")
    ap.add_argument("--camera", default="both")
    ap.add_argument("--eo", default=None,
                    help="Pointing the frames were resampled with (default: from "
                         "chelsea_setups.csv per frame)")
    ap.add_argument("--per-frame", type=int, default=200,
                    help="Points used per frame, spread along the line (default 200)")
    ap.add_argument("--force", action="store_true",
                    help="Offer flagged fits to --write-setups too -- only after checking them on "
                         "the photo (project_survey.py)")
    ap.add_argument("--write-setups", action="store_true",
                    help="Point the group's rows in calibration/chelsea_setups.csv at the new EO")
    args = ap.parse_args()

    survey = read_survey(args.survey)
    setups = load_setups()
    groups = defaultdict(lambda: defaultdict(list))      # (cam, eo) -> frame -> rows
    with open(args.contours, newline="") as f:
        for r in csv.DictReader(f):
            day = r.get("capture_time_utc", "")[:10]
            if (args.start_date and day < args.start_date) or (args.end_date and day > args.end_date):
                continue
            cam = r["camera"]
            if args.camera != "both" and cam != args.camera:
                continue
            eo = args.eo or setup_of(cam, int(float(r["capture_epoch"])), setups)
            groups[(cam, eo)][r["source_file"]].append(r)
    if not groups:
        sys.exit("no waterline points in that date range")

    new_setups = {}
    for (cam, eo_name), frames in sorted(groups.items()):
        io = load_intrinsics(CAL / f"CACO05_{cam}_20240801_IO.yaml")
        eo_path = Path(eo_name) if Path(eo_name).exists() else CAL / eo_name
        eo0 = load_extrinsics(eo_path)
        eo_new = load_extrinsics(CAL / f"CACO05_{cam}_{EO_NEW}_EO-CV.yaml")
        U, V, Z, F = [], [], [], []
        for k, (name, rows) in enumerate(sorted(frames.items())):
            idx = np.linspace(0, len(rows) - 1, min(args.per_frame, len(rows))).astype(int)
            for i in idx:
                U.append(float(rows[i]["pixel_column"])); V.append(float(rows[i]["pixel_row"]))
                Z.append(float(rows[i].get("beach_elevation_navd88") or rows[i]["tide_elevation_navd88"]))
                F.append(k)
        U, V, Z, F = map(np.array, (U, V, Z, F))
        # undo the resampling: current-view pixel -> original pixel under the assumed pointing
        u0, v0, ok = dst_to_src_points(U, V, Z, io, eo0, eo_new)
        u0, v0, Z, F = u0[ok], v0[ok], Z[ok], F[ok]
        print(f"\n{cam}, resampled with {eo_path.name}: {len(frames)} frames, {len(Z)} points")
        cost = cost_fn(io, eo0, u0, v0, Z, survey)
        a0, r0 = cost(np.zeros(3), per_point=True)
        on0 = np.isfinite(r0)
        print(f"  before: mean |survey - z| {a0.mean():.3f} m (capped at {CAP:g}), "
              f"{100 * on0.mean():.0f}% of points on the survey, median offset "
              f"{np.median(r0[on0]) if on0.any() else np.nan:+.2f} m")
        best, c = search(cost)
        a1, r1 = cost(best, per_point=True)
        on1 = np.isfinite(r1)
        print(f"  after : mean |survey - z| {c:.3f} m, {100 * on1.mean():.0f}% on the survey, "
              f"median offset {np.median(r1[on1]) if on1.any() else np.nan:+.2f} m")
        print(f"  change: azimuth {best[0]:+.2f} deg, tilt {best[1]:+.2f} deg, roll {best[2]:+.2f} deg")
        # sharpness: cost increase for a 0.25 deg move along each axis
        sharp = []
        for k, nm in enumerate(("azimuth", "tilt", "roll")):
            d = np.zeros(3); d[k] = 0.25
            sharp.append(min(cost(best + d), cost(best - d)) - c)
            flag = "  <-- FLAT: not determined by these data" if sharp[-1] < 0.01 else ""
            print(f"    {nm:8s}: cost +{sharp[-1]:.3f} m for 0.25 deg{flag}")
        # Is the answer believable? Tilt and roll are fixed to a few tenths of a
        # degree by the sea horizon (1 deg moves the beach 1-5 m here); pan is
        # free, and moves lines by about range x angle (1 deg = ~1-2.5 m at 60-150 m).
        E0, N0 = pixel_to_ground(u0, v0, Z, io, eo0)
        e1 = eo0.copy(); e1[3:6] += np.deg2rad(best)
        E1, N1 = pixel_to_ground(u0, v0, Z, io, e1)
        moved = np.hypot(E1 - E0, N1 - N0)
        print(f"  the fit moves the lines by median {np.nanmedian(moved):.1f} m "
              f"(p90 {np.nanpercentile(moved, 90):.1f} m)")
        problems = []
        for k, nm in enumerate(("azimuth", "tilt", "roll")):
            if abs(best[k]) >= SPAN[k] - 0.02:
                problems.append(f"{nm} hit the search limit ({best[k]:+.2f} deg): no real optimum"
                                + (" -- tilt and roll cannot move further without contradicting "
                                   "the sea horizon" if k else ""))
        if on1.mean() < on0.mean() - 0.10:
            problems.append(f"the fit pushed lines off the survey ({100 * on0.mean():.0f}% -> "
                            f"{100 * on1.mean():.0f}% on it): it improves by escaping, not matching")
        if problems:
            print("  NOT TRUSTWORTHY:")
            for pr in problems:
                print(f"    - {pr}")
            print("    Check the overlays: if the lines are on the water, the camera position or "
                  "lens may be off too (this fit keeps both); if not, the detections are wrong.")
        per = []
        for k in range(len(frames)):
            m = (F == k) & on1
            if m.sum() >= 10:
                per.append(np.median(r1[m]))
        if per:
            per = np.array(per)
            print(f"  per-frame offsets after: median {np.median(per):+.3f} m, p10-p90 "
                  f"{np.percentile(per, 10):+.2f} to {np.percentile(per, 90):+.2f} m "
                  f"({len(per)} frames)  <- what remains is datum + waves, not pointing")
        ang = np.degrees(eo0[3:6]) + best
        days = sorted({r[0]["capture_time_utc"][:10] for r in frames.values()})
        out = CAL / f"CACO05_{cam}_{days[0]}_to_{days[-1]}_lidar_EO.yaml"
        write_eo(out, eo0, ang, f"fit_eo_to_survey.py: {eo_path.name} with azimuth {best[0]:+.3f}, "
                 f"tilt {best[1]:+.3f}, roll {best[2]:+.3f} deg fitted to {Path(args.survey).name} "
                 f"using {len(frames)} waterline frames {days[0]} to {days[-1]}; position kept.")
        print(f"  wrote {out}")
        if problems and not args.force:
            print("  (not offered to --write-setups; --force after checking with project_survey.py)")
        else:
            new_setups[(cam, eo_path.name)] = out.name

    if args.write_setups:
        sc = CAL / "chelsea_setups.csv"
        if not sc.exists():
            sys.exit(f"{sc} not found; add the rows by hand")
        rows = list(csv.DictReader(open(sc, newline="")))
        changed = 0
        for r in rows:
            k = (r["camera"], r["eo_file"])
            if k in new_setups:
                r["eo_file"] = new_setups[k]; changed += 1
        with open(sc, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["camera", "first_date", "last_date", "eo_file"])
            w.writeheader(); w.writerows(rows)
        print(f"\nchelsea_setups.csv: {changed} row(s) now use the lidar-fitted pointing. Rerun "
              f"process_chelsea.py, then compare_dem_survey.py to confirm.")
        if any(k[1].endswith("20250219_EO.yaml") for k in new_setups):
            print("  NOTE: frames using the default 2025-02-19 pointing have no row to change; add "
                  "a row for their dates by hand if the fit holds.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

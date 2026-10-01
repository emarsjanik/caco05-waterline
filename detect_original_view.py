#!/usr/bin/env python3
"""
Waterlines Detected In The Original Photos (Historical Frames)
=================================================================
For a camera whose historical pointing is too far from today's for the
frames to be redrawn into today's view, detects the waterlines in the
ORIGINAL photos and georectifies them with that period's own pointing.

WHY. process_chelsea.py redraws every historical frame into the
2025-11-13 view, so the detector's hand-tuned crop, search envelope and
bias corrections apply. Once the Jan 2025 pointing was solved against the
23 Jan 2025 lidar (fit_eo_to_survey.py), c2 turned out to have looked
~25 deg further east than today: most of a redrawn frame is filled margin
and the detector found no waterline on any Jan 18-23 c2 frame. The
pointing is now known, so the redraw is not needed.

HOW.
  1. Search envelope from geometry, not tuning: every survey cell between
     --band elevations (default -1.5 to +2.5 m NAVD88, the intertidal zone
     with a margin) is projected into the photo with the period's pointing;
     per column, the rows it spans (padded) become the envelope. The band
     is broad -- about +/-20 m across the beach -- so the survey bounds the
     search without steering the line onto any particular contour.
  2. waterline_detector_v5.py runs on the original photos with that
     envelope, the whole frame as crop, and the view-specific bias
     correction switched off (--profile-json).
  3. extract_elevation_contours.py adds the ADCP water level and waves,
     georectify.py maps the points with the period's pointing.
  4. The points are merged with the main contour file (this camera's
     rows for the dates replaced), ready for dem_from_contours.py and
     compare_dem_survey.py.

Usage:
    python3 detect_original_view.py --camera c2 \\
        --eo calibration/CACO05_c2_2025-01-18_to_2025-01-23_lidar_EO.yaml \\
        --survey /mnt/I2Rgus_Data/Chelsea_calibration/2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif \\
        --start-date 2025-01-18 --end-date 2025-01-23
"""

import os
import re
import sys
import csv
import json
import shutil
import argparse
import subprocess
from pathlib import Path
from datetime import datetime, timezone

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from georectify import load_extrinsics, load_intrinsics      # noqa: E402
from view_reproject import ground_to_pixel                   # noqa: E402
from compare_dem_survey import read_survey                   # noqa: E402

CAL = HERE / "calibration"
CHELSEA = Path("/mnt/I2Rgus_Data/Chelsea_calibration")


def envelope_from_survey(io, eo, survey, zlo, zhi, max_range, n_cols=40,
                         pad_sea=0.06, pad_land=0.03):
    """
    Per-column [min_row, max_row] (fractions of the full image) spanned by
    the survey's zlo..zhi band, padded: more on the seaward (upper) side,
    where the lidar has no data below the water at flight time.
    Also returns which image side is nearer the camera.
    """
    grid, x0, y0, cell = survey
    cx, cy = eo[0], eo[1]
    c0 = max(int((cx - max_range - x0) / cell), 0)
    c1 = min(int((cx + max_range - x0) / cell), grid.shape[1])
    r0 = max(int((y0 - (cy + max_range)) / cell), 0)
    r1 = min(int((y0 - (cy - max_range)) / cell), grid.shape[0])
    sub = grid[r0:r1:2, c0:c1:2]                    # every other cell: 0.5 m is plenty
    rr, cc = np.nonzero(np.isfinite(sub) & (sub >= zlo) & (sub <= zhi))
    E = x0 + (c0 + 2 * cc + 0.5) * cell
    N = y0 - (r0 + 2 * rr + 0.5) * cell
    Z = sub[rr, cc]
    U, V, ok = ground_to_pixel(E, N, Z, io, eo)
    if ok.sum() < 100:
        sys.exit("the survey band hardly appears in this camera's view with this pointing")
    U, V = U[ok], V[ok]
    rng = np.hypot(E[ok] - cx, N[ok] - cy)
    W, H = io[0], io[1]
    edges = np.linspace(0, W, n_cols + 1)
    xs, lo, hi = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (U >= a) & (U < b)
        if m.sum() < 5:
            continue
        xs.append(((a + b) / 2) / W)
        lo.append(np.percentile(V[m], 1) / H - pad_sea)
        hi.append(np.percentile(V[m], 99) / H + pad_land)
    if len(xs) < 3:
        sys.exit("too few image columns see the survey band")
    # extend to the frame edges so the detector's interpolation covers them
    xs = [0.0] + xs + [1.0]
    lo = [lo[0]] + lo + [lo[-1]]
    hi = [hi[0]] + hi + [hi[-1]]
    env = [[round(x, 4), round(float(np.clip(a, 0, 1)), 4), round(float(np.clip(b, 0, 1)), 4)]
           for x, a, b in zip(xs, lo, hi)]
    left_far = np.median(rng[U < W / 3]) > np.median(rng[U > 2 * W / 3]) \
        if (U < W / 3).any() and (U > 2 * W / 3).any() else False
    return env, left_far, float(ok.mean())


def run(cmd, log):
    with open(log, "a") as f:
        f.write("\n$ " + " ".join(str(c) for c in cmd) + "\n")
        p = subprocess.run([str(c) for c in cmd], stdout=f, stderr=subprocess.STDOUT)
    return p.returncode == 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--camera", required=True, choices=["c1", "c2"])
    ap.add_argument("--eo", required=True, help="the period's pointing (e.g. the lidar-fitted EO)")
    ap.add_argument("--survey", required=True, help="lidar DSM for the search envelope")
    ap.add_argument("--start-date", required=True); ap.add_argument("--end-date", required=True)
    ap.add_argument("--originals", default=str(CHELSEA / "work" / "original"),
                    help="folder of original timex photos")
    ap.add_argument("--water-level", default=str(CHELSEA / "adcp_water_level_navd88.csv"))
    ap.add_argument("--waves", default=str(CHELSEA / "adcp_waves.csv"))
    ap.add_argument("--merge-into", default=str(CHELSEA / "contour_points_ground.csv"),
                    help="main contour file; this camera's rows for the dates are replaced")
    ap.add_argument("--merged-output", default=str(CHELSEA / "contour_points_ground_merged.csv"))
    ap.add_argument("--band", nargs=2, type=float, default=[-1.5, 2.5],
                    help="elevations (m NAVD88) bounding the search envelope (default -1.5 2.5)")
    ap.add_argument("--range", type=float, default=300.0, help="survey within this many m")
    ap.add_argument("--min-signal-fraction", type=float, default=0.15)
    ap.add_argument("--utc-hours", default="13.5-18",
                    help="daylight window used for the frames (default 13.5-18, as process_chelsea)")
    ap.add_argument("--output", default=None,
                    help="work folder (default <Chelsea>/original_view_<camera>)")
    args = ap.parse_args()

    out = Path(args.output or CHELSEA / f"original_view_{args.camera}")
    for sub in ("src", "in", "detections", "debug", "overlays"):
        shutil.rmtree(out / sub, ignore_errors=True)
        (out / sub).mkdir(parents=True)
    log = out / "processing.log"
    log.write_text("")

    io = load_intrinsics(CAL / f"CACO05_{args.camera}_20240801_IO.yaml")
    eo_path = Path(args.eo) if Path(args.eo).exists() else CAL / args.eo
    eo = load_extrinsics(eo_path)
    survey = read_survey(args.survey)

    # 1. envelope
    env, left_far, frac_ok = envelope_from_survey(io, eo, survey, *args.band, args.range)
    prof = {"crop": [0.0, 1.0, 0.0, 1.0], "envelope": env, "no_bias": True}
    if left_far:      # profiles assume column 0 = near the camera; swap the per-column ramps
        prof["step"] = [4, 18]
        prof["profile_scale"] = [0.4, 1.0]
    (out / "profile.json").write_text(json.dumps({args.camera: prof}, indent=1))
    print(f"envelope          : {len(env)} control points from the survey's "
          f"{args.band[0]:+g}..{args.band[1]:+g} m band ({100 * frac_ok:.0f}% of its cells in view); "
          f"far side {'left' if left_far else 'right'}")

    # 2. frames
    h0, h1 = (float(v) for v in args.utc_hours.split("-"))
    n = 0
    for p in sorted(Path(args.originals).glob(f"*.{args.camera}.timex.jpg")):
        m = re.match(r"^(\d{9,11})\.", p.name)
        if not m:
            continue
        t = datetime.fromtimestamp(int(m.group(1)), tz=timezone.utc)
        if not (args.start_date <= t.strftime("%Y-%m-%d") <= args.end_date):
            continue
        if not (h0 <= t.hour + t.minute / 60 <= h1):
            continue
        dst = out / "src" / p.name
        try:
            os.link(p, dst)
        except OSError:
            shutil.copy2(p, dst)
        n += 1
    print(f"frames            : {n} {args.camera} original photo(s), {args.start_date} to "
          f"{args.end_date}, {args.utc_hours} UTC")
    if not n:
        sys.exit(f"no {args.camera} timex frames in {args.originals} for those dates")

    # 3. detect
    print("detecting         : (roughly 3-4 s per frame)")
    if not run([sys.executable, HERE / "waterline_detector_v5.py",
                "--profile-json", out / "profile.json",
                "--min-signal-fraction", args.min_signal_fraction,
                "--image-suffix", "timex.jpg",
                "--source-dir", out / "src", "--input-dir", out / "in",
                "--output-dir", out / "detections", "--debug-dir", out / "debug"], log):
        sys.exit(f"detector failed; see {log}")
    dets = sorted((out / "detections").glob("*.csv"))
    print(f"                    {len(dets)} frame(s) passed the detector's quality filters")
    if not dets:
        sys.exit(f"no detections; see {log} and the debug images in {out / 'debug'}")

    # 4. water level + georectify
    contours, ground = out / "contour_points.csv", out / "contour_points_ground.csv"
    if not run([sys.executable, HERE / "extract_elevation_contours.py", args.water_level,
                "--time-col", "time", "--level-col", "water_level_navd88",
                "--min-coverage", args.min_signal_fraction,
                "--waves", args.waves, "--wave-max-gap-minutes", 60,
                "--processed-dir", out / "detections", "--output", contours], log) \
            or not contours.exists():
        sys.exit(f"contour extraction failed; see {log}")
    if not run([sys.executable, HERE / "georectify.py", contours, ground,
                f"--io-{args.camera}", CAL / f"CACO05_{args.camera}_20240801_IO.yaml",
                f"--eo-{args.camera}", eo_path], log):
        sys.exit(f"georectification failed; see {log}")
    new = list(csv.DictReader(open(ground, newline="")))
    frames = sorted({r["source_file"] for r in new})
    print(f"waterlines        : {len(new)} point(s) in {len(frames)} frame(s) after the contour filters")

    # 5. overlays: the detected line on each original photo
    try:
        import cv2
        by = {}
        for r in new:
            by.setdefault(r["source_file"], []).append((float(r["pixel_column"]), float(r["pixel_row"]),
                                                        float(r["tide_elevation_navd88"])))
        for f, pts in by.items():
            img = cv2.imread(str(next(iter((out / "src").glob(f + "*")), "")))
            if img is None:
                continue
            pts = sorted(pts)
            poly = np.array([[u, v] for u, v, _ in pts], np.int32).reshape(-1, 1, 2)
            cv2.polylines(img, [poly], False, (0, 0, 0), 9, cv2.LINE_AA)
            cv2.polylines(img, [poly], False, (0, 255, 255), 5, cv2.LINE_AA)
            # envelope bounds, thin
            W, H = img.shape[1], img.shape[0]
            for k in (1, 2):
                line = np.array([[e[0] * W, e[k] * H] for e in env], np.int32).reshape(-1, 1, 2)
                cv2.polylines(img, [line], False, (255, 0, 255), 2, cv2.LINE_AA)
            cv2.putText(img, f"{pts[0][2]:+.2f} m NAVD88", (40, 90), cv2.FONT_HERSHEY_SIMPLEX, 2.4,
                        (0, 0, 0), 10, cv2.LINE_AA)
            cv2.putText(img, f"{pts[0][2]:+.2f} m NAVD88", (40, 90), cv2.FONT_HERSHEY_SIMPLEX, 2.4,
                        (0, 255, 255), 4, cv2.LINE_AA)
            cv2.imwrite(str(out / "overlays" / (f + ".jpg")),
                        cv2.resize(img, (W // 2, H // 2)), [cv2.IMWRITE_JPEG_QUALITY, 85])
        print(f"overlays          : {out / 'overlays'} (yellow = waterline, magenta = search envelope)")
    except Exception as e:                       # overlays are a convenience
        print(f"overlays skipped: {e}")

    # 6. merge into the main contour file
    if args.merge_into and Path(args.merge_into).exists() and new:
        base = list(csv.DictReader(open(args.merge_into, newline="")))
        keep = [r for r in base if not (r["camera"] == args.camera and
                                        args.start_date <= r.get("capture_time_utc", "")[:10] <= args.end_date)]
        fields = list(base[0].keys()) if base else []
        for k in new[0].keys():
            if k not in fields:
                fields.append(k)
        with open(args.merged_output, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(keep)
            w.writerows(new)
        print(f"merged            : {args.merged_output} ({len(base) - len(keep)} old {args.camera} "
              f"row(s) for those dates replaced by {len(new)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())

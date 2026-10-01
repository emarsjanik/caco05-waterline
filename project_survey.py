#!/usr/bin/env python3
"""
Draw A Survey's Contours Onto A Camera Photo
===============================================
Projects elevation contours of a lidar DSM (or any survey grid) into a
camera image with one or more candidate pointings (EO files), each in its
own colour. If a pointing is right, the contours land on what they
describe: the dune toe and berm edge on the dune toe and berm, the 0 m
line near the water's edge at mid tide. A wrong pointing misses them by
an amount anyone can see.

WHY. Waterline comparisons against the 23 Jan 2025 lidar put the Jan
2025 lines 15-30 m off, and a pointing fit to the lidar wanted the same
~22 deg pan change for both cameras -- yet a Jan 21 / Feb 19 blend shows
no such turn. This checks the pointing on the photo itself, with no
waterline detection or water level involved.

Usage:
    python3 project_survey.py --image <original photo> --camera c2 \\
        --survey 2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif \\
        --eo calibration/CACO05_c2_2025-01-08_to_2025-01-23_EO.yaml \\
             calibration/CACO05_c2_20250219_EO.yaml \\
        --levels 0 1 2 3 --output ~/projected_c2.jpg
"""

import sys
import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from georectify import load_extrinsics, load_intrinsics          # noqa: E402
from view_reproject import ground_to_pixel                       # noqa: E402
from compare_dem_survey import read_survey                       # noqa: E402

CAL = HERE / "calibration"
# BGR, one per EO file, readable on sand and water
COLOURS = [(0, 255, 255), (255, 255, 0), (255, 0, 255), (0, 165, 255), (0, 255, 0)]


def contour_segments(grid, x0, y0, cell, levels, bbox):
    """World-coordinate contour polylines {level: [array(n,2), ...]}, via matplotlib."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xa, xb, ya, yb = bbox
    c0, c1 = max(int((xa - x0) / cell), 0), min(int((xb - x0) / cell), grid.shape[1])
    r0, r1 = max(int((y0 - yb) / cell), 0), min(int((y0 - ya) / cell), grid.shape[0])
    sub = grid[r0:r1, c0:c1]
    xs = x0 + (np.arange(c0, c1) + 0.5) * cell
    ys = y0 - (np.arange(r0, r1) + 0.5) * cell
    cs = plt.contour(xs, ys, np.ma.masked_invalid(sub), levels=sorted(levels))
    out = {lv: [s for s in segs if len(s) > 5] for lv, segs in zip(cs.levels, cs.allsegs)}
    plt.close("all")
    return out


def main():
    import cv2
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--image", required=True, help="original (not resampled) photo")
    ap.add_argument("--camera", required=True, choices=["c1", "c2"])
    ap.add_argument("--survey", required=True)
    ap.add_argument("--eo", nargs="+", required=True, help="EO file(s) to compare, one colour each")
    ap.add_argument("--levels", nargs="+", type=float, default=[0.0, 1.0, 2.0, 3.0])
    ap.add_argument("--range", type=float, default=250.0,
                    help="Only survey within this distance of the camera, m (default 250)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    img = cv2.imread(args.image)
    if img is None:
        sys.exit(f"cannot read {args.image}")
    io = load_intrinsics(CAL / f"CACO05_{args.camera}_20240801_IO.yaml")
    eos = [(Path(p).name, load_extrinsics(p if Path(p).exists() else CAL / p)) for p in args.eo]
    grid, x0, y0, cell = read_survey(args.survey)
    cx, cy = eos[0][1][0], eos[0][1][1]
    segs = contour_segments(grid, x0, y0, cell, args.levels,
                            (cx - args.range, cx + args.range, cy - args.range, cy + args.range))
    print(f"contours: " + ", ".join(f"{lv:+g} m ({len(s)} pieces)" for lv, s in segs.items()))
    for k, (name, eo) in enumerate(eos):
        col = COLOURS[k % len(COLOURS)]
        drawn = 0
        for lv, pieces in segs.items():
            for s in pieces:
                U, V, ok = ground_to_pixel(s[:, 0], s[:, 1], np.full(len(s), lv), io, eo)
                # draw only runs of consecutive in-frame points
                idx = np.where(ok)[0]
                if len(idx) < 2:
                    continue
                runs = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
                for r in runs:
                    if len(r) < 2:
                        continue
                    pts = np.c_[U[r], V[r]].astype(np.int32).reshape(-1, 1, 2)
                    cv2.polylines(img, [pts], False, (0, 0, 0), 7, cv2.LINE_AA)
                    cv2.polylines(img, [pts], False, col, 3, cv2.LINE_AA)
                    mid = r[len(r) // 2]
                    cv2.putText(img, f"{lv:+g}", (int(U[mid]) + 6, int(V[mid]) - 6),
                                cv2.FONT_HERSHEY_SIMPLEX, 1.2, col, 3, cv2.LINE_AA)
                    drawn += 1
        cv2.putText(img, name, (30, 60 + 55 * k), cv2.FONT_HERSHEY_SIMPLEX, 1.4, (0, 0, 0), 8,
                    cv2.LINE_AA)
        cv2.putText(img, name, (30, 60 + 55 * k), cv2.FONT_HERSHEY_SIMPLEX, 1.4, col, 3, cv2.LINE_AA)
        print(f"  {name}: {drawn} contour piece(s) in the frame")
    small = cv2.resize(img, (img.shape[1] // 2, img.shape[0] // 2), interpolation=cv2.INTER_AREA)
    cv2.imwrite(args.output, small, [cv2.IMWRITE_JPEG_QUALITY, 85])
    print(f"wrote {args.output} (half size)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Sweep profile_scale_far Against Ground Truth
-----------------------------------------------
Re-detects ONE frame at several values of profile_scale_far and reports
the error against a hand-traced ground truth, broken down by column
band.

WHY: on a well-lit C2 midday frame the detector's line goes flat from
about column 1900 to the right edge, while the water/sand boundary is
plainly visible there and the ground truth follows it. The far-field
error reached 216 px. The profile matcher samples vertical zones
(water -40:-15, foam -15:+5, wet sand +5:+25, dry sand +25:+55 px)
scaled by profile_scale_far at the far end; at 0.4 the water zone
shrinks to roughly 16-6 px. The hypothesis under test is that those
windows collapse faster than the feature does in the foreshortened far
field, so the matcher stops responding to a transition that is still
there.

This tests that hypothesis by measurement rather than by eye. A value
that reduces far-band error without inflating near-band error supports
it; no improvement at any value refutes it, and the cause lies
elsewhere.

Run from the directory holding waterline_detector_v5.py.

Usage:
    python3 sweep_profile_scale.py <image> <ground_truth_csv> <camera c1|c2>
        [--values 0.4,0.5,0.6,0.7,0.8,1.0]
"""

import sys
import csv
import shutil
import tempfile
import argparse
from pathlib import Path
from dataclasses import replace

import numpy as np

import waterline_detector_v5 as detector


CAMERA_KEYS = {"c1": "CACO05_C1", "c2": "CACO05_C2"}
BANDS = [("near 0-800", 0, 800), ("mid 800-1600", 800, 1600),
         ("far 1600-2448", 1600, 2448)]


def load_two_column_csv(path):
    cols, rows = [], []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        precise = reader.fieldnames and "Row_Precise" in reader.fieldnames
        for r in reader:
            cols.append(int(r["Column"]))
            if precise and r.get("Row_Precise") not in ("", None):
                rows.append(float(r["Row_Precise"]))
            else:
                rows.append(float(r["Row"]))
    return np.array(cols), np.array(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("image")
    ap.add_argument("ground_truth_csv")
    ap.add_argument("camera", choices=["c1", "c2"])
    ap.add_argument("--param", default="profile_scale_far",
                    help="Which parameter to sweep: profile_scale_far, step_far, step_near "
                         "(CameraProfile fields), or vertical_penalty (a CONFIG field).")
    ap.add_argument("--values", default="0.4,0.5,0.6,0.7,0.8,1.0")
    args = ap.parse_args()

    values = [float(v) for v in args.values.split(",")]
    key = CAMERA_KEYS[args.camera]
    original_profile = detector.CAMERAS[key]
    original_penalty = detector.CONFIG.vertical_penalty

    gt_cols, gt_rows = load_two_column_csv(args.ground_truth_csv)

    print("=" * 84)
    print(f"{args.param} sweep -- {Path(args.image).name}")
    print(f"ground truth: {len(gt_cols)} columns from {Path(args.ground_truth_csv).name}")
    print("=" * 84)
    print(f"{'scale':>6} | " + " | ".join(f"{b[0]:>14}" for b in BANDS) + " |    overall")
    print("-" * 84)

    results = []
    for value in values:
        workdir = Path(tempfile.mkdtemp())
        try:
            src = workdir / "src"; src.mkdir()
            shutil.copy2(args.image, src / Path(args.image).name)

            if args.param == "vertical_penalty":
                detector.CONFIG.vertical_penalty = value
            else:
                cast = int if args.param.startswith("step_") else float
                detector.CAMERAS[key] = replace(original_profile, **{args.param: cast(value)})
            detector.SOURCE_FOLDER = str(src)
            detector.INPUT_FOLDER = str(workdir / "in")
            detector.OUTPUT_FOLDER = str(workdir / "out")
            detector.DEBUG_FOLDER = str(workdir / "dbg")
            for d in ("in", "out", "dbg"):
                (workdir / d).mkdir(exist_ok=True)

            # Silence per-frame chatter; we only want the numbers.
            real_print = detector.print if hasattr(detector, "print") else None
            import io, contextlib
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                shutil.copy2(args.image, workdir / "in" / Path(args.image).name)
                result, status = detector.process_image(workdir / "in" / Path(args.image).name)

            out_csv = workdir / "out" / (Path(args.image).stem + ".csv")
            if not out_csv.exists():
                print(f"{value:6.2f} | frame not saved (status: {status})")
                continue

            det_cols, det_rows = load_two_column_csv(out_csv)
            gt_interp = np.interp(det_cols, gt_cols, gt_rows)
            error = det_rows - gt_interp

            cells = []
            for _, lo, hi in BANDS:
                m = (det_cols >= lo) & (det_cols < hi)
                cells.append(f"{np.sqrt((error[m]**2).mean()):14.1f}" if m.any()
                             else f"{'--':>14}")
            overall = np.sqrt((error ** 2).mean())
            print(f"{value:6.2f} | " + " | ".join(cells) + f" | {overall:10.1f}")
            results.append((value, overall, error))
        finally:
            shutil.rmtree(workdir, ignore_errors=True)

    detector.CAMERAS[key] = original_profile
    detector.CONFIG.vertical_penalty = original_penalty

    print("-" * 84)
    print("values are RMS error in pixels against ground truth (lower is better)")
    print()
    if results:
        best = min(results, key=lambda r: r[1])
        defaults = {"profile_scale_far": 0.4, "step_far": 4, "step_near": 18,
                    "vertical_penalty": 0.35}
        dv = defaults.get(args.param)
        baseline = next((r for r in results if dv is not None and abs(r[0] - dv) < 1e-9), None)
        print(f"Best: profile_scale_far={best[0]:.2f} at {best[1]:.1f} px overall")
        if baseline and best[0] != baseline[0]:
            change = 100 * (best[1] - baseline[1]) / baseline[1]
            print(f"  vs current {baseline[0]:g} ({baseline[1]:.1f} px): {change:+.0f}%")
            print()
            print("  ONE FRAME ONLY. A value that wins here may lose on other tide states,")
            print("  lighting, or on c1 -- whose geometry differs. Before changing the")
            print("  detector, re-run this on several frames spanning the tidal range, and")
            print("  check the near band did not get worse to buy the far-band gain.")
        elif baseline:
            print(f"  The current value is already best among those tested, so {args.param}")
            print("  is not the cause and it lies elsewhere.")


if __name__ == "__main__":
    main()

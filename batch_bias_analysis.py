#!/usr/bin/env python3
"""
Batch Bias Analysis Tool
--------------------------
Automatically finds every ground truth CSV, matches it to the
corresponding detector output CSV (same base filename) in the
processed/ folder, and runs RMS error + signed bias analysis on each
pair. Prints a per-image summary and an overall average at the end.

Usage:
    python3 batch_bias_analysis.py
"""

import csv
import numpy as np
from pathlib import Path


PROCESSED_FOLDER = "/mnt/I2Rgus_Data/waterline/processed"
GROUND_TRUTH_FOLDER = "/mnt/I2Rgus_Data/waterline/ground_truth"


def load_csv_rows(path, row_col_index):
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        next(reader)
        for line in reader:
            rows.append(int(line[row_col_index]))
    return np.array(rows, dtype=np.float32)


N_SEGMENTS = 5


def detect_camera(name):
    name = name.lower()
    if ".c1." in name:
        return "C1"
    if ".c2." in name:
        return "C2"
    return "UNKNOWN"


def analyze_pair(detector_csv, ground_truth_csv):
    det_rows = load_csv_rows(detector_csv, row_col_index=1)
    gt_rows = load_csv_rows(ground_truth_csv, row_col_index=1)

    n = min(len(det_rows), len(gt_rows))
    det_rows = det_rows[:n]
    gt_rows = gt_rows[:n]

    signed_error = det_rows - gt_rows
    abs_error = np.abs(signed_error)

    rms = float(np.sqrt(np.mean(signed_error ** 2)))
    mae = float(np.mean(abs_error))
    max_dev = float(np.max(abs_error))
    mean_signed = float(np.mean(signed_error))
    within_10 = 100.0 * np.sum(abs_error <= 10) / n
    within_25 = 100.0 * np.sum(abs_error <= 25) / n

    # Per-segment mean signed error, same 5-way column split as
    # bias_analysis.py, so results from both tools line up directly.
    # This is what lets us tell a genuine per-column (perspective)
    # effect apart from single-frame noise: a real effect should show
    # the same segment shape averaged across many frames, not just on
    # one.
    segment_size = n // N_SEGMENTS
    segment_means = []
    for i in range(N_SEGMENTS):
        lo = i * segment_size
        hi = n if i == N_SEGMENTS - 1 else (i + 1) * segment_size
        segment_means.append(float(np.mean(signed_error[lo:hi])))

    return {
        "rms": rms,
        "mae": mae,
        "max_dev": max_dev,
        "mean_signed": mean_signed,
        "within_10": within_10,
        "within_25": within_25,
        "n": n,
        "segment_means": segment_means,
    }


def main():
    gt_dir = Path(GROUND_TRUTH_FOLDER)
    processed_dir = Path(PROCESSED_FOLDER)

    gt_files = sorted(gt_dir.glob("*.csv"))

    if not gt_files:
        print(f"No ground truth CSVs found in {GROUND_TRUTH_FOLDER}")
        return

    print(f"Found {len(gt_files)} ground truth file(s).")
    print()

    results = []
    missing = []

    for gt_path in gt_files:
        base_name = gt_path.stem  # e.g. "1763666100...c2.snap"
        detector_path = processed_dir / (base_name + ".csv")

        if not detector_path.exists():
            missing.append(base_name)
            continue

        stats = analyze_pair(detector_path, gt_path)
        stats["name"] = base_name
        stats["camera"] = detect_camera(base_name)
        results.append(stats)

        print(f"{base_name}")
        print(f"  RMS: {stats['rms']:7.2f}px   MAE: {stats['mae']:7.2f}px   "
              f"Max: {stats['max_dev']:7.2f}px   Mean signed: {stats['mean_signed']:+7.2f}px   "
              f"Within10px: {stats['within_10']:5.1f}%   Within25px: {stats['within_25']:5.1f}%")

    if missing:
        print()
        print(f"WARNING: {len(missing)} ground truth file(s) had no matching "
              f"detector output in {PROCESSED_FOLDER} (skipped):")
        for name in missing:
            print(f"  - {name}")

    if results:
        print()
        print("=" * 100)
        print("OVERALL SUMMARY (average across all matched images)")
        print("=" * 100)

        avg_rms = np.mean([r["rms"] for r in results])
        avg_mae = np.mean([r["mae"] for r in results])
        avg_max = np.mean([r["max_dev"] for r in results])
        avg_signed = np.mean([r["mean_signed"] for r in results])
        avg_within_10 = np.mean([r["within_10"] for r in results])
        avg_within_25 = np.mean([r["within_25"] for r in results])

        print(f"Images compared     : {len(results)}")
        print(f"Average RMS         : {avg_rms:.2f} px")
        print(f"Average MAE         : {avg_mae:.2f} px")
        print(f"Average Max Dev     : {avg_max:.2f} px")
        print(f"Average Mean Signed : {avg_signed:+.2f} px")
        print(f"Average Within 10px : {avg_within_10:.1f}%")
        print(f"Average Within 25px : {avg_within_25:.1f}%")
        print("=" * 100)

        # Sort and show worst offenders, since averages can hide a few
        # very bad frames -- worth knowing which specific images are
        # dragging down the average.
        worst = sorted(results, key=lambda r: -r["rms"])[:5]
        print()
        print("Worst 5 by RMS error:")
        for r in worst:
            print(f"  {r['name']}: RMS {r['rms']:.2f}px")

        # Per-camera, per-segment average signed error, aggregated
        # across every matched frame for that camera. This is the
        # tool for telling a genuine perspective/column-position
        # effect apart from single-frame noise: a real effect shows
        # the same segment shape averaged over many frames; noise
        # averages toward flat/zero across frames even if any single
        # frame looks dramatic.
        cameras = sorted(set(r["camera"] for r in results))
        print()
        print("=" * 100)
        print("PER-CAMERA SEGMENT-WISE MEAN SIGNED ERROR (averaged across all matched frames)")
        print("=" * 100)
        for camera in cameras:
            camera_results = [r for r in results if r["camera"] == camera]
            n_frames = len(camera_results)
            segment_avgs = np.mean(
                [r["segment_means"] for r in camera_results], axis=0
            )
            segment_str = "  ".join(f"{v:+8.2f}px" for v in segment_avgs)
            print(f"{camera} (n={n_frames} frames): {segment_str}")
        print("=" * 100)
        print("Columns are left-to-right segments (near camera -> far camera).")
        print("A genuine per-column/perspective effect should show a consistent")
        print("shape here across frames; if it looks flat/noisy instead, the")
        print("per-frame segment numbers you were looking at were likely just")
        print("single-frame variation, not a structural effect worth correcting.")


if __name__ == "__main__":
    main()

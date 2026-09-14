#!/usr/bin/env python3
"""
Ground Truth Comparison Tool
-----------------------------
Compares the detector's output CSV against a hand-collected ground
truth CSV, and reports RMS error, mean absolute error, max deviation,
and percentage of columns within +/-5 pixels.

Usage:
    python3 compare_to_ground_truth.py <detector_csv> <ground_truth_csv>

Example:
    python3 compare_to_ground_truth.py \
        /mnt/I2Rgus_Data/waterline/processed/1763666100.Thu.Nov.20_19_15_00.GMT.2025.CACO05.c1.snap.csv \
        /mnt/I2Rgus_Data/waterline/ground_truth/c1_19_15_00.csv
"""

import sys
import csv
import numpy as np
from pathlib import Path


def load_csv_rows(path, has_confidence):
    columns = []
    rows = []

    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)  # skip header row

        for line in reader:
            columns.append(int(line[0]))
            rows.append(int(line[1]))

    return np.array(columns), np.array(rows)


def compute_rms_error(reference, detected):
    reference = np.asarray(reference, dtype=np.float32)
    detected = np.asarray(detected, dtype=np.float32)

    n = min(len(reference), len(detected))
    if n == 0:
        return np.nan

    error = reference[:n] - detected[:n]
    return float(np.sqrt(np.mean(error * error)))


def compute_mae(reference, detected):
    reference = np.asarray(reference, dtype=np.float32)
    detected = np.asarray(detected, dtype=np.float32)

    n = min(len(reference), len(detected))
    error = np.abs(reference[:n] - detected[:n])
    return float(np.mean(error))


def compute_max_deviation(reference, detected):
    reference = np.asarray(reference, dtype=np.float32)
    detected = np.asarray(detected, dtype=np.float32)

    n = min(len(reference), len(detected))
    error = np.abs(reference[:n] - detected[:n])
    return float(np.max(error))


def compute_within_tolerance(reference, detected, tolerance=5):
    reference = np.asarray(reference, dtype=np.float32)
    detected = np.asarray(detected, dtype=np.float32)

    n = min(len(reference), len(detected))
    error = np.abs(reference[:n] - detected[:n])
    within = np.sum(error <= tolerance)
    return 100.0 * within / n


def rating_from_rms(rms):
    if rms < 2:
        return "Excellent"
    elif rms < 5:
        return "Very Good"
    elif rms < 10:
        return "Good"
    else:
        return "Needs Improvement"


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 compare_to_ground_truth.py <detector_csv> <ground_truth_csv>")
        sys.exit(1)

    detector_csv = Path(sys.argv[1])
    ground_truth_csv = Path(sys.argv[2])

    if not detector_csv.exists():
        print(f"ERROR: detector CSV not found: {detector_csv}")
        sys.exit(1)

    if not ground_truth_csv.exists():
        print(f"ERROR: ground truth CSV not found: {ground_truth_csv}")
        sys.exit(1)

    det_columns, det_rows = load_csv_rows(detector_csv, has_confidence=True)
    gt_columns, gt_rows = load_csv_rows(ground_truth_csv, has_confidence=False)

    if len(det_columns) != len(gt_columns):
        print(
            f"WARNING: column count mismatch "
            f"(detector={len(det_columns)}, ground_truth={len(gt_columns)}). "
            f"Comparing over the overlapping range only."
        )

    rms = compute_rms_error(gt_rows, det_rows)
    mae = compute_mae(gt_rows, det_rows)
    max_dev = compute_max_deviation(gt_rows, det_rows)
    within_5 = compute_within_tolerance(gt_rows, det_rows, tolerance=5)
    within_10 = compute_within_tolerance(gt_rows, det_rows, tolerance=10)

    print()
    print("=" * 60)
    print("Ground Truth Comparison")
    print("=" * 60)
    print(f"Detector CSV     : {detector_csv.name}")
    print(f"Ground truth CSV : {ground_truth_csv.name}")
    print()
    print(f"RMS Error              : {rms:.2f} pixels")
    print(f"Mean Absolute Error     : {mae:.2f} pixels")
    print(f"Max Deviation           : {max_dev:.2f} pixels")
    print(f"Within +/-5 pixels      : {within_5:.1f}%")
    print(f"Within +/-10 pixels     : {within_10:.1f}%")
    print()
    print(f"Quality Rating          : {rating_from_rms(rms)}")
    print("=" * 60)
    print()


if __name__ == "__main__":
    main()

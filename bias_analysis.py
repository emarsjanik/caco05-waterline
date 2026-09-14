#!/usr/bin/env python3
"""
Bias Analysis Tool
-------------------
For a detector/ground-truth CSV pair, reports not just RMS error but
the SIGNED mean error (positive = detector is BELOW true line / too
far toward dry sand, negative = detector is ABOVE / too far toward
open water), plus a per-column breakdown to see if the offset is
truly constant or varies across the frame.

Usage:
    python3 bias_analysis.py <detector_csv> <ground_truth_csv>
"""

import sys
import csv
import numpy as np
from pathlib import Path


def load_csv_rows(path):
    columns, rows = [], []
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        next(reader)
        for line in reader:
            columns.append(int(line[0]))
            rows.append(int(line[1]))
    return np.array(columns), np.array(rows)


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 bias_analysis.py <detector_csv> <ground_truth_csv>")
        sys.exit(1)

    det_path = Path(sys.argv[1])
    gt_path = Path(sys.argv[2])

    _, det_rows = load_csv_rows(det_path)
    _, gt_rows = load_csv_rows(gt_path)

    n = min(len(det_rows), len(gt_rows))
    det_rows = det_rows[:n].astype(np.float32)
    gt_rows = gt_rows[:n].astype(np.float32)

    signed_error = det_rows - gt_rows

    print()
    print("=" * 60)
    print("Bias Analysis")
    print("=" * 60)
    print(f"Detector CSV     : {det_path.name}")
    print(f"Ground truth CSV : {gt_path.name}")
    print()
    print(f"Mean signed error   : {np.mean(signed_error):+.2f} px")
    print(f"Median signed error : {np.median(signed_error):+.2f} px")
    print(f"Std dev of error    : {np.std(signed_error):.2f} px")
    print()
    print("Interpretation:")
    print("  Positive mean -> detector consistently BELOW true line")
    print("                   (toward dry sand / berm)")
    print("  Negative mean -> detector consistently ABOVE true line")
    print("                   (toward open water)")
    print()

    segment_size = n // 5
    print("Signed error by horizontal segment (left to right):")
    for i in range(5):
        lo = i * segment_size
        hi = n if i == 4 else (i + 1) * segment_size
        seg_mean = np.mean(signed_error[lo:hi])
        print(f"  Segment {i+1} (cols {lo}-{hi}): {seg_mean:+.2f} px")

    print("=" * 60)
    print()


if __name__ == "__main__":
    main()

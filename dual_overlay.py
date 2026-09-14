#!/usr/bin/env python3
"""
Dual-Line Overlay
--------------------
Draws BOTH the detector's line (red) and the ground truth line (lime)
on the same image, so we can see directly how/where they diverge,
rather than inferring it from numbers alone.

Usage:
    python3 dual_overlay.py <image_path> <detector_csv> <ground_truth_csv> <output_path>
"""

import sys
import csv
import cv2
import numpy as np
from pathlib import Path


def load_csv_rows(path):
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        next(reader)
        for line in reader:
            rows.append(int(line[1]))
    return np.array(rows, dtype=np.int32)


def draw_line(image, rows, color, thickness=3):
    points = [(x, int(y)) for x, y in enumerate(rows)]
    for i in range(len(points) - 1):
        cv2.line(image, points[i], points[i + 1], color, thickness, cv2.LINE_AA)
    return image


def main():
    if len(sys.argv) != 5:
        print("Usage: python3 dual_overlay.py <image_path> <detector_csv> <ground_truth_csv> <output_path>")
        sys.exit(1)

    image_path = Path(sys.argv[1])
    det_path = Path(sys.argv[2])
    gt_path = Path(sys.argv[3])
    output_path = Path(sys.argv[4])

    image = cv2.imread(str(image_path))
    if image is None:
        print(f"Could not load image: {image_path}")
        sys.exit(1)

    det_rows = load_csv_rows(det_path)
    gt_rows = load_csv_rows(gt_path)

    n = min(len(det_rows), len(gt_rows), image.shape[1])
    det_rows = det_rows[:n]
    gt_rows = gt_rows[:n]

    # Ground truth (lime, drawn first so red is on top and easy to see)
    draw_line(image, gt_rows, (0, 255, 0), thickness=4)
    # Detector (red)
    draw_line(image, det_rows, (0, 0, 255), thickness=3)

    # Label markers at both ends so left/right is unambiguous
    cv2.putText(image, "LEFT (col 0)", (20, 60), cv2.FONT_HERSHEY_SIMPLEX, 1.2, (255, 255, 0), 3)
    cv2.putText(image, f"RIGHT (col {n-1})", (image.shape[1] - 500, 60),
                cv2.FONT_HERSHEY_SIMPLEX, 1.2, (255, 255, 0), 3)

    cv2.imwrite(str(output_path), image)
    print(f"Saved: {output_path}")
    print(f"Green = ground truth, Red = detector")
    print(f"Column 0 (LEFT side of image) to column {n-1} (RIGHT side of image)")


if __name__ == "__main__":
    main()

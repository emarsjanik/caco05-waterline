#!/usr/bin/env python3
"""
Crop Coverage Check
---------------------
Shows where the ground truth line falls relative to the camera's
current crop_top/crop_bottom, in both fractional and pixel terms --
to check whether the crop region is even well-positioned for this
camera's true shoreline.

Usage:
    python3 check_crop_coverage.py <ground_truth_csv> <image_path> <crop_top_frac> <crop_bottom_frac>
"""

import sys
import csv
import numpy as np
import cv2
from pathlib import Path


def main():
    if len(sys.argv) != 5:
        print("Usage: python3 check_crop_coverage.py <ground_truth_csv> <image_path> <crop_top_frac> <crop_bottom_frac>")
        sys.exit(1)

    gt_path = Path(sys.argv[1])
    image_path = Path(sys.argv[2])
    crop_top_frac = float(sys.argv[3])
    crop_bottom_frac = float(sys.argv[4])

    image = cv2.imread(str(image_path))
    height = image.shape[0]

    crop_top_px = int(crop_top_frac * height)
    crop_bottom_px = int(crop_bottom_frac * height)

    rows = []
    with open(gt_path, "r", newline="") as f:
        reader = csv.reader(f)
        next(reader)
        for line in reader:
            rows.append(int(line[1]))

    rows = np.array(rows, dtype=np.float32)

    print()
    print("=" * 60)
    print("Crop Coverage Check")
    print("=" * 60)
    print(f"Image height        : {height}px")
    print(f"Current crop region  : {crop_top_px}px to {crop_bottom_px}px "
          f"({crop_top_frac:.2f} to {crop_bottom_frac:.2f})")
    print()
    print(f"Ground truth min row : {rows.min():.0f}px  ({rows.min()/height:.3f} of image height)")
    print(f"Ground truth max row : {rows.max():.0f}px  ({rows.max()/height:.3f} of image height)")
    print(f"Ground truth mean row: {rows.mean():.0f}px  ({rows.mean()/height:.3f} of image height)")
    print()

    below_crop = np.sum(rows < crop_top_px)
    above_crop = np.sum(rows > crop_bottom_px)
    inside_crop = len(rows) - below_crop - above_crop

    print(f"Ground truth points ABOVE crop_top    : {below_crop} ({100*below_crop/len(rows):.1f}%)")
    print(f"Ground truth points INSIDE crop region: {inside_crop} ({100*inside_crop/len(rows):.1f}%)")
    print(f"Ground truth points BELOW crop_bottom : {above_crop} ({100*above_crop/len(rows):.1f}%)")
    print("=" * 60)
    print()


if __name__ == "__main__":
    main()

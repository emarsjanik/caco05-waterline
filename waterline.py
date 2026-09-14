#!/usr/bin/env python3
"""
Waterline Detection Script (fixed-camera, curved dune mask)
-------------------------------------------------------------
Detects the waterline by finding the strongest brightness gradient
in each image column, restricted to a hand-set curved region that
excludes the dune/vegetation area (since this is a fixed webcam,
that region never moves between frames).

Includes a debug mode that draws the current dune-mask boundary
as a BLUE line so you can visually tune DUNE_EDGE_POINTS.

Requirements:
    pip install opencv-python numpy

Usage:
    python3 waterline_detector.py
"""

import cv2
import numpy as np
import os
import glob

# ---------------- CONFIG ----------------
IMAGE_DIR = "/mnt/I2Rgus_Data/waterline"
INPUT_FILENAME = "beach_2026-07-21_1300.jpg"   # only used by main_single() and debug mode

TOP_CROP_FRAC = 0.20     # skip sky/timestamp
BOTTOM_CROP_FRAC = 0.10  # skip close dune toe / watermark

# --- DUNE EXCLUSION CURVE (piecewise, 3+ points) ---
# Each tuple is (height_frac_in_roi, x_frac_of_width): at that fractional
# height down the cropped ROI, everything LEFT of that x-fraction is
# excluded from the waterline search (treated as dune/vegetation).
# Tune these using DEBUG MODE below -- once set correctly for this fixed
# camera, they should not need to change again.
DUNE_EDGE_POINTS = [
    (0.00, 0.22),   # top of ROI
    (0.55, 0.52),   # mid-height -- covers the sandy dune-apron bulge
    (1.00, 0.65),   # bottom of ROI
]

BLUR_KERNEL = (21, 21)
SMOOTHING_WINDOW = 41
MEDIAN_WINDOW = 15

LINE_COLOR = (0, 0, 255)   # red in BGR
LINE_THICKNESS = 3


# ---------------- CORE FUNCTIONS ----------------

def build_dune_mask(height, width):
    """
    Returns a boolean mask, shape (height, width), True = excluded
    (dune/vegetation region). The boundary is a piecewise-linear curve
    through DUNE_EDGE_POINTS, so it can bulge to match the real dune
    shape instead of assuming a single straight diagonal.
    """
    mask = np.zeros((height, width), dtype=bool)

    heights = np.array([p[0] for p in DUNE_EDGE_POINTS])
    x_fracs = np.array([p[1] for p in DUNE_EDGE_POINTS])

    for row in range(height):
        frac = row / max(1, height - 1)
        boundary_x_frac = np.interp(frac, heights, x_fracs)
        boundary_x = int(boundary_x_frac * width)
        mask[row, :boundary_x] = True

    return mask


def find_waterline(L_roi, exclude_mask):
    """
    For each column, find the row with the strongest vertical brightness
    gradient, ignoring excluded (dune) rows. Columns that are entirely
    excluded return -1 (invalid, not drawn).
    """
    blur = cv2.GaussianBlur(L_roi, BLUR_KERNEL, 0).astype(np.float32)
    height, width = blur.shape

    gradient = np.abs(np.diff(blur, axis=0))
    exclude_grad = exclude_mask[:-1, :]

    gradient_masked = gradient.copy()
    gradient_masked[exclude_grad] = 0

    waterline_rows = np.argmax(gradient_masked, axis=0)
    col_has_signal = gradient_masked.max(axis=0) > 0
    waterline_rows = np.where(col_has_signal, waterline_rows, -1)

    return waterline_rows


def median_filter_1d(rows, window):
    """Rolling median, operating only on valid (>=0) points; -1 stays -1."""
    n = len(rows)
    half = window // 2
    out = rows.copy()
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        window_vals = rows[lo:hi]
        valid_vals = window_vals[window_vals >= 0]
        if len(valid_vals) > 0:
            out[i] = int(np.median(valid_vals))
    return out


def smooth_line(rows, window):
    """Moving average, but only across valid points, leaving -1 as -1."""
    n = len(rows)
    half = window // 2
    out = rows.copy().astype(np.float32)
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        window_vals = rows[lo:hi]
        valid_vals = window_vals[window_vals >= 0]
        if len(valid_vals) > 0:
            out[i] = np.mean(valid_vals)
        else:
            out[i] = -1
    return out.astype(np.int32)


def draw_waterline(image, rows, top_offset, color=LINE_COLOR, thickness=LINE_THICKNESS):
    """
    Draw the line only across contiguous runs of valid (>=0) points,
    so invalid/excluded columns leave a gap instead of a spike.
    """
    segment = []
    for x, y in enumerate(rows):
        if y >= 0:
            segment.append((x, y + top_offset))
        else:
            if len(segment) >= 2:
                pts = np.array(segment, dtype=np.int32)
                cv2.polylines(image, [pts], isClosed=False, color=color, thickness=thickness)
            segment = []
    if len(segment) >= 2:
        pts = np.array(segment, dtype=np.int32)
        cv2.polylines(image, [pts], isClosed=False, color=color, thickness=thickness)
    return image


def process_image(input_path, output_path):
    """Run full waterline detection + overlay on a single image file."""
    image = cv2.imread(input_path)
    if image is None:
        print(f"Skipping unreadable file: {input_path}")
        return False

    h, w = image.shape[:2]
    top = int(h * TOP_CROP_FRAC)
    bottom = int(h * (1 - BOTTOM_CROP_FRAC))

    bgr_roi = image[top:bottom, :]
    lab_roi = cv2.cvtColor(bgr_roi, cv2.COLOR_BGR2LAB)
    L_roi = lab_roi[:, :, 0]

    roi_h, roi_w = L_roi.shape
    dune_mask = build_dune_mask(roi_h, roi_w)

    raw_rows = find_waterline(L_roi, dune_mask)
    median_rows = median_filter_1d(raw_rows, MEDIAN_WINDOW)
    smooth_rows = smooth_line(median_rows, SMOOTHING_WINDOW)

    result = draw_waterline(image.copy(), smooth_rows, top_offset=top)

    cv2.imwrite(output_path, result)
    print(f"Saved: {output_path}")
    return True


# ---------------- DEBUG MODE ----------------

def debug_draw_mask_boundary(input_path, output_path):
    """
    Draws the current DUNE_EDGE_POINTS boundary as a BLUE line directly
    on the image, so you can see exactly where the mask sits relative
    to the real dune edge and adjust the config with certainty.
    """
    image = cv2.imread(input_path)
    if image is None:
        raise FileNotFoundError(f"Could not read image at {input_path}")

    h, w = image.shape[:2]
    top = int(h * TOP_CROP_FRAC)
    bottom = int(h * (1 - BOTTOM_CROP_FRAC))
    roi_h = bottom - top

    heights = np.array([p[0] for p in DUNE_EDGE_POINTS])
    x_fracs = np.array([p[1] for p in DUNE_EDGE_POINTS])

    pts = []
    for row in range(roi_h):
        frac = row / max(1, roi_h - 1)
        boundary_x_frac = np.interp(frac, heights, x_fracs)
        boundary_x = int(boundary_x_frac * w)
        pts.append((boundary_x, row + top))

    pts = np.array(pts, dtype=np.int32)
    cv2.polylines(image, [pts], isClosed=False, color=(255, 0, 0), thickness=4)  # blue

    cv2.imwrite(output_path, image)
    print(f"Debug mask boundary saved to: {output_path}")


# ---------------- ENTRY POINTS ----------------

def main_single():
    """Process just one specific file (INPUT_FILENAME)."""
    input_path = os.path.join(IMAGE_DIR, INPUT_FILENAME)
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Could not find image at {input_path}")

    base, ext = os.path.splitext(input_path)
    output_path = f"{base}_waterline{ext}"
    process_image(input_path, output_path)


def main_batch():
    """Process every .jpg in IMAGE_DIR, skipping already-processed outputs."""
    input_paths = glob.glob(os.path.join(IMAGE_DIR, "*.jpg"))

    if not input_paths:
        print(f"No .jpg files found in {IMAGE_DIR}")
        return

    for input_path in input_paths:
        if "_waterline" in input_path or "debug_mask_boundary" in input_path:
            continue  # skip files we already generated

        base, ext = os.path.splitext(input_path)
        output_path = f"{base}_waterline{ext}"
        process_image(input_path, output_path)


def main_debug():
    """Draw the current dune-mask boundary on INPUT_FILENAME for tuning."""
    input_path = os.path.join(IMAGE_DIR, INPUT_FILENAME)
    output_path = os.path.join(IMAGE_DIR, "debug_mask_boundary.jpg")
    debug_draw_mask_boundary(input_path, output_path)


if __name__ == "__main__":
    # ---- CHOOSE ONE MODE BELOW ----

    main_debug()   # STEP 1: run this first, check debug_mask_boundary.jpg,
                   #         adjust DUNE_EDGE_POINTS above until the blue
                   #         line hugs the real dune edge in your image.

    # main_batch()   # STEP 2: once DUNE_EDGE_POINTS looks right in the debug
                     #         image, comment out main_debug() above and
                     #         uncomment this line to process all images.

    # main_single()  # (optional) process just one file instead of the whole folder

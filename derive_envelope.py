#!/usr/bin/env python3
"""
Derive Camera Envelope From Ground Truth
-------------------------------------------
Uses ALL collected ground truth CSVs for a given camera to compute the
actual observed min/max shoreline row per column fraction, with a
safety margin added, and converts directly to CROP-RELATIVE fractions
ready to paste into waterline_detector_v5.py's envelope_points.

Ground truth CSVs store rows in FULL-IMAGE pixel coordinates (see
collect_ground_truth.py / collect_all_ground_truth.py). The detector's
envelope_points are fractions of the CROPPED frame. Converting between
the two by hand is exactly the kind of step that silently produces a
wrong envelope if the arithmetic is off by one crop parameter -- so
this script imports the real CameraProfile (crop_top/crop_bottom) from
waterline_detector_v5.py directly, rather than asking the user to
redo that math.

Requires waterline_detector_v5.py to be importable (same directory,
or on PYTHONPATH).

Usage:
    python3 derive_envelope.py <camera_suffix e.g. c1 or c2>
"""

import sys
import csv
import numpy as np
from pathlib import Path

try:
    import waterline_detector_v5 as detector
except ImportError:
    print("ERROR: could not import waterline_detector_v5.py.")
    print("Run this script from the same directory as waterline_detector_v5.py,")
    print("or add it to PYTHONPATH.")
    sys.exit(1)


GROUND_TRUTH_FOLDER = "/mnt/I2Rgus_Data/waterline/ground_truth"
MARGIN_PX = 60   # extra buffer added above/below observed range, to allow
                 # for tide/conditions not represented in this sample

CAMERA_SUFFIX_TO_KEY = {
    "c1": "CACO05_C1",
    "c2": "CACO05_C2",
}


def load_csv_rows(path):
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        next(reader)
        for line in reader:
            rows.append(int(line[1]))
    return np.array(rows, dtype=np.float32)


def find_matching_image(gt_path):
    """
    Ground truth CSVs are saved with the same stem as their source
    image (see collect_all_ground_truth.py: output_name =
    image_path.stem + ".csv"). Look the source image up in
    INPUT_FOLDER so we can read its real full-image height, instead of
    asking the user to supply it.
    """
    candidate = Path(detector.INPUT_FOLDER) / (gt_path.stem + ".jpg")
    if candidate.exists():
        return candidate
    matches = list(Path(detector.INPUT_FOLDER).glob(gt_path.stem + "*"))
    return matches[0] if matches else None


def parse_margins(argv):
    """
    Returns (margin_above_px, margin_below_px). Defaults to MARGIN_PX
    for both. Allows an ASYMMETRIC margin because the two directions
    guard against different things: the margin ABOVE the observed
    range must cover higher waterlines than the (possibly narrow)
    ground-truth sample happened to capture -- e.g. higher tide, or a
    season not represented in the sample -- and a too-tight ceiling
    there makes the true waterline structurally unreachable by the DP
    solver, which shows up as the line being wrong or 'missing' at
    whichever edge of the frame the waterline sits highest.
    """
    margin_above = margin_below = MARGIN_PX
    forced_height = None
    rest = []
    i = 0
    while i < len(argv):
        if argv[i] == "--margin" and i + 1 < len(argv):
            margin_above = margin_below = float(argv[i + 1]); i += 2
        elif argv[i] == "--margin-above" and i + 1 < len(argv):
            margin_above = float(argv[i + 1]); i += 2
        elif argv[i] == "--margin-below" and i + 1 < len(argv):
            margin_below = float(argv[i + 1]); i += 2
        elif argv[i] == "--image-height" and i + 1 < len(argv):
            forced_height = int(argv[i + 1]); i += 2
        else:
            rest.append(argv[i]); i += 1
    return margin_above, margin_below, forced_height, rest


def main():
    margin_above, margin_below, forced_height, argv_rest = parse_margins(sys.argv[1:])
    sys.argv = [sys.argv[0]] + argv_rest
    if len(sys.argv) != 2:
        print("Usage: python3 derive_envelope.py [--margin N | --margin-above N --margin-below N] "
          "[--image-height N] <camera_suffix e.g. c1 or c2>")
        sys.exit(1)

    camera_suffix = sys.argv[1].lower()
    camera_key = CAMERA_SUFFIX_TO_KEY.get(camera_suffix)
    if camera_key is None or camera_key not in detector.CAMERAS:
        print(f"ERROR: unknown camera suffix '{camera_suffix}'. "
              f"Known: {sorted(CAMERA_SUFFIX_TO_KEY.keys())}")
        sys.exit(1)

    profile = detector.CAMERAS[camera_key]
    gt_dir = Path(GROUND_TRUTH_FOLDER)

    matching_files = sorted([
        f for f in gt_dir.glob("*.csv")
        if f".{camera_suffix}." in f.name.lower() or f.stem.lower().endswith(camera_suffix)
    ])

    if not matching_files:
        print(f"No ground truth files found matching camera '{camera_suffix}'")
        sys.exit(1)

    print(f"Camera profile      : {camera_key} "
          f"(crop_top={profile.crop_top}, crop_bottom={profile.crop_bottom})")
    print(f"Found {len(matching_files)} ground truth file(s) for camera '{camera_suffix}':")
    for f in matching_files:
        print(f"  {f.name}")
    print()

    all_rows = []
    n_columns = None
    crop_top_px = None
    crop_bottom_px = None
    full_height = None

    for f in matching_files:
        rows = load_csv_rows(f)

        if forced_height is not None:
            height = forced_height
        else:
            image_path = find_matching_image(f)
            if image_path is None:
                print(f"  WARNING: no source image found for {f.name} in "
                      f"{detector.INPUT_FOLDER} -- skipping. The source image is only "
                      f"needed to read the full-image HEIGHT; since input/ is cleared "
                      f"every detector run and ImageProducts is wiped nightly, this "
                      f"lookup fails for essentially all historical ground truth. "
                      f"Re-run with --image-height <pixels> (e.g. --image-height 2048) "
                      f"to use this file anyway.")
                continue
            image = detector.load_image(image_path)
            height = image.shape[0]

        this_top_px = int(profile.crop_top * height)
        this_bottom_px = int(profile.crop_bottom * height)

        if full_height is None:
            full_height, crop_top_px, crop_bottom_px = height, this_top_px, this_bottom_px
        elif height != full_height:
            print(f"  WARNING: {f.name}'s source image height ({height}px) differs "
                  f"from earlier images ({full_height}px). Skipping to avoid mixing "
                  f"inconsistent scales.")
            continue

        if n_columns is None:
            n_columns = len(rows)
        elif len(rows) != n_columns:
            print(f"  WARNING: {f.name} has {len(rows)} columns, expected {n_columns}. Skipping.")
            continue

        all_rows.append(rows)

    if not all_rows:
        print("No usable ground truth (with matching source images) found. Aborting.")
        sys.exit(1)

    all_rows = np.array(all_rows)  # shape (n_images, n_columns), full-image px
    n_frames_used = all_rows.shape[0]

    print()
    print(f"FRAMES ACTUALLY USED: {n_frames_used} of {len(matching_files)} ground truth file(s)")
    if n_frames_used < 4:
        print("=" * 100)
        print(f"WARNING: only {n_frames_used} frame(s) contributed to this envelope.")
        print("An envelope derived from very few frames has NO observed spread -- both the")
        print("upper and lower bounds collapse onto whatever those frames happened to show,")
        print("plus the margins. The result can be NARROWER than your existing envelope in")
        print("one direction while widening in the other, silently excluding conditions that")
        print("previously worked. Do not paste these values in without comparing them against")
        print("the envelope currently in waterline_detector_v5.py, bound by bound.")
        print("=" * 100)
    print()

    crop_height_px = crop_bottom_px - crop_top_px
    if crop_height_px <= 0:
        print(f"ERROR: computed crop height <= 0 (top={crop_top_px}, bottom={crop_bottom_px}). "
              f"Check crop_top/crop_bottom for {camera_key}.")
        sys.exit(1)

    n_control_points = 6
    control_x_indices = np.linspace(0, n_columns - 1, n_control_points).astype(int)

    print("=" * 100)
    print(f"Derived envelope for camera '{camera_suffix}' "
          f"(full image height={full_height}px, crop {crop_top_px}-{crop_bottom_px}px, "
          f"margin: -{margin_above:.0f}px above / +{margin_below:.0f}px below)")
    print("=" * 100)
    print(f"{'x_frac':>8} | {'min_row(full)':>14} | {'max_row(full)':>14} | "
          f"{'min_frac(crop)':>16} | {'max_frac(crop)':>16}")

    suggested_points = []

    for idx in control_x_indices:
        column_values = all_rows[:, idx]
        min_row = float(np.min(column_values))
        max_row = float(np.max(column_values))

        min_with_margin = max(0.0, min_row - margin_above)
        max_with_margin = max_row + margin_below

        # Convert from full-image pixel rows to crop-relative fractions,
        # since envelope_points are interpreted relative to the CROPPED
        # frame in build_envelope_mask().
        min_crop_frac = (min_with_margin - crop_top_px) / crop_height_px
        max_crop_frac = (max_with_margin - crop_top_px) / crop_height_px
        min_crop_frac = float(np.clip(min_crop_frac, 0.0, 1.0))
        max_crop_frac = float(np.clip(max_crop_frac, 0.0, 1.0))

        x_frac = idx / (n_columns - 1)

        print(f"{x_frac:8.3f} | {min_with_margin:14.1f} | {max_with_margin:14.1f} | "
              f"{min_crop_frac:16.3f} | {max_crop_frac:16.3f}")

        suggested_points.append((x_frac, min_crop_frac, max_crop_frac))

    print("=" * 100)
    print()
    print("Ready-to-paste anchors for CameraProfile.envelope_points "
          "(already crop-relative -- wrap in make_envelope_points(...)):")
    print("envelope_points=make_envelope_points((")
    for x_frac, min_f, max_f in suggested_points:
        print(f"    ({x_frac:.3f}, {min_f:.3f}, {max_f:.3f}),")
    print(")),")

    if any(f == 0.0 for _, f, _ in suggested_points) or any(f == 1.0 for _, _, f in suggested_points):
        print()
        print("NOTE: one or more fractions clipped to 0.0 or 1.0 -- the observed "
              "ground-truth range (plus margin) extends outside the current crop "
              "region for this camera. Consider widening crop_top/crop_bottom, or "
              "re-run check_crop_coverage.py to confirm.")


if __name__ == "__main__":
    main()


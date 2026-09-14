#!/usr/bin/env python3
"""
Derive Bias Correction From Ground Truth
-------------------------------------------
Aggregates the per-column signed error (detector - ground_truth)
across EVERY matched ground-truth/processed-output pair for a given
camera -- not just one batch, one day, or one small sample -- and
prints ready-to-paste bias_correction_points for that camera's
CameraProfile in waterline_detector_v5.py.

WHY THIS EXISTS: a bias_correction_points fit to a handful of frames
from a single day was shown to only partially generalize -- when run
on a different day's frames, the correction's own residual came back
in nearly the same shape it started with (see project history). The
fix isn't to keep re-deriving 5 numbers from whatever small batch is
on hand; it's to aggregate across as much varied ground truth
(different days, tides, seasons) as exists, and to re-run this
whenever a meaningfully larger/more varied set of ground truth
accumulates, rather than trusting one snapshot indefinitely.

This uses more resolution (default 10 segments, vs. bias_analysis.py's
5) since it's meant to produce the actual correction curve, not just a
diagnostic.

Usage:
    python3 derive_bias_correction.py <camera_suffix e.g. c1 or c2>
"""

import sys
import csv
import warnings
import numpy as np
from pathlib import Path

try:
    import waterline_detector_v5 as detector
except ImportError:
    detector = None  # only needed for display purposes here, not required


GROUND_TRUTH_FOLDER = "/mnt/I2Rgus_Data/waterline/ground_truth"
PROCESSED_FOLDER = "/mnt/I2Rgus_Data/waterline/processed"
N_SEGMENTS = 10


def detect_camera(name):
    name = name.lower()
    if ".c1." in name:
        return "c1"
    if ".c2." in name:
        return "c2"
    return None


def load_csv_rows(path, row_col_index, exclude_no_signal=True):
    """
    Loads a shoreline. Prefers sub-pixel "Row_Precise" when present.

    Columns flagged Has_Signal=0 are returned as NaN rather than
    dropped, so the array stays aligned with ground truth by column
    index. Segment means below use nanmean and therefore skip them.
    This matters: those columns are where the detector found no usable
    evidence, and fitting a bias correction to them would be fitting
    to the fabricated far-field rows -- exactly the data the signal
    filter exists to exclude.
    """
    rows = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        names = reader.fieldnames or []
        has_precise = "Row_Precise" in names
        has_signal = "Has_Signal" in names
        for line in reader:
            if has_precise and line.get("Row_Precise") not in ("", None):
                value = float(line["Row_Precise"])
            else:
                value = float(line["Row"])
            if (exclude_no_signal and has_signal
                    and line.get("Has_Signal") not in ("", None)
                    and int(float(line["Has_Signal"])) == 0):
                value = float("nan")
            rows.append(value)
    return np.array(rows, dtype=np.float32)


def parse_options(argv):
    """
    Directory and filtering options. Defaults preserve the original
    behaviour (snap detections in processed/, ground truth in
    ground_truth/), so existing usage is unchanged.
    """
    opts = {"gt_dir": GROUND_TRUTH_FOLDER, "proc_dir": PROCESSED_FOLDER,
            "suffix": None, "segments": N_SEGMENTS, "exclude_no_signal": True,
            "uncorrected": False}
    rest, i = [], 0
    while i < len(argv):
        a = argv[i]
        if a == "--ground-truth-dir" and i + 1 < len(argv):
            opts["gt_dir"] = argv[i + 1]; i += 2
        elif a == "--processed-dir" and i + 1 < len(argv):
            opts["proc_dir"] = argv[i + 1]; i += 2
        elif a == "--suffix" and i + 1 < len(argv):
            opts["suffix"] = argv[i + 1]; i += 2
        elif a == "--segments" and i + 1 < len(argv):
            opts["segments"] = int(argv[i + 1]); i += 2
        elif a == "--include-no-signal":
            opts["exclude_no_signal"] = False; i += 1
        elif a == "--uncorrected":
            opts["uncorrected"] = True; i += 1
        else:
            rest.append(a); i += 1
    return opts, rest


def main():
    opts, argv_rest = parse_options(sys.argv[1:])
    sys.argv = [sys.argv[0]] + argv_rest
    if len(sys.argv) != 2:
        print("Usage: python3 derive_bias_correction.py [options] <camera_suffix c1|c2>")
        print()
        print("  --ground-truth-dir DIR   default: " + GROUND_TRUTH_FOLDER)
        print("  --processed-dir DIR      default: " + PROCESSED_FOLDER)
        print("  --suffix STR             only use files containing STR, e.g. '.timex.'")
        print("  --segments N             number of horizontal segments (default 10)")
        print("  --include-no-signal      keep columns flagged Has_Signal=0 (default: drop)")
        print("  --uncorrected            detector output had NO bias correction applied,")
        print("                           so the measured error IS the correction. REQUIRED")
        print("                           for timex, whose cron runs --no-bias-correction;")
        print("                           without it the existing snap correction is added")
        print("                           to a residual that never had it subtracted.")
        print()
        print("For timex: --processed-dir archive/processed_timex --suffix .timex.")
        sys.exit(1)

    camera_suffix = sys.argv[1].lower()

    # Look up the camera profile so the correction ALREADY applied by
    # the detector can be read and composed with the measured residual.
    camera_key = {"c1": "CACO05_C1", "c2": "CACO05_C2"}.get(camera_suffix)
    if detector is None or camera_key is None or camera_key not in getattr(detector, "CAMERAS", {}):
        print(f"ERROR: could not load camera profile for '{camera_suffix}'. This script must "
              f"be run from the same directory as waterline_detector_v5.py, because it needs "
              f"to read the bias correction already applied in order to compose the new one "
              f"correctly.")
        sys.exit(1)
    profile = detector.CAMERAS[camera_key]

    gt_dir = Path(opts["gt_dir"])
    processed_dir = Path(opts["proc_dir"])
    n_segments = opts["segments"]

    gt_files = sorted([
        f for f in gt_dir.glob("*.csv")
        if detect_camera(f.name) == camera_suffix
        and (opts["suffix"] is None or opts["suffix"] in f.name)
    ])
    if opts["suffix"]:
        print(f"Filtering to files containing '{opts['suffix']}'")
    print(f"Ground truth dir : {gt_dir}")
    print(f"Processed dir    : {processed_dir}")

    if not gt_files:
        print(f"No ground truth files found matching camera '{camera_suffix}'")
        sys.exit(1)

    print(f"Found {len(gt_files)} ground truth file(s) for camera '{camera_suffix}':")

    all_segment_means = []
    n_columns = None
    used_frames = []
    skipped_frames = []

    for gt_path in gt_files:
        detector_path = processed_dir / (gt_path.stem + ".csv")
        if not detector_path.exists():
            skipped_frames.append((gt_path.name, "no matching processed CSV"))
            continue

        det_rows = load_csv_rows(detector_path, row_col_index=1,
                                 exclude_no_signal=opts["exclude_no_signal"])
        gt_rows = load_csv_rows(gt_path, row_col_index=1, exclude_no_signal=False)

        n = min(len(det_rows), len(gt_rows))
        if n == 0:
            skipped_frames.append((gt_path.name, "empty CSV"))
            continue
        if n_columns is None:
            n_columns = n
        elif n != n_columns:
            skipped_frames.append(
                (gt_path.name, f"column count {n} != expected {n_columns}, skipped to keep segments aligned")
            )
            continue

        signed_error = det_rows[:n] - gt_rows[:n]

        segment_size = n // n_segments
        segment_means = []
        for i in range(n_segments):
            lo = i * segment_size
            hi = n if i == n_segments - 1 else (i + 1) * segment_size
            seg = signed_error[lo:hi]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                value = float(np.nanmean(seg)) if not np.all(np.isnan(seg)) else float("nan")
            segment_means.append(value)

        all_segment_means.append(segment_means)
        used_frames.append(gt_path.name)

    for name in used_frames:
        print(f"  {name}")

    if skipped_frames:
        print()
        print(f"Skipped {len(skipped_frames)} file(s):")
        for name, reason in skipped_frames:
            print(f"  {name}: {reason}")

    if not all_segment_means:
        print()
        print("No usable frame pairs found. Aborting.")
        sys.exit(1)

    all_segment_means = np.array(all_segment_means)  # (n_frames, N_SEGMENTS)
    n_frames = all_segment_means.shape[0]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        segment_avg = np.nanmean(all_segment_means, axis=0)
        segment_std = np.nanstd(all_segment_means, axis=0)
    segment_n = np.sum(~np.isnan(all_segment_means), axis=0)
    segment_avg = np.nan_to_num(segment_avg, nan=0.0)

    segment_size = n_columns // n_segments
    x_fracs = []
    for i in range(n_segments):
        lo = i * segment_size
        hi = n_columns if i == n_segments - 1 else (i + 1) * segment_size
        mid = (lo + hi) / 2.0
        x_fracs.append(mid / (n_columns - 1))

    # COMPOSE with the correction already applied by the detector.
    #
    # This is essential and was previously missing: the detector output
    # these measurements came from ALREADY had
    # profile.bias_correction_points subtracted from it. So
    # segment_avg is the RESIDUAL error remaining AFTER correction,
    # not a correction in its own right. Pasting the residual in as
    # the new bias_correction_points would silently DISCARD whatever
    # the existing correction was contributing and re-introduce that
    # error.
    #
    # Since apply_bias_correction() computes
    #     corrected = raw - correction
    # the corrected error is
    #     err = (raw - correction) - truth
    # and to drive err to zero the correction must become
    #     correction_new = correction_existing + err
    existing = () if opts["uncorrected"] else profile.bias_correction_points
    if opts["uncorrected"]:
        print()
        print("--uncorrected: treating detector output as having NO correction applied.")
        print("  The measured error is therefore the correction directly, with nothing")
        print("  composed in. This is the correct mode for timex detections produced by")
        print("  waterline_timex_cron.sh, which passes --no-bias-correction.")
    if existing:
        ex_x = np.array([pt[0] for pt in existing])
        ex_v = np.array([pt[1] for pt in existing])
        existing_at_x = np.interp(x_fracs, ex_x, ex_v)
    else:
        existing_at_x = np.zeros(len(x_fracs))

    composed = existing_at_x + segment_avg
    raw_error = composed  # identical by construction; the detector's uncorrected error

    print()
    print("=" * 104)
    print(f"Bias correction for camera '{camera_suffix}' "
          f"({n_frames} frame(s) used, {n_segments} segments)")
    print("=" * 104)
    if existing:
        print("NOTE: this camera already has bias_correction_points applied, so the measured")
        print("      error below is the RESIDUAL after that correction. The ready-to-paste")
        print("      block composes existing + residual, which is what you actually want.")
    else:
        print("NOTE: this camera has no existing bias_correction_points, so residual == raw error.")
    print()
    print(f"{'x_frac':>8} | {'existing':>10} | {'residual':>10} | {'=> RAW err':>11} | "
          f"{'NEW corr':>10} | {'std':>8} | {'n':>4}")
    for x, ex, res, comp, std_v, nseg in zip(x_fracs, existing_at_x, segment_avg,
                                             composed, segment_std, segment_n):
        flag = ""
        if nseg < max(3, 0.5 * n_frames):
            flag = "  <-- FEW FRAMES"
        elif std_v > abs(res) and std_v > 10:
            flag = "  <-- HIGH VARIANCE"
        print(f"{x:8.3f} | {ex:10.2f} | {res:10.2f} | {comp:11.2f} | {comp:10.2f} | "
              f"{std_v:8.2f} | {int(nseg):4d}{flag}")
    print()
    print("  n = frames contributing to that segment. Segments where the no-signal filter")
    print("  removed most columns have a low n and their correction is correspondingly")
    print("  less trustworthy -- for c2 this is typically the far field.")
    print("=" * 104)

    spread = float(np.max(raw_error) - np.min(raw_error))
    print(f"RAW error spread across the frame: {spread:.1f}px "
          f"(mean {float(np.mean(raw_error)):+.1f}px)")
    if spread < 25:
        print("  -> The raw error is close to FLAT. A tilted correction would ADD a trend")
        print("     that is not in the data. Prefer a near-constant correction here.")
    else:
        print("  -> The raw error varies substantially across the frame, so a per-column")
        print("     (rather than constant) correction is justified.")
    print()
    print("Ready-to-paste bias_correction_points "
          "(existing + residual; applied via subtraction):")
    print("bias_correction_points=(")
    for x, comp in zip(x_fracs, composed):
        print(f"    ({x:.3f}, {comp:.2f}),")
    print("),")
    print()
    print(f"NOTE: derived from {n_frames} frame(s) currently in {GROUND_TRUTH_FOLDER}. "
          f"Re-run this whenever a meaningfully larger or more varied (different "
          f"days/tides/seasons) set of ground truth becomes available -- a "
          f"correction fit to a narrow sample was already shown to only "
          f"partially generalize to new days.")


if __name__ == "__main__":
    main()

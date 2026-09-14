#!/usr/bin/env python3
"""
Envelope Pinning Check
-------------------------
Reports, per column, how often the detected waterline sits ON the
search envelope boundary rather than somewhere inside it.

WHY: the DP solver can only select rows the envelope allows. If the
true waterline lies outside that band, the solver cannot reach it and
instead returns the best-scoring row it CAN reach -- which is the
boundary itself. Visually this shows up as unnaturally straight
segments and hard right-angle steps, usually hugging one edge of the
frame. It looks like a detection failure but is actually a constraint
failure, and no amount of energy-function tuning will fix it.

This distinguishes the two. A column reported as "pinned" means the
detection is within --tolerance pixels of the envelope edge, i.e. the
solver had no freedom there.

Run it from the directory holding waterline_detector_v5.py.

Usage:
    python3 envelope_pinning_check.py <processed_dir> <camera c1|c2>
        [--tolerance 2.0] [--image-height 2048] [--segments 10]
"""

import sys
import csv
import argparse
from pathlib import Path
from collections import defaultdict

import numpy as np

try:
    import waterline_detector_v5 as detector
except ImportError:
    print("ERROR: cannot import waterline_detector_v5.py. Run this from the same directory.")
    sys.exit(1)


CAMERA_KEYS = {"c1": "CACO05_C1", "c2": "CACO05_C2"}


def load_rows(path):
    """Returns (columns, full_image_rows), preferring sub-pixel Row_Precise."""
    cols, rows = [], []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        precise = reader.fieldnames and "Row_Precise" in reader.fieldnames
        for r in reader:
            cols.append(int(r["Column"]))
            if precise and r["Row_Precise"] not in ("", None):
                rows.append(float(r["Row_Precise"]))
            else:
                rows.append(float(r["Row"]))
    return np.array(cols), np.array(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("processed_dir")
    ap.add_argument("camera", choices=["c1", "c2"])
    ap.add_argument("--tolerance", type=float, default=2.0,
                    help="Pixels from the envelope edge still counted as pinned (default 2).")
    ap.add_argument("--image-height", type=int, default=2048)
    ap.add_argument("--segments", type=int, default=10)
    args = ap.parse_args()

    profile = detector.CAMERAS[CAMERA_KEYS[args.camera]]
    H = args.image_height
    crop_top = int(profile.crop_top * H)
    crop_bottom = int(profile.crop_bottom * H)
    crop_h = crop_bottom - crop_top

    files = sorted(Path(args.processed_dir).glob("*.csv"))
    files = [f for f in files if f".{args.camera}." in f.name.lower()]
    if not files:
        print(f"No {args.camera} CSVs in {args.processed_dir}")
        sys.exit(1)

    # Envelope mask, built once at the width of the first file.
    cols0, _ = load_rows(files[0])
    width = len(cols0)
    mask = detector.build_envelope_mask(crop_h, width, profile)

    lo_edge = np.full(width, np.nan)
    hi_edge = np.full(width, np.nan)
    for c in range(width):
        allowed = np.where(mask[:, c])[0]
        if allowed.size:
            lo_edge[c] = allowed.min()
            hi_edge[c] = allowed.max()

    pinned_lo = np.zeros(width)
    pinned_hi = np.zeros(width)
    counted = np.zeros(width)

    for f in files:
        cols, rows_full = load_rows(f)
        if len(cols) != width:
            continue
        rows_crop = rows_full - crop_top
        near_lo = np.abs(rows_crop - lo_edge) <= args.tolerance
        near_hi = np.abs(rows_crop - hi_edge) <= args.tolerance
        pinned_lo += near_lo
        pinned_hi += near_hi
        counted += 1

    n = int(counted.max()) if counted.size else 0
    if n == 0:
        print("No comparable frames found.")
        sys.exit(1)

    print("=" * 88)
    print(f"ENVELOPE PINNING CHECK -- camera {args.camera}, {len(files)} file(s), "
          f"tolerance {args.tolerance} px")
    print("=" * 88)
    print(f"crop rows {crop_top}-{crop_bottom} (full-image), envelope built at width {width}")
    print()
    print(f"{'segment':>18} | {'cols':>13} | {'pinned TOP':>11} | {'pinned BOTTOM':>13} | "
          f"{'either':>8}")

    seg = width // args.segments
    worst = None
    for i in range(args.segments):
        a = i * seg
        b = width if i == args.segments - 1 else (i + 1) * seg
        frac_lo = pinned_lo[a:b].sum() / (counted[a:b].sum() or 1)
        frac_hi = pinned_hi[a:b].sum() / (counted[a:b].sum() or 1)
        frac_any = min(1.0, frac_lo + frac_hi)
        flag = "  <-- PINNED" if frac_any > 0.25 else ""
        print(f"{i+1:>10} ({a:5d}-{b:5d}) | {b-a:13d} | {100*frac_lo:10.1f}% | "
              f"{100*frac_hi:12.1f}% | {100*frac_any:7.1f}%{flag}")
        if worst is None or frac_any > worst[1]:
            worst = (i + 1, frac_any, frac_lo, frac_hi)

    overall = min(1.0, (pinned_lo.sum() + pinned_hi.sum()) / (counted.sum() or 1))
    print("=" * 88)
    print(f"Overall pinned fraction: {100*overall:.1f}%")
    print()

    if overall > 0.20:
        side = "TOP (envelope ceiling)" if worst[2] > worst[3] else "BOTTOM (envelope floor)"
        print("VERDICT: the envelope is constraining the solution, not just guiding it.")
        print(f"  Worst region is segment {worst[0]} ({100*worst[1]:.0f}% pinned), mostly against")
        print(f"  the {side}.")
        print()
        print("  The true waterline there is outside the allowed band, so the DP returns the")
        print("  boundary instead. This produces straight runs and right-angle steps that look")
        print("  like detection failure but are a constraint failure -- tuning the energy")
        print("  function will not help.")
        print()
        print("  Fix: re-derive the envelope from ground truth collected on the SAME image")
        print("  product being detected. An envelope fitted to snap imagery will not fit timex,")
        print("  because the two put the waterline in different places.")
    else:
        print("VERDICT: the envelope is not materially constraining the solution.")
        print("  Artifacts here are more likely from the energy function or the imagery")
        print("  itself than from the search bounds.")


if __name__ == "__main__":
    main()


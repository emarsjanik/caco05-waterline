#!/usr/bin/env python3
"""
Score Training Images For Usability
=====================================
Measures every image in a pairing manifest for the two defects that a
crop preview showed dominate this training set, and writes a filtered
manifest.

WHY THIS EXISTS: a first OWG model scored 0.51 m RMSE, barely better
than predicting the mean, with strong shrinkage toward the average.
A crop preview then showed that of four frames sampled across the
wave-height range, two were entirely black night frames and the storm
frame was nothing but water droplets on the lens. Each carried a real,
measured wave height. The network was being asked to read waves from
images that contain none -- and for such an image the only rational
answer is the average, which is exactly the shrinkage observed.

An earlier check looked for a few bad days dominating the error and
found none. That check was blind to this defect: night falls on every
day, so dark frames degrade all days evenly and look like a weak model
rather than bad data.

TWO MEASURES:

  * BRIGHTNESS -- mean grey level of the frame. Night and deep dusk
    fall low. The default threshold of 45 is not new: it is the value
    the waterline detector already uses, raised from 20 after three
    dusk frames at 13.5-32.9 passed a lower one.

  * SHARPNESS -- variance of the Laplacian, a standard focus measure.
    Water on the lens, fog and spray remove fine detail and drive it
    down while leaving brightness intact, so a bright but blurred frame
    passes the first test and fails this one. This threshold is NOT
    given a default: it depends on scene content and must be chosen
    from the measured distribution and a look at the frames near it,
    which this script supplies.

Filtering the obscured frames will remove storm conditions, because
storms are when the lens is wet. That thins an already thin tail, and
the script reports it per wave-height bin rather than hiding it. It is
the right trade: a frame the camera cannot see through teaches nothing
about waves -- and it marks a real limit of optical wave gauging, which
cannot measure through a wet lens.

Usage:
    python3 score_image_quality.py --manifest manifest_c2.csv \\
        --image-dir ~/owg_marconi/images --output manifest_c2_scored.csv
    # inspect the report and contact sheets, then:
    python3 score_image_quality.py --manifest manifest_c2.csv \\
        --image-dir ~/owg_marconi/images --output manifest_c2_scored.csv \\
        --min-sharpness 40 --write-clean manifest_c2_clean.csv
"""

import os
import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import cv2


def score(path, work_width=612, crop=None):
    """
    Returns (brightness, sharpness, saturation) or NaNs if unreadable.

    SATURATION is the fraction of pixels at 250 or above. It catches a
    failure the other two miss entirely: low winter sun glinting off
    the water drives a maximum-intensity (bright) composite to pure
    white over large areas. Such a frame is bright, and the edge of the
    white region is crisp so it scores WELL on sharpness, yet it
    carries no usable structure. Measured over the crop region rather
    than the whole frame, because a blown-out sky above the analysis
    band does not matter.

    Scored at reduced resolution: Laplacian variance at full 2448 px is
    dominated by sensor noise and JPEG blocking, and a 4x reduction
    keeps the scene structure that distinguishes a clear frame from an
    obscured one while running far faster.
    """
    img = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
    if img is None:
        return np.nan, np.nan, np.nan
    h, w = img.shape
    small = cv2.resize(img, (work_width, int(h * work_width / w)),
                       interpolation=cv2.INTER_AREA)
    brightness = float(small.mean())
    sharpness = float(cv2.Laplacian(small, cv2.CV_64F).var())
    band = small
    if crop is not None:
        t, b, l, r = crop
        hh, ww = small.shape
        band = small[int(t * hh):int(b * hh), int(l * ww):int(r * ww)]
        if band.size == 0:
            band = small
    saturation = float((band >= 250).mean())
    return brightness, sharpness, saturation


def contact_sheet(rows, image_dir, out_path, title_fn, cols=4, tile=360):
    tiles = []
    for _, r in rows.iterrows():
        img = cv2.imread(os.path.join(image_dir, r["filename"]))
        if img is None:
            continue
        t = cv2.resize(img, (tile, int(tile * img.shape[0] / img.shape[1])))
        cv2.rectangle(t, (0, 0), (tile, 30), (0, 0, 0), -1)
        cv2.putText(t, title_fn(r), (6, 21), cv2.FONT_HERSHEY_SIMPLEX, 0.55,
                    (255, 255, 255), 1, cv2.LINE_AA)
        tiles.append(t)
    if not tiles:
        return False
    th = max(t.shape[0] for t in tiles)
    tiles = [cv2.copyMakeBorder(t, 0, th - t.shape[0], 0, 0, cv2.BORDER_CONSTANT,
                                value=(30, 30, 30)) for t in tiles]
    while len(tiles) % cols:
        tiles.append(np.full_like(tiles[0], 30))
    grid = np.vstack([np.hstack(tiles[i:i + cols]) for i in range(0, len(tiles), cols)])
    cv2.imwrite(out_path, grid, [cv2.IMWRITE_JPEG_QUALITY, 82])
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--image-dir", required=True)
    ap.add_argument("--output", required=True,
                    help="Manifest with brightness and sharpness columns added.")
    ap.add_argument("--min-brightness", type=float, default=45.0,
                    help="Default 45, the waterline detector's proven night threshold.")
    ap.add_argument("--max-saturation", type=float, default=None,
                    help="Reject frames where more than this FRACTION of pixels in the crop "
                         "are at 250 or above (e.g. 0.25). Sun glare blows out "
                         "maximum-intensity composites; such frames pass the brightness and "
                         "sharpness tests while carrying no structure. No default -- choose "
                         "it from the reported distribution and the contact sheet.")
    ap.add_argument("--crop", default=None,
                    help="Crop as 'top,bottom,left,right' fractions, matching the one used "
                         "for training, so saturation is measured where the network looks.")
    ap.add_argument("--min-sharpness", type=float, default=None,
                    help="No default. Choose from the reported distribution and the "
                         "borderline contact sheet.")
    ap.add_argument("--write-clean", default=None,
                    help="Write a manifest of only the frames passing both tests.")
    ap.add_argument("--sheets", default=None,
                    help="Stem for contact sheets (default: next to --output).")
    args = ap.parse_args()

    image_dir = os.path.expanduser(args.image_dir)
    m = pd.read_csv(args.manifest)
    print("=" * 70)
    print("IMAGE QUALITY")
    print("=" * 70)
    print(f"frames            : {len(m)}")

    # Reuse scores from a previous run: the second invocation, with a
    # sharpness threshold chosen, should not re-decode every image.
    crop = None
    if args.crop:
        crop = tuple(float(v) for v in args.crop.split(","))
        if len(crop) != 4:
            sys.exit("--crop needs four values: top,bottom,left,right")
        print(f"crop for saturation: {crop}")

    cached = None
    if Path(args.output).exists():
        prev = pd.read_csv(args.output)
        cols = {"filename", "brightness", "sharpness", "saturation"}
        if cols.issubset(prev.columns) and set(prev["filename"]) >= set(m["filename"]):
            cached = prev.set_index("filename")[["brightness", "sharpness", "saturation"]]
            print("scores            : reused from previous run")

    if cached is None:
        print("scoring ...")
        b, sh, sa = [], [], []
        for i, f in enumerate(m["filename"], 1):
            bb, ss, tt = score(os.path.join(image_dir, f), crop=crop)
            b.append(bb); sh.append(ss); sa.append(tt)
            if i % 250 == 0:
                print(f"  {i}/{len(m)}")
        m["brightness"], m["sharpness"], m["saturation"] = b, sh, sa
    else:
        m = m.join(cached, on="filename")

    unreadable = int(m["brightness"].isna().sum())
    if unreadable:
        print(f"  WARNING: {unreadable} unreadable frame(s) -- excluded")
    m.to_csv(args.output, index=False)
    print(f"wrote             : {args.output}")

    # ---- brightness ------------------------------------------------------
    print()
    print("Brightness (mean grey level)")
    edges = [0, 10, 20, 30, 45, 60, 80, 120, 256]
    counts, _ = np.histogram(m["brightness"].dropna(), bins=edges)
    for i, c in enumerate(counts):
        bar = "#" * int(50 * c / max(counts.max(), 1))
        flag = "  <- below threshold" if edges[i + 1] <= args.min_brightness else ""
        print(f"  {edges[i]:>3}-{edges[i+1]:<3} {c:5d} {bar}{flag}")
    dark = m["brightness"] < args.min_brightness
    print(f"  dark (< {args.min_brightness:g}): {int(dark.sum())} of {len(m)} "
          f"({100*dark.mean():.0f}%)")

    # ---- sharpness, among frames that are bright enough ------------------
    lit = m[~dark & m["brightness"].notna()]
    print()
    print("Sharpness (Laplacian variance), among frames bright enough to keep")
    pct = [1, 5, 10, 25, 50, 75, 90]
    vals = np.percentile(lit["sharpness"], pct)
    print("  " + "  ".join(f"p{p}={v:.0f}" for p, v in zip(pct, vals)))

    # ---- contact sheets for choosing the sharpness threshold -------------
    stem = args.sheets or str(Path(args.output).with_suffix(""))
    lo = lit.nsmallest(16, "sharpness")
    contact_sheet(lo, image_dir, stem + "_least_sharp.jpg",
                  lambda r: f"sharp {r['sharpness']:.0f}  H {r['wave_height_m']:.2f}")
    # frames around the median, as a reference for what "clear" looks like
    mid = lit.iloc[(lit["sharpness"] - lit["sharpness"].median()).abs().argsort()[:8]]
    contact_sheet(mid, image_dir, stem + "_typical.jpg",
                  lambda r: f"sharp {r['sharpness']:.0f}  H {r['wave_height_m']:.2f}")
    print()
    print("Saturation (fraction of crop at 250+), among frames bright enough to keep")
    sv = np.percentile(lit["saturation"], [50, 75, 90, 95, 99])
    print("  " + "  ".join(f"p{p}={v:.3f}" for p, v in zip([50, 75, 90, 95, 99], sv)))
    hi = lit.nlargest(16, "saturation")
    contact_sheet(hi, image_dir, stem + "_most_saturated.jpg",
                  lambda r: f"sat {r['saturation']:.2f}  H {r['wave_height_m']:.2f}")
    print(f"wrote {stem}_most_saturated.jpg (16 most blown-out of the bright frames)")

    print()
    print(f"wrote {stem}_least_sharp.jpg   (16 least sharp of the bright frames)")
    print(f"wrote {stem}_typical.jpg       (8 frames near the median, for reference)")
    print("  Choose --min-sharpness just above the frames that show no beach detail.")

    if args.min_sharpness is None:
        print()
        print("No --min-sharpness given: nothing filtered on sharpness yet. Inspect the")
        print("sheets, then rerun with a threshold and --write-clean.")
        return 0

    # ---- apply both tests, and report what the filter costs -------------
    blurred = ~dark & (m["sharpness"] < args.min_sharpness)
    blown = pd.Series(False, index=m.index)
    if args.max_saturation is not None:
        blown = ~dark & ~blurred & (m["saturation"] > args.max_saturation)
    keep = ~dark & ~blurred & ~blown & m["brightness"].notna()
    print()
    print(f"obscured (sharpness < {args.min_sharpness:g}): {int(blurred.sum())}")
    if args.max_saturation is not None:
        print(f"blown out (saturation > {args.max_saturation:g}): {int(blown.sum())}")
    print(f"KEPT              : {int(keep.sum())} of {len(m)} ({100*keep.mean():.0f}%)")
    print()
    print("What filtering removes, by wave height:")
    print(f"  {'bin (m)':<12}{'total':>7}{'dark':>7}{'obscured':>10}{'blown':>7}{'kept':>7}")
    bins = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 4.0]
    cat = pd.cut(m["wave_height_m"], bins)
    for b in cat.cat.categories:
        sel = cat == b
        print(f"  {b.left:.1f}-{b.right:.1f} m   {int(sel.sum()):>7}{int((sel & dark).sum()):>7}"
              f"{int((sel & blurred).sum()):>10}{int((sel & blown).sum()):>7}"
              f"{int((sel & keep).sum()):>7}")

    storm = m["wave_height_m"] >= 2.0
    if storm.sum():
        frac = (storm & blurred).sum() / storm.sum()
        print()
        print(f"  {100*frac:.0f}% of frames at 2 m and above are obscured, against "
              f"{100*(blurred & ~storm).sum() / max((~storm & ~dark).sum(),1):.0f}% below.")
        print("  Storms are when the lens is wet. This is a limit of optical wave gauging")
        print("  itself -- no model can read waves through water on the lens.")

    if args.write_clean:
        m[keep].drop(columns=["brightness", "sharpness", "saturation"]).to_csv(
            args.write_clean, index=False)
        print()
        print(f"wrote {args.write_clean}   ({int(keep.sum())} usable frames)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

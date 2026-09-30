#!/usr/bin/env python3
"""
Optical Wave Gauge Training Photos In Today's Camera View
============================================================
Redraws every Dec 2024 - Mar 2025 training photo as the CURRENT camera
(calibration 2025-11-13) would have seen it, on the sea surface at that
moment's measured water level, using that day's own camera pointing.

WHY. The training period spans several camera pointings -- horizon_check.py
found c2 re-aimed by 1-2 degrees around Jan 6 and again around Jan 24, and
steady from Jan 24 to Mar 10 but 0.4 deg off the 2025-02-19 calibration.
The network reads a fixed crop of the frame, so each pointing change moved
the sea within that crop by tens to a hundred pixels: the model partly
learned the camera's pointing. Resampled, every photo shows the same sea
in the same pixels -- and in the pixels of TODAY's camera, so the trained
model can read live frames directly, with no conversion at run time.

HOW. Each frame's pointing comes from calibration/chelsea_setups.csv
(written by horizon_check.py --write-setups; days not listed use the
2025-02-19 EO), its water level from the ADCP record on NAVD88. The
resampling is view_reproject.reproject_image, as for the Chelsea waterlines:
exact on the plane z = water level, i.e. the sea surface.

COVERAGE AND CROP. Parts of today's view the old camera never saw are
filled and blurred; a network must not be trained on those. The script
intersects what every pointing used actually saw, and writes the largest
crop rectangle inside it (fractions top,bottom,left,right, the format of
train_marconi_owg.py --crop) to <output-dir>/coverage.json, with a preview
image showing it on a resampled frame.

Usage (on the station):
    python3 resample_owg_images.py --manifest ~/owg_marconi/marconi-c2-clean.csv \\
        --image-dir ~/owg_marconi/images --output-dir ~/owg_marconi/images_current_view
    (manifests with an 'id' or a 'filename' column; --skip-days to leave out
     days whose pointing is in doubt, e.g. 2025-02-16)
"""

import csv
import sys
import json
import argparse
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone

import numpy as np

HERE = Path(__file__).resolve().parent
CAL = HERE / "calibration"
EO_NEW = "20251113"


def names_from_manifest(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return []
    col = "filename" if "filename" in rows[0] else "id" if "id" in rows[0] else None
    if col is None:
        sys.exit(f"{path}: needs a 'filename' or 'id' column")
    out = []
    for r in rows:
        n = r[col].strip()
        out.append(n)
    return sorted(set(out))


def largest_rectangle(mask):
    """Largest all-True axis-aligned rectangle: (r0, r1, c0, c1), end-exclusive."""
    h, w = mask.shape
    heights = np.zeros(w, dtype=int)
    best = (0, 0, 0, 0, 0)
    for r in range(h):
        heights = np.where(mask[r], heights + 1, 0)
        stack = []
        for c in range(w + 1):
            cur = heights[c] if c < w else 0
            start = c
            while stack and stack[-1][1] >= cur:
                s, hh = stack.pop()
                area = hh * (c - s)
                if area > best[0]:
                    best = (area, r - hh + 1, r + 1, s, c)
                start = s
            stack.append((start, cur))
    return best[1:]


def main():
    import cv2
    from georectify import load_extrinsics, load_intrinsics
    from view_reproject import reproject_image, _tables
    from process_chelsea import load_setups, setup_of, load_water_level, level_at

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--manifest", required=True, action="append",
                    help="CSV of training/validation frames (repeatable)")
    ap.add_argument("--image-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--water-level",
                    default="/mnt/I2Rgus_Data/Chelsea_calibration/adcp_water_level_navd88.csv",
                    help="ADCP water level on NAVD88 (from adcp_to_navd88.py)")
    ap.add_argument("--max-gap-minutes", type=float, default=60.0)
    ap.add_argument("--skip-days", nargs="*", default=[],
                    help="UTC dates to leave out, e.g. a day whose pointing is in doubt")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    names = sorted({n for m in args.manifest for n in names_from_manifest(m)})
    src_dir, out_dir = Path(args.image_dir).expanduser(), Path(args.output_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    ep_wl, lv_wl = load_water_level(args.water_level)
    setups = load_setups()
    new_eo, io_by_cam, eo_cache = {}, {}, {}

    done, skipped = 0, defaultdict(int)
    used = defaultdict(list)                     # (cam, eo) -> water levels used
    sample = None
    for i, n in enumerate(names, 1):
        fname = n if n.lower().endswith(".jpg") else n + ".jpg"
        src = src_dir / fname
        if not src.exists():
            skipped["image not found"] += 1
            continue
        parts = fname.split(".")
        cam = "c1" if ".c1." in fname else "c2" if ".c2." in fname else None
        try:
            epoch = int(parts[0])
        except ValueError:
            epoch = None
        if cam is None or epoch is None:
            skipped["name not understood"] += 1
            continue
        day = datetime.fromtimestamp(epoch, tz=timezone.utc).strftime("%Y-%m-%d")
        if day in args.skip_days:
            skipped["--skip-days"] += 1
            continue
        z = level_at(ep_wl, lv_wl, epoch, args.max_gap_minutes * 60)
        if z is None:
            skipped["no water level"] += 1
            continue
        eo_name = setup_of(cam, epoch, setups)
        if cam not in io_by_cam:
            io_by_cam[cam] = load_intrinsics(CAL / f"CACO05_{cam}_20240801_IO.yaml")
            new_eo[cam] = load_extrinsics(CAL / f"CACO05_{cam}_{EO_NEW}_EO-CV.yaml")
        if eo_name not in eo_cache:
            if not (CAL / eo_name).exists():
                sys.exit(f"{eo_name} (for {day}) is not in {CAL}")
            eo_cache[eo_name] = load_extrinsics(CAL / eo_name)
        used[(cam, eo_name)].append(z)
        dst = out_dir / fname
        if dst.exists() and not args.overwrite:
            done += 1
            continue
        img = cv2.imread(str(src))
        if img is None:
            skipped["unreadable"] += 1
            continue
        out = reproject_image(img, io_by_cam[cam], eo_cache[eo_name], new_eo[cam], z)
        cv2.imwrite(str(dst), out, [cv2.IMWRITE_JPEG_QUALITY, 95])
        sample = sample or dst
        done += 1
        if i % 200 == 0:
            print(f"  {i}/{len(names)} ...", flush=True)

    print(f"Resampled {done} of {len(names)} frame(s) into the {EO_NEW} view -> {out_dir}")
    for why, k in sorted(skipped.items()):
        print(f"  skipped {k}: {why}")
    if not used:
        return

    # What every pointing saw, at the lowest and highest water level it was used at.
    coverage = {}
    for cam in sorted({c for c, _ in used}):
        io = io_by_cam[cam]
        mask = None
        lines = []
        for (c, eo_name), zs in sorted(used.items()):
            if c != cam:
                continue
            for z in (min(zs), max(zs)):
                _, _, ok = _tables(tuple(io), tuple(eo_cache[eo_name]), tuple(new_eo[cam]),
                                   int(round(z * 100)))
                ok = ok > 0.5
                mask = ok if mask is None else (mask & ok)
            lines.append(f"    {eo_name}: {len(zs)} frames, water level {min(zs):+.2f} to {max(zs):+.2f} m")
        h, w = mask.shape
        r0, r1, c0, c1 = largest_rectangle(mask)
        crop = [round(r0 / h, 3), round(r1 / h, 3), round(c0 / w, 3), round(c1 / w, 3)]
        coverage[cam] = {"crop_top_bottom_left_right": crop,
                         "seen_by_every_pointing_fraction": round(float(mask.mean()), 3),
                         "pointings": {e: len(zs) for (c, e), zs in used.items() if c == cam}}
        print(f"\n{cam}: {len(lines)} pointing(s)")
        print("\n".join(lines))
        print(f"  seen by every pointing: {100 * mask.mean():.0f}% of today's view")
        print(f"  largest crop inside it (top,bottom,left,right): "
              f"{crop[0]},{crop[1]},{crop[2]},{crop[3]}")
        # preview: a resampled frame with the unseen area darkened and the crop outlined
        ex = next((out_dir / (n if n.endswith('.jpg') else n + '.jpg') for n in names
                   if f".{cam}." in n and (out_dir / (n if n.endswith('.jpg') else n + '.jpg')).exists()), None)
        if ex is not None:
            img = cv2.imread(str(ex))
            H, W = img.shape[:2]
            m = cv2.resize(mask.astype(np.uint8), (W, H), interpolation=cv2.INTER_NEAREST).astype(bool)
            img[~m] = (img[~m] * 0.35).astype(np.uint8)
            cv2.rectangle(img, (int(crop[2] * W), int(crop[0] * H)),
                          (int(crop[3] * W) - 1, int(crop[1] * H) - 1), (0, 255, 255), 6)
            prev = out_dir / f"coverage_{cam}.jpg"
            cv2.imwrite(str(prev), cv2.resize(img, (W // 2, H // 2)))
            print(f"  preview: {prev} (dark = not seen by every pointing; yellow = crop)")
    (out_dir / "coverage.json").write_text(json.dumps(coverage, indent=1) + "\n")
    print(f"\nWrote {out_dir / 'coverage.json'}")


if __name__ == "__main__":
    main()

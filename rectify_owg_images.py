#!/usr/bin/env python3
"""
Wave-Gauge Training Photos Projected Onto A Fixed Patch Of Sea
=================================================================
Projects every Dec 2024 - Mar 2025 training photo onto one rectangle of sea
surface (sea_patch.py), at that moment's measured water level, using that
day's corrected camera pointing. owg_live.py projects live frames onto the
same rectangle with today's calibration, so the network sees the same sea
at the same scale in training and in use.

WHY. resample_owg_images.py redrew the photos into today's camera view
with the horizon-fitted pointings, which before 24 Jan 2025 were ~22 deg
off in pan: the cameras were re-aimed that day, after the Jan 2025 lidar
and before the 2025-02-19 calibration (fit_eo_to_survey.py). With the corrected
pointing, c2 before 24 Jan 2025 saw only ~32% of today's view, so those
redraws put sea from the wrong place in most of the frame; from 24 Jan the
pointing kept the calibration's pan, a few degrees from today's. A patch on the
ground has no such problem -- it is the same patch whatever the pointing,
provided every pointing saw it. For c2 the common patch of SEA is small
(about 50 x 35 m at 195-250 m from the camera, Jan 24 - Mar 2025 plus
today), which is the price of matching the geometry exactly.

The patch is chosen here as the largest near-square rectangle of sea
(--min-seaward past the 0 m contour) that every pointing used -- each
training period's AND today's -- sees at every water level involved, and
written to <output-dir>/patch.json. Copy that file next to the trained model
as <model>.patch.json; owg_live.py then projects live frames onto it.

Usage (on the station, after apply_pointing_correction.py):
    python3 rectify_owg_images.py \\
        --manifest ~/owg_marconi/marconi-c2-bright2-train.csv \\
        --manifest ~/owg_marconi/marconi-c2-bright2-val.csv \\
        --image-dir ~/owg_marconi/images_bright --output-dir ~/owg_marconi/images_patch \\
        --skip-days 2025-02-16 --first-date 2025-01-24
"""

import sys
import json
import argparse
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from resample_owg_images import names_from_manifest          # noqa: E402

CAL = HERE / "calibration"
EO_NEW = "20251113"


def main():
    import cv2
    from georectify import load_extrinsics, load_intrinsics
    from process_chelsea import load_setups, setup_of, load_water_level, level_at
    from sea_patch import Patch, choose_patch, rectify

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--manifest", required=True, action="append")
    ap.add_argument("--image-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--camera", default="c2")
    ap.add_argument("--water-level",
                    default="/mnt/I2Rgus_Data/Chelsea_calibration/adcp_water_level_navd88.csv")
    ap.add_argument("--max-gap-minutes", type=float, default=60.0)
    ap.add_argument("--skip-days", nargs="*", default=[])
    ap.add_argument("--first-date", default=None,
                    help="only frames from this UTC date on, e.g. 2025-01-24: periods pointed far "
                         "from today's camera shrink the patch every pointing sees")
    ap.add_argument("--last-date", default=None)
    ap.add_argument("--patch", default=None, help="use this patch.json instead of choosing one")
    ap.add_argument("--live-levels", nargs=2, type=float, default=[-1.5, 2.0],
                    help="water levels (m NAVD88) today's camera must see the patch at")
    ap.add_argument("--min-seaward", type=float, default=10.0,
                    help="patch must lie this far (m) seaward of the 0 m contour (default 10): "
                         "without it the common view is mostly beach and bluff")
    ap.add_argument("--max-range", type=float, default=250.0,
                    help="farthest patch point from a camera, m (default 250)")
    ap.add_argument("--res", type=float, default=0.25, help="patch resolution, m (default 0.25)")
    ap.add_argument("--min-coverage", type=float, default=0.99,
                    help="Skip frames whose photo shows less than this share of the patch "
                         "(default 0.99): black margins would teach the network something false")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    names = sorted({n for m in args.manifest for n in names_from_manifest(m)})
    src_dir, out_dir = Path(args.image_dir).expanduser(), Path(args.output_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    ep_wl, lv_wl = load_water_level(args.water_level)
    setups = load_setups()
    if not any("_corr_EO" in r[3] or "_lidar_EO" in r[3] for r in setups if r[0] == args.camera):
        print("WARNING: no corrected pointings in chelsea_setups.csv for this camera -- run "
              "apply_pointing_correction.py first, or the frames are projected with the "
              "uncorrected pointing.")
    io = load_intrinsics(CAL / f"CACO05_{args.camera}_20240801_IO.yaml")
    eo_now = load_extrinsics(CAL / f"CACO05_{args.camera}_{EO_NEW}_EO-CV.yaml")

    # frames -> (pointing, level)
    todo, skipped, used = [], defaultdict(int), defaultdict(list)
    eos = {}
    for n in names:
        fname = n if n.lower().endswith(".jpg") else n + ".jpg"
        if f".{args.camera}." not in fname:
            skipped["other camera"] += 1
            continue
        src = src_dir / fname
        if not src.exists():
            skipped["image not found"] += 1
            continue
        try:
            epoch = int(fname.split(".")[0])
        except ValueError:
            skipped["name not understood"] += 1
            continue
        day = datetime.fromtimestamp(epoch, tz=timezone.utc).strftime("%Y-%m-%d")
        if day in args.skip_days:
            skipped["--skip-days"] += 1
            continue
        if (args.first_date and day < args.first_date) or (args.last_date and day > args.last_date):
            skipped["outside --first-date/--last-date"] += 1
            continue
        z = level_at(ep_wl, lv_wl, epoch, args.max_gap_minutes * 60)
        if z is None:
            skipped["no water level"] += 1
            continue
        eo_name = setup_of(args.camera, epoch, setups)
        if eo_name not in eos:
            eos[eo_name] = load_extrinsics(CAL / eo_name)
        used[eo_name].append(z)
        todo.append((fname, src, eo_name, z))

    print(f"frames            : {len(todo)} of {len(names)}")
    for why, k in sorted(skipped.items()):
        print(f"  skipped {k}: {why}")
    for e, zs in sorted(used.items()):
        print(f"  {e}: {len(zs)} frames, water level {min(zs):+.2f} to {max(zs):+.2f} m")

    if args.patch:
        patch = Patch.from_json(args.patch)
    else:
        views = [(io, eos[e], [min(zs), max(zs)]) for e, zs in used.items()]
        views.append((io, eo_now, list(args.live_levels)))
        patch, area = choose_patch(views, res=args.res, min_seaward=args.min_seaward,
                                   max_range=args.max_range)
        print(f"patch             : {patch.along[1] - patch.along[0]:.0f} x "
              f"{patch.cross[1] - patch.cross[0]:.0f} m ({area:.0f} m2) seen by all "
              f"{len(views)} pointings incl. today's; {patch.shape[1]} x {patch.shape[0]} px")
    c = patch.corners()
    rng = np.hypot(c[:, 0] - eo_now[0], c[:, 1] - eo_now[1])
    from sea_patch import DEFAULT_ORIGIN, DEFAULT_SHORE_BEARING
    sb = np.radians(DEFAULT_SHORE_BEARING + 90.0)
    sea = (c[:, 0] - DEFAULT_ORIGIN[0]) * np.sin(sb) + (c[:, 1] - DEFAULT_ORIGIN[1]) * np.cos(sb)
    print(f"                    {sea.min():.0f}-{sea.max():.0f} m seaward of the Jan 2025 0 m contour")
    print(f"                    corners (UTM 19N): " +
          "; ".join(f"{e:.0f},{n:.0f}" for e, n in c) + f"  -- {rng.min():.0f}-{rng.max():.0f} m from the camera")
    (out_dir / "patch.json").write_text(json.dumps(patch.to_json(), indent=1) + "\n")

    done, low = 0, 0
    sample = {}
    for i, (fname, src, eo_name, z) in enumerate(todo, 1):
        dst = out_dir / fname
        if dst.exists() and not args.overwrite:
            done += 1
            continue
        img = cv2.imread(str(src))
        if img is None:
            continue
        out, cov = rectify(img, io, eos[eo_name], z, patch)
        if cov < args.min_coverage:
            low += 1
            if dst.exists():
                dst.unlink()
            continue
        cv2.imwrite(str(dst), out, [cv2.IMWRITE_JPEG_QUALITY, 95])
        sample.setdefault(eo_name, (src, z))
        done += 1
        if i % 200 == 0:
            print(f"  {i}/{len(todo)} ...", flush=True)
    print(f"projected         : {done} frame(s) -> {out_dir} ({patch.shape[1]} x {patch.shape[0]} px)"
          + (f"; {low} SKIPPED, the photo shows less than {args.min_coverage:.0%} of the patch"
             if low else ""))

    # check picture: the patch outline on one photo per pointing
    from view_reproject import ground_to_pixel
    for eo_name, (src, z) in sample.items():
        img = cv2.imread(str(src))
        U, V, ok = ground_to_pixel(np.r_[c[:, 0], c[0, 0]], np.r_[c[:, 1], c[0, 1]], z, io, eos[eo_name])
        if ok.all():
            pts = np.c_[U, V].astype(np.int32).reshape(-1, 1, 2)
            cv2.polylines(img, [pts], False, (0, 0, 0), 9, cv2.LINE_AA)
            cv2.polylines(img, [pts], False, (0, 255, 255), 4, cv2.LINE_AA)
        prev = out_dir / f"patch_on_{Path(eo_name).stem}.jpg"
        cv2.imwrite(str(prev), cv2.resize(img, (img.shape[1] // 2, img.shape[0] // 2)),
                    [cv2.IMWRITE_JPEG_QUALITY, 85])
    print(f"check pictures    : {out_dir}/patch_on_*.jpg (yellow = the patch on each pointing's photo)")

    # training lists with only the frames that were projected
    import csv
    kept = {f for f, *_ in todo if (out_dir / f).exists()}
    for m in args.manifest:
        rows = list(csv.DictReader(open(Path(m).expanduser(), newline="")))
        if not rows:
            continue
        col = "filename" if "filename" in rows[0] else "id"
        keep = [r for r in rows if (r[col] if r[col].lower().endswith(".jpg") else r[col] + ".jpg") in kept]
        dst = out_dir / (Path(m).stem + "_patch.csv")
        with open(dst, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(keep)
        print(f"training list     : {dst} ({len(keep)} of {len(rows)} rows)")
    print(f"\nNext: train on {out_dir} with --img-size {patch.shape[1]} --img-height {patch.shape[0]} "
          f"and no --crop; copy patch.json next to the model as <model>.patch.json.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

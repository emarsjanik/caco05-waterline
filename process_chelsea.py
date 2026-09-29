#!/usr/bin/env python3
"""
Process The ADCP-Deployment Imagery (CACO03/CACO04, pre-Nov-2025 geometry)
============================================================================
Runs the waterline chain on imagery from before the Nov 2025 camera move,
using the Signature 1000 ADCP as the water level, and writes everything to
one output folder:

    <output>/adcp_water_level_navd88.csv   ADCP water level on NAVD88
    <output>/geometry_check.csv            did the cameras move during the period?
    <output>/contour_points.csv            waterlines with elevation, CURRENT-view pixels
    <output>/contour_points_original_view.csv   same, pixels of the ORIGINAL photos
    <output>/contour_points_ground.csv     same, UTM 19N easting/northing
    <output>/maps/<date>_<cam>.png         each day's waterlines on the photo, by elevation
    <output>/maps/all_<cam>.png            every waterline in the period on one photo
    <output>/overlays/<cam>/<frame>.jpg    each frame with its waterline (half size)
    <output>/dem/<period>_dem.*            DEM per image period, and for all periods
    <output>/processing.log                every stage's full output

HOW, AND WHY IT DIFFERS FROM process_historical.py

  * Geometry. These frames were taken with the 2025-02-19 extrinsics
    (calibration/CACO05_c?_20250219_EO.yaml); the detector is tuned in
    pixel positions of the 2025-11-13 view, ~300-350 px away. Each frame
    is first resampled into the 2025-11-13 view on the plane z = its
    water level (view_reproject.py). The waterline lies on that plane,
    so detecting on the resampled frame and georectifying with the
    2025-11-13 calibration at that z gives the true ground position.
    Detections are mapped back for the photo overlays.

  * Water level. Signature 1000 ADCP (sig1000_waves_ALL.csv), converted
    to NAVD88 by adcp_to_navd88.py (run automatically if the file is not
    in the output folder yet). It covers 2024-12-09 18:00 to 2025-03-10
    15:00 UTC; frames outside it, or more than --max-gap-minutes from an
    ADCP reading on either side, are skipped and counted.

  * Camera stability. The 2025-02-19 solve is only right for frames taken
    with the same pointing. geometry_check.csv gives each day's image
    shift against a 2025-02-19 frame (phase correlation). Shifts beyond
    GEOMETRY_WARN_PX mean that day needs its own calibration; they are
    reported, not silently used.

  * Search envelope. The detector only looks inside a per-column band
    tuned on the 2026 beach. If the winter 2025 beach put the waterline
    well outside it, frames are rejected as "no usable signal" (counted
    in processing.log) rather than mis-detected. Many such rejections on
    otherwise clear days mean the envelope needs widening for this period.

  * Input layout. Downloaded folders may hold duplicate copies (the same
    period with and without a leading '__'), c1/c2 subfolders and
    quarantine_* subfolders. Every *.timex.jpg below --images is used
    once, quarantine folders excepted.

Usage (on the station):
    cd /mnt/I2Rgus_Data/waterline
    python3 process_chelsea.py --images /mnt/I2Rgus_Data/Chelsea_pics \\
        --output /mnt/I2Rgus_Data/Chelsea_calibration
"""

import os
import re
import csv
import sys
import shutil
import argparse
import subprocess
from pathlib import Path
from datetime import datetime, timezone
from collections import defaultdict

import numpy as np

HERE = Path(__file__).resolve().parent
CAL = HERE / "calibration"
EO_OLD = "20250219"
EO_NEW = "20251113"
GEOMETRY_WARN_PX = 10.0
GEOMETRY_MIN_CORRELATION = 0.2   # weaker matches say nothing either way
PERIOD_RE = re.compile(r"(\d{4})_(\d{2})_(\d{2})_to_(\d{4})_(\d{2})_(\d{2})")
EPOCH_RE = re.compile(r"^(\d{9,11})\.")


def say(msg=""):
    print(msg, flush=True)


def run(cmd, log_path, label):
    say(f"  {label} ...")
    with open(log_path, "a") as log:
        log.write(f"\n{'=' * 70}\n{label}\n{' '.join(str(c) for c in cmd)}\n{'=' * 70}\n")
        log.flush()
        proc = subprocess.run([str(c) for c in cmd], stdout=log, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        say(f"    FAILED (exit {proc.returncode}) -- see {log_path}")
        return False
    return True


def camera_of(name):
    n = name.lower()
    return "c1" if ".c1." in n else "c2" if ".c2." in n else None


def period_of(path, images_root):
    """'2025_02_13_to_2025_02_24' for the top-level folder a frame came from."""
    top = path.relative_to(images_root).parts[0]
    m = PERIOD_RE.search(top)
    return f"{m.group(1)}_{m.group(2)}_{m.group(3)}_to_{m.group(4)}_{m.group(5)}_{m.group(6)}" if m else top.lstrip("_")


def collect_images(images_root):
    frames = {}
    for p in sorted(images_root.rglob("*.timex.jpg")):
        if any(part.lower().startswith("quarantine") for part in p.parts):
            continue
        m = EPOCH_RE.match(p.name)
        cam = camera_of(p.name)
        if not m or not cam or p.name in frames:
            continue
        frames[p.name] = {"path": p, "epoch": int(m.group(1)), "camera": cam,
                          "period": period_of(p, images_root)}
    return frames


def load_water_level(path):
    ep, lv = [], []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            t = datetime.strptime(r["time"], "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
            ep.append(t.timestamp())
            lv.append(float(r["water_level_navd88"]))
    order = np.argsort(ep)
    return np.array(ep)[order], np.array(lv)[order]


def level_at(ep_wl, lv_wl, epoch, max_gap_s):
    """Linear interpolation; None if either bracketing reading is more than max_gap_s away."""
    i = np.searchsorted(ep_wl, epoch)
    if i == 0 or i == len(ep_wl):
        return None
    if max(epoch - ep_wl[i - 1], ep_wl[i] - epoch) > max_gap_s:
        return None
    w = (epoch - ep_wl[i - 1]) / (ep_wl[i] - ep_wl[i - 1])
    return float(lv_wl[i - 1] + w * (lv_wl[i] - lv_wl[i - 1]))


def geometry_check(frames, out_csv):
    """
    Image shift of each day's midday frame against the 2025-02-19 frame
    nearest 17:00 UTC, per camera. Gradient images at quarter size and a
    window keep the comparison on fixed edges (dune line, horizon,
    structures) more than on moving water.
    """
    import cv2
    ref_day = datetime(2025, 2, 19, 17, tzinfo=timezone.utc).timestamp()
    rows = []

    def prep(path):
        g = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
        if g is None:
            return None
        g = cv2.resize(g, (g.shape[1] // 4, g.shape[0] // 4), interpolation=cv2.INTER_AREA)
        g = cv2.GaussianBlur(g.astype(np.float32), (0, 0), 1.5)
        mag = cv2.magnitude(cv2.Sobel(g, cv2.CV_32F, 1, 0), cv2.Sobel(g, cv2.CV_32F, 0, 1))
        return mag * cv2.createHanningWindow(mag.shape[::-1], cv2.CV_32F)

    for cam in ("c1", "c2"):
        mine = [f for f in frames.values() if f["camera"] == cam]
        if not mine:
            continue
        ref = min(mine, key=lambda f: abs(f["epoch"] - ref_day))
        ref_img = prep(ref["path"])
        by_day = defaultdict(list)
        for f in mine:
            by_day[datetime.fromtimestamp(f["epoch"], tz=timezone.utc).date()].append(f)
        for day, fs in sorted(by_day.items()):
            noon = datetime(day.year, day.month, day.day, 17, tzinfo=timezone.utc).timestamp()
            f = min(fs, key=lambda f: abs(f["epoch"] - noon))
            img = prep(f["path"])
            if img is None or ref_img is None or img.shape != ref_img.shape:
                rows.append([cam, day.isoformat(), f["path"].name, "", "", "", "unreadable or different size"])
                continue
            (dx, dy), resp = cv2.phaseCorrelate(ref_img, img)
            shift = 4 * float(np.hypot(dx, dy))
            if resp < GEOMETRY_MIN_CORRELATION:
                note = "inconclusive (weak match: fog, glare or a very different scene)"
            elif shift > GEOMETRY_WARN_PX:
                note = "MOVED? needs its own calibration"
            else:
                note = ""
            rows.append([cam, day.isoformat(), f["path"].name, round(4 * dx, 1), round(4 * dy, 1),
                         round(resp, 3), note])
        say(f"  {cam}: reference {ref['path'].name}")
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["camera", "date", "frame", "shift_x_px", "shift_y_px", "correlation", "note"])
        w.writerows(rows)
    moved = [r for r in rows if r[6].startswith("MOVED")]
    unsure = [r for r in rows if r[6].startswith("inconclusive")]
    if unsure:
        say(f"  {len(unsure)} camera-day(s) inconclusive (weak image match) -- see geometry_check.csv")
    return rows, moved


def back_to_original_view(contours_in, contours_out, levels, calib, current_view_out):
    """
    Rewrite pixel_column/pixel_row from the current view into the original
    photo's pixels (contours_out). Points the old camera never saw -- i.e.
    detected in the nearest-filled margin of a resampled frame -- are not
    real and are dropped from BOTH outputs; current_view_out keeps the
    surviving points in current-view pixels for georectification.
    """
    from view_reproject import dst_to_src_points
    with open(contours_in, newline="") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames
        rows = list(reader)
    current = [dict(r) for r in rows]
    groups = defaultdict(list)
    for i, r in enumerate(rows):
        groups[r["source_file"]].append(i)
    keep = np.ones(len(rows), bool)
    for src, idx in groups.items():
        cam = rows[idx[0]]["camera"]
        z = levels.get(src, float(rows[idx[0]]["tide_elevation_navd88"]))
        io, eo_old, eo_new = calib[cam]
        u = np.array([float(rows[i]["pixel_column"]) for i in idx])
        v = np.array([float(rows[i]["pixel_row"]) for i in idx])
        su, sv, ok = dst_to_src_points(u, v, z, io, eo_old, eo_new)
        for k, i in enumerate(idx):
            if ok[k]:
                rows[i]["pixel_column"] = f"{su[k]:.1f}"
                rows[i]["pixel_row"] = f"{sv[k]:.1f}"
            else:
                keep[i] = False
    for path, table in ((contours_out, rows), (current_view_out, current)):
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for r, k in zip(table, keep):
                if k:
                    w.writerow(r)
    return int(keep.sum()), int((~keep).sum())


def draw_overlays(contours_csv, originals, out_dir, max_frames=None):
    """Each frame with its waterline, coloured by elevation, at half size."""
    import cv2
    import matplotlib
    matplotlib.use("Agg")
    try:
        cmap = matplotlib.colormaps["turbo"]
    except AttributeError:                     # matplotlib < 3.5
        from matplotlib import cm
        cmap = cm.get_cmap("turbo")
    frames = defaultdict(lambda: {"u": [], "v": [], "z": None, "cam": None})
    with open(contours_csv, newline="") as f:
        for r in csv.DictReader(f):
            d = frames[r["source_file"]]
            d["u"].append(float(r["pixel_column"]))
            d["v"].append(float(r["pixel_row"]))
            d["z"] = float(r["tide_elevation_navd88"])
            d["cam"] = r["camera"]
    if not frames:
        return 0
    zs = [d["z"] for d in frames.values()]
    zmin, zmax = min(zs), max(zs)
    n = 0
    for src, d in sorted(frames.items())[:max_frames]:
        img_path = next(iter(originals.glob(src + "*")), None)
        img = cv2.imread(str(img_path)) if img_path else None
        if img is None:
            continue
        order = np.argsort(d["u"])
        pts = np.c_[np.array(d["u"])[order], np.array(d["v"])[order]].astype(np.int32)
        frac = 0.5 if zmax == zmin else (d["z"] - zmin) / (zmax - zmin)
        r, g, b, _ = cmap(frac)
        colour = (int(b * 255), int(g * 255), int(r * 255))
        cv2.polylines(img, [pts.reshape(-1, 1, 2)], False, (0, 0, 0), 9, cv2.LINE_AA)
        cv2.polylines(img, [pts.reshape(-1, 1, 2)], False, colour, 5, cv2.LINE_AA)
        label = f"{d['z']:+.2f} m NAVD88"
        cv2.putText(img, label, (40, 90), cv2.FONT_HERSHEY_SIMPLEX, 2.4, (0, 0, 0), 10, cv2.LINE_AA)
        cv2.putText(img, label, (40, 90), cv2.FONT_HERSHEY_SIMPLEX, 2.4, colour, 4, cv2.LINE_AA)
        img = cv2.resize(img, (img.shape[1] // 2, img.shape[0] // 2), interpolation=cv2.INTER_AREA)
        target = out_dir / d["cam"]
        target.mkdir(parents=True, exist_ok=True)
        cv2.imwrite(str(target / (src + ".jpg")), img, [cv2.IMWRITE_JPEG_QUALITY, 85])
        n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--images", default="/mnt/I2Rgus_Data/Chelsea_pics")
    ap.add_argument("--output", default="/mnt/I2Rgus_Data/Chelsea_calibration")
    ap.add_argument("--adcp", default=str(HERE / "sig1000_waves_ALL.csv"))
    ap.add_argument("--max-gap-minutes", type=float, default=60.0,
                    help="Largest distance to the ADCP reading on either side (ADCP is hourly).")
    ap.add_argument("--cell", type=float, default=2.0, help="DEM cell size, m.")
    ap.add_argument("--max-hs", type=float, default=1.5,
                    help="Leave frames with ADCP Hs above this out of the DEMs (default 1.5, as "
                         "the station cron). Waves push the timex waterline up the beach; "
                         "validate_waves.py part B measures by how much. 0 keeps every frame.")
    ap.add_argument("--max-day-offset", type=float, default=0.15,
                    help="Leave out whole days whose elevations differ from the other days' "
                         "DEM by more than this (m, default 0.15, as the station cron): catches "
                         "days where the detector followed snow or shadow. 0 disables.")
    ap.add_argument("--envelope-pad", type=float, default=0.1,
                    help="Widen the detector's search envelope by this fraction of the crop "
                         "height (default 0.1). The mass rejection of the winter frames was "
                         "caused by the black border of the resampled frames, not the envelope; "
                         "at 0.3 the detector wandered onto the bluff's shadow on the dry beach.")
    ap.add_argument("--utc-hours", default="13-18",
                    help="Only use frames from this UTC hour range, start inclusive, end "
                         "exclusive (default 13-18, i.e. 8 am - 1 pm EST). In winter the low "
                         "afternoon sun throws the bluff's shadow across the beach, and the "
                         "detector followed its edge (and snow patches) instead of the water. "
                         "'0-24' uses every frame.")
    ap.add_argument("--min-signal-fraction", type=float, default=0.35,
                    help="Share of a frame's columns that must show the waterline (default "
                         "0.35; the station uses 0.60). On overcast winter frames the far field "
                         "is too hazy: a real 21 Jan frame had signal in 95-100%% of the near and "
                         "middle columns but 20-40%% of the far third, 44%% overall. Columns "
                         "without signal are dropped one by one, so only visible parts are used.")
    ap.add_argument("--no-overlays", action="store_true", help="Skip the per-frame overlay images.")
    ap.add_argument("--skip-geometry-check", action="store_true")
    args = ap.parse_args()

    import cv2
    from georectify import load_extrinsics, load_intrinsics
    from view_reproject import reproject_image

    images_root = Path(args.images)
    out = Path(args.output)
    work = out / "work"
    # Products are rebuilt every run; leftovers from a run with other
    # settings (e.g. afternoon frames) would otherwise sit among them.
    for d in (out / "maps", out / "dem", out / "overlays"):
        shutil.rmtree(d, ignore_errors=True)
    for d in (out, work / "src", work / "in", work / "detections", work / "debug",
              work / "original", out / "maps", out / "dem"):
        d.mkdir(parents=True, exist_ok=True)
    log = out / "processing.log"
    log.write_text(f"process_chelsea.py started {datetime.now().isoformat()}\n")

    calib = {}
    for cam in ("c1", "c2"):
        calib[cam] = (load_intrinsics(CAL / f"CACO05_{cam}_20240801_IO.yaml"),
                      load_extrinsics(CAL / f"CACO05_{cam}_{EO_OLD}_EO.yaml"),
                      load_extrinsics(CAL / f"CACO05_{cam}_{EO_NEW}_EO-CV.yaml"))

    # --- 1. images --------------------------------------------------------
    say("1. Collecting images")
    frames = collect_images(images_root)
    if not frames:
        sys.exit(f"No *.timex.jpg found below {images_root}")
    by_period = defaultdict(lambda: defaultdict(int))
    for f in frames.values():
        by_period[f["period"]][f["camera"]] += 1
    for p in sorted(by_period):
        say(f"  {p}: c1 {by_period[p]['c1']}, c2 {by_period[p]['c2']}")
    say(f"  {len(frames)} unique frames (duplicate copies and quarantine folders skipped)")

    # --- 2. water level ----------------------------------------------------
    say()
    say("2. ADCP water level")
    wl_csv = out / "adcp_water_level_navd88.csv"
    if not wl_csv.exists():
        if not run([sys.executable, HERE / "adcp_to_navd88.py", "--adcp", args.adcp,
                    "--output", wl_csv], log, "converting the ADCP record to NAVD88"):
            sys.exit("Could not make the ADCP water level. Run adcp_to_navd88.py by hand "
                     "(add --offset 0.09 if the NOAA gauge cannot be reached).")
    ep_wl, lv_wl = load_water_level(wl_csv)
    # The ADCP's own wave record, in the format fetch_buoy_waves.py writes,
    # so every frame is tagged with the waves measured AT Marconi.
    waves_csv = out / "adcp_waves.csv"
    import pandas as pd
    wv = pd.read_csv(args.adcp)
    wt = pd.to_datetime(wv["time"], utc=True)
    pd.DataFrame({"time_utc": wt.dt.strftime("%Y-%m-%dT%H:%M:%SZ"),
                  "epoch": (wt - pd.Timestamp("1970-01-01", tz="UTC")) // pd.Timedelta(seconds=1),
                  "wvht_m": wv["wh_4061"].round(3), "dpd_s": wv["wp_peak"].round(2),
                  "mwd_deg": wv["wvdir"]}).to_csv(waves_csv, index=False)
    say(f"  {len(ep_wl)} readings, "
        f"{datetime.fromtimestamp(ep_wl[0], tz=timezone.utc):%Y-%m-%d %H:%M} to "
        f"{datetime.fromtimestamp(ep_wl[-1], tz=timezone.utc):%Y-%m-%d %H:%M} UTC")
    h0, h1 = (int(x) for x in args.utc_hours.split("-"))
    wrong_hour = [n for n, f in frames.items()
                  if not h0 <= datetime.fromtimestamp(f["epoch"], tz=timezone.utc).hour < h1]
    for n in wrong_hour:
        del frames[n]
    if wrong_hour:
        say(f"  {len(wrong_hour)} frame(s) outside {h0:02d}:00-{h1:02d}:00 UTC skipped "
            f"(afternoon shadow of the bluff; --utc-hours)")
    levels, no_level = {}, defaultdict(int)
    for name, f in frames.items():
        z = level_at(ep_wl, lv_wl, f["epoch"], args.max_gap_minutes * 60)
        if z is None:
            no_level[f["period"]] += 1
        else:
            levels[name] = z
    say(f"  {len(levels)} frames have an ADCP water level")
    for p, n in sorted(no_level.items()):
        say(f"  {p}: {n} frame(s) skipped, outside the ADCP record")
    if not levels:
        sys.exit("No frame falls inside the ADCP record.")

    # --- 3. camera stability ---------------------------------------------------
    say()
    if args.skip_geometry_check:
        say("3. Camera stability check skipped")
        moved = []
    else:
        say("3. Camera stability against the 2025-02-19 calibration")
        _, moved = geometry_check({n: f for n, f in frames.items() if n in levels},
                                  out / "geometry_check.csv")
        if moved:
            say(f"  WARNING: {len(moved)} camera-day(s) shifted more than {GEOMETRY_WARN_PX:.0f} px "
                f"from the 2025-02-19 view -- see geometry_check.csv. Their ground positions "
                f"(and the DEM) are only as good as the 2025-02-19 calibration is for them.")
        else:
            say(f"  all days within {GEOMETRY_WARN_PX:.0f} px of the 2025-02-19 view")

    # --- 4. reproject into the current view -------------------------------------
    say()
    say("4. Resampling frames into the 2025-11-13 view on the water surface")
    from view_reproject import REPROJECT_VERSION
    marker = work / "src" / ".reproject_version"
    if not marker.exists() or marker.read_text().strip() != REPROJECT_VERSION:
        shutil.rmtree(work / "src", ignore_errors=True)       # made by an older resampling
        (work / "src").mkdir(parents=True)
        marker.write_text(REPROJECT_VERSION + "\n")
    # The detector reads every frame in work/src, so frames cached by an
    # earlier run with other settings (e.g. --utc-hours) must go.
    for sub in ("src", "original"):
        for p in (work / sub).glob("*.jpg"):
            if p.name not in levels:
                p.unlink()
    with open(work / "frame_levels.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "camera", "water_level_navd88"])
        for i, name in enumerate(sorted(levels), 1):
            f = frames[name]
            orig = work / "original" / name
            if not orig.exists():
                try:
                    os.link(f["path"], orig)
                except OSError:
                    shutil.copy2(f["path"], orig)
            dst = work / "src" / name
            if not dst.exists():
                img = cv2.imread(str(f["path"]))
                if img is None:
                    continue
                io, eo_old, eo_new = calib[f["camera"]]
                cv2.imwrite(str(dst), reproject_image(img, io, eo_old, eo_new, levels[name]),
                            [cv2.IMWRITE_JPEG_QUALITY, 95])
            w.writerow([name, f["camera"], f"{levels[name]:.4f}"])
            if i % 200 == 0:
                say(f"  {i} / {len(levels)}")
    say(f"  {len(levels)} frames resampled")

    # --- 5-7. detect, water level, georectify ---------------------------------
    say()
    say("5. Detecting waterlines (roughly 3-4 s per frame)")
    # A rerun must not mix in detections from an earlier run's settings.
    for d in (work / "in", work / "detections", work / "debug"):
        shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True)
    if args.envelope_pad:
        say(f"  search envelope widened by {args.envelope_pad:.2f} of the crop height")
    if not run([sys.executable, HERE / "waterline_detector_v5.py",
                "--envelope-pad", args.envelope_pad,
                "--min-signal-fraction", args.min_signal_fraction,
                "--image-suffix", "timex.jpg",
                "--source-dir", work / "src", "--input-dir", work / "in",
                "--output-dir", work / "detections", "--debug-dir", work / "debug"],
               log, "detecting waterlines"):
        sys.exit(f"Detection failed; see {log}")
    say(f"    {len(list((work / 'detections').glob('*.csv')))} frame(s) passed the detector's quality filters")

    all_contours = work / "contour_points_all.csv"
    if not run([sys.executable, HERE / "extract_elevation_contours.py", wl_csv,
                "--time-col", "time", "--level-col", "water_level_navd88",
                "--max-gap-minutes", args.max_gap_minutes,
                "--min-coverage", args.min_signal_fraction,
                "--waves", waves_csv, "--wave-max-gap-minutes", 60,
                "--processed-dir", work / "detections", "--output", all_contours],
               log, "matching water levels and filtering contours") or not all_contours.exists():
        sys.exit(f"Contour extraction failed; see {log}")

    # --- 6. drop detections outside the original photo, map the rest back ----------------
    say()
    say("6. Keeping only waterline points the original photo actually shows")
    level_by_stem = {Path(n).name.rsplit(".jpg", 1)[0]: z for n, z in levels.items()}
    contours = out / "contour_points.csv"
    orig_contours = out / "contour_points_original_view.csv"
    n_pts, n_out = back_to_original_view(all_contours, orig_contours, level_by_stem, calib, contours)
    say(f"  {n_pts} point(s) kept; {n_out} dropped (in the filled margin the old camera did not see)")

    # --- 7. georectify ------------------------------------------------------------------
    ground = out / "contour_points_ground.csv"
    if not run([sys.executable, HERE / "georectify.py", contours, ground,
                "--io-c1", CAL / "CACO05_c1_20240801_IO.yaml",
                "--eo-c1", CAL / f"CACO05_c1_{EO_NEW}_EO-CV.yaml",
                "--io-c2", CAL / "CACO05_c2_20240801_IO.yaml",
                "--eo-c2", CAL / f"CACO05_c2_{EO_NEW}_EO-CV.yaml"],
               log, "georectifying to UTM Zone 19 (2025-11-13 geometry = the resampled view)"):
        sys.exit(f"Georectification failed; see {log}")

    # --- 9. maps and overlays --------------------------------------------------------
    say()
    say("7. Maps")
    dates = defaultdict(set)
    with open(orig_contours, newline="") as fh:
        for r in csv.DictReader(fh):
            dates[r["camera"]].add(r["capture_time_utc"][:10])
    for cam in sorted(dates):
        run([sys.executable, HERE / "daily_elevation_map.py", orig_contours, work / "original",
             cam, out / "maps" / f"all_{cam}.png", "--days", "36500"],
            log, f"all waterlines, {cam}")
        for d in sorted(dates[cam]):
            run([sys.executable, HERE / "daily_elevation_map.py", orig_contours, work / "original",
                 cam, out / "maps" / f"{d}_{cam}.png", "--date", d],
                log, f"waterlines {d} {cam}")
    say(f"  {len(list((out / 'maps').glob('*.png')))} map(s) in {out / 'maps'}")
    if not args.no_overlays:
        n = draw_overlays(orig_contours, work / "original", out / "overlays")
        say(f"  {n} per-frame overlay(s) in {out / 'overlays'}")

    # --- 10. DEMs ---------------------------------------------------------------------
    say()
    say("8. DEMs")
    dem_filters = []
    if args.max_hs:
        dem_filters += ["--max-hs", args.max_hs]
    if args.max_day_offset:
        dem_filters += ["--max-day-offset", args.max_day_offset]
    if dem_filters:
        say(f"  filters: {' '.join(str(x) for x in dem_filters)}")
    periods = sorted({f["period"] for n, f in frames.items() if n in levels})
    made = []
    for p in periods:
        m = PERIOD_RE.search(p)
        if not m:
            continue
        start = f"{m.group(1)}-{m.group(2)}-{m.group(3)}"
        end = f"{m.group(4)}-{m.group(5)}-{m.group(6)}"
        if run([sys.executable, HERE / "dem_from_contours.py", ground, out / "dem" / p,
                "--cell", args.cell, "--start-date", start, "--end-date", end] + dem_filters,
               log, f"DEM {start} to {end}"):
            made.append(p)
    if run([sys.executable, HERE / "dem_from_contours.py", ground, out / "dem" / "all_periods",
            "--cell", args.cell] + dem_filters, log, "DEM, all periods together"):
        made.append("all_periods")

    say()
    say("=" * 70)
    say(f"DONE -> {out}")
    say("=" * 70)
    say(f"  DEMs: {', '.join(made) if made else 'none (see processing.log)'}")
    if moved:
        say(f"  NOTE: {len(moved)} camera-day(s) flagged in geometry_check.csv")
    say(f"  Full log: {log}")


if __name__ == "__main__":
    main()

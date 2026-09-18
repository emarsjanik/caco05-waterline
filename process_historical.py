#!/usr/bin/env python3
"""
Historical Imagery Processing -- CACO05 Marconi Beach
========================================================
Interactive driver for reprocessing archived TIMEX imagery on a machine
other than the station NUC.

Prompts for three folders -- images, tide model, output -- validates
each before doing any work, and then runs the full chain:

    detect -> water level -> contours -> georectify -> DEM + maps

WHAT IT CHECKS BEFORE PROCESSING, and why each check exists:

  * Calibration validity. The extrinsic solve is dated 2025-11-13 and
    followed a physical camera move. Imagery from before that date has
    a different geometry, and the current calibration would turn it
    into ground coordinates that look entirely plausible and are wrong.
    The script refuses rather than warns.

  * Image product. Snap and timex place the waterline in different
    places -- C1's raw error is -49.8 px on snap against -11.0 px on
    timex -- so the bias corrections are product-specific. This script
    processes TIMEX ONLY and says so if it finds anything else.

  * Tide coverage against image dates. A mismatch here produces an
    empty result with no error: every frame is silently skipped for
    having no nearby water level. Checked up front instead.

  * Dependencies and companion scripts. Failing at the first prompt
    with a clear message beats failing twenty minutes into a run.

VERTICAL DATUM. A tide MODEL is not the same measurement as GNSS-R.
Compared over a 47-day overlap at this site the model sat 0.2574 m
above the GNSS-R water level, so a -0.25 m offset is applied by
default to bring modelled levels onto the same datum as the
operational product. Without it every elevation is about a quarter of
a metre too high. Override with care.

Usage:
    python3 process_historical.py                 (prompts for everything)
    python3 process_historical.py --images DIR --tides DIR --output DIR
"""

import os
import re
import sys
import glob
import shutil
import argparse
import subprocess
from pathlib import Path
from datetime import datetime, timezone

# The camera move this calibration follows. Imagery before this date
# needs a different extrinsic solve.
CALIBRATION_VALID_FROM = "2025-11-13"

# Tide model minus GNSS-R, measured over a 47-day overlap at this site.
DEFAULT_DATUM_OFFSET = -0.25

REQUIRED_SCRIPTS = [
    "waterline_detector_v5.py",
    "extract_elevation_contours.py",
    "georectify.py",
    "dem_from_contours.py",
    "daily_elevation_map.py",
]

REQUIRED_CALIBRATION = [
    "CACO05_c1_20240801_IO.yaml",
    "CACO05_c2_20240801_IO.yaml",
    "CACO05_c1_20251113_EO-CV.yaml",
    "CACO05_c2_20251113_EO-CV.yaml",
]

HERE = Path(__file__).resolve().parent


# ----------------------------------------------------------------- utils

def say(msg=""):
    print(msg, flush=True)


def fail(msg):
    say()
    say("STOPPED: " + msg)
    sys.exit(1)


def ask_folder(prompt, must_exist=True, create=False):
    while True:
        raw = input(f"{prompt}\n  > ").strip().strip('"').strip("'")
        if not raw:
            say("  (nothing entered)")
            continue
        path = Path(os.path.expanduser(raw)).resolve()
        if create:
            try:
                path.mkdir(parents=True, exist_ok=True)
                return path
            except Exception as exc:
                say(f"  cannot create that folder: {exc}")
                continue
        if must_exist and not path.is_dir():
            say(f"  not a folder: {path}")
            continue
        return path


def confirm(prompt):
    while True:
        r = input(f"{prompt} [y/n] > ").strip().lower()
        if r in ("y", "yes"):
            return True
        if r in ("n", "no"):
            return False


def parse_stem_date(name):
    """Date from an Argus-style filename, or None."""
    m = re.search(r"\.([A-Z][a-z]{2})\.(\d{1,2})_\d{2}_\d{2}_\d{2}\.GMT\.(\d{4})\.", name)
    if m:
        month, day, year = m.group(1), int(m.group(2)), int(m.group(3))
        try:
            return datetime.strptime(f"{year} {month} {day}", "%Y %b %d").date().isoformat()
        except ValueError:
            return None
    m = re.match(r"^(\d{9,11})\.", name)
    if m:
        try:
            return datetime.fromtimestamp(int(m.group(1)), tz=timezone.utc).date().isoformat()
        except (ValueError, OSError):
            return None
    return None


# ----------------------------------------------------- preflight checks

def check_environment():
    say("Checking environment...")
    missing = []
    for mod, why in [("numpy", "arrays"), ("cv2", "image reading (opencv-python)"),
                     ("matplotlib", "map and DEM plots"), ("pandas", "tide spreadsheets"),
                     ("openpyxl", "reading .xlsx")]:
        try:
            __import__(mod)
        except ImportError:
            missing.append(f"    {mod:<12} -- {why}")
    if missing:
        say("  Missing Python packages:")
        for m in missing:
            say(m)
        say()
        say("  Install with:")
        say("    pip install numpy opencv-python matplotlib pandas openpyxl")
        fail("required packages are not installed.")

    missing_scripts = [s for s in REQUIRED_SCRIPTS if not (HERE / s).exists()]
    if missing_scripts:
        say(f"  These scripts must sit beside this one in {HERE}:")
        for s in missing_scripts:
            say(f"    {s}")
        fail("companion scripts not found.")

    cal_dir = HERE / "calibration"
    missing_cal = [c for c in REQUIRED_CALIBRATION if not (cal_dir / c).exists()]
    if missing_cal:
        say(f"  These calibration files must sit in {cal_dir}:")
        for c in missing_cal:
            say(f"    {c}")
        fail("calibration not found -- georectification and the DEM cannot run.")

    say("  OK: packages, scripts and calibration all present.")
    return cal_dir


def scan_images(folder):
    """Reports what is in the image folder and returns the timex files."""
    all_jpg = sorted(glob.glob(str(folder / "*.jpg")))
    if not all_jpg:
        fail(f"no .jpg files found in {folder}")

    timex, other_products = [], {}
    for f in all_jpg:
        name = Path(f).name.lower()
        if ".timex." in name:
            timex.append(f)
        else:
            for product in ("snap", "bright", "dark", "var", "ras"):
                if f".{product}." in name:
                    other_products[product] = other_products.get(product, 0) + 1
                    break

    cams = {}
    dates = []
    for f in timex:
        n = Path(f).name
        for c in ("c1", "c2"):
            if f".{c}." in n.lower():
                cams[c] = cams.get(c, 0) + 1
        d = parse_stem_date(n)
        if d:
            dates.append(d)

    say(f"  files found        : {len(all_jpg)}")
    say(f"  timex              : {len(timex)}")
    if other_products:
        listed = ", ".join(f"{k} {v}" for k, v in sorted(other_products.items()))
        say(f"  other products     : {listed}  (IGNORED -- this script does timex only,")
        say(f"                       because bias corrections differ by product)")
    if not timex:
        fail("no timex images in that folder. This script processes timex only.")
    say(f"  cameras            : " + ", ".join(f"{k} {v}" for k, v in sorted(cams.items())))
    if dates:
        say(f"  date range         : {min(dates)} to {max(dates)}  ({len(set(dates))} day(s))")
    else:
        say("  date range         : could not be read from filenames")

    return timex, sorted(set(dates))


def check_calibration_dates(dates):
    if not dates:
        say()
        say("  WARNING: no dates could be read from the filenames, so the calibration")
        say("  validity check could not run. The extrinsics date from "
            f"{CALIBRATION_VALID_FROM}")
        say("  and followed a camera move; earlier imagery would be georectified WRONG,")
        say("  and wrong in a way that looks entirely plausible.")
        return confirm("  Continue anyway?")
    too_old = [d for d in dates if d < CALIBRATION_VALID_FROM]
    if too_old:
        say()
        say(f"  {len(too_old)} day(s) fall BEFORE {CALIBRATION_VALID_FROM}, "
            f"earliest {min(too_old)}.")
        say("  The extrinsic calibration follows a physical camera move on that date, so")
        say("  those frames have a different geometry. Georectifying them with this")
        say("  calibration produces coordinates that are wrong but look reasonable.")
        say("  A separate extrinsic solve is needed for the earlier period.")
        fail("imagery predates the calibration.")
    say(f"  calibration        : valid (all frames on or after {CALIBRATION_VALID_FROM})")


def merge_tide_files(folder, out_path):
    """
    Combines every tide spreadsheet in the folder into one file.

    The extraction step takes a single water-level file, but coverage
    here is typically split across several spreadsheets by month.
    Merging on 'time' and dropping duplicates gives one continuous
    record and avoids running the extraction repeatedly and stitching
    the results.
    """
    import pandas as pd

    files = sorted(glob.glob(str(folder / "*.xlsx")) + glob.glob(str(folder / "*.csv")))
    files = [f for f in files if not Path(f).name.startswith("~$")]
    if not files:
        fail(f"no .xlsx or .csv tide files found in {folder}")

    say(f"  files found        : {len(files)}")
    frames = []
    for f in files:
        try:
            df = pd.read_excel(f) if f.lower().endswith(".xlsx") else pd.read_csv(f)
        except Exception as exc:
            say(f"    skipped {Path(f).name}: {exc}")
            continue
        cols = {c.lower(): c for c in df.columns}
        if "time" not in cols:
            say(f"    skipped {Path(f).name}: no 'time' column")
            continue
        height_cols = [c for c in df.columns if "heightm" in c.lower()]
        if not height_cols:
            say(f"    skipped {Path(f).name}: no '<model>_heightm' column")
            continue
        say(f"    {Path(f).name}: {len(df)} row(s), "
            f"{len(height_cols)} model(s) [{', '.join(height_cols)}]")
        frames.append(df)

    if not frames:
        fail("no usable tide files -- each needs a 'time' column and at least one "
             "'<model>_heightm' column.")

    merged = pd.concat(frames, ignore_index=True)
    time_col = [c for c in merged.columns if c.lower() == "time"][0]
    merged[time_col] = pd.to_datetime(merged[time_col])
    merged = merged.sort_values(time_col).drop_duplicates(subset=[time_col], keep="first")
    merged.to_excel(out_path, index=False)

    lo = merged[time_col].min().date().isoformat()
    hi = merged[time_col].max().date().isoformat()
    say(f"  merged             : {len(merged)} row(s), {lo} to {hi}")
    return lo, hi


def check_overlap(image_dates, tide_lo, tide_hi):
    if not image_dates:
        return
    covered = [d for d in image_dates if tide_lo <= d <= tide_hi]
    uncovered = [d for d in image_dates if not (tide_lo <= d <= tide_hi)]
    if not covered:
        say()
        say(f"  Image dates   : {min(image_dates)} to {max(image_dates)}")
        say(f"  Tide coverage : {tide_lo} to {tide_hi}")
        fail("tide data does not cover any of the imagery. Every frame would be "
             "skipped for having no nearby water level, and the run would finish "
             "reporting success with an empty result.")
    if uncovered:
        say(f"  coverage           : {len(covered)} of {len(image_dates)} day(s) covered; "
            f"{len(uncovered)} outside the tide record will be skipped")
    else:
        say(f"  coverage           : all {len(covered)} image day(s) covered")


# ------------------------------------------------------------- pipeline

def run(cmd, log_path, label):
    say(f"  {label} ...")
    with open(log_path, "a") as log:
        log.write(f"\n{'='*70}\n{label}\n{' '.join(str(c) for c in cmd)}\n{'='*70}\n")
        log.flush()
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        say(f"    FAILED (exit {proc.returncode}) -- see {log_path}")
        return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--images"); ap.add_argument("--tides"); ap.add_argument("--output")
    ap.add_argument("--datum-offset", type=float, default=DEFAULT_DATUM_OFFSET)
    ap.add_argument("--cell", type=float, default=2.0)
    ap.add_argument("--yes", action="store_true", help="Skip the confirmation prompt.")
    args = ap.parse_args()

    say("=" * 70)
    say("CACO05 HISTORICAL IMAGERY PROCESSING")
    say("=" * 70)
    say()
    cal_dir = check_environment()
    say()

    # --- images -------------------------------------------------------
    if args.images:
        img_dir = Path(os.path.expanduser(args.images)).resolve()
        if not img_dir.is_dir():
            fail(f"not a folder: {img_dir}")
    else:
        say("Where are the images?")
        img_dir = ask_folder("  Full path to the folder containing the timex .jpg files")
    say()
    say(f"Scanning {img_dir}")
    timex_files, image_dates = scan_images(img_dir)
    check_calibration_dates(image_dates)
    say()

    # --- output (needed before merging tides, which writes there) -----
    if args.output:
        out_dir = Path(os.path.expanduser(args.output)).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
    else:
        say("Where should the results go?")
        say("  A new folder is best -- this writes several files and a working")
        say("  subfolder, and keeping runs separate makes them easy to compare.")
        out_dir = ask_folder("  Full path to the output folder (created if absent)",
                             must_exist=False, create=True)
    say(f"  output folder      : {out_dir}")
    say()

    # --- tides --------------------------------------------------------
    if args.tides:
        tide_dir = Path(os.path.expanduser(args.tides)).resolve()
        if not tide_dir.is_dir():
            fail(f"not a folder: {tide_dir}")
    else:
        say("Where is the tide model data?")
        say("  A folder of .xlsx files with a 'time' column and one or more")
        say("  '<model>_heightm' columns. Several files covering different")
        say("  periods is fine -- they will be merged.")
        tide_dir = ask_folder("  Full path to the tide model folder")
    say()
    say(f"Scanning {tide_dir}")
    merged_tide = out_dir / "_tides_merged.xlsx"
    tide_lo, tide_hi = merge_tide_files(tide_dir, merged_tide)
    check_overlap(image_dates, tide_lo, tide_hi)
    say()

    # --- summary ------------------------------------------------------
    say("-" * 70)
    say("READY TO PROCESS")
    say("-" * 70)
    say(f"  images        : {len(timex_files)} timex frames from {img_dir}")
    if image_dates:
        say(f"                  {min(image_dates)} to {max(image_dates)}")
    say(f"  water level   : tide MODEL, {tide_lo} to {tide_hi}")
    say(f"  datum offset  : {args.datum_offset:+.2f} m")
    say( "                  (a tide model sits about 0.26 m above GNSS-R at this site;")
    say( "                   without this every elevation would be that much too high)")
    say(f"  output        : {out_dir}")
    say(f"  DEM cell      : {args.cell} m")
    say()
    say("  Will produce: elevation maps with waterlines for c1 and c2, an")
    say("  intertidal DEM with repeatability and sample-count grids, and the")
    say("  intermediate contour files.")
    say()
    if not args.yes and not confirm("Proceed?"):
        say("Cancelled.")
        return 0
    say()

    # --- run ----------------------------------------------------------
    work = out_dir / "work"
    src = work / "src"
    for d in (work, src, work / "in", work / "detections", work / "debug"):
        d.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "processing.log"
    log_path.write_text(f"CACO05 historical processing\nstarted {datetime.now().isoformat()}\n")

    say("Staging timex images...")
    for f in timex_files:
        target = src / Path(f).name
        if not target.exists():
            try:
                os.link(f, target)          # hard link: no copy, no extra disk
            except OSError:
                shutil.copy2(f, target)
    say(f"  {len(timex_files)} file(s) staged")
    say()

    say("Processing (this takes roughly 3-4 seconds per image)...")

    ok = run([sys.executable, str(HERE / "waterline_detector_v5.py"),
              "--image-suffix", "timex.jpg",
              "--source-dir", str(src),
              "--input-dir", str(work / "in"),
              "--output-dir", str(work / "detections"),
              "--debug-dir", str(work / "debug")],
             log_path, "detecting waterlines")
    if not ok:
        fail(f"detection failed. See {log_path}")

    n_det = len(list((work / "detections").glob("*.csv")))
    say(f"    {n_det} detection(s) passed the quality filters")
    if n_det == 0:
        fail("no frames produced a usable waterline. The log will say which filter "
             "rejected them -- most often darkness, fog, or no usable signal.")

    contours = out_dir / "contour_points.csv"
    ok = run([sys.executable, str(HERE / "extract_elevation_contours.py"),
              str(merged_tide),
              "--processed-dir", str(work / "detections"),
              "--output", str(contours),
              "--vertical-datum-offset", str(args.datum_offset)],
             log_path, "matching water levels and filtering contours")
    if not ok or not contours.exists():
        fail(f"contour extraction failed. See {log_path}")

    ground = out_dir / "contour_points_ground.csv"
    ok = run([sys.executable, str(HERE / "georectify.py"), str(contours), str(ground),
              "--io-c1", str(cal_dir / "CACO05_c1_20240801_IO.yaml"),
              "--eo-c1", str(cal_dir / "CACO05_c1_20251113_EO-CV.yaml"),
              "--io-c2", str(cal_dir / "CACO05_c2_20240801_IO.yaml"),
              "--eo-c2", str(cal_dir / "CACO05_c2_20251113_EO-CV.yaml")],
             log_path, "georectifying to UTM Zone 19")
    if not ok:
        fail(f"georectification failed. See {log_path}")

    made_map = []
    for cam in ("c1", "c2"):
        out_png = out_dir / f"elevation_map_{cam}.png"
        if run([sys.executable, str(HERE / "daily_elevation_map.py"),
                str(contours), str(src), cam, str(out_png),
                "--days", "36500"],
               log_path, f"drawing waterline map for {cam}"):
            if out_png.exists():
                made_map.append(out_png.name)

    made_dem = False
    if run([sys.executable, str(HERE / "dem_from_contours.py"),
            str(ground), str(out_dir / "dem"), "--cell", str(args.cell)],
           log_path, "building the DEM"):
        made_dem = (out_dir / "dem_dem.asc").exists()

    # --- report -------------------------------------------------------
    say()
    say("=" * 70)
    say("DONE")
    say("=" * 70)
    say(f"  {out_dir}")
    say()
    for name in made_map:
        say(f"    {name:<28} waterlines drawn on a photo, coloured by elevation")
    if made_dem:
        say(f"    {'dem_dem.png':<28} DEM and repeatability, side by side")
        say(f"    {'dem_dem.asc':<28} elevation grid (opens in QGIS)")
        say(f"    {'dem_spread.asc':<28} repeatability grid")
        say(f"    {'dem_count.asc':<28} samples per cell")
    say(f"    {'contour_points.csv':<28} waterlines in pixel coordinates")
    say(f"    {'contour_points_ground.csv':<28} the same, in UTM Zone 19")
    say(f"    {'processing.log':<28} full output of every stage")
    say()
    say("  The 'work' subfolder holds staged images and per-frame detections.")
    say("  It can be deleted once you are happy with the results.")
    say()
    say("  REMEMBER: elevations come from a tide MODEL, not a measurement, so they")
    say("  carry the model's own error on top of everything else. The detected edge")
    say("  is also the wave RUNUP limit, which sits above still water by an amount")
    say("  that grows with wave height -- a systematic high bias, not corrected here.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        say()
        say("Interrupted.")
        sys.exit(130)

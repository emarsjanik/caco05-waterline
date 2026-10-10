#!/usr/bin/env python3
"""
Survey-Date Products: Waterlines, DEM, Maps And An Honest Comparison
======================================================================
One command per survey date builds, from the station photos around that
date, the waterlines, an intertidal DEM (ASCII grids and GeoTIFFs), the
DEM page, waterline maps, a comparison with that date's survey labelled
by how independent it really is, and a provenance file that says exactly
what went in. --summary puts every date built so far in one table and
one figure.

WHY. A DEM from the cameras is worth something only when it agrees with
a survey that played no part in making it. Every number that builds the
DEM -- the camera pointing, the water level, the waves and the setup
coefficient, the detection search envelope -- has a source, and some of
those sources are surveys. This script runs the whole chain the same way
for every date, from one configuration table (surveys/survey_dates.csv),
writes down every source, and labels each comparison by its weakest step:

  INDEPENDENT      nothing in the chain was fitted to or placed with
                   this survey
  CROSS-VALIDATED  fitted on other data from the same source, held out
                   here
  PARTLY-CIRCULAR  some step used this survey (e.g. the search envelope
                   was placed with it, or the calibration was solved from
                   these very points)
  CIRCULAR         pointing or offsets fitted to this survey

The label in surveys/surveys.csv is a starting point; the script then
looks at the chain it actually ran and can only make it worse: a search
envelope placed with the same survey -> PARTLY-CIRCULAR; a point survey
of the calibration's own day (maybe its GCPs) -> PARTLY-CIRCULAR; a
camera pointing fitted to it (fit_eo_to_survey.py) -> CIRCULAR; a setup
coefficient fitted to it -> CIRCULAR (C = 0.037 was fitted to the
2026-09-29 RTK shots, so with it the 2026-09-29 RTK comparison is
CIRCULAR; the shots are recognised by their content, not their file
name). survey_compare.py's own checks run before the comparison and the
worst of all is the headline label, with every reason. A check that can
misfire (a date or name rule) is overruled only by a reason written in
surveys.csv (column label_override_reason), which is printed next to
the label. The setup coefficient is read from the waterlines the DEM was
built from, never from the command line of a later --steps compare.
A --setup-coef other than the known fits must say where it was fitted
(--setup-fitted-to <survey file> or 'none:<how>'), or it is refused.

HOW, per date (--steps, in this order; each step is skipped when its
outputs exist, are newer than its inputs and were made with the same
settings, the same scripts (sha256) and the same upstream builds, unless
--force; every skip is printed; a rebuilt step makes every step after it
stale, and nothing stale is ever reported as current):
  forcing   water level and waves. Oct 2024 - Mar 2025 (era adcp or
            chatham): historical_forcing.py (the ADCP where it measured,
            Chatham 8447435 transferred to Marconi elsewhere; ADCP waves, or
            WIS/NDBC converted to the ADCP's Hs). Live dates (era live): the
            GNSS-R levels and waves the station's own archive rows carry.
  pointing  every day of photos, per camera: the sea horizon against the
            calibration (horizon_check.py's fit) and the rotation from a
            reference frame (pointing_check.py's feature matching; the
            reference is a frame from the calibration day when the window
            has one). A day whose pointing differs is LEFT OUT of the
            detection and listed (--keep-moved-days keeps it); if most days
            differ from the reference nothing is left out and a warning says
            to look. Nothing is dropped silently.
  detect    historical eras: detect_original_view.py per camera with the
            station ID, the date's calibration, the search envelope from the
            OTHER survey, the forcing files and --setup-coef (the detector's
            debug images, ~2 GB a two-camera week, are then removed unless
            --keep-detector-debug). Live era: the
            cron's own detections (contour_points_timex.csv) for the window,
            with the setup correction computed here from each row's
            offshore_hs_m / offshore_tp_s by the formula of
            extract_elevation_contours.py --setup-coef (same rounding, same
            columns), then georectify.py with the live calibration. Then the
            cameras are merged into waterlines/contour_points_ground.csv,
            each row tagged with its water-level and wave source.
  filter    waterline_consistency.py (its defaults, as the live cron).
            If it fails, the unfiltered file is used and that is said.
  dem       dem_from_contours.py over the window (the cron's DEM settings;
            its day test is not run on fewer than 3 camera-days, where it
            would reject both of two days that disagree, and the README says
            so), asc_to_geotiff.py (EPSG:32619; NAVD88 EPSG:5703 vertical keys
            on the DEM), and the DEM page by dem_figure.py with the date's own
            camera positions.
  maps      daily_elevation_map.py per camera on the window's photos, and a
            plan view of the filtered waterlines coloured by elevation over
            the survey (or the DEM).
  compare   survey_compare.py for every survey of the date in
            surveys/surveys.csv. A survey whose file is missing is skipped
            with a message saying which file to provide and where.
Then provenance.json (calibration, water level and wave sources, setup
coefficient and where it was fitted, filter and DEM settings, search
envelope, photo list, pointing results, and the commit and script
sha256s that built each step's outputs) and README.txt (what this is,
its label and why, its caveats, the exact commands).

STATUS AND EXIT CODE. complete, exit 0. partial, exit 1: a step failed
(a failed comparison included). partial, exit 3: a camera contributed
nothing, a part is missing (a camera's pointing check or its waterline
map on the photos: no photo found, e.g. a wrong photo root; a survey file
surveys.csv names that was not found, or has no current comparison), or
the outputs on disk are not one build (e.g. a step rebuilt without the
steps after it). A disabled date exits 2. A build marks the product 'in
progress' (provenance.json, and a banner over the previous README) before
any step runs and rewrites both at the end; one interrupted on the way
(Ctrl-C, a full disk) writes status 'failed' with the reason, or leaves
the mark if even that cannot be written: the previous build's numbers are
never presented as current. --summary re-checks every comparison against
the stamps on disk. One build of a date at a time: a build holds
<date>/.build.lock while it runs; a second one refuses to start (exit 4).
--all returns the code of its worst date by severity: 1 if any date failed,
else 4 if any was busy, else 3 if any was partial, else 0.

WHAT A LATER RUN CANNOT CHANGE. The outputs on disk keep the settings they
were built with: the DEM settings, the setup coefficient, the photo window
and hours, and the input files. A later run that does not rebuild the steps
using an option (e.g. --steps compare --window ..., or a --live-contours
whose step fails) says NOT APPLIED, and the README, provenance and the
"rebuild everything" line report what built the outputs. --window and
--utc-hours apply only to a run of every step from the forcing to the DEM.
--synthetic NOTE (or a SYNTHETIC_FIXTURE file in a photo root) marks a build
from synthetic test inputs on every output.

DISABLED DATES. A date whose rows in survey_dates.csv all have enabled=0
is refused with the reason, unless --force-disabled.

RUN TIME: forcing 10-30 s; pointing ~3 s per day and camera; detection
~25 s per photo (detect_original_view.py on the whole frame: measured 22-26
s per photo on a 2.8 GHz Xeon, not yet on the station's NUC, whose
progress lines print the real rate), i.e. ~50 min for a two-camera week
of daytime photos (~120); filter, DEM and GeoTIFFs under a minute; maps
~20 s per camera; each lidar comparison 15-40 s. The live era reads the
station's whole contour file (~1-2 min) instead of detecting. The run
lowers its own priority (--nice 10) so the station's cron keeps the two
cores; it says when it would overlap the cron (12:30, 19:25, 20:55 local)
and cleanup.sh (19:45, 21:00): a start after ~21:30 local keeps clear.

Usage:
    python3 survey_products.py --date 2025-01-23
    python3 survey_products.py --date 2025-01-23 --dry-run
    python3 survey_products.py --date 2025-01-23 --steps compare --force
    python3 survey_products.py --date 2025-01-23 --setup-coef 0.034 \\
        --setup-fitted-to 'none:repeat crossings, dem_from_contours.py --fit-setup on 2025-01-18..23'
    python3 survey_products.py --all
    python3 survey_products.py --summary
    (station defaults; override the inputs with --photo-roots, --survey-dirs,
     --adcp, --adcp-navd88, --chatham, --wis, --ndbc, --live-contours, ...)
"""

import io
import contextlib
import os
import re
import sys
import csv
import json
import time
import shlex
import shutil
import hashlib
import argparse
import subprocess
from pathlib import Path
from datetime import datetime, timezone, date as date_cls, timedelta

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Station (NUC) defaults; every one can be overridden on the command line.
DEFAULT_OUTPUT_ROOT = Path("/mnt/I2Rgus_Data/survey_products")
DEFAULT_CONFIG = HERE / "surveys" / "survey_dates.csv"
DEFAULT_SURVEYS = HERE / "surveys" / "surveys.csv"
CHELSEA = Path("/mnt/I2Rgus_Data/Chelsea_calibration")
WATERLINE = Path("/mnt/I2Rgus_Data/waterline")
DEFAULT_PHOTO_ROOTS = ["/mnt/I2Rgus_Data/Chelsea_pics", "/mnt/I2Rgus_Data/waterline/composite",
                       "/home/argus_user/owg_marconi/images",
                       "/mnt/I2Rgus_Data/waterline/archive/images_timex"]
DEFAULT_SURVEY_DIRS = [str(CHELSEA), str(WATERLINE), str(WATERLINE / "calibration"),
                       str(HERE / "calibration"), str(HERE / "surveys")]
DEFAULT_LIVE_CONTOURS = WATERLINE / "contour_points_timex.csv"
LIVE_CRON = HERE / "waterline_timex_cron.sh"

STEPS = ("forcing", "pointing", "detect", "filter", "dem", "maps", "compare")
ERAS = ("adcp", "chatham", "live")
LABELS = ("INDEPENDENT", "CROSS-VALIDATED", "PARTLY-CIRCULAR", "CIRCULAR")   # best -> worst

# Wave setup, as the live station since 6 Oct 2026 (waterline_timex_cron.sh SETUP_COEF).
SETUP_COEF = 0.037
G = 9.81
# Where each setup coefficient in use was fitted. A comparison with that
# survey cannot test the setup (nor the overall level the setup sets). The
# survey is recognised by its CONTENT (sha256), not only its file name: a
# copy saved under another name is still the same shots; and a points
# survey of the same date is treated as possibly the same shots (a
# re-export changes the bytes) unless surveys.csv says otherwise
# (label_override_reason). Any other C must say where it came from
# (--setup-fitted-to), or the run is refused.
SETUP_FITS = {
    0.037: {"survey_stem": "2026-09-29_Marconi_Checkshots",
            "survey_sha256": "49cc01bfa7a41f2ebe34a678167e47376cbd4e07d9b6a62e708647ead7b23beb",
            "survey_type": "points", "survey_date": "2026-09-29",
            "still_water": "GNSS-R (tide_elevation_navd88 of the live rows)",
            "how": "median of the per-frame C = -(waterline - RTK)/sqrt(Hs*L0) of the waterlines "
                   "lying on the 2026-09-29 RTK transects (14 frames, 29 Sep - 5 Oct 2026, +0.8 to "
                   "+1.3 m); repeat crossings (dem_from_contours.py --fit-setup) agree best at "
                   "0.03-0.04 (waterline_timex_cron.sh). The per-frame fit took the RTK elevation where "
                   "the line lay WITHOUT setup. On the 18 frames of 29 Sep - 5 Oct that lie on the transects "
                   "both ways, RTK - waterline is +0.33 m without setup and still +0.13 m with C = 0.037: "
                   "~0.07 m of that because a line given the setup is re-projected landward onto higher "
                   "beach (~0.3 x the setup, by the slope and the cross-shore distance), ~0.06 m because "
                   "C = 0.037 is below these 18 frames' own per-frame C (0.046 with no re-projection; the "
                   "14-frame fit gave 0.037)",
            "wave_currency": "offshore_hs_m = hs_best of archive/waves_marconi.csv (ADCP wh_4061 "
                             "currency), offshore_tp_s = NDBC 44008 peak period",
            # The 29 Sep - 5 Oct 2026 frames on the RTK transects both ways (Oct 2026: survey_products.py
            # --date 2026-09-29 --window 2026-09-29 2026-10-05 on the station's rows, then survey_compare's
            # transect check of each frame without setup (re-projected at still water), with the setup
            # added in place, and with it re-projected): medians of waterline - RTK over the 18 frames
            # (all lines: compare_rtk.py's set), the median setup applied, and per frame
            # (with - without) / setup ('keep_paired': the share of its setup a line keeps once
            # re-projected). The consistency-filtered lines (17 of those frames) give 0.317 / 0.113 m and
            # the same 0.71.
            "rtk_check": {"without_setup_m": 0.333, "with_setup_m": 0.127, "setup_applied_m": 0.257,
                          "with_setup_in_place_m": 0.058, "keep_paired": 0.71, "c_per_frame_in_place": 0.046,
                          "frames": 18, "window": "2026-09-29 .. 2026-10-05"}},
}

# Earlier measurements of a date's REAL photos, with the same calibrations, against its survey,
# made before this script: until the date is built from those photos (on the station computer,
# where they are) they are the only real numbers for it. The README states them, per camera,
# next to what this build should then give. Matched by date, survey and every camera's
# calibration file.
#
# SIGN. Each entry keeps the quantity its tool printed ('reported', 'reported_as', and what a
# positive value means) and the same numbers as waterline - survey ('waterline_minus_survey',
# the sign of DEM - survey), written out so nothing is flipped in the head. pointing_fix_figure.py
# prints 'lidar elevation - water level at the waterline': positive = the lidar beach lies ABOVE
# the level the line was given, i.e. the line reads LOW and waterline - lidar is NEGATIVE.
PRIOR_CHECKS = [
    {"date": "2025-01-23", "survey_stem": "2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm",
     "survey_text": "the Jan 2025 lidar",
     "eo": {"c1": "CACO03_c1_20250123_EO.yaml", "c2": "CACO03_c2_20250123_EO.yaml"},
     "window": "2025-01-18 .. 2025-01-23",
     "tool": "pointing_fix_figure.py",
     # what pointing_fix_figure.py printed on the real Jan 18-23 photos: median and NMAD per camera
     "reported_as": "lidar elevation - water level at the waterline",
     "reported": {"c1": [0.46, 0.19], "c2": [0.15, 0.23]},
     "reported_meaning": "positive there = the lidar beach lies ABOVE the level each line was given: the lines "
                         "read LOW",
     # the same numbers as waterline - lidar (negative = the lines read LOW)
     "waterline_minus_survey": {"c1": [-0.46, 0.19], "c2": [-0.15, 0.23]},
     "setup_coef": 0.0, "predicted_setup_m": 0.28,
     "how": "unfiltered waterlines, NO wave setup (C = 0), the detection search envelope from the earlier "
            "lidar-fitted pointing (so its detections may differ from this build's, whose envelope is placed "
            "with the Mar 2025 lidar)",
     # by this script's own rule (a search envelope placed with the compared survey): that envelope was
     # the Jan 2025 lidar's band projected through the pointing fitted to the Jan 2025 lidar
     "label": "PARTLY-CIRCULAR",
     "label_why": "its search envelope was placed with this same Jan 2025 lidar, through the lidar-fitted "
                  "pointing; the calibrations are the GCP ones and the elevations used no setup"},
]

# Storms the beach is known to have changed in, measured on the station's own data. A product
# whose window reaches into one grids frames of a changing beach: its README says so, with what
# was measured and where.
STORMS = [
    {"name": "the 25-26 Sep 2026 storm", "first": "2026-09-25", "last": "2026-09-26",
     "change": "the station's like-for-like DEM rebuild (17-22 Sep against 1-6 Oct 2026, the same processing: "
               "dem_change.py, Oct 2026) found real, non-uniform change across it: median -0.12 m, p10 -0.50, "
               "p90 +0.32 m over 354 cells (level of detection 0.156 m), the largest loss (to ~-1 m) at "
               "N 4638400-4638470, near the c1/c2 seam",
     "loss_northing": (4638400.0, 4638470.0)},
]

# DEM settings of the live cron (waterline_timex_cron.sh), so a survey-date
# DEM is built like the station's own.
DEM_CELL, DEM_MIN_POINTS, DEM_MAX_SPREAD, DEM_MAX_HS, DEM_MAX_DAY_OFFSET = 2.0, 3, 0.5, 1.5, 0.15

# Pointing check thresholds (pointing_check.py --warn/--fail).
POINT_WARN_DEG, POINT_FAIL_DEG = 0.1, 0.5

# Fewer compared cells or points than this is not an estimate (as survey_compare.py's MIN_BAND_N).
FEW_N = 10

# A product built from synthetic test inputs (the sandbox fixture: photos rendered from the lidar with the
# water level and setup planted) says so on every output: --synthetic NOTE, or a file of this name in a
# photo root.
SYNTHETIC_NOTE = "tests the pipeline, not the beach"
SYNTHETIC_MARKER = "SYNTHETIC_FIXTURE"
RUN = {"synthetic": None}

# Expected run times, for the printed estimates. Detection: detect_original_view.py runs the
# detector on the whole 2448 x 2048 frame, measured at 22-26 s per photo on a 2.8 GHz Xeon in
# Oct 2026 (60 photos in 1358-1505 s; the 3.5 s used before was process_chelsea.py's cropped,
# resampled frames). NOT yet measured on the station's NUC10i3: its per-10% progress lines
# print the real rate; put it here once seen.
SEC_PER_PHOTO_DETECT = 25.0
SEC_PER_DAY_POINTING = 3.0

# Disk a detection needs while it runs: the detector's debug images (~8 MB a photo) and full-size
# overlays (~4 MB), removed at the end of each camera (prune_detector_scratch), plus the products
# (~0.3 GB a two-camera week). A full disk mid-detection leaves half-written files; checked first.
DETECT_MB_PER_PHOTO = 12.0
PRODUCT_MB = 300.0

# The station's own jobs (local time; crontab lines in waterline_timex_cron.sh and owg.sh). A long
# detection on the 2-core NUC should not run across the waterline cron, which reads the same
# disks and needs both cores; this run lowers its own priority (--nice) and says when it would
# overlap them.
STATION_JOBS_LOCAL = [("12:30", "waterline_timex_cron.sh"), ("19:25", "waterline_timex_cron.sh"),
                      ("19:45", "cleanup.sh"), ("20:55", "waterline_timex_cron.sh"), ("21:00", "cleanup.sh")]


def station_job_overlap(seconds, now=None):
    """The station jobs (STATION_JOBS_LOCAL) that start within the next `seconds` of local time.
    -> list of 'HH:MM name'."""
    now = now or datetime.now()
    hits = []
    for day in (0, 1, 2):
        base = (now + timedelta(days=day)).replace(second=0, microsecond=0)
        for hm, name in STATION_JOBS_LOCAL:
            t = base.replace(hour=int(hm[:2]), minute=int(hm[3:]))
            if now <= t <= now + timedelta(seconds=seconds):
                hits.append(f"{hm} {name}")
    return hits


def overlap_note(seconds):
    hits = station_job_overlap(seconds)
    return ("" if not hits else
            f"it would overlap the station's own jobs ({', '.join(hits)} local): it runs at a lower priority "
            f"(--nice), but a start after ~21:30 or before ~06:00 local keeps clear of them (owg.sh also runs "
            f"hourly at :40, 06:40-20:40)")

# The scripts each step runs. Their sha256 is part of the step's signature,
# so a changed script rebuilds the step (another job revises some of these),
# and each stamp records the commit and the sha256s that built its outputs:
# provenance.json reports those, not the code present when it is written.
# Steps that run inside this script (pointing, live forcing and setup, the
# merge, the plan view) carry a version number in their signature instead
# (bumped when their code changes) and record this script's sha256.
STEP_SCRIPTS = {
    "forcing": ["historical_forcing.py"],
    "pointing": ["horizon_check.py", "pointing_check.py", "estimate_eo_rotation.py"],
    "detect": ["detect_original_view.py", "waterline_detector_v5.py", "extract_elevation_contours.py",
               "georectify.py", "view_reproject.py", "compare_dem_survey.py"],
    "detect_live": ["georectify.py", "extract_elevation_contours.py"],
    "merge": [],
    "filter": ["waterline_consistency.py"],
    "dem": ["dem_from_contours.py", "asc_to_geotiff.py", "dem_figure.py"],
    "map": ["daily_elevation_map.py"],
    "map_plan": ["compare_dem_survey.py", "survey_compare.py"],
    "compare": ["survey_compare.py", "compare_dem_survey.py", "compare_rtk.py", "asc_to_geotiff.py"],
}

# Horizontal CRS written into the GeoTIFFs (--geotiff-epsg). The grids are in
# the frame of the calibration's GCPs, assumed to be that of the surveys
# (NAD83(2011) / UTM 19N, EPSG:6348: the GCP files state no CRS); 32619 (WGS 84 /
# UTM 19N) was asked for, so it is the default, and every README says the tag is
# nominal.
GEOTIFF_EPSG = 32619

# Colours (as dem_figure.py / survey_compare.py): neutral ink, data in colour.
INK, INK2, MUTED, GRID = "#1f1e1c", "#52514e", "#898781", "#e6e5df"
SAND = ["#483413", "#674d21", "#876834", "#a7844d", "#c5a36e", "#e0c399", "#f7e5cb"]   # low -> high
SERIES_1, SERIES_2 = "#2a78d6", "#eb6834"      # categorical slots 1 and 2


# ---------------------------------------------------------------------
# Printing, running, stamps
# ---------------------------------------------------------------------

def say(label, value=""):
    print(f"{label:<18}: {value}" if value != "" else label, flush=True)


def warn(text):
    say("WARNING", text)


def rule(text):
    print(f"\n=== {text} " + "=" * max(3, 70 - len(text)), flush=True)


def cmd_text(cmd):
    return " ".join(shlex.quote(str(c)) for c in cmd)


def run_cmd(cmd, log_path, cwd=None):
    """Runs cmd, streaming its output to the console (indented) and to log_path. -> exit code."""
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    say("running", cmd_text(cmd))
    t0 = time.time()
    with open(log_path, "a") as lf:
        lf.write(f"\n$ {cmd_text(cmd)}\n")
        lf.flush()
        p = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True, bufsize=1, cwd=str(cwd or HERE))
        for line in p.stdout:
            lf.write(line)
            print("    | " + line.rstrip(), flush=True)
        rc = p.wait()
        lf.write(f"[exit {rc}, {time.time() - t0:.0f} s]\n")
    say("finished", f"exit {rc} in {time.time() - t0:.0f} s (log {log_path})")
    return rc


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(1 << 20)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def stamp_path(out, key):
    return Path(out) / "logs" / "stamps" / f"{key}.json"


def normal(obj):
    """JSON round trip, so a signature compares equal to the one read back."""
    return json.loads(json.dumps(obj, default=str))


def is_fresh(out, key, sig, outputs, inputs=()):
    """(True, '') when step `key` last ran with this signature, its outputs all exist and
    none of its inputs is newer than them; else (False, why)."""
    st = stamp_path(out, key)
    if not st.exists():
        return False, "not built yet"
    try:
        old = json.loads(st.read_text())
    except (OSError, ValueError):
        return False, "its stamp is unreadable"
    if old.get("sig") != normal(sig):
        return False, "its settings changed"
    missing = [o for o in outputs if not Path(o).exists()]
    if missing:
        return False, f"{Path(missing[0]).name} is missing"
    t_out = min(Path(o).stat().st_mtime for o in outputs) if outputs else 0
    newer = [i for i in inputs if i and Path(i).exists() and Path(i).stat().st_mtime > t_out + 1]
    if newer:
        return False, f"{Path(newer[0]).name} is newer than its outputs"
    return True, ""


def write_stamp(out, key, sig, cmds=(), extra=None, scripts=None):
    """Records a finished step: its signature, commands, the software that ran it
    (commit, sha256 of each script) and a token later steps put in their own signature."""
    st = stamp_path(out, key)
    st.parent.mkdir(parents=True, exist_ok=True)
    d = {"sig": normal(sig), "commands": [cmd_text(c) if not isinstance(c, str) else c for c in cmds],
         "finished_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")}
    sw = software()
    used = list(scripts or []) + ["survey_products.py"]
    d["software"] = {"commit": sw.get("commit"),
                     "scripts": {s: sw["script_sha256"].get(s) or script_sha(s) for s in used},
                     "not_as_committed": [s for s in used if s in sw.get("scripts_not_as_committed", [])]}
    if extra:
        d.update(normal(extra))
    if RUN.get("synthetic"):
        d["synthetic"] = RUN["synthetic"]           # built from synthetic test inputs: every report says so
    d["token"] = hashlib.sha1((json.dumps(d["sig"], sort_keys=True) + d["finished_utc"]).encode()).hexdigest()[:16]
    tmp = Path(str(st) + ".tmp")
    tmp.write_text(json.dumps(d, indent=1) + "\n")
    os.replace(str(tmp), str(st))


def update_stamp(out, key, **kv):
    """Sets fields of a finished step's stamp in place, KEEPING its token (nothing downstream
    rebuilds): for what describes the same outputs (e.g. the folders the same photos were found in)."""
    st = read_stamp(out, key)
    if not st or all(st.get(k) == normal(v) for k, v in kv.items()):
        return
    st.update(normal(kv))
    sp_ = stamp_path(out, key)
    tmp = Path(str(sp_) + ".tmp")
    tmp.write_text(json.dumps(st, indent=1) + "\n")
    os.replace(str(tmp), str(sp_))


def read_stamp(out, key):
    try:
        return json.loads(stamp_path(out, key).read_text())
    except (OSError, ValueError):
        return None


def token(out, key):
    """The stamp token of step `key` (None when it has not finished): a later step puts it in
    its signature, so a rebuilt upstream step makes it stale."""
    d = read_stamp(out, key)
    return (d or {}).get("token")


def drop_stamp(out, key):
    """A step about to rebuild loses its stamp first: if it fails, nothing downstream can look
    current against outputs it may have half rewritten."""
    try:
        stamp_path(out, key).unlink()
    except OSError:
        pass


_SHA_CACHE = {}


def script_sha(name):
    """sha256 (16 hex) of a repository script, or None if absent."""
    if name not in _SHA_CACHE:
        p = HERE / name
        _SHA_CACHE[name] = sha256(p)[:16] if p.exists() else None
    return _SHA_CACHE[name]


def scripts_sig(step):
    return {s: script_sha(s) for s in STEP_SCRIPTS[step]}


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

CONFIG_FIELDS = ["date", "camera", "enabled", "station", "first_day", "last_day", "eo_file",
                 "envelope_survey", "utc_hours", "era", "notes"]
SURVEY_FIELDS = ["date", "name", "path", "survey_type", "survey_date", "label", "why"]


def read_table(path, fields):
    path = Path(path)
    if not path.exists():
        sys.exit(f"configuration not found: {path}")
    with open(path, newline="") as f:
        rd = csv.DictReader(f)
        missing = [k for k in fields if k not in (rd.fieldnames or [])]
        if missing:
            sys.exit(f"{path}: missing column(s) {', '.join(missing)} (expected {', '.join(fields)})")
        rows = []
        for r in rd:
            r = {k: (v or "").strip() for k, v in r.items() if k}
            if not r.get(fields[0]) or r[fields[0]].startswith("#"):
                continue
            rows.append(r)
    return rows


def parse_hours(text):
    """'13.5-18' -> (13.5, 18.0); '' or 'all' -> (0, 24)."""
    t = (text or "").strip().lower()
    if t in ("", "all"):
        return 0.0, 24.0
    try:
        a, b = (float(v) for v in t.split("-"))
    except ValueError:
        sys.exit(f"utc_hours {text!r}: expected e.g. 13.5-18")
    return a, b


def day_add(day, n):
    return (date_cls.fromisoformat(day) + timedelta(days=n)).isoformat()


def days_between(first, last):
    a, b = date_cls.fromisoformat(first), date_cls.fromisoformat(last)
    return [(a + timedelta(days=i)).isoformat() for i in range((b - a).days + 1)]


def resolve_file(text, dirs):
    """A path as given, else the first of `dirs` holding that file name; globs allowed.
    -> list of Paths (empty when nothing matches)."""
    if not text:
        return []
    p = Path(text)
    if any(ch in text for ch in "*?["):
        hits = sorted(Path(p.parent).glob(p.name)) if p.is_absolute() else []
        for d in dirs:
            if not hits:
                hits = sorted(Path(d).glob(text))
        return hits
    if p.exists():
        return [p]
    for d in dirs:
        q = Path(d) / p.name
        if q.exists():
            return [q]
    return []


def build_plan(date, cfg_rows, survey_rows, args):
    """Everything the steps need for one date, from the tables and the options."""
    rows = [r for r in cfg_rows if r["date"] == date]
    if not rows:
        known = sorted({r["date"] for r in cfg_rows})
        sys.exit(f"{date} is not in {args.config} (dates: {', '.join(known)})")
    enabled = [r for r in rows if r["enabled"] in ("1", "yes", "true", "True")]
    if not enabled and not args.force_disabled:
        notes = "; ".join(sorted({r["notes"] for r in rows if r["notes"]}))
        print(f"\n{date} is DISABLED in {args.config}: {notes or '(no reason given)'}\n"
              f"Not built. To build it anyway: --date {date} --force-disabled", flush=True)
        return None
    use = rows if args.force_disabled else enabled
    if args.force_disabled and len(enabled) < len(rows):
        warn(f"--force-disabled: building {date} with its disabled row(s) "
             f"({', '.join(r['camera'] for r in rows if r not in enabled)})")
    eras = {r["era"] for r in use}
    if len(eras) != 1 or not eras <= set(ERAS):
        sys.exit(f"{date}: era must be one of {ERAS} and the same for every camera (got {sorted(eras)})")
    era = eras.pop()
    cams = {}
    for r in sorted(use, key=lambda r: r["camera"]):
        cam = r["camera"]
        if cam not in ("c1", "c2"):
            sys.exit(f"{date}: camera {cam!r} (expected c1 or c2)")
        first, last = (args.window or (r["first_day"], r["last_day"]))
        hours_text = args.utc_hours if args.utc_hours is not None else r["utc_hours"]
        eo = resolve_file(r["eo_file"], [args.calibration])
        env = resolve_file(r["envelope_survey"], args.survey_dirs) if r["envelope_survey"] else []
        hours = parse_hours(hours_text)
        th = parse_hours(r["utc_hours"])
        cams[cam] = {"station": r["station"], "first": first, "last": last,
                     "hours": hours, "hours_text": f"{hours[0]:g}-{hours[1]:g}",
                     # the table's own window and hours (what no --window / --utc-hours means)
                     "table_first": r["first_day"], "table_last": r["last_day"],
                     "table_hours_text": f"{th[0]:g}-{th[1]:g}",
                     "eo_file": r["eo_file"], "eo_path": eo[0] if eo else None,
                     "io_path": Path(args.calibration) / f"CACO05_{cam}_20240801_IO.yaml",
                     "envelope_survey": r["envelope_survey"], "envelope_path": env[0] if env else None,
                     "notes": r["notes"], "enabled": r in enabled}
    surveys = [s for s in survey_rows if s["date"] == date]
    plan = {"date": date, "era": era, "cams": cams,
            "first": min(c["first"] for c in cams.values()),
            "last": max(c["last"] for c in cams.values()),
            "surveys": surveys, "out": Path(args.output_root) / date,
            "disabled_cams": {r["camera"]: r["notes"] for r in rows if r not in use},
            "forced_disabled": bool(args.force_disabled and len(enabled) < len(rows))}
    return plan


# ---------------------------------------------------------------------
# Photos
# ---------------------------------------------------------------------

def epoch_of(name):
    m = re.match(r"^(\d{9,11})\.", name)
    return int(m.group(1)) if m else None


def utc_of(name):
    e = epoch_of(name)
    return datetime.fromtimestamp(e, tz=timezone.utc) if e is not None else None


_SCAN = {}


def scan_photos(plan, cam, args, skip_days=()):
    """The camera's photos in the window (detect_original_view.collect_photos: each name once, only
    copies of the lens file's frame size). Cached per run: the photo folders are large."""
    from detect_original_view import collect_photos
    from georectify import load_intrinsics
    c = plan["cams"][cam]
    size = None
    if Path(c["io_path"]).exists():
        io_ = load_intrinsics(c["io_path"])
        size = (int(io_[0]), int(io_[1]))
    key = (cam, c["first"], c["last"], tuple(c["hours"]), c["station"], tuple(args.photo_roots),
           tuple(sorted(skip_days or ())), size)
    if key not in _SCAN:
        _SCAN[key] = collect_photos(args.photo_roots, cam, c["first"], c["last"], c["hours"],
                                    station=c["station"] or None, skip_days=skip_days, frame_size=size)
        cnt = _SCAN[key][1]
        if cnt.get("wrong_size") or cnt.get("unreadable"):
            warn(f"{cam}: photo copies rejected: {cnt['wrong_size']} of another size than the lens file's "
                 f"{size or '?'} (reduced copies), {cnt['unreadable']} unreadable (broken links?)"
                 + (f"; {cnt['no_good_copy']} photo(s) with no good copy left out" if cnt.get("no_good_copy")
                    else ""))
    return _SCAN[key]


def per_day(photos):
    out = {}
    for p in photos:
        t = utc_of(p.name)
        out.setdefault(t.strftime("%Y-%m-%d"), []).append(p)
    return out


def link_photos(paths, folder):
    """Hard link (else symbolic link, else copy) each photo into folder."""
    from detect_original_view import place
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    n = 0
    for p in paths:
        dst = folder / Path(p).name
        if not dst.exists():
            place(Path(p), dst)
            n += 1
    return n


# ---------------------------------------------------------------------
# Step: forcing
# ---------------------------------------------------------------------

FORCING_INPUTS = ("adcp", "adcp_navd88", "chatham", "wis", "ndbc")


def forcing_cmd(plan, args):
    cmd = [sys.executable, HERE / "historical_forcing.py", "--start", plan["first"], "--end", plan["last"],
           "--output-dir", plan["out"] / "forcing", "--setup-coef", args.setup_coef]
    for opt in FORCING_INPUTS:
        v = getattr(args, opt)
        if v:
            cmd += ["--" + opt.replace("_", "-"), v]
    if getattr(args, "gnssr_spline", None):
        cmd += ["--gnssr-spline", args.gnssr_spline]
    if getattr(args, "gauge_csv", None):
        cmd += ["--gauge-archive", args.gauge_csv]
    if args.no_download:
        cmd.append("--no-download")
    return cmd


def forcing_inputs(args):
    """{key: Path or None}: each input as historical_forcing.py will find it (the file given, else
    the first station folder holding its default name), so a file that arrives in one of those
    folders, or is replaced, rebuilds the forcing."""
    from historical_forcing import find_input, ForcingError
    found = {}
    for k in FORCING_INPUTS:
        try:
            found[k] = find_input(k, getattr(args, k))
        except ForcingError as exc:
            raise StepFailed(str(exc))
    return found


def file_sig(p):
    """[path, size, mtime] of a file, or None."""
    if not p or not Path(p).exists():
        return None
    st = Path(p).stat()
    return [str(p), st.st_size, int(st.st_mtime)]


def step_forcing(plan, args, state):
    out = plan["out"]
    fdir = out / "forcing"
    outputs = [fdir / "water_level.csv", fdir / "waves.csv", fdir / "forcing.json"]
    if plan["era"] == "live":
        return live_forcing(plan, args, state)
    cmd = forcing_cmd(plan, args)
    found = forcing_inputs(args)
    # The 2026 record the share of setup GNSS-R sees is measured from (README, provenance
    # setup.still_water_reference). NOT in the signature: archive/waves_marconi.csv grows every
    # day with the cron, and a new forcing token would redo an hour of detection for a number
    # only the README reports. A file that arrives or changes refreshes that number in place.
    share_sig = setup_share_sig(args)
    sig = {"cmd": cmd_text(cmd), "inputs": {k: file_sig(p) for k, p in found.items()},
           "scripts": scripts_sig("forcing")}
    fresh, why = is_fresh(out, "forcing", sig, outputs, list(found.values()))
    if fresh and not args.force:
        if (read_stamp(out, "forcing") or {}).get("setup_share_inputs") != share_sig:
            return refresh_setup_share(out, args, share_sig)
        say("forcing", f"SKIPPED: {fdir} is up to date (--force to rebuild)")
        return "skipped"
    say("forcing", f"{plan['first']} .. {plan['last']} ({why}); expected 10-30 s on the NUC")
    for k, p in found.items():
        say(f"  {k}", str(p) if p else "not found (historical_forcing.py looks in the station folders)")
    drop_stamp(out, "forcing")
    rc = run_cmd(cmd, out / "logs" / "forcing.log")
    if rc != 0 or not all(o.exists() for o in outputs):
        raise StepFailed(f"historical_forcing.py failed (exit {rc}); see {out / 'logs' / 'forcing.log'}"
                         + ("" if rc != 2 else " (its ERROR line says why: e.g. too few wave hours for the setup)"))
    n_waves = max(0, sum(1 for _ in open(fdir / "waves.csv")) - 1)
    if args.setup_coef and not n_waves:
        # never an hour of detection that then drops every frame for want of a setup
        raise StepFailed(f"{fdir / 'waves.csv'} has no wave hours, and the setup (C = {args.setup_coef}) needs "
                         f"them: give the wave files (--wis, --ndbc, --adcp) or --setup-coef 0")
    write_stamp(out, "forcing", sig, [cmd], scripts=STEP_SCRIPTS["forcing"],
                extra={"setup_share_inputs": share_sig})
    return "built"


SHARE_KEYS = ("quantified", "note", "share", "k", "k_se", "c", "n", "r", "typical_x", "typical_m", "how", "inputs")


def setup_share_sig(args):
    """[[path, size, mtime], ...] of the three files historical_forcing.gnssr_setup_share() reads."""
    from historical_forcing import setup_share_inputs
    share = setup_share_inputs(getattr(args, "gnssr_spline", None), getattr(args, "gauge_csv", None))
    return normal([[p, file_sig(p)] for p in share])


def refresh_setup_share(out, args, share_sig):
    """The forcing is current but the files the GNSS-R setup share is measured from changed (or
    arrived): measure it again into forcing.json and the stamp, KEEPING the stamp's token, so
    nothing downstream rebuilds (water level and waves are unchanged). forcing_report.txt keeps
    the old text; the README and provenance read forcing.json."""
    from historical_forcing import gnssr_setup_share
    fj = out / "forcing" / "forcing.json"
    d = load_json(fj) or {}
    wv = d.get("waves") or {}
    coef = wv.get("setup_coef")
    swr = dict(wv.get("still_water_reference") or {})
    if coef:
        for k in SHARE_KEYS:
            swr.pop(k, None)
        swr.update(gnssr_setup_share(coef, getattr(args, "gnssr_spline", None), getattr(args, "gauge_csv", None)))
        wv["still_water_reference"] = swr
        d["waves"] = wv
        tmp = Path(str(fj) + ".tmp")
        tmp.write_text(json.dumps(d, indent=1) + "\n")
        os.replace(str(tmp), str(fj))
    st = read_stamp(out, "forcing") or {}
    st["setup_share_inputs"] = share_sig
    sp = stamp_path(out, "forcing")
    tmp = Path(str(sp) + ".tmp")
    tmp.write_text(json.dumps(st, indent=1) + "\n")
    os.replace(str(tmp), str(sp))
    say("forcing", "up to date; GNSS-R setup share measured again (its 2026 input files changed): "
        + (f"GNSS-R sees {swr['share']:.2f} of the setup" if swr.get("quantified") and swr.get("share") is not None
           else (swr.get("note") or "not measured")))
    return "skipped (up to date; setup share refreshed)"


class StepFailed(Exception):
    pass


def live_rows_file(plan, args, state):
    """The live archive's rows for the window (and the enabled cameras and hours), unmodified,
    in waterlines/archive_rows.csv: read once, used by the forcing and detect steps."""
    out = plan["out"]
    dst = out / "waterlines" / "archive_rows.csv"
    src = Path(args.live_contours)
    if not src.exists():
        raise StepFailed(f"live contour file not found: {src} (--live-contours)")
    st = src.stat()
    sig = {"src": str(src), "size": st.st_size, "mtime": int(st.st_mtime), "v": 2,
           "cams": {c: [v["first"], v["last"], v["hours"], v["station"]] for c, v in plan["cams"].items()}}
    fresh, why = is_fresh(out, "archive_rows", sig, [dst])
    # read at most once per run: the forcing and detect steps both ask for the rows, and with --force
    # a second read would give the detection another archive_rows build than the forcing's (the
    # outputs would then never be one build)
    if fresh and (not args.force or state.get("archive_rows_built")):
        return dst
    say("archive rows", f"reading {src} ({st.st_size / 1e6:.0f} MB; ~1 min per GB on the NUC) for "
        f"{plan['first']} .. {plan['last']}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    # Read into a new file first: the rows of the earlier build (and their stamp, which says which
    # window built the outputs on disk) are replaced only once the new window has rows. A window with
    # none, or a read cut short, leaves them, so the README and the 'rebuild everything' line keep
    # naming the window the DEM and waterlines on disk came from.
    new = dst.parent / (dst.name + ".new")
    # The cron rewrites this file in place (opened with 'w') at 12:30, 19:25 and 20:55 local:
    # a read that overlaps a rewrite sees a short or cut file. Size and time are compared
    # before and after; on a change the file is read again, once.
    try:
        for attempt in (1, 2):
            before = (src.stat().st_size, src.stat().st_mtime)
            t0 = time.time()
            res = read_live_rows(src, new, plan)
            after = (src.stat().st_size, src.stat().st_mtime)
            if before == after:
                break
            if attempt == 1:
                warn(f"{src} changed while it was read (the cron rewrites it): reading it again")
            else:
                raise StepFailed(f"{src} changed twice while it was read (the cron is rewriting it): "
                                 f"run again in a few minutes")
        n_in, n_out, frames, other, bad = res
        say("archive rows", f"{n_out} of {n_in} rows, {len(frames)} frames ({time.time() - t0:.0f} s)")
        if other:
            warn("rows of frames named for another station than the date's left out: "
                 + ", ".join(f"{c} {sum(v.values())} rows ({', '.join(f'{k} {n}' for k, n in sorted(v.items()))})"
                             for c, v in sorted(other.items())))
        if bad:
            warn(f"{bad} row(s) with an unreadable capture time skipped (a cut-off line?)")
        if not n_out:
            kept = read_stamp(out, "archive_rows")
            raise StepFailed(f"no rows in {src} for {plan['first']} .. {plan['last']}"
                             + (f"; the archive rows of the earlier build ({kept.get('rows_out')} rows) and the "
                                f"outputs made from them are kept" if kept and dst.exists() else ""))
        drop_stamp(out, "archive_rows")
        os.replace(str(new), str(dst))
    finally:
        if new.exists():
            try:
                new.unlink()
            except OSError:
                pass
    say("archive rows", f"-> {dst}")
    write_stamp(out, "archive_rows", sig, extra={"rows_in": n_in, "rows_out": n_out, "frames": len(frames),
                                                 "other_station_rows": other, "malformed_rows": bad})
    state["archive_rows_built"] = True
    return dst


def read_live_rows(src, dst, plan):
    """Copies the rows of the window, cameras, hours and station from src to dst.
    -> (rows read, rows kept, frames kept, {cam: {station: rows left out}}, malformed rows)."""
    from detect_original_view import station_of
    n_in = n_out = bad = 0
    frames, other = set(), {}
    with open(src, newline="") as fi, open(str(dst) + ".tmp", "w", newline="") as fo:
        rd = csv.reader(fi)
        head = next(rd, None)
        if not head:
            raise StepFailed(f"{src} is empty (being rewritten by the cron?)")
        ix = {c: i for i, c in enumerate(head)}
        for k in ("source_file", "camera", "capture_time_utc", "tide_elevation_navd88"):
            if k not in ix:
                raise StepFailed(f"{src}: no {k} column")
        wr = csv.writer(fo)
        wr.writerow(head)
        ic, it, isrc = ix["camera"], ix["capture_time_utc"], ix["source_file"]
        win = {c: (v["first"], v["last"], v["hours"], v["station"]) for c, v in plan["cams"].items()}
        need = max(ic, it, isrc)
        for r in rd:
            n_in += 1
            if len(r) <= need:
                bad += 1
                continue
            w = win.get(r[ic])
            if w is None:
                continue
            ts = r[it]
            if not (w[0] <= ts[:10] <= w[1]):
                continue
            try:
                hh = int(ts[11:13]) + int(ts[14:16]) / 60.0
            except ValueError:
                bad += 1
                continue
            if not (w[2][0] <= hh <= w[2][1]):
                continue
            # the station ID in the frame name: a window spanning a re-set must not mix two
            # pointings (the historical eras guard this with detect_original_view.py --station)
            stn = station_of(r[isrc] + ".jpg") if not r[isrc].endswith(".jpg") else station_of(r[isrc])
            if w[3] and stn != w[3]:
                o = other.setdefault(r[ic], {})
                o[stn or "?"] = o.get(stn or "?", 0) + 1
                continue
            wr.writerow(r)
            n_out += 1
            frames.add(r[isrc])
    os.replace(str(dst) + ".tmp", dst)
    return n_in, n_out, frames, other, bad


def cron_settings(path=LIVE_CRON):
    """The live cron's wave and setup settings (VAR=value lines), as text."""
    keep = ("WAVE_BUOY", "USE_MARCONI_WAVES", "SETUP_COEF", "DEM_MAX_HS", "DEM_CELL", "DEM_MIN_POINTS",
            "DEM_MAX_SPREAD", "DEM_MAX_DAY_OFFSET")
    out = {}
    try:
        for line in Path(path).read_text().splitlines():
            m = re.match(r"^\s*([A-Z_]+)=\"?([^\"#]*)\"?", line)
            if m and m.group(1) in keep:
                out[m.group(1)] = m.group(2).strip()
    except OSError:
        pass
    return out


def live_forcing(plan, args, state):
    """Per-frame water level and waves from the live archive rows (GNSS-R and the cron's waves)."""
    out = plan["out"]
    fdir = out / "forcing"
    outputs = [fdir / "water_level.csv", fdir / "waves.csv", fdir / "forcing.json"]
    rows_file = live_rows_file(plan, args, state)
    sig = {"rows": str(rows_file), "after": {"archive_rows": token(out, "archive_rows")},
           "spline": file_sig(args.gnssr_spline) or str(args.gnssr_spline), "setup_coef": args.setup_coef, "v": 4}
    fresh, why = is_fresh(out, "forcing", sig, outputs, [rows_file])
    if fresh and not args.force:
        say("forcing", f"SKIPPED: {fdir} is up to date (--force to rebuild)")
        return "skipped"
    say("forcing", f"live era: the archive rows' own GNSS-R levels and waves ({why})")
    drop_stamp(out, "forcing")
    frames = {}
    with open(rows_file, newline="") as f:
        for r in csv.DictReader(f):
            k = r["source_file"]
            if k in frames:
                continue
            frames[k] = (int(float(r["capture_epoch"])) if r.get("capture_epoch") else epoch_of(k),
                         r["tide_elevation_navd88"], r.get("water_level_source", ""),
                         r.get("offshore_hs_m", ""), r.get("offshore_tp_s", ""), r["camera"])
    fdir.mkdir(parents=True, exist_ok=True)
    items = sorted(frames.items(), key=lambda kv: kv[1][0])
    ep = np.array([v[0] for _, v in items], float)
    lv = np.array([float(v[1]) for _, v in items])
    with open(fdir / "water_level.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time", "water_level_navd88", "source", "sigma_m", "source_file"])
        for k, v in items:
            w.writerow([iso_utc(v[0]), v[1], v[2] or "live archive", "", k])
    hs_ok = [(v[0], float(v[3]), float(v[4])) for _, v in items if v[3] != "" and v[4] != ""]
    with open(fdir / "waves.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_utc", "epoch", "wvht_m", "dpd_s", "source", "hs_sigma_m", "source_file"])
        for k, v in items:
            w.writerow([iso_utc(v[0]), v[0], v[3], v[4], "live archive (offshore_hs_m, offshore_tp_s)", "", k])
    cron = cron_settings(args.live_cron)
    srcs = {}
    for _, v in items:
        srcs[v[2] or "?"] = srcs.get(v[2] or "?", 0) + 1
    # The archive is rebuilt every cron run from the CURRENT spline: check it still agrees.
    spline_check = None
    if args.gnssr_spline and Path(args.gnssr_spline).exists():
        try:
            from marconi_water_level import WaterLevel
            wl = WaterLevel(spline=str(args.gnssr_spline), gauge_csv=args.gauge_csv)
            now, src_now = wl.at(ep)
            d = now - lv
            ok = np.isfinite(d)
            spline_check = {"frames": int(ok.sum()), "median_m": round(float(np.median(d[ok])), 4) if ok.any() else None,
                            "max_abs_m": round(float(np.max(np.abs(d[ok]))), 4) if ok.any() else None,
                            "sources_now": {s: int((src_now == s).sum()) for s in set(src_now.tolist())},
                            "describe": wl.describe()}
            say("GNSS-R check", f"archive level minus today's spline at {ok.sum()} frames: median "
                f"{spline_check['median_m']} m, max |d| {spline_check['max_abs_m']} m")
        except Exception as exc:                       # a cross-check, never a blocker
            spline_check = {"error": str(exc)}
            warn(f"GNSS-R cross-check failed: {exc}")
    else:
        say("GNSS-R check", f"spline {args.gnssr_spline} not here: the archive rows' levels are used as they are")
    setups = [setup_m(h, t, args.setup_coef) for _, h, t in hs_ok]
    use_marconi = cron.get("USE_MARCONI_WAVES") == "1"
    info = {"period": {"start": plan["first"], "end": plan["last"]},
            "water_level": {"method": "GNSS-R measured at the station (gnssrefl spline on NAVD88, QC'd "
                                      "against the Chatham gauge by the cron: gnssr_qc.py)",
                            "sources_used": srcs, "frames": len(items),
                            "cv_rms_m": None, "extrapolation_rms_m": None,
                            "inputs": {"live_contours": str(args.live_contours),
                                       "gnssr_spline": str(args.gnssr_spline)},
                            "spline_check": spline_check,
                            "note": "the levels are those the station's contour file carries; it is rebuilt "
                                    "from the whole archive at every cron run with the spline of that run"},
            "waves": {"method": "the live archive's per-frame offshore_hs_m / offshore_tp_s",
                      "cron_settings": cron,
                      "source": ("archive/waves_marconi.csv (USE_MARCONI_WAVES=1): hs_best = camera wave "
                                 "models trained on the ADCP wh_4061, blended with buoys converted to it "
                                 "(buoy_transfer.py); tp_s = NDBC 44008 peak period"
                                 if use_marconi else
                                 f"NDBC {cron.get('WAVE_BUOY', '44008')} raw (USE_MARCONI_WAVES=0)"),
                      "hs_currency": ("the setup coefficient's own (C was fitted on these same live rows)"
                                      if setup_fit_info(args.setup_coef, args).get("kind") == "known" else
                                      "the live rows' (offshore_hs_m, offshore_tp_s); the currency this C was "
                                      "fitted in is not known here"),
                      "frames_with_waves": len(hs_ok), "frames_without_waves": len(items) - len(hs_ok),
                      "setup_coef": args.setup_coef,
                      "setup_in_window": ({"median_m": round(float(np.median(setups)), 3),
                                           "min_m": round(float(np.min(setups)), 3),
                                           "max_m": round(float(np.max(setups)), 3)} if setups else None)},
            "software": {"script": "survey_products.py (live era)"}}
    (fdir / "forcing.json").write_text(json.dumps(info, indent=1) + "\n")
    lines = [f"Live-era forcing for {plan['date']} ({plan['first']} .. {plan['last']})", "",
             "Water level: the GNSS-R level each archive row carries (water_level_source):",
             "  " + ", ".join(f"{k} {v} frames" for k, v in srcs.items()),
             f"  cross-check against today's spline: {spline_check}",
             "", f"Waves: {info['waves']['source']}",
             f"  {len(hs_ok)} of {len(items)} frames have Hs and Tp; setup C = {args.setup_coef}: "
             + (f"median {np.median(setups):.3f} m ({np.min(setups):.3f}-{np.max(setups):.3f})" if setups else "none"),
             f"  cron settings: {cron}"]
    (fdir / "forcing_report.txt").write_text("\n".join(lines) + "\n")
    try:
        live_forcing_figure(fdir / "forcing.png", ep, lv, hs_ok, args.setup_coef, plan)
    except Exception as exc:
        warn(f"forcing figure failed: {exc}")
    say("forcing", f"{len(items)} frames: water level {lv.min():+.2f} .. {lv.max():+.2f} m NAVD88; "
        f"{len(hs_ok)} with waves -> {fdir}")
    write_stamp(out, "forcing", sig, ["survey_products.py, in-process (live era): per-frame water level and "
                                      "waves from waterlines/archive_rows.csv"], scripts=["marconi_water_level.py"])
    return "built"


def iso_utc(ep):
    return datetime.fromtimestamp(float(ep), tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def setup_m(hs, tp, coef):
    """Wave setup C*sqrt(Hs*L0), L0 = g*Tp^2/(2*pi): extract_elevation_contours.py --setup-coef."""
    return coef * np.sqrt(float(hs) * 9.81 * float(tp) ** 2 / (2 * np.pi))


def live_forcing_figure(path, ep, lv, hs_ok, coef, plan):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    t = [datetime.fromtimestamp(e, tz=timezone.utc) for e in ep]
    fig, axes = plt.subplots(2, 1, figsize=(9, 5.4), sharex=True)
    ax = axes[0]
    ax.plot(t, lv, "o", ms=2.5, color=SERIES_1)
    ax.set_ylabel("water level (m NAVD88)", color=INK2)
    ax.set_title(f"{plan['date']}: still-water level at each frame (GNSS-R, live archive)",
                 loc="left", fontsize=10, color=INK)
    ax = axes[1]
    if hs_ok:
        ts = [datetime.fromtimestamp(e, tz=timezone.utc) for e, _, _ in hs_ok]
        ax.plot(ts, [setup_m(h, p, coef) for _, h, p in hs_ok], "o", ms=2.5, color=SERIES_1)
    ax.set_ylabel(f"setup (m), C = {coef}", color=INK2)
    ax.set_title("wave setup added to each frame", loc="left", fontsize=10, color=INK)
    for a in axes:
        a.grid(color=GRID, lw=0.6)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    loc = mdates.AutoDateLocator()          # hours within a day, days across a week: no repeated labels
    # the whole window, so days without frames show as gaps
    axes[1].set_xlim(datetime.fromisoformat(plan["first"]).replace(tzinfo=timezone.utc),
                     datetime.fromisoformat(plan["last"]).replace(tzinfo=timezone.utc) + timedelta(days=1))
    axes[1].xaxis.set_major_locator(loc)
    axes[1].xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
    axes[1].set_xlabel("UTC", color=INK2)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


# ---------------------------------------------------------------------
# Step: pointing
# ---------------------------------------------------------------------

def eo_date(eo_file):
    m = re.search(r"_(\d{8})_", Path(str(eo_file)).name)
    return f"{m.group(1)[:4]}-{m.group(1)[4:6]}-{m.group(1)[6:]}" if m else None


def survey_day(plan):
    days = [s["survey_date"] for s in plan["surveys"] if s.get("survey_date")]
    return days[0] if days else plan["date"]


def check_pointing(plan, cam, photos, args):
    """
    Per day: the sea horizon against the calibration (tilt and roll; horizon_check.py's fit)
    and the rotation from a reference frame (pan too; pointing_check.py's matching).
    -> (rows, info).
    """
    import cv2
    from georectify import load_intrinsics, load_extrinsics
    from estimate_eo_rotation import horizon_rows, observed_horizon, land_mask
    from horizon_check import fit_tilt_roll, MOVED_PX, FIT_OK_PX, RUN_TOL_DEG
    import pointing_check as pc

    c = plan["cams"][cam]
    io, eo = load_intrinsics(c["io_path"]), load_extrinsics(c["eo_path"])
    days = per_day(photos)
    for d in days:
        days[d].sort(key=lambda p: abs(utc_of(p.name).hour + utc_of(p.name).minute / 60 - 17))
    cal_day = eo_date(c["eo_file"])
    if cal_day in days:
        ref_day, ref_kind = cal_day, "the calibration day: pan, tilt and roll are against the calibration"
    else:
        target = date_cls.fromisoformat(survey_day(plan))
        ref_day = min(days, key=lambda d: abs((date_cls.fromisoformat(d) - target).days))
        ref_kind = (f"the day nearest the survey ({cal_day or 'calibration day'} has no photos here): "
                    f"changes are against that day, not the calibration")
    ref_src = days[ref_day][0]
    rdir = plan["out"] / "pointing" / f"reference_{cam}"
    shutil.rmtree(rdir, ignore_errors=True)
    rdir.mkdir(parents=True)
    ref_img = rdir / ref_src.name
    shutil.copy2(ref_src, ref_img)
    ref_img.with_suffix(".json").write_text(json.dumps({"eo": eo.tolist(), "eo_file": c["eo_file"]}))
    img0 = cv2.imread(str(ref_img))
    if img0 is None:
        raise StepFailed(f"cannot read {ref_src}")
    cv2.imwrite(str(ref_img.with_name(ref_img.stem + "_mask.png")), land_mask(img0.shape[:2], io, eo))
    ref = pc.Reference(ref_img)
    cols = np.arange(150, int(io[0]) - 150, 100)
    pred = horizon_rows(io, eo, cols)
    say(f"pointing {cam}", f"{len(days)} day(s); reference {ref_src.name} ({ref_kind}); "
        f"~{SEC_PER_DAY_POINTING:.0f} s per day on the NUC")
    rows = []
    for day in sorted(days):
        cands = days[day][:args.pointing_frames]
        offs, fits = [], []
        for p in cands:
            g = cv2.imread(str(p), cv2.IMREAD_GRAYSCALE)
            if g is None:
                continue
            obs = observed_horizon(g, cols, pred, window=250)
            if np.isfinite(obs).sum() < 3:
                continue
            (dt, dr), res = fit_tilt_roll(io, eo, cols, obs)
            offs.append(float(np.nanmedian(np.abs(obs - pred))))
            fits.append((float(dt), float(dr), float(np.median(res))))
        r = {"camera": cam, "date": day, "photos": len(days[day]), "horizon_frames": len(offs),
             "horizon_offset_px": round(float(np.median(offs)), 1) if offs else "",
             "horizon_dtilt_deg": round(float(np.median([f[0] for f in fits])), 3) if fits else "",
             "horizon_droll_deg": round(float(np.median([f[1] for f in fits])), 3) if fits else "",
             "horizon_fit_px": round(float(np.median([f[2] for f in fits])), 1) if fits else ""}
        r["horizon"] = ("no horizon" if not fits else
                        "unclear" if r["horizon_fit_px"] > FIT_OK_PX else
                        "off the calibration" if r["horizon_offset_px"] > MOVED_PX else "matches")
        probe = next((p for p in cands if p.name != ref_src.name), None)
        if probe is None and len(days[day]) > 1:
            probe = days[day][1]
        m = None
        if probe is not None:
            img = cv2.imread(str(probe))
            m = pc.measure(img, [ref], io) if img is not None else None
        if m is None:
            r.update(feature_frame=probe.name if probe else "", inliers="", d_azimuth_deg="", d_tilt_deg="",
                     d_roll_deg="", features="not measured" if probe is None else "unreadable")
        elif m.get("status") == "unmatched":
            r.update(feature_frame=probe.name, inliers=m.get("inliers", 0), d_azimuth_deg="", d_tilt_deg="",
                     d_roll_deg="", features="unmatched")             # fog, rain, night, few features
        else:
            big = max(abs(m["d_azimuth"]), abs(m["d_tilt"]), abs(m["d_roll"]))
            r.update(feature_frame=probe.name, inliers=m["inliers"], d_azimuth_deg=round(m["d_azimuth"], 3),
                     d_tilt_deg=round(m["d_tilt"], 3), d_roll_deg=round(m["d_roll"], 3),
                     features="MOVED" if big > args.pointing_fail else "moved" if big > POINT_WARN_DEG else "ok")
        rows.append(r)
    # Day-to-day horizon changes: a day whose tilt or roll differs from the window's
    # median by more than horizon_check's run tolerance has another pointing.
    clear = [r for r in rows if r["horizon"] in ("matches", "off the calibration")]
    med_t = float(np.median([r["horizon_dtilt_deg"] for r in clear])) if clear else 0.0
    med_r = float(np.median([r["horizon_droll_deg"] for r in clear])) if clear else 0.0
    for r in rows:
        horizon_jump = (r in clear and (abs(r["horizon_dtilt_deg"] - med_t) > RUN_TOL_DEG
                                        or abs(r["horizon_droll_deg"] - med_r) > RUN_TOL_DEG))
        if r["features"] == "MOVED" or horizon_jump:
            r["verdict"] = "DIFFERENT POINTING"
            r["why"] = "; ".join(x for x in (
                f"rotation from the reference {r['d_azimuth_deg']:+.2f}/{r['d_tilt_deg']:+.2f}/{r['d_roll_deg']:+.2f} "
                f"deg (pan/tilt/roll) > {args.pointing_fail}" if r["features"] == "MOVED" else "",
                f"horizon tilt/roll {r['horizon_dtilt_deg']:+.2f}/{r['horizon_droll_deg']:+.2f} deg vs the "
                f"window's {med_t:+.2f}/{med_r:+.2f} (> {RUN_TOL_DEG})" if horizon_jump else "") if x)
        elif r["features"] == "moved":
            r["verdict"], r["why"] = "small change", f"rotation {POINT_WARN_DEG}-{args.pointing_fail} deg from the reference"
        elif r["features"] in ("unmatched", "not measured", "unreadable"):
            if r["horizon"] in ("unclear", "no horizon"):
                r["verdict"], r["why"] = "unchecked", "neither the features nor the horizon could be measured"
            else:
                r["verdict"], r["why"] = "horizon only", ("tilt and roll agree with the other days; the pan "
                                                          "could not be measured (features unmatched)")
        else:
            r["verdict"], r["why"] = "matches", ""
    info = {"reference_frame": ref_src.name, "reference_day": ref_day, "reference_kind": ref_kind,
            "calibration": c["eo_file"], "window_horizon_dtilt_deg": round(med_t, 3),
            "window_horizon_droll_deg": round(med_r, 3), "clear_horizon_days": len(clear)}
    # where the calibration puts the horizon: at the top edge of the frame (or off it) the lens model
    # is extrapolated far beyond its GCPs, and a constant offset there says little about the pointing
    fin = np.isfinite(pred)
    if fin.any():
        inside = fin & (pred >= 0) & (pred < float(io[1]))
        info.update(horizon_pred_rows=[round(float(np.min(pred[fin]))), round(float(np.max(pred[fin])))],
                    horizon_cols_in_frame=[int(inside.sum()), int(len(cols))], frame_rows=int(io[1]))
    offs_all = [r["horizon_offset_px"] for r in clear]
    info["window_horizon_offset_px"] = round(float(np.median(offs_all)), 1) if offs_all else None
    off = pointing_offset(info, c["eo_file"])
    if off:
        info["calibration_offset_deg"], info["calibration_offset"] = [off[0], off[1]], off[2]
    return rows, info


# A constant tilt or roll of the whole window against the calibration larger than this (deg) is
# reported with its DEM sensitivity (0.05 deg of tilt is ~0.07 m of DEM at 250-350 m). It is
# never used to leave days out: horizon_check.py's MOVED_PX and RUN_TOL_DEG judge a CHANGE of
# pointing from day to day; this judges what a constant offset may do to the DEM.
# The predicted horizon (estimate_eo_rotation.horizon_rows) carries the full Earth-curvature dip
# and standard refraction (k = 0.13): a perfectly calibrated camera reads 0.000 deg, and
# +/-0.010 deg for k from 0 to 0.25 (true horizon projected with each station EO of 2025-26 and
# fitted with horizon_check.fit_tilt_roll), so 0.05 deg is 5x the refraction uncertainty. Before
# Oct 2026 the prediction had half the dip and read +0.06-0.07 deg on a perfect camera.
OFFSET_REPORT_DEG = 0.05


def pointing_offset(info, eo_file, rate=True):
    """(tilt, roll, text) of the window's constant horizon offset against the calibration when
    either exceeds OFFSET_REPORT_DEG (from check_pointing's info; an older result without the
    numbers keeps its own text), else None. rate=False leaves out the rough tilt-only rate (for a
    README that gives the per-band table from the date's own DEM, which counts the roll too)."""
    if not info:
        return None
    dt, dr = info.get("window_horizon_dtilt_deg"), info.get("window_horizon_droll_deg")
    if (dt is None or dr is None) and info.get("calibration_offset") and info.get("calibration_offset_deg"):
        a, b = info["calibration_offset_deg"]
        return a, b, info["calibration_offset"]
    if dt is None or dr is None or not info.get("clear_horizon_days") or \
            max(abs(dt), abs(dr)) <= OFFSET_REPORT_DEG:
        return None
    px = info.get("window_horizon_offset_px")
    where = ""
    pr_, ic_, H_ = info.get("horizon_pred_rows"), info.get("horizon_cols_in_frame"), info.get("frame_rows")
    if pr_ and ic_ and H_:
        off_frame = ic_[0] < ic_[1]
        top = min(pr_) < 0.1 * H_
        where = (f". The calibration puts the horizon at rows {max(min(pr_), 0):.0f}-{max(pr_):.0f} of the {H_}-row "
                 f"frame" + (f", inside the frame on only {ic_[0]} of the {ic_[1]} columns checked" if off_frame else "")
                 + (": partly OFF the frame, where the lens model is extrapolated, so the offset there may be the "
                    "lens model's rather than the pointing's" if off_frame else
                    ": near the top edge, far from the GCPs, where the lens model is least constrained" if top else ""))
    text = (f"the sea horizon sits " + (f"{px:.0f} px " if px is not None else "")
            + f"from where {eo_file} puts it over the window ({info.get('clear_horizon_days')} day(s) with a clear "
            f"horizon; as a tilt/roll change: {dt:+.2f}/{dr:+.2f} deg). A constant offset like this is either "
            f"a calibration error or the lens model's error near the top of the frame (the horizon lies far "
            f"outside the GCPs the calibration was solved from); only day-to-day CHANGES are used to leave "
            f"days out. If it is pointing, it matters: the view grazes the beach, so 0.1 deg of tilt moves a "
            f"waterline ~6-9 m along a range of 250-350 m, ~0.13-0.17 m of DEM there"
            + (f" (rough rate, tilt only, roll ignored: ~{abs(dt) / 0.1 * 0.15:.2f} m for this offset at ~300 m; "
               f"the README gives it per range band from this date's own DEM, roll included)" if rate else "")
            + where)
    return dt, dr, text


POINTING_VERSION = 5


def pointing_sig(plan, cam, args, photos):
    c = plan["cams"][cam]
    return {"eo": c["eo_file"], "eo_file": file_sig(c["eo_path"]),
            "photos": hashlib.sha1("\n".join(p.name for p in photos).encode()).hexdigest(),
            "frames": args.pointing_frames, "fail": args.pointing_fail, "v": POINTING_VERSION,
            "scripts": scripts_sig("pointing")}


def step_pointing(plan, args, state):
    out = plan["out"]
    pdir = out / "pointing"
    results = {}
    any_built = False
    for cam, c in plan["cams"].items():
        photos, cnt = scan_photos(plan, cam, args)
        state.setdefault("photo_counts", {})[cam] = cnt
        sig = pointing_sig(plan, cam, args, photos)
        csv_out, js = pdir / f"pointing_{cam}.csv", pdir / f"pointing_{cam}.json"
        fresh, why = is_fresh(out, f"pointing_{cam}", sig, [csv_out, js], [c["eo_path"]])
        if fresh and not args.force:
            say(f"pointing {cam}", f"SKIPPED: {csv_out} is up to date (--force to rebuild)")
            results[cam] = json.loads(js.read_text())
            # the same photos, found in these roots now: the roots the result stands for
            update_stamp(out, f"pointing_{cam}", photo_roots=list(args.photo_roots))
            continue
        drop_stamp(out, f"pointing_{cam}")
        pdir.mkdir(parents=True, exist_ok=True)
        if not photos:
            warn(f"pointing {cam}: no photos found (roots: {', '.join(args.photo_roots)}); not checked")
            results[cam] = {"rows": [], "info": {"note": "no photos found: not checked"}, "different_days": [],
                            "leave_out": False,
                            "station_log": station_pointing_log(plan, cam, args) if plan["era"] == "live" else None}
            for p in (csv_out, js):
                if p.exists():
                    p.unlink()
            js.write_text(json.dumps(results[cam], indent=1) + "\n")   # read back by the README's caveats
            csv_out.write_text("camera,date,verdict\n")
            write_stamp(out, f"pointing_{cam}", sig, [f"check_pointing() for {cam}: no photos"],
                        scripts=STEP_SCRIPTS["pointing"], extra={"photo_roots": list(args.photo_roots)})
            continue
        if not c["eo_path"]:
            raise StepFailed(f"{cam}: calibration {c['eo_file']} not found in {args.calibration}")
        rows, info = check_pointing(plan, cam, photos, args)
        moved = [r["date"] for r in rows if r["verdict"] == "DIFFERENT POINTING"]
        measured = [r for r in rows if r["verdict"] != "unchecked"]
        # Only what was MEASURED is stored; whether those days are left out is decided when the
        # detection runs (--keep-moved-days), so changing that option takes effect on a rerun.
        leave_out = True
        if moved and len(moved) > len(measured) / 2.0:
            info["warning"] = (f"{len(moved)} of {len(measured)} measured days differ from the reference "
                               f"{info['reference_day']}: the reference itself may be the odd one out. Nothing "
                               f"is left out automatically: look at {pdir} and pass --skip-days if needed")
            leave_out = False
        fields = ["camera", "date", "photos", "verdict", "why", "horizon", "horizon_offset_px",
                  "horizon_dtilt_deg", "horizon_droll_deg", "horizon_fit_px", "horizon_frames", "features",
                  "feature_frame", "inliers", "d_azimuth_deg", "d_tilt_deg", "d_roll_deg"]
        with open(csv_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        res = {"rows": rows, "info": info, "different_days": moved, "leave_out": leave_out,
               "station_log": station_pointing_log(plan, cam, args) if plan["era"] == "live" else None}
        js.write_text(json.dumps(res, indent=1) + "\n")
        write_stamp(out, f"pointing_{cam}", sig, [f"survey_products.py, in-process: check_pointing() for {cam} "
                                                  f"(horizon_check.fit_tilt_roll + pointing_check.measure)"],
                    scripts=STEP_SCRIPTS["pointing"], extra={"photo_roots": list(args.photo_roots)})
        results[cam] = res
        any_built = True
    for cam, res in results.items():
        print_pointing(cam, res, args)
    state["pointing"] = results
    if results and all(not r.get("rows") for r in results.values()):
        return "not checked (no photos found)"
    return "built" if any_built else "skipped"


def skip_days_of(res, args):
    """The days the detection leaves out, from a pointing result and this run's options."""
    if not res or args.keep_moved_days or not res.get("leave_out", True):
        return []
    return sorted(res.get("different_days") or [])


def station_pointing_log(plan, cam, args):
    """Per day of the window: the statuses in the station's own pointing monitor log
    (pointing_check.py, archive/pointing_<cam>.csv), when there is one."""
    log = Path(args.pointing_log_dir) / f"pointing_{cam}.csv"
    if not log.exists():
        return None
    c = plan["cams"][cam]
    days = {}
    with open(log, newline="") as f:
        for r in csv.DictReader(f):
            t = utc_of(r.get("filename", "")) or (datetime.fromtimestamp(float(r["epoch"]), tz=timezone.utc)
                                                  if r.get("epoch") else None)
            if t is None:
                continue
            d = t.strftime("%Y-%m-%d")
            if not (c["first"] <= d <= c["last"]):
                continue
            e = days.setdefault(d, {"frames": 0, "status": {}, "d_azimuth": []})
            e["frames"] += 1
            e["status"][r.get("status", "?")] = e["status"].get(r.get("status", "?"), 0) + 1
            try:
                e["d_azimuth"].append(float(r["d_azimuth"]))
            except (KeyError, ValueError):
                pass
    for e in days.values():
        e["median_d_azimuth_deg"] = round(float(np.median(e["d_azimuth"])), 3) if e["d_azimuth"] else None
        del e["d_azimuth"]
    return {"log": str(log), "days": days}


def print_pointing(cam, res, args):
    info = res.get("info", {})
    sl = res.get("station_log")
    if sl:
        print(f"\n  {cam}: the station's own pointing monitor ({sl['log']}):")
        for d, e in sorted(sl["days"].items()):
            print(f"  {d}  {e['frames']:3d} frames  {e['status']}  median pan change {e['median_d_azimuth_deg']}")
        if any(k in ("MOVED", "moved") for e in sl["days"].values() for k in e["status"]):
            warn(f"{cam}: the station's pointing monitor flagged frames in the window (see above)")
    if not res.get("rows"):
        return
    print(f"\n  {cam}: reference {info.get('reference_frame')} -- {info.get('reference_kind')}")
    print("  date        photos  horizon(px)  tilt/roll chg (deg)  pan/tilt/roll vs ref (deg)  inliers  verdict")
    for r in res["rows"]:
        ang = (f"{r['d_azimuth_deg']:+.2f}/{r['d_tilt_deg']:+.2f}/{r['d_roll_deg']:+.2f}"
               if r["d_azimuth_deg"] != "" else r["features"][:24])
        hz = (f"{r['horizon_dtilt_deg']:+.2f}/{r['horizon_droll_deg']:+.2f}"
              if r["horizon_dtilt_deg"] != "" else r["horizon"])
        off = f"{r['horizon_offset_px']:.1f}" if r["horizon_offset_px"] != "" else "-"
        print(f"  {r['date']}  {r['photos']:5d}  {off:>10s}  {hz:>19s}  {ang:>26s}  {str(r['inliers']):>7s}  "
              f"{r['verdict']}")
    off = pointing_offset(info, info.get("calibration"))
    if off:
        warn(f"{cam}: {off[2]}")
    if info.get("warning"):
        warn(f"{cam}: {info['warning']}")
    if res.get("different_days"):
        skip = skip_days_of(res, args)
        if skip:
            warn(f"{cam}: LEFT OUT of the detection (different pointing): {', '.join(skip)} "
                 f"(--keep-moved-days keeps them)")
        else:
            warn(f"{cam}: days with a different pointing KEPT: {', '.join(res['different_days'])}"
                 + (" (--keep-moved-days)" if args.keep_moved_days else " (most days differ: see above)"))


def pointing_skip_days(plan, state, cam, args, photos):
    """
    -> (days the detection leaves out, note). The pointing results must have been made for
    these photos and this calibration; results of another window or calibration stop the
    detection (they would leave the wrong days out) rather than being used.
    """
    res = (state.get("pointing") or {}).get(cam)
    out = plan["out"]
    js = out / "pointing" / f"pointing_{cam}.json"
    if res is None:
        if not js.exists():
            return [], "pointing not checked (the pointing step has not run): no day left out"
        st = read_stamp(out, f"pointing_{cam}")
        if not st or st.get("sig") != normal(pointing_sig(plan, cam, args, photos)):
            raise StepFailed(f"{cam}: the pointing check in {js} was made for other photos or another "
                             f"calibration; run the pointing step again (--steps pointing,detect,...)")
        res = json.loads(js.read_text())
    return skip_days_of(res, args), ""


# ---------------------------------------------------------------------
# Step: detect (historical eras) / select + setup (live era), then merge
# ---------------------------------------------------------------------

def detect_cmd(plan, cam, args, skip_days):
    c = plan["cams"][cam]
    f = plan["out"] / "forcing"
    cmd = [sys.executable, HERE / "detect_original_view.py", "--camera", cam,
           "--eo", c["eo_path"], "--envelope-survey", c["envelope_path"],
           "--start-date", c["first"], "--end-date", c["last"], "--utc-hours", c["hours_text"],
           "--originals"] + list(args.photo_roots) + \
          ["--water-level", f / "water_level.csv", "--waves", f / "waves.csv",
           "--setup-coef", args.setup_coef, "--merge-into", "", "--output", plan["out"] / "waterlines" / cam]
    if c["station"]:
        cmd += ["--station", c["station"]]
    if skip_days:
        cmd += ["--skip-days"] + sorted(skip_days)
    return cmd


def step_detect(plan, args, state):
    out = plan["out"]
    wdir = out / "waterlines"
    f = out / "forcing"
    if not (f / "water_level.csv").exists():
        raise StepFailed(f"no forcing in {f}: run the forcing step first")
    fit = setup_fit_info(args.setup_coef, args)
    if plan["era"] == "live":
        parts = live_detect(plan, args, state, fit)
    else:
        parts = []
        for cam, c in plan["cams"].items():
            ground = wdir / cam / "contour_points_ground.csv"
            if not c["eo_path"]:
                raise StepFailed(f"{cam}: calibration {c['eo_file']} not found in {args.calibration}")
            if not c["envelope_path"]:
                raise StepFailed(f"{cam}: envelope survey {c['envelope_survey']} not found in "
                                 f"{', '.join(args.survey_dirs)} (--survey-dirs)")
            all_photos, cnt = scan_photos(plan, cam, args)
            skip, note = pointing_skip_days(plan, state, cam, args, all_photos)
            if note:
                warn(f"{cam}: {note}")
                state.setdefault("pointing_notes", {})[cam] = note
            photos = [p for p in all_photos if utc_of(p.name).strftime("%Y-%m-%d") not in set(skip)]
            cmd = detect_cmd(plan, cam, args, skip)
            sig = {"cmd": cmd_text(cmd), "photos": hashlib.sha1("\n".join(p.name for p in photos).encode()).hexdigest(),
                   "after": {"forcing": token(out, "forcing")}, "eo_file": file_sig(c["eo_path"]),
                   "envelope": file_sig(c["envelope_path"]), "setup_fit": fit_sig(fit), "scripts": scripts_sig("detect")}
            fresh, why = is_fresh(out, f"detect_{cam}", sig, [ground],
                                  [f / "water_level.csv", f / "waves.csv", c["eo_path"]])
            if fresh and not args.force:
                say(f"detect {cam}", f"SKIPPED: {ground} is up to date (--force to rebuild)")
                parts.append((cam, ground))
                continue
            drop_stamp(out, f"detect_{cam}")
            if not photos:
                warn(f"detect {cam}: no photos for {c['first']} .. {c['last']} (station {c['station']}) in "
                     f"{', '.join(args.photo_roots)}" + (f"; missing folders: {', '.join(cnt['missing_roots'])}"
                                                          if cnt["missing_roots"] else "") + ": camera left out")
                state.setdefault("failed_cams", {})[cam] = ("no photos" + (
                    f" ({cnt['wrong_size']} copies of another size than the lens file's)" if cnt.get("wrong_size")
                    else ""))
                continue
            say(f"detect {cam}", f"{len(photos)} photos ({why}); expected ~{len(photos) * SEC_PER_PHOTO_DETECT / 60:.0f} "
                f"min (~{SEC_PER_PHOTO_DETECT:.0f} s per photo as measured off the station; the progress lines "
                f"below print the real rate; follow with tail -f {wdir / cam / 'processing.log'})")
            note = overlap_note(len(photos) * SEC_PER_PHOTO_DETECT)
            if note:
                warn(f"detect {cam}: {note}")
            check_disk(out, len(photos) * DETECT_MB_PER_PHOTO + PRODUCT_MB,
                       f"the detection of {len(photos)} {cam} photos (~{DETECT_MB_PER_PHOTO:.0f} MB a photo of debug "
                       f"images and overlays until the camera is done, + ~{PRODUCT_MB / 1e3:.1f} GB of products)")
            rc = run_cmd(cmd, out / "logs" / f"detect_{cam}.log")
            if rc != 0 or not ground.exists():
                warn(f"detect {cam} FAILED (exit {rc}); see {out / 'logs' / f'detect_{cam}.log'} and "
                     f"{wdir / cam / 'processing.log'}: camera left out")
                state.setdefault("failed_cams", {})[cam] = f"detection failed (exit {rc})"
                continue
            write_stamp(out, f"detect_{cam}", sig, [cmd], scripts=STEP_SCRIPTS["detect"],
                        extra={"setup_coef": args.setup_coef, "setup_fit": fit, "skip_days": skip,
                               "photos": len(photos)})
            if not args.keep_detector_debug:
                prune_detector_scratch(wdir / cam)
            parts.append((cam, ground))
    if not parts:
        raise StepFailed("no camera produced waterlines")
    merged = wdir / "contour_points_ground.csv"
    sig = {"parts": [[c, str(p)] for c, p in parts], "era": plan["era"], "v": 2,
           "after": dict({f"detect_{c}": token(out, f"detect_{c}") for c, _ in parts},
                         forcing=token(out, "forcing")),
           "drop_no_setup": bool(args.setup_coef)}
    fresh, why = is_fresh(out, "merge", sig, [merged, wdir / "photos"],
                          [p for _, p in parts] + [f / "water_level.csv", f / "waves.csv"])
    if fresh and not args.force:
        say("merge", f"SKIPPED: {merged.name} is up to date")
        return "skipped"
    drop_stamp(out, "merge")
    # the filtered file and everything built on the old merged file are now stale
    drop_stamp(out, "filter")
    for p in (wdir / "contour_points_ground_filtered.csv",):
        if p.exists():
            p.unlink()
    stats = merge_cameras(plan, parts, args)
    link_all_photos(plan)
    setup = {}
    for c, _ in parts:
        st = read_stamp(out, f"detect_{c}") or {}
        setup[c] = {"coef": st.get("setup_coef"), "fit": st.get("setup_fit"),
                    "implied_coef": (stats["implied_coef"] or {}).get(c)}
        imp = setup[c]["implied_coef"]
        if imp is not None and st.get("setup_coef") is not None and abs(imp - st["setup_coef"]) > 2e-4:
            warn(f"{c}: the rows' own setup implies C = {imp:.4f}, but the detection stamp says "
                 f"{st['setup_coef']}: the rows are taken as the truth")
            setup[c]["coef"] = imp
            setup[c]["fit"] = setup_fit_info(imp, None)
    write_stamp(out, "merge", sig, extra={"setup": setup, "rows": stats["rows"],
                                          "no_setup_frames": stats["no_setup_frames"],
                                          "no_setup_rows": stats["no_setup_rows"],
                                          "cameras": [c for c, _ in parts],
                                          "failed_cams": state.get("failed_cams") or {}})
    return "built"


def check_disk(path, need_mb, what):
    """StepFailed unless the disk holding `path` has need_mb free (MB)."""
    try:
        free = shutil.disk_usage(str(path)).free / 1e6
    except OSError:
        return
    if free < need_mb:
        raise StepFailed(f"only {free / 1e3:.2f} GB free on the disk of {path}; {what} needs ~{need_mb / 1e3:.2f} GB "
                         f"while it runs: free some space or give --output-root on another disk (nothing was rebuilt "
                         f"by this step)")


def prune_detector_scratch(cdir):
    """Deletes the detector's debug images (debug/, ~8 MB a frame), its resized inputs (in/) and
    its full-size PNG overlays (detections/*_overlay.png, ~4 MB a frame) once the waterlines are
    written: ~2 GB for a two-camera week otherwise. The detections' CSVs, the half-size overlays
    of detect_original_view.py (overlays/: the final line and the search envelope on each photo)
    and the photo links stay. --keep-detector-debug keeps everything."""
    freed = 0
    for sub in ("debug", "in"):
        d = Path(cdir) / sub
        if d.is_dir():
            freed += sum(f.stat().st_size for f in d.rglob("*") if f.is_file() and not f.is_symlink())
            shutil.rmtree(d, ignore_errors=True)
    for f in (Path(cdir) / "detections").glob("*_overlay.png"):
        freed += f.stat().st_size
        f.unlink()
    if freed:
        say("disk", f"removed the detector's debug images and resized inputs in {cdir} ({freed / 1e6:.0f} MB; "
            f"--keep-detector-debug keeps them)")


LIVE_DETECT_VERSION = 3


def live_detect(plan, args, state, fit):
    """The cron's detections for the window: setup applied, georectified with the live calibration."""
    out = plan["out"]
    wdir = out / "waterlines"
    rows_file = live_rows_file(plan, args, state)
    parts = []
    for cam, c in plan["cams"].items():
        if not c["eo_path"]:
            raise StepFailed(f"{cam}: calibration {c['eo_file']} not found in {args.calibration}")
        all_photos, _ = scan_photos(plan, cam, args)
        skip, note = pointing_skip_days(plan, state, cam, args, all_photos)
        if note:
            warn(f"{cam}: {note}")
            state.setdefault("pointing_notes", {})[cam] = note
        cdir = wdir / cam
        contours, ground = cdir / "contour_points.csv", cdir / "contour_points_ground.csv"
        geo = [sys.executable, HERE / "georectify.py", contours, ground,
               f"--io-{cam}", c["io_path"], f"--eo-{cam}", c["eo_path"]]
        sig = {"rows": str(rows_file), "after": {"archive_rows": token(out, "archive_rows")},
               "setup_coef": args.setup_coef, "setup_fit": fit_sig(fit), "skip": sorted(skip), "cmd": cmd_text(geo),
               "eo_file": file_sig(c["eo_path"]), "v": LIVE_DETECT_VERSION, "scripts": scripts_sig("detect_live")}
        fresh, why = is_fresh(out, f"detect_{cam}", sig, [ground], [rows_file, c["eo_path"]])
        if fresh and not args.force:
            say(f"detect {cam}", f"SKIPPED: {ground} is up to date (--force to rebuild)")
            parts.append((cam, ground))
            continue
        drop_stamp(out, f"detect_{cam}")
        cdir.mkdir(parents=True, exist_ok=True)
        stats = apply_setup(rows_file, contours, cam, args.setup_coef, skip)
        say(f"setup {cam}", f"{stats['rows']} rows, {stats['frames']} frames; C = {args.setup_coef}: "
            + (f"{stats['setup_min']:.3f}-{stats['setup_max']:.3f} m (median {stats['setup_median']:.3f})"
               if stats["frames_with_setup"] else "none")
            + (f"; {stats['frames_without_waves']} frame(s) without Hs+Tp left uncorrected"
               if stats["frames_without_waves"] else "")
            + (f"; {stats['skipped_rows']} rows on left-out days" if stats["skipped_rows"] else ""))
        if stats.get("had_setup"):
            say(f"setup {cam}", f"the archive rows already had a setup column: it differs from this run's by "
                f"up to {stats['max_diff_existing']:.4f} m (this run's C = {args.setup_coef} is used)")
        if not stats["rows"]:
            warn(f"{cam}: no archive rows in the window: camera left out")
            state.setdefault("failed_cams", {})[cam] = "no archive rows"
            continue
        rc = run_cmd(geo, out / "logs" / f"detect_{cam}.log")
        if rc != 0 or not ground.exists():
            warn(f"georectify {cam} FAILED (exit {rc}): camera left out")
            state.setdefault("failed_cams", {})[cam] = f"georectification failed (exit {rc})"
            continue
        # photos of the frames, for the maps (the maps step links them again from the photo roots of
        # its own run, so a photo that turns up later is used without a new georectification)
        link_frame_photos(plan, cam, args)
        write_stamp(out, f"detect_{cam}", sig, [geo], scripts=STEP_SCRIPTS["detect_live"],
                    extra={"setup": {k: v for k, v in stats.items() if k != "frame_names"},
                           "setup_coef": args.setup_coef, "setup_fit": fit, "skip_days": skip})
        parts.append((cam, ground))
    return parts


def frame_names(path):
    """The source_file names in a contour file (as photo file names, .jpg)."""
    names = set()
    try:
        with open(path, newline="") as f:
            for r in csv.DictReader(f):
                n = r.get("source_file") or ""
                if n:
                    names.add(n if n.endswith(".jpg") else n + ".jpg")
    except OSError:
        pass
    return names


def link_frame_photos(plan, cam, args):
    """
    Live era: waterlines/<cam>/src holds the photos of the camera's frames found in the photo roots
    NOW (scan_photos of this run), relinked when that set changed, e.g. after a corrected
    --photo-roots or once a mount is back; photos.txt lists them. -> (photos linked, frames).
    """
    cdir = plan["out"] / "waterlines" / cam
    frames = frame_names(cdir / "contour_points.csv")
    all_photos, _ = scan_photos(plan, cam, args)
    found = sorted((p for p in all_photos if p.name in frames), key=lambda p: p.name)
    src = cdir / "src"
    have = sorted(q.name for q in src.glob("*.jpg")) if src.exists() else []
    if have != [p.name for p in found]:
        shutil.rmtree(src, ignore_errors=True)
        link_photos(found, src)
        say(f"photos {cam}", f"{len(found)} of {len(frames)} frames' photos found in the photo roots -> {src}")
    (cdir / "photos.txt").write_text("".join(f"{p}\n" for p in found))
    return found, len(frames)


def apply_setup(src, dst, cam, coef, skip_days=()):
    """
    Writes the camera's rows of `src` to `dst` with setup_correction_m and
    beach_elevation_navd88 computed exactly as extract_elevation_contours.py
    --setup-coef does: per frame C*sqrt(Hs*g*Tp^2/(2*pi)) from the row's
    offshore_hs_m and offshore_tp_s (already rounded to 0.01 m and 0.1 s
    there), rounded to 0.1 mm; the beach elevation is the tide elevation plus
    the setup, rounded to 0.1 mm; a frame without both keeps a blank setup
    and its tide elevation. Existing setup columns (rows the cron wrote after
    6 Oct 2026) are recomputed with this C and compared.
    """
    skip = set(skip_days or ())
    st = {"rows": 0, "frames": 0, "frames_with_setup": 0, "frames_without_waves": 0, "skipped_rows": 0,
          "had_setup": False, "max_diff_existing": 0.0, "frame_names": set()}
    setups = {}
    with open(src, newline="") as fi, open(str(dst) + ".tmp", "w", newline="") as fo:
        rd = csv.DictReader(fi)
        fields = [k for k in rd.fieldnames if k not in ("setup_correction_m", "beach_elevation_navd88")]
        st["had_setup"] = "setup_correction_m" in rd.fieldnames
        if "offshore_hs_m" not in rd.fieldnames or "offshore_tp_s" not in rd.fieldnames:
            raise StepFailed(f"{src}: no offshore_hs_m / offshore_tp_s columns (the cron ran without --waves?)")
        w = csv.DictWriter(fo, fieldnames=fields + ["setup_correction_m", "beach_elevation_navd88"],
                           extrasaction="ignore")
        w.writeheader()
        for r in rd:
            if r["camera"] != cam:
                continue
            if r["capture_time_utc"][:10] in skip:
                st["skipped_rows"] += 1
                continue
            hs, tp = r["offshore_hs_m"], r["offshore_tp_s"]
            k = r["source_file"]
            if hs != "" and tp != "":
                full = float(setup_m(hs, tp, coef))          # as there: the beach uses the unrounded setup
                s = round(full, 4)
                beach = round(round(float(r["tide_elevation_navd88"]), 4) + full, 4)
                setups[k] = s
            else:
                s, beach = "", r["tide_elevation_navd88"]
                setups.setdefault(k, None)
            if st["had_setup"] and r.get("setup_correction_m", "") != "" and s != "":
                st["max_diff_existing"] = max(st["max_diff_existing"], abs(float(r["setup_correction_m"]) - s))
            r["setup_correction_m"], r["beach_elevation_navd88"] = s, beach
            w.writerow(r)
            st["rows"] += 1
            st["frame_names"].add(k)
    os.replace(str(dst) + ".tmp", dst)
    vals = [v for v in setups.values() if v is not None]
    st["frames"] = len(setups)
    st["frames_with_setup"] = len(vals)
    st["frames_without_waves"] = len(setups) - len(vals)
    if vals:
        st.update(setup_min=min(vals), setup_max=max(vals), setup_median=float(np.median(vals)))
    return st


def forcing_lookup(fdir):
    """(level epochs, sources, sigmas, wave epochs, wave sources, hs sigmas) from the forcing files."""
    import pandas as pd
    wl = pd.read_csv(fdir / "water_level.csv")
    wv = pd.read_csv(fdir / "waves.csv")
    t = pd.to_datetime(wl["time"], utc=True)
    wl_ep = ((t - pd.Timestamp("1970-01-01", tz="UTC")) // pd.Timedelta(seconds=1)).to_numpy(float)
    o = np.argsort(wl_ep)
    wv_ep = wv["epoch"].to_numpy(float)
    p = np.argsort(wv_ep)
    col = lambda df, k: df[k].astype(str).to_numpy() if k in df else np.array([""] * len(df))  # noqa: E731
    return (wl_ep[o], col(wl, "source")[o], col(wl, "sigma_m")[o],
            wv_ep[p], col(wv, "source")[p], col(wv, "hs_sigma_m")[p])


def nearest(ep, t, max_gap):
    if not len(ep):
        return None
    i = int(np.searchsorted(ep, t))
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(ep) and abs(ep[j] - t) <= max_gap and (best is None or abs(ep[j] - t) < abs(ep[best] - t)):
            best = j
    return best


def implied_c(setup, hs, tp):
    """C = setup / sqrt(Hs*L0) of one row, or None."""
    try:
        s, h, p = float(setup), float(hs), float(tp)
    except (TypeError, ValueError):
        return None
    x = np.sqrt(h * G * p ** 2 / (2 * np.pi))
    return s / x if x > 0.05 else None


def merge_cameras(plan, parts, args):
    """
    waterlines/contour_points_ground.csv: every camera's rows, each tagged with its water-level and
    wave source (forcing_level_source, forcing_level_sigma_m, forcing_waves_source).
    With a setup coefficient (C > 0) a frame whose setup could not be computed (no Hs/Tp within
    the wave record's gap limit) is LEFT OUT: its elevation would be the still water only, about
    one setup (0.2-0.4 m) below its neighbours, and dem_from_contours.py would keep it. Those
    frames are listed in waterlines/no_setup_frames.csv and counted per camera.
    -> stats: rows, no_setup_frames {cam: n}, no_setup_rows, implied_coef {cam: median C of the rows}.
    """
    out = plan["out"]
    dst = out / "waterlines" / "contour_points_ground.csv"
    heads = []
    for cam, p in parts:
        with open(p, newline="") as f:
            heads.append(next(csv.reader(f)))
    fields = []
    for h in heads:
        fields += [k for k in h if k not in fields]
    tag = plan["era"] != "live"
    if tag:
        lk = forcing_lookup(out / "forcing")
        fields += ["forcing_level_source", "forcing_level_sigma_m", "forcing_waves_source"]
    drop = bool(args.setup_coef) and "setup_correction_m" in fields
    if args.setup_coef and "setup_correction_m" not in fields:
        warn("the waterlines have no setup_correction_m column although C > 0: no setup was applied")
    n = 0
    cache = {}
    no_setup = {}                 # (cam, frame) -> (time, rows)
    cs = {}                       # cam -> {frame: implied C}
    with open(str(dst) + ".tmp", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=fields, extrasaction="ignore", restval="")
        w.writeheader()
        for cam, p in parts:
            with open(p, newline="") as fi:
                for r in csv.DictReader(fi):
                    k = r["source_file"]
                    if drop and (r.get("setup_correction_m") or "") == "":
                        e = no_setup.setdefault((r.get("camera") or cam, k), [r.get("capture_time_utc", ""), 0])
                        e[1] += 1
                        continue
                    fc = cs.setdefault(r.get("camera") or cam, {})
                    if k not in fc and r.get("setup_correction_m"):
                        fc[k] = implied_c(r["setup_correction_m"], r.get("offshore_hs_m"), r.get("offshore_tp_s"))
                    if tag:
                        if k not in cache:
                            t = float(r["capture_epoch"])
                            i = nearest(lk[0], t, 3600)
                            j = nearest(lk[3], t, 5400)
                            cache[k] = ((lk[1][i] if i is not None else "", lk[2][i] if i is not None else "",
                                         lk[4][j] if j is not None else ""))
                        r["forcing_level_source"], r["forcing_level_sigma_m"], r["forcing_waves_source"] = cache[k]
                    w.writerow(r)
                    n += 1
    os.replace(str(dst) + ".tmp", dst)
    ns = out / "waterlines" / "no_setup_frames.csv"
    with open(ns, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["camera", "source_file", "capture_time_utc", "rows", "why"])
        for (cam, k), (t, nr) in sorted(no_setup.items(), key=lambda kv: kv[1][0]):
            wr.writerow([cam, k, t, nr, "no wave height and period within the wave record's gap limit: "
                                        "setup unknown, left out"])
    per_cam = {}
    rows_out = 0
    for (cam, _), (_, nr) in no_setup.items():
        per_cam[cam] = per_cam.get(cam, 0) + 1
        rows_out += nr
    implied = {}
    for cam, fc in cs.items():
        v = [x for x in fc.values() if x is not None]
        implied[cam] = round(float(np.median(v)), 5) if v else None
    say("merged", f"{n} waterline points from {', '.join(c for c, _ in parts)} -> {dst}")
    if per_cam:
        warn("frames WITHOUT a wave record (no setup) left out of the waterlines: "
             + ", ".join(f"{c} {v}" for c, v in sorted(per_cam.items())) + f" (listed in {ns})")
    return {"rows": n, "no_setup_frames": per_cam, "no_setup_rows": rows_out, "implied_coef": implied}


def link_all_photos(plan):
    """waterlines/photos: every camera's photos in one folder (the consistency plot's backdrop)."""
    wdir = plan["out"] / "waterlines"
    dst = wdir / "photos"
    shutil.rmtree(dst, ignore_errors=True)
    dst.mkdir(parents=True)
    n = 0
    for cam in plan["cams"]:
        for p in sorted((wdir / cam / "src").glob("*.jpg")) if (wdir / cam / "src").exists() else []:
            try:
                os.link(p, dst / p.name)
            except OSError:
                try:
                    os.symlink(p.resolve(), dst / p.name)
                except OSError:
                    shutil.copy2(p, dst / p.name)
            n += 1
    return n


# ---------------------------------------------------------------------
# Step: filter
# ---------------------------------------------------------------------

def step_filter(plan, args, state):
    out = plan["out"]
    wdir = out / "waterlines"
    src, dst = wdir / "contour_points_ground.csv", wdir / "contour_points_ground_filtered.csv"
    if not src.exists():
        raise StepFailed(f"no {src}: run the detect step first")
    ndays = len(days_between(plan["first"], plan["last"]))
    cmd = [sys.executable, HERE / "waterline_consistency.py", src, "--output", dst,
           "--report", wdir / "consistency_report.csv", "--plot", wdir / "consistency",
           "--image-dir", wdir / "photos", "--plot-days", ndays]
    sig = {"cmd": cmd_text(cmd), "after": {"merge": token(out, "merge")}, "scripts": scripts_sig("filter")}
    fresh, why = is_fresh(out, "filter", sig, [dst], [src])
    if fresh and not args.force:
        say("filter", f"SKIPPED: {dst} is up to date (--force to rebuild)")
        return "skipped"
    say("filter", f"waterline_consistency.py, its defaults as the live cron ({why}); ~10-60 s on the NUC")
    drop_stamp(out, "filter")
    if dst.exists():
        dst.unlink()                       # never leave an old filtered file looking current
    log = out / "logs" / "filter.log"
    if log.exists():
        log.unlink()                       # its settings are read back from this run's log
    rc = run_cmd(cmd, log)
    summary, settings = "", {}
    try:
        lines = log.read_text().splitlines()
        summary = [ln for ln in lines if ln.startswith("CONSISTENCY ")][-1][12:]
        for ln in lines:                   # 'Window : ...', 'Thresholds : ...' as the script prints them
            m = re.match(r"^(Window|Thresholds|Reference|Bins?|Min\w*)\s*:\s*(.+)$", ln.strip())
            if m:
                settings[m.group(1).lower()] = m.group(2).strip()
    except (OSError, IndexError):
        pass
    if rc != 0 or not dst.exists() or dst.stat().st_size == 0:
        if dst.exists():
            dst.unlink()
        warn(f"waterline consistency filter FAILED (exit {rc}): the DEM, maps and comparison use the "
             f"UNFILTERED {src.name} (as the live cron does); see {log}")
        state["filter_failed"] = True
        return "failed"
    write_stamp(out, "filter", sig, [cmd], scripts=STEP_SCRIPTS["filter"],
                extra={"summary": summary, "settings": settings})
    return "built"


def waterline_file(plan, quiet_=False):
    """The waterlines the DEM, maps and comparisons use: the filtered file, but only when the filter
    ran on the CURRENT merged file (its stamp names the merge it followed); else the unfiltered
    merged file, with a warning, never a filtered file left from earlier settings."""
    out = plan["out"]
    wdir = out / "waterlines"
    f = wdir / "contour_points_ground_filtered.csv"
    m = wdir / "contour_points_ground.csv"
    if not f.exists():
        return m
    st = read_stamp(out, "filter")
    ok = (st is not None and (st.get("sig") or {}).get("after", {}).get("merge") == token(out, "merge")
          and token(out, "merge") is not None and (not m.exists() or f.stat().st_mtime + 1 >= m.stat().st_mtime))
    if ok:
        return f
    if not quiet_:
        warn(f"{f.name} was not made from the current {m.name} (older settings): it is NOT used; "
             f"run the filter step (--steps filter,dem,maps,compare)")
    return m


def waterline_token(plan):
    """The stamp token of the waterline file in use (filter or merge)."""
    src = waterline_file(plan, quiet_=True)
    return token(plan["out"], "filter" if src.name.endswith("_filtered.csv") else "merge")


# ---------------------------------------------------------------------
# Step: dem
# ---------------------------------------------------------------------

def view_calibration(plan):
    """dem/view_calibration: the date's own IO/EO under the names dem_figure.py reads, so the page
    draws the cameras where they stood then (dem_figure.py --calibration)."""
    d = plan["out"] / "dem" / "view_calibration"
    d.mkdir(parents=True, exist_ok=True)
    lines = []
    for cam, c in plan["cams"].items():
        if not c["eo_path"]:
            continue
        shutil.copy2(c["io_path"], d / f"CACO05_{cam}_20240801_IO.yaml")
        shutil.copy2(c["eo_path"], d / f"CACO05_{cam}_20251113_EO-CV.yaml")
        lines.append(f"{cam}: CACO05_{cam}_20251113_EO-CV.yaml here = {c['eo_file']}; "
                     f"CACO05_{cam}_20240801_IO.yaml = {Path(c['io_path']).name}")
    (d / "README.txt").write_text("Copies of this date's calibration under the file names dem_figure.py "
                                  "looks for, so the DEM page shows the cameras where they stood:\n"
                                  + "\n".join(lines) + "\n")
    return d


# dem_from_contours.py's day test compares each camera-day with the DEM built WITHOUT it. With fewer
# than 3 camera-days that DEM is the other day alone: two days that differ by more than the limit get
# equal and opposite offsets and BOTH are rejected, leaving no points (the script then stops with a
# traceback). With fewer than this many camera-days the test is not run, and the README says so.
DEM_DAY_TEST_MIN = 3


def camera_days(plan, args, src=None):
    """The (camera, local day) pairs dem_from_contours.py's day test would see: frames of the window
    (UTC day of capture_time_utc) with ground coordinates, not left out by its wave filter (Hs >
    --dem-max-hs; a frame with no wave record is kept, as it does), keyed by the station's local day
    (America/New_York) as it keys them. -> sorted list of 'cam YYYY-MM-DD'."""
    src = Path(src or waterline_file(plan, quiet_=True))
    max_hs = float(args.dem_max_hs) if args.dem_max_hs else None
    epochs = {}
    try:
        with open(src, newline="") as f:
            for r in csv.DictReader(f):
                if not r.get("easting_utm19") or not r.get("northing_utm19"):
                    continue
                day = (r.get("capture_time_utc") or "")[:10]
                if day < plan["first"] or day > plan["last"]:
                    continue
                k = (r.get("camera", ""), r.get("source_file", ""))
                if k in epochs:
                    continue
                hs = r.get("offshore_hs_m", "")
                if max_hs is not None and hs not in ("", None) and float(hs) > max_hs:
                    continue
                try:
                    epochs[k] = float(r.get("capture_epoch") or "nan")
                except ValueError:
                    epochs[k] = float("nan")
    except (OSError, ValueError):
        return []
    import pandas as pd
    keys = set()
    for (cam, _), ep in epochs.items():
        if ep == ep:
            keys.add(f"{cam} " + pd.Timestamp(int(ep), unit="s", tz="UTC").tz_convert("America/New_York")
                     .strftime("%Y-%m-%d"))
    return sorted(keys)


def day_test_off(plan, args, src=None):
    """Why the DEM's day test is not run for this window ('' when it is, or when it is off by option)."""
    if not args.dem_max_day_offset:
        return ""
    cd = camera_days(plan, args, src)
    if len(cd) >= DEM_DAY_TEST_MIN:
        return ""
    return (f"only {len(cd)} camera-day(s) in the window ({', '.join(cd) or 'none'}): with fewer than "
            f"{DEM_DAY_TEST_MIN}, the DEM without one day is the other day alone, so two days that differ by more "
            f"than {args.dem_max_day_offset} m would get equal and opposite offsets and both be rejected, leaving no "
            f"points. Every camera-day is kept; a difference between them stays in the DEM's spread")


def dem_cmds(plan, args, day_test=True):
    out = plan["out"]
    stem = out / "dem" / plan["date"]
    cmd = [sys.executable, HERE / "dem_from_contours.py", waterline_file(plan), stem,
           "--start-date", plan["first"], "--end-date", plan["last"],
           "--cell", args.dem_cell, "--min-points", args.dem_min_points, "--max-spread", args.dem_max_spread,
           "--no-plot"]
    if args.dem_max_hs:
        cmd += ["--max-hs", args.dem_max_hs]
    if args.dem_max_day_offset and day_test:
        cmd += ["--max-day-offset", args.dem_max_day_offset]
    tif = [sys.executable, HERE / "asc_to_geotiff.py"] + [f"{stem}_{k}.asc" for k in ("dem", "spread", "count")] + \
          ["--epsg", args.geotiff_epsg]
    station = "/".join(sorted({c["station"] for c in plan["cams"].values()}))
    cal = ", ".join(sorted({eo_date(c["eo_file"]) or "?" for c in plan["cams"].values()}))
    page = [sys.executable, HERE / "dem_figure.py", stem, "--calibration", out / "dem" / "view_calibration",
            "--station", f"{station} Marconi Beach (calibration {cal})"]
    return stem, cmd, tif, page


def step_dem(plan, args, state):
    out = plan["out"]
    src = waterline_file(plan)
    if not src.exists():
        raise StepFailed(f"no waterlines in {src.parent}: run the detect step first")
    if src.name == "contour_points_ground.csv":
        warn("the DEM is built from the UNFILTERED waterlines (no filtered file)")
    no_day_test = day_test_off(plan, args, src)
    if no_day_test:
        warn(f"DEM day test NOT run: {no_day_test}")
    stem, cmd, tif, page = dem_cmds(plan, args, day_test=not no_day_test)
    outs = [Path(f"{stem}_{k}.{e}") for k in ("dem", "spread", "count") for e in ("asc", "tif")] + \
           [Path(f"{stem}_dem.png")]
    src_key = "filter" if src.name.endswith("_filtered.csv") else "merge"
    sig = {"cmd": cmd_text(cmd), "tif": cmd_text(tif), "page": cmd_text(page),
           "after": {src_key: token(out, src_key)}, "scripts": scripts_sig("dem")}
    fresh, why = is_fresh(out, "dem", sig, outs, [src])
    if fresh and not args.force:
        say("dem", f"SKIPPED: {stem}_dem.asc/.tif/.png are up to date (--force to rebuild)")
        return "skipped"
    say("dem", f"dem_from_contours.py over {plan['first']} .. {plan['last']} ({why}); ~10-40 s on the NUC")
    drop_stamp(out, "dem")
    (out / "dem").mkdir(parents=True, exist_ok=True)
    for o in outs + [Path(f"{stem}_info.json")]:        # no stale grid or page from an earlier build
        if o.exists():
            o.unlink()
    clear_comparisons(plan, "the DEM is being rebuilt")
    log = out / "logs" / "dem.log"
    if log.exists():
        log.unlink()                        # its cell counts are read back from this run's log
    rc = run_cmd(cmd, log)
    if rc != 0 or not Path(f"{stem}_dem.asc").exists():
        raise StepFailed(f"dem_from_contours.py failed (exit {rc}): {dem_failure_reason(log)}; see {log}")
    cells = dem_log_counts(log)
    rc = run_cmd(tif, log)
    if rc != 0:
        raise StepFailed(f"asc_to_geotiff.py failed (exit {rc}); see {log}")
    view_calibration(plan)
    rc = run_cmd(page, log)
    extra = {"cells": cells, "source": str(src),
             "settings": dict(dem_settings_requested(args), day_test_off=no_day_test or None)}
    if rc != 0 or not Path(f"{stem}_dem.png").exists():
        warn(f"the DEM page (dem_figure.py) failed (exit {rc}); the grids are fine. See {log}")
        write_stamp(out, "dem", dict(sig, page_failed=True), [cmd, tif, page], scripts=STEP_SCRIPTS["dem"],
                    extra=extra)
        return "built (page failed)"
    write_stamp(out, "dem", sig, [cmd, tif, page], scripts=STEP_SCRIPTS["dem"], extra=extra)
    return "built"


DEM_SETTING_OPTS = {"cell_m": "--cell", "min_points": "--min-points", "max_spread_m": "--max-spread",
                    "max_hs_m": "--max-hs", "max_day_offset_m": "--max-day-offset"}


def dem_settings_requested(args):
    """The DEM and GeoTIFF settings this run asks for (0 = off -> None, as dem_cmds passes them)."""
    return {"cell_m": float(args.dem_cell), "min_points": int(args.dem_min_points),
            "max_spread_m": float(args.dem_max_spread), "max_hs_m": float(args.dem_max_hs) or None,
            "max_day_offset_m": float(args.dem_max_day_offset) or None, "geotiff_epsg": int(args.geotiff_epsg)}


def dem_settings_built(out):
    """The settings the DEM ON DISK was built with: from its stamp (recorded since Oct 2026, else read
    back from the dem_from_contours.py and asc_to_geotiff.py commands the stamp holds). None when no
    DEM has finished. A later run with other options (e.g. --steps compare --dem-cell 1) does not
    change them: the provenance, README and labels report these."""
    st = read_stamp(out, "dem")
    if not st:
        return None
    if st.get("settings"):
        return dict(st["settings"])
    sig = st.get("sig") or {}
    try:
        cmd, tif = shlex.split(sig.get("cmd", "")), shlex.split(sig.get("tif", ""))
    except ValueError:
        return None

    def opt(words, name, cast):
        return cast(words[words.index(name) + 1]) if name in words[:-1] else None
    d = {k: opt(cmd, o, int if k == "min_points" else float) for k, o in DEM_SETTING_OPTS.items()}
    d["geotiff_epsg"] = opt(tif, "--epsg", int)
    return d


def dem_settings_differ(built, wanted):
    """The settings that differ between the DEM on disk and this run's options, as text ('' if none)."""
    if not built:
        return ""
    diff = []
    for k, v in wanted.items():
        b = built.get(k)
        if (b is None) != (v is None) or (b is not None and abs(float(b) - float(v)) > 1e-9):
            diff.append(f"{k}: asked {v}, DEM on disk {b}")
    return ", ".join(diff)


def dem_failure_reason(log):
    """Why dem_from_contours.py stopped, from its log: every camera-day rejected by the day test (no
    points left), no points in the window, or its last error line. -> text"""
    try:
        lines = [ln.rstrip() for ln in Path(log).read_text().splitlines()]
    except OSError:
        return "no log"
    tested, rejected, in_day = 0, 0, False
    for ln in lines:
        if ln.startswith("Day consistency"):
            in_day, tested, rejected = True, 0, 0
            continue
        if in_day:
            m = re.match(r"^\s+(\S+ \d{4}-\d\d-\d\d)\s", ln)
            if m:
                tested += 1
                rejected += "REJECTED" in ln
            elif not ln.strip().startswith("left out"):
                in_day = False
    if tested and rejected == tested:
        return (f"the day-consistency test (--dem-max-day-offset) rejected every camera-day ({rejected} of "
                f"{tested}): no points were left to grid. The days disagree by more than the limit; with so few "
                f"days the test cannot tell which one is off: rerun with --dem-max-day-offset 0 (test off) or a "
                f"wider --window")
    for ln in lines:
        if ln.startswith("No georectified points matched"):
            return "no waterline points in the window (dem_from_contours.py: 'No georectified points matched')"
    err = [ln for ln in lines if re.match(r"^(ERROR|[A-Za-z]+Error)\b", ln.strip())]
    return err[-1].strip() if err else "see its log"


def dem_log_counts(log):
    """The cell counts dem_from_contours.py prints (cells with points / filled / blanked, and why).
    -> dict (empty if the lines are not there)."""
    out = {}
    try:
        for ln in Path(log).read_text().splitlines():
            m = re.match(r"^cells (with points|filled|blanked)\s*:\s*(\d+)\s*(.*)$", ln.strip())
            if m:
                k = {"with points": "with_points", "filled": "filled", "blanked": "blanked"}[m.group(1)]
                out[k] = int(m.group(2))
                if k == "blanked" and m.group(3):
                    out["blanked_why"] = m.group(3).strip().strip("()")
            m = re.match(r"^Wave filter\s*:\s*(.+)$", ln.strip())
            if m:
                out["wave_filter"] = m.group(1).strip()
    except OSError:
        pass
    return out


def comparison_files(plan, s):
    """Every file survey_compare.py writes for survey row s (and its merged points)."""
    name = f"{plan['date']}_{s['name']}"
    cdir = plan["out"] / "compare"
    return [cdir / f"{name}_{k}" for k in ("comparison.json", "comparison.png", "comparison.csv",
                                           "comparison.txt", "points.csv", "dem_minus_survey.asc",
                                           "dem_minus_survey.tif", "waterline_frames.csv")]


def clear_comparisons(plan, why):
    """Removes the date's comparison outputs (and their stamps): they compared another DEM."""
    n = 0
    for s in plan["surveys"]:
        for p in comparison_files(plan, s):
            if p.exists():
                p.unlink()
                n += 1
        drop_stamp(plan["out"], f"compare_{s['name']}")
    if n:
        say("compare", f"{n} comparison file(s) of the previous build removed ({why})")


# ---------------------------------------------------------------------
# Step: maps
# ---------------------------------------------------------------------

MAP_STATUS = "photo_maps.json"


def map_status(plan):
    """{camera: {"drawn": bool, "why": text, "photos": n}} recorded by the maps step (maps/photo_maps.json)."""
    return load_json(plan["out"] / "maps" / MAP_STATUS) or {}


def step_maps(plan, args, state):
    out = plan["out"]
    mdir = out / "maps"
    src = waterline_file(plan)
    if not src.exists():
        raise StepFailed(f"no waterlines in {src.parent}: run the detect step first")
    mdir.mkdir(parents=True, exist_ok=True)
    built = False
    bg_day = survey_day(plan)
    if not (plan["first"] <= bg_day <= plan["last"]):
        bg_day = plan["last"]
    status = map_status(plan)
    relinked = False
    for cam in plan["cams"]:
        png = mdir / f"{plan['date']}_waterlines_{cam}.png"
        img_dir = out / "waterlines" / cam / "src"
        if plan["era"] == "live" and (out / "waterlines" / cam / "contour_points.csv").exists():
            before = sorted(q.name for q in img_dir.glob("*.jpg")) if img_dir.exists() else []
            link_frame_photos(plan, cam, args)
            relinked |= before != (sorted(q.name for q in img_dir.glob("*.jpg")) if img_dir.exists() else [])
        names = sorted(q.name for q in img_dir.glob("*.jpg")) if img_dir.exists() else []
        cmd = [sys.executable, HERE / "daily_elevation_map.py", src, img_dir, cam, png,
               "--start-date", plan["first"], "--end-date", plan["last"], "--background-date", bg_day]
        # the photos drawn on are part of the map: a photo that arrives (or goes) redraws it
        sig = {"cmd": cmd_text(cmd), "after": {"waterlines": waterline_token(plan)}, "scripts": scripts_sig("map"),
               "photos": hashlib.sha1("\n".join(names).encode()).hexdigest()}
        fresh, why = is_fresh(out, f"map_{cam}", sig, [png], [src])
        if fresh and not args.force:
            say(f"map {cam}", f"SKIPPED: {png.name} is up to date (--force to rebuild)")
            status[cam] = {"drawn": True, "photos": len(names)}
            continue
        drop_stamp(out, f"map_{cam}")
        if png.exists():
            png.unlink()
        if not names:
            why_not = (f"no photo of its frames was found in the photo roots "
                       f"({', '.join(map(str, args.photo_roots))})" if plan["era"] == "live" else
                       f"no photo in {img_dir} (the detection's photo links)")
            warn(f"map {cam}: {why_not}: photo map not drawn")
            status[cam] = {"drawn": False, "why": why_not, "kind": "no photos", "photos": 0}
            continue
        log = out / "logs" / f"map_{cam}.log"
        rc = run_cmd(cmd, log)
        if rc != 0 or not png.exists():
            # the background day may have no frame of this camera: let the script choose
            say(f"map {cam}", "retrying without --background-date")
            rc = run_cmd(cmd[:-2], log)
            if rc != 0 or not png.exists():
                why_not = (f"daily_elevation_map.py FAILED (exit {rc}; see logs/{log.name}); the photos were "
                           f"found ({len(names)} in waterlines/{cam}/src)")
                warn(f"map {cam}: {why_not}")
                status[cam] = {"drawn": False, "why": why_not, "kind": "script failed", "photos": len(names),
                               "exit": rc, "log": str(log)}
                continue
            write_stamp(out, f"map_{cam}", sig, [cmd[:-2]], scripts=STEP_SCRIPTS["map"])
        else:
            write_stamp(out, f"map_{cam}", sig, [cmd], scripts=STEP_SCRIPTS["map"])
        status[cam] = {"drawn": True, "photos": len(names)}
        built = True
    if relinked:
        link_all_photos(plan)             # the consistency plots' backdrop folder, for the next filter run
    state["map_status"] = status
    tmp = mdir / (MAP_STATUS + ".tmp")
    tmp.write_text(json.dumps(status, indent=1) + "\n")
    os.replace(str(tmp), str(mdir / MAP_STATUS))
    png = mdir / f"{plan['date']}_waterlines_plan.png"
    surveys = [s for s in plan["surveys"] if resolve_survey(s, plan, args)[0]]
    sig = {"src": str(src), "after": {"waterlines": waterline_token(plan)}, "v": 10,
           "surveys": [[s["name"], file_sig(resolve_survey(s, plan, args)[0][0])] for s in surveys],
           "scripts": scripts_sig("map_plan")}
    fresh, why = is_fresh(out, "map_plan", sig, [png], [src])
    if fresh and not args.force:
        say("plan map", f"SKIPPED: {png.name} is up to date (--force to rebuild)")
    else:
        drop_stamp(out, "map_plan")
        try:
            plan_view_map(plan, src, surveys, args, png)
            write_stamp(out, "map_plan", sig, ["survey_products.py (plan view, in-process)"],
                        scripts=STEP_SCRIPTS["map_plan"])
            say("plan map", str(png))
            built = True
        except Exception as exc:
            warn(f"plan-view map failed: {exc!r}")
    missing = [cam for cam in plan["cams"] if not (mdir / f"{plan['date']}_waterlines_{cam}.png").exists()]
    if missing:
        warn(f"waterline map(s) on the photos NOT drawn: "
             + "; ".join(f"{c} ({(status.get(c) or {}).get('kind', '?')})" for c in missing) + " (status partial)")
    return ("built" if built else "skipped") + (f" (photo maps missing: {', '.join(missing)})" if missing else "")


def map_missing_reason(plan, state, cam):
    """Why camera cam has no waterline map on the photos: this run's maps step, else the one
    recorded (maps/photo_maps.json), else what the photo folder shows."""
    st = ((state or {}).get("map_status") or map_status(plan)).get(cam) or {}
    if st.get("why") and not st.get("drawn"):
        return st["why"]
    src = plan["out"] / "waterlines" / cam / "src"
    if not src.exists() or not any(src.glob("*.jpg")):
        return "no photo of the window's frames was found in the photo roots"
    return "the maps step has not drawn it (rerun --steps maps)"


def read_ground_points(path, first, last, max_points=150000, days_out=None):
    """-> E, N, Z, camera, total points, points whose elevation includes a wave setup. days_out (a
    set) gets the capture days of the points read."""
    E, N, Z, C = [], [], [], []
    n_setup = 0
    with open(path, newline="") as f:
        rd = csv.reader(f)
        head = next(rd)
        ix = {c: i for i, c in enumerate(head)}
        ie, inn, it = ix["easting_utm19"], ix["northing_utm19"], ix["tide_elevation_navd88"]
        ib, ic, itm = ix.get("beach_elevation_navd88"), ix["camera"], ix["capture_time_utc"]
        isu = ix.get("setup_correction_m")
        for r in rd:
            if len(r) <= max(ie, inn) or not r[ie] or not (first <= r[itm][:10] <= last):
                continue
            z = r[ib] if ib is not None and r[ib] != "" else r[it]
            if isu is not None and r[isu] not in ("", "0", "0.0", "0.0000"):
                n_setup += 1
            E.append(float(r[ie])); N.append(float(r[inn])); Z.append(float(z)); C.append(r[ic])
            if days_out is not None:
                days_out.add(r[itm][:10])
    E, N, Z, C = np.array(E), np.array(N), np.array(Z), np.array(C)
    if len(E) > max_points:                  # evenly thinned for drawing
        k = np.linspace(0, len(E) - 1, max_points).astype(int)
        return E[k], N[k], Z[k], C[k], len(E), n_setup
    return E, N, Z, C, len(E), n_setup


def signed(v, fmt="{:+.1f}"):
    """A signed number with a true minus sign (U+2212): a hyphen vanishes on a busy map."""
    t = fmt.format(v)
    if t.lstrip("+-").strip("0.") == "":            # 0 shown as 0, not -0.00 / +0.00 (float noise)
        t = t.lstrip("+-")
    return t.replace("-", "−")


def plan_view_map(plan, src, surveys, args, png):
    """Filtered waterlines in plan view coloured by elevation (one sequential hue ramp, as the DEM
    page), over the survey: a DSM as thin grey contours at the same elevations, survey points as
    markers in the same colours. Where a waterline and the survey agree, the colours match."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from georectify import load_extrinsics
    import matplotlib.patheffects as pe_
    pdays = set()
    E, N, Z, C, n_all, n_setup = read_ground_points(src, plan["first"], plan["last"], days_out=pdays)
    if not len(E):
        raise StepFailed("no georectified waterline points in the window")
    cmap = LinearSegmentedColormap.from_list("sand", SAND)
    zlo, zhi = np.percentile(Z, [1, 99])
    zlo, zhi = np.floor(zlo * 4) / 4, np.ceil(zhi * 4) / 4
    norm = Normalize(zlo, zhi)
    pad = 25.0
    x0, x1, y0, y1 = E.min() - pad, E.max() + pad, N.min() - pad, N.max() + pad
    cams_xy = {}
    for cam, c in plan["cams"].items():
        if c["eo_path"]:
            eo = load_extrinsics(c["eo_path"])
            if x0 - 200 < eo[0] < x1 + 200 and y0 - 200 < eo[1] < y1 + 200:
                cams_xy[cam] = (float(eo[0]), float(eo[1]))
                x0, x1 = min(x0, eo[0] - 15), max(x1, eo[0] + 15)
                y0, y1 = min(y0, eo[1] - 15), max(y1, eo[1] + 15)
    # equal scale in E and N: size the page to the area's shape
    h_over_w = (y1 - y0) / max(x1 - x0, 1.0)
    fig_w = 9.0 if h_over_w < 1.6 else 7.5
    fig, ax = plt.subplots(figsize=(fig_w, float(np.clip(fig_w * 0.8 * h_over_w + 1.6, 5.0, 11.5))))
    legend = []
    notes = []
    dsm = next((s for s in surveys if s["survey_type"] == "dsm"), None)
    pts = next((s for s in surveys if s["survey_type"] == "points"), None)
    if dsm:
        from compare_dem_survey import read_survey
        path = resolve_survey(dsm, plan, args)[0][0]
        g, gx0, gy0, cell = read_survey(path)
        c0, c1 = max(int((x0 - gx0) / cell), 0), min(int((x1 - gx0) / cell) + 1, g.shape[1])
        r0, r1 = max(int((gy0 - y1) / cell), 0), min(int((gy0 - y0) / cell) + 1, g.shape[0])
        step = max(1, int(round(1.0 / cell)))           # 1 m is plenty for contours
        sub = g[r0:r1:step, c0:c1:step]
        if sub.size and np.isfinite(sub).any():
            xs = gx0 + (c0 + np.arange(sub.shape[1]) * step + 0.5) * cell
            ys = gy0 - (r0 + np.arange(sub.shape[0]) * step + 0.5) * cell
            # Contours at whole half-metres (a level of -0.25 m would be labelled "-0.2"),
            # in the waterlines' own colours and on top of them, outlined in ink: where
            # the two agree, a waterline lies along the contour of its own colour.
            levels = np.arange(np.ceil(zlo * 2) / 2, zhi + 1e-6, 0.5)
            halo = [pe_.Stroke(linewidth=2.6, foreground=INK), pe_.Normal()]
            with matplotlib.rc_context({"contour.negative_linestyle": "solid"}):
                cs = ax.contour(xs, ys, sub, levels=levels, cmap=cmap, norm=norm, linewidths=1.3, zorder=3)
            for coll in getattr(cs, "collections", []):          # matplotlib <= 3.7
                coll.set_path_effects(halo)
            if not getattr(cs, "collections", None):
                cs.set_path_effects(halo)                          # matplotlib >= 3.8
            # labels with a true minus and a white halo, above the waterline scatter: a hyphen in
            # dark ink vanishes on the dark lines and a -0.5 contour would read 0.5
            labels = ax.clabel(cs, fmt=lambda v: signed(v), fontsize=7.5, colors=INK, inline_spacing=2)
            for t in labels or []:
                t.set_path_effects([pe_.withStroke(linewidth=2.8, foreground="white")])
                t.set_zorder(6)
            legend.append(plt.Line2D([], [], color=SAND[4], lw=1.3, path_effects=halo,
                                     label=f"{dsm['name']} ({dsm['survey_date']}): contours every 0.5 m, "
                                           f"same colours"))
            notes.append(f"outlined lines: {dsm['name']} contours (m NAVD88)")
    ax.scatter(E, N, c=Z, cmap=cmap, norm=norm, s=0.6, lw=0, rasterized=True, zorder=2)
    # the colour range is the lines' p1..p99: the ~2% beyond it take the end colours, said by the
    # colour bar's arrows (with the survey points' own, below)
    wl_hi, wl_lo = float((Z > zhi).mean()), float((Z < zlo).mean())
    legend.append(plt.Line2D([], [], marker="o", ls="", color=SAND[2], ms=4,
                             label=f"waterlines of {min(pdays)} .. {max(pdays)} ({len(pdays)} days, {n_all:,} points)"
                             if pdays else f"waterlines ({n_all:,} points)"))
    if pts:
        from survey_compare import read_points
        paths = resolve_survey(pts, plan, args)[0]
        p = read_points(paths[0])
        inside = (p["E"] > x0) & (p["E"] < x1) & (p["N"] > y0) & (p["N"] < y1)
        pe, pn, pz = p["E"][inside], p["N"][inside], p["Z"][inside]
        # Only points inside the waterlines' elevation range can be read on its colour
        # scale; a point above it would take the end colour and look like a waterline
        # match. Those are drawn as open markers instead and counted in the legend.
        rng = (pz >= zlo) & (pz <= zhi)
        ax.scatter(pe[rng], pn[rng], c=pz[rng], cmap=cmap, norm=norm, s=34,
                   edgecolors=INK, linewidths=0.8, zorder=4)
        legend.append(plt.Line2D([], [], marker="o", ls="", mfc=SAND[4], mec=INK, ms=6,
                                 label=f"{pts['name']} ({pts['survey_date']}): {int(rng.sum())} points in "
                                       f"the colour range"))
        hi, lo = pz > zhi, pz < zlo
        if hi.any():
            ax.scatter(pe[hi], pn[hi], s=22, facecolors="white", edgecolors=MUTED, linewidths=0.8, zorder=3)
            legend.append(plt.Line2D([], [], marker="o", ls="", mfc="white", mec=MUTED, ms=5,
                                     label=f"{int(hi.sum())} survey points above {signed(zhi, '{:+.2f}')} m"))
        if lo.any():
            ax.scatter(pe[lo], pn[lo], s=22, marker="s", facecolors="white", edgecolors=MUTED,
                       linewidths=0.8, zorder=3)
            legend.append(plt.Line2D([], [], marker="s", ls="", mfc="white", mec=MUTED, ms=5,
                                     label=f"{int(lo.sum())} survey points below {signed(zlo, '{:+.2f}')} m"))
        sv_hi, sv_lo = bool(hi.any()), bool(lo.any())
        notes.append("circles: survey points, same colours (open: outside the range)")
    else:
        sv_hi = sv_lo = False
    up, down = sv_hi or wl_hi > 0, sv_lo or wl_lo > 0
    extend = "both" if (up and down) else "max" if up else "min" if down else "neither"
    if cams_xy:                     # the cameras stand metres apart: one marker, one label
        xy = np.array(list(cams_xy.values()))
        ax.plot(xy[:, 0], xy[:, 1], ls="", marker="^", color=INK, ms=7, zorder=5)
        ax.annotate(("cameras " if len(cams_xy) > 1 else "camera ") + " + ".join(cams_xy),
                    (xy[:, 0].mean(), xy[:, 1].mean()), xytext=(8, -3),
                    textcoords="offset points", color=INK, fontsize=9)
    from matplotlib.ticker import MaxNLocator, FuncFormatter
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    for a in (ax.xaxis, ax.yaxis):
        a.set_major_locator(MaxNLocator(4))
        a.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_xlabel("Easting (m, UTM 19N)", color=INK2)
    ax.set_ylabel("Northing (m, UTM 19N)", color=INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    for s in ax.spines.values():
        s.set_color("#c3c2b7")
    ax.grid(color=GRID, lw=0.5, zorder=0)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.02, extend=extend)
    # '+ wave setup' only when the points drawn carry one
    elev = ("still water + wave setup" if n_setup == n_all else
            f"still water + wave setup on {100 * n_setup / max(n_all, 1):.0f}% of points" if n_setup else
            "still water only, no wave setup")
    beyond = "; ".join(x for x in (f"{100 * wl_lo:.1f}% of line points below (end colour)" if wl_lo else "",
                                    f"{100 * wl_hi:.1f}% of line points above (end colour)" if wl_hi else "",
                                    "survey points (open markers)" if (sv_hi or sv_lo) else "") if x)
    cb.set_label(f"waterline elevation (m NAVD88, {elev})"
                 + (f"\narrows = beyond the colour range: {beyond}" if extend != "neither" else ""),
                 color=INK2)
    cb.ax.tick_params(colors=INK2, labelsize=8)
    cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: signed(v, "{:+.2f}")))
    # the legend below the map, so it never covers the ends of the contours and waterlines
    ax.legend(handles=legend, loc="upper left", bbox_to_anchor=(0.0, -0.09), fontsize=8, frameon=False,
              ncol=1)
    ax.set_title(f"{plan['date']}: filtered waterlines in plan view\n" + ("\n".join(notes) or "no survey drawn"),
                 loc="left", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(png, dpi=130, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


# ---------------------------------------------------------------------
# Step: compare, and the label
# ---------------------------------------------------------------------

def resolve_survey(s, plan, args):
    """-> (list of Paths, note). An empty path is looked for by the survey date (GCP files arrive later)."""
    if s.get("path"):
        hits = resolve_file(s["path"], args.survey_dirs)
        return hits, ("" if hits else f"{s['path']} not found in {', '.join(args.survey_dirs)}")
    if s["survey_type"] == "points" and s.get("survey_date"):
        compact = s["survey_date"].replace("-", "")
        hits = []
        for pat in (f"{s['survey_date']}*xyz*.csv", f"{s['survey_date']}*GCP*.csv", f"*{compact}*xyz*.csv"):
            for d in args.survey_dirs:
                hits += [p for p in sorted(Path(d).glob(pat)) if p not in hits]
        if hits:
            return hits, f"no path in surveys.csv; found {', '.join(str(h) for h in hits)} by its date"
    return [], "no survey file given in surveys.csv (column path)"


def missing_survey_message(s, plan, args, note):
    disabled = plan.get("forced_disabled") or not any(c["enabled"] for c in plan["cams"].values())
    if s.get("path"):
        # a file was named: say exactly which, and where it was looked for
        expect = (f"the file {s['path']} named in {args.surveys} (looked for as given, then by its name in "
                  f"{', '.join(args.survey_dirs)}; --survey-dirs adds folders), or correct its path there")
    elif s["survey_type"] == "points":
        expect = (f"{CHELSEA}/{s.get('survey_date') or plan['date']}_Marconi_Extrinsic_Targets_<cam>_xyz.csv "
                  f"(no header: num,E,N,Z in UTM 19N / NAVD88, like "
                  f"calibration/2025-11-13_Marconi_Extrinsic_Targets_c1_xyz.csv), or write its path in "
                  f"{args.surveys} (column path)")
    else:
        expect = f"the survey file (lidar DSM .tif or .asc) and write its name in {args.surveys} (column path)"
    return (f"SKIPPED {s['name']}: {note}. Provide {expect}, then rerun: python3 survey_products.py "
            f"--date {plan['date']} --steps compare" + (" --force-disabled" if disabled else ""))


def merged_points(paths, dst):
    """Several CIRN target files (e.g. one per camera) -> one, duplicates (same E,N,Z) once."""
    seen, rows = set(), []
    for p in paths:
        with open(p, newline="") as f:
            for r in csv.reader(f):
                if len(r) < 4:
                    continue
                key = tuple(r[1:4])
                if key in seen:
                    continue
                seen.add(key)
                rows.append(r[:4])
    buf = io.StringIO()
    csv.writer(buf).writerows(rows)
    if not (Path(dst).exists() and Path(dst).read_text() == buf.getvalue()):   # keep its age if unchanged
        Path(dst).write_text(buf.getvalue())
    return dst


def setup_fit_info(coef, args=None):
    """
    Where a setup coefficient came from. -> dict with coef, kind and how:
      none      C = 0, no setup
      known     a fit listed in SETUP_FITS (survey, its sha256, type and date)
      declared  --setup-fitted-to: a survey file (then recognised by content and name) or
                'none:<how>' (fitted to no survey, e.g. repeat crossings)
      unknown   neither: where it was fitted is not known
    """
    coef = float(coef or 0.0)
    if not coef:
        return {"coef": 0.0, "kind": "none", "how": "no wave-setup correction (C = 0): elevations are "
                                                    "the still-water level only"}
    known = next((dict(v) for k, v in SETUP_FITS.items() if abs(k - coef) < 1e-9), None)
    if known:
        known.pop("rtk_check", None)
        known.update(coef=coef, kind="known")
        return known
    given = (getattr(args, "setup_fitted_to", None) or "").strip() if args is not None else ""
    if given.lower().startswith("none:"):
        return {"coef": coef, "kind": "declared", "survey_stem": None,
                "how": f"{given[5:].strip()} (declared with --setup-fitted-to: fitted to no survey)"}
    if given:
        hits = resolve_file(given, getattr(args, "survey_dirs", []) or [])
        p = hits[0] if hits else None
        return {"coef": coef, "kind": "declared", "survey_stem": Path(given).stem, "survey_name": Path(given).name,
                "survey_sha256": file_sha(p) if p else None,
                "how": f"fitted to {Path(given).name} (declared with --setup-fitted-to"
                       + ("" if p else "; the file was not found, so it is recognised by name only") + ")"}
    return {"coef": coef, "kind": "unknown", "survey_stem": None,
            "how": "given on the command line without --setup-fitted-to: where it was fitted is NOT known"}


FIT_KEYS = ("coef", "kind", "survey_stem", "survey_name", "survey_sha256", "survey_type", "survey_date")


def fit_sig(fit):
    """What identifies a setup fit in a step's signature: the coefficient and the survey it was fitted
    to, not the prose describing it (rewording SETUP_FITS must not redo an hour of detection)."""
    fit = fit or {}
    return {k: fit.get(k) for k in FIT_KEYS if fit.get(k) is not None}


_FILE_SHA = {}


def file_sha(p):
    p = str(p)
    if p not in _FILE_SHA:
        _FILE_SHA[p] = sha256(p) if Path(p).is_file() else None
    return _FILE_SHA[p]


def setup_in_use(plan, args):
    """
    -> (fits, source): the setup coefficient(s) the waterlines ON DISK were built with, from the
    merge stamp (each camera's detection stamp, checked against the C the rows' own setup
    implies), each with where it was fitted and the cameras it applies to. Before anything is
    built (a dry run): this run's --setup-coef. Never the command line once waterlines exist:
    a comparison is labelled by what built the DEM it scores.
    """
    st = read_stamp(plan["out"], "merge") if plan.get("out") else None
    if st and st.get("setup"):
        fits = []
        for cam, v in sorted(st["setup"].items()):
            f = dict(v.get("fit") or setup_fit_info(v.get("coef"), None))
            if v.get("coef") is not None and abs(float(f.get("coef") or 0) - float(v["coef"])) > 1e-9:
                f = setup_fit_info(v["coef"], None)
            elif f.get("kind") in ("known", "none"):
                f = setup_fit_info(f.get("coef"), None)       # the same fit, described as the code now does
            same = next((x for x in fits if {k: x[k] for k in x if k != "cameras"} == f), None)
            if same:
                same["cameras"].append(cam)
            else:
                fits.append(dict(f, cameras=[cam]))
        return fits, "the waterlines on disk (merge stamp)"
    return [dict(setup_fit_info(args.setup_coef, args), cameras=sorted(plan["cams"]))], \
        "this run's --setup-coef (nothing built yet)"


_FIT_POINTS = {}


def fit_survey_points(fit, args):
    """(E, N) of the shots a setup coefficient was fitted to, when that file is on this computer
    (found by its name in the survey folders, with the fit's sha256 when one is known); else None."""
    names = [n for n in (fit.get("survey_name"), (fit.get("survey_stem") or "") + ".csv") if n and n != ".csv"]
    dirs = [str(x) for x in (getattr(args, "survey_dirs", None) or [])] if args is not None else []
    key = (tuple(names), fit.get("survey_sha256"), tuple(dirs))
    if key not in _FIT_POINTS:
        _FIT_POINTS[key] = None
        for n in names:
            hits = resolve_file(n, dirs)
            if hits and (not fit.get("survey_sha256") or file_sha(hits[0]) == fit["survey_sha256"]):
                try:
                    from survey_compare import read_points
                    p = read_points(hits[0])
                    _FIT_POINTS[key] = (p["E"], p["N"], hits[0].name)
                except (SystemExit, Exception):
                    pass
                break
    return _FIT_POINTS[key]


def file_dates(path):
    """The days in an Emlid file's own 'Averaging start' column (local dates), or set()."""
    try:
        with open(path, newline="") as f:
            rd = csv.DictReader(f)
            if "Averaging start" not in (rd.fieldnames or []):
                return set()
            return {r["Averaging start"][:10] for r in rd if re.match(r"\d{4}-\d\d-\d\d", r.get("Averaging start") or "")}
    except (OSError, UnicodeDecodeError, csv.Error):
        return set()


def survey_matches_fit(fit, s, survey_path, args=None):
    """Is survey s the one setup coefficient `fit` was fitted to? By content (sha256), file name,
    the shots' own coordinates (more than half within 0.05 m of the fit's shots, when that file is
    here: a re-export changes the bytes, not the points), then by type and date: the table's
    survey_date or, in an Emlid file, its own 'Averaging start' dates.
    -> ('same' | 'maybe' | None, why)."""
    if not fit or fit.get("kind") in ("none", "unknown") or not (fit.get("survey_stem") or fit.get("survey_sha256")):
        return None, ""
    p = Path(survey_path)
    sha = fit.get("survey_sha256")
    if sha and p.is_file() and (not fit.get("survey_type") or fit["survey_type"] == s.get("survey_type")):
        if file_sha(p) == sha:
            return "same", f"{p.name} has the same content (sha256) as {fit.get('survey_stem')}"
    names = {fit.get("survey_stem"), fit.get("survey_name")} - {None}
    if names & {p.stem, p.name}:
        return "same", f"file name {p.name}"
    if s.get("survey_type") == "points" and p.is_file() and fit.get("survey_type", "points") == "points":
        ref = fit_survey_points(fit, args)
        if ref is not None:
            try:
                from survey_compare import read_points
                q = read_points(p)
                d = np.full(len(q["E"]), np.inf)
                for k in range(0, len(ref[0]), 500):
                    d = np.minimum(d, np.min(np.hypot(q["E"][:, None] - ref[0][None, k:k + 500],
                                                      q["N"][:, None] - ref[1][None, k:k + 500]), axis=1))
                near = int((d <= 0.05).sum())
                if len(d) and near > 0.5 * len(d):
                    return "same", (f"{near} of its {len(d)} points lie within 0.05 m of the shots C was fitted "
                                    f"to ({ref[2]}): the same shots, re-exported")
            except (SystemExit, Exception):
                pass
    if fit.get("survey_type") and fit.get("survey_date") and s.get("survey_type") == fit["survey_type"]:
        if s.get("survey_date") == fit["survey_date"]:
            return "maybe", (f"a {fit['survey_type']} survey of {fit['survey_date']}, the type and date of the survey "
                             f"C was fitted to: a copy or re-export of it cannot be ruled out by content")
        if p.is_file() and fit["survey_date"] in file_dates(p):
            return "maybe", (f"its own 'Averaging start' dates include {fit['survey_date']}, the date of the survey "
                             f"C was fitted to (whatever surveys.csv says): a copy or re-export cannot be ruled out")
    return None, ""


class Verdict:
    """A label that can only get worse, with every reason. A 'heuristic' downgrade (a date or file-name
    rule that can misfire) can be overruled only by the survey row's label_override_reason, which is
    then printed next to it; a downgrade from what the chain actually did cannot."""

    def __init__(self, s):
        self.table = s["label"]
        self.label = s["label"] if s["label"] in LABELS else "CIRCULAR"
        self.reasons = []
        self.overridden = []
        self.override = (s.get("label_override_reason") or "").strip()
        if s["label"] not in LABELS:
            self.reasons.append(f"CIRCULAR: label {s['label']!r} is not one of {', '.join(LABELS)}")

    def worse(self, new, why, heuristic=False):
        if heuristic and self.override and LABELS.index(new) > LABELS.index(self.label):
            self.overridden.append(f"{new}: {why}")
            return
        self.reasons.append(f"{new}: {why}")
        if LABELS.index(new) > LABELS.index(self.label):
            self.label = new


def honest_label(s, plan, args, survey_path, fits=None, verdict=None):
    """The table's label, made worse (never better) by what the chain actually used.
    -> (label, reasons); with verdict=Verdict(s) given, that object is filled (overrides too)."""
    v = verdict or Verdict(s)
    sname, sstem = Path(survey_path).name, Path(survey_path).stem
    for cam, c in plan["cams"].items():
        env = c.get("envelope_path")
        if env and (Path(env).name == sname or Path(env).stem == sstem or
                    (Path(env).is_file() and Path(survey_path).is_file() and file_sha(env) == file_sha(survey_path))):
            v.worse("PARTLY-CIRCULAR", f"the {cam} search envelope was placed with this survey ({Path(env).name})")
        eo = c.get("eo_path")
        if eo and Path(eo).exists():
            # A pointing fitted to a survey (fit_eo_to_survey.py, *_lidar_EO.yaml) or carrying such a
            # fit's change to another period (apply_pointing_correction.py, *_corr_EO.yaml): the
            # lineage is followed to the survey each fit names, recognised by name or content.
            from survey_compare import pointing_lineage, same_survey, lineage_text
            look = eo_dirs(args)
            lin = pointing_lineage(eo, look)
            via = f" (lineage {lineage_text(lin)})" if len(lin["chain"]) > 1 else ""
            hit = [f for f in lin["fitted_to"] if same_survey(f, survey_path, look)]
            named = sname in lin["notes"] or (len(sstem) > 6 and sstem in lin["notes"])
            if hit or (lin["survey_fit"] and named):
                v.worse("CIRCULAR", f"the {cam} pointing {Path(eo).name} "
                                    + ("carries a correction fitted to" if len(lin["chain"]) > 1 else "was fitted to")
                                    + f" this survey{via}")
            elif named:
                v.worse("CIRCULAR", f"the {cam} calibration {Path(eo).name} names this survey in its notes "
                                    f"(was it fitted to it?)", heuristic=True)
            elif lin["survey_fit"] and lin["unresolved"]:
                v.worse("CIRCULAR", f"the {cam} pointing {Path(eo).name} is survey-derived{via} and which survey "
                                    f"cannot be resolved ({'; '.join(lin['unresolved'])}): it may be this one",
                        heuristic=True)
        # any point survey made on the calibration's own day may be the GCPs it was solved from,
        # whatever the file is called
        if s["survey_type"] == "points" and eo_date(c["eo_file"]) == s.get("survey_date"):
            v.worse("PARTLY-CIRCULAR", f"the {cam} calibration {c['eo_file']} was solved on the survey day: "
                                       f"these points may be the GCPs it was solved from", heuristic=True)
    if fits is None:
        fits, _ = setup_in_use(plan, args)
    for f in fits:
        cams = f"; cameras {', '.join(f['cameras'])}" if f.get("cameras") and len(plan["cams"]) > 1 else ""
        if f.get("kind") == "unknown":
            v.worse("PARTLY-CIRCULAR", f"setup coefficient C = {f['coef']} of unknown origin (no "
                                       f"--setup-fitted-to): it may have been fitted to this survey{cams}")
            continue
        how, why = survey_matches_fit(f, s, survey_path, args)
        if how:
            v.worse("CIRCULAR", f"the setup coefficient C = {f['coef']} was fitted to this survey ({why}{cams}), "
                                f"so this comparison cannot test the setup or the overall level; its spread "
                                f"over the elevations the survey covers still says something",
                    heuristic=(how == "maybe"))
    return v.label, v.reasons


def eo_dirs(args):
    """Folders where the files an EO's notes name (the fit a *_corr_EO.yaml was carried from) are
    looked for: the calibration folder given, the repository's, and the survey folders."""
    dirs = [str(getattr(args, "calibration", None) or HERE / "calibration"), str(HERE / "calibration")]
    dirs += [str(d) for d in (getattr(args, "survey_dirs", None) or [])]
    return list(dict.fromkeys(dirs))


def survey_compare_findings(plan, s, survey_path, label, dem, args=None):
    """survey_compare.py's own label checks (audit_label), run here BEFORE the comparison so the
    headline label already includes them. -> list of (suggested label or None, text, kind), kind
    'chain' (what the chain did: never overruled), 'rule' (a date or name rule that can misfire:
    overruled only by surveys.csv label_override_reason) or 'note'."""
    from survey_compare import audit_label, load_cameras
    eos = {cam: Path(c["eo_path"]) for cam, c in plan["cams"].items()
           if c.get("eo_path") and (plan["out"] / "waterlines" / cam / "contour_points_ground.csv").exists()}
    envs = sorted({str(c["envelope_path"]) for c in plan["cams"].values() if c.get("envelope_path")})
    try:
        cams = load_cameras(eos, {}) if eos else {}
        return audit_label(label, survey_path, s["survey_type"], s.get("survey_date"), cams,
                           envs[0] if len(envs) == 1 else None, dem, dirs=eo_dirs(args) if args else ())
    except SystemExit as exc:                      # its loaders exit on a missing file
        return [(None, f"survey_compare.py's checks could not run: {exc}", "note")]


def frames_used_dates(plan):
    """(first, last) day of the frames the DEM gridded (dem_from_contours.py's info), else of the
    waterlines in the window, else the window."""
    info = load_json(plan["out"] / "dem" / f"{plan['date']}_info.json") or {}
    if info.get("first_date") and info.get("last_date"):
        return info["first_date"], info["last_date"], "the frames the DEM gridded"
    days = sorted({d for v in frames_per_day(waterline_file(plan, quiet_=True)).values() for d in v
                   if plan["first"] <= d <= plan["last"]})
    if days:
        return days[0], days[-1], "the waterline frames in the window"
    return plan["first"], plan["last"], "the configured window (no frame found)"


def step_compare(plan, args, state):
    out = plan["out"]
    cdir = out / "compare"
    stem = out / "dem" / plan["date"]
    dem = Path(f"{stem}_dem.asc")
    if not plan["surveys"]:
        say("compare", f"no survey listed for {plan['date']} in {args.surveys}: nothing to compare")
        return "skipped"
    results = state.setdefault("comparisons", {})
    todo = []
    for s in plan["surveys"]:
        paths, note = resolve_survey(s, plan, args)
        if paths:
            todo.append((s, paths, note))
            continue
        msg = missing_survey_message(s, plan, args, note)
        say("compare", msg)
        for p in comparison_files(plan, s):              # nothing of an earlier file may look current
            if p.exists():
                p.unlink()
        drop_stamp(out, f"compare_{s['name']}")
        results[s["name"]] = {"status": "skipped", "reason": msg, "label": s["label"]}
    if not todo:
        return "skipped (no survey file)"
    if not dem.exists() or not token(out, "dem"):
        raise StepFailed(f"no current DEM {dem}: run the dem step first")
    write_provenance(plan, args, state, final=False)     # survey_compare.py audits it
    fits, fits_from = setup_in_use(plan, args)
    used = {round(float(f["coef"]), 6) for f in fits}
    if fits_from.startswith("the waterlines") and used != {round(float(args.setup_coef), 6)}:
        warn(f"this run's --setup-coef {args.setup_coef} is NOT what the waterlines were built with "
             f"(C = {', '.join(str(c) for c in sorted(used))}): the label, README and provenance use the C that "
             f"built them. To apply {args.setup_coef}, rebuild from the detection (--steps detect,filter,dem,"
             f"maps,compare or all steps)")
    dd = dem_settings_differ(dem_settings_built(out), dem_settings_requested(args))
    if dd:
        warn(f"this run's DEM options are NOT what the DEM on disk was built with ({dd}): the comparison, README "
             f"and provenance use the DEM on disk and its settings. To apply them, rebuild: --steps dem,maps,compare")
    first_used, last_used, used_how = frames_used_dates(plan)
    built = False
    failed_any = 0
    for s, paths, note in todo:
        name = f"{plan['date']}_{s['name']}"
        if note:
            say("compare", note)
        survey_path = paths[0]
        if len(paths) > 1:
            cdir.mkdir(parents=True, exist_ok=True)
            survey_path = merged_points(paths, cdir / f"{name}_survey_points.csv")
            say("compare", f"{len(paths)} files merged into {survey_path}")
        v = Verdict(s)
        # every file of the survey (several point files are merged): the worst of them
        for p_ in paths:
            honest_label(s, plan, args, p_, fits=fits, verdict=v)
        for sug, txt, kind in survey_compare_findings(plan, s, survey_path, v.label, dem, args):
            if sug and LABELS.index(sug) > LABELS.index(v.label):
                v.worse(sug, f"survey_compare.py check: {txt}", heuristic=(kind == "rule"))
        label = v.label
        why = s["why"] + ("" if not v.reasons else " | chain checks -> " + "; ".join(v.reasons))
        if v.overridden:
            why += (f" | {len(v.overridden)} downgrade(s) to "
                    + "/".join(sorted({o.split(':')[0] for o in v.overridden}))
                    + f" overruled in surveys.csv: {v.override}")
        if label != s["label"]:
            warn(f"{s['name']}: labelled {label}, not {s['label']} as in {Path(args.surveys).name}: "
                 + "; ".join(v.reasons))
        for o in v.overridden:
            warn(f"{s['name']}: downgrade NOT applied, overruled in {Path(args.surveys).name} "
                 f"({v.override}): {o}")
        envs = sorted({str(c["envelope_path"]) for c in plan["cams"].values() if c.get("envelope_path")})
        cmd = [sys.executable, HERE / "survey_compare.py", "--dem", dem, "--survey", survey_path,
               "--survey-type", s["survey_type"], "--name", name, "--output-dir", cdir,
               "--label", label, "--why", why, "--photo-dates", first_used, last_used,
               "--contours", waterline_file(plan),
               "--spread", f"{stem}_spread.asc", "--count", f"{stem}_count.asc",
               # the difference GeoTIFF carries the same tag as the DEM's
               "--epsg", (dem_settings_built(out) or {}).get("geotiff_epsg") or args.geotiff_epsg]
        if s.get("survey_date"):
            cmd += ["--survey-date", s["survey_date"]]
        eos = [f"{cam}={c['eo_path']}" for cam, c in plan["cams"].items() if c["eo_path"] and
               (out / "waterlines" / cam / "contour_points_ground.csv").exists()]
        if eos:
            cmd += ["--camera-eo"] + eos + ["--eo-dirs"] + eo_dirs(args)
        if len(envs) == 1:
            cmd += ["--envelope-source", envs[0]]
        syn = product_synthetic(plan)
        if syn:
            cmd += ["--synthetic", syn]
        sig = {"cmd": cmd_text(cmd), "after": {"dem": token(out, "dem"), "waterlines": waterline_token(plan)},
               "survey": file_sig(survey_path), "scripts": scripts_sig("compare")}
        js = cdir / f"{name}_comparison.json"
        fresh, why_f = is_fresh(out, f"compare_{s['name']}", sig, [js], [dem, survey_path, waterline_file(plan)])
        info = {"label": label, "table_label": s["label"], "reasons": v.reasons, "overridden": v.overridden,
                "override_reason": v.override or None, "photo_dates_used": [first_used, last_used, used_how]}
        if fresh and not args.force:
            # the label recorded when it ran may be worse (survey_compare.py's checks after the run)
            old = (read_stamp(out, f"compare_{s['name']}") or {}).get("label") or {}
            if old.get("label") in LABELS and LABELS.index(old["label"]) > LABELS.index(info["label"]):
                info = dict(old, photo_dates_used=info["photo_dates_used"])
            say("compare", f"SKIPPED {name}: up to date (--force to rebuild); {info['label']}")
            results[s["name"]] = dict(info, status="skipped (up to date)", json=str(js))
            continue
        drop_stamp(out, f"compare_{s['name']}")
        for p in comparison_files(plan, s):        # a failed run must not leave the last one looking current
            if p.exists():
                p.unlink()
        say("compare", f"{name}: {label} ({why_f or '--force'}); expected 15-40 s on the NUC for a lidar")
        log = out / "logs" / f"compare_{s['name']}.log"
        rc = run_cmd(cmd, log)
        if rc != 0 or not js.exists():
            warn(f"comparison {name} FAILED (exit {rc}); see {log}")
            results[s["name"]] = dict(info, status="failed")
            state["failures"].append(f"compare {s['name']}: survey_compare.py failed (exit {rc}); see {log}")
            failed_any += 1
            continue
        head = load_json(js) or {}
        sug = head.get("suggested_label")
        if sug in LABELS and LABELS.index(sug) > LABELS.index(label):
            # its checks still found something worse than the checks run before it (they should be
            # the same): the worse label wins (a rule that can misfire only if not overruled), and
            # the page is drawn again with it
            for f_ in head.get("findings") or [{"suggested": sug, "text": "; ".join(head.get("warnings") or []),
                                                 "kind": "rule"}]:
                fs = f_.get("suggested")
                if fs in LABELS and LABELS.index(fs) > LABELS.index(label) and \
                        not any(o.startswith(f"{fs}: survey_compare.py check") for o in v.overridden):
                    v.worse(fs, "survey_compare.py check after the run: " + str(f_.get("text")),
                            heuristic=(f_.get("kind") != "chain"))
            info.update(label=v.label, reasons=v.reasons, overridden=v.overridden)
            if v.label != label:
                warn(f"{s['name']}: survey_compare.py's checks make it {v.label}: drawn again with that label")
                why = why + " | survey_compare.py after the run -> " + v.reasons[-1]
                k = cmd.index("--label")
                cmd[k + 1], cmd[cmd.index("--why") + 1] = v.label, why
                rc = run_cmd(cmd, log)
                if rc != 0 or not js.exists():
                    results[s["name"]] = dict(info, status="failed")
                    state["failures"].append(f"compare {s['name']}: survey_compare.py failed (exit {rc}); see {log}")
                    failed_any += 1
                    continue
        write_stamp(out, f"compare_{s['name']}", sig, [cmd], scripts=STEP_SCRIPTS["compare"], extra={"label": info})
        results[s["name"]] = dict(info, status="built", json=str(js))
        built = True
    if failed_any:
        return f"built ({failed_any} failed)" if built else "failed"
    return "built" if built else "skipped"


# ---------------------------------------------------------------------
# Provenance and README
# ---------------------------------------------------------------------

SCRIPTS = ["survey_products.py", "historical_forcing.py", "detect_original_view.py", "waterline_detector_v5.py",
           "extract_elevation_contours.py", "georectify.py", "waterline_consistency.py", "dem_from_contours.py",
           "dem_figure.py", "asc_to_geotiff.py", "daily_elevation_map.py", "survey_compare.py",
           "pointing_check.py", "horizon_check.py", "estimate_eo_rotation.py", "view_reproject.py",
           "compare_dem_survey.py", "compare_rtk.py", "marconi_water_level.py"]

_SOFTWARE = {}


def software():
    """The code present now (commit, scripts not as committed, sha256s); computed once per run."""
    if not _SOFTWARE:
        _SOFTWARE.update(_software_now())
    return _SOFTWARE


def _software_now():
    def git(*a):
        try:
            r = subprocess.run(["git", "-C", str(HERE)] + list(a), stdout=subprocess.PIPE,
                               stderr=subprocess.DEVNULL, universal_newlines=True, timeout=20)
            return r.stdout if r.returncode == 0 else None
        except Exception:
            return None
    commit = (git("rev-parse", "HEAD") or "").strip() or None
    status = git("status", "--porcelain", "--", *SCRIPTS) or ""
    changed = [ln[3:] for ln in status.splitlines() if ln.strip()]
    return {"commit": commit, "branch": (git("rev-parse", "--abbrev-ref", "HEAD") or "").strip() or None,
            "scripts_not_as_committed": changed,
            "note": ("the scripts listed in scripts_not_as_committed differ from the commit (modified or "
                     "untracked): their sha256 below identifies them" if changed else
                     "every script used is as committed"),
            "script_sha256": {s: script_sha(s) for s in SCRIPTS if (HERE / s).exists()},
            "python": sys.version.split()[0]}


def software_of_outputs(out):
    """What built the outputs on disk, from their stamps: per step, the commit and the sha256
    of each script it ran; and whether that is the code present now."""
    now = software()
    steps, commits, dirty, differs = {}, set(), set(), set()
    ran = {}                                      # script -> set of sha256s it ran with
    sdir = Path(out) / "logs" / "stamps"
    for p in sorted(sdir.glob("*.json")) if sdir.exists() else []:
        sw = (load_json(p) or {}).get("software")
        if not sw:
            steps[p.stem] = {"commit": None, "note": "built before the stamps recorded the software"}
            continue
        steps[p.stem] = sw
        commits.add(sw.get("commit"))
        dirty.update(sw.get("not_as_committed") or [])
        for s, h in (sw.get("scripts") or {}).items():
            ran.setdefault(s, set()).add(h)
            if h and now["script_sha256"].get(s) and h != now["script_sha256"][s]:
                differs.add(s)
    old = [k for k, v in steps.items() if v.get("commit") is None]
    # a script uncommitted when it ran, committed since without a change: say which commit holds it
    same_now = sorted(s for s in dirty if s not in (now.get("scripts_not_as_committed") or [])
                      and ran.get(s) == {now["script_sha256"].get(s)})
    if not steps:
        note = "nothing built yet"
    elif len(commits) == 1 and not dirty and not old:
        note = f"every output was built at commit {next(iter(commits))} with the scripts as committed"
    else:
        note = ("outputs built at " + (f"{len(commits)} commits ({', '.join(sorted(str(c)[:10] for c in commits))})"
                                       if len(commits) > 1 else f"commit {str(next(iter(commits)))[:10]}")
                + (f"; scripts not as committed when they ran: {', '.join(sorted(dirty))} (their sha256 is in "
                   f"each step's entry)" if dirty else "")
                + (f"; of those, {', '.join(sorted(same_now))} are byte for byte the files committed in "
                   f"{str(now.get('commit'))[:10]}" if same_now else "")
                + (f"; {len(old)} step(s) built before the software was recorded" if old else ""))
    ext = sorted(differs - {"survey_products.py"})
    if ext:
        note += (f". The code present now differs from what built them in: {', '.join(ext)} "
                 f"(a rerun rebuilds the steps that use them)")
    if "survey_products.py" in differs:
        note += (". survey_products.py has changed since some steps ran (its own in-process steps rebuild "
                 "when their version number changes)")
    return {"by_step": steps, "commits": sorted(str(c) for c in commits if c), "note": note,
            "scripts_changed_since": sorted(differs),
            "now": {"commit": now.get("commit"), "branch": now.get("branch"),
                    "scripts_not_as_committed": now.get("scripts_not_as_committed"), "python": now.get("python")}}


def calibration_history(eo_file, cal_dir):
    try:
        with open(Path(cal_dir) / "calibration_history.csv", newline="") as f:
            for r in csv.DictReader(f):
                if r.get("file") == Path(str(eo_file)).name:
                    return r
    except OSError:
        pass
    return None


def load_json(p):
    try:
        return json.loads(Path(p).read_text())
    except (OSError, ValueError):
        return None


def photo_section(plan, cam):
    pt = plan["out"] / "waterlines" / cam / "photos.txt"
    if not pt.exists():
        return None
    paths = [ln for ln in pt.read_text().splitlines() if ln.strip()]
    names = [Path(p).name for p in paths]
    days = {}
    for n in names:
        t = utc_of(n)
        if t:
            days[t.strftime("%Y-%m-%d")] = days.get(t.strftime("%Y-%m-%d"), 0) + 1
    return {"count": len(names), "list_file": str(pt), "per_day": days,
            "folders": sorted({str(Path(p).parent) for p in paths}), "names": names}


def frames_per_day(path):
    """{camera: {day: frames}} of a waterline file (the one the DEM used), or {} if none yet."""
    out = {}
    if not Path(path).exists():
        return out
    seen = set()
    with open(path, newline="") as f:
        rd = csv.reader(f)
        head = next(rd, [])
        ix = {c: i for i, c in enumerate(head)}
        if not all(k in ix for k in ("source_file", "camera", "capture_time_utc")):
            return out
        i_s, i_c, i_t = ix["source_file"], ix["camera"], ix["capture_time_utc"]
        for r in rd:
            if len(r) <= max(i_s, i_c, i_t) or r[i_s] in seen:
                continue
            seen.add(r[i_s])
            d = out.setdefault(r[i_c], {})
            d[r[i_t][:10]] = d.get(r[i_t][:10], 0) + 1
    return {c: dict(sorted(v.items())) for c, v in sorted(out.items())}


def failed_cameras(plan, state):
    """{camera: why} of the enabled cameras with no waterlines in the merged file: this run's
    failures, else those recorded when the merged file was built (a later --steps compare must
    not call a one-camera product complete)."""
    failed = dict(state.get("failed_cams") or {})
    m = read_stamp(plan["out"], "merge")
    if m and m.get("cameras") is not None:
        for cam in plan["cams"]:
            if cam not in m["cameras"] and cam not in failed:
                failed[cam] = (m.get("failed_cams") or {}).get(cam) or "no waterlines in the merged file"
    return failed


def missing_parts(plan, state, args=None):
    """The parts a complete product has that are missing, per camera that has waterlines: the
    pointing check (it needs the photos) and the waterline map on the photos; and per survey whose
    file surveys.csv names (or that was found by its date), a current comparison. On the station a
    wrong or empty photo root, or a misplaced or misnamed survey file (or a wrong --survey-dirs),
    would otherwise look like a finished product that validated nothing. A survey row with no file
    yet (e.g. the Oct 2024 GCPs, not delivered) is not a missing part. -> list of text."""
    out = plan["out"]
    failed = failed_cameras(plan, state)
    miss = []
    for s in plan.get("surveys") or []:
        found = None
        if args is not None and getattr(args, "survey_dirs", None) is not None:
            found = bool(resolve_survey(s, plan, args)[0])
        r = (state.get("comparisons") or {}).get(s["name"]) or {}
        if s.get("path") and found is False:
            miss.append(f"{s['name']}: survey file {s['path']} NOT FOUND, not compared (surveys.csv names it; "
                        f"looked in {', '.join(args.survey_dirs)}; --survey-dirs adds folders)")
        elif (s.get("path") or found) and not comparison_current(plan, s)[0] and not str(r.get("status", "")).startswith("failed"):
            miss.append(f"{s['name']}: no current comparison ({comparison_current(plan, s)[1]}; rerun --steps compare)")
    for cam in plan["cams"]:
        if cam in failed:
            continue                              # already a reason of its own
        res = (state.get("pointing") or {}).get(cam) or load_json(out / "pointing" / f"pointing_{cam}.json")
        if res is None:
            miss.append(f"{cam}: pointing not checked (the pointing step has not run)")
        elif not res.get("rows"):
            miss.append(f"{cam}: pointing NOT checked (no photo of the window in the photo roots)")
        if not (out / "maps" / f"{plan['date']}_waterlines_{cam}.png").exists():
            miss.append(f"{cam}: no waterline map on the photos ({map_missing_reason(plan, state, cam)})")
    return miss


def chain_current(plan, maps=True):
    """Do the outputs on disk form ONE build? Each stamp names the upstream stamps it was built
    after (its token); a stamp dropped by a failed or interrupted rebuild breaks the chain.
    -> (True, '') or (False, why)."""
    out = plan["out"]
    if not token(out, "forcing"):
        return False, "the forcing step has no finished build (it failed, or has not run, with these settings)"
    m = read_stamp(out, "merge")
    if not m:
        return False, "the waterlines have no finished build"
    keys = [("merge", m)]
    for k in (m.get("sig") or {}).get("after") or {}:
        if k.startswith("detect_"):
            keys.append((k, read_stamp(out, k) or {}))
    keys.append(("forcing", read_stamp(out, "forcing") or {}))
    d = read_stamp(out, "dem")
    if not d:
        return False, "the DEM has no finished build"
    keys.append(("dem", d))
    # the grids on disk are the ones that stamp's settings built (dem_from_contours.py writes its own
    # settings into <date>_info.json): a grid rebuilt by hand or by another run is not this build
    built = dem_settings_built(out) or {}
    info = load_json(out / "dem" / f"{plan['date']}_info.json") or {}
    if info and built:
        filt = info.get("filters") or {}
        on_disk = {"cell_m": info.get("cell"), "min_points": info.get("min_points"),
                   "max_spread_m": info.get("max_spread"), "max_hs_m": filt.get("max_hs"),
                   "max_day_offset_m": filt.get("max_day_offset")}
        for k, v in on_disk.items():
            b = built.get(k)
            if k == "max_day_offset_m" and built.get("day_test_off"):
                b = None                     # asked for, but not run on so few camera-days (step_dem)
            if (b is None) != (v is None) or (b is not None and abs(float(b) - float(v)) > 1e-9):
                return False, (f"the DEM grids on disk ({plan['date']}_info.json: {k} {v}) are not the ones the DEM "
                               f"stamp built ({k} {b}); rerun --steps dem,maps,compare")
    for k, st in keys:
        for up, tk in ((st.get("sig") or {}).get("after") or {}).items():
            if token(out, up) != tk:
                return False, f"{k} was built from another {up} build than the one on disk"
    # the maps drawn from the waterlines: a map of earlier waterlines is not part of this build
    wt = waterline_token(plan)
    sdir = out / "logs" / "stamps"
    for p in sorted(sdir.glob("map_*.json")) if (maps and sdir.exists()) else []:
        st = load_json(p) or {}
        if ((st.get("sig") or {}).get("after") or {}).get("waterlines") != wt:
            return False, f"{p.stem} was drawn from earlier waterlines (rerun --steps maps)"
    return True, ""


def comparison_current(plan, s):
    """Is the comparison on disk for survey s one of the DEM and waterlines on disk now?"""
    st = read_stamp(plan["out"], f"compare_{s['name']}")
    if not st:
        return False, "no finished comparison"
    after = (st.get("sig") or {}).get("after") or {}
    if after.get("dem") != token(plan["out"], "dem"):
        return False, "made with another DEM"
    if after.get("waterlines") != waterline_token(plan):
        return False, "made with other waterlines"
    ok, why = chain_current(plan, maps=False)          # a comparison does not use the maps
    return ok, why


def tilt_sensitivity(plan, cam, dtilt, droll):
    """
    The DEM error a pointing error of (dtilt, droll) deg would cause, per range band, from this
    date's own DEM and calibration. Each DEM cell is projected into the photo with the
    calibration, mapped back to the ground at its own elevation with the pointing changed, and
    the DEM's local slope times that horizontal shift is the elevation error. The view grazes the
    beach, so a small tilt moves a line metres along the ground (~1/sin of the depression angle).
    -> list of dicts (range band, cells, median shift m, median error m), or None without a DEM.
    """
    from georectify import load_intrinsics, load_extrinsics, pixel_to_ground
    from view_reproject import ground_to_pixel
    from asc_to_geotiff import read_asc
    c = plan["cams"].get(cam) or {}
    demp = plan["out"] / "dem" / f"{plan['date']}_dem.asc"
    if not demp.exists() or not c.get("eo_path") or not Path(c["io_path"]).exists():
        return None
    g, h = read_asc(demp)
    cell = h["cellsize"]
    xll, ytop = h["xllcorner"], h["yllcorner"] + g.shape[0] * cell
    d_row, d_col = np.gradient(g, cell)
    dzde, dzdn = d_col, -d_row                   # rows run north to south
    rr, cc = np.nonzero(np.isfinite(g) & np.isfinite(dzde) & np.isfinite(dzdn))
    if len(rr) < 20:
        return None
    E, N, Z = xll + (cc + 0.5) * cell, ytop - (rr + 0.5) * cell, g[rr, cc]
    io, eo = load_intrinsics(c["io_path"]), load_extrinsics(c["eo_path"])
    U, V, ok = ground_to_pixel(E, N, Z, io, eo)
    if ok.sum() < 20:
        return None
    e2 = eo.copy()
    e2[4] += np.deg2rad(dtilt)
    e2[5] += np.deg2rad(droll)
    X2, Y2 = pixel_to_ground(U[ok], V[ok], Z[ok], io, e2)
    dE, dN = X2 - E[ok], Y2 - N[ok]
    err = dzde[rr[ok], cc[ok]] * dE + dzdn[rr[ok], cc[ok]] * dN
    shift = np.hypot(dE, dN)
    rng = np.hypot(E[ok] - eo[0], N[ok] - eo[1])
    rows = []
    for a, b in ((0, 100), (100, 150), (150, 200), (200, 250), (250, 300), (300, 350), (350, 450), (450, 800)):
        m = (rng >= a) & (rng < b) & np.isfinite(err)
        if m.sum() >= 15:
            rows.append({"range_m": f"{a}-{b}", "cells": int(m.sum()),
                         "shift_m": round(float(np.median(shift[m])), 2),
                         "dem_error_m": round(float(np.median(err[m])), 3)})
    return rows or None


def survey_path_for(s, plan, args):
    """The survey file of row s as found (or its name when it is not here): for the label checks."""
    try:
        hits = resolve_survey(s, plan, args)[0] if getattr(args, "survey_dirs", None) is not None else []
    except Exception:
        hits = []
    return hits[0] if hits else Path(s.get("path") or s["name"])


def prior_checks(plan):
    """The PRIOR_CHECKS entries for this date: same date, one of its surveys, the same calibration
    for every camera the entry measured that this date builds."""
    out = []
    stems = {Path(s.get("path") or "").stem for s in plan.get("surveys") or []}
    for pc in PRIOR_CHECKS:
        if pc["date"] != plan.get("date") or pc["survey_stem"] not in stems:
            continue
        cams = [c for c in pc["eo"] if c in plan["cams"]]
        if cams and all(plan["cams"][c].get("eo_file") == pc["eo"][c] for c in cams):
            out.append(dict(pc, cameras=cams))
    return out


def build_setup_effect(prov, survey_stems):
    """This build's own paired C = 0 sensitivity: per camera, the median rise the setup gave the lines
    once re-projected (with - without, the same frames or points), from the current comparison(s) with
    the given survey stem(s) on disk. A camera's value rests on at least FEW_N frames or survey points
    (else it is not used); of several surveys, the one with the most is used. -> {camera: (effect m,
    what it rests on)} ({} when there is none)."""
    if isinstance(survey_stems, str):
        survey_stems = [survey_stems]
    best = {}
    for c in (prov or {}).get("comparisons") or []:
        if not c.get("current") or Path(c.get("survey_path") or "").stem not in survey_stems:
            continue
        pr = ((((c.get("headline") or {}).get("waterlines") or {}).get("without_setup") or {}).get("paired") or {})
        for cam, v in (pr.get("by_camera") or {}).items():
            if v.get("setup_effect_median") is None or not v.get("n"):
                continue
            k, what = independent_count(v)
            if k >= FEW_N and k > best.get(cam, (None, -1))[1]:
                best[cam] = (float(v["setup_effect_median"]), k, what)
    return {cam: (e, f"{k} {what}") for cam, (e, k, what) in best.items()}


def prior_check_caveat(plan, pc, fits, forcing, prov=None):
    """The README caveat for an earlier check of the real photos: what its tool printed, the same
    numbers as waterline - survey (per camera, signs written out), and what this build should then
    give with its setup: the earlier lines plus the rise THIS build's setup gives each camera's lines
    once re-projected (its own paired C = 0 sensitivity, per camera), else a range over the shares
    kept measured so far (SHARE_KEPT_RANGE x the window's setup). The share differs by date, camera
    and geometry (2026 RTK 0.71, Jan 2025 ~0.6, Mar 2025 ~0.8), so no single share is carried over."""
    wv = (forcing or {}).get("waves", {})
    su = wv.get("setup_in_window") or {}
    med = su.get("daytime_median_m") if su.get("daytime_median_m") is not None else su.get("median_m")
    setup = med if med is not None else pc["predicted_setup_m"]
    C = fits[0]["coef"] if fits else None
    eff = build_setup_effect(prov, pc["survey_stem"]) if C else {}
    cams = sorted(pc["cameras"])
    per, rep = pc["waterline_minus_survey"], pc.get("reported") or {}
    low = [c for c in cams if per[c][0] < -0.05]
    high = [c for c in cams if per[c][0] > 0.05]
    wms = ", ".join(f"{c} {per[c][0]:+.2f} m (NMAD {per[c][1]:.2f})" for c in cams)
    text = (f"EARLIER CHECK (real photos)" + (f", itself {pc['label']}: {pc['label_why']}." if pc.get("label") else ":")
            + f" Before this script, {pc.get('tool', 'an earlier check')} on the real "
            f"photos of {pc['window']} with these same calibrations ({', '.join(pc['eo'][c] for c in pc['cameras'])}), "
            f"{pc['how']}, ")
    if rep and pc.get("reported_as"):
        text += (f"printed {pc['reported_as']} of "
                 + ", ".join(f"{c} {rep[c][0]:+.2f} m (NMAD {rep[c][1]:.2f})" for c in cams if c in rep)
                 + f" ({pc.get('reported_meaning', '')}). As waterline - {pc['survey_text']} (the sign of DEM - "
                   f"survey in this README) that is {wms}. ")
    else:
        text += f"gave waterline - {pc['survey_text']} of {wms}. "
    nosetup = not pc.get("setup_coef")
    if low and not high and nosetup:
        text += ("LOW is the sign a line given no setup must have (a timex line is marked where the swash reaches, "
                 "above still water), the sign of the biases listed below and the sign this build's C = 0 "
                 "sensitivity should show: there is no opposite-sign discrepancy. ")
    elif high and not low and nosetup:
        text += ("HIGH without any setup is NOT what the swash explains (it would make the lines read low): the "
                 "setup will make them read higher still, so look for a pointing or water-level error. ")
    if C and nosetup:
        lo_k, hi_k = SHARE_KEPT_RANGE
        rng = ", ".join(f"{c} {per[c][0] + setup * lo_k:+.2f} to {per[c][0] + setup * hi_k:+.2f}" for c in cams)
        full = ", ".join(f"{c} {per[c][0] + setup:+.2f}" for c in cams)
        if all(c in eff for c in cams):
            own = ", ".join(f"{c} ~{per[c][0] + eff[c][0]:+.2f} m ({per[c][0]:+.2f} {eff[c][0]:+.3f}, on "
                            f"{eff[c][1]})" for c in cams)
            how_ = (f"this build's own paired C = 0 sensitivity gives each camera's lines a rise of "
                    + ", ".join(f"{c} {eff[c][0]:+.3f} m ({eff[c][0] / setup:.2f} x the setup median)" for c in cams)
                    + f": expect waterline - {pc['survey_text']} of about {own}")
        else:
            how_ = (f"once its landward re-projection is counted a line keeps ~{lo_k}-{hi_k} x its setup (the shares "
                    f"measured so far: it depends on the date, the camera and the beach; this build's own paired C = 0 "
                    f"sensitivity, in compare/*_comparison.txt, gives it per camera): expect waterline - "
                    f"{pc['survey_text']} of about {rng} m")
        text += (f"This build adds the setup (C = {C}: median {setup:.2f} m in the window's daytime"
                 f"{'' if med is not None else ', predicted'}), which raises every line; {how_} unless this build's "
                 f"envelope changes the detections ({rng} m over the shares {lo_k}-{hi_k}; {full} m if the whole "
                 f"setup counted; the earlier check's own NMAD, {', '.join(f'{c} {per[c][1]:.2f}' for c in cams)} m, "
                 f"is the scale of agreement to expect). That expectation already holds the setup's "
                 f"under-correction listed below: do not subtract it again"
                 + (". It holds the still-water-reference effect listed below too: the earlier check was made on "
                    "lines at the same ADCP still water as this build's, so that effect is inside its numbers. Those "
                    "caveats say why the lines may differ from the lidar, not from this expectation. "
                    if plan.get("era") in ("adcp", "chatham") else ". "))
    vals = [per[c][0] for c in cams]
    if len(vals) > 1 and max(vals) - min(vals) > 0.1:
        text += (f"The difference between the cameras ({' vs '.join(f'{c} {per[c][0]:+.2f}' for c in cams)} m) "
                 f"cannot come from the water level or the setup as far as both cameras were measured over the same "
                 f"hours of {pc['window']}: it points at the pointing or the lens model; see the horizon offsets per "
                 f"camera in this README. ")
    text += (f"Until this product is built from the real photos of {pc['window']} (on the station computer only) "
             f"these are the only real numbers for this date: a build from other photos (e.g. the synthetic test "
             f"fixture, rendered with the setup planted) says nothing about the real beach. Read the per-camera lines "
             f"of the headline and the C = 0 sensitivity in compare/*_comparison.txt.")
    return text


def setup_caveats(plan, args, state, forcing, fits, era, prov=None):
    """Caveats about the setup coefficient actually used, its origin and its reference frame."""
    cav = []
    wv = (forcing or {}).get("waves", {})
    for f in fits:
        C = f["coef"]
        cams = f" (cameras {', '.join(f['cameras'])})" if len(fits) > 1 else ""
        if f.get("kind") == "none":
            cav.append(f"No wave-setup correction{cams}: each waterline is given the still-water level only, so "
                       f"lines marked by the swash read LOW, by roughly {SHARE_KEPT_RANGE[0]}-{SHARE_KEPT_RANGE[1]} x "
                       f"the setup (the setup is ~0.2-0.4 m in moderate waves; a line placed at still water also "
                       f"moves seaward onto lower beach, which takes back the rest).")
            continue
        if f.get("kind") == "unknown":
            cav.append(f"The setup coefficient C = {C}{cams} was given on the command line WITHOUT saying where it "
                       f"was fitted (--setup-fitted-to): no comparison of this product can be called INDEPENDENT.")
        elif f.get("kind") == "declared":
            cav.append(f"The setup coefficient C = {C}{cams}: {f.get('how')}. Its wave currency (which Hs and "
                       f"Tp it was fitted with) and still-water reference are as declared, not checked here.")
        else:
            cav.append(f"The setup coefficient C = {C}{cams} was fitted on this beach in Sep-Oct 2026 "
                       f"(waterline_timex_cron.sh) with {f.get('wave_currency')}; setup scales with the "
                       f"beach-face slope, which may have differed on this date.")
            if era in ("adcp", "chatham") and wv.get("tp_currency_setup_bias_m") is not None:
                cav.append(f"Wave-period currency: C = {C} was fitted with NDBC 44008 peak periods; this date "
                           f"uses the ADCP's (or WIS converted to it): expected setup bias "
                           f"{wv['tp_currency_setup_bias_m']:+.3f} m.")
            if era in ("adcp", "chatham"):
                cav.append(still_water_caveat(C, wv, era))
        same = [s for s in plan.get("surveys") or []
                if survey_matches_fit(f, s, survey_path_for(s, plan, args), args)[0]]
        if same:
            cav.append(f"{', '.join(s['name'] for s in same)}: CIRCULAR for as long as C = {C} is used, because C "
                       f"was fitted to those very shots (see 'C fitted' under WHAT WAS USED). That comparison cannot test "
                       f"the setup or the overall level of the lines; it still shows their spread over the elevations "
                       f"the survey covers (WHERE, under the comparison). A build with C = 0 (or a C fitted to other "
                       f"data, --setup-fitted-to) compares INDEPENDENTLY. An independent check of the whole "
                       f"intertidal needs another survey that fitted nothing, e.g. a calm, low-tide RTK across it.")
        fit = SETUP_FITS.get(C) if f.get("kind") == "known" else None
        rc = (fit or {}).get("rtk_check")
        if rc:
            su = (wv.get("setup_in_window") or {})
            med = su.get("daytime_median_m") if era != "live" else None
            med = med if med is not None else su.get("median_m")
            keep = rc.get("keep_paired") or 0.7
            own = build_setup_effect(prov, [Path(s_.get("path") or "").stem for s_ in plan.get("surveys") or []])
            inplace = rc.get("with_setup_in_place_m")
            reproj = rc["with_setup_m"] - (inplace or 0.0)
            prior = bool(prior_checks(plan))
            # what the re-projection takes back on THIS date: from this build's own paired C = 0 rise per
            # camera when there is one; on the live beach and days of the RTK, that RTK's own share; else
            # a range over the shares measured so far (they differ by date, camera and geometry)
            if not med:
                reproj_here = ""
            elif own and era != "live":
                reproj_here = (" Here the re-projection alone leaves the lines ~"
                               + ", ".join(f"{c} {med - e:.2f} m (on {what})" for c, (e, what) in sorted(own.items()))
                               + f" low (the setup median {med:.2f} m less this build's own paired rise per camera, "
                                 f"compare/*_comparison.txt), plus whatever C mismatch this date has")
            elif era == "live":
                reproj_here = (f" On this beach and these days the re-projection alone leaves this DEM "
                               f"~{(1 - keep) * med:.2f} m low (~{1 - keep:.1f} x the setup median of {med:.2f} m), "
                               f"plus the C mismatch")
            else:
                reproj_here = (f" The share kept depends on the date, camera and beach ({SHARE_KEPT_RANGE[0]}-"
                               f"{SHARE_KEPT_RANGE[1]} in the builds so far): the re-projection alone may leave this "
                               f"DEM ~{(1 - SHARE_KEPT_RANGE[1]) * med:.2f}-{(1 - SHARE_KEPT_RANGE[0]) * med:.2f} m "
                               f"low ({1 - SHARE_KEPT_RANGE[1]:.1f}-{1 - SHARE_KEPT_RANGE[0]:.1f} x the setup median of "
                               f"{med:.2f} m), plus whatever C mismatch this date has")
            cav.append(f"C = {C} under-corrects. On the 2026-09-29 RTK transects ({rc['frames']} frames of "
                       f"{rc['window']}, the same frames both ways) the lines read {rc['without_setup_m']:.2f} m low "
                       f"without setup and still {rc['with_setup_m']:.2f} m low with it (median setup applied "
                       f"{rc['setup_applied_m']:.2f} m). Two causes: (1) a line given the setup is re-projected "
                       f"landward onto higher beach, which takes back ~{reproj:.2f} m there (each line keeps "
                       f"~{keep:.2f} x its setup; the rest depends on the beach slope and the cross-shore distance); "
                       + (f"(2) C = {C} is below those frames' own per-frame C ({rc['c_per_frame_in_place']:.3f} with "
                          f"no re-projection), ~{inplace:.2f} m (C was fitted on 14 of them, each line taken where it "
                          f"lay WITHOUT setup)." if inplace is not None else "(2) the rest is a C mismatch on those frames.")
                       + reproj_here
                       + (" -- already inside the EARLIER CHECK expectation above: do not add it to that." if prior
                          else ".")
                       + (" That residual was measured in the GNSS-R frame (lines at the GNSS-R level); against the "
                          "lidar, this date's still-water-reference effect (above) adds to it"
                          + (" (the EARLIER CHECK expectation holds both)." if prior else ".")
                          if era in ("adcp", "chatham") else ""))
    return cav


# GNSS-R sits this much below the Chatham gauge's mean level after its +0.349 m datum fix (OPUS,
# +/-0.061 m): a comparison of MEAN levels, so it holds the GNSS-R's mean share of the setup.
GNSSR_BELOW_CHATHAM_M = 0.02


def still_water_caveat(C, wv, era):
    """
    C's still-water reference (GNSS-R, surf-zone footprint) against this date's (ADCP or Chatham,
    no setup), in words and, when measured, in metres. GNSS-R(t) = still water + s x setup(t) + its
    datum error; that datum error was set from mean levels (GNSS-R ~0.02 m below Chatham), so it holds
    -s x the mean setup of that comparison; the ADCP datum is tied to Chatham's mean. A line given
    C's setup then reads about +0.02 - s x (setup(t) - mean setup) against the frame C was fitted in:
    only the WAVE-DEPENDENT part of the share shows; its mean part is in the datum chain.
    """
    sw = wv.get("still_water_reference") or {}
    q = sw.get("share") if sw.get("quantified") else None
    su = wv.get("setup_in_window") or {}
    med = su.get("daytime_median_m") if su.get("daytime_median_m") is not None else su.get("median_m")
    text = ("Still-water reference of C: C was fitted with lines placed at the GNSS-R water level, and the GNSS-R "
            "footprint is the surf zone (gnssr_qc.py, gnssir_reflection_audit.py), where breaking waves raise the "
            "mean level: GNSS-R contains part of the setup (a share s), so C carries only the rest. This date's "
            "still water is the " + ("ADCP at 21 m depth" if era == "adcp" else "Chatham harbour gauge transferred")
            + ", outside the surf zone, with no setup. The MEAN part of that share is already in the datum chain: "
            f"GNSS-R sits ~{GNSSR_BELOW_CHATHAM_M:.2f} m below Chatham after its +0.349 m datum fix (OPUS +/-0.061 "
            "m), a comparison of mean levels that holds the GNSS-R's mean setup share, and the ADCP datum is tied to "
            "Chatham's mean. So, against the frame C was fitted in, a frame of this date is expected to read about "
            f"+{GNSSR_BELOW_CHATHAM_M:.2f} m - s x (its setup - the mean setup of the GNSS-R/Chatham comparison "
            "period): ~0.02 m HIGH in average waves, LOW only by the share of the setup ABOVE that average.")
    if q is not None:
        typ = sw.get("typical_m") or 0.0
        text += (f" Measured on the 2026 record ({sw.get('how')}): s = {q:.2f} (GNSS-R - Chatham transfer = "
                 f"{sw.get('k', 0):+.4f} x sqrt(Hs*L0) + const), a wave-dependent effect; at the record's median "
                 f"sqrt(Hs*L0) ({sw.get('typical_x', 0):.1f} m) the GNSS-R holds ~{typ:.2f} m of setup, the mean "
                 f"part already in the datum chain.")
        if med is not None:
            exp = GNSSR_BELOW_CHATHAM_M - (q * med - typ)
            text += (f" For this window's daytime median setup ({med:.2f} m) that gives "
                     f"{GNSSR_BELOW_CHATHAM_M:+.2f} - ({q:.2f} x {med:.2f} - {typ:.2f}) = {exp:+.2f} m.")
    else:
        text += (" (s is not measured here: " + (sw.get("note") or "the 2026 GNSS-R record is not on this computer")
                 + "; on the station historical_forcing.py measures it.)")
    text += (" To take the question out, C can be fitted in this date's own frame: rebuild with --setup-coef 0, run "
             "dem_from_contours.py --fit-setup on waterlines/contour_points_ground.csv (repeat crossings of the same "
             "cells at different wave heights: no survey), then rebuild with --setup-coef <that C> --setup-fitted-to "
             "'none:repeat crossings, <window>'; the label stays INDEPENDENT.")
    return text


def collect_caveats(plan, args, state, forcing, prov):
    cav = []
    era = plan["era"]
    out = plan["out"]
    wl = (forcing or {}).get("water_level", {})
    wv = (forcing or {}).get("waves", {})
    fits = (prov.get("setup") or {}).get("fits") or setup_in_use(plan, args)[0]
    if era in ("adcp", "chatham"):
        d = wl.get("adcp_datum") or {}
        cav.append("ADCP datum: the ADCP's NAVD88 offset was derived from the Chatham gauge's mean level over "
                   "the deployment (adcp_to_navd88.py), i.e. Marconi's mean level is assumed equal to "
                   f"Chatham's{' (' + str(d.get('source')) + ')' if d.get('source') else ''}. Every elevation "
                   "shares this datum error.")
        srcs = wl.get("sources_used") or {}
        if any("chatham" in str(k) for k in srcs):
            cav.append(f"Some water levels are the Chatham gauge transferred to Marconi ({wl.get('method')}; "
                       f"CV RMS {wl.get('cv_rms_m')} m, extrapolation RMS {wl.get('extrapolation_rms_m')} m): "
                       f"sources used {srcs}.")
        for c in wv.get("caveats") or []:
            cav.append("Waves: " + c)
    else:
        had = [c for c in plan["cams"] if (((read_stamp(out, f"detect_{c}") or {}).get("setup") or {}).get("had_setup"))]
        C = fits[0]["coef"] if fits else args.setup_coef
        if had:
            diffs = [((read_stamp(out, f"detect_{c}") or {}).get("setup") or {}).get("max_diff_existing") or 0.0
                     for c in had]
            cav.append(f"Live era: the station's rows already carried a setup (the cron applies it since 6 Oct 2026); "
                       + (f"this product recomputed it with C = {C} from each row's offshore_hs_m / offshore_tp_s "
                          f"(extract_elevation_contours.py's formula); the largest difference from the cron's value was "
                          f"{max(diffs):.4f} m." if C else
                          "this product set it to zero (C = 0): every line keeps its still-water (GNSS-R) level."))
        elif C:
            cav.append("Live era: the station's own detections were made before the setup correction was switched "
                       f"on (6 Oct 2026); this product adds the setup (C = {C}) to those rows with the "
                       "formula of extract_elevation_contours.py --setup-coef, from each row's own offshore_hs_m / "
                       "offshore_tp_s.")
        else:
            cav.append("Live era: the station's own detections, with NO wave setup added (C = 0): each line keeps "
                       "its still-water (GNSS-R) level.")
        cav.append("GNSS-R datum: the live water level is the gnssrefl spline put on NAVD88 through ONE survey of the "
                   "antenna, NGS OPUS with GEOID18 (+0.349 m from the earlier CSRS-PPP CGVD2013 height; ellipsoid "
                   "heights agree to 2 mm; +/-0.061 m, mostly the geoid model). Every elevation of this DEM shares "
                   "that offset. It largely cancels against RTK shots on NAVD88 (GEOID18), the same geoid model, but "
                   "not in the absolute NAVD88 heights of the GeoTIFFs.")
    # an earlier real measurement of these photos first: it can contradict the expectations below
    cav[:0] = [prior_check_caveat(plan, pc, fits, forcing, prov) for pc in prior_checks(plan)]
    cav += setup_caveats(plan, args, state, forcing, fits, era, prov)
    # frames that got no setup: left out of the waterlines (and so the DEM), counted
    m = read_stamp(out, "merge") or {}
    ns = m.get("no_setup_frames") or {}
    unc = wv.get("uncovered_spans") or []
    if ns:
        cav.append("Frames with NO wave record (no Hs/Tp within the gap limit, so no setup) were LEFT OUT of the "
                   "waterlines and the DEM: " + ", ".join(f"{c} {n}" for c, n in sorted(ns.items()))
                   + " frame(s) (waterlines/no_setup_frames.csv); kept, they would read ~one setup (0.2-0.4 m) low.")
    if unc:
        cav.append("Hours with no wave record in the forcing: " + "; ".join(unc[:6])
                   + (f" (+{len(unc) - 6} more)" if len(unc) > 6 else "") + ".")
    for cam, days in (prov.get("frames_per_day") or {}).items():
        window = days_between(plan["cams"][cam]["first"], plan["cams"][cam]["last"]) if cam in plan["cams"] else []
        empty = [d for d in window if not days.get(d)]
        if window and empty:
            cav.append(f"{cam}: waterlines on {len(window) - len(empty)} of the window's {len(window)} days "
                       f"({', '.join(d for d in window if days.get(d)) or 'none'}); none on {', '.join(empty)}.")
    # Pointing results of this run, else of the run that built them (a --steps compare rerun
    # must not lose these caveats).
    pointing = {}
    for cam in plan["cams"]:
        res = (state.get("pointing") or {}).get(cam) or load_json(out / "pointing" / f"pointing_{cam}.json")
        if res is not None:
            pointing[cam] = res
        if res is not None and not res.get("rows"):
            cav.append(f"{cam}: pointing NOT checked: no photo of the window was found in the photo roots.")
        note = (state.get("pointing_notes") or {}).get(cam)
        if note:
            cav.append(f"{cam}: {note}.")
        if (out / "maps").exists() and \
                not (out / "maps" / f"{plan['date']}_waterlines_{cam}.png").exists():
            cav.append(f"{cam}: no waterline map on the photos: {map_missing_reason(plan, state, cam)} (the "
                       f"plan-view map does not need the photos).")
    for cam, res in pointing.items():
        info = res.get("info", {})
        skip = (read_stamp(out, f"detect_{cam}") or {}).get("skip_days")
        if skip:
            cav.append(f"{cam}: days left out for a different pointing: {', '.join(skip)}.")
        elif res.get("different_days"):
            cav.append(f"{cam}: days with a different pointing KEPT in the detection: "
                       f"{', '.join(res['different_days'])}"
                       + (" (--keep-moved-days)." if getattr(args, "keep_moved_days", False) else "."))
        eo_ = info.get("calibration") or plan["cams"].get(cam, {}).get("eo_file")
        off = pointing_offset(info, eo_)
        if off:
            dt, dr = off[0], off[1]
            sens = tilt_sensitivity(plan, cam, dt, dr)
            if sens:
                off = pointing_offset(info, eo_, rate=False) or off
            text = f"{cam}: {off[2]}"
            if sens:
                text += (f". From this date's own DEM and calibration, a pointing error of {dt:+.2f}/{dr:+.2f} deg "
                         f"(tilt/roll) would shift the DEM by: "
                         + ", ".join(f"{r['dem_error_m']:+.2f} m at {r['range_m']} m (lines moved "
                                     f"{r['shift_m']:.1f} m)" for r in sens)
                         + " -- if the offset is the pointing and not the lens model")
                test = pointing_vs_survey(plan, cam, (dt, dr), sens, prov, fits)
                if test:
                    text += ". " + test
            cav.append(text + ".")
        if info.get("warning"):
            cav.append(f"{cam}: {info['warning']}.")
        rows_ = res.get("rows") or []
        no_pan = [r["date"] for r in rows_ if r.get("verdict") in ("horizon only", "unchecked")]
        if rows_ and len(no_pan) == len(rows_):
            cav.append(f"{cam}: the pan (azimuth) could not be measured on any day (too few matched features "
                       f"against the reference frame): the days were checked by the sea horizon only (tilt and "
                       f"roll), so a turn of the camera would not have been seen.")
        elif no_pan:
            cav.append(f"{cam}: pan not measured on {', '.join(no_pan)} (horizon only).")
        if info.get("reference_kind", "").startswith("the day nearest"):
            cav.append(f"{cam}: no photo of the calibration day in the window, so the feature check only shows "
                       f"that the days agree with each other ({info.get('reference_day')}), not with the calibration.")
    if plan.get("disabled_cams"):
        for cam, why in plan["disabled_cams"].items():
            cav.append(f"{cam} not built (disabled row): {why}")
    for cam, why in failed_cameras(plan, state).items():
        cav.append(f"{cam} has NO waterlines in this product (PARTIAL build): {why}.")
    if state.get("filter_failed"):
        cav.append("The consistency filter failed: the DEM and maps use the UNFILTERED waterlines.")
    elif read_stamp(out, "merge") and not waterline_file(plan, quiet_=True).name.endswith("_filtered.csv"):
        cav.append("No current consistency-filtered waterlines (the filter step failed or has not run on the "
                   "current merged file): the DEM, maps and comparison use the UNFILTERED waterlines.")
    epsg = ((dem_settings_built(out) or {}).get("geotiff_epsg")      # what the GeoTIFFs on disk carry
            or getattr(args, "geotiff_epsg", GEOTIFF_EPSG))
    cav.append(f"Coordinates: the grids are in the frame of the calibration's GCPs, ASSUMED to be that of the "
               f"surveys, NAD83(2011) / UTM zone 19N (EPSG:6348): the GCP files state no CRS. Were the GCPs in "
               f"WGS 84 / ITRF instead, the grid would sit ~1-1.5 m off the surveys, up to ~0.1-0.15 m of "
               f"elevation where that shift runs across a 1:10 foreshore. The GeoTIFFs are tagged EPSG:{epsg}"
               + (" (WGS 84 / UTM 19N) as asked, a NOMINAL tag: read them as NAD83(2011) (WGS 84 differs by "
                  "~1-1.5 m here; nothing was transformed; --geotiff-epsg 6348 writes the true code)"
                  if epsg == 32619 else "")
               + ". The surveys are sampled at the DEM's coordinates with no shift.")
    first_used, last_used, _ = frames_used_dates(plan)
    for s in plan["surveys"]:
        if s["survey_type"] == "dsm":
            cav.append(f"{s['name']}: lidar sees only the beach that was DRY at flight time; DEM cells below "
                       f"the water line then have no comparison.")
        if s.get("survey_date"):
            sd = date_cls.fromisoformat(s["survey_date"])
            a, b = date_cls.fromisoformat(first_used), date_cls.fromisoformat(last_used)
            gap = 0 if a <= sd <= b else min(abs((sd - a).days), abs((sd - b).days))
            dp = (sd - date_cls.fromisoformat(plan["date"])).days
            if gap or dp:
                cav.append(f"{s['name']}: surveyed {s['survey_date']}"
                           + (f", {abs(dp)} day(s) {'after' if dp > 0 else 'before'} the product date {plan['date']}"
                              if dp else "")
                           + (f"; {gap} day(s) outside the days of the frames used ({first_used} .. {last_used})"
                              if gap else f"; inside the days of the frames used ({first_used} .. {last_used})")
                           + ": the beach may have moved in between.")
        js = load_json(out / "compare" / f"{plan['date']}_{s['name']}_comparison.json") or {}
        if js.get("vertical_datum_note"):
            cav.append(f"{s['name']}: {js['vertical_datum_note']}.")
    cav += storm_caveats(plan, args, prov, first_used, last_used)
    cav += wave_event_caveats(plan, prov, first_used, last_used)
    return cav


# A setup coefficient fitted to NO survey: repeat crossings of the same DEM cells agree best at C =
# 0.03-0.04 (dem_from_contours.py --fit-setup on lines without setup; waterline_timex_cron.sh). A
# constant pointing error moves every frame's line alike, so it does not change how well repeat
# crossings agree: this C is blind to it, unlike the C fitted to the 2026 RTK (which absorbs any level
# error at the transects, of the setup or of the pointing). Stockdon et al. (2006): setup = 0.35 x
# foreshore slope x sqrt(H0 L0), so C ~0.035-0.045 on a 1:10-1:8 beach face.
SETUP_C_INTERNAL = (0.03, 0.04)
# The share of its setup a line keeps once re-projected landward at still water + setup, as measured
# in the paired C = 0 sensitivities of the builds (Jan 2025 ~0.6, the 2026 RTK 0.71, Mar 2025 ~0.8):
# used only where a build has no paired measurement of its own (a C = 0 build).
SHARE_KEPT_RANGE = (0.6, 0.8)


def compared_frames_scale(plan, name, cam):
    """Median sqrt(Hs*L0) (m) over the frames of camera `cam` compared on the survey's transects
    (compare/<name>_waterline_frames.csv), from the waterlines' own offshore_hs_m / offshore_tp_s.
    -> (scale, frames) or (None, 0)."""
    fr = plan["out"] / "compare" / f"{name}_waterline_frames.csv"
    try:
        with open(fr, newline="") as f:
            want = {r["frame"] for r in csv.DictReader(f) if r.get("camera") == cam}
    except (OSError, KeyError):
        return None, 0
    if not want:
        return None, 0
    want |= {w + ".jpg" for w in want}
    x = {}
    try:
        with open(waterline_file(plan, quiet_=True), newline="") as f:
            for r in csv.DictReader(f):
                k = r.get("source_file", "")
                if k in want and k not in x:
                    try:
                        hs, tp = float(r.get("offshore_hs_m") or "nan"), float(r.get("offshore_tp_s") or "nan")
                    except ValueError:
                        continue
                    if hs == hs and tp == tp:
                        x[k] = float(np.sqrt(hs * G * tp ** 2 / (2 * np.pi)))
    except OSError:
        return None, 0
    return (float(np.median(list(x.values()))), len(x)) if x else (None, 0)


def survey_range_from(plan, name, cam, zlo, zhi):
    """Median distance (m) from camera `cam` of the survey's transect points in its view with
    elevations in [zlo, zhi] (where the lines were compared), from compare/<name>_points.csv."""
    try:
        with open(plan["out"] / "compare" / f"{name}_points.csv", newline="") as f:
            rows = list(csv.DictReader(f))
    except OSError:
        return None
    d = []
    for tr in (True, False):
        for r in rows:
            try:
                z, dist = float(r.get("survey_z_navd88") or "nan"), float(r.get("distance_from_camera_m") or "nan")
            except ValueError:
                continue
            if r.get("camera") != cam or not (zlo <= z <= zhi) or dist != dist:
                continue
            if tr and "transect" not in (r.get("description", "") + r.get("name", "")).lower():
                continue
            d.append(dist)
        if d:
            break
    return float(np.median(d)) if d else None


def pointing_vs_survey(plan, cam, offset, sens, prov, fits):
    """
    Can a survey's waterline check tell a pointing error from a lens-model one? Only on evidence that
    fitted nothing to that survey: the lines WITHOUT setup (C = 0) against the survey, minus the setup a
    pointing-blind C (SETUP_C_INTERNAL, repeat crossings) explains, leaves what pointing, swash and
    detector bias share; compared with what the full offset would do at the survey's range from the
    camera (tilt_sensitivity), it bounds the part of the offset that cannot be ruled out. A check with
    a C fitted to this same survey cannot say (that C absorbs any level error there). -> text or None.
    """
    C = float(fits[0]["coef"]) if fits and fits[0].get("coef") is not None else None
    for c in prov.get("comparisons") or []:
        h = c.get("headline") or {}
        wl = h.get("waterlines") or {}
        # only a check at a known range from the camera (frames on RTK transects): a DSM comparison
        # pools every range, where the offset's effect varies from nothing to metres
        if not c.get("current") or not wl.get("n") or not wl.get("by_camera") or wl.get("unit") != "frames":
            continue
        name = f"{plan['date']}_{c['name']}"
        no_setup = str(wl.get("setup", "")).startswith("NO setup") or not C
        if no_setup:
            v = (wl.get("by_camera") or {}).get(cam) or {}
            obs0, n0 = v.get("median"), v.get("frames") or v.get("n")
            scale, nsc = compared_frames_scale(plan, name, cam)
            if obs0 is None or not scale:
                continue
            per_c = (SHARE_KEPT_RANGE[0] * scale, SHARE_KEPT_RANGE[1] * scale)
            how = (f"its frames' median sqrt(Hs L0) of {scale:.1f} m x {SETUP_C_INTERNAL[0]}-{SETUP_C_INTERNAL[1]} "
                   f"x the share of its setup a line keeps once re-projected, {SHARE_KEPT_RANGE[0]}-"
                   f"{SHARE_KEPT_RANGE[1]} as measured in the paired builds")
        else:
            v = (((wl.get("without_setup") or {}).get("paired") or {}).get("by_camera") or {}).get(cam) or {}
            obs0, n0, eff = v.get("without_setup_median"), v.get("frames") or v.get("n"), v.get("setup_effect_median")
            if obs0 is None or eff is None:
                continue
            per_c = (eff / C, eff / C)
            how = (f"this build's own paired effect, {eff:+.3f} m at C = {C:g} with the re-projection, scaled to "
                   f"C = {SETUP_C_INTERNAL[0]}-{SETUP_C_INTERNAL[1]}")
        lo = obs0 + per_c[0] * SETUP_C_INTERNAL[0]
        hi = obs0 + per_c[1] * SETUP_C_INTERNAL[1]
        zc = wl.get("line_elevation_compared") or {}
        rng = survey_range_from(plan, name, cam, (zc.get("p5") or 0.0) - 0.2, (zc.get("max") or 9.0) + 0.2)
        if rng is None:
            continue
        band = None
        for r in sens:
            a, b = (float(x) for x in r["range_m"].split("-"))
            if a <= rng < b:
                band = r
        if band is None:
            band = min(sens, key=lambda r: min(abs(rng - float(x)) for x in r["range_m"].split("-")))
        pred = band["dem_error_m"]
        circ = (" The check of these lines WITH the setup cannot say: C = %g was fitted to this same survey, so any "
                "level error on its transects, of the setup or of the pointing, is absorbed into C (that comparison "
                "is CIRCULAR for the level)." % C) if (C and any(f.get("kind") == "known" for f in fits)
                                                       and c.get("label") == "CIRCULAR") else ""
        few_ = f" ({n0} frames: too few for an estimate)" if (n0 or 0) < FEW_N else f" ({n0} frames)"
        head = (f"Is it the pointing?{circ} Judged on what fitted nothing to {c['name']}: repeat crossings of the same "
                f"DEM cells agree best at C = {SETUP_C_INTERNAL[0]}-{SETUP_C_INTERNAL[1]} (dem_from_contours.py "
                f"--fit-setup on lines without setup), blind to a constant pointing error, which moves every frame's "
                f"line alike, and Stockdon's 0.35 x a ~1:10 beach face expects about that. Without setup the {cam} "
                f"lines read {obs0:+.2f} m against the survey on its transects{few_}; that setup raises them by "
                f"{lo - obs0:.2f}-{hi - obs0:.2f} m ({how}), leaving {lo:+.2f} to {hi:+.2f} m")
        where = f" at the transects (~{rng:.0f} m from {cam}, the {band['range_m']} m band)"
        if pred >= -0.02 and pred <= 0.02:
            return head + where + f", where the offset would barely move the DEM ({pred:+.2f} m): this survey cannot test it"
        f_lo, f_hi = sorted((lo / pred, hi / pred))
        f_lo, f_hi = max(f_lo, 0.0), max(f_hi, 0.0)
        if f_hi <= 0.0:
            verdict = (f", of the opposite sign to the {pred:+.2f} m the full offset would cause there: the "
                       f"transects show no sign of it being pointing (on {n0} frames)")
        elif f_hi >= 1.0:
            verdict = (f", as much as the {pred:+.2f} m the full offset would cause there: the transects cannot rule "
                       f"out that the whole offset is pointing; swash and detector bias at the waterline are "
                       f"confounded with it")
        else:
            verdict = (f", which is {f_lo * 100:.0f}-{f_hi * 100:.0f}% of the {pred:+.2f} m the full offset would "
                       f"cause there: the full offset is unlikely to be pointing, but a part of it, up to ~"
                       f"{f_hi * offset[0]:+.2f}/{f_hi * offset[1]:+.2f} deg tilt/roll, cannot be ruled out; swash "
                       f"and detector bias at the waterline are confounded with it")
        return head + where + verdict + ". Read the shifts above as what the offset WOULD do if it were all pointing"
    return None


def storm_caveats(plan, args, prov, first_used, last_used):
    """A caveat per known storm (STORMS) that reaches into the window or ended within 3 days of it:
    which frame days are the storm's tail before the survey, what changed and where, and whether a
    point survey lies on the part that changed most."""
    cav = []
    days = sorted({d for v in (prov.get("frames_per_day") or {}).values() for d in v}) or \
        days_between(first_used, last_used)
    for sev in STORMS:
        a, b = date_cls.fromisoformat(sev["first"]), date_cls.fromisoformat(sev["last"])
        w0, w1 = date_cls.fromisoformat(plan["first"]), date_cls.fromisoformat(plan["last"])
        if b < w0 - timedelta(days=3) or a > w1:
            continue
        where = (f"reaches into this window ({plan['first']} .. {plan['last']})" if b >= w0 else
                 f"ended {(w0 - b).days} day(s) before this window")
        sdays = sorted({x["survey_date"] for x in plan.get("surveys") or [] if x.get("survey_date")})
        sd = sdays[0] if sdays else plan["date"]
        tail = [d for d in days if sev["first"] <= d < sd]
        text = (f"STORM: {sev['name']} ({sev['first']} .. {sev['last']}) {where}"
                + (f"; the frames of {', '.join(tail)} (the storm's tail, before the {sd} survey) are gridded with "
                   f"the rest" if tail else "")
                + f". On this beach {sev['change']}.")
        lo, hi = sev.get("loss_northing") or (None, None)
        for s in plan.get("surveys") or []:
            if s["survey_type"] != "points" or lo is None:
                continue
            try:
                from survey_compare import read_points, is_transect
                q = read_points(survey_path_for(s, plan, args))
                tr = np.array([is_transect(x) for x in q["desc"]], bool)
                n = q["N"][tr] if tr.sum() >= 3 else q["N"]
                if len(n) and n.min() <= hi and n.max() >= lo:
                    text += (f" The {s['name']} {'transects' if tr.sum() >= 3 else 'points'} span N "
                             f"{n.min():.0f}-{n.max():.0f}: they overlap that loss.")
            except (SystemExit, Exception):
                pass
        text += " " + PER_DAY_ROWS_NOTE
        cav.append(text)
    return cav


# What a step or trend in the per-day rows of waterlines - survey can and cannot say. On the sandbox's
# synthetic fixture (every day rendered from ONE lidar: a static beach) this method's per-day medians
# drift by ~0.09 m over Jan 18-23 (+0.01 to -0.08 m) and ~0.17 m over Mar 3-9 (+0.06 to -0.10 m), with a
# ~0.07 m step across the 20 Jan wave event (survey_build/fixture/README.md, 'By day'); 5 and 8 Mar, both
# at a median survey elevation of -0.07 m, read +0.023 and -0.103 m (Oct 2026 fixture build): the error
# depends on the elevation, the waves and setup and the tide phase of the hours sampled, not on one alone.
PER_DAY_ROWS_NOTE = (
    "The per-day rows of waterlines - survey (compare/*_comparison.txt) show whether the lines moved over the "
    "window, but a step or trend there may be change of the beach OR a method error that depends on the "
    "conditions (the elevation, the waves and setup, the tide phase of the hours sampled, which moves from day to "
    "day): on a static synthetic beach this method's per-day medians drift by up to ~0.17 m over a week (Jan "
    "-0.09 m, Mar -0.17 m, with a ~0.07 m step across a wave event), and days at the same median elevation still "
    "differ by up to ~0.13 m. A difference between days below ~0.15-0.2 m is not by itself evidence of change; "
    "compare days at a similar survey elevation and setup (the rows' 'at z' and 'setup' columns, or the elevation "
    "bands) before calling it change.")


def frame_epochs(plan):
    """{camera: {source_file: capture epoch}} of the waterlines in use (for placing frames against an
    event by time, not by calendar day)."""
    res = {}
    try:
        with open(waterline_file(plan, quiet_=True), newline="") as f:
            for r in csv.DictReader(f):
                k, cam = r.get("source_file", ""), r.get("camera", "?")
                d = res.setdefault(cam, {})
                if k in d:
                    continue
                try:
                    d[k] = float(r.get("capture_epoch") or "nan")
                except ValueError:
                    d[k] = float("nan")
                if d[k] != d[k]:
                    e = epoch_of(k)
                    d[k] = float(e) if e is not None else float("nan")
    except OSError:
        return {}
    return res


# Offshore Hs (the ADCP currency of the forcing) from which a span of the window counts as a wave
# event that may have reshaped the beach face between the frames gridded together (the DEM leaves
# out frames above its own max Hs, 1.5 m, but not the days on either side of the event).
EVENT_HS_M = 2.0


def wave_events(fdir, first, last, threshold=EVENT_HS_M, max_gap_h=6.0):
    """Spans of the forcing's waves (forcing/waves.csv) in the window with Hs >= threshold.
    -> list of dicts (start, end epochs, peak Hs and its epoch, source)."""
    ep, hs, src = [], [], []
    try:
        with open(Path(fdir) / "waves.csv", newline="") as f:
            for r in csv.DictReader(f):
                try:
                    e, h = float(r["epoch"]), float(r["wvht_m"])
                except (KeyError, TypeError, ValueError):
                    continue
                ep.append(e)
                hs.append(h)
                src.append(r.get("source") or "")
    except OSError:
        return []
    if not ep:
        return []
    t0 = datetime.fromisoformat(first).replace(tzinfo=timezone.utc).timestamp()
    t1 = (datetime.fromisoformat(last) + timedelta(days=1)).replace(tzinfo=timezone.utc).timestamp()
    o = np.argsort(ep)
    ep, hs, src = np.array(ep)[o], np.array(hs)[o], np.array(src, dtype=object)[o]
    m = (ep >= t0) & (ep < t1) & np.isfinite(hs) & (hs >= threshold)
    events = []
    for e, h, s_ in zip(ep[m], hs[m], src[m]):
        if events and e - events[-1]["end"] <= max_gap_h * 3600:
            ev = events[-1]
            ev["end"] = e
            if h > ev["peak"]:
                ev["peak"], ev["peak_at"] = h, e
            ev["sources"].add(s_)
        else:
            events.append({"start": e, "end": e, "peak": h, "peak_at": e, "sources": {s_}})
    return events


def rough_days(path, max_hs):
    """{camera: [days]} on which EVERY frame of the waterline file has offshore Hs > max_hs (the DEM
    leaves such a day out entirely)."""
    per = {}
    try:
        with open(path, newline="") as f:
            rd = csv.DictReader(f)
            seen = set()
            for r in rd:
                k = r.get("source_file")
                if k in seen:
                    continue
                seen.add(k)
                d = per.setdefault(r.get("camera", "?"), {}).setdefault((r.get("capture_time_utc") or "")[:10],
                                                                         [0, 0])
                d[0] += 1
                try:
                    d[1] += float(r.get("offshore_hs_m") or "nan") > max_hs
                except ValueError:
                    pass
    except OSError:
        return {}
    return {c: sorted(day for day, (n, r_) in v.items() if n and r_ == n) for c, v in per.items()
            if any(n and r_ == n for n, r_ in v.values())}


def wave_event_caveats(plan, prov, first_used, last_used):
    """A caveat per span of the window with offshore Hs >= EVENT_HS_M in the date's own forcing (the
    hand-kept STORMS table covers only the events measured on the station): the frame days before and
    after it, where the survey lies, the days the DEM left out entirely for their waves, and where
    to look for change. An event inside a STORMS entry is left to that caveat."""
    out = plan["out"]
    cav = []
    days = sorted({d for v in (prov.get("frames_per_day") or {}).values() for d in v}) or \
        days_between(first_used, last_used)
    sdays = sorted({x["survey_date"] for x in plan.get("surveys") or [] if x.get("survey_date")})
    sd = sdays[0] if sdays else plan["date"]
    max_hs = (dem_settings_built(out) or {}).get("max_hs_m")
    dropped = rough_days(waterline_file(plan, quiet_=True), float(max_hs)) if max_hs else {}
    events = wave_events(out / "forcing", plan["first"], plan["last"])
    fep = frame_epochs(plan) if events else {}
    for ev in events:
        a = datetime.fromtimestamp(ev["start"], tz=timezone.utc)
        b = datetime.fromtimestamp(ev["end"], tz=timezone.utc)
        da, db = a.strftime("%Y-%m-%d"), b.strftime("%Y-%m-%d")
        if any(not (db < sv["first"] or da > sv["last"]) for sv in STORMS):
            continue
        pk = datetime.fromtimestamp(ev["peak_at"], tz=timezone.utc)
        # frames placed by their capture TIME against the event (a wave record stands for the half hour
        # on either side of it), not by calendar day: a one-hour event is not 'during' a whole day
        t_a, t_b = ev["start"] - 1800.0, ev["end"] + 1800.0
        cls = {"before": {}, "during": {}, "after": {}}
        for cam, fr in (fep or {}).items():
            for k, e in fr.items():
                if e != e:
                    continue
                c_ = "before" if e < t_a else "after" if e > t_b else "during"
                dd = datetime.fromtimestamp(e, tz=timezone.utc).strftime("%Y-%m-%d")
                cls[c_][dd] = cls[c_].get(dd, 0) + 1
        if not any(cls.values()):                     # no capture times: by day, as before
            cls = {"before": {d: 0 for d in days if d < da}, "during": {d: 0 for d in days if da <= d <= db},
                   "after": {d: 0 for d in days if d > db}}

        def frames_text(k, what):
            v = cls[k]
            if not v:
                return ""
            n = sum(v.values())
            return (f"{n} frame(s) {what} it ({', '.join(sorted(v))})" if n else
                    f"the frames of {', '.join(sorted(v))} {what} it")
        # the survey's time of day is not known: on an event day it is 'the same day', not before/after
        rel = ("after it" if sd > db else "before it" if sd < da else
               "on the same day as the event (its time is not known)" if da == db else
               "within the event's days (its time is not known)")
        span = (f"{a:%Y-%m-%d %H:%M} .. {b:%Y-%m-%d %H:%M} UTC" if ev["end"] > ev["start"] else
                f"{a:%Y-%m-%d %H:%M} UTC")
        srcs = ", ".join(sorted(x for x in ev["sources"] if x)) or "the forcing"
        text = (f"HIGH WAVES in the window: offshore Hs >= {EVENT_HS_M:.1f} m over {span} (peak {ev['peak']:.2f} m "
                f"at {pk:%Y-%m-%d %H:%M} UTC; {srcs}, forcing/waves.csv). The beach may have changed across it, and "
                f"the DEM grids the frames of both sides together: "
                + "; ".join(x for x in (frames_text("before", "before"), frames_text("during", "during"),
                                        frames_text("after", "after")) if x)
                + f"; the survey ({sd}) is {rel}.")
        gone = [f"{c} {d}" for c, ds in sorted(dropped.items()) for d in ds if da <= d <= db or d in cls["during"]]
        if gone:
            text += (f" The DEM left out every frame of {', '.join(gone)} (offshore Hs > {float(max_hs):g} m, its "
                     f"wave filter).")
        text += " " + PER_DAY_ROWS_NOTE
        cav.append(text)
    return cav


def few(n):
    return n is None or n < FEW_N


# Per-camera medians further apart than this (m) are flagged: on a two-camera date a pooled median
# near zero can hide each camera being off by half of it, with opposite signs (the Jan 2025 earlier
# check: c1 and c2 0.3 m apart on the same water level).
CAMERA_SPLIT_M = 0.10


def camera_medians(bc):
    """{camera: (median, n)} of the real cameras in a by_camera dict (survey_compare.py's 'neither'
    group, outside every view, left out)."""
    return {k: (float(v["median"]), int(v.get("n") or 0)) for k, v in sorted((bc or {}).items())
            if k != "neither" and v and v.get("n") and v.get("median") is not None}


def _camera_groups(h):
    """(what, by_camera dict, unit) of a comparison's per-camera statistics: DEM - survey (DEM cells or
    points: n itself) and waterlines - survey (frames or distinct survey points)."""
    wl = h.get("waterlines") or {}
    return (("DEM - survey", h.get("by_camera") or {}, None),
            ("waterlines - survey", wl.get("by_camera") or {}, wl.get("unit") or "points"))


def _split_of(bc, unit, limit):
    """-> (medians {cam: (median, n)}, spread m, rests-on {cam: (count, what)}) or None for < 2 cameras."""
    m = camera_medians(bc)
    if len(m) < 2:
        return None
    vals = [v[0] for v in m.values()]
    rests = {c: (independent_count(bc[c], unit) if unit else (m[c][1], "values")) for c in m}
    return m, max(vals) - min(vals), rests


def _split_cause(what, bc, wl):
    """Why two cameras may differ, from what their compared values share. The DEM cells of two cameras
    are different parts of the beach; waterlines are the cameras' own only when both were compared at
    the same times (same water level, waves and setup)."""
    if what.startswith("DEM"):
        return ("the cameras' cells are different parts of the beach, made from different frames: pointing, "
                "lens model, or the beach and the frames themselves")
    cams = sorted(c for c in bc if c != "neither")
    tm = {c: (bc[c].get("times"), bc[c].get("shared_times")) for c in cams}
    if all(t is not None and s_ is not None for t, s_ in tm.values()):
        if all(t and s_ / float(t) >= 0.8 for t, s_ in tm.values()):
            return ("compared at the same times (" + ", ".join(f"{c} {s_} of {t}" for c, (t, s_) in tm.items())
                    + "): the water level and the setup are the same for both, so pointing or lens model")
        txt = ("compared at different times (" + ", ".join(f"{c} {t} time(s), {s_} shared" for c, (t, s_) in
                                                           tm.items()) + "): the water level, waves and setup differ "
               "too, so this need not be the cameras")
    else:
        txt = ("not known to be compared at the same times: the water level, waves and setup may differ too, so "
               "this need not be the cameras")
    su = {c: bc[c].get("setup_median") for c in cams if bc[c].get("setup_median") is not None}
    if su and any(abs(v) > 1e-6 for v in su.values()):
        txt += "; median setup of the frames compared " + ", ".join(f"{c} {v:.2f} m" for c, v in su.items())
    hs = {c: bc[c].get("hs_median") for c in cams if bc[c].get("hs_median") is not None}
    if hs:
        txt += "; their median offshore Hs " + ", ".join(f"{c} {v:.2f} m" for c, v in hs.items())
    pr = (((wl or {}).get("without_setup") or {}).get("paired") or {}).get("by_camera") or {}
    pp = {c: pr[c] for c in cams if c in pr and pr[c].get("with_setup_median") is not None
          and pr[c].get("without_setup_median") is not None}
    if len(pp) > 1:
        txt += ("; with / without the setup " + ", ".join(
            f"{c} {v['with_setup_median']:+.3f} / {v['without_setup_median']:+.3f}" for c, v in pp.items()))
    return txt


def camera_split(h, limit=CAMERA_SPLIT_M):
    """Text when the per-camera medians of a comparison (DEM - survey, or waterlines - survey) are more
    than `limit` apart AND every camera's median rests on at least FEW_N values (DEM cells or points;
    frames or distinct survey points for the waterlines), else ''. The cause is named only from what
    the values share (_split_cause)."""
    out = []
    for what, bc, unit in _camera_groups(h):
        r = _split_of(bc, unit, limit)
        if not r or r[1] <= limit or any(k < FEW_N for k, _ in r[2].values()):
            continue
        m, spread, rests = r
        out.append(f"{what} by camera " + ", ".join(f"{c} {v[0]:+.3f} ({rests[c][0]} {rests[c][1]})"
                                                    for c, v in m.items())
                   + f": {spread:.2f} m apart, more than {limit:.2f} m: the pooled median hides a camera "
                     f"difference ({_split_cause(what, bc, h.get('waterlines'))})")
    return "; ".join(out)


def camera_split_note(h, limit=CAMERA_SPLIT_M):
    """Text when per-camera medians are more than `limit` apart but a camera rests on fewer than FEW_N
    values: too few to compare the cameras (no flag), else ''."""
    out = []
    for what, bc, unit in _camera_groups(h):
        r = _split_of(bc, unit, limit)
        if not r or r[1] <= limit or not any(k < FEW_N for k, _ in r[2].values()):
            continue
        m, spread, rests = r
        out.append(f"{what}: too few to compare the cameras (" + ", ".join(
            f"{c} {v[0]:+.3f} on {rests[c][0]} {rests[c][1]}" for c, v in m.items())
                   + f"; {spread:.2f} m apart, but below {FEW_N} a median is not an estimate)")
    return "; ".join(out)


def independent_count(v, unit=None):
    """(count, what) a waterline statistic rests on (survey_compare.independent_n): frames on RTK
    transects, the DISTINCT survey points of a points survey, frames on a DSM; for DEM cells or
    points, n itself."""
    v = v or {}
    if v.get("n_independent") is not None:
        return int(v["n_independent"]), v.get("n_independent_unit") or "values"
    if unit == "frames":
        return int(v.get("n") or 0), "frames"
    if v.get("survey_points") is not None:
        return int(v["survey_points"]), "survey points"
    if unit == "points" and v.get("frames") is not None:
        return int(v["frames"]), "frames"
    return int(v.get("n") or 0), "values"


def per_camera_text(by_cam, unit=None):
    """'c1 +0.120 m (n 345), c2 ...' from a {camera: {n, median, ...}} dict ('' if none). For the
    waterlines (unit 'frames' or 'points') the frames and the distinct survey points are given too,
    and 'too few' is judged on what the statistic rests on (frames, survey points), not on n."""
    parts = []
    for k, v in sorted((by_cam or {}).items()):
        n = (v or {}).get("n") or 0
        if n and v.get("median") is not None:
            # survey_compare.py's 'neither': cells outside every camera's footprint (from the calibrations)
            name = "outside the camera views" if k == "neither" else k
            if unit is None:
                parts.append(f"{name} {v['median']:+.3f} m (n {n}{', too few' if n < FEW_N else ''})")
                continue
            k_, what = independent_count(v, unit)
            det = (f"n {n} frames" if unit == "frames" else
                   f"n {n}, {v.get('frames')} frames" + (f", {v['survey_points']} survey points"
                                                         if v.get("survey_points") is not None else ""))
            parts.append(f"{name} {v['median']:+.3f} m ({det}{', too few ' + what if k_ < FEW_N else ''})")
    return ", ".join(parts)


def stat_text(h, prefix="DEM - survey"):
    """'median ..., NMAD ..., RMSE ..., n ...' with '-' where n < 3 and a 'too few' marker below FEW_N
    (for the waterlines: below FEW_N frames or distinct survey points, whatever n is)."""
    n = h.get("n") or 0
    if not n:
        return f"{prefix}: nothing compared"
    nm = f"{h['nmad']:.3f} m" if n >= 3 and h.get("nmad") is not None and np.isfinite(h["nmad"]) else "-"
    k, what = (independent_count(h, h.get("unit")) if h.get("n_independent") is not None else (n, "n"))
    return (f"{prefix}: median {h['median']:+.3f} m, NMAD {nm}, RMSE {h['rmse']:.3f} m, n {n}"
            + (f" -- TOO FEW (n < {FEW_N}): not an estimate" if k < FEW_N and what == "n" else
               f" -- TOO FEW ({k} {what} < {FEW_N}): not an estimate" if k < FEW_N else ""))


def write_provenance(plan, args, state, final=True):
    out = plan["out"]
    out.mkdir(parents=True, exist_ok=True)
    forcing = load_json(out / "forcing" / "forcing.json")
    cams = {}
    for cam, c in plan["cams"].items():
        hist = calibration_history(c["eo_file"], args.calibration)
        pres = (state.get("pointing") or {}).get(cam) or load_json(out / "pointing" / f"pointing_{cam}.json")
        dst = read_stamp(out, f"detect_{cam}")
        cams[cam] = {
            "station": c["station"], "window": [c["first"], c["last"]], "utc_hours": c["hours_text"],
            "calibration": {"eo_file": c["eo_file"], "eo_path": str(c["eo_path"]),
                            "eo_sha256": sha256(c["eo_path"])[:16] if c["eo_path"] else None,
                            "io_file": Path(c["io_path"]).name, "history": hist},
            "search_envelope": ({"survey": str(c["envelope_path"] or c["envelope_survey"]),
                                 "how": "detect_original_view.py: the survey's -1.5..+2.5 m band projected "
                                        "into the photo with this calibration"}
                                if plan["era"] != "live" else
                                {"survey": None, "how": "the live detector's own envelope "
                                 "(waterline_detector_v5.py CameraProfile, from hand-clicked ground truth: "
                                 "derive_envelope.py), not a survey"}),
            "photos": photo_section(plan, cam),
            "pointing_check": ({"reference": pres.get("info"),
                                "left_out_days": (dst or {}).get("skip_days"),
                                "different_days": pres.get("different_days"),
                                "keep_moved_days": bool(args.keep_moved_days),
                                "station_monitor_log": pres.get("station_log"),
                                "per_day": [{k: r[k] for k in ("date", "verdict", "why", "horizon_offset_px",
                                                                "d_azimuth_deg", "d_tilt_deg", "d_roll_deg")}
                                            for r in pres.get("rows", [])]} if pres else None),
            "detection": (dst or {}).get("setup") if plan["era"] == "live" else
            ({"photos": (dst or {}).get("photos"), "skip_days": (dst or {}).get("skip_days")} if dst else None),
            "built": bool(dst) and cam not in failed_cameras(plan, state),
            "notes": c["notes"]}
    fits, fits_from = setup_in_use(plan, args)
    flt = read_stamp(out, "filter")
    mrg = read_stamp(out, "merge") or {}
    dstamp = read_stamp(out, "dem") or {}
    chain_ok, chain_why = chain_current(plan)
    comps = []
    for s in plan["surveys"]:
        name = f"{plan['date']}_{s['name']}"
        r = (state.get("comparisons") or {}).get(s["name"], {})
        cst = read_stamp(out, f"compare_{s['name']}") or {}
        lab = r.get("label") and r or (cst.get("label") or {})
        cur, why_not = comparison_current(plan, s)
        head = load_json(out / "compare" / f"{name}_comparison.json") if cur else None
        status = r.get("status") or ("built" if cur else ("not built" if not cst else f"stale ({why_not})"))
        if r.get("status", "").startswith("skipped (up to date)") and not cur:
            status = f"stale ({why_not})"
        comps.append({"name": s["name"], "survey_path": s["path"], "survey_type": s["survey_type"],
                      "why_not_current": None if cur else why_not,
                      "survey_date": s["survey_date"], "table_label": s["label"],
                      "label": (lab.get("label") or s["label"]) if (cur or r) else s["label"],
                      "label_reasons": lab.get("reasons"),
                      "overruled": lab.get("overridden") or None,
                      "override_reason": lab.get("override_reason"),
                      "suggested_by_survey_compare": (head or {}).get("suggested_label"),
                      "why": (head or {}).get("why") or s["why"], "status": status, "reason": r.get("reason"),
                      "current": bool(cur),
                      "photo_dates_used": lab.get("photo_dates_used"),
                      "headline": ({k: head.get(k) for k in ("n", "median", "nmad", "rmse", "mean", "p5", "p95",
                                                             "dem_cells", "dem_cells_compared", "points",
                                                             "points_compared", "time_gap", "by_camera",
                                                             "waterlines", "days_outside_photo_window",
                                                             "days_from_photo_window_middle", "photo_dates",
                                                             "dem_cells_compared_partial_cover",
                                                             "points_outside_dem", "points_in_empty_dem_cells",
                                                             "vertical_datum_note")}
                                   if head else None)})
    stamps = {}
    for p in sorted((out / "logs" / "stamps").glob("*.json")) if (out / "logs" / "stamps").exists() else []:
        d = load_json(p) or {}
        if d.get("commands"):
            stamps[p.stem] = d["commands"]
    # options of this run that the outputs on disk do not reflect: a window kept from the build
    # (keep_built_window), input files whose step did not run or failed
    not_applied = list(state.get("not_applied") or [])
    win_diff = []
    if final:
        ina = inputs_not_applied(plan, args)
        if ina:
            not_applied.append("this run's input file(s) NOT applied (the step that reads them failed or did not "
                               "run): " + "; ".join(ina) + ". The outputs on disk come from the files in brackets; "
                               "the 'rebuild everything' command names those")
        win_diff = window_differences(plan, built_settings(plan))
    commands = {"survey_products": cmd_text([sys.executable, Path(__file__).resolve()]
                                            + state_argv(built_args(args, plan), plan["date"]))}
    if fits_from.startswith("the waterlines") and \
            {round(float(f["coef"]), 6) for f in fits} != {round(float(args.setup_coef), 6)}:
        commands["note"] = (f"this run asked for C = {args.setup_coef}, but the waterlines on disk were built with "
                            f"C = {', '.join(str(f['coef']) for f in fits)} (see the detect commands below); the "
                            f"first line rebuilds this product as it is (C = "
                            f"{', '.join(str(f['coef']) for f in fits)}); to apply C = {args.setup_coef}, rebuild "
                            f"from the detection with it")
    order = {st: i for i, st in enumerate(("forcing", "pointing", "detect", "filter", "dem", "map", "compare"))}
    for k in sorted(stamps, key=lambda k: (order.get(k.split("_")[0], 99), k)):
        commands[k if not k.startswith("compare_") else "comparison_" + k[8:]] = stamps[k]
    dem_want = dem_settings_requested(args)
    dem_built = dem_settings_built(out)
    dem_set = dem_built or dem_want
    dem_diff = dem_settings_differ(dem_built, dem_want)
    if dem_diff:
        commands["note_dem"] = (f"this run asked for other DEM settings than the DEM on disk was built with ({dem_diff}): "
                                f"the grids, README and provenance are those of the build; rerun --steps dem,maps,"
                                f"compare to apply this run's")
    if not_applied:
        commands["note_not_applied"] = ("the first line rebuilds the product ON DISK (its window, hours and input "
                                        "files); this run's options that were not applied: " + " | ".join(not_applied))
    failed_cams = failed_cameras(plan, state)
    missing = missing_parts(plan, state, args) if final else []
    status = ("in progress" if not final else
              "failed" if state.get("interrupted") else
              "complete" if not state.get("failures") and not failed_cams and chain_ok and not missing
              and not win_diff else "partial")
    info = load_json(out / "dem" / f"{plan['date']}_info.json")
    prov = {
        "date": plan["date"], "era": plan["era"], "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "status": status, "synthetic": product_synthetic(plan),
        "status_why": ([f"build INTERRUPTED ({state['interrupted']}): the steps after it did not run, and the "
                        f"one it was in may be half done; rerun the same command"] if state.get("interrupted") else [])
                      + ([f"{c}: {w}" for c, w in failed_cams.items()] + list(state.get("failures", []))
                       + ([] if chain_ok else [f"outputs on disk are not one build: {chain_why}"]) + missing
                       + ([f"the outputs on disk were not all built for this window: {'; '.join(win_diff)}"]
                          if win_diff else [])),
        "not_applied": not_applied or None,
        "steps_this_run": state.get("steps", {}), "failures": state.get("failures", []),
        "cameras_failed": failed_cams or None,
        "window": {"first_day": plan["first"], "last_day": plan["last"]},
        "cameras": cams, "cameras_not_built": plan.get("disabled_cams"),
        "water_level": (forcing or {}).get("water_level"),
        "waves": (forcing or {}).get("waves"),
        "forcing_files": {"dir": str(out / "forcing")},
        "setup": {"coef": (fits[0]["coef"] if len(fits) == 1 else [f["coef"] for f in fits]),
                  "from": fits_from,
                  "requested_this_run": args.setup_coef,
                  "fits": fits,
                  "formula": "beach elevation = still-water level + C*sqrt(Hs*L0), L0 = g*Tp^2/(2*pi)",
                  "applied_by": ("extract_elevation_contours.py --setup-coef (inside detect_original_view.py)"
                                 if plan["era"] != "live" else
                                 "survey_products.py apply_setup(): the same formula and rounding as "
                                 "extract_elevation_contours.py, on the archive rows' offshore_hs_m/offshore_tp_s"),
                  "fitted_to_survey": next((f.get("survey_stem") for f in fits if f.get("survey_stem")), None),
                  "fitted_how": "; ".join(f"C = {f['coef']}: {f.get('how')}" for f in fits),
                  "wave_currency": next((f.get("wave_currency") for f in fits if f.get("wave_currency")), None),
                  "still_water_reference": (
                      {"c_fitted_with": next((f.get("still_water") for f in fits if f.get("still_water")), None),
                       "this_date": ("GNSS-R (the same)" if plan["era"] == "live" else
                                     "ADCP at 21 m / Chatham transfer: outside the surf zone, no setup"),
                       "expected": (None if plan["era"] == "live" else
                                    "about +0.02 m - share x (setup - the mean setup of the GNSS-R/Chatham comparison "
                                    "period): the mean part of the share is in the datum chain, only the "
                                    "wave-dependent part shows"),
                       "measured": ((forcing or {}).get("waves") or {}).get("still_water_reference")}),
                  "frames_left_out_without_setup": mrg.get("no_setup_frames") or {},
                  "implied_by_rows": {c: v.get("implied_coef") for c, v in (mrg.get("setup") or {}).items()}},
        "filter": {"script": "waterline_consistency.py", "settings": dict({"options": "its defaults (as the live cron)"},
                                                                           **((flt or {}).get("settings") or {})),
                   "result": (flt or {}).get("summary"), "failed": bool(state.get("filter_failed")),
                   "file": str(waterline_file(plan, quiet_=True))},
        "dem": dict(dem_set, script="dem_from_contours.py",
                    settings_from=("the DEM on disk (its build stamp)" if dem_built else
                                   "this run's options (no DEM built yet)"),
                    requested_this_run=dem_want,
                    files={k: f"{out / 'dem' / plan['date']}_{k}" for k in ("dem", "spread", "count")},
                    cells=dstamp.get("cells"),
                    built_from=info,
                    geotiff=(f"asc_to_geotiff.py: tagged EPSG:{dem_set.get('geotiff_epsg')}"
                             + (" (WGS 84 / UTM 19N), NOMINAL" if dem_set.get("geotiff_epsg") == 32619 else "")
                             + "; the grid coordinates are in the GCP frame, assumed NAD83(2011) / UTM zone 19N "
                               "as the surveys (the GCP files state no CRS) "
                               "(EPSG:6348), untransformed; only the _dem GeoTIFF carries vertical keys (NAVD88, "
                               "EPSG:5703); spread, count and difference grids are not heights")),
        "chain": {"current": chain_ok, "why": chain_why or None},
        "comparisons": comps,
        "commands": commands,
        "software": software_of_outputs(out),
    }
    prov["frames_per_day"] = frames_per_day(waterline_file(plan, quiet_=True))
    prov["caveats"] = collect_caveats(plan, args, state, forcing, prov)
    tmp = out / "provenance.json.tmp"
    tmp.write_text(json.dumps(clean_json(prov), indent=1, default=str, allow_nan=False) + "\n")
    os.replace(str(tmp), str(out / "provenance.json"))
    if final:
        write_readme(plan, args, prov)
    return prov


def clean_json(o):
    """NaN and infinities -> None (JSON has no NaN)."""
    if isinstance(o, dict):
        return {k: clean_json(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean_json(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, np.integer):
        return int(o)
    return o


def waterline_readme_lines(wl_):
    """The README lines of a comparison's waterlines - survey headline: the statistic (with what it
    rests on), WHERE on the beach it was taken, by survey elevation, by camera, and the C = 0
    sensitivity read on the same frames or points both ways."""
    L = []
    unit = wl_.get("unit")
    L.append("    " + stat_text(wl_, "waterlines - survey")
             + (f" ({wl_.get('method')}: one value per frame)" if unit == "frames" else
                f", {wl_.get('frames')} frames"
                + (f", {wl_['survey_points']} distinct survey points" if wl_.get("survey_points") is not None else "")
                + f" ({wl_.get('method')})"))
    zu, za = wl_.get("line_elevation_compared"), wl_.get("line_elevation_all")
    if zu:
        L.append(f"      WHERE: the line points compared lie at {zu['p5']:+.2f} to {zu['p95']:+.2f} m (p5..p95; "
                 f"median {zu['median']:+.2f})"
                 + (f", the window's lines at {za['p5']:+.2f} to {za['p95']:+.2f} m (p5..p95)" if za else "")
                 + ": this check covers only that part of the beach")
    be = [b for b in wl_.get("by_elevation") or [] if b.get("n")]
    if len(be) > 1:
        def band(b):
            k, what = independent_count(b, unit)
            det = (f"n {b['n']} frames" if unit == "frames" else
                   f"n {b['n']}, {b.get('frames')} frames" + (f", {b['survey_points']} survey points"
                                                              if b.get("survey_points") is not None else ""))
            return f"{b['group']} m {b['median']:+.3f} ({det}{', too few ' + what if k < FEW_N else ''})"
        L.append("      by survey elevation: " + ", ".join(band(b) for b in be))
    line = per_camera_text(wl_.get("by_camera"), unit)
    if line:
        L.append("      by camera: " + line)
    fs = wl_.get("frame_set") or {}
    if fs:
        lo, al = fs.get("left_out_on_survey") or {}, fs.get("all_frames") or {}

        def one(st_):
            k, what = independent_count(st_)
            return (f"{st_['median']:+.3f} m (n {st_['n']}" + ("" if unit == "frames" else f", {st_.get('frames')} frames")
                    + (f", {st_['survey_points']} survey points" if st_.get("survey_points") is not None else "")
                    + (f", too few {what}" if k < FEW_N else "") + ")")
        L.append(f"      FRAMES: these are the {fs.get('dem_frames')} frames the DEM gridded"
                 + (f"; it left out {fs['left_out_frames']} ("
                    + ", ".join(f"{v} with {k}" for k, v in (fs.get("left_out_why") or {}).items()) + ")"
                    + (f", which read {one(lo)} on the survey; all frames of the window {one(al)}" if lo.get("n") else
                       ", none of them on the survey") if fs.get("left_out_frames") else "; it left out none")
                 + (f" [{fs['note']}]" if fs.get("note") else ""))
    ws = wl_.get("without_setup") or {}
    pr, unp = ws.get("paired") or {}, ws.get("unpaired") or {}
    if ws.get("n") and pr.get("n"):
        cams = [(c, v) for c, v in sorted((pr.get("by_camera") or {}).items()) if v.get("n")]

        def few_mark(v):
            k, what = independent_count(v)
            return f", too few {what}" if k < FEW_N else ""
        kp, whatp = independent_count(pr)
        L.append(f"      without the wave setup (C = 0, a sensitivity), on the SAME {pr['n']} {pr.get('unit')} both "
                 f"ways" + (f" ({kp} {whatp})" if whatp != pr.get("unit") and pr.get("n_independent") is not None else "")
                 + f": median {pr['without_setup_median']:+.3f} m without the setup, {pr['with_setup_median']:+.3f} m "
                 f"with it (effect of the setup {pr['setup_effect_median']:+.3f} m)"
                 + (f" -- TOO FEW ({kp} {whatp} < {FEW_N}): not an estimate" if kp < FEW_N else "")
                 + (("; by camera " + ", ".join(f"{c} {v['without_setup_median']:+.3f} without / "
                                                f"{v['with_setup_median']:+.3f} with (n {v['n']}{few_mark(v)})"
                                                for c, v in cams))
                    if len(cams) > 1 else ""))
        if unp.get("with_setup_only") or unp.get("without_setup_only"):
            L.append(f"        {unp.get('with_setup_only', 0)} {pr.get('unit')} of the headline have no C = 0 value "
                     f"({unp.get('why')})"
                     + (f"; their median with the setup {unp['with_setup_only_median']:+.3f} m"
                        if unp.get("with_setup_only_median") is not None else "")
                     + (f"; {unp['without_setup_only']} have one only at C = 0" if unp.get("without_setup_only") else "")
                     + f". All C = 0 values: median {ws['median']:+.3f} m, n {ws['n']}.")
        L.append(f"        ({ws.get('how')})")
    elif ws.get("n"):                    # a comparison made before the pairing (Oct 2026)
        L.append(f"      without the wave setup (C = 0, a sensitivity): median {ws['median']:+.3f} m"
                 + (f"; by camera {per_camera_text(ws.get('by_camera'), unit)}" if ws.get("by_camera") else "")
                 + f" ({ws.get('how')}; NOT on the same frames as the headline)")
    return L


def write_readme(plan, args, prov):
    out = plan["out"]
    d = plan["date"]
    L = []
    w = L.append
    w(f"SURVEY-DATE PRODUCT {d}  (built {prov['built_utc']}, status: {prov['status']})")
    w("=" * 72)
    if prov.get("synthetic"):
        w(f"SYNTHETIC FIXTURE: {SYNTHETIC_NOTE}. Built from synthetic test inputs ({prov['synthetic']}): every "
          f"number below tests the code, none is a result about the real beach.")
    if prov["status"] != "complete":
        for x in prov.get("status_why") or []:
            w(f"  NOT COMPLETE: {x}")
    for x in prov.get("not_applied") or []:
        w(f"  NOT APPLIED: {x}")
    w("")
    w("WHAT THIS IS")
    fpd = prov.get("frames_per_day") or {}
    used = sorted(c for c, v in fpd.items() if v)
    days = sorted({dd for v in fpd.values() for dd in v})
    w(f"An intertidal beach DEM of Marconi Beach made only from the station's camera photos "
      f"({', '.join(used) or 'no camera'}"
      + (f"; {', '.join(c for c in plan['cams'] if c not in used)} contributed nothing" if
         [c for c in plan["cams"] if c not in used] else "") + f"), to check against the survey(s) of {d} "
      f"listed below. Photo window {plan['first']} .. {plan['last']}; waterline frames on "
      + (f"{days[0]} .. {days[-1]} ({len(days)} days)" if days else "no day") + ".")
    fits_ = (prov.get("setup") or {}).get("fits") or []
    no_setup = bool(fits_) and all(f.get("kind") == "none" or not f.get("coef") for f in fits_)
    w("Each photo's waterline is where the water met the sand; its elevation is the measured water level "
      + ("only: NO wave setup (C = 0), so a line marked by the swash reads LOW, by roughly "
         f"{SHARE_KEPT_RANGE[0]}-{SHARE_KEPT_RANGE[1]} x the setup (the line also moves seaward onto lower beach). "
         if no_setup else
         f"plus the wave setup (C = {(prov.get('setup') or {}).get('coef')}). ")
      + "Many waterlines at different tides make the DEM.")
    w("")
    w("COMPARISON AND ITS LABEL (how independent the check is)")
    for c in prov["comparisons"]:
        h = c.get("headline") or {}
        if c["status"].startswith("skipped") and not h:
            w(f"  {c['name']}: NOT COMPARED. {c.get('reason') or ''}")
            continue
        if not c.get("current"):
            wn = c.get("why_not_current") or ""
            w(f"  {c['name']}: NO CURRENT COMPARISON ({c['status']}"
              + (f"; {wn}" if wn and wn not in c["status"] else "") + ").")
            continue
        w(f"  {c['name']} ({c['survey_type']}, {c['survey_date']}): {c['label']}")
        w(f"    why: {c['why']}")
        if c.get("overruled"):
            w(f"    downgrades OVERRULED in surveys.csv ({c.get('override_reason')}): {'; '.join(c['overruled'])}")
        if c.get("suggested_by_survey_compare") and c["suggested_by_survey_compare"] != c["label"]:
            w(f"    survey_compare.py's own checks suggest {c['suggested_by_survey_compare']}")
        if h.get("n") is not None:
            w("    " + stat_text(h) + " (positive = DEM too high)")
            line = per_camera_text(h.get("by_camera"))
            if line:
                w("      by camera: " + line)
        if camera_split(h):
            w("    CAMERAS DIFFER: " + camera_split(h))
        if camera_split_note(h):
            w("    cameras: " + camera_split_note(h))
        wl_ = h.get("waterlines") or {}
        if wl_.get("n"):
            for ln in waterline_readme_lines(wl_):
                w(ln)
        if h.get("time_gap"):
            w(f"    {h['time_gap']} (frames used; product date {d})")
    w("  Labels: INDEPENDENT = nothing in the chain was fitted to or placed with this survey; "
      "CROSS-VALIDATED = fitted on other data from the same source, held out here; PARTLY-CIRCULAR = some "
      "step used this survey; CIRCULAR = pointing or offsets fitted to this survey. A comparison is never "
      "labelled better than its weakest step.")
    w("")
    w("KNOWN CAVEATS")
    for c in prov["caveats"]:
        w("  - " + c)
    w("")
    w("WHAT WAS USED")
    for cam, c in prov["cameras"].items():
        cal = c["calibration"]
        hist = (cal.get("history") or {}).get("note") or ""
        w(f"  {cam}: station {c['station']}, photos {c['window'][0]} .. {c['window'][1]} {c['utc_hours']} UTC"
          + (f", {c['photos']['count']} photos (list: {c['photos']['list_file']})" if c.get("photos") else ""))
        fc = fpd.get(cam) or {}
        w("      waterline frames used, per day: "
          + (", ".join(f"{dd[5:]} {n}" for dd, n in fc.items()) + f" (total {sum(fc.values())})" if fc else "none"))
        w(f"      calibration {cal['eo_file']}" + (f" ({hist})" if hist else ""))
        w(f"      search envelope: {c['search_envelope']['survey'] or c['search_envelope']['how']}")
        pc_ = c.get("pointing_check") or {}
        if pc_.get("per_day"):
            verdicts = {}
            for r in pc_["per_day"]:
                verdicts.setdefault(r["verdict"], []).append(r["date"])
            w("      pointing check: " + "; ".join(f"{k} {', '.join(v)}" for k, v in verdicts.items()))
        ref = pc_.get("reference") or {}
        if ref.get("window_horizon_dtilt_deg") is not None and ref.get("clear_horizon_days"):
            w(f"      sea horizon vs the calibration over the window: tilt/roll {ref['window_horizon_dtilt_deg']:+.2f}/"
              f"{ref['window_horizon_droll_deg']:+.2f} deg"
              + (f" ({ref['window_horizon_offset_px']:.1f} px)" if ref.get("window_horizon_offset_px") is not None else "")
              + f" on {ref['clear_horizon_days']} clear day(s)"
              + (f" (above {OFFSET_REPORT_DEG} deg: DEM sensitivity under KNOWN CAVEATS)"
                 if max(abs(ref["window_horizon_dtilt_deg"]), abs(ref["window_horizon_droll_deg"])) > OFFSET_REPORT_DEG
                 else ""))
    wl, wv = prov.get("water_level") or {}, prov.get("waves") or {}
    srcs = wl.get("sources_used") or {}
    if plan["era"] != "live" and srcs and set(srcs) == {"adcp"}:
        w(f"  water level: the ADCP (Signature 1000) measured every hour used ({srcs['adcp']} rows); the "
          f"Chatham transfer ({wl.get('method')}) was not needed")
    else:
        w(f"  water level: {wl.get('method')}; rows by source {srcs}"
          + (f"; transfer CV RMS {wl.get('cv_rms_m')} m, extrapolation RMS {wl.get('extrapolation_rms_m')} m"
             if any('chatham' in str(k) for k in srcs) else ""))
    w(f"  waves: {wv.get('method') or wv.get('hs_currency')}; sources {wv.get('sources_used') or wv.get('source')}")
    st = prov["setup"]
    w(f"  wave setup: C = {st['coef']} ({st['formula']}); applied by {st['applied_by']}; from {st['from']}"
      + (f" (this run asked for {st['requested_this_run']})"
         if st["from"].startswith("the waterlines") and st["coef"] != st["requested_this_run"] else ""))
    w(f"      C fitted: {st['fitted_how']}")
    swr = st.get("still_water_reference") or {}
    if swr.get("c_fitted_with") and plan["era"] != "live":
        m = swr.get("measured") or {}
        w(f"      C's still-water reference: {swr['c_fitted_with']}; this date's: {swr.get('this_date')}; "
          f"expected {swr.get('expected')}"
          + (f"; measured share of the setup in GNSS-R {m['share']:.2f} (wave-dependent; ~{m.get('typical_m', 0):.2f} m "
             f"at the 2026 median waves, the part already in the datum chain)"
             if m.get("quantified") else "; share not measured here (see KNOWN CAVEATS)"))
    if st.get("frames_left_out_without_setup"):
        w("      frames left out, no wave record: " + ", ".join(f"{k} {v}" for k, v in
                                                               st["frames_left_out_without_setup"].items()))
    fs = prov["filter"]["settings"] or {}
    w(f"  filter: {prov['filter']['script']} (" + "; ".join(f"{k} {v}" for k, v in fs.items()) + "): "
      + ("FAILED: the unfiltered waterlines were used" if prov["filter"]["failed"] else
         (prov["filter"]["result"] or "not run yet")))
    dm = prov["dem"]
    w(f"  DEM: {dm['script']} cell {dm['cell_m']} m, min {dm['min_points']} frames per cell, max spread "
      f"{dm['max_spread_m']} m, max Hs {dm['max_hs_m']} m, max day offset {dm['max_day_offset_m']} m"
      + (" (day test NOT run)" if dm.get("day_test_off") else "") + f" (from {dm.get('settings_from')})"
      + (f"; this run asked for other DEM settings ({dem_settings_differ(dm, dm.get('requested_this_run') or {})}): "
         f"NOT applied, rerun --steps dem,maps,compare" if dem_settings_differ(dm, dm.get("requested_this_run") or {})
         else ""))
    if dm.get("day_test_off"):
        w(f"      day test not run: {dm['day_test_off']}")
    cl = dm.get("cells") or {}
    if cl:
        w(f"      cells the waterlines crossed {cl.get('with_points')}: filled {cl.get('filled')}, blanked "
          f"{cl.get('blanked')}" + (f" ({cl['blanked_why']})" if cl.get("blanked_why") else ""))
    bf = dm.get("built_from") or {}
    if bf.get("frames"):
        flt_ = bf.get("filters") or {}
        w(f"      gridded {sum(bf['frames'].values())} frames ("
          + ", ".join(f"{k} {v}" for k, v in sorted(bf["frames"].items())) + f") on {bf.get('days')} days; "
          f"left out before gridding: {flt_.get('rough_frames', 0)} frames with Hs > {flt_.get('max_hs')} m, "
          + (f"{len(flt_.get('camera_days_rejected') or [])} camera-days off the rest by > "
             f"{flt_.get('max_day_offset')} m" if flt_.get("max_day_offset") else "no day test"))
    w(f"      GeoTIFFs: {dm['geotiff']}")
    sw = prov["software"]
    w(f"  software: {sw.get('note')} (per step in provenance.json software.by_step)")
    w("")
    w("FILES")
    w("  forcing/      water_level.csv, waves.csv, forcing_report.txt, forcing.png, forcing.json")
    w("  pointing/     per-day pointing check per camera")
    w("  waterlines/   <cam>/ (detection work folders), contour_points_ground.csv (all), "
      "contour_points_ground_filtered.csv (after the consistency filter), consistency_report.csv, "
      "no_setup_frames.csv (frames left out: no wave record)")
    w(f"  dem/          {d}_dem/_spread/_count .asc and .tif, {d}_dem.png (the DEM page)")
    w("  maps/         waterlines drawn on the photos, per camera (daily_elevation_map.py, as the station's own "
      "maps: its colour bar says which elevation colours the lines); plan view over the survey (coloured by "
      "still water + setup, the elevation the DEM uses)")
    w("  compare/      DEM - survey grids/points, report, figure, json")
    w("  provenance.json  every source above, machine-readable")
    w("")
    w("EXACT COMMANDS (rebuild everything: the first line; each step's own command below it)")
    for k, v in prov["commands"].items():
        if isinstance(v, list):
            for c in v:
                w(f"  [{k}] {c}")
        else:
            w(f"  [{k}] {v}")
    (out / "README.txt").write_text("\n".join(L) + "\n")


# ---------------------------------------------------------------------
# One date, all dates, summary
# ---------------------------------------------------------------------

def dry_run(plan, args, steps):
    rule(f"{plan['date']} (DRY RUN: nothing is written)")
    out = plan["out"]
    say("era", plan["era"])
    say("output", str(out))
    say("window", f"{plan['first']} .. {plan['last']}")
    for cam, why in (plan.get("disabled_cams") or {}).items():
        say(f"{cam}", f"not built (disabled): {why}")
    total = 0
    for cam, c in plan["cams"].items():
        say(f"{cam}", f"station {c['station']}, {c['first']} .. {c['last']} {c['hours_text']} UTC")
        say("  calibration", f"{c['eo_file']} -> {c['eo_path'] or 'NOT FOUND in ' + str(args.calibration)}")
        if plan["era"] != "live":
            say("  envelope", f"{c['envelope_survey']} -> "
                f"{c['envelope_path'] or 'NOT FOUND in ' + ', '.join(args.survey_dirs)}")
        photos, cnt = scan_photos(plan, cam, args)
        for r in cnt["missing_roots"]:
            say("  photo folder", f"not found: {r}")
        days = per_day(photos)
        say("  photos", f"{len(photos)} in the window and hours"
            + (f" ({cnt['duplicates']} duplicate names skipped, {cnt['other_station']} of another station)"
               if cnt["duplicates"] or cnt["other_station"] else ""))
        for dd in days_between(c["first"], c["last"]):
            print(f"      {dd}: {len(days.get(dd, []))}")
        total += len(photos)
    if plan["era"] == "live":
        p = Path(args.live_contours)
        say("live contours", f"{p} " + (f"({p.stat().st_size / 1e6:.0f} MB)" if p.exists() else "NOT FOUND"))
        say("GNSS-R spline", f"{args.gnssr_spline} " + ("" if Path(args.gnssr_spline).exists() else "(not here)"))
    else:
        say("forcing", cmd_text(forcing_cmd(plan, args)))
        # the files historical_forcing.py will use, so a missing one shows now, not after the detection
        try:
            found = forcing_inputs(args)
        except StepFailed as exc:
            found = {}
            say("  inputs", f"ERROR: {exc}")
        need = ({"chatham": "the Chatham gauge (else downloaded from NOAA when the period needs it)",
                 "wis": "WIS ST63064 waves (no download: required for the setup unless NDBC 44013 covers)",
                 "ndbc": "NDBC 44013 waves"} if plan["era"] == "chatham" else
                {"adcp": "the ADCP (water level and waves)", "adcp_navd88": "the ADCP on NAVD88"})
        for k, p in found.items():
            say(f"  {k}", (str(p) if p else "NOT FOUND") + (f"  <- {need[k]}" if k in need and not p else ""))
        if plan["era"] == "chatham" and args.setup_coef and not (found.get("wis") or found.get("ndbc")):
            warn("no WIS or NDBC wave file: the forcing step will FAIL (no waves for the setup); give --wis / "
                 "--ndbc or put the files in the Chelsea_calibration folder")
        from historical_forcing import setup_share_inputs
        sh = setup_share_inputs(getattr(args, "gnssr_spline", None), getattr(args, "gauge_csv", None))
        say("  setup share", "2026 record for the share of setup GNSS-R sees: "
            + ", ".join(f"{Path(p).name} {'here' if p and Path(p).exists() else 'NOT here'}" for p in sh))
    for s in plan["surveys"]:
        paths, note = resolve_survey(s, plan, args)
        v = Verdict(s)
        for p_ in (paths or [Path(s["path"] or s["name"])]):
            honest_label(s, plan, args, p_, fits=[setup_fit_info(args.setup_coef, args)], verdict=v)
        label, reasons = v.label, v.reasons
        say(f"survey {s['name']}", f"{paths[0] if paths else 'MISSING'}{' (' + note + ')' if note else ''}")
        if not paths:
            say("  then", missing_survey_message(s, plan, args, note))
        say("  label", label + ("" if label == s["label"] else f" (table: {s['label']}; " + "; ".join(reasons) + ")")
            + (f" [overruled in surveys.csv ({v.override}): {'; '.join(v.overridden)}]" if v.overridden else ""))
    est = 30 + (total * SEC_PER_PHOTO_DETECT if plan["era"] != "live" else 120) + 60 + 40 * len(plan["cams"]) + \
        40 * len(plan["surveys"])
    say("steps", ", ".join(steps))
    for st in steps:
        key = {"forcing": "forcing", "filter": "filter", "dem": "dem"}.get(st)
        if key:
            d = read_stamp(out, key)
            say(f"  {st}", "built before (" + d["finished_utc"] + ")" if d else "not built yet")
    say("expected", f"~{est / 60:.0f} min (detection ~{SEC_PER_PHOTO_DETECT:.0f} s per photo, measured off the "
        f"station; the NUC's own rate is printed by the detection's progress lines)")
    if overlap_note(est):
        say("  note", f"started now, {overlap_note(est)}")


IN_PROGRESS_MARK = "BUILD IN PROGRESS"


def mark_in_progress(plan, args):
    """Before any step can rebuild anything: provenance.json says status 'in progress' (every
    comparison not current) and README.txt starts with a BUILD IN PROGRESS / INTERRUPTED banner over
    the previous build's text. The end of the build rewrites both; a build that dies on the way
    (Ctrl-C, a full disk, a power cut) leaves the mark, so neither the README nor --summary can
    present the previous build's numbers as current over files that have since changed."""
    out = plan["out"]
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    cmd = cmd_text(["python3", "survey_products.py"] + state_argv(args, plan["date"]))
    note = (f"{IN_PROGRESS_MARK} since {now} (pid {os.getpid()} on {os.uname()[1] if hasattr(os, 'uname') else '?'}): "
            f"if no build of this date is running now, it was INTERRUPTED and the files on disk may be part "
            f"old, part new; rerun: {cmd}")
    prov = load_json(out / "provenance.json") or {}
    if prov:
        if not str(prov.get("status", "")).startswith("in progress"):
            prov["previous_status"] = prov.get("status")
        for c in prov.get("comparisons") or []:
            c["current"], c["why_not_current"] = False, "a build of this date is in progress or was interrupted"
    else:
        prov = {"date": plan["date"], "era": plan["era"], "comparisons": [], "synthetic": RUN.get("synthetic"),
                "window": {"first_day": plan["first"], "last_day": plan["last"]}}
    # every survey of the date has a row (not current), so --summary lists the date while it builds
    have = {c.get("name") for c in prov.get("comparisons") or []}
    for sv in plan.get("surveys") or []:
        if sv["name"] not in have:
            prov.setdefault("comparisons", []).append(
                {"name": sv["name"], "survey_path": sv.get("path"), "survey_type": sv["survey_type"],
                 "survey_date": sv.get("survey_date"), "table_label": sv["label"], "label": sv["label"],
                 "why": sv.get("why"), "status": "not built yet", "current": False,
                 "why_not_current": "a build of this date is in progress or was interrupted"})
    prov.update(status="in progress", status_why=[note], in_progress={"since_utc": now, "pid": os.getpid(),
                                                                       "command": cmd})
    tmp = out / "provenance.json.tmp"
    tmp.write_text(json.dumps(clean_json(prov), indent=1, default=str) + "\n")
    os.replace(str(tmp), str(out / "provenance.json"))
    rd = out / "README.txt"
    old = rd.read_text() if rd.exists() else ""
    k = old.find("\n" + "=" * 72 + " PREVIOUS BUILD\n")
    if old.startswith(IN_PROGRESS_MARK) and k >= 0:
        old = old[k + 1 + 72 + len(" PREVIOUS BUILD\n"):]          # one banner, not one per attempt
    banner = (f"{IN_PROGRESS_MARK} / INTERRUPTED -- {plan['date']}\n{note}\n"
              f"Below: the PREVIOUS build's README, for reference only. Its numbers and status may not match the "
              f"files on disk now; nothing in it is current until this build ends and rewrites this file.\n")
    tmp = out / "README.txt.tmp"
    tmp.write_text(banner + (("\n" + "=" * 72 + " PREVIOUS BUILD\n" + old) if old else ""))
    os.replace(str(tmp), str(rd))


LOCK_NAME = ".build.lock"


def host_name():
    try:
        import socket
        return socket.gethostname()
    except Exception:
        return os.uname()[1] if hasattr(os, "uname") else "?"


def build_running(pid):
    """Is process `pid` (on this computer) alive and a survey_products.py build? A pid reused by
    another program after a power cut is not one (Linux: /proc/<pid>/cmdline says)."""
    try:
        os.kill(int(pid), 0)
    except PermissionError:
        return True
    except (OSError, ValueError, TypeError):
        return False
    try:
        return b"survey_products" in Path(f"/proc/{int(pid)}/cmdline").read_bytes()
    except OSError:
        return True                        # no /proc: alive is all that can be told


def take_lock(out, args, date):
    """
    One build of a date at a time: two would share every work file (the archive rows, the
    detector's folders, the stamps) and mark each other 'failed'. The lock is a file created
    atomically (O_CREAT | O_EXCL) in the date folder, holding the pid, the computer and the command.
    A lock left by a build that is no longer running on this computer (killed, power cut) is taken
    over, and said so. -> (True, note) when this run holds it, (False, why) when another build does.
    """
    lock = Path(out) / LOCK_NAME
    me = {"pid": os.getpid(), "host": host_name(),
          "since_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
          "command": cmd_text(["python3", "survey_products.py"] + state_argv(args, date))}
    note = ""
    for _ in range(3):
        try:
            fd = os.open(str(lock), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            other = load_json(lock)
            try:
                age = time.time() - lock.stat().st_mtime
            except OSError:
                continue                   # released meanwhile: try again
            if other is None and age < 30:
                time.sleep(1.0)            # being written by a build that started this instant
                other = load_json(lock)
            if other is None:
                note = f"an unreadable lock file {lock} ({age / 60:.0f} min old) was removed"
            elif other.get("host") == me["host"] and not build_running(other.get("pid")):
                note = (f"the lock of an earlier build (pid {other.get('pid')}, since {other.get('since_utc')}) "
                        f"was left behind: that build is not running (killed or interrupted); taken over")
            elif other.get("host") != me["host"]:
                return False, (f"another computer ({other.get('host')}) holds {lock} (pid {other.get('pid')} since "
                               f"{other.get('since_utc')}: {other.get('command')}). If no build of {date} runs "
                               f"there, delete that file and run again")
            else:
                return False, (f"another build of {date} is RUNNING (pid {other.get('pid')} since "
                               f"{other.get('since_utc')}: {other.get('command')}). Two builds of one date "
                               f"would overwrite each other's files: wait for it to end (its README.txt says "
                               f"BUILD IN PROGRESS until then), or stop it, then run again")
            try:
                lock.unlink()
            except OSError:
                pass
            continue
        with os.fdopen(fd, "w") as f:
            f.write(json.dumps(me) + "\n")
        return True, note
    return False, f"could not take {lock} (another build keeps taking it)"


def release_lock(out):
    lock = Path(out) / LOCK_NAME
    if (load_json(lock) or {}).get("pid") == os.getpid():
        try:
            lock.unlink()
        except OSError:
            pass


def build_date(date, cfg, surveys, args, steps):
    args = argparse.Namespace(**vars(args))          # this date's own copy (built settings go in it)
    plan = build_plan(date, cfg, surveys, args)
    if plan is None:
        return "disabled"
    state = {"steps": {}, "failures": [], "not_applied": []}
    keep_built_window(plan, args, steps, state)
    if args.dry_run:
        dry_run(plan, args, steps)
        return "dry-run"
    rule(f"{date}  era {plan['era']}  cameras {', '.join(plan['cams'])}  photos {plan['first']} .. {plan['last']}")
    say("output", str(plan["out"]))
    state["argv"] = state_argv(args, date)
    try:
        plan["out"].mkdir(parents=True, exist_ok=True)
        ok, note = take_lock(plan["out"], args, date)
    except OSError as exc:
        warn(f"{date}: cannot write in {plan['out']} ({exc}): NOTHING rebuilt")
        return "failed"
    if not ok:
        warn(f"{date}: NOT built: {note}")
        return "busy"
    if note:
        warn(f"{date}: {note}")
    try:
        return run_date(plan, args, steps, state)
    finally:
        release_lock(plan["out"])


def run_synthetic(args):
    """This run's synthetic-input note: --synthetic, or a SYNTHETIC_FIXTURE file in a photo root."""
    if getattr(args, "synthetic", None):
        return str(args.synthetic)
    for r in getattr(args, "photo_roots", None) or []:
        m = Path(r) / SYNTHETIC_MARKER
        if m.exists():
            try:
                t = m.read_text().strip().splitlines()
            except (OSError, UnicodeDecodeError):
                t = []
            return (t[0] if t else "") or f"marker file {m}"
    return None


def product_synthetic(plan):
    """The product is (partly) synthetic when this run is, or any step on disk was built in a synthetic run."""
    if RUN.get("synthetic"):
        return RUN["synthetic"]
    sdir = plan["out"] / "logs" / "stamps"
    for p_ in sorted(sdir.glob("*.json")) if sdir.exists() else []:
        v = (load_json(p_) or {}).get("synthetic")
        if v:
            return v
    return None


def run_date(plan, args, steps, state):
    date = plan["date"]
    RUN["synthetic"] = run_synthetic(args)
    if RUN["synthetic"]:
        warn(f"{date}: SYNTHETIC inputs ({RUN['synthetic']}): this product {SYNTHETIC_NOTE}; every output says so")
    try:
        check_disk(plan["out"], 50, "marking the build in progress and writing its provenance")
        mark_in_progress(plan, args)
    except (OSError, StepFailed) as exc:
        warn(f"{date}: cannot mark the build in progress in {plan['out']} ({exc}): NOTHING rebuilt")
        return "failed"
    funcs = {"forcing": step_forcing, "pointing": step_pointing, "detect": step_detect, "filter": step_filter,
             "dem": step_dem, "maps": step_maps, "compare": step_compare}
    t_all = time.time()
    st = None
    try:
        for st in STEPS:
            if st not in steps:
                state["steps"][st] = "not requested"
                continue
            rule(f"{date}: {st}")
            t0 = time.time()
            try:
                state["steps"][st] = funcs[st](plan, args, state)
            except StepFailed as exc:
                state["steps"][st] = "FAILED"
                state["failures"].append(f"{st}: {exc}")
                warn(f"step {st} FAILED: {exc}")
                if st in ("forcing", "detect", "dem"):
                    warn(f"later steps need {st}: stopping {date} here")
                    break
            say(f"{st} done", f"{state['steps'][st]} ({time.time() - t0:.0f} s)")
    except BaseException as exc:
        # Ctrl-C, a full disk, a bug: the forcing or waterlines on disk may already be new while the
        # README and provenance.json (marked BUILD IN PROGRESS above) are not. Say so in both, as
        # 'failed', then stop as before. If even that cannot be written (disk full), the in-progress
        # mark stays, and --summary shows nothing of this date as current.
        state["steps"][st or "?"] = "INTERRUPTED"
        state["interrupted"] = f"{type(exc).__name__} in step {st}: {exc}".rstrip(": ")
        state["failures"].append(f"{st}: interrupted by {type(exc).__name__}: {exc}")
        warn(f"{date}: build INTERRUPTED in step {st} ({type(exc).__name__}: {exc})")
        # a second Ctrl-C (or the signal a process group gets twice from 'timeout') must not cut this
        # short: SIGINT is ignored for the few seconds it takes
        import signal
        try:
            old_int = signal.signal(signal.SIGINT, signal.SIG_IGN)
        except (ValueError, OSError):
            old_int = None
        try:
            write_provenance(plan, args, state, final=True)
            warn(f"README.txt and provenance.json say 'status: failed' (interrupted); rerun the same command")
        except BaseException as exc2:
            warn(f"could not write the final README/provenance ({exc2!r}): they still say BUILD IN PROGRESS / "
                 f"INTERRUPTED")
        finally:
            if old_int is not None:
                signal.signal(signal.SIGINT, old_int)
        raise
    if state["failures"]:
        keep_built_window(plan, args, (), state, why="a step failed before the outputs were rebuilt for it")
    prov = write_provenance(plan, args, state, final=True)
    rule(f"{date}: summary of this run ({time.time() - t_all:.0f} s)")
    for st in STEPS:
        say(st, state["steps"].get(st, "not reached"))
    for c in prov["comparisons"]:
        h = c.get("headline") or {}
        if c.get("current") and h.get("n") is not None:
            say(c["name"], f"{c['label']}: " + stat_text(h))
            if per_camera_text(h.get("by_camera")):
                say("  by camera", per_camera_text(h.get("by_camera")))
            if camera_split(h):
                warn(f"{c['name']}: {camera_split(h)}")
            if camera_split_note(h):
                say("  cameras", camera_split_note(h))
            wl_ = h.get("waterlines") or {}
            if wl_.get("n"):
                say("  waterlines", stat_text(wl_, "waterlines - survey") + "; by camera "
                    + (per_camera_text(wl_.get("by_camera"), wl_.get("unit")) or "-"))
                zu = wl_.get("line_elevation_compared")
                if zu:
                    za = wl_.get("line_elevation_all") or {}
                    say("  where", f"line points compared at {zu['p5']:+.2f} .. {zu['p95']:+.2f} m only (p5..p95)"
                        + (f"; the window's lines {za['p5']:+.2f} .. {za['p95']:+.2f} m" if za.get("p95") is not None
                           else ""))
            ws = wl_.get("without_setup") or {}
            pr = ws.get("paired") or {}
            if pr.get("n"):
                kp, whatp = independent_count(pr)
                say("  without setup", f"C = 0 sensitivity on the same {pr['n']} {pr.get('unit')}"
                    + (f" (TOO FEW: {kp} {whatp})" if kp < FEW_N else "") + ": "
                    f"{pr['without_setup_median']:+.3f} m without, {pr['with_setup_median']:+.3f} m with the setup; "
                    + ", ".join(f"{c} {v['without_setup_median']:+.3f} / {v['with_setup_median']:+.3f}"
                                + (f" (too few: {independent_count(v)[0]} {independent_count(v)[1]})"
                                   if independent_count(v)[0] < FEW_N else "")
                                for c, v in sorted((pr.get("by_camera") or {}).items()) if v.get("n"))
                    + f"; {(ws.get('unpaired') or {}).get('with_setup_only', 0)} drop out at C = 0")
            elif str(wl_.get("setup", "")).startswith("NO setup"):
                say("  without setup", "no C = 0 sensitivity: these lines carry no setup (C = 0)")
            fs = wl_.get("frame_set") or {}
            if fs.get("left_out_frames"):
                lo = fs.get("left_out_on_survey") or {}
                say("  DEM frames", f"lines on the {fs.get('dem_frames')} frames the DEM gridded; the "
                    f"{fs['left_out_frames']} it left out " + (f"read {lo['median']:+.3f} m (n {lo['n']})" if lo.get("n")
                                                              else "are not on the survey"))
        else:
            say(c["name"], f"{c['label']}: {c['status']}")
    # what an operator must not miss, again at the end
    for cam, why in failed_cameras(plan, state).items():
        warn(f"{cam} contributed NOTHING to this product: {why} (status partial)")
    ns = (prov.get("setup") or {}).get("frames_left_out_without_setup") or {}
    if ns:
        warn("frames without a wave record (no setup) left out: " + ", ".join(f"{k} {v}" for k, v in ns.items()))
    unc = ((prov.get("waves") or {}).get("uncovered_spans")) or []
    if unc:
        warn("no wave record for: " + "; ".join(unc[:4]) + (" ..." if len(unc) > 4 else ""))
    if not prov["chain"]["current"]:
        warn(f"the outputs on disk are not one build ({prov['chain']['why']}): status partial"
             + ("" if any(c.get("current") for c in prov["comparisons"]) else "; no comparison is reported as current"))
    say("status", prov["status"] + ("" if prov["status"] == "complete" else
                                    " -- " + "; ".join(prov.get("status_why") or [])))
    say("provenance", str(plan["out"] / "provenance.json"))
    say("README", str(plan["out"] / "README.txt"))
    if state["failures"]:
        return "failed"
    return "built" if prov["status"] == "complete" else "partial"


def opt_after(words, name, many=False):
    """The value after option `name` in a command split into words (many: every value up to the next
    option), or None when the option is not there."""
    if name not in words:
        return None
    i = words.index(name)
    if not many:
        return words[i + 1] if i + 1 < len(words) else None
    vals = []
    for w in words[i + 1:]:
        if w.startswith("--"):
            break
        vals.append(w)
    return vals


def stamp_words(out, key, field="cmd"):
    try:
        return shlex.split(((read_stamp(out, key) or {}).get("sig") or {}).get(field) or "")
    except ValueError:
        return []


def built_settings(plan):
    """
    What the outputs ON DISK were built with, from their stamps (never this run's options): each
    camera's photo window and UTC hours (the detection's command; live: the archive rows read), the
    forcing period, and the input files (the forcing's command, the live contour file, the photo
    roots). -> {"cams": {cam: {first, last, hours}}, "forcing": (first, last) or None,
    "inputs": {option: value}}; empty where nothing is built.
    """
    out = plan["out"]
    res = {"cams": {}, "forcing": None, "inputs": {}}
    if plan.get("era") == "live":
        sig = (read_stamp(out, "archive_rows") or {}).get("sig") or {}
        for cam, v in (sig.get("cams") or {}).items():
            try:
                res["cams"][cam] = {"first": v[0], "last": v[1], "hours": (float(v[2][0]), float(v[2][1]))}
            except (TypeError, ValueError, IndexError):
                pass
        if not res["cams"]:
            # no archive-rows stamp (e.g. a product of an older version whose read failed): the DEM's
            # own window (dem_from_contours.py --start-date/--end-date) still says what the grids on
            # disk came from; the hours are not recorded there, so this run's are assumed
            w = stamp_words(out, "dem")
            a, b = opt_after(w, "--start-date"), opt_after(w, "--end-date")
            if a and b:
                for cam, c in plan["cams"].items():
                    if read_stamp(out, f"detect_{cam}"):
                        res["cams"][cam] = {"first": a, "last": b, "hours": tuple(float(x) for x in c["hours"])}
        if res["cams"]:
            res["forcing"] = (min(v["first"] for v in res["cams"].values()),
                              max(v["last"] for v in res["cams"].values()))
        if sig.get("src"):
            res["inputs"]["live_contours"] = sig["src"]
    else:
        w = stamp_words(out, "forcing")
        if w:
            a, b = opt_after(w, "--start"), opt_after(w, "--end")
            if a and b:
                res["forcing"] = (a, b)
            for k in FORCING_INPUTS:
                res["inputs"][k] = opt_after(w, "--" + k.replace("_", "-"))
            res["inputs"]["gnssr_spline"] = opt_after(w, "--gnssr-spline")
            res["inputs"]["gauge_csv"] = opt_after(w, "--gauge-archive")
            res["inputs"]["no_download"] = "--no-download" in w
        for cam in plan["cams"]:
            w = stamp_words(out, f"detect_{cam}")
            a, b, h = opt_after(w, "--start-date"), opt_after(w, "--end-date"), opt_after(w, "--utc-hours")
            if a and b and h:
                res["cams"][cam] = {"first": a, "last": b, "hours": parse_hours(h)}
                roots = opt_after(w, "--originals", many=True)
                if roots and "photo_roots" not in res["inputs"]:
                    res["inputs"]["photo_roots"] = roots
    if "photo_roots" not in res["inputs"]:
        for cam in plan["cams"]:
            r = (read_stamp(out, f"pointing_{cam}") or {}).get("photo_roots")
            if r:
                res["inputs"]["photo_roots"] = list(r)
                break
    return res


WINDOW_STEPS = ("forcing", "pointing", "detect", "filter", "dem")


def window_differences(plan, bs):
    """This run's photo window and hours against those the outputs on disk were built with. -> list
    of text ('' when they agree or nothing is built)."""
    diff = []
    for cam, c in plan["cams"].items():
        b = bs["cams"].get(cam)
        if not b:
            continue
        if (c["first"], c["last"]) != (b["first"], b["last"]):
            diff.append(f"{cam} window {c['first']} .. {c['last']} (built: {b['first']} .. {b['last']})")
        if tuple(float(x) for x in c["hours"]) != tuple(b["hours"]):
            diff.append(f"{cam} UTC hours {c['hours_text']} (built: {b['hours'][0]:g}-{b['hours'][1]:g})")
    f = bs.get("forcing")
    if f and not bs["cams"] and (plan["first"], plan["last"]) != tuple(f):
        diff.append(f"window {plan['first']} .. {plan['last']} (the forcing on disk: {f[0]} .. {f[1]})")
    return diff


def keep_built_window(plan, args, steps, state, why=""):
    """
    A product's photo window and hours are fixed by its waterlines (and forcing). A run asking for
    another --window / --utc-hours applies them only when it rebuilds from the forcing through the
    DEM (WINDOW_STEPS); otherwise (e.g. --steps compare --window ...) they are NOT applied: this run
    uses the window and hours the outputs on disk were built with, says so, and the README,
    provenance and 'rebuild everything' line report those.
    """
    bs = built_settings(plan)
    diff = window_differences(plan, bs)
    if not diff or set(WINDOW_STEPS) <= set(steps):
        return
    text = (f"this run's window/hours NOT applied ({'; '.join(diff)})" + (f" ({why})" if why else "")
            + f": the outputs on disk were built for the window shown; to apply them, rebuild from the forcing "
              f"(all steps: drop --steps, or --steps {','.join(STEPS)})")
    warn(text)
    state["not_applied"].append(text)
    f = bs.get("forcing")
    for cam, c in plan["cams"].items():
        b = bs["cams"].get(cam) or ({"first": f[0], "last": f[1], "hours": c["hours"]} if f else None)
        if not b:
            continue
        c["first"], c["last"] = b["first"], b["last"]
        c["hours"] = tuple(float(x) for x in b["hours"])
        c["hours_text"] = f"{c['hours'][0]:g}-{c['hours'][1]:g}"
    plan["first"] = min(c["first"] for c in plan["cams"].values())
    plan["last"] = max(c["last"] for c in plan["cams"].values())
    window_args(plan, args)


def window_args(plan, args):
    """--window / --utc-hours that give the plan's window and hours (None where the table's do)."""
    cams = list(plan["cams"].values())
    same_w = all((c["first"], c["last"]) == (c["table_first"], c["table_last"]) for c in cams)
    one_w = len({(c["first"], c["last"]) for c in cams}) == 1
    args.window = None if same_w else ([cams[0]["first"], cams[0]["last"]] if one_w else args.window)
    same_h = all(c["hours_text"] == c["table_hours_text"] for c in cams)
    one_h = len({c["hours_text"] for c in cams}) == 1
    args.utc_hours = None if same_h else (cams[0]["hours_text"] if one_h else args.utc_hours)


def inputs_not_applied(plan, args):
    """This run's input files (--live-contours, the forcing files, --photo-roots) that differ from
    those the outputs on disk were built from (the step that reads them did not run, or failed).
    -> list of text."""
    bs = built_settings(plan)
    out = []
    for k, v in bs["inputs"].items():
        if k == "no_download" or v is None and getattr(args, k, None) is None:
            continue
        now = getattr(args, k, None)
        if k == "photo_roots":
            if v and list(map(str, now or [])) != list(map(str, v)):
                out.append(f"--photo-roots {' '.join(map(str, now or []))} (outputs built from photo roots "
                           f"{' '.join(map(str, v))})")
        elif str(now or "") != str(v or ""):
            out.append(f"--{k.replace('_', '-')} {now or '(default)'} (outputs built from {v or 'the default'})")
    return out


def built_args(args, plan):
    """args with the settings the outputs ON DISK were built with in place of this run's: the DEM
    settings and GeoTIFF EPSG (dem stamp) and the setup coefficient (merge stamp). The 'rebuild
    everything' line must rebuild THIS product, not what a later --steps compare asked for."""
    b = argparse.Namespace(**vars(args))
    dem = dem_settings_built(plan["out"])
    if dem:
        for k, opt in (("cell_m", "dem_cell"), ("min_points", "dem_min_points"), ("max_spread_m", "dem_max_spread"),
                       ("max_hs_m", "dem_max_hs"), ("max_day_offset_m", "dem_max_day_offset")):
            if k in dem:
                setattr(b, opt, dem[k] if dem[k] is not None else 0.0)
        if dem.get("geotiff_epsg"):
            b.geotiff_epsg = int(dem["geotiff_epsg"])
    if read_stamp(plan["out"], "merge"):
        fits, _ = setup_in_use(plan, args)
        coefs = {f.get("coef") for f in fits}
        if len(coefs) == 1 and None not in coefs and abs(float(fits[0]["coef"]) - float(args.setup_coef)) > 1e-9:
            f = fits[0]
            b.setup_coef = float(f["coef"])
            b.setup_fitted_to = (None if f.get("kind") != "declared" else
                                 f.get("survey_name") or "none:" + str(f.get("how", "")).split(" (declared with")[0])
    # the window, hours and input files the outputs were built with (their stamps), not a later
    # run's that did not apply (a failed step, or one not run)
    bs = built_settings(plan)
    if bs["cams"] and len(bs["cams"]) == len(plan["cams"]):
        p2 = {"cams": {c: dict(v, first=bs["cams"][c]["first"], last=bs["cams"][c]["last"],
                               hours_text=f"{bs['cams'][c]['hours'][0]:g}-{bs['cams'][c]['hours'][1]:g}")
                       for c, v in plan["cams"].items()}}
        window_args(p2, b)
    for k, v in bs["inputs"].items():
        if k == "no_download":
            continue
        elif k == "photo_roots":
            if v:
                b.photo_roots = list(v)
        elif k in ("gnssr_spline", "gauge_csv"):
            if v:
                setattr(b, k, v)
        else:
            setattr(b, k, v)
    if "no_download" in bs["inputs"]:
        b.no_download = bool(bs["inputs"]["no_download"])
    return b


def state_argv(args, date):
    """The command line that rebuilds this date: every option that changes the result and is not
    at its default (the inputs, window and hours, setup, calibration, GNSS-R, pointing and DEM
    settings). --force, --steps, --nice and --keep-detector-debug change no result and are left out."""
    a = ["--date", date, "--output-root", str(args.output_root)]
    if Path(args.config).resolve() != DEFAULT_CONFIG.resolve():
        a += ["--config", str(args.config)]
    if Path(args.surveys).resolve() != DEFAULT_SURVEYS.resolve():
        a += ["--surveys", str(args.surveys)]
    if list(args.photo_roots) != DEFAULT_PHOTO_ROOTS:
        a += ["--photo-roots"] + [str(p) for p in args.photo_roots]
    if list(args.survey_dirs) != DEFAULT_SURVEY_DIRS:
        a += ["--survey-dirs"] + [str(p) for p in args.survey_dirs]
    for k in ("adcp", "adcp_navd88", "chatham", "wis", "ndbc"):
        if getattr(args, k):
            a += ["--" + k.replace("_", "-"), str(getattr(args, k))]
    if str(args.live_contours) != str(DEFAULT_LIVE_CONTOURS):
        a += ["--live-contours", str(args.live_contours)]
    if args.window:
        a += ["--window"] + list(args.window)
    if args.utc_hours is not None:
        a += ["--utc-hours", args.utc_hours]
    if args.setup_coef != SETUP_COEF:
        a += ["--setup-coef", str(args.setup_coef)]
    if getattr(args, "setup_fitted_to", None):
        a += ["--setup-fitted-to", args.setup_fitted_to]
    if args.geotiff_epsg != GEOTIFF_EPSG:
        a += ["--geotiff-epsg", str(args.geotiff_epsg)]
    if Path(args.calibration).resolve() != (HERE / "calibration").resolve():
        a += ["--calibration", str(args.calibration)]
    try:
        from marconi_water_level import GNSSR_SPLINE, GAUGE_CSV
    except Exception:
        GNSSR_SPLINE = GAUGE_CSV = None
    if args.gnssr_spline and str(args.gnssr_spline) != str(GNSSR_SPLINE):
        a += ["--gnssr-spline", str(args.gnssr_spline)]
    if getattr(args, "gauge_csv", None) and str(args.gauge_csv) != str(GAUGE_CSV):
        a += ["--gauge-csv", str(args.gauge_csv)]
    if str(args.live_cron) != str(LIVE_CRON):
        a += ["--live-cron", str(args.live_cron)]
    if str(args.pointing_log_dir) != str(HERE / "archive"):
        a += ["--pointing-log-dir", str(args.pointing_log_dir)]
    for opt, default in (("pointing_frames", 3), ("pointing_fail", POINT_FAIL_DEG), ("dem_cell", DEM_CELL),
                         ("dem_min_points", DEM_MIN_POINTS), ("dem_max_spread", DEM_MAX_SPREAD),
                         ("dem_max_hs", DEM_MAX_HS), ("dem_max_day_offset", DEM_MAX_DAY_OFFSET)):
        if getattr(args, opt) != default:
            a += ["--" + opt.replace("_", "-"), str(getattr(args, opt))]
    if args.keep_moved_days:
        a.append("--keep-moved-days")
    if args.force_disabled:
        a.append("--force-disabled")
    if args.no_download:
        a.append("--no-download")
    if getattr(args, "synthetic", None):
        a += ["--synthetic", str(args.synthetic)]
    return a


def strip_brackets(t):
    while re.search(r"\([^()]*\)", t):              # nested brackets: innermost first
        t = re.sub(r"\s*\([^()]*\)", "", t)
    return t


def short_reason(c):
    """One clause saying why a comparison has its label (for the summary figure and table): the
    downgrade that set it, else the first sentence of the table's reason."""
    lab = c.get("label") or ""
    for r in c.get("label_reasons") or []:
        if r.startswith(lab + ":"):
            t = strip_brackets(r[len(lab) + 1:].strip())
            return t.split(", so ")[0].split(": ")[0].strip()
    why = strip_brackets((c.get("why") or "").split(" | ")[0])
    return re.split(r"(?<=[a-z0-9)])\. |; ", why)[0].strip().rstrip(".")


def setup_short(prov):
    """The setup coefficient a product was built with and where it was fitted, in a few words."""
    st = prov.get("setup") or {}
    fits = st.get("fits") or []
    if not fits:
        return f"C = {st.get('coef')}" if st.get("coef") is not None else "C not known"
    parts = []
    for f in fits:
        k, c = f.get("kind"), f.get("coef")
        cams = f" ({'+'.join(f['cameras'])})" if len(fits) > 1 and f.get("cameras") else ""
        if k == "none" or not c:
            parts.append(f"C = 0: no wave setup{cams}")
        elif k == "known":
            parts.append(f"C = {c:g}, fitted to {f.get('survey_stem')}{cams}")
        elif k == "declared":
            parts.append(f"C = {c:g}, declared " + (f"fitted to {f.get('survey_name') or f.get('survey_stem')}"
                                                    if f.get("survey_stem") else
                                                    "fitted to no survey: " + str(f.get("how", "")).split(" (declared")[0])
                         + cams)
        else:
            parts.append(f"C = {c:g} of UNKNOWN origin{cams}")
    return "; ".join(parts)


def summary_reason(cc, prov):
    """Why a summary row has its label, from the build: the downgrade that set it, else (the table's
    label kept) what the setup coefficient was fitted to -- the part of the chain that decides the
    label of a survey C could have been fitted to -- then the table's own reason."""
    lab = cc.get("label") or ""
    if any(r.startswith(lab + ":") for r in cc.get("label_reasons") or []):
        return short_reason(cc)
    fits = (prov.get("setup") or {}).get("fits") or []
    if fits and all(f.get("kind") == "none" or not f.get("coef") for f in fits):
        lead = "C = 0: no setup, nothing fitted to this survey"
    elif fits and all(f.get("kind") == "known" for f in fits):
        lead = f"C fitted to {', '.join(sorted({str(f.get('survey_stem')) for f in fits}))}, not this survey"
    else:
        lead = ""
    t = short_reason(cc)
    return f"{lead}; {t}" if lead else t


def setup_line(cc, prov):
    """The summary figure's setup line: setup_short, and, for a label the table gave (no downgrade),
    whether C was fitted to this survey."""
    t = setup_short(prov)
    lab = cc.get("label") or ""
    if any(r.startswith(lab + ":") for r in cc.get("label_reasons") or []):
        return t
    fits = (prov.get("setup") or {}).get("fits") or []
    if fits and all(f.get("kind") == "none" or not f.get("coef") for f in fits):
        return t + ", nothing fitted to this survey"
    if fits and all(f.get("kind") == "known" for f in fits):
        return t + ", not this survey"
    return t


def wrap_phrases(text, width, max_lines):
    """text in lines of <= width, at most max_lines; cut at a phrase boundary (', ' '; ' ': ') with
    '...' when longer, never mid-word."""
    import textwrap
    lines = textwrap.wrap(text, width)
    if len(lines) <= max_lines:
        return lines
    head = " ".join(lines[:max_lines])
    cut = max(head.rfind(", "), head.rfind("; "), head.rfind(": "))
    head = (head[:cut] if cut > len(head) // 2 else head[:head.rfind(" ")]) + " ..."
    out = textwrap.wrap(head, width)
    while len(out) > max_lines:                       # the '...' can spill over: shorten a word at a time
        head = head[:head[:-4].rfind(" ")] + " ..."
        out = textwrap.wrap(head, width)
    return out


def num(v):
    """A finite float or None."""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if np.isfinite(f) else None


def summary_comparison(prov, pdir, c):
    """(current, why not, headline, label info) of one comparison of a product, from what is on DISK
    now: the stamps' chain (comparison_current) and the comparison's own json, never the 'current'
    flag a provenance.json recorded -- a build interrupted after it rewrote the forcing or the
    waterlines leaves an old provenance.json behind. A product whose provenance says a build is in
    progress (or was interrupted) or failed shows nothing as current."""
    date = prov.get("date") or pdir.name
    st = str(prov.get("status") or "")
    if st.startswith(("in progress", "failed")):
        return False, f"the product's build is {st} (README.txt says why)", {}, {}
    plan_min = {"date": date, "out": Path(pdir)}
    cur, why = comparison_current(plan_min, {"name": c["name"]})
    if not cur:
        return False, why, {}, {}
    head = load_json(Path(pdir) / "compare" / f"{date}_{c['name']}_comparison.json")
    if not head:
        return False, "its comparison.json is missing or unreadable", {}, {}
    lab = (read_stamp(Path(pdir), f"compare_{c['name']}") or {}).get("label") or {}
    return True, "", head, lab


def summary(args):
    """One table and one figure across every date built so far under --output-root."""
    root = Path(args.output_root)
    rows = []
    for pj in sorted(root.glob("*/provenance.json")):
        prov = load_json(pj) or {}
        for c in prov.get("comparisons", []):
            cur, why_not, h, lab = summary_comparison(prov, pj.parent, c)
            wl = h.get("waterlines") or {}
            ws = wl.get("without_setup") or {}
            pr = ws.get("paired") or {}
            ph = sum((v.get("photos") or {}).get("count", 0) for v in (prov.get("cameras") or {}).values())
            # the cameras and frames that made the DEM (the live era has no photo list of its own)
            fpd = prov.get("frames_per_day") or {}
            cams = sorted(k for k, v in fpd.items() if v) or sorted(
                k for k, v in (prov.get("cameras") or {}).items() if (v.get("photos") or {}).get("count"))
            total = h.get("dem_cells") if c["survey_type"] == "dsm" else h.get("points")
            done = h.get("dem_cells_compared") if c["survey_type"] == "dsm" else h.get("points_compared")
            cells = ((prov.get("dem") or {}).get("cells") or {})
            pd_ = h.get("photo_dates") or (lab.get("photo_dates_used") or c.get("photo_dates_used") or [None, None])[:2]
            mid = None
            if pd_ and pd_[0] and c.get("survey_date"):
                a, b = date_cls.fromisoformat(pd_[0]), date_cls.fromisoformat(pd_[1])
                mid = date_cls.fromisoformat(c["survey_date"]).toordinal() - (a.toordinal() + b.toordinal()) / 2.0
            dprod = ((date_cls.fromisoformat(c["survey_date"]) - date_cls.fromisoformat(prov["date"])).days
                     if c.get("survey_date") and prov.get("date") else None)
            n = h.get("n")
            # no label without a current comparison: a label is a statement about numbers shown
            label = (lab.get("label") or c.get("label") or "") if cur else "not current"
            cc = dict(c, label=label, label_reasons=lab.get("reasons") if cur else c.get("label_reasons"),
                      why=h.get("why") or c.get("why"))
            dem_cam = camera_medians(h.get("by_camera"))
            wl_cam = camera_medians(wl.get("by_camera"))
            wk, wunit = independent_count(wl, wl.get("unit")) if wl.get("n") else (0, "")
            zu = wl.get("line_elevation_compared") or {}
            row = {"date": prov.get("date"), "era": prov.get("era"), "status": prov.get("status"),
                   "cameras": "+".join(cams) or "none",
                   "frames": sum(sum(v.values()) for v in fpd.values()) if fpd else None,
                   "photos": ph, "window": f"{prov['window']['first_day']}..{prov['window']['last_day']}",
                   "frame_days": f"{pd_[0]}..{pd_[1]}" if pd_ and pd_[0] else None,
                   "survey": c["name"], "survey_type": c["survey_type"], "survey_date": c["survey_date"],
                   "label": label, "table_label": c.get("table_label") or "",
                   "label_reason": summary_reason(cc, prov) if cur else short_reason(cc),
                   "overruled": "; ".join((lab.get("overridden") if cur else c.get("overruled")) or []) or None,
                   "survey_compare_suggests": h.get("suggested_label") or "",
                   "comparison": (c.get("status") if cur else f"not current: {why_not}"),
                   "current": cur, "n": n, "median_m": num(h.get("median")),
                   "nmad_m": num(h.get("nmad")) if (n or 0) >= 3 else None,
                   "rmse_m": num(h.get("rmse")),
                   "p5_m": num(h.get("p5")) if (n or 0) >= 3 else None,
                   "p95_m": num(h.get("p95")) if (n or 0) >= 3 else None,
                   "too_few": bool(n is not None and n < FEW_N),
                   "compared": done, "of": total,
                   "coverage": (f"{done}/{total} {'DEM cells' if c['survey_type'] == 'dsm' else 'points'}"
                                if total else None),
                   "dem_cells_crossed": cells.get("with_points"), "dem_cells_filled": cells.get("filled"),
                   "waterlines_median_m": num(wl.get("median")),
                   "waterlines_nmad_m": num(wl.get("nmad")) if (wl.get("n") or 0) >= 3 else None,
                   "waterlines_n": wl.get("n"), "waterlines_frames": wl.get("frames"),
                   "waterlines_survey_points": wl.get("survey_points"),
                   "waterlines_n_independent": wk if wl.get("n") else None,
                   "waterlines_n_independent_unit": wunit or None,
                   "waterlines_too_few": bool(wl.get("n")) and wk < FEW_N,
                   "waterlines_method": wl.get("method"),
                   "waterlines_elevation_p5_m": num(zu.get("p5")), "waterlines_elevation_p95_m": num(zu.get("p95")),
                   "waterlines_elevation_max_m": num(zu.get("max")),
                   "c0_paired_n": pr.get("n"), "c0_paired_unit": pr.get("unit"),
                   "c0_without_setup_median_m": num(pr.get("without_setup_median")),
                   "c0_with_setup_median_m": num(pr.get("with_setup_median")),
                   "c0_dropped": (ws.get("unpaired") or {}).get("with_setup_only"),
                   "camera_split": camera_split(h) or None,
                   "camera_note": camera_split_note(h) or None,
                   "days_outside_frames": h.get("days_outside_photo_window"),
                   "days_survey_minus_frames_middle": round(mid, 1) if mid is not None else None,
                   "days_survey_minus_product_date": dprod,
                   "setup_coef": (prov.get("setup") or {}).get("coef"), "setup": setup_short(prov),
                   "setup_line": setup_line(cc, prov) if cur else setup_short(prov),
                   "label_reason_table": short_reason(cc),
                   "synthetic": prov.get("synthetic") or None,
                   "built_utc": prov.get("built_utc"), "why": cc.get("why")}
            for cam in ("c1", "c2"):
                row[f"median_{cam}_m"], row[f"n_{cam}"] = (dem_cam[cam] if cam in dem_cam else (None, None))
                row[f"waterlines_median_{cam}_m"], row[f"waterlines_n_{cam}"] = (wl_cam[cam] if cam in wl_cam
                                                                                  else (None, None))
                # what a camera's line statistic rests on (frames, or distinct survey points)
                v_ = (wl.get("by_camera") or {}).get(cam)
                k_, u_ = independent_count(v_, wl.get("unit")) if v_ else (None, None)
                row[f"waterlines_n_independent_{cam}"], row[f"waterlines_n_independent_unit_{cam}"] = k_, u_
                pc_ = (pr.get("by_camera") or {}).get(cam)
                row[f"c0_without_setup_median_{cam}_m"] = num((pc_ or {}).get("without_setup_median"))
                row[f"c0_with_setup_median_{cam}_m"] = num((pc_ or {}).get("with_setup_median"))
                row[f"c0_n_independent_{cam}"] = independent_count(pc_)[0] if pc_ and pc_.get("n") else None
            kp_, up_ = independent_count(pr) if pr.get("n") else (None, None)
            row["c0_n_independent"], row["c0_n_independent_unit"] = kp_, up_
            # the label of the C = 0 values: what the comparison would be without a setup coefficient
            row["c0_label"] = (label_without_setup(c.get("table_label") or label, cc.get("label_reasons"))
                               if cur and pr.get("n") else None)
            row["c0_too_few"] = bool(pr.get("n")) and kp_ < FEW_N
            row["dem_by_camera"] = {k: v[0] for k, v in dem_cam.items()}
            row["waterlines_by_camera"] = {k: v[0] for k, v in wl_cam.items()}
            rows.append(row)
    # the dates the configuration leaves out (the user asked to skip them): listed, so the summary
    # never reads as if they had been forgotten
    disabled, not_yet = [], []
    try:
        cfg = read_table(args.config, CONFIG_FIELDS)
        srv = read_table(args.surveys, SURVEY_FIELDS)
        built_dates = {r["date"] for r in rows}
        for d in sorted({r["date"] for r in cfg}):
            rr = [r for r in cfg if r["date"] == d]
            if d in built_dates:
                continue
            if any(r["enabled"] in ("1", "yes", "true", "True") for r in rr):
                # enabled, but no product (or no comparison row) under this root yet
                st_ = (load_json(root / d / "provenance.json") or {}).get("status")
                for sv in ([x for x in srv if x["date"] == d] or [{"name": "-", "survey_type": "", "survey_date": "",
                                                                      "label": ""}]):
                    lw, lw_why = label_when_built(d, sv, cfg, srv, args)
                    not_yet.append({"date": d, "survey": sv["name"], "survey_type": sv["survey_type"],
                                    "survey_date": sv["survey_date"], "table_label": sv["label"],
                                    "label_when_built": lw, "label_when_built_reason": lw_why,
                                    "comparison": (f"not built yet (no {d}/ under this root)" if not st_ else
                                                   f"no comparison row yet (product status: {st_})")})
                continue
            why = "; ".join(sorted({r["notes"].split(". ")[0] for r in rr if r["notes"]})) or "disabled"
            for sv in ([x for x in srv if x["date"] == d] or [{"name": "-", "survey_type": "", "survey_date": "",
                                                                  "label": ""}]):
                disabled.append({"date": d, "survey": sv["name"], "survey_type": sv["survey_type"],
                                 "survey_date": sv["survey_date"], "table_label": sv["label"],
                                 "comparison": f"not built: disabled in {Path(args.config).name} ({why})"})
    except SystemExit:
        pass
    if not rows and not disabled and not not_yet:
        print(f"no product with a comparison under {root} yet")
        return 1
    fields = [k for k in (rows[0].keys() if rows else ["date", "survey", "survey_type", "survey_date", "comparison"])
              if k not in ("dem_by_camera", "waterlines_by_camera")]
    for k in ("table_label", "label_when_built", "label_when_built_reason"):
        if k not in fields:
            fields.append(k)
    with open(root / "summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
        w.writerows(not_yet)
        w.writerows(disabled)

    def fmt(v, s="{:+.3f}"):
        if v is None or v == "" or (isinstance(v, float) and not np.isfinite(v)):
            return "-"
        return s.format(v)
    lines = [f"Survey-date products under {root} ({datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC)",
             "DEM - survey, positive = DEM too high. Label = how independent the check is (reason below).", "",
             "n       DEM cells (lidar) or survey points compared; '*' = fewer than 10: NOT an estimate",
             "NMAD    robust spread; '-' when n < 3",
             "lines   waterline - survey: median over the waterline values on the survey (frames on RTK",
             "        transects: one value per frame), with n/frames (/survey points of a points survey);",
             "        '*' = fewer than 10 frames, or 10 distinct survey points: NOT an estimate",
             "coverage compared / DEM cells with a value (dsm) or / survey points; cells: DEM cells the",
             "        waterlines crossed -> filled (the rest blanked: too few frames or spread)",
             "out d   days the survey lies OUTSIDE the days of the frames used (0 = inside)",
             "mid d   survey date minus the middle of the frames used; prod d = survey date minus product date",
             "status  of the product build (complete / partial / in progress / failed); every comparison is",
             "        re-checked against the files on disk: one that is not current shows no numbers",
             "C       the wave-setup coefficient the product was built with (0 = no setup); where it was",
             "        fitted is under 'Why each label'",
             "under each row: the medians per camera ('!' = cameras more than "
             f"{CAMERA_SPLIT_M:.2f} m apart, each on >= {FEW_N} values: the pooled median hides it),",
             "        with what each rests on (DEM cells; frames or distinct survey points for the lines;",
             "        'too few' below 10), where on the beach the lines were compared, and the C = 0",
             "        sensitivity on the same frames", "",
             f"{'date':<11} {'survey':<16} {'label':<16} {'C':>5} {'n':>6} {'median':>8} {'NMAD':>6} {'RMSE':>6} "
             f"{'lines':>7} {'n/frames':>12} {'out d':>5} {'mid d':>6} {'prod d':>6}  {'coverage':<20} "
             f"{'cells':<9} cameras frames status"]
    for r in rows:
        star = "*" if r["too_few"] else " "
        wst = "*" if r["waterlines_too_few"] else " "
        wn = (f"{r['waterlines_n']}/{r['waterlines_frames']}"
              + (f"/{r['waterlines_survey_points']}p" if r.get("waterlines_survey_points") is not None else "")
              if r["waterlines_n"] else "-")
        cells = (f"{r['dem_cells_crossed']}->{r['dem_cells_filled']}" if r["dem_cells_crossed"] else "-")
        cval = r.get("setup_coef")
        ctxt = ("/".join(f"{float(x):g}" for x in cval) if isinstance(cval, list) else
                f"{float(cval):g}" if cval is not None else "?")
        lines.append(f"{r['date']:<11} {r['survey'][:16]:<16} {(r['label'] or r['comparison'] or '')[:16]:<16} "
                     f"{ctxt:>5} "
                     f"{fmt(r['n'], '{:d}'):>5}{star} {fmt(r['median_m']):>8} {fmt(r['nmad_m'], '{:.3f}'):>6} "
                     f"{fmt(r['rmse_m'], '{:.3f}'):>6} {fmt(r['waterlines_median_m']):>7} {wn:>11}{wst} "
                     f"{fmt(r['days_outside_frames'], '{:d}'):>5} {fmt(r['days_survey_minus_frames_middle'], '{:+.1f}'):>6} "
                     f"{fmt(r['days_survey_minus_product_date'], '{:+d}'):>6}  {r['coverage'] or '-':<20} "
                     f"{cells:<9} {r['cameras']:<7} {fmt(r['frames'], '{:d}'):>6} {r['status']}")
        if r.get("synthetic"):
            lines.append(f"{'':<11} SYNTHETIC FIXTURE: tests the pipeline, not the beach ({r['synthetic']})")
        if not r["current"] or (r["comparison"] and not str(r["comparison"]).startswith(("built", "skipped (up to date)"))):
            lines.append(f"{'':<11} comparison: {r['comparison']}")
        bc = []
        if len(r["dem_by_camera"]) > 1:
            bc.append("DEM " + ", ".join(f"{k} {v:+.3f} (n {r.get('n_' + k)}"
                                         + (", too few" if (r.get("n_" + k) or 0) < FEW_N else "") + ")"
                                         for k, v in r["dem_by_camera"].items()))
        if len(r["waterlines_by_camera"]) > 1:
            def wl_cam(k):
                ki, ku = r.get("waterlines_n_independent_" + k), r.get("waterlines_n_independent_unit_" + k)
                det = (f"{ki} {ku}" if ki is not None and ku and ku != "values" else f"n {r.get('waterlines_n_' + k)}")
                return det + (", too few" if ki is not None and ki < FEW_N else "")
            bc.append("lines " + ", ".join(f"{k} {v:+.3f} ({wl_cam(k)})" for k, v in r["waterlines_by_camera"].items()))
        if bc:
            lines.append(f"{'':<11} {'!' if r['camera_split'] else ' '} by camera: " + "; ".join(bc))
        if r.get("camera_note"):
            lines.append(f"{'':<11}   cameras: {r['camera_note']}")
        if r.get("waterlines_elevation_p5_m") is not None:
            lines.append(f"{'':<11}   lines compared at {r['waterlines_elevation_p5_m']:+.2f} .. "
                         f"{r['waterlines_elevation_p95_m']:+.2f} m only (p5..p95)")
        if r.get("c0_paired_n"):
            cams0 = [k for k in ("c1", "c2") if r.get(f"c0_n_independent_{k}")]
            lines.append(f"{'':<11}   C = 0 sensitivity, the same {r['c0_paired_n']} {r['c0_paired_unit']}"
                         + (f" ({r['c0_n_independent']} {r['c0_n_independent_unit']})"
                            if r.get("c0_n_independent_unit") not in (None, r["c0_paired_unit"]) else "")
                         + f": {r['c0_without_setup_median_m']:+.3f} m without the setup vs "
                         f"{r['c0_with_setup_median_m']:+.3f} m with it"
                         + (" -- TOO FEW: not an estimate" if r.get("c0_too_few") else "")
                         + (f"; the C = 0 values are {r['c0_label']} (no setup coefficient fitted to this survey in them)"
                            if r.get("c0_label") and r.get("label") in LABELS
                            and LABELS.index(r["c0_label"]) < LABELS.index(r["label"]) else "")
                         + (f" ({r['c0_dropped']} drop out at C = 0)" if r.get("c0_dropped") else "")
                         + (("; by camera " + ", ".join(
                             f"{k} {r[f'c0_without_setup_median_{k}_m']:+.3f} / {r[f'c0_with_setup_median_{k}_m']:+.3f}"
                             + (" (too few)" if r[f"c0_n_independent_{k}"] < FEW_N else "") for k in cams0))
                            if len(cams0) > 1 else ""))
        if r["survey_compare_suggests"] and r["survey_compare_suggests"] != r["label"]:
            lines.append(f"{'':<11} survey_compare.py's checks suggest {r['survey_compare_suggests']}"
                         + (f" (overruled: {r['overruled']})" if r["overruled"] else ""))
    if not_yet:
        lines += ["", "Enabled in the configuration, no comparison under this root yet (its label if built with "
                      f"this run's options, C = {args.setup_coef:g}):"]
        for r in not_yet:
            lw, tl = r.get("label_when_built"), r.get("table_label")
            lines.append(f"  {r['date']} {r['survey']}: {r['comparison']}"
                         + (f"; table label {tl}" if tl else "")
                         + (f"; built with C = {args.setup_coef:g} it becomes {lw} ({r['label_when_built_reason']})"
                            if lw and lw != tl and r.get("label_when_built_reason") else
                            f"; its label when built: {lw}" if lw else ""))
    if disabled:
        lines += ["", "Not built (disabled at the user's request in the configuration):"]
        for r in disabled:
            lines.append(f"  {r['date']} {r['survey']}: {r['comparison']}"
                         + (f"; table label {r['table_label']} (the build may make it worse, never better)"
                            if r.get("table_label") else ""))
    lines += ["", "Why each label:"]
    for r in rows:
        if not r["current"]:
            lines.append(f"  {r['date']} {r['survey']}: no current comparison ({r['comparison']})")
            continue
        lines.append(f"  {r['date']} {r['survey']}: {r['label']} -- {r['label_reason']}")
        lines.append(f"      setup: {r.get('setup')}")
        if r.get("synthetic"):
            lines.append(f"      SYNTHETIC FIXTURE: tests the pipeline, not the beach ({r['synthetic']})")
        lines.append(f"      full reason: {r['why']}")
    (root / "summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    try:
        summary_figure(rows, root / "summary.png", disabled, not_yet)
        print(f"\nwrote {root / 'summary.csv'}, {root / 'summary.txt'}, {root / 'summary.png'}")
    except Exception as exc:
        warn(f"summary figure failed: {exc!r}")
    return 0


def label_when_built(date, sv, cfg, srv, args):
    """The label a not-yet-built comparison would get if built with this run's options (the C of
    --setup-coef, the table's envelopes and calibrations), as --dry-run reports it: never the table's
    label when the chain would make it worse. -> (label, short reason) ('', '' for no survey row)."""
    if not sv.get("label") or sv.get("name") in (None, "", "-"):
        return "", ""
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            plan = build_plan(date, cfg, srv, args)
            if not plan:
                return sv["label"], ""
            paths, _ = resolve_survey(sv, plan, args)
            v = Verdict(sv)
            for p_ in (paths or [Path(sv["path"] or sv["name"])]):
                honest_label(sv, plan, args, p_, fits=[setup_fit_info(args.setup_coef, args)], verdict=v)
    except (SystemExit, Exception) as exc:          # never let the summary fail on a table problem
        return sv["label"], f"not checked: {exc}"
    worst = [r for r in v.reasons if r.startswith(v.label + ":")]
    return v.label, (worst[0][len(v.label) + 1:].strip() if worst else "")


def label_without_setup(table_label, reasons):
    """The label a comparison would have with no setup coefficient (C = 0): the table's label made
    worse by every reason except the setup coefficient's (a paired C = 0 value of a comparison that is
    CIRCULAR only through C fitted to its survey is INDEPENDENT)."""
    lab = table_label if table_label in LABELS else "CIRCULAR"
    for r in reasons or []:
        k = str(r).split(":", 1)[0]
        if k in LABELS and "setup coefficient" not in str(r) and LABELS.index(k) > LABELS.index(lab):
            lab = k
    return lab


def summary_figure(rows, path, disabled=(), not_yet=()):
    """Dot plot: median DEM - survey per comparison, NMAD as the bar (n >= 3 only), the label and the
    reason for it written beside it; on a two-camera date each camera's median as a small labelled
    marker (a pooled median near zero can hide two cameras off in opposite directions); the dates
    not built (disabled) listed under it."""
    import textwrap
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rs = [r for r in rows if r["median_m"] is not None]
    if not rs:
        raise ValueError("no comparison with numbers")
    fig, ax = plt.subplots(figsize=(11.5, 2.6 + 1.2 * len(rs)))
    y = np.arange(len(rs))[::-1]
    med = np.array([r["median_m"] for r in rs], float)
    nm = np.array([np.nan if r["nmad_m"] is None else r["nmad_m"] for r in rs], float)
    lim = max(0.5, float(np.nanmax(np.abs(med) + np.nan_to_num(nm))) * 1.15)
    wl_med = np.array([np.nan if r["waterlines_median_m"] is None else r["waterlines_median_m"] for r in rs], float)
    wl_nm = np.array([np.nan if r["waterlines_nmad_m"] is None else r["waterlines_nmad_m"] for r in rs], float)
    if np.isfinite(wl_med).any():
        lim = max(lim, float(np.nanmax(np.abs(wl_med) + np.nan_to_num(wl_nm))) * 1.15)
    camv = [v for r in rs for v in list(r["dem_by_camera"].values()) + list(r["waterlines_by_camera"].values())]
    if camv:
        lim = max(lim, max(abs(v) for v in camv) * 1.15)
    ax.axvline(0, color=MUTED, lw=1, zorder=1)
    # A median of a handful of cells or points is not an estimate: drawn faint and hollow, and
    # said so beside it, so it cannot sit level with a comparison of a thousand cells. No
    # spread bar below 3 values (the NMAD of one value is 0, which would read as perfect).
    few_ = np.array([(r["n"] or 0) < FEW_N for r in rs])
    wl_few = np.array([bool(r.get("waterlines_too_few")) for r in rs])
    # A CIRCULAR or PARTLY-CIRCULAR comparison is not a test of what was fitted to (or placed with) its
    # survey: drawn grey and hollow, never in the colours of an independent result. Its paired C = 0
    # values (the same frames without the setup), when nothing else made it circular, are drawn
    # beside it as their own, independent mark. Rows of synthetic test inputs carry a hatched band and
    # a SYNTHETIC tag on their name.
    circ = np.array([r.get("label") in ("CIRCULAR", "PARTLY-CIRCULAR") for r in rs])
    c0_ok = np.array([bool(circ[k] and r.get("c0_label") and r.get("c0_without_setup_median_m") is not None
                           and LABELS.index(r["c0_label"]) < LABELS.index(r["label"]))
                      for k, r in enumerate(rs)])
    if c0_ok.any():
        lim = max(lim, max(abs(r["c0_without_setup_median_m"]) for k, r in enumerate(rs) if c0_ok[k]) * 1.15)
    synth = np.array([bool(r.get("synthetic")) for r in rs])
    for k, yy in enumerate(y):
        if synth[k]:
            ax.axhspan(yy - 0.45, yy + 0.45, facecolor="none", edgecolor=GRID, hatch="////", lw=0, zorder=0)
    any_cam = False
    for k, (yy, m_, e_) in enumerate(zip(y, med, nm)):
        c1_, c2_ = (MUTED, MUTED) if circ[k] else (SERIES_1, SERIES_2)
        ax.errorbar([m_], [yy + 0.12], xerr=None if not np.isfinite(e_) else [e_], fmt="o", color=c1_,
                    ecolor=c1_, capsize=0, elinewidth=1.2 if (few_[k] or circ[k]) else 2, ms=8,
                    mfc="white" if (few_[k] or circ[k]) else c1_, mew=1.6, alpha=0.45 if few_[k] else 1.0, zorder=3)
        if np.isfinite(wl_med[k]):
            ax.errorbar([wl_med[k]], [yy - 0.12], xerr=None if not np.isfinite(wl_nm[k]) else [wl_nm[k]], fmt="D",
                        mfc="white", mec=c2_, color=c2_, ecolor=c2_, elinewidth=1.5, capsize=0,
                        ms=6, mew=1.6, alpha=0.45 if wl_few[k] else 1.0, zorder=3)
        if c0_ok[k]:
            r = rs[k]
            ax.plot([r["c0_without_setup_median_m"]], [yy - 0.34], marker="s", ms=6.5, ls="", color=SERIES_2,
                    mfc="white" if r.get("c0_too_few") else SERIES_2, mec=SERIES_2, mew=1.6,
                    alpha=0.45 if r.get("c0_too_few") else 1.0, zorder=4)
            ax.annotate(f"C = 0: {r['c0_label']}", (r["c0_without_setup_median_m"], yy - 0.34), xytext=(7, 0),
                        textcoords="offset points", va="center", ha="left", fontsize=7, color=INK2, zorder=5)
        # each camera's median, small, named, on the same row (only where there are two, and not on a
        # row that is itself too few to be an estimate); faint where that camera has fewer than FEW_N
        for bc, yo, col, pre, skip in ((rs[k]["dem_by_camera"], 0.12, c1_, "n_", few_[k]),
                                       (rs[k]["waterlines_by_camera"], -0.12, c2_, "waterlines_n_", wl_few[k])):
            if len(bc) > 1 and not skip:
                any_cam = True
                # names to the outer side of each mark (left of the lower, right of the higher), so two
                # cameras that agree do not print their names over each other
                for j, (cam, v) in enumerate(sorted(bc.items(), key=lambda kv: kv[1])):
                    faint = (rs[k].get(pre + cam) or 0) < FEW_N
                    ax.plot([v], [yy + yo], marker="|", ms=11, mew=1.6, color=col, ls="", zorder=4,
                            alpha=0.45 if faint else 1.0)
                    ax.annotate(cam, (v, yy + yo), xytext=(-2 if j == 0 else 2, 5), textcoords="offset points",
                                ha="right" if j == 0 else "left", va="bottom", fontsize=7, color=INK2, zorder=5)
    handles = [plt.Line2D([], [], marker="o", ls="-", color=SERIES_1, mfc=SERIES_1, ms=7, lw=2,
                          label="DEM - survey (cells or points on the DEM)")]
    if np.isfinite(wl_med).any():
        handles.append(plt.Line2D([], [], marker="D", ls="-", color=SERIES_2, mfc="white", mec=SERIES_2, ms=6,
                                  lw=1.5, label="waterlines - survey (one value per frame on RTK transects, else per point)"))
    if any_cam:
        handles.append(plt.Line2D([], [], marker="|", ls="", color=INK2, ms=10, mew=1.6,
                                  label="one camera's median (named above it)"))
    if few_.any() or wl_few[np.isfinite(wl_med)].any():
        handles.append(plt.Line2D([], [], marker="o", ls="", color=SERIES_1, mfc="white", alpha=0.45, ms=7,
                                  label=f"faint: fewer than {FEW_N} compared (frames or survey points for the "
                                        f"lines), not an estimate"))
    if circ.any():
        handles.append(plt.Line2D([], [], marker="o", ls="-", color=MUTED, mfc="white", mec=MUTED, ms=7, lw=1.2,
                                  label="grey, hollow: CIRCULAR / PARTLY-CIRCULAR, fitted to (or placed with) this "
                                        "survey: not a test of what was fitted"))
    if c0_ok.any():
        handles.append(plt.Line2D([], [], marker="s", ls="", color=SERIES_2, mfc=SERIES_2, ms=6.5,
                                  label="waterlines - survey with C = 0 (no setup), the same frames: its own label "
                                        "beside it"))
    if synth.any():
        from matplotlib.patches import Patch
        handles.append(Patch(facecolor="white", edgecolor=GRID, hatch="////",
                             label="hatched: SYNTHETIC test fixture: tests the pipeline, not the beach"))
    leg = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.42 if len(rs) < 3 else -0.18),
                    fontsize=8, frameon=False, ncol=2)
    for k, (yy, r) in enumerate(zip(y, rs)):
        txt = r["label"] + (f"  (n {r['n']}{', too few' if few_[k] else ''})" if r["n"] is not None else "")
        if r.get("survey_compare_suggests") and r["survey_compare_suggests"] != r["label"]:
            txt += (f"; {r['survey_compare_suggests']} overruled in surveys.csv" if r.get("overruled") else
                    f"; checks suggest {r['survey_compare_suggests']}")
        if r.get("synthetic"):
            txt = "SYNTHETIC FIXTURE -- " + txt
        ax.annotate(txt, (lim, yy + 0.2), xytext=(8, 0), textcoords="offset points", va="center", ha="left",
                    fontsize=8.5, color=INK, fontweight="bold", annotation_clip=False)
        # the setup it was built with, why the label (whole phrases, two lines), then what to watch
        reason = wrap_phrases(r.get("setup_line") or r.get("setup") or "", 68, 2)
        reason += wrap_phrases(r.get("label_reason_table") or r.get("label_reason") or "", 68, 2)
        extra = []
        if r.get("synthetic"):
            extra.append("tests the pipeline, not the beach")
        if r.get("waterlines_too_few"):
            extra.append(f"lines: {r['waterlines_n_independent']} {r['waterlines_n_independent_unit']}, too few")
        if c0_ok[k]:
            extra.append(f"C = 0, the same {r.get('c0_paired_n')} {r.get('c0_paired_unit')}: "
                         f"{r['c0_without_setup_median_m']:+.2f} m, {r['c0_label']}"
                         + (" (too few)" if r.get("c0_too_few") else ""))
        if r.get("camera_split"):
            extra.append(f"cameras disagree by more than {CAMERA_SPLIT_M:.1f} m: see the per-camera marks")
        reason += extra
        ax.annotate("\n".join(reason), (lim, yy + 0.04), xytext=(8, 0), textcoords="offset points", va="top",
                    ha="left", fontsize=7.6, color=INK2, annotation_clip=False, linespacing=1.15)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['date']}  {r['survey']}" + ("\nSYNTHETIC FIXTURE" if r.get("synthetic") else "")
                        for r in rs], color=INK)
    ax.set_ylim(-0.75, len(rs) - 0.35)
    ax.tick_params(axis="x", colors=INK2)
    ax.set_xlim(-lim, lim)
    from matplotlib.ticker import FuncFormatter
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: signed(v, "{:+.2f}")))
    ax.set_xlabel("camera elevation - survey (m): dot = median, bar = ±NMAD; positive = camera too high", color=INK2)
    ax.grid(axis="x", color=GRID, lw=0.6)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_title("Camera DEMs against the surveys, with how independent each check is and why", loc="left",
                 fontsize=10.5, color=INK)
    fig.tight_layout()
    fig.subplots_adjust(right=0.64)
    if disabled or not_yet:
        dd = sorted({(r["date"], r["comparison"].split("(", 1)[-1].rstrip(")")) for r in disabled})
        txt = ("Not built, at the user's request: " + "; ".join(
            f"{d} ({textwrap.shorten(w, 60, placeholder='...')})" for d, w in dd)) if dd else ""
        ny = sorted({r["date"] for r in not_yet})
        if ny:
            txt = (txt + ". " if txt else "") + "Enabled, no comparison under this root yet: " + ", ".join(ny)
        # under the legend, never over it: placed from the legend's drawn extent (the saved figure
        # grows to hold it, bbox_inches='tight')
        fig.canvas.draw()
        bb = leg.get_window_extent().transformed(fig.transFigure.inverted())
        fig.text(bb.x0, bb.y0 - 0.015, "\n".join(textwrap.wrap(txt, 150)), fontsize=8, color=INK2, va="top",
                 ha="left")
    fig.savefig(path, dpi=130, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    what = ap.add_mutually_exclusive_group(required=True)
    what.add_argument("--date", help="survey date to build, YYYY-MM-DD (a date in --config)")
    what.add_argument("--all", action="store_true", help="every enabled date")
    what.add_argument("--summary", action="store_true", help="table and figure across the dates built")
    ap.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT),
                    help=f"products go in <root>/<date>/ (default {DEFAULT_OUTPUT_ROOT})")
    ap.add_argument("--config", default=str(DEFAULT_CONFIG), help="dates and cameras (surveys/survey_dates.csv)")
    ap.add_argument("--surveys", default=str(DEFAULT_SURVEYS), help="surveys per date (surveys/surveys.csv)")
    ap.add_argument("--steps", default=",".join(STEPS),
                    help=f"comma-separated subset of {','.join(STEPS)} (default all, in that order)")
    ap.add_argument("--force", action="store_true", help="rebuild the steps asked for even if up to date")
    ap.add_argument("--force-disabled", action="store_true", help="build a date (or camera) disabled in --config")
    ap.add_argument("--dry-run", action="store_true",
                    help="print what would be done: photos per day, files found, labels, run time")
    ap.add_argument("--photo-roots", nargs="+", default=DEFAULT_PHOTO_ROOTS,
                    help="folders searched recursively for the photos (default: the station's)")
    ap.add_argument("--survey-dirs", nargs="+", default=DEFAULT_SURVEY_DIRS,
                    help="folders where survey and envelope files named in the tables are looked for")
    ap.add_argument("--calibration", default=str(HERE / "calibration"), help="calibration folder")
    ap.add_argument("--adcp", help="sig1000_waves_ALL.csv (historical_forcing.py default: the repo's)")
    ap.add_argument("--adcp-navd88", help="adcp_water_level_navd88.csv")
    ap.add_argument("--chatham", help="Chatham 8447435 6-min NAVD88 CSV (time, level)")
    ap.add_argument("--wis", help="WIS ST63064 hourly CSV")
    ap.add_argument("--ndbc", help="NDBC 44013 CSV")
    ap.add_argument("--no-download", action="store_true", help="never download the Chatham record")
    ap.add_argument("--live-contours", default=str(DEFAULT_LIVE_CONTOURS),
                    help="live era: the station's contour file (contour_points_timex.csv)")
    ap.add_argument("--live-cron", default=str(LIVE_CRON), help="waterline_timex_cron.sh (its wave settings)")
    ap.add_argument("--gnssr-spline", default=None,
                    help="live era: GNSS-R spline for a cross-check (default marconi_water_level.py's)")
    ap.add_argument("--gauge-csv", default=None, help="live era: Chatham archive for that cross-check")
    ap.add_argument("--window", nargs=2, metavar=("FIRST", "LAST"), default=None,
                    help="photo days to use instead of the table's (all cameras)")
    ap.add_argument("--utc-hours", default=None, help="instead of the table's, e.g. 13.5-18 or all")
    ap.add_argument("--setup-coef", type=float, default=SETUP_COEF,
                    help=f"wave-setup coefficient C (default {SETUP_COEF}, the live station's; 0 = none)")
    ap.add_argument("--setup-fitted-to", default=None, metavar="SURVEY or none:HOW",
                    help="where a --setup-coef other than the known fits was fitted: the survey file it was "
                         "fitted to (comparisons with that survey become CIRCULAR), or 'none:<how>' for a C fitted "
                         "to no survey (e.g. 'none:repeat crossings, dem_from_contours.py --fit-setup on "
                         "2025-01-18..23'). Required for such a C: a C of unknown origin is refused")
    ap.add_argument("--geotiff-epsg", type=int, default=GEOTIFF_EPSG,
                    help=f"horizontal EPSG code written into the GeoTIFFs (default {GEOTIFF_EPSG}, as asked; the "
                         f"grids are in the frame of the GCPs, assumed NAD83(2011) / UTM 19N, EPSG:6348, as the "
                         f"surveys)")
    ap.add_argument("--keep-detector-debug", action="store_true",
                    help="keep the detector's debug images (~8 MB a photo; default: removed once the "
                         "waterlines are written)")
    ap.add_argument("--keep-moved-days", action="store_true",
                    help="keep days whose pointing differs from the calibration (default: left out)")
    ap.add_argument("--pointing-frames", type=int, default=3, help="frames per day for the horizon (default 3)")
    ap.add_argument("--pointing-log-dir", default=str(HERE / "archive"),
                    help="live era: folder of the station's pointing_<cam>.csv monitor logs (summarised)")
    ap.add_argument("--pointing-fail", type=float, default=POINT_FAIL_DEG,
                    help=f"deg: a day rotated more than this from the reference is a different pointing "
                         f"(default {POINT_FAIL_DEG}, as pointing_check.py)")
    ap.add_argument("--nice", type=int, default=10,
                    help="lower this run's priority by this much (os.nice; inherited by every script it runs) so "
                         "the station's cron keeps the two cores (default 10; 0 = not lowered)")
    ap.add_argument("--synthetic", default=None, metavar="NOTE",
                    help="the inputs are SYNTHETIC test data (e.g. the sandbox fixture): the README, comparison "
                         "pages and summary rows say 'SYNTHETIC FIXTURE: tests the pipeline, not the beach' with this "
                         f"note (also set by a file named {SYNTHETIC_MARKER} in a photo root)")
    ap.add_argument("--dem-cell", type=float, default=DEM_CELL)
    ap.add_argument("--dem-min-points", type=int, default=DEM_MIN_POINTS)
    ap.add_argument("--dem-max-spread", type=float, default=DEM_MAX_SPREAD)
    ap.add_argument("--dem-max-hs", type=float, default=DEM_MAX_HS, help="0 = keep every frame")
    ap.add_argument("--dem-max-day-offset", type=float, default=DEM_MAX_DAY_OFFSET, help="0 = off")
    args = ap.parse_args()

    if args.summary:
        return summary(args)
    # a shutdown or kill (SIGTERM) ends the run like Ctrl-C: the build writes 'status: failed' and why
    # instead of leaving the previous README looking current
    import signal

    def _terminated(signum, frame):
        raise SystemExit(f"terminated by signal {signum} (SIGTERM)")
    try:
        signal.signal(signal.SIGTERM, _terminated)
    except (ValueError, OSError):
        pass
    if args.nice and not args.dry_run:
        try:
            level = os.nice(int(args.nice))
            say("priority", f"niceness {level} (--nice {args.nice}: the station's cron keeps priority)")
        except (OSError, AttributeError) as exc:
            warn(f"--nice {args.nice} not applied: {exc}")
    steps = [s.strip() for s in args.steps.split(",") if s.strip()]
    bad = [s for s in steps if s not in STEPS]
    if bad:
        ap.error(f"--steps: unknown {', '.join(bad)} (expected {','.join(STEPS)})")
    if args.window:
        try:
            a_, b_ = (date_cls.fromisoformat(x) for x in args.window)
        except ValueError:
            ap.error(f"--window {' '.join(args.window)}: expected two dates YYYY-MM-DD")
        if b_ < a_:
            ap.error(f"--window {' '.join(args.window)}: the first day is after the last")
    if args.utc_hours is not None and args.utc_hours.strip().lower() not in ("", "all"):
        try:
            h0, h1 = (float(v) for v in args.utc_hours.split("-"))
        except ValueError:
            ap.error(f"--utc-hours {args.utc_hours!r}: expected e.g. 13.5-18 or all")
        if not (0 <= h0 < h1 <= 24):
            ap.error(f"--utc-hours {args.utc_hours!r}: need 0 <= first < last <= 24 (UTC)")
    if args.setup_coef < 0:
        ap.error("--setup-coef must be >= 0")
    if args.setup_coef and not any(abs(k - args.setup_coef) < 1e-9 for k in SETUP_FITS) \
            and not args.setup_fitted_to and not args.dry_run:
        ap.error(f"--setup-coef {args.setup_coef}: where was this C fitted? Give --setup-fitted-to <survey file> "
                 f"(the survey it was fitted to) or --setup-fitted-to 'none:<how>' (fitted to no survey). A C of "
                 f"unknown origin could be fitted to the very survey it is compared with")
    if args.setup_fitted_to and any(abs(k - args.setup_coef) < 1e-9 for k in SETUP_FITS):
        warn(f"--setup-fitted-to ignored: C = {args.setup_coef} is a known fit (SETUP_FITS)")
    if args.gnssr_spline is None:
        from marconi_water_level import GNSSR_SPLINE, GAUGE_CSV
        args.gnssr_spline = GNSSR_SPLINE
        args.gauge_csv = args.gauge_csv or GAUGE_CSV
    args.photo_roots = [str(p) for p in args.photo_roots]
    args.survey_dirs = [str(p) for p in args.survey_dirs]
    if args.setup_coef == 0:
        args.setup_coef = 0.0
    cfg = read_table(args.config, CONFIG_FIELDS)
    srv = read_table(args.surveys, SURVEY_FIELDS)
    say("configuration", f"{args.config} ({len(cfg)} rows), {args.surveys} ({len(srv)} surveys)")
    say("setup coef", f"C = {args.setup_coef}")
    if args.date:
        res = build_date(args.date, cfg, srv, args, steps)
        return EXIT_CODES.get(res, 1)
    results = {}
    for d in sorted({r["date"] for r in cfg}):
        results[d] = build_date(d, cfg, srv, args, steps)
    rule("all dates")
    for d, r in results.items():
        say(d, r + {"partial": " (PARTIAL: see its README)", "failed": " (FAILED: see its log)"}.get(r, ""))
    if not args.dry_run:
        summary(args)
    return worst_exit(results.values())


# Exit codes: 0 built (or dry run), 1 a step failed, 2 the date is disabled, 3 built but PARTIAL
# (a camera contributed nothing, a part is missing: a camera's pointing check or photo map, or the
# outputs on disk are not one build), 4 not started: another build of the date is running.
EXIT_CODES = {"built": 0, "dry-run": 0, "disabled": 2, "failed": 1, "partial": 3, "busy": 4}
# --all returns the code of its WORST date by severity, not the largest number: a failed date (1)
# outranks one not started because another build held it (4), which outranks a partial one (3).
# Disabled dates (2) are skipped by --all and do not count.
SEVERITY = ("failed", "busy", "partial")


def worst_exit(results):
    """Exit code of an --all run: the most severe per-date result (failed > busy > partial > built)."""
    res = [r for r in results if r != "disabled"]
    unknown = [r for r in res if r not in EXIT_CODES]
    if unknown:
        return 1                                  # an unexpected result is a failure, never a pass
    for r in SEVERITY:
        if r in res:
            return EXIT_CODES[r]
    return 0


if __name__ == "__main__":
    sys.exit(main())

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
(--setup-fitted-to <survey file> or 'none:<how>'), or it is refused (exit 5).

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
            has one). Features matched on 40+ points decide (good to ~0.1
            deg): a turn over 0.5 deg is a DIFFERENT POINTING, no turn
            overrules a horizon jump; a horizon reading counts only on a
            clear sea horizon (not a foam, spray or fog edge). A day measured to
            differ, or whose change one method reported and nothing could
            confirm (UNCONFIRMED), is LEFT OUT of the detection and listed
            with what was measured (--keep-moved-days keeps both); if most
            days differ from the reference nothing is left out and a warning
            says to look. Nothing is dropped silently. OPT-IN: with
            --correct-pointing (an evaluation, in its own root
            <root>_pointing-corrected) a day whose matched features pass
            pointing_daily.py's gates (a turn over MATCH_ERR_DEG, held
            through the day on its first and last photos) is detected and
            its waterlines mapped with its own pointing (the calibration
            turned by the change measured that day) instead; every other
            day as without it.
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
            cannot tell which day is off: every day is kept, and the README
            says so), asc_to_geotiff.py (EPSG:6348, NAD83(2011) / UTM 19N,
            the surveys' frame; NAVD88 EPSG:5703 vertical keys on the DEM),
            and the DEM page by dem_figure.py with the date's own camera
            positions. --dem-loess (off by default) also writes the
            scale-controlled loess surface and its prediction error
            (dem_loess.py; <stem>_loess*.asc/.tif): a second surface beside
            the DEM, which no figure, comparison or label uses.
  maps      daily_elevation_map.py per camera on the window's photos, and a
            plan view of the filtered waterlines coloured by elevation over
            the survey (or the DEM).
  compare   survey_compare.py for every survey of the date in
            surveys/surveys.csv. A survey whose file is missing is skipped
            with a message saying which file to provide and where. Never on
            a BLIND date (below).
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
A usage error (an unknown option or a bad value) exits 6, never 2.
ONE C PER ROOT: a run that would rebuild a date's waterlines with a C other
than the one they carry (it would overwrite that product) is refused, exit
5, and told the root that C belongs in (the default root for the station's
C, <default root>_c<C> for any other; a root whose name gives it to that
very C although it holds another, e.g. a default root whose first build
used --setup-coef 0, is said to hold the other C, and the run is given a
root beside it, <root>_c<C>); --replace overwrites on purpose. A
--setup-coef of unknown origin without --setup-fitted-to is refused, exit 5
(a --dry-run says so). --all returns the code of its worst date by
severity: 1 if any date failed, else 5 if any was refused, else 4 if any
was busy, else 3 if any was partial, else 0.

WHAT A LATER RUN CANNOT CHANGE. The outputs on disk keep the settings they
were built with: the DEM settings, the setup coefficient, the photo window
and hours, and the input files. A later run that does not rebuild the steps
using an option (e.g. --steps compare --window ..., or a --live-contours
whose step fails) says NOT APPLIED, and the README, provenance and the
"rebuild everything" line report what built the outputs. --window and
--utc-hours apply only to a run of every step from the forcing to the DEM;
such a run replaces the date's product built for another window (it says
so: build a recheck with another window in its own --output-root).
--synthetic NOTE (or a SYNTHETIC_FIXTURE file in a photo root) marks a build
from synthetic test inputs on every output.

BLIND DATES (BLIND_DATES below). A date whose survey is HELD by someone
else for an independent check of the DEM (2024-10-23: the GCPs of 23 Oct
2024, held by Dr Chris Sherwood, USGS) is built as a DEM only: forcing,
pointing, detection, filter, DEM, GeoTIFFs, maps and README. Its survey
rows never reach a file look-up: build_plan withholds them (names, types
and dates only), resolve_survey refuses any row of the date or of the held
survey's date (BlindSurvey, before any folder is searched), the compare
step does not run and the plan view draws no survey, whatever the tables
say and whatever lies in the survey folders. As a backstop, a runtime
guard (an audit hook) refuses, in this script's own process, a Python-level
open for reading of a file named like a survey of the held date (its date
and xyz/GCP/target/... in a .csv/.txt/... name; by its name or the name of
the file a link points to), and no step this script starts is handed one
(run_cmd checks every argument; the photo folders a step is given were
scanned first in this process). Readers that open files in C (pandas'
default CSV engine, OpenCV) and the steps' own processes are not seen by
the hook: they are safe because no look-up ever hands them such a path.
The README and provenance say 'blind: the survey is held by ... for an
independent check; nothing in this build read it', and what that check can
test given the calibration used (PARTLY-CIRCULAR for the camera geometry
when the calibration rests on the held survey; INDEPENDENT when it traces
to other GCPs only). Status complete when everything else was built;
--summary lists the date as 'blind (DEM only)'. A folder holding outputs an
earlier build made with a survey of the date is refused (exit 5): build it
in a fresh --output-root. A date is unblinded by removing its BLIND_DATES
entry (a commit), after the holder has reported.

ONE CALIBRATION PER ROOT. --eo c1=PATH c2=PATH --calibration-variant NAME
builds a date with another calibration than the table's eo_file (e.g. one
made without the date's held GCPs), into a root named for it,
<root>_<NAME> (any other root is refused, exit 5); a table build is refused
in a root that holds a variant. The README, the DEM page, provenance and
--summary say what each calibration was solved from (CALIBRATION_ORIGINS,
by name and sha256; a calibration turned from another by matching photos
is traced through the reference(s) its notes name: estimate_eo_rotation.py's
'X.yaml rotated by N deg', the first named giving the camera position, or
gcpfree_eo.py's 'Reference EO: X.yaml ... sha256 H', 'Position (x, y, z)
kept from X.yaml' and 'pan change N deg from X.yaml', each claim checked
against the files: a claim that is not true makes the origin NOT KNOWN,
never independent) and carry its pre-registered acceptance (gcpfree_eo.py's
'CALIBRATION B ...: DELIVERED | DELIVERED-WIDER (<why>)', or a line
'# acceptance: PASS | PASS-WIDER <uncertainty> | FAIL <why>'): FAIL is
refused (exit 5), and so, on a BLIND date, is a file that states none; a
--eo file that is missing is refused, never replaced by the table's. As
with C, a run that would rebuild a date's pointing or waterlines with
another calibration (by content, sha256) than they carry is refused, exit
5; --replace overwrites on purpose. A run that does not rebuild them says
NOT APPLIED and reports the calibration that built the outputs. Every
freshness stamp of the pointing and detection records the calibration's
content (sha256): an edited calibration rebuilds. A variant's DEM files
carry its name (dem/<date>_<NAME>_dem.tif) and every GeoTIFF's description
names the calibration (table or variant, EO files, sha256 prefix, and on a
blind date what the holder's check can test), so files separated from
their folders still say which DEM they are.

CALIBRATION B OF A BLIND DATE (GCPFREE_PINS). A turned calibration is
called INDEPENDENT of the held survey only when the PRE-REGISTERED GCP-free
procedure made it: its 'Made by gcpfree_eo.py (sha256 H)' line names a
pinned gcpfree_eo.py, its report beside it (same run folder) records the
same B, verdict and file sha256 for a delivery run with the pinned modules,
its run.json records the file, and its 'X.yaml rotated by N deg' line is
true of its numbers. Anything else is 'NOT the pre-registered procedure',
PARTLY-CIRCULAR for the holder's check (a turn by hand, by another or an
edited procedure, could have been chosen after seeing calibration A). The
procedure's variant names (gcpfree; gcpfree-aug for the second product of
its one allowed second run) are its own: a B is built as the variant its
report names, all cameras from one run, and no other calibration under
those names (exit 5). A camera the procedure did NOT DELIVER (outcome rule
per-camera) may be left out of --eo: its report in the same run folder must
say so (a delivered camera is never left out, exit 5); the README names it
NOT DELIVERED with that report, and the table's calibration is never used
in its place. The README gives B's two stated uncertainty terms (rotation,
inherited camera position: neither a validated bound on real photos), its
run, outcome rule and whether its file is committed to git.

SHARING. --share --date D [--output-root R] packs one product into
<R>/share/<name>.tgz (+ .sha256): README, provenance, the GeoTIFFs and DEM
page, the maps, the calibration files used and, for B, the procedure's
reports, run record and pre-registration, with a MANIFEST of every file's
sha256. On a blind date the name says which DEM it is:
<D>_A_official_calibration.tgz (the table's calibration) and
<D>_B_gcpfree.tgz (<D>_B2_gcpfree_aug.tgz for a second GCP-free product).
A B that is no longer the pre-registered procedure's is refused; an
archive is never written over (--replace).

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
        --setup-fitted-to 'none:repeat crossings, dem_from_contours.py --fit-setup on 2025-01-18..23' \\
        --output-root /mnt/I2Rgus_Data/survey_products_c0.034
        (another C in its own --output-root: in the default root it would replace the C = 0.037
         product of that date, so it is refused there without --replace; a 2025 date's detection
         runs again, ~50 min a two-camera week)
    python3 survey_products.py --date 2024-10-23           # BLIND: DEM only (calibration A, the table's)
    python3 survey_products.py --date 2024-10-23 \\
        --eo c1=calibration/gcpfree/runs/<run id>/CACO03_c1_gcpfree-oct2024_EO.yaml \\
             c2=calibration/gcpfree/runs/<run id>/CACO03_c2_gcpfree-oct2024_EO.yaml \\
        --calibration-variant gcpfree --output-root /mnt/I2Rgus_Data/survey_products_gcpfree
        (calibration B, made without the held GCPs by gcpfree_eo.py, in its own root; give only the
         cameras it delivered)
    python3 survey_products.py --date 2024-10-23 --share   # -> <root>/share/2024-10-23_A_official_calibration.tgz
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
from setup_model import G, setup_c, wave_scale          # noqa: E402,F401  the setup formula, once (G: as before)

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

# Wave setup, as the live station since 6 Oct 2026 (waterline_timex_cron.sh SETUP_COEF). ONE value
# for the station and these products: historical_forcing.py's SETUP_COEF and the cron's must equal it
# (test_survey_products.py checks; a run warns when the cron's differs), and it must be a key of
# SETUP_FITS (where it was fitted decides which comparisons are CIRCULAR). A new C is changed in all
# three places and given its own SETUP_FITS entry.
SETUP_COEF = 0.037
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
            "how": "fitted 8 Oct 2026 as the median of the per-frame C = -(waterline - RTK)/sqrt(Hs*L0) of the "
                   "waterlines lying on the 2026-09-29 RTK transects (14 frames, 29 Sep - 5 Oct 2026; 90% "
                   "0.028-0.058), each line taken where it lay WITHOUT setup; repeat crossings (no survey): "
                   "dem_from_contours.py --fit-setup on those lines built with C = 0 prints a near-best band (C "
                   "within 5 mm of its best spread) of +0.006 to +0.034 (filtered, --max-hs 1.5; a flat minimum) "
                   "and +0.026 to +0.040 on the unfiltered rows, also with --max-hs 1.5 (the same spread with each "
                   "line re-projected: 0.029, 90% ~0.00-0.04), so on one week it does not pin C by itself. "
                   "Rechecked Oct 2026 with "
                   "each line re-projected at still water + setup: the transects stop at their lowest shot "
                   "(+1.17 to +1.80 m), below which a line can only read low and above which a high-reading line "
                   "drops out, so a fit over every crossing is biased by where the survey stops (0.067 over all "
                   "crossings); over only the lines at least 0.1-0.2 m above the lowest shots the fit gives "
                   "~0.043-0.046 (day-block 90% ~0.02-0.055), and the within-cell slope of repeat crossings (no "
                   "survey, re-projected) 0.043 (0.026-0.053). 0.037 lies inside each 90% interval of these "
                   "(~0.00-0.04, 0.026-0.053, ~0.02-0.055; not the biased all-crossings value) and the unfiltered "
                   "rows' --fit-setup band (the filtered lines' band stops just below it): Stockdon et al. "
                   "(2006) setup 0.35 x beta_f for a ~1:10 beach face (beta_f ~0.106). The RTK, which stops above "
                   "the low-tide beach, cannot pin C more closely",
            "wave_currency": "offshore_hs_m = hs_best of archive/waves_marconi.csv (ADCP wh_4061 "
                             "currency), offshore_tp_s = NDBC 44008 peak period",
            # The 29 Sep - 5 Oct 2026 frames on the RTK transects both ways (Oct 2026: survey_products.py
            # --date 2026-09-29 --window 2026-09-29 2026-10-05 --output-root <root>_recheck on the station's
            # rows -- its own root: in the default one it would replace the 26 Sep - 2 Oct product -- then
            # survey_compare's transect check of each frame without setup (re-projected at still water) and with it
            # re-projected): medians of waterline - RTK over the 18 frames (all lines: compare_rtk.py's
            # set), the median setup applied, and per frame (with - without) / setup ('keep_paired': the
            # share of its setup a line keeps once re-projected). These all-crossings medians are biased LOW
            # by where the survey stops (SURVEY_FLOOR_MARGIN in survey_compare.py). 'above_floor': the
            # setup refit's adversarial check (Oct 2026), over only the crossings whose line lies at least
            # 0.1 / 0.2 m above the transect's lowest shot: waterline - RTK at C = 0.037, frames, and the
            # re-projected C fitted there (day-block 90%); 'c_internal_slope': the within-cell slope of
            # repeat crossings on sqrt(Hs*L0), re-projected, no survey (frame-bootstrap 90%).
            "rtk_check": {"without_setup_m": 0.333, "with_setup_m": 0.127, "setup_applied_m": 0.257,
                          "keep_paired": 0.71, "frames": 18, "window": "2026-09-29 .. 2026-10-05",
                          "above_floor": {"margins_m": (0.1, 0.2), "with_setup_m": (-0.037, -0.028),
                                          "frames": (10, 7), "c_fit": (0.046, 0.043),
                                          "c_fit_dayblock90": (0.02, 0.055)},
                          "c_internal_slope": 0.043, "c_internal_slope90": (0.026, 0.053),
                          "c_all_crossings": 0.067}},
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

# Offshore Hs (m) over which a day counts as a storm for the pointing check (storm_days_of): the live
# cron's DEM_MAX_HS, from which the DEM leaves frames out for their waves (fixed here, not --dem-max-hs,
# whose 0 keeps every frame). A camera can be turned in a storm (c1 turned -1.2 deg over 26-28 Sep 2026,
# offshore Hs 1.95 / 3.71 / 2.11 m at its peak on those days), while foam and spray hide the sea horizon:
# a day UNCONFIRMED on a misread horizon in a storm passes UNCONFIRMED to the days next to it that
# nothing could read (pointing_verdicts), where one misread fog edge on a calm day does not.
STORM_HS_M = DEM_MAX_HS

# Pointing check thresholds (pointing_check.py --warn/--fail).
POINT_WARN_DEG, POINT_FAIL_DEG = 0.1, 0.5

# How far each of the pointing step's two methods is trusted (check_pointing -> pointing_verdicts).
# Measured in Oct 2026 on the only real frames there are: the 17:00 UTC c1/c2 timex of 26, 27 and 28
# Sep 2026 (the storm's tail, sent from the NUC; for c2 26/27 Sep the very frames the NUC's feature
# match used), 7 Oct 2026 16:00 c1/c2, CACO03 c1 21 Jan 2025 15:00; planted rotations of them and
# synthetic weather on them (a scratch investigation, not code here; docs/VALIDATION.md 4.6).
#
# FEATURES. pointing_check.measure calls a frame unmatched below MATCH_INLIERS inliers (its
# min_inliers). At or above it the rotation is good to ~0.1 deg: on real pairs with a known answer,
# 69 results at >= 40 inliers erred by <= 0.063 deg; on 402 planted 0.05-2.5 deg rotations by
# <= 0.053 deg (a frame against itself) and <= 0.101 deg (against another day; 40-59 inliers: <= 0.10
# deg on 45 trials, 29 of them against another day); no matched trial with a true rotation under 0.45 deg
# was called MOVED (0 of 229)
# and every one of 0.55 deg or more was (157 of 157). Below it errors grow fast (20-39 inliers: up
# to 0.18 deg; 10-19: 0.58 deg; under 10: anything), and two scenes that do not match share at most
# 4. So a MATCHED result is trusted both ways: a turn above POINT_FAIL_DEG is a different pointing
# (c1 27 Sep 2026: -1.21 deg of pan on 48 inliers, confirmed by tracked dune landmarks: -1.20 to
# -1.25 deg on 26-28 Sep), and one of POINT_WARN_DEG or less means the camera did not turn, whatever
# a horizon reading that is NOT a clear horizon says (c2 26/27 Sep 2026: 73 and 375 inliers, <= 0.07
# deg, while the horizon read a 2.6-3.3 deg tilt change off a foam band).
# The one known way a confident match is wrong: dirt or drops on the housing window at the same
# pixels in the reference and the new frame match at ZERO turn whatever the camera did
# (pointing_check.runner_up). A near-zero winner with a runner-up consensus of AMBIGUOUS_INLIERS or
# more, turned more than AMBIGUOUS_DEG from it, is 'ambiguous' and decides nothing (without dirt the
# runner-up had <= 35 inliers on 551 matched trials and the real pairs; 64-561 in each of the 17
# dirt failures). The runner-up shows the dirt only while the scene itself still matches: with the
# land fogged or dark the dirt is the only consensus (the real 7 Oct frames with the same drops on
# every day and a 0.6 deg tilt on a fogged day: 68 inliers at 0.00 deg, runner-up 6), and a dirt
# winner can hold 388-1573 inliers, so no inlier count tells it from the scene. A zero turn
# therefore never overrules a CLEAR horizon edge that jumps: the day is UNCONFIRMED. A turned
# winner (POINT_WARN_DEG - POINT_FAIL_DEG) is the scene, not dirt, and measures the day's change from
# the reference directly: its tilt or roll is judged by the step's tilt/roll tolerance RUN_TOL_DEG
# (over it: DIFFERENT POINTING, whatever the horizon reads; 0.1 deg of tilt is ~0.13-0.17 m of DEM),
# its pan by POINT_FAIL_DEG (no horizon sees a pan); a smaller turn is a 'small change', kept. A match
# reads tilt and roll to MATCH_ERR_DEG only (the planted rotations against another day: <= 0.101 deg at
# 40-59 inliers, 29 trials; <= 0.094 at 60-99, 41; <= 0.060 at 100-199, 11 -- every run of that study; its
# first summary, without one of them, read 28 / 35 / 10), so a change of 0.2-0.3 deg
# it reads may be one over RUN_TOL_DEG (a clear horizon scatters about as much from day to day at rest: the
# NUC's clear days of 29 Sep - 2 Oct 2026, their features within 0.02 deg, read c1 +0.90 to +0.96 deg of
# tilt and +0.16 to +0.26 of roll): a CLEAR horizon that reads the change over RUN_TOL_DEG and agrees with
# the features within MATCH_ERR_DEG in tilt and in roll makes the day a DIFFERENT POINTING, whatever the
# inlier count (review 3: 0.28 deg on 50 inliers against a clear horizon's 0.37 was a 'small change',
# kept). A clear horizon further from the features than their own error is not the change they measured,
# and the features decide (review 4: with review 3's cut at 100 inliers and agreement within RUN_TOL_DEG,
# 0.34 confirmed a 0.05 deg tilt on 50 inliers while 0.40 did not, and 0.36 never confirmed 0.29 on 150).
# Dirt also matches when the scene does NOT match at all (fog, spray over the land), and then no
# runner-up shows it: its consensus sits at exactly zero turn, as it moves with the camera, while a
# match of the scene between two days reads the camera's own small turns. pointing_check.turn_sigma
# measures how far the fitted turn is from none in units of its own standard error (from the
# inliers' residual scatter; a chi of 3 degrees of freedom for noise alone: 99.9% below ~4.0).
# Measured in Oct 2026 (review 2; the station-floor venv): the consensus of synthetic drops alone on
# the real storm and clear frames (the scene turned or lost) 0.5-2.7 on 21 matches of 113-483
# inliers; the real 26-28 Sep storm frames of c2 (the camera at rest) against the real 7 Oct 16:00
# frame and a copy of it 4.9-16 on 6 matches of 25-102 inliers (lowest: 28 Sep, -0.06/-0.07 deg of
# roll, which the match constrains least); copies of one frame (no true turn: the sandbox's
# stand-ins) 1.2-5.0. A second measurement with other drops: drops alone 0.4-3.3 (18 consensuses of
# 16-117 inliers), the real c2 storm frames against 7 Oct 6.5-13.4 and against each other (26-28 Sep,
# 64-292 inliers) 15.6-18.5. A zero-turn match below TURN_SIGMA_MIN is therefore no reading
# (blind_match, features_trusted); the reference day's own frames are exempt.
MATCH_INLIERS = 40
MATCH_ERR_DEG = 0.1
# What MATCH_ERR_DEG rests on (review r1): the planted rotations against ANOTHER day are of CACO05 c2 2026 only (29
# trials at 40-59 inliers, 41 at 60-99, 11 at 100-199); the matches of a frame against itself add CACO05 c1 2026
# and CACO03 c1 Jan 2025; no CACO04 (Mar 2025) or CACO03 c2 view was tried. Every match against another day had
# its inliers spread 180 px or more along the minor axis of their scatter (pointing_check.inlier_spread; 180-292 px,
# 292-618 along the major one). A match bunched in a smaller part of the frame reads roll
# (and tilt) less well (a review model: 50 rays in a 300 px patch, roll errors up to 0.17 deg): the spread is
# written for every match (inlier_spread_px) and a deciding match under MATCH_SPREAD_TESTED_PX says so in its why
# and in the README; it does not change the verdict.
MATCH_SPREAD_TESTED_PX = 180
MATCH_GOOD = (f"a match of {MATCH_INLIERS} or more features spread over the frame is good to ~{MATCH_ERR_DEG:g} deg on "
              f"the views tested")
AMBIGUOUS_INLIERS = 40
TURN_SIGMA_MIN = 4.0
AMBIGUOUS_DEG = 0.3
RUNNER_UP_SHOWN_MIN = 10           # a runner-up below this is chance matches: its rotation is not written
#
# HORIZON. observed_horizon takes the strongest horizontal edge within +/-HORIZON_WINDOW_PX of the
# predicted horizon in each column checked, and fit_tilt_roll fits tilt and roll to it. On CACO05 the
# horizon is inside the frame on 9 of the 22 columns only (an 800 px segment), where tilt and roll are
# an offset and a slope of a short line: ANY smooth edge fits, so the fit residual cannot tell a
# misread (single storm frames 1.0-7.5 px, the NUC's 3-frame medians of the storm days 1.6-5.0 px,
# against 0.2-0.4 px on its clear days but 0.3-1.4 px on the fixture's rendered horizon): FIT_OK_PX
# only says the fit failed. What tells a sea horizon from a foam, spray or fog edge is the edge
# itself (horizon_edge; per column, then the median over the columns and over the day's frames):
#   prominence   |d/dy| of the picked edge / the median |d/dy| of its search window: real clear
#                frames 61-85, exact rotations of them (47 true rotations) >= 49.8, real hazy frames
#                19.8-75; the six storm frames 2.2-6.7
#   rival        the strongest other edge >= 20 px away / the picked one: clear 0.06-0.24; storm
#                0.57-0.96
#   sky above    share of the columns where the edge is bright above and dark below (sky over sea):
#                clear >= 0.89; the storm picks of 26-27 Sep 0.00-0.22 (the dark top of a bright foam
#                or fog band)
#   to limit     px from the edge to the nearer end of its search window, where that end is not the
#                frame's edge: clear 174-214; storm 26-27 Sep 15-62 (the window cut a larger gradient
#                off; a horizon this near the limit may also lie beyond it)
# A reading is CLEAR when prominence >= HORIZON_PROMINENCE_MIN, rival <= HORIZON_RIVAL_MAX, sky above
# >= HORIZON_SKY_MIN and to limit > HORIZON_LIMIT_PX. Every storm frame fails at least two of these;
# every real clear frame passes, and so do exact rotations of them by up to 2 deg of tilt or roll
# (not clear: c2 tilted -2 deg, its horizon off the top of the frame, and +3 deg of tilt on c2 and
# the Jan 2025 c1, the horizon at the end of the search; matched features measure such turns). On
# synthetic weather put on the clear frames the tests remove 135 of 206 false jumps; a horizon-like
# fog front with a sharp lower edge can still pass, which is why a clear horizon is checked against
# the matched features (pointing_verdicts: a jump the features contradict is no clear day, never
# confirms another day's, and its own day is decided by its turned features, or UNCONFIRMED at zero
# turn). Each day is judged
# against the REFERENCE's tilt/roll as the horizon reads it (the clear days' horizon less the turn
# their matched features measured), not the window's median, which falls between two pointings when
# a camera moved on half the days. The sandbox fixture's RENDERED horizon sits near
# the limits: Jan c1 and Mar c2 pass (prominence 14.9-22, rival 0.49-0.62), Jan c2 does not
# (prominence 6.5-8.0): a rendering artefact, not a camera. A reading that is not clear never enters
# the window's median, offset or clear-day count, and alone never shows a day moved.
HORIZON_WINDOW_PX = 250
HORIZON_PROMINENCE_MIN = 12.0
HORIZON_RIVAL_MAX = 0.65
HORIZON_SKY_MIN = 2.0 / 3.0
HORIZON_LIMIT_PX = 50
# The verdicts that leave a day out of the detection (--keep-moved-days keeps them): DIFFERENT
# POINTING is MEASURED (matched features, or a clear horizon whose change two or more days share);
# UNCONFIRMED is a change one method reported that nothing could confirm or rule out (a horizon jump
# with no matched features, on one day or on an edge that is not a clear horizon; a zero-turn match
# against a clear horizon jump; an ambiguous match), or a day nothing could read, or only a clear
# horizon next to a measured pan, next to a day left out on firm evidence or in a storm (pointing_verdicts:
# nothing shows the camera back). When the days left out are most of the window,
# the step and the README say so with the count (pointing_majority_note).
LEFT_OUT_VERDICTS = ("DIFFERENT POINTING", "UNCONFIRMED")
# OPT-IN per-day pointing CORRECTION (--correct-pointing; pointing_daily.py, whose docstring has the method and the
# gates). Instead of leaving a turned day out, a day whose matched features pass pointing_daily's gates (MATCH_INLIERS
# or more, TURN_SIGMA_MIN or more, inliers spread MATCH_SPREAD_TESTED_PX or more, a turn of more than MATCH_ERR_DEG
# (a smaller one is the match's own error: the day is at rest and gridded as by default) and at most
# pointing_daily.CORRECT_MAX_DEG, the same turn within MATCH_ERR_DEG on the day's first and last photos (the pointing
# held through the day: check_pointing matches them too), not UNCONFIRMED, not the reference day) is detected like the
# others and its waterlines are mapped with its own pointing: the calibration turned by the change the features
# measured that day, the coastal-imaging practice (Rodriguez-Padilla et al. 2020; Arriaga et al. 2022: ~40% less
# shoreline error at Sisal). Every other day is as in the default build (a day measured to differ that fails a gate
# is still left out).
# OFF by default, so the products stay as they are; a corrected build goes in its own root, <root>_pointing-corrected
# (pointing_root_problem), so the two can be compared date by date against the same survey (docs/VALIDATION.md
# 4.6b). No survey enters the correction, so the labels do not change. Not on a BLIND date (its A and B builds are
# fixed by its protocol) and not with --eo (a GCP-free B is the pre-registered procedure's product as made).
POINTING_CORRECTED = "pointing-corrected"

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
    "forcing": ["historical_forcing.py", "setup_model.py"],
    "pointing": ["horizon_check.py", "pointing_check.py", "estimate_eo_rotation.py"],
    "detect": ["detect_original_view.py", "waterline_detector_v5.py", "extract_elevation_contours.py",
               "georectify.py", "view_reproject.py", "compare_dem_survey.py", "timex_timing.py",
               "setup_model.py"],
    "detect_live": ["georectify.py", "extract_elevation_contours.py", "timex_timing.py", "setup_model.py"],
    "merge": [],
    "filter": ["waterline_consistency.py", "waterline_track.py"],
    "dem": ["dem_from_contours.py", "waterline_track.py", "asc_to_geotiff.py", "dem_figure.py", "dem_map.py",
            "setup_model.py"],
    "map": ["daily_elevation_map.py"],
    "map_plan": ["compare_dem_survey.py", "survey_compare.py", "setup_model.py"],
    "compare": ["survey_compare.py", "compare_dem_survey.py", "compare_rtk.py", "asc_to_geotiff.py", "setup_model.py"],
}

# Horizontal CRS written into the GeoTIFFs (--geotiff-epsg). The grids are in
# the frame of the calibration's GCPs, assumed to be that of the surveys: the 2025
# lidar DSMs are tagged EPSG:6348 and the 2026 RTK file says NAD83(2011) / UTM 19N
# (the GCP files state no CRS). Until Oct 2026 the default was 32619 (WGS 84 / UTM
# 19N), a NOMINAL tag taken from an old .prj, which a GIS places ~1-1.5 m off the
# surveys; --geotiff-epsg 32619 still writes it, and the README then says it is
# nominal.
GEOTIFF_EPSG = 6348

# BLIND DATES: a date whose survey is HELD by someone else for an independent check of our DEM.
# The build makes the DEM, its GeoTIFFs, maps and README only. Its survey rows never reach a file
# look-up (build_plan withholds them; resolve_survey raises BlindSurvey for any row of the date or of
# the held survey's date before any folder is searched), whatever surveys.csv (or --surveys) says and
# whatever lies in the survey folders; a runtime guard (arm_blind_guard: an audit hook) refuses, in this
# script's own process, a Python-level open for reading of a file named like a survey of the held date
# (not pandas' C reader, OpenCV or child processes: _blind_audit), and run_cmd hands no step such a
# file. A date is unblinded by removing its entry here -- a commit, so when and why stays in the
# history -- after the holder has reported.
BLIND_DATES = {
    "2024-10-23": {
        "held_by": "Dr Chris Sherwood (USGS)",
        "survey": "the GCPs surveyed 23 Oct 2024 (surveys.csv oct_gcps)",
        "survey_date": "2024-10-23",
        "since": "2026-10-13",
        "why": "the owner's request (13 Oct 2026): the 23 Oct 2024 DEM only, so that Chris can verify its accuracy "
               "independently; the GCP file is kept secret, so no survey of that date can have been used to make "
               "or to compare it"},
}
# The runtime guard refuses to open for reading a file whose name holds a held survey's date (2024-10-23
# or 20241023), one of these words and one of these extensions: the names resolve_survey would look for
# (<date>*xyz*.csv, <date>*GCP*.csv, *<compact date>*xyz*.csv) and the like. Never a calibration (.yaml),
# a product of this script (dem/<date>_info.json, .asc, .tif, .png) or a photo.
BLIND_GUARD_WORDS = ("xyz", "gcp", "target", "extrinsic", "control", "checkshot", "rtk", "survey", "points")
BLIND_GUARD_EXT = (".csv", ".txt", ".xyz", ".dat", ".xls", ".xlsx", ".las", ".laz", ".shp", ".gpkg", ".geojson",
                   ".kml", ".zip")

# What each GCP-solved calibration was solved from, for the README and provenance (and, on a BLIND date,
# for what the holder's check can test). Keyed by file name, with its sha256 (a file of that name with
# other content is 'not known'). CACO03 stood on one mount from Aug 2024 to 24 Jan 2025 and was solved
# from ground control three times (calibration_history.csv). A calibration turned from another by matching
# photos (estimate_eo_rotation.py, or a procedure built on it) is traced through the reference(s) its
# notes name (ROTATED_REF), down to these.
CALIBRATION_ORIGINS = {
    "CACO03_c1_20240827_EO.yaml": ("b942f527d17281f57087975ae7e44aa2de8af079d4b0f912cd87c5c38e435eb2", "2024-08-27",
                                   "the CACO03 GCPs of Aug 2024 (the solution of 27 Aug 2024)"),
    "CACO03_c2_20240827_EO.yaml": ("aab35e9d42c88b89adc2268a837e1b68baa3d7c993c1690a55524aae5c1fb0b2", "2024-08-27",
                                   "the CACO03 GCPs of Aug 2024 (the solution of 27 Aug 2024)"),
    "CACO03_c1_20241023_EO.yaml": ("f0f8a262a6f4dfd956158b74262e9a24440f15054a8c6c30b8f66fd746ed5c98", "2024-10-23",
                                   "the GCPs surveyed 23 Oct 2024 (surveys.csv oct_gcps)"),
    "CACO03_c2_20241023_EO.yaml": ("3a128aee6743e5625e48a0c82b1975b36c021a56a9e6bf41fe6bbb3827dffc2a", "2024-10-23",
                                   "the GCPs surveyed 23 Oct 2024 (surveys.csv oct_gcps)"),
    "CACO03_c1_20250123_EO.yaml": ("c31e0f461e845a394ae62c45a1f524b691b25be64e4348a15599a478e37d5c77", "2025-01-23",
                                   "the GCPs surveyed 23 Jan 2025 (not the Jan 2025 lidar)"),
    "CACO03_c2_20250123_EO.yaml": ("53a021843f9d11f5fe7e9fae5f303e5354b55849f28f8e3603715927e46a035e", "2025-01-23",
                                   "the GCPs surveyed 23 Jan 2025 (not the Jan 2025 lidar)"),
}
# A calibration turned from another by matching photos names its reference(s) in a comment line, as
# estimate_eo_rotation.py writes it: '<reference>.yaml rotated by <angle> deg ...'. Several references
# (e.g. Jan 2025 and, once photos of it arrive, Aug 2024) may be named; the FIRST one named is the one
# whose camera position (x, y, z) the file keeps. A chain (Jan -> Dec -> Oct) is followed file by file:
# each intermediate file must lie next to the one naming it (or in the calibration folder).
ROTATED_REF = re.compile(r"(\S+?\.ya?ml) (?:rotated|turned) by (-?[0-9.]+) deg")
# gcpfree_eo.py (the pre-registered GCP-free procedure for 23 Oct 2024) names its reference as
# 'Reference EO: <file>.yaml (...) sha256 <hex>', the position it keeps as 'Position (x, y, z) kept from
# <file>.yaml', its turn as 'pan change <+deg> deg from <file>.yaml' and how it got it as 'Rotation source:
# <how>'. Each claim is checked against the files here: the reference's content (the stated sha256), the
# file's azimuth against the reference's (the stated pan change), its position against the reference's.
REFERENCE_EO = re.compile(r"Reference EO:\s*(\S+?\.ya?ml)\b[^;]*?\bsha256\s+([0-9a-fA-F]{64})")
PAN_FROM = re.compile(r"pan change ([+-]?[0-9.]+) deg from (\S+?\.ya?ml)")
POSITION_FROM = re.compile(r"Position \(x, y, z\) kept from (\S+?\.ya?ml)")
ROTATION_SOURCE = re.compile(r"Rotation source:\s*([^;]+?)\s*(?:;|$)")
# The pre-registered acceptance of such a calibration, one comment line: '# acceptance: PASS',
# '# acceptance: PASS-WIDER <the larger uncertainty, stated>' or '# acceptance: FAIL <why>'. FAIL is never
# built (exit 5); PASS-WIDER is built and its stated uncertainty goes on every output.
ACCEPTANCE_LINE = re.compile(r"^acceptance\b[^:]*:\s*(PASS-WIDER|PASS|FAIL)\b[\s;:,.-]*(.*)$", re.IGNORECASE)
# gcpfree_eo.py states its pre-registered verdict in its first comment line instead: 'CALIBRATION B (GCP-FREE)
# for the 23 Oct 2024 DEM, camera c1: DELIVERED (all checks passed)' | 'DELIVERED-WIDER (<why: the larger
# uncertainty>)' | 'NOT DELIVERED' (it writes no file then): read as PASS, PASS-WIDER and FAIL.
DELIVERY_LINE = re.compile(r"^CALIBRATION\b[^:]*:\s*(NOT DELIVERED|DELIVERED-WIDER|DELIVERED)\b[\s;:,.-]*(.*)$")
DELIVERY_VERDICT = {"DELIVERED": "PASS", "DELIVERED-WIDER": "PASS-WIDER", "NOT DELIVERED": "FAIL"}
# --calibration-variant: a short name, used in the root's name (<root>_<NAME>) and on every output
VARIANT_NAME = re.compile(r"^[a-z0-9][a-z0-9-]{0,31}$")

# The PRE-REGISTERED GCP-free procedure (gcpfree_eo.py) that makes calibration B for a BLIND date, pinned by
# content: the sha256 of gcpfree_eo.py, of the five modules it computes with, and of its pre-registered constants
# as a run records them in effect (constants_sha256), as committed. A turned calibration is called INDEPENDENT of a
# held survey only when procedure_check finds all of this true: it was made by exactly this code (its 'Made by
# gcpfree_eo.py (sha256 H)' line) with these constants; its numbers are those of the report written beside it (same
# angles, same verdict, the report's record of the file's sha256); its 'X.yaml rotated by N deg', 'pan change' and
# 'Position (x, y, z) kept from' lines are true of its numbers and its references' (the references' content the one
# the report recorded); its run folder is THIS repository's calibration/gcpfree/runs/<id> (GCPFREE_DIR), the run
# was made under THIS repository's PREREGISTERED.txt, unchanged and not committed again since the run, from the
# station's photo folders and with the photos its pre-registration pins, and it is the first finished delivery run
# there (or the one allowed second run with Aug photos); every run's record there (a trial, a run cut off) is
# committed; no commit outside this branch's history (another branch, a stash, the reflog) wrote a pre-registration
# or run record this branch does not hold; and the B file, its report, run.json and run_report.txt are committed in
# this repository, unchanged, each touched by one commit. Otherwise it is 'not the pre-registered procedure': never
# INDEPENDENT, never built or shared under the procedure's variant names. On a BLIND date the B DEM must also be the
# pre-registered B BUILD (blind_b_build_problem: built as A was, one root). What no check here can see: a B edited
# BEFORE its first commit with its report and run.json rewritten to agree (an edit that keeps the 'rotated by'
# magnitude, e.g. a small tilt); a copy of the procedure run outside this repository (another clone or a copy of the
# script) before the pre-registration, or a whole re-clone of the repository after a first run; a branch or stash
# removed together with its reflog entries. Only the sha256s sent to the holder beforehand (the procedure's, the
# pre-registration's, each B file's when its run ends) answer those.
# test_survey_products.py checks that the LAST pin is the files as committed: change gcpfree_eo.py or one of its
# modules and a new pin is appended in the same commit, BEFORE any B is made with it. Once the station has run a
# pinned version (pre-registration or run), its pin is never removed or edited (a B made by it stays traceable);
# until then (Oct 2026: no run outside tests) the single pin follows the committed files.
GCPFREE_PINS = [
    {"gcpfree_eo.py": "e97c24f045c56f85284499c0145dddca8761c0f94ec3612d92eea1bebf52ec2a",
     "modules": {
         "estimate_eo_rotation.py": "2fbb27b69db3604498b5cfd63739e4ab5ebdbb061e58e2397bc1323956406e72",
         "georectify.py": "96f3f52e5adc4ca033c444ae2affb7377f7bb7eb66f1384f00ec14a8e1108a25",
         "horizon_check.py": "a7eeac44fe21ab1414b98568f235cb6c6e1e1892cee965f9ee88c7923f88c23f",
         "detect_original_view.py": "c807b2f907b4856cda08b21ddcf46c1c0470bd23b7559f990fc21e831a19036a",
         "view_reproject.py": "aba55fb29a9104103a2884fbf23efb2ed7bdba7d0462773eede8a701e1abca0c"},
     "constants_sha256": "6259cec4a81e4e0d7da96fec2a6ec4430e589b5a5935401ab274f069feee9f96",
     "pinned": "2026-10-20 (two modules: georectify.py and view_reproject.py test 'in front of the camera' by the "
               "point's depth along the optical axis, no longer by the sign of the projective scale, which depends "
               "on the camera's absolute UTM position and lost every point for azimuths of ~95-275 deg; "
               "georectify.undistort_uv returns NaN where its passes do not converge or beyond the lens model's "
               "fold. Their outputs are bit-identical for every IO and EO file in calibration/ (calibration A not "
               "opened), over every pixel of the frame, so nothing the procedure computes changes; gcpfree_eo.py "
               "and its constants are the same files). It replaces the pin of 2026-10-14 (integration review round "
               "1: the same gcpfree_eo.py, georectify.py e9f1d2c18726... and view_reproject.py 2d63650ec558...), "
               "which replaced the pins of review rounds 0-2 (2387f2d3857d... the last). Nothing was run outside "
               "tests under any of them. The pin is replaced, not appended: a second pin with the same "
               "gcpfree_eo.py sha256 would never be looked up (gcpfree_pin_of takes the first)"},
]
# gcpfree_eo.py's ROOTS: the station's photo folders a delivery run must use (a B made from other folders is a test of
# the procedure, never the delivered B)
GCPFREE_STATION_ROOTS = ["/mnt/I2Rgus_Data/Chelsea_pics", "/mnt/I2Rgus_Data/waterline/composite",
                         "/home/argus_user/owg_marconi/images"]
# gcpfree_eo.py's OUT_DIR in this repository: PREREGISTERED.txt and runs/<run id>/. A B from anywhere else (another
# pre-registration, another folder or repository) is never the delivered B.
GCPFREE_DIR = HERE / "calibration" / "gcpfree"
# The --calibration-variant a pinned B is built with, fixed by the procedure (gcpfree_eo.py RUNS AND HISTORY): the
# delivered B, or the second product of the one allowed second run (Aug 2024 photos) when the first delivered a
# camera. Any other name for such a file is refused (exit 5).
GCPFREE_VARIANTS = {"gcpfree": "calibration B, the GCP-free calibration delivered for the blind check",
                    "gcpfree-aug": "a SECOND GCP-free calibration, from the procedure's second run with Aug 2024 photos "
                                   "(the first run's B stays the delivered one)"}
MADE_BY = re.compile(r"Made by gcpfree_eo\.py \(sha256 ([0-9a-fA-F]{64})\)[^;]*;\s*report (\S+?\.json)")
# A pinned B's 'X.yaml rotated by N deg' line must be true of its numbers and X's to this (gcpfree_eo.py writes N to
# 0.001 deg and the angles to 0.0001 deg: a hand edit of either shows; turned calibrations of other procedures are
# checked to 0.05 deg in calibration_origin)
GCPFREE_TURN_TOL_DEG = 0.01
VARIANT_LINE = re.compile(r"--calibration-variant (\S+?)[\s.(]")
# What a blind date's products are called when they are shared (--share): the table's calibration (A) and the
# GCP-free ones; the archive's name says which, so files separated from their folders cannot be mistaken.
SHARE_NAMES = {"table": "A_official_calibration", "gcpfree": "B_gcpfree", "gcpfree-aug": "B2_gcpfree_aug"}

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
    """Runs cmd, streaming its output to the console (indented) and to log_path. -> exit code. A step is never handed
    a file named like a held survey (BLIND_DATES), by its own name or the name of the file a link points to: the
    runtime guard (an audit hook) does not reach into the child processes, so their arguments are checked here
    (BlindSurvey)."""
    for c in cmd:
        if not isinstance(c, (str, Path)) or not BLIND_DATES:
            continue
        # a path (a Path, or a string with no white space) by its name and the name of the file it really is; any
        # other argument (a sentence, e.g. a GeoTIFF description) by its last component only, with no file look-up
        path_like = isinstance(c, Path) or not any(ch.isspace() for ch in str(c))
        b = held_file_name(c) if path_like else held_name(os.path.basename(str(c)))
        if b:
            raise BlindSurvey(f"refused to start {Path(str(cmd[1] if len(cmd) > 1 else cmd[0])).name}: it would be "
                              f"handed {c}, named like a survey of {b['survey_date']}, and {blind_statement(b)} "
                              f"(survey_products.py BLIND_DATES {b['date']})")
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


def eo_same_content(old, new):
    """The calibration (EO) file of a pointing or detection signature whose stamp differs from it only in that file's
    time, while its content is the one the stamp recorded (eo_sha256): a calibration restored byte for byte (git
    checkout) or only touched needs no rebuild. -> its path, or None."""
    osig = old.get("sig") if isinstance(old, dict) else None
    if not isinstance(new, dict) or not isinstance(osig, dict):
        return None
    ne, oe = new.get("eo_file"), osig.get("eo_file")
    if not (isinstance(ne, list) and isinstance(oe, list) and len(ne) == 3 and len(oe) == 3 and ne[:2] == oe[:2]):
        return None
    if not old.get("eo_sha256") or file_sha(ne[0]) != old["eo_sha256"]:
        return None
    if dict(osig, eo_file=None) != dict(new, eo_file=None):
        return None
    return str(ne[0])


def is_fresh(out, key, sig, outputs, inputs=()):
    """(True, '') when step `key` last ran with this signature, its outputs all exist and
    none of its inputs is newer than them; else (False, why). The calibration file counts by its
    content (the stamp's eo_sha256), not its time: restored byte for byte, it rebuilds nothing."""
    st = stamp_path(out, key)
    if not st.exists():
        return False, "not built yet"
    try:
        old = json.loads(st.read_text())
    except (OSError, ValueError):
        return False, "its stamp is unreadable"
    new = normal(sig)
    same_eo = eo_same_content(old, new) if old.get("sig") != new else None
    if old.get("sig") != new and not same_eo:
        return False, "its settings changed"
    missing = [o for o in outputs if not Path(o).exists()]
    if missing:
        return False, f"{Path(missing[0]).name} is missing"
    t_out = min(Path(o).stat().st_mtime for o in outputs) if outputs else 0
    eo_now = (new.get("eo_file") or [None])[0] if isinstance(new, dict) and isinstance(new.get("eo_file"), list) \
        else None
    unchanged_eo = str(eo_now) if eo_now and old.get("eo_sha256") and file_sha(eo_now) == old["eo_sha256"] else None
    newer = [i for i in inputs if i and str(i) != unchanged_eo and Path(i).exists()
             and Path(i).stat().st_mtime > t_out + 1]
    if newer:
        return False, f"{Path(newer[0]).name} is newer than its outputs"
    if same_eo:
        update_stamp(out, key, sig=new)         # the file's time now, the same token: nothing downstream rebuilds
    return True, ""


def eo_content_fresh(out, key, eo_path):
    """Is step `key` (a pointing or detection) built with the calibration content this run uses? Its
    stamp records the EO file's sha256 (eo_sha256). (True, '') when the same; (False, why) when the
    content differs. A stamp written before the content was recorded has the file's path, size and
    mtime in its signature, unchanged if it got here: the content is recorded now, keeping its token."""
    st = read_stamp(out, key)
    if not st:
        return False, "not built yet"
    now = file_sha(eo_path) if eo_path else None
    if "eo_sha256" not in st:
        update_stamp(out, key, eo_sha256=now)
        return True, ""
    if st["eo_sha256"] != now:
        return False, "its calibration's content (sha256) changed"
    return True, ""


FIRST_BUILT = "first_built.json"


def first_built(out, key=None):
    """When each step of a product FIRST finished ({key: time}, or one key's time): logs/first_built.json, kept across
    rebuilds (a stamp is dropped and rewritten when its step reruns, so its finished_utc is the LAST build). A blind
    date's ORDER states when A's DEM first existed from this, not from its last rebuild."""
    d = load_json(Path(out) / "logs" / FIRST_BUILT) or {}
    return d if key is None else d.get(key)


def note_first_built(out, key, when):
    """Records `when` as step `key`'s first finish unless an earlier one is recorded."""
    if not when:
        return
    f = Path(out) / "logs" / FIRST_BUILT
    d = load_json(f) or {}
    if key in d and str(d[key]) <= str(when):
        return
    d[key] = when
    try:
        f.parent.mkdir(parents=True, exist_ok=True)
        tmp = Path(str(f) + ".tmp")
        tmp.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n")
        os.replace(str(tmp), str(f))
    except OSError:
        pass


def write_stamp(out, key, sig, cmds=(), extra=None, scripts=None):
    """Records a finished step: its signature, commands, the software that ran it
    (commit, sha256 of each script) and a token later steps put in their own signature. The step's first
    finish is kept apart (first_built), whatever rebuilds follow."""
    st = stamp_path(out, key)
    st.parent.mkdir(parents=True, exist_ok=True)
    d = {"sig": normal(sig), "commands": [cmd_text(c) if not isinstance(c, str) else c for c in cmds],
         "finished_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")}
    note_first_built(out, key, d["finished_utc"])
    d["first_finished_utc"] = first_built(out, key) or d["finished_utc"]
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
    current against outputs it may have half rewritten. When the step first finished is kept
    (first_built; a stamp written before that record existed gives its own time)."""
    old = read_stamp(out, key)
    if old:
        note_first_built(out, key, old.get("first_finished_utc") or old.get("finished_utc"))
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


_TABLE_WARNED = set()


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
            # csv puts the fields past the header under None; read as they are, the text after a comma
            # that is not inside double quotes would be dropped without a word. Empty ones (a row ending
            # in a comma, as a spreadsheet may save it) drop nothing and are not warned about
            extra = [x for x in (r.get(None) or []) if (x or "").strip()]
            if extra and (str(path), rd.line_num) not in _TABLE_WARNED:
                _TABLE_WARNED.add((str(path), rd.line_num))
                warn(f"{path} line {rd.line_num}: {len(r[None])} field(s) more than its {len(rd.fieldnames)} "
                     f"columns (a comma in a text that is not in double quotes?): the text after it is dropped")
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


def resolve_file(text, dirs, typed=False):
    """A file the tables name -> list of Paths (empty when nothing matches); globs allowed.
    A bare file name (no folder in it) is looked for by that name in `dirs` only, never in the
    folder the command is run from: what a build finds must not depend on where it was started
    (a stray copy there would win over the station's file). A path with a folder (absolute, or
    relative such as ./x.csv) is taken as given when it exists, else by its name in `dirs`.
    typed=True: a file named on the command line, looked for as given (the current folder) first."""
    if not text:
        return []
    p = Path(text)
    if any(ch in text for ch in "*?["):
        hits = sorted(Path(p.parent).glob(p.name)) if p.is_absolute() else []
        for d in dirs:
            if not hits:
                hits = sorted(Path(d).glob(text))
        return hits
    bare = "/" not in text and os.sep not in text
    if (typed or not bare) and p.exists():
        return [p]
    for d in dirs:
        q = Path(d) / p.name
        if q.exists():
            return [q]
    return []


class BlindSurvey(BaseException):
    """Raised BEFORE any folder is searched when a survey of a BLIND date (BLIND_DATES) is about to be
    looked for. A BaseException on purpose: no 'except Exception' on the way may swallow it and carry
    on; a build that reaches it stops, its README saying failed and why."""


def blind_of(day):
    """The BLIND_DATES entry (with its 'date') for a product date or a held survey's date, else None."""
    if not day:
        return None
    for d, b in BLIND_DATES.items():
        if day in (d, b.get("survey_date")):
            return dict(b, date=d)
    return None


def held_survey(s):
    """The BLIND_DATES entry when survey row s belongs to a blind date or is a survey of a held
    survey's date (another row naming the same survey), else None."""
    return blind_of(s.get("date")) or blind_of(s.get("survey_date"))


def blind_statement(b):
    """The sentence every output of a blind product carries."""
    return (f"blind: the survey is held by {b['held_by']} for an independent check; nothing in this build "
            f"read it")


def held_name(n):
    """The BLIND_DATES entry when the file NAME n is like a survey of a held date (its date, either form, with
    one of BLIND_GUARD_WORDS, and one of BLIND_GUARD_EXT), else None."""
    n = n.lower()
    if not n.endswith(BLIND_GUARD_EXT):
        return None
    words = n.replace("gcpfree", "")          # the GCP-free procedure's own reports name no survey
    if not any(w in words for w in BLIND_GUARD_WORDS):
        return None
    for d, b in BLIND_DATES.items():
        sd = b.get("survey_date") or d
        if sd in n or sd.replace("-", "") in n:
            return dict(b, date=d)
    return None


def held_file_name(path):
    """The BLIND_DATES entry when `path` is named like a survey of a held date, by its own name or by the name
    of the file it really is (a symbolic link, os.path.realpath), else None."""
    try:
        p = os.fsdecode(path)
    except TypeError:
        return None
    hit = held_name(os.path.basename(p))
    if hit or not BLIND_DATES:
        return hit
    try:
        real = os.path.basename(os.path.realpath(p))
    except (OSError, ValueError):
        return None
    return held_name(real) if real != os.path.basename(p) else None


def _blind_audit(event, args):
    """Audit hook (sys.addaudithook): an open for READING of a file named like a held survey (by its name or the
    name of the file a link points to) raises BlindSurvey, when the open goes through Python's own open/io/os IN
    THIS PROCESS: this script's reads, numpy's, PIL's and other Python-level readers'. Writing is not refused.
    NOT covered: readers that open files in C (pandas' default C CSV engine, OpenCV) and the scripts this one
    starts (child processes have no hook; run_cmd checks the files it hands them by name instead). Those are safe
    because no look-up of this script ever returns a held survey's path (resolve_survey refuses before any folder
    is searched); this hook is the backstop for this process's Python reads, not a guarantee for every reader."""
    if event != "open" or not BLIND_DATES or not args:
        return
    path = args[0]
    if path is None or isinstance(path, int):
        return
    mode = args[1] if len(args) > 1 else None
    flags = args[2] if len(args) > 2 else None
    if mode is not None:
        reading = "r" in str(mode) or "+" in str(mode)
    else:
        reading = not isinstance(flags, int) or (flags & 3) != os.O_WRONLY
    b = held_file_name(path) if reading else None
    if b:
        raise BlindSurvey(f"refused to open {os.fsdecode(path)}: it is named like a survey of {b['survey_date']}, "
                          f"and {blind_statement(b)} (survey_products.py BLIND_DATES {b['date']})")


_BLIND_GUARD = {"armed": False}


def arm_blind_guard():
    """Installs _blind_audit once per process (Python 3.8+; the NUC runs 3.8). The hook reads BLIND_DATES
    at each open, so it guards exactly the dates blind at that moment. -> True when armed."""
    if not _BLIND_GUARD["armed"] and hasattr(sys, "addaudithook"):
        sys.addaudithook(_blind_audit)
        _BLIND_GUARD["armed"] = True
    return _BLIND_GUARD["armed"]


def build_plan(date, cfg_rows, survey_rows, args):
    """Everything the steps need for one date, from the tables and the options. On a BLIND date
    (BLIND_DATES) the date's survey rows are withheld from the plan (plan['held']: their names, types,
    dates and table labels only) and plan['surveys'] is empty: no step can look for them."""
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
        # --eo CAM=PATH: this build's calibration instead of the table's
        eo_given = (getattr(args, "eo", None) or {}).get(cam)
        eo_name = Path(eo_given).name if eo_given else r["eo_file"]
        table_eo = resolve_file(r["eo_file"], [args.calibration])
        # an --eo file is taken as given, never looked up by its name elsewhere: a wrong path is refused
        # (eo_option_problems), never replaced by a file of the same name or by the table's calibration
        eo = ([Path(eo_given)] if Path(eo_given).is_file() else []) if eo_given else table_eo
        env = resolve_file(r["envelope_survey"], args.survey_dirs) if r["envelope_survey"] else []
        hours = parse_hours(hours_text)
        th = parse_hours(r["utc_hours"])
        cams[cam] = {"station": r["station"], "first": first, "last": last,
                     "hours": hours, "hours_text": f"{hours[0]:g}-{hours[1]:g}",
                     # the table's own window and hours (what no --window / --utc-hours means)
                     "table_first": r["first_day"], "table_last": r["last_day"],
                     "table_hours_text": f"{th[0]:g}-{th[1]:g}",
                     "eo_file": eo_name, "eo_path": eo[0] if eo else None,
                     "table_eo_file": r["eo_file"], "table_eo_path": table_eo[0] if table_eo else None,
                     "eo_given": eo_given or None,
                     "variant": getattr(args, "calibration_variant", None) if eo_given else None,
                     "io_path": Path(args.calibration) / f"CACO05_{cam}_20240801_IO.yaml",
                     "envelope_survey": r["envelope_survey"], "envelope_path": env[0] if env else None,
                     "notes": r["notes"], "enabled": r in enabled}
    # --eo for some cameras only: a camera the pre-registered GCP-free procedure did NOT DELIVER (its report beside
    # the given B files) is left out of this build, said on every output; any other omission is refused
    # (eo_option_problems), never built with the table's calibration
    not_delivered, omitted_why = ({}, {})
    if getattr(args, "eo", None):
        not_delivered, omitted_why = not_delivered_cameras(list(cams), args.eo, eo_dirs(args))
        for cam in not_delivered:
            cams.pop(cam)
    surveys = [s for s in survey_rows if s["date"] == date]
    blind = blind_of(date) or next((held_survey(s) for s in surveys if held_survey(s)), None)
    held = surveys if blind else []
    surveys = [] if blind else surveys
    plan = {"date": date, "era": era, "cams": cams,
            "first": min(c["first"] for c in cams.values()),
            "last": max(c["last"] for c in cams.values()),
            "surveys": surveys, "out": Path(args.output_root) / date,
            "disabled_cams": {r["camera"]: r["notes"] for r in rows if r not in use},
            "not_delivered": not_delivered, "omitted_why": omitted_why,
            "forced_disabled": bool(args.force_disabled and len(enabled) < len(rows))}
    if blind:
        plan["blind"] = blind
        plan["held"] = [{k: s.get(k) for k in ("name", "survey_type", "survey_date", "label")} for s in held]
        guard = arm_blind_guard()
        say("BLIND", f"{date}: DEM only. {blind_statement(blind)} ({blind['survey']}); no survey of this date is "
                     f"looked for, opened or compared (survey_products.py BLIND_DATES"
                     + ("; a Python-level open of a file named like it is refused in this process, and no step is "
                        "started with one)" if guard else ")"))
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
    """Wave setup C*sqrt(Hs*L0), L0 = g*Tp^2/(2*pi): extract_elevation_contours.py --setup-coef
    (setup_model.setup_c; hs, tp may be the text of a contour row)."""
    return setup_c(float(hs), float(tp), coef)


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


def storm_days_of(plan, days):
    """{day: why} of the days (of `days`, 'YYYY-MM-DD') that count as a storm for the pointing check
    (pointing_verdicts): inside a known storm (STORMS), or with an offshore Hs over STORM_HS_M in the
    forcing's waves (forcing/waves.csv, which the forcing step writes before the pointing step; without
    it, the STORMS days only). Days by UTC, as the photos are named."""
    days = sorted(set(days))
    out = {}
    for d in days:
        for sev in STORMS:
            if sev["first"] <= d <= sev["last"]:
                out[d] = f"inside {sev['name']}"
    peak = {}
    try:
        with open(Path(plan["out"]) / "forcing" / "waves.csv", newline="") as f:
            for r in csv.DictReader(f):
                try:
                    e, h = float(r["epoch"]), float(r["wvht_m"])
                except (KeyError, TypeError, ValueError):
                    continue
                if h == h:
                    d = datetime.fromtimestamp(e, tz=timezone.utc).strftime("%Y-%m-%d")
                    peak[d] = max(h, peak.get(d, h))
    except OSError:
        pass
    for d in days:
        if peak.get(d, 0.0) > STORM_HS_M:
            t = f"offshore Hs up to {peak[d]:.2f} m, over {STORM_HS_M:g} m (forcing/waves.csv)"
            out[d] = f"{out[d]}; {t}" if d in out else t
    return out


# check_pointing's reference for a calibration the GCP-free procedure turned to the window's own photos
# (turned_to_window): it has no solve day, so the reference is the day nearest the survey, and that is no lack
TURNED_REF_KIND = "the day nearest the survey (this calibration was turned to the window's own photos"


def turned_to_window(c, args):
    """Was this camera's calibration made by the GCP-free procedure (gcpfree_eo.py's format: its comment lines say
    so; procedure_check), which turns a GCP solution to the DEM window's own photos pooled? Such a calibration has
    no solve day of its own (its file name holds none)."""
    try:
        o = calibration_origin(c.get("eo_path") or c.get("eo_file"), eo_dirs(args))
    except Exception:                          # what it was made from cannot be read: not claimed
        return False
    return bool((o.get("procedure") or {}).get("format"))


def check_pointing(plan, cam, photos, args):
    """
    Per day: the sea horizon against the calibration (tilt and roll; horizon_check.py's fit) with
    how clearly its edge is a sea horizon (horizon_edge), and the rotation from a reference frame
    (pan too; pointing_check.py's matching, with its runner-up consensus). The verdicts:
    pointing_verdicts. -> (rows, info).
    """
    import cv2
    from georectify import load_intrinsics, load_extrinsics
    from estimate_eo_rotation import horizon_rows, observed_horizon, land_mask
    from horizon_check import fit_tilt_roll, MOVED_PX, FIT_OK_PX
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
        if cal_day is None and turned_to_window(c, args):
            ref_kind = (f"{TURNED_REF_KIND}, pooled: it has no solve day): changes are against that day, the "
                        f"date's own")
        else:
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
        offs, fits, edges = [], [], []
        for p in cands:
            g = cv2.imread(str(p), cv2.IMREAD_GRAYSCALE)
            if g is None:
                continue
            obs = observed_horizon(g, cols, pred, window=HORIZON_WINDOW_PX)
            if np.isfinite(obs).sum() < 3:
                continue
            (dt, dr), res = fit_tilt_roll(io, eo, cols, obs)
            offs.append(float(np.nanmedian(np.abs(obs - pred))))
            fits.append((float(dt), float(dr), float(np.median(res))))
            edges.append(horizon_edge(g, cols, pred, obs, window=HORIZON_WINDOW_PX))
        r = {"camera": cam, "date": day, "photos": len(days[day]), "horizon_frames": len(offs),
             "horizon_offset_px": round(float(np.median(offs)), 1) if offs else "",
             "horizon_dtilt_deg": round(float(np.median([f[0] for f in fits])), 3) if fits else "",
             "horizon_droll_deg": round(float(np.median([f[1] for f in fits])), 3) if fits else "",
             "horizon_fit_px": round(float(np.median([f[2] for f in fits])), 1) if fits else ""}
        r["horizon"] = ("no horizon" if not fits else
                        "unclear" if r["horizon_fit_px"] > FIT_OK_PX else
                        "off the calibration" if r["horizon_offset_px"] > MOVED_PX else "matches")
        r.update(day_edge(edges))
        probe = next((p for p in cands if p.name != ref_src.name), None)
        if probe is None and len(days[day]) > 1:
            probe = days[day][1]
        m = None
        if probe is not None:
            img = cv2.imread(str(probe))
            m = pc.measure(img, [ref], io, min_inliers=MATCH_INLIERS, second=True) if img is not None else None
        r.update(feature_record(m, probe, args.pointing_fail))
        if getattr(args, "correct_pointing", False):
            r.update(pointing_day_check(days[day], probe, ref_src, ref, io, r, ref_day, args))
        rows.append(r)
    storms = storm_days_of(plan, days)
    info = {"reference_frame": ref_src.name, "reference_day": ref_day, "reference_kind": ref_kind,
            "calibration": c["eo_file"], "storm_days": storms}
    info.update(pointing_verdicts(rows, fail=args.pointing_fail, ref_day=ref_day, storm_days=storms))
    # where the calibration puts the horizon: at the top edge of the frame (or off it) the lens model
    # is extrapolated far beyond its GCPs, and a constant offset there says little about the pointing
    fin = np.isfinite(pred)
    if fin.any():
        inside = fin.copy()                # compared only where finite: no NaN warning on numpy 1.17
        inside[fin] = (pred[fin] >= 0) & (pred[fin] < float(io[1]))
        info.update(horizon_pred_rows=[round(float(np.min(pred[fin]))), round(float(np.max(pred[fin])))],
                    horizon_cols_in_frame=[int(inside.sum()), int(len(cols))], frame_rows=int(io[1]))
    off = pointing_offset(info, c["eo_file"])
    if off:
        info["calibration_offset_deg"], info["calibration_offset"] = [off[0], off[1]], off[2]
    return rows, info


def horizon_edge(gray, cols, near_rows, obs, window=HORIZON_WINDOW_PX):
    """
    How clearly the edge estimate_eo_rotation.observed_horizon picked in one frame is a sea horizon,
    from the same gradient (blur 2.5, Sobel 5 in y, a 12-column strip), signed. Per column with an
    edge:
      prominence  |d/dy| at the edge / the median |d/dy| of its search window
      rival       the strongest |d/dy| at least 20 px away / the edge's
      sky above   whether the edge is bright above and dark below (d/dy < 0), as sky over sea
      to limit    px from the edge to the nearer end of its search window, where that end is not
                  the frame's edge (no such end: none)
    -> {"prominence", "rival" (medians), "sky_above" (share of the columns), "to_limit_px" (median,
    or None)}, or {} when no column has an edge.
    """
    import cv2
    fin = np.isfinite(near_rows) & np.isfinite(obs)
    if not fin.any():
        return {}
    H = gray.shape[0]
    # only the rows the windows use, plus a margin wider than the blur and Sobel kernels (the rows
    # inside the windows then get the same gradient as observed_horizon's whole-frame one)
    top = int(min(H, np.max(near_rows[fin]) + window + 32))
    g = cv2.GaussianBlur(gray[:max(top, 1)].astype(np.float32), (0, 0), 2.5)
    gy = cv2.Sobel(g, cv2.CV_32F, 0, 1, ksize=5)
    prom, rival, sky, lim = [], [], [], []
    for c, r0, o in zip(cols, near_rows, obs):
        if not (np.isfinite(r0) and np.isfinite(o)):
            continue
        a, b = int(max(0, r0 - window)), int(min(H, r0 + window))
        prof = gy[a:b, max(0, int(c) - 6):int(c) + 6].mean(axis=1)
        ap = np.abs(prof)
        k = int(round(float(o))) - a
        if not 0 <= k < ap.size or ap[k] <= 0:
            continue
        rest = ap.copy()
        rest[max(0, k - 20):k + 21] = 0
        prom.append(float(ap[k]) / max(float(np.median(ap)), 1e-6))
        rival.append(float(rest.max()) / float(ap[k]))
        sky.append(bool(prof[k] < 0))
        ends = ([k] if r0 - window > 0 else []) + ([ap.size - 1 - k] if r0 + window < H else [])
        if ends:
            lim.append(min(ends))
    if not prom:
        return {}
    return {"prominence": float(np.median(prom)), "rival": float(np.median(rival)),
            "sky_above": float(np.mean(sky)), "to_limit_px": float(np.median(lim)) if lim else None}


EDGE_KEYS = ("horizon_prominence", "horizon_rival", "horizon_sky_above", "horizon_to_limit_px")


def day_edge(edges):
    """A day's horizon-edge columns: the medians of its frames' horizon_edge(), and 'horizon_edge':
    'clear' or what is wrong with the edge (horizon_doubts). Blank without a fitted frame."""
    e = [x for x in edges if x]
    if not e:
        return dict({k: "" for k in EDGE_KEYS}, horizon_edge="")
    lim = [x["to_limit_px"] for x in e if x.get("to_limit_px") is not None]
    r = {"horizon_prominence": round(float(np.median([x["prominence"] for x in e])), 1),
         "horizon_rival": round(float(np.median([x["rival"] for x in e])), 2),
         "horizon_sky_above": round(float(np.median([x["sky_above"] for x in e])), 2),
         "horizon_to_limit_px": int(round(float(np.median(lim)))) if lim else ""}
    r["horizon_edge"] = "; ".join(horizon_doubts(r)) or "clear"
    return r


def _num(v):
    """A number from a row's column ('' or None, as an older build's row without it: None)."""
    try:
        return None if v in ("", None) else float(v)
    except (TypeError, ValueError):
        return None


def horizon_doubts(r):
    """What makes a day's horizon reading something other than a clear sea horizon (the HORIZON_*
    thresholds above): 'faint', 'rival edge', 'inverted', 'at search limit'. A row without the edge
    numbers (made before Oct 2026, when they were not measured) has none."""
    pr, rv, sk, tl = (_num(r.get(k)) for k in EDGE_KEYS)
    out = []
    if pr is not None and pr < HORIZON_PROMINENCE_MIN:
        out.append("faint")
    if rv is not None and rv > HORIZON_RIVAL_MAX:
        out.append("rival edge")
    if sk is not None and sk < HORIZON_SKY_MIN:
        out.append("inverted")
    if tl is not None and tl <= HORIZON_LIMIT_PX:
        out.append("at search limit")
    return out


def doubts_text(r):
    """horizon_doubts in words, with the day's numbers and the clear horizon's."""
    pr, rv, sk, tl = (_num(r.get(k)) for k in EDGE_KEYS)
    t = []
    for d in horizon_doubts(r):
        if d == "faint":
            t.append(f"faint: {pr:.1f}x the median gradient of its search window, where a clear sea horizon is "
                     f"{HORIZON_PROMINENCE_MIN:g}x or more")
        elif d == "rival edge":
            t.append(f"a rival edge {rv:.2f} as strong in the window, where a clear horizon has {HORIZON_RIVAL_MAX:g} "
                     f"or less")
        elif d == "inverted":
            t.append(f"dark above and bright below on {100 * (1 - sk):.0f}% of the columns, as the top of a bright "
                     f"foam or fog band, not sky over sea")
        elif d == "at search limit":
            t.append(f"{tl:.0f} px from the end of the +/-{HORIZON_WINDOW_PX} px search, where an edge is cut off "
                     f"or a horizon lies beyond")
    return "; ".join(t)


def feature_record(m, probe, fail=POINT_FAIL_DEG, warn=POINT_WARN_DEG):
    """A day's feature columns from pointing_check.measure(..., second=True) on its probe frame.
    features: MOVED (turn > fail), moved (warn - fail), ok (<= warn), ambiguous (a winner at <= warn
    with a runner-up consensus of AMBIGUOUS_INLIERS or more, turned more than AMBIGUOUS_DEG from it:
    window dirt?), unmatched (fewer than MATCH_INLIERS inliers), not measured (no probe frame),
    unreadable."""
    rec = {"feature_frame": probe.name if probe is not None else "", "inliers": "", "d_azimuth_deg": "",
           "d_tilt_deg": "", "d_roll_deg": "", "turn_sigma": "", "runner_up_inliers": "", "runner_up_deg": "",
           "runner_up_d": "", "inlier_spread_px": ""}
    if m is None:
        return dict(rec, features="not measured" if probe is None else "unreadable")
    if m.get("status") == "unmatched":                       # fog, rain, night, few features
        return dict(rec, inliers=m.get("inliers", 0), features="unmatched")
    big = max(abs(m["d_azimuth"]), abs(m["d_tilt"]), abs(m["d_roll"]))
    n2, a2, d2 = int(m.get("runner_up_inliers") or 0), m.get("runner_up_deg"), m.get("runner_up_d")
    # dirt on the window can only add a consensus at ZERO turn: a winner that turned is the camera
    amb = big <= warn and n2 >= AMBIGUOUS_INLIERS and a2 is not None and a2 > AMBIGUOUS_DEG
    # a runner-up of a few inliers is two or three chance matches (scenes that do not match share up to 4):
    # its rotation means nothing and is not written
    shown = n2 >= RUNNER_UP_SHOWN_MIN
    ts = m.get("turn_sigma")
    spr = m.get("inlier_spread_px")
    rec.update(inlier_spread_px="/".join(f"{x:.0f}" for x in spr) if spr else "")
    rec.update(inliers=m["inliers"], d_azimuth_deg=round(m["d_azimuth"], 3), d_tilt_deg=round(m["d_tilt"], 3),
               d_roll_deg=round(m["d_roll"], 3), turn_sigma=round(float(ts), 1) if ts is not None else "",
               runner_up_inliers=n2,
               runner_up_deg=round(float(a2), 3) if a2 is not None and shown else "",
               runner_up_d="/".join(f"{x:+.3f}" for x in d2) if d2 and shown else "",
               features="MOVED" if big > fail else "ambiguous" if amb else "moved" if big > warn else "ok")
    return rec


def _rot(r):
    return f"{float(r['d_azimuth_deg']):+.2f}/{float(r['d_tilt_deg']):+.2f}/{float(r['d_roll_deg']):+.2f}"


def _turn(r):
    return max(abs(float(r["d_azimuth_deg"])), abs(float(r["d_tilt_deg"])), abs(float(r["d_roll_deg"])))


def _inliers(r):
    try:
        return int(float(r.get("inliers")))
    except (TypeError, ValueError):
        return 0


def _sigma(r):
    """A day's turn_sigma (pointing_check.turn_sigma), or None (not measured: an older build's row)."""
    try:
        return float(r.get("turn_sigma"))
    except (TypeError, ValueError):
        return None


def spread_of(r):
    """A row's inlier_spread_px as [minor, major] (px), or None (not measured: an older build's row)."""
    try:
        a, b = (float(x) for x in str(r.get("inlier_spread_px") or "").split("/"))
        return [a, b]
    except ValueError:
        return None


def bunched(r):
    """A match whose inliers spread less than MATCH_SPREAD_TESTED_PX along the minor axis of their scatter: smaller
    than any match MATCH_ERR_DEG rests on."""
    sp_ = spread_of(r)
    return sp_ is not None and sp_[0] < MATCH_SPREAD_TESTED_PX


def bunched_text(r):
    sp_ = spread_of(r) or [0.0, 0.0]
    return (f"its {_inliers(r)} inliers are bunched (spread {sp_[0]:.0f} x {sp_[1]:.0f} px; every match the "
            f"~{MATCH_ERR_DEG:g} deg rests on spread {MATCH_SPREAD_TESTED_PX} px or more): its tilt and roll may be read "
            f"worse than that")


def no_feature_reading(days, rows):
    """Why the features gave no reading on `days` (rows: the pointing rows), in words, grouped: the days whose
    MATCH_INLIERS or more features matched within their own noise (blind_match), those with too few matched
    features, those with no second photo or an unreadable one. Never 'too few matched features' of a day whose
    CSV row shows hundreds."""
    by = {r.get("date"): r for r in rows}
    groups = {}
    for d in sorted(days):
        r = by.get(d) or {}
        k = ("blind" if blind_match(r) and _inliers(r) >= MATCH_INLIERS else
             "none" if r.get("features") == "not measured" else
             "unreadable" if r.get("features") == "unreadable" else "few")
        groups.setdefault(k, []).append(d)
    text = {"blind": (f"the features matched the reference, but within their own noise (turn_sigma under "
                      f"{TURN_SIGMA_MIN:g}: what dirt on the window, or a camera that did not move at all, gives), "
                      f"so the match is no reading"),
            "few": f"too few matched features against the reference frame (under {MATCH_INLIERS})",
            "none": "no second photo of the day to match features on",
            "unreadable": "the photo for the feature match could not be read"}
    if len(groups) == 1:
        k, ds = next(iter(groups.items()))
        return ("on each " if k == "blind" and len(ds) > 1 else "") + text[k]
    return "; ".join(f"on {', '.join(ds)} {text[k]}" for k, ds in groups.items())


def blind_match(r):
    """A zero-turn match ('ok') no further from no turn than its own noise (turn_sigma under
    TURN_SIGMA_MIN): what dirt on the housing window gives whatever the camera did (and a camera that did not
    move at all), so no reading."""
    ts = _sigma(r)
    return r.get("features") == "ok" and ts is not None and ts < TURN_SIGMA_MIN


def features_trusted(r, ref=None):
    """Whether a day's feature match is a reading (pointing_verdicts): matched on MATCH_INLIERS or more,
    not ambiguous and, at zero turn, more than its own noise from none (blind_match; ref, the reference
    day's row, is exempt: two frames of one day match at zero turn by nature)."""
    return (r.get("features") in ("ok", "moved", "MOVED") and _inliers(r) >= MATCH_INLIERS
            and (r is ref or not blind_match(r)))


def horizon_alone(r):
    """An UNCONFIRMED day whose only evidence is a horizon reading that is not a clear horizon (no
    feature reading, not even an ambiguous one): such a reading alone never shows a day moved. Not a day
    taken in for the day next to it (unread_next_to)."""
    return (r.get("verdict") == "UNCONFIRMED" and not features_trusted(r) and r.get("features") != "ambiguous"
            and r.get("horizon") in ("matches", "off the calibration") and bool(horizon_doubts(r))
            and not r.get("next_to"))


def unread_next_to(r):
    """An UNCONFIRMED day taken in for the day next to it (pointing_verdicts): nothing could read it (no
    feature reading, and no horizon fitted or one that is not a clear horizon), or only a clear horizon,
    which cannot see the pan its neighbour turned by; its 'next_to' names the day its run starts at. A row
    of an older build without that column: UNCONFIRMED with no horizon fitted and no feature reading."""
    if "next_to" in r:
        return r.get("verdict") == "UNCONFIRMED" and bool(r.get("next_to"))
    return (r.get("verdict") == "UNCONFIRMED" and not features_trusted(r) and r.get("features") != "ambiguous"
            and r.get("horizon") not in ("matches", "off the calibration"))


def pointing_verdicts(rows, fail=POINT_FAIL_DEG, warn=POINT_WARN_DEG, ref_day=None, storm_days=None):
    """
    Each day's verdict and why from what check_pointing measured: a pure function of the rows'
    columns (completed in place; rows of an older build without the edge numbers count as clear).
    ref_day: the day of the reference frame the features are measured against (None: not known).
    storm_days: {day: why} of the window's storm days (storm_days_of; None: none known).
    -> the window's horizon numbers for check_pointing's info.

    The two methods are put on one footing first. The features measure each day's turn from the
    reference frame; the horizon, each day's tilt/roll against the calibration. The ANCHOR is the
    reference's tilt/roll as the horizon reads it:
      1. from the days with a clear horizon AND a feature reading (features_trusted): each one's
         horizon less the tilt/roll its features measured (a turn of warn or less counts as none);
         the anchor is the median of the group of them that agree within RUN_TOL_DEG: the reference
         day's group when another day confirms its reading (a steady fog front on more days must not
         outvote the reference day's own clear horizon), else the largest (ties: the reference day's,
         then the earlier one; a lone reference reading the others outnumber is the odd one, and the
         step warns that the anchor leaves it out). A day of this kind outside the group is
         CONTRADICTED: its horizon and its features disagree about its change; it is no clear day,
         and a day without features whose jump it reads the same is UNCONFIRMED.
      2. else, the reference day's own clear reading;
      3. else, the largest group of clear readings that agree within RUN_TOL_DEG (nothing ties it to
         the reference: the why says when other clear days read otherwise; when that group is not
         most of the clear days, no day of another group counts as measured to differ);
      4. else (no clear day), the median of the fitted readings; no offset is reported then.
    A JUMP: tilt or roll more than RUN_TOL_DEG from the anchor. The window's offset against the
    calibration and its clear-day count are those of the clear days at the anchor (neither jumping
    nor contradicted). A reading that is not clear never enters them, and alone never shows a day
    moved.

      features (vs the reference)     horizon jump                            no jump
      MOVED (> fail)                  DIFFERENT POINTING                      DIFFERENT POINTING
      ok (<= warn)                    not a clear edge: matches (the          matches
                                      features overrule a misread horizon);
                                      a clear edge: UNCONFIRMED (the two
                                      disagree: a fog front can look like a
                                      horizon, and dirt on the housing window
                                      matches at zero turn when the scene
                                      does not); on the reference day itself:
                                      matches (its reading is the odd one)
      moved (warn - fail)             the features' own change of tilt or roll more than RUN_TOL_DEG:
                                      DIFFERENT POINTING, whatever the horizon reads (a turned match is
                                      the scene, not dirt, and measures the change from the reference
                                      directly); also a clear horizon jump that agrees with the features
                                      within MATCH_ERR_DEG (their own error); else small change (a
                                      horizon jump is overruled)
      ambiguous (window dirt?)        as 'none' (a runner-up turn over fail   small change (UNCONFIRMED
                                      in pan or RUN_TOL_DEG in tilt/roll:     if the runner-up is
                                      UNCONFIRMED)                            turned that much)
      none (unmatched, not            clear edge with the same change on      clear horizon: horizon only;
      measured, unreadable, or ok     2+ days (none of them contradicted):    not clear, or no fit:
      within its own noise:           DIFFERENT POINTING; one day only, or    unchecked
      blind_match)                    not a clear edge: UNCONFIRMED

    A zero-turn match ('ok') whose turn is within its own noise (turn_sigma under TURN_SIGMA_MIN) is
    what dirt on the housing window gives whatever the camera did: it is no reading (features_trusted),
    except on the reference day, whose own frames match at zero turn by nature. A row of an older
    build without turn_sigma counts as a reading.

    A turned match ('moved') reads tilt and roll to ~MATCH_ERR_DEG only: its tilt/roll within
    RUN_TOL_DEG, a CLEAR horizon that reads more than RUN_TOL_DEG and agrees with the features within
    MATCH_ERR_DEG in tilt and in roll makes the day DIFFERENT POINTING too, at any inlier count.

    A day nothing could read ('unchecked': no feature reading, and no horizon fitted or one that is not
    a clear horizon) next to a FIRM day is UNCONFIRMED too: nothing shows the camera back. Firm: a
    DIFFERENT POINTING day, an UNCONFIRMED one with more behind it than a horizon reading that is not a
    clear horizon, or one UNCONFIRMED on such a reading alone (horizon_alone) in a storm: a run of
    consecutive days each left out on its own readings or read by nothing that holds a firm day of the
    first kinds, such misread horizons on two days or more, or a day of storm_days ({day: why}: inside
    STORMS or offshore Hs over STORM_HS_M; storm_days_of). A day UNCONFIRMED on a misread horizon alone
    on a calm day never reaches past itself. A 'horizon only' day (a clear horizon, no feature reading)
    is taken in only when the firm day measured a pan over fail (its horizon cannot see the pan). Runs of
    such days take it from their neighbour, either side, as long as the run starts at a firm day; the
    row's 'next_to' names that day (unread_next_to). Never the reference day itself.

    The window's other pointing ('other_pointing'): when every day measured to differ agrees with the
    others (within RUN_TOL_DEG in tilt and roll, fail in pan) and they outnumber the days measured at
    the reference's pointing together with the days that read nothing else whose features matched the
    reference at zero turn within their own noise ('zero_turn_days': blind_match, MATCH_INLIERS or more),
    the reference itself may be the odd one out: the step and the README say so, and the build is not
    complete (missing_parts), unless the reference day is the calibration's own solve day (its pointing is
    then the calibration's: reference_is_calibration_day). Those days are left out all the same (a product
    of one pointing, the reference's; before Oct 2026 they were all kept, UNCONFIRMED days with them,
    which mixed two pointings).

    The rows are judged in date order whatever order they come in (a tie between two groups goes to the
    earliest day).

    DIFFERENT POINTING and UNCONFIRMED days are left out (LEFT_OUT_VERDICTS; pointing_leave_out says
    which are, as the detection decides it, so no why here says so); every why gives the numbers
    measured.
    """
    from horizon_check import RUN_TOL_DEG as tol
    rows = sorted(rows, key=lambda r: str(r.get("date") or ""))     # the same dicts, completed in place
    fitted = [r for r in rows if r.get("horizon") in ("matches", "off the calibration")]
    clear = [r for r in fitted if not horizon_doubts(r)]
    ref = next((r for r in rows if ref_day and r.get("date") == ref_day), None)

    def trusted(r):
        return features_trusted(r, ref)

    def hz_of(r):
        return float(r["horizon_dtilt_deg"]), float(r["horizon_droll_deg"])

    def feat_of(r):
        """the tilt/roll change the features measured (a turn of warn or less is none)"""
        if trusted(r) and r.get("features") in ("moved", "MOVED"):
            return float(r["d_tilt_deg"]), float(r["d_roll_deg"])
        return 0.0, 0.0

    def est(r):
        """the reference's tilt/roll as this day's horizon and features put it"""
        (a, b), (c, e) = hz_of(r), feat_of(r)
        return a - c, b - e

    def apart(p, q):
        return abs(p[0] - q[0]) > tol or abs(p[1] - q[1]) > tol

    def med(ps):
        return ((float(np.median([p[0] for p in ps])), float(np.median([p[1] for p in ps]))) if ps else (0.0, 0.0))

    def largest_group(rs, val, prefer=None):
        """the group of rs whose values agree within tol of one member's: one holding prefer and another
        member first (prefer's reading confirmed), then the largest (ties: the one holding prefer, then the
        one of the earliest member)"""
        best = None
        for r in rs:
            g = [o for o in rs if not apart(val(o), val(r))]
            has = any(o is prefer for o in g)
            key = (has and len(g) >= 2, len(g), has)
            if best is None or key > best[0]:
                best = (key, g)
        return best[1] if best else []

    tied = [r for r in clear if trusted(r)]
    others = []
    if tied:
        grp = largest_group(tied, est, prefer=ref)
        anchor, how = med([est(r) for r in grp]), "tied"
    elif ref is not None and any(ref is r for r in clear):
        grp = [ref]
        anchor, how = hz_of(ref), "reference"
    elif clear:
        grp = largest_group(clear, hz_of)
        anchor, how = med([hz_of(r) for r in grp]), "clear"
        others = [r for r in clear if not any(r is o for o in grp)]
    else:
        grp = list(fitted)
        anchor, how = med([hz_of(r) for r in fitted]), "fitted"
    at, ar = anchor
    contra = [r for r in tied if apart(est(r), anchor)]
    is_contra = set(id(r) for r in contra)
    base = [r for r in clear if not apart(hz_of(r), anchor) and id(r) not in is_contra]
    is_clear = set(id(r) for r in clear)
    is_fitted = set(id(r) for r in fitted)
    jumped = set(id(r) for r in fitted if apart(hz_of(r), anchor))
    if (how == "tied" and ref is not None and len(grp) == 1 and grp[0] is ref) or how == "reference":
        where = f"where the reference day {ref_day} reads {at:+.2f}/{ar:+.2f}"
        short = f"the reference day {ref_day} ({at:+.2f}/{ar:+.2f} deg)"
    elif how == "tied":
        where = (f"where the reference's pointing reads {at:+.2f}/{ar:+.2f} (the clear days' horizon less the "
                 f"turn their matched features measured)")
        short = f"the reference's pointing ({at:+.2f}/{ar:+.2f} deg)"
    elif how == "clear":
        where = f"where the window's clear days read {at:+.2f}/{ar:+.2f}"
        short = f"the window's clear days ({at:+.2f}/{ar:+.2f} deg)"
    else:
        where = f"where the window's days read {at:+.2f}/{ar:+.2f} (none has a clear horizon)"
        short = f"the window's days ({at:+.2f}/{ar:+.2f} deg)"
    two = (f"; {len(others)} other clear day(s) ({', '.join(o['date'] for o in others)}) read more than {tol} deg "
           f"from these, and nothing ties either group to the reference: the window may hold two pointings"
           if others else "")
    # without a tie to the reference, a group of clear days that is not most of them is not the window's
    # pointing for sure: a day of another group is then never MEASURED to differ
    split = how == "clear" and len(grp) <= len(clear) / 2.0
    overruled = []
    for r in rows:
        f, n = r.get("features"), _inliers(r)
        fit_, jump, clear_ = id(r) in is_fitted, id(r) in jumped, id(r) in is_clear
        disagree = id(r) in is_contra
        h0, h1 = hz_of(r) if fit_ else (0.0, 0.0)
        hz = (f"the horizon read {h0:+.2f}/{h1:+.2f} deg (tilt/roll against the calibration) {where}, more than "
              f"{tol} deg apart" if fit_ else "")
        dz = doubts_text(r) if fit_ else ""
        if not fit_:
            hstate = (f"the horizon could not be fitted (residual {r.get('horizon_fit_px')} px)"
                      if r.get("horizon") == "unclear" else "no horizon was found")
        elif not clear_:
            hstate = (f"{hz}, but that is not a clear horizon ({dz}): a misread (haze, spray or foam), not used"
                      if jump else
                      f"the horizon reading is within {tol} deg of {short} but is not a clear horizon ({dz})")
        elif disagree:
            ft, fr = feat_of(r)
            hstate = (f"the horizon read {h0:+.2f}/{h1:+.2f} deg (tilt/roll against the calibration) {where}: a "
                      f"change of {h0 - at:+.2f}/{h1 - ar:+.2f} deg, where the features measured {ft:+.2f}/{fr:+.2f}, "
                      f"more than {tol} deg apart, on a clear-looking edge: the horizon reading is not the pointing "
                      f"(a fog front can look like a horizon), not used")
        elif jump:
            hstate = f"{hz}, on a clear horizon edge" + (f": the change the features measured, within {tol} deg"
                                                          if trusted(r) else "")
        else:
            hstate = f"the horizon agrees with {short} in tilt and roll"
        kind = (f if trusted(r) else "ambiguous" if f == "ambiguous" else "none")
        rot = _rot(r) if kind in ("ok", "moved", "MOVED", "ambiguous") else ""
        if kind == "none" and f == "ok" and _inliers(r) >= MATCH_INLIERS:
            ftxt = (f"{n} features matched the reference at {_rot(r)} deg (pan/tilt/roll), no further from no turn "
                    f"than their own noise ({_sigma(r):.1f} standard errors, under {TURN_SIGMA_MIN:g}): dirt or drops "
                    f"on the housing window match so, at exactly zero turn, whatever the camera did, when fog or "
                    f"spray hides the scene (a camera that did not move at all reads so too), so the match is no "
                    f"reading")
        elif kind == "none":
            ftxt = {"not measured": "there was no second photo of the day to match features on",
                    "unreadable": "the photo for the feature match could not be read",
                    "ok": f"only {n} features matched the reference ({MATCH_INLIERS} are needed for a reading)"}.get(
                "ok" if f in ("ok", "moved", "MOVED") else f,
                f"the features could not be matched ({n} inliers; {MATCH_INLIERS} are needed for a reading)")
        elif kind == "ambiguous":
            ftxt = (f"the feature match is ambiguous: {n} features at {rot} deg (pan/tilt/roll) and a second "
                    f"consensus of {r.get('runner_up_inliers')} at {r.get('runner_up_d')} deg from the reference "
                    f"(dirt or drops on the housing window match at zero turn whatever the camera did)")
        else:
            ftxt = f"{n} matched features show a turn of {rot} deg (pan/tilt/roll) from the reference"
        if kind == "MOVED":
            r["verdict"] = "DIFFERENT POINTING"
            r["why"] = (f"{ftxt}, more than {fail} deg ({MATCH_GOOD}); {hstate}"
                        + (" (it cannot see a pan)" if fit_ and clear_ and not jump and not disagree else ""))
        elif kind == "ok":
            r["verdict"] = "matches"
            r["why"] = ""
            ts = _sigma(r)
            if jump and not clear_:
                overruled.append(r["date"])
                r["why"] = (f"{hz}, on an edge that is not a clear horizon ({dz}); but {n} matched features show at "
                            f"most {_turn(r):.2f} deg of turn ({rot} deg pan/tilt/roll) from the reference, {warn} deg "
                            f"or less ({MATCH_GOOD})"
                            + (f", and that small turn is {ts:.0f} times its own standard error from none: the scene "
                               f"matched (dirt on the housing window matches at exactly zero turn)"
                               if ts is not None and r is not ref else "")
                            + ": the horizon reading is not the pointing (haze, spray or foam read as the horizon); "
                              "the features decide, day kept")
            elif jump and r is ref:
                # the reference day: the other days' features are measured against its own frame, and
                # its own match (two frames of one day) shows no turn: its horizon reading is the odd one
                overruled.append(r["date"])
                r["why"] = (f"{hz}, on a clear-looking edge; but this is the reference day (the other days' "
                            f"features are measured against its frame) and {n} matched features of its own frames "
                            f"show at most {_turn(r):.2f} deg of turn: its horizon reading is not the pointing (a "
                            f"fog front can look like a horizon); day kept")
            elif jump:
                ru = r.get("runner_up_inliers")
                r["verdict"] = "UNCONFIRMED"
                r["why"] = (f"{hz}, on a clear-looking horizon edge; but {n} matched features show at most "
                            f"{_turn(r):.2f} deg of turn ({rot} deg pan/tilt/roll) from the reference: the two "
                            f"disagree. A fog front can look like a horizon, and dirt or drops on the housing window "
                            f"match at zero turn when the scene itself no longer matches (a second consensus shows "
                            f"that only while the scene still matches"
                            + (f"; here it has {ru} inliers" if ru not in ("", None) else "")
                            + "): the pointing could not be confirmed either way")
        elif kind == "moved":
            # a turned match is the scene (dirt matches at zero turn): its tilt/roll is the day's change from
            # the reference, measured directly and to ~0.1 deg, and is judged by the step's own tilt/roll
            # tolerance, whatever the horizon reads (the pan, which no horizon sees, by pointing_check's fail)
            ft, fr = feat_of(r)
            if max(abs(ft), abs(fr)) > tol:
                r["verdict"] = "DIFFERENT POINTING"
                r["why"] = (f"{ftxt}: a change of tilt or roll of more than {tol} deg ({MATCH_GOOD}, and a turned "
                            f"match is the scene, not dirt on the window); {hstate}")
            elif jump and clear_ and abs(h0 - at - ft) <= MATCH_ERR_DEG and abs(h1 - ar - fr) <= MATCH_ERR_DEG:
                # a match reads tilt/roll to ~MATCH_ERR_DEG: a clear horizon the features confirm within that
                # error tells whether the change is over tol
                r["verdict"] = "DIFFERENT POINTING"
                r["why"] = (f"{ftxt}: tilt and roll within {tol} deg by the features, which read them to "
                            f"~{MATCH_ERR_DEG:g} deg only; the horizon read {h0:+.2f}/{h1:+.2f} deg (tilt/roll against "
                            f"the calibration) {where}, on a clear horizon edge: a change of {h0 - at:+.2f}/"
                            f"{h1 - ar:+.2f} deg, more than {tol} deg, which the features confirm within "
                            f"{MATCH_ERR_DEG:g} deg: a change of tilt or roll of more than {tol} deg")
            else:
                r["verdict"] = "small change"
                r["why"] = (f"{ftxt}: a turn of {warn}-{fail} deg with tilt and roll within {tol} deg (the features "
                            f"measure the change from the reference)" + (f"; {hstate}" if jump or disagree else "")
                            + (f"; that is {h0 - at:+.2f}/{h1 - ar:+.2f} deg from the reference's, more than "
                               f"{MATCH_ERR_DEG:g} deg (the features' own error) from their {ft:+.2f}/{fr:+.2f}: not "
                               f"the change they measured, so the features decide"
                               if jump and clear_ and not disagree else ""))
                if jump:
                    overruled.append(r["date"])
        else:
            rup = ([abs(float(x)) for x in str(r.get("runner_up_d") or "").split("/") if x.strip()]
                   if kind == "ambiguous" else [])
            def alike(o):
                return (o is not r and abs(float(o["horizon_dtilt_deg"]) - h0) <= tol
                        and abs(float(o["horizon_droll_deg"]) - h1) <= tol)
            # the clear days that read the same jump: confirmers unless matched features contradicted them
            same = ([o for o in fitted if id(o) in jumped and id(o) in is_clear and id(o) not in is_contra
                     and alike(o)] if jump and clear_ else [])
            gainsaid = [o["date"] for o in contra if alike(o)] if jump and clear_ else []
            # the runner-up may be the scene's own turn: judged as a turned match is (pan by fail, tilt/roll by tol)
            if rup and (rup[0] > fail or max(rup[1:] or [0.0]) > tol):
                r["verdict"] = "UNCONFIRMED"
                r["why"] = (f"{ftxt}: the camera may have turned by up to {max(rup):.2f} deg (more than {fail} deg of "
                            f"pan or {tol} deg of tilt or roll); {hstate}: the pointing could not be confirmed either "
                            f"way")
            elif same and not split and (not gainsaid or any(trusted(o) for o in same)):
                r["verdict"] = "DIFFERENT POINTING"
                r["why"] = (f"{hz}, on a clear horizon edge, and {', '.join(o['date'] for o in same)} read the same "
                            f"(within {tol} deg): a change of tilt/roll seen on more than one day; {ftxt}, so the pan "
                            f"was not measured")
            elif same or gainsaid:
                r["verdict"] = "UNCONFIRMED"
                r["why"] = (f"{hz}, on a clear-looking horizon edge"
                            + (f", and {', '.join(o['date'] for o in same)} read the same" if same else "")
                            + (f"; but {', '.join(gainsaid)} read the same on a clear-looking edge while matched "
                               f"features there measured no such change (a fog front can look like a horizon)"
                               if gainsaid else
                               f"; but the window's clear days fall in groups more than {tol} deg apart and nothing "
                               f"ties either to the reference")
                            + f"; {ftxt}: the pointing could not be confirmed either way")
            elif jump and clear_:
                r["verdict"] = "UNCONFIRMED"
                r["why"] = (f"{hz}, on a clear-looking horizon edge, on this day only; {ftxt}: nothing confirms the "
                            f"change or rules it out (a fog front can look like a horizon)")
            elif jump:
                r["verdict"] = "UNCONFIRMED"
                r["why"] = (f"{hz}, but that is not a clear horizon ({dz}): it may be a haze, spray or foam edge "
                            f"rather than the horizon; {ftxt}: the pointing could not be confirmed either way")
            elif kind == "ambiguous":
                r["verdict"] = "small change"
                r["why"] = f"{ftxt}: a turn of {max(rup) if rup else 0.0:.2f} deg at most; {hstate}"
            elif fit_ and clear_:
                r["verdict"] = "horizon only"
                r["why"] = (f"tilt and roll within {tol} deg of {short}; the pan could not be measured ("
                            + {"not measured": "no second photo of the day to match features on",
                               "unreadable": "the photo for the feature match could not be read"}.get(
                                   f, ftxt if blind_match(r) else "features unmatched") + ")" + two)
            elif fit_:
                r["verdict"] = "unchecked"
                r["why"] = f"{hstate}; {ftxt}: nothing measured the pointing (kept)"
            else:
                r["verdict"] = "unchecked"
                r["why"] = f"neither the features nor the horizon gave a reading ({ftxt}; {hstate})"
    # a deciding match whose inliers are bunched: its why says so (not the verdict: MATCH_SPREAD_TESTED_PX)
    for r in rows:
        if trusted(r) and r.get("features") in ("ok", "moved", "MOVED") and bunched(r):
            note = bunched_text(r)
            r["why"] = (r["why"] + "; " + note) if r.get("why") else note[:1].upper() + note[1:]
    # Days nothing could read are taken in with the days left out next to them: nothing shows the camera back,
    # so such a day is not kept as if it matched (review 2, simR: c1 turned -1.2 deg on 26-28 Sep, 28 Sep's
    # 3-frame horizon fit failed at 14.6 px and its features did not match; it was 'unchecked', kept, and its
    # frames gridded).
    #   UNREAD ('unchecked'): no feature reading, and no horizon fitted or a fitted one that is not a clear
    #     horizon (one that falls within tol of the anchor by chance reads no more than a failed fit; review 4).
    #   PAN UNREAD ('horizon only'): a clear horizon at the anchor and no feature reading: it shows tilt and roll
    #     back, never the pan (review 4, P10_30_FOG01: a 1.0 deg pan measured on 30 Sep, 1 Oct fogged over the
    #     land and 'horizon only', kept with the pan on both cameras).
    #   FIRM, the days a run starts at: a day measured to differ; UNCONFIRMED on more than a horizon reading that
    #     is not a clear horizon (a clear-looking jump, an ambiguous match); or UNCONFIRMED on such a reading
    #     alone (horizon_alone) in a STORM: a run of consecutive days each left out on its own readings or unread
    #     that holds a firm day of the first kinds, misread horizons on two days or more, or a storm day
    #     (storm_days: inside STORMS, or offshore Hs over STORM_HS_M in the forcing's waves). A camera can turn in
    #     a storm (c1 turned -1.2 deg over 26-28 Sep 2026) while foam or spray hides the horizon and the land, and
    #     its days then come misread or unread in any order (review 4: real storm frames turned 1.0 deg, ST10_S3
    #     misread, unread, unread; S2728_T10_FOG misread, unread, misread: the unread days were kept, turned).
    #     One misread fog edge on a calm day never reaches past itself (review 3: it left out the run of four
    #     fog days after it).
    # An unread day next to a firm day, or to a day taken in, is UNCONFIRMED; a pan-unread one only when the
    # day its run starts at measured a pan over fail (or an ambiguous match put one forward). Runs go either
    # side (the step cannot tell when the camera turned or came back). Never the reference day: the others are
    # measured against its frame, so a neighbour that differs from it says nothing of it.
    order = sorted(rows, key=lambda r: str(r.get("date") or ""))
    storm_days = dict(storm_days or {})
    for r in order:
        r["next_to"] = ""
    own_out = set(id(o) for o in order if o.get("verdict") in LEFT_OUT_VERDICTS)     # on their own readings
    unread = set(id(o) for o in order if o.get("verdict") == "unchecked" and o is not ref)
    pan_unread = set(id(o) for o in order if o.get("verdict") == "horizon only" and o is not ref)

    def pan_of(o):
        """the pan a day's features measured, or its ambiguous match's runner-up put forward (None: neither)"""
        if trusted(o) and o.get("features") in ("moved", "MOVED"):
            return float(o["d_azimuth_deg"])
        if o.get("features") == "ambiguous" and o.get("verdict") == "UNCONFIRMED":
            p = [x for x in str(o.get("runner_up_d") or "").split("/") if x.strip()]
            return float(p[0]) if p else None
        return None

    def pan_over(o):
        p = pan_of(o)
        return p is not None and abs(p) > fail

    storm = {}                   # id(a day UNCONFIRMED on a misread horizon in a storm) -> what makes it one
    k = 0
    while k < len(order):
        j = k
        while j < len(order) and (id(order[j]) in own_out or id(order[j]) in unread):
            j += 1
        run = order[k:j]
        mis = [o for o in run if id(o) in own_out and horizon_alone(o)]
        hard = [o for o in run if id(o) in own_out and not horizon_alone(o)]
        wav = [o["date"] for o in run if o.get("date") in storm_days]
        if mis and (hard or wav or len(mis) >= 2):
            why_ = ("; ".join(f"{d_}: {storm_days[d_]}" for d_ in wav[:2]) if wav else
                    f"misread horizons on {', '.join(o['date'] for o in mis)}, in one run of days left out or "
                    f"unread" if len(mis) >= 2 else
                    f"in one run of days left out or unread with {hard[0]['date']}, left out on more than a "
                    f"misread horizon")
            for o in mis:
                storm[id(o)] = why_
        k = max(j, k + 1)

    def firm(o):
        return (o.get("verdict") == "DIFFERENT POINTING"
                or (o.get("verdict") == "UNCONFIRMED" and (not horizon_alone(o) or id(o) in storm)))

    def firm_text(o):
        p = pan_of(o)
        pan_ = (f" (a pan of {p:+.2f} deg, which no horizon sees)" if p is not None and abs(p) > fail else "")
        if o.get("verdict") == "DIFFERENT POINTING":
            return "was measured to differ from the reference" + pan_
        if horizon_alone(o):
            return (f"is UNCONFIRMED on a misread horizon in a storm ({storm[id(o)]}; a camera can turn in one: c1 "
                    f"turned -1.2 deg of pan on 27 Sep 2026, measured; on 26 and 28 Sep by tracked landmarks, "
                    f"outside this code)")
        return (f"is UNCONFIRMED on an ambiguous feature match (the camera may have turned{pan_})"
                if o.get("features") == "ambiguous" else "is UNCONFIRMED on a horizon jump on a clear-looking edge")

    src_of, pan_src = {}, {}     # id(a day taken in) -> the firm day its run starts at; -> one with a pan over fail
    grown = True
    while grown:
        grown = False
        for i, r in enumerate(order):
            if id(r) not in unread and id(r) not in pan_unread:
                continue
            got, gotp = [], []
            for k in (i - 1, i + 1):
                if not 0 <= k < len(order):
                    continue
                o = order[k]
                if id(o) in src_of:
                    got.append(src_of[id(o)])
                    if id(o) in pan_src:
                        gotp.append(pan_src[id(o)])
                elif firm(o):
                    got.append(o)
                    if pan_over(o):
                        gotp.append(o)
            if id(r) in pan_unread:
                got = gotp           # its clear horizon shows tilt and roll back, never the pan
            if not got:
                continue
            if id(r) not in src_of:
                src_of[id(r)] = sorted(got, key=lambda o: o.get("verdict") != "DIFFERENT POINTING")[0]
                grown = True
            if gotp and id(r) not in pan_src:
                pan_src[id(r)] = gotp[0]
                grown = True
    for i, r in enumerate(order):
        src = src_of.get(id(r))
        if src is None:
            continue
        nb = any(0 <= k < len(order) and order[k] is src for k in (i - 1, i + 1))
        hz_only = id(r) in pan_unread
        r["verdict"], r["next_to"] = "UNCONFIRMED", src["date"]
        r["why"] = (r["why"].replace(": nothing measured the pointing (kept)", ": nothing measured the pointing")
                    + "; and "
                    + (f"the day next to it, {src['date']}, {firm_text(src)}" if nb else
                       f"it is in a run of days nothing could read{', or only the horizon,' if hz_only else ''} "
                       f"next to {src['date']}, which {firm_text(src)}")
                    + (": its clear horizon shows tilt and roll back but cannot see the pan, so nothing shows the "
                       "camera back; the pointing could not be confirmed either way" if hz_only else
                       ": nothing shows the camera back, so the pointing could not be confirmed either way"))
    # the window's other pointing: every day measured to differ agrees with the others, and they outnumber the
    # days measured at the reference's pointing (the reference itself may be the odd one out)
    def change(r):
        """a DIFFERENT POINTING day's change from the reference: [pan or None (not measured), tilt, roll]"""
        if trusted(r) and r.get("features") in ("moved", "MOVED"):
            return [float(r["d_azimuth_deg"]), float(r["d_tilt_deg"]), float(r["d_roll_deg"])]
        h0, h1 = hz_of(r)
        return [None, h0 - at, h1 - ar]

    def agree(p, q):
        return (abs(p[1] - q[1]) <= tol and abs(p[2] - q[2]) <= tol
                and (p[0] is None or q[0] is None or abs(p[0] - q[0]) <= fail))

    moved_ = [(r, change(r)) for r in rows if r.get("verdict") == "DIFFERENT POINTING"]
    same_ = [r["date"] for r in rows if r.get("verdict") in ("matches", "small change", "horizon only")]
    # days that read nothing else, whose MATCH_INLIERS or more features matched the reference at zero turn within
    # their own noise (blind_match): no reading of a turn, but what the reference's own pointing gives (or dirt on
    # the window): kept 'unchecked', or UNCONFIRMED only for the day next to them. They are counted against the
    # other pointing's days (review r1: the Oct 2024 fixture's c2, two planted +1.0 deg pan days against the
    # reference day alone, while three days matched it at 0.00 deg on 222-236 features: 'the reference's is the
    # smaller group' and a PARTIAL build). A day whose clear horizon jumps is not one of them.
    zero_ = [r["date"] for r in rows if r is not ref and blind_match(r) and _inliers(r) >= MATCH_INLIERS
             and (r.get("verdict") == "unchecked" or (r.get("verdict") == "UNCONFIRMED" and r.get("next_to")))]
    other = None
    if moved_ and len(moved_) > len(same_) + len(zero_) and all(agree(p, q) for _, p in moved_ for _, q in moved_):
        pans = [p[0] for _, p in moved_ if p[0] is not None]
        other = {"days": [r["date"] for r, _ in moved_], "reference_days": same_, "zero_turn_days": zero_,
                 "change_deg": [round(float(np.median(pans)), 3) if pans else None,
                                round(float(np.median([p[1] for _, p in moved_])), 3),
                                round(float(np.median([p[2] for _, p in moved_])), 3)]}
    anchor_note = None
    if how == "tied" and ref is not None and any(ref is r for r in tied) and not any(ref is r for r in grp):
        anchor_note = (f"the reference day {ref_day}'s own clear horizon reads {hz_of(ref)[0]:+.2f}/{hz_of(ref)[1]:+.2f} "
                       f"deg, more than {tol} deg from the {len(grp)} other clear day(s) with matched features "
                       f"({', '.join(r['date'] for r in grp)}: {at:+.2f}/{ar:+.2f} deg), which the days are judged "
                       f"against; if those are a fog front rather than the horizon, the verdicts are wrong: look at "
                       f"the photos")
    if clear:
        mt, mr = med([hz_of(r) for r in base]) if base else anchor
        offs = [float(r["horizon_offset_px"]) for r in base]
    else:
        (mt, mr), offs = anchor, []
    offs_f = [float(r["horizon_offset_px"]) for r in fitted]
    return {"window_horizon_dtilt_deg": round(mt, 3), "window_horizon_droll_deg": round(mr, 3),
            "clear_horizon_days": len(base) if clear else 0,
            "window_horizon_offset_px": round(float(np.median(offs)), 1) if offs else None,
            "fitted_horizon_days": len(fitted),
            # with no clear day: the fitted readings' median offset, reported as such, never as the horizon's
            "fitted_horizon_offset_px": round(float(np.median(offs_f)), 1) if offs_f and not clear else None,
            "horizon_window_from": ("clear" if clear else "fitted, none clear" if fitted else None),
            # what each day's horizon was judged against (the anchor above) and where it came from
            "horizon_anchor_deg": [round(at, 3), round(ar, 3)],
            "horizon_anchor_from": {"tied": "clear days with matched features (horizon less their measured turn)",
                                    "reference": f"the reference day {ref_day}'s clear horizon",
                                    "clear": "the largest group of clear days that agree (nothing ties it to the "
                                             "reference)",
                                    "fitted": "the fitted days (none clear)"}[how],
            "horizon_other_clear_days": [r["date"] for r in others],
            "horizon_anchor_note": anchor_note,
            "horizon_not_clear_days": [r["date"] for r in fitted if id(r) not in is_clear],
            "horizon_contradicted_days": [r["date"] for r in contra],
            "horizon_overruled_days": overruled,
            # every day measured to differ at one other pointing that outnumbers the reference's (None otherwise)
            "other_pointing": other}


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
    text = ("the sea horizon sits " + (f"{px:.0f} px " if px is not None else "")
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


# The pointing step's rule (pointing_verdicts) is part of its result, so a result made by another
# version of the rule is rebuilt: 5 (before Oct 2026), 6 (matched features decide; edge numbers),
# 7 (each day judged against the reference; a zero turn never overrules a clear horizon jump),
# 8 (a turned match's tilt/roll judged by RUN_TOL_DEG; a zero turn within its own noise shows nothing;
# the reference day's confirmed group anchors; a day nothing read next to a left-out day is UNCONFIRMED),
# 9 (that day only next to a firm one, never a horizon reading alone; a turned match under 100 inliers
# with a clear horizon over RUN_TOL_DEG; 'most days differ' no longer keeps days: the other pointing named),
# 10 (storms: a run of days left out or unread with two misread days, a firm day or a storm day by the
# waves or STORMS; a fitted horizon that is not clear reads nothing; a 'horizon only' day next to a
# measured pan; a turned match confirmed by a clear horizon within MATCH_ERR_DEG at any inlier count),
# 11 (the other pointing is counted against the zero-turn days too; on the calibration's own day it is no
# reason for a partial build and no rebuild is advised).
POINTING_VERSION = 11


def pointing_sig(plan, cam, args, photos):
    c = plan["cams"][cam]
    sig = {"eo": c["eo_file"], "eo_file": file_sig(c["eo_path"]),
           "photos": hashlib.sha1("\n".join(p.name for p in photos).encode()).hexdigest(),
           "frames": args.pointing_frames, "fail": args.pointing_fail, "v": POINTING_VERSION,
           # the window's storm days (from the forcing's waves): a new wave record reruns the verdicts only
           # when it changes them
           "storm_days": storm_days_of(plan, per_day(photos)),
           "scripts": scripts_sig("pointing")}
    if getattr(args, "correct_pointing", False):
        # --correct-pointing also matches the first and last photos of a day it would correct (pointing_day_check);
        # a default build's signature is as it was
        import pointing_daily
        sig["day_check"] = {"v": pointing_daily.CORRECTION_VERSION, "rules": pointing_correction_rules(),
                            "script": script_sha("pointing_daily.py")}
    return sig


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
        if fresh:
            fresh, why = eo_content_fresh(out, f"pointing_{cam}", c["eo_path"])
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
                        scripts=STEP_SCRIPTS["pointing"], extra={"photo_roots": list(args.photo_roots),
                                                                  "eo_sha256": calibration_record(c)["eo_sha256"],
                                                                  "calibration": calibration_record(c)})
            continue
        if not c["eo_path"]:
            raise StepFailed(f"{cam}: calibration {c['eo_file']} not found in {args.calibration}")
        rows, info = check_pointing(plan, cam, photos, args)
        lo = pointing_leave_out(rows)
        # Only what was MEASURED is stored; whether those days are left out is decided when the
        # detection runs (--keep-moved-days), so changing that option takes effect on a rerun.
        warn_ = other_pointing_warning(info, lo, pdir, plan["date"])
        if warn_:
            info["warning"] = warn_
        fields = ["camera", "date", "photos", "verdict", "why", "horizon", "horizon_offset_px",
                  "horizon_dtilt_deg", "horizon_droll_deg", "horizon_fit_px", "horizon_frames", "horizon_edge",
                  "horizon_prominence", "horizon_rival", "horizon_sky_above", "horizon_to_limit_px", "features",
                  "feature_frame", "inliers", "d_azimuth_deg", "d_tilt_deg", "d_roll_deg", "turn_sigma",
                  "inlier_spread_px", "runner_up_inliers", "runner_up_d", "next_to"]
        with open(csv_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        res = {"rows": rows, "info": info, "different_days": lo["different_days"], "leave_out": lo["leave_out"],
               "measured_different_days": lo["measured_different_days"], "unconfirmed_days": lo["unconfirmed_days"],
               "horizon_alone_days": lo["horizon_alone_days"], "next_to_days": lo["next_to_days"],
               "window_days": lo["window_days"],
               "horizon_overruled_days": info.get("horizon_overruled_days") or [],
               "station_log": station_pointing_log(plan, cam, args) if plan["era"] == "live" else None}
        js.write_text(json.dumps(res, indent=1) + "\n")
        write_stamp(out, f"pointing_{cam}", sig, [f"survey_products.py, in-process: check_pointing() for {cam} "
                                                  f"(horizon_check.fit_tilt_roll + pointing_check.measure)"],
                    scripts=STEP_SCRIPTS["pointing"], extra={"photo_roots": list(args.photo_roots),
                                                              "eo_sha256": calibration_record(c)["eo_sha256"],
                                                              "calibration": calibration_record(c)})
        results[cam] = res
        any_built = True
    for cam, res in results.items():
        print_pointing(cam, res, args, detect_now="detect" in (state.get("requested") or STEPS),
                       detected=read_stamp(out, f"detect_{cam}"))
    state["pointing"] = results
    if results and all(not r.get("rows") for r in results.values()):
        return "not checked (no photos found)"
    return "built" if any_built else "skipped"


def pointing_leave_out(rows):
    """Which days the detection leaves out by default, from the verdicts (pointing_verdicts):
    'different_days' (the key keeps its name: older builds and code read it) holds the measured
    DIFFERENT POINTING and the UNCONFIRMED days alike, whatever share of the window they are:
    pointing_majority_note says so when they are most of it, and other_pointing_warning when the days
    measured to differ agree with each other and outnumber those at the reference's pointing (the
    reference itself may be the odd one out). Until Oct 2026 more than half of the measured days
    differing from the reference kept EVERY day ('leave_out' False: nothing left out), the UNCONFIRMED
    storm days and days that differ from each other included, which mixed two pointings or more in
    the product (review 3: an exact +0.6 deg tilt on c1 30 Sep - 2 Oct kept all 7 days, 26-28 Sep
    with them). 'leave_out' is now always True; an older result's False still means what its
    detection did (skip_days_of)."""
    differ = [r["date"] for r in rows if r.get("verdict") == "DIFFERENT POINTING"]
    measured = [r for r in rows if r.get("verdict") not in ("unchecked", "UNCONFIRMED")]
    return {"different_days": [r["date"] for r in rows if r.get("verdict") in LEFT_OUT_VERDICTS],
            "measured_different_days": differ,
            "unconfirmed_days": [r["date"] for r in rows if r.get("verdict") == "UNCONFIRMED"],
            "horizon_alone_days": [r["date"] for r in rows if horizon_alone(r)],
            "next_to_days": [r["date"] for r in rows if unread_next_to(r)],
            "measured_days": len(measured),
            "window_days": len(rows),
            "leave_out": True}


def other_pointing_window(op, ref_day):
    """(FIRST, LAST) of a --window that holds the other pointing's days (pointing_verdicts' 'other_pointing')
    on one side of the reference day and leaves that day out (the side with more of them; ties: after it),
    or None. A build with that window takes the day nearest the survey in it as its reference."""
    if not op or not ref_day:
        return None
    before = sorted(d for d in op.get("days") or [] if d < ref_day)
    after = sorted(d for d in op.get("days") or [] if d > ref_day)
    side = after if len(after) >= len(before) else before
    return (side[0], side[-1]) if side else None


def other_pointing_rebuild(op, ref_day, date=None):
    """The concrete way to build on the window's other pointing instead of the reference's, or ''."""
    w = other_pointing_window(op, ref_day)
    if not w:
        return ""
    return ("to build on the other pointing instead, leave the reference day out: survey_products.py"
            + (f" --date {date}" if date else "") + f" --window {w[0]} {w[1]} --output-root <another root> (all "
            f"steps; the day nearest the survey in that window becomes the reference; in this root it would "
            f"replace this build)")


def reference_is_calibration_day(info):
    """Is the pointing step's reference day the calibration's own solve day (the date in its file name, e.g.
    CACO03_c2_20241023_EO.yaml on 2024-10-23; check_pointing takes that day as the reference whenever it has
    photos)? Its pointing is then the calibration's: a window's other pointing, however many days hold it, can
    never make the reference the wrong pointing FOR THIS CALIBRATION (its days would need a calibration of their
    own). From check_pointing's info (an older result's too: it records the calibration and the reference day)."""
    rd = (info or {}).get("reference_day")
    return bool(rd) and eo_date((info or {}).get("calibration") or "") == rd


def other_pointing_groups(op):
    """The days at the reference's pointing against the other pointing's, in words: those measured there, and
    those whose features matched the reference at zero turn within their own noise (no reading, counted with
    them: pointing_verdicts' 'zero_turn_days')."""
    ref, zero = op.get("reference_days") or [], op.get("zero_turn_days") or []
    return ((", ".join(ref) or "no day measured")
            + (f"; with {', '.join(zero)}, whose features matched it at zero turn within their own noise (no reading)"
               if zero else ""))


def other_pointing_warning(info, lo, where, date=None):
    """The step's (and the README's) warning when every day measured to differ from the reference
    agrees with the others and they outnumber the days at the reference's pointing (pointing_verdicts'
    'other_pointing': the days measured there and the zero-turn days): the window holds two pointings and
    the reference's is the smaller group, so the reference itself may be the odd one out. When the reference
    is the calibration's own day (reference_is_calibration_day) it cannot be: the warning says the product
    rests on the calibration's pointing and advises no rebuild. None otherwise."""
    op = (info or {}).get("other_pointing")
    if not op:
        return None
    pan, dt, dr = op["change_deg"]
    ref = op.get("reference_days") or []
    head = (f"{len(op['days'])} of the {lo.get('measured_days', len(op['days']) + len(ref))} days measured differ "
            f"from the reference {info.get('reference_day')} and agree with each other ("
            + (f"{pan:+.2f}/{dt:+.2f}/{dr:+.2f} deg pan/tilt/roll" if pan is not None else
               f"{dt:+.2f}/{dr:+.2f} deg tilt/roll, pan not measured")
            + f" from it, median; {', '.join(op['days'])}): the window holds two pointings")
    if reference_is_calibration_day(info):
        return (head + f", and more of its days at the other one than at the reference's ({other_pointing_groups(op)})"
                f". The reference day "
                f"is the calibration's own day ({info.get('calibration')}), so its pointing is the calibration's (only "
                f"a turn on that very day, between the calibration's photos and the reference frame "
                f"{info.get('reference_frame') or ''}, would make it otherwise): "
                f"the detection leaves the other pointing's days out by default (with this calibration they would "
                f"be gridded with that turn as a pointing error), and the product rests on the calibration's "
                f"pointing; look at the photos ({where}). --keep-moved-days keeps every day (both pointings, the "
                f"UNCONFIRMED days too)")
    return (head + ", and "
            + (f"the reference's ({other_pointing_groups(op)}) is the smaller group"
               if op.get("reference_days") or op.get("zero_turn_days") else
               "no day was measured at the reference's pointing")
            + ". The detection leaves the other "
            f"pointing's days out by default, so the product rests on the reference's pointing; if the "
            f"reference day itself is the odd one out, that is the wrong one: look at the photos ({where}). "
            + ((other_pointing_rebuild(op, info.get("reference_day"), date)[:1].upper()
                + other_pointing_rebuild(op, info.get("reference_day"), date)[1:] + ". ")
               if other_pointing_window(op, info.get("reference_day")) else "")
            + "--keep-moved-days keeps every day (both pointings, the UNCONFIRMED days too), and --window FIRST "
            "LAST narrows the photo days; there is no option that keeps only some days")


def skip_days_of(res, args):
    """The days the detection leaves out, from a pointing result and this run's options (None: the
    defaults). A result of the rule before Oct 2026 whose 'most days differ' guard kept every day
    ('leave_out' False) gives none, as its detection left out none. With --correct-pointing the days it
    corrects (pointing_corrections) are detected, not left out."""
    if not res or getattr(args, "keep_moved_days", False) or not res.get("leave_out", True):
        return []
    skip = sorted(res.get("different_days") or [])
    if getattr(args, "correct_pointing", False):
        fixed = pointing_corrections(res)[0]
        skip = [d for d in skip if d not in fixed]
    return skip


def pointing_correction_rules():
    """pointing_daily's gates (correction_rules) with this script's thresholds: MATCH_INLIERS, TURN_SIGMA_MIN and
    MATCH_SPREAD_TESTED_PX, its own CORRECT_MAX_DEG, and the gates of the day at MATCH_ERR_DEG (a turn over it; the
    day's other frames within it of the probe's turn)."""
    import pointing_daily
    return pointing_daily.correction_rules(MATCH_INLIERS, TURN_SIGMA_MIN, MATCH_SPREAD_TESTED_PX,
                                           min_turn_deg=MATCH_ERR_DEG, max_day_spread_deg=MATCH_ERR_DEG)


def pointing_day_check(day_photos, probe, ref_src, ref, io, row, ref_day, args):
    """
    --correct-pointing (pointing_daily.day_held): did the pointing the probe measured hold through the day? For a day
    whose probe passes every other gate, its first and last photos of the detection (other than the probe and the
    reference frame) are matched against the same reference; the other days are not corrected, so they are not
    matched again. -> the columns added to the day's row, feature_time_utc and day_check
    ([feature_record + time_utc]; [] when the probe is the day's only photo), or {} for a day not checked.
    """
    import cv2
    import pointing_check as pc
    import pointing_daily
    if probe is None:
        return {}
    rules = {k: v for k, v in pointing_correction_rules().items() if k != "max_day_spread_deg"}
    if pointing_daily.day_correction(row, rules, ref_day)[0] is None:
        return {}
    others = sorted((p for p in day_photos if p.name not in (probe.name, ref_src.name)), key=lambda p: utc_of(p.name))
    chk = []
    for p in (others[:1] + others[-1:]) if len(others) > 1 else others:
        img = cv2.imread(str(p))
        m = pc.measure(img, [ref], io, min_inliers=MATCH_INLIERS, second=True) if img is not None else None
        chk.append(dict(feature_record(m, p, args.pointing_fail), time_utc=f"{utc_of(p.name):%H:%M}"))
    return {"feature_time_utc": f"{utc_of(probe.name):%H:%M}", "day_check": chk}


def pointing_corrections(res):
    """--correct-pointing: the days of a pointing result that pass pointing_daily's gates, {day: the change measured
    and its match}, and why each other day does not, {day: [why]}; judged against the result's reference day.
    ({}, {}) without rows."""
    import pointing_daily
    if not res or not res.get("rows"):
        return {}, {}
    return pointing_daily.window_corrections(res["rows"], pointing_correction_rules(),
                                             (res.get("info") or {}).get("reference_day"))


def corrected_days(dst):
    """The days a detection mapped with their own pointing (its stamp's pointing_correction; none without it)."""
    return sorted(((dst or {}).get("pointing_correction") or {}).get("days") or {})


def pointing_correction_plan(plan, state, cam, args):
    """--correct-pointing: what the detection of this camera corrects, from the pointing results it uses
    (pointing_skip_days has checked they are this window's and this calibration's): {"days": {day: change},
    "not_corrected": {day: [why]}, "res": the results}; None without the option."""
    if not getattr(args, "correct_pointing", False):
        return None
    res = (state.get("pointing") or {}).get(cam) or load_json(plan["out"] / "pointing" / f"pointing_{cam}.json")
    corr, why_not = pointing_corrections(res)
    return {"days": corr, "not_corrected": why_not, "res": res}


def pointing_correction_sig(pc_):
    """The part of a corrected detection's signature that the correction adds: its version, gates, the days and
    their changes, pointing_daily.py's sha256. Only with --correct-pointing: a default build's signature is as it
    was, so nothing rebuilds for this code."""
    import pointing_daily
    return {"v": pointing_daily.CORRECTION_VERSION, "rules": pointing_correction_rules(),
            "days": {d: v["d_deg"] for d, v in sorted(pc_["days"].items())}, "script": script_sha("pointing_daily.py")}


def correct_pointing_rows(plan, cam, ground, pc_):
    """--correct-pointing, once the camera's waterlines are georectified with the calibration: each corrected day's
    rows of `ground` mapped again with that day's own pointing (pointing_daily.apply_corrections; its EO files in
    waterlines/<cam>/pointing_corrected_eo/). -> the record the detection's stamp keeps, which the README and
    provenance read."""
    import pointing_daily
    from georectify import load_intrinsics, load_extrinsics
    c = plan["cams"][cam]
    info = (pc_.get("res") or {}).get("info") or {}
    days = pointing_daily.apply_corrections(ground, cam, load_intrinsics(c["io_path"]), load_extrinsics(c["eo_path"]),
                                            c["eo_file"], pc_["days"], Path(ground).parent / "pointing_corrected_eo",
                                            reference=info.get("reference_frame") or "")
    say(f"pointing correction {cam}", (", ".join(
        f"{d} {v['d_deg'][0]:+.2f}/{v['d_deg'][1]:+.2f}/{v['d_deg'][2]:+.2f} deg ({v['rows']} points"
        + (f", moved {v['moved_m_median']:.2f} m median" if v.get("moved_m_median") is not None else "") + ")"
        for d, v in sorted(days.items())) or "no day passes the gates")
        + " (--correct-pointing: mapped with each day's own pointing)")
    return {"v": pointing_daily.CORRECTION_VERSION, "rules": pointing_correction_rules(),
            "reference_day": info.get("reference_day"), "reference_frame": info.get("reference_frame"),
            "days": days, "not_corrected": pc_["not_corrected"]}


def pointing_eo_day_specs(plan, cams):
    """--correct-pointing: survey_compare.py's --camera-eo-day CAM:DAY=PATH for each day the detection of one of
    `cams` mapped with its own pointing (its stamp's record; the EO files in waterlines/<cam>/pointing_corrected_eo/),
    so the C = 0 sensitivity re-projects those days' lines with the pointing they were mapped with, not the
    calibration's. [] for a default build (its compare command is as it was)."""
    specs = []
    for cam in sorted(cams):
        pcr = (read_stamp(plan["out"], f"detect_{cam}") or {}).get("pointing_correction") or {}
        for d, v in sorted((pcr.get("days") or {}).items()):
            specs.append(f"{cam}:{d}={plan['out'] / 'waterlines' / cam / 'pointing_corrected_eo' / v['eo_file']}")
    return specs


def pointing_root_problem(plan, args):
    """--correct-pointing builds in its own root, named for it (<root>_pointing-corrected, optionally + _c<C>), and no
    other build goes in such a root, so a corrected product never replaces a default one and a root's name says
    which it holds; not on a BLIND date, not with --eo. -> text ('' when fine). Refused, exit 5."""
    on = bool(getattr(args, "correct_pointing", False))
    root = Path(str(getattr(args, "output_root", "") or ""))
    named = bool(re.search(rf"_{re.escape(POINTING_CORRECTED)}(_c\d+(\.\d+)?)?$", root.name))
    if on and plan.get("blind"):
        return (f"--correct-pointing is not for a BLIND date: its builds (A, the table's calibration, and the GCP-free "
                f"B) are fixed by its protocol")
    if on and any(c.get("eo_given") for c in plan["cams"].values()):
        return "--correct-pointing does not go with --eo: a calibration given there is built as it was made"
    if on and not named:
        return (f"a --correct-pointing build goes in its own root, named for it: --output-root "
                f"{root.parent / (root.name + '_' + POINTING_CORRECTED)} (not {root}: a corrected product never "
                f"replaces a default one, and the root's name says which it holds)")
    if named and not on:
        return (f"{root} is named for pointing-corrected products: build there with --correct-pointing, or build this "
                f"one in another --output-root")
    return ""


def pointing_majority_note(res, args, skip=None, ran=True):
    """A warning when the days left out for their pointing are more than half of the window's (the
    product then rests on few days): the count, how many of them rest on a horizon reading alone that
    is not a clear horizon or are left out only for the day next to them, and what to do. None
    otherwise. (Days measured to differ that agree with each other have their own warning,
    other_pointing_warning, in info['warning'].) skip: the days the detection left out (its stamp;
    None: as this run's options leave them out); ran=False: the detection has not run with these
    results, and the days are those it would leave out."""
    if not res or not res.get("rows"):
        return None
    n = res.get("window_days") or len(res["rows"])
    skip = sorted(skip_days_of(res, args) if skip is None else skip)
    if len(skip) <= n / 2.0:
        return None
    alone = [d for d in (res.get("horizon_alone_days") or [r["date"] for r in res["rows"] if horizon_alone(r)])
             if d in skip]
    nxt = [d for d in (res.get("next_to_days") or [r["date"] for r in res["rows"] if unread_next_to(r)])
           if d in skip]
    return (f"{len(skip)} of the window's {n} days "
            + ("are left out of the detection" if ran else
               "would be left out of the detection (it has not run with these pointing results yet)")
            + f" for their pointing ({', '.join(skip)}): the product "
            + ("rests" if ran else "would rest") + f" on the {n - len(skip)} other day(s)"
            + (f". {len(alone)} of them ({', '.join(alone)}) {'is' if len(alone) == 1 else 'are'} UNCONFIRMED on a "
               f"horizon reading alone that is not "
               f"a clear horizon, with no feature reading: such a reading never shows a day moved, but in a storm "
               f"the camera can turn on those very days (c1, 27 Sep 2026: -1.2 deg of pan, measured; 26 and 28 Sep by "
               f"tracked landmarks, outside this code)" if alone else "")
            + (f". {len(nxt)} of them ({', '.join(nxt)}) {'is' if len(nxt) == 1 else 'are'} UNCONFIRMED only because "
               f"nothing could read {'it' if len(nxt) == 1 else 'them'} (no feature reading, and no clear horizon; "
               f"or only a clear horizon, next to a measured pan it cannot see) next to a day left out (measured to "
               f"differ, unconfirmed on more than a misread horizon, or misread in a storm: two misread days, or "
               f"high waves): nothing shows the camera back" if nxt else "")
            + ". Look at the photos and the pointing CSV; --keep-moved-days keeps every such day (there is no "
              "option that keeps only some; --window FIRST LAST narrows the photo days)")


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


def print_pointing(cam, res, args, detect_now=True, detected=None):
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
    # why, for every day that is not plainly the window's pointing (and every day a horizon
    # reading was overruled): what was measured, and what the step concluded from it
    for r in res["rows"]:
        if r.get("why") and r.get("verdict") != "horizon only":
            print(f"    {r['date']} {r['verdict']}: {r['why']}")
    off = pointing_offset(info, info.get("calibration"))
    if off:
        warn(f"{cam}: {off[2]}")
    if info.get("warning"):
        warn(f"{cam}: {info['warning']}")
    if info.get("horizon_anchor_note"):
        warn(f"{cam}: {info['horizon_anchor_note']}")
    note = pointing_majority_note(res, args)
    if note:
        warn(f"{cam}: {note}")
    fixed = pointing_corrections(res)[0] if getattr(args, "correct_pointing", False) else {}
    if fixed:
        verdict = {r["date"]: r.get("verdict", "?") for r in res["rows"]}
        say(f"{cam} corrected", ("" if detect_now else "when the detection runs again: ")
            + ", ".join(f"{d} ({verdict.get(d, '?')}, {v['d_deg'][0]:+.2f}/{v['d_deg'][1]:+.2f}/{v['d_deg'][2]:+.2f} "
                        f"deg)" for d, v in sorted(fixed.items()))
            + " (--correct-pointing: detected, and mapped with the calibration turned by the change matched features "
              "measured that day)")
    if res.get("different_days"):
        skip = skip_days_of(res, args)
        kept = [d for d in res["different_days"] if d not in skip and d not in fixed]
        verdict = {r["date"]: r.get("verdict", "?") for r in res["rows"]}
        # what the detection does with these results: this run's when it runs the detection too; a pointing-only
        # run (--steps pointing) says what the next detection will do and what the one on disk did
        later = "" if detect_now else (
            " when the detection runs again with these results (--steps detect and the steps after it; "
            + (f"the detection on disk left out {', '.join(sorted((detected or {}).get('skip_days') or [])) or 'no day'})"
               if detected else "no detection has run yet)"))
        if skip:
            warn(f"{cam}: " + ("LEFT OUT of the detection: " if detect_now else "to be LEFT OUT of the detection: ")
                 + ", ".join(f"{d} ({verdict.get(d, '?')})" for d in skip)
                 + " (why: above; --keep-moved-days keeps them)" + later)
        if kept:
            warn(f"{cam}: days with a different or unconfirmed pointing " + ("KEPT: " if detect_now else "to be KEPT: ")
                 + ", ".join(f"{d} ({verdict.get(d, '?')})" for d in kept) + f" ({kept_reason(res, args)})" + later)


def kept_reason(res, args, detect=None):
    """Why days with a different or unconfirmed pointing are kept in the detection (skip_days_of):
    from this run's options, or, given the detection's stamp (detect), from the options it ran with.
    The option is never inferred: a stamp of an older build did not record --keep-moved-days, and
    then only an older rule's 'most days differ' guard (the pointing result's 'leave_out' False) is
    named; else its record does not say (the README calls such a detection stale before this)."""
    kmd = getattr(args, "keep_moved_days", False) if detect is None else detect.get("keep_moved_days")
    if kmd:
        return "--keep-moved-days"
    if not res.get("leave_out", True):
        return "most days differ from the reference, which kept every day by the rule before Oct 2026: see the warning"
    if kmd is None:
        return ("--keep-moved-days, or other pointing results: the detection's record, of an older build, does "
                "not say which")
    return "not by its options: run the detection again"


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
        want = normal(pointing_sig(plan, cam, args, photos))
        if not st or st.get("sig") != want:
            older = bool(st) and (st.get("sig") or {}).get("v") != want.get("v")
            storms_ = bool(st) and (st.get("sig") or {}).get("storm_days") != want.get("storm_days")
            raise StepFailed(f"{cam}: the pointing check in {js} was made "
                             + (f"by an older version of the check (v{(st.get('sig') or {}).get('v')}; this code "
                                f"decides by v{want.get('v')}: which days are left out may differ)" if older else
                                "for other storm days (the forcing's waves changed since: they decide whether a "
                                "day nothing read next to a misread one is left out)" if storms_ else
                                "for other photos or another calibration")
                             + "; run the pointing step again (--steps pointing,detect,...)")
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
            # The water-level timing (timex_timing.py) reaches the detection through
            # detect_original_view.py's detector and contour subprocesses (TIMING_MODE). A mode
            # other than legacy is part of the signature, so a switch rebuilds the detection;
            # legacy adds nothing to it (no key), so going back to legacy matches legacy stamps.
            from timex_timing import timing_mode
            try:
                timing = timing_mode()
            except ValueError as e:
                raise StepFailed(str(e))
            if timing != "legacy":
                sig["timing_mode"] = timing
            pc_ = pointing_correction_plan(plan, state, cam, args)
            if pc_ is not None:
                sig["pointing_correction"] = pointing_correction_sig(pc_)
            fresh, why = is_fresh(out, f"detect_{cam}", sig, [ground],
                                  [f / "water_level.csv", f / "waves.csv", c["eo_path"]])
            if fresh:
                fresh, why = eo_content_fresh(out, f"detect_{cam}", c["eo_path"])
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
            extra = {"setup_coef": args.setup_coef, "setup_fit": fit, "skip_days": skip,
                     "keep_moved_days": bool(getattr(args, "keep_moved_days", False)),
                     "photos": len(photos), "eo_sha256": calibration_record(c)["eo_sha256"],
                     "calibration": calibration_record(c)}
            if pc_ is not None:
                try:
                    extra["pointing_correction"] = correct_pointing_rows(plan, cam, ground, pc_)
                except (OSError, ValueError, KeyError) as exc:
                    warn(f"detect {cam}: the pointing correction FAILED ({exc.__class__.__name__}: {exc}): camera "
                         f"left out")
                    state.setdefault("failed_cams", {})[cam] = f"pointing correction failed ({exc})"
                    continue
            write_stamp(out, f"detect_{cam}", sig, [cmd], scripts=STEP_SCRIPTS["detect"]
                        + (["pointing_daily.py"] if pc_ is not None else []), extra=extra)
            # the search envelope's content as this detection used it (a GCP-free B of a blind date must use A's)
            update_stamp(out, f"detect_{cam}", envelope_sha256=file_sha(c["envelope_path"]))
            if not args.keep_detector_debug:
                prune_detector_scratch(wdir / cam)
            parts.append((cam, ground))
    if not parts:
        raise StepFailed("no camera produced waterlines")
    merged = wdir / "contour_points_ground.csv"
    # v 3 (historical eras): water_level_source / tide_model_spread now name the forcing's level source;
    # the live era's merge is unchanged (v 2), so a built live date does not rebuild for it
    sig = {"parts": [[c, str(p)] for c, p in parts], "era": plan["era"], "v": 2 if plan["era"] == "live" else 3,
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
        pc_ = pointing_correction_plan(plan, state, cam, args)
        if pc_ is not None:
            sig["pointing_correction"] = pointing_correction_sig(pc_)
        fresh, why = is_fresh(out, f"detect_{cam}", sig, [ground], [rows_file, c["eo_path"]])
        if fresh:
            fresh, why = eo_content_fresh(out, f"detect_{cam}", c["eo_path"])
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
        extra = {"setup": {k: v for k, v in stats.items() if k != "frame_names"},
                 "setup_coef": args.setup_coef, "setup_fit": fit, "skip_days": skip,
                 "keep_moved_days": bool(getattr(args, "keep_moved_days", False)),
                 "eo_sha256": calibration_record(c)["eo_sha256"], "calibration": calibration_record(c)}
        if pc_ is not None:
            try:
                extra["pointing_correction"] = correct_pointing_rows(plan, cam, ground, pc_)
            except (OSError, ValueError, KeyError) as exc:
                warn(f"georectify {cam}: the pointing correction FAILED ({exc.__class__.__name__}: {exc}): camera "
                     f"left out")
                state.setdefault("failed_cams", {})[cam] = f"pointing correction failed ({exc})"
                continue
        # photos of the frames, for the maps (the maps step links them again from the photo roots of
        # its own run, so a photo that turns up later is used without a new georectification)
        link_frame_photos(plan, cam, args)
        write_stamp(out, f"detect_{cam}", sig, [geo], scripts=STEP_SCRIPTS["detect_live"]
                    + (["pointing_daily.py"] if pc_ is not None else []), extra=extra)
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
        # the slope-aware model's columns (the cron's SETUP_MODEL="stockdon") go too: these rows get C's setup
        fields = [k for k in rd.fieldnames if k not in ("setup_correction_m", "beach_elevation_navd88",
                                                        "setup_beta_f", "setup_model")]
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
    x = wave_scale(h, p)
    return s / x if x > 0.05 else None


def merge_cameras(plan, parts, args):
    """
    waterlines/contour_points_ground.csv: every camera's rows, each tagged with its water-level and
    wave source (forcing_level_source, forcing_level_sigma_m, forcing_waves_source); outside the live
    era water_level_source is the forcing's level source too (adcp / chatham_transfer, not the
    extractor's 'tide_model') and tide_model_spread is blank (no tide model was used).
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
                        # extract_elevation_contours.py calls every --water-level CSV a 'tide_model' (with a
                        # 0.0 spread across model columns): here the level is the forcing's (ADCP or Chatham
                        # transfer), so the row says that, not that a tide model gave it
                        if "water_level_source" in r:
                            r["water_level_source"] = cache[k][0] or "historical_forcing"
                        if "tide_model_spread" in r:
                            r["tide_model_spread"] = ""
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

def remove_filter_outputs(wdir, dst):
    """Everything waterline_consistency.py writes for a build: the filtered file, its report and the
    per-camera diagnostic plots (and their .tmp copies). Removed before the filter runs and again
    when it fails, as the cron does, so a camera this run leaves out, or a failed run, can never
    leave last build's files next to this one's."""
    old = [dst, Path(str(dst) + ".tmp"), wdir / "consistency_report.csv",
           wdir / "consistency_report.csv.tmp"] + sorted(wdir.glob("consistency_*.png"))
    for f in old:
        try:
            f.unlink()
        except FileNotFoundError:
            pass


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
    remove_filter_outputs(wdir, dst)       # never leave an old filtered file, report or plot looking current
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
        remove_filter_outputs(wdir, dst)
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
# than 3 camera-days that DEM is the other day alone, so the test cannot tell which day is off (two
# days that differ get equal and opposite offsets). dem_from_contours.py keeps every day itself then
# (and exits 4, 'no points', if a test ever rejects every camera-day); here the option is not passed at
# all with fewer than this many camera-days, so the README can say the test was not run, and why.
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
            f"{DEM_DAY_TEST_MIN} the test cannot tell which day is off (the DEM without one day is the other day "
            f"alone; dem_from_contours.py keeps them all in that case), so the option ({args.dem_max_day_offset} m) "
            f"is not passed. Every camera-day is kept; a difference between them stays in the DEM's spread")


def plan_variant(plan):
    """The calibration variant this plan builds with (--calibration-variant), None for the table's (also for a plan
    that names no cameras)."""
    return next((c.get("variant") for c in (plan.get("cams") or {}).values()
                 if c.get("variant") and c.get("variant") != "table"), None)


def product_plan(prov, pdir):
    """The plan a BUILT product's files are read back with (--summary): its date, its folder and each camera's
    calibration variant from its provenance.json ('table' = none), so a variant's DEM is read as
    dem/<date>_<variant>_* (dem_stem), never as the table calibration's."""
    cams = {}
    for cam, v in (prov.get("cameras") or {}).items():
        var = ((v or {}).get("calibration") or {}).get("variant")
        cams[cam] = {"variant": None if var in (None, "", "table") else var}
    return {"date": prov.get("date") or Path(pdir).name, "out": Path(pdir), "cams": cams}


def dem_stem(plan):
    """The DEM files' stem: dem/<date> with the table's calibration, dem/<date>_<variant> with another (e.g.
    2024-10-23_gcpfree_dem.tif), so a DEM separated from its folder still says which calibration built it."""
    v = plan_variant(plan)
    return plan["out"] / "dem" / (plan["date"] + (f"_{v}" if v else ""))


def map_png(plan, what):
    """A map's file (what: a camera, or 'plan'): maps/<date>_waterlines_<what>.png with the table's calibration,
    maps/<date>_<variant>_waterlines_<what>.png with another (as dem_stem), so a map separated from its folder still
    says which calibration drew it."""
    v = plan_variant(plan)
    return plan["out"] / "maps" / f"{plan['date']}{'_' + v if v else ''}_waterlines_{what}.png"


def calibration_identity(plan, args, blind_note=True):
    """One ASCII line naming this build's calibration, written into every GeoTIFF's description: the table's or
    the variant, each camera's EO file and the first 12 hex of its sha256, and on a BLIND date what the holder's
    check of the DEM can test with it (held_check). A camera left out because the GCP-free procedure did not
    deliver it is named."""
    parts = []
    for cam, c in sorted(plan["cams"].items()):
        sha = file_sha(c["eo_path"])[:12] if c.get("eo_path") else "missing"
        lab = ""
        if plan.get("blind"):
            o_ = calibration_origin(c.get("eo_path") or c.get("eo_file"), eo_dirs(args))
            lab = " " + held_check(o_, plan["blind"], envelope_check(plan, cam, args))[0]
            if lab == " INDEPENDENT" and (o_.get("procedure") or {}).get("pinned"):
                # the label separated from the README: it still points to what the procedure's makers had seen of A
                lab += " (procedure made with stated exposure to A: README EXPOSURE)"
        parts.append(f"{cam} {c.get('eo_file')} sha256 {sha}{lab}")
    for cam in sorted(plan.get("not_delivered") or {}):
        parts.append(f"{cam} NOT DELIVERED by the GCP-free procedure (no data)")
    v = plan_variant(plan)
    text = (f"calibration {('variant ' + v + ' (--eo)') if v else 'of the table (survey_dates.csv)'}: "
            + "; ".join(parts)
            + (f"; BLIND: survey held by {plan['blind']['held_by']}, labels = his check of this DEM"
               if plan.get("blind") and blind_note else ""))
    return text.encode("ascii", "replace").decode("ascii")


def calibration_label(plan, args):
    """What the DEM page subtitle says of the calibration: the day a GCP solution was solved, or for a turned one
    the GCP surveys it rests on and the variant ('Jan 2025 GCPs, turned'), never a date it was not solved on."""
    out = []
    for cam, c in sorted(plan["cams"].items()):
        o = calibration_origin(c.get("eo_path") or c.get("eo_file"), eo_dirs(args))
        if o.get("kind") == "gcp" and o.get("gcp_dates"):
            out.append(o["gcp_dates"][0])
        elif o.get("kind") == "rotated" and o.get("gcp_dates"):
            out.append(f"{'/'.join(o['gcp_dates'])} GCPs, turned")
        elif len(o.get("chain") or []) > 1 or (o.get("procedure") or {}).get("format"):
            out.append("a turned calibration of unknown origin")
        else:                                      # a GCP solution not in CALIBRATION_ORIGINS: its file's date
            out.append(eo_date(c.get("eo_file")) or "?")
    return ", ".join(sorted(set(out)))


# The grids dem_from_contours.py writes that the dem step converts to GeoTIFF; with --dem-loess also the loess
# surface and its error (dem_loess.LOESS_FILES), none of which a figure, comparison or label reads.
DEM_GRIDS = ("dem", "spread", "count", "sigma")
DEM_LOESS_GRIDS = ("loess", "loess_err", "loess_nmse", "loess_skill")


def dem_loess_on(args):
    """--dem-loess (off by default; False for callers whose options predate it)."""
    return bool(getattr(args, "dem_loess", False))


def dem_grids(args):
    """The grids the dem step converts and expects: DEM_GRIDS, + DEM_LOESS_GRIDS with --dem-loess."""
    return DEM_GRIDS + (DEM_LOESS_GRIDS if dem_loess_on(args) else ())


def dem_loess_settings():
    """What the DEM's stamp, provenance.json and README record of a --dem-loess build."""
    from dem_loess import LOESS_SCALES, LOESS_ORDER, LOESS_MAX_NMSE
    return {"script": "dem_loess.py", "half_widths_m": {"alongshore": LOESS_SCALES[0], "cross_shore": LOESS_SCALES[1]},
            "order": LOESS_ORDER, "max_nmse": LOESS_MAX_NMSE, "grids": [f"_{k}" for k in DEM_LOESS_GRIDS],
            "note": "a second surface beside the DEM (quadratic loess with explicit scales, Plant, Holland & Puleo "
                    "2002) and its prediction error; no figure, comparison or label uses it; check its error with "
                    "python3 dem_loess.py evaluate <dem stem> <survey>"}


def dem_cmds(plan, args, day_test=True):
    out = plan["out"]
    stem = dem_stem(plan)
    cmd = [sys.executable, HERE / "dem_from_contours.py", waterline_file(plan), stem,
           "--start-date", plan["first"], "--end-date", plan["last"],
           "--cell", args.dem_cell, "--min-points", args.dem_min_points, "--max-spread", args.dem_max_spread,
           "--no-plot"]
    if args.dem_max_hs:
        cmd += ["--max-hs", args.dem_max_hs]
    if args.dem_max_day_offset and day_test:
        # anchored on the survey date: on a camera whose days drift through the window the days kept
        # are those around the survey, not whichever end weighs more (dem_from_contours.py
        # --anchor-date; 23 Oct 2024, c2: 24 and 26 Oct kept, the survey's own day left out, without it)
        cmd += ["--max-day-offset", args.dem_max_day_offset, "--anchor-date", plan["date"]]
    if dem_loess_on(args):
        cmd += ["--loess"]                  # opt-in: the loess surface and its error beside the DEM (dem_loess.py)
    tif = [sys.executable, HERE / "asc_to_geotiff.py"] + [f"{stem}_{k}.asc" for k in dem_grids(args)] + \
          ["--epsg", args.geotiff_epsg, "--note", calibration_identity(plan, args)]
    station = "/".join(sorted({c["station"] for c in plan["cams"].values()}))
    cal = calibration_label(plan, args)
    var = sorted({c["variant"] for c in plan["cams"].values() if c.get("variant")})   # an --eo build says which
    held_day = (plan.get("blind") or {}).get("survey_date")
    page = [sys.executable, HERE / "dem_figure.py", stem, "--calibration", out / "dem" / "view_calibration",
            "--station", (f"{station} Marconi Beach (calibration {'/'.join(var)}: {cal})" if var else
                          # a blind date's table calibration solved from the held survey's GCPs: said on the page
                          f"{station} Marconi Beach (table calibration {cal}: solved from the held GCPs)"
                          if held_day and held_day in cal.split(", ") else
                          f"{station} Marconi Beach (calibration {cal})")]
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
    outs = [Path(f"{stem}_{k}.{e}") for k in dem_grids(args) for e in ("asc", "tif")] + \
           [Path(f"{stem}_dem.png")]
    src_key = "filter" if src.name.endswith("_filtered.csv") else "merge"
    scripts = STEP_SCRIPTS["dem"] + (["dem_loess.py"] if dem_loess_on(args) else [])
    sig = {"cmd": cmd_text(cmd), "tif": cmd_text(tif), "page": cmd_text(page),
           "after": {src_key: token(out, src_key)}, "scripts": {s: script_sha(s) for s in scripts}}
    fresh, why = is_fresh(out, "dem", sig, outs, [src])
    if fresh and not args.force:
        say("dem", f"SKIPPED: {stem}_dem.asc/.tif/.png are up to date (--force to rebuild)")
        return "skipped"
    say("dem", f"dem_from_contours.py over {plan['first']} .. {plan['last']} ({why}); ~10-40 s on the NUC")
    drop_stamp(out, "dem")
    (out / "dem").mkdir(parents=True, exist_ok=True)
    # no stale grid or page from an earlier build (a --dem-loess build's loess grids included, whatever this one asks)
    for o in outs + [Path(f"{stem}_info.json"), Path(f"{stem}_map.png")] + \
            [Path(f"{stem}_{k}.{e}") for k in DEM_LOESS_GRIDS for e in ("asc", "tif")]:
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
    # the plain north-up map (dem_map.py, <stem>_map.png): what Dr Sherwood found easier to read than the
    # page below; a failure is said, never fatal (the grids are the product)
    if run_cmd([sys.executable, HERE / "dem_map.py", stem, "--title", f"{plan['date']}: {plan['first']} .. "
                f"{plan['last']}, cameras {', '.join(sorted(plan['cams']))}"], log) != 0:
        warn(f"the plain DEM map (dem_map.py) failed; the grids are fine. See {log}")
    rc = run_cmd(page, log)
    extra = {"cells": cells, "source": str(src),
             "settings": dict(dem_settings_requested(args), day_test_off=no_day_test or None,
                              **({"loess": dem_loess_settings()} if dem_loess_on(args) else {}))}
    if rc != 0 or not Path(f"{stem}_dem.png").exists():
        warn(f"the DEM page (dem_figure.py) failed (exit {rc}); the grids are fine. See {log}")
        write_stamp(out, "dem", dict(sig, page_failed=True), [cmd, tif, page], scripts=scripts, extra=extra)
        return "built (page failed)"
    write_stamp(out, "dem", sig, [cmd, tif, page], scripts=scripts, extra=extra)
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
        if ln.startswith("No georectified points"):
            return (f"no waterline points in the window (dem_from_contours.py, exit 4: "
                    f"'{ln.strip()[:60]}')")
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
        png = map_png(plan, cam)
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
    png = map_png(plan, "plan")
    surveys = [s for s in plan["surveys"] if resolve_survey(s, plan, args)[0]]
    # every file of a survey (the GCPs come as one file per camera): all are drawn, so all sign it
    sig = {"src": str(src), "after": {"waterlines": waterline_token(plan)}, "v": 11,
           "surveys": [[s["name"]] + [file_sig(p) for p in resolve_survey(s, plan, args)[0]] for s in surveys],
           "scripts": scripts_sig("map_plan")}
    if plan.get("blind"):
        sig["blind"] = plan["blind"]["date"]          # drawn over no survey, and titled so
    if plan.get("blind") or plan_variant(plan):
        sig["calibration"] = calibration_identity(plan, args, blind_note=False)   # in its title
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
        except SystemExit as exc:                 # a loader's sys.exit (survey_compare.py's readers)
            if is_termination(exc):
                raise
            warn(f"plan-view map failed: {exc}")
    missing = [cam for cam in plan["cams"] if not map_png(plan, cam).exists()]
    if missing:
        warn("waterline map(s) on the photos NOT drawn: "
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


def plan_survey_points(s, plan, args, folder):
    """The points of point survey s for the plan view (survey_compare.read_points' dict): every file
    of it (the GCPs come as one file per camera), merged as the comparison merges them (a target in
    two files once, into folder/<date>_<name>_plan_points.csv), not only the first file found."""
    paths = resolve_survey(s, plan, args)[0]
    if len(paths) == 1:
        return read_survey_points(paths[0])
    return read_survey_points(merged_points(paths, Path(folder) / f"{plan['date']}_{s['name']}_plan_points.csv"))


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
    p = None
    if pts:
        try:
            p = plan_survey_points(pts, plan, args, Path(png).parent)
        except SurveyFileError as exc:
            warn(f"plan view: the points of {pts['name']} are NOT drawn: {exc}")
            notes.append(f"{pts['name']}: points NOT drawn (survey file unreadable, see the build log)")
    if p is not None:
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
    if plan.get("blind"):
        notes = [f"BLIND: no survey drawn (the survey is held by {plan['blind']['held_by']} for an independent check)"]
    if plan.get("blind") or plan_variant(plan):
        # the calibration that placed these lines (A and B of a blind date differ; their maps must not be mixed up)
        import textwrap
        notes = textwrap.wrap(calibration_identity(plan, args, blind_note=False), 100) + notes
    ax.set_title(f"{plan['date']}: filtered waterlines in plan view\n" + ("\n".join(notes) or "no survey drawn"),
                 loc="left", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(png, dpi=130, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


# ---------------------------------------------------------------------
# Step: compare, and the label
# ---------------------------------------------------------------------

def resolve_survey(s, plan, args):
    """-> (list of Paths, note). An empty path is looked for by the survey date (GCP files arrive later).
    A survey of a BLIND date (BLIND_DATES: the row's date or survey date, or a blind plan) is never looked
    for: BlindSurvey, before any folder is searched."""
    b = held_survey(s) or (plan or {}).get("blind")
    if b:
        raise BlindSurvey(f"{s.get('name')}: {blind_statement(b)} ({b['survey']}): survey_products.py never looks "
                          f"for it (BLIND_DATES {b['date']})")
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
        expect = (f"the file {s['path']} named in {args.surveys} (looked for by its name in "
                  f"{', '.join(args.survey_dirs)}, never in the folder the command runs from; a path with a "
                  f"folder in it as given first; --survey-dirs adds folders), or correct its path there")
    elif s["survey_type"] == "points":
        # a folder this run looks in: the station's Chelsea folder when it is one of them (the
        # default), else the first --survey-dirs folder (a file put in a folder not searched would
        # never be found)
        dirs = [str(d) for d in (args.survey_dirs or [])]
        where = str(CHELSEA) if (not dirs or str(CHELSEA) in dirs) else dirs[0]
        others = [d for d in dirs if d != where]
        expect = (f"{where}/{s.get('survey_date') or plan['date']}_Marconi_Extrinsic_Targets_<cam>_xyz.csv "
                  + (f"(or in {', '.join(others)}: the folders this run searches, --survey-dirs) " if others else "")
                  + f"(UTM 19N / NAVD88; num,E,N,Z rows without a header like "
                  f"calibration/2025-11-13_Marconi_Extrinsic_Targets_c1_xyz.csv, any target names, or an Emlid "
                  f"Flow CSV; one file per camera is fine), or write its path in {args.surveys} (column path)")
    else:
        expect = f"the survey file (lidar DSM .tif or .asc) and write its name in {args.surveys} (column path)"
    # maps too: the plan view draws the survey (a compare-only run left it drawn without the points)
    return (f"SKIPPED {s['name']}: {note}. Provide {expect}, then rerun: python3 survey_products.py "
            f"--date {plan['date']} --steps maps,compare" + (" --force-disabled" if disabled else "")
            + " (the photo maps are kept if current; the plan view and the comparison are made, ~1 min)")


class SurveyFileError(Exception):
    """A survey file survey_compare.read_points cannot read (empty, no coordinate columns, ...)."""


def is_termination(exc):
    """True for the SystemExit main()'s SIGTERM handler raises: a kill must end the run, never be
    taken for a loader's sys.exit on a bad file."""
    return str(exc).startswith("terminated by signal")


def read_survey_points(path):
    """survey_compare.read_points, with its sys.exit on a file it cannot read turned into
    SurveyFileError: a bad survey file costs that survey's comparison or the plan view's points, not
    the build (a SystemExit escaped step_maps' 'except Exception' and failed the product)."""
    from survey_compare import read_points
    try:
        return read_points(path)
    except SystemExit as exc:
        if is_termination(exc):
            raise
        raise SurveyFileError(str(exc) or f"{path}: not a point survey survey_compare.py can read")


def merged_points(paths, dst):
    """Several files of one point survey (the GCPs come as one file per camera) -> one file that
    survey_compare.py reads. Each file is read as survey_compare.py reads it (CIRN num,E,N,Z rows with
    or without numeric target names, Emlid Flow CSV, a CSV with E/N/Z headers; CRLF or not) and the
    points are merged on their parsed coordinates (a target in two files, to the mm, once): copying
    the first four columns of every row merged two Emlid files into a header-less mess survey_compare
    .py could not read, and dropped a second file whose columns came in another order. Written as
    CIRN rows (name,E,N,Z) when every file is CIRN, else with Name,Code,Easting,Northing,Elevation,
    Description,CRS headers (code, description and CRS kept; compare_rtk.py's transect reader reads
    it too). The file is rewritten only when its content changes (its age keeps the comparison
    current). Raises SurveyFileError naming a file that cannot be read. -> dst"""
    parts = [(p_, read_survey_points(p_)) for p_ in paths]
    cirn = all(q["format"].startswith("CIRN") for _, q in parts)
    seen, rows = set(), []
    for _, q in parts:
        for k in range(len(q["E"])):
            key = tuple(int(round(float(q[c][k]) * 1000)) for c in ("E", "N", "Z"))
            if key in seen:
                continue
            seen.add(key)
            xyz = [f"{float(q[c][k]):.4f}".rstrip("0").rstrip(".") for c in ("E", "N", "Z")]
            rows.append([q["name"][k]] + xyz if cirn else
                        [q["name"][k], q["code"][k], xyz[0], xyz[1], xyz[2], q["desc"][k], q["crs"][k]])
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    if not cirn:
        w.writerow(["Name", "Code", "Easting", "Northing", "Elevation", "Description", "CRS"])
    w.writerows(rows)
    data = buf.getvalue().encode()
    # compared as bytes: a text read turns csv's default '\r\n' into '\n', so the file never
    # matched, was rewritten on every run, and its new time re-ran the comparison every time
    if not (Path(dst).exists() and Path(dst).read_bytes() == data):   # keep its age if unchanged
        Path(dst).write_bytes(data)
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
        hits = resolve_file(given, getattr(args, "survey_dirs", []) or [], typed=True)
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


def setup_fit_known(coef):
    """True when C is a key of SETUP_FITS (where it was fitted is known)."""
    return any(abs(k - float(coef or 0.0)) < 1e-9 for k in SETUP_FITS)


def replaced_coefs(plan, args, steps):
    """The C(s) of the waterlines on disk, sorted, when this run would rebuild them with another C in
    the same folder (e.g. [0.037]), else None."""
    if "detect" not in steps or not plan.get("out") or not read_stamp(plan["out"], "merge"):
        return None
    fits, _ = setup_in_use(plan, args)
    built = sorted({round(float(f["coef"]), 6) for f in fits if f.get("coef") is not None})
    if not built or built == [round(float(args.setup_coef or 0.0), 6)]:
        return None
    return built


def root_for_c(root, coef, held=None):
    """The output root a build with this C belongs in: the C = SETUP_COEF products in the base root
    (the default /mnt/I2Rgus_Data/survey_products), any other C in <base root>_c<C> (e.g. ..._c0).
    The base is this root without a '_c<number>' suffix, so the advice from a ..._c0 root points
    back at the station-C root, never at ..._c0_c0.037.
    held: the C(s) this root's waterlines carry. A name says only which C a root is MEANT for (a date's
    first build is never refused: a root without suffix may hold C = 0, a ..._c0 root C = 0.037). When
    the name gives this very root to coef although it holds another C, the root given is one beside
    it, <root>_c<C>, never the root that refused the run; root_name_note says what this root holds."""
    root = Path(root)
    base = re.sub(r"_c\d+(\.\d+)?$", "", root.name) or root.name
    want = round(float(coef or 0.0), 6)
    alt = root.parent / (base if abs(want - SETUP_COEF) < 1e-9 else f"{base}_c{want:g}")
    if alt == root and any(abs(round(float(h), 6) - want) > 1e-9 for h in (held or ())):
        alt = root.parent / f"{root.name}_c{want:g}"
    return alt


def root_name_note(root, coef, held, how=None):
    """When the name of this root gives it to C = coef but its waterlines carry another C: what it holds,
    where that C belongs and how to free the root for coef. -> text, or '' when the name does not
    give this root to coef or the root holds coef. how: how it came to hold that C (default: a first
    build, which is never refused)."""
    root = Path(root)
    want = round(float(coef or 0.0), 6)
    held = sorted({round(float(h), 6) for h in (held or ())})
    if not held or held == [want] or root_for_c(root, want) != root:
        return ""
    hs = " / ".join(f"{h:g}" for h in held)
    if how is None:
        how = "a date's first build in a root is never refused, whatever its C"
        if root_for_c(root, SETUP_COEF) == root:          # a root meant for the station's C
            how += f", e.g. --setup-coef {held[0]:g} without --output-root builds into the default root"
    text = f"the name of {root} gives it to C = {want:g}, but it holds the C = {hs} product ({how})"
    if len(held) == 1:
        home = root_for_c(root, held[0])
        text += (f". To keep this root for C = {want:g}: build C = {hs} again in its own root (--setup-coef {hs} "
                 f"--output-root {home}), then rerun here with --replace")
    return text


def setup_replace_note(plan, args, steps):
    """A run that would build a C other than the one the waterlines on disk carry, into the SAME
    folder: what it would replace, what reruns, and the root that C belongs in (never this root:
    root_for_c with what it holds). -> text or None.
    Such a run is refused unless it passes --replace (build_date)."""
    built = replaced_coefs(plan, args, steps)
    if built is None:
        return None
    old = " / ".join(f"{c:g}" for c in built)
    want = round(float(args.setup_coef or 0.0), 6)
    root = Path(args.output_root)
    alt = root_for_c(root, want, built)
    text = (f"the C = {old} product in {plan['out']} would be REPLACED: this run builds C = {want:g} into the same "
            f"folder (waterlines, filter, DEM, maps and comparisons rebuilt, the C = {old} ones overwritten, and "
            f"--summary on {root} then shows only the C = {want:g} build). Build each C in its own root: "
            f"--output-root {alt} for C = {want:g} (then --summary --output-root on each root); give --replace "
            f"only to replace this product on purpose.")
    misnamed = root_name_note(root, want, built)
    if misnamed:
        text += f" NOTE: {misnamed}."
    if plan.get("era") != "live":
        text += (f" For this date the photo detection runs again (~{SEC_PER_PHOTO_DETECT:.0f} s a photo, ~50 min "
                 f"a two-camera week): each line is placed at the elevation its setup gives it.")
    return text


def window_replace_note(plan, args, steps):
    """A run that rebuilds a built date with another photo window or UTC hours, in the SAME folder:
    what it replaces and where to build it instead. -> text or None (only when this run rebuilds from
    the forcing through the DEM: otherwise keep_built_window keeps the built window)."""
    if not set(WINDOW_STEPS) <= set(steps):
        return None
    diff = window_differences(plan, built_settings(plan))
    if not diff:
        return None
    root = Path(args.output_root)
    return (f"replacing the product in {plan['out']} built for another window/hours ({'; '.join(diff)}): this "
            f"run rebuilds the date for its own window into the same folder, and --summary on {root} then shows "
            f"only this build. To keep both, build it in its own --output-root, e.g. "
            f"{root.parent / (root.name + '_recheck')}.")


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
        # each reason once: a survey of several files (one GCP file per camera) is checked file by
        # file, and a rule about the cameras' calibrations fires once per file
        if heuristic and self.override and LABELS.index(new) > LABELS.index(self.label):
            if f"{new}: {why}" not in self.overridden:
                self.overridden.append(f"{new}: {why}")
            return
        if f"{new}: {why}" not in self.reasons:
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
            elif named and s["survey_type"] == "points" and re.search(r"target|gcp", sname, re.I):
                # calibration targets the calibration names (solved from them?): PARTLY-CIRCULAR, as the
                # GCP-day rule below and survey_compare.audit_label say
                v.worse("PARTLY-CIRCULAR", f"the {cam} calibration {Path(eo).name} names these calibration targets "
                                           f"in its notes (solved from them?)", heuristic=True)
            elif named:
                v.worse("CIRCULAR", f"the {cam} calibration {Path(eo).name} names this survey in its notes "
                                    f"(was it fitted to it?)", heuristic=True)
            elif lin["survey_fit"] and lin["unresolved"]:
                v.worse("CIRCULAR", f"the {cam} pointing {Path(eo).name} is survey-derived{via} and which survey "
                                    f"cannot be resolved ({'; '.join(lin['unresolved'])}): it may be this one",
                        heuristic=True)
        # any point survey made on the day of the GCPs the calibration rests on may be those GCPs, whatever
        # the file is called; a calibration turned from another rests on that one's GCPs (CALIBRATION_ORIGINS),
        # whatever date its own name carries
        if s["survey_type"] == "points":
            org = calibration_origin(c.get("eo_path") or c["eo_file"], eo_dirs(args))
            days = org.get("gcp_dates") if org.get("kind") in ("gcp", "rotated") else [eo_date(c["eo_file"])]
            if s.get("survey_date") in (days or []):
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
        # the setup fits are applied here by survey_products' own rule (setup_caveats, Verdict), with
        # surveys.csv's label_override_reason; survey_compare.py's run applies them again on its own
        return audit_label(label, survey_path, s["survey_type"], s.get("survey_date"), cams,
                           envs[0] if len(envs) == 1 else None, dem, dirs=eo_dirs(args) if args else (),
                           setup_fits=False)
    except SystemExit as exc:                      # its loaders exit on a missing file
        return [(None, f"survey_compare.py's checks could not run: {exc}", "note")]


def frames_used_dates(plan):
    """(first, last) day of the frames the DEM gridded (dem_from_contours.py's info), else of the
    waterlines in the window, else the window."""
    info = load_json(Path(f"{dem_stem(plan)}_info.json")) or {}
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
    stem = dem_stem(plan)
    dem = Path(f"{stem}_dem.asc")
    if plan.get("blind"):
        b = plan["blind"]
        say("compare", f"NOT RUN: {blind_statement(b)} ({b['survey']}); this is a DEM-only product and no survey "
                       f"of {plan['date']} is looked for (BLIND_DATES)")
        state["comparisons"] = {}
        return "not run (blind: DEM only)"
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
        misnamed = root_name_note(args.output_root, args.setup_coef, used)
        warn(f"this run's --setup-coef {args.setup_coef} is NOT what the waterlines were built with "
             f"(C = {', '.join(str(c) for c in sorted(used))}): the label, README and provenance use the C that "
             f"built them. To apply {args.setup_coef}, build it from the detection (all steps) in its own "
             f"--output-root ({root_for_c(args.output_root, args.setup_coef, used)}): in this root it would replace "
             f"this product (refused without --replace)" + (f". NOTE: {misnamed}" if misnamed else ""))
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
            try:
                survey_path = merged_points(paths, cdir / f"{name}_survey_points.csv")
            except SurveyFileError as exc:
                warn(f"comparison {name} FAILED: {exc}")
                results[s["name"]] = {"status": "failed", "label": s["label"], "reason": str(exc)}
                state["failures"].append(f"compare {s['name']}: a survey file cannot be read ({exc})")
                for p in comparison_files(plan, s):     # nothing of an earlier file may look current
                    if p.exists():
                        p.unlink()
                drop_stamp(out, f"compare_{s['name']}")
                failed_any += 1
                continue
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
        eod = pointing_eo_day_specs(plan, [e.split("=", 1)[0] for e in eos])
        if eod:
            cmd += ["--camera-eo-day"] + eod
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
           "compare_dem_survey.py", "compare_rtk.py", "marconi_water_level.py", "timex_timing.py",
           "dem_loess.py", "setup_model.py", "pointing_daily.py"]

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


def calibration_record(c):
    """What identifies a camera's calibration in a stamp: file name, path, content, variant."""
    return {"eo_file": c.get("eo_file"), "eo_path": str(c["eo_path"]) if c.get("eo_path") else None,
            "eo_sha256": file_sha(c["eo_path"]) if c.get("eo_path") else None, "variant": c.get("variant")}


def eo_comment_lines(path):
    """The comment lines of an EO file, '#' and spaces stripped, in order ([] when it cannot be read)."""
    try:
        text = Path(path).read_text(errors="replace")
    except OSError:
        return []
    return [ln.strip().lstrip("#").strip() for ln in text.splitlines() if ln.strip().startswith("#")]


def calibration_acceptance(lines):
    """The pre-registered acceptance a calibration file states, from its comment lines: an ACCEPTANCE_LINE
    ('acceptance: PASS ...'), or gcpfree_eo.py's DELIVERY_LINE ('CALIBRATION B ...: DELIVERED ...', its words
    kept in the text). -> {"verdict": 'PASS' | 'PASS-WIDER' | 'FAIL', "text": what follows} or None when it
    states none."""
    for ln in lines:
        m = ACCEPTANCE_LINE.match(ln)
        said = None
        if not m:
            m = DELIVERY_LINE.match(ln)
            said = m.group(1).upper() if m else None
        if m:
            t = m.group(2).strip()
            if t.startswith("(") and t.endswith(")"):
                t = t[1:-1].strip()
            if said:
                return {"verdict": DELIVERY_VERDICT[said], "text": f"gcpfree_eo.py: {said}" + (f"; {t}" if t else "")}
            return {"verdict": m.group(1).upper(), "text": t}
    return None


def rotation_claims(notes, own_name):
    """What a turned calibration's notes claim, reference by reference in the order they are named: the turn
    ('X.yaml rotated by N deg', estimate_eo_rotation.py), the pan change ('pan change N deg from X.yaml') and
    the reference's sha256 ('Reference EO: X.yaml ... sha256 H', gcpfree_eo.py); the reference whose camera
    position the file keeps ('Position (x, y, z) kept from X.yaml', else the first named); how the turn was
    measured ('Rotation source: ...'). -> {"refs": [{"name", "angle", "pan", "sha256"}], "position_ref",
    "source"}."""
    found = [(m.start(), Path(m.group(1)).name, "angle", m.group(2)) for m in ROTATED_REF.finditer(notes)]
    found += [(m.start(), Path(m.group(2)).name, "pan", m.group(1)) for m in PAN_FROM.finditer(notes)]
    found += [(m.start(), Path(m.group(1)).name, "sha256", m.group(2).lower()) for m in REFERENCE_EO.finditer(notes)]
    pm = POSITION_FROM.search(notes)
    if pm:
        found.append((pm.start(), Path(pm.group(1)).name, None, None))
    refs = {}
    for _, n, k, v in sorted(found, key=lambda t: t[0]):
        if n == own_name:
            continue
        r = refs.setdefault(n, {"name": n, "angle": None, "pan": None, "sha256": None})
        if k and r[k] is None:
            r[k] = v
    order = list(refs.values())
    pos = Path(pm.group(1)).name if pm else (order[0]["name"] if order else None)
    sm = ROTATION_SOURCE.search(notes)
    return {"refs": order, "position_ref": pos if pos != own_name else None,
            "source": sm.group(1).strip() if sm else None}


def git_file_state(path):
    """Is `path` committed to git and unchanged since? -> {"committed": bool, "commit", "date" (the last commit that
    touched it), "commits": how many commits touched it, "history": [[commit, date], ...] newest first, "toplevel":
    the repository it is in, "why"}."""
    q = Path(os.path.realpath(str(path)))

    def git(*a):
        return subprocess.run(["git", "-C", str(q.parent)] + list(a), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True)
    try:
        top = git("rev-parse", "--show-toplevel")
        lg = git("log", "--format=%H %cI", "--", q.name)
        st = git("status", "--porcelain", "--", q.name)
    except OSError as exc:
        return {"committed": False, "commits": 0, "why": f"git cannot be run ({exc})"}
    hist = [ln.split(" ", 1) for ln in lg.stdout.splitlines() if ln.strip()] if lg.returncode == 0 else []
    out = {"commits": len(hist), "history": hist,
           "toplevel": os.path.realpath(top.stdout.strip()) if top.returncode == 0 and top.stdout.strip() else None}
    if not hist:
        return dict(out, committed=False, why="NOT committed to git")
    out.update(commit=hist[0][0], date=hist[0][1] if len(hist[0]) > 1 else None)
    if st.returncode != 0 or st.stdout.strip():
        return dict(out, committed=False, why="changed since its commit")
    return dict(out, committed=True, why="")


_REPO_TOP = {}


def repo_toplevel():
    """The git repository this script is in (its top folder, real path), or None."""
    if "top" not in _REPO_TOP:
        try:
            r = subprocess.run(["git", "-C", str(HERE), "rev-parse", "--show-toplevel"], stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, universal_newlines=True)
            _REPO_TOP["top"] = os.path.realpath(r.stdout.strip()) if r.returncode == 0 and r.stdout.strip() else None
        except OSError:
            _REPO_TOP["top"] = None
    return _REPO_TOP["top"]


def gcpfree_pin_of(sha):
    """The GCPFREE_PINS entry whose gcpfree_eo.py has this sha256, else None."""
    return next((p for p in GCPFREE_PINS if p["gcpfree_eo.py"] == str(sha or "").lower()), None)


def gcpfree_delivery_runs(runs_dir):
    """The finished DELIVERY runs (eligible) in a gcpfree runs folder, oldest first (as gcpfree_eo.py orders them):
    [(run.json dict, folder)]."""
    found = []
    for d in sorted(q for q in Path(runs_dir).glob("*") if q.is_dir()) if Path(runs_dir).is_dir() else []:
        js = load_json(d / "run.json") or {}
        if js.get("status") == "finished" and js.get("eligible"):
            found.append((js, d))
    return sorted(found, key=lambda t: (t[0].get("started_utc") or t[0].get("run_id") or "", t[1].name))


def gcpfree_run_records(runs_dir):
    """The run records (runs/<id>/run.json) of a gcpfree runs folder that are not as first committed: committed
    and gone now (deleted in a later commit, or from the working tree), or committed and then changed (in a later
    commit, or in the working tree). Every run's record is kept as written: an earlier run removed, or rewritten
    (e.g. as a trial), would let a later run pass for the first. -> {"removed": [run id], "changed": [run id]}
    (git not available, or the folder not in git: both empty)."""
    rdir = Path(str(runs_dir))
    out = {"removed": set(), "changed": set()}
    if not rdir.is_dir():
        return {k: [] for k in out}

    def git(*a):
        return subprocess.run(["git", "-C", str(rdir)] + list(a), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True)
    try:
        lg = git("log", "--format=", "--name-status", "--", ".")
        st = git("status", "--porcelain", "--", ".")
    except OSError:
        return {k: [] for k in out}
    rows = [ln.split("\t") for ln in lg.stdout.splitlines() if "\t" in ln] if lg.returncode == 0 else []
    rows += [[ln[:2].strip(), ln[3:]] for ln in st.stdout.splitlines() if len(ln) > 3] if st.returncode == 0 else []
    for row in rows:
        code, q = row[0][:1], Path(row[-1].strip().strip('"'))
        if q.name != "run.json" or not q.parent.name:
            continue
        if code == "D":
            out["removed"].add(q.parent.name)
        elif code in ("M", "R", "T"):
            out["changed"].add(q.parent.name)
    out["removed"] = {r for r in out["removed"] if not (rdir / r / "run.json").is_file()}
    return {k: sorted(v) for k, v in out.items()}


GCPFREE_RECORD = re.compile(r"(^|/)(PREREGISTERED\.txt|runs/[^/]+/run\.json)$")


def gcpfree_foreign_records(gdir):
    """Commits OUTSIDE this branch's history (another branch, a stash, the reflog: a branch made from before the
    pre-registration, a reset) that wrote a PREREGISTERED.txt or a runs/<id>/run.json in the gcpfree folder whose
    content this branch's history does not hold (a rebased copy holds the same content: not counted). The first run,
    or an earlier pre-registration, kept there would be invisible from this branch (gcpfree_eo.foreign_records, the
    same rule). -> [(commit, path)] ([] when git cannot tell)."""
    anc = Path(str(gdir))
    while not anc.exists() and anc != anc.parent:
        anc = anc.parent
    spec = os.path.relpath(str(Path(str(gdir))), str(anc))

    def blobs(*refs):
        r = subprocess.run(["git", "-C", str(anc), "log"] + list(refs) + ["--no-abbrev", "--format=@%H", "--raw",
                                                                           "--", spec],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        out, commit = [], None
        for ln in r.stdout.splitlines() if r.returncode == 0 else []:
            if ln.startswith("@"):
                commit = ln[1:].strip()
            elif ln.startswith(":") and "\t" in ln:
                meta, path = ln.split("\t", 1)
                f = meta.split()
                if len(f) >= 4 and GCPFREE_RECORD.search(path.strip()) and set(f[3]) != {"0"}:
                    out.append((commit, path.strip(), f[3]))
        return out, r.returncode
    try:
        head, _ = blobs("HEAD")
        every, rc = blobs("--all", "--reflog")
    except OSError:
        return []
    have = {b for _, _, b in head}
    found = []
    for c, path, b in every if rc == 0 else []:
        if b not in have and (c, path) not in found:
            found.append((c, path))
    return found


def gcpfree_uncommitted_runs(runs_dir):
    """The run folders of a gcpfree runs folder whose record (run.json) is not committed as it is (untracked, or
    changed), or that have none: every run, a trial or one cut off, is committed (its record is what B's README lists).
    -> [run id]."""
    rdir = Path(str(runs_dir))
    if not rdir.is_dir():
        return []
    try:
        st = subprocess.run(["git", "-C", str(rdir), "status", "--porcelain", "--untracked-files=all", "--", "."],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    except OSError:
        return []
    if st.returncode != 0:
        return []
    top = subprocess.run(["git", "-C", str(rdir), "rev-parse", "--show-toplevel"], stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, universal_newlines=True).stdout.strip()
    bad = set()
    for ln in st.stdout.splitlines():
        q = Path(top) / ln[3:].strip().strip('"')
        if q.name == "run.json" and q.parent.parent == Path(os.path.realpath(str(rdir))):
            bad.add(q.parent.name)
    bad |= {d.name for d in rdir.iterdir() if d.is_dir() and not (d / "run.json").is_file()}
    return sorted(bad)


def prereg_photos(prereg):
    """The photos a PREREGISTERED.txt pins ('photo: <cam> <window> <day> <name> <size> <sha256>' lines, as gcpfree_eo.py
    reads them) and its 'photos_sha256:' line. -> (rows, sha256 line or None)."""
    rows, sha = [], None
    try:
        for ln in Path(str(prereg)).read_text().splitlines():
            s = ln.strip()
            if s.startswith("photos_sha256:"):
                sha = s.split(":", 1)[1].strip()
            elif s.startswith("photo:"):
                w = s.split(":", 1)[1].split()
                rows.append(w[:4] + [int(w[4]) if w[4].isdigit() else w[4]] + w[5:] if len(w) == 6 else w)
    except OSError:
        pass
    return rows, sha


def photos_sha256(rows):
    """gcpfree_eo.photos_sha256: the sha256 of a list of photo rows."""
    return hashlib.sha256(json.dumps([list(r) for r in rows]).encode()).hexdigest()


def prereg_rule_history(prereg):
    """The outcome rule of every commit of a PREREGISTERED.txt, oldest first: [{"commit", "date", "rule"}] ([] when
    it is not in git). Only that file is read (git show)."""
    q = Path(os.path.realpath(str(prereg)))

    def git(*a):
        return subprocess.run(["git", "-C", str(q.parent)] + list(a), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True)
    try:
        lg = git("log", "--reverse", "--format=%H %cI", "--", q.name)
    except OSError:
        return []
    out = []
    for ln in lg.stdout.splitlines() if lg.returncode == 0 else []:
        h, _, when = ln.strip().partition(" ")
        if not h:
            continue
        sh = git("show", f"{h}:./{q.name}")
        rule = next((s.split(":", 1)[1].strip() for s in sh.stdout.splitlines()
                     if s.strip().startswith("outcome_if_one_camera:")), None) if sh.returncode == 0 else None
        out.append({"commit": h, "date": when, "rule": rule})
    return out


def angle_between_deg(e1, e2):
    """The rotation (deg) between the pointings of two EOs."""
    from georectify import cirn_angles_to_R
    R1, R2 = cirn_angles_to_R(*e1[3:6]), cirn_angles_to_R(*e2[3:6])
    return float(np.degrees(np.arccos(np.clip((np.trace(R1 @ R2.T) - 1.0) / 2.0, -1.0, 1.0))))


def procedure_check(eo_path, dirs=()):
    """
    Was this calibration made by the PRE-REGISTERED GCP-free procedure (GCPFREE_PINS), and is it what that
    procedure DELIVERED here? Reads the file (its real path: a link is followed to its run folder), its comment
    lines and numbers, the report and run.json written beside it, its references (the report's; looked for beside
    it, in `dirs`, in the calibration folder its run folder sits in, and in the repository's), this repository's
    gcpfree folder (GCPFREE_DIR: the pre-registration and every run in runs/) and git. Every one of these must hold
    (the GCPFREE_PINS comment lists them), else it is 'NOT the pre-registered procedure' and `why` says what failed:
    made by a pinned gcpfree_eo.py with its modules and constants; the report beside it agrees (numbers, verdict,
    delivery, sha256); its 'rotated by', 'pan change' and 'Position ... kept from' lines are true of its numbers and
    of its references (whose content is the one the report recorded); its run folder is GCPFREE_DIR/runs/<id>, made
    under GCPFREE_DIR/PREREGISTERED.txt as it is now (same sha256, no commit since the run), from the station's photo
    folders and the photos it pins, as the first finished delivery run there (or the one allowed second run, with Aug
    photos); every run record there committed; nothing outside this branch's history (gcpfree_foreign_records); its
    file, report, run.json and run_report.txt committed in this repository, unchanged, each touched by ONE commit
    (the same one). What it cannot see: the GCPFREE_PINS comment's list (a B edited BEFORE that commit with its report
    and run.json rewritten to agree; a copy of the procedure run elsewhere before the pre-registration; a re-clone; a
    branch removed with its reflog): only the sha256s sent to the holder beforehand show those.
    -> {"format": it says it is gcpfree_eo.py's, "pinned": every check passed, "why": [what failed], and from its
    report: camera, run_id, run_no, variant, first_run, verdict, delivery, source, solved_in_run, outcome_rule,
    preregistration {sha256, commit, commits, decided_by, path, git}, roots, station_roots, aug_dir, uncertainty,
    exposure, run_started_utc, run_finished_utc, report, report_sha256, b_sha256, run_dir, committed (git state of
    the B file: commit, date, commits)}.
    """
    p = Path(os.path.realpath(str(eo_path)))
    out = {"format": False, "pinned": False, "why": []}
    lines = eo_comment_lines(p)
    notes = " ".join(lines)
    made = MADE_BY.search(notes)
    said = next((m.group(1) for m in (DELIVERY_LINE.match(ln) for ln in lines) if m), None)
    if not (made or said or "gcpfree_eo.py" in notes):
        out["why"].append("not made by gcpfree_eo.py (no 'Made by gcpfree_eo.py (sha256 ...)' line): turned by hand "
                          "or by another procedure")
        return out
    out["format"] = True
    why = out["why"]
    if not made:
        why.append("no 'Made by gcpfree_eo.py (sha256 ...); report ....json' line")
        return out
    pin = gcpfree_pin_of(made.group(1))
    if pin is None:
        why.append(f"made by a gcpfree_eo.py (sha256 {made.group(1)[:12]}...) that is not the pre-registered one "
                   f"({', '.join(x['gcpfree_eo.py'][:12] + '...' for x in GCPFREE_PINS)})")
    rd = p.parent
    rp = rd / Path(made.group(2)).name
    out["report"] = str(rp)
    try:
        rep = json.loads(rp.read_text())
    except (OSError, ValueError) as exc:
        why.append(f"its report {rp.name} is not beside it or cannot be read ({exc.__class__.__name__})")
        return out
    run = rep.get("run") or {}
    pr = run.get("preregistration") or {}
    rj = load_json(rd / "run.json")
    sha = sha256(p)                        # read now (not file_sha's cache): the file as it is
    out.update(camera=rep.get("camera"), run_id=run.get("run_id"), run_no=run.get("run_no"), variant=run.get("variant"),
               first_run=run.get("first_run"), verdict=rep.get("verdict"), delivery=rep.get("delivery"),
               source=rep.get("source"), solved_in_run=rep.get("solved_in_run"),
               outcome_rule=pr.get("outcome_if_one_camera"),
               preregistration={"sha256": pr.get("sha256"), "commit": pr.get("commit"),
                                "commits": len(pr.get("git_history") or []), "decided_by": pr.get("decided_by"),
                                "decided_on": pr.get("decided_on"), "path": pr.get("path"),
                                "commit_date": ((pr.get("git_history") or [{}])[0] or {}).get("date")},
               roots=run.get("roots"), aug_dir=run.get("aug_dir"),
               station_roots=[str(r) for r in (run.get("roots") or [])] == GCPFREE_STATION_ROOTS,
               uncertainty=rep.get("uncertainty"), exposure=rep.get("exposure"), report_sha256=sha256(rp),
               b_sha256=sha, run_dir=str(rd), procedure_sha256=made.group(1).lower(),
               run_started_utc=(rj or {}).get("started_utc") or run.get("started_utc"),
               run_finished_utc=(rj or {}).get("finished_utc"),
               preregistration_here=str(Path(os.path.realpath(str(GCPFREE_DIR / "PREREGISTERED.txt")))))
    # the report beside it: written by this run, for this file
    if rep.get("b_file_sha256") != sha:
        why.append("its report records another sha256 for the B file (changed after the run)")
    if rep.get("delivery") != "DELIVERED":
        why.append(f"its report says it was not delivered ({rep.get('delivery')})")
    if said != rep.get("verdict"):
        why.append(f"its verdict line ({said}) is not its report's ({rep.get('verdict')})")
    if run.get("procedure_sha256") != made.group(1).lower():
        why.append("its report was written by another gcpfree_eo.py than its 'Made by' line names")
    if pin is not None and (run.get("modules_sha256") or {}) != pin["modules"]:
        bad = sorted(k for k in set(pin["modules"]) | set(run.get("modules_sha256") or {})
                     if (run.get("modules_sha256") or {}).get(k) != pin["modules"].get(k))
        why.append(f"its run used other modules than the pre-registered ones ({', '.join(bad)})")
    if pin is not None and run.get("constants_sha256") != pin.get("constants_sha256"):
        why.append(f"its run's pre-registered constants are not the pinned ones (sha256 "
                   f"{str(run.get('constants_sha256'))[:12]}..., pinned {str(pin.get('constants_sha256'))[:12]}...): "
                   f"changed in the running process (e.g. other photo folders patched in)")
    if not run.get("eligible"):
        why.append("a TRIAL run of the procedure, not a delivery run"
                   + (f" ({'; '.join(run['trial_why'])})" if run.get("trial_why") else ""))
    if not out["station_roots"]:
        why.append(f"its photo folders are not the station's ({', '.join(str(r) for r in run.get('roots') or [])}): "
                   f"a test of the procedure, not the delivered B")
    if [str(r) for r in (pr.get("roots") or [])] != GCPFREE_STATION_ROOTS:
        why.append("its pre-registration names other photo folders than the station's")
    if run.get("variant") not in GCPFREE_VARIANTS:
        why.append(f"its run names no known variant ({run.get('variant')})")
    vm = VARIANT_LINE.search(notes)
    if not vm or vm.group(1) != run.get("variant"):
        why.append("its 'Variant' line is not its report's")
    # its numbers, and its lines about its references, against the files
    from georectify import load_extrinsics
    e = None
    try:
        e = load_extrinsics(p)
        b_ = [float(v) for v in (rep.get("B") or [])]
        ang = np.degrees(np.asarray(e[3:6], float))
        same = (len(b_) == 6 and np.allclose(np.asarray(e[:3], float), b_[:3], atol=2e-3)
                and all(abs((a - b + 180.0) % 360.0 - 180.0) <= 2e-3 for a, b in zip(ang, b_[3:6])))
    except Exception:
        same = False
    if not same:
        why.append("its numbers are not the B of its report")
    places = [rd] + [Path(str(x)) for x in dirs] + ([p.parents[3]] if len(p.parents) > 3 else []) \
        + [HERE / "calibration"]

    def ref_file(name, want_sha, what):
        """The reference `name` here, with the content the report recorded (and, for a GCP solution known here, its
        archived content), else None (said)."""
        q = next((Path(c) / name for c in places if name and (Path(c) / name).is_file()), None)
        if q is None:
            why.append(f"its {what} {name or '(none named)'} is not here (beside it or in a calibration folder): what it "
                       f"says of it cannot be checked")
            return None
        if file_sha(q) != str(want_sha or "").lower():
            why.append(f"the {name} here ({q}) is not the {what} its report says it was turned from (sha256)")
            return None
        known = CALIBRATION_ORIGINS.get(name)
        if known and known[0] != file_sha(q):
            why.append(f"the {name} here is not the archived {name} (sha256)")
            return None
        try:
            return load_extrinsics(q)
        except Exception:
            why.append(f"its {what} {name} cannot be read")
            return None
    ref_name = str(rep.get("b_reference_eo") or "")
    turns = list(ROTATED_REF.finditer(notes))
    if not turns:
        why.append("no '<reference>.yaml rotated by <angle> deg' line")
    e_ref = ref_file(ref_name, rep.get("b_reference_eo_sha256"), "reference")
    for m in turns:
        if Path(m.group(1)).name != ref_name:
            why.append(f"its 'rotated by' line names {Path(m.group(1)).name}, its report's reference is "
                       f"{ref_name or 'none'}")
        elif e is not None and e_ref is not None:
            turn = angle_between_deg(e, e_ref)
            if abs(turn - abs(float(m.group(2)))) > GCPFREE_TURN_TOL_DEG:
                why.append(f"its 'rotated by' line says {m.group(2)} deg from {ref_name}, but its numbers are turned "
                           f"{turn:.3f} deg from it")
    pm = POSITION_FROM.search(notes)
    if not pm or Path(pm.group(1)).name != ref_name:
        why.append(f"its 'Position (x, y, z) kept from' line does not name its report's reference ({ref_name})")
    elif e is not None and e_ref is not None:
        dpos = float(np.linalg.norm(np.asarray(e[:3], float) - np.asarray(e_ref[:3], float)))
        if dpos > 0.001:
            why.append(f"its camera position is {dpos:.3f} m from that of {ref_name}, whose position it says it keeps "
                       f"(changed after the run)")
    pn = PAN_FROM.search(notes)
    jan_name = str(rep.get("reference_eo") or "")
    if not pn or Path(pn.group(2)).name != jan_name:
        why.append(f"its 'pan change ... deg from' line does not name its report's Jan 2025 reference ({jan_name})")
    else:
        e_jan = e_ref if jan_name == ref_name else ref_file(jan_name, rep.get("reference_eo_sha256"), "Jan 2025 reference")
        if e is not None and e_jan is not None:
            dpan = float((np.degrees(e[3] - e_jan[3]) + 180.0) % 360.0 - 180.0)
            if abs(dpan - float(pn.group(1))) > GCPFREE_TURN_TOL_DEG or \
                    abs(dpan - float(rep.get("pan_change_deg") or 0.0)) > GCPFREE_TURN_TOL_DEG:
                why.append(f"its 'pan change' line says {pn.group(1)} deg from {jan_name} (its report "
                           f"{rep.get('pan_change_deg')}), but its azimuth differs from that of {jan_name} by "
                           f"{dpan:+.3f} deg")
    if not rj:
        why.append("no run.json beside it")
    else:
        cam_ = (rj.get("cameras") or {}).get(rep.get("camera")) or {}
        if rj.get("status") != "finished" or cam_.get("b_file_sha256") != sha:
            why.append("the run.json beside it does not record this B file (a finished run, its sha256)")
        if rj.get("run_id") != run.get("run_id") or rj.get("run_id") != rd.name:
            why.append(f"its run folder ({rd.name}), its run.json ({rj.get('run_id')}) and its report "
                       f"({run.get('run_id')}) name different runs")
    # where it was made: this repository's runs folder, under this repository's pre-registration, as it is now
    runs_here = Path(os.path.realpath(str(GCPFREE_DIR / "runs")))
    pr_here = Path(out["preregistration_here"])
    if rd.parent != runs_here:
        why.append(f"its run folder ({rd}) is not in this repository's runs folder ({runs_here}): a run under another "
                   f"pre-registration, in another folder or repository, is never the delivered B")
    if os.path.realpath(str(pr.get("path") or "")) != str(pr_here):
        why.append(f"its run was made under the pre-registration {pr.get('path')}, not this repository's ({pr_here})")
    if not pr_here.is_file():
        why.append(f"this repository's pre-registration ({pr_here}) does not exist")
    else:
        gp = git_file_state(pr_here)
        out["preregistration"]["git"] = {k: gp.get(k) for k in ("committed", "commit", "date", "commits", "why")}
        if sha256(pr_here) != str(pr.get("sha256") or "").lower():
            why.append("this repository's pre-registration is not the file its run recorded (sha256): changed after "
                       "the run")
        if not gp.get("committed"):
            why.append(f"this repository's pre-registration is {gp.get('why')}")
        elif gp.get("commits") != len(pr.get("git_history") or []):
            # by the number of commits, not their hashes: a history rewritten by a rebase keeps the count (and the
            # content is checked by its sha256 above)
            why.append(f"this repository's pre-registration was committed again after the run ({gp.get('commits')} "
                       f"commits now, {len(pr.get('git_history') or [])} when the run was made; the run recorded "
                       f"{str(pr.get('commit'))[:10]})")
        elif gp.get("toplevel") != repo_toplevel():
            why.append("this repository's pre-registration is committed in another git repository than this script's")
    # the first finished delivery run there (or the one allowed second run, with Aug photos)
    order = gcpfree_delivery_runs(runs_here)
    ids = [d.name for _, d in order]
    if rd.parent == runs_here:
        if rd.name not in ids:
            why.append(f"its run ({rd.name}) is not a finished delivery run of this repository")
        else:
            k = ids.index(rd.name)
            first = order[0][0]
            delivered_first = any((v or {}).get("delivery") == "DELIVERED"
                                  for v in (first.get("cameras") or {}).values())
            want = list(GCPFREE_VARIANTS)[0 if k == 0 or not delivered_first else 1]
            if k >= 2:
                why.append(f"it is delivery run {k + 1} of this repository: the pre-registered rule allows two at most, "
                           f"and B is the first's ({ids[0]})")
            elif k == 1 and not (not first.get("aug_dir") and run.get("aug_dir")
                                 and (first.get("preregistration") or {}).get("sha256") == pr.get("sha256")):
                why.append(f"it is the second delivery run of this repository ({rd.name}; the first is {ids[0]}): one "
                           f"is allowed only with Aug photos, after a first run without them, under the same "
                           f"pre-registration")
            elif run.get("run_no") != k + 1 or run.get("variant") != want:
                why.append(f"its report calls it delivery run {run.get('run_no')} (variant {run.get('variant')}); "
                           f"it is run {k + 1} of this repository, whose B is variant {want}")
    out["delivery_runs_here"] = ids
    if rd.parent == runs_here:
        # every run's record committed (a trial and a run cut off too: B's README lists them all)
        unc = gcpfree_uncommitted_runs(runs_here)
        if unc:
            why.append(f"the record of run(s) {', '.join(unc)} in this repository's runs folder is not committed (a "
                       f"trial, a run cut off: every run is committed, and B's README lists them): commit it")
        # nothing outside this branch's history: another branch, a stash or the reflog holding an earlier
        # pre-registration or run (e.g. a branch made from before the pre-registration)
        fr = gcpfree_foreign_records(GCPFREE_DIR)
        if fr:
            why.append(f"{len(fr)} commit(s) outside this branch's history wrote a pre-registration or a run record that "
                       f"this branch does not hold (" + ", ".join(f"{c[:10]} {q}" for c, q in fr[:3])
                       + (" ..." if len(fr) > 3 else "") + "): an earlier pre-registration or run may sit there, so "
                       "which run came first cannot be told")
        # the photos its run solved with are the ones the pre-registration pins (for the cameras it solved)
        rows_, sha_ = prereg_photos(pr_here)
        solved_ = set((rj or {}).get("solved") or [])
        want_ph = photos_sha256([r for r in rows_ if r and r[0] in solved_])
        if not rows_ or sha_ != photos_sha256(rows_) or not (rj or {}).get("photos_preregistered") \
                or (rj or {}).get("photos_used_sha256") != want_ph:
            why.append("its run did not use the photos its pre-registration pins (run.json photos_used_sha256 "
                       f"{str((rj or {}).get('photos_used_sha256'))[:12]}..., the pre-registration's for "
                       f"{', '.join(sorted(solved_)) or 'no camera'} {want_ph[:12]}...)")
    recs = gcpfree_run_records(runs_here) if rd.parent == runs_here else {}
    if recs.get("removed"):
        why.append(f"the record of run(s) {', '.join(recs['removed'])} was committed in this repository's runs folder "
                   f"and is gone now: every run is kept, so which delivery run came first cannot be told")
    other = [r for r in recs.get("changed") or [] if r != rd.name]
    if other:
        why.append(f"the record (run.json) of run(s) {', '.join(other)} was changed after its commit: every run is kept "
                   f"as written, so which delivery run came first cannot be told")
    # committed: the B file, its report and run.json, unchanged, each touched by one commit (the same)
    gs = {}
    for f, what in ((p, "B file"), (rp, "report"), (rd / "run.json", "run.json"), (rd / "run_report.txt",
                                                                                  "run_report.txt")):
        g = git_file_state(f)
        gs[what] = g
        name = f"its {what}" + ("" if f.name == what else f" {f.name}")
        if not g.get("committed"):
            why.append(f"{name} is {g.get('why')}: commit the run folder, in one commit, before the B DEM is built")
        elif g.get("commits") != 1:
            why.append(f"{name} was touched by {g.get('commits')} commits (one, the run's, is the rule: a file changed "
                       f"and committed again is not what the run wrote)")
        elif g.get("toplevel") != repo_toplevel():
            why.append(f"{name} is committed in another git repository than this script's")
    if all(g.get("committed") for g in gs.values()) and len({g.get("commit") for g in gs.values()}) > 1:
        why.append("its B file, report, run.json and run_report.txt were committed in different commits (the run "
                   "folder is committed once)")
    out["committed"] = {k: gs["B file"].get(k) for k in ("committed", "commit", "date", "commits", "why")}
    out["b_commit"] = gs["B file"].get("commit") if gs["B file"].get("committed") else None
    out["pinned"] = not why
    return out


def procedure_text(pc):
    """One line on a calibration B's procedure, for the README and the dry run."""
    if not pc or not pc.get("format"):
        return "not made by the pre-registered procedure (gcpfree_eo.py)"
    if not pc.get("pinned"):
        return "NOT the pre-registered procedure: " + "; ".join(pc.get("why") or [])
    pr = pc.get("preregistration") or {}
    cm = pc.get("committed") or {}
    return (f"made by the pre-registered procedure (gcpfree_eo.py sha256 {pc['procedure_sha256'][:12]}...), delivery "
            f"run {pc.get('run_no')} ({pc.get('run_id')}; " + (f"solved in run {pc['solved_in_run']}; "
                                                              if pc.get("solved_in_run") not in (None, pc.get("run_id"))
                                                              else "")
            + f"variant {pc.get('variant')}), verdict {pc.get('verdict')}, under the outcome rule "
            f"'{pc.get('outcome_rule')}' (pre-registration {pc.get('preregistration_here')} sha256 "
            f"{str(pr.get('sha256'))[:12]}..., {pr.get('commits')} commit(s) before the run, the last "
            f"{str(pr.get('commit'))[:10]}, decided by {pr.get('decided_by')}); photo folders the station's three"
            + (f", Aug photos {pc['aug_dir']}" if pc.get("aug_dir") else "")
            + f"; run folder {pc.get('run_dir')}; B file sha256 {str(pc.get('b_sha256'))[:12]}..., committed to git "
              f"once ({str(cm.get('commit'))[:10]}, {cm.get('date')}) with its report and run.json")


def uncertainty_terms(pc):
    """The two stated uncertainty terms of a calibration B (its report), as text, or ''."""
    u = (pc or {}).get("uncertainty") or {}
    if not u:
        return ""
    pos = u.get("position") or {}
    return (f"rotation +/- {u.get('rotation_px')} px on the beach part of the frame (~{u.get('rotation_deg')} deg; "
            f"{u.get('rotation_rule')})"
            + (f", {u['beach_m']} m on the beach" if u.get("beach_m") else "")
            + f"; inherited camera position: in the largest of five synthetic position-inheritance cases, DEM |dz| "
              f"median {pos.get('dem_abs_median_m')} m, p90 {pos.get('dem_p90_m')} m, horizontal "
              f"{pos.get('horizontal_median_m')} m (not a bound; the GCP solutions of Aug 2024 and Jan 2025 put this "
              f"camera {pos.get('reference_solution_spread')} apart); NEITHER term is a validated bound on real "
              f"photos: both come from synthetic tests")


def tilt_effect_text(rt):
    """The README line on a GCP-free B's rotation term taken as a pure tilt (write_provenance: rotation_term_as_tilt)."""
    if not rt or not rt.get("bands"):
        return (f"its DEM effect as a tilt: not computed ({(rt or {}).get('why_none') or 'no DEM yet'}); a rotation of "
                f"that size as tilt moves waterlines metres at 200-400 m from the camera")
    return (f"the rotation term as a pure TILT of {rt['tilt_deg']} deg, from this date's DEM and calibration (median per "
            f"range band from the camera): "
            + "; ".join(f"{b['range_m']} m: waterlines move {b['shift_m']} m, DEM |dz| "
                        f"{float(b.get('dem_abs_error_m', abs(float(b['dem_error_m'])))):.2f} m" for b in rt["bands"])
            + (" (local slopes from one-sided differences: a thin strip of DEM)" if any(b.get("slope") for b in
                                                                                     rt["bands"]) else "")
            + " (the DEM error is the local slope times the move; the rotation term is not a validated bound, nor "
              "is this)")


def utc_time(text):
    """A time as this project writes them (git's ISO with an offset, a stamp's '...Z', a run's 'YYYY-MM-DD HH:MM:SS'
    UTC) as an aware UTC datetime, else None."""
    try:
        t = str(text or "").strip().replace(" ", "T")
        if t.endswith("Z"):
            t = t[:-1] + "+00:00"
        d = datetime.fromisoformat(t)
    except ValueError:
        return None
    return (d if d.tzinfo else d.replace(tzinfo=timezone.utc)).astimezone(timezone.utc)


def utc_text(d):
    return d.strftime("%Y-%m-%d %H:%M UTC") if d else "not known"


def pin_commit(procedure_sha):
    """The first commit of this repository whose survey_products.py holds this gcpfree_eo.py sha256 (its pin), and
    its date. -> (commit, date) or (None, None)."""
    try:
        r = subprocess.run(["git", "-C", str(HERE), "log", "--reverse", "--format=%H %cI", "-S", str(procedure_sha),
                            "--", "survey_products.py"], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           universal_newlines=True)
    except OSError:
        return None, None
    first = next((ln.split(" ", 1) for ln in r.stdout.splitlines() if ln.strip()), None) if r.returncode == 0 else None
    return (first[0], first[1] if len(first) > 1 else None) if first else (None, None)


def argv_options(words):
    """{option: [values]} of a command line split into words (state_argv's form; the words before the first option,
    e.g. the interpreter and the script, are left out)."""
    out, k = {}, None
    for w in words or []:
        if str(w).startswith("--"):
            k = str(w)
            out[k] = []
        elif k is not None:
            out[k].append(str(w))
    return out


def a_product(args, date, variant):
    """The table-calibration (A) product of a blind date, from its root (--a-root, else the root this variant's root
    is named after): when its DEM FIRST existed and when it was last rebuilt, its build options (provenance: the
    command that rebuilds it), and per camera the search envelope and forcing inputs its steps used (stamps). Only
    these records are read: never A's calibration, grids or waterlines. -> dict."""
    root = Path(getattr(args, "a_root", None) or root_for_variant(args.output_root, None, variant))
    out = root / date
    st, prov = read_stamp(out, "dem") or {}, load_json(out / "provenance.json") or {}
    var = sorted({(v.get("calibration") or {}).get("variant") or "table" for v in (prov.get("cameras") or {}).values()})
    times = [x for x in (first_built(out, "dem"), st.get("first_finished_utc"), st.get("finished_utc")) if x]
    cmd = (prov.get("commands") or {}).get("survey_products")
    try:
        words = shlex.split(cmd) if cmd else None
    except ValueError:
        words = None
    env, cams = {}, sorted((prov.get("cameras") or {}).keys())
    for cam in cams:
        ds = read_stamp(out, f"detect_{cam}") or {}
        sig = (ds.get("sig") or {}).get("envelope")
        if ds:
            env[cam] = {"path": sig[0] if sig else None, "sig": sig, "sha256": ds.get("envelope_sha256")}
    return {"root": str(root), "found": bool(st or prov), "variant": "/".join(var) or None,
            "dem_built_utc": min(times) if times else None, "dem_last_built_utc": st.get("finished_utc"),
            "provenance_built_utc": prov.get("built_utc"), "options": argv_options(words) if words else None,
            "command": cmd, "envelope": env, "cameras": cams,
            "forcing_inputs": ((read_stamp(out, "forcing") or {}).get("sig") or {}).get("inputs"),
            "forcing_sha256": {k: file_sha(out / "forcing" / k) for k in ("water_level.csv", "waves.csv")},
            "software": prov.get("software")}


def a_product_times(args, date, variant):
    """a_product() (its times are what B's ORDER states)."""
    return a_product(args, date, variant)


# What makes a GCP-free B independent of the held survey, said plainly to its holder and to the owner (B's README
# BLIND section and the share MANIFEST): what the code enforces, and what rests on the owner following the
# procedure -- software on a computer the owner controls cannot prove that a rule was fixed before anyone looked.
B_PROTOCOL_ENFORCED = (
    "WHAT THE CODE ENFORCES (survey_products.py and gcpfree_eo.py, on this computer): B is made only by the pinned "
    "procedure (gcpfree_eo.py, its five modules and its constants, by sha256: GCPFREE_PINS); every run of it, a trial "
    "too, is refused unless this repository's calibration/gcpfree/PREREGISTERED.txt (the outcome rule, and the photos "
    "a run uses, by sha256) is committed, unchanged and pushed (its commit on a remote-tracking branch: this "
    "computer's record that a push was made, not of when GitHub received it), and every run is recorded in "
    "calibration/gcpfree/runs/; the "
    "first delivery run fixes B (one more only with Aug 2024 photos, for a camera not delivered); B's file, report "
    "and run record must be committed once and unchanged since; B's DEM is built as A was (the committed tables, A's "
    "build options, the search envelope A's detection used, A's forcing), in one root; the procedure never opens "
    "calibration A, a GCP or target file or the calibration history, and no build looks for the held survey; the "
    "ORDER above is read from git and from the run records")
B_PROTOCOL_OWNER = (
    "WHAT RESTS ON THE OWNER FOLLOWING THE PROCEDURE (no software on a computer the owner controls can prove it): that "
    "PREREGISTERED.txt was committed AND pushed to GitHub before any run of gcpfree_eo.py, a trial included (the "
    "procedure refuses a run until a remote-tracking branch holds the commit, but only GitHub's receipt time shows "
    "when the push reached it); that no "
    "copy of the procedure was run anywhere else (another clone, a copied script) before it; that no run folder, "
    "branch or reflog entry was deleted and the repository not re-cloned; that no B file was edited before its commit; "
    "that this computer's clock, which dates every git commit and run above, was right. What answers these comes from "
    "outside this computer: the sha256s of the procedure and of the pre-registration, sent to the holder before A's "
    "DEM was shared; each B file's sha256, sent the moment its run ended; and the time GitHub received the "
    "pre-registration commit (below)")


def commit_remote_state(path, commit):
    """Where the commit `commit` of the repository holding `path` is known to be pushed: the remote-tracking branches
    holding it as this computer last fetched them, and origin's URL (any credentials left out). -> {"branches": [...],
    "url": str or None, "why": ''}."""
    q = Path(os.path.realpath(str(path)))

    def git(*a):
        return subprocess.run(["git", "-C", str(q.parent)] + list(a), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True)
    if not commit:
        return {"branches": [], "url": None, "why": "no commit"}
    try:
        br = git("branch", "-r", "--contains", str(commit))
        url = git("remote", "get-url", "origin")
    except OSError as exc:
        return {"branches": [], "url": None, "why": f"git cannot be run ({exc})"}
    names = [x.strip() for x in br.stdout.splitlines() if x.strip() and "->" not in x] if br.returncode == 0 else []
    u = url.stdout.strip() if url.returncode == 0 and url.stdout.strip() else None
    return {"branches": names, "url": re.sub(r"//[^/@]*@", "//", u) if u else None, "why": ""}


def prereg_commit_text(pc):
    """The pre-registration's commit as B's holder is told it: its hash, its date by this computer's clock, whether a
    remote branch this computer has fetched holds it, and where its public time is."""
    pr = (pc or {}).get("preregistration") or {}
    c = pr.get("commit")
    rs = commit_remote_state(pc.get("preregistration_here") or GCPFREE_DIR / "PREREGISTERED.txt", c)
    return (f"PRE-REGISTRATION COMMIT: {c or '(none found)'} ({pr.get('commit_date')} by this computer's clock; the "
            f"file's sha256 {str(pr.get('sha256'))[:12]}...); "
            + (f"held by {', '.join(rs['branches'])} of {rs['url'] or 'the remote'} as this computer last fetched it"
               if rs["branches"] else
               "on NO remote branch this computer has fetched (not pushed, or not fetched since)")
            + ". A commit's own date is set by the computer that made it: the public timestamp is the time GitHub "
              "received the push, which the holder can check on the repository's GitHub page (its Activity, or the "
              "GitHub events API) against this commit hash and the sha256 he was sent")


def gcpfree_disclosure(plan, args, pcs):
    """What the holder of a blind date's survey is told of a GCP-free B, besides its procedure line: the EXPOSURE
    statement of its procedure (from its report), and the ORDER of what decides B and of A: the procedure's pin, the
    pre-registration of the outcome rule, A's DEM, the B run and B's commit, with their dates and what came first.
    pcs: {camera: procedure_check() of a pinned B}. -> {"exposure", "order": [lines], "times": {...}}."""
    pc = next(iter(pcs.values()))
    pr = pc.get("preregistration") or {}
    pin_c, pin_d = pin_commit(pc.get("procedure_sha256"))
    a = a_product(args, plan["date"], plan_variant(plan))
    # the rule's own history: when the rule in force was first committed (an earlier rule, or an edit of who decided,
    # is said); the run recorded the pre-registration's last commit
    hist = prereg_rule_history(pc.get("preregistration_here") or GCPFREE_DIR / "PREREGISTERED.txt")
    since = None
    for h in reversed(hist):                       # newest first: back to the first commit of the rule in force
        if h.get("rule") != pc.get("outcome_rule"):
            break
        since = h
    t = {"pin": utc_time(pin_d), "rule": utc_time((since or {}).get("date") or pr.get("commit_date")),
         "prereg": utc_time(pr.get("commit_date")), "a_dem": utc_time(a.get("dem_built_utc")),
         "a_dem_last": utc_time(a.get("dem_last_built_utc")),
         "b_run": utc_time(pc.get("run_started_utc")), "b_end": utc_time(pc.get("run_finished_utc")),
         "b_commit": utc_time((pc.get("committed") or {}).get("date"))}
    # every run of the procedure in this repository's runs folder (a trial, a run cut off), with its start: a trial
    # is a full result, so one made before B's run, or any run made before the rule was committed, is said
    runs = []
    rdir = Path(pc.get("run_dir") or GCPFREE_DIR / "runs" / "x").parent
    for d_ in sorted(q for q in rdir.glob("*") if q.is_dir()) if rdir.is_dir() else []:
        js = load_json(d_ / "run.json") or {}
        runs.append({"run_id": d_.name, "kind": "delivery" if js.get("eligible") else "TRIAL" if js else "unreadable",
                     "status": js.get("status") or "no run.json", "started": utc_time(js.get("started_utc")),
                     "finished": utc_time(js.get("finished_utc")), "why": js.get("trial_why"),
                     "this": d_.name == Path(str(pc.get("run_dir") or "")).name})
    order = [f"procedure pinned (gcpfree_eo.py sha256 {str(pc.get('procedure_sha256'))[:12]}... in survey_products.py "
             f"GCPFREE_PINS): " + (f"commit {pin_c[:10]}, {utc_text(t['pin'])}" if pin_c else
                                   "its commit not found in this repository's history"),
             f"outcome rule '{pc.get('outcome_rule')}' pre-registered (PREREGISTERED.txt, decided by "
             f"{pr.get('decided_by')} on {pr.get('decided_on')}): first committed with this rule "
             f"{str((since or {}).get('commit') or pr.get('commit'))[:10]}, {utc_text(t['rule'])}"
             + ("; the file's commits: " + ", ".join(f"{str(h['commit'])[:10]} {utc_text(utc_time(h['date']))} "
                                                       f"'{h['rule']}'" for h in hist) if len(hist) > 1 else ""),
             (f"A (the table's calibration) DEM first built {utc_text(t['a_dem'])}"
              + (f" (last rebuilt {utc_text(t['a_dem_last'])})" if t["a_dem_last"] and t["a_dem_last"] != t["a_dem"]
                 else "") + f" ({a['root']}/{plan['date']}"
              + (f"; variant {a['variant']}" if a.get("variant") not in (None, "table") else "") + ")"
              if a["found"] and t["a_dem"] else
              f"A (the table's calibration): its product at {a['root']}/{plan['date']} has no finished DEM yet (being "
              f"built, or failed): its build time is not stated here (build this README again once it has one)"
              if a["found"] else
              f"A (the table's calibration): no product found at {a['root']}/{plan['date']} (its build time is not "
              f"stated here; give --a-root)"),
             f"B run {pc.get('run_id')}: {utc_text(t['b_run'])} to {utc_text(t['b_end'])}",
             f"B committed: {', '.join(sorted({str((p.get('committed') or {}).get('commit'))[:10] for p in pcs.values()}))}, "
             f"{utc_text(t['b_commit'])}",
             f"every run of the procedure in {rdir} ({len(runs)}): " + "; ".join(
                 f"{r['run_id']} {r['kind']}" + (" (B's run)" if r["this"] else "")
                 + f", started {utc_text(r['started'])}, "
                 + (f"finished {utc_text(r['finished'])}" if r["status"] == "finished" else f"status {r['status']}")
                 for r in runs)]
    flags = []
    early = [r for r in runs if r["started"] and t["rule"] and r["started"] < t["rule"]]
    if early:
        flags.append(f"run(s) {', '.join(r['run_id'] for r in early)} started BEFORE the outcome rule in force was "
                     f"committed: a result could have been seen when the rule was chosen")
    trials = [r for r in runs if r["kind"] != "delivery" and r["started"] and t["b_run"] and r["started"] < t["b_run"]]
    if trials:
        flags.append(f"TRIAL run(s) {', '.join(r['run_id'] for r in trials)} (a full result: every verdict and angle) "
                     f"were made before B's run {pc.get('run_id')}, under the same committed rule")
    cut = [r for r in runs if r["status"] != "finished"]
    if cut:
        flags.append(f"run(s) {', '.join(r['run_id'] for r in cut)} never finished (cut off, or still running): what "
                     f"they had computed was printed only")
    if not a["found"]:
        flags.append(f"the order of B relative to A is not known: no table-calibration (A) product was found at "
                     f"{a['root']}/{plan['date']} (give --a-root)")
    if len({h.get("rule") for h in hist}) > 1:
        flags.append("the outcome rule of PREREGISTERED.txt was changed before the run ("
                     + " -> ".join(f"'{h['rule']}' {utc_text(utc_time(h['date']))}" for h in hist) + ")")
    if t["a_dem"] and t["rule"] and t["a_dem"] < t["rule"]:
        flags.append("A's DEM was built BEFORE the outcome rule in force was pre-registered (it could have been seen "
                     "when the rule was chosen; the procedure's code was fixed by its pin)")
    if t["a_dem"] and t["pin"] and t["a_dem"] < t["pin"]:
        flags.append("A's DEM was built BEFORE the procedure was pinned")
    if t["b_commit"] and t["a_dem"] and t["b_commit"] > t["a_dem"]:
        flags.append("B was committed after A's DEM was built: only B's sha256, written down (and sent) when its run "
                     "ended, shows it was not changed in between")
    if t["a_dem_last"] and t["b_run"] and t["a_dem_last"] > t["b_run"]:
        # B's DEM is built with A's options as A was LAST built: a build of A after B's run could have set them with
        # B's calibration known (a rebuild after a code change keeps them: A's provenance gives its commands). A DEM
        # built for the FIRST time after B's run is not a rebuild (review r1: 'rebuilt' next to 'first built' then)
        first_after = not t["a_dem"] or t["a_dem"] >= t["a_dem_last"] or t["a_dem"] > t["b_run"]
        flags.append((f"A's DEM was first built ({utc_text(t['a_dem'] or t['a_dem_last'])}) after B's run began: its "
                      f"options could have been chosen with B known, and B's DEM is built with them; A's provenance "
                      f"gives the command of that build"
                      + (f" (and of its last rebuild, {utc_text(t['a_dem_last'])})"
                         if t["a_dem"] and t["a_dem_last"] != t["a_dem"] else "")) if first_after else
                     f"A's DEM was rebuilt ({utc_text(t['a_dem_last'])}) after B's run began: B's DEM is built with "
                     f"A's options as last built, which could then have been chosen with B known; A's provenance gives "
                     f"the command of that rebuild (a rebuild after a code change keeps A's options)")
    def ran(sw):                                   # {script: {sha256 it ran with}} of a product's steps
        out_ = {}
        for v in ((sw or {}).get("by_step") or {}).values() if isinstance(sw, dict) else []:
            for s, h in ((v or {}).get("scripts") or {}).items():
                out_.setdefault(s, set()).add(h)
        return out_
    sw_a, sw_b = ran(a.get("software")), ran(software_of_outputs(plan["out"]) if plan.get("out") else None)
    other_code = sorted(k for k in set(sw_a) & set(sw_b) if sw_a[k] != sw_b[k])
    if other_code:
        flags.append(f"B was built with other code than A ({', '.join(other_code)}: see the software of both "
                     f"provenance files)")
    return {"exposure": pc.get("exposure") or "the procedure's report states no exposure",
            "order": order, "order_flags": flags,
            "protocol": [B_PROTOCOL_ENFORCED, B_PROTOCOL_OWNER, prereg_commit_text(pc)],
            "times": {k: (v.strftime("%Y-%m-%dT%H:%M:%SZ") if v else None) for k, v in t.items()},
            "a_product": {k: v for k, v in a.items() if k not in ("options", "software")},
            "runs": [dict(r, started=utc_text(r["started"]), finished=utc_text(r["finished"])) for r in runs],
            "rule_history": hist}


def not_delivered_cameras(date_cams, given, dirs=()):
    """--eo for some of a date's cameras only (a one-camera B): each camera left out must be one the
    PRE-REGISTERED procedure did NOT DELIVER in the very run the given B files come from (its report beside them
    says so); a delivered camera is never left out, and only the procedure's own B files may leave one out.
    -> ({camera: what the procedure said}, {camera: why leaving it out is not accepted})."""
    missing = [c for c in date_cams if c not in given]
    ok, bad = {}, {}
    if not missing:
        return ok, bad
    procs = {c: procedure_check(pth, dirs) for c, pth in given.items() if c in date_cams and Path(pth).is_file()}
    runs = {pc.get("run_dir") for pc in procs.values()}
    for cam in missing:
        if not procs or not all(pc.get("pinned") for pc in procs.values()):
            notb = "; ".join(f"{c}: {'; '.join(pc.get('why') or [])}" for c, pc in sorted(procs.items())
                             if not pc.get("pinned"))
            bad[cam] = ("only calibration B files of the pre-registered procedure (gcpfree_eo.py, as pinned here), as "
                        "it delivered them here, may leave out a camera the procedure did not deliver"
                        + (f"; these are not (or not yet): {notb}" if notb else ""))
            continue
        if len(runs) != 1:
            bad[cam] = "the --eo files come from different runs of the procedure"
            continue
        rd = Path(next(iter(runs)))
        found = None
        for f in sorted(rd.glob("*_report.json")):
            try:
                rep = json.loads(f.read_text())
            except (OSError, ValueError):
                continue
            if rep.get("camera") == cam:
                found = (f, rep)
                break
        if found is None:
            bad[cam] = f"no report of {cam} in the run folder of the given B files ({rd})"
            continue
        f, rep = found
        if rep.get("delivery") == "DELIVERED":
            bad[cam] = (f"the procedure DELIVERED {cam} in that run ({f.name}): a delivered camera is never left out "
                        f"of its B")
            continue
        ok[cam] = {"verdict": rep.get("verdict"), "delivery": rep.get("delivery"), "reasons": rep.get("reasons") or [],
                   "report": str(f), "report_sha256": file_sha(f), "run_id": rd.name}
    return ok, bad


def calibration_origin(eo_path, dirs=(), _depth=0, name_only=False):
    """
    What a calibration (EO file) was solved from: CACO03's GCP solutions are known by name and content
    (CALIBRATION_ORIGINS); a calibration turned from others by matching photos is traced through the
    reference(s) its comment lines name (ROTATED_REF: looked for next to it, then in `dirs`), down to
    those. Only comment lines are read for the trace, never a calibration's numbers.
    -> {"file", "kind": 'gcp' | 'rotated' | 'unknown', "gcp_dates": [dates of the GCP surveys its
    geometry rests on], "text", "chain": [file names], "position_from": the solution whose camera
    position it keeps, "acceptance": calibration_acceptance() of this file, "notes"}.
    'rotated' only when every reference traces to a known GCP solution; else 'unknown'.
    """
    p = Path(str(eo_path)) if eo_path else None
    name = p.name if p else "?"
    out = {"file": name, "kind": "unknown", "gcp_dates": [], "chain": [name], "position_from": None,
           "acceptance": None, "notes": ""}
    reg = CALIBRATION_ORIGINS.get(name)
    if name_only:                          # a reference named in notes and not found: known by its name only
        p = None
    sha = file_sha(p) if p is not None and p.is_file() else None
    if reg and (sha is None or sha == reg[0]):
        out.update(kind="gcp", gcp_dates=[reg[1]], position_from=name,
                   text=f"solved from {reg[2]}" + ("" if sha else " (known by its name: the file is not here)"))
        return out
    if reg:
        out["text"] = (f"NOT KNOWN: its content differs from the archived {name} (sha256), so what it was solved "
                       f"from is not known")
        return out
    if p is None or not p.is_file():
        out["text"] = f"NOT KNOWN: {name} is not here and its name is not in CALIBRATION_ORIGINS"
        return out
    lines = eo_comment_lines(p)
    out["notes"] = " ".join(lines)
    out["acceptance"] = calibration_acceptance(lines)
    out["procedure"] = procedure_check(p, dirs)  # made by the pre-registered GCP-free procedure? (held_check)
    claims = rotation_claims(out["notes"], name)
    refs, pos_ref = claims["refs"], claims["position_ref"]
    if not refs or _depth >= 4:
        out["text"] = ("NOT KNOWN: not in survey_products.py CALIBRATION_ORIGINS, and its notes name no reference "
                       "it was turned from ('<file>.yaml rotated by <angle> deg', or gcpfree_eo.py's 'Reference EO: "
                       "<file>.yaml ... sha256 ...')" if not refs else
                       "NOT KNOWN: a chain of turned calibrations longer than 4")
        return out
    subs, parts, wrong = [], [], []
    try:
        from georectify import load_extrinsics, cirn_angles_to_R
        e_self = load_extrinsics(p)
    except Exception:                      # not a readable EO: nothing it says can be checked
        e_self = None
        wrong.append("its numbers cannot be read")
    for r in refs:
        ref = r["name"]
        q = next((Path(d) / ref for d in [p.parent] + [str(x) for x in dirs] if (Path(d) / ref).is_file()), None)
        sub = calibration_origin(q or ref, dirs, _depth + 1, name_only=q is None)
        subs.append(sub)
        how = (f"turned {r['angle']} deg" if r["angle"] is not None else
               f"turned {float(r['pan']):+.3f} deg in pan" if r["pan"] is not None else "turned")
        parts.append(f"{ref} ({sub['text']}) {how}")
        # what the notes say must be true of the files: the reference's content, the turn (or pan change) from
        # each reference, and the camera position of the one it keeps; a position from elsewhere (e.g. the held
        # survey's solution) is caught here
        if e_self is None:
            continue
        if q is None:
            wrong.append(f"{ref} is not here, so the turn and position it claims cannot be checked")
            continue
        if r["sha256"] and file_sha(q) != r["sha256"]:
            wrong.append(f"the {ref} here is not the one it was turned from (its sha256 is not the one its notes "
                         f"state)")
        try:
            e_ref = load_extrinsics(q)
            Rs, Rr = cirn_angles_to_R(*e_self[3:6]), cirn_angles_to_R(*e_ref[3:6])
            turn = float(np.degrees(np.arccos(np.clip((np.trace(Rs @ Rr.T) - 1.0) / 2.0, -1.0, 1.0))))
            dpan = float((np.degrees(e_self[3] - e_ref[3]) + 180.0) % 360.0 - 180.0)
        except Exception:
            wrong.append(f"{ref} cannot be read")
            continue
        if r["angle"] is not None and abs(turn - abs(float(r["angle"]))) > 0.05:
            wrong.append(f"its notes say {r['angle']} deg from {ref}, but it is turned {turn:.2f} deg from it")
        if r["pan"] is not None and abs(dpan - float(r["pan"])) > 0.05:
            wrong.append(f"its notes say a pan change of {float(r['pan']):+.3f} deg from {ref}, but its azimuth "
                         f"differs from that of {ref} by {dpan:+.3f} deg")
        dpos = float(np.linalg.norm(np.asarray(e_self[:3]) - np.asarray(e_ref[:3])))
        if ref == pos_ref and dpos > 0.001:
            wrong.append(f"its camera position is {dpos:.3f} m from that of {ref}, whose position it claims to "
                         f"keep: it comes from elsewhere")
    known = all(x["kind"] in ("gcp", "rotated") for x in subs) and not wrong
    pos_sub = next((x for x, r in zip(subs, refs) if r["name"] == pos_ref), subs[0])
    pos_from = pos_sub.get("position_from") or pos_ref
    out.update(kind="rotated" if known else "unknown", refs=[r["name"] for r in refs], checked=wrong or None,
               gcp_dates=sorted({dd for x in subs for dd in x["gcp_dates"]}),
               chain=[name] + [c for x in subs for c in x["chain"]],
               position_from=pos_from, rotation_source=claims["source"],
               text=("made WITHOUT ground control of its own: " + "; and ".join(parts)
                     + " to this pointing by matching photos"
                     + (f" (rotation source: {claims['source']})" if claims["source"] else "")
                     + ", its camera position (x, y, z) that of " + str(pos_from)
                     + ("" if known else
                        "; NOT KNOWN: " + "; ".join(wrong) if wrong else "; a reference's origin is NOT KNOWN")))
    return out


def held_check(origin, b, envelope=None):
    """What the holder's check of a BLIND product can test, given one camera's calibration origin and (envelope:
    envelope_check()) the search envelope that placed its detection. -> (label, text). The label is that of a
    comparison of this DEM with the held survey: the build never read it, so what is left is the camera geometry
    (the calibration) and the envelope (a DEM built with calibration A, or a file not known, carries A's geometry
    into the search band: PARTLY-CIRCULAR)."""
    lab, text = held_check_calibration(origin, b)
    if envelope and not envelope.get("ok") and lab == "INDEPENDENT":
        return "PARTLY-CIRCULAR", (f"{text}; BUT {envelope.get('text')}: the envelope may rest on the held survey's "
                                   f"calibration, so his check of this DEM is not called independent")
    if envelope and envelope.get("ok") and lab == "INDEPENDENT" and envelope.get("text"):
        text += f"; {envelope['text']}"
    return lab, text


def held_check_calibration(origin, b):
    """held_check by the calibration alone."""
    who = b["held_by"].split(" (")[0]
    head = (f"{origin['file']} {origin['text']}" if origin["kind"] == "gcp" else f"{origin['file']}: {origin['text']}")
    pc_ = origin.get("procedure") or {}
    if origin["kind"] == "unknown":
        return "PARTLY-CIRCULAR", (f"{head}: what its geometry rests on is not known, so it may rest on "
                                   f"{b['survey']}: a check of this DEM cannot be called independent of the camera "
                                   f"geometry"
                                   + (f"; and it is NOT the pre-registered procedure's B ({'; '.join(pc_['why'])})"
                                      if pc_.get("format") and not pc_.get("pinned") and pc_.get("why") else ""))
    if b["survey_date"] in origin["gcp_dates"]:
        return "PARTLY-CIRCULAR", (f"{head}: that is the survey {who} holds, so the camera geometry is fitted to his "
                                   f"points and his check of this DEM tests the water level, the detection and the "
                                   f"wave setup, not the geometry (PARTLY-CIRCULAR for the geometry)")
    pc = origin.get("procedure") or {}
    if origin["kind"] == "rotated" and not pc.get("pinned"):
        return "PARTLY-CIRCULAR", (f"{head}: NOT the pre-registered procedure ({'; '.join(pc.get('why') or [])}): a "
                                   f"calibration turned by hand, by another procedure or by an edited one could have "
                                   f"been chosen after seeing the held survey's calibration, so his check of this DEM "
                                   f"is not called independent of the camera geometry")
    acc = origin.get("acceptance") or {}
    return "INDEPENDENT", (f"{head}: nothing in it comes from {b['survey']}, so his check tests the camera geometry "
                           f"too (INDEPENDENT"
                           + (f"; {procedure_text(pc)}" if pc.get("pinned") else "") + ")"
                           + (f"; its pre-registered acceptance: {acc['verdict']}"
                              + (f" ({acc['text']})" if acc.get("text") else "") if acc else
                              "; it states no pre-registered acceptance (no '# acceptance: ...' line, no "
                              "gcpfree_eo.py 'CALIBRATION B ...: DELIVERED' line)"
                              if origin["kind"] == "rotated" else ""))


def built_calibration(plan):
    """{camera: {"eo_file", "eo_path", "eo_sha256", "variant"}} of the pointing and waterlines ON DISK
    (each camera's detection stamp, else its pointing stamp; one written before Oct 2026 gives the path
    in its signature only). {} where nothing is built."""
    res = {}
    for cam in plan["cams"]:
        for key in (f"detect_{cam}", f"pointing_{cam}"):
            st = read_stamp(plan["out"], key) if plan.get("out") else None
            if not st:
                continue
            cal = dict(st.get("calibration") or {})
            sig_eo = (st.get("sig") or {}).get("eo_file") or [None]
            cal.setdefault("eo_path", sig_eo[0])
            if not cal.get("eo_file") and cal.get("eo_path"):
                cal["eo_file"] = Path(cal["eo_path"]).name
            cal["eo_sha256"] = cal.get("eo_sha256") or st.get("eo_sha256")
            cal["step"] = key
            res[cam] = cal
            break
    return res


def calibration_differences(plan, built=None):
    """This run's calibration against the one the outputs on disk were built with, per camera (by
    content; by path for a stamp that has no content). -> list of text (empty when they agree)."""
    built = built_calibration(plan) if built is None else built
    diff = []
    for cam, c in plan["cams"].items():
        b = built.get(cam)
        if not b:
            continue
        now = file_sha(c["eo_path"]) if c.get("eo_path") else None
        same = (b["eo_sha256"] == now) if b.get("eo_sha256") and now else \
            (str(b.get("eo_path")) == str(c.get("eo_path")))
        if not same:
            diff.append(f"{cam} {c['eo_file']}" + (f" (variant {c['variant']})" if c.get("variant") else "")
                        + f", built with {b.get('eo_file')}" + (f" (variant {b['variant']})" if b.get("variant") else ""))
    return diff


def root_for_variant(root, variant, held_variant=None):
    """The output root a build with this calibration variant belongs in: the table's calibration in the
    base root, a variant in <base>_<variant> (before any _c<C> suffix, so root_for_c still finds its base).
    held_variant: the variant this root holds (its suffix is taken off first). Never this very root."""
    root = Path(root)
    m = re.match(r"^(.*?)(_c\d+(?:\.\d+)?)?$", root.name)
    stem, csuf = m.group(1), m.group(2) or ""
    if held_variant and stem.endswith("_" + held_variant):
        stem = stem[:-len(held_variant) - 1]
    alt = root.parent / (stem + (f"_{variant}" if variant else "") + csuf)
    if alt == root:
        alt = root.parent / f"{root.name}_{variant or 'table'}"
    return alt


def root_variants(root):
    """{variant ('table' for the table's calibration): [dates]} of the products under `root`, from their
    pointing and detection stamps (a stamp written before Oct 2026 records none: 'table')."""
    found = {}
    for d in sorted(Path(root).iterdir()) if Path(root).is_dir() else []:
        sdir = d / "logs" / "stamps"
        if not sdir.is_dir():
            continue
        for st in sorted(sdir.glob("pointing_c*.json")) + sorted(sdir.glob("detect_c*.json")):
            v = ((load_json(st) or {}).get("calibration") or {}).get("variant") or "table"
            if d.name not in found.setdefault(v, []):
                found[v].append(d.name)
    return found


def variant_root_problem(plan, args):
    """One calibration variant per output root, told by its name: a --calibration-variant NAME build goes in
    a root named <...>_NAME (optionally + _c<C>), and no root mixes variants across its dates (a table
    build is refused in a root holding a variant, and the reverse). -> text ('' when fine). Refused, exit
    5, whatever --replace says: the remedy is the right --output-root, not overwriting."""
    variant = next((c.get("variant") for c in plan["cams"].values() if c.get("variant")), None)
    root = Path(args.output_root)
    if variant and not re.search(rf"_{re.escape(variant)}(_c\d+(\.\d+)?)?$", root.name):
        return (f"a --calibration-variant {variant} build goes in its own root, named for it: --output-root "
                f"{root_for_variant(root, variant)} (not {root}: one calibration variant per root, so a "
                f"product's root says which calibration built it)")
    other = {v: ds for v, ds in root_variants(root).items() if v != (variant or "table")}
    if other:
        held = ", ".join(f"{v} ({', '.join(ds)})" for v, ds in sorted(other.items()))
        return (f"{root} holds products built with another calibration variant: {held}. One calibration variant "
                f"per root: build this one "
                + (f"(--calibration-variant {variant}) in --output-root {root_for_variant(root, variant)}" if variant
                   else f"(the table's calibration) in --output-root "
                        f"{root_for_variant(root, None, next(iter(sorted(other))))}"))
    return ""


def calibration_replace_note(plan, args, steps):
    """A run that would rebuild the pointing or the waterlines of a built date with another calibration
    in the SAME folder: what it would replace and the root this calibration belongs in. -> text or None.
    Refused unless --replace (build_date), as another C is."""
    if not ({"pointing", "detect"} & set(steps)):
        return None
    built = built_calibration(plan)
    diff = calibration_differences(plan, built)
    if not diff:
        return None
    variant = next((c.get("variant") for c in plan["cams"].values() if c.get("variant")), None)
    held = next((b.get("variant") for b in built.values() if b.get("variant")), None)
    same_file = all(str((built.get(cam) or {}).get("eo_path")) == str(c.get("eo_path"))
                    for cam, c in plan["cams"].items() if built.get(cam))
    if variant and variant == held and same_file and \
            re.search(rf"_{re.escape(variant)}(_c\d+(\.\d+)?)?$", Path(args.output_root).name):
        # the same variant and file in its own root: the calibration FILE changed since the build (another root
        # would not help)
        return (f"the product in {plan['out']} was built with another calibration ({'; '.join(diff)}): the same --eo "
                f"file, changed since this product was built (its sha256). Restore the file as it was, or give "
                f"--replace to rebuild this product with it on purpose"
                + (" (a GCP-free B edited after its run is not the pre-registered procedure's: it is refused)"
                   if variant in GCPFREE_VARIANTS else "")
                + (f". For this date the photo detection runs again (~{SEC_PER_PHOTO_DETECT:.0f} s a photo)."
                   if plan.get("era") != "live" else ""))
    alt = root_for_variant(args.output_root, variant, held)
    return (f"the product in {plan['out']} was built with another calibration ({'; '.join(diff)}): this run would "
            f"rebuild its pointing, waterlines, DEM and maps with "
            + (f"the --eo calibration (variant {variant})" if variant else "the table's calibration")
            + f" into the same folder, replacing it. One calibration per root: --output-root {alt} for this one "
              f"(then --summary --output-root on each root); give --replace only to replace this product on "
              f"purpose."
            + (f" For this date the photo detection runs again (~{SEC_PER_PHOTO_DETECT:.0f} s a photo)."
               if plan.get("era") != "live" else ""))


def keep_built_calibration(plan, args, steps, state):
    """A run that rebuilds neither the pointing nor the waterlines cannot change the calibration they were
    built with: another --eo (or none, on a variant product) is NOT APPLIED, said, and the plan carries
    the calibration on disk, so the README, provenance, DEM page and maps name the one that built them."""
    if {"pointing", "detect"} & set(steps):
        return
    built = built_calibration(plan)
    diff = calibration_differences(plan, built)
    if not diff:
        return
    text = (f"this run's calibration NOT applied ({'; '.join(diff)}): the outputs on disk were built with the "
            f"calibration shown; to apply it, rebuild from the pointing (all steps), in its own --output-root")
    warn(text)
    state["not_applied"].append(text)
    for cam, b in built.items():
        c = plan["cams"][cam]
        c["eo_file"], c["variant"] = b.get("eo_file"), b.get("variant")
        c["eo_path"] = Path(b["eo_path"]) if b.get("eo_path") and Path(b["eo_path"]).is_file() else None
        c["eo_given"] = None if str(b.get("eo_path")) == str(c.get("table_eo_path")) else b.get("eo_path")
    args.eo = {cam: c["eo_given"] for cam, c in plan["cams"].items() if c.get("eo_given")} or None
    args.calibration_variant = next((c.get("variant") for c in plan["cams"].values() if c.get("variant")), None)


def eo_option_problems(plan, args):
    """--eo must give every camera the date builds, and each file must be found. -> text ('' if fine)."""
    given = getattr(args, "eo", None) or {}
    if not given:
        return ""
    probs = []
    extra = sorted(set(given) - set(plan["cams"]))
    if extra:
        probs.append(f"--eo names {', '.join(extra)}, which {plan['date']} does not build")
    missing = [cam for cam in plan["cams"] if cam not in given]
    if missing:
        why = plan.get("omitted_why") or {}
        probs.append(f"--eo gives no calibration for {', '.join(missing)}: one variant is one calibration for "
                     f"every camera of the date (a camera is left out only when the pre-registered GCP-free procedure "
                     f"did not deliver it"
                     + ("; here: " + "; ".join(f"{c}: {w}" for c, w in sorted(why.items()) if c in missing)
                        if any(c in why for c in missing) else "") + ")")
    for cam, c in plan["cams"].items():
        if cam in given and not c.get("eo_path"):
            probs.append(f"--eo {cam}={given[cam]}: no such file (an --eo path is taken as given; the table's "
                         f"calibration is never used in its place)")
    return "; ".join(probs)


def eo_acceptance_problem(plan, args):
    """An --eo calibration whose pre-registered acceptance FAILED is not built (exit 5): the procedure's
    outcome for a failed check is 'not delivered', and the table's calibration is never used in its
    place. On a BLIND date an --eo calibration must state its acceptance ('# acceptance: ...'): a
    delivery for the holder's check is never built from a calibration whose checks were not run.
    -> text ('' when fine)."""
    bad, unstated = [], []
    for cam, c in plan["cams"].items():
        if c.get("eo_given") and c.get("eo_path"):
            acc = calibration_origin(c["eo_path"], eo_dirs(args)).get("acceptance") or {}
            if acc.get("verdict") == "FAIL":
                bad.append(f"{cam} {c['eo_file']}: '{acc['verdict']}'" + (f" ({acc['text']})" if acc.get("text") else ""))
            elif not acc and plan.get("blind"):
                unstated.append(f"{cam} {c['eo_file']}")
    if bad:
        return ("the pre-registered acceptance of the --eo calibration FAILED: " + "; ".join(bad)
                + ". Not built: a failed calibration is not delivered, and the table's calibration is never built in "
                  "its place under this variant's name")
    if unstated:
        return ("on a BLIND date an --eo calibration must state its pre-registered acceptance, and "
                + ", ".join(unstated) + " states none: run the procedure's acceptance checks first (gcpfree_eo.py "
                "writes its file only when they pass, with a first comment line 'CALIBRATION B ...: DELIVERED' or "
                "'DELIVERED-WIDER (<why>)'; another procedure, a comment line '# acceptance: PASS', '# acceptance: "
                "PASS-WIDER <the larger uncertainty>' or '# acceptance: FAIL <why>'; FAIL is never built)")
    return ""


def variant_lookalike(variant):
    """A --calibration-variant name that reads as the GCP-free procedure's B without being one of its names
    (GCPFREE_VARIANTS): 'gcpfree...' (gcpfree-oct, gcpfree2) or 'b', 'b2', 'b-...'. Its files and archive would be
    named like B's (2024-10-23_gcpfree-oct_dem.tif, ..._gcpfree-oct_calibration.tgz) with only the text inside saying
    what it is (review r1). -> text ('' when fine)."""
    v = str(variant or "")
    if not v or v in GCPFREE_VARIANTS:
        return ""
    if v.startswith("gcpfree") or re.match(r"^b\d*(-|$)", v):
        return (f"--calibration-variant {v!r} reads as the pre-registered GCP-free procedure's calibration B "
                f"({', '.join(GCPFREE_VARIANTS)}), whose names belong to the B files it delivered: another "
                f"calibration takes a name that cannot be read as B (e.g. turned-oct, handfit-2)")
    return ""


def procedure_variant_problem(plan, args):
    """The GCP-free procedure's variant names (GCPFREE_VARIANTS) belong to its B files: a B it made is built as the
    variant it names (gcpfree for the delivered B; gcpfree-aug for the second product of its one allowed second
    run), every camera from ONE of its runs; and no other calibration is built under those names, so an archive
    called B_gcpfree is always the pre-registered B. -> text ('' when fine)."""
    variant = getattr(args, "calibration_variant", None)
    given = {cam: c for cam, c in plan["cams"].items() if c.get("eo_given") and c.get("eo_path")}
    if not given:
        return ""
    look = variant_lookalike(variant)
    if look:
        return look
    procs = {cam: procedure_check(c["eo_path"], eo_dirs(args)) for cam, c in given.items()}
    probs = []
    if variant in GCPFREE_VARIANTS:
        bad = [f"{cam} {Path(given[cam]['eo_path']).name} ({'; '.join(pc.get('why') or [])})"
               for cam, pc in procs.items() if not pc.get("pinned")]
        if bad:
            probs.append(f"--calibration-variant {variant} is the name of the pre-registered GCP-free procedure's "
                         f"calibration B, and these are not (or not yet) the B files it delivered here: "
                         f"{'; '.join(bad)}. A B whose run folder is not yet committed is built once it is (one "
                         f"commit); another calibration goes under another name (it is then never called independent "
                         f"on a blind date)")
    made = {cam: pc for cam, pc in procs.items() if pc.get("pinned")}
    if len({pc.get("run_dir") for pc in made.values()}) > 1:
        probs.append("the --eo B files come from different runs of the GCP-free procedure ("
                     + ", ".join(f"{cam}: {pc.get('run_id')}" for cam, pc in sorted(made.items()))
                     + "): a B is one run's (a camera kept from the first run is written again into the second's)")
    for cam, pc in sorted(made.items()):
        if pc.get("variant") != variant:
            probs.append(f"{cam} {Path(given[cam]['eo_path']).name} was made by delivery run {pc.get('run_no')} "
                         f"({pc.get('run_id')}) of the pre-registered procedure as --calibration-variant "
                         f"{pc.get('variant')} ({GCPFREE_VARIANTS.get(pc.get('variant'), '?')}): build it as that "
                         f"variant, in its own root")
    return "; ".join(probs)


# A GCP-free B of a BLIND date is built exactly as A was, only the calibration different: the committed default
# tables, A's build options as built (provenance), the search envelope A's detection used (path and sha256) and A's
# forcing. These options name only the product, its calibration or its label, and are not compared:
B_BUILD_NOT_COMPARED = ("--date", "--output-root", "--eo", "--calibration-variant", "--synthetic", "--force-disabled")
# The one root and command a B is built with, written into its run folder at its first build (commit it): a B of
# that run is built and shared from that root only, with that command.
B_BUILD_FILE = "B_BUILD.txt"


def gcpfree_b_of(plan, args):
    """{camera: procedure_check()} of a pinned B (a GCP-free variant's --eo files that the pre-registered procedure
    delivered here) on a BLIND date, else {}."""
    if not plan.get("blind") or getattr(args, "calibration_variant", None) not in GCPFREE_VARIANTS:
        return {}
    pcs = {cam: procedure_check(c["eo_path"], eo_dirs(args)) for cam, c in plan["cams"].items()
           if c.get("eo_given") and c.get("eo_path")}
    return pcs if pcs and all(pc.get("pinned") for pc in pcs.values()) else {}


def tiff_description(path):
    """The ImageDescription (tag 270) of a classic little-endian TIFF, as asc_to_geotiff.py writes them, read from
    the file's first IFD only (no image data), else None (another layout, e.g. a BigTIFF, or not a TIFF)."""
    import struct
    try:
        with open(str(path), "rb") as f:
            head = f.read(8)
            if head[:4] != b"II*\0":
                return None
            off = struct.unpack("<I", head[4:8])[0]
            f.seek(off)
            n = struct.unpack("<H", f.read(2))[0]
            ents = f.read(12 * min(n, 512))
            for k in range(len(ents) // 12):
                tag, typ, cnt, val = struct.unpack("<HHI4s", ents[12 * k:12 * k + 12])
                if tag == 270 and typ == 2:
                    if cnt <= 4:
                        return val[:cnt].rstrip(b"\0").decode("ascii", "replace")
                    f.seek(struct.unpack("<I", val)[0])
                    return f.read(min(cnt, 1 << 16)).rstrip(b"\0").decode("ascii", "replace")
    except (OSError, struct.error, ValueError):
        return None
    return None


def product_dem_note(path):
    """The description of a GeoTIFF this script wrote (asc_to_geotiff.py --note: '... calibration ...'), else None:
    a search envelope that is one of these is a DEM of this pipeline, not a survey."""
    q = Path(str(path or ""))
    if q.suffix.lower() not in (".tif", ".tiff") or not q.is_file():
        return None
    d = tiff_description(q) or ""
    return d if ("calibration" in d and ("survey_products" in d or "EPSG:" in d or "BLIND" in d)) else None


def envelope_check(plan, cam, args, a=None):
    """On a BLIND date, what a camera's search envelope rests on (it places the detection's search band, so an
    envelope built with calibration A carries A's geometry into the waterlines). The envelope of the detection on
    disk (its stamp), else the one a build would use. For a calibration variant it counts as the date's envelope
    survey only when it is the file A's detection used (same path and sha256) and is no DEM of this pipeline.
    -> {"ok": bool, "text", "path", "sha256"}."""
    c = plan["cams"].get(cam) or {}
    st = read_stamp(plan["out"], f"detect_{cam}") or {}
    sig = (st.get("sig") or {}).get("envelope")
    path = (sig[0] if sig else None) or (str(c["envelope_path"]) if c.get("envelope_path") else None)
    sha = st.get("envelope_sha256") if sig else None
    sha = sha or file_sha(path)
    out = {"ok": True, "path": path, "sha256": sha,
           "text": f"search envelope {Path(str(path)).name if path else c.get('envelope_survey')}"
                   + (f" (sha256 {sha[:12]}...)" if sha else "")}
    if plan.get("era") == "live":
        return dict(out, text="the live detector's own envelope (not a survey)")
    note = product_dem_note(path)
    if note:
        held = any(n in note for n in ("20241023_EO", "PARTLY-CIRCULAR", "CIRCULAR"))
        return dict(out, ok=False, text=f"its search envelope {path} is a DEM of this pipeline ({note[:120]})"
                    + ("; built with the held survey's calibration (A)" if held else "")
                    + ": the detection's search band then carries that DEM's camera geometry")
    if not c.get("eo_given") or not plan.get("blind"):
        return out
    a = a or a_product(args, plan["date"], plan_variant(plan))
    ae = (a.get("envelope") or {}).get(cam)
    if not a.get("found") or not ae:
        return dict(out, ok=False, text=f"its search envelope {path} cannot be checked against A's (no A product with "
                                        f"a {cam} detection at {a.get('root')}/{plan['date']}; give --a-root): it "
                                        f"may not be the date's envelope survey")
    if str(path) != str(ae.get("path")):
        return dict(out, ok=False, text=f"its search envelope {path} is not the one A's detection used "
                                        f"({ae.get('path')})")
    if ae.get("sha256") and sha != ae["sha256"]:
        return dict(out, ok=False, text=f"its search envelope {path} (sha256 {str(sha)[:12]}...) is not the file A's "
                                        f"detection used (sha256 {ae['sha256'][:12]}...)")
    if not ae.get("sha256") and file_sig(path) != ae.get("sig"):
        return dict(out, ok=False, text=f"its search envelope {path} is not the file A's detection used (size or "
                                        f"time; A's stamp records no sha256)")
    return dict(out, text=out["text"] + ", the file A's detection used")


def b_build_record(pcs):
    """The B_BUILD.txt of a pinned B's run folder: {path, exists, root, command, options_sha256, git} ({} for
    no B)."""
    if not pcs:
        return {}
    rd = Path(str(next(iter(pcs.values())).get("run_dir")))
    f = rd / B_BUILD_FILE
    out = {"path": str(f), "exists": f.is_file()}
    if f.is_file():
        for ln in f.read_text().splitlines():
            if ":" in ln and not ln.startswith("#"):
                k, v = ln.split(":", 1)
                out.setdefault(k.strip(), v.strip())
        out["git"] = git_file_state(f)
    return out


def b_build_options(args, date):
    """The result-changing options of this build (state_argv), less those that name only the product
    (B_BUILD_NOT_COMPARED). -> {option: [values]}."""
    return {k: v for k, v in argv_options(state_argv(args, date)).items() if k not in B_BUILD_NOT_COMPARED}


def blind_b_build_problem(plan, args):
    """A GCP-free B of a BLIND date must be the pre-registered B BUILD, not only the pre-registered calibration: the
    committed default tables (surveys/survey_dates.csv, surveys.csv); A's product found (--a-root) and every
    result-changing option equal to A's as built (its provenance: the command that rebuilds it), the same search
    envelope (path and sha256 of A's detection) and the same forcing inputs; and, once its first build wrote
    B_BUILD.txt into the run folder, that root and command. Anything else is refused (exit 5): a B built with other
    options, another envelope (e.g. A's DEM) or in a second root could be chosen after comparing it with A.
    -> text ('' when fine, or not such a build)."""
    pcs = gcpfree_b_of(plan, args)
    if not pcs:
        return ""
    probs = []
    for opt, given, default in (("--config", args.config, DEFAULT_CONFIG), ("--surveys", args.surveys,
                                                                             DEFAULT_SURVEYS)):
        if os.path.realpath(str(given)) != os.path.realpath(str(default)):
            probs.append(f"{opt} {given}: a B is built with the committed {Path(default).name} of this repository")
            continue
        g = git_file_state(default)
        if not g.get("committed"):
            probs.append(f"{default} is {g.get('why')}: a B is built with the table as committed")
    a = a_product(args, plan["date"], plan_variant(plan))
    if not a.get("found") or a.get("options") is None:
        probs.append(f"no table-calibration (A) product of {plan['date']} with its build record at {a.get('root')} "
                     f"(give --a-root <A's output root>): a B is built as A was, and checked against A's record")
    else:
        mine, theirs = b_build_options(args, plan["date"]), {k: v for k, v in a["options"].items()
                                                             if k not in B_BUILD_NOT_COMPARED}
        def shown(o, k):
            return "(default)" if k not in o else (" ".join(o[k]) or "on")
        diff = [f"{k} {shown(mine, k)} (A: {shown(theirs, k)})"
                for k in sorted(set(mine) | set(theirs)) if mine.get(k) != theirs.get(k)]
        if diff:
            probs.append("its build options are not A's as built (" + "; ".join(diff) + "): a B is built with "
                         "every result-changing option as A was, only the calibration different")
        for cam in plan["cams"]:
            ec = envelope_check(plan, cam, args, a)
            if not ec["ok"]:
                probs.append(f"{cam}: {ec['text']}")
        if a.get("forcing_inputs"):
            try:
                now = {k: file_sig(q) for k, q in forcing_inputs(args).items()}
            except StepFailed:
                now = None
            if now is not None and normal(now) != normal(a["forcing_inputs"]):
                bad = sorted(k for k in set(now) | set(a["forcing_inputs"])
                             if normal(now.get(k)) != normal(a["forcing_inputs"].get(k)))
                probs.append(f"its forcing inputs are not the files A's forcing was built from ({', '.join(bad)})")
    rec = b_build_record(pcs)
    if rec.get("exists"):
        if os.path.realpath(str(args.output_root)) != rec.get("root"):
            probs.append(f"this B (run {next(iter(pcs.values())).get('run_id')}) is built in {rec.get('root')} only "
                         f"({rec['path']}: its first build): one B root, so no second B can be built and compared")
        elif rec.get("options_sha256") != options_sha256(b_build_options(args, plan["date"])):
            probs.append(f"this B's build command is fixed by its first build ({rec['path']}: "
                         f"{rec.get('command')}): rebuild it with that command")
    return "; ".join(probs)


def options_sha256(opts):
    return hashlib.sha256(json.dumps(opts, sort_keys=True).encode()).hexdigest()


def write_b_build_record(plan, args, pcs):
    """At a B's first build: its root and command, into its run folder (B_BUILD.txt; to be committed). -> path or
    None (already there)."""
    rec = b_build_record(pcs)
    if not pcs or rec.get("exists"):
        return None
    a = a_product(args, plan["date"], plan_variant(plan))
    f = Path(rec["path"])
    opts = b_build_options(args, plan["date"])
    f.write_text(
        "# B BUILD: the one root and command this run's GCP-free B is built with (survey_products.py, at its first\n"
        "# build). COMMIT THIS FILE now (git add + git commit): the B is built and shared from this root only, with\n"
        "# this command; another root or command is refused.\n"
        f"date: {plan['date']}\n"
        f"variant: {args.calibration_variant}\n"
        f"run: {next(iter(pcs.values())).get('run_id')}\n"
        f"root: {os.path.realpath(str(args.output_root))}\n"
        f"command: {cmd_text(['python3', 'survey_products.py'] + state_argv(args, plan['date']))}\n"
        f"options_sha256: {options_sha256(opts)}\n"
        f"a_root: {a.get('root')}\n"
        f"a_dem_first_built_utc: {a.get('dem_built_utc')}\n"
        + "".join(f"b_file_{cam}: {pc.get('b_sha256')}\n" for cam, pc in sorted(pcs.items()))
        + "".join(f"envelope_{cam}: {envelope_check(plan, cam, args, a).get('sha256')}\n" for cam in sorted(plan['cams']))
        + f"written_utc: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}\n")
    return f


def b_forcing_problem(plan, args):
    """After a pinned B's forcing step: its water level and waves must be A's, byte for byte (the forcing does not
    depend on the calibration). -> text or ''."""
    if not gcpfree_b_of(plan, args):
        return ""
    a = a_product(args, plan["date"], plan_variant(plan))
    bad = []
    for k, want in (a.get("forcing_sha256") or {}).items():
        q = plan["out"] / "forcing" / k
        got = sha256(q) if q.is_file() else None
        if want and got != want:
            bad.append(k)
    return (f"its forcing ({', '.join(bad)}) is not A's ({a.get('root')}/{plan['date']}/forcing): a B is built with A's "
            f"water level and waves" if bad else "")


def option_refusal(plan, args):
    """What refuses this run before anything is written (exit 5): a wrong --eo, a variant in the wrong root,
    an --eo calibration whose acceptance failed, a GCP-free B under another name than the procedure gave it,
    a --correct-pointing build outside its own root (or another build in that root).
    -> text ('' when none)."""
    return (eo_option_problems(plan, args) or variant_root_problem(plan, args) or eo_acceptance_problem(plan, args)
            or procedure_variant_problem(plan, args) or pointing_root_problem(plan, args))


def blind_independence_warnings(plan, args):
    """On a BLIND date, an --eo calibration whose geometry rests (or may rest) on the held survey: said
    loudly before the build (its README labels the holder's check PARTLY-CIRCULAR). -> list of text."""
    if not plan.get("blind"):
        return []
    out = []
    for cam, c in plan["cams"].items():
        if c.get("eo_given") and c.get("eo_path"):
            lab, why = held_check(calibration_origin(c["eo_path"], eo_dirs(args)), plan["blind"],
                                  envelope_check(plan, cam, args))
            if lab != "INDEPENDENT":
                out.append(f"{cam}: the --eo calibration (variant {c.get('variant')}) is NOT independent of the held "
                           f"survey: {why}")
    return out


def blind_leftovers(plan):
    """Outputs in this date's folder that an earlier, not blind, build made with a survey of the date: its
    comparisons and their stamps, a plan view drawn over a survey. A blind product is never built over
    them (build_date refuses). Only this product's own folder is listed. -> list of text."""
    out = plan["out"]
    found = []
    cdir = out / "compare"
    if cdir.is_dir() and any(cdir.iterdir()):
        found.append(f"{cdir}/ (comparisons)")
    sdir = out / "logs" / "stamps"
    if sdir.is_dir():
        found += [f"{p} (a comparison's stamp)" for p in sorted(sdir.glob("compare_*.json"))]
        mp = load_json(sdir / "map_plan.json") or {}
        if (mp.get("sig") or {}).get("surveys"):
            found.append(f"{map_png(plan, 'plan')} (drawn over a survey)")
    return found


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
    yet (e.g. the Oct 2024 GCPs, not delivered) is not a missing part. A window whose days measured to
    differ agree with each other and outnumber the reference's (pointing_verdicts' 'other_pointing')
    is not complete either: its product rests on the reference's pointing, the smaller group, which
    is the wrong one if the reference day is the odd one out (review 3, REF_T06: the reference day
    tilted 0.6 deg, every other day at the calibration's pointing), or mixes the two with
    --keep-moved-days; only the photos tell. When the reference is the calibration's own solve day
    (reference_is_calibration_day) its pointing is the calibration's, and only a product that mixes the
    two (or a detection not run with these results) is not complete. -> list of text."""
    out = plan["out"]
    failed = failed_cameras(plan, state)
    miss = []
    for s in plan.get("surveys") or []:
        if held_survey(s) or plan.get("blind"):
            continue                              # held for an independent check: never looked for here
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
        elif (res.get("info") or {}).get("other_pointing") and reference_is_calibration_day(res.get("info")):
            # the reference is the calibration's own day: a product on its pointing is the right one for this
            # calibration (KNOWN CAVEATS says so); only one that holds both pointings, or a detection not run with
            # these results, is not complete
            op = res["info"]["other_pointing"]
            dst = read_stamp(out, f"detect_{cam}")
            inside = [d for d in op["days"] if d not in ((dst or {}).get("skip_days") or [])
                      and d not in corrected_days(dst)]
            if dst is None or inside:
                miss.append(f"{cam}: the window holds two pointings: {', '.join(op['days'])} agree with each other, "
                            f"away from the calibration's own day {res['info'].get('reference_day')}; "
                            + ("the detection has not run with these pointing results" if dst is None else
                               f"the product holds both (the other's {', '.join(inside)} kept: with this "
                               f"calibration they are gridded with that turn as a pointing error)")
                            + f" (pointing/pointing_{cam}.csv)")
        elif (res.get("info") or {}).get("other_pointing"):
            op = res["info"]["other_pointing"]
            dst = read_stamp(out, f"detect_{cam}")
            inside = [d for d in op["days"] if d not in ((dst or {}).get("skip_days") or [])
                      and d not in corrected_days(dst)]
            miss.append(f"{cam}: the window holds two pointings: {', '.join(op['days'])} agree with each other, "
                        f"away from the reference's ({other_pointing_groups(op)}); "
                        + ("the detection has not run with these pointing results" if dst is None else
                           f"the product holds both (the other's {', '.join(inside)} kept)" if inside else
                           "the product rests on the reference's, the smaller group")
                        + f"; if the reference day is the odd one out, that is the wrong one: look at the photos "
                          f"(pointing/pointing_{cam}.csv)"
                        + (f"; {other_pointing_rebuild(op, (res.get('info') or {}).get('reference_day'), plan['date'])}"
                           if other_pointing_window(op, (res.get("info") or {}).get("reference_day")) else ""))
        if not map_png(plan, cam).exists():
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
    info = load_json(Path(f"{dem_stem(plan)}_info.json")) or {}
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


def nan_gradient(g, cell):
    """(d/drow, d/dcol) of a grid with empty (NaN) cells: a central difference where both neighbours are there, else
    a one-sided one, else NaN (np.gradient leaves NaN wherever a neighbour is empty: a thin strip of cells then has
    no slope at all)."""
    p = np.pad(np.asarray(g, float), 1, mode="constant", constant_values=np.nan)
    out = []
    for prev, nxt in ((p[:-2, 1:-1], p[2:, 1:-1]), (p[1:-1, :-2], p[1:-1, 2:])):
        d = np.where(np.isfinite(prev) & np.isfinite(nxt), (nxt - prev) / (2.0 * cell),
                     np.where(np.isfinite(nxt), (nxt - g) / cell, np.where(np.isfinite(prev), (g - prev) / cell,
                                                                           np.nan)))
        out.append(np.where(np.isfinite(g), d, np.nan))
    return out[0], out[1]


TILT_MIN_CELLS, TILT_MIN_BAND = 20, 15


def tilt_sensitivity(plan, cam, dtilt, droll, why=None):
    """
    The DEM error a pointing error of (dtilt, droll) deg would cause, per range band, from this
    date's own DEM and calibration. Each DEM cell is projected into the photo with the
    calibration, mapped back to the ground at its own elevation with the pointing changed, and
    the DEM's local slope times that horizontal shift is the elevation error. The view grazes the
    beach, so a small tilt moves a line metres along the ground (~1/sin of the depression angle).
    The local slope is np.gradient's; where that leaves too few cells (a thin strip: every cell has an
    empty neighbour), one-sided differences (nan_gradient).
    -> list of dicts (range band, cells, median shift m, median error m, median |error| m), or None; why
    (a dict, if given) gets "text": why there is no band.
    """
    from georectify import load_intrinsics, load_extrinsics, pixel_to_ground
    from view_reproject import ground_to_pixel
    from asc_to_geotiff import read_asc
    why = why if why is not None else {}
    c = plan["cams"].get(cam) or {}
    demp = Path(f"{dem_stem(plan)}_dem.asc")
    if not demp.exists() or not c.get("eo_path") or not Path(c["io_path"]).exists():
        why["text"] = ("no DEM yet" if not demp.exists() else f"no calibration or lens file for {cam}")
        return None
    g, h = read_asc(demp)
    cell = h["cellsize"]
    xll, ytop = h["xllcorner"], h["yllcorner"] + g.shape[0] * cell
    io, eo = load_intrinsics(c["io_path"]), load_extrinsics(c["eo_path"])
    rows, n_view, n_slope = [], 0, 0
    for how, grad in (("", lambda: np.gradient(g, cell)), (" (one-sided slopes)", lambda: nan_gradient(g, cell))):
        d_row, d_col = grad()
        dzde, dzdn = d_col, -d_row                   # rows run north to south
        rr, cc = np.nonzero(np.isfinite(g) & np.isfinite(dzde) & np.isfinite(dzdn))
        n_slope = max(n_slope, len(rr))
        if len(rr) < TILT_MIN_CELLS:
            continue
        E, N, Z = xll + (cc + 0.5) * cell, ytop - (rr + 0.5) * cell, g[rr, cc]
        U, V, ok = ground_to_pixel(E, N, Z, io, eo)
        n_view = int(ok.sum())
        if ok.sum() < TILT_MIN_CELLS:
            continue
        rows = tilt_bands(E, N, Z, U, V, ok, rr, cc, dzde, dzdn, io, eo, dtilt, droll, pixel_to_ground)
        if rows:
            for r in rows if how else []:
                r["slope"] = "one-sided differences where a neighbour cell is empty"
            return rows
    n_cells = int(np.isfinite(g).sum())
    why["text"] = (f"the DEM has {n_cells} cells, {n_slope} with a local slope (even from one-sided differences), "
                   f"fewer than the {TILT_MIN_CELLS} needed" if n_slope < TILT_MIN_CELLS else
                   f"the DEM has {n_cells} cells, but only {n_view} with a local slope in {cam}'s view (fewer than "
                   f"{TILT_MIN_CELLS})" if n_view < TILT_MIN_CELLS else
                   f"no range band holds {TILT_MIN_BAND} or more of the {n_view} DEM cells with a local slope in "
                   f"{cam}'s view")
    return None


def tilt_bands(E, N, Z, U, V, ok, rr, cc, dzde, dzdn, io, eo, dtilt, droll, pixel_to_ground):
    """tilt_sensitivity's per-band medians for the cells in view (ok). -> list of dicts (empty: no band)."""
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
        if m.sum() >= TILT_MIN_BAND:
            rows.append({"range_m": f"{a}-{b}", "cells": int(m.sum()),
                         "shift_m": round(float(np.median(shift[m])), 2),
                         "dem_error_m": round(float(np.median(err[m])), 3),
                         "dem_abs_error_m": round(float(np.median(np.abs(err[m]))), 3)})
    return rows


def survey_path_for(s, plan, args):
    """The survey file of row s as found (or its name when it is not here): for the label checks.
    A file not found keeps its name for the name rule, placed in the first survey folder (where it
    was looked for and is not), so no check reads a file of that name in the current folder."""
    dirs = getattr(args, "survey_dirs", None)
    try:
        hits = resolve_survey(s, plan, args)[0] if dirs is not None else []
    except Exception:                # BlindSurvey is not one: a held survey is never given a path here
        hits = []
    if hits:
        return hits[0]
    name = s.get("path") or s["name"]
    if dirs and "/" not in name and os.sep not in name and not any(ch in name for ch in "*?["):
        return Path(dirs[0]) / name
    return Path(name)


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
    text = ("EARLIER CHECK (real photos)" + (f", itself {pc['label']}: {pc['label_why']}." if pc.get("label") else ":")
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
            how_ = ("this build's own paired C = 0 sensitivity gives each camera's lines a rise of "
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
                 f"is the scale of agreement to expect). That expectation already holds the re-projection (the "
                 f"share of its setup a line keeps) and any mismatch of C on this date, listed below: do not "
                 f"subtract them again"
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
    held = [f["coef"] for f in fits if f.get("coef") is not None]     # what this build's folder holds
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
                           f"{wv['tp_currency_setup_bias_m']:+.3f} m (the DEM about {share_text()} of it).")
            if era in ("adcp", "chatham"):
                cav.append(still_water_caveat(C, wv, era))
        same = [s for s in plan.get("surveys") or []
                if survey_matches_fit(f, s, survey_path_for(s, plan, args), args)[0]]
        if same:
            cav.append(f"{', '.join(s['name'] for s in same)}: CIRCULAR for as long as C = {C} is used, because C "
                       f"was fitted to those very shots (see 'C fitted' under WHAT WAS USED). That comparison cannot test "
                       f"the setup or the overall level of the lines; it still shows their spread over the elevations "
                       f"the survey covers (WHERE, under the comparison). A build with C = 0 (or a C fitted to other "
                       f"data, --setup-fitted-to) compares INDEPENDENTLY; build it in its own folder, e.g. "
                       f"survey_products.py --date {plan['date']} --setup-coef 0 --output-root "
                       f"{root_for_c(plan['out'].parent, 0.0, held)} (in this folder it would "
                       f"replace this product, so it is refused there without --replace), and read it with --summary "
                       f"--output-root on that root"
                       + (f" (NOTE: {root_name_note(plan['out'].parent, 0.0, held)})"
                          if root_name_note(plan['out'].parent, 0.0, held) else "")
                       + ". An independent check of the whole intertidal needs another survey that fitted nothing, "
                         "e.g. a calm, low-tide RTK across it.")
        fit = SETUP_FITS.get(C) if f.get("kind") == "known" else None
        rc = (fit or {}).get("rtk_check")
        if rc:
            su = (wv.get("setup_in_window") or {})
            med = su.get("daytime_median_m") if era != "live" else None
            med = med if med is not None else su.get("median_m")
            keep = rc.get("keep_paired") or 0.7
            own = build_setup_effect(prov, [Path(s_.get("path") or "").stem for s_ in plan.get("surveys") or []])
            prior = bool(prior_checks(plan))
            # what a C between the survey-covered fits would change on THIS date: a re-projected line keeps
            # a share of the extra setup (this build's own paired C = 0 rise per camera when there is one,
            # else the shares measured so far), at this window's sqrt(Hs L0) (setup median / C)
            af = rc.get("above_floor") or {}
            cf = af.get("c_fit")
            reproj_here = ""
            if med and cf and C:
                shares = sorted(min(1.0, max(0.0, e / med)) for c, (e, what) in own.items()) \
                    if (own and era != "live") else []
                lo_k, hi_k = (shares[0], shares[-1]) if shares else SHARE_KEPT_RANGE
                phi = med / C
                lo_d, hi_d = lo_k * (min(cf) - C) * phi, hi_k * (max(cf) - C) * phi
                reproj_here = (f" Were C the ~{min(cf):.3f}-{max(cf):.3f} of the recheck's re-projected fit, the lines here (setup "
                               f"median {med:.2f} m, sqrt(Hs L0) ~{phi:.1f} m) would sit ~{lo_d:.2f}-{hi_d:.2f} m "
                               f"higher than with C = {C} ("
                               + (f"{lo_k:.2f}" if abs(hi_k - lo_k) < 0.005 else f"{lo_k:.2f}-{hi_k:.2f}")
                               + " x the extra setup, the share a "
                               "re-projected line keeps"
                               + (": this build's own paired C = 0 rise per camera)" if shares else
                                  ", as measured in the builds so far)"))
            # this build's own check of the same kind (survey_compare.py's WHERE THE SURVEY STOPS), next to
            # the recheck's numbers, which are held in SETUP_FITS and come from other frames
            own_floor = ""
            for c_ in (prov or {}).get("comparisons") or []:
                af_ = ((c_.get("headline") or {}).get("waterlines") or {}).get("above_survey_floor") or {}
                if c_.get("current") and af_.get("frames_all"):
                    own_floor = (f" THIS BUILD's own check of that kind ({c_['name']}, its own window, frames and "
                                 f"filter: WHERE THE SURVEY STOPS, above): over only the lines at least "
                                 f"{af_.get('margin_m', 0.1):.2f} m above the lowest shots, "
                                 + (f"median {af_['median']:+.3f} m on {af_['frames']} frames"
                                    + (" (too few frames: not an estimate)" if af_["frames"] < FEW_N else "")
                                    if af_.get("frames") else "no frame")
                                 + f" ({af_.get('frames_below_floor', 0)} of {af_['frames_all']} frames below the "
                                   f"lowest shot); where it differs from the recheck, read this build's.")
                    break
            cav.append(f"C = {C} is not pinned by the 2026 RTK. The Oct 2026 recheck (the station's rows of "
                       f"{rc['window']}, {rc['frames']} frames on its transects, the same frames both ways, with the "
                       f"filter of that time; numbers held in this script, SETUP_FITS, not this build's) found the "
                       f"lines {rc['without_setup_m']:.2f} m low "
                       f"without setup and {rc['with_setup_m']:.2f} m low with it (median setup applied "
                       f"{rc['setup_applied_m']:.2f} m), but both medians are biased LOW by where that survey stops "
                       f"(the transects end at +1.17 to +1.80 m): a line below a transect's lowest shot can only read "
                       f"low there, and a line that would read high at that water level lands seaward of the survey "
                       f"and is not compared."
                       + (f" Over only the lines at least {af['margins_m'][0]:.1f}-{af['margins_m'][1]:.1f} m above "
                          f"the lowest shots, that recheck (a scratch analysis, outside this repository's code) read "
                          f"the lines with C = {C} at {min(af['with_setup_m']):+.2f} to "
                          f"{max(af['with_setup_m']):+.2f} m ({min(af['frames'])}-{max(af['frames'])} frames: too "
                          f"few to be an estimate), and a re-projected fit there gave C ~{min(af['c_fit']):.3f}-"
                          f"{max(af['c_fit']):.3f} (day-block 90% ~{af['c_fit_dayblock90'][0]:.2f}-"
                          f"{af['c_fit_dayblock90'][1]:.3f})" if af else "")
                       + (f"; the within-cell slope of repeat crossings (no survey, the same analysis) gave "
                          f"{rc['c_internal_slope']:.3f} "
                          f"({rc['c_internal_slope90'][0]:.3f}-{rc['c_internal_slope90'][1]:.3f})"
                          if rc.get("c_internal_slope") else "")
                       + (f". C = {C} lies inside each of these 90% intervals" if all(
                           lo_ <= C <= hi_ for lo_, hi_ in ([af["c_fit_dayblock90"]] if af else [])
                           + ([rc["c_internal_slope90"]] if rc.get("c_internal_slope90") else []))
                          else ". Not every 90% interval here holds C = {}".format(C))
                       + "; a survey reaching the low-tide line on days with different waves would pin it."
                       + own_floor
                       + f" A line given the setup is re-projected landward onto higher "
                         f"beach, so it keeps only part of its setup (~{keep:.2f} x on the recheck's frames; the rest depends "
                         f"on the beach slope and the cross-shore distance)."
                       + reproj_here
                       + ((f". The EARLIER CHECK expectation above is made with C = {C} and already holds the "
                           f"re-projection: do not add it to that"
                           + ("; with such a C that expectation would rise by as much" if reproj_here else "") + ".")
                          if prior else ".")
                       + (" The RTK numbers were measured in the GNSS-R frame (lines at the GNSS-R level); against the "
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
    text += (" To take the question out, C can be fitted in this date's own frame, each build in its OWN "
             "--output-root (in this folder a build with another C would replace this product, so it is refused "
             "without --replace; its detection runs again, ~50 min a two-camera week): build with --setup-coef 0 "
             "--output-root <root>_c0, run "
             "dem_from_contours.py <root>_c0/<date>/waterlines/contour_points_ground_filtered.csv <scratch stem> "
             "--fit-setup --max-hs 1.5 --no-plot (repeat crossings of the same cells at different wave heights: no "
             "survey) and read what it prints. It offers a C to apply only when C = 0 and the C in use both lie "
             "outside its near-best band (C within 5 mm of its best spread) and the minimum is not flat. One week's "
             "curve is noisy: on the C = 0 lines of 29 Sep - 5 Oct 2026 (filtered, --max-hs 1.5) the band was "
             "+0.006 to +0.034 (flat), and on 29 Sep - 2 Oct alone -0.020 to +0.010 ('not distinguishable from "
             "C = 0'). Adopt a C only if it holds over weeks of different waves, by building with --setup-coef "
             "<that C> --setup-fitted-to 'none:repeat crossings, <window>' in a third root; the label stays "
             "INDEPENDENT. Otherwise keep the station's C.")
    return text


def transfer_caveat(wl, period, n_tr, n_all):
    """The README's caveat on water levels that are the Chatham gauge transferred to Marconi.
    What the extrapolation tests say about this period (the test predicting on its side of the fit
    window, its reach, RMS and bias, and whether the period lies further out) is worked out here from
    forcing.json's own record of the tests (historical_forcing.transfer_reach), not read from fields
    only a newer historical_forcing.py writes: a README rewritten after a pull (--steps maps,compare,
    a minute) from a forcing.json built before it says it too, without rebuilding the forcing step
    and so the hour of detection after it."""
    fw = wl.get("fit_window") or []
    every = n_tr == n_all
    head = (("Every water level" if every else f"{n_tr} of {n_all} water levels")
            + f" of this date {'is' if every else 'are'} the Chatham gauge 8447435 transferred to Marconi "
            f"({wl.get('method')} transfer" + (f", fitted on the ADCP {fw[0]} .. {fw[1]} UTC" if len(fw) == 2 else "")
            + f"; leave-one-week-out RMS {wl.get('cv_rms_m')} m)"
            + (f". This period lies {wl['period_vs_adcp']}" if wl.get("period_vs_adcp") else ""))
    parts = [head]
    rch = None
    if period.get("start_utc") and period.get("end_utc") and wl.get("methods") and len(fw) == 2:
        try:
            from historical_forcing import transfer_reach
            rch = transfer_reach(wl, period["start_utc"], period["end_utc"])
        except Exception:                          # an unreadable record: the fields below still say it
            rch = None
    sig = wl.get("transfer_sigma_m")
    if rch and rch["tests"]:
        side, t = rch["side"], rch["tests"]
        parts.append(f"On that side the transfer was tested only by the {'backward' if side == 'before' else 'forward'} "
                     "extrapolation test, predicting "
                     + ", ".join(f"{x['reach_days'][0]}-{x['reach_days'][1]}" for x in t)
                     + f" days {side} its fit window: RMS "
                     + ", ".join(f"{x['rms_m']:.3f} m (daytime "
                                 + ("n/a" if x.get("rms_daytime_m") is None else f"{x['rms_daytime_m']:.3f}") + ")"
                                 for x in t)
                     + ", bias " + ", ".join(f"{x['bias_m']:+.3f}" for x in t)
                     + " m (transfer minus ADCP: a positive bias reads the water, and every elevation, that much high)")
        if rch["beyond_tested"]:
            near, far = rch["beyond_days"]
            parts.append(f"This date lies {near:.0f}-{far:.0f} days {side} the fit window: FURTHER THAN TESTED, so its "
                         "error is not measured and may exceed that RMS and bias (the harbour and inlet change; the "
                         "2026 GNSS-R fit of the same transfer differs, forcing_report.txt)")
        if sig is not None and abs(float(sig) - rch["sigma_m"]) > 0.0005:
            parts.append(f"The sigma_m written in forcing/water_level.csv, {sig} m, is the pooled RMS of both tests "
                         f"(an older historical_forcing.py): this side's test gives {rch['sigma_m']:.3f} m (rebuilding "
                         "the forcing step writes it, and rebuilds everything after it)")
        else:
            parts.append(f"Each transferred level carries {sig} m as its sigma_m (forcing/water_level.csv)")
    else:
        # hours the ADCP lacks inside its record (the pooled tests stand), or a forcing.json without them
        reach = wl.get("extrapolation_reach_days") or []
        parts.append(f"Extrapolation RMS {wl.get('extrapolation_rms_m')} m (pooled over both directions); each "
                     f"transferred level carries "
                     f"{sig} m as its sigma_m"
                     + (f"; the extrapolation tests predict {' and '.join(f'{a}-{b}' for a, b in reach)} days "
                        f"beyond their fit windows" if reach else ""))
        if wl.get("beyond_tested"):
            parts.append("FURTHER THAN TESTED: the transfer's error for this date is not measured and may exceed "
                         "that RMS (the harbour and inlet change; the 2026 GNSS-R fit of the same transfer differs, "
                         "forcing_report.txt)")
    if every:
        wm = wl.get("week_mean_error")
        parts.append("The error is shared, not averaged away: every frame of the window carries the window's one "
                     "mean transfer error"
                     + (f" (over the ADCP overlap, the mean error of a whole week left out and predicted from the "
                        f"rest, 13.5-18 UTC, was {wm['min_m']:+.3f} .. {wm['max_m']:+.3f} m, RMS {wm['rms_m']:.3f} m, "
                        f"over {wm['weeks']} weeks)" if wm else ""))
    return (". ".join(parts) + ". Every line moves with it, and the DEM by about " + share_text()
            + " of it once re-projected (a higher water level also places the line further landward, on higher "
              "beach: the share of a level shift a line keeps, as of the setup in the paired C = 0 builds).")


def share_text():
    """'0.6-0.8 x' (SHARE_KEPT_RANGE): how much of a shift of the water level (or setup, or datum) the DEM keeps."""
    return f"{SHARE_KEPT_RANGE[0]}-{SHARE_KEPT_RANGE[1]} x"


def datum_size_text(d):
    """The size the ADCP datum assumption could have, from the evidence there is (forcing.json adcp_datum; the 2026
    GNSS-R comparison): a few cm, not measured."""
    a, b = d.get("offset_from_chatham_mean_m"), d.get("offset_from_msl_constant_m")
    alt = (f"the two ways of putting the ADCP on NAVD88 differ by {abs(float(a) - float(b)):.3f} m (Chatham's mean "
           f"level {float(a):+.4f} m, the tide model's MSL-to-NAVD88 constant {float(b):+.4f} m)"
           if a is not None and b is not None else "")
    return ("Its size is not measured; the evidence bounds it to a few cm: " + (alt + "; " if alt else "")
            + "on the live station the GNSS-R level sits ~0.02 m below Chatham's after its datum fix (itself +/-0.061 m, "
              "NGS OPUS with GEOID18), a comparison of mean levels of the open coast and the harbour")


def wave_conversion_caveat(wv):
    """The README's line for hours whose waves are WIS / NDBC converted to the ADCP's currency: the setup error the
    conversion causes as cross-validated (forcing.json waves.combination.setup_blend_cv, or a lone source's
    setup_cv_on_common_m), how much of it a DEM keeps, that it was tested on December hours only, and that a seasonal
    bias would be shared by every frame (review r1: one of the two largest per-frame error terms of 2024-10-23, given
    no size). None when no hour is converted."""
    used = wv.get("sources_used") or {}
    n_conv = sum(int(v) for k, v in used.items() if k != "adcp")
    n_all = sum(int(v) for v in used.values())
    if not n_conv:
        return None
    comb = wv.get("combination") or {}
    cv = comb.get("setup_blend_cv") or {}
    one = comb.get("setup_cv_on_common_m") or {}
    labels = "; ".join(v for k, v in (wv.get("source_labels") or {}).items() if k != "adcp")
    sw = (wv.get("setup_in_window") or {}).get("daytime_median_m")
    head = (f"Waves converted: {'every' if n_conv == n_all else f'{n_conv} of {n_all}'} wave hour(s) of this date "
            f"{'are' if n_conv != 1 else 'is'} {labels or 'a converted source'}, put into the ADCP's currency (the "
            f"setup coefficient's)")
    if cv.get("rms_m") is None:
        return head + (". The setup error of the conversion was not cross-validated on hours all sources share "
                       "(forcing_report.txt gives each source's own); a seasonal bias of it would be shared by every "
                       "frame of the window.")
    lo, hi = SHARE_KEPT_RANGE
    cp = comb.get("common_period") or ["?", "?"]
    return (head + f". The setup they give was cross-validated against the ADCP's own setup on the "
            f"{comb.get('common_hours')} h all sources share ({cp[0]} .. {cp[1]} UTC) only: RMS {cv['rms_m']:.3f} m, "
            f"bias {cv.get('bias_m', 0.0):+.3f} m"
            + (" (" + ", ".join(f"{k} alone {v:.3f} m" for k, v in sorted(one.items())) + ")" if one else "")
            + (f", against a setup of {sw:.2f} m (daytime median) in this window" if sw is not None else "")
            + f". A line keeps about {lo}-{hi} x a setup error once re-projected, so ~{cv['rms_m'] * lo:.2f}-"
              f"{cv['rms_m'] * hi:.2f} m RMS per frame. Those hours are of December only: a seasonal bias of the "
              f"conversion (e.g. October swell directions, not tested) would be shared by every frame of the window, "
              f"not averaged away (forcing_report.txt, 'Combination for hours without the ADCP').")


def mount_solutions(cam, station, cal_dir, leave_out_dates=()):
    """The GCP solutions of camera `cam` on mount `station` in calibration_history.csv whose files are in
    CALIBRATION_ORIGINS, without those solved on leave_out_dates (a held survey's: their numbers are not
    parsed). -> [(file, date, x, y, z)], by date."""
    sols = []
    try:
        with open(Path(cal_dir) / "calibration_history.csv", newline="") as f:
            for r in csv.DictReader(f):
                reg = CALIBRATION_ORIGINS.get(r.get("file") or "")
                if not reg or reg[1] in leave_out_dates or r.get("camera") != cam or r.get("setup") != station:
                    continue
                try:
                    sols.append((r["file"], reg[1], float(r["x"]), float(r["y"]), float(r["z"])))
                except (KeyError, TypeError, ValueError):
                    continue
    except OSError:
        return []
    return sorted(sols, key=lambda t: t[1])


def calibration_caveats(plan, args, prov):
    """A calibration turned from another by matching photos keeps that one's camera position: said, with
    how far the mount's GCP solutions put the camera (calibration_history.csv; on a blind date the
    solution from the held survey is left out, unread), and its pre-registered acceptance."""
    cav = []
    held_dates = {plan["blind"]["survey_date"]} if plan.get("blind") else set()
    cal_dir = getattr(args, "calibration", None) or HERE / "calibration"
    for cam, v in (prov.get("cameras") or {}).items():
        o = (v.get("calibration") or {}).get("origin") or {}
        if len(o.get("chain") or []) < 2:            # not turned from another calibration
            continue
        pos = o.get("position_from")
        sols = mount_solutions(cam, v.get("station"), cal_dir, held_dates)
        spread = ""
        if len(sols) > 1:
            zs, xy = [t[4] for t in sols], [(t[2], t[3]) for t in sols]
            dxy = max(float(np.hypot(a[0] - b[0], a[1] - b[1])) for a in xy for b in xy)
            spread = (". The mount did not move, yet its GCP solutions "
                      + ", ".join(f"{t[1]} ({t[0]})" for t in sols)
                      + f" put this camera {max(zs) - min(zs):.2f} m apart in height and {dxy:.2f} m horizontally "
                        f"(calibration_history.csv"
                      + ("; the solution from the held survey is left out: nothing of it is used by this build"
                         if held_dates else "")
                      + "): each solution trades its position against its angles over its own GCPs, so this DEM "
                        f"carries the position error of {pos} away from where its GCPs lay")
        acc = o.get("acceptance") or {}
        acc_t = ({"PASS": "its pre-registered acceptance checks passed (PASS"
                          + (f": {acc['text']}" if acc.get("text") else "") + ")",
                  "PASS-WIDER": f"its pre-registered acceptance checks passed only with a LARGER UNCERTAINTY, stated "
                                f"with it: {acc.get('text') or '(no uncertainty stated)'}",
                  "FAIL": f"its pre-registered acceptance FAILED ({acc.get('text') or 'no reason stated'}): it should "
                          f"not have been built"}.get(acc.get("verdict"), "")
                 if acc else "it states no pre-registered acceptance (no '# acceptance: ...' line, no "
                             "gcpfree_eo.py 'CALIBRATION B ...: DELIVERED' line): its checks are not shown here")
        cav.append(f"{cam} calibration {o.get('file')}: made without ground control of its own; its angles were "
                   f"turned from {' and '.join(o.get('refs') or [])} by matching photos"
                   + (f" (the whole chain: {' <- '.join(o.get('chain') or [])})"
                      if len(o.get("chain") or []) > len(o.get("refs") or []) + 1 else "")
                   + f", and its camera position "
                   f"(x, y, z) is that of {pos}{spread}. The turn adds its own error (the matching residual); the "
                   f"pointing step's sea-horizon check of this product's photos against this calibration (WHAT WAS "
                   f"USED: 'sea horizon vs the calibration') is independent of the matches; {acc_t}."
                   + (f" Its STATED UNCERTAINTY (the procedure's pre-registered rule): {uncertainty_terms(o['procedure'])}."
                      if (o.get("procedure") or {}).get("pinned") and uncertainty_terms(o["procedure"]) else "")
                   + (f" It is NOT the pre-registered procedure's: {'; '.join(o['procedure'].get('why') or [])}."
                      if (o.get("procedure") or {}).get("format") and not o["procedure"].get("pinned") else ""))
    for cam, v in sorted((plan.get("not_delivered") or {}).items()):
        cav.append(f"{cam} is NOT in this product: the pre-registered GCP-free procedure did not deliver it "
                   f"({v.get('verdict')}; {'; '.join(v.get('reasons') or []) or 'no reason given'}; report "
                   f"{v.get('report')}). The table's calibration is never used in its place.")
    return cav


def pointing_correction_caveat(cam, pcr, byday=None):
    """The README's line for a --correct-pointing detection, from its stamp's record (correct_pointing_rows): each day
    mapped with its own pointing, with the change, the match and how far its lines moved; the days measured to differ
    (or unconfirmed) that a gate kept from a correction, and why; what the correction is against and what it is not."""
    byday = byday or {}
    ru = pcr.get("rules") or {}
    parts = []
    for d, v in sorted((pcr.get("days") or {}).items()):
        a, t, r = v["d_deg"]
        mv = (f", its {v.get('georectified')} points moved {v['moved_m_median']:.2f} m median ({v['moved_m_max']:.2f} m "
              f"max) from where the calibration put them" if v.get("moved_m_median") is not None else
              f", {v.get('georectified', 0)} points")
        chk = v.get("day_check")
        held = ("" if chk is None else
                f" (the {v.get('probe_time_utc') or '?'} UTC frame; held through the day: its "
                f"{', '.join((f.get('time_utc') or '?') for f in chk)} UTC frames read it within "
                f"{v.get('day_spread_deg') or 0:.3f} deg)" if chk else
                f" (the {v.get('probe_time_utc') or '?'} UTC frame, the day's only photo)")
        parts.append(f"{d} ({(byday.get(d) or {}).get('verdict') or v.get('verdict') or '?'}): {a:+.3f}/{t:+.3f}/"
                     f"{r:+.3f} deg pan/tilt/roll on {v.get('inliers')} inliers, turn_sigma {v.get('turn_sigma')}"
                     + held + mv)
    left = [f"{d} ({'; '.join(w)})" for d, w in sorted((pcr.get("not_corrected") or {}).items())
            if (byday.get(d) or {}).get("verdict") in LEFT_OUT_VERDICTS]
    return (f"{cam}: POINTING CORRECTED PER DAY (--correct-pointing: an evaluation product, in its own root): "
            + ("; ".join(parts) if parts else "no day passed the gates")
            + f". Each such day was detected like the others and its waterlines mapped with the calibration turned by "
              f"the change matched features measured that day against the reference frame "
              f"{pcr.get('reference_frame') or '?'} ({pcr.get('reference_day') or '?'}) (pointing_daily.py; the EO "
              f"files in waterlines/{cam}/pointing_corrected_eo/); the default product leaves a day measured to differ "
              f"out instead, and grids a smaller change uncorrected. Gates: {ru.get('min_inliers')} or more inliers, "
              f"turn_sigma {ru.get('min_turn_sigma', 0):g} or more, inliers spread {ru.get('min_spread_px', 0):g} px or "
              f"more, a turn of at most {ru.get('max_turn_deg', 0):g} deg"
            + (f" and of more than {ru['min_turn_deg']:g} deg (the match's own error: a day at rest is gridded as by "
               f"default)" if "min_turn_deg" in ru else "")
            + (f", the same turn within {ru['max_day_spread_deg']:g} deg on the day's first and last photos (the "
               f"pointing held through the day)" if "max_day_spread_deg" in ru else "")
            + ", not UNCONFIRMED, not the reference day"
            + (f". Measured to differ or unconfirmed and NOT corrected (left out or kept as the lines above say): "
               f"{'; '.join(left)}" if left else "")
            + ". Only the georectification is corrected: the detector searched each day with the calibration's "
              "envelope or crop (padded ~1-2.7 deg of tilt), so a day tilted ~1 deg or more may have lines outside "
              "it. The window's constant offset against the calibration (the sea horizon's, above) is not corrected "
              "(docs/VALIDATION.md 4.6b).")


def small_change_caveat(plan, cam, rows, skip=()):
    """The README's line for a camera's KEPT 'small change' days (pointing_verdicts: matched features measured a turn
    of POINT_WARN_DEG - POINT_FAIL_DEG with tilt and roll within RUN_TOL_DEG, or an ambiguous match): each day's
    measured turn and what that tilt/roll does to its lines and to the DEM where they are gridded, per range band
    (tilt_sensitivity, from this date's DEM and calibration). '' when there is none. Such a day is kept with the
    calibration's pointing, so its frames carry that turn as a pointing error (review r1: up to ~0.4 m at range)."""
    days = [r for r in rows if r.get("verdict") == "small change" and r.get("date") not in set(skip or [])]
    if not days:
        return ""
    parts = []
    for r in days:
        if r.get("features") == "ambiguous":
            rup = [abs(float(x)) for x in str(r.get("runner_up_d") or "").split("/") if x.strip()]
            parts.append(f"{r['date']}: an ambiguous match ({_inliers(r)} features at zero turn and a second consensus "
                         f"of {r.get('runner_up_inliers')} at {r.get('runner_up_d')} deg: a turn of up to "
                         f"{max(rup) if rup else 0.0:.2f} deg)")
            continue
        try:
            a, t, rl = (float(r[k]) for k in ("d_azimuth_deg", "d_tilt_deg", "d_roll_deg"))
        except (KeyError, TypeError, ValueError):
            parts.append(f"{r['date']}: {r.get('why') or 'a small change'}")
            continue
        txt = f"{r['date']}: {a:+.2f}/{t:+.2f}/{rl:+.2f} deg pan/tilt/roll on {_inliers(r)} features"
        sens = tilt_sensitivity(plan, cam, t, rl) if max(abs(t), abs(rl)) > 0 else None
        if sens:
            txt += (" (its tilt/roll moves that day's lines, and the DEM where they are gridded, by "
                    + ", ".join(f"~{x['dem_abs_error_m']:.2f} m at {x['range_m']} m" for x in sens) + ")")
        parts.append(txt)
    return (f"{cam}: kept with a small change of pointing, measured by matched features against the reference (a turn "
            f"of {POINT_WARN_DEG}-{POINT_FAIL_DEG} deg, tilt and roll within the step's tolerance; a match reads them to "
            f"~{MATCH_ERR_DEG:g} deg): " + "; ".join(parts) + ". These frames are gridded with the calibration's "
            f"pointing, so each carries its day's turn as a pointing error (median over the cells of each range band, "
            f"from this date's DEM and calibration; a pan moves the lines along the beach and is not in these "
            f"numbers); the DEM holds it in proportion to that day's share of the frames (pointing/pointing_{cam}.csv).")


def collect_caveats(plan, args, state, forcing, prov):
    cav = []
    era = plan["era"]
    out = plan["out"]
    wl = (forcing or {}).get("water_level", {})
    wv = (forcing or {}).get("waves", {})
    fits = (prov.get("setup") or {}).get("fits") or setup_in_use(plan, args)[0]
    if era in ("adcp", "chatham"):
        d = wl.get("adcp_datum") or {}
        srcs = wl.get("sources_used") or {}
        n_tr = sum(int(v) for k, v in srcs.items() if "chatham" in str(k))
        n_all = sum(int(v) for v in srcs.values())
        cav.append("ADCP datum: the ADCP's NAVD88 offset was derived from the Chatham gauge's mean level over "
                   "the deployment (adcp_to_navd88.py), i.e. Marconi's mean level is assumed equal to "
                   f"Chatham's{' (' + str(d.get('source')) + ')' if d.get('source') else ''}. Every line "
                   "shares this datum error (the DEM about " + share_text() + " of it, once re-projected)"
                   + (" (no ADCP level is used here, but the Chatham transfer was fitted to the ADCP's NAVD88 "
                      "levels, so it carries the same datum)" if n_all and n_tr == n_all else "")
                   + ". " + datum_size_text(d) + ".")
        if n_tr:
            cav.append(transfer_caveat(wl, (forcing or {}).get("period") or {}, n_tr, n_all))
        for c in wv.get("caveats") or []:
            cav.append("Waves: " + c)
        wc_ = wave_conversion_caveat(wv)
        if wc_:
            cav.append(wc_)
    else:
        had = [c for c in plan["cams"] if (((read_stamp(out, f"detect_{c}") or {}).get("setup") or {}).get("had_setup"))]
        C = fits[0]["coef"] if fits else args.setup_coef
        if had:
            diffs = [((read_stamp(out, f"detect_{c}") or {}).get("setup") or {}).get("max_diff_existing") or 0.0
                     for c in had]
            cav.append("Live era: the station's rows already carried a setup (the cron applies it since 6 Oct 2026); "
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
                   "heights agree to 2 mm; +/-0.061 m, mostly the geoid model). Every line shares that offset (the "
                   "DEM about " + share_text() + " of it, once re-projected). It largely cancels against RTK shots on "
                   "NAVD88 (GEOID18), the same geoid model, but not in the absolute NAVD88 heights of the GeoTIFFs.")
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
    # the DEM's day test, from its own record: a camera-day it could not judge or kept untested, and a
    # camera whose day offsets run one way with time (its WARNINGs, as dem_from_contours printed them)
    flt_ = (load_json(Path(f"{dem_stem(plan)}_info.json")) or {}).get("filters") or {}
    if flt_.get("max_day_offset") is not None:
        from dem_figure import day_test_warnings
        for ln in day_test_warnings(flt_):
            cav.append("DEM day test " + ln + (" (WHAT WAS USED: each camera-day it left out, how and why)"
                                               if flt_.get("camera_days_rejected") else
                                               " (it left out no camera-day)"))
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
                not map_png(plan, cam).exists():
            cav.append(f"{cam}: no waterline map on the photos: {map_missing_reason(plan, state, cam)} (the "
                       f"plan-view map does not need the photos).")
    for cam, res in pointing.items():
        info = res.get("info", {})
        # what the detection did is in its own stamp (the days it left out, and whether it ran with
        # --keep-moved-days), never this run's options: a later --steps compare or pointing run without
        # the option rewrites this README too
        dst = read_stamp(out, f"detect_{cam}")
        skip = sorted((dst or {}).get("skip_days") or [])
        byday = {r.get("date"): r for r in res.get("rows") or []}

        def with_why(ds):
            return "; ".join(f"{d_} ({byday[d_].get('verdict')}: {byday[d_]['why']})" if byday.get(d_, {}).get("why")
                             else d_ for d_ in ds)
        diff_ = list(res.get("different_days") or [])
        # a detection made with other pointing results (the pointing step was run again on its own since,
        # e.g. after a change of its rule): the days it left out are not those these results leave out, so
        # neither its list nor these results' 'kept' days describe the product. A stamp of an older build
        # did not record --keep-moved-days: whether it ran with it is never inferred (review 3: a deployed
        # stamp that left out no day, under results that now leave days out, read 'KEPT in the detection
        # (--keep-moved-days)'); its days are compared with these results' defaults, and a difference has
        # either reason
        kmd_ = (dst or {}).get("keep_moved_days")
        # a --correct-pointing detection (its stamp's record; its root holds no other) is judged by that option
        cp_ = "pointing_correction" in (dst or {})
        dflt = argparse.Namespace(correct_pointing=True) if cp_ else None
        stale = dst is not None and (skip != skip_days_of(res, argparse.Namespace(keep_moved_days=bool(kmd_),
                                                                                  correct_pointing=cp_))
                                     or (cp_ and corrected_days(dst) != sorted(pointing_corrections(res)[0])))
        if dst is None or stale:
            if stale:
                cav.append(f"{cam}: the detection left out " + (", ".join(skip) if skip else "no day")
                           + (f" and corrected {', '.join(corrected_days(dst)) or 'no day'}" if cp_ else "")
                           + ", where these pointing results leave out "
                           + (", ".join(skip_days_of(res, dflt)) or "no day")
                           + (f" and correct {', '.join(sorted(pointing_corrections(res)[0])) or 'no day'}" if cp_
                              else "")
                           # --keep-moved-days leaves out no day: it can explain a detection that left out
                           # none, never one that left days out (review 4)
                           + (" by default: it ran with other pointing results (the pointing step was run again "
                              "since) or with --keep-moved-days (its record, of an older build, does not say which)"
                              if kmd_ is None and not skip else
                              ": it ran with other pointing results (the pointing step was run again since)")
                           + "; run the detection again (--steps detect and the steps after it) to apply these.")
            if diff_:
                cav.append(f"{cam}: days with a different or unconfirmed pointing: {with_why(diff_)}. The detection "
                           f"has not run with these pointing results"
                           + (" (not requested, or it failed)" if dst is None else "")
                           + ": it leaves such days out by default"
                           + ("" if res.get("leave_out", True) else
                              ", but these results come from the rule before Oct 2026, whose 'most days differ' guard "
                              "kept every day (see the warning); a detection now needs the pointing step run again")
                           + "; --keep-moved-days keeps them.")
        else:
            if skip:
                cav.append(f"{cam}: days left out of the detection for their pointing: {with_why(skip)}.")
            kept_ = [d_ for d_ in diff_ if d_ not in skip and d_ not in corrected_days(dst)]
            if kept_:
                cav.append(f"{cam}: days with a different or unconfirmed pointing KEPT in the detection ("
                           f"{kept_reason(res, args, dst)}): {with_why(kept_)}.")
            if cp_:
                cav.append(pointing_correction_caveat(cam, dst["pointing_correction"], byday))
        ran_ = dst is not None and not stale
        note_ = (pointing_majority_note(res, args, skip=skip if ran_ else skip_days_of(res, dflt), ran=ran_)
                 if res.get("rows") else None)
        if note_:
            cav.append(f"{cam}: {note_}.")
        over = [d_ for d_, r_ in sorted(byday.items()) if r_.get("verdict") not in LEFT_OUT_VERDICTS
                and d_ in (res.get("horizon_overruled_days") or info.get("horizon_overruled_days") or [])]
        if over and ran_:
            cav.append(f"{cam}: kept, although the sea-horizon check alone would have left them out (matched "
                       f"features measured no turn, or not the one the horizon read): {with_why(over)}.")
        elif over:
            # the detection on disk did not run with these results: what it did is said above, and these days
            # are kept only once it runs again (review 4: a pointing-only rerun on a deployed build listed days
            # the detection on disk had left out as 'kept')
            gone = [d_ for d_ in over if d_ in skip] if dst is not None else []
            cav.append(f"{cam}: these pointing results keep {', '.join(over)}, although the sea-horizon check alone "
                       f"would leave them out (matched features measured no turn, or not the one the horizon read)"
                       + (f"; the detection on disk left out {', '.join(gone)}, so they are kept once it runs again"
                          if gone else "")
                       + f": {with_why(over)}.")
        small_ = small_change_caveat(plan, cam, res.get("rows") or [], (skip + corrected_days(dst)) if ran_ else
                                     skip_days_of(res, dflt) + (sorted(pointing_corrections(res)[0]) if cp_ else []))
        if small_:
            cav.append(small_)
        bun_ = [r_ for r_ in res.get("rows") or [] if r_.get("verdict") not in LEFT_OUT_VERDICTS and bunched(r_)
                and features_trusted(r_)]
        if bun_:
            cav.append(f"{cam}: kept on matched features whose inliers are bunched in a small part of the frame ("
                       + ", ".join(f"{r_['date']}: {_inliers(r_)} inliers, spread {spread_of(r_)[0]:.0f} x "
                                   f"{spread_of(r_)[1]:.0f} px" for r_ in bun_)
                       + f"): every match the ~{MATCH_ERR_DEG:g} deg of a match rests on spread "
                         f"{MATCH_SPREAD_TESTED_PX} px or more, so these days' tilt and roll may be read worse than "
                         f"that (pointing/pointing_{cam}.csv, inlier_spread_px).")
        if info.get("fitted_horizon_days") and not info.get("clear_horizon_days"):
            px_ = info.get("fitted_horizon_offset_px")
            cav.append(f"{cam}: no day of the window showed a clear sea horizon (the edge read as the horizon was "
                       f"faint, inverted, matched by a rival edge or at the end of its search, as haze, spray or foam "
                       f"is; a rendered horizon can be too): its readings (median "
                       f"{info.get('window_horizon_dtilt_deg', 0.0):+.2f}/{info.get('window_horizon_droll_deg', 0.0):+.2f}"
                       f" deg tilt/roll against the calibration"
                       + (f", {px_:.1f} px" if px_ is not None else "")
                       + f" over {info['fitted_horizon_days']} fitted day(s)) may not be the horizon, so no offset "
                       f"against the calibration is reported; they were compared with each other only.")
        eo_ = info.get("calibration") or plan["cams"].get(cam, {}).get("eo_file")
        off = pointing_offset(info, eo_)
        if off:
            dt, dr = off[0], off[1]
            sens_why = {}
            sens = tilt_sensitivity(plan, cam, dt, dr, why=sens_why)
            if sens:
                off = pointing_offset(info, eo_, rate=False) or off
            text = f"{cam}: {off[2]}"
            if not sens:                     # the rough rate's promise of a per-band table is kept only when it is made
                text = text.replace("the README gives it per range band from this date's own DEM, roll included",
                                    f"not given per range band from this date's DEM: {sens_why.get('text') or 'no DEM'}")
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
        if info.get("horizon_anchor_note"):
            cav.append(f"{cam}: {info['horizon_anchor_note']}.")
        rows_ = res.get("rows") or []
        hz_only = [r["date"] for r in rows_ if r.get("verdict") == "horizon only"]
        unchk = [r["date"] for r in rows_ if r.get("verdict") == "unchecked"]
        # a day whose features matched, but within their own noise (blind_match), counts as no reading: said so
        # (never 'too few matched features'), or the inlier count in the table looks like a contradiction
        if rows_ and len(unchk) == len(rows_):
            cav.append(f"{cam}: nothing measured the pointing on any day ({no_feature_reading(unchk, rows_)}, and "
                       f"no clear sea horizon): a turn or tilt of the camera would not have been seen.")
        elif rows_ and len(hz_only) + len(unchk) == len(rows_):
            cav.append(f"{cam}: the pan (azimuth) could not be measured on any day ("
                       f"{no_feature_reading(hz_only + unchk, rows_)}): the days were checked by the sea horizon "
                       f"only (tilt and roll), so a turn of the camera would not have been seen"
                       + (f"; on {', '.join(unchk)} not even that (no clear sea horizon)." if unchk else "."))
        else:
            if hz_only:
                cav.append(f"{cam}: pan not measured on {', '.join(hz_only)} (horizon only: "
                           f"{no_feature_reading(hz_only, rows_)}).")
            if unchk:
                cav.append(f"{cam}: pointing not measured on {', '.join(unchk)} (unchecked: "
                           f"{no_feature_reading(unchk, rows_)}, and no clear sea horizon); kept.")
        if info.get("reference_kind", "").startswith("the day nearest") and \
                not info["reference_kind"].startswith(TURNED_REF_KIND):
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
    cav += calibration_caveats(plan, args, prov)
    epsg = ((dem_settings_built(out) or {}).get("geotiff_epsg")      # what the GeoTIFFs on disk carry
            or getattr(args, "geotiff_epsg", GEOTIFF_EPSG))
    cav.append(f"Coordinates: the grids are in the frame of the calibration's GCPs, ASSUMED to be that of the "
               f"surveys, NAD83(2011) / UTM zone 19N (EPSG:6348): the GCP files state no CRS. Were the GCPs in "
               f"WGS 84 / ITRF instead, the grid would sit ~1-1.5 m off the surveys, up to ~0.1-0.15 m of "
               f"elevation where that shift runs across a 1:10 foreshore. The GeoTIFFs are tagged EPSG:{epsg}"
               + (" (WGS 84 / UTM 19N), a NOMINAL tag: read them as NAD83(2011) (WGS 84 differs by "
                  "~1-1.5 m here; nothing was transformed; --geotiff-epsg 6348, the default since Oct 2026, writes "
                  "the true code)" if epsg == 32619 else
                  ", the frame assumed above: a GIS places them with the surveys" if epsg == 6348 else "")
               + (". Nothing was shifted: a survey is to be compared at the DEM's coordinates as they are (this blind "
                  "product made no comparison)." if plan.get("blind") else
                  ". The surveys are sampled at the DEM's coordinates with no shift."))
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
    # the storm and high-wave caveats share one note on what the per-day rows can say: given once, after them
    sw = storm_caveats(plan, args, prov, first_used, last_used, note=False) + \
        wave_event_caveats(plan, prov, first_used, last_used, note=False)
    if len(sw) == 1:
        sw[-1] += " " + PER_DAY_ROWS_NOTE
    elif sw:
        sw.append("STORM and HIGH WAVES (the caveats above): " + PER_DAY_ROWS_NOTE)
    cav += sw
    return cav


# A setup coefficient fitted to NO survey, from repeat crossings of the same DEM cells, each line
# re-projected at still water + setup (the pointing argument below applies it to re-projected lines):
# the Oct 2026 recheck of the filtered 29 Sep - 5 Oct 2026 lines (a scratch analysis made while
# reviewing C, not code in this repository) gave 0.029 for the spread that
# dem_from_contours.py --fit-setup minimises (a flat minimum, 90% ~0.00-0.04; --fit-setup itself, which
# does not re-project, prints a near-best band of +0.006 to +0.034 there on the C = 0 lines, --max-hs 1.5)
# and 0.043 (90% 0.026-0.053) for the within-cell slope
# of elevation on sqrt(Hs*L0). The range is the two estimates, not an interval of either. A
# constant pointing error moves every frame's line alike, so it does not change how well repeat
# crossings agree: this C is blind to it, unlike the C fitted to the 2026 RTK (which absorbs any level
# error at the transects, of the setup or of the pointing). Stockdon et al. (2006): setup = 0.35 x
# foreshore slope x sqrt(H0 L0), so C ~0.035-0.045 on a 1:10-1:8 beach face.
SETUP_C_INTERNAL = (0.029, 0.043)
SETUP_C_INTERNAL_HOW = ("Oct 2026 recheck of the 29 Sep - 5 Oct 2026 lines, each re-projected, outside this "
                        "repository's code: the repeat-crossing spread that dem_from_contours.py --fit-setup "
                        "minimises, 0.029 (a flat minimum, 90% ~0.00-0.04; --fit-setup itself, not re-projecting, "
                        "prints a near-best band of +0.006 to +0.034 on the C = 0 lines), and the within-cell slope "
                        "of elevation on sqrt(Hs L0), 0.043 (0.026-0.053)")
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
                        x[k] = float(wave_scale(hs, tp))
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
                f"DEM cells (no survey) give C = {SETUP_C_INTERNAL[0]}-{SETUP_C_INTERNAL[1]} ({SETUP_C_INTERNAL_HOW}), "
                f"blind to a constant pointing error, which moves every frame's "
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
        # the C = 0 lines sit low on the beach: those below where the survey stops can only read low there
        fl = (wl.get("above_survey_floor") if no_setup else (wl.get("without_setup") or {}).get("above_survey_floor")) or {}
        floor_txt = ""
        if fl.get("frames_below_floor"):
            floor_txt = (f". But {fl['frames_below_floor']} of {fl['frames_all']} of those C = 0 frames lie below "
                         f"their transect's lowest shot, where a line can only read low (one that would read high "
                         f"lands seaward of the survey): the lines without setup read low by construction there, so "
                         f"what is left, and the share of the offset it allows, is OVERSTATED")
        return (head + where + verdict + floor_txt
                + ". Read the shifts above as what the offset WOULD do if it were all pointing")
    return None


def dem_gridded_frames(plan):
    """
    Which waterline frames the DEM gridded, from the file it read (waterline_file) and its own
    record (dem/<date>_info.json, dem/<date>_<variant>_info.json with another calibration (dem_stem):
    its window, its wave filter (a frame with offshore Hs > max_hs is left out; one with no wave record
    is kept, as there) and the camera-days its day test rejected, local days as it counts them, each
    with its own reason (dem_figure.day_test_reason)), as survey_compare.dem_frame_selection.
    -> {"gridded": {camera: {frame: (epoch, utc day)}}, "left_out": {camera: {frame: (epoch, utc day,
    why)}}, "max_hs": float or None}, or None without the DEM's record or the file.
    """
    info = load_json(Path(f"{dem_stem(plan)}_info.json"))        # a variant's DEM is <date>_<variant>
    if not info:
        return None
    filt = info.get("filters") or {}
    mx = filt.get("max_hs")
    rej = set(filt.get("camera_days_rejected") or [])
    w0, w1 = info.get("window_start"), info.get("window_end")
    res = {"gridded": {}, "left_out": {}, "max_hs": float(mx) if mx is not None else None}
    seen = set()
    try:
        with open(waterline_file(plan, quiet_=True), newline="") as f:
            for r in csv.DictReader(f):
                k, cam = r.get("source_file", ""), r.get("camera", "?")
                if (cam, k) in seen:
                    continue
                seen.add((cam, k))
                day = (r.get("capture_time_utc") or "")[:10]
                if (w0 and day < w0) or (w1 and day > w1):
                    continue
                e = _num(r.get("capture_epoch"))
                if e is None:
                    e = epoch_of(k)
                e = float(e) if e is not None else float("nan")
                hs = _num(r.get("offshore_hs_m"))
                why = None
                if mx is not None and hs is not None and hs > float(mx):
                    why = f"offshore Hs > {float(mx):g} m (the DEM's wave filter)"
                elif rej and e == e:
                    from dem_from_contours import local_day
                    if f"{cam} {local_day(e)}" in rej:
                        # by the day's own reason (a day the iterative test left out within the limit
                        # of the kept days is not 'off the rest'), as the DEM page words it
                        from dem_figure import day_test_reason
                        why = day_test_reason(filt, f"{cam} {local_day(e)}") + " (the DEM's day test)"
                if why:
                    res["left_out"].setdefault(cam, {})[k] = (e, day, why)
                else:
                    res["gridded"].setdefault(cam, {})[k] = (e, day)
    except OSError:
        return None
    return res


def dem_days_text(plan, fpd):
    """'; the DEM grids N of their M frames, on FIRST .. LAST (k days)' (dem_gridded_frames), when the
    DEM's own record shows it gridded fewer frames or days than the waterlines hold (fpd: their
    frames_per_day; the DEM's window, its wave filter, its day test); '' otherwise, or without that
    record."""
    dg = dem_gridded_frames(plan)
    days = sorted({d for v in (fpd or {}).values() for d in v})
    if not dg or not days:
        return ""
    n_all = sum(n for v in fpd.values() for n in v.values())
    g = [v for c in dg["gridded"].values() for v in c.values()]
    gd = sorted({d for _, d in g})
    if not g:
        return "; the DEM grids none of them"
    if len(g) >= n_all and gd == days:
        return ""
    return (f"; the DEM grids {len(g)} of their {n_all} frames, on "
            + (f"{gd[0]} .. {gd[-1]} ({len(gd)} days)" if len(gd) > 1 else gd[0]))


def storm_caveats(plan, args, prov, first_used, last_used, note=True):
    """A caveat per known storm (STORMS) that reaches into the window or ended within 3 days of it:
    the window's days of the storm's tail before the survey and, from the DEM's own record, whether
    the DEM grids frames of them, or why not (left out by the pointing check, by the DEM's wave
    filter or day test, no waterline); what changed and where, and whether a point survey lies on the
    part that changed most. note: end with PER_DAY_ROWS_NOTE (a README gives it once for all the
    storm and high-wave caveats)."""
    cav = []
    out = plan["out"]
    fpd = prov.get("frames_per_day") or {}
    days = sorted({d for v in fpd.values() for d in v}) or days_between(first_used, last_used)
    dg = dem_gridded_frames(plan)
    cams = sorted(set(plan.get("cams") or {}) | set(fpd))
    for sev in STORMS:
        a, b = date_cls.fromisoformat(sev["first"]), date_cls.fromisoformat(sev["last"])
        w0, w1 = date_cls.fromisoformat(plan["first"]), date_cls.fromisoformat(plan["last"])
        if b < w0 - timedelta(days=3) or a > w1:
            continue
        where = (f"reaches into this window ({plan['first']} .. {plan['last']})" if b >= w0 else
                 f"ended {(w0 - b).days} day(s) before this window")
        sdays = sorted({x["survey_date"] for x in plan.get("surveys") or [] if x.get("survey_date")})
        sd = sdays[0] if sdays else plan["date"]
        if dg is None:
            # no DEM record (the DEM step has not run or failed): what the waterlines hold
            tail = [d for d in days if sev["first"] <= d < sd]
            tail_txt = (f"; {', '.join(tail)} (the storm's tail, before the {sd} survey) "
                        f"{'is' if len(tail) == 1 else 'are'} in the window (no DEM record says which of their "
                        f"frames it gridded)" if tail else "")
        else:
            tail = [d for d in days_between(plan["first"], plan["last"]) if sev["first"] <= d < sd]
            grid, why_not = [], {}
            for cam in cams:
                skip = set((read_stamp(out, f"detect_{cam}") or {}).get("skip_days") or [])
                verdict = {r.get("date"): r.get("verdict") for r in
                           ((load_json(out / "pointing" / f"pointing_{cam}.json") or {}).get("rows") or [])}
                g_ = (dg["gridded"].get(cam) or {}).values()
                lo_ = (dg["left_out"].get(cam) or {}).values()
                for d in tail:
                    n = sum(1 for _, dd in g_ if dd == d)
                    if n:
                        grid.append(f"{cam} {d} ({n} frame{'s' if n != 1 else ''})")
                        continue
                    whys = sorted({w_ for _, dd, w_ in lo_ if dd == d})
                    if d in skip and verdict.get(d) and verdict[d] not in LEFT_OUT_VERDICTS:
                        # the detection ran with other pointing results: today's verdict keeps the day
                        k_ = ("left out by an earlier pointing result (the pointing step was run again since; see "
                              "the pointing caveats)")
                        lab = f"{cam} {d} ({verdict[d]})"
                    elif d in skip:
                        k_ = "left out by the pointing check (see the pointing caveats)"
                        lab = f"{cam} {d}" + (f" ({verdict[d]})" if verdict.get(d) else "")
                    elif whys:
                        k_, lab = "every frame left out: " + ", ".join(whys), f"{cam} {d}"
                    else:
                        k_, lab = "no waterline", f"{cam} {d}"
                    why_not.setdefault(k_, []).append(lab)
            nots = "; ".join(f"{', '.join(v)} {k}" for k, v in why_not.items())
            if not tail:
                tail_txt = ""
            elif grid:
                tail_txt = (f"; the DEM grids frames of the storm's tail before the {sd} survey with the rest: "
                            f"{', '.join(grid)}" + (f"; not gridded: {nots}" if nots else ""))
            else:
                tail_txt = (f"; the DEM grids no frame of the storm's tail before the {sd} survey "
                            f"({tail[0]}{' .. ' + tail[-1] if len(tail) > 1 else ''})" + (f": {nots}" if nots else ""))
        text = f"STORM: {sev['name']} ({sev['first']} .. {sev['last']}) {where}{tail_txt}. On this beach {sev['change']}."
        lo, hi = sev.get("loss_northing") or (None, None)
        for s in plan.get("surveys") or []:
            if s["survey_type"] != "points" or lo is None:
                continue
            try:
                from survey_compare import read_points, is_transect
                hits = resolve_survey(s, plan, args)[0] if getattr(args, "survey_dirs", None) is not None else []
                if not hits:                     # not here: nothing to read (never a file of its name elsewhere)
                    continue
                q = read_points(hits[0])
                tr = np.array([is_transect(x) for x in q["desc"]], bool)
                n = q["N"][tr] if tr.sum() >= 3 else q["N"]
                if len(n) and n.min() <= hi and n.max() >= lo:
                    text += (f" The {s['name']} {'transects' if tr.sum() >= 3 else 'points'} span N "
                             f"{n.min():.0f}-{n.max():.0f}: they overlap that loss.")
            except (SystemExit, Exception):
                pass
        cav.append(text)
    if cav and note:
        cav[-1] += " " + PER_DAY_ROWS_NOTE
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


# The same on a BLIND product, which has no comparison: what the holder's own check may see by day.
PER_DAY_ROWS_NOTE_BLIND = (
    "This product has no comparison (blind). The DEM merges the days; a check that splits its error by day needs the "
    "waterline points themselves (waterlines/contour_points_ground_filtered.csv, or contour_points_ground.csv when no "
    "filtered file was used: each point with its frame's capture time, still water and setup; in the --share archive "
    "too; the photo days: provenance.json frames_per_day). Such a check may show a step or trend between days: it may "
    "be change of the beach OR a method error that "
    "depends on the conditions (the elevation, the waves and setup, the tide phase of the hours sampled, which moves "
    "from day to day): on a static synthetic beach this method's per-day medians drift by up to ~0.17 m over a week "
    "(Jan -0.09 m, Mar -0.17 m, with a ~0.07 m step across a wave event), and days at the same median elevation still "
    "differ by up to ~0.13 m. A difference between days below ~0.15-0.2 m is not by itself evidence of change; compare "
    "days at a similar elevation and setup before calling it change.")


def per_day_rows_note(plan):
    """PER_DAY_ROWS_NOTE, or on a BLIND product (no comparison) PER_DAY_ROWS_NOTE_BLIND."""
    return PER_DAY_ROWS_NOTE_BLIND if plan.get("blind") else PER_DAY_ROWS_NOTE


def blind_caveats(plan, cav):
    """The caveats as a BLIND product says them: PER_DAY_ROWS_NOTE (which points at the comparison's per-day
    rows) becomes PER_DAY_ROWS_NOTE_BLIND wherever a caveat carries it. Unchanged otherwise."""
    if not plan.get("blind"):
        return cav
    return [c.replace(PER_DAY_ROWS_NOTE, PER_DAY_ROWS_NOTE_BLIND) if isinstance(c, str) else c for c in cav]


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
# Events closer than this (hours, end of one to start of the next) are one high-wave caveat: the
# beach has no time to settle between them (27-28 Sep 2026: three events in 33 h, the tail of the
# 25-26 Sep storm).
EVENT_MERGE_H = 48.0


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


def wave_event_caveats(plan, prov, first_used, last_used, note=True):
    """One caveat per group of spans of the window with offshore Hs >= EVENT_HS_M in the date's own
    forcing (spans less than EVENT_MERGE_H apart are one group; the hand-kept STORMS table covers only
    the events measured on the station): the frames the DEM gridded before, between, during and after
    them (from its own record, dem_gridded_frames; without it, the waterlines' frames), the frames it
    left out for their waves, where the survey lies, and where to look for change. An event inside a
    STORMS entry is left to that caveat. note: end with PER_DAY_ROWS_NOTE (a README gives it once)."""
    out = plan["out"]
    cav = []
    days = sorted({d for v in (prov.get("frames_per_day") or {}).values() for d in v}) or \
        days_between(first_used, last_used)
    sdays = sorted({x["survey_date"] for x in plan.get("surveys") or [] if x.get("survey_date")})
    sd = sdays[0] if sdays else plan["date"]
    events = [ev for ev in wave_events(out / "forcing", plan["first"], plan["last"])
              if not any(not (datetime.fromtimestamp(ev["end"], tz=timezone.utc).strftime("%Y-%m-%d") < sv["first"]
                              or datetime.fromtimestamp(ev["start"], tz=timezone.utc).strftime("%Y-%m-%d") > sv["last"])
                         for sv in STORMS)]
    groups = []
    for ev in events:
        if groups and ev["start"] - groups[-1][-1]["end"] < EVENT_MERGE_H * 3600:
            groups[-1].append(ev)
        else:
            groups.append([ev])
    dg = dem_gridded_frames(plan) if groups else None
    left = {}                                  # {why: {camera: {frame: epoch}}} the DEM left out
    if dg is not None:
        grid = {c: {k: v[0] for k, v in fr.items()} for c, fr in dg["gridded"].items()}
        for c, fr in dg["left_out"].items():
            for k, v in fr.items():
                left.setdefault(v[2], {}).setdefault(c, {})[k] = v[0]
    else:
        grid = frame_epochs(plan) if groups else {}
    for grp in groups:
        it = "it" if len(grp) == 1 else "them"
        da = datetime.fromtimestamp(grp[0]["start"], tz=timezone.utc).strftime("%Y-%m-%d")
        db = datetime.fromtimestamp(grp[-1]["end"], tz=timezone.utc).strftime("%Y-%m-%d")
        # frames placed by their capture TIME against the events (a wave record stands for the half
        # hour on either side of it), not by calendar day: a one-hour event is not 'during' a whole day
        spans = [(ev["start"] - 1800.0, ev["end"] + 1800.0) for ev in grp]

        def classify(fr):
            cls = {"before": {}, "during": {}, "between": {}, "after": {}}
            for cam, f_ in (fr or {}).items():
                for k, e in f_.items():
                    if e != e:
                        continue
                    c_ = ("before" if e < spans[0][0] else "after" if e > spans[-1][1] else
                          "during" if any(s0 <= e <= s1 for s0, s1 in spans) else "between")
                    dd = datetime.fromtimestamp(e, tz=timezone.utc).strftime("%Y-%m-%d")
                    cls[c_][dd] = cls[c_].get(dd, 0) + 1
            return cls
        cls = classify(grid)
        if dg is None and not any(cls.values()):           # no capture times: by day, as before
            cls = {"before": {d: 0 for d in days if d < da}, "during": {d: 0 for d in days if da <= d <= db},
                   "between": {}, "after": {d: 0 for d in days if d > db}}

        def frames_text(cl, k, what):
            v = cl[k]
            if not v:
                return ""
            n = sum(v.values())
            return (f"{n} frame(s) {what} {it} ({', '.join(sorted(v))})" if n else
                    f"the frames of {', '.join(sorted(v))} {what} {it}")
        parts = [x for x in (frames_text(cls, "before", "before"), frames_text(cls, "between", "between"),
                             frames_text(cls, "during", "during"), frames_text(cls, "after", "after")) if x]
        sides = sum(bool(cls[k]) for k in ("before", "between", "during", "after"))
        # the survey's time of day is not known: on an event day it is 'the same day', not before/after
        rel = (f"after {it}" if sd > db else f"before {it}" if sd < da else
               "on the same day as the event (its time is not known)" if da == db and len(grp) == 1 else
               "within the events' days (its time is not known)" if len(grp) > 1 else
               "within the event's days (its time is not known)")
        evs = []
        for ev in grp:
            a = datetime.fromtimestamp(ev["start"], tz=timezone.utc)
            b = datetime.fromtimestamp(ev["end"], tz=timezone.utc)
            pk = datetime.fromtimestamp(ev["peak_at"], tz=timezone.utc)
            evs.append((f"{a:%Y-%m-%d %H:%M} .. {b:%Y-%m-%d %H:%M} UTC" if ev["end"] > ev["start"] else
                        f"{a:%Y-%m-%d %H:%M} UTC") + f" (peak {ev['peak']:.2f} m at {pk:%Y-%m-%d %H:%M} UTC)")
        srcs = ", ".join(sorted({x for ev in grp for x in ev["sources"] if x})) or "the forcing"
        span = evs[0] if len(evs) == 1 else ", ".join(evs[:-1]) + " and " + evs[-1]
        text = (f"HIGH WAVES in the window: offshore Hs >= {EVENT_HS_M:.1f} m over {span} ({srcs}, forcing/waves.csv"
                + (f"; {len(grp)} spans less than {EVENT_MERGE_H:.0f} h apart, taken together" if len(grp) > 1 else "")
                + f"). The beach may have changed across {it}")
        if dg is None:
            text += (", and the waterlines hold frames of both sides: " if sides > 1 else
                     ": the waterlines hold " if parts else "") + "; ".join(parts)
            text += " (no DEM record says which of them it gridded)" if parts else ""
        elif not parts:
            text += "; the DEM grids no frame of the window"
        else:
            text += ("; the DEM grids the frames of both sides together: " if sides > 1 else
                     "; the DEM grids only ") + "; ".join(parts)
        for why_, fr in sorted(left.items()):
            lcls = classify(fr)
            lparts = [x for x in (frames_text(lcls, "before", "before"), frames_text(lcls, "between", "between"),
                                  frames_text(lcls, "during", "during"), frames_text(lcls, "after", "after")) if x]
            if lparts:
                text += (f"; it left out {sum(sum(v.values()) for v in lcls.values())} frame(s) "
                         f"{'of' if why_.startswith('camera-days') else 'with'} {why_}: " + "; ".join(lparts))
        text += f"; the survey ({sd}) is {rel}."
        cav.append(text)
    if cav and note:
        cav[-1] += " " + PER_DAY_ROWS_NOTE
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
                            "io_file": Path(c["io_path"]).name, "history": hist,
                            "variant": c.get("variant") or "table",
                            "given_with": (f"--eo {cam}={c['eo_given']} --calibration-variant {c.get('variant')}"
                                           if c.get("eo_given") else "survey_dates.csv eo_file"),
                            "table_eo_file": c.get("table_eo_file") or c["eo_file"],
                            "origin": calibration_origin(c["eo_path"] or c["eo_file"], eo_dirs(args))},
            "search_envelope": ({"survey": str(c["envelope_path"] or c["envelope_survey"]),
                                 "sha256": envelope_check(plan, cam, args).get("sha256"),
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
                                "unconfirmed_days": pres.get("unconfirmed_days"),
                                "horizon_overruled_days": pres.get("horizon_overruled_days"),
                                # as the detection ran (its stamp; None: not run, or an older build's stamp)
                                "keep_moved_days": (dst or {}).get("keep_moved_days"),
                                "station_monitor_log": pres.get("station_log"),
                                # the columns of an older build's rows that it did not have are left out
                                "per_day": [{k: r[k] for k in ("date", "verdict", "why", "horizon_offset_px",
                                                                "horizon_dtilt_deg", "horizon_droll_deg",
                                                                "horizon_edge", "features", "inliers",
                                                                "d_azimuth_deg", "d_tilt_deg", "d_roll_deg",
                                                                "turn_sigma", "runner_up_inliers", "runner_up_d")
                                             if k in r}
                                            for r in pres.get("rows", [])]} if pres else None),
            "detection": (dst or {}).get("setup") if plan["era"] == "live" else
            ({"photos": (dst or {}).get("photos"), "skip_days": (dst or {}).get("skip_days")} if dst else None),
            "built": bool(dst) and cam not in failed_cameras(plan, state),
            "notes": c["notes"]}
        if cams[cam]["pointing_check"] is not None and "pointing_correction" in (dst or {}):
            # --correct-pointing (its own root): the days mapped with their own pointing (absent from a default build)
            cams[cam]["pointing_check"]["pointing_correction"] = dst["pointing_correction"]
    # a GCP-free B made by the pre-registered procedure: its stated rotation term's DEM effect taken as a pure TILT,
    # per range band, from this date's DEM and calibration (the procedure's synthetic errors were almost all roll)
    pinned = {cam: v["calibration"]["origin"]["procedure"] for cam, v in cams.items()
              if ((v["calibration"].get("origin") or {}).get("procedure") or {}).get("pinned")}
    for cam, pc_ in pinned.items():
        u_ = pc_.get("uncertainty") or {}
        if u_.get("rotation_deg"):
            why_ = {}
            try:
                bands = tilt_sensitivity(plan, cam, float(u_["rotation_deg"]), 0.0, why=why_)
            except Exception as exc:             # no readable DEM yet: said as not computed
                bands, why_["text"] = None, f"the DEM cannot be read ({exc.__class__.__name__})"
            cams[cam]["calibration"]["rotation_term_as_tilt"] = {
                "tilt_deg": u_["rotation_deg"], "bands": bands, "why_none": None if bands else why_.get("text"),
                "how": "each DEM cell projected into the photo with this calibration and back to the ground at its own "
                       "elevation with the tilt changed by the rotation term; DEM error = local slope x that move; "
                       "median per range band from the camera"}
    fits, fits_from = setup_in_use(plan, args)
    flt = read_stamp(out, "filter")
    mrg = read_stamp(out, "merge") or {}
    dstamp = read_stamp(out, "dem") or {}
    chain_ok, chain_why = chain_current(plan)
    blind = None
    if plan.get("blind"):
        b = plan["blind"]
        per_cam = {cam: dict(zip(("label", "why"), held_check(v["calibration"]["origin"], b,
                                                             envelope_check(plan, cam, args))))
                   for cam, v in cams.items()}
        other = [{"label": "PARTLY-CIRCULAR", "why": f"setup coefficient C = {f['coef']} of unknown origin (no "
                                                     f"--setup-fitted-to): it may have been fitted to the held survey"}
                 for f in fits if f.get("kind") == "unknown"]
        worst = max([v["label"] for v in per_cam.values()] + [o["label"] for o in other], key=LABELS.index,
                    default="PARTLY-CIRCULAR")
        blind = {"statement": blind_statement(b), "held_by": b["held_by"], "survey": b["survey"],
                 "survey_date": b["survey_date"], "since": b.get("since"), "why": b.get("why"),
                 "held_rows": plan.get("held") or [],
                 "how": "survey_products.py BLIND_DATES: no survey of this date, nor of the held survey's date, was "
                        "looked for or opened (resolve_survey refuses before any folder is searched); the compare "
                        "step did not run; the plan view draws no survey",
                 "check": {"label": worst, "per_camera": per_cam, "other": other,
                           "what": "how independent the holder's check of THIS DEM is, by the calibration it was "
                                   "built with and the search envelope that placed its detection (the build never "
                                   "read the held survey; a GCP-free B is built as A was, only the calibration "
                                   "different)"}}
        if pinned:                  # a GCP-free B: what its makers had seen of A, and the order of pin, rule, A and B
            blind["gcpfree"] = gcpfree_disclosure(plan, args, pinned)
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
                      # the waterlines without their 'strata' (survey_compare's report-only tables, which
                      # stay in NAME_comparison.json beside its 'intervals', neither copied here)
                      "headline": ({k: ({kk: vv for kk, vv in head[k].items() if kk != "strata"}
                                        if k == "waterlines" and isinstance(head.get(k), dict) else head.get(k))
                                    for k in ("n", "median", "nmad", "rmse", "mean", "p5", "p95",
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
    info = load_json(Path(f"{dem_stem(plan)}_info.json"))
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
        "cameras_not_delivered": plan.get("not_delivered") or None,
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
                    files=dict({k: f"{dem_stem(plan)}_{k}" for k in ("dem", "spread", "count")},
                               **({k: f"{dem_stem(plan)}_{k}" for k in DEM_LOESS_GRIDS} if dem_set.get("loess") else {})),
                    cells=dstamp.get("cells"),
                    built_from=info,
                    geotiff=(f"asc_to_geotiff.py: tagged EPSG:{dem_set.get('geotiff_epsg')}"
                             + (" (WGS 84 / UTM 19N), NOMINAL" if dem_set.get("geotiff_epsg") == 32619 else "")
                             + "; the grid coordinates are in the GCP frame, assumed NAD83(2011) / UTM zone 19N "
                               "as the surveys (the GCP files state no CRS) "
                               "(EPSG:6348), untransformed; only the _dem GeoTIFF carries vertical keys (NAVD88, "
                               "EPSG:5703); spread, count, sigma and difference grids are not heights"
                             + ("; with --dem-loess the _loess GeoTIFF (a second surface, NAVD88) carries them too; "
                                "its _loess_err, _loess_nmse and _loess_skill are not heights"
                                if dem_set.get("loess") else ""))),
        "chain": {"current": chain_ok, "why": chain_why or None},
        "blind": blind,
        "comparisons": comps,
        "commands": commands,
        "software": software_of_outputs(out),
    }
    prov["frames_per_day"] = frames_per_day(waterline_file(plan, quiet_=True))
    prov["caveats"] = blind_caveats(plan, collect_caveats(plan, args, state, forcing, prov))
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
    from dem_figure import frames_left_out_text
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
    af = wl_.get("above_survey_floor") or {}
    if af.get("frames_all"):
        fl = ", ".join(f"{k.split(' ')[0]} {v:+.2f}" for k, v in (af.get("transect_floors_m") or {}).items())
        L.append(f"      WHERE THE SURVEY STOPS: the transects end at their lowest shot ({fl} m); a line below that "
                 f"can only read LOW there and one that would read high lands seaward of the survey, so the value "
                 f"above is biased LOW by where the survey stops ({af.get('frames_below_floor', 0)} of "
                 f"{af['frames_all']} frames lie below their transect's lowest shot). Over only the points whose line "
                 f"lies at least {af.get('margin_m', 0.1):.2f} m above that shot: "
                 + (f"median {af['median']:+.3f} m on {af['frames']} frames"
                    + (" (too few frames: not an estimate)" if af["frames"] < FEW_N else "")
                    if af.get("frames") else "no frame"))
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
                 + (f"; it left out {fs['left_out_frames']} (" + frames_left_out_text(fs.get("left_out_why")) + ")"
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
        af0 = ws.get("above_survey_floor") or {}
        if af0.get("frames_all"):
            L.append(f"        the C = 0 lines sit lower, so more lie below where the survey stops "
                     f"({af0.get('frames_below_floor', 0)} of {af0['frames_all']} frames) and these C = 0 medians are "
                     f"biased LOW by more; above the transects' lowest shot + {af0.get('margin_m', 0.1):.2f} m: "
                     + (f"median {af0['median']:+.3f} m on {af0['frames']} frames"
                        + (" (too few frames: not an estimate)" if af0["frames"] < FEW_N else "")
                        if af0.get("frames") else "no frame"))
        L.append(f"        ({ws.get('how')})")
    elif ws.get("n"):                    # a comparison made before the pairing (Oct 2026)
        L.append(f"      without the wave setup (C = 0, a sensitivity): median {ws['median']:+.3f} m"
                 + (f"; by camera {per_camera_text(ws.get('by_camera'), unit)}" if ws.get("by_camera") else "")
                 + f" ({ws.get('how')}; NOT on the same frames as the headline)")
    return L


def level_words(prov):
    """The README's words for the water level each line is put at: 'the measured water level' was
    printed for the October 2024 date too, whose every level is the Chatham gauge transferred (and
    extrapolated 44-51 days beyond the ADCP record), none measured at Marconi."""
    srcs = (prov.get("water_level") or {}).get("sources_used") or {}
    n_tr = sum(int(v) for k, v in srcs.items() if "chatham" in str(k))
    n_all = sum(int(v) for v in srcs.values())
    if n_all and n_tr == n_all:
        return ("the Chatham gauge's water level transferred to Marconi (no level of this window was measured "
                "at Marconi; see KNOWN CAVEATS) ")
    if n_tr:
        return (f"the measured water level ({n_tr} of {n_all} levels are the Chatham gauge transferred to "
                "Marconi; see KNOWN CAVEATS) ")
    return "the measured water level "


def dem_built_from_lines(info):
    """The README's lines (WHAT WAS USED, under the DEM) on what the DEM gridded and what its filters
    left out, from its own record (<stem>_info.json): the frames and days, then the day test's record --
    each camera-day it left out, how and why (only a day 'left out' is off the rest by more than the
    limit: dem_figure's wording, as on the DEM page) and its WARNINGs (also among the KNOWN CAVEATS)."""
    from dem_figure import day_test_left_out, day_test_rejected_lines, day_test_warnings
    flt_ = info.get("filters") or {}
    out = [f"      gridded {sum(info['frames'].values())} frames ("
           + ", ".join(f"{k} {v}" for k, v in sorted(info["frames"].items())) + f") on {info.get('days')} days; "
           f"left out before gridding: {flt_.get('rough_frames', 0)} frames with Hs > {flt_.get('max_hs')} m, "
           + (day_test_left_out(flt_, "+/-") if flt_.get("max_day_offset") else "no day test")]
    if flt_.get("max_day_offset"):
        out += [f"        day test: {ln}" for ln in day_test_rejected_lines(flt_)]
        out += [f"        day test {ln}" for ln in day_test_warnings(flt_)]
    return out


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
    bl = prov.get("blind")
    if bl:
        w(f"BLIND PRODUCT, DEM ONLY: {bl['statement']}.")
    cals = {cam: (v.get("calibration") or {}) for cam, v in (prov.get("cameras") or {}).items()}
    var = sorted({c.get("variant") or "table" for c in cals.values()})
    if bl or var != ["table"]:
        w("CALIBRATION: " + (f"variant {'/'.join(v for v in var if v != 'table')} (--eo), NOT the table's"
                             if var != ["table"] else "the table's (surveys/survey_dates.csv eo_file)")
          + ": " + "; ".join(f"{cam} {c.get('eo_file')}" for cam, c in sorted(cals.items()))
          + " (what each was solved from: WHAT WAS USED)")
    for cam, c in sorted(cals.items()):
        pc = (c.get("origin") or {}).get("procedure") or {}
        if pc.get("format"):
            w(f"  {cam} {c.get('eo_file')}: {procedure_text(pc)}")
            if pc.get("pinned") and uncertainty_terms(pc):
                w(f"    STATED UNCERTAINTY of this calibration (pre-registered rule): {uncertainty_terms(pc)}")
                w(f"    {tilt_effect_text(c.get('rotation_term_as_tilt'))}")
    for cam, v in sorted((prov.get("cameras_not_delivered") or {}).items()):
        w(f"NOT DELIVERED: {cam}: the pre-registered GCP-free procedure did not deliver it ({v.get('verdict')}; "
          f"{'; '.join(v.get('reasons') or []) or 'no reason given'}). Its report: {v.get('report')} (sha256 "
          f"{str(v.get('report_sha256'))[:12]}..., run {v.get('run_id')}). This DEM has no {cam} data: the table's "
          f"calibration is never used in its place.")
    w("")
    w("WHAT THIS IS")
    fpd = prov.get("frames_per_day") or {}
    used = sorted(c for c, v in fpd.items() if v)
    days = sorted({dd for v in fpd.values() for dd in v})
    w(f"An intertidal beach DEM of Marconi Beach made only from the station's camera photos "
      f"({', '.join(used) or 'no camera'}"
      + (f"; {', '.join(c for c in plan['cams'] if c not in used)} contributed nothing" if
         [c for c in plan["cams"] if c not in used] else "")
      + (f"), for an independent check by {bl['held_by']}, who holds the survey of {d} ({bl['survey']}). "
         if bl else f"), to check against the survey(s) of {d} listed below. ")
      + f"Photo window {plan['first']} .. {plan['last']}; waterline frames on "
      + (f"{days[0]} .. {days[-1]} ({len(days)} days)" if days else "no day") + dem_days_text(plan, fpd) + ".")
    fits_ = (prov.get("setup") or {}).get("fits") or []
    no_setup = bool(fits_) and all(f.get("kind") == "none" or not f.get("coef") for f in fits_)
    w("Each photo's waterline is where the water met the sand; its elevation is " + level_words(prov)
      + ("only: NO wave setup (C = 0), so a line marked by the swash reads LOW, by roughly "
         f"{SHARE_KEPT_RANGE[0]}-{SHARE_KEPT_RANGE[1]} x the setup (the line also moves seaward onto lower beach). "
         if no_setup else
         f"plus the wave setup (C = {(prov.get('setup') or {}).get('coef')}). ")
      + "Many waterlines at different tides make the DEM.")
    w("")
    if bl:
        w("BLIND: NO COMPARISON (the survey is held for an independent check)")
        w(f"  {bl['statement']}.")
        w(f"  The survey: {bl['survey']}, held by {bl['held_by']} since {bl.get('since')}. This build never looked "
          f"for a survey file of {d}" + (f" or of {bl['survey_date']}" if bl["survey_date"] != d else "")
          + ", under any name or in any folder (survey_products.py "
          "BLIND_DATES), made no comparison and drew no survey on the maps: no survey of that date, whatever "
          "lies on this computer, went into this product.")
        ck = bl.get("check") or {}
        w(f"  How independent his check of this DEM is: {ck.get('label')} (by the calibration it was built with and "
          f"the search envelope that placed its detection, each named below with its sha256; a GCP-free B is also "
          f"built exactly as A was, only the calibration different: the committed tables, A's options, envelope "
          f"and forcing)")
        for cam, v in sorted((prov.get("cameras") or {}).items()):
            se = v.get("search_envelope") or {}
            if se.get("survey"):
                w(f"    {cam} search envelope: {se.get('survey')}" + (f" sha256 {se['sha256']}" if se.get("sha256")
                                                                       else ""))
        for cam, v in sorted((ck.get("per_camera") or {}).items()):
            w(f"    {cam}: {v['label']}: {v['why']}")
        for o in ck.get("other") or []:
            w(f"    {o['label']}: {o['why']}")
        gf = bl.get("gcpfree") or {}
        if gf:
            w(f"  {gf.get('exposure')}")
            w("  ORDER (what decides B, and A): " + "; ".join(gf.get("order") or []) + ".")
            for x in gf.get("order_flags") or []:
                w(f"    NOTE: {x}.")
            for x in gf.get("protocol") or []:
                w(f"  {x}.")
        if ck.get("label") != "INDEPENDENT":
            w("  A DEM of this date built with a calibration that used none of his points (--eo ... "
              "--calibration-variant, in its own --output-root) lets his check test the camera geometry too.")
    else:
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
        org = cal.get("origin") or {}
        w(f"      calibration {cal['eo_file']}"
          + (f" (--eo, variant {cal.get('variant')}; the table names {cal.get('table_eo_file')})"
             if str(cal.get("given_with", "")).startswith("--eo") else "")
          + (f": {org['text']}" if org.get("kind") in ("gcp", "rotated") or cal.get("variant") not in (None, "table")
             or bl else f" ({hist})" if hist else ""))
        if org.get("kind") == "gcp" and hist:
            w(f"        ({hist})")
        if org.get("acceptance"):
            w(f"        pre-registered acceptance: {org['acceptance']['verdict']}"
              + (f" ({org['acceptance']['text']})" if org["acceptance"].get("text") else ""))
        w(f"      search envelope: {c['search_envelope']['survey'] or c['search_envelope']['how']}")
        pc_ = c.get("pointing_check") or {}
        if pc_.get("per_day"):
            verdicts = {}
            for r in pc_["per_day"]:
                verdicts.setdefault(r["verdict"], []).append(r["date"])
            w("      pointing check: " + "; ".join(f"{k} {', '.join(v)}" for k, v in verdicts.items()))
        pcr = pc_.get("pointing_correction")
        if pcr is not None:
            w("      pointing CORRECTED per day (--correct-pointing): "
              + ("; ".join(f"{d} {v['d_deg'][0]:+.2f}/{v['d_deg'][1]:+.2f}/{v['d_deg'][2]:+.2f} deg pan/tilt/roll"
                           for d, v in sorted((pcr.get("days") or {}).items())) or "no day passed the gates")
              + " (KNOWN CAVEATS)")
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
          + (f"; transfer CV RMS {wl.get('cv_rms_m')} m, extrapolation RMS {wl.get('extrapolation_rms_m')} m pooled "
             f"over both directions"
             + (f" (each level's sigma_m {wl.get('transfer_sigma_m')} m: {wl.get('transfer_sigma_basis')})"
                if wl.get("transfer_sigma_m") is not None and wl.get("transfer_sigma_basis") else
                "; KNOWN CAVEATS gives the test of this date's side")
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
        for ln in dem_built_from_lines(bf):
            w(ln)
    if dm.get("loess"):
        lo = dm["loess"]
        w(f"      loess (--dem-loess): {lo.get('script')} half-widths {lo['half_widths_m']['alongshore']:g} m "
          f"alongshore x {lo['half_widths_m']['cross_shore']:g} m cross-shore: {lo.get('note')}")
    w(f"      GeoTIFFs: {dm['geotiff']}")
    sw = prov["software"]
    w(f"  software: {sw.get('note')} (per step in provenance.json software.by_step)")
    w("")
    w("FILES")
    w("  forcing/      water_level.csv, waves.csv, forcing_report.txt, forcing.png, forcing.json")
    w("  pointing/     per-day pointing check per camera")
    w("  waterlines/   <cam>/ (detection work folders), contour_points_ground.csv (all), "
      "contour_points_ground_filtered.csv (after the consistency filter), consistency_report.csv and "
      "consistency_<cam>.png (what the filter dropped, this build's cameras only), "
      "no_setup_frames.csv (frames left out: no wave record)")
    ds = dem_stem(plan).name
    w(f"  dem/          {ds}_dem/_spread/_count/_sigma .asc and .tif, {ds}_dem.png (the DEM page)"
      + (f"; {ds}_loess/_loess_err/_loess_nmse/_loess_skill .asc and .tif (--dem-loess: the loess surface and its "
         f"prediction error, dem_loess.py)" if (prov.get("dem") or {}).get("loess") else "")
      + (f" (the calibration variant in the file names: {plan_variant(plan)})" if plan_variant(plan) else ""))
    w(f"  maps/         {map_png(plan, '<cam>').name}: waterlines drawn on the photos, per camera "
      f"(daily_elevation_map.py, as the station's own maps: its colour bar says which elevation colours the lines); "
      f"{map_png(plan, 'plan').name}: plan view "
      + ("(BLIND: no survey drawn; " if bl else "over the survey (")
      + "coloured by still water + setup, the elevation the DEM uses; its title names the calibration)"
      + (f" (the calibration variant in the file names: {plan_variant(plan)})" if plan_variant(plan) else ""))
    w("  compare/      " + (f"none: a blind product ({bl['statement']})" if bl else
                           "DEM - survey grids/points, report, figure, json"))
    w("  provenance.json  every source above, machine-readable")
    w("  (survey_products.py --share packs all of these into one archive with a MANIFEST, but not waterlines/<cam>/, "
      "pointing/reference_<cam>/ or logs/stamps/, and of the two waterline files only the one the DEM used)")
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

def forcing_input_notes(era, found, setup_coef):
    """What the dry run says of each forcing input historical_forcing.py will use ({key: Path or None},
    forcing_inputs()): -> ([(key, text)], [warnings]). The chatham era needs the ADCP too: the transfer
    is fitted on its overlap with Chatham and the WIS / NDBC waves are converted to its Hs, so without
    it the forcing step fails; it was listed as 'NOT FOUND' with nothing to say it mattered."""
    need = ({"adcp": "REQUIRED: the Chatham transfer is fitted on it and the WIS/NDBC waves are converted to "
                     "its Hs",
             "adcp_navd88": "the ADCP on NAVD88 (else its datum is derived from the Chatham mean, as "
                            "adcp_to_navd88.py does)",
             "chatham": "the Chatham gauge (else downloaded from NOAA when the period needs it)",
             "wis": "WIS ST63064 waves (no download: required for the setup unless NDBC 44013 covers)",
             "ndbc": "NDBC 44013 waves"} if era == "chatham" else
            {"adcp": "the ADCP (water level and waves)", "adcp_navd88": "the ADCP on NAVD88"})
    lines = [(k, (str(p) if p else "NOT FOUND") + (f"  <- {need[k]}" if k in need and not p else ""))
             for k, p in found.items()]
    warns = []
    if era == "chatham" and found and not found.get("adcp"):
        warns.append("no ADCP record (sig1000_waves_ALL.csv): the forcing step will FAIL (no Chatham transfer can "
                     "be fitted and the WIS/NDBC waves cannot be converted); give --adcp")
    if era == "chatham" and setup_coef and not (found.get("wis") or found.get("ndbc")):
        warns.append("no WIS or NDBC wave file: the forcing step will FAIL (no waves for the setup); give --wis / "
                     "--ndbc or put the files in the Chelsea_calibration folder")
    return lines, warns


def dry_run(plan, args, steps):
    rule(f"{plan['date']} (DRY RUN: nothing is written)")
    out = plan["out"]
    say("era", plan["era"])
    say("output", str(out))
    say("window", f"{plan['first']} .. {plan['last']}")
    for cam, why in (plan.get("disabled_cams") or {}).items():
        say(f"{cam}", f"not built (disabled): {why}")
    for cam, v in (plan.get("not_delivered") or {}).items():
        say(f"{cam}", f"NOT in this build: the pre-registered GCP-free procedure did not deliver it ({v.get('verdict')}; "
                      f"report {v.get('report')})")
    total = 0
    for cam, c in plan["cams"].items():
        say(f"{cam}", f"station {c['station']}, {c['first']} .. {c['last']} {c['hours_text']} UTC")
        say("  calibration", f"{c['eo_file']} -> {c['eo_path'] or 'NOT FOUND in ' + str(args.calibration)}"
            + (f" (--eo, variant {c.get('variant')}; the table names {c.get('table_eo_file')})" if c.get("eo_given") else ""))
        org = calibration_origin(c["eo_path"] or c["eo_file"], eo_dirs(args))
        say("  solved from", org["text"])
        if (org.get("procedure") or {}).get("format"):
            say("  procedure", procedure_text(org["procedure"]))
        if org.get("acceptance"):
            say("  acceptance", org["acceptance"]["verdict"]
                + (f" ({org['acceptance']['text']})" if org["acceptance"].get("text") else ""))
        if plan.get("blind"):
            say("  held check", " -- ".join(held_check(org, plan["blind"], envelope_check(plan, cam, args))))
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
        lines, warns = forcing_input_notes(plan["era"], found, args.setup_coef)
        for k, t in lines:
            say(f"  {k}", t)
        for t in warns:
            warn(t)
        from historical_forcing import setup_share_inputs
        sh = setup_share_inputs(getattr(args, "gnssr_spline", None), getattr(args, "gauge_csv", None))
        say("  setup share", "2026 record for the share of setup GNSS-R sees: "
            + ", ".join(f"{Path(p).name} {'here' if p and Path(p).exists() else 'NOT here'}" for p in sh))
    if plan.get("blind"):
        b = plan["blind"]
        say("BLIND", f"{blind_statement(b)} ({b['survey']}): DEM only; the compare step will not run and no survey "
                     f"of {plan['date']} is looked for")
    for s in plan["surveys"]:
        paths, note = resolve_survey(s, plan, args)
        v = Verdict(s)
        # a survey not found is labelled by its name in the survey folders (survey_path_for), never by a
        # file of that name in the folder the command is run from: the label the build would give
        for p_ in (paths or [survey_path_for(s, plan, args)]):
            honest_label(s, plan, args, p_, fits=[setup_fit_info(args.setup_coef, args)], verdict=v)
        label, reasons = v.label, v.reasons
        say(f"survey {s['name']}", f"{paths[0] if paths else 'MISSING'}{' (' + note + ')' if note else ''}")
        if not paths:
            say("  then", missing_survey_message(s, plan, args, note))
        say("  label", label + ("" if label == s["label"] else f" (table: {s['label']}; " + "; ".join(reasons) + ")")
            + (f" [overruled in surveys.csv ({v.override}): {'; '.join(v.overridden)}]" if v.overridden else ""))
    if plan["era"] == "live":
        # no detection: the rows are read from the station's contour file (~1 min per GB), then the
        # pointing check, georectification, filter, DEM, maps and comparison (79 s for 67 MB of rows
        # off the station, 2 cameras, 1 survey; the NUC's 2 cores take longer)
        p = Path(args.live_contours)
        gb = p.stat().st_size / 1e9 if p.exists() else 1.0
        est = 20 + max(10.0, 60.0 * gb) + 40 * len(plan["cams"]) + 50 + 15 * len(plan["surveys"])
        how = ("no detection: reads the live contour file, ~1 min per GB, then pointing, filter, DEM, maps and "
               "comparison")
    else:
        est = 30 + total * SEC_PER_PHOTO_DETECT + 60 + 40 * len(plan["cams"]) + 40 * len(plan["surveys"])
        # the NUC10i3's rate has never been measured (no date with detection was built there): the
        # estimate is the sandbox's rate, so a range up to twice it is given, and where to read the real one
        how = (f"{total} photos at ~{SEC_PER_PHOTO_DETECT:.0f} s each as measured OFF the station; the NUC10i3's rate "
               f"is NOT measured yet: allow up to ~{max(1.0, (est + total * SEC_PER_PHOTO_DETECT) / 60):.0f} min "
               f"(twice the detection) and read the real rate from the first 'detecting : photo k of N ... "
               f"(S s per photo)' line")
    if getattr(args, "correct_pointing", False):
        say("pointing", "CORRECTED per day (--correct-pointing): each day whose matched features pass "
            "pointing_daily.py's gates is detected and mapped with its own pointing; the others as by default")
    say("steps", ", ".join(steps))
    for st in steps:
        key = {"forcing": "forcing", "filter": "filter", "dem": "dem"}.get(st)
        if key:
            d = read_stamp(out, key)
            say(f"  {st}", "built before (" + d["finished_utc"] + ")" if d else "not built yet")
    rep = setup_replace_note(plan, args, steps)
    if rep:
        warn(rep + (" --replace given: a real run replaces it." if args.replace else
                    " A real run would be REFUSED (exit 5) without --replace."))
    rep_w = window_replace_note(plan, args, steps)
    if rep_w:
        warn(rep_w)
    rep_c = calibration_replace_note(plan, args, steps)
    if rep_c:
        warn(rep_c + (" --replace given: a real run replaces it." if args.replace else
                      " A real run would be REFUSED (exit 5) without --replace."))
    bad = option_refusal(plan, args)
    if bad:
        warn(f"a real run would be REFUSED (exit 5): {bad}")
    else:
        bb = blind_b_build_problem(plan, args)
        if bb:
            warn(f"a real run would be REFUSED (exit 5): NOT the pre-registered B build: {bb}")
        elif gcpfree_b_of(plan, args):
            rec = b_build_record(gcpfree_b_of(plan, args))
            say("B build", ("the pre-registered B build (the committed tables, A's options, envelope and forcing "
                            "inputs); " + (f"its root and command as recorded in {rec['path']}" if rec.get("exists") else
                                           f"a real run records this root and command in {rec.get('path')} (commit "
                                           f"it then)")))
    for w_ in blind_independence_warnings(plan, args):
        warn(w_)
    if plan.get("blind") and blind_leftovers(plan):
        warn(f"a real run would be REFUSED (exit 5): {plan['out']} holds outputs made with a survey of this date "
             f"({'; '.join(blind_leftovers(plan))}): build the blind product in a fresh --output-root")
    if args.setup_coef and not args.setup_fitted_to and not setup_fit_known(args.setup_coef):
        warn(f"a real run would be REFUSED (exit 5): --setup-coef {args.setup_coef} is not a known fit; give "
             f"--setup-fitted-to <survey file> or --setup-fitted-to 'none:<how>'")
    say("expected", f"~{max(1.0, est / 60):.0f} min ({how})")
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
    if plan.get("blind"):                # --summary lists it as blind (DEM only) while it builds
        b = plan["blind"]
        prov["blind"] = dict(prov.get("blind") or {}, statement=blind_statement(b), held_by=b["held_by"],
                             survey=b["survey"], survey_date=b["survey_date"], held_rows=plan.get("held") or [])
        prov["comparisons"] = []
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
    another program after a power cut is not one (Linux: /proc/<pid>/cmdline says): a build is a
    process one of whose command-line arguments (NUL-separated there) is a path whose last part is
    survey_products.py. 'survey_products' anywhere else in the command line (an interpreter in a
    folder of that name, a file such as survey_products_notes.txt, a test folder) is not one."""
    try:
        os.kill(int(pid), 0)
    except PermissionError:
        return True
    except (OSError, ValueError, TypeError):
        return False
    try:
        argv = Path(f"/proc/{int(pid)}/cmdline").read_bytes().split(b"\0")
    except OSError:
        return True                        # no /proc: alive is all that can be told
    return any(os.path.basename(a) == b"survey_products.py" for a in argv)


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


def photo_scan_refusal(plan, args, steps):
    """The photo folders scanned before any step runs (scan_photos is cached: the steps reuse it): a photo-named file
    there that is, or links to, a file named like a held survey is refused by this process's guard (BlindSurvey, before
    it is opened) while the scan reads its size. Said as a refusal, before anything is written (no forcing step, no
    folder). -> text or ''."""
    if plan.get("era") == "live" or not (args.dry_run or {"pointing", "detect"} & set(steps)):
        return ""
    try:
        for cam in plan["cams"]:
            scan_photos(plan, cam, args)
    except BlindSurvey as exc:
        return f"a file in the photo folders is (or links to) a file named like the held survey: {exc}"
    return ""


def build_date(date, cfg, surveys, args, steps):
    args = argparse.Namespace(**vars(args))          # this date's own copy (built settings go in it)
    plan = build_plan(date, cfg, surveys, args)
    if plan is None:
        return "disabled"
    state = {"steps": {}, "failures": [], "not_applied": [], "requested": list(steps)}
    keep_built_window(plan, args, steps, state)
    bad = option_refusal(plan, args)
    if not bad:
        keep_built_calibration(plan, args, steps, state)
        bad = blind_b_build_problem(plan, args)
        if bad:
            bad = f"NOT the pre-registered B build: {bad}"
    ps = photo_scan_refusal(plan, args, steps)
    if ps:                                           # before anything is written, a dry run too
        rule(f"{date}  era {plan['era']}  cameras {', '.join(plan['cams'])}")
        warn(f"{date}: REFUSED (exit 5): {ps}. NOTHING built")
        return "refused"
    if args.dry_run:
        dry_run(plan, args, steps)
        return "dry-run"
    rule(f"{date}  era {plan['era']}  cameras {', '.join(plan['cams'])}  photos {plan['first']} .. {plan['last']}")
    say("output", str(plan["out"]))
    if bad:                                          # before anything is written: no folder, no lock
        warn(f"{date}: NOT built: {bad}")
        return "refused"
    for w_ in blind_independence_warnings(plan, args):
        warn(w_)
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
        if plan.get("blind"):
            left = blind_leftovers(plan)
            if left:
                warn(f"{date}: NOT built: {plan['out']} holds outputs an earlier build made with a survey of this "
                     f"date ({'; '.join(left)}). A blind product ({blind_statement(plan['blind'])}) is never built "
                     f"over them: build it in a fresh --output-root, or move them away yourself")
                return "refused"
        rep_c = calibration_replace_note(plan, args, steps)
        if rep_c and not args.replace:
            warn(f"{date}: NOT built: {rep_c}")
            return "refused"
        if rep_c:
            warn(rep_c + " (--replace given)")
            state["replaced_calibration"] = rep_c
        rep = setup_replace_note(plan, args, steps)
        if rep and not args.replace:
            warn(f"{date}: NOT built: {rep}")
            return "refused"
        if rep:
            warn(rep + " (--replace given)")
            state["replaced_c"] = replaced_coefs(plan, args, steps)
        rep_w = window_replace_note(plan, args, steps)
        if rep_w:
            warn(rep_w)
            state["replaced_window"] = rep_w
        rec = write_b_build_record(plan, args, gcpfree_b_of(plan, args))
        if rec:
            warn(f"{date}: this B's root and command are recorded in {rec} (its first build): COMMIT IT NOW (git -C "
                 f"{rec.parent} add {rec.name} && git -C {rec.parent} commit -m 'gcpfree: B build'); the B is built "
                 f"and shared from {os.path.realpath(str(args.output_root))} only, with this command")
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
                if st == "forcing":
                    fb = b_forcing_problem(plan, args)
                    if fb:
                        raise StepFailed(f"NOT the pre-registered B build: {fb}")
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
            warn("README.txt and provenance.json say 'status: failed' (interrupted); rerun the same command")
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
    if prov.get("blind"):
        bl = prov["blind"]
        say("BLIND", f"DEM only: {bl['statement']}; his check of this DEM: {(bl.get('check') or {}).get('label')}")
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
    if state.get("replaced_c"):
        gone = state["replaced_c"]
        now = round(float(args.setup_coef or 0.0), 6)
        gone_t = " / ".join(f"{c:g}" for c in gone)
        keep = (f"a C = {gone_t} build goes in its own --output-root, {root_for_c(args.output_root, gone[0], [now])}"
                if len(gone) == 1 else "another C goes in its own --output-root")
        misnamed = (root_name_note(args.output_root, gone[0], [now], how="this run's --replace put it there")
                    if len(gone) == 1 else "")
        warn(f"this run REPLACED the C = {gone_t} product in {plan['out']} with a C = {now:g} build (--replace; "
             f"{keep})" + (f". NOTE: {misnamed}" if misnamed else ""))
    if state.get("replaced_window"):
        warn(f"this run REPLACED the product built for another window/hours in {plan['out']} (a recheck with "
             f"another window goes in its own --output-root)")
    if state.get("replaced_calibration"):
        warn(f"this run REPLACED the product built with another calibration in {plan['out']} (--replace; another "
             f"calibration goes in its own --output-root)")
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
        b.dem_loess = bool(dem.get("loess"))     # --dem-loess: whether the DEM on disk was built with the loess
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
    # the calibration the pointing and waterlines on disk were built with (their stamps)
    bc = built_calibration(plan)
    if bc and len(bc) == len(plan["cams"]):
        given = {cam: v.get("eo_path") for cam, v in bc.items()
                 if v.get("eo_path") and str(v["eo_path"]) != str(plan["cams"][cam].get("table_eo_path"))}
        b.eo = given or None
        b.calibration_variant = next((v.get("variant") for v in bc.values() if v.get("variant")), None) if given else None
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
    if getattr(args, "eo", None):
        a += ["--eo"] + [f"{c}={p}" for c, p in sorted(args.eo.items())]
    if getattr(args, "calibration_variant", None):
        a += ["--calibration-variant", str(args.calibration_variant)]
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
    if getattr(args, "dem_loess", False):
        a.append("--dem-loess")
    if args.keep_moved_days:
        a.append("--keep-moved-days")
    if getattr(args, "correct_pointing", False):
        a.append("--correct-pointing")
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
    cur, why = comparison_current(product_plan(prov, pdir), {"name": c["name"]})
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
    blind_rows = []
    for pj in sorted(root.glob("*/provenance.json")):
        prov = load_json(pj) or {}
        if prov.get("blind"):
            bl = prov["blind"]
            ck = bl.get("check") or {}
            cal = sorted({(v.get("calibration") or {}).get("variant") or "table"
                          for v in (prov.get("cameras") or {}).values()})
            blind_rows.append({"date": prov.get("date"), "era": prov.get("era"), "status": prov.get("status"),
                               "survey": ", ".join(h.get("name") or "" for h in bl.get("held_rows") or []) or "-",
                               "survey_type": ", ".join(sorted({h.get("survey_type") or "" for h in bl.get("held_rows") or []})),
                               "survey_date": bl.get("survey_date"), "label": "blind (DEM only)",
                               "comparison": f"blind (DEM only): {bl.get('statement')}",
                               "calibration_variant": "/".join(cal),
                               "held_check_label": ck.get("label"),
                               "held_check_why": "; ".join(f"{k}: {v.get('why')}"
                                                           for k, v in sorted((ck.get("per_camera") or {}).items())),
                               "synthetic": prov.get("synthetic") or None, "built_utc": prov.get("built_utc")})
            continue
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
            # no label without a current comparison: a label is a statement about numbers shown. A
            # survey whose file has not arrived (the Oct 2024 GCPs) was never compared: not 'stale'
            nofile = not cur and c.get("status") == "skipped" and not c.get("survey_path")
            if nofile:
                why_not = "not compared: no survey file yet (" + (c.get("reason") or "").split(". Provide")[0] + ")"
            label = (lab.get("label") or c.get("label") or "") if cur else ("no survey file" if nofile
                                                                           else "not current")
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
                   "comparison": (c.get("status") if cur else why_not if nofile else f"not current: {why_not}"),
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
                   # transect checks: the same lines over only those above where the survey stops
                   "waterlines_above_floor_median_m": num((wl.get("above_survey_floor") or {}).get("median"))
                   if (wl.get("above_survey_floor") or {}).get("frames") else None,
                   "waterlines_above_floor_frames": (wl.get("above_survey_floor") or {}).get("frames"),
                   "waterlines_below_floor_frames": (wl.get("above_survey_floor") or {}).get("frames_below_floor"),
                   "waterlines_floor_frames_all": (wl.get("above_survey_floor") or {}).get("frames_all"),
                   "c0_below_floor_frames": (ws.get("above_survey_floor") or {}).get("frames_below_floor"),
                   "c0_floor_frames_all": (ws.get("above_survey_floor") or {}).get("frames_all"),
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
        built_dates = {r["date"] for r in rows} | {r["date"] for r in blind_rows}
        for d in sorted({r["date"] for r in cfg}):
            rr = [r for r in cfg if r["date"] == d]
            if d in built_dates:
                continue
            if blind_of(d):
                # held for an independent check: listed, never looked for (no label_when_built)
                on = any(r["enabled"] in ("1", "yes", "true", "True") for r in rr)
                st_ = (load_json(root / d / "provenance.json") or {}).get("status")
                blind_rows.append({"date": d, "survey": ", ".join(x["name"] for x in srv if x["date"] == d) or "-",
                                   "label": "blind (DEM only)", "status": st_,
                                   "comparison": "blind (DEM only): " + (
                                       f"product status {st_}" if st_ else
                                       f"not built yet (no {d}/ under this root)" if on else
                                       f"disabled in {Path(args.config).name}")})
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
    if not rows and not disabled and not not_yet and not blind_rows:
        print(f"no product with a comparison under {root} yet")
        return 1
    fields = [k for k in (rows[0].keys() if rows else ["date", "era", "status", "survey", "survey_type", "survey_date",
                                                       "label", "comparison"])
              if k not in ("dem_by_camera", "waterlines_by_camera")]
    for k in ("table_label", "label_when_built", "label_when_built_reason", "calibration_variant",
              "held_check_label", "held_check_why", "synthetic", "built_utc"):
        if k not in fields:
            fields.append(k)
    with open(root / "summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
        w.writerows(blind_rows)
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
        if r.get("waterlines_floor_frames_all"):
            af_n = r.get("waterlines_above_floor_frames") or 0
            lines.append(f"{'':<11}   where the survey stops: {r.get('waterlines_below_floor_frames') or 0} of "
                         f"{r['waterlines_floor_frames_all']} frames lie below their transect's lowest shot, where a "
                         f"line can only read low (the lines median is biased LOW); the lines above it: "
                         + (f"{r['waterlines_above_floor_median_m']:+.3f} m on {af_n} frames"
                            + (" (too few: not an estimate)" if af_n < FEW_N else "") if af_n else "none"))
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
                         + (f"; at C = 0 {r['c0_below_floor_frames']} of {r['c0_floor_frames_all']} frames lie below "
                            f"where the survey stops: biased LOW" if r.get("c0_below_floor_frames") else "")
                         + (("; by camera " + ", ".join(
                             f"{k} {r[f'c0_without_setup_median_{k}_m']:+.3f} / {r[f'c0_with_setup_median_{k}_m']:+.3f}"
                             + (" (too few)" if r[f"c0_n_independent_{k}"] < FEW_N else "") for k in cams0))
                            if len(cams0) > 1 else ""))
        if r["survey_compare_suggests"] and r["survey_compare_suggests"] != r["label"]:
            lines.append(f"{'':<11} survey_compare.py's checks suggest {r['survey_compare_suggests']}"
                         + (f" (overruled: {r['overruled']})" if r["overruled"] else ""))
    if blind_rows:
        lines += ["", "Blind (DEM only): the survey is held by someone else for an independent check of the DEM; "
                      "nothing here looked for it or compared it:"]
        for r in blind_rows:
            lines.append(f"  {r['date']} {r['survey']}: {r['comparison']}"
                         + (f"; status {r['status']}" if r.get("status") else "")
                         + (f"; calibration {r['calibration_variant']}" if r.get("calibration_variant") else "")
                         + (f"; the holder's check of this DEM is {r['held_check_label']} ({r['held_check_why']})"
                            if r.get("held_check_label") else ""))
            if r.get("synthetic"):
                lines.append(f"      SYNTHETIC FIXTURE: tests the pipeline, not the beach ({r['synthetic']})")
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
    if not any(r["median_m"] is not None for r in rows):
        # nothing to plot (a root of blind products, DEMs only; a date whose survey file has not arrived):
        # an older figure would show numbers that are no longer the table's
        if (root / "summary.png").exists():
            (root / "summary.png").unlink()
        print(f"\nwrote {root / 'summary.csv'}, {root / 'summary.txt'} (no summary figure: no comparison has "
              f"numbers yet under this root" + ("; a blind product is a DEM only" if blind_rows else "") + ")")
        return 0
    try:
        summary_figure(rows, root / "summary.png", disabled, not_yet, blind_rows)
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
    if blind_of(date) or held_survey(sv):
        return "blind (DEM only)", "the survey is held for an independent check: never looked for here"
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            plan = build_plan(date, cfg, srv, args)
            if not plan:
                return sv["label"], ""
            paths, _ = resolve_survey(sv, plan, args)
            v = Verdict(sv)
            for p_ in (paths or [survey_path_for(sv, plan, args)]):     # never a file of the current folder
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


def summary_figure(rows, path, disabled=(), not_yet=(), blind=()):
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
        if r.get("waterlines_below_floor_frames"):
            extra.append(f"lines: {r['waterlines_below_floor_frames']} of {r['waterlines_floor_frames_all']} frames below "
                         f"where the survey stops (read low by construction)")
        if c0_ok[k]:
            extra.append(f"C = 0, the same {r.get('c0_paired_n')} {r.get('c0_paired_unit')}: "
                         f"{r['c0_without_setup_median_m']:+.2f} m, {r['c0_label']}"
                         + (" (too few)" if r.get("c0_too_few") else "")
                         + (f", biased low ({r['c0_below_floor_frames']} of {r['c0_floor_frames_all']} below the "
                            f"survey's floor)" if r.get("c0_below_floor_frames") else ""))
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
    if disabled or not_yet or blind:
        dd = sorted({(r["date"], r["comparison"].split("(", 1)[-1].rstrip(")")) for r in disabled})
        txt = ("Not built, at the user's request: " + "; ".join(
            f"{d} ({textwrap.shorten(w, 60, placeholder='...')})" for d, w in dd)) if dd else ""
        ny = sorted({r["date"] for r in not_yet})
        if ny:
            txt = (txt + ". " if txt else "") + "Enabled, no comparison under this root yet: " + ", ".join(ny)
        bd = sorted({r["date"] for r in blind})
        if bd:
            txt = ((txt + ". " if txt else "") + "Blind, DEM only (the survey is held for an independent check, "
                   "never compared here): " + ", ".join(bd))
        # under the legend, never over it: placed from the legend's drawn extent (the saved figure
        # grows to hold it, bbox_inches='tight')
        fig.canvas.draw()
        bb = leg.get_window_extent().transformed(fig.transFigure.inverted())
        fig.text(bb.x0, bb.y0 - 0.015, "\n".join(textwrap.wrap(txt, 150)), fontsize=8, color=INK2, va="top",
                 ha="left")
    fig.savefig(path, dpi=130, facecolor="white", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


class UsageParser(argparse.ArgumentParser):
    """argparse with this script's own exit code for a usage error (an unknown option, a bad value):
    EXIT_CODES['usage'] (6), not argparse's 2, which here means a disabled date."""

    def error(self, message):
        self.print_usage(sys.stderr)
        self.exit(EXIT_CODES["usage"], f"{self.prog}: error: {message} (NOTHING built)\n")


def main():
    ap = UsageParser(description=__doc__.split("\n")[1])
    what = ap.add_mutually_exclusive_group(required=True)
    what.add_argument("--date", help="survey date to build, YYYY-MM-DD (a date in --config)")
    what.add_argument("--all", action="store_true", help="every enabled date")
    what.add_argument("--summary", action="store_true", help="table and figure across the dates built")
    ap.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT),
                    help=f"products go in <root>/<date>/ (default {DEFAULT_OUTPUT_ROOT})")
    ap.add_argument("--share", action="store_true",
                    help="with --date: pack the date's product under --output-root into one archive to send "
                         "(<share-dir>/<name>.tgz + .sha256; on a BLIND date named by its calibration: "
                         "<date>_A_official_calibration.tgz, <date>_B_gcpfree.tgz); nothing is built")
    ap.add_argument("--share-dir", default=None, help="where --share writes (default <output-root>/share)")
    ap.add_argument("--a-root", default=None,
                    help="a GCP-free (B) build of a blind date: the output root of the table-calibration (A) product, "
                         "whose DEM build time B's README states (default: the root B's root is named after, e.g. "
                         f"{DEFAULT_OUTPUT_ROOT} for {DEFAULT_OUTPUT_ROOT}_gcpfree); only times are read")
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
    ap.add_argument("--eo", nargs="+", default=None, metavar="CAM=PATH",
                    help="build with these calibrations (EO files, each path taken as given) instead of the table's "
                         "eo_file, one per camera the date builds (a camera the GCP-free procedure did NOT DELIVER may "
                         "be left out: its report beside the given B files says so), e.g. "
                         "c2=calibration/gcpfree/runs/<run id>/CACO03_c2_gcpfree-oct2024_EO.yaml; --date only, with "
                         "--calibration-variant NAME, into a root named <root>_NAME, e.g. "
                         f"{DEFAULT_OUTPUT_ROOT}_gcpfree (any other root is refused, exit 5). A file whose "
                         "stated acceptance is FAIL ('# acceptance: FAIL', 'NOT DELIVERED') is refused; a missing one too (never "
                         "the table's instead)")
    ap.add_argument("--calibration-variant", default=None, metavar="NAME",
                    help="short name of the --eo calibration (e.g. gcpfree: lower case, digits, '-'): said in the "
                         "README, provenance, the DEM page and --summary, and names its root (<root>_<NAME>)")
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
                    help=f"wave-setup coefficient C (default {SETUP_COEF}, the live station's; 0 = none). A C other "
                         f"than the one a date was built with would REPLACE that build in <root>/<date>/ (and rerun a "
                         f"2025 date's detection), so such a run is refused without --replace: give each C its own "
                         f"--output-root, e.g. {DEFAULT_OUTPUT_ROOT}_c0")
    ap.add_argument("--replace", action="store_true",
                    help="build into a date folder whose waterlines carry another C and REPLACE that product. "
                         "Without it such a run is refused (exit 5): another C goes in its own --output-root "
                         f"({DEFAULT_OUTPUT_ROOT}_c<C>; the C = {SETUP_COEF} products in {DEFAULT_OUTPUT_ROOT})")
    ap.add_argument("--setup-fitted-to", default=None, metavar="SURVEY or none:HOW",
                    help="where a --setup-coef other than the known fits was fitted: the survey file it was "
                         "fitted to (comparisons with that survey become CIRCULAR), or 'none:<how>' for a C fitted "
                         "to no survey (e.g. 'none:repeat crossings, dem_from_contours.py --fit-setup on "
                         "2025-01-18..23'). Required for such a C: a C of unknown origin is refused")
    ap.add_argument("--geotiff-epsg", type=int, default=GEOTIFF_EPSG,
                    help=f"horizontal EPSG code written into the GeoTIFFs (default {GEOTIFF_EPSG}: the grids are in "
                         f"the frame of the GCPs, assumed NAD83(2011) / UTM 19N as the surveys; 32619, WGS 84 / UTM "
                         f"19N, was the default until Oct 2026 and is a NOMINAL tag, ~1-1.5 m off in a GIS)")
    ap.add_argument("--keep-detector-debug", action="store_true",
                    help="keep the detector's debug images (~8 MB a photo; default: removed once the "
                         "waterlines are written)")
    ap.add_argument("--keep-moved-days", action="store_true",
                    help="keep the days the pointing check marks DIFFERENT POINTING or UNCONFIRMED (default: "
                         "left out of the detection)")
    ap.add_argument("--correct-pointing", action="store_true",
                    help="OPT-IN, an evaluation: a day whose matched features pass pointing_daily.py's gates (a "
                         f"turn over {MATCH_ERR_DEG:g} deg, held through the day) is detected and its waterlines "
                         "mapped with its own pointing (the calibration turned by the change measured that day) "
                         "instead of being left out; every other day as by default. Into "
                         f"its own root only, <root>_{POINTING_CORRECTED} (e.g. {DEFAULT_OUTPUT_ROOT}_"
                         f"{POINTING_CORRECTED}; any other root is refused, exit 5); not on a BLIND date, not with --eo")
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
    ap.add_argument("--dem-loess", action="store_true",
                    help="also write the scale-controlled loess surface of the DEM's points and its prediction error "
                         "(dem_from_contours.py --loess, dem_loess.py: <stem>_loess, _loess_err, _loess_nmse, "
                         "_loess_skill .asc and .tif). Off by default; the DEM, its page, maps and comparisons are "
                         "the same with it. Check its error against the date's survey with python3 dem_loess.py "
                         "evaluate dem/<date> <survey file>")
    args = ap.parse_args()

    arm_blind_guard()          # no file named like a held survey (BLIND_DATES) is opened by this run, whatever asks
    if args.eo or args.calibration_variant:
        if not args.date:
            ap.error("--eo / --calibration-variant: with --date only (one date, one calibration)")
        if not (args.eo and args.calibration_variant):
            ap.error("--eo and --calibration-variant go together: the calibrations, and the name of that variant")
        if not VARIANT_NAME.match(args.calibration_variant) or args.calibration_variant in ("table", "station"):
            ap.error(f"--calibration-variant {args.calibration_variant!r}: a short name of lower-case letters, digits "
                     f"and '-' (not 'table' or 'station'), e.g. gcpfree")
        eo = {}
        for item in args.eo:
            cam, sep, path = item.partition("=")
            if not sep or cam not in ("c1", "c2") or not path or cam in eo:
                ap.error(f"--eo {item!r}: expected c1=PATH and/or c2=PATH, each camera once")
            eo[cam] = os.path.abspath(path)     # recorded as an absolute path: the stamps and the rebuild line
        args.eo = eo
    if args.summary:
        return summary(args)
    if args.share:
        if not args.date:
            ap.error("--share: with --date (one date's product, from --output-root)")
        return share(args)
    # a shutdown or kill (SIGTERM) ends the run like Ctrl-C: the build writes 'status: failed' and why
    # instead of leaving the previous README looking current
    import signal

    def _terminated(signum, frame):
        raise SystemExit(f"terminated by signal {signum} (SIGTERM)")
    try:
        signal.signal(signal.SIGTERM, _terminated)
    except (ValueError, OSError):
        pass
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
    if args.setup_coef and not setup_fit_known(args.setup_coef) and not args.setup_fitted_to and not args.dry_run:
        # refused with exit 5 ('refused'): the option is valid, the run's settings are what is refused
        print(f"{ap.prog}: error: --setup-coef {args.setup_coef}: where was this C fitted? Give --setup-fitted-to "
              f"<survey file> (the survey it was fitted to) or --setup-fitted-to 'none:<how>' (fitted to no "
              f"survey). A C of unknown origin could be fitted to the very survey it is compared with. NOTHING "
              f"built", file=sys.stderr)
        return EXIT_CODES["refused"]
    if args.setup_fitted_to and setup_fit_known(args.setup_coef):
        warn(f"--setup-fitted-to ignored: C = {args.setup_coef} is a known fit (SETUP_FITS)")
    if args.nice and not args.dry_run:
        try:
            level = os.nice(int(args.nice))
            say("priority", f"niceness {level} (--nice {args.nice}: the station's cron keeps priority)")
        except (OSError, AttributeError) as exc:
            warn(f"--nice {args.nice} not applied: {exc}")
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
    cron_c = cron_settings(args.live_cron).get("SETUP_COEF")
    try:
        cron_c = float(cron_c) if cron_c not in (None, "") else 0.0
    except ValueError:
        cron_c = None
    if cron_c is None or abs(cron_c - SETUP_COEF) > 1e-9:
        warn(f"the station's cron ({args.live_cron}) has SETUP_COEF {cron_settings(args.live_cron).get('SETUP_COEF')!r}, "
             f"not this script's default C = {SETUP_COEF}: the live waterlines and these products would carry "
             f"different setups. Change SETUP_COEF in survey_products.py, historical_forcing.py and the cron "
             f"together (with a SETUP_FITS entry saying where the new C was fitted)")
    # the cron's opt-in slope-aware setup (SETUP_MODEL="stockdon", setup_model.py): read here only, not in
    # cron_settings (whose keys go into forcing.json); these products keep C's setup either way
    try:
        cron_model = next((m.group(1) for m in (re.match(r'^\s*SETUP_MODEL="?([A-Za-z]*)"?', ln) for ln in
                                                Path(args.live_cron).read_text().splitlines()) if m), "")
    except OSError:
        cron_model = ""
    if cron_model == "stockdon":
        warn(f"the station's cron ({args.live_cron}) runs the slope-aware setup (SETUP_MODEL=\"stockdon\"): its live "
             f"DEMs carry Stockdon's setup, these products C = {args.setup_coef}'s (a live date's rows are recomputed "
             f"with C), so the two differ by the models' difference")
    if args.date:
        res = build_date(args.date, cfg, srv, args, steps)
        return EXIT_CODES.get(res, 1)
    results = {}
    for d in sorted({r["date"] for r in cfg}):
        results[d] = build_date(d, cfg, srv, args, steps)
    return finish_all(results, args)


def finish_all(results, args):
    """The end of an --all run: each date's result, then the summary table across the root. The exit code is the
    builds' worst (worst_exit): a summary that cannot be made is said (with its traceback), never the run's code."""
    rule("all dates")
    for d, r in results.items():
        say(d, r + {"partial": " (PARTIAL: see its README)", "failed": " (FAILED: see its log)",
                    "refused": " (NOT built: see the WARNING above)"}.get(r, ""))
    if not args.dry_run:
        try:
            summary(args)
        except Exception:                          # the builds are done and on disk: their result stands
            import traceback
            traceback.print_exc()
            warn("the summary table across the root could not be made (traceback above); every date's product is "
                 "as its own README says. Run --summary again after a fix")
    return worst_exit(results.values())


# Exit codes: 0 built (or dry run), 1 a step failed (or a configuration table cannot be read), 2 the
# date is disabled, 3 built but PARTIAL (a camera contributed nothing, a part is missing: a camera's
# pointing check or photo map, or the outputs on disk are not one build), 4 not started: another build
# of the date is running, 5 not started: refused by the run's own settings (a C other than the one the
# date's waterlines carry, without --replace; a --setup-coef of unknown origin, without
# --setup-fitted-to), 6 not started: a usage error (an unknown option or a bad value, e.g. --steps,
# --window, --utc-hours, a negative --setup-coef; UsageParser), so a caller can tell it from a
# disabled date.
def share_name(date, prov):
    """The archive name of a date's product: on a BLIND date the calibration says which DEM it is
    (<date>_A_official_calibration for the table's, <date>_B_gcpfree for the GCP-free one: SHARE_NAMES), else
    <date>_survey_product[_<variant>]. -> name (no extension)."""
    var = sorted({(v.get("calibration") or {}).get("variant") or "table" for v in (prov.get("cameras") or {}).values()})
    v = var[0] if len(var) == 1 else "/".join(var)
    if prov.get("blind"):
        return f"{date}_" + SHARE_NAMES.get(v, f"{v}_calibration")
    return f"{date}_survey_product" + ("" if v == "table" else f"_{v}")


def lens_note(io_, cam, station, files, seen):
    """The detector reads every date's lens as CACO05_<cam>_20240801_IO.yaml (the lens of the camera since Aug 2024,
    on every mount). A product of another station (CACO03, CACO04) also gets that station's file of the same lens
    (e.g. CACO03_<cam>_20240801_IO.yaml, which gcpfree_eo.py solved with), and the MANIFEST says whether the two
    have the same numbers. -> MANIFEST line or ''."""
    if not station or station == "CACO05" or not str(io_.name).startswith("CACO05_"):
        return ""
    other = io_.parent / io_.name.replace("CACO05_", f"{station}_", 1)
    if not other.is_file():
        return (f"lens {cam}: {io_.name} (the detector's lens file; no {other.name} here to compare with)")
    try:
        from georectify import load_intrinsics
        dmax = float(np.max(np.abs(np.asarray(load_intrinsics(io_), float) - np.asarray(load_intrinsics(other), float))))
    except Exception:
        dmax = None
    if other not in seen:
        seen.add(other)
        files.append((other, f"calibration/{other.name}"))
    return (f"lens {cam}: {io_.name} is the file the detector reads (every station's lens under its CACO05 name); "
            f"{other.name}, this {station} product's own lens file, is in the archive too: "
            + (f"the same lens (the CACO05 file holds its numbers rounded: the largest difference {dmax:.1e}, in px "
               f"or a distortion coefficient)" if dmax is not None and dmax < 1e-5 else
               f"THEY DIFFER (by up to {dmax:.3g}): this DEM was built with the CACO05 file" if dmax is not None else
               "not compared (not readable)"))


def reasons_text(reasons):
    """Reasons joined with '; ', each once, in their first order (two cameras of one B share its build record)."""
    return "; ".join(dict.fromkeys(r for r in reasons if r))


def share(args):
    """
    One archive of a built date to send (e.g. to the holder of a blind date's survey): README, provenance, the
    DEM GeoTIFFs and page, the maps, the calibration files used (with, for a GCP-free calibration, its
    procedure's reports, run record and pre-registration, and the report of a camera it did not deliver), and
    MANIFEST.txt (what the archive is, the calibration identity, the holder's-check label, whether the
    pre-registration is still the file the run recorded, every file's sha256). Named so the two DEMs of a blind
    date cannot be mixed up once separated from their folders: 2024-10-23_A_official_calibration.tgz and
    2024-10-23_B_gcpfree.tgz. Never written over (--replace to rebuild it). Refused (exit 5) for a product in
    progress or failed, and for a GCP-free variant whose calibration is no longer the pre-registered procedure's.
    -> exit code.
    """
    import tarfile
    out = Path(args.output_root) / args.date
    prov = load_json(out / "provenance.json")
    if not prov:
        print(f"--share: no product of {args.date} under {args.output_root} (no provenance.json): nothing shared")
        return EXIT_CODES["refused"]
    if prov.get("status") not in ("complete", "partial"):
        print(f"--share: the product of {args.date} is {prov.get('status')}: only a finished build is shared")
        return EXIT_CODES["refused"]
    name = share_name(args.date, prov)
    cams = prov.get("cameras") or {}
    var = {(v.get("calibration") or {}).get("variant") or "table" for v in cams.values()}
    problems = []
    for cam, v in sorted(cams.items()):
        cal = v.get("calibration") or {}
        # the calibration file packed must be the one this DEM was built with (provenance; the detection's stamp)
        q = Path(str(cal.get("eo_path") or ""))
        now = sha256(q) if cal.get("eo_path") and q.is_file() else None
        st = read_stamp(out, f"detect_{cam}") or read_stamp(out, f"pointing_{cam}") or {}
        built = str((st.get("calibration") or {}).get("eo_sha256") or st.get("eo_sha256") or "")
        rec = str(cal.get("eo_sha256") or "")
        if now is None:
            problems.append(f"{cam} {cal.get('eo_file')}: the calibration file it was built with is not here ({q})")
        elif (rec and not now.startswith(rec)) or (built and built != now):
            problems.append(f"{cam} {cal.get('eo_file')}: the file there now (sha256 {now[:16]}) is not the one this "
                            f"product was built with (provenance {rec[:16] or '?'}, detection {built[:16] or '?'}): an "
                            f"archive never pairs a DEM with another calibration file (restore the file, or rebuild)")
        if cal.get("variant") in GCPFREE_VARIANTS:
            pc = procedure_check(cal.get("eo_path"), eo_dirs(args)) if cal.get("eo_path") else {}
            if not pc.get("pinned"):
                problems.append(f"{cam} {cal.get('eo_file')}: NOT the pre-registered procedure's B "
                                f"({'; '.join(pc.get('why') or ['its file is not here'])})")
            elif prov.get("blind"):
                # the one B root its first build recorded (B_BUILD.txt, committed once): no other root is shared
                rec = b_build_record({cam: pc})
                g = rec.get("git") or {}
                if not rec.get("exists"):
                    problems.append(f"{cam}: no {B_BUILD_FILE} in its run folder ({rec.get('path')}): a B is shared "
                                    f"only from the root its first build recorded there")
                elif rec.get("root") != os.path.realpath(str(args.output_root)):
                    problems.append(f"{cam}: this B is built in {rec.get('root')} ({rec['path']}), not in "
                                    f"{args.output_root}: only that root's B is shared")
                elif not g.get("committed") or g.get("commits") != 1:
                    problems.append(f"{rec['path']} is {g.get('why') or 'committed more than once'}: commit it (once) "
                                    f"before the B is shared")
    for v_ in sorted(var):
        if variant_lookalike(v_):
            problems.append(variant_lookalike(v_))
    if problems:
        print(f"--share: {name} REFUSED: " + reasons_text(problems))
        return EXIT_CODES["refused"]
    sdir = Path(args.share_dir or (Path(args.output_root) / "share"))
    tgz = sdir / f"{name}.tgz"
    if tgz.exists() and not args.replace:
        print(f"--share: {tgz} exists: an archive is never written over (--replace to make it again)")
        return EXIT_CODES["refused"]
    files = []                                     # (path on disk, path in the archive)
    prereg = []                                    # what the archive says of the GCP-free pre-registration
    lens = []                                      # what the archive says of the lens files
    for f in ("README.txt", "provenance.json"):
        if (out / f).is_file():
            files.append((out / f, f))
    stem = next(iter(sorted((out / "dem").glob("*_dem.asc"))), None) if (out / "dem").is_dir() else None
    if stem is not None:
        st = str(stem)[:-len("_dem.asc")]
        for suf in ("_dem.tif", "_spread.tif", "_count.tif", "_sigma.tif", "_dem.png", "_map.png", "_info.json") + \
                tuple(f"_{k}.tif" for k in DEM_LOESS_GRIDS):          # the loess ones only when --dem-loess built them
            if Path(st + suf).is_file():
                files.append((Path(st + suf), f"dem/{Path(st + suf).name}"))
    for f in sorted((out / "maps").glob("*.png")) if (out / "maps").is_dir() else []:
        files.append((f, f"maps/{f.name}"))
    # what the README's FILES section and caveats point to, so the holder can check the forcing, the pointing and
    # the DEM day by day from the archive alone (review r1): the forcing, each camera's pointing table, the
    # waterlines the DEM used with the filter's report, the comparison (none on a blind date) and the step logs.
    # Not packed: waterlines/<cam>/ (the detection's work folders, with copies of the photos), the other waterline
    # file, pointing/reference_<cam>/ (a copy of the reference photo) and logs/stamps/ (build records)
    for sub, pats in (("forcing", ("forcing_report.txt", "forcing.json", "water_level.csv", "waves.csv", "forcing.png")),
                      ("pointing", tuple(f"pointing_{c}.{e}" for c in sorted(cams) for e in ("csv", "json"))),
                      ("waterlines", ("consistency_report.csv", "no_setup_frames.csv")
                       + tuple(f"consistency_{c}.png" for c in sorted(cams))),
                      ("compare", ("*",)), ("logs", ("*.log",))):
        for pat in pats:
            for f in sorted((out / sub).glob(pat)) if (out / sub).is_dir() else []:
                if f.is_file() and (f, f"{sub}/{f.name}") not in files:
                    files.append((f, f"{sub}/{f.name}"))
    wl = waterline_file(product_plan(prov, out), quiet_=True)
    if wl.is_file():
        files.append((wl, f"waterlines/{wl.name}"))
    # a --correct-pointing product: the per-day pointings its waterlines were mapped with, with the calibration files
    # (none in a default build)
    for f in sorted((out / "waterlines").glob("c*/pointing_corrected_eo/*.yaml")) if (out / "waterlines").is_dir() else []:
        files.append((f, f"calibration/pointing_corrected_eo/{f.name}"))
    seen = set()
    for cam, v in sorted(cams.items()):
        cal = v.get("calibration") or {}
        for key in ("eo_path",):
            q = Path(str(cal.get(key) or ""))
            if q.is_file() and q not in seen:
                seen.add(q)
                files.append((q, f"calibration/{q.name}"))
        io_ = Path(args.calibration) / str(cal.get("io_file") or "")
        if io_.is_file() and io_ not in seen:
            seen.add(io_)
            files.append((io_, f"calibration/{io_.name}"))
            lens.append(lens_note(io_, cam, (v.get("station") or ""), files, seen))
        pc = ((cal.get("origin") or {}).get("procedure") or {})
        if pc.get("pinned") and pc.get("run_dir"):
            rd = Path(pc["run_dir"])
            for q in sorted(rd.glob("*_report.json")) + sorted(rd.glob("*_report.txt")) + [rd / "run.json",
                                                                                           rd / "run_report.txt",
                                                                                           rd / B_BUILD_FILE]:
                if q.is_file() and q not in seen:
                    seen.add(q)
                    files.append((q, f"calibration/gcpfree_run_{rd.name}/{q.name}"))
            # the GCP solution(s) B was turned from (its report's b_reference_eo and Jan 2025 reference_eo), so the
            # holder can check B's 'rotated by', 'pan change' and 'Position kept from' lines from the archive alone
            rep_ = load_json(Path(str(pc.get("report") or ""))) or {}
            for key in ("b_reference_eo", "reference_eo"):
                nm, want_ = str(rep_.get(key) or ""), str(rep_.get(key + "_sha256") or "").lower()
                q = next((Path(c_) / nm for c_ in (rd, Path(args.calibration), HERE / "calibration")
                          if nm and (Path(c_) / nm).is_file()), None)
                if nm and (q is None or q in seen):
                    if q is None:
                        prereg.append(f"reference EO {nm} (B's report: sha256 {want_[:16]}): not found here, not packed")
                    continue
                if q is not None:
                    seen.add(q)
                    files.append((q, f"calibration/{q.name}"))
                    prereg.append(f"reference EO {nm} (what B was turned from, {key}): sha256 {file_sha(q)[:16]}, "
                                  + ("the file B's report recorded" if file_sha(q) == want_ else
                                     f"NOT the file B's report recorded ({want_[:16]})"))
            # the pre-registration (the outcome rule, fixed before any run): calibration/gcpfree/PREREGISTERED.txt,
            # two folders up from runs/<run id>/; said whether it is still the file the run recorded
            q = rd.parent.parent / "PREREGISTERED.txt"
            rec = str((pc.get("preregistration") or {}).get("sha256") or "")
            if q.is_file() and q not in seen:
                seen.add(q)
                files.append((q, f"calibration/{q.name}"))
                prereg.append(f"pre-registration: {q.name} sha256 {file_sha(q)[:16]}, "
                              + ("the file the procedure's run recorded" if file_sha(q) == rec else
                                 f"NOT the file the run recorded (sha256 {rec[:16]}): changed after the run"))
            elif not q.is_file():
                prereg.append(f"pre-registration: {q} not found (the run recorded sha256 {rec[:16]})")
    for cam, v in sorted((prov.get("cameras_not_delivered") or {}).items()):
        q = Path(str(v.get("report") or ""))
        if q.is_file() and q not in seen:
            seen.add(q)
            files.append((q, f"calibration/gcpfree_run_{q.parent.name}/{q.name}"))
    bl = prov.get("blind") or {}
    ident = []
    for cam, v in sorted(cams.items()):
        cal = v.get("calibration") or {}
        ident.append(f"  {cam}: {cal.get('eo_file')} (variant {cal.get('variant') or 'table'}, sha256 "
                     f"{file_sha(cal['eo_path'])[:16] if cal.get('eo_path') and Path(cal['eo_path']).is_file() else '?'})"
                     + (f": {((bl.get('check') or {}).get('per_camera') or {}).get(cam, {}).get('label')}" if bl else ""))
    lines = [f"{name}: the survey-date product of {args.date}, as built under {args.output_root}",
             f"status: {prov.get('status')}" + (" (PARTIAL: see README NOT COMPLETE)" if prov.get("status") == "partial"
                                                else ""),
             ("SYNTHETIC FIXTURE: tests the pipeline, not the beach" if prov.get("synthetic") else ""),
             f"calibration: {'/'.join(sorted(var))}" + (" (the table's, survey_dates.csv)" if var == {"table"} else
                                                         " (--eo; NOT the table's)")]
    if bl:
        lab_ = (bl.get("check") or {}).get("label")
        lines += [f"BLIND: {bl.get('statement')}",
                  f"his check of this DEM: {lab_}"
                  + (" (procedure made with stated exposure to A: README EXPOSURE, and below)"
                     if lab_ == "INDEPENDENT" and (bl.get("gcpfree") or {}).get("exposure") else "")
                  + " (per camera below; README: BLIND)"]
    lines += ["cameras:"] + ident
    for cam, v in sorted(cams.items()):
        cal = v.get("calibration") or {}
        pc = (cal.get("origin") or {}).get("procedure") or {}
        if pc.get("pinned"):
            lines.append(f"  {cam} procedure: {procedure_text(pc)}")
            lines.append(f"  {cam} STATED UNCERTAINTY (pre-registered rule): {uncertainty_terms(pc)}")
            lines.append(f"  {cam} {tilt_effect_text(cal.get('rotation_term_as_tilt'))}")
    for cam, v in sorted((prov.get("cameras_not_delivered") or {}).items()):
        lines.append(f"  {cam}: NOT DELIVERED by the pre-registered GCP-free procedure (report "
                     f"{Path(str(v.get('report'))).name}); no {cam} data in this DEM")
    gf = bl.get("gcpfree") or {}
    if gf:
        lines += [gf.get("exposure") or "", "ORDER (what decides B, and A): " + "; ".join(gf.get("order") or [])]
        lines += [f"NOTE: {x}" for x in gf.get("order_flags") or []]
        lines += [f"{x}." for x in gf.get("protocol") or []]
    lines += sorted(set(prereg))
    lines += [x for x in lens if x]
    lines += ["The DEM and map files carry the calibration variant in their names (the table's: <date>_dem.tif, "
              "maps/<date>_waterlines_plan.png; a variant: <date>_<variant>_dem.tif, maps/<date>_<variant>_waterlines_"
              "plan.png); every GeoTIFF's description (and, on a blind date or with a variant, the plan view's title) "
              "names the calibration (EO files, sha256) and, on a blind date, the GeoTIFFs the holder's-check label.",
              "What of the build's folder is here (README: FILES): README.txt, provenance.json, dem/, maps/, forcing/ "
              "(water level, waves, the forcing report), pointing/<cam>.csv and .json (each day's verdict and why), "
              f"waterlines/{wl.name} (the waterline points the DEM used, each with its frame's capture time: a check "
              "by day or by frame starts there) with the consistency filter's report, compare/ (none on a blind "
              "date), logs/*.log and calibration/. Not here: waterlines/<cam>/ (the detection's work folders, with "
              "copies of the photos), pointing/reference_<cam>/ (a copy of the reference photo), logs/stamps/.",
              "", "files (sha256, bytes, path):"]
    lines += [f"  {file_sha(src)}  {src.stat().st_size:>10}  {arc}" for src, arc in files]
    sdir.mkdir(parents=True, exist_ok=True)
    man = sdir / f"{name}_MANIFEST.txt"
    man.write_text("\n".join(x for x in lines if x is not None) + "\n")
    tmp = sdir / f".{name}.tgz.partial"
    with tarfile.open(str(tmp), "w:gz") as tf:
        tf.add(str(man), arcname=f"{name}/MANIFEST.txt")
        for src, arc in files:
            tf.add(str(src), arcname=f"{name}/{arc}")
    os.replace(str(tmp), str(tgz))
    man.unlink()
    h = hashlib.sha256(tgz.read_bytes()).hexdigest()
    (sdir / f"{name}.tgz.sha256").write_text(f"{h}  {tgz.name}\n")
    print(f"wrote {tgz} ({len(files) + 1} files, {tgz.stat().st_size} bytes)\nsha256 {h}  (also in {tgz.name}.sha256)")
    if prov.get("status") == "partial":
        warn(f"{args.date}: the product is PARTIAL (README: NOT COMPLETE): said in the archive's MANIFEST")
    return 0


EXIT_CODES = {"built": 0, "dry-run": 0, "disabled": 2, "failed": 1, "partial": 3, "busy": 4, "refused": 5,
              "usage": 6}
# --all returns the code of its WORST date by severity, not the largest number: a failed date (1)
# outranks one refused (5) or not started because another build held it (4), which outrank a partial
# one (3). Disabled dates (2) are skipped by --all and do not count.
SEVERITY = ("failed", "refused", "busy", "partial")


def worst_exit(results):
    """Exit code of an --all run: the most severe per-date result
    (failed > refused > busy > partial > built)."""
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

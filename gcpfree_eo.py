#!/usr/bin/env python3
"""
GCP-Free Pointing For 23 Oct 2024 (Calibration B), Pre-Registered
==================================================================
Writes calibration B for the 23 Oct 2024 DEM that Dr Chris Sherwood (USGS)
checks BLIND against his 23 Oct 2024 GCP survey. The station's official
calibration A (CACO03_<cam>_20241023_EO.yaml) was solved from those same
GCPs; B must owe nothing to them. This script never opens calibration A, a
GCP or target file, or calibration_history.csv (guard() refuses them, by
name and by the name a link points to; a photo-named link to one is never
even probed), and nothing in it was set from A's numbers by intent (what
its makers had seen of A: EXPOSURE below): the procedure, photo rules and
thresholds were fixed on 13-14 Oct 2026 from synthetic tests and one real
Jan 2025 frame, BEFORE any run on the real Oct 2024 photos. Do not change a
constant after seeing a result; a changed procedure is a new, separately
labelled product (and survey_products.py, which pins this file's sha256,
its modules' and its pre-registered constants', calls it 'not the
pre-registered procedure').

WHAT B IS. The 23 Jan 2025 GCP solution of the same CACO03 mount
(CACO03_<cam>_20250123_EO.yaml) turned by the rotation measured between Jan
2025 photos and Oct 2024 photos (estimate_eo_rotation.py's method: SIFT on
land, rays through the lens model, a rotation fitted with RANSAC). The
camera POSITION is kept from that solution: B inherits it. The CACO03 GCP
solutions of Aug 2024 and Jan 2025 place the camera 0.36 m apart for c1
(0.34 m horizontally, 0.11 m in height) and 1.23 m for c2 (1.08 m, 0.59 m)
on a mount that did not move.

PRE-REGISTRATION FILE (required). calibration/gcpfree/PREREGISTERED.txt next
to this script (PREREG; there is no other), written by `--preregister` with
the OWNER's choice of the outcome rule (below) and the PHOTO MANIFEST (name,
size and sha256 of every Oct, Jan and, for c2, Dec photo a run selects from
the station's three folders), COMMITTED to git and PUSHED (git push) before
any run. EVERY run, a trial too, is made under this file and refuses (exit 2)
while it is missing, not committed, changed since its last commit, its commit
on no remote-tracking branch (not pushed: the time GitHub receives the push is
the rule's only public timestamp), names another sha256
of this script or of its modules (estimate_eo_rotation.py, georectify.py,
horizon_check.py, detect_original_view.py, view_reproject.py), names photo
folders other than the station's three (ROOTS), pins no photos, or states no
valid rule; and once any run has been made under it, under any other version
of it. So the rule is fixed, committed, before this code solves any photo,
and a run on the station's folders that selects other photos than the
manifest (a photo removed, replaced or added after a look) is refused. Each
run records the file's sha256 and its git history (every commit that touched
it). Nothing is pre-registered or run unless this file and its modules are
the version survey_products.py pins last (its GCPFREE_PINS, read without
importing it): a B of any other version would be 'NOT the pre-registered
procedure' there, and the one delivery run would be spent on it. A run is
also refused while a commit OUTSIDE this branch's history (another branch, a
stash, the reflog: a branch made from before the pre-registration, a reset)
wrote a PREREGISTERED.txt or a run record (runs/*/run.json) whose content
this branch's history does not hold: an earlier pre-registration or run there
would otherwise be invisible. What no check here can see: a copy of this
procedure run outside this repository (another clone, a copy of the script)
before the pre-registration, and a branch or reflog entry removed with its
history; only the sha256s the holder is sent beforehand (the deploy notes)
answer those. survey_products.py accepts a B only from runs/ next to this
repository's file, as the first delivery run there, and its README lists
every run in runs/ (a trial, a killed run) with its start time.

STEPS (per camera)
  1. Photos (CACO03 timex, full frame as the lens file, each name once):
     Oct = 20-26 Oct 2024 (the DEM window); Jan = 18-23 Jan 2025 (the
     calibration day and the five days before it; 24 Jan is the re-set day);
     Dec 2024 (c2, all days); Aug = 25-29 Aug 2024, only from --aug-dir. Per
     day up to 3 photos: the nearest to 14, 17 and 20 UTC within 1 h; mean
     grey 40-245 (no night, no white-out).
  2. Reference check: every Jan day is matched to 23 Jan (the anchor); a day
     MEASURED to point more than REF_DAY_TOL_PX away is left out; a day the
     features cannot measure (snow, fog) stays unless its horizon differs
     from the anchor's by more than REF_HZ_TOL_DEG.
     Then, BEFORE any Oct photo is read, the reference EO is fitted to the
     visible sea horizon of its own kept photos and that offset is logged
     (the c1 Jan 2025 solution puts the horizon ~1.6 deg of tilt off the
     visible one on a real frame; c2's offset is measured here first).
  3. DIRECT solve: Oct photos x kept Jan photos of the SAME target hour
     (similar sun), pooled. Gate F: >= MIN_INLIERS inliers, from >=
     MIN_NEW_DAYS Oct days and >= MIN_REF_PHOTOS reference photos, spread
     over >= MIN_EXTENT_PX across and down, residual median <=
     MAX_RESID_MED_PX, leave-one-day-out change <= JACKKNIFE_MAX_PX, and no
     second consensus of >= AMBIGUITY_RATIO of the first.
  4. Per-day solves (each Oct day alone). Gate S: >= MIN_DAYS_MEASURED days
     with >= DAY_MIN_INLIERS inliers, each within DAY_TOL_PX of the pooled
     rotation.
  5. Horizon, DIFFERENTIAL, judged in a BAND (tilt and roll only, never pan).
     d = fit(B on the Oct photos) - fit(reference EO on its own photos).
     A perfect B does not give d = 0: it carries the reference EO's own
     horizon error turned by the rotation Q. So the expected bias is
     measured: e = fit(B on the kept reference photos warped by Q) - fit(ref
     EO on the same photos unwarped). The check passes when d lies in the
     band [min(0, e) - tol, max(0, e) + tol] for tilt (tol HZ_TILT_TOL) and
     roll (tol HZ_ROLL_TOL). On the turned photos the horizon is looked for
     near where the reference EO, corrected by its own measured offset and
     turned by Q, puts it (exact for a rotation), and only in frame columns
     where the turned photos have reference content there (the side of a
     turned frame has none); the Oct fit uses the same columns. The band, the
     bias and a geometric prediction of it (the reference offset as a
     rotation, turned by Q) are reported. 'unavailable' when fewer than
     HZ_MIN_CLEAR clear horizons are found on any of the three sets.
  6. c2 CHAIN Jan -> Dec -> Oct (a cross-check, or B when the direct solve
     fails): each link is judged as a route (below); B_chain within
     CHAIN_TOL_PX of B_direct.
  7. Aug (only with --aug-dir): Aug -> Oct with CACO03_<cam>_20240827_EO.yaml
     (Aug 2024 GCPs, not the Oct ones), anchor 27 Aug, the same rules and
     gates; B_aug keeps the Aug 2024 position. Compared with B on the beach.
  8. c1/c2: the two pan changes (B azimuth - Jan azimuth) within PAN_TOL_DEG.
  (Pointing differences in px: the largest image displacement over the
  lower two thirds of the frame, where the beach is.)

VERDICT (per camera; never a silent fallback to A). A ROUTE counts when its
solve passes gates F and S and its horizon check is not FAIL (pass or
unavailable; for the chain: both links). Routes in order: the direct
Jan->Oct solve; the Aug->Oct solve (only with --aug-dir); c2's chain.
  DELIVERED          the route used counts with its horizon PASSED, and every
                     cross-check that counts agrees: c2's chain within
                     CHAIN_TOL_PX of the direct B, B_aug within AUG_TOL_M on
                     the beach.
  DELIVERED-WIDER    the route used counts but its horizon is unavailable
                     (tilt and roll unchecked); or a cross-check disagrees
                     (the uncertainty is then at least the disagreement); or
                     B = c2's chain (two links: never better than WIDER).
  NOT DELIVERED      no route counts. No B file. A is never used instead.
  If both cameras pass and their pan changes differ by more than
  PAN_TOL_DEG, both are NOT DELIVERED.

OUTCOME RULE when only one camera passes (the owner's choice, in
PREREGISTERED.txt, outcome_if_one_camera):
  both-or-none   B is delivered only when both cameras are; otherwise no B
                 file is written for either (the passing camera's result is
                 WITHHELD: kept in its report), and Chris gets A alone with
                 the reports.
  per-camera     each passing camera's B is delivered; survey_products.py
                 builds a one-camera B DEM whose README names the other
                 camera NOT DELIVERED, with this procedure's report for it.

RUNS AND HISTORY. Every run writes into its own folder runs/<UTC run id>/
next to PREREGISTERED.txt (never over an earlier one): per camera the B
file (only when delivered), its report .json/.txt, and run.json (argv, photo
folders, camera, Aug folder, the pre-registration and every sha256; the
pre-registered constants as they were in effect, with their sha256, which
survey_products.py pins, so a run whose constants were changed in the
running process is not the pre-registered procedure; and every repository
module the run had loaded, with its sha256: the five MODULES are pinned,
the others -- compare_dem_survey.py, asc_to_geotiff.py and
dem_from_contours.py, imported by detect_original_view.py for its survey
envelope -- compute nothing of B and are recorded only), and the photos it
selected (name, size, sha256; for Aug photos from --aug-dir too, which no
pre-registration can list before they arrive). A DELIVERY run uses the
station's three photo folders (ROOTS, the default --roots), exactly the
pre-registered photos, and --camera both; any other run (one camera, other
--roots) is a TRIAL (folder suffix _trial), made under the same committed
pre-registration and recorded in the same runs/: reports only, no B file. A
trial is a FULL result (every verdict, every angle): it is never run on the
station (survey_products.py lists every run in B's README, and a trial or a
killed run left uncommitted stops the B DEM build). A run (trial or
delivery) is refused (exit 2) when a camera it solves has no Oct 2024 or no
Jan 2025 photo, and a delivery run when a photo folder is missing (an
unmounted disk would change the photos used): bad input never becomes a
verdict (--dry-run says beforehand whether a run would be refused). A photo
whose name, or the name of the file a link points to, is forbidden (guard)
is never opened, not even for its size. The first finished delivery run
fixes B. Only ONE more delivery run is ever allowed, and only to bring new
evidence: with --aug-dir (photos of 25-29 Aug 2024 from Chris) when the
first run used none and left a camera NOT DELIVERED. It solves only those
cameras; a camera that passed in the first run keeps that run's verdict and
B (written again into the second run's folder, 'solved in run <id>').
WHICH B IS DELIVERED (fixed now, before any run): the first delivery run's
B is "B" (built with --calibration-variant gcpfree); when it delivered no
camera at all, the second run's B is "B" instead; otherwise the second
run's B is a SECOND, separately labelled product (--calibration-variant
gcpfree-aug, in its own root) and the first B stays the delivered one. The
B file and its report name the variant; survey_products.py refuses any
other. Commit each run's folder (B file, reports, run.json: ONE commit, never
amended), push it, and write down its B sha256s AT ONCE, before the B DEM is built and
before anything is compared with the A DEM: survey_products.py builds and
shares a B only when its file, report and run.json are committed, unchanged,
each touched by one commit. (The A DEM may be built before this procedure
runs: nothing of it, and nothing of calibration A, reaches this procedure,
whose code is pinned by sha256 in survey_products.py and whose outcome rule
is pre-registered; when A was built relative to the pin, the
pre-registration and the run is stated in the B DEM's README.)

STATED UNCERTAINTY (two terms, on every delivered B; NEITHER is a validated
bound on real photos):
  rotation   max(ROT_SPREAD_FACTOR x the source solve's spread, ROT_FLOOR_PX)
             px on the beach part of the frame (spread: the larger of its
             leave-one-day-out change and its median per-day difference; the
             chain: the sum of its links), at least the chain-direct
             disagreement when that made B WIDER. The 4 px floor covers the
             true error of every delivered B in the synthetic sets (1.8-3.6
             px), but those errors were almost all roll (tilt 0.008-0.02 deg),
             which moves the beach little: the same px as TILT moves waterlines
             metres at 200-400 m. So no DEM figure is given for it here;
             survey_products.py states the DEM effect of this term taken as a
             tilt, per range band, from the date's own DEM and calibration.
  position   the inherited camera position: the largest DEM effect over the
             five synthetic position-inheritance cases (camera height
             +/-0.25 and +/-0.55 m, 0.5 m east and north, and the Aug 2024 -
             Jan 2025 solution difference): POSITION_TERM below. The largest
             of five cases, not a bound. These cases were chosen AFTER A's
             rows had been seen, and three of them match c1's A minus Jan
             2025 offset (EXPOSURE): c1's DEM term comes from the -0.55 m
             case.

EXPECTED OUTCOME (stated before the run, from the tests of 13 Oct 2026):
  c1: NOT DELIVERED from Jan photos (its Jan and Oct views share no fixed
      feature after a ~12.5 deg turn: only the horizon, which gives tilt and
      roll, never pan). c2: uncertain (the direct solve rests on a corner of
      bluff/dune that the real Jan frames show under snow). Aug 2024 photos
      (--aug-dir) are the route most likely to work.

AUG PHOTOS FROM CHRIS (--aug-dir). Named as the station names them,
'<epoch>.<Day>.<Mon>.<dd>_<HH>_<MM>_<SS>.GMT.<YYYY>.CACO03.<cam>.timex.jpg'
(e.g. 1724770800.Tue.Aug.27_15_00_00.GMT.2024.CACO03.c1.timex.jpg), full
frame (the lens file's size: a reduced copy is never used), 25-29 Aug 2024,
any folder layout. A folder that yields no such photo for a camera solved is
refused (exit 2), never run without it.

EXPOSURE (what the people who made this procedure had seen of calibration
A, and when; written into every report and B file, see EXPOSURE below).
  13 Oct 2026, 15:11 UTC, before the procedure was designed: a grep printed
    the 2024-10-23 rows of calibration_history.csv (calibration A). The
    repository's table notes then also said the Aug 2024 calibration differs
    from A by ~2.3-2.5 deg of azimuth.
  The synthetic test turns that validated the gates (pan -12.5 deg for c1,
    -13.4 deg for c2) are the sizes the brief gave: A minus the Jan 2025
    solution.
  13 Oct 2026, 15:16 UTC, five minutes after A's rows were printed: the five
    position-inheritance cases that size the stated POSITION term were
    written. Three of them match c1's A minus Jan 2025 camera offset
    (-0.39 m east, +0.41 m north, +0.55 m up): the +/-0.55 m height cases
    and the 0.5 m east-and-north case; the design log gives no other source
    for 0.55 m (c1's Aug 2024 - Jan 2025 height difference is 0.11 m; c2's
    0.59 m has its own case). c1's stated DEM position term (POSITION_TERM:
    |dz| median 0.011 m, p90 0.026 m) comes from the -0.55 m case: it was
    sized after A had been seen, and is in effect the DEM effect of A's own
    height offset from the Jan solution. It is a stated uncertainty only: no
    gate, threshold or rule uses it.
  14 Oct 2026: the implementer read a synthetic test fixture's notes that
    give A's azimuths to 0.1 deg; the same day (the order is not recorded)
    the horizon pass rule was changed from |d| <= tolerance (the 13 Oct
    design) to the band [min(0, e) - tol, max(0, e) + tol] (the design's
    adversarial check asked for it), and its two window parameters were set
    (HZ_EXPECTED_WINDOW_PX = 60, HZ_EDGE_MARGIN_PX = 25). Later changes of 14
    Oct 2026 (review) touched texts, records and refusals only: no gate,
    threshold, window or rule of the computation.
  Apart from those position cases, no constant, threshold or rule was set
  from A's numbers by intent (this is a statement, it cannot be proven), and
  the computation never reads A.
  The 23 Oct 2024 GCPs were never on this side at all. On the station, the
  deploy notes have this procedure's sha256 (survey_products.py
  GCPFREE_PINS) sent to the holder before A's DEM is built; A's DEM may be
  built and seen before the outcome rule is pre-registered and before B is
  run: the B DEM's README states the actual order, with dates.

Usage (on the station, from the repo folder):
    python3 gcpfree_eo.py --dry-run       # the pre-registration, the photos (would a run be refused? are they
                                          # the pinned ones?), the reference horizon; writes nothing
    python3 gcpfree_eo.py --preregister --outcome-if-one-camera <both-or-none|per-camera> --decided-by "NAME"
        (the rule and the photos a run will use; then commit calibration/gcpfree/PREREGISTERED.txt, before
        anything else, and send its sha256 to Chris)
    python3 gcpfree_eo.py                 # the delivery run, both cameras (~10-30 min each)
        (then commit calibration/gcpfree/runs/<run id>/ in ONE commit, and write down each B file's sha256
        and send it to Chris)
    python3 gcpfree_eo.py --aug-dir /path/to/aug_2024_photos      # the one allowed second run
    (--camera c1|c2 or other --roots make a TRIAL run: a full result, never run on the station)
Exit status: 0 B delivered for both cameras, each DELIVERED; 3 B delivered
in part, or a camera DELIVERED-WIDER; 4 nothing delivered (also a TRIAL
run); 2 bad input or refused.
"""

import os
import re
import sys
import json
import shlex
import hashlib
import argparse
import subprocess
from pathlib import Path
from datetime import datetime, timezone

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from georectify import (load_intrinsics, load_extrinsics, cirn_angles_to_R, pixel_to_ground,  # noqa: E402
                        undistort_uv)
import estimate_eo_rotation as EER                                                       # noqa: E402

CAL = HERE / "calibration"
OUT_DIR = CAL / "gcpfree"
PREREG = OUT_DIR / "PREREGISTERED.txt"
STATION = "CACO03"
ROOTS = ["/mnt/I2Rgus_Data/Chelsea_pics", "/mnt/I2Rgus_Data/waterline/composite",
         "/home/argus_user/owg_marconi/images"]
IO_FILE = "CACO03_{cam}_20240801_IO.yaml"
REF_EO_FILE = {"jan": "CACO03_{cam}_20250123_EO.yaml", "aug": "CACO03_{cam}_20240827_EO.yaml"}
# B's name carries no solve date (it is not a solution of 23 Oct 2024): 'oct2024' is the pointing it is for
OUT_EO = "CACO03_{cam}_gcpfree-oct2024_EO.yaml"
OUT_REPORT = "CACO03_{cam}_gcpfree-oct2024_report"
RUN_FILE = "run.json"
# the --calibration-variant survey_products.py builds a B with (see RUNS AND HISTORY): the delivered B, or a
# second, separately labelled product from the one allowed second run (Aug photos)
B_VARIANT, B2_VARIANT = "gcpfree", "gcpfree-aug"
# the code the procedure runs: its sha256s are pre-registered with this file's
MODULES = ("estimate_eo_rotation.py", "georectify.py", "horizon_check.py", "detect_original_view.py",
           "view_reproject.py")
# the archived content of every calibration file the procedure reads (a file of that name with other content,
# e.g. a link to calibration A, is refused before anything is solved)
ARCHIVED_SHA256 = {
    "CACO03_c1_20250123_EO.yaml": "c31e0f461e845a394ae62c45a1f524b691b25be64e4348a15599a478e37d5c77",
    "CACO03_c2_20250123_EO.yaml": "53a021843f9d11f5fe7e9fae5f303e5354b55849f28f8e3603715927e46a035e",
    "CACO03_c1_20240827_EO.yaml": "b942f527d17281f57087975ae7e44aa2de8af079d4b0f912cd87c5c38e435eb2",
    "CACO03_c2_20240827_EO.yaml": "aab35e9d42c88b89adc2268a837e1b68baa3d7c993c1690a55524aae5c1fb0b2",
    "CACO03_c1_20240801_IO.yaml": "83c7484e20d6dcdffbc6f79ba44b98071c9957e3b5a3bbef94e57025348584b1",
    "CACO03_c2_20240801_IO.yaml": "a7bac4285fd813b921451cac7bc461a76c9c17c58b4b7f5770297fe9b5d92417",
}
OUTCOME_RULES = ("both-or-none", "per-camera")

# ---- pre-registered: photos ---------------------------------------------------------------
WINDOWS = {"oct": ("2024-10-20", "2024-10-26"),     # the 23 Oct 2024 DEM window (surveys/survey_dates.csv)
           "jan": ("2025-01-18", "2025-01-23"),     # calibration day and the 5 days before; 24 Jan = re-set
           "dec": ("2024-12-01", "2024-12-31"),     # c2 only (c1 has no Dec photos)
           "aug": ("2024-08-25", "2024-08-29")}     # only from --aug-dir
ANCHOR_DAY = {"jan": "2025-01-23", "aug": "2024-08-27"}
TARGET_HOURS = (14.0, 17.0, 20.0)                   # UTC; nearest photo within MAX_DT_H of each
MAX_DT_H = 1.0
GREY_RANGE = (40.0, 245.0)

# ---- pre-registered: solve and gates ------------------------------------------------------
THRESH_PX = 3.0                 # inlier: ray angle below 3 px (as estimate_eo_rotation.py)
RANSAC_ITERS = 10000            # 2-ray samples (enough at 5 % inliers)
SEED = 20241023
MIN_INLIERS = 60
MIN_NEW_DAYS = 3
MIN_REF_PHOTOS = 2
MIN_EXTENT_PX = 50.0            # inliers must span (P10-P90) this much of the new frame across AND down: the
                                # synthetic c2 solves spanned 183 x 120 px (all on the bluff) and were right to
                                # 1.8-3.6 px; the spurious c1 consensus spanned 8 x 2 to 61 x 4 px
MAX_RESID_MED_PX = 1.0
JACKKNIFE_MAX_PX = 5.0
AMBIGUITY_RATIO = 0.5
AMBIGUITY_SEP_DEG = 0.1
DAY_MIN_INLIERS = 25
DAY_TOL_PX = 10.0
MIN_DAYS_MEASURED = 3
REF_DAY_TOL_PX = 4.0
REF_HZ_TOL_DEG = 0.10          # a reference day the features cannot measure: horizon within this of the anchor's
HZ_FIT_OK_PX = 8.0              # horizon_check.FIT_OK_PX: a clear horizon
HZ_MIN_CLEAR = 3
HZ_TILT_TOL = 0.25              # noise widening of the horizon band (deg), each side
HZ_ROLL_TOL = 0.35
CHAIN_TOL_PX = 10.0
SAME_SLOT_PAIRS = True          # only photos of the same target hour are matched (similar sun): a 14 h vs 20 h
                                # pair biased the synthetic c2 solve by 0.9 deg of roll (29 px)
PAN_TOL_DEG = 3.0
AUG_TOL_M = 2.0                 # B_aug vs B: median ground distance of beach points (z = 0.5 m). Gross check only:
                                # the Aug 2024 and Jan 2025 GCP solutions put the c2 camera 1.23 m apart and their
                                # angles absorb that, so two correct B's can differ by ~0.8 m there

# ---- pre-registered: the stated uncertainty -----------------------------------------------
ROT_SPREAD_FACTOR = 2.0
ROT_FLOOR_PX = 4.0              # >= the true error of every delivered B of the synthetic sets (1.84, 2.34, 3.0, 3.6 px;
                                # almost all roll: their tilt errors were 0.008-0.02 deg, so their tiny DEM errors say
                                # nothing of this term as a tilt: no DEM figure is stated for it here)
# the largest DEM effect over the five synthetic position-inheritance cases (height +/-0.25, +/-0.55 m; 0.5 m
# east and north; the Aug 2024 - Jan 2025 solution difference), B re-solved in each (Oct view, m): the largest of
# five cases, not a bound. The cases were written after A's rows had been seen (EXPOSURE): the +/-0.55 m and 0.5 m
# cases match c1's A minus Jan 2025 offset, and c1's DEM numbers here come from the -0.55 m case
POSITION_TERM = {"c1": {"dem_abs_median_m": 0.011, "dem_p90_m": 0.026, "horizontal_median_m": 0.29},
                 "c2": {"dem_abs_median_m": 0.018, "dem_p90_m": 0.086, "horizontal_median_m": 0.83}}
POSITION_SPREAD = {"c1": "0.36 m (0.34 m horizontally, 0.11 m in height)",
                   "c2": "1.23 m (1.08 m horizontally, 0.59 m in height)"}
UNCERTAINTY_NOT_VALIDATED = ("neither term is a validated bound on real photos: both come from synthetic tests (the "
                             "DEM effect of the rotation term taken as a tilt, per range band, is stated with the B "
                             "DEM by survey_products.py, from the date's own geometry)")

# ---- what the people who made the procedure had seen of calibration A (docstring: EXPOSURE) ------------------
EXPOSURE = (
    "EXPOSURE (stated by the people who made this procedure; it cannot be proven): on 13 Oct 2026 15:11 UTC, before "
    "the procedure was designed, a grep printed the 2024-10-23 rows of calibration_history.csv (calibration A, solved "
    "from the 23 Oct 2024 GCPs), and the repository's table notes said the Aug 2024 calibration differs from A by "
    "~2.3-2.5 deg of azimuth; the synthetic test turns that validated the gates (pan -12.5 deg for c1, -13.4 deg for "
    "c2) are the sizes the brief gave, i.e. A minus the Jan 2025 solution; at 15:16 UTC the same day, five minutes after "
    "A's rows were printed, the five position-inheritance cases that size the stated POSITION term were written, and "
    "three of them match c1's A minus Jan 2025 camera offset (-0.39 m east, +0.41 m north, +0.55 m up): the +/-0.55 m "
    "height cases and the 0.5 m east-and-north case (the design log gives no other source for 0.55 m; c1's Aug 2024 - "
    "Jan 2025 height difference is 0.11 m), and c1's stated DEM position term (|dz| median 0.011 m, p90 0.026 m) comes "
    "from the -0.55 m case: it was sized after A had been seen (a stated uncertainty only: no gate, threshold or rule "
    "uses it); on 14 Oct 2026 the implementer read a test "
    "fixture's notes giving A's azimuths to 0.1 deg, and the same day (the order is not recorded) the horizon pass "
    "rule was changed from |d| <= tolerance (the 13 Oct design) to the band [min(0, e) - tol, max(0, e) + tol] (asked "
    "for by the design's adversarial check) and its two window parameters were set (60 px, 25 px); later changes of "
    "14 Oct 2026 touched texts, records and refusals only. Apart from those position cases, no constant, threshold or "
    "rule was set from A's numbers by intent, and the computation never reads A; the 23 Oct 2024 GCPs were never on "
    "this side. On the station, the "
    "deploy notes have this procedure's sha256 sent to the holder before A's DEM is built; A's DEM may be built and seen "
    "before the outcome rule is pre-registered and before B is run: the B DEM's README states the actual order, with "
    "dates.")

FORBIDDEN = (re.compile(r"20241023_EO", re.I), re.compile(r"gcp(?!free)", re.I),
             re.compile(r"calibration_history", re.I), re.compile(r"Extrinsic_Targets", re.I),
             re.compile(r"2024-?10-?23.*(xyz|target|control|checkshot|rtk|survey|points)", re.I))

AUG_NAMING = ("'<epoch>.<Day>.<Mon>.<dd>_<HH>_<MM>_<SS>.GMT.<YYYY>.CACO03.<cam>.timex.jpg' (e.g. "
              "1724770800.Tue.Aug.27_15_00_00.GMT.2024.CACO03.c1.timex.jpg), full frame as the lens file, "
              "25-29 Aug 2024")


class Refused(Exception):
    """A run that must not go ahead (exit 2): the reason is printed."""


def forbidden_name(path):
    """The name that makes `path` forbidden (its own, or the name of the file a link points to), else None."""
    p = Path(path)
    for name in (p.name, Path(os.path.realpath(str(p))).name):
        if any(f.search(name) for f in FORBIDDEN):
            return name
    return None


def guard(path):
    """Every file this script reads goes through here: never calibration A, a GCP or target file, a survey of
    23 Oct 2024 or the calibration history -- by its name and by the name it really has (a link). Refused
    (exit 2) before it is opened."""
    p = Path(path)
    bad = forbidden_name(p)
    if bad:
        raise Refused(f"REFUSED to read {p.name}" + (f" (a link to {bad})" if bad != p.name else "")
                      + ": calibration B must not use calibration A, the Oct 2024 GCPs or the calibration history")
    return p


def sha256(path):
    h = hashlib.sha256()
    with open(guard(path), "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def own_sha256():
    return sha256(Path(__file__).resolve())


def module_sha256():
    return {m: sha256(HERE / m) for m in MODULES}


# The pre-registered constants (every number, name, folder and rule the computation and the verdict use), recorded
# in run.json as they were in effect in the running process, with their sha256, which survey_products.py pins with
# this file's (GCPFREE_PINS 'constants_sha256'): a run whose constants were changed in its process (other photo
# folders patched in, a threshold moved) is not the pre-registered procedure there.
PREREGISTERED_CONSTANTS = (
    "STATION", "ROOTS", "IO_FILE", "REF_EO_FILE", "OUT_EO", "OUT_REPORT", "B_VARIANT", "B2_VARIANT", "MODULES",
    "ARCHIVED_SHA256", "OUTCOME_RULES", "FORBIDDEN", "WINDOWS", "ANCHOR_DAY", "TARGET_HOURS", "MAX_DT_H", "GREY_RANGE",
    "THRESH_PX", "RANSAC_ITERS", "SEED", "MIN_INLIERS", "MIN_NEW_DAYS", "MIN_REF_PHOTOS", "MIN_EXTENT_PX",
    "MAX_RESID_MED_PX", "JACKKNIFE_MAX_PX", "AMBIGUITY_RATIO", "AMBIGUITY_SEP_DEG", "DAY_MIN_INLIERS", "DAY_TOL_PX",
    "MIN_DAYS_MEASURED", "REF_DAY_TOL_PX", "REF_HZ_TOL_DEG", "HZ_FIT_OK_PX", "HZ_MIN_CLEAR", "HZ_TILT_TOL",
    "HZ_ROLL_TOL", "HZ_WINDOW_PX", "HZ_EXPECTED_WINDOW_PX", "HZ_EDGE_MARGIN_PX", "CHAIN_TOL_PX", "SAME_SLOT_PAIRS",
    "PAN_TOL_DEG", "AUG_TOL_M", "ROT_SPREAD_FACTOR", "ROT_FLOOR_PX", "POSITION_TERM", "POSITION_SPREAD")


def constants():
    """{name: value} of PREREGISTERED_CONSTANTS as they are in effect now (JSON values: tuples as lists, the
    FORBIDDEN patterns as text)."""
    g = globals()
    return {n: json.loads(json.dumps([f.pattern for f in g[n]] if n == "FORBIDDEN" else g[n])) for n in
            PREREGISTERED_CONSTANTS}


def constants_sha256(c=None):
    return hashlib.sha256(json.dumps(constants() if c is None else c, sort_keys=True).encode()).hexdigest()


def loaded_modules():
    """{file name: sha256} of every repository module loaded in this process (the five MODULES are pinned; the
    others are recorded: none computes B)."""
    out = {}
    for m in list(sys.modules.values()):
        f = getattr(m, "__file__", None)
        if f and f.endswith(".py") and Path(f).resolve().parent == HERE:
            out[Path(f).name] = sha256(Path(f).resolve())
    return dict(sorted(out.items()))


PINS_FILE = HERE / "survey_products.py"     # its GCPFREE_PINS: the versions of this procedure whose B it accepts


def pinned_versions(path=None):
    """survey_products.py's GCPFREE_PINS, read with ast (never imported). -> (list or None, why not)."""
    import ast
    p = Path(path or PINS_FILE)
    try:
        tree = ast.parse(Path(guard(p)).read_text())
    except (OSError, SyntaxError, ValueError) as exc:
        return None, f"{p} cannot be read ({exc.__class__.__name__})"
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "GCPFREE_PINS" for t in node.targets):
            try:
                return ast.literal_eval(node.value), ""
            except ValueError:
                break
    return None, f"{p} has no GCPFREE_PINS"


def pin_problem():
    """'' when this gcpfree_eo.py and its modules are the version survey_products.py pins last (GCPFREE_PINS[-1]);
    else why not: a B made by any other version is 'NOT the pre-registered procedure' there, and the one delivery
    run would be spent on it, so nothing is pre-registered or run."""
    pins, why = pinned_versions()
    if not pins:
        return f"the versions of this procedure survey_products.py accepts cannot be read: {why or 'none pinned'}"
    own, mods, last = own_sha256(), module_sha256(), pins[-1]
    if last.get("gcpfree_eo.py") == own and last.get("modules") == mods:
        return ""
    bad = [m for m in MODULES if (last.get("modules") or {}).get(m) != mods[m]]
    return (f"this gcpfree_eo.py (sha256 {own[:12]}...)" + (f" or its modules ({', '.join(bad)})" if bad else "")
            + f" is not the version survey_products.py pins last (GCPFREE_PINS[-1]: gcpfree_eo.py "
              f"{str(last.get('gcpfree_eo.py'))[:12]}...): a B it made would be 'NOT the pre-registered procedure'. "
              f"Bring the checkout to the committed version (git status; git pull) before anything else")


def constants_problem():
    """'' when the pre-registered constants in effect in this process are the pinned ones (GCPFREE_PINS[-1]
    'constants_sha256'); else why not. Not a refusal (tests change the photo folders in their process): recorded in
    run.json, and survey_products.py calls such a run's B 'NOT the pre-registered procedure'."""
    pins, _ = pinned_versions()
    want = ((pins or [{}])[-1] or {}).get("constants_sha256")
    now = constants_sha256()
    if want == now:
        return ""
    return (f"the pre-registered constants in effect (sha256 {now[:12]}...) are not the pinned ones "
            f"({str(want)[:12]}...): a B of this run would be 'NOT the pre-registered procedure'")


def archived(path):
    """The calibration file at path, refused unless its content is the archived one (ARCHIVED_SHA256)."""
    p = guard(path)
    want = ARCHIVED_SHA256.get(p.name)
    got = sha256(p) if p.is_file() else None
    if want is None or got != want:
        raise Refused(f"{p}: " + ("not found" if got is None else
                                  "not one of the calibration files this procedure may read" if want is None else
                                  f"its content (sha256 {got[:12]}...) is not the archived {p.name} "
                                  f"({want[:12]}...)") + ": refused before anything is solved")
    return p


# ------------------------------------------------------------------ pre-registration

def read_preregistration(path):
    """{key: value} of PREREGISTERED.txt ('key: value' lines; '#' lines are comments; roots and
    modules_sha256 are JSON; every 'photo: <cam> <window> <day> <name> <size> <sha256>' line goes, in
    order, into 'photos' as a list)."""
    out = {"photos": []}
    for ln in Path(guard(path)).read_text().splitlines():
        s = ln.strip()
        if not s or s.startswith("#") or ":" not in s:
            continue
        k, v = s.split(":", 1)
        k, v = k.strip(), v.strip()
        if k == "photo":
            w = v.split()
            out["photos"].append(w[:4] + [int(w[4]) if w[4].isdigit() else w[4]] + w[5:] if len(w) == 6 else w)
            continue
        if k in ("roots", "modules_sha256"):
            try:
                v = json.loads(v)
            except ValueError:
                v = None
        out[k] = v
    return out


# ------------------------------------------------------------------ the pre-registered photos

def photo_manifest(roots, cams=("c1", "c2")):
    """The photos a run selects from these folders (select_photos: Oct and Jan for every camera, Dec for c2), each
    with its size and sha256: [[cam, window, day, name, size, sha256]], in the order a run reads them. Only names,
    sizes and bytes: nothing is computed from a photo here. The pre-registration lists them; a run on the station's
    folders that selects other photos is refused (a photo removed, replaced or added after a look)."""
    rows = []
    for cam in cams:
        io = load_intrinsics(archived(CAL / IO_FILE.format(cam=cam)))
        for k in ("oct", "jan") + (("dec",) if cam == "c2" else ()):
            days, _ = select_photos(roots, cam, WINDOWS[k], io)
            for day, paths in sorted(days.items()):
                for q in paths:
                    rows.append([cam, k, day, Path(q).name, os.path.getsize(str(guard(q))), sha256(q)])
    return rows


def photos_sha256(rows):
    return hashlib.sha256(json.dumps([list(r) for r in rows]).encode()).hexdigest()


def manifest_difference(want, now):
    """'' when the photos selected now are the pre-registered ones, else what differs (a few examples)."""
    w = {(r[0], r[1], r[3]): tuple(r[4:6]) for r in want}
    n = {(r[0], r[1], r[3]): tuple(r[4:6]) for r in now}
    gone = sorted(set(w) - set(n))
    new = sorted(set(n) - set(w))
    changed = sorted(k for k in set(w) & set(n) if w[k] != n[k])
    if not (gone or new or changed) and [list(r) for r in want] == [list(r) for r in now]:
        return ""
    parts = []
    for what, ks in (("no longer selected (removed, or replaced by another)", gone),
                     ("selected now but not pre-registered", new), ("of other content (size or sha256)", changed)):
        if ks:
            parts.append(f"{len(ks)} {what}: " + ", ".join(f"{c} {k} {nm}" for c, k, nm in ks[:3])
                         + (" ..." if len(ks) > 3 else ""))
    return "; ".join(parts) or "the same photos in another order"


# ------------------------------------------------------------------ git: records outside this branch's history

RECORD_NAME = re.compile(r"(^|/)(PREREGISTERED\.txt|runs/[^/]+/run\.json)$")


def foreign_records(base):
    """Commits OUTSIDE this branch's history (another branch, a stash, the reflog: a branch made from before the
    pre-registration, a reset, a rebase) that wrote a PREREGISTERED.txt or a runs/<id>/run.json under `base` whose
    content this branch's history does not hold (a rebased copy of a commit holds the same content and does not
    count). An earlier pre-registration or run there would otherwise be invisible from this branch.
    -> ([(commit, path)], why git could not tell, or '')."""
    anc = Path(base)
    while not anc.exists() and anc != anc.parent:
        anc = anc.parent
    spec = os.path.relpath(str(Path(base)), str(anc))
    env = dict(os.environ, GIT_TERMINAL_PROMPT="0")

    def blobs(*refs):
        r = subprocess.run(["git", "-C", str(anc), "log"] + list(refs) + ["--no-abbrev", "--format=@%H", "--raw",
                                                                           "--", spec],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, env=env)
        out, commit = [], None
        for ln in r.stdout.splitlines() if r.returncode == 0 else []:
            if ln.startswith("@"):
                commit = ln[1:].strip()
            elif ln.startswith(":") and "\t" in ln:
                meta, path = ln.split("\t", 1)
                f = meta.split()
                if len(f) >= 4 and RECORD_NAME.search(path.strip()) and set(f[3]) != {"0"}:
                    out.append((commit, path.strip(), f[3]))
        return out, r.returncode
    try:
        head, rc_h = blobs("HEAD")
        every, rc_a = blobs("--all", "--reflog")
    except OSError as exc:
        return [], f"git cannot be run ({exc})"
    if rc_a != 0:
        return [], "" if rc_h != 0 else "git log --all --reflog failed"
    have = {b for _, _, b in head}
    seen, found = set(), []
    for c, path, b in every:
        if b not in have and (c, path) not in seen:
            seen.add((c, path))
            found.append((c, path))
    return found, ""


def foreign_problem(base):
    """A run is refused while foreign_records() finds one. -> text or ''."""
    found, _ = foreign_records(base)
    if not found:
        return ""
    return (f"{len(found)} commit(s) outside this branch's history wrote a pre-registration or a run record that this "
            f"branch does not hold (" + ", ".join(f"{c[:10]} {p}" for c, p in found[:4])
            + (" ..." if len(found) > 4 else "") + "): another branch, a stash, or the reflog (a branch made from "
            "before the pre-registration, a reset) holds an earlier pre-registration or run, which would be invisible "
            "from here. Bring that history into this branch (git merge), never around it")


def git_state(path):
    """(committed and unchanged since its last commit?, the commit that last changed it, why not, its history:
    [{"commit", "date"}] of every commit that touched it, newest first)."""
    p = Path(path).resolve()
    env = dict(os.environ, GIT_TERMINAL_PROMPT="0")

    def git(*a):
        return subprocess.run(["git", "-C", str(p.parent)] + list(a), stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, universal_newlines=True, env=env)
    try:
        tracked = git("ls-files", "--error-unmatch", "--", p.name)
    except OSError as exc:
        return False, None, f"git cannot be run ({exc})", []
    if tracked.returncode != 0:
        return False, None, "it is not committed to git (git add + git commit + git push it first)", []
    log = git("log", "--format=%H %cI", "--", p.name)
    hist = [dict(zip(("commit", "date"), ln.split(" ", 1))) for ln in log.stdout.splitlines() if ln.strip()]
    st = git("status", "--porcelain", "--", p.name)
    if st.returncode != 0 or st.stdout.strip():
        return False, None, "it was changed after its last commit (commit the change, or restore it)", hist
    return True, (hist[0]["commit"] if hist else None), "", hist


def remote_branches(path, commit):
    """The remote-tracking branches (as this computer last pushed to or fetched them) whose history holds `commit`.
    -> ([branch names], why none: '' when there is one)."""
    p = Path(path).resolve()
    env = dict(os.environ, GIT_TERMINAL_PROMPT="0")
    try:
        r = subprocess.run(["git", "-C", str(p.parent), "branch", "-r", "--contains", commit], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, universal_newlines=True, env=env)
    except OSError as exc:
        return [], f"not checked: git cannot be run ({exc})"
    if r.returncode != 0:
        return [], f"not checked: git branch -r --contains failed ({r.stderr.strip()[:200]})"
    br = [ln.strip() for ln in r.stdout.splitlines() if ln.strip() and "->" not in ln]
    return br, ("" if br else "on no remote-tracking branch (not pushed)")


def preregistration(path, roots=None):
    """The pre-registration a run rests on, checked. -> (info, [problems]); problems empty when it may run."""
    p = Path(path)
    info = {"path": str(p.resolve()) if p.exists() else str(p)}
    if not p.is_file():
        return info, [f"{p} does not exist: the owner writes it with --preregister (choosing the outcome rule) "
                      f"and commits it before any run"]
    probs = []
    pr = read_preregistration(p)
    info.update(sha256=sha256(p), outcome_if_one_camera=pr.get("outcome_if_one_camera"),
                decided_by=pr.get("decided_by"), decided_on=pr.get("decided_on"), roots=pr.get("roots"))
    ok, commit, why, hist = git_state(p)
    info["commit"], info["git_history"] = commit, hist
    if not ok:
        probs.append(f"{p}: {why}")
    elif commit:
        br, why_r = remote_branches(p, commit)
        info["remote_branches"] = br
        if not br:
            probs.append(f"{p}: its commit {commit[:10]} is {why_r}: push it (git push) BEFORE any run, a trial too. The "
                         f"time GitHub receives the push is the rule's only public timestamp (a commit's date is set "
                         f"by this computer's clock), and a rule never pushed could be chosen again after a look (a "
                         f"local reset). After the push, git branch -r --contains {commit[:10]} names the remote branch")
    if pr.get("procedure_sha256") != own_sha256():
        probs.append(f"{p} pre-registers gcpfree_eo.py sha256 {str(pr.get('procedure_sha256'))[:12]}..., this "
                     f"gcpfree_eo.py is {own_sha256()[:12]}...: not the pre-registered procedure")
    mods, now = pr.get("modules_sha256") or {}, module_sha256()
    bad = [m for m in MODULES if mods.get(m) != now[m]]
    if bad:
        probs.append(f"the procedure's modules differ from the pre-registered ones: {', '.join(bad)}")
    pp = pin_problem()
    if pp:
        probs.append(pp)
    if pr.get("outcome_if_one_camera") not in OUTCOME_RULES:
        probs.append(f"outcome_if_one_camera is {pr.get('outcome_if_one_camera')!r}: it must be one of "
                     f"{' | '.join(OUTCOME_RULES)} (the owner's choice)")
    if not pr.get("decided_by"):
        probs.append("decided_by is empty: who chose the outcome rule")
    if not isinstance(pr.get("roots"), list) or [str(r) for r in pr["roots"]] != [str(r) for r in ROOTS]:
        probs.append(f"roots must be the station's three photo folders ({', '.join(ROOTS)}), as a delivery run "
                     f"uses them; it says {pr.get('roots')}")
    rows = pr.get("photos") or []
    info.update(photos=rows, photos_sha256=pr.get("photos_sha256"))
    if not rows:
        probs.append(f"{p} pins no photos (no 'photo:' lines): write it with --preregister of this version")
    elif any(len(r) != 6 for r in rows) or photos_sha256(rows) != pr.get("photos_sha256"):
        probs.append(f"{p}: its 'photo:' lines are not those its photos_sha256 line was written for (edited)")
    elif {r[0] for r in rows} != {"c1", "c2"}:
        probs.append(f"{p} pins photos of {', '.join(sorted({r[0] for r in rows}))} only: a delivery run solves c1 "
                     f"and c2")
    return info, probs


def write_preregistration(path, outcome, decided_by, roots, photos=None):
    """Writes the pre-registration once: the owner's rule, who decided, the procedure's and its modules' sha256, the
    station's photo folders and the photos a run selects from them (photo_manifest; `photos` given only by tests
    that have no such folders). Refused (nothing written) unless every folder is there and every camera has Oct and
    Jan photos (check_photos), as a delivery run needs them."""
    p = Path(path)
    if p.exists():
        raise Refused(f"{p} exists: a pre-registration is written once (before any run it may still be edited and "
                      f"committed again; once a run was made under it, never)")
    if [str(r) for r in roots] != [str(r) for r in ROOTS]:
        raise Refused(f"--roots: a pre-registration names the station's three photo folders ({', '.join(ROOTS)}); "
                      f"other folders make TRIAL runs only")
    if outcome not in OUTCOME_RULES:
        raise Refused(f"--outcome-if-one-camera must be one of {' | '.join(OUTCOME_RULES)} (the owner's choice)")
    if not (decided_by or "").strip():
        raise Refused("--decided-by: who chose the outcome rule")
    pp = pin_problem()
    if pp:
        raise Refused(f"nothing pre-registered: {pp}")
    fp = foreign_problem(p.parent)
    if fp:
        raise Refused(f"nothing pre-registered: {fp}")
    if photos is None:
        check_photos([str(r) for r in roots], ["c1", "c2"], None, True)
        photos = photo_manifest([str(r) for r in roots])
    p.parent.mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc)
    p.write_text(
        "# PRE-REGISTRATION of calibration B (gcpfree_eo.py) for the BLIND 23 Oct 2024 DEM\n"
        f"# Written by gcpfree_eo.py --preregister on {now:%Y-%m-%d %H:%M} UTC. COMMIT AND PUSH THIS FILE before any\n"
        f"# run:  git add {p.name} && git commit -m 'gcpfree: pre-registration' && git push\n"
        "# gcpfree_eo.py refuses to run while it is missing, uncommitted, changed since its commit or not pushed\n"
        "# (its commit on no remote-tracking branch), or while the procedure or its modules differ from the\n"
        "# sha256s below.\n"
        "procedure: gcpfree_eo.py\n"
        f"procedure_sha256: {own_sha256()}\n"
        f"modules_sha256: {json.dumps(module_sha256(), sort_keys=True)}\n"
        f"roots: {json.dumps([str(r) for r in roots])}\n"
        f"outcome_if_one_camera: {outcome}\n"
        f"decided_by: {decided_by.strip()}\n"
        f"decided_on: {now:%Y-%m-%d}\n"
        "# outcome_if_one_camera: what B is when the procedure's gates pass for one camera only:\n"
        "#   both-or-none  B is delivered only when both cameras are; otherwise no B (Chris gets A alone, with\n"
        "#                 the reports; the passing camera's result is WITHHELD in its report)\n"
        "#   per-camera    each passing camera is delivered (a one-camera B DEM; its README names the other\n"
        "#                 camera NOT DELIVERED, with the procedure's report)\n"
        "# Fixed by the procedure, not a choice: a delivery run uses the roots above and both cameras; the first\n"
        "# finished delivery run fixes B; one more delivery run is allowed only with --aug-dir (Aug 2024 photos)\n"
        "# when the first used none and left a camera NOT DELIVERED, and it solves only those cameras. Its B is\n"
        "# the delivered B only when the first run delivered no camera; otherwise it is a second product\n"
        "# (--calibration-variant gcpfree-aug) and the first run's B stays the delivered one.\n"
        "# Once any run is made under this file, it can never change (gcpfree_eo.py refuses).\n"
        "# THE PHOTOS: every photo a run selects from the folders above (Oct and Jan for each camera, Dec for c2),\n"
        "# with its size and sha256. A run on these folders that selects other photos is refused.\n"
        f"photos_sha256: {photos_sha256(photos)}\n"
        + "".join(f"photo: {' '.join(str(x) for x in r)}\n" for r in photos))
    return p


# ------------------------------------------------------------------ photos

def epoch_of(p):
    m = re.match(r"^(\d{9,11})\.", Path(p).name)
    return int(m.group(1)) if m else None


def select_photos(roots, cam, window, io):
    """{day: [paths]} -- up to len(TARGET_HOURS) photos per day, the nearest to each target hour. The scan reads
    each candidate's image header for its size: a photo-named file whose name, or the name of the file a link points
    to, is forbidden (guard) is never opened for it, and counts as unreadable (counts['forbidden_not_opened'])."""
    import detect_original_view as DOV
    real_size, skipped = DOV.jpeg_size, []

    def size_unless_forbidden(path):
        if forbidden_name(path):
            skipped.append(Path(path).name)
            return None
        return real_size(path)
    DOV.jpeg_size = size_unless_forbidden
    try:
        found, counts = DOV.collect_photos(roots, cam, window[0], window[1], (0.0, 24.0), station=STATION,
                                           frame_size=(int(io[0]), int(io[1])))
    finally:
        DOV.jpeg_size = real_size
    counts["forbidden_not_opened"] = skipped
    by_day = {}
    for p in found:
        t = datetime.fromtimestamp(epoch_of(p), tz=timezone.utc)
        by_day.setdefault(t.strftime("%Y-%m-%d"), []).append((t.hour + t.minute / 60.0, p))
    out = {}
    for day in sorted(by_day):
        chosen = []
        for h in TARGET_HOURS:
            cands = sorted((abs(hh - h), str(p), p) for hh, p in by_day[day]
                           if abs(hh - h) <= MAX_DT_H and p not in chosen)
            if cands:
                chosen.append(cands[0][2])
        if chosen:
            out[day] = chosen
    return out, counts


class Photo:
    def __init__(self, path, day, mask):
        import cv2
        self.path, self.day = Path(path), day
        t = datetime.fromtimestamp(epoch_of(path), tz=timezone.utc)
        h = t.hour + t.minute / 60.0
        self.slot = int(np.argmin([abs(h - x) for x in TARGET_HOURS]))
        img = cv2.imread(str(guard(path)), cv2.IMREAD_GRAYSCALE)
        self.grey = float(img.mean()) if img is not None else -1.0
        self.ok = img is not None and GREY_RANGE[0] <= self.grey <= GREY_RANGE[1]
        self.pts, self.des = (EER.features(path, mask)[1:] if self.ok else (np.empty((0, 2)), None))


def load_set(days, mask):
    photos = []
    for day, paths in sorted(days.items()):
        for p in paths:
            ph = Photo(p, day, mask)
            if ph.ok:
                photos.append(ph)
    return photos


# ------------------------------------------------------------------ solving

def angle_deg(Q):
    return float(np.degrees(np.arccos(np.clip((np.trace(Q) - 1) / 2, -1, 1))))


def dist_deg(Q1, Q2):
    return angle_deg(Q1 @ Q2.T)


def disp_px(Q1, Q2, io, step=64):
    """Largest image displacement (px) between two pointings of the same camera over the lower two
    thirds of the frame (the beach): roll counts by its effect at the frame edge, pan and tilt fully."""
    U, V = np.meshgrid(np.arange(step / 2, io[0], step), np.arange(io[1] / 3, io[1], step))
    r = EER.rays(U.ravel(), V.ravel(), io)
    a = np.arccos(np.clip(np.sum((r @ (Q1 @ Q2.T).T) * r, axis=1), -1, 1))
    return float(a.max() * io[4])


def ray_err(Q, ra, rb):
    return np.arccos(np.clip(np.sum((ra @ Q.T) * rb, axis=1), -1, 1))


def ransac2(ra, rb, thresh, iters=RANSAC_ITERS, seed=SEED):
    """Rotation by RANSAC over 2-ray samples, then refits on the inliers (as estimate_eo_rotation.py)."""
    n = len(ra)
    if n < 3:
        return None, np.zeros(n, bool)
    rng = np.random.default_rng(seed)
    best = np.zeros(n, bool)
    for _ in range(iters):
        i, j = rng.choice(n, 2, replace=False)
        if abs(float(ra[i] @ ra[j])) > 0.99996:          # nearly the same ray: no rotation from it
            continue
        Q = EER.kabsch(ra[[i, j]], rb[[i, j]])
        inl = ray_err(Q, ra, rb) < thresh
        if inl.sum() > best.sum():
            best = inl
    if best.sum() < 3:
        return None, best
    Q = EER.kabsch(ra[best], rb[best])
    for _ in range(3):
        best = ray_err(Q, ra, rb) < thresh
        Q = EER.kabsch(ra[best], rb[best])
    return Q, best


class Matches:
    """Pooled matches between a set of new photos and a set of reference photos."""

    def __init__(self, new, ref, io):
        A, B, ni, ri = [], [], [], []
        for a, pn in enumerate(new):
            for b, pr in enumerate(ref):
                if SAME_SLOT_PAIRS and pn.slot != pr.slot:
                    continue
                for i, j in EER.match(pr.des, pn.des):
                    A.append(pr.pts[i]); B.append(pn.pts[j]); ni.append(a); ri.append(b)
        self.new, self.ref = new, ref
        self.A, self.B = np.array(A).reshape(-1, 2), np.array(B).reshape(-1, 2)
        self.ni, self.ri = np.array(ni, int), np.array(ri, int)
        self.ra = EER.rays(self.A[:, 0], self.A[:, 1], io) if len(A) else np.empty((0, 3))
        self.rb = EER.rays(self.B[:, 0], self.B[:, 1], io) if len(A) else np.empty((0, 3))


def solve(M, io, label):
    """Pooled solve with all the gate-F metrics."""
    thresh = THRESH_PX / io[4]
    out = {"label": label, "n_new_photos": len(M.new), "n_ref_photos": len(M.ref), "matches": int(len(M.A))}
    Q, inl = ransac2(M.ra, M.rb, thresh)
    if Q is None:
        out.update(inliers=int(inl.sum()), gate_F=False, why="no consistent rotation")
        return None, out
    err = ray_err(Q, M.ra[inl], M.rb[inl]) * io[4]
    new_days = sorted({M.new[k].day for k in M.ni[inl]})
    ref_photos = sorted({M.ref[k].path.name for k in M.ri[inl]})
    ext_u = float(np.percentile(M.B[inl, 0], 90) - np.percentile(M.B[inl, 0], 10))
    ext_v = float(np.percentile(M.B[inl, 1], 90) - np.percentile(M.B[inl, 1], 10))
    jack = []
    for d in new_days:                                    # leave one new day out (no new RANSAC)
        keep = inl & np.array([M.new[k].day != d for k in M.ni])
        if keep.sum() >= 3:
            jack.append(disp_px(EER.kabsch(M.ra[keep], M.rb[keep]), Q, io))
    rest = ~inl                                           # a second consensus? (dirt, overlays, repeats)
    Q2, inl2 = ransac2(M.ra[rest], M.rb[rest], thresh, iters=RANSAC_ITERS // 2, seed=SEED + 1)
    second = int(inl2.sum()) if Q2 is not None else 0
    ambiguous = Q2 is not None and second >= AMBIGUITY_RATIO * inl.sum() and dist_deg(Q2, Q) > AMBIGUITY_SEP_DEG
    out.update(inliers=int(inl.sum()), inlier_pct=round(100.0 * inl.mean(), 1),
               resid_med_px=round(float(np.median(err)), 3), resid_p90_px=round(float(np.percentile(err, 90)), 3),
               rotation_deg=round(angle_deg(Q), 4), inlier_new_days=new_days, inlier_ref_photos=len(ref_photos),
               extent_px=[round(ext_u), round(ext_v)], jackknife_max_px=round(max(jack), 2) if jack else None,
               second_consensus=second, ambiguous=bool(ambiguous))
    why = []
    if inl.sum() < MIN_INLIERS:
        why.append(f"{int(inl.sum())} inliers < {MIN_INLIERS}")
    if len(new_days) < MIN_NEW_DAYS:
        why.append(f"inliers on {len(new_days)} new-photo days < {MIN_NEW_DAYS}")
    if len(ref_photos) < MIN_REF_PHOTOS:
        why.append(f"inliers on {len(ref_photos)} reference photos < {MIN_REF_PHOTOS}")
    if min(ext_u, ext_v) < MIN_EXTENT_PX:
        why.append(f"inliers span only {ext_u:.0f} x {ext_v:.0f} px (< {MIN_EXTENT_PX:.0f})")
    if np.median(err) > MAX_RESID_MED_PX:
        why.append(f"residual {np.median(err):.2f} px > {MAX_RESID_MED_PX}")
    if not jack or max(jack) > JACKKNIFE_MAX_PX:
        why.append(f"leave-one-day-out {max(jack) if jack else float('nan'):.1f} px > {JACKKNIFE_MAX_PX}")
    if ambiguous:
        why.append(f"a second consensus of {second} inliers")
    out["gate_F"] = not why
    out["why"] = "; ".join(why)
    return Q, out


def per_day(M, io, Q):
    """Each new day alone (its own RANSAC). Gate S."""
    thresh = THRESH_PX / io[4]
    rows = []
    for d in sorted({p.day for p in M.new}):
        sel = np.array([M.new[k].day == d for k in M.ni], bool)
        if sel.sum() < 3:
            rows.append({"day": d, "inliers": 0})
            continue
        Qd, inl = ransac2(M.ra[sel], M.rb[sel], thresh, iters=RANSAC_ITERS // 2, seed=SEED + 2)
        r = {"day": d, "inliers": int(inl.sum())}
        if Qd is not None and inl.sum() >= DAY_MIN_INLIERS:
            r["diff_px"] = round(disp_px(Qd, Q, io), 2)
        rows.append(r)
    measured = [r for r in rows if "diff_px" in r]
    bad = [r["day"] for r in measured if r["diff_px"] > DAY_TOL_PX]
    ok = len(measured) >= MIN_DAYS_MEASURED and not bad
    why = ("" if ok else f"{len(measured)} days measured (need {MIN_DAYS_MEASURED})"
           + (f"; days off the pooled pointing: {bad}" if bad else ""))
    return {"days": rows, "measured": len(measured), "off_days": bad, "gate_S": ok, "why": why}


def eo_from(Q, eo_ref):
    az, tilt, roll = EER.R_to_angles(Q @ cirn_angles_to_R(*eo_ref[3:6]))
    return np.r_[eo_ref[:3], az, tilt, roll]


def angles_deg(eo):
    a = np.degrees(np.asarray(eo, float)[3:6])
    a[0] = a[0] % 360.0
    return a


def eo_deg(eo):
    """[x, y, z, azimuth, tilt, roll] (m, deg), unrounded, as the report keeps B."""
    return [float(v) for v in np.r_[np.asarray(eo, float)[:3], angles_deg(eo)]]


def eo_rad(v):
    v = np.asarray(v, float)
    return np.r_[v[:3], np.radians(v[3:6])]


# ------------------------------------------------------------------ horizon

_WARP = {}


def warp_maps(io, Q):
    """cv2.remap maps of the frame a camera turned by Q (rays_new = Q rays_ref, same centre) sees, built from
    a reference frame: exact for a pure rotation. Pixels outside the old view are marked (-1e4)."""
    key = (tuple(np.round(np.asarray(io, float), 9)), tuple(np.round(np.asarray(Q, float).ravel(), 12)))
    if key not in _WARP:
        from view_reproject import distort_uv
        nu, nv = int(io[0]), int(io[1])
        U, V = np.meshgrid(np.arange(nu, dtype=float), np.arange(nv, dtype=float))
        x, y = undistort_uv(U.ravel(), V.ravel(), io)
        K = EER.K_of(io)
        r_ref = Q.T @ np.linalg.solve(K, np.vstack([x, y, np.ones_like(x)]))
        p = K @ r_ref
        with np.errstate(invalid="ignore", divide="ignore"):
            us, vs = p[0] / p[2], p[1] / p[2]
        bad = (p[2] <= 0) | ~np.isfinite(us) | ~np.isfinite(vs)
        us, vs = np.where(bad, io[2], us), np.where(bad, io[3], vs)
        bad |= np.hypot((us - io[2]) / io[4], (vs - io[3]) / io[5]) > 1.2   # beyond the lens model's reach
        ud, vd = distort_uv(us, vs, io)
        ud, vd = np.where(bad, -1e4, ud), np.where(bad, -1e4, vd)
        _WARP.clear()                                     # one link at a time: keep memory down
        _WARP[key] = (ud.reshape(nv, nu).astype(np.float32), vd.reshape(nv, nu).astype(np.float32))
    return _WARP[key]


def warp_image(grey, maps):
    """The frame a camera turned by Q would see, from a reference frame. Where the reference frame saw nothing,
    each column repeats its nearest reference pixel above or below: no new content and no horizontal edge (the
    horizon is read from vertical brightness changes) where the old view ends."""
    import cv2
    mx, my = maps
    out = cv2.remap(grey, mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT, borderValue=0)
    inside = cv2.remap(np.full(grey.shape[:2], 255, np.uint8), mx, my, cv2.INTER_NEAREST,
                       borderMode=cv2.BORDER_CONSTANT, borderValue=0) > 0
    h = out.shape[0]
    rows = np.arange(h)[:, None]
    above = np.maximum.accumulate(np.where(inside, rows, -1), axis=0)              # last inside row at or above
    below = np.minimum.accumulate(np.where(inside, rows, h)[::-1], axis=0)[::-1]  # first inside row at or below
    src = np.clip(np.where(above >= 0, above, below), 0, h - 1)
    return out[src, np.arange(out.shape[1])[None, :]]


HZ_WINDOW_PX = 250                  # observed_horizon's search window around the predicted row
# The turned reference photos (the expected bias): their horizon lies where the reference EO, corrected by its own
# measured offset and turned by Q, puts it (exact for a rotation, to the reference fit's residual): looked for
# within this many px of it, in columns where the turned photos have reference content over that window and a
# margin for the edge filter (an implementation of the 14 Oct 2026 band, not a gate: no threshold changed).
HZ_EXPECTED_WINDOW_PX = 60
HZ_EDGE_MARGIN_PX = 25


def horizon_cols(io):
    return np.arange(150, int(io[0]) - 150, 100)


def corrected_eo(eo, h):
    """eo with the tilt and roll change that puts its horizon on the photos h was measured on."""
    e = np.asarray(eo, float).copy()
    e[4] += np.radians(h["dtilt_med"])
    e[5] += np.radians(h["droll_med"])
    return e


def covered_cols(io, Q, near):
    """The horizon columns (and their rows `near`) where the frame turned by Q has reference content within
    HZ_EDGE_MARGIN_PX of the row its horizon should be on: the part of a turned frame that the reference frame
    never saw holds no horizon (warp_image repeats the nearest seen pixel there)."""
    cols = horizon_cols(io)
    mx, my = warp_maps(io, Q)
    nu, nv = int(io[0]), int(io[1])
    keep, rows_kept = [], []
    for c, r in zip(cols, near):
        if not np.isfinite(r) or r < 0 or r > nv - 1:
            continue
        rows = np.clip(np.round(r + np.linspace(-HZ_EDGE_MARGIN_PX, HZ_EDGE_MARGIN_PX, 11)).astype(int), 0, nv - 1)
        u, v = mx[rows, int(c)], my[rows, int(c)]
        if ((u >= 0) & (u <= nu - 1) & (v >= 0) & (v <= nv - 1)).all():
            keep.append(c)
            rows_kept.append(r)
    return np.array(keep, dtype=cols.dtype), np.array(rows_kept, float)


def horizon_fits(io, eo, photos, warp=None, cols=None, near=None, window=HZ_WINDOW_PX):
    """(tilt, roll) changes that move eo's horizon onto each photo's (horizon_check.fit_tilt_roll); with warp=Q,
    onto each photo turned by Q (the frame a camera turned by Q would see); over the given columns (default
    every 100 px), the horizon looked for within `window` px of eo's rows (or of `near`)."""
    import cv2
    from horizon_check import fit_tilt_roll
    cols = horizon_cols(io) if cols is None else np.asarray(cols)
    pred = EER.horizon_rows(io, eo, cols) if near is None else np.asarray(near, float)
    maps = warp_maps(io, warp) if warp is not None else None
    rows = []
    for ph in photos:
        g = cv2.imread(str(guard(ph.path)), cv2.IMREAD_GRAYSCALE)
        if g is None:
            rows.append({"photo": ph.path.name, "clear": False})
            continue
        if maps is not None:
            g = warp_image(g, maps)
        obs = EER.observed_horizon(g, cols, pred, window=window)
        if np.isfinite(obs).sum() < 6:
            rows.append({"photo": ph.path.name, "clear": False})
            continue
        (dt, dr), res = fit_tilt_roll(io, eo, cols, obs)
        rows.append({"photo": ph.path.name, "dtilt": round(float(dt), 4), "droll": round(float(dr), 4),
                     "fit_px": round(float(np.median(res)), 2), "clear": bool(np.median(res) <= HZ_FIT_OK_PX)})
    clear = [r for r in rows if r["clear"]]
    med = (float(np.median([r["dtilt"] for r in clear])), float(np.median([r["droll"] for r in clear]))) \
        if clear else (None, None)
    return {"photos": rows, "n_clear": len(clear), "dtilt_med": med[0], "droll_med": med[1]}


def horizon_prediction(eo_ref, href, Q, eo_b):
    """The bias a perfect B would show, from geometry alone: the reference EO's horizon offset taken as a
    rotation of the reference camera, turned by Q (reported only; the band uses the measured bias)."""
    if href.get("dtilt_med") is None:
        return None
    fix = np.asarray(eo_ref, float).copy()
    fix[4] += np.radians(href["dtilt_med"])
    fix[5] += np.radians(href["droll_med"])
    a_true = np.degrees(EER.R_to_angles(Q @ cirn_angles_to_R(*fix[3:6])))
    d = a_true - np.degrees(np.asarray(eo_b, float)[3:6])
    d[0] = (d[0] + 180.0) % 360.0 - 180.0
    return {"tilt_deg": round(float(d[1] - href["dtilt_med"]), 4), "roll_deg": round(float(d[2] - href["droll_med"]), 4),
            "pan_carried_deg": round(float(d[0]), 4)}


def band(bias, tol):
    return [round(min(0.0, bias) - tol, 4), round(max(0.0, bias) + tol, 4)]


def horizon_check(io, eo_ref, ref_photos, eo_b, new_photos, Q, href=None):
    """The differential horizon check of B, judged in the band [0, expected bias] widened by the noise
    tolerance (see the docstring). -> dict with status 'pass' | 'FAIL' | 'unavailable'."""
    href = href or horizon_fits(io, eo_ref, ref_photos)
    if href["n_clear"] < HZ_MIN_CLEAR:
        return {"reference": href, "status": "unavailable",
                "why": f"clear horizons on the reference photos: {href['n_clear']} (need {HZ_MIN_CLEAR})"}
    # where the turned reference photos show their horizon: the reference EO as its photos see it, turned by Q
    near = EER.horizon_rows(io, eo_from(Q, corrected_eo(eo_ref, href)), horizon_cols(io))
    cols, near = covered_cols(io, Q, near)             # the same columns for both B fits
    hexp = horizon_fits(io, eo_b, ref_photos, warp=Q, cols=cols, near=near, window=HZ_EXPECTED_WINDOW_PX)
    hnew = horizon_fits(io, eo_b, new_photos, cols=cols)
    out = {"reference": href, "new": hnew, "expected": hexp, "prediction": horizon_prediction(eo_ref, href, Q, eo_b),
           "columns": {"used": int(len(cols)), "of": int(len(horizon_cols(io)))}}
    n = (href["n_clear"], hnew["n_clear"], hexp["n_clear"])
    if min(n) < HZ_MIN_CLEAR:
        out.update(status="unavailable", why=f"clear horizons: reference {n[0]}, new {n[1]}, reference turned by "
                                             f"the solved rotation {n[2]} (need {HZ_MIN_CLEAR} each)")
        return out
    dt, dr = hnew["dtilt_med"] - href["dtilt_med"], hnew["droll_med"] - href["droll_med"]
    et, er = hexp["dtilt_med"] - href["dtilt_med"], hexp["droll_med"] - href["droll_med"]
    bt, br = band(et, HZ_TILT_TOL), band(er, HZ_ROLL_TOL)
    ok = bt[0] <= dt <= bt[1] and br[0] <= dr <= br[1]
    out.update(diff_tilt_deg=round(dt, 4), diff_roll_deg=round(dr, 4), bias_tilt_deg=round(et, 4),
               bias_roll_deg=round(er, 4), band_tilt_deg=bt, band_roll_deg=br, status="pass" if ok else "FAIL",
               why="" if ok else f"horizon differs from the reference by tilt {dt:+.3f}, roll {dr:+.3f} deg, outside "
                                 f"the band tilt [{bt[0]:+.3f}, {bt[1]:+.3f}], roll [{br[0]:+.3f}, {br[1]:+.3f}] "
                                 f"(expected bias {et:+.3f} / {er:+.3f})")
    return out


def horizon_text(h):
    if not h:
        return "n/a"
    if "diff_tilt_deg" not in h:
        return f"{h.get('status')} ({h.get('why')})"
    return (f"{h['status']}: tilt {h['diff_tilt_deg']:+.3f} in [{h['band_tilt_deg'][0]:+.3f}, {h['band_tilt_deg'][1]:+.3f}], "
            f"roll {h['diff_roll_deg']:+.3f} in [{h['band_roll_deg'][0]:+.3f}, {h['band_roll_deg'][1]:+.3f}] deg "
            f"(expected bias {h['bias_tilt_deg']:+.3f} / {h['bias_roll_deg']:+.3f}"
            + (f"; geometric prediction {h['prediction']['tilt_deg']:+.3f} / {h['prediction']['roll_deg']:+.3f}"
               if h.get("prediction") else "") + ")")


def ground_diff_m(io, eo1, eo2, z=0.5, step=48):
    """Median / P90 distance (m) between where two EOs put the same beach pixels on the plane z (50-400 m)."""
    U, V = np.meshgrid(np.arange(step / 2, io[0], step), np.arange(step / 2, io[1], step))
    X1, Y1 = pixel_to_ground(U.ravel(), V.ravel(), z, io, eo1, max_range_m=400.0)
    X2, Y2 = pixel_to_ground(U.ravel(), V.ravel(), z, io, eo2, max_range_m=400.0)
    r = np.hypot(X1 - eo1[0], Y1 - eo1[1])
    d = np.hypot(X1 - X2, Y1 - Y2)
    ok = np.isfinite(d) & np.isfinite(r)
    ok[ok] = r[ok] > 50
    return (float(np.median(d[ok])), float(np.percentile(d[ok], 90))) if ok.any() else (None, None)


# ------------------------------------------------------------------ one camera

def reference_set(cam, io, eo_ref, days, anchor, mask_ref, mask_new, log):
    """Reference photos (Jan, or Aug): the anchor (calibration) day and every other day of the window
    whose pointing is not MEASURED to differ from the anchor's (see the docstring)."""
    if not days:
        return [], {"why": "no reference photos"}
    anchor = anchor if anchor in days else max(days)
    anc = load_set({anchor: days[anchor]}, mask_ref)
    h_anc = horizon_fits(io, eo_ref, anc)
    kept, rows = list(anc), []
    for d in sorted(days):
        if d == anchor:
            continue
        ph_ref = load_set({d: days[d]}, mask_ref)
        ph_new = load_set({d: days[d]}, mask_new)
        M = Matches(ph_new, anc, io)
        Q, inl = ransac2(M.ra, M.rb, THRESH_PX / io[4], iters=RANSAC_ITERS // 2, seed=SEED + 3)
        r = {"day": d, "inliers": int(inl.sum())}
        if Q is not None and inl.sum() >= DAY_MIN_INLIERS:
            r["rotation_deg"] = round(angle_deg(Q), 4)
            r["diff_px"] = round(disp_px(Q, np.eye(3), io), 2)
            r["kept"] = r["diff_px"] <= REF_DAY_TOL_PX
            r["by"] = "features"
        else:
            h = horizon_fits(io, eo_ref, ph_ref)
            if h["n_clear"] and h_anc["n_clear"]:
                dt, dr = h["dtilt_med"] - h_anc["dtilt_med"], h["droll_med"] - h_anc["droll_med"]
                r.update(horizon_dtilt=round(dt, 4), horizon_droll=round(dr, 4),
                         kept=abs(dt) <= REF_HZ_TOL_DEG and abs(dr) <= REF_HZ_TOL_DEG, by="horizon (pan not measured)")
            else:
                r.update(kept=True, by="not measurable: kept (no evidence of a move)")
        if r["kept"]:
            kept += ph_ref
        rows.append(r)
        log(f"    reference day {d} vs {anchor}: {r['inliers']} inliers"
            + (f", {r['diff_px']:.1f} px" if "diff_px" in r else "")
            + (" kept" if r["kept"] else " LEFT OUT") + f" ({r['by']})")
    return kept, {"anchor": anchor, "anchor_photos": [p.path.name for p in anc], "days": rows}


def run_link(cam, io, eo_ref, ref_days, anchor, new_days, label, log):
    """reference set -> new set: the reference horizon (logged before any new photo is read), pooled solve (F),
    per-day (S), horizon band. -> (Q, B, info)"""
    shape = (int(io[1]), int(io[0]))
    mask_ref = EER.land_mask(shape, io, eo_ref)
    mask_new = EER.land_mask(shape, io, eo_ref, margin=-250)
    ref, refinfo = reference_set(cam, io, eo_ref, ref_days, anchor, mask_ref, mask_new, log)
    info = {"label": label, "reference_check": refinfo, "ref_photos": [p.path.name for p in ref]}
    href = horizon_fits(io, eo_ref, ref) if ref else None
    if href is not None:
        info["reference_horizon"] = {k: href[k] for k in ("n_clear", "dtilt_med", "droll_med")}
        log(f"  {label}: the reference EO against the visible horizon of its own {len(ref)} kept photos (measured "
            f"before any new photo is read): "
            + (f"tilt {href['dtilt_med']:+.3f}, roll {href['droll_med']:+.3f} deg ({href['n_clear']} clear)"
               if href["n_clear"] else "no clear horizon"))
    new = load_set(new_days, mask_new)
    info["new_photos"] = [p.path.name for p in new]
    if not new or not ref:
        info.update(gate_F=False, gate_S=False, why="no usable photos on one side")
        return None, None, info
    M = Matches(new, ref, io)
    Q, sol = solve(M, io, label)
    info["solve"] = sol
    if Q is None:
        info.update(gate_F=False, gate_S=False, why=sol["why"])
        return None, None, info
    info["per_day"] = per_day(M, io, Q)
    B = eo_from(Q, eo_ref)
    info["eo"] = [round(v, 6) for v in eo_deg(B)]
    info["gate_F"], info["gate_S"] = sol["gate_F"], info["per_day"]["gate_S"]
    info["horizon"] = horizon_check(io, eo_ref, ref, B, new, Q, href)
    log(f"  {label}: {sol['inliers']} inliers of {sol['matches']} matches, residual {sol['resid_med_px']} px, "
        f"turn {sol['rotation_deg']:.2f} deg, days {sol['inlier_new_days']}, gate F "
        f"{'pass' if sol['gate_F'] else 'FAIL: ' + sol['why']}")
    log("    per day: " + ", ".join(f"{r['day'][5:]} {r['inliers']}" + (f"/{r['diff_px']:.1f}px" if "diff_px" in r else "")
                                    for r in info["per_day"]["days"])
        + f" -> gate S {'pass' if info['gate_S'] else 'FAIL: ' + info['per_day']['why']}")
    log(f"    horizon (differential, band): {horizon_text(info['horizon'])}")
    return Q, B, info


def link_counts(link):
    """A route (one solve link) counts when it passes gates F and S and its horizon check is not FAIL."""
    return (bool(link) and bool(link.get("gate_F")) and bool(link.get("gate_S"))
            and (link.get("horizon") or {}).get("status") != "FAIL")


def link_hz(link):
    return ((link or {}).get("horizon") or {}).get("status")


def link_why(lk):
    """Why a link failed: its photos, gate F, gate S and the horizon, as far as it got."""
    lk = lk or {}
    parts = [lk.get("why", ""), (lk.get("solve") or {}).get("why", ""), (lk.get("per_day") or {}).get("why", "")]
    h = lk.get("horizon") or {}
    if h.get("status") == "FAIL":
        parts.append("horizon FAIL: " + h.get("why", ""))
    elif h.get("status") == "unavailable":
        parts.append("horizon unavailable: " + h.get("why", ""))
    out = []
    for x in parts:                                       # each reason once (a failed solve's 'why' is the link's too)
        if x and x not in out:
            out.append(x)
    return "; ".join(out) or "passed"


def decide(res):
    """The pre-registered verdict of one camera from its routes (see the docstring: VERDICT). Pure: reads only
    res (direct, chain, aug with has_B and the comparisons run_camera measured). -> (verdict, source, reasons,
    notes)."""
    direct, chain, aug = res.get("direct") or {}, res.get("chain"), res.get("aug")
    reasons, notes = [], []
    d_ok = bool(direct.get("has_B")) and link_counts(direct)
    chain_ok = (chain is not None and bool(chain.get("has_B")) and link_counts(chain.get("link1"))
                and link_counts(chain.get("link2")))
    aug_ok = aug is not None and bool(aug.get("has_B")) and link_counts(aug.get("link"))
    if chain is not None and not chain_ok:
        notes.append("chain Jan->Dec->Oct does not count: " + "; ".join(
            f"{k}: {link_why(chain.get(k))}" for k in ("link1", "link2") if k in chain)
            + ("" if "link2" in chain else "; link 2 not run (link 1 failed)"))
    if aug is not None and not aug_ok:
        notes.append("Aug->Oct does not count: " + link_why(aug.get("link")))
    if d_ok:
        verdict, source = "DELIVERED", "direct Jan->Oct"
        if link_hz(direct) != "pass":
            verdict = "DELIVERED-WIDER"
            reasons.append("horizon check unavailable: tilt and roll not checked independently")
        if chain_ok and chain.get("vs_direct_px") is not None and chain["vs_direct_px"] > CHAIN_TOL_PX:
            verdict = "DELIVERED-WIDER"
            reasons.append(f"chain differs from direct by {chain['vs_direct_px']:.1f} px "
                           f"({chain['vs_direct_deg']:.3f} deg; > {CHAIN_TOL_PX} px): uncertainty at least that")
        m = (aug or {}).get("vs_direct_ground_m") or {}
        if aug_ok and m.get("median") is not None and m["median"] > AUG_TOL_M:
            verdict = "DELIVERED-WIDER"
            reasons.append(f"B_aug differs from B by {m['median']:.2f} m on the beach (> {AUG_TOL_M} m)")
    elif aug_ok:
        verdict, source = "DELIVERED", "Aug->Oct"
        notes.append("direct Jan->Oct does not count (" + link_why(direct) + "): B = the Aug 2024 solve")
        if link_hz(aug["link"]) != "pass":
            verdict = "DELIVERED-WIDER"
            reasons.append("horizon check unavailable on the Aug->Oct solve: tilt and roll not checked independently")
        m = aug.get("vs_chain_ground_m")
        if chain_ok and m is not None and m > AUG_TOL_M:
            verdict = "DELIVERED-WIDER"
            reasons.append(f"the c2 chain differs from B_aug by {m:.2f} m on the beach (> {AUG_TOL_M} m)")
    elif chain_ok:
        verdict, source = "DELIVERED-WIDER", "chain Jan->Dec->Oct"
        reasons.append("direct: " + link_why(direct))
        reasons.append("B = the chain (direct failed its gates; two links: never better than WIDER)")
        if any(link_hz(chain.get(k)) != "pass" for k in ("link1", "link2")):
            reasons.append("a chain link's horizon check is unavailable: tilt and roll not checked independently")
    else:
        verdict, source = "NOT DELIVERED", None
        reasons.append("direct: " + link_why(direct))
        if chain is not None:
            reasons.append("chain: " + "; ".join(f"{k}: {link_why(chain.get(k))}" for k in ("link1", "link2")
                                                 if k in chain) + ("" if "link2" in chain else "; link 2 not run"))
        if aug is not None:
            reasons.append("Aug->Oct: " + link_why(aug.get("link")))
    return verdict, source, reasons, notes


def run_camera(cam, roots, aug_dir, log):
    """Every route of one camera, then decide(). -> (res, B)"""
    io_p, ref_p = archived(CAL / IO_FILE.format(cam=cam)), archived(CAL / REF_EO_FILE["jan"].format(cam=cam))
    io, eo_jan = load_intrinsics(io_p), load_extrinsics(ref_p)
    sel = {}
    for k in ("oct", "jan") + (("dec",) if cam == "c2" else ()):
        sel[k], _ = select_photos(roots, cam, WINDOWS[k], io)
    if aug_dir:
        sel["aug"], _ = select_photos([aug_dir], cam, WINDOWS["aug"], io)
    res = {"camera": cam, "io": io_p.name, "io_sha256": sha256(io_p), "reference_eo": ref_p.name,
           "reference_eo_sha256": sha256(ref_p),
           "photos_selected": {k: {d: [Path(p).name for p in v] for d, v in s.items()} for k, s in sel.items()}}
    log(f"\n{cam}: " + ", ".join(f"{k} {sum(len(v) for v in s.values())} photos on {len(s)} days" for k, s in sel.items()))
    Qd, Bd, direct = run_link(cam, io, eo_jan, sel["jan"], ANCHOR_DAY["jan"], sel["oct"], "direct Jan->Oct", log)
    direct["has_B"] = Bd is not None
    res["direct"] = direct
    Bs = {"direct Jan->Oct": Bd}
    if cam == "c2" and sel.get("dec"):
        Q1, Bdec, l1 = run_link(cam, io, eo_jan, sel["jan"], ANCHOR_DAY["jan"], sel["dec"], "chain Jan->Dec", log)
        chain = {"link1": l1, "has_B": False}
        if Bdec is not None and link_counts(l1):
            dec_days = {d: v for d, v in sel["dec"].items()
                        if not any(r["day"] == d for r in l1["per_day"]["days"] if r.get("diff_px", 0) > DAY_TOL_PX)}
            Q2, Bch, l2 = run_link(cam, io, Bdec, dec_days, max(dec_days), sel["oct"], "chain Dec->Oct", log)
            chain["link2"] = l2
            if Bch is not None:
                chain.update(has_B=True, eo=l2["eo"])
                Bs["chain Jan->Dec->Oct"] = Bch
                if Bd is not None:
                    chain["vs_direct_px"] = round(disp_px(cirn_angles_to_R(*Bch[3:6]), cirn_angles_to_R(*Bd[3:6]), io), 2)
                    chain["vs_direct_deg"] = round(dist_deg(cirn_angles_to_R(*Bch[3:6]), cirn_angles_to_R(*Bd[3:6])), 4)
        res["chain"] = chain
    if aug_dir:
        aug_p = archived(CAL / REF_EO_FILE["aug"].format(cam=cam))
        eo_aug = load_extrinsics(aug_p)
        Qa, Ba, la = run_link(cam, io, eo_aug, sel["aug"], ANCHOR_DAY["aug"], sel["oct"], "Aug->Oct", log)
        aug = {"link": la, "reference_eo": aug_p.name, "reference_eo_sha256": sha256(aug_p), "has_B": Ba is not None}
        if Ba is not None:
            Bs["Aug->Oct"] = Ba
            if Bd is not None:
                m, p90 = ground_diff_m(io, Bd, Ba)
                aug["vs_direct_ground_m"] = {"median": None if m is None else round(m, 3),
                                             "p90": None if p90 is None else round(p90, 3)}
                aug["vs_direct_px"] = round(disp_px(cirn_angles_to_R(*Ba[3:6]), cirn_angles_to_R(*Bd[3:6]), io), 2)
            if "chain Jan->Dec->Oct" in Bs:
                m, _ = ground_diff_m(io, Ba, Bs["chain Jan->Dec->Oct"])
                aug["vs_chain_ground_m"] = None if m is None else round(m, 3)
        res["aug"] = aug
    verdict, source, reasons, notes = decide(res)
    res.update(verdict=verdict, source=source, reasons=reasons, notes=notes)
    B = Bs.get(source) if source else None
    if source == "Aug->Oct":
        res["b_reference_eo"], res["b_reference_eo_sha256"] = res["aug"]["reference_eo"], res["aug"]["reference_eo_sha256"]
    elif source is not None:
        res["b_reference_eo"], res["b_reference_eo_sha256"] = res["reference_eo"], res["reference_eo_sha256"]
    if B is not None:
        res["B"] = eo_deg(B)
        res["pan_change_deg"] = round(float(((np.degrees(B[3] - eo_jan[3]) + 180) % 360) - 180), 4)
        res["uncertainty"] = stated_uncertainty(res, io)
    return res, B


def link_u_px(lk):
    """A solve's own spread (px): the larger of its leave-one-day-out change and the median per-day difference."""
    s, pd = (lk or {}).get("solve", {}) or {}, (lk or {}).get("per_day", {}) or {}
    diffs = [r["diff_px"] for r in pd.get("days", []) if "diff_px" in r]
    return float(max([s.get("jackknife_max_px") or 0.0] + ([float(np.median(diffs))] if diffs else [])))


def stated_uncertainty(res, io):
    """The two pre-registered uncertainty terms B carries (see the docstring). Neither is a validated bound on
    real photos."""
    src, cam = res["source"], res["camera"]
    if src == "chain Jan->Dec->Oct":
        spread = link_u_px(res["chain"]["link1"]) + link_u_px(res["chain"]["link2"])
    else:
        spread = link_u_px({"direct Jan->Oct": res.get("direct"), "Aug->Oct": (res.get("aug") or {}).get("link")}[src])
    rot = max(ROT_SPREAD_FACTOR * spread, ROT_FLOOR_PX)
    why = [f"max({ROT_SPREAD_FACTOR:g} x the {src} solve's spread {spread:.2f} px, floor {ROT_FLOOR_PX:g} px)"]
    ch = res.get("chain") or {}
    if (src == "direct Jan->Oct" and link_counts(ch.get("link1")) and link_counts(ch.get("link2"))
            and (ch.get("vs_direct_px") or 0) > CHAIN_TOL_PX and ch["vs_direct_px"] > rot):
        rot = ch["vs_direct_px"]
        why.append("the chain-direct disagreement")
    out = {"rotation_px": round(rot, 2), "rotation_deg": round(float(np.degrees(rot / io[4])), 3),
           "spread_px": round(spread, 2), "rotation_rule": " and ".join(why),
           "position": dict(POSITION_TERM[cam], reference_solution_spread=POSITION_SPREAD[cam],
                            rule="the largest DEM effect over the five synthetic position-inheritance cases (not a "
                                 "bound)"),
           "not_validated": UNCERTAINTY_NOT_VALIDATED}
    aug = res.get("aug") or {}
    for k in ("vs_direct_ground_m", "vs_chain_ground_m"):
        v = aug.get(k)
        m = v.get("median") if isinstance(v, dict) else v
        if m is not None and m > AUG_TOL_M:
            out["beach_m"] = round(m, 2)
    if link_hz({"direct Jan->Oct": res.get("direct"), "Aug->Oct": aug.get("link"),
                "chain Jan->Dec->Oct": ch.get("link2")}[src]) != "pass":
        out["note"] = "tilt and roll not checked by the horizon"
    return out


def uncertainty_text(u, cam):
    if not u:
        return "n/a"
    p = u["position"]
    return (f"rotation +/- {u['rotation_px']} px on the beach part of the frame (~{u['rotation_deg']} deg; "
            f"{u['rotation_rule']})"
            + (f", {u['beach_m']} m on the beach" if "beach_m" in u else "") + (f"; {u['note']}" if "note" in u else "")
            + f"; inherited camera position: in the largest of five synthetic position-inheritance cases, DEM |dz| "
              f"median {p['dem_abs_median_m']} m, p90 {p['dem_p90_m']} m, horizontal {p['horizontal_median_m']} m (not "
              f"a bound; the {cam} GCP solutions of Aug 2024 and Jan 2025 are {p['reference_solution_spread']} apart); "
              f"{u['not_validated']}")


# ------------------------------------------------------------------ runs and outcome

def prior_runs(out_dir, finished=True):
    """Runs under out_dir/runs (finished ones only, or every one started), oldest first: [(run.json dict,
    folder)]. A folder whose run.json cannot be read counts as a run of unknown pre-registration."""
    found = []
    for d in sorted(p for p in (Path(out_dir) / "runs").glob("*") if p.is_dir()):
        try:
            js = json.loads((d / RUN_FILE).read_text())
        except (OSError, ValueError):
            js = {"run_id": d.name, "status": "unreadable"}
        if not finished or js.get("status") == "finished":
            found.append((js, d))
    return sorted(found, key=lambda t: t[0].get("started_utc") or t[0].get("run_id") or "")


def earlier_preregistration_problem(out_dir, prereg_info):
    """A run under a pre-registration that differs from the one any earlier run (trial or delivery, finished or
    not) was made under: the rule would have been changed after a result could be seen. -> text or ''."""
    other = [(js, d) for js, d in prior_runs(out_dir, finished=False)
             if (js.get("preregistration") or {}).get("sha256") != prereg_info.get("sha256")]
    if not other:
        return ""
    return (f"{len(other)} earlier run(s) were made under another version of the pre-registration (or one that "
            f"cannot be read): {', '.join(d.name for _, d in other[:4])}. A pre-registration cannot change once a run "
            f"was made under it: restore it as it was (git checkout), or the outcome rule would be chosen after "
            f"seeing a result")


def plan_run(out_dir, eligible, cams, aug_dir, prereg_info):
    """What this run may do, by the pre-registered rerun rule (docstring: RUNS AND HISTORY).
    -> {"run_no", "solve": [cams], "keep": {cam: (report dict, folder, run id)}}. Raises Refused."""
    if not eligible:
        return {"run_no": None, "solve": list(cams), "keep": {}, "variant": None}
    done = [(js, d) for js, d in prior_runs(out_dir) if js.get("eligible")]
    if not done:
        return {"run_no": 1, "solve": list(cams), "keep": {}, "variant": B_VARIANT}
    first, fdir = done[0]
    head = f"a delivery run already finished: {first.get('run_id')} ({fdir})"
    if len(done) >= 2:
        raise Refused(f"{head}, and a second one, {done[1][0].get('run_id')}: the pre-registered rule allows no "
                      f"third delivery run (a TRIAL run, with other --roots or one --camera, writes reports only)")
    if first.get("aug_dir"):
        raise Refused(f"{head}, with Aug photos: the pre-registered rule allows no other delivery run")
    if not aug_dir:
        raise Refused(f"{head}: B is fixed by it. The only other delivery run allowed brings new evidence: "
                      f"--aug-dir with photos of 25-29 Aug 2024, for a camera it did NOT DELIVER")
    if (first.get("preregistration") or {}).get("sha256") != prereg_info.get("sha256"):
        raise Refused(f"{head}, under another pre-registration ({(first.get('preregistration') or {}).get('sha256', '')[:12]}"
                      f"..., now {str(prereg_info.get('sha256'))[:12]}...): a pre-registration cannot change after "
                      f"the first delivery run")
    verdicts = {c: (first.get("cameras") or {}).get(c, {}).get("verdict") for c in cams}
    solve_ = [c for c in cams if verdicts.get(c) in (None, "NOT DELIVERED")]
    if not solve_:
        raise Refused(f"{head}, and it passed the gates for {', '.join(cams)}: nothing is left for the Aug photos "
                      f"to solve (B stays the first run's)")
    keep = {}
    for c in cams:
        if c in solve_:
            continue
        rp = fdir / (OUT_REPORT.format(cam=c) + ".json")
        try:
            keep[c] = (json.loads(rp.read_text()), fdir, first.get("run_id"))
        except (OSError, ValueError) as exc:
            raise Refused(f"the first delivery run's report {rp} cannot be read ({exc}): its {c} result cannot be "
                          f"kept, so no second run")
    first_delivered = [c for c, v in (first.get("cameras") or {}).items() if v.get("delivery") == "DELIVERED"]
    return {"run_no": 2, "solve": solve_, "keep": keep, "first_run": first.get("run_id"),
            "first_delivered": first_delivered,
            # the second run's B is THE B only when the first delivered no camera; else a second product
            "variant": B2_VARIANT if first_delivered else B_VARIANT}


def apply_outcome(results, rule, eligible):
    """Per camera 'delivery' from the gate verdicts and the outcome rule (docstring: OUTCOME RULE)."""
    passed = [c for c, r in results.items() if r["verdict"] != "NOT DELIVERED"]
    for c, r in results.items():
        if not eligible:
            r["delivery"] = "TRIAL (not delivered)"
        elif r["verdict"] == "NOT DELIVERED":
            r["delivery"] = "NOT DELIVERED"
        elif rule == "both-or-none" and len(passed) < len(results):
            r["delivery"] = (f"WITHHELD (outcome rule both-or-none: "
                             f"{', '.join(x for x in results if x not in passed)} NOT DELIVERED)")
        else:
            r["delivery"] = "DELIVERED"
    return results


def pan_check(results, kept=()):
    """c1/c2: both passed -> their pan changes must agree within PAN_TOL_DEG, else both NOT DELIVERED. In the
    second delivery run a camera kept from the first keeps its verdict: only the camera solved now fails."""
    if len(results) != 2:
        for r in results.values():
            r["notes"].append("c1/c2 pan check not run (one camera only)")
        return
    if all(r["verdict"] != "NOT DELIVERED" for r in results.values()):
        dp = abs(results["c1"]["pan_change_deg"] - results["c2"]["pan_change_deg"])
        for r in results.values():
            r["c1_c2_pan_difference_deg"] = round(dp, 4)
        if dp > PAN_TOL_DEG:
            for c, r in results.items():
                if c in kept:
                    r["notes"].append(f"c1 and c2 pan changes differ by {dp:.2f} deg (> {PAN_TOL_DEG}): this camera's "
                                      f"verdict is the first run's and stands; the camera solved now is not delivered")
                    continue
                r["verdict"] = "NOT DELIVERED"
                r["reasons"].append(f"c1 and c2 pan changes differ by {dp:.2f} deg (> {PAN_TOL_DEG})")
    else:
        for r in results.values():
            r["notes"].append("c1/c2 pan check unavailable (one camera not delivered)")


def write_b(run_dir, cam, res, run_info):
    """The B file of a delivered camera (comment lines: everything survey_products.py checks)."""
    B = eo_rad(res["B"])
    a = angles_deg(B)
    rep = OUT_REPORT.format(cam=cam)
    src = res["source"]
    d = {"direct Jan->Oct": res.get("direct"), "chain Jan->Dec->Oct": (res.get("chain") or {}).get("link2"),
         "Aug->Oct": (res.get("aug") or {}).get("link")}[src] or {}
    pr = run_info["preregistration"]
    notes = [
        f"CALIBRATION B (GCP-FREE) for the 23 Oct 2024 DEM, camera {cam}: {res['verdict']} "
        f"({'; '.join(res['reasons']) or 'all checks passed'})",
        f"Delivery: {res['delivery']} under the pre-registered outcome rule '{pr.get('outcome_if_one_camera')}' "
        f"(run {run_info['run_id']}, delivery run {run_info['run_no']}"
        + (f"; solved in run {res['solved_in_run']}" if res.get("solved_in_run") != run_info["run_id"] else "") + ").",
        f"No Oct 2024 GCP used: CACO03_{cam}_20241023_EO.yaml (calibration A, solved from them) was not read.",
        f"Reference EO: {res['b_reference_eo']} ({'Aug 2024' if src == 'Aug->Oct' else 'Jan 2025'} GCPs, not the "
        f"Oct 2024 ones) sha256 {res['b_reference_eo_sha256']}; lens {res['io']}.",
        f"Position (x, y, z) kept from {res['b_reference_eo']}: B inherits that solution's camera position (the "
        f"CACO03 GCP solutions of Aug 2024 and Jan 2025 place the {cam} camera {POSITION_SPREAD[cam]} apart on one "
        f"mount).",
        f"Rotation source: {src}; pan change {res['pan_change_deg']:+.3f} deg from {res['reference_eo']}.",
        f"Solve ({d.get('label')}): {(d.get('solve') or {}).get('inliers')} inliers of "
        f"{(d.get('solve') or {}).get('matches')} matches, residual median {(d.get('solve') or {}).get('resid_med_px')} px, "
        f"inlier days {(d.get('solve') or {}).get('inlier_new_days')}, leave-one-day-out "
        f"{(d.get('solve') or {}).get('jackknife_max_px')} px.",
        f"Per-day: {(d.get('per_day') or {}).get('measured')} days measured, off days {(d.get('per_day') or {}).get('off_days')}.",
        "Horizon (differential, in the band [0, expected bias] widened by the noise tolerance): "
        + horizon_text(d.get("horizon")) + ".",
    ]
    if "chain" in res:
        notes.append(f"c2 chain Jan->Dec->Oct: counts={link_counts((res['chain'] or {}).get('link1')) and link_counts((res['chain'] or {}).get('link2'))}, "
                     f"vs direct {res['chain'].get('vs_direct_px')} px ({res['chain'].get('vs_direct_deg')} deg).")
    if "aug" in res:
        notes.append(f"Aug 2024 solve: counts={link_counts(res['aug'].get('link'))}, vs B "
                     f"{res['aug'].get('vs_direct_ground_m')} m.")
    ref_eo = load_extrinsics(archived(CAL / res["b_reference_eo"]))
    turn = dist_deg(cirn_angles_to_R(*B[3:6]), cirn_angles_to_R(*ref_eo[3:6]))
    route = {"direct Jan->Oct": "Jan 2025 -> Oct 2024 photos", "chain Jan->Dec->Oct": "Jan 2025 -> Dec 2024 -> Oct 2024 photos",
             "Aug->Oct": "Aug 2024 -> Oct 2024 photos"}[src]
    notes.append(f"Solved by gcpfree_eo.py: {res['b_reference_eo']} rotated by {turn:.3f} deg ({route}); position unchanged.")
    u_txt = uncertainty_text(res.get("uncertainty"), cam)
    if res["verdict"] == "DELIVERED":
        notes.append(f"acceptance: PASS (pre-registered gates F, S and the horizon band passed; {u_txt})")
    else:
        notes.append(f"acceptance: PASS-WIDER {u_txt}: {'; '.join(res['reasons'])}")
    notes.append(f"Uncertainty (stated; not a validated bound on real photos): {u_txt}.")
    notes.append(f"Oct photos: {', '.join(d.get('new_photos', []))}")
    notes.append(f"Reference photos: {', '.join(d.get('ref_photos', []))}")
    notes.append(f"Run {run_info['run_id']}: {run_info['argv_text']}; photo folders {', '.join(run_info['roots'])}"
                 + (f"; Aug photos {run_info['aug_dir']}" if run_info.get("aug_dir") else "") + ".")
    var = run_info.get("variant")
    notes.append(f"Variant: build this B with survey_products.py --calibration-variant {var} ("
                 + ("the delivered B" if var == B_VARIANT else
                    f"a SECOND product, from delivery run 2 with Aug 2024 photos; the B of run "
                    f"{run_info.get('first_run')} stays the delivered one") + ").")
    notes.append(f"Made by gcpfree_eo.py (sha256 {run_info['procedure_sha256']}) on {run_info['started_utc']} UTC; "
                 f"report {rep}.json; pre-registration {pr.get('path')} sha256 {pr.get('sha256')} "
                 f"(outcome_if_one_camera: {pr.get('outcome_if_one_camera')}); pre-registered constants sha256 "
                 f"{run_info.get('constants_sha256')}.")
    notes.append(EXPOSURE)
    yml = run_dir / OUT_EO.format(cam=cam)
    yml.write_text(
        f"x: {B[0]:.4f}\ny: {B[1]:.4f}\nz: {B[2]:.4f}\n"
        f"azimuth: {a[0]:.4f}\ntilt: {a[1]:.4f}\nroll: {a[2]:.4f}\n"
        + "".join(f"# {n}\n" for n in notes)
        + "#x, y - camera Easting / Northing (m, UTM zone 19); z - camera elevation (m, NAVD88)\n"
          "#azimuth, tilt, roll - degrees, CIRN convention\n")
    return yml


def new_run_dir(out_dir, eligible):
    rid = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + ("" if eligible else "_trial")
    d = Path(out_dir) / "runs" / rid
    try:
        d.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        raise Refused(f"{d} exists: a run never writes over an earlier one (wait a second and run again)")
    return rid, d


def finish_run(run_dir, run_info, results, lines, log):
    """Reports for every camera, B files for the delivered ones, run.json. -> exit status."""
    run_info["loaded_modules"] = loaded_modules()        # every repository module the run had loaded (sha256)
    for cam, r in results.items():
        r["run"] = {k: run_info.get(k) for k in ("run_id", "run_no", "eligible", "trial_why", "variant", "first_run",
                                                 "argv", "roots", "aug_dir", "procedure_sha256", "modules_sha256",
                                                 "constants_sha256", "preregistration", "started_utc")}
        r.setdefault("solved_in_run", run_info["run_id"])
        r["exposure"] = EXPOSURE
        if r["delivery"] == "DELIVERED":
            r["b_file"] = write_b(run_dir, cam, r, run_info).name
            r["b_file_sha256"] = sha256(run_dir / r["b_file"])
            log(f"  wrote {run_dir / r['b_file']} (sha256 {r['b_file_sha256']})")
        (run_dir / (OUT_REPORT.format(cam=cam) + ".json")).write_text(json.dumps(r, indent=1, default=str) + "\n")
        (run_dir / (OUT_REPORT.format(cam=cam) + ".txt")).write_text("\n".join(lines) + "\n")
    run_info.update(status="finished", finished_utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
                    cameras={c: {"verdict": r["verdict"], "delivery": r["delivery"], "source": r.get("source"),
                                 "solved_in_run": r.get("solved_in_run"), "b_file": r.get("b_file"),
                                 "b_file_sha256": r.get("b_file_sha256")} for c, r in results.items()})
    (run_dir / RUN_FILE).write_text(json.dumps(run_info, indent=1, default=str) + "\n")
    (run_dir / "run_report.txt").write_text("\n".join(lines) + "\n")
    dl = [r["delivery"] for r in results.values()]
    if run_info["eligible"] and all(x == "DELIVERED" for x in dl) and len(dl) == 2 \
            and all(r["verdict"] == "DELIVERED" for r in results.values()):
        return 0
    return 3 if any(x == "DELIVERED" for x in dl) else 4


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--camera", choices=["c1", "c2", "both"], default="both",
                    help="both (default; a delivery run needs both); one camera makes a TRIAL run")
    ap.add_argument("--roots", nargs="+", default=None,
                    help="photo folders (default: the station's three, ROOTS; any other folders make a TRIAL run)")
    ap.add_argument("--aug-dir", default=None, help="CACO03 photos of 25-29 Aug 2024 from Chris: " + AUG_NAMING)
    ap.add_argument("--dry-run", action="store_true",
                    help="the pre-registration, the photos, the reference horizon; writes nothing")
    ap.add_argument("--preregister", action="store_true", help="write the pre-registration file (then commit it)")
    ap.add_argument("--outcome-if-one-camera", choices=OUTCOME_RULES, default=None,
                    help="with --preregister: the owner's choice (see the docstring)")
    ap.add_argument("--decided-by", default=None, help="with --preregister: who chose the outcome rule")
    args = ap.parse_args(argv)
    try:
        return run(args, list(sys.argv if argv is None else ["gcpfree_eo.py"] + list(argv)))
    except Refused as exc:
        print(f"REFUSED: {exc}")
        return 2


def check_photos(roots, cams, aug_dir, eligible):
    """Bad input is refused before a run starts, never turned into a verdict: a camera solved with no Oct 2024
    or no Jan 2025 photo (or, with --aug-dir, no Aug 2024 photo); for a delivery run, a missing photo folder.
    -> {cam: {window: photo count}}. Raises Refused."""
    found = {}
    for cam in cams:
        io = load_intrinsics(archived(CAL / IO_FILE.format(cam=cam)))
        found[cam] = {}
        for k in ("oct", "jan") + (("aug",) if aug_dir else ()):
            days, counts = select_photos([aug_dir] if k == "aug" else roots, cam, WINDOWS[k], io)
            n = sum(len(v) for v in days.values())
            found[cam][k] = n
            if counts.get("missing_roots") and (eligible or k == "aug"):
                raise Refused(f"photo folder(s) missing: {', '.join(counts['missing_roots'])}: "
                              + ("a delivery run uses every one of the station's folders (an unmounted disk would "
                                 "change the photos it uses)" if k != "aug" else "--aug-dir is not a folder"))
            if not n:
                what = {"oct": "20-26 Oct 2024 (the DEM window)", "jan": "18-23 Jan 2025 (the reference days)",
                        "aug": "25-29 Aug 2024"}[k]
                raise Refused(f"{cam}: no usable {STATION} {cam} photo of {what} in "
                              f"{aug_dir if k == 'aug' else ', '.join(roots)} ({counts.get('files', 0)} "
                              f"'{cam}.timex.jpg' files with an epoch name, {counts.get('other_station', 0)} of another "
                              f"station, {counts.get('outside_window', 0)} outside the days, {counts.get('wrong_size', 0)} "
                              f"reduced or of another size). Photos must be named {AUG_NAMING.replace('25-29 Aug 2024', 'the days above')}. "
                              f"Not run: bad input is never a verdict")
    return found


def trial_reasons(args, roots):
    """Why a run is a TRIAL run (reports only), [] for a DELIVERY run: one camera, or photo folders other than the
    station's. (Every run, a trial too, is made under this repository's committed pre-registration, PREREG.)"""
    why = []
    if args.camera != "both":
        why.append("one camera only")
    if roots != [str(r) for r in ROOTS]:
        why.append(f"photo folders other than the station's ({', '.join(ROOTS)})")
    return why


def photos_problem(prereg_info, roots, cams):
    """On the station's folders, the photos a run of these cameras selects must be the pre-registered ones.
    -> (text or '', the manifest selected now)."""
    now = photo_manifest(roots, cams)
    if roots != [str(r) for r in ROOTS]:
        return "", now
    want = [r for r in (prereg_info.get("photos") or []) if r[0] in cams]
    diff = manifest_difference(want, now)
    return ((f"the photos a run selects now from the station's folders are not the pre-registered ones ({diff}): a run "
             f"uses exactly the photos its pre-registration lists (a photo removed, replaced or added after the "
             f"pre-registration would change the result)") if diff else ""), now


def run(args, argv):
    roots = [str(r) for r in (args.roots if args.roots is not None else ROOTS)]
    for p in roots + ([args.aug_dir] if args.aug_dir else []):
        guard(p)
    if args.preregister:
        p = write_preregistration(PREREG, args.outcome_if_one_camera, args.decided_by, roots)
        print(f"wrote {p} (the outcome rule '{args.outcome_if_one_camera}' and "
              f"{len(read_preregistration(p)['photos'])} photos)\nNow commit AND push it before any run (a trial too; "
              f"every run is refused until its commit is on a remote branch): git -C {p.parent} add {p.name} && "
              f"git -C {p.parent} commit -m 'gcpfree: pre-registration ({args.outcome_if_one_camera})' && "
              f"git -C {p.parent} push\nThen send Chris its sha256 (sha256sum {p}) and the commit hash")
        return 0
    import cv2
    if not hasattr(cv2, "SIFT_create"):
        print(f"OpenCV {cv2.__version__} has no SIFT (needs 4.4 or later): the pre-registered procedure "
              f"cannot run here (no other detector may replace it).")
        return 2
    cams = ["c1", "c2"] if args.camera == "both" else [args.camera]
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    prereg_info, probs = preregistration(PREREG)        # every run, a trial too: this repository's, committed
    out_dir = Path(PREREG).parent                        # runs/<run id>/ next to it: every run is recorded there
    trial_why = trial_reasons(args, roots)
    eligible = not trial_why
    if args.dry_run:
        return dry_run(args, cams, roots, prereg_info, probs, eligible, out_dir, log)
    if probs:
        raise Refused("the pre-registration does not allow a run (a trial neither: the rule is fixed and committed "
                      "before this code solves any photo): " + "; ".join(probs))
    bad = earlier_preregistration_problem(out_dir, prereg_info) or foreign_problem(out_dir)
    if bad:
        raise Refused(bad)
    for cam in cams:                                     # every file read is checked before anything is solved
        archived(CAL / IO_FILE.format(cam=cam))
        archived(CAL / REF_EO_FILE["jan"].format(cam=cam))
        if args.aug_dir:
            archived(CAL / REF_EO_FILE["aug"].format(cam=cam))
    plan = plan_run(out_dir, eligible, cams, args.aug_dir, prereg_info)
    photos = check_photos(roots, plan["solve"], args.aug_dir, eligible)
    bad, manifest = photos_problem(prereg_info, roots, plan["solve"])
    if bad:
        raise Refused(bad)
    aug_manifest = []
    if args.aug_dir:
        for cam in plan["solve"]:
            io = load_intrinsics(archived(CAL / IO_FILE.format(cam=cam)))
            days, _ = select_photos([args.aug_dir], cam, WINDOWS["aug"], io)
            aug_manifest += [[cam, "aug", dd, Path(q).name, os.path.getsize(str(guard(q))), sha256(q)]
                             for dd, qs in sorted(days.items()) for q in qs]
    rid, run_dir = new_run_dir(out_dir, eligible)
    started = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    run_info = {"run_id": rid, "run_no": plan["run_no"], "eligible": eligible, "trial_why": trial_why or None,
                "variant": plan.get("variant"), "first_run": plan.get("first_run"), "status": "running",
                "started_utc": started, "argv": argv, "argv_text": " ".join(shlex.quote(a) for a in argv),
                "roots": roots, "camera": args.camera, "aug_dir": args.aug_dir, "photos_found": photos,
                "procedure_sha256": own_sha256(), "modules_sha256": module_sha256(), "preregistration": prereg_info,
                "constants": constants(), "constants_sha256": constants_sha256(),
                "constants_problem": constants_problem() or None, "exposure": EXPOSURE,
                "solved": plan["solve"], "kept": {c: v[2] for c, v in plan["keep"].items()},
                # the photos this run selected (name, size, sha256): on the station's folders the pre-registered ones
                "photos_used": manifest, "photos_used_sha256": photos_sha256(manifest),
                "photos_preregistered": roots == [str(r) for r in ROOTS],
                "aug_photos_used": aug_manifest or None,
                "aug_photos_sha256": photos_sha256(aug_manifest) if aug_manifest else None}
    (run_dir / RUN_FILE).write_text(json.dumps(run_info, indent=1, default=str) + "\n")
    log(f"gcpfree_eo.py run {rid} ({'DELIVERY run ' + str(plan['run_no']) if eligible else 'TRIAL run: reports only'}"
        f"; pre-registration {prereg_info.get('path')} ({len(prereg_info.get('git_history') or [])} commit(s)), "
        f"outcome_if_one_camera {prereg_info.get('outcome_if_one_camera')})")
    if eligible:
        log(f"  B of this run, if any: --calibration-variant {plan['variant']}"
            + ("" if plan["variant"] == B_VARIANT else
               f" (a second product: the B of run {plan.get('first_run')} stays the delivered one)"))
    else:
        log("  TRIAL: " + "; ".join(trial_why) + ": nothing it finds is delivered")
    if run_info["constants_problem"]:
        log(f"  NOTE: {run_info['constants_problem']}")
    results = {}
    for cam in cams:
        if cam in plan["keep"]:
            r, fdir, frid = plan["keep"][cam]
            r = {k: v for k, v in r.items() if k not in ("delivery", "run", "b_file", "b_file_sha256")}
            r["notes"] = list(r.get("notes") or []) + [f"kept from the first delivery run {frid} (not solved again)"]
            r["solved_in_run"] = r.get("solved_in_run") or frid
            log(f"\n{cam}: kept from run {frid}: {r['verdict']}")
            results[cam] = r
            continue
        res, B = run_camera(cam, roots, args.aug_dir, log)
        res["solved_in_run"] = rid
        results[cam] = res
    pan_check(results, kept=set(plan["keep"]))
    apply_outcome(results, prereg_info.get("outcome_if_one_camera"), eligible)
    log(f"\nVERDICT (pre-registered procedure; calibration A not read; outcome rule "
        f"{prereg_info.get('outcome_if_one_camera')})")
    for cam, r in results.items():
        log(f"  {cam}: {r['verdict']} -> {r['delivery']}"
            + (f"  B = {r['B'][3]:.4f} / {r['B'][4]:.4f} / {r['B'][5]:.4f} deg (source {r['source']})" if r.get("B") else "")
            + ("" if not r["reasons"] else "  -- " + "; ".join(r["reasons"])))
    rc = finish_run(run_dir, run_info, results, lines, log)
    if eligible:
        log(f"\nrun folder {run_dir}: commit it now, in ONE commit, and push it (git add {run_dir} && git commit && "
            f"git push), and write down the sha256 of each B file (and send them to Chris), BEFORE the B DEM is built "
            f"and before anything is "
            f"compared with the A DEM: survey_products.py builds no B whose file, report or run.json is uncommitted, "
            f"changed or touched by more than one commit")
    else:
        log(f"\nTRIAL run folder {run_dir}: reports only, nothing delivered. Commit and push it too, at once (git add "
            f"{run_dir} && git commit && git push: the record of every run: survey_products.py lists it in B's "
            f"README, and builds no B DEM while a run folder is uncommitted)")
    (run_dir / "run_report.txt").write_text("\n".join(lines) + "\n")
    return rc


def dry_run(args, cams, roots, prereg_info, probs, eligible, out_dir, log):
    """Writes nothing. The pre-registration, what a run would be now, the photos, and the reference EO's
    horizon offset on its own anchor-day photos (no Oct photo is read for it)."""
    pp = pin_problem()
    log(f"procedure: gcpfree_eo.py sha256 {own_sha256()}: " + ("the version survey_products.py pins (GCPFREE_PINS)"
                                                               if not pp else "REFUSED: " + pp))
    cp = constants_problem()
    log(f"pre-registered constants: sha256 {constants_sha256()}: " + ("the pinned ones" if not cp else "NOT: " + cp))
    log(f"pre-registration {prereg_info.get('path')}: " + ("OK" if not probs else "a run (a trial too) is REFUSED: "
                                                              + "; ".join(probs)))
    if not probs:
        log(f"  outcome_if_one_camera {prereg_info.get('outcome_if_one_camera')} (decided by "
            f"{prereg_info.get('decided_by')} on {prereg_info.get('decided_on')}); commit {prereg_info.get('commit')}, "
            f"pushed: on {', '.join(prereg_info.get('remote_branches') or [])}; "
            f"{len(prereg_info.get('photos') or [])} photos pinned")
    log(f"this run would be a {'DELIVERY' if eligible else 'TRIAL'} run"
        + ("" if eligible else f" (reports only: {'; '.join(trial_reasons(args, roots))}; a trial is a FULL result: "
                               f"never run one on the station)"))
    bad = (earlier_preregistration_problem(out_dir, prereg_info) if not probs else "") or foreign_problem(out_dir)
    if bad:
        log(f"  REFUSED: {bad}")
    runs_ = prior_runs(out_dir, finished=False)
    if runs_:
        log(f"  runs already in {out_dir / 'runs'}: " + ", ".join(
            f"{d.name} ({'delivery' if js.get('eligible') else 'trial'}, {js.get('status')})" for js, d in runs_))
    if eligible:
        try:
            pl = plan_run(out_dir, eligible, cams, args.aug_dir, prereg_info)
            log(f"  delivery run {pl['run_no']}: solves {', '.join(pl['solve'])}"
                + (f"; keeps {', '.join(pl['keep'])} from the first run" if pl["keep"] else "")
                + f"; its B: --calibration-variant {pl['variant']}")
        except Refused as exc:
            log(f"  REFUSED: {exc}")
    try:                                                 # the photo checks a run makes before it starts
        check_photos(roots, cams, args.aug_dir, eligible)
        log(f"photos: every folder there and Oct and Jan photos for {', '.join(cams)}"
            + (" and Aug photos" if args.aug_dir else "") + ": a run would not be refused for its photos")
        if not probs:
            pb, now = photos_problem(prereg_info, roots, cams)
            log("  " + (f"a run would be REFUSED (exit 2): {pb}" if pb else
                        f"the {len(now)} photos selected now are the pre-registered ones (names, sizes, sha256)"
                        if roots == [str(r) for r in ROOTS] else
                        f"{len(now)} photos selected (other folders: not compared with the pre-registered ones)"))
        else:
            log(f"  {len(photo_manifest(roots, cams))} photos would be pinned by --preregister now (names, sizes, "
                f"sha256; nothing is computed from them)")
    except Refused as exc:
        log(f"photos: a {'delivery' if eligible else 'trial'} run would be REFUSED (exit 2): {exc}")
    for cam in cams:
        io = load_intrinsics(archived(CAL / IO_FILE.format(cam=cam)))
        eo_jan = load_extrinsics(archived(CAL / REF_EO_FILE["jan"].format(cam=cam)))
        for k in ("oct", "jan") + (("dec",) if cam == "c2" else ()) + (("aug",) if args.aug_dir else ()):
            days, counts = select_photos([args.aug_dir] if k == "aug" else roots, cam, WINDOWS[k], io)
            log(f"{cam} {k} {WINDOWS[k][0]}..{WINDOWS[k][1]}: {sum(len(v) for v in days.values())} photos on "
                f"{len(days)} days (missing folders: {counts.get('missing_roots')})"
                + (f"; NOT opened (named like a forbidden file, or a link to one): "
                   f"{', '.join(counts['forbidden_not_opened'])}" if counts.get("forbidden_not_opened") else ""))
            for d, ps in days.items():
                log(f"   {d}: " + " ".join(Path(p).name.split('.')[3][-8:] for p in ps))
            if k == "aug" and not days:
                log(f"   NO usable Aug photo: a run would be REFUSED. Names must be {AUG_NAMING}")
            if k == "jan" and days.get(ANCHOR_DAY["jan"]):
                shape = (int(io[1]), int(io[0]))
                anc = load_set({ANCHOR_DAY["jan"]: days[ANCHOR_DAY["jan"]]}, EER.land_mask(shape, io, eo_jan))
                h = horizon_fits(io, eo_jan, anc)
                log(f"   the Jan 2025 EO against the visible horizon on {ANCHOR_DAY['jan']} (no Oct photo read): "
                    + (f"tilt {h['dtilt_med']:+.3f}, roll {h['droll_med']:+.3f} deg ({h['n_clear']} clear of {len(anc)})"
                       if h["n_clear"] else f"no clear horizon ({len(anc)} photos)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())

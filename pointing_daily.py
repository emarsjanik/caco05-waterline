#!/usr/bin/env python3
"""
Camera Pointing Day By Day: A Log, And An Opt-In Correction
==============================================================
WHY. Every waterline is georectified with ONE calibration per period:
the live cron uses CACO05_<cam>_20251113_EO-CV.yaml for the whole
archive, a survey-date product the date's calibration for its week. A
camera that turns moves its lines with it: 0.1 deg of tilt moves a line
~6-9 m at 250-350 m range, ~0.13-0.17 m of DEM (docs/VALIDATION.md 4.6),
and c1 turned -1.2 deg of pan over 26-28 Sep 2026, in a storm, which only
a survey-date build found. survey_products.py measures each day's turn
(pointing_check.measure: land features matched against a reference
frame) and LEAVES OUT a day that turned; nothing measures the live
cameras day by day.

The coastal-imaging practice is to measure the turn on each image from
static land features and to correct the camera's extrinsics with it:
keypoints and sub-pixel correlation in static land sub-images, then the
pointing updated (Rodriguez-Padilla et al. 2020: movements of up to
10 px of pan, 30 px of tilt and 0.4 deg of roll at Anglet, up to 400 m
of rectification error at 2.5 km when left uncorrected; Pearre & Puleo
2009 corrected small azimuth and tilt changes by image correlation at
Rehoboth Beach). At Sisal, stabilisation cut the shoreline position
error by ~40% and the error of its trend from 0.25 to 0.01 m/yr (Arriaga
et al. 2022). View angles also drift with the sun on sunny days (Bouvier
et al. 2019), which a record of the turn frame by frame can test here.

THREE PARTS, none of which changes a product unless asked:

  1. THE GATES (day_correction). A day's feature match corrects that
     day only when it is a reading the survey products already trust,
     and more:
       * matched: features 'ok', 'moved' or 'MOVED', not 'ambiguous' (a
         second consensus at another turn: dirt or drops on the housing
         window match at zero turn whatever the camera did), not
         unmatched, not measured, unreadable;
       * on MATCH_INLIERS (40) or more inliers: good to ~0.1 deg
         (survey_products.py MATCH_ERR_DEG, from 402 planted turns and 69
         real pairs with a known answer);
       * turn_sigma TURN_SIGMA_MIN (4) or more: the turn is that many of
         its own standard errors from none, which a dirt consensus (at
         exactly zero) does not reach; a turn within its noise is no
         correction either;
       * its inliers spread MATCH_SPREAD_TESTED_PX (180 px) or more along
         the minor axis of their scatter: every match the 0.1 deg rests
         on did; a bunched match reads roll and tilt worse;
       * a turn of at most CORRECT_MAX_DEG (2.5 deg) in pan, tilt and
         roll: the largest planted turn the matching was tested on; the
         larger turns seen here were re-sets that also MOVED the camera
         (24 Jan 2025: ~15 deg and 5.6-5.8 m), which no rotation
         describes;
       * not the reference day itself (the others are measured against
         its frame), and, in the survey products, not UNCONFIRMED (the
         pointing step's two methods disagree, an ambiguous match, or a
         day nothing could read next to one left out:
         survey_products.pointing_verdicts).
     The survey products (--correct-pointing) add two gates of the day:
       * a turn of MORE than MATCH_ERR_DEG (0.1 deg) in pan, tilt or
         roll: a smaller one is within the match's own error against
         another day (<= 0.101 deg at 40-59 inliers), so a day at rest
         is left as the default build grids it, not moved by the
         matcher's noise (turn_sigma is a distance from no turn in the
         inliers' own scatter, not this error: real at-rest matches
         read 4.9-18.5). The corrected days are then the turned ones,
         those the default leaves out or keeps as a 'small change';
       * the pointing HELD through the day: the day's first and last
         photos of the detection (other than its probe frame and the
         reference frame) are matched against the same reference, and
         each must pass the gates above and read the probe's turn
         within MATCH_ERR_DEG in pan, tilt and roll. The correction is
         measured on one probe frame (nearest 17:00 UTC) and maps every
         row of the UTC day; on a day the camera turned or came back
         (wind load, a slip) the frames before the change would be
         moved by the whole turn. A frame that cannot be read (fog,
         rain) shows nothing held: the day is not corrected, and stays
         left out when the default leaves it out.
     The first three thresholds and MATCH_ERR_DEG are survey_products.py's
     own, passed in (correction_rules), so there is one definition of
     each.

  2. THE CORRECTED POINTING (corrected_eo, regeorectify). The camera
     only rotates, so its pointing on the day is the reference's turned
     by the measured change. pointing_check.measure fits the day's
     rotation R = Q R_ref (Kabsch inside RANSAC) and reports R's CIRN
     angles less the reference's, so the day's EO is the calibration's
     position with its angles plus that change: exact, to the 0.001 deg
     the pointing rows keep (~6-9 cm of line position in range, ~1.5 mm
     of DEM, at 250-350 m). The day's waterline points
     are mapped again with it (georectify.pixel_to_ground, each on its
     own elevation plane, written as georectify.py writes them); every
     other row is left as it was, byte for byte. The pinned files
     (georectify.py, estimate_eo_rotation.py, view_reproject.py, ...)
     are used, not changed: the corrected EO is a new object made here,
     written next to the waterlines it mapped. The correction aligns each
     day with the REFERENCE frame's pointing, taken as the calibration's
     (as the pointing step does): a constant offset of the whole window
     against the calibration (the sea horizon's, reported in the README)
     is not corrected.
     survey_products.py --correct-pointing corrects so every day that
     passes, in its own output root (<root>_pointing-corrected), keeps
     every other day as the default build does (a day measured to differ
     that fails a gate is still left out), and says which days, by how
     much, on which frames and how far their lines moved, in its README.
     Its compare step re-projects a corrected day's lines for the C = 0
     sensitivity with that day's EO (survey_compare.py --camera-eo-day).
     Only the georectification is corrected: the detector's search
     (the envelope of a historical era, placed with --eo; the live
     crop) is the calibration's, padded by 0.03 / 0.06 of the frame
     height (~61 / ~123 px: ~1.4 / ~2.7 deg of tilt for c1, fy ~2573
     px, ~1.0 / ~1.9 deg for c2, fy ~3681 px), so a corrected day
     tilted ~1 deg or more (or rolled so the band's ends move as far)
     may have part of its waterline band outside the search, where the
     detector can follow another edge (docs/VALIDATION.md 4.6b).

  3. THE LIVE LOG (main; waterline_timex_cron.sh). For each camera and
     each complete UTC day not yet logged (the last --lookback-days, at
     most --max-days a run): up to --frames archived timex spread over
     the day (--utc-hours), each matched against the camera's reference,
     gated as above (the feature gates only: the survey products also
     weigh the sea horizon and the neighbouring days), and the day's
     median turn over the frames that pass, the spread between them (a
     drift through the day), and the EO a correction WOULD use ->
     archive/pointing_daily_<cam>.json, plus one line per day in the
     cron's log. LOG ONLY: the live products keep the calibration as it
     is ('applied_live': false), and nothing is written where another
     job reads (not calibration/pointing/, owg.sh's monitor bank).
     The reference, in this order:
       * the frames of pointing_check.py's bank (calibration/pointing/
         <cam>/, the owner's) taken at the live calibration;
       * else this log's own (archive/pointing_daily_reference/<cam>/),
         made once from the frame nearest 17:00 UTC of --reference-day,
         or of the calibration's own day (the date in its file name) when
         the archive still holds it: then the changes are against the
         calibration's pointing, as survey_products' reference on the
         calibration day;
       * else nothing is measured, and the log says how to give one.
     A camera whose reference cannot be matched day after day (a foggy
     reference frame) shows 'no reading' every day: remove
     archive/pointing_daily_reference/<cam>/ and give --reference-day.

EVALUATION ON THE NUC: docs/VALIDATION.md section 4.6b (the corrected
and the default product of the same date against its survey; the log's
days against the survey products' pointing tables).

Usage:
    python3 pointing_daily.py                          the cron's call: new days, both cameras
    python3 pointing_daily.py --camera c1 --days 2026-09-25 2026-10-02     these days again
    python3 pointing_daily.py --reference-day 2026-10-07 --camera c2       a reference of that day
    python3 pointing_daily.py --self-test [--image <c2 timex.jpg>]
        planted turns of a synthetic scene (or of a real c2 frame) recovered through the gates, the
        corrected EO and the waterline points mapped again; no scipy needed (the NUC has none)

References:
  Arriaga, J., Medellin, G., Ojeda, E., Salles, P., 2022. Shoreline detection accuracy from video
    monitoring systems. Journal of Marine Science and Engineering 10(1), 95. doi:10.3390/jmse10010095
  Bouvier, C., Balouin, Y., Castelle, B., Holman, R., 2019. Modelling camera viewing angle deviation to
    improve nearshore video monitoring. Coastal Engineering 147, 99-106.
  Pearre, N.S., Puleo, J.A., 2009. Quantifying seasonal shoreline variability at Rehoboth Beach, Delaware,
    using automated imaging techniques. Journal of Coastal Research 25(4), 900-914. doi:10.2112/08-1029.1
  Rodriguez-Padilla, I., Castelle, B., Marieu, V., Morichon, D., 2020. A simple and efficient image
    stabilization method for coastal monitoring video systems. Remote Sensing 12(1), 70.
    doi:10.3390/rs12010070
"""

import os
import re
import sys
import csv
import json
import shutil
import hashlib
import argparse
from pathlib import Path
from datetime import datetime, timezone, timedelta

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from georectify import load_intrinsics, load_extrinsics, pixel_to_ground   # noqa: E402

CAL = HERE / "calibration"

# The largest turn (deg; the largest of pan, tilt and roll) a correction is made from. The feature match was
# measured on planted turns of 0.05-2.5 deg (survey_products.py MATCH_INLIERS/MATCH_ERR_DEG; pointing_check.py
# --self-test up to 2.5 deg of pan); a larger turn is outside what was tested, and the larger turns seen here
# were re-sets that also MOVED the camera (24 Jan 2025: ~15 deg and 5.6-5.8 m), which no rotation describes.
CORRECT_MAX_DEG = 2.5
# The correction's version (its gates, and how the corrected pointing is made and applied): part of
# survey_products' detection signature with --correct-pointing, bumped when either changes, so a corrected
# product is rebuilt.
CORRECTION_VERSION = 2
# The live log: the cron's calibrations and archive (waterline_timex_cron.sh: step 6 and ARCHIVE_IMG).
LIVE_EO = {c: CAL / f"CACO05_{c}_20251113_EO-CV.yaml" for c in ("c1", "c2")}
ARCHIVE = HERE / "archive"
FRAME_EPOCH = re.compile(r"^(\d{9,11})\.")


def _num(v):
    """A number from a row's column ('' or None, as an older build's row without it: None)."""
    try:
        return None if v in ("", None) else float(v)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------
# 1. The gates
# ---------------------------------------------------------------------

def correction_rules(min_inliers, min_turn_sigma, min_spread_px, max_turn_deg=CORRECT_MAX_DEG, min_turn_deg=None,
                     max_day_spread_deg=None):
    """The gates of day_correction, as numbers. The first three are the survey products' own (survey_products.py
    MATCH_INLIERS, TURN_SIGMA_MIN, MATCH_SPREAD_TESTED_PX), passed in so each has one definition; max_turn_deg is
    this module's CORRECT_MAX_DEG. The gates of the day, given by the survey products only (MATCH_ERR_DEG; the
    live log reads frames, not days to correct): min_turn_deg, the turn a day needs to be corrected at all, and
    max_day_spread_deg, how far the day's other frames (its row's day_check) may read from its probe."""
    r = {"min_inliers": int(min_inliers), "min_turn_sigma": float(min_turn_sigma),
         "min_spread_px": float(min_spread_px), "max_turn_deg": float(max_turn_deg)}
    if min_turn_deg is not None:
        r["min_turn_deg"] = float(min_turn_deg)
    if max_day_spread_deg is not None:
        r["max_day_spread_deg"] = float(max_day_spread_deg)
    return r


def station_rules():
    """correction_rules with survey_products.py's thresholds (the live log has no survey product to pass them)."""
    import survey_products as sp
    return correction_rules(sp.MATCH_INLIERS, sp.TURN_SIGMA_MIN, sp.MATCH_SPREAD_TESTED_PX)


def day_correction(row, rules, ref_day=None):
    """
    Whether one day's feature match may correct that day's pointing, and by how much. row: a day as the survey
    products' pointing step writes it (survey_products.feature_record's columns, and its 'verdict' once
    pointing_verdicts judged it; the live log's rows have no verdict). rules: correction_rules. ref_day: the day of
    the reference frame. -> (d_deg [pan, tilt, roll] or None, [why not, ...]): the change measured from the
    reference, when every gate (the module's docstring) passes; a row of an older build without turn_sigma or the
    inliers' spread fails those gates (not measured is not passed). With the rules' gates of the day: a turn of
    min_turn_deg or less is at rest, and the row's day_check must show the pointing held (day_held).
    """
    why = []
    if ref_day and row.get("date") == ref_day:
        return None, ["the reference day (the other days are measured against its frame)"]
    if row.get("verdict") == "UNCONFIRMED":
        why.append("UNCONFIRMED by the pointing step (its methods disagree, or nothing could read it)")
    f = row.get("features") or "not measured"
    if f not in ("ok", "moved", "MOVED"):
        why.append(f"no feature reading ({f})")
    n = _num(row.get("inliers"))
    if n is None or n < rules["min_inliers"]:
        why.append(f"{int(n or 0)} inliers, under {rules['min_inliers']}")
    ts = _num(row.get("turn_sigma"))
    if ts is None or ts < rules["min_turn_sigma"]:
        why.append("turn_sigma not measured" if ts is None else
                   f"turn_sigma {ts:.1f}, under {rules['min_turn_sigma']:g} (the turn is within its own noise)")
    try:
        minor = float(str(row.get("inlier_spread_px") or "").split("/")[0])
    except ValueError:
        minor = None
    if minor is None or minor < rules["min_spread_px"]:
        why.append("inlier spread not measured" if minor is None else
                   f"inliers bunched ({minor:.0f} px along the minor axis, under {rules['min_spread_px']:g})")
    d = [_num(row.get(k)) for k in ("d_azimuth_deg", "d_tilt_deg", "d_roll_deg")]
    if any(x is None for x in d):
        if f in ("ok", "moved", "MOVED"):
            why.append("no measured turn")
        return None, why
    turn = max(abs(x) for x in d)
    if turn > rules["max_turn_deg"]:
        why.append(f"a turn of {turn:.2f} deg, over {rules['max_turn_deg']:g} (a re-set? it may have moved too)")
    if "min_turn_deg" in rules and turn <= rules["min_turn_deg"]:
        why.append(f"a turn of {turn:.3f} deg, within the match's own error against another day "
                   f"({rules['min_turn_deg']:g} deg): at rest, as the default build grids it")
    if not why and "max_day_spread_deg" in rules:
        why = day_held(row, d, rules)
    return (None if why else d), why


def frame_gates(rules):
    """The gates of one frame's reading (correction_rules without the gates of the day)."""
    return {k: rules[k] for k in ("min_inliers", "min_turn_sigma", "min_spread_px", "max_turn_deg")}


def day_held(row, d, rules):
    """
    Did the pointing a day's probe measured (d) hold through the day? row["day_check"]: the feature records
    (survey_products.feature_record, with 'time_utc') of the day's other frames, its first and last photos of the
    detection (survey_products.check_pointing; [] when the probe is the day's only photo). Each must pass the gates
    of a frame and read d within rules["max_day_spread_deg"] in pan, tilt and roll. -> [why not] ([]: it held).
    """
    chk = row.get("day_check")
    if chk is None:
        return ["the pointing through the day not checked (only the probe frame was matched)"]
    lim = rules["max_day_spread_deg"]
    why = []
    for fr in chk:
        dd, w = day_correction(fr, frame_gates(rules))
        t = f"the {fr.get('time_utc') or '?'} UTC frame"
        if dd is None:
            why.append(f"{t} gives no reading ({'; '.join(w)}): the pointing is not shown to hold through the day")
        elif max(abs(a - b) for a, b in zip(dd, d)) > lim:
            why.append(f"{t} reads {dd[0]:+.3f}/{dd[1]:+.3f}/{dd[2]:+.3f} deg, more than {lim:g} deg from the "
                       f"probe's {d[0]:+.3f}/{d[1]:+.3f}/{d[2]:+.3f}: the pointing did not hold through the day "
                       f"(turned or came back)")
    return why


def window_corrections(rows, rules, ref_day=None):
    """day_correction for each day of a window's pointing rows. -> ({day: {d_deg, inliers, turn_sigma,
    inlier_spread_px, verdict, and with the gates of the day: probe_frame, day_check [{frame, time_utc, d_deg}],
    day_spread_deg (the largest difference from the probe's turn)}} of the days that pass, {day: [why not]} of the
    others), by date."""
    corr, why_not = {}, {}
    for r in sorted(rows or [], key=lambda r_: str(r_.get("date") or "")):
        d, why = day_correction(r, rules, ref_day)
        if d is None:
            why_not[r.get("date")] = why
            continue
        corr[r["date"]] = {"d_deg": [round(float(x), 3) for x in d], "inliers": int(_num(r.get("inliers"))),
                           "turn_sigma": _num(r.get("turn_sigma")),
                           "inlier_spread_px": r.get("inlier_spread_px"), "verdict": r.get("verdict")}
        if "max_day_spread_deg" in rules:
            chk = [{"frame": f.get("feature_frame"), "time_utc": f.get("time_utc"),
                    "d_deg": [_num(f.get(k)) for k in ("d_azimuth_deg", "d_tilt_deg", "d_roll_deg")]}
                   for f in r.get("day_check") or []]
            corr[r["date"]].update(probe_frame=r.get("feature_frame"), probe_time_utc=r.get("feature_time_utc"),
                                   day_check=chk, day_spread_deg=round(max(
                                       [abs(a - b) for f in chk for a, b in zip(f["d_deg"], d)] or [0.0]), 3))
    return corr, why_not


# ---------------------------------------------------------------------
# 2. The corrected pointing and the waterlines mapped with it
# ---------------------------------------------------------------------

def corrected_eo(eo, d_deg):
    """The EO of a camera that turned by d_deg ([pan, tilt, roll], deg) from eo: the same position, the angles plus
    the change. Exact for pointing_check.measure's change, which is R's CIRN angles less the reference's
    (R_to_angles, the inverse of georectify.cirn_angles_to_R), R the day's fitted rotation."""
    out = np.array(eo, float).copy()
    out[3:6] = out[3:6] + np.radians(np.asarray(d_deg, float))
    return out


def write_day_eo(path, eo, notes):
    """An EO file in the station's format (key: value, degrees), with more decimals than a calibration's (its
    angles are the calibration's plus a change), and '#' note lines saying what it is. -> the EO as read back by
    georectify.load_extrinsics: what the rows are mapped with, so georectify.py with this file gives them again."""
    a = np.degrees(np.asarray(eo, float)[3:6])
    Path(path).write_text(f"x: {eo[0]:.6f}\ny: {eo[1]:.6f}\nz: {eo[2]:.6f}\n"
                          f"azimuth: {a[0]:.8f}\ntilt: {a[1]:.8f}\nroll: {a[2]:.8f}\n"
                          + "".join(f"# {n}\n" for n in notes))
    return load_extrinsics(path)


def row_day(r):
    """A waterline row's UTC day: from capture_epoch (the photo name's epoch, as the pointing step's days), else the
    first 10 characters of capture_time_utc."""
    try:
        return datetime.fromtimestamp(float(r["capture_epoch"]), tz=timezone.utc).strftime("%Y-%m-%d")
    except (KeyError, TypeError, ValueError, OverflowError):
        return str(r.get("capture_time_utc") or "")[:10]


def regeorectify(ground_csv, cam, io, eo_by_day):
    """
    Maps again, in place, the rows of `cam` on the days of eo_by_day ({day: EO}) of a georectified waterline file
    (georectify.py's output): each point on its own plane (beach_elevation_navd88, else tide_elevation_navd88) with
    that day's EO, by georectify.pixel_to_ground as georectify.py does (the same 1000 m range limit, 3 decimals,
    blank where the ray does not reach the plane). Every other row is written back as it was read, which for
    georectify.py's files (csv's own quoting and line ends) is the same bytes; with no row of those days the file is
    not touched. -> {day: {rows, georectified, moved_m_median, moved_m_max}}: how far the day's points moved from
    where the calibration had put them (points mapped both ways).
    """
    with open(ground_csv, newline="") as f:
        rd = csv.DictReader(f)
        fields = list(rd.fieldnames or [])
        rows = list(rd)
    by = {}
    for i, r in enumerate(rows):
        if r.get("camera") == cam:
            d = row_day(r)
            if d in eo_by_day:
                by.setdefault(d, []).append(i)
    stats = {}
    for d, idx in sorted(by.items()):
        g = [rows[i] for i in idx]
        U = np.array([float(r["pixel_column"]) for r in g])
        V = np.array([float(r["pixel_row"]) for r in g])
        Z = np.array([float(r.get("beach_elevation_navd88") or r["tide_elevation_navd88"]) for r in g])
        E, N = pixel_to_ground(U, V, Z, io, eo_by_day[d])
        moved, n_ok = [], 0
        for r, e, n in zip(g, E, N):
            e0, n0 = _num(r.get("easting_utm19")), _num(r.get("northing_utm19"))
            if np.isfinite(e) and np.isfinite(n):
                r["easting_utm19"], r["northing_utm19"] = f"{e:.3f}", f"{n:.3f}"
                n_ok += 1
                if e0 is not None and n0 is not None:
                    moved.append(float(np.hypot(e - e0, n - n0)))
            else:
                r["easting_utm19"], r["northing_utm19"] = "", ""
        stats[d] = {"rows": len(g), "georectified": n_ok,
                    "moved_m_median": round(float(np.median(moved)), 3) if moved else None,
                    "moved_m_max": round(float(np.max(moved)), 3) if moved else None}
    if by:
        tmp = str(ground_csv) + ".tmp"
        with open(tmp, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        os.replace(tmp, str(ground_csv))
    return stats


def apply_corrections(ground_csv, cam, io, eo, eo_name, corrections, eo_dir, reference=""):
    """
    survey_products.py --correct-pointing, for one camera: for each corrected day ({day: window_corrections'
    record}) the EO file of its corrected pointing in eo_dir (<cam>_<day>_pointing_corrected_EO.yaml; the folder is
    emptied first, so it holds this build's only) and its rows of ground_csv mapped again with it (regeorectify).
    -> {day: the record, with the file, its sha256, and the rows and how far they moved}.
    """
    eo_dir = Path(eo_dir)
    shutil.rmtree(eo_dir, ignore_errors=True)
    eo_dir.mkdir(parents=True)
    eos, out = {}, {}
    for day, c in sorted(corrections.items()):
        a, t, r = c["d_deg"]
        p = eo_dir / f"{cam}_{day}_pointing_corrected_EO.yaml"
        eos[day] = write_day_eo(p, corrected_eo(eo, c["d_deg"]), [
            f"Made by pointing_daily.py (survey_products.py --correct-pointing): the pointing of {cam} on {day} = "
            f"{eo_name} with the turn matched features measured that day from the reference frame "
            f"{reference or '(the pointing step)'} added: {a:+.3f}/{t:+.3f}/{r:+.3f} deg pan/tilt/roll "
            f"({c['inliers']} inliers, turn_sigma {c['turn_sigma']}); position unchanged.",
            "Not a calibration: valid for that day's photos only."])
        out[day] = dict(c, eo_file=p.name, eo_sha256=hashlib.sha256(p.read_bytes()).hexdigest()[:16])
    st = regeorectify(ground_csv, cam, io, eos)
    for day in out:
        out[day].update(st.get(day) or {"rows": 0, "georectified": 0, "moved_m_median": None, "moved_m_max": None})
    return out


# ---------------------------------------------------------------------
# A camera that only turned, for the checks (no scipy)
# ---------------------------------------------------------------------

def rotated_view(img, io, eo_src, eo_dst, step=4):
    """The exact image a camera that only turned, from eo_src to eo_dst, would take of what img (taken at eo_src)
    shows: each pixel's ray rotated into the source frame and projected with the lens model (view_reproject's
    distort_uv, georectify's undistort_uv); OpenCV only, as the NUC has no scipy (view_reproject.reproject_image,
    which pointing_check.py --self-test uses, needs it)."""
    import cv2
    import pointing_check as pc
    from georectify import cirn_angles_to_R, undistort_uv
    from view_reproject import distort_uv
    nu, nv = int(io[0]), int(io[1])
    UU, VV = np.meshgrid(np.arange(0, nu + step, step, dtype=float), np.arange(0, nv + step, step, dtype=float))
    uu, vv = undistort_uv(UU.ravel(), VV.ravel(), io)
    r1 = cirn_angles_to_R(*eo_src[3:6]) @ (cirn_angles_to_R(*eo_dst[3:6]).T @ (pc.K_inv(io) @ np.vstack(
        [uu, vv, np.ones_like(uu)])))
    q = np.linalg.inv(pc.K_inv(io)) @ r1
    su, sv = distort_uv(q[0] / q[2], q[1] / q[2], io)
    mx, my = np.meshgrid(np.arange(nu, dtype=np.float32) / step, np.arange(nv, dtype=np.float32) / step)
    mapx = cv2.remap(su.reshape(UU.shape).astype(np.float32), mx, my, cv2.INTER_LINEAR)
    mapy = cv2.remap(sv.reshape(UU.shape).astype(np.float32), mx, my, cv2.INTER_LINEAR)
    return cv2.remap(img, mapx, mapy, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)


def textured_scene(io, contrast=40.0, seed=1):
    """A frame-sized grey texture (blurred noise): features everywhere, for the synthetic checks."""
    import cv2
    rng = np.random.default_rng(seed)
    t = cv2.GaussianBlur(rng.normal(0, 1, (int(io[1]), int(io[0]))).astype(np.float32), (0, 0), 3.0)
    t = 128 + contrast * t / t.std()
    return np.clip(np.repeat(t[:, :, None], 3, axis=2), 0, 255).astype(np.uint8)


def beach_points(io, eo, z=0.5, n=400, seed=3, max_range=400.0):
    """Ground points on the plane z that the camera at eo sees up to max_range m away (the waterlines lie at ~100-350
    m): (E, N), for the checks."""
    rng = np.random.default_rng(seed)
    u = rng.uniform(0.05, 0.95, 4 * n) * io[0]
    v = rng.uniform(0.08, 0.6, 4 * n) * io[1]
    E, N = pixel_to_ground(u, v, z, io, eo)
    ok = np.isfinite(E) & np.isfinite(N)
    ok[ok] = np.hypot(E[ok] - eo[0], N[ok] - eo[1]) <= max_range
    return E[ok][:n], N[ok][:n]


def self_test(image=None, quiet=False):
    """
    Planted turns recovered end to end: a synthetic scene (or a real c2 frame, --image) seen by c2 at its live
    calibration is turned by known angles (rotated_view: exact for a camera that only turns); pointing_check.measure
    (second=True) against the unturned frame, survey_products.feature_record, the gates (day_correction), the
    corrected EO (corrected_eo) and waterline points mapped again (regeorectify) must give the planted turn back
    to 0.02 deg and the points where they lie on the ground to 0.5 m; an unturned frame, whose match is within its
    own noise, must correct nothing and leave its rows as they were. -> 0 (PASS) or 1.
    """
    import cv2
    import tempfile
    import pointing_check as pc
    import survey_products as sp
    from view_reproject import ground_to_pixel
    say = (lambda *a: None) if quiet else print
    io = load_intrinsics(CAL / "CACO05_c2_20240801_IO.yaml")
    eo = load_extrinsics(LIVE_EO["c2"])
    img = cv2.imread(str(image)) if image else textured_scene(io)
    if img is None or img.shape[:2] != (int(io[1]), int(io[0])):
        say(f"pointing_daily.py --self-test: {image} is not a {int(io[0])} x {int(io[1])} c2 frame")
        return 1
    rules = station_rules()
    tmp = Path(tempfile.mkdtemp(prefix="pointing_daily_selftest_"))
    ok = True
    try:
        ref = tmp / "ref.jpg"
        cv2.imwrite(str(ref), img)
        ref.with_suffix(".json").write_text(json.dumps({"eo": eo.tolist()}))
        refs = [pc.Reference(ref)]
        E0, N0 = beach_points(io, eo)
        say(f"pointing_daily.py self-test ({'real frame ' + Path(image).name if image else 'synthetic scene'}, c2 "
            f"at {LIVE_EO['c2'].name}; {len(E0)} waterline points at z = +0.5 m, "
            f"{np.hypot(E0 - eo[0], N0 - eo[1]).min():.0f}-{np.hypot(E0 - eo[0], N0 - eo[1]).max():.0f} m away)")
        for turn in ((0.0, 0.0, 0.0), (0.30, -0.20, 0.10), (-1.0, 0.5, -0.3), (2.0, 0.0, 0.0), (0.05, 0.0, 0.0)):
            eo2 = corrected_eo(eo, turn)
            m = pc.measure(rotated_view(img, io, eo, eo2), refs, io, min_inliers=rules["min_inliers"], second=True)
            row = sp.feature_record(m, Path("probe.jpg"))
            d, why = day_correction(row, rules)
            # the points the turned camera sees, as its detector would report them; georectified with the
            # calibration (what the default products do), then mapped again
            U, V, vis = ground_to_pixel(E0, N0, np.full_like(E0, 0.5), io, eo2)
            ec, nc = pixel_to_ground(U[vis], V[vis], 0.5, io, eo)
            g = tmp / "ground.csv"
            with open(g, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["source_file", "camera", "capture_time_utc", "capture_epoch", "pixel_column", "pixel_row",
                            "tide_elevation_navd88", "beach_elevation_navd88", "easting_utm19", "northing_utm19"])
                for u_, v_, e_, n_ in zip(U[vis], V[vis], ec, nc):
                    w.writerow(["f", "c2", "2026-10-07T16:00:00Z", 1791388800, f"{u_:.2f}", f"{v_:.2f}", "0.5", "0.5",
                                f"{e_:.3f}" if np.isfinite(e_) else "", f"{n_:.3f}" if np.isfinite(n_) else ""])
            before = g.read_bytes()
            err0 = np.hypot(ec - E0[vis], nc - N0[vis])
            if d is None:
                good = max(abs(x) for x in turn) < 0.1 and before == g.read_bytes()
                say(f"  turn {turn[0]:+.2f}/{turn[1]:+.2f}/{turn[2]:+.2f} deg: not corrected ({'; '.join(why)}); "
                    f"{row.get('inliers')} inliers; the uncorrected lines lie {np.nanmedian(err0):.2f} m (median) off "
                    f"-> {'ok' if good else 'FAIL (a planted turn was not corrected)'}")
                ok &= good
                continue
            regeorectify(g, "c2", io, {"2026-10-07": corrected_eo(eo, d)})
            back = list(csv.DictReader(open(g, newline="")))
            e1 = np.array([_num(r["easting_utm19"]) if r["easting_utm19"] else np.nan for r in back])
            n1 = np.array([_num(r["northing_utm19"]) if r["northing_utm19"] else np.nan for r in back])
            err1 = np.hypot(e1 - E0[vis], n1 - N0[vis])
            derr = max(abs(a - b) for a, b in zip(d, turn))
            good = derr < 0.02 and np.nanmax(err1) < 0.5 and np.isfinite(err1).mean() > 0.95
            say(f"  turn {turn[0]:+.2f}/{turn[1]:+.2f}/{turn[2]:+.2f} deg -> measured {d[0]:+.3f}/{d[1]:+.3f}/"
                f"{d[2]:+.3f} ({row.get('inliers')} inliers, turn_sigma {row.get('turn_sigma')}; error {derr:.3f} "
                f"deg); waterline points off by {np.nanmedian(err0):.2f} m median ({np.nanmax(err0):.2f} max) with "
                f"the calibration, {np.nanmedian(err1):.3f} m ({np.nanmax(err1):.3f} max) corrected "
                f"-> {'ok' if good else 'FAIL'}")
            ok &= bool(good)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    say("PASS" if ok else "FAIL")
    return 0 if ok else 1


# ---------------------------------------------------------------------
# 3. The live log
# ---------------------------------------------------------------------

def eo_day(eo_file):
    """The day in a calibration's file name (CACO05_c1_20251113_EO-CV.yaml -> 2025-11-13), or None."""
    m = re.search(r"_(\d{8})_", Path(str(eo_file)).name)
    return f"{m.group(1)[:4]}-{m.group(1)[4:6]}-{m.group(1)[6:]}" if m else None


def eo_station(eo_file):
    """The station ID a calibration's file name starts with (CACO05_...), or None."""
    m = re.match(r"^([A-Z]+\d+)_", Path(str(eo_file)).name)
    return m.group(1) if m else None


def archive_frames(image_dir, cam, station=None):
    """{UTC day: [(epoch, Path), ...] by time} of a folder's <cam> timex photos (the cron's archive/images_timex,
    flat), of `station` only when given (the station ID in the name, detect_original_view.station_of; a name
    without one is kept)."""
    from detect_original_view import station_of
    out = {}
    for p in Path(image_dir).glob(f"*.{cam}.timex.jpg"):
        m = FRAME_EPOCH.match(p.name)
        if not m:
            continue
        st = station_of(p.name)
        if station and st and st != station:
            continue
        e = int(m.group(1))
        out.setdefault(datetime.fromtimestamp(e, tz=timezone.utc).strftime("%Y-%m-%d"), []).append((e, p))
    for v in out.values():
        v.sort()
    return out


def spread_frames(frames, n, hours):
    """Up to n of a day's frames inside the UTC hours (h0, h1), spread evenly over the day (first, middle, last)."""
    h0, h1 = hours
    inside = [(e, p) for e, p in frames
              if h0 <= datetime.fromtimestamp(e, tz=timezone.utc).hour
              + datetime.fromtimestamp(e, tz=timezone.utc).minute / 60.0 <= h1]
    if len(inside) <= n:
        return inside
    return [inside[i] for i in sorted(set(np.linspace(0, len(inside) - 1, n).round().astype(int).tolist()))]


def same_eo(meta_path, eo):
    """Was the reference whose .json this is taken at this EO?"""
    try:
        meta = json.loads(Path(meta_path).read_text())
        return bool(np.allclose(np.asarray(meta.get("eo"), float), eo, rtol=0, atol=1e-9))
    except (OSError, ValueError, TypeError):
        return False


def live_references(cam, io, eo, eo_file, frames, args):
    """The reference frame(s) of a camera's log (the order of the module's docstring), made when needed. ->
    (paths, info {kind, frames, text}) or ([], info with 'why' and 'how')."""
    bank = Path(args.bank_dir) / cam
    refs = [p for p in sorted(bank.glob("*.jpg")) if not p.stem.endswith("_mask") and same_eo(p.with_suffix(".json"), eo)]
    if refs:
        return refs, {"kind": "bank", "frames": [p.name for p in refs],
                      "text": f"pointing_check.py's reference bank {bank}, taken at {eo_file}: the changes are against "
                              f"the calibration's pointing as those frames show it"}
    own = Path(args.log_dir) / "pointing_daily_reference" / cam
    have = [p for p in sorted(own.glob("*.jpg")) if not p.stem.endswith("_mask") and same_eo(p.with_suffix(".json"), eo)]
    if have:
        info = json.loads(have[0].with_suffix(".json").read_text()).get("info") or {"kind": "own"}
        if args.reference_day and args.reference_day != info.get("day"):
            print(f"pointing {cam}: --reference-day {args.reference_day} not used: the log keeps its reference "
                  f"{have[0].name} ({info.get('day') or '?'}) at this calibration; to replace it, remove {own} and "
                  f"run again with --reference-day")
        return have[:1], info
    day = args.reference_day or eo_day(eo_file)
    kind = "reference day" if args.reference_day else "calibration day"
    cands = (frames.get(day) or []) if day else []
    if not cands:
        what = "--reference-day" if args.reference_day else f"the calibration day of {eo_file}"
        return [], {"kind": None,
                    "why": (f"no reference: no {cam} frame of {day} ({what}) in {args.image_dir}, and no frame in "
                            f"{bank} taken at {eo_file}"),
                    "how": (f"python3 pointing_daily.py --camera {cam} --reference-day <a clear day at the "
                            f"calibration's pointing>, or python3 pointing_check.py --camera {cam} --add-reference "
                            f"<image.jpg> --eo {eo_file} (the hourly monitor in owg.sh then uses it too)")}
    import cv2
    from estimate_eo_rotation import land_mask
    e, src = min(cands, key=lambda ep: abs(datetime.fromtimestamp(ep[0], tz=timezone.utc).hour
                                           + datetime.fromtimestamp(ep[0], tz=timezone.utc).minute / 60.0 - 17))
    img = cv2.imread(str(src))
    if img is None:
        return [], {"kind": None, "why": f"no reference: {src} cannot be read", "how": "give another --reference-day"}
    own.mkdir(parents=True, exist_ok=True)
    dst = own / src.name
    shutil.copy2(str(src), str(dst))
    info = {"kind": kind, "frames": [src.name], "day": day,
            "text": (f"the frame nearest 17:00 UTC of the calibration's own day {day}: the changes are against the "
                     f"calibration's pointing" if kind == "calibration day" else
                     f"the frame nearest 17:00 UTC of {day} (--reference-day), taken as at {eo_file}: the changes are "
                     f"against that day's pointing, not the calibration's"),
            "made_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
    dst.with_suffix(".json").write_text(json.dumps({"eo": np.asarray(eo, float).tolist(), "eo_file": str(eo_file),
                                                    "info": info}, indent=1) + "\n")
    cv2.imwrite(str(dst.with_name(dst.stem + "_mask.png")), land_mask(img.shape[:2], io, eo))
    return [dst], info


def measure_day(cam, day, frames, refs, ref_names, io, eo, rules, args):
    """One camera-day of the log: each chosen frame matched (pointing_check.measure, second=True), its feature record
    (survey_products.feature_record) and its gates (day_correction); the day's median turn over the frames that pass,
    the spread between them and the EO a correction would use. -> the day's entry."""
    import cv2
    import pointing_check as pc
    import survey_products as sp
    hours = tuple(float(x) for x in args.utc_hours.split("-"))
    chosen = spread_frames(frames, args.frames, hours)
    ent = {"checked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "frames_in_day": len(frames),
           "frames": []}
    for e, p in chosen:
        fr = {"frame": p.name, "time_utc": datetime.fromtimestamp(e, tz=timezone.utc).strftime("%H:%M")}
        if p.name in ref_names:
            fr.update(passes=False, why_not=["the reference frame itself"])
            ent["frames"].append(fr)
            continue
        img = cv2.imread(str(p))
        m = pc.measure(img, refs, io, min_inliers=rules["min_inliers"], second=True) if img is not None else None
        rec = sp.feature_record(m, p)
        d, why = day_correction(rec, rules)
        fr.update({k: rec.get(k) for k in ("features", "inliers", "turn_sigma", "inlier_spread_px",
                                           "runner_up_inliers", "runner_up_d")})
        if m is not None and m.get("reference"):
            fr["reference"] = m["reference"]
        if rec.get("d_azimuth_deg") not in ("", None):
            fr["d_deg"] = [rec["d_azimuth_deg"], rec["d_tilt_deg"], rec["d_roll_deg"]]
        fr.update(passes=d is not None, why_not=why)
        ent["frames"].append(fr)
    ok = [f["d_deg"] for f in ent["frames"] if f.get("passes")]
    ent["frames_passing"] = len(ok)
    if not chosen:
        ent.update(status="no frames", d_deg=None)
        return ent
    if not ok:
        why = sorted({w for f in ent["frames"] for w in f.get("why_not") or []})
        ent.update(status="no reading", d_deg=None, why_not=why)
        return ent
    a = np.array(ok, float)
    d = [round(float(x), 3) for x in np.median(a, axis=0)]
    turn = max(abs(x) for x in d)
    ent.update(d_deg=d, range_deg=[round(float(x), 3) for x in a.max(axis=0) - a.min(axis=0)],
               status="MOVED" if turn > sp.POINT_FAIL_DEG else "moved" if turn > sp.POINT_WARN_DEG else "ok",
               would_correct_eo_deg=[round(float(x), 4) for x in np.degrees(corrected_eo(eo, d)[3:6])])
    return ent


def load_log(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return None


def save_log(path, log):
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(log, indent=1, sort_keys=False) + "\n")
    os.replace(str(tmp), str(path))


def day_line(cam, day, ent, eo_name):
    """The cron log's line for one camera-day (never the word WARNING: station_status.py counts those)."""
    n = ent.get("frames_in_day", 0)
    if ent.get("status") == "no frames":
        return (f"pointing {cam} {day}: no frames" + (f" ({n} in the archive that day, none in the UTC hours)" if n else
                                                     " in the archive") + "; log only")
    if ent.get("status") == "no reading":
        return f"pointing {cam} {day}: no reading ({n} frame(s) that day; {'; '.join(ent.get('why_not') or [])}); log only"
    d, rg = ent["d_deg"], ent.get("range_deg") or [0, 0, 0]
    return (f"pointing {cam} {day}: {ent['status']} {d[0]:+.3f}/{d[1]:+.3f}/{d[2]:+.3f} deg pan/tilt/roll from the "
            f"reference ({ent['frames_passing']} of {len(ent['frames'])} frame(s) pass the gates"
            + (f", spread {max(rg):.3f} deg between them" if ent["frames_passing"] > 1 else "")
            + f"); log only: the live products keep {eo_name}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--camera", nargs="+", default=["c1", "c2"], choices=["c1", "c2"])
    ap.add_argument("--image-dir", default=str(ARCHIVE / "images_timex"),
                    help="the archived timex photos (default the cron's archive/images_timex)")
    ap.add_argument("--log-dir", default=str(ARCHIVE), help="where pointing_daily_<cam>.json goes (default archive/)")
    ap.add_argument("--bank-dir", default=str(CAL / "pointing"),
                    help="pointing_check.py's reference banks, <dir>/<cam>/ (read only; default calibration/pointing)")
    ap.add_argument("--eo-c1", default=str(LIVE_EO["c1"]), help="c1's live calibration (default the cron's)")
    ap.add_argument("--eo-c2", default=str(LIVE_EO["c2"]), help="c2's live calibration (default the cron's)")
    ap.add_argument("--days", nargs=2, metavar=("FIRST", "LAST"), default=None,
                    help="measure these UTC days (again, if logged); default: the complete days not yet logged")
    ap.add_argument("--lookback-days", type=int, default=7, help="how far back new days are looked for (default 7)")
    ap.add_argument("--max-days", type=int, default=2,
                    help="days measured per camera and run, newest first (default 2: a backlog fills over a few runs)")
    ap.add_argument("--frames", type=int, default=3, help="frames per day, spread over it (default 3)")
    ap.add_argument("--utc-hours", default="13-20", help="the frames' UTC hours (default 13-20: daylight all year)")
    ap.add_argument("--reference-day", default=None,
                    help="make the log's own reference from this day's frame nearest 17:00 UTC (when it has none)")
    ap.add_argument("--self-test", action="store_true", help="planted turns recovered (no scipy needed)")
    ap.add_argument("--image", default=None, help="with --self-test: plant the turns on this real c2 frame")
    args = ap.parse_args()
    if args.self_test:
        return self_test(args.image)
    try:
        h0, h1 = (float(x) for x in args.utc_hours.split("-"))
    except ValueError:
        ap.error(f"--utc-hours {args.utc_hours!r}: expected e.g. 13-20")
    rules = station_rules()
    today = datetime.now(timezone.utc).date()
    for cam in args.camera:
        eo_path = Path(getattr(args, f"eo_{cam}"))
        io = load_intrinsics(CAL / f"CACO05_{cam}_20240801_IO.yaml")
        eo = load_extrinsics(eo_path)
        log_path = Path(args.log_dir) / f"pointing_daily_{cam}.json"
        log = load_log(log_path) or {}
        frames = archive_frames(args.image_dir, cam, eo_station(eo_path))
        refs, rinfo = live_references(cam, io, eo, eo_path.name, frames, args)
        log.update(camera=cam, what=(
            "pointing_daily.py: each complete UTC day's turn of the camera from its reference frame, by matched land "
            "features on up to a few archived timex spread over the day, gated as survey_products.py --correct-pointing "
            "gates a day (the feature gates only); d_deg is pan/tilt/roll (deg) from the reference, the median of the "
            "frames that pass; would_correct_eo_deg the EO angles a correction would use. LOG ONLY: the live products "
            "are georectified with the calibration below as it is"),
            calibration={"eo_file": eo_path.name, "eo_sha256": hashlib.sha256(eo_path.read_bytes()).hexdigest()[:16]},
            reference=rinfo, rules=rules, applied_live=False)
        log.setdefault("days", {})
        if not refs:
            print(f"pointing {cam}: {rinfo['why']}: nothing measured. To start: {rinfo['how']}")
            Path(args.log_dir).mkdir(parents=True, exist_ok=True)
            save_log(log_path, log)
            continue
        if args.days:
            a, b = (datetime.strptime(x, "%Y-%m-%d").date() for x in args.days)
            todo = [(a + timedelta(days=k)).isoformat() for k in range((b - a).days + 1)]
        else:
            new = [(today - timedelta(days=k)).isoformat() for k in range(1, max(args.lookback_days, 1) + 1)]
            todo = sorted([d for d in new if d not in log["days"]], reverse=True)[:max(args.max_days, 0)]
        import pointing_check as pc
        loaded = [pc.Reference(p) for p in refs]
        names = set(p.name for p in refs)
        for day in sorted(todo):
            ent = measure_day(cam, day, frames.get(day) or [], loaded, names, io, eo, rules, args)
            # what this day was measured against (the header's calibration and reference are the latest run's)
            ent.update(calibration=dict(log["calibration"]),
                       reference={k: rinfo.get(k) for k in ("kind", "day", "frames")})
            log["days"][day] = ent
            log["days"] = dict(sorted(log["days"].items()))
            save_log(log_path, log)
            print(day_line(cam, day, ent, eo_path.name))
        if not todo:
            print(f"pointing {cam}: no new complete day to log (the last {args.lookback_days} are in {log_path})")
    return 0


if __name__ == "__main__":
    sys.exit(main())

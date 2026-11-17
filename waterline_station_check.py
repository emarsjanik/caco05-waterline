#!/usr/bin/env python3
"""
Checks For A New Station's Waterline Settings (coastcam-station-setup)
=======================================================================
    python3 waterline_station_check.py capabilities
    python3 waterline_station_check.py check-overlay --conf waterline_station.conf
    python3 waterline_station_check.py pair --xyz T.csv --uv UV.csv --io IO1 [IO2 ..] --eo EO1 [EO2 ..] [--image WxH]
    python3 waterline_station_check.py profile --camera c1 --io IO --eo EO --image-size WxH --output P.json
                                       (--survey-csv S.csv | --dem D.tif | --geometry [--plane-z MSL])
                                       [--mllw M --mhhw M] [--merge]
    python3 waterline_station_check.py trial --camera c1 --io IO --eo EO --profile-json P.json
                                       --spline SPLINE --frames T1 [.. T5] --antenna-navd88 M [--work-dir DIR]
    python3 waterline_station_check.py pointing --io IO --eo EO --calibration-picture IMG --timex T1 [T2 ..]

WHAT. coastcam-station-setup (./setup.sh --waterline) sets this repository up
on a station other than CACO-05. It asks this script what this checkout can
do (capabilities) and has it check what it is about to install. CACO-05 itself
never needs it: it has no waterline_station.conf.

  capabilities     prints the feature tags of this checkout, on one line:
                   overlay-1 storm-alerts-1 station-setup-1
  check-overlay    checks a waterline_station.conf before waterline_timex_cron.sh
                   sources it (the cron runs it on every run, and stops with a
                   FATAL line when it refuses). The file is shell, so only this
                   is accepted, one item a line:
                     * a comment (# ...) or a blank line;
                     * KEY="value": a name of ALLOWED_KEYS at the start of the
                       line, '=', the value in double quotes and nothing after
                       them. Inside the quotes anything but $, `, \\, " and
                       control characters (so ; ( ) = are fine: the shell
                       expands nothing there). Each key once;
                     * one line 'export NAME ...' naming only EXPORTED_KEYS.
                   Then each value is checked for its key (numbers, 0/1
                   switches, absolute paths without spaces, file names, the
                   cameras, the time zone, detector options that would move
                   the cron's own folders), the two cameras' calibration names
                   against CAMS, and, unless WATERLINE_STAGE="archive" (stage 1,
                   which runs no Python), that the file names every station
                   fact CACO-05's default would otherwise silently supply
                   (STAGE2_REQUIRED) and exports the names the Python scripts
                   read from the environment. A forbidden key (FORBIDDEN_KEYS:
                   the cron's own settings, which the storm watch and the
                   survey products read from the cron itself) or an unknown one
                   is refused.
  pair             (design 6.4) which lens (IO), position (EO), target and picks
                   files belong together, for ONE target-coordinate file (header-
                   less num,E,N,Z) and ONE picks file (num,U,V): every IO x EO goes
                   through georectify.self_test, unchanged (its PASS: mean under
                   1.0 m), and is reported with its mean, max and n. A set is
                   chosen when it passes and the next best POSITION is at least 3
                   times worse (the ratio test is for the choice of EO, so two
                   cameras with the same lens model are not refused): two lenses
                   that pass with it within 20% of each other are an "io_tie" for
                   the caller to decide (a name hint, or asking), never a refusal.
                   Then: the lens's NU x NV against --image, picks inside the
                   frame, targets 10-1000 m from the camera and at -5 m .. camera
                   z + 5 m. Fewer than 4 shared targets: not paired.
  profile          (design 6.5) the detector's search band, written as
                   detect_original_view.py writes it ({"<cam>": {"crop": [0,1,0,1],
                   "envelope": ..., "no_bias": true}}, plus "step" [4,18] and
                   "profile_scale" [0.4,1.0] when the far side is on the left), by
                   detect_original_view.envelope_from_survey, unchanged: from a beach
                   survey's points (a 0.5 m grid by matplotlib.tri) or a beach DEM
                   (.tif/.asc, read by compare_dem_survey.read_survey), band MLLW
                   - 0.5 to MHHW + 1.0 m (else -1.5 to +2.5) within 300 m; or, with
                   no survey, a level plane at the gauge's MSL 40-600 m from the
                   camera: PROVISIONAL. The detector's other CACO05 tunings
                   (band_contrast_floor, snap settings) stay; the JSON says so.
  trial            (design 6.6) the detector (--no-bias-correction --profile-json,
                   as the station file's DETECTOR_ARGS), extract_elevation_contours.py
                   and georectify.py on up to 5 timex of one camera, every folder a
                   scratch one (never input/, processed_timex/ or archive/; the
                   pictures are copied). A frame passes with a line in >= 50% of
                   the columns, >= 80% of its points placed, all 30-1000 m from the
                   camera. Exit 0 when one frame passes. The antenna height comes
                   from --antenna-navd88 (or GNSSR_ANTENNA_NAVD88_M already in the
                   environment), never CACO-05's 19.014 m: without either, usage
                   (exit 2). The JSON says which height was used.
  pointing         (correction 7) a picture taken on the calibration day against
                   recent timex of the same camera, by pointing_check's matcher:
                   the rotation in degrees, PASS under 0.3, WARN 0.3-1.0, FAIL over
                   1.0; UNMATCHED (nothing checked) when the pictures do not match.

OUTPUT. Subcommands that report results print JSON on stdout (setup reads it)
and a summary for people on stderr. Exit status: 0 ok, 1 refused or failed,
2 usage.

Python 3.8 (CACO-05 and the NUCs: Ubuntu 20.04/22.04). capabilities and
check-overlay use the standard library only (check-overlay also asks this
checkout's station_env.py, itself standard library, for its time-zone rule,
so it accepts exactly the zones the Python steps accept); pair, profile, trial
and pointing load numpy and this repository's modules when they run (profile
also matplotlib, pointing OpenCV), and never change them.
"""

import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

CAPABILITIES = "overlay-1 storm-alerts-1 station-setup-1"

# ---------------------------------------------------------------------
# check-overlay: the names a waterline_station.conf may set
# ---------------------------------------------------------------------

# Every name the station file may set (docs: design 7.1, plus RAS_SOURCE).
ALLOWED_KEYS = (
    # stage, cameras and calibration
    "WATERLINE_STAGE", "CAMS", "CAL", "CAL_IO_C1", "CAL_EO_C1", "CAL_IO_C2", "CAL_EO_C2",
    # water levels and waves
    "GNSSR_SPLINE", "TIDE_GAUGE", "GAUGE_CSV", "WAVE_BUOY", "WAVES_CSV", "DEM_MAX_HS", "SETUP_COEF",
    # switches
    "DETECTOR_ARGS", "OWG_ENABLE", "RUNUP_ENABLE", "POINTING_LOG_ENABLE", "ARCHIVE_PRODUCTS", "RAS_SOURCE",
    # settings the Python scripts read from the environment (exported)
    "GNSSR_ANTENNA_NAVD88_M", "STATION_TZ", "GEOTIFF_EPSG", "STATION_NAME", "STATION_PROVISIONAL",
    "STATION_WAVE_CURRENCY", "GNSSR_QC_GROSS_M",
)

# The names the one export line may list: what the Python scripts read from the environment.
EXPORTED_KEYS = ("GNSSR_ANTENNA_NAVD88_M", "STATION_TZ", "GEOTIFF_EPSG", "STATION_NAME", "STATION_PROVISIONAL",
                 "STATION_WAVE_CURRENCY", "GNSSR_QC_GROSS_M", "WAVE_BUOY")

# Names that are refused with a reason (any other unknown name is refused too). The cron's
# own settings stay the cron's: the storm watch, survey products and timing recalibration read
# them from waterline_timex_cron.sh itself, so a station value would split the station in two.
FORBIDDEN_KEYS = {
    "TIMING_MODE": "the timing mode is the cron's one switch (recalibrate_timing.py and the station status read it)",
    "DEM_ENABLE": "the DEM settings are the cron's own (the storm watch builds its DEMs with them)",
    "DEM_CELL": "the DEM settings are the cron's own (the storm watch builds its DEMs with them)",
    "DEM_MIN_POINTS": "the DEM settings are the cron's own (the storm watch builds its DEMs with them)",
    "DEM_MAX_SPREAD": "the DEM settings are the cron's own (the storm watch builds its DEMs with them)",
    "DEM_MAX_DAY_OFFSET": "the DEM settings are the cron's own (the storm watch builds its DEMs with them)",
    "DEM_WINDOW_DAYS": "the DEM settings are the cron's own (the storm watch builds its DEMs with them)",
    "GROUND": "a file of the cron's own (the storm watch reads it)",
    "LOCK": "the cron's lock (the storm watch takes the same one)",
    "GNSSR_START_FILE": "a file of the cron's own (the storm watch reads it)",
    "BASE": "the folder the cron runs in",
    "LOG": "the cron's log",
    "SCRATCH": "a scratch folder the detector CLEARS on every run",
    "SCRATCH_INPUT": "a scratch folder the detector CLEARS on every run",
    "DEBUG_DIR": "a folder of the cron's own",
    "ARCHIVE_CSV": "the cron's archive",
    "ARCHIVE_IMG": "the cron's archive",
    "ARCHIVE_RAS": "the cron's archive",
    "PATH": "the shell's search path",
    "IFS": "the shell's word splitting",
    "PYTHONPATH": "the Python path",
    "PYTHONHOME": "the Python installation",
    "LD_PRELOAD": "the dynamic linker",
    "LD_LIBRARY_PATH": "the dynamic linker",
    "BASH_ENV": "a file bash runs",
    "ENV": "a file the shell runs",
    "SHELLOPTS": "the shell's options",
}

# A stage-2 file (WATERLINE_STAGE not "archive") must name these: without them the cron would use
# CACO-05's value without a word (its calibration, GNSS-R file, gauge, buoy, antenna height, setup
# coefficient, time zone, owg.sh with CACO-05's recipients ...). These are the names of setup's own
# full file (design 7.1), except STATION_PROVISIONAL, ARCHIVE_PRODUCTS and RAS_SOURCE.
STAGE2_REQUIRED = (
    "CAMS", "CAL", "CAL_IO_C1", "CAL_EO_C1", "CAL_IO_C2", "CAL_EO_C2",
    "GNSSR_SPLINE", "TIDE_GAUGE", "GAUGE_CSV", "WAVE_BUOY", "WAVES_CSV", "DEM_MAX_HS", "SETUP_COEF",
    "DETECTOR_ARGS", "OWG_ENABLE", "RUNUP_ENABLE", "POINTING_LOG_ENABLE",
    "GNSSR_ANTENNA_NAVD88_M", "STATION_TZ", "GEOTIFF_EPSG", "STATION_NAME",
)

# Detector options the cron sets itself: given again in DETECTOR_ARGS, argparse would take the
# later one (an --output-dir elsewhere would be CLEARED by the detector), and --help would print
# the help and detect nothing. argparse also takes any unambiguous prefix (--out), so prefixes
# are refused too.
DETECTOR_OWN_OPTIONS = ("--image-suffix", "--timing-mode", "--input-dir", "--output-dir", "--debug-dir",
                        "--source-dir", "--help")

KEY_LINE = re.compile(r'^([A-Za-z_][A-Za-z0-9_]*)="(.*)"$')
EXPORT_LINE = re.compile(r"^export(?:[ \t]+[A-Za-z_][A-Za-z0-9_]*)+[ \t]*$")
BAD_IN_QUOTES = re.compile(r'[$`\\"\x00-\x1f\x7f]')


def _float(v):
    """float(v) when v is a finite decimal number (no 'nan', 'inf', '1e', ' 1'), else None."""
    if not re.match(r"^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$", v):
        return None
    x = float(v)
    return x if math.isfinite(x) else None


def _is_abs_path(v):
    """An absolute path the cron can use unquoted: no spaces, no glob characters, no '..' parts."""
    return (v.startswith("/") and not re.search(r"[\s*?\[\]~]", v)
            and ".." not in v.split("/"))


def _known_zone(name):
    """True/False when the zone is known to the Python steps here, None when nothing can look it up.
    station_env's own rule (zoneinfo, else pytz, else dateutil), so this check and the contour and
    DEM steps the cron runs with the same python3 cannot disagree: on Python 3.8 dateutil alone
    would accept POSIX strings (EST5, UTC+3) and files such as 'localtime' that pytz, which
    station_env asks first there, does not know. station_env imports only the standard library."""
    _repo_modules()
    import station_env
    return station_env._zone_known(name)


def _shell_words(v):
    """The words the cron's unquoted $VAR gives: bash splits on its default IFS (space, tab, newline)
    only, never on Unicode spaces (a pasted no-break space keeps 'c1\u00a0c2' one word). Tab and
    newline cannot occur in an accepted value (BAD_IN_QUOTES), so splitting on ' ' is the same."""
    return [w for w in v.split(" ") if w]


def check_value(key, v):
    """The problem with value v of `key`, or None. Only the key's own format (see check_overlay)."""
    if key == "WATERLINE_STAGE":
        return None if v in ("", "archive") else 'must be "archive" (stage 1) or ""'
    if key == "CAMS":
        cams = _shell_words(v)
        if cams not in (["c1"], ["c1", "c2"]):
            return ('must be "c1" or "c1 c2" (the cron maps c1 and, on a two-camera station, c2)'
                    + ("; it holds a character that is not plain ASCII (a pasted no-break space? the cron "
                       "splits on plain spaces only)" if re.search(r"[^\x00-\x7f]", v) else ""))
        return None
    if key in ("CAL", "GNSSR_SPLINE", "GAUGE_CSV", "WAVES_CSV"):
        return None if _is_abs_path(v) else "must be an absolute path without spaces, wildcards or '..'"
    if key == "RAS_SOURCE":
        return None if v == "" or _is_abs_path(v) else \
            "must be \"\" (no timestacks kept) or an absolute path without spaces, wildcards or '..'"
    if key in ("CAL_IO_C1", "CAL_EO_C1", "CAL_IO_C2", "CAL_EO_C2"):
        if v == "" and key.endswith("C2"):
            return None
        if not re.match(r"^[A-Za-z0-9_][A-Za-z0-9_.+-]*$", v):
            return "must be a file name in CAL (letters, digits, _ . + -; no folder)" + \
                   (', or "" for a one-camera station' if key.endswith("C2") else "")
        return None
    if key == "TIDE_GAUGE":
        return None if re.match(r"^\d{7}$", v) else "must be a 7-digit NOAA station number"
    if key == "WAVE_BUOY":
        return None if re.match(r"^[A-Za-z0-9]{4,7}$", v) else "must be an NDBC station id (e.g. 44008)"
    if key == "DEM_MAX_HS":
        x = _float(v)
        return None if v == "" or (x is not None and x > 0) else \
            'must be a wave height in m above 0, or "" (keep every frame)'
    if key == "SETUP_COEF":
        x = _float(v)
        return None if v == "" or (x is not None and 0 <= x < 1) else \
            'must be a setup coefficient C from 0 to below 1, or "" (no wave setup)'
    if key in ("OWG_ENABLE", "RUNUP_ENABLE", "POINTING_LOG_ENABLE"):
        return None if v in ("0", "1") else 'must be "0" or "1"'
    if key == "ARCHIVE_PRODUCTS":
        return None if re.match(r"^([a-z]+( [a-z]+)*)?$", v) else \
            'must be product names in lower case separated by spaces (e.g. "bright dark"), or ""'
    if key == "DETECTOR_ARGS":
        if re.search(r"[*?\[\]~]", v):
            return "must not hold wildcards or '~' (the cron splits it into words unquoted)"
        for tok in v.split():
            if tok == "-h":
                return "must not hold -h"
            if tok.startswith("--"):
                opt = tok.split("=", 1)[0]
                own = [o for o in DETECTOR_OWN_OPTIONS if len(opt) > 2 and o.startswith(opt)]
                if own:
                    return f"must not hold {opt} ({own[0]} is the cron's own option)"
        return None
    if key == "GNSSR_ANTENNA_NAVD88_M":
        x = _float(v)
        return None if x is not None and -100 < x < 5000 else \
            "must be the antenna height in m NAVD88 (a number)"
    if key == "STATION_TZ":
        if not re.match(r"^[A-Za-z][A-Za-z0-9_+-]*(/[A-Za-z0-9_+-]+)*$", v):
            return "must be an IANA time zone name (e.g. America/New_York)"
        return None if _known_zone(v) is not False else "is not a time zone known here"
    if key == "GEOTIFF_EPSG":
        return None if re.match(r"^\d{4,6}$", v) else "must be an EPSG code (e.g. 6348)"
    if key == "GNSSR_QC_GROSS_M":
        if v == "":
            return None
        parts = v.split(",")
        lo, hi = (_float(parts[0]), _float(parts[1])) if len(parts) == 2 else (None, None)
        return None if lo is not None and hi is not None and lo < hi else \
            'must be "lo,hi" in m NAVD88 with lo below hi, or ""'
    if key == "STATION_WAVE_CURRENCY":
        return None if v == "" or re.match(r"^[A-Za-z0-9_]+:[A-Za-z0-9_]+/[A-Za-z0-9_]+$", v) else \
            'must be "<kind>:<height>/<period>" (e.g. buoy:hs/tp), or ""'
    if key == "STATION_NAME":
        return None if v.strip() else "must not be empty"
    return None                                          # STATION_PROVISIONAL: any text


def check_overlay(text):
    """
    The findings for the text of a waterline_station.conf:
    -> dict(ok, stage, keys {name: value}, exported [names], errors [{line, key, message}], notes [text]).
    """
    errors, notes = [], []
    keys, key_line, exported, export_lines = {}, {}, [], []

    def err(line, key, message):
        errors.append({"line": line, "key": key, "message": message})

    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()                                      # the newline that ends the last line
    for n, line in enumerate(lines, 1):
        if line.strip(" \t") == "":
            continue
        if line.lstrip(" \t").startswith("#"):
            continue
        if "\r" in line:
            err(n, None, "carriage return (a file saved with Windows line ends): save it with Unix ones")
            continue
        if line.startswith("export") and (len(line) == 6 or line[6] in " \t"):
            if not EXPORT_LINE.match(line):
                err(n, None, "an export line may only name variables: export NAME NAME ...")
                continue
            export_lines.append(n)
            for name in line.split()[1:]:
                if name not in EXPORTED_KEYS:
                    err(n, name, f"export {name}: only these may be exported: {' '.join(EXPORTED_KEYS)}")
                elif name not in exported:
                    exported.append(name)
            continue
        m = KEY_LINE.match(line)
        if not m:
            m2 = re.match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=", line)
            if m2 and not line.startswith(m2.group(1) + "="):
                err(n, m2.group(1), f'no spaces around "=", none before the name: {m2.group(1)}="value"')
            elif m2:
                err(n, m2.group(1), f'the value must be in double quotes with nothing after them: '
                                    f'{m2.group(1)}="value"')
            else:
                err(n, None, 'not a comment, a KEY="value" line or the export line')
            continue
        key, value = m.group(1), m.group(2)
        bad = BAD_IN_QUOTES.search(value)
        if bad:
            what = {"$": "$", "`": "a backtick", "\\": "a backslash", '"': 'a "'}.get(
                bad.group(0), "a control character")
            err(n, key, f"{key}: {what} inside the quotes (the shell would expand or end it)")
            continue
        if key in FORBIDDEN_KEYS:
            err(n, key, f"{key} may not be set here: {FORBIDDEN_KEYS[key]}")
            continue
        if key not in ALLOWED_KEYS:
            err(n, key, f"{key} is not a name this file may set")
            continue
        if key in keys:
            err(n, key, f"{key} is set twice (lines {key_line[key]} and {n})")
            continue
        keys[key], key_line[key] = value, n
        problem = check_value(key, value)
        if problem:
            err(n, key, f'{key}="{value}" {problem}')
        if key == "STATION_TZ" and not problem and _known_zone(value) is None:
            notes.append(f"STATION_TZ={value}: not looked up (no zoneinfo, pytz or dateutil here)")

    if len(export_lines) > 1:
        err(export_lines[1], None, f"more than one export line (lines {', '.join(map(str, export_lines))})")

    stage = "archive" if keys.get("WATERLINE_STAGE") == "archive" else "full"

    # the cameras and their calibration names
    if "CAMS" in keys and check_value("CAMS", keys["CAMS"]) is None:
        two = "c2" in _shell_words(keys["CAMS"])
        for k in ("CAL_IO_C2", "CAL_EO_C2"):
            if two and keys.get(k) == "":
                err(key_line[k], k, f'{k}="" but CAMS="{keys["CAMS"]}" maps c2: name its file')
            if not two and keys.get(k, None) != "":
                err(key_line.get(k, key_line["CAMS"]), k,
                    f'CAMS="{keys["CAMS"]}" (one camera): set {k}="" (else CACO-05\'s c2 file is looked for)')
    io2, eo2 = keys.get("CAL_IO_C2"), keys.get("CAL_EO_C2")
    if (io2 is None) != (eo2 is None) or (io2 == "") != (eo2 == ""):
        if not any(e["key"] in ("CAL_IO_C2", "CAL_EO_C2") for e in errors):
            err(key_line.get("CAL_IO_C2", key_line.get("CAL_EO_C2")), "CAL_IO_C2",
                "CAL_IO_C2 and CAL_EO_C2 go together: both named, or both \"\"")

    # stage 2: every station fact named, and exported where the Python scripts read it
    if stage == "full":
        missing = [k for k in STAGE2_REQUIRED if k not in keys]
        if missing:
            err(None, None, "missing (CACO-05's value would be used without a word): " + " ".join(missing)
                + '; a stage-1 file says WATERLINE_STAGE="archive"')
    unexported = [k for k in EXPORTED_KEYS if k in keys and k not in exported]
    if unexported:
        msg = ("set but not exported, so the Python scripts would not see them: " + " ".join(unexported)
               + " (add them to the export line)")
        if stage == "full":
            err(None, None, msg)
        else:
            notes.append(msg + "; stage 1 runs no Python")
    errors.sort(key=lambda e: (e["line"] is None, e["line"] or 0))
    return {"ok": not errors, "stage": stage, "keys": keys, "exported": exported, "errors": errors,
            "notes": notes}


# ---------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------

def say(text):
    """The summary for people: stderr (stdout is setup's JSON)."""
    print(text, file=sys.stderr, flush=True)


def emit(obj):
    print(json.dumps(obj, indent=1, sort_keys=True), flush=True)


def cmd_capabilities(args):
    print(CAPABILITIES, flush=True)
    return 0


def cmd_check_overlay(args):
    path = Path(args.conf)
    name = path.name
    try:
        data = path.read_bytes()
    except OSError as exc:
        emit({"command": "check-overlay", "conf": str(path), "ok": False,
              "errors": [{"line": None, "key": None, "message": f"cannot read: {exc.strerror or exc}"}]})
        say(f"check-overlay: {name} REFUSED: cannot read it ({exc.strerror or exc})")
        return 1
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        emit({"command": "check-overlay", "conf": str(path), "ok": False,
              "errors": [{"line": None, "key": None, "message": f"not UTF-8 text ({exc.reason})"}]})
        say(f"check-overlay: {name} REFUSED: not UTF-8 text")
        return 1
    res = check_overlay(text)
    out = {"command": "check-overlay", "conf": str(path), "ok": res["ok"], "stage": res["stage"],
           "keys": sorted(res["keys"]), "exported": res["exported"], "errors": res["errors"],
           "notes": res["notes"]}
    emit(out)
    if res["ok"]:
        say(f"check-overlay: {name} OK (stage {'1, archive only' if res['stage'] == 'archive' else 2}: "
            f"{len(res['keys'])} setting(s), {len(res['exported'])} exported)")
    else:
        say(f"check-overlay: {name} REFUSED ({len(res['errors'])} problem(s)):")
        for e in res["errors"]:
            say(f"  {'line %d' % e['line'] if e['line'] else 'file'}: {e['message']}")
    for n in res["notes"]:
        say(f"  note: {n}")
    return 0 if res["ok"] else 1


# ---------------------------------------------------------------------
# pair, profile, trial, pointing (design 6.4-6.6; corrections C5, C6, C8 and 7)
# ---------------------------------------------------------------------
# These import numpy and the repository's own modules only when they run, and use them unchanged
# (georectify.self_test and its readers; detect_original_view.envelope_from_survey, which uses
# view_reproject.ground_to_pixel; compare_dem_survey.read_survey; pointing_check's matcher), so
# capabilities and check-overlay, which the cron runs every time, stay standard library.

HERE = Path(__file__).resolve().parent

# pair (design 6.4; C5, C6)
PAIR_RATIO = 3.0              # the next best POSITION (EO) for the picks must be this many times worse
IO_TIE = 0.20                 # two lenses passing with the chosen EO within 20% of each other: an io_tie
PAIR_MIN_SHARED = 4           # targets the picks and the coordinates share (CHECKLIST: 4 or more)
PAIR_GOOD_SHARED = 6
UNPLACED_M = 1000.0           # in the ranking only: a target georectify cannot place counts as 1000 m off
TARGET_RANGE_M = (10.0, 1000.0)
TARGET_Z_LO_M = -5.0          # target Z from -5 m to camera z + 5 m NAVD88 (C6: a dune or bluff top is fine)
TARGET_Z_ABOVE_CAMERA_M = 5.0

# profile (design 6.5; C8)
SURVEY_BAND_M = (-1.5, 2.5)   # without the gauge's datums (detect_original_view's default --band)
DATUM_PAD_M = (0.5, 1.0)      # with them: MLLW - 0.5 .. MHHW + 1.0
SURVEY_RANGE_M = 300.0
SURVEY_CELL_M = 0.5
SURVEY_MAX_DIST_M = 1500.0    # survey points farther from the camera are ignored (design 6.2)
TRI_MAX_EDGE_M = 100.0        # no surface across a gap wider than this between survey points
GEOMETRY_RING_M = (40.0, 600.0)
GEOMETRY_CELL_M = 1.0
GEOMETRY_HALF_BAND_M = 0.5
PROFILE_OVERRIDES = "crop, envelope, bias, step and profile_scale"
TUNINGS_STAY = ("detector tunings other than " + PROFILE_OVERRIDES
                + " stay CACO05's (band_contrast_floor, snap settings)")

# trial (design 6.6)
TRIAL_MAX_FRAMES = 5
TRIAL_MIN_COLUMNS = 0.50      # a line in at least half of the picture's columns
TRIAL_MIN_PLACED = 0.80       # 80% or more of its points placed by georectify
TRIAL_RANGE_M = (30.0, 1000.0)
TRIAL_MARKER = ".waterline_station_check_trial"
STATION_DATA = Path("/mnt/I2Rgus_Data")
# the station's own folders (the cron's, the detector's defaults, the calibration): never a trial's
NEVER_WRITE = ("input", "input_timex", "processed", "processed_timex", "archive", "debug",
               "calibration", "station_calibration")

# pointing (correction 7)
POINTING_PASS_DEG = 0.3
POINTING_FAIL_DEG = 1.0
POINTING_MIN_INLIERS = 40     # pointing_check.measure's own default

GCP_LINE = re.compile(r"^\s*(\d+) +(\S+) +(\S+) +(\S+)\s*$")
SELF_TEST_SUMMARY = re.compile(r"^n=(\d+)\s+mean (\S+) m\s+median (\S+) m\s+max (\S+) m\s+RMS (\S+) m")
FRAME_NAME = re.compile(r"^(\d{9,10})\.")


def _repo_modules():
    """This checkout's folder first on sys.path, for its own modules."""
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))


def _num(x, nd=3):
    """round(x, nd) for a finite number, else None (JSON has no NaN)."""
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return round(x, nd) if math.isfinite(x) else None


def _clean(obj):
    """obj for strict JSON: numpy scalars and arrays as Python ones, NaN and infinities as None."""
    if type(obj).__module__ == "numpy" and hasattr(obj, "tolist"):
        obj = obj.tolist()
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    return obj


def _emit(obj):
    emit(_clean(obj))


def _size_arg(text):
    """'2448x2048' -> (2448, 2048) (argparse type)."""
    m = re.match(r"^\s*(\d+)\s*[xX]\s*(\d+)\s*$", text)
    if not m or int(m.group(1)) == 0 or int(m.group(2)) == 0:
        raise argparse.ArgumentTypeError(f"{text!r}: give the picture size as WIDTHxHEIGHT, e.g. 2448x2048")
    return int(m.group(1)), int(m.group(2))


def _usage(command, message):
    """A usage problem found after argparse (a missing input file ...): JSON, a line, exit 2."""
    _emit({"command": command, "ok": False, "error": message})
    say(f"{command}: {message}")
    return 2


def _missing(paths):
    return [str(p) for p in paths if p and not Path(p).is_file()]


def _image_size(path):
    """(width, height) of a picture from its header, or None when it cannot be read."""
    try:
        from PIL import Image
        with Image.open(str(path)) as im:
            return tuple(int(v) for v in im.size)
    except Exception:
        pass
    try:
        import cv2
        img = cv2.imread(str(path))
        return (int(img.shape[1]), int(img.shape[0])) if img is not None else None
    except Exception:
        return None


# --- pair ------------------------------------------------------------

def read_targets(path):
    """Target coordinates as georectify.self_test reads them (georectify.py:331-339): a header-less
    num,E,N,Z file, positional (C6); rows that are not four numbers are skipped. -> {num: (E, N, Z)}."""
    xyz = {}
    with open(path, "r", newline="") as f:
        for row in csv.reader(f):
            if len(row) < 4:
                continue
            try:
                xyz[int(float(row[0]))] = (float(row[1]), float(row[2]), float(row[3]))
            except ValueError:
                continue
    return xyz


def read_picks(path):
    """Target picks as georectify.self_test reads them (georectify.py:341-344): num,U,V with that header.
    -> {num: (U, V)}. ValueError when the header is not num,U,V or a value is not a number."""
    uv = {}
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        if not {"num", "U", "V"} <= set(reader.fieldnames or ()):
            raise ValueError(f"the header must be num,U,V (it is {','.join(reader.fieldnames or ()) or 'empty'})")
        for row in reader:
            uv[int(float(row["num"]))] = (float(row["U"]), float(row["V"]))
    return uv


def self_test_stats(io_path, eo_path, uv_path, xyz_path):
    """
    georectify.self_test (pinned, called unchanged) for one lens + position + picks + targets, its
    printout captured and read back: -> dict(n, mean_m, median_m, max_m, rms_m, unplaced, pass,
    error, score). pass is self_test's own verdict (its exit 0: mean under 1.0 m). score ranks the
    combinations: the mean, with a target georectify cannot place (beyond its 1000 m, or above the
    horizon) counted as UNPLACED_M off; inf when nothing could be tested.
    """
    _repo_modules()
    import io as _io
    import contextlib
    import warnings
    import numpy as np
    import georectify
    buf = _io.StringIO()
    try:
        with contextlib.redirect_stdout(buf), warnings.catch_warnings(), np.errstate(all="ignore"):
            warnings.simplefilter("ignore")
            rc = georectify.self_test(str(io_path), str(eo_path), str(uv_path), str(xyz_path))
    except Exception as exc:                              # a file georectify cannot read
        return {"n": 0, "mean_m": None, "median_m": None, "max_m": None, "rms_m": None, "unplaced": 0,
                "pass": False, "error": f"{type(exc).__name__}: {exc}", "score": math.inf}
    dist, summary = [], None
    for line in buf.getvalue().splitlines():
        m = GCP_LINE.match(line)
        if m:
            try:
                dist.append(float(m.group(4)))
            except ValueError:
                pass
            continue
        m = SELF_TEST_SUMMARY.match(line.strip())
        if m:
            summary = m
    n = len(dist)
    fin = [d for d in dist if math.isfinite(d)]
    unplaced = n - len(fin)
    stats = [float(summary.group(k)) for k in (2, 3, 4, 5)] if summary else [math.nan] * 4
    if not n:
        score = math.inf
    elif unplaced or not math.isfinite(stats[0]):
        score = (sum(fin) + UNPLACED_M * unplaced) / n
    else:
        score = stats[0]
    return {"n": n, "mean_m": _num(stats[0]), "median_m": _num(stats[1]), "max_m": _num(stats[2]),
            "rms_m": _num(stats[3]), "unplaced": unplaced, "pass": rc == 0 and n > 0,
            "error": None if n else "no targets in common between the picks and the coordinates",
            "score": score}


def _camera_of(eo_path):
    """The EO file's position and angles (degrees, as written) for the report, or {}."""
    _repo_modules()
    from georectify import read_simple_yaml
    try:
        v = read_simple_yaml(eo_path)
    except OSError:
        return {}
    return {k: _num(v.get(k), 4) for k in ("x", "y", "z", "azimuth", "tilt", "roll")}


def _lens_of(io_path):
    _repo_modules()
    from georectify import read_simple_yaml
    try:
        v = read_simple_yaml(io_path)
    except OSError:
        return {}
    return {k: _num(v.get(k), 4) for k in ("NU", "NV", "fx", "fy")}


def pair_sanity(io_path, eo_path, xyz, uv, shared, image=None):
    """The checks after a choice (design 6.4.4, C6): the lens file's NU x NV against the picture size
    (when given), every pick inside the frame, every target 10-1000 m from the camera and at
    -5 m .. camera z + 5 m. -> [{check, ok (True/False, None = not checked), detail}]."""
    lens, cam = _lens_of(io_path), _camera_of(eo_path)
    nu, nv = lens.get("NU"), lens.get("NV")
    out = []
    if image is None:
        out.append({"check": "image_size", "ok": None, "detail": "not checked (no --image)"})
    else:
        ok = nu is not None and nv is not None and (int(nu), int(nv)) == tuple(image)
        out.append({"check": "image_size", "ok": ok,
                    "detail": f"lens file {nu and int(nu)} x {nv and int(nv)}, pictures {image[0]} x {image[1]}"
                              + ("" if ok else ": the lens file is for another picture size")})
    if nu is not None and nv is not None:
        bad = [n for n in shared if not (0 <= uv[n][0] <= nu and 0 <= uv[n][1] <= nv)]
        out.append({"check": "picks_in_frame", "ok": not bad,
                    "detail": f"every pick within 0..{int(nu)} x 0..{int(nv)}" if not bad else
                              f"pick(s) outside the {int(nu)} x {int(nv)} frame: {', '.join(map(str, bad))}"})
    x, y, z = cam.get("x"), cam.get("y"), cam.get("z")
    if x is None or y is None or z is None:
        out.append({"check": "target_range", "ok": False, "detail": "the EO file has no x, y and z"})
        return out
    dist = {n: math.hypot(xyz[n][0] - x, xyz[n][1] - y) for n in shared}
    lo, hi = TARGET_RANGE_M
    bad = [n for n in shared if not lo <= dist[n] <= hi]
    rng = f"{min(dist.values()):.0f}-{max(dist.values()):.0f} m" if dist else "none"
    out.append({"check": "target_range", "ok": bool(shared) and not bad,
                "detail": f"targets {rng} from the camera" + (
                    f"; outside {lo:g}-{hi:g} m: " + ", ".join(f"{n} ({dist[n]:.0f} m)" for n in bad) if bad else "")})
    zhi = z + TARGET_Z_ABOVE_CAMERA_M
    bad = [n for n in shared if not TARGET_Z_LO_M <= xyz[n][2] <= zhi]
    zs = [xyz[n][2] for n in shared]
    out.append({"check": "target_z", "ok": bool(shared) and not bad,
                "detail": (f"target Z {min(zs):.2f} to {max(zs):.2f} m" if zs else "no targets")
                          + f" (allowed {TARGET_Z_LO_M:g} to {zhi:.2f} m NAVD88: camera z + {TARGET_Z_ABOVE_CAMERA_M:g})"
                          + ("; outside: " + ", ".join(f"{n} ({xyz[n][2]:.2f} m)" for n in bad) if bad else "")})
    return out


def cmd_pair(args):
    ios = list(dict.fromkeys(args.io))                    # each file once, in the order given
    eos = list(dict.fromkeys(args.eo))
    missing = _missing([args.xyz, args.uv] + ios + eos)
    if missing:
        return _usage("pair", "not found: " + ", ".join(missing))
    out = {"command": "pair", "xyz": args.xyz, "uv": args.uv, "io": ios, "eo": eos, "ok": False,
           "chosen": "none", "io_tie": None, "eo_tie": None, "ratio_next_eo": None, "next_eo": None,
           "sanity": None, "warnings": [], "reason": None,
           "rule": {"pass_mean_m": 1.0, "ratio_next_eo": PAIR_RATIO, "io_tie_within": IO_TIE,
                    "min_shared": PAIR_MIN_SHARED, "target_range_m": list(TARGET_RANGE_M),
                    "target_z_m": [TARGET_Z_LO_M, "camera z + %g" % TARGET_Z_ABOVE_CAMERA_M]}}
    xyz = read_targets(args.xyz)
    try:
        uv = read_picks(args.uv)
    except (ValueError, KeyError) as exc:
        out["reason"] = f"{Path(args.uv).name}: not a picks file ({exc})"
        _emit(out)
        say(f"pair: REFUSED: {out['reason']}")
        return 1
    shared = sorted(set(xyz) & set(uv))
    out.update(n_xyz=len(xyz), n_uv=len(uv), n_shared=len(shared))

    combos = []
    for eo in eos:
        for io_ in ios:
            combos.append(dict(io=io_, eo=eo, **self_test_stats(io_, eo, args.uv, args.xyz)))
    out["combinations"] = [{k: v for k, v in c.items() if k != "score"} for c in combos]

    # the ratio test is for the choice of POSITION (C5): each EO's best lens, the EOs ranked
    best_of_eo = {}
    for c in combos:
        if math.isfinite(c["score"]) and (c["eo"] not in best_of_eo or c["score"] < best_of_eo[c["eo"]]["score"]):
            best_of_eo[c["eo"]] = c
    ranked = sorted(best_of_eo.values(), key=lambda c: c["score"])
    best = ranked[0] if ranked else None
    nxt = ranked[1] if len(ranked) > 1 else None
    ratio = (nxt["score"] / max(best["score"], 1e-9)) if nxt else None
    if nxt:
        out["ratio_next_eo"] = _num(ratio, 2)
        out["next_eo"] = {"eo": nxt["eo"], "io": nxt["io"], "mean_m": nxt["mean_m"], "unplaced": nxt["unplaced"]}

    def mean_txt(c):
        return f"{c['mean_m']:.2f} m" if c["mean_m"] is not None else f"{c['unplaced']} target(s) not placed"

    if len(shared) < PAIR_MIN_SHARED:
        out["reason"] = (f"only {len(shared)} target(s) shared by the picks and the coordinates "
                         f"({PAIR_MIN_SHARED} or more needed): not paired")
    elif best is None:
        out["reason"] = "no lens and position file could be tested with these picks (see the errors)"
    elif not best["pass"]:
        out["reason"] = (f"the calibration does not pass its test (best mean {mean_txt(best)}): "
                         f"the lens, position and picks do not belong together")
    elif nxt is not None and ratio < PAIR_RATIO:
        out["eo_tie"] = [{"eo": c["eo"], "io": c["io"], "mean_m": c["mean_m"]} for c in ranked if c["pass"]]
        out["reason"] = (f"two positions fit these picks almost equally: {Path(best['eo']).name} "
                         f"({mean_txt(best)}) and {Path(nxt['eo']).name} ({mean_txt(nxt)}), only "
                         f"{ratio:.1f} times worse ({PAIR_RATIO:g} needed): keep one position file")
    if out["reason"] is None:
        # the lens for that position: two passing within IO_TIE of each other is an io_tie, never a refusal (C5)
        same = sorted((c for c in combos if c["eo"] == best["eo"] and c["pass"]), key=lambda c: c["score"])
        tie = [c for c in same if c["score"] <= best["score"] * (1.0 + IO_TIE)]
        out["chosen"] = {"io": best["io"], "eo": best["eo"], "n": best["n"], "mean_m": best["mean_m"],
                         "max_m": best["max_m"], "io_decided": len(tie) == 1,
                         "camera": _camera_of(best["eo"]), "lens": _lens_of(best["io"])}
        if len(tie) > 1:
            entries = []
            for c in tie:
                lens = _lens_of(c["io"])
                fits = None if args.image is None or lens.get("NU") is None else \
                    (int(lens["NU"]), int(lens.get("NV") or 0)) == tuple(args.image)
                entries.append({"io": c["io"], "mean_m": c["mean_m"], "max_m": c["max_m"],
                                "nu_nv": [lens.get("NU"), lens.get("NV")], "image_size_ok": fits})
            out["io_tie"] = {"eo": best["eo"], "ios": entries,
                             "message": f"{len(tie)} lens files pass with this position within "
                                        f"{IO_TIE:.0%} of each other: decide by a name hint (c1, cam1, "
                                        f"CameraA ...) or ask; the chosen io is only the lowest mean"}
        out["sanity"] = pair_sanity(best["io"], best["eo"], xyz, uv, shared, args.image)
        failed = [s for s in out["sanity"] if s["ok"] is False]
        if failed:
            out["reason"] = "; ".join(f"{s['check']}: {s['detail']}" for s in failed)
        out["ok"] = not failed
        if len(shared) < PAIR_GOOD_SHARED:
            out["warnings"].append(f"only {len(shared)} targets shared: {PAIR_GOOD_SHARED} or more is better")

    say(f"pair: {Path(args.uv).name} + {Path(args.xyz).name}: {len(shared)} target(s) shared, "
        f"{len(combos)} lens x position combination(s)")
    for c in combos:
        res = ("ERROR " + c["error"]) if c["error"] and not c["n"] else \
            (f"n={c['n']} mean {mean_txt(c)} max {c['max_m'] if c['max_m'] is not None else '-'} m "
             + ("PASS" if c["pass"] else "FAIL"))
        say(f"  {Path(c['io']).name} + {Path(c['eo']).name}: {res}")
    if out["chosen"] != "none":
        ch = out["chosen"]
        say(f"  chosen: {Path(ch['io']).name} + {Path(ch['eo']).name}, mean {ch['mean_m']:.2f} m"
            + (f"; next best position {Path(nxt['eo']).name} {ratio:.1f} times worse" if nxt else
               "; no other position to compare with"))
        if out["io_tie"]:
            say("  io_tie: " + ", ".join(f"{Path(e['io']).name} ({e['mean_m']:.2f} m)" for e in out["io_tie"]["ios"])
                + ": decide by a name hint or ask")
        for s in out["sanity"]:
            say(f"  {'ok  ' if s['ok'] else 'FAIL' if s['ok'] is False else '--  '} {s['check']}: {s['detail']}")
    for w in out["warnings"]:
        say(f"  warning: {w}")
    say(f"pair: {'OK' if out['ok'] else 'REFUSED: ' + out['reason']}")
    _emit(out)
    return 0 if out["ok"] else 1


# --- profile ---------------------------------------------------------

def _east_north_z_columns(header):
    """Column indexes (E, N, Z) named in a survey file's header, or None."""
    names = [h.strip().strip('"').strip().lower() for h in header]

    def find(test):
        for i, h in enumerate(names):
            if test(h):
                return i
        return None
    e = find(lambda h: h.startswith("east") or h in ("e", "x"))
    n = find(lambda h: h.startswith("north") or h in ("n", "y"))
    z = find(lambda h: h.startswith("elev") or h.startswith("ortho") or h in ("z", "h", "height")
             or h.startswith("z_") or h.startswith("z("))
    return (e, n, z) if None not in (e, n, z) else None


def read_survey_points(path):
    """
    A beach survey's points -> (E, N, Z) numpy arrays (UTM m, NAVD88 m). Comma, semicolon, tab or
    space separated. With a header, the columns it names easting/northing/elevation (E/N/Z, x/y/z;
    other columns, such as a text ID, are not read); without one, positional as the normalised
    files: three numbers E,N,Z, four or more num,E,N,Z (C6). Rows that are not numbers are skipped.
    ValueError when no UTM points are left.
    """
    import numpy as np
    rows, cols, seen = [], None, False
    with open(path, "r", newline="", encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = [p.strip() for p in (s.split(",") if "," in s else s.split(";") if ";" in s else s.split())]
            first, seen = not seen, True
            if cols is not None:
                try:
                    rows.append([float(parts[c]) for c in cols])
                except (ValueError, IndexError):
                    pass
                continue
            try:
                vals = [float(p) for p in parts]
            except ValueError:
                if first:                                 # a header: the columns it names
                    cols = _east_north_z_columns(parts)
                    if cols is None:
                        raise ValueError(f"its header ({s[:80]}) names no easting, northing and elevation columns")
                continue
            if len(vals) == 3:
                rows.append(vals)
            elif len(vals) >= 4:
                rows.append(vals[1:4])
    if not rows:
        raise ValueError("no rows of numbers")
    a = np.array(rows, float)
    ok = (a[:, 0] >= 166000) & (a[:, 0] <= 834000) & (a[:, 1] >= 1e6) & (a[:, 1] <= 1e7) & np.isfinite(a[:, 2])
    if not ok.any():
        raise ValueError("no UTM easting/northing (E 166000-834000, N 1,000,000-10,000,000): "
                         "other columns or another coordinate system")
    return a[ok, 0], a[ok, 1], a[ok, 2]


def survey_grid(E, N, Z, cx, cy, max_range=SURVEY_RANGE_M, cell=SURVEY_CELL_M, max_edge=TRI_MAX_EDGE_M):
    """
    Survey points -> the north-up grid envelope_from_survey wants, (grid rows north to south, x_left,
    y_top, cell), NaN where there is no survey: the points triangulated (matplotlib.tri) and
    interpolated linearly at the centres of `cell` cells over the square +-max_range box round the
    camera, no surface across gaps wider than max_edge. ValueError when nothing is left.
    """
    import numpy as np
    import matplotlib.tri as mtri
    margin = max_range + max_edge
    keep = (np.abs(E - cx) <= margin) & (np.abs(N - cy) <= margin)
    E, N, Z = E[keep], N[keep], Z[keep]
    if len(E):
        _, idx = np.unique(np.round(np.c_[E, N], 3), axis=0, return_index=True)   # one point per mm
        E, N, Z = E[idx], N[idx], Z[idx]
    if len(E) < 3:
        raise ValueError(f"fewer than 3 survey points within {max_range:g} m of the camera")
    try:
        tri = mtri.Triangulation(E, N)
    except (RuntimeError, ValueError) as exc:
        raise ValueError(f"the survey points cannot be triangulated ({exc}): all on one line?")
    t = tri.triangles
    edge = np.max([np.hypot(E[t[:, i]] - E[t[:, j]], N[t[:, i]] - N[t[:, j]])
                   for i, j in ((0, 1), (1, 2), (2, 0))], axis=0)
    if not (edge <= max_edge).any():
        raise ValueError(f"the survey points are too sparse (no triangle with sides under {max_edge:g} m)")
    tri.set_mask(edge > max_edge)
    x0, x1 = max(E.min(), cx - max_range), min(E.max(), cx + max_range)
    y0, y1 = max(N.min(), cy - max_range), min(N.max(), cy + max_range)
    if x0 >= x1 or y0 >= y1:
        raise ValueError(f"the survey is not within {max_range:g} m of the camera")
    x_left = math.floor(x0 / cell) * cell
    y_top = math.ceil(y1 / cell) * cell
    ncol = int(math.ceil((x1 - x_left) / cell))
    nrow = int(math.ceil((y_top - y0) / cell))
    XX, YY = np.meshgrid(x_left + (np.arange(ncol) + 0.5) * cell, y_top - (np.arange(nrow) + 0.5) * cell)
    zi = mtri.LinearTriInterpolator(tri, Z)(XX, YY)
    grid = np.ma.filled(np.ma.masked_invalid(np.ma.asarray(zi).astype(float)), np.nan)
    if not np.isfinite(grid).any():
        raise ValueError(f"the survey covers no ground within {max_range:g} m of the camera")
    return grid, x_left, y_top, cell


def geometry_grid(cx, cy, z, ring=GEOMETRY_RING_M, cell=GEOMETRY_CELL_M):
    """A level plane at z over the ring ring[0]..ring[1] m round the camera (NaN elsewhere), as the
    north-up grid envelope_from_survey wants: the search band from geometry alone (C8)."""
    import numpy as np
    n = int(math.ceil(2 * ring[1] / cell))
    x_left, y_top = cx - n * cell / 2.0, cy + n * cell / 2.0
    XX, YY = np.meshgrid(x_left + (np.arange(n) + 0.5) * cell, y_top - (np.arange(n) + 0.5) * cell)
    r = np.hypot(XX - cx, YY - cy)
    return np.where((r >= ring[0]) & (r <= ring[1]), float(z), np.nan), x_left, y_top, cell


def _grid_covers(grid, x_left, y_top, cell, cx, cy, max_range):
    """True when the grid overlaps the square +-max_range box round the camera."""
    h, w = grid.shape
    return (x_left < cx + max_range and x_left + w * cell > cx - max_range
            and y_top > cy - max_range and y_top - h * cell < cy + max_range)


def cmd_profile(args):
    out = {"command": "profile", "camera": args.camera, "ok": False, "status": "REFUSED", "output": args.output,
           "io": args.io, "eo": args.eo, "provisional_reason": None, "note": TUNINGS_STAY}
    missing = _missing([args.io, args.eo, args.survey_csv, args.dem])
    if missing:
        return _usage("profile", "not found: " + ", ".join(missing))
    if (args.mllw is None) != (args.mhhw is None):
        return _usage("profile", "--mllw and --mhhw go together")
    if args.mllw is not None and not args.mllw < args.mhhw:
        return _usage("profile", f"--mllw {args.mllw:g} must be below --mhhw {args.mhhw:g}")
    if args.plane_z is not None and not args.geometry:
        return _usage("profile", "--plane-z goes with --geometry (the plane's height, the gauge's MSL)")

    def refuse(reason):
        out["reason"] = reason
        _emit(out)
        say(f"profile {args.camera}: REFUSED: {reason}")
        return 1

    _repo_modules()
    import numpy as np
    from georectify import load_intrinsics, load_extrinsics
    try:
        io_, eo = load_intrinsics(args.io), load_extrinsics(args.eo)
    except (OSError, ValueError) as exc:
        return refuse(f"calibration: {exc}")
    nu, nv = int(io_[0]), int(io_[1])
    if (nu, nv) != tuple(args.image_size):
        return refuse(f"the lens file is {nu} x {nv} but the pictures are {args.image_size[0]} x "
                      f"{args.image_size[1]}: another lens file, or reduced copies")
    cx, cy = float(eo[0]), float(eo[1])
    notes = []
    if args.geometry:
        z = 0.0 if args.plane_z is None else args.plane_z
        band = [z - GEOMETRY_HALF_BAND_M, z + GEOMETRY_HALF_BAND_M]
        max_range = GEOMETRY_RING_M[1]
        survey = geometry_grid(cx, cy, z)
        out.update(source="geometry", source_file=None, plane_z_m=z, ring_m=list(GEOMETRY_RING_M),
                   band_from=f"a level plane at {z:+.2f} m NAVD88 (the gauge's MSL)")
        if args.mllw is not None:
            notes.append("--mllw/--mhhw not used: the band from geometry is the plane alone")
        reason = (f"search band from the camera geometry (a level plane at {z:+.2f} m NAVD88, "
                  f"{GEOMETRY_RING_M[0]:g}-{GEOMETRY_RING_M[1]:g} m from the camera), not from a beach survey or DEM")
        out["provisional_reason"] = reason
        out["station_provisional"] = (f"search band from camera geometry (level plane at {z:+.2f} m, "
                                      f"{GEOMETRY_RING_M[0]:g}-{GEOMETRY_RING_M[1]:g} m); {TUNINGS_STAY}")
    else:
        if args.mllw is not None:
            band = [args.mllw - DATUM_PAD_M[0], args.mhhw + DATUM_PAD_M[1]]
            out["band_from"] = (f"the gauge's MLLW {args.mllw:+.2f} - {DATUM_PAD_M[0]:g} m to MHHW "
                                f"{args.mhhw:+.2f} + {DATUM_PAD_M[1]:g} m (NAVD88)")
        else:
            band = list(SURVEY_BAND_M)
            out["band_from"] = f"the default {SURVEY_BAND_M[0]:+g} to {SURVEY_BAND_M[1]:+g} m NAVD88 (no gauge datums given)"
        max_range = SURVEY_RANGE_M
        if args.survey_csv:
            out.update(source="survey", source_file=args.survey_csv)
            try:
                E, N, Z = read_survey_points(args.survey_csv)
                d = np.hypot(E - cx, N - cy)
                near = d <= SURVEY_MAX_DIST_M
                out["survey_points"] = int(near.sum())
                if (~near).any():
                    notes.append(f"{int((~near).sum())} survey point(s) more than {SURVEY_MAX_DIST_M:g} m "
                                 f"from the camera ignored")
                survey = survey_grid(E[near], N[near], Z[near], cx, cy)
            except (OSError, ValueError) as exc:
                return refuse(f"{Path(args.survey_csv).name}: {exc}")
        else:
            out.update(source="dem", source_file=args.dem)
            from compare_dem_survey import read_survey
            try:
                survey = read_survey(args.dem)
            except SystemExit as exc:                       # read_survey says what it cannot read
                return refuse(f"{Path(args.dem).name}: {exc.code}")
            except (OSError, ValueError, KeyError) as exc:
                return refuse(f"{Path(args.dem).name}: cannot read ({exc})")
            if not _grid_covers(*survey, cx, cy, max_range):
                return refuse(f"{Path(args.dem).name} does not cover the ground within {max_range:g} m of the "
                              f"camera (another UTM zone or coordinate system?)")
        out["station_provisional"] = TUNINGS_STAY
    out.update(band_m=[_num(band[0]), _num(band[1])], range_m=max_range)

    from detect_original_view import envelope_from_survey
    try:
        with np.errstate(invalid="ignore"):
            env, left_far, frac_ok = envelope_from_survey(io_, eo, survey, band[0], band[1], max_range)
    except SystemExit as exc:                               # "the survey band hardly appears ...": C8
        return refuse(f"{exc.code} ({out['source']}, band {band[0]:+.2f} to {band[1]:+.2f} m, "
                      f"within {max_range:g} m)")
    # exactly as detect_original_view.py:343-347 writes it
    prof = {"crop": [0.0, 1.0, 0.0, 1.0], "envelope": env, "no_bias": True}
    if left_far:      # profiles assume column 0 = near the camera; swap the per-column ramps
        prof["step"] = [4, 18]
        prof["profile_scale"] = [0.4, 1.0]
    data = {}
    target = Path(args.output)
    if args.merge and target.is_file():
        try:
            data = json.loads(target.read_text())
        except ValueError as exc:
            return refuse(f"--merge: {target} is not JSON ({exc})")
        if not isinstance(data, dict):
            return refuse(f"--merge: {target} is not a profile file")
        data.pop(args.camera, None)
    data[args.camera] = prof
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        tmp = target.with_name(target.name + ".tmp")
        tmp.write_text(json.dumps(data, indent=1))
        tmp.replace(target)
    except OSError as exc:
        return refuse(f"cannot write {target}: {exc}")
    out.update(ok=True, status="PROVISIONAL" if args.geometry else "OK", profile=prof,
               cameras_in_file=sorted(data), envelope_points=len(env),
               far_side="left" if left_far else "right", in_view_fraction=_num(frac_ok, 3), notes=notes)
    _emit(out)
    say(f"profile {args.camera}: {out['status']}: {len(env)} envelope points from {out['source']} "
        f"(band {band[0]:+.2f} to {band[1]:+.2f} m, within {max_range:g} m; {100 * frac_ok:.0f}% of its cells "
        f"in view); far side {out['far_side']} -> {target}")
    if out["provisional_reason"]:
        say(f"  PROVISIONAL: {out['provisional_reason']}")
    say(f"  note: {TUNINGS_STAY}")
    for n in notes:
        say(f"  note: {n}")
    return 0


# --- trial -----------------------------------------------------------

def _inside(path, folder):
    try:
        path.relative_to(folder)
        return True
    except ValueError:
        return False


def trial_dir_problem(work):
    """Why `work` may not hold the trial's scratch folders, or None: never inside the station's
    own folders (input/, processed_timex/, archive/ ... of this checkout or of the station's
    waterline folder), nor the camera's ImageProducts."""
    w = Path(os.path.abspath(str(work)))
    try:
        w = w.resolve()
    except OSError:
        pass
    roots = [HERE / name for name in NEVER_WRITE] + [STATION_DATA / "waterline" / name for name in NEVER_WRITE]
    roots.append(STATION_DATA / "ImageProducts")
    for r in roots:
        for f in {r, Path(os.path.realpath(str(r)))}:
            if _inside(w, f):
                return f"{w} is inside {f}, one of the station's own folders: give a scratch folder"
    return None


def _run_logged(cmd, log, cwd, env=None):
    """Runs cmd with its output appended to log; -> exit status."""
    with open(log, "a") as f:
        f.write("\n$ " + " ".join(str(c) for c in cmd) + "\n")
        f.flush()
        p = subprocess.run([str(c) for c in cmd], stdout=f, stderr=subprocess.STDOUT, cwd=str(cwd), env=env)
    return p.returncode


def _read_csv(path):
    try:
        with open(path, newline="") as f:
            return list(csv.DictReader(f))
    except OSError:
        return []


def _frame_result(reasons=None):
    """A frame's entry in the trial's JSON, before it is judged (every frame has the same keys)."""
    return {"detected": False, "signal_columns_fraction": None, "points": 0, "placed": 0,
            "placed_fraction": None, "range_m": None, "water_level_navd88": None, "pass": False,
            "reasons": list(reasons or [])}


def judge_frame(stem, det_csv, ground_rows, nu, cam_xy):
    """The trial's verdict on one frame (design 6.6): a line in >= TRIAL_MIN_COLUMNS of the picture's
    columns, >= TRIAL_MIN_PLACED of its points placed, all of them TRIAL_RANGE_M from the camera."""
    res = _frame_result()
    det = _read_csv(det_csv) if det_csv.is_file() else []
    if not det:
        res["reasons"].append("the detector kept no line (no usable signal, too dark or low confidence: "
                              "see the log and the debug images)")
        return res
    res["detected"] = True
    sig = sum(1 for r in det if (r.get("Has_Signal") or "1").strip() in ("1", "1.0", "True"))
    frac = sig / float(nu)
    res["signal_columns_fraction"] = _num(frac, 3)
    if frac < TRIAL_MIN_COLUMNS:
        res["reasons"].append(f"a line in {frac:.0%} of the columns ({TRIAL_MIN_COLUMNS:.0%} needed)")
    rows = [r for r in ground_rows if r.get("source_file") == stem]
    res["points"] = len(rows)
    if not rows:
        res["reasons"].append("no contour points (extract_elevation_contours.py left the frame out: the "
                              "water level does not cover it, or its line filters; see the log)")
        return res
    lev = rows[0].get("beach_elevation_navd88") or rows[0].get("tide_elevation_navd88")
    res["water_level_navd88"] = _num(lev)
    dist = []
    for r in rows:
        try:
            e, n = float(r.get("easting_utm19") or "nan"), float(r.get("northing_utm19") or "nan")
        except ValueError:
            continue
        if math.isfinite(e) and math.isfinite(n):
            dist.append(math.hypot(e - cam_xy[0], n - cam_xy[1]))
    res["placed"] = len(dist)
    pf = len(dist) / float(len(rows))
    res["placed_fraction"] = _num(pf, 3)
    if pf < TRIAL_MIN_PLACED:
        res["reasons"].append(f"{pf:.0%} of its points placed ({TRIAL_MIN_PLACED:.0%} needed)")
    if dist:
        res["range_m"] = [_num(min(dist), 1), _num(max(dist), 1)]
        lo, hi = TRIAL_RANGE_M
        out_of = sum(1 for d in dist if not lo <= d <= hi)
        if out_of:
            res["reasons"].append(f"{out_of} point(s) outside {lo:g}-{hi:g} m from the camera "
                                  f"(they lie {min(dist):.0f}-{max(dist):.0f} m)")
    res["pass"] = not res["reasons"]
    return res


def cmd_trial(args):
    cam = args.camera
    frames = list(dict.fromkeys(args.frames))
    if len(frames) > TRIAL_MAX_FRAMES:
        return _usage("trial", f"{len(frames)} frames given: {TRIAL_MAX_FRAMES} at most")
    missing = _missing([args.io, args.eo, args.profile_json, args.spline, args.waves, args.gauge_csv] + frames)
    if missing:
        return _usage("trial", "not found: " + ", ".join(missing))
    bad = [f for f in frames if not (FRAME_NAME.match(Path(f).name) and f".{cam}." in Path(f).name.lower()
                                     and Path(f).name.endswith(".timex.jpg"))]
    if bad:
        return _usage("trial", f"not a {cam} timex named <epoch>....{cam}.timex.jpg: " + ", ".join(bad))
    if args.work_dir:
        problem = trial_dir_problem(args.work_dir)
        if problem:
            return _usage("trial", problem)
    # the antenna height: extract_elevation_contours.py would otherwise use CACO-05's 19.014 m and every
    # water level would be off by the difference while the trial still passes (its range check cannot see it)
    if args.antenna_navd88 is not None:
        antenna, antenna_from, raw = args.antenna_navd88, "--antenna-navd88", repr(args.antenna_navd88)
    else:
        raw = os.environ.get("GNSSR_ANTENNA_NAVD88_M", "").strip()
        if not raw:
            return _usage("trial", "--antenna-navd88 M is needed (the station's antenna height in m NAVD88, the "
                                   "OPUS height): without it extract_elevation_contours.py would use CACO-05's "
                                   "19.014 m and put every water level off by the difference")
        antenna, antenna_from = _float(raw), "GNSSR_ANTENNA_NAVD88_M in the environment"
    if antenna is None or not math.isfinite(antenna) or not -100 < antenna < 5000:
        return _usage("trial", f"{antenna_from} {raw}: not an antenna height in m NAVD88")
    scripts = Path(args.scripts_dir).resolve()
    out = {"command": "trial", "camera": cam, "ok": False, "frames": [], "passed": 0, "tried": len(frames),
           "best": None, "steps": {}, "work_dir": None, "log": None,
           "antenna_navd88_m": antenna, "antenna_from": antenna_from,
           "rule": {"min_columns_fraction": TRIAL_MIN_COLUMNS, "min_placed_fraction": TRIAL_MIN_PLACED,
                    "range_m": list(TRIAL_RANGE_M)}}

    def refuse(reason):
        out["reason"] = reason
        _emit(out)
        say(f"trial {cam}: REFUSED: {reason}")
        return 1

    try:
        prof = json.loads(Path(args.profile_json).read_text())
    except (OSError, ValueError) as exc:
        return refuse(f"{args.profile_json}: not a profile file ({exc})")
    if not isinstance(prof, dict) or cam not in prof:
        return refuse(f"{args.profile_json} has no entry for {cam} (the detector would use CACO05's own "
                      f"crop and search band): make it with the profile command")
    _repo_modules()
    from georectify import load_intrinsics, load_extrinsics
    try:
        io_, eo = load_intrinsics(args.io), load_extrinsics(args.eo)
    except (OSError, ValueError) as exc:
        return refuse(f"calibration: {exc}")
    nu, nv = int(io_[0]), int(io_[1])

    import shutil
    import tempfile
    # absolute: every step runs with cwd=t and is given paths inside t (a relative one would mean t/t/...)
    work = Path(os.path.abspath(args.work_dir)) if args.work_dir else Path(tempfile.mkdtemp(prefix="waterline_trial_"))
    t = work / f"trial_{cam}"
    if t.exists() and not (t / TRIAL_MARKER).is_file():
        return refuse(f"{t} exists and was not made by this command: not touching it")
    shutil.rmtree(str(t), ignore_errors=True)
    for sub in ("src", "in", "detections", "debug"):
        (t / sub).mkdir(parents=True)
    (t / TRIAL_MARKER).write_text("scratch folder of waterline_station_check.py trial; replaced on each run\n")
    log = t / "trial.log"
    log.write_text("")
    out.update(work_dir=str(t), log=str(log))

    results, usable = {}, []
    for f in frames:
        name = Path(f).name
        stem = name[:-len(".jpg")]
        size = _image_size(f)
        if size is None:
            results[name] = dict(frame=name, **_frame_result(["cannot read the picture"]))
            continue
        if size != (nu, nv):
            results[name] = dict(frame=name, **_frame_result(
                [f"the picture is {size[0]} x {size[1]}, the lens file {nu} x {nv}"]))
            continue
        shutil.copy2(f, str(t / "src" / name))            # a copy: nothing in the station's folders is touched
        usable.append((name, stem))

    if usable:
        rc = _run_logged([sys.executable, scripts / "waterline_detector_v5.py",
                          "--image-suffix", "timex.jpg", "--timing-mode", args.timing_mode,
                          "--source-dir", t / "src", "--input-dir", t / "in",
                          "--output-dir", t / "detections", "--debug-dir", t / "debug",
                          "--no-bias-correction", "--profile-json", Path(args.profile_json).resolve()], log, t)
        out["steps"]["detector"] = rc
        contours, ground = t / "contour_points.csv", t / "contour_points_ground.csv"
        if any((t / "detections").glob("*.csv")):
            env = dict(os.environ)
            env["GNSSR_ANTENNA_NAVD88_M"] = repr(float(antenna))
            extra = (["--waves", Path(args.waves).resolve()] if args.waves else []) + \
                    (["--gnssr-qc-reference", Path(args.gauge_csv).resolve()] if args.gauge_csv else [])
            rc = _run_logged([sys.executable, scripts / "extract_elevation_contours.py", Path(args.spline).resolve(),
                              "--processed-dir", t / "detections", "--output", contours,
                              "--timing-mode", args.timing_mode] + extra, log, t, env=env)
            out["steps"]["contours"] = rc
            if rc == 0 and contours.is_file():
                rc = _run_logged([sys.executable, scripts / "georectify.py", contours, ground,
                                  f"--io-{cam}", Path(args.io).resolve(), f"--eo-{cam}", Path(args.eo).resolve()],
                                 log, t)
                out["steps"]["georectify"] = rc
        ground_rows = _read_csv(ground) if out["steps"].get("georectify") == 0 else []
        for name, stem in usable:
            res = judge_frame(stem, t / "detections" / (stem + ".csv"), ground_rows, nu, (float(eo[0]), float(eo[1])))
            if res["detected"] and not res["points"] and out["steps"].get("georectify") not in (None, 0):
                res["reasons"].append("georectify.py failed (see the log)")
            results[name] = dict(frame=name, **res)

    for f in frames:
        name = Path(f).name
        r = results[name]
        m = FRAME_NAME.match(name)
        r["time_utc"] = datetime.fromtimestamp(int(m.group(1)), tz=timezone.utc).strftime("%Y-%m-%d %H:%M")
        out["frames"].append(r)
    passed = [r for r in out["frames"] if r["pass"]]
    out["passed"] = len(passed)
    out["ok"] = bool(passed)
    if passed:
        b = max(passed, key=lambda r: (r["points"], r["placed_fraction"] or 0))
        out["best"] = {"frame": b["frame"], "time_utc": b["time_utc"], "points": b["points"],
                       "placed_percent": _num(100 * (b["placed_fraction"] or 0), 0), "range_m": b["range_m"]}
    else:
        out["reason"] = f"no frame passed ({len(frames)} tried): the maps may stay empty until the detector is tuned"
    _emit(out)
    say(f"trial {cam}: {len(passed)} of {len(frames)} frame(s) pass (work folder {t}; antenna "
        f"{antenna:.3f} m NAVD88 from {antenna_from})")
    for r in out["frames"]:
        if r["pass"]:
            say(f"  ok    {r['time_utc']}: {r['points']} points, {100 * r['placed_fraction']:.0f}% placed, "
                f"{r['range_m'][0]:.0f}-{r['range_m'][1]:.0f} m from the camera")
        else:
            say(f"  FAIL  {r['time_utc']}: " + "; ".join(r["reasons"]))
    return 0 if passed else 1


# --- pointing --------------------------------------------------------

def pointing_grade(rotation_deg):
    """PASS under 0.3 deg, WARN 0.3-1.0, FAIL over 1.0 (correction 7)."""
    if rotation_deg < POINTING_PASS_DEG:
        return "PASS"
    return "WARN" if rotation_deg <= POINTING_FAIL_DEG else "FAIL"


def cmd_pointing(args):
    out = {"command": "pointing", "ok": False, "verdict": None, "calibration_picture": args.calibration_picture,
           "io": args.io, "eo": args.eo, "frames": [], "rotation_deg": None,
           "rule": {"pass_below_deg": POINTING_PASS_DEG, "fail_above_deg": POINTING_FAIL_DEG,
                    "min_inliers": args.min_inliers}}
    timex = list(dict.fromkeys(args.timex))
    missing = _missing([args.io, args.eo, args.calibration_picture] + timex)
    if missing:
        return _usage("pointing", "not found: " + ", ".join(missing))

    def stop(verdict, reason):
        out.update(verdict=verdict, reason=reason)
        _emit(out)
        say(f"pointing: {verdict}: {reason}")
        return 1

    _repo_modules()
    try:
        import numpy as np
        import cv2
        import pointing_check as pc
        from georectify import load_intrinsics, load_extrinsics, cirn_angles_to_R
    except ImportError as exc:
        return stop("UNAVAILABLE", f"the picture matcher needs OpenCV (python3-opencv), which does not import "
                                   f"here ({exc}): nothing checked")
    try:
        io_, eo = load_intrinsics(args.io), load_extrinsics(args.eo)
    except (OSError, ValueError) as exc:
        return stop("ERROR", f"calibration: {exc}")
    nu, nv = int(io_[0]), int(io_[1])
    ref_img = cv2.imread(str(args.calibration_picture))
    if ref_img is None:
        return stop("ERROR", f"cannot read the calibration-day picture {args.calibration_picture}")
    if (ref_img.shape[1], ref_img.shape[0]) != (nu, nv):
        return stop("ERROR", f"the calibration-day picture is {ref_img.shape[1]} x {ref_img.shape[0]}, "
                             f"the lens file {nu} x {nv}")
    import tempfile
    with tempfile.TemporaryDirectory(prefix="waterline_pointing_") as td:
        # pointing_check.Reference reads <picture> and <picture>.json (its pointing): a link to the
        # picture in a scratch folder, so nothing is written next to it
        src = Path(args.calibration_picture).resolve()
        ref_path = Path(td) / ("reference" + (src.suffix or ".jpg"))
        try:
            os.symlink(str(src), str(ref_path))
        except OSError:
            import shutil
            shutil.copy2(str(src), str(ref_path))
        ref_path.with_suffix(".json").write_text(json.dumps({"eo": eo.tolist(), "eo_file": Path(args.eo).name}))
        try:
            ref = pc.Reference(ref_path)
        except Exception as exc:                          # an OpenCV error on an odd picture
            return stop("ERROR", f"cannot use the calibration-day picture ({type(exc).__name__}: {exc})")
        for p in timex:
            fr = {"frame": Path(p).name}
            out["frames"].append(fr)
            img = cv2.imread(str(p))
            if img is None:
                fr.update(status="unreadable")
                continue
            if (img.shape[1], img.shape[0]) != (nu, nv):
                fr.update(status="wrong_size", size=[img.shape[1], img.shape[0]])
                continue
            try:
                r = pc.measure(img, [ref], io_, min_inliers=args.min_inliers)
            except Exception as exc:                      # too few points for RANSAC and the like
                fr.update(status="unmatched", matches=0, inliers=0, error=f"{type(exc).__name__}: {exc}")
                continue
            if r.get("status") == "unmatched":
                fr.update(status="unmatched", matches=int(r.get("matches", 0)), inliers=int(r.get("inliers", 0)))
                continue
            rot = pc.rotation_angle_deg(ref.R, cirn_angles_to_R(*r["eo"][3:6]))
            fr.update(status="matched", rotation_deg=_num(rot, 3), d_azimuth_deg=_num(r["d_azimuth"], 3),
                      d_tilt_deg=_num(r["d_tilt"], 3), d_roll_deg=_num(r["d_roll"], 3),
                      matches=int(r["matches"]), inliers=int(r["inliers"]),
                      residual_deg=_num(r["residual_deg"], 4), turn_sigma=_num(r.get("turn_sigma"), 2),
                      verdict=pointing_grade(rot))
    matched = [f for f in out["frames"] if f["status"] == "matched"]
    compared = [f for f in out["frames"] if f["status"] in ("matched", "unmatched")]
    if not compared:
        return stop("ERROR", "no timex could be compared: " + ", ".join(
            f"{f['frame']} ({f['status']}" + (f" {f['size'][0]} x {f['size'][1]}, the lens file {nu} x {nv}"
                                              if "size" in f else "") + ")" for f in out["frames"]))
    if not matched:
        best = max(compared, key=lambda f: f["inliers"])
        return stop("UNMATCHED", f"not enough matches between the calibration-day picture and the timex (best "
                                 f"{best['inliers']} inliers of {best['matches']} matches, {args.min_inliers} "
                                 f"needed): fog, night, rain on the lens, another camera's picture, or a camera "
                                 f"that turned far; nothing checked")
    rot = float(np.median([f["rotation_deg"] for f in matched]))
    verdict = pointing_grade(rot)
    out.update(verdict=verdict, ok=verdict in ("PASS", "WARN"), rotation_deg=_num(rot, 3),
               matched=len(matched),
               shift_at_300m_m=_num(300.0 * math.tan(math.radians(rot)), 2))
    out["reason"] = {
        "PASS": f"the camera has not turned since the calibration day ({rot:.2f} deg)",
        "WARN": f"the camera turned {rot:.2f} deg since the calibration day (about "
                f"{out['shift_at_300m_m']:.1f} m at 300 m): the maps will be off by that much",
        "FAIL": f"the camera turned {rot:.2f} deg since the calibration day (about "
                f"{out['shift_at_300m_m']:.1f} m at 300 m): the calibration no longer fits; calibrate again"}[verdict]
    _emit(out)
    say(f"pointing: {verdict}: {out['reason']} (median of {len(matched)} matched frame(s))")
    for f in out["frames"]:
        if f["status"] == "matched":
            say(f"  {f['frame']}: {f['rotation_deg']:.3f} deg (azimuth {f['d_azimuth_deg']:+.3f}, tilt "
                f"{f['d_tilt_deg']:+.3f}, roll {f['d_roll_deg']:+.3f}), {f['inliers']}/{f['matches']} inliers")
        else:
            say(f"  {f['frame']}: {f['status']}")
    return 0 if out["ok"] else 1


def add_station_setup_parsers(sub):
    """The parsers of pair, profile, trial and pointing."""
    p = sub.add_parser("pair", help="which lens (IO), position (EO), target and picks files belong together "
                                    "(design 6.4): georectify's self-test on every IO x EO")
    p.add_argument("--xyz", required=True, help="target coordinates: header-less num,E,N,Z (UTM m, NAVD88 m)")
    p.add_argument("--uv", required=True, help="the picks of those targets in one camera: num,U,V with that header")
    p.add_argument("--io", nargs="+", required=True, metavar="IO", help="candidate lens files (key: value)")
    p.add_argument("--eo", nargs="+", required=True, metavar="EO", help="candidate position files (key: value)")
    p.add_argument("--image", type=_size_arg, default=None, metavar="WxH",
                   help="the camera's picture size (its newest timex): checked against the chosen lens's NU x NV")
    p.set_defaults(func=cmd_pair)

    p = sub.add_parser("profile", help="the detector's search band from a beach survey, a DEM or the camera "
                                       "geometry (design 6.5): writes the --profile-json file")
    p.add_argument("--camera", required=True, choices=("c1", "c2"))
    p.add_argument("--io", required=True)
    p.add_argument("--eo", required=True)
    p.add_argument("--image-size", required=True, type=_size_arg, metavar="WxH")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--survey-csv", default=None, metavar="CSV",
                     help="beach survey points (E,N,Z; num,E,N,Z; or a header naming them), NAVD88")
    src.add_argument("--dem", default=None, metavar="TIF_OR_ASC",
                     help="a beach DEM, GeoTIFF or Arc ASCII grid, in the camera's UTM zone, NAVD88")
    src.add_argument("--geometry", action="store_true",
                     help="no survey: a level plane at --plane-z, 40-600 m from the camera (PROVISIONAL)")
    p.add_argument("--plane-z", type=float, default=None, metavar="M",
                   help="with --geometry: the plane's height, the gauge's MSL in m NAVD88 (default 0)")
    p.add_argument("--mllw", type=float, default=None, metavar="M", help="the gauge's MLLW, m NAVD88")
    p.add_argument("--mhhw", type=float, default=None, metavar="M", help="the gauge's MHHW, m NAVD88")
    p.add_argument("--output", required=True, help="the profile JSON to write (detector --profile-json)")
    p.add_argument("--merge", action="store_true",
                   help="keep the other camera's entry of an existing --output file")
    p.set_defaults(func=cmd_profile)

    p = sub.add_parser("trial", help="the detector, the contours and georectify on up to 5 kept timex of one "
                                     "camera, in scratch folders only (design 6.6)")
    p.add_argument("--camera", required=True, choices=("c1", "c2"))
    p.add_argument("--io", required=True)
    p.add_argument("--eo", required=True)
    p.add_argument("--profile-json", required=True, help="the profile command's file")
    p.add_argument("--spline", required=True, help="the GNSS-R water levels (gnssrefl spline_out.txt)")
    p.add_argument("--frames", nargs="+", required=True, metavar="TIMEX",
                   help=f"up to {TRIAL_MAX_FRAMES} timex pictures (<epoch>....<camera>.timex.jpg); only read")
    p.add_argument("--work-dir", default=None,
                   help="scratch folder (default: a new temporary one); trial_<camera>/ in it is replaced")
    p.add_argument("--antenna-navd88", type=float, default=None, metavar="M",
                   help="GNSSR_ANTENNA_NAVD88_M for extract_elevation_contours.py (the OPUS height); needed "
                        "unless GNSSR_ANTENNA_NAVD88_M is set in the environment (never CACO-05's 19.014)")
    p.add_argument("--waves", default=None, help="wave archive for extract_elevation_contours.py --waves")
    p.add_argument("--gauge-csv", default=None, help="tide-gauge archive for --gnssr-qc-reference")
    p.add_argument("--timing-mode", default="legacy", choices=("legacy", "mid"))
    p.add_argument("--scripts-dir", default=str(HERE), help=argparse.SUPPRESS)    # tests: stub scripts
    p.set_defaults(func=cmd_trial)

    p = sub.add_parser("pointing", help="a calibration-day picture against recent timex of the same camera: "
                                        "has the camera turned? (PASS < 0.3 deg, WARN 0.3-1.0, FAIL > 1.0)")
    p.add_argument("--io", required=True)
    p.add_argument("--eo", required=True, help="the calibration's position file (the picture's pointing)")
    p.add_argument("--calibration-picture", required=True, metavar="IMAGE",
                   help="a picture this camera took on the calibration day")
    p.add_argument("--timex", nargs="+", required=True, metavar="TIMEX", help="recent timex of the same camera")
    p.add_argument("--min-inliers", type=int, default=POINTING_MIN_INLIERS, help=argparse.SUPPRESS)
    p.set_defaults(func=cmd_pointing)


def build_parser():
    ap = argparse.ArgumentParser(
        prog="waterline_station_check.py",
        description="Checks for a new station's waterline settings (coastcam-station-setup). "
                    "Exit 0 ok, 1 refused or failed, 2 usage.")
    sub = ap.add_subparsers(dest="command", metavar="COMMAND")
    sub.required = True

    p = sub.add_parser("capabilities", help=f"print this checkout's feature tags ({CAPABILITIES})")
    p.set_defaults(func=cmd_capabilities)

    p = sub.add_parser("check-overlay", help="check a waterline_station.conf before the cron sources it")
    p.add_argument("--conf", required=True, help="the file (waterline_station.conf)")
    p.set_defaults(func=cmd_check_overlay)

    add_station_setup_parsers(sub)
    return ap


def main(argv=None):
    ap = build_parser()
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

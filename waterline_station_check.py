#!/usr/bin/env python3
"""
Checks For A New Station's Waterline Settings (coastcam-station-setup)
=======================================================================
    python3 waterline_station_check.py capabilities
    python3 waterline_station_check.py check-overlay --conf waterline_station.conf
    python3 waterline_station_check.py pair | profile | trial | pointing ...   (not implemented yet)

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

OUTPUT. Subcommands that report results print JSON on stdout (setup reads it)
and a summary for people on stderr. Exit status: 0 ok, 1 refused or failed,
2 usage (including a subcommand not implemented yet).

Python 3.8 (CACO-05 and the NUCs: Ubuntu 20.04/22.04); standard library only.
"""

import argparse
import json
import math
import re
import sys
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
    """True/False when the zone can be looked up here (zoneinfo, else dateutil), None when it cannot."""
    try:
        import zoneinfo                                  # Python 3.9+
    except ImportError:
        zoneinfo = None
    if zoneinfo is not None:
        try:
            zoneinfo.ZoneInfo(name)
            return True
        except Exception:                                # ZoneInfoNotFoundError, ValueError on odd names
            return False
    try:
        from dateutil import tz                          # Python 3.8: pandas brings dateutil
    except ImportError:
        return None
    return tz.gettz(name) is not None


def check_value(key, v):
    """The problem with value v of `key`, or None. Only the key's own format (see check_overlay)."""
    if key == "WATERLINE_STAGE":
        return None if v in ("", "archive") else 'must be "archive" (stage 1) or ""'
    if key == "CAMS":
        cams = v.split()
        if cams not in (["c1"], ["c1", "c2"]):
            return 'must be "c1" or "c1 c2" (the cron maps c1 and, on a two-camera station, c2)'
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
            notes.append(f"STATION_TZ={value}: not looked up (no zoneinfo or dateutil here)")

    if len(export_lines) > 1:
        err(export_lines[1], None, f"more than one export line (lines {', '.join(map(str, export_lines))})")

    stage = "archive" if keys.get("WATERLINE_STAGE") == "archive" else "full"

    # the cameras and their calibration names
    if "CAMS" in keys and check_value("CAMS", keys["CAMS"]) is None:
        two = "c2" in keys["CAMS"].split()
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


# --- Subcommands still to come (design 6.4-6.6; correction 7). Each will replace its entry in
# --- NOT_YET with its own add_<name>_parser(sub) and cmd_<name>(args), print JSON on stdout and a
# --- summary on stderr, and exit 0 ok / 1 refused or failed / 2 usage.
NOT_YET = (
    ("pair", "which IO, EO and target files belong together (design 6.4)"),
    ("profile", "the detector's search band from a survey, a DEM or the camera geometry (design 6.5)"),
    ("trial", "a trial run of the detector and the georectification on the kept pictures (design 6.6)"),
    ("pointing", "a calibration-day picture against a recent timex: has the camera turned?"),
)


def cmd_not_yet(args):
    say(f"waterline_station_check.py {args.command}: not implemented yet")
    return 2


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

    for name, what in NOT_YET:
        p = sub.add_parser(name, help=f"{what} -- not implemented yet")
        p.set_defaults(func=cmd_not_yet)
    return ap


def main(argv=None):
    ap = build_parser()
    args, extra = ap.parse_known_args(argv)
    if extra and args.func is not cmd_not_yet:          # a stub says "not implemented yet" whatever it is given
        ap.error("unrecognized arguments: " + " ".join(extra))
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

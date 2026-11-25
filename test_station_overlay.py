#!/usr/bin/env python3
"""
Tests For The Station Overlay Of The Cron (waterline_station.conf)
====================================================================
Plain python3 (no pytest needed on the station):

    python3 test_station_overlay.py [--keep DIR] [--baseline REV|FILE]

WHAT. On a station other than CACO-05, coastcam-station-setup writes
waterline_station.conf beside waterline_timex_cron.sh; the cron has
waterline_station_check.py check it, sources it after its own settings, and
takes from it the cameras, calibration names, detector options and switches.
CACO-05 has no such file, so there the cron must run exactly as before:
  a. without a station file the cron runs the same commands, writes the same
     log and the same files as the cron before the overlay (REV below, from
     git), apart from station_status.py's --buoy-csv/--gauge-csv (the same
     files as its defaults); names set in the job's environment change nothing;
  b. with a station file, for one and for two cameras: no owg.sh, no runup, no
     pointing log, the station's DETECTOR_ARGS, CAL and calibration names, its
     GNSS-R file, gauge, buoy, DEM_MAX_HS and setup; the exported settings
     reach the Python scripts; the pointing log, when on, uses the station's EO;
  c. a stage-1 file (WATERLINE_STAGE="archive") keeps the timex pictures and
     runs no Python but the check;
  d. the five readers of the cron's VAR= lines (survey_products.py,
     recalibrate_timing.py, storm_watch.py, station_status.py,
     dem_from_contours.py) read the same values from the cron as before, and
     the cron's column-0 NAME= lines are the same lines;
  e. check-overlay accepts setup's files and refuses a bad one (an unquoted
     value, $, a backtick, a backslash, a quote, a control character, an
     unknown or forbidden key, a bad number, a bad export line, a missing
     station fact ...); the cron then exits 1 with a FATAL line in its log and
     runs nothing else;
  f. a fresh folder without logs/ works (the check runs after the cron makes
     logs/ and defines log(), design critique C1);
  g. the weekly frames step (7c, weekly_frames.py daily --base BASE --have-lock,
     docs/WEEKLY_FRAMES.md): called once, after the pointing log and before
     station_status.py, only while DEM_ENABLE=1; a failure is the cron's
     "WARNING: weekly frames failed (see above)" and the run goes on;
  h. the extra timex (:15/:45, /mnt/I2Rgus_Data/extra_timex/YYYY-MM-DD/; the
     cron's EXTRA TIMEX): with that folder, a second detector call right after
     the first, the same options plus --source-dir, its source folder holding
     only the window's (3 days) timex pictures with no detection yet; their
     CSVs archived with cp -pn, the pictures kept in archive/images_timex; a
     frame already detected skipped, a later run with nothing new makes no
     call; the station's DETECTOR_ARGS on it too; stage 1 keeps the pictures;
     --extra-timex-backfill takes every day, one call per day, and nothing
     else; any other argument is refused. Without the folder (parts a-g) the
     cron is as before.

HOW. The cron is copied into a temporary folder ($BASE, with ImageProducts,
the GNSS-R file and the extra timex folder there too) and run by bash with a
stand-in python3 that records each call (and runs waterline_station_check.py for real) and a
stand-in owg.sh. REV (default BASELINE below: the commit before the overlay)
gives the cron "as before"; update it when the cron changes on purpose. Without
git or REV, part (a)'s comparison and part (d) are skipped with a SKIP line.
The cron changed on purpose after REV in one way only: step 7c, the weekly
frames call (WEEKLY_CALL below). Part (a) expects exactly that one call more,
just before station_status.py, and otherwise the same commands, log and files.

STATION FILES. None are read. The calibration names are created as empty
files in the temporary folder; no folder of the station is listed.
"""

import argparse
import ast
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import waterline_station_check as wsc   # noqa: E402

CRON = HERE / "waterline_timex_cron.sh"
CHECK = HERE / "waterline_station_check.py"
BASELINE = "0b5e98ed44f96fab55687b12d8a1008b0b4b5e12"      # origin/main before the overlay (#104)
CACO05_CAL = ("CACO05_c1_20240801_IO.yaml", "CACO05_c1_20251113_EO-CV.yaml",
              "CACO05_c2_20240801_IO.yaml", "CACO05_c2_20251113_EO-CV.yaml")
STATION_CAL = ("CACO07_c1_20261201_IO.yaml", "CACO07_c1_20261201_EO.yaml",
               "CACO07_c2_20261201_IO.yaml", "CACO07_c2_20261201_EO.yaml")
TIMEX = ("1792500000.Fri.Oct.20_14_00_00.GMT.2026.CACO07.c1.timex.jpg",
         "1792500000.Fri.Oct.20_14_00_00.GMT.2026.CACO07.c2.timex.jpg")
failures = []


def check(name, cond, detail=""):
    print(f"  {'ok ' if cond else 'FAIL'} {name}" + (f"  ({detail})" if detail != "" and not cond else ""))
    if not cond:
        failures.append(name)


# ---------------------------------------------------------------------
# The harness
# ---------------------------------------------------------------------

STUB = r'''#!/bin/bash
echo "PY $*" >> "$CALLS"
case "$1" in
    */waterline_station_check.py) exec "$REAL_PY" "$@" ;;
esac
echo "ENV ${1##*/} ANT=${GNSSR_ANTENNA_NAVD88_M-unset} TZ=${STATION_TZ-unset} EPSG=${GEOTIFF_EPSG-unset} NAME=${STATION_NAME-unset} BUOY=${WAVE_BUOY-unset} PROV=${STATION_PROVISIONAL-unset}" >> "$ENVS"
case "$1" in
    */weekly_frames.py)
        [ -n "$WEEKLY_SAY" ] && echo "$WEEKLY_SAY"
        exit "${WEEKLY_RC:-0}" ;;
esac
prev=""
for a in "$@"; do
    [ "$prev" = "--output" ] && echo "h" > "$a"
    prev="$a"
done
case "$1" in */georectify.py) printf 'h\nr\n' > "$3" ;; esac
exit 0
'''
OWG_STUB = '#!/bin/bash\necho "OWG $*" >> "$CALLS"\nexit 0\n'
# the one call the cron makes after REV on purpose (step 7c), just before station_status.py
WEEKLY_CALL = "PY /B/weekly_frames.py daily --base /B --have-lock"


def with_weekly(calls):
    """REV's calls with the weekly frames call put in just before station_status.py (exactly one more call)."""
    out, put = [], False
    for c in calls:
        if c.startswith("PY /B/station_status.py ") and not put:
            out.append(WEEKLY_CALL)
            put = True
        out.append(c)
    return out


def localise(text, b):
    """The cron's text with $BASE, ImageProducts, the extra timex folder and the GNSS-R file in folder b."""
    n_base = len(re.findall(r"^BASE=/mnt/I2Rgus_Data/waterline$", text, re.M))
    n_spl = len(re.findall(r"^GNSSR_SPLINE=.*$", text, re.M))
    if n_base != 1 or n_spl != 1:
        raise SystemExit(f"the cron's BASE= / GNSSR_SPLINE= lines not found once ({n_base}, {n_spl})")
    text = re.sub(r"^BASE=/mnt/I2Rgus_Data/waterline$", f"BASE={b}", text, flags=re.M)
    text = re.sub(r"^GNSSR_SPLINE=.*$", f"GNSSR_SPLINE={b}/gnss/usgs_spline_out.txt", text, flags=re.M)
    # the extra timex folder too, so a station that has one (CACO-05 once it is installed) never shows it to
    # a run (parts a-g have none; part h makes one in b)
    return text.replace("/mnt/I2Rgus_Data/ImageProducts", f"{b}/ImageProducts").replace(
        "/mnt/I2Rgus_Data/extra_timex", f"{b}/extra_timex")


def make_base(root, name, cron_text, caco05_cal=True, station_cal=()):
    """A station folder: no logs/ (the cron makes it), stand-ins, pictures, a GNSS-R file 2 days old."""
    b = Path(root) / name
    for sub in ("bin", "ImageProducts/products", "gnss"):
        (b / sub).mkdir(parents=True)
    if caco05_cal:
        (b / "calibration").mkdir()
        for f in CACO05_CAL:
            (b / "calibration" / f).write_text("")
    if station_cal:
        (b / "station_calibration").mkdir()
        for f in station_cal:
            (b / "station_calibration" / f).write_text("")
    for f in TIMEX:
        (b / "ImageProducts" / f).write_bytes(b"jpg")
    (b / "ImageProducts" / TIMEX[1].replace("timex", "bright")).write_bytes(b"jpg")
    (b / "ImageProducts/products/1792500000.c2.ras.tiff").write_bytes(b"ras")
    day = subprocess.run(["date", "-u", "-d", "2 days ago", "+%Y %m %d"], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout.strip()
    spline = f"% Hortho 5.0\n60000 10.0 {day} 12 00 00 0.5\n"
    for s in ("usgs_spline_out.txt", "cac7_spline_out.txt"):
        (b / "gnss" / s).write_text(spline)
    (b / "bin/python3").write_text(STUB)
    (b / "bin/python3").chmod(0o755)
    (b / "owg.sh").write_text(OWG_STUB)
    (b / "owg.sh").chmod(0o755)
    shutil.copy(str(CHECK), str(b / CHECK.name))
    shutil.copy(str(HERE / "station_env.py"), str(b / "station_env.py"))   # the check's time zone rule
    (b / "waterline_timex_cron.sh").write_text(localise(cron_text, b))
    return b


class Run:
    def __init__(self, b, rc, out):
        self.b, self.rc, self.out = b, rc, out
        norm = lambda t: t.replace(str(b), "/B")        # noqa: E731
        rd = lambda p: p.read_text() if p.exists() else ""   # noqa: E731
        self.calls = [norm(x) for x in rd(b / "calls.txt").splitlines()]
        self.envs = [norm(x) for x in rd(b / "envs.txt").splitlines()]
        self.log_raw = norm(rd(b / "logs/timex_cron.log"))
        self.log = re.sub(r"^\[\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ\] ", "", self.log_raw, flags=re.M)

    def py(self, script):
        return [c for c in self.calls if c.startswith("PY ") and script in c.split(" ")[1]]

    def files(self):
        """{relative path: bytes} of every file the run left, but the harness's own and the log."""
        out = {}
        skip = {"calls.txt", "envs.txt", "logs/timex_cron.log", "waterline_timex_cron.sh"}
        for p in sorted(self.b.rglob("*")):
            rel = str(p.relative_to(self.b))
            if p.is_file() and rel not in skip:
                out[rel] = p.read_bytes().replace(str(self.b).encode(), b"/B")
            elif p.is_dir():
                out[rel + "/"] = b""
        return out


def run_cron(b, extra_env=None, args=()):
    env = {"PATH": f"{b}/bin:{os.environ.get('PATH', '/usr/bin:/bin')}", "HOME": str(b), "LC_ALL": "C",
           "CALLS": str(b / "calls.txt"), "ENVS": str(b / "envs.txt"), "REAL_PY": sys.executable}
    if os.environ.get("TZ"):
        env["TZ"] = os.environ["TZ"]
    env.update(extra_env or {})
    r = subprocess.run(["bash", str(b / "waterline_timex_cron.sh")] + list(args), stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, env=env, cwd=str(b))
    return Run(b, r.returncode, r.stdout)


def station_conf(b, cams="c1 c2", stage=None, drop=(), extra=None, **over):
    """setup's full station file (design 7.1, with RAS_SOURCE), paths in folder b; or the stage-1 file."""
    if stage == "archive":
        kv = {"WATERLINE_STAGE": "archive", "OWG_ENABLE": "0", "ARCHIVE_PRODUCTS": "", "STATION_NAME": "CACO07"}
        export = "export STATION_NAME"
    else:
        two = "c2" in cams.split()
        kv = {"CAMS": cams, "CAL": f"{b}/station_calibration",
              "CAL_IO_C1": STATION_CAL[0], "CAL_EO_C1": STATION_CAL[1],
              "CAL_IO_C2": STATION_CAL[2] if two else "", "CAL_EO_C2": STATION_CAL[3] if two else "",
              "GNSSR_SPLINE": f"{b}/gnss/cac7_spline_out.txt", "TIDE_GAUGE": "8447930",
              "GAUGE_CSV": f"{b}/archive/gauge_8447930.csv", "WAVE_BUOY": "44097",
              "WAVES_CSV": f"{b}/archive/waves_44097.csv", "DEM_MAX_HS": "1.0", "SETUP_COEF": "",
              "DETECTOR_ARGS": f"--no-bias-correction --profile-json {b}/station_calibration/detector_profile.json",
              "OWG_ENABLE": "0", "RUNUP_ENABLE": "0", "POINTING_LOG_ENABLE": "0", "ARCHIVE_PRODUCTS": "",
              "RAS_SOURCE": "", "GNSSR_ANTENNA_NAVD88_M": "12.345", "STATION_TZ": "America/New_York",
              "GEOTIFF_EPSG": "6348", "STATION_NAME": "CACO07",
              "STATION_PROVISIONAL": "search band from camera geometry; no wave setup (C = 0)"}
        export = ("export GNSSR_ANTENNA_NAVD88_M STATION_TZ GEOTIFF_EPSG STATION_NAME STATION_PROVISIONAL "
                  "STATION_WAVE_CURRENCY GNSSR_QC_GROSS_M WAVE_BUOY")
    kv.update(over)
    lines = ["# waterline_station.conf -- CACO07's settings for waterline_timex_cron.sh.",
             "# Written by coastcam-station-setup (./setup.sh --waterline); not tracked by git.",
             "# Do not edit: change station.conf and run ./setup.sh --waterline again."]
    lines += [f'{k}="{v}"' for k, v in kv.items() if k not in drop]
    lines += list(extra or [])
    lines.append(export)
    return "\n".join(lines) + "\n"


def baseline_text(rev):
    """The cron before the overlay: git revision `rev` of it, or the file `rev`; None when neither is there."""
    if os.path.isfile(rev):
        return Path(rev).read_text()
    try:
        r = subprocess.run(["git", "-C", str(HERE), "show", f"{rev}:waterline_timex_cron.sh"],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    except OSError:
        return None
    return r.stdout if r.returncode == 0 and r.stdout.strip() else None


# ---------------------------------------------------------------------
# a. no station file: as before
# ---------------------------------------------------------------------

def test_no_station_file(d, old):
    print("a. no station file: the same commands, log and files as the cron before the overlay")
    new = run_cron(make_base(d, "a_new", CRON.read_text()))
    check("the run ends cleanly", new.rc == 0 and "=== run end ===" in new.log, new.out[-300:])
    check("no check-overlay call, owg.sh run, both cameras mapped, CACO-05's calibration and its pointing log",
          not new.py("waterline_station_check.py") and "OWG --no-fetch --no-report" in new.calls
          and len(new.py("daily_elevation_map.py")) == 4 and len(new.py("pointing_daily.py")) == 1
          and len(new.py("runup_from_timestack.py")) == 2, new.calls)
    geo = new.py("georectify.py")
    check("georectify: CACO-05's four calibration files",
          len(geo) == 1 and geo[0].endswith("--io-c1 /B/calibration/CACO05_c1_20240801_IO.yaml "
                                            "--eo-c1 /B/calibration/CACO05_c1_20251113_EO-CV.yaml "
                                            "--io-c2 /B/calibration/CACO05_c2_20240801_IO.yaml "
                                            "--eo-c2 /B/calibration/CACO05_c2_20251113_EO-CV.yaml"), geo)
    st = new.py("station_status.py")
    check("station_status.py: CACO-05's buoy and gauge files (its own defaults)",
          st == ["PY /B/station_status.py --base /B --gnssr /B/gnss/usgs_spline_out.txt "
                 "--buoy-csv /B/archive/waves_44008.csv --gauge-csv /B/archive/gauge_8447435.csv"], st)
    check("no station setting reaches a Python script",
          new.envs and all("ANT=unset TZ=unset EPSG=unset NAME=unset BUOY=unset PROV=unset" in e
                           for e in new.envs), new.envs[:2])
    if old is None:
        print("  SKIP the comparison with the cron before the overlay (no git, or the revision is not here)")
        return
    gold = run_cron(make_base(d, "a_old", old))
    want = with_weekly([c + " --buoy-csv /B/archive/waves_44008.csv --gauge-csv /B/archive/gauge_8447435.csv"
                        if c.startswith("PY /B/station_status.py ") else c for c in gold.calls])
    check("the golden run (the cron before the overlay) ran its chain",
          gold.rc == 0 and len(gold.calls) > 15 and "OWG --no-fetch --no-report" in gold.calls, gold.calls)
    check("the golden run has no weekly frames call; the new one exactly one, just before station_status.py",
          not [c for c in gold.calls if "weekly_frames.py" in c] and len(want) == len(gold.calls) + 1
          and [c for c in new.calls if "weekly_frames.py" in c] == [WEEKLY_CALL]
          and new.calls.index(WEEKLY_CALL) + 1 == next(i for i, c in enumerate(new.calls)
                                                       if c.startswith("PY /B/station_status.py ")),
          [c for c in new.calls if "weekly_frames.py" in c or "station_status.py" in c])
    check("the same commands in the same order (station_status.py: + --buoy-csv/--gauge-csv; + the weekly frames "
          "call before it)", new.calls == want, "\n".join(x for x in set(new.calls) ^ set(want)))
    check("the same log, line for line (timestamps aside)", new.log == gold.log,
          "\n".join(set(new.log.splitlines()) ^ set(gold.log.splitlines())))
    check("the same files, with the same contents", new.files() == gold.files(),
          sorted(set(new.files().items()) ^ set(gold.files().items()))[:6])
    check("the same exit status", new.rc == gold.rc)
    env = {"WATERLINE_STAGE": "archive", "CAMS": "c1", "OWG_ENABLE": "0", "DETECTOR_ARGS": "--no-bias-correction",
           "CAL_IO_C1": "x.yaml", "CAL_EO_C1": "y.yaml", "CAL_IO_C2": "", "CAL_EO_C2": ""}
    leak = run_cron(make_base(d, "a_env", CRON.read_text()), env)
    check("the same names in the job's environment change nothing (the cron unsets them)",
          leak.calls == want and leak.log == gold.log, "\n".join(set(leak.calls) ^ set(want)))


# ---------------------------------------------------------------------
# b. a station file
# ---------------------------------------------------------------------

def test_station_file(d):
    print("b. a station file: the station's cameras, calibration, detector options, sources and switches")
    for cams in ("c1 c2", "c1"):
        two = cams == "c1 c2"
        tag = "two cameras" if two else "one camera"
        b = make_base(d, "b_" + cams.replace(" ", ""), CRON.read_text(), caco05_cal=False,
                      station_cal=STATION_CAL if two else STATION_CAL[:2])
        (b / "waterline_station.conf").write_text(station_conf(b, cams))
        r = run_cron(b)
        check(f"{tag}: the run ends cleanly, the check logged OK, no ERROR or FATAL",
              r.rc == 0 and "=== run end ===" in r.log and "check-overlay: waterline_station.conf OK" in r.log
              and "ERROR" not in r.log and "FATAL" not in r.log, r.log[-600:])
        check(f"{tag}: the file checked first, then the detector",
              r.calls[:2] == ["PY /B/waterline_station_check.py check-overlay --conf /B/waterline_station.conf",
                              "PY waterline_detector_v5.py --image-suffix timex.jpg --timing-mode legacy "
                              "--input-dir /B/input_timex --output-dir /B/processed_timex --debug-dir /B/debug "
                              "--no-bias-correction --profile-json /B/station_calibration/detector_profile.json"],
              r.calls[:2])
        check(f"{tag}: no owg.sh, no runup, no pointing log (OWG/RUNUP/POINTING_LOG_ENABLE=0)",
              not any(c.startswith("OWG") for c in r.calls) and not r.py("runup_from_timestack.py")
              and not r.py("pointing_daily.py") and "runup stage disabled (RUNUP_ENABLE=0)" in r.log, r.calls)
        geo = r.py("georectify.py")
        want_geo = ("--io-c1 /B/station_calibration/CACO07_c1_20261201_IO.yaml "
                    "--eo-c1 /B/station_calibration/CACO07_c1_20261201_EO.yaml")
        if two:
            want_geo += (" --io-c2 /B/station_calibration/CACO07_c2_20261201_IO.yaml "
                         "--eo-c2 /B/station_calibration/CACO07_c2_20261201_EO.yaml")
        check(f"{tag}: georectify with the station's CAL and calibration names only",
              len(geo) == 1 and geo[0] == "PY /B/georectify.py /B/contour_points_timex_qc.csv "
                                          "/B/contour_points_ground.csv " + want_geo, geo)
        maps = [c.split()[4] for c in r.py("daily_elevation_map.py")]
        check(f"{tag}: maps of CAMS only", maps == (["c1", "c1", "c2", "c2"] if two else ["c1", "c1"]), maps)
        check(f"{tag}: the station's buoy, gauge and GNSS-R file; its DEM_MAX_HS; no setup (SETUP_COEF=\"\")",
              r.py("fetch_buoy_waves.py") == ["PY /B/fetch_buoy_waves.py --station 44097 "
                                              "--output /B/archive/waves_44097.csv"]
              and r.py("fetch_tide_gauge.py") == ["PY /B/fetch_tide_gauge.py --station 8447930 --days 120 "
                                                  "--output /B/archive/gauge_8447930.csv"]
              and r.py("extract_elevation_contours.py") == [
                  "PY /B/extract_elevation_contours.py /B/gnss/cac7_spline_out.txt --processed-dir "
                  "/B/archive/processed_timex --waves /B/archive/waves_44097.csv --gnssr-qc-reference "
                  "/B/archive/gauge_8447930.csv "
                  "--timing-mode legacy --output /B/contour_points_timex.csv"]
              and all("--max-hs 1.0 " in c for c in r.py("dem_from_contours.py"))
              and r.py("station_status.py") == ["PY /B/station_status.py --base /B --gnssr "
                                                "/B/gnss/cac7_spline_out.txt --buoy-csv /B/archive/waves_44097.csv "
                                                "--gauge-csv /B/archive/gauge_8447930.csv"],
              r.calls)
        ext = [e for e in r.envs if e.startswith("ENV extract_elevation_contours.py ")]
        check(f"{tag}: the exported settings reach the Python scripts",
              ext == ["ENV extract_elevation_contours.py ANT=12.345 TZ=America/New_York EPSG=6348 NAME=CACO07 "
                      "BUOY=44097 PROV=search band from camera geometry; no wave setup (C = 0)"], ext)
        check(f"{tag}: no timestacks (RAS_SOURCE=\"\"), no other products (ARCHIVE_PRODUCTS=\"\")",
              not list((b / "archive/ras_c2").iterdir()) and not (b / "archive/images_bright").exists()
              and len(list((b / "archive/images_timex").iterdir())) == 2, sorted(os.listdir(b / "archive")))
        check(f"{tag}: CACO-05's calibration folder never needed", not (b / "calibration").exists())
    check('RAS_SOURCE="" copies nothing: the copy is guarded (not a glob of /*.c2.ras.tiff)',
          '[ -n "$RAS_SOURCE" ] && cp -pn "$RAS_SOURCE"/*.c2.ras.tiff "$ARCHIVE_RAS"/ 2>/dev/null' in CRON.read_text())

    # the pointing log on, one camera: the station's EO, c1 only
    b = make_base(d, "b_point", CRON.read_text(), caco05_cal=False, station_cal=STATION_CAL[:2])
    (b / "waterline_station.conf").write_text(station_conf(b, "c1", POINTING_LOG_ENABLE="1", RUNUP_ENABLE="1"))
    r = run_cron(b)
    check("POINTING_LOG_ENABLE=1 on one camera: pointing_daily.py with the station's c1 EO only",
          r.py("pointing_daily.py") == ["PY /B/pointing_daily.py --image-dir /B/archive/images_timex --log-dir "
                                        "/B/archive --eo-c1 /B/station_calibration/CACO07_c1_20261201_EO.yaml"],
          r.py("pointing_daily.py"))
    check("RUNUP_ENABLE=1 from the station file: the runup stage runs", len(r.py("runup_from_timestack.py")) == 2)

    # a station calibration file missing: the cron's own ERROR, not a silent CACO-05 file
    b = make_base(d, "b_nocal", CRON.read_text(), caco05_cal=True, station_cal=STATION_CAL[1:])
    (b / "waterline_station.conf").write_text(station_conf(b, "c1 c2"))
    r = run_cron(b)
    check("the station's c1 IO missing: 'calibration not found', no georectification (CACO-05's files unused)",
          "ERROR: calibration not found in /B/station_calibration" in r.log and not r.py("georectify.py"),
          r.log[-400:])
    b = make_base(d, "b_nocal2", CRON.read_text(), caco05_cal=False, station_cal=STATION_CAL[:3])
    (b / "waterline_station.conf").write_text(station_conf(b, "c1 c2"))
    r = run_cron(b)
    check("two cameras, c2's IO there but not its EO: georectify is still called (georectify.py reports it)",
          len(r.py("georectify.py")) == 1)
    b = make_base(d, "b_nocal3", CRON.read_text(), caco05_cal=False, station_cal=STATION_CAL[:2])
    (b / "waterline_station.conf").write_text(station_conf(b, "c1 c2"))
    r = run_cron(b)
    check("two cameras, c2's IO missing: 'calibration not found'",
          "ERROR: calibration not found" in r.log and not r.py("georectify.py"))


# ---------------------------------------------------------------------
# c. stage 1
# ---------------------------------------------------------------------

def test_stage1(d):
    print("c. stage 1 (WATERLINE_STAGE=\"archive\"): the timex kept, no Python but the check")
    b = make_base(d, "c_stage1", CRON.read_text(), caco05_cal=False)
    (b / "waterline_station.conf").write_text(station_conf(b, stage="archive"))
    r = run_cron(b)
    kept = sorted(os.listdir(b / "archive/images_timex"))
    check("exit 0, the run logged start to end", r.rc == 0 and "=== run start ===" in r.log
          and r.log.rstrip().endswith("=== run end ==="), r.log)
    check("the only Python call is the check; no owg.sh", r.calls == [
        "PY /B/waterline_station_check.py check-overlay --conf /B/waterline_station.conf"], r.calls)
    check("the timex pictures kept in archive/images_timex", kept == sorted(TIMEX), kept)
    check("logged: stage 1, the number kept", "stage 1 (WATERLINE_STAGE=archive): timex kept: 2; no detection, "
          "maps or DEM until setup's continue step" in r.log, r.log)
    check("nothing else archived or written (no detections, products, timestacks, status)",
          not (b / "archive/images_bright").exists() and not list((b / "archive/ras_c2").iterdir())
          and not (b / "station_status.txt").exists() and "detections archived" not in r.log, r.log)
    (b / "ImageProducts" / TIMEX[0]).write_bytes(b"other")
    (b / "ImageProducts" / TIMEX[0].replace("1792500000", "1792503600")).write_bytes(b"jpg")
    r2 = run_cron(b)
    check("a second run adds the new picture, never overwrites one kept",
          r2.rc == 0 and len(os.listdir(b / "archive/images_timex")) == 3
          and (b / "archive/images_timex" / TIMEX[0]).read_bytes() == b"jpg" and "timex kept: 3" in r2.log)


# ---------------------------------------------------------------------
# d. the readers of the cron's VAR= lines
# ---------------------------------------------------------------------

def test_readers(d, old):
    print("d. the five readers of the cron's VAR= lines: the same values as before")
    new = CRON.read_text()
    col0 = lambda t: [ln for ln in t.splitlines() if re.match(r"^[A-Za-z_][A-Za-z0-9_]*=", ln)]   # noqa: E731
    check('no new NAME= line: the defaults are : "${NAME=default}" lines',
          'unset WATERLINE_STAGE CAMS OWG_ENABLE DETECTOR_ARGS CAL_IO_C1 CAL_EO_C1 CAL_IO_C2 CAL_EO_C2' in new
          and ': "${WATERLINE_STAGE=}" "${CAMS=c1 c2}" "${OWG_ENABLE=1}" "${DETECTOR_ARGS=}"' in new
          and ': "${CAL_IO_C1=CACO05_c1_20240801_IO.yaml}" "${CAL_EO_C1=CACO05_c1_20251113_EO-CV.yaml}"' in new
          and ': "${CAL_IO_C2=CACO05_c2_20240801_IO.yaml}" "${CAL_EO_C2=CACO05_c2_20251113_EO-CV.yaml}"' in new
          and not [ln for ln in col0(new) if ln.split("=")[0] in ("WATERLINE_STAGE", "CAMS", "OWG_ENABLE",
                                                                   "DETECTOR_ARGS") or ln.startswith("CAL_")])
    if old is None:
        print("  SKIP the readers' comparison (no git, or the revision is not here)")
        return
    check("every column-0 NAME= line is the line it was, in the same order", col0(new) == col0(old),
          "\n".join(set(col0(new)) ^ set(col0(old))))
    import survey_products as sp
    import recalibrate_timing as rt
    import storm_watch as sw
    import station_status as ss
    import dem_from_contours as dc
    got = {}
    for name, text in (("old", old), ("new", new)):
        f = Path(d) / f"d_{name}"
        f.mkdir()
        (f / "waterline_timex_cron.sh").write_text(text)
        p = f / "waterline_timex_cron.sh"
        nobase = Path(d) / "d_base"
        nobase.mkdir(exist_ok=True)
        got[name] = {"survey_products.cron_settings": sp.cron_settings(p),
                     "recalibrate_timing.cron_settings": rt.cron_settings(p),
                     "storm_watch.cron_settings": sw.cron_settings(nobase, cron_path=p),
                     "storm_watch.var_lines": sw.var_lines(text),
                     "station_status.cron_timing_mode": ss.cron_timing_mode(f),
                     "dem_from_contours.station_setup_coef": dc.station_setup_coef(p)}
    for k in got["old"]:
        check(f"{k}: the same", got["old"][k] == got["new"][k], (got["old"][k], got["new"][k]))
    g = got["new"]
    check("... and CACO-05's values (not an empty read)",
          g["survey_products.cron_settings"].get("SETUP_COEF") == "0.037"
          and g["survey_products.cron_settings"].get("WAVE_BUOY") == "44008"
          and g["recalibrate_timing.cron_settings"]["TIDE_GAUGE"] == "8447435"
          and g["storm_watch.cron_settings"][0]["DEM_CELL"] == "2.0" and not g["storm_watch.cron_settings"][1]
          and g["station_status.cron_timing_mode"] == "legacy"
          and g["dem_from_contours.station_setup_coef"] == 0.037, g)


# ---------------------------------------------------------------------
# e. check-overlay
# ---------------------------------------------------------------------

def test_check_overlay(d):
    print("e. check-overlay: setup's files accepted, bad ones refused; the cron stops on a refusal")
    b = "/mnt/I2Rgus_Data/waterline"
    good = {
        "setup's full file, two cameras (STATION_PROVISIONAL with ; ( ) =)": station_conf(b),
        "one camera": station_conf(b, "c1"),
        "the stage-1 file": station_conf(b, stage="archive"),
        "the stage-1 file without its export line (stage 1 runs no Python: a note)":
            station_conf(b, stage="archive").replace("export STATION_NAME\n", ""),
        "the design's example (no RAS_SOURCE), comments with $ and quotes, blank lines, no final newline":
            station_conf(b, drop=("RAS_SOURCE",), extra=["", "# a comment may say $BASE or `x` or \"q\"",
                                                          "   ", "\t# indented comment"]).rstrip("\n"),
        "optional settings set": station_conf(b, extra=['GNSSR_QC_GROSS_M="-3.0,4.5"',
                                                        'STATION_WAVE_CURRENCY="buoy:hs/tp"',
                                                        'ARCHIVE_PRODUCTS="bright dark"']).replace(
                                                            'ARCHIVE_PRODUCTS=""\n', ""),
        "the stage-1 file with WATERLINE_STAGE first and more settings":
            station_conf(b, stage="archive", CAMS="c1", CAL_IO_C2="", CAL_EO_C2=""),
        "CAMS with extra plain spaces (bash splits it into c1 and c2 all the same)": station_conf(b, CAMS=" c1  c2"),
    }
    for what, text in good.items():
        res = wsc.check_overlay(text)
        check(f"accepted: {what}", res["ok"], res["errors"])
    full = station_conf(b)
    bad = [
        ("an unquoted value", full.replace('CAMS="c1 c2"', "CAMS=c1"), "double quotes"),
        ("single quotes", full.replace('CAMS="c1 c2"', "CAMS='c1 c2'"), "double quotes"),
        ("$ inside the quotes", full.replace(f'CAL="{b}/station_calibration"', 'CAL="$BASE/station_calibration"'),
         "$ inside"),
        ("$( ) inside the quotes", full.replace('STATION_NAME="CACO07"', 'STATION_NAME="$(id)"'), "$ inside"),
        ("a backtick", full.replace('STATION_NAME="CACO07"', 'STATION_NAME="`id`"'), "backtick"),
        ("a backslash", full.replace('STATION_NAME="CACO07"', 'STATION_NAME="CA\\CO"'), "backslash"),
        ("a quote inside", full.replace('STATION_NAME="CACO07"', 'STATION_NAME="CA"CO"'), 'a "'),
        ("a tab inside", full.replace('STATION_NAME="CACO07"', 'STATION_NAME="CA\tCO"'), "control"),
        ("a command after the value", full.replace('CAMS="c1 c2"', 'CAMS="c1 c2"; touch /tmp/x'), "nothing after"),
        ("a comment after the value", full.replace('CAMS="c1 c2"', 'CAMS="c1 c2"  # two'), "nothing after"),
        ("spaces around =", full.replace('CAMS="c1 c2"', 'CAMS = "c1 c2"'), "no spaces"),
        ("an indented setting", full.replace('CAMS="c1 c2"', '  CAMS="c1 c2"'), "no spaces"),
        ("a command line", full + "rm -rf /tmp/x\n", "not a comment"),
        ("a source line", full + ". /tmp/other.conf\n", "not a comment"),
        ("Windows line ends", full.replace("\n", "\r\n"), "carriage return"),
        ("an unknown key", full + 'FOO="1"\n', "not a name"),
        ("a lower-case key", full + 'cams="c1"\n', "not a name"),
        ("forbidden: TIMING_MODE", full + 'TIMING_MODE="mid"\n', "may not be set"),
        ("forbidden: DEM_CELL", full + 'DEM_CELL="1.0"\n', "may not be set"),
        ("forbidden: BASE", full + 'BASE="/tmp"\n', "may not be set"),
        ("forbidden: LOG", full + 'LOG="/tmp/x.log"\n', "may not be set"),
        ("forbidden: PATH", full + 'PATH="/tmp"\n', "may not be set"),
        ("a key set twice", full + 'CAMS="c1 c2"\n', "set twice"),
        ("a bad number: DEM_MAX_HS 1,5", full.replace('DEM_MAX_HS="1.0"', 'DEM_MAX_HS="1,5"'), "DEM_MAX_HS"),
        ("a bad number: antenna height with a unit", full.replace('"12.345"', '"12.3m"'), "GNSSR_ANTENNA"),
        ("a bad number: antenna height empty", full.replace('"12.345"', '""'), "GNSSR_ANTENNA"),
        ("a bad number: nan", full.replace('"12.345"', '"nan"'), "GNSSR_ANTENNA"),
        ("a bad number: SETUP_COEF", full.replace('SETUP_COEF=""', 'SETUP_COEF="0.0.37"'), "SETUP_COEF"),
        ("a bad switch: OWG_ENABLE yes", full.replace('OWG_ENABLE="0"', 'OWG_ENABLE="yes"'), "OWG_ENABLE"),
        ("a bad EPSG", full.replace('GEOTIFF_EPSG="6348"', 'GEOTIFF_EPSG="EPSG:6348"'), "GEOTIFF_EPSG"),
        ("a bad gauge", full.replace('TIDE_GAUGE="8447930"', 'TIDE_GAUGE="844793"'), "TIDE_GAUGE"),
        ("a bad gross range", full + 'GNSSR_QC_GROSS_M="4,-3"\n', "GNSSR_QC_GROSS_M"),
        ("an unknown time zone", full.replace("America/New_York", "America/Nowhere_Town"), "STATION_TZ"),
        ("a bad stage", full + 'WATERLINE_STAGE="stage1"\n', "WATERLINE_STAGE"),
        ("cameras c2 alone", full.replace('CAMS="c1 c2"', 'CAMS="c2"'), "CAMS"),
        ("CAMS with a no-break space (bash keeps 'c1\u00a0c2' one word)",
         full.replace('CAMS="c1 c2"', 'CAMS="c1\u00a0c2"'), "not plain ASCII"),
        ("CAMS with an em space", full.replace('CAMS="c1 c2"', 'CAMS="c1\u2003c2"'), "CAMS"),
        ("CAMS with an ideographic space", full.replace('CAMS="c1 c2"', 'CAMS="c1\u3000c2"'), "CAMS"),
        ("a relative path", full.replace(f'GAUGE_CSV="{b}/', 'GAUGE_CSV="archive/'), "absolute path"),
        ("a path with a space", full.replace(f'WAVES_CSV="{b}/archive/', f'WAVES_CSV="{b}/my archive/'),
         "absolute path"),
        ("a calibration name with a folder", full.replace('CAL_IO_C1="CACO07', 'CAL_IO_C1="../CACO07'), "CAL_IO_C1"),
        ("DETECTOR_ARGS moving the detector's output folder",
         full.replace("--no-bias-correction", "--output-dir /tmp/x"), "--output-dir"),
        ("DETECTOR_ARGS with an abbreviation of --input-dir", full.replace("--no-bias-correction", "--inp=/x"),
         "--input-dir"),
        ("DETECTOR_ARGS with --help", full.replace("--no-bias-correction", "--help"), "--help"),
        ("DETECTOR_ARGS with a wildcard", full.replace("--no-bias-correction", "--profile-json /x/*.json"),
         "wildcards"),
        ("one camera with c2's calibration still named", station_conf(b, "c1", CAL_IO_C2=STATION_CAL[2],
                                                                        CAL_EO_C2=STATION_CAL[3]), "one camera"),
        ("one camera without CAL_IO_C2/CAL_EO_C2", station_conf(b, "c1", drop=("CAL_IO_C2", "CAL_EO_C2")),
         "one camera"),
        ("two cameras with c2's IO empty", station_conf(b, CAL_IO_C2=""), "maps c2"),
        ("c2's IO named, its EO empty", station_conf(b, "c1", CAL_IO_C2=STATION_CAL[2]), "CAL_IO_C2"),
        ("a stage-2 file without the antenna height", station_conf(b, drop=("GNSSR_ANTENNA_NAVD88_M",)),
         "missing"),
        ("a stage-2 file without SETUP_COEF (CACO-05's 0.037 would apply)", station_conf(b, drop=("SETUP_COEF",)),
         "missing"),
        ("an empty file (CACO-05's facts all round)", "", "missing"),
        ("a stage-2 file whose antenna height is not exported",
         full.replace("export GNSSR_ANTENNA_NAVD88_M ", "export "), "not exported"),
        ("an export line exporting PATH", full.replace("export GNSSR", "export PATH GNSSR"), "export PATH"),
        ("an export line with an assignment", full.replace("export GNSSR_ANTENNA_NAVD88_M ",
                                                           "export GNSSR_ANTENNA_NAVD88_M=1 "), "export line"),
        ("two export lines", full + "export STATION_NAME\n", "more than one export"),
        ("export -n", full.replace("export GNSSR", "export -n GNSSR"), "export line"),
    ]
    for what, text, why in bad:
        res = wsc.check_overlay(text)
        msgs = " | ".join(e["message"] for e in res["errors"])
        check(f"refused: {what}", not res["ok"] and why in msgs, msgs or "accepted")

    # CAMS is split as the cron's unquoted $CAMS is after it sources the file (bash, default IFS): a value
    # check-overlay takes for c1 and c2 is c1 and c2 in the map loop, never one word with a Unicode space
    for v in ("c1 c2", "c1", " c1  c2", "c1\u00a0c2", "c1\u2003c2", "c1\u3000c2"):
        conf = Path(d) / "e_cams.conf"
        conf.write_bytes(f'CAMS="{v}"\n'.encode("utf-8"))
        r = subprocess.run(["bash", "-c", '. "$1"; for cam in $CAMS; do printf "[%s]" "$cam"; done', "x", str(conf)],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        words = re.findall(r"\[(.*?)\]", r.stdout.decode("utf-8"))
        accepted = wsc.check_value("CAMS", v) is None
        check(f"CAMS={v!r}: split as bash splits it ({words}); accepted only when that is c1 or c1 c2",
              r.returncode == 0 and wsc._shell_words(v) == words and accepted == (words in (["c1"], ["c1", "c2"])),
              (words, wsc._shell_words(v), accepted))

    # STATION_TZ: check-overlay asks station_env's own lookup, so a zone it accepts never stops the contour
    # and DEM steps at import (station_env.station_tz(), run by the cron's same python3)
    import station_env
    zones = ("America/New_York", "Europe/Lisbon", "America/Los_Angeles", "EST5", "CST6", "UTC+3", "GMT-5",
             "localtime", "posixrules", "Factory", "America/Nowhere_Town")
    saved = os.environ.get("STATION_TZ")
    try:
        for v in zones:
            os.environ["STATION_TZ"] = v
            try:
                station_env.station_tz()
                steps_take_it = True
            except SystemExit:
                steps_take_it = False
            gate = wsc.check_value("STATION_TZ", v) is None
            check(f"STATION_TZ={v}: check-overlay {'accepts' if gate else 'refuses'} it, as the Python steps "
                  f"{'do' if steps_take_it else 'refuse it'} (the same lookup: {wsc._known_zone(v)})",
                  gate == steps_take_it and wsc._known_zone(v) == station_env._zone_known(v),
                  (gate, steps_take_it, wsc._known_zone(v), station_env._zone_known(v)))
    finally:
        if saved is None:
            os.environ.pop("STATION_TZ", None)
        else:
            os.environ["STATION_TZ"] = saved
    # ... and where zoneinfo is missing (Python 3.8), with pytz there: pytz decides for both (dateutil alone
    # would take EST5 or localtime, which the steps then refuse on every run)
    probe = ("import sys, json\nsys.modules['zoneinfo'] = None\nsys.path.insert(0, %r)\n"
             "import waterline_station_check as w, station_env as s\n"
             "try:\n    import pytz\n    have = True\nexcept ImportError:\n    have = False\n"
             "print(json.dumps([have, {v: [w.check_value('STATION_TZ', v) is None, s._zone_known(v)] for v in %r}]))"
             % (str(HERE), list(zones)))
    r = subprocess.run([sys.executable, "-c", probe], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       universal_newlines=True)
    try:
        have_pytz, got = json.loads(r.stdout)
    except ValueError:
        have_pytz, got = None, {}
    check("without zoneinfo: check-overlay and station_env give the same answer for every zone",
          r.returncode == 0 and got and all(a == (k is not False) for a, k in got.values()), (r.stdout, r.stderr[-400:]))
    if have_pytz:
        check("without zoneinfo, with pytz: EST5, UTC+3 and localtime refused by check-overlay (pytz does not know "
              "them; the steps would stop at import)", got and not any(got[v][0] for v in ("EST5", "UTC+3", "localtime"))
              and got["America/New_York"][0] and got["Europe/Lisbon"][0], got)
    else:
        print("  SKIP without zoneinfo, with pytz: pytz does not import here (the Python 3.8 run checks it)")

    # the command line: JSON on stdout, the summary on stderr, exit 0 / 1 / 2
    f = Path(d) / "e_cli.conf"

    def cli(*args):
        r = subprocess.run([sys.executable, str(CHECK)] + list(args), stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, universal_newlines=True)
        return r.returncode, r.stdout, r.stderr

    rc, out, err = cli("capabilities")
    check("capabilities prints exactly 'overlay-1 storm-alerts-1 station-setup-1'",
          rc == 0 and out == "overlay-1 storm-alerts-1 station-setup-1\n" and err == "", (rc, out, err))
    f.write_text(full)
    rc, out, err = cli("check-overlay", "--conf", str(f))
    try:
        js = json.loads(out)
    except ValueError:
        js = {}
    check("check-overlay OK: exit 0, JSON on stdout (ok, stage, keys), the summary on stderr",
          rc == 0 and js.get("ok") is True and js.get("stage") == "full" and "CAMS" in js.get("keys", [])
          and js.get("errors") == [] and "e_cli.conf OK" in err and "{" not in err, (rc, out[:200], err))
    f.write_text(full + 'FOO="1"\n')
    rc, out, err = cli("check-overlay", "--conf", str(f))
    js = json.loads(out) if out.startswith("{") else {}
    check("check-overlay refused: exit 1, the problem in the JSON and on stderr with its line",
          rc == 1 and js.get("ok") is False and js["errors"][0]["key"] == "FOO" and "REFUSED" in err
          and f"line {len(full.splitlines()) + 1}: FOO is not a name" in err, (rc, out[:300], err))
    rc, out, err = cli("check-overlay", "--conf", str(Path(d) / "no_such.conf"))
    check("a file that cannot be read: exit 1", rc == 1 and json.loads(out)["ok"] is False, (rc, err))
    f.write_bytes(b'STATION_NAME="\xff"\n')
    rc, out, err = cli("check-overlay", "--conf", str(f))
    check("a file that is not UTF-8: exit 1", rc == 1 and "not UTF-8" in err, (rc, err))
    rc, out, err = cli("check-overlay")
    check("no --conf: usage, exit 2", rc == 2, (rc, err))
    # pair, profile, trial and pointing are real commands (test_waterline_station_check.py tests what they
    # do): here only that each has its own parser (a usage error naming its required options, exit 2,
    # nothing on stdout) and that an input file that is not there is a usage error too, said in the JSON
    nx = str(Path(d) / "no_such_file")
    not_there = {
        "pair": ("--xyz", nx, "--uv", nx, "--io", nx, "--eo", nx),
        "profile": ("--camera", "c1", "--io", nx, "--eo", nx, "--image-size", "2448x2048", "--geometry",
                    "--output", nx + ".json"),
        "trial": ("--camera", "c1", "--io", nx, "--eo", nx, "--profile-json", nx, "--spline", nx, "--frames", nx),
        "pointing": ("--io", nx, "--eo", nx, "--calibration-picture", nx, "--timex", nx),
    }
    for sub in ("pair", "profile", "trial", "pointing"):
        rc, out, err = cli(sub, "--anything", "x")
        check(f"{sub} with an unknown option and none of its own: its usage, the required options named, exit 2, "
              f"nothing on stdout", rc == 2 and out == "" and f"usage: waterline_station_check.py {sub}" in err
              and f"waterline_station_check.py {sub}: error: the following arguments are required" in err,
              (rc, out[:200], err))
        rc, out, err = cli(sub, *not_there[sub])
        try:
            js = json.loads(out)
        except ValueError:
            js = {}
        check(f"{sub} with input files that are not there: usage, exit 2, 'not found' in the JSON and on stderr",
              rc == 2 and js.get("command") == sub and js.get("ok") is False and "not found" in js.get("error", "")
              and f"{sub}: not found" in err and not Path(nx + ".json").exists(), (rc, out[:300], err))
    rc, out, err = cli()
    check("no subcommand: usage, exit 2", rc == 2)

    # the cron: a refused file stops it with a FATAL line, nothing else run
    for what, text, why in (("an unquoted value", full.replace('CAMS="c1 c2"', "CAMS=c1"), "double quotes"),
                            ("$", full.replace('STATION_NAME="CACO07"', 'STATION_NAME="$(touch pwned)"'), "$ inside"),
                            ("a backtick", full.replace('STATION_NAME="CACO07"', 'STATION_NAME="`touch pwned`"'),
                             "backtick"),
                            ("an unknown key", full + 'FOO="1"\n', "FOO is not a name"),
                            ("a forbidden key", full + 'TIMING_MODE="mid"\n', "TIMING_MODE may not be set"),
                            ("a bad number", full.replace('DEM_MAX_HS="1.0"', 'DEM_MAX_HS="1,0"'), "DEM_MAX_HS"),
                            ("a command line", full + "touch pwned\n", "not a comment")):
        bb = make_base(d, "e_cron_" + re.sub(r"\W+", "_", what), CRON.read_text())
        (bb / "waterline_station.conf").write_text(text.replace(b, str(bb)))
        r = run_cron(bb)
        check(f"the cron with {what}: exit 1, the reason and a FATAL line in its log, nothing run",
              r.rc == 1 and "FATAL: /B/waterline_station.conf refused (see above)" in r.log and why in r.log
              and "REFUSED" in r.log and r.calls == ["PY /B/waterline_station_check.py check-overlay "
                                                     "--conf /B/waterline_station.conf"]
              and "=== run start ===" not in r.log and not (bb / "pwned").exists(), (r.rc, r.log, r.calls))


# ---------------------------------------------------------------------
# f. a fresh folder
# ---------------------------------------------------------------------

def test_fresh_base(d):
    print("f. a fresh folder without logs/ (critique C1): the check logs once log() and logs/ exist")
    text = CRON.read_text()
    m, h, lk = (text.index('mkdir -p "$SCRATCH"'), text.index('if [ -f "$BASE/waterline_station.conf" ]'),
                text.index('exec 9>"$LOCK"'))
    check("the hook comes after the mkdir of logs/ and the log() definition, before the lock",
          m < text.index("log() {") < h < lk)
    s1, det, cd = (text.index('if [ "$WATERLINE_STAGE" = "archive" ]'), text.index("# 1. Detect"),
                   text.index('cd "$BASE" ||'))
    check("the stage-1 block comes after the lock and the cd, just before step 1", lk < cd < s1 < det)
    r = subprocess.run(["bash", "-n", str(CRON)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True)
    check("bash -n", r.returncode == 0, r.stdout)
    for what, conf, rc_want, line in (("a good station file", station_conf, 0, "=== run end ==="),
                                      ("a refused station file", lambda bb: station_conf(bb) + 'FOO="1"\n', 1,
                                       "FATAL: /B/waterline_station.conf refused")):
        b = make_base(d, "f_" + str(rc_want), text, caco05_cal=False, station_cal=STATION_CAL)
        (b / "waterline_station.conf").write_text(conf(b))
        before = (b / "logs").exists()
        r = run_cron(b)
        check(f"{what}: no logs/ before; the run makes it, exit {rc_want}, the check's line and '{line}' in the log",
              not before and r.rc == rc_want and "check-overlay: waterline_station.conf" in r.log
              and line in r.log, (r.rc, r.log[-500:], r.out[-300:]))


# ---------------------------------------------------------------------
# g. the weekly frames step (7c)
# ---------------------------------------------------------------------

def test_weekly_step(d):
    print("g. the weekly frames step: once, after the pointing log, before station_status.py, only with DEM_ENABLE=1")
    text = CRON.read_text()
    call = 'python3 "$BASE/weekly_frames.py" daily --base "$BASE" --have-lock >> "$LOG" 2>&1'
    i7b, i7c, i8 = (text.find("# 7b. Camera pointing log"), text.find(call), text.find("# 8. Station health summary"))
    guard = text.rfind('if [ "$DEM_ENABLE" = "1" ]; then', 0, i7c)
    check("the cron calls weekly_frames.py once, with --have-lock, between step 7b and step 8, inside "
          "'if [ \"$DEM_ENABLE\" = \"1\" ]' at column 0, with its WARNING line",
          text.count('"$BASE/weekly_frames.py"') == 1 and text.count(call) == 1 and 0 < i7b < guard < i7c < i8
          and text[guard - 1] == "\n" and "fi" not in text[guard:i7c].split()
          and text[i7c:i8].count('|| log "WARNING: weekly frames failed (see above)"') == 1, (i7b, guard, i7c, i8))
    b = make_base(d, "g_ok", text)
    r = run_cron(b, {"WEEKLY_SAY": "weekly frames: done today already (2026-11-23 UTC)"})
    wk = r.py("weekly_frames.py")
    pos = [i for i, c in enumerate(r.calls) if "pointing_daily.py" in c or "weekly_frames.py" in c
           or "station_status.py" in c]
    check("DEM_ENABLE=1: one call, 'daily --base BASE --have-lock', after pointing_daily.py and before "
          "station_status.py; its own line in the log, no WARNING for it",
          wk == [WEEKLY_CALL] and [r.calls[i].split()[1].split("/")[-1] for i in pos]
          == ["pointing_daily.py", "weekly_frames.py", "station_status.py"]
          and "weekly frames: done today already (2026-11-23 UTC)" in r.log
          and "weekly frames failed" not in r.log and r.rc == 0, (wk, r.log[-400:]))
    b = make_base(d, "g_fail", text)
    r = run_cron(b, {"WEEKLY_SAY": "WARNING: weekly frames: week 2026-11-09..2026-11-15: dem_from_contours.py failed",
                     "WEEKLY_RC": "1"})
    lines = r.log.splitlines()
    k = next((i for i, ln in enumerate(lines) if ln.startswith("WARNING: weekly frames: week")), -1)
    check("weekly_frames.py exit 1: its line, then the cron's 'WARNING: weekly frames failed (see above)'; "
          "station_status.py still runs and the run ends cleanly",
          k >= 0 and lines[k + 1] == "WARNING: weekly frames failed (see above)" and r.rc == 0
          and len(r.py("station_status.py")) == 1 and "=== run end ===" in r.log, r.log[-500:])
    off = re.sub(r"^DEM_ENABLE=1$", "DEM_ENABLE=0", text, count=1, flags=re.M)
    b = make_base(d, "g_off", off)
    r = run_cron(b)
    check("DEM_ENABLE=0: no weekly frames call (and the DEM stage logged as disabled)",
          off != text and not r.py("weekly_frames.py") and "DEM stage disabled (DEM_ENABLE=0)" in r.log
          and len(r.py("station_status.py")) == 1, r.calls[-4:])


# ---------------------------------------------------------------------
# h. the extra timex (:15/:45)
# ---------------------------------------------------------------------

# the stand-in python3, but the detector also writes one CSV per timex picture of its source folder (--source-dir,
# else ImageProducts) into --output-dir (cleared first, as the detector clears it) and lists its source folder
DET_STUB = STUB.replace('prev=""\n', r"""case "$1" in
    *waterline_detector_v5.py)
        src="$SRCDEF"; out=""; prev=""
        for a in "$@"; do
            [ "$prev" = "--source-dir" ] && src="$a"
            [ "$prev" = "--output-dir" ] && out="$a"
            prev="$a"
        done
        echo "SRC $src: $(cd "$src" && find . -type f | sed 's|^\./||' | sort | tr '\n' ' ')" >> "$SRCS"
        rm -f "$out"/*
        for f in "$src"/*.timex.jpg; do
            [ -f "$f" ] || continue
            n="${f##*/}"
            echo "det $n" > "$out/${n%.jpg}.csv"
            echo "$(stat -c %Y "$f") $n" >> "$SRCS.mtimes"
        done
        exit "${DET_RC:-0}" ;;
esac
prev=""
""", 1)
DET_CALL = ("PY waterline_detector_v5.py --image-suffix timex.jpg --timing-mode legacy --input-dir /B/input_timex "
            "--output-dir /B/processed_timex --debug-dir /B/debug")
XT_CALL = ("PY waterline_detector_v5.py --image-suffix timex.jpg --timing-mode legacy --source-dir "
           "/B/extra_timex_source --input-dir /B/input_timex --output-dir /B/processed_timex --debug-dir /B/debug")


def local_day(days_ago):
    """The local date days_ago days before today, as the cron's date reads it."""
    return subprocess.run(["date", "-d", f"today -{days_ago} days", "+%F"], stdout=subprocess.PIPE,
                          universal_newlines=True).stdout.strip()


def xt_name(day, cam, minute, hour=14):
    """An extra timex picture's name (Argus style) at <hour>:<minute> local on <day>."""
    e = int(subprocess.run(["date", "-d", f"{day} {hour}:{minute:02d}", "+%s"], stdout=subprocess.PIPE,
                           universal_newlines=True).stdout.strip())
    return f"{e}.Mon.Jan.01_{hour}_{minute:02d}_00.GMT.2026.CACO05.{cam}.timex.jpg"


def make_extra(b, detected=True):
    """b/extra_timex as extra_timex.sh leaves it: dated folders of today and 1, 2, 3 days ago (3 is outside the
    window), a snap, a raw/ run folder and ArgusEyes' products/ folder; one frame (1 day ago, c1) already detected
    when `detected`. -> {name: day} of every timex picture in a dated folder."""
    xt = b / "extra_timex"
    days = [local_day(k) for k in range(4)]
    pics = {}
    for k, cams_minutes in enumerate(([("c1", 15), ("c2", 15)], [("c1", 45), ("c2", 45)], [("c1", 15)],
                                      [("c1", 15)])):
        (xt / days[k]).mkdir(parents=True)
        for cam, m in cams_minutes:
            n = xt_name(days[k], cam, m)
            (xt / days[k] / n).write_bytes(b"xt " + n.encode())
            pics[n] = days[k]
    snap = xt_name(days[0], "c1", 15).replace("timex", "snap")
    (xt / days[0] / snap).write_bytes(b"snap")
    (xt / days[0] / "raw/1792500000000").mkdir(parents=True)
    (xt / days[0] / "raw/1792500000000" / xt_name(days[0], "c1", 16)).write_bytes(b"raw")
    (xt / "products").mkdir()
    (xt / "products" / xt_name(days[0], "c2", 17)).write_bytes(b"in progress")
    (xt / "raw").mkdir()
    (xt / "notes").mkdir()
    (xt / "notes" / xt_name(days[0], "c2", 18)).write_bytes(b"not a dated folder")
    os.utime(str(xt / days[0] / xt_name(days[0], "c1", 15)), (1792500000, 1792500000))
    if detected:
        (b / "archive/processed_timex").mkdir(parents=True)
        (b / "archive/processed_timex" / xt_name(days[1], "c1", 45).replace(".jpg", ".csv")).write_text("old\n")
    return days, pics


def xt_run(b, args=(), **env):
    for f in ("calls.txt", "envs.txt", "srcs.txt", "srcs.txt.mtimes"):
        if (b / f).exists():
            (b / f).unlink()
    r = run_cron(b, dict({"SRCDEF": str(b / "ImageProducts"), "SRCS": str(b / "srcs.txt")}, **env), args)
    r.srcs = [x.replace(str(b), "/B") for x in
              ((b / "srcs.txt").read_text().splitlines() if (b / "srcs.txt").exists() else [])]
    mt = b / "srcs.txt.mtimes"
    r.mtimes = dict(ln.split(" ", 1)[::-1] for ln in (mt.read_text().splitlines() if mt.exists() else []))
    return r


def src_line(folder, names):
    return f"SRC {folder}: " + "".join(n + " " for n in sorted(names))


CP_FAIL_STUB = """#!/bin/bash
# cp, failing (as on a full disk) when an argument holds $CP_FAIL
if [ -n "${CP_FAIL:-}" ]; then
    for a in "$@"; do
        case "$a" in *"$CP_FAIL"*) echo "cp: cannot create regular file '$a': No space left on device" >&2; exit 1 ;; esac
    done
fi
exec %s "$@"
"""


def test_extra_timex(d):
    print("h. the extra timex (:15/:45): a second detector call on the window's undetected pictures, archived as the others")
    # A time zone whose date differs from UTC's right now (UTC-12 before 12:00 UTC, UTC+14 after), so the window
    # is proven to use the LOCAL date (the 20:55 run on the NUC is the next day in UTC): `date -u` fails here.
    saved = os.environ.get("TZ")
    os.environ["TZ"] = "XST12" if time.gmtime().tm_hour < 12 else "XST-14"
    try:
        utc = subprocess.run(["date", "-u", "+%F"], stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
        check(f"part h runs in TZ={os.environ['TZ']}, where today is not UTC's today", local_day(0) != utc,
              (local_day(0), utc))
        extra_timex_cases(d)
    finally:
        if saved is None:
            os.environ.pop("TZ", None)
        else:
            os.environ["TZ"] = saved


def extra_timex_cases(d):
    text = CRON.read_text()
    lines = [ln for ln in text.splitlines() if "/mnt/I2Rgus_Data/extra_timex" in ln]
    check("the extra timex folder is a local of extra_timex_step (no column-0 NAME= line for it)",
          [ln.strip() for ln in lines if not ln.lstrip().startswith("#")]
          == ['local mode="$1" xt=/mnt/I2Rgus_Data/extra_timex cut="" d days=()'], lines)
    b = make_base(d, "h_extra", text)
    (b / "bin/python3").write_text(DET_STUB)
    days, pics = make_extra(b)
    old = xt_name(days[1], "c1", 45)
    new = sorted(n for n, day in pics.items() if day != days[3] and n != old)
    r = xt_run(b)
    dets = r.py("waterline_detector_v5.py")
    check("the run ends cleanly, no WARNING about the extra timex",
          r.rc == 0 and "=== run end ===" in r.log and "WARNING: extra timex" not in r.log, r.log[-600:])
    check("two detector calls, the second right after the first: the same options plus --source-dir "
          "/B/extra_timex_source", dets == [DET_CALL, XT_CALL] and r.calls[:2] == dets, r.calls[:3])
    check("its source folder: only the window's (today and the 2 days before) timex pictures with no detection yet "
          "(not the one detected, not 3 days ago, no snap, nothing from raw/, products/ or a folder not dated)",
          len(r.srcs) == 2 and r.srcs[0].startswith("SRC /B/ImageProducts: ")
          and r.srcs[1] == src_line("/B/extra_timex_source", new) and len(new) == 4, (r.srcs, new))
    arch = b / "archive/processed_timex"
    got = sorted(p.name for p in arch.iterdir())
    want = sorted([n.replace(".jpg", ".csv") for n in list(TIMEX) + new] + [old.replace(".jpg", ".csv")])
    check("their CSVs archived beside the first run's; the frame already detected keeps its CSV (cp -pn); none for "
          "the day outside the window", got == want and (arch / old.replace(".jpg", ".csv")).read_text() == "old\n"
          and all((arch / n.replace(".jpg", ".csv")).read_text() == f"det {n}\n" for n in new), got)
    img = b / "archive/images_timex"
    kept = sorted(p.name for p in img.iterdir())
    check("the window's extra timex pictures kept in archive/images_timex (the one detected before too), with "
          "their times (cp -p); not 3 days ago's, no snap",
          kept == sorted(list(TIMEX) + new + [old])
          and int((img / xt_name(days[0], "c1", 15)).stat().st_mtime) == 1792500000, kept)
    check("the detector saw each picture with its own time (cp -p into the source folder: its file-name/file-time "
          "check)", r.mtimes.get(xt_name(days[0], "c1", 15)) == "1792500000", r.mtimes)
    check("logged: the extra run's count, then the archive's images with them; the source folder removed",
          "extra timex detections archived: +4 of 4 new picture(s) (total 7)" in r.log
          and "images in archive: 7" in r.log and not (b / "extra_timex_source").exists(), r.log)
    check("the extra timex folder left as it was (nothing moved or removed)",
          sorted(p.name for p in (b / "extra_timex" / days[0]).iterdir())
          == sorted([xt_name(days[0], "c1", 15), xt_name(days[0], "c2", 15),
                     xt_name(days[0], "c1", 15).replace("timex", "snap"), "raw"]))
    r = xt_run(b)
    check("a second run with nothing new: no second detector call, a line in the log, no WARNING",
          r.py("waterline_detector_v5.py") == [DET_CALL] and r.rc == 0
          and "extra timex: nothing new to detect (3 dated folder(s), every picture detected)" in r.log
          and "WARNING: extra timex" not in r.log and not (b / "extra_timex_source").exists(),
          (r.calls[:3], r.log[-500:]))
    late = xt_name(days[0], "c1", 45)
    seed = xt_name(days[0], "c1", 15)
    seed_csv = (arch / seed.replace(".jpg", ".csv")).read_text()
    (b / "extra_timex" / days[0] / late).write_bytes(b"late")
    r = xt_run(b, DET_RC="1")
    check("a new picture: the second call's source folder holds it, and the same camera's frame 30 min before it "
          "(detected by an earlier run) as the blend's previous frame, whose CSV is kept; the detector failing is "
          "a WARNING and what it wrote is still archived",
          r.srcs[1:] == [src_line("/B/extra_timex_source", [late, seed])]
          and "extra timex detections archived: +1 of 1 new picture(s) (total 8; 1 frame(s) detected before given "
              "as the blend's previous frame, their CSVs kept)" in r.log
          and (arch / seed.replace(".jpg", ".csv")).read_text() == seed_csv
          and "WARNING: extra timex: detector exited 1 (archiving whatever it produced)" in r.log
          and (arch / late.replace(".jpg", ".csv")).exists() and (img / late).exists(), (r.srcs, r.log[-500:]))
    # the blend's previous frame only within 35 min (TEMPORAL_MAX_GAP_MINUTES), of the same camera; an empty
    # picture (an extra run killed while writing it) skipped with a line, never given to the detector
    c1_new, c2_new = xt_name(days[0], "c1", 15, 15), xt_name(days[0], "c2", 15, 15)
    c2_empty = xt_name(days[0], "c2", 45)
    for n, data in ((c1_new, b"new"), (c2_new, b"new"), (c2_empty, b"")):
        (b / "extra_timex" / days[0] / n).write_bytes(data)
    r = xt_run(b)
    last = r.log[r.log.rfind("=== run start ==="):]     # the log holds every run of this base
    check("15:15 c1 gets its 14:45 frame (30 min, detected before); 15:15 c2 none (its last detected frame 60 min "
          "before, the empty 14:45 one never detected); the empty one skipped, logged, no ERROR, no extra timex WARNING",
          r.srcs[1:] == [src_line("/B/extra_timex_source", [c1_new, late, c2_new])]
          and "extra timex: 1 empty picture(s) (0 bytes: an extra run stopped while writing them) skipped, never "
              "detected" in last
          and "extra timex detections archived: +2 of 2 new picture(s) (total 10; 1 frame(s) detected before" in last
          and "ERROR" not in last and "WARNING: extra timex" not in last
          and not (arch / c2_empty.replace(".jpg", ".csv")).exists(), (r.srcs, last))
    r = xt_run(b)
    last = r.log[r.log.rfind("=== run start ==="):]     # the log holds every run of this base
    check("run again: the empty picture skipped again, nothing new to detect, no second call",
          r.py("waterline_detector_v5.py") == [DET_CALL] and "1 empty picture(s)" in last
          and "extra timex: nothing new to detect (3 dated folder(s), every picture detected)" in last, last)
    # a picture that cannot be copied into the source folder (a full disk): a WARNING, never "every picture
    # detected"
    stuck = xt_name(days[0], "c1", 45, 15)
    (b / "extra_timex" / days[0] / stuck).write_bytes(b"new")
    (b / "bin/cp").write_text(CP_FAIL_STUB % shutil.which("cp"))
    (b / "bin/cp").chmod(0o755)
    r = xt_run(b, CP_FAIL=stuck)
    last = r.log[r.log.rfind("=== run start ==="):]
    check("a picture cp cannot copy into the source folder: a WARNING naming it, cp's message in the log, no second "
          "call, no 'every picture detected', the source folder removed",
          r.py("waterline_detector_v5.py") == [DET_CALL]
          and f"WARNING: extra timex: 1 picture(s) could not be copied into /B/extra_timex_source (e.g. "
              f"/B/extra_timex/{days[0]}/{stuck}; see above): not detected this run" in last
          and "No space left on device" in last and "nothing new to detect" not in last
          and not (b / "extra_timex_source").exists(), last)
    (b / "bin/cp").unlink()

    # the station's DETECTOR_ARGS on the second call too
    b = make_base(d, "h_station", text, caco05_cal=False, station_cal=STATION_CAL)
    (b / "bin/python3").write_text(DET_STUB)
    (b / "waterline_station.conf").write_text(station_conf(b))
    make_extra(b)
    r = xt_run(b)
    args = " --no-bias-correction --profile-json /B/station_calibration/detector_profile.json"
    check("with a station file: both detector calls carry its DETECTOR_ARGS",
          r.py("waterline_detector_v5.py") == [DET_CALL + args, XT_CALL + args] and r.rc == 0, r.calls[:3])

    # stage 1: the pictures kept, no detection
    b = make_base(d, "h_stage1", text, caco05_cal=False)
    (b / "bin/python3").write_text(DET_STUB)
    (b / "waterline_station.conf").write_text(station_conf(b, stage="archive"))
    days, pics = make_extra(b, detected=False)
    r = xt_run(b)
    win = [n for n, day in pics.items() if day != days[3]]
    check("stage 1: the only Python call is the check; the window's extra timex kept with the others",
          r.calls == ["PY /B/waterline_station_check.py check-overlay --conf /B/waterline_station.conf"]
          and sorted(os.listdir(b / "archive/images_timex")) == sorted(list(TIMEX) + win)
          and "timex kept: 7;" in r.log, (r.calls, r.log))
    r = xt_run(b, ["--extra-timex-backfill"])
    check("stage 1 with --extra-timex-backfill: every day's extra timex kept, still nothing detected",
          r.calls == ["PY /B/waterline_station_check.py check-overlay --conf /B/waterline_station.conf"]
          and sorted(os.listdir(b / "archive/images_timex")) == sorted(list(TIMEX) + list(pics))
          and "timex kept: 8;" in r.log, (r.calls, r.log))

    # the backfill: every day, one detector call per day, oldest first, nothing else run
    b = make_base(d, "h_backfill", text)
    (b / "bin/python3").write_text(DET_STUB)
    days, pics = make_extra(b)
    old = xt_name(days[1], "c1", 45)
    r = xt_run(b, ["--extra-timex-backfill"])
    bf = XT_CALL.replace("--debug-dir /B/debug", "--debug-dir /B/extra_timex_backfill_debug")
    by_day = [sorted(n for n, day in pics.items() if day == days[k] and n != old) for k in (3, 2, 1, 0)]
    check("--extra-timex-backfill: one detector call per dated folder, oldest first, each on that day's undetected "
          "pictures, its debug images in a scratch folder removed after; no other Python call, no owg.sh",
          r.calls == [bf] * 4 and r.srcs == [src_line("/B/extra_timex_source", x) for x in by_day]
          and not (b / "extra_timex_backfill_debug").exists() and not (b / "extra_timex_source").exists()
          and r.rc == 0, (r.calls, r.srcs))
    check("... every day's CSVs and pictures archived (3 days ago's too), the one detected before kept; logged",
          sorted(p.name for p in (b / "archive/processed_timex").iterdir())
          == sorted(n.replace(".jpg", ".csv") for n in pics)
          and (b / "archive/processed_timex" / old.replace(".jpg", ".csv")).read_text() == "old\n"
          and sorted(os.listdir(b / "archive/images_timex")) == sorted(pics)
          and "extra timex backfill: 4 dated folder(s) in /B/extra_timex" in r.log
          and f"extra timex backfill: {days[3]}" in r.log and "extra timex backfill: done;" in r.log
          and r.log.rstrip().endswith("=== run end ===") and "detections archived: +" not in
          r.log.replace("extra timex detections archived: +", ""), r.log)
    r = xt_run(b, ["--extra-timex-backfill"])
    check("run again: nothing to detect, no call", r.calls == [] and r.rc == 0
          and r.log.count("extra timex: nothing new to detect") >= 4, (r.calls, r.log[-400:]))
    r = xt_run(make_base(d, "h_backfill_none", text), ["--extra-timex-backfill"])
    check("--extra-timex-backfill without the extra timex folder: says so, nothing run",
          r.calls == [] and r.rc == 0 and "extra timex backfill: no /B/extra_timex (extra timex not installed): "
          "nothing to do" in r.log, r.log)
    bb = make_base(d, "h_badarg", text)
    r = xt_run(bb, ["--extra-timex-backfil"])
    check("any other argument: usage, exit 2, nothing run, nothing written",
          r.rc == 2 and "usage:" in r.out and r.calls == [] and not (bb / "logs").exists(), (r.rc, r.out))

    # the folder there, no dated folder in the window: a line, no call
    b = make_base(d, "h_old_only", text)
    (b / "bin/python3").write_text(DET_STUB)
    days, pics = make_extra(b)
    for k in (0, 1, 2):
        shutil.rmtree(str(b / "extra_timex" / days[k]))
    r = xt_run(b)
    check("only folders older than the window: one detector call, a line in the log, no WARNING, nothing kept",
          r.py("waterline_detector_v5.py") == [DET_CALL]
          and "extra timex: no dated folder of the last 3 days in /B/extra_timex: nothing to detect" in r.log
          and "WARNING: extra timex" not in r.log
          and sorted(os.listdir(b / "archive/images_timex")) == sorted(TIMEX), r.log[-500:])


# ---------------------------------------------------------------------
# Python 3.8
# ---------------------------------------------------------------------

def test_py38():
    print("Python 3.8 syntax")
    for f in ("waterline_station_check.py", "test_station_overlay.py"):
        try:
            ast.parse((HERE / f).read_text(), filename=f, feature_version=(3, 8))
            ok = True
        except SyntaxError as e:
            ok = e
        check(f"{f} parses as Python 3.8", ok is True, ok)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="keep the run folders, in a new folder inside this one")
    ap.add_argument("--baseline", default=BASELINE,
                    help=f"the cron before the overlay: a git revision (default {BASELINE[:7]}) or a copy of the file")
    args = ap.parse_args()
    if args.keep:
        Path(args.keep).mkdir(parents=True, exist_ok=True)
    d = os.path.realpath(tempfile.mkdtemp(prefix="station_overlay_", dir=args.keep))
    print(f"test folder: {d}")
    old = baseline_text(args.baseline)
    try:
        test_no_station_file(d, old)
        test_station_file(d)
        test_stage1(d)
        test_readers(d, old)
        test_check_overlay(d)
        test_fresh_base(d)
        test_weekly_step(d)
        test_extra_timex(d)
        test_py38()
    finally:
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    if failures:
        print(f"\n{len(failures)} FAILED:")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    print("\nall station overlay tests passed")


if __name__ == "__main__":
    main()

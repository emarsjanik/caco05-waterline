#!/usr/bin/env python3
"""
Tests: a new station's settings from the environment (station_env.py)
=======================================================================
A new station's waterline_station.conf (coastcam-station-setup writes it;
CACO-05 has none) exports GNSSR_ANTENNA_NAVD88_M, STATION_TZ,
STATION_WAVE_CURRENCY, STATION_NAME, STATION_PROVISIONAL, GEOTIFF_EPSG,
GNSSR_QC_GROSS_M and WAVE_BUOY to the cron's Python steps. This checks:

  a. unset: every default is CACO-05's value (19.014, America/New_York,
     marconi:hs_best/tp_s, 6348, -2.5..3.5, 44008, the CACO-05 titles), and
     asc_to_geotiff.py still says 'default=6348' (test_survey_products.py
     reads it there);
  b. each one takes effect where it is read: extract_elevation_contours.py
     (the antenna height on the levels, local days, the wave currency),
     dem_from_contours.py (local days, the station file's SETUP_COEF,
     STATION_PROVISIONAL), gnssr_qc.py (gross range, waves_<buoy>.csv),
     asc_to_geotiff.py (the GeoTIFF's EPSG), dem_figure.py / dem_map.py /
     dem_change.py (titles, frame labels);
  c. a value that does not parse stops the script with an ERROR naming the
     variable (exit 1, no traceback), never CACO-05's value in its place;
     above all a bad or empty antenna height;
  d. the raw-buoy line of extract_elevation_contours.py: a WARNING with
     --setup-coef, a NOTE without one; none when the station's C was
     fitted on its raw buoy (STATION_WAVE_CURRENCY);
  e. a one-camera chain on a synthetic c1 fixture (rows from a planar beach
     through the CACO05 c1 2025-11-13 calibration, named explicitly):
     extract -> georectify (c1 only) -> DEM with its page -> map, GeoTIFF,
     change. With the variables unset and with CACO-05's values set
     explicitly, every file is byte-identical; with a new station's values
     the antenna height reaches the levels and the DEM, info.json and the
     page carry STATION_PROVISIONAL, the page and map the station's name,
     the GeoTIFF and the axes the station's frame.

    python3 test_station_env.py [--keep DIR] [--before REF]

--before REF also runs the chain (variables unset) with the seven edited
scripts as of git revision REF (e.g. the commit before station_env.py) and
requires every output file byte-identical to today's: the no-change check
of the edit itself. Not run by default: later, unrelated changes to those
scripts would make it fail.

Python 3.8 compatible; no network; no pytest. Calibration files are named
one by one (calibration/CACO05_c1_20240801_IO.yaml and
calibration/CACO05_c1_20251113_EO-CV.yaml); the folder is never listed.
"""

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
CAL_IO_C1 = HERE / "calibration" / "CACO05_c1_20240801_IO.yaml"
CAL_EO_C1 = HERE / "calibration" / "CACO05_c1_20251113_EO-CV.yaml"

STATION_VARS = ("GNSSR_ANTENNA_NAVD88_M", "STATION_TZ", "STATION_WAVE_CURRENCY", "STATION_NAME",
                "STATION_PROVISIONAL", "GEOTIFF_EPSG", "GNSSR_QC_GROSS_M", "WAVE_BUOY")
# CACO-05's values, set explicitly: everything must come out as with none set (no STATION_NAME or
# STATION_PROVISIONAL: CACO-05 has neither, and its two pages have different built-in titles)
CACO05_EXPLICIT = {"GNSSR_ANTENNA_NAVD88_M": "19.014", "STATION_TZ": "America/New_York",
                   "STATION_WAVE_CURRENCY": "marconi:hs_best/tp_s", "GEOTIFF_EPSG": "6348",
                   "GNSSR_QC_GROSS_M": "-2.5,3.5", "WAVE_BUOY": "44008"}
NEW_STATION = {"GNSSR_ANTENNA_NAVD88_M": "20.014", "STATION_TZ": "America/Los_Angeles",
               "STATION_WAVE_CURRENCY": "buoy:wvht_m/dpd_s", "GEOTIFF_EPSG": "6339",
               "GNSSR_QC_GROSS_M": "-3.0,4.5", "WAVE_BUOY": "46029", "STATION_NAME": "CACO07 Test Beach",
               "STATION_PROVISIONAL": "search band from camera geometry; no wave setup (C = 0)"}
EDITED = ("extract_elevation_contours.py", "dem_from_contours.py", "gnssr_qc.py", "asc_to_geotiff.py",
          "dem_figure.py", "dem_map.py", "dem_change.py")

FAILED, SKIPPED = [], []


def check(cond, msg, detail=None):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILED.append(msg)
        if detail is not None:
            print("       " + str(detail)[-1500:].replace("\n", "\n       "))


def skip(msg):
    print("  skip " + msg)
    SKIPPED.append(msg)


def clean_env(extra=None, pythonpath=None):
    """This process's environment without any station variable, plus `extra`."""
    e = {k: v for k, v in os.environ.items() if k not in STATION_VARS and k not in ("TIMING_MODE", "CACO_LOCAL_WAVES")}
    e["PYTHONPATH"] = str(pythonpath or HERE)
    e["MPLBACKEND"] = "Agg"
    e.update(extra or {})
    return e


def run(args, cwd, env=None, pythonpath=None):
    """-> (exit code, stdout, stderr) of python3 args (sys.executable)."""
    r = subprocess.run([sys.executable] + [str(a) for a in args], stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True, cwd=str(cwd),
                       env=clean_env(env, pythonpath))
    return r.returncode, r.stdout, r.stderr


def py(code, env=None):
    """Runs python -c code with the station variables `env` only; -> (rc, stdout, stderr)."""
    return run(["-c", code], HERE, env)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def folder_digest(d):
    return {p.name: sha(p) for p in sorted(Path(d).iterdir()) if p.is_file()}


# ---------------------------------------------------------------------
# a. defaults
# ---------------------------------------------------------------------

def test_defaults():
    print("a. unset: CACO-05's values everywhere")
    code = ("import json, station_env as s, extract_elevation_contours as x, dem_from_contours as d, gnssr_qc as g, "
            "dem_figure as f\n"
            "print(json.dumps([x.GNSSR_ANTENNA_NAVD88_M, x.LOCAL_TZ, x.WAVE_CURRENCY_FITTED, d.LOCAL_TZ, "
            "[g.GROSS_MIN, g.GROSS_MAX], f.STATION, s.geotiff_epsg(), s.wave_buoy(), s.station_provisional(), "
            "s.frame_label(), s.frame_label(short=True), "
            "str(g.default_waves_path('/nonexistent/gauge.csv')), "
            "x.wave_currency_note(x.WAVE_CURRENCY_FITTED)[1]]))")
    rc, out, err = py(code)
    got = json.loads(out) if rc == 0 else None
    check(rc == 0 and got[:12] == [19.014, "America/New_York", "marconi:hs_best/tp_s", "America/New_York",
                                   [-2.5, 3.5], "CACO05 Marconi Beach", 6348, "44008", None,
                                   "NAD83(2011) / UTM 19N, EPSG:6348", "UTM 19N", "None"],
          "antenna 19.014, America/New_York (both scripts), marconi:hs_best/tp_s, gross -2.5..3.5, the page's "
          "'CACO05 Marconi Beach', EPSG 6348 (NAD83(2011) / UTM 19N), buoy 44008, nothing provisional",
          (rc, out, err))
    check(rc == 0 and got[12] == "marconi:hs_best/tp_s -- the currency the setup coefficient was fitted in "
                                 "(hs_best: ADCP wh_4061 Hs; tp_s: NDBC 44008 peak period)",
          "the fitted-currency note is CACO-05's text, word for word", got and got[12])
    check(sum(1 for line in (HERE / "asc_to_geotiff.py").read_text().splitlines()
              if 'add_argument("--epsg", type=int, default=6348,' in line) == 1,
          "asc_to_geotiff.py still declares --epsg default=6348 (test_survey_products.py reads it there)")
    # the station variables set to "" (an empty export) are unset, except the antenna height
    rc2, out2, err2 = py(code, {k: "" for k in STATION_VARS if k != "GNSSR_ANTENNA_NAVD88_M"})
    check(rc2 == 0 and out2 == out, "every other variable set to \"\": the same as unset", (rc2, out2, err2))


# ---------------------------------------------------------------------
# b. each one takes effect; c. bad values refused
# ---------------------------------------------------------------------

def refused(code, name, value, label):
    rc, out, err = py(code, {name: value})
    check(rc == 1 and f"ERROR: {name}=" in err and "Traceback" not in err and "unset " + name in err,
          f"{name}={value!r} ({label}): refused, exit 1, an ERROR naming it, no traceback", (rc, out, err))


def write_spline(path, epochs, levels, hortho=19.014):
    """gnssrefl subdaily spline_out.txt layout, Hortho in the header (as test_timing_mode.py)."""
    with open(path, "w") as f:
        f.write("% station test gnssrefl v4.1.5\n")
        f.write(f"% orthometric height minus RH, where Hortho (m) is    {hortho:.3f} \n")
        f.write("% MJD RH(m) YYYY MM DD HH MM SS quasi-sea-level(m)\n")
        for t, z in zip(epochs, levels):
            u = datetime.fromtimestamp(int(t), tz=timezone.utc)
            f.write(f"{t / 86400.0 + 40587.0:15.7f} {hortho - z:10.3f} {u.year:4d} {u.month:3d} {u.day:3d} "
                    f"{u.hour:3d} {u.minute:3d} {u.second:3d} {z:10.3f}\n")


def test_antenna(d):
    print("b/c. GNSSR_ANTENNA_NAVD88_M: the levels follow it; a bad or empty value is refused")
    t0 = int(datetime(2026, 10, 1, 12, 0, tzinfo=timezone.utc).timestamp())
    ep = np.arange(t0, t0 + 6 * 3600 + 1, 1800)
    write_spline(d / "spline.txt", ep, np.round(0.5 * np.sin(np.arange(len(ep)) / 3.0), 3))
    code = ("import json, extract_elevation_contours as x\n"
            f"e, lv, _, h = x.load_gnssr_spline({str(d / 'spline.txt')!r})\n"
            "print(json.dumps([x.GNSSR_ANTENNA_NAVD88_M, h, [round(float(v), 4) for v in lv]]))")
    rc0, out0, err0 = py(code)
    rc1, out1, err1 = py(code, {"GNSSR_ANTENNA_NAVD88_M": "12.345"})
    rc2, out2, _ = py(code, {"GNSSR_ANTENNA_NAVD88_M": "19.014"})
    a, b = (json.loads(out0.splitlines()[-1]), json.loads(out1.splitlines()[-1])) if rc0 == rc1 == 0 else (None, None)
    check(a is not None and a[0] == 19.014 and a[1] == 19.014 and "already NAVD88" in out0,
          "unset: 19.014, the file's Hortho 19.014: no correction", (out0, err0))
    check(b is not None and b[0] == 12.345 and all(abs(y - x - (12.345 - 19.014)) < 1e-6 for x, y in zip(a[2], b[2]))
          and "-6.669 m added to every level" in out1,
          "GNSSR_ANTENNA_NAVD88_M=12.345: every level moves by 12.345 - 19.014 = -6.669 m", (out1, err1))
    check(rc2 == 0 and out2 == out0, "GNSSR_ANTENNA_NAVD88_M=19.014 (CACO-05's, explicit): as unset")
    rc3, out3, err3 = py("import extract_elevation_contours as x; print(x.GNSSR_ANTENNA_NAVD88_M)",
                         {"GNSSR_ANTENNA_NAVD88_M": " 7.5 "})
    check(rc3 == 0 and out3.strip() == "7.5", "surrounding blanks are ignored", (rc3, out3, err3))
    imp = "import extract_elevation_contours"
    for value, label in (("abc", "not a number"), ("", "empty"), ("  ", "blank"), ("12.3m", "with a unit"),
                         ("nan", "nan"), ("inf", "inf"), ("1e9", "out of range"), ("19,014", "decimal comma")):
        refused(imp, "GNSSR_ANTENNA_NAVD88_M", value, label)
    # the scripts that import the constant follow it
    rc4, out4, err4 = py("from extract_elevation_contours import GNSSR_ANTENNA_NAVD88_M as A; print(A)",
                         {"GNSSR_ANTENNA_NAVD88_M": "12.345"})
    check(rc4 == 0 and out4.strip() == "12.345",
          "'from extract_elevation_contours import GNSSR_ANTENNA_NAVD88_M' (compare_gnssr_to_tidemodel.py, "
          "gnssir_reflection_audit.py) gets the station's height", (rc4, out4, err4))


def test_tz():
    print("b/c. STATION_TZ: local days in the station's zone (both scripts); an unknown zone is refused")
    # 2026-10-01 04:30 UTC: 00:30 on 1 Oct in New York, 21:30 on 30 Sep in Los Angeles
    e = int(datetime(2026, 10, 1, 4, 30, tzinfo=timezone.utc).timestamp())
    code = ("import extract_elevation_contours as x, dem_from_contours as d\n"
            f"print(x.LOCAL_TZ, x.local_date({e}), d.LOCAL_TZ, d.local_day({e}))")
    for env, want in ((None, "America/New_York 2026-10-01 America/New_York 2026-10-01"),
                      ({"STATION_TZ": "America/Los_Angeles"},
                       "America/Los_Angeles 2026-09-30 America/Los_Angeles 2026-09-30"),
                      ({"STATION_TZ": "UTC"}, "UTC 2026-10-01 UTC 2026-10-01")):
        rc, out, err = py(code, env)
        check(rc == 0 and out.strip() == want, f"STATION_TZ={(env or {}).get('STATION_TZ', 'unset')}: {want}",
              (rc, out, err))
    for value in ("America/Nowhere_Town", "EST5EDT; rm", "../etc"):
        refused("import extract_elevation_contours", "STATION_TZ", value, "not a zone")
        refused("import dem_from_contours", "STATION_TZ", value, "not a zone (dem_from_contours.py)")


def test_wave_currency():
    print("b/c. STATION_WAVE_CURRENCY: the currency C was fitted in")
    code = ("import extract_elevation_contours as x, dem_from_contours as d, numpy as np\n"
            "print(x.WAVE_CURRENCY_FITTED)\n"
            "print(x.wave_currency_note('buoy:wvht_m/dpd_s'))\n"
            "rec, lines = d.wave_currency_record(np.array(['f1', 'f2']), "
            "{'f1': 'buoy:wvht_m/dpd_s', 'f2': 'buoy:wvht_m/dpd_s'}, 1.0)\n"
            "print(rec['all_fitted_currency'], any(l.startswith('WARNING') for l in lines))")
    rc, out, err = py(code)
    lines = out.splitlines()
    check(rc == 0 and len(lines) == 3 and lines[0] == "marconi:hs_best/tp_s"
          and lines[1].startswith("(False, 'buoy:wvht_m/dpd_s -- NOT the currency") and lines[2] == "False True",
          "unset: marconi:hs_best/tp_s; raw-buoy frames with a setup are NOT the currency, a WARNING (as before)",
          (rc, out, err))
    rc, out, err = py(code, {"STATION_WAVE_CURRENCY": "buoy:wvht_m/dpd_s"})
    lines = out.splitlines()
    check(rc == 0 and lines[0] == "buoy:wvht_m/dpd_s"
          and lines[1] == "(True, 'buoy:wvht_m/dpd_s -- the currency the setup coefficient was fitted in')"
          and lines[2] == "True False",
          "STATION_WAVE_CURRENCY=buoy:wvht_m/dpd_s: the raw buoy IS the station's currency (no CACO-05 wording, "
          "no WARNING in the DEM's record)", (rc, out, err))
    for value in ("buoy", "buoy:hs", "buoy:hs/tp;x", "a b:c/d"):
        refused("import extract_elevation_contours", "STATION_WAVE_CURRENCY", value, "not <kind>:<h>/<p>")


def test_gnssr_qc(d):
    print("b/c. gnssr_qc.py: GNSSR_QC_GROSS_M and WAVE_BUOY")
    code = ("import gnssr_qc as g, numpy as np\n"
            "ep = 1790000000 + 600 * np.arange(7)\n"
            "lv = np.array([0.0, 0.1, 4.0, 0.2, -2.8, 0.1, 0.0])\n"
            "f, why, _, _ = g.run_qc(ep, lv)\n"
            "print(g.GROSS_MIN, g.GROSS_MAX, sorted(i for i, w in enumerate(why) if 'gross_range' in w))")
    rc, out, err = py(code)
    check(rc == 0 and out.strip() == "-2.5 3.5 [2, 4]", "unset: -2.5..3.5, so 4.0 and -2.8 m fail gross range",
          (rc, out, err))
    rc, out, err = py(code, {"GNSSR_QC_GROSS_M": "-3.0,4.5"})
    check(rc == 0 and out.strip() == "-3.0 4.5 []", "GNSSR_QC_GROSS_M=-3.0,4.5: neither fails gross range",
          (rc, out, err))
    for value in ("4,-3", "1", "a,b", "1,2,3", "nan,1"):
        refused("import gnssr_qc", "GNSSR_QC_GROSS_M", value, "not lo,hi with lo < hi")
    (d / "gauge.csv").write_text("epoch,level_navd88\n")
    for name in ("waves_44008.csv", "waves_46029.csv"):
        (d / name).write_text("epoch,wvht_m,dpd_s,apd_s\n")
    code = f"import gnssr_qc as g; p = g.default_waves_path({str(d / 'gauge.csv')!r}); print(p.name if p else None)"
    rc0, out0, _ = py(code)
    rc1, out1, err1 = py(code, {"WAVE_BUOY": "46029"})
    check(rc0 == 0 and out0.strip() == "waves_44008.csv" and rc1 == 0 and out1.strip() == "waves_46029.csv",
          "the buoy archive beside the gauge: waves_44008.csv unset, waves_<WAVE_BUOY>.csv with WAVE_BUOY=46029",
          (out0, out1, err1))
    rc2, out2, _ = py(code, {"WAVE_BUOY": "44097"})
    check(rc2 == 0 and out2.strip() == "None", "WAVE_BUOY=44097 and no waves_44097.csv: no wave file "
                                              "(never CACO-05's waves_44008.csv)", out2)
    for value in ("../x", "44 08", "a/b"):
        refused("import gnssr_qc; gnssr_qc.default_waves_path('/x/gauge.csv')", "WAVE_BUOY", value, "not an id")


def test_frames_and_names():
    print("b/c. GEOTIFF_EPSG frame labels; STATION_NAME titles")
    code = ("import station_env as s, dem_figure as f\n"
            "print(s.frame_label(), '|', s.frame_label(short=True), '|', f.STATION)")
    for env, want in ((None, "NAD83(2011) / UTM 19N, EPSG:6348 | UTM 19N | CACO05 Marconi Beach"),
                      ({"GEOTIFF_EPSG": "6339", "STATION_NAME": "CACO07 Test Beach"},
                       "NAD83(2011) / UTM 10N, EPSG:6339 | UTM 10N | CACO07 Test Beach"),
                      ({"GEOTIFF_EPSG": "26918"}, "NAD83 / UTM 18N, EPSG:26918 | UTM 18N | CACO05 Marconi Beach"),
                      ({"GEOTIFF_EPSG": "32619"}, "WGS 84 / UTM 19N, EPSG:32619 | UTM 19N | CACO05 Marconi Beach"),
                      ({"GEOTIFF_EPSG": "32756"}, "WGS 84 / UTM 56S, EPSG:32756 | UTM 56S | CACO05 Marconi Beach"),
                      ({"GEOTIFF_EPSG": "2249"}, "EPSG:2249 | EPSG:2249 | CACO05 Marconi Beach")):
        rc, out, err = py(code, env)
        check(rc == 0 and out.strip() == want, f"{env or 'unset'}: {want}", (rc, out, err))
    import station_env
    check(station_env.frame_name(6330) == "NAD83(2011) / UTM 1N" and station_env.frame_name(6328)
          == "NAD83(2011) / UTM 59N" and station_env.frame_name(6329) == "NAD83(2011) / UTM 60N"
          and station_env.frame_name(6349) is None, "NAD83(2011) UTM codes 6328-6348 (59N, 60N, 1N-19N) only")
    for value in ("EPSG:6348", "63", "6348.0"):
        refused("import station_env; station_env.frame_label()", "GEOTIFF_EPSG", value, "not a code")


# ---------------------------------------------------------------------
# d. WARNING vs NOTE
# ---------------------------------------------------------------------

def test_warning_note(d):
    print("d. extract_elevation_contours.py: the raw-buoy line is a WARNING with --setup-coef, a NOTE without")
    if not (CAL_IO_C1.is_file() and CAL_EO_C1.is_file()):
        skip(f"WARNING vs NOTE: {CAL_IO_C1.name} / {CAL_EO_C1.name} not in calibration/ (the fixture uses them)")
        return
    fx = make_fixture(d / "wn", days=1)
    hours = list(range(fx["t0"] - 3600, fx["t_end"] + 4 * 3600, 3600))
    buoy = d / "wn" / "waves_44008.csv"
    buoy.write_text("time_utc,epoch,wvht_m,dpd_s,apd_s\n"
                    + "".join(f"x,{e},{0.8 + 0.01 * k:.2f},8.0,6.0\n" for k, e in enumerate(hours)))

    def extract(name, extra, env=None):
        rc, out, err = run([HERE / "extract_elevation_contours.py", fx["spline"], "--processed-dir", fx["det"],
                            "--output", d / "wn" / f"{name}.csv", "--waves", buoy] + extra, d / "wn", env)
        lines = [ln for ln in out.splitlines() if "tagged with the raw buoy waves" in ln]
        return rc, lines, out + err

    rc, lines, log = extract("with_c", ["--setup-coef", "0.037"])
    check(rc == 0 and len(lines) == 1 and lines[0].startswith("WARNING: the frames are tagged with the raw buoy")
          and "C = 0.037 is NOT the station's calibrated setup" in lines[0],
          "--setup-coef 0.037: WARNING (the setup is not the calibrated one), as before", log[-800:])
    rc, lines, log = extract("no_c", [])
    check(rc == 0 and len(lines) == 1 and lines[0].startswith("NOTE: the frames are tagged with the raw buoy")
          and "no setup is applied here" in lines[0] and "WARNING" not in "".join(lines),
          "no --setup-coef: a NOTE (no setup applied; only a --max-hs filter judges the buoy Hs)", log[-800:])
    rc, lines, log = extract("c0", ["--setup-coef", "0"])
    check(rc == 0 and len(lines) == 1 and lines[0].startswith("NOTE:"), "--setup-coef 0 (no setup): a NOTE",
          log[-800:])
    rc, lines, log = extract("station_c", ["--setup-coef", "0.02"], {"STATION_WAVE_CURRENCY": "buoy:wvht_m/dpd_s"})
    check(rc == 0 and not lines and "Wave currency            : buoy:wvht_m/dpd_s -- the currency the setup "
                                    "coefficient was fitted in\n" in log,
          "STATION_WAVE_CURRENCY=buoy:wvht_m/dpd_s and --setup-coef: the station's own currency, no WARNING",
          log[-800:])


# ---------------------------------------------------------------------
# e. the one-camera chain on a synthetic fixture
# ---------------------------------------------------------------------

BEACH_R0, BEACH_Z0, BEACH_SLOPE = 38.0, 0.3, 0.1    # a planar beach: z = Z0 - slope * (range - R0)


def tide(t, t0):
    return 0.3 + 1.1 * np.sin(2 * np.pi * (np.asarray(t, float) - t0) / (12.42 * 3600))


def stamp(epoch):
    t = datetime.fromtimestamp(epoch, tz=timezone.utc)
    return f"{epoch}.{t:%a.%b.%d_%H_%M_%S}.GMT.{t:%Y}.CACO05"


def make_fixture(d, days=2):
    """c1 timex detections every 30 min, 11:00-20:30 UTC on `days` days, whose rows are where a planar beach
    meets the water level through the CACO05 c1 2025-11-13 calibration; the spline they were made from; a
    gauge record (gauge.csv) with the two buoy archives beside it (waves_44008.csv, waves_46029.csv) and
    the Marconi waves (waves_marconi.csv)."""
    import georectify as G
    d.mkdir(parents=True, exist_ok=True)
    det = d / "det"
    det.mkdir(exist_ok=True)
    io, eo = G.load_intrinsics(CAL_IO_C1), G.load_extrinsics(CAL_EO_C1)
    t0 = int(datetime(2026, 10, 1, 11, 0, tzinfo=timezone.utc).timestamp())
    frames = [t0 + 86400 * day + 1800 * k for day in range(days) for k in range(20)]
    t_end = frames[-1]
    sp = np.arange(t0 - 2 * 3600, t_end + 2 * 3600 + 1, 1800)
    write_spline(d / "spline.txt", sp, np.round(tide(sp, t0), 3))
    g = np.arange(sp[0], sp[-1] + 1, 360)
    (d / "gauge.csv").write_text("epoch,level_navd88\n" + "".join(
        f"{int(e)},{z:.3f}\n" for e, z in zip(g, tide(g - 2400, t0) * 0.8 - 0.05)))
    hours = np.arange(sp[0] - 3600, sp[-1] + 4 * 3600, 3600)
    for name, hs in (("waves_44008.csv", 0.9), ("waves_46029.csv", 2.5)):
        (d / name).write_text("time_utc,epoch,wvht_m,dpd_s,apd_s\n" + "".join(
            f"x,{int(e)},{hs + 0.01 * (k % 7):.2f},9.0,6.5\n" for k, e in enumerate(hours)))
    (d / "waves_marconi.csv").write_text("time_utc,epoch,hs_best,hs_sigma,tp_s\n" + "".join(
        f"x,{int(e)},{0.5 + 0.01 * (k % 5):.2f},0.1,8.0\n" for k, e in enumerate(hours)))
    cols = np.arange(500, 1900, 4).astype(float)
    vs = np.arange(1050.0, 2040.0, 2.0)
    U, V = np.meshgrid(cols, vs)
    rng = np.random.RandomState(7)
    for i, ep in enumerate(frames):
        z = float(np.interp(ep, sp, np.round(tide(sp, t0), 3)))
        X, Y = G.pixel_to_ground(U.ravel(), V.ravel(), z, io, eo)
        r = np.hypot(X - eo[0], Y - eo[1]).reshape(U.shape)          # falls as the row grows
        want = BEACH_R0 + (BEACH_Z0 - z) / BEACH_SLOPE
        lines = ["Column,Row,Confidence,Row_Precise,Peak_Sharpness,Has_Signal"]
        for j, c in enumerate(cols):
            row = float(np.interp(-want, -r[:, j], vs)) + 1.5 * np.sin(c / 9.0) + rng.normal(0, 0.3)
            lines.append(f"{int(c)},{int(round(row))},0.95,{row:.3f},0.5,1")
        (det / f"{stamp(ep)}.c1.timex.csv").write_text("\n".join(lines) + "\n")
    return {"det": det, "spline": d / "spline.txt", "t0": t0, "t_end": t_end, "frames": frames,
            "gauge": d / "gauge.csv", "waves": d / "waves_marconi.csv"}


def chain(fx, out, env=None, scripts=None):
    """GNSS-R QC report; extract (Marconi waves, C = 0.037, the QC on) -> georectify (c1 only) -> DEM (+ page)
    -> map -> GeoTIFF -> change, in folder `out`.
    scripts: a folder holding other copies of the edited scripts (--before), else this checkout's.
    -> {step: (rc, log)}"""
    out.mkdir(parents=True, exist_ok=True)
    src = Path(scripts) if scripts else HERE
    pp = os.pathsep.join([str(src), str(HERE)])
    logs = {}

    def step(name, args):
        rc, o, e = run(args, out, env, pythonpath=pp)
        logs[name] = (rc, (o + e).replace(str(out), "<out>").replace(str(src), "<src>"))
        return rc

    step("gnssr_qc", [src / "gnssr_qc.py", fx["spline"], "--reference", fx["gauge"], "--output", out / "qc.csv",
                      "--plot", out / "qc.png"])
    step("extract", [src / "extract_elevation_contours.py", fx["spline"], "--processed-dir", fx["det"],
                     "--output", out / "contours.csv", "--waves", fx["waves"], "--setup-coef", "0.037",
                     "--gnssr-qc-reference", fx["gauge"]])
    step("georectify", [HERE / "georectify.py", out / "contours.csv", out / "ground.csv",
                        "--io-c1", CAL_IO_C1, "--eo-c1", CAL_EO_C1])
    step("dem", [src / "dem_from_contours.py", out / "ground.csv", out / "dem", "--min-points", "2"])
    d1 = datetime.fromtimestamp(fx["t0"], tz=timezone.utc).strftime("%Y-%m-%d")
    d2 = datetime.fromtimestamp(fx["t_end"], tz=timezone.utc).strftime("%Y-%m-%d")
    for name, day in (("day1", d1), ("day2", d2)):
        step(f"dem_{name}", [src / "dem_from_contours.py", out / "ground.csv", out / name, "--min-points", "2",
                             "--start-date", day, "--end-date", day, "--no-plot"])
    step("map", [src / "dem_map.py", out / "dem"])
    tif = out / "tif"
    tif.mkdir(exist_ok=True)
    if (out / "dem_dem.asc").exists():
        shutil.copy(str(out / "dem_dem.asc"), str(tif / "dem_dem.asc"))
    step("geotiff", [src / "asc_to_geotiff.py", tif])
    # dem_change.compare with its figure; the labels are read off the figure before it is saved
    step("change", ["-c", "import sys, json, matplotlib; matplotlib.use('Agg')\n"
                          "import matplotlib.figure as F\n"
                          "seen = []\n"
                          "orig = F.Figure.savefig\n"
                          "def keep(self, *a, **k):\n"
                          "    seen.append([ax.get_xlabel() for ax in self.axes if ax.get_xlabel()])\n"
                          "    return orig(self, *a, **k)\n"
                          "F.Figure.savefig = keep\n"
                          "import dem_change\n"
                          f"dem_change.compare({str(out / 'day1')!r}, {str(out / 'day2')!r}, 0.05, "
                          f"out_dir={str(out / 'change')!r}, min_count=2, min_uniform_cells=5)\n"
                          "print('LABELS ' + json.dumps(seen))"])
    return logs


def page_texts(stem, info_path, env=None):
    """The heading texts and title of the DEM page and the map, as drawn (read off the figures)."""
    code = ("import json, matplotlib; matplotlib.use('Agg')\n"
            "import matplotlib.figure as F\n"
            "seen = []\n"
            "orig = F.Figure.savefig\n"
            "def keep(self, *a, **k):\n"
            "    seen.append({'texts': [t.get_text() for t in self.texts], "
            "'titles': [ax.get_title() for ax in self.axes], 'xlabels': [ax.get_xlabel() for ax in self.axes]})\n"
            "    return orig(self, *a, **k)\n"
            "F.Figure.savefig = keep\n"
            "import dem_figure, dem_map\n"
            f"info = json.load(open({str(info_path)!r}))\n"
            f"dem_figure.draw_page({str(stem)!r}, output={str(stem) + '_redrawn.png'!r}, info=info)\n"
            f"dem_figure.draw_message_page({str(stem)!r}, {str(stem) + '_short.png'!r}, info, dem_figure.STATION, "
            "'Nothing to map.', ['detail'])\n"
            f"dem_map.draw({str(stem)!r}, output={str(stem) + '_map_redrawn.png'!r})\n"
            "print('SEEN ' + json.dumps(seen))")
    rc, out, err = py(code, env)
    line = next((ln for ln in out.splitlines() if ln.startswith("SEEN ")), None)
    return (json.loads(line[5:]) if rc == 0 and line else None), out + err


def read_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def geotiff_tags(path):
    """(ImageDescription, ProjectedCSType) of a baseline little-endian GeoTIFF as asc_to_geotiff.py writes it."""
    import struct
    b = Path(path).read_bytes()
    off = struct.unpack("<I", b[4:8])[0]
    n = struct.unpack("<H", b[off:off + 2])[0]
    tags = {}
    for i in range(n):
        t, ty, c, v = struct.unpack("<HHI4s", b[off + 2 + 12 * i: off + 14 + 12 * i])
        size = {2: 1, 3: 2, 4: 4, 12: 8}.get(ty, 1) * c
        raw = v[:size] if size <= 4 else b[struct.unpack("<I", v)[0]: struct.unpack("<I", v)[0] + size]
        if ty == 2:
            tags[t] = raw.rstrip(b"\0").decode("latin-1")
        elif ty == 3:
            tags[t] = struct.unpack(f"<{c}H", raw)
    keys = tags.get(34735, ())
    pcs = next((keys[4 + 4 * k + 3] for k in range(keys[3] if keys else 0) if keys[4 + 4 * k] == 3072), None)
    return tags.get(270), pcs


def test_chain(d, before=None):
    print("e. one camera (c1), synthetic fixture: extract -> georectify -> DEM + page -> map, GeoTIFF, change")
    if not (CAL_IO_C1.is_file() and CAL_EO_C1.is_file()):
        skip(f"the chain: {CAL_IO_C1.name} / {CAL_EO_C1.name} not in calibration/")
        return
    fx = make_fixture(d / "fixture")
    runs = {"unset": None, "caco05": dict(CACO05_EXPLICIT), "station": dict(NEW_STATION)}
    logs = {k: chain(fx, d / k, env) for k, env in runs.items()}
    bad = {(k, s): v for k, lg in logs.items() for s, v in lg.items() if v[0] != 0}
    check(not bad, "every step ran (exit 0) in the three runs", {k: v[1][-600:] for k, v in bad.items()})
    if bad:
        return
    u = read_rows(d / "unset" / "ground.csv")
    frames = {r["source_file"] for r in u}
    check(len(frames) == 40 and all(r["camera"] == "c1" for r in u) and all(r["easting_utm19"] for r in u),
          f"40 c1 frames extracted and georectified on c1's calibration alone ({len(u)} points)",
          logs["unset"]["extract"][1][-600:])
    info_u = json.loads((d / "unset" / "dem_info.json").read_text())
    dem_u = np.loadtxt(d / "unset" / "dem_dem.asc", skiprows=6)
    cells_u = dem_u[dem_u != -9999.0]
    check(len(cells_u) >= 50 and info_u["frames"] == {"c1": 40},
          f"a DEM from the one camera: {len(cells_u)} cells, median {np.median(cells_u):+.2f} m (the beach "
          f"{BEACH_Z0 - BEACH_SLOPE * 11:+.1f} to {BEACH_Z0 + BEACH_SLOPE * 11:+.1f} m)")

    # unset == CACO-05's values set explicitly, file for file
    for sub in ("", "tif", "change"):
        a, b = folder_digest(d / "unset" / sub), folder_digest(d / "caco05" / sub)
        check(a and a == b, f"CACO-05's values set explicitly: every file in {sub or 'the run folder'} "
                            f"byte-identical to unset ({len(a)} files)",
              sorted(n for n in set(a) | set(b) if a.get(n) != b.get(n)))
    check(all(logs["unset"][s][1] == logs["caco05"][s][1] for s in logs["unset"]),
          "... and every step's log the same", [s for s in logs["unset"] if logs["unset"][s][1] != logs["caco05"][s][1]])
    check("station_provisional" not in info_u, "unset: info.json has no station_provisional (as before)")
    if before:
        test_before(d, fx, logs["unset"], before)

    # a new station's values
    s = read_rows(d / "station" / "contours.csv")
    u0 = read_rows(d / "unset" / "contours.csv")
    check(len(s) == len(u0) and all(abs(float(b["tide_elevation_navd88"]) - float(a["tide_elevation_navd88"]) - 1.0)
                                    < 1e-6 for a, b in zip(u0, s)),
          "GNSSR_ANTENNA_NAVD88_M=20.014: every contour level +1.000 m (the file's Hortho 19.014)")
    dem_s = np.loadtxt(d / "station" / "dem_dem.asc", skiprows=6)
    cells_s = dem_s[dem_s != -9999.0]
    shift = float(np.median(cells_s) - np.median(cells_u)) if len(cells_s) else np.nan
    check(len(cells_s) >= 50 and 0.7 < shift < 1.3,
          f"... and the DEM: median {shift:+.2f} m higher (the lines georectified on the higher plane)")
    info_s = json.loads((d / "station" / "dem_info.json").read_text())
    check(info_s.get("station_provisional") == NEW_STATION["STATION_PROVISIONAL"]
          and {k: v for k, v in info_s.items() if k != "station_provisional"}.keys() == info_u.keys(),
          "STATION_PROVISIONAL: info.json carries station_provisional, nothing else new")
    check(info_s["frames"] == {"c1": 40} and info_s["days"] == info_u["days"],
          "STATION_TZ=America/Los_Angeles: the same frames and days (11:00-20:30 UTC is the same date there)")
    seen, log = page_texts(d / "station" / "dem", d / "station" / "dem_info.json", NEW_STATION)
    if seen is None:
        check(False, "the page and the map redrawn with the station's values", log[-800:])
    else:
        page, short, mp = seen
        prov = [t for t in page["texts"] if t.startswith("PROVISIONAL: ")]
        check(page["texts"][0].startswith("Intertidal DEM, CACO07 Test Beach: ")
              and prov and prov[0].startswith("PROVISIONAL: search band from camera geometry;"),
              "the DEM page: 'Intertidal DEM, CACO07 Test Beach: ...' and a 'PROVISIONAL: ...' line", page["texts"][:6])
        check(short["texts"][0].startswith("Intertidal DEM, CACO07 Test Beach: ")
              and any(t.startswith("PROVISIONAL: ") for t in short["texts"]),
              "the short page too", short["texts"][:6])
        check(mp["titles"][0].startswith("CACO07 Test Beach: ")
              and mp["xlabels"][0] == "easting (m, NAD83(2011) / UTM 10N, EPSG:6339)",
              "the map: title 'CACO07 Test Beach: ...', axis 'easting (m, NAD83(2011) / UTM 10N, EPSG:6339)'",
              (mp["titles"], mp["xlabels"]))
    seen_u, log_u = page_texts(d / "unset" / "dem", d / "unset" / "dem_info.json")
    if seen_u is not None:
        page, short, mp = seen_u
        check(page["texts"][0].startswith("Intertidal DEM, CACO05 Marconi Beach: ")
              and not any(t.startswith("PROVISIONAL") for t in page["texts"] + short["texts"])
              and mp["titles"][0].startswith("Marconi Beach: ")
              and mp["xlabels"][0] == "easting (m, NAD83(2011) / UTM 19N, EPSG:6348)",
              "unset: CACO-05's titles and frame, no PROVISIONAL line", (page["texts"][:3], mp["titles"]))
    else:
        check(False, "the page and the map redrawn, unset", log_u[-800:])
    desc_u, pcs_u = geotiff_tags(d / "unset" / "tif" / "dem_dem.tif")
    desc_s, pcs_s = geotiff_tags(d / "station" / "tif" / "dem_dem.tif")
    check(pcs_u == 6348 and "EPSG:6348+5703" in (desc_u or "") and pcs_s == 6339 and "EPSG:6339+5703" in (desc_s or ""),
          "the GeoTIFF: EPSG 6348 unset, 6339 with GEOTIFF_EPSG=6339 (tag and description)", (pcs_u, pcs_s, desc_s))
    rc, o, e = run([HERE / "asc_to_geotiff.py", d / "station" / "tif" / "dem_dem.asc", "--epsg", "32619"],
                   d, NEW_STATION)
    check(rc == 0 and geotiff_tags(d / "station" / "tif" / "dem_dem.tif")[1] == 32619,
          "--epsg on the command line still wins over GEOTIFF_EPSG", o + e)
    rc, o, e = run([HERE / "asc_to_geotiff.py", d / "station" / "tif" / "dem_dem.asc"], d, {"GEOTIFF_EPSG": "UTM19"})
    check(rc == 1 and "ERROR: GEOTIFF_EPSG=" in e and "Traceback" not in e, "GEOTIFF_EPSG=UTM19: refused", (rc, e))
    lab = {k: next((ln for ln in logs[k]["change"][1].splitlines() if ln.startswith("LABELS ")), "") for k in logs}
    check('"easting (m, UTM 19N)"' in lab["unset"] and '"easting (m, UTM 10N)"' in lab["station"],
          "the change map: 'easting (m, UTM 19N)' unset, 'easting (m, UTM 10N)' with GEOTIFF_EPSG=6339", lab)


def test_before(d, fx, now_logs, ref):
    print(f"e'. --before {ref}: the chain with the scripts as of {ref}, variables unset: every file the same")
    old = d / "old_scripts"
    old.mkdir()
    for name in EDITED:
        r = subprocess.run(["git", "-C", str(HERE), "show", f"{ref}:{name}"], stdout=subprocess.PIPE,
                           stderr=subprocess.DEVNULL)
        if r.returncode != 0:
            skip(f"--before {ref}: git cannot read {name} there")
            return
        (old / name).write_bytes(r.stdout)
    # dem_figure.py draws the cameras from the calibration beside it: the four CACO-05 files it reads, by name
    (old / "calibration").mkdir()
    for cam in ("c1", "c2"):
        for name in (f"CACO05_{cam}_20240801_IO.yaml", f"CACO05_{cam}_20251113_EO-CV.yaml"):
            if (HERE / "calibration" / name).is_file():
                (old / "calibration" / name).symlink_to(HERE / "calibration" / name)
    logs = chain(fx, d / "before", None, scripts=old)
    for sub in ("", "tif", "change"):
        a, b = folder_digest(d / "unset" / sub), folder_digest(d / "before" / sub)
        check(a and a == b, f"{sub or 'the run folder'}: every file byte-identical to the scripts of {ref} "
                            f"({len(a)} files)", sorted(n for n in set(a) | set(b) if a.get(n) != b.get(n)))
    differ = [s for s in logs if logs[s][1] != now_logs[s][1]]
    check(not differ, f"every step's log the same as with the scripts of {ref}",
          {s: (logs[s][1][-400:], now_logs[s][1][-400:]) for s in differ})


# ---------------------------------------------------------------------
# dem_from_contours.station_setup_coef: the station file last
# ---------------------------------------------------------------------

def test_setup_coef(d):
    print("b. dem_from_contours.station_setup_coef: waterline_station.conf's SETUP_COEF after the cron's")
    import dem_from_contours as dc
    f = d / "sc"
    f.mkdir()
    cron = f / "waterline_timex_cron.sh"
    cron.write_text('BASE=/x\nSETUP_COEF="0.037"\n')
    check(dc.station_setup_coef(cron) == 0.037
          and dc.station_setup_coef(cron, source=True) == (0.037, "waterline_timex_cron.sh"),
          "no station file (CACO-05): the cron's 0.037")
    st = f / "waterline_station.conf"
    for text, want, label in (('CAMS="c1"\nSETUP_COEF="0.02"\n', (0.02, "waterline_station.conf"), "its 0.02"),
                              ('SETUP_COEF=""\n', (0.0, "waterline_station.conf"), '"" (no setup): 0'),
                              ('CAMS="c1"\n', (0.037, "waterline_timex_cron.sh"), "none set there: the cron's"),
                              ('SETUP_COEF="abc"\n', (None, "waterline_station.conf"), "a bad value: None, "
                                                                                      "never the cron's")):
        st.write_text(text)
        check(dc.station_setup_coef(cron, source=True) == want, f"station file {text.strip()!r}: {label}",
              dc.station_setup_coef(cron, source=True))
    st.unlink()
    check(dc.station_setup_coef(f / "missing.sh") is None, "no cron and no station file: None (as before)")


def main():
    ap = argparse.ArgumentParser(description="station_env.py and the scripts that read it")
    ap.add_argument("--keep", default=None, help="keep the work in this folder")
    ap.add_argument("--before", default=None, metavar="REF",
                    help="also run the chain with the edited scripts as of git revision REF and require "
                         "byte-identical outputs (variables unset)")
    a = ap.parse_args()
    print(f"test_station_env.py (Python {sys.version.split()[0]})")
    work = Path(a.keep) if a.keep else Path(tempfile.mkdtemp(prefix="station_env_"))
    if a.keep:
        if work.exists():
            shutil.rmtree(str(work))
        work.mkdir(parents=True)
    try:
        test_defaults()
        test_antenna(work)
        test_tz()
        test_wave_currency()
        test_gnssr_qc(work)
        test_frames_and_names()
        test_setup_coef(work)
        test_warning_note(work)
        test_chain(work, a.before)
    finally:
        if not a.keep:
            shutil.rmtree(str(work), ignore_errors=True)
    print()
    if SKIPPED:
        print(f"{len(SKIPPED)} skipped")
    if FAILED:
        print(f"{len(FAILED)} FAILED:")
        for m in FAILED:
            print("  - " + m)
        sys.exit(1)
    print("all station_env checks passed")


if __name__ == "__main__":
    main()

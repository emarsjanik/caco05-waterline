#!/usr/bin/env python3
"""
Tests For survey_products.py And The detect_original_view.py Photo Options
=============================================================================
Plain python3 (no pytest needed on the station):

    python3 test_survey_products.py [--keep DIR]

WHY. The survey products are judged by their labels and their numbers, so
the parts that decide them are tested by construction:
  * the live era adds the wave setup to the station's own rows itself; it
    must give EXACTLY what extract_elevation_contours.py --setup-coef writes
    (same values, same rounding, byte for byte);
  * the photo search must use each photo once when the same file sits in
    several folders (and the full-size copy, not a reduced one), keep only
    the station ID asked for, leave out the days asked for, and not hang on
    a symbolic-link loop;
  * a label can only get worse: an envelope placed with the compared
    survey, a pointing fitted to it, GCPs on the calibration day, or a setup
    coefficient fitted to it must each downgrade it;
  * the configuration must refuse the disabled dates (the user asked to skip
    23 Oct 2024 for now and 19 Mar 2025) and never place a search envelope
    with the survey of the same date;
  * a step is skipped only when its outputs exist, are newer than its inputs
    and were made with the same settings.

HOW. Everything is written into a temporary folder (--keep DIR to look) and
asserted on (--keep DIR: kept, in a new folder inside DIR). Runs in ~5 s.

STATION FILES. The checks run inside their own temporary folder (also the
current folder for the run, so a file name never finds a file of the folder
the test was started from), on fixture tables, survey files and folders of
their own: the outcome is the same on the station NUC as in a clean
checkout. Four kinds of check read the station's own state, and say so:
  * the station's surveys/survey_dates.csv and surveys/surveys.csv as they
    are now, which the owner edits (test_station_tables, run last with the
    live cron so a table problem cannot hide the code's results, every
    problem listed; a table the product cannot read at all, e.g. not UTF-8,
    is one of them): only what
    survey_products needs of them or does with them -- read with no warning
    (no text cut by a comma outside quotes), every date switched off refused
    and planned with --force-disabled, every date switched on planned, the
    calibration of every camera switched on there, and the product's own
    label rules applied to their rows (an envelope placed with the compared
    survey, a points survey of its calibration's day unless a
    label_override_reason is written: no better than PARTLY-CIRCULAR), and
    every label_override_reason they hold is listed. An edit the product
    reads with no warning and accepts never fails it (a text with a comma
    goes inside double quotes: surveys/README.md);
  * the live cron, waterline_timex_cron.sh as the station runs it (the
    owner may edit it, e.g. SETUP_COEF="" for a refit): its C, read as
    survey_products reads it, must be the products' (test_live_cron, last,
    with the tables: each runs even when the other fails);
  * the two checks on the real 2026-09-29 RTK file, run when it is found
    where surveys.csv says (the product's own rule: a path with a folder as
    given, a bare name in the survey folders);
  * the commit a stamp records, checked when the code is a git checkout.
What was asked for in the tables (the dates switched off, windows, the 29 Sep
RTK row and its text, envelopes, calibrations in the repository, the Oct GCP
label) is checked on the tables as the repository shipped them with this
test (git, as of the last commit of this file), so the owner's edits, even
committed on the station, never reach it; when the tables at HEAD differ
from those, it says so (a commit that changes the shipped tables changes
this test with them). A check that cannot run is skipped, counted and said
why.
"""

import os
import io
import re
import sys
import csv
import json
import shutil
import argparse
import time
import tempfile
import calendar
import contextlib
import subprocess
from pathlib import Path
from types import SimpleNamespace
from datetime import datetime, timezone, timedelta

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import survey_products as sp                                  # noqa: E402
from detect_original_view import collect_photos, station_of   # noqa: E402


def quiet(fn, *a, **k):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = fn(*a, **k)
    return out, buf.getvalue()


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)
    print(f"  ok  {msg}")


# where the real 2026-09-29 RTK file is looked for: the station's survey folders, as survey_products
# looks (main() adds --survey-dirs); checks that need it are skipped, and counted, when it is not found
RTK_DIRS = list(sp.DEFAULT_SURVEY_DIRS)
SKIPPED = []          # (what, why)
# main(): the test's own folders, used wherever a check needs an output root or survey folders
# ('out': an output root where nothing is built; 'surveys': an empty survey folder)
TEST_DIRS = {}


def real_rtk():
    """The real RTK file the station's C was fitted to, looked for as survey_products looks (the path
    surveys.csv gives, by the product's own rule, survey_products.resolve_file: a path with a folder in it
    as given, a bare name in the survey folders; then that name in every survey folder) and confirmed by
    its sha256. -> (path, "") or (None, why it cannot be used)."""
    fit = sp.SETUP_FITS[sp.SETUP_COEF]
    try:
        rows, _ = quiet(sp.read_table, sp.DEFAULT_SURVEYS, sp.SURVEY_FIELDS)
    except (SystemExit, Exception):     # a table survey_products cannot read (refused, not UTF-8, bad csv):
        rows = []                       # test_station_tables, run last, says so
    given = [r["path"] for r in rows if Path(r.get("path") or "").stem == fit["survey_stem"]]
    given = given or [fit["survey_stem"] + ".csv"]
    found = []
    for g in given:
        for f in sp.resolve_file(g, RTK_DIRS) + [Path(x) / Path(g).name for x in RTK_DIRS]:
            if f not in found and f.is_file():
                found.append(f)
    differ, unreadable = [], []
    for f in found:
        try:
            if sp.sha256(f) == fit["survey_sha256"]:
                return f, ""
            differ.append(str(f))
        except OSError as exc:
            unreadable.append(f"{f} ({exc.strerror or exc})")
    if unreadable:
        return None, f"the real RTK file is here but cannot be read: {', '.join(unreadable)}"
    if differ:
        return None, (f"found {', '.join(differ)}, but its sha256 is not the one C = {sp.SETUP_COEF} was fitted to "
                      f"(SETUP_FITS): re-saved or edited? The file as surveyed is needed")
    g = given[0]
    if "/" in g or os.sep in g:
        return None, (f"the real RTK file {g} (surveys.csv) is not there, nor by its name in {', '.join(RTK_DIRS)}; "
                      f"--survey-dirs adds folders")
    return None, f"the real RTK file {g} is not in {', '.join(RTK_DIRS)}; --survey-dirs adds folders"


def skip(what, why):
    SKIPPED.append((what, why))
    print(f"  skip  {what} ({why})")


def command_line_shown(proc, timeout=30.0):
    """Wait until /proc shows the command line of a process just started: when Popen returns, the
    kernel may not have set it yet, and /proc/<pid>/cmdline then reads empty, so
    survey_products.build_running() would take a build just started for another program and the lock
    checks would fail by timing alone (seen here with Python 3.11, whose subprocess starts children with
    vfork, in 33 of 240 tries with six runs at once; not seen with 3.8). A build that holds a lock has long
    had its own. -> True when shown, or when there is no /proc (alive is then all that can be told)."""
    f = Path(f"/proc/{proc.pid}/cmdline")
    if not Path("/proc/self/cmdline").exists():
        return True
    t_end = time.time() + timeout
    while time.time() < t_end:
        try:
            if f.read_bytes():
                return True
        except OSError:
            return False
        time.sleep(0.01)
    return False


def photo_name(t, station, cam):
    e = calendar.timegm(t.timetuple())
    return f"{e}.{t:%a.%b.%d_%H_%M_%S}.GMT.{t.year}.{station}.{cam}.timex.jpg"


# ---------------------------------------------------------------------

def test_setup_identical(d):
    print("live-era setup = extract_elevation_contours.py --setup-coef")
    d = Path(d) / "setup"
    (d / "proc").mkdir(parents=True)
    rng = np.random.default_rng(1)
    t0 = datetime(2026, 9, 29, 12, 0, tzinfo=timezone.utc)
    with open(d / "wl.csv", "w") as f:
        f.write("time,water_level_navd88\n")
        for k in range(48 * 6):
            t = t0 + timedelta(minutes=10 * k)
            f.write(f"{t:%Y-%m-%dT%H:%M:%SZ},{0.9 * np.sin(k / 37.0) + rng.normal(0, 0.01):.4f}\n")
    with open(d / "waves.csv", "w") as f:
        f.write("time_utc,epoch,wvht_m,dpd_s\n")
        for k in range(48):
            t = t0 + timedelta(hours=k)
            if k == 3:                      # an hour without waves: that frame keeps no setup
                continue
            f.write(f"{t:%Y-%m-%dT%H:%M:%SZ},{calendar.timegm(t.timetuple())},"
                    f"{rng.uniform(0.5, 2.5):.3f},{rng.uniform(4, 13):.2f}\n")
    for k in range(12):
        t = t0 + timedelta(minutes=30 * k + 7)
        name = photo_name(t, "CACO05", "c1").replace(".jpg", ".csv")
        cols = np.arange(0, 2448)
        rows = 900 + 40 * np.sin(cols / 90.0) + 25 * np.sin(cols / 13.0) + k * 3 + rng.normal(0, 2, cols.size)
        with open(d / "proc" / name, "w") as f:
            f.write("Column,Row,Row_Precise,Peak_Sharpness,Has_Signal\n")
            for c, r in zip(cols, rows):
                f.write(f"{c},{int(r)},{r:.3f},0.5,1\n")
    base = [sys.executable, str(HERE / "extract_elevation_contours.py"), str(d / "wl.csv"), "--time-col", "time",
            "--level-col", "water_level_navd88", "--processed-dir", str(d / "proc"), "--waves",
            str(d / "waves.csv"), "--wave-max-gap-minutes", "40", "--min-tide-agreement", "0"]
    for out, extra in (("with.csv", ["--setup-coef", "0.037"]), ("without.csv", [])):
        r = subprocess.run(base + extra + ["--output", str(d / out)], stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, universal_newlines=True)
        check(r.returncode == 0, f"extract_elevation_contours.py {' '.join(extra) or '(no setup)'} ran")
    st = sp.apply_setup(d / "without.csv", d / "mine.csv", "c1", 0.037)
    a = list(csv.DictReader(open(d / "with.csv")))
    b = list(csv.DictReader(open(d / "mine.csv")))
    check(len(a) == len(b) and len(a) > 1000, f"same rows ({len(a)})")
    same = all(x["setup_correction_m"] == y["setup_correction_m"] and
               x["beach_elevation_navd88"] == y["beach_elevation_navd88"] for x, y in zip(a, b))
    check(same, "setup_correction_m and beach_elevation_navd88 identical, as text, on every row")
    blank = [y for y in b if y["setup_correction_m"] == ""]
    check(st["frames_without_waves"] >= 1 and all(y["beach_elevation_navd88"] == y["tide_elevation_navd88"]
                                                  for y in blank),
          f"{st['frames_without_waves']} frame(s) without waves keep a blank setup and the still-water level")
    # rows that already carry a setup (the cron after 6 Oct 2026) are recomputed and compared
    st2 = sp.apply_setup(d / "with.csv", d / "mine2.csv", "c1", 0.040)
    check(st2["had_setup"] and abs(st2["max_diff_existing"] - (st2["setup_max"] * (1 - 0.037 / 0.040))) < 2e-4,
          f"existing setup columns recomputed with the run's C (max change {st2['max_diff_existing']:.4f} m)")


def test_collect_photos(d):
    print("photo search: several folders, duplicates, station ID, skipped days")
    d = Path(d) / "photos"
    for sub in ("A", "B/sub/deeper", "C", "B/loop"):
        (d / sub).mkdir(parents=True)
    big, small = b"x" * 5000, b"x" * 100
    names = []
    for day in (21, 22, 23, 24):
        for hh, mm in ((14, 0), (15, 30), (17, 0), (19, 0)):
            t = datetime(2025, 1, day, hh, mm, tzinfo=timezone.utc)
            st = "CACO04" if (day == 24 and hh >= 15) else "CACO03"
            n = photo_name(t, st, "c1")
            (d / "A" / n).write_bytes(big)
            (d / "B" / "sub" / "deeper" / n).write_bytes(big)
            names.append(n)
    reduced = photo_name(datetime(2025, 1, 22, 15, 30, tzinfo=timezone.utc), "CACO03", "c1")
    (d / "C" / reduced).write_bytes(small)                      # a reduced copy, in the FIRST folder
    (d / "A" / reduced.replace("timex", "bright")).write_bytes(small)
    os.symlink(d / "B", d / "B" / "loop" / "back")             # a loop: must not hang
    roots = [str(d / "C"), str(d / "A"), str(d / "B"), str(d / "missing")]
    ph, n = collect_photos(roots, "c1", "2025-01-21", "2025-01-24", (13.5, 18), station="CACO03")
    check(len({p.name for p in ph}) == len(ph), "each file name once")
    check(n["duplicates"] == 11 and n["duplicates_differing"] == 1,
          f"duplicates counted ({n['duplicates']}, {n['duplicates_differing']} of another size)")
    kept = [p for p in ph if p.name == reduced][0]
    check(kept.stat().st_size == len(big), "the full-size copy kept, not the reduced one listed first")
    check(all(station_of(p.name) == "CACO03" for p in ph) and n["other_station"] == 6,
          "only CACO03 photos (the CACO04 ones of 24 Jan left out)")
    check(len(ph) == 3 * 3 + 1, f"window and 13.5-18 UTC hours ({len(ph)} photos)")
    check(n["missing_roots"] == [str(d / "missing")], "a missing folder is reported")
    ph2, n2 = collect_photos(roots, "c1", "2025-01-21", "2025-01-24", (13.5, 18), station="CACO03",
                             skip_days=["2025-01-22"])
    check(not any("Jan.22" in p.name for p in ph2) and len(ph2) == len(ph) - 3, "--skip-days leaves the day out")
    ph3, _ = collect_photos([str(d / "A")], "c1", "2025-01-21", "2025-01-24", (0, 24))
    check(len(ph3) == 16, "no --station: every station, as before")
    # real JPEGs: a reduced copy is rejected against the lens file's frame size, a broken link never kept
    import cv2
    j = d / "jpg"
    for sub in ("a", "b", "c"):
        (j / sub).mkdir(parents=True)
    t = datetime(2025, 1, 18, 15, 0, tzinfo=timezone.utc)
    full, half = photo_name(t, "CACO03", "c1"), photo_name(t + timedelta(hours=1), "CACO03", "c1")
    cv2.imwrite(str(j / "b" / full), np.zeros((64, 80, 3), np.uint8))
    os.symlink(j / "nowhere.jpg", j / "a" / full)                       # broken, in the first root
    cv2.imwrite(str(j / "c" / half), np.zeros((32, 40, 3), np.uint8))   # only a reduced copy
    ph4, n4 = collect_photos([str(j / "a"), str(j / "b"), str(j / "c")], "c1", "2025-01-18", "2025-01-18",
                             (0, 24), frame_size=(80, 64))
    check([p.name for p in ph4] == [full] and ph4[0].exists(), "the good copy kept, not the broken link listed first")
    check(n4["unreadable"] == 1 and n4["wrong_size"] == 1 and n4["no_good_copy"] == 1,
          f"broken link and reduced copy counted ({n4['unreadable']}, {n4['wrong_size']}); the reduced-only photo left out")


def plan_for(d, era="adcp", env="2025005FA_Marconi_Mar_YSMP_Lidar_DSM_25cm.tif",
             eo="CACO03_c1_20250123_EO.yaml", eo_text=None):
    d = Path(d)
    eo_path = d / eo
    eo_path.write_text(eo_text or (HERE / "calibration" / "CACO03_c1_20250123_EO.yaml").read_text())
    return {"date": "2025-01-23", "era": era, "first": "2025-01-18", "last": "2025-01-23", "surveys": [],
            "out": d, "cams": {"c1": {"eo_file": eo, "eo_path": eo_path,
                                      "envelope_path": d / env if env else None, "envelope_survey": env}}}


def named(d, name):
    """A survey known by its name only: a path in the test's folder where no such file is (a bare name
    would be looked for in the current folder, the station's when the test runs from there)."""
    p = Path(d) / "not_here" / name
    assert not p.exists()
    return p


def test_labels(d):
    print("labels can only get worse")
    d = Path(d) / "labels"
    d.mkdir()
    args = SimpleNamespace(setup_coef=0.037)
    jan = {"name": "jan_lidar", "survey_type": "dsm", "survey_date": "2025-01-23", "label": "INDEPENDENT", "why": ""}
    janf = named(d, "2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif")
    rtkf = named(d, "2026-09-29_Marconi_Checkshots.csv")
    lab, why = sp.honest_label(jan, plan_for(d), args, janf)
    check(lab == "INDEPENDENT", "Jan lidar, envelope from the Mar lidar, GCP calibration: INDEPENDENT")
    lab, why = sp.honest_label(jan, plan_for(d, env="2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif"), args, janf)
    check(lab == "PARTLY-CIRCULAR", "envelope placed with the compared lidar: PARTLY-CIRCULAR")
    fitted = ("# fit_eo_to_survey.py: fitted to 2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif\n"
              + (HERE / "calibration" / "CACO03_c1_20250123_EO.yaml").read_text())
    lab, why = sp.honest_label(jan, plan_for(d, eo="CACO05_c1_2025-01-18_to_2025-01-23_lidar_EO.yaml",
                                             eo_text=fitted), args, janf)
    check(lab == "CIRCULAR", "pointing fitted to the compared lidar: CIRCULAR")
    rtk = {"name": "rtk_2026-09-29", "survey_type": "points", "survey_date": "2026-09-29",
           "label": "INDEPENDENT", "why": ""}
    p = plan_for(d, era="live", env=None, eo="CACO05_c1_20251113_EO-CV.yaml")
    lab, why = sp.honest_label(rtk, p, args, rtkf)
    check(lab == "CIRCULAR" and "setup" in why[0] and "file name" in why[0],
          "RTK with C = 0.037 (fitted to it): CIRCULAR (by its name)")
    # a C of unknown origin is never INDEPENDENT; a C declared fitted to no survey keeps the label
    for c in (0.035, 0.05, 0.0371):
        lab, why = sp.honest_label(rtk, p, SimpleNamespace(setup_coef=c), rtkf)
        check(lab == "PARTLY-CIRCULAR" and "unknown origin" in why[0], f"RTK with C = {c} of unknown origin: PARTLY-CIRCULAR")
    lab, why = sp.honest_label(rtk, p, SimpleNamespace(setup_coef=0.035, setup_fitted_to="none:repeat crossings"),
                               rtkf)
    check(lab == "INDEPENDENT", "RTK with a C declared fitted to no survey (--setup-fitted-to none:...): INDEPENDENT")
    # the real RTK file under another name is recognised by its content (station-independent: test_fit_recognised)
    real = Path(d) / "rtk_copy.csv"
    rtk_src, rtk_why = real_rtk()
    if rtk_src:
        shutil.copy(str(rtk_src), str(real))
        lab, why = sp.honest_label(dict(rtk, survey_date=""), p, args, real)
        check(lab == "CIRCULAR" and "same content" in why[0],
              "the real RTK shots renamed (rtk_copy.csv): still CIRCULAR, by sha256 (station file)")
    else:
        skip("renamed-RTK check", rtk_why)
    other = Path(d) / "another_rtk.csv"
    other.write_text("Name,Easting,Northing,Elevation,Description\n1,420100,4638400,1.2,Transect\n")
    v = sp.Verdict(rtk)
    sp.honest_label(rtk, p, args, other, verdict=v)
    check(v.label == "CIRCULAR" and "cannot be ruled out" in v.reasons[0],
          "another points survey of the fit's date: possibly the same shots -> CIRCULAR")
    v = sp.Verdict(dict(rtk, label_override_reason="a different survey: the 2026-09-29 afternoon RTK, not the fit's shots"))
    sp.honest_label(rtk, p, args, other, verdict=v)
    check(v.label == "INDEPENDENT" and v.overridden and "CIRCULAR" in v.overridden[0],
          "... overruled only by label_override_reason, and the overruled downgrade is kept")
    v = sp.Verdict(dict(rtk, label_override_reason="anything"))
    sp.honest_label(rtk, p, args, real if real.exists() else rtkf, verdict=v)
    check(v.label == "CIRCULAR", "the same shots (content or name) cannot be overruled")
    # a C declared fitted to a survey file: CIRCULAR against that file
    lab, why = sp.honest_label(rtk, p, SimpleNamespace(setup_coef=0.05, setup_fitted_to=str(other),
                                                       survey_dirs=[]), other)
    check(lab == "CIRCULAR", "a C declared fitted to this very file: CIRCULAR")
    gcp = {"name": "oct_gcps", "survey_type": "points", "survey_date": "2024-10-23", "label": "INDEPENDENT", "why": ""}
    p = plan_for(d, eo="CACO03_c1_20241023_EO.yaml")
    lab, why = sp.honest_label(gcp, p, args, named(d, "2024-10-23_Marconi_Extrinsic_Targets_c1_xyz.csv"))
    check(lab == "PARTLY-CIRCULAR", "GCPs of the calibration day, even if the table said INDEPENDENT")
    lab, why = sp.honest_label(gcp, p, args, named(d, "2024-10-23_control.csv"))
    check(lab == "PARTLY-CIRCULAR", "a point survey of the calibration day under any name ('_control.csv'): PARTLY-CIRCULAR")


def test_round2(d):
    print("carried pointings, earlier checks, horizon offsets, missing parts, DEM settings of the build")
    # The lidar fit a correction is carried from: a name unique to this run. Besides the EO's own folder the
    # label looks for it in the repository's calibration/ (sp.eo_dirs), where the station keeps its own fits,
    # so the 'fit not on this computer' case below must use a name no station file can have.
    fit_name = f"CACO05_c1_selftest_{Path(d).name}_lidar_EO.yaml"
    d = Path(d) / "round2"
    d.mkdir()
    args = SimpleNamespace(setup_coef=0.037)
    jan = {"name": "jan_lidar", "survey_type": "dsm", "survey_date": "2025-01-23", "label": "INDEPENDENT", "why": ""}
    mar = dict(jan, name="mar_lidar", survey_date="2025-03-06")
    janf, marf = (named(d, "2025005FA_Marconi_%s_YSMP_Lidar_DSM_25cm.tif" % m) for m in ("Jan", "Mar"))
    base = (HERE / "calibration" / "CACO03_c1_20250123_EO.yaml").read_text()
    # apply_pointing_correction.py carries a lidar fit's change to another period (*_corr_EO.yaml)
    (d / fit_name).write_text("# fit_eo_to_survey.py: fitted to 2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif\n" + base)
    corr = ("# apply_pointing_correction.py: CACO04_c1_20250219_EO.yaml + (+0.120, -0.310, +0.020) deg from "
            f"{fit_name}; 2024-10-01..2025-01-17\n" + base)
    p = plan_for(d, eo="CACO05_c1_20250219_corr_EO.yaml", eo_text=corr)
    v = sp.Verdict(dict(jan, label_override_reason="anything"))
    sp.honest_label(jan, p, args, janf, verdict=v)
    check(v.label == "CIRCULAR" and any("carries a correction fitted to" in r for r in v.reasons),
          "a *_corr_EO.yaml carried from the Jan lidar fit, against the Jan lidar: CIRCULAR, not overrulable")
    pm = plan_for(d, eo="CACO05_c1_20250219_corr_EO.yaml", eo_text=corr, env=janf.name)   # Mar: envelope from Jan
    lab, why = sp.honest_label(mar, pm, args, marf)
    check(lab == "INDEPENDENT", "... against the Mar lidar (another survey): INDEPENDENT")
    (d / fit_name).unlink()
    v = sp.Verdict(dict(mar, label_override_reason="the fit was to the Jan lidar (operator's notes)"))
    sp.honest_label(mar, pm, args, marf, verdict=v)
    check(v.label == "INDEPENDENT" and v.overridden and "cannot be resolved" in v.overridden[0]
          and not any((Path(x) / fit_name).exists() for x in [d] + sp.eo_dirs(args)),
          "the fit it was carried from is not here: CIRCULAR as a rule, overruled only by a written reason")
    lab, why = sp.honest_label(mar, pm, args, marf)
    check(lab == "CIRCULAR", "... and without that reason: CIRCULAR")
    # the earlier check of 2025-01-23 on the REAL photos leads the caveats, per camera, with what to
    # expect now. pointing_fix_figure.py printed lidar - water level = c1 +0.46, c2 +0.15: the lines
    # read LOW (waterline - lidar c1 -0.46, c2 -0.15), the sign a line without setup must have.
    plan = {"era": "adcp", "date": "2025-01-23", "first": "2025-01-18", "last": "2025-01-23", "out": d,
            "surveys": [dict(jan, path=janf.name)],
            "cams": {c: {"eo_file": "CACO03_%s_20250123_EO.yaml" % c, "first": "2025-01-18", "last": "2025-01-23"}
                     for c in ("c1", "c2")}}
    pc0 = sp.PRIOR_CHECKS[0]
    check(all(abs(pc0["waterline_minus_survey"][c][0] + pc0["reported"][c][0]) < 1e-12 for c in ("c1", "c2"))
          and pc0["waterline_minus_survey"]["c1"][0] < 0 and "lidar elevation - water level" in pc0["reported_as"],
          "the earlier check's sign: lidar - water level +0.46/+0.15 = waterline - lidar -0.46/-0.15 (lines LOW)")
    forcing = {"waves": {"setup_in_window": {"median_m": 0.29, "daytime_median_m": 0.28}}}
    cav = sp.collect_caveats(plan, SimpleNamespace(setup_coef=sp.SETUP_COEF), {}, forcing, {})
    c0 = cav[0]
    check(c0.startswith("EARLIER CHECK (real photos)") and "OPPOSITE" not in c0 and "HIGHER" not in c0
          and "+0.46 m (NMAD 0.19)" in c0 and "+0.15 m (NMAD 0.23)" in c0
          and "c1 -0.46 m (NMAD 0.19)" in c0 and "c2 -0.15 m (NMAD 0.23)" in c0
          and "read LOW" in c0 and "no opposite-sign discrepancy" in c0,
          "2025-01-23: the earlier result comes first, printed as measured (+0.46/+0.15 lidar - level) and as "
          "waterline - lidar (-0.46/-0.15: LOW, the same sign as no setup)")
    check("c1 -0.29 to -0.24, c2 +0.02 to +0.07 m" in c0 and "c1 -0.18, c2 +0.13 m if the" in c0
          and "0.71 x its setup" not in c0 and "+/-0.01" not in c0 and "0.6-0.8 x its setup" in c0,
          "... without a paired C = 0 effect of its own: c1 -0.29..-0.24 / c2 +0.02..+0.07 m over the shares kept "
          "measured so far (0.6-0.8), no carried-over 0.71 and no false +/-0.01 (-0.18 / +0.13 m if all counted)")
    # with this build's own paired C = 0 sensitivity, the expectation is the earlier lines plus that rise
    prov_ = {"comparisons": [{"current": True, "survey_path": janf.name, "headline": {"waterlines": {
        "without_setup": {"paired": {"by_camera": {
            "c1": {"n": 40, "frames": 40, "setup_effect_median": 0.176, "n_independent": 40,
                   "n_independent_unit": "frames"},
            "c2": {"n": 30, "frames": 30, "setup_effect_median": 0.160, "n_independent": 30,
                   "n_independent_unit": "frames"}}}}}}}]}
    c0b = sp.prior_check_caveat(plan, dict(pc0, cameras=["c1", "c2"]), [{"coef": sp.SETUP_COEF}], forcing, prov_)
    check("c1 ~-0.28 m (-0.46 +0.176, on 40 frames)" in c0b and "c2 ~+0.01 m (-0.15 +0.160, on 30 frames)" in c0b
          and "c1 +0.176 m (0.63 x the setup median)" in c0b,
          "... with this build's own paired effect per camera: c1 -0.46 + 0.176 ~ -0.28, c2 -0.15 + 0.160 ~ +0.01")
    # a second survey's effect on a handful of points never replaces it (nor is used alone)
    few_pts = {"current": True, "survey_path": "2025-01-23_points_xyz.csv", "headline": {"waterlines": {
        "without_setup": {"paired": {"by_camera": {
            "c1": {"n": 13, "setup_effect_median": 0.28, "n_independent": 3, "n_independent_unit": "survey points"},
            "c2": {"n": 13, "setup_effect_median": 0.28, "n_independent": 4,
                   "n_independent_unit": "survey points"}}}}}}}
    eff2 = sp.build_setup_effect({"comparisons": prov_["comparisons"] + [few_pts]},
                                 [janf.stem, "2025-01-23_points_xyz"])
    check(eff2 == {"c1": (0.176, "40 frames"), "c2": (0.16, "30 frames")}
          and sp.build_setup_effect({"comparisons": [few_pts]}, ["2025-01-23_points_xyz"]) == {},
          "... the effect on the most frames is used; one on 3-4 survey points is too few to use")
    check("PARTLY-CIRCULAR: its search envelope was placed with this same Jan 2025 lidar" in c0,
          "... the earlier check carries its own label: PARTLY-CIRCULAR (envelope placed with the same lidar)")
    under = [c for c in cav if c.startswith("C = 0.037 is not pinned by the 2026 RTK")]
    check("do not subtract them again" in c0 and under and "EARLIER CHECK expectation above is made with C = 0.037 "
          "and already holds the re-projection" in under[0] and "would sit ~0.03-0.05 m higher than with C = 0.037" in under[0]
          and "biased LOW by where that survey stops" in under[0] and "C ~0.043-0.046" in under[0]
          and "-0.04 to -0.03 m (7-10 frames: too few" in under[0] and "0.043 (0.026-0.053)" in under[0]
          and "under-corrects" not in under[0],
          "... the 2026 RTK numbers are said to be biased LOW by where that survey stops; above it the lines read "
          "-0.04..-0.03 m and C ~0.043-0.046, internal 0.043: C = 0.037 is not shown to under-correct; the "
          "re-projection is inside the expectation, never counted twice")
    check("station computer only" in c0 and "synthetic test fixture" in c0 and "nothing about the real beach" in c0,
          "... and says the real photos are on the station only; a fixture build says nothing about the real beach")
    cav0 = sp.collect_caveats(plan, SimpleNamespace(setup_coef=0.0), {}, forcing, {})
    check(cav0[0].startswith("EARLIER CHECK (real photos)") and "This build adds the setup" not in cav0[0],
          "... with C = 0: no expectation from a setup that is not applied")
    other = dict(plan, cams={"c1": dict(plan["cams"]["c1"], eo_file="CACO05_c1_20250219_corr_EO.yaml")})
    check(not sp.prior_checks(other), "... not claimed for another calibration")
    # a constant horizon offset is reported from 0.05 deg (not only above the 12 px day-to-day limit)
    info = {"window_horizon_dtilt_deg": 0.12, "window_horizon_droll_deg": 0.01, "clear_horizon_days": 5,
            "window_horizon_offset_px": 5.6}          # made-up values: a unit test, not a measurement
    off = sp.pointing_offset(info, "CACO03_c1_20250123_EO.yaml")
    check(off and abs(off[0] - 0.12) < 1e-9 and "6 px" in off[2], "a constant 0.12 deg tilt under 12 px: reported")
    check(sp.pointing_offset(dict(info, window_horizon_dtilt_deg=0.03), "x") is None, "0.03 deg: not reported")
    # a camera whose pointing was not checked or whose photo map is missing: status partial
    miss = sp.missing_parts({"out": d, "date": "2025-01-23", "cams": {"c1": {}}},
                            {"pointing": {"c1": {"rows": [], "info": {}}}})
    check(any("pointing NOT checked" in m for m in miss) and any("no waterline map" in m for m in miss),
          "no photos for the pointing check or the photo map: listed as missing parts")
    # a survey surveys.csv names whose file is not found: a missing part (status partial, exit 3);
    # a row with no file yet (the Oct GCPs) is not
    empty = d / "no_surveys_here"
    empty.mkdir()
    named_row = dict(jan, path=janf.name)
    nofile = {"name": "oct_gcps", "survey_type": "points", "survey_date": "2024-10-23", "label": "PARTLY-CIRCULAR",
              "why": "", "path": ""}
    miss = sp.missing_parts({"out": d, "date": "2025-01-23", "cams": {}, "surveys": [named_row, nofile]}, {},
                            SimpleNamespace(survey_dirs=[str(empty)]))
    check(any(m.startswith("jan_lidar: survey file") and "NOT FOUND" in m for m in miss)
          and not any(m.startswith("oct_gcps") for m in miss),
          "a named survey file not found: a missing part; a row with no file yet: not")
    # the DEM settings reported are those the DEM on disk was built with
    sp.write_stamp(d, "dem", {"cmd": "python3 dem_from_contours.py x.csv out --cell 2.0 --min-points 3 --max-spread 0.5 "
                                     "--max-hs 1.5 --max-day-offset 0.15", "tif": "python3 asc_to_geotiff.py a --epsg 32619"})
    built = sp.dem_settings_built(d)
    check(built == {"cell_m": 2.0, "min_points": 3, "max_spread_m": 0.5, "max_hs_m": 1.5, "max_day_offset_m": 0.15,
                    "geotiff_epsg": 32619}, "DEM settings read back from the build's own commands")
    want = sp.dem_settings_requested(SimpleNamespace(dem_cell=1.0, dem_min_points=3, dem_max_spread=0.5, dem_max_hs=3.0,
                                                     dem_max_day_offset=0.15, geotiff_epsg=32619))
    dd = sp.dem_settings_differ(built, want)
    check("cell_m: asked 1.0, DEM on disk 2.0" in dd and "max_hs_m: asked 3.0, DEM on disk 1.5" in dd
          and "min_points" not in dd,
          "a later command line's other DEM options are named as NOT applied, asked vs on disk the right way round")
    later = SimpleNamespace(dem_cell=1.0, dem_min_points=3, dem_max_spread=0.5, dem_max_hs=3.0, dem_max_day_offset=0.15,
                            geotiff_epsg=6348, setup_coef=sp.SETUP_COEF)
    b = sp.built_args(later, {"out": d, "date": "2025-01-23", "cams": {"c1": {}}})
    check(b.dem_cell == 2.0 and b.dem_max_hs == 1.5 and b.geotiff_epsg == 32619 and later.dem_cell == 1.0,
          "the 'rebuild everything' line takes the DEM settings the DEM was built with, not this run's")
    # a re-export of the real RTK the station's C was fitted to, entered under another date: same shots
    # (station-independent: test_fit_recognised)
    rtk_src, rtk_why = real_rtk()
    if rtk_src:
        lines = rtk_src.read_text().splitlines()
        lines[1] = lines[1].replace(",", ", ", 1)
        rex = d / "rtk_reexport.csv"
        rex.write_text("\n".join(lines) + "\n")
        rtk = {"name": "rtk_other", "survey_type": "points", "survey_date": "2026-09-30", "label": "INDEPENDENT",
               "why": ""}
        pl = plan_for(d, era="live", env=None, eo="CACO05_c1_20251113_EO-CV.yaml")
        lab, why = sp.honest_label(rtk, pl, SimpleNamespace(setup_coef=0.037, survey_dirs=[str(rtk_src.parent)]), rex)
        check(lab == "CIRCULAR" and "re-exported" in why[0], "a byte-changed re-export of the real RTK under another "
                                                               "date: CIRCULAR by its coordinates (station file)")
    else:
        skip("re-exported-RTK check", rtk_why)


def test_fit_recognised(d):
    print("the survey a setup coefficient was fitted to: by content, by its coordinates, by name (a test fit)")
    d = Path(d) / "fitrec"
    station, elsewhere = d / "station", d / "elsewhere"
    for x in (station, elsewhere):
        x.mkdir(parents=True)
    rng = np.random.RandomState(29)
    E, N = 420050.0 + rng.uniform(0, 80, 40), 4638300.0 + rng.uniform(0, 150, 40)

    def shots(path, dE=0.0, dated=True):
        """Emlid-like RTK shots; dE moves every point east; dated: with its 'Averaging start' column."""
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Name", "Code", "Code description", "Easting", "Northing", "Elevation", "Description",
                        "CS name"] + (["Averaging start"] if dated else []))
            for k in range(len(E)):
                w.writerow([str(k + 1), "T", "Transect", f"{E[k] + dE:.3f}", f"{N[k]:.3f}",
                            f"{0.5 + 0.03 * (420130.0 - E[k]):.3f}", "Transect",
                            "NAD83(2011) / UTM zone 19N + NAVD88(GEOID18) height"]
                           + (["2026-10-15 14:00:00 UTC-04:00"] if dated else []))
        return path
    fitted = shots(station / "selftest_fit_shots.csv")
    fake = 0.0123                                 # a fit made for this test only, to these shots
    sp.SETUP_FITS[fake] = {"survey_stem": "selftest_fit_shots", "survey_sha256": sp.sha256(fitted),
                           "survey_type": "points", "survey_date": "2026-10-15", "how": "a fit for this test only"}
    try:
        p = plan_for(d, era="live", env=None, eo="CACO05_c1_20251113_EO-CV.yaml")
        here = SimpleNamespace(setup_coef=fake, survey_dirs=[str(station)])
        fits = [dict(sp.setup_fit_info(fake, here), cameras=["c1"])]
        s = {"name": "shots", "survey_type": "points", "survey_date": "", "label": "INDEPENDENT", "why": ""}

        def label(path, args=here):
            return sp.honest_label(s, p, args, path, fits=fits)
        renamed = elsewhere / "renamed_export.csv"
        shutil.copy(str(fitted), str(renamed))
        lab, why = label(renamed)
        check(lab == "CIRCULAR" and "same content (sha256)" in why[0],
              "the fitted shots under another name: CIRCULAR, by their content")
        rex = shots(elsewhere / "reexport.csv", dE=0.004, dated=False)
        lab, why = label(rex)
        check(lab == "CIRCULAR" and "40 of its 40 points lie within 0.05 m" in why[0] and "re-exported" in why[0],
              "a byte-changed re-export (+4 mm east, another name, no date in it): CIRCULAR by its coordinates")
        lab, why = label(shots(elsewhere / "another_survey.csv", dE=1.0, dated=False))
        check(lab == "INDEPENDENT" and not why, "another survey (every point 1 m away), no date: INDEPENDENT")
        lab, why = label(named(d, "selftest_fit_shots.csv"))
        check(lab == "CIRCULAR" and "file name selftest_fit_shots.csv" in why[0],
              "the fit's file name, the file itself not here: CIRCULAR by its name")
        # the fit's own shots are looked for in the survey folders only: a copy in the current folder is not used
        start = os.getcwd()
        os.chdir(str(elsewhere))
        try:
            shutil.copy(str(fitted), "selftest_fit_shots.csv")
            lab, why = label(rex, SimpleNamespace(setup_coef=fake, survey_dirs=[str(TEST_DIRS["surveys"])]))
        finally:
            os.chdir(start)
        check(lab == "INDEPENDENT", "... the same re-export, the fit's shots in no survey folder (a copy in the "
                                    "current folder): its coordinates cannot be compared, nothing is recognised")
    finally:
        del sp.SETUP_FITS[fake]


def test_cwd_never_searched(d):
    print("file names in the tables: looked for in the survey folders, never in the current folder")
    d = Path(d) / "cwd"
    run_dir, sdir, empty = d / "run_here", d / "surveys", d / "empty"
    for x in (run_dir, sdir, empty):
        x.mkdir(parents=True)
    start = os.getcwd()
    os.chdir(str(run_dir))
    try:
        Path("stray.csv").write_text("x\n")
        (sdir / "found.csv").write_text("x\n")
        check(sp.resolve_file("stray.csv", [str(empty)]) == []
              and sp.resolve_file("found.csv", [str(empty), str(sdir)]) == [sdir / "found.csv"],
              "a bare file name: found in the folders given, not in the current folder")
        check(sp.resolve_file("./stray.csv", [str(empty)]) == [Path("stray.csv")]
              and sp.resolve_file("stray.csv", [str(empty)], typed=True) == [Path("stray.csv")],
              "a path with a folder (./stray.csv), or a name typed on the command line: as given")
        # a survey whose file sits only in the current folder: NOT FOUND (a missing part), its label goes by its
        # name, and the storm caveat reads nothing of it; once in a survey folder it is read
        row = {"name": "rtk", "survey_type": "points", "survey_date": "2026-09-29", "label": "INDEPENDENT", "why": "",
               "path": "probe_shots.csv"}
        Path("probe_shots.csv").write_text("Name,Easting,Northing,Elevation,Description\n" + "".join(
            f"{k},420100.0,{4638380.0 + 30 * k:.1f},1.0,Transect\n" for k in range(5)))
        a = SimpleNamespace(survey_dirs=[str(empty)])
        miss = sp.missing_parts({"out": d, "date": "2026-09-29", "cams": {}, "surveys": [row]}, {}, a)
        found = sp.survey_path_for(row, None, a)
        pl = {"date": "2026-09-29", "era": "live", "first": "2026-09-26", "last": "2026-10-02", "out": d, "cams": {},
              "surveys": [row]}
        cav = sp.storm_caveats(pl, a, {}, "2026-09-26", "2026-10-02")
        check(any(m.startswith("rtk: survey file") and "NOT FOUND" in m for m in miss) and not found.is_file()
              and found.name == "probe_shots.csv" and len(cav) == 1 and "overlap that loss" not in cav[0],
              "a survey file only in the current folder: NOT FOUND, labelled by its name, the storm caveat reads nothing")
        shutil.copy("probe_shots.csv", str(sdir / "probe_shots.csv"))
        cav = sp.storm_caveats(pl, SimpleNamespace(survey_dirs=[str(sdir)]), {}, "2026-09-26", "2026-10-02")
        check(len(cav) == 1 and "The rtk transects span N 4638380-4638500: they overlap that loss." in cav[0],
              "... once in a survey folder: found, and the storm caveat reads it")
        # --dry-run and --summary (a date not built yet) label the survey as the build would: by its name in the
        # survey folders, never by the current folder's copy. A setup fit of the test's own whose shots have
        # this file's content: CIRCULAR only where the file is found.
        fake = 0.0123
        while any(abs(k - fake) < 1e-6 for k in sp.SETUP_FITS):
            fake += 0.0001
        sp.SETUP_FITS[fake] = {"survey_stem": "selftest_fit_shots", "survey_sha256": sp.sha256(Path("probe_shots.csv")),
                               "survey_type": "points", "survey_date": "2026-10-15", "how": "a test fit"}
        try:
            row_t = dict(row, date="2026-09-29", label_override_reason="")
            cfg_t = [{"date": "2026-09-29", "camera": "c1", "enabled": "1", "station": "CACO05",
                      "first_day": "2026-09-26", "last_day": "2026-10-02", "eo_file": "CACO05_c1_20251113_EO-CV.yaml",
                      "envelope_survey": "", "utc_hours": "0-24", "era": "live", "notes": ""}]
            for dirs, want in (([str(empty)], "INDEPENDENT"), ([str(sdir)], "CIRCULAR")):
                a = SimpleNamespace(setup_coef=fake, setup_fitted_to=None, survey_dirs=dirs, calibration=str(empty),
                                    photo_roots=[str(empty)], live_contours=str(empty / "none.csv"),
                                    gnssr_spline=str(empty / "none.txt"), gauge_csv=str(empty / "none.csv"),
                                    surveys="surveys.csv", config="survey_dates.csv", force_disabled=False,
                                    window=None, utc_hours=None, output_root=str(d / "no_output_root"))
                plan_t = quiet(sp.build_plan, "2026-09-29", cfg_t, [row_t], a)[0]
                _, out = quiet(sp.dry_run, plan_t, a, ["compare"])
                lab = [ln.split(":", 1)[1].strip() for ln in out.splitlines() if ln.strip().startswith("label")]
                when, _ = sp.label_when_built("2026-09-29", row_t, cfg_t, [row_t], a)
                where = "only in the current folder" if want == "INDEPENDENT" else "in a survey folder"
                check(len(lab) == 1 and lab[0].split(" ")[0] == want and when == want
                      and (want == "INDEPENDENT" or "same content (sha256)" in lab[0]),
                      f"--dry-run and --summary, the survey file {where}: {want} (a setup fit whose shots have its "
                      f"content: CIRCULAR only where the build finds it)")
        finally:
            del sp.SETUP_FITS[fake]
    finally:
        os.chdir(start)


def test_setup_in_use(d):
    print("the label and provenance use the C the waterlines were built with")
    d = Path(d) / "inuse"
    (d / "waterlines").mkdir(parents=True)
    p = plan_for(d, era="live", env=None, eo="CACO05_c1_20251113_EO-CV.yaml")
    rtk = {"name": "rtk_2026-09-29", "survey_type": "points", "survey_date": "2026-09-29",
           "label": "INDEPENDENT", "why": ""}
    built = SimpleNamespace(setup_coef=0.037)
    sp.write_stamp(d, "merge", {"x": 1}, extra={"setup": {"c1": {"coef": 0.037,
                                                                  "fit": sp.setup_fit_info(0.037, built),
                                                                  "implied_coef": 0.0370}}})
    later = SimpleNamespace(setup_coef=0.02, setup_fitted_to="none:test")
    fits, src = sp.setup_in_use(p, later)
    check([f["coef"] for f in fits] == [0.037] and src.startswith("the waterlines"),
          "a later --setup-coef 0.02 --steps compare: the C in use is the one the waterlines carry (0.037)")
    lab, why = sp.honest_label(rtk, p, later, named(d, "2026-09-29_Marconi_Checkshots.csv"))
    check(lab == "CIRCULAR", "... so the RTK comparison stays CIRCULAR")
    check(abs(sp.implied_c(0.3541, 0.65, 9.5) - 0.037) < 2e-4, "C implied by a row's own setup, Hs and Tp")


def test_merge_drops_no_setup(d):
    print("frames without a wave record (no setup) are left out when C > 0")
    d = Path(d) / "nosetup"
    (d / "waterlines" / "c1").mkdir(parents=True)
    f = d / "waterlines" / "c1" / "contour_points_ground.csv"
    f.write_text("source_file,camera,capture_time_utc,capture_epoch,tide_elevation_navd88,offshore_hs_m,"
                 "offshore_tp_s,setup_correction_m,beach_elevation_navd88,easting_utm19,northing_utm19\n"
                 "a,c1,2026-09-29T14:00:00Z,1790690400,0.5,1.0,9.0,0.3,0.8,420100,4638400\n"
                 "a,c1,2026-09-29T14:00:00Z,1790690400,0.5,1.0,9.0,0.3,0.8,420101,4638401\n"
                 "b,c1,2026-09-29T15:00:00Z,1790694000,0.6,,,,0.6,420100,4638400\n")
    plan = {"date": "2026-09-29", "era": "live", "out": d, "first": "2026-09-29", "last": "2026-09-29",
            "cams": {"c1": {}}, "surveys": []}
    st, out = quiet(sp.merge_cameras, plan, [("c1", f)], SimpleNamespace(setup_coef=0.037))
    rows = list(csv.DictReader(open(d / "waterlines" / "contour_points_ground.csv")))
    check(len(rows) == 2 and all(r["source_file"] == "a" for r in rows), "the frame without Hs/Tp is left out")
    check(st["no_setup_frames"] == {"c1": 1} and "WARNING" in out, "counted per camera and warned")
    listed = list(csv.DictReader(open(d / "waterlines" / "no_setup_frames.csv")))
    check(len(listed) == 1 and listed[0]["source_file"] == "b", "listed in no_setup_frames.csv")
    st, _ = quiet(sp.merge_cameras, plan, [("c1", f)], SimpleNamespace(setup_coef=0.0))
    check(st["rows"] == 3, "C = 0: nothing left out (no frame has a setup)")


def cfg_args(**kw):
    """Options as survey_products' command line gives them, with the test's own output root (nothing is
    built there) and survey folder (empty) instead of the station's."""
    a = dict(config=str(sp.DEFAULT_CONFIG), surveys=str(sp.DEFAULT_SURVEYS), force_disabled=False, window=None,
             utc_hours=None, calibration=str(HERE / "calibration"), survey_dirs=[str(TEST_DIRS["surveys"])],
             output_root=str(TEST_DIRS["out"]), setup_coef=sp.SETUP_COEF)
    a.update(kw)
    return SimpleNamespace(**a)


# Fixture tables: the behaviour of the configuration (a date switched off, a camera switched off, the
# window, the labels, GCPs that arrive later) is tested on these, not on surveys/*.csv, which the owner
# edits on the station as surveys/README.md says (switch a date on, give the GCP file's path).
FIXTURE_DATES = """date,camera,enabled,station,first_day,last_day,eo_file,envelope_survey,utc_hours,era,notes
2024-10-23,c1,0,CACO03,2024-10-20,2024-10-26,CACO03_c1_20241023_EO.yaml,2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif,13.5-18,chatham,"switched off (a test row), with a comma"
2024-10-23,c2,0,CACO03,2024-10-20,2024-10-26,CACO03_c2_20241023_EO.yaml,2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif,13.5-18,chatham,"switched off (a test row), with a comma"
2025-01-23,c1,1,CACO03,2025-01-18,2025-01-23,CACO03_c1_20250123_EO.yaml,2025005FA_Marconi_Mar_YSMP_Lidar_DSM_25cm.tif,13.5-18,adcp,
2025-01-23,c2,1,CACO03,2025-01-18,2025-01-23,CACO03_c2_20250123_EO.yaml,2025005FA_Marconi_Mar_YSMP_Lidar_DSM_25cm.tif,13.5-18,adcp,
2025-03-06,c1,0,CACO04,2025-03-03,2025-03-09,CACO04_c1_20250219_EO.yaml,2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif,13.5-18,adcp,no c1 photos (a test row)
2025-03-06,c2,1,CACO04,2025-03-03,2025-03-09,CACO04_c2_20250219_EO.yaml,2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif,13.5-18,adcp,
2026-09-29,c1,1,CACO05,2026-09-26,2026-10-02,CACO05_c1_20251113_EO-CV.yaml,,0-24,live,
2026-09-29,c2,1,CACO05,2026-09-26,2026-10-02,CACO05_c2_20251113_EO-CV.yaml,,0-24,live,
"""
FIXTURE_SURVEYS = """date,name,path,survey_type,survey_date,label,why,label_override_reason
2024-10-23,oct_gcps,,points,2024-10-23,PARTLY-CIRCULAR,"GCPs of the calibration day, file not here yet (a test row)",
2025-01-23,jan_lidar,2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif,dsm,2025-01-23,INDEPENDENT,a test row,
2026-09-29,rtk_2026-09-29,2026-09-29_Marconi_Checkshots.csv,points,2026-09-29,INDEPENDENT,a test row,
"""


def fixture_tables(d):
    """-> (cfg rows, survey rows, cfg_args keywords naming the two fixture tables written into d/tables)."""
    t = Path(d) / "tables"
    t.mkdir(parents=True, exist_ok=True)
    (t / "survey_dates.csv").write_text(FIXTURE_DATES)
    (t / "surveys.csv").write_text(FIXTURE_SURVEYS)
    return (sp.read_table(t / "survey_dates.csv", sp.CONFIG_FIELDS), sp.read_table(t / "surveys.csv", sp.SURVEY_FIELDS),
            {"config": str(t / "survey_dates.csv"), "surveys": str(t / "surveys.csv")})


def is_enabled(row):
    return row["enabled"] in ("1", "yes", "true", "True")


def calibration_day(eo_file):
    """The day a calibration was solved, from its name (CACO03_c1_20241023_EO.yaml -> 2024-10-23), read
    here by the test itself, so a check of the product's calibration-day rule does not use the product's."""
    m = re.search(r"_(\d{4})(\d\d)(\d\d)_", Path(str(eo_file)).name)
    return "-".join(m.groups()) if m else None


def test_config(d):
    print("the configuration: its behaviour on fixture tables")
    d = Path(d) / "config"
    d.mkdir()
    cfg, srv, tabs = fixture_tables(d)
    plan, out = quiet(sp.build_plan, "2024-10-23", cfg, srv, cfg_args(**tabs))
    check(plan is None and "DISABLED" in out and "switched off (a test row), with a comma" in out,
          "a date with no enabled row is refused, with the reason from its notes")
    plan, out = quiet(sp.build_plan, "2024-10-23", cfg, srv, cfg_args(force_disabled=True, **tabs))
    check(plan is not None and plan["forced_disabled"] and list(plan["cams"]) == ["c1", "c2"],
          "... and builds with --force-disabled")
    plan, _ = quiet(sp.build_plan, "2025-03-06", cfg, srv, cfg_args(**tabs))
    check(list(plan["cams"]) == ["c2"] and plan["cams"]["c2"]["station"] == "CACO04",
          "a camera switched off: the date is built with the other only, and its station ID (2025-03-06: c2, CACO04)")
    plan, _ = quiet(sp.build_plan, "2026-09-29", cfg, srv, cfg_args(**tabs))
    check(plan["first"] == "2026-09-26" and plan["last"] == "2026-10-02" and plan["era"] == "live"
          and [s["name"] for s in plan["surveys"]] == ["rtk_2026-09-29"],
          "the window, era and surveys of a date come from its rows (2026-09-29: 26 Sep - 2 Oct, live, the 29 Sep RTK)")
    plan, _ = quiet(sp.build_plan, "2025-01-23", cfg, srv, cfg_args(**tabs))
    check(plan["last"] == "2025-01-23" and all(c["station"] == "CACO03" for c in plan["cams"].values()),
          "2025-01-23: CACO03 only, 24 Jan (the re-set) outside the window")
    # the RTK the station's C was fitted to, known here by its name only (content: test_fit_recognised)
    rtk_row = [s for s in srv if s["name"] == "rtk_2026-09-29"][0]
    rtk_file = Path(TEST_DIRS["surveys"]) / rtk_row["path"]
    pl, _ = quiet(sp.build_plan, "2026-09-29", cfg, srv, cfg_args(**tabs))
    lab, why = sp.honest_label(rtk_row, pl, SimpleNamespace(setup_coef=sp.SETUP_COEF), rtk_file)
    check(lab == "CIRCULAR" and "setup coefficient C = 0.037 was fitted to this survey (file name" in why[0],
          "... a default build (C = 0.037): CIRCULAR by the setup rule (here by the file's name)")
    lab, why = sp.honest_label(rtk_row, pl, SimpleNamespace(setup_coef=0.0), rtk_file)
    check(lab == "INDEPENDENT" and not why, "... a C = 0 build (fitted to nothing): INDEPENDENT")
    # The GCP files arrive later, one per camera, under the survey date's name: picked up by
    # themselves and merged, a target both cameras saw counted once.
    oct_ = [s for s in srv if s["name"] == "oct_gcps"][0]
    sdir = Path(d) / "gcp_arrive"
    sdir.mkdir()
    a = cfg_args(survey_dirs=[str(sdir)], **tabs)
    paths, note = sp.resolve_survey(oct_, None, a)
    check(not paths and "no survey file" in note, "Oct GCPs not there yet: skipped, with the reason")
    (sdir / "2024-10-23_Marconi_Extrinsic_Targets_c1_xyz.csv").write_text(
        "1,420100.0,4638400.0,1.20\n2,420110.0,4638410.0,0.80\n")
    (sdir / "2024-10-23_Marconi_Extrinsic_Targets_c2_xyz.csv").write_text(
        "2,420110.0,4638410.0,0.80\n3,420120.0,4638500.0,0.50\n")
    paths, note = sp.resolve_survey(oct_, None, a)
    check(len(paths) == 2 and "found" in note, "Oct GCPs arrive: both camera files found by their date")
    rows = [ln for ln in Path(sp.merged_points(paths, sdir / "merged.csv")).read_text().splitlines() if ln]
    check(len(rows) == 3, f"Oct GCPs merged: 3 targets, the one both cameras saw once (got {len(rows)})")


def test_table_commas(d):
    print("a comma in a table's text: inside double quotes, or said")
    d = Path(d) / "commas"
    d.mkdir()
    bad = d / "bad_dates.csv"
    bad.write_text("date,camera,enabled,station,first_day,last_day,eo_file,envelope_survey,utc_hours,era,notes\n"
                   "2025-03-06,c2,1,CACO04,2025-03-03,2025-03-09,CACO04_c2_20250219_EO.yaml,,13.5-18,adcp,"
                   "c2 only, a comma not in quotes\n")
    rows, out = quiet(sp.read_table, bad, sp.CONFIG_FIELDS)
    check("WARNING" in out and "line 2: 1 field(s) more than its 11 columns" in out and rows[0]["notes"] == "c2 only",
          "a table row with more fields than columns (a comma in a text not in quotes): a WARNING naming the line")
    # a row that only ends in a comma (or two, or one and a space) has an empty extra field: nothing is dropped
    trail = d / "trailing_dates.csv"
    trail.write_text("date,camera,enabled,station,first_day,last_day,eo_file,envelope_survey,utc_hours,era,notes\n"
                     "2025-03-06,c2,1,CACO04,2025-03-03,2025-03-09,CACO04_c2_20250219_EO.yaml,,13.5-18,adcp,"
                     "\"c2 only, CACO04\",\n"
                     "2025-03-07,c2,1,CACO04,2025-03-03,2025-03-09,CACO04_c2_20250219_EO.yaml,,13.5-18,adcp,"
                     "c2 again,, \n")
    rows, out = quiet(sp.read_table, trail, sp.CONFIG_FIELDS)
    check("WARNING" not in out and [r["notes"] for r in rows] == ["c2 only, CACO04", "c2 again"],
          "a row ending in empty extra fields (a trailing comma, no text after it): read whole, no warning")
    short = d / "short_surveys.csv"           # the optional last column left off: read as empty, no warning
    short.write_text("date,name,path,survey_type,survey_date,label,why,label_override_reason\n"
                     "2026-09-29,rtk_b,2026-09-29_b.csv,points,2026-09-29,INDEPENDENT,second set\n")
    rows, out = quiet(sp.read_table, short, sp.SURVEY_FIELDS)
    check("WARNING" not in out and rows[0]["why"] == "second set" and rows[0]["label_override_reason"] == "",
          "a row without its optional last column (label_override_reason): read as empty, no warning")
    # surveys/README.md: in such a row a comma outside double quotes ends `why` with no warning, and the rest
    # becomes the row's label_override_reason (test_station_tables lists every one the station's tables hold)
    for why, reason, said in (("second GCP set, surveyed by Chris", "surveyed by Chris", "outside double quotes"),
                              ('"second GCP set, surveyed by Chris"', "", "inside double quotes")):
        short.write_text("date,name,path,survey_type,survey_date,label,why,label_override_reason\n"
                         f"2024-10-23,oct_gcps_b,,points,2024-10-23,INDEPENDENT,{why}\n")
        rows, out = quiet(sp.read_table, short, sp.SURVEY_FIELDS)
        check("WARNING" not in out and rows[0]["label_override_reason"] == reason
              and rows[0]["why"] == ("second GCP set" if reason else "second GCP set, surveyed by Chris"),
              f"that row with a comma in `why` {said}: " + (f"no warning, the rest read as label_override_reason "
                                                            f"{reason!r} (why the README asks for the quotes)"
                                                            if reason else "read whole, no override"))


def git_out(*words, repo=None):
    """git run on the code's checkout (or `repo`) -> (stdout text, "") or (None, why not)."""
    repo = HERE if repo is None else repo
    try:
        r = subprocess.run(["git", "-C", str(repo)] + list(words), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           timeout=20)
    except (OSError, subprocess.SubprocessError) as exc:
        return None, f"git cannot be run here ({exc})"
    if r.returncode != 0:
        err = r.stderr.decode("utf-8", "replace").strip().splitlines()
        return None, f"git cannot read a checkout at {repo}" + (f" ({err[0][:120]})" if err else "")
    return r.stdout.decode("utf-8", "replace"), ""


def shipped_commit(repo=None):
    """The commit the repository shipped this test in: the last one that changed test_survey_products.py.
    -> (commit, "") or (None, why not). The tables are read as of that commit, not of HEAD: on the station
    the owner edits them, and may commit the edits to pull past them; neither reaches the checks of what
    was asked for (a merge that brings upstream code in leaves this file as upstream wrote it)."""
    out, why = git_out("log", "-1", "--format=%H", "--", "test_survey_products.py", repo=repo)
    if out is None or not out.strip():
        return None, why or f"git has no commit of test_survey_products.py in a checkout at {repo or HERE}"
    return out.strip(), ""


TABLES = ("surveys/survey_dates.csv", "surveys/surveys.csv")


def tables_changed_after(commit, repo=None):
    """The tables at HEAD against those of `commit` -> (names that differ, the later commits that changed
    them, newest first), or (None, why not). A commit that changes the shipped tables without this test leaves
    test_shipped_tables checking the older ones: it says so."""
    names, why = git_out("diff", "--name-only", commit, "HEAD", "--", *TABLES, repo=repo)
    if names is None:
        return None, why
    later, why = git_out("log", "--format=%h", f"{commit}..HEAD", "--", *TABLES, repo=repo)
    if later is None:
        return None, why
    return sorted(names.split()), later.split()


def test_tables_changed_note(d):
    """tables_changed_after on a git repository of the test's own (the code's checkout is not touched): a
    commit that changes only the tables after the test's last commit is found, and one that changes the
    test with them is not."""
    print("the shipped-tables checks say when the tables changed after this test (a git repository of its own)")
    repo = Path(d) / "tables_git"
    (repo / "surveys").mkdir(parents=True)
    home = Path(d) / "tables_git_home"            # git's own settings: none of the computer's (hooks, signing)
    home.mkdir()
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update(HOME=str(home), XDG_CONFIG_HOME=str(home), GIT_CONFIG_NOSYSTEM="1")

    def commit(msg, **files):
        for name, text in files.items():
            (repo / name).write_text(text)
        for words in (["add", "-A"], ["-c", "user.name=test", "-c", "user.email=test@invalid", "-c",
                                      "commit.gpgsign=false", "commit", "-q", "-m", msg]):
            subprocess.run(["git", "-C", str(repo)] + words, env=env, check=True, timeout=20,
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    try:
        subprocess.run(["git", "init", "-q", str(repo)], env=env, check=True, timeout=20,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        commit("shipped", **{"test_survey_products.py": "v1\n", "surveys/survey_dates.csv": "a\n",
                             "surveys/surveys.csv": "b\n"})
        shipped = shipped_commit(repo)[0]
        commit("code only", **{"survey_products.py": "x\n"})
        same = tables_changed_after(shipped, repo)
        commit("tables only", **{"surveys/survey_dates.csv": "a2\n"})
        later = tables_changed_after(shipped, repo)
        commit("test with them", **{"test_survey_products.py": "v2\n", "surveys/surveys.csv": "b2\n"})
        shipped2 = shipped_commit(repo)[0]
        after = tables_changed_after(shipped2, repo)
    except (OSError, subprocess.SubprocessError) as exc:
        err = getattr(exc, "stderr", None)
        skip("tables-changed note check", f"git cannot make a repository here ({exc}"
             + (f": {err.decode('utf-8', 'replace').strip()[:120]}" if err else "") + ")")
        return

    def got(r):                                   # without the commit ids, which differ from run to run
        return r if r[0] is None else f"{r[0]}, {len(r[1])} commit(s)"
    check(same == ([], []), f"a later commit of the code only: nothing said (got {got(same)})")
    check(later[0] == ["surveys/survey_dates.csv"] and len(later[1]) == 1,
          f"then a commit of the tables only: the table named, with that one commit (got {got(later)})")
    check(shipped2 != shipped and after == ([], []),
          f"the test changed with the tables: its new commit is the shipped one, nothing said (got {got(after)})")


def test_shipped_tables(d):
    """surveys/*.csv as the repository shipped them with this test (git, at the last commit of this file),
    not as the station has them: the dates, windows, envelopes, calibrations and labels this branch was
    asked for. On the station or in a clean checkout the same files are checked; the owner's later edits,
    committed or not, are test_station_tables' business."""
    print("the tables as the repository shipped them with this test (git): what was asked for")
    d = Path(d) / "shipped"
    d.mkdir()
    commit, why = shipped_commit()
    if commit is None:
        skip("shipped-tables checks", why)
        return
    print(f"  (surveys/*.csv and calibration/ as of commit {commit[:12]}, the last commit of this test)")
    changed, later = tables_changed_after(commit)
    if changed is None:
        print(f"  (whether the tables changed after it is not known: {later})")
    elif changed:
        print(f"  ({', '.join(changed)} changed after it, in {len(later)} later commit(s) ({', '.join(later[:6])}"
              f"{', ...' if len(later) > 6 else ''}): what was asked for is checked below on the tables as of "
              f"{commit[:12]}, not on those; a commit that changes the shipped tables changes this test with them. "
              f"The tables as they are now are checked for what the product needs, last)")
    texts = {}
    for name in ("survey_dates.csv", "surveys.csv"):
        out, why = git_out("show", f"{commit}:./surveys/{name}")
        if out is None:
            skip("shipped-tables checks", why)
            return
        texts[name] = out
        (d / name).write_text(out, encoding="utf-8")
    tracked, why = git_out("ls-tree", "--name-only", commit, "calibration/")
    if tracked is None:
        skip("shipped-tables checks", why)
        return
    in_repo = {Path(x).name for x in tracked.splitlines()}
    for name in texts:
        rows = list(csv.reader(io.StringIO(texts[name], newline="")))
        rows = [r for r in rows if r and not r[0].startswith("#")]
        bad_rows = [i + 1 for i, r in enumerate(rows) if len(r) != len(rows[0])]
        check(not bad_rows, f"surveys/{name}: every row has its {len(rows[0])} columns (a comma in a text is "
                            f"inside double quotes)" + (f": NOT rows {bad_rows}" if bad_rows else ""))
    cfg = sp.read_table(d / "survey_dates.csv", sp.CONFIG_FIELDS)
    srv = sp.read_table(d / "surveys.csv", sp.SURVEY_FIELDS)
    tabs = {"config": str(d / "survey_dates.csv"), "surveys": str(d / "surveys.csv")}
    for date in ("2024-10-23", "2025-03-19"):
        plan, out = quiet(sp.build_plan, date, cfg, srv, cfg_args(**tabs))
        forced, _ = quiet(sp.build_plan, date, cfg, srv, cfg_args(force_disabled=True, **tabs))
        check(plan is None and "DISABLED" in out and forced is not None,
              f"{date}: switched off (the user asked to skip it): refused, and builds with --force-disabled")
    plan, _ = quiet(sp.build_plan, "2025-01-23", cfg, srv, cfg_args(**tabs))
    check(plan["last"] == "2025-01-23" and all(c["station"] == "CACO03" for c in plan["cams"].values()),
          "2025-01-23: CACO03 only, 24 Jan (the re-set) outside the window")
    plan, _ = quiet(sp.build_plan, "2025-03-06", cfg, srv, cfg_args(**tabs))
    check(list(plan["cams"]) == ["c2"] and plan["cams"]["c2"]["station"] == "CACO04", "2025-03-06: c2 only, CACO04")
    plan, _ = quiet(sp.build_plan, "2026-09-29", cfg, srv, cfg_args(**tabs))
    check(plan["first"] == "2026-09-26" and plan["last"] == "2026-10-02" and plan["era"] == "live"
          and [s["name"] for s in plan["surveys"]] == ["rtk_2026-09-29"]
          and not [r for r in cfg if r["date"] == "2026-09-27"],
          "2026-09-29 (the user first called it 27 Sep): window 26 Sep - 2 Oct, the 29 Sep RTK; no 2026-09-27 product")
    rtk_row = [s for s in srv if s["name"] == "rtk_2026-09-29"][0]
    check(rtk_row["date"] == rtk_row["survey_date"] == "2026-09-29" and rtk_row["label"] == "INDEPENDENT"
          and "CIRCULAR" in rtk_row["why"] and "C = 0.037" in rtk_row["why"]
          and "calm low-tide RTK" in rtk_row["why"] and "+1.2 m" in rtk_row["why"],
          "the 2026 RTK row: surveyed 29 Sep, the 2026-09-29 product; INDEPENDENT as a table label, saying the "
          "setup rule makes it CIRCULAR with C = 0.037, where the check is and what an independent check needs")
    for r in cfg:
        if r["era"] == "live":
            continue
        for s in srv:
            if s["date"] == r["date"] and s["path"]:
                check(Path(s["path"]).name != Path(r["envelope_survey"]).name,
                      f"{r['date']} {r['camera']}: envelope {r['envelope_survey'][:28]}... is not the compared survey")
    missing = sorted({r["eo_file"] for r in cfg if r["eo_file"] not in in_repo})
    check(not missing, f"every calibration the {len(cfg)} rows name is in the repository's calibration/"
                       + (f": NOT {', '.join(missing)}" if missing else ""))
    # points surveyed on a calibration's own day may be the GCPs it was solved from: the table says so itself
    gcps = [s for s in srv if s["survey_type"] == "points" and not s.get("label_override_reason")
            and any(calibration_day(r["eo_file"]) == s["survey_date"] for r in cfg if r["date"] == s["date"])]
    check(all(s["label"] in ("PARTLY-CIRCULAR", "CIRCULAR") for s in gcps),
          f"GCP rows (points of their calibration's day: {', '.join(s['name'] for s in gcps) or 'none'}): "
          f"PARTLY-CIRCULAR in the table, whether or not their file has arrived")
    oct_ = [s for s in srv if s["name"] == "oct_gcps"][0]
    check(oct_["label"] == "PARTLY-CIRCULAR" and not oct_["path"], "Oct GCPs: PARTLY-CIRCULAR, path empty until the "
                                                                   "file arrives")


def test_station_tables(d):
    """The station's own tables as they are on this computer, which the owner edits (surveys/README.md: add a
    date, switch one on or off, give a survey file that has arrived). Only what survey_products itself needs
    of them, or does with them, is asserted, so an edit the product reads with no warning and accepts never
    fails here (surveys/README.md: a text with a comma goes inside double quotes, or the product cuts it); what
    was asked for in the tables the repository ships is test_shipped_tables'. Every problem is listed, not just the
    first, and this group runs last: a table problem cannot hide the results of the code's checks."""
    print("the station's own tables (surveys/*.csv as they are now): what survey_products needs of them")
    problems = []

    def tcheck(cond, msg):
        if cond:
            print(f"  ok  {msg}")
        else:
            print(f"  FAIL  {msg}")
            problems.append(msg)
    tables = {}
    for path, fields in ((sp.DEFAULT_CONFIG, sp.CONFIG_FIELDS), (sp.DEFAULT_SURVEYS, sp.SURVEY_FIELDS)):
        # read as the product reads them (its warning is printed once per line and run: forget earlier reads)
        sp._TABLE_WARNED.difference_update({k for k in sp._TABLE_WARNED if k[0] == str(Path(path))})
        err = ""
        try:
            tables[path], out = quiet(sp.read_table, path, fields)
        except SystemExit as exc:           # the product's own refusal (no such file, a column missing)
            tables[path], out, err = None, "", str(exc)
        except Exception as exc:            # e.g. not UTF-8 (a spreadsheet saving Latin-1), or a csv.Error: the
            tables[path], out, err = None, "", f"{type(exc).__name__}: {exc}"     # product cannot read it either
        if err:
            tcheck(False, f"surveys/{Path(path).name}: survey_products cannot read it: {err}")
            continue
        warned = [ln.split(": ", 1)[-1] for ln in out.splitlines() if "WARNING" in ln]
        tcheck(not warned, f"surveys/{Path(path).name}: read by survey_products with no warning"
                           + (f": {'; '.join(warned)}" if warned else ""))
    cfg, srv = tables[sp.DEFAULT_CONFIG], tables[sp.DEFAULT_SURVEYS]
    if cfg is None or srv is None:
        check(False, f"the station's tables: {len(problems)} problem(s): {' | '.join(problems)}")
    # a label_override_reason overrules the downgrades that can misfire, and in a row without its last column
    # the end of a `why` cut at a comma outside double quotes lands there unwarned: each one is shown
    over = [s for s in srv if s.get("label_override_reason")]
    print(f"  ({len(over)} survey row(s) with a label_override_reason, which overrules the calibration-day, "
          f"setup-fit-date and survey_compare name/date downgrades"
          + (": each should be meant, not the end of a `why` cut at a comma outside double quotes)" if over else ")"))
    for s in over:
        print(f"      {s['date']} {s['name']}: label_override_reason {s['label_override_reason']!r}; "
              f"its why ends {s['why'][-60:]!r}")
    dates = sorted({r["date"] for r in cfg})
    off = [x for x in dates if not any(is_enabled(r) for r in cfg if r["date"] == x)]
    plans = {}
    for x in dates:
        try:
            plan, out = quiet(sp.build_plan, x, cfg, srv, cfg_args())
            forced = quiet(sp.build_plan, x, cfg, srv, cfg_args(force_disabled=True))[0] if x in off else plan
        except SystemExit as exc:
            tcheck(False, f"{x}: its rows cannot be planned ({exc})")
            continue
        plans[x] = forced
        if x in off:
            tcheck(plan is None and "DISABLED" in out and forced is not None,
                   f"{x}: switched off: refused, and builds with --force-disabled")
        else:
            tcheck(plan is not None and plan["cams"], f"{x}: switched on: its plan is built")
    for r in cfg:
        if is_enabled(r):
            tcheck(bool(sp.resolve_file(r["eo_file"], [str(HERE / "calibration")])),
                   f"{r['date']} {r['camera']} (switched on): its calibration {r['eo_file']} is in calibration/")
    # The product's label rules on these rows (survey_products.honest_label): a survey the search envelope was
    # placed with, and a points survey made on its calibration's day (maybe its GCPs), can be no better than
    # PARTLY-CIRCULAR; the second rule is overruled only by a written label_override_reason. Stand-in files
    # (each its own content) let the envelopes and surveys be found by their names alone.
    stand = Path(d) / "station_tables"
    stand.mkdir()
    for n in {Path(r["envelope_survey"]).name for r in cfg if r["envelope_survey"]} | \
             {Path(s["path"]).name for s in srv if s["path"] and not any(ch in s["path"] for ch in "*?[")}:
        (stand / n).write_text(f"stand-in for {n}\n")
    n_rule = 0
    for x, plan in plans.items():
        if plan is None:
            continue
        # the cameras a build of the date uses: those switched on, or every one of a date switched off (forced)
        a = cfg_args(survey_dirs=[str(stand)], force_disabled=x in off)
        plan = quiet(sp.build_plan, x, cfg, srv, a)[0]
        for s in plan["surveys"]:
            v = sp.Verdict(s)
            sp.honest_label(s, plan, a, sp.survey_path_for(s, plan, a), fits=[], verdict=v)
            worse = sp.LABELS.index(v.label) >= sp.LABELS.index("PARTLY-CIRCULAR")
            env = [cam for cam, c in plan["cams"].items() if s["path"] and c["envelope_survey"]
                   and Path(c["envelope_survey"]).name == Path(s["path"]).name]
            if env:
                n_rule += 1
                tcheck(worse and any("search envelope was placed with this survey" in t for t in v.reasons),
                       f"{x} {s['name']}: the {'/'.join(env)} search envelope was placed with it: "
                       f"{v.label} (PARTLY-CIRCULAR or worse)")
            if s["survey_type"] == "points" and any(calibration_day(c["eo_file"]) == s.get("survey_date")
                                                    for c in plan["cams"].values()):
                n_rule += 1
                said = [t for t in v.reasons + v.overridden if "was solved on the survey day" in t]
                if v.override:
                    tcheck(said, f"{x} {s['name']}: points of its calibration's day, label_override_reason "
                                 f"written: the calibration-day downgrade overruled and listed ({v.label})")
                else:
                    tcheck(worse and said, f"{x} {s['name']}: points of its calibration's day (maybe its GCPs): "
                                           f"{v.label} (PARTLY-CIRCULAR or worse)")
    print(f"  ({n_rule} survey row(s) of these tables concerned by the envelope and calibration-day rules)")
    check(not problems, f"the station's tables: {len(problems)} problem(s): {' | '.join(problems)}"
          if problems else f"the station's tables: {len(cfg)} date rows and {len(srv)} survey rows, no problem")


def test_frames_per_day(d):
    print("frames per day (README caveats)")
    f = Path(d) / "fpd.csv"
    f.write_text("source_file,camera,capture_time_utc,x\n"
                 "a,c1,2025-01-18T14:00:00Z,1\na,c1,2025-01-18T14:00:00Z,2\nb,c1,2025-01-20T15:00:00Z,1\n"
                 "c,c2,2025-01-20T15:00:00Z,1\n")
    check(sp.frames_per_day(f) == {"c1": {"2025-01-18": 1, "2025-01-20": 1}, "c2": {"2025-01-20": 1}},
          "frames counted once per photo, per camera and day")
    check(sp.frames_per_day(Path(d) / "missing.csv") == {}, "no waterline file yet: nothing")
    # the detector's scratch goes once the waterlines are written; the waterline CSVs stay
    c = Path(d) / "prune" / "c1"
    for sub in ("debug/x", "in", "detections", "overlays"):
        (c / sub).mkdir(parents=True)
    (c / "debug" / "x" / "a.png").write_bytes(b"0" * 10)
    (c / "in" / "a.jpg").write_bytes(b"0" * 10)
    (c / "detections" / "a.csv").write_text("Column,Row\n")
    (c / "detections" / "a_overlay.png").write_bytes(b"0" * 10)
    (c / "overlays" / "a.jpg").write_bytes(b"0" * 10)
    quiet(sp.prune_detector_scratch, c)
    left = sorted(str(p.relative_to(c)) for p in c.rglob("*") if p.is_file())
    check(left == ["detections/a.csv", "overlays/a.jpg"], f"detector scratch removed, waterlines kept ({left})")


def test_caveats(d):
    print("README caveats")
    out = Path(d) / "cav"
    out.mkdir()
    a = SimpleNamespace(setup_coef=sp.SETUP_COEF)
    for era in ("live", "adcp"):
        plan = {"era": era, "date": "2026-09-29", "first": "2026-09-26", "last": "2026-10-02", "out": out,
                "surveys": [], "cams": {"c1": {"first": "2026-09-26", "last": "2026-10-02"}}}
        forcing = {"waves": {"setup_in_window": {"median_m": 0.3, "daytime_median_m": 0.28}}}
        cav = sp.collect_caveats(plan, a, {}, forcing, {"frames_per_day": {"c1": {"2026-09-29": 3}}})
        text = "\n".join(cav)
        first = "Live era" if era == "live" else "ADCP datum"
        check(first in text, f"{era}: its own era caveat is there")
        check("is not pinned by the 2026 RTK" in text and "under-corrects" not in text,
              f"{era}: the caveat on what the 2026 RTK can say about C is there too")
        check("none on 2026-09-26" in text, f"{era}: the window days without waterlines are listed")
        check(("GNSS-R" in text and "surf zone" in text) == (era == "adcp"),
              f"{era}: C's still-water reference (GNSS-R, surf zone) stated for the ADCP era only")
        check("NOMINAL" in text and "6348" in text, f"{era}: the GeoTIFF CRS tag is said to be nominal")
    cav = sp.collect_caveats(plan, SimpleNamespace(setup_coef=0.05), {}, forcing, {})
    check(not any("not pinned by the 2026 RTK" in c for c in cav) and any("WITHOUT saying where" in c for c in cav),
          "a C of unknown origin: no RTK caveat, and said to be of unknown origin")
    cav = sp.collect_caveats(plan, SimpleNamespace(setup_coef=0.05, setup_fitted_to="none:crossings"), {}, forcing, {})
    check(not any("fitted on this beach in Sep-Oct 2026" in c for c in cav),
          "a declared C is not described as the 2026 fit")
    check(sp.stat_text({"n": 1, "median": -0.3, "nmad": 0.0, "rmse": 0.3}).count("NMAD -") == 1
          and "TOO FEW" in sp.stat_text({"n": 1, "median": -0.3, "nmad": 0.0, "rmse": 0.3}),
          "one value: no NMAD, and marked too few")


def test_stamps(d):
    print("resumable steps")
    d = Path(d) / "stamps"
    d.mkdir()
    inp, outp = d / "in.csv", d / "out.csv"
    inp.write_text("a\n")
    outp.write_text("b\n")
    sig = {"cmd": "x --y 1"}
    check(not sp.is_fresh(d, "k", sig, [outp], [inp])[0], "never built: not fresh")
    sp.write_stamp(d, "k", sig, [["x", "--y", 1]])
    check(sp.is_fresh(d, "k", sig, [outp], [inp])[0], "built with these settings: fresh")
    check(not sp.is_fresh(d, "k", {"cmd": "x --y 2"}, [outp], [inp])[0], "other settings: rebuilt")
    t = outp.stat().st_mtime
    os.utime(inp, (t + 10, t + 10))
    fresh, why = sp.is_fresh(d, "k", sig, [outp], [inp])
    check(not fresh and "newer" in why, "an input newer than the output: rebuilt")
    outp.unlink()
    check(not sp.is_fresh(d, "k", sig, [outp], [])[0], "a missing output: rebuilt")
    # chained: a rebuilt upstream step makes the filtered file stale
    plan = {"out": d, "date": "x"}
    w = d / "waterlines"
    w.mkdir()
    (w / "contour_points_ground.csv").write_text("a\n")
    sp.write_stamp(d, "merge", {"m": 1})
    (w / "contour_points_ground_filtered.csv").write_text("a\n")
    sp.write_stamp(d, "filter", {"after": {"merge": sp.token(d, "merge")}})
    check(sp.waterline_file(plan).name == "contour_points_ground_filtered.csv", "filtered file of the current merge: used")
    sp.write_stamp(d, "merge", {"m": 2})
    f, out = quiet(sp.waterline_file, plan)
    check(f.name == "contour_points_ground.csv" and "NOT used" in out,
          "filtered file of an earlier merge: NOT used, and said so")
    check("waterline_consistency.py" in sp.scripts_sig("filter") and sp.scripts_sig("filter")["waterline_consistency.py"],
          "the script's sha256 is part of the step's signature")
    st = sp.read_stamp(d, "merge")
    check(st["software"]["scripts"].get("survey_products.py") == sp.script_sha("survey_products.py"),
          "each stamp records the sha256 of the scripts that built it")
    # the commit: whether git can read this code's checkout is the computer's business (a copy that is not
    # a checkout, a user git refuses: 'dubious ownership'), so it is checked when it can be, and said when not
    try:
        r = subprocess.run(["git", "-C", str(HERE), "rev-parse", "HEAD"], stdout=subprocess.PIPE,
                           stderr=subprocess.DEVNULL, universal_newlines=True, timeout=20)
        head = r.stdout.strip() if r.returncode == 0 else ""
    except (OSError, subprocess.SubprocessError):
        head = ""
    if len(head) == 40:
        check(st["software"]["commit"] == head, "each stamp records the commit of the code that built it")
    else:
        check(st["software"]["commit"] is None, "no commit here: the stamp says none (not a made-up one)")
        skip("commit check", f"git cannot read a checkout at {HERE}: the stamps record no commit, only the sha256s")


def test_keep_moved_days(d):
    print("--keep-moved-days is applied when the detection runs, not stored")
    res = {"different_days": ["2025-01-24"], "leave_out": True}
    check(sp.skip_days_of(res, SimpleNamespace(keep_moved_days=False)) == ["2025-01-24"], "default: the day is left out")
    check(sp.skip_days_of(res, SimpleNamespace(keep_moved_days=True)) == [], "--keep-moved-days on a rerun: kept")
    check(sp.skip_days_of(dict(res, leave_out=False), SimpleNamespace(keep_moved_days=False)) == [],
          "most days differ: nothing left out automatically")


# -- pointing verdicts (Oct 2026): what the step concludes from the two methods' measurements --------
# Edge numbers (horizon_prominence, horizon_rival, horizon_sky_above, horizon_to_limit_px) measured with
# survey_products.horizon_edge on the REAL frames: the 17:00 UTC timex of 26-28 Sep 2026 sent from the
# NUC, and 7 Oct 2026 16:00 (clear). The fixture's rendered horizon: prominence 14.7-24, rival 0.49-0.62.
EDGE_CLEAR = {"c1": (73.1, 0.08, 1.0, 207), "c2": (84.9, 0.06, 1.0, 214)}
EDGE_STORM = {("c1", "2026-09-26"): (6.3, 0.70, 0.00, 43), ("c1", "2026-09-27"): (6.0, 0.68, 0.00, 15),
              ("c1", "2026-09-28"): (2.2, 0.96, 0.89, 141), ("c2", "2026-09-26"): (6.1, 0.72, 0.00, 62),
              ("c2", "2026-09-27"): (6.7, 0.57, 0.22, 25), ("c2", "2026-09-28"): (2.7, 0.92, 1.00, 172)}
EDGE_FIXTURE = (18.0, 0.55, 1.0, 240)
# The NUC's own 2026-09-29 pointing table (python3 survey_products.py --date 2026-09-29, Oct 2026):
# date: (horizon, offset px, tilt, roll, fit px, features, inliers, pan, tilt, roll vs the reference)
NUC_TABLE = {
    "c1": {"2026-09-26": ("off the calibration", 205.8, 4.227, 0.061, 1.6, "unmatched", 0, "", "", ""),
           "2026-09-27": ("off the calibration", 230.8, 4.818, 0.113, 5.0, "MOVED", 48, -1.214, 0.004, 0.03),
           "2026-09-28": ("off the calibration", 83.4, 1.716, 0.952, 3.5, "unmatched", 3, "", "", ""),
           "2026-09-29": ("off the calibration", 44.7, 0.956, 0.158, 0.3, "ok", 650, 0.0, -0.0, -0.005),
           "2026-09-30": ("off the calibration", 43.7, 0.949, 0.245, 0.4, "ok", 423, -0.007, -0.002, 0.001),
           "2026-10-01": ("off the calibration", 41.4, 0.908, 0.235, 0.3, "ok", 364, -0.003, -0.004, 0.005),
           "2026-10-02": ("off the calibration", 41.4, 0.899, 0.255, 0.4, "ok", 134, -0.005, -0.008, 0.019)},
    "c2": {"2026-09-26": ("off the calibration", 186.6, 2.624, 3.843, 4.3, "ok", 73, -0.02, 0.063, 0.071),
           "2026-09-27": ("off the calibration", 214.2, 3.248, 0.741, 3.3, "ok", 375, -0.009, 0.003, 0.017),
           "2026-09-28": ("unclear", 76.8, 0.731, -1.29, 40.3, "ok", 393, 0.008, -0.002, -0.01),
           "2026-09-29": ("off the calibration", 35.0, 0.567, -0.62, 0.2, "ok", 2723, -0.0, 0.002, -0.0),
           "2026-09-30": ("off the calibration", 34.0, 0.557, -0.66, 0.2, "ok", 1400, -0.004, 0.002, -0.009),
           "2026-10-01": ("off the calibration", 31.0, 0.521, -0.7, 0.2, "ok", 1341, -0.001, 0.002, -0.006),
           "2026-10-02": ("off the calibration", 31.0, 0.526, -0.69, 0.2, "ok", 610, -0.002, -0.003, -0.005)}}


def _pt_row(cam, date, horizon="matches", off=8.0, dt=0.15, dr=0.0, fit=0.4, features="ok", inliers=500,
            d=(0.0, 0.0, 0.0), edge=None):
    """A pointing row as check_pointing writes it (edge: the four edge numbers, or "" for none)."""
    r = {"camera": cam, "date": date, "photos": 31, "horizon": horizon, "horizon_offset_px": off,
         "horizon_dtilt_deg": dt, "horizon_droll_deg": dr, "horizon_fit_px": fit, "horizon_frames": 3,
         "features": features, "inliers": inliers, "d_azimuth_deg": "", "d_tilt_deg": "", "d_roll_deg": "",
         "runner_up_inliers": "", "runner_up_d": ""}
    if features in ("ok", "moved", "MOVED", "ambiguous"):
        r.update(d_azimuth_deg=d[0], d_tilt_deg=d[1], d_roll_deg=d[2])
    if edge != "":
        e = edge or EDGE_CLEAR.get(cam, EDGE_FIXTURE)
        r.update(zip(sp.EDGE_KEYS, e))
    return r


def _nuc_rows(cam):
    rows = []
    for day, (hz, off, dt, dr, fit, feat, n, a, t, rl) in NUC_TABLE[cam].items():
        rows.append(_pt_row(cam, day, hz, off, dt, dr, fit, feat, n, (a, t, rl),
                            edge=EDGE_STORM.get((cam, day), EDGE_CLEAR[cam])))
    return rows


def test_pointing_verdicts(d):
    print("pointing verdicts: matched features decide, a misread horizon is no clear day, UNCONFIRMED left out")
    by = lambda rows: {r["date"]: r for r in rows}          # noqa: E731
    # -- the NUC's 2026-09-29 table, c1: the storm days are left out with an honest reason
    c1 = _nuc_rows("c1")
    old_med = float(np.median([r["horizon_dtilt_deg"] for r in c1]))          # the old 'clear': every fitted day
    i1 = sp.pointing_verdicts(c1)
    v = by(c1)
    check(v["2026-09-26"]["verdict"] == "UNCONFIRMED" and v["2026-09-28"]["verdict"] == "UNCONFIRMED"
          and all("could not be confirmed" in v[x]["why"] and "not a clear horizon" in v[x]["why"]
                  and "DIFFERENT" not in v[x]["why"] for x in ("2026-09-26", "2026-09-28"))
          and "(0 inliers; 40 are needed" in v["2026-09-26"]["why"] and "(3 inliers" in v["2026-09-28"]["why"],
          "NUC c1 26/28 Sep: a horizon jump on a misread edge with no matched features: UNCONFIRMED (not 'moved')")
    w27 = v["2026-09-27"]["why"]
    check(v["2026-09-27"]["verdict"] == "DIFFERENT POINTING" and w27.startswith("48 matched features show a turn of "
                                                                                 "-1.21/+0.00/+0.03 deg")
          and "+4.82/+0.11" in w27 and "not a clear horizon" in w27 and "not used" in w27,
          "NUC c1 27 Sep: -1.21 deg of pan on 48 inliers (>= 40: good to ~0.1 deg) is DIFFERENT POINTING, from the "
          "features; the horizon's +4.82 deg is named a misread")
    check(all(v[x]["verdict"] == "matches" and not v[x]["why"] for x in ("2026-09-29", "2026-09-30", "2026-10-01",
                                                                           "2026-10-02")),
          "NUC c1 29 Sep - 2 Oct: matches")
    check(i1["clear_horizon_days"] == 4 and abs(i1["window_horizon_dtilt_deg"] - 0.928) < 0.0011
          and abs(i1["window_horizon_droll_deg"] - 0.240) < 0.0011 and abs(i1["window_horizon_offset_px"] - 42.55) < 0.06
          and i1["horizon_not_clear_days"] == ["2026-09-26", "2026-09-27", "2026-09-28"] and abs(old_med - 0.956) < 1e-9,
          f"NUC c1 window: the 4 clear days only (+0.93/+0.24 deg, 42.5 px), not the old 7 with the storm misreads "
          f"(median tilt {old_med:+.3f})")
    lo = sp.pointing_leave_out(c1)
    check(lo["leave_out"] and lo["different_days"] == ["2026-09-26", "2026-09-27", "2026-09-28"]
          and lo["measured_different_days"] == ["2026-09-27"] and lo["unconfirmed_days"] == ["2026-09-26", "2026-09-28"],
          "NUC c1: 26-28 Sep left out (1 measured, 2 unconfirmed); the 'most days differ' guard stays off")
    # -- c2: the features overrule the misread horizon; nothing left out
    c2 = _nuc_rows("c2")
    i2 = sp.pointing_verdicts(c2)
    v = by(c2)
    w26, w27 = v["2026-09-26"]["why"], v["2026-09-27"]["why"]
    check(v["2026-09-26"]["verdict"] == "matches" and v["2026-09-27"]["verdict"] == "matches"
          and "+2.62/+3.84" in w26 and "+0.54/-0.68" in w26 and "73 matched features show at most 0.07 deg" in w26
          and "+3.25/+0.74" in w27 and "375 matched features show at most 0.02 deg" in w27
          and all("not a clear horizon" in w and "the features decide, day kept" in w for w in (w26, w27)),
          "NUC c2 26/27 Sep: 73 and 375 matched features (<= 0.07 deg) overrule the horizon's 2.6-3.3 deg; the why "
          "gives both numbers and the window's")
    check(all(r["verdict"] == "matches" for r in c2) and i2["horizon_overruled_days"] == ["2026-09-26", "2026-09-27"]
          and not sp.pointing_leave_out(c2)["different_days"],
          "NUC c2: every day matches (28 Sep: its horizon fit failed, 393 features decide); none left out")
    check(i2["clear_horizon_days"] == 4 and abs(i2["window_horizon_dtilt_deg"] - 0.5415) < 0.0011
          and abs(i2["window_horizon_droll_deg"] + 0.675) < 0.0011 and abs(i2["window_horizon_offset_px"] - 32.5) < 0.06,
          "NUC c2 window: the 4 clear days (+0.54/-0.68 deg, 32.5 px), not the old 6 (+0.56/-0.64)")
    # -- printed: every day's why below the table, the left-out days with their verdicts
    res1 = dict(sp.pointing_leave_out(c1), rows=c1, info=dict(i1, calibration="CACO05_c1_20251113_EO-CV.yaml"))
    _, txt = quiet(sp.print_pointing, "c1", res1, SimpleNamespace(keep_moved_days=False))
    check("LEFT OUT of the detection: 2026-09-26 (UNCONFIRMED), 2026-09-27 (DIFFERENT POINTING), 2026-09-28 "
          "(UNCONFIRMED)" in txt and "    2026-09-26 UNCONFIRMED: the horizon read +4.23/+0.06" in txt,
          "print: the left-out days with their verdicts, and each day's why")
    res2 = dict(sp.pointing_leave_out(c2), rows=c2, info=dict(i2, calibration="CACO05_c2_20251113_EO-CV.yaml"))
    _, txt = quiet(sp.print_pointing, "c2", res2, SimpleNamespace(keep_moved_days=False))
    check("    2026-09-27 matches: the horizon read +3.25/+0.74" in txt and "LEFT OUT" not in txt,
          "print: an overruled horizon is shown with its why; nothing left out")
    _, txt = quiet(sp.print_pointing, "c1", res1, SimpleNamespace(keep_moved_days=True))
    check("days with a different or unconfirmed pointing KEPT: 2026-09-26 (UNCONFIRMED)" in txt
          and sp.skip_days_of(res1, SimpleNamespace(keep_moved_days=True)) == []
          and sp.skip_days_of(res1, SimpleNamespace(keep_moved_days=False)) == ["2026-09-26", "2026-09-27", "2026-09-28"],
          "--keep-moved-days keeps DIFFERENT POINTING and UNCONFIRMED days alike")
    # -- c1 27 Sep as a weak match: the same turn on 9 inliers (against the 7 Oct frame) or a reading on 39
    rows = _nuc_rows("c1")
    by(rows)["2026-09-27"].update(features="unmatched", inliers=9, d_azimuth_deg="", d_tilt_deg="", d_roll_deg="")
    sp.pointing_verdicts(rows)
    w = by(rows)["2026-09-27"]
    check(w["verdict"] == "UNCONFIRMED" and "(9 inliers; 40 are needed for a reading)" in w["why"]
          and "could not be confirmed" in w["why"],
          "c1 27 Sep matched on 9 features only (its -1.20 deg not trusted) + a misread horizon jump: UNCONFIRMED")
    rows = _nuc_rows("c1")
    by(rows)["2026-09-27"].update(inliers=39)
    sp.pointing_verdicts(rows)
    w = by(rows)["2026-09-27"]
    check(w["verdict"] == "UNCONFIRMED" and "only 39 features matched the reference" in w["why"],
          "a turn read on 39 inliers (below the 40 a reading needs) + a horizon jump: UNCONFIRMED, never measured")
    # -- a window with MORE misread days than clear ones: the clear days are not flagged, the offset stays theirs
    rows = [_pt_row("c1", f"2026-09-2{k}", "off the calibration", 200.0 + 5 * k, t, 0.1, 3.0, "unmatched", 2,
                    edge=(6.0, 0.7, 0.0, 30)) for k, t in zip(range(3, 7), (4.2, 4.8, 3.5, 4.5))]
    rows += [_pt_row("c1", f"2026-09-{dd}", "off the calibration", o, t, 0.24, 0.3, "unmatched", 5)
             for dd, o, t in (("27", 44.7, 0.95), ("28", 41.4, 0.91), ("29", 43.0, 0.90))]
    old = float(np.median([r["horizon_dtilt_deg"] for r in rows]))
    i = sp.pointing_verdicts(rows)
    v = by(rows)
    check(old > 3.0 and all(v[x]["verdict"] == "horizon only" for x in ("2026-09-27", "2026-09-28", "2026-09-29"))
          and all(v[f"2026-09-2{k}"]["verdict"] == "UNCONFIRMED" for k in range(3, 7))
          and i["clear_horizon_days"] == 3 and abs(i["window_horizon_dtilt_deg"] - 0.91) < 1e-9
          and abs(i["window_horizon_offset_px"] - 43.0) < 1e-9 and sp.pointing_leave_out(rows)["leave_out"],
          f"4 misread days against 3 clear: the clear days are not flagged (the old median {old:+.2f} would have "
          f"flagged them), the window is theirs (+0.91 deg, 43.0 px), the 4 are UNCONFIRMED and left out")
    # a clear-LOOKING jump (a fog front with a sharp lower edge) that 500 matched features at no turn
    # contradict: no clear day; the two methods disagree, and a zero turn may be dirt on the window (the
    # scene unmatched: no runner-up shows it), so the day is UNCONFIRMED; on the reference day itself
    # (the other days' features are measured against its frame) it is the horizon that is the odd one
    rows = [_pt_row("c2", f"2026-09-{dd}", "off the calibration", 33.0, t, -0.67, 0.2) for dd, t in
            (("29", 0.567), ("30", 0.557), ("01", 0.521), ("02", 0.526))]
    rows.append(_pt_row("c2", "2026-10-03", "off the calibration", 130.0, 2.0, -0.6, 0.6))
    i = sp.pointing_verdicts(rows)
    check(by(rows)["2026-10-03"]["verdict"] == "UNCONFIRMED" and "the two disagree" in by(rows)["2026-10-03"]["why"]
          and i["horizon_contradicted_days"] == ["2026-10-03"] and i["clear_horizon_days"] == 4
          and abs(i["window_horizon_offset_px"] - 33.0) < 1e-9 and not i["horizon_overruled_days"]
          and sp.pointing_leave_out(rows)["different_days"] == ["2026-10-03"],
          "a clear-looking horizon jump with 500 matched features at no turn: UNCONFIRMED (left out), no clear day")
    rows = [_pt_row("c2", f"2026-09-{dd}", "off the calibration", 33.0, t, -0.67, 0.2) for dd, t in
            (("30", 0.557), ("01", 0.521), ("02", 0.526))]
    rows.append(_pt_row("c2", "2026-09-29", "off the calibration", 130.0, 2.0, -0.6, 0.6))
    i = sp.pointing_verdicts(rows, ref_day="2026-09-29")
    check(by(rows)["2026-09-29"]["verdict"] == "matches" and "this is the reference day" in by(rows)["2026-09-29"]["why"]
          and i["horizon_overruled_days"] == ["2026-09-29"] and i["clear_horizon_days"] == 3,
          "... on the reference day itself: matches (its own reading is the odd one), and no clear day")
    # -- window dirt with the scene unmatched (review, Oct 2026: the real 7 Oct frames, the same drops on every
    # day, a 0.6 deg tilt on a fogged day: 68 inliers at 0.00 deg, a runner-up of 6, the horizon clear and
    # +0.60 deg off): the zero turn does not overrule the clear horizon; the old rule left the day out too
    rows = [_pt_row("c1", f"2026-09-{dd}", "off the calibration", 41.8, 0.91, 0.29, 0.4, "ok", 2840) for dd in
            ("26", "27", "28", "29", "30")] + [_pt_row("c1", "2026-10-01", "off the calibration", 41.8, 0.91, 0.29,
                                                       0.4, "ok", 2851)]
    rows.append(_pt_row("c1", "2026-10-02", "off the calibration", 70.8, 1.514, 0.31, 0.4, "ok", 68, (-0.0, 0.0, 0.0),
                        edge=(63.7, 0.08, 0.89, 178)))
    by(rows)["2026-10-02"].update(runner_up_inliers=6)
    i = sp.pointing_verdicts(rows, ref_day="2026-09-29")
    w = by(rows)["2026-10-02"]
    check(w["verdict"] == "UNCONFIRMED" and "here it has 6 inliers" in w["why"] and "+1.51/+0.31" in w["why"]
          and sp.skip_days_of(dict(sp.pointing_leave_out(rows), rows=rows), SimpleNamespace(keep_moved_days=False))
          == ["2026-10-02"] and i["clear_horizon_days"] == 6,
          "dirt with the scene unmatched (68 inliers at 0.00 deg, runner-up 6) against a clear +0.60 deg jump: "
          "UNCONFIRMED, left out (was 'matches', kept)")
    # -- a TURNED match (0.1-0.5 deg: the scene, not dirt) decides against a clear-looking horizon whose change
    # is not its own (review probe: 1500 inliers at +0.15 deg tilt, the horizon +2.00 deg): small change, kept;
    # when the two agree in size (0.35 / +0.35 below), both see the change
    base4 = [_pt_row("c2", f"2026-10-0{k}", "off the calibration", 33.0, 0.55, -0.67, 0.2, "ok", 1400) for k in range(1, 5)]
    for feat, n_, dtilt, want, kept in (("moved", 1500, 0.15, "small change", True), ("moved", 55, 0.11, "small change", True),
                                        ("ok", 1500, 0.09, "UNCONFIRMED", False)):
        rows = [dict(r) for r in base4] + [_pt_row("c2", "2026-10-05", "off the calibration", 165.0, 2.55, -0.6, 0.6,
                                                   feat, n_, (0.0, dtilt, 0.0))]
        i = sp.pointing_verdicts(rows)
        r = by(rows)["2026-10-05"]
        check(r["verdict"] == want and (("2026-10-05" in i["horizon_overruled_days"]) == kept)
              and i["horizon_contradicted_days"] == ["2026-10-05"] and i["clear_horizon_days"] == 4
              and (not kept or (f"where the features measured +{dtilt:.2f}/+0.00" in r["why"]
                                and "a change of +2.00/+0.07 deg" in r["why"])),
              f"{feat} {dtilt} deg on {n_} features against a clear-looking +2.00 deg jump: {want}"
              + (" (the turned match decides; the why gives both changes)" if kept else
                 " (a zero turn may be dirt: it cannot overrule a clear edge)"))
    # -- a real 0.5 deg move with strong features (exact rotations of real frames read 0.30/0.39/0.50 deg on
    # 64-106 inliers): measured, never overruled
    base = [_pt_row("c2", f"2026-10-0{k}", "off the calibration", 33.0, 0.55, -0.67, 0.2, "ok", 1400) for k in range(1, 5)]
    for feat, dtilt, hz_dt, edge, want in (
            ("MOVED", 0.51, 1.05, None, "DIFFERENT POINTING"),           # read 0.51: MOVED, a clear horizon too
            ("MOVED", 0.51, 1.05, (6.0, 0.7, 0.0, 30), "DIFFERENT POINTING"),   # ... a misread horizon
            ("MOVED", 0.51, 0.55, "", "DIFFERENT POINTING"),             # ... no horizon numbers at all
            ("moved", 0.49, 1.04, None, "DIFFERENT POINTING"),           # read 0.49: both methods see it
            ("moved", 0.49, 1.04, (6.0, 0.7, 0.0, 30), "small change"),  # the horizon misread: features decide
            ("moved", 0.35, 0.90, None, "DIFFERENT POINTING"),           # 0.35 deg on a clear horizon
            ("moved", 0.29, 0.84, None, "small change")):                # 0.29 deg: within the horizon's 0.3
        rows = [dict(r) for r in base] + [_pt_row("c2", "2026-10-05", "off the calibration", 33.0 + 66 * (hz_dt - 0.55),
                                                  hz_dt, -0.67, 0.2, feat, 1500, (0.0, dtilt, 0.0), edge=edge)]
        sp.pointing_verdicts(rows)
        r = by(rows)["2026-10-05"]
        check(r["verdict"] == want and ("1500 matched features" in r["why"]),
              f"a true tilt read {dtilt} deg ({feat}) on 1500 features, horizon "
              f"{'without edge numbers' if edge == '' else 'misread' if edge else 'clear'} at {hz_dt - 0.55:+.2f} deg: {want}")
    # -- a camera tilted 0.6 deg on half the clear days, features unmatched there (review scenario E: real
    # frames, the land fogged): the window's median falls between the two pointings (0.91/1.51/1.51/0.90 ->
    # 1.21) and saw nothing; judged against the reference (29 Sep, its own frames matched) both days are
    # measured, the one-day roll unconfirmed
    E = [_pt_row("c1", "2026-09-29", "off the calibration", 41.8, 0.908, 0.255, 0.4, "ok", 2804),
         _pt_row("c1", "2026-09-30", "off the calibration", 70.8, 1.506, 0.255, 0.4, "unmatched", 3),
         _pt_row("c1", "2026-10-01", "off the calibration", 70.8, 1.506, 0.255, 0.4, "unmatched", 3),
         _pt_row("c1", "2026-10-02", "off the calibration", 45.0, 0.895, 1.185, 0.4, "unmatched", 2)]
    storm3 = [_pt_row("c1", f"2026-09-2{k}", "off the calibration", 200.0, 4.2 + 0.3 * k, 0.1, 3.0, "unmatched", 2,
                      edge=(6.0, 0.7, 0.0, 30)) for k in (6, 7, 8)]
    rows = [dict(r) for r in storm3 + E]
    i = sp.pointing_verdicts(rows, ref_day="2026-09-29")
    v = by(rows)
    lo = sp.pointing_leave_out(rows)
    check([v[x]["verdict"] for x in ("2026-09-29", "2026-09-30", "2026-10-01", "2026-10-02")]
          == ["matches", "DIFFERENT POINTING", "DIFFERENT POINTING", "UNCONFIRMED"]
          and "where the reference day 2026-09-29 reads +0.91/+0.26" in v["2026-09-30"]["why"]
          and i["horizon_anchor_deg"] == [0.908, 0.255] and i["clear_horizon_days"] == 1
          and lo["measured_different_days"] == ["2026-09-30", "2026-10-01"] and not lo["leave_out"],
          "a 0.6 deg tilt on half the clear days (features unmatched): measured against the reference day, not "
          "missed by a median that falls between the pointings (2 of the 3 days measured differ: the 'most days "
          "differ' guard then keeps them, with its warning)")
    rows = [dict(r) for r in storm3 + E] + [_pt_row("c1", "2026-10-03", "off the calibration", 41.8, 0.91, 0.25, 0.4,
                                                    "ok", 900)]
    sp.pointing_verdicts(rows, ref_day="2026-09-29")
    res = dict(sp.pointing_leave_out(rows), rows=rows, info={})
    check(sp.skip_days_of(res, SimpleNamespace(keep_moved_days=False))
          == ["2026-09-26", "2026-09-27", "2026-09-28", "2026-09-30", "2026-10-01", "2026-10-02"]
          and (sp.pointing_majority_note(res, SimpleNamespace(keep_moved_days=False)) or "").startswith(
              "6 of the window's 8 days are left out"),
          "... with one more matched day: left out, and the warning counts the 6 of 8 days left out")
    # ... with nothing tying the clear days to the reference (no feature reading, the reference day not clear):
    # two equal groups; neither is called measured, and 'horizon only' says the window may hold two pointings
    rows = [dict(r) for r in storm3] + [dict(r, features="unmatched", inliers=3, d_azimuth_deg="", d_tilt_deg="",
                                             d_roll_deg="") for r in E]
    by(rows)["2026-10-02"].update(horizon_droll_deg=0.25)
    sp.pointing_verdicts(rows)
    v = by(rows)
    check(not any(r["verdict"] == "DIFFERENT POINTING" for r in rows)
          and v["2026-09-30"]["verdict"] == "UNCONFIRMED" and "nothing ties either to the reference" in v["2026-09-30"]["why"]
          and v["2026-09-29"]["verdict"] == "horizon only" and "the window may hold two pointings" in v["2026-09-29"]["why"]
          and "agree with the other days" not in v["2026-09-29"]["why"],
          "two equal groups of clear days, nothing tying either to the reference: UNCONFIRMED, never 'measured'; "
          "'horizon only' does not claim the days agree")
    # -- a fog front with a sharp lower edge passes the edge tests on 26-28 Sep (review simulation sB); 27 Sep's
    # 154 matched features measure its real pan but no tilt, where the front reads +2.9 deg: 27 Sep is no clear
    # day and confirms nothing, and 26/28 Sep, reading the same front, are UNCONFIRMED (not 'a lasting change')
    FOG = (40.0, 0.3, 0.9, 60)
    rows = [_pt_row("c1", "2026-09-26", "off the calibration", 186.0, 3.85, 0.2, 1.0, "unmatched", 3, edge=FOG),
            _pt_row("c1", "2026-09-27", "off the calibration", 186.0, 3.86, 0.25, 1.0, "MOVED", 154, (-1.21, 0.0, 0.03),
                    edge=FOG),
            _pt_row("c1", "2026-09-28", "off the calibration", 188.0, 3.9, 0.2, 1.0, "unmatched", 3, edge=FOG)]
    rows += [_pt_row("c1", f"2026-{md}", "off the calibration", 42.0, 0.93, 0.24, 0.4, "ok", 600)
             for md in ("09-29", "09-30", "10-01", "10-02")]
    i = sp.pointing_verdicts(rows, ref_day="2026-09-29")
    v = by(rows)
    check([v[x]["verdict"] for x in ("2026-09-26", "2026-09-27", "2026-09-28")]
          == ["UNCONFIRMED", "DIFFERENT POINTING", "UNCONFIRMED"]
          and "2026-09-27 read the same on a clear-looking edge while matched features there measured no such change"
          in v["2026-09-26"]["why"] and "where the features measured +0.00/+0.03" in v["2026-09-27"]["why"]
          and i["horizon_contradicted_days"] == ["2026-09-27"] and i["clear_horizon_days"] == 4
          and abs(i["window_horizon_offset_px"] - 42.0) < 1e-9,
          "a horizon-like fog front: the day whose features contradict it is no clear day and confirms nothing; the "
          "days reading the same front are UNCONFIRMED; the window's 4 clear days only")
    # -- a hazy winter week (review probe): features unmatched, one clear day, six hazy fitted readings 0.33-0.55
    # deg off it: UNCONFIRMED on a non-clear reading alone, most of the window. Left out (a storm can turn the
    # camera on such days: c1 26-28 Sep 2026), and the warning says how many and why
    rows = [_pt_row("c1", "2025-01-18", "matches", 7.0, 0.15, 0.0, 0.4, "unmatched", 3)]
    rows += [_pt_row("c1", f"2025-01-{19 + k}", "matches", 9.0, t, 0.02, 0.6, "unmatched", 3,
                     edge=(10.0, 0.4, 1.0, 200)) for k, t in enumerate((0.5, 0.62, 0.48, 0.7, 0.55, 0.66))]
    sp.pointing_verdicts(rows)
    res = dict(sp.pointing_leave_out(rows), rows=rows, info={})
    no = SimpleNamespace(keep_moved_days=False)
    note = sp.pointing_majority_note(res, no) or ""
    hazy = [f"2025-01-{d}" for d in range(19, 25)]
    check(sum(r["verdict"] == "UNCONFIRMED" for r in rows) == 6 and res["horizon_alone_days"] == hazy
          and sp.skip_days_of(res, no) == hazy and note.startswith("6 of the window's 7 days are left out")
          and "6 of them (2025-01-19" in note and "not a clear horizon" in note and "--keep-moved-days" in note,
          "a hazy week: 6 of 7 days UNCONFIRMED on a non-clear horizon reading alone: left out, and the warning "
          "counts them and says why")
    _, txt = quiet(sp.print_pointing, "c1", res, no)
    check("6 of the window's 7 days are left out" in txt and "LEFT OUT of the detection: 2025-01-19 (UNCONFIRMED)" in txt,
          "print: the warning and the left-out days")
    _, txt = quiet(sp.print_pointing, "c1", res, SimpleNamespace(keep_moved_days=True))
    check("6 of the window's 7 days" not in txt and "KEPT: 2025-01-19 (UNCONFIRMED)" in txt and "(--keep-moved-days)" in txt,
          "... with --keep-moved-days: kept, no warning about days left out")
    # ... and days left out on clear evidence that are most of the window: left out, with a warning and the count
    rows = [_pt_row("c1", f"2025-01-{18 + k}", "matches", 7.0, 0.15, 0.0, 0.4, "ok", 300) for k in range(3)]
    rows += [_pt_row("c1", f"2025-01-{21 + k}", "off the calibration", 30.0, 0.15 + 0.4 * (k + 1), 0.0, 0.4, "ok", 300)
             for k in range(4)]
    sp.pointing_verdicts(rows)
    res = dict(sp.pointing_leave_out(rows), rows=rows, info={})
    note = sp.pointing_majority_note(res, no) or ""
    check(sp.skip_days_of(res, no) == [f"2025-01-{d}" for d in range(21, 25)]
          and note.startswith("4 of the window's 7 days are left out of the detection") and "--keep-moved-days" in note
          and sp.pointing_majority_note(res, SimpleNamespace(keep_moved_days=True)) is None,
          "4 of 7 days left out (each a clear jump that its zero-turn features contradict): left out, with the count "
          "and the options")
    # -- no saved why says 'left out' (the detection decides that, and --keep-moved-days keeps them)
    allrows, more = _nuc_rows("c1") + _nuc_rows("c2"), [dict(r) for r in storm3 + E]
    sp.pointing_verdicts(allrows[:7])
    sp.pointing_verdicts(allrows[7:])
    sp.pointing_verdicts(more, ref_day="2026-09-29")
    check(not any("left out" in r["why"] for r in allrows + more),
          "no saved why says 'left out' (with --keep-moved-days such a day is kept)")
    # -- the 24 Jan 2025 re-aim (the cameras moved 5.7 m and turned ~15 deg): always left out
    jan = [_pt_row("c1", f"2025-01-{18 + k}", "matches", 7.5, t, r_, 0.4, "unmatched", 4, edge=EDGE_FIXTURE)
           for k, (t, r_) in enumerate(((0.153, 0.012), (0.162, 0.007), (0.164, 0.002), (0.171, 0.0), (0.175, 0.006),
                                        (0.18, -0.003)))]
    for extra, want in (
            ([_pt_row("c1", "2025-01-24", "off the calibration", 90.0, 0.16, -4.08, 1.0, "unmatched", 3,
                      edge=EDGE_FIXTURE)], {"2025-01-24": "UNCONFIRMED"}),
            ([_pt_row("c1", "2025-01-24", "off the calibration", 90.0, 0.16, -4.08, 1.0, "unmatched", 3, edge=EDGE_FIXTURE),
              _pt_row("c1", "2025-01-25", "off the calibration", 91.0, 0.18, -4.05, 1.0, "unmatched", 2, edge=EDGE_FIXTURE)],
             {"2025-01-24": "DIFFERENT POINTING", "2025-01-25": "DIFFERENT POINTING"}),
            ([_pt_row("c1", "2025-01-24", "off the calibration", 90.0, 0.16, -4.08, 1.0, "MOVED", 300, (-15.3, 0.2, -4.0),
                      edge=EDGE_FIXTURE)], {"2025-01-24": "DIFFERENT POINTING"}),
            ([_pt_row("c1", "2025-01-24", "off the calibration", 90.0, 0.16, -4.08, 1.0, "unmatched", 3,
                      edge=(8.0, 0.8, 1.0, 120))], {"2025-01-24": "UNCONFIRMED"})):
        rows = [dict(r) for r in jan] + extra
        sp.pointing_verdicts(rows)
        got = {r["date"]: r["verdict"] for r in rows if r["date"] in want}
        check(got == want and all(r["verdict"] == "horizon only" for r in rows if r["date"] not in want)
              and set(sp.pointing_leave_out(rows)["different_days"]) == set(want),
              f"re-aim: {got} (one day on a clear horizon: unconfirmed; the same change on two days, or the features: "
              f"measured; a faint horizon: unconfirmed): left out either way")
    # -- the Jan fixture: every day unmatched, a clear (rendered) horizon: horizon only, as before, nothing left out
    rows = [dict(r) for r in jan]
    i = sp.pointing_verdicts(rows)
    check(all(r["verdict"] == "horizon only" and r["why"].startswith("tilt and roll within 0.3 deg of the window's "
                                                                       "clear days (+0.17/+0.00 deg)")
              and "features unmatched" in r["why"] and "two pointings" not in r["why"] for r in rows)
          and i["clear_horizon_days"] == 6
          and abs(i["window_horizon_dtilt_deg"] - 0.1675) < 0.0011 and abs(i["window_horizon_droll_deg"] - 0.004) < 1e-9
          and not sp.pointing_leave_out(rows)["different_days"],
          "Jan fixture (all unmatched, rendered horizon prominence ~15-24): 'horizon only' with the same why, nothing "
          "left out, the window +0.17/+0.00 over 6 clear days")
    # -- no clear day at all (the Jan 2025 fixture's c2 render: prominence 6.5-8.0, rival 0.60-0.64): the
    # readings are compared with each other only, no offset against the calibration, nothing claimed measured
    faint = [_pt_row("c2", f"2025-01-{18 + k}", "matches", 9.0 + 0.2 * k, 0.13 + 0.004 * k, -0.02, 0.9, "unmatched", 7,
                     edge=(7.0, 0.62, 1.0, 239)) for k in range(6)]
    i = sp.pointing_verdicts(faint)
    check(all(r["verdict"] == "unchecked" and "not a clear horizon (faint: 7.0x" in r["why"] for r in faint)
          and i["clear_horizon_days"] == 0 and i["window_horizon_offset_px"] is None
          and abs(i["fitted_horizon_offset_px"] - 9.5) < 1e-9 and sp.pointing_offset(i, "x.yaml") is None
          and not sp.pointing_leave_out(faint)["different_days"],
          "no clear horizon on any day, no matched features: 'unchecked' (kept), no offset against the calibration")
    rows = [dict(r) for r in faint]
    rows[2].update(horizon_dtilt_deg=0.6)
    sp.pointing_verdicts(rows)
    check(rows[2]["verdict"] == "UNCONFIRMED" and "(none has a clear horizon)" in rows[2]["why"],
          "... and a jump among such readings: UNCONFIRMED, judged against the window's fitted days")
    # -- rows of an older build (no edge numbers) count as clear, as before Oct 2026, and a one-day horizon jump
    # without matched features is unconfirmed (left out, as it was); the window's offset is now that of the
    # clear days at the reference's pointing (the storm days' readings jump from it: 4 days, not the old 7)
    rows = [{k: v for k, v in r.items() if k not in sp.EDGE_KEYS} for r in _nuc_rows("c1")]
    i = sp.pointing_verdicts(rows)
    check([by(rows)[x]["verdict"] for x in ("2026-09-26", "2026-09-27", "2026-09-28", "2026-09-29")]
          == ["UNCONFIRMED", "DIFFERENT POINTING", "UNCONFIRMED", "matches"] and not i["horizon_not_clear_days"]
          and i["clear_horizon_days"] == 4 and abs(i["window_horizon_dtilt_deg"] - 0.928) < 0.0011,
          "rows without the edge numbers (an older build's): every fitted day clear, 26-28 Sep still left out; the "
          "window +0.93 deg over the 4 days at the reference's pointing")
    # -- an older build's pointing result on disk: read as it is when current; when made by the older check
    # (POINTING_VERSION 5), the detection is stopped with that reason, not 'other photos'
    pd_ = Path(d) / "ptold"
    (pd_ / "pointing").mkdir(parents=True)
    pl = plan_for(pd_)
    a = SimpleNamespace(pointing_frames=3, pointing_fail=sp.POINT_FAIL_DEG, keep_moved_days=False)
    photos = [Path("1737219600.Sat.Jan.18_17_00_00.GMT.2025.CACO03.c1.timex.jpg")]
    (pd_ / "pointing" / "pointing_c1.json").write_text(json.dumps(
        {"rows": [{"date": "2025-01-24", "verdict": "DIFFERENT POINTING", "why": "horizon tilt/roll ..."}],
         "info": {}, "different_days": ["2025-01-24"], "leave_out": True}))
    sig = sp.pointing_sig(pl, "c1", a, photos)
    sp.write_stamp(pd_, "pointing_c1", sig)
    check(sp.pointing_skip_days(pl, {}, "c1", a, photos) == (["2025-01-24"], ""),
          "an older build's pointing JSON (its own keys) still gives the days to leave out")
    sp.write_stamp(pd_, "pointing_c1", dict(sig, v=5))
    try:
        sp.pointing_skip_days(pl, {}, "c1", a, photos)
        msg = ""
    except sp.StepFailed as exc:
        msg = str(exc)
    check("older version of the check (v5" in msg and "--steps pointing" in msg,
          "a pointing result of the older check (v5) is not used by the detection, and the reason says so")
    # -- window dirt: a confident zero turn with a large second consensus at a real turn
    m = {"inliers": 900, "d_azimuth": 0.001, "d_tilt": 0.0, "d_roll": 0.002, "runner_up_inliers": 145,
         "runner_up_deg": 1.0, "runner_up_d": [1.0, 0.0, 0.0]}
    fr = sp.feature_record(m, Path("x.jpg"))
    check(fr["features"] == "ambiguous" and fr["runner_up_d"] == "+1.000/+0.000/+0.000"
          and sp.feature_record(dict(m, runner_up_inliers=21), Path("x.jpg"))["features"] == "ok"
          and sp.feature_record(dict(m, d_azimuth=1.0, runner_up_deg=1.0, runner_up_d=[0.0, 0.0, 0.0]),
                                Path("x.jpg"))["features"] == "MOVED"
          and sp.feature_record({"status": "unmatched", "inliers": 9}, Path("x.jpg"))["features"] == "unmatched"
          and sp.feature_record(None, None)["features"] == "not measured",
          "features: a zero turn with a runner-up of 145 inliers 1 deg away is ambiguous (21: ok; a turned winner: "
          "MOVED, dirt can only add zero)")
    rows = [dict(r) for r in base] + [_pt_row("c2", "2026-10-05", "no horizon", "", "", "", "", "ambiguous", 900,
                                              (0.001, 0.0, 0.002), edge="")]
    by(rows)["2026-10-05"].update(runner_up_inliers=145, runner_up_d="+1.000/+0.000/+0.000")
    sp.pointing_verdicts(rows)
    check(by(rows)["2026-10-05"]["verdict"] == "UNCONFIRMED" and "ambiguous" in by(rows)["2026-10-05"]["why"],
          "an ambiguous match whose runner-up turned 1 deg: UNCONFIRMED (left out), never 'matches'")
    # -- horizon_edge on drawn frames: one crisp sky-over-sea edge, an inverted edge, two equal edges, an edge
    # at the end of the search, a ramp
    rng = np.random.default_rng(1)
    cols, near = np.array([100, 200, 300]), np.array([200.0, 200.0, 200.0])
    from estimate_eo_rotation import observed_horizon

    def edge_of(prof):
        g = np.clip(np.repeat(np.asarray(prof, float)[:, None], 400, axis=1) + rng.normal(0, 1.5, (600, 400)), 0,
                    255).astype(np.uint8)
        obs = observed_horizon(g, cols, near, window=sp.HORIZON_WINDOW_PX)
        return sp.day_edge([sp.horizon_edge(g, cols, near, obs)])
    y = np.arange(600)
    crisp = edge_of(np.where(y < 210, 200, 60))
    inv = edge_of(np.where(y < 210, 60, 200))
    two = edge_of(np.where(y < 150, 200, np.where(y < 300, 60, 200)))
    lim = edge_of(np.where(y < 430, 200, 60))
    ramp = edge_of(np.clip(200 - (y - 50) * 140.0 / 400, 60, 200))
    check(crisp["horizon_edge"] == "clear" and crisp["horizon_prominence"] > 50 and crisp["horizon_rival"] < 0.2
          and crisp["horizon_sky_above"] == 1.0 and 230 <= crisp["horizon_to_limit_px"] <= 245,
          f"horizon_edge: one crisp sky-over-sea edge is clear ({crisp['horizon_prominence']}x, rival "
          f"{crisp['horizon_rival']}, {crisp['horizon_to_limit_px']} px from the search's end)")
    check("inverted" in inv["horizon_edge"] and "rival edge" in two["horizon_edge"]
          and "at search limit" in lim["horizon_edge"] and "faint" in ramp["horizon_edge"],
          f"horizon_edge: dark over bright '{inv['horizon_edge']}', two equal edges '{two['horizon_edge']}', an edge "
          f"{lim['horizon_to_limit_px']} px from the end '{lim['horizon_edge']}', a ramp '{ramp['horizon_edge']}'")


def test_storm_and_wave_caveats(d):
    print("README: pointing whys, the storm tail and high waves from the frames the DEM gridded, one per-day note")
    out = Path(d) / "stormcav"
    for sub in ("waterlines", "forcing", "dem", "pointing"):
        (out / sub).mkdir(parents=True)
    # the NUC's 2026-09-29 build: the pointing results as this step writes them (c1, the new rule) and as an
    # OLDER build wrote them (c2: the old rule's rows, no edge numbers, no overruled list)
    c1 = _nuc_rows("c1")
    i1 = sp.pointing_verdicts(c1)
    (out / "pointing" / "pointing_c1.json").write_text(json.dumps(dict(sp.pointing_leave_out(c1), rows=c1, info=i1)))
    old_c2 = [{"camera": "c2", "date": day, "verdict": "DIFFERENT POINTING" if day < "2026-09-28" else "matches",
               "why": (f"horizon tilt/roll {t[2]:+.2f}/{t[3]:+.2f} deg vs the window's +0.56/-0.64 (> 0.3)"
                       if day < "2026-09-28" else ""), "horizon_offset_px": t[1], "d_azimuth_deg": t[7],
               "d_tilt_deg": t[8], "d_roll_deg": t[9]} for day, t in NUC_TABLE["c2"].items()]
    (out / "pointing" / "pointing_c2.json").write_text(json.dumps(
        {"rows": old_c2, "info": {}, "different_days": ["2026-09-26", "2026-09-27"], "leave_out": True}))
    sp.write_stamp(out, "detect_c1", {"x": 1}, extra={"skip_days": ["2026-09-26", "2026-09-27", "2026-09-28"]})
    sp.write_stamp(out, "detect_c2", {"x": 2}, extra={"skip_days": ["2026-09-26", "2026-09-27"]})
    # the forcing: three high-wave spans in 33 h (27 Sep 11:30-12:30, 27 Sep 20:30, 28 Sep 20:00-22:00)
    t0 = datetime(2026, 9, 26, tzinfo=timezone.utc)
    high = {datetime(2026, 9, 27, 11, 30): 2.4, datetime(2026, 9, 27, 12, 0): 3.71, datetime(2026, 9, 27, 12, 30): 2.2,
            datetime(2026, 9, 27, 20, 30): 3.08, datetime(2026, 9, 28, 20, 0): 2.11, datetime(2026, 9, 28, 21, 0): 2.05,
            datetime(2026, 9, 28, 22, 0): 2.02}
    with open(out / "forcing" / "waves.csv", "w") as f:
        f.write("time_utc,epoch,wvht_m,dpd_s,source\n")
        for k in range(7 * 48):
            t = t0 + timedelta(minutes=30 * k)
            f.write(f"x,{int(t.timestamp())},{high.get(t.replace(tzinfo=None), 1.2)},9,live archive\n")

    def frames(spec, dem=True):
        """{(cam, day, hour): hs} -> the waterline file and (dem) the DEM's record of what it gridded"""
        with open(out / "waterlines" / "contour_points_ground.csv", "w") as f:
            f.write("source_file,camera,capture_time_utc,capture_epoch,offshore_hs_m\n")
            for (cam, day, hh), hs in sorted(spec.items()):
                t = datetime.fromisoformat(day).replace(tzinfo=timezone.utc) + timedelta(hours=hh)
                f.write(f"{photo_name(t, 'CACO05', cam)},{cam},{t:%Y-%m-%dT%H:%M:%SZ},{int(t.timestamp())},{hs}\n")
        info = out / "dem" / "2026-09-29_info.json"
        if dem:
            info.write_text(json.dumps(
                {"window_start": "2026-09-26", "window_end": "2026-10-02", "filters": {
                    "max_hs": 1.5, "rough_frames": sum(h > 1.5 for h in spec.values()), "camera_days_rejected": [],
                    "max_day_offset": 0.15}}))
        elif info.exists():
            info.unlink()
        fpd = {}
        for (cam, day, hh) in spec:
            fpd.setdefault(cam, {})[day] = fpd.setdefault(cam, {}).get(day, 0) + 1
        return {"frames_per_day": fpd}
    after = {(c, f"2026-{dd}", hh): 1.2 for c in ("c1", "c2") for dd in ("09-29", "09-30", "10-01", "10-02")
             for hh in ((14, 15, 16, 17) if c == "c1" else (14, 15, 16))}
    nuc = dict(after)
    nuc.update({("c2", "2026-09-28", 20 + h / 2.0): 2.1 for h in range(4)})     # 4 during the 28 Sep span
    nuc.update({("c1", "2026-09-30", 18): 1.7, ("c1", "2026-09-30", 19): 1.7})  # 2 rough frames after it all
    plan = {"era": "live", "date": "2026-09-29", "first": "2026-09-26", "last": "2026-10-02", "out": out,
            "surveys": [{"name": "rtk", "survey_type": "points", "survey_date": "2026-09-29",
                         "path": str(named(d, "rtk_shots.csv"))}],
            "cams": {c: {"first": "2026-09-26", "last": "2026-10-02", "eo_file": f"CACO05_{c}_20251113_EO-CV.yaml"}
                     for c in ("c1", "c2")}}
    a = cfg_args(keep_moved_days=False)
    forcing = {"waves": {"setup_in_window": {"median_m": 0.3, "daytime_median_m": 0.28}}}
    cav = sp.collect_caveats(plan, a, {}, forcing, frames(nuc))
    text = "\n".join(cav)
    storm = [c for c in cav if c.startswith("STORM:")]
    waves = [c for c in cav if c.startswith("HIGH WAVES")]
    # (pointing) the left-out days with their verdicts and whys; an older build's rows read as they are
    pc1 = next(c for c in cav if c.startswith("c1: days left out of the detection for their pointing:"))
    check("2026-09-26 (UNCONFIRMED: the horizon read +4.23/+0.06" in pc1 and "2026-09-27 (DIFFERENT POINTING: 48 "
          "matched features" in pc1 and "for a different pointing" not in text,
          "README: c1's left-out days with their verdicts and the measured reason (not 'a different pointing')")
    check("c2: days left out of the detection for their pointing: 2026-09-26 (DIFFERENT POINTING: horizon tilt/roll "
          "+2.62/+3.84" in text and "c2: kept, although" not in text,
          "README: an older build's pointing JSON (old rows) is read as it is")
    # (a) the storm's tail: what the DEM gridded, not the waterline days
    check(len(storm) == 1 and "gridded with the rest" not in storm[0]
          and "the DEM grids no frame of the storm's tail before the 2026-09-29 survey (2026-09-26 .. 2026-09-28)"
          in storm[0] and "c1 2026-09-26 (UNCONFIRMED), c1 2026-09-27 (DIFFERENT POINTING), c1 2026-09-28 "
          "(UNCONFIRMED), c2 2026-09-26 (DIFFERENT POINTING), c2 2026-09-27 (DIFFERENT POINTING) left out by the "
          "pointing check" in storm[0] and "c2 2026-09-28 every frame left out: offshore Hs > 1.5 m (the DEM's wave "
          "filter)" in storm[0],
          "(a) STORM: the DEM grids no storm-tail frame; which days the pointing check left out, and that the DEM's "
          "wave filter took c2 28 Sep")
    # (b) + (c) the three spans in one caveat, the DEM's own frames, one per-day note
    check(len(waves) == 1 and "2026-09-27 11:30 .. 2026-09-27 12:30 UTC (peak 3.71 m at 2026-09-27 12:00 UTC), "
          "2026-09-27 20:30 UTC (peak 3.08 m at 2026-09-27 20:30 UTC) and 2026-09-28 20:00 .. 2026-09-28 22:00 UTC"
          in waves[0] and "3 spans less than 48 h apart" in waves[0],
          "(c) three high-wave spans within 33 h: ONE caveat that lists them")
    check("the DEM grids only 28 frame(s) after them (2026-09-29, 2026-09-30, 2026-10-01, 2026-10-02)" in waves[0]
          and "both sides" not in waves[0] and "34 frame(s)" not in waves[0]
          and "it left out 6 frame(s) with offshore Hs > 1.5 m (the DEM's wave filter): 4 frame(s) during them "
          "(2026-09-28); 2 frame(s) after them (2026-09-30)" in waves[0] and "the survey (2026-09-29) is after them"
          in waves[0],
          "(b) HIGH WAVES counts the 28 frames the DEM gridded (all after), not the 34 waterlines, and names the 6 it "
          "left out (4 during, 2 after)")
    check(text.count(sp.PER_DAY_ROWS_NOTE) == 1
          and any(c.startswith("STORM and HIGH WAVES (the caveats above): ") for c in cav)
          and not any(sp.PER_DAY_ROWS_NOTE in c for c in storm + waves),
          "(c) the per-day note is given once for the storm and high-wave caveats (was 4 copies on the NUC)")
    # the new rule's c2 (26/27 Sep kept): their rows reach the DEM unless their waves are too high
    c2 = _nuc_rows("c2")
    i2 = sp.pointing_verdicts(c2)
    (out / "pointing" / "pointing_c2.json").write_text(json.dumps(dict(sp.pointing_leave_out(c2), rows=c2, info=i2,
                                                                       horizon_overruled_days=i2[
                                                                           "horizon_overruled_days"])))
    sp.write_stamp(out, "detect_c2", {"x": 3}, extra={"skip_days": []})
    new = dict(nuc)
    new.update({("c2", "2026-09-26", 15): 1.9, ("c2", "2026-09-27", 15): 1.4, ("c2", "2026-09-27", 16): 1.4})
    cav = sp.collect_caveats(plan, a, {}, forcing, frames(new))
    storm = [c for c in cav if c.startswith("STORM:")]
    waves = [c for c in cav if c.startswith("HIGH WAVES")]
    check(len(storm) == 1 and "the DEM grids frames of the storm's tail before the 2026-09-29 survey with the rest: "
          "c2 2026-09-27 (2 frames); not gridded: c1 2026-09-26 (UNCONFIRMED), c1 2026-09-27 (DIFFERENT "
          "POINTING), c1 2026-09-28 (UNCONFIRMED) left out by the pointing check" in storm[0]
          and "c2 2026-09-26, c2 2026-09-28 every frame left out: offshore Hs > 1.5 m" in storm[0],
          "(a) c2 26/27 Sep kept by the new rule: 27 Sep gridded (said, with its count), 26 and 28 Sep taken by the "
          "DEM's wave filter")
    check(len(waves) == 1 and "the DEM grids the frames of both sides together: 2 frame(s) between them (2026-09-27); "
          "28 frame(s) after them" in waves[0] and "1 frame(s) before them (2026-09-26)" in waves[0],
          "(b) frames between the spans are named; the frame the DEM left out before them too")
    over = next((c for c in cav if c.startswith("c2: kept, although the sea-horizon check alone")), "")
    check("2026-09-26 (matches: the horizon read +2.62/+3.84" in over and "2026-09-27 (matches:" in over,
          "README: the days the features kept against the horizon are listed with their why")
    # without the DEM's record (the DEM step has not run): the waterline days and frames, said to be that
    cav = sp.collect_caveats(plan, a, {}, forcing, frames(new, dem=False))
    storm = [c for c in cav if c.startswith("STORM:")]
    waves = [c for c in cav if c.startswith("HIGH WAVES")]
    check(len(storm) == 1 and "2026-09-26, 2026-09-27, 2026-09-28 (the storm's tail, before the 2026-09-29 survey) are "
          "in the window (no DEM record says which of their frames it gridded)" in storm[0]
          and len(waves) == 1 and "the waterlines hold frames of both sides" in waves[0]
          and "(no DEM record says which of them it gridded)" in waves[0] and "30 frame(s) after them" in waves[0],
          "without the DEM's record: the waterline days and frames, never called gridded")
    # a camera with no clear sea horizon on any day and no matched features (as the Jan 2025 fixture's c2)
    faint = [_pt_row("c2", f"2026-09-{26 + k}", "matches", 9.0, 0.13, -0.02, 0.9, "unmatched", 7,
                     edge=(7.0, 0.62, 1.0, 239)) for k in range(5)]
    i2 = sp.pointing_verdicts(faint)
    (out / "pointing" / "pointing_c2.json").write_text(json.dumps(dict(sp.pointing_leave_out(faint), rows=faint,
                                                                       info=i2)))
    cav = sp.collect_caveats(plan, a, {}, forcing, frames(new))
    check(any(c.startswith("c2: no day of the window showed a clear sea horizon") and "+0.13/-0.02 deg" in c
              and "no offset against the calibration is reported" in c for c in cav)
          and any(c.startswith("c2: nothing measured the pointing on any day") for c in cav)
          and not any(c.startswith("c2: the sea horizon sits") for c in cav),
          "README: no clear horizon on any day: its readings given as such, no offset caveat, 'nothing measured'")
    # --keep-moved-days: the README says the days were KEPT, and no why it quotes says they were left out
    sp.write_stamp(out, "detect_c1", {"x": 4}, extra={"skip_days": []})
    cav = sp.collect_caveats(plan, cfg_args(keep_moved_days=True), {}, forcing, frames(new))
    kept = next((c for c in cav if c.startswith("c1: days with a different or unconfirmed pointing KEPT in the "
                                                 "detection (--keep-moved-days): 2026-09-26 (UNCONFIRMED: ")), "")
    check(kept and "left out" not in kept and "2026-09-27 (DIFFERENT POINTING: 48 matched features" in kept,
          "README with --keep-moved-days: the days KEPT, with whys that do not say 'left out'")


def _rotated(img, io, eo_src, eo_dst, step=4):
    """The exact image of a camera that only turned, from eo_src to eo_dst (each pixel's ray rotated into
    the source frame and projected with the lens model; scipy-free, as the NUC has no scipy)."""
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


def test_feature_runner_up(d):
    print("pointing_check.measure(second=True) on a drawn scene turned 1.0 deg under dirt fixed on the window")
    import cv2
    import pointing_check as pc
    from georectify import load_intrinsics, load_extrinsics
    io = load_intrinsics(HERE / "calibration" / "CACO05_c2_20240801_IO.yaml")
    eo = load_extrinsics(HERE / "calibration" / "CACO05_c2_20251113_EO-CV.yaml")
    nu, nv = int(io[0]), int(io[1])
    rng = np.random.default_rng(1)
    eo2 = eo.copy()
    eo2[3:6] += np.radians([1.0, 0.0, 0.0])
    wd = Path(d) / "runner_up"
    wd.mkdir()

    def scene(contrast):
        t = cv2.GaussianBlur(rng.normal(0, 1, (nv, nu)).astype(np.float32), (0, 0), 3.0)
        t = 128 + contrast * t / t.std()
        return np.clip(np.repeat(t[:, :, None], 3, axis=2), 0, 255).astype(np.uint8)

    def dirt(img, n):
        r, out = np.random.default_rng(7), img.copy()
        for _ in range(n):
            x, y, rad = int(r.uniform(50, nu - 50)), int(r.uniform(50, nv - 50)), int(r.uniform(3, 8))
            cv2.circle(out, (x, y), rad, (20, 20, 20), -1)
            cv2.circle(out, (x, y), rad + 3, (235, 235, 235), 2)
        return out
    got = {}
    for tag, contrast, n in (("scene wins", 40, 600), ("dirt wins", 18, 1100), ("scene lost", 15, 1200)):
        tex = scene(contrast)
        ref_img = wd / f"ref_{n}.jpg"
        cv2.imwrite(str(ref_img), dirt(tex, n))
        ref_img.with_suffix(".json").write_text(json.dumps({"eo": eo.tolist()}))
        m = pc.measure(dirt(_rotated(tex, io, eo, eo2), n), [pc.Reference(ref_img)], io, min_inliers=40, second=True)
        plain = pc.measure(dirt(_rotated(tex, io, eo, eo2), n), [pc.Reference(ref_img)], io, min_inliers=40)
        got[tag] = (m, sp.feature_record(m, Path("x.jpg")), plain)
    m, fr, plain = got["scene wins"]
    check(fr["features"] == "MOVED" and abs(m["d_azimuth"] - 1.0) < 0.02 and m["runner_up_inliers"] >= 10
          and abs(m["runner_up_d"][0]) < 0.02 and abs(m["runner_up_deg"] - 1.0) < 0.02
          and plain["inliers"] == m["inliers"] and "runner_up_inliers" not in plain,
          f"the scene's 1.0 deg pan wins ({m['inliers']} inliers); the runner-up is the dirt at zero turn "
          f"({m['runner_up_inliers']}); second=False gives the same winner, without a runner-up")
    m, fr, _ = got["dirt wins"]
    check(fr["features"] == "ambiguous" and max(abs(m["d_azimuth"]), abs(m["d_tilt"]), abs(m["d_roll"])) < 0.02
          and m["runner_up_inliers"] >= sp.AMBIGUOUS_INLIERS and abs(m["runner_up_d"][0] - 1.0) < 0.02,
          f"heavier dirt wins at zero turn ({m['inliers']} inliers) and the scene's 1.0 deg is the runner-up "
          f"({m['runner_up_inliers']}): 'ambiguous'")
    m, fr, _ = got["scene lost"]
    check(fr["features"] == "ok" and m["runner_up_inliers"] < sp.RUNNER_UP_SHOWN_MIN and fr["runner_up_d"] == "",
          f"the scene too faint to match: the dirt is the only consensus ({m['inliers']} inliers at zero turn, a "
          f"runner-up of {m['runner_up_inliers']}), and the match reads 'ok'")
    rows = [_pt_row("c2", f"2026-10-0{k}", "off the calibration", 33.0, 0.55, -0.67, 0.2, "ok", 1400) for k in range(1, 5)]
    r = _pt_row("c2", "2026-10-05", "off the calibration", 66.0, 1.15, -0.67, 0.2, edge=None)
    r.update({k: v for k, v in fr.items() if k != "feature_frame"})
    sp.pointing_verdicts(rows + [r])
    check(r["verdict"] == "UNCONFIRMED" and "UNCONFIRMED" in sp.LEFT_OUT_VERDICTS,
          "... so against a clear horizon that moved 0.6 deg it decides nothing: UNCONFIRMED, left out")


def test_live_rows(d):
    print("live rows: window, hours, station ID; malformed rows counted")
    d = Path(d) / "live"
    d.mkdir()
    src = d / "contours.csv"
    rows = ["source_file,camera,capture_time_utc,tide_elevation_navd88,offshore_hs_m,offshore_tp_s"]
    for st, h in (("CACO05", 14), ("CACO05", 15), ("CACO04", 16)):
        t = datetime(2026, 9, 29, h, 0, tzinfo=timezone.utc)
        n = photo_name(t, st, "c2").replace(".jpg", "")
        rows += [f"{n},c2,{t:%Y-%m-%dT%H:%M:%SZ},0.5,1.0,9.0"] * 3
    rows.append("x,c2,2026-09-29T1")                  # a line cut off by the cron's rewrite
    src.write_text("\n".join(rows) + "\n")
    plan = {"cams": {"c2": {"first": "2026-09-29", "last": "2026-09-29", "hours": (0.0, 24.0), "station": "CACO05"}}}
    n_in, n_out, frames, other, bad = sp.read_live_rows(src, d / "out.csv", plan)
    check(n_out == 6 and len(frames) == 2, "only the date's station ID kept (CACO05)")
    check(other == {"c2": {"CACO04": 3}}, "rows of frames named for another station counted, per camera")
    check(bad == 1, "a cut-off row skipped and counted, no traceback")
    # with --force the rows are read once per run, not once per step that asks for them (the forcing
    # and the detection): two reads would leave them two different builds, never one product
    pl = dict(plan, out=d / "prod", first="2026-09-29", last="2026-09-29", date="2026-09-29")
    a = SimpleNamespace(live_contours=str(src), force=True)
    state = {}
    quiet(sp.live_rows_file, pl, a, state)
    t1 = sp.token(pl["out"], "archive_rows")
    quiet(sp.live_rows_file, pl, a, state)
    check(t1 and sp.token(pl["out"], "archive_rows") == t1, "--force: the live rows read once per run")
    quiet(sp.live_rows_file, pl, a, {})
    check(sp.token(pl["out"], "archive_rows") != t1, "... and again by the next forced run")


def test_fix1(d):
    print("horizon dip, in-progress mark and summary re-check, storm, camera split, independent counts")
    d = Path(d) / "fix1"
    d.mkdir()
    # the predicted sea horizon holds the full Earth-curvature dip: the true horizon of a PERFECT camera
    # (geometric dip with refraction k, projected with the station EO) fits as no tilt
    from georectify import load_extrinsics, load_intrinsics
    from view_reproject import ground_to_pixel
    from horizon_check import fit_tilt_roll
    import estimate_eo_rotation as er
    for eo_name, cam in (("CACO03_c1_20250123_EO.yaml", "c1"), ("CACO05_c2_20251113_EO-CV.yaml", "c2")):
        io_ = load_intrinsics(HERE / "calibration" / f"CACO05_{cam}_20240801_IO.yaml")
        eo = load_extrinsics(HERE / "calibration" / eo_name)
        cols = np.arange(150, int(io_[0]) - 150, 100)
        fits = []
        for k in (0.0, er.REFRACTION_K, 0.25):
            dip, dd_ = er.horizon_dip(eo[2], k), 20000.0
            az = eo[3] + np.deg2rad(np.linspace(-80, 80, 6000))
            U, V, ok = ground_to_pixel(eo[0] + dd_ * np.sin(az), eo[1] + dd_ * np.cos(az), eo[2] - dd_ * np.tan(dip),
                                       io_, eo)
            o = np.argsort(U[ok])
            obs = np.interp(cols, U[ok][o], V[ok][o], left=np.nan, right=np.nan)
            fits.append(float(fit_tilt_roll(io_, eo, cols, obs)[0][0]))
        check(abs(fits[1]) < 0.003 and max(abs(f) for f in fits) < 0.015 < sp.OFFSET_REPORT_DEG,
              f"{eo_name}: a perfect camera fits tilt {fits[1]:+.3f} deg (k 0.13), {fits[0]:+.3f}/{fits[2]:+.3f} for "
              f"k 0/0.25 (was +0.06-0.07 with half the dip), well under the {sp.OFFSET_REPORT_DEG} deg report limit")
    # a build marks the product in progress first: an interrupted build never leaves the previous
    # README / provenance looking current, and --summary re-checks against the stamps on disk
    out = d / "2025-01-23"
    out.mkdir()
    (out / "README.txt").write_text("SURVEY-DATE PRODUCT 2025-01-23  (built x, status: complete)\nold numbers\n")
    (out / "provenance.json").write_text('{"date": "2025-01-23", "status": "complete", "window": {"first_day": '
                                         '"2025-01-18", "last_day": "2025-01-23"}, "comparisons": [{"name": '
                                         '"jan_lidar", "current": true, "headline": {"n": 5, "median": 0.1}}]}')
    plan = {"date": "2025-01-23", "era": "adcp", "first": "2025-01-18", "last": "2025-01-23", "out": out,
            "cams": {}, "surveys": []}
    a = cfg_args(output_root=str(d), photo_roots=sp.DEFAULT_PHOTO_ROOTS, adcp=None, adcp_navd88=None, chatham=None,
                 wis=None, ndbc=None, live_contours=str(sp.DEFAULT_LIVE_CONTOURS), setup_fitted_to=None,
                 geotiff_epsg=sp.GEOTIFF_EPSG, gnssr_spline=None, gauge_csv=None, live_cron=str(sp.LIVE_CRON),
                 pointing_log_dir=str(HERE / "archive"), pointing_frames=3, pointing_fail=sp.POINT_FAIL_DEG,
                 dem_cell=sp.DEM_CELL, dem_min_points=sp.DEM_MIN_POINTS, dem_max_spread=sp.DEM_MAX_SPREAD,
                 dem_max_hs=sp.DEM_MAX_HS, dem_max_day_offset=sp.DEM_MAX_DAY_OFFSET, keep_moved_days=False,
                 no_download=False)
    sp.mark_in_progress(plan, a)
    sp.mark_in_progress(plan, a)                       # a second attempt: one banner, the old text kept once
    rd = (out / "README.txt").read_text()
    prov = sp.load_json(out / "provenance.json")
    check(rd.startswith("BUILD IN PROGRESS / INTERRUPTED") and rd.count("BUILD IN PROGRESS / INTERRUPTED") == 1
          and rd.count("old numbers") == 1 and prov["status"] == "in progress"
          and prov["previous_status"] == "complete" and not prov["comparisons"][0]["current"],
          "the build marks README and provenance in progress (once), every comparison not current")
    cur, why, h, _ = sp.summary_comparison(prov, out, prov["comparisons"][0])
    check(not cur and "in progress" in why and not h, "--summary shows nothing of a product in progress as current")
    prov2 = dict(prov, status="complete")
    cur, why, h, _ = sp.summary_comparison(prov2, out, prov2["comparisons"][0])
    check(not cur and "no finished comparison" in why,
          "--summary re-checks the stamps on disk: a 'complete' provenance with no comparison stamp is not current")
    # the storm of 25-26 Sep 2026 reaches into the live window: said, with the storm-tail days
    pl = {"date": "2026-09-29", "era": "live", "first": "2026-09-26", "last": "2026-10-02", "out": d,
          "cams": {}, "surveys": [{"name": "rtk", "survey_type": "points", "survey_date": "2026-09-29", "path": ""}]}
    cav = sp.storm_caveats(pl, cfg_args(survey_dirs=[]), {"frames_per_day": {"c1": {"2026-09-26": 3, "2026-09-27": 2,
                                                                                    "2026-09-30": 4}}},
                           "2026-09-26", "2026-10-02")
    check(len(cav) == 1 and cav[0].startswith("STORM: the 25-26 Sep 2026 storm") and "2026-09-26, 2026-09-27" in cav[0]
          and "-0.12 m" in cav[0] and "per-day rows" in cav[0],
          "2026-09-29: the storm caveat names the storm-tail frame days before the survey and the measured change")
    pl_jan = dict(pl, first="2025-01-18", last="2025-01-23")
    check(not sp.storm_caveats(pl_jan, cfg_args(survey_dirs=[]), {}, "2025-01-18", "2025-01-23"),
          "... and nothing for a window far from it")
    # two cameras off in opposite directions: flagged, whatever the pooled median
    h = {"median": 0.0, "by_camera": {"c1": {"n": 300, "median": -0.15}, "c2": {"n": 226, "median": 0.15},
                                      "neither": {"n": 3, "median": 2.0}}}
    check("0.30 m apart" in sp.camera_split(h) and "pooled median hides" in sp.camera_split(h)
          and not sp.camera_split({"by_camera": {"c1": {"n": 9, "median": 0.02}, "c2": {"n": 9, "median": 0.05}}}),
          "cameras 0.30 m apart flagged (the 'outside the views' group ignored); 0.03 m apart: not")
    # a points survey: 33 values from 7 survey points is 7, TOO FEW
    wl = {"n": 33, "median": -0.03, "nmad": 0.09, "rmse": 0.1, "unit": "points", "frames": 32, "survey_points": 7,
          "n_independent": 7, "n_independent_unit": "survey points",
          "by_camera": {"c1": {"n": 18, "median": -0.02, "frames": 18, "survey_points": 3}}}
    t = sp.stat_text(wl, "waterlines - survey")
    check("TOO FEW (7 survey points < 10)" in t and "c1 -0.020 m (n 18, 18 frames, 3 survey points, too few survey points)"
          in sp.per_camera_text(wl["by_camera"], "points"),
          "waterlines on a points survey: judged on the distinct survey points (7: TOO FEW), per camera too")
    t = sp.stat_text({"n": 16, "median": -0.127, "nmad": 0.098, "rmse": 0.19, "unit": "frames", "n_independent": 16,
                      "n_independent_unit": "frames"}, "waterlines - survey")
    check("TOO FEW" not in t, "16 frames on RTK transects: an estimate")


def test_fix2(d):
    print("one build at a time, the window and inputs of the build, map reasons, still water, waves, summary")
    d = Path(d) / "fix2"
    d.mkdir()
    a = cfg_args(output_root=str(d), photo_roots=["/nowhere"], adcp=None, adcp_navd88=None, chatham=None, wis=None,
                 ndbc=None, live_contours="/x/contour_points_timex.csv", setup_fitted_to=None,
                 geotiff_epsg=sp.GEOTIFF_EPSG, gnssr_spline=None, gauge_csv=None, live_cron=str(sp.LIVE_CRON),
                 pointing_log_dir=str(HERE / "archive"), pointing_frames=3, pointing_fail=sp.POINT_FAIL_DEG,
                 dem_cell=sp.DEM_CELL, dem_min_points=sp.DEM_MIN_POINTS, dem_max_spread=sp.DEM_MAX_SPREAD,
                 dem_max_hs=sp.DEM_MAX_HS, dem_max_day_offset=sp.DEM_MAX_DAY_OFFSET, keep_moved_days=False,
                 no_download=False)
    # -- the lock: a second build of the same date is refused while the first runs; a lock left by a
    # build that is no longer running is taken over, and said so. 'The first' is a process of its own whose
    # command line is a survey_products build's (as /proc shows it), whatever runs this test
    out = d / "2026-09-29"
    out.mkdir()
    ok1, _ = sp.take_lock(out, a, "2026-09-29")
    held = sp.load_json(out / sp.LOCK_NAME) or {}
    first = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(300)", "survey_products.py", "--date",
                              "2026-09-29"])
    other = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(300)"])
    try:
        shown = [command_line_shown(proc) for proc in (first, other)]
        check(all(shown), "the two stand-in processes have their command lines in /proc"
                          + ("" if all(shown) else f" (shown: {shown}, after 30 s)"))
        lock = dict(held, pid=first.pid)
        (out / sp.LOCK_NAME).write_text(json.dumps(lock))
        ok2, why2 = sp.take_lock(out, a, "2026-09-29")
        check(ok1 and held.get("pid") == os.getpid() and not ok2 and "RUNNING" in why2,
              "a second build of the date is refused while the first runs" + ("" if ok1 and not ok2 else
                                                                               f" (got: {why2[:120]})"))
        if Path(f"/proc/{other.pid}/cmdline").exists():
            (out / sp.LOCK_NAME).write_text(json.dumps(dict(held, pid=other.pid)))
            ok5, note5 = sp.take_lock(out, a, "2026-09-29")
            check(ok5 and "not running" in note5,
                  "a lock whose pid is now another program (not a survey_products build): taken over, and said so")
        else:
            skip("pid-reused lock check", "no /proc here: a live pid is all that can be told")
    finally:
        for proc in (first, other):
            proc.kill()
            proc.wait()
    (out / sp.LOCK_NAME).write_text(json.dumps(held))      # this run's own lock again
    sp.release_lock(out)
    check(not (out / sp.LOCK_NAME).exists(), "the lock is released at the end")
    dead = subprocess.Popen([sys.executable, "-c", "pass"])
    dead.wait()
    (out / sp.LOCK_NAME).write_text('{"pid": %d, "host": "%s", "since_utc": "x", "command": "y"}'
                                    % (dead.pid, sp.host_name()))
    ok3, note3 = sp.take_lock(out, a, "2026-09-29")
    check(ok3 and "not running" in note3, "a lock left by a build that died is taken over, and said so")
    sp.release_lock(out)
    (out / sp.LOCK_NAME).write_text('{"pid": 1, "host": "another-computer", "since_utc": "x", "command": "y"}')
    ok4, why4 = sp.take_lock(out, a, "2026-09-29")
    check(not ok4 and "another computer" in why4, "a lock held on another computer is not taken over")
    (out / sp.LOCK_NAME).unlink()
    # -- the window: a --steps compare --window run does not apply the window; the plan keeps the built one
    cams = {c: {"first": "2026-09-30", "last": "2026-10-01", "hours": (0.0, 24.0), "hours_text": "0-24",
                "table_first": "2026-09-26", "table_last": "2026-10-02", "table_hours_text": "0-24"}
            for c in ("c1", "c2")}
    plan = {"date": "2026-09-29", "era": "live", "out": out, "first": "2026-09-30", "last": "2026-10-01",
            "cams": cams, "surveys": []}
    sp.write_stamp(out, "archive_rows", {"src": "/station/contour_points_timex.csv",
                                        "cams": {c: ["2026-09-26", "2026-10-02", [0.0, 24.0], "CACO05"]
                                                 for c in ("c1", "c2")}})
    st = {"not_applied": []}
    a2 = cfg_args(**dict(vars(a), window=["2026-09-30", "2026-10-01"]))
    quiet(sp.keep_built_window, plan, a2, ["compare"], st)
    check(plan["first"] == "2026-09-26" and plan["last"] == "2026-10-02" and cams["c1"]["first"] == "2026-09-26"
          and a2.window is None and st["not_applied"] and "NOT applied" in st["not_applied"][0],
          "--steps compare --window: NOT applied; the plan, README and rebuild line keep the built window")
    plan2 = dict(plan, first="2026-09-30", last="2026-10-01",
                 cams={c: dict(v, first="2026-09-30", last="2026-10-01") for c, v in cams.items()})
    st2 = {"not_applied": []}
    quiet(sp.keep_built_window, plan2, a2, list(sp.STEPS), st2)
    check(plan2["first"] == "2026-09-30" and not st2["not_applied"], "... a run of every step applies it")
    # -- the inputs of the build: a --live-contours whose step failed is not the build's
    b = sp.built_args(cfg_args(**dict(vars(a), live_contours="/nonexistent.csv")), plan)
    na = sp.inputs_not_applied(plan, cfg_args(**dict(vars(a), live_contours="/nonexistent.csv")))
    check(b.live_contours == "/station/contour_points_timex.csv" and na and "/nonexistent.csv" in na[0],
          "a failed --live-contours: the rebuild line names the file the outputs came from, and says it")
    # -- why a photo map is missing: the maps step's own record
    (out / "maps").mkdir()
    (out / "maps" / sp.MAP_STATUS).write_text('{"c2": {"drawn": false, "kind": "script failed", "why": '
                                              '"daily_elevation_map.py FAILED (exit 1; see logs/map_c2.log); the photos '
                                              'were found (70 in waterlines/c2/src)"}}')
    check("FAILED (exit 1" in sp.map_missing_reason(plan, {}, "c2")
          and "photo roots" in sp.map_missing_reason(plan, {}, "c1"),
          "a missing photo map says why: the script failed (with its log), or no photo was found")
    # -- the still-water reference: +0.02 - s x (setup - its mean), not 'LOW by the share'
    wv = {"still_water_reference": {"quantified": True, "share": 0.4, "k": 0.0148, "typical_x": 7.0,
                                    "typical_m": 0.1036, "how": "2026-07-09 .. 2026-10-06"},
          "setup_in_window": {"daytime_median_m": 0.28}}
    t = sp.still_water_caveat(0.037, wv, "adcp")
    check("+0.02 m - s x (its setup - the mean setup" in t and "already in the datum chain" in t
          and "+0.02 - (0.40 x 0.28 - 0.10) = +0.01 m" in t and "lack" not in t,
          "still water: the mean part of the GNSS-R share is in the datum chain; this window +0.01 m")
    pc = dict(sp.PRIOR_CHECKS[0], cameras=["c1", "c2"])
    t = sp.prior_check_caveat({"era": "adcp"}, pc, [{"coef": 0.037}],
                              {"waves": {"setup_in_window": {"daytime_median_m": 0.28}}})
    check("comes on top" not in t and "holds the still-water-reference effect" in t,
          "the earlier check's expectation already holds the still-water effect (ADCP-frame lines)")
    # -- high waves in the window from the date's own forcing
    fdir = out / "forcing"
    fdir.mkdir()
    with open(fdir / "waves.csv", "w") as f:
        f.write("time_utc,epoch,wvht_m,dpd_s,source\n")
        for h in range(0, 24 * 6):
            ep = 1737158400 + h * 3600                 # 2025-01-18 00:00 UTC on
            hs = 2.8 if 50 <= h <= 60 else 1.0          # 20 Jan 02:00-12:00
            f.write(f"x,{ep},{hs},9,adcp\n")
    pj = dict(plan, era="adcp", first="2025-01-18", last="2025-01-23", date="2025-01-23",
              surveys=[{"name": "jan_lidar", "survey_type": "dsm", "survey_date": "2025-01-23"}])
    ev = sp.wave_events(fdir, "2025-01-18", "2025-01-23")
    cav = sp.wave_event_caveats(pj, {"frames_per_day": {"c1": {"2025-01-18": 9, "2025-01-19": 9, "2025-01-21": 9,
                                                               "2025-01-23": 9}}}, "2025-01-18", "2025-01-23")
    check(len(ev) == 1 and abs(ev[0]["peak"] - 2.8) < 1e-9 and len(cav) == 1 and "HIGH WAVES" in cav[0]
          and "2025-01-18, 2025-01-19 before it" in cav[0] and "2025-01-21, 2025-01-23 after it" in cav[0]
          and "the survey (2025-01-23) is after it" in cav[0] and "per-day rows" in cav[0],
          "a high-wave event in the window: the frame days before and after it, the survey after it")
    check(not sp.wave_event_caveats(dict(pj, first="2025-01-22", last="2025-01-23"), {}, "2025-01-22", "2025-01-23"),
          "... and no caveat for a window without one")
    # -- the summary: the reason of a C = 0 INDEPENDENT row comes from the build, wrapped at phrases
    prov = {"setup": {"coef": 0.0, "fits": [{"coef": 0.0, "kind": "none"}]}}
    cc = {"label": "INDEPENDENT", "label_reasons": [], "why": "RTK check shots of 29 Sep 2026 (the user first called "
          "the date 27 Sep): not used by the calibration (GCPs of 13 Nov 2025), the search envelope or the water "
          "level (GNSS-R). INDEPENDENT only for a build whose ..."}
    r = sp.summary_reason(cc, prov)
    check(r.startswith("C = 0: no setup, nothing fitted to this survey; RTK check shots of 29 Sep 2026: not used by")
          and "INDEPENDENT only" not in r and sp.setup_short(prov) == "C = 0: no wave setup",
          "summary: a C = 0 INDEPENDENT row says why from the build (C = 0), then the table's first sentence")
    w = sp.wrap_phrases(r, 60, 2)
    check(len(w) == 2 and w[-1].endswith("...") and all(len(x) <= 60 for x in w)
          and " ".join(w).replace(" ...", "") in r.replace("\n", " "),
          "the figure's reason: two lines, cut at a phrase boundary, never mid-word")
    # -- paired C = 0 numbers on few frames are marked
    wl = {"n": 16, "median": -0.127, "nmad": 0.098, "rmse": 0.19, "unit": "frames", "n_independent": 16,
          "n_independent_unit": "frames", "method": "frames on RTK transects",
          "without_setup": {"n": 12, "median": -0.33, "how": "re-projected",
                            "paired": {"n": 12, "unit": "frames", "n_independent": 12, "n_independent_unit": "frames",
                                       "without_setup_median": -0.33, "with_setup_median": -0.108,
                                       "setup_effect_median": 0.175,
                                       "by_camera": {"c1": {"n": 7, "n_independent": 7, "n_independent_unit": "frames",
                                                            "without_setup_median": -0.389, "with_setup_median": -0.113},
                                                     "c2": {"n": 5, "n_independent": 5, "n_independent_unit": "frames",
                                                            "without_setup_median": -0.272,
                                                            "with_setup_median": -0.102}}}}}
    txt = "\n".join(sp.waterline_readme_lines(wl))
    check("(n 7, too few frames)" in txt and "(n 5, too few frames)" in txt and "TOO FEW" not in txt.split("by camera c1")[0],
          "the C = 0 sensitivity marks each camera's few frames (12 overall: an estimate)")
    # -- synthetic inputs: a marker file in a photo root
    root = d / "roots"
    root.mkdir()
    (root / sp.SYNTHETIC_MARKER).write_text("rendered test photos\n")
    check(sp.run_synthetic(cfg_args(photo_roots=[str(root)], synthetic=None)) == "rendered test photos"
          and sp.run_synthetic(cfg_args(photo_roots=["/nowhere"], synthetic=None)) is None,
          "a SYNTHETIC_FIXTURE file in a photo root marks the build synthetic")


def test_fix3(d):
    print("exit code by severity, DEM day test, not-yet-built label, kept live rows, pointing, cameras, shares")
    d = Path(d) / "fix3"
    d.mkdir()
    # -- --all: the worst date by severity, not the largest number
    check(sp.worst_exit(["built", "failed", "partial"]) == 1 and sp.worst_exit(["partial", "busy"]) == 4
          and sp.worst_exit(["built", "partial", "disabled"]) == 3 and sp.worst_exit(["built", "disabled"]) == 0
          and sp.worst_exit([]) == 0 and sp.worst_exit(["built", "odd"]) == 1,
          "--all: failed (1) > busy (4) > partial (3) > built (0); disabled dates do not count")
    # -- the DEM's day test on two camera-days rejects both: not run, and said why
    wf = d / "w.csv"
    head = "source_file,camera,capture_time_utc,capture_epoch,easting_utm19,northing_utm19,offshore_hs_m\n"
    rows = []
    for day, h, hs in (("2025-03-06", 15, 1.0), ("2025-03-07", 15, 1.0), ("2025-03-08", 15, 2.5)):
        t = datetime.fromisoformat(f"{day}T{h:02d}:00:00+00:00")
        rows.append(f"f{day},c2,{day}T{h:02d}:00:00Z,{int(t.timestamp())},420000,4638000,{hs}")
    wf.write_text(head + "\n".join(rows) + "\n")
    plan = {"first": "2025-03-06", "last": "2025-03-08", "out": d}
    a = SimpleNamespace(dem_max_hs=1.5, dem_max_day_offset=0.15)
    check(sp.camera_days(plan, a, wf) == ["c2 2025-03-06", "c2 2025-03-07"],
          "camera-days as the DEM's day test sees them (the Hs 2.5 m day left out by its wave filter)")
    why = sp.day_test_off(plan, a, wf)
    check("only 2 camera-day(s)" in why and "cannot tell which day is off" in why and "both be rejected" not in why
          and "keeps them all" in why,
          "fewer than 3 camera-days: the day test is not run, and why (dem_from_contours.py would keep them all)")
    check(sp.day_test_off(plan, SimpleNamespace(dem_max_hs=0, dem_max_day_offset=0.15), wf) == ""
          and sp.day_test_off(plan, SimpleNamespace(dem_max_hs=1.5, dem_max_day_offset=0), wf) == "",
          "3 camera-days (no wave filter), or the test off by option: nothing to say")
    cmd = sp.dem_cmds(dict(plan, date="2025-03-06", cams={"c2": {"station": "CACO04", "eo_file": "x_20250219_EO.yaml"}}),
                      SimpleNamespace(dem_cell=2.0, dem_min_points=3, dem_max_spread=0.5, dem_max_hs=1.5,
                                      dem_max_day_offset=0.15, geotiff_epsg=32619), day_test=False)[1]
    check("--max-day-offset" not in [str(c) for c in cmd], "... and dem_from_contours.py is run without it")
    log = d / "dem.log"
    log.write_text("Day consistency (median offset from the DEM without that day, limit +/-0.15 m):\n"
                   "   c2 2025-03-06   +0.190 m   (514 samples)  <-- REJECTED\n"
                   "   c2 2025-03-07   -0.190 m   (408 samples)  <-- REJECTED\n"
                   "   left out 2 camera-day(s), 38167 point(s)\n\nTraceback (most recent call last):\n"
                   "ValueError: zero-size array to reduction operation minimum which has no identity\n")
    check("rejected every camera-day (2 of 2)" in sp.dem_failure_reason(log),
          "a DEM that failed because every camera-day was rejected says so, not a traceback")
    # -- a not-yet-built date: its label with THIS run's C, as --dry-run gives it (fixture tables)
    cfg, srv, tabs = fixture_tables(d)
    rtk = [s_ for s_ in srv if s_["date"] == "2026-09-29"][0]
    lab, why = sp.label_when_built("2026-09-29", rtk, cfg, srv, cfg_args(survey_dirs=[str(d)], **tabs))
    check(rtk["label"] == "INDEPENDENT" and lab == "CIRCULAR" and "fitted to this survey" in why,
          "2026-09-29 not built yet, C = 0.037: CIRCULAR when built (the table says INDEPENDENT)")
    lab0, _ = sp.label_when_built("2026-09-29", rtk, cfg, srv, cfg_args(survey_dirs=[str(d)], setup_coef=0.0, **tabs))
    check(lab0 == "INDEPENDENT", "... with C = 0: INDEPENDENT")
    # -- live era: a window with no rows keeps the rows (and stamp) of the earlier build
    src = d / "contours.csv"
    lines_ = ["source_file,camera,capture_time_utc,tide_elevation_navd88,offshore_hs_m,offshore_tp_s"]
    for h in (14, 15):
        t = datetime(2026, 9, 29, h, 0, tzinfo=timezone.utc)
        n = photo_name(t, "CACO05", "c2").replace(".jpg", "")
        lines_ += [f"{n},c2,{t:%Y-%m-%dT%H:%M:%SZ},0.5,1.0,9.0"] * 2
    src.write_text("\n".join(lines_) + "\n")
    cams = {"c2": {"first": "2026-09-29", "last": "2026-09-29", "hours": (0.0, 24.0), "station": "CACO05"}}
    pl = {"cams": cams, "out": d / "live", "first": "2026-09-29", "last": "2026-09-29", "date": "2026-09-29",
          "era": "live"}
    a = SimpleNamespace(live_contours=str(src), force=False)
    quiet(sp.live_rows_file, pl, a, {})
    t1 = sp.token(pl["out"], "archive_rows")
    before = (pl["out"] / "waterlines" / "archive_rows.csv").read_text()
    cams2 = {"c2": dict(cams["c2"], first="2026-11-01", last="2026-11-03")}
    pl2 = dict(pl, cams=cams2, first="2026-11-01", last="2026-11-03")
    try:
        quiet(sp.live_rows_file, pl2, a, {})
        failed = False
    except sp.StepFailed as exc:
        failed = "earlier build" in str(exc) and "kept" in str(exc)
    check(failed and (pl["out"] / "waterlines" / "archive_rows.csv").read_text() == before
          and sp.token(pl["out"], "archive_rows") == t1
          and not (pl["out"] / "waterlines" / "archive_rows.csv.new").exists(),
          "a live window with no rows fails WITHOUT wiping the earlier rows or their stamp")
    bs = sp.built_settings(pl2)
    check(bs["cams"]["c2"]["first"] == "2026-09-29", "... so the window on disk is still the one that built them")
    # -- per-day rows: a step may be change OR a method error; frames placed against an event by TIME
    check("change of the beach OR a method error" in sp.PER_DAY_ROWS_NOTE and "static synthetic beach"
          in sp.PER_DAY_ROWS_NOTE and "not an error" not in sp.PER_DAY_ROWS_NOTE,
          "per-day rows: a trend may be change OR a method error (the static fixture drifts too)")
    import survey_compare as sc_
    cont_ = {"cam": np.array(["c2"] * 4), "day": np.array(["2025-03-05", "2025-03-05", "2025-03-08", "2025-03-08"]),
             "frame": np.array(["a", "a", "b", "b"]), "Z": np.array([0.3, 0.2, 0.7, 0.5]),
             "Zt": np.array([0.1, 0.0, 0.1, -0.1])}
    wr = sc_.waterline_rows(cont_, np.array([0.02, 0.04, -0.10, -0.12]), zref=np.array([-0.1, -0.05, -0.1, -0.05]))
    dd = {r["group"]: r for r in wr if r["group_type"] == "waterline_day"}
    check(abs(dd["2025-03-05 c2"]["survey_z_median"] + 0.075) < 1e-9 and abs(dd["2025-03-05 c2"]["setup_median"]
                                                                         - 0.2) < 1e-9
          and abs(dd["2025-03-08 c2"]["setup_median"] - 0.6) < 1e-9,
          "per-day rows carry the survey elevation and the setup they were taken at (to compare days like for like)")
    pe = d / "ev"
    (pe / "forcing").mkdir(parents=True)
    (pe / "waterlines").mkdir()
    ep0 = int(datetime(2025, 1, 20, 3, 0, tzinfo=timezone.utc).timestamp())
    with open(pe / "forcing" / "waves.csv", "w") as f:
        f.write("time_utc,epoch,wvht_m,dpd_s,source\n")
        for h in range(-24, 48):
            f.write(f"x,{ep0 + h * 3600},{2.6 if h == 0 else 1.0},9,adcp\n")
    with open(pe / "waterlines" / "contour_points_ground.csv", "w") as f:
        f.write("source_file,camera,capture_time_utc,capture_epoch,offshore_hs_m\n")
        for k, hh in enumerate((14, 15, 16)):          # same day as the 03:00 event, all hours after it
            t = datetime(2025, 1, 20, hh, 0, tzinfo=timezone.utc)
            f.write(f"fr{k},c1,{t:%Y-%m-%dT%H:%M:%SZ},{int(t.timestamp())},1.0\n")
    pj = {"out": pe, "era": "adcp", "first": "2025-01-19", "last": "2025-01-21", "date": "2025-01-20",
          "surveys": [{"name": "s", "survey_type": "dsm", "survey_date": "2025-01-20"}], "cams": {"c1": {}}}
    cv = sp.wave_event_caveats(pj, {"frames_per_day": {"c1": {"2025-01-20": 3}}}, "2025-01-19", "2025-01-21")
    check(len(cv) == 1 and "3 frame(s) after it" in cv[0] and "during it" not in cv[0]
          and "on the same day as the event (its time is not known)" in cv[0],
          "a one-hour event: frames taken hours after it are 'after it', not 'during'; the survey 'the same day'")
    # -- cameras: no flag on too few values; the cause named only for the same times
    few_h = {"waterlines": {"unit": "frames", "by_camera": {
        "c1": {"n": 7, "median": -0.389, "frames": 7}, "c2": {"n": 5, "median": -0.272, "frames": 5}}}}
    check(sp.camera_split(few_h) == "" and "too few to compare the cameras" in sp.camera_split_note(few_h),
          "7 and 5 frames, 0.12 m apart: not flagged, 'too few to compare the cameras'")
    many = {"waterlines": {"unit": "frames", "by_camera": {
        "c1": {"n": 20, "median": -0.30, "frames": 20, "times": 20, "shared_times": 4, "setup_median": 0.40},
        "c2": {"n": 15, "median": -0.15, "frames": 15, "times": 15, "shared_times": 4, "setup_median": 0.24}}}}
    t_ = sp.camera_split(many)
    check("different times" in t_ and "need not be the cameras" in t_ and "c1 0.40 m" in t_
          and "the water level and setup are the same" not in t_,
          "different frames: the split is not put on the cameras; each camera's median setup is given")
    same = {"waterlines": {"unit": "frames", "by_camera": {
        "c1": {"n": 20, "median": -0.30, "frames": 20, "times": 20, "shared_times": 19},
        "c2": {"n": 20, "median": -0.15, "frames": 20, "times": 20, "shared_times": 19}}}}
    check("same times" in sp.camera_split(same) and "pointing or lens model" in sp.camera_split(same),
          "the same times: pointing or lens model")
    check(sp.camera_split({"by_camera": {"c1": {"n": 1, "median": -0.36}, "c2": {"n": 2, "median": -0.20}}}) == "",
          "DEM by camera on 1 vs 2 values: no flag")
    # -- the pointing argument rests on what fitted nothing to the RTK
    cp = d / "pt"
    (cp / "compare").mkdir(parents=True)
    with open(cp / "compare" / "2026-09-29_rtk_points.csv", "w") as f:
        f.write("name,code,description,easting,northing,survey_z_navd88,camera,distance_from_camera_m\n")
        for k, (z, dist) in enumerate(((1.18, 75.3), (1.26, 56.7), (1.68, 50.7), (1.17, 55.1), (3.0, 20.0))):
            f.write(f"{k},P,Transect,0,0,{z},c1,{dist}\n")
    prov_ = {"comparisons": [{"name": "rtk", "current": True, "label": "CIRCULAR", "headline": {"waterlines": {
        "n": 16, "unit": "frames", "setup": "setup-corrected", "by_camera": {"c1": {"n": 10}},
        "line_elevation_compared": {"p5": 1.01, "p95": 1.56, "max": 1.56},
        "without_setup": {"paired": {"by_camera": {"c1": {"n": 7, "frames": 7, "without_setup_median": -0.389,
                                                          "setup_effect_median": 0.2295}}}}}}}]}
    sens = [{"range_m": "0-100", "cells": 300, "shift_m": 4.7, "dem_error_m": -0.41}]
    txt = sp.pointing_vs_survey({"date": "2026-09-29", "out": cp}, "c1", (0.90, 0.23), sens, prov_,
                                [{"coef": 0.037, "kind": "known"}])
    check(txt and "absorbed into C" in txt and "repeat crossings" in txt and "-0.21 to -0.12 m" in txt
          and "raises them by 0.18-0.27 m" in txt and "give C = 0.029-0.043" in txt and "agree best at" not in txt
          and "outside this repository's code" in txt and "+0.006 to +0.034" in txt and "0.024-0.026" not in txt
          and "30-51%" in txt and "cannot be ruled out" in txt and "which the RTK does not show" not in txt
          and "the 0-100 m band" in txt and "too few" in txt,
          "pointing vs lens: C = 0 lines minus the setup of the repeat-crossing C (0.029-0.043, both estimates, "
          "said where they come from) leave -0.21..-0.12 m at ~55 m, 30-51% of the full offset's -0.41 m: a part "
          "cannot be ruled out (and the C fitted to the RTK cannot say)")
    check("OVERSTATED" not in txt, "... no survey-floor note when the build records none")
    wl_ = prov_["comparisons"][0]["headline"]["waterlines"]
    wl_["without_setup"]["above_survey_floor"] = {"frames_below_floor": 12, "frames_all": 12, "frames": 0}
    txt2 = sp.pointing_vs_survey({"date": "2026-09-29", "out": cp}, "c1", (0.90, 0.23), sens, prov_,
                                 [{"coef": 0.037, "kind": "known"}])
    check(txt2 and "12 of 12 of those C = 0 frames lie below their transect's lowest shot" in txt2
          and "OVERSTATED" in txt2,
          "... C = 0 lines below where the survey stops: what is left (and the share of the offset) is OVERSTATED")
    del wl_["without_setup"]["above_survey_floor"]
    dsm = {"comparisons": [dict(prov_["comparisons"][0], headline={"waterlines": dict(
        prov_["comparisons"][0]["headline"]["waterlines"], unit="points")})]}
    check(sp.pointing_vs_survey({"date": "2026-09-29", "out": cp}, "c1", (0.90, 0.23), sens, dsm,
                                [{"coef": 0.037, "kind": "known"}]) is None,
          "... not on a DSM comparison, which pools every range from the camera")
    # -- the summary figure: a CIRCULAR row grey, its C = 0 values as their own mark, synthetic rows tagged
    check(sp.label_without_setup("INDEPENDENT", ["CIRCULAR: the setup coefficient C = 0.037 was fitted to this "
                                                 "survey"]) == "INDEPENDENT"
          and sp.label_without_setup("INDEPENDENT", ["PARTLY-CIRCULAR: the c1 search envelope was placed with this "
                                                     "survey", "CIRCULAR: the setup coefficient C = 0.037 ..."])
          == "PARTLY-CIRCULAR", "the C = 0 values' label: every reason but the setup coefficient's")
    base = {"n": 2, "median_m": -0.19, "nmad_m": None, "waterlines_median_m": -0.127, "waterlines_nmad_m": 0.098,
            "dem_by_camera": {}, "waterlines_by_camera": {"c1": -0.16, "c2": -0.12}, "waterlines_too_few": False,
            "survey_compare_suggests": "", "label_reason": "x", "date": "2026-09-29", "survey": "rtk",
            "c0_label": "INDEPENDENT", "c0_without_setup_median_m": -0.33, "c0_paired_n": 12,
            "c0_paired_unit": "frames", "c0_too_few": False}
    rows_ = [dict(base, label="CIRCULAR"), dict(base, label="INDEPENDENT", date="2025-01-23", survey="jan",
                                                 synthetic="fixture", c0_label=None)]
    sp.summary_figure(rows_, d / "summary.png")
    check((d / "summary.png").stat().st_size > 10000, "the summary figure draws CIRCULAR, C = 0 and SYNTHETIC rows")
    # -- C = 0 texts: LOW by a share of the setup, not by the setup
    cv0 = sp.setup_caveats({"surveys": []}, None, {}, {}, [{"coef": 0.0, "kind": "none"}], "live")
    check("0.6-0.8 x the setup" in cv0[0] and "by about the setup" not in cv0[0],
          "C = 0: lines read LOW by roughly 0.6-0.8 x the setup, not 'by about the setup'")


def test_dem_two_days(d):
    print("dem_from_contours.py: two camera-days that disagree are both kept, never a crash")
    d = Path(d) / "dem2"
    d.mkdir()
    g = d / "ground.csv"
    lines_ = ["source_file,camera,capture_time_utc,capture_epoch,easting_utm19,northing_utm19,"
              "tide_elevation_navd88,beach_elevation_navd88"]
    for day, dz in (("2025-03-06", 0.0), ("2025-03-07", 0.3)):        # the second day reads 0.3 m higher
        for h in (14, 15, 16, 17):
            t = datetime.fromisoformat(f"{day}T{h:02d}:00:00+00:00")
            f = f"f{day}_{h}"
            for i in range(40):
                e, n = 420000.0 + 2.0 * (i % 8) + 1.0, 4638000.0 + 2.0 * (i // 8) + 1.0
                z = 1.0 - 0.05 * (e - 420000.0) + dz + 0.01 * (h - 15)
                lines_.append(f"{f},c2,{t:%Y-%m-%dT%H:%M:%S}+00:00,{int(t.timestamp())},{e:.2f},{n:.2f},{z:.3f},{z:.3f}")
    g.write_text("\n".join(lines_) + "\n")
    r = subprocess.run([sys.executable, str(HERE / "dem_from_contours.py"), str(g), str(d / "dem"),
                        "--max-day-offset", "0.15", "--no-plot"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True, cwd=str(d))
    out = r.stdout
    check(r.returncode == 0 and (d / "dem_dem.asc").exists() and "REJECTED" not in out
          and "kept: fewer than 3 testable days" in out and "every day kept" in out,
          f"two camera-days 0.3 m apart with --max-day-offset 0.15: both kept, DEM built (exit {r.returncode})")


def test_setup_share(d):
    print("a C = 0 build carries no setup: the DEM's share, its REMINDER and the maps say so")
    import daily_elevation_map as dm
    import dem_from_contours as dc
    d = Path(d) / "share"
    d.mkdir()
    hdr = ("source_file,camera,capture_time_utc,capture_epoch,pixel_column,pixel_row,easting_utm19,northing_utm19,"
           "tide_elevation_navd88,setup_correction_m,beach_elevation_navd88")
    res = {}
    for name, setup in (("c0", "0.0"), ("c037", "0.21"), ("mixed", None)):
        lines_ = [hdr]
        for h in range(14, 18):
            t = datetime(2026, 9, 29, h, tzinfo=timezone.utc)
            s_ = setup if setup is not None else ("" if h == 14 else "0.21")   # one frame with no wave record
            for i in range(40):
                e, n = 420000.0 + 2.0 * (i % 8) + 1.0, 4638000.0 + 2.0 * (i // 8) + 1.0
                z = 1.0 - 0.05 * (e - 420000.0) + 0.01 * (h - 15)
                lines_.append(f"f{h},c2,{t:%Y-%m-%dT%H:%M:%S}+00:00,{int(t.timestamp())},{100 + i},{500 + h},"
                              f"{e:.2f},{n:.2f},{z:.3f},{s_},{z + float(s_ or 0):.3f}")
        g = d / f"{name}.csv"
        g.write_text("\n".join(lines_) + "\n")
        st = {}
        quiet(dc.load_points, str(g), stats=st)
        r = subprocess.run([sys.executable, str(HERE / "dem_from_contours.py"), str(g), str(d / name), "--no-plot"],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
        res[name] = (st["setup_share"], r.stdout)
    check(res["c0"][0] == 0.0 and "no setup in" in res["c0"][1] and "include the wave setup" not in res["c0"][1],
          "C = 0 (setup_correction_m 0.0 on every row): setup share 0, the REMINDER says the DEM has no setup")
    check(res["c037"][0] == 1.0 and "include the wave setup" in res["c037"][1] and "reads LOW" not in res["c037"][1],
          "C > 0 on every row: share 1, the REMINDER says the elevations include the setup (not 'reads LOW')")
    check(abs(res["mixed"][0] - 0.75) < 1e-9 and "on 75% of the points (the rest: no wave record" in res["mixed"][1],
          "one frame of four with no wave record: share 0.75, said so")
    check(not dm.carries_setup("0.0") and not dm.carries_setup("") and not dm.carries_setup(None)
          and dm.carries_setup("0.2134"), "elevation maps: a 0.0 (C = 0) or blank setup is water level only")
    import waterline_consistency as wc
    su = {k: wc.read_points(str(d / f"{k}.csv"))["setup_used"] for k in ("c0", "c037", "mixed")}
    check(not su["c0"].any() and su["c037"].all() and su["mixed"].sum() == 3 and su["mixed"].size == 4,
          "consistency diagnostic: its 'wave setup' label follows the lines' setup (none at C = 0, 3 of 4 mixed)")


def test_fit_setup_flat(d):
    print("dem_from_contours.py --fit-setup offers a C only when it beats C = 0 and the C in use")
    d = Path(d) / "fitflat"
    d.mkdir()
    out = {}
    wide, same = (5.0, 7.0, 9.0, 11.0, 13.0, 15.0), (9.9, 9.95, 10.0, 10.05, 10.1, 10.0)
    # (name, sqrt(Hs L0) per frame, planted C, --setup-coef-in-use, setup_correction_m carried with this C)
    cases = (("sharp", wide, 0.04, ["--setup-coef-in-use", "0.02"], None),
             ("in_use", wide, 0.037, ["--setup-coef-in-use", "0.037"], None),
             ("zero", wide, 0.0, ["--setup-coef-in-use", "0.037"], None),
             ("flat", same, 0.04, ["--setup-coef-in-use", "0.08"], None),
             ("carried", wide, 0.04, [], 0.06))
    for name, phis, c_true, extra, carried in cases:
        lines_ = ["source_file,camera,capture_time_utc,capture_epoch,easting_utm19,northing_utm19,"
                  "tide_elevation_navd88,offshore_hs_m,offshore_tp_s" + (",setup_correction_m" if carried else "")]
        for k, phi in enumerate(phis):          # every frame's line at beach 1.0 m = tide + C sqrt(Hs L0)
            t = datetime(2026, 9, 29, 12 + k, tzinfo=timezone.utc)
            hs = (phi / 8.0) ** 2 / (9.81 / (2 * np.pi))
            for i in range(40):
                e, n = 420000.0 + 2.0 * (i % 8) + 1.0, 4638000.0 + 2.0 * (i // 8) + 1.0
                lines_.append(f"f{k},c2,{t:%Y-%m-%dT%H:%M:%S}+00:00,{int(t.timestamp())},{e:.2f},{n:.2f},"
                              f"{1.0 - c_true * phi:.4f},{hs:.4f},8.0" + (f",{carried * phi:.4f}" if carried else ""))
        g = d / f"{name}.csv"
        g.write_text("\n".join(lines_) + "\n")
        r = subprocess.run([sys.executable, str(HERE / "dem_from_contours.py"), str(g), str(d / name), "--fit-setup",
                            "--no-plot"] + extra, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           universal_newlines=True, cwd=str(d))
        out[name] = r.stdout
    check("FLAT" not in out["sharp"] and "--setup-coef 0.04" in out["sharp"]
          and "Beats C = 0 and the C in use (0.02)" in out["sharp"],
          "waves that differ, planted C = 0.04, C in use 0.02: a sharp minimum at 0.04 that beats both, 'To apply' "
          "given")
    check("Not distinguishable from the C in use (0.037)" in out["in_use"] and "To apply" not in out["in_use"],
          "planted C = the C in use: within the near-best band, no other C offered")
    check("Not distinguishable from C = 0" in out["zero"] and "To apply" not in out["zero"],
          "planted C = 0: C = 0 within the near-best band, no C offered")
    check("a FLAT minimum: this fit alone does not pin C" in out["flat"] and "To apply" not in out["flat"],
          "nearly the same waves in every frame: the fit says FLAT and does not offer a C to apply")
    check("these contours carry a setup (C ~0.06)" in out["carried"]
          and "C in use 0.06 (the contours' setup_correction_m)" in out["carried"]
          and "--setup-coef 0.04" in out["carried"],
          "contours that carry a setup: said so, and the C they carry is the C in use compared with")


def test_c_change(d):
    print("another C in the same folder is refused without --replace, and the docs give each C its own root")
    d = Path(d) / "cchange"
    out = d / "root" / "2025-01-23"
    out.mkdir(parents=True)
    plan = {"date": "2025-01-23", "era": "adcp", "out": out, "cams": {"c1": {}, "c2": {}}}
    a0 = SimpleNamespace(setup_coef=0.0, setup_fitted_to=None, output_root=str(d / "root"))
    check(sp.setup_replace_note(plan, a0, sp.STEPS) is None, "nothing built yet: nothing to replace")
    sp.write_stamp(out, "merge", {"x": 1}, extra={"setup": {
        c: {"coef": sp.SETUP_COEF, "fit": sp.fit_sig(sp.setup_fit_info(sp.SETUP_COEF))} for c in ("c1", "c2")}})
    txt = sp.setup_replace_note(plan, a0, sp.STEPS)
    check(txt and txt.startswith(f"the C = {sp.SETUP_COEF:g} product in {out} would be REPLACED")
          and f"--output-root {d / 'root_c0'} for C = 0" in txt and "detection runs again" in txt
          and "--replace" in txt,
          "C = 0 over a C = 0.037 build: 'would be REPLACED', its own root named, detection reruns, --replace to force")
    # the root a C belongs in: the base root for the station's C, <base>_c<C> otherwise -- also from a _c root
    r = Path("/mnt/I2Rgus_Data/survey_products")
    check(sp.root_for_c(r, 0) == Path(str(r) + "_c0") and sp.root_for_c(r, sp.SETUP_COEF) == r
          and sp.root_for_c(Path(str(r) + "_c0"), sp.SETUP_COEF) == r
          and sp.root_for_c(Path(str(r) + "_c0"), 0.05) == Path(str(r) + "_c0.05")
          and sp.root_for_c(Path(str(r) + "_c0.034"), 0) == Path(str(r) + "_c0"),
          "root for a C: survey_products -> _c0; survey_products_c0 + the station's C -> survey_products (never "
          "survey_products_c0_c0.037)")
    out_c0 = d / "root_c0" / "2026-09-29"
    out_c0.mkdir(parents=True)
    sp.write_stamp(out_c0, "merge", {"x": 1}, extra={"setup": {
        c: {"coef": 0.0, "fit": sp.fit_sig(sp.setup_fit_info(0.0))} for c in ("c1", "c2")}})
    a37 = SimpleNamespace(setup_coef=sp.SETUP_COEF, setup_fitted_to=None, output_root=str(d / "root_c0"))
    t37 = sp.setup_replace_note(dict(plan, date="2026-09-29", era="live", out=out_c0), a37, sp.STEPS)
    check(t37 and f"--output-root {d / 'root'} for C = {sp.SETUP_COEF:g}" in t37 and "root_c0_c" not in t37
          and "NOTE" not in t37,
          "a refresh without --setup-coef over a C = 0 root: the station-C product's root is the base root")
    # a root whose name does not say what it holds: the advice never points back at the root that refused
    check(sp.root_for_c(r, sp.SETUP_COEF, [0.0]) == Path(str(r) + "_c0.037")
          and sp.root_for_c(Path(str(r) + "_c0"), 0, [sp.SETUP_COEF]) == Path(str(r) + "_c0_c0")
          and sp.root_for_c(Path(str(r) + "_c0"), sp.SETUP_COEF, [0.0]) == r
          and sp.root_for_c(r, 0, [sp.SETUP_COEF]) == Path(str(r) + "_c0"),
          "root for a C, given what the root holds: a root named for that C but holding another -> a root "
          "beside it; a root named for what it holds -> the name rule as before")
    out_p = d / "plain" / "2026-09-29"                       # no suffix, its first build used --setup-coef 0
    out_p.mkdir(parents=True)
    sp.write_stamp(out_p, "merge", {"x": 1}, extra={"setup": {
        c: {"coef": 0.0, "fit": sp.fit_sig(sp.setup_fit_info(0.0))} for c in ("c1", "c2")}})
    ap_ = SimpleNamespace(setup_coef=sp.SETUP_COEF, setup_fitted_to=None, output_root=str(d / "plain"))
    tp = sp.setup_replace_note(dict(plan, date="2026-09-29", era="live", out=out_p), ap_, sp.STEPS)
    check(tp and f"--output-root {d / 'plain_c0.037'} for C = {sp.SETUP_COEF:g}" in tp
          and f"--output-root {d / 'plain'} for" not in tp
          and f"the name of {d / 'plain'} gives it to C = {sp.SETUP_COEF:g}, but it holds the C = 0 product" in tp
          and f"--setup-coef 0 --output-root {d / 'plain_c0'}), then rerun here with --replace" in tp,
          "a root without suffix holding C = 0, default C: refused, a root beside it named (not the same root), "
          "and what it holds and how to free it")
    out_m = d / "misplaced_c0" / "2026-09-29"                # a ..._c0 root holding C = 0.037
    out_m.mkdir(parents=True)
    sp.write_stamp(out_m, "merge", {"x": 1}, extra={"setup": {
        c: {"coef": sp.SETUP_COEF, "fit": sp.fit_sig(sp.setup_fit_info(sp.SETUP_COEF))} for c in ("c1", "c2")}})
    am_ = SimpleNamespace(setup_coef=0.0, setup_fitted_to=None, output_root=str(d / "misplaced_c0"))
    tm = sp.setup_replace_note(dict(plan, date="2026-09-29", era="live", out=out_m), am_, sp.STEPS)
    check(tm and f"--output-root {d / 'misplaced_c0_c0'} for C = 0" in tm
          and f"--output-root {d / 'misplaced_c0'} for" not in tm
          and f"but it holds the C = {sp.SETUP_COEF:g} product" in tm
          and f"--output-root {d / 'misplaced'}), then rerun here with --replace" in tm,
          "a ..._c0 root holding C = 0.037, --setup-coef 0: a root beside it, and where its C = 0.037 belongs")
    cav_m = sp.setup_caveats({"date": "2026-09-29", "out": out_m, "surveys": [
        {"name": "rtk", "path": "2026-09-29_Marconi_Checkshots.csv", "survey_type": "points",
         "survey_date": "2026-09-29"}]}, SimpleNamespace(survey_dirs=[str(TEST_DIRS["surveys"])], setup_fitted_to=None),
        {}, {}, [dict(sp.setup_fit_info(sp.SETUP_COEF), cameras=["c1", "c2"])], "live")
    circ_m = [c for c in cav_m if c.startswith("rtk: CIRCULAR")]
    check(circ_m and f"--setup-coef 0 --output-root {d / 'misplaced_c0_c0'}" in circ_m[0]
          and f"but it holds the C = {sp.SETUP_COEF:g} product" in circ_m[0],
          "the CIRCULAR caveat in a ..._c0 root holding C = 0.037: the C = 0 build goes beside it, not here")
    check(sp.root_name_note(d / "plain", 0.0, [sp.SETUP_COEF]) == ""
          and sp.root_name_note(d / "root_c0", sp.SETUP_COEF, [0.0]) == ""
          and sp.root_name_note(d / "plain", sp.SETUP_COEF, [sp.SETUP_COEF]) == "",
          "no note where the name gives the root to what it holds, or the root holds the C asked for")
    # exit codes: refused is 5 (not 2, a disabled date), and ranks after failed in --all
    check(sp.EXIT_CODES["refused"] == 5 and sp.worst_exit(["built", "refused", "busy", "partial"]) == 5
          and sp.worst_exit(["refused", "failed"]) == 1, "refused: exit 5; --all: failed > refused > busy > partial")
    # a usage error has its own code (6), never 2 (a disabled date) and never a date result's
    codes = {k: v for k, v in sp.EXIT_CODES.items() if k != "usage"}
    check(sp.EXIT_CODES["usage"] == 6 and 6 not in codes.values(), "usage error: exit 6, used by no date result")
    for bad_args in (["--date", "2025-01-23", "--steps", "bogus"], ["--date", "2025-01-23", "--setup-coef", "-1"],
                     ["--date", "2025-01-23", "--window", "2025-01-23", "2025-01-18"],
                     ["--date", "2025-01-23", "--utc-hours", "x"], ["--no-such-option"]):
        r_ = subprocess.run([sys.executable, str(HERE / "survey_products.py")] + bad_args
                            + ["--output-root", str(d / "usage")], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True, cwd=str(d))
        check(r_.returncode == 6 and "NOTHING built" in r_.stdout and not (d / "usage").exists(),
              f"usage error {' '.join(bad_args)}: exit 6, nothing written (rc {r_.returncode})")
    r_ = subprocess.run([sys.executable, str(HERE / "survey_products.py"), "--date", "2025-01-23", "--setup-coef",
                         "0.05", "--output-root", str(d / "unk")], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                        universal_newlines=True, cwd=str(d))
    check(r_.returncode == 5 and "where was this C fitted" in r_.stdout and not (d / "unk").exists(),
          f"a C of unknown origin is refused with exit 5 (not argparse's 2), nothing written (rc {r_.returncode})")
    # a full rebuild with another window in the same root: says it replaces, names a recheck root
    st = d / "root_w" / "2025-01-23"
    st.mkdir(parents=True)
    sp.write_stamp(st, "detect_c1", {"cmd": "detect_original_view.py --start-date 2025-01-18 --end-date 2025-01-23 "
                                            "--utc-hours 13.5-18"})
    pw = {"date": "2025-01-23", "era": "adcp", "out": st, "first": "2025-01-20", "last": "2025-01-23",
          "cams": {"c1": {"first": "2025-01-20", "last": "2025-01-23", "hours": (13.5, 18.0),
                          "hours_text": "13.5-18"}}}
    aw = SimpleNamespace(output_root=str(d / "root_w"))
    tw = sp.window_replace_note(pw, aw, sp.STEPS)
    check(tw and tw.startswith(f"replacing the product in {st} built for another window/hours")
          and "c1 window 2025-01-20 .. 2025-01-23 (built: 2025-01-18 .. 2025-01-23)" in tw
          and str(d / "root_w_recheck") in tw and sp.window_replace_note(pw, aw, ("compare",)) is None,
          "another window, all steps: 'replacing the product ... built for another window', a recheck root named; "
          "--steps compare keeps the built window (no note)")
    check(sp.setup_replace_note(dict(plan, era="live"), a0, sp.STEPS) and
          "detection runs again" not in sp.setup_replace_note(dict(plan, era="live"), a0, sp.STEPS),
          "... the live era has no detection to rerun")
    same = SimpleNamespace(setup_coef=sp.SETUP_COEF, setup_fitted_to=None, output_root=str(d / "root"))
    check(sp.setup_replace_note(plan, same, sp.STEPS) is None
          and sp.setup_replace_note(plan, a0, ("compare",)) is None,
          "the same C, or a run that does not rebuild the waterlines (--steps compare): no warning")
    # the product README's advice and the fit recipe build another C in its own root
    sw = sp.still_water_caveat(sp.SETUP_COEF, {}, "adcp")
    check("OWN --output-root" in sw and "--setup-coef 0 --output-root <root>_c0" in sw and "curve is noisy" in sw
          and "both lie outside its near-best band" in sw and "+0.006 to +0.034" in sw
          and "-0.004 to +0.008" not in sw and "Otherwise keep the station's C" in sw,
          "the fit-C recipe: each build in its own root; --fit-setup offers a C only outside its near-best band, "
          "with the bands it printed")
    cav = sp.setup_caveats({"date": "2026-09-29", "out": out, "surveys": [
        {"name": "rtk", "path": "2026-09-29_Marconi_Checkshots.csv", "survey_type": "points",
         "survey_date": "2026-09-29"}]}, SimpleNamespace(survey_dirs=[str(TEST_DIRS["surveys"])], setup_fitted_to=None),
        {}, {}, [dict(sp.setup_fit_info(sp.SETUP_COEF), cameras=["c1", "c2"])], "live",
        prov={"comparisons": [{"name": "rtk", "current": True, "headline": {"waterlines": {"above_survey_floor": {
            "frames_all": 16, "frames": 7, "median": -0.069, "margin_m": 0.1, "frames_below_floor": 6}}}}]})
    circ = [c for c in cav if c.startswith("rtk: CIRCULAR")]
    check(circ and "--setup-coef 0 --output-root" in circ[0] and "root_c0" in circ[0]
          and "would replace this product" in circ[0],
          "the CIRCULAR caveat: a C = 0 build, with the command and its own root")
    under = [c for c in cav if c.startswith(f"C = {sp.SETUP_COEF} is not pinned")]
    check(under and "The Oct 2026 recheck" in under[0] and "not this build's" in under[0]
          and "outside this repository's code" in under[0]
          and "THIS BUILD's own check of that kind (rtk" in under[0] and "median -0.069 m on 7 frames" in under[0]
          and "6 of 16 frames below the lowest shot" in under[0] and "read this build's" in under[0],
          "the recheck's fixed numbers say where they come from, next to this build's own floor statistic")


def test_filter_leaves_nothing_old(d):
    print("filter step: last build's report and diagnostic plots never sit next to this build's")
    d = Path(d) / "filt"
    w = d / "waterlines"
    w.mkdir(parents=True)
    # a camera the new run leaves out, a report and a filtered copy from last build, and their .tmp copies
    old = [w / "consistency_c1.png", w / "consistency_c2.png", w / "consistency_report.csv",
           w / "consistency_report.csv.tmp", w / "contour_points_ground_filtered.csv",
           w / "contour_points_ground_filtered.csv.tmp"]
    for f in old:
        f.write_text("last build\n")
    keep = w / "contour_points_ground.csv"
    keep.write_text("camera,capture_epoch\n")              # no pixel columns: the filter exits 1
    plan = {"out": d, "first": "2026-09-26", "last": "2026-10-02"}
    state = {}
    res, _ = quiet(sp.step_filter, plan, SimpleNamespace(force=False), state)
    gone = [f.name for f in old if f.exists()]
    check(res == "failed" and state.get("filter_failed") and not gone and keep.exists(),
          f"a failed filter run leaves none of last build's report, plots or filtered file (left: {gone})")
    for f in old[:3]:
        f.write_text("last build\n")
    sp.remove_filter_outputs(w, w / "contour_points_ground_filtered.csv")
    check(not any(f.exists() for f in old) and keep.exists(),
          "before the filter runs: the old report and every consistency_<cam>.png are removed, the input kept")


def test_one_setup_coef():
    print("one setup coefficient for the station and the products")
    import historical_forcing as hf
    check(abs(hf.SETUP_COEF - sp.SETUP_COEF) < 1e-12,
          f"survey_products {sp.SETUP_COEF}, historical_forcing {hf.SETUP_COEF}: the same C "
          f"(the live cron's: test_live_cron, with the station's state)")
    check(any(abs(k - sp.SETUP_COEF) < 1e-12 for k in sp.SETUP_FITS),
          "the station's C has a SETUP_FITS entry (where it was fitted: which comparisons are CIRCULAR)")


def test_live_cron():
    """The cron the station runs, as it is on this computer (waterline_timex_cron.sh, which the owner may edit:
    its own comments say SETUP_COEF="" for a refit), read as survey_products.main() reads it (no value or ''
    is C = 0; a value that is not a number is named): the live waterlines and the products must carry the
    same C, or main() warns. Run with the station's tables, last."""
    print("the live cron's setup coefficient (waterline_timex_cron.sh as it is now)")
    raw = sp.cron_settings(sp.LIVE_CRON).get("SETUP_COEF")
    try:
        c = float(raw) if raw not in (None, "") else 0.0
    except ValueError:
        c = None
    check(c is not None and abs(c - sp.SETUP_COEF) < 1e-12,
          f"{sp.LIVE_CRON.name} SETUP_COEF {raw!r} "
          + ("(not a number)" if c is None else f"(C = {c:g})") + f", survey_products {sp.SETUP_COEF}: the same C")


def test_station_state(d):
    """The station's own state, last, so it cannot hide the code's results: its tables (test_station_tables)
    and its live cron (test_live_cron), each run even when the other fails, every failure raised at the end."""
    failed = []
    for fn, a in ((test_station_tables, (d,)), (test_live_cron, ())):
        try:
            fn(*a)
        except AssertionError as exc:
            print(f"  FAIL  {exc}")
            failed.append(str(exc))
    if failed:
        raise AssertionError(" || ".join(failed))


def test_py38():
    print("Python 3.8 syntax")
    import ast
    for f in ("survey_products.py", "detect_original_view.py", "test_survey_products.py", "estimate_eo_rotation.py",
              "survey_compare.py"):
        src = (HERE / f).read_text(encoding="utf-8")
        try:
            ast.parse(src, feature_version=(3, 8))
        except TypeError:                        # Python 3.7 has no feature_version
            ast.parse(src)
        check(True, f"{f} parses as Python 3.8")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="keep the outputs, in a new folder inside this one")
    ap.add_argument("--survey-dirs", nargs="+", default=[],
                    help="more folders to look for the real 2026-09-29 RTK file in (after the station's)")
    args = ap.parse_args()
    RTK_DIRS.extend(os.path.abspath(x) for x in args.survey_dirs)
    if args.keep:                       # a new folder each run: never the files an earlier run left
        Path(args.keep).mkdir(parents=True, exist_ok=True)
    d = os.path.realpath(tempfile.mkdtemp(prefix="survey_products_test_", dir=args.keep))
    TEST_DIRS.update(out=Path(d) / "no_output_root", surveys=Path(d) / "no_surveys")
    TEST_DIRS["surveys"].mkdir()
    start = os.getcwd()
    os.chdir(d)                         # a file name never finds a file of the folder the test was started from
    print(f"test folder: {d} (the current folder during the run)")
    try:
        test_setup_identical(d)
        test_collect_photos(d)
        test_labels(d)
        test_setup_in_use(d)
        test_merge_drops_no_setup(d)
        test_keep_moved_days(d)
        test_pointing_verdicts(d)
        test_storm_and_wave_caveats(d)
        test_feature_runner_up(d)
        test_live_rows(d)
        test_config(d)
        test_table_commas(d)
        test_frames_per_day(d)
        test_caveats(d)
        test_stamps(d)
        test_round2(d)
        test_fit_recognised(d)
        test_cwd_never_searched(d)
        test_fix1(d)
        test_fix2(d)
        test_fix3(d)
        test_dem_two_days(d)
        test_filter_leaves_nothing_old(d)
        test_setup_share(d)
        test_c_change(d)
        test_fit_setup_flat(d)
        test_one_setup_coef()
        test_py38()
        # the tables last: the repository's as shipped, then the station's as they are now (with its live
        # cron), so a table or cron the owner has edited cannot stop the checks of the code before they have run
        test_tables_changed_note(d)
        test_shipped_tables(d)
        test_station_state(d)
    finally:
        os.chdir(start)
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    why = {}
    for what, reason in SKIPPED:
        why.setdefault(reason, []).append(what)
    print("all survey_products tests passed" + (f" ({len(SKIPPED)} skipped: " + "; ".join(
        f"{', '.join(w)}: {r}" for r, w in why.items()) + ")" if SKIPPED else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())

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
    surveys.csv gives -- a row of the fit's file name, or the points row of the fit's survey date under
    whatever name the owner gave it --, by the product's own rule, survey_products.resolve_file: a path
    with a folder in it as given, a bare name in the survey folders; then that name, and the fit's own
    file name, in every survey folder) and confirmed by its sha256, the product's rule for 'the same
    shots'. -> (path, "") or (None, why it cannot be used)."""
    fit = sp.SETUP_FITS[sp.SETUP_COEF]
    try:
        rows, _ = quiet(sp.read_table, sp.DEFAULT_SURVEYS, sp.SURVEY_FIELDS)
    except (SystemExit, Exception):     # a table survey_products cannot read (refused, not UTF-8, bad csv):
        rows = []                       # test_station_tables, run last, says so
    given = []
    for r in rows:
        g = (r.get("path") or "").strip()
        if g and g not in given and (Path(g).stem == fit["survey_stem"] or (
                r.get("survey_type") == fit.get("survey_type") and r.get("survey_date") == fit.get("survey_date"))):
            given.append(g)
    if fit["survey_stem"] + ".csv" not in given:
        given.append(fit["survey_stem"] + ".csv")
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


@contextlib.contextmanager
def unblinded(date="2024-10-23"):
    """The behaviour of a BLIND date once it is unblinded (its BLIND_DATES entry removed, after the holder
    has reported): the survey looked for, merged, compared. Restored afterwards, whatever happens."""
    held = sp.BLIND_DATES.pop(date, None)
    try:
        yield
    finally:
        if held is not None:
            sp.BLIND_DATES[date] = held


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
    # the same shots cannot be overruled: by their name, always (the file need not be here); by their
    # content on the test's own fit (test_fit_recognised), and on the real shots when they are here
    v = sp.Verdict(dict(rtk, label_override_reason="anything"))
    sp.honest_label(rtk, p, args, rtkf, verdict=v)
    check(v.label == "CIRCULAR" and not v.overridden,
          "the same shots by name (2026-09-29_Marconi_Checkshots.csv) cannot be overruled by label_override_reason")
    if real.exists():
        v = sp.Verdict(dict(rtk, label_override_reason="anything"))
        sp.honest_label(dict(rtk, survey_date=""), p, args, real, verdict=v)
        check(v.label == "CIRCULAR" and not v.overridden,
              "... nor the real shots renamed, by their content (station file)")
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
    # the DEM settings reported are those the DEM on disk was built with (here a DEM built before Oct 2026,
    # tagged EPSG:32619, the default then; the default since is 6348)
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
                            geotiff_epsg=sp.GEOTIFF_EPSG, setup_coef=sp.SETUP_COEF)
    b = sp.built_args(later, {"out": d, "date": "2025-01-23", "cams": {"c1": {}}})
    check(b.dem_cell == 2.0 and b.dem_max_hs == 1.5 and b.geotiff_epsg == 32619 and later.dem_cell == 1.0,
          "the 'rebuild everything' line takes the DEM settings the DEM was built with, not this run's")
    import survey_compare
    defaults = {"survey_products.GEOTIFF_EPSG": sp.GEOTIFF_EPSG, "survey_compare.DEFAULT_EPSG": survey_compare.DEFAULT_EPSG}
    for name in ("asc_to_geotiff.py", "process_chelsea.py"):          # their --epsg default, read, not imported
        m_ = re.search(r'add_argument\("--epsg", type=int, default=(\d+)', (HERE / name).read_text())
        defaults[f"{name} --epsg"] = int(m_.group(1)) if m_ else None
    check(set(defaults.values()) == {6348},
          f"GeoTIFFs are tagged EPSG:6348 (NAD83(2011) / UTM 19N, the surveys' frame) by default everywhere ({defaults})")
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
        v = sp.Verdict(dict(s, label_override_reason="anything"))
        sp.honest_label(s, p, here, renamed, fits=fits, verdict=v)
        check(v.label == "CIRCULAR" and not v.overridden,
              "... and label_override_reason cannot overrule it: the same shots by content (test_labels: by name)")
        # real_rtk(), which finds the real shots for the checks that need them, on this fit: the points row of
        # the fit's date under another name is looked up by the product's rule and confirmed by its sha256
        saved = (sp.SETUP_FITS[sp.SETUP_COEF], sp.DEFAULT_SURVEYS, list(RTK_DIRS))
        tbl = d / "surveys_renamed.csv"
        tbl.write_text("date,name,path,survey_type,survey_date,label,why\n"
                       "2026-10-15,shots,renamed_export.csv,points,2026-10-15,INDEPENDENT,the owner's own name\n")
        try:
            sp.SETUP_FITS[sp.SETUP_COEF], sp.DEFAULT_SURVEYS, RTK_DIRS[:] = sp.SETUP_FITS[fake], tbl, [str(elsewhere)]
            got, _ = real_rtk()
            check(got == elsewhere / "renamed_export.csv",
                  f"real_rtk: the fit's shots under the name surveys.csv's points row of its date gives, confirmed by "
                  f"sha256 ({got})")
            (elsewhere / "renamed_export.csv").write_text("Name,Easting,Northing,Elevation\n1,1,1,1\n")
            got, why_ = real_rtk()
            check(got is None and "sha256 is not the one" in why_,
                  "... a file of that name with other content: not used, and said so")
        finally:
            sp.SETUP_FITS[sp.SETUP_COEF], sp.DEFAULT_SURVEYS = saved[0], saved[1]
            RTK_DIRS[:] = saved[2]
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
    # themselves and merged, a target both cameras saw counted once. 2024-10-23 is BLIND (its GCPs are
    # held for an independent check: test_blind_and_calibration), so this is the behaviour once it is
    # unblinded (its BLIND_DATES entry removed).
    oct_ = [s for s in srv if s["name"] == "oct_gcps"][0]
    sdir = Path(d) / "gcp_arrive"
    sdir.mkdir()
    a = cfg_args(survey_dirs=[str(sdir)], **tabs)
    held = sp.BLIND_DATES.pop("2024-10-23", None)
    try:
        paths, note = sp.resolve_survey(oct_, None, a)
        check(not paths and "no survey file" in note, "Oct GCPs not there yet: skipped, with the reason")
        (sdir / "2024-10-23_Marconi_Extrinsic_Targets_c1_xyz.csv").write_text(
            "1,420100.0,4638400.0,1.20\n2,420110.0,4638410.0,0.80\n")
        (sdir / "2024-10-23_Marconi_Extrinsic_Targets_c2_xyz.csv").write_text(
            "2,420110.0,4638410.0,0.80\n3,420120.0,4638500.0,0.50\n")
        paths, note = sp.resolve_survey(oct_, None, a)
        check(len(paths) == 2 and "found" in note, "Oct GCPs arrive (once unblinded): both camera files found by "
                                                   "their date")
        rows = [ln for ln in Path(sp.merged_points(paths, sdir / "merged.csv")).read_text().splitlines() if ln]
        check(len(rows) == 3, f"Oct GCPs merged: 3 targets, the one both cameras saw once (got {len(rows)})")
    finally:
        if held:
            sp.BLIND_DATES["2024-10-23"] = held


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
    for date in ("2025-03-19",):
        plan, out = quiet(sp.build_plan, date, cfg, srv, cfg_args(**tabs))
        forced, _ = quiet(sp.build_plan, date, cfg, srv, cfg_args(force_disabled=True, **tabs))
        check(plan is None and "DISABLED" in out and forced is not None,
              f"{date}: switched off (no photos, no ADCP): refused, and builds with --force-disabled")
    plan, out = quiet(sp.build_plan, "2024-10-23", cfg, srv, cfg_args(**tabs))
    check(plan is not None and plan.get("blind") and plan["surveys"] == [] and not plan["forced_disabled"]
          and [h["name"] for h in plan["held"]] == ["oct_gcps"]
          and all(c["eo_file"] == f"CACO03_{k}_20241023_EO.yaml" for k, c in plan["cams"].items()),
          "2024-10-23: switched on (the user asked for it, 13 Oct 2026), BLIND: DEM only, its survey withheld, "
          "calibration A (the table's) unless --eo")
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
    check(oct_["label"] == "PARTLY-CIRCULAR" and not oct_["path"] and "BLIND" in oct_["why"]
          and "Sherwood" in oct_["why"],
          "Oct GCPs: PARTLY-CIRCULAR (the label of a check of the DEM built with calibration A), no path: held by "
          "Chris Sherwood (BLIND), never on the NUC")


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
        check("6348" in text and ("NOMINAL" in text) == (sp.GEOTIFF_EPSG == 32619),
              f"{era}: the GeoTIFF CRS tag (EPSG:{sp.GEOTIFF_EPSG}) is said, nominal only if 32619")
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
          "a result of the rule before Oct 2026 whose 'most days differ' guard kept every day: none (as its "
          "detection did)")
    check(all(sp.pointing_leave_out([{"date": "2025-01-2%d" % k, "verdict": v} for k, v in enumerate(vs)])["leave_out"]
              for vs in (["DIFFERENT POINTING"] * 3 + ["matches"], ["UNCONFIRMED"] * 4, [])),
          "the rule of Oct 2026 never keeps days automatically ('leave_out' always True)")


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
    # -- a day nothing could read (horizon fit failed, features unmatched) next to a left-out day (review 2,
    # simR: c1 28 Sep's 3-frame fit failed at 14.6 px; 'unchecked', kept and gridded while c1 was turned):
    # UNCONFIRMED, and runs of such days with it; one with no left-out neighbour stays 'unchecked'
    blindrows = _nuc_rows("c1")
    by(blindrows)["2026-09-28"].update(horizon="unclear", horizon_fit_px=14.6)
    by(blindrows)["2026-10-01"].update(horizon="unclear", horizon_fit_px=20.0, features="unmatched", inliers=5,
                                       d_azimuth_deg="", d_tilt_deg="", d_roll_deg="")
    sp.pointing_verdicts(blindrows, ref_day="2026-09-29")
    vb = by(blindrows)
    check(vb["2026-09-28"]["verdict"] == "UNCONFIRMED" and "could not be fitted (residual 14.6 px)" in vb["2026-09-28"]["why"]
          and "the day next to it, 2026-09-27, was measured to differ from the reference" in vb["2026-09-28"]["why"]
          and vb["2026-10-01"]["verdict"] == "unchecked" and "left out" not in vb["2026-09-28"]["why"],
          "a day nothing could read next to a left-out day (c1 28 Sep, its fit failed): UNCONFIRMED; one between "
          "matching days (1 Oct) stays 'unchecked'")
    by(blindrows)["2026-09-29"].update(horizon="unclear", horizon_fit_px=20.0, features="unmatched", inliers=5,
                                       d_azimuth_deg="", d_tilt_deg="", d_roll_deg="")
    sp.pointing_verdicts(blindrows)
    check([vb[x]["verdict"] for x in ("2026-09-28", "2026-09-29", "2026-09-30")]
          == ["UNCONFIRMED", "UNCONFIRMED", "matches"]
          and "it is in a run of days nothing could read next to 2026-09-27, which was measured to differ" in vb["2026-09-29"]["why"]
          and sp.pointing_leave_out(blindrows)["next_to_days"] == ["2026-09-28", "2026-09-29"],
          "... a run of such days after a day measured to differ: each UNCONFIRMED, its why naming the day the run "
          "starts at")
    # ... but never the reference day: the other days are measured against its frame, so a neighbour that
    # differs from it says nothing of it (its run stops there)
    sp.pointing_verdicts(blindrows, ref_day="2026-09-29")
    check([vb[x]["verdict"] for x in ("2026-09-28", "2026-09-29", "2026-09-30")]
          == ["UNCONFIRMED", "unchecked", "matches"],
          "... the reference day itself, nothing read on it and next to a left-out day: 'unchecked', not UNCONFIRMED "
          f"({vb['2026-09-29']['verdict']})")
    # ... but never from a day UNCONFIRMED on a horizon reading alone that is not a clear horizon (review 3: one
    # misread fog edge with unmatched features, 19 Jan, took the three fog days after it, whose fits failed, with
    # it: 4 of 7 days left out, where the same days next to a matching day stay 'unchecked'): such a reading
    # alone never shows a day moved, so it does not reach past its own day
    STORMY = (6.1, 0.72, 0.00, 62)
    fog = [dict(_pt_row("c2", "2025-01-18", "off the calibration", 34.0, 0.56, -0.68, 0.3, "ok", 2000,
                        (0.004, 0.003, 0.01)), turn_sigma=30.0),
           _pt_row("c2", "2025-01-19", "off the calibration", 34.0, 2.6, 1.2, 3.0, "unmatched", 4, edge=STORMY)]
    fog += [_pt_row("c2", f"2025-01-{dd}", "unclear", 34.0, t, r_, f_, "unmatched", n_, edge=STORMY)
            for dd, t, r_, f_, n_ in (("20", 0.9, -0.2, 25.0, 6), ("21", 0.4, -0.9, 30.0, 3), ("22", 1.4, 0.2, 18.0, 7))]
    fog += [dict(_pt_row("c2", "2025-01-23", "off the calibration", 34.0, 0.56, -0.68, 0.3, "ok", 3000), turn_sigma=2.0),
            _pt_row("c2", "2025-01-24", "unclear", 34.0, 0.7, -0.5, 22.0, "unmatched", 2, edge=STORMY)]
    sp.pointing_verdicts(fog, ref_day="2025-01-23")
    vf = by(fog)
    lof = dict(sp.pointing_leave_out(fog), rows=fog)
    check([vf[f"2025-01-{dd}"]["verdict"] for dd in range(18, 25)]
          == ["matches", "UNCONFIRMED", "unchecked", "unchecked", "unchecked", "matches", "unchecked"]
          and sp.skip_days_of(lof, SimpleNamespace(keep_moved_days=False)) == ["2025-01-19"]
          and lof["horizon_alone_days"] == ["2025-01-19"] and not lof["next_to_days"],
          "a fog week (review 3): the day UNCONFIRMED on a misread horizon alone is left out, the days nothing read "
          "next to it stay 'unchecked', kept (were UNCONFIRMED: 4 of 7 days left out)")
    # ... while a day UNCONFIRMED on more (a clear-looking one-day jump with no feature reading) starts a run
    vf["2025-01-19"].update(**dict(zip(sp.EDGE_KEYS, EDGE_CLEAR["c2"])))
    sp.pointing_verdicts(fog, ref_day="2025-01-23")
    lof = dict(sp.pointing_leave_out(fog), rows=fog)
    check([vf[f"2025-01-{dd}"]["verdict"] for dd in range(19, 23)] == ["UNCONFIRMED"] * 4
          and "the day next to it, 2025-01-19, is UNCONFIRMED on a horizon jump on a clear-looking edge"
          in vf["2025-01-20"]["why"] and "a run of days nothing could read next to 2025-01-19" in vf["2025-01-22"]["why"]
          and vf["2025-01-24"]["verdict"] == "unchecked" and lof["next_to_days"] == ["2025-01-20", "2025-01-21",
                                                                                      "2025-01-22"]
          and "3 of them (2025-01-20, 2025-01-21, 2025-01-22) are UNCONFIRMED only because nothing could read them"
          in (sp.pointing_majority_note(lof, SimpleNamespace(keep_moved_days=False)) or ""),
          "... the same days next to a clear-looking one-day jump: UNCONFIRMED, a run from it; the warning counts "
          "them as left out for the day next to them")
    # ... and next to misread horizons on two days running, a storm (review 3's own S_* scenarios: real frames,
    # c1 turned 0.3-1.0 deg on all three storm days; 26/27 Sep misread with unmatched features, 28 Sep's fit
    # failed: without this 28 Sep was 'unchecked', kept, its frames turned)
    vf["2025-01-19"].update(**dict(zip(sp.EDGE_KEYS, STORMY)))
    vf["2025-01-20"].update(horizon="off the calibration", horizon_fit_px=2.0, horizon_dtilt_deg=3.1)
    sp.pointing_verdicts(fog, ref_day="2025-01-23")
    lof = dict(sp.pointing_leave_out(fog), rows=fog)
    check([vf[f"2025-01-{dd}"]["verdict"] for dd in range(19, 23)] == ["UNCONFIRMED"] * 4
          and lof["horizon_alone_days"] == ["2025-01-19", "2025-01-20"] and lof["next_to_days"] == ["2025-01-21",
                                                                                                    "2025-01-22"]
          and "the day next to it, 2025-01-20, is UNCONFIRMED on a misread horizon in a storm (misread horizons on "
              "2025-01-19, 2025-01-20, in one run of days left out or unread; a camera can turn in one" in vf["2025-01-21"]["why"]
          and vf["2025-01-24"]["verdict"] == "unchecked",
          "... misread horizons on two days running (a storm): the days nothing read after them UNCONFIRMED too")
    # -- review 4: a storm's days come misread or unread in any order, and the storm run is built over the days
    # left out on their own readings OR read by nothing (v9 needed two misread days side by side; v8 left these
    # days out). The NUC's c1 rows with 27 Sep unread (its fit failed at 33.0 px, no features): misread, unread,
    # misread: 27 Sep is UNCONFIRMED (v9: 'unchecked', kept), with or without a wave record
    nu = _nuc_rows("c1")
    by(nu)["2026-09-27"].update(horizon="unclear", horizon_fit_px=33.0, features="unmatched", inliers=0,
                                d_azimuth_deg="", d_tilt_deg="", d_roll_deg="")
    sp.pointing_verdicts(nu, ref_day="2026-09-29")
    lon = sp.pointing_leave_out(nu)
    check([by(nu)[x]["verdict"] for x in ("2026-09-26", "2026-09-27", "2026-09-28")] == ["UNCONFIRMED"] * 3
          and lon["different_days"] == ["2026-09-26", "2026-09-27", "2026-09-28"] and lon["next_to_days"] == ["2026-09-27"]
          and lon["horizon_alone_days"] == ["2026-09-26", "2026-09-28"]
          and "the day next to it, 2026-09-26, is UNCONFIRMED on a misread horizon in a storm (misread horizons on "
              "2026-09-26, 2026-09-28, in one run of days left out or unread; a camera can turn" in by(nu)["2026-09-27"]["why"],
          "review 4: the NUC's c1 rows with 27 Sep unread between the misread 26 and 28 Sep: a storm, 27 Sep "
          "UNCONFIRMED (v9 kept it)")
    # ... ST10_S3 (review 4: the real c1 storm frames tilted +1.0 deg on 26-28 Sep; 26 Sep misread, 27/28 Sep fits
    # failed at 33.0 / 10.7 px with 6 / 2 inliers): ONE misread day, then unread days, the shape of review 3's calm
    # fog week. The wave record tells them apart: inside STORMS or offshore Hs over STORM_HS_M (1.5 m, the DEM's
    # own limit) on a day of the run, the unread days are UNCONFIRMED; on calm days they stay 'unchecked'
    st = _nuc_rows("c1")
    by(st)["2026-09-26"].update(horizon_dtilt_deg=5.20, horizon_droll_deg=0.15, horizon_offset_px=240.0,
                                **dict(zip(sp.EDGE_KEYS, (5.5, 0.55, 0.0, 12))))
    for day, fit_, n_ in (("2026-09-27", 33.0, 6), ("2026-09-28", 10.7, 2)):
        by(st)[day].update(horizon="unclear", horizon_fit_px=fit_, features="unmatched", inliers=n_,
                           d_azimuth_deg="", d_tilt_deg="", d_roll_deg="")
    got = {}
    for name, storms in (("calm", None), ("STORMS", {"2026-09-26": "inside the 25-26 Sep 2026 storm"}),
                         ("waves 27", {"2026-09-27": "offshore Hs up to 3.71 m, over 1.5 m (forcing/waves.csv)"}),
                         ("waves 30", {"2026-09-30": "offshore Hs up to 1.60 m, over 1.5 m (forcing/waves.csv)"})):
        rows = [dict(r) for r in st]
        sp.pointing_verdicts(rows, ref_day="2026-09-29", storm_days=storms)
        got[name] = ([by(rows)[x]["verdict"] for x in ("2026-09-26", "2026-09-27", "2026-09-28")], by(rows))
    check(got["calm"][0] == ["UNCONFIRMED", "unchecked", "unchecked"]
          and got["STORMS"][0] == ["UNCONFIRMED"] * 3 and got["waves 27"][0] == ["UNCONFIRMED"] * 3
          and got["waves 30"][0] == ["UNCONFIRMED", "unchecked", "unchecked"]
          and "the day next to it, 2026-09-26, is UNCONFIRMED on a misread horizon in a storm (2026-09-26: inside the "
              "25-26 Sep 2026 storm; a camera can turn in one" in got["STORMS"][1]["2026-09-27"]["why"]
          and "it is in a run of days nothing could read next to 2026-09-26" in got["STORMS"][1]["2026-09-28"]["why"]
          and "(2026-09-27: offshore Hs up to 3.71 m" in got["waves 27"][1]["2026-09-28"]["why"],
          "review 4, ST10_S3: one misread day, then days nothing read: UNCONFIRMED when the run holds a storm day "
          "(STORMS, or the waves over 1.5 m on any of its days), 'unchecked' on calm days (as review 3's fog week) "
          "and when the high waves fall on a day outside the run")
    # ... and the forcing's waves give the storm days: STORMS, and offshore Hs over 1.5 m (UTC days)
    wd = Path(d) / "stormdays"
    (wd / "forcing").mkdir(parents=True)
    with open(wd / "forcing" / "waves.csv", "w") as f:
        f.write("time_utc,epoch,wvht_m,dpd_s,source\n")
        for day, hs in ((26, 1.2), (27, 3.71), (28, 1.5), (29, 1.52), (30, float("nan"))):
            e = datetime(2026, 9, day, 17, tzinfo=timezone.utc).timestamp()
            f.write(f"x,{e:.0f},{hs},8.0,buoy\n")
            f.write(f"x,{e - 3600:.0f},0.4,8.0,buoy\n")
    sd_ = sp.storm_days_of({"out": wd}, ["2026-09-%02d" % k for k in range(25, 31)])
    sd0 = sp.storm_days_of({"out": Path(d) / "no_forcing"}, ["2026-09-25", "2026-09-26", "2026-09-27"])
    check(sorted(sd_) == ["2026-09-25", "2026-09-26", "2026-09-27", "2026-09-29"]
          and sd_["2026-09-27"] == "offshore Hs up to 3.71 m, over 1.5 m (forcing/waves.csv)"
          and sd_["2026-09-26"] == "inside the 25-26 Sep 2026 storm" and sorted(sd0) == ["2026-09-25", "2026-09-26"],
          "storm_days_of: the STORMS days and the days with offshore Hs over 1.5 m (28 Sep's 1.50 is not over it); "
          "with no forcing, the STORMS days only")
    # -- review 4, P10_30_FOG01 (real frames: a 1.0 deg pan on 30 Sep and 1 Oct, 1 Oct fogged over the land): 30
    # Sep measured (+1.00 deg on 1851 inliers), 1 Oct 'horizon only' at the anchor, kept with the pan on both
    # cameras: a horizon cannot see a pan, so nothing shows the camera back. The tilt version (T06_30_FOG01, 1 Oct
    # reading the same jump) was caught; a horizon that shows a tilt back is believed (it sees tilt)
    pf = _nuc_rows("c1")
    by(pf)["2026-09-30"].update(features="MOVED", inliers=1851, d_azimuth_deg=1.0, d_tilt_deg=0.01, d_roll_deg=0.0,
                                turn_sigma=5613.3)
    by(pf)["2026-10-01"].update(features="unmatched", inliers=0, d_azimuth_deg="", d_tilt_deg="", d_roll_deg="")
    variants = {}
    for name, upd in (("horizon only", {}),
                      ("misread near the anchor", dict(zip(sp.EDGE_KEYS, (6.0, 0.7, 0.0, 30)))),
                      ("fit failed", {"horizon": "unclear", "horizon_fit_px": 30.0}),
                      ("tilt, horizon only", "tilt")):
        rows = [dict(r) for r in pf]
        if upd == "tilt":
            by(rows)["2026-09-30"].update(d_azimuth_deg=0.0, d_tilt_deg=0.6, horizon_dtilt_deg=1.55)
        else:
            by(rows)["2026-10-01"].update(upd)
        sp.pointing_verdicts(rows, ref_day="2026-09-29")
        variants[name] = by(rows)
    w1 = variants["horizon only"]["2026-10-01"]["why"]
    check(all(variants[x]["2026-10-01"]["verdict"] == "UNCONFIRMED" for x in ("horizon only", "misread near the anchor",
                                                                              "fit failed"))
          and variants["horizon only"]["2026-10-02"]["verdict"] == "matches"
          and w1.startswith("tilt and roll within 0.3 deg of the reference's pointing")
          and "the day next to it, 2026-09-30, was measured to differ from the reference (a pan of +1.00 deg, which no "
              "horizon sees): its clear horizon shows tilt and roll back but cannot see the pan" in w1
          and "nothing measured the pointing (kept)" not in variants["misread near the anchor"]["2026-10-01"]["why"]
          and variants["tilt, horizon only"]["2026-10-01"]["verdict"] == "horizon only"
          and variants["tilt, horizon only"]["2026-09-30"]["verdict"] == "DIFFERENT POINTING",
          "review 4, P10_30_FOG01: next to a measured 1.0 deg pan, a 'horizon only' day (and one whose fitted horizon "
          "is not clear, or failed) is UNCONFIRMED; next to a measured tilt with no pan, a clear horizon at the "
          "reference's shows the camera back: 'horizon only', kept")
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
    # -- the zero-turn guard (review 2: window dirt matches at exactly zero turn when fog hides the scene). The
    # NUC's rows predate turn_sigma (an older row counts as the scene, as above); with the values the real storm
    # frames project for them (the real c2 26-28 Sep frames against the real 7 Oct frame, their standard errors
    # scaled to the NUC's inlier counts: ~17, ~9, ~9) nothing changes; at a dirt-like 2.0 the days whose horizon
    # jumps on a misread edge are UNCONFIRMED, and 28 Sep (its fit failed) next to them, two such days running
    # (a storm), the reference day and the clear days are not
    # (review 4: 28 Sep follows only when its own match is no reading either: at its projected ~9 its 393-inlier
    # match decides, and it matches)
    for sig, want in (({"2026-09-26": 17.5, "2026-09-27": 8.7, "2026-09-28": 8.4}, ["matches"] * 3),
                      ({"2026-09-26": 3.0, "2026-09-27": 3.0, "2026-09-28": 8.4},
                       ["UNCONFIRMED", "UNCONFIRMED", "matches"]),
                      ({"2026-09-26": 2.0, "2026-09-27": 2.0, "2026-09-28": 2.0}, ["UNCONFIRMED"] * 3)):
        c2s = _nuc_rows("c2")
        for r in c2s:
            r["turn_sigma"] = sig.get(r["date"], 2.0)          # 29 Sep - 2 Oct: dirt-like too
        sp.pointing_verdicts(c2s, ref_day="2026-09-29")
        vs = by(c2s)
        got = [vs[x]["verdict"] for x in ("2026-09-26", "2026-09-27", "2026-09-28")]
        check(got == want and vs["2026-09-29"]["verdict"] == "matches"
              and all(vs[x]["verdict"] == "horizon only" and "no further from no turn than their own noise" in
                      vs[x]["why"] for x in ("2026-09-30", "2026-10-01", "2026-10-02"))
              and (want[2] == "matches" or ("than their own noise (2.0 standard errors, under 4)" in vs["2026-09-27"]["why"]
                                            and "the day next to it, 2026-09-27, is UNCONFIRMED on a misread horizon "
                                                "in a storm (misread horizons on 2026-09-26, 2026-09-27" in
                                            vs["2026-09-28"]["why"])),
              f"NUC c2 with turn_sigma {sig['2026-09-26']}/{sig['2026-09-27']}/{sig['2026-09-28']} on 26-28 Sep: {got} "
              f"(28 Sep, its fit failed: " + ("its own 393-inlier match decides" if want[2] == "matches" else
                                             "a day nothing read, next to two days misread, a storm")
              + "); the reference day matches whatever its turn_sigma, the clear "
              "days at a dirt-like 2.0 are 'horizon only' (the match is no reading), kept")
    # ... the reviewer's DIRT2 (the real c1 storm frames, c1 turned -1.2 deg, the same synthetic drops on every
    # frame: 113 and 165 inliers at 0.00 deg, turn_sigma 1.1 and 2.7, the horizon misread) and DIRT2_ST10 (c2
    # tilted +1.0 deg on 28 Sep, its fit failed: 244 inliers at 0.00 deg, turn_sigma 2.1): UNCONFIRMED, left out
    # (the first version of the rule kept all three; 28 Sep through 26/27 Sep, misread two days running)
    dr = [_pt_row("c1", "2026-09-26", "off the calibration", 205.4, 4.32, 0.23, 1.9, "ok", 113, (0.001, 0.0, -0.001),
                  edge=EDGE_STORM[("c1", "2026-09-26")]),
          _pt_row("c1", "2026-09-27", "off the calibration", 246.2, 4.87, 0.17, 4.0, "ok", 165, (0.001, 0.0, 0.0),
                  edge=EDGE_STORM[("c1", "2026-09-27")]),
          _pt_row("c1", "2026-09-28", "unclear", 102.5, 2.56, 1.55, 34.1, "ok", 268, (0.0, 0.001, 0.0), edge="")]
    for r, ts in zip(dr, (1.1, 2.7, 2.7)):
        r["turn_sigma"] = ts
    dr += [dict(_pt_row("c1", day, "off the calibration", 41.4, 0.90, 0.24, 0.4, "ok", 3500), turn_sigma=2.5)
           for day in ("2026-09-29", "2026-09-30", "2026-10-01", "2026-10-02")]
    sp.pointing_verdicts(dr, ref_day="2026-09-29")
    vd = by(dr)
    check([vd[x]["verdict"] for x in ("2026-09-26", "2026-09-27", "2026-09-28")] == ["UNCONFIRMED"] * 3
          and vd["2026-09-29"]["verdict"] == "matches"
          and all(vd[x]["verdict"] == "horizon only" for x in ("2026-09-30", "2026-10-01", "2026-10-02"))
          and "could not be fitted" in vd["2026-09-28"]["why"],
          "DIRT2: a dirt-only consensus (turn_sigma 1.1-2.7) against a misread horizon, or with no horizon fit next "
          "to two such days: UNCONFIRMED, left out; the clear days 'horizon only', kept")
    # -- printed: every day's why below the table, the left-out days with their verdicts
    res1 = dict(sp.pointing_leave_out(c1), rows=c1, info=dict(i1, calibration="CACO05_c1_20251113_EO-CV.yaml"))
    _, txt = quiet(sp.print_pointing, "c1", res1, SimpleNamespace(keep_moved_days=False))
    check("LEFT OUT of the detection: 2026-09-26 (UNCONFIRMED), 2026-09-27 (DIFFERENT POINTING), 2026-09-28 "
          "(UNCONFIRMED)" in txt and "    2026-09-26 UNCONFIRMED: the horizon read +4.23/+0.06" in txt,
          "print: the left-out days with their verdicts, and each day's why")
    # a pointing-only run (--steps pointing): the detection does not run with these results now
    _, txt = quiet(sp.print_pointing, "c1", res1, SimpleNamespace(keep_moved_days=False), detect_now=False,
                   detected={"skip_days": ["2026-09-27"]})
    check("to be LEFT OUT of the detection: 2026-09-26 (UNCONFIRMED), 2026-09-27 (DIFFERENT POINTING)" in txt
          and "when the detection runs again with these results (--steps detect and the steps after it; the "
              "detection on disk left out 2026-09-27)" in txt and "c1: LEFT OUT of the detection" not in txt,
          "print, --steps pointing alone: 'to be LEFT OUT' when the detection runs again, and what the detection on "
          "disk left out (it said 'LEFT OUT of the detection' as if done)")
    _, txt = quiet(sp.print_pointing, "c1", res1, SimpleNamespace(keep_moved_days=False), detect_now=False)
    check("no detection has run yet)" in txt, "... and with no detection on disk, says so")
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
          and i["horizon_overruled_days"] == ["2026-09-29"] and i["clear_horizon_days"] == 3
          and "the reference day 2026-09-29's own clear horizon reads +2.00/-0.60 deg" in (i["horizon_anchor_note"] or ""),
          "... on the reference day itself: matches (its own reading is the odd one), and no clear day; the step warns "
          "that the anchor leaves the reference day's own clear reading out")
    # -- a steady fog front that passes the edge tests on MORE days than the true horizon (review 2, simX2: a
    # planted clear edge at +130 px on c2 26, 27, 28 and 30 Sep, features matched at zero turn; 29 Sep the
    # reference, 1 and 2 Oct real clear frames): the reference day's group, confirmed by another day, is the
    # anchor, not the larger group; it was the front (1 and 2 Oct UNCONFIRMED and left out, the README's
    # offset 131 px, +1.96/+0.02 deg, where the true one is 34 px, +0.56/-0.67)
    rows = [_pt_row("c2", f"2026-09-{dd}", "off the calibration", 131.0, 1.96, 0.02, 0.4, "ok", 300, (0.0, 0.01, 0.0))
            for dd in ("26", "27", "28", "30")]
    rows += [_pt_row("c2", day, "off the calibration", 34.0, 0.56, -0.67, 0.2, "ok", 4000) for day in
             ("2026-09-29", "2026-10-01", "2026-10-02")]
    i = sp.pointing_verdicts(rows, ref_day="2026-09-29")
    v = by(rows)
    check(all(v[x]["verdict"] == "UNCONFIRMED" for x in ("2026-09-26", "2026-09-27", "2026-09-28", "2026-09-30"))
          and all(v[x]["verdict"] == "matches" for x in ("2026-09-29", "2026-10-01", "2026-10-02"))
          and i["horizon_anchor_deg"] == [0.56, -0.67] and i["clear_horizon_days"] == 3
          and i["window_horizon_offset_px"] == 34.0 and not i["horizon_anchor_note"],
          "a clear-looking front on 4 days against the reference day and 2 more: the reference's confirmed group is the "
          "anchor (34 px, +0.56/-0.67), the front's days UNCONFIRMED")
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
    # -- a real move with strong features (exact rotations of real frames read 0.30/0.39/0.50 deg on 64-106
    # inliers): measured, never overruled. A turned match is the scene and measures the change from the
    # reference directly: its tilt or roll is judged by the step's tilt/roll tolerance (0.3 deg), its pan by
    # pointing_check's 0.5, whatever the horizon reads (review 2: the real c2 27 Sep storm frame turned
    # +0.45 deg of tilt read +0.446 on 90 inliers against a misread horizon and was a 'small change', kept;
    # the deployed rule left it out)
    base = [_pt_row("c2", f"2026-10-0{k}", "off the calibration", 33.0, 0.55, -0.67, 0.2, "ok", 1400) for k in range(1, 5)]
    MISREAD = (6.0, 0.7, 0.0, 30)
    for feat, d_, hz_dt, edge, want in (
            ("MOVED", (0.0, 0.51, 0.0), 1.05, None, "DIFFERENT POINTING"),       # read 0.51: MOVED, a clear horizon too
            ("MOVED", (0.0, 0.51, 0.0), 1.05, MISREAD, "DIFFERENT POINTING"),    # ... a misread horizon
            ("MOVED", (0.0, 0.51, 0.0), 0.55, "", "DIFFERENT POINTING"),         # ... no horizon numbers at all
            ("moved", (0.0, 0.49, 0.0), 1.04, None, "DIFFERENT POINTING"),       # read 0.49: over 0.3 of tilt
            ("moved", (0.0, 0.49, 0.0), 1.04, MISREAD, "DIFFERENT POINTING"),    # ... the horizon misread: features
            ("moved", (-0.02, 0.446, -0.05), 4.22, MISREAD, "DIFFERENT POINTING"),   # review 2 ST045, c2 27 Sep
            ("moved", (-0.02, -0.005, 0.388), 3.89, MISREAD, "DIFFERENT POINTING"),  # review 2 SR045: a roll
            ("moved", (0.0, 0.35, 0.0), 0.90, None, "DIFFERENT POINTING"),       # 0.35 deg on a clear horizon
            ("moved", (0.0, 0.35, 0.0), 0.55, None, "DIFFERENT POINTING"),       # ... the clear horizon reads none
            ("moved", (0.0, 0.31, 0.0), 0.55, "", "DIFFERENT POINTING"),         # 0.31: no horizon numbers
            ("moved", (0.0, 0.29, 0.0), 0.84, None, "small change"),             # 0.29 deg: within 0.3
            ("moved", (0.0, 0.29, 0.0), 0.86, None, "DIFFERENT POINTING"),       # ... a clear 0.31 confirms (T03)
            ("moved", (0.0, 0.19, 0.0), 0.86, None, "small change"),             # ... 0.12 from 0.19: not confirmed
            ("moved", (0.0, 0.29, 0.0), 3.2, MISREAD, "small change"),           # ... a misread horizon overruled
            ("moved", (0.45, 0.05, 0.0), 0.60, None, "small change")):           # a 0.45 deg pan: pan by 0.5
        rows = [dict(r) for r in base] + [_pt_row("c2", "2026-10-05", "off the calibration", 33.0 + 66 * (hz_dt - 0.55),
                                                  hz_dt, -0.67, 0.2, feat, 1500, d_, edge=edge)]
        sp.pointing_verdicts(rows)
        r = by(rows)["2026-10-05"]
        check(r["verdict"] == want and ("1500 matched features" in r["why"]),
              f"a turn read {'/'.join(f'{x:+.3f}' for x in d_)} deg ({feat}) on 1500 features, horizon "
              f"{'without edge numbers' if edge == '' else 'misread' if edge else 'clear'} at {hz_dt - 0.55:+.2f} deg: "
              f"{want}" + ("" if want != "DIFFERENT POINTING" or feat == "MOVED" else
                           " (a tilt or roll over 0.3 deg from turned features)"))
    # -- a turned match (tilt/roll to ~0.1 deg: MATCH_ERR_DEG) that reads a change within 0.3 deg, on a day whose
    # CLEAR horizon (a few hundredths) reads one over 0.3 deg and agrees with the features within 0.1 deg in tilt
    # and roll (review 3: 0.28 on 50 inliers against a clear 0.37 was a 'small change', kept): DIFFERENT POINTING,
    # at any inlier count; with the horizon misread, or further than 0.1 deg from the features, the features
    # decide. Review 4: with a cut at 100 inliers the rule was not monotonic (0.34 confirmed a 0.05 deg tilt on
    # 50 inliers, 0.40 did not; 0.36 never confirmed 0.29 on 150): a larger jump never gives a weaker verdict now
    for n_, d_, hz_dt, edge, want in ((50, (0.0, 0.28, 0.0), 0.92, None, "DIFFERENT POINTING"),
                                      (99, (0.0, 0.25, 0.0), 0.86, None, "DIFFERENT POINTING"),
                                      (50, (0.0, 0.28, 0.0), 0.84, None, "small change"),     # the horizon reads 0.29
                                      (100, (0.0, 0.28, 0.0), 0.92, None, "DIFFERENT POINTING"),
                                      (150, (0.0, 0.29, 0.0), 0.91, None, "DIFFERENT POINTING"),  # review 4: 0.36
                                      (50, (0.0, 0.28, 0.0), 0.92, MISREAD, "small change"),
                                      (50, (0.3, 0.05, 0.0), 0.89, None, "small change"),     # review 4: 0.34 vs 0.05
                                      (50, (0.3, 0.05, 0.0), 0.95, None, "small change"),     # ... 0.40 vs 0.05
                                      (50, (0.0, 0.15, 0.0), 2.55, None, "small change")):    # contradicted: 2.0
        rows = [dict(r) for r in base] + [_pt_row("c2", "2026-10-05", "off the calibration", 33.0 + 66 * (hz_dt - 0.55),
                                                  hz_dt, -0.67, 0.2, "moved", n_, d_, edge=edge)]
        sp.pointing_verdicts(rows)
        r = by(rows)["2026-10-05"]
        check(r["verdict"] == want and (want == "small change" or "which the features confirm within 0.1 deg" in r["why"])
              and (want != "small change" or edge or hz_dt - 0.55 <= 0.3 or "not the change they measured" in r["why"]
                   or "not used" in r["why"]),
              f"a turn of {'/'.join(f'{x:+.2f}' for x in d_)} deg read on {n_} features, horizon "
              f"{'misread' if edge else 'clear'} at {hz_dt - 0.55:+.2f} deg: {want}")
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
    res = dict(lo, rows=rows, info=dict(i, reference_day="2026-09-29"))
    warn_ = sp.other_pointing_warning(res["info"], lo, "pointing/") or ""
    check([v[x]["verdict"] for x in ("2026-09-29", "2026-09-30", "2026-10-01", "2026-10-02")]
          == ["matches", "DIFFERENT POINTING", "DIFFERENT POINTING", "UNCONFIRMED"]
          and "where the reference day 2026-09-29 reads +0.91/+0.26" in v["2026-09-30"]["why"]
          and i["horizon_anchor_deg"] == [0.908, 0.255] and i["clear_horizon_days"] == 1
          and lo["measured_different_days"] == ["2026-09-30", "2026-10-01"] and lo["leave_out"]
          and i["other_pointing"] == {"days": ["2026-09-30", "2026-10-01"], "reference_days": ["2026-09-29"],
                                      "zero_turn_days": [],
                                      "change_deg": [None, 0.598, 0.0]}
          and sp.skip_days_of(res, SimpleNamespace(keep_moved_days=False))
          == ["2026-09-26", "2026-09-27", "2026-09-28", "2026-09-30", "2026-10-01", "2026-10-02"]
          and warn_.startswith("2 of the 3 days measured differ from the reference 2026-09-29 and agree with each "
                               "other (+0.60/+0.00 deg tilt/roll, pan not measured")
          and "the window holds two pointings" in warn_ and "--keep-moved-days keeps every day" in warn_,
          "a 0.6 deg tilt on half the clear days (features unmatched): measured against the reference day, not "
          "missed by a median that falls between the pointings; 2 of the 3 days measured differ and agree: the "
          "window's other pointing is named, and its days left out with the rest (the 'most days differ' guard "
          "kept all 7, UNCONFIRMED storm days included)")
    rows = [dict(r) for r in storm3 + E] + [_pt_row("c1", "2026-10-03", "off the calibration", 41.8, 0.91, 0.25, 0.4,
                                                    "ok", 900)]
    sp.pointing_verdicts(rows, ref_day="2026-09-29")
    res = dict(sp.pointing_leave_out(rows), rows=rows, info={})
    check(sp.skip_days_of(res, SimpleNamespace(keep_moved_days=False))
          == ["2026-09-26", "2026-09-27", "2026-09-28", "2026-09-30", "2026-10-01", "2026-10-02"]
          and (sp.pointing_majority_note(res, SimpleNamespace(keep_moved_days=False)) or "").startswith(
              "6 of the window's 8 days are left out"),
          "... with one more matched day: left out, and the warning counts the 6 of 8 days left out")
    check(not sp.pointing_verdicts([dict(r) for r in rows], ref_day="2026-09-29")["other_pointing"],
          "... and no other pointing named: the 2 days that differ do not outnumber the 2 at the reference's")
    # -- review 3, H_T06 (real frames: c1 tilted exactly +0.6 deg on 30 Sep - 2 Oct, read +0.592-0.598 on ~1680
    # inliers, its clear horizon +1.49-1.51 against +0.91; 26-28 Sep UNCONFIRMED, as the sandbox's storm frames
    # did not match its stand-in reference): 3 of the 4 days measured differ, and agree. The 'most days differ'
    # guard kept all 7 days (two pointings and the storm days in one product); now every one of them is left
    # out and the step names the other pointing
    ht = _nuc_rows("c1")
    by(ht)["2026-09-27"].update(features="unmatched", inliers=9, d_azimuth_deg="", d_tilt_deg="", d_roll_deg="")
    for day, t_ in (("2026-09-30", 0.598), ("2026-10-01", 0.597), ("2026-10-02", 0.592)):
        by(ht)[day].update(features="MOVED", inliers=1680, d_azimuth_deg=0.001, d_tilt_deg=t_, d_roll_deg=0.002,
                           horizon_dtilt_deg=0.91 + t_, horizon_offset_px=70.8, turn_sigma=80.0)
    iht = sp.pointing_verdicts(ht, ref_day="2026-09-29")
    loht = sp.pointing_leave_out(ht)
    resht = dict(loht, rows=ht, info=dict(iht, reference_day="2026-09-29"))
    no = SimpleNamespace(keep_moved_days=False)
    wht = sp.other_pointing_warning(resht["info"], loht, "pointing/") or ""
    check([r["verdict"] for r in ht] == ["UNCONFIRMED"] * 3 + ["matches"] + ["DIFFERENT POINTING"] * 3
          and loht["leave_out"] and sp.skip_days_of(resht, no) == [r["date"] for r in ht if r["date"] != "2026-09-29"]
          and iht["other_pointing"]["days"] == ["2026-09-30", "2026-10-01", "2026-10-02"]
          and iht["other_pointing"]["change_deg"] == [0.001, 0.597, 0.002]
          and wht.startswith("3 of the 4 days measured differ from the reference 2026-09-29 and agree with each other "
                             "(+0.00/+0.60/+0.00 deg pan/tilt/roll")
          and "the reference's (2026-09-29) is the smaller group" in wht
          and "the photos (pointing/). To build on the other pointing instead, leave the reference day out: "
              "survey_products.py --window 2026-09-30 2026-10-02 --output-root <another root>" in wht
          and (sp.pointing_majority_note(resht, no) or "").startswith("6 of the window's 7 days are left out"),
          "H_T06 (review 3): a +0.6 deg tilt on 3 of the 4 days measured: all 6 days that differ or are unconfirmed "
          "left out (were all kept), the other pointing named, and the count of days left out")
    _, txt = quiet(sp.print_pointing, "c1", dict(resht, info=dict(resht["info"], warning=wht)), no)
    check("the window holds two pointings" in txt and "LEFT OUT of the detection: 2026-09-26 (UNCONFIRMED)" in txt
          and "KEPT" not in txt,
          "print: the other pointing's warning and the days left out, none KEPT")
    # ... and the build is not complete: its product rests on the reference's pointing, the smaller group, the
    # wrong one if the reference day is the odd one out (review 3, REF_T06: the reference day tilted 0.6 deg)
    hd = Path(d) / "two_pointings"
    (hd / "pointing").mkdir(parents=True)
    (hd / "pointing" / "pointing_c1.json").write_text(json.dumps(resht))
    plan_h = {"out": hd, "date": "2026-09-29", "cams": {"c1": {}}}
    tp = []
    for stamp in (None, {"skip_days": sp.skip_days_of(resht, no), "keep_moved_days": False},
                  {"skip_days": [], "keep_moved_days": True}):
        if stamp is not None:
            sp.write_stamp(hd, "detect_c1", {"x": len(tp)}, extra=stamp)
        tp.append(next((m for m in sp.missing_parts(plan_h, {}) if "two pointings" in m), ""))
    check(tp[0].startswith("c1: the window holds two pointings: 2026-09-30, 2026-10-01, 2026-10-02 agree with each "
                           "other, away from the reference's (2026-09-29); the detection has not run")
          and "the product rests on the reference's, the smaller group; if the reference day is the odd one out"
          in tp[1] and "the product holds both (the other's 2026-09-30, 2026-10-01, 2026-10-02 kept)" in tp[2]
          and "survey_products.py --date 2026-09-29 --window 2026-09-30 2026-10-02 --output-root <another root>" in tp[1],
          "two pointings in the window: a missing part (status partial), whatever the detection did, with the "
          "--window that builds on the other pointing (review 4)")
    # ... that window is the other pointing's days on one side of the reference day, the side with more of them
    check(sp.other_pointing_window({"days": ["2026-09-27", "2026-09-28", "2026-09-30", "2026-10-01", "2026-10-02"]},
                                   "2026-09-29") == ("2026-09-30", "2026-10-02")
          and sp.other_pointing_window({"days": ["2026-09-26", "2026-09-27", "2026-09-28", "2026-09-30"]},
                                       "2026-09-29") == ("2026-09-26", "2026-09-28")
          and sp.other_pointing_window({"days": ["2026-09-30"]}, None) is None,
          "other_pointing_window: the side of the reference day with more of the other pointing's days")
    # -- review 3, turn_full: c1 as on the NUC (26/28 Sep UNCONFIRMED, 27 Sep -1.21 deg of pan) and an exact
    # +0.45 deg roll on 30 Sep and 1 Oct (2286/2103 inliers, the clear horizon +0.69 against +0.24): 3 of the 5
    # days measured differ, but not in one way: no other pointing, every one of them left out (the guard kept
    # all 7, the rolled days gridded: ~0.6 m of DEM at 250-350 m)
    tf = _nuc_rows("c1")
    for day, n_ in (("2026-09-30", 2286), ("2026-10-01", 2103)):
        by(tf)[day].update(features="moved", inliers=n_, d_azimuth_deg=0.002, d_tilt_deg=0.001, d_roll_deg=0.449,
                           horizon_droll_deg=0.69, turn_sigma=60.0)
    itf = sp.pointing_verdicts(tf, ref_day="2026-09-29")
    restf = dict(sp.pointing_leave_out(tf), rows=tf, info=dict(itf, reference_day="2026-09-29"))
    check([r["verdict"] for r in tf] == ["UNCONFIRMED", "DIFFERENT POINTING", "UNCONFIRMED", "matches",
                                          "DIFFERENT POINTING", "DIFFERENT POINTING", "matches"]
          and itf["other_pointing"] is None and sp.other_pointing_warning(restf["info"], restf, "x") is None
          and sp.skip_days_of(restf, no) == ["2026-09-26", "2026-09-27", "2026-09-28", "2026-09-30", "2026-10-01"],
          "turn_full (review 3): days that differ from the reference but not from each other in one way (a -1.21 deg "
          "pan; a +0.45 deg roll): no other pointing named, all left out (were all kept)")
    check(not [m for m in sp.missing_parts(plan_h, {"pointing": {"c1": restf}}) if "two pointings" in m],
          "... and no 'two pointings' missing part")
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
          and "--skip-days" not in note and "--window FIRST LAST" in note
          and sp.pointing_majority_note(res, SimpleNamespace(keep_moved_days=True)) is None,
          "4 of 7 days left out (each a clear jump that its zero-turn features contradict): left out, with the count "
          "and the options this script has (no --skip-days: that is detect_original_view.py's)")
    note2 = sp.pointing_majority_note(res, no, skip=sp.skip_days_of(res, None), ran=False) or ""
    check(note2.startswith("4 of the window's 7 days would be left out of the detection (it has not run with these "
                           "pointing results yet)") and "would rest on the 3 other day(s)" in note2
          and sp.pointing_majority_note(res, no, skip=[]) is None,
          "the README of a build whose detection has not run: the days 'would be left out'; the detection's own "
          "stamp (nothing left out) gives no warning")
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
        sp.pointing_verdicts(rows, ref_day="2025-01-23")       # the calibration day, as check_pointing takes it
        got = {r["date"]: r["verdict"] for r in rows if r["date"] in want}
        check(got == want and all(r["verdict"] == "horizon only" for r in rows if r["date"] not in want)
              and set(sp.pointing_leave_out(rows)["different_days"]) == set(want),
              f"re-aim: {got} (one day on a clear horizon: unconfirmed; the same change on two days, or the features: "
              f"measured; a faint horizon: unconfirmed): left out either way; the reference day next to it is never "
              f"taken in")
    # ... a pan measured next to days whose clear horizon shows tilt and roll back but not the pan ('horizon only':
    # no feature reading): with no reference day between, nothing shows the pan back on any of them (review 4,
    # P10_30_FOG01: a 1.0 deg pan on 30 Sep and 1 Oct, 1 Oct fogged over the land, 'horizon only', kept)
    rows = [dict(r) for r in jan] + [_pt_row("c1", "2025-01-24", "off the calibration", 90.0, 0.16, -4.08, 1.0,
                                             "MOVED", 300, (-15.3, 0.2, -4.0), edge=EDGE_FIXTURE)]
    sp.pointing_verdicts(rows)
    check(all(r["verdict"] == "UNCONFIRMED" and r["next_to"] for r in rows[:6])
          and "the day next to it, 2025-01-24, was measured to differ from the reference (a pan of -15.30 deg, which "
              "no horizon sees): its clear horizon shows tilt and roll back but cannot see the pan" in rows[5]["why"]
          and "it is in a run of days nothing could read, or only the horizon, next to 2025-01-24" in rows[0]["why"]
          and sp.pointing_leave_out(rows)["next_to_days"] == [r["date"] for r in rows[:6]],
          "... with no reference day between: every 'horizon only' day next to the measured -15.3 deg pan, and its "
          "run, UNCONFIRMED (its horizon cannot see a pan)")
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
    # ... and one made before the forcing's waves showed a storm day (they decide whether a day nothing read next
    # to a misread one is left out): stale, and said so
    (pd_ / "forcing").mkdir(exist_ok=True)
    (pd_ / "forcing" / "waves.csv").write_text("time_utc,epoch,wvht_m,dpd_s,source\n"
                                                "x,1737219600,2.4,9.0,buoy\n")       # 18 Jan 2025 17:00 UTC
    sig2 = sp.pointing_sig(pl, "c1", a, photos)
    sp.write_stamp(pd_, "pointing_c1", sig)
    try:
        sp.pointing_skip_days(pl, {}, "c1", a, photos)
        msg = ""
    except sp.StepFailed as exc:
        msg = str(exc)
    check(sig2["storm_days"] == {"2025-01-18": "offshore Hs up to 2.40 m, over 1.5 m (forcing/waves.csv)"}
          and not sig["storm_days"] and "for other storm days (the forcing's waves changed since" in msg,
          "the pointing step's signature holds the window's storm days: a result made before the waves showed one is "
          "not used by the detection, and the reason says so")
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
    got = []
    for rud in ("+0.000/+0.400/+0.000", "+0.000/+0.000/-0.350", "+0.450/+0.000/+0.000"):
        by(rows)["2026-10-05"].update(runner_up_d=rud)
        sp.pointing_verdicts(rows)
        got.append(by(rows)["2026-10-05"]["verdict"])
    check(got == ["UNCONFIRMED", "UNCONFIRMED", "small change"],
          f"an ambiguous match's runner-up is judged as a turned match: 0.40 deg of tilt or 0.35 of roll "
          f"UNCONFIRMED, a 0.45 deg pan a small change ({got})")
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


def test_oct_pointing_pattern(d):
    print("the Oct 2024 fixture's pointing pattern: two planted pan days against zero-turn days; the calibration's own "
          "day as the reference; the README's why for unread days")
    by = lambda rows: {r["date"]: r for r in rows}          # noqa: E731
    faint = (8.3, 0.4, 1.0, 200)                            # prominence 8.3: not a clear horizon
    # review r1 (runner-builds, review-blind): c2 of the synthetic Oct 2024 fixture at HEAD 123224b. 20-21 Oct carry
    # a planted +1.0 deg pan (181/206 inliers); every other day is rendered at the reference's pointing and matches
    # it at 0.00 deg within its own noise (no reading: blind_match), each horizon faint
    def oct_rows():
        rows = []
        for day, feat, n, dd, ts in (("2024-10-20", "MOVED", 181, (1.0, 0.0, 0.0), 60.0),
                                     ("2024-10-21", "MOVED", 206, (1.0, 0.0, 0.0), 70.0),
                                     ("2024-10-22", "ok", 263, (0.0, 0.0, 0.0), 0.5),
                                     ("2024-10-23", "ok", 1500, (0.0, 0.0, 0.0), 0.3),
                                     ("2024-10-24", "ok", 234, (0.0, 0.0, 0.01), 1.9),
                                     ("2024-10-25", "ok", 236, (0.0, 0.0, 0.0), 1.6),
                                     ("2024-10-26", "ok", 222, (0.0, 0.0, 0.01), 2.1)):
            r = _pt_row("c2", day, "matches", 9.0, 0.14, -0.01, 0.6, feat, n, dd, edge=faint)
            r["turn_sigma"] = ts
            rows.append(r)
        return rows
    rows = oct_rows()
    import random
    shuffled = [dict(r) for r in rows]
    random.Random(7).shuffle(shuffled)
    i = sp.pointing_verdicts(rows, ref_day="2024-10-23")
    i_s = sp.pointing_verdicts(shuffled, ref_day="2024-10-23")
    v, vs = by(rows), by(shuffled)
    check([v[x]["verdict"] for x in sorted(v)] == ["DIFFERENT POINTING", "DIFFERENT POINTING", "UNCONFIRMED",
                                                    "matches", "unchecked", "unchecked", "unchecked"]
          and v["2024-10-22"]["next_to"] == "2024-10-21"
          and all(vs[x]["verdict"] == v[x]["verdict"] and vs[x]["why"] == v[x]["why"] for x in v) and i_s == i,
          "Oct fixture c2: 20-21 Oct DIFFERENT POINTING (the planted pan), 22 Oct UNCONFIRMED next to 21 Oct (its "
          "zero-turn match is no reading: days lost, not wrong ones), 24-26 Oct unchecked, kept; the same in any "
          "row order")
    check(i["other_pointing"] is None,
          "... and NO other pointing: the 2 moved days do not outnumber the reference day with the 4 days whose "
          "features matched it at zero turn within their noise (was: 'the reference's is the smaller group', PARTIAL)")
    # a shorter window: 20-23 Oct and 24 Oct unmatched. 22 Oct, UNCONFIRMED only for the day next to it, read
    # nothing but zero turn and counts: 2 moved against 1 + 1, no other pointing
    rows = oct_rows()[:4] + [dict(oct_rows()[4], features="unmatched", inliers=12, turn_sigma="")]
    i = sp.pointing_verdicts(rows, ref_day="2024-10-23")
    check(i["other_pointing"] is None and by(rows)["2024-10-22"]["verdict"] == "UNCONFIRMED",
          "2 moved days against the reference day and 22 Oct (UNCONFIRMED next to 21 Oct, its features at zero turn "
          "within their noise): not outnumbered, no other pointing")
    # 22 Oct unmatched too: 2 moved against the reference day alone, the other pointing is named
    rows = (oct_rows()[:2] + [dict(oct_rows()[2], features="unmatched", inliers=8, turn_sigma="")] + oct_rows()[3:4]
            + [dict(oct_rows()[4], features="unmatched", inliers=12, turn_sigma="")])
    i = sp.pointing_verdicts(rows, ref_day="2024-10-23")
    op = i["other_pointing"] or {}
    lo = sp.pointing_leave_out(rows)
    for cal, want_cal in (("CACO03_c2_20241023_EO.yaml", True), ("CACO05_c2_20251113_EO-CV.yaml", False)):
        info = dict(i, reference_day="2024-10-23", calibration=cal, reference_frame="1729702800.c2.timex.jpg")
        w = sp.other_pointing_warning(info, lo, "pointing/", "2024-10-23") or ""
        check(sp.reference_is_calibration_day(info) == want_cal and op.get("days") == ["2024-10-20", "2024-10-21"]
              and ("the calibration's own day (CACO03_c2_20241023_EO.yaml), so its pointing is the calibration's" in w
                   and "leave the reference day out" not in w and "is the smaller group" not in w
                   if want_cal else
                   "the reference's (2024-10-23) is the smaller group" in w and "leave the reference day out" in w),
              f"other pointing against a reference that {'IS' if want_cal else 'is NOT'} the calibration's own solve "
              f"day: " + ("the warning says the product rests on the calibration's pointing, advises no rebuild"
                          if want_cal else "the reference may be the odd one out, the rebuild without it is given"))
        hd = Path(d) / f"oct_two_pointings_{want_cal}"
        (hd / "pointing").mkdir(parents=True)
        res = dict(lo, rows=rows, info=dict(info, warning=w))
        (hd / "pointing" / "pointing_c2.json").write_text(json.dumps(res))
        plan_h = {"out": hd, "date": "2024-10-23", "cams": {"c2": {}}}
        tp = []
        for stamp in (None, {"skip_days": sp.skip_days_of(res, SimpleNamespace(keep_moved_days=False)),
                             "keep_moved_days": False}, {"skip_days": [], "keep_moved_days": True}):
            if stamp is not None:
                sp.write_stamp(hd, "detect_c2", {"x": len(tp)}, extra=stamp)
            tp.append(next((m for m in sp.missing_parts(plan_h, {}) if "two pointings" in m), ""))
        if want_cal:
            check(tp[0] and not tp[1] and "the product holds both" in tp[2] and "leave the reference day out" not in
                  "".join(tp),
                  "... on the calibration's own day: NOT a missing part once the detection left the other pointing out "
                  "(status complete); still one while the detection has not run, or with --keep-moved-days (both "
                  "pointings in the product)")
        else:
            check(all(tp) and "the smaller group" in tp[1],
                  "... on another reference: a missing part whatever the detection did (as before)")
    # a 'small change' day kept: the README gives its measured turn and the DEM effect per range band (review r1)
    from georectify import load_intrinsics, load_extrinsics, pixel_to_ground
    cal = HERE / "calibration"
    io_f, eo_f = cal / "CACO05_c2_20240801_IO.yaml", cal / "CACO03_c2_20250123_EO.yaml"
    io_, e_ = load_intrinsics(io_f), load_extrinsics(eo_f)
    U, V = np.meshgrid(np.linspace(300, 2100, 40), np.linspace(280, 1100, 40))
    X, Y = pixel_to_ground(U.ravel(), V.ravel(), np.full(U.size, 0.5), io_, e_)
    cell = 2.0
    xll, yll = float(np.floor(X.min())) - 4, float(np.floor(Y.min())) - 4
    nc, nr = int((X.max() - xll) / cell) + 4, int((Y.max() - yll) / cell) + 4
    yy = yll + (nr - 0.5 - np.arange(nr))[:, None] * cell
    g = 0.5 + 0.05 * (yy - yy.mean()) + 0.0 * np.arange(nc)[None, :]
    plan_s = {"date": "2025-01-23", "era": "adcp", "out": Path(d) / "small_change" / "2025-01-23",
              "cams": {"c2": {"eo_path": eo_f, "io_path": io_f}}}
    stem = sp.dem_stem(plan_s)
    stem.parent.mkdir(parents=True)
    Path(f"{stem}_dem.asc").write_text(f"ncols {nc}\nnrows {nr}\nxllcorner {xll}\nyllcorner {yll}\ncellsize {cell}\n"
                                       f"NODATA_value -9999\n" + "\n".join(" ".join(f"{v:.3f}" for v in row)
                                                                          for row in g) + "\n")
    sc = [_pt_row("c2", "2025-01-22", "matches", 8.0, 0.15, 0.0, 0.4, "moved", 180, (0.32, 0.12, -0.05),
                  edge=EDGE_CLEAR["c2"]), _pt_row("c2", "2025-01-23", features="ok", inliers=900, edge=EDGE_CLEAR["c2"])]
    for r in sc:
        r["turn_sigma"] = 40.0
    sp.pointing_verdicts(sc, ref_day="2025-01-23")
    txt = sp.small_change_caveat(plan_s, "c2", sc)
    bands = sp.tilt_sensitivity(plan_s, "c2", 0.12, -0.05) or []
    check(sc[0]["verdict"] == "small change" and txt.startswith("c2: kept with a small change of pointing")
          and "2025-01-22: +0.32/+0.12/-0.05 deg pan/tilt/roll on 180 features" in txt and len(bands) >= 2
          and all(f"~{b['dem_abs_error_m']:.2f} m at {b['range_m']} m" in txt for b in bands)
          and not sp.small_change_caveat(plan_s, "c2", sc, skip=["2025-01-22"]),
          f"a kept 'small change' day: the README gives its measured turn and its DEM effect per range band "
          f"({len(bands)} bands: " + ", ".join(f"{b['dem_abs_error_m']:.2f} m at {b['range_m']} m" for b in bands)
          + "); nothing for a day left out")
    # where the ~0.1 deg of a match rests: the inliers' spread is recorded, and a deciding match bunched in a small
    # part of the frame says so (review r1: a 300 px patch reads roll up to 0.17 deg; every tested match spread
    # 180 px or more)
    import pointing_check as pc_
    box = np.array([[x, y] for x in np.linspace(1000, 1300, 31) for y in np.linspace(500, 600, 11)])
    spr = pc_.inlier_spread(box)
    fr = sp.feature_record({"d_azimuth": 0.8, "d_tilt": 0.02, "d_roll": -0.01, "inliers": 341, "turn_sigma": 60.0,
                            "runner_up_inliers": 3, "inlier_spread_px": spr}, Path("x.jpg"))
    brow = [dict(_pt_row("c2", "2024-10-21", "matches", 9.0, 0.14, -0.01, 0.6, edge=faint), **fr),
            dict(oct_rows()[3])]
    sp.pointing_verdicts(brow, ref_day="2024-10-23")
    want = [float(np.std(box[:, 1], ddof=1)), float(np.std(box[:, 0], ddof=1))]      # an axis-aligned box
    check(np.allclose(spr, want, atol=0.01) and fr["inlier_spread_px"] == "32/90"
          and brow[0]["verdict"] == "DIFFERENT POINTING" and sp.bunched(brow[0])
          and "its 341 inliers are bunched (spread 32 x 90 px; every match the ~0.1 deg rests on spread 180 px or more)"
          in brow[0]["why"] and "spread over the frame is good to ~0.1 deg on the views tested" in brow[0]["why"]
          and not sp.bunched(dict(brow[0], inlier_spread_px="250/520")) and not sp.bunched({"inliers": 50}),
          "a match's inlier spread is recorded (inlier_spread_px, minor/major px); a deciding match bunched under "
          "180 px says so in its why (the verdict unchanged), and the claim reads 'spread over the frame ... on the "
          "views tested'")
    # the README's why for days nothing read: matched within noise is said so, never 'too few matched features'
    rows = oct_rows()
    sp.pointing_verdicts(rows, ref_day="2024-10-23")
    unchk = [r["date"] for r in rows if r["verdict"] == "unchecked"]
    t = sp.no_feature_reading(unchk, rows)
    t2 = sp.no_feature_reading(unchk + ["2024-10-27"], rows + [_pt_row("c2", "2024-10-27", features="unmatched",
                                                                       inliers=7, edge=faint)])
    check(t.startswith("on each the features matched the reference, but within their own noise") and "too few" not in t
          and "a camera that did not move at all" in t
          and "on 2024-10-24, 2024-10-25, 2024-10-26 the features matched" in t2
          and "on 2024-10-27 too few matched features against the reference frame (under 40)" in t2,
          "the README's why for an unchecked day: matched within their own noise where 40+ features matched, too few "
          "matched only where they did not")


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
    # the DEM's day test, in the same sentence (review 2: 'with the camera-day is off the rest ...')
    frames(nuc)
    inf_ = out / "dem" / "2026-09-29_info.json"
    j_ = json.loads(inf_.read_text())
    j_["filters"].update(camera_days_rejected=["c1 2026-10-01"])
    inf_.write_text(json.dumps(j_))
    w2 = next(c for c in sp.collect_caveats(plan, a, {}, forcing, {"frames_per_day": {}}) if c.startswith("HIGH WAVES"))
    check("it left out 4 frame(s) of camera-days off the rest by > 0.15 m (the DEM's day test): 4 frame(s) after "
          "them (2026-10-01)" in w2 and "it left out 6 frame(s) with offshore Hs > 1.5 m" in w2,
          "HIGH WAVES: the frames of a camera-day the DEM's day test rejected, in a sentence that reads")
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
    # ... a pointing-only rerun over a deployed build (review 4: its detection left c2 26/27 Sep out on the old
    # rule's horizon jump; these results keep them): the README says what the detection on disk did, never
    # '--keep-moved-days' (which leaves out no day), never 'kept' for days the detection on disk left out, and
    # the STORM caveat names the earlier pointing result
    sp.write_stamp(out, "detect_c2", {"x": 8}, extra={"skip_days": ["2026-09-26", "2026-09-27"], "setup_coef": 0.037})
    cav = sp.collect_caveats(plan, a, {}, forcing, frames(nuc))
    c2_ = [c for c in cav if c.startswith("c2: ")]
    storm = next(c for c in cav if c.startswith("STORM:"))
    check(any(c.startswith("c2: the detection left out 2026-09-26, 2026-09-27, where these pointing results leave out "
                           "no day: it ran with other pointing results (the pointing step was run again since); run the "
                           "detection again") for c in c2_)
          and not any("--keep-moved-days" in c for c in c2_) and not any(c.startswith("c2: kept, although") for c in c2_)
          and any(c.startswith("c2: these pointing results keep 2026-09-26, 2026-09-27, although the sea-horizon check "
                               "alone would leave them out") and "the detection on disk left out 2026-09-26, 2026-09-27, "
                               "so they are kept once it runs again" in c for c in c2_)
          and "c2 2026-09-26 (matches), c2 2026-09-27 (matches) left out by an earlier pointing result (the pointing "
              "step was run again since" in storm,
          "README after a pointing-only rerun over a deployed build that left c2 26/27 Sep out: the detection's own "
          "days, no '--keep-moved-days', no 'kept' for them, the STORM caveat names the earlier pointing result")
    sp.write_stamp(out, "detect_c2", {"x": 3}, extra={"skip_days": []})
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
    sp.write_stamp(out, "detect_c1", {"x": 4}, extra={"skip_days": [], "keep_moved_days": True})
    cav = sp.collect_caveats(plan, cfg_args(keep_moved_days=True), {}, forcing, frames(new))
    kept = next((c for c in cav if c.startswith("c1: days with a different or unconfirmed pointing KEPT in the "
                                                 "detection (--keep-moved-days): 2026-09-26 (UNCONFIRMED: ")), "")
    check(kept and "left out" not in kept and "2026-09-27 (DIFFERENT POINTING: 48 matched features" in kept,
          "README with --keep-moved-days: the days KEPT, with whys that do not say 'left out'")
    # ... and a later run without the option (--steps compare) that rewrites the README says what the
    # detection did, from its stamp, not what this run's options would do (review 2)
    cav2 = sp.collect_caveats(plan, cfg_args(), {}, forcing, frames(new))
    check(any(c.startswith("c1: days with a different or unconfirmed pointing KEPT in the detection "
                           "(--keep-moved-days): 2026-09-26") for c in cav2)
          and not any("most days differ from the reference: see the warning" in c for c in cav2)
          and not any(c.startswith("c1: days left out of the detection") for c in cav2),
          "README rewritten by a later run without --keep-moved-days: still 'KEPT (--keep-moved-days)', from the "
          "detection's stamp")
    # an older stamp, as the deployed code writes it (no option recorded), that left out no day, under results
    # that leave days out (review 3: the README inferred '--keep-moved-days', which was never passed, while the
    # step said 'LEFT OUT of the detection'): the option is never inferred; either reason is named
    sp.write_stamp(out, "detect_c1", {"x": 5}, extra={"skip_days": [], "setup_coef": 0.037})
    cav2 = [c for c in sp.collect_caveats(plan, cfg_args(), {}, forcing, frames(new)) if c.startswith("c1: ")]
    check(any(c.startswith("c1: the detection left out no day, where these pointing results leave out 2026-09-26, "
                           "2026-09-27, 2026-09-28 by default: it ran with other pointing results (the pointing step was "
                           "run again since) or with --keep-moved-days (its record, of an older build, does not say "
                           "which); run the detection again") for c in cav2)
          and not any("KEPT in the detection" in c for c in cav2)
          and any(c.startswith("c1: days with a different or unconfirmed pointing: 2026-09-26 (UNCONFIRMED: ")
                  and "The detection has not run with these pointing results" in c for c in cav2),
          "an older detection stamp without the option that left out no day, under results that leave days out: "
          "'other pointing results or --keep-moved-days', never 'KEPT (--keep-moved-days)'")
    check(sp.kept_reason({"leave_out": True}, cfg_args(), {"skip_days": []}).startswith(
              "--keep-moved-days, or other pointing results: the detection's record, of an older build, does not say")
          and sp.kept_reason({"leave_out": False}, cfg_args(), {"skip_days": []}).startswith("most days differ")
          and sp.kept_reason({"leave_out": True}, cfg_args(), {"keep_moved_days": True}) == "--keep-moved-days",
          "kept_reason never infers --keep-moved-days from the pointing result")
    # ... and with no detection at all (a pointing-only or partial build, or a failed detection): neither
    # KEPT nor 'left out of the detection' (review 2: the README said both, for the same days)
    sp.drop_stamp(out, "detect_c1")
    cav3 = sp.collect_caveats(plan, cfg_args(), {}, forcing, frames(new))
    c1_ = [c for c in cav3 if c.startswith("c1: ")]
    nd = next((c for c in c1_ if c.startswith("c1: days with a different or unconfirmed pointing: 2026-09-26 "
                                               "(UNCONFIRMED: ")), "")
    check(nd and "The detection has not run with these pointing results" in nd and "leaves such days out by "
          "default" in nd and not any("KEPT in the detection" in c or c.startswith("c1: days left out of the "
                                                                                   "detection") for c in c1_)
          and not any("are left out of the detection" in c for c in c1_),
          "README with no detection stamp: the days and their whys, 'the detection has not run', never KEPT or left "
          "out")
    # ... and a detection made with OTHER pointing results (the pointing step run again on its own since, as
    # after a change of its rule): what it left out is said as such, and these results' days are the ones it
    # would leave out, never 'KEPT' (it left out 30 Sep only; these results leave out 26-28 Sep)
    sp.write_stamp(out, "detect_c1", {"x": 6}, extra={"skip_days": ["2026-09-30"], "keep_moved_days": False})
    c1_ = [c for c in sp.collect_caveats(plan, cfg_args(), {}, forcing, frames(new)) if c.startswith("c1: ")]
    check(any(c.startswith("c1: the detection left out 2026-09-30, where these pointing results leave out 2026-09-26, "
                           "2026-09-27, 2026-09-28: it ran with other pointing results (the pointing step was run "
                           "again since); run the detection again") and "--steps detect" in c for c in c1_)
          and any(c.startswith("c1: days with a different or unconfirmed pointing: 2026-09-26 (UNCONFIRMED: ")
                  and "The detection has not run with these pointing results: it leaves such days out" in c for c in c1_)
          and not any("KEPT in the detection" in c or c.startswith("c1: days left out of the detection") for c in c1_),
          "README after a pointing-only rerun over an older detection: the detection's own left-out days said as "
          "such, these results' days 'not run with these', never KEPT")
    sp.write_stamp(out, "detect_c1", {"x": 7}, extra={"skip_days": ["2026-09-26", "2026-09-27", "2026-09-28"],
                                                      "keep_moved_days": False})
    c1_ = [c for c in sp.collect_caveats(plan, cfg_args(), {}, forcing, frames(new)) if c.startswith("c1: ")]
    check(any(c.startswith("c1: days left out of the detection for their pointing: 2026-09-26 (UNCONFIRMED: ")
              for c in c1_) and not any("other pointing results" in c or "has not run" in c for c in c1_),
          "... a detection made with these results: its left-out days, nothing about other results")
    sp.drop_stamp(out, "detect_c1")
    # WHAT THIS IS: the DEM's own frames and days next to the waterlines' (review 3: '7 days' where the DEM
    # gridded 4, its wave filter having taken the storm's frames)
    check(sp.dem_days_text(plan, frames(nuc)["frames_per_day"])
          == "; the DEM grids 28 of their 34 frames, on 2026-09-29 .. 2026-10-02 (4 days)"
          and sp.dem_days_text(plan, frames(after)["frames_per_day"]) == ""
          and sp.dem_days_text(plan, frames(nuc, dem=False)["frames_per_day"]) == "",
          "WHAT THIS IS: the DEM's 28 of 34 frames on 4 days, when it gridded fewer than the waterlines hold; nothing "
          "added when it gridded them all, or without its record")
    # 'horizon only' days whose features matched within their own noise (the Mar 2025 fixture's c2: 51-59 inliers,
    # turn_sigma 0.7-3.4): the README says why the pan was not measured (review 3)
    mar = [dict(_pt_row("c2", f"2026-09-{dd}", "matches", 8.0, 0.15, 0.0, 0.4, "ok", n_, (0.001, 0.002, -0.004),
                        edge=EDGE_FIXTURE), turn_sigma=ts_)
           for dd, n_, ts_ in (("26", 55, 3.0), ("27", 51, 2.4), ("28", 59, 1.3), ("29", 57, 0.7))]
    mar.append(_pt_row("c2", "2026-09-30", "matches", 8.0, 0.15, 0.0, 0.4, "unmatched", 12, edge=EDGE_FIXTURE))
    im = sp.pointing_verdicts(mar, ref_day="2026-09-29")
    (out / "pointing" / "pointing_c2.json").write_text(json.dumps(dict(sp.pointing_leave_out(mar), rows=mar, info=im)))
    hz = next((c for c in sp.collect_caveats(plan, a, {}, forcing, frames(new)) if c.startswith("c2: pan not measured")),
              "")
    check(hz.startswith("c2: pan not measured on 2026-09-26, 2026-09-27, 2026-09-28, 2026-09-30 (horizon only: on "
                        "2026-09-26, 2026-09-27, 2026-09-28 the features matched the reference, but within their own "
                        "noise (turn_sigma under 4") and "on 2026-09-30 too few matched features" in hz,
          "README: 'horizon only' days whose features matched within their own noise say so (the unmatched day: too "
          "few matched)")
    mar[-1].update(features="ok", inliers=52, d_azimuth_deg=0.0, d_tilt_deg=0.001, d_roll_deg=0.0, turn_sigma=1.9)
    im = sp.pointing_verdicts(mar, ref_day="2026-09-29")
    (out / "pointing" / "pointing_c2.json").write_text(json.dumps(dict(sp.pointing_leave_out(mar), rows=mar, info=im)))
    hz = next((c for c in sp.collect_caveats(plan, a, {}, forcing, frames(new)) if c.startswith("c2: pan not measured")),
              "")
    check(hz.startswith("c2: pan not measured on 2026-09-26, 2026-09-27, 2026-09-28, 2026-09-30 (horizon only: on each "
                        "the features matched the reference, but within their own noise"),
          "... and when every such day's match was within its noise: 'on each', not the list twice")


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


def test_day_test_record(d):
    print("the DEM's day test record reaches the README, the caveats and the per-frame reasons (a variant's DEM too)")
    import dem_from_contours as dfc
    import dem_figure
    out = Path(d) / "daytest_rec"
    for sub in ("dem", "waterlines", "pointing", "logs"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    # a record as dem_from_contours writes it: c1 19 Jan left out 1st, then refused take-back (within the
    # limit of the kept days, but taking it back would put c2 20 Jan over); c2 20 Jan left out 2nd; c1 22 Jan
    # not judged (over the limit in the first pass, its DEM one other day); c2's offsets run one way with time,
    # the trend resting on c2 23 Jan, the one day left out at its end
    lim = 0.15
    days = {"c1 2025-01-19": 0.21, "c2 2025-01-20": -0.30, "c1 2025-01-22": 0.40, "c1 2025-01-21": 0.02,
            "c2 2025-01-21": 0.01, "c2 2025-01-22": -0.05}
    dt = {"first": {k: (v, 300) for k, v in days.items()},
          "final": dict({k: (v, 300) for k, v in days.items()}, **{"c1 2025-01-19": (0.08, 300)}),
          "steps": [("c1 2025-01-19", 0.21, 300), ("c2 2025-01-20", -0.30, 300)],
          "not_judged": [{"day": "c1 2025-01-22", "offset": 0.40, "first": 0.40, "left_out": True,
                          "why": "its DEM is one other day (c2 2025-01-22) that the rest do not confirm"}],
          "refused": [("c1 2025-01-19", 0.08, ["c2 2025-01-20"])], "taken_back": [],
          "rejected": ["c1 2025-01-19", "c1 2025-01-22", "c2 2025-01-20"], "untested": [], "over_kept": [],
          "too_few": False, "testable": 6}
    trend = {"camera": "c2", "days": 5, "spearman_rho": -1.0, "chance": 1 / 120.0, "range_m": 0.31,
             "first": ("2025-01-19", 0.012), "last": ("2025-01-23", -0.298),
             "without": {"day": "c2 2025-01-23", "days": 4, "spearman_rho": -1.0, "chance": 1 / 24.0,
                         "range_m": 0.09, "alone": False}}
    rec = dfc.day_test_record(dt, lim, [trend])
    filt = dict(max_hs=1.5, rough_frames=0, max_day_offset=lim, camera_days_rejected=list(dt["rejected"]), **rec)
    how = {r["day"]: r["how"] for r in filt["day_test"]["rejected"]}
    check(how == {"c1 2025-01-19": "refused take-back", "c1 2025-01-22": "not judged", "c2 2025-01-20": "left out"}
          and filt["day_test"]["warning"] and filt["day_test"]["trend_warnings"][0]["without"]["alone"] is False,
          f"the fixture is dem_from_contours' own record: one day of each 'how', its WARNING, a trend ({how})")
    # one wording: the DEM page's
    check(dem_figure.day_test_left_out(filt) == "3 camera-days left out by the day test (±0.15 m; 1 within it but "
                                                "at odds with a kept day)"
          and dem_figure.day_test_left_out(filt, "+/-").startswith("3 camera-days left out by the day test (+/-0.15")
          and dem_figure.day_test_left_out(dict(filt, day_test=None)) == "3 camera-days off the rest by > 0.15 m",
          "dem_figure.day_test_left_out: 'off the rest' only without a day left out within the limit")
    r_ = {k: dem_figure.day_test_reason(filt, k) for k in how}
    check(r_["c2 2025-01-20"] == "camera-days off the rest by > 0.15 m"
          and r_["c1 2025-01-19"].startswith("camera-days within +/-0.15 m of the kept days")
          and "off the rest" not in r_["c1 2025-01-19"] and "could not judge" in r_["c1 2025-01-22"]
          and dem_figure.day_test_reason(dict(filt, day_test=None), "c1 2025-01-19") ==
          "camera-days off the rest by > 0.15 m",
          "dem_figure.day_test_reason: each day by its 'how'; a record of the one-pass test: 'off the rest'")
    # the README's lines under the DEM (WHAT WAS USED)
    info = {"window_start": "2025-01-18", "window_end": "2025-01-23", "frames": {"c1": 10, "c2": 10}, "days": 6,
            "filters": filt}
    rl = "\n".join(sp.dem_built_from_lines(info))
    check("3 camera-days left out by the day test (+/-0.15 m; 1 within it but at odds with a kept day)" in rl
          and "off the rest by > 0.15 m" not in rl.split("\n")[0]
          and "day test: c1 2025-01-19: refused take-back, +0.080 m from the DEM of the kept days (+0.210 m when it "
              "was left out); why: left out 1st; within the limit of the kept days, but taking it back would put "
              "c2 2025-01-20 over the limit" in rl
          and "day test: c2 2025-01-20: left out, -0.300 m from the DEM of the kept days; why: left out 2nd" in rl
          and "day test: c1 2025-01-22: not judged" in rl
          and "day test WARNING: the day test could not judge every camera-day -- left out without being judged: "
              "c1 2025-01-22." in rl
          and "day test WARNING: c2's camera-day offsets run one way with time over 5 days (Spearman rho -1.00, "
              "which days in random order reach in that direction 1 time in 120 (either way 1 in 60)" in rl
          and "The method alone drifts so: on a static synthetic beach" in rl
          and "it rests on c2 2025-01-23, the one day left out, at an end: the other 4 days give rho -1.00" in rl
          and all(ord(ch) < 128 for ch in rl),
          "README: every camera-day left out with how, its offsets and why; the WARNING; the trend 'rests on' its "
          "trimmed end day (ASCII)")
    alone = json.loads(json.dumps(filt))
    alone["day_test"]["trend_warnings"][0]["without"]["alone"] = True
    tw = dem_figure.day_test_warnings(alone)[-1]
    check("without c2 2025-01-23, the one day left out, at an end, the other 4 days give rho -1.00" in tw
          and "they pass the rule on their own" in tw and "rests on" not in tw,
          "a trend whose other days pass the rule alone does not 'rest on' the trimmed day")
    # a variant's DEM (B: dem/<date>_gcpfree_*): its record is the one read (was dem/<date>_info.json only)
    plan = {"era": "adcp", "date": "2025-01-23", "first": "2025-01-18", "last": "2025-01-23", "out": out,
            "surveys": [], "cams": {c: {"eo_file": f"CACO03_{c}_20250123_EO.yaml", "first": "2025-01-18",
                                        "last": "2025-01-23", "variant": "gcpfree"} for c in ("c1", "c2")}}
    (out / "dem" / "2025-01-23_gcpfree_info.json").write_text(json.dumps(info))
    with open(out / "waterlines" / "contour_points_ground.csv", "w") as f:
        f.write("source_file,camera,capture_time_utc,capture_epoch,offshore_hs_m\n")
        for cam, day in (("c1", 19), ("c2", 20), ("c1", 21), ("c1", 22)):
            t = datetime(2025, 1, day, 16, tzinfo=timezone.utc)
            f.write(f"{photo_name(t, 'CACO03', cam)},{cam},{t:%Y-%m-%dT%H:%M:%SZ},{int(t.timestamp())},0.8\n")
    old_wf = sp.waterline_file
    sp.waterline_file = lambda plan_, quiet_=False: out / "waterlines" / "contour_points_ground.csv"
    try:
        dg = sp.dem_gridded_frames(plan)
    finally:
        sp.waterline_file = old_wf
    lw = {cam: sorted({v[2] for v in fr.values()}) for cam, fr in ((dg or {}).get("left_out") or {}).items()}
    check(dg is not None and lw.get("c1") == sorted([
        "camera-days within +/-0.15 m of the kept days, left out as taking them back would put a kept day over "
        "(the DEM's day test)", r_["c1 2025-01-22"] + " (the DEM's day test)"])
          and lw.get("c2") == ["camera-days off the rest by > 0.15 m (the DEM's day test)"]
          and len(dg["gridded"].get("c1") or {}) == 1,
          f"dem_gridded_frames reads the variant DEM's record; each left-out frame by its day's own reason ({lw})")
    forcing = {"waves": {"setup_in_window": {"median_m": 0.29, "daytime_median_m": 0.28}}}
    cav = [c for c in sp.collect_caveats(plan, SimpleNamespace(setup_coef=sp.SETUP_COEF), {}, forcing, {})
           if c.startswith("DEM day test")]
    check(len(cav) == 2 and cav[0].startswith("DEM day test WARNING: the day test could not judge every camera-day")
          and "it rests on c2 2025-01-23" in cav[1] and all("WHAT WAS USED" in c for c in cav),
          "KNOWN CAVEATS: the day test's WARNING and its trend WARNING, from the (variant) DEM's own record")
    (out / "dem" / "2025-01-23_gcpfree_info.json").write_text(json.dumps(dict(info, filters=dict(
        filt, camera_days_rejected=[], day_test=dict(filt["day_test"], warning=None, rejected=[], trend_warnings=[
            dict(filt["day_test"]["trend_warnings"][0], without=None)])))))
    cav = [c for c in sp.collect_caveats(plan, SimpleNamespace(setup_coef=sp.SETUP_COEF), {}, forcing, {})
           if c.startswith("DEM day test")]
    check(len(cav) == 1 and "c2's camera-day offsets run one way" in cav[0] and "rests on" not in cav[0]
          and cav[0].endswith("(it left out no camera-day)"),
          "a trend WARNING with no day left out (the Mar 2025 fixture's case): the caveat says none was left out")
    (out / "dem" / "2025-01-23_gcpfree_info.json").write_text(json.dumps(dict(info, filters=dict(
        filt, day_test=dict(filt["day_test"], warning=None, trend_warnings=[])))))
    check(not [c for c in sp.collect_caveats(plan, SimpleNamespace(setup_coef=sp.SETUP_COEF), {}, forcing, {})
               if c.startswith("DEM day test")], "no WARNING in the record: no day-test caveat")


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
    check(m["turn_sigma"] is not None and m["turn_sigma"] < sp.TURN_SIGMA_MIN
          and fr["turn_sigma"] == round(m["turn_sigma"], 1),
          f"... and its turn is within its own noise: turn_sigma {m['turn_sigma']:.1f} (under {sp.TURN_SIGMA_MIN:g}): "
          f"dirt moves with the camera and matches at exactly zero")
    # the scene itself, the camera turned by a small 0.03 deg of pan (no dirt): 'ok' (0.1 deg or less), but
    # the turn is many times its standard error, which no dirt consensus gives
    eo3 = eo.copy()
    eo3[3:6] += np.radians([0.03, 0.0, 0.0])
    tex = scene(40)
    ref_img = wd / "ref_small.jpg"
    cv2.imwrite(str(ref_img), tex)
    ref_img.with_suffix(".json").write_text(json.dumps({"eo": eo.tolist()}))
    ms = pc.measure(_rotated(tex, io, eo, eo3), [pc.Reference(ref_img)], io, min_inliers=40, second=True)
    frs = sp.feature_record(ms, Path("x.jpg"))
    check(frs["features"] == "ok" and abs(ms["d_azimuth"] - 0.03) < 0.01 and ms["turn_sigma"] > 3 * sp.TURN_SIGMA_MIN,
          f"the scene turned 0.03 deg ({ms['inliers']} inliers, read {ms['d_azimuth']:+.3f}): 'ok', turn_sigma "
          f"{ms['turn_sigma']:.0f}, far above {sp.TURN_SIGMA_MIN:g}: the scene matched")
    rows = [_pt_row("c2", f"2026-10-0{k}", "off the calibration", 33.0, 0.55, -0.67, 0.2, "ok", 1400) for k in range(1, 5)]
    r = _pt_row("c2", "2026-10-05", "off the calibration", 66.0, 1.15, -0.67, 0.2, edge=None)
    r.update({k: v for k, v in fr.items() if k != "feature_frame"})
    sp.pointing_verdicts(rows + [r])
    check(r["verdict"] == "UNCONFIRMED" and "UNCONFIRMED" in sp.LEFT_OUT_VERDICTS,
          "... so against a clear horizon that moved 0.6 deg it decides nothing: UNCONFIRMED, left out")
    r.update(horizon_prominence=6.0, horizon_rival=0.7, horizon_sky_above=0.0, horizon_to_limit_px=30)
    sp.pointing_verdicts(rows + [r])
    check(r["verdict"] == "UNCONFIRMED" and "no further from no turn than their own noise" in r["why"],
          "... and against a horizon that is NOT clear (a storm's foam edge) too: UNCONFIRMED (the first version of the "
          "rule kept it: 'matches')")


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


def summary_product(root, date, variant=None, name="lidar"):
    """A BUILT product as --summary reads it back: the stamps' chain (forcing -> merge -> dem -> compare), the DEM's
    record dem/<stem>_info.json (the variant's own stem, dem_stem) with the settings the dem stamp holds, the
    comparison's json and provenance.json naming each camera's calibration variant ('table' = none)."""
    out = Path(root) / date
    (out / "waterlines").mkdir(parents=True)
    (out / "dem").mkdir()
    (out / "compare").mkdir()
    (out / "waterlines" / "contour_points_ground.csv").write_text("camera,x\n")
    sp.write_stamp(out, "forcing", {"forcing": "test"})
    sp.write_stamp(out, "merge", {"after": {"forcing": sp.token(out, "forcing")}})
    settings = {"cell_m": 1.0, "min_points": 3, "max_spread_m": 0.5, "max_hs_m": 1.5, "max_day_offset_m": 0.15}
    sp.write_stamp(out, "dem", {"after": {"merge": sp.token(out, "merge")}}, extra={"settings": settings})
    info = {"cell": 1.0, "min_points": 3, "max_spread": 0.5, "filters": {"max_hs": 1.5, "max_day_offset": 0.15}}
    (out / "dem" / f"{date}{'_' + variant if variant else ''}_info.json").write_text(json.dumps(info))
    if variant:          # the table calibration's DEM record left in a variant's folder (other settings): not B's
        (out / "dem" / f"{date}_info.json").write_text(json.dumps(dict(info, cell=2.0)))
    sp.write_stamp(out, f"compare_{name}", {"after": {"dem": sp.token(out, "dem"), "waterlines": sp.token(out, "merge")}},
                   extra={"label": {"label": "INDEPENDENT", "reasons": []}})
    (out / "compare" / f"{date}_{name}_comparison.json").write_text(json.dumps(
        {"n": 120, "median": -0.034, "nmad": 0.057, "rmse": 0.08, "p5": -0.12, "p95": 0.07, "dem_cells": 900,
         "dem_cells_compared": 120, "photo_dates": [date, date], "days_outside_photo_window": 0}))
    prov = {"date": date, "era": "adcp", "status": "complete", "built_utc": "2026-10-14T21:00:00Z",
            "window": {"first_day": date, "last_day": date}, "setup": {"coef": sp.SETUP_COEF},
            "frames_per_day": {"c2": {date: 12}},
            "cameras": {"c2": {"calibration": {"variant": variant or "table"}, "photos": {"count": 12}}},
            "comparisons": [{"name": name, "survey_type": "dsm", "survey_date": date, "label": "INDEPENDENT",
                             "table_label": "INDEPENDENT", "status": "built", "current": True}]}
    (out / "provenance.json").write_text(json.dumps(prov))
    return out, prov


def test_summary_built(d):
    print("--summary and the end of --all over BUILT products with a comparison: a table root and a variant root")
    for variant in (None, "test-var"):
        root = Path(d) / f"summary_{variant or 'table'}"
        out, prov = summary_product(root, "2025-03-06", variant)
        pl = sp.product_plan(prov, out)
        check(sp.plan_variant(pl) == variant and sp.dem_stem(pl) == out / "dem" / ("2025-03-06" + (f"_{variant}" if
                                                                                                     variant else "")),
              f"a built product read back ({variant or 'table'}): its plan names the variant from provenance.json, so "
              f"its DEM record is dem/{sp.dem_stem(pl).name}_info.json")
        cur, why, h, lab = sp.summary_comparison(prov, out, prov["comparisons"][0])
        check(cur and not why and h.get("n") == 120 and lab.get("label") == "INDEPENDENT",
              f"summary_comparison on a complete product with a DEM and a comparison ({variant or 'table'} root): "
              f"current, no exception (was KeyError 'cams')")
        a = cfg_args(output_root=str(root))
        rc, txt = quiet(sp.summary, a)
        with open(root / "summary.csv", newline="") as f:
            rows = [r for r in csv.DictReader(f) if r.get("survey") == "lidar"]
        check(rc == 0 and len(rows) == 1 and rows[0]["current"] == "True" and "-0.034" in txt,
              f"--summary over that root: exit 0, the comparison current with its numbers ({variant or 'table'})")
        if variant:
            (out / "dem" / f"2025-03-06_{variant}_info.json").write_text(json.dumps(
                {"cell": 2.0, "min_points": 3, "max_spread": 0.5, "filters": {"max_hs": 1.5, "max_day_offset": 0.15}}))
            cur, why, _, _ = sp.summary_comparison(prov, out, prov["comparisons"][0])
            check(not cur and "2025-03-06_info.json: cell_m 2.0" in why,
                  "... and the variant's OWN DEM record is the one checked: another grid there makes it not current")
        # the end of --all: the per-date results, then the summary; its exit code is the builds' worst
        a = cfg_args(output_root=str(root), dry_run=False)
        rc, txt = quiet(sp.finish_all, {"2025-03-06": "built", "2024-10-23": "disabled"}, a)
        check(rc == 0 and "Traceback" not in txt and "2025-03-06" in txt,
              f"the end of --all ({variant or 'table'} root): the summary is made, exit 0 (the builds')")
    real = sp.summary

    def broken(args):
        raise KeyError("cams")
    sp.summary = broken
    try:
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            rc, txt = quiet(sp.finish_all, {"2025-03-06": "partial"}, cfg_args(output_root=str(root), dry_run=False))
    finally:
        sp.summary = real
    check(rc == sp.EXIT_CODES["partial"] and "KeyError" in err.getvalue() and "could not be made" in
          (txt + err.getvalue()),
          "a summary that fails at the end of --all is said, with its traceback; the exit code stays the builds' "
          "(3, partial), not 1")


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
    # 'another program': a command line the test controls (not this Python, whose path the test cannot
    # choose: an interpreter in a folder named survey_products would read as a build)
    other = subprocess.Popen(["sleep", "300"] if shutil.which("sleep") else
                             [sys.executable, "-c", "import time; time.sleep(300)"])
    # 'survey_products' inside other arguments: this test's own folder (survey_products_test_...), a file
    # named like the script but not it, the script's name with a suffix
    lookalike = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(300)", str(Path(d) / "x.txt"),
                                  str(Path(d) / "survey_products_notes.txt"), "survey_products.py.bak",
                                  "--survey_products"])
    abs_path = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(300)",
                                 str(HERE / "survey_products.py"), "--date", "2026-09-29"])
    try:
        shown = [command_line_shown(proc) for proc in (first, other, lookalike, abs_path)]
        check(all(shown), "the four stand-in processes have their command lines in /proc"
                          + ("" if all(shown) else f" (shown: {shown}, after 30 s)"))
        if Path(f"/proc/{other.pid}/cmdline").exists():
            check(sp.build_running(first.pid) and sp.build_running(abs_path.pid)
                  and not sp.build_running(other.pid) and not sp.build_running(lookalike.pid),
                  "build_running: an argument whose last part is survey_products.py is a build (as given, or a full "
                  "path); 'sleep 300', or 'survey_products' inside other arguments (the test's folder "
                  f"{Path(d).parent.name}, survey_products_notes.txt, survey_products.py.bak), is not")
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
        for proc in (first, other, lookalike, abs_path):
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
    wlf = dict(wl, frame_set={"dem_frames": 40, "left_out_frames": 9,
                              "left_out_why": {"offshore Hs > 1.5 m": 4, "camera-days off the rest by > 0.15 m": 5}})
    txt = " ".join("\n".join(sp.waterline_readme_lines(wlf)).split())
    check("it left out 9 (4 frame(s) with offshore Hs > 1.5 m, 5 frame(s) of camera-days off the rest by > 0.15 m)"
          in txt and "with camera-days" not in txt,
          "README FRAMES line: the left-out frames by reason in the comparison's own wording ('N frame(s) of "
          "camera-days ...', never 'N with camera-days')")
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
                                      dem_max_day_offset=0.15, geotiff_epsg=sp.GEOTIFF_EPSG), day_test=False)[1]
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


def test_oct_path(d):
    print("the October 2024 (chatham era) path, as built end to end on the synthetic October fixture")
    cfg = sp.read_table(sp.DEFAULT_CONFIG, sp.CONFIG_FIELDS)
    srv = sp.read_table(sp.DEFAULT_SURVEYS, sp.SURVEY_FIELDS)
    oct_ = [s for s in srv if s["name"] == "oct_gcps"][0]
    # The GCP advice names a folder the run searches: with --survey-dirs that leave out the station's
    # Chelsea folder, a file put there would never be found.
    sdir = Path(d) / "oct_surveys"
    sdir.mkdir()
    a = cfg_args(survey_dirs=[str(sdir)], force_disabled=True)
    plan, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, a)
    paths, note = sp.resolve_survey(oct_, plan, a)
    msg = sp.missing_survey_message(oct_, plan, a, note)
    want = f"{sdir}/2024-10-23_Marconi_Extrinsic_Targets_<cam>_xyz.csv"
    check(want in msg and str(sp.CHELSEA) not in msg and ("--force-disabled" in msg) == bool(plan["forced_disabled"]),
          "Oct GCPs missing, --survey-dirs given: the advice names that folder, not the station's Chelsea folder "
          "(--force-disabled in the rerun line only while the date is switched off)")
    check("--steps maps,compare" in msg,
          "... and reruns the maps too (a compare-only rerun left the plan view drawn without the GCPs)")
    (sdir / "2024-10-23_Marconi_Extrinsic_Targets_c1_xyz.csv").write_text("1,420100.0,4638400.0,1.20\r\n")
    paths, note = sp.resolve_survey(oct_, plan, a)
    check(len(paths) == 1, "... and a file put where the advice says is found")
    # the station's folders, Chelsea first (named only: missing_survey_message looks nothing up)
    a2 = cfg_args(force_disabled=True, survey_dirs=[str(x) for x in sp.DEFAULT_SURVEY_DIRS])
    msg = sp.missing_survey_message(oct_, plan, a2, "no survey file given")
    check(f"{sp.CHELSEA}/2024-10-23_Marconi_Extrinsic_Targets_<cam>_xyz.csv" in msg,
          "... and with the station's own folders, the Chelsea folder (searched first)")
    check("not on the NUC" not in oct_["why"],
          "the Oct GCP label reason says nothing about the file being absent (it is printed once the file is there)")

    # The GCPs arrive as one CRLF file per camera (as the 2025-11-13 target files), a target seen by
    # both cameras in both files.
    (sdir / "2024-10-23_Marconi_Extrinsic_Targets_c1_xyz.csv").write_bytes(
        b"101,420131.051,4638335.088,1.317\r\n301,420135.051,4638385.088,0.816\r\n")
    (sdir / "2024-10-23_Marconi_Extrinsic_Targets_c2_xyz.csv").write_bytes(
        b"205,420126.801,4638440.088,1.034\r\n301,420135.051,4638385.088,0.816\r\n")
    paths, _ = sp.resolve_survey(oct_, plan, a)
    check(len(paths) == 2, "both camera files found by the survey date")
    # (1) the merged file is rewritten only when its content changes: it was rewritten on every run
    # (csv's '\r\n' never equalled a text read), so the comparison was redone on every run
    dst = sdir / "merged.csv"
    sp.merged_points(paths, dst)
    os.utime(dst, (1.0e9, 1.0e9))
    sp.merged_points(paths, dst)
    check(dst.stat().st_mtime == 1.0e9 and len(dst.read_text().splitlines()) == 3,
          "merged GCP file: 3 targets (301 once), not rewritten when nothing changed (keeps its age)")
    # (2) the label reasons: the calibration-day rule ran once per FILE and listed every reason twice
    v = sp.Verdict(oct_)
    for p_ in paths:
        sp.honest_label(oct_, plan, a, p_, fits=[], verdict=v)
    check(v.label == "PARTLY-CIRCULAR" and len(v.reasons) == len(set(v.reasons)) == 2,
          f"two GCP files: PARTLY-CIRCULAR, each camera's calibration-day reason once ({len(v.reasons)} reasons)")
    # (3) the plan view draws every file's points, not only the first file's
    pts = sp.plan_survey_points(oct_, plan, a, sdir)
    check(len(pts["E"]) == 3 and np.isclose(sorted(pts["Z"])[-1], 1.317) and np.isclose(sorted(pts["Z"])[0], 0.816),
          f"plan view: the targets of both camera files, 301 once ({len(pts['E'])} points)")

    # (4) the summary of a date whose survey file has not arrived: never compared, not 'not current'
    root = Path(d) / "oct_root"
    pdir = root / "2024-10-23"
    pdir.mkdir(parents=True)
    (root / "summary.png").write_text("an old figure")
    reason = sp.missing_survey_message(oct_, plan, a, "no survey file given in surveys.csv (column path)")
    (pdir / "provenance.json").write_text(sp.json.dumps({
        "date": "2024-10-23", "era": "chatham", "status": "complete",
        "window": {"first_day": "2024-10-20", "last_day": "2024-10-26"}, "cameras": {}, "frames_per_day": {},
        "comparisons": [{"name": "oct_gcps", "survey_path": "", "survey_type": "points", "survey_date": "2024-10-23",
                         "table_label": "PARTLY-CIRCULAR", "label": "PARTLY-CIRCULAR", "why": oct_["why"],
                         "status": "skipped", "reason": reason, "current": False}]}))
    _, out = quiet(sp.summary, SimpleNamespace(output_root=str(root), config=str(sp.DEFAULT_CONFIG),
                                               surveys=str(sp.DEFAULT_SURVEYS), setup_coef=sp.SETUP_COEF,
                                               survey_dirs=[str(sdir)], calibration=str(HERE / "calibration")))
    rows = list(csv.DictReader(open(root / "summary.csv")))
    row = [r for r in rows if r["survey"] == "oct_gcps"][0]
    check(row["label"] == "no survey file" and "no survey file yet" in row["comparison"]
          and "not current" not in row["comparison"],
          f"summary: the Oct GCPs NOT COMPARED for want of the file, not 'not current' ({row['label']!r})")
    check(not (root / "summary.png").exists() and "no comparison has numbers yet" in out,
          "summary: no figure when no comparison has numbers (and no old one left to look current)")

    # (5) the README's water-level caveat for the chatham era: EVERY level is transferred (it said
    # 'Some'), the period lies 44-51 days before the ADCP record, further than the extrapolation tests
    # reach (0-23 and 0-29 days), and the transfer still carries the ADCP's datum
    forcing = {"water_level": {
        "sources_used": {"chatham_transfer": 1681}, "method": "harmonic", "cv_rms_m": 0.0756,
        "extrapolation_rms_m": 0.0851, "fit_window": ["2024-12-09 18:00", "2025-03-10 15:00"],
        "period_vs_adcp": "44-51 days BEFORE the ADCP record (transfer extrapolated)",
        "extrapolation_reach_days": [[0, 23], [0, 29]],
        "beyond_tested": "this period lies 44-51 days from the ADCP record, further than the extrapolation tests reach",
        "adcp_datum": {"source": "adcp_water_level_navd88.csv (adcp_to_navd88.py)"}}, "waves": {}}
    cav = sp.collect_caveats(plan, SimpleNamespace(setup_coef=sp.SETUP_COEF), {}, forcing, {})
    wlc = [c for c in cav if "Chatham gauge 8447435 transferred" in c]
    check(len(wlc) == 1 and wlc[0].startswith("Every water level") and "44-51 days BEFORE" in wlc[0]
          and "0-23 and 0-29 days" in wlc[0] and "FURTHER THAN TESTED" in wlc[0] and "Some water" not in " ".join(cav),
          "README: every level of the Oct date is the Chatham transfer, 44-51 days out, FURTHER THAN TESTED")
    check(any(c.startswith("ADCP datum") and "carries the same datum" in c for c in cav),
          "README: no ADCP level used, but the transfer carries the ADCP's datum")
    forcing["water_level"].update(sources_used={"adcp": 300, "chatham_transfer": 12}, period_vs_adcp="96% covered",
                                  beyond_tested=None)
    cav = sp.collect_caveats(plan, SimpleNamespace(setup_coef=sp.SETUP_COEF), {}, forcing, {})
    check(any(c.startswith("12 of 312 water levels") and "FURTHER THAN TESTED" not in c for c in cav),
          "README: a few transferred hours inside the ADCP period are counted, not 'every'")


# forcing.json's water_level of the real Oct 2024 build, as historical_forcing.py of f4ba895 (main, on
# the NUC) writes it: the tests' record and the fit window, none of the fields a newer one adds
OCT_WL_OLD = {
    "sources_used": {"chatham_transfer": 1681}, "method": "harmonic", "cv_rms_m": 0.0756,
    "extrapolation_rms_m": 0.0851, "extrapolation_rms_daytime_m": 0.0806, "extrapolation_bias_m": 0.0139,
    "transfer_sigma_m": 0.0851, "fit_window": ["2024-12-09 18:00", "2025-03-10 15:00"],
    "period_vs_adcp": "44 days BEFORE the ADCP record (transfer extrapolated)",
    "adcp_datum": {"source": "adcp_water_level_navd88.csv (adcp_to_navd88.py)"},
    "methods": {"harmonic": {"extrapolation_tests": [
        {"n": 534, "rms_m": 0.095, "bias_m": 0.0375, "rms_daytime_m": 0.0849,
         "test": "fit 2025-01-01..2025-03-11 -> predict 2024-12-09..2025-01-01"},
        {"n": 688, "rms_m": 0.0765, "bias_m": -0.0045, "rms_daytime_m": 0.0772,
         "test": "fit 2024-12-09..2025-02-10 -> predict 2025-02-10..2025-03-11"}]}}}
OCT_PERIOD = {"start": "2024-10-20", "end": "2024-10-26", "start_utc": "2024-10-20 00:00",
              "end_utc": "2024-10-27 00:00"}


def test_oct_transfer_caveat(d):
    print("the Oct 2024 water-level caveat, worked out from forcing.json itself (old or new)")
    sdir = Path(d) / "oct_caveat_surveys"
    sdir.mkdir()
    cfg = sp.read_table(sp.DEFAULT_CONFIG, sp.CONFIG_FIELDS)
    srv = sp.read_table(sp.DEFAULT_SURVEYS, sp.SURVEY_FIELDS)
    plan, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, cfg_args(survey_dirs=[str(sdir)], force_disabled=True))
    a = SimpleNamespace(setup_coef=sp.SETUP_COEF)
    # (1) a forcing.json built by the code on the NUC now (f4ba895): the README rewritten by this code
    # (--steps maps,compare after a pull, no forcing rebuild) still gives the reach, the backward test
    # and FURTHER THAN TESTED: it said only '44 days BEFORE the ADCP record' (review, finding 1)
    forcing = {"period": dict(OCT_PERIOD), "water_level": dict(OCT_WL_OLD), "waves": {}}
    cav = [c for c in sp.collect_caveats(plan, a, {}, forcing, {}) if "Chatham gauge 8447435 transferred" in c]
    c = cav[0] if cav else ""
    check(len(cav) == 1 and c.startswith("Every water level") and "backward extrapolation test" in c
          and "0-23 days before" in c and "0.095 m" in c and "+0.037" in c and "FURTHER THAN TESTED" in c
          and "44-51 days before" in c, "old forcing.json: backward test 0-23 days, RMS 0.095, bias +0.037, "
                                        "44-51 days out, FURTHER THAN TESTED")
    check("0-29" not in c and "0-23 and 0-29" not in c,
          "... judged by the backward test's reach only, not the farthest of both tests")
    check("0.0851 m, is the pooled RMS of both tests (an older historical_forcing.py)" in c
          and "this side's test gives 0.095 m" in c,
          "... and the 0.0851 m sigma_m an older forcing step wrote is named as the pooled RMS, not this side's")
    check("not averaged away" in c and "week left out" not in c,
          "... the shared (window-mean) error said; no week-mean numbers an old forcing.json lacks")
    # (2) a forcing.json of this code: sigma_m is the backward test's, the week-mean error is given
    wl = dict(OCT_WL_OLD, transfer_sigma_m=0.095, transfer_sigma_measured=False,
              period_vs_adcp="44-51 days BEFORE the ADCP record (transfer extrapolated)",
              week_mean_error={"weeks": 13, "rms_m": 0.0506, "min_m": -0.0892, "max_m": 0.1016},
              beyond_tested="this period lies 44-51 days before the transfer's fit window, ...")
    forcing = {"period": dict(OCT_PERIOD), "water_level": wl, "waves": {}}
    c = [x for x in sp.collect_caveats(plan, a, {}, forcing, {}) if "Chatham gauge 8447435 transferred" in x][0]
    check("carries 0.095 m as its sigma_m" in c and "older" not in c and "-0.089 .. +0.102 m" in c
          and "RMS 0.051 m" in c and "13 weeks" in c,
          "new forcing.json: each level carries 0.095 m; a whole week's mean error -0.089..+0.102 m (RMS 0.051)")
    check(c.endswith("Every line moves with it, and the DEM by about 0.6-0.8 x of it once re-projected (a higher water "
                     "level also places the line further landward, on higher beach: the share of a level shift a line "
                     "keeps, as of the setup in the paired C = 0 builds)."),
          "... a level error moves every line, and the DEM by about 0.6-0.8 x of it (not 1:1)")
    # the ADCP datum: a size from the evidence there is; converted waves: the setup error they cause (review r1)
    wl2 = dict(wl, adcp_datum={"offset_from_chatham_mean_m": 0.101, "offset_from_msl_constant_m": 0.1127,
                               "offset_used_m": 0.1012, "source": "adcp_to_navd88.py"})
    wv = {"sources_used": {"hs=wis+ndbc44013;tp=wis": 169},
          "source_labels": {"hs=wis+ndbc44013;tp=wis": "Hs: WIS 63064 + NDBC 44013 converted; Tp: WIS 63064"},
          "setup_in_window": {"daytime_median_m": 0.361},
          "combination": {"common_hours": 532, "common_period": ["2024-12-09 18:00", "2024-12-31 23:00"],
                          "setup_cv_on_common_m": {"wis": 0.0964, "ndbc44013": 0.0791},
                          "setup_blend_cv": {"n": 532, "rms_m": 0.079, "bias_m": 0.0041, "p95_abs_m": 0.1844}}}
    cav = sp.collect_caveats(plan, a, {}, {"period": dict(OCT_PERIOD), "water_level": wl2, "waves": wv}, {})
    dat = next((x for x in cav if x.startswith("ADCP datum")), "")
    wcv = next((x for x in cav if x.startswith("Waves converted")), "")
    check("the DEM about 0.6-0.8 x of it" in dat and "differ by 0.012 m" in dat and "a few cm" in dat
          and "~0.02 m below Chatham" in dat,
          "ADCP datum: a size from the evidence (0.012 m between the two datum routes; GNSS-R ~0.02 m below Chatham): "
          "a few cm, not measured; the DEM about 0.6-0.8 x of it")
    check(wcv.startswith("Waves converted: every wave hour(s) of this date are Hs: WIS 63064 + NDBC 44013 converted")
          and "532 h all sources share (2024-12-09 18:00 .. 2024-12-31 23:00 UTC) only: RMS 0.079 m, bias +0.004 m"
          in wcv and "ndbc44013 alone 0.079 m, wis alone 0.096 m" in wcv and "setup of 0.36 m" in wcv
          and "~0.05-0.06 m RMS per frame" in wcv and "shared by every frame of the window" in wcv
          and sp.wave_conversion_caveat({"sources_used": {"adcp": 171}}) is None,
          "converted waves: the setup error they cause (CV RMS and bias, December only), what a line keeps of it, and "
          "a seasonal bias shared by every frame; nothing for an ADCP-only date")
    # (3) an ADCP-era date with a few transferred hours inside the record: the pooled RMS, no side
    wl = dict(OCT_WL_OLD, sources_used={"adcp": 300, "chatham_transfer": 12}, period_vs_adcp="96% covered")
    forcing = {"period": {"start_utc": "2025-01-16 00:00", "end_utc": "2025-01-24 00:00"}, "water_level": wl,
               "waves": {}}
    c = [x for x in sp.collect_caveats(plan, a, {}, forcing, {}) if "Chatham gauge 8447435 transferred" in x][0]
    check(c.startswith("12 of 312 water levels") and "backward" not in c and "FURTHER THAN TESTED" not in c
          and "Extrapolation RMS 0.0851 m (pooled over both directions)" in c and "not averaged away" not in c,
          "a hole inside the record: the pooled RMS, no side, no FURTHER THAN TESTED, no window-mean sentence")
    # (4) the README's WHAT THIS IS: no 'measured water level' for a date none of whose levels is
    w = sp.level_words({"water_level": {"sources_used": {"chatham_transfer": 1681}}})
    check("measured water level" not in w and "Chatham gauge's water level transferred" in w,
          "README: the Oct levels are called the Chatham gauge transferred, not measured")
    check(sp.level_words({"water_level": {"sources_used": {"adcp": 1500}}}) == "the measured water level "
          and "12 of 312" in sp.level_words({"water_level": {"sources_used": {"adcp": 300, "chatham_transfer": 12}}}),
          "... the ADCP's levels are measured; a mixed date says how many are transferred")


def test_gcp_files_any_format(d):
    print("GCP files in any format survey_compare.py reads, merged on their coordinates")
    import survey_compare as sc
    d = Path(d) / "gcp_formats"
    d.mkdir()
    hdr = "Name,Code,Easting,Northing,Elevation,Description,Longitude,Latitude,Ellipsoidal height,CS name\r\n"
    crs = "NAD83(2011) / UTM zone 19N + NAVD88(GEOID18) height"
    e1 = d / "emlid_c1.csv"
    e1.write_text(hdr + f"T1,gcp,420134.298,4638335.154,2.274,target,-69.9,41.9,-26.0,{crs}\r\n"
                        f"T2,gcp,420140.100,4638380.200,1.105,target,-69.9,41.9,-27.1,{crs}\r\n")
    e2 = d / "emlid_c2.csv"
    e2.write_text(hdr + f"T2,gcp,420140.100,4638380.200,1.105,target,-69.9,41.9,-27.1,{crs}\r\n"
                        f"T9,gcp,420126.801,4638440.088,0.834,target,-69.9,41.9,-27.4,{crs}\r\n")
    check(sc.read_points(e1)["format"] == "Emlid Flow CSV", "an Emlid file alone reads as Emlid")
    m = sp.merged_points([e1, e2], d / "m_emlid.csv")
    q = sc.read_points(m)
    check(len(q["E"]) == 3 and sorted(q["name"]) == ["T1", "T2", "T9"] and set(q["crs"]) == {crs}
          and np.isclose(q["Z"][list(q["name"]).index("T9")], 0.834),
          "two Emlid files merged: 3 targets (T2 once), names and CRS kept, readable (it failed: 'no "
          "Easting/Northing/Elevation columns')")
    from compare_rtk import load_transects
    check(load_transects(str(m), 5.0) == [], "... and compare_rtk's transect reader reads the merged file")
    # headers in another order: every target kept (the second file's were dropped)
    a1 = d / "hdr_a.csv"
    a1.write_text("Name,Easting,Northing,Elevation\nA1,420100.0,4638400.0,1.20\nA2,420110.0,4638410.0,0.80\n")
    a2 = d / "hdr_b.csv"
    a2.write_text("Elevation,Northing,Easting,Name\n0.50,4638500.0,420120.0,B1\n0.40,4638510.0,420125.0,B2\n")
    q = sc.read_points(sp.merged_points([a1, a2], d / "m_hdr.csv"))
    k = list(q["name"]).index("B1")
    check(len(q["E"]) == 4 and np.isclose(q["E"][k], 420120.0) and np.isclose(q["N"][k], 4638500.0)
          and np.isclose(q["Z"][k], 0.50), "two header files, columns in another order: all 4 targets, right axes")
    # CIRN rows with target names that are not numbers ('T1'): read alone (it failed) and merged
    c1 = d / "cirn_alnum.csv"
    c1.write_bytes(b"T1,420134.298,4638335.154,2.274\r\nGCP-2,420140.100,4638380.200,1.105\r\n")
    q = sc.read_points(c1)
    check(q["format"].startswith("CIRN") and list(q["name"]) == ["T1", "GCP-2"] and np.isclose(q["Z"][0], 2.274),
          "CIRN rows named T1, GCP-2: read as CIRN, names kept (it exited: 'no Easting/Northing/Elevation')")
    c2 = d / "cirn_num.csv"
    c2.write_bytes(b"301,420135.051,4638385.088,0.816\r\nGCP-2,420140.100,4638380.200,1.105\r\n")
    m = sp.merged_points([c1, c2], d / "m_cirn.csv")
    q = sc.read_points(m)
    check(q["format"].startswith("CIRN") and sorted(q["name"]) == ["301", "GCP-2", "T1"]
          and not Path(m).read_text().startswith("Name"),
          "CIRN files with mixed names merged: still CIRN (no header), 3 targets")
    q = sc.read_points(sp.merged_points([c1, e2], d / "m_mixed.csv"))
    check(len(q["E"]) == 3 and sorted(q["name"]) == ["GCP-2", "T1", "T9"] and crs in q["crs"],
          "a CIRN file and an Emlid file merged: one file, the CRS of the Emlid targets kept, GCP-2/T2 once")
    # a file survey_compare.py cannot read: SurveyFileError naming it, never SystemExit (which escaped
    # step_maps' 'except Exception' and failed the whole product)
    bad = d / "empty_c2.csv"
    bad.write_text("\n")
    try:
        sp.merged_points([c1, bad], d / "m_bad.csv")
        check(False, "an unreadable file raises")
    except sp.SurveyFileError as exc:
        check(str(bad) in str(exc), f"an unreadable file: SurveyFileError naming it ({exc})")
    words = d / "words.csv"
    words.write_text("Name,Code\nT1,gcp\n")
    try:
        sp.read_survey_points(words)
        check(False, "no coordinate columns raises")
    except sp.SurveyFileError as exc:
        check("Easting/Northing/Elevation" in str(exc), "no coordinate columns: SurveyFileError, the reason kept")
    # a kill during the read still ends the run (the SIGTERM handler raises SystemExit)
    real = sc.read_points
    try:
        def killed(path):
            raise SystemExit("terminated by signal 15 (SIGTERM)")
        sc.read_points = killed
        try:
            sp.read_survey_points(c1)
            check(False, "a termination passes through")
        except sp.SurveyFileError:
            check(False, "a termination was taken for a bad file")
        except SystemExit as exc:
            check(sp.is_termination(exc), "a SIGTERM during the read is not taken for a bad file: the run ends")
    finally:
        sc.read_points = real

    # the plan view with a bad survey file: drawn without its points, with a warning (it failed the build)
    src = d / "ground.csv"
    src.write_text("source_file,camera,capture_time_utc,tide_elevation_navd88,beach_elevation_navd88,"
                   "easting_utm19,northing_utm19\n"
                   + "".join(f"f{k},c1,2024-10-2{k % 5}T15:00:00Z,0.1,{0.2 + 0.1 * k:.2f},{420100 + k},{4638380 + k}\n"
                             for k in range(20)))
    row = {"name": "oct_gcps", "survey_type": "points", "survey_date": "2024-10-23", "label": "PARTLY-CIRCULAR",
           "path": ""}
    sdir = d / "surveys"
    sdir.mkdir()
    shutil.copy(c1, sdir / "2024-10-23_Marconi_Extrinsic_Targets_c1_xyz.csv")
    (sdir / "2024-10-23_Marconi_Extrinsic_Targets_c2_xyz.csv").write_text("\n")
    a = cfg_args(survey_dirs=[str(sdir)], force_disabled=True)
    plan = {"date": "2024-10-23", "era": "chatham", "first": "2024-10-20", "last": "2024-10-26", "out": d,
            "cams": {}, "surveys": [row]}
    check(len(sp.resolve_survey(row, plan, a)[0]) == 2, "both camera files found")
    png = d / "plan.png"
    _, out = quiet(sp.plan_view_map, plan, src, [row], a, png)
    check(png.exists() and "points of oct_gcps are NOT drawn" in out and "2024-10-23_Marconi_Extrinsic_Targets_c2"
          in out, "plan view: an unreadable GCP file costs its points (warned), not the map or the build")
    (sdir / "2024-10-23_Marconi_Extrinsic_Targets_c2_xyz.csv").write_bytes(
        b"T9,420126.801,4638440.088,0.834\r\n")
    png.unlink()
    _, out = quiet(sp.plan_view_map, plan, src, [row], a, png)
    check(png.exists() and "NOT drawn" not in out, "... and with both files readable, drawn with the points")


def test_merge_level_source(d):
    print("the merged waterlines name the forcing's level source, not 'tide_model'")
    d = Path(d) / "lvlsrc"
    (d / "waterlines" / "c1").mkdir(parents=True)
    (d / "forcing").mkdir()
    (d / "forcing" / "water_level.csv").write_text(
        "time,water_level_navd88,source,sigma_m\n2024-10-23T15:00:00Z,0.40,chatham_transfer,0.095\n"
        "2024-10-23T15:06:00Z,0.41,chatham_transfer,0.095\n")
    (d / "forcing" / "waves.csv").write_text(
        "time_utc,epoch,wvht_m,dpd_s,mwd_deg,source,hs_sigma_m,tp_sigma_s\n"
        "2024-10-23T15:00:00Z,1729695600,1.0,9.0,90,wis+ndbc44013,0.3,1.0\n")
    f = d / "waterlines" / "c1" / "contour_points_ground.csv"
    f.write_text("source_file,camera,capture_time_utc,capture_epoch,tide_elevation_navd88,tide_model_spread,"
                 "water_level_source,offshore_hs_m,offshore_tp_s,setup_correction_m,beach_elevation_navd88,"
                 "easting_utm19,northing_utm19\n"
                 "a,c1,2024-10-23T15:00:00Z,1729695600,0.40,0.0,tide_model,1.0,9.0,0.3,0.70,420100,4638400\n")
    plan = {"date": "2024-10-23", "era": "chatham", "out": d, "first": "2024-10-20", "last": "2024-10-26",
            "cams": {"c1": {}}, "surveys": []}
    quiet(sp.merge_cameras, plan, [("c1", f)], SimpleNamespace(setup_coef=0.037))
    r = list(csv.DictReader(open(d / "waterlines" / "contour_points_ground.csv")))[0]
    check(r["water_level_source"] == "chatham_transfer" and r["forcing_level_source"] == "chatham_transfer"
          and r["tide_model_spread"] == "" and r["forcing_level_sigma_m"] == "0.095",
          "chatham era: water_level_source = chatham_transfer (it said tide_model), no tide-model spread")
    plan["era"] = "live"
    f.write_text(f.read_text().replace(",tide_model,", ",gnssr_measured,"))
    quiet(sp.merge_cameras, plan, [("c1", f)], SimpleNamespace(setup_coef=0.037))
    r = list(csv.DictReader(open(d / "waterlines" / "contour_points_ground.csv")))[0]
    check(r["water_level_source"] == "gnssr_measured" and "forcing_level_source" not in r,
          "live era: the station's own source kept, no forcing columns")


def test_dry_run_needs_adcp():
    print("the dry run marks the ADCP as required for the October (chatham) date")
    found = {"adcp": None, "adcp_navd88": None, "chatham": Path("chatham.csv"), "wis": Path("wis.csv"), "ndbc": None}
    lines, warns = sp.forcing_input_notes("chatham", found, sp.SETUP_COEF)
    t = dict(lines)
    check("REQUIRED" in t["adcp"] and "derived" in t["adcp_navd88"],
          "chatham era, no ADCP: listed as REQUIRED (it was a bare 'NOT FOUND'); its NAVD88 file optional")
    check(any("no ADCP record" in w and "FAIL" in w for w in warns), "... and a warning that the forcing step fails")
    found["adcp"] = Path("sig1000_waves_ALL.csv")
    lines, warns = sp.forcing_input_notes("chatham", found, sp.SETUP_COEF)
    check(not warns, "with the ADCP and WIS: no warning")
    found["wis"] = None
    _, warns = sp.forcing_input_notes("chatham", found, sp.SETUP_COEF)
    check(any("no WIS or NDBC" in w for w in warns), "no WIS and no NDBC: the wave warning stays")
    _, warns = sp.forcing_input_notes("adcp", {"adcp": None, "adcp_navd88": None}, sp.SETUP_COEF)
    check(not warns, "ADCP era: no chatham-era warnings")


def copy_tracked_calibration(dst, repo=HERE):
    """The repository's calibration folder as COMMITTED (git ls-files), into dst: never the station's own untracked
    files there (review r1: copytree opened all 49 of them on the simulated NUC), and never a file named like a held
    survey (survey_products.held_file_name, by its name: such a file is not even opened), whatever git holds. Without
    git (a copy made by git archive) the folder as it is, by the same rule. -> how the files were listed."""
    src = Path(repo) / "calibration"
    r = subprocess.run(["git", "-C", str(repo), "ls-files", "-z", "--", "calibration"], stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
    if r.returncode == 0 and r.stdout:
        rels = [Path(x).relative_to("calibration") for x in r.stdout.decode("utf-8", "replace").split("\0") if x]
        how = "the committed files (git ls-files)"
    else:
        rels = [q.relative_to(src) for q in sorted(src.rglob("*")) if not q.is_dir()]
        how = "the folder as it is (no git here)"
    for rel in rels:
        if any(sp.held_name(part) for part in rel.parts):
            continue
        q = src / rel
        if not q.is_file():                         # tracked but deleted in this checkout
            continue
        (dst / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(q), str(dst / rel))
    dst.mkdir(parents=True, exist_ok=True)
    return how


def test_blind_and_calibration(d):
    print("a BLIND date (DEM only: its survey is held for an independent check), its two calibrations, one per root")
    import builtins
    import pathlib
    import estimate_eo_rotation as er
    from georectify import load_extrinsics
    d = Path(d) / "blind"
    d.mkdir()
    cfg, srv, tabs = fixture_tables(d)
    oct_row = [s for s in srv if s["name"] == "oct_gcps"][0]
    held = sp.BLIND_DATES["2024-10-23"]
    # canaries: files named like the held October GCPs, in every folder a build could search (the survey
    # folders, the calibration folder, the output roots, the current folder). Never looked for, never opened.
    sdirs = [d / "survey_a", d / "survey_b"]
    cal = d / "calibration"
    how_cal = copy_tracked_calibration(cal)
    check((cal / "calibration_history.csv").is_file() and (cal / "CACO03_c1_20241023_EO.yaml").is_file()
          and not any(sp.held_name(q.name) for q in cal.rglob("*")),
          f"the test's calibration folder: {how_cal} of the repository's, none named like the held survey (the "
          f"station's own untracked files there are never opened)")
    # ... on a repository of its own: a tracked file is copied, an untracked one (the station's own) is not
    rp = d / "cal_repo"
    (rp / "calibration" / "sub").mkdir(parents=True)
    (rp / "calibration" / "CACO03_c9_20250123_EO.yaml").write_text("x: 1\n")
    (rp / "calibration" / "sub" / "kept.txt").write_text("k\n")
    git_in(rp, "init", "-q")
    git_in(rp, "add", "-A")
    git_in(rp, "commit", "-q", "-m", "tracked")
    (rp / "calibration" / "station_own_EO.yaml.bak").write_text("untracked\n")
    how2 = copy_tracked_calibration(d / "cal_copy", rp)
    check(how2.startswith("the committed files") and (d / "cal_copy" / "CACO03_c9_20250123_EO.yaml").is_file()
          and (d / "cal_copy" / "sub" / "kept.txt").is_file() and not (d / "cal_copy" / "station_own_EO.yaml.bak").exists(),
          "copy_tracked_calibration: the committed files only (an untracked station file is never copied or opened)")
    # the calibration history of the test's copy: the rows of the solution from the held GCPs get a height
    # no camera has (999 m): it must never reach a README (those rows are not read on a blind date)
    hist = cal / "calibration_history.csv"
    rows_h = list(csv.DictReader(open(hist, newline="")))
    for r in rows_h:
        if "20241023" in r["file"]:
            r["z"] = "999.000"
    with open(hist, "w", newline="") as f:
        w_ = csv.DictWriter(f, fieldnames=list(rows_h[0].keys()))
        w_.writeheader()
        w_.writerows(rows_h)
    roots = {k: d / k for k in ("out", "out_turned", "outb_turned", "fresh", "fresh_turned", "proc_proc",
                                "fresh_proc")}
    for f in sdirs + list(roots.values()):
        f.mkdir()
    canary_names = ["2024-10-23_Marconi_Extrinsic_Targets_c1_xyz.csv", "2024-10-23_Marconi_Extrinsic_Targets_c2_xyz.csv",
                    "2024-10-23_GCP_points.csv", "targets_20241023_xyz.csv"]
    canaries = set()
    for f in sdirs + [cal, roots["out"], Path.cwd()]:
        for n in canary_names:
            (f / n).write_text("1,420100.0,4638400.0,1.20\n")          # writing is allowed: only reading is refused
            canaries.add(os.path.realpath(str(f / n)))
    (sdirs[0] / "2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif").write_text("not a real lidar")
    forcing = {}
    for k in ("adcp", "adcp_navd88", "chatham", "wis", "ndbc"):
        forcing[k] = str(d / f"{k}.csv")
        Path(forcing[k]).write_text("")
    a = cfg_args(output_root=str(roots["out"]), survey_dirs=[str(x) for x in sdirs], calibration=str(cal),
                 force_disabled=True, photo_roots=[str(d / "no_photos")], live_contours="/x/contour_points_timex.csv",
                 setup_fitted_to=None, geotiff_epsg=sp.GEOTIFF_EPSG, gnssr_spline=str(d / "no_spline.txt"),
                 gauge_csv=str(d / "no_gauge.csv"), live_cron=str(sp.LIVE_CRON), pointing_log_dir=str(d / "no_log"),
                 pointing_frames=3, pointing_fail=sp.POINT_FAIL_DEG, dem_cell=sp.DEM_CELL,
                 dem_min_points=sp.DEM_MIN_POINTS, dem_max_spread=sp.DEM_MAX_SPREAD, dem_max_hs=sp.DEM_MAX_HS,
                 dem_max_day_offset=sp.DEM_MAX_DAY_OFFSET, keep_moved_days=False, no_download=True, synthetic=None,
                 replace=False, dry_run=False, force=False, eo=None, calibration_variant=None, **forcing, **tabs)

    # -- the tracer: every look-up (a glob of the held date, a stat or open of a file named like it) through
    # pathlib, os and open(), and every open() of a canary seen by an audit hook (also C-level opens: pandas,
    # numpy), and every command run (survey_compare.py, or a canary named on a command line)
    seen = []
    trace = {"on": False}

    def held_name(x):          # a survey-like file (never a calibration .yaml) named with the held date
        n = os.path.basename(str(x)).lower()
        return (("2024-10-23" in n or "20241023" in n) and any(k in n for k in ("xyz", "gcp", "extrinsic", "target"))
                and not n.endswith((".yaml", ".yml")))

    def held_pattern(pat):
        return any(k in str(pat) for k in ("2024-10-23", "20241023"))
    orig = {"glob": pathlib.Path.glob, "rglob": pathlib.Path.rglob, "stat": pathlib.Path.stat,
            "open": pathlib.Path.open, "io_open": io.open, "b_open": builtins.open, "os_stat": os.stat}

    def wrap(name, fn, kind):
        def inner(*args, **kw):
            if kind == "pattern" and len(args) > 1 and held_pattern(args[1]):
                seen.append(f"{name}({args[0]}, {args[1]!r})")
            elif kind == "path" and args and held_name(args[0]):
                seen.append(f"{name}({args[0]})")
            return fn(*args, **kw)
        return inner

    watch_a = {"on": False, "paths": {os.path.realpath(str(cal / f"CACO03_{c}_20241023_EO.yaml")) for c in ("c1", "c2")},
               "seen": []}

    def audit(event, args_):
        if event == "open" and watch_a["on"] and args_ and isinstance(args_[0], (str, bytes, os.PathLike)):
            if os.path.realpath(os.fsdecode(args_[0])) in watch_a["paths"]:
                watch_a["seen"].append(os.fsdecode(args_[0]))
        if not trace["on"]:
            return
        if event == "open" and args_ and isinstance(args_[0], (str, bytes, os.PathLike)):
            if os.path.realpath(os.fsdecode(args_[0])) in canaries:
                seen.append(f"audit open({os.fsdecode(args_[0])})")
        elif event == "subprocess.Popen":
            argv_ = [str(x) for x in (args_[1] or [])]
            words = " ".join(argv_)
            if (len(argv_) > 1 and argv_[1].endswith("survey_compare.py")) or any(n in words for n in canary_names):
                seen.append(f"audit subprocess({words[:120]})")       # a comparison run, or a canary named
    if hasattr(sys, "addaudithook"):
        sys.addaudithook(audit)

    @contextlib.contextmanager
    def traced():
        pathlib.Path.glob = wrap("Path.glob", orig["glob"], "pattern")
        pathlib.Path.rglob = wrap("Path.rglob", orig["rglob"], "pattern")
        pathlib.Path.stat = wrap("Path.stat", orig["stat"], "path")
        pathlib.Path.open = wrap("Path.open", orig["open"], "path")
        io.open = builtins.open = wrap("open", orig["io_open"], "path")
        os.stat = wrap("os.stat", orig["os_stat"], "path")
        trace["on"] = True
        try:
            yield
        finally:
            trace["on"] = False
            pathlib.Path.glob, pathlib.Path.rglob = orig["glob"], orig["rglob"]
            pathlib.Path.stat, pathlib.Path.open = orig["stat"], orig["open"]
            io.open, builtins.open, os.stat = orig["io_open"], orig["b_open"], orig["os_stat"]

    # the tracer catches what it should: an UNBLINDED look-up of the same date is seen (glob, stat, open)
    sp.BLIND_DATES.pop("2024-10-23")
    try:
        with traced():
            hits = sp.resolve_survey(oct_row, None, a)[0]
            sp.merged_points(hits, d / "merged_unblinded.csv")
    finally:
        sp.BLIND_DATES["2024-10-23"] = held
    check(len(hits) >= 2 and any(x.startswith("Path.glob") for x in seen)
          and (any(x.startswith("audit open") for x in seen) or not hasattr(sys, "addaudithook")),
          f"the tracer sees an (unblinded) look-up and read of the October GCPs ({len(seen)} calls)")
    seen.clear()
    (d / "merged_unblinded.csv").unlink()

    def fake_product(plan_, args_):
        """Stamps of every step, waterlines and maps of a two-camera build (no detection is run)."""
        outd = plan_["out"]
        wdir = outd / "waterlines"
        wdir.mkdir(parents=True)
        rows = ["source_file,camera,capture_time_utc,easting_utm19,northing_utm19,tide_elevation_navd88,"
                "beach_elevation_navd88,setup_correction_m"]
        for k in range(40):
            for cam in ("c1", "c2"):
                rows.append(f"f{k}_{cam},{cam},2024-10-2{k % 7} 15:00:00,{420150 + k},{4638350 + 3 * k},"
                            f"{0.1 * (k % 9) - 0.3:.2f},{0.1 * (k % 9) - 0.1:.2f},0.2")
        (wdir / "contour_points_ground.csv").write_text("\n".join(rows) + "\n")
        sp.write_stamp(outd, "forcing", {"x": 1})
        fit = sp.setup_fit_info(sp.SETUP_COEF)
        for cam, c in plan_["cams"].items():
            rec = sp.calibration_record(c)
            sp.write_stamp(outd, f"pointing_{cam}", {"eo_file": sp.file_sig(c["eo_path"])},
                           extra={"eo_sha256": rec["eo_sha256"], "calibration": rec})
            sp.write_stamp(outd, f"detect_{cam}", {"eo_file": sp.file_sig(c["eo_path"])},
                           extra={"setup_coef": sp.SETUP_COEF, "setup_fit": fit, "calibration": rec,
                                  "eo_sha256": rec["eo_sha256"]})
            (outd / "pointing").mkdir(exist_ok=True)
            (outd / "pointing" / f"pointing_{cam}.json").write_text(json.dumps(
                {"rows": [{"date": "2024-10-23", "verdict": "matches", "why": "", "horizon_offset_px": 2.0,
                           "d_azimuth_deg": 0.01, "d_tilt_deg": 0.0, "d_roll_deg": 0.0}], "info": {}}))
        sp.write_stamp(outd, "merge", {"after": {"detect_c1": sp.token(outd, "detect_c1"),
                                                 "detect_c2": sp.token(outd, "detect_c2"),
                                                 "forcing": sp.token(outd, "forcing")}},
                       extra={"setup": {c: {"coef": sp.SETUP_COEF, "fit": fit} for c in ("c1", "c2")},
                              "cameras": ["c1", "c2"]})
        sp.write_stamp(outd, "dem", {"after": {"merge": sp.token(outd, "merge")}},
                       extra={"settings": sp.dem_settings_requested(args_)})
        st_ = {"steps": {}, "failures": [], "not_applied": []}
        res_, out_ = quiet(sp.step_maps, plan_, args_, st_)
        for cam in ("c1", "c2"):
            sp.map_png(plan_, cam).write_bytes(b"png")          # <date>[_<variant>]_waterlines_<cam>.png
            sp.write_stamp(outd, f"map_{cam}", {"after": {"waterlines": sp.waterline_token(plan_)}})
        st_["steps"] = {s: "built" for s in sp.STEPS}
        res_c, out_c = quiet(sp.step_compare, plan_, args_, st_)
        st_["steps"]["compare"] = res_c
        prov_, _ = quiet(sp.write_provenance, plan_, args_, st_, True)
        return prov_, (outd / "README.txt").read_text(), res_c, out_c

    # -- a calibration made WITHOUT the held GCPs but NOT by the pre-registered procedure: synthetic files here,
    # turned from the 23 Jan 2025 solution (c2 through a Dec 2024 step), as estimate_eo_rotation.py writes them:
    # the reference(s) in a comment line ('X.yaml rotated by N deg'), an acceptance line added by hand. Their angles
    # are made up. Traced to the Jan GCPs, but never INDEPENDENT: a turn not made by the pinned procedure could
    # have been chosen after seeing calibration A (the procedure's own B: test_gcpfree_b).
    gf = d / "gcpfree"
    gf.mkdir()

    from georectify import cirn_angles_to_R

    def turn_deg(e1, e2):
        q = cirn_angles_to_R(*e1[3:6]) @ cirn_angles_to_R(*e2[3:6]).T
        return float(np.degrees(np.arccos(np.clip((np.trace(q) - 1.0) / 2.0, -1.0, 1.0))))

    def turned(path, ref_names, acceptance=None, ref_dir=cal, d_angles=(-11.0, 0.2, 0.1), shift=(0.0, 0.0, 0.0),
               note_turn=None):
        """A synthetic turned calibration: the first reference's position (+ shift), its angles + d_angles,
        and the notes estimate_eo_rotation.py writes, with the true turn from each reference (or note_turn)."""
        refs_e = [load_extrinsics(ref_dir / n) if (ref_dir / n).is_file() else load_extrinsics(gf / n)
                  for n in ref_names]
        ang = np.degrees(refs_e[0][3:6]) + np.array(d_angles)
        e_new = np.r_[refs_e[0][:3] + np.array(shift), np.radians(ang)]
        er.write_eo(path, e_new, ang,
                    "Solved by gcpfree_eo: " + "; ".join(
                        f"{n} rotated by {note_turn if note_turn is not None else turn_deg(e_new, r):.2f} deg from "
                        f"140 feature matches (6 new / 6 reference photos)" for n, r in zip(ref_names, refs_e))
                    + "; position unchanged.")
        if acceptance:
            with open(path, "a") as f:
                f.write(f"# acceptance: {acceptance}\n")
        return path
    b_eo = {"c1": turned(gf / "CACO03_c1_20241023_gcpfree_EO.yaml", ["CACO03_c1_20250123_EO.yaml"],
                         "PASS (inliers 140 >= 40, residual 0.8 px, horizon on the Oct photos 3 px)")}
    turned(gf / "CACO03_c2_20241205_gcpfree_EO.yaml", ["CACO03_c2_20250123_EO.yaml"])
    b_eo["c2"] = turned(gf / "CACO03_c2_20241023_gcpfree_EO.yaml", ["CACO03_c2_20241205_gcpfree_EO.yaml"],
                        "PASS-WIDER: the chain through Dec and the direct turn differ by 0.12 deg of tilt: +/-0.12 deg, "
                        "~0.2 m of DEM at 300 m")
    b_eo = {k: str(v) for k, v in b_eo.items()}
    ab = argparse.Namespace(**dict(vars(a), eo=b_eo, calibration_variant="turned",
                                   output_root=str(roots["out_turned"])))

    with traced():
        # 1. the plan: the date's survey rows withheld
        plan, out = quiet(sp.build_plan, "2024-10-23", cfg, srv, a)
        check(plan["blind"]["held_by"].startswith("Dr Chris Sherwood") and plan["surveys"] == []
              and [h["name"] for h in plan["held"]] == ["oct_gcps"] and "path" not in plan["held"][0]
              and "BLIND" in out,
              "2024-10-23 is BLIND: its survey rows are withheld from the plan (names, types, dates only), said")
        # 2. every look-up of a held survey is refused before any folder is searched
        for what, call in (("resolve_survey, no plan", lambda: sp.resolve_survey(oct_row, None, a)),
                           ("resolve_survey, the blind plan", lambda: sp.resolve_survey(oct_row, plan, a)),
                           ("survey_path_for", lambda: sp.survey_path_for(oct_row, plan, a)),
                           ("another date's row of the held survey's date",
                            lambda: sp.resolve_survey(dict(oct_row, date="2024-10-25", path=""), None, a)),
                           ("a row of the date with a path",
                            lambda: sp.resolve_survey(dict(oct_row, path=canary_names[0]), None, a)),
                           ("any survey row in a blind plan",
                            lambda: sp.resolve_survey([s for s in srv if s["name"] == "jan_lidar"][0], plan, a))):
            try:
                call()
                check(False, f"{what}: NOT refused")
            except sp.BlindSurvey as exc:
                check("held by Dr Chris Sherwood" in str(exc), f"{what}: refused (BlindSurvey, never swallowed by "
                                                               f"'except Exception')")
        check(not issubclass(sp.BlindSurvey, Exception), "BlindSurvey is not an Exception: no broad except swallows it")
        # 3. the compare step does not run; nothing is a missing part; no label is looked up
        state = {"steps": {}, "failures": [], "not_applied": []}
        res, out = quiet(sp.step_compare, plan, a, state)
        check(res == "not run (blind: DEM only)" and "NOT RUN" in out and state["comparisons"] == {},
              "compare: not run on a blind date")
        lw, why = sp.label_when_built("2024-10-23", oct_row, cfg, srv, a)
        check(lw == "blind (DEM only)", "--summary's label of a not-yet-built blind date: 'blind (DEM only)', never "
                                        "looked up")
        # 4. calibration A (the table's), a whole (fake) product: status complete, README and provenance blind
        prov, rd, res_c, _ = fake_product(plan, a)
        mp = sp.read_stamp(plan["out"], "map_plan") or {}
        check((plan["out"] / "maps" / "2024-10-23_waterlines_plan.png").exists() and mp["sig"]["surveys"] == []
              and mp["sig"].get("blind") == "2024-10-23",
              "maps: the plan view is drawn over no survey, its stamp says blind")
        check(prov["status"] == "complete" and prov["comparisons"] == [] and not prov["status_why"]
              and res_c == "not run (blind: DEM only)" and not (plan["out"] / "compare").exists(),
              "A: a blind product with everything else built: status complete (not partial), no comparison made")
        stmt = ("blind: the survey is held by Dr Chris Sherwood (USGS) for an independent check; nothing in this "
                "build read it")
        check(prov["blind"]["statement"] == stmt and stmt in rd and "BLIND: NO COMPARISON" in rd
              and "compare/      none: a blind product" in rd and "DEM - survey grids" not in rd,
              "README and provenance: 'blind: the survey is held by Dr Chris Sherwood ...; nothing in this build read it'")
        ck = prov["blind"]["check"]
        check(ck["label"] == "PARTLY-CIRCULAR"
              and all("not the geometry (PARTLY-CIRCULAR for the geometry)" in v["why"] for v in ck["per_camera"].values())
              and "CALIBRATION: the table's (surveys/survey_dates.csv eo_file): c1 CACO03_c1_20241023_EO.yaml" in rd
              and "c1: PARTLY-CIRCULAR: CACO03_c1_20241023_EO.yaml solved from the GCPs surveyed 23 Oct 2024" in rd
              and "tests the water level, the detection and the wave setup, not the geometry" in rd,
              "A: solved from the held GCPs: his check is PARTLY-CIRCULAR for the geometry (water level, detection, "
              "setup only), said in the README")
        check("calibration CACO03_c1_20241023_EO.yaml: solved from the GCPs surveyed 23 Oct 2024" in rd
              and "999" not in rd, "A: WHAT WAS USED says what the calibration was solved from")
        # 5. --summary lists it as 'blind (DEM only)', also before it is built
        rc, out = quiet(sp.summary, a)
        sm = (roots["out"] / "summary.txt").read_text()
        row = [r for r in csv.DictReader(open(roots["out"] / "summary.csv")) if r["date"] == "2024-10-23"]
        check(rc == 0 and "Blind (DEM only)" in sm and "2024-10-23 oct_gcps: blind (DEM only)" in sm
              and row and row[0]["label"] == "blind (DEM only)" and row[0]["held_check_label"] == "PARTLY-CIRCULAR"
              and row[0]["calibration_variant"] == "table",
              "--summary: 'blind (DEM only)', with the calibration and what the holder's check can test")
        on = Path(tabs["config"]).with_name("survey_dates_on.csv")
        on.write_text(FIXTURE_DATES.replace("2024-10-23,c1,0,", "2024-10-23,c1,1,").replace("2024-10-23,c2,0,",
                                                                                            "2024-10-23,c2,1,"))
        rc, out = quiet(sp.summary, argparse.Namespace(**dict(vars(a), output_root=str(roots["fresh"]), config=str(on))))
        sm = (roots["fresh"] / "summary.txt").read_text()
        check("2024-10-23 oct_gcps: blind (DEM only): not built yet" in sm,
              "--summary before it is built: the enabled blind date is listed 'blind (DEM only): not built yet'")
        # 6. a dry run of the blind date says so and looks for nothing
        _, out = quiet(sp.dry_run, plan, a, list(sp.STEPS))
        check("BLIND" in out and "survey oct_gcps" not in out and "held check" in out,
              "--dry-run: BLIND, the survey not listed, what the check can test per camera")
        # 7. a turned calibration in its own root: what it was solved from, NOT independent, its acceptance, position
        pb, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, ab)
        check(all(c["eo_file"] == f"CACO03_{k}_20241023_gcpfree_EO.yaml" and c["variant"] == "turned"
                  and c["table_eo_file"] == f"CACO03_{k}_20241023_EO.yaml" for k, c in pb["cams"].items())
              and sp.option_refusal(pb, ab) == "",
              "--eo: the plan builds with the given calibrations, named variant, the table's kept for the record")
        ag = argparse.Namespace(**dict(vars(ab), calibration_variant="gcpfree", output_root=str(d / "x_gcpfree")))
        pg, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, ag)
        check("these are not (or not yet) the B files it delivered here" in sp.option_refusal(pg, ag),
              "... and under the procedure's own name (gcpfree): REFUSED (only its B is called gcpfree)")
        watch_a["on"] = True
        try:
            provb, rdb, _, _ = fake_product(pb, ab)
        finally:
            watch_a["on"] = False
        if hasattr(sys, "addaudithook"):
            check(not watch_a["seen"], "B's build never opens calibration A (the table's, solved from the held GCPs): "
                                       + (", ".join(watch_a["seen"]) or "none opened"))
        ckb = provb["blind"]["check"]
        c1o = provb["cameras"]["c1"]["calibration"]["origin"]
        c2o = provb["cameras"]["c2"]["calibration"]["origin"]
        check(provb["status"] == "complete" and ckb["label"] == "PARTLY-CIRCULAR"
              and all("NOT the pre-registered procedure" in v["why"] for v in ckb["per_camera"].values())
              and c1o["kind"] == "rotated" and c1o["gcp_dates"] == ["2025-01-23"]
              and c1o["position_from"] == "CACO03_c1_20250123_EO.yaml"
              and c2o["chain"] == ["CACO03_c2_20241023_gcpfree_EO.yaml", "CACO03_c2_20241205_gcpfree_EO.yaml",
                                   "CACO03_c2_20250123_EO.yaml"] and c2o["gcp_dates"] == ["2025-01-23"],
              "turned by hand: traced to the 23 Jan 2025 GCPs (c2 through its Dec 2024 step), position of the Jan "
              "solution, but NOT the pre-registered procedure: PARTLY-CIRCULAR, never INDEPENDENT")
        check("CALIBRATION: variant turned (--eo), NOT the table's" in rdb
              and "made WITHOUT ground control of its own: CACO03_c1_20250123_EO.yaml (solved from the GCPs surveyed "
                  "23 Jan 2025" in rdb
              and "NOT the pre-registered procedure (not made by gcpfree_eo.py" in rdb and "(INDEPENDENT" not in rdb
              and "pre-registered acceptance: PASS (inliers 140" in rdb
              and "pre-registered acceptance: PASS-WIDER (the chain through Dec" in rdb,
              "its README: the variant, what it was solved from, NOT the pre-registered procedure, its acceptance")
        cav_b = [c for c in provb["caveats"] if c.startswith("c2 calibration CACO03_c2_20241023_gcpfree_EO.yaml")]
        check(cav_b and "camera position (x, y, z) is that of CACO03_c2_20250123_EO.yaml" in cav_b[0]
              and "2024-08-27 (CACO03_c2_20240827_EO.yaml), 2025-01-23 (CACO03_c2_20250123_EO.yaml) put this camera "
                  "0.59 m apart in height" in cav_b[0]
              and "nothing of it is used by this build" in cav_b[0] and "LARGER UNCERTAINTY" in cav_b[0]
              and "999" not in rdb and "999" not in json.dumps(provb),
              "B: its camera position is the Jan solution's, with the mount's GCP solutions' spread (Aug and Jan: "
              "0.59 m in height for c2; the held survey's solution never read)")
        _, out = quiet(sp.dry_run, pb, ab, list(sp.STEPS))
        check("held check" in out and "PARTLY-CIRCULAR" in out and "acceptance" in out and "REFUSED" not in out
              and "NOT independent of the held survey" in out,
              "--dry-run: what it was solved from, PARTLY-CIRCULAR (warned), its acceptance, not refused")
        pf, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, argparse.Namespace(**dict(vars(a),
                                                                                      output_root=str(roots["fresh"]))))
        pf["out"].mkdir(parents=True)
        sp.mark_in_progress(pf, argparse.Namespace(**dict(vars(a), output_root=str(roots["fresh"]))))
        rc, out = quiet(sp.summary, argparse.Namespace(**dict(vars(a), output_root=str(roots["fresh"]))))
        sm = (roots["fresh"] / "summary.txt").read_text()
        check(sp.load_json(pf["out"] / "provenance.json")["blind"]["statement"].startswith("blind: the survey is held")
              and "2024-10-23 oct_gcps: blind (DEM only)" in sm and "status in progress" in sm,
              "a blind build in progress (or interrupted): provenance says blind, --summary lists it so, in progress")
        shutil.rmtree(str(pf["out"]))
        rc, out = quiet(sp.summary, ab)
        rowb = [r for r in csv.DictReader(open(roots["out_turned"] / "summary.csv")) if r["date"] == "2024-10-23"]
        check(rowb and rowb[0]["calibration_variant"] == "turned" and rowb[0]["held_check_label"] == "PARTLY-CIRCULAR",
              "--summary of its root: blind (DEM only), calibration turned, PARTLY-CIRCULAR")
    check(not seen, "the held survey was never looked for or opened (no glob of its date, no stat or open of a file "
                    "named like it, no comparison run): " + ("; ".join(seen[:5]) or "none"))

    # -- the runtime guard: whatever tries to read a file named like the held survey is refused
    if sp.arm_blind_guard():
        try:
            with open(sorted(canaries)[0]) as f_:
                f_.read()
            check(False, "a canary read: NOT refused")
        except sp.BlindSurvey as exc:
            check("named like a survey of 2024-10-23" in str(exc), "the runtime guard refuses to open a file named "
                                                                   "like the held survey (any reader, any folder)")
        for ok_name in ("CACO03_c1_20241023_gcpfree_EO.yaml", "2024-10-23_info.json", "2024-10-23_dem.asc"):
            check(sp.held_file_name(d / ok_name) is None, f"... never a calibration or a product ({ok_name})")
        sp.BLIND_DATES.pop("2024-10-23")
        try:
            check(open(sorted(canaries)[0]).read().startswith("1,"), "... and only while the date is blind")
        finally:
            sp.BLIND_DATES["2024-10-23"] = held
    else:
        skip("runtime guard", "no sys.addaudithook (Python < 3.8)")

    # -- refusals, before anything is written
    outd_a = plan["out"]
    (outd_a / "compare").mkdir()
    (outd_a / "compare" / "2024-10-23_oct_gcps_comparison.json").write_text("{}")
    check(any("compare" in x for x in sp.blind_leftovers(plan)), "earlier comparisons in the folder: listed")
    res, out = quiet(sp.build_date, "2024-10-23", cfg, srv, a, list(sp.STEPS))
    check(res == "refused" and "fresh --output-root" in out, "... and a blind build over them is REFUSED (exit 5)")
    shutil.rmtree(str(outd_a / "compare"))

    def refused(args_, words, what):
        root_ = Path(args_.output_root)
        before = sorted(p_.name for p_ in root_.iterdir()) if root_.exists() else []
        res_, out_ = quiet(sp.build_date, "2024-10-23", cfg, srv, args_, list(sp.STEPS))
        after = sorted(p_.name for p_ in root_.iterdir()) if root_.exists() else []
        check(res_ == "refused" and all(w in out_ for w in words) and before == after,
              f"{what}: REFUSED (exit 5), nothing written")
    refused(argparse.Namespace(**dict(vars(ab), output_root=str(roots["out"]))), ["own root", "out_turned"],
            "a variant into A's root (the default root's name)")
    refused(argparse.Namespace(**dict(vars(a), output_root=str(roots["out_turned"]))),
            ["another calibration variant", "turned (2024-10-23)"], "A (the table's) into a variant's root")
    refused(argparse.Namespace(**dict(vars(ab), output_root=str(roots["fresh_turned"]),
                                      eo=dict(b_eo, c2=str(gf / "missing.yaml")))),
            ["no such file", "never used in its place"], "an --eo file that is missing (never the table's instead)")
    p_m, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv,
                   argparse.Namespace(**dict(vars(ab), eo=dict(b_eo, c2=str(gf / "missing.yaml")))))
    check(p_m["cams"]["c2"]["eo_path"] is None and p_m["cams"]["c2"]["table_eo_path"] is not None,
          "... the plan has no calibration for it, not the table's")
    turned(gf / "failed.yaml", ["CACO03_c1_20250123_EO.yaml"], "FAIL: c1 and c2 disagree by 0.4 deg of roll")
    refused(argparse.Namespace(**dict(vars(ab), output_root=str(roots["fresh_turned"]),
                                      eo=dict(b_eo, c1=str(gf / "failed.yaml")))),
            ["acceptance of the --eo calibration FAILED", "0.4 deg of roll", "never built in its place"],
            "an --eo calibration whose pre-registered acceptance FAILED (B not delivered)")
    turned(gf / "unstated.yaml", ["CACO03_c1_20250123_EO.yaml"])
    refused(argparse.Namespace(**dict(vars(ab), output_root=str(roots["fresh_turned"]),
                                      eo=dict(b_eo, c1=str(gf / "unstated.yaml")))),
            ["must state its pre-registered acceptance", "unstated.yaml"],
            "on a blind date, an --eo calibration that states no acceptance (its checks not run)")
    check(sp.root_for_variant("/x/survey_products", "gcpfree") == Path("/x/survey_products_gcpfree")
          and sp.root_for_variant("/x/survey_products_gcpfree", None, "gcpfree") == Path("/x/survey_products")
          and sp.root_for_variant("/x/survey_products_c0", "gcpfree") == Path("/x/survey_products_gcpfree_c0"),
          "a variant's root: <root>_<variant> (before any _c<C>); the table's: the base root")
    # the same B root, B's file edited since: another calibration by content, refused without --replace
    b2 = dict(b_eo, c1=str(turned(gf / "CACO03_c1_20241023_gcpfree-v2_EO.yaml", ["CACO03_c1_20250123_EO.yaml"],
                                   "PASS")))
    ab2 = argparse.Namespace(**dict(vars(ab), eo=b2))
    pb2, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, ab2)
    note = sp.calibration_replace_note(pb2, ab2, list(sp.STEPS))
    check(note and "built with another calibration" in note and "--replace" in note,
          "B's root, another B (by content): a rebuild with it is refused without --replace")
    res, out = quiet(sp.build_date, "2024-10-23", cfg, srv, ab2, list(sp.STEPS))
    check(res == "refused" and "built with another calibration" in out, "... build_date refuses it (exit 5)")
    st0 = {"steps": {}, "failures": [], "not_applied": []}
    pa2, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, ab2)
    _, out = quiet(sp.keep_built_calibration, pa2, ab2, ["dem", "maps"], st0)
    check("NOT applied" in out and pa2["cams"]["c1"]["eo_file"] == "CACO03_c1_20241023_gcpfree_EO.yaml",
          "--steps dem,maps with another --eo: NOT APPLIED; the product keeps (and reports) the calibration that "
          "built it")
    words = sp.state_argv(ab, "2024-10-23")
    check("--eo" in words and f"c1={b_eo['c1']}" in words and words[words.index("--calibration-variant") + 1] == "turned",
          "the 'rebuild everything' line carries --eo and --calibration-variant")

    # -- origins: other references, unknown ones, one turned from the held survey's solution
    both = turned(gf / "both.yaml", ["CACO03_c1_20250123_EO.yaml", "CACO03_c1_20240827_EO.yaml"], "PASS")
    o = sp.calibration_origin(both, [str(cal)])
    check(o["kind"] == "rotated" and o["gcp_dates"] == ["2024-08-27", "2025-01-23"]
          and o["position_from"] == "CACO03_c1_20250123_EO.yaml"
          and sp.held_check(o, sp.blind_of("2024-10-23"))[0] == "PARTLY-CIRCULAR",
          "photos of Aug 2024 as a second reference: both traced, position from the first named (by hand: not "
          "INDEPENDENT)")
    moved = turned(gf / "moved.yaml", ["CACO03_c1_20250123_EO.yaml"], "PASS", shift=(0.0, 0.0, 0.112))
    o = sp.calibration_origin(moved, [str(cal)])
    check(o["kind"] == "unknown" and "comes from elsewhere" in o["text"]
          and sp.held_check(o, sp.blind_of("2024-10-23"))[0] == "PARTLY-CIRCULAR",
          "a 'turned' file whose camera position is not its reference's (e.g. taken from another solution): NOT KNOWN, "
          "never called independent")
    lie = turned(gf / "lie.yaml", ["CACO03_c1_20250123_EO.yaml"], "PASS", note_turn=3.0)
    o = sp.calibration_origin(lie, [str(cal)])
    check(o["kind"] == "unknown" and "but it is turned" in o["text"],
          "a file whose notes give another turn than its numbers: NOT KNOWN")
    from_a = turned(gf / "turned_from_oct.yaml", ["CACO03_c1_20241023_EO.yaml"])
    check(sp.held_check(sp.calibration_origin(from_a, [str(cal)]), sp.blind_of("2024-10-23"))[0] == "PARTLY-CIRCULAR",
          "one turned from the held survey's solution: PARTLY-CIRCULAR")
    (gf / "noted_nothing.yaml").write_text((cal / "CACO03_c1_20250123_EO.yaml").read_text())
    (gf / "CACO03_c1_20240827_EO.yaml").write_text((cal / "CACO03_c1_20240827_EO.yaml").read_text() + "# edited\n")
    lost = d / "lost"
    lost.mkdir()
    shutil.copy2(str(b_eo["c2"]), str(lost / "CACO03_c2_20241023_gcpfree_EO.yaml"))   # its Dec step not beside it
    for p_, what in ((gf / "noted_nothing.yaml", "no reference in its notes"),
                     (gf / "CACO03_c1_20240827_EO.yaml", "an archived name with other content"),
                     (lost / "CACO03_c2_20241023_gcpfree_EO.yaml", "a chain whose Dec step is not found")):
        o = sp.calibration_origin(p_, [])
        check(o["kind"] == "unknown" and sp.held_check(o, sp.blind_of("2024-10-23"))[0] == "PARTLY-CIRCULAR",
              f"a calibration of unknown origin ({what}): never called independent")
    (gf / "CACO03_c1_20240827_EO.yaml").unlink()

    # -- a file in the GCP-free procedure's format, written by hand (no report beside it, its 'Made by' sha256 not
    # the pinned one): its verdict in the first comment line, its reference with that file's sha256, the position
    # it keeps, the pan change and the rotation source. Traced and its claims checked as before, but NOT the
    # pre-registered procedure: never INDEPENDENT, and never built under the procedure's name (gcpfree). The
    # angles are made up (from the Jan 2025 solution: calibration A is never read for them). The procedure's own
    # files: test_gcpfree_b.
    def procedure_file(path, cam, verdict="DELIVERED", reasons="all checks passed", d_angles=(-12.5, 0.15, -0.1),
                       numbers_from=None, notes_ref=None, sha=None, pan=None, source="direct Jan->Oct"):
        numbers_from = numbers_from or f"CACO03_{cam}_20250123_EO.yaml"    # the solution its numbers come from
        notes_ref = notes_ref or numbers_from                              # the one its notes name
        e_num, e_note = load_extrinsics(cal / numbers_from), load_extrinsics(cal / notes_ref)
        ang = np.degrees(e_num[3:6]) + np.array(d_angles)
        pan_ = ((ang[0] - np.degrees(e_note[3]) + 180.0) % 360.0) - 180.0 if pan is None else pan
        notes = [f"CALIBRATION B (GCP-FREE) for the 23 Oct 2024 DEM, camera {cam}: {verdict} ({reasons})",
                 f"No Oct 2024 GCP used: CACO03_{cam}_20241023_EO.yaml (calibration A, solved from them) was not read.",
                 f"Reference EO: {notes_ref} (Jan 2025 GCPs) sha256 {sha or sp.file_sha(cal / notes_ref)}; "
                 f"lens CACO03_{cam}_20240801_IO.yaml.",
                 f"Position (x, y, z) kept from {notes_ref}: B inherits the Jan 2025 solution's camera position.",
                 f"Rotation source: {source}; pan change {pan_:+.3f} deg from {notes_ref}.",
                 "Solve (Jan->Oct): 212 inliers of 1830 matches, residual median 0.61 px, inlier days [20, 21, 22, 24, "
                 "25], leave-one-day-out 1.9 px.",
                 "Per-day: 5 days measured, off days [].",
                 "Horizon (differential, vs the reference EO on the reference photos): tilt +0.031, roll -0.012 deg.",
                 "Oct photos: 1729782000.Thu.Oct.24_15_00_00.GMT.2024.CACO03.c1.timex.jpg",
                 "Reference photos: 1737471600.Tue.Jan.21_15_00_00.GMT.2025.CACO03.c1.timex.jpg",
                 f"Made by gcpfree_eo.py (sha256 {'0' * 64}) on 2026-10-14 01:00 UTC; report "
                 f"CACO03_{cam}_gcpfree-oct2024_report.json"]
        path.write_text(f"x: {e_num[0]:.4f}\ny: {e_num[1]:.4f}\nz: {e_num[2]:.4f}\n"
                        f"azimuth: {ang[0]:.4f}\ntilt: {ang[1]:.4f}\nroll: {ang[2]:.4f}\n"
                        + "".join(f"# {n}\n" for n in notes)
                        + "#x, y - camera Easting / Northing (m, UTM zone 19); z - camera elevation (m, NAVD88)\n"
                          "#azimuth, tilt, roll - degrees, CIRN convention\n")
        return path
    pe = {"c1": str(procedure_file(cal / "CACO03_c1_20241023_gcpfree_EO.yaml", "c1")),
          "c2": str(procedure_file(cal / "CACO03_c2_20241023_gcpfree_EO.yaml", "c2", "DELIVERED-WIDER",
                                   "chain differs from direct by 12.3 px (0.150 deg; > 10.0 px): uncertainty at least "
                                   "that", d_angles=(-13.4, 0.1, 0.05), source="direct Jan->Oct"))}
    ap_ = argparse.Namespace(**dict(vars(a), eo=pe, calibration_variant="proc",
                                    output_root=str(roots["proc_proc"])))
    with traced():
        pp, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, ap_)
        check(sp.option_refusal(pp, ap_) == ""
              and all("NOT independent" in x for x in sp.blind_independence_warnings(pp, ap_))
              and len(sp.blind_independence_warnings(pp, ap_)) == 2,
              "files in gcpfree_eo.py's format written by hand (DELIVERED, DELIVERED-WIDER): built as another variant, "
              "each WARNED as not independent")
        watch_a["on"] = True
        try:
            provp, rdp, _, _ = fake_product(pp, ap_)
        finally:
            watch_a["on"] = False
    check(not seen, "... and the held survey was never looked for: " + ("; ".join(seen[:5]) or "none"))
    if hasattr(sys, "addaudithook"):
        check(not watch_a["seen"], "... nor calibration A opened: " + (", ".join(watch_a["seen"]) or "none opened"))
    o1, o2 = (provp["cameras"][k]["calibration"]["origin"] for k in ("c1", "c2"))
    check(provp["status"] == "complete" and provp["blind"]["check"]["label"] == "PARTLY-CIRCULAR"
          and all(o["kind"] == "rotated" and o["gcp_dates"] == ["2025-01-23"] and o["rotation_source"] == "direct Jan->Oct"
                  and o["procedure"]["format"] and not o["procedure"]["pinned"]
                  and any("not the pre-registered one" in w for w in o["procedure"]["why"]) for o in (o1, o2))
          and o1["position_from"] == "CACO03_c1_20250123_EO.yaml" and o2["position_from"] == "CACO03_c2_20250123_EO.yaml"
          and o1["acceptance"]["verdict"] == "PASS" and o2["acceptance"]["verdict"] == "PASS-WIDER",
          "gcpfree_eo.py's format by hand: traced to the 23 Jan 2025 GCPs (reference sha256, pan change, position "
          "checked), DELIVERED = PASS, DELIVERED-WIDER = PASS-WIDER, but its 'Made by' sha256 is not the pinned "
          "procedure's: PARTLY-CIRCULAR")
    check("CACO03_c1_20250123_EO.yaml (solved from the GCPs surveyed 23 Jan 2025 (not the Jan 2025 lidar)) turned "
          "-12.500 deg in pan to this pointing by matching photos (rotation source: direct Jan->Oct), its camera "
          "position (x, y, z) that of CACO03_c1_20250123_EO.yaml" in rdp
          and "pre-registered acceptance: PASS (gcpfree_eo.py: DELIVERED; all checks passed)" in rdp
          and "pre-registered acceptance: PASS-WIDER (gcpfree_eo.py: DELIVERED-WIDER; chain differs from direct by "
              "12.3 px" in rdp and "NOT the pre-registered procedure" in rdp and "(INDEPENDENT" not in rdp
          and "999" not in rdp,
          "... its README: what it was solved from, its turn, its position, the verdict word for word, NOT the "
          "pre-registered procedure")
    cav_p = [c for c in provp["caveats"] if c.startswith("c2 calibration CACO03_c2_20241023_gcpfree_EO.yaml")]
    check(cav_p and "LARGER UNCERTAINTY" in cav_p[0] and "chain differs from direct by 12.3 px" in cav_p[0]
          and "0.59 m apart in height" in cav_p[0],
          "... DELIVERED-WIDER: the larger uncertainty, as stated, and the mount's position spread under KNOWN CAVEATS")
    for n, kw, words, what in (
            ("nd", dict(verdict="NOT DELIVERED", reasons="c1 and c2 pan changes differ by 4.10 deg (> 3.0)"), None,
             "NOT DELIVERED"),
            ("sha", dict(sha="0" * 64), "not the one it was turned from", "a reference whose sha256 is not the stated one"),
            ("pan", dict(pan=2.0), "pan change of +2.000 deg", "a pan change its numbers do not have"),
            ("aug", dict(numbers_from="CACO03_c1_20240827_EO.yaml", notes_ref="CACO03_c1_20250123_EO.yaml",
                         source="Aug->Oct"), "comes from elsewhere",
             "numbers from the Aug 2024 solution, notes naming the Jan one")):
        f_ = procedure_file(cal / f"CACO03_c1_20241023_gcpfree-{n}_EO.yaml", "c1", **kw)
        o = sp.calibration_origin(f_, [str(cal)])
        if words is None:
            check(o["acceptance"]["verdict"] == "FAIL", f"gcpfree_eo.py's {what}: read as FAIL")
            refused(argparse.Namespace(**dict(vars(ap_), output_root=str(roots["fresh_proc"]), eo=dict(pe, c1=str(f_)))),
                    ["acceptance of the --eo calibration FAILED", "NOT DELIVERED", "pan changes differ"],
                    f"gcpfree_eo.py's {what}")
            continue
        check(o["kind"] == "unknown" and words in o["text"]
              and sp.held_check(o, sp.blind_of("2024-10-23"))[0] == "PARTLY-CIRCULAR",
              f"gcpfree_eo.py-style file, {what}: NOT KNOWN, never called independent")
    pa_ = argparse.Namespace(**dict(vars(ap_), eo=dict(pe, c1=str(cal / "CACO03_c1_20241023_gcpfree-aug_EO.yaml"))))
    pq, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, pa_)
    check(any("NOT independent" in x for x in sp.blind_independence_warnings(pq, pa_)),
          "... a WARNING before such a build (its README would call the holder's check PARTLY-CIRCULAR)")
    for n in ("CACO03_c1_20241023_gcpfree_report.txt", "CACO03_c1_20241023_gcpfree_report.json"):
        check(sp.held_file_name(cal / n) is None, f"the runtime guard never refuses the procedure's report ({n})")
    check(sp.held_file_name(cal / "2024-10-23_gcpfree_GCP_xyz.csv") is not None,
          "... but a survey-like name with 'gcpfree' in it is still refused")

    # -- freshness: the calibration's content is in every pointing/detection stamp
    e = d / "eo_edit.yaml"
    e.write_text((cal / "CACO03_c1_20250123_EO.yaml").read_text())
    sp.write_stamp(d, "pointing_c1", {"x": 1}, extra={"eo_sha256": sp.file_sha(e)})
    check(sp.eo_content_fresh(d, "pointing_c1", e) == (True, ""), "the same calibration content: fresh")
    e.write_text(e.read_text().replace("azimuth", "azimuth ", 1))
    sp._FILE_SHA.clear()
    ok, why = sp.eo_content_fresh(d, "pointing_c1", e)
    check(not ok and "content" in why, "an edited calibration (same name, same size or not): stale, rebuilt")
    sp.write_stamp(d, "detect_c1", {"x": 1})
    tk = sp.token(d, "detect_c1")
    check(sp.eo_content_fresh(d, "detect_c1", e) == (True, "") and sp.read_stamp(d, "detect_c1")["eo_sha256"]
          == sp.file_sha(e) and sp.token(d, "detect_c1") == tk,
          "a stamp from before Oct 2026 (no content): its content recorded now, its token kept (nothing rebuilds)")
    # -- the GeoTIFFs carry the surveys' CRS
    check(sp.GEOTIFF_EPSG == 6348 and "--epsg 6348" in sp.cmd_text(sp.dem_cmds(plan, a)[2]),
          "GeoTIFFs tagged EPSG:6348 by default (NAD83(2011) / UTM 19N, the surveys' frame)")
    # -- the command line: --eo needs --date and a variant name
    for argv, what in ((["--all", "--eo", "c1=x.yaml", "--calibration-variant", "gcpfree"], "--eo with --all"),
                       (["--date", "2024-10-23", "--eo", "c1=x.yaml"], "--eo without --calibration-variant"),
                       (["--date", "2024-10-23", "--eo", "c3=x.yaml", "--calibration-variant", "gcpfree"], "camera c3"),
                       (["--date", "2024-10-23", "--eo", "c1=x.yaml", "--calibration-variant", "GCP free"],
                        "a variant name with capitals and a space")):
        r = subprocess.run([sys.executable, str(HERE / "survey_products.py")] + argv + ["--dry-run", "--output-root",
                                                                                       str(d / "usage")],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, cwd=str(d))
        check(r.returncode == sp.EXIT_CODES["usage"] and not (d / "usage").exists(),
              f"{what}: a usage error (exit {r.returncode}), nothing written")
    # -- the station's own tables: 2024-10-23 is blind, built when asked
    scfg, ssrv = sp.read_table(sp.DEFAULT_CONFIG, sp.CONFIG_FIELDS), sp.read_table(sp.DEFAULT_SURVEYS, sp.SURVEY_FIELDS)
    p_, _ = quiet(sp.build_plan, "2024-10-23", scfg, ssrv, cfg_args())
    check(p_ is not None and p_.get("blind") and p_["surveys"] == [] and not p_["forced_disabled"],
          "the station's tables: 2024-10-23 enabled and BLIND (DEM only), no --force-disabled needed")
    check(all(not s["path"] for s in ssrv if s["date"] == "2024-10-23") and
          all("BLIND" in s["why"] and "Sherwood" in s["why"] for s in ssrv if s["date"] == "2024-10-23"),
          "surveys.csv: the held survey has no path, and its row says BLIND, held by Chris Sherwood")


def tiff_description_and_epsg(path):
    """(ImageDescription, ProjectedCSType) of a baseline little-endian GeoTIFF as asc_to_geotiff.py writes it."""
    import struct
    b = Path(path).read_bytes()
    off = struct.unpack("<I", b[4:8])[0]
    n = struct.unpack("<H", b[off:off + 2])[0]
    desc, epsg = None, None
    for k in range(n):
        tag, typ, cnt, val = struct.unpack("<HHI4s", b[off + 2 + 12 * k: off + 14 + 12 * k])
        data_off = struct.unpack("<I", val)[0]
        if tag == 270:
            raw = b[data_off:data_off + cnt] if cnt > 4 else val[:cnt]
            desc = raw.rstrip(b"\0").decode("ascii")
        if tag == 34735:
            keys = struct.unpack("<" + "H" * cnt, b[data_off:data_off + 2 * cnt])
            for i in range(4, len(keys), 4):
                if keys[i] == 3072:
                    epsg = keys[i + 3]
    return desc, epsg


def fake_built_product(plan_, args_):
    """Stamps of every step, waterlines, maps and DEM files of a build (no detection is run). -> (prov, README)."""
    outd = plan_["out"]
    wdir = outd / "waterlines"
    wdir.mkdir(parents=True)
    rows = ["source_file,camera,capture_time_utc,easting_utm19,northing_utm19,tide_elevation_navd88,"
            "beach_elevation_navd88,setup_correction_m"]
    for k in range(40):
        for cam in plan_["cams"]:
            rows.append(f"f{k}_{cam},{cam},2024-10-2{k % 7} 15:00:00,{420150 + k},{4638350 + 3 * k},"
                        f"{0.1 * (k % 9) - 0.3:.2f},{0.1 * (k % 9) - 0.1:.2f},0.2")
    (wdir / "contour_points_ground.csv").write_text("\n".join(rows) + "\n")
    sp.write_stamp(outd, "forcing", {"x": 1})
    fit = sp.setup_fit_info(sp.SETUP_COEF)
    for cam, c in plan_["cams"].items():
        rec = sp.calibration_record(c)
        sp.write_stamp(outd, f"pointing_{cam}", {"eo_file": sp.file_sig(c["eo_path"])},
                       extra={"eo_sha256": rec["eo_sha256"], "calibration": rec})
        sp.write_stamp(outd, f"detect_{cam}", {"eo_file": sp.file_sig(c["eo_path"]),
                                               "envelope": sp.file_sig(c["envelope_path"])},
                       extra={"setup_coef": sp.SETUP_COEF, "setup_fit": fit, "calibration": rec,
                              "eo_sha256": rec["eo_sha256"], "envelope_sha256": sp.file_sha(c["envelope_path"])})
        (outd / "pointing").mkdir(exist_ok=True)
        (outd / "pointing" / f"pointing_{cam}.json").write_text(json.dumps(
            {"rows": [{"date": "2024-10-23", "verdict": "matches", "why": "", "horizon_offset_px": 2.0,
                       "d_azimuth_deg": 0.01, "d_tilt_deg": 0.0, "d_roll_deg": 0.0}], "info": {}}))
    sp.write_stamp(outd, "merge", {"after": {**{f"detect_{c}": sp.token(outd, f"detect_{c}") for c in plan_["cams"]},
                                             "forcing": sp.token(outd, "forcing")}},
                   extra={"setup": {c: {"coef": sp.SETUP_COEF, "fit": fit} for c in plan_["cams"]},
                          "cameras": sorted(plan_["cams"])})
    stem = sp.dem_stem(plan_)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for suf in ("_dem.asc", "_dem.tif", "_spread.tif", "_count.tif", "_dem.png"):
        Path(f"{stem}{suf}").write_bytes(b"grid")
    sp.write_stamp(outd, "dem", {"after": {"merge": sp.token(outd, "merge")}},
                   extra={"settings": sp.dem_settings_requested(args_)})
    st_ = {"steps": {}, "failures": [], "not_applied": []}
    quiet(sp.step_maps, plan_, args_, st_)
    for cam in plan_["cams"]:
        sp.map_png(plan_, cam).write_bytes(b"png")
        sp.write_stamp(outd, f"map_{cam}", {"after": {"waterlines": sp.waterline_token(plan_)}})
    st_["steps"] = {s: "built" for s in sp.STEPS}
    res_c, _ = quiet(sp.step_compare, plan_, args_, st_)
    st_["steps"]["compare"] = res_c
    prov_, _ = quiet(sp.write_provenance, plan_, args_, st_, True)
    return prov_, (outd / "README.txt").read_text()


def git_in(repo, *a):
    return subprocess.run(["git", "-C", str(repo), "-c", "user.name=test", "-c", "user.email=test@example.org",
                           "-c", "commit.gpgsign=false"] + [str(x) for x in a],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)


# The photos a test pre-registration pins (gcpfree_eo.py --preregister scans the station's folders for them; a test
# has none): one Oct and one Jan photo of each camera, and a Dec photo of c2, made-up sizes and sha256s.
FAKE_PHOTOS = [["c1", "oct", "2024-10-23", "1729695600.Wed.Oct.23_15_00_00.GMT.2024.CACO03.c1.timex.jpg", 1200000,
                "1" * 64],
               ["c1", "jan", "2025-01-23", "1737644400.Thu.Jan.23_15_00_00.GMT.2025.CACO03.c1.timex.jpg", 1100000,
                "2" * 64],
               ["c2", "oct", "2024-10-23", "1729695600.Wed.Oct.23_15_00_00.GMT.2024.CACO03.c2.timex.jpg", 1300000,
                "3" * 64],
               ["c2", "jan", "2025-01-23", "1737644400.Thu.Jan.23_15_00_00.GMT.2025.CACO03.c2.timex.jpg", 1150000,
                "4" * 64],
               ["c2", "dec", "2024-12-12", "1734015600.Thu.Dec.12_15_00_00.GMT.2024.CACO03.c2.timex.jpg", 1250000,
                "5" * 64]]


class GcpfreeRepo:
    """A scratch git repository standing in for this one while a calibration B is checked: its calibration/gcpfree/
    is survey_products' GCPFREE_DIR and its top folder the repository the B's commits are checked in (while 'with'
    it). gcpfree_eo.py writes its pre-registration (committed at once); runs are written by gcpfree_eo.py's own
    finish_run from given results and committed in one commit, as on the station."""

    def __init__(self, root, rule="per-camera", push=False):
        import gcpfree_eo as G
        self.G, self.root, self.rule = G, Path(root), rule
        self.root.mkdir(parents=True)
        git_in(self.root, "init", "-q")
        self.dir = self.root / "calibration" / "gcpfree"
        self.prereg = G.write_preregistration(self.dir / "PREREGISTERED.txt", rule, "the owner (a test)", G.ROOTS,
                                              photos=FAKE_PHOTOS)
        self.commit("gcpfree: pre-registration")
        if push:
            self.push()
        self.info, probs = G.preregistration(self.prereg)
        # gcpfree_eo.py refuses a run until the commit is pushed; these runs are written by its finish_run directly
        # (no run is refused here), so an unpushed pre-registration stands for the B README's 'on NO remote branch'
        assert not [x for x in probs if push or "on no remote-tracking branch" not in x], probs
        self.saved = []

    def push(self):
        """The owner's `git push`: to a bare repository standing in for GitHub (origin/<branch> then holds it)."""
        bare = Path(str(self.root) + "_github.git")
        if not bare.exists():
            git_in(bare.parent, "init", "-q", "--bare", bare.name)
            git_in(self.root, "remote", "add", "origin", str(bare))
        return git_in(self.root, "-c", "push.negotiate=false", "push", "-q", "origin", "HEAD")

    def commit(self, msg):
        git_in(self.root, "add", "-A")
        return git_in(self.root, "commit", "-q", "-m", msg)

    def __enter__(self):
        self.saved.append((sp.GCPFREE_DIR, dict(sp._REPO_TOP)))
        sp.GCPFREE_DIR = self.dir
        sp._REPO_TOP["top"] = os.path.realpath(str(self.root))
        sp._FILE_SHA.clear()
        return self

    def __exit__(self, *exc):
        sp.GCPFREE_DIR, top = self.saved.pop()
        sp._REPO_TOP.clear()
        sp._REPO_TOP.update(top)
        sp._FILE_SHA.clear()

    def run(self, run_id, results, started, run_no=1, variant="gcpfree", first_run=None, aug_dir=None, eligible=True,
            commit=True, edit_info=None, folder=None, prereg=None):
        """A run of the procedure with these per-camera results (gcpfree_eo.run_camera's), written by its own
        finish_run. -> {camera: its B file's path (written only when delivered)}."""
        G = self.G
        folder = Path(folder or self.dir / "runs" / run_id)
        folder.mkdir(parents=True)
        info = {"run_id": run_id, "run_no": run_no if eligible else None, "eligible": eligible,
                "trial_why": None if eligible else ["a test"], "variant": variant if eligible else None,
                "first_run": first_run, "status": "running", "started_utc": started,
                "argv": ["gcpfree_eo.py"] + (["--aug-dir", aug_dir] if aug_dir else []),
                "argv_text": "gcpfree_eo.py" + (f" --aug-dir {aug_dir}" if aug_dir else ""), "roots": list(G.ROOTS),
                "camera": "both", "aug_dir": aug_dir, "procedure_sha256": G.own_sha256(),
                "modules_sha256": G.module_sha256(), "preregistration": prereg or self.info,
                "constants": G.constants(), "constants_sha256": G.constants_sha256(), "exposure": G.EXPOSURE,
                "solved": sorted(results), "kept": {}}
        used = [r for r in FAKE_PHOTOS if r[0] in results]
        info.update(photos_used=used, photos_used_sha256=G.photos_sha256(used), photos_preregistered=True)
        if edit_info:
            edit_info(info)
        G.apply_outcome(results, (prereg or self.info).get("outcome_if_one_camera"), eligible)
        with contextlib.redirect_stdout(io.StringIO()):
            G.finish_run(folder, info, results, ["a test run"], print)
        if commit:
            self.commit(f"gcpfree: run {run_id}")
        return {c: folder / G.OUT_EO.format(cam=c) for c in results}


def gcpfree_res(cam, d_angles, verdict="DELIVERED", ref="jan"):
    """A camera's result as gcpfree_eo.run_camera returns it: made-up angles turned from the Jan 2025 (or Aug 2024)
    solution, nothing solved here (the procedure's own tests do that); calibration A is never read."""
    import gcpfree_eo as G
    from georectify import load_extrinsics, load_intrinsics
    cal = HERE / "calibration"
    io_f, jan_f, ref_f = (G.IO_FILE.format(cam=cam), G.REF_EO_FILE["jan"].format(cam=cam),
                          G.REF_EO_FILE[ref].format(cam=cam))
    e, e_jan = load_extrinsics(cal / ref_f), load_extrinsics(cal / jan_f)
    io_ = load_intrinsics(cal / io_f)
    link = {"label": "direct Jan->Oct", "gate_F": verdict != "NOT DELIVERED", "gate_S": True, "has_B": True,
            "solve": {"inliers": 300, "matches": 2000, "resid_med_px": 0.6, "jackknife_max_px": 0.8,
                      "inlier_new_days": ["2024-10-21", "2024-10-23", "2024-10-25"], "why": ""},
            "per_day": {"days": [{"day": "2024-10-21", "inliers": 90, "diff_px": 0.7}], "measured": 3,
                        "off_days": [], "gate_S": True, "why": ""},
            "horizon": {"status": "pass", "why": "a test"}, "new_photos": ["oct.jpg"], "ref_photos": ["jan.jpg"]}
    res = {"camera": cam, "io": io_f, "io_sha256": sp.file_sha(cal / io_f), "reference_eo": jan_f,
           "reference_eo_sha256": sp.file_sha(cal / jan_f), "direct": link, "verdict": verdict,
           "source": None if verdict == "NOT DELIVERED" else "direct Jan->Oct", "notes": [],
           "reasons": ["direct: 14 inliers < 60 (a test)"] if verdict == "NOT DELIVERED" else []}
    if ref == "aug":
        res["aug"] = {"link": dict(link, label="Aug->Oct"), "reference_eo": ref_f,
                      "reference_eo_sha256": sp.file_sha(cal / ref_f), "has_B": True}
    if verdict != "NOT DELIVERED":
        B = e.copy()
        B[3:6] += np.radians(d_angles)
        res.update(b_reference_eo=ref_f, b_reference_eo_sha256=sp.file_sha(cal / ref_f), B=G.eo_deg(B),
                   pan_change_deg=round(float(((np.degrees(B[3] - e_jan[3]) + 180) % 360) - 180), 4),
                   source="Aug->Oct" if ref == "aug" else "direct Jan->Oct")
        res["uncertainty"] = G.stated_uncertainty(res, io_)
    return res


def reseal(run_dir, cam, edit_b=None, edit_report=None):
    """Edit a B file (and its report) and rewrite what is checked against it -- its report's b_file_sha256 and
    run.json's -- to agree: a consistent forgery. -> the B file."""
    import hashlib
    import gcpfree_eo as G
    b = Path(run_dir) / G.OUT_EO.format(cam=cam)
    if edit_b:
        b.write_text(edit_b(b.read_text()))
    sha = hashlib.sha256(b.read_bytes()).hexdigest()
    rp = Path(run_dir) / (G.OUT_REPORT.format(cam=cam) + ".json")
    js = json.loads(rp.read_text())
    js["b_file_sha256"] = sha
    if edit_report:
        edit_report(js)
    rp.write_text(json.dumps(js, indent=1) + "\n")
    rj = json.loads((Path(run_dir) / "run.json").read_text())
    rj["cameras"][cam]["b_file_sha256"] = sha
    (Path(run_dir) / "run.json").write_text(json.dumps(rj, indent=1) + "\n")
    sp._FILE_SHA.clear()
    return b


def shifted(key, delta, fmt="{:.4f}"):
    """An edit of a B file's number `key` (x, y, z, azimuth, tilt, roll) by delta."""
    def edit(t):
        v = float(re.search(rf"^{key}: (\S+)$", t, re.M).group(1))
        return re.sub(rf"^{key}: \S+$", f"{key}: {fmt.format(v + delta)}", t, flags=re.M)
    return edit


def test_gcpfree_b(d):
    print("calibration B of the pre-registered GCP-free procedure: pinned, checked against its report, made and "
          "committed here, one camera, its own names, shared apart from A")
    import tarfile
    import gcpfree_eo as G
    d = Path(d) / "gcpfree_b"
    d.mkdir()
    cfg, srv, tabs = fixture_tables(d)
    # the pin is the procedure as it is here (and, once committed, as committed)
    pin = sp.GCPFREE_PINS[-1]
    check(pin["gcpfree_eo.py"] == G.own_sha256() and pin["modules"] == G.module_sha256()
          and pin["constants_sha256"] == G.constants_sha256(),
          f"GCPFREE_PINS[-1] is gcpfree_eo.py ({G.own_sha256()[:12]}...), its modules and its pre-registered "
          f"constants as they are")
    r = subprocess.run(["git", "-C", str(HERE), "diff", "--quiet", "HEAD", "--", "gcpfree_eo.py", *sorted(G.MODULES)],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if r.returncode == 1:
        print("  (note: gcpfree_eo.py or a module differs from its last commit: the pin must match the COMMITTED files "
              "before any B is made)")
    elif r.returncode != 0:                 # 128/129: not a git checkout here (e.g. a copy made by git archive)
        print("  (note: git cannot read a checkout here: the committed state of the pinned files is not checked)")
    check(sp.GCPFREE_STATION_ROOTS == G.ROOTS and set(sp.GCPFREE_VARIANTS) == {G.B_VARIANT, G.B2_VARIANT}
          and list(sp.GCPFREE_VARIANTS) == [G.B_VARIANT, G.B2_VARIANT]
          and os.path.realpath(str(sp.GCPFREE_DIR)) == os.path.realpath(str(G.OUT_DIR))
          and sp.CALIBRATION_ORIGINS["CACO03_c1_20250123_EO.yaml"][0] == G.ARCHIVED_SHA256["CACO03_c1_20250123_EO.yaml"]
          and sp.CALIBRATION_ORIGINS["CACO03_c2_20240827_EO.yaml"][0] == G.ARCHIVED_SHA256["CACO03_c2_20240827_EO.yaml"],
          "survey_products and gcpfree_eo agree: the station's photo folders, the variant names, the gcpfree folder, the "
          "archived calibrations' sha256")
    held = sp.BLIND_DATES["2024-10-23"]
    forcing = {}
    for k in ("adcp", "adcp_navd88", "chatham", "wis", "ndbc"):
        forcing[k] = str(d / f"{k}.csv")
        Path(forcing[k]).write_text("")

    def full_args(**kw):
        """Every option a build reads (the station's defaults where nothing is looked up)."""
        base = dict(survey_dirs=[str(d / "no_surveys")], calibration=str(HERE / "calibration"), force_disabled=True,
                    photo_roots=[str(d / "no_photos")], live_contours="/x/contour_points_timex.csv",
                    setup_fitted_to=None, geotiff_epsg=sp.GEOTIFF_EPSG, gnssr_spline=str(d / "no_spline.txt"),
                    gauge_csv=str(d / "no_gauge.csv"), live_cron=str(sp.LIVE_CRON), pointing_log_dir=str(d / "no_log"),
                    pointing_frames=3, pointing_fail=sp.POINT_FAIL_DEG, dem_cell=sp.DEM_CELL,
                    dem_min_points=sp.DEM_MIN_POINTS, dem_max_spread=sp.DEM_MAX_SPREAD, dem_max_hs=sp.DEM_MAX_HS,
                    dem_max_day_offset=sp.DEM_MAX_DAY_OFFSET, keep_moved_days=False, no_download=True,
                    synthetic=None, replace=False, dry_run=False, force=False, eo=None, calibration_variant=None)
        base.update(forcing)
        base.update(tabs)
        base.update(kw)
        return cfg_args(**base)

    check(sp.reasons_text(["x is NOT committed", "y", "x is NOT committed", ""]) == "x is NOT committed; y",
          "a refusal gives each reason once (two cameras of one B share its build record)")
    look = {v: sp.variant_lookalike(v) for v in ("gcpfree-oct", "gcpfree2", "b", "b2", "b-oct", "gcpfree",
                                                 "gcpfree-aug", "turned", "bird", "handfit-2")}
    check(all(look[v] for v in ("gcpfree-oct", "gcpfree2", "b", "b2", "b-oct"))
          and not any(look[v] for v in ("gcpfree", "gcpfree-aug", "turned", "bird", "handfit-2"))
          and "reads as the pre-registered GCP-free procedure's calibration B" in look["gcpfree-oct"]
          and sp.procedure_variant_problem({"cams": {"c2": {"eo_given": True, "eo_path": "/x/hand_turned_EO.yaml"}}},
                                           SimpleNamespace(calibration_variant="gcpfree-oct")) == look["gcpfree-oct"],
          "a variant name that reads as B (gcpfree-oct, gcpfree2, b2, b-oct) is refused unless it is one of the "
          "procedure's own names (gcpfree, gcpfree-aug): no lookalike file or archive name")

    def share_args(rt, sdir):
        return SimpleNamespace(date="2024-10-23", output_root=str(d / rt), share_dir=str(d / sdir), replace=False,
                               calibration=str(HERE / "calibration"))

    repo = GcpfreeRepo(d / "repo")                      # stands for this repository: calibration/gcpfree/ ...
    with repo:
        runs = repo.dir / "runs"
        prereg_sha = repo.info["sha256"]
        # delivery run 1 (per-camera): c1 NOT DELIVERED, c2 DELIVERED; delivery run 2 (Aug photos) solves c1 and
        # keeps c2: a second product (gcpfree-aug), the first run's B stays the delivered one
        r1 = repo.run("20261014T010000Z", {"c1": gcpfree_res("c1", (0, 0, 0), "NOT DELIVERED"),
                                           "c2": gcpfree_res("c2", (-13.4, -0.5, 0.3))}, "2026-10-14 01:00:00")
        kept = dict(gcpfree_res("c2", (-13.4, -0.5, 0.3)), solved_in_run="20261014T010000Z")
        r2 = repo.run("20261014T020000Z", {"c1": gcpfree_res("c1", (-10.2, 0.6, -0.4), ref="aug"), "c2": kept},
                      "2026-10-14 02:00:00", run_no=2, variant="gcpfree-aug", first_run="20261014T010000Z",
                      aug_dir="/aug_photos")
        check(r1["c2"].exists() and not r1["c1"].exists() and r2["c1"].exists() and r2["c2"].exists(),
              "gcpfree_eo.py's own writer: run 1 (per-camera) a c2 B only; run 2 (Aug photos) c1, and c2 kept")
        pc = sp.procedure_check(r1["c2"])
        check(pc["format"] and pc["pinned"] and not pc["why"] and pc["variant"] == "gcpfree" and pc["run_no"] == 1
              and pc["station_roots"] and pc["outcome_rule"] == "per-camera" and pc["b_commit"]
              and pc["committed"]["commits"] == 1 and pc["delivery_runs_here"] == ["20261014T010000Z", "20261014T020000Z"]
              and pc["uncertainty"]["rotation_px"] == G.ROT_FLOOR_PX and "EXPOSURE" in pc["exposure"],
              "its B: made by the pinned procedure, its report beside it agrees, the first delivery run of this "
              "repository, under its committed pre-registration, committed once")
        pa2 = sp.procedure_check(r2["c1"])
        check(pa2["pinned"] and pa2["variant"] == "gcpfree-aug" and pa2["run_no"] == 2,
              "run 2's B (Aug photos, after a run 1 without them): the one allowed second run, variant gcpfree-aug"
              + ("" if pa2["pinned"] else f" -- NOT: {pa2['why']}"))
        org = sp.calibration_origin(r1["c2"], [str(HERE / "calibration")])
        lab, why = sp.held_check(org, sp.blind_of("2024-10-23"))
        txt = sp.procedure_text(pc)
        check(org["kind"] == "rotated" and org["gcp_dates"] == ["2025-01-23"] and lab == "INDEPENDENT"
              and "made by the pre-registered procedure (gcpfree_eo.py sha256" in why
              and "committed to git once" in txt and str(runs / "20261014T010000Z") in txt
              and str(repo.dir / "PREREGISTERED.txt") in txt,
              "the holder's check of a DEM built with it: INDEPENDENT; its run folder, pre-registration and commit "
              "named")

        # forgeries copied out of the runs folder: never INDEPENDENT
        def forged(name, edit_b=None, edit_report=None, drop=None):
            f = d / "forged" / name
            shutil.copytree(str(runs / "20261014T010000Z"), str(f))
            b = f / G.OUT_EO.format(cam="c2")
            if edit_b:
                b.write_text(edit_b(b.read_text()))
            rp = f / (G.OUT_REPORT.format(cam="c2") + ".json")
            if edit_report:
                js = json.loads(rp.read_text())
                edit_report(js)
                rp.write_text(json.dumps(js))
            if drop:
                (f / drop).unlink()
            return b
        cases = [("copy", {}, "is not in this repository's runs folder"),
                 ("tilt", dict(edit_b=shifted("tilt", 0.7)), "its numbers are not the B of its report"),
                 ("sha", dict(edit_b=lambda t: t.replace(G.own_sha256(), "1" * 64)), "not the pre-registered one"),
                 ("noreport", dict(drop=G.OUT_REPORT.format(cam="c2") + ".json"), "is not beside it"),
                 ("trial", dict(edit_report=lambda js: js["run"].update(eligible=False)), "a TRIAL run"),
                 ("modules", dict(edit_report=lambda js: js["run"]["modules_sha256"].update({"georectify.py": "0" * 64})),
                  "other modules than the pre-registered ones (georectify.py)"),
                 ("norot", dict(edit_b=lambda t: t.replace(" rotated by ", " moved by ")), "rotated by <angle> deg"),
                 ("withheld", dict(edit_report=lambda js: js.update(delivery="WITHHELD (both-or-none)")),
                  "it was not delivered")]
        for name, kw, words in cases:
            b = forged(name, **kw)
            pcx = sp.procedure_check(b)
            o = sp.calibration_origin(b, [str(HERE / "calibration")])
            check(not pcx["pinned"] and any(words in w for w in pcx["why"])
                  and sp.held_check(o, sp.blind_of("2024-10-23"))[0] == "PARTLY-CIRCULAR",
                  f"B {name}: '{words}': NOT the pre-registered procedure, never INDEPENDENT")
        # A first (its DEM time goes into B's README; B's search envelope is checked against A's detection)
        aa = full_args(output_root=str(d / "out"))
        pa, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, aa)
        fake_built_product(pa, aa)
        # a camera the procedure did not deliver: left out, said; a delivered one never
        a = full_args(output_root=str(d / "out_gcpfree"), eo={"c2": str(r1["c2"])}, calibration_variant="gcpfree")
        pb, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, a)
        check(list(pb["cams"]) == ["c2"] and pb["not_delivered"]["c1"]["delivery"] == "NOT DELIVERED"
              and pb["not_delivered"]["c1"]["report"].endswith("CACO03_c1_gcpfree-oct2024_report.json")
              and sp.option_refusal(pb, a) == "",
              "--eo c2 only, from a run that did NOT DELIVER c1: c1 left out of the plan (its report kept), not refused")
        check(sp.dem_stem(pb).name == "2024-10-23_gcpfree" and "--epsg 6348" in sp.cmd_text(sp.dem_cmds(pb, a)[2])
              and sp.map_png(pb, "plan").name == "2024-10-23_gcpfree_waterlines_plan.png"
              and sp.map_png(pb, "c2").name == "2024-10-23_gcpfree_waterlines_c2.png",
              "its DEM and map files: 2024-10-23_gcpfree_dem.tif, 2024-10-23_gcpfree_waterlines_plan.png ... (the "
              "variant in the name), tagged EPSG:6348")
        ident = sp.calibration_identity(pb, a)
        check(ident.startswith("calibration variant gcpfree (--eo): c2 CACO03_c2_gcpfree-oct2024_EO.yaml sha256 "
                               + sp.file_sha(r1["c2"])[:12] + " INDEPENDENT (procedure made with stated exposure to "
                                 "A: README EXPOSURE)")
              and "c1 NOT DELIVERED by the GCP-free procedure" in ident and "BLIND: survey held by Dr Chris Sherwood" in ident
              and ident.encode("ascii"),
              f"the GeoTIFF description: '{ident[:110]}...'")
        page = sp.dem_cmds(pb, a)[3]
        check("CACO03 Marconi Beach (calibration gcpfree: 2025-01-23 GCPs, turned)" in sp.cmd_text(page),
              "its DEM page names the variant first ('calibration gcpfree: 2025-01-23 GCPs, turned'; dem_figure.py "
              "shrinks a heading that would run off the page)")
        a2 = full_args(output_root=str(d / "out2_gcpfree-aug"), eo={"c2": str(r2["c2"])},
                       calibration_variant="gcpfree-aug")
        p2, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, a2)
        msg = sp.option_refusal(p2, a2)
        check("c1" in p2["cams"] and "the procedure DELIVERED c1" in msg,
              "--eo c2 only, from a run that DELIVERED c1: REFUSED (a delivered camera is never left out)")
        a3 = full_args(output_root=str(d / "out3_gcpfree"), eo={"c1": str(r2["c1"]), "c2": str(r1["c2"])},
                       calibration_variant="gcpfree")
        p3, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, a3)
        check("different runs" in sp.option_refusal(p3, a3), "B files from two different runs: REFUSED")
        a4 = full_args(output_root=str(d / "out_gcpfree-aug"), eo={"c2": str(r1["c2"])}, calibration_variant="gcpfree-aug")
        p4, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, a4)
        check("as --calibration-variant gcpfree" in sp.option_refusal(p4, a4),
              "the first run's B under the second product's name (gcpfree-aug): REFUSED")
        # A, then B: README, provenance, the GeoTIFF description, the maps
        ida = sp.calibration_identity(pa, aa)
        check(ida.startswith("calibration of the table (survey_dates.csv): c1 CACO03_c1_20241023_EO.yaml sha256 ")
              and "PARTLY-CIRCULAR" in ida and "stated exposure" not in ida and sp.dem_stem(pa).name == "2024-10-23"
              and sp.map_png(pa, "plan").name == "2024-10-23_waterlines_plan.png",
              "A's GeoTIFF description: the table's calibration, its EO files and sha256, PARTLY-CIRCULAR; stem "
              "2024-10-23, maps 2024-10-23_waterlines_*")
        check("CACO03 Marconi Beach (table calibration 2024-10-23: solved from the held GCPs)"
              in sp.cmd_text(sp.dem_cmds(pa, aa)[3]),
              "A's DEM page heading says its calibration is the table's, solved from the held GCPs (PARTLY-CIRCULAR)")
        provb, rdb = fake_built_product(pb, a)
        unc = next((ln for ln in rdb.splitlines() if "STATED UNCERTAINTY" in ln), "")
        check(provb["status"] == "complete" and provb["blind"]["check"]["label"] == "INDEPENDENT"
              and "NOT DELIVERED: c1: the pre-registered GCP-free procedure did not deliver it" in rdb
              and "rotation +/- 4.0 px" in unc and "in the largest of five synthetic position-inheritance cases" in unc
              and "NEITHER term is a validated bound on real photos" in unc and "<=" not in unc
              and "real-texture" not in rdb and "|dz| median of 3-7 mm" not in rdb
              and "its DEM effect as a tilt: not computed" in rdb
              and "dem/          2024-10-23_gcpfree_dem/_spread/_count" in rdb
              and "maps/         2024-10-23_gcpfree_waterlines_<cam>.png" in rdb
              and provb["cameras_not_delivered"]["c1"]["verdict"] == "NOT DELIVERED"
              and any(c.startswith("c1 is NOT in this product") for c in provb["caveats"]),
              "B's README: c1 NOT DELIVERED with the procedure's report, B's two stated uncertainty terms (the "
              "position term the largest of five cases, neither a validated bound; no synthetic DEM figure for the "
              "rotation), the DEM and map files with the variant")
        dfiles = (provb.get("dem") or {}).get("files") or {}
        check(sorted(dfiles) == ["count", "dem", "spread"]
              and all(v == f"{pb['out'] / 'dem' / '2024-10-23_gcpfree'}_{k}" for k, v in dfiles.items()),
              "B's provenance.json names B's own DEM grids (dem/2024-10-23_gcpfree_dem, _spread, _count), never A's "
              f"names ({dfiles.get('dem')})")
        check("dem/          2024-10-23_gcpfree_dem/_spread/_count/_sigma .asc and .tif" in rdb,
              "B's README FILES names the _sigma grid too (dem_from_contours.sigma_grid), which dem/ and --share "
              "carry")
        gf = provb["blind"]["gcpfree"]
        check("EXPOSURE (stated by the people who made this procedure" in rdb and "A minus the Jan 2025 solution" in rdb
              and "ORDER (what decides B, and A): procedure pinned" in rdb and "outcome rule 'per-camera' pre-registered"
              in rdb and "A (the table's calibration) DEM first built 20" in rdb and "B run 20261014T010000Z" in rdb
              and gf["a_product"]["found"] and gf["times"]["prereg"] and gf["times"]["b_commit"],
              "B's README, BLIND: the procedure's EXPOSURE statement, and the ORDER of the pin, the pre-registration, "
              "A's DEM, the B run and its commit, with dates")
        check("NOTE: A's DEM was built BEFORE" not in rdb and "NOTE: B was committed after A's DEM" not in rdb
              and gf["times"]["a_dem"] >= gf["times"]["b_commit"] >= gf["times"]["prereg"],
              "... and what came first (here: the rule, B's commit, then A's DEM: no NOTE)")
        pr_c = pc["preregistration"]["commit"]
        rdb1 = " ".join(rdb.split())
        check("WHAT THE CODE ENFORCES (survey_products.py and gcpfree_eo.py, on this computer)" in rdb1
              and "WHAT RESTS ON THE OWNER FOLLOWING THE PROCEDURE (no software on a computer the owner controls can "
                  "prove it): that PREREGISTERED.txt was committed AND pushed to GitHub before any run of "
                  "gcpfree_eo.py, a trial included" in rdb1
              and f"PRE-REGISTRATION COMMIT: {pr_c}" in rdb1 and "the time GitHub received the push" in rdb1
              and "on NO remote branch this computer has fetched (not pushed, or not fetched since)" in rdb1
              and len(gf["protocol"]) == 3,
              "B's README, BLIND: what the code enforces, what rests on the owner (committed AND pushed before any "
              "run), and the pre-registration commit with where its public time is (here no remote holds it: said)")
        rpush = GcpfreeRepo(d / "pushed", push=True)
        rs_ = sp.commit_remote_state(rpush.prereg, rpush.info["commit"])
        txt_ = sp.prereg_commit_text({"preregistration": rpush.info, "preregistration_here": str(rpush.prereg)})
        check(rpush.info.get("remote_branches") and rs_["branches"] == rpush.info["remote_branches"]
              and rs_["url"].endswith("pushed_github.git") and f"held by {rs_['branches'][0]} of " in txt_
              and "NO remote branch" not in txt_ and "pushed (its commit on a remote-tracking branch" in sp.B_PROTOCOL_ENFORCED,
              f"a pushed pre-registration: gcpfree_eo.py allows a run (its commit on {rs_['branches']}), and B's README "
              f"says 'held by {rs_['branches'][0]} of <the remote>'; the code's part says it refuses until pushed")
        info_ = {"procedure_sha256": pc["procedure_sha256"], "preregistration": pc["preregistration"],
                 "outcome_rule": "per-camera", "run_id": "x", "run_started_utc": "2026-10-14 00:00:00",
                 "committed": dict(pc["committed"], date="2099-01-01T00:00:00+00:00"),
                 "preregistration_here": str(d / "no_git" / "PREREGISTERED.txt")}
        g2 = sp.gcpfree_disclosure(pb, a, {"c2": dict(info_, preregistration=dict(
            pc["preregistration"], commit_date="2099-01-01T00:00:00+00:00"))})
        check(any("A's DEM was built BEFORE the outcome rule in force was pre-registered" in x
                  for x in g2["order_flags"])
              and any("B was committed after A's DEM was built" in x for x in g2["order_flags"]),
              "... A's DEM before the rule's pre-registration, or before B's commit: said (NOTE lines)")
        # the rule changed (before any run) after A's DEM: the rule's history is said, from when the rule in force
        # was first committed
        rr = GcpfreeRepo(d / "rule_history", rule="both-or-none")
        rr.prereg.write_text(rr.prereg.read_text().replace("outcome_if_one_camera: both-or-none",
                                                           "outcome_if_one_camera: per-camera"))
        git_in(rr.root, "add", "-A")
        subprocess.run(["git", "-C", str(rr.root), "-c", "user.name=test", "-c", "user.email=test@example.org", "-c",
                        "commit.gpgsign=false", "commit", "-q", "-m", "pre-registration: the other rule"],
                       env=dict(os.environ, GIT_COMMITTER_DATE="2099-01-02T00:00:00+00:00"),
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        hist = sp.prereg_rule_history(rr.prereg)
        g3 = sp.gcpfree_disclosure(pb, a, {"c2": dict(info_, preregistration_here=str(rr.prereg))})
        check([h["rule"] for h in hist] == ["both-or-none", "per-camera"]
              and any("the outcome rule of PREREGISTERED.txt was changed before the run ('both-or-none'" in x
                      for x in g3["order_flags"])
              and any("A's DEM was built BEFORE the outcome rule in force was pre-registered" in x
                      for x in g3["order_flags"])
              and f"first committed with this rule {hist[1]['commit'][:10]}" in g3["order"][1]
              and "the file's commits: " in g3["order"][1],
              "... the outcome rule changed between the pre-registration's commits (after A's DEM): its history and "
              "NOTE lines; the rule in force dated from its own first commit")
        check("plan view (BLIND: no survey drawn; coloured by still water + setup" in rdb
              and "its title names the calibration" in rdb,
              "B's README FILES: the blind plan view's line reads as a sentence; its title names the calibration")
        mp = sp.read_stamp(pb["out"], "map_plan")
        check(mp and mp["sig"].get("calibration", "").startswith("calibration variant gcpfree (--eo): c2 "
                                                                   "CACO03_c2_gcpfree-oct2024_EO.yaml")
              and sp.map_png(pb, "plan").exists(),
              "B's plan view: drawn as 2024-10-23_gcpfree_waterlines_plan.png, titled with the calibration identity")
        check(sp.per_day_rows_note(pb).startswith("This product has no comparison (blind)")
              and sp.per_day_rows_note(pa) == sp.PER_DAY_ROWS_NOTE_BLIND
              and "compare/" not in sp.PER_DAY_ROWS_NOTE_BLIND and "show whether the lines moved" not in rdb
              and "Nothing was shifted: a survey is to be compared at the DEM's coordinates" in rdb
              and "The surveys are sampled at the DEM's coordinates" not in rdb,
              "a blind README never points at a comparison it does not have (per-day rows, the surveys sampled)")
        asc = d / "geo" / "t_dem.asc"
        asc.parent.mkdir()
        asc.write_text("ncols 2\nnrows 2\nxllcorner 420000\nyllcorner 4638000\ncellsize 1\nNODATA_value -9999\n"
                       "1 2\n3 -9999\n")
        rr = subprocess.run([sys.executable, str(HERE / "asc_to_geotiff.py"), str(asc), "--note", ident],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        desc, epsg = tiff_description_and_epsg(asc.with_suffix(".tif"))
        check(rr.returncode == 0 and epsg == 6348 and desc.endswith(ident) and "EPSG:6348+5703" in desc,
              "asc_to_geotiff.py (its default now EPSG:6348) --note: the GeoTIFF reads back EPSG 6348 and the "
              "calibration identity")
        # B is shared only from the root its first build recorded in its run folder (B_BUILD.txt), committed once
        rc, out = quiet(sp.share, share_args("out_gcpfree", "share_nob"))
        check(rc == sp.EXIT_CODES["refused"] and f"no {sp.B_BUILD_FILE} in its run folder" in out
              and not (d / "share_nob").exists(), "--share of a B whose first build recorded no root: REFUSED")
        pcs_b = sp.gcpfree_b_of(pb, a)
        recf = sp.write_b_build_record(pb, a, pcs_b)
        rec = sp.b_build_record(pcs_b)
        check(recf == runs / "20261014T010000Z" / sp.B_BUILD_FILE and rec["root"] == os.path.realpath(str(d / "out_gcpfree"))
              and "--calibration-variant gcpfree" in rec["command"] and sp.write_b_build_record(pb, a, pcs_b) is None,
              "a B's first build records its root and command in its run folder (B_BUILD.txt), once")
        rc, out = quiet(sp.share, share_args("out_gcpfree", "share_nob"))
        check(rc == sp.EXIT_CODES["refused"] and "commit it (once) before the B is shared" in out,
              "... not yet committed: --share REFUSED")
        repo.commit("gcpfree: B build")
        sp._FILE_SHA.clear()
        a_other = full_args(output_root=str(d / "out_other_gcpfree"), eo={"c2": str(r1["c2"])},
                            calibration_variant="gcpfree")
        pbo, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, a_other)
        check("is built in " + rec["root"] + " only" in sp.blind_b_build_problem(pbo, a_other),
              "a second root for the same B: REFUSED (one B root: no second B can be built and compared)")
        for rt in ("out", "out_gcpfree"):
            rc, out = quiet(sp.share, share_args(rt, "share"))
            check(rc == 0, f"--share of {rt}: written" + ("" if rc == 0 else f" -- {out}"))
        ta, tb = d / "share" / "2024-10-23_A_official_calibration.tgz", d / "share" / "2024-10-23_B_gcpfree.tgz"
        with tarfile.open(str(ta)) as t:
            na = t.getnames()
            man_a = t.extractfile("2024-10-23_A_official_calibration/MANIFEST.txt").read().decode()
        with tarfile.open(str(tb)) as t:
            nb = t.getnames()
            man_b = t.extractfile("2024-10-23_B_gcpfree/MANIFEST.txt").read().decode()
        check(ta.exists() and tb.exists() and (d / "share" / "2024-10-23_B_gcpfree.tgz.sha256").exists()
              and "2024-10-23_A_official_calibration/dem/2024-10-23_dem.tif" in na
              and "2024-10-23_B_gcpfree/dem/2024-10-23_gcpfree_dem.tif" in nb
              and "2024-10-23_A_official_calibration/maps/2024-10-23_waterlines_plan.png" in na
              and "2024-10-23_B_gcpfree/maps/2024-10-23_gcpfree_waterlines_plan.png" in nb
              and "2024-10-23_B_gcpfree/maps/2024-10-23_gcpfree_waterlines_c2.png" in nb
              and "2024-10-23_A_official_calibration/calibration/CACO03_c1_20241023_EO.yaml" in na
              and "2024-10-23_B_gcpfree/calibration/gcpfree_run_20261014T010000Z/CACO03_c1_gcpfree-oct2024_report.json" in nb
              and "2024-10-23_B_gcpfree/calibration/gcpfree_run_20261014T010000Z/run.json" in nb
              and "2024-10-23_B_gcpfree/calibration/gcpfree_run_20261014T010000Z/B_BUILD.txt" in nb
              and "2024-10-23_B_gcpfree/calibration/CACO03_c2_20240801_IO.yaml" in nb
              and "2024-10-23_A_official_calibration/calibration/CACO03_c1_20240801_IO.yaml" in na
              and not any("2024-10-23_A" in n for n in nb) and not any("2024-10-23_B" in n for n in na),
              "two archives with distinct names: A (the table's calibration) and B (gcpfree, with the procedure's reports, "
              "c1's NOT DELIVERED report among them); their maps and DEM files named apart too")
        check("2024-10-23_B_gcpfree/calibration/PREREGISTERED.txt" in nb and not any("PREREGISTERED" in n for n in na)
              and f"pre-registration: PREREGISTERED.txt sha256 {prereg_sha[:16]}, the file the procedure's run recorded"
              in man_b, "B's archive holds the pre-registration (the outcome rule fixed before the run), its MANIFEST "
                        "saying it is the file the run recorded; A's holds none")
        check("his check of this DEM: PARTLY-CIRCULAR" in man_a and "his check of this DEM: INDEPENDENT" in man_b
              and "c1: NOT DELIVERED by the pre-registered GCP-free procedure" in man_b and "BLIND" in man_a
              and "c2 procedure: made by the pre-registered procedure" in man_b
              and str(runs / "20261014T010000Z") in man_b and "EXPOSURE (stated by" in man_b
              and "ORDER (what decides B, and A)" in man_b and "STATED UNCERTAINTY" in man_b
              and "WHAT THE CODE ENFORCES" in man_b and "WHAT RESTS ON THE OWNER" in man_b
              and "PRE-REGISTRATION COMMIT: " in man_b and "WHAT RESTS ON THE OWNER" not in man_a
              and "EXPOSURE" not in man_a
              and "lens c2: CACO05_c2_20240801_IO.yaml is the file the detector reads" in man_b
              and "CACO03_c2_20240801_IO.yaml, this CACO03 product's own lens file, is in the archive too: the same "
                  "lens" in man_b,
              "their MANIFESTs: what each is, its calibration and the holder's-check label, every file's sha256; B's "
              "also its procedure (run folder, pre-registration, commit), EXPOSURE, ORDER and stated uncertainty")
        ref_b = "CACO03_c2_20250123_EO.yaml"
        check("2024-10-23_B_gcpfree/waterlines/contour_points_ground.csv" in nb
              and "2024-10-23_A_official_calibration/waterlines/contour_points_ground.csv" in na
              and "2024-10-23_B_gcpfree/pointing/pointing_c2.json" in nb
              and f"2024-10-23_B_gcpfree/calibration/{ref_b}" in nb and f"calibration/{ref_b}" not in " ".join(na)
              and f"reference EO {ref_b} (what B was turned from, b_reference_eo): sha256 " in man_b
              and "the file B's report recorded" in man_b
              and "What of the build's folder is here (README: FILES)" in man_b
              and "waterlines/contour_points_ground.csv (the waterline points the DEM used" in man_b
              and not any("/waterlines/c2/" in n or "reference_c2" in n or "/stamps/" in n for n in nb),
              "the archives hold what their README points to: the waterlines the DEM used (a check by day or by "
              "frame), the pointing tables, the forcing and logs; B's also the GCP solution it was turned from, with "
              "the sha256 its report recorded; never the detection's work folders or the reference photo copy")
        check("his check of this DEM: INDEPENDENT (procedure made with stated exposure to A: README EXPOSURE" in man_b
              and "stated exposure" not in man_a,
              "B's MANIFEST label line points to the EXPOSURE statement (A's has none)")
        rc, out = quiet(sp.share, share_args("out_gcpfree", "share"))
        check(rc == sp.EXIT_CODES["refused"] and "never written over" in out, "--share again: never written over")
        # the pre-registration changed after the run (not yet committed): B is no longer the pre-registered one
        prereg_txt = repo.prereg.read_text()
        repo.prereg.write_text(prereg_txt.replace("outcome_if_one_camera: per-camera", "outcome_if_one_camera: both-or-none"))
        rc, out = quiet(sp.share, share_args("out_gcpfree", "share_p"))
        check(rc == sp.EXIT_CODES["refused"] and "not the file its run recorded" in out
              and not (d / "share_p" / "2024-10-23_B_gcpfree.tgz").exists(),
              "the pre-registration changed after the run: --share REFUSED (B is not the pre-registered procedure's)")
        repo.prereg.write_text(prereg_txt)
        # the calibration file there now is not the one the product was built with: never packed with its DEM
        provp_ = json.loads((d / "out" / "2024-10-23" / "provenance.json").read_text())
        provp_["cameras"]["c1"]["calibration"]["eo_sha256"] = "0" * 16
        (d / "out" / "2024-10-23" / "provenance.json").write_text(json.dumps(provp_))
        rc, out = quiet(sp.share, share_args("out", "share_a2"))
        check(rc == sp.EXIT_CODES["refused"] and "is not the one this product was built with" in out
              and not (d / "share_a2").exists(),
              "A: a calibration file that is not the one the product was built with (its sha256): --share REFUSED")
        reseal(runs / "20261014T010000Z", "c2", edit_b=shifted("azimuth", 0.4))      # B edited, report resealed
        rc, out = quiet(sp.share, share_args("out_gcpfree", "share2"))
        check(rc == sp.EXIT_CODES["refused"] and "NOT the pre-registered procedure's B" in out
              and "is not the one this product was built with" in out and "changed since its commit" in out
              and not (d / "share2" / "2024-10-23_B_gcpfree.tgz").exists(),
              "--share of a B edited after it was built (report and run.json rewritten to agree): REFUSED, no archive")
    # the runtime guard sees through a link; a step is never handed such a file
    real = d / "held" / "2024-10-23_Marconi_GCP_xyz.csv"
    real.parent.mkdir()
    real.write_text("1,2,3,4\n")
    os.symlink(str(real), str(d / "held" / "data.csv"))
    check(sp.held_file_name(d / "held" / "data.csv") is not None and sp.held_file_name(d / "held" / "other.csv") is None,
          "the runtime guard knows a file by the name it really has (a link 'data.csv' to a held survey's file)")
    if sp.arm_blind_guard():
        try:
            open(str(d / "held" / "data.csv")).read()
            check(False, "a link to a held survey: NOT refused")
        except sp.BlindSurvey:
            check(True, "... and refuses to open it through the link (a Python-level open in this process)")
        # a photo-named link to it in a photo folder: the photo scan of this process (it runs before any step that
        # is handed the photo folders) probes its size with PIL, a Python-level open: refused, nothing read
        from detect_original_view import collect_photos
        pics = d / "held" / "pics"
        pics.mkdir()
        os.symlink(str(real), str(pics / "1729699200.Wed.Oct.23_16_00_00.GMT.2024.CACO03.c1.timex.jpg"))
        try:
            collect_photos([str(pics)], "c1", "2024-10-20", "2024-10-26", (0.0, 24.0), station="CACO03",
                           frame_size=(2448, 2048))
            check(False, "a photo-named link to a held survey: NOT refused by the photo scan")
        except sp.BlindSurvey:
            check(True, "a photo-named link to a held survey in a photo folder: the in-process photo scan is refused "
                        "(BlindSurvey, before any step is started on the folder)")
    try:
        quiet(sp.run_cmd, [sys.executable, str(HERE / "asc_to_geotiff.py"), str(d / "held" / "data.csv")],
              d / "held" / "x.log")
        check(False, "a step handed a link to a held survey: NOT refused")
    except sp.BlindSurvey as exc:
        check("refused to start asc_to_geotiff.py" in str(exc) and not (d / "held" / "x.log").exists(),
              "a step (a child process, beyond the hook) is never handed such a file: refused before it starts")
    check(held is sp.BLIND_DATES["2024-10-23"], "BLIND_DATES untouched")


def test_gcpfree_b_build(d):
    print("a GCP-free B of a blind date is the pre-registered B BUILD: the committed tables, A's options, envelope and "
          "forcing, one root; A's first DEM time kept; every run listed; an uncommitted trial stops it")
    d = Path(d) / "gcpfree_b_build"
    d.mkdir()
    repo = GcpfreeRepo(d / "repo")
    # this repository's committed tables (the station's) stand in for surveys/*.csv
    tdir = repo.root / "surveys"
    tdir.mkdir()
    (tdir / "survey_dates.csv").write_text(FIXTURE_DATES)
    (tdir / "surveys.csv").write_text(FIXTURE_SURVEYS)
    repo.commit("the tables")
    cfg, srv = sp.read_table(tdir / "survey_dates.csv", sp.CONFIG_FIELDS), sp.read_table(tdir / "surveys.csv",
                                                                                          sp.SURVEY_FIELDS)
    sdir = d / "surveys_dir"
    sdir.mkdir()
    env = sdir / "2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif"
    env.write_bytes(b"II*\0a stand-in for the Jan 2025 lidar")
    forcing = {}
    for k in ("adcp", "adcp_navd88", "chatham", "wis", "ndbc"):
        forcing[k] = str(d / f"{k}.csv")
        Path(forcing[k]).write_text("")

    def args_(**kw):
        base = dict(survey_dirs=[str(sdir)], calibration=str(HERE / "calibration"), force_disabled=True,
                    photo_roots=[str(d / "no_photos")], live_contours="/x/contour_points_timex.csv",
                    setup_fitted_to=None, geotiff_epsg=sp.GEOTIFF_EPSG, gnssr_spline=str(d / "no_spline.txt"),
                    gauge_csv=str(d / "no_gauge.csv"), live_cron=str(sp.LIVE_CRON), pointing_log_dir=str(d / "no_log"),
                    pointing_frames=3, pointing_fail=sp.POINT_FAIL_DEG, dem_cell=sp.DEM_CELL,
                    dem_min_points=sp.DEM_MIN_POINTS, dem_max_spread=sp.DEM_MAX_SPREAD, dem_max_hs=sp.DEM_MAX_HS,
                    dem_max_day_offset=sp.DEM_MAX_DAY_OFFSET, keep_moved_days=False, no_download=True,
                    synthetic=None, replace=False, dry_run=False, force=False, eo=None, calibration_variant=None,
                    config=str(tdir / "survey_dates.csv"), surveys=str(tdir / "surveys.csv"), a_root=None)
        base.update(forcing)
        base.update(kw)
        return cfg_args(**base)
    real_cfg, real_srv = sp.DEFAULT_CONFIG, sp.DEFAULT_SURVEYS
    sp.DEFAULT_CONFIG, sp.DEFAULT_SURVEYS = tdir / "survey_dates.csv", tdir / "surveys.csv"
    try:
        with repo:
            r1 = repo.run("20261014T010000Z", {"c1": gcpfree_res("c1", (0, 0, 0), "NOT DELIVERED"),
                                               "c2": gcpfree_res("c2", (-13.4, -0.5, 0.3))}, "2026-10-14 01:00:00")
            # A first: its DEM first built at 12:00, rebuilt at 14:00 (a later merge reruns its DEM step)
            aa = args_(output_root=str(d / "out"))
            pa, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, aa)
            fake_built_product(pa, aa)
            for k in ("water_level.csv", "waves.csv"):
                (pa["out"] / "forcing").mkdir(exist_ok=True)
                (pa["out"] / "forcing" / k).write_text("time,level\n2024-10-23 15:00,0.5\n")
            st = sp.read_stamp(pa["out"], "dem")
            st["finished_utc"] = st["first_finished_utc"] = "2026-10-01T12:00:00.000000Z"
            sp.stamp_path(pa["out"], "dem").write_text(json.dumps(st))
            (pa["out"] / "logs" / sp.FIRST_BUILT).unlink()        # a product of the code before the record existed
            sp.drop_stamp(pa["out"], "dem")                         # its DEM step reruns...
            sp.write_stamp(pa["out"], "dem", st["sig"], extra={"settings": st.get("settings")})   # ... now
            check(sp.first_built(pa["out"], "dem") == "2026-10-01T12:00:00.000000Z"
                  and sp.read_stamp(pa["out"], "dem")["first_finished_utc"] == "2026-10-01T12:00:00.000000Z"
                  and sp.read_stamp(pa["out"], "dem")["finished_utc"] > "2026-10-01T12:00:00",
                  "A's DEM rebuilt: the first finish (1 Oct 12:00) is kept apart (logs/first_built.json, the stamp's "
                  "first_finished_utc), the stamp's finished_utc is the rebuild's")
            ab = args_(output_root=str(d / "out_gcpfree"), eo={"c2": str(r1["c2"])}, calibration_variant="gcpfree")
            pb, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, ab)
            pcs = sp.gcpfree_b_of(pb, ab)
            check(pcs and sp.blind_b_build_problem(pb, ab) == "",
                  "B built as A was (the committed tables, A's options, the envelope A's detection used, A's "
                  "forcing inputs): the pre-registered B build" + ("" if not pcs else
                                                                   f" -- {sp.blind_b_build_problem(pb, ab)}"))
            pc = pcs["c2"]
            info_ = dict(pc, preregistration=dict(pc["preregistration"], commit_date="2026-10-01T13:00:00+00:00"),
                         preregistration_here=str(d / "no_git" / "PREREGISTERED.txt"),
                         run_started_utc="2026-10-01 13:10:00", run_finished_utc="2026-10-01 13:20:00",
                         committed=dict(pc["committed"], date="2026-10-01T13:25:00+00:00"))
            g = sp.gcpfree_disclosure(pb, ab, {"c2": info_})
            a_line = next(x for x in g["order"] if x.startswith("A "))
            check("DEM first built 2026-10-01 12:00 UTC (last rebuilt " in a_line
                  and any("A's DEM was built BEFORE the outcome rule in force was pre-registered" in x
                          for x in g["order_flags"])
                  and any("B was committed after A's DEM was built" in x for x in g["order_flags"]),
                  f"B's ORDER after A's DEM step was rebuilt: '{a_line[:70]}...', and the NOTEs stay (A came first)")
            check(any("A's DEM was rebuilt (" in x and "after B's run began" in x for x in g["order_flags"]),
                  "A's DEM rebuilt after B's run began (B is built with A's options as last built): a NOTE says so")
            g_first = sp.gcpfree_disclosure(pb, ab, {"c2": dict(info_, run_started_utc="2026-10-01 11:00:00",
                                                                run_finished_utc="2026-10-01 11:10:00")})
            fl = [x for x in g_first["order_flags"] if "after B's run began" in x]
            check(len(fl) == 1 and fl[0].startswith("A's DEM was first built (2026-10-01 12:00 UTC) after B's run began")
                  and "rebuilt (" not in fl[0] and "and of its last rebuild" in fl[0],
                  "B's run began before A's DEM was FIRST built: the NOTE says 'first built' (not 'rebuilt', which "
                  "contradicted the ORDER's 'first built'), and names the last rebuild too")
            g_early = sp.gcpfree_disclosure(pb, ab, {"c2": dict(info_, run_started_utc="2099-01-01 00:00:00",
                                                                run_finished_utc="2099-01-01 00:10:00")})
            check(not any("after B's run began" in x for x in g_early["order_flags"]),
                  "... and no such NOTE when A's last build came before B's run")
            g0 = sp.gcpfree_disclosure(pb, args_(output_root=str(d / "out_gcpfree"), eo={"c2": str(r1["c2"])},
                                                 calibration_variant="gcpfree", a_root=str(d / "nowhere")),
                                       {"c2": info_})
            check(any("the order of B relative to A is not known" in x for x in g0["order_flags"]),
                  "no A product found: a NOTE says the order relative to A is not known")
            check(any("every run of the procedure in" in x and "20261014T010000Z delivery (B's run)" in x
                      for x in g["order"]), "B's ORDER lists every run of the procedure there, B's marked")

            def problem(**kw):
                a_ = args_(**dict(dict(output_root=str(d / "out_gcpfree"), eo={"c2": str(r1["c2"])},
                                       calibration_variant="gcpfree"), **kw))
                p_, _ = quiet(sp.build_plan, "2024-10-23", cfg, srv, a_)
                return sp.blind_b_build_problem(p_, a_), p_, a_
            msg, _, _ = problem(dem_cell=1.0)
            check("its build options are not A's as built (--dem-cell 1.0 (A: (default)))" in msg,
                  "B with another DEM cell than A's: REFUSED (" + msg[:80] + ")")
            msg, _, _ = problem(keep_moved_days=True, window=["2024-10-21", "2024-10-25"])
            check("--keep-moved-days" in msg and "--window 2024-10-21 2024-10-25 (A: (default))" in msg,
                  "B with --keep-moved-days and another --window: REFUSED")
            other_tab = d / "tables_other" / "survey_dates.csv"
            other_tab.parent.mkdir()
            other_tab.write_text(FIXTURE_DATES)
            msg, _, _ = problem(config=str(other_tab))
            check("--config " + str(other_tab) + ": a B is built with the committed survey_dates.csv" in msg,
                  "B with another --config (e.g. its envelope A's DEM): REFUSED")
            msg, _, _ = problem(a_root=str(d / "nowhere"))
            check("no table-calibration (A) product of 2024-10-23" in msg, "B without A's product: REFUSED (--a-root)")
            # the envelope: changed since A's detection, or a DEM of this pipeline (A's DEM under the lidar's name)
            keep = env.read_bytes()
            env.write_bytes(keep + b" edited")
            sp._FILE_SHA.clear()
            msg, p_, a_ = problem()
            ec = sp.envelope_check(p_, "c2", a_)
            lab = sp.held_check(sp.calibration_origin(r1["c2"], [str(HERE / "calibration")]), p_["blind"], ec)[0]
            check("is not the file A's detection used (sha256" in msg and not ec["ok"] and lab == "PARTLY-CIRCULAR",
                  "the envelope's content changed since A's detection: REFUSED, and his check PARTLY-CIRCULAR")
            env.write_bytes(keep)
            sp._FILE_SHA.clear()
            asc = d / "geo" / "2024-10-23_dem.asc"
            asc.parent.mkdir()
            asc.write_text("ncols 2\nnrows 2\nxllcorner 420000\nyllcorner 4638000\ncellsize 1\nNODATA_value -9999\n"
                           "1 2\n3 -9999\n")
            subprocess.run([sys.executable, str(HERE / "asc_to_geotiff.py"), str(asc), "--note",
                            sp.calibration_identity(pa, aa)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            adem = asc.with_suffix(".tif")
            fake_env = {"envelope_path": adem}
            p2 = dict(p_, cams={"c2": dict(p_["cams"]["c2"], **fake_env)}, out=d / "out_none" / "2024-10-23")
            ec2 = sp.envelope_check(p2, "c2", a_)
            lab2 = sp.held_check(sp.calibration_origin(r1["c2"], [str(HERE / "calibration")]), p2["blind"], ec2)[0]
            check(not ec2["ok"] and "is a DEM of this pipeline" in ec2["text"] and "held survey's calibration (A)" in
                  ec2["text"] and lab2 == "PARTLY-CIRCULAR",
                  "an envelope that is A's DEM (built with calibration A): PARTLY-CIRCULAR, never INDEPENDENT")
            # the forcing: A's water level and waves, byte for byte
            (pb["out"] / "forcing").mkdir(parents=True)
            for k in ("water_level.csv", "waves.csv"):
                (pb["out"] / "forcing" / k).write_text("time,level\n2024-10-23 15:00,0.5\n")
            ok_f = sp.b_forcing_problem(pb, ab)
            (pb["out"] / "forcing" / "water_level.csv").write_text("time,level\n2024-10-23 15:00,0.6\n")
            bad_f = sp.b_forcing_problem(pb, ab)
            check(ok_f == "" and "its forcing (water_level.csv) is not A's" in bad_f,
                  "after B's forcing step: A's water level and waves byte for byte, else the build stops")
            # a trial (or a run cut off) left uncommitted: not the pre-registered procedure until it is committed
            repo.run("20261014T020000Z_trial", {"c2": gcpfree_res("c2", (-13.4, -0.5, 0.3))},
                          "2026-10-14 02:00:00", eligible=False, commit=False)
            sp._FILE_SHA.clear()
            pct = sp.procedure_check(r1["c2"])
            check(not pct["pinned"] and any("the record of run(s) 20261014T020000Z_trial" in w and "is not committed" in w
                                            for w in pct["why"]),
                  "an uncommitted TRIAL run in runs/: B is not (yet) the pre-registered procedure's (commit it)")
            repo.commit("gcpfree: trial 20261014T020000Z_trial")
            pct = sp.procedure_check(r1["c2"])
            g2 = sp.gcpfree_disclosure(pb, ab, {"c2": dict(pct, run_started_utc="2026-10-14 03:00:00")})
            check(pct["pinned"] and any("20261014T020000Z_trial TRIAL" in x for x in g2["order"])
                  and any("TRIAL run(s) 20261014T020000Z_trial (a full result" in x for x in g2["order_flags"]),
                  "... committed: B pinned again; its ORDER lists the trial, and a NOTE says a full result was seen "
                  "before B's run" + ("" if pct["pinned"] else f" -- {pct['why']}"))
            # a reset away from the first run, a new 'first' run: the old one is still in the reflog
            rr = GcpfreeRepo(d / "reset")
            with rr:
                ra = rr.run("20261014T010000Z", {"c2": gcpfree_res("c2", (-13.4, -0.5, 0.3))}, "2026-10-14 01:00:00")
                git_in(rr.root, "reset", "-q", "--hard", "HEAD~1")
                rb_ = rr.run("20261014T030000Z", {"c2": gcpfree_res("c2", (-13.4, -0.4, 0.3))}, "2026-10-14 03:00:00")
                pcr = sp.procedure_check(rb_["c2"])
                check(not ra["c2"].exists() and not pcr["pinned"] and pcr["delivery_runs_here"] == ["20261014T030000Z"]
                      and any("outside this branch's history" in w for w in pcr["why"]),
                      "the first run reset away (git reset --hard), a new run committed: the old run is still in the "
                      "reflog: NOT the pre-registered procedure")
    finally:
        sp.DEFAULT_CONFIG, sp.DEFAULT_SURVEYS = real_cfg, real_srv


def test_review2_nits(d):
    print("a restored calibration rebuilds nothing; a photo-named link to the held survey is refused (exit 5) before "
          "any step; the tilt table of a thin DEM strip, and why none; the same B file edited; the lens")
    from georectify import load_intrinsics, load_extrinsics, pixel_to_ground
    d = Path(d) / "review2"
    d.mkdir()
    # (1) the calibration file counted by content: touched (git checkout) -> fresh, its token kept; edited -> rebuilt
    out = d / "fresh" / "2024-10-23"
    eo = d / "fresh" / "b.yaml"
    out.mkdir(parents=True)
    shutil.copy(str(HERE / "calibration" / "CACO03_c2_20250123_EO.yaml"), str(eo))
    o = out / "pointing" / "pointing_c2.json"
    o.parent.mkdir(parents=True)
    o.write_text("{}")
    sp._FILE_SHA.clear()
    sp.write_stamp(out, "pointing_c2", {"eo_file": sp.file_sig(eo), "x": 1}, extra={"eo_sha256": sp.file_sha(eo)})
    tok = sp.token(out, "pointing_c2")
    t0 = os.stat(str(eo)).st_mtime
    os.utime(str(eo), (t0 + 100, t0 + 100))
    fresh, why = sp.is_fresh(out, "pointing_c2", {"eo_file": sp.file_sig(eo), "x": 1}, [o], [eo])
    check(fresh and sp.token(out, "pointing_c2") == tok
          and sp.read_stamp(out, "pointing_c2")["sig"]["eo_file"] == sp.normal(sp.file_sig(eo)),
          "a calibration file touched (restored by git, the same bytes): the step is fresh, its token kept (nothing "
          "downstream rebuilds), the stamp records the file's time now" + ("" if fresh else f" -- {why}"))
    eo.write_text(eo.read_text().replace("tilt:", "tilt:  ", 1)[:-2] + "\n\n")
    os.utime(str(eo), (t0 + 200, t0 + 200))
    sp._FILE_SHA.clear()
    fresh2, why2 = sp.is_fresh(out, "pointing_c2", {"eo_file": sp.file_sig(eo), "x": 1}, [o], [eo])
    check(not fresh2 and why2 == "its settings changed", f"... its content edited: rebuilt ({why2})")
    # (2) a photo-named link to a held survey's file in a photo folder: refused (exit 5) before any step, a dry run too
    cfg, srv, tabs = fixture_tables(d)
    held = d / "held" / "2024-10-23_Marconi_GCP_points.csv"
    held.parent.mkdir()
    held.write_text("1,2,3,4\n")
    pics = d / "pics"
    pics.mkdir()
    os.symlink(str(held), str(pics / "1729699200.Wed.Oct.23_16_00_00.GMT.2024.CACO03.c1.timex.jpg"))
    if sp.arm_blind_guard():
        for dry in (True, False):
            a_ = cfg_args(output_root=str(d / "trap_out"), photo_roots=[str(pics)], force_disabled=True, dry_run=dry,
                          eo=None, calibration_variant=None, survey_dirs=[str(d)], **tabs)
            for k, v in dict(live_contours="/x/c.csv", setup_fitted_to=None, geotiff_epsg=sp.GEOTIFF_EPSG,
                             gnssr_spline="/x/s.txt", gauge_csv="/x/g.csv", live_cron=str(sp.LIVE_CRON),
                             pointing_log_dir="/x", pointing_frames=3, pointing_fail=sp.POINT_FAIL_DEG,
                             dem_cell=sp.DEM_CELL, dem_min_points=sp.DEM_MIN_POINTS, dem_max_spread=sp.DEM_MAX_SPREAD,
                             dem_max_hs=sp.DEM_MAX_HS, dem_max_day_offset=sp.DEM_MAX_DAY_OFFSET, keep_moved_days=False,
                             no_download=True, synthetic=None, replace=False, force=False, adcp=None,
                             adcp_navd88=None, chatham=None, wis=None, ndbc=None, keep_detector_debug=False,
                             a_root=None).items():
                setattr(a_, k, v)
            sp._SCAN.clear()
            res, txt = quiet(sp.build_date, "2024-10-23", cfg, srv, a_, list(sp.STEPS))
            check(res == "refused" and "REFUSED (exit 5): a file in the photo folders is (or links to) a file named "
                                       "like the held survey" in txt and not (d / "trap_out").exists(),
                  f"a photo-named link to the held survey's file ({'dry run' if dry else 'build'}): REFUSED (exit 5) "
                  f"before any step, nothing written")
    # (3) the rotation term's tilt table on a thin DEM strip (two rows: no central slope), and why none on a few cells
    cal = HERE / "calibration"
    io_f, eo_f = cal / "CACO05_c2_20240801_IO.yaml", cal / "CACO03_c2_20250123_EO.yaml"
    io_, e_ = load_intrinsics(io_f), load_extrinsics(eo_f)
    U = np.linspace(300, 2100, 40)
    X, Y = pixel_to_ground(U, np.full_like(U, 1500.0), np.full_like(U, 0.5), io_, e_)
    cell = 2.0
    xll, yll = float(np.floor(X.min())) - 10, float(np.floor(Y.min())) - 10
    nc, nr = int((X.max() - xll) / cell) + 12, int((Y.max() - yll) / cell) + 12
    for name, rows_used, cols in (("strip", 2, slice(None)), ("tiny", 2, slice(0, 5))):
        g = np.full((nr, nc), -9999.0)
        r0 = nr - 1 - int((float(np.median(Y)) - yll) / cell)
        cc = np.arange(nc)[cols]
        for k in range(rows_used):
            g[r0 + k, cc] = 0.5 + 0.03 * cc + 0.1 * k
        plan_ = {"date": "2024-10-23", "era": "chatham", "out": d / name / "2024-10-23",
                 "cams": {"c2": {"eo_path": eo_f, "io_path": io_f, "variant": "gcpfree"}}}
        stem = sp.dem_stem(plan_)
        stem.parent.mkdir(parents=True)
        Path(f"{stem}_dem.asc").write_text(f"ncols {nc}\nnrows {nr}\nxllcorner {xll}\nyllcorner {yll}\ncellsize {cell}\n"
                                           f"NODATA_value -9999\n" + "\n".join(" ".join(f"{v:.3f}" for v in row)
                                                                             for row in g) + "\n")
        why_ = {}
        bands = sp.tilt_sensitivity(plan_, "c2", 0.06, 0.0, why=why_)
        if name == "strip":
            d_row, d_col = np.gradient(np.where(g > -9000, g, np.nan), cell)
            check(not np.isfinite(d_row).any() and bands and all(b.get("slope") for b in bands)
                  and "one-sided differences" in sp.tilt_effect_text({"tilt_deg": 0.06, "bands": bands}),
                  f"a DEM strip two cells wide (np.gradient finds no slope): the tilt table from one-sided slopes "
                  f"({len(bands or [])} band(s)), said")
        else:
            txt = sp.tilt_effect_text({"tilt_deg": 0.06, "bands": None, "why_none": why_.get("text")})
            check(bands is None and "10 with a local slope (even from one-sided differences), fewer than the 20 "
                                    "needed" in why_.get("text", "")
                  and "not computed (the DEM has 10 cells" in txt and "no DEM yet" not in txt,
                  f"a DEM of 10 cells: no tilt table, and the README says why ('{txt[:90]}...'), not 'no DEM yet'")
    # (4) the same B file edited after its build, in its own root: said as a changed file, not as another variant
    rootg = d / "sp_gcpfree"
    f = d / "runs" / "CACO03_c2_gcpfree-oct2024_EO.yaml"
    f.parent.mkdir()
    shutil.copy(str(eo_f), str(f))
    plan_ = {"date": "2024-10-23", "era": "chatham", "out": rootg / "2024-10-23",
             "cams": {"c2": {"eo_path": f, "eo_file": f.name, "variant": "gcpfree"}}}
    sp.write_stamp(plan_["out"], "detect_c2", {"eo_file": sp.file_sig(f)},
                   extra={"calibration": {"eo_file": f.name, "eo_path": str(f), "eo_sha256": "0" * 64,
                                          "variant": "gcpfree"}})
    sp._FILE_SHA.clear()
    note = sp.calibration_replace_note(plan_, SimpleNamespace(output_root=str(rootg), replace=False),
                                       ["pointing", "detect"])
    check(note and "the same --eo file, changed since this product was built" in note and "_gcpfree_gcpfree" not in note,
          "the same B file edited after its build: 'the same --eo file, changed since' (no new root advised)")


def test_gcpfree_where(d):
    print("a calibration B is the delivered one only from this repository's runs folder, under its pre-registration "
          "as committed, from the station's folders, the first delivery run, committed once, its lines true")
    import gcpfree_eo as G
    d = Path(d) / "gcpfree_where"
    c2 = (-13.4, -0.5, 0.3)

    def res_nd():
        return {"c1": gcpfree_res("c1", (0, 0, 0), "NOT DELIVERED"), "c2": gcpfree_res("c2", c2)}

    def why_of(b, *words):
        pc = sp.procedure_check(b)
        miss = [w for w in words if not any(w in x for x in pc["why"])]
        lab = sp.held_check(sp.calibration_origin(b, [str(HERE / "calibration")]), sp.blind_of("2024-10-23"))[0]
        return (not pc["pinned"] and not miss and lab == "PARTLY-CIRCULAR"), (pc["why"] if miss else "")
    # (1) a second pre-registration (the other rule), made after the official run's result was seen, in another
    # folder of the same repository, its run written as a delivery run would be: never the delivered B
    off = GcpfreeRepo(d / "official", rule="both-or-none")
    with off:
        r = off.run("20261014T010000Z", res_nd(), "2026-10-14 01:00:00")
        check(not r["c2"].exists(), "the official run, both-or-none, c1 not delivered: no B at all (c2 WITHHELD)")
        sh = G.write_preregistration(off.root / "shadow" / "PREREGISTERED.txt", "per-camera", "someone else", G.ROOTS,
                                     photos=FAKE_PHOTOS)
        off.commit("a shadow pre-registration")
        sh_info, _ = G.preregistration(sh)
        rs = off.run("20261014T030000Z", res_nd(), "2026-10-14 03:00:00", folder=sh.parent / "runs" / "20261014T030000Z",
                     prereg=sh_info)
        ok, w = why_of(rs["c2"], "is not in this repository's runs folder", "not this repository's")
        check(rs["c2"].exists() and ok, f"(1) a shadow pre-registration's 'delivery run' (per-camera, c2 B): NOT the "
                                         f"pre-registered procedure, never INDEPENDENT {w}")
        nd, bad = sp.not_delivered_cameras(["c1", "c2"], {"c2": str(rs["c2"])})
        check(not nd and "c1" in bad, "... and its c2-only B may not leave c1 out")
    # (2) a B of another repository (another clone, any rule): never this repository's
    other = GcpfreeRepo(d / "other")
    ro = other.run("20261014T010000Z", res_nd(), "2026-10-14 01:00:00")
    with off:
        ok, w = why_of(ro["c2"], "is not in this repository's runs folder", "not this repository's")
        check(ok, f"(2) a B made in another repository: NOT the pre-registered procedure {w}")
    # (3) edited after its commit, report and run.json rewritten: 'changed'; committed again: two commits
    rp = GcpfreeRepo(d / "recommit")
    with rp:
        rr = rp.run("20261014T010000Z", res_nd(), "2026-10-14 01:00:00")
        check(sp.procedure_check(rr["c2"])["pinned"], "(3) the B as the run wrote and committed it: pinned")
        reseal(rr["c2"].parent, "c2", edit_b=shifted("azimuth", 0.5))
        ok, w = why_of(rr["c2"], "its B file CACO03_c2_gcpfree-oct2024_EO.yaml is changed since its commit",
                       "its report CACO03_c2_gcpfree-oct2024_report.json is changed since its commit")
        check(ok, f"... azimuth +0.5 deg, report and run.json rewritten to agree, not committed: NOT the "
                  f"pre-registered procedure {w}")
        rp.commit("B nudged")
        ok, w = why_of(rr["c2"], "was touched by 2 commits")
        check(ok, f"... and committed again: two commits of the B file: NOT the pre-registered procedure {w}")
        (rr["c2"].parent / "run_report.txt").unlink()
        rp.commit("run_report.txt removed")
        ok, w = why_of(rr["c2"], "its run_report.txt was touched by 2 commits")
        check(ok, f"... its run_report.txt deleted in a later commit: NOT the pre-registered procedure {w}")
    # (4) not committed at all; then in another git repository nested inside
    nc = GcpfreeRepo(d / "nocommit")
    with nc:
        rn = nc.run("20261014T010000Z", res_nd(), "2026-10-14 01:00:00", commit=False)
        ok, w = why_of(rn["c2"], "is NOT committed to git: commit the run folder, in one commit")
        check(ok, f"(4) a run folder not committed yet: NOT (yet) the pre-registered procedure {w}")
        git_in(rn["c2"].parent, "init", "-q")
        git_in(rn["c2"].parent, "add", "-A")
        git_in(rn["c2"].parent, "commit", "-q", "-m", "inside")
        ok, w = why_of(rn["c2"], "committed in another git repository than this script's")
        check(ok, f"... committed in a repository of its own inside this one: NOT the pre-registered procedure {w}")
    # (5) photo folders other than the station's, and constants changed in the running process
    ro5 = GcpfreeRepo(d / "roots")

    def other_roots(info):
        info["roots"] = ["/tmp/rendered/Chelsea_pics"]
        info["constants"]["ROOTS"] = info["roots"]
        info["constants_sha256"] = G.constants_sha256(info["constants"])
    with ro5:
        r5 = ro5.run("20261014T010000Z", res_nd(), "2026-10-14 01:00:00", edit_info=other_roots)
        ok, w = why_of(r5["c2"], "its photo folders are not the station's", "constants are not the pinned ones")
        check(ok, f"(5) photo folders not the station's (ROOTS patched in the running process): a test of the "
                  f"procedure, never INDEPENDENT {w}")
    # (6) the camera position, (7) the pan change line, (8) a small tilt: edited BEFORE the first commit, report and
    # run.json rewritten to agree
    def report_b(k, delta):
        def edit(js):
            js["B"][k] += delta
        return edit

    def pan_line(t):
        return re.sub(r"pan change ([+-][0-9.]+) deg", lambda m: f"pan change {float(m.group(1)) + 1:+.3f} deg", t)
    for name, edit, words in (("position", dict(edit_b=shifted("x", 0.5), edit_report=report_b(0, 0.5)),
                               "its camera position is 0.500 m from"),
                              ("pan", dict(edit_b=pan_line), "its 'pan change' line says"),
                              ("tilt", dict(edit_b=shifted("tilt", 0.05), edit_report=report_b(4, 0.05)), None)):
        rx = GcpfreeRepo(d / name)
        with rx:
            r_ = rx.run("20261014T010000Z", res_nd(), "2026-10-14 01:00:00", commit=False)
            reseal(r_["c2"].parent, "c2", **edit)
            rx.commit("run (edited before its commit)")
            if words:
                ok, w = why_of(r_["c2"], words)
                check(ok, f"({name}) edited before the first commit, report and run.json rewritten: '{words}': NOT the "
                          f"pre-registered procedure {w}")
            else:
                check(sp.procedure_check(r_["c2"])["pinned"],
                      "(tilt) a 0.05 deg tilt edit BEFORE the first commit, report and run.json rewritten: NOT seen "
                      "here (the 'rotated by' magnitude barely moves): only the B sha256 written down when the run "
                      "ended shows it (stated in the docs)")
    # (9) the pre-registration committed again after the run; (10) a second delivery run without Aug photos
    rq = GcpfreeRepo(d / "prereg_again")
    with rq:
        r9 = rq.run("20261014T010000Z", res_nd(), "2026-10-14 01:00:00")
        r10 = rq.run("20261014T020000Z", res_nd(), "2026-10-14 02:00:00")
        ok, w = why_of(r10["c2"], "it is the second delivery run of this repository")
        check(ok and sp.procedure_check(r9["c2"])["pinned"],
              f"(10) a second delivery run without Aug photos: NOT the delivered B (the first run's is) {w}")
        rq.prereg.write_text(rq.prereg.read_text() + "# a note added after the run\n")
        rq.commit("pre-registration: a note")
        ok, w = why_of(r9["c2"], "not the file its run recorded", "was committed again after the run")
        check(ok, f"(9) the pre-registration changed and committed again after the run: NOT the pre-registered "
                  f"procedure {w}")
    # (11) the first delivery run's record removed (from the working tree, then in a commit), a later run then
    # passing for the first: never the delivered B
    rg = GcpfreeRepo(d / "removed")
    with rg:
        r11 = rg.run("20261014T010000Z", {"c1": gcpfree_res("c1", (0, 0, 0), "NOT DELIVERED"),
                                          "c2": gcpfree_res("c2", c2, "NOT DELIVERED")}, "2026-10-14 01:00:00")
        shutil.rmtree(str(r11["c2"].parent))
        r12 = rg.run("20261014T020000Z", res_nd(), "2026-10-14 02:00:00", commit=False)
        git_in(rg.root, "add", "--", str(r12["c2"].parent))         # the new run only: the removal not committed
        git_in(rg.root, "commit", "-q", "-m", "gcpfree: run 20261014T020000Z")
        ok, w = why_of(r12["c2"], "the record of run(s) 20261014T010000Z was committed")
        check(ok, f"(11) the first run's folder deleted, a second run then the only one there: NOT the delivered B {w}")
        rg.commit("the first run removed")
        ok2, w2 = why_of(r12["c2"], "the record of run(s) 20261014T010000Z was committed")
        check(ok2 and sp.gcpfree_run_records(rg.dir / "runs") == {"removed": ["20261014T010000Z"], "changed": []},
              f"... and with its removal committed: still NOT (git keeps the record) {w2}")
    # (12) the first delivery run's record rewritten as a trial (so a later run would be the first): never the
    # delivered B either
    rh = GcpfreeRepo(d / "rewritten")
    with rh:
        r13 = rh.run("20261014T010000Z", {"c1": gcpfree_res("c1", (0, 0, 0), "NOT DELIVERED"),
                                          "c2": gcpfree_res("c2", c2, "NOT DELIVERED")}, "2026-10-14 01:00:00")
        rj = r13["c2"].parent / "run.json"
        rjs = json.loads(rj.read_text())
        rjs["eligible"] = False
        rj.write_text(json.dumps(rjs, indent=1) + "\n")
        rh.commit("run 1 made a trial")
        r14 = rh.run("20261014T020000Z", res_nd(), "2026-10-14 02:00:00")
        pc14 = sp.procedure_check(r14["c2"])
        ok, w = why_of(r14["c2"], "the record (run.json) of run(s) 20261014T010000Z was changed after its commit")
        check(ok and pc14["delivery_runs_here"] == ["20261014T020000Z"],
              f"(12) the first run's record rewritten as a trial and committed, a later run then the first delivery run "
              f"there: NOT the delivered B {w}")


def test_gcpfree_turn_line(d):
    print("a pinned B's 'X.yaml rotated by N deg' line: X is its report's reference, here with the content the report "
          "recorded, and N is true of the two files' numbers; else NOT the pre-registered procedure")
    import gcpfree_eo as G
    d = Path(d) / "gcpfree_turn"
    cal = HERE / "calibration"
    ref_f = G.REF_EO_FILE["jan"].format(cam="c2")
    repo = GcpfreeRepo(d / "repo")
    with repo:
        b0 = repo.run("20261014T030000Z", {"c2": gcpfree_res("c2", (-13.4, -0.5, 0.3))}, "2026-10-14 03:00:00",
                      edit_info=lambda info: info.update(camera="both"))["c2"]
        run = b0.parent
        m = sp.ROTATED_REF.search(" ".join(sp.eo_comment_lines(b0)))
        pc = sp.procedure_check(b0, [str(cal)])
        check(pc["pinned"] and m and m.group(1) == ref_f and 13.0 < float(m.group(2)) < 14.0,
              f"as gcpfree_eo.py wrote it: '{m.group(0) if m else '?'}', true of the numbers: pinned"
              + ("" if pc["pinned"] else f" -- NOT: {pc['why']}"))

        def forged(name, edit_b=None, ref_beside=None, edit_report=None):
            """A copy of the run, edited, its report and run.json made to agree with the edited B (a consistent
            forgery: only the turn line, the reference named, or the reference's content tells)."""
            f = d / "forged" / name
            shutil.copytree(str(run), str(f))
            if ref_beside:
                (f / ref_f).write_text(ref_beside)
            return reseal(f, "c2", edit_b=edit_b, edit_report=edit_report)
        turn = m.group(2)
        other = "CACO03_c2_20240827_EO.yaml"
        jan_tilted = re.sub(r"^tilt: (\S+)$", lambda t: f"tilt: {float(t.group(1)) + 1.0:.4f}",
                            (cal / ref_f).read_text(), flags=re.M)
        cases = [("angle", dict(edit_b=lambda t: t.replace(f"rotated by {turn} deg",
                                                           f"rotated by {float(turn) + 1:.3f} deg")),
                  "but its numbers are turned"),
                 ("other_ref", dict(edit_b=lambda t: t.replace(f"{ref_f} rotated by", f"{other} rotated by")),
                  f"names {other}, its report's reference is {ref_f}"),
                 ("ref_content", dict(ref_beside=jan_tilted), "is not the reference its report says it was turned from"),
                 ("ref_missing", dict(edit_b=lambda t: t.replace(f"{ref_f} rotated by",
                                                                 "CACO03_c2_20250199_EO.yaml rotated by"),
                                      edit_report=lambda js: js.update(b_reference_eo="CACO03_c2_20250199_EO.yaml")),
                  "is not here")]
        for name, kw, words in cases:
            b = forged(name, **kw)
            pcx = sp.procedure_check(b, [str(cal)])
            lab = sp.held_check(sp.calibration_origin(b, [str(cal)]), sp.blind_of("2024-10-23"))[0]
            msg = sp.procedure_variant_problem({"cams": {"c2": {"eo_given": str(b), "eo_path": b}}},
                                               SimpleNamespace(calibration_variant="gcpfree", calibration=str(cal)))
            check(not pcx["pinned"] and any(words in w for w in pcx["why"]) and lab == "PARTLY-CIRCULAR"
                  and "these are not (or not yet) the B files it delivered here" in msg,
                  f"B {name} (report and run.json resealed): '{words}': NOT the pre-registered procedure, never "
                  f"INDEPENDENT, not built as gcpfree")


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


def end_on_signals():
    """SIGTERM and SIGHUP end the run like Ctrl-C (SystemExit), so main()'s finally removes the test folder: a
    killed run would otherwise leave its canary files, named like the held survey, in /tmp (review r1)."""
    import signal

    def _end(signum, frame):
        raise SystemExit(f"test run ended by signal {signum}: its folder is removed")
    for sig in ("SIGTERM", "SIGHUP"):
        if hasattr(signal, sig):
            try:
                signal.signal(getattr(signal, sig), _end)
            except (ValueError, OSError):
                pass


def test_pointing_selftest_skips(d):
    print("pointing_check.py --self-test without a working scipy: SKIPPED, exit 0 (also with a scipy that is installed "
          "but does not import)")
    stub = Path(d) / "broken_scipy" / "scipy"
    stub.mkdir(parents=True)
    (stub / "__init__.py").write_text("raise ImportError('a broken scipy (test stub)')\n")
    r = subprocess.run([sys.executable, str(HERE / "pointing_check.py"), "--self-test"], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d),
                       env=dict(os.environ, PYTHONPATH=str(stub.parent)))
    check(r.returncode == 0 and "SKIPPED (nothing checked)" in r.stdout and "Traceback" not in r.stdout,
          f"a scipy found but broken: SKIPPED, exit {r.returncode} (find_spec alone would have found it and then "
          f"failed in view_reproject)")


def test_killed_run_cleans_up(d):
    print("a test run ended by SIGTERM removes its own folder (its canaries named like the held survey with it)")
    if os.environ.get("SP_TEST_CHILD"):
        return
    import signal
    pr = subprocess.Popen([sys.executable, str(HERE / "test_survey_products.py")], stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, universal_newlines=True, env=dict(os.environ, SP_TEST_CHILD="1"),
                          cwd=str(d))
    first, folder = "", None
    try:
        for _ in range(40):
            first = pr.stdout.readline()
            m = re.match(r"test folder: (\S+)", first)
            if m:
                folder = Path(m.group(1))
                break
        time.sleep(1.0)                          # it has started its checks (and written files there)
        pr.send_signal(signal.SIGTERM)
        rest = pr.communicate(timeout=120)[0]
    finally:
        if pr.poll() is None:
            pr.kill()
    check(folder is not None and pr.returncode != 0 and not folder.exists() and "signal" in rest,
          f"SIGTERM to a run of this test: it stops (exit {pr.returncode}) and its folder {folder} is gone")


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
    end_on_signals()                    # a kill (SIGTERM/SIGHUP) ends the run through the finally below: the folder,
    #                                     with its canaries named like the held survey, is removed (review r1)
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
        test_oct_pointing_pattern(d)
        test_storm_and_wave_caveats(d)
        test_day_test_record(d)
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
        test_summary_built(d)
        test_fix2(d)
        test_fix3(d)
        test_dem_two_days(d)
        test_filter_leaves_nothing_old(d)
        test_setup_share(d)
        test_c_change(d)
        test_blind_and_calibration(d)
        test_gcpfree_b(d)
        test_gcpfree_b_build(d)
        test_review2_nits(d)
        test_gcpfree_where(d)
        test_gcpfree_turn_line(d)
        test_fit_setup_flat(d)
        with unblinded():            # the GCP path as it runs once Chris has reported
            test_oct_path(d)
        test_oct_transfer_caveat(d)
        with unblinded():
            test_gcp_files_any_format(d)
        test_merge_level_source(d)
        test_dry_run_needs_adcp()
        test_one_setup_coef()
        test_pointing_selftest_skips(d)
        test_killed_run_cleans_up(d)
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

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
asserted on. Runs in ~5 s.
"""

import os
import io
import sys
import csv
import shutil
import argparse
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
                                      "envelope_path": Path(env) if env else None, "envelope_survey": env}}}


def test_labels(d):
    print("labels can only get worse")
    d = Path(d) / "labels"
    d.mkdir()
    args = SimpleNamespace(setup_coef=0.037)
    jan = {"name": "jan_lidar", "survey_type": "dsm", "survey_date": "2025-01-23", "label": "INDEPENDENT", "why": ""}
    lab, why = sp.honest_label(jan, plan_for(d), args, Path("2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif"))
    check(lab == "INDEPENDENT", "Jan lidar, envelope from the Mar lidar, GCP calibration: INDEPENDENT")
    lab, why = sp.honest_label(jan, plan_for(d, env="2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif"), args,
                               Path("2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif"))
    check(lab == "PARTLY-CIRCULAR", "envelope placed with the compared lidar: PARTLY-CIRCULAR")
    fitted = ("# fit_eo_to_survey.py: fitted to 2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif\n"
              + (HERE / "calibration" / "CACO03_c1_20250123_EO.yaml").read_text())
    lab, why = sp.honest_label(jan, plan_for(d, eo="CACO05_c1_2025-01-18_to_2025-01-23_lidar_EO.yaml",
                                             eo_text=fitted), args,
                               Path("2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif"))
    check(lab == "CIRCULAR", "pointing fitted to the compared lidar: CIRCULAR")
    rtk = {"name": "rtk_2026-09-29", "survey_type": "points", "survey_date": "2026-09-29",
           "label": "INDEPENDENT", "why": ""}
    p = plan_for(d, era="live", env=None, eo="CACO05_c1_20251113_EO-CV.yaml")
    lab, why = sp.honest_label(rtk, p, args, Path("2026-09-29_Marconi_Checkshots.csv"))
    check(lab == "CIRCULAR" and "setup" in why[0], "RTK with C = 0.037 (fitted to it): CIRCULAR")
    # a C of unknown origin is never INDEPENDENT; a C declared fitted to no survey keeps the label
    for c in (0.035, 0.05, 0.0371):
        lab, why = sp.honest_label(rtk, p, SimpleNamespace(setup_coef=c), Path("2026-09-29_Marconi_Checkshots.csv"))
        check(lab == "PARTLY-CIRCULAR" and "unknown origin" in why[0], f"RTK with C = {c} of unknown origin: PARTLY-CIRCULAR")
    lab, why = sp.honest_label(rtk, p, SimpleNamespace(setup_coef=0.035, setup_fitted_to="none:repeat crossings"),
                               Path("2026-09-29_Marconi_Checkshots.csv"))
    check(lab == "INDEPENDENT", "RTK with a C declared fitted to no survey (--setup-fitted-to none:...): INDEPENDENT")
    # the RTK file under another name is recognised by its content
    real = Path(d) / "rtk_copy.csv"
    rtk_src = Path("/tmp/claude-0/-home-user-caco05-waterline/1e668ac9-7b7a-5204-b62b-1a0cde16501d/scratchpad/"
                   "survey_inputs/2026-09-29_Marconi_Checkshots.csv")
    if rtk_src.exists():
        shutil.copy(rtk_src, real)
        lab, why = sp.honest_label(dict(rtk, survey_date=""), p, args, real)
        check(lab == "CIRCULAR" and "same content" in why[0], "the RTK shots renamed (rtk_copy.csv): still CIRCULAR, by sha256")
    else:
        print("  skip  renamed-RTK check (the real RTK file is not on this computer)")
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
    sp.honest_label(rtk, p, args, real if real.exists() else Path("2026-09-29_Marconi_Checkshots.csv"), verdict=v)
    check(v.label == "CIRCULAR", "the same shots (content or name) cannot be overruled")
    # a C declared fitted to a survey file: CIRCULAR against that file
    lab, why = sp.honest_label(rtk, p, SimpleNamespace(setup_coef=0.05, setup_fitted_to=str(other),
                                                       survey_dirs=[]), other)
    check(lab == "CIRCULAR", "a C declared fitted to this very file: CIRCULAR")
    gcp = {"name": "oct_gcps", "survey_type": "points", "survey_date": "2024-10-23", "label": "INDEPENDENT", "why": ""}
    p = plan_for(d, eo="CACO03_c1_20241023_EO.yaml")
    lab, why = sp.honest_label(gcp, p, args, Path("2024-10-23_Marconi_Extrinsic_Targets_c1_xyz.csv"))
    check(lab == "PARTLY-CIRCULAR", "GCPs of the calibration day, even if the table said INDEPENDENT")
    lab, why = sp.honest_label(gcp, p, args, Path("2024-10-23_control.csv"))
    check(lab == "PARTLY-CIRCULAR", "a point survey of the calibration day under any name ('_control.csv'): PARTLY-CIRCULAR")


def test_round2(d):
    print("carried pointings, earlier checks, horizon offsets, missing parts, DEM settings of the build")
    d = Path(d) / "round2"
    d.mkdir()
    args = SimpleNamespace(setup_coef=0.037)
    jan = {"name": "jan_lidar", "survey_type": "dsm", "survey_date": "2025-01-23", "label": "INDEPENDENT", "why": ""}
    mar = dict(jan, name="mar_lidar", survey_date="2025-03-06")
    janf, marf = (Path("2025005FA_Marconi_%s_YSMP_Lidar_DSM_25cm.tif" % m) for m in ("Jan", "Mar"))
    base = (HERE / "calibration" / "CACO03_c1_20250123_EO.yaml").read_text()
    # apply_pointing_correction.py carries a lidar fit's change to another period (*_corr_EO.yaml)
    (d / "CACO05_c1_2025-01-18_to_2025-01-22_lidar_EO.yaml").write_text(
        "# fit_eo_to_survey.py: fitted to 2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif\n" + base)
    corr = ("# apply_pointing_correction.py: CACO04_c1_20250219_EO.yaml + (+0.120, -0.310, +0.020) deg from "
            "CACO05_c1_2025-01-18_to_2025-01-22_lidar_EO.yaml; 2024-10-01..2025-01-17\n" + base)
    p = plan_for(d, eo="CACO05_c1_20250219_corr_EO.yaml", eo_text=corr)
    v = sp.Verdict(dict(jan, label_override_reason="anything"))
    sp.honest_label(jan, p, args, janf, verdict=v)
    check(v.label == "CIRCULAR" and any("carries a correction fitted to" in r for r in v.reasons),
          "a *_corr_EO.yaml carried from the Jan lidar fit, against the Jan lidar: CIRCULAR, not overrulable")
    pm = plan_for(d, eo="CACO05_c1_20250219_corr_EO.yaml", eo_text=corr, env=janf.name)   # Mar: envelope from Jan
    lab, why = sp.honest_label(mar, pm, args, marf)
    check(lab == "INDEPENDENT", "... against the Mar lidar (another survey): INDEPENDENT")
    (d / "CACO05_c1_2025-01-18_to_2025-01-22_lidar_EO.yaml").unlink()
    v = sp.Verdict(dict(mar, label_override_reason="the fit was to the Jan lidar (operator's notes)"))
    sp.honest_label(mar, pm, args, marf, verdict=v)
    check(v.label == "INDEPENDENT" and v.overridden and "cannot be resolved" in v.overridden[0],
          "the fit it was carried from is not here: CIRCULAR as a rule, overruled only by a written reason")
    lab, why = sp.honest_label(mar, pm, args, marf)
    check(lab == "CIRCULAR", "... and without that reason: CIRCULAR")
    # the earlier check of 2025-01-23 on the REAL photos leads the caveats, per camera, with what to
    # expect now. pointing_fix_figure.py printed lidar - water level = c1 +0.46, c2 +0.15: the lines
    # read LOW (waterline - lidar c1 -0.46, c2 -0.15), the sign a line without setup must have.
    plan = {"era": "adcp", "date": "2025-01-23", "first": "2025-01-18", "last": "2025-01-23", "out": d,
            "surveys": [dict(jan, path=str(janf))],
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
    check("c1 ~-0.26 m" in c0 and "c2 ~+0.05 m" in c0 and "c1 -0.18, c2 +0.13 m if the" in c0
          and "0.71 x its setup" in c0,
          "... expecting c1 ~-0.26 / c2 ~+0.05 m with C = 0.037: each line keeps 0.71 x its setup, paired per frame "
          "on the 2026 RTK (-0.18 / +0.13 m if the whole setup counted)")
    check("PARTLY-CIRCULAR: its search envelope was placed with this same Jan 2025 lidar" in c0,
          "... the earlier check carries its own label: PARTLY-CIRCULAR (envelope placed with the same lidar)")
    under = [c for c in cav if c.startswith("C = 0.037 under-corrects")]
    check("do not subtract it again" in c0 and under and "already inside the EARLIER CHECK expectation" in under[0]
          and "~0.07 m" in under[0] and "0.046" in under[0],
          "... the under-correction (re-projection ~0.07 m + C mismatch, per-frame C 0.046) is said to be inside "
          "that expectation, never to be counted twice")
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
    named = dict(jan, path=str(janf.name))
    nofile = {"name": "oct_gcps", "survey_type": "points", "survey_date": "2024-10-23", "label": "PARTLY-CIRCULAR",
              "why": "", "path": ""}
    miss = sp.missing_parts({"out": d, "date": "2025-01-23", "cams": {}, "surveys": [named, nofile]}, {},
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
    # a re-export of the RTK the station's C was fitted to, entered under another date: same shots
    rtk_src = Path("/tmp/claude-0/-home-user-caco05-waterline/1e668ac9-7b7a-5204-b62b-1a0cde16501d/scratchpad/"
                   "survey_inputs/2026-09-29_Marconi_Checkshots.csv")
    if rtk_src.exists():
        lines = rtk_src.read_text().splitlines()
        lines[1] = lines[1].replace(",", ", ", 1)
        rex = d / "rtk_reexport.csv"
        rex.write_text("\n".join(lines) + "\n")
        rtk = {"name": "rtk_other", "survey_type": "points", "survey_date": "2026-09-30", "label": "INDEPENDENT",
               "why": ""}
        pl = plan_for(d, era="live", env=None, eo="CACO05_c1_20251113_EO-CV.yaml")
        lab, why = sp.honest_label(rtk, pl, SimpleNamespace(setup_coef=0.037, survey_dirs=[str(rtk_src.parent)]), rex)
        check(lab == "CIRCULAR" and "re-exported" in why[0], "a byte-changed re-export under another date: CIRCULAR "
                                                               "by its coordinates")
    else:
        print("  skip  re-exported-RTK check (the real RTK file is not on this computer)")


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
    lab, why = sp.honest_label(rtk, p, later, Path("2026-09-29_Marconi_Checkshots.csv"))
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
    a = dict(config=str(sp.DEFAULT_CONFIG), surveys=str(sp.DEFAULT_SURVEYS), force_disabled=False, window=None,
             utc_hours=None, calibration=str(HERE / "calibration"), survey_dirs=sp.DEFAULT_SURVEY_DIRS,
             output_root="/tmp/none", setup_coef=sp.SETUP_COEF)
    a.update(kw)
    return SimpleNamespace(**a)


def test_config(d):
    print("the configuration tables")
    cfg = sp.read_table(sp.DEFAULT_CONFIG, sp.CONFIG_FIELDS)
    srv = sp.read_table(sp.DEFAULT_SURVEYS, sp.SURVEY_FIELDS)
    for date in ("2024-10-23", "2025-03-19"):
        plan, out = quiet(sp.build_plan, date, cfg, srv, cfg_args())
        check(plan is None and "DISABLED" in out, f"{date} is refused (disabled at the user's request)")
        plan, out = quiet(sp.build_plan, date, cfg, srv, cfg_args(force_disabled=True))
        check(plan is not None, f"{date} builds with --force-disabled")
    plan, _ = quiet(sp.build_plan, "2025-03-06", cfg, srv, cfg_args())
    check(list(plan["cams"]) == ["c2"] and plan["cams"]["c2"]["station"] == "CACO04", "2025-03-06: c2 only, CACO04")
    plan, _ = quiet(sp.build_plan, "2026-09-29", cfg, srv, cfg_args())
    check(plan["first"] == "2026-09-26" and plan["last"] == "2026-10-02" and plan["era"] == "live"
          and [s["name"] for s in plan["surveys"]] == ["rtk_2026-09-29"],
          "2026-09-29 (the user first called it 27 Sep): window 26 Sep - 2 Oct, the 29 Sep RTK")
    check(not [r for r in cfg if r["date"] == "2026-09-27"], "no 2026-09-27 product left in the table")
    rtk_row = [s for s in srv if s["name"] == "rtk_2026-09-29"][0]
    check(rtk_row["label"] == "INDEPENDENT" and "CIRCULAR" in rtk_row["why"] and "C = 0.037" in rtk_row["why"]
          and "calm low-tide RTK" in rtk_row["why"] and "+1.2 m" in rtk_row["why"],
          "the 2026 RTK row: INDEPENDENT as a table label, saying the setup rule makes it CIRCULAR with C = 0.037, "
          "where the check is and what an independent check needs")
    rtk_file = Path(sp.DEFAULT_SURVEY_DIRS[0]) / rtk_row["path"]
    pl, _ = quiet(sp.build_plan, "2026-09-29", cfg, srv, cfg_args())
    lab, why = sp.honest_label(rtk_row, pl, SimpleNamespace(setup_coef=sp.SETUP_COEF), rtk_file)
    check(lab == "CIRCULAR" and "setup coefficient C = 0.037 was fitted to this survey" in why[0],
          "... a default build (C = 0.037): CIRCULAR by the setup rule (the file's name, or its content when here)")
    lab, why = sp.honest_label(rtk_row, pl, SimpleNamespace(setup_coef=0.0), rtk_file)
    check(lab == "INDEPENDENT" and not why, "... a C = 0 build (fitted to nothing): INDEPENDENT")
    plan, _ = quiet(sp.build_plan, "2025-01-23", cfg, srv, cfg_args())
    check(plan["last"] == "2025-01-23" and all(c["station"] == "CACO03" for c in plan["cams"].values()),
          "2025-01-23: CACO03 only, 24 Jan (the re-set) outside the window")
    for r in cfg:
        if r["era"] == "live":
            continue
        for s in srv:
            if s["date"] == r["date"] and s["path"]:
                check(Path(s["path"]).name != Path(r["envelope_survey"]).name,
                      f"{r['date']} {r['camera']}: envelope {r['envelope_survey'][:28]}... is not the compared survey")
    for r in cfg:
        eo = HERE / "calibration" / r["eo_file"]
        check(eo.exists(), f"{r['date']} {r['camera']}: calibration {r['eo_file']} exists")
    labels = {s["name"]: s["label"] for s in srv}
    check(labels["oct_gcps"] == "PARTLY-CIRCULAR" and not [s for s in srv if s["name"] == "oct_gcps"][0]["path"],
          "Oct GCPs: PARTLY-CIRCULAR, path empty until the file arrives")
    # The GCP files arrive later, one per camera, under the survey date's name: picked up by
    # themselves and merged, a target both cameras saw counted once.
    oct_ = [s for s in srv if s["name"] == "oct_gcps"][0]
    sdir = Path(d) / "gcp_arrive"
    sdir.mkdir()
    a = cfg_args(survey_dirs=[str(sdir)])
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
        check("under-corrects" in text, f"{era}: the setup under-correction caveat is there too")
        check("none on 2026-09-26" in text, f"{era}: the window days without waterlines are listed")
        check(("GNSS-R" in text and "surf zone" in text) == (era == "adcp"),
              f"{era}: C's still-water reference (GNSS-R, surf zone) stated for the ADCP era only")
        check("NOMINAL" in text and "6348" in text, f"{era}: the GeoTIFF CRS tag is said to be nominal")
    cav = sp.collect_caveats(plan, SimpleNamespace(setup_coef=0.05), {}, forcing, {})
    check(not any("under-corrects" in c for c in cav) and any("WITHOUT saying where" in c for c in cav),
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
    check(st["software"]["scripts"].get("survey_products.py") and "commit" in st["software"],
          "each stamp records the commit and the sha256 of the scripts that built it")


def test_keep_moved_days(d):
    print("--keep-moved-days is applied when the detection runs, not stored")
    res = {"different_days": ["2025-01-24"], "leave_out": True}
    check(sp.skip_days_of(res, SimpleNamespace(keep_moved_days=False)) == ["2025-01-24"], "default: the day is left out")
    check(sp.skip_days_of(res, SimpleNamespace(keep_moved_days=True)) == [], "--keep-moved-days on a rerun: kept")
    check(sp.skip_days_of(dict(res, leave_out=False), SimpleNamespace(keep_moved_days=False)) == [],
          "most days differ: nothing left out automatically")


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


def test_py38():
    print("Python 3.8 syntax")
    import ast
    for f in ("survey_products.py", "detect_original_view.py", "test_survey_products.py", "estimate_eo_rotation.py",
              "survey_compare.py"):
        src = (HERE / f).read_text()
        try:
            ast.parse(src, feature_version=(3, 8))
        except TypeError:                        # Python 3.7 has no feature_version
            ast.parse(src)
        check(True, f"{f} parses as Python 3.8")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="write the outputs here and keep them")
    args = ap.parse_args()
    d = args.keep or tempfile.mkdtemp(prefix="survey_products_test_")
    Path(d).mkdir(parents=True, exist_ok=True)
    try:
        test_setup_identical(d)
        test_collect_photos(d)
        test_labels(d)
        test_setup_in_use(d)
        test_merge_drops_no_setup(d)
        test_keep_moved_days(d)
        test_live_rows(d)
        test_config(d)
        test_frames_per_day(d)
        test_caveats(d)
        test_stamps(d)
        test_round2(d)
        test_fix1(d)
        test_py38()
    finally:
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    print("all survey_products tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

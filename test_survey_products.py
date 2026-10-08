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
    lab, why = sp.honest_label(rtk, p, SimpleNamespace(setup_coef=0.035), Path("2026-09-29_Marconi_Checkshots.csv"))
    check(lab == "INDEPENDENT", "RTK with a C fitted elsewhere: the table's label")
    gcp = {"name": "oct_gcps", "survey_type": "points", "survey_date": "2024-10-23", "label": "INDEPENDENT", "why": ""}
    p = plan_for(d, eo="CACO03_c1_20241023_EO.yaml")
    lab, why = sp.honest_label(gcp, p, args, Path("2024-10-23_Marconi_Extrinsic_Targets_c1_xyz.csv"))
    check(lab == "PARTLY-CIRCULAR", "GCPs of the calibration day, even if the table said INDEPENDENT")


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
        plan = {"era": era, "date": "2026-09-27", "first": "2026-09-24", "last": "2026-09-30", "out": out,
                "surveys": [], "cams": {"c1": {"first": "2026-09-24", "last": "2026-09-30"}}}
        forcing = {"waves": {"setup_in_window": {"median_m": 0.3, "daytime_median_m": 0.28}}}
        cav = sp.collect_caveats(plan, a, {}, forcing, {"frames_per_day": {"c1": {"2026-09-29": 3}}})
        text = "\n".join(cav)
        first = "Live era" if era == "live" else "ADCP datum"
        check(first in text, f"{era}: its own era caveat is there")
        check("under-corrects" in text, f"{era}: the setup under-correction caveat is there too")
        check("none on 2026-09-24" in text, f"{era}: the window days without waterlines are listed")
    cav = sp.collect_caveats(plan, SimpleNamespace(setup_coef=0.05), {}, forcing, {})
    check(not any("under-corrects" in c for c in cav), "a C fitted elsewhere: no RTK under-correction caveat")


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


def test_py38():
    print("Python 3.8 syntax")
    import ast
    for f in ("survey_products.py", "detect_original_view.py", "test_survey_products.py"):
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
        test_config(d)
        test_frames_per_day(d)
        test_caveats(d)
        test_stamps(d)
        test_py38()
    finally:
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    print("all survey_products tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Tests of weekly_frames.py (the weekly DEM and waterline frames and their videos).

    python3 test_weekly_frames.py

Plain script like the other tests (no pytest); exits 1 when a check fails. Python 3.8+ (the NUC's 3.8:
numpy 1.17, matplotlib 3.3). Everything runs in a temporary BASE: a copy of the real waterline_timex_cron.sh
(its DEM options), a synthetic ground file, and the REAL dem_from_contours.py for the weeks of T2; a fake
runner elsewhere (speed). ffmpeg: a fake one on PATH, and, when a real ffmpeg is installed, one real tiny MP4
whose frame count ffprobe checks. No network. The blind date's week is built from rows in the ground file
only; nothing named for that date is created or opened (an audit hook watches every open).

  T1  weeks (UTC Monday..Sunday), titles, statuses, settings, the blind list, read_points
  T2  the REAL dem_from_contours.py: the exact cron options, the week window, grids never in archive/dems
  T3  statuses over days: current -> provisional -> final (coverage, max wait), the empty week, a failed
      build, the daily step's once-per-day stamp (a second run changes nothing)
  T4  the frozen look: a changed setting is a NOTE (nothing redrawn), respec redraws every frame
  T5  backfill: dry run writes nothing, the blind week is skipped (no file named for it opened or made),
      resumable after a stop
  T6  the videos: concat order (weeks and storm frames by date, change after DEM), ffmpeg's arguments,
      rebuilt only on change, ffmpeg missing (a NOTE, frames still made, rc 0), a real tiny MP4
  T7  locks: the cron's lock held, a cron run due, another weekly run busy; the CLI end to end
Every command's writes are checked to land under archive/weekly/ or logs/.
"""

import os
import sys
import csv
import json
import time
import shutil
import tempfile
import subprocess
from pathlib import Path
from datetime import datetime, timedelta, timezone

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import weekly_frames as wf             # noqa: E402
import storm_watch as sw               # noqa: E402

UTC = timezone.utc
FAILED = []
REAL_FFMPEG = shutil.which("ffmpeg")
REAL_FFPROBE = shutil.which("ffprobe")
GHEAD = ("source_file,camera,capture_time_utc,capture_epoch,pixel_column,easting_utm19,northing_utm19,"
         "tide_elevation_navd88,offshore_hs_m,offshore_tp_s,setup_correction_m,beach_elevation_navd88")
W_TIDE = 2 * np.pi / (12.42 * 3600)
BLIND = "2024-10-23"
AUDIT = {"on": False, "hits": []}


def check(cond, msg):
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        FAILED.append(msg)


def audit(event, args):
    if AUDIT["on"] and event == "open" and args and isinstance(args[0], (str, bytes, os.PathLike)):
        p = os.fsdecode(args[0])
        if wf.blind_named(Path(p).name) or any(wf.blind_named(x) for x in Path(p).parts):
            AUDIT["hits"].append(p)


if hasattr(sys, "addaudithook"):
    sys.addaudithook(audit)


def T(s):
    return datetime.fromisoformat(s.replace("Z", "+00:00"))


def make_base(tmp, name="base", settings=None):
    base = Path(tmp) / name
    (base / "archive").mkdir(parents=True)
    (base / "logs").mkdir()
    shutil.copy(str(HERE / "waterline_timex_cron.sh"), str(base))
    (base / ".timex_cron.lock").write_text("")
    if settings is not None:
        (base / "weekly_frames_settings.conf").write_text(settings)
    return base


def rows_for(d, rng, hs=0.8, cams=("c1", "c2"), frames=25, shift=0.0, first_hour=11):
    """A day of a small station: frames every 30 min from first_hour UTC, 20 points per camera line."""
    out = []
    t0 = datetime.fromisoformat(d + "T00:00:00+00:00") + timedelta(hours=first_hour)
    for k in range(frames):
        t = t0 + timedelta(minutes=30 * k)
        e_ = int(t.timestamp())
        tide = 0.6 + 1.2 * np.sin(W_TIDE * e_)
        for cam in cams:
            x0 = 420001.0 if cam == "c1" else 420021.0
            for col in range(20):
                e = x0 + 2 * col
                xr = int(round((tide + 0.05 + 1.0) / 0.2))
                z = tide + 0.05 + shift + rng.normal(0, 0.02)
                h = hs(d, k) if callable(hs) else hs
                out.append(f"{e_}.{cam}.timex,{cam},{t:%Y-%m-%dT%H:%M:%S}+00:00,{e_},{col},{e:.2f},"
                           f"{4638001.0 + 2 * xr:.2f},{z:.4f},{h:.2f},8.0,0.05,{z:.4f}")
    return out


def dates(a, b):
    d = datetime.fromisoformat(a)
    out = []
    while d <= datetime.fromisoformat(b):
        out.append(d.strftime("%Y-%m-%d"))
        d += timedelta(days=1)
    return out


def write_ground(base, days, seed=1, **kw):
    rng = np.random.RandomState(seed)
    rows = []
    for d in days:
        rows += rows_for(d, rng, **kw)
    (base / "contour_points_ground.csv").write_text(GHEAD + "\n" + "\n".join(rows) + "\n")


def tree(base):
    out = {}
    for root, dirs, files in os.walk(str(base)):
        for f in files:
            p = Path(root) / f
            st = p.stat()
            out[str(p.relative_to(base))] = (st.st_size, st.st_mtime_ns)
    return out


def writes_ok(before, after):
    """New or changed files, and whether every one lies under archive/weekly/ or logs/ (and none removed
    elsewhere)."""
    changed = [k for k in after if before.get(k) != after[k]]
    gone = [k for k in before if k not in after]
    bad = [k for k in changed + gone if not (k.startswith("archive/weekly/") or k.startswith("logs/"))
           and k != "contour_points_ground.csv"]          # the test's own ground file
    return changed, bad


def fake_runner(calls, rc_for=None, cells=(12, 25), raise_for=None):
    """dem_from_contours.py: a 12 x 25 grid (or rc_for[monday]); every other command: the real runner."""
    def run(argv, log_path=None):
        if Path(argv[1]).name == "dem_from_contours.py":
            calls.append(list(argv))
            stem = argv[3]
            a, b = argv[argv.index("--start-date") + 1], argv[argv.index("--end-date") + 1]
            if raise_for and a in raise_for:
                raise RuntimeError("the NUC lost power")
            if rc_for and a in rc_for:
                return rc_for[a], "No georectified points matched.\n"
            import dem_change
            g = np.add.outer(np.linspace(-0.6, 1.4, cells[0]), np.zeros(cells[1]))
            for suf, grid in (("_dem", g), ("_spread", np.full(cells, 0.1)), ("_count", np.full(cells, 9.0)),
                              ("_sigma", np.full(cells, 0.02))):
                dem_change.write_grid(f"{stem}{suf}.asc", grid, 420000.0, 4638004.0, 2.0)
            Path(f"{stem}_info.json").write_text(json.dumps({"first_date": a, "last_date": b, "days": 5,
                                                             "window_start": a, "window_end": b,
                                                             "frames": {"c1": 40, "c2": 41}}))
            return 0, "NOTE: fake build\n"
        return wf.subprocess_runner(argv, log_path)
    return run


def run(base, argv, now, runner=None):
    import contextlib
    import io
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = wf.main(list(argv) + ["--base", str(base)], runner=runner, now=T(now))
    out = buf.getvalue()
    return rc, out


def all_pngs(base):
    return sorted(p for p in (base / "archive/weekly").rglob("*.png"))


def wk(base, monday):
    return json.loads((base / "archive/weekly/weeks" / monday / "week.json").read_text())


# --------------------------------------------------------------------------------------------
# T1
# --------------------------------------------------------------------------------------------

def t1():
    print("T1 weeks, titles, statuses, settings, the blind list, read_points")
    check(wf.d_iso(wf.monday_of("2026-11-23")) == "2026-11-23" and wf.d_iso(wf.monday_of("2026-11-29")) == "2026-11-23"
          and wf.sunday_of("2026-11-23") == "2026-11-29" and wf.d_iso(wf.monday_of("2027-01-01")) == "2026-12-28",
          "a week is Monday..Sunday, named by its Monday (Fri 1 Jan 2027 -> Mon 28 Dec 2026)")
    check(wf.week_title("2026-11-16", "CACO-05 Marconi Beach")
          == "CACO-05 Marconi Beach: week of Mon 16 Nov to Sun 22 Nov 2026 (UTC dates)"
          and wf.week_title("2026-12-28", "X") == "X: week of Mon 28 Dec 2026 to Sun 3 Jan 2027 (UTC dates)",
          "titles: 'week of Mon 16 Nov to Sun 22 Nov 2026 (UTC dates)'; across a new year both years")
    D = wf.d_of
    check(wf.decide_status("2026-11-16", D("2026-11-16"), "2026-11-14", 7) == ("current", "")
          and wf.decide_status("2026-11-16", D("2026-11-22"), "2026-11-20", 7) == ("current", "")
          and wf.decide_status("2026-11-09", D("2026-11-16"), "2026-11-14", 7) == ("provisional", "")
          and wf.decide_status("2026-11-09", D("2026-11-16"), "2026-11-15", 7) == ("provisional", "")
          and wf.decide_status("2026-11-09", D("2026-11-17"), "2026-11-16", 7) == ("final", "covered")
          and wf.decide_status("2026-11-09", D("2026-11-21"), "2026-11-14", 7) == ("provisional", "")
          and wf.decide_status("2026-11-09", D("2026-11-22"), "2026-11-14", 7) == ("final", "waited")
          and wf.decide_status("2026-11-09", D("2026-11-22"), None, 7) == ("final", "waited")
          and wf.decide_status("2026-11-09", D("2026-11-16"), None, 0) == ("final", "waited"),
          "statuses: current (today in the week), provisional (ground through the Sunday or before), final when a "
          "later date is covered, or WEEKLY_MAX_WAIT_DAYS after the Sunday")
    import survey_products
    check(sorted(wf.BLIND_DATES) == sorted(survey_products.BLIND_DATES),
          f"BLIND_DATES {list(wf.BLIND_DATES)} equal survey_products.BLIND_DATES")
    forms = ["2024-10-23", "20241023", "2024-Oct-23", "10-23-2024", "Oct_23_2024", "October 23 2024", "23Oct2024",
             "1729699200.Wed.Oct.23_16_00_00.GMT.2024.caco03.c1.timex.jpg"]
    check(all(wf.blind_named("x_" + f + "_y") for f in forms) and not wf.blind_named("2026-11-23_week.png")
          and not wf.blind_named("2024-10-27_week.png"), "blind_named: every name form of the blind date, no other")
    check(wf.blind_in_week("2024-10-21") == BLIND and wf.blind_in_week("2024-10-14") is None
          and wf.blind_in_week("2024-10-28") is None, "blind_in_week: only Mon 21..Sun 27 Oct 2024")
    check(wf.storm_frame_name("2026-11-14", "2026-11-09_fc", "post3") == "2026-11-14_storm_2026-11-09_fc_post3.png",
          "storm_frame_name: <date>_storm_<event>_<part>.png (sorts by its first 10 characters)")
    # settings
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        cfg, raw, notes = wf.load_settings(base)
        check(cfg == {"WEEKLY_ENABLE": True, "WEEKLY_FPS": 2.0, "WEEKLY_MAX_WAIT_DAYS": 7, "WEEKLY_EXTENT": None,
                      "WEEKLY_Z_LIMITS": None, "WEEKLY_CHANGE_LIMIT": 0.5, "WEEKLY_SIZE": [1920, 1080]}
              and not raw and not notes, "no settings file: the defaults (1, 2 fps, 7 days, auto, auto, 0.5 m, 1920x1080)")
        (base / "weekly_frames_settings.conf").write_text(
            "WEEKLY_FPS=4\nWEEKLY_SIZE=1921x1080\nWEEKLY_EXTENT=\"5,1,0,9\"\nWEEKLY_BOGUS=1\n"
            "WEEKLY_Z_LIMITS=\"-1,2.5\"\nWEEKLY_MAX_WAIT_DAYS=x\n")
        cfg, raw, notes = wf.load_settings(base)
        check(cfg["WEEKLY_FPS"] == 4.0 and cfg["WEEKLY_SIZE"] == [1920, 1080] and cfg["WEEKLY_EXTENT"] is None
              and cfg["WEEKLY_Z_LIMITS"] == [-1.0, 2.5] and cfg["WEEKLY_MAX_WAIT_DAYS"] == 7
              and any("WEEKLY_BOGUS" in n and "not known" in n for n in notes)
              and any("WEEKLY_SIZE" in n and "even" in n for n in notes)
              and any("WEEKLY_EXTENT" in n for n in notes) and any("WEEKLY_MAX_WAIT_DAYS" in n for n in notes),
              f"settings: valid values used; an odd size, a reversed extent, a bad number and an unknown key are "
              f"notes with the default kept ({len(notes)} notes)")
    ex = (HERE / "weekly_frames_settings.conf.example").read_text()
    live = sw.var_lines(ex)
    check(not live and all(f"#{k}=" in ex for k in wf.SETTINGS),
          "weekly_frames_settings.conf.example: every key, commented out (no live line)")
    defaults_in_example = {k: sw.var_value(l.split("=", 1)[1]) for l in ex.splitlines() for k in wf.SETTINGS
                           if l.startswith(f"#{k}=")}
    check(all(defaults_in_example[k] == wf.SETTINGS[k][0] for k in wf.SETTINGS),
          "the example shows every default as the script has it")
    # read_points: the UTC date window, blank coordinates, --max-hs, beach else tide
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "g.csv"
        p.write_text(GHEAD + "\n"
                     "a,c1,2026-11-22T23:59:00+00:00,1,0,420001.00,4638001.00,0.5,0.8,8,0.05,0.55\n"
                     "b,c1,2026-11-23T00:30:00+00:00,2,0,420001.00,4638001.00,0.5,0.8,8,0.05,0.56\n"
                     "c,c2,2026-11-24T12:00:00+00:00,3,0,,,0.5,0.8,8,0.05,0.57\n"
                     "d,c2,2026-11-24T13:00:00+00:00,4,0,420003.00,4638003.00,0.7,2.4,8,0.05,0.77\n"
                     "e,c2,2026-11-24T14:00:00+00:00,5,0,420005.00,4638005.00,0.4,,8,,\n"
                     "f,c1,2026-11-30T00:00:00+00:00,6,0,420001.00,4638001.00,0.5,0.8,8,0.05,0.58\n")
        pts = wf.read_points(p, "2026-11-23", "2026-11-29", 1.5)
        check([x["source_file"] for x in pts] == ["b", "e"] and abs(pts[0]["elevation"] - 0.56) < 1e-9
              and abs(pts[1]["elevation"] - 0.4) < 1e-9,
              "read_points: Sunday 22 Nov 23:59Z is last week, Mon 23 Nov 00:30Z (Sun 19:30 EST) this week; blank "
              "coordinates skipped; Hs 2.4 > 1.5 dropped, blank Hs kept; beach elevation, else the tide's")
        pts = wf.read_points(p, "2026-11-23", "2026-11-29", None)
        check([x["source_file"] for x in pts] == ["b", "d", "e"], "read_points without --max-hs keeps the rough frame")


# --------------------------------------------------------------------------------------------
# T2 the real dem_from_contours.py
# --------------------------------------------------------------------------------------------

def t2():
    print("T2 the REAL dem_from_contours.py: the exact cron options, the week window, never archive/dems")
    t0 = time.time()
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        # 2-14 Nov 2026; 4-5 Nov rough (Hs 2.2: left out by --max-hs 1.5)
        write_ground(base, dates("2026-11-02", "2026-11-14"),
                     hs=lambda d, k: 2.2 if d in ("2026-11-04", "2026-11-05") else 0.8)
        (base / "archive/dems").mkdir()
        calls = []

        def recording(argv, log_path=None):
            if Path(argv[1]).name == "dem_from_contours.py":
                calls.append(list(argv))
            return wf.subprocess_runner(argv, log_path)
        before = tree(base)
        rc, out = run(base, ["backfill"], "2026-11-16T13:00:00Z", recording)
        after = tree(base)
        changed, bad = writes_ok(before, after)
        check(rc == 0 and "week 2026-11-02..2026-11-08 FINAL" in out and "week 2026-11-09..2026-11-15 PROVISIONAL" in out,
              f"backfill on Mon 16 Nov: week 2-8 Nov final (the file covers 9 Nov), 9-15 Nov provisional "
              f"(rc {rc})")
        cron, _ = sw.cron_settings(base)
        expect_opts = ["--max-hs", "1.5", "--max-day-offset", "0.15", "--trend-note", "--cell", "2.0", "--min-points",
                       "3", "--max-spread", "0.5"]
        ground = str((base / "archive/weekly/work/ground_window.csv").resolve())
        check(len(calls) == 2 and all(c[0] == sys.executable and c[1] == str(HERE / "dem_from_contours.py")
                                      and c[2] == ground for c in calls),
              "two builds, dem_from_contours.py on the working copy of the ground rows")
        m = "2026-11-02"
        stem = str((base / "archive/weekly/weeks" / m / f"wk_{m}_2026-11-08").resolve())
        check(calls and calls[0][3:] == [stem, "--start-date", m, "--end-date", "2026-11-08"] + expect_opts
              + ["--no-plot"],
              "the week's build: stem wk_<mon>_<sun>, --start-date Monday --end-date Sunday, then exactly the "
              "cron's window-DEM options (--max-hs 1.5 --max-day-offset 0.15 --trend-note --cell 2.0 --min-points 3 "
              "--max-spread 0.5), --no-plot")
        check(all(c[8:-1] == sw.dem_args(cron, True) for c in calls),
              "... the options are storm_watch.dem_args of the cron's VAR= lines (one source, no drift)")
        check(not any("--series-dir" in c or "--last-days" in c for c in calls)
              and not os.listdir(str(base / "archive/dems")) and not bad,
              f"no --series-dir / --last-days; archive/dems stays empty; every write under archive/weekly or logs "
              f"({len(changed)} files; outside: {bad[:3]})")
        info = json.loads(Path(stem + "_info.json").read_text())
        W = wk(base, m)
        check(info["window_start"] == m and info["window_end"] == "2026-11-08" and m <= info["first_date"]
              and info["last_date"] <= "2026-11-08" and info["days"] == 5,
              f"the DEM's own window and dates lie in Monday..Sunday; 5 usable days (4-5 Nov rough) "
              f"({info['first_date']}..{info['last_date']}, {info['days']} days)")
        check(W["status"] == "final" and W["cells"] >= 100 and W["dem_rc"] == 0 and W["usable_days"] == 5
              and W["points"] == 5 * 25 * 2 * 20 and W["lines"] == 5 * 25 * 2,
              f"week.json: final, {W['cells']} cells, 5 usable days; 5000 waterline points (the rough days' frames "
              "left out, as the DEM leaves them out), 250 lines")
        spec = wf.load_spec(base)
        sizes = {wf.png_size(p) for p in all_pngs(base)}
        hashes = {wf.png_text(p).get("weekly_frames_spec") for p in all_pngs(base)}
        check(sizes == {(1920, 1080)} and hashes == {spec["hash"]} and len(all_pngs(base)) == 6,
              f"every frame 1920x1080 with the spec's hash {spec['hash']} in its PNG metadata ({len(all_pngs(base))} "
              "frames: 2 weeks, 2 final copies)")
        check(spec["extent"] == [419980.0, 420080.0, 4637980.0, 4638050.0] and spec["cell"] == 2.0
              and spec["made_from"]["rows"] == 13 * 25 * 2 * 20 and spec["extent_source"] == "auto"
              and spec["z_limits"][0] % 0.25 == 0 and spec["z_limits"][1] % 0.25 == 0
              and spec["z_limits"][0] < -0.5 and spec["z_limits"][1] > 1.75,
              f"the frozen spec: extent padded 20 m and rounded to 10 m, z limits on 0.25 m ({spec['z_limits']}), "
              "cell = the cron's DEM_CELL")
        check(not (base / "archive/weekly/work").exists(), "the working copy is removed after the run")
    print(f"  ({time.time() - t0:.1f} s)")


# --------------------------------------------------------------------------------------------
# T3 statuses over days
# --------------------------------------------------------------------------------------------

def t3():
    print("T3 current -> provisional -> final (coverage, max wait), the empty week, a failed build, once per day")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        calls = []
        runner = fake_runner(calls)
        write_ground(base, dates("2026-11-02", "2026-11-10"))
        before = tree(base)
        rc, out = run(base, ["daily", "--have-lock"], "2026-11-11T17:00:00Z", runner)
        st = json.loads((base / "archive/weekly/state.json").read_text())
        check(rc == 0 and st["weeks"] == {"2026-11-09": "current"} and st["last_daily"] == "2026-11-11"
              and (base / "archive/weekly/current/dem_frame.png").exists()
              and not (base / "archive/weekly/frames/dem").exists(),
              "Wed 11 Nov, first run: only the current week (9-15 Nov), in current/; no final frame")
        cur1 = (base / "archive/weekly/current/dem_frame.png").read_bytes()
        # a second run the same day: nothing
        t_before = tree(base)
        n = len(calls)
        rc, out = run(base, ["daily", "--have-lock"], "2026-11-11T23:55:00Z", runner)
        check(rc == 0 and "done today already" in out and tree(base) == t_before and len(calls) == n,
              "a second run the same UTC date: 'done today already', no build, no file touched")
        # Monday 16 Nov: the past week provisional at once, a brand new current week
        write_ground(base, dates("2026-11-02", "2026-11-14"))
        rc, out = run(base, ["daily", "--have-lock"], "2026-11-16T17:00:00Z", runner)
        W = wk(base, "2026-11-09")
        cur2 = (base / "archive/weekly/current/dem_frame.png").read_bytes()
        check(W["status"] == "provisional" and W["data_through"] == "2026-11-14"
              and wf.frame_texts(wf.load_spec(base), W)[3] == "PROVISIONAL: data through 2026-11-14"
              and wk(base, "2026-11-16")["status"] == "current" and cur2 != cur1
              and cur2 == (base / "archive/weekly/weeks/2026-11-16/dem_frame.png").read_bytes()
              and not (base / "archive/weekly/frames/dem").exists(),
              "Mon 16 Nov: 9-15 Nov PROVISIONAL (stamp 'PROVISIONAL: data through 2026-11-14'), a brand new "
              "current week 16-22 Nov in current/; nothing final yet")
        W16 = wk(base, "2026-11-16")
        check(W16["cells"] == 0 and W16["points"] == 0 and "no waterline points" in W16["why"]
              and "EMPTY" in wf.frame_texts(wf.load_spec(base), W16)[3],
              "the new week has no data yet (GNSS-R lag): its frames say EMPTY and why")
        # Tue 17 Nov: the file covers 16 Nov -> 9-15 Nov final
        write_ground(base, dates("2026-11-02", "2026-11-16"))
        n = len(calls)
        rc, out = run(base, ["daily", "--have-lock"], "2026-11-17T17:00:00Z", runner)
        W = wk(base, "2026-11-09")
        fd = base / "archive/weekly/frames/dem/2026-11-15_week.png"
        fw = base / "archive/weekly/frames/waterlines/2026-11-15_week.png"
        check(W["status"] == "final" and W["final_reason"] == "covered" and W.get("finalised_utc")
              and fd.exists() and fw.exists()
              and fd.read_bytes() == (base / "archive/weekly/weeks/2026-11-09/dem_frame.png").read_bytes()
              and wf.frame_texts(wf.load_spec(base), W)[3] == "",
              "Tue 17 Nov: the ground file covers 16 Nov: 9-15 Nov FINAL, copied to frames/{dem,waterlines}/"
              "2026-11-15_week.png, no stamp")
        # max wait: the ground file stops at 20 Nov
        write_ground(base, dates("2026-11-02", "2026-11-20"))
        for day in range(18, 30):
            rc, out = run(base, ["daily", "--have-lock"], f"2026-11-{day:02d}T17:00:00Z", runner)
            s = wk(base, "2026-11-16")["status"]
            if day == 28:
                check(s == "provisional", "Sat 28 Nov (6 days after Sun 22 Nov, ground through 20 Nov): 16-22 Nov "
                                          "still provisional")
        W = wk(base, "2026-11-16")
        check(W["status"] == "final" and W["final_reason"] == "waited" and W["data_through"] == "2026-11-20"
              and wf.frame_texts(wf.load_spec(base), W)[3] == "data through 2026-11-20"
              and (base / "archive/weekly/frames/dem/2026-11-22_week.png").exists(),
              "Sun 29 Nov (WEEKLY_MAX_WAIT_DAYS=7 after Sun 22 Nov): 16-22 Nov FINAL, stamped 'data through "
              "2026-11-20'")
        # the week 23-29 Nov has no data at all -> empty when final (max wait)
        for day in (30, 31, 32, 33, 34, 35, 36, 37):
            now = (datetime(2026, 11, 30, 17, tzinfo=UTC) + timedelta(days=day - 30)).strftime("%Y-%m-%dT%H:%M:%SZ")
            rc, out = run(base, ["daily", "--have-lock"], now, runner)
        W = wk(base, "2026-11-23")
        check(W["status"] == "empty" and W["cells"] == 0 and W["points"] == 0 and "no waterline points" in W["why"]
              and (base / "archive/weekly/frames/dem/2026-11-29_week.png").exists()
              and "EMPTY" in wf.frame_texts(wf.load_spec(base), W)[3],
              "a week with no points: EMPTY when final (frames say why; the timeline stays continuous)")
        built_mondays = sorted({c[c.index("--start-date") + 1] for c in calls})
        check(built_mondays == ["2026-11-09", "2026-11-16"],
              f"builds only for weeks with rows (no build for an empty week): {built_mondays}")
        n9 = sum(1 for c in calls if c[c.index("--start-date") + 1] == "2026-11-09")
        check(n9 == 3, f"9-15 Nov built 3 times (11 Nov, 16 Nov, 17 Nov: its rows changed), never after final ({n9})")
        after = tree(base)
        changed, bad = writes_ok(before, after)
        check(not bad, f"every write under archive/weekly or logs ({len(changed)} files; outside: {bad[:3]})")
        sizes = {wf.png_size(p) for p in all_pngs(base)}
        check(sizes == {(1920, 1080)}, f"every frame the same pixel size: {sizes}")
    # rc 4 (the day test rejected all) and a failed build
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        write_ground(base, dates("2026-11-02", "2026-11-17"))
        calls = []
        rc, out = run(base, ["backfill"], "2026-11-18T03:00:00Z",
                      fake_runner(calls, rc_for={"2026-11-02": 4, "2026-11-09": 1}))
        W2, W9 = wk(base, "2026-11-02"), wk(base, "2026-11-09")
        check(W2["status"] == "empty" and "day test" in W2["why"] and W2["dem_rc"] == 4,
              "exit 4 with points in the week: an EMPTY final week ('... the day test rejected every camera-day')")
        check(rc == 1 and W9["status"] == "provisional" and "WARNING: weekly frames: week 2026-11-09" in out
              and "failed (exit 1)" in out,
              "exit 1: a WARNING, rc 1, and the week is not final although covered (built again next run)")
        calls2 = []
        rc, out = run(base, ["daily", "--have-lock"], "2026-11-18T17:00:00Z", fake_runner(calls2))
        check(wk(base, "2026-11-09")["status"] == "final" and rc == 0
              and [c[c.index("--start-date") + 1] for c in calls2] == ["2026-11-09", "2026-11-16"],
              "the next daily run builds it again: final")
        cells_few = []
        with tempfile.TemporaryDirectory() as tmp2:
            b2 = make_base(tmp2)
            write_ground(b2, dates("2026-11-02", "2026-11-10"))
            run(b2, ["backfill"], "2026-11-11T03:00:00Z", fake_runner(cells_few, cells=(2, 3)))
            W = wk(b2, "2026-11-02")
            check(W["status"] == "empty" and "too few cells" in W["why"], f"6 cells (< 10): an EMPTY week ({W['why']})")


# --------------------------------------------------------------------------------------------
# T4 the frozen look, respec
# --------------------------------------------------------------------------------------------

def t4():
    print("T4 the frozen look: a changed setting is a NOTE; respec redraws every frame")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        write_ground(base, dates("2026-11-02", "2026-11-17"))
        calls = []
        run(base, ["backfill"], "2026-11-18T03:00:00Z", fake_runner(calls))
        run(base, ["daily", "--have-lock"], "2026-11-18T17:00:00Z", fake_runner(calls))
        sp = base / "archive/weekly/frame_spec.json"
        spec_bytes = sp.read_bytes()
        old = wf.load_spec(base)
        (base / "weekly_frames_settings.conf").write_text('WEEKLY_SIZE=1280x720\nWEEKLY_Z_LIMITS="-1,3"\n')
        before = {str(p): p.read_bytes() for p in all_pngs(base)}
        rc, out = run(base, ["daily", "--have-lock", "--force"], "2026-11-18T18:00:00Z", fake_runner(calls))
        notes = [l for l in out.splitlines() if l.startswith("NOTE:") and "respec" in l]
        check(rc == 0 and len(notes) == 2 and any("WEEKLY_SIZE=1280x720" in l for l in notes)
              and any("WEEKLY_Z_LIMITS" in l for l in notes) and sp.read_bytes() == spec_bytes,
              "settings changed after the freeze: one NOTE each naming `respec`; frame_spec.json unchanged")
        sizes = {wf.png_size(p) for p in all_pngs(base)}
        hashes = {wf.png_text(p).get("weekly_frames_spec") for p in all_pngs(base)}
        check(sizes == {(1920, 1080)} and hashes == {old["hash"]},
              "every frame (the current week's redrawn too) keeps the frozen size and look")
        final_before = {k: v for k, v in before.items() if "/frames/" in k}
        check(all(Path(k).read_bytes() == v for k, v in final_before.items()), "final frames untouched")
        n = len(calls)
        rc, out = run(base, ["respec"], "2026-11-18T19:00:00Z", fake_runner(calls))
        new = wf.load_spec(base)
        pngs = all_pngs(base)
        sizes = {wf.png_size(p) for p in pngs}
        hashes = {wf.png_text(p).get("weekly_frames_spec") for p in pngs}
        check(rc == 0 and new["size"] == [1280, 720] and new["z_limits"] == [-1.0, 3.0] and new["hash"] != old["hash"]
              and (base / "archive/weekly/frame_spec.prev.json").exists() and "re-frozen" in out,
              f"respec: the look re-frozen from the settings ({old['hash']} -> {new['hash']}), the old one kept")
        check(sizes == {(1280, 720)} and hashes == {new["hash"]} and len(pngs) == len(before),
              f"respec redrew every frame ({len(pngs)}: weeks, current, final copies) at 1280x720 with the new hash")
        check(n == len(calls) and "3 week(s) redrawn" in out,
              "respec rebuilt no DEM (3 weeks redrawn from the kept grids)")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        (base / "contour_points_ground.csv").write_text(GHEAD + "\n")
        rc, out = run(base, ["daily", "--have-lock"], "2026-11-18T17:00:00Z", fake_runner([]))
        st = wf.read_json(base / "archive/weekly/state.json") or {}
        check(rc == 0 and "no frame spec yet" in out and wf.load_spec(base) is None and not st.get("last_daily"),
              "an empty ground file: a NOTE (no frame spec yet), nothing drawn, no stamp (tried again next run)")
        notes = []
        check(wf.ensure_spec(base, base / "contour_points_ground.csv", log=notes.append) is None and notes,
              "ensure_spec on it: None and a NOTE")
        (base / "weekly_frames_settings.conf").write_text('WEEKLY_EXTENT="419000,421000,4637000,4639000"\n'
                                                          'WEEKLY_Z_LIMITS="-1,2"\n')
        spec = wf.ensure_spec(base, base / "contour_points_ground.csv", log=notes.append)
        check(spec and spec["extent"] == [419000.0, 421000.0, 4637000.0, 4639000.0] and spec["z_limits"] == [-1.0, 2.0]
              and spec["extent_source"] == "WEEKLY_EXTENT" and wf.load_spec(base) == spec,
              "the settings' extent and limits make a spec without points (written to frame_spec.json)")


# --------------------------------------------------------------------------------------------
# T5 backfill
# --------------------------------------------------------------------------------------------

def t5():
    print("T5 backfill: dry run, the blind week skipped, resumable")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        write_ground(base, dates("2024-10-14", "2024-11-04"))
        AUDIT["hits"] = []
        AUDIT["on"] = True
        try:
            before = tree(base)
            rc, out = run(base, ["backfill", "--dry-run"], "2024-11-20T03:00:00Z", fake_runner([]))
            check(rc == 0 and tree(base) == before and "would build week 2024-10-14..2024-10-20 (final, covered" in out
                  and "would build week 2024-10-28..2024-11-03" in out and "skipped: it contains the blind date "
                  "2024-10-23" in out, "--dry-run: says what it would build and skip; writes nothing")
            calls, seen_rows = [], []

            def runner(argv, log_path=None):
                if Path(argv[1]).name == "dem_from_contours.py":
                    with open(argv[2]) as f:
                        seen_rows.extend(r["capture_time_utc"][:10] for r in csv.DictReader(f))
                return fake_runner(calls, raise_for={"2024-10-28"})(argv, log_path)
            rc, out = run(base, ["backfill"], "2024-11-20T03:00:00Z", runner)
            check(rc == 1 and "stopped (RuntimeError" in out and wk(base, "2024-10-14")["status"] == "final",
                  "a stop during the 3rd week (power lost): rc 1 with a WARNING; the first week is final already")
            calls2 = []
            rc, out = run(base, ["backfill"], "2024-11-20T03:00:00Z", fake_runner(calls2))
            mondays = [c[c.index("--start-date") + 1] for c in calls2]
            check(rc == 0 and mondays == ["2024-10-28", "2024-11-04"] and "2024-10-14..2024-10-20 already final: skipped" in out,
                  f"resumed: only the unfinished week is built, the final one skipped ({mondays})")
            check(BLIND not in seen_rows and "2024-10-22" in seen_rows and "2024-10-24" in seen_rows,
                  "the working copy of the ground rows never holds a row of the blind date")
            names = [k for k in tree(base)]
            check(not (base / "archive/weekly/weeks/2024-10-21").exists()
                  and not any(wf.blind_named(n) for n in names)
                  and not (base / "archive/weekly/frames/dem/2024-10-27_week.png").exists(),
                  "the blind week 21-27 Oct 2024: no folder, no frame, no file named for the date")
            rc, out = run(base, ["daily", "--have-lock"], "2024-10-23T17:00:00Z", fake_runner(calls2))
            check("skipped: it contains the blind date 2024-10-23" in out
                  and not (base / "archive/weekly/weeks/2024-10-21").exists(),
                  "daily skips the blind week too (were it ever the current week)")
            frames, notes = wf.collect_frames(base, "dem")
            check([n for n, _ in frames] == ["2024-10-20_week.png", "2024-11-03_week.png", "2024-11-10_week.png"],
                  "the DEM video's frames: the two weeks around the blind week and the next one, no gap filler "
                  "for the blind week")
        finally:
            AUDIT["on"] = False
        check(not AUDIT["hits"], f"the audit hook saw no open of a path named for the blind date ({AUDIT['hits'][:2]})")
        check(not [h for h in calls + calls2 if "2024-10-21" in h], "no DEM build of the blind week")


# --------------------------------------------------------------------------------------------
# T6 the videos
# --------------------------------------------------------------------------------------------

FAKE_FFMPEG = """#!/bin/sh
printf '%s\\n' "$@" > "{log}.$$"
cat "{log}.$$" >> "{log}"
printf '%s\\n' '--' >> "{log}"
rm -f "{log}.$$"
for a in "$@"; do last="$a"; done
printf 'fake mp4' > "$last"
"""


def storm_event(base, spec, ev="2026-11-09_fc", entries=None):
    d = base / "archive/storm_watch/events" / ev / "frames"
    d.mkdir(parents=True, exist_ok=True)
    entries = entries if entries is not None else [
        {"part": "pre", "series": "dem", "date": "2026-11-08", "window": ["2026-11-02", "2026-11-08"],
         "png": "dem_pre.png"},
        {"part": "pre", "series": "waterlines", "date": "2026-11-08", "window": ["2026-11-02", "2026-11-08"],
         "png": "waterlines_pre.png"},
        {"part": "post3", "series": "dem", "date": "2026-11-14", "window": ["2026-11-12", "2026-11-14"],
         "png": "dem_post3.png"},
        {"part": "post3", "series": "change", "date": "2026-11-14", "window": ["2026-11-12", "2026-11-14"],
         "png": "change_post3.png"},
        {"part": "post7", "series": "dem", "date": "2026-11-18", "window": ["2026-11-12", "2026-11-18"],
         "png": "dem_post7.png"},
        {"part": "x", "series": "dem", "date": "2026-11-18", "window": [], "png": "../../../weekly/state.png"},
        {"part": "y", "series": "dem", "date": "2026-11-18", "window": [], "png": "missing.png"},
    ]
    for e in entries:
        if "/" not in e["png"] and e["png"] != "missing.png":
            wf.draw_dem_frame(None, d / e["png"], spec, f"storm {e['part']}", banner=f"STORM {ev}: {e['part']}")
    (d / "frames.json").write_text(json.dumps(entries))
    return d


def t6():
    print("T6 the videos: order, ffmpeg's arguments, rebuilt on change only, ffmpeg missing, a real MP4")
    old_path = os.environ.get("PATH", "")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        write_ground(base, dates("2026-11-02", "2026-11-17"))
        fakebin = Path(tmp) / "fakebin"
        fakebin.mkdir()
        flog = Path(tmp) / "ffmpeg_calls.txt"
        (fakebin / "ffmpeg").write_text(FAKE_FFMPEG.format(log=flog))
        (fakebin / "ffmpeg").chmod(0o755)
        os.environ["PATH"] = f"{fakebin}{os.pathsep}{old_path}"
        try:
            calls = []
            rc, out = run(base, ["backfill"], "2026-11-18T03:00:00Z", fake_runner(calls))
            spec = wf.load_spec(base)
            storm_event(base, spec)
            flog.unlink()
            rc, out = run(base, ["daily", "--have-lock"], "2026-11-18T17:00:00Z", fake_runner(calls))
            lst = (base / "archive/weekly/video/dem_frames.txt").read_text()
            names = [l[2:] for l in lst.splitlines() if l.startswith("# ") and l.endswith(".png")]
            expect = ["2026-11-08_storm_2026-11-09_fc_pre.png", "2026-11-08_week.png",
                      "2026-11-14_storm_2026-11-09_fc_post3.png", "2026-11-14_storm_2026-11-09_fc_change_post3.png",
                      "2026-11-15_week.png", "2026-11-18_storm_2026-11-09_fc_post7.png"]
            check(names == expect, f"DEM video order: by date, then name; the change frame after that date's DEM "
                                   f"frame; a bad or missing png left out ({names})")
            files = [l for l in lst.splitlines() if l.startswith("file ")]
            durs = [l for l in lst.splitlines() if l.startswith("duration ")]
            check(lst.startswith("ffconcat version 1.0") and len(files) == 6 and set(durs) == {"duration 0.5"}
                  and files[1] == f"file '{(base / 'archive/weekly/frames/dem/2026-11-08_week.png').resolve()}'"
                  and "events/2026-11-09_fc/frames/dem_pre.png'" in files[0],
                  "the concat list: absolute paths, 1/fps = 0.5 s each; storm frames read in place")
            wl = (base / "archive/weekly/video/waterlines_frames.txt").read_text()
            wnames = [l[2:] for l in wl.splitlines() if l.startswith("# ") and l.endswith(".png")]
            check(wnames == ["2026-11-08_storm_2026-11-09_fc_pre.png", "2026-11-08_week.png", "2026-11-15_week.png"],
                  f"waterline video: its weeks and the storm's waterline frame ({wnames})")
            check("left out" in out and "missing.png" in out, "the bad entries are NOTEs")
            args = flog.read_text().split("--\n")
            args = [a.split("\n")[:-1] for a in args if a.strip()]
            a = args[0] if args else []
            check(len(args) == 2 and a[a.index("-f") + 1] == "concat" and a[a.index("-safe") + 1] == "0"
                  and a[a.index("-i") + 1] == str(base / "archive/weekly/video/dem_frames.txt")
                  and a[a.index("-r") + 1] == "2" and a[a.index("-c:v") + 1] == "libx264"
                  and a[a.index("-pix_fmt") + 1] == "yuv420p" and "scale=1920:1080" in a[a.index("-vf") + 1]
                  and "setpts=N/(2*TB)" in a[a.index("-vf") + 1] and a[-1].endswith(".dem_weekly.mp4.tmp"),
                  "ffmpeg: concat demuxer, -safe 0, one frame each 1/fps (setpts, -r 2), H.264 yuv420p, the frame "
                  "size, a temp file")
            mp4 = base / "archive/weekly/video/dem_weekly.mp4"
            check(mp4.read_bytes() == b"fake mp4" and not list((base / "archive/weekly/video").glob(".*.tmp")),
                  "the temp file replaced the video")
            st = json.loads((base / "archive/weekly/state.json").read_text())
            check(st["videos"]["dem"]["frames"] == 6 and st["videos"]["dem"]["storm_frames"] == 4,
                  "state.json: the last video build (6 frames, 4 storm)")
            n = len(args)
            run(base, ["daily", "--have-lock"], "2026-11-19T17:00:00Z", fake_runner(calls))
            args2 = [x for x in flog.read_text().split("--\n") if x.strip()]
            check(len(args2) == n, "the next day, no frame changed: no video rebuilt")
            p = base / "archive/storm_watch/events/2026-11-09_fc/frames/dem_post7.png"
            os.utime(str(p), (time.time() + 5, time.time() + 5))
            run(base, ["daily", "--have-lock"], "2026-11-20T17:00:00Z", fake_runner(calls))
            args3 = [x for x in flog.read_text().split("--\n") if x.strip()]
            check(len(args3) == n + 1, "a storm frame redrawn (its mtime): the DEM video alone is rebuilt")
            # ffmpeg missing
            empty = Path(tmp) / "nobin"
            empty.mkdir()
            os.environ["PATH"] = str(empty)
            (base / "archive/weekly/frames/dem/2026-11-15_week.png").touch()
            rc, out = run(base, ["daily", "--have-lock"], "2026-11-21T17:00:00Z", fake_runner(calls))
            rc2, out2 = run(base, ["daily", "--have-lock", "--force"], "2026-11-21T18:00:00Z", fake_runner(calls))
            check(rc == 0 and rc2 == 0 and out.count("ffmpeg not found") == 1 and "sudo apt-get install -y ffmpeg" in out
                  and "ffmpeg not found" not in out2 and "WARNING" not in out + out2
                  and (base / "archive/weekly/current/dem_frame.png").exists(),
                  "ffmpeg missing: one NOTE a day (with the install command), no WARNING, rc 0, frames still made")
        finally:
            os.environ["PATH"] = old_path
    if not (REAL_FFMPEG and REAL_FFPROBE):
        print("  (no ffmpeg/ffprobe here: the real MP4 check is skipped)")
        return
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, settings="WEEKLY_SIZE=320x240\nWEEKLY_FPS=4\n")
        write_ground(base, dates("2026-11-02", "2026-11-24"))
        rc, out = run(base, ["backfill"], "2026-11-25T03:00:00Z", fake_runner([]))
        spec = wf.load_spec(base)
        storm_event(base, spec, entries=[{"part": "pre", "series": "dem", "date": "2026-11-08",
                                          "window": ["2026-11-02", "2026-11-08"], "png": "dem_pre.png"}])
        rc, out = run(base, ["daily", "--have-lock"], "2026-11-25T17:00:00Z", fake_runner([]))
        mp4 = base / "archive/weekly/video/dem_weekly.mp4"
        p = subprocess.run([REAL_FFPROBE, "-v", "error", "-count_frames", "-select_streams", "v:0", "-show_entries",
                            "stream=nb_read_frames,width,height,pix_fmt,codec_name", "-of", "json", str(mp4)],
                           stdout=subprocess.PIPE, universal_newlines=True)
        s = (json.loads(p.stdout or "{}").get("streams") or [{}])[0]
        check(rc == 0 and int(s.get("nb_read_frames") or 0) == 4 and s.get("width") == 320 and s.get("height") == 240
              and s.get("pix_fmt") == "yuv420p" and s.get("codec_name") == "h264",
              f"the real ffmpeg: a 320x240 H.264 yuv420p MP4 of 4 frames (3 weeks + 1 storm frame): {s}")


# --------------------------------------------------------------------------------------------
# T7 locks and the CLI
# --------------------------------------------------------------------------------------------

def hold(path):
    p = subprocess.Popen(["flock", str(path), "sleep", "60"], start_new_session=True)
    time.sleep(0.5)
    return p


def t7():
    print("T7 locks: the cron's lock, a cron run due, another weekly run; the CLI")
    saved = (wf.LOCK_WAIT_S, wf.local_now)
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        write_ground(base, dates("2026-11-02", "2026-11-10"))
        wf.LOCK_WAIT_S = 0.3
        try:
            h = hold(base / ".timex_cron.lock")
            try:
                rc, out = run(base, ["daily"], "2026-11-11T17:00:00Z", fake_runner([]))
                st = wf.read_json(base / "archive/weekly/state.json") or {}
                check(rc == 0 and "lock (.timex_cron.lock) is held" in out and not st.get("last_daily"),
                      "the cron's lock held (a manual run): a NOTE, nothing done, no stamp")
                rc, out = run(base, ["backfill"], "2026-11-11T17:00:00Z", fake_runner([]))
                check(rc == 1 and "is held" in out, "backfill with the cron's lock held: an error, rc 1")
                rc, out = run(base, ["daily", "--have-lock"], "2026-11-11T17:00:00Z", fake_runner([]))
                st = wf.read_json(base / "archive/weekly/state.json") or {}
                check(rc == 0 and st.get("last_daily") == "2026-11-11",
                      "--have-lock (the cron holds it): the ground file read directly")
            finally:
                h.kill()
            wf.local_now = lambda: datetime(2026, 11, 12, 12, 20)
            rc, out = run(base, ["daily", "--force"], "2026-11-12T17:00:00Z", fake_runner([]))
            check(rc == 0 and "cron run is due (12:30 waterline_timex_cron.sh local)" in out,
                  "a waterline cron run due within 15 min: the ground file is not read (a NOTE)")
            wf.local_now = lambda: datetime(2026, 11, 12, 3, 0)
            h = hold(base / "logs/weekly_frames.lock")
            try:
                rc, out = run(base, ["daily", "--have-lock"], "2026-11-12T17:00:00Z", fake_runner([]))
                check(rc == 0 and "another weekly_frames.py run is busy: skipped" in out,
                      "another weekly_frames.py run (logs/weekly_frames.lock): skipped")
            finally:
                h.kill()
        finally:
            wf.LOCK_WAIT_S, wf.local_now = saved
        (base / "weekly_frames_settings.conf").write_text("WEEKLY_ENABLE=0\n")
        rc, out = run(base, ["daily", "--have-lock"], "2026-11-13T17:00:00Z", fake_runner([]))
        check(rc == 0 and "off (WEEKLY_ENABLE=0)" in out, "WEEKLY_ENABLE=0: does nothing")
        (base / "weekly_frames_settings.conf").unlink()
        # the CLI, end to end in its own process
        env = dict(os.environ, MPLBACKEND="Agg")
        before = tree(base)
        p = subprocess.run([sys.executable, str(HERE / "weekly_frames.py"), "status", "--base", str(base)],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, env=env)
        check(p.returncode == 0 and "frame spec" in p.stdout and "current     : 1 week(s): 2026-11-09" in p.stdout
              and tree(base) == before, "status: the spec, the weeks by status; read only (nothing written)")
        p = subprocess.run([sys.executable, str(HERE / "weekly_frames.py"), "daily", "--base", str(base),
                            "--have-lock", "--now", "2026-11-11T20:00:00Z"],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, env=env)
        check(p.returncode == 0 and p.stdout.strip() == "weekly frames: done today already (2026-11-11 UTC)",
              f"the cron's call when there is nothing to do: one line, rc 0 ({p.stdout.strip()[:80]})")
        p = subprocess.run([sys.executable, str(HERE / "weekly_frames.py"), "backfill", "--base", str(base),
                            "--from", "2026-13-01"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           universal_newlines=True, env=env)
        check(p.returncode == 2 and "not a date" in p.stdout, "a bad --from: refused (argparse)")


def main():
    os.environ["MPLCONFIGDIR"] = tempfile.mkdtemp(prefix="mpl_wf_")
    os.environ.setdefault("MPLBACKEND", "Agg")
    wf.local_now = lambda: datetime(2026, 1, 1, 3, 0)       # far from every station job
    t_all = time.time()
    for t in (t1, t2, t3, t4, t5, t6, t7):
        t()
    shutil.rmtree(os.environ["MPLCONFIGDIR"], ignore_errors=True)
    print(f"\n{str(len(FAILED)) + ' check(s) FAILED' if FAILED else 'all weekly frames tests passed'} "
          f"({time.time() - t_all:.0f} s)")
    for f in FAILED:
        print("  - " + f)
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()

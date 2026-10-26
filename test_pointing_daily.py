#!/usr/bin/env python3
"""
Tests For The Per-Day Pointing Log And Its Opt-In Correction (pointing_daily.py)
===================================================================================
Plain python3 (no pytest needed on the station; no scipy either):

    python3 test_pointing_daily.py [--keep DIR]

WHY. A day whose camera turned is left out of a survey-date product
(survey_products.py) and is gridded with the wrong pointing in the live
one. pointing_daily.py measures the turn on matched land features and,
when asked (survey_products.py --correct-pointing, or never, live),
maps that day's waterlines with the calibration turned by it. These
checks plant the answer and ask for it back:
  * the gates: every way a day's match can fail one (unmatched,
    ambiguous, under 40 inliers, within its own noise, bunched, a turn
    over 2.5 deg, UNCONFIRMED, the reference day, an older row without
    the numbers), each said; a day that passes is corrected by exactly
    what was measured; the thresholds are survey_products.py's own;
    the survey products' gates of the day: a day at rest (a turn within
    MATCH_ERR_DEG) is not corrected, nor one whose pointing did not
    hold through the day (a turn in the middle of the day, its first
    or last photo unread, or not checked); the pointing step matches a
    day's first and last photos only with --correct-pointing, and a
    camera that turned at 15:00 is not corrected (real matches);
  * the corrected pointing is exact: the CIRN angles plus the measured
    change give the rotation pointing_check.measure fitted, also after
    the EO file is written and read back;
  * synthetic rotation recovery end to end (pointing_daily.self_test):
    planted turns of 0.05-2 deg on a camera that only turned, measured
    by pointing_check.measure to 0.02 deg, and waterline points that sat
    metres off (tens at range) put back to under 0.5 m; an unturned frame
    corrects nothing; dirt on the window that wins the match, or a scene
    too faint to match, corrects nothing either;
  * the waterlines mapped again: the corrected days' rows equal what
    georectify.py gives with the day's EO file (string for string),
    every other row is unchanged byte for byte, and a file with no
    corrected day is not touched;
  * survey_products.py: the default build is as before (the same days
    left out, no correction in the stamp or the signature); with
    --correct-pointing a passing day is detected and mapped with its own
    pointing (the live era: the station's rows; the historical eras: the
    detection with that day no longer skipped), a failing one is still
    left out, and the README, provenance and the rebuild line say so;
    the root rule (its own <root>_pointing-corrected, never a blind date
    or --eo); the compare step's C = 0 lines of a corrected day are
    re-projected with that day's EO (survey_compare.py --camera-eo-day),
    the default compare command as before;
  * the live log: complete days not yet logged, the reference made from
    --reference-day, a turned day, a day within its noise, a turn over
    the limit, the reference frame itself, another station's photos,
    a rerun that measures nothing again (and says a new --reference-day
    is not used), no reference at all; each day records its calibration
    and reference; nothing is applied, nothing written where another job
    reads;
  * the cron calls it log-only, with the calibrations of its own step 6.

STATION FILES. Only the repository's calibration files are read; every
photo, log and product is made in a temporary folder.
"""

import os
import io
import re
import sys
import csv
import json
import shutil
import argparse
import calendar
import tempfile
import contextlib
import subprocess
from pathlib import Path
from types import SimpleNamespace
from datetime import datetime, timezone, timedelta

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pointing_daily as pdl                                 # noqa: E402
import survey_products as sp                                 # noqa: E402
from georectify import load_intrinsics, load_extrinsics, cirn_angles_to_R, pixel_to_ground   # noqa: E402

CAL = HERE / "calibration"
IO2 = load_intrinsics(CAL / "CACO05_c2_20240801_IO.yaml")
EO2 = load_extrinsics(CAL / "CACO05_c2_20251113_EO-CV.yaml")


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


def prow(date, verdict="matches", features="ok", inliers=400, sigma=12.0, spread="250/520", d=(0.0, 0.0, 0.0),
         day_check=None):
    """A pointing row as survey_products' pointing step writes it (the columns the gates read); day_check: the
    records of the day's first and last photos (--correct-pointing; None: not checked)."""
    r = {"camera": "c2", "date": date, "verdict": verdict, "features": features, "inliers": inliers,
         "turn_sigma": sigma, "inlier_spread_px": spread, "d_azimuth_deg": "", "d_tilt_deg": "", "d_roll_deg": ""}
    if features in ("ok", "moved", "MOVED", "ambiguous"):
        r.update(d_azimuth_deg=d[0], d_tilt_deg=d[1], d_roll_deg=d[2])
    if day_check is not None:
        r.update(feature_time_utc="17:00", day_check=day_check)
    return r


def frec(time_utc, d, features="MOVED", inliers=60, sigma=9.0, spread="230/480"):
    """One of a day's other frames as check_pointing records it (pointing_day_check)."""
    return dict(prow("x", None, features, inliers, sigma, spread, d), feature_frame=f"f{time_utc}.jpg",
                time_utc=time_utc)


STORM_D = (-1.214, 0.004, 0.03)


# ---------------------------------------------------------------------

def test_gates():
    print("the gates: a day is corrected only on a match the survey products trust, and more")
    rules = sp.pointing_correction_rules()
    frame_rules = pdl.station_rules()
    check(frame_rules == {"min_inliers": sp.MATCH_INLIERS, "min_turn_sigma": sp.TURN_SIGMA_MIN,
                          "min_spread_px": float(sp.MATCH_SPREAD_TESTED_PX), "max_turn_deg": pdl.CORRECT_MAX_DEG}
          and (rules["min_inliers"], rules["min_turn_sigma"], rules["min_spread_px"], rules["max_turn_deg"])
          == (40, 4.0, 180.0, 2.5)
          and rules == dict(frame_rules, min_turn_deg=sp.MATCH_ERR_DEG, max_day_spread_deg=sp.MATCH_ERR_DEG)
          and pdl.frame_gates(rules) == frame_rules and sp.MATCH_ERR_DEG == 0.1,
          "the thresholds are survey_products.py's own (40 inliers, turn_sigma 4, 180 px) and CORRECT_MAX_DEG 2.5; "
          "the survey products add the gates of the day at MATCH_ERR_DEG 0.1, the live log's frames have none")
    held = [frec("14:05", (-1.20, 0.01, 0.02)), frec("19:55", (-1.23, 0.0, 0.05))]
    storm = prow("2026-09-27", "DIFFERENT POINTING", "MOVED", 48, 9.1, "210/450", STORM_D, day_check=held)
    d, why = pdl.day_correction(storm, rules, "2026-09-29")
    check(d == [-1.214, 0.004, 0.03] and why == [],
          "c1's measured -1.21 deg storm pan (48 inliers, turn_sigma 9.1, spread 210 px), the same on the day's "
          "14:05 and 19:55 frames: corrected by exactly the probe's turn")
    rest = prow("2026-09-30", sigma=12.0, d=(0.012, -0.004, 0.0), day_check=[])
    d, why = pdl.day_correction(rest, rules, "2026-09-29")
    check(d is None and why[0].startswith("a turn of 0.012 deg, within the match's own error") and "at rest" in why[0],
          f"a 'matches' day whose small turn is 12 of its standard errors: at rest, NOT corrected ({why[0][:70]})")
    check(pdl.day_correction(rest, frame_rules, "2026-09-29")[0] == [0.012, -0.004, 0.0],
          "... while the live log's frame gates still read it (its record of the drift)")
    small = prow("2026-09-30", "small change", "moved", 300, 15.0, d=(0.0, 0.18, 0.0), day_check=[])
    check(pdl.day_correction(small, rules, "2026-09-29")[0] == [0.0, 0.18, 0.0],
          "a 'small change' kept by default (0.18 deg of tilt), its only photo: corrected")
    print("  a turn in the middle of the day, and a day not shown to hold:")
    mid = {
        "came back by 14:05 (at rest, within its noise)":
            ([frec("14:05", (0.0, 0.001, 0.0), "ok", 300, 1.6), held[1]], "the 14:05 UTC frame gives no reading"),
        "14:05 read at rest on a confident match":
            ([frec("14:05", (0.02, 0.01, 0.0), "ok", 300, 8.0), held[1]], "the 14:05 UTC frame reads +0.020"),
        "turned further by 19:55":
            ([held[0], frec("19:55", (-1.40, 0.0, 0.03))], "the 19:55 UTC frame reads -1.400"),
        "fog at 14:05": ([frec("14:05", None, "unmatched", 12), held[1]], "the 14:05 UTC frame gives no reading"),
        "not checked (an older row)": (None, "the pointing through the day not checked"),
    }
    for name, (chk, want) in mid.items():
        d, why = pdl.day_correction(prow("2026-09-27", "DIFFERENT POINTING", "MOVED", 48, 9.1, "210/450", STORM_D,
                                         day_check=chk), rules, "2026-09-29")
        check(d is None and any(w.startswith(want) for w in why), f"{name}: not corrected ({'; '.join(why)[:110]})")
    cases = {
        "unmatched": (prow("x", "UNCONFIRMED", "unmatched", 12), "no feature reading (unmatched)"),
        "ambiguous": (prow("x", "small change", "ambiguous", 300, 20.0, d=(0.0, 0.0, 0.0)),
                      "no feature reading (ambiguous)"),
        "39 inliers": (prow("x", "DIFFERENT POINTING", "MOVED", 39, 9.0, d=(1.0, 0, 0)), "39 inliers, under 40"),
        "noise": (prow("x", sigma=2.7, d=(0.0, 0.001, 0.0)), "turn_sigma 2.7, under 4"),
        "bunched": (prow("x", "DIFFERENT POINTING", "MOVED", 90, 15.0, "150/400", (0.8, 0, 0)), "inliers bunched"),
        "over the limit": (prow("x", "DIFFERENT POINTING", "MOVED", 300, 50.0, d=(3.1, 0.0, 0.0)),
                           "a turn of 3.10 deg, over 2.5"),
        "UNCONFIRMED": (prow("x", "UNCONFIRMED", "ok", 300, 9.0, d=(0.02, 0, 0)), "UNCONFIRMED by the pointing step"),
        "older row": (dict(prow("x", "DIFFERENT POINTING", "MOVED", 300, "", "", (1.0, 0, 0))),
                      "turn_sigma not measured"),
    }
    for name, (row, want) in cases.items():
        d, why = pdl.day_correction(row, rules, "2026-09-29")
        check(d is None and any(w.startswith(want) for w in why), f"{name}: not corrected ({'; '.join(why)})")
    d, why = pdl.day_correction(prow("2026-09-29", sigma=30.0, d=(0.001, 0, 0)), rules, "2026-09-29")
    check(d is None and why[0].startswith("the reference day"), "the reference day itself: never corrected")
    corr, not_ = pdl.window_corrections([prow("2026-09-30", sigma=2.0), storm, prow("2026-09-29")], rules,
                                        "2026-09-29")
    check(list(corr) == ["2026-09-27"] and corr["2026-09-27"]["inliers"] == 48 and sorted(not_) == ["2026-09-29",
                                                                                                    "2026-09-30"],
          "window_corrections: the days that pass, by date, and why each other one does not")
    c27 = corr["2026-09-27"]
    check([f["time_utc"] for f in c27["day_check"]] == ["14:05", "19:55"] and c27["day_spread_deg"] == 0.02
          and c27["probe_time_utc"] == "17:00" and c27["day_check"][0]["d_deg"] == [-1.2, 0.01, 0.02],
          "... with the frames the day was checked on and their largest difference from the probe (0.02 deg)")


def test_corrected_eo_exact():
    print("the corrected pointing is the rotation the match fitted: exact, and after the EO file round trip")
    import pointing_check as pc
    rng = np.random.default_rng(5)
    worst = worst_file = 0.0
    d = Path(tempfile.mkdtemp(dir="."))
    for k in range(200):
        eo = np.r_[420000.0, 4638000.0, 20.0, np.radians([rng.uniform(0, 360), rng.uniform(55, 80),
                                                           rng.uniform(-10, 10)])]
        w = rng.normal(0, 1, 3)
        w *= np.radians(rng.uniform(0.01, 2.5)) / np.linalg.norm(w)
        K = np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])
        th = np.linalg.norm(w)
        Q = np.eye(3) + np.sin(th) / th * K + (1 - np.cos(th)) / th ** 2 * K @ K       # Rodrigues
        R_new = Q @ cirn_angles_to_R(*eo[3:6])
        # what pointing_check.measure reports: R's CIRN angles less the reference's, the pan wrapped
        dd = np.degrees(np.array(pc.R_to_angles(R_new)) - eo[3:6])
        dd[0] = (dd[0] + 180) % 360 - 180
        got = cirn_angles_to_R(*pdl.corrected_eo(eo, dd)[3:6])
        worst = max(worst, float(np.abs(got - R_new).max()))
        back = pdl.write_day_eo(d / "x.yaml", pdl.corrected_eo(eo, dd), ["test"])
        worst_file = max(worst_file, float(np.abs(cirn_angles_to_R(*back[3:6]) - R_new).max()),
                         float(np.abs(back[:3] - eo[:3]).max()))
    check(worst < 1e-12, f"200 random turns of 0.01-2.5 deg: angles + change give the fitted rotation ({worst:.1e})")
    check(worst_file < 1e-8, f"... and so does the EO file as georectify.load_extrinsics reads it back ({worst_file:.1e})")
    shutil.rmtree(d)


def test_self_test():
    print("synthetic rotation recovery end to end (pointing_daily.self_test, a scene of c2's size, no scipy)")
    rc, out = quiet(pdl.self_test)
    lines = [l for l in out.splitlines() if l.strip().startswith("turn ")]
    check(rc == 0 and out.rstrip().endswith("PASS") and len(lines) == 5 and all(l.endswith("-> ok") for l in lines),
          "PASS: " + " | ".join(re.sub(r"\s+", " ", l.strip())[:120] for l in lines[1:3]))
    check("not corrected (turn_sigma" in lines[0] and "under 4" in lines[0],
          "the unturned frame matches within its own noise: nothing corrected, its rows untouched")


def test_dirt_and_fog_gated():
    print("a match the window's dirt wins, or a scene too faint to match, corrects nothing (real matches)")
    import cv2
    import pointing_check as pc
    nu, nv = int(IO2[0]), int(IO2[1])
    eo2 = pdl.corrected_eo(EO2, (1.0, 0.0, 0.0))
    d = Path(tempfile.mkdtemp(dir="."))

    def dirt(img, n):
        r, out = np.random.default_rng(7), img.copy()
        for _ in range(n):
            x, y, rad = int(r.uniform(50, nu - 50)), int(r.uniform(50, nv - 50)), int(r.uniform(3, 8))
            cv2.circle(out, (x, y), rad, (20, 20, 20), -1)
            cv2.circle(out, (x, y), rad + 3, (235, 235, 235), 2)
        return out
    rules = sp.pointing_correction_rules()
    for tag, contrast, n, want in (("dirt wins", 18, 1100, "no feature reading (ambiguous)"),
                                   ("scene lost", 15, 1200, "turn_sigma")):
        tex = pdl.textured_scene(IO2, contrast=contrast, seed=1)
        ref = d / f"ref_{n}.jpg"
        cv2.imwrite(str(ref), dirt(tex, n))
        ref.with_suffix(".json").write_text(json.dumps({"eo": EO2.tolist()}))
        m = pc.measure(dirt(pdl.rotated_view(tex, IO2, EO2, eo2), n), [pc.Reference(ref)], IO2, min_inliers=40,
                       second=True)
        row = sp.feature_record(m, Path("x.jpg"))
        dd, why = pdl.day_correction(row, rules)
        check(dd is None and any(w.startswith(want) for w in why),
              f"{tag}: the camera turned 1.0 deg, the match reads {row.get('features')} at "
              f"{row.get('d_azimuth_deg')} deg ({row.get('inliers')} inliers, turn_sigma {row.get('turn_sigma')}): "
              f"not corrected ({'; '.join(why)})")
    shutil.rmtree(d)


def _contours(path, rows):
    fields = ["source_file", "camera", "capture_time_utc", "capture_epoch", "pixel_column", "pixel_row",
              "tide_elevation_navd88", "offshore_hs_m", "offshore_tp_s"]
    if any("beach_elevation_navd88" in r for r in rows):
        fields.insert(fields.index("tide_elevation_navd88") + 1, "beach_elevation_navd88")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def _lines_for(cam, t, eo, n=40, seed=0, z=0.6, setup=None):
    """Waterline rows of one frame of `cam` at time t, as its detector reports them for a camera at eo (the true
    pointing): pixels of ground points on the plane z (all inside the frame). setup: the rows also carry
    beach_elevation_navd88 = z, their still water (tide_elevation_navd88) that much lower, as a C > 0 build's."""
    from view_reproject import ground_to_pixel
    E, N = pdl.beach_points(IO2, EO2, z=z, n=3 * n, seed=seed)
    U, V, ok = ground_to_pixel(E, N, np.full_like(E, z), IO2, eo)
    E, N, U, V = E[ok][:n], N[ok][:n], U[ok][:n], V[ok][:n]
    e = calendar.timegm(t.timetuple())
    name = photo_name(t, "CACO05", cam)[:-4]
    rows = [{"source_file": name, "camera": cam, "capture_time_utc": f"{t:%Y-%m-%dT%H:%M:%SZ}", "capture_epoch": e,
             "pixel_column": f"{u:.2f}", "pixel_row": f"{v:.2f}", "tide_elevation_navd88": f"{z:.4f}",
             "offshore_hs_m": "1.0", "offshore_tp_s": "8.0"} for u, v in zip(U, V)]
    if setup is not None:
        for r in rows:
            r.update(tide_elevation_navd88=f"{z - setup:.4f}", beach_elevation_navd88=f"{z:.4f}")
    return rows, E, N


def test_regeorectify(d):
    print("the corrected days' rows are what georectify.py gives with the day's EO (on the beach plane of a C > 0 "
          "build, not the still water); every other row is unchanged")
    d = Path(d) / "regeo"
    d.mkdir()
    t1, t2 = datetime(2026, 9, 27, 17, tzinfo=timezone.utc), datetime(2026, 9, 29, 17, tzinfo=timezone.utc)
    turn = (-1.214, 0.004, 0.03)
    eo_t = pdl.corrected_eo(EO2, turn)
    # each line on its beach plane (still water + a 0.25 m setup, as C = 0.037 gives in ~1.5 m waves), as
    # georectify.py places it; mapped on the still water instead it would sit metres off
    r1, E1, N1 = _lines_for("c2", t1, eo_t, seed=1, setup=0.25)
    r2, _, _ = _lines_for("c2", t2, EO2, seed=2, setup=0.25)
    r3 = [dict(r, camera="c1") for r in _lines_for("c2", t1, EO2, seed=3, setup=0.25)[0]]   # another camera
    _contours(d / "contours.csv", r1 + r2 + r3)
    io1 = CAL / "CACO05_c1_20240801_IO.yaml"
    eo1 = CAL / "CACO05_c1_20251113_EO-CV.yaml"
    geo = [sys.executable, str(HERE / "georectify.py"), str(d / "contours.csv"), str(d / "ground.csv"),
           "--io-c1", str(io1), "--eo-c1", str(eo1), "--io-c2", str(CAL / "CACO05_c2_20240801_IO.yaml"),
           "--eo-c2", str(CAL / "CACO05_c2_20251113_EO-CV.yaml")]
    subprocess.run(geo, check=True, stdout=subprocess.DEVNULL)
    before = (d / "ground.csv").read_bytes()
    st = pdl.regeorectify(d / "ground.csv", "c2", IO2, {"2026-09-28": eo_t})
    check(st == {} and (d / "ground.csv").read_bytes() == before, "no row on the corrected days: the file is not touched")
    day_eo = pdl.write_day_eo(d / "day.yaml", eo_t, ["test"])
    st = pdl.regeorectify(d / "ground.csv", "c2", IO2, {"2026-09-27": day_eo})
    after = (d / "ground.csv").read_bytes().split(b"\r\n")
    old = before.split(b"\r\n")
    # the reference: georectify.py itself with the day's EO file, on the same rows
    subprocess.run(geo[:-1] + [str(d / "day.yaml")], check=True, stdout=subprocess.DEVNULL)
    ref = (d / "ground.csv").read_bytes().split(b"\r\n")
    hdr = old[0].decode().split(",")
    ic, iday = hdr.index("camera"), hdr.index("capture_time_utc")
    moved = [k for k in range(1, len(old)) if old[k] and old[k].decode().split(",")[ic] == "c2"
             and old[k].decode().split(",")[iday].startswith("2026-09-27")]
    same = [k for k in range(len(old)) if k not in moved]
    check(len(moved) == len(r1) and all(after[k] == ref[k] for k in moved) and all(after[k] != old[k] for k in moved),
          f"the {len(moved)} c2 rows of 2026-09-27: exactly georectify.py's with the day's EO file")
    check(all(after[k] == old[k] for k in same) and len(after) == len(old),
          f"the other {len(same) - 1} rows (c2 2026-09-29, c1 the same day) and the header: byte for byte as they were")
    back = list(csv.DictReader(open(d / "ground.csv", newline="")))
    # georectify.py rewrote the file above: what regeorectify wrote is `after`; its numbers against the truth
    rows_after = list(csv.DictReader(io.StringIO(b"\r\n".join(after).decode())))
    err = [np.hypot(float(r["easting_utm19"]) - e, float(r["northing_utm19"]) - n)
           for r, e, n in zip([r for r in rows_after if r["camera"] == "c2" and r["capture_time_utc"].startswith(
               "2026-09-27")], E1, N1)]
    check(st["2026-09-27"]["rows"] == len(r1) and st["2026-09-27"]["georectified"] == len(r1)
          and st["2026-09-27"]["moved_m_median"] > 1.0 and max(err) < 0.01 and back,
          f"they now lie where the turned camera saw them (to {max(err) * 1000:.1f} mm), "
          f"{st['2026-09-27']['moved_m_median']:.1f} m (median) from where the calibration put them")


def test_apply_corrections(d):
    print("apply_corrections: one EO file per corrected day (the folder emptied first) and the record")
    d = Path(d) / "apply"
    d.mkdir()
    t1 = datetime(2026, 9, 27, 17, tzinfo=timezone.utc)
    r1, _, _ = _lines_for("c2", t1, pdl.corrected_eo(EO2, (-1.214, 0.004, 0.03)), seed=1)
    _contours(d / "contours.csv", r1)
    subprocess.run([sys.executable, str(HERE / "georectify.py"), str(d / "contours.csv"), str(d / "ground.csv"),
                    "--io-c2", str(CAL / "CACO05_c2_20240801_IO.yaml"),
                    "--eo-c2", str(CAL / "CACO05_c2_20251113_EO-CV.yaml")], check=True, stdout=subprocess.DEVNULL)
    (d / "eo").mkdir()
    (d / "eo" / "stale_EO.yaml").write_text("x: 0\n")
    corr = {"2026-09-27": {"d_deg": [-1.214, 0.004, 0.03], "inliers": 48, "turn_sigma": 9.1,
                           "inlier_spread_px": "210/450", "verdict": "DIFFERENT POINTING"},
            "2026-09-30": {"d_deg": [0.01, 0.0, 0.0], "inliers": 400, "turn_sigma": 12.0,
                           "inlier_spread_px": "250/520", "verdict": "matches"}}
    rec = pdl.apply_corrections(d / "ground.csv", "c2", IO2, EO2, "CACO05_c2_20251113_EO-CV.yaml", corr, d / "eo",
                                reference="ref.jpg")
    files = sorted(p.name for p in (d / "eo").iterdir())
    txt = (d / "eo" / "c2_2026-09-27_pointing_corrected_EO.yaml").read_text()
    check(files == ["c2_2026-09-27_pointing_corrected_EO.yaml", "c2_2026-09-30_pointing_corrected_EO.yaml"]
          and "CACO05_c2_20251113_EO-CV.yaml with the turn" in txt and "-1.214/+0.004/+0.030" in txt
          and "Not a calibration" in txt and " rotated by " not in txt and " turned by " not in txt,
          "one file per day, the stale one gone; its notes say what it is (never read as a turned calibration: no "
          "'X.yaml rotated/turned by N deg' line)")
    check(rec["2026-09-27"]["rows"] == len(r1) and rec["2026-09-30"]["rows"] == 0
          and rec["2026-09-27"]["eo_sha256"] and rec["2026-09-27"]["moved_m_median"] > 1,
          "the record: rows mapped per day (none on a day without rows), the file's sha256, how far they moved")


# --- survey_products.py ---------------------------------------------------------------------------------------

def _res():
    rows = [prow("2026-09-27", "DIFFERENT POINTING", "MOVED", 48, 9.1, "210/450", STORM_D,
                 day_check=[frec("14:05", (-1.20, 0.01, 0.02)), frec("19:55", (-1.23, 0.0, 0.05))]),
            prow("2026-09-28", "UNCONFIRMED", "unmatched", 3),
            prow("2026-09-29", sigma=30.0),
            prow("2026-09-30", sigma=2.0, d=(0.0, 0.001, 0.0))]
    for r in rows:
        r["why"] = "test"
    return {"rows": rows, "info": {"reference_day": "2026-09-29", "reference_frame": "ref.c2.timex.jpg"},
            "different_days": ["2026-09-27", "2026-09-28"], "leave_out": True}


def test_skip_days():
    print("survey_products.skip_days_of: the default build leaves out what it did; --correct-pointing detects the "
          "days it corrects")
    res = _res()
    check(sp.skip_days_of(res, SimpleNamespace(keep_moved_days=False)) == ["2026-09-27", "2026-09-28"]
          and sp.skip_days_of(res, None) == ["2026-09-27", "2026-09-28"],
          "default (and an args without the option): both left out, as before")
    check(sp.skip_days_of(res, SimpleNamespace(keep_moved_days=False, correct_pointing=True)) == ["2026-09-28"],
          "--correct-pointing: 27 Sep (passes the gates) detected; 28 Sep (UNCONFIRMED, unmatched) still left out")
    check(sp.skip_days_of(res, SimpleNamespace(keep_moved_days=True, correct_pointing=True)) == [],
          "--keep-moved-days still keeps every day")
    check(sp.pointing_correction_plan({}, {}, "c2", SimpleNamespace()) is None
          and sp.pointing_correction_plan({}, {}, "c2", SimpleNamespace(correct_pointing=False)) is None,
          "without the option nothing is planned (the detection's signature gets no new key)")


def test_root_rule():
    print("survey_products: a corrected build only in its own root, never a blind date or with --eo")
    plan = {"cams": {"c2": {}}, "blind": None}
    root = "/x/survey_products"
    a = lambda **k: SimpleNamespace(output_root=k.pop("root", root), **k)               # noqa: E731
    check(sp.pointing_root_problem(plan, a(correct_pointing=False)) == "" and
          sp.pointing_root_problem(plan, a()) == "", "a default build in a default root: fine")
    t = sp.pointing_root_problem(plan, a(correct_pointing=True))
    check("--output-root /x/survey_products_pointing-corrected" in t, f"--correct-pointing elsewhere: refused ({t[:90]})")
    check(sp.pointing_root_problem(plan, a(correct_pointing=True, root=root + "_pointing-corrected")) == ""
          and sp.pointing_root_problem(plan, a(correct_pointing=True, root=root + "_pointing-corrected_c0.034")) == "",
          "in <root>_pointing-corrected (and with a C suffix): fine")
    t = sp.pointing_root_problem(plan, a(correct_pointing=False, root=root + "_pointing-corrected"))
    check("named for pointing-corrected products" in t, "a default build in a corrected root: refused")
    t = sp.pointing_root_problem(dict(plan, blind={"x": 1}), a(correct_pointing=True, root=root + "_pointing-corrected"))
    check("BLIND" in t, "a blind date: refused")
    t = sp.pointing_root_problem({"cams": {"c2": {"eo_given": "x.yaml"}}},
                                 a(correct_pointing=True, root=root + "_pointing-corrected"))
    check("--eo" in t, "with --eo: refused")
    check(sp.option_refusal.__doc__ and "pointing_root_problem" in open(HERE / "survey_products.py").read()
          .split("def option_refusal")[1].split("def ")[0], "option_refusal checks it (exit 5, before anything is written)")


def _live_plan(out):
    c = {"eo_path": str(CAL / "CACO05_c2_20251113_EO-CV.yaml"), "io_path": str(CAL / "CACO05_c2_20240801_IO.yaml"),
         "eo_file": "CACO05_c2_20251113_EO-CV.yaml", "first": "2026-09-27", "last": "2026-09-30", "hours": (13.5, 18.0),
         "hours_text": "13.5-18", "station": "CACO05", "variant": None, "notes": ""}
    return {"date": "2026-09-29", "era": "live", "first": "2026-09-27", "last": "2026-09-30", "out": Path(out),
            "cams": {"c2": c}, "surveys": []}


def _day_rows():
    """c2's archive rows for 27-30 Sep: 27 Sep seen by a camera turned by c1's measured storm pan, the others at the
    calibration. -> (rows, {day: (E, N)} true ground points)"""
    rows, truth = [], {}
    for k, day in enumerate(("2026-09-27", "2026-09-28", "2026-09-29", "2026-09-30")):
        t = datetime.strptime(day, "%Y-%m-%d").replace(hour=17, tzinfo=timezone.utc)
        eo = pdl.corrected_eo(EO2, (-1.214, 0.004, 0.03)) if day == "2026-09-27" else EO2
        r, E, N = _lines_for("c2", t, eo, n=30, seed=10 + k)
        rows += r
        truth[day] = (E, N)
    return rows, truth


def test_live_detect(d):
    print("survey_products live era: the default detection as before; --correct-pointing maps 27 Sep with its own "
          "pointing")
    rows, truth = _day_rows()
    src = Path(d) / "live_rows.csv"
    _contours(src, rows)
    saved = (sp.live_rows_file, sp.scan_photos, sp.link_frame_photos)
    sp.live_rows_file = lambda plan, args, state: src
    sp.scan_photos = lambda plan, cam, args, skip_days=(): ([], {"missing_roots": [], "duplicates": 0})
    sp.link_frame_photos = lambda plan, cam, args: ([], 0)
    try:
        got = {}
        for tag, on in (("default", False), ("corrected", True)):
            out = Path(d) / f"live_{tag}"
            plan = _live_plan(out)
            # C = 0: each line's plane is its still water, where the test put the points
            args = SimpleNamespace(setup_coef=0.0, force=False, keep_moved_days=False, correct_pointing=on,
                                   setup_fitted_to=None)
            state = {"pointing": {"c2": _res()}}
            fit = sp.setup_fit_info(args.setup_coef, args)
            parts, log = quiet(sp.live_detect, plan, args, state, fit)
            st = sp.read_stamp(out, "detect_c2")
            got[tag] = (parts, st, list(csv.DictReader(open(out / "waterlines" / "c2" / "contour_points_ground.csv",
                                                             newline=""))), log)
    finally:
        sp.live_rows_file, sp.scan_photos, sp.link_frame_photos = saved
    _, st0, g0, _ = got["default"]
    _, st1, g1, log1 = got["corrected"]
    days0 = sorted({r["capture_time_utc"][:10] for r in g0})
    check(days0 == ["2026-09-29", "2026-09-30"] and st0["skip_days"] == ["2026-09-27", "2026-09-28"]
          and "pointing_correction" not in st0 and "pointing_correction" not in st0["sig"]
          and "pointing_daily.py" not in st0["software"]["scripts"],
          "default: 27 and 28 Sep left out as before; no correction in the stamp, its signature or its scripts")
    days1 = sorted({r["capture_time_utc"][:10] for r in g1})
    pcr = st1.get("pointing_correction") or {}
    check(days1 == ["2026-09-27", "2026-09-29", "2026-09-30"] and st1["skip_days"] == ["2026-09-28"]
          and list(pcr.get("days") or {}) == ["2026-09-27"] and st1["sig"]["pointing_correction"]["days"] ==
          {"2026-09-27": [-1.214, 0.004, 0.03]} and "pointing_daily.py" in st1["software"]["scripts"],
          "--correct-pointing: 27 Sep detected and corrected (stamp, signature, scripts); 28 Sep still left out")
    E, N = truth["2026-09-27"]
    r27 = [r for r in g1 if r["capture_time_utc"].startswith("2026-09-27")]
    err = np.hypot(np.array([float(r["easting_utm19"]) for r in r27]) - E,
                   np.array([float(r["northing_utm19"]) for r in r27]) - N)
    E2, N2 = pixel_to_ground(np.array([float(r["pixel_column"]) for r in r27]),
                             np.array([float(r["pixel_row"]) for r in r27]),
                             np.array([float(r["beach_elevation_navd88"]) for r in r27]), IO2, EO2)
    off = np.hypot(E2 - E, N2 - N)
    check(err.max() < 0.02 and np.median(off) > 1.0,
          f"27 Sep's lines lie where the turned camera saw them (to {err.max() * 1000:.0f} mm), not "
          f"{np.median(off):.1f} m (median) off with the calibration")
    check([r for r in g1 if not r["capture_time_utc"].startswith("2026-09-27")] == g0,
          "every other row as the default build's, value for value")
    eo_f = Path(d) / "live_corrected" / "waterlines" / "c2" / "pointing_corrected_eo" / \
        "c2_2026-09-27_pointing_corrected_EO.yaml"
    check(eo_f.exists() and "pointing correction c2" in log1 and "2026-09-27 -1.21/+0.00/+0.03 deg" in log1,
          "the day's EO file is next to the waterlines, and the step says what it corrected")
    # what the README and provenance say, from the stamp
    byday = {r["date"]: r for r in _res()["rows"]}
    cav_ = sp.pointing_correction_caveat("c2", pcr, byday)
    check(cav_.startswith("c2: POINTING CORRECTED PER DAY (--correct-pointing")
          and "2026-09-27 (DIFFERENT POINTING): -1.214/+0.004/+0.030 deg pan/tilt/roll on 48 inliers" in cav_
          and "(the 17:00 UTC frame; held through the day: its 14:05, 19:55 UTC frames read it within 0.020 deg)"
          in cav_ and "of more than 0.1 deg (the match's own error" in cav_ and "held through the day)" in cav_
          and "Only the georectification is corrected" in cav_
          and "NOT corrected" in cav_ and "2026-09-28 (" in cav_ and "2026-09-30" not in cav_.split("NOT corrected")[1]
          and all(ord(ch) < 128 for ch in cav_),
          "README: the corrected day with its change, the frames it held on and how far its lines moved; the gates of "
          "the day; the day measured to differ that a gate kept from it, with why (not the kept days); the detector's "
          "search is the calibration's; ASCII")
    cav = {}
    for tag in ("default", "corrected"):
        out = Path(d) / f"live_{tag}"
        (out / "pointing").mkdir(exist_ok=True)
        (out / "pointing" / "pointing_c2.json").write_text(json.dumps(_res()))
        cav[tag] = [c for c in sp.collect_caveats(_live_plan(out), cfg_args(correct_pointing=tag == "corrected"), {},
                                                  {}, {"frames_per_day": {}}) if c.startswith("c2")]
    check(cav["default"] == ["c2: days left out of the detection for their pointing: 2026-09-27 (DIFFERENT POINTING: "
                             "test); 2026-09-28 (UNCONFIRMED: test)."]
          and cav["corrected"][0] == "c2: days left out of the detection for their pointing: 2026-09-28 (UNCONFIRMED: "
                                     "test)."
          and len(cav["corrected"]) == 2 and cav["corrected"][1] == cav_
          and not any("KEPT" in c or "stale" in c or "run the detection again" in c for c in cav["corrected"]),
          "KNOWN CAVEATS (collect_caveats): the default product's as before; the corrected one's left-out day, then "
          "the correction, nothing 'KEPT' or stale")
    check("--correct-pointing" in sp.state_argv(SimpleNamespace(**dict(vars(cfg_args()), correct_pointing=True)),
                                                "2026-09-29")
          and "--correct-pointing" not in sp.state_argv(cfg_args(), "2026-09-29"),
          "the 'rebuild everything' line carries --correct-pointing, only when given")


def cfg_args(**kw):
    a = dict(config=str(sp.DEFAULT_CONFIG), surveys=str(sp.DEFAULT_SURVEYS), force_disabled=False, window=None,
             utc_hours=None, calibration=str(HERE / "calibration"), survey_dirs=list(sp.DEFAULT_SURVEY_DIRS),
             output_root="/x", setup_coef=sp.SETUP_COEF, photo_roots=list(sp.DEFAULT_PHOTO_ROOTS), adcp=None,
             adcp_navd88=None, chatham=None, wis=None, ndbc=None, live_contours=str(sp.DEFAULT_LIVE_CONTOURS),
             setup_fitted_to=None, geotiff_epsg=sp.GEOTIFF_EPSG, gnssr_spline=None, gauge_csv=None,
             live_cron=str(sp.LIVE_CRON), pointing_log_dir=str(HERE / "archive"), pointing_frames=3,
             pointing_fail=sp.POINT_FAIL_DEG, dem_cell=sp.DEM_CELL, dem_min_points=sp.DEM_MIN_POINTS,
             dem_max_spread=sp.DEM_MAX_SPREAD, dem_max_hs=sp.DEM_MAX_HS, dem_max_day_offset=sp.DEM_MAX_DAY_OFFSET,
             keep_moved_days=False, no_download=False)
    a.update(kw)
    return SimpleNamespace(**a)


def test_historical_detect(d):
    print("survey_products historical eras: --correct-pointing detects 27 Sep (no longer --skip-days) and maps it "
          "with its own pointing")
    rows, truth = _day_rows()
    allc = Path(d) / "hist_contours.csv"
    _contours(allc, rows)
    seen = {}
    real_run = sp.run_cmd

    def fake_run(cmd, log_path, cwd=None):
        """detect_original_view.py as far as these checks need: the window's rows minus --skip-days, georectified
        with --eo by georectify.py itself (the detector is not what is tested here)"""
        words = [str(c) for c in cmd]
        if not words[1].endswith("detect_original_view.py"):
            return real_run(cmd, log_path, cwd)
        skip = sp.opt_after(words, "--skip-days", many=True) or []
        outd = Path(sp.opt_after(words, "--output"))
        outd.mkdir(parents=True, exist_ok=True)
        seen.setdefault("skip", []).append(skip)
        keep = [r for r in rows if r["capture_time_utc"][:10] not in skip]
        _contours(outd / "contour_points.csv", keep)
        return subprocess.run([sys.executable, str(HERE / "georectify.py"), str(outd / "contour_points.csv"),
                               str(outd / "contour_points_ground.csv"), "--io-c2", str(CAL / "CACO05_c2_20240801_IO.yaml"),
                               "--eo-c2", sp.opt_after(words, "--eo")], stdout=subprocess.DEVNULL).returncode
    photos = [Path("/nowhere") / photo_name(datetime.strptime(day, "%Y-%m-%d").replace(hour=17, tzinfo=timezone.utc),
                                            "CACO05", "c2")
              for day in ("2026-09-27", "2026-09-28", "2026-09-29", "2026-09-30")]
    saved = (sp.run_cmd, sp.scan_photos, sp.link_all_photos)
    sp.run_cmd = fake_run
    sp.scan_photos = lambda plan, cam, args, skip_days=(): (list(photos), {"missing_roots": [], "duplicates": 0})
    sp.link_all_photos = lambda plan: 0
    got = {}
    try:
        for tag, on in (("default", False), ("corrected", True)):
            out = Path(d) / f"hist_{tag}"
            (out / "forcing").mkdir(parents=True)
            (out / "forcing" / "water_level.csv").write_text("time,water_level_navd88,source,sigma_m\n"
                                                             "2026-09-27T00:00:00Z,0.6,adcp,0.05\n")
            (out / "forcing" / "waves.csv").write_text("epoch,wvht_m,dpd_s,source,hs_sigma_m\n0,1.0,8.0,adcp,0.1\n")
            env = out / "envelope.tif"
            env.write_text("x")
            plan = _live_plan(out)
            plan["era"] = "adcp"
            plan["cams"]["c2"].update(envelope_path=str(env), envelope_survey="envelope.tif")
            args = cfg_args(force=False, correct_pointing=on, setup_coef=0.0, keep_detector_debug=True)
            state = {"pointing": {"c2": _res()}}
            res, log = quiet(sp.step_detect, plan, args, state)
            got[tag] = (res, sp.read_stamp(out, "detect_c2"),
                        list(csv.DictReader(open(out / "waterlines" / "c2" / "contour_points_ground.csv", newline=""))),
                        log)
    finally:
        sp.run_cmd, sp.scan_photos, sp.link_all_photos = saved
    check(seen["skip"] == [["2026-09-27", "2026-09-28"], ["2026-09-28"]],
          "the detection's --skip-days: 27 and 28 Sep by default; only 28 Sep with --correct-pointing")
    _, st0, g0, _ = got["default"]
    _, st1, g1, _ = got["corrected"]
    check("pointing_correction" not in st0 and "pointing_correction" not in st0["sig"]
          and list(st1["pointing_correction"]["days"]) == ["2026-09-27"],
          "the stamps: none by default; 27 Sep corrected with the option")
    E, N = truth["2026-09-27"]
    r27 = [r for r in g1 if r["capture_time_utc"].startswith("2026-09-27")]
    err = np.hypot(np.array([float(r["easting_utm19"]) for r in r27]) - E,
                   np.array([float(r["northing_utm19"]) for r in r27]) - N)
    check(err.max() < 0.02 and [r for r in g1 if not r["capture_time_utc"].startswith("2026-09-27")] == g0,
          f"27 Sep mapped with its own pointing (to {err.max() * 1000:.0f} mm); every other row as by default")


def test_check_pointing_day_held(d):
    print("survey_products.check_pointing --correct-pointing: the day's first and last photos matched too; a camera "
          "that turned at 15:00 is not corrected (real matches)")
    import cv2
    d = Path(d) / "dayheld"
    img_dir = d / "photos"
    img_dir.mkdir(parents=True)
    tex = pdl.textured_scene(IO2, seed=6)
    turned = pdl.rotated_view(tex, IO2, EO2, pdl.corrected_eo(EO2, STORM_D))
    photos = []
    # 27 Sep turned all day; 28 Sep at rest at 14:00, turned from 15:00; 29 Sep (the reference day) at rest
    for day, plan_ in (("2026-09-27", (1, 1, 1)), ("2026-09-28", (0, 1, 1)), ("2026-09-29", (0, 0, 0))):
        for hour, on in zip((14, 17, 19), plan_):
            t = datetime.strptime(day, "%Y-%m-%d").replace(hour=hour, tzinfo=timezone.utc)
            img = turned if on else tex
            # each frame its own sensor noise (two copies of one image match at no noise at all)
            noise = np.random.default_rng(calendar.timegm(t.timetuple())).normal(0, 3.0, img.shape)
            p = img_dir / photo_name(t, "CACO05", "c2")
            cv2.imwrite(str(p), np.clip(img + noise, 0, 255).astype(np.uint8))
            photos.append(p)
    got = {}
    for tag, on in (("default", False), ("corrected", True)):
        plan = _live_plan(d / tag)
        args = cfg_args(correct_pointing=on)
        (rows, info), _ = quiet(sp.check_pointing, plan, "c2", photos, args)
        got[tag] = ({r["date"]: r for r in rows}, info, sp.pointing_sig(plan, "c2", args, photos))
    rows0, _, sig0 = got["default"]
    rows1, info1, sig1 = got["corrected"]
    check(not any("day_check" in r or "feature_time_utc" in r for r in rows0.values()) and "day_check" not in sig0
          and sig1["day_check"]["rules"] == sp.pointing_correction_rules(),
          "default: one probe frame a day, as before (no new column, no new key in the step's signature)")
    check(info1["reference_day"] == "2026-09-29" and "day_check" not in rows1["2026-09-29"]
          and [f["time_utc"] for f in rows1["2026-09-27"]["day_check"]] == ["14:00", "19:00"]
          and [f["time_utc"] for f in rows1["2026-09-28"]["day_check"]] == ["14:00", "19:00"],
          "--correct-pointing: the 14:00 and 19:00 photos of each day the 17:00 probe would correct are matched "
          "against the same reference (not the reference day)")
    # the features' gates (a synthetic scene has no sea horizon: the verdicts are not what is tested here)
    corr, why = pdl.window_corrections([{k: v for k, v in r.items() if k != "verdict"} for r in rows1.values()],
                                       sp.pointing_correction_rules(), info1["reference_day"])
    c27 = corr.get("2026-09-27") or {}
    check(list(corr) == ["2026-09-27"] and max(abs(a - b) for a, b in zip(c27["d_deg"], STORM_D)) < 0.02
          and c27["day_spread_deg"] < 0.02,
          f"27 Sep (turned all day): corrected by {c27.get('d_deg')}, its 14:00 and 19:00 frames within "
          f"{c27.get('day_spread_deg')} deg")
    w28 = "; ".join(why.get("2026-09-28") or [])
    check(rows1["2026-09-28"]["features"] == "MOVED" and "the 14:00 UTC frame gives no reading" in w28
          and "not shown to hold" in w28,
          f"28 Sep (turned at 15:00): its 17:00 probe reads the turn, but not corrected ({w28[:100]})")


def test_compare_c0(d):
    print("survey_compare C = 0: a corrected day's lines re-projected with that day's EO; the default command as "
          "before")
    import survey_compare as sc
    out = Path(d) / "live_corrected"
    g = out / "waterlines" / "c2" / "contour_points_ground.csv"
    eo_f = out / "waterlines" / "c2" / "pointing_corrected_eo" / "c2_2026-09-27_pointing_corrected_EO.yaml"
    cont = sc.read_contours(str(g))
    cams = sc.load_cameras({"c2": CAL / "CACO05_c2_20251113_EO-CV.yaml"}, {})
    plain, how0 = sc.without_setup(cont, cams)
    days = sc.parse_cam_day_specs([f"c2:2026-09-27={eo_f}"], "--camera-eo-day")
    cams["c2"]["eo_by_day"] = {k: load_extrinsics(v) for k, v in days["c2"].items()}
    c0, how = sc.without_setup(cont, cams)
    on = np.array([str(x).startswith("2026-09-27") for x in cont["day"]])
    E, N = pixel_to_ground(cont["U"][on], cont["V"][on], cont["Zt"][on], IO2, load_extrinsics(eo_f))
    E0, N0 = pixel_to_ground(cont["U"][on], cont["V"][on], cont["Zt"][on], IO2, EO2)
    check(days == {"c2": {"2026-09-27": eo_f}} and np.abs(c0["E"][on] - E).max() < 1e-9
          and np.abs(c0["N"][on] - N).max() < 1e-9 and np.median(np.hypot(plain["E"][on] - E, plain["N"][on] - N)) > 1,
          f"27 Sep's C = 0 lines are pixel_to_ground with the day's EO (with the calibration they would sit "
          f"{np.median(np.hypot(E0 - E, N0 - N)):.1f} m off: the correction undone)")
    check(np.array_equal(c0["E"][~on], plain["E"][~on]) and np.array_equal(c0["N"][~on], plain["N"][~on])
          and "(c2 2026-09-27) with that day's EO" in how and "--camera-eo-day" not in how0,
          "every other day's C = 0 lines as with the calibration alone; the how-text says which days")
    plan = _live_plan(out)
    check(sp.pointing_eo_day_specs(plan, ["c2"]) == [f"c2:2026-09-27={eo_f}"]
          and sp.pointing_eo_day_specs(_live_plan(Path(d) / "live_default"), ["c2"]) == []
          and sp.pointing_eo_day_specs(plan, []) == [],
          "survey_products passes --camera-eo-day for a corrected detection's days only (none for a default build, "
          "so its compare command is as it was)")


def test_provenance_readme_lines(d):
    print("survey_products: provenance and the README's camera line carry the correction only when it was made")
    src = open(HERE / "survey_products.py").read()
    blk = src.split("def write_provenance")[1].split("\ndef ")[0]
    check('if cams[cam]["pointing_check"] is not None and "pointing_correction" in (dst or {}):' in blk,
          "provenance: pointing_check.pointing_correction added only when the detection's stamp has it")
    rd = src.split("def write_readme")[1].split("\ndef ")[0]
    check('pcr = pc_.get("pointing_correction")' in rd and "if pcr is not None:" in rd,
          "README: the 'pointing CORRECTED per day' line only for such a product")


# --- the live log ---------------------------------------------------------------------------------------------

def test_live_log(d):
    print("the live log: complete days, its own reference, the gates per frame, nothing applied")
    import cv2
    d = Path(d) / "live_log"
    img_dir, log_dir, bank = d / "images_timex", d / "archive", d / "bank"
    for p in (img_dir, log_dir, bank / "c2"):
        p.mkdir(parents=True)
    tex = pdl.textured_scene(IO2, seed=4)
    today = datetime.now(timezone.utc).date()
    D = [(today - timedelta(days=k)).isoformat() for k in (5, 4, 3, 2)]       # ref day, turned, still, too far
    plant = {D[1]: (0.30, -0.20, 0.10), D[2]: (0.0, 0.0, 0.0), D[3]: (3.0, 0.0, 0.0)}

    def put(day, hour, turn, station="CACO05"):
        """a frame of the scene, turned by `turn`, with its own sensor noise (two copies of one image match at no
        noise at all, which no camera gives)"""
        t = datetime.strptime(day, "%Y-%m-%d").replace(hour=hour, tzinfo=timezone.utc)
        img = tex if turn is None else pdl.rotated_view(tex, IO2, EO2, pdl.corrected_eo(EO2, turn))
        noise = np.random.default_rng(calendar.timegm(t.timetuple()) + len(station)).normal(0, 3.0, img.shape)
        cv2.imwrite(str(img_dir / photo_name(t, station, "c2")), np.clip(img + noise, 0, 255).astype(np.uint8))
    put(D[0], 17, None)
    put(D[1], 14, plant[D[1]])
    put(D[1], 19, plant[D[1]])
    put(D[2], 16, None)
    put(D[3], 16, plant[D[3]])
    put(D[3], 15, (1.0, 0, 0), station="CACO04")              # another station: never used
    put(D[1], 3, (2.0, 0, 0))                                # outside the UTC hours: never used
    base = [sys.executable, str(HERE / "pointing_daily.py"), "--camera", "c2", "--image-dir", str(img_dir),
            "--log-dir", str(log_dir), "--bank-dir", str(bank)]
    r = subprocess.run(base + ["--days", D[0], D[3], "--reference-day", D[0]], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True)
    log = json.loads((log_dir / "pointing_daily_c2.json").read_text())
    days = log["days"]
    e1 = days[D[1]]
    check(r.returncode == 0 and log["applied_live"] is False and log["reference"]["kind"] == "reference day"
          and log["calibration"]["eo_file"] == "CACO05_c2_20251113_EO-CV.yaml" and log["rules"] == pdl.station_rules(),
          "the log says what it is: LOG ONLY (applied_live false), its reference (--reference-day), calibration, gates")
    check(e1["status"] == "moved" and e1["frames_passing"] == 2 and len(e1["frames"]) == 2
          and max(abs(a - b) for a, b in zip(e1["d_deg"], plant[D[1]])) < 0.02
          and max(abs(a - b) for a, b in zip(e1["would_correct_eo_deg"],
                                             np.degrees(EO2[3:6]) + np.array(plant[D[1]]))) < 0.02,
          f"{D[1]}: the planted 0.30/-0.20/0.10 deg read {e1['d_deg']} on both frames (the 03:00 frame and the other "
          f"station's left out), 'moved', with the EO angles a correction would use")
    check(days[D[2]]["status"] == "no reading" and any("within its own noise" in w for w in days[D[2]]["why_not"]),
          f"{D[2]}: the camera did not turn, the match is within its own noise: no reading, no correction")
    f3 = days[D[3]]["frames"][0]
    check(days[D[3]]["status"] == "no reading" and abs(f3["d_deg"][0] - 3.0) < 0.05
          and any("over 2.5" in w for w in f3["why_not"]),
          f"{D[3]}: a 3.0 deg turn is measured ({f3['d_deg'][0]:+.2f}) and logged, but over the 2.5 deg limit: no "
          f"reading a correction would use")
    check(days[D[0]]["frames"][0]["why_not"] == ["the reference frame itself"],
          f"{D[0]}: the reference frame itself is never measured against itself")
    check((log_dir / "pointing_daily_reference" / "c2").is_dir() and not list(bank.glob("c2/*"))
          and f"pointing c2 {D[1]}: moved" in r.stdout and "log only" in r.stdout and "WARNING" not in r.stdout,
          "its own reference in archive/pointing_daily_reference/c2/ (not the monitor's bank); one line a day in the "
          "log, no WARNING")
    check(e1["calibration"]["eo_file"] == "CACO05_c2_20251113_EO-CV.yaml" and e1["calibration"]["eo_sha256"]
          and e1["reference"] == {"kind": "reference day", "day": D[0], "frames": log["reference"]["frames"]},
          "each day records the calibration and the reference it was measured against")
    stamp = e1["checked_utc"]
    r = subprocess.run(base + ["--max-days", "2", "--reference-day", D[2]], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True)
    check(f"--reference-day {D[2]} not used: the log keeps its reference" in r.stdout and "remove " in r.stdout,
          "a new --reference-day with a reference already made: said, with how to replace it (not silently ignored)")
    log = json.loads((log_dir / "pointing_daily_c2.json").read_text())
    newd = sorted(set(log["days"]) - set(days))
    want = sorted([(today - timedelta(days=1)).isoformat(), (today - timedelta(days=6)).isoformat()])
    check(r.returncode == 0 and newd == want and all(log["days"][x]["status"] == "no frames" for x in newd)
          and log["days"][D[1]]["checked_utc"] == stamp,
          f"the cron's call: the newest 2 complete days not yet logged ({', '.join(newd)}: no frames), the logged "
          f"days not measured again, today never")
    r = subprocess.run(base[:-4] + ["--log-dir", str(d / "other"), "--bank-dir", str(bank)], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True)
    lg = json.loads((d / "other" / "pointing_daily_c2.json").read_text())
    check(r.returncode == 0 and lg["reference"]["kind"] is None and "no reference" in lg["reference"]["why"]
          and "--reference-day" in lg["reference"]["how"] and lg["days"] == {} and "nothing measured" in r.stdout,
          "no reference (no bank frame at this calibration, no photo of its own day 2025-11-13): nothing measured, the "
          "log says how to give one, exit 0")


def test_cron():
    print("waterline_timex_cron.sh: the log runs log-only, with the calibrations its own georectification uses")
    src = (HERE / "waterline_timex_cron.sh").read_text()
    call = src.split('python3 "$BASE/pointing_daily.py"')[1].split("fi\n")[0]
    geo = src.split('python3 "$BASE/georectify.py" "$MAP_CONTOURS" "$GROUND"')[1].split(">>")[0]
    eos = lambda t: re.findall(r'--eo-(c\d) "(\$CAL/[^"]+)"', t)               # noqa: E731
    check(eos(call) == eos(geo) and len(eos(call)) == 2, f"the same EO files as step 6 ({eos(call)})")
    check("POINTING_LOG_ENABLE=1" in src and 'if [ "$POINTING_LOG_ENABLE" = "1" ]' in src and "--correct" not in call
          and "LOG ONLY" in src and "WARNING: camera pointing log failed" in call,
          "on by default, LOG ONLY (no correction switch), its failure a WARNING")
    r = subprocess.run(["bash", "-n", str(HERE / "waterline_timex_cron.sh")])
    check(r.returncode == 0, "bash -n waterline_timex_cron.sh")


def test_py38():
    print("Python 3.8 syntax")
    import ast
    for f in ("pointing_daily.py", "test_pointing_daily.py", "survey_products.py"):
        src = (HERE / f).read_text(encoding="utf-8")
        try:
            ast.parse(src, feature_version=(3, 8))
        except TypeError:                        # Python 3.7 has no feature_version
            ast.parse(src)
        check(True, f"{f} parses as Python 3.8")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="keep the outputs, in a new folder inside this one")
    args = ap.parse_args()
    if args.keep:
        Path(args.keep).mkdir(parents=True, exist_ok=True)
    d = os.path.realpath(tempfile.mkdtemp(prefix="pointing_daily_test_", dir=args.keep))
    start = os.getcwd()
    os.chdir(d)
    print(f"test folder: {d} (the current folder during the run)")
    try:
        test_gates()
        test_corrected_eo_exact()
        test_self_test()
        test_dirt_and_fog_gated()
        test_regeorectify(d)
        test_apply_corrections(d)
        test_skip_days()
        test_root_rule()
        test_live_detect(d)
        test_historical_detect(d)
        test_compare_c0(d)
        test_check_pointing_day_held(d)
        test_provenance_readme_lines(d)
        test_live_log(d)
        test_cron()
        test_py38()
    finally:
        os.chdir(start)
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    print("all pointing_daily tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

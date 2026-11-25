#!/usr/bin/env python3
"""
Tests For pair, profile, trial and pointing (waterline_station_check.py)
==========================================================================
Plain python3 (no pytest; the station's Python 3.8 runs it too):

    python3 test_waterline_station_check.py [--keep DIR]

WHAT. coastcam-station-setup (./setup.sh --waterline) reads these commands'
JSON on stdout (a summary for people goes to stderr; exit 0 ok, 1 refused or
failed, 2 usage). These checks run the script as setup does:
  * pair (design 6.4; C5, C6): CACO-05's own 2025-11-13 calibration, its
    lenses, positions, targets and picks named one by one: the right sets pass
    (mean under 1 m) with the next best position at least 3 times worse; the
    swapped position fails clearly; a duplicate lens is an io_tie, never a
    refusal; a duplicate position (the picks cannot tell them apart), a wrong
    picture size, too few shared targets and a picks file without its header
    are refused; the target range and height checks (camera z + 5 m);
  * profile (6.5; C8): the search band from a synthetic beach survey (points,
    with or without a header, and as a DEM .asc), written as
    detect_original_view.py writes it; from the camera geometry alone
    (PROVISIONAL, with a STATION_PROVISIONAL text check-overlay accepts); the
    pinned envelope_from_survey's sys.exit ("the band hardly appears")
    caught; a DEM that does not cover the camera; --merge; usage errors;
    geometry_grid's 40-600 m ring itself;
  * trial (6.6): a stub detector (it records its options) with the real
    extract_elevation_contours.py and georectify.py, five frames, one of each
    outcome; the scratch-folder rules (never input/, processed_timex/,
    archive/ ...; the pictures only read); a relative --work-dir; the antenna
    height from --antenna-navd88 or the environment, never CACO-05's 19.014
    without a word; then the real detector on one synthetic frame (needs
    OpenCV);
  * pointing (correction 7): the PASS/WARN/FAIL limits (pointing_grade, with
    or without OpenCV); planted turns of a synthetic scene give PASS,
    WARN and FAIL; unrelated pictures give UNMATCHED (nothing checked), a
    picture of another size an ERROR, all as JSON with exit 1 (needs OpenCV:
    a SKIP line without it).

STATION FILES. Only these calibration files are read, each named:
calibration/CACO05_c{1,2}_20240801_IO.yaml, CACO05_c{1,2}_20251113_EO-CV.yaml,
2025-11-13_Marconi_Extrinsic_Targets_c{1,2}_xyz.csv and gcp_uv_c{1,2}.csv.
Everything else is made in a temporary folder.
"""

import os
import sys
import json
import math
import shutil
import argparse
import tempfile
import subprocess
from pathlib import Path
from datetime import datetime, timezone

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import waterline_station_check as wsc                              # noqa: E402
from georectify import load_intrinsics, load_extrinsics, pixel_to_ground   # noqa: E402

CHECK = HERE / "waterline_station_check.py"
CAL = HERE / "calibration"
IO = {c: CAL / f"CACO05_{c}_20240801_IO.yaml" for c in ("c1", "c2")}
EO = {c: CAL / f"CACO05_{c}_20251113_EO-CV.yaml" for c in ("c1", "c2")}
XYZ = {c: CAL / f"2025-11-13_Marconi_Extrinsic_Targets_{c}_xyz.csv" for c in ("c1", "c2")}
UV = {c: CAL / f"gcp_uv_{c}.csv" for c in ("c1", "c2")}

FAILED = []
try:
    import cv2                                                      # noqa: F401
    HAVE_CV2 = True
except ImportError:
    HAVE_CV2 = False


def check(msg, cond, detail=None):
    print(("  ok    " if cond else "  FAIL  ") + msg)
    if not cond:
        FAILED.append(msg)
        if detail is not None:
            print("        " + str(detail)[:1500].replace("\n", "\n        "))


def strict_json(text):
    """JSON as setup reads it: NaN or Infinity anywhere is an error."""
    def no(c):
        raise ValueError(f"non-standard JSON constant {c}")
    return json.loads(text, parse_constant=no)


def cli(*args, cwd=None, env=None):
    """-> (exit status, JSON dict or None, stdout, stderr)."""
    p = subprocess.run([sys.executable, str(CHECK)] + [str(a) for a in args], stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True, cwd=str(cwd or HERE), env=env)
    try:
        js = strict_json(p.stdout) if p.stdout.strip().startswith("{") else None
    except ValueError:
        js = None
    return p.returncode, js, p.stdout, p.stderr


def frame_name(epoch, cam, station="TEST07"):
    u = datetime.fromtimestamp(epoch, tz=timezone.utc)
    return f"{epoch}.{u.strftime('%a.%b.%d_%H_%M_%S')}.GMT.{u.year}.{station}.{cam}.timex.jpg"


T0 = 1789000200                       # 2026-09-10 00:30 UTC


def write_spline(path, t0=T0, hours=6, hortho=19.014):
    """gnssrefl subdaily spline_out.txt layout, Hortho 19.014 by default (no datum shift at CACO-05's
    antenna height)."""
    ts = np.arange(t0 - hours * 3600, t0 + hours * 3600 + 1, 1800)
    with open(path, "w") as f:
        f.write(f"% station test gnssrefl\n% orthometric height minus RH, where Hortho (m) is {hortho:9.3f} \n"
                "% MJD RH(m) YYYY MM DD HH MM SS quasi-sea-level(m)\n")
        for t in ts:
            z = 0.3 + 1.0 * math.sin(2 * math.pi * (t - ts[0]) / (12.42 * 3600))
            u = datetime.fromtimestamp(int(t), tz=timezone.utc)
            f.write(f"{t / 86400 + 40587:15.7f} {hortho - z:10.3f} {u.year:4d} {u.month:3d} {u.day:3d} "
                    f"{u.hour:3d} {u.minute:3d} {u.second:3d} {z:10.3f}\n")


def beach_points(cam="c1", behind=False, step_along=10.0, step_across=4.0):
    """A synthetic beach in front of (or behind) the camera: z = 4 - 0.04 (d - 30) m NAVD88 at d m along
    the camera's azimuth, profiles every step_along m across 300 m. -> (E, N, Z) arrays."""
    eo = load_extrinsics(EO[cam])
    cx, cy, az = eo[0], eo[1], eo[3]
    sign = -1.0 if behind else 1.0
    E, N, Z = [], [], []
    for along in np.arange(-150, 151, step_along):
        for d in np.arange(20, 260, step_across):
            E.append(cx + sign * d * math.sin(az) + along * math.cos(az))
            N.append(cy + sign * d * math.cos(az) - along * math.sin(az))
            Z.append(4.0 - 0.04 * (d - 30))
    return np.array(E), np.array(N), np.array(Z)


# ---------------------------------------------------------------------
# pair
# ---------------------------------------------------------------------

def test_pair(d):
    print("pair: CACO-05's 2025-11-13 calibration (files named one by one)")
    for cam, other in (("c1", "c2"), ("c2", "c1")):
        rc, js, out, err = cli("pair", "--xyz", XYZ[cam], "--uv", UV[cam], "--io", IO["c1"], IO["c2"],
                               "--eo", EO["c1"], EO["c2"], "--image", "2448x2048")
        ok = js is not None and rc == 0 and js["ok"] is True and isinstance(js["chosen"], dict)
        ch = js["chosen"] if ok else {}
        check(f"{cam}: the right set is chosen (exit 0): {IO[cam].name} + {EO[cam].name}, mean "
              f"{ch.get('mean_m')} m over {ch.get('n')} targets",
              ok and ch["io"] == str(IO[cam]) and ch["eo"] == str(EO[cam]) and ch["mean_m"] < 1.0
              and ch["n"] >= 13 and ch["io_decided"] is True and js["io_tie"] is None, (rc, out[:600], err))
        if not ok:
            continue
        check(f"{cam}: the next best position ({Path(js['next_eo']['eo']).name}) is {js['ratio_next_eo']} times "
              f"worse (3 or more needed)", js["ratio_next_eo"] >= 3 and js["next_eo"]["eo"] == str(EO[other]))
        combos = {(Path(c["io"]).name, Path(c["eo"]).name): c for c in js["combinations"]}
        check(f"{cam}: all 4 lens x position combinations reported with mean, max and n; only the right one passes",
              len(combos) == 4 and all(set(c) >= {"mean_m", "max_m", "n", "pass"} for c in combos.values())
              and [k for k, c in combos.items() if c["pass"]] == [(IO[cam].name, EO[cam].name)], combos)
        wrong_io = combos[(IO[other].name, EO[cam].name)]
        check(f"{cam}: the other lens with the right position fails ({wrong_io['mean_m']} m) but is no part of "
              f"the ratio test (C5: the ratio is for the position)", not wrong_io["pass"] and wrong_io["mean_m"] > 1)
        san = {s["check"]: s for s in js["sanity"]}
        check(f"{cam}: sanity after the choice: picture size, picks in frame, targets 10-1000 m, target Z up to "
              f"camera z + 5 m: {'; '.join(s['detail'] for s in js['sanity'])}",
              set(san) == {"image_size", "picks_in_frame", "target_range", "target_z"}
              and all(s["ok"] is True for s in san.values()))
        check(f"{cam}: the camera's position and angles in the JSON (for naming the cameras)",
              set(ch["camera"]) == {"x", "y", "z", "azimuth", "tilt", "roll"} and ch["lens"]["NU"] == 2448)
        check(f"{cam}: a summary for people on stderr, JSON only on stdout",
              "pair: OK" in err and "chosen:" in err and out.lstrip().startswith("{"), err)

    rc, js, out, err = cli("pair", "--xyz", XYZ["c1"], "--uv", UV["c1"], "--io", IO["c1"], IO["c2"], "--eo", EO["c2"])
    best = min((c["mean_m"] for c in js["combinations"]), default=None) if js else None
    check(f"the swapped position (c1's picks, c2's EO) fails clearly: exit 1, chosen 'none', best mean {best} m",
          rc == 1 and js and js["chosen"] == "none" and js["ok"] is False and best > 10
          and "does not pass its test" in js["reason"] and "REFUSED" in err, (rc, out[:400], err))
    check("... with one position given there is no ratio to report", js and js["ratio_next_eo"] is None)

    # an io_tie: the same lens twice (a copy under another name): never a refusal (C5)
    dup = d / "c1_lens_copy_IO.yaml"
    dup.write_text("# a copy of c1's lens under another name (test)\n" + IO["c1"].read_text())
    rc, js, out, err = cli("pair", "--xyz", XYZ["c1"], "--uv", UV["c1"], "--io", IO["c1"], dup, IO["c2"],
                           "--eo", EO["c1"], EO["c2"], "--image", "2448x2048")
    tie = (js or {}).get("io_tie") or {}
    check("two lenses passing with one position within 20%: an io_tie listing both, exit 0, not refused",
          rc == 0 and js["ok"] is True and sorted(e["io"] for e in tie.get("ios", [])) == sorted([str(IO["c1"]), str(dup)])
          and tie.get("eo") == str(EO["c1"]) and js["chosen"]["io_decided"] is False and "io_tie" in err,
          (rc, out[:800], err))
    check("... the ratio test still on the positions (the next position many times worse)",
          js and js["ratio_next_eo"] >= 3)
    check("... each tied lens says whether it fits the picture size",
          all(e["image_size_ok"] is True and e["nu_nv"] == [2448.0, 2048.0] for e in tie.get("ios", [])), tie)

    # two positions the picks cannot tell apart: refused, both listed
    dup_eo = d / "c1_position_copy_EO.yaml"
    dup_eo.write_text(EO["c1"].read_text().replace("azimuth: 54.4904", "azimuth: 54.4954"))
    rc, js, out, err = cli("pair", "--xyz", XYZ["c1"], "--uv", UV["c1"], "--io", IO["c1"], "--eo", EO["c1"], dup_eo)
    check(f"two positions fitting the picks equally (ratio {js and js['ratio_next_eo']}): refused, exit 1, both "
          f"in eo_tie", rc == 1 and js["chosen"] == "none" and js["ratio_next_eo"] < 3
          and len(js["eo_tie"] or []) == 2 and "keep one position file" in js["reason"], (rc, out[:600]))

    rc, js, out, err = cli("pair", "--xyz", XYZ["c1"], "--uv", UV["c1"], "--io", IO["c1"], IO["c2"],
                           "--eo", EO["c1"], EO["c2"], "--image", "1920x1080")
    san = {s["check"]: s for s in (js or {}).get("sanity") or []}
    check("a lens for another picture size: the set is chosen, the size check fails, exit 1",
          rc == 1 and js["ok"] is False and isinstance(js["chosen"], dict) and san["image_size"]["ok"] is False
          and "image_size" in js["reason"], (rc, out[:600]))

    few = d / "uv_three.csv"
    few.write_text("".join(UV["c1"].read_text().splitlines(True)[:4]))
    rc, js, out, err = cli("pair", "--xyz", XYZ["c1"], "--uv", few, "--io", IO["c1"], "--eo", EO["c1"])
    check("3 shared targets: not paired (4 or more needed), exit 1",
          rc == 1 and js["n_shared"] == 3 and js["chosen"] == "none" and "not paired" in js["reason"], (rc, out[:400]))

    noheader = d / "uv_noheader.csv"
    noheader.write_text("".join(UV["c1"].read_text().splitlines(True)[1:]))
    rc, js, out, err = cli("pair", "--xyz", XYZ["c1"], "--uv", noheader, "--io", IO["c1"], "--eo", EO["c1"])
    check("a picks file without its num,U,V header: refused with the reason, exit 1",
          rc == 1 and js and "header must be num,U,V" in js["reason"], (rc, out[:300], err))

    # the sanity checks themselves (C6): a target on a bluff top is fine up to camera z + 5 m
    xyz = wsc.read_targets(XYZ["c1"])
    uv = wsc.read_picks(UV["c1"])
    shared = sorted(set(xyz) & set(uv))
    z_cam = load_extrinsics(EO["c1"])[2]
    n0 = shared[0]
    mod = dict(xyz)
    mod[n0] = (xyz[n0][0], xyz[n0][1], z_cam + 4.0)
    s = {c["check"]: c for c in wsc.pair_sanity(IO["c1"], EO["c1"], mod, uv, shared)}
    check(f"a target at camera z + 4 m ({z_cam + 4:.1f} m, a bluff top): accepted (not a fixed +15 m cap)",
          s["target_z"]["ok"] is True and s["image_size"]["ok"] is None, s["target_z"])
    mod[n0] = (xyz[n0][0], xyz[n0][1], z_cam + 6.0)
    s = {c["check"]: c for c in wsc.pair_sanity(IO["c1"], EO["c1"], mod, uv, shared)}
    check("a target above camera z + 5 m: refused, named", s["target_z"]["ok"] is False
          and f"{n0} (" in s["target_z"]["detail"], s["target_z"])
    mod = dict(xyz)
    mod[n0] = (xyz[n0][0] + 1500.0, xyz[n0][1], xyz[n0][2])
    s = {c["check"]: c for c in wsc.pair_sanity(IO["c1"], EO["c1"], mod, uv, shared)}
    check("a target 1500 m from the camera: refused (10-1000 m)", s["target_range"]["ok"] is False, s["target_range"])
    mod_uv = dict(uv)
    mod_uv[n0] = (3000.0, uv[n0][1])
    s = {c["check"]: c for c in wsc.pair_sanity(IO["c1"], EO["c1"], xyz, mod_uv, shared)}
    check("a pick outside the 2448 x 2048 frame: refused", s["picks_in_frame"]["ok"] is False)

    rc, js, out, err = cli("pair", "--xyz", d / "no_such.csv", "--uv", UV["c1"], "--io", IO["c1"], "--eo", EO["c1"])
    check("a file that is not there: usage, exit 2, said in the JSON", rc == 2 and js and "not found" in js["error"])
    rc, js, out, err = cli("pair", "--xyz", XYZ["c1"], "--uv", UV["c1"], "--io", IO["c1"], "--eo", EO["c1"],
                           "--image", "big")
    check("--image not WxH: usage, exit 2", rc == 2 and "WIDTHxHEIGHT" in err, err)
    rc, js, out, err = cli("pair", "--xyz", XYZ["c1"])
    check("--uv, --io, --eo missing: usage, exit 2", rc == 2)


# ---------------------------------------------------------------------
# profile
# ---------------------------------------------------------------------

def _profile_shape(path, cam):
    """The file as detect_original_view.py:343-347 writes it -> (ok, profile)."""
    data = strict_json(Path(path).read_text())
    prof = data.get(cam, {})
    keys = set(prof)
    env = prof.get("envelope") or []
    ok = (keys in ({"crop", "envelope", "no_bias"}, {"crop", "envelope", "no_bias", "step", "profile_scale"})
          and prof["crop"] == [0.0, 1.0, 0.0, 1.0] and prof["no_bias"] is True and len(env) >= 5
          and all(len(p) == 3 and 0 <= p[1] <= p[2] <= 1 for p in env)
          and env[0][0] == 0.0 and env[-1][0] == 1.0 and [p[0] for p in env] == sorted(p[0] for p in env)
          and (("step" in prof) == ("profile_scale" in prof))
          and prof.get("step", [4, 18]) == [4, 18] and prof.get("profile_scale", [0.4, 1.0]) == [0.4, 1.0])
    return ok, prof


def test_profile(d):
    print("profile: the detector's search band")
    common = ["--camera", "c1", "--io", IO["c1"], "--eo", EO["c1"], "--image-size", "2448x2048"]
    E, N, Z = beach_points("c1")
    survey = d / "survey.csv"                        # header-less num,E,N,Z (the normalised form)
    survey.write_text("".join(f"{k},{e:.3f},{n:.3f},{z:.3f}\n" for k, (e, n, z) in enumerate(zip(E, N, Z), 1)))
    out_json = d / "profile_survey.json"
    rc, js, out, err = cli("profile", *common, "--survey-csv", survey, "--mllw", "-0.9", "--mhhw", "1.1",
                           "--output", out_json)
    ok, prof = _profile_shape(out_json, "c1") if out_json.exists() else (False, {})
    check(f"from a beach survey ({len(E)} points): exit 0, status OK, written as detect_original_view.py writes it "
          f"({len(prof.get('envelope', []))} envelope points)",
          rc == 0 and js["ok"] is True and js["status"] == "OK" and ok and js["source"] == "survey"
          and js["provisional_reason"] is None and js["profile"] == prof, (rc, out[:600], err))
    check(f"band from the gauge: MLLW - 0.5 to MHHW + 1.0 = {js and js['band_m']}, within 300 m",
          js and js["band_m"] == [-1.4, 2.1] and js["range_m"] == 300.0)
    check("the file has only this camera, with JSON's indent=1 like detect_original_view.py",
          list(strict_json(out_json.read_text())) == ["c1"] and out_json.read_text().startswith('{\n "c1": {'))
    check("the note says the other CACO05 tunings stay (band_contrast_floor, snap settings)",
          js and "band_contrast_floor" in js["note"] and "snap settings" in js["station_provisional"])

    # the band follows the survey: its +2.5 m and -1.5 m contours (67.5 and 167.5 m along the view) lie
    # inside the envelope at their own columns
    io_, eo = load_intrinsics(IO["c1"]), load_extrinsics(EO["c1"])
    from view_reproject import ground_to_pixel
    env = np.array(prof.get("envelope") or [[0, 0, 0], [1, 0, 0]], float)
    worst, n_in = 0.0, 0
    for dd, zz in ((67.5, 2.5), (167.5, -1.5)):
        al = np.arange(-100, 101, 20.0)
        E2 = eo[0] + dd * math.sin(eo[3]) + al * math.cos(eo[3])
        N2 = eo[1] + dd * math.cos(eo[3]) - al * math.sin(eo[3])
        U, V, okv = ground_to_pixel(E2, N2, np.full(len(al), zz), io_, eo)
        x, v = U[okv] / io_[0], V[okv] / io_[1]
        lo, hi = np.interp(x, env[:, 0], env[:, 1]), np.interp(x, env[:, 0], env[:, 2])
        worst = max([worst] + list(lo - v) + list(v - hi))
        n_in += int(okv.sum())
    check(f"the band's +2.5 and -1.5 m contours ({n_in} points in view) lie inside the envelope at their columns "
          f"(worst {worst:+.3f} of the height outside; 0.02 allowed)", n_in >= 10 and worst <= 0.02, worst)

    # with a header and the default band
    hdr = d / "survey_header.csv"
    hdr.write_text("Name,Easting,Northing,Elevation,CS name\n" + "".join(
        f"P{k},{e:.3f},{n:.3f},{z:.3f},NAD83(2011) / UTM 19N + NAVD88\n" for k, (e, n, z) in enumerate(zip(E, N, Z), 1)))
    rc, js, out, err = cli("profile", *common, "--survey-csv", hdr, "--output", d / "profile_header.json")
    check("a survey with a header (Name,Easting,Northing,Elevation): read by name; default band -1.5 to +2.5",
          rc == 0 and js["band_m"] == [-1.5, 2.5] and js["survey_points"] == len(E), (rc, out[:400], err))

    # a DEM (.asc), as compare_dem_survey.read_survey reads it
    cell = 1.0
    x0, y1 = math.floor(E.min()), math.ceil(N.max())
    nc, nr = int(math.ceil(E.max() - x0)) + 1, int(math.ceil(y1 - N.min())) + 1
    XX, YY = np.meshgrid(x0 + (np.arange(nc) + 0.5) * cell, y1 - (np.arange(nr) + 0.5) * cell)
    along = (XX - eo[0]) * math.sin(eo[3]) + (YY - eo[1]) * math.cos(eo[3])
    G = np.where((along > 20) & (along < 260), 4.0 - 0.04 * (along - 30), -9999.0)
    asc = d / "beach_dem.asc"
    with open(asc, "w") as f:
        f.write(f"ncols {nc}\nnrows {nr}\nxllcorner {x0}\nyllcorner {y1 - nr * cell}\ncellsize {cell}\n"
                f"NODATA_value -9999\n")
        np.savetxt(f, G, fmt="%.3f")
    rc, js, out, err = cli("profile", *common, "--dem", asc, "--mllw", "-0.9", "--mhhw", "1.1",
                           "--output", d / "profile_dem.json")
    ok, prof_dem = _profile_shape(d / "profile_dem.json", "c1") if (d / "profile_dem.json").exists() else (False, {})
    check("from a beach DEM (.asc): exit 0, the same shape", rc == 0 and js["source"] == "dem" and ok,
          (rc, out[:400], err))
    far = d / "far_dem.asc"
    with open(far, "w") as f:
        f.write(f"ncols 10\nnrows 10\nxllcorner {x0 + 50000}\nyllcorner {y1}\ncellsize 1\nNODATA_value -9999\n")
        np.savetxt(f, np.zeros((10, 10)), fmt="%.1f")
    rc, js, out, err = cli("profile", *common, "--dem", far, "--output", d / "profile_far.json")
    check("a DEM that does not cover the camera's ground: refused, exit 1, nothing written",
          rc == 1 and "does not cover" in js["reason"] and not (d / "profile_far.json").exists(), (rc, out[:300]))

    # the pinned envelope_from_survey calls sys.exit when the band hardly appears (C8): caught
    Eb, Nb, Zb = beach_points("c1", behind=True)
    behind = d / "survey_behind.csv"
    behind.write_text("".join(f"{e:.3f},{n:.3f},{z:.3f}\n" for e, n, z in zip(Eb, Nb, Zb)))
    rc, js, out, err = cli("profile", *common, "--survey-csv", behind, "--output", d / "profile_behind.json")
    check("a survey behind the camera: envelope_from_survey's SystemExit caught, refused with its words, exit 1",
          rc == 1 and js and js["ok"] is False and "hardly appears" in js["reason"] and "Traceback" not in err
          and not (d / "profile_behind.json").exists(), (rc, out[:400], err))

    # from geometry alone: PROVISIONAL
    geo = d / "detector_profile.json"
    rc, js, out, err = cli("profile", *common, "--geometry", "--plane-z", "0.12", "--output", geo)
    ok, prof_geo = _profile_shape(geo, "c1") if geo.exists() else (False, {})
    check("from the camera geometry (a level plane at the gauge's MSL, 40-600 m): exit 0, PROVISIONAL with "
          "the reason", rc == 0 and ok and js["status"] == "PROVISIONAL" and js["source"] == "geometry"
          and "camera geometry" in js["provisional_reason"] and "+0.12 m" in js["provisional_reason"]
          and js["ring_m"] == [40.0, 600.0] and "PROVISIONAL" in err, (rc, out[:600], err))
    # the ring itself (C8): the plane only 40-600 m from the camera, NaN nearer and farther. The envelope
    # above changes when the ring goes (up to 0.1 of the picture's height), but nothing else pins it.
    cx, cy, zp = 412345.6, 4645678.9, 0.12
    g, xl, yt, cell = wsc.geometry_grid(cx, cy, zp)
    XX, YY = np.meshgrid(xl + (np.arange(g.shape[1]) + 0.5) * cell, yt - (np.arange(g.shape[0]) + 0.5) * cell)
    rr = np.hypot(XX - cx, YY - cy)
    lo, hi = wsc.GEOMETRY_RING_M
    inside = (rr >= lo) & (rr <= hi)
    check(f"geometry_grid: a level plane at z only {lo:g}-{hi:g} m from the camera, NaN nearer and farther, "
          f"over a north-up box round the camera",
          (lo, hi) == (40.0, 600.0) and inside.any() and (rr < lo).any() and (rr > hi).any()
          and bool(np.isnan(g[~inside]).all()) and bool((g[inside] == zp).all())
          and abs(xl + g.shape[1] * cell / 2.0 - cx) < cell and abs(yt - g.shape[0] * cell / 2.0 - cy) < cell
          and g.shape[0] * cell >= 2 * hi, (g.shape, xl, yt, cell))
    sp = js["station_provisional"] if js else ""
    res = wsc.check_overlay(f'STATION_PROVISIONAL="{sp}"\nWATERLINE_STAGE="archive"\n')
    check(f"its STATION_PROVISIONAL text passes check-overlay: \"{sp[:70]}...\"",
          not [e for e in res["errors"] if e["key"] == "STATION_PROVISIONAL"]
          and "band_contrast_floor" in sp, res["errors"])
    rc, js, out, err = cli("profile", "--camera", "c2", "--io", IO["c2"], "--eo", EO["c2"], "--image-size",
                           "2448x2048", "--geometry", "--output", geo, "--merge")
    data = strict_json(geo.read_text()) if geo.exists() else {}
    check("--merge: c2 added to the same file, c1 kept as it was (one --profile-json for the detector)",
          rc == 0 and sorted(data) == ["c1", "c2"] and data["c1"] == prof_geo and js["plane_z_m"] == 0.0,
          (rc, sorted(data)))
    sky = d / "c1_looking_up_EO.yaml"
    sky.write_text(EO["c1"].read_text().replace("tilt: 66.7572", "tilt: 179.0"))
    rc, js, out, err = cli("profile", "--camera", "c1", "--io", IO["c1"], "--eo", sky, "--image-size", "2448x2048",
                           "--geometry", "--output", d / "p_up.json")
    check("a camera looking straight up: no ground in view from geometry either, SystemExit caught, exit 1",
          rc == 1 and "hardly appears" in js["reason"] and "Traceback" not in err
          and not (d / "p_up.json").exists(), (rc, out[:300], err))

    rc, js, out, err = cli("profile", "--camera", "c1", "--io", IO["c1"], "--eo", EO["c1"], "--image-size",
                           "1224x1024", "--geometry", "--output", d / "p_size.json")
    check("pictures of another size than the lens file: refused, exit 1",
          rc == 1 and "lens file is 2448 x 2048" in js["reason"], (rc, out[:300]))
    rc, js, out, err = cli("profile", *common, "--survey-csv", survey, "--mllw", "-0.9", "--output", d / "x.json")
    check("--mllw without --mhhw: usage, exit 2", rc == 2 and "go together" in err, err)
    rc, js, out, err = cli("profile", *common, "--survey-csv", survey, "--geometry", "--output", d / "x.json")
    check("a survey and --geometry together: usage, exit 2", rc == 2, err)
    rc, js, out, err = cli("profile", *common, "--output", d / "x.json")
    check("no band source: usage, exit 2", rc == 2, err)
    rc, js, out, err = cli("profile", *common, "--survey-csv", survey, "--plane-z", "0.1", "--output", d / "x.json")
    check("--plane-z without --geometry: usage, exit 2", rc == 2 and "--plane-z goes with --geometry" in err, err)
    rubbish = d / "not_utm.csv"
    rubbish.write_text("1,2,3\n4,5,6\n7,8,9\n")
    rc, js, out, err = cli("profile", *common, "--survey-csv", rubbish, "--output", d / "x.json")
    check("a survey that is not UTM: refused with the reason, exit 1", rc == 1 and "UTM" in js["reason"], out[:300])


# ---------------------------------------------------------------------
# trial
# ---------------------------------------------------------------------

STUB_DETECTOR = r'''#!/usr/bin/env python3
"""Stub of waterline_detector_v5.py for test_waterline_station_check.py: records its options and
writes the detections plan.json asks for, in the detector's CSV form."""
import argparse, csv, json, sys
from pathlib import Path
ap = argparse.ArgumentParser()
for o in ("--image-suffix", "--timing-mode", "--source-dir", "--input-dir", "--output-dir", "--debug-dir",
          "--profile-json"):
    ap.add_argument(o)
ap.add_argument("--no-bias-correction", action="store_true")
a = ap.parse_args()
here = Path(__file__).resolve().parent
(here / "detector_argv.json").write_text(json.dumps(sys.argv[1:]))
plan = json.loads((here / "plan.json").read_text())
for p in sorted(Path(a.source_dir).glob("*" + a.image_suffix)):
    (Path(a.input_dir) / p.name).write_bytes(p.read_bytes())
    stem = p.name[:-len(".jpg")]
    spec = plan.get(stem)
    if not spec:
        print("no usable waterline:", p.name)
        continue
    n, sig = spec["columns"], int(round(spec["signal"] * spec["columns"]))
    with open(Path(a.output_dir) / (stem + ".csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Column", "Row", "Confidence", "Row_Precise", "Peak_Sharpness", "Has_Signal"])
        for c in range(n):
            r = spec["row"] + 0.6 * (((c * 7919) % 11) - 5) / 5.0
            w.writerow([c, int(round(r)), 0.95, "%.3f" % r, 0.5, 1 if c < sig else 0])
'''


def _jpeg(path, w, h, row=900):
    """A synthetic timex: dark sea above `row`, bright sand below (blurred)."""
    if HAVE_CV2:
        img = np.full((h, w, 3), 90, np.uint8)
        img[row:, :] = 170
        cv2.imwrite(str(path), cv2.GaussianBlur(img, (0, 0), 6))
    else:
        from PIL import Image
        a = np.full((h, w, 3), 90, np.uint8)
        a[row:, :] = 170
        Image.fromarray(a).save(str(path), quality=90)


def _sha(path):
    import hashlib
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_trial(d):
    print("trial: a stub detector, the real contours and georectify, scratch folders only")
    src = d / "kept_timex"                           # stands for archive/images_timex: only read
    src.mkdir()
    names = [frame_name(T0 + k * 1800, "c1") for k in range(5)]
    for k, n in enumerate(names):
        _jpeg(src / n, *((1224, 1024) if k == 4 else (2448, 2048)))
    before = {n: _sha(src / n) for n in names}
    spline = d / "spline.txt"
    write_spline(spline)
    stubs = d / "stub_scripts"
    stubs.mkdir()
    (stubs / "waterline_detector_v5.py").write_text(STUB_DETECTOR)
    for s in ("extract_elevation_contours.py", "georectify.py"):
        os.symlink(str(HERE / s), str(stubs / s))
    io_, eo = load_intrinsics(IO["c1"]), load_extrinsics(EO["c1"])
    # frame 2's line so low in the frame that it lies under 30 m from the camera
    e, n = pixel_to_ground(np.array([1224.0]), np.array([2000.0]), 0.5, io_, eo)
    near_m = float(np.hypot(e[0] - eo[0], n[0] - eo[1]))
    stems = [x[:-len(".jpg")] for x in names]
    plan = {stems[0]: {"row": 900, "signal": 1.0, "columns": 2448},      # passes
            stems[1]: {"row": 900, "signal": 0.4, "columns": 2448},      # a line in 40% of the columns
            stems[2]: {"row": 2000, "signal": 1.0, "columns": 2448}}     # under 30 m from the camera
    (stubs / "plan.json").write_text(json.dumps(plan))               # frame 3: the detector keeps no line
    geo = d / "profile_c1.json"
    rc, js, out, err = cli("profile", "--camera", "c1", "--io", IO["c1"], "--eo", EO["c1"], "--image-size",
                           "2448x2048", "--geometry", "--output", geo)
    check("(a profile for the trial)", rc == 0, (rc, err))
    work = d / "candidate"
    args = ["trial", "--camera", "c1", "--io", IO["c1"], "--eo", EO["c1"], "--profile-json", geo,
            "--spline", spline, "--antenna-navd88", "19.014", "--work-dir", work, "--scripts-dir", stubs]
    rc, js, out, err = cli(*args, "--frames", *[src / x for x in names])
    fr = {f["frame"]: f for f in (js or {}).get("frames", [])}
    check("five frames, one passes: exit 0, passed 1 of 5, the best frame summarised",
          rc == 0 and js["ok"] is True and js["passed"] == 1 and js["tried"] == 5
          and js["best"]["frame"] == names[0] and js["best"]["points"] == 2448
          and js["best"]["placed_percent"] == 100, (rc, out[:1500], err))
    if not js:
        return
    f0 = fr[names[0]]
    check(f"the good frame: a line in 100% of the columns, {f0['placed']} of {f0['points']} points placed, "
          f"{f0['range_m']} m from the camera, its water level {f0['water_level_navd88']} m",
          f0["pass"] and f0["signal_columns_fraction"] == 1.0 and f0["placed_fraction"] == 1.0
          and 30 <= f0["range_m"][0] <= f0["range_m"][1] <= 1000 and f0["water_level_navd88"] is not None)
    check("a line in 40% of the columns: fails, said (and the contours drop it)",
          not fr[names[1]]["pass"] and any("40% of the columns" in r for r in fr[names[1]]["reasons"]),
          fr[names[1]])
    check(f"a line {near_m:.0f} m from the camera: fails the 30-1000 m rule, said",
          near_m < 30 and not fr[names[2]]["pass"] and any("outside 30-1000 m" in r for r in fr[names[2]]["reasons"]),
          fr[names[2]])
    check("no line kept by the detector: fails, said", not fr[names[3]]["pass"]
          and fr[names[3]]["detected"] is False and "kept no line" in fr[names[3]]["reasons"][0], fr[names[3]])
    check("a picture of another size than the lens file: fails, never given to the detector",
          not fr[names[4]]["pass"] and "1224 x 1024" in fr[names[4]]["reasons"][0], fr[names[4]])
    argv = json.loads((stubs / "detector_argv.json").read_text())
    opt = dict(zip(argv[::1], argv[1:] + [""]))
    t = work / "trial_c1"
    folders = [Path(opt[o]) for o in ("--source-dir", "--input-dir", "--output-dir", "--debug-dir")]
    check("the detector ran with --no-bias-correction --profile-json <the profile>, timex, and all four folders "
          "in the scratch folder", "--no-bias-correction" in argv and Path(opt["--profile-json"]) == geo.resolve()
          and opt["--image-suffix"] == "timex.jpg" and opt["--timing-mode"] == "legacy"
          and all(wsc._inside(f.resolve(), t.resolve()) for f in folders), argv)
    check("the steps' exit codes, the scratch folder and its log in the JSON",
          js["steps"] == {"detector": 0, "contours": 0, "georectify": 0} and js["work_dir"] == str(t)
          and Path(js["log"]).is_file() and "extract_elevation_contours.py" in Path(js["log"]).read_text())
    check("the kept pictures only read: unchanged, nothing added to their folder",
          {n: _sha(src / n) for n in names} == before and sorted(os.listdir(str(src))) == sorted(names))
    check("a summary for people on stderr", "1 of 5 frame(s) pass" in err and "FAIL" in err, err)

    rc, js, out, err = cli(*args, "--frames", src / names[1], src / names[3])
    check("every frame failing: exit 1, ok false, the reason", rc == 1 and js["ok"] is False and js["passed"] == 0
          and "no frame passed" in js["reason"], (rc, out[:400]))

    rc, js, out, err = cli(*args, "--frames", *[src / x for x in names], src / names[0])
    check("a frame named twice counts once", rc == 0 and js["tried"] == 5, (rc, out[:200]))
    extra = src / frame_name(T0 + 9 * 1800, "c1")
    shutil.copy2(str(src / names[0]), str(extra))
    rc, js, out, err = cli(*args, "--frames", *[src / x for x in names], extra)
    check("six frames: usage, exit 2 (5 at most)", rc == 2 and "5 at most" in err, err)
    extra.unlink()
    c2 = src / frame_name(T0, "c2")
    shutil.copy2(str(src / names[0]), str(c2))
    rc, js, out, err = cli(*args, "--frames", c2)
    check("another camera's timex: usage, exit 2", rc == 2 and "not a c1 timex" in err, err)
    c2.unlink()

    for bad in (HERE / "archive" / "trial_test_scratch", HERE / "processed_timex" / "x", HERE / "input",
                Path("/mnt/I2Rgus_Data/waterline/archive/x"), HERE / "extra_timex_source" / "x",
                Path("/mnt/I2Rgus_Data/waterline/extra_timex_backfill_debug"), Path("/mnt/I2Rgus_Data/extra_timex/x")):
        existed = bad.exists()
        rc, js, out, err = cli("trial", "--camera", "c1", "--io", IO["c1"], "--eo", EO["c1"], "--profile-json", geo,
                               "--spline", spline, "--frames", src / names[0], "--work-dir", bad,
                               "--scripts-dir", stubs)
        check(f"--work-dir {bad}: refused (the station's own folder), exit 2, nothing made",
              rc == 2 and "station's own folders" in err and bad.exists() == existed, (rc, err))

    prof_c2 = d / "profile_c2_only.json"
    prof_c2.write_text(json.dumps({"c2": strict_json(geo.read_text())["c1"]}))
    rc, js, out, err = cli("trial", "--camera", "c1", "--io", IO["c1"], "--eo", EO["c1"], "--profile-json", prof_c2,
                           "--spline", spline, "--antenna-navd88", "19.014", "--frames", src / names[0],
                           "--work-dir", work, "--scripts-dir", stubs)
    check("a profile without this camera: refused, exit 1 (CACO05's own band would be used)",
          rc == 1 and "no entry for c1" in js["reason"], (rc, out[:300]))
    other = d / "someone_elses"
    (other / "trial_c1").mkdir(parents=True)
    (other / "trial_c1" / "keep.txt").write_text("x")
    rc, js, out, err = cli(*args[:-4], "--work-dir", other, "--scripts-dir", stubs, "--frames", src / names[0])
    check("a trial_c1 folder this command did not make: not touched, refused, exit 1",
          rc == 1 and "not touching it" in js["reason"] and (other / "trial_c1" / "keep.txt").is_file(), (rc, out[:300]))
    rc, js, out, err = cli(*args, "--frames", src / names[0])
    check("run again in the same scratch folder: its own trial_c1 replaced, exit 0", rc == 0 and js["passed"] == 1,
          (rc, out[:300]))

    # a relative --work-dir (docs/NEW_STATION.md: --work-dir work/trial): every step runs inside the scratch
    # folder, so the folder must be made absolute first (else the detector reads trial_c1/<relative>/src)
    rel_args = [a if a != work else "relwork/trial" for a in args]
    rc, js, out, err = cli(*rel_args, "--frames", src / names[0], cwd=d)
    t_rel = d / "relwork" / "trial" / "trial_c1"
    check("a relative --work-dir: resolved against the current folder; the detector, contours and georectify "
          "all find their files, exit 0, the absolute folder in the JSON, no nested copy",
          rc == 0 and js["passed"] == 1 and js["steps"] == {"detector": 0, "contours": 0, "georectify": 0}
          and js["work_dir"] == str(t_rel) and Path(js["log"]).is_file() and not (t_rel / "relwork").exists(),
          (rc, out[:600], err[-400:]))

    # the antenna height: never CACO-05's 19.014 m without a word
    spline12 = d / "spline_hortho12.txt"
    write_spline(spline12, hortho=12.3)
    no_ant = [a for a in args if a not in ("--antenna-navd88", "19.014")]
    env_clean = {k: v for k, v in os.environ.items() if k != "GNSSR_ANTENNA_NAVD88_M"}
    i = no_ant.index("--spline") + 1
    no_ant[i] = spline12
    rc, js, out, err = cli(*no_ant, "--frames", src / names[0], env=env_clean)
    check("no --antenna-navd88 and no GNSSR_ANTENNA_NAVD88_M: usage, exit 2, said (CACO-05's 19.014 would be used)",
          rc == 2 and js and js["ok"] is False and "--antenna-navd88" in js["error"] and "19.014" in err, (rc, out[:300], err))
    levels = {}
    for what, extra, env in (("--antenna-navd88 12.3", ["--antenna-navd88", "12.3"], env_clean),
                             ("GNSSR_ANTENNA_NAVD88_M=12.3 in the environment", [],
                              dict(env_clean, GNSSR_ANTENNA_NAVD88_M="12.3")),
                             ("--antenna-navd88 19.014 (another height)", ["--antenna-navd88", "19.014"], env_clean)):
        rc, js, out, err = cli(*no_ant, *extra, "--frames", src / names[0], env=env)
        lev = js and js["frames"][0]["water_level_navd88"]
        levels[what] = lev
        check(f"{what}: the trial runs, the height used and where it came from in the JSON (water level {lev} m)",
              rc == 0 and js["antenna_navd88_m"] == (19.014 if "19.014" in what else 12.3) and lev is not None
              and js["antenna_from"] == ("--antenna-navd88" if extra else "GNSSR_ANTENNA_NAVD88_M in the environment")
              and "m NAVD88 from" in err, (rc, out[:400], err[-300:]))
    l12, lenv, l19 = levels.values()
    check(f"... the height reaches extract_elevation_contours.py: the same level from the option and the "
          f"environment ({l12}, {lenv}), 19.014 - 12.3 = 6.714 m higher with 19.014 ({l19})",
          None not in (l12, lenv, l19) and l12 == lenv and abs(l19 - l12 - 6.714) < 0.01, levels)
    rc, js, out, err = cli(*no_ant, "--frames", src / names[0], env=dict(env_clean, GNSSR_ANTENNA_NAVD88_M="12.3m"))
    check("GNSSR_ANTENNA_NAVD88_M not a number (and no option): usage, exit 2",
          rc == 2 and "not an antenna height" in err, (rc, err))
    rc, js, out, err = cli(*no_ant, "--antenna-navd88", "nan", "--frames", src / names[0], env=env_clean)
    check("--antenna-navd88 nan: usage, exit 2", rc == 2 and "not an antenna height" in err, (rc, err))

    if not HAVE_CV2:
        print("  SKIP  the real detector on a synthetic frame (OpenCV does not import here)")
        return
    real = d / "real_frame"
    real.mkdir()
    _jpeg(real / names[0], 2448, 2048, row=900)
    rc, js, out, err = cli("trial", "--camera", "c1", "--io", IO["c1"], "--eo", EO["c1"], "--profile-json", geo,
                           "--spline", spline, "--antenna-navd88", "19.014", "--work-dir", d / "real_work",
                           "--frames", real / names[0])
    f0 = (js or {}).get("frames", [{}])[0]
    check(f"the real detector, contours and georectify on a synthetic c1 frame (sea above row 900): passes, "
          f"{f0.get('points')} points {f0.get('range_m')} m from the camera",
          rc == 0 and js["steps"] == {"detector": 0, "contours": 0, "georectify": 0} and f0.get("pass") is True
          and f0["signal_columns_fraction"] >= 0.5, (rc, out[:800], err[-800:]))
    t = d / "real_work" / "trial_c1"
    log = (t / "trial.log").read_text() if (t / "trial.log").is_file() else ""
    check("... with the profile and the bias off (the detector's own log lines), its debug images in scratch",
          "PROFILE OVERRIDE for CACO05_C1" in log and "BIAS CORRECTION DISABLED" in log
          and any((t / "debug").iterdir()), log[-600:])


# ---------------------------------------------------------------------
# pointing
# ---------------------------------------------------------------------

POINTING_RULE = {"pass_below_deg": 0.3, "fail_above_deg": 1.0, "min_inliers": 40}


def test_pointing(d):
    print("pointing: a calibration-day picture against recent timex")
    # the limits themselves (pure Python: checked with or without OpenCV): PASS < 0.3, WARN 0.3-1.0, FAIL > 1.0
    for x, want in ((0.0, "PASS"), (0.29, "PASS"), (0.2999, "PASS"), (0.30, "WARN"), (0.6, "WARN"),
                    (1.00, "WARN"), (1.0001, "FAIL"), (1.01, "FAIL"), (5.0, "FAIL")):
        check(f"pointing_grade({x}) is {want} (PASS under 0.3 deg, WARN 0.3-1.0, FAIL over 1.0)",
              wsc.pointing_grade(x) == want, wsc.pointing_grade(x))
    if not HAVE_CV2:
        print("  SKIP  pointing (OpenCV does not import here; the command then says UNAVAILABLE)")
        rc, js, out, err = cli("pointing", "--io", IO["c2"], "--eo", EO["c2"], "--calibration-picture", __file__,
                               "--timex", __file__)
        check("without OpenCV: UNAVAILABLE, exit 1, JSON, the rule it would grade by",
              rc == 1 and js and js["verdict"] == "UNAVAILABLE" and js["rule"] == POINTING_RULE, (rc, out[:300]))
        return
    import pointing_daily as pdl                                     # synthetic scene and exact turns (OpenCV only)
    # c2's lens at half size (the same angles; faster matching)
    v = {}
    for line in IO["c2"].read_text().splitlines():
        if ":" in line and not line.lstrip().startswith("#"):
            k, _, rest = line.partition(":")
            v[k.strip()] = float(rest.split("#")[0])
    for k in ("NU", "NV", "coU", "coV", "fx", "fy"):
        v[k] = v[k] / 2.0
    io_path = d / "half_c2_IO.yaml"
    io_path.write_text("".join(f"{k}: {x!r}\n" for k, x in v.items()))
    io_, eo = load_intrinsics(io_path), load_extrinsics(EO["c2"])
    scene = pdl.textured_scene(io_, contrast=40.0, seed=1)
    cal = d / "calibration_day_c2.jpg"
    cv2.imwrite(str(cal), scene, [cv2.IMWRITE_JPEG_QUALITY, 95])
    for turn, want in ((0.1, "PASS"), (0.6, "WARN"), (1.6, "FAIL")):
        moved = pdl.rotated_view(scene, io_, eo, pdl.corrected_eo(eo, (turn, 0.0, 0.0)))
        tx = d / f"timex_turned_{turn}.jpg"
        cv2.imwrite(str(tx), moved, [cv2.IMWRITE_JPEG_QUALITY, 95])
        rc, js, out, err = cli("pointing", "--io", io_path, "--eo", EO["c2"], "--calibration-picture", cal,
                               "--timex", tx)
        rot = js and js["rotation_deg"]
        check(f"a camera turned {turn} deg: measured {rot} deg, {want}, exit {0 if want != 'FAIL' else 1}",
              js and js["verdict"] == want and abs(rot - turn) < 0.05 and rc == (1 if want == "FAIL" else 0)
              and js["frames"][0]["status"] == "matched" and js["ok"] is (want != "FAIL")
              and js["rule"] == POINTING_RULE and f"pointing: {want}" in err, (rc, out[:600], err))
    check("nothing written next to the calibration-day picture (its pointing goes to a scratch folder)",
          not cal.with_suffix(".json").exists())

    other = d / "unrelated.jpg"
    cv2.imwrite(str(other), pdl.textured_scene(io_, contrast=40.0, seed=7), [cv2.IMWRITE_JPEG_QUALITY, 95])
    rc, js, out, err = cli("pointing", "--io", io_path, "--eo", EO["c2"], "--calibration-picture", cal,
                           "--timex", other)
    check(f"an unrelated picture: UNMATCHED (best {js and js['frames'][0].get('inliers')} inliers), nothing checked, "
          f"exit 1, no traceback", rc == 1 and js and js["verdict"] == "UNMATCHED" and js["ok"] is False
          and js["rotation_deg"] is None and "nothing checked" in js["reason"] and "Traceback" not in err,
          (rc, out[:500], err[-500:]))
    blank = d / "fog.jpg"
    cv2.imwrite(str(blank), np.full((int(io_[1]), int(io_[0]), 3), 128, np.uint8))
    rc, js, out, err = cli("pointing", "--io", io_path, "--eo", EO["c2"], "--calibration-picture", cal,
                           "--timex", blank)
    check("a featureless (fog) frame: UNMATCHED, exit 1, no traceback",
          rc == 1 and js and js["verdict"] == "UNMATCHED" and "Traceback" not in err, (rc, out[:300], err[-300:]))
    small = d / "small.jpg"
    cv2.imwrite(str(small), cv2.resize(scene, (320, 256)))
    rc, js, out, err = cli("pointing", "--io", io_path, "--eo", EO["c2"], "--calibration-picture", cal,
                           "--timex", small)
    check("a timex of another size: ERROR (nothing compared), exit 1",
          rc == 1 and js["verdict"] == "ERROR" and "320 x 256" in js["reason"], (rc, out[:300]))
    rc, js, out, err = cli("pointing", "--io", IO["c2"], "--eo", EO["c2"], "--calibration-picture", cal,
                           "--timex", other)
    check("a calibration-day picture of another size than the lens file: ERROR, exit 1",
          rc == 1 and js["verdict"] == "ERROR" and "lens file 2448 x 2048" in js["reason"], (rc, out[:300]))
    rc, js, out, err = cli("pointing", "--io", io_path, "--eo", EO["c2"], "--calibration-picture", d / "none.jpg",
                           "--timex", other)
    check("a picture that is not there: usage, exit 2", rc == 2 and "not found" in err, err)


def test_cli():
    print("the command line")
    rc, js, out, err = cli("capabilities")
    check("capabilities unchanged: 'overlay-1 storm-alerts-1 station-setup-1'",
          rc == 0 and out == "overlay-1 storm-alerts-1 station-setup-1\n", out)
    for sub in ("pair", "profile", "trial", "pointing"):
        rc, js, out, err = cli(sub)
        check(f"{sub} with no options: usage, exit 2", rc == 2 and "required" in err, err[-200:])
    import ast
    top = set()
    for node in ast.parse(CHECK.read_text()).body:                 # the module's own top-level imports
        if isinstance(node, ast.Import):
            top.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            top.add((node.module or "").split(".")[0])
    std = {"argparse", "csv", "json", "math", "os", "re", "subprocess", "sys", "datetime", "pathlib"}
    check(f"capabilities and check-overlay stay standard library: the top-level imports are {sorted(top)}",
          top <= std, sorted(top - std))
    for f in ("waterline_station_check.py", "test_waterline_station_check.py"):
        try:
            ast.parse((HERE / f).read_text(), filename=f, feature_version=(3, 8))
            ok = True
        except SyntaxError as e:
            ok = e
        check(f"{f} parses as Python 3.8", ok is True, ok)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="keep the test folder, made inside this one")
    args = ap.parse_args()
    print(f"Python {sys.version.split()[0]}, numpy {np.__version__}, OpenCV {'yes' if HAVE_CV2 else 'no'}")
    d = Path(tempfile.mkdtemp(prefix="test_wsc_", dir=args.keep))
    try:
        test_cli()
        test_pair(d)
        test_profile(d)
        test_trial(d)
        test_pointing(d)
    finally:
        if args.keep:
            print(f"kept: {d}")
        else:
            shutil.rmtree(str(d), ignore_errors=True)
    print()
    if FAILED:
        print(f"FAILED: {len(FAILED)} check(s)")
        for m in FAILED:
            print("  " + m)
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())

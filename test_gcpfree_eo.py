#!/usr/bin/env python3
"""
Tests For gcpfree_eo.py (calibration B, 23 Oct 2024)
======================================================
Plain python3 (no pytest needed on the station):

    python3 test_gcpfree_eo.py [--keep DIR]

WHY. Calibration B goes to a blind check: it must never read calibration A
(CACO03_<cam>_20241023_EO.yaml), a GCP or target file or the calibration
history; it must recover a known turn when the photos allow it and say NOT
DELIVERED (no B file) when they do not, never falling back to anything; its
outcome rule must be the owner's, pre-registered, committed and pushed before any
run; and nothing it delivers may be chosen after seeing a result (run
folders, no overwrite, one more run only with Aug photos).

HOW. A synthetic world that depends only on the viewing DIRECTION (a sky
above the true sea horizon, a dark sea, textured land below 4 deg of
depression), so any turn of the camera is rendered exactly through the lens
model of CACO03_<cam>_20240801_IO.yaml. Jan photos are rendered with
CACO03_<cam>_20250123_EO.yaml, Oct photos with it turned by a known
rotation (the brief's test sizes), Dec photos (c2) half way, Aug photos with
CACO03_<cam>_20240827_EO.yaml. A decoy calibration A, a decoy calibration
history and a decoy GCP file sit in the test's calibration folder: an audit
hook (sys.addaudithook: every open made through Python, by any module:
builtins, io, pathlib, os, numpy, PIL) and a watch on cv2.imread (OpenCV
reads files in C++, which no audit hook sees) record any read of a file named
like them, and one fails the test. The run logic
(pre-registration, run folders, outcome rules, the second run) is tested
with the real procedure's results of the rendered sets, replayed.

STATION FILES. None: the photos, the calibration folder, the
pre-registration (a git repository of its own) and the outputs are in a
temporary folder; the lens, Jan and Aug EO files are copied there.
"""

import os
import re
import sys
import copy
import json
import time
import shutil
import argparse
import tempfile
import subprocess
from pathlib import Path
from types import SimpleNamespace
from datetime import datetime, timezone

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gcpfree_eo as G                                        # noqa: E402
import estimate_eo_rotation as EER                            # noqa: E402
from georectify import load_intrinsics, load_extrinsics, cirn_angles_to_R, undistort_uv  # noqa: E402

TURN = {"c1": (-12.5, 0.6, -0.4), "c2": (-13.4, -0.5, 0.3)}  # test sizes only (deg: pan, tilt, roll)
N_OK = [0]

# ---- 'never opened': an audit hook sees every open made through Python, whatever module makes it (builtins, io,
# pathlib, os, numpy, PIL); OpenCV reads its files in C++, unseen by it, so cv2.imread is watched as well (watch_cv2)
HELD = re.compile(r"(20241023_EO|calibration_history|Extrinsic_Targets|2024-?10-?23[^/]*(xyz|gcp|target|survey))",
                  re.I)
WATCH = {"on": False, "opened": []}


def _audit(event, args):
    if not WATCH["on"] or event != "open" or not args or args[0] is None or isinstance(args[0], int):
        return
    try:
        p = os.fsdecode(args[0])
    except TypeError:
        return
    for n in (os.path.basename(p), os.path.basename(os.path.realpath(p))):
        if HELD.search(n):
            WATCH["opened"].append(p)
            return


if hasattr(sys, "addaudithook"):
    sys.addaudithook(_audit)


def watch_cv2():
    """cv2.imread through _audit too (the procedure reads every photo with it). -> a function that undoes it."""
    import cv2
    real = cv2.imread

    def imread(path, *a, **kw):
        _audit("open", (path, "r", None))
        return real(path, *a, **kw)
    cv2.imread = imread

    def undo():
        cv2.imread = real
    return undo


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)
    N_OK[0] += 1
    print(f"  ok  {msg}")


def turned(eo, d):
    e = np.array(eo, float).copy()
    e[3:6] += np.radians(d)
    return e


class World:
    """Brightness as a function of the world viewing direction only (exact under any rotation)."""

    def __init__(self, seed, h):
        import cv2
        rng = np.random.default_rng(seed)
        self.res = 0.04                                         # deg per texel
        n_az, n_el = int(360 / self.res), int(40 / self.res)    # elevation -40..0 deg
        tex = np.zeros((n_el, n_az), np.float32)
        for s, a in ((3, 1.0), (8, 0.8), (20, 0.6)):
            t = cv2.GaussianBlur(rng.standard_normal((n_el, n_az)).astype(np.float32), (0, 0), s)
            tex += a * t / (t.std() + 1e-6)
        self.land = tex
        self.dip = np.degrees(EER.horizon_dip(h))
        self._rays, self._dirs = None, {}

    def directions(self, io, eo):
        """(azimuth, elevation) in deg of every pixel's ray; the lens is undistorted once per world."""
        key = tuple(np.round(eo, 12))
        if key not in self._dirs:
            nu, nv = int(io[0]), int(io[1])
            if self._rays is None:
                U, V = np.meshgrid(np.arange(nu, dtype=float), np.arange(nv, dtype=float))
                x, y = undistort_uv(U.ravel(), V.ravel(), io)
                K = np.array([[-io[4], 0, io[2]], [0, -io[5], io[3]], [0, 0, 1.0]])
                self._rays = np.linalg.solve(K, np.vstack([x, y, np.ones_like(x)]))
            d = cirn_angles_to_R(*eo[3:6]).T @ self._rays
            d /= np.linalg.norm(d, axis=0)
            az = (np.degrees(np.arctan2(d[0], d[1])) % 360.0).reshape(nv, nu).astype(np.float32)
            el = np.degrees(np.arcsin(np.clip(d[2], -1, 1))).reshape(nv, nu).astype(np.float32)
            if len(self._dirs) > 6:
                self._dirs.clear()
            self._dirs[key] = (az, el)
        return self._dirs[key]

    def render(self, io, eo, day_seed, beach_change=0.0, gain=1.0, fog=False):
        import cv2
        az, el = self.directions(io, eo)
        rng = np.random.default_rng(day_seed)
        if fog:                                                 # nothing to see: no horizon, no texture
            img = np.full(az.shape, 150.0) + rng.normal(0, 2.0, az.shape)
            return np.clip(img, 0, 255).astype(np.uint8)
        land = cv2.remap(self.land, az / self.res, -el / self.res, cv2.INTER_LINEAR, borderMode=cv2.BORDER_WRAP)
        if beach_change:                                        # a beach that changed: its own texture
            ch = cv2.GaussianBlur(rng.standard_normal(self.land.shape).astype(np.float32), (0, 0), 6)
            ch = cv2.remap(ch / (ch.std() + 1e-6), az / self.res, -el / self.res, cv2.INTER_LINEAR,
                           borderMode=cv2.BORDER_WRAP)
            near = np.clip((-el - 4.0) / 6.0, 0, 1)             # the change on the nearer beach only
            land = land * (1 - beach_change * near) + ch * beach_change * near
        img = np.where(el > -self.dip, 205.0 - 0.8 * el, 70.0 + 3.0 * np.clip(-el, 0, 4))   # sky / sea
        img = np.where(el < -4.0, 120.0 + 28.0 * land, img)                                     # land
        img = img * gain + rng.normal(0, 2.0, img.shape)
        img = cv2.GaussianBlur(img.astype(np.float32), (0, 0), 0.8)
        return np.clip(img, 0, 255).astype(np.uint8)


def name(epoch, cam, station="CACO03"):
    t = datetime.fromtimestamp(epoch, tz=timezone.utc)
    return f"{epoch}.{t.strftime('%a.%b.%d_%H_%M_%S')}.GMT.{t.year}.{station}.{cam}.timex.jpg"


def epoch(day, hour, minute=0):
    y, m, d = (int(v) for v in day.split("-"))
    return int(datetime(y, m, d, hour, minute, tzinfo=timezone.utc).timestamp())


def write(path, img):
    import cv2
    path.parent.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(path), img, [cv2.IMWRITE_JPEG_QUALITY, 92])


def make_photos(root, cam, io, eo_jan, world, oct_change=0.5, oct_world=None, dec=True, jan_eo=None):
    """Jan 22-23 Jan 2025 (jan_eo, default the Jan EO), Oct 21/23/25 Oct 2024 (turned by TURN), c2's Dec
    10/12/14 Dec 2024 (half the turn), at 14, 17 and 20 UTC; distractors: another station, a half-size copy."""
    import cv2
    jan_eo = eo_jan if jan_eo is None else jan_eo
    eo_oct = turned(jan_eo, TURN[cam])
    eo_dec = turned(jan_eo, tuple(0.5 * v for v in TURN[cam]))
    plan = [(d, jan_eo, 0.0) for d in ("2025-01-22", "2025-01-23")]
    plan += [(d, eo_oct, oct_change) for d in ("2024-10-21", "2024-10-23", "2024-10-25")]
    if cam == "c2" and dec:
        plan += [(d, eo_dec, 0.3) for d in ("2024-12-10", "2024-12-12", "2024-12-14")]
    k, img = 0, None
    for day, eo, ch in plan:
        w = oct_world if (oct_world is not None and day.startswith("2024-10")) else world
        for h in (14, 17, 20):
            k += 1
            img = w.render(io, eo, day_seed=1000 + k, beach_change=ch, gain=0.9 + 0.05 * (k % 4))
            write(root / "pics" / name(epoch(day, h), cam), img)
    write(root / "pics" / name(epoch("2024-10-23", 17, 30), cam, "CACO04"), img)      # another station
    write(root / "owg" / name(epoch("2024-10-23", 17), cam), cv2.resize(img, (1224, 1024)))  # reduced copy
    return eo_oct, eo_dec


def make_aug(folder, cam, io, eo_aug, world, seed=500):
    for day in ("2024-08-25", "2024-08-27", "2024-08-29"):
        for h in (14, 17, 20):
            seed += 1
            write(folder / name(epoch(day, h), cam), world.render(io, eo_aug, seed, 0.4, 0.95))


def roots_of(root):
    return [str(root / "pics"), str(root / "owg")]


def cal_eos(cam):
    io = load_intrinsics(G.CAL / G.IO_FILE.format(cam=cam))
    return (io, load_extrinsics(G.CAL / G.REF_EO_FILE["jan"].format(cam=cam)),
            load_extrinsics(G.CAL / G.REF_EO_FILE["aug"].format(cam=cam)))


# ------------------------------------------------------------------ pure parts

def test_watch_sees_reads(d):
    print("the 'never opened' watch itself: it sees a Python open and an OpenCV read of a decoy")
    import cv2
    decoy = G.CAL / "CACO03_c1_20241023_EO.yaml"
    WATCH["on"] = True
    try:
        with open(str(decoy)) as f:
            f.read()
        python_seen = list(WATCH["opened"])
        del WATCH["opened"][:]
        cv2.imread(str(decoy), cv2.IMREAD_GRAYSCALE)          # not an image: None, but the read was asked for
        cv2_seen = list(WATCH["opened"])
        del WATCH["opened"][:]
    finally:
        WATCH["on"] = False
    if hasattr(sys, "addaudithook"):
        check(python_seen == [str(decoy)], "a Python open() of the decoy calibration A is seen (audit hook)")
    check(cv2_seen == [str(decoy)], "a cv2.imread of the decoy calibration A is seen (the cv2 watch)")


def test_guard(d):
    print("guard: never calibration A, a GCP or target file, the calibration history; by name and by link")
    for bad in ("CACO03_c1_20241023_EO.yaml", "2024-10-23_GCPs.csv", "calibration_history.csv",
                "2025-11-13_Marconi_Extrinsic_Targets_c1_xyz.csv", "2024-10-23_Marconi_points_xyz.csv"):
        try:
            G.guard(d / bad)
            check(False, f"guard refuses {bad}")
        except G.Refused:
            check(True, f"guard refuses {bad} (Refused: exit 2)")
    for ok in ("CACO03_c1_20250123_EO.yaml", "gcpfree_eo.py", "CACO03_c1_gcpfree-oct2024_EO.yaml",
               "CACO03_c1_gcpfree-oct2024_report.json", "PREREGISTERED.txt"):
        check(G.guard(d / ok).name == ok, f"guard lets {ok} through")
    link = d / "links" / "CACO03_c2_20250123_EO.yaml"
    link.parent.mkdir()
    (d / "links" / "CACO03_c2_20241023_EO.yaml").write_text("x: 0\n# DECOY\n")
    os.symlink(str(d / "links" / "CACO03_c2_20241023_EO.yaml"), str(link))
    try:
        G.guard(link)
        check(False, "a link named like the Jan EO that points to calibration A: refused")
    except G.Refused as exc:
        check("a link to CACO03_c2_20241023_EO.yaml" in str(exc),
              "a link named like the Jan EO that points to calibration A: refused (os.path.realpath)")


def test_scan_never_opens(d):
    print("the photo scan: a photo-named link to a forbidden file is never opened (not even for its size)")
    import detect_original_view as DOV
    io, _, _ = cal_eos("c1")
    root = d / "scan" / "pics"
    root.mkdir(parents=True)
    target = G.CAL / "2024-10-23_Marconi_Extrinsic_Targets_c1_xyz.csv"          # the test's decoy GCP file
    link = root / name(epoch("2024-10-23", 16), "c1")
    os.symlink(str(target), str(link))
    WATCH["on"] = True
    try:
        DOV.jpeg_size(link)                                  # the unguarded probe: the watch sees it
        probe_seen = list(WATCH["opened"])
        del WATCH["opened"][:]
        days, counts = G.select_photos([str(root)], "c1", G.WINDOWS["oct"], io)
        scan_seen = list(WATCH["opened"])
        del WATCH["opened"][:]
    finally:
        WATCH["on"] = False
    if hasattr(sys, "addaudithook"):
        check(probe_seen, "an unguarded size probe of the link opens the decoy (the watch sees it)")
        check(not scan_seen, "gcpfree_eo's scan never opens it")
    check(not days and counts["forbidden_not_opened"] == [link.name] and DOV.jpeg_size is not None
          and DOV.jpeg_size.__name__ == "jpeg_size",
          f"... it is left out and counted as not opened ({link.name}); the module's own probe is restored")


def test_archived(d):
    print("the calibration files read are the archived ones (sha256), else refused before anything is solved")
    for f in G.ARCHIVED_SHA256:
        check(G.archived(G.CAL / f).name == f, f"{f}: its archived content")
    other = d / "edited"
    other.mkdir()
    (other / "CACO03_c1_20250123_EO.yaml").write_text((G.CAL / "CACO03_c1_20250123_EO.yaml").read_text() + "# edit\n")
    (other / "CACO03_c1_20240101_EO.yaml").write_text("x: 0\n")
    for p, words in ((other / "CACO03_c1_20250123_EO.yaml", "is not the archived"),
                     (other / "CACO03_c1_20240101_EO.yaml", "not one of the calibration files"),
                     (other / "CACO03_c2_20250123_EO.yaml", "not found")):
        try:
            G.archived(p)
            check(False, f"{p.name}: refused")
        except G.Refused as exc:
            check(words in str(exc), f"{p.name} ({words}): refused before anything is solved")


def test_ransac():
    print("2-ray RANSAC")
    rng = np.random.default_rng(1)
    Q = cirn_angles_to_R(*np.radians([60.0, 72.0, -7.0])) @ cirn_angles_to_R(*np.radians([72.57, 71.6, -6.67])).T
    a = rng.normal(size=(400, 3))
    a /= np.linalg.norm(a, axis=1)[:, None]
    b = a @ Q.T
    b[40:] = rng.normal(size=(360, 3))
    b /= np.linalg.norm(b, axis=1)[:, None]                    # 90 % outliers
    Qe, inl = G.ransac2(a, b, 3.0 / 2575.0)
    # 1e-4 deg: the arccos of a trace near 3 is rounding noise below ~1e-5 deg, and that noise differs with the
    # LAPACK under numpy (0 on numpy 1.17, 1.7e-6 deg on numpy 2.4)
    check(Qe is not None and G.dist_deg(Qe, Q) < 1e-4 and inl[:40].all() and inl.sum() < 45,
          f"2-ray RANSAC finds a 13 deg turn among 90 % outliers ({int(inl.sum())} inliers)")


def link(gate_F=True, gate_S=True, hz="pass", px=1.0):
    return {"gate_F": gate_F, "gate_S": gate_S, "horizon": {"status": hz, "why": f"horizon {hz}"},
            "solve": {"jackknife_max_px": px, "why": "" if gate_F else "60 inliers < 60"},
            "per_day": {"days": [{"day": "2024-10-21", "diff_px": px}], "why": "" if gate_S else "2 days measured"}}


def test_verdict_rules():
    print("the verdict rules: the code does what the docstring says, also when the horizon is unavailable")
    doc = " ".join(G.__doc__.split())
    for words in ("A ROUTE counts when its solve passes gates F and S and its horizon check is not FAIL (pass or "
                  "unavailable; for the chain: both links)",
                  "DELIVERED-WIDER the route used counts but its horizon is unavailable (tilt and roll unchecked)",
                  "B = c2's chain (two links: never better than WIDER)", "NOT DELIVERED no route counts"):
        check(words in doc, f"docstring: '{words[:60]}...'")
    cases = [
        ("direct, horizon pass", {"direct": dict(link(), has_B=True)}, ("DELIVERED", "direct Jan->Oct")),
        ("direct, horizon unavailable", {"direct": dict(link(hz="unavailable"), has_B=True)},
         ("DELIVERED-WIDER", "direct Jan->Oct")),
        ("direct, horizon FAIL, nothing else", {"direct": dict(link(hz="FAIL"), has_B=True)}, ("NOT DELIVERED", None)),
        ("direct horizon FAIL, Aug passes", {"direct": dict(link(hz="FAIL"), has_B=True),
                                             "aug": {"link": link(), "has_B": True}}, ("DELIVERED", "Aug->Oct")),
        ("direct fails F, Aug horizon unavailable", {"direct": dict(link(gate_F=False), has_B=True),
                                                     "aug": {"link": link(hz="unavailable"), "has_B": True}},
         ("DELIVERED-WIDER", "Aug->Oct")),
        ("direct fails S, chain counts (link 2 horizon unavailable)",
         {"direct": dict(link(gate_S=False), has_B=True),
          "chain": {"link1": link(), "link2": link(hz="unavailable"), "has_B": True}},
         ("DELIVERED-WIDER", "chain Jan->Dec->Oct")),
        ("direct fails, chain link 2 horizon FAIL", {"direct": dict(link(gate_F=False), has_B=True),
                                                     "chain": {"link1": link(), "link2": link(hz="FAIL"),
                                                               "has_B": True}}, ("NOT DELIVERED", None)),
        ("direct passes, chain differs by 12 px", {"direct": dict(link(), has_B=True),
                                                   "chain": {"link1": link(), "link2": link(), "has_B": True,
                                                             "vs_direct_px": 12.0, "vs_direct_deg": 0.3}},
         ("DELIVERED-WIDER", "direct Jan->Oct")),
        ("direct passes, B_aug 3 m off on the beach", {"direct": dict(link(), has_B=True),
                                                       "aug": {"link": link(), "has_B": True,
                                                               "vs_direct_ground_m": {"median": 3.0}}},
         ("DELIVERED-WIDER", "direct Jan->Oct")),
        ("no B at all", {"direct": {"has_B": False, "gate_F": False, "why": "no usable photos on one side"}},
         ("NOT DELIVERED", None)),
    ]
    for what, res, want in cases:
        v, src, reasons, notes = G.decide(res)
        check((v, src) == want, f"{what}: {v} ({src})" + (f" -- {reasons[0][:70]}" if reasons else ""))
    v, _, reasons, _ = G.decide({"direct": dict(link(hz="unavailable"), has_B=True)})
    check(any("horizon check unavailable" in r for r in reasons), "... WIDER says the horizon was unavailable")


def test_uncertainty():
    print("the stated uncertainty: two terms, neither a validated bound on real photos")
    io = load_intrinsics(G.CAL / G.IO_FILE.format(cam="c2"))
    res = {"camera": "c2", "source": "direct Jan->Oct", "direct": link(px=1.2)}
    u = G.stated_uncertainty(res, io)
    check(u["rotation_px"] == G.ROT_FLOOR_PX and "floor 4 px" in u["rotation_rule"],
          f"rotation: max(2 x spread 1.2 px, floor 4 px) = {u['rotation_px']} px")
    res = {"camera": "c2", "source": "direct Jan->Oct", "direct": link(px=3.0)}
    check(G.stated_uncertainty(res, io)["rotation_px"] == 6.0, "rotation: 2 x a spread of 3 px = 6 px")
    res = {"camera": "c2", "source": "chain Jan->Dec->Oct", "chain": {"link1": link(px=1.5), "link2": link(px=1.5)}}
    check(G.stated_uncertainty(res, io)["rotation_px"] == 6.0, "the chain: the sum of its links' spreads, doubled")
    check(u["position"]["dem_abs_median_m"] == 0.018 and u["position"]["dem_p90_m"] == 0.086
          and "1.23 m" in u["position"]["reference_solution_spread"]
          and G.POSITION_SPREAD["c1"].startswith("0.36 m"),
          "position: c2 DEM |dz| median 0.018 m, p90 0.086 m; the Aug-Jan solutions 0.36 m (c1) / 1.23 m (c2) apart")
    txt = G.uncertainty_text(u, "c2")
    check("neither term is a validated bound on real photos" in txt and "inherited camera position" in txt
          and "in the largest of five synthetic position-inheritance cases" in txt and "<=" not in txt
          and "3-7 mm" not in txt and "real-texture" not in txt and "0.001-0.012" not in txt
          and "rotation_synthetic" not in u,
          "the text: both terms ('neither term is a validated bound on real photos'; the position term the largest of "
          "five cases, not a bound); no synthetic DEM figure for the rotation (those errors were roll), no "
          "real-texture figure")


def test_outcome_rules():
    print("the outcome rule (the owner's, pre-registered): both-or-none | per-camera")
    base = {"c1": {"verdict": "NOT DELIVERED"}, "c2": {"verdict": "DELIVERED"}}
    r = G.apply_outcome(copy.deepcopy(base), "both-or-none", True)
    check(r["c1"]["delivery"] == "NOT DELIVERED" and r["c2"]["delivery"].startswith("WITHHELD"),
          "both-or-none, only c2 passed: c2 WITHHELD, c1 NOT DELIVERED (no B)")
    r = G.apply_outcome(copy.deepcopy(base), "per-camera", True)
    check(r["c2"]["delivery"] == "DELIVERED" and r["c1"]["delivery"] == "NOT DELIVERED",
          "per-camera, only c2 passed: c2 DELIVERED, c1 NOT DELIVERED")
    r = G.apply_outcome({"c1": {"verdict": "DELIVERED-WIDER"}, "c2": {"verdict": "DELIVERED"}}, "both-or-none", True)
    check(all(x["delivery"] == "DELIVERED" for x in r.values()), "both pass: both delivered under either rule")
    r = G.apply_outcome(copy.deepcopy(base), "per-camera", False)
    check(all(x["delivery"].startswith("TRIAL") for x in r.values()), "a TRIAL run delivers nothing")


# ------------------------------------------------------------------ rendered sets

def test_horizon_band(d):
    print("the horizon check in a band: a perfect B carries the reference EO's own horizon error, turned")
    cam = "c2"
    io, eo_jan, _ = cal_eos(cam)
    truth_jan = turned(eo_jan, (0.0, 3.0, 0.0))        # the 'Jan solution' 3 deg of tilt off the world's truth
    truth_oct = turned(truth_jan, TURN[cam])
    world = World(31, eo_jan[2])
    ref, new = [], []
    for k, h in enumerate((14, 17, 20)):
        p = d / "hz" / name(epoch("2025-01-23", h), cam)
        write(p, world.render(io, truth_jan, 300 + k))
        ref.append(SimpleNamespace(path=p))
        p = d / "hz" / name(epoch("2024-10-23", h), cam)
        write(p, world.render(io, truth_oct, 400 + k))
        new.append(SimpleNamespace(path=p))
    R = cirn_angles_to_R
    Q = R(*truth_oct[3:6]) @ R(*truth_jan[3:6]).T        # what the photos measure
    B = G.eo_from(Q, eo_jan)                              # a perfect B: the measured turn on the Jan solution
    h = G.horizon_check(io, eo_jan, ref, B, new, Q)
    check(h["status"] == "pass" and abs(h["diff_tilt_deg"] - h["bias_tilt_deg"]) < 0.05
          and abs(h["diff_roll_deg"] - h["bias_roll_deg"]) < 0.05,
          f"perfect B: differential {h['diff_tilt_deg']:+.3f} / {h['diff_roll_deg']:+.3f} deg = the expected bias "
          f"{h['bias_tilt_deg']:+.3f} / {h['bias_roll_deg']:+.3f} (the Jan photos turned by Q): PASS in the band")
    old = abs(h["diff_tilt_deg"]) <= G.HZ_TILT_TOL and abs(h["diff_roll_deg"]) <= G.HZ_ROLL_TOL
    check(not old, f"... which a gate around 0 (+/-{G.HZ_TILT_TOL} / {G.HZ_ROLL_TOL} deg) would have FAILED")
    check(h["band_tilt_deg"][0] <= min(0, h["bias_tilt_deg"]) - G.HZ_TILT_TOL + 1e-9
          and h["band_roll_deg"][1] >= max(0, h["bias_roll_deg"]) + G.HZ_ROLL_TOL - 1e-9
          and h.get("prediction") is not None,
          f"band [min(0, bias) - tol, max(0, bias) + tol] reported (tilt {h['band_tilt_deg']}, roll "
          f"{h['band_roll_deg']}), with the geometric prediction")
    for dd, what in (((0.0, 0.0, -0.8), "roll -0.8 deg"), ((0.0, 0.6, 0.0), "tilt +0.6 deg")):
        Bw = turned(B, dd)
        Qw = R(*Bw[3:6]) @ R(*eo_jan[3:6]).T
        hw = G.horizon_check(io, eo_jan, ref, Bw, new, Qw)
        check(hw["status"] == "FAIL", f"a B off by {what}: FAIL ({hw.get('why', '')[:90]})")
    fog = []
    for k, hh in enumerate((14, 17, 20)):
        p = d / "hz" / ("fog_" + name(epoch("2024-10-24", hh), cam))
        write(p, world.render(io, truth_oct, 500 + k, fog=True))
        fog.append(SimpleNamespace(path=p))
    hf = G.horizon_check(io, eo_jan, ref, B, fog, Q)
    check(hf["status"] == "unavailable", f"fog: no clear horizon -> 'unavailable' ({hf.get('why')})")


def run_cam(cam, roots, aug_dir=None):
    lines = []
    res, B = G.run_camera(cam, roots, aug_dir, lines.append)
    print("\n".join("      | " + x for x in lines))
    return res, B


def test_c2_direct_and_chain(d):
    print("c2: Jan -> Oct (and the chain through Dec): DELIVERED, the true pointing, calibration A never opened")
    io, eo_jan, _ = cal_eos("c2")
    eo_oct, _ = make_photos(d / "roots", "c2", io, eo_jan, World(7, eo_jan[2]))
    WATCH["on"] = True
    try:
        res, B = run_cam("c2", roots_of(d / "roots"))
    finally:
        WATCH["on"] = False
    check(not WATCH["opened"], "calibration A, the history and the GCP decoys were never opened (audit hook: "
                               + (", ".join(WATCH["opened"]) or "none") + ")")
    sel = res["photos_selected"]
    check(sum(len(v) for v in sel["oct"].values()) == 9 and not any("CACO04" in n for v in sel["oct"].values() for n in v),
          "9 Oct photos, the other station left out, the half-size copy never chosen")
    check(res["direct"]["gate_F"] and res["direct"]["gate_S"], "direct solve passes gates F and S")
    err = G.dist_deg(cirn_angles_to_R(*B[3:6]), cirn_angles_to_R(*eo_oct[3:6]))
    check(err < 0.01, f"B within 0.01 deg of the true Oct pointing ({err:.4f} deg)")
    check(np.allclose(B[:3], eo_jan[:3]), "position kept from the Jan EO")
    h = res["direct"]["horizon"]
    check(h["status"] == "pass" and "band_tilt_deg" in h, f"horizon in the band: {G.horizon_text(h)[:100]}")
    ch = res["chain"]
    check(G.link_counts(ch["link1"]) and G.link_counts(ch["link2"]) and ch["vs_direct_deg"] < 0.02,
          f"c2 chain Jan->Dec->Oct counts and agrees with direct ({ch['vs_direct_deg']} deg)")
    check(res["verdict"] == "DELIVERED" and res["source"] == "direct Jan->Oct", f"verdict DELIVERED ({res['reasons']})")
    u = res["uncertainty"]
    check(u["rotation_px"] >= G.ROT_FLOOR_PX and "position" in u, f"stated uncertainty: {u['rotation_px']} px + position")
    return res


def test_c1_with_aug(d):
    print("c1 with Aug photos (--aug-dir): Aug is a second reference, the Jan-based B stays primary")
    io, eo_jan, eo_aug = cal_eos("c1")
    world = World(11, eo_jan[2])
    make_photos(d / "roots", "c1", io, eo_jan, world)
    make_aug(d / "aug_c1", "c1", io, eo_aug, world)
    WATCH["on"] = True
    try:
        res, B = run_cam("c1", roots_of(d / "roots"), str(d / "aug_c1"))
    finally:
        WATCH["on"] = False
    check(not WATCH["opened"], "nothing held opened")
    check(sum(len(v) for v in res["photos_selected"]["aug"].values()) == 9, "9 Aug photos taken from --aug-dir")
    a = res["aug"]
    check(G.link_counts(a["link"]) and a["reference_eo"] == "CACO03_c1_20240827_EO.yaml", "Aug->Oct counts (Aug EO)")
    check(res["verdict"] == "DELIVERED" and res["source"] == "direct Jan->Oct" and "vs_direct_ground_m" in a,
          f"the Jan-based B stays primary, B_aug compared on the beach ({a.get('vs_direct_ground_m')} m)")
    return res


def test_c1_aug_only(d):
    print("c1: Jan shares nothing with Oct, Aug does: B = the Aug 2024 solve, its position the Aug solution's")
    io, eo_jan, eo_aug = cal_eos("c1")
    autumn = World(23, eo_jan[2])
    root = d / "roots_augonly"
    eo_oct, _ = make_photos(root, "c1", io, eo_jan, World(5, eo_jan[2]), oct_change=0.0, oct_world=autumn)
    make_aug(d / "aug_augonly", "c1", io, eo_aug, autumn, seed=700)
    res, B = run_cam("c1", roots_of(root), str(d / "aug_augonly"))
    check(not (res["direct"].get("gate_F") and res["direct"].get("gate_S")), "direct Jan->Oct fails (nothing in common)")
    check(res["verdict"] == "DELIVERED" and res["source"] == "Aug->Oct",
          f"B = the Aug 2024 solve, DELIVERED ({res['verdict']}, {res['source']})")
    err = G.dist_deg(cirn_angles_to_R(*B[3:6]), cirn_angles_to_R(*eo_oct[3:6]))
    check(err < 0.01 and np.allclose(B[:3], eo_aug[:3]) and res["b_reference_eo"] == "CACO03_c1_20240827_EO.yaml",
          f"B_aug on the true Oct pointing ({err:.4f} deg), the Aug position, the Aug EO named as its reference")
    return res


def test_not_delivered(d):
    print("c2: Oct photos of ANOTHER world (nothing in common): NOT DELIVERED, never a fallback")
    io, eo_jan, _ = cal_eos("c2")
    root = d / "roots_other"
    make_photos(root, "c2", io, eo_jan, World(7, eo_jan[2]), oct_change=0.0, oct_world=World(99, 20.0), dec=False)
    res, B = run_cam("c2", roots_of(root))
    check(res["verdict"] == "NOT DELIVERED" and B is None and "B" not in res,
          f"no common features -> NOT DELIVERED, no B ({res['reasons'][0][:90]})")
    return res


# ------------------------------------------------------------------ runs: pre-registration, history, outcome

def git(repo, *a):
    return subprocess.run(["git", "-C", str(repo), "-c", "user.name=test", "-c", "user.email=test@example.org"]
                          + list(a), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)


def push(repo):
    """The owner's `git push` after a commit: to a bare repository of the test standing in for GitHub (made at the
    first push), which puts the commit on a remote-tracking branch (origin/<branch>), as on the station."""
    bare = Path(str(repo) + "_github.git")
    if not bare.exists():
        subprocess.run(["git", "init", "-q", "--bare", str(bare)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        git(repo, "remote", "add", "origin", str(bare))
    return git(repo, "-c", "push.negotiate=false", "push", "-q", "origin", "HEAD")


def main_q(argv):
    """gcpfree_eo.main, its output captured. -> (exit status, output)"""
    import io as io_
    import contextlib
    buf = io_.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = G.main(argv)
    return rc, buf.getvalue()


def nd_of(res):
    """The same camera with nothing that counts: its direct link failed (as the rendered c1 aug-only set)."""
    r = copy.deepcopy(res)
    for k in ("aug", "B", "uncertainty", "pan_change_deg", "b_reference_eo", "b_reference_eo_sha256", "chain"):
        r.pop(k, None)
    r["verdict"], r["source"], r["reasons"], r["notes"] = G.decide(r)
    return r


def test_pinned_version(d):
    print("only the version survey_products.py pins is pre-registered or run (the one delivery run is never spent "
          "on a B that would not count)")
    pins, why = G.pinned_versions()
    check(pins and pins[-1]["gcpfree_eo.py"] == G.own_sha256() and pins[-1]["modules"] == G.module_sha256()
          and pins[-1]["constants_sha256"] == G.constants_sha256() and G.pin_problem() == ""
          and G.constants_problem() == "",
          f"survey_products.py's GCPFREE_PINS[-1] is this gcpfree_eo.py ({G.own_sha256()[:12]}...), its modules and "
          f"its pre-registered constants")
    real_roots = list(G.ROOTS)
    G.ROOTS = ["/elsewhere"]
    try:
        cp = G.constants_problem()
    finally:
        G.ROOTS = real_roots
    check("not the pinned ones" in cp and G.constants_problem() == "",
          "a constant changed in the running process (ROOTS): its constants are not the pinned ones (recorded in "
          "run.json; survey_products.py calls such a B not the pre-registered procedure)")
    other = d / "pins_other" / "survey_products.py"
    other.parent.mkdir()
    other.write_text("GCPFREE_PINS = [\n    {'gcpfree_eo.py': '" + "1" * 64 + "', 'modules': "
                     + repr(G.module_sha256()) + "},\n]\n")
    repo = d / "prereg_pin"
    repo.mkdir()
    git(repo, "init", "-q")
    pr = repo / "PREREGISTERED.txt"
    real, real_roots, real_prereg = G.PINS_FILE, list(G.ROOTS), G.PREREG
    G.ROOTS = roots_of(d / "roots")                     # the 'station's' folders of this test (never the station's)
    G.PREREG = pr                                       # this test's repository stands for the station's
    try:
        rc, out = main_q(["--preregister", "--outcome-if-one-camera", "per-camera", "--decided-by", "the owner"])
        G.PINS_FILE = other
        rc2, out2 = main_q(["--dry-run"])
        git(repo, "add", "-A")
        git(repo, "commit", "-q", "-m", "prereg")
        rc3, out3 = main_q(["--camera", "c2"])
        pr2 = d / "prereg_pin2" / "PREREGISTERED.txt"
        G.PREREG = pr2
        rc4, out4 = main_q(["--preregister", "--outcome-if-one-camera", "per-camera", "--decided-by", "the owner"])
        G.PINS_FILE, G.PREREG = d / "pins_none" / "survey_products.py", pr
        rc5, out5 = main_q(["--camera", "c2"])
    finally:
        G.PINS_FILE, G.ROOTS, G.PREREG = real, real_roots, real_prereg
    check(rc == 0 and pr.exists(), "the pinned version: --preregister writes the file")
    check(rc2 == 0 and "REFUSED: this gcpfree_eo.py" in out2 and "not the version survey_products.py pins" in out2,
          "--dry-run under another pin: says the procedure is not the pinned version")
    check(rc3 == 2 and "not the version survey_products.py pins" in out3 and not (repo / "runs").exists(),
          "a run under another pin: REFUSED (exit 2), no run folder")
    check(rc4 == 2 and not pr2.exists() and "nothing pre-registered" in out4,
          "--preregister under another pin: REFUSED, nothing written")
    check(rc5 == 2 and "cannot be read" in out5 and not (repo / "runs").exists(),
          "no survey_products.py to read the pins from: REFUSED")


def test_runs(d, canned):
    print("runs: the committed pre-registration (a trial too), its pinned photos, run folders never overwritten, one "
          "more run only with Aug photos")
    real_run_camera, real_roots = G.run_camera, list(G.ROOTS)
    scen, solved = {}, []

    def fake_run_camera(cam, roots, aug_dir, log):
        solved.append(cam)
        r = copy.deepcopy(scen[(cam, bool(aug_dir))])
        log(f"{cam}: replayed ({r['verdict']})")
        return r, (G.eo_rad(r["B"]) if r.get("B") else None)
    G.run_camera = fake_run_camera
    G.ROOTS = roots_of(d / "roots")                     # the 'station's' folders of this test
    WATCH["on"] = True
    repo = d / "prereg"
    repo.mkdir()
    git(repo, "init", "-q")
    pr = repo / "calibration" / "gcpfree" / "PREREGISTERED.txt"
    real_prereg, G.PREREG = G.PREREG, pr                # this test's repository stands for the station's
    base = []
    try:
        scen[("c1", False)] = nd_of(canned["c1_aug_only"])
        scen[("c2", False)] = canned["c2"]
        rc, out = main_q(base)
        check(rc == 2 and "does not exist" in out, "no pre-registration: every run REFUSED (exit 2)")
        # the review's peek: a full run before the rule exists, as a trial (one camera, or other folders): refused,
        # nothing solved, nothing written (the rule is fixed and committed before this code solves any photo)
        for argv, what in ((["--camera", "c2"], "one camera"), (["--roots", str(d / "roots" / "pics")], "other --roots")):
            rc, out = main_q(argv)
            check(rc == 2 and "does not exist" in out and "a trial neither" in out and not solved
                  and not (pr.parent / "runs").exists(),
                  f"no pre-registration: a TRIAL run ({what}) REFUSED too, nothing solved, no run folder")
        import io as io_
        import contextlib
        try:
            with contextlib.redirect_stderr(io_.StringIO()):
                G.main(["--preregistered", str(d / "shadow" / "PREREGISTERED.txt"), "--dry-run"])
            gone = False
        except SystemExit as exc:
            gone = exc.code == 2
        check(gone, "no --preregistered option any more: every run is made under this repository's file (PREREG)")
        rc, out = main_q(base + ["--preregister", "--decided-by", "the owner"])
        check(rc == 2 and "must be one of" in out and not pr.exists(), "--preregister without the owner's rule: refused")
        rc, out = main_q(base + ["--preregister", "--outcome-if-one-camera", "per-camera", "--decided-by", "the owner",
                                 "--roots", str(d / "elsewhere")])
        check(rc == 2 and "station's three photo folders" in out and not pr.exists(),
              "--preregister naming other photo folders: refused")
        G.ROOTS = roots_of(d / "roots") + [str(d / "not_mounted")]
        rc, out = main_q(base + ["--preregister", "--outcome-if-one-camera", "per-camera", "--decided-by", "the owner"])
        G.ROOTS = roots_of(d / "roots")
        check(rc == 2 and "photo folder(s) missing" in out and not pr.exists(),
              "--preregister with a station folder missing: refused (the photos it pins would not be a run's)")
        rc, out = main_q(base + ["--preregister", "--outcome-if-one-camera", "per-camera", "--decided-by", "the owner"])
        txt = pr.read_text()
        pinned = G.read_preregistration(pr)["photos"]
        n_c1 = sum(1 for r in pinned if r[0] == "c1")
        check(rc == 0 and "outcome_if_one_camera: per-camera" in txt and f"procedure_sha256: {G.own_sha256()}" in txt
              and "COMMIT AND PUSH THIS FILE" in txt and "&& git push" in txt and "commit AND push it" in out
              and f"-C {pr.parent} push" in out
              and pinned == G.photo_manifest(G.ROOTS) and n_c1 == 15 and len(pinned) == 39
              and f"photos_sha256: {G.photos_sha256(pinned)}" in txt
              and all(len(r[5]) == 64 and isinstance(r[4], int) for r in pinned),
              "--preregister: the owner's rule, the procedure's and its modules' sha256, the photo folders and the 39 "
              "photos a run selects (c1 15: Oct and Jan; c2 24: and Dec), each with its size and sha256; its header "
              "and the console say commit AND push")
        rc, out = main_q(base + ["--preregister", "--outcome-if-one-camera", "both-or-none", "--decided-by", "x"])
        check(rc == 2 and "exists" in out, "... written once (a second --preregister is refused)")
        rc, out = main_q(base)
        check(rc == 2 and "not committed" in out and not (pr.parent / "runs").exists(),
              "not committed yet: REFUSED, nothing written")
        rc, out = main_q(base + ["--camera", "c2"])
        check(rc == 2 and "not committed" in out and not solved, "... a trial too")
        git(repo, "add", "-A")
        git(repo, "commit", "-q", "-m", "prereg")
        pr.write_text(txt.replace("decided_by: the owner", "decided_by: the owner (13 Oct)"))
        rc, out = main_q(base)
        check(rc == 2 and "changed after its last commit" in out, "changed after its commit: REFUSED")
        git(repo, "commit", "-q", "-am", "prereg: who decided")
        head = git(repo, "rev-parse", "HEAD").stdout.strip()
        rc, out = main_q(base)
        rc2, out2 = main_q(base + ["--camera", "c2"])
        rc3, out3 = main_q(base + ["--dry-run"])
        check(rc == 2 and rc2 == 2 and f"its commit {head[:10]} is on no remote-tracking branch (not pushed)" in out
              and "push it (git push) BEFORE any run, a trial too" in out and "on no remote-tracking branch" in out2
              and not solved and not (pr.parent / "runs").exists()
              and "a run (a trial too) is REFUSED" in out3 and "git push" in out3,
              "committed but NOT pushed: every run REFUSED, a trial too, nothing solved (the time GitHub receives the "
              "push is the rule's only public timestamp); the dry run says so")
        push(repo)
        rc, out = main_q(base + ["--dry-run"])
        check(rc == 0 and "OK" in out and "pushed: on origin/" in out and "DELIVERY" in out
              and not (pr.parent / "runs").exists()
              and "a run would not be refused for its photos" in out
              and "the 39 photos selected now are the pre-registered ones" in out,
              "--dry-run: the pre-registration OK (2 commits, no run yet), a DELIVERY run, its photos there and the "
              "pinned ones, nothing written")
        G.ROOTS = roots_of(d / "roots") + [str(d / "not_mounted")]
        rc, out = main_q(base + ["--dry-run"])
        G.ROOTS = roots_of(d / "roots")
        check(rc == 0 and "a delivery run would be REFUSED (exit 2): photo folder(s) missing" in out
              and not (pr.parent / "runs").exists(),
              "--dry-run with a station folder missing: says a delivery run would be REFUSED (exit 2)")
        # the pinned photos: one removed, or one of other content, after the pre-registration: refused
        oct_one = Path(G.ROOTS[0]) / pinned[0][3]
        keep = oct_one.read_bytes()
        try:
            oct_one.rename(oct_one.with_name(oct_one.name + ".away"))
            rc, out = main_q(base)
            rc2, out2 = main_q(base + ["--dry-run"])
            oct_one.with_name(oct_one.name + ".away").rename(oct_one)
            oct_one.write_bytes(keep + b"\0")
            rc3, out3 = main_q(base)
        finally:
            oct_one.write_bytes(keep)
        check(rc == 2 and "not the pre-registered ones" in out and "no longer selected" in out and not solved
              and not (pr.parent / "runs").exists() and "a run would be REFUSED (exit 2): the photos" in out2,
              "a pinned photo removed after the pre-registration: a run REFUSED before anything is solved (the dry run "
              "says so)")
        check(rc3 == 2 and "of other content (size or sha256)" in out3 and not solved,
              "a pinned photo of other content (one byte added): REFUSED")
        # a TRIAL run (one camera), under this repository's pre-registration: reports only, recorded in its runs/
        rc, out = main_q(base + ["--camera", "c2"])
        runs = sorted((pr.parent / "runs").iterdir())
        tj = json.loads((runs[0] / "run.json").read_text())
        check(rc == 4 and len(runs) == 1 and runs[0].name.endswith("_trial") and not list(runs[0].glob("*_EO.yaml"))
              and tj["status"] == "finished" and tj["preregistration"]["path"] == str(pr.resolve())
              and tj["photos_used"] == [r for r in pinned if r[0] == "c2"]
              and "Commit and push it too, at once" in out and "git push" in out,
              "--camera c2: a TRIAL run (exit 4) under this repository's pre-registration, recorded in its runs/ "
              "(its report, its photos, no B file); commit it too")
        git(repo, "add", "-A")
        git(repo, "commit", "-q", "-m", "the trial")
        push(repo)
        # the rule cannot change once a run (a trial too) was made under it
        pr.write_text(pr.read_text().replace("outcome_if_one_camera: per-camera", "outcome_if_one_camera: both-or-none"))
        git(repo, "commit", "-q", "-am", "prereg: other rule")
        push(repo)
        rc, out = main_q(base)
        check(rc == 2 and "another version of the pre-registration" in out,
              "the rule changed after a trial was made: REFUSED (a rule is never chosen after seeing a result)")
        git(repo, "checkout", "-q", "HEAD~1", "--", str(pr.relative_to(repo)))
        git(repo, "commit", "-q", "-am", "prereg: as it was")
        push(repo)
        # delivery run 1: c1 NOT DELIVERED, c2 DELIVERED (per-camera): a c2-only B
        scen[("c1", False)] = nd_of(canned["c1_aug_only"])
        rc, out = main_q(base)
        r1 = sorted(p for p in (pr.parent / "runs").iterdir() if not p.name.endswith("_trial"))
        check(rc == 3 and len(r1) == 1 and "and push it (git add" in out and "&& git push)" in out,
              f"delivery run 1 (exit {rc}): c2 delivered, c1 not; it says commit AND push the run folder")
        r1 = r1[0]
        b2 = r1 / G.OUT_EO.format(cam="c2")
        check(b2.exists() and not (r1 / G.OUT_EO.format(cam="c1")).exists()
              and (r1 / (G.OUT_REPORT.format(cam="c1") + ".json")).exists(),
              "per-camera: c2's B file written, c1 none, c1's report kept beside it")
        check(not re.search(r"\d{8}", b2.name), f"B's name carries no solve date ({b2.name})")
        txt = b2.read_text()
        rep = json.loads((r1 / (G.OUT_REPORT.format(cam="c2") + ".json")).read_text())
        e = load_extrinsics(b2)
        ang = np.degrees(e[3:6])
        ang[0] %= 360.0
        check(np.allclose(e[:3], rep["B"][:3], atol=1e-4) and np.allclose(ang, rep["B"][3:6], atol=1e-4)
              and rep["delivery"] == "DELIVERED" and rep["b_file_sha256"] == G.sha256(b2)
              and rep["run"]["variant"] == "gcpfree" and rep["run"]["eligible"] is True,
              "its report beside it: the same B, DELIVERED, the B file's sha256, variant gcpfree, a delivery run")
        for words in (f"Made by gcpfree_eo.py (sha256 {G.own_sha256()})", "CACO03_c2_20250123_EO.yaml rotated by",
                      "Reference EO: CACO03_c2_20250123_EO.yaml", "No Oct 2024 GCP used", "acceptance: PASS",
                      "--calibration-variant gcpfree (the delivered B)", "not a validated bound on real photos",
                      "inherited camera position", "per-camera", "photo folders",
                      f"pre-registered constants sha256 {G.constants_sha256()}",
                      "EXPOSURE (stated by the people who made this procedure", "A minus the Jan 2025 solution",
                      "the horizon pass rule was changed from |d| <= tolerance",
                      "five position-inheritance cases that size the stated POSITION term were written",
                      "c1's stated DEM position term (|dz| median 0.011 m, p90 0.026 m) comes from the -0.55 m case"):
            check(words in txt, f"B file: '{words}'")
        rj = json.loads((r1 / "run.json").read_text())
        check(rj["argv"][1:] == base and rj["roots"] == G.ROOTS and len(rj["preregistration"]["git_history"]) == 4
              and rj["constants_sha256"] == G.constants_sha256() and rj["constants"]["ROOTS"] == G.ROOTS
              and "not the pinned ones" in rj["constants_problem"] and rj["exposure"] == G.EXPOSURE
              and rep["exposure"] == G.EXPOSURE and rj["photos_used"] == pinned
              and rj["photos_used_sha256"] == rj["preregistration"]["photos_sha256"] == G.photos_sha256(pinned)
              and all(m in rj["loaded_modules"] for m in G.MODULES + ("gcpfree_eo.py",)),
              "run.json: argv, the photo folders, the pre-registration and its git history (4 commits), the "
              "constants in effect (this test's photo folders: not the pinned ones, said), the exposure statement, "
              "the photos it used (the pinned ones), every repository module loaded with its sha256")
        git(repo, "add", "-A")
        git(repo, "commit", "-q", "-m", "delivery run 1")
        # no second run without new evidence; never over an earlier run
        rc, out = main_q(base)
        check(rc == 2 and "B is fixed by it" in out, "a second delivery run without Aug photos: REFUSED")
        class FixedNow(datetime):
            @classmethod
            def now(cls, tz=None):
                return datetime(2026, 10, 14, 1, 2, 3, tzinfo=timezone.utc)
        real_dt, G.datetime = G.datetime, FixedNow
        try:
            G.new_run_dir(d / "overwrite", True)
            G.new_run_dir(d / "overwrite", True)
            check(False, "a run folder is never written over")
        except G.Refused as exc:
            check("never writes over an earlier one" in str(exc), "a run folder is never written over (same id: refused)")
        finally:
            G.datetime = real_dt
        (d / "aug_empty").mkdir()
        rc, out = main_q(base + ["--aug-dir", str(d / "aug_empty")])
        check(rc == 2 and "no usable CACO03 c1 photo of 25-29 Aug 2024" in out and "<epoch>" in out,
              "--aug-dir with no usable photo: REFUSED (exit 2), the naming Chris's photos need")
        # delivery run 2 with Aug photos: c1 solved, c2 kept; a second product (gcpfree-aug)
        scen[("c1", True)] = canned["c1_aug_only"]
        time.sleep(1.1)                                 # run ids are UTC seconds
        rc, out = main_q(base + ["--aug-dir", str(d / "aug_c1")])
        r2 = sorted(p for p in (pr.parent / "runs").iterdir() if not p.name.endswith("_trial") and p != r1)
        check(rc == 0 and len(r2) == 1 and (r2[0] / G.OUT_EO.format(cam="c1")).exists()
              and (r2[0] / G.OUT_EO.format(cam="c2")).exists(), f"delivery run 2 (Aug photos, exit {rc}): c1 and c2")
        t2 = (r2[0] / G.OUT_EO.format(cam="c2")).read_text()
        t1 = (r2[0] / G.OUT_EO.format(cam="c1")).read_text()
        rj2 = json.loads((r2[0] / "run.json").read_text())
        check(f"solved in run {r1.name}" in t2 and "--calibration-variant gcpfree-aug (a SECOND product" in t1
              and "Reference EO: CACO03_c1_20240827_EO.yaml" in t1 and len(rj2["aug_photos_used"]) == 9
              and rj2["aug_photos_sha256"] == G.photos_sha256(rj2["aug_photos_used"]),
              "c2 kept from run 1 ('solved in run ...'); c1 from the Aug solve; a SECOND product: gcpfree-aug; the 9 "
              "Aug photos it used recorded (names, sizes, sha256)")
        rc, out = main_q(base + ["--aug-dir", str(d / "aug_c1")])
        check(rc == 2 and "no third delivery run" in out, "a third delivery run: REFUSED")
        # both-or-none: c2 passes, c1 does not -> nothing delivered, c2 WITHHELD in its report
        repo2 = d / "prereg2"
        repo2.mkdir()
        git(repo2, "init", "-q")
        pr2 = repo2 / "PREREGISTERED.txt"
        G.PREREG = pr2
        main_q(["--preregister", "--outcome-if-one-camera", "both-or-none", "--decided-by", "the owner"])
        git(repo2, "add", "-A")
        git(repo2, "commit", "-q", "-m", "prereg")
        push(repo2)
        rc, out = main_q([])
        rr = sorted((repo2 / "runs").iterdir())[0]
        rep2 = json.loads((rr / (G.OUT_REPORT.format(cam="c2") + ".json")).read_text())
        check(rc == 4 and not list(rr.glob("*_EO.yaml")) and rep2["delivery"].startswith("WITHHELD")
              and rep2["verdict"] == "DELIVERED",
              "both-or-none, only c2 passed: no B file at all, c2 WITHHELD (its result kept in its report)")
        # bad input is never a verdict
        repo3 = d / "prereg3"
        repo3.mkdir()
        git(repo3, "init", "-q")
        pr3 = repo3 / "PREREGISTERED.txt"
        G.PREREG = pr3
        main_q(["--preregister", "--outcome-if-one-camera", "per-camera", "--decided-by", "the owner"])
        git(repo3, "add", "-A")
        git(repo3, "commit", "-q", "-m", "prereg")
        push(repo3)
        owg = Path(G.ROOTS[1])
        owg.rename(owg.with_name("owg_unmounted"))
        try:
            n0 = len(solved)
            rc, out = main_q([])
        finally:
            owg.with_name("owg_unmounted").rename(owg)
        check(rc == 2 and "photo folder(s) missing" in out and not (repo3 / "runs").exists() and len(solved) == n0,
              "a delivery run with a photo folder missing (unmounted): REFUSED, no run")
        rc, out = main_q(["--camera", "c1", "--roots", str(d / "roots_other" / "pics")])
        check(rc == 2 and "no usable CACO03 c1 photo of 20-26 Oct 2024" in out,
              "no Oct photo for a camera: REFUSED (exit 2), never NOT DELIVERED")
        # a link to calibration A under the Jan EO's name: refused before anything is solved
        jan = G.CAL / "CACO03_c2_20250123_EO.yaml"
        keep = jan.read_bytes()
        jan.unlink()
        os.symlink(str(G.CAL / "CACO03_c2_20241023_EO.yaml"), str(jan))
        try:
            G.ROOTS, G.PREREG = roots_of(d / "roots"), pr2
            rc, out = main_q([])
            check(rc == 2 and "a link to CACO03_c2_20241023_EO.yaml" in out,
                  "the Jan EO replaced by a link to calibration A: REFUSED before anything is solved")
        finally:
            jan.unlink()
            jan.write_bytes(keep)
        test_foreign_history(d, scen)
    finally:
        WATCH["on"] = False
        G.run_camera, G.ROOTS, G.PREREG = real_run_camera, real_roots, real_prereg


def test_foreign_history(d, scen):
    print("records outside this branch's history (a branch made from before the pre-registration, a stash, the "
          "reflog): refused; a rebased copy of the same records is not")
    repo = d / "branch_away"
    repo.mkdir()
    git(repo, "init", "-q")
    (repo / "README").write_text("the station's repository (a test)\n")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "before the pre-registration")
    first = git(repo, "rev-parse", "HEAD").stdout.strip()
    pr = repo / "calibration" / "gcpfree" / "PREREGISTERED.txt"
    G.PREREG = pr
    rc, _ = main_q(["--preregister", "--outcome-if-one-camera", "both-or-none", "--decided-by", "the owner"])
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "gcpfree: pre-registration")
    push(repo)
    rc, out = main_q([])
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "gcpfree: delivery run 1")
    main_branch = git(repo, "rev-parse", "--abbrev-ref", "HEAD").stdout.strip()
    check(rc == 4 and G.foreign_records(pr.parent) == ([], ""),
          "the first delivery run, committed on its branch: nothing outside the branch's history")
    git(repo, "checkout", "-q", "-b", "alt", first)
    check(not pr.exists(), "a branch made from before the pre-registration: no pre-registration in its tree")
    rc, out = main_q(["--preregister", "--outcome-if-one-camera", "per-camera", "--decided-by", "the owner, again"])
    found, _ = G.foreign_records(pr.parent)
    check(rc == 2 and not pr.exists() and "outside this branch's history" in out and len(found) == 2,
          f"--preregister there: REFUSED (the other branch's pre-registration and run records: {len(found)} found)")
    git(repo, "branch", "-q", "-D", main_branch)
    found, _ = G.foreign_records(pr.parent)
    check(len(found) == 2, "... still seen with that branch deleted (the reflog keeps its commits)")
    # a rebase-like copy of the same records onto this branch: the same content, not foreign
    for c in git(repo, "log", "--reflog", "--topo-order", "--reverse", "--format=%H", "--", "calibration").stdout.split():
        git(repo, "cherry-pick", c)
    found, why = G.foreign_records(pr.parent)
    rc, out = main_q(["--dry-run"])
    check(not found and pr.exists() and "outside this branch's history" not in out,
          "the same records cherry-picked onto the branch (as a rebase copies them): nothing foreign")


def test_py38():
    import ast
    for f in ("gcpfree_eo.py", "test_gcpfree_eo.py"):
        ast.parse((HERE / f).read_text(), feature_version=(3, 8))
        check(True, f"{f} parses as Python 3.8")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="keep the test folder under this folder")
    args = ap.parse_args()
    if args.keep:
        Path(args.keep).mkdir(parents=True, exist_ok=True)
    d = Path(os.path.realpath(tempfile.mkdtemp(prefix="gcpfree_test_", dir=args.keep)))
    cal = d / "calibration"
    cal.mkdir()
    for f in G.ARCHIVED_SHA256:
        shutil.copy(str(HERE / "calibration" / f), str(cal / f))
    # decoys: calibration A, the history, a GCP file: an open of any fails the test (audit hook)
    for c in ("c1", "c2"):
        (cal / f"CACO03_{c}_20241023_EO.yaml").write_text("x: 0\ny: 0\nz: 0\nazimuth: 0\ntilt: 0\nroll: 0\n# DECOY\n")
    (cal / "calibration_history.csv").write_text("file,camera,x,y,z\nDECOY,c1,0,0,999\n")
    (cal / "2024-10-23_Marconi_Extrinsic_Targets_c1_xyz.csv").write_text("1,0,0,999\n")
    station_cal, G.CAL = G.CAL, cal
    unwatch = watch_cv2()
    start = os.getcwd()
    os.chdir(str(d))
    print(f"test folder: {d}")
    t0 = datetime.now()
    try:
        test_py38()
        test_watch_sees_reads(d)
        test_guard(d)
        test_scan_never_opens(d)
        test_archived(d)
        test_ransac()
        test_verdict_rules()
        test_uncertainty()
        test_outcome_rules()
        test_horizon_band(d)
        canned = {"c2": test_c2_direct_and_chain(d)}
        test_c1_with_aug(d)
        canned["c1_aug_only"] = test_c1_aug_only(d)
        test_not_delivered(d)
        test_pinned_version(d)
        test_runs(d, canned)
        if hasattr(sys, "addaudithook"):
            check(not WATCH["opened"], "nothing named like calibration A, a GCP file or the history was ever opened")
    finally:
        os.chdir(start)
        G.CAL = station_cal
        unwatch()
        if not args.keep:
            shutil.rmtree(str(d), ignore_errors=True)
    print(f"all gcpfree_eo tests passed ({N_OK[0]} checks, {(datetime.now() - t0).total_seconds():.0f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

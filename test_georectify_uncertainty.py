#!/usr/bin/env python3
"""
Tests For The Georectification Uncertainty
=============================================
Plain python3 (no pytest needed on the station):

    python3 test_georectify_uncertainty.py [--keep DIR] [--quick]

WHAT.
  * georectify_uncertainty.py on a SYNTHETIC camera, where the answers
    are known: the solve recovers the camera exactly from exact picks and
    leaves every target out at ~0 m; with pick noise of known size, the
    Jacobian covariance matches the scatter of the solves over many noise
    draws; the Monte Carlo spread matches the geometry by hand (an azimuth
    error moves a point range x angle across, a tilt error (h^2 + r^2)/h
    x angle along); the pixel footprint matches the pinhole by hand; the
    table writes, reads back and interpolates; a file named for 23 Oct
    2024 is refused and never opened.
  * dem_from_contours.py --sigma-obs on a synthetic planar beach:
    sqrt(SE^2 + sigma_wl^2 + sigma_setup^2 + g' S g) cell by cell, the
    beach's own slope recovered, each cell's geometry the mean over its
    FRAMES, a camera without a table, the frame and day filters; without
    the switch nothing is written or changed.
  * --evaluate on a synthetic DEM and survey: about 68 / 95 % within 1 / 2
    sigma when the residuals are drawn from sigma_obs.
  * The REAL calibration (smoke test): the in-sample numbers are
    georectify.py --self-test's, leave-one-out RMS c1 0.59 m (6
    parameters), c2 0.41 m (3 angles), the error grows with range, and
    the committed calibration/<cam>_geom_sigma.csv are what the script
    writes now (--quick skips that last, slow, check).

STATION FILES. calibration/ (the 2025-11-13 targets and picks, the IO/EO
YAMLs and the two tables) is read by the smoke test; nothing is written
outside a temporary folder.
"""

import argparse
import csv
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

import georectify_uncertainty as GU          # noqa: E402
from georectify import pixel_to_ground, load_intrinsics, load_extrinsics   # noqa: E402

FAILED = []
WATCH = {"on": False, "opened": []}


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILED.append(msg)


def run(args, d):
    r = subprocess.run([sys.executable] + [str(a) for a in args], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
    return r.returncode, r.stdout


def _audit(event, args):
    if WATCH["on"] and event == "open" and args and isinstance(args[0], (str, bytes, os.PathLike)):
        p = os.fsdecode(args[0])
        for n in (os.path.basename(p), os.path.basename(os.path.realpath(p))):
            if GU.HELD_NAME.search(n):
                WATCH["opened"].append(p)


# ---------------------------------------------------------------------
# A synthetic camera
# ---------------------------------------------------------------------

def synthetic_camera(distorted=True, roll_deg=1.0):
    """A camera like the station's: 2448 x 2048, f 2500 px, 20 m up, looking along the shore normal."""
    io = np.array([2448, 2048, 1224.0, 1024.0, 2500.0, 2500.0] +
                  ([-0.2, 0.2, -0.05, 0.001, 0.002] if distorted else [0.0] * 5))
    eo = np.array([1000.0, 2000.0, 20.0, np.radians(75.8), np.radians(72.0), np.radians(roll_deg)])
    return io, eo


def synthetic_targets(io, eo, n=15, seed=1):
    """Targets about 20-200 m out, across the lower two thirds of the frame, on planes 0-4 m up."""
    rng = np.random.RandomState(seed)
    u = rng.uniform(150, io[0] - 150, n)
    v = rng.uniform(560, io[1] - 80, n)
    z = rng.uniform(0.0, 4.0, n)
    E, N = pixel_to_ground(u, v, z, io, eo)
    xyz = np.c_[E, N, z]
    U, V = GU.project(xyz, io, eo)
    return xyz, np.c_[U, V]


def test_fit():
    print("fit_eo / leave_one_out on a synthetic camera")
    io, eo = synthetic_camera()
    xyz, uv = synthetic_targets(io, eo)
    rngs = np.hypot(xyz[:, 0] - eo[0], xyz[:, 1] - eo[1])
    check(np.isfinite(xyz).all() and rngs.min() > 20 and rngs.max() < 200,
          f"15 targets {rngs.min():.0f}-{rngs.max():.0f} m out, all in view")
    check(rngs.max() > 100, "some of them past 100 m")
    start = eo + np.r_[0.3, -0.3, 0.2, np.radians([0.3, -0.3, 0.3])]
    f = GU.fit_eo(io, start, xyz, uv, [0, 1, 2, 3, 4, 5])
    d = f["eo"] - eo
    check(f["converged"] and np.abs(d[:3]).max() < 1e-4 and np.degrees(np.abs(d[3:])).max() < 1e-5,
          f"exact picks, a start 0.3 m / 0.3 deg off: the solve returns the camera (|dx| {np.abs(d[:3]).max():.1e} m, "
          f"|dangle| {np.degrees(np.abs(d[3:])).max():.1e} deg)")
    fa = GU.fit_eo(io, eo + np.r_[0, 0, 0, np.radians([0.2, 0.2, -0.2])], xyz, uv, [3, 4, 5])
    check(np.degrees(np.abs(fa["eo"][3:] - eo[3:])).max() < 1e-5 and np.array_equal(fa["eo"][:3], eo[:3]),
          "angles only: the angles come back, the position is not touched")
    free, label, why, _ = GU.choose_free(io, eo, xyz, uv)
    check(free == [0, 1, 2, 3, 4, 5] and "reproduces" in why, f"auto: all six reproduce an exact camera ({label})")
    dE, dN, dpx, sols = GU.leave_one_out(io, eo, xyz, uv, [0, 1, 2, 3, 4, 5])
    check(np.hypot(dE, dN).max() < 1e-3 and dpx.max() < 1e-2 and sols.shape == (15, 6),
          f"exact picks: every target left out lands {np.hypot(dE, dN).max():.1e} m from where it is")
    # a calibration that held a position 0.5 m off and solved the angles: six parameters would move it
    # back, the angles' refit stays put, so auto solves the angles
    held = GU.fit_eo(io, eo + np.r_[0.5, 0, 0, 0, 0, 0], xyz, uv, [3, 4, 5])["eo"]
    free2, label2, why2, refits = GU.choose_free(io, held, xyz, uv)
    check(free2 == [3, 4, 5] and refits["all"][0] > 0.4 and refits["angles"][1] < 1e-4,
          f"a solve that held the position: auto picks {label2} (six would move {refits['all'][0]:.3f} m)")
    check(GU.choose_free(io, eo, xyz, uv, "angles")[0] == [3, 4, 5] and
          "--free angles" in GU.choose_free(io, eo, xyz, uv, "angles")[2], "--free angles is obeyed")


def test_covariance():
    print("covariances: the Jacobian against the scatter of solves over noise draws")
    io, eo = synthetic_camera()
    xyz, uv = synthetic_targets(io, eo)
    sigma_px = 2.0
    free = [0, 1, 2, 3, 4, 5]
    J = GU.fit_eo(io, eo, xyz, uv, free)["J"]
    c_true = sigma_px ** 2 * np.linalg.inv(J.T @ J)
    rs = np.random.RandomState(7)
    sols = []
    for _ in range(150):
        sols.append(GU.fit_eo(io, eo, xyz, uv + rs.normal(0, sigma_px, uv.shape), free)["eo"][free])
    emp = np.std(np.array(sols), axis=0, ddof=1)
    ratio = emp / np.sqrt(np.diag(c_true))
    check(np.all((ratio > 0.8) & (ratio < 1.25)),
          "150 draws of 2 px pick noise: the solves scatter as sigma^2 (J'J)^-1 says (ratio per parameter "
          + " ".join(f"{x:.2f}" for x in ratio) + ")")
    one = uv + np.random.RandomState(11).normal(0, sigma_px, uv.shape)
    _, _, _, loo = GU.leave_one_out(io, eo, xyz, one, free)
    covs = GU.covariances(io, eo, xyz, one, free, loo)
    rj = np.sqrt(np.diag(covs["jacobian"])) / np.sqrt(np.diag(c_true))
    check(set(covs) == {"jacobian", "jackknife"} and np.all((rj > 0.6) & (rj < 1.6))
          and np.all(np.linalg.eigvalsh(covs["jackknife"]) > 0),
          "on one draw: the Jacobian estimate (s^2 from its residuals) within 0.6-1.6 of the truth ("
          + " ".join(f"{x:.2f}" for x in rj) + "); the jackknife is a covariance (positive definite)")


def test_monte_carlo():
    print("monte_carlo / footprints against the geometry by hand (pinhole, no roll)")
    io, eo = synthetic_camera(distorted=False, roll_deg=0.0)
    h = eo[2]
    u, v, r, z = GU.table_nodes(io, eo, columns=[io[2]], ranges=(100, 200), zs=(0.0,))
    check(len(u) == 2 and np.allclose(r, [100, 200]), "the principal column reaches 100 and 200 m")
    E, N = pixel_to_ground(u, v, z, io, eo)
    check(np.allclose(np.hypot(E - eo[0], N - eo[1]), r, atol=0.02), "table_nodes finds the rows of those ranges")
    sa = np.radians(0.1)
    sxx, syy, sxy, share = GU.monte_carlo(io, eo, [3], np.array([[sa ** 2]]), u, v, z, samples=4000, seed=3)
    want = r * sa
    got = np.sqrt(sxx + syy)
    check(np.all(np.abs(got / want - 1) < 0.04) and np.all(share == 1),
          f"azimuth 0.1 deg: {got[0]:.3f} / {got[1]:.3f} m, by hand range x angle {want[0]:.3f} / {want[1]:.3f} m")
    st = np.radians(0.1)
    sxx, syy, sxy, _ = GU.monte_carlo(io, eo, [4], np.array([[st ** 2]]), u, v, z, samples=4000, seed=3)
    want = (h ** 2 + r ** 2) / h * st
    got = np.sqrt(sxx + syy)
    # along the view (azimuth 75.8 deg), so x and y are correlated +1
    corr = sxy / np.sqrt(sxx * syy)
    check(np.all(np.abs(got / want - 1) < 0.04) and np.all(corr > 0.99),
          f"tilt 0.1 deg: {got[0]:.3f} / {got[1]:.3f} m, by hand (h^2 + r^2)/h x angle {want[0]:.3f} / "
          f"{want[1]:.3f} m, along the view (corr {corr[0]:.3f})")
    normal, _ = GU.shore_normal()
    cross, frow, fcol = GU.footprints(io, eo, u, v, z, normal)
    theta = np.arctan(r / h)
    tilt = eo[4]
    want = h * np.cos(theta - tilt) ** 2 / (io[5] * np.cos(theta) ** 2)
    check(np.all(np.abs(frow / want - 1) < 0.005) and np.all(np.abs(cross / frow - 1) < 0.01),
          f"pixel footprint down the image {frow[0]:.3f} / {frow[1]:.3f} m, by hand h cos^2(theta - tilt) / "
          f"(f cos^2 theta) {want[0]:.3f} / {want[1]:.3f} m; looking along the shore normal it is the cross-shore one")
    want = np.hypot(r, h) * np.cos(theta - tilt) / io[4]
    check(np.all(np.abs(fcol / want - 1) < 0.005),
          f"across the image: slant range x cos(theta - tilt) / f ({fcol[0]:.3f} / {fcol[1]:.3f} m, by hand "
          f"{want[0]:.3f} / {want[1]:.3f} m)")


def test_table(d):
    print("analyse_camera / write_table / read_geom_sigma / geom_cov")
    d = Path(d) / "table"
    d.mkdir(parents=True, exist_ok=True)
    io, eo = synthetic_camera()
    xyz, uv = synthetic_targets(io, eo)
    uv = uv + np.random.RandomState(5).normal(0, 3.0, uv.shape)
    paths = write_camera_files(d, io, eo, xyz, uv)
    rows, head, res = GU.analyse_camera("c1", *paths, samples=400, ranges=(30, 50, 100, 150, 200, 300),
                                        zs=(-1.0, 1.0), log=lambda *a: None)
    out = d / "c1_geom_sigma.csv"
    GU.write_table(out, rows, head)
    t = GU.read_geom_sigma(out)
    keys = ("solve", "in_sample", "leave_one_out", "covariance", "eo_sd", "monte_carlo", "script_sha256", "io", "eo",
            "targets_xyz", "targets_uv", "grid_columns", "grid_ranges_m", "grid_z_navd88")
    check(all(k in t["provenance"] for k in keys) and abs(t["x"] - eo[0]) < 1e-3 and t["camera"] == "c1"
          and t["sxx"].shape == (2, 6, len(GU.grid_columns(io))),
          "the header carries the provenance and the grid; read back as z x range x column")
    with open(out, newline="") as f:
        data = list(csv.DictReader(line for line in f if not line.startswith("#")))
    check(list(data[0]) == GU.COLUMNS and len(data) == len(rows), f"{len(rows)} rows with the documented columns")
    r0 = next(r for r in data if r["sigma_x_m"] and float(r["range_m"]) == 100 and r["z_navd88"] == "1")
    u0, rg, z0 = float(r0["pixel_column"]), float(r0["range_m"]), float(r0["z_navd88"])
    sxx, syy, sxy = GU.geom_cov(t, u0, rg, z0)
    sx, sy, c = float(r0["sigma_x_m"]), float(r0["sigma_y_m"]), float(r0["corr_xy"])
    check(abs(sxx - sx * sx) < 1e-9 and abs(syy - sy * sy) < 1e-9 and abs(sxy - c * sx * sy) < 1e-9,
          f"at a node geom_cov returns the node ({u0:.0f} px, {rg:g} m, z {z0:g})")
    iz, ir, ic = (int(np.searchsorted(t[k], x)) for k, x in (("zs", z0), ("ranges", rg), ("columns", u0)))
    if ir + 1 < len(t["ranges"]) and np.isfinite(t["sxx"][iz, ir + 1, ic]):
        mid = GU.geom_cov(t, u0, (t["ranges"][ir] + t["ranges"][ir + 1]) / 2, z0)[0]
        check(abs(mid - (t["sxx"][iz, ir, ic] + t["sxx"][iz, ir + 1, ic]) / 2) < 1e-9,
              "half way between two ranges: the mean of the two")
    lo = GU.geom_cov(t, u0, rg, -5.0)[0]
    check(abs(lo - t["sxx"][0, ir, ic]) < 1e-9 and np.isnan(GU.geom_cov(t, np.nan, rg, z0)[0]),
          "below the lowest plane: that plane (clamped); no pixel column: NaN")
    sh = np.array([float(r["sigma_h_m"]) for r in data if r["z_navd88"] == "1" and r["sigma_h_m"]])
    rr = np.array([float(r["range_m"]) for r in data if r["z_navd88"] == "1" and r["sigma_h_m"]])
    check(np.median(sh[rr == 300]) > np.median(sh[rr == 100]) > np.median(sh[rr == 30]),
          f"the error grows with range (median {np.median(sh[rr == 30]):.2f} / {np.median(sh[rr == 100]):.2f} / "
          f"{np.median(sh[rr == 300]):.2f} m at 30 / 100 / 300 m)")
    # the CLI on the same files writes the same rows (same seed, same samples)
    (d / "cli").mkdir()
    rc, log = run([HERE / "georectify_uncertainty.py", "--camera", "c1", "--io", paths[0], "--eo", paths[1],
                   "--targets-xyz", paths[2], "--targets-uv", paths[3], "--samples", "300",
                   "--output-dir", d / "cli"], d)
    check(rc == 0 and (d / "cli" / "c1_geom_sigma.csv").exists() and "leave-one-out" in log and "wrote" in log,
          "the command line writes <cam>_geom_sigma.csv to --output-dir")
    rc2, log2 = run([HERE / "georectify_uncertainty.py", "--io", paths[0]], d)
    check(rc2 != 0 and "for one camera" in log2, "--io without --camera: refused (it is for one camera)")
    return paths


def write_camera_files(d, io, eo, xyz, uv, stem="synth"):
    iop, eop = Path(d) / f"{stem}_IO.yaml", Path(d) / f"{stem}_EO.yaml"
    iop.write_text("".join(f"{k}: {v!r}\n" for k, v in zip(["NU", "NV", "coU", "coV", "fx", "fy", "d1", "d2", "d3",
                                                             "t1", "t2"], io.tolist())))
    x, y, z, az, tilt, roll = eo[:3].tolist() + np.degrees(eo[3:]).tolist()
    eop.write_text(f"x: {x!r}\ny: {y!r}\nz: {z!r}\nazimuth: {az!r}\ntilt: {tilt!r}\nroll: {roll!r}\n")
    xp, up = Path(d) / f"{stem}_targets_xyz.csv", Path(d) / f"{stem}_picks_uv.csv"
    xp.write_text("".join(f"{100 + k},{x!r},{y!r},{z!r}\n" for k, (x, y, z) in enumerate(xyz.tolist())))
    up.write_text("num,U,V\n" + "".join(f"{100 + k},{a!r},{b!r}\n" for k, (a, b) in enumerate(uv.tolist())))
    return iop, eop, xp, up


def test_refused(d, paths):
    print("files named for 23 Oct 2024: refused, never opened")
    d = Path(d) / "held"
    d.mkdir(parents=True, exist_ok=True)
    held = d / "2024-10-23_Marconi_Extrinsic_Targets_c1_xyz.csv"
    shutil.copy(str(paths[2]), str(held))
    link = d / "innocent_targets.csv"
    os.symlink(str(held), str(link))
    survey = d / "2024-10-23_survey.asc"
    survey.write_text("ncols 1\n")
    WATCH["on"] = True
    try:
        outcomes = []
        for p in (held, link):
            try:
                GU.load_targets(p, paths[3])
                outcomes.append(False)
            except GU.Refused:
                outcomes.append(True)
        try:
            GU.evaluate(str(d / "nostem"), survey, log=lambda *a: None)
            outcomes.append(False)
        except GU.Refused:
            outcomes.append(True)
    finally:
        WATCH["on"] = False
    check(all(outcomes) and not WATCH["opened"],
          "targets by name, through a link with another name, a survey by name: Refused, and none was opened")
    rc, log = run([HERE / "georectify_uncertainty.py", "--camera", "c1", "--io", paths[0], "--eo", paths[1],
                   "--targets-xyz", held, "--targets-uv", paths[3], "--no-write"], d)
    check(rc == 2 and "REFUSED" in log, "the command line: exit 2, REFUSED")


# ---------------------------------------------------------------------
# dem_from_contours.py --sigma-obs
# ---------------------------------------------------------------------

HEADER = ("source_file,camera,capture_time_utc,capture_epoch,pixel_column,easting_utm19,northing_utm19,"
          "tide_elevation_navd88,offshore_hs_m,offshore_tp_s,setup_correction_m,beach_elevation_navd88,wave_source")
CAM_XY = (420000.0, 4638000.0)
GRAD = np.array([0.06, 0.03])                  # dz/dE, dz/dN of the synthetic beach: tan(beta) 0.067


def beach_rows(days=6, seed=2, bad_frame=None):
    """Both cameras every 30 min, 13-21 UTC: each frame's line is the contour of a planar beach
    (gradient GRAD, z 0 m at 420130 E, 4638130 N, 130-230 m from CAM_XY) at the frame's level, 60 m long;
    c1 the west part, c2 the east, overlapping. bad_frame: (day, k) of a frame given 0.5 m too much."""
    rng = np.random.RandomState(seed)
    w = 2 * np.pi / (12.42 * 3600)
    g2 = float(GRAD @ GRAD)
    t = np.array([-GRAD[1], GRAD[0]]) / np.sqrt(g2)
    c0 = np.array([420130.0, 4638130.0])
    rows = [HEADER]
    for day in range(days):
        for k in range(17):
            ep = int(datetime(2026, 10, 1 + day, 13, 0, tzinfo=timezone.utc).timestamp()) + 1800 * k
            level = 0.2 + 0.9 * np.sin(w * ep)
            when = datetime.fromtimestamp(ep, tz=timezone.utc)
            for cam, s_lo, s_hi in (("c1", -30, 6), ("c2", -6, 30)):
                dz = rng.normal(0, 0.02) + (0.5 if bad_frame == (day, k) and cam == "c1" else 0.0)
                base = c0 + level / g2 * GRAD
                for s in np.arange(s_lo, s_hi + 0.1, 0.5):
                    e, n = base + s * t
                    rows.append(f"{ep}.{cam}.timex,{cam},{when:%Y-%m-%dT%H:%M:%S}+00:00,{ep},{1000 + 20 * s:.1f},"
                                f"{e:.3f},{n:.3f},{level + dz - 0.05:.4f},1.0,8.0,0.05,{level + dz:.4f},"
                                f"marconi:hs_best/tp_s")
    return rows


def constant_table(path, cam, sx, sy, corr):
    """A <cam>_geom_sigma.csv with the same error at every node (the camera at CAM_XY)."""
    cols, rngs, zs = (0, 2447), (50, 100, 200, 300, 400), (-1, 2)
    head = [f"test table for {cam}", f"camera_x_utm19: {CAM_XY[0]}", f"camera_y_utm19: {CAM_XY[1]}",
            "covariance: constant -- test", "solve: none -- test",
            "grid_columns: " + " ".join(map(str, cols)), "grid_ranges_m: " + " ".join(map(str, rngs)),
            "grid_z_navd88: " + " ".join(map(str, zs))]
    rows = [[cam, u, r, z, "", "", "", sx, sy, corr, "", "", "", "", "", "1.000"]
            for z in zs for u in cols for r in rngs]
    GU.write_table(path, rows, head)


def read_asc(path):
    a = np.loadtxt(path, skiprows=6)
    return np.flipud(np.where(a == -9999, np.nan, a))            # row 0 = south, as build_grid


def test_local_slope():
    print("dem_from_contours.local_slope: the least-squares plane of 3 x 3 cells")
    import dem_from_contours as dfc
    yy, xx = np.mgrid[0:6, 0:7] * 2.0
    plane = 1.0 + 0.06 * xx + 0.03 * yy
    gE, gN = dfc.local_slope(plane, 2.0)
    check(np.allclose(gE, 0.06) and np.allclose(gN, 0.03), "a plane: its gradient everywhere, edges and corners too")
    holes = plane.copy()
    holes[2, :] = np.nan
    holes[:, 4] = np.nan
    lone = np.full((5, 5), np.nan)
    lone[2, 2] = 1.0
    line = np.full((5, 5), np.nan)
    line[2, :] = np.arange(5.0)
    gE2, _ = dfc.local_slope(holes, 2.0)
    check(np.allclose(gE2[np.isfinite(holes)], 0.06) and np.isnan(gE2[~np.isfinite(holes)]).all()
          and np.isnan(dfc.local_slope(lone, 2.0)[0]).all() and np.isnan(dfc.local_slope(line, 2.0)[1]).all(),
          "around gaps: from the cells there; a lone cell, or cells all in one line: NaN")


def test_frame_mean():
    print("dem_from_contours.observation_sigma_grid: a cell's geometry is the mean over its FRAMES")
    import dem_from_contours as dfc
    yy, xx = np.mgrid[0:3, 0:3] * 2.0
    dem = 0.1 * xx
    se = np.full((3, 3), 0.01)
    E = np.r_[np.full(10, 3.0), 3.2]
    N = np.full(11, 3.0)
    frames = np.array(["a"] * 10 + ["b"])
    geo = np.r_[np.tile([1.0, 0.0, 0.0], (10, 1)), [[3.0, 0.0, 0.0]]]
    out, rec, _ = dfc.observation_sigma_grid(dem, se, E, N, frames, geo, 0.0, 0.0, 2.0, 0.0, 0.0)
    want = np.sqrt(0.01 ** 2 + 0.1 ** 2 * 2.0)
    check(abs(out[1, 1] - want) < 1e-12 and np.isnan(out[0, 0]),
          f"frame a 10 points of sxx 1, frame b one of 3: sxx 2 (not 1.18 by points); sigma_obs {out[1, 1]:.4f} m")
    out2, _, _ = dfc.observation_sigma_grid(dem, se, E, N, frames, None, 0.0, 0.0, 2.0, 0.12, 0.21)
    check(np.isnan(out2).all(), "no georectification term at all: no value")


def build(d, name, rows, extra=()):
    g = Path(d) / f"{name}.csv"
    g.write_text("\n".join(rows) + "\n")
    rc, log = run([HERE / "dem_from_contours.py", g, Path(d) / name, "--no-plot"] + [str(x) for x in extra], d)
    (Path(d) / f"{name}.log").write_text(log)
    info = json.loads((Path(d) / f"{name}_info.json").read_text()) if (Path(d) / f"{name}_info.json").exists() else {}
    return rc, log, info


def test_dem_sigma_obs(d):
    print("dem_from_contours.py --sigma-obs on a synthetic planar beach")
    d = Path(d) / "dem"
    tables = d / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    sx, sy, corr = 2.0, 1.0, 0.5
    for cam in ("c1", "c2"):
        constant_table(tables / f"{cam}_geom_sigma.csv", cam, sx, sy, corr)
    rows = beach_rows()
    rc0, log0, info0 = build(d, "stemplain", rows)
    rc, log, info = build(d, "stemobs", rows, ["--sigma-obs", "--geom-sigma-dir", tables])
    check(rc0 == 0 and rc == 0 and (d / "stemobs_sigma_obs.asc").exists()
          and not (d / "stemplain_sigma_obs.asc").exists(),
          "with --sigma-obs: <stem>_sigma_obs.asc; without: none")
    same = all((d / f"stemplain_{k}.asc").read_bytes() == (d / f"stemobs_{k}.asc").read_bytes()
               for k in ("dem", "spread", "count", "sigma"))
    i0 = dict(info)
    rec = i0.pop("sigma_obs", None)
    check(same and "sigma_obs" not in info0 and i0 == info0,
          "the dem, spread, count and sigma grids are byte for byte the same, and _info.json but for sigma_obs")
    plain = log0.replace("stemplain", "STEM").splitlines()
    extra = [ln for ln in log.replace("stemobs", "STEM").splitlines() if ln not in plain]
    ours = ("Georectification  : ", "observation sigma : ", "                    georectification term",
            "                    no value on")
    check(len(extra) >= 4 and all(ln.startswith(ours) or not ln.strip() or (ln.startswith("wrote ") and
                                                                           "STEM_sigma_obs.asc (" in ln)
                                  for ln in extra) and "WARNING" not in log,
          f"the log differs only by the --sigma-obs lines ({len(extra)}), no WARNING")
    dem, se, so = (read_asc(d / f"stemobs_{k}.asc") for k in ("dem", "sigma", "sigma_obs"))
    import dem_from_contours as dfc
    gE, gN = dfc.local_slope(dem, 2.0)
    S = np.array([[sx * sx, corr * sx * sy], [corr * sx * sy, sy * sy]])
    var = gE * gE * S[0, 0] + gN * gN * S[1, 1] + 2 * gE * gN * S[0, 1]
    want = np.sqrt(se ** 2 + dfc.SIGMA_WL_M ** 2 + dfc.SIGMA_SETUP_M ** 2 + var)
    has = np.isfinite(so)
    check(has.sum() > 100 and np.array_equal(has, np.isfinite(se) & np.isfinite(gE))
          and np.allclose(so[has], want[has], atol=3e-4),
          f"sqrt(SE^2 + 0.12^2 + 0.21^2 + g' S g) on every cell with an SE and a local slope ({int(has.sum())} "
          f"cells), NaN elsewhere")
    tb = rec["tan_beta_median"] if rec else np.nan
    gt = np.sqrt(GRAD @ S @ GRAD)
    check(rec and abs(tb / np.hypot(*GRAD) - 1) < 0.15 and abs(rec["geometry_median_m"] / gt - 1) < 0.2,
          f"the beach's own slope: tan(beta) median {tb} (planted {np.hypot(*GRAD):.4f}), geometry term "
          f"{rec and rec['geometry_median_m']} m (planted sqrt(g' S g) {gt:.4f})")
    check(rec and rec["file"] == "_sigma_obs.asc" and rec["sigma_wl_m"] == dfc.SIGMA_WL_M
          and set(rec["tables"]) == {"c1", "c2"} and len(rec["tables"]["c1"]["sha256"]) == 64,
          "_info.json sigma_obs: the file, the defaults, the tables by sha256")
    # zero observation terms and a zero table: the SE itself
    for cam in ("c1", "c2"):
        constant_table(tables / f"{cam}_geom_sigma.csv", cam, 0.0, 0.0, 0.0)
    rcz, _, _ = build(d, "zero", rows, ["--sigma-obs", "--geom-sigma-dir", tables, "--sigma-wl", "0",
                                        "--sigma-setup", "0"])
    soz, sez = read_asc(d / "zero_sigma_obs.asc"), read_asc(d / "zero_sigma.asc")
    hz = np.isfinite(soz)
    check(rcz == 0 and np.allclose(soz[hz], sez[hz], atol=1e-4), "no observation error at all: sigma_obs = SE")
    # a camera without a table
    (tables / "c2_geom_sigma.csv").unlink()
    rc2, log2, info2 = build(d, "noc2", rows, ["--sigma-obs", "--geom-sigma-dir", tables])
    so2 = read_asc(d / "noc2_sigma_obs.asc")
    check(rc2 == 0 and "no " in log2 and "c2_geom_sigma.csv" in log2 and "NOTE: --sigma-obs" in log2
          and np.isfinite(so2).sum() < has.sum() and info2["sigma_obs"]["tables"].keys() == {"c1"},
          f"no table for c2: a NOTE, its cells without c1 points have no value ({int(np.isfinite(so2).sum())} of "
          f"{int(has.sum())} cells left)")
    # the frame and day filters keep each point's term beside it
    constant_table(tables / "c2_geom_sigma.csv", "c2", sx, sy, corr)
    rcf, logf, infof = build(d, "filters", beach_rows(bad_frame=(2, 8)),
                             ["--sigma-obs", "--geom-sigma-dir", tables, "--max-frame-offset", "0.2",
                              "--max-day-offset", "0.15"])
    check(rcf == 0 and "REJECTED" in logf and (d / "filters_sigma_obs.asc").exists()
          and infof["filters"]["frames_rejected"] >= 1,
          "--max-frame-offset rejects the planted frame, --max-day-offset runs: sigma_obs still built")
    # the series copy
    rcs, logs, _ = build(d, "ser", rows, ["--sigma-obs", "--geom-sigma-dir", tables, "--last-days", "6",
                                          "--series-dir", d / "series"])
    check(rcs == 0 and len(list((d / "series").glob("*_sigma_obs.asc"))) == 1, "the dated series copy too")


def test_evaluate(d):
    print("georectify_uncertainty.py --evaluate on a synthetic DEM and survey")
    import dem_from_contours as dfc
    d = Path(d) / "evaluate"
    d.mkdir(parents=True, exist_ok=True)
    rs = np.random.RandomState(9)
    yy, xx = np.mgrid[0:60, 0:60] * 1.0
    dem = 0.5 + 0.05 * xx
    so = 0.1 + 0.3 * xx / 59.0
    r = rs.normal(0, 1, dem.shape) * so
    for name, grid in (("t_dem", dem), ("t_sigma", np.full(dem.shape, 0.02)), ("t_sigma_obs", so),
                       ("survey", dem - r)):
        dfc.write_ascii_grid(str(d / f"{name}.asc"), grid, 420000.0, 4638000.0, 1.0)
    rc, log = run([HERE / "georectify_uncertainty.py", "--evaluate", d / "t", "--survey", d / "survey.asc"], d)
    rec = GU.evaluate(str(d / "t"), d / "survey.asc", log=lambda *a: None)
    w1, w2 = rec["raw|sigma_obs"]["within_1"], rec["raw|sigma_obs"]["within_2"]
    thirds = [t["rms_m"] for t in rec["thirds"]]
    check(rc == 0 and abs(w1 - 0.683) < 0.04 and abs(w2 - 0.954) < 0.02 and thirds == sorted(thirds)
          and rec["raw|sigma (SE)"]["within_1"] < 0.2 and "within 1 sigma" in log,
          f"residuals drawn from sigma_obs: {100 * w1:.1f} % within 1, {100 * w2:.1f} % within 2; RMS rises over the "
          f"thirds ({', '.join(f'{x:.3f}' for x in thirds)}); the SE alone covers {100 * rec['raw|sigma (SE)']['within_1']:.0f} %")
    # a point survey (an RTK file): each shot against the cell it falls in
    pr = np.random.RandomState(10)
    ce, cn = pr.randint(0, 60, 1500), pr.randint(0, 60, 1500)
    E, N = 420000.0 + ce + pr.uniform(0.05, 0.95, 1500), 4638000.0 + cn + pr.uniform(0.05, 0.95, 1500)
    zs = dem[cn, ce] - pr.normal(0, 1, 1500) * so[cn, ce]
    (d / "shots.csv").write_text("Name,Easting,Northing,Elevation\n" + "".join(
        f"p{k},{e:.3f},{n:.3f},{z:.4f}\n" for k, (e, n, z) in enumerate(zip(E, N, zs))) + "far,1,1,9\n")
    rp = GU.evaluate(str(d / "t"), d / "shots.csv", log=lambda *a: None)
    p1 = rp["raw|sigma_obs"]["within_1"]
    check(rp["unit"] == "shot" and rp["cells"] == 1500 and abs(p1 - 0.683) < 0.05,
          f"an RTK-style CSV: {rp['cells']} shots in DEM cells (one outside left out), {100 * p1:.1f} % within 1 sigma")
    rc2, log2 = run([HERE / "georectify_uncertainty.py", "--evaluate", d / "nothere", "--survey", d / "survey.asc"], d)
    check(rc2 != 0 and "not found" in log2, "a stem without _sigma_obs.asc: an error naming the file")


# ---------------------------------------------------------------------
# The real calibration
# ---------------------------------------------------------------------

def test_real(quick):
    print("the real 2025-11-13 calibration (smoke test)")
    cal = HERE / "calibration"
    want = {"c1": (0.328, 0.416, 1.181, 0.589, [0, 1, 2, 3, 4, 5]), "c2": (0.330, 0.416, 0.994, 0.410, [3, 4, 5])}
    for cam in ("c1", "c2"):
        paths = (cal / f"CACO05_{cam}_20240801_IO.yaml", cal / f"CACO05_{cam}_20251113_EO-CV.yaml",
                 cal / GU.TARGETS.format(cam=cam), cal / GU.PICKS.format(cam=cam))
        rows, head, res = GU.analyse_camera(cam, *paths, samples=300, ranges=(50, 100, 200, 300, 400), zs=(0.0,),
                                            log=lambda *a: None)
        mean, rms, mx, loo, free = want[cam]
        si = res["in_sample"]
        check(abs(si["mean"] - mean) < 0.002 and abs(si["rms"] - rms) < 0.002 and abs(si["max"] - mx) < 0.002,
              f"{cam}: in-sample mean {si['mean']:.3f} / RMS {si['rms']:.3f} / max {si['max']:.3f} m, as "
              f"georectify.py --self-test")
        check(res["free"] == free and abs(res["loo"]["rms"] - loo) < 0.005,
              f"{cam}: solve {res['solve'].split(' (')[0]}, leave-one-out RMS {res['loo']['rms']:.3f} m "
              f"(known_verdicts C-B6: {loo})")
        sh = {}
        for r in rows:
            if r[10]:
                sh.setdefault(int(r[1]), []).append((float(r[2]), float(r[10]), float(r[12])))
        grows = all(np.all(np.diff([x[1] for x in v]) > 0) and np.all(np.diff([x[2] for x in v]) > 0)
                    for v in sh.values() if len(v) > 1)
        check(grows and len(sh) >= 8, f"{cam}: in every column the error and the footprint grow with range")
        t = GU.read_geom_sigma(cal / f"{cam}_geom_sigma.csv")
        check(t["camera"] == cam and np.isfinite(t["sxx"]).sum() > 1000
              and t["provenance"]["eo"].startswith(f"CACO05_{cam}_20251113_EO-CV.yaml sha256 ")
              and t["provenance"]["eo"].endswith(GU.sha256(paths[1])),
              f"{cam}: the committed table reads, and names this EO by its sha256")
        if quick:
            continue
        rows2, head2, _ = GU.analyse_camera(cam, *paths, log=lambda *a: None)
        with open(cal / f"{cam}_geom_sigma.csv", newline="") as f:
            lines = f.read().splitlines()
        data = list(csv.reader(line for line in lines if not line.startswith("#")))[1:]
        same = len(data) == len(rows2) and all(
            a[:4] == b[:4] and all((x == "" and y == "") or (x != "" and y != "" and abs(float(x) - float(y)) <= 2e-3)
                                   for x, y in zip(a[4:], b[4:])) for a, b in zip(data, rows2))
        skip = ("written_utc", "git_head", "script_sha256")
        hc = [ln[2:] for ln in lines if ln.startswith("#") and not ln[2:].startswith(skip)]
        hn = [ln for ln in head2 if not ln.startswith(skip)]
        check(same and hc == hn, f"{cam}: the committed table is what the script writes now ({len(rows2)} nodes)")


def test_py38():
    print("Python 3.8 syntax")
    import ast
    for f in ("georectify_uncertainty.py", "test_georectify_uncertainty.py", "dem_from_contours.py"):
        src = (HERE / f).read_text(encoding="utf-8")
        try:
            ast.parse(src, feature_version=(3, 8))
        except TypeError:                        # Python 3.7 has no feature_version
            ast.parse(src)
        check(True, f"{f} parses as Python 3.8")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="keep the outputs, in a new folder inside this one")
    ap.add_argument("--quick", action="store_true", help="skip regenerating the committed tables (~20 s)")
    args = ap.parse_args()
    if args.keep:
        Path(args.keep).mkdir(parents=True, exist_ok=True)
    sys.addaudithook(_audit)
    d = os.path.realpath(tempfile.mkdtemp(prefix="geom_uncertainty_", dir=args.keep))
    print(f"test folder: {d}")
    try:
        test_fit()
        test_covariance()
        test_monte_carlo()
        paths = test_table(d)
        test_refused(d, paths)
        test_local_slope()
        test_frame_mean()
        test_dem_sigma_obs(d)
        test_evaluate(d)
        test_real(args.quick)
        test_py38()
    finally:
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED")
        return 1
    print("all georectification uncertainty tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

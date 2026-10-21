#!/usr/bin/env python3
"""
Tests For The Georectification Geometry
==========================================
Plain python3 (no pytest needed on the station):

    python3 test_georectify.py [--keep DIR]

WHAT.
  * "In front of the camera" (georectify.camera_depth, used by
    pixel_to_ground and view_reproject.ground_to_pixel). A synthetic
    camera at the station's UTM position (420,000 E, 4,638,000 N, 20 m
    up, 70 deg tilt, and again at a steep 25 deg, where the vertical part
    of the depth matters) is turned through every azimuth, 0-359 deg. At every
    azimuth each pixel must land where a ray cast independently lands,
    always in front of the camera; the principal point must land on the
    azimuth, at h tan(tilt) from the camera; ground_to_pixel must take
    those points back to their pixels as valid, and refuse a point behind
    the camera. The same camera moved to a local frame (camera at
    E = N = 0) must give the same ground points, shifted. Until Oct 2026
    the test was the sign of the projective scale. That sign depends on
    the camera's absolute position. With the code of 42e3e86, this camera
    lost every in-view pixel at azimuths of 96-275 deg, and
    ground_to_pixel refused every point there. It took sky pixels as
    ground points behind the camera, and accepted points behind the
    camera as in view. In the local frame it failed at every azimuth
    tried (0-300 deg). The checks here FAIL on that code.
  * The lens inversion's guards (georectify.undistort_uv,
    lens_fold_r2): on every IO file in calibration/ nothing is removed and
    the result is bit-identical to the inversion without guards; on a
    synthetic lens that folds inside the frame, a pixel beyond the fold or
    one the passes do not converge on comes back NaN, never a finite wrong
    point, while the inversion without guards returns finite points that
    miss their pixel. On a lens whose r_d rises again past its fold, the
    passes converge onto that outer branch; only the fold guard makes those
    pixels NaN. None of this, a NaN pixel included, raises a numpy warning.
  * georectify.py's command line on a camera aimed at 150 deg: every
    waterline point georectified (it used to write none and blame the
    horizon).

STATION FILES. Only the committed lens files calibration/*_IO.yaml are
read, for the lens checks. The georectification checks use a synthetic
camera, written to a temporary folder for the command line.
"""

import argparse
import ast
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import georectify as G          # noqa: E402
import view_reproject as VR     # noqa: E402

FAILED = []

# The synthetic camera: the CACO05 frame and focal length, no lens
# distortion (so a ray can be cast here without the lens inversion), at the
# station's position in UTM 19N.
NU, NV, F = 2448, 2048, 2000.0
PINHOLE = np.array([NU, NV, NU / 2, NV / 2, F, F, 0.0, 0.0, 0.0, 0.0, 0.0])
CAM_E, CAM_N, CAM_Z = 420000.0, 4638000.0, 20.0
TILT, ROLL = 70.0, -8.0
STEEP_TILT = 25.0           # near nadir: the depth's vertical part outweighs its horizontal part
Z_PLANE = 0.5


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILED.append(msg)


def eo_at(az_deg, e=CAM_E, n=CAM_N, tilt=TILT):
    return np.array([e, n, CAM_Z, np.radians(az_deg), np.radians(tilt), np.radians(ROLL)])


def grid(step=32):
    """Pixel centres every `step` px over the whole frame, off its very edge (a
    round trip may come back 1e-9 px outside the frame there)."""
    U, V = np.meshgrid(np.arange(step / 2.0, NU, step), np.arange(step / 2.0, NV, step))
    return U.ravel(), V.ravel()


def cast_rays(U, V, z, eo, max_range_m=1000.0):
    """The ground point of each (undistorted) pixel by casting its ray: the
    viewing direction R^T K^-1 [u v 1] from the camera, to the plane z, kept
    only when the plane is met IN FRONT (t > 0) and within max_range_m."""
    K = np.array([[-F, 0.0, NU / 2], [0.0, -F, NV / 2], [0.0, 0.0, 1.0]])
    R = G.cirn_angles_to_R(*eo[3:6])
    w = R.T @ np.linalg.solve(K, np.vstack([U, V, np.ones_like(U)]))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (z - eo[2]) / w[2]
    X, Y = eo[0] + t * w[0], eo[1] + t * w[1]
    ok = np.isfinite(t) & (t > 0) & (np.hypot(X - eo[0], Y - eo[1]) <= max_range_m)
    return np.where(ok, X, np.nan), np.where(ok, Y, np.nan)


def forward(eo):
    """The camera's viewing direction from its azimuth and tilt (tilt 0 =
    straight down, 90 = level): the axis a point must lie ahead on."""
    az, ti = eo[3], eo[4]
    return np.array([np.sin(ti) * np.sin(az), np.sin(ti) * np.cos(az), -np.cos(ti)])


# ---------------------------------------------------------------------
# In front of the camera, at every azimuth
# ---------------------------------------------------------------------

def test_every_azimuth(tilt=TILT):
    print(f"pixel_to_ground and ground_to_pixel at every azimuth, camera at the station's UTM position, "
          f"tilt {tilt:g} deg")
    U, V = grid()
    n_ref = None
    worst_xy = worst_px = worst_pp = 0.0
    lost, ghosts, wrong_mask, refused, behind_ok = [], [], [], [], []
    for az in range(0, 360):
        eo = eo_at(az, tilt=tilt)
        X, Y = G.pixel_to_ground(U, V, Z_PLANE, PINHOLE, eo)
        Xr, Yr = cast_rays(U, V, Z_PLANE, eo)
        fin, fin_r = np.isfinite(X), np.isfinite(Xr)
        if n_ref is None:
            n_ref = int(fin_r.sum())
        if fin_r.any() and not fin[fin_r].any():
            lost.append(az)
        if not np.array_equal(fin, fin_r):
            wrong_mask.append(az)
        both = fin & fin_r
        if both.any():
            worst_xy = max(worst_xy, float(np.max(np.hypot(X[both] - Xr[both], Y[both] - Yr[both]))))
        depth = (np.vstack([X[fin] - eo[0], Y[fin] - eo[1], np.full(int(fin.sum()), Z_PLANE - eo[2])]).T
                 @ forward(eo))
        if np.any(depth <= 0):
            ghosts.append(az)
        # the principal point: on the azimuth, at h tan(tilt)
        pe, pn = G.pixel_to_ground(NU / 2, NV / 2, Z_PLANE, PINHOLE, eo)
        if np.isfinite(pe[0]):
            d_bear = (np.degrees(np.arctan2(pe[0] - eo[0], pn[0] - eo[1])) - az + 180.0) % 360.0 - 180.0
            h_tan = (CAM_Z - Z_PLANE) * np.tan(np.radians(tilt))
            d_dist = np.hypot(pe[0] - eo[0], pn[0] - eo[1]) - h_tan
            worst_pp = max(worst_pp, abs(np.radians(d_bear)) * h_tan, abs(d_dist))   # across and along, in m
        else:
            worst_pp = np.inf
        # the ray cast's ground points, back to their pixels
        u2, v2, ok = VR.ground_to_pixel(Xr[fin_r], Yr[fin_r], Z_PLANE, PINHOLE, eo)
        if not ok.all():
            refused.append(az)
        elif ok.any():
            worst_px = max(worst_px, float(np.max(np.hypot(u2 - U[fin_r], v2 - V[fin_r]))))
        # a point 100 m BEHIND the camera on its optical axis: it projects onto the
        # principal point (the mirror image of a point ahead), but is never in view
        bx, by, bz = eo[0:3] - 100.0 * forward(eo)
        ub, vb, okb = VR.ground_to_pixel(bx, by, bz, PINHOLE, eo)
        if bool(okb) or abs(float(ub) - NU / 2) > 1e-6 or abs(float(vb) - NV / 2) > 1e-6:
            behind_ok.append(az)
    span = lambda a: f"{len(a)} azimuths" + (f" ({a[0]}-{a[-1]} deg)" if a else "")   # noqa: E731
    check(n_ref and n_ref > 1000, f"the reference ray cast meets the plane at {n_ref} of {U.size} pixels")
    check(not lost, f"no azimuth loses every in-view pixel (lost: {span(lost)})")
    check(not wrong_mask, f"at every azimuth exactly the pixels whose ray meets the plane in front, within 1000 m, "
                          f"are kept (not: {span(wrong_mask)})")
    check(worst_xy < 1e-6, f"... and land where the ray cast lands (worst {worst_xy:.1e} m)")
    check(not ghosts, f"no point returned behind the camera (at: {span(ghosts)})")
    check(worst_pp < 1e-5, f"the principal point lands on the azimuth, at h tan(tilt) from the camera "
                           f"(worst {worst_pp:.1e} m across or along)")
    check(not refused, f"ground_to_pixel takes every ground point the ray cast found back as valid "
                       f"(refused at: {span(refused)})")
    check(worst_px < 1e-6, f"... to its own pixel (worst {worst_px:.1e} px)")
    check(not behind_ok, f"a point behind the camera on its axis: projected onto the principal point, never valid "
                         f"(valid at: {span(behind_ok)})")


def test_local_frame():
    print("a camera at the origin of a local frame: the same ground points, shifted")
    U, V = grid()
    bad = []
    for az in (0, 54, 96, 150, 180, 270, 300):
        Xu, Yu = G.pixel_to_ground(U, V, Z_PLANE, PINHOLE, eo_at(az))
        Xl, Yl = G.pixel_to_ground(U, V, Z_PLANE, PINHOLE, eo_at(az, 0.0, 0.0))
        same = np.array_equal(np.isfinite(Xu), np.isfinite(Xl)) and np.isfinite(Xl).sum() > 1000
        if not same or np.nanmax(np.hypot(Xl + CAM_E - Xu, Yl + CAM_N - Yu)) > 1e-6:
            bad.append(az)
        _, _, ok = VR.ground_to_pixel(Xl[np.isfinite(Xl)], Yl[np.isfinite(Xl)], Z_PLANE, PINHOLE, eo_at(az, 0.0, 0.0))
        if not ok.all():
            bad.append(az)
    check(not bad, f"camera at E = N = 0, azimuths 0-300 deg: the UTM points minus the camera's UTM position, "
                   f"each valid in ground_to_pixel (not at: {bad})")


def test_station_lens_every_azimuth():
    print("the live c1 lens turned through every azimuth: the same pixels kept, the ground turned with it")
    io = G.load_intrinsics(HERE / "calibration" / "CACO05_c1_20240801_IO.yaml")
    U, V = grid(48)
    X0, Y0 = G.pixel_to_ground(U, V, Z_PLANE, io, eo_at(0))
    keep0 = np.isfinite(X0)
    # pixels whose ray meets the plane close to the 1000 m range cut could go either way on rounding
    with np.errstate(invalid="ignore"):
        near_cut = np.abs(np.hypot(X0 - CAM_E, Y0 - CAM_N) - 1000.0) < 1e-3
    bad, worst, back = [], 0.0, []
    for az in range(0, 360, 3):
        X, Y = G.pixel_to_ground(U, V, Z_PLANE, io, eo_at(az))
        a = np.radians(az)
        # turning the camera by az about the vertical turns the ground points by az about the camera
        Xe = CAM_E + (X0 - CAM_E) * np.cos(a) + (Y0 - CAM_N) * np.sin(a)
        Ye = CAM_N - (X0 - CAM_E) * np.sin(a) + (Y0 - CAM_N) * np.cos(a)
        if not np.array_equal(np.isfinite(X)[~near_cut], keep0[~near_cut]):
            bad.append(az)
        # the turned az-0 ground points, back to their pixels
        u2, v2, ok = VR.ground_to_pixel(Xe[keep0], Ye[keep0], Z_PLANE, io, eo_at(az))
        if not ok.all() or np.max(np.hypot(u2 - U[keep0], v2 - V[keep0])) > 1e-4:
            back.append(az)
        both = np.isfinite(X) & keep0
        worst = (max(worst, float(np.max(np.hypot(X[both] - Xe[both], Y[both] - Ye[both])))) if both.any()
                 else np.inf)
    check(keep0.sum() > 1000 and not bad, f"{int(keep0.sum())} of {U.size} pixels kept at every azimuth "
                                          f"(other pixels at: {bad})")
    check(worst < 1e-6, f"the ground points turn with the camera (worst {worst:.1e} m)")
    check(not back, f"ground_to_pixel takes them back to their pixels (not at: {back})")


# ---------------------------------------------------------------------
# The lens inversion's guards
# ---------------------------------------------------------------------

def unguarded(Ud, Vd, io, iterations=20):
    """undistort_uv as it was before its guards (to Oct 2026), for comparison."""
    c0U, c0V, fx, fy = io[2], io[3], io[4], io[5]
    d1, d2, d3, t1, t2 = io[6], io[7], io[8], io[9], io[10]
    xd = (np.asarray(Ud, dtype=float) - c0U) / fx
    yd = (np.asarray(Vd, dtype=float) - c0V) / fy
    x, y = xd.copy(), yd.copy()
    for _ in range(iterations):
        r2 = x * x + y * y
        radial = 1.0 + d1 * r2 + d2 * r2 * r2 + d3 * r2 * r2 * r2
        dx = 2.0 * t1 * x * y + t2 * (r2 + 2.0 * x * x)
        dy = t1 * (r2 + 2.0 * y * y) + 2.0 * t2 * x * y
        x = (xd - dx) / radial
        y = (yd - dy) / radial
    return x * fx + c0U, y * fy + c0V


def test_station_lenses():
    print("the lens guards on the station's lens files: nothing removed, bit-identical")
    files = sorted((HERE / "calibration").glob("*_IO.yaml"))
    check(len(files) >= 8, f"{len(files)} lens files in calibration/")
    U, V = grid(4)
    U = np.r_[U, NU - 1.0, 0.0, NU - 1.0]                       # the frame's far corners too
    V = np.r_[V, 0.0, NV - 1.0, NV - 1.0]
    for p in files:
        io = G.load_intrinsics(p)
        x, y = G.undistort_uv(U, V, io)
        x0, y0 = unguarded(U, V, io)
        ud, vd = VR.distort_uv(x, y, io)
        res = float(np.max(np.hypot(ud - U, vd - V)))
        fold = G.lens_fold_r2(io)
        r2 = float(np.max(((x - io[2]) / io[4]) ** 2 + ((y - io[3]) / io[5]) ** 2))
        check(np.isfinite(x).all() and np.array_equal(x, x0) and np.array_equal(y, y0) and res < 1e-9
              and r2 < fold,
              f"{p.name}: every pixel inverted, the same bits as without the guards, residual {res:.0e} px, "
              f"corner r^2 {r2:.3f} < fold {fold:.3f}")


def test_folding_lens():
    print("the lens guards on a synthetic lens that folds inside the frame")
    io = PINHOLE.copy()
    io[6] = -0.55                                    # d1: dr_d/dr = 1 - 1.65 r^2 folds at r^2 = 1/1.65
    fold = G.lens_fold_r2(io)
    check(abs(fold - 1.0 / 1.65) < 1e-12, f"lens_fold_r2 {fold:.6f} = 1/(3 x 0.55)")
    no_lens = G.lens_fold_r2(PINHOLE)
    io2 = PINHOLE.copy()
    io2[6], io2[7], io2[8] = 0.1, 0.0, -0.02          # folds only through d3: 1 + 0.3 s - 0.14 s^3 = 0
    s = G.lens_fold_r2(io2)
    check(no_lens == np.inf and abs(1 + 0.3 * s - 0.14 * s ** 3) < 1e-9 and s > 0,
          f"no distortion: no fold (inf); d1 0.1, d3 -0.02: fold at r^2 {s:.4f}, a root of 1 + 0.3 s - 0.14 s^3")
    U, V = grid(8)
    # the guards add no numpy warning to those of the passes themselves (which could overflow beyond the fold)
    with warnings.catch_warnings(record=True) as w_guarded:
        warnings.simplefilter("always")
        x, y = G.undistort_uv(U, V, io)
    with warnings.catch_warnings(record=True) as w_passes:
        warnings.simplefilter("always")
        x0, y0 = unguarded(U, V, io)
    quiet = sorted(str(w.message) for w in w_guarded) == sorted(str(w.message) for w in w_passes)
    with np.errstate(all="ignore"):
        ud, vd = VR.distort_uv(x, y, io)
        ud0, vd0 = VR.distort_uv(x0, y0, io)
    fin = np.isfinite(x)
    rd2 = ((U - io[2]) / io[4]) ** 2 + ((V - io[3]) / io[5]) ** 2
    rd_fold2 = fold * (1 + io[6] * fold) ** 2        # the largest distorted radius the lens reaches
    beyond = rd2 > rd_fold2
    res = np.hypot(ud - U, vd - V)
    res0 = np.hypot(ud0 - U, vd0 - V)
    r2 = ((x - io[2]) / io[4]) ** 2 + ((y - io[3]) / io[5]) ** 2
    garbage = np.isfinite(x0) & ~(res0 <= G.UNDISTORT_TOL_PX)
    check(quiet, f"the guards raise no numpy warning of their own ({len(w_guarded)} from the passes, as without them)")
    check(beyond.sum() > 100 and not fin[beyond].any(),
          f"{int(beyond.sum())} pixels lie beyond the largest radius the lens reaches: all NaN")
    check(fin.sum() > U.size // 2 and np.all(res[fin] <= G.UNDISTORT_TOL_PX) and np.all(r2[fin] < fold),
          f"{int(fin.sum())} of {U.size} pixels kept, each maps back onto its pixel (worst {np.max(res[fin]):.1e} px) "
          f"and lies inside the fold")
    check(garbage.sum() > 100 and not fin[garbage].any(),
          f"without the guards {int(garbage.sum())} pixels came back finite and missed their pixel by more than "
          f"{G.UNDISTORT_TOL_PX} px (up to {np.max(res0[garbage]):.0f} px): now NaN")
    ok_both = fin & np.isfinite(x0)
    check(np.array_equal(x[ok_both], x0[ok_both]) and np.array_equal(y[ok_both], y0[ok_both]),
          "where a pixel is kept, its value is the unguarded one, bit for bit")
    # a lens whose r_d rises again past its fold: there the passes can converge onto that outer branch, onto a
    # point that maps back onto its pixel, which only the fold guard catches
    io3 = PINHOLE.copy()
    io3[6], io3[7] = -0.9, 0.35                      # dr_d/dr = 1 - 2.7 s + 1.75 s^2: 0 at s = 0.618 and 0.925
    fold3 = G.lens_fold_r2(io3)
    with warnings.catch_warnings(record=True) as w_outer:
        warnings.simplefilter("always")
        x3, y3 = G.undistort_uv(U, V, io3)
    with np.errstate(all="ignore"):
        x30, y30 = unguarded(U, V, io3)
        ud30, vd30 = VR.distort_uv(x30, y30, io3)
    r230 = ((x30 - io3[2]) / io3[4]) ** 2 + ((y30 - io3[3]) / io3[5]) ** 2
    outer = np.isfinite(x30) & (np.hypot(ud30 - U, vd30 - V) <= G.UNDISTORT_TOL_PX) & (r230 >= fold3)
    check(abs(fold3 - (2.7 - 0.29 ** 0.5) / 3.5) < 1e-12 and outer.sum() > 100
          and not np.isfinite(x3[outer]).any() and not np.isfinite(y3[outer]).any(),
          f"d1 -0.9, d2 0.35 (fold at r^2 {fold3:.4f}): without the guards {int(outer.sum())} pixels converged "
          f"onto the branch beyond the fold, each within {G.UNDISTORT_TOL_PX} px of its pixel: now NaN")
    check(not w_outer, f"... with no numpy warning ({len(w_outer)} raised)")
    xs, ys = G.undistort_uv(1300.0, 1100.0, PINHOLE)
    # numpy 1.17 (the station's) warns on comparing a NaN outside np.errstate
    with warnings.catch_warnings(record=True) as w_nan:
        warnings.simplefilter("always")
        xn, yn = G.undistort_uv(np.nan, 1100.0, io)
        xg, yg = G.undistort_uv(np.r_[np.nan, 1300.0], np.r_[1100.0, 1100.0], io)
    check(isinstance(xs, float) and abs(xs - 1300.0) < 1e-9 and abs(ys - 1100.0) < 1e-9 and np.isnan(xn)
          and np.isnan(yn) and np.isnan(xg[0]) and np.isnan(yg[0]) and np.isfinite(xg[1]),
          "a scalar pixel gives a scalar back; a NaN pixel gives NaN, alone or beside a good one")
    check(not w_nan, f"a NaN pixel raises no numpy warning ({len(w_nan)} raised)")


# ---------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------

def write_yaml(path, keys, values):
    with open(path, "w") as f:
        for k, v in zip(keys, values):
            f.write(f"{k}: {float(v)!r}\n")


def test_cli(d):
    print("georectify.py on a camera aimed at 150 deg: every waterline point georectified")
    io_p, eo_p = Path(d) / "synthetic_IO.yaml", Path(d) / "synthetic_EO.yaml"
    write_yaml(io_p, G.IO_KEYS, list(PINHOLE))
    write_yaml(eo_p, G.EO_KEYS, [CAM_E, CAM_N, CAM_Z, 150.0, TILT, ROLL])
    cols = np.linspace(100.0, NU - 100.0, 23)
    contour, out = Path(d) / "contour.csv", Path(d) / "ground.csv"
    with open(contour, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["camera", "pixel_column", "pixel_row", "tide_elevation_navd88"])
        for c in cols:
            w.writerow(["c1", f"{c:.3f}", "1300.000", f"{Z_PLANE}"])
    r = subprocess.run([sys.executable, str(HERE / "georectify.py"), str(contour), str(out),
                        "--io-c1", str(io_p), "--eo-c1", str(eo_p)],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    rows = []
    if out.exists():
        with open(out, newline="") as fh:
            rows = list(csv.DictReader(fh))
    E = np.array([float(x["easting_utm19"]) if x["easting_utm19"] else np.nan for x in rows])
    N = np.array([float(x["northing_utm19"]) if x["northing_utm19"] else np.nan for x in rows])
    Er, Nr = cast_rays(cols, np.full(cols.size, 1300.0), Z_PLANE, G.load_extrinsics(eo_p))
    check(r.returncode == 0 and "georectified          : 23" in r.stdout and "no intersection" not in r.stdout,
          "exit 0, 'georectified : 23', none 'at or above the horizon'")
    check(len(rows) == 23 and np.isfinite(E).all() and np.max(np.hypot(E - Er, N - Nr)) < 1e-3,
          "the 23 rows carry the ray cast's easting and northing (to the file's 1 mm)")


def test_py38():
    print("Python 3.8 syntax")
    for f in ("georectify.py", "view_reproject.py", "test_georectify.py"):
        src = (HERE / f).read_text(encoding="utf-8")
        try:
            try:
                ast.parse(src, feature_version=(3, 8))
            except TypeError:                    # Python 3.7 has no feature_version
                ast.parse(src)
            ok = True
        except SyntaxError as e:
            print(f"       {e}")
            ok = False
        check(ok, f"{f} parses as Python 3.8")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="keep the outputs, in a new folder inside this one")
    args = ap.parse_args()
    if args.keep:
        Path(args.keep).mkdir(parents=True, exist_ok=True)
    d = os.path.realpath(tempfile.mkdtemp(prefix="georectify_", dir=args.keep))
    print(f"test folder: {d}")
    try:
        test_every_azimuth()
        test_every_azimuth(STEEP_TILT)
        test_local_frame()
        test_station_lens_every_azimuth()
        test_station_lenses()
        test_folding_lens()
        test_cli(d)
        test_py38()
    finally:
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED")
        return 1
    print("all georectify tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Reproject Imagery Between Two Camera Geometries
==================================================
Makes an image taken with one camera pointing (an older extrinsic
calibration) look as the SAME camera would have seen it with another
pointing (the current calibration), on a chosen horizontal plane.

WHY. The detector's search strips, per-column envelope, step sizes and
bias corrections are all tuned in pixel positions of the 2025-11-13
view. Imagery from before the Nov 2025 camera move (e.g. the Dec 2024 -
Mar 2025 ADCP deployment, EO solved 2025-02-19) sees the beach several
hundred pixels away from where those settings look. Rather than retune
every setting, each old frame is resampled into the current view.

EXACT ON THE WATER SURFACE. The mapping goes through the horizontal
plane z = water level at the time of the frame: for every pixel of the
current view, cast its ray to that plane and read the old image where
the old camera saw that ground point. The waterline lies on that plane,
so it lands exactly where the current camera would have seen it, and
georectifying the detection with the current calibration at that same
z returns the true ground position. Anything off the plane (dune, sky,
wave crests) is displaced, which does not matter for the waterline.

The lens intrinsics (IO, 2024-08-01) are the same for both geometries:
the camera was re-aimed, not re-lensed.

Functions:
    ground_to_pixel(X, Y, Z, io, eo)      -> distorted (U, V), valid
    remap_tables(io, eo_src, eo_dst, z)   -> map_x, map_y for cv2.remap
    reproject_image(img, io, eo_src, eo_dst, z)
    dst_to_src_points(u, v, z, io, eo_src, eo_dst)   (detections back to the original frame)

Self-test (round trip and surveyed control points):
    python3 view_reproject.py --self-test
"""

import argparse
from functools import lru_cache
from pathlib import Path

import numpy as np

from georectify import build_P, load_extrinsics, load_intrinsics, pixel_to_ground

HERE = Path(__file__).resolve().parent
MAP_STEP = 4            # remap tables are computed every MAP_STEP px, then interpolated


def distort_uv(U, V, intrinsics):
    """Ideal -> distorted pixel coordinates (CIRN distortUV.m forward model)."""
    c0U, c0V, fx, fy = intrinsics[2], intrinsics[3], intrinsics[4], intrinsics[5]
    d1, d2, d3, t1, t2 = intrinsics[6], intrinsics[7], intrinsics[8], intrinsics[9], intrinsics[10]
    x = (np.asarray(U, dtype=float) - c0U) / fx
    y = (np.asarray(V, dtype=float) - c0V) / fy
    r2 = x * x + y * y
    radial = 1.0 + d1 * r2 + d2 * r2 * r2 + d3 * r2 * r2 * r2
    dx = 2.0 * t1 * x * y + t2 * (r2 + 2.0 * x * x)
    dy = t1 * (r2 + 2.0 * y * y) + 2.0 * t2 * x * y
    return (x * radial + dx) * fx + c0U, (y * radial + dy) * fy + c0V


def ground_to_pixel(X, Y, Z, intrinsics, extrinsics):
    """
    World (UTM easting, northing, NAVD88 z) -> distorted image (U, V).
    `valid` is False for points behind the camera or outside the frame.
    """
    P = build_P(intrinsics, extrinsics)
    X, Y, Z = np.broadcast_arrays(np.asarray(X, float), np.asarray(Y, float), np.asarray(Z, float))
    s = P[2, 0] * X + P[2, 1] * Y + P[2, 2] * Z + P[2, 3]
    with np.errstate(invalid="ignore", divide="ignore"):
        U = (P[0, 0] * X + P[0, 1] * Y + P[0, 2] * Z + P[0, 3]) / s
        V = (P[1, 0] * X + P[1, 1] * Y + P[1, 2] * Z + P[1, 3]) / s
    Ud, Vd = distort_uv(U, V, intrinsics)
    nu, nv = intrinsics[0], intrinsics[1]
    # s < 0 in front of the camera under this convention (see pixel_to_ground)
    valid = np.isfinite(Ud) & np.isfinite(Vd) & (s < 0) & \
        (Ud >= 0) & (Ud <= nu - 1) & (Vd >= 0) & (Vd <= nv - 1)
    return Ud, Vd, valid


def dst_to_src_points(u, v, z, intrinsics, eo_src, eo_dst):
    """Pixels in the destination (current) view at height z -> pixels in the source (old) view."""
    X, Y = pixel_to_ground(np.asarray(u, float), np.asarray(v, float), z, intrinsics, eo_dst)
    return ground_to_pixel(X, Y, z, intrinsics, eo_src)


@lru_cache(maxsize=8)
def _tables(io_key, src_key, dst_key, z_cm):
    intrinsics, eo_src, eo_dst = np.array(io_key), np.array(src_key), np.array(dst_key)
    nu, nv = int(intrinsics[0]), int(intrinsics[1])
    gu = np.arange(0, nu + MAP_STEP, MAP_STEP, dtype=float)
    gv = np.arange(0, nv + MAP_STEP, MAP_STEP, dtype=float)
    UU, VV = np.meshgrid(gu, gv)
    su, sv, ok = dst_to_src_points(UU.ravel(), VV.ravel(), z_cm / 100.0, intrinsics, eo_src, eo_dst)
    su = np.where(ok, su, -1e4).reshape(UU.shape).astype(np.float32)
    sv = np.where(ok, sv, -1e4).reshape(UU.shape).astype(np.float32)
    return su, sv


def remap_tables(intrinsics, eo_src, eo_dst, z):
    """Full-resolution cv2.remap tables: destination pixel -> source pixel, at plane z."""
    import cv2
    su, sv = _tables(tuple(intrinsics), tuple(eo_src), tuple(eo_dst), int(round(z * 100)))
    nu, nv = int(intrinsics[0]), int(intrinsics[1])
    # The tables are smooth, so a coarse grid interpolated up is exact to
    # a small fraction of a pixel and ~16x faster than every pixel.
    gw, gh = su.shape[1], su.shape[0]
    xs = np.arange(nu, dtype=np.float32) / MAP_STEP
    ys = np.arange(nv, dtype=np.float32) / MAP_STEP
    mx, my = np.meshgrid(xs, ys)
    map_x = cv2.remap(su, mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)
    map_y = cv2.remap(sv, mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)
    assert gw * MAP_STEP >= nu and gh * MAP_STEP >= nv
    return map_x, map_y


def reproject_image(img, intrinsics, eo_src, eo_dst, z):
    """Old-view image -> current-view image on plane z. Unseen areas are black."""
    import cv2
    map_x, map_y = remap_tables(intrinsics, eo_src, eo_dst, z)
    return cv2.remap(img, map_x, map_y, cv2.INTER_LINEAR,
                     borderMode=cv2.BORDER_CONSTANT, borderValue=0)


def self_test():
    cal = HERE / "calibration"
    worst = 0.0
    for cam in ("c1", "c2"):
        io = load_intrinsics(cal / f"CACO05_{cam}_20240801_IO.yaml")
        eo_new = load_extrinsics(cal / f"CACO05_{cam}_20251113_EO-CV.yaml")
        eo_old = load_extrinsics(cal / f"CACO05_{cam}_20250219_EO.yaml")

        # 1. Surveyed control: projecting the targets must land on the picked pixels.
        uv = np.genfromtxt(cal / f"gcp_uv_{cam}.csv", delimiter=",", names=True)
        xyz = np.genfromtxt(cal / f"2025-11-13_Marconi_Extrinsic_Targets_{cam}_xyz.csv", delimiter=",")
        by_num = {int(r[0]): r[1:] for r in xyz}
        d = []
        for n, u, v in zip(uv["num"], uv["U"], uv["V"]):
            if int(n) in by_num:
                pu, pv, ok = ground_to_pixel(*by_num[int(n)], io, eo_new)
                if ok:
                    d.append(np.hypot(pu - u, pv - v))
        print(f"{cam}: control targets reprojected to {len(d)} picked pixels, "
              f"median {np.median(d):.1f} px, max {np.max(d):.1f} px")

        # 2. Round trip pixel -> ground -> pixel, same geometry, at z = 0.5 m.
        U, V = np.meshgrid(np.linspace(100, io[0] - 100, 25), np.linspace(900, io[1] - 100, 20))
        X, Y = pixel_to_ground(U.ravel(), V.ravel(), 0.5, io, eo_new)
        u2, v2, ok = ground_to_pixel(X, Y, 0.5, io, eo_new)
        err = np.hypot(u2 - U.ravel(), v2 - V.ravel())[ok]
        worst = max(worst, float(err.max()))
        print(f"{cam}: round trip over {ok.sum()} pixels, max error {err.max():.3f} px")

        # 3. Identity: same geometry both sides -> tables are the identity.
        mx, my = remap_tables(io, eo_new, eo_new, 0.5)
        sel = (slice(1000, 1900, 50), slice(100, 2300, 50))
        idn = np.hypot(mx[sel] - np.arange(io[0])[None, 100:2300:50],
                       my[sel] - np.arange(io[1])[1000:1900:50, None])
        print(f"{cam}: identity remap, max {np.nanmax(idn):.3f} px")

        # 4. How far the Feb 2025 view is from the Nov 2025 view at the waterline.
        su, sv, ok = dst_to_src_points(U.ravel(), V.ravel(), 0.5, io, eo_old, eo_new)
        sh = np.hypot(su - U.ravel(), sv - V.ravel())[ok]
        print(f"{cam}: Feb 2025 vs Nov 2025 view at z=0.5 m: median shift {np.median(sh):.0f} px "
              f"({ok.sum()} of {ok.size} sample points visible in both)")
    return worst


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--self-test", action="store_true")
    if ap.parse_args().self_test:
        self_test()

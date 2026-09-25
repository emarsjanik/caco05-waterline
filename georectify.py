#!/usr/bin/env python3
"""
Pixel -> Ground Georectification (CIRN convention)
-----------------------------------------------------
Converts detected waterline pixel coordinates into real-world UTM
Zone 19 easting/northing, using the camera intrinsics (IO) and
extrinsics (EO) supplied by USGS Woods Hole.

HOW IT WORKS: a pixel defines a ray from the camera through the scene,
which alone does not fix a 3D point -- one more constraint is needed.
For a waterline that constraint comes free: every point on the line
lies at the water surface, whose elevation is already known from
GNSS-R. So each contour point is intersected with a horizontal plane
at ITS OWN measured water elevation. That is why this works for
shorelines when it would not for arbitrary image features.

CONVENTION: CIRN, as used by the supplied calibration scripts
(C_singleExtrinsicSolution_CACO05_*.m). Extrinsics are
[x y z azimuth tilt swing] in radians internally; the EO-CV YAML files
name the third angle "roll", and the solved value (-8.3432) matches the
YAML's roll (-8.3346) to within the fit uncertainty, so they are
treated as the same quantity. Note K uses -fx, not +fx: that sign is
part of the CIRN convention and flipping it mirrors the scene.

Intrinsics are dated 2024-08-01 while extrinsics are 2025-11-13. That
is deliberate per the supplied notes: the camera was physically moved
in Nov 2025 but not re-lensed, so the intrinsics carry over. IO and EO
therefore need to be selected by DIFFERENT dates.

Usage (validate against surveyed ground control points):
    python3 georectify.py --self-test \\
        --io CACO05_c1_20240801_IO.yaml \\
        --eo CACO05_c1_20251113_EO-CV.yaml \\
        --gcp-uv gcp_uv_c1.csv --gcp-xyz targets_c1_xyz.csv

Usage (convert a contour file):
    python3 georectify.py contour_points_timex.csv out_ground.csv \\
        --io-c1 ... --eo-c1 ... --io-c2 ... --eo-c2 ...
"""

import sys
import csv
import argparse
from pathlib import Path

import numpy as np


IO_KEYS = ["NU", "NV", "coU", "coV", "fx", "fy", "d1", "d2", "d3", "t1", "t2"]
EO_KEYS = ["x", "y", "z", "azimuth", "tilt", "roll"]


def read_simple_yaml(path):
    """Reads the flat 'key: value  # comment' files used for IO and EO."""
    values = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" not in line:
                continue
            key, _, rest = line.partition(":")
            rest = rest.split("#")[0].strip()
            try:
                values[key.strip()] = float(rest)
            except ValueError:
                continue
    return values


def load_intrinsics(path):
    v = read_simple_yaml(path)
    missing = [k for k in IO_KEYS if k not in v]
    if missing:
        raise ValueError(f"{path}: missing intrinsic key(s) {missing}")
    return np.array([v[k] for k in IO_KEYS], dtype=float)


def load_extrinsics(path):
    v = read_simple_yaml(path)
    missing = [k for k in EO_KEYS if k not in v]
    if missing:
        raise ValueError(f"{path}: missing extrinsic key(s) {missing}")
    return np.array([v["x"], v["y"], v["z"],
                     np.deg2rad(v["azimuth"]), np.deg2rad(v["tilt"]),
                     np.deg2rad(v["roll"])], dtype=float)


def cirn_angles_to_R(azimuth, tilt, swing):
    """Rotation matrix, CIRN convention (matches CIRNangles2R.m)."""
    R = np.empty((3, 3))
    R[0, 0] = -np.cos(azimuth) * np.cos(swing) - np.sin(azimuth) * np.cos(tilt) * np.sin(swing)
    R[0, 1] = np.cos(swing) * np.sin(azimuth) - np.sin(swing) * np.cos(tilt) * np.cos(azimuth)
    R[0, 2] = -np.sin(swing) * np.sin(tilt)
    R[1, 0] = -np.sin(swing) * np.cos(azimuth) + np.cos(swing) * np.cos(tilt) * np.sin(azimuth)
    R[1, 1] = np.sin(swing) * np.sin(azimuth) + np.cos(swing) * np.cos(tilt) * np.cos(azimuth)
    R[1, 2] = np.cos(swing) * np.sin(tilt)
    R[2, 0] = np.sin(tilt) * np.sin(azimuth)
    R[2, 1] = np.sin(tilt) * np.cos(azimuth)
    R[2, 2] = -np.cos(tilt)
    return R


def build_P(intrinsics, extrinsics):
    """Camera projection matrix P = K R [I | -C], CIRN convention."""
    c0U, c0V, fx, fy = intrinsics[2], intrinsics[3], intrinsics[4], intrinsics[5]
    # Both focal terms are NEGATIVE. Verified against the surveyed
    # ground control: with +fy the projected V came out mirrored about
    # the principal point (V_projected + V_picked = 2048 = image
    # height, on every single GCP) while U matched to a few pixels.
    # That signature is a flipped image row axis, not a bad solve.
    K = np.array([[-fx, 0.0, c0U],
                  [0.0, -fy, c0V],
                  [0.0, 0.0, 1.0]])
    R = cirn_angles_to_R(extrinsics[3], extrinsics[4], extrinsics[5])
    C = extrinsics[0:3]
    IC = np.hstack([np.eye(3), -C.reshape(3, 1)])
    P = K @ R @ IC
    return P / P[2, 3]


def undistort_uv(Ud, Vd, intrinsics, iterations=20):
    """
    Removes lens distortion from image coordinates.

    The calibration provides the FORWARD model (ideal -> distorted, per
    distortUV.m). Recovering ideal coordinates requires inverting it,
    which has no closed form, so this iterates: start from the distorted
    point and repeatedly subtract the distortion implied by the current
    estimate. Converges in a few passes for normal lens coefficients.
    """
    c0U, c0V, fx, fy = intrinsics[2], intrinsics[3], intrinsics[4], intrinsics[5]
    d1, d2, d3, t1, t2 = intrinsics[6], intrinsics[7], intrinsics[8], intrinsics[9], intrinsics[10]

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


def pixel_to_ground(Ud, Vd, Z, intrinsics, extrinsics, max_range_m=1000.0):
    """
    Intersects the ray through each pixel with the horizontal plane at
    elevation Z, returning (easting, northing).

    Z may be a scalar or per-point, which is what lets each waterline
    contour use its own measured water elevation.

    Points whose ray is near-parallel to the plane (looking at or above
    the horizon) are returned as NaN rather than a huge extrapolated
    coordinate: the intersection is real but meaningless, and silently
    emitting a value kilometres offshore would be worse than a gap.
    """
    P = build_P(intrinsics, extrinsics)
    U, V = undistort_uv(Ud, Vd, intrinsics)
    U = np.atleast_1d(U); V = np.atleast_1d(V)
    Z = np.broadcast_to(np.atleast_1d(np.asarray(Z, dtype=float)), U.shape)

    # s[U,V,1]' = P[X,Y,Z,1]'  ->  eliminate s, solve 2x2 for X,Y.
    a11 = U * P[2, 0] - P[0, 0]
    a12 = U * P[2, 1] - P[0, 1]
    a21 = V * P[2, 0] - P[1, 0]
    a22 = V * P[2, 1] - P[1, 1]
    b1 = (P[0, 2] * Z + P[0, 3]) - U * (P[2, 2] * Z + P[2, 3])
    b2 = (P[1, 2] * Z + P[1, 3]) - V * (P[2, 2] * Z + P[2, 3])

    det = a11 * a22 - a12 * a21
    bad = np.abs(det) < 1e-12
    safe = np.where(bad, 1.0, det)
    X = (b1 * a22 - b2 * a12) / safe
    Y = (a11 * b2 - a21 * b1) / safe
    X = np.where(bad, np.nan, X)
    Y = np.where(bad, np.nan, Y)

    # Reject rays that do not reach the plane in front of the camera.
    #
    # The scale factor is NEGATIVE for valid in-view points under this
    # convention, because both focal terms in K are negative. Measured
    # on the surveyed ground control, every in-view point gave
    # s between -1.2e-5 and -2.6e-5. An earlier version tested s > 0
    # and rejected every point in the dataset -- the sign here is not
    # something to reason about, it is something to check.
    s = P[2, 0] * X + P[2, 1] * Y + P[2, 2] * Z + P[2, 3]
    invalid = ~np.isfinite(s) | (s >= 0)

    # Also drop absurd ranges. Near the horizon the ray becomes almost
    # parallel to the water plane, so a pixel one row higher can move
    # the intersection by kilometres. Such a point is arithmetically
    # real but physically meaningless.
    if max_range_m and max_range_m > 0:
        C = extrinsics[0:3]
        rng = np.hypot(X - C[0], Y - C[1])
        invalid = invalid | ~np.isfinite(rng) | (rng > max_range_m)

    X = np.where(invalid, np.nan, X)
    Y = np.where(invalid, np.nan, Y)
    return X, Y


def self_test(io_path, eo_path, gcp_uv_path, gcp_xyz_path):
    """
    Projects surveyed ground control points from their picked image
    coordinates and compares against their surveyed positions.

    This is the only honest check of the whole chain -- convention,
    sign of K, distortion inversion, angle units. Synthetic data would
    only confirm the code matches itself.
    """
    intrinsics = load_intrinsics(io_path)
    extrinsics = load_extrinsics(eo_path)

    xyz = {}
    with open(gcp_xyz_path, "r") as f:
        for row in csv.reader(f):
            if len(row) < 4:
                continue
            try:
                xyz[int(float(row[0]))] = (float(row[1]), float(row[2]), float(row[3]))
            except ValueError:
                continue

    uv = {}
    with open(gcp_uv_path, "r") as f:
        for row in csv.DictReader(f):
            uv[int(float(row["num"]))] = (float(row["U"]), float(row["V"]))

    shared = sorted(set(xyz) & set(uv))
    if not shared:
        print("No ground control points in common between the UV and XYZ files.")
        return 1

    print("=" * 72)
    print(f"GEORECTIFICATION SELF-TEST  ({Path(io_path).name} + {Path(eo_path).name})")
    print("=" * 72)
    print(f"{'GCP':>6} {'east err':>10} {'north err':>10} {'dist':>8}")

    errors = []
    for num in shared:
        E, N, Z = xyz[num]
        U, V = uv[num]
        Xe, Yn = pixel_to_ground(U, V, Z, intrinsics, extrinsics)
        de, dn = float(Xe[0]) - E, float(Yn[0]) - N
        dist = float(np.hypot(de, dn))
        errors.append(dist)
        print(f"{num:>6} {de:+10.3f} {dn:+10.3f} {dist:8.3f}")

    errors = np.array(errors)
    print("=" * 72)
    print(f"n={len(errors)}   mean {errors.mean():.3f} m   median "
          f"{np.median(errors):.3f} m   max {errors.max():.3f} m   "
          f"RMS {np.sqrt((errors**2).mean()):.3f} m")
    print()
    if errors.mean() < 1.0:
        print("PASS -- consistent with the reprojection errors in the supplied")
        print("calibration notes (mostly < 0.3 m, one GCP at 1.17 m).")
    else:
        print("FAIL -- errors far exceed the supplied calibration's own residuals.")
        print("Something in the chain is wrong: angle convention, the sign of fx")
        print("in K, degrees vs radians, or swing/roll not being equivalent.")
        print("Do NOT georectify production data until this passes.")
    return 0 if errors.mean() < 1.0 else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("contour_csv", nargs="?")
    ap.add_argument("output_csv", nargs="?")
    ap.add_argument("--io-c1"); ap.add_argument("--eo-c1")
    ap.add_argument("--io-c2"); ap.add_argument("--eo-c2")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--io"); ap.add_argument("--eo")
    ap.add_argument("--gcp-uv"); ap.add_argument("--gcp-xyz")
    args = ap.parse_args()

    if args.self_test:
        sys.exit(self_test(args.io, args.eo, args.gcp_uv, args.gcp_xyz))

    if not args.contour_csv or not args.output_csv:
        ap.error("contour_csv and output_csv are required unless --self-test is given")

    cams = {}
    if args.io_c1 and args.eo_c1:
        cams["c1"] = (load_intrinsics(args.io_c1), load_extrinsics(args.eo_c1))
    if args.io_c2 and args.eo_c2:
        cams["c2"] = (load_intrinsics(args.io_c2), load_extrinsics(args.eo_c2))
    if not cams:
        ap.error("supply at least one camera's --io-cN and --eo-cN")

    rows = list(csv.DictReader(open(args.contour_csv)))
    if not rows:
        print("Empty contour file."); sys.exit(1)

    fieldnames = list(rows[0].keys()) + ["easting_utm19", "northing_utm19"]
    written = skipped_cam = skipped_nan = 0

    by_cam = {}
    for r in rows:
        by_cam.setdefault(r["camera"], []).append(r)

    out = []
    for cam, group in by_cam.items():
        if cam not in cams:
            skipped_cam += len(group)
            print(f"WARNING: no calibration supplied for '{cam}' -- {len(group)} point(s) "
                  f"written without ground coordinates.")
            for r in group:
                r["easting_utm19"] = ""; r["northing_utm19"] = ""
                out.append(r)
            continue
        intr, extr = cams[cam]
        U = np.array([float(r["pixel_column"]) for r in group])
        V = np.array([float(r["pixel_row"]) for r in group])
        # The plane each waterline point lies on: the beach elevation
        # (water level + wave setup) when extract_elevation_contours.py
        # was run with --setup-coef, otherwise the water level.
        Z = np.array([float(r.get("beach_elevation_navd88") or r["tide_elevation_navd88"])
                      for r in group])
        E, N = pixel_to_ground(U, V, Z, intr, extr)
        for r, e, n in zip(group, E, N):
            if np.isfinite(e) and np.isfinite(n):
                r["easting_utm19"] = f"{e:.3f}"; r["northing_utm19"] = f"{n:.3f}"
                written += 1
            else:
                r["easting_utm19"] = ""; r["northing_utm19"] = ""
                skipped_nan += 1
            out.append(r)

    with open(args.output_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader(); w.writerows(out)

    print(f"Wrote {len(out)} row(s) to {args.output_csv}")
    print(f"  georectified          : {written}")
    if skipped_nan:
        print(f"  no intersection       : {skipped_nan}  (ray at or above the horizon)")
    if skipped_cam:
        print(f"  no calibration for cam: {skipped_cam}")

    if written:
        E = np.array([float(r["easting_utm19"]) for r in out if r["easting_utm19"]])
        N = np.array([float(r["northing_utm19"]) for r in out if r["northing_utm19"]])
        print(f"  easting  range: {E.min():.1f} to {E.max():.1f} m")
        print(f"  northing range: {N.min():.1f} to {N.max():.1f} m")
        print(f"  extent: {E.max()-E.min():.1f} m E-W by {N.max()-N.min():.1f} m N-S")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Camera Pointing Of A Re-Aimed Camera, From Its Photos
========================================================
Solves the extrinsics (EO) of a camera setup that has no calibration of
its own -- e.g. the CACO03 frames of Jan 2025 -- from photos of it and
photos of a setup that HAS one (CACO04, EO 2025-02-19).

WHY THIS WORKS. The camera stayed on the same mount; only its pointing
changed. For a camera that only rotates, every scene point -- near or
far, on the plane or not -- moves between the two photos by the same
3-D rotation of its viewing ray. So matched features (dune edge, bluff,
fence posts, the horizon) give that rotation Q exactly, and

    R_new = Q R_ref        (position, lens and height unchanged)

HOW. SIFT features on land (sky and moving water masked out with the
reference geometry) are matched between each new-setup photo and each
reference photo; the rays are undistorted with the lens model; a
rotation is fitted with RANSAC (3 matches per trial, Kabsch), then
refined on all inliers pooled over every pair.

CHECKS (printed). Inlier count and ray residual in pixels; and the
HORIZON, independent of the matches: where the new EO puts the sea
horizon against where each new-setup photo shows it, across the frame.
A good solve agrees within a few pixels; the reference EO applied to the
same photo shows how far off it was.

Usage (dates as printed by horizon_check.py):
    python3 estimate_eo_rotation.py --camera c1 \\
        --images-dir /mnt/I2Rgus_Data/Chelsea_calibration/work/original \\
        --first 2025-01-18 --last 2025-01-22 --add-to-setups
    -> calibration/CACO05_c1_2025-01-18_to_2025-01-22_EO.yaml, and a line in
       calibration/chelsea_setups.csv so process_chelsea.py uses it for those days.
    Reference photos are the other days (not in any listed setup).
    Or give photos directly: --new A.jpg B.jpg --ref C.jpg D.jpg --output X.yaml
    (up to --max-images photos of each, one per day nearest 17:00 UTC, are used)
"""

import re
import sys
import argparse
from pathlib import Path
from datetime import datetime, timezone

import numpy as np

from georectify import cirn_angles_to_R, load_extrinsics, load_intrinsics, undistort_uv
from view_reproject import ground_to_pixel

HERE = Path(__file__).resolve().parent
CAL = HERE / "calibration"
EARTH_R = 6371000.0
SCALE = 0.5            # features are found at half size


def K_of(io):
    # Same (negative-focal) convention as georectify.build_P.
    return np.array([[-io[4], 0.0, io[2]], [0.0, -io[5], io[3]], [0.0, 0.0, 1.0]])


def rays(u, v, io):
    """Distorted pixels -> unit viewing rays in the camera frame."""
    x, y = undistort_uv(u, v, io)
    r = np.linalg.solve(K_of(io), np.vstack([x, y, np.ones_like(x)]))
    return (r / np.linalg.norm(r, axis=0)).T


def kabsch(a, b):
    """Rotation Q minimising |b - Q a| for row-vector sets a, b."""
    H = a.T @ b
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def R_to_angles(R):
    """Inverse of georectify.cirn_angles_to_R: (azimuth, tilt, roll) in radians."""
    tilt = np.arccos(np.clip(-R[2, 2], -1, 1))
    az = np.arctan2(R[2, 0], R[2, 1])
    roll = np.arctan2(-R[0, 2], R[1, 2])
    return az, tilt, roll


def horizon_rows(io, eo, cols):
    """Rows of the sea horizon (with Earth curvature) at the given columns, or NaN."""
    from georectify import build_P
    h = eo[2]
    d = np.sqrt(2 * EARTH_R * h)
    az = eo[3] + np.deg2rad(np.linspace(-80, 80, 6000))
    X, Y = eo[0] + d * np.sin(az), eo[1] + d * np.cos(az)
    U, V, ok = ground_to_pixel(X, Y, 0.0, io, eo)
    # Keep only points whose UNDISTORTED position is near the frame: far outside
    # it the polynomial lens model folds back and puts them inside the image.
    P = build_P(io, eo)
    q = P @ np.vstack([X, Y, np.zeros_like(X), np.ones_like(X)])
    ui, vi = q[0] / q[2], q[1] / q[2]
    ok &= (ui > -0.1 * io[0]) & (ui < 1.1 * io[0]) & (vi > -0.1 * io[1]) & (vi < 1.1 * io[1])
    if ok.sum() < 2:
        return np.full(len(cols), np.nan)
    o = np.argsort(U[ok])
    return np.interp(cols, U[ok][o], V[ok][o], left=np.nan, right=np.nan)


def observed_horizon(gray, cols, near_rows, window=200):
    """Strongest horizontal edge within +/-window px of an expected row, per column."""
    import cv2
    g = cv2.GaussianBlur(gray.astype(np.float32), (0, 0), 2.5)
    gy = np.abs(cv2.Sobel(g, cv2.CV_32F, 0, 1, ksize=5))
    out = []
    for c, r0 in zip(cols, near_rows):
        if not np.isfinite(r0):
            out.append(np.nan)
            continue
        a, b = int(max(0, r0 - window)), int(min(gray.shape[0], r0 + window))
        prof = gy[a:b, max(0, c - 6):c + 6].mean(axis=1)
        out.append(a + float(np.argmax(prof)) if prof.size and prof.max() > 0 else np.nan)
    return np.array(out)


def land_mask(shape, io, eo, margin=40):
    """Pixels below the reference horizon + margin (no sky, whose clouds move)."""
    rows = horizon_rows(io, eo, np.arange(shape[1]))
    rows = np.where(np.isfinite(rows), rows, 0)
    yy = np.arange(shape[0])[:, None]
    return (yy > rows[None, :] + margin).astype(np.uint8) * 255


def pick(paths, n):
    """Up to n photos, one per day, nearest 17:00 UTC."""
    by_day = {}
    for p in paths:
        m = re.match(r"^(\d{9,11})\.", Path(p).name)
        if not m:
            continue
        t = datetime.fromtimestamp(int(m.group(1)), tz=timezone.utc)
        score = abs(t.hour + t.minute / 60 - 17)
        k = t.date()
        if k not in by_day or score < by_day[k][0]:
            by_day[k] = (score, p)
    days = sorted(by_day)
    if len(days) > n:                              # spread over the period
        days = [days[i] for i in np.linspace(0, len(days) - 1, n).round().astype(int)]
    return [by_day[d][1] for d in days]


def features(path, mask):
    import cv2
    img = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
    if img is None:
        return None, None, None
    small = cv2.resize(img, None, fx=SCALE, fy=SCALE, interpolation=cv2.INTER_AREA)
    m = cv2.resize(mask, (small.shape[1], small.shape[0]), interpolation=cv2.INTER_NEAREST)
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
    kp, des = cv2.SIFT_create(nfeatures=6000).detectAndCompute(clahe.apply(small), m)
    pts = np.array([k.pt for k in kp], float) / SCALE if kp else np.empty((0, 2))
    return img, pts, des


def match(des_a, des_b, ratio=0.75):
    import cv2
    if des_a is None or des_b is None or len(des_a) < 3 or len(des_b) < 3:
        return []
    knn = cv2.BFMatcher(cv2.NORM_L2).knnMatch(des_a, des_b, k=2)
    return [(m.queryIdx, m.trainIdx) for m, *rest in knn if rest and m.distance < ratio * rest[0].distance]


def ransac_rotation(a, b, thresh_rad, iters=3000, seed=0):
    rng = np.random.default_rng(seed)
    best = np.zeros(len(a), bool)
    for _ in range(iters):
        i = rng.choice(len(a), 3, replace=False)
        Q = kabsch(a[i], b[i])
        inl = np.arccos(np.clip(np.sum((a @ Q.T) * b, axis=1), -1, 1)) < thresh_rad
        if inl.sum() > best.sum():
            best = inl
    Q = kabsch(a[best], b[best]) if best.sum() >= 3 else None
    return Q, best


def write_eo(path, eo, angles_deg, note):
    Path(path).write_text(
        f"x: {eo[0]:.4f}\ny: {eo[1]:.4f}\nz: {eo[2]:.4f}\n"
        f"azimuth: {angles_deg[0]:.4f}\ntilt: {angles_deg[1]:.4f}\nroll: {angles_deg[2]:.4f}\n"
        f"# {note}\n"
        "#x, y - camera Easting / Northing (m, UTM zone 19); z - camera elevation (m, NAVD88)\n"
        "#azimuth, tilt, roll - degrees, CIRN convention (see the 2025-02-19 EO files)\n")


def main():
    import cv2
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--camera", required=True, choices=["c1", "c2"])
    ap.add_argument("--images-dir", help="folder of *.timex.jpg; with --first/--last")
    ap.add_argument("--first", help="first UTC date (YYYY-MM-DD) of the setup to solve")
    ap.add_argument("--last", help="last UTC date of the setup to solve")
    ap.add_argument("--add-to-setups", action="store_true",
                    help="append the solved EO to calibration/chelsea_setups.csv for --first..--last "
                         "(only if the horizon check passes)")
    ap.add_argument("--new", nargs="+", help="photos of the setup to solve")
    ap.add_argument("--ref", nargs="+", help="photos of the calibrated setup")
    ap.add_argument("--ref-eo", default=None, help="EO of --ref (default calibration/CACO05_<cam>_20250219_EO.yaml)")
    ap.add_argument("--io", default=None, help="default calibration/CACO05_<cam>_20240801_IO.yaml")
    ap.add_argument("--output", default=None,
                    help="default calibration/CACO05_<cam>_<first>_to_<last>_EO.yaml")
    ap.add_argument("--max-images", type=int, default=5)
    ap.add_argument("--thresh-px", type=float, default=3.0, help="RANSAC inlier threshold (px)")
    args = ap.parse_args()

    if args.images_dir:
        if not (args.first and args.last):
            ap.error("--images-dir needs --first and --last")
        setups = []
        sc = CAL / "chelsea_setups.csv"
        if sc.exists():
            import csv
            with open(sc, newline="") as f:
                setups = [r for r in csv.DictReader(f) if r.get("camera") == args.camera]
        new, ref = [], []
        for p in sorted(Path(args.images_dir).glob(f"*.{args.camera}.timex.jpg")):
            m = re.match(r"^(\d{9,11})\.", p.name)
            if not m:
                continue
            day = datetime.fromtimestamp(int(m.group(1)), tz=timezone.utc).strftime("%Y-%m-%d")
            if args.first <= day <= args.last:
                new.append(p)
            elif not any(r["first_date"] <= day <= r["last_date"] for r in setups):
                ref.append(p)
        args.new, args.ref = args.new or new, args.ref or ref
        args.output = args.output or str(CAL / f"CACO05_{args.camera}_{args.first}_to_{args.last}_EO.yaml")
    if not args.new or not args.ref or not args.output:
        ap.error("give --images-dir/--first/--last, or --new, --ref and --output")

    io = load_intrinsics(args.io or CAL / f"CACO05_{args.camera}_20240801_IO.yaml")
    ref_eo_path = args.ref_eo or CAL / f"CACO05_{args.camera}_20250219_EO.yaml"
    eo_ref = load_extrinsics(ref_eo_path)
    R_ref = cirn_angles_to_R(*eo_ref[3:6])
    new, ref = pick(args.new, args.max_images), pick(args.ref, args.max_images)
    if not new or not ref:
        sys.exit("Need at least one photo of each setup (file names must start with the epoch).")
    print(f"{len(new)} new-setup and {len(ref)} reference photo(s):")
    for p in new + ref:
        print("   ", Path(p).name)

    shape = (int(io[1]), int(io[0]))
    ref_mask = land_mask(shape, io, eo_ref)
    # The new setup is rotated by an unknown amount: keep sky out with a looser margin.
    new_mask = land_mask(shape, io, eo_ref, margin=-250)
    A, B = [], []
    feats_ref = [features(p, ref_mask) for p in ref]
    for p in new:
        _, pn, dn = features(p, new_mask)
        for _, pr, dr in feats_ref:
            for i, j in match(dr, dn):
                A.append(pr[i]); B.append(pn[j])
    if len(A) < 20:
        sys.exit(f"Only {len(A)} feature matches -- too few (snow, fog, or very different light). "
                 f"Try other photos, e.g. clear midday frames close in date to the setup change.")
    A, B = np.array(A), np.array(B)
    ra, rb = rays(A[:, 0], A[:, 1], io), rays(B[:, 0], B[:, 1], io)
    thresh = args.thresh_px / io[4]
    Q, inl = ransac_rotation(ra, rb, thresh)
    if Q is None or inl.sum() < 15:
        sys.exit(f"No consistent rotation ({int(inl.sum())} inliers of {len(A)} matches).")
    for _ in range(3):                                  # refine on all inliers
        err = np.arccos(np.clip(np.sum((ra @ Q.T) * rb, axis=1), -1, 1))
        inl = err < thresh
        Q = kabsch(ra[inl], rb[inl])
    err_px = np.arccos(np.clip(np.sum((ra[inl] @ Q.T) * rb[inl], axis=1), -1, 1)) * io[4]
    angle = np.degrees(np.arccos(np.clip((np.trace(Q) - 1) / 2, -1, 1)))

    R_new = Q @ R_ref
    az, tilt, roll = R_to_angles(R_new)
    assert np.allclose(cirn_angles_to_R(az, tilt, roll), R_new, atol=1e-6)
    eo_new = np.r_[eo_ref[:3], az, tilt, roll]
    ang = np.degrees([az, tilt, roll])
    ref_ang = np.degrees(eo_ref[3:6])

    print(f"\nMatches {len(A)}, inliers {int(inl.sum())} ({100 * inl.mean():.0f}%), "
          f"residual median {np.median(err_px):.2f} px, P90 {np.percentile(err_px, 90):.2f} px")
    print(f"Rotation between the setups: {angle:.2f} deg")
    print(f"             reference ({Path(ref_eo_path).name})   new")
    for k, a0, a1 in zip(("azimuth", "tilt", "roll"), ref_ang, ang):
        print(f"  {k:8s}   {a0:9.4f}   {a1:9.4f}   ({a1 - a0:+.3f})")

    # Independent check: the sea horizon in the new-setup photos.
    cols = np.arange(150, shape[1] - 150, 200)
    pred_new, pred_ref = horizon_rows(io, eo_new, cols), horizon_rows(io, eo_ref, cols)
    print("\nHorizon check (row of the sea horizon; observed = strongest edge near the new EO's row)")
    print("  column   observed   new EO   reference EO")
    diffs_new, diffs_ref = [], []
    for p in new:
        g = cv2.imread(str(p), cv2.IMREAD_GRAYSCALE)
        obs = observed_horizon(g, cols, pred_new)
        diffs_new.append(obs - pred_new); diffs_ref.append(obs - pred_ref)
    obs_med = np.nanmedian(np.array(diffs_new), axis=0) + pred_new
    for c, o, pn_, pr_ in zip(cols, obs_med, pred_new, pred_ref):
        if np.isfinite(o):
            print(f"  {c:6d}   {o:8.0f}   {pn_:6.0f}   {pr_:6.0f}")
    dn, dr = np.abs(np.array(diffs_new)), np.abs(np.array(diffs_ref))
    print(f"  median |observed - EO|: new {np.nanmedian(dn):.1f} px, reference {np.nanmedian(dr):.1f} px")
    horizon_ok = np.nanmedian(dn) <= 10
    if not horizon_ok:
        print("  WARNING: the new EO does not put the horizon where the photos show it; "
              "do not use it (check the photos: fog or haze can hide the horizon).")

    write_eo(args.output, eo_new, ang,
             f"Solved by estimate_eo_rotation.py: {Path(ref_eo_path).name} rotated by {angle:.2f} deg "
             f"from {int(inl.sum())} feature matches ({len(new)} new / {len(ref)} reference photos); "
             f"position unchanged.")
    print(f"\nWrote {args.output}")
    if args.add_to_setups:
        if not (args.first and args.last):
            print("--add-to-setups needs --first and --last; chelsea_setups.csv not changed.")
        elif not horizon_ok:
            print("Horizon check failed: chelsea_setups.csv NOT changed.")
        else:
            import csv
            sc = CAL / "chelsea_setups.csv"
            rows = []
            if sc.exists():
                with open(sc, newline="") as f:
                    rows = [r for r in csv.DictReader(f)
                            if not (r["camera"] == args.camera and r["first_date"] == args.first
                                    and r["last_date"] == args.last)]
            rows.append({"camera": args.camera, "first_date": args.first, "last_date": args.last,
                         "eo_file": Path(args.output).name})
            with open(sc, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=["camera", "first_date", "last_date", "eo_file"])
                w.writeheader(); w.writerows(rows)
            print(f"Added {args.camera} {args.first}..{args.last} -> {Path(args.output).name} "
                  f"to {sc}. process_chelsea.py will use it for those days.")


if __name__ == "__main__":
    main()

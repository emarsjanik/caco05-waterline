#!/usr/bin/env python3
"""
Camera Pointing Monitor: Has The Camera Turned?
=================================================
Measures each camera's pointing (azimuth, tilt, roll) on every new image
by matching it against reference frames whose pointing is known, and
logs the change. A ~22 deg re-aim (24 Jan 2025) went unnoticed for months
and put the Jan 2025 DEM 1.3 m off; at 300 m range on this beach a 0.1 deg pan
error is already ~4 cm of DEM elevation and 0.5 deg is ~21 cm. The
horizon checks (horizon_check.py) see tilt and roll only -- this sees pan.

HOW. The camera stays where it is and can only rotate, so the motion
between a reference frame and a new one is a pure rotation, whatever the
depth of each feature. Features (SIFT, else ORB) are matched between the
two frames; each match gives a viewing ray in the new frame and, through
the reference's known pointing, the same ray in the world. The rotation
that best maps world rays onto new-frame rays (Kabsch, inside RANSAC so
moving water and sand are rejected) IS the new pointing. Lens distortion
is removed first, so the answer is exact rather than a flat-image
approximation.

Timex images are used by default: the sea is blurred, so the features
that match are the bluff, dune and horizon, which do not move.

REFERENCES. A bank per camera in calibration/pointing/<cam>/: each
reference is an image plus a .json with the pointing (EO) it was taken
at. Several references (morning, afternoon, overcast, winter) make the
matching robust to light; each new frame uses the one it matches best.
Start a new bank after every re-aim / recalibration.

    python3 pointing_check.py --camera c2 --add-reference <image.jpg> \\
        --eo calibration/CACO05_c2_20251113_EO-CV.yaml

LOG. archive/pointing_<cam>.csv, one row per image:
    epoch, filename, reference, matches, inliers, residual_deg,
    d_azimuth, d_tilt, d_roll (deg, new minus reference), status
status: ok (|change| <= --warn), moved (> --warn), MOVED (> --fail),
or unmatched (too few features: fog, night, rain on the lens).

Usage:
    python3 pointing_check.py --camera c2                 new images on disk
    python3 pointing_check.py --camera c2 --image x.jpg   one image, printed
    python3 pointing_check.py --self-test --image x.jpg   rotate a frame by known
                                                          angles and recover them
"""

import re
import sys
import json
import argparse
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import cv2

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from georectify import load_intrinsics, load_extrinsics, cirn_angles_to_R, undistort_uv   # noqa: E402

CAL = HERE / "calibration"
FIELDS = ["epoch", "time_utc", "filename", "reference", "matches", "inliers", "residual_deg",
          "d_azimuth", "d_tilt", "d_roll", "status"]


# ------------------------------------------------------------- geometry

def K_inv(io):
    c0U, c0V, fx, fy = io[2], io[3], io[4], io[5]
    K = np.array([[-fx, 0.0, c0U], [0.0, -fy, c0V], [0.0, 0.0, 1.0]])   # georectify.build_P
    return np.linalg.inv(K)


def rays(u, v, io):
    """Distorted pixels -> unit viewing rays in the camera frame (CIRN sign convention)."""
    uu, vv = undistort_uv(u, v, io)
    x = np.vstack([uu, vv, np.ones_like(uu)])
    r = K_inv(io) @ x
    return (r / np.linalg.norm(r, axis=0)).T


def R_to_angles(R):
    """Inverse of georectify.cirn_angles_to_R -> (azimuth, tilt, roll) in radians."""
    tilt = np.arccos(np.clip(-R[2, 2], -1, 1))
    az = np.arctan2(R[2, 0], R[2, 1])
    roll = np.arctan2(-R[0, 2], R[1, 2])
    return az, tilt, roll


def kabsch(world, cam):
    """Rotation R minimising sum |cam_i - R world_i|^2 (both n x 3 unit vectors)."""
    H = world.T @ cam
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1, 1, d]) @ U.T


def ang_err(R, world, cam):
    c = np.clip(np.sum((world @ R.T) * cam, axis=1), -1, 1)
    return np.degrees(np.arccos(c))


def solve_rotation(world, cam, thresh_deg=0.05, iters=2000, rng=None):
    """RANSAC over 2-ray samples, then a least-squares refit on the inliers."""
    rng = rng or np.random.default_rng(0)
    n = len(world)
    best = None
    for _ in range(iters):
        i = rng.choice(n, 2, replace=False)
        R = kabsch(world[i], cam[i])
        inl = ang_err(R, world, cam) < thresh_deg
        if best is None or inl.sum() > best.sum():
            best = inl
    for _ in range(3):                                    # refit and re-select
        R = kabsch(world[best], cam[best])
        best = ang_err(R, world, cam) < thresh_deg
    return R, best


# ------------------------------------------------------------- matching

def prep(img):
    g = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY) if img.ndim == 3 else img
    return cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8)).apply(g)


def detector():
    if hasattr(cv2, "SIFT_create"):
        return cv2.SIFT_create(nfeatures=6000), cv2.NORM_L2
    return cv2.ORB_create(nfeatures=10000), cv2.NORM_HAMMING


def features(img, mask=None):
    det, norm = detector()
    kp, des = det.detectAndCompute(prep(img), mask)
    return kp, des, norm


def match(f1, f2, ratio=0.75):
    (kp1, d1, norm), (kp2, d2, _) = f1, f2
    if d1 is None or d2 is None or len(kp1) < 10 or len(kp2) < 10:
        return np.empty((0, 2)), np.empty((0, 2))
    bf = cv2.BFMatcher(norm)
    good = []
    for m in bf.knnMatch(d1, d2, k=2):
        if len(m) == 2 and m[0].distance < ratio * m[1].distance:
            good.append(m[0])
    p1 = np.array([kp1[m.queryIdx].pt for m in good]).reshape(-1, 2)
    p2 = np.array([kp2[m.trainIdx].pt for m in good]).reshape(-1, 2)
    return p1, p2


class Reference:
    def __init__(self, img_path):
        self.path = Path(img_path)
        meta = json.loads(self.path.with_suffix(".json").read_text())
        self.eo = np.array(meta["eo"], float)                  # x y z az tilt roll (rad)
        self.eo_name = meta.get("eo_file", "")
        img = cv2.imread(str(self.path))
        mask = None
        mp = self.path.with_name(self.path.stem + "_mask.png")
        if mp.exists():
            mask = cv2.imread(str(mp), cv2.IMREAD_GRAYSCALE)
        self.feat = features(img, mask)
        self.R = cirn_angles_to_R(*self.eo[3:6])


def measure(img, refs, io, min_inliers=40):
    """-> dict with the best reference's result."""
    f = features(img)
    best = None
    for ref in refs:
        p_ref, p_new = match(ref.feat, f)
        if len(p_ref) < 8:
            continue
        world = rays(p_ref[:, 0], p_ref[:, 1], io) @ ref.R      # R^T applied row-wise
        cam = rays(p_new[:, 0], p_new[:, 1], io)
        R, inl = solve_rotation(world, cam)
        n_in = int(inl.sum())
        if best is None or n_in > best["inliers"]:
            az, tilt, roll = R_to_angles(R)
            d = np.degrees(np.array([az, tilt, roll]) - ref.eo[3:6])
            d[0] = (d[0] + 180) % 360 - 180
            best = {"reference": ref.path.name, "matches": len(p_ref), "inliers": n_in,
                    "residual_deg": float(np.median(ang_err(R, world[inl], cam[inl]))) if n_in else np.nan,
                    "d_azimuth": d[0], "d_tilt": d[1], "d_roll": d[2],
                    "eo": np.r_[ref.eo[:3], az, tilt, roll], "eo_ref": ref.eo_name}
    if best is None or best["inliers"] < min_inliers:
        return {"status": "unmatched", **(best or {"matches": 0, "inliers": 0})}
    return best


# ------------------------------------------------------------- commands

def self_test(image=None):
    """Rotate a frame by known angles (exact for a camera that only rotates) and recover them."""
    from view_reproject import reproject_image
    io = load_intrinsics(CAL / "CACO05_c2_20240801_IO.yaml")
    eo = load_extrinsics(CAL / "CACO05_c2_20251113_EO-CV.yaml")
    if not image:
        sys.exit("self-test needs a c2 frame: pointing_check.py --self-test --image <c2 image.jpg>")
    img = cv2.imread(str(image))
    import tempfile
    tmp = Path(tempfile.mkdtemp(prefix="pointing_selftest_"))
    ref_img = tmp / "ref.jpg"; cv2.imwrite(str(ref_img), img)
    ref_img.with_suffix(".json").write_text(json.dumps({"eo": eo.tolist()}))
    refs = [Reference(ref_img)]
    worst = 0.0
    for daz, dtilt, droll in ((0, 0, 0), (0.30, -0.20, 0.10), (-1.0, 0.5, -0.3), (2.5, 0.0, 0.0)):
        eo2 = eo.copy(); eo2[3:6] += np.radians([daz, dtilt, droll])
        # what camera eo2 sees, built from the frame taken at eo (exact: same centre)
        moved = reproject_image(img, io, eo, eo2, 0.0)
        r = measure(moved, refs, io)
        if r.get("status") == "unmatched":
            print(f"  true ({daz:+.2f}, {dtilt:+.2f}, {droll:+.2f}): UNMATCHED"); worst = 99; continue
        err = np.abs(np.array([r["d_azimuth"], r["d_tilt"], r["d_roll"]]) - [daz, dtilt, droll])
        worst = max(worst, err.max())
        print(f"  true az/tilt/roll ({daz:+.2f}, {dtilt:+.2f}, {droll:+.2f}) deg -> measured "
              f"({r['d_azimuth']:+.3f}, {r['d_tilt']:+.3f}, {r['d_roll']:+.3f}); "
              f"{r['inliers']}/{r['matches']} inliers, error {err.max():.3f} deg")
    print(("PASS" if worst < 0.03 else "FAIL") + f": worst error {worst:.3f} deg")
    return 0 if worst < 0.03 else 1


def epoch_of(name):
    m = re.match(r"^(\d{9,11})\.", name)
    return int(m.group(1)) if m else None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--camera", default="c2")
    ap.add_argument("--image", help="measure this one image and print the result")
    ap.add_argument("--add-reference", metavar="IMAGE",
                    help="add IMAGE to the camera's reference bank, taken at --eo")
    ap.add_argument("--eo", help="EO yaml of the reference image (with --add-reference)")
    ap.add_argument("--image-dir", nargs="+",
                    default=["/mnt/I2Rgus_Data/ImageProducts/products", "/mnt/I2Rgus_Data/ImageProducts"])
    ap.add_argument("--product", default="timex", help="image product to check (default timex)")
    ap.add_argument("--log", default=None, help="default archive/pointing_<cam>.csv")
    ap.add_argument("--warn", type=float, default=0.1, help="deg; 'moved' above this (default 0.1)")
    ap.add_argument("--fail", type=float, default=0.5, help="deg; 'MOVED' above this (default 0.5)")
    ap.add_argument("--every", type=int, default=3600,
                    help="seconds between checked frames (default one per hour)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test(args.image)

    bank = CAL / "pointing" / args.camera
    io = load_intrinsics(CAL / f"CACO05_{args.camera}_20240801_IO.yaml")

    if args.add_reference:
        if not args.eo:
            sys.exit("--add-reference needs --eo (the pointing the image was taken at)")
        bank.mkdir(parents=True, exist_ok=True)
        src = Path(args.add_reference)
        dst = bank / src.name
        dst.write_bytes(src.read_bytes())
        eo = load_extrinsics(Path(args.eo) if Path(args.eo).exists() else CAL / args.eo)
        dst.with_suffix(".json").write_text(json.dumps({"eo": eo.tolist(), "eo_file": Path(args.eo).name,
                                                        "added": datetime.now(timezone.utc).isoformat()}, indent=1))
        n = len(list(bank.glob("*.jpg")))
        print(f"reference added: {dst} (bank now {n} image(s)). Add a few more taken in different light.")
        return 0

    refs = [Reference(p) for p in sorted(bank.glob("*.jpg")) if not p.stem.endswith("_mask")]
    if not refs:
        sys.exit(f"no reference images in {bank} -- add one with --add-reference IMAGE --eo EO.yaml")

    if args.image:
        img = cv2.imread(args.image)
        r = measure(img, refs, io)
        if r.get("status") == "unmatched":
            print(f"unmatched: {r.get('inliers', 0)} inliers (fog, night, rain, or a new view?)")
            return 1
        print(f"reference {r['reference']}: {r['inliers']}/{r['matches']} inliers, residual "
              f"{r['residual_deg']:.3f} deg")
        print(f"change: azimuth {r['d_azimuth']:+.3f}, tilt {r['d_tilt']:+.3f}, roll {r['d_roll']:+.3f} deg")
        az, tilt, roll = np.degrees(r["eo"][3:6])
        print(f"pointing now: azimuth {az:.4f}, tilt {tilt:.4f}, roll {roll:.4f} deg")
        return 0

    import pandas as pd
    log = Path(args.log or HERE / "archive" / f"pointing_{args.camera}.csv")
    log.parent.mkdir(parents=True, exist_ok=True)
    done = pd.read_csv(log) if log.exists() else pd.DataFrame(columns=FIELDS)
    seen = set(done["filename"])
    last = done["epoch"].max() if len(done) else 0
    files = {}
    for d in args.image_dir:
        for p in Path(d).glob(f"*.{args.camera}.{args.product}.jpg"):
            files.setdefault(p.name, p)
    todo = []
    for name in sorted(files):
        e = epoch_of(name)
        if e is None or name in seen or e < last + args.every:
            continue
        todo.append((e, name)); last = e
    rows = []
    for e, name in todo:
        img = cv2.imread(str(files[name]))
        if img is None:
            continue
        r = measure(img, refs, io)
        if r.get("status") == "unmatched":
            status = "unmatched"
        else:
            big = max(abs(r["d_azimuth"]), abs(r["d_tilt"]), abs(r["d_roll"]))
            status = "MOVED" if big > args.fail else "moved" if big > args.warn else "ok"
        rows.append({"epoch": e, "time_utc": datetime.fromtimestamp(e, tz=timezone.utc).isoformat(),
                     "filename": name, "reference": r.get("reference", ""),
                     "matches": r.get("matches", 0), "inliers": r.get("inliers", 0),
                     "residual_deg": round(r.get("residual_deg", np.nan), 4),
                     **{k: round(r[k], 4) if k in r else np.nan for k in ("d_azimuth", "d_tilt", "d_roll")},
                     "status": status})
        print(f"{name}: {status}" + (f"  az {r['d_azimuth']:+.3f} tilt {r['d_tilt']:+.3f} "
                                     f"roll {r['d_roll']:+.3f} deg ({r['inliers']} inliers)"
                                     if status != "unmatched" else ""))
    if rows:
        out = pd.concat([done, pd.DataFrame(rows, columns=FIELDS)], ignore_index=True)
        out.to_csv(log, index=False)
    print(f"pointing {args.camera}: {len(rows)} new frame(s) checked -> {log}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

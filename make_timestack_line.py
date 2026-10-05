#!/usr/bin/env python3
"""
Surf-Zone Timestack Lines For The Argus Pixel List
====================================================
The c2 timestack lines in arguseyes/build/c2_timestack.pix were laid out
for RUNUP: they cross the swash and upper beach. On the station 216 of 310
stacks had no pixel seaward of the swash, so the wave period could only
be measured as a mean (Tm01), not a dominant period. This makes
cross-shore lines that run from the dry beach out through the surf zone,
in image pixels for the camera's current calibration, ready to append to
the .pix file.

Each line runs along the seaward normal of the beach (81 deg) from --from
to --to metres past the Jan 2025 0 m NAVD88 contour, at --along metres
north along the beach (bearing 351 deg) from the reference point, sampled
every --step metres on the water plane --z. Points the camera cannot see
are dropped; consecutive points falling on the same pixel are merged.

c2 coverage (z = 0): 150 m along -> 0 m contour to ~110 m seaward, 10-17 px/m;
                     220 m along -> to ~170 m seaward,              8-13 px/m.

OUTPUT
  --output      the new lines only, "u v" per row, one line after another
  --combined    the existing .pix (--existing) with the new lines APPENDED, so the
                runup lines keep their numbers (1-4) and the new ones follow (5, 6)
  --csv         per point: line, u, v, easting, northing, metres past the 0 m contour
  --preview     the lines drawn on a camera frame (--image), to check before use

THE CAMERA MUST BE TOLD. The .pix file is read by the Argus collection
software (arguseyes) to decide which pixels to record; replacing it is a
station configuration change -- keep the old file, and restart collection
as that software requires. Stacks recorded before and after have
different column counts: runup_from_timestack.py and
timestack_wave_period.py skip stacks that do not match the current .pix.

Usage (on the station):
    python3 make_timestack_line.py --camera c2 --along 150 220 \\
        --existing /home/argus_user/arguseyes/build/c2_timestack.pix \\
        --combined c2_timestack_surf.pix --preview c2_surf_lines.jpg \\
        --image $(ls /mnt/I2Rgus_Data/ImageProducts/*.c2.timex.jpg | tail -1)
"""

import sys
import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
SHORE_ORIGIN = np.array([420150.0, 4638400.0])     # 0 m NAVD88 contour, Jan 2025 lidar (sea_patch.py)
ALONG_DEG, SEAWARD_DEG = 351.0, 81.0


def line_pixels(io, eo, along, d_from, d_to, step, z):
    """Ground line -> (u, v, E, N, d) integer pixels the camera sees, consecutive duplicates merged."""
    from view_reproject import ground_to_pixel
    d = np.arange(d_from, d_to + step / 2, step)
    p0 = SHORE_ORIGIN + along * np.array([np.sin(np.radians(ALONG_DEG)), np.cos(np.radians(ALONG_DEG))])
    E = p0[0] + d * np.sin(np.radians(SEAWARD_DEG))
    N = p0[1] + d * np.cos(np.radians(SEAWARD_DEG))
    U, V, ok = ground_to_pixel(E, N, z, io, eo)
    u, v = np.round(U[ok]).astype(int), np.round(V[ok]).astype(int)
    keep = np.r_[True, (np.diff(u) != 0) | (np.diff(v) != 0)]
    return u[keep], v[keep], E[ok][keep], N[ok][keep], d[ok][keep]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--camera", default="c2")
    ap.add_argument("--along", nargs="+", type=float, default=[150.0, 220.0],
                    help="metres north along the beach of each line (default 150 220)")
    ap.add_argument("--from", dest="d_from", type=float, default=-20.0,
                    help="start, metres past the 0 m contour (negative = up the beach; default -20)")
    ap.add_argument("--to", dest="d_to", type=float, default=250.0, help="end, metres seaward (default 250)")
    ap.add_argument("--step", type=float, default=1.0, help="ground spacing, m (default 1)")
    ap.add_argument("--z", type=float, default=0.0, help="water plane, m NAVD88 (default 0)")
    ap.add_argument("--eo", default=None, help="default calibration/CACO05_<cam>_20251113_EO-CV.yaml")
    ap.add_argument("--output", default=None, help="new lines only (default <cam>_surf_lines.pix)")
    ap.add_argument("--existing", default=None, help="the camera's current .pix file")
    ap.add_argument("--combined", default=None, help="existing lines + new lines, for the camera")
    ap.add_argument("--csv", default=None, help="per-point table (default <cam>_surf_lines.csv)")
    ap.add_argument("--image", default=None, help="a frame from this camera, for --preview")
    ap.add_argument("--preview", default=None, help="JPEG with the lines drawn on --image")
    args = ap.parse_args()

    from georectify import load_intrinsics, load_extrinsics
    io = load_intrinsics(HERE / "calibration" / f"CACO05_{args.camera}_20240801_IO.yaml")
    eo = load_extrinsics(Path(args.eo) if args.eo else
                         HERE / "calibration" / f"CACO05_{args.camera}_20251113_EO-CV.yaml")
    nu, nv = int(io[0]), int(io[1])

    lines, rows = [], []
    for k, along in enumerate(args.along, 1):
        u, v, E, N, d = line_pixels(io, eo, along, args.d_from, args.d_to, args.step, args.z)
        if len(u) < 10:
            print(f"line at {along:.0f} m along: the camera sees only {len(u)} point(s) -- skipped")
            continue
        rng = np.hypot(E - eo[0], N - eo[1])
        print(f"line at {along:.0f} m along: {len(u)} pixels, {d.min():+.0f} to {d.max():+.0f} m past the "
              f"0 m contour, {rng.min():.0f}-{rng.max():.0f} m from the camera, "
              f"u {u.min()}-{u.max()}, v {v.min()}-{v.max()}")
        lines.append((along, u, v))
        rows += [(k, a, b, e, n, dd) for a, b, e, n, dd in zip(u, v, E, N, d)]
    if not lines:
        sys.exit("no line is visible to this camera -- change --along / --to")

    out = Path(args.output or f"{args.camera}_surf_lines.pix")
    out.write_text("".join(f"{a} {b}\n" for _, u, v in lines for a, b in zip(u, v)))
    print(f"new lines         : {out} ({sum(len(u) for _, u, _ in lines)} pixels)")
    csv = Path(args.csv or f"{args.camera}_surf_lines.csv")
    csv.write_text("line,u,v,easting,northing,m_past_0m_contour\n" +
                   "".join(f"{k},{a},{b},{e:.2f},{n:.2f},{dd:.1f}\n" for k, a, b, e, n, dd in rows))
    print(f"point table       : {csv}")

    if args.combined:
        if not args.existing or not Path(args.existing).exists():
            sys.exit("--combined needs --existing (the camera's current .pix file)")
        old = Path(args.existing).read_text()
        if old and not old.endswith("\n"):
            old += "\n"
        n_old = len([l for l in old.splitlines() if l.strip()])
        from runup_from_timestack import split_lines
        old_pix = np.loadtxt(args.existing)[:, :2]
        n_lines_old = len(split_lines(old_pix))
        new = out.read_text()
        cols = len(old.splitlines()[0].split())
        if cols > 2:                     # keep the existing file's column layout
            print(f"NOTE: {args.existing} has {cols} columns per row; the new rows have 2 (u v). "
                  "Check the Argus documentation before installing.")
        Path(args.combined).write_text(old + new)
        first = n_lines_old + 1
        print(f"combined          : {args.combined} = {n_old} existing pixels ({n_lines_old} lines) "
              f"+ the new lines as line(s) {first}..{first + len(lines) - 1}")
        # the line splitter must see the new lines as separate lines
        chk = split_lines(np.loadtxt(args.combined)[:, :2])
        print(f"                    line splitter finds {len(chk)} lines "
              f"({'as expected' if len(chk) == n_lines_old + len(lines) else 'CHECK: not as expected'})")

    if args.preview:
        import cv2
        img = cv2.imread(args.image) if args.image else None
        if img is None:
            img = np.full((nv, nu, 3), 40, np.uint8)
        for k, (along, u, v) in enumerate(lines):
            pts = np.c_[u, v].reshape(-1, 1, 2).astype(np.int32)
            cv2.polylines(img, [pts], False, (0, 0, 0), 9, cv2.LINE_AA)
            cv2.polylines(img, [pts], False, (0, 255, 255), 4, cv2.LINE_AA)
            cv2.putText(img, f"{along:.0f} m", (int(u[-1]) + 10, int(v[-1])), cv2.FONT_HERSHEY_SIMPLEX,
                        1.6, (0, 255, 255), 4, cv2.LINE_AA)
        cv2.imwrite(args.preview, cv2.resize(img, (img.shape[1] // 2, img.shape[0] // 2)),
                    [cv2.IMWRITE_JPEG_QUALITY, 88])
        print(f"preview           : {args.preview} (yellow: the new lines; the seaward end is labelled)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

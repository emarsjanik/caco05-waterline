#!/usr/bin/env python3
"""
Waterlines Against A Survey: Datum Offset Or Waves?
======================================================
Compares every georectified waterline point with a surveyed surface
(lidar DSM, RTK grid) and splits the difference into its likely causes.

WHY. C. Sherwood compared the Jan 18-23 2025 waterline DEM with the
23 Jan 2025 lidar (Oct 2026): right place, right shape, but ~1.2-1.4 m
LOW. Two explanations predict different patterns, so the data can
decide between them:

  * A DATUM error (water level not on the survey's NAVD88, or the survey
    on another datum) shifts EVERY frame by the same amount, whatever the
    waves.
  * WAVE SETUP AND SWASH: a timex waterline marks where the swash
    reaches on average, above the still-water line, but is given the
    still-water elevation. That error GROWS with wave height (and
    period), and on calm days it should be small.

So for each frame: offset = survey elevation at its waterline points -
the elevation the waterline was given (median over the frame). Fitting
offset = a + b * sqrt(Hs * L0) across frames, the intercept a is the
part waves do not explain (datum, and anything else common to all
frames) and b is the setup/swash coefficient (Stockdon et al. 2006: mean
setup ~0.35 * beach slope * sqrt(Hs L0)). If a ~ 0 and b fits, it is
waves; if b ~ 0 and a ~ the whole offset, it is the datum.

It also reports the offset against water level (a dependence there
means a horizontal or slope error, not a datum one) and writes, with
--dem, the survey-minus-DEM grid.

SURVEY INPUT. An ESRI ASCII grid (.asc; in Global Mapper: File > Export >
Elevation Grid > Arc ASCII Grid, cropped to the beach) or a float
GeoTIFF readable by PIL. Same horizontal coordinates as the waterlines
(UTM zone 19N, metres); the elevations should be NAVD88 -- if the survey
is on another datum, that difference lands in the intercept.

Usage:
    python3 compare_dem_survey.py /mnt/I2Rgus_Data/Chelsea_calibration/contour_points_ground.csv \\
        --survey lidar_2025-01-23.asc --start-date 2025-01-18 --end-date 2025-01-23 \\
        --dem /mnt/I2Rgus_Data/Chelsea_calibration/dem/jan18_23_fine_dem.asc
"""

import sys
import csv
import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from asc_to_geotiff import read_asc, write_geotiff          # noqa: E402
from dem_from_contours import write_ascii_grid               # noqa: E402


def read_tiff(path):
    """
    Minimal GeoTIFF / BigTIFF reader, numpy + zlib only: one band, 32-bit float, strips or
    tiles, no compression or Deflate, predictor 1 (none) or 3 (floating point).
    The station's Pillow cannot open the 2025 YSMP lidar (tiled, Deflate,
    predictor 3). Returns (array rows north to south, {tag: values}).
    """
    import struct
    import zlib
    raw = Path(path).read_bytes()
    bo = {b"II": "<", b"MM": ">"}.get(raw[:2])
    version = struct.unpack(bo + "H", raw[2:4])[0] if bo else 0
    if version == 42:                                  # classic TIFF
        off = struct.unpack(bo + "I", raw[4:8])[0]
        n = struct.unpack(bo + "H", raw[off:off + 2])[0]
        head, esize, cfmt, ofmt, inline = 2, 12, "I", "I", 4
    elif version == 43:                                # BigTIFF (the 2025 YSMP lidar)
        off = struct.unpack(bo + "Q", raw[8:16])[0]
        n = struct.unpack(bo + "Q", raw[off:off + 8])[0]
        head, esize, cfmt, ofmt, inline = 8, 20, "Q", "Q", 8
    else:
        raise ValueError("not a TIFF")
    size = {1: 1, 2: 1, 3: 2, 4: 4, 11: 4, 12: 8, 16: 8}
    fmt = {1: "B", 3: "H", 4: "I", 11: "f", 12: "d", 16: "Q"}
    tags = {}
    for i in range(n):
        e = raw[off + head + esize * i: off + head + esize * (i + 1)]
        tag, typ = struct.unpack(bo + "HH", e[:4])
        cnt = struct.unpack(bo + cfmt, e[4:4 + struct.calcsize(cfmt)])[0]
        vo = 4 + struct.calcsize(cfmt)
        nb = size.get(typ, 1) * cnt
        data = e[vo:vo + nb] if nb <= inline else \
            raw[struct.unpack(bo + ofmt, e[vo:vo + struct.calcsize(ofmt)])[0]:][:nb]
        if typ == 2:
            tags[tag] = data.rstrip(b"\0").decode("latin-1")
        elif typ in fmt:
            tags[tag] = struct.unpack(bo + fmt[typ] * cnt, data)
    w, h = tags[256][0], tags[257][0]
    if tags.get(258, (32,))[0] != 32 or tags.get(339, (3,))[0] != 3 or tags.get(277, (1,))[0] != 1:
        raise ValueError("only single-band 32-bit float supported")
    comp, pred = tags.get(259, (1,))[0], tags.get(317, (1,))[0]
    if comp not in (1, 8, 32946) or pred not in (1, 3):
        raise ValueError(f"compression {comp} / predictor {pred} not supported")
    if 322 in tags:
        tw, th = tags[322][0], tags[323][0]
        offs, cnts = tags[324], tags[325]
        across = -(-w // tw)
        blocks = [((k // across) * th, (k % across) * tw) for k in range(len(offs))]
    else:
        tw, th = w, tags.get(278, (h,))[0]
        offs, cnts = tags[273], tags[279]
        blocks = [(k * th, 0) for k in range(len(offs))]
    out = np.full((h, w), np.nan, dtype=np.float32)
    for (r0, c0), o, c in zip(blocks, offs, cnts):
        buf = raw[o:o + c]
        if comp != 1:
            buf = zlib.decompress(buf)
        rows = len(buf) // (tw * 4)
        b = np.frombuffer(buf[:rows * tw * 4], dtype=np.uint8).reshape(rows, tw * 4)
        if pred == 3:
            # floating-point predictor: bytes differenced along the row, then
            # stored as byte planes, most significant first
            b = np.cumsum(b, axis=1, dtype=np.uint8)
            vals = b.reshape(rows, 4, tw).transpose(0, 2, 1).copy().view(">f4")[..., 0]
        else:
            vals = b.copy().view(bo + "f4")
        rr, cc = min(rows, h - r0), min(tw, w - c0)
        out[r0:r0 + rr, c0:c0 + cc] = vals[:rr, :cc]
    return out, tags


def read_survey(path):
    """-> (grid rows north to south, x_left, y_top, cell). NaN = nodata."""
    p = Path(path)
    if p.suffix.lower() == ".asc":
        d, h = read_asc(p)
        return d, h["xllcorner"], h["yllcorner"] + d.shape[0] * h["cellsize"], h["cellsize"]
    try:
        d, tags = read_tiff(p)
        d = d.astype(np.float64)
    except (ValueError, KeyError) as err:
        from PIL import Image                          # anything the reader above cannot do
        try:
            im = Image.open(p)
        except Exception:
            sys.exit(f"{p}: cannot read ({err}); export it as an Arc ASCII grid (.asc)")
        tags = dict(im.tag_v2)
        d = np.array(im, dtype=np.float64)
    scale, tie = tags.get(33550), tags.get(33922)
    if not scale or not tie:
        sys.exit(f"{p}: no GeoTIFF georeferencing tags; export it as an Arc ASCII grid (.asc)")
    nd = tags.get(42113)
    if nd is not None:
        try:
            d[d == float(str(nd).strip("\0 "))] = np.nan
        except ValueError:
            pass
    with np.errstate(invalid="ignore"):      # the nodata cells are NaN already (numpy 1.17 warns)
        d[d < -1e30] = np.nan
    if abs(scale[0] - scale[1]) > 1e-9:
        sys.exit(f"{p}: non-square pixels; export as .asc")
    x0, y0 = tie[3] - tie[0] * scale[0], tie[4] + tie[1] * scale[1]
    # GTRasterType (GeoKey 1025) 2 = PixelIsPoint: the tie point is the CENTRE of
    # the first pixel, not its outer corner (the 2025 YSMP lidar is written so).
    keys = tags.get(34735) or ()
    for i in range(4, len(keys) - 3, 4):
        if keys[i] == 1025 and keys[i + 3] == 2:
            x0, y0 = x0 - scale[0] / 2, y0 + scale[1] / 2
    return d, x0, y0, scale[0]


def sample(grid, x0, y0, cell, E, N):
    """Bilinear sample at cell centres; NaN outside or next to nodata."""
    fc = (E - x0) / cell - 0.5
    fr = (y0 - N) / cell - 0.5
    c0, r0 = np.floor(fc).astype(int), np.floor(fr).astype(int)
    tc, tr = fc - c0, fr - r0
    out = np.full(E.shape, np.nan)
    ok = (c0 >= 0) & (r0 >= 0) & (c0 + 1 < grid.shape[1]) & (r0 + 1 < grid.shape[0])
    c, r, a, b = c0[ok], r0[ok], tc[ok], tr[ok]
    out[ok] = (grid[r, c] * (1 - a) * (1 - b) + grid[r, c + 1] * a * (1 - b)
               + grid[r + 1, c] * (1 - a) * b + grid[r + 1, c + 1] * a * b)
    return out


def waves_vs_datum(off, hs, tp, tide):
    """
    Fits per-frame offsets as a + b*sqrt(Hs*L0) and prints the split:
    a = the part waves do not explain (datum, common bias), b = setup/swash.
    Also the trend against water level. Returns (a, b, phi, has_waves).
    """
    phi = np.sqrt(hs * 9.81 * tp ** 2 / (2 * np.pi))         # sqrt(Hs L0), m
    w = np.isfinite(phi)
    a = b = np.nan
    spread_phi = float(np.ptp(phi[w])) if w.any() else 0.0
    if w.sum() >= 5 and spread_phi < 2:
        print(f"WAVES vs DATUM    : not fitted -- the frames' wave forcing is nearly the same "
              f"(sqrt(Hs*L0) spans {spread_phi:.1f} m; Hs {np.nanmin(hs):.2f}-{np.nanmax(hs):.2f} m), "
              f"so a datum shift and a wave effect cannot be told apart. Frames from days with "
              f"different seas are needed.")
    elif w.sum() >= 5:
        X = np.column_stack([np.ones(w.sum()), phi[w]])
        (a, b), *_ = np.linalg.lstsq(X, off[w], rcond=None)
        res = off[w] - X @ np.array([a, b])
        r2 = 1 - (res ** 2).sum() / max(((off[w] - off[w].mean()) ** 2).sum(), 1e-12)
        print(f"WAVES vs DATUM    : offset = {a:+.3f} + {b:.4f} * sqrt(Hs*L0)   "
              f"(R2 {r2:.2f}, {w.sum()} frames with Hs and Tp)")
        print(f"   Hs {np.nanmin(hs):.2f}-{np.nanmax(hs):.2f} m, Tp {np.nanmin(tp):.1f}-"
              f"{np.nanmax(tp):.1f} s -> wave part {b * phi[w].min():+.2f} to {b * phi[w].max():+.2f} m")
        print(f"   intercept {a:+.3f} m = the part waves do not explain (datum, or a common bias)")
        print(f"   slope {b:.4f}: Stockdon mean setup alone would be ~0.35 x beach slope, "
              f"i.e. {0.35 * 0.1:.3f} for a 0.1 slope; swash adds to it")
        for lo, hi in ((0, 0.75), (0.75, 1.25), (1.25, 2), (2, 9)):
            m = (hs >= lo) & (hs < hi)
            if m.any():
                print(f"   Hs {lo:.2f}-{hi:.2f} m: {m.sum():3d} frames, offset median "
                      f"{np.median(off[m]):+.3f} m")
    else:
        print("WAVES vs DATUM    : fewer than 5 frames with Hs and Tp -- cannot separate "
              "(run extract_elevation_contours.py with --waves)")
    if len(off) >= 3:
        X = np.column_stack([np.ones(len(off)), tide])
        (a2, b2), *_ = np.linalg.lstsq(X, off, rcond=None)
        print(f"vs WATER LEVEL    : offset = {a2:+.3f} + {b2:+.3f} * level   "
              f"(a trend here means a slope or cross-shore position error, not a datum shift)")
    return a, b, phi, w


def horizontal_shift(grid, x0, y0, cell, E, N, Z, max_shift=60.0, step=0.25):
    """
    Cross-shore distance (m) a frame's line must move so the survey elevation
    under it equals the elevation the line was given; + = seaward. The beach
    normal is perpendicular to the line's own principal axis, pointing downhill
    in the survey. NaN if no shift within +/-max_shift fits better than 0.15 m.
    """
    if len(E) < 5:
        return np.nan
    c = np.array([E.mean(), N.mean()])
    _, vecs = np.linalg.eigh(np.cov(np.stack([E - c[0], N - c[1]])))
    nrm = vecs[:, 0]
    up = sample(grid, x0, y0, cell, np.array([c[0] + 5 * nrm[0]]), np.array([c[1] + 5 * nrm[1]]))[0]
    dn = sample(grid, x0, y0, cell, np.array([c[0] - 5 * nrm[0]]), np.array([c[1] - 5 * nrm[1]]))[0]
    if np.isfinite(up) and np.isfinite(dn) and up > dn:
        nrm = -nrm                                      # point seaward (downhill)
    best, best_err = np.nan, 0.15
    for dlt in np.arange(-max_shift, max_shift + step, step):
        z = sample(grid, x0, y0, cell, E + dlt * nrm[0], N + dlt * nrm[1])
        ok = np.isfinite(z)
        if ok.sum() < max(5, 0.5 * len(E)):
            continue
        err = float(np.median(np.abs(z[ok] - Z[ok])))
        if err < best_err:
            best, best_err = float(dlt), err
    return best


def draw_map(path, grid, x0, y0, cell, mapped):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cams = sorted({m[0] for m in mapped})
    allE = np.concatenate([m[1] for m in mapped]); allN = np.concatenate([m[2] for m in mapped])
    pad = 40
    xa, xb, ya, yb = allE.min() - pad, allE.max() + pad, allN.min() - pad, allN.max() + pad
    c0, c1 = max(int((xa - x0) / cell), 0), min(int((xb - x0) / cell), grid.shape[1])
    r0, r1 = max(int((y0 - yb) / cell), 0), min(int((y0 - ya) / cell), grid.shape[0])
    sub = grid[r0:r1, c0:c1]
    ext = [x0 + c0 * cell, x0 + c1 * cell, y0 - r1 * cell, y0 - r0 * cell]
    xs = x0 + (np.arange(c0, c1) + 0.5) * cell
    ys = y0 - (np.arange(r0, r1) + 0.5) * cell
    fig, axes = plt.subplots(1, len(cams), figsize=(6.5 * len(cams), 9), dpi=100, squeeze=False)
    for ax, cam in zip(axes[0], cams):
        ax.imshow(sub, extent=ext, cmap="gray", vmin=-1.5, vmax=6, origin="upper")
        cs = ax.contour(xs, ys, sub, levels=np.arange(-1, 4.01, 0.5), colors="#6b6a64",
                        linewidths=0.6)
        ax.clabel(cs, fontsize=7, fmt="%.1f")
        pts = [m for m in mapped if m[0] == cam]
        E = np.concatenate([m[1] for m in pts]); N = np.concatenate([m[2] for m in pts])
        D = np.concatenate([m[3] for m in pts])
        k = slice(None, None, max(1, len(E) // 20000))
        sc = ax.scatter(E[k], N[k], c=D[k], cmap="RdBu_r", vmin=-2, vmax=2, s=3, lw=0)
        ax.set_title(f"{cam}: survey - waterline elevation (m)\nred = waterline too LOW "
                     f"(line sits too far landward), blue = too high", loc="left", fontsize=9)
        ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal")
        ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
        fig.colorbar(sc, ax=ax, shrink=0.6)
    fig.tight_layout()
    fig.savefig(path)
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("contours", help="contour_points_ground.csv (georectified waterline points)")
    ap.add_argument("--survey", required=True, help="Survey surface: .asc or float GeoTIFF")
    ap.add_argument("--start-date"); ap.add_argument("--end-date")
    ap.add_argument("--camera", default="both")
    ap.add_argument("--dem", help="Our DEM (.asc) to difference against the survey")
    ap.add_argument("--min-points", type=int, default=20,
                    help="Points a frame needs on the survey to count (default 20)")
    ap.add_argument("--output", default=None,
                    help="Per-frame table CSV (default <survey stem>_frames.csv)")
    ap.add_argument("--plot", default=None, help="PNG of offset vs wave forcing")
    ap.add_argument("--map", default=None,
                    help="PNG map: the survey with every compared waterline point coloured by "
                         "its offset, one panel per camera")
    args = ap.parse_args()

    grid, x0, y0, cell = read_survey(args.survey)
    print(f"survey            : {args.survey}  {grid.shape[1]} x {grid.shape[0]} at {cell:g} m, "
          f"z {np.nanmin(grid):+.2f} to {np.nanmax(grid):+.2f}")

    frames = {}
    with open(args.contours, newline="") as f:
        for r in csv.DictReader(f):
            if not r.get("easting_utm19"):
                continue
            day = r.get("capture_time_utc", "")[:10]
            if (args.start_date and day < args.start_date) or (args.end_date and day > args.end_date):
                continue
            if args.camera != "both" and r.get("camera") != args.camera:
                continue
            fr = frames.setdefault(r["source_file"], {
                "time": r.get("capture_time_utc", ""), "cam": r.get("camera", ""),
                "E": [], "N": [], "Z": [],
                "tide": float(r["tide_elevation_navd88"]),
                "hs": float(r["offshore_hs_m"]) if r.get("offshore_hs_m") else np.nan,
                "tp": float(r["offshore_tp_s"]) if r.get("offshore_tp_s") else np.nan})
            fr["E"].append(float(r["easting_utm19"]))
            fr["N"].append(float(r["northing_utm19"]))
            fr["Z"].append(float(r.get("beach_elevation_navd88") or r["tide_elevation_navd88"]))
    if not frames:
        sys.exit("no georectified waterline points in that date range")

    rows, all_d, mapped = [], [], []
    for name, fr in sorted(frames.items(), key=lambda kv: kv[1]["time"]):
        E, N, Z = np.array(fr["E"]), np.array(fr["N"]), np.array(fr["Z"])
        s = sample(grid, x0, y0, cell, E, N)
        ok = np.isfinite(s)
        d = (s - Z)[ok]
        if len(d) < args.min_points:
            continue
        shift = horizontal_shift(grid, x0, y0, cell, E[ok], N[ok], Z[ok])
        all_d.append(d)
        mapped.append((fr["cam"], E[ok], N[ok], d))
        rows.append({"frame": name, "time_utc": fr["time"], "camera": fr["cam"],
                     "tide_navd88": fr["tide"],
                     "hs_m": fr["hs"], "tp_s": fr["tp"], "n_points": len(d),
                     "offset_m": float(np.median(d)),
                     "spread_m": float(np.percentile(d, 84) - np.percentile(d, 16)),
                     "shift_m": shift})
    if not rows:
        sys.exit("no frame has enough points on the survey -- do the extents overlap? "
                 "(check the survey's coordinate system: UTM 19N metres)")
    all_d = np.concatenate(all_d)
    off = np.array([r["offset_m"] for r in rows])
    print(f"frames compared   : {len(rows)} of {len(frames)} ({len(all_d)} points)")
    print(f"OFFSET survey - waterline elevation: median {np.median(all_d):+.3f} m, "
          f"frame medians {np.percentile(off, 10):+.2f} to {np.percentile(off, 90):+.2f} m (p10-p90)")
    print("   positive = our waterline elevations are LOW")
    print("   per camera (shift = how far the line would have to move seaward(+) across the "
          "beach to sit on its own elevation in the survey):")
    for cam in sorted({r["camera"] for r in rows}):
        rc = [r for r in rows if r["camera"] == cam]
        o = np.array([r["offset_m"] for r in rc])
        sh = np.array([r["shift_m"] for r in rc], float)
        sh = sh[np.isfinite(sh)]
        print(f"   {cam or '?':3s}: {len(rc):3d} frames, offset median {np.median(o):+.3f} m "
              f"(p10-p90 {np.percentile(o, 10):+.2f} to {np.percentile(o, 90):+.2f}), shift median "
              + (f"{np.median(sh):+.1f} m (p10-p90 {np.percentile(sh, 10):+.1f} to "
                 f"{np.percentile(sh, 90):+.1f})" if len(sh) else "--"))
    print("   A DATUM error gives a roughly constant vertical offset (the shift then varies with "
          "the beach slope);\n   a GEOMETRY (camera pointing) error gives a roughly constant "
          "shift, and huge vertical offsets where the beach is flat (berm).")

    tide = np.array([r["tide_navd88"] for r in rows])
    hs = np.array([r["hs_m"] for r in rows])
    tp = np.array([r["tp_s"] for r in rows])
    print()
    a, b, phi, w = waves_vs_datum(off, hs, tp, tide)

    if args.map:
        draw_map(args.map, grid, x0, y0, cell, mapped)

    out = args.output or str(Path(args.survey).with_suffix("")) + "_frames.csv"
    with open(out, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    print(f"\nwrote {out}")

    if args.dem:
        dem, h = read_asc(args.dem)
        c = h["cellsize"]
        rr, cc = np.nonzero(np.isfinite(dem))
        top = h["yllcorner"] + dem.shape[0] * c
        E = h["xllcorner"] + (cc + 0.5) * c
        N = top - (rr + 0.5) * c
        s = sample(grid, x0, y0, cell, E, N)
        diff = np.full(dem.shape, np.nan)
        diff[rr, cc] = s - dem[rr, cc]
        ok = np.isfinite(diff)
        print(f"DEM vs survey     : {int(ok.sum())} cells, survey - DEM median "
              f"{np.nanmedian(diff):+.3f} m, p16-p84 {np.nanpercentile(diff, 16):+.2f} to "
              f"{np.nanpercentile(diff, 84):+.2f} m")
        # accuracy by distance from the camera tower: errors that cancel in one
        # median are visible here (far cells have coarse pixels and grazing views)
        from georectify import load_extrinsics
        cam = load_extrinsics(Path(__file__).resolve().parent / "calibration" / "CACO05_c2_20251113_EO-CV.yaml")
        rng_m = np.hypot(E - cam[0], N - cam[1])
        dv = diff[rr, cc]
        print(f"  {'distance':>12} {'cells':>7} {'median':>8} {'NMAD':>6} {'RMSE':>6} {'p95 |e|':>8}")
        bins = [0, 100, 150, 200, 250, 300, 400, 600, 1e9]
        acc = []
        for lo, hi in zip(bins, bins[1:]):
            m = np.isfinite(dv) & (rng_m >= lo) & (rng_m < hi)
            if m.sum() < 20:
                continue
            e = dv[m]
            med = float(np.median(e)); nmad = float(1.4826 * np.median(np.abs(e - med)))
            rmse = float(np.sqrt(np.mean(e ** 2))); p95 = float(np.percentile(np.abs(e), 95))
            label = f"{lo:.0f}-{hi:.0f} m" if hi < 1e8 else f">{lo:.0f} m"
            acc.append({"distance": label, "cells": int(m.sum()), "median": med, "nmad": nmad,
                        "rmse": rmse, "p95_abs": p95})
            print(f"  {label:>12} {int(m.sum()):>7} {med:>+8.3f} {nmad:>6.3f} {rmse:>6.3f} {p95:>8.3f}")
        e = dv[np.isfinite(dv)]
        print(f"  {'all':>12} {len(e):>7} {np.median(e):>+8.3f} "
              f"{1.4826 * np.median(np.abs(e - np.median(e))):>6.3f} {np.sqrt(np.mean(e ** 2)):>6.3f} "
              f"{np.percentile(np.abs(e), 95):>8.3f}   (survey - DEM, m; NMAD = robust sd)")
        if acc:
            with open(str(Path(args.dem).with_suffix("")) + "_accuracy_by_distance.csv", "w", newline="") as f:
                wr = csv.DictWriter(f, fieldnames=list(acc[0]))
                wr.writeheader(); wr.writerows(acc)
        stem = str(Path(args.dem).with_suffix("")) + "_minus_survey"
        write_ascii_grid(stem + ".asc", np.flipud(-diff), h["xllcorner"], h["yllcorner"], c)
        print(f"wrote {stem}.asc (DEM - survey, m; negative = DEM low)")

    if args.plot and w.sum() >= 5:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), dpi=110)
        sc = axes[0].scatter(phi[w], off[w], c=tide[w], cmap="viridis", s=18)
        xs = np.linspace(0, np.nanmax(phi) * 1.05, 2)
        axes[0].plot(xs, a + b * xs, color="#1f1f1e", lw=1.5)
        axes[0].axhline(0, color="#6b6a64", lw=0.8)
        axes[0].set_xlabel("sqrt(Hs x L0)  (m)")
        axes[0].set_ylabel("survey - waterline elevation (m)")
        axes[0].set_title(f"intercept {a:+.2f} m (not waves), slope {b:.3f}", loc="left", fontsize=10)
        fig.colorbar(sc, ax=axes[0], label="water level (m NAVD88)")
        axes[1].scatter(tide, off, s=18, color="#2a78d6")
        axes[1].axhline(0, color="#6b6a64", lw=0.8)
        axes[1].set_xlabel("water level (m NAVD88)")
        axes[1].set_ylabel("survey - waterline elevation (m)")
        for ax in axes:
            ax.grid(True, color="#e4e3dc"); ax.set_axisbelow(True)
        fig.tight_layout()
        fig.savefig(args.plot)
        print(f"wrote {args.plot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

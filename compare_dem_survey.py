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


def read_survey(path):
    """-> (grid rows north to south, x_left, y_top, cell). NaN = nodata."""
    p = Path(path)
    if p.suffix.lower() == ".asc":
        d, h = read_asc(p)
        return d, h["xllcorner"], h["yllcorner"] + d.shape[0] * h["cellsize"], h["cellsize"]
    from PIL import Image
    im = Image.open(p)
    tags = im.tag_v2
    scale, tie = tags.get(33550), tags.get(33922)
    if not scale or not tie:
        sys.exit(f"{p}: no GeoTIFF georeferencing tags; export it as an Arc ASCII grid (.asc)")
    d = np.array(im, dtype=np.float64)
    nd = tags.get(42113)
    if nd is not None:
        d[d == float(str(nd).strip("\0"))] = np.nan
    d[d < -1e30] = np.nan
    if abs(scale[0] - scale[1]) > 1e-9:
        sys.exit(f"{p}: non-square pixels; export as .asc")
    return d, tie[3] - tie[0] * scale[0], tie[4] + tie[1] * scale[1], scale[0]


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
                "time": r.get("capture_time_utc", ""), "E": [], "N": [], "Z": [],
                "tide": float(r["tide_elevation_navd88"]),
                "hs": float(r["offshore_hs_m"]) if r.get("offshore_hs_m") else np.nan,
                "tp": float(r["offshore_tp_s"]) if r.get("offshore_tp_s") else np.nan})
            fr["E"].append(float(r["easting_utm19"]))
            fr["N"].append(float(r["northing_utm19"]))
            fr["Z"].append(float(r.get("beach_elevation_navd88") or r["tide_elevation_navd88"]))
    if not frames:
        sys.exit("no georectified waterline points in that date range")

    rows, all_d = [], []
    for name, fr in sorted(frames.items(), key=lambda kv: kv[1]["time"]):
        s = sample(grid, x0, y0, cell, np.array(fr["E"]), np.array(fr["N"]))
        d = s - np.array(fr["Z"])
        d = d[np.isfinite(d)]
        if len(d) < args.min_points:
            continue
        all_d.append(d)
        rows.append({"frame": name, "time_utc": fr["time"], "tide_navd88": fr["tide"],
                     "hs_m": fr["hs"], "tp_s": fr["tp"], "n_points": len(d),
                     "offset_m": float(np.median(d)),
                     "spread_m": float(np.percentile(d, 84) - np.percentile(d, 16))})
    if not rows:
        sys.exit("no frame has enough points on the survey -- do the extents overlap? "
                 "(check the survey's coordinate system: UTM 19N metres)")
    all_d = np.concatenate(all_d)
    off = np.array([r["offset_m"] for r in rows])
    print(f"frames compared   : {len(rows)} of {len(frames)} ({len(all_d)} points)")
    print(f"OFFSET survey - waterline elevation: median {np.median(all_d):+.3f} m, "
          f"frame medians {np.percentile(off, 10):+.2f} to {np.percentile(off, 90):+.2f} m (p10-p90)")
    print("   positive = our waterline elevations are LOW")

    tide = np.array([r["tide_navd88"] for r in rows])
    hs = np.array([r["hs_m"] for r in rows])
    tp = np.array([r["tp_s"] for r in rows])
    phi = np.sqrt(hs * 9.81 * tp ** 2 / (2 * np.pi))         # sqrt(Hs L0), m
    w = np.isfinite(phi)
    print()
    if w.sum() >= 5:
        X = np.column_stack([np.ones(w.sum()), phi[w]])
        (a, b), *_ = np.linalg.lstsq(X, off[w], rcond=None)
        res = off[w] - X @ np.array([a, b])
        r2 = 1 - (res ** 2).sum() / ((off[w] - off[w].mean()) ** 2).sum()
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
        print("WAVES vs DATUM    : no Hs/Tp in the contour file (run extract_elevation_contours.py "
              "with --waves); cannot separate")
    X = np.column_stack([np.ones(len(off)), tide])
    (a2, b2), *_ = np.linalg.lstsq(X, off, rcond=None)
    print(f"vs WATER LEVEL    : offset = {a2:+.3f} + {b2:+.3f} * level   "
          f"(a trend here means a slope or cross-shore position error, not a datum shift)")

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

#!/usr/bin/env python3
"""
Intertidal DEM From Waterline Contours
-----------------------------------------
Grids georectified waterline points into a beach elevation model.

THE IDEA: each detected waterline point is a place where the water
surface met the sand, so the BEACH elevation there equals the WATER
elevation, which GNSS-R measured independently. Over a tidal cycle the
waterline sweeps across the intertidal zone and traces out the beach
profile. This grids those traces into a surface.

WHAT THIS DEM IS AND IS NOT:
  * It covers the INTERTIDAL ZONE ONLY -- between the lowest and
    highest water levels observed. There is no information above the
    high-water line or below the low-water line, and none is invented.
  * Cells are filled ONLY where waterlines actually passed. Gaps are
    left as nodata rather than interpolated, because a smooth surface
    across an unmeasured gap looks like data and is not. Use
    --fill-gaps if you need a continuous surface and accept that.
  * Every cell carries a point count and an elevation spread, written
    as companion grids. A cell built from 2 points spanning 0.4 m is
    not the same measurement as one from 200 points spanning 0.03 m,
    and the DEM alone cannot show that difference.

SOURCES OF ERROR, roughly in order of size:
  * Waterline detection: the largest term, and it varies across the
    frame. The far field is worse than the near field.
  * Georectification: 0.33 m horizontal RMS against surveyed control,
    as of the Nov 2025 calibration. Unverified since.
  * Water level: GNSS-R, a measurement rather than a model.
  * Wave runup: the detected edge is the runup limit, which sits ABOVE
    the still-water line by an amount that grows with wave height. This
    is a systematic positive bias in elevation and is NOT corrected
    here -- it needs wave data the station does not currently provide.

Usage:
    python3 dem_from_contours.py contour_points_ground.csv dem_out \\
        [--cell 2.0] [--min-points 3] [--camera c1|c2|both]
        [--max-spread 0.5] [--fill-gaps] [--start-date ...] [--end-date ...]
"""

import sys
import csv
import argparse
from pathlib import Path

import numpy as np


def load_points(path, camera=None, start_date=None, end_date=None):
    """Reads georectified contour points. Rows without ground coordinates are skipped."""
    E, N, Z, cams, dates = [], [], [], [], []
    missing_ground = 0
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        if "easting_utm19" not in (reader.fieldnames or []):
            print("ERROR: no 'easting_utm19' column. Run georectify.py on the contour file first.")
            sys.exit(1)
        for r in reader:
            if not r.get("easting_utm19") or not r.get("northing_utm19"):
                missing_ground += 1
                continue
            if camera and camera != "both" and r["camera"] != camera:
                continue
            day = r.get("capture_time_utc", "")[:10]
            if start_date and day < start_date:
                continue
            if end_date and day > end_date:
                continue
            E.append(float(r["easting_utm19"]))
            N.append(float(r["northing_utm19"]))
            Z.append(float(r["tide_elevation_navd88"]))
            cams.append(r["camera"])
            dates.append(day)
    return (np.array(E), np.array(N), np.array(Z),
            np.array(cams), np.array(dates), missing_ground)


def build_grid(E, N, Z, cell, min_points, max_spread):
    """
    Bins points into cells and takes the MEDIAN elevation of each.

    Median rather than mean: a cell can catch an outlier from a single
    bad detection, and one wild value would drag a mean while barely
    moving a median. With repeat tidal crossings most cells hold many
    samples, so the median is well determined.
    """
    e0 = np.floor(E.min() / cell) * cell
    n0 = np.floor(N.min() / cell) * cell
    ncols = int(np.ceil((E.max() - e0) / cell)) + 1
    nrows = int(np.ceil((N.max() - n0) / cell)) + 1

    col = ((E - e0) / cell).astype(int)
    row = ((N - n0) / cell).astype(int)
    flat = row * ncols + col

    order = np.argsort(flat)
    flat_s, Z_s = flat[order], Z[order]
    edges = np.flatnonzero(np.diff(flat_s)) + 1
    starts = np.concatenate([[0], edges])
    ends = np.concatenate([edges, [len(flat_s)]])

    dem = np.full(nrows * ncols, np.nan)
    count = np.zeros(nrows * ncols, dtype=int)
    spread = np.full(nrows * ncols, np.nan)

    for s, e in zip(starts, ends):
        idx = flat_s[s]
        vals = Z_s[s:e]
        count[idx] = len(vals)
        if len(vals) >= min_points:
            dem[idx] = np.median(vals)
            # Spread as the 16th-84th percentile range -- a robust
            # stand-in for +/-1 sigma that a couple of outliers cannot
            # inflate the way a standard deviation can.
            spread[idx] = (np.percentile(vals, 84) - np.percentile(vals, 16)) if len(vals) > 2 \
                else float(vals.max() - vals.min())

    if max_spread and max_spread > 0:
        too_noisy = np.isfinite(spread) & (spread > max_spread)
        dem[too_noisy] = np.nan

    return (dem.reshape(nrows, ncols), count.reshape(nrows, ncols),
            spread.reshape(nrows, ncols), e0, n0, ncols, nrows)


def fill_small_gaps(dem, max_iterations=3):
    """
    Fills isolated nodata cells from their immediate neighbours.

    Deliberately limited: it closes pinholes inside measured areas, not
    voids between them. Each pass fills only cells with at least five
    of eight neighbours present, so a gap wider than a few cells stays
    open. Filling a real gap would manufacture terrain.
    """
    out = dem.copy()
    for _ in range(max_iterations):
        holes = np.isnan(out)
        if not holes.any():
            break
        padded = np.pad(out, 1, constant_values=np.nan)
        stack = np.stack([padded[a:a + out.shape[0], b:b + out.shape[1]]
                          for a in range(3) for b in range(3)
                          if not (a == 1 and b == 1)])
        neighbours = np.sum(np.isfinite(stack), axis=0)
        with np.errstate(invalid="ignore"):
            mean_nb = np.nanmean(stack, axis=0)
        fill = holes & (neighbours >= 5)
        if not fill.any():
            break
        out[fill] = mean_nb[fill]
    return out


def write_ascii_grid(path, grid, e0, n0, cell, nodata=-9999.0):
    """ESRI ASCII grid -- readable by QGIS, ArcGIS, GDAL, MATLAB."""
    flipped = np.flipud(grid)  # ASCII grids run north to south
    with open(path, "w") as f:
        f.write(f"ncols {grid.shape[1]}\n")
        f.write(f"nrows {grid.shape[0]}\n")
        f.write(f"xllcorner {e0:.3f}\n")
        f.write(f"yllcorner {n0:.3f}\n")
        f.write(f"cellsize {cell}\n")
        f.write(f"NODATA_value {nodata}\n")
        for r in flipped:
            f.write(" ".join("%.4f" % (v if np.isfinite(v) else nodata) for v in r) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("contour_csv")
    ap.add_argument("output_stem")
    ap.add_argument("--cell", type=float, default=2.0,
                    help="Grid cell size in metres (default 2.0). Below the georectification's "
                         "own 0.33 m accuracy there is nothing to gain; well above it the "
                         "beach profile gets smoothed away.")
    ap.add_argument("--min-points", type=int, default=3,
                    help="Minimum samples for a cell to be filled (default 3). One point is a "
                         "single detection with no way to tell whether it was a good one.")
    ap.add_argument("--max-spread", type=float, default=0.5,
                    help="Blank cells whose 16-84 percentile elevation range exceeds this, in "
                         "metres (default 0.5). A cell where repeat crossings disagree by more "
                         "than half a metre is not measuring one surface.")
    ap.add_argument("--camera", default="both", choices=["c1", "c2", "both"])
    ap.add_argument("--start-date"); ap.add_argument("--end-date")
    ap.add_argument("--fill-gaps", action="store_true",
                    help="Close isolated single-cell holes from neighbours. Off by default: an "
                         "interpolated cell looks identical to a measured one in the output.")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    E, N, Z, cams, dates, missing = load_points(
        args.contour_csv, args.camera, args.start_date, args.end_date)

    if len(E) == 0:
        print("No georectified points matched. Check --camera and the date range.")
        sys.exit(1)

    print("=" * 74)
    print("INTERTIDAL DEM")
    print("=" * 74)
    print(f"points            : {len(E)}")
    if missing:
        print(f"  (skipped {missing} row(s) with no ground coordinates)")
    for c in sorted(set(cams)):
        print(f"    {c}: {int((cams == c).sum())}")
    print(f"dates             : {min(dates)} to {max(dates)}  "
          f"({len(set(dates))} day(s))")
    print(f"elevation range   : {Z.min():+.2f} to {Z.max():+.2f} m NAVD88")
    print(f"extent            : {E.max()-E.min():.1f} m E-W by {N.max()-N.min():.1f} m N-S")
    print()

    dem, count, spread, e0, n0, ncols, nrows = build_grid(
        E, N, Z, args.cell, args.min_points, args.max_spread)

    total_cells = dem.size
    with_any = int((count > 0).sum())
    filled = int(np.isfinite(dem).sum())
    blanked = with_any - filled

    print(f"grid              : {nrows} x {ncols} cells at {args.cell} m")
    print(f"cells with points : {with_any} ({100*with_any/total_cells:.1f}% of grid)")
    print(f"cells filled      : {filled}")
    if blanked:
        print(f"cells blanked     : {blanked}  (fewer than {args.min_points} points, "
              f"or spread > {args.max_spread} m)")

    occupied = count[count > 0]
    if len(occupied):
        print(f"points per cell   : median {int(np.median(occupied))}, "
              f"p10 {int(np.percentile(occupied,10))}, p90 {int(np.percentile(occupied,90))}")
    valid_spread = spread[np.isfinite(spread) & np.isfinite(dem)]
    if len(valid_spread):
        print(f"elevation spread  : median {np.median(valid_spread):.3f} m, "
              f"p90 {np.percentile(valid_spread,90):.3f} m")
        print("                    (repeat crossings of the same cell; this is the DEM's")
        print("                     own repeatability, not its accuracy)")

    if args.fill_gaps:
        before = filled
        dem = fill_small_gaps(dem)
        after = int(np.isfinite(dem).sum())
        print(f"gap fill          : +{after-before} cell(s) interpolated from neighbours")

    stem = Path(args.output_stem)
    write_ascii_grid(str(stem) + "_dem.asc", dem, e0, n0, args.cell)
    write_ascii_grid(str(stem) + "_spread.asc", spread, e0, n0, args.cell)
    write_ascii_grid(str(stem) + "_count.asc",
                     count.astype(float), e0, n0, args.cell, nodata=0.0)
    print()
    print(f"wrote {stem}_dem.asc     (elevation, m NAVD88)")
    print(f"wrote {stem}_spread.asc  (16-84 percentile range, m)")
    print(f"wrote {stem}_count.asc   (samples per cell)")

    if not args.no_plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, axes = plt.subplots(1, 2, figsize=(16, 7), dpi=110)
            extent = [e0, e0 + ncols * args.cell, n0, n0 + nrows * args.cell]

            im0 = axes[0].imshow(dem, origin="lower", extent=extent,
                                 cmap="terrain", aspect="equal")
            axes[0].set_title(f"Intertidal DEM  ({min(dates)} to {max(dates)})\n"
                              f"{filled} cells at {args.cell} m", fontsize=10)
            plt.colorbar(im0, ax=axes[0], label="elevation (m NAVD88)", shrink=0.8)

            im1 = axes[1].imshow(spread, origin="lower", extent=extent,
                                 cmap="magma", aspect="equal", vmin=0,
                                 vmax=args.max_spread if args.max_spread else None)
            axes[1].set_title("Repeatability\n16-84 percentile elevation range", fontsize=10)
            plt.colorbar(im1, ax=axes[1], label="spread (m)", shrink=0.8)

            for ax in axes:
                ax.set_xlabel("easting (m, UTM 19N)")
                ax.set_ylabel("northing (m, UTM 19N)")
                ax.ticklabel_format(useOffset=False, style="plain")
            plt.tight_layout()
            plt.savefig(str(stem) + "_dem.png", bbox_inches="tight")
            print(f"wrote {stem}_dem.png")
        except Exception as exc:
            print(f"(plot skipped: {exc})")

    print()
    print("REMINDER: intertidal zone only, between the lowest and highest water")
    print("levels observed. No data above or below, and none invented. The detected")
    print("edge is the wave RUNUP limit, which sits above still water by an amount")
    print("that grows with wave height -- a systematic positive bias not corrected here.")


if __name__ == "__main__":
    main()

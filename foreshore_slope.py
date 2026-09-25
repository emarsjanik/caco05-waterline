#!/usr/bin/env python3
"""
Foreshore Slope From The Intertidal DEM
==========================================
Extracts the local foreshore beach slope needed by the Stockdon et al.
(2006) runup parameterization, and writes it as a grid alongside the
DEM.

WHY A GRID RATHER THAN ONE NUMBER: Stockdon et al. showed that using a
single alongshore-averaged slope where the foreshore is
three-dimensional produces a runup error equal to 51 per cent of the
fractional slope variability -- underestimates of 12 per cent on steep
cusp horns and overestimates of up to 38 per cent in embayments. A DEM
already resolves that variation, so there is no reason to discard it.

HOW THE SLOPE IS DEFINED: Stockdon define the foreshore slope as the
average over the region of significant swash activity, <eta> +/- S/2,
rather than over the whole profile. That definition is circular at
first pass -- the swash region depends on the wave conditions, which
is what the slope is needed to compute -- so this script offers both:

  * a fixed elevation band (--z-min/--z-max), the simple case; and
  * an iterative solve (--iterate), which computes slope over a first
    guess, evaluates the Stockdon swash height, re-centres the band on
    the resulting swash region, and repeats. Two or three passes
    converge.

The slope is computed along the DIRECTION OF STEEPEST DESCENT of the
surface, not along a grid axis. That distinction is not cosmetic here:
this shoreline runs obliquely to the UTM grid, so an easting-aligned
profile cuts the beach at an angle, overstating its width and
understating its slope. Profiles taken that way gave 1:15 and 1:19,
which are therefore lower bounds.

Usage:
    python3 foreshore_slope.py dem_intertidal --output slope
    python3 foreshore_slope.py dem_intertidal --output slope \\
        --z-min -0.5 --z-max 0.5 --iterate --wave-height 0.9 --wave-period 6.8
"""

import sys
import argparse
from pathlib import Path

import numpy as np

G = 9.81


def read_ascii_grid(path):
    path = Path(path)
    if not path.exists():
        sys.exit(f"not found: {path}")
    hdr = {}
    with open(path) as f:
        for _ in range(6):
            parts = f.readline().split()
            if len(parts) != 2:
                sys.exit(f"{path}: malformed header")
            hdr[parts[0].lower()] = float(parts[1])
    grid = np.loadtxt(path, skiprows=6)
    grid = np.where(grid == hdr.get("nodata_value", -9999.0), np.nan, grid)
    return np.flipud(grid), hdr          # flip: ASCII runs north to south


def write_ascii_grid(path, grid, hdr, nodata=-9999.0):
    with open(path, "w") as f:
        f.write(f"ncols {grid.shape[1]}\n")
        f.write(f"nrows {grid.shape[0]}\n")
        f.write(f"xllcorner {hdr['xllcorner']:.3f}\n")
        f.write(f"yllcorner {hdr['yllcorner']:.3f}\n")
        f.write(f"cellsize {hdr['cellsize']}\n")
        f.write(f"NODATA_value {nodata}\n")
        for row in np.flipud(grid):
            f.write(" ".join("%.6f" % (v if np.isfinite(v) else nodata) for v in row) + "\n")


def write_prj(path):
    """
    Writes the projection sidecar.

    The ASCII grid format carries corner coordinates and cell size but
    no coordinate system, so anything reading it has to be told the
    projection by hand -- an easy thing to get wrong or forget. A .prj
    beside the .asc closes that gap for QGIS, ArcGIS and GDAL alike.
    """
    wkt = ('PROJCS["WGS_1984_UTM_Zone_19N",GEOGCS["GCS_WGS_1984",'
           'DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],'
           'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],'
           'PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],'
           'PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-69.0],'
           'PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],'
           'UNIT["Meter",1.0]]')
    Path(path).write_text(wkt + "\n")


def surface_slope(dem, cell, window):
    """
    Slope magnitude along the direction of steepest descent.

    A plane is fitted by least squares to the elevations in a window
    about each cell; the slope is the magnitude of that plane's
    gradient. Fitting a plane rather than differencing neighbours is
    deliberate: the DEM has gaps and per-cell noise of order 0.23 m
    (its measured repeatability), and a two-cell difference would
    amplify both. The fit uses whatever cells are present, so it
    degrades gracefully at the edges of coverage rather than
    propagating NaN.

    Taking the magnitude of the gradient, rather than a component
    along a grid axis, matters because this shoreline runs obliquely
    to the UTM grid: an axis-aligned profile crosses the beach at an
    angle and understates its steepness.
    """
    rows, cols = dem.shape
    half = max(1, int(window // 2))
    out = np.full_like(dem, np.nan, dtype=float)

    # Local coordinates in metres, shared by every window.
    yy, xx = np.mgrid[-half:half + 1, -half:half + 1]
    xx = xx.astype(float) * cell
    yy = yy.astype(float) * cell

    for r in range(rows):
        r0, r1 = max(0, r - half), min(rows, r + half + 1)
        for c in range(cols):
            if not np.isfinite(dem[r, c]):
                continue
            c0, c1 = max(0, c - half), min(cols, c + half + 1)
            z = dem[r0:r1, c0:c1]
            ok = np.isfinite(z)
            # A plane has three unknowns; fewer than six points makes
            # the fit unstable against the DEM's own noise.
            if ok.sum() < 6:
                continue
            X = xx[r0 - (r - half):r1 - (r - half), c0 - (c - half):c1 - (c - half)][ok]
            Y = yy[r0 - (r - half):r1 - (r - half), c0 - (c - half):c1 - (c - half)][ok]
            Z = z[ok]
            A = np.column_stack([X, Y, np.ones(len(Z))])
            try:
                coef, *_ = np.linalg.lstsq(A, Z, rcond=None)
            except np.linalg.LinAlgError:
                continue
            out[r, c] = float(np.hypot(coef[0], coef[1]))
    return out


def stockdon_swash(H0, T0, beta):
    """Significant swash height, Stockdon et al. (2006) Eqs. 11-12."""
    L0 = G * T0 ** 2 / (2 * np.pi)
    Sinc = 0.75 * beta * np.sqrt(H0 * L0)
    SIG = 0.06 * np.sqrt(H0 * L0)
    return np.sqrt(Sinc ** 2 + SIG ** 2)


def stockdon_setup(H0, T0, beta):
    """Shoreline setup, Stockdon et al. (2006) Eq. 10."""
    L0 = G * T0 ** 2 / (2 * np.pi)
    return 0.35 * beta * np.sqrt(H0 * L0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dem_stem", help="Stem of the DEM, e.g. 'dem_intertidal' for "
                                     "dem_intertidal_dem.asc")
    ap.add_argument("--output", default=None, help="Output stem (default: <dem_stem>_slope)")
    ap.add_argument("--window", type=int, default=5,
                    help="Window width in cells for the plane fit (default 5). At 2 m cells "
                         "that is a 10 m footprint, comparable to the swash zone width here.")
    ap.add_argument("--z-min", type=float, default=None,
                    help="Restrict the slope band to elevations above this (m). Stockdon "
                         "define foreshore slope over the swash region, not the whole "
                         "profile; without --iterate this sets that band by hand.")
    ap.add_argument("--z-max", type=float, default=None)
    ap.add_argument("--iterate", action="store_true",
                    help="Solve the circular definition: slope determines the swash region, "
                         "which determines where slope should be measured. Requires "
                         "--wave-height and --wave-period.")
    ap.add_argument("--wave-height", type=float, help="Representative H0 (m), for --iterate")
    ap.add_argument("--wave-period", type=float, help="Representative T0 (s), for --iterate")
    ap.add_argument("--still-water", type=float, default=None,
                    help="Still-water level (m NAVD88) the swash region is centred on, for "
                         "--iterate. Setup and swash are measured from still water, not "
                         "from 0 m NAVD88; without this the band assumed still water at "
                         "0 m. Default: the middle of the DEM's elevation range, i.e. the "
                         "middle of the water levels that built it.")
    ap.add_argument("--passes", type=int, default=3)
    args = ap.parse_args()

    stem = args.output or (args.dem_stem + "_slope")
    dem, hdr = read_ascii_grid(args.dem_stem + "_dem.asc")
    cell = hdr["cellsize"]

    print("=" * 68)
    print("FORESHORE SLOPE")
    print("=" * 68)
    n_dem = int(np.isfinite(dem).sum())
    print(f"DEM cells         : {n_dem} at {cell} m")
    print(f"elevation range   : {np.nanmin(dem):+.2f} to {np.nanmax(dem):+.2f} m")

    band = dem.copy()
    if args.z_min is not None:
        band = np.where(band >= args.z_min, band, np.nan)
    if args.z_max is not None:
        band = np.where(band <= args.z_max, band, np.nan)
    if args.z_min is not None or args.z_max is not None:
        print(f"slope band        : {int(np.isfinite(band).sum())} cell(s) within the "
              f"requested elevation range")

    print(f"plane-fit window  : {args.window} cells ({args.window*cell:.0f} m)")
    slope = surface_slope(band, cell, args.window)

    if args.iterate:
        if args.wave_height is None or args.wave_period is None:
            sys.exit("--iterate needs --wave-height and --wave-period")
        swl = args.still_water
        if swl is None:
            swl = 0.5 * (float(np.nanmin(dem)) + float(np.nanmax(dem)))
        print()
        print("Iterating: slope sets the swash region, which sets where slope is measured.")
        print(f"  still-water level {swl:+.2f} m NAVD88"
              + ("" if args.still_water is not None else " (middle of the DEM range; "
                 "set --still-water)"))
        for i in range(args.passes):
            med = float(np.nanmedian(slope))
            if not np.isfinite(med) or med <= 0:
                print("  slope did not converge to a usable value -- stopping")
                break
            setup = stockdon_setup(args.wave_height, args.wave_period, med)
            S = stockdon_swash(args.wave_height, args.wave_period, med)
            lo, hi = swl + setup - S / 2.0, swl + setup + S / 2.0
            band = np.where((dem >= lo) & (dem <= hi), dem, np.nan)
            n_band = int(np.isfinite(band).sum())
            print(f"  pass {i+1}: slope 1:{1/med:.1f}  ->  swash region "
                  f"{lo:+.2f} to {hi:+.2f} m  ({n_band} cell(s))")
            if n_band < 20:
                print("    too few cells in that band -- keeping the previous slope")
                break
            slope = surface_slope(band, cell, args.window)

    valid = slope[np.isfinite(slope)]
    if len(valid) == 0:
        sys.exit("No slope could be computed. Widen the elevation band or the window.")

    print()
    print(f"slope cells       : {len(valid)}")
    print(f"  median          : {np.median(valid):.4f}  (1:{1/np.median(valid):.1f})")
    print(f"  p10 - p90       : {np.percentile(valid,10):.4f} - "
          f"{np.percentile(valid,90):.4f}  "
          f"(1:{1/np.percentile(valid,10):.0f} to 1:{1/np.percentile(valid,90):.0f})")

    # Alongshore variability is what Stockdon warn about: a single
    # averaged slope costs 51% of the fractional variability in runup
    # error, so the spread here is a direct measure of what using one
    # number would cost.
    spread = (np.percentile(valid, 90) - np.percentile(valid, 10)) / np.median(valid)
    print(f"  p10-p90 spread  : {100*spread:.0f}% of the median")
    print(f"  -> using a single averaged slope would cost roughly "
          f"{51*spread:.0f}% error in runup")
    print("     (Stockdon et al. 2006: runup error is 51% of fractional slope variability)")

    write_ascii_grid(stem + ".asc", slope, hdr)
    write_prj(stem + ".prj")
    write_prj(args.dem_stem + "_dem.prj")      # the DEM lacked one too
    print()
    print(f"wrote {stem}.asc   (foreshore slope, dimensionless)")
    print(f"wrote {stem}.prj   (UTM 19N projection)")
    print(f"wrote {args.dem_stem}_dem.prj   (the DEM had no projection defined)")

    if args.wave_height is not None and args.wave_period is not None:
        med = float(np.median(valid))
        L0 = G * args.wave_period ** 2 / (2 * np.pi)
        xi = med / np.sqrt(args.wave_height / L0)
        setup = stockdon_setup(args.wave_height, args.wave_period, med)
        S = stockdon_swash(args.wave_height, args.wave_period, med)
        R2 = 1.1 * (setup + S / 2)
        print()
        print(f"At H0={args.wave_height} m, T0={args.wave_period} s and the median slope:")
        print(f"  Iribarren xi0   : {xi:.2f}  "
              f"({'dissipative' if xi < 0.3 else 'reflective' if xi > 1.25 else 'intermediate'})")
        print(f"  setup <eta>     : {setup:.3f} m")
        print(f"  swash S         : {S:.3f} m")
        print(f"  R2 (2% exceed)  : {R2:.3f} m")
        print()
        print("  These bound the elevation correction the DEM currently omits. Which value")
        print("  applies depends on where within the swash zone the detector places the")
        print("  waterline -- setup alone if it tracks the time-mean shoreline, closer to R2")
        print("  if it tracks the upper swash limit. That has not been established for this")
        print("  detector and must be calibrated, not assumed.")


if __name__ == "__main__":
    main()

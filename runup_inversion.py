#!/usr/bin/env python3
"""
Wave Height By Inverting The Runup Relation
==============================================
Estimates offshore significant wave height from the elevation of the
detected waterline, using Stockdon et al. (2006) Eq. 19 solved
backwards for H0.

THE IDEA. The visible waterline is not the still-water line: it sits
at the runup limit, above still water by an amount that grows with
wave height and beach slope. That offset is normally a nuisance -- the
largest uncorrected term in the DEM error budget. Here it is the
signal. Given the still-water level from GNSS-IR or a tide record, the
local foreshore slope from the DEM, and a wave period, the height of
the waterline above still water determines H0.

    R2 - z0 = 1.1 [ 0.35 b sqrt(H0 L0)
                    + sqrt( H0 L0 (0.563 b^2 + 0.004) ) / 2 ]

with L0 = g T0^2 / 2pi. Every term on the right is proportional to
sqrt(H0), so the relation inverts in closed form -- no iteration, no
ambiguity about which root to take.

WHY IT IS WORTH HAVING ALONGSIDE THE NEURAL MODEL. Its errors come
from different places: waterline geometry and beach slope, rather than
image texture and surf-zone width. Two estimators that fail
independently can be combined usefully, which is not true of averaging
one estimator over time. It may also degrade more gracefully in
storms, where the CNN saturates because the surf zone becomes a
uniform white wash while the runup limit remains a detectable
boundary.

WHAT IT CANNOT DO. Stockdon's own scatter is 38 cm on R2, which maps
to a large relative error in H0 when waves are small. It needs a wave
period, which this project could not predict from imagery (R2 = 0.10)
and which must come from a buoy or model. And it inherits whatever
vertical error the waterline elevation carries. Treat it as a second
opinion, not a replacement.

Usage:
    python3 runup_inversion.py --contours contour_points.csv \\
        --slope dem_intertidal_slope.asc --period 8.0 \\
        --output h_from_runup.csv --compare owg_vs_adcp_scored.csv
"""

import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

G = 9.81


def read_ascii_grid(path):
    hdr = {}
    with open(path) as f:
        for _ in range(6):
            k, v = f.readline().split()
            hdr[k.lower()] = float(v)
    grid = np.loadtxt(path, skiprows=6)
    grid = np.where(grid == hdr.get("nodata_value", -9999.0), np.nan, grid)
    return np.flipud(grid), hdr


def slope_at(grid, hdr, easting, northing):
    """Nearest-cell slope, NaN outside the grid."""
    col = ((easting - hdr["xllcorner"]) / hdr["cellsize"]).astype(int)
    row = ((northing - hdr["yllcorner"]) / hdr["cellsize"]).astype(int)
    ok = ((row >= 0) & (row < grid.shape[0]) & (col >= 0) & (col < grid.shape[1]))
    out = np.full(len(easting), np.nan)
    out[ok] = grid[row[ok], col[ok]]
    return out


def h0_from_runup(dz, beta, T0):
    """
    Solve Stockdon Eq. 19 for H0 given the runup elevation above still
    water.

    Both terms scale as sqrt(H0 L0), so writing k for the bracket
    coefficient the relation is dz = 1.1 k sqrt(H0 L0), and

        H0 = (dz / (1.1 k))^2 / L0

    exactly. Negative dz -- a waterline below still water, which is
    physically impossible and indicates an error in the elevation or
    the datum -- returns NaN rather than a number.
    """
    L0 = G * np.asarray(T0, float) ** 2 / (2 * np.pi)
    beta = np.asarray(beta, float)
    k = 0.35 * beta + 0.5 * np.sqrt(0.563 * beta ** 2 + 0.004)
    dz = np.asarray(dz, float)
    out = np.full(dz.shape, np.nan)
    good = (dz > 0) & np.isfinite(beta) & (k > 0) & np.isfinite(L0)
    out[good] = (dz[good] / (1.1 * k[good])) ** 2 / L0[good]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contours", required=True,
                    help="Georectified contour points with easting, northing, elevation "
                         "and the still-water level used.")
    ap.add_argument("--slope", default=None,
                    help="Foreshore slope grid from foreshore_slope.py. Without it, "
                         "--fixed-slope is used everywhere.")
    ap.add_argument("--fixed-slope", type=float, default=None,
                    help="Single slope, e.g. 0.141 for 1:7.1 measured here.")
    ap.add_argument("--period", type=float, default=None,
                    help="Wave period T0 in seconds, if the contour file has none. "
                         "Period could not be predicted from imagery in this project, so it "
                         "must come from a buoy or wave model.")
    ap.add_argument("--elev-col", default="elevation_m")
    ap.add_argument("--still-col", default="water_level_m",
                    help="Still-water level at the time of capture.")
    ap.add_argument("--time-col", default="time_utc")
    ap.add_argument("--output", required=True)
    ap.add_argument("--compare", default=None,
                    help="Scored CNN predictions, to compare the two estimates.")
    args = ap.parse_args()

    df = pd.read_csv(args.contours)
    print("=" * 70)
    print("WAVE HEIGHT FROM RUNUP INVERSION")
    print("=" * 70)
    print(f"contour points    : {len(df)}")

    missing = [c for c in (args.elev_col, args.still_col) if c not in df.columns]
    if missing:
        print(f"  columns present : {list(df.columns)}")
        sys.exit(f"contour file lacks {missing}. Pass --elev-col / --still-col to match it.")

    # ---- slope ----------------------------------------------------------
    if args.slope and Path(args.slope).exists():
        grid, hdr = read_ascii_grid(args.slope)
        for c in ("easting", "northing"):
            if c not in df.columns:
                sys.exit(f"--slope needs an '{c}' column in the contour file")
        beta = slope_at(grid, hdr, df["easting"].to_numpy(), df["northing"].to_numpy())
        n_ok = int(np.isfinite(beta).sum())
        print(f"slope             : from grid, {n_ok} of {len(df)} points inside coverage")
        if args.fixed_slope:
            beta = np.where(np.isfinite(beta), beta, args.fixed_slope)
            print(f"                    {len(df)-n_ok} point(s) fell back to "
                  f"{args.fixed_slope:.3f}")
    elif args.fixed_slope:
        beta = np.full(len(df), args.fixed_slope)
        print(f"slope             : fixed at {args.fixed_slope:.3f} (1:{1/args.fixed_slope:.1f})")
    else:
        sys.exit("give --slope or --fixed-slope")

    # ---- period ---------------------------------------------------------
    if "wave_period_s" in df.columns and df["wave_period_s"].notna().any():
        T0 = df["wave_period_s"].to_numpy(float)
        print("period            : from the contour file")
    elif args.period:
        T0 = np.full(len(df), args.period)
        print(f"period            : fixed at {args.period:.1f} s")
        print("  A single period for the whole record is a strong assumption: H0 scales as")
        print("  1/T0^2 here, so a period wrong by 25% biases H0 by about 60%.")
    else:
        sys.exit("no wave period available; give --period or add wave_period_s")

    # ---- invert ----------------------------------------------------------
    dz = df[args.elev_col].to_numpy(float) - df[args.still_col].to_numpy(float)
    df["runup_above_still_m"] = dz
    df["H0_from_runup_m"] = h0_from_runup(dz, beta, T0)
    df["beta_used"] = beta

    good = np.isfinite(df["H0_from_runup_m"])
    print()
    print(f"runup above still : median {np.nanmedian(dz):.3f} m, "
          f"p10 {np.nanpercentile(dz,10):+.3f}, p90 {np.nanpercentile(dz,90):+.3f}")
    n_neg = int((dz <= 0).sum())
    if n_neg:
        print(f"  {n_neg} point(s) ({100*n_neg/len(df):.0f}%) sit at or below still water, "
              f"which is physically impossible.")
        print("  That points to an elevation or datum error rather than to small waves, and")
        print("  those points are dropped rather than forced to a value.")
    print(f"H0 estimated      : {int(good.sum())} point(s)")
    if good.sum():
        h = df.loc[good, "H0_from_runup_m"]
        print(f"                    median {h.median():.2f} m, "
              f"p10 {h.quantile(0.1):.2f}, p90 {h.quantile(0.9):.2f}")

    df.to_csv(args.output, index=False)
    print(f"wrote             : {args.output}")

    # ---- aggregate per capture and compare ------------------------------
    if args.compare and Path(args.compare).exists() and args.time_col in df.columns:
        per = (df[good].groupby(args.time_col)["H0_from_runup_m"]
               .median().rename("H0_runup").reset_index())
        cnn = pd.read_csv(args.compare)
        tcol = "time" if "time" in cnn.columns else args.time_col
        if tcol not in cnn.columns:
            print("comparison file has no time column; skipping comparison")
            return 0
        per["_t"] = pd.to_datetime(per[args.time_col]).dt.floor("30min")
        cnn["_t"] = pd.to_datetime(cnn[tcol]).dt.floor("30min")
        j = cnn.merge(per[["_t", "H0_runup"]], on="_t", how="inner")
        if len(j) < 10:
            print(f"only {len(j)} matching capture(s); too few to compare")
            return 0
        print()
        print(f"COMPARISON, {len(j)} capture(s) with both estimates")
        for name, col in (("CNN", "predicted"), ("runup inversion", "H0_runup")):
            e = j[col] - j["observed"]
            print(f"  {name:<16} RMSE {np.sqrt((e**2).mean()):.3f} m   "
                  f"bias {e.mean():+.3f} m")
        both = (j["predicted"] + j["H0_runup"]) / 2
        e = both - j["observed"]
        print(f"  {'mean of the two':<16} RMSE {np.sqrt((e**2).mean()):.3f} m   "
              f"bias {e.mean():+.3f} m")
        r = float(np.corrcoef(j["predicted"] - j["observed"],
                              j["H0_runup"] - j["observed"])[0, 1])
        print(f"  error correlation between the two: {r:+.2f}")
        if r < 0.4:
            print("  The two fail largely independently, which is what makes combining them")
            print("  worthwhile -- unlike averaging one estimator over time.")
        else:
            print("  The two err together, so combining them gains little.")
        j.to_csv(args.output.replace(".csv", "_compared.csv"), index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())

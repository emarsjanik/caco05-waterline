#!/usr/bin/env python3
"""
Georectified Waterlines On The Ground
========================================
Every waterline point after georectification, one dot per point,
coloured by its water level -- before any gridding, frame check or
blanking. Shows what the DEM is built from:

  * a sound beach gives nested contours, low water seaward, high water
    landward, each camera's lines parallel to the shore;
  * a bad frame shows as a line crossing the others;
  * lines of different levels landing on top of each other mean the
    DEM cells there will be blanked (spread > 0.5 m).

One panel per camera, plus both together. Frames the DEM's frame check
rejected can be listed with --rejected (the log lines work as is).

Usage:
    python3 plot_ground_points.py /mnt/I2Rgus_Data/Chelsea_calibration/contour_points_ground.csv \\
        /mnt/I2Rgus_Data/Chelsea_calibration/ground_points.png
"""

import sys
import argparse

import numpy as np
import pandas as pd


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("contour_csv", help="contour_points_ground.csv (with easting_utm19/northing_utm19)")
    ap.add_argument("output_png")
    ap.add_argument("--every", type=int, default=5, help="plot every Nth point per frame (default 5)")
    args = ap.parse_args()

    d = pd.read_csv(args.contour_csv)
    zcol = next((c for c in ("tide_elevation_navd88", "water_level_navd88") if c in d.columns), None)
    if zcol is None or "easting_utm19" not in d.columns:
        sys.exit("Need easting_utm19, northing_utm19 and tide_elevation_navd88 columns.")
    d = d.dropna(subset=["easting_utm19", "northing_utm19"])
    d = d[d.groupby("source_file").cumcount() % args.every == 0]
    lo, hi = d[zcol].min(), d[zcol].max()

    cams = sorted(d["camera"].unique())
    fig, axes = plt.subplots(1, len(cams) + 1, figsize=(6 * (len(cams) + 1), 9), squeeze=False)
    for ax, cam in zip(axes[0], cams + ["both"]):
        sub = d if cam == "both" else d[d["camera"] == cam]
        sub = sub.sort_values(zcol)
        sc = ax.scatter(sub["easting_utm19"], sub["northing_utm19"], c=sub[zcol], s=2,
                        cmap="turbo", vmin=lo, vmax=hi, linewidths=0)
        ax.set_aspect("equal")
        ax.set_title(f"{cam}: {sub['source_file'].nunique()} frames, {len(sub)} points shown")
        ax.set_xlabel("easting (m, UTM 19N)")
        ax.set_ylabel("northing (m, UTM 19N)")
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.grid(alpha=0.3)
    fig.colorbar(sc, ax=axes[0].tolist(), shrink=0.7, label="water level (m NAVD88)")
    fig.suptitle("Georectified waterlines before gridding (low water should be seaward, "
                 "high water landward, lines not crossing)")
    fig.savefig(args.output_png, dpi=110, bbox_inches="tight")
    print(f"wrote {args.output_png}")

    # Numbers the picture cannot give: extent of each camera's points.
    for cam in cams:
        s = d[d["camera"] == cam]
        print(f"{cam}: easting {s.easting_utm19.min():.0f}-{s.easting_utm19.max():.0f}, "
              f"northing {s.northing_utm19.min():.0f}-{s.northing_utm19.max():.0f}, "
              f"levels {s[zcol].min():+.2f} to {s[zcol].max():+.2f} m")


if __name__ == "__main__":
    main()

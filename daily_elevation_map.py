#!/usr/bin/env python3
"""
Daily Beach Elevation Map
----------------------------
Takes one day's detected waterlines from contour_points.csv and draws
them all onto a single background image (by default the ~12:15
capture), shaded by water elevation.

WHAT THE FILLED BANDS MEAN -- read this before interpreting the figure:
  Over a day the tide crosses each elevation TWICE, once rising and
  once falling. Where two or more detected waterlines fall in the same
  elevation bin, the area between them is filled with that elevation's
  colour. On a stable beach those lines should coincide, so:

      band width  ~=  REPEATABILITY of the measurement at that
                      elevation, NOT a morphological feature.

  A wide band means the two crossings disagreed. That can come from
  detection error, from wave runup differing between rising and
  falling tide, or from genuine morphological change during the day --
  and this figure alone cannot separate those. Treat a wide band as
  "uncertain here", not as "the beach is this shape here".

  Bands are drawn semi-transparent so the underlying image stays
  visible, which lets you judge by eye whether a line sits on the real
  water's edge.

Usage:
    python3 daily_elevation_map.py <contour_points.csv> <image_dir> <camera> <output.png>
        [--date YYYY-MM-DD] [--background-hour 12] [--background-minute 15]
        [--elevation-bin 0.05] [--line-alpha 0.85] [--fill-alpha 0.30]

Example:
    python3 daily_elevation_map.py contour_points.csv \\
        /mnt/I2Rgus_Data/waterline/archive/images_timex c1 \\
        elevation_map_c1_20260911.png --date 2026-09-11
"""

import sys
import csv
import argparse
from collections import defaultdict
from datetime import datetime, timezone, timedelta, date as date_cls
from pathlib import Path

import numpy as np
import cv2

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D


ELEVATION_COLUMNS = ("tide_elevation_navd88", "tide_elevation")


def load_contours(path, camera, date_filter=None):
    """
    Groups contour points by source frame. Returns a dict keyed by
    source_file, each with sorted column/row arrays, the elevation,
    and the capture time.
    """
    frames = defaultdict(lambda: {"columns": [], "rows": [],
                                  "elevation": None, "capture": None})

    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        elev_col = next((c for c in ELEVATION_COLUMNS if c in reader.fieldnames), None)
        if elev_col is None:
            raise KeyError(
                f"No elevation column in {path}. Expected one of {ELEVATION_COLUMNS}; "
                f"found: {reader.fieldnames}")

        for row in reader:
            if row.get("camera") != camera:
                continue
            capture = row["capture_time_utc"]
            if date_filter and not capture.startswith(date_filter):
                continue
            key = row["source_file"]
            frames[key]["columns"].append(float(row["pixel_column"]))
            frames[key]["rows"].append(float(row["pixel_row"]))
            frames[key]["elevation"] = float(row[elev_col])
            frames[key]["capture"] = capture

    for data in frames.values():
        order = np.argsort(data["columns"])
        data["columns"] = np.array(data["columns"])[order]
        data["rows"] = np.array(data["rows"])[order]

    return dict(frames), elev_col


def pick_background(frames, image_dir, target_hour, target_minute):
    """
    Chooses the frame closest to the requested time of day and returns
    its image path. Falls back to any readable image from the day if
    the preferred one is missing from disk.
    """
    image_dir = Path(image_dir)
    target = target_hour * 60 + target_minute

    def minutes_from_target(item):
        capture = datetime.fromisoformat(item[1]["capture"])
        return abs((capture.hour * 60 + capture.minute) - target)

    for key, data in sorted(frames.items(), key=minutes_from_target):
        candidate = image_dir / (key + ".jpg")
        if candidate.exists():
            return candidate, data["capture"]
        matches = list(image_dir.glob(key + "*"))
        if matches:
            return matches[0], data["capture"]
    return None, None


def main():
    parser = argparse.ArgumentParser(
        description="Draw one day's waterlines on a single image, shaded by elevation.")
    parser.add_argument("contour_csv")
    parser.add_argument("image_dir")
    parser.add_argument("camera", choices=["c1", "c2"])
    parser.add_argument("output_png")
    parser.add_argument("--date", default=None,
                        help="Single YYYY-MM-DD to plot. Shorthand for --start-date X "
                             "--end-date X.")
    parser.add_argument("--days", type=int, default=1,
                        help="Number of days to include, ending at --end-date (or at the most "
                             "recent date present). Use --days 7 for a rolling week. Ignored "
                             "if --date or --start-date is given.")
    parser.add_argument("--start-date", default=None, help="First YYYY-MM-DD to include.")
    parser.add_argument("--end-date", default=None, help="Last YYYY-MM-DD to include.")
    parser.add_argument("--background-date", default=None,
                        help="Which day's image to use as the backdrop. Default: the last day "
                             "in the range, so the picture sits on the most recent view of the "
                             "beach.")
    parser.add_argument("--background-hour", type=int, default=12)
    parser.add_argument("--background-minute", type=int, default=15)
    parser.add_argument("--elevation-bin", type=float, default=0.10,
                        help="Elevation bin width in metres for grouping lines treated as the "
                             "same level (default 0.10). Lines never match exactly so some "
                             "tolerance is needed. A day of 30-minute captures typically "
                             "spaces levels ~0.1 m apart, so a 0.05 m bin rarely catches the "
                             "rising and falling crossing of the same level and little gets "
                             "filled; too wide and genuinely different levels merge, making "
                             "bands look artificially large. Check the bin-occupancy summary "
                             "printed at the end and adjust.")
    parser.add_argument("--max-gap-columns", type=float, default=40,
                        help="Leave a break in the drawn line wherever the nearest column with "
                             "actual data is further than this (default 40). Prevents the plot "
                             "from implying a continuous waterline across regions the detector "
                             "flagged as having no signal. Set 0 to interpolate across "
                             "everything (the old behaviour).")
    parser.add_argument("--line-alpha", type=float, default=0.85)
    parser.add_argument("--line-width", type=float, default=None,
                        help="Line width. Default scales with how many lines are drawn, so a "
                             "week's worth stays legible.")
    parser.add_argument("--fill-alpha", type=float, default=0.30)
    parser.add_argument("--colormap", default="turbo")
    parser.add_argument("--dpi", type=int, default=130)
    args = parser.parse_args()

    frames, elev_col = load_contours(args.contour_csv, args.camera)
    if not frames:
        print(f"No contour points for camera '{args.camera}' in {args.contour_csv}")
        sys.exit(1)

    all_dates = sorted({d["capture"][:10] for d in frames.values()})

    if args.date:
        start_date = end_date = args.date
    elif args.start_date or args.end_date:
        start_date = args.start_date or all_dates[0]
        end_date = args.end_date or all_dates[-1]
    else:
        end_date = all_dates[-1]
        start_date = (date_cls.fromisoformat(end_date)
                      - timedelta(days=max(1, args.days) - 1)).isoformat()

    frames = {k: v for k, v in frames.items()
              if start_date <= v["capture"][:10] <= end_date}
    if not frames:
        print(f"No frames for camera '{args.camera}' between {start_date} and {end_date}.")
        print(f"Dates available: {all_dates}")
        sys.exit(1)

    dates_used = sorted({d["capture"][:10] for d in frames.values()})
    date_label = (dates_used[0] if len(dates_used) == 1
                  else f"{dates_used[0]} to {dates_used[-1]}  ({len(dates_used)} days)")

    # The backdrop comes from ONE day -- by default the most recent in
    # range, so accumulated lines are drawn over the latest view of the
    # beach rather than a week-old one.
    bg_day = args.background_date or dates_used[-1]
    bg_candidates = {k: v for k, v in frames.items()
                     if v["capture"].startswith(bg_day)} or frames
    bg_path, bg_capture = pick_background(
        bg_candidates, args.image_dir, args.background_hour, args.background_minute)
    if bg_path is None:
        print(f"No background image found in {args.image_dir} for any frame on {bg_day}.")
        print("Frames needing an image:")
        for k in sorted(frames):
            print(f"  {k}.jpg")
        sys.exit(1)

    image = cv2.imread(str(bg_path))
    if image is None:
        print(f"Could not read background image: {bg_path}")
        sys.exit(1)
    image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
    height, width = image_rgb.shape[:2]

    elevations = np.array([d["elevation"] for d in frames.values()])
    norm = Normalize(vmin=elevations.min(), vmax=elevations.max())
    colormap = matplotlib.colormaps[args.colormap]

    # Resample every line onto a common column grid so lines in the
    # same elevation bin can be filled between directly.
    grid = np.arange(width, dtype=float)
    resampled = {}
    for key, data in frames.items():
        if len(data["columns"]) < 2:
            continue
        y = np.interp(grid, data["columns"], data["rows"], left=np.nan, right=np.nan)

        # Do NOT interpolate across wide gaps. Columns absent from the
        # contour file were dropped as no-signal, so bridging them
        # invents a waterline where the detector explicitly found none.
        # On real C2 data a frame reduced to 13 of 447 far-field
        # columns was still being drawn as a continuous line spanning
        # the whole frame. Gaps wider than --max-gap-columns are left
        # as NaN, which matplotlib renders as a break.
        if args.max_gap_columns > 0 and len(data["columns"]) > 1:
            present = np.zeros(width, dtype=bool)
            idx = np.clip(data["columns"].astype(int), 0, width - 1)
            present[idx] = True
            # Distance to the nearest column that actually has data.
            pos = np.where(present)[0]
            if pos.size:
                nearest = np.abs(grid[:, None] - pos[None, :]).min(axis=1) \
                    if pos.size <= 4000 else None
                if nearest is None:
                    # Large series: compute distance by forward/backward fill.
                    fwd = np.full(width, np.inf); last = -np.inf
                    for i in range(width):
                        if present[i]: last = i
                        fwd[i] = i - last
                    bwd = np.full(width, np.inf); nxt = np.inf
                    for i in range(width - 1, -1, -1):
                        if present[i]: nxt = i
                        bwd[i] = nxt - i
                    nearest = np.minimum(fwd, bwd)
                y[nearest > args.max_gap_columns] = np.nan
        resampled[key] = y

    # Group by elevation bin.
    bins = defaultdict(list)
    for key, data in frames.items():
        if key not in resampled:
            continue
        bin_index = int(round(data["elevation"] / args.elevation_bin))
        bins[bin_index].append(key)

    # With a week of captures this can be ~120 lines. Drawn at
    # single-day weight they stack into an opaque mass and the
    # background disappears, which defeats the point of shading. Scale
    # down automatically unless the user has overridden.
    n_lines = len(resampled)
    line_alpha = args.line_alpha
    line_width = args.line_width
    if args.line_width is None:
        line_width = 1.6 if n_lines <= 25 else (1.1 if n_lines <= 60 else 0.8)
    fill_alpha = args.fill_alpha
    if n_lines > 100:
        # A week of captures saturates the intertidal zone into a solid
        # ramp and the photo behind it disappears. Keep it translucent
        # so the line positions can still be checked against the image.
        line_alpha = min(line_alpha, 0.35)
        fill_alpha = min(fill_alpha, 0.12)
    elif n_lines > 60:
        line_alpha = min(line_alpha, 0.45)
        fill_alpha = min(fill_alpha, 0.18)
    elif n_lines > 25:
        line_alpha = min(line_alpha, 0.70)
        fill_alpha = min(fill_alpha, 0.25)

    fig, ax = plt.subplots(figsize=(width / args.dpi, height / args.dpi), dpi=args.dpi)
    ax.imshow(image_rgb)

    filled_bins = 0
    for bin_index, keys in sorted(bins.items(), key=lambda kv: frames[kv[1][0]]["elevation"]):
        bin_elev = float(np.mean([frames[k]["elevation"] for k in keys]))
        colour = colormap(norm(bin_elev))

        if len(keys) >= 2:
            stack = np.vstack([resampled[k] for k in keys])
            # Columns covered by no line in this bin are all-NaN; that
            # is expected, so silence the warning rather than let it
            # clutter every run.
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                lo = np.nanmin(stack, axis=0)
                hi = np.nanmax(stack, axis=0)
            valid = ~np.isnan(lo) & ~np.isnan(hi)
            if valid.any():
                ax.fill_between(grid[valid], lo[valid], hi[valid],
                                color=colour, alpha=fill_alpha, linewidth=0)
                filled_bins += 1

        for k in keys:
            y = resampled[k]
            valid = ~np.isnan(y)
            ax.plot(grid[valid], y[valid], color=colour,
                    alpha=line_alpha, linewidth=line_width)

    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.axis("off")

    ax.set_title(
        f"{args.camera.upper()}  {date_label}   {len(frames)} waterlines, "
        f"{elevations.min():+.2f} to {elevations.max():+.2f} m NAVD88\n"
        f"background: {bg_day} {bg_capture[11:16]} UTC   |   "
        + ("shaded bands = spread between same-elevation crossings "
           "(repeatability, not morphology)"
           if len(dates_used) == 1 else
           f"shaded bands = spread across {len(dates_used)} days at each elevation "
           "(repeatability + real shoreline movement)"),
        fontsize=9)

    scalar_map = matplotlib.cm.ScalarMappable(cmap=colormap, norm=norm)
    scalar_map.set_array([])
    cbar = fig.colorbar(scalar_map, ax=ax, fraction=0.030, pad=0.015)
    cbar.set_label("water elevation (m, NAVD88)")

    fig.tight_layout()
    fig.savefig(args.output_png, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    print(f"Date range        : {date_label}")
    print(f"Camera            : {args.camera}")
    print(f"Waterlines drawn  : {len(resampled)}")
    print(f"Elevation range   : {elevations.min():+.3f} to {elevations.max():+.3f} m NAVD88")
    print(f"Elevation bins    : {len(bins)}  ({filled_bins} had 2+ crossings and were filled)")
    print(f"Background image  : {bg_path.name}  ({bg_capture} UTC)")
    print(f"Saved             : {args.output_png}")

    if len(dates_used) > 1:
        print("Lines per day     :")
        for d in dates_used:
            n = sum(1 for v in frames.values() if v["capture"].startswith(d))
            print(f"    {d}: {n}")

    occupancy = sorted((len(v) for v in bins.values()), reverse=True)
    print(f"Lines per bin     : {occupancy}")
    print()
    if len(dates_used) == 1:
        print("Reminder: filled band width = spread between repeat crossings of the same")
        print("elevation within ONE day. It measures repeatability (detection error, runup")
        print("differences, or change during the day) -- not a morphology feature.")
    else:
        print(f"Reminder: this spans {len(dates_used)} days, so a filled band now mixes TWO")
        print("things: measurement repeatability, AND genuine movement of the shoreline at")
        print("that elevation over the period. A band that widens across the week is the")
        print("interesting case, but this figure alone cannot separate real change from")
        print("detection scatter -- compare against the single-day bands to judge which.")

    if filled_bins == 0:
        print()
        print("NOTE: no elevation bin contained two or more waterlines, so nothing was")
        print("filled. That happens when the day's captures never revisit the same level")
        print(f"within +/-{args.elevation_bin} m -- try a wider --elevation-bin, or a day")
        print("whose captures span both a rising and a falling tide.")


if __name__ == "__main__":
    main()

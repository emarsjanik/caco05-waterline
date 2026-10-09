#!/usr/bin/env python3
"""
Daily Beach Elevation Map
----------------------------
Takes one day's detected waterlines from contour_points.csv and draws
them all onto a single background image (by default the ~12:15
capture), shaded by the elevation each line marks on the beach.

WHICH ELEVATION: the one the consistency filter and the DEM use --
beach_elevation_navd88, the water level plus wave setup, when the
contour file has it (the cron applies setup, SETUP_COEF), else the
water level (tide_elevation_navd88), row by row the way
dem_from_contours.py does. Binned and coloured by the water level alone,
one 0.1 m bin held lines whose real levels differed by up to ~0.5 m on
stormy days (setup 0.15-0.88 m a frame in the station archive), the band
stretched across them, and the colours disagreed with the filter's report
and the DEM (30-day synthetic station through the cron, Oct 2026: c1's
7-day bands covered 629k px2 binned by water level, 276k px2 by water
level + setup; a line coloured +1.05 m here was +1.34 m in the report).
The colour bar says which elevation it is.

WHAT THE FILLED BANDS MEAN -- read this before interpreting the figure:
  Lines are grouped into elevation bins (--elevation-bin, 0.10 m; the
  elevation as above). In
  each image column covered by at least 3 lines of a bin that AGREE,
  the band spans the 16th to 84th percentile of their rows (about
  +/-1 sigma). Lines agree in a column when they form a group of at
  least 3: sorted by row, a group breaks wherever the next line is
  further than 12 px, or 4x the usual spacing of that bin's lines in
  that column if larger -- but never further than the rows half a
  metre of water level moves the line in that column (from all the
  map's lines). On a stable beach the lines of one bin should
  coincide, so:

      band width  ~=  REPEATABILITY of the measurement at that
                      elevation, NOT a morphological feature.

  Over several days the band also takes in real movement of the
  shoreline at that elevation: lines that moved together form a group
  of their own, or stay linked to the others, so they stay in, and the
  band spans both positions.

  A stray line -- a detection error, a line out on the water -- is
  still drawn, but cannot widen the band; nor can two or three strays
  that happen to lie close to each other. Until Oct 2026 the band ran
  from the lowest to the highest line in the bin (nanmin/nanmax), so
  one stray line turned it into a large translucent rectangle reaching
  that line, with vertical edges where the stray line started and
  stopped, and the picture looked far worse than the data. Where fewer
  than 3 lines agree there is no band: two lines cannot tell which of
  them is off. Gaps of fewer than 20 columns between stretches of band
  are bridged (no notch where a line or two drop out for a few
  columns), the edges are smoothed over 15 columns, stretches shorter
  than 20 columns are not drawn, and each stretch narrows to a point at
  its ends by at most 2 px per column on each side -- a tall band over
  a proportionally long stretch -- instead of stopping in a wall.

  A wide band means the lines that agree still disagree by that much:
  detection error, wave runup differing between rising and falling
  tide, or genuine change -- this figure alone cannot separate those.
  Treat a wide band as "uncertain here", not as "the beach is this
  shape here".

  Bands are drawn semi-transparent so the underlying image stays
  visible, which lets you judge by eye whether a line sits on the real
  water's edge.

  The cron draws these maps from the consistency-filtered contours
  (waterline_consistency.py, contour_points_timex_qc.csv), so lines out
  of order with the rest are already removed; the file used is named
  at the foot of the figure.

Usage:
    python3 daily_elevation_map.py <contour_points.csv> <image_dir> <camera> <output.png>
        [--date YYYY-MM-DD] [--background-hour 12] [--background-minute 15]
        [--elevation-bin 0.05] [--line-alpha 0.85] [--fill-alpha 0.30]
    python3 daily_elevation_map.py --self-test     (the band's shape, band_self_test())

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
BEACH_COLUMN = "beach_elevation_navd88"       # water level + wave setup (--setup-coef)


def load_contours(path, camera, date_filter=None, only=None):
    """
    Groups contour points by source frame. Returns a dict keyed by
    source_file, each with sorted column/row arrays, the elevation,
    whether that elevation includes wave setup, and the capture time;
    and the name of the water-level column.

    The elevation is beach_elevation_navd88 (water level + wave setup)
    where the file has it and the row's value is not blank, else the
    water level -- the same choice, row by row, as dem_from_contours.py,
    georectify.py and waterline_consistency.py make.
    """
    frames = defaultdict(lambda: {"columns": [], "rows": [],
                                  "elevation": None, "setup": False, "capture": None})

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
            if only is not None and key not in only:
                continue
            frames[key]["columns"].append(float(row["pixel_column"]))
            frames[key]["rows"].append(float(row["pixel_row"]))
            beach = row.get(BEACH_COLUMN)
            frames[key]["elevation"] = float(beach or row[elev_col])
            # setup_correction_m is blank for a frame with no wave record:
            # its beach_elevation_navd88 is then the water level itself
            frames[key]["setup"] = bool(beach) and bool(row.get("setup_correction_m"))
            frames[key]["capture"] = capture

    for data in frames.values():
        order = np.argsort(data["columns"])
        data["columns"] = np.array(data["columns"])[order]
        data["rows"] = np.array(data["rows"])[order]

    return dict(frames), elev_col


def agreeing_band(stack, min_lines=3, alone_px=12.0, alone_factor=4.0,
                  smooth_columns=15, max_edge_slope=2.0, min_run=20,
                  rows_per_m=None, link_m=0.5):
    """
    Band of the lines that agree, per column of `stack` (lines x columns,
    NaN where a line has no data). Returns (lo, hi, valid).

    Lines agree in a column when they form a GROUP of at least
    `min_lines`: sorted by row, a group breaks wherever the gap to the next
    line exceeds the link distance, max(`alone_px`, `alone_factor` x the
    column's median nearest-neighbour distance), capped at `link_m` (0.5 m)
    of water level in that column's rows (`rows_per_m`, per column, from all
    the map's lines; NaN or None: no cap). Smaller groups are left out of
    the band there. A stray line -- out on the water, or on the wrong
    feature -- is a group of one, so it can never stretch the band; nor
    can two strays lying close together (a group of two). Lines that move
    together, as after real beach change over a multi-day window, are a
    group of their own, or stay linked to the others: they stay in, and
    the band then shows that movement, as the map's title says. (A
    median/MAD rule was tried first: with two groups of lines it switched
    between showing both and only the larger one from column to column as
    lines started and stopped, which drew a blotchy band. Groups are kept
    or not on their own size, so nothing switches.)

    WHY groups and the cap: the earlier rule left out only a line with no
    other line near it. On the station's 7-day c2 map of 29 Sep - 5 Oct
    2026 two -0.4 m lines out on the water lay 1 px apart, so neither was
    alone, and the band ran 100 px out from the shore lines to them; in
    the -0.5 m bin one stray sat 115 px from two shore lines 31 px apart,
    under 4x their spacing. Together they drew a translucent 'bottle' on
    the water. 115 px there is ~1-2 m of water level; half a metre is
    already more than a 0.1 m bin, setup scatter and a week of real
    change put between neighbouring lines.

    The band is the 16th-84th percentile of the remaining rows -- about
    +/-1 sigma -- drawn only where at least `min_lines` remain. The 12 px
    floor is a few times the detector's ~2 px repeatability and below the
    spread of one 0.1 m bin in the near field (10-25 px), so honest lines at
    the edge of a bin are not mistaken for strays.

    Shape, in this order:
      * gaps shorter than `min_run` columns between two runs are bridged,
        the edges joined in straight lines across: where one or two
        columns drop below `min_lines` (a line starting, a detector gap)
        the band used to fall to nothing and back, a full-height V-notch
        (c2, bin +1.4 m on the synthetic week: a 3-column gap in a 239 px
        band);
      * edges smoothed over `smooth_columns`;
      * runs shorter than `min_run` columns dropped;
      * each edge moves away from the band's core by at most
        `max_edge_slope` px per column: from a point at each end, so a
        run narrows to a point with its half-width shrinking by at most
        2 px per column -- a 240 px band over 60 columns, a 20 px one
        over 5 -- and inside a run, where a group of lines starts or
        stops part way along, so the edge eases out instead of stepping.
        A fixed 15-column taper left tall bands ending in near-vertical
        walls (80% of the height within 14 columns), and steps where
        lines joined, which read as translucent rectangles. 2 px per
        column is several times the waterlines' own slope in the photo
        (c1 ~0.1, c2 up to ~0.6 px per column), so it bends the band's
        edges, not the lines it follows.
    """
    import warnings
    n_lines, width = stack.shape
    srt = np.sort(stack, axis=0)                     # NaN sort to the end
    d_prev = np.full_like(srt, np.inf); d_next = np.full_like(srt, np.inf)
    if n_lines > 1:
        gaps = np.diff(srt, axis=0)
        d_prev[1:] = gaps; d_next[:-1] = gaps
    nn = np.fmin(np.where(np.isnan(d_prev), np.inf, d_prev),
                 np.where(np.isnan(d_next), np.inf, d_next))
    nn[np.isnan(srt)] = np.nan
    with warnings.catch_warnings():
        # all-NaN columns are expected (no line of the bin there)
        warnings.simplefilter("ignore", category=RuntimeWarning)
        typical = np.nanmedian(np.where(np.isfinite(nn), nn, np.nan), axis=0)
        limit = alone_factor * np.nan_to_num(typical, nan=alone_px)
        if rows_per_m is not None:
            cap = link_m * np.abs(np.asarray(rows_per_m, dtype=float))
            limit = np.where(np.isfinite(cap), np.minimum(limit, cap), limit)
        limit = np.maximum(alone_px, limit)
        # groups: consecutive sorted rows no further apart than the limit
        idx = np.broadcast_to(np.arange(n_lines)[:, None], srt.shape)
        brk = np.ones(srt.shape, bool)
        if n_lines > 1:
            brk[1:] = ~(np.diff(srt, axis=0) <= limit)      # a NaN gap breaks too
        ends = np.ones(srt.shape, bool)
        ends[:-1] = brk[1:]
        first = np.maximum.accumulate(np.where(brk, idx, 0), axis=0)
        last = np.minimum.accumulate(np.where(ends, idx, n_lines - 1)[::-1], axis=0)[::-1]
        size = last - first + 1
        kept = np.where(np.isfinite(srt) & (size >= min_lines), srt, np.nan)
        n = np.isfinite(kept).sum(axis=0)
        lo = np.nanpercentile(kept, 16, axis=0)
        hi = np.nanpercentile(kept, 84, axis=0)
    valid = n >= min_lines
    if valid.any():
        # bridge short gaps between runs, edges joined straight across
        edges = np.flatnonzero(np.diff(np.r_[0, valid.astype(int), 0]))
        for e, s2 in zip(edges[1::2][:-1], edges[::2][1:]):     # gap = columns e .. s2-1
            if s2 - e < min_run:
                t = (np.arange(e, s2) - (e - 1)) / float(s2 - (e - 1))
                lo[e:s2] = lo[e - 1] + t * (lo[s2] - lo[e - 1])
                hi[e:s2] = hi[e - 1] + t * (hi[s2] - hi[e - 1])
                valid[e:s2] = True
    if smooth_columns > 1 and valid.any():
        kernel = np.ones(int(smooth_columns))
        weight = np.convolve(valid.astype(float), kernel, mode="same")
        with np.errstate(invalid="ignore", divide="ignore"):
            lo = np.convolve(np.where(valid, lo, 0.0), kernel, mode="same") / weight
            hi = np.convolve(np.where(valid, hi, 0.0), kernel, mode="same") / weight
    if valid.any():
        edges = np.flatnonzero(np.diff(np.r_[0, valid.astype(int), 0]))
        lo, hi = lo.copy(), hi.copy()
        for a, b in zip(edges[::2], edges[1::2]):
            if b - a < min_run:
                valid[a:b] = False
                continue
            # Each edge may move AWAY from the band's core by at most
            # max_edge_slope px per column -- from a point just beyond each
            # end, and from every neighbouring column (a cumulative min/max
            # of edge +/- slope x distance, each way). Edges only ever move
            # towards the core: no band is drawn where there was none.
            x = max_edge_slope * (np.arange(b - a) + 0.5)
            m0, m1 = 0.5 * (lo[a] + hi[a]), 0.5 * (lo[b - 1] + hi[b - 1])
            h, l = hi[a:b], lo[a:b]
            h_lim = np.minimum.reduce([
                h, m0 + x, m1 + x[::-1],
                x + np.minimum.accumulate(h - x),
                (x + np.minimum.accumulate(h[::-1] - x))[::-1]])
            l_lim = np.maximum.reduce([
                l, m0 - x, m1 - x[::-1],
                -x + np.maximum.accumulate(l + x),
                (-x + np.maximum.accumulate(l[::-1] + x))[::-1]])
            cross = l_lim > h_lim                     # the line itself moved faster than that
            c = 0.5 * (l_lim + h_lim)
            lo[a:b] = np.where(cross, c, l_lim)
            hi[a:b] = np.where(cross, c, h_lim)
    return lo, hi, valid


def rows_per_metre(stack, elevations, min_lines=8, min_range=0.3):
    """
    How many rows the waterline moves per metre of water level, in each
    column: the least-squares slope of row against elevation over all the
    lines that reach it (NaN where fewer than `min_lines` do, or their
    elevations span less than `min_range` m). A few strays among the ~100
    lines of a week barely move it. agreeing_band() uses it to say how far
    apart, in metres, two lines of one bin are.
    """
    ok = np.isfinite(stack)
    z = np.where(ok, np.asarray(elevations, dtype=float)[:, None], np.nan)
    n = ok.sum(axis=0)
    out = np.full(stack.shape[1], np.nan)
    use = n >= min_lines
    if not use.any():
        return out
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        zm = np.nanmean(z[:, use], axis=0)
        ym = np.nanmean(stack[:, use], axis=0)
        dz = z[:, use] - zm
        sxx = np.nansum(dz * dz, axis=0)
        sxy = np.nansum(dz * (stack[:, use] - ym), axis=0)
        zr = np.nanmax(z[:, use], axis=0) - np.nanmin(z[:, use], axis=0)
        b = np.where((sxx > 0) & (zr >= min_range), sxy / np.where(sxx > 0, sxx, 1.0), np.nan)
    out[use] = np.abs(b)
    return out


def crop_bounds_for(camera, image_height):
    """
    Returns (top, bottom) rows of the region the detector analyses, so
    clarity is judged on the beach rather than on sky and dune. Reads
    the real crop from the detector when it can be imported; otherwise
    falls back to the middle half of the frame, which is cruder but
    still excludes most sky.
    """
    try:
        import waterline_detector_v5 as detector
        key = {"c1": "CACO05_C1", "c2": "CACO05_C2"}[camera]
        profile = detector.CAMERAS[key]
        return (int(profile.crop_top * image_height),
                int(profile.crop_bottom * image_height))
    except Exception:
        return (int(0.25 * image_height), int(0.75 * image_height))


def resolve_image(image_dir, key):
    candidate = image_dir / (key + ".jpg")
    if candidate.exists():
        return candidate
    matches = list(image_dir.glob(key + "*"))
    return matches[0] if matches else None


def pick_background(frames, image_dir, target_hour, target_minute,
                    camera=None, select="clearest"):
    """
    Chooses the backdrop image.

    select="clearest" (default) picks the frame with the highest
    contrast INSIDE the detector's crop region. select="time" picks the
    frame nearest a requested time of day.

    Clearest is the default because time-of-day selection repeatedly
    produced unusable backdrops: a midday frame on a foggy day is still
    fog, and three separate maps came out drawn over a grey wash in
    which no waterline was visible, making the overlay impossible to
    check by eye. Contrast inside the crop separates these decisively
    -- measured on this station, fogged frames scored 5-7 while clear
    ones scored 45-58, an order of magnitude apart.

    Contrast is measured on the crop, not the whole frame, because a
    fogged beach under a bright sky can still show high whole-image
    contrast: one fogged frame scored 31.8 globally and 5.3 in-band.

    This affects ONLY which photo sits underneath. The contours drawn
    on top are unchanged.
    """
    image_dir = Path(image_dir)

    if select == "time":
        target = target_hour * 60 + target_minute

        def minutes_from_target(item):
            capture = datetime.fromisoformat(item[1]["capture"])
            return abs((capture.hour * 60 + capture.minute) - target)

        for key, data in sorted(frames.items(), key=minutes_from_target):
            path = resolve_image(image_dir, key)
            if path is not None:
                return path, data["capture"], None
        return None, None, None

    best = None
    for key, data in frames.items():
        path = resolve_image(image_dir, key)
        if path is None:
            continue
        image = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
        if image is None:
            continue
        top, bottom = crop_bounds_for(camera, image.shape[0])
        top = max(0, min(top, image.shape[0] - 1))
        bottom = max(top + 1, min(bottom, image.shape[0]))
        score = float(image[top:bottom, :].std())
        if best is None or score > best[0]:
            best = (score, path, data["capture"])

    if best is None:
        # Nothing readable -- fall back to time-based rather than failing.
        return pick_background(frames, image_dir, target_hour, target_minute,
                               camera=camera, select="time")
    return best[1], best[2], best[0]


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
    parser.add_argument("--frames-list", default=None,
                        help="File with one source_file name per line: draw only those frames "
                             "(e.g. one camera setup's, as lines from a differently aimed setup "
                             "land in the wrong place on this setup's background photo).")
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
    parser.add_argument("--max-lines", type=int, default=0,
                        help="Cap how many waterlines are drawn, subsampling EVENLY IN TIME "
                             "across the range (0 = draw all, the default). A 7-day map draws "
                             "about 80 lines and stays readable; a 30-day map would draw 500+ "
                             "and saturate the intertidal zone into a solid mass with the "
                             "photo invisible behind it. Subsampling in time rather than "
                             "taking the first N keeps the tidal range and the whole period "
                             "represented.")
    parser.add_argument("--background-select", choices=["clearest", "time"],
                        default="clearest",
                        help="How to choose the backdrop photo. 'clearest' (default) picks the "
                             "highest-contrast frame within the detector's crop region; 'time' "
                             "picks the one nearest --background-hour/minute. Time-based "
                             "selection repeatedly chose fogged frames, since midday on a "
                             "foggy day is still fog.")
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

    only = None
    if args.frames_list:
        with open(args.frames_list) as fh:
            only = {line.strip() for line in fh if line.strip()}
    frames, elev_col = load_contours(args.contour_csv, args.camera, only=only)
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

    if args.max_lines and len(frames) > args.max_lines:
        # Even in TIME, not every Nth file: captures are irregular
        # (nights and bad-weather frames are missing), so index-based
        # thinning would over-sample dense periods and under-sample
        # sparse ones.
        ordered = sorted(frames.items(), key=lambda kv: kv[1]["capture"])
        pick = np.linspace(0, len(ordered) - 1, args.max_lines).round().astype(int)
        pick = sorted(set(pick.tolist()))
        kept = {ordered[i][0]: ordered[i][1] for i in pick}
        print(f"Subsampled        : {len(frames)} waterlines -> {len(kept)} "
              f"(--max-lines {args.max_lines}), evenly spaced in time")
        frames = kept

    dates_used = sorted({d["capture"][:10] for d in frames.values()})
    date_label = (dates_used[0] if len(dates_used) == 1
                  else f"{dates_used[0]} to {dates_used[-1]}  ({len(dates_used)} days)")

    # The backdrop comes from ONE day -- by default the most recent in
    # range, so accumulated lines are drawn over the latest view of the
    # beach rather than a week-old one.
    bg_day = args.background_date or dates_used[-1]
    bg_candidates = {k: v for k, v in frames.items()
                     if v["capture"].startswith(bg_day)} or frames
    bg_path, bg_capture, bg_score = pick_background(
        bg_candidates, args.image_dir, args.background_hour, args.background_minute,
        camera=args.camera, select=args.background_select)
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
    n_setup = sum(bool(d["setup"]) for d in frames.values())
    if n_setup == len(frames):
        elev_what = "water level + wave setup"
    elif n_setup:
        elev_what = (f"water level + wave setup ({len(frames) - n_setup} of {len(frames)} "
                     f"lines: water level only, no wave record)")
    else:
        elev_what = "water level"
    norm = Normalize(vmin=elevations.min(), vmax=elevations.max())
    try:
        colormap = matplotlib.colormaps[args.colormap]      # matplotlib >= 3.5
    except AttributeError:
        colormap = plt.get_cmap(args.colormap)              # 3.3 / 3.4 (the station's floor)

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
    rows_per_m = rows_per_metre(np.vstack([resampled[k] for k in resampled]),
                                np.array([frames[k]["elevation"] for k in resampled])) \
        if resampled else None

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
    if n_lines > 250:
        # Beyond roughly 250 lines the band goes opaque whatever the
        # alpha, so this is a floor rather than a fix -- prefer
        # --max-lines to thin the set instead.
        line_alpha = min(line_alpha, 0.18)
        fill_alpha = min(fill_alpha, 0.05)
        line_width = min(line_width, 0.5)
    elif n_lines > 100:
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

        if len(keys) >= 3:
            stack = np.vstack([resampled[k] for k in keys])
            lo, hi, valid = agreeing_band(stack, rows_per_m=rows_per_m)
            if valid.any():
                # where= (not indexing) so a gap stays a gap instead of being
                # bridged by a straight edge across columns with no data.
                ax.fill_between(grid, lo, hi, where=valid,
                                color=colour, alpha=fill_alpha, linewidth=0)
                filled_bins += 1

        for k in keys:
            # NaNs break the line: plotting only the valid columns joined the
            # two sides of every gap with a straight segment, which looked like
            # a detected (flat) waterline where there was none.
            ax.plot(grid, resampled[k], color=colour,
                    alpha=line_alpha, linewidth=line_width)

    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.axis("off")

    ax.set_title(
        f"{args.camera.upper()}  {date_label}   {len(frames)} waterlines, "
        f"{elevations.min():+.2f} to {elevations.max():+.2f} m NAVD88 "
        f"({'water level + wave setup' if n_setup else 'water level'})\n"
        f"background: {bg_day} {bg_capture[11:16]} UTC   |   "
        + ("shaded bands = 16-84% spread of the 3+ same-elevation lines that agree "
           "(repeatability, not morphology)"
           if len(dates_used) == 1 else
           f"shaded bands = 16-84% spread of the 3+ lines that agree, per elevation, over "
           f"{len(dates_used)} days (repeatability + real shoreline movement)"),
        fontsize=9)

    scalar_map = matplotlib.cm.ScalarMappable(cmap=colormap, norm=norm)
    scalar_map.set_array([])
    cbar = fig.colorbar(scalar_map, ax=ax, fraction=0.030, pad=0.015)
    cbar.set_label(f"{'water level + wave setup' if n_setup else 'water level'} (m NAVD88)")

    fig.text(0.01, 0.005, f"lines: {Path(args.contour_csv).name}", fontsize=7, color="0.35",
             ha="left", va="bottom")
    fig.tight_layout()
    fig.savefig(args.output_png, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    print(f"Date range        : {date_label}")
    print(f"Camera            : {args.camera}")
    print(f"Waterlines drawn  : {len(resampled)}")
    print(f"Elevation range   : {elevations.min():+.3f} to {elevations.max():+.3f} m NAVD88 "
          f"({elev_what})")
    print(f"Elevation bins    : {len(bins)}  ({filled_bins} had 3+ agreeing lines somewhere and "
          f"were filled)")
    if bg_score is not None:
        print(f"Background image  : {bg_path.name}")
        print(f"                    {bg_capture} UTC, in-crop contrast {bg_score:.1f} "
              f"(picked as clearest of {len(bg_candidates)} candidate(s))")
    else:
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
        print("Reminder: filled band = 16-84% spread of the same-elevation lines that agree,")
        print("within ONE day. It measures repeatability (detection error, runup differences,")
        print("or change during the day) -- not a morphology feature. Stray lines are drawn")
        print("but do not widen it.")
    else:
        print(f"Reminder: this spans {len(dates_used)} days, so a filled band now mixes TWO")
        print("things: measurement repeatability, AND genuine movement of the shoreline at")
        print("that elevation over the period. A band that widens across the week is the")
        print("interesting case, but this figure alone cannot separate real change from")
        print("detection scatter -- compare against the single-day bands to judge which.")

    if filled_bins == 0:
        print()
        print("NOTE: no elevation bin had three or more agreeing waterlines, so nothing was")
        print("filled. That happens when the day's captures never revisit the same level")
        print(f"within +/-{args.elevation_bin} m -- try a wider --elevation-bin, or a day")
        print("whose captures span both a rising and a falling tide.")


def band_self_test():
    """
    python3 daily_elevation_map.py --self-test

    The band's shape on synthetic stacks of one bin's lines (rows; 600
    columns), each case one way the station's maps went wrong:
      1. a stray 115 px from two shore lines 31 px apart (0.5 m = 30 px
         there): no band where only those two agree (the c2 'bottle');
      2. two strays 1 px apart, 100 px off eight shore lines: the band
         stays on the shore lines;
      3. real change, four lines before and four 35 px further after
         (0.35 m): both groups stay, the band spans both;
      4. a ~240 px tall band (lines spread evenly): its edges move at
         most 2 px per column, so it is at 80% of its height no sooner
         than 0.8 x 240 / 4 = 48 columns from its ends -- no wall;
      5. three columns where only two of eight lines reach: bridged, no
         notch.
    Returns 0 if all pass, 1 otherwise.
    """
    rng = np.random.default_rng(1)
    W = 600
    cols = np.arange(W)
    fails = []

    def lines(rows, noise=1.0):
        return np.array([np.full(W, r) + rng.normal(0, noise, W) for r in rows])

    # 1
    st = lines([686, 717, 571])
    lo, hi, ok = agreeing_band(st, rows_per_m=np.full(W, 60.0))
    print(f"  1. stray + two shore lines: band in {ok.sum()} of {W} columns")
    if ok.any():
        fails.append("1: a band where only two lines agree")
    # 2
    st = lines([700, 702, 705, 707, 709, 711, 714, 716, 604, 605])
    lo, hi, ok = agreeing_band(st, rows_per_m=np.full(W, 60.0))
    top = float(np.nanmin(np.where(ok, lo, np.nan)))
    print(f"  2. two strays 1 px apart: band top row {top:.0f} (shore lines from 700)")
    if not ok.any() or top < 690:
        fails.append(f"2: band reaches the strays (top {top:.0f})")
    # 3
    st = lines([700, 703, 705, 708, 735, 738, 740, 743])
    lo, hi, ok = agreeing_band(st, rows_per_m=np.full(W, 100.0))
    mid = slice(200, 400)
    print(f"  3. real change (0.35 m): band {np.nanmean(lo[mid]):.0f}-{np.nanmean(hi[mid]):.0f}")
    if not ok[mid].all() or np.nanmean(lo[mid]) > 706 or np.nanmean(hi[mid]) < 737:
        fails.append("3: the band does not span both groups")
    # 4
    st = lines(np.linspace(560, 920, 9), noise=0.5)
    st[:, :50] = np.nan
    st[:, 550:] = np.nan
    lo, hi, ok = agreeing_band(st)
    h = np.where(ok, hi - lo, 0.0)
    H = h.max()
    run = np.flatnonzero(ok)
    to80 = min(int(np.argmax(h[run] >= 0.8 * H)), int(np.argmax(h[run][::-1] >= 0.8 * H)))
    slope = max(np.abs(np.diff(lo[run])).max(), np.abs(np.diff(hi[run])).max())
    print(f"  4. {H:.0f} px band: 80% of its height {to80} columns in from its ends, "
          f"edges move <= {slope:.1f} px/col")
    if to80 < 0.8 * H / 4.0 - 2 or slope > 2.01:
        fails.append(f"4: wall at the end (80% in {to80} columns, {slope:.1f} px/col)")
    # 5
    st = lines([700, 704, 707, 710, 713, 716, 720, 723])
    st[2:, 300:303] = np.nan
    lo, hi, ok = agreeing_band(st)
    h = np.where(ok, hi - lo, 0.0)
    dip = h[298:305].min() / np.median(h[100:250])
    print(f"  5. three columns with two lines: band {'bridged' if ok[300:303].all() else 'BROKEN'}, "
          f"narrowest {dip:.0%} of its usual height")
    if not ok[300:303].all() or dip < 0.6:
        fails.append("5: notch where two of eight lines remain")
    print("SELF-TEST " + ("PASSED" if not fails else "FAILED: " + "; ".join(fails)))
    return 0 if not fails else 1


if __name__ == "__main__":
    if "--self-test" in sys.argv[1:]:
        sys.exit(band_self_test())
    main()

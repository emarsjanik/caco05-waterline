#!/usr/bin/env python3
"""
Extract Elevation Contours
-----------------------------
Matches each detected-shoreline CSV (from waterline_detector_v5.py's
processed/ output) to a water elevation from a tide model, using the
timestamp already embedded in the filename, and produces one combined
"contour points" file: every (column, row, tide_elevation) triple
across every processed image.

WHY THIS SHAPE: this is the core idea behind video-based intertidal
beach mapping (Plant & Holman 1997; Aarninkhof et al.) -- each
detected shoreline is a snapshot of where the water's edge sits at one
instant, so if you know the water elevation at that instant, the
ENTIRE shoreline is a contour line at that elevation. Collect shorelines
across many different tide levels and you get a set of elevation-
tagged contour lines that can be gridded into a DEM.

WHAT THIS SCRIPT DOES NOT DO (yet): output columns are still in PIXEL
space (column, row within the cropped/full image), not real-world
ground coordinates (x, y in meters). Converting pixel -> ground
requires camera georectification (intrinsic + extrinsic calibration),
which this project does not have yet (see project history: CACO-05's
calibration status is still being confirmed with USGS Woods Hole).
Once that exists, insert the pixel->ground mapping as a step between
this script's output and the DEM-gridding step -- everything here
(tide matching, contour point structure) stays the same either way.

TIDE MODEL FORMAT: tide model file formats vary a lot (NOAA CO-OPS
exports, custom prediction tables, etc.), so this script tries to
auto-detect a timestamp column and a water-level column by common
names. If auto-detection fails, pass --time-col and --level-col
explicitly to name them.

Usage:
    python3 extract_elevation_contours.py <tide_model_csv> [options]

Options:
    --processed-dir DIR   default: /mnt/I2Rgus_Data/waterline/processed
    --output FILE         default: /mnt/I2Rgus_Data/waterline/contour_points.csv
    --time-col NAME       explicit timestamp column name in tide model CSV
    --level-col NAME      explicit water-level column name in tide model CSV
    --max-gap-minutes N   warn (and by default skip) if the nearest tide
                          reading is more than N minutes from the image
                          timestamp (default: 60)
    --include-stale       include matches beyond --max-gap-minutes anyway
                          (flagged in output) instead of skipping them
"""

import argparse
import csv
import re
import sys
from pathlib import Path
from datetime import datetime, timezone

from collections import defaultdict
import numpy as np
import pandas as pd


DEFAULT_PROCESSED_DIR = "/mnt/I2Rgus_Data/waterline/processed"
DEFAULT_OUTPUT = "/mnt/I2Rgus_Data/waterline/contour_points.csv"

# Common column name variants seen across tide model / NOAA CO-OPS
# exports and simple custom prediction tables.
TIME_COL_CANDIDATES = [
    "date time", "datetime", "date_time", "time", "timestamp", "date", "t",
]
LEVEL_COL_CANDIDATES = [
    "prediction", "predicted", "predicted (ft)", "predicted (m)",
    "water level", "water_level", "wl", "verified (ft)", "verified (m)",
    "level", "elevation", "value",
]
# Global ocean tide models (pyTMD-style output) commonly appear as
# "<MODEL>_heightm" columns, e.g. EOT20_heightm, GOT5.5_heightm,
# FES2022_heightm -- often several at once as a cross-check between
# models, rather than one single "the" water level column.
HEIGHT_COL_PATTERN = re.compile(r"_height\s*m?$", re.IGNORECASE)


def detect_column(header, candidates):
    header_lower = [h.strip().lower() for h in header]
    for cand in candidates:
        for i, h in enumerate(header_lower):
            if h == cand or cand in h:
                return header[i]
    return None


def detect_height_columns(header):
    return [h for h in header if HEIGHT_COL_PATTERN.search(h.strip())]


def load_tide_model(path, time_col=None, level_col=None, tide_timezone="UTC"):
    """
    Loads a tide model from .xlsx/.xls (via pandas) or .csv, and
    returns (timestamps_epoch_utc, levels, spreads). `spreads` is the
    std-dev across multiple model columns at each timestamp when more
    than one height column is present and none was explicitly chosen
    via level_col (0.0 everywhere otherwise) -- a rough per-point
    cross-model agreement metric, carried through to the output as a
    QA column rather than silently discarded.

    TIMEZONE: tide model timestamps are assumed to be `tide_timezone`
    (default UTC, matching the "GMT" already in this project's image
    filenames) if the timestamps are timezone-naive. This is an
    explicit assumption, not a silent guess -- pass --tide-timezone if
    your tide model uses local time instead.
    """
    suffix = Path(path).suffix.lower()
    if suffix in (".xlsx", ".xls"):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)

    header = list(df.columns)

    if time_col is None:
        time_col = detect_column(header, TIME_COL_CANDIDATES)
    if time_col is None:
        raise ValueError(
            f"Could not auto-detect a timestamp column in tide model "
            f"header: {header}. Pass --time-col explicitly."
        )

    level_cols = []
    if level_col is not None:
        if level_col not in header:
            raise ValueError(f"--level-col '{level_col}' not found in header: {header}")
        level_cols = [level_col]
    else:
        height_cols = detect_height_columns(header)
        if height_cols:
            level_cols = height_cols
        else:
            single = detect_column(header, LEVEL_COL_CANDIDATES)
            if single is None:
                raise ValueError(
                    f"Could not auto-detect a water-level column in tide model "
                    f"header: {header}. Pass --level-col explicitly."
                )
            level_cols = [single]

    print(f"Time column   : '{time_col}'")
    if len(level_cols) > 1:
        print(f"Level columns : {level_cols} (averaging {len(level_cols)} models; "
              f"per-point std-dev across models reported as tide_model_spread)")
    else:
        print(f"Level column  : '{level_cols[0]}'")

    timestamps = pd.to_datetime(df[time_col], errors="coerce")
    valid = timestamps.notna()
    for col in level_cols:
        valid &= df[col].notna()
    df = df[valid].copy()
    timestamps = timestamps[valid]

    if timestamps.dt.tz is None:
        timestamps = timestamps.dt.tz_localize(tide_timezone)
        print(f"Tide timestamps were timezone-naive -- localized as '{tide_timezone}'. "
              f"Pass --tide-timezone if this is wrong for your tide model.")
    timestamps = timestamps.dt.tz_convert("UTC")

    levels_matrix = df[level_cols].to_numpy(dtype=float)
    levels = levels_matrix.mean(axis=1)
    spreads = levels_matrix.std(axis=1) if len(level_cols) > 1 else np.zeros(len(levels))

    epoch = np.array([t.timestamp() for t in timestamps])

    order = np.argsort(epoch)
    epoch = epoch[order]
    levels = levels[order]
    spreads = spreads[order]

    print(f"Loaded {len(epoch)} tide readings, spanning "
          f"{datetime.fromtimestamp(epoch[0], tz=timezone.utc)} to "
          f"{datetime.fromtimestamp(epoch[-1], tz=timezone.utc)}")
    if len(level_cols) > 1:
        print(f"Cross-model spread: mean={spreads.mean():.4f}, max={spreads.max():.4f} "
              f"(same units as tide levels)")
    return epoch, levels, spreads


GNSSR_GAP_SENTINEL = 999.0


def looks_like_gnssr_spline(path):
    """
    gnssrefl subdaily spline output is a plain text file whose header
    lines begin with '%'. Detect it that way rather than by extension,
    since tide models also arrive as .csv/.txt.
    """
    try:
        with open(path, "r") as f:
            for _ in range(5):
                line = f.readline()
                if not line:
                    break
                if line.lstrip().startswith("%"):
                    return True
    except (OSError, UnicodeDecodeError):
        return False
    return False


def load_wave_archive(path):
    """
    Reads the offshore wave archive written by fetch_buoy_waves.py.
    Returns (epochs_hs, hs, epochs_tp, tp), each sorted by time, with
    rows lacking a value dropped separately for height and period.
    """
    df = pd.read_csv(path)
    hs = df.dropna(subset=["wvht_m"]).sort_values("epoch")
    tp = df.dropna(subset=["dpd_s"]).sort_values("epoch")
    return (hs["epoch"].to_numpy(float), hs["wvht_m"].to_numpy(float),
            tp["epoch"].to_numpy(float), tp["dpd_s"].to_numpy(float))


def nearest_within(epochs, values, epoch, max_gap_s):
    """Value of the record nearest `epoch`, or None if none within max_gap_s."""
    if len(epochs) == 0:
        return None
    i = int(np.searchsorted(epochs, epoch))
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(epochs):
            gap = abs(epochs[j] - epoch)
            if gap <= max_gap_s and (best is None or gap < best[0]):
                best = (gap, values[j])
    return None if best is None else float(best[1])


def load_gnssr_spline(path):
    """
    Loads gnssrefl subdaily spline output as a water-level source.

    Columns (1-indexed): 1 MJD, 2 RH(m), 3 YYYY, 4 MM, 5 DD, 6 HH,
    7 MM, 8 SS, 9 quasi-sea-level(m). Column 9 is computed upstream as
    (station orthometric height - reflector height), and orthometric
    height in the US means NAVD88 -- the SAME datum as the camera EO
    z-coordinates. So GNSS-R needs NO vertical datum offset at all,
    unlike the global tide models.

    It is also a MEASUREMENT at this station, so it includes storm
    surge, wind setup and local bias that an astronomical tide model
    cannot predict. Against ~47 days of overlap the model differed
    from GNSS-R by about -0.25 m in the mean with a 0.08 m standard
    deviation; using GNSS-R removes both of those terms.

    Returns (epochs, levels, spreads, hortho). `spreads` is zeros:
    there is no cross-model spread for a single measured series, and
    fabricating one would misrepresent it. Note the spline is a
    SMOOTHED fit, not raw observations (the file header says so), so
    it will flatten genuinely short-period signal.
    """
    times, levels = [], []
    hortho = None

    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("%"):
                if "orthometric height" in s.lower() or "Hortho" in s:
                    for tok in s.replace(",", " ").split():
                        try:
                            val = float(tok)
                            if 0.0 < val < 1000.0:
                                hortho = val
                        except ValueError:
                            continue
                continue
            parts = s.split()
            if len(parts) < 9:
                continue
            try:
                dt = datetime(int(parts[2]), int(parts[3]), int(parts[4]),
                              int(parts[5]), int(parts[6]), int(parts[7]),
                              tzinfo=timezone.utc)
                level = float(parts[8])
            except ValueError:
                continue
            if abs(level - GNSSR_GAP_SENTINEL) < 1e-6:
                continue
            times.append(dt)
            levels.append(level)

    if not times:
        raise ValueError(f"No usable GNSS-R rows parsed from {path}")

    epoch = np.array([t.timestamp() for t in times])
    order = np.argsort(epoch)
    epoch = epoch[order]
    levels = np.array(levels)[order]

    print(f"Water-level source: GNSS-R MEASURED (gnssrefl spline)")
    if hortho is not None:
        print(f"  station orthometric height (NAVD88): {hortho:.3f} m")
        print(f"  water elevation = {hortho:.3f} - reflector height  (already NAVD88)")
    print(f"  Loaded {len(epoch)} readings, spanning "
          f"{datetime.fromtimestamp(epoch[0], tz=timezone.utc)} to "
          f"{datetime.fromtimestamp(epoch[-1], tz=timezone.utc)}")
    print("  No vertical datum offset required -- this source is already NAVD88.")
    return epoch, levels, np.zeros(len(epoch)), hortho


def interpolate_tide(epoch, tide_timestamps, tide_levels, tide_spreads):
    """
    Linear interpolation between the two nearest tide readings.
    Returns (elevation, spread, gap_minutes). gap_minutes is the distance
    from `epoch` to the FARTHER of the two readings it is interpolated
    between -- i.e. how far the straight line is stretched. An earlier
    version used the NEARER reading, so an image 10 min after the last
    reading before a 6-hour GNSS-R outage passed a 60-min check while
    its level was a straight line drawn across the whole outage.
    Outside the record the end value is held, and the gap is the
    distance to that end.
    """
    if epoch <= tide_timestamps[0]:
        gap_minutes = abs(epoch - tide_timestamps[0]) / 60.0
        return float(tide_levels[0]), float(tide_spreads[0]), gap_minutes
    if epoch >= tide_timestamps[-1]:
        gap_minutes = abs(epoch - tide_timestamps[-1]) / 60.0
        return float(tide_levels[-1]), float(tide_spreads[-1]), gap_minutes

    idx = np.searchsorted(tide_timestamps, epoch)
    t0, t1 = tide_timestamps[idx - 1], tide_timestamps[idx]
    l0, l1 = tide_levels[idx - 1], tide_levels[idx]
    s0, s1 = tide_spreads[idx - 1], tide_spreads[idx]
    frac = (epoch - t0) / (t1 - t0)
    elevation = l0 + frac * (l1 - l0)
    spread = s0 + frac * (s1 - s0)
    gap_minutes = max(abs(epoch - t0), abs(epoch - t1)) / 60.0
    return float(elevation), float(spread), gap_minutes


def extract_epoch_from_filename(name):
    """
    Filenames in this project start with a Unix epoch timestamp, e.g.
    "1785148200.Mon.Jul.27_10_30_00.GMT.2026.CACO05.c1.snap.csv".
    Returns None (rather than guessing) if no leading epoch is found.
    """
    match = re.match(r"^(\d{9,10})\.", name)
    if not match:
        return None
    return int(match.group(1))


def detect_camera(name):
    name = name.lower()
    if ".c1." in name:
        return "c1"
    if ".c2." in name:
        return "c2"
    return None


# Fewest pixel columns two consecutive frames must share for their
# tide-direction comparison to count.
MIN_COMMON_COLUMNS = 50

# The station's local time zone. "Day" filters are judged on local days:
# UTC dates split each evening's captures (20:00 EDT = 00:00 UTC) into
# the next day.
LOCAL_TZ = "America/New_York"


def local_date(epoch):
    """YYYY-MM-DD of `epoch` in the station's local time zone."""
    return pd.Timestamp(int(epoch), unit="s", tz="UTC").tz_convert(LOCAL_TZ).strftime("%Y-%m-%d")


def tide_agreement_by_day(output_rows):
    """
    Measures, per camera-day, how often the detected waterline moves in
    the SAME DIRECTION as the tide between consecutive frames.

    A real waterline has no choice about this: if the water rises, the
    line must move one way; if it falls, the other. The magnitude can
    vary with beach slope, but the SIGN cannot. So this tests the one
    thing that must hold, without needing to judge whether an image
    looks degraded.

    That distinction matters here. Four image-quality statistics --
    brightness, whole-frame contrast, in-band contrast, and detector
    confidence -- all failed to separate a day of torrential rain and
    fog from clear days: the bad day scored 0.950-0.998 confidence and
    its contrast overlapped the good days. Measured agreement on the
    same three days was 81%, 89% and 56%, ranking them in the same
    order as their correlation with tide (+0.97, +0.73, -0.17). A day
    at 56% is a coin flip: those detections carry no information about
    the water level at all.

    Returns {(camera, date): (agreeing, total)}. Transitions are only
    counted between frames less than ~90 min apart and where the tide
    actually moved more than 0.05 m, since a near-stationary tide gives
    the sign no meaning.

    Consecutive frames are compared over the pixel columns BOTH kept.
    Averaging each frame over whichever columns survived its own
    filtering (as an earlier version did) compares different stretches
    of beach whenever far-field coverage varies, which alone can flip
    the sign of a transition. Days are LOCAL days (see local_date).
    """
    per_frame = {}
    for row in output_rows:
        stem, camera, epoch = row[0], row[1], row[3]
        key = (camera, stem)
        if key not in per_frame:
            per_frame[key] = {"epoch": epoch, "date": local_date(epoch),
                              "elev": row[6], "rows": {}}
        per_frame[key]["rows"][row[4]] = row[5]

    by_camera = defaultdict(list)
    for (camera, _stem), v in per_frame.items():
        by_camera[camera].append((v["epoch"], v["date"], v["elev"], v["rows"]))

    agreement = defaultdict(lambda: [0, 0])
    for camera, entries in by_camera.items():
        entries.sort(key=lambda e: e[0])
        previous = None
        for epoch, date, elev, rows in entries:
            if (previous is not None
                    and epoch - previous[0] <= 5400
                    and abs(elev - previous[1]) > 0.05):
                common = rows.keys() & previous[2].keys()
                if len(common) >= MIN_COMMON_COLUMNS:
                    now = float(np.mean([rows[c] for c in common]))
                    before = float(np.mean([previous[2][c] for c in common]))
                    agreement[(camera, date)][1] += 1
                    if np.sign(elev - previous[1]) == np.sign(now - before):
                        agreement[(camera, date)][0] += 1
            previous = (epoch, elev, rows)
    return agreement


def straightness_mask(rows, window, min_residual):
    """
    Flags columns where the detected line is implausibly straight.

    A dynamic-programming solver with no signal to follow falls back on
    its smoothness penalty and emits a dead-straight run. A real
    waterline never does: even after sub-pixel refinement it retains
    small column-to-column wiggle from swash, foam and beach texture.

    Measured on 40 archived C2 frames with a 31-column window, the RMS
    deviation from a local straight-line fit was:
        columns the detector kept     median 0.173 px
        columns it flagged no-signal  median 0.004 px
    a 43x separation, with 83% of flagged columns below 0.05 px against
    only 6.9% of kept ones.

    This is INDEPENDENT of the Has_Signal test, which is why it is
    worth having as well: Has_Signal asks whether the energy was strong
    enough, this asks whether the answer looks like a waterline. Some
    straight segments pass the energy test and are caught only here.

    Returns True for columns to KEEP.
    """
    n = len(rows)
    if window < 5 or n < window or min_residual <= 0:
        return np.ones(n, dtype=bool)

    half = window // 2
    x = np.arange(window, dtype=float)
    A = np.vstack([x, np.ones(window)]).T
    # Least-squares projection matrix, computed once: the residual for
    # every window is then a single matrix multiply rather than a fit.
    pinv = np.linalg.pinv(A)

    keep = np.ones(n, dtype=bool)
    for i in range(half, n - half):
        seg = rows[i - half:i - half + window]
        coef = pinv @ seg
        residual = np.sqrt(np.mean((seg - (A @ coef)) ** 2))
        if residual < min_residual:
            keep[i] = False

    # Window edges are never evaluated, so inherit from the nearest
    # evaluated column rather than being kept by default -- otherwise a
    # straight run reaching the frame edge keeps a fringe of survivors.
    if n > 2 * half:
        keep[:half] = keep[half]
        keep[n - half:] = keep[n - half - 1]
    return keep


def local_coverage_mask(has_signal, window, min_fraction):
    """
    Marks columns that sit in a NEIGHBOURHOOD with adequate signal,
    rather than judging each column alone.

    Why this is needed: a column can pass the detector's own signal
    test while being an isolated survivor in a region that is
    otherwise dead. Measured on real C2 data, low-tide frames retained
    2.7-12% of their far-field columns -- a scatter of fragments that
    the per-column filter kept, and whose mean row then broke the
    physical elevation ordering (-1.4 m contours plotting ABOVE -1.0 m
    ones). The near field of those same frames was perfectly good, so
    discarding the whole frame would throw away the lowest elevations
    in the record, which are the hardest part of the tidal range to
    sample.

    Truncating instead keeps the near field and cuts the fragmentary
    far field. Operating on a sliding window makes this geometry-
    agnostic: it trims whichever end (or interior region) has gone
    sparse, without assuming which side of the frame is far from the
    camera.
    """
    n = len(has_signal)
    k = int(window)
    if k < 3 or n == 0:
        return np.ones(n, dtype=bool)
    kernel = np.ones(k, dtype=float)
    counts = np.convolve(has_signal.astype(float), kernel, mode="same")
    denom = np.convolve(np.ones(n, dtype=float), kernel, mode="same")
    fraction = counts / np.maximum(denom, 1.0)
    return fraction >= min_fraction


def load_shoreline_csv(path):
    """
    Returns (columns, rows, peak_sharpness). Prefers the sub-pixel
    "Row_Precise" column when present (detector v6.1+) over the
    rounded-integer "Row" column, since fractional-pixel precision
    matters once this feeds into an elevation estimate -- rounding to
    whole pixels here would throw away exactly the precision the
    sub-pixel refinement was added to provide. Falls back to "Row"
    (as int) for older detector output that predates this column, and
    "Peak_Sharpness" defaults to NaN if absent, rather than a fabricated
    value.
    """
    columns, rows, sharpness, has_signal = [], [], [], []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        has_precise = reader.fieldnames and "Row_Precise" in reader.fieldnames
        has_sharpness = reader.fieldnames and "Peak_Sharpness" in reader.fieldnames
        has_signal_col = reader.fieldnames and "Has_Signal" in reader.fieldnames
        for line in reader:
            columns.append(int(line["Column"]))
            if has_precise and line["Row_Precise"] not in ("", None):
                rows.append(float(line["Row_Precise"]))
            else:
                rows.append(float(line["Row"]))
            if has_sharpness and line["Peak_Sharpness"] not in ("", None):
                sharpness.append(float(line["Peak_Sharpness"]))
            else:
                sharpness.append(float("nan"))
            # Has_Signal=0 marks a column where the detector found no
            # usable energy, so the row there was chosen without
            # evidence. Absent in pre-v6.7 output: assume 1 rather
            # than 0, since discarding all older data would be worse
            # than carrying it forward unflagged.
            if has_signal_col and line["Has_Signal"] not in ("", None):
                has_signal.append(int(float(line["Has_Signal"])))
            else:
                has_signal.append(1)
    return (np.array(columns), np.array(rows), np.array(sharpness),
            np.array(has_signal, dtype=bool))


def main():
    parser = argparse.ArgumentParser(description="Match detected shorelines to tide elevations, "
                                                   "producing elevation-tagged contour points.")
    parser.add_argument("tide_model_csv", metavar="WATER_LEVEL_FILE",
                        help="Either a tide model (.xlsx/.csv with '<model>_heightm' columns) "
                             "or gnssrefl subdaily spline output (usgs_spline_out.txt). The "
                             "type is auto-detected. GNSS-R is preferred where it covers the "
                             "imagery: it is measured rather than predicted, already in "
                             "NAVD88, and includes surge and setup.")
    parser.add_argument("--processed-dir", default=DEFAULT_PROCESSED_DIR)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--time-col", default=None)
    parser.add_argument("--level-col", default=None,
                         help="Explicit water-level column name. If omitted and multiple "
                              "'<model>_heightm' columns are found (e.g. several global tide "
                              "models), they are averaged and their spread reported as a QA metric.")
    parser.add_argument("--tide-timezone", default="UTC",
                         help="Timezone to assume for timezone-naive tide model timestamps "
                              "(default: UTC, matching this project's GMT-labeled image filenames).")
    parser.add_argument("--max-gap-minutes", type=float, default=60.0)
    parser.add_argument("--include-stale", action="store_true")
    parser.add_argument("--min-tide-agreement", type=float, default=0.70,
                        help="Reject a whole camera-day if the detected waterline moves in "
                             "the same direction as the tide less often than this (default "
                             "0.70). Measured on three days: 81%%, 89%% and 56%%, ranking them "
                             "as their tide correlation did (+0.97, +0.73, -0.17). The 56%% "
                             "day was torrential rain and fog, which brightness, contrast, "
                             "in-band contrast and confidence all failed to detect. Set 0 to "
                             "disable. NOTE: this judges a whole DAY, so it discards usable "
                             "frames within a mixed day -- appropriate for day-long weather, "
                             "less so for an afternoon that fogs in.")
    parser.add_argument("--min-agreement-transitions", type=int, default=5,
                        help="Minimum consecutive-frame transitions before the agreement test "
                             "is applied to a day (default 5). Below this the fraction is too "
                             "noisy to act on, so the day is kept and reported rather than "
                             "judged.")
    parser.add_argument("--straightness-window", type=int, default=31,
                        help="Columns used to judge local straightness (default 31). Measured "
                             "on 40 archived frames this window separated known-bad from "
                             "known-good columns by 43x in median residual -- better than 61 "
                             "or 121.")
    parser.add_argument("--straightness-residual-c1", type=float, default=0.0,
                        help="Straightness threshold for c1, in pixels. DISABLED by default "
                             "(0.0) on measurement: c1's flagged and kept columns separate by "
                             "only 1.8x in median residual (0.083 vs 0.150), so a threshold "
                             "catches just 34%% of known-bad while removing 13.2%% of "
                             "known-good. c1's failures are fog days where the detector still "
                             "found something with natural-looking wiggle -- wrong, but not "
                             "straight. The tide-agreement filter is what catches those.")
    parser.add_argument("--straightness-residual-c2", type=float, default=0.05,
                        help="Straightness threshold for c2, in pixels (default 0.05). c2's "
                             "populations separate by 43x (0.004 vs 0.173), because its "
                             "failure mode is the far-field solver flat-lining, which is dead "
                             "straight. Catches 83%% of known-bad for 6.9%% of known-good.")
    parser.add_argument("--min-straightness-residual", type=float, default=None,
                        help="Drop columns whose deviation from a local straight-line fit is "
                             "below this, in pixels (default 0.05). Catches 83%% of columns the "
                             "detector already flags as signal-free, and some it does not -- a "
                             "DP with nothing to follow emits dead-straight runs. Costs about "
                             "6.9%% of good columns, which are the flattest and least "
                             "informative. Overrides BOTH per-camera values above when given "
                             "-- use it to sweep a single value, not in production, since the "
                             "cameras demonstrably need different thresholds.")
    parser.add_argument("--max-straight-fraction", type=float, default=0.50,
                        help="If more than this fraction of a frame's columns are straight, "
                             "drop the frame entirely rather than trimming (default 0.50). A "
                             "mostly-straight line is not a waterline with a bad patch.")
    parser.add_argument("--truncate-window", type=int, default=120,
                        help="Width in columns of the sliding window used to judge local "
                             "signal density (default 120). Columns whose neighbourhood is "
                             "sparser than --truncate-min-coverage are dropped even if they "
                             "individually passed the detector's signal test.")
    parser.add_argument("--truncate-min-coverage", type=float, default=0.60,
                        help="Minimum local signal fraction for a column to be kept (default "
                             "0.60). On real C2 data, far-field coverage was 92-99.9%% on good "
                             "frames and 2.7-47%% on bad ones, so anything in 0.5-0.8 separates "
                             "them; 0.60 sits mid-gap. Set 0 to disable truncation.")
    parser.add_argument("--min-coverage", type=float, default=0.50,
                        help="Minimum fraction of a frame's columns that must survive the "
                             "no-signal filter for the frame to be used at all (default 0.50). "
                             "A frame reduced to a handful of scattered columns is not a "
                             "waterline: measured on real C2 data, frames retaining 13-54 of "
                             "447 far-field columns were the ones breaking the physical "
                             "elevation ordering, because a mean over unrepresentative "
                             "fragments is meaningless. Set 0 to disable.")
    parser.add_argument("--include-no-signal", action="store_true",
                        help="Keep columns the detector flagged Has_Signal=0. By default they "
                             "are DROPPED: the detector found no usable energy there, so the "
                             "row was chosen without evidence and the elevation assigned to it "
                             "is fabricated. C2's far field runs 15-45%% such columns, so this "
                             "materially affects that camera. Pre-v6.7 output has no flag and "
                             "is treated as signal-bearing either way.")
    parser.add_argument("--waves", default=None,
                        help="Offshore wave archive from fetch_buoy_waves.py. Adds "
                             "offshore_hs_m and offshore_tp_s to every point, so rough-water "
                             "frames (whose waterline sits above still water by the wave "
                             "setup) can be left out of the DEM.")
    parser.add_argument("--wave-max-gap-minutes", type=float, default=90.0,
                        help="Largest gap to the nearest wave record before a frame's wave "
                             "columns are left blank (default 90; the buoy reports every "
                             "10-60 min).")
    parser.add_argument("--setup-coef", type=float, default=None,
                        help="Wave-setup correction coefficient C (needs --waves). Each "
                             "frame's beach elevation becomes water level + C*sqrt(Hs*L0), "
                             "L0 = g*Tp^2/(2*pi), from the offshore buoy record: the timex "
                             "waterline is the MEAN swash position, which sits above still "
                             "water by the setup. Stockdon et al. (2006) give C = 0.35*slope; "
                             "fit it for this site with dem_from_contours.py --fit-setup. "
                             "Adds setup_correction_m and beach_elevation_navd88 columns, "
                             "which georectify.py and dem_from_contours.py then use; "
                             "tide_elevation_navd88 stays the water level.")
    parser.add_argument("--vertical-datum-offset", type=float, default=0.0,
                         help="Meters to ADD to every tide elevation, converting the tide "
                              "model's vertical datum (mean sea level or a geoid, for global "
                              "ocean models like EOT20/GOT/FES) into NAVD88 -- the datum the "
                              "camera EO z-coordinates use. The preferred way to obtain this "
                              "is compare_gnssr_to_tidemodel.py, which measures it against "
                              "GNSS-R water levels that are already in NAVD88 (about -0.25 m "
                              "at this station, +/- 0.04 m depending on the comparison "
                              "window). NOAA's VDatum tool is an alternative but gives a "
                              "modelled rather than measured value. Defaults to 0.0, which "
                              "is almost certainly WRONG -- see the warning printed if unset.")
    args = parser.parse_args()

    using_gnssr = looks_like_gnssr_spline(args.tide_model_csv)

    if using_gnssr:
        tide_timestamps, tide_levels, tide_spreads, _hortho = load_gnssr_spline(
            args.tide_model_csv)
        source_label = "gnssr_measured"
        if args.vertical_datum_offset != 0.0:
            print()
            print("=" * 90)
            print(f"WARNING: --vertical-datum-offset {args.vertical_datum_offset} was supplied, but the")
            print("water-level source is GNSS-R, which is ALREADY in NAVD88. Applying an offset")
            print("here would shift correct elevations into being wrong. Ignoring it.")
            print("=" * 90)
            print()
            args.vertical_datum_offset = 0.0
    else:
        tide_timestamps, tide_levels, tide_spreads = load_tide_model(
            args.tide_model_csv, time_col=args.time_col, level_col=args.level_col,
            tide_timezone=args.tide_timezone,
        )
        source_label = "tide_model"

    if (not using_gnssr) and args.vertical_datum_offset == 0.0:
        print()
        print("=" * 90)
        print("WARNING: --vertical-datum-offset is 0.0 (the default).")
        print("=" * 90)
        print("The tide values in this model are almost certainly referenced to mean sea")
        print("level or a geoid (standard for global ocean tide models like EOT20/GOT/FES),")
        print("NOT to NAVD88 -- the datum this station's camera EO z-coordinates use. Using")
        print("tide_elevation as-is means every elevation this script produces carries a")
        print("constant, uncorrected offset of unknown size (commonly tens of centimeters).")
        print()
        print("To fix this properly, derive the offset from GNSS-R water levels, which are")
        print("already referenced to NAVD88:")
        print("    python3 compare_gnssr_to_tidemodel.py <gnssr_spline_file> <tide_file>")
        print("then re-run with:")
        print("    --vertical-datum-offset <value_in_meters>")
        print()
        print("At this station that came out near -0.25 m, though it varied by about 0.04 m")
        print("between comparison windows -- so treat it as approximate, not exact. Where")
        print("GNSS-R covers the imagery directly it is better still: no offset is needed.")
        print("=" * 90)
        print()

    processed_dir = Path(args.processed_dir)
    shoreline_files = sorted(processed_dir.glob("*.csv"))

    if not shoreline_files:
        print(f"No shoreline CSVs found in {processed_dir}")
        sys.exit(1)

    print(f"Found {len(shoreline_files)} processed shoreline file(s) in {processed_dir}")
    print()

    output_rows = []
    used, skipped_no_epoch, skipped_no_camera, skipped_stale = 0, 0, 0, 0
    dropped_no_signal = 0
    skipped_all_no_signal = 0
    skipped_low_coverage = 0
    low_coverage_detail = []
    truncated_columns = 0
    straight_trimmed = 0
    skipped_mostly_straight = 0

    for path in shoreline_files:
        epoch = extract_epoch_from_filename(path.name)
        if epoch is None:
            skipped_no_epoch += 1
            continue

        camera = detect_camera(path.name)
        if camera is None:
            skipped_no_camera += 1
            continue

        elevation, spread, gap_minutes = interpolate_tide(
            epoch, tide_timestamps, tide_levels, tide_spreads
        )
        elevation_navd88 = elevation + args.vertical_datum_offset

        is_stale = gap_minutes > args.max_gap_minutes
        if is_stale and not args.include_stale:
            skipped_stale += 1
            continue

        columns, rows, sharpness, col_has_signal = load_shoreline_csv(path)
        capture_time = datetime.fromtimestamp(epoch, tz=timezone.utc).isoformat()

        # Straightness is judged on the FULL row series before any
        # columns are removed: dropping columns first would leave gaps
        # that make a straight run look wiggly.
        rows_full_for_straightness = rows.copy()

        if not args.include_no_signal:
            keep = col_has_signal
            total_columns = len(keep)
            if args.truncate_min_coverage > 0.0:
                dense = local_coverage_mask(col_has_signal, args.truncate_window,
                                            args.truncate_min_coverage)
                truncated_columns += int((col_has_signal & ~dense).sum())
                keep = col_has_signal & dense
            dropped_no_signal += int((~keep).sum())

            # Per-camera straightness threshold. The two cameras fail
            # in different ways and one value cannot serve both --
            # measured, not assumed (see the argument help above).
            #
            # This mask MUST be folded into `keep` BEFORE the arrays
            # are subset below. An earlier version applied it after,
            # which counted every trimmed column while removing none:
            # the log reported 56426 trimmed and the output was
            # unchanged. Order matters here, and the point counts are
            # what reveal it -- the counter alone looked correct.
            if args.min_straightness_residual is not None:
                straight_threshold = args.min_straightness_residual
            else:
                straight_threshold = {"c1": args.straightness_residual_c1,
                                      "c2": args.straightness_residual_c2}.get(camera, 0.0)

            if straight_threshold > 0.0 and len(keep) >= args.straightness_window:
                straight_ok = straightness_mask(rows_full_for_straightness,
                                                args.straightness_window,
                                                straight_threshold)
                n_straight = int((keep & ~straight_ok).sum())
                straight_trimmed += n_straight
                if n_straight / max(int(keep.sum()), 1) > args.max_straight_fraction:
                    skipped_mostly_straight += 1
                    continue
                keep = keep & straight_ok

            columns, rows, sharpness = columns[keep], rows[keep], sharpness[keep]
            col_has_signal = col_has_signal[keep]
            if len(columns) == 0:
                skipped_all_no_signal += 1
                continue

            coverage = len(columns) / max(total_columns, 1)
            if args.min_coverage > 0.0 and coverage < args.min_coverage:
                skipped_low_coverage += 1
                low_coverage_detail.append((path.stem, coverage, len(columns), total_columns))
                continue

        for col, row, sharp in zip(columns, rows, sharpness):
            output_rows.append([
                path.stem, camera, capture_time, epoch,
                col, row, round(elevation_navd88, 4), round(spread, 4),
                round(gap_minutes, 1), is_stale,
                round(sharp, 4) if not np.isnan(sharp) else "",
                source_label,
            ])
        used += 1

    if not output_rows:
        print("No contour points produced -- check water-level source coverage and file naming.")
        sys.exit(1)

    # Tide-direction agreement, per camera-day. Computed AFTER all
    # frames are gathered because it needs consecutive frames -- it
    # cannot judge a frame in isolation, unlike every other filter here.
    agreement = tide_agreement_by_day(output_rows)
    rejected_days = set()
    if args.min_tide_agreement > 0.0 and agreement:
        print()
        print("Tide-direction agreement (detected line moves with the tide):")
        for (camera, date) in sorted(agreement):
            ok, total = agreement[(camera, date)]
            frac = ok / total if total else float("nan")
            if total < args.min_agreement_transitions:
                note = f"  (only {total} transitions -- not judged)"
            elif frac < args.min_tide_agreement:
                rejected_days.add((camera, date))
                note = "  <-- REJECTED, no directional signal"
            else:
                note = ""
            print(f"   {camera}  {date}   {ok:2d}/{total:2d}  {100*frac:5.0f}%{note}")

        if rejected_days:
            before = len(output_rows)
            output_rows = [r for r in output_rows
                           if (r[1], local_date(r[3])) not in rejected_days]
            print(f"   dropped {before - len(output_rows)} point(s) from "
                  f"{len(rejected_days)} camera-day(s)")
            if not output_rows:
                print("No contour points left after the agreement filter.")
                sys.exit(1)
        print()

    header = [
        "source_file", "camera", "capture_time_utc", "capture_epoch",
        "pixel_column", "pixel_row", "tide_elevation_navd88", "tide_model_spread",
        "tide_gap_minutes", "tide_is_stale", "peak_sharpness",
        "water_level_source",
    ]

    if args.setup_coef is not None and not args.waves:
        print("ERROR: --setup-coef needs --waves (the wave height and period per frame).")
        sys.exit(1)

    # Offshore waves per frame, appended as the LAST columns so every
    # positional index used above (r[1], r[2], ...) is unchanged.
    if args.waves:
        ep_hs, hs, ep_tp, tp = load_wave_archive(args.waves)
        gap_s = args.wave_max_gap_minutes * 60.0
        per_frame = {}
        for r in output_rows:
            if r[0] not in per_frame:
                h = nearest_within(ep_hs, hs, r[3], gap_s)
                t = nearest_within(ep_tp, tp, r[3], gap_s)
                per_frame[r[0]] = ("" if h is None else round(h, 2),
                                   "" if t is None else round(t, 1))
            r.extend(per_frame[r[0]])
        header += ["offshore_hs_m", "offshore_tp_s"]

        if args.setup_coef is not None:
            # Wave setup per frame; frames without a height AND period
            # get no correction and a blank setup column, so they can be
            # told apart from a genuine zero.
            no_waves = set()
            for r in output_rows:
                hs, tp = r[-2], r[-1]
                if hs != "" and tp != "":
                    setup = args.setup_coef * np.sqrt(float(hs) * 9.81 * float(tp) ** 2 / (2 * np.pi))
                    r.extend([round(float(setup), 4), round(r[6] + float(setup), 4)])
                else:
                    no_waves.add(r[0])
                    r.extend(["", r[6]])
            header += ["setup_correction_m", "beach_elevation_navd88"]
            setups = [r[-2] for r in output_rows if r[-2] != ""]
            print(f"Setup correction         : C = {args.setup_coef}, "
                  + (f"{min(setups):.3f}-{max(setups):.3f} m" if setups else "none applied")
                  + (f"; {len(no_waves)} frame(s) without wave height+period left uncorrected"
                     if no_waves else ""))
        tagged = [v[0] for v in per_frame.values() if v[0] != ""]
        print(f"Offshore waves           : {len(tagged)}/{len(per_frame)} frame(s) tagged"
              + (f", Hs {min(tagged):.2f}-{max(tagged):.2f} m, "
                 f"median {float(np.median(tagged)):.2f} m" if tagged else ""))

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(output_rows)

    used = len({r[0] for r in output_rows})
    print(f"Frames used              : {used}")
    print(f"Skipped (no epoch prefix): {skipped_no_epoch}")
    print(f"Skipped (no c1/c2 in name): {skipped_no_camera}")
    if args.include_no_signal:
        print(f"No-signal columns        : KEPT (--include-no-signal)")
    else:
        print(f"Dropped (no signal)      : {dropped_no_signal} column(s) across all frames")
        if straight_trimmed or skipped_mostly_straight:
            if args.min_straightness_residual is not None:
                thresh_note = f"{args.min_straightness_residual} px (both cameras)"
            else:
                thresh_note = (f"c1 {args.straightness_residual_c1} px, "
                               f"c2 {args.straightness_residual_c2} px")
            print(f"Trimmed (implausibly straight): {straight_trimmed} column(s) -- below "
                  f"{thresh_note}, i.e. the solver was drawing, not detecting")
            if skipped_mostly_straight:
                print(f"Skipped (mostly straight) : {skipped_mostly_straight} frame(s) -- over "
                      f"{100*args.max_straight_fraction:.0f}% of columns straight")
        if args.truncate_min_coverage > 0.0:
            print(f"Truncated (sparse region): {truncated_columns} column(s) -- passed the "
                  f"per-column test but sat in a neighbourhood below "
                  f"{100*args.truncate_min_coverage:.0f}% coverage")
        if skipped_all_no_signal:
            print(f"Skipped (all columns flagged no-signal): {skipped_all_no_signal} frame(s)")
        if skipped_low_coverage:
            print(f"Skipped (coverage < {100*args.min_coverage:.0f}%): "
                  f"{skipped_low_coverage} frame(s) -- too fragmentary to be a waterline")
            for stem, cov, kept_n, tot_n in sorted(low_coverage_detail,
                                                   key=lambda t: t[1])[:8]:
                print(f"      {100*cov:5.1f}%  ({kept_n:5d}/{tot_n:5d})  {stem}")
            if len(low_coverage_detail) > 8:
                print(f"      ... and {len(low_coverage_detail) - 8} more")
    print(f"Skipped (tide gap > {args.max_gap_minutes:.0f} min): {skipped_stale}"
          + ("" if not args.include_stale else " -- included anyway (see --include-stale)"))
    print()
    print(f"Wrote {len(output_rows)} contour point(s) from {used} frame(s) to {output_path}")
    print()
    print("NOTE: pixel_column/pixel_row are still in PIXEL space, not real-world")
    print("coordinates -- georectification (camera calibration) is still needed")
    print("before this can be gridded into an actual elevation map.")


if __name__ == "__main__":
    main()

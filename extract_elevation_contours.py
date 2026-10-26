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
    --timing-mode MODE    legacy (default) or mid: when in the exposure the
                          level is taken (timex_timing.py)

WATER-LEVEL TIMING (timex_timing.py). In the default, legacy, mode a
frame's level is taken at its file-name epoch, as always. With
--timing-mode mid (or TIMING_MODE=mid) the level, the waves and the
stale-gap check are taken at the middle of the exposure: epoch +
TIMEX_MID_OFFSET_S (300 s) for a timex, + 0 s for a snap. capture_time_utc
and capture_epoch stay the file-name time, so nothing that joins on them
changes; each row then also says, in a last column level_time_offset_s,
how many seconds later its level was taken. A detection is only used in
the mode it was detected in (its Timing_Mode column; none = legacy): a
blended legacy line given a mid-exposure level is 11.5 min behind it,
worse than either mode, so such files are left out and counted.
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

import timex_timing
from setup_model import (SETUP_MODELS, STOCKDON_SETUP_K, EvolvingSlope, deshoal_depth, model_tag,
                         setup_c, stockdon_frame)


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


# THE WAVE "CURRENCY" OF THE SETUP COEFFICIENT. C (waterline_timex_cron.sh
# SETUP_COEF = 0.037) was fitted on contours tagged from
# archive/waves_marconi.csv (marconi_waves.py, USE_MARCONI_WAVES=1):
# offshore_hs_m = hs_best, which is ADCP wh_4061 currency (camera wave
# models trained on that ADCP, buoys converted to it), offshore_tp_s = tp_s
# (NDBC 44008 peak period). The same formula on another Hs is another
# setup: the raw 44008 wvht_m overstates sheltered south-west seas, and the
# DEM moves ~0.06 m per 0.01 of effective C. The cron falls back to the raw
# buoy archive whenever waves_marconi.csv is missing (owg.sh stopped, or a
# copy of the repository without it), and nothing said so; every contour
# row now names its wave columns (wave_source) and a WARNING is printed
# when they are the raw buoy's (a NOTE for the other kinds that are not
# this currency). The kinds, told apart by the columns:
#   marconi     hs_best / tp_s             -- the currency C was fitted in
#   adcp_equiv  wvht_m / dpd_s + hs_sigma_m -- historical_forcing.py's
#               waves.csv: ADCP-equivalent Hs (C's Hs currency); its Tp is
#               ADCP-equivalent, not 44008's (historical_forcing.py's report
#               gives the setup bias that implies)
#   buoy        wvht_m / dpd_s + apd_s      -- fetch_buoy_waves.py, the raw
#               offshore buoy (its archive always has apd_s): NOT the
#               currency C was fitted in
#   unknown     wvht_m / dpd_s, neither     -- a file the columns do not
#               name, e.g. process_chelsea.py's adcp_waves.csv (the ADCP's
#               own wh_4061 / peak period, which an earlier version of this
#               rule called the raw buoy and warned about): a NOTE, no
#               WARNING
WAVE_CURRENCY_FITTED = "marconi:hs_best/tp_s"


def wave_columns(columns):
    """
    The (height, period) columns load_wave_archive reads from a wave file
    with these column names, and the file's currency tag
    ('<kind>:<height>/<period>', see WAVE_CURRENCY_FITTED). The raw buoy is
    recognised by fetch_buoy_waves.py's apd_s column, not by wvht_m / dpd_s
    alone: other writers use those two names too (process_chelsea.py's
    ADCP record), and calling them the raw buoy put a false WARNING in
    their logs.
    """
    columns = set(columns)
    if "hs_best" in columns:          # marconi_waves.py's local best estimate
        return "hs_best", "tp_s", "marconi:hs_best/tp_s"
    kind = ("adcp_equiv" if "hs_sigma_m" in columns else
            "buoy" if "apd_s" in columns else "unknown")
    return "wvht_m", "dpd_s", f"{kind}:wvht_m/dpd_s"


def wave_currency_note(tag):
    """One line on a wave_source tag against the currency C was fitted in:
    (is_fitted_currency, text)."""
    if tag == WAVE_CURRENCY_FITTED:
        return True, (f"{tag} -- the currency the setup coefficient was fitted in "
                      f"(hs_best: ADCP wh_4061 Hs; tp_s: NDBC 44008 peak period)")
    if tag.startswith("adcp_equiv:"):
        return False, (f"{tag} -- ADCP-equivalent Hs (the setup coefficient's Hs currency) but "
                       f"ADCP-equivalent Tp, not the NDBC 44008 peak period C was fitted with "
                       f"(historical_forcing.py's report gives the setup bias that implies)")
    if tag.startswith("unknown:"):
        return False, (f"{tag} -- a wave file its columns do not name (not fetch_buoy_waves.py's raw "
                       f"buoy archive, no apd_s; not historical_forcing.py's waves.csv, no hs_sigma_m), "
                       f"e.g. process_chelsea.py's adcp_waves.csv (the ADCP's own Hs and peak period): "
                       f"not known to be the currency the setup coefficient was fitted in "
                       f"({WAVE_CURRENCY_FITTED})")
    return False, (f"{tag} -- NOT the currency the setup coefficient was fitted in "
                   f"({WAVE_CURRENCY_FITTED}): raw offshore buoy Hs/Tp overstate sheltered "
                   f"south-west seas, so the setup (and the --max-hs filter) differ")


def load_wave_archive(path):
    """
    Reads the offshore wave archive written by fetch_buoy_waves.py, or the
    hourly Marconi best estimate written by marconi_waves.py.
    Returns (epochs_hs, hs, epochs_tp, tp), each sorted by time, with
    rows lacking a value dropped separately for height and period.
    Which columns are read: wave_columns().
    """
    df = pd.read_csv(path)
    # marconi_waves.py's local best estimate (hs_best, tp_s) or the raw buoy
    hcol, tcol, _tag = wave_columns(df.columns)
    hs = df.dropna(subset=[hcol]).sort_values("epoch")
    tp = df.dropna(subset=[tcol]).sort_values("epoch")
    return (hs["epoch"].to_numpy(float), hs[hcol].to_numpy(float),
            tp["epoch"].to_numpy(float), tp[tcol].to_numpy(float))


def stockdon_setups(output_rows, args, wave_tag):
    """
    The slope-aware setup of every frame with a wave height and period (--setup-model stockdon;
    setup_model.py has the model, the deshoaling question and the fallback): Stockdon et al. (2006)
    k * beta_f * sqrt(H0 * L0), beta_f measured on the DEM nearest the frame's day over its swash band
    (still water + setup +/- S/2, iterated), else C / k. Called after the wave columns are appended,
    so a row's last two values are its offshore_hs_m and offshore_tp_s.
    -> {'frames': {source_file: (setup m, beta_f, source)}, 'tag', 'how', 'lines'}.
    """
    k = args.setup_stockdon_k
    depth, how = deshoal_depth(args.setup_deshoal_depth, wave_tag)
    fallback = args.setup_coef / k
    slopes = EvolvingSlope(args.setup_slope_dem)
    frames, why = {}, defaultdict(int)
    for r in output_rows:
        # a Tp of 0 (or a setup that is not a number) is a frame without waves, left uncorrected
        if r[0] in frames or r[-2] == "" or r[-1] == "" or not float(r[-1]) > 0:
            continue
        fs, why_none = slopes.for_epoch(r[3])
        st = stockdon_frame(r[6], r[-2], r[-1], fs, fallback, k, depth, why_none=why_none)
        if not np.isfinite(st[0]):
            continue
        frames[r[0]] = st
        why[st[2]] += 1
    betas = [v[1] for v in frames.values() if v[2] == "dem"]
    lines = [f"Foreshore slope          : from {slopes.what}: "
             + (f"{len(betas)} frame(s), beta_f median {np.median(betas):.3f} (p10-p90 "
                f"{np.percentile(betas, 10):.3f}-{np.percentile(betas, 90):.3f})" if betas else "no frame")]
    fell = {w: n for w, n in why.items() if w != "dem"}
    if fell:
        lines.append(f"{'WARNING' if not betas else 'NOTE'}: {sum(fell.values())} frame(s) took the fallback "
                     f"slope beta_f = C/k = {fallback:.3f} (the C model's setup): "
                     + ", ".join(f"{w} {n}" for w, n in sorted(fell.items())))
    return {"frames": frames, "tag": model_tag(k, depth), "how": how, "lines": lines}


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


# NAVD88 height of the GNSS-R antenna. gnssrefl computes each water
# level as (station orthometric height - reflector height), so any error
# in that height shifts EVERY water level -- and every waterline and DEM
# elevation -- by the same amount.
#
# The height in the gnssrefl station setup, 18.665 m, came from a
# CSRS-PPP solution (NRCan, 24 h of RINEX, 2026-07-11). CSRS-PPP reports
# orthometric heights in CGVD2013 (Canadian datum, geoid CGG2013a), not
# NAVD88. NGS OPUS on the same RINEX (2026-09-28; NAD83(2011) epoch
# 2010.0, GEOID18) gives 19.014 m NAVD88 (+/-0.061 m, mostly geoid model
# uncertainty); the two solutions agree on the ITRF2020 ellipsoid height
# to 2 mm (-10.025 vs -10.023 m), so the antenna position is not in
# question -- only the height system. GNSS-R levels were therefore
# 0.349 m LOW, consistent with the Chatham gauge (-0.37 m, including a
# harbour effect) and the tide-model offset (0.26 m).
#
# The correction applied is GNSSR_ANTENNA_NAVD88_M minus the height
# written in the spline file's header, so if the gnssrefl station file is
# later corrected to 19.014 m, the correction falls to zero by itself
# rather than being applied twice.
GNSSR_ANTENNA_NAVD88_M = 19.014


def load_gnssr_spline(path, qc_reference=None):
    """
    Loads gnssrefl subdaily spline output as a water-level source.

    Columns (1-indexed): 1 MJD, 2 RH(m), 3 YYYY, 4 MM, 5 DD, 6 HH,
    7 MM, 8 SS, 9 quasi-sea-level(m). Column 9 is computed upstream as
    (station orthometric height - reflector height). The levels are
    shifted onto NAVD88 -- the datum of the camera EO z-coordinates --
    using GNSSR_ANTENNA_NAVD88_M (see there: the header height was a
    CGVD2013 value, 0.349 m below NAVD88).

    It is also a MEASUREMENT at this station, so it includes storm
    surge, wind setup and local bias that an astronomical tide model
    cannot predict. Against ~47 days of overlap a global tide model
    differed from the (uncorrected) GNSS-R by about 0.26 m in the mean
    with a 0.08 m standard deviation.

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

    print("Water-level source: GNSS-R MEASURED (gnssrefl spline)")
    print(f"  Loaded {len(epoch)} readings, spanning "
          f"{datetime.fromtimestamp(epoch[0], tz=timezone.utc)} to "
          f"{datetime.fromtimestamp(epoch[-1], tz=timezone.utc)}")
    if hortho is None:
        print("  WARNING: no antenna orthometric height in the file header, so the NAVD88")
        print("  correction cannot be checked; levels used as written.")
    else:
        correction = GNSSR_ANTENNA_NAVD88_M - hortho
        print(f"  antenna height in file {hortho:.3f} m; NAVD88 (OPUS/GEOID18) "
              f"{GNSSR_ANTENNA_NAVD88_M:.3f} m")
        if abs(correction) >= 0.001:
            levels = levels + correction
            print(f"  -> {correction:+.3f} m added to every level to put it on NAVD88")
        else:
            print("  -> already NAVD88, no correction")
    if qc_reference:
        # Drop readings that fail quality control (gnssr_qc.py); they
        # become gaps, so nothing near them gets a wrong water level.
        from gnssr_qc import qc_filter
        epoch, levels, summary = qc_filter(epoch, levels, qc_reference)
        print("  " + summary)
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


def cliff_segments(columns, rows, window, cliff_px):
    """
    Label runs of a detected line separated by CLIFFS: places where the row
    changes by more than cliff_px within `window` columns. A waterline is
    smooth across the image (well under 1 px per column); a jump of tens of
    pixels in a few columns is the solver leaving the water for another
    edge -- on the winter imagery, the edges of snow patches, which gave
    staircase lines climbing the dry beach. Across a column gap wider than
    `window` nothing is compared (the gap is not a cliff, and does not
    split the run). Returns (labels, number of cliffs).
    """
    n = len(rows)
    labels = np.zeros(n, dtype=int)
    if n < 2:
        return labels, 0
    cliffs = 0
    run = 0
    j = 0
    for i in range(1, n):
        if columns[i] - columns[i - 1] > window:
            j = i
        else:
            while columns[i] - columns[j] > window:
                j += 1
            if abs(rows[i] - rows[j:i].min()) > cliff_px or abs(rows[i] - rows[j:i].max()) > cliff_px:
                if labels[i - 1] == run:          # first column past the jump: new run
                    run += 1
                    cliffs += 1
                j = i
        labels[i] = run
    return labels, cliffs


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
    parser.add_argument("--cliff-px", type=float, default=40.0,
                        help="A row change of more than this within --cliff-window columns is a "
                             "cliff: the line left the water for another edge (default 40; 0 "
                             "disables). A waterline changes by well under 1 px per column.")
    parser.add_argument("--cliff-window", type=int, default=10)
    parser.add_argument("--cliff-min-run", type=float, default=0.6,
                        help="With cliffs, keep the longest smooth run if it holds at least this "
                             "share of the frame's columns (default 0.6), else drop the frame.")
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
                        help="Offshore wave archive from fetch_buoy_waves.py (or "
                             "marconi_waves.py's waves_marconi.csv). Adds "
                             "offshore_hs_m and offshore_tp_s to every point, so rough-water "
                             "frames (whose waterline sits above still water by the wave "
                             "setup) can be left out of the DEM, and wave_source: which "
                             "columns of which kind of file tagged them (a WARNING when it "
                             "is the raw buoy, a NOTE for another currency than the one the "
                             "setup coefficient was fitted in; see WAVE_CURRENCY_FITTED).")
    parser.add_argument("--wave-max-gap-minutes", type=float, default=90.0,
                        help="Largest gap to the nearest wave record before a frame's wave "
                             "columns are left blank (default 90; the buoy reports every "
                             "10-60 min).")
    parser.add_argument("--gnssr-qc-reference", default=None,
                        help="Tide-gauge archive (fetch_tide_gauge.py) for GNSS-R quality control "
                             "(gnssr_qc.py): readings that fail are removed and become gaps, so "
                             "frames at those times are skipped as stale instead of getting a "
                             "wrong level. Without this option no QC is applied.")
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
    parser.add_argument("--setup-model", choices=SETUP_MODELS, default="c",
                        help="Which setup --setup-coef applies (setup_model.py). 'c' (default, the "
                             "station's): C*sqrt(Hs*L0) as above. 'stockdon' (opt-in; the station "
                             "stays on 'c' until evaluate_setup_model.py favours it on the surveys): "
                             "Stockdon et al. (2006) k*beta_f*sqrt(H0*L0), with the foreshore slope "
                             "beta_f measured for each frame on the evolving DEM "
                             "(--setup-slope-dem) over its swash band, and H0 the deep-water height "
                             "(--setup-deshoal-depth). A frame the DEM cannot serve takes beta_f = "
                             "C/k, i.e. the C model's setup. Adds the columns setup_beta_f and "
                             "setup_model.")
    parser.add_argument("--setup-slope-dem", default=None,
                        help="With --setup-model stockdon (required there): the DEM the foreshore "
                             "slope is measured on. A folder of dated DEMs (dem_from_contours.py "
                             "--series-dir, the cron's archive/dems): each frame takes the one whose "
                             "window is centred nearest its day, within 7 days. Or one DEM "
                             "(<stem> or <stem>_dem.asc) for every frame.")
    parser.add_argument("--setup-stockdon-k", type=float, default=STOCKDON_SETUP_K,
                        help="With --setup-model stockdon: the setup coefficient k (default 0.35, "
                             "Stockdon et al. 2006 Eq. 10; evaluate_setup_model.py refits it).")
    parser.add_argument("--setup-deshoal-depth", default="auto",
                        help="With --setup-model stockdon: the depth (m) the wave heights are "
                             "deshoaled from to Stockdon's deep-water H0 (linear theory, H0 = Hs/Ks). "
                             "'auto' (default): 21 m for the ADCP currency (waves_marconi.csv "
                             "hs_best, historical_forcing.py's waves), none for a raw buoy file; "
                             "'none': H0 = Hs.")
    parser.add_argument("--vertical-datum-offset", type=float, default=0.0,
                         help="Meters to ADD to every tide elevation, converting the tide "
                              "model's vertical datum (mean sea level or a geoid, for global "
                              "ocean models like EOT20/GOT/FES) into NAVD88 -- the datum the "
                              "camera EO z-coordinates use. The preferred way to obtain this "
                              "is compare_gnssr_to_tidemodel.py, which measures it against "
                              "GNSS-R water levels on NAVD88 (about +0.09 m at this station "
                              "after the GNSS-R datum correction; it was -0.25 m against the "
                              "uncorrected GNSS-R). NOAA's VDatum tool is an alternative but gives a "
                              "modelled rather than measured value. Defaults to 0.0, which "
                              "is almost certainly WRONG -- see the warning printed if unset.")
    parser.add_argument("--timing-mode", default=None, choices=timex_timing.TIMING_MODES,
                        help="When in the exposure each frame's water level is taken (timex_timing.py). "
                             "legacy (the default): at the file-name epoch, the exposure's start, as "
                             "always. mid: the level, the waves and the stale-gap check at epoch + "
                             f"{timex_timing.TIMEX_MID_OFFSET_S} s for a timex (mid-exposure; + 0 s for a "
                             "snap), for detections made with waterline_detector_v5.py --timing-mode mid "
                             "(no temporal blend). Detections of the other mode are left out and counted. "
                             "Default: the TIMING_MODE environment variable, else legacy.")
    args = parser.parse_args()
    try:
        timing = timex_timing.timing_mode(args.timing_mode)
    except ValueError as e:
        parser.error(str(e))
    if args.setup_model == "stockdon" and (args.setup_coef is None or not args.waves
                                           or not args.setup_slope_dem):
        print("ERROR: --setup-model stockdon needs --setup-coef (the C of the frames the DEM cannot "
              "serve, and the switch that turns the setup on), --waves and --setup-slope-dem.")
        sys.exit(1)
    if args.setup_model == "stockdon":
        # checked before any frame is read, not after the whole archive has been
        try:
            deshoal_depth(args.setup_deshoal_depth, "")
        except ValueError:
            print(f"ERROR: --setup-deshoal-depth {args.setup_deshoal_depth!r}: 'auto', 'none' or metres.")
            sys.exit(1)
        if not 0 < args.setup_stockdon_k < float("inf"):
            print(f"ERROR: --setup-stockdon-k {args.setup_stockdon_k}: must be a finite number above 0.")
            sys.exit(1)

    using_gnssr = looks_like_gnssr_spline(args.tide_model_csv)

    if using_gnssr:
        tide_timestamps, tide_levels, tide_spreads, _hortho = load_gnssr_spline(
            args.tide_model_csv, args.gnssr_qc_reference)
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
        print("At this station that comes out near +0.09 m with GNSS-R on NAVD88 (it was")
        print("-0.25 m before the GNSS-R datum correction of +0.349 m) and varied by about")
        print("0.04 m between comparison windows -- so treat it as approximate. Where")
        print("GNSS-R covers the imagery directly it is better still: no offset is needed.")
        print("=" * 90)
        print()

    processed_dir = Path(args.processed_dir)
    shoreline_files = sorted(processed_dir.glob("*.csv"))

    if not shoreline_files:
        print(f"No shoreline CSVs found in {processed_dir}")
        sys.exit(1)

    print(f"Found {len(shoreline_files)} processed shoreline file(s) in {processed_dir}")
    print(f"Water-level time         : {timex_timing.describe(timing, level_only=True)}")
    print()

    output_rows = []
    # Per frame: the epoch its water level (and waves) were taken at --
    # the file-name epoch in legacy mode -- and its offset from it.
    level_epoch_of, level_offset_of = {}, {}
    skipped_other_mode = defaultdict(int)
    used, skipped_no_epoch, skipped_no_camera, skipped_stale = 0, 0, 0, 0
    dropped_no_signal = 0
    skipped_all_no_signal = 0
    skipped_low_coverage = 0
    low_coverage_detail = []
    truncated_columns = 0
    straight_trimmed = 0
    skipped_mostly_straight = 0
    skipped_staircase = 0
    cliff_trimmed = 0

    for path in shoreline_files:
        epoch = extract_epoch_from_filename(path.name)
        if epoch is None:
            skipped_no_epoch += 1
            continue

        camera = detect_camera(path.name)
        if camera is None:
            skipped_no_camera += 1
            continue

        # A line is given a level only in the timing mode it was detected
        # in: the level time and the blend go together (timex_timing.py).
        det_mode = timex_timing.detection_timing_mode(path)
        if det_mode != timing:
            skipped_other_mode[det_mode] += 1
            continue

        level_offset = timex_timing.level_time_offset_s(path.name, timing)
        level_epoch = epoch + level_offset
        elevation, spread, gap_minutes = interpolate_tide(
            level_epoch, tide_timestamps, tide_levels, tide_spreads
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

            # Staircase check: keep the longest smooth run, or drop the frame.
            if args.cliff_px > 0 and len(columns) > 1:
                labels, n_cliffs = cliff_segments(columns, rows, args.cliff_window, args.cliff_px)
                if n_cliffs:
                    sizes = np.bincount(labels)
                    best = int(np.argmax(sizes))
                    if sizes[best] < args.cliff_min_run * len(columns):
                        skipped_staircase += 1
                        continue
                    run_keep = labels == best
                    cliff_trimmed += int((~run_keep).sum())
                    columns, rows, sharpness = columns[run_keep], rows[run_keep], sharpness[run_keep]
                    col_has_signal = col_has_signal[run_keep]

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
        level_epoch_of[path.stem] = level_epoch
        level_offset_of[path.stem] = level_offset
        used += 1

    if skipped_other_mode:
        for mode, n in sorted(skipped_other_mode.items()):
            print(f"WARNING: {n} detection(s) made in timing mode {mode!r} left out: this run takes the "
                  f"level in mode {timing!r}, and a line gets a level only in the mode it was detected "
                  f"in (timex_timing.py). Point --processed-dir at that mode's detections, or use "
                  f"--timing-mode {mode}.")

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
        hcol, tcol, wave_tag = wave_columns(pd.read_csv(args.waves, nrows=0).columns)
        fitted, note = wave_currency_note(wave_tag)
        print(f"Wave file                : {args.waves} (columns {hcol}, {tcol})")
        print(f"Wave currency            : {note}")
        if not fitted and wave_tag.startswith("buoy:"):
            print(f"WARNING: the frames are tagged with the raw buoy waves ({wave_tag}), not the "
                  f"currency the setup coefficient was fitted in ({WAVE_CURRENCY_FITTED})"
                  + (f": the setup with C = {args.setup_coef} is NOT the station's calibrated setup"
                     if args.setup_coef else "")
                  + ". Was archive/waves_marconi.csv (owg.sh / marconi_waves.py) missing?")
        elif not fitted and wave_tag.startswith("unknown:"):
            print(f"NOTE: the wave file's kind is not known from its columns, so not known to be the "
                  f"currency C was fitted in ({WAVE_CURRENCY_FITTED}); see the line above.")
        elif not fitted:
            print(f"NOTE: wave period currency differs from the one C was fitted in "
                  f"({WAVE_CURRENCY_FITTED}); see the line above.")
        gap_s = args.wave_max_gap_minutes * 60.0
        per_frame = {}
        for r in output_rows:
            if r[0] not in per_frame:
                # the waves at the time the level was taken (the file-name
                # epoch in legacy mode, mid-exposure in mid mode)
                t_level = level_epoch_of.get(r[0], r[3])
                h = nearest_within(ep_hs, hs, t_level, gap_s)
                t = nearest_within(ep_tp, tp, t_level, gap_s)
                per_frame[r[0]] = ("" if h is None else round(h, 2),
                                   "" if t is None else round(t, 1))
            r.extend(per_frame[r[0]])
        header += ["offshore_hs_m", "offshore_tp_s"]

        if args.setup_coef is not None:
            # Wave setup per frame; frames without a height AND period
            # get no correction and a blank setup column, so they can be
            # told apart from a genuine zero.
            no_waves = set()
            slope_aware = None
            if args.setup_model == "stockdon":
                slope_aware = stockdon_setups(output_rows, args, wave_tag)
            for r in output_rows:
                hs, tp = r[-2], r[-1]
                if slope_aware is not None:
                    # the slope-aware model's setup and the slope it used (setup_model.py)
                    s = slope_aware["frames"].get(r[0])
                    if s is not None:
                        r.extend([round(s[0], 4), round(r[6] + s[0], 4), round(s[1], 4), slope_aware["tag"]])
                    else:
                        no_waves.add(r[0])
                        r.extend(["", r[6], "", slope_aware["tag"]])
                elif hs != "" and tp != "":
                    setup = setup_c(float(hs), float(tp), args.setup_coef)
                    r.extend([round(float(setup), 4), round(r[6] + float(setup), 4)])
                else:
                    no_waves.add(r[0])
                    r.extend(["", r[6]])
            header += ["setup_correction_m", "beach_elevation_navd88"]
            if slope_aware is not None:
                header += ["setup_beta_f", "setup_model"]
            isu = header.index("setup_correction_m")
            setups = [r[isu] for r in output_rows if r[isu] != ""]
            what = (f"C = {args.setup_coef}" if slope_aware is None else
                    f"slope-aware (Stockdon et al. 2006, k = {args.setup_stockdon_k:g}; {slope_aware['how']})")
            print(f"Setup correction         : {what}, "
                  + (f"{min(setups):.3f}-{max(setups):.3f} m" if setups else "none applied")
                  + (f"; {len(no_waves)} frame(s) without wave height+period left uncorrected"
                     if no_waves else ""))
            if slope_aware is not None:
                for line in slope_aware["lines"]:
                    print(line)
        # Which wave columns tagged the frames, on every row (as
        # water_level_source names the level): the last column, so every
        # positional index above is unchanged. A '#' header line is not an
        # option -- georectify.py and the other readers take the first line
        # as the column names -- and the rows carry it on through the
        # consistency filter and georectify.py to dem_from_contours.py,
        # which records it in <stem>_info.json.
        for r in output_rows:
            r.append(wave_tag)
        header += ["wave_source"]
        tagged = [v[0] for v in per_frame.values() if v[0] != ""]
        print(f"Offshore waves           : {len(tagged)}/{len(per_frame)} frame(s) tagged"
              + (f", Hs {min(tagged):.2f}-{max(tagged):.2f} m, "
                 f"median {float(np.median(tagged)):.2f} m" if tagged else ""))

    # Outside legacy mode every row says when its level was taken, as the
    # last column (so every positional index above, and wave_source's
    # place before it, are unchanged); legacy rows are exactly as before.
    if timing != "legacy":
        for r in output_rows:
            r.append(level_offset_of[r[0]])
        header += [timex_timing.LEVEL_OFFSET_COLUMN]
        offsets = sorted(set(level_offset_of[r[0]] for r in output_rows))
        print(f"Water-level time         : {timing}: level (and waves) taken at capture_epoch + "
              + ", ".join(f"{o} s" for o in offsets) + f" ({timex_timing.LEVEL_OFFSET_COLUMN})")

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
        if skipped_staircase or cliff_trimmed:
            print(f"Staircase lines          : {skipped_staircase} frame(s) skipped, {cliff_trimmed} "
                  f"column(s) trimmed -- jumps over {args.cliff_px:g} px within "
                  f"{args.cliff_window} columns (the solver leaving the water, e.g. for snow edges)")
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
    if skipped_other_mode:
        print(f"Skipped (detected in another timing mode than {timing}): "
              + ", ".join(f"{n} {m}" for m, n in sorted(skipped_other_mode.items())))
    print()
    print(f"Wrote {len(output_rows)} contour point(s) from {used} frame(s) to {output_path}")
    print()
    print("NOTE: pixel_column/pixel_row are still in PIXEL space, not real-world")
    print("coordinates -- georectification (camera calibration) is still needed")
    print("before this can be gridded into an actual elevation map.")


if __name__ == "__main__":
    main()

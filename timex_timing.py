#!/usr/bin/env python3
"""
When A Timex Waterline's Water Level Is Taken (the timing switch)
===================================================================
One place for the timing of the image products: the middle of a timex
exposure, and the switch (TIMING_MODE) between the timing the station has
always used and the corrected one.

THE COLLECTION. An Argus/CoastCam collection is a 10-minute exposure: the
snap is one frame, the timex the mean of every frame, the ras.tiff
timestack 1200 rows at 2 Hz (Holman & Stanley 2007). The file name carries
the collection's START: the names fall on :00:00 and :30:00, which fits a
scheduled start and not an end. So the middle of a timex is the file-name
epoch + TIMEX_MID_OFFSET_S. runup_from_timestack.py, timestack_wave_period.py,
relabel_manifests.py and reaim_figure.py already took their water levels
and waves there; they read the constant from here. NEITHER the start nor
the 10 minutes has been read off the station's collection settings: both
are the repository's assumption (see docs/VALIDATION.md section 2.7 for the
check on the NUC). If the exposure is not 600 s, this constant is the one
number to change.

WHY THE WATERLINES NEED A SWITCH. Since v6.8 the timex waterlines were
given the water level at the file-name epoch, the START of the exposure,
while the detector blended 30% of the previous saved frame (30 min older)
into each line (waterline_detector_v5.py temporal_filter). The blend is
RECURSIVE -- the stored previous line is itself blended -- so a line
trails mid-exposure by 0, 9, 11.7, ... up to 12.9 min along a chain of
frames, and is 4 to 7.9 min OLDER than the level it is given: high on the
flood, low on the ebb, by that lag times dz/dt (up to ~0.1 m at spring
rates). The fix is three changes that only work TOGETHER (the review's
reproduction, known finding C-B1):
  * the level (and the waves, and the stale-gap check) at mid-exposure;
  * the blend off -- a causal blend cannot be made lag-free;
  * timex bias corrections re-derived with the blend off: the built-in
    ones were fitted on 11 Sep 2026 frames that were all on the flood, so
    they absorbed ~2.8 px (c1) / ~2.6 px (c2) of blend lag.
Any one alone leaves an error or makes it worse: the level at mid-exposure
with the blend still on puts the line 11.5 min behind its level.

THE SWITCH. TIMING_MODE, one value for every step:
  legacy  (the default) the level at the file-name epoch, the blend on,
          the timex bias of calibration/timex_bias_legacy.json (the
          built-in v6.8 values): the station's products as they have
          always been made, bit for bit.
  mid     the level, waves and stale-gap check at epoch + 300 s for a
          timex (+0 s for a snap, a single frame), the blend off, the
          timex bias of calibration/timex_bias_mid.json, which only
          recalibrate_timing.py writes (on the NUC, from the hand-clicked
          ground truth re-detected with the blend off). No such file, no
          mid detection with the bias correction on.
The detector and extract_elevation_contours.py both read it: from their
--timing-mode option, else from the TIMING_MODE environment variable, else
legacy. The cron sets it in one line, passes it to both and exports it;
process_historical.py and process_chelsea.py resolve it once
(--timing-mode, else TIMING_MODE) and give the same mode to both;
detect_original_view.py (pinned by the GCP-free pre-registration, so not
edited) runs both as subprocesses that inherit TIMING_MODE, which
survey_products.py passes on. A mid-mode detection CSV says
so in its Timing_Mode column, and the contour step leaves out (and
counts) every detection of the other mode, so the level time and the
blend can never be applied one without the other. recalibrate_timing.py
is the path to mid on the station (docs/VALIDATION.md section 2.7): it
checks the exposure, derives timex_bias_mid.json, compares the two modes
on a window and backfills the mid archive.

Reference: Holman, R.A., Stanley, J., 2007. The history and technical
capabilities of Argus. Coastal Engineering 54(6-7), 477-491.
"""

import csv
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent

# The middle of the 10-minute collection, in seconds after the file-name
# epoch (the collection's start). The ONE definition: everything that
# takes a level or waves at mid-exposure reads it from here
# (recalibrate_timing.py --steps exposure checks it on the station;
# test_timing_mode.py pins 300, so change it there too).
TIMEX_MID_OFFSET_S = 300

TIMING_MODES = ("legacy", "mid")
DEFAULT_TIMING_MODE = "legacy"
# The environment variable both the detector and the contour step read
# when no --timing-mode is given (the cron exports it).
TIMING_ENV = "TIMING_MODE"

# A detection CSV made in a non-legacy mode carries this column (the mode
# on every row); legacy detections have none, as before.
TIMING_COLUMN = "Timing_Mode"
# A contour row made in a non-legacy mode carries this last column: how
# many seconds after capture_epoch its water level (and waves) were taken.
LEVEL_OFFSET_COLUMN = "level_time_offset_s"

CALIBRATION_DIR = HERE / "calibration"


def timing_mode(value=None, environ=None):
    """
    The timing mode in force: `value` (a --timing-mode option) if given,
    else the TIMING_MODE environment variable, else 'legacy'. An unknown
    value is an error, never a guess: a typo must not silently run the
    other mode.
    """
    env = os.environ if environ is None else environ
    mode = value if value not in (None, "") else env.get(TIMING_ENV, "") or DEFAULT_TIMING_MODE
    mode = str(mode).strip().lower()
    if mode not in TIMING_MODES:
        raise ValueError(f"timing mode {mode!r} is not one of {', '.join(TIMING_MODES)} "
                         f"(--timing-mode, or the {TIMING_ENV} environment variable)")
    return mode


def level_time_offset_s(name, mode):
    """
    Seconds after the file-name epoch at which a frame's water level is
    taken: 0 in legacy mode (the file-name epoch, as always); in mid mode
    TIMEX_MID_OFFSET_S for a timex (the middle of its exposure) and 0 for a
    snap, which is a single frame at the start.
    """
    if mode == "legacy":
        return 0
    return TIMEX_MID_OFFSET_S if ".timex." in Path(name).name.lower() else 0


def describe(mode, level_only=False):
    """One line on what the mode does, for the logs. level_only: only the
    water-level time (the contour step, which sets neither the blend nor
    the bias)."""
    if mode == "legacy":
        return ("legacy: water level at the file-name epoch (the exposure's start)"
                + ("" if level_only else ", temporal blend on, timex bias calibration/timex_bias_legacy.json"))
    return (f"{mode}: water level, waves and stale-gap check at the file-name epoch + "
            f"{TIMEX_MID_OFFSET_S} s for a timex (mid-exposure; + 0 s for a snap)"
            + ("" if level_only else f", temporal blend off, timex bias calibration/timex_bias_{mode}.json"))


def timex_bias_path(mode, calibration_dir=None):
    """calibration/timex_bias_<mode>.json."""
    return Path(calibration_dir or CALIBRATION_DIR) / f"timex_bias_{mode}.json"


def read_timex_bias(path, mode=None):
    """
    {'c1': ((x_fraction, correction_px), ...), 'c2': (...)} from a timex
    bias file, the points exactly as written (floats parsed from the JSON
    text, so -11.03 is the same double as the literal -11.03 in the
    detector). The file's 'mode' must be `mode` when one is given: a file
    fitted for one timing (blend on or off) is never applied in the other.
    Its 'about' is provenance. Raises ValueError if the mode differs, a
    camera is missing or a point is malformed.
    """
    data = json.loads(Path(path).read_text())
    if mode is not None and data.get("mode") != mode:
        raise ValueError(f"{path} is the timex bias of timing mode {data.get('mode')!r}, "
                         f"not of {mode!r}")
    out = {}
    for cam in ("c1", "c2"):
        pts = data.get(cam)
        if not pts:
            raise ValueError(f"{path}: no '{cam}' points")
        out[cam] = tuple((float(x), float(c)) for x, c in pts)
        xs = [p[0] for p in out[cam]]
        if xs != sorted(xs) or not all(0.0 <= x <= 1.0 for x in xs):
            raise ValueError(f"{path}: '{cam}' x fractions must rise within 0..1")
    return out


def load_timex_bias(mode, path=None, calibration_dir=None):
    """(the mode's points as read_timex_bias, or None if the file does not
    exist; the path). `path` overrides calibration/timex_bias_<mode>.json."""
    p = Path(path) if path else timex_bias_path(mode, calibration_dir)
    if not p.exists():
        return None, p
    return read_timex_bias(p, mode), p


def detection_timing_mode(path):
    """
    The timing mode a detection CSV was made in: its Timing_Mode column
    (first row), or 'legacy' when the column is absent (every detection
    made before the switch existed, and every legacy-mode one since).
    """
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        if TIMING_COLUMN not in (reader.fieldnames or []):
            return "legacy"
        for row in reader:
            return (row.get(TIMING_COLUMN) or "legacy").strip().lower()
    return "legacy"


def archive_suffix(mode):
    """The station archive of a mode's detections is archive/processed_timex<suffix>:
    '' for legacy (the archive as it has always been), '_mid' for mid."""
    return "" if mode == "legacy" else f"_{mode}"

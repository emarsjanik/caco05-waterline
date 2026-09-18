#!/usr/bin/env python3
# =====================================================================
#
#  USGS ARGUS SHORELINE DETECTOR
#  Version 6.9
#
#  CHANGES IN THIS VERSION:
#    - PER-CAMERA BAND CONTRAST FLOOR (CameraProfile.band_contrast_floor).
#      Whole-image contrast does not detect fog: a fogged C1 frame
#      scored global std 31.8 -- higher than many clear frames --
#      while the cropped band it is actually detected on scored 5.3.
#      Sky and dune keep enough variation to mask a flattened beach.
#      Measured across 108 archived frames, the eight lowest in-band
#      values were ALL from Sep 13 before any other day appeared, and
#      Sep 13 is the day whose C1 detections correlate with tide at
#      -0.01 against Sep 11's +0.97.
#    - C1 floor set to 8.0 (rejects 8 of 25 Sep 13 frames, 1 of 25
#      from Sep 12, none of the other 58). C2 left at 0.0/disabled:
#      at the same threshold it rejects nothing at all, which is
#      consistent EITHER with its crop genuinely holding contrast in
#      fog OR with 8.0 being wrong for that crop. No fogged C2 frame
#      has been identified to calibrate against, so no value is set.
#    - This is a PARTIAL fix. The surviving Sep 13 frames sit just
#      above the threshold and are likely degraded too, so expect C1's
#      pooled tide correlation to improve from +0.38 rather than reach
#      Sep 11 levels.
#
#  CHANGES IN V6.8:
#    - TIMEX BIAS CORRECTIONS, and they are now PRODUCT-AWARE. Snap and
#      timex place the waterline differently, so one correction cannot
#      serve both: C1's raw error is -49.8 px on snap but -11.0 px on
#      timex, so cross-applying would inject ~39 px. CameraProfile now
#      carries bias_correction_points (snap) and
#      bias_correction_points_timex, and apply_bias_correction() picks
#      by IMAGE_SUFFIX. If a timex run finds no timex points it warns
#      once and falls back, rather than silently using the wrong set.
#    - C1 timex: CONSTANT -11.03 px (spread/SEM 2.7 -- shape is noise).
#      C2 timex: per-column curve (spread/SEM 12.0 -- shape is real),
#      with the outer two segments clamped to the last trustworthy
#      value because their measured error bars exceeded their values.
#    - Thresholds tuned on station data: NIGHT_BRIGHTNESS_THRESHOLD
#      20 -> 45 (three dusk frames at 13.5/19.0/32.9 were passing while
#      the next-brightest usable frame was 74.0), IMAGE_CONTRAST_FLOOR
#      12 -> 6 (C1 legitimately runs 9.3-38.6 and was being rejected).
#
#  CHANGES IN V6.7:
#    - ABSOLUTE SIGNAL FLOOR. The confidence metric is computed per
#      column as (score - min) / (max - min) and is therefore scale-
#      invariant: it can say which row is best, never whether any row
#      is meaningful. A C2 frame at -1.487 m tide scored 0.99 mean
#      confidence while the right 35% of its energy field peaked at
#      12/255 -- the DP drew a line through noise and it was saved as
#      valid data. Three checks now guard this:
#        * per-column REGIONAL test (peak vs frame peak) -- catches a
#          dead far-field beside a live near-field;
#        * per-column PROMINENCE test (peak vs that column's own
#          median, over envelope rows only) -- catches featureless
#          columns that the regional test passes;
#        * whole-frame IMAGE CONTRAST test -- catches fog/glare/flat
#          scenes. Verified necessary: a synthetic fog frame scored
#          100% signal and 181x prominence on the relative tests and
#          was caught only by this one.
#    - New "Has_Signal" CSV column marks per-column whether the row was
#      chosen with real evidence. Columns with 0 should be EXCLUDED
#      from elevation contours, not weighted equally.
#    - New batch counter "Rejected (no usable signal)".
#
#    TWO THINGS FOUND WRONG DURING IMPLEMENTATION, recorded so they are
#    not re-attempted:
#      * Thresholding the energy magnitude cannot detect whole-frame
#        loss. Every input to the energy (profile, foam, gradient) is
#        already per-frame normalised, so its peak is ~1.0 regardless
#        of image quality. Absolute scale must be measured on the
#        IMAGE, not on anything derived from it.
#      * Column prominence must use the median over ENVELOPE rows only.
#        Including the zeroed region outside the envelope drives the
#        median to ~0 and the ratio to ~1e8, which silently passes
#        everything -- including a featureless grey frame.
#
#  CHANGES IN V6.6:
#    - OPTIONAL CLI OVERRIDES (--image-suffix, --source-dir,
#      --input-dir, --output-dir, --debug-dir, --no-bias-correction).
#      With no arguments the behavior is unchanged. These exist to run
#      an alternative configuration into SEPARATE folders without
#      touching the working snap-based setup or its tuned parameters.
#    - Motivation: the remaining error is now dominated by frame-to-
#      frame scatter rather than bias (bias is corrected to within a
#      few px; day-to-day variation is ~29px). A "snap" is a single
#      instantaneous frame, so its waterline is wherever the swash
#      happened to be at that millisecond -- which is precisely that
#      scatter. A "timex" is a ~10-minute time exposure whose bright
#      band represents the MEAN swash excursion, and is the standard
#      product for this measurement (Plant & Holman 1997; Aarninkhof
#      et al. 2003). Comparing them fairly requires --no-bias-
#      correction, since the existing corrections were fit against
#      snap imagery and would otherwise confound the comparison.
#
#  CHANGES IN V6.5:
#    - BIAS CORRECTIONS RE-DERIVED from real ground truth on current
#      imagery (17 frames C1, 6 frames C2), and -- importantly --
#      COMPOSED with the correction already applied, rather than
#      replacing it with the residual. derive_bias_correction.py
#      previously measured error from already-corrected output and
#      emitted that residual as if it were a new correction; pasting
#      it would have discarded the existing correction's contribution.
#    - C1 IS NOW A CONSTANT -49.8px. Its measured raw error varies
#      only 15.6px across the frame against an 11.6px standard error
#      (1.4 SEM = noise), while the offset itself is 4.3 SEM = real.
#      The OLD C1 curve was fit to 8 frames from a single day, and
#      that fitted tilt was itself generating the monotonic "residual
#      trend" we then spent time diagnosing -- the raw error was flat
#      all along. Fitting shape into noise actively degraded accuracy.
#    - C2 KEEPS A PER-COLUMN CURVE, because its shape genuinely is
#      real: +36.7px to -143.0px, a 179.7px spread at 7.4 SEM. The two
#      cameras needed opposite treatments; applying one policy to both
#      would have been wrong for one of them.
#
#  CHANGES IN V6.4:
#    - NIGHT/DARK FRAME REJECTION. Images whose mean brightness falls
#      below NIGHT_BRIGHTNESS_THRESHOLD are now rejected before
#      detection runs, and counted separately in the batch summary.
#      This closes a real data-quality hole: on 2026-09-10 the 02:00
#      UTC (22:00 local) night captures were completely black yet
#      scored 0.936 / 0.952 mean confidence and were SAVED as valid
#      detections, feeding meaningless shorelines into
#      contour_points.csv. The confidence metric structurally cannot
#      catch this -- it is normalised within each column, so it
#      measures "is the chosen row better than others in its column",
#      not "is there any real signal here", and a black frame passes
#      easily. An explicit brightness gate is the right check.
#    - process_image() now returns a status string ("saved",
#      "low_confidence", "dark", "failed") rather than a bool, so the
#      batch summary can distinguish deliberate dark rejections from
#      confidence discards and from genuine errors.
#
#  CHANGES IN V6.3:
#    - C2 ENVELOPE REVERTED to its previous (Jul 24, 8-frame) values.
#      The v6.2 widening (--margin-above 250) regressed this camera:
#      RMS on the May 21 frame went 174px -> 225px, with mean signed
#      error reaching -391px in the far segment, because the ceiling
#      opened from row 842 to 608 and the detector drifted up to it
#      while the true line sat at the floor (row ~999).
#    - C1 ENVELOPE KEPT WIDENED. The identical change HELPED C1
#      (RMS 15px on the same frame, vs a ~24-34px prior baseline).
#      The cameras' geometry differs enough that envelope tuning must
#      be done independently per camera, not applied symmetrically.
#
#  CHANGES IN V6.2:
#    - BAND_MEAN EDGE BUG FIXED: band_mean() returned 0.0 when a band
#      fell entirely off-frame, which is the MAXIMUM penalty and the
#      exact opposite of its own documented intent ("bands near the
#      top/bottom edge aren't unfairly penalized"). For the water zone
#      (-40,-15) this hard-zeroed the "open water above" evidence for
#      every candidate row within 15px of the crop top, discounting the
#      profile score 3.3x and the final energy (profile squared) ~11x,
#      purely for being near the top of the frame. It now falls back to
#      the nearest in-frame row. This was a prerequisite for the
#      envelope change below: widening a ceiling to crop_top is
#      pointless if rows near crop_top carry an 11x energy penalty.
#    - ENVELOPES RE-DERIVED FROM 14 FRAMES (was 8, single day). Both
#      cameras' ceilings open substantially (C1 up to +166px at the
#      left edge, C2 up to +251px) while floors stay unchanged.
#      Addresses the waterline going missing/wrong at the top-left,
#      which the previous ceiling made structurally unreachable.
#
#  CHANGES IN V6.1:
#    - SUB-PIXEL SHORELINE LOCALIZATION: refine_subpixel() fits a
#      parabola to the smoothed energy around each column's integer DP
#      row and refines it to a fractional-pixel estimate, instead of
#      reporting only whole-pixel rows. The full downstream pipeline
#      (horizontal smoothing, temporal blending, bias correction,
#      clipping) now carries float32 precision throughout instead of
#      rounding to int at each stage -- computing sub-pixel precision
#      and then immediately discarding it a few lines later would have
#      been pointless. CSV output gains "Row_Precise" (float) alongside
#      the existing integer "Row" (kept for backward compatibility with
#      every existing consumer script).
#    - PEAK SHARPNESS: new "Peak_Sharpness" column (0-1), derived from
#      the same parabola fit's curvature -- a RELATIVE per-column data
#      quality indicator (sharp, well-defined energy peak = high;
#      flat/ambiguous = low), not a calibrated pixel-uncertainty value.
#      Worth validating against batch_bias_analysis.py once enough
#      ground truth exists, to see whether low-sharpness columns really
#      do correlate with larger error.
#    - IMAGE PRODUCT TYPE documented explicitly (IMAGE_SUFFIX): the
#      detector still defaults to "snap.jpg" (single-instant frames)
#      for continuity with existing tuning, but ImageProducts also
#      contains "timex.jpg" (time-exposure) and other Argus/CoastCam
#      product types. timex is the scientifically standard choice for
#      shoreline/elevation mapping specifically (reduces swash-cycle
#      noise), and switching is a one-line change -- but requires
#      re-deriving the envelope/profile-scale/bias-correction
#      parameters against timex ground truth before trusting results,
#      since all of those were tuned against snap images.
#
#  CHANGES IN V6.0:
#    - VERTICAL PROFILE MATCHER: the core shoreline energy no longer
#      scores a single pixel in isolation. compute_vertical_profile_
#      score() evaluates a whole vertical neighborhood around each
#      candidate row -- open water well above, foam near the row, wet
#      sand just below, dry sand further below -- and scores how well
#      that neighborhood matches the expected shoreline sequence. This
#      is the single highest-leverage change identified in review: a
#      shoreline is a transition region spanning tens of pixels, not
#      one pixel, and scoring in isolation was the main source of
#      false positives from locally shoreline-like single pixels.
#    - WET SAND LIKELIHOOD: new compute_wet_sand_likelihood() fills
#      the missing class between water/foam and dry sand. The
#      foam-to-wet-sand transition is often a cleaner signal than
#      foam-to-dry-sand, and the profile matcher now uses it as its
#      own zone.
#    - Everything else (camera handling, envelope, DP solver,
#      temporal filter, evaluation tools) is UNCHANGED from V5.2, on
#      purpose -- isolating this one change lets its impact be
#      measured directly with the existing ground-truth tools rather
#      than being confounded with other simultaneous changes.
#
#  CHANGES IN V5.2:
#    - CAMERA-SPECIFIC SEARCH ENVELOPE: instead of a single rectangular
#      crop, each camera now has a per-column [min_row, max_row] band
#      built by interpolating between a few hand-set control points.
#      Energy outside this envelope is zeroed before the DP solver,
#      preventing the detector from ever considering physically
#      impossible shoreline locations (e.g. open ocean far offshore).
#    - WATER WEIGHTING REPLACED: water_gate (hard threshold) replaced
#      with a smooth power-law weight (water ** 2.5), so calm offshore
#      water (which still passes a hard threshold easily) contributes
#      much less than genuinely strong, shoreline-adjacent water
#      signal, instead of being treated as equally valid.
#    - HORIZONTAL CONTINUITY BLENDING: energy is now blended with its
#      immediate left/right neighboring columns (0.6/0.2/0.2) before
#      the DP solver runs. This reinforces genuine shoreline signal
#      that persists across columns and suppresses isolated anomalous
#      spikes in a single column.
#    - Reverted to the simpler single-state DP solver (dropped the
#      multi-state curvature-tracking version), since the envelope +
#      horizontal blending should address the "jump and stay" failure
#      more directly, without the ~11x slowdown from tracking per-step
#      state.
#
#  Requirements
#  ------------
#      Python 3.8+
#      OpenCV
#      NumPy
#
# =====================================================================

import os
import re
import cv2
import warnings
import argparse
import csv
import time
import shutil

import numpy as np

from pathlib import Path
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone


VERSION = "6.9"
PROGRAM_NAME = "USGS Waterline Detector"

SOURCE_FOLDER = "/mnt/I2Rgus_Data/ImageProducts"
INPUT_FOLDER = "/mnt/I2Rgus_Data/waterline/input"
OUTPUT_FOLDER = "/mnt/I2Rgus_Data/waterline/processed"
DEBUG_FOLDER = "/mnt/I2Rgus_Data/waterline/debug"

# WHICH IMAGE PRODUCT TO DETECT ON -- this matters more than it looks.
# ImageProducts contains several Argus/CoastCam product types per
# capture, not just one:
#   snap.jpg   a single instantaneous video frame. The shoreline in a
#              snap is wherever the swash happened to be at that exact
#              moment -- noisy, and biased by whichever wave phase was
#              captured. THIS IS THE CURRENT DEFAULT, for continuity
#              with everything already tuned (envelope, profile-matcher
#              scaling, bias correction) against snap images.
#   timex.jpg  a time-exposure (typically ~10 min average). The bright
#              band directly represents the MEAN swash excursion zone
#              over that window -- this is the actual basis of the
#              published video-shoreline-elevation method (Plant &
#              Holman 1997; Aarninkhof et al. 2003), and is the
#              scientifically preferred choice for elevation mapping
#              specifically. Switching to it is a one-line change
#              below, but the envelope/profile-matcher/bias-correction
#              parameters were all tuned against snap images and would
#              need to be re-derived against timex ground truth before
#              trusting results at the same level -- don't switch this
#              without also re-running the full ground-truth validation
#              loop (collect_all_ground_truth.py -> batch_bias_analysis.py)
#              against timex images specifically.
#   bright.jpg / dark.jpg / var.jpg  brightest/darkest/variance-of-pixel
#              images over the same averaging window. var.jpg in
#              particular highlights the swash zone as a band of high
#              variance (constantly wetting/drying), which is sometimes
#              used as a complementary/alternative shoreline indicator
#              -- not currently used by this detector at all.
IMAGE_SUFFIX = "snap.jpg"

IMAGES_PER_RUN = None  # None = process every matching image in SOURCE_FOLDER, no cap

# Any processed image whose mean confidence falls below this is
# considered unreliable and is NOT saved to OUTPUT_FOLDER (no CSV, no
# overlay PNG) -- tracked separately from failures/successes in the
# batch summary so it's clear these were deliberately discarded, not
# errors.
CONFIDENCE_DISCARD_THRESHOLD = 0.9

# Mean grayscale brightness (0-255) below which an image is treated as
# a night/dark frame containing no usable waterline, and is rejected
# BEFORE detection runs.
#
# This exists because the confidence metric cannot catch this case. On
# 2026-09-10 the 02:00 UTC (22:00 local) night captures -- completely
# black frames -- were scored at 0.936 and 0.952 mean confidence and
# saved as valid detections. Confidence is normalised within each
# column, so it stays high even when a column carries no signal at
# all; it measures "how much better is the chosen row than the rest of
# its column", not "is there anything here". A black frame therefore
# passes the confidence gate easily while producing a meaningless
# shoreline that then flows into contour_points.csv as if it were real
# data. An explicit brightness gate is the correct check.
NIGHT_BRIGHTNESS_THRESHOLD = 45.0

# ABSOLUTE SIGNAL FLOOR ------------------------------------------------
# The confidence metric cannot detect the ABSENCE of signal. It is
# computed per column as (score - col_min) / (col_max - col_min), so it
# is scale-invariant by construction: a column containing nothing but
# noise still yields a high value, because some row is always nominally
# "best". A C2 frame at -1.487 m tide scored 0.99 mean confidence while
# the entire right 35% of its energy field peaked at 12/255 -- the DP
# was drawing a line through noise, and the result was saved as valid.
# The night-frame brightness gate closed one specific instance of this;
# these two checks address the general case.
#
# (1) REGIONAL loss. build_shoreline_energy() ends with
#     normalize_float(), which rescales each frame to [0,1] using that
#     frame's own min/max. Absolute scale is therefore gone, but
#     RELATIVE structure within the frame survives -- which is why the
#     dead right-hand side reads near zero while the frame max reads 1.
#     A column whose peak falls below this fraction of the frame peak is
#     treated as carrying no usable signal.
COLUMN_SIGNAL_FLOOR_FRACTION = 0.05

#     The fraction test above compares each column to the FRAME peak,
#     so it detects regional loss -- a dead far-field beside a live
#     near-field. It does NOT detect uniform degradation: pure noise
#     has no regional contrast, so every column passes. (Verified: a
#     synthetic all-noise field scored 100% by that test alone.)
#
#     This second test asks whether a column's peak stands out from
#     that column's OWN distribution. A real energy ridge sits far
#     above the column median; noise does not. Both tests must pass.
COLUMN_PEAK_PROMINENCE_MIN = 5.0

#     If fewer than this fraction of columns carry signal, the frame as
#     a whole is rejected rather than saved with a partly-fabricated
#     line.
MIN_SIGNAL_COLUMN_FRACTION = 0.60

# (2) WHOLE-FRAME loss. Uniform degradation (fog, glare, a dirty lens)
#     weakens every column together, so the relative tests above see
#     nothing wrong.
#
#     NOTE ON AN APPROACH THAT DOES NOT WORK: measuring the energy
#     field's magnitude is useless here, because every input to it
#     (profile, foam, gradient) is already per-frame normalised, so the
#     energy peak is ~1.0 by construction regardless of image quality.
#     Absolute scale is destroyed long before the final normalisation.
#     Detecting whole-frame loss requires measuring the IMAGE, not the
#     energy derived from it.
#
#     Standard deviation of the cropped greyscale is such a measure: it
#     is absolute, cheap, and low for fog/flat/washed-out scenes while
#     high for a beach with water, foam and sand. A featureless grey
#     test frame scores ~3; a normal beach frame scores ~50+.
IMAGE_CONTRAST_FLOOR = 0.0

SAVE_DEBUG_IMAGES = True
EXPORT_CSV = True
VERBOSE = True


# =====================================================================
# CONFIG
# =====================================================================

@dataclass
class DetectorConfig:
    max_vertical_step: int = 15
    vertical_penalty: float = 0.35
    smoothing_kernel: int = 15
    blur_sizes: tuple = (3, 5, 9)

    # Water weighting: replaced hard threshold with a smooth power-law
    # curve. Calm offshore water still scores highly on raw
    # "water_likelihood", but this makes only VERY water-like pixels
    # contribute strongly, sharply discounting marginal/ambiguous ones.
    water_weight_power: float = 2.5

    sand_gate_threshold: float = 0.45
    transition_power: float = 1.5
    gradient_power: float = 1.2

    # Horizontal continuity blending weights (must sum to 1.0).
    horiz_center_weight: float = 0.9
    horiz_neighbor_weight: float = 0.05

    confidence_threshold: float = 0.55


CONFIG = DetectorConfig()


@dataclass
class DetectorState:
    # Keyed by camera profile name (e.g. "CACO05_C1" / "CACO05_C2"),
    # NOT a single global previous shoreline. Images are processed in
    # alphabetical filename order (see discover_source_images), which
    # interleaves cameras (t1.c1, t1.c2, t2.c1, t2.c2, ...) since "c1"
    # sorts before "c2" for the same timestamp. A single shared
    # previous_shoreline meant temporal_filter was blending each
    # frame with the PREVIOUS CAMERA'S completely different geometry
    # 30% of the time it ran -- both cameras happen to crop to the
    # same width (2448px), so the length-mismatch guard never caught
    # it. This was a real bug, not a tuning issue: it explains the
    # smooth, opposite-signed, far-column-growing error seen in
    # bias_analysis.py on both cameras.
    previous_shoreline: dict = field(default_factory=dict)
    previous_confidence: dict = field(default_factory=dict)
    processed_images: int = 0
    warned_timex_bias: bool = False


STATE = DetectorState()


def log(message):
    if VERBOSE:
        print(message, flush=True)


def normalize(image):
    image = image.astype(np.float32)
    minimum, maximum = np.min(image), np.max(image)
    if maximum - minimum < 1e-8:
        return np.zeros(image.shape, dtype=np.uint8)
    image = (image - minimum) / (maximum - minimum)
    return (image * 255.0).astype(np.uint8)


def normalize_float(image):
    image = image.astype(np.float32)
    minimum, maximum = image.min(), image.max()
    if maximum - minimum < 1e-8:
        return np.zeros_like(image)
    return (image - minimum) / (maximum - minimum)


# =====================================================================
# CAMERA PROFILES + SEARCH ENVELOPES
# =====================================================================

def make_envelope_points(anchors, num_points=25):
    """
    Takes a small set of hand-set anchors
    [(x_fraction, min_row_fraction, max_row_fraction), ...]
    and densifies them into `num_points` evenly-spaced control points
    by linear interpolation. This is equivalent to interpolating from
    the anchors directly, but gives build_envelope_mask a dense set of
    points to work with (per review feedback: prefer 20-30 points over
    3-5), which matters once individual points are hand-tuned against
    ground truth and no longer fall on a straight line between anchors.
    """
    anchors = sorted(anchors, key=lambda p: p[0])
    xs = np.array([a[0] for a in anchors])
    min_fracs = np.array([a[1] for a in anchors])
    max_fracs = np.array([a[2] for a in anchors])

    dense_x = np.linspace(0.0, 1.0, num_points)
    dense_min = np.interp(dense_x, xs, min_fracs)
    dense_max = np.interp(dense_x, xs, max_fracs)
    return tuple(zip(dense_x.tolist(), dense_min.tolist(), dense_max.tolist()))


@dataclass
class CameraProfile:
    name: str
    station: str
    camera: str

    # Outer crop, same as before -- coarse bounding box applied first.
    crop_top: float
    crop_bottom: float
    crop_left: float
    crop_right: float

    # NEW: per-column search envelope, defined as control points
    # (x_fraction, min_row_fraction, max_row_fraction) within the
    # CROPPED image. Interpolated across columns. Energy outside this
    # band is zeroed before the DP solver runs, so the detector can
    # never consider physically impossible shoreline locations (e.g.
    # open ocean far offshore, or dry sand/dune).
    #
    # These are STARTING ESTIMATES based on the one ground-truth frame
    # collected so far -- they should be refined with more ground truth
    # across different tide states before being trusted as final.
    envelope_points: tuple = field(default_factory=tuple)

    max_vertical_step: int = 5
    vertical_penalty: float = 0.35
    tracking_window: int = 40

    # Per-column adaptive DP step: perspective means the shoreline can
    # jump many pixels between adjacent columns near the camera, but
    # only a few pixels far from the camera. step_near/step_far define
    # the allowed range, linearly interpolated across columns (column
    # 0 = near camera/left edge, last column = far/right edge). If
    # step_far is None, falls back to the flat max_vertical_step.
    step_near: int = None
    step_far: int = None

    # Per-column scaling of the vertical profile matcher's band
    # widths (see compute_vertical_profile_score). Same motivation as
    # step_near/step_far: perspective compresses real-world vertical
    # extent into fewer pixels toward the far edge of an oblique
    # camera, so a fixed-pixel-height zone that's the right physical
    # scale near the camera can span multiple real classes at once
    # far from it. 1.0 = no scaling (band widths as written). Falls
    # back to no scaling if either is unset.
    profile_scale_near: float = None
    profile_scale_far: float = None

    # Per-column signed-error correction, in raw pixels: (x_fraction,
    # correction_px) points, interpolated across columns and
    # SUBTRACTED from the raw detector row. Derived empirically from
    # batch_bias_analysis.py's per-camera segment-wise mean signed
    # error, averaged across multiple ground-truth frames (validated
    # as a real, consistent effect, not single-frame noise, before
    # being added here -- see derive_bias_correction.py). None/empty
    # means no correction applied (default, unchanged behavior).
    bias_correction_points: tuple = field(default_factory=tuple)

    # Snap and timex put the waterline in DIFFERENT places, so they
    # need different corrections -- measured, not assumed: C1's raw
    # error is -49.8 px on snap but only -11.0 px on timex. Applying
    # one product's correction to the other would inject roughly 39 px
    # of error. apply_bias_correction() selects between these on
    # IMAGE_SUFFIX, so a snap run and a timex run each get their own.
    bias_correction_points_timex: tuple = field(default_factory=tuple)

    # Minimum standard deviation of the CROPPED band for a frame to be
    # used. 0.0 disables.
    #
    # This exists because whole-image contrast does not detect fog.
    # Measured on 108 archived frames: a fogged C1 frame scored global
    # std 31.8 -- higher than many clear frames -- while its cropped
    # band scored 5.3. Fog flattens the beach region specifically,
    # while sky and dune keep enough brightness variation to inflate
    # the global figure. Measuring where the detector actually looks
    # separates them; measuring the whole frame does not.
    #
    # The value is PER CAMERA and does not transfer: each crop covers a
    # different region with different content. C1 is set from data
    # (see below); C2 is left disabled because no fogged C2 frame has
    # yet been identified to calibrate against, and guessing a
    # threshold risks discarding good data.
    band_contrast_floor: float = 0.0


CAMERAS = {
    "CACO05_C1": CameraProfile(
        "CACO05_C1", "CACO05", "C1",
        crop_top=0.26, crop_bottom=0.86, crop_left=0.00, crop_right=1.00,
        max_vertical_step=15,   # fallback if step_near/step_far unset
        step_near=18,
        step_far=4,
        # First-pass hypothesis based on bias_analysis.py's per-segment
        # breakdown (error grew from -33px at col 0 to +183px at the
        # far edge on 2026-07-24 test frame): shrink profile-matcher
        # band widths toward the far column, leave near-camera side
        # unchanged since it already showed the smallest error.
        # VALIDATE by re-running batch_bias_analysis.py after this
        # change -- if the segment-wise error flattens out, this was
        # the right fix; if not, revert and look elsewhere.
        profile_scale_near=1.0,
        profile_scale_far=0.4,
        # Bias correction, derived from batch_bias_analysis.py's
        # per-camera segment-wise mean signed error averaged across 8
        # ground-truth frames (2026-07-24). NOTE: this was temporarily
        # DISABLED over a concern that a single-day fit wouldn't
        # generalize -- but an A/B test disabling it on July 27 data
        # (a different day) showed RMS getting WORSE (34.06px ->
        # 37.93px) and the segment-wise error roughly DOUBLING in
        # magnitude with the same shape. So despite being derived from
        # only one (now-unrecoverable) day, this correction empirically
        # still helps out-of-sample. Re-enabled on that direct evidence
        # rather than the theoretical generalization concern. Still
        # worth replacing with a derive_bias_correction.py fit against
        # the growing archive (see archive_run.py) once enough varied
        # days accumulate there -- but "no correction" is not the safer
        # default it was assumed to be.
        bias_correction_points=(
            # CONSTANT correction, derived from 17 ground-truth frames
            # (Nov 20 + May 21 + Sep 10) via derive_bias_correction.py.
            #
            # WHY CONSTANT AND NOT A CURVE: the measured raw error for
            # this camera averages -49.8px and varies only 15.6px
            # across the frame, while the standard error of the mean at
            # each point is 11.6px. The apparent shape is therefore
            # 1.4 SEM -- indistinguishable from noise. The previous
            # correction WAS a curve (fit to 8 frames from one day) and
            # that fitted tilt was itself creating the monotonic
            # residual we then measured: raw error was flat, the
            # correction added a 35px slope, and the "trend" we chased
            # was an artifact of our own correction.
            #
            # The offset, by contrast, is real at 4.3 SEM. Removing
            # just the mean should take RMS from ~69px to ~48px.
            #
            # Do not replace this with a per-column curve unless
            # derive_bias_correction.py reports spread/SEM > 3.
            (0.00, -49.80),
            (1.00, -49.80),
        ),
        # Set from 108 archived frames across Sep 11-15. The eight
        # lowest band_std values were ALL from Sep 13 (5.3, 5.3, 6.0,
        # 6.1, 6.6, 7.3, 7.4, 7.7) before any other day appeared, and
        # Sep 13 is the day whose detections correlate with tide at
        # -0.01 while Sep 11 manages +0.97. A threshold of 8.0 rejects
        # 8 of 25 Sep 13 frames, 1 of 25 from Sep 12, and none of the
        # other 58.
        #
        # PARTIAL FIX, not a complete one: the remaining Sep 13 frames
        # sit just above the threshold and are probably degraded too.
        # Expect this to move C1's pooled correlation up from +0.38,
        # not all the way to Sep 11 levels.
        band_contrast_floor=8.0,
        bias_correction_points_timex=(
            # TIMEX correction, from 12 ground-truth frames (Sep 11)
            # via derive_bias_correction.py --suffix .timex. --uncorrected.
            #
            # CONSTANT, not a curve: measured spread across the frame is
            # 11.2 px against a 4.2 px standard error -- spread/SEM 2.7,
            # below the threshold of 3, so the apparent shape is not
            # distinguishable from noise. The offset itself is real at
            # 2.6 SEM. Fitting the wiggle would repeat the snap mistake,
            # where a curve fitted to one day encoded a 35 px tilt that
            # was not in the data.
            #
            # Note how much smaller this is than the snap correction
            # (-49.8 px): the timex time-average removes most of the
            # swash-phase offset rather than requiring it to be
            # corrected away.
            #
            # ONE DAY of ground truth. Re-derive as the cron accumulates
            # more; every single-day fit in this project has generalised
            # poorly.
            (0.00, -11.03),
            (1.00, -11.03),
        ),
        envelope_points=make_envelope_points((
            # Derived from 14 ground-truth frames (May 21 + Jul 24 +
            # Jul 27) via:
            #   derive_envelope.py --margin-above 250 --margin-below 60
            #                      --image-height 2048 c1
            # Replaces an earlier 8-frame/single-day (Jul 24) envelope
            # whose ceiling at the left edge sat at full-image row 698.
            # That ceiling made any waterline higher in the frame than
            # row 698 structurally unreachable by the DP solver -- the
            # suspected cause of the waterline going missing/wrong at
            # the top-left. Ceiling now opens up to 166px higher (to
            # crop_top); floors are unchanged from the previous
            # envelope, so no previously-working condition is excluded.
            #
            # TRADEOFF: min_row_fraction is 0.000 at every column, so
            # the envelope no longer constrains the TOP at all for this
            # camera -- only crop_top does. That removes a guardrail
            # against the detector wandering upward into open water. If
            # C1's RMS degrades vs the ~24-34px baseline, or the line
            # drifts offshore, re-derive with a tighter --margin-above
            # (e.g. 120) rather than keeping 250.
            (0.000, 0.000, 0.278),
            (0.200, 0.000, 0.236),
            (0.400, 0.000, 0.207),
            (0.600, 0.000, 0.182),
            (0.800, 0.000, 0.158),
            (1.000, 0.000, 0.154),
        )),
    ),
    "CACO05_C2": CameraProfile(
        "CACO05_C2", "CACO05", "C2",
        crop_top=0.03, crop_bottom=0.50, crop_left=0.00, crop_right=1.00,
        max_vertical_step=15,   # fallback if step_near/step_far unset
        step_near=18,
        step_far=4,
        # Same hypothesis as C1 (see comment there): error also grew
        # sharply toward the far column on this camera (+7px -> -280px
        # across segments on the same test frame). Same starting scale
        # applied; validate independently since C2's geometry/sweep is
        # very different from C1's.
        profile_scale_near=1.0,
        profile_scale_far=0.4,
        # Bias correction, re-enabled -- see C1 comment for the full
        # rationale (disabling it was tested and made results worse,
        # not better, on out-of-sample July 27 data).
        bias_correction_points=(
            # PER-COLUMN correction, derived from 6 ground-truth frames
            # (May 21 + Sep 10) via derive_bias_correction.py.
            #
            # Unlike C1, this camera's shape IS real: the raw error
            # sweeps +36.7px to -143.0px across the frame, a 179.7px
            # spread against a 24.2px standard error -- 7.4 SEM. A
            # constant correction would leave almost all of that
            # uncorrected.
            #
            # CAVEATS, both material:
            #  * n=6, because ~70% of this camera's daylight frames are
            #    discarded by the confidence gate and so produce no
            #    output to compare against. The 6 scored frames are
            #    therefore a biased sample of C2's EASIER images.
            #  * The per-point std grows from 38px at the near edge to
            #    119px at the far edge, so the large far-field values
            #    (-103, -143) are the least certain part of the curve
            #    despite being the largest.
            # Re-derive once C2's discard rate is understood and more
            # of its frames can be scored.
            (0.050, 36.74),
            (0.150, 22.96),
            (0.249, 14.89),
            (0.349, 4.86),
            (0.449, -13.35),
            (0.548, -24.42),
            (0.648, -39.80),
            (0.748, -68.77),
            (0.848, -102.88),
            (0.949, -143.00),
        ),
        # Disabled deliberately. At C1's threshold of 8.0 this camera
        # rejects ZERO frames across all 108, including on Sep 13 when
        # C2 was also hazy -- so either C2's crop genuinely holds more
        # contrast in fog (more water and dune texture in view), or 8.0
        # is simply not calibrated for this crop. Both are consistent
        # with the data, so no value is set until a fogged C2 frame
        # exists to calibrate against.
        band_contrast_floor=6.0,
        bias_correction_points_timex=(
            # TIMEX correction, from 13 ground-truth frames (Sep 11)
            # via derive_bias_correction.py --suffix .timex. --uncorrected.
            #
            # PER-COLUMN here, unlike C1: measured spread is 61.6 px
            # against a 5.1 px standard error -- spread/SEM 12.0, so the
            # shape is unambiguous. The two cameras genuinely differ and
            # need different treatments, as they did on snap.
            #
            # THE LAST TWO SEGMENTS ARE CLAMPED, not measured. Their
            # measured values were -41.28 +/- 48.26 and -57.78 +/- 67.00
            # -- error bars LARGER than the values. That is the far
            # field, where the energy peak collapses ~15x (0.1168 ->
            # 0.0076), where the detector was shown to hold a flat line
            # 53 px RMS from hand-traced truth regardless of
            # profile_scale_far, step_far or vertical_penalty, and where
            # the Has_Signal filter removes most columns from the
            # contours anyway. Fitting a correction to those rows would
            # be fitting to output already known to be wrong, so they
            # instead inherit the last trustworthy value (-17.22).
            #
            # This clamp is a judgement call, not a measurement. If the
            # far field is ever made reliable, re-derive rather than
            # assuming the correction stays flat out there.
            (0.050,   3.87),
            (0.150,  -1.19),
            (0.249,  -8.15),
            (0.349,  -8.49),
            (0.449,  -5.26),
            (0.548,  -4.55),
            (0.648,  -8.44),
            (0.748, -17.22),
            (0.848, -17.22),
            (0.949, -17.22),
        ),
        envelope_points=make_envelope_points((
            # REVERTED to the 8-frame (Jul 24) envelope after a widened
            # 14-frame version measurably REGRESSED this camera.
            #
            # Evidence: with --margin-above 250 the ceiling at x=1.0
            # opened from full-image row 842 to 608. On the May 21
            # frame the true waterline there sits at row ~999 (the
            # envelope FLOOR), and the detector drifted to the ceiling
            # -- mean signed error -391px in the far segment, RMS for
            # the frame going 174px -> 225px. The old 842 ceiling
            # capped that same drift at 157px. C1's identical widening
            # HELPED (RMS 15px), so this is a C2-specific outcome, not
            # a general one -- the two cameras' geometry differs enough
            # that they need independent envelope tuning.
            #
            # Do NOT re-widen this camera without first collecting
            # ground truth on CURRENT imagery (see note below) --
            # widening was validated on one out-of-season frame, which
            # is not a sufficient basis either way.
            (0.00, 0.046, 0.191),
            (0.20, 0.158, 0.306),
            (0.40, 0.352, 0.509),
            (0.60, 0.511, 0.666),
            (0.80, 0.658, 0.826),
            (1.00, 0.812, 0.975),
        )),
    ),
}

def identify_camera(filename):
    """
    Identifies camera from filename by looking for ".C1." / ".C2."
    (dot-delimited, so e.g. "...GMT.2026.CACO05.c1.snap.jpg" matches
    but a bare "c1.snap.jpg" -- no leading dot before "c1" -- does
    NOT). Previously this silently fell back to CACO05_C1 for
    anything that didn't match either pattern, which meant oddly-named
    files (e.g. "c2.snap.jpg", missing its timestamp prefix) got
    silently processed with the WRONG camera's crop/envelope/bias
    correction -- confirmed as the cause of a 256px-RMS outlier in a
    real batch (a "c2" file processed as if it were C1). Raising here
    instead surfaces the problem immediately (process_image already
    catches exceptions per-image and reports them, so this just turns
    a silent wrong answer into a visible, skippable error).
    """
    name = Path(filename).name.upper()
    if ".C1." in name:
        return CAMERAS["CACO05_C1"]
    if ".C2." in name:
        return CAMERAS["CACO05_C2"]
    raise ValueError(
        f"Could not identify camera from filename '{filename}' -- "
        f"expected '.c1.' or '.c2.' somewhere in the name (dot-delimited). "
        f"Refusing to guess, since guessing wrong silently corrupts results."
    )


def crop_to_camera(image, profile):
    h, w = image.shape[:2]
    top = int(profile.crop_top * h)
    bottom = int(profile.crop_bottom * h)
    left = int(profile.crop_left * w)
    right = int(profile.crop_right * w)
    return image[top:bottom, left:right], top, left


def build_envelope_mask(height, width, profile):
    """
    Builds a boolean mask (True = allowed) over the CROPPED image,
    based on interpolating profile.envelope_points across columns.
    Each control point is (x_fraction, min_row_fraction, max_row_fraction).
    """
    mask = np.zeros((height, width), dtype=bool)

    xs = np.array([p[0] for p in profile.envelope_points])
    min_fracs = np.array([p[1] for p in profile.envelope_points])
    max_fracs = np.array([p[2] for p in profile.envelope_points])

    col_x_fracs = np.linspace(0.0, 1.0, width)
    min_row_frac = np.interp(col_x_fracs, xs, min_fracs)
    max_row_frac = np.interp(col_x_fracs, xs, max_fracs)

    min_rows = np.clip((min_row_frac * height).astype(np.int32), 0, height - 1)
    max_rows = np.clip((max_row_frac * height).astype(np.int32), 0, height - 1)

    for col in range(width):
        lo = min(min_rows[col], max_rows[col])
        hi = max(min_rows[col], max_rows[col])
        mask[lo:hi + 1, col] = True

    return mask


# =====================================================================
# FILE MANAGEMENT
# =====================================================================

def ensure_directory(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def initialize_directories():
    ensure_directory(INPUT_FOLDER)
    ensure_directory(OUTPUT_FOLDER)
    ensure_directory(DEBUG_FOLDER)


def clear_directory(directory):
    directory = Path(directory)
    if not directory.exists():
        return
    removed = 0
    for item in directory.iterdir():
        if item.is_file():
            item.unlink()
            removed += 1
    log(f"Cleared {removed} file(s) from {directory}")


def extract_epoch_from_filename(name):
    """
    Filenames in this project start with a Unix epoch timestamp, e.g.
    "1763666100.Thu.Nov.20_19_15_00.GMT.2025.CACO-05.c1.snap.jpg".
    Returns None (rather than guessing) if no leading epoch is found.
    """
    match = re.match(r"^(\d{9,10})\.", name)
    if not match:
        return None
    return int(match.group(1))


# If a source image's actual filesystem modification time differs from
# what its filename claims by more than this many days, flag it as a
# mismatch. Set well above normal upload/transfer slop (which is
# seconds, not days) so this only fires on genuinely suspicious cases
# -- e.g. a file dated November 2025 that's still sitting in the
# source folder in August 2026, which happened for real on this
# project and could only be confirmed after the fact via a copy that
# happened to preserve the original mtime. Logging this at discovery
# time means that evidence is captured immediately, before any nightly
# wipe or cleanup can remove it.
TIMESTAMP_MISMATCH_WARNING_DAYS = 1


def check_timestamp_mismatch(image_path):
    """
    Compares the epoch embedded in a source image's filename to its
    actual filesystem modification time. Returns None if they agree
    (within TIMESTAMP_MISMATCH_WARNING_DAYS) or if the filename has no
    parseable epoch; otherwise returns a human-readable warning string
    describing the mismatch, for the caller to log and/or count.
    """
    claimed_epoch = extract_epoch_from_filename(image_path.name)
    if claimed_epoch is None:
        return None

    try:
        actual_mtime = image_path.stat().st_mtime
    except OSError:
        return None

    gap_days = abs(actual_mtime - claimed_epoch) / 86400.0
    if gap_days <= TIMESTAMP_MISMATCH_WARNING_DAYS:
        return None

    claimed_dt = datetime.fromtimestamp(claimed_epoch, tz=timezone.utc)
    actual_dt = datetime.fromtimestamp(actual_mtime, tz=timezone.utc)
    return (f"TIMESTAMP MISMATCH: {image_path.name} -- filename claims "
            f"{claimed_dt.isoformat()}, but filesystem modification time is "
            f"{actual_dt.isoformat()} ({gap_days:.1f} days apart). This file "
            f"may be stale/mislabeled -- worth checking with whoever manages "
            f"the camera upload/sync process.")


def discover_source_images(limit=None):
    images = sorted(Path(SOURCE_FOLDER).rglob(f"*{IMAGE_SUFFIX}"))
    log(f"Discovered {len(images)} source images.")

    mismatches = []
    for image in images:
        warning = check_timestamp_mismatch(image)
        if warning is not None:
            log(warning)
            mismatches.append(image.name)
    if mismatches:
        log(f"WARNING: {len(mismatches)} source image(s) have a filename/filesystem "
            f"timestamp mismatch of more than {TIMESTAMP_MISMATCH_WARNING_DAYS} day(s) "
            f"-- see above for details. This does not block processing, but the "
            f"resulting shoreline/date association for those files may not be trustworthy.")

    if limit is not None:
        images = images[:limit]
        log(f"Limiting to first {len(images)} images (IMAGES_PER_RUN={limit}).")
    else:
        log(f"Processing all {len(images)} discovered images (no IMAGES_PER_RUN cap).")
    return images


def copy_input_images(limit=IMAGES_PER_RUN):
    initialize_directories()
    images = discover_source_images(limit=limit)
    copied = 0
    for image in images:
        destination = Path(INPUT_FOLDER) / image.name
        shutil.copy2(image, destination)
        copied += 1
    log(f"Copied {copied} images into {INPUT_FOLDER}.")
    return copied


def list_input_images():
    return sorted(Path(INPUT_FOLDER).glob("*.jpg"))


def load_image(filename):
    image = cv2.imread(str(filename))
    if image is None:
        raise RuntimeError(f"Could not read image: {filename}")
    return image


def debug_directory(filename):
    directory = Path(DEBUG_FOLDER) / Path(filename).stem
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def save_image(filename, image):
    cv2.imwrite(str(filename), image)


# =====================================================================
# IMAGE PREPROCESSING
# =====================================================================

@dataclass
class ImageData:
    filename: str = ""
    original: np.ndarray = None
    cropped: np.ndarray = None
    normalized: np.ndarray = None
    gray: np.ndarray = None
    hsv: np.ndarray = None
    lab: np.ndarray = None
    red: np.ndarray = None
    green: np.ndarray = None
    blue: np.ndarray = None
    saturation: np.ndarray = None
    value: np.ndarray = None
    lab_b: np.ndarray = None
    crop_top: int = 0
    crop_left: int = 0
    profile: CameraProfile = None


def normalize_color(image):
    channels = cv2.split(image)
    return cv2.merge([normalize(c) for c in channels])


def correct_illumination(image):
    image = image.astype(np.float32)
    background = cv2.GaussianBlur(image, (0, 0), sigmaX=50, sigmaY=50)
    corrected = np.clip(image - background + 128.0, 0, 255)
    return corrected.astype(np.uint8)


def prepare_color_spaces(data):
    image = data.normalized
    data.gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    data.hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
    data.lab = cv2.cvtColor(image, cv2.COLOR_BGR2LAB)
    data.blue = image[:, :, 0]
    data.green = image[:, :, 1]
    data.red = image[:, :, 2]
    data.saturation = data.hsv[:, :, 1]
    data.value = data.hsv[:, :, 2]
    data.lab_b = data.lab[:, :, 2]


def prepare_image(filename):
    profile = identify_camera(filename)
    image = load_image(filename)
    cropped, top, left = crop_to_camera(image, profile)
    corrected = correct_illumination(cropped)
    normalized = normalize_color(corrected)

    data = ImageData()
    data.filename = str(filename)
    data.original = image
    data.cropped = cropped
    data.normalized = normalized
    data.crop_top = top
    data.crop_left = left
    data.profile = profile
    prepare_color_spaces(data)
    return data


# =====================================================================
# PHYSICS-BASED LIKELIHOOD MODELS
# =====================================================================

def compute_water_likelihood(data):
    blue = data.blue.astype(np.float32)
    red = data.red.astype(np.float32)
    green = data.green.astype(np.float32)
    blue_ratio = normalize_float(blue / (red + green + 1.0))

    lab_b = data.lab_b.astype(np.float32)
    coolness = normalize_float(255.0 - lab_b)

    sat = normalize_float(data.saturation.astype(np.float32))

    gray = data.gray.astype(np.float32)
    mean = cv2.blur(gray, (9, 9))
    texture = normalize_float(np.abs(gray - mean))
    smoothness = 1.0 - texture

    water = (0.40 * blue_ratio + 0.30 * coolness + 0.15 * sat + 0.15 * smoothness)
    return normalize_float(water)


def compute_dry_sand_likelihood(data):
    brightness = normalize_float(data.value.astype(np.float32))
    lab_b = data.lab_b.astype(np.float32)
    warmth = normalize_float(lab_b)
    dryness = normalize_float(255.0 - data.saturation.astype(np.float32))

    gray = data.gray.astype(np.float32)
    mean = cv2.blur(gray, (9, 9))
    mean2 = cv2.blur(gray * gray, (9, 9))
    variance = np.maximum(mean2 - mean * mean, 0.0)
    texture = normalize_float(np.sqrt(variance))

    dry_sand = (0.40 * brightness + 0.25 * warmth + 0.20 * dryness + 0.15 * texture)
    return normalize_float(dry_sand)


def compute_wet_sand_likelihood(data):
    """
    Wet sand sits between water and dry sand: darker than dry sand
    (it reflects sky/water and stays saturated with moisture), but
    less saturated/blue than open water, with a distinctive muted,
    slightly warm tone and a smoother (less granular) texture than
    dry sand. This is a first-pass heuristic -- worth validating
    against ground truth alongside the other likelihood maps, same as
    water/dry_sand/foam.
    """
    brightness = normalize_float(data.value.astype(np.float32))
    sat = normalize_float(data.saturation.astype(np.float32))
    lab_b = data.lab_b.astype(np.float32)
    warmth = normalize_float(lab_b)

    darkness = 1.0 - brightness

    # Peaks at moderate saturation (not washed-out dry sand, not
    # richly saturated open water); falls off toward either extreme.
    moderate_sat = np.clip(1.0 - np.abs(sat - 0.40) * 2.5, 0.0, 1.0)

    gray = data.gray.astype(np.float32)
    smoothness = 1.0 - normalize_float(compute_local_std(gray, ksize=7))

    wet_sand = (0.35 * darkness + 0.30 * moderate_sat + 0.20 * warmth + 0.15 * smoothness)
    return normalize_float(wet_sand)


def compute_local_std(gray, ksize=7):
    mean = cv2.blur(gray, (ksize, ksize))
    mean_sq = cv2.blur(gray * gray, (ksize, ksize))
    variance = np.maximum(mean_sq - mean * mean, 0.0)
    return np.sqrt(variance)


def compute_local_entropy_surrogate(gray, ksize=9, bins=16):
    """
    Cheap texture-entropy stand-in: bins local intensities coarsely,
    then measures how "spread out" the local histogram is via the
    local standard deviation of the *quantized* image. True Shannon
    entropy over a sliding window is expensive; this captures the
    same intuition (foam is locally busy/varied, not flat) far more
    cheaply and is good enough as one of several texture cues.
    """
    quantized = np.floor(gray / (256.0 / bins))
    return compute_local_std(quantized.astype(np.float32), ksize=ksize)


def compute_foam_likelihood(data):
    brightness = normalize_float(data.value.astype(np.float32))
    sat = normalize_float(data.saturation.astype(np.float32))

    lab_b = data.lab_b.astype(np.float32)
    coolness = normalize_float(255.0 - lab_b)

    gray = data.gray.astype(np.float32)

    # Texture cues: foam has a distinctive busy/turbulent texture that
    # plain color/brightness misses on thin, sunlit, or partially
    # submerged foam. Combine local std, Laplacian magnitude, and an
    # entropy surrogate into a single texture score.
    local_std = normalize_float(compute_local_std(gray, ksize=7))

    laplacian = cv2.Laplacian(gray, cv2.CV_32F, ksize=3)
    laplacian_mag = normalize_float(np.abs(laplacian))

    entropy_surrogate = normalize_float(compute_local_entropy_surrogate(gray))

    texture = (0.4 * local_std + 0.35 * laplacian_mag + 0.25 * entropy_surrogate)
    texture = normalize_float(texture)

    color_score = brightness * (1.0 - sat) * coolness

    # Blend color-based foam score with texture, rather than requiring
    # both -- this lets thin/sunlit/partially-submerged foam (weak
    # color signal, strong texture) still register, instead of only
    # catching obvious bright-white breakers.
    foam = 0.55 * color_score + 0.45 * texture
    return normalize_float(np.power(foam, 1.4))


def band_mean(arr, lo, hi):
    """
    For every row r (per column), the mean of arr over rows
    [r+lo, r+hi) -- i.e. a band offset from each candidate row, with
    hi exclusive. lo/hi can be negative (band above r) or positive
    (band below r). Rows of the band that fall outside the image are
    simply excluded from the average (not treated as zero), so bands
    near the top/bottom edge aren't unfairly penalized just for being
    partially off-frame. Vectorized across columns via a cumulative
    sum; the outer loop is only over rows (not rows*cols).
    """
    rows, cols = arr.shape
    csum = np.cumsum(arr, axis=0)
    csum = np.vstack([np.zeros((1, cols), dtype=np.float32), csum])  # csum[i] = sum(arr[0:i])

    result = np.zeros((rows, cols), dtype=np.float32)
    for r in range(rows):
        a = int(np.clip(r + lo, 0, rows))
        b = int(np.clip(r + hi, 0, rows))
        if b <= a:
            # The band is ENTIRELY off-frame. Previously this returned
            # 0.0, which is the maximum possible penalty and the exact
            # opposite of this function's documented intent -- for the
            # water zone (-40,-15) it silently zeroed the "open water
            # above" evidence for every candidate row within 15px of
            # the crop top, discounting the profile score by 3.3x and
            # the final energy (profile squared) by ~11x, purely for
            # being near the top of the frame. That systematically
            # biased the detector AWAY from the top of the crop, and
            # would have sabotaged any attempt to widen a camera's
            # envelope upward to reach a high waterline.
            #
            # Instead, fall back to the nearest in-frame row -- i.e.
            # extrapolate the closest available evidence, which is
            # what "don't penalize for being off-frame" actually
            # means. A band above the frame uses row 0; a band below
            # uses the last row.
            edge = 0 if (r + hi) <= 0 else rows - 1
            result[r, :] = arr[edge, :]
        else:
            result[r, :] = (csum[b] - csum[a]) / (b - a)
    return result


def band_mean_scaled(arr, lo, hi, scale_near, scale_far):
    """
    Perspective-aware version of band_mean: instead of one fixed band
    width for every column, the band width scales per column between
    scale_near (at column 0, near the camera) and scale_far (at the
    last column, far from the camera) -- same near/far convention as
    build_step_schedule. Computed as a blend of two full band_mean
    passes (at the near and far extremes) rather than a true
    per-column variable-width band, since that keeps cost bounded
    (2x band_mean instead of a per-column loop) at the expense of
    being an approximation between the two extremes rather than an
    exact interpolation.

    If either scale is None or both are 1.0, behaves exactly like the
    unscaled band_mean (no behavior change for cameras that haven't
    set profile_scale_near/profile_scale_far).
    """
    if scale_near is None or scale_far is None or (scale_near == 1.0 and scale_far == 1.0):
        return band_mean(arr, lo, hi)

    near_lo, near_hi = int(round(lo * scale_near)), int(round(hi * scale_near))
    far_lo, far_hi = int(round(lo * scale_far)), int(round(hi * scale_far))

    near_map = band_mean(arr, near_lo, near_hi)
    far_map = band_mean(arr, far_lo, far_hi)

    cols = arr.shape[1]
    weight_far = np.linspace(0.0, 1.0, cols, dtype=np.float32)
    return near_map * (1.0 - weight_far) + far_map * weight_far


def compute_vertical_profile_score(water, foam, wet_sand, dry_sand, profile=None):
    """
    V6: replaces single-pixel transition scoring with a VERTICAL
    PROFILE match. A shoreline isn't one pixel -- it's a transition
    region spanning tens of pixels: open water well above, foam near
    the candidate row, wet sand just below, dry sand further below.
    Scoring each candidate row by how well its whole vertical
    neighborhood matches that expected sequence uses far more context
    than judging the pixel at (row, col) in isolation, and is much
    less prone to a single anomalous pixel creating a false positive.

    Band offsets (in pixels, image rows increase downward so "above"
    is negative and "below" is positive):
      water zone     : rows-40 to rows-15 above the candidate (open water)
      foam zone      : rows-15 to rows+5   around the candidate (breaking wave)
      wet sand zone  : rows+5  to rows+25  below the candidate
      dry sand zone  : rows+25 to rows+55  further below the candidate

    These base widths are scaled per-column via profile.
    profile_scale_near/profile_scale_far when set (see
    band_mean_scaled) to account for perspective: bias_analysis.py
    showed error growing sharply toward the far column on both
    cameras, consistent with these fixed-pixel zones being too tall
    (spanning multiple real classes at once) far from the camera.
    """
    scale_near = getattr(profile, "profile_scale_near", None) if profile else None
    scale_far = getattr(profile, "profile_scale_far", None) if profile else None

    water_zone = band_mean_scaled(water, -40, -15, scale_near, scale_far)
    foam_zone = band_mean_scaled(foam, -15, 5, scale_near, scale_far)
    wet_sand_zone = band_mean_scaled(wet_sand, 5, 25, scale_near, scale_far)
    dry_sand_zone = band_mean_scaled(dry_sand, 25, 55, scale_near, scale_far)

    # Additive floor on each zone (same rationale as the hybrid energy
    # formula): a single weak zone should discount the profile match,
    # not collapse it to near-zero the way a pure product would.
    profile_score = (
        (0.3 + 0.7 * water_zone)
        * (0.3 + 0.7 * foam_zone)
        * (0.3 + 0.7 * wet_sand_zone)
        * (0.3 + 0.7 * dry_sand_zone)
    )
    return normalize_float(profile_score)


def compute_negative_gradient(data):
    gray = data.gray.astype(np.float32)
    blurred = cv2.GaussianBlur(gray, (5, 5), 0)
    gy = cv2.Sobel(blurred, cv2.CV_32F, 0, 1, ksize=3)
    return normalize_float(np.abs(gy))


def apply_horizontal_continuity(energy):
    """
    Blends each column's energy with its immediate left/right
    neighbors, using CONFIG.horiz_center_weight /
    CONFIG.horiz_neighbor_weight. This reinforces genuine shoreline
    signal that persists across nearby columns, and suppresses a
    single anomalous column's isolated spike, since its neighbors
    (which won't share that spike) pull the blended value back down.
    """
    center_w = CONFIG.horiz_center_weight
    neighbor_w = CONFIG.horiz_neighbor_weight

    left = np.roll(energy, 1, axis=1)
    left[:, 0] = energy[:, 0]  # no wraparound at the edge

    right = np.roll(energy, -1, axis=1)
    right[:, -1] = energy[:, -1]  # no wraparound at the edge

    blended = center_w * energy + neighbor_w * left + neighbor_w * right
    return normalize_float(blended)


def build_shoreline_energy(data):
    water = compute_water_likelihood(data)
    dry_sand = compute_dry_sand_likelihood(data)
    wet_sand = compute_wet_sand_likelihood(data)
    foam = compute_foam_likelihood(data)
    profile = compute_vertical_profile_score(water, foam, wet_sand, dry_sand, profile=data.profile)
    gradient = normalize_float(np.power(compute_negative_gradient(data), CONFIG.gradient_power))

    # Smooth power-law water weighting, replacing the old hard
    # threshold gate. Calm offshore water can still score reasonably
    # on raw "water_likelihood", but only VERY water-like pixels
    # contribute strongly here, sharply discounting marginal/ambiguous
    # regions rather than treating any value above a cutoff as equally
    # valid.
    water_weight = np.power(np.clip(water, 0.0, 1.0), CONFIG.water_weight_power)

    sand_suppression = np.clip(1.0 - (dry_sand - CONFIG.sand_gate_threshold) * 3.0, 0.0, 1.0)
    sand_suppression = np.where(dry_sand > CONFIG.sand_gate_threshold, sand_suppression, 1.0)

    # V6: profile (vertical neighborhood match) replaces the old
    # single-pixel transition term as the primary discriminating
    # signal. foam/gradient remain secondary boosts with an additive
    # floor so a locally weak one only discounts, not zeroes, the
    # score.
    energy = (
        np.power(profile, 2.0)
        * np.power(water_weight, 1.5)
        * (0.4 + 0.6 * foam)
        * (0.4 + 0.6 * gradient)
        * sand_suppression
    )

    # Camera-specific per-column search envelope: zero out anything
    # outside the allowed band BEFORE horizontal blending and the DP
    # solver ever see it, so physically impossible shoreline locations
    # (e.g. far offshore water) can never be selected, regardless of
    # how strong their raw energy happens to be.
    envelope = build_envelope_mask(energy.shape[0], energy.shape[1], data.profile)
    energy = energy * envelope.astype(np.float32)

    # Horizontal continuity blending across neighboring columns.
    energy = apply_horizontal_continuity(energy)

    # Capture the RAW peak before normalisation destroys absolute
    # scale. Without this there is no way to distinguish a strong frame
    # from a uniformly weak one after the fact.
    raw_energy_peak = float(np.max(energy))

    return normalize_float(energy), {
        "water": water, "dry_sand": dry_sand, "wet_sand": wet_sand, "foam": foam,
        "profile": profile, "gradient": gradient, "envelope": envelope,
        "raw_energy_peak": raw_energy_peak,
    }


def save_energy_debug(data, energy, components):
    if not SAVE_DEBUG_IMAGES:
        return
    # components now also carries scalar diagnostics; only the arrays
    # are image-saveable.
    directory = debug_directory(data.filename)
    save_image(directory / "energy_water.png", normalize(components["water"]))
    save_image(directory / "energy_dry_sand.png", normalize(components["dry_sand"]))
    save_image(directory / "energy_wet_sand.png", normalize(components["wet_sand"]))
    save_image(directory / "energy_foam.png", normalize(components["foam"]))
    save_image(directory / "energy_profile.png", normalize(components["profile"]))
    save_image(directory / "energy_gradient.png", normalize(components["gradient"]))
    save_image(directory / "energy_envelope.png", (components["envelope"].astype(np.uint8) * 255))
    save_image(directory / "energy_final.png", normalize(energy))


def smooth_energy(energy):
    kernel = CONFIG.smoothing_kernel
    if kernel % 2 == 0:
        kernel += 1
    return normalize_float(cv2.GaussianBlur(energy, (kernel, kernel), 0))


# =====================================================================
# DYNAMIC PROGRAMMING SOLVER (simple single-state version)
# =====================================================================

@dataclass
class DPResult:
    shoreline: np.ndarray = None          # integer row per column (backward compatible)
    subpixel_shoreline: np.ndarray = None  # float row per column, sub-pixel refined
    peak_sharpness: np.ndarray = None      # 0-1 normalized, higher = better-localized peak
    confidence: np.ndarray = None
    image_height: int = 0


def build_step_schedule(cols, profile):
    """
    Per-column allowed vertical step. Perspective in an oblique Argus
    camera means the shoreline can move many pixels between adjacent
    columns near the camera, but only a few pixels far away. Column 0
    is treated as "near" (left edge) and the last column as "far"
    (right edge); step_near/step_far are interpolated linearly across
    columns. Falls back to the flat max_vertical_step if either is
    unset, preserving old behavior for any profile that hasn't been
    tuned yet.
    """
    if profile.step_near is None or profile.step_far is None:
        return np.full(cols, profile.max_vertical_step, dtype=np.int32)
    steps = np.linspace(profile.step_near, profile.step_far, cols)
    return np.round(steps).astype(np.int32)


def refine_subpixel(energy, shoreline_int):
    """
    Refines each column's integer DP row to a sub-pixel estimate via
    parabolic interpolation of the energy profile around that row --
    the standard peak-interpolation formula used for correlation/edge
    localization (fit a parabola through (row-1, row, row+1) energy
    values, offset = 0.5*(E[-1]-E[+1]) / (E[-1] - 2*E[0] + E[+1])).

    Also returns a 0-1 "peak sharpness" score per column, derived from
    the same parabola's curvature: a sharp, well-defined energy peak
    (strongly curved) gets a high score; a flat/ambiguous local energy
    landscape gets a low one. This is a RELATIVE quality indicator,
    not a calibrated pixel-uncertainty number -- turning it into an
    actual +/- pixel confidence interval would need to be validated
    against ground truth (e.g. checking whether low-sharpness columns
    really do show larger error in batch_bias_analysis.py) before it
    should be trusted quantitatively.

    At the top/bottom image edge (no row-1 or row+1 neighbor), the
    missing neighbor is replaced with the row itself (zero gradient
    that direction), which correctly yields zero offset and low
    sharpness there rather than crashing or extrapolating.
    """
    rows, cols = energy.shape
    r = shoreline_int.astype(np.int64)

    r_minus = np.clip(r - 1, 0, rows - 1)
    r_plus = np.clip(r + 1, 0, rows - 1)
    col_idx = np.arange(cols)

    e_minus = energy[r_minus, col_idx]
    e0 = energy[r, col_idx]
    e_plus = energy[r_plus, col_idx]

    denom = (e_minus - 2.0 * e0 + e_plus)
    safe_denom = np.where(np.abs(denom) > 1e-8, denom, 1.0)
    raw_offset = 0.5 * (e_minus - e_plus) / safe_denom
    offset = np.where(np.abs(denom) > 1e-8, raw_offset, 0.0)
    offset = np.clip(offset, -0.5, 0.5)

    subpixel_shoreline = r.astype(np.float32) + offset.astype(np.float32)

    # Sharpness: positive curvature magnitude at a true local max
    # (denom < 0 for a peak); anything else (denom >= 0, i.e. not
    # actually a local max -- can happen at envelope edges or flat
    # regions) scores zero rather than a misleading negative value.
    sharpness_raw = np.where(denom < 0, -denom, 0.0)
    peak_sharpness = normalize_float(sharpness_raw)

    return subpixel_shoreline, peak_sharpness.astype(np.float32)


def assess_signal(energy, envelope=None):
    """
    Measures where the energy field actually carries usable signal.

    Returns (column_has_signal, stats). A column counts as having
    signal when its peak reaches COLUMN_SIGNAL_FLOOR_FRACTION of the
    frame's peak. Only rows INSIDE the envelope are considered, since
    energy outside it is zeroed and cannot be selected anyway --
    including it would understate the peak and wrongly condemn good
    columns.

    This is deliberately a RELATIVE test. After normalize_float() the
    frame is scaled to [0,1], so an absolute threshold here would be
    meaningless; what survives normalisation is the contrast BETWEEN
    regions of the same frame, and that is what distinguishes a dead
    far-field from a live one.
    """
    if envelope is not None:
        masked = np.where(envelope, energy, 0.0)
    else:
        masked = energy

    column_peak = masked.max(axis=0)
    frame_peak = float(column_peak.max())

    if frame_peak <= 1e-9:
        return np.zeros(energy.shape[1], dtype=bool), {
            "frame_peak": 0.0, "signal_fraction": 0.0,
            "median_column_peak": 0.0, "median_prominence": 0.0,
        }

    # Test 1 -- regional: is this column comparable to the strongest
    # part of the frame, or is it a dead zone beside a live one?
    passes_fraction = column_peak >= (COLUMN_SIGNAL_FLOOR_FRACTION * frame_peak)

    # Test 2 -- prominence: does the peak stand out from this column's
    # own distribution, or is the column featureless? Uses the median
    # rather than the mean so a single bright row cannot inflate the
    # baseline it is being compared against.
    # Median over rows INSIDE the envelope only. Using `masked` here
    # would average in the zeroed-out region beyond the envelope,
    # driving the median to ~0 and the ratio to ~1e8 -- which silently
    # made this test pass everything on real data, including a
    # featureless grey frame.
    if envelope is not None:
        inside = np.where(envelope, energy, np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            column_median = np.nanmedian(inside, axis=0)
        column_median = np.nan_to_num(column_median, nan=0.0)
    else:
        column_median = np.median(masked, axis=0)
    prominence = column_peak / np.maximum(column_median, 1e-9)
    passes_prominence = prominence >= COLUMN_PEAK_PROMINENCE_MIN

    column_has_signal = passes_fraction & passes_prominence
    return column_has_signal, {
        "frame_peak": frame_peak,
        "signal_fraction": float(column_has_signal.mean()),
        "median_column_peak": float(np.median(column_peak) / frame_peak),
        "median_prominence": float(np.median(prominence)),
    }


def solve_shoreline_vectorized(energy, profile):
    rows, cols = energy.shape
    step_schedule = build_step_schedule(cols, profile)
    max_step = int(step_schedule.max())
    vertical_penalty = CONFIG.vertical_penalty

    score = np.full((rows, cols), -1e30, dtype=np.float32)
    prev_idx = np.full((rows, cols), -1, dtype=np.int32)
    score[:, 0] = energy[:, 0]

    offsets = np.arange(-max_step, max_step + 1)

    for col in range(1, cols):
        col_step = step_schedule[col]
        prev_col_score = score[:, col - 1]
        best_score = np.full(rows, -1e30, dtype=np.float32)
        best_prev = np.full(rows, -1, dtype=np.int32)

        for offset in offsets:
            if abs(offset) > col_step:
                continue  # outside this column's allowed step range

            candidate_rows = np.arange(rows) - offset
            valid = (candidate_rows >= 0) & (candidate_rows < rows)

            candidate_score = np.full(rows, -1e30, dtype=np.float32)
            candidate_score[valid] = (
                prev_col_score[candidate_rows[valid]]
                + energy[np.arange(rows)[valid], col]
                - abs(offset) * vertical_penalty
            )

            improved = candidate_score > best_score
            best_score = np.where(improved, candidate_score, best_score)
            best_prev = np.where(improved, candidate_rows, best_prev)

        score[:, col] = best_score
        prev_idx[:, col] = best_prev

    shoreline = np.zeros(cols, dtype=np.int32)
    row = int(np.argmax(score[:, cols - 1]))
    shoreline[cols - 1] = row
    for col in range(cols - 1, 0, -1):
        row = prev_idx[row, col]
        if row < 0:
            row = shoreline[col]
        shoreline[col - 1] = row

    col_max = score.max(axis=0)
    col_min = score.min(axis=0)
    span = np.maximum(col_max - col_min, 1e-6)
    confidence = np.array([
        (score[shoreline[c], c] - col_min[c]) / span[c] for c in range(cols)
    ], dtype=np.float32)

    subpixel_shoreline, peak_sharpness = refine_subpixel(energy, shoreline)

    result = DPResult()
    result.shoreline = shoreline
    result.subpixel_shoreline = subpixel_shoreline
    result.peak_sharpness = peak_sharpness
    result.confidence = confidence
    result.image_height = rows
    return result


def smooth_shoreline(shoreline, kernel_size):
    """
    Averages each column with its neighbors across a horizontal
    window. Operates on and returns float32 (NOT rounded to int) so a
    sub-pixel-refined shoreline stays sub-pixel through this step
    rather than being truncated back to whole pixels immediately after
    the refinement that computed it.
    """
    if kernel_size < 3:
        return shoreline.astype(np.float32)
    if kernel_size % 2 == 0:
        kernel_size += 1
    radius = kernel_size // 2
    padded = np.pad(shoreline.astype(np.float32), radius, mode="edge")
    smoothed = np.array([np.mean(padded[i:i + kernel_size]) for i in range(len(shoreline))])
    return smoothed.astype(np.float32)


def temporal_filter(shoreline, camera_name):
    """
    Blends with the previous frame's (pre-correction) shoreline for
    the SAME camera. Operates on and returns float32 to preserve
    sub-pixel precision through the blend.
    """
    previous = STATE.previous_shoreline.get(camera_name)
    if previous is None or len(previous) != len(shoreline):
        return shoreline.astype(np.float32)
    alpha = 0.70
    blended = alpha * shoreline.astype(np.float32) + (1.0 - alpha) * previous.astype(np.float32)
    return blended.astype(np.float32)


def apply_bias_correction(shoreline, profile):
    """
    Subtracts a per-column empirical correction derived from
    batch_bias_analysis.py's per-camera segment-wise mean signed
    error (see profile.bias_correction_points). If no points are set
    (default), returns shoreline unchanged -- this is opt-in per
    camera, not a behavior change for unconfigured profiles.
    """
    # Select by image product. Falls back to the snap points if no
    # timex-specific set exists, with a warning, rather than silently
    # applying a correction fitted to a different product.
    if "timex" in IMAGE_SUFFIX.lower():
        points = profile.bias_correction_points_timex
        if not points:
            points = profile.bias_correction_points
            if points and not STATE.warned_timex_bias:
                print("WARNING: no timex-specific bias correction for this camera; falling "
                      "back to the SNAP correction, which was fitted to a different image "
                      "product and will be wrong. Derive one with "
                      "derive_bias_correction.py --suffix .timex. --uncorrected")
                STATE.warned_timex_bias = True
    else:
        points = profile.bias_correction_points

    if not points:
        return shoreline.astype(np.float32)

    cols = len(shoreline)
    xs = np.array([p[0] for p in points])
    corrections = np.array([p[1] for p in points])

    col_x_fracs = np.linspace(0.0, 1.0, cols)
    correction_per_col = np.interp(col_x_fracs, xs, corrections)

    corrected = shoreline.astype(np.float32) - correction_per_col.astype(np.float32)
    return corrected.astype(np.float32)


def clip_shoreline(shoreline, image_height):
    """Clips to valid image bounds. Preserves float32 precision -- only
    the final CSV export rounds to int, for the legacy integer column."""
    return np.clip(shoreline, 0, image_height - 1).astype(np.float32)


# =====================================================================
# OUTPUT
# =====================================================================

def draw_shoreline(image, shoreline, offset_y, offset_x, color=(0, 0, 255), thickness=2):
    overlay = image.copy()
    points = [(int(x + offset_x), int(y + offset_y)) for x, y in enumerate(shoreline)]
    for i in range(len(points) - 1):
        cv2.line(overlay, points[i], points[i + 1], color, thickness, cv2.LINE_AA)
    return overlay


def export_csv(data, shoreline, confidence, peak_sharpness=None, column_has_signal=None):
    """
    Writes the detected shoreline. "Row" stays a rounded integer for
    backward compatibility with every existing consumer script
    (compare_to_ground_truth.py, bias_analysis.py, batch_bias_analysis.py,
    dual_overlay.py, extract_elevation_contours.py all read Column/Row
    positionally by index -- adding trailing columns doesn't break any
    of them, verified by re-checking each one's parser).

    "Row_Precise" is the SAME position at full sub-pixel float
    precision -- this is what should be used for anything where
    fractional-pixel accuracy actually matters (e.g. elevation
    contour extraction), rather than "Row".

    "Peak_Sharpness" (0-1) reflects how well-defined the ORIGINAL raw
    per-column energy peak was, before smoothing/temporal blending/
    bias correction -- a relative data-quality indicator, not a
    calibrated pixel-uncertainty value (see refine_subpixel()'s
    docstring).
    """
    if not EXPORT_CSV:
        return
    csv_name = Path(OUTPUT_FOLDER) / (Path(data.filename).stem + ".csv")
    with open(csv_name, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Column", "Row", "Confidence", "Row_Precise", "Peak_Sharpness",
                         "Has_Signal"])
        for x in range(len(shoreline)):
            precise_row = float(shoreline[x]) + data.crop_top
            full_image_row = int(round(precise_row))
            sharpness = float(peak_sharpness[x]) if peak_sharpness is not None else ""
            # Has_Signal=0 marks a column where the energy peak was
            # below the floor, i.e. the row here was chosen with no
            # real evidence. Such columns should be EXCLUDED from
            # elevation contours rather than trusted equally.
            has_sig = int(bool(column_has_signal[x])) if column_has_signal is not None else ""
            writer.writerow([x, full_image_row, float(confidence[x]), precise_row,
                             sharpness, has_sig])


def export_overlay(data, shoreline):
    overlay = draw_shoreline(data.original, shoreline, data.crop_top, data.crop_left)
    filename = Path(OUTPUT_FOLDER) / (Path(data.filename).stem + "_overlay.png")
    save_image(filename, overlay)


# =====================================================================
# MAIN PIPELINE
# =====================================================================

def measure_brightness(image):
    """Mean grayscale brightness (0-255) of a BGR image."""
    if image.ndim == 3:
        gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    else:
        gray = image
    return float(np.mean(gray))


def process_image(image_path):
    print()
    print("=" * 70)
    print("Processing:", image_path.name)
    print("=" * 70)

    try:
        data = prepare_image(image_path)

        # Reject night/dark frames BEFORE detection. See
        # NIGHT_BRIGHTNESS_THRESHOLD -- the confidence metric cannot
        # catch these (black frames have scored >0.93), so this gate is
        # what stops meaningless night detections reaching the CSVs.
        brightness = measure_brightness(data.original)
        if brightness < NIGHT_BRIGHTNESS_THRESHOLD:
            print(f"Mean brightness : {brightness:.2f}")
            print(f"DISCARDED (too dark: mean brightness {brightness:.2f} < "
                  f"{NIGHT_BRIGHTNESS_THRESHOLD} threshold) -- night/no-light frame, "
                  f"no usable waterline. Not saved to {OUTPUT_FOLDER}.")
            return None, "dark"

        energy, components = build_shoreline_energy(data)
        save_energy_debug(data, energy, components)

        smoothed_energy = smooth_energy(energy)

        result = solve_shoreline_vectorized(smoothed_energy, data.profile)

        # Signal checks. These run BEFORE the confidence gate because
        # confidence cannot see absence of signal -- a frame can be
        # simultaneously 0.99 "confident" and entirely fabricated.
        column_has_signal, sig = assess_signal(smoothed_energy, components.get("envelope"))
        # Absolute image contrast, measured on the CROPPED region the
        # detector actually analyses (not the whole frame, which
        # includes sky and dune).
        crop_gray = cv2.cvtColor(data.cropped, cv2.COLOR_BGR2GRAY) \
            if data.cropped.ndim == 3 else data.cropped
        crop_contrast = float(np.std(crop_gray))
        print("Image contrast  :", f"{crop_contrast:.1f}")
        print("Signal columns  :", f"{100*sig['signal_fraction']:.1f}%"
              f"  (median column peak {100*sig['median_column_peak']:.1f}% of frame peak,"
              f" median prominence {sig['median_prominence']:.1f}x)")

        if IMAGE_CONTRAST_FLOOR > 0.0 and crop_contrast < IMAGE_CONTRAST_FLOOR:
            print(f"DISCARDED (image contrast {crop_contrast:.1f} < {IMAGE_CONTRAST_FLOOR} "
                  f"floor) -- whole-frame signal loss (fog, glare, obstruction, flat scene). "
                  f"Not saved to {OUTPUT_FOLDER}.")
            return result, "no_signal"

        # Per-camera band contrast. Distinct from the check above:
        # that one measures the crop the detector analyses against a
        # global floor, this one applies a CAMERA-SPECIFIC floor to the
        # same region. Fog can leave global contrast looking normal
        # (a fogged C1 frame scored 31.8 globally, 5.3 in-band) so this
        # catches what the global test cannot.
        band_floor = getattr(data.profile, "band_contrast_floor", 0.0)
        if band_floor > 0.0 and crop_contrast < band_floor:
            print(f"DISCARDED (band contrast {crop_contrast:.1f} < {band_floor} floor for "
                  f"{data.profile.name}) -- the cropped region is too flat to contain a "
                  f"resolvable waterline, typically fog. Not saved to {OUTPUT_FOLDER}.")
            return result, "no_signal"

        if sig["signal_fraction"] < MIN_SIGNAL_COLUMN_FRACTION:
            print(f"DISCARDED (only {100*sig['signal_fraction']:.1f}% of columns carry signal, "
                  f"need {100*MIN_SIGNAL_COLUMN_FRACTION:.0f}%) -- the solver would be drawing "
                  f"a line through noise across most of the frame. "
                  f"Not saved to {OUTPUT_FOLDER}.")
            return result, "no_signal"

        shoreline = smooth_shoreline(result.subpixel_shoreline, CONFIG.smoothing_kernel)
        shoreline = temporal_filter(shoreline, data.profile.name)

        # Store the PRE-correction shoreline for next frame's temporal
        # blending, not the corrected one. Storing the corrected value
        # here would mean each frame's correction partially leaks into
        # the next frame via temporal_filter's blend, then gets
        # corrected AGAIN on top of that -- a compounding bug that
        # produces a smaller, sign-flipped residual instead of the
        # correction cleanly cancelling out, exactly what
        # batch_bias_analysis.py caught after the first attempt at
        # this. Bias correction should be applied fresh, once, per
        # frame, right before export -- never baked into temporal
        # state.
        STATE.previous_shoreline[data.profile.name] = shoreline.copy()
        STATE.previous_confidence[data.profile.name] = result.confidence.copy()

        shoreline = apply_bias_correction(shoreline, data.profile)
        shoreline = clip_shoreline(shoreline, result.image_height)

        mean_confidence = float(np.mean(result.confidence))
        print("Mean confidence :", mean_confidence)

        if mean_confidence < CONFIDENCE_DISCARD_THRESHOLD:
            print(f"DISCARDED (mean confidence {mean_confidence:.4f} < "
                  f"{CONFIDENCE_DISCARD_THRESHOLD} threshold) -- not saved to {OUTPUT_FOLDER}.")
            return result, "low_confidence"

        export_csv(data, shoreline, result.confidence, result.peak_sharpness,
                   column_has_signal)
        export_overlay(data, shoreline)

        print("Finished.")

        return result, "saved"

    except Exception as e:
        print()
        print("ERROR")
        print(image_path)
        print(str(e))
        return None, "failed"


def process_directory():
    initialize_directories()

    print()
    print("Clearing previous input/output files...")
    clear_directory(INPUT_FOLDER)
    clear_directory(OUTPUT_FOLDER)

    copy_input_images(limit=IMAGES_PER_RUN)

    images = list_input_images()

    print()
    print("Images found:", len(images))

    successful, failed, discarded, dark, nosignal = 0, 0, 0, 0, 0
    start_time = time.time()

    for image in images:
        result, status = process_image(image)
        if status == "failed":
            failed += 1
        elif status == "dark":
            dark += 1
        elif status == "no_signal":
            nosignal += 1
            STATE.processed_images += 1
        elif status == "low_confidence":
            discarded += 1
            STATE.processed_images += 1
        else:
            successful += 1
            STATE.processed_images += 1

    elapsed = time.time() - start_time

    print()
    print("=" * 70)
    print("Batch Complete")
    print("=" * 70)
    print("Successful (saved)         :", successful)
    print(f"Discarded (confidence < {CONFIDENCE_DISCARD_THRESHOLD}):", discarded)
    print(f"Rejected (too dark / night) :", dark)
    print(f"Rejected (no usable signal) :", nosignal)
    print("Failed                     :", failed)
    total_attempted = successful + discarded + dark + nosignal + failed
    if total_attempted > 0:
        print("Average/Image :", round(elapsed / total_attempted, 2), "seconds")
    print("Elapsed :", round(elapsed, 2), "seconds")


def self_test():
    print()
    print("Running detector self-test...")
    initialize_directories()
    images = list_input_images()
    if len(images) == 0:
        print("No images found yet (expected right after clearing input/).")
        return True
    try:
        prepare_image(images[0])
        print("Image loading........PASS")
    except Exception:
        print("Image loading........FAIL")
        return False
    print("Self-test complete.")
    return True


def print_banner():
    print()
    print("=" * 70)
    print(PROGRAM_NAME)
    print(VERSION)
    print("=" * 70)
    print("OpenCV :", cv2.__version__)
    print("NumPy  :", np.__version__)
    print()
    print("Image product:", IMAGE_SUFFIX)
    print("Source Folder:", SOURCE_FOLDER)
    print("Input Folder :", INPUT_FOLDER)
    print("Output Folder:", OUTPUT_FOLDER)
    print("Debug Folder :", DEBUG_FOLDER)
    print("Images per run:", "ALL (no cap)" if IMAGES_PER_RUN is None else IMAGES_PER_RUN)
    print("Confidence discard threshold:", CONFIDENCE_DISCARD_THRESHOLD)
    print()


def apply_cli_overrides():
    """
    Optional command-line overrides. With no arguments the behavior is
    EXACTLY as before -- this exists so an alternative configuration
    (e.g. detecting on timex instead of snap imagery) can be run into
    separate folders without disturbing the working setup or its
    tuned parameters.
    """
    global SOURCE_FOLDER, INPUT_FOLDER, OUTPUT_FOLDER, DEBUG_FOLDER, IMAGE_SUFFIX

    parser = argparse.ArgumentParser(
        description="USGS Argus shoreline detector. All arguments are optional; "
                    "with none supplied the built-in configuration is used unchanged."
    )
    parser.add_argument("--image-suffix", default=None,
                        help='Which Argus image product to detect on, e.g. "snap.jpg" '
                             '(default, single instantaneous frame) or "timex.jpg" '
                             '(time-exposure, whose bright band represents the MEAN swash '
                             'excursion). NOTE: the envelope, profile-matcher scaling and '
                             'bias corrections were all tuned against snap imagery, so a '
                             'timex run needs its own ground-truth validation before its '
                             'numbers can be compared like-for-like.')
    parser.add_argument("--source-dir", default=None)
    parser.add_argument("--input-dir", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--debug-dir", default=None)
    parser.add_argument("--no-bias-correction", action="store_true",
                        help="Run with bias_correction_points disabled for every camera. "
                             "Use this when comparing image products: the existing "
                             "corrections were fit against snap imagery, so leaving them on "
                             "would confound an otherwise clean comparison.")
    args = parser.parse_args()

    if args.image_suffix:
        IMAGE_SUFFIX = args.image_suffix
    if args.source_dir:
        SOURCE_FOLDER = args.source_dir
    if args.input_dir:
        INPUT_FOLDER = args.input_dir
    if args.output_dir:
        OUTPUT_FOLDER = args.output_dir
    if args.debug_dir:
        DEBUG_FOLDER = args.debug_dir

    if args.no_bias_correction:
        for key, profile in CAMERAS.items():
            CAMERAS[key] = replace(profile, bias_correction_points=(),
                                   bias_correction_points_timex=())
        print("BIAS CORRECTION DISABLED for all cameras (--no-bias-correction).")

    return args


def main():
    args = apply_cli_overrides()
    start_time = time.time()
    print_banner()

    try:
        process_directory()
    except KeyboardInterrupt:
        print()
        print("Interrupted by user.")
    finally:
        elapsed = time.time() - start_time
        print()
        print(f"Total elapsed: {elapsed:.2f} seconds")
        if STATE.processed_images > 0:
            print(f"Images processed: {STATE.processed_images}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Wave Height From The Runup Timestacks (Camera Only, Opt-In)
=============================================================
Estimates the significant wave height off Marconi from each 10-minute
c2 runup timestack (archive/ras_c2, 2 Hz): the swash edge that
runup_from_timestack.py tracks is put on the beach, its runup, setup and
swash are measured as Stockdon et al. (2006) defined them, and their
parameterization is inverted for the deep-water wave height H0 in three
ways. OPT-IN: it writes only archive/hs_timestack_c2.csv. The cron,
every plot and every published value are untouched; whether it is ever
shown or blended is decided by docs/VALIDATION.md 4.11. Of the three,
only the infragravity swash (S_IG) is independent of the station's own
wave heights (CIRCULARITY below): hs_sig_m is the camera-only candidate.

THE RULE
 1. EDGE. runup_from_timestack.swash_edge with that script's own settings
    (foam threshold, 15 px tongues, 2.5 s median), on lines 1 and 2 (the
    c2 lines that cross the swash). One change: its backwash limit, 3 m/s,
    is meant on the beach, but it is applied to distances along the line
    on the still-water plane, the sum of the pixel-to-pixel steps (its
    dist). On a synthetic cross-shore line of single-pixel steps, a
    1:7-1:14 beach 230 m away and the real c2 calibration, that path is
    2.0-3.0 m per beach metre over the swash: 1.1-1.8 of it the plane
    geometry (the plane's straight-line span of the swash: c2 is 19.8 m
    up and sees the swash at a grazing angle, so a point above still water
    lands further along its ray), the rest the zigzag of the staircase (a
    line drawn from 1 m ground samples, as make_timestack_line.py draws
    its surf lines, has none; how the runup lines were drawn is not known
    here). So on such a line the limit is 1.0-1.5 m/s on the beach,
    slower than real backwash. It is therefore scaled by the line's own
    path per beach metre over the swash (plane_stretch, which measures both
    parts on the actual line): 3 m/s on the beach (--max-retreat-speed).
    On the synthetic scene the plane limit read the setup +0.09 to +0.20 m
    (rescaled: +0.06 to +0.10 m) and S_inc 10-21% low (rescaled: -3 to
    +14%), and lost the 3 m case entirely (no spectrum, no estimate); but
    for H0 0.8 m the rescale made S_IG worse (+6.5% -> +16.6%) and the
    combined H0 (+9% -> +20%). runup_from_timestack.py's *_horiz_m are
    those path metres (that script is unchanged).
 2. ELEVATION. Each pixel of the line is a ray from the camera, followed
    down to where it meets the station's 7-day DEM (archive/dems: the DEM
    setup_model.EvolvingSlope.for_epoch picks, its window centre nearest
    the stack and within 7 days) in 1 cm steps, read bilinearly
    (survey_compare.bilinear) and only between cell centres whose weighted
    cells all exist. Past the DEM's ends the profile goes on as the
    straight line through its last 10 m, if that rises landward, for at
    most 10 m (--max-extrapolate; z_extrapolated_fraction is the share of
    the edge read there); gaps up to 6 m are bridged. Without a DEM,
    --beta gives a plane through the Jan 2025 0 m contour
    (timestack_wave_period.py's beach geometry; profile 'planar'). The
    edge's elevation in each 2 Hz row is the profile at its (fractional)
    pixel, minus still water (marconi_water_level.py, GPS first,
    interpolated across the burst; --level-csv another record). Plane
    positions times a slope would not do: on the synthetic scene the plane
    path spread the swash 2.1-2.7 times wider (sd) than the beach did.
    THE LINE'S END. An edge in find_edge's land reference (the line's
    landward 15%, less 2 px: line_end) is capped by the line and its foam
    would sit in the detector's no-foam reference; it counts as past the
    line, like an edge past the profile. (Foam there also lifts the row's
    reference brightness, which can pull the detected edge seaward: that
    is not detectable here.)
    STILL WATER IS NOT SETUP-FREE. Stockdon's R2 and <eta> are measured
    above the tide. The GNSS-R footprint is the surf zone, where breaking
    waves raise the mean level: its level already holds part of the setup
    (historical_forcing.py; its report measures the share), and the
    Chatham fallback, fitted to the GPS, holds the share's mean. Setup and
    R2 therefore read low by share x setup, more the bigger the waves: at
    H0 2 m, T0 10 s, beta_f 0.10 (setup 0.62 m) a share of 20% reads H0
    from setup 1.28 m and from R2 1.71 m, 40% 0.72 m and 1.44 m. S_IG is
    unaffected (the level is smooth over 10 minutes). --compare splits the
    scores by water_level_source; --no-gps (the Chatham transfer: only the
    share's mean, a constant) or a --level-csv from a transfer fitted on
    calm hours gives R2 and setup a level nearer setup-free.
 3. STATISTICS, Stockdon's definitions. R2: the 98th percentile of the
    uprush maxima above still water (the edge's local maxima, found as
    runup_from_timestack.py finds them); an uprush past the landward end of
    the profile, or of the line, counts as higher than any value on it,
    and if the 98th percentile falls among those, no R2 ('R2 above the
    profile' / 'R2 above the line'). Setup <eta>: the mean (none when over
    2% of the edge is off the profile). Swash: 4 x the standard deviation
    of the series above 0.05 Hz (S_inc) and below (S_IG), Stockdon's split,
    with runup_from_timestack.band_sig (gaps up to 5 s bridged; none with
    10% missing), after the straight line through the series is removed
    (band_sig removes only the mean, and a drift of 0.5 m across the record
    adds 10% or more to S_IG; the tide, at most ~0.1 m in 10 minutes, goes
    with still water anyway). Samples off the profile are bridged too,
    which flattens the uprush tops: S_IG reads low by 0.65-0.85 x the share
    off (Stockdon-generated swash, its top 2 / 5 / 8% taken off: -1.3 /
    -3.5 / -5.8% at H0 3 m, -1.4 / -4.0 / -6.8% at 2 m), so the share is
    written (off_profile_fraction), noted and added to S_IG's error (7.).
 4. FORESHORE SLOPE beta_f: Stockdon's, the mean slope over <eta> +/- 2
    sigma of the swash. The median DEM slope (setup_model.slope_grid: the
    plane through each cell's 5 x 5 window, steepest descent, so the
    cross-shore slope and not the slope along an oblique line) of the cells
    under the line in that band, if 3 or more; else of the whole DEM in the
    band (setup_model.ForeshoreSlope.band); else --beta. Outside 0.01-0.30:
    none (a scarp or a gap edge, setup_model.BETA_RANGE).
 5. PERIOD T0: the camera's Tp for the same stack from timestack_wave_period
    .py (the median of its lines with status ok), else the buoy peak period
    in archive/waves_marconi.csv (44008 DPD) or waves_44008.csv within 1.5 h;
    period_source says which. H0 goes as 1/T0^2: a period 25% too long reads
    H0 36% low, 25% too short 78% high (runup_inversion.py's warning), and
    the camera's Tp is the less stable of its two periods
    (timestack_wave_period.py reports Tm01 for that reason): --period-source
    buoy is the check. A buoy period ties the estimate to the buoy: the
    converted buoy (buoy_transfer.py) uses the same Tp.
 6. INVERSION, closed form (each statistic is proportional to sqrt(H0 L0),
    L0 = g T0^2 / 2 pi, so H0 = (statistic / coefficient)^2 / L0):
      (a) R2    = 1.1 [0.35 b sqrt(H0 L0) + sqrt(H0 L0 (0.563 b^2 + 0.004)) / 2]
                  Eq. 19, solved by runup_inversion.h0_from_runup;
      (b) S_IG  = 0.06 sqrt(H0 L0), Eq. 12: no slope in it, and the swash
                  band that keeps growing with the waves when the incident
                  swash saturates in storms (Guza & Thornton 1982; Ruessink et
                  al. 1998), so the estimate meant for big seas;
      (c) <eta> = 0.35 b sqrt(H0 L0), Eq. 10.
    An R2 or setup at or below still water gives none (an error of level,
    not small waves). These are Stockdon's forms for xi0 = b / sqrt(H0/L0)
    >= 0.3. Below it (dissipative) they give R2 = 0.043 sqrt(H0 L0) (Eq.
    18) and, beside it, setup 0.016 and S_IG 0.046 sqrt(H0 L0) (not
    checkable here against the paper, whose full text was not retrievable,
    but consistent with Eq. 18: 1.1 (0.016 + 0.046 / 2) = 0.043); the
    intermediate forms there would read H0 from S_IG 41% low and, at
    beta_f 0.03, from setup 2.3 x high. Both sets are inverted and combined
    (estimate_h0), each giving its own xi0; the set whose own xi0 agrees
    with it is kept, and when both or neither do (near the threshold) the
    one whose estimates agree better, noted 'regime unsure'. On Marconi's
    0.07-0.14 foreshore xi0 < 0.3 needs a deep-water steepness H0/L0 above
    0.05-0.22, steeper than storm seas usually are; on a flatter DEM slope
    (BETA_RANGE allows 0.01) it happens.
 7. COMBINATION (combine): the inverse-variance mean of log H0. A relative
    error e of a statistic is 2e in H0. The error of each statistic: R2 26%
    (Stockdon's own scatter, 38 cm rms on a mean observed R2 of 144 cm) and
    0.28 m of level (still water now, 0.12 m: dem_from_contours.SIGMA_WL_M;
    the DEM's own level, its median sigma_obs 0.258 m: VALIDATION 4.10);
    S_IG 26% (Stockdon's per-band scatter was not retrievable here, so the R2
    figure is assumed), its sampling error and the share of the edge off
    the profile, no level term; setup 0.21 m (Stockdon's setup scatter
    re-evaluated by Dalinghaus, Coco & Higuera 2023: dem_from_contours
    .SIGMA_SETUP_M) and the 0.28 m of level. The SAMPLING error of S_IG is
    sqrt(2 / nu) / 2 with nu the record's own effective degrees of freedom
    (ig_dof_eff: from its periodogram, corrected for the periodogram's own
    scatter). Ten minutes hold at most 58 (a flat infragravity spectrum:
    9%), but swash spectra are peaked or red: on Gaussian records of 10
    minutes the sd of S_IG was 0.096 flat, 0.116 for a 1/f spectrum and
    0.174 peaked at 0.01 Hz (17 minutes: 0.076 / 0.089 / 0.137), and
    ig_dof_eff predicted 0.094 / 0.109 / 0.162. A SLOPE term: 15% of the
    slope (SLOPE_REL, PROVISIONAL: a 5 x 5 plane fit with the DEM's 0.258 m
    per cell, cells independent, has a slope error of 0.018, 18% of 0.10;
    the 29 Sep - 5 Oct 2026 DEM's foreshore spans 0.089-0.125 at four
    positions, +/-17% about 0.106), times each estimate's sensitivity
    (slope_sensitivity): H0 from setup goes as beta_f^-2, from R2 as
    beta_f^-1.5 (2 beta k'/k = 1.52 at 0.10), from S_IG as the DEM's slope
    squared (its elevations scale with it); a 20% slope error is about 30%
    in H0 from R2 and 40% from setup. The errors are taken at the COMBINED H0,
    iterated, not at each estimate's own (that would give a reading that
    is too high a smaller error and more weight). The estimates are NOT
    independent and the combination treats them as if they were: R2 is
    about 1.1 (<eta> + S / 2) of the same edge series, so its errors move
    with setup's and S_IG's (on Stockdon-generated swash a 30% excess of
    infragravity swash moved H0 from S_IG +69% and from R2 +9 to +15%);
    R2 and setup share the level error; all three share the slope (when
    beta_f is read under the same line, R2's and setup's own elevations
    scale with it too, which partly cancels their beta_f term; that is not
    modelled). h0_sigma_rel (the sigma of log H0) is therefore a floor. A
    planar profile does not know where still water meets the beach: there
    only S_IG is inverted (h0_r2_m and h0_setup_m stay blank). --scatter
    replaces the error model once the NUC has measured each estimate
    (VALIDATION 4.11).
 8. CURRENCY. H0 is a DEEP-WATER height, Stockdon's (their nearshore wave
    records taken back to deep water with linear theory, as setup_model.py
    records from how the paper is cited). One shoaling step is applied, the
    other way: hs_m = H0 x Ks(T0, 21 m) (setup_model.shoaling_coefficient,
    0.91-0.97 for 6-14 s) is the height at the ADCP's 21 m, the currency of
    hs_best and of the OWG (both fitted to the ADCP's wh_4061), so hs_m and
    hs_sig_m (H0 from S_IG alone x Ks) are the numbers to compare with
    them. Refraction is not undone either way (setup_model's choice). Being
    inverted from the swash on this beach, H0 is the height that reaches
    Marconi, sheltering included: the raw buoy 44008 (offshore, any
    direction) compares only through buoy_transfer.py (waves_marconi.csv
    hs_buoy_marconi).

STATUS (first match; note says more)
  dark           mean grey level < 25 (timestack_wave_period.py's 'dark')
  low contrast   the c2 bright frame of the same collection has owg_live.py
                 sharpness < 10 (station_status.py's blurred window: 20-24 Oct
                 2026 read 3.4-4.7, clear days 49-88), or edge_contrast < 0.05.
                 Its numbers are still written, to measure what blur does, but
                 they are not trusted.
  no swash edge  the runup script's own test failed (an edge in < 80% of the
                 rows or < 20 uprushes), or no estimate survived (edge gaps)
  no slope       no DEM within 7 days covers the line, and no --beta
  no period      neither the camera's Tp nor a buoy period
  ok
  WAITING. A stack younger than --wait-days (3) is left for a later run,
  not written, while any input that arrives late is missing: the camera's
  period (its file has no row for the stack yet), a DEM, any period, the
  water level, or -- where the GPS is the level source in use -- the GPS
  level (it arrives 1-2 days late; the cron's runup stage waits for it
  too). That holds whatever the contrast label; only 'dark' and a failed
  edge (whose numbers no input would change) are written at once.
  EDGE CONTRAST is the Michelson contrast (foam - sand) / (foam + sand) of
  the mean grey 2-8 px either side of the tip, median over the uprush maxima
  (where the tip stands still). Blur spreads the tip; spray haze lifts both
  sides. On the synthetic scene 0.21-0.23 clear; with 50% haze and a 2 / 4 /
  6 / 10 px blur 0.096 / 0.065 / 0.047 / 0.027. 0.05 is PROVISIONAL, and on
  that scene too lenient: the 2 and 4 px blurs (0.096, 0.065) stay 'ok'
  while moving hs_m 16-17% and hs_sig_m 34% below the clear stack's (the
  same swash). No real stack was available here; the NUC sets it from 20-24
  Oct 2026 against clear days, where the heights start to drift
  (VALIDATION 4.11).

WHAT IT CANNOT DO
  * CIRCULARITY. The DEM is built from the timex waterlines with the
    station's setup C sqrt(Hs L0) added: Hs = hs_best on the live DEMs,
    which blends the OWG and the converted buoy (in storms the buoy carries
    the blend: marconi_waves.py), and the ADCP's own wh_4061 on the
    survey-date DEMs of Dec 2024 - Mar 2025 (historical_forcing.py). R2 and
    setup read the DEM's LEVEL, so they carry that setup: a DEM made of
    waterlines at the mean swash position is high or low by (model setup -
    true setup) there, and for a stack inside the DEM's own week the
    measured setup tends to C sqrt(Hs_ref L0) itself, i.e. H0 from setup to
    Hs_ref x (0.106 / beta_f)^2 (C = 0.037 = 0.35 x 0.106), and the setup
    inside R2 (<eta> is ~38% of R2 at beta_f 0.10) likewise. So R2, setup and the combined hs_m are NOT independent
    of hs_best, of hs_buoy_marconi or of the ADCP (on a DEM forced by it);
    their agreement with any of them is partly built in, and their bias,
    the very thing VALIDATION 4.11's criterion 3 scores, most. S_IG reads
    only the DEM's SLOPE (elevation differences a few metres apart: a setup
    error common to a frame cancels; a setup that grows with elevation
    over the week is second-order): H0 from S_IG (hs_sig_m) is the only
    camera-only candidate against hs_buoy_marconi or the ADCP (with a buoy
    period it shares the buoy's Tp, 5.). --compare lists the others as
    '[not indep.]'.
  * THE DETECTOR on a synthetic scene (the image the real c2 calibration
    would record of a planar 1:7-1:11 beach whose DEM covers -1.0 to +1.6
    m, with Stockdon-shaped swash whose backwash is at most 2 m/s; four
    cases, H0 0.8-3.0 m: test_timestack_wave_height.py prints the table)
    [SYNTH]: R2 within 0.6% of the swash's own, H0 from R2 within 1% (the
    3 m case ran above the DEM's extension: no R2); setup +0.06 to +0.10 m
    (the edge trails a fast tip: the detector's 3 x 3 blur reaches the
    neighbouring rows), H0 from setup +38 to +59% (part of it the synthetic
    swash's own: its capped backwash lifts its mean above Eq. 10); S_IG +4
    to +17% (the edge drops seaward at the deepest rundowns), H0 from S_IG
    -8 to +32%; combined H0 +11 to +20%. Real foam is not synthetic foam:
    only the NUC comparison measures the real errors.
  * Ten minutes is short for infragravity statistics: S_IG's sampling error
    alone is 9% for a flat spectrum (19% in H0 from it), about 12% for a
    1/f spectrum and 17% for one peaked at 0.01 Hz (23% and 35% in H0);
    ig_dof_eff measures it per stack.
  * Still water holds part of the setup (2.): R2 and setup read low, most
    in storms.
  * The period (5. above).
  * The DEM is a week's waterlines, i.e. the mean swash of the week's tides:
    at high water and in big waves the swash runs above it (10 m of
    extension, then no R2), and a berm built or cut within the week is not in
    it. The line's landward 15% is not read (2.): at high water in storms
    the uprush may reach it, and then there is no R2 ('R2 above the line').
  * Stockdon's relations are fits over ten experiments on other beaches. A
    site coefficient could be fitted to this beach once truth exists, which
    is what VALIDATION 4.11 is for.

OUTPUT archive/hs_timestack_c2.csv, one row per stack and line, APPENDED
(every 50 line-stacks, so a crash loses at most those): rows already there
are never rewritten except by --reprocess, and a file with other columns is
refused before any stack is processed.
  epoch (burst middle), time_utc, filename, line, status, note, brightness,
  edge_contrast, window_sharpness, valid_fraction, n_peaks, sea_at, profile,
  z_extrapolated_fraction, off_profile_fraction (share of the edge off the
  profile or past the line), water_level_m, water_level_source, beta_f,
  beta_source, iribarren, period_s, period_source, setup_m, r2_m, s_inc_m,
  s_ig_m, s_m (m, above still water), h0_r2_m, h0_sig_m, h0_setup_m (each
  estimate, deep water), h0_m, h0_sigma_rel, combined_from (the
  combination), ks_21m, hs_m (H0 x Ks at 21 m), hs_sig_m (H0 from S_IG alone
  x Ks: the camera-only candidate)

References: Stockdon, Holman, Howd & Sallenger 2006, Coastal Engineering
53, 573-588 (its full text was not retrievable here; equations checked
against py-wave-runup's Stockdon2006, which holds Eqs. 10-12, 18 and 19);
USGS Data Series 602 (its runup statistics); Holman & Sallenger 1985, J.
Geophys. Res. 90(C1), 945-953 (setup and swash from video); Guza & Thornton
1982, J. Geophys. Res. 87(C1), 483-491; Ruessink, Kleinhans & van den
Beukel 1998, J. Geophys. Res. 103(C2), 3111-3118; Dalinghaus, Coco &
Higuera 2023, Nat. Hazards Earth Syst. Sci. 23, 2157-2169.

Usage:
    python3 timestack_wave_height.py                  new stacks in archive/ras_c2
    python3 timestack_wave_height.py --reprocess      every archived stack again
    python3 timestack_wave_height.py --compare        against the ADCP / converted buoy, the OWG on the same stacks
    python3 timestack_wave_height.py --compare --since YYYY-MM-DD      only the held-out stacks
    python3 timestack_wave_height.py --self-test
Tests: python3 test_timestack_wave_height.py
"""

import sys
import csv
import copy
import argparse
from pathlib import Path
from datetime import datetime, timezone

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import runup_from_timestack as rft                                          # noqa: E402
from runup_from_timestack import SAMPLE_HZ, IG_CUTOFF_HZ, LAND_REF_FRACTION, band_sig, bridge_gaps, peak_mask  # noqa: E402
from setup_model import (deep_water_wavelength, stockdon_setup, shoaling_coefficient,  # noqa: E402
                         slope_grid, ForeshoreSlope, EvolvingSlope, STOCKDON_SETUP_K,
                         STOCKDON_SWASH_IG, DESHOAL_DEPTH_BY_CURRENCY, BETA_RANGE, SLOPE_MAX_DAYS)
from dem_from_contours import SIGMA_WL_M, SIGMA_SETUP_M                     # noqa: E402
from station_status import LENS_BLURRED                                     # noqa: E402
from timestack_wave_period import SHORE_ORIGIN, SEAWARD_DEG                 # noqa: E402
from survey_compare import bilinear                                         # noqa: E402

FIELDS = ["epoch", "time_utc", "filename", "line", "status", "note",
          "brightness", "edge_contrast", "window_sharpness", "valid_fraction", "n_peaks", "sea_at",
          "profile", "z_extrapolated_fraction", "off_profile_fraction", "water_level_m", "water_level_source",
          "beta_f", "beta_source", "iribarren", "period_s", "period_source",
          "setup_m", "r2_m", "s_inc_m", "s_ig_m", "s_m",
          "h0_r2_m", "h0_sig_m", "h0_setup_m", "h0_m", "h0_sigma_rel", "combined_from",
          "ks_21m", "hs_m", "hs_sig_m"]
STATUSES = ("ok", "low contrast", "no swash edge", "no slope", "no period", "dark")

DEFAULT_LINES = (1, 2)            # the c2 lines that cross the swash (waterline_timex_cron.sh RUNUP_LINES)
DARK_MEAN = 25.0                  # grey level: timestack_wave_period.period_from_stack's "dark"
MIN_CONTRAST = 0.05               # Michelson contrast across the swash tip: provisional, see EDGE CONTRAST
CONTRAST_NEAR, CONTRAST_FAR = 2, 8    # pixels either side of the tip compared by edge_contrast
OWG_MAX_GAP_S = 600.0             # the bright frame of the same 10-min collection
PERIOD_MAX_GAP_S = 5400.0         # a buoy hour within 1.5 h of the burst's middle
WAIT_DAYS = 3.0                   # inputs that arrive late (period, DEM, water level) are waited for
CHUNK_ROWS = 50                   # line-stacks per append to the output (a crash loses at most these)
DESHOAL_DEPTH_M = DESHOAL_DEPTH_BY_CURRENCY["marconi"]     # 21 m: the ADCP / hs_best / OWG currency

# The beach profile along a line
DZ_STEP = 0.01                    # m: the ray is marched down the DEM in these steps
EXT_FIT_M = 10.0                  # m of DEM profile at each end the extension is fitted to
MAX_EXTRAPOLATE_M = 10.0          # m beyond the DEM's edge (along the line) the profile is carried on
GAP_FILL_M = 6.0                  # m: DEM gaps inside the profile up to this long are bridged
MISSING_ELEV_MAX = 0.02           # share of edge samples off the profile above which no setup is given
MIN_SLOPE_CELLS = 3               # distinct DEM cells under the line inside the swash band
LINE_END_PX = 2                   # an edge this close to find_edge's land reference counts as past the line

# Stockdon et al. (2006) below xi0 = 0.3 (dissipative): R2 = 0.043 sqrt(H0 L0) (Eq. 18), setup 0.016 and
# S_IG 0.046 sqrt(H0 L0) (the dissipative forms beside it; 1.1 (0.016 + 0.046 / 2) = 0.043)
XI_DISSIPATIVE = 0.3
DISSIPATIVE_K = {"r2": 0.043, "setup": 0.016, "sig": 0.046}

# Error model of the combination (log H0, inverse variance): per estimate (relative, absolute m,
# slope relative): the scatter of the measured quantity about Stockdon's relation, and the relative
# error of the slope it depends on, see THE COMBINATION.
DEM_SIGMA_M = 0.258               # median per-cell sigma_obs of the station's DEM, VALIDATION 4.10 [REAL]
LEVEL_SIGMA_M = float(np.hypot(SIGMA_WL_M, DEM_SIGMA_M))   # still water now and the DEM's level: 0.28 m
R2_REL = 0.38 / 1.44              # Stockdon et al. (2006): R2 rms 38 cm on a mean observed 144 cm
SLOPE_REL = 0.15                  # beta_f and the DEM's slope, relative: PROVISIONAL (THE COMBINATION)
SCATTER = {"r2": (R2_REL, LEVEL_SIGMA_M, SLOPE_REL),
           "sig": (R2_REL, 0.0, SLOPE_REL),
           "setup": (0.0, float(np.hypot(SIGMA_SETUP_M, LEVEL_SIGMA_M)), SLOPE_REL)}


# ---------------------------------------------------------------------
# Stockdon et al. (2006), forward and inverted
# ---------------------------------------------------------------------

def stockdon_r2(h0, t0, beta, regime="intermediate"):
    """R2 = 1.1 [0.35 beta sqrt(H0 L0) + sqrt(H0 L0 (0.563 beta^2 + 0.004)) / 2]: Stockdon et al.
    (2006) Eq. 19 as published, the relation runup_inversion.h0_from_runup inverts. Its swash
    term rounds Eqs. 11-12 (0.75^2 = 0.5625, 0.06^2 = 0.0036), so 1.1 (<eta> + S/2) from
    setup_model.stockdon_setup / stockdon_swash is 0-3% lower. regime 'dissipative': Eq. 18,
    R2 = 0.043 sqrt(H0 L0)."""
    h0, t0, beta = (np.asarray(a, float) for a in (h0, t0, beta))
    x = h0 * deep_water_wavelength(t0)
    if regime == "dissipative":
        return DISSIPATIVE_K["r2"] * np.sqrt(x)
    return 1.1 * (STOCKDON_SETUP_K * beta * np.sqrt(x) + np.sqrt(x * (0.563 * beta ** 2 + 0.004)) / 2.0)


def stockdon_sig(h0, t0, regime="intermediate"):
    """Infragravity swash S_IG = 0.06 sqrt(H0 L0), Stockdon et al. (2006) Eq. 12 (no slope in it);
    0.046 sqrt(H0 L0) for regime 'dissipative'."""
    k = DISSIPATIVE_K["sig"] if regime == "dissipative" else STOCKDON_SWASH_IG
    return k * np.sqrt(h0 * deep_water_wavelength(t0))


def stockdon_eta(h0, t0, beta, regime="intermediate"):
    """Setup <eta> = 0.35 beta sqrt(H0 L0) (Eq. 10, setup_model.stockdon_setup); 0.016 sqrt(H0 L0)
    for regime 'dissipative'."""
    if regime == "dissipative":
        return DISSIPATIVE_K["setup"] * np.sqrt(np.asarray(h0, float) * deep_water_wavelength(t0))
    return stockdon_setup(h0, t0, beta)


def _period_ok(t0):
    return np.isfinite(t0) and t0 > 0


def h0_from_r2(r2, beta, t0, regime="intermediate"):
    """H0 from R2 above still water: runup_inversion.h0_from_runup (closed form; NaN for R2 <= 0);
    regime 'dissipative': (R2 / 0.043)^2 / L0, no slope."""
    if regime == "dissipative":
        r2, t0 = float(r2), float(t0)
        if not (np.isfinite(r2) and r2 > 0 and _period_ok(t0)):
            return np.nan
        return (r2 / DISSIPATIVE_K["r2"]) ** 2 / float(deep_water_wavelength(t0))
    from runup_inversion import h0_from_runup
    return float(h0_from_runup(np.asarray([r2], float), np.asarray([beta], float),
                               np.asarray([t0], float))[0])


def h0_from_sig(s_ig, t0, regime="intermediate"):
    """H0 = (S_IG / 0.06)^2 / L0: Eq. 12 inverted (0.046 for regime 'dissipative'). NaN for
    S_IG <= 0 or no period."""
    s_ig, t0 = float(s_ig), float(t0)
    if not (np.isfinite(s_ig) and s_ig > 0 and _period_ok(t0)):
        return np.nan
    k = DISSIPATIVE_K["sig"] if regime == "dissipative" else STOCKDON_SWASH_IG
    return (s_ig / k) ** 2 / float(deep_water_wavelength(t0))


def h0_from_setup(setup, beta, t0, regime="intermediate"):
    """H0 = (setup / (0.35 beta))^2 / L0: Eq. 10 inverted ((setup / 0.016)^2 / L0, no slope, for
    regime 'dissipative'). NaN for setup <= 0 (below still water: an error in the level, not small
    waves), no slope or no period."""
    setup, beta, t0 = float(setup), float(beta), float(t0)
    if not (np.isfinite(setup) and setup > 0 and _period_ok(t0)):
        return np.nan
    if regime == "dissipative":
        return (setup / DISSIPATIVE_K["setup"]) ** 2 / float(deep_water_wavelength(t0))
    if not (np.isfinite(beta) and beta > 0):
        return np.nan
    return (setup / (STOCKDON_SETUP_K * beta)) ** 2 / float(deep_water_wavelength(t0))


def iribarren(beta, h0, t0):
    """xi0 = beta / sqrt(H0 / L0); Stockdon's Eq. 19 is for xi0 >= 0.3 (Eq. 18 below)."""
    if not (np.isfinite(beta) and np.isfinite(h0) and h0 > 0 and np.isfinite(t0)):
        return np.nan
    return float(beta / np.sqrt(h0 / deep_water_wavelength(t0)))


def predicted(kind, h0, t0, beta, regime="intermediate"):
    """The swash statistic an estimate was inverted from, as Stockdon predicts it for H0."""
    if kind == "r2":
        return float(stockdon_r2(h0, t0, beta, regime))
    if kind == "sig":
        return float(stockdon_sig(h0, t0, regime))
    if kind == "setup":
        return float(stockdon_eta(h0, t0, beta, regime))
    raise ValueError(kind)


def slope_sensitivity(kind, beta, regime="intermediate"):
    """
    |d ln H0 / d ln slope| of an estimate. 'sig': 2 (S_IG has no slope in Eq. 12, but its
    elevations are read off the DEM, so they scale with the DEM's slope, and H0 with S_IG^2).
    'setup': 2 (H0 goes as beta_f^-2). 'r2': 2 beta k'/k, k = 0.35 beta + sqrt(0.563 beta^2 +
    0.004) / 2 the bracket of Eq. 19: 1.52 at beta 0.10, i.e. H0 ~ beta_f^-1.5. Dissipative: 0 for
    R2 and setup (no slope in Stockdon's dissipative forms). Each is beta_f's effect on the inversion alone:
    R2's and setup's own elevations scale with the DEM's slope too, which, when beta_f is read
    under the same line, partly cancels it (setup's to first order); that is not modelled.
    """
    if kind == "sig":
        return 2.0
    if regime == "dissipative" or not (np.isfinite(beta) and beta > 0):
        return 0.0
    if kind == "setup":
        return 2.0
    root = np.sqrt(0.563 * beta ** 2 + 0.004)
    k = STOCKDON_SETUP_K * beta + root / 2.0
    dk = STOCKDON_SETUP_K + 0.563 * beta / (2.0 * root)
    return float(2.0 * beta * dk / k)


def estimate_sigma(kind, h0, t0, beta, scatter=None, sig_extra=0.0, regime="intermediate"):
    """
    The error (sigma of ln H0) of one estimate, evaluated at H0: its statistic q's error
    sqrt((rel q)^2 + abs^2) (scatter[kind] = (rel, abs m, slope rel); for 'sig' rel also takes
    sig_extra, its sampling and off-profile terms, in quadrature) is 2 s_q / q in ln H0, since q
    goes as sqrt(H0); the slope term slope_sensitivity x slope rel adds in quadrature. inf when
    q <= 0.
    """
    scatter = scatter or SCATTER
    rel, ab, sl = (tuple(scatter[kind]) + (0.0, 0.0))[:3]
    q = predicted(kind, h0, t0, beta, regime)
    if not q > 0:
        return np.inf
    if kind == "sig":
        rel = float(np.hypot(rel, sig_extra))
    s_q = float(np.hypot(rel * q, ab))
    return float(np.hypot(2.0 * s_q / q, slope_sensitivity(kind, beta, regime) * sl))


def combine(estimates, t0, beta, scatter=None, sig_extra=0.0, regime="intermediate", passes=20):
    """
    One H0 from the separate estimates {'r2' | 'sig' | 'setup': H0 m}: the inverse-variance mean
    of log H0, the variance of each from estimate_sigma. It is evaluated at the COMBINED H0, not at
    each estimate's own (that would give an estimate that reads high a smaller relative error and
    so more weight), iterated to a fixed point from the geometric mean. The estimates are treated
    as independent, which they are not (THE COMBINATION): the sigma returned is a floor.
    -> (H0, sigma of log H0, '+'-joined names used); (NaN, NaN, '') with no usable estimate.
    """
    scatter = scatter or SCATTER
    use = {k: float(v) for k, v in estimates.items() if v is not None and np.isfinite(v) and v > 0}
    if not use:
        return np.nan, np.nan, ""
    names = sorted(use)
    logs = np.array([np.log(use[k]) for k in names])
    h = float(np.exp(logs.mean()))
    w = np.ones(len(names))
    for _ in range(passes):
        for i, k in enumerate(names):
            sd = estimate_sigma(k, h, t0, beta, scatter, sig_extra, regime)
            w[i] = 1.0 / sd ** 2 if 0 < sd < np.inf else 0.0
        if w.sum() <= 0:
            return np.nan, np.nan, ""
        h_new = float(np.exp((w * logs).sum() / w.sum()))
        done = abs(h_new - h) <= 1e-9 * h
        h = h_new
        if done:
            break
    return h, float(1.0 / np.sqrt(w.sum())), "+".join(names)


def parse_scatter(text):
    """'r2=0.26:0.28:0.15,sig=0.26:0,setup=0:0.36' (relative:absolute m[:slope relative]; a slope
    term left out keeps the default's) -> a SCATTER dict."""
    out = dict(SCATTER)
    if not text:
        return out
    for part in str(text).split(","):
        k, _, v = part.partition("=")
        k = k.strip()
        if k not in SCATTER:
            raise ValueError(f"--scatter: unknown estimate {k!r} (r2, sig, setup)")
        bits = v.split(":")
        rel = float(bits[0])
        ab = float(bits[1]) if len(bits) > 1 and bits[1] != "" else 0.0
        sl = float(bits[2]) if len(bits) > 2 and bits[2] != "" else SCATTER[k][2]
        out[k] = (rel, ab, sl)
    return out


# ---------------------------------------------------------------------
# The beach under the line: where each pixel's ray meets it
# ---------------------------------------------------------------------

class Dem:
    """An ESRI ASCII DEM (dem_from_contours.py): bilinear elevations, and the foreshore slope grid of
    setup_model.slope_grid (the plane fitted to each cell's 5 x 5 window) on demand."""

    def __init__(self, path=None, grid=None, xll=None, yll=None, cell=None, name=""):
        if path is not None:
            from asc_to_geotiff import read_asc
            grid, head = read_asc(path)
            cell = float(head["cellsize"])
            xll = head["xllcorner"] if "xllcorner" in head else head["xllcenter"] - cell / 2.0
            yll = head["yllcorner"] if "yllcorner" in head else head["yllcenter"] - cell / 2.0
            name = Path(path).name
        self.z = np.asarray(grid, float)                   # rows north to south
        self.x0, self.cell, self.name = float(xll), float(cell), name
        self.ytop = float(yll) + self.z.shape[0] * self.cell
        fin = self.z[np.isfinite(self.z)]
        self.zmin, self.zmax = (float(fin.min()), float(fin.max())) if len(fin) else (np.nan, np.nan)
        self._slope = self._foreshore = None

    def _rc(self, E, N):
        return (self.ytop - np.asarray(N, float)) / self.cell - 0.5, \
               (np.asarray(E, float) - self.x0) / self.cell - 0.5

    def _cell(self, grid, r, c):
        r, c = np.asarray(r, float), np.asarray(c, float)
        fin = np.isfinite(r) & np.isfinite(c)
        rr = np.where(fin, np.floor(np.where(fin, r, 0.0)), -1).astype(int)
        cc = np.where(fin, np.floor(np.where(fin, c, 0.0)), -1).astype(int)
        ok = (rr >= 0) & (rr < grid.shape[0]) & (cc >= 0) & (cc < grid.shape[1])
        out = np.full(np.shape(r), np.nan)
        out[ok] = grid[rr[ok], cc[ok]]
        return out

    def sample(self, E, N):
        """Bilinear between cell centres: survey_compare.bilinear (NaN where a neighbour that carries
        weight is nodata), and NaN in the outer half cell too, where it would hold the edge value."""
        r, c = self._rc(E, N)
        with np.errstate(invalid="ignore"):
            z = bilinear(self.z, self.x0, self.ytop, self.cell, E, N)
            between = (r >= 0) & (r <= self.z.shape[0] - 1) & (c >= 0) & (c <= self.z.shape[1] - 1)
        return np.where(between, z, np.nan)

    def slope_and_cell(self, E, N):
        """Foreshore slope of the nearest cell, and that cell's index (NaN / -1 outside)."""
        if self._slope is None:
            self._slope = slope_grid(self.z, self.cell)
        r, c = self._rc(E, N)
        rr, cc = np.floor(r + 0.5), np.floor(c + 0.5)
        s = self._cell(self._slope, rr, cc)
        idx = np.where(np.isfinite(s), np.nan_to_num(rr) * self.z.shape[1] + np.nan_to_num(cc), -1)
        return s, idx.astype(np.int64)

    def foreshore(self):
        """setup_model.ForeshoreSlope of this DEM (the whole grid), for the band fallback: the one
        EvolvingSlope already holds when it picked this DEM, else computed here."""
        if self._foreshore is None:
            self._foreshore = ForeshoreSlope(self.z, self.cell, self.name)
        return self._foreshore


class DemSource:
    """The DEM for a stack: a dated series (dem_from_contours.py --series-dir, the cron's archive/dems)
    or one DEM file for every stack, picked by setup_model.EvolvingSlope.for_epoch itself -- the DEM
    whose window centre is nearest, within max_days -- so the two can never pick differently."""

    def __init__(self, source, max_days=SLOPE_MAX_DAYS):
        self.ev = EvolvingSlope(source, max_days) if source else None
        self._cache = {}

    def path_for(self, epoch):
        """The DEM file EvolvingSlope picks for epoch (None: none within max_days, or unreadable)."""
        if self.ev is None:
            return None
        fs, _why = self.ev.for_epoch(epoch)
        if fs is None:
            return None
        if self.ev.single is not None:
            return self.ev.single
        return next((p for _, p in self.ev.series if p.name == fs.name), None)

    def for_epoch(self, epoch):
        p = self.path_for(epoch)
        if p is None:
            return None
        if p not in self._cache:
            try:
                self._cache[p] = Dem(p)
                self._cache[p]._foreshore = self.ev.for_epoch(epoch)[0]      # the same ForeshoreSlope
            except (OSError, ValueError, KeyError) as exc:
                print(f"WARNING: cannot read the DEM {p} ({exc})")
                self._cache[p] = None
        return self._cache[p]


def seaward_distance(E, N):
    """Metres seaward of the Jan 2025 0 m NAVD88 contour along the 81 deg beach normal (the beach
    geometry of timestack_wave_period.py and make_timestack_line.py)."""
    a = np.radians(SEAWARD_DEG)
    return (np.asarray(E, float) - SHORE_ORIGIN[0]) * np.sin(a) + (np.asarray(N, float) - SHORE_ORIGIN[1]) * np.cos(a)


def line_profile(u, v, io, eo, dem=None, beta=None, max_extrapolate=MAX_EXTRAPOLATE_M):
    """
    The beach elevation seen by each pixel of a line (sea end first): where the pixel's ray meets
    the surface. A ray's ground point on the horizontal plane at height z is exactly linear in z
    (georectify.pixel_to_ground), so the ray is P(z) = P(0) + z (P(1) - P(0)).

    dem: the ray is marched down from above the DEM's top in DZ_STEP steps; the first crossing of
    the surface (ray above it, then at or below) is refined linearly. The DEM is read bilinearly and
    only where all four surrounding cells exist (at its edges the nearest cell would be up to a cell
    off: 0.1-0.2 m on a 1:10 foreshore). Where the line runs off the DEM, the profile is carried on
    as the straight line (elevation against horizontal distance along the line) fitted to the last
    EXT_FIT_M of DEM profile at that end, for at most max_extrapolate metres (the pixels from the
    DEM's end until the line first reaches that far), rising landward or not at all; gaps inside up
    to GAP_FILL_M are bridged. kind 'dem'.
    beta (no dem): the plane z = -beta x (metres seaward of the Jan 2025 0 m contour). kind 'planar'.
    -> dict z, E, N, d (m along the line on the beach, from the first pixel with an elevation),
       extrap (bool: extended or bridged), slope, cell (DEM foreshore slope and cell under each
       DEM pixel), kind, name.
    """
    import georectify
    u, v = np.asarray(u, float), np.asarray(v, float)
    n = len(u)
    E0, N0 = georectify.pixel_to_ground(u, v, 0.0, io, eo)
    E1, N1 = georectify.pixel_to_ground(u, v, 1.0, io, eo)
    E0, N0 = np.asarray(E0, float), np.asarray(N0, float)
    dE, dN = np.asarray(E1, float) - E0, np.asarray(N1, float) - N0
    out = {"z": np.full(n, np.nan), "E": np.full(n, np.nan), "N": np.full(n, np.nan),
           "d": np.full(n, np.nan), "extrap": np.zeros(n, bool), "slope": np.full(n, np.nan),
           "cell": np.full(n, -1, np.int64), "kind": "none", "name": ""}
    z = np.full(n, np.nan)
    extrap = np.zeros(n, bool)

    def axis_of(zz):
        """The line's horizontal direction ON THE BEACH (its ground points at their own height: the
        plane at one height would turn and stretch it), and d(z) = d0 + d1 z on every ray."""
        ii = np.flatnonzero(np.isfinite(zz))
        Ea, Na = E0[ii[0]] + zz[ii[0]] * dE[ii[0]], N0[ii[0]] + zz[ii[0]] * dN[ii[0]]
        Eb, Nb = E0[ii[-1]] + zz[ii[-1]] * dE[ii[-1]], N0[ii[-1]] + zz[ii[-1]] * dN[ii[-1]]
        ax = np.array([Eb - Ea, Nb - Na]) / max(np.hypot(Eb - Ea, Nb - Na), 1e-9)
        d0 = (E0 - Ea) * ax[0] + (N0 - Na) * ax[1]
        d1 = dE * ax[0] + dN * ax[1]
        return d0, d1

    if dem is not None and np.isfinite(dem.zmin):
        lev = np.arange(dem.zmax + DZ_STEP, dem.zmin - 2 * DZ_STEP, -DZ_STEP)
        Ek = E0[None, :] + lev[:, None] * dE[None, :]
        Nk = N0[None, :] + lev[:, None] * dN[None, :]
        f = lev[:, None] - dem.sample(Ek, Nk)
        with np.errstate(invalid="ignore"):
            cross = (f[:-1] > 0) & (f[1:] <= 0)
        has = cross.any(axis=0)
        k = np.argmax(cross, axis=0)
        cols = np.flatnonzero(has)
        fa, fb = f[k[cols], cols], f[k[cols] + 1, cols]
        la, lb = lev[k[cols]], lev[k[cols] + 1]
        z[cols] = la + fa * (lb - la) / (fa - fb)
        covered = np.isfinite(z)
        if covered.sum() >= 2:
            ci = np.flatnonzero(covered)
            d0, d1 = axis_of(z)
            dcov = d0[ci] + d1[ci] * z[ci]
            zcov = z[ci].copy()
            for end in ("land", "sea"):
                # targets in order away from the DEM's end
                if end == "land":
                    sel = dcov >= dcov.max() - EXT_FIT_M
                    targets = np.arange(ci[-1] + 1, n)
                    d_end, sgn = dcov.max(), 1.0
                else:
                    sel = dcov <= dcov.min() + EXT_FIT_M
                    targets = np.arange(ci[0] - 1, -1, -1)
                    d_end, sgn = dcov.min(), -1.0
                if sel.sum() < 10 or np.ptp(dcov[sel]) < 3.0 or not len(targets):
                    continue
                b, a = np.polyfit(dcov[sel], zcov[sel], 1)
                if not b > 0:
                    continue
                with np.errstate(invalid="ignore", divide="ignore"):
                    zt = (a + b * d0[targets]) / (1.0 - b * d1[targets])
                    dt = d0[targets] + d1[targets] * zt
                    # how far past the DEM's end the line has reached so far: a staircase pixel line
                    # steps back a little at each turn, which must not leave holes in the extension
                    reach = np.maximum.accumulate(np.where(np.isfinite(dt), sgn * (dt - d_end), -np.inf))
                    ok = np.isfinite(zt) & (reach <= max_extrapolate)
                z[targets[ok]] = zt[ok]
                extrap[targets[ok]] = True
            # bridge short gaps inside the DEM profile
            inner = np.arange(ci[0], ci[-1] + 1)
            gap = inner[~np.isfinite(z[inner])]
            if len(gap):
                for r in np.split(gap, np.flatnonzero(np.diff(gap) > 1) + 1):
                    lo, hi = r[0] - 1, r[-1] + 1
                    dl, dh = d0[lo] + d1[lo] * z[lo], d0[hi] + d1[hi] * z[hi]
                    if abs(dh - dl) <= GAP_FILL_M:
                        z[r] = np.interp(r, [lo, hi], [z[lo], z[hi]])
                        extrap[r] = True
            out["kind"], out["name"] = "dem", dem.name
            s, cell = dem.slope_and_cell(np.where(covered, E0 + z * dE, np.nan),
                                         np.where(covered, N0 + z * dN, np.nan))
            out["slope"], out["cell"] = np.where(covered, s, np.nan), np.where(covered, cell, -1)
        else:
            z[:] = np.nan
    elif beta is not None and np.isfinite(beta) and beta > 0:
        s0 = seaward_distance(E0, N0)
        s1 = seaward_distance(E0 + dE, N0 + dN) - s0
        with np.errstate(invalid="ignore", divide="ignore"):
            z = -beta * s0 / (1.0 + beta * s1)
        if np.isfinite(z).sum() >= 2:
            out["kind"], out["name"] = "planar", f"beta {beta:g}"
        else:
            z[:] = np.nan
    out["z"], out["extrap"] = z, extrap & np.isfinite(z)
    out["E"], out["N"] = E0 + z * dE, N0 + z * dN
    if np.isfinite(z).sum() >= 2:
        d0, d1 = axis_of(z)
        out["d"] = d0 + d1 * z
    return out


def line_end(n_px, margin=LINE_END_PX):
    """The first column (from the sea end) of the line's landward end as find_edge sees it: its land
    reference, the last 15% (runup_from_timestack.LAND_REF_FRACTION), less `margin` pixels. An edge
    there is capped by the line, and its foam would sit in find_edge's no-foam reference."""
    return max(int(LAND_REF_FRACTION * int(n_px)) - int(margin), 1)


def elevation_at(edge, z_line, land_end=None):
    """
    Elevation of the edge (fractional column index, sea end first) on the profile z_line, linear
    between pixels. -> (z, above, below): NaN where the edge is NaN or off the profile; above /
    below flag edge samples landward / seaward of the profile's ends (an uprush above the profile
    is still known to be above every elevation the profile has). land_end (a column, line_end):
    edge samples at or past it are censored the same way (NaN, above), whatever the profile there:
    the line caps them, so their elevation would read low.
    """
    edge = np.asarray(edge, float)
    z_line = np.asarray(z_line, float)
    n = len(z_line)
    z = np.full(len(edge), np.nan)
    above = np.zeros(len(edge), bool)
    below = np.zeros(len(edge), bool)
    fin = np.flatnonzero(np.isfinite(z_line))
    ok = np.isfinite(edge)
    if not len(fin):
        return z, above, below
    e = np.clip(edge[ok], 0, n - 1)
    i0 = np.clip(np.floor(e).astype(int), 0, max(n - 2, 0))
    i1 = np.minimum(i0 + 1, n - 1)
    fr = e - i0
    za, zb = z_line[i0], z_line[i1]
    val = np.where(fr == 0, za, za + fr * (zb - za))
    idx = np.flatnonzero(ok)
    if land_end is not None:
        val = np.where(edge[idx] >= land_end, np.nan, val)
    z[idx] = val
    miss = idx[~np.isfinite(val)]
    above[miss] = edge[miss] > fin[-1]
    if land_end is not None:
        above[miss] |= edge[miss] >= land_end
    below[miss] = (edge[miss] < fin[0]) & ~above[miss]
    return z, above, below


# ---------------------------------------------------------------------
# Swash statistics (Stockdon et al. 2006 definitions)
# ---------------------------------------------------------------------

def detrend(x):
    """x minus its least-squares straight line in time (finite samples); NaNs kept."""
    x = np.asarray(x, float)
    t = np.arange(len(x), dtype=float)
    ok = np.isfinite(x)
    if ok.sum() < 2:
        return x.copy()
    b, a = np.polyfit(t[ok], x[ok], 1)
    return x - (a + b * t)


def ig_dof(n, fs=SAMPLE_HZ, cutoff=IG_CUTOFF_HZ):
    """Degrees of freedom of the infragravity variance of an n-sample record with a FLAT spectrum:
    2 per frequency bin in (0, cutoff) -- 58 for 10 min at 2 Hz. The most a record can have."""
    f = np.fft.rfftfreq(int(n), 1.0 / fs)
    return 2 * int(((f > 0) & (f < cutoff)).sum())


def ig_dof_eff(series, fs=SAMPLE_HZ, cutoff=IG_CUTOFF_HZ):
    """
    The effective degrees of freedom of the infragravity variance of this record, from its own
    periodogram P_k over the bins in (0, cutoff): nu = 2 (sum S_k)^2 / sum S_k^2 for the spectrum
    S_k, which a peaked or red swash spectrum makes far fewer than the flat-spectrum 2 per bin. The
    raw periodogram's own scatter (each bin a chi-square of 2 degrees of freedom) inflates sum P_k^2
    twofold and would halve it, so the moments are corrected: nu = 4 (sum P_k)^2 / sum P_k^2 - 2
    (2 per bin on average for a flat spectrum), within [2, ig_dof(n)]. The series is taken as
    band_sig takes it (gaps up to 5 s bridged, the mean removed). NaN when band_sig gives no value.
    """
    x = bridge_gaps(series)
    if x is None or len(x) < 64:
        return np.nan
    x = x - x.mean()
    f = np.fft.rfftfreq(len(x), 1.0 / fs)
    P = np.abs(np.fft.rfft(x)) ** 2
    p = P[(f > 0) & (f < cutoff)]
    if not len(p) or not (p ** 2).sum() > 0:
        return np.nan
    nu = 4.0 * p.sum() ** 2 / (p ** 2).sum() - 2.0
    return float(np.clip(nu, 2.0, ig_dof(len(x), fs, cutoff)))


def sig_sampling_rel(nu):
    """Relative sampling error of S_IG (4 x its standard deviation) for nu degrees of freedom of its
    variance: sqrt(2 / nu) / 2."""
    return float(np.sqrt(2.0 / nu) / 2.0) if np.isfinite(nu) and nu > 0 else np.nan


def swash_stats(z_abs, swl, edge, above, below, valid, min_peaks, min_valid, past_line=None):
    """
    z_abs: elevation of the swash edge (m NAVD88) per 2 Hz sample; swl: still-water level per
    sample (NaN if unknown); edge: the edge (fractional column index, for the uprush peaks);
    above / below: elevation_at's flags; past_line: the samples it censored at the line's landward
    end (for the note). -> dict setup_m (mean of z - swl), r2_m (98th percentile of the uprush
    maxima of z - swl, the maxima found on the edge series as runup_from_timestack.py finds them:
    peak_mask), s_inc_m, s_ig_m, s_m (4 x std of the detrended z - swl in > 0.05 Hz, < 0.05 Hz,
    all: band_sig), ig_dof (ig_dof_eff of that series), n_peaks, eta and sigma (mean and std of
    z_abs, for the foreshore band), missing (share of edge samples with no elevation: off the
    profile or past the line), notes.
    """
    z_abs, swl, edge = (np.asarray(a, float) for a in (z_abs, swl, edge))
    notes = []
    rel = z_abs - swl
    mask = peak_mask(edge)
    n_peaks = int(mask.sum())
    have_edge = np.isfinite(edge)
    missing = float((have_edge & ~np.isfinite(z_abs)).sum() / max(have_edge.sum(), 1))
    r2 = np.nan
    if np.isfinite(swl).any():
        pz = rel[mask].copy()
        pz[above[mask]] = np.inf
        pz[below[mask]] = -np.inf
        keep = ~np.isnan(pz)
        line_cut = (np.zeros(len(pz), bool) if past_line is None else np.asarray(past_line, bool)[mask])[keep]
        pz = pz[keep]
        if n_peaks >= min_peaks and valid >= min_valid and len(pz) >= min_peaks:
            with np.errstate(invalid="ignore"):
                r2 = float(np.percentile(pz, 98))
            if not np.isfinite(r2):
                r2 = np.nan
                notes.append("R2 above the line" if line_cut.any() else "R2 above the profile")
        setup = float(np.nanmean(rel)) if missing <= MISSING_ELEV_MAX and np.isfinite(rel).any() else np.nan
        if missing > MISSING_ELEV_MAX:
            notes.append(f"{missing:.0%} of the edge off the profile: no setup")
        series = rel
    else:
        setup = np.nan
        notes.append("no water level: no R2 or setup")
        series = z_abs
    zd = detrend(series)
    s_ig = float(band_sig(zd, 1e-9, IG_CUTOFF_HZ))
    s_inc = float(band_sig(zd, IG_CUTOFF_HZ, SAMPLE_HZ))
    s_all = float(band_sig(zd, 0.0, SAMPLE_HZ))
    dof = ig_dof_eff(zd) if np.isfinite(s_ig) else np.nan
    if not np.isfinite(s_ig):
        notes.append("edge gaps: no swash spectrum")
    elif missing > 0:
        notes.append(f"S_IG with {missing:.1%} of the edge off the profile (bridged: reads low)")
    fin = z_abs[np.isfinite(z_abs)]
    return {"setup_m": setup, "r2_m": r2, "s_inc_m": s_inc, "s_ig_m": s_ig, "s_m": s_all,
            "ig_dof": dof, "n_peaks": n_peaks, "eta": float(fin.mean()) if len(fin) else np.nan,
            "sigma": float(fin.std()) if len(fin) else np.nan, "missing": missing, "notes": notes}


def foreshore_beta(profile, dem, eta, sigma, fixed=None):
    """
    Stockdon's foreshore slope: the mean slope over <eta> +/- 2 sigma of the swash (their
    definition; setup_model.stockdon_frame uses the same band). The median DEM slope
    (setup_model.slope_grid) of the cells under the line's DEM pixels in that band, if they are at
    least MIN_SLOPE_CELLS distinct cells; else the band over the whole DEM (setup_model.
    ForeshoreSlope.band, >= MIN_BAND_CELLS cells); else `fixed` (--beta). A slope outside
    BETA_RANGE is no foreshore (a scarp, a gap edge). -> (beta or NaN, source).
    """
    lo, hi = eta - 2 * sigma, eta + 2 * sigma
    if profile["kind"] == "dem" and dem is not None and np.isfinite(lo):
        z, s, c = profile["z"], profile["slope"], profile["cell"]
        with np.errstate(invalid="ignore"):
            sel = np.isfinite(z) & ~profile["extrap"] & np.isfinite(s) & (z >= lo) & (z <= hi)
        cells = np.unique(c[sel])
        if len(cells) >= MIN_SLOPE_CELLS:
            b = float(np.median(s[sel]))
            if BETA_RANGE[0] <= b <= BETA_RANGE[1]:
                return b, f"dem under the line ({len(cells)} cells)"
        b, nc = dem.foreshore().band(lo, hi)
        if b is not None and BETA_RANGE[0] <= b <= BETA_RANGE[1]:
            return float(b), f"dem band, whole grid ({nc} cells)"
    if fixed is not None and np.isfinite(fixed):
        return float(fixed), "fixed (--beta)"
    return np.nan, "none"


def edge_contrast(gray, edge, near=CONTRAST_NEAR, far=CONTRAST_FAR, min_rows=10):
    """
    Michelson contrast across the swash tip, (foam - sand) / (foam + sand), at the uprush maxima
    (peak_mask of the edge: there the tip stands still and is a tip; in backwash the edge may be
    held behind the water by the retreat limit): foam the mean grey level `near`..`far` pixels
    seaward of the edge, sand the same landward; the median over those rows (gray time x pixels,
    sea end first, as the edge). Blur spreads the tip over many pixels and haze lifts both sides;
    either drives it toward 0. NaN with fewer than min_rows maxima clear of the line's ends.
    """
    g = np.asarray(gray, float)
    edge = np.asarray(edge, float)
    n = g.shape[1]
    rows = np.flatnonzero(peak_mask(edge))
    e = np.round(edge[rows]).astype(int)
    ok = (e - far >= 0) & (e + far < n)
    rows, e = rows[ok], e[ok]
    if len(rows) < min_rows:
        return np.nan
    c = np.concatenate([np.zeros((len(rows), 1)), np.cumsum(g[rows], axis=1)], axis=1)
    r = np.arange(len(rows))
    w = far - near + 1
    foam = (c[r, e - near + 1] - c[r, e - far]) / w
    sand = (c[r, e + far + 1] - c[r, e + near]) / w
    with np.errstate(invalid="ignore", divide="ignore"):
        m = (foam - sand) / (foam + sand)
    m = m[np.isfinite(m)]
    return float(np.median(m)) if len(m) else np.nan


def classify(brightness, contrast, window_sharpness, valid, n_peaks, have_profile, period,
             min_contrast=MIN_CONTRAST, min_valid=0.8, min_peaks=20):
    """The status of a stack and line, first match wins: dark / low contrast / no swash edge /
    no slope / no period / ok."""
    if not (np.isfinite(brightness) and brightness >= DARK_MEAN):
        return "dark"
    if (np.isfinite(window_sharpness) and window_sharpness < LENS_BLURRED) or \
            (np.isfinite(contrast) and contrast < min_contrast):
        return "low contrast"
    if not (valid >= min_valid and n_peaks >= min_peaks):
        return "no swash edge"
    if not have_profile:
        return "no slope"
    if not (np.isfinite(period) and period > 0):
        return "no period"
    return "ok"


# ---------------------------------------------------------------------
# Inputs: period, water level, the bright frame's sharpness
# ---------------------------------------------------------------------

def _num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return np.nan
    return v


def _read_rows(path):
    p = Path(path) if path else None
    if p is None or not p.exists():
        return None
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def _nearest(ep, val, t, max_gap):
    if ep is None or not len(ep):
        return np.nan
    i = int(np.argmin(np.abs(ep - t)))
    return float(val[i]) if abs(ep[i] - t) <= max_gap and np.isfinite(val[i]) else np.nan


class Periods:
    """
    T0 for a stack. 'camera': the median Tp of the stack's lines whose status is ok in
    archive/wave_period_<cam>.csv (timestack_wave_period.py, the same stacks). 'buoy': the hourly
    buoy peak period of archive/waves_marconi.csv (tp_s; marconi_waves.py: NDBC 44008 DPD), else
    archive/waves_44008.csv dpd_s, within 1.5 h. 'auto': camera, else buoy. for_stack -> (T0,
    source, pending): pending when the camera file has no row for the stack yet (it runs hourly).
    """

    def __init__(self, camera_csv=None, waves_csv=None, buoy_csv=None, mode="auto"):
        self.mode = mode
        self.cam = None
        rows = _read_rows(camera_csv)
        if rows is not None:
            self.cam = {}
            for r in rows:
                self.cam.setdefault(r.get("filename", ""), []).append((r.get("status", ""), _num(r.get("tp_s"))))
        self.buoys = []
        for path, col, label in ((waves_csv, "tp_s", "buoy Tp (waves_marconi)"),
                                 (buoy_csv, "dpd_s", "buoy 44008 DPD")):
            rows = _read_rows(path)
            if rows:
                ep = np.array([_num(r.get("epoch")) for r in rows])
                tp = np.array([_num(r.get(col)) for r in rows])
                ok = np.isfinite(ep)
                self.buoys.append((ep[ok], tp[ok], label))

    def for_stack(self, filename, mid_epoch):
        pending = False
        if self.mode in ("auto", "camera"):
            if self.cam is not None and filename not in self.cam:
                pending = True                     # the period script has not reached this stack yet
            elif self.cam is not None:
                tps = [tp for st, tp in self.cam[filename] if st == "ok" and np.isfinite(tp) and tp > 0]
                if tps:
                    return float(np.median(tps)), "camera Tp", False
            if self.mode == "camera":
                return np.nan, "", pending
        for ep, tp, label in self.buoys:
            t = _nearest(ep, tp, float(mid_epoch), PERIOD_MAX_GAP_S)
            if np.isfinite(t) and t > 0:
                return t, label, pending
        return np.nan, "", pending


class LevelSeries:
    """A water-level record (epoch, level_navd88 -- or the formats historical_forcing.read_level_csv
    reads) as a level source with WaterLevel's .at(): linear, NaN across gaps over 1 h."""

    def __init__(self, path, max_gap_s=3600.0, label="level csv"):
        from historical_forcing import read_level_csv
        self.ep, self.lv = read_level_csv(path)
        self.max_gap, self.label = max_gap_s, label

    def at(self, epochs):
        t = np.atleast_1d(np.asarray(epochs, float))
        out = np.full(len(t), np.nan)
        if len(self.ep) >= 2:
            i = np.clip(np.searchsorted(self.ep, t), 1, len(self.ep) - 1)
            ok = (t >= self.ep[0]) & (t <= self.ep[-1]) & (self.ep[i] - self.ep[i - 1] <= self.max_gap)
            out[ok] = np.interp(t[ok], self.ep, self.lv)
        return out, np.where(np.isfinite(out), self.label, "").astype(object)


class WindowSharpness:
    """owg_live.py's sharpness of the c2 bright frame of the same collection (archive/owg_c2_H.csv):
    the window check of station_status.py, per stack."""

    def __init__(self, path):
        rows = _read_rows(path)
        self.ep = self.sh = None
        if rows:
            self.ep = np.array([_num(r.get("epoch")) for r in rows])
            self.sh = np.array([_num(r.get("sharpness")) for r in rows])

    def at(self, epoch):
        return _nearest(self.ep, self.sh, float(epoch), OWG_MAX_GAP_S)


# ---------------------------------------------------------------------
# One stack and line
# ---------------------------------------------------------------------

def _r(x, nd):
    return "" if x is None or not np.isfinite(x) else round(float(x), nd)


def plane_stretch(dist, d_beach, edge):
    """How many metres of path along the line on the still-water plane (runup_from_timestack's dist,
    the sum of its pixel-to-pixel steps) one metre on the beach is, over the pixels the swash edge
    covers (its 5-95th percentile). Two things make it: the plane geometry (c2 sees the swash at a
    grazing angle from 19.8 m up, so a point above still water lands further along its ray on that
    plane) and, on a line drawn as single-pixel steps, the staircase's zigzag; both are measured on
    the actual line. 1.0 where it cannot be measured."""
    dist, d_beach, edge = (np.asarray(a, float) for a in (dist, d_beach, edge))
    e = edge[np.isfinite(edge)]
    if len(e) < 10:
        return 1.0
    lo, hi = np.percentile(e, [5, 95])
    i = np.arange(int(np.floor(lo)), int(np.ceil(hi)) + 1)
    i = i[(i >= 0) & (i < len(d_beach))]
    i = i[np.isfinite(d_beach[i]) & np.isfinite(dist[i])]
    if len(i) < 10 or d_beach[i[-1]] - d_beach[i[0]] < 0.5:
        return 1.0
    r = (dist[i[-1]] - dist[i[0]]) / (d_beach[i[-1]] - d_beach[i[0]])
    return float(np.clip(r, 0.2, 20.0)) if np.isfinite(r) and r > 0 else 1.0


def sig_extra_rel(st, n):
    """S_IG's own relative error beyond Stockdon's scatter: its sampling error (sig_sampling_rel of
    the record's effective degrees of freedom, st['ig_dof']; the flat-spectrum count of an n-sample
    record where that is unknown) and the share of the edge off the profile (st['missing']: bridged
    samples flatten the uprush tops, so S_IG reads low by about that much), in quadrature."""
    samp = sig_sampling_rel(st["ig_dof"])
    return float(np.hypot(samp if np.isfinite(samp) else sig_sampling_rel(ig_dof(n)), st["missing"]))


def late_inputs(pending, have_profile, t0, wl_mid, wl_src, want_gps):
    """The inputs that arrive late and are missing for a stack: the camera period (pending), a DEM,
    any period, the water level, and the GPS level where the GPS is the level source in use."""
    miss = []
    if pending:
        miss.append("camera period")
    if not have_profile:
        miss.append("DEM")
    if not (np.isfinite(t0) and t0 > 0):
        miss.append("period")
    if not np.isfinite(wl_mid):
        miss.append("water level")
    elif want_gps and wl_src != "gps":
        miss.append("GPS level")
    return miss


def invert_all(st, beta, t0, dem_level, regime):
    """The three H0 estimates of one stack in one regime. R2 and setup only on a DEM (dem_level: a
    planar profile does not know where still water meets the beach) and, but in the dissipative
    regime, with a foreshore slope."""
    est = {"sig": h0_from_sig(st["s_ig_m"], t0, regime), "r2": np.nan, "setup": np.nan}
    if dem_level and (np.isfinite(beta) or regime == "dissipative"):
        est["r2"] = h0_from_r2(st["r2_m"], beta, t0, regime)
        est["setup"] = h0_from_setup(st["setup_m"], beta, t0, regime)
    return est


def misfit(est, h0, t0, beta, scatter=None, sig_extra=0.0, regime="intermediate"):
    """How badly the estimates agree with their combination: sum of ((ln H0_k - ln H0) / sigma_k)^2,
    sigma_k at the combined H0 (0 for a single estimate)."""
    out = 0.0
    for k, v in est.items():
        if v is not None and np.isfinite(v) and v > 0 and np.isfinite(h0) and h0 > 0:
            sd = estimate_sigma(k, h0, t0, beta, scatter, sig_extra, regime)
            if 0 < sd < np.inf:
                out += ((np.log(v) - np.log(h0)) / sd) ** 2
    return float(out)


def estimate_h0(st, beta, t0, dem_level, scatter=None, sig_extra=0.0):
    """
    H0 from one stack's statistics (swash_stats' dict), inverted and combined with Stockdon's
    intermediate forms (Eqs. 10, 12, 19) and with the dissipative ones (Eq. 18 and the setup and
    S_IG beside it), each giving its own xi0. The regime is the one whose own H0 agrees with it
    (intermediate xi0 >= 0.3, dissipative below); when both do or neither does (H0 near the
    threshold), the one whose estimates agree better (misfit; the intermediate on a tie), noted as
    unsure. Without a foreshore slope, the intermediate. -> dict est (the three H0s), h0, sigma (of
    ln H0), used, regime, xi0 (at the H0 kept), note.
    """
    b_comb = beta if np.isfinite(beta) else 0.1     # S_IG alone: the slope does not enter
    res = {}
    for regime in ("intermediate", "dissipative"):
        est = invert_all(st, beta, t0, dem_level, regime)
        h0, sig_ln, used = combine(est, t0, b_comb, scatter, sig_extra, regime)
        res[regime] = {"est": est, "h0": h0, "sigma": sig_ln, "used": used, "regime": regime,
                       "xi0": iribarren(beta, h0, t0), "note": "",
                       "misfit": misfit(est, h0, t0, b_comb, scatter, sig_extra, regime)}
    i, d = res["intermediate"], res["dissipative"]
    if not (np.isfinite(i["xi0"]) and np.isfinite(d["xi0"])):
        return i
    ok_i, ok_d = i["xi0"] >= XI_DISSIPATIVE, d["xi0"] < XI_DISSIPATIVE
    if ok_i and not ok_d:
        return i
    pick = d if (ok_d and not ok_i) or d["misfit"] < i["misfit"] else i
    unsure = ok_i == ok_d
    if pick is d or unsure:
        pick["note"] = (f"xi0 {i['xi0']:.2f} with the intermediate forms, {d['xi0']:.2f} with the dissipative "
                        f"(Stockdon's, below xi0 {XI_DISSIPATIVE:g}: R2 0.043, setup 0.016, S_IG 0.046 x sqrt(H0 L0)): "
                        f"{pick['regime']} kept" + ("; regime unsure" if unsure else ""))
    return pick


def estimate_stack(ras_path, line, ctx, now=None):
    """The row for one stack and line, None if the stack cannot be read, or 'wait' when an input that
    arrives late (camera period, DEM, water level, the GPS level) is missing, the stack is younger
    than ctx.wait_days and its numbers would be computed (not dark, an edge found)."""
    rargs = copy.copy(ctx.rargs)
    rargs.line = int(line)
    # georectified (for the retreat limit) on the plane at still water, as the cron's runup is
    e0 = rft.epoch_from_name(Path(ras_path).name)
    if e0 is not None:
        z0 = ctx.wl.at([e0 + rft.BURST_MID_OFFSET_S])[0][0]
        rargs.z_plane = float(z0) if np.isfinite(z0) else ctx.rargs.z_plane
    try:
        s = rft.swash_edge(ras_path, rargs, ctx.pix_cache, None)
    except SystemExit as exc:                      # a line the stack's .pix does not have
        print(f"  {Path(ras_path).name}: {exc} -- skipped")
        return None
    if s is None:
        return None
    io, eo = ctx.pix_cache[s["cam"]]
    mid = float(s["mid_epoch"])
    dem = ctx.dems.for_epoch(mid) if ctx.dems is not None else None
    key = (dem.name if dem is not None else "", ctx.beta, s["line_uv"].tobytes())
    if key not in ctx.profiles:
        prof = line_profile(s["line_uv"][:, 0], s["line_uv"][:, 1], io, eo, dem=dem,
                            max_extrapolate=ctx.max_extrapolate)
        if not np.isfinite(prof["z"]).any() and ctx.beta is not None:
            prof = line_profile(s["line_uv"][:, 0], s["line_uv"][:, 1], io, eo, beta=ctx.beta)
        ctx.profiles[key] = prof
    prof = ctx.profiles[key]
    have_profile = bool(np.isfinite(prof["z"]).any())
    # the retreat limit, on the beach: the runup script applies it along its still-water-plane path
    stretch = plane_stretch(s["dist"], prof["d"], s["edge"]) if have_profile else 1.0
    if abs(stretch - 1.0) > 0.02:
        rargs.max_retreat_speed = ctx.retreat_speed * stretch
        s = rft.swash_edge(ras_path, rargs, ctx.pix_cache, None)
    gray, edge = s["gray"], s["edge"]
    n = len(edge)
    brightness = float(gray.mean())
    contrast = edge_contrast(gray, edge)
    sharp = ctx.window.at(s["epoch"]) if ctx.window is not None else np.nan
    valid = float(s["valid"])
    n_peaks = int(peak_mask(edge).sum())
    t0, psrc, pending = ctx.periods.for_stack(s["name"], mid)

    t = s["epoch"] + np.arange(n) / SAMPLE_HZ
    lv, src = ctx.wl.at([t[0], mid, t[-1]])
    if np.isfinite(lv).all():
        swl = np.interp(t, [t[0], mid, t[-1]], lv)
        wl_mid, wl_src = float(lv[1]), str(src[1])
    elif np.isfinite(lv[1]):
        swl, wl_mid, wl_src = np.full(n, float(lv[1])), float(lv[1]), str(src[1])
    else:
        swl, wl_mid, wl_src = np.full(n, np.nan), np.nan, ""

    status = classify(brightness, contrast, sharp, valid, n_peaks, have_profile, t0,
                      ctx.min_contrast, rargs.min_valid, rargs.min_peaks)
    now = datetime.now(timezone.utc).timestamp() if now is None else float(now)
    young = now - mid < ctx.wait_days * 86400.0
    if young and status not in ("dark", "no swash edge") and \
            late_inputs(pending, have_profile, t0, wl_mid, wl_src, getattr(ctx, "want_gps", False)):
        return "wait"

    row = {"epoch": int(mid), "time_utc": datetime.fromtimestamp(mid, tz=timezone.utc).isoformat(),
           "filename": s["name"], "line": int(line), "status": status, "note": "",
           "brightness": _r(brightness, 1), "edge_contrast": _r(contrast, 3),
           "window_sharpness": _r(sharp, 2), "valid_fraction": _r(valid, 3), "n_peaks": n_peaks,
           "sea_at": s["sea_at"],
           "profile": (f"{prof['kind']}:{prof['name']}" if prof["kind"] != "none" else ""),
           "water_level_m": _r(wl_mid, 3), "water_level_source": wl_src,
           "period_s": _r(t0, 2), "period_source": psrc}
    notes = []
    if status in ("ok", "low contrast") and have_profile and np.isfinite(t0):
        n_px = len(prof["z"])
        cut = line_end(n_px)
        z_abs, above, below = elevation_at(edge, prof["z"], cut)
        with np.errstate(invalid="ignore"):
            past_line = np.isfinite(edge) & (edge >= cut)
        st = swash_stats(z_abs, swl, edge, above, below, valid, rargs.min_peaks, rargs.min_valid, past_line)
        notes += st["notes"]
        with np.errstate(invalid="ignore"):
            on = np.isfinite(z_abs)
            xf = float(prof["extrap"][np.clip(np.round(edge[on]).astype(int), 0, n_px - 1)].mean()) \
                if on.any() else np.nan
        beta, bsrc = foreshore_beta(prof, dem if prof["kind"] == "dem" else None, st["eta"], st["sigma"],
                                    ctx.beta)
        dem_level = prof["kind"] == "dem"           # R2 and setup need the beach's true level
        if not np.isfinite(beta):
            notes.append("no foreshore slope: no R2 or setup estimate")
        if not dem_level:
            notes.append("planar profile: R2 and setup not inverted")
        sig_extra = sig_extra_rel(st, n)
        res = estimate_h0(st, beta, t0, dem_level, ctx.scatter, sig_extra)
        est, h0, sig_ln, used, xi = res["est"], res["h0"], res["sigma"], res["used"], res["xi0"]
        if res["note"]:
            notes.append(res["note"])
        ks = float(shoaling_coefficient(t0, DESHOAL_DEPTH_M))
        row.update({"z_extrapolated_fraction": _r(xf, 3), "off_profile_fraction": _r(st["missing"], 3),
                    "beta_f": _r(beta, 4), "beta_source": bsrc,
                    "iribarren": _r(xi, 2) if np.isfinite(beta) else "",
                    "setup_m": _r(st["setup_m"], 3), "r2_m": _r(st["r2_m"], 3),
                    "s_inc_m": _r(st["s_inc_m"], 3), "s_ig_m": _r(st["s_ig_m"], 3), "s_m": _r(st["s_m"], 3),
                    "h0_r2_m": _r(est["r2"], 3), "h0_sig_m": _r(est["sig"], 3),
                    "h0_setup_m": _r(est["setup"], 3), "h0_m": _r(h0, 3), "h0_sigma_rel": _r(sig_ln, 3),
                    "combined_from": used, "ks_21m": _r(ks, 4), "hs_m": _r(h0 * ks, 3),
                    "hs_sig_m": _r(est["sig"] * ks, 3)})
        if status == "ok" and not np.isfinite(h0):
            row["status"] = "no swash edge"
    if status == "low contrast":
        why = ([f"window sharpness {sharp:.1f} (< {LENS_BLURRED:g})"] if np.isfinite(sharp) and sharp < LENS_BLURRED
               else []) + ([f"edge contrast {contrast:.3f} (< {ctx.min_contrast:g})"]
                           if np.isfinite(contrast) and contrast < ctx.min_contrast else [])
        notes.insert(0, "; ".join(why))
    row["note"] = "; ".join(notes)
    return row


# ---------------------------------------------------------------------
# The output file: appended, never rewritten except by --reprocess
# ---------------------------------------------------------------------

def read_done(path):
    """{(filename, line)} already in the output file (empty when there is none)."""
    rows = _read_rows(path)
    return set() if not rows else {(r["filename"], str(r["line"])) for r in rows}


def check_header(path, reprocess=False):
    """A non-empty output file whose header is not FIELDS is refused (SystemExit) unless reprocess,
    so a column change never shifts old rows. -> whether the file exists with content."""
    path = Path(path)
    exists = path.exists() and path.stat().st_size > 0
    if exists:
        with open(path, newline="") as f:
            header = next(csv.reader(f), [])
        if header != FIELDS and not reprocess:
            raise SystemExit(f"{path}: its columns are not this version's; rerun with --reprocess")
    return exists


def write_rows(path, rows, reprocess=False):
    """
    Appends rows to the CSV (header first if new). Rows whose (filename, line) are already in it
    are not appended again. With reprocess, the file is rewritten: its rows that these rows redo
    are replaced, the others kept, all in (epoch, line) order. A file with another header is
    refused (check_header): rerun with --reprocess. -> the number of rows written.
    """
    path = Path(path)
    if not rows:
        return 0
    rows = sorted(rows, key=lambda r: (int(r["epoch"]), int(r["line"])))
    exists = check_header(path, reprocess)
    if reprocess and exists:
        redo = {(r["filename"], str(r["line"])) for r in rows}
        kept = [r for r in (_read_rows(path) or []) if (r.get("filename"), str(r.get("line"))) not in redo]
        merged = sorted(kept + rows, key=lambda r: (int(float(r["epoch"])), int(float(r["line"]))))
        tmp = path.with_suffix(path.suffix + ".tmp")
        with open(tmp, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS, restval="", extrasaction="ignore")
            w.writeheader()
            w.writerows(merged)
        tmp.replace(path)
        return len(rows)
    done = read_done(path) if exists else set()
    new = [r for r in rows if (r["filename"], str(r["line"])) not in done]
    seen, uniq = set(), []
    for r in new:
        k = (r["filename"], str(r["line"]))
        if k not in seen:
            seen.add(k)
            uniq.append(r)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, restval="", extrasaction="ignore")
        if not exists:
            w.writeheader()
        w.writerows(uniq)
    return len(uniq)


# ---------------------------------------------------------------------
# Self-test: a synthetic swash series with a known H0
# ---------------------------------------------------------------------

def synthetic_swash(h0, t0, beta, n=1200, seed=0):
    """
    A 2 Hz swash elevation series (m above still water) whose setup, infragravity swash and R2 are
    what Stockdon et al. (2006) give for (h0, t0, beta): mean exactly Eq. 10; an infragravity part
    (random phases, flat 0.005-0.05 Hz) scaled so its detrended 4-sigma is exactly Eq. 12; an
    incident part (Gaussian spectrum about 1/t0) scaled, by bisection, until the 98th percentile of
    the uprush maxima is Eq. 19. -> (series, targets dict). ValueError when Eq. 19's R2 is below
    what the infragravity swash alone reaches (very low slopes: Eq. 19 is not for them).
    """
    rng = np.random.default_rng(seed)
    f = np.fft.rfftfreq(n, 1.0 / SAMPLE_HZ)
    eta = float(stockdon_setup(h0, t0, beta))
    s_ig = float(stockdon_sig(h0, t0))
    r2 = float(stockdon_r2(h0, t0, beta))

    def part(amp):
        X = amp * np.exp(2j * np.pi * rng.random(len(f)))
        return np.fft.irfft(X, n)

    ig = part(((f >= 0.005) & (f < IG_CUTOFF_HZ)).astype(float))
    inc = part(np.exp(-0.5 * ((f - 1.0 / t0) / (0.2 / t0)) ** 2) * (f >= IG_CUTOFF_HZ))
    inc /= np.std(inc)

    def series(k):
        a = 1.0
        for _ in range(4):                  # detrending is linear: a few passes make S_IG exact
            a *= s_ig / band_sig(detrend(a * ig + k * inc), 1e-9, IG_CUTOFF_HZ)
        x = a * ig + k * inc
        return eta + x - x.mean()

    def r2_of(k):
        z = series(k)
        return float(np.percentile(z[peak_mask(z)], 98)), z

    lo, hi = 0.0, 4.0 * s_ig + 1.0
    if r2_of(lo)[0] > r2:
        raise ValueError(f"R2 {r2:.2f} m is below the infragravity swash alone")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if r2_of(mid)[0] < r2:
            lo = mid
        else:
            hi = mid
    z = r2_of(0.5 * (lo + hi))[1]
    return z, {"setup": eta, "s_ig": s_ig, "r2": r2}


def self_test():
    """Stockdon-generated swash on a curved synthetic profile, with a rising tide: H0 recovered."""
    n_px = 900
    p = np.arange(n_px, dtype=float)
    z_line = -0.8 + 3.6 * (p / (n_px - 1)) ** 1.3          # sea end first, uneven metres per pixel
    worst = {"exact": 0.0, "noisy": 0.0}
    rng = np.random.default_rng(3)
    for h0, t0, beta in ((0.8, 8.0, 0.10), (1.5, 10.0, 0.12), (3.0, 12.0, 0.09), (2.2, 7.0, 0.14)):
        z_rel, tg = synthetic_swash(h0, t0, beta)
        tt = np.arange(len(z_rel)) / SAMPLE_HZ
        swl = 0.3 + 0.12 * tt / tt[-1]                     # rising 0.12 m over the 10 minutes
        for case in ("exact", "noisy"):
            z_abs = swl + z_rel
            if case == "noisy":                             # 3 cm noise, 4% of rows lost in short runs
                z_abs = z_abs + rng.normal(0.0, 0.03, len(z_abs))
                for s0 in rng.choice(len(z_abs) - 4, 12, replace=False):
                    z_abs[s0:s0 + 4] = np.nan
            edge = np.interp(z_abs, z_line, p)              # what the camera would see
            edge[~np.isfinite(z_abs)] = np.nan
            zz, above, below = elevation_at(edge, z_line)
            st = swash_stats(zz, swl, edge, above, below, float(np.isfinite(edge).mean()), 20, 0.8)
            est = {"r2": h0_from_r2(st["r2_m"], beta, t0), "sig": h0_from_sig(st["s_ig_m"], t0),
                   "setup": h0_from_setup(st["setup_m"], beta, t0)}
            h, sl, used = combine(est, t0, beta)
            errs = [v / h0 - 1 for v in est.values()] + [h / h0 - 1]
            worst[case] = max(worst[case], max(abs(e) for e in errs))
            print(f"  {case:<5} H0 {h0:.1f} m, T0 {t0:4.1f} s, beta {beta:.2f}: R2 {st['r2_m']:.3f} "
                  f"(Eq. 19 {tg['r2']:.3f}), S_IG {st['s_ig_m']:.3f} ({tg['s_ig']:.3f}), setup "
                  f"{st['setup_m']:.3f} ({tg['setup']:.3f}) m -> H0 "
                  + ", ".join(f"{k} {est[k]:.2f}" for k in ("r2", "sig", "setup"))
                  + f", combined {h:.3f} m (+/-{sl:.0%})")
    ok = worst["exact"] < 0.02 and worst["noisy"] < 0.05
    print(("PASS" if ok else "FAIL") + f": worst H0 error {worst['exact']:.2%} exact (tolerance 2%), "
          f"{worst['noisy']:.1%} with noise and gaps (tolerance 5%)")
    return 0 if ok else 1


# ---------------------------------------------------------------------
# Running over the archive
# ---------------------------------------------------------------------

class Context:
    pass


def make_context(args, wl=None):
    ctx = Context()
    rargs = rft.build_parser().parse_args(["-", "--output", "-", "--camera", args.camera])
    rargs.pix, rargs.io, rargs.eo, rargs.plot_dir = args.pix, args.io, args.eo, None
    ctx.rargs = rargs
    if args.max_retreat_speed is None:
        args.max_retreat_speed = rargs.max_retreat_speed
    ctx.pix_cache, ctx.profiles = {}, {}
    ctx.dems = DemSource(args.dem) if args.dem and Path(args.dem).exists() else None
    ctx.beta = args.beta
    ctx.max_extrapolate = args.max_extrapolate
    ctx.retreat_speed = args.max_retreat_speed
    ctx.min_contrast = args.min_contrast
    ctx.wait_days = args.wait_days
    ctx.scatter = parse_scatter(args.scatter)
    ctx.periods = Periods(args.period_csv, args.waves_csv, args.buoy_csv, args.period_source)
    ctx.window = WindowSharpness(args.owg_csv)
    if wl is not None:
        ctx.wl = wl
    elif args.level_csv:
        ctx.wl = LevelSeries(args.level_csv)
    else:
        from marconi_water_level import WaterLevel
        ctx.wl = WaterLevel(None if args.no_gps else args.gnssr_spline, args.gauge_csv)
    # young stacks wait for the GPS level where the GPS is in use (it arrives 1-2 days late; the
    # cron's runup stage waits for it too, --require-water-level)
    ctx.want_gps = getattr(ctx.wl, "s_ep", None) is not None
    return ctx


def _day_epoch(s, end=False):
    t = datetime.strptime(s, "%Y-%m-%d").replace(tzinfo=timezone.utc).timestamp()
    return t + 86400.0 if end else t


def run(args, wl=None, now=None):
    out = Path(args.output or HERE / "archive" / f"hs_timestack_{args.camera}.csv")
    check_header(out, args.reprocess)               # refuse a file of other columns before any work
    ctx = make_context(args, wl)
    done = set() if args.reprocess else read_done(out)
    paths = sorted(Path(args.ras_dir).glob(f"*.{args.camera}.ras.tiff"), key=lambda q: q.name)
    lo = _day_epoch(args.since) if args.since else -np.inf
    hi = _day_epoch(args.until, end=True) if args.until else np.inf
    print(f"timestack wave height {args.camera}: {len(paths)} stack(s) in {args.ras_dir}; DEM: "
          + (ctx.dems.ev.what if ctx.dems is not None else f"none at {args.dem}")
          + (f"; planar beta {ctx.beta:g} where no DEM" if ctx.beta is not None else
             ("" if ctx.dems is not None else " (every line 'no slope' without --beta)")))
    if hasattr(ctx.wl, "describe"):
        print("  " + ctx.wl.describe())
    chunk, waiting, n, counts, latest = [], 0, 0, {}, None
    try:
        for p in paths:
            e = rft.epoch_from_name(p.name)
            if e is None or not lo <= e < hi:
                continue
            for ln in args.lines:
                if (p.name, str(ln)) in done:
                    continue
                r = estimate_stack(p, ln, ctx, now=now)
                if r is None:
                    break
                if r == "wait":
                    waiting += 1
                    continue
                chunk.append(r)
                counts[r["status"]] = counts.get(r["status"], 0) + 1
                if r["status"] == "ok" and r["hs_m"] != "":
                    latest = r
                if len(chunk) >= CHUNK_ROWS:            # appended as it goes: a crash loses at most these
                    n += write_rows(out, chunk, reprocess=args.reprocess)
                    chunk = []
    finally:
        if chunk:
            n += write_rows(out, chunk, reprocess=args.reprocess)
    print(f"timestack wave height {args.camera}: {n} new row(s) -> {out}"
          + (" (" + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())) + ")" if counts else "")
          + (f"; {waiting} line-stack(s) wait for a period, DEM, water level or the GPS level" if waiting else "")
          + (f"; latest Hs {latest['hs_m']} m ({latest['time_utc'][:16]}Z)" if latest else ""))
    return 0


# ---------------------------------------------------------------------
# --compare
# ---------------------------------------------------------------------

def scores(est, ref):
    """(n, bias, RMSE, r) of est - ref over the pairs where both are finite."""
    est, ref = np.asarray(est, float), np.asarray(ref, float)
    ok = np.isfinite(est) & np.isfinite(ref)
    n = int(ok.sum())
    if n == 0:
        return 0, np.nan, np.nan, np.nan
    d = est[ok] - ref[ok]
    r = float(np.corrcoef(est[ok], ref[ok])[0, 1]) if n > 2 and np.std(est[ok]) > 0 and np.std(ref[ok]) > 0 else np.nan
    return n, float(d.mean()), float(np.sqrt(np.mean(d ** 2))), r


def with_owg(est, owg, ref):
    """The estimate and the OWG on the SAME stacks (all three finite): -> dict n, bias, rmse (the
    estimate's), owg_bias, owg_rmse, err_corr (of the two errors), blend_rmse (their 50:50 mean)."""
    est, owg, ref = (np.asarray(a, float) for a in (est, owg, ref))
    ok = np.isfinite(est) & np.isfinite(owg) & np.isfinite(ref)
    out = {"n": int(ok.sum())}
    if not out["n"]:
        return out
    e1, e2 = est[ok] - ref[ok], owg[ok] - ref[ok]
    eb = 0.5 * (est[ok] + owg[ok]) - ref[ok]
    out.update(bias=float(e1.mean()), rmse=float(np.sqrt(np.mean(e1 ** 2))), owg_bias=float(e2.mean()),
               owg_rmse=float(np.sqrt(np.mean(e2 ** 2))), blend_rmse=float(np.sqrt(np.mean(eb ** 2))),
               err_corr=(float(np.corrcoef(e1, e2)[0, 1]) if out["n"] > 2 and np.std(e1) > 0 and np.std(e2) > 0
                         else np.nan))
    return out


HS_BANDS = ((0.0, 1.0), (1.0, 2.0), (2.0, 99.0))
STORM_HS_M = 2.0                  # a storm stack: the reference Hs 2 m or more (VALIDATION 4.11)
STORM_GAP_S = 86400.0             # storm stacks 24 h or more apart are separate storms (a night is shorter)


def compare_subsets(df, ref_col, by_source=True):
    """[(name, mask)]: all clear (status ok), clear by Hs band of the reference, clear by water-level
    source and by period source (when they differ), and blurred (status low contrast)."""
    clear = df["status"] == "ok"
    out = [("all clear", clear)]
    for a, b in HS_BANDS:
        out.append((f"clear, {a:g}-{b:g} m" if b < 99 else f"clear, >= {a:g} m",
                    clear & (df[ref_col] >= a) & (df[ref_col] < b)))
    if by_source:
        for col, short in (("water_level_source", "level"), ("period_kind", "T0")):
            vals = sorted(v for v in df.loc[clear, col].dropna().unique() if str(v))
            if len(vals) > 1:
                out += [(f"clear, {short} {v}", clear & (df[col] == v)) for v in vals]
    out.append(("blurred", df["status"] == "low contrast"))
    return out


def storm_count(epochs, gap_s=STORM_GAP_S):
    """Separate storms among storm-stack times: a new one wherever 24 h or more pass without one."""
    t = np.sort(np.asarray(epochs, float))
    return int(len(t) > 0) + int((np.diff(t) >= gap_s).sum()) if len(t) else 0


def compare_counts(df, ref_col, col="hs_sig"):
    """The counts of decision rule 1 for an estimate column against a reference: clear stacks with
    both, distinct UTC days, storm stacks (reference >= 2 m) and separate storms."""
    m = (df["status"] == "ok") & df[col].notna() & df[ref_col].notna()
    st = m & (df[ref_col] >= STORM_HS_M)
    days = {datetime.fromtimestamp(float(t), tz=timezone.utc).strftime("%Y-%m-%d") for t in df.loc[m, "start"]}
    return {"stacks": int(m.sum()), "days": len(days), "storm_stacks": int(st.sum()),
            "storms": storm_count(df.loc[st, "start"])}


def compare_table(df, ref_col, ref_name, cols, owg_col=None, result=None):
    """
    Prints, and stores in result[ref_col] when given: (1) each estimate column against ref_col on
    its own stacks (n / bias / RMSE / r) for compare_subsets; (2) with owg_col, each estimate and the
    OWG on the stacks where the estimate, the OWG and the reference all exist (the estimate's bias
    and RMSE there, the OWG's, the correlation of their errors and the RMSE of their 50:50 mean):
    all clear, by Hs band, blurred.
    """
    print(f"\n  against {ref_name}:")
    if any("[not indep.]" in label for label, _ in cols):
        print(f"    [not indep.]: not independent of {ref_col} (the DEM's level carries setup computed from it)")
    print(f"    {'estimate':<34} {'subset':<20} {'n':>5} {'bias m':>8} {'RMSE m':>8} {'r':>6}")
    alone, paired = [], []
    for label, col in cols:
        for name, m in compare_subsets(df, ref_col):
            n, bias, rmse, r = scores(df.loc[m, col], df.loc[m, ref_col])
            if n:
                alone.append({"estimate": label, "col": col, "subset": name, "n": n, "bias": bias,
                              "rmse": rmse, "r": r})
                print(f"    {label:<34} {name:<20} {n:>5} {bias:>+8.3f} {rmse:>8.3f} {r:>6.2f}")
    if owg_col is not None and owg_col in df and df[owg_col].notna().any():
        print(f"    with the OWG, on the stacks where the estimate, the OWG and {ref_col} all exist:")
        print(f"    {'estimate':<34} {'subset':<20} {'n':>5} {'bias m':>8} {'RMSE m':>8} {'OWG bias':>9} "
              f"{'OWG RMSE':>9} {'err corr':>9} {'50:50 RMSE':>11}")
        for label, col in cols:
            for name, m in compare_subsets(df, ref_col, by_source=False):
                w = with_owg(df.loc[m, col], df.loc[m, owg_col], df.loc[m, ref_col])
                if w["n"]:
                    w.update(estimate=label, col=col, subset=name)
                    paired.append(w)
                    print(f"    {label:<34} {name:<20} {w['n']:>5} {w['bias']:>+8.3f} {w['rmse']:>8.3f} "
                          f"{w['owg_bias']:>+9.3f} {w['owg_rmse']:>9.3f} {w['err_corr']:>+9.2f} "
                          f"{w['blend_rmse']:>11.3f}")
    c = compare_counts(df, ref_col)
    print(f"    rule 1 counts (from S_IG, clear): {c['stacks']} stack(s) on {c['days']} day(s); "
          f"{c['storm_stacks']} storm stack(s) ({ref_col} >= {STORM_HS_M:g} m) in {c['storms']} storm(s) "
          f"(24 h without one starts a new storm)")
    if result is not None:
        result[ref_col] = {"alone": alone, "with_owg": paired, "counts": c}


def load_compare(args):
    """The output file as --compare scores it: one row per stack and status (the median of its lines),
    with hs_sig / hs_r2 / hs_setup (each estimate x Ks, R2 and setup only from DEM profiles), the
    stack's start (its filename epoch), only stacks inside --since / --until (as the run selects
    them), and the references: the OWG's ok frame of the same collection, the ADCP within 30 min,
    and hs_best / hs_buoy_marconi of the hour the stack starts in (as marconi_waves.py files frames).
    -> (DataFrame, whether any stack falls in the ADCP record)."""
    import pandas as pd
    out = Path(args.output or HERE / "archive" / f"hs_timestack_{args.camera}.csv")
    if not out.exists():
        sys.exit(f"no {out} yet: run without --compare first")
    d = pd.read_csv(out, dtype={"water_level_source": str, "period_source": str, "profile": str})
    d = d[d["status"].isin(["ok", "low contrast"]) & d["hs_m"].notna()].copy()
    d["start"] = d["filename"].astype(str).str.split(".").str[0].astype(float)
    lo = _day_epoch(args.since) if args.since else -np.inf
    hi = _day_epoch(args.until, end=True) if args.until else np.inf
    d = d[(d["start"] >= lo) & (d["start"] < hi)]
    if not len(d):
        sys.exit("no stack with a wave height yet" + (" in --since / --until" if args.since or args.until else ""))
    on_dem = d["profile"].fillna("").str.startswith("dem:")
    for k in ("r2", "sig", "setup"):
        d[f"hs_{k}"] = d[f"h0_{k}_m"] * d["ks_21m"]
        if k != "sig":
            d.loc[~on_dem, f"hs_{k}"] = np.nan             # a planar profile has no level for them
    d["period_kind"] = np.where(d["period_source"].fillna("").str.startswith("camera"), "camera", "buoy")
    d["water_level_source"] = d["water_level_source"].fillna("")
    g = d.groupby(["filename", "status"], as_index=False).agg(
        epoch=("epoch", "first"), start=("start", "first"), water_level_source=("water_level_source", "first"),
        period_kind=("period_kind", "first"), hs_m=("hs_m", "median"), hs_r2=("hs_r2", "median"),
        hs_sig=("hs_sig", "median"), hs_setup=("hs_setup", "median"))
    owg = _read_rows(args.owg_csv)
    g["owg"] = np.nan
    if owg:
        oe = np.array([_num(r.get("epoch")) for r in owg])
        oh = np.array([_num(r.get("hs_m")) if r.get("status") == "ok" else np.nan for r in owg])
        good = np.isfinite(oe) & np.isfinite(oh)
        oe, oh = oe[good], oh[good]
        g["owg"] = [_nearest(oe, oh, t, OWG_MAX_GAP_S) for t in g["start"]]
    g["adcp"] = np.nan
    adcp_path = Path(args.adcp)
    if adcp_path.exists():
        from historical_forcing import read_adcp
        a = read_adcp(adcp_path)
        g["adcp"] = [_nearest(a["ep"].to_numpy(float), a["hs"].to_numpy(float), t, 1800.0) for t in g["epoch"]]
    w = _read_rows(args.waves_csv)
    for c in ("hs_best", "hs_buoy_marconi"):
        g[c] = np.nan
        if w:
            we = np.array([_num(r.get("epoch")) for r in w])
            vals = np.array([_num(r.get(c)) for r in w])
            # the hour the stack starts in, as marconi_waves.py files its frames
            g[c] = [_nearest(we, vals, np.floor(t / 3600.0) * 3600.0, 0.0) for t in g["start"]]
    return g, bool(g["adcp"].notna().any()), out


def compare(args, result=None):
    """--compare: the output file against the references (load_compare), printed; result (a dict),
    when given, receives each table (compare_table)."""
    g, used_adcp, out = load_compare(args)
    span = ""
    if args.since or args.until:
        span = f" from {args.since or 'the start'} to {args.until or 'the end'}"
    print(f"timestack wave height: {len(g)} stack(s) with a value{span} ({(g['status'] == 'ok').sum()} clear, "
          f"{(g['status'] == 'low contrast').sum()} blurred) in {out}")

    def cols_for(ref):
        tag = " [not indep.]"
        return [("from S_IG (camera only)", "hs_sig"), ("hs_m (combined)" + tag, "hs_m"),
                ("from R2" + tag, "hs_r2"), ("from setup" + tag, "hs_setup")]
    refs = []
    if used_adcp:
        refs.append(("adcp", "the ADCP (wh_4061, 21 m; the truth the OWG was trained on)"))
    else:
        print("  no stack falls in the ADCP record (Dec 2024 - Mar 2025): the converted buoy is the reference")
    if g["hs_buoy_marconi"].notna().any():
        refs.append(("hs_buoy_marconi", "the buoy converted to Marconi (hs_buoy_marconi: no camera in it)"))
    for ref, name in refs:
        compare_table(g, ref, name, cols_for(ref), "owg", result)
    if g["hs_best"].notna().any():
        compare_table(g, "hs_best", "hs_best (a consistency check only: it CONTAINS the OWG, and the DEM's "
                      "setup was computed from it)", cols_for("hs_best"), None, result)
    if g["owg"].notna().any():
        compare_table(g, "owg", "the OWG's ok frames of the same collections", cols_for("owg"), None, result)
    print("\n  hs_* is the deep-water H0 shoaled to 21 m (Ks, linear theory): the currency of the ADCP,")
    print("  hs_best and the OWG. [not indep.]: R2 and setup read the DEM's LEVEL, which carries the")
    print("  station's setup C x sqrt(Hs L0) with Hs = hs_best (in storms mostly hs_buoy_marconi) or, on a")
    print("  survey-date DEM, the ADCP's: their agreement with that reference is partly built in. S_IG")
    print("  reads only the DEM's slope: 'from S_IG' is the camera-only candidate (with a buoy period it")
    print("  shares the buoy's Tp). Decision rule: docs/VALIDATION.md 4.11.")
    return 0


# ---------------------------------------------------------------------

def build_parser():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--ras-dir", default=str(HERE / "archive" / "ras_c2"))
    ap.add_argument("--camera", default="c2")
    ap.add_argument("--lines", nargs="+", type=int, default=list(DEFAULT_LINES),
                    help="timestack lines (default 1 2: the c2 lines that cross the swash)")
    ap.add_argument("--pix", default=None, help="default: <cam>_timestack.pix in arguseyes/build, or a saved "
                    "copy matching each stack's width (runup_from_timestack.pix_for_width)")
    ap.add_argument("--io", default=None, help="default calibration/CACO05_<cam>_20240801_IO.yaml")
    ap.add_argument("--eo", default=None, help="default calibration/CACO05_<cam>_20251113_EO-CV.yaml "
                    "(stacks before Nov 2025 need their own: calibration/calibration_history.csv)")
    ap.add_argument("--output", default=None, help="default archive/hs_timestack_<cam>.csv")
    ap.add_argument("--dem", default=str(HERE / "archive" / "dems"),
                    help="folder of dated DEMs (the cron's archive/dems; nearest within 7 days) or one "
                         "DEM .asc for every stack")
    ap.add_argument("--beta", type=float, default=None,
                    help="planar beach of this slope through the Jan 2025 0 m contour where no DEM covers "
                         "the line, and the foreshore slope where the DEM gives none (default: none)")
    ap.add_argument("--max-extrapolate", type=float, default=MAX_EXTRAPOLATE_M,
                    help=f"m the DEM profile is carried on past its edge (default {MAX_EXTRAPOLATE_M:g})")
    ap.add_argument("--max-retreat-speed", type=float, default=None,
                    help="fastest the swash edge may move seaward, m/s ON THE BEACH (default: "
                         "runup_from_timestack.py's 3, which it applies along the line's path on the "
                         "still-water plane)")
    ap.add_argument("--period-source", choices=["auto", "camera", "buoy"], default="auto",
                    help="T0: camera Tp when ok, else the buoy (auto, default); or only one of them")
    ap.add_argument("--period-csv", default=str(HERE / "archive" / "wave_period_c2.csv"))
    ap.add_argument("--waves-csv", default=str(HERE / "archive" / "waves_marconi.csv"))
    ap.add_argument("--buoy-csv", default=str(HERE / "archive" / "waves_44008.csv"))
    ap.add_argument("--owg-csv", default=str(HERE / "archive" / "owg_c2_H.csv"),
                    help="owg_live.py's c2 frames: the window sharpness per stack, and --compare")
    ap.add_argument("--adcp", default=str(HERE / "sig1000_waves_ALL.csv"), help="for --compare")
    ap.add_argument("--min-contrast", type=float, default=MIN_CONTRAST,
                    help=f"edge contrast below which a stack is 'low contrast' (default {MIN_CONTRAST:g}, "
                         "provisional)")
    ap.add_argument("--scatter", default=None,
                    help="error model of the combination, e.g. 'r2=0.26:0.28:0.15,sig=0.26:0,setup=0:0.36' "
                         "(relative:absolute m[:slope relative] per estimate; default the module's SCATTER)")
    ap.add_argument("--wait-days", type=float, default=WAIT_DAYS,
                    help=f"leave stacks younger than this without a period, DEM, water level or (where the "
                         f"GPS is in use) the GPS level for a later run (default {WAIT_DAYS:g})")
    ap.add_argument("--level-csv", default=None,
                    help="water level (epoch, level_navd88, or the formats historical_forcing.py reads) "
                         "instead of marconi_water_level.py")
    from marconi_water_level import GNSSR_SPLINE, GAUGE_CSV
    ap.add_argument("--gnssr-spline", default=GNSSR_SPLINE)
    ap.add_argument("--gauge-csv", default=GAUGE_CSV)
    ap.add_argument("--no-gps", action="store_true", help="water level from the Chatham transfer only")
    ap.add_argument("--since", default=None, help="only stacks from this UTC day (YYYY-MM-DD); with "
                    "--compare, only these are scored (the held-out stacks: --since the day after D)")
    ap.add_argument("--until", default=None, help="only stacks up to this UTC day (YYYY-MM-DD), runs and "
                    "--compare alike")
    ap.add_argument("--reprocess", action="store_true",
                    help="redo the stacks (of --since/--until) and replace their rows")
    ap.add_argument("--compare", action="store_true",
                    help="against the ADCP (Dec 2024 - Mar 2025 stacks) and the converted buoy, with the OWG "
                         "on the same stacks; hs_best and the OWG as consistency checks")
    ap.add_argument("--self-test", action="store_true")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.self_test:
        return self_test()
    if args.compare:
        return compare(args)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())

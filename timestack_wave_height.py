#!/usr/bin/env python3
"""
Wave Height From The Runup Timestacks (Camera Only, Opt-In)
=============================================================
Estimates the significant wave height off Marconi from each 10-minute
c2 runup timestack (archive/ras_c2, 2 Hz): the swash edge that
runup_from_timestack.py tracks is put on the beach, its runup, setup and
swash are measured as Stockdon et al. (2006) defined them, and their
parameterization is inverted for the deep-water wave height H0 in three
independent ways. OPT-IN: it writes only archive/hs_timestack_c2.csv. The
cron, every plot and every published value are untouched; whether it is
ever shown or blended is decided by docs/VALIDATION.md 4.11.

THE RULE
 1. EDGE. runup_from_timestack.swash_edge with that script's own settings
    (foam threshold, 15 px tongues, 2.5 s median), on lines 1 and 2 (the
    c2 lines that cross the swash). One change: its backwash limit, 3 m/s,
    is meant on the beach, but it is applied on the still-water plane, and
    c2 (19.8 m up, the swash 150-300 m away at a grazing angle) stretches
    the line on that plane: 2.1-2.8 plane metres per beach metre on a
    1:7-1:11 beach 230 m away (the real c2 calibration, a synthetic line and
    beach), so the limit is 1.1-1.4 m/s on the beach, slower than real
    backwash. It is therefore scaled by the line's own stretch over the
    swash (plane_stretch): 3 m/s on the beach (--max-retreat-speed). On the
    synthetic scene the plane limit made the setup error +0.09 to +0.11 m
    instead of +0.06 m and cut S_inc by 10%. For the same reason
    runup_from_timestack.py's *_horiz_m are plane metres, 2-3 of them to a
    beach metre on such a line (that script is unchanged).
 2. ELEVATION. Each pixel of the line is a ray from the camera, followed
    down to where it meets the station's 7-day DEM (archive/dems: the DEM
    whose window centre is nearest the stack, within 7 days, picked as
    setup_model.EvolvingSlope picks it) in 1 cm steps, read bilinearly and
    only where all four cells exist. Past the DEM's ends the profile goes on
    as the straight line through its last 10 m, for at most 10 m
    (--max-extrapolate; z_extrapolated_fraction is the share of the edge
    read there); gaps up to 6 m are bridged. Without a DEM, --beta gives a
    plane through the Jan 2025 0 m contour (timestack_wave_period.py's beach
    geometry; profile 'planar'). The edge's elevation in each 2 Hz row is the
    profile at its (fractional) pixel, minus still water
    (marconi_water_level.py, GPS first, interpolated across the burst).
    Plane positions times a slope would not do: on the synthetic scene the
    plane spread the swash 2.1-2.7 times wider (sd) than the beach did.
 3. STATISTICS, Stockdon's definitions. R2: the 98th percentile of the
    uprush maxima above still water (the edge's local maxima, found as
    runup_from_timestack.py finds them); an uprush past the landward end of
    the profile counts as higher than any value on it, and if the 98th
    percentile falls among those, no R2 ('R2 above the profile'). Setup
    <eta>: the mean (none when over 2% of the edge is off the profile).
    Swash: 4 x the standard deviation of the series above 0.05 Hz (S_inc)
    and below (S_IG), Stockdon's split, with runup_from_timestack.band_sig
    (gaps up to 5 s bridged; none with 10% missing), after the straight
    line through the series is removed (band_sig removes only the mean,
    and a drift of 0.5 m across the record adds 10% or more to S_IG; the
    tide, at most ~0.1 m in 10 minutes, goes with still water anyway).
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
    buoy is the check.
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
    not small waves). Eq. 19 is Stockdon's for xi0 = b / sqrt(H0/L0) >= 0.3,
    with R2 = 0.043 sqrt(H0 L0) (Eq. 18) below; on Marconi's 0.07-0.14
    foreshore xi0 < 0.3 needs a deep-water steepness H0/L0 above 0.05-0.22,
    steeper than storm seas usually are, so Eq. 19 is used throughout and
    iribarren is written for the check.
 7. COMBINATION (combine): the inverse-variance mean of log H0. A relative
    error e of a statistic is 2e in H0. The error of each statistic: R2 26%
    (Stockdon's own scatter, 38 cm rms on a mean observed R2 of 144 cm) and
    0.28 m of level (still water now, 0.12 m: dem_from_contours.SIGMA_WL_M;
    the DEM's own level, its median sigma_obs 0.258 m: VALIDATION 4.10);
    S_IG 26% (Stockdon's per-band scatter was not retrievable here, so the R2
    figure is assumed) and its sampling error, sqrt(2/dof)/2 = 9% for 10
    minutes (58 degrees of freedom), no level term; setup 0.21 m (Stockdon's
    setup scatter re-evaluated by Dalinghaus, Coco & Higuera 2023:
    dem_from_contours.SIGMA_SETUP_M) and the 0.28 m of level. The errors are
    taken at the COMBINED H0, iterated, not at each estimate's own (that would
    give a reading that is too high a smaller error and more weight). R2 and
    setup share the level error, so h0_sigma_rel (the sigma of log H0) is a
    floor. A planar profile does not know where still water meets the beach:
    there only S_IG is combined. --scatter replaces the error model once the
    NUC has measured each estimate (VALIDATION 4.11).
 8. CURRENCY. H0 is a DEEP-WATER height, Stockdon's (their nearshore wave
    records taken back to deep water with linear theory, as setup_model.py
    records from how the paper is cited). One shoaling step is applied, the
    other way: hs_m = H0 x Ks(T0, 21 m) (setup_model.shoaling_coefficient,
    0.91-0.97 for 6-14 s) is the height at the ADCP's 21 m, the currency of
    hs_best and of the OWG (both fitted to the ADCP's wh_4061), so hs_m is
    the number to compare with them. Refraction is not undone either way
    (setup_model's choice). Being
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
  A stack younger than --wait-days (3) still lacking a camera period, a DEM
  or a water level is left for a later run, not written.
  EDGE CONTRAST is the Michelson contrast (foam - sand) / (foam + sand) of
  the mean grey 2-8 px either side of the tip, median over the uprush maxima
  (where the tip stands still). Blur spreads the tip; spray haze lifts both
  sides. On the synthetic scene 0.21-0.23 clear; with 50% haze and a 2 / 4 /
  6 / 10 px blur 0.096 / 0.065 / 0.047 / 0.027. 0.05 is PROVISIONAL: no real
  stack was available here, and the NUC sets it from 20-24 Oct 2026 against
  clear days (VALIDATION 4.11).

WHAT IT CANNOT DO
  * CIRCULARITY. The DEM is built from the timex waterlines with the
    station's setup C sqrt(Hs L0) added, Hs = hs_best, which contains the
    OWG. R2 and setup read the DEM's LEVEL, so the setup model's error on the
    DEM's frames moves both, and part of their agreement with hs_best is
    built in. S_IG reads only the DEM's SLOPE (elevation differences a few
    metres apart: a setup error common to a frame cancels). Fair references
    are the ADCP (Dec 2024 - Mar 2025) and the converted buoy only.
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
    -8 to +32%; combined H0 +9 to +20%. Real foam is not synthetic foam:
    only the NUC comparison measures the real errors.
  * Ten minutes is short for infragravity statistics: sampling alone is 9%
    on S_IG, 19% on H0 from it.
  * The period (5. above).
  * The DEM is a week's waterlines, i.e. the mean swash of the week's tides:
    at high water and in big waves the swash runs above it (10 m of
    extension, then no R2), and a berm built or cut within the week is not in
    it.
  * Stockdon's relations are fits over ten experiments on other beaches. A
    site coefficient could be fitted to this beach once truth exists, which
    is what VALIDATION 4.11 is for.

OUTPUT archive/hs_timestack_c2.csv, one row per stack and line, APPENDED:
rows already there are never rewritten except by --reprocess, and a file
with other columns is refused.
  epoch (burst middle), time_utc, filename, line, status, note, brightness,
  edge_contrast, window_sharpness, valid_fraction, n_peaks, sea_at, profile,
  z_extrapolated_fraction, water_level_m, water_level_source, beta_f,
  beta_source, iribarren, period_s, period_source, setup_m, r2_m, s_inc_m,
  s_ig_m, s_m (m, above still water), h0_r2_m, h0_sig_m, h0_setup_m (each
  estimate, deep water), h0_m, h0_sigma_rel, combined_from (the
  combination), ks_21m, hs_m (H0 x Ks at 21 m: compare this one)

References: Stockdon, Holman, Howd & Sallenger 2006, Coastal Engineering
53, 573-588 (its full text was not retrievable here; equations checked
against py-wave-runup's Stockdon2006); USGS Data Series 602 (its runup
statistics); Holman & Sallenger 1985, J. Geophys. Res. 90(C1), 945-953
(setup and swash from video); Guza & Thornton 1982, J. Geophys. Res.
87(C1), 483-491; Ruessink, Kleinhans & van den Beukel 1998, J. Geophys.
Res. 103(C2), 3111-3118; Dalinghaus, Coco & Higuera 2023, Nat. Hazards
Earth Syst. Sci. 23, 2157-2169.

Usage:
    python3 timestack_wave_height.py                  new stacks in archive/ras_c2
    python3 timestack_wave_height.py --reprocess      every archived stack again
    python3 timestack_wave_height.py --compare        against ADCP / converted buoy / hs_best / OWG
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
from runup_from_timestack import SAMPLE_HZ, IG_CUTOFF_HZ, band_sig, peak_mask  # noqa: E402
from setup_model import (deep_water_wavelength, stockdon_setup, shoaling_coefficient,  # noqa: E402
                         slope_grid, ForeshoreSlope, EvolvingSlope, STOCKDON_SETUP_K,
                         STOCKDON_SWASH_IG, DESHOAL_DEPTH_BY_CURRENCY, BETA_RANGE, SLOPE_MAX_DAYS)
from dem_from_contours import SIGMA_WL_M, SIGMA_SETUP_M                     # noqa: E402
from station_status import LENS_BLURRED                                     # noqa: E402
from timestack_wave_period import SHORE_ORIGIN, SEAWARD_DEG                 # noqa: E402

FIELDS = ["epoch", "time_utc", "filename", "line", "status", "note",
          "brightness", "edge_contrast", "window_sharpness", "valid_fraction", "n_peaks", "sea_at",
          "profile", "z_extrapolated_fraction", "water_level_m", "water_level_source",
          "beta_f", "beta_source", "iribarren", "period_s", "period_source",
          "setup_m", "r2_m", "s_inc_m", "s_ig_m", "s_m",
          "h0_r2_m", "h0_sig_m", "h0_setup_m", "h0_m", "h0_sigma_rel", "combined_from",
          "ks_21m", "hs_m"]
STATUSES = ("ok", "low contrast", "no swash edge", "no slope", "no period", "dark")

DEFAULT_LINES = (1, 2)            # the c2 lines that cross the swash (waterline_timex_cron.sh RUNUP_LINES)
DARK_MEAN = 25.0                  # grey level: timestack_wave_period.period_from_stack's "dark"
MIN_CONTRAST = 0.05               # Michelson contrast across the swash tip: provisional, see EDGE CONTRAST
CONTRAST_NEAR, CONTRAST_FAR = 2, 8    # pixels either side of the tip compared by edge_contrast
OWG_MAX_GAP_S = 600.0             # the bright frame of the same 10-min collection
PERIOD_MAX_GAP_S = 5400.0         # a buoy hour within 1.5 h of the burst's middle
WAIT_DAYS = 3.0                   # inputs that arrive late (period, DEM, water level) are waited for
DESHOAL_DEPTH_M = DESHOAL_DEPTH_BY_CURRENCY["marconi"]     # 21 m: the ADCP / hs_best / OWG currency

# The beach profile along a line
DZ_STEP = 0.01                    # m: the ray is marched down the DEM in these steps
EXT_FIT_M = 10.0                  # m of DEM profile at each end the extension is fitted to
MAX_EXTRAPOLATE_M = 10.0          # m beyond the DEM's edge (along the line) the profile is carried on
GAP_FILL_M = 6.0                  # m: DEM gaps inside the profile up to this long are bridged
MISSING_ELEV_MAX = 0.02           # share of edge samples off the profile above which no setup is given
MIN_SLOPE_CELLS = 3               # distinct DEM cells under the line inside the swash band

# Error model of the combination (log H0, inverse variance): per estimate (relative, absolute m)
# scatter of the measured quantity about Stockdon's relation, see THE COMBINATION.
DEM_SIGMA_M = 0.258               # median per-cell sigma_obs of the station's DEM, VALIDATION 4.10 [REAL]
LEVEL_SIGMA_M = float(np.hypot(SIGMA_WL_M, DEM_SIGMA_M))   # still water now and the DEM's level: 0.28 m
R2_REL = 0.38 / 1.44              # Stockdon et al. (2006): R2 rms 38 cm on a mean observed 144 cm
SCATTER = {"r2": (R2_REL, LEVEL_SIGMA_M),
           "sig": (R2_REL, 0.0),
           "setup": (0.0, float(np.hypot(SIGMA_SETUP_M, LEVEL_SIGMA_M)))}


# ---------------------------------------------------------------------
# Stockdon et al. (2006), forward and inverted
# ---------------------------------------------------------------------

def stockdon_r2(h0, t0, beta):
    """R2 = 1.1 [0.35 beta sqrt(H0 L0) + sqrt(H0 L0 (0.563 beta^2 + 0.004)) / 2]: Stockdon et al.
    (2006) Eq. 19 as published, the relation runup_inversion.h0_from_runup inverts. Its swash
    term rounds Eqs. 11-12 (0.75^2 = 0.5625, 0.06^2 = 0.0036), so 1.1 (<eta> + S/2) from
    setup_model.stockdon_setup / stockdon_swash is 0-3% lower."""
    h0, t0, beta = (np.asarray(a, float) for a in (h0, t0, beta))
    x = h0 * deep_water_wavelength(t0)
    return 1.1 * (STOCKDON_SETUP_K * beta * np.sqrt(x) + np.sqrt(x * (0.563 * beta ** 2 + 0.004)) / 2.0)


def stockdon_sig(h0, t0):
    """Infragravity swash S_IG = 0.06 sqrt(H0 L0), Stockdon et al. (2006) Eq. 12 (no slope in it)."""
    return STOCKDON_SWASH_IG * np.sqrt(h0 * deep_water_wavelength(t0))


def h0_from_r2(r2, beta, t0):
    """H0 from R2 above still water: runup_inversion.h0_from_runup (closed form; NaN for R2 <= 0)."""
    from runup_inversion import h0_from_runup
    return float(h0_from_runup(np.asarray([r2], float), np.asarray([beta], float),
                               np.asarray([t0], float))[0])


def h0_from_sig(s_ig, t0):
    """H0 = (S_IG / 0.06)^2 / L0: Eq. 12 inverted. NaN for S_IG <= 0 or no period."""
    s_ig, t0 = float(s_ig), float(t0)
    if not (np.isfinite(s_ig) and s_ig > 0 and np.isfinite(t0) and t0 > 0):
        return np.nan
    return (s_ig / STOCKDON_SWASH_IG) ** 2 / float(deep_water_wavelength(t0))


def h0_from_setup(setup, beta, t0):
    """H0 = (setup / (0.35 beta))^2 / L0: Eq. 10 inverted. NaN for setup <= 0 (below still water:
    an error in the level, not small waves), no slope or no period."""
    setup, beta, t0 = float(setup), float(beta), float(t0)
    if not (np.isfinite(setup) and setup > 0 and np.isfinite(beta) and beta > 0
            and np.isfinite(t0) and t0 > 0):
        return np.nan
    return (setup / (STOCKDON_SETUP_K * beta)) ** 2 / float(deep_water_wavelength(t0))


def iribarren(beta, h0, t0):
    """xi0 = beta / sqrt(H0 / L0); Stockdon's Eq. 19 is for xi0 >= 0.3 (Eq. 18 below)."""
    if not (np.isfinite(beta) and np.isfinite(h0) and h0 > 0 and np.isfinite(t0)):
        return np.nan
    return float(beta / np.sqrt(h0 / deep_water_wavelength(t0)))


def predicted(kind, h0, t0, beta):
    """The swash statistic an estimate was inverted from, as Stockdon predicts it for H0."""
    if kind == "r2":
        return float(stockdon_r2(h0, t0, beta))
    if kind == "sig":
        return float(stockdon_sig(h0, t0))
    if kind == "setup":
        return float(stockdon_setup(h0, t0, beta))
    raise ValueError(kind)


def combine(estimates, t0, beta, scatter=None, sig_sampling=0.0, passes=20):
    """
    One H0 from the separate estimates {'r2' | 'sig' | 'setup': H0 m}: the inverse-variance mean
    of log H0. Each estimate inverts a statistic q that Stockdon puts proportional to sqrt(H0), so a
    relative error e in q is 2e in H0. The error of q is sqrt((rel q)^2 + abs^2) (scatter, per
    kind), plus for 'sig' its sampling error (sig_sampling, relative). It is evaluated at the
    COMBINED H0, not at each estimate's own (that would give an estimate that reads high a smaller
    relative error and so more weight), iterated to a fixed point from the geometric mean.
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
            rel, ab = scatter[k]
            q = predicted(k, h, t0, beta)
            if k == "sig":
                rel = float(np.hypot(rel, sig_sampling))
            s_q = np.hypot(rel * q, ab)
            w[i] = (q / (2.0 * s_q)) ** 2 if q > 0 and s_q > 0 else 0.0
        if w.sum() <= 0:
            return np.nan, np.nan, ""
        h_new = float(np.exp((w * logs).sum() / w.sum()))
        done = abs(h_new - h) <= 1e-9 * h
        h = h_new
        if done:
            break
    return h, float(1.0 / np.sqrt(w.sum())), "+".join(names)


def parse_scatter(text):
    """'r2=0.26:0.28,sig=0.26:0,setup=0:0.36' (relative:absolute m) -> a SCATTER dict."""
    out = dict(SCATTER)
    if not text:
        return out
    for part in str(text).split(","):
        k, _, v = part.partition("=")
        k = k.strip()
        if k not in SCATTER:
            raise ValueError(f"--scatter: unknown estimate {k!r} (r2, sig, setup)")
        rel, _, ab = v.partition(":")
        out[k] = (float(rel), float(ab or 0.0))
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
        """Bilinear between cell centres; NaN where any of the four cells is missing."""
        r, c = self._rc(E, N)
        with np.errstate(invalid="ignore"):
            r0, c0 = np.floor(r), np.floor(c)
        fr, fc = r - r0, c - c0
        z00 = self._cell(self.z, r0, c0)
        z01 = self._cell(self.z, r0, c0 + 1)
        z10 = self._cell(self.z, r0 + 1, c0)
        z11 = self._cell(self.z, r0 + 1, c0 + 1)
        return z00 * (1 - fr) * (1 - fc) + z01 * (1 - fr) * fc + z10 * fr * (1 - fc) + z11 * fr * fc

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
        """setup_model.ForeshoreSlope of this DEM (the whole grid), for the band fallback."""
        if self._foreshore is None:
            self._foreshore = ForeshoreSlope(self.z, self.cell, self.name)
        return self._foreshore


class DemSource:
    """The DEM for a stack: a dated series (dem_from_contours.py --series-dir, the cron's archive/dems)
    picked as setup_model.EvolvingSlope picks it -- the DEM whose window centre is nearest, within
    max_days -- or one DEM file for every stack."""

    def __init__(self, source, max_days=SLOPE_MAX_DAYS):
        self.ev = EvolvingSlope(source, max_days) if source else None
        self._cache = {}

    def path_for(self, epoch):
        ev = self.ev
        if ev is None:
            return None
        if ev.single is not None:
            return ev.single
        if not ev.series:
            return None
        i = int(np.argmin(np.abs(ev.centres - float(epoch))))
        if abs(ev.centres[i] - float(epoch)) > ev.max_days * 86400.0:
            return None
        return ev.series[i][1]

    def for_epoch(self, epoch):
        p = self.path_for(epoch)
        if p is None:
            return None
        if p not in self._cache:
            try:
                self._cache[p] = Dem(p)
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
    EXT_FIT_M of DEM profile at that end, for at most max_extrapolate metres, rising landward or not
    at all; gaps inside up to GAP_FILL_M are bridged. kind 'dem'.
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
                if end == "land":
                    sel = dcov >= dcov.max() - EXT_FIT_M
                    targets = np.arange(ci[-1] + 1, n)
                    d_end = dcov.max()
                else:
                    sel = dcov <= dcov.min() + EXT_FIT_M
                    targets = np.arange(0, ci[0])
                    d_end = dcov.min()
                if sel.sum() < 10 or np.ptp(dcov[sel]) < 3.0 or not len(targets):
                    continue
                b, a = np.polyfit(dcov[sel], zcov[sel], 1)
                if not b > 0:
                    continue
                with np.errstate(invalid="ignore", divide="ignore"):
                    zt = (a + b * d0[targets]) / (1.0 - b * d1[targets])
                    dt = d0[targets] + d1[targets] * zt
                    ok = np.isfinite(zt) & ((dt - d_end) * (1 if end == "land" else -1) > 0) & \
                        (np.abs(dt - d_end) <= max_extrapolate)
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


def elevation_at(edge, z_line):
    """
    Elevation of the edge (fractional column index, sea end first) on the profile z_line, linear
    between pixels. -> (z, above, below): NaN where the edge is NaN or off the profile; above /
    below flag edge samples landward / seaward of the profile's ends (an uprush above the profile
    is still known to be above every elevation the profile has).
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
    z[idx] = val
    miss = idx[~np.isfinite(val)]
    above[miss] = edge[miss] > fin[-1]
    below[miss] = edge[miss] < fin[0]
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
    """Degrees of freedom of the infragravity variance of an n-sample record: 2 per frequency bin
    in (0, cutoff) -- 58 for 10 min at 2 Hz."""
    f = np.fft.rfftfreq(int(n), 1.0 / fs)
    return 2 * int(((f > 0) & (f < cutoff)).sum())


def swash_stats(z_abs, swl, edge, above, below, valid, min_peaks, min_valid):
    """
    z_abs: elevation of the swash edge (m NAVD88) per 2 Hz sample; swl: still-water level per
    sample (NaN if unknown); edge: the edge (fractional column index, for the uprush peaks).
    -> dict setup_m (mean of z - swl), r2_m (98th percentile of the uprush maxima of z - swl, the
    maxima found on the edge series as runup_from_timestack.py finds them: peak_mask), s_inc_m,
    s_ig_m, s_m (4 x std of the detrended z - swl in > 0.05 Hz, < 0.05 Hz, all: band_sig), n_peaks,
    eta and sigma (mean and std of z_abs, for the foreshore band), missing (share of edge samples
    with no elevation), notes.
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
        pz = pz[~np.isnan(pz)]
        if n_peaks >= min_peaks and valid >= min_valid and len(pz) >= min_peaks:
            with np.errstate(invalid="ignore"):
                r2 = float(np.percentile(pz, 98))
            if not np.isfinite(r2):
                r2 = np.nan
                notes.append("R2 above the profile")
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
    if not np.isfinite(s_ig):
        notes.append("edge gaps: no swash spectrum")
    fin = z_abs[np.isfinite(z_abs)]
    return {"setup_m": setup, "r2_m": r2, "s_inc_m": s_inc, "s_ig_m": s_ig, "s_m": s_all,
            "n_peaks": n_peaks, "eta": float(fin.mean()) if len(fin) else np.nan,
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
    """How many metres along the line on the still-water plane (runup_from_timestack's dist) one
    metre on the beach is, over the pixels the swash edge covers (its 5-95th percentile): c2 sees
    the swash at a grazing angle from 19.8 m up, so a point above still water lands further along
    its ray on that plane. 1.0 where it cannot be measured."""
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


def estimate_stack(ras_path, line, ctx, now=None):
    """The row for one stack and line, None if the stack cannot be read, or 'wait' when an input that
    arrives late (camera period, DEM, water level) is missing and the stack is younger than
    ctx.wait_days."""
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
    # the retreat limit, on the beach: the runup script applies it on the still-water plane
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
    if young and (pending or status in ("no period", "no slope")
                  or (status == "ok" and not np.isfinite(wl_mid))):
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
        z_abs, above, below = elevation_at(edge, prof["z"])
        st = swash_stats(z_abs, swl, edge, above, below, valid, rargs.min_peaks, rargs.min_valid)
        notes += st["notes"]
        with np.errstate(invalid="ignore"):
            on = np.isfinite(z_abs)
            xf = float(prof["extrap"][np.clip(np.round(edge[on]).astype(int), 0, len(prof["z"]) - 1)].mean()) \
                if on.any() else np.nan
        beta, bsrc = foreshore_beta(prof, dem if prof["kind"] == "dem" else None, st["eta"], st["sigma"],
                                    ctx.beta)
        est = {"sig": h0_from_sig(st["s_ig_m"], t0)}
        absolute = prof["kind"] == "dem"            # R2 and setup need the beach's true level
        est["r2"] = h0_from_r2(st["r2_m"], beta, t0) if np.isfinite(beta) else np.nan
        est["setup"] = h0_from_setup(st["setup_m"], beta, t0) if np.isfinite(beta) else np.nan
        if not np.isfinite(beta):
            notes.append("no foreshore slope: no R2 or setup estimate")
        use = dict(est) if absolute else {"sig": est["sig"]}
        if not absolute:
            notes.append("planar profile: R2 and setup not combined")
        samp = np.sqrt(2.0 / max(ig_dof(n), 1)) / 2.0
        h0, sig_ln, used = combine(use, t0, beta if np.isfinite(beta) else 0.1, ctx.scatter, samp)
        ks = float(shoaling_coefficient(t0, DESHOAL_DEPTH_M))
        row.update({"z_extrapolated_fraction": _r(xf, 3), "beta_f": _r(beta, 4), "beta_source": bsrc,
                    "iribarren": _r(iribarren(beta, h0, t0), 2) if np.isfinite(beta) else "",
                    "setup_m": _r(st["setup_m"], 3), "r2_m": _r(st["r2_m"], 3),
                    "s_inc_m": _r(st["s_inc_m"], 3), "s_ig_m": _r(st["s_ig_m"], 3), "s_m": _r(st["s_m"], 3),
                    "h0_r2_m": _r(est["r2"], 3), "h0_sig_m": _r(est["sig"], 3),
                    "h0_setup_m": _r(est["setup"], 3), "h0_m": _r(h0, 3), "h0_sigma_rel": _r(sig_ln, 3),
                    "combined_from": used, "ks_21m": _r(ks, 4), "hs_m": _r(h0 * ks, 3)})
        xi = iribarren(beta, h0, t0)
        if np.isfinite(xi) and xi < 0.3:
            notes.append(f"xi0 {xi:.2f} < 0.3: dissipative, where Stockdon use Eq. 18 for R2")
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


def write_rows(path, rows, reprocess=False):
    """
    Appends rows to the CSV (header first if new). Rows whose (filename, line) are already in it
    are not appended again. With reprocess, the file is rewritten once: its rows that this run
    redid are replaced, the others kept, all in (epoch, line) order. A file with another header is
    refused (SystemExit), so a column change never shifts old rows: rerun with --reprocess.
    -> the number of rows written.
    """
    path = Path(path)
    if not rows:
        return 0
    rows = sorted(rows, key=lambda r: (int(r["epoch"]), int(r["line"])))
    exists = path.exists() and path.stat().st_size > 0
    if exists:
        with open(path, newline="") as f:
            header = next(csv.reader(f), [])
        if header != FIELDS:
            if not reprocess:
                raise SystemExit(f"{path}: its columns are not this version's; rerun with --reprocess")
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
    return ctx


def _day_epoch(s, end=False):
    t = datetime.strptime(s, "%Y-%m-%d").replace(tzinfo=timezone.utc).timestamp()
    return t + 86400.0 if end else t


def run(args, wl=None, now=None):
    out = Path(args.output or HERE / "archive" / f"hs_timestack_{args.camera}.csv")
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
    rows, waiting = [], 0
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
            rows.append(r)
    n = write_rows(out, rows, reprocess=args.reprocess)
    counts = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    ok = [r for r in rows if r["status"] == "ok" and r["hs_m"] != ""]
    print(f"timestack wave height {args.camera}: {n} new row(s) -> {out}"
          + (" (" + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())) + ")" if counts else "")
          + (f"; {waiting} line-stack(s) wait for a period, DEM or water level" if waiting else "")
          + (f"; latest Hs {ok[-1]['hs_m']} m ({ok[-1]['time_utc'][:16]}Z)" if ok else ""))
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


HS_BANDS = ((0.0, 1.0), (1.0, 2.0), (2.0, 99.0))


def compare_table(df, ref_col, ref_name, cols):
    """Prints n / bias / RMSE of each estimate column against ref_col: all, by Hs band of the
    reference, and clear (status ok) against blurred (status low contrast)."""
    print(f"\n  against {ref_name}:")
    print(f"    {'estimate':<26} {'subset':<16} {'n':>5} {'bias m':>8} {'RMSE m':>8} {'r':>6}")
    subsets = [("all clear", df["status"] == "ok")]
    for a, b in HS_BANDS:
        subsets.append((f"clear, {a:g}-{b:g} m" if b < 99 else f"clear, >= {a:g} m",
                        (df["status"] == "ok") & (df[ref_col] >= a) & (df[ref_col] < b)))
    subsets.append(("blurred", df["status"] == "low contrast"))
    for label, col in cols:
        for name, m in subsets:
            n, bias, rmse, r = scores(df.loc[m, col], df.loc[m, ref_col])
            if n:
                print(f"    {label:<26} {name:<16} {n:>5} {bias:>+8.3f} {rmse:>8.3f} {r:>6.2f}")


def compare(args):
    import pandas as pd
    out = Path(args.output or HERE / "archive" / f"hs_timestack_{args.camera}.csv")
    if not out.exists():
        sys.exit(f"no {out} yet: run without --compare first")
    d = pd.read_csv(out)
    d = d[d["status"].isin(["ok", "low contrast"]) & d["hs_m"].notna()].copy()
    if not len(d):
        sys.exit("no stack with a wave height yet")
    for k in ("r2", "sig", "setup"):
        d[f"hs_{k}"] = d[f"h0_{k}_m"] * d["ks_21m"]
    # one value per stack: the median of its lines (of the same status)
    g = d.groupby(["filename", "status"], as_index=False).agg(
        epoch=("epoch", "first"), hs_m=("hs_m", "median"), hs_r2=("hs_r2", "median"),
        hs_sig=("hs_sig", "median"), hs_setup=("hs_setup", "median"))
    g["start"] = g["filename"].str.split(".").str[0].astype(float)
    print(f"timestack wave height: {len(g)} stack(s) with a value ({(g['status'] == 'ok').sum()} clear, "
          f"{(g['status'] == 'low contrast').sum()} blurred) in {out}")
    cols = [("timestack hs_m (combined)", "hs_m"), ("  from R2", "hs_r2"), ("  from S_IG", "hs_sig"),
            ("  from setup", "hs_setup")]

    owg = _read_rows(args.owg_csv)
    if owg:
        oe = np.array([_num(r.get("epoch")) for r in owg])
        oh = np.array([_num(r.get("hs_m")) if r.get("status") == "ok" else np.nan for r in owg])
        good = np.isfinite(oe) & np.isfinite(oh)
        oe, oh = oe[good], oh[good]
        g["owg"] = [_nearest(oe, oh, t, OWG_MAX_GAP_S) for t in g["start"]]
    else:
        g["owg"] = np.nan

    adcp_path = Path(args.adcp)
    used_adcp = False
    if adcp_path.exists():
        from historical_forcing import read_adcp
        a = read_adcp(adcp_path)
        g["adcp"] = [_nearest(a["ep"].to_numpy(float), a["hs"].to_numpy(float), t, 1800.0) for t in g["epoch"]]
        if g["adcp"].notna().any():
            used_adcp = True
            compare_table(g, "adcp", "the ADCP (wh_4061, 21 m; the truth the OWG was trained on)",
                          cols + [("OWG (same stacks)", "owg")])
    if not used_adcp:
        print("  no stack falls in the ADCP record (Dec 2024 - Mar 2025): compared with the station's "
              "own estimates instead, which are NOT independent truth")
        w = _read_rows(args.waves_csv)
        if w:
            we = np.array([_num(r.get("epoch")) for r in w])
            for c in ("hs_best", "hs_buoy_marconi"):
                vals = np.array([_num(r.get(c)) for r in w])
                # the hour the stack starts in, as marconi_waves.py files its frames
                g[c] = [_nearest(we, vals, np.floor(t / 3600.0) * 3600.0, 0.0) for t in g["start"]]
            compare_table(g, "hs_buoy_marconi", "the buoy converted to Marconi (hs_buoy_marconi: no camera in it)",
                          cols + [("OWG (same stacks)", "owg")])
            compare_table(g, "hs_best", "hs_best (the hourly blend: it CONTAINS the OWG, and the DEM's setup "
                          "was computed from it)", cols)
        if g["owg"].notna().any():
            compare_table(g, "owg", "the OWG's ok frames of the same collections", cols)
    # does it complement the OWG? error correlation and the 50:50 blend against the reference
    ref = "adcp" if used_adcp else ("hs_buoy_marconi" if "hs_buoy_marconi" in g else None)
    if ref is not None:
        m = (g["status"] == "ok") & g[ref].notna() & g["owg"].notna() & g["hs_m"].notna()
        if m.sum() >= 5:
            e1, e2 = g.loc[m, "hs_m"] - g.loc[m, ref], g.loc[m, "owg"] - g.loc[m, ref]
            blend = 0.5 * (g.loc[m, "hs_m"] + g.loc[m, "owg"])
            n, _, rb, _ = scores(blend, g.loc[m, ref])
            print(f"\n  with the OWG on {n} clear stack(s), against {ref}: error correlation "
                  f"{np.corrcoef(e1, e2)[0, 1]:+.2f}; RMSE timestack {np.sqrt(np.mean(e1 ** 2)):.3f}, "
                  f"OWG {np.sqrt(np.mean(e2 ** 2)):.3f}, 50:50 blend {rb:.3f} m")
    print("\n  hs_m is the deep-water H0 shoaled to 21 m (Ks, linear theory): the currency of the ADCP,")
    print("  hs_best and the OWG. R2 and setup read the DEM's level, which carries the station's setup")
    print("  (C x sqrt(hs_best L0)); S_IG uses only the DEM's slope. Decision rule: docs/VALIDATION.md 4.11.")
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
                         "runup_from_timestack.py's 3, which it applies on the still-water plane)")
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
                    help="error model of the combination, e.g. 'r2=0.26:0.28,sig=0.26:0,setup=0:0.36' "
                         "(relative:absolute m per estimate; default the module's SCATTER)")
    ap.add_argument("--wait-days", type=float, default=WAIT_DAYS,
                    help=f"leave stacks younger than this without a period, DEM or water level for a "
                         f"later run (default {WAIT_DAYS:g})")
    ap.add_argument("--level-csv", default=None,
                    help="water level (epoch, level_navd88, or the formats historical_forcing.py reads) "
                         "instead of marconi_water_level.py")
    from marconi_water_level import GNSSR_SPLINE, GAUGE_CSV
    ap.add_argument("--gnssr-spline", default=GNSSR_SPLINE)
    ap.add_argument("--gauge-csv", default=GAUGE_CSV)
    ap.add_argument("--no-gps", action="store_true", help="water level from the Chatham transfer only")
    ap.add_argument("--since", default=None, help="only stacks from this UTC day (YYYY-MM-DD)")
    ap.add_argument("--until", default=None, help="only stacks up to this UTC day (YYYY-MM-DD)")
    ap.add_argument("--reprocess", action="store_true",
                    help="redo the stacks (of --since/--until) and replace their rows")
    ap.add_argument("--compare", action="store_true",
                    help="against the ADCP (Jan-Mar 2025 stacks), else hs_best, the converted buoy and the OWG")
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

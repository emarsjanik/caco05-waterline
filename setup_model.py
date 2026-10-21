#!/usr/bin/env python3
"""
Wave Setup On The Waterline Elevation: One Formula, Two Models
=================================================================
A timex waterline is the time-averaged edge of the water, which sits
ABOVE still water by the wave setup (and part of the swash). Each line
is therefore given the still-water level plus a setup computed from the
offshore waves of its frame. This module holds that formula, once.

THE STATION'S MODEL ('c', the default everywhere):

    setup = C * sqrt(Hs * L0),     L0 = g * Tp^2 / (2 * pi)

with C = 0.037 (waterline_timex_cron.sh SETUP_COEF; survey_products.py
and historical_forcing.py carry the same value, and a test keeps the
three equal). It is the setup of Stockdon et al. (2006), Eq. 10,
<eta> = 0.35 * beta_f * sqrt(H0 * L0), with the foreshore slope folded
into ONE site coefficient (C = 0.35 * beta_f: beta_f ~0.106, about 1:9.5,
for C = 0.037) and Hs in the currency C was fitted in (hs_best of
archive/waves_marconi.csv: ADCP Hs in 21 m of water off Marconi; Tp the
NDBC 44008 peak period) instead of Stockdon's deep-water H0. How C was
fitted and every refit since: docs/VALIDATION.md section 4.1.

ONE COPY. Until Oct 2026 this formula was written out wherever it was
needed: extract_elevation_contours.py (the rows' setup), survey_products.py
(setup_m, which apply_setup uses; implied_c; compared_frames_scale),
historical_forcing.py (setup_m; the GNSS-R setup share), dem_from_contours.py
(the C the contours carry; --fit-setup), compare_dem_survey.py,
validate_waves.py, gnssr_qc.py (its Hs * L0) and runup_inversion.py (its
L0), and Stockdon's setup and swash in foreshore_slope.py. They now call the
functions here. The arithmetic is unchanged to the last bit: the same
operations in the same order -- (Hs * g * Tp^2) / (2 pi), its square root,
times C; and for Stockdon L0 = (g * T^2) / (2 pi), sqrt(H0 * L0), times
0.35 * beta. A reordering would move the last bit of some values (floating
point is not associative), and the setup is rounded to 0.1 mm into the
contour files, so test_setup_model.py checks every function and call site
against the expressions as they were written, over a grid of Hs and Tp, bit
for bit. Left as they are: waterline_consistency.py's self-test, which
plants synthetic lines with 0.037 * sqrt(hs * 9.81 * tp * tp / (2 pi)) (test
data; its product is grouped differently, so this module's would not be the
same to the bit), gnssr_qc.py's own setup and runup allowance (a fixed slope
of 0.10, for the GNSS-R quality control only), and the tests. The modules
the GCP-free procedure pins by sha256 (gcpfree_eo.py MODULES) hold no copy:
detect_original_view.py calls the extractor.

THE SLOPE-AWARE MODEL ('stockdon', opt-in: extract_elevation_contours.py
--setup-model stockdon, waterline_timex_cron.sh SETUP_MODEL="stockdon"):

    setup = k * beta_f * sqrt(H0 * L0),     k = 0.35

Stockdon, Holman, Howd & Sallenger (2006), Empirical parameterization of
setup, swash, and runup, Coastal Engineering 53, 573-588,
doi:10.1016/j.coastaleng.2005.12.005 (Eq. 10; the swash, Eqs. 11-12,
sets where the slope is measured). What changes from 'c':

  * beta_f is MEASURED for each frame on the evolving DEM (the station's
    dated 7-day DEMs, dem_from_contours.py --series-dir; or one DEM given
    by the caller) instead of being fixed at ~0.106 inside C. Stockdon
    define the foreshore slope as the mean slope over the region of
    significant swash, <eta> +/- 2 sigma about the setup, i.e. still water
    + setup +/- S/2 with S the significant swash height (4 sigma). That is
    circular (the band needs the setup, the setup needs the slope), so it
    is solved as foreshore_slope.py --iterate does: start from the
    fallback slope, measure the slope of the DEM cells in the band, recentre
    the band, three passes. The slope of a cell is the gradient of a plane
    fitted to its 5 x 5 neighbourhood (foreshore_slope.surface_slope; here
    vectorised, slope_grid), the band's slope the MEDIAN over its cells. One
    difference: the plane is fitted to every DEM cell of the window, not
    only those inside the band (foreshore_slope.py masks the DEM to the band
    first, which on a band a few cells wide leaves the windows at its edges
    half a plane).
  * WHY: one C cannot follow a beach face that is not one plane and that
    changes with the storms. On the real 29 Sep - 5 Oct 2026 DEM the
    foreshore slope over -0.5 to +0.5 m is 0.089-0.125 at four profiles
    (VALIDATION 2.3) [REAL]; the fixture DEMs shaped by the 2025 lidar give
    0.125-0.135 (Jan 2025) and 0.067-0.075 (Mar 2025); a within-cell C of
    0.017 below -0.2 m, 0.045 at -0.2 to 0.6 m and 0.035 above 0.6 m
    [SCRATCH] (VALIDATION 4.1, Known limits). Stockdon's form follows the
    slope. The literature recommends a local, iteratively updated
    foreshore slope from the evolving DEM rather than an alongshore mean
    (Stockdon et al. 2006 warn that an alongshore-averaged slope adds
    error on three-dimensional, intermediate or reflective foreshores --
    foreshore_slope.py quotes their figure, 51% of the fractional slope
    variability in the runup; the review behind this repository,
    best-practice step 5). Whether it is better HERE is for
    evaluate_setup_model.py to say, on the surveys, with leave-one-survey-out
    cross-validation; until it does, the station stays on 'c'.
  * The slope is one value per frame, over the DEM cells of the band (both
    cameras' cells: the DEM is one grid), not per pixel column: the
    contours are tagged in pixel space, before georectification, so a
    point's ground position (and the slope under it) is not known yet. It
    follows the cross-shore shape (the elevation band moves with the tide),
    not the alongshore variation within the view.
  * The slope comes from DEMs built from the waterlines themselves (with the
    setup of their day). A slope is the difference of elevations a few
    metres apart, so a setup error common to a frame cancels in it; what is
    left (a setup that grows with elevation) is a second-order effect.

THE DESHOALING QUESTION (H0). Stockdon's H0 is the DEEP-WATER significant
wave height; their wave records were taken back to deep water with linear
theory (as the paper is usually cited: its full text was not retrieved for
this repository's literature review). The station's Hs is NOT a deep-water
height: hs_best is ADCP Hs at 21 m depth (the camera wave models were
trained on the ADCP's wh_4061 and the buoys converted to it;
historical_forcing.py, buoy_transfer.py), so for Stockdon it should be
deshoaled: H0 = Hs / Ks, Ks = sqrt(Cg0 / Cg) the linear-theory shoaling
coefficient at 21 m (shoaling_coefficient). At 21 m Ks is 0.97 at Tp 6 s,
0.92 at 8 s, 0.913 near 9 s (its minimum), 0.94 at 12 s and 0.97 at 14 s,
so over Tp 5-14 s H0 is 0.7-9.5% above Hs and sqrt(H0 L0) 0.3-4.7% above
sqrt(Hs L0) (swell longer than ~15 s has Ks > 1 there: H0 below Hs); the
tide (+/-1.2 m of depth) moves Ks by under 0.006 for Tp 6-12 s, 0.010 at
16 s (test_setup_model.py checks these). Refraction is NOT undone: that needs the wave direction and the
offshore bathymetry, and whether an 'unrefracted' deep-water height is
what the parameterisation wants is not settled. The raw NDBC 44008 currency
(fetch_buoy_waves.py wvht_m) is taken as deep water (no deshoaling): it is
not the currency C was fitted in anyway. DESHOAL_DEPTH_BY_CURRENCY holds
the depths; the extractor's --setup-deshoal-depth overrides them ('none' =
H0 = Hs, as the 'c' model uses it). The ratio H0/Hs is nearly constant
(1.06-1.10 over Tp 7-12 s), so a refitted k absorbs most of it:
evaluate_setup_model.py scores both H0 = Hs / Ks and H0 = Hs.

THE FALLBACK. A frame gets beta_f = C / k (C the --setup-coef in use:
0.106 for 0.037 at k = 0.35) when no DEM lies within SLOPE_MAX_DAYS of its
day, the swash band holds fewer than MIN_BAND_CELLS DEM cells, or the slope
found lies outside BETA_RANGE (a scarp, a gap edge, a flat patch: not a
foreshore). Its setup is then C * sqrt(H0 * L0): the station's setup, up to
the deshoaling. Every row the slope-aware model writes names the slope it
used (setup_beta_f) and the model (setup_model), and the extractor prints
how many frames fell back and why.

Self-test of this module: python3 test_setup_model.py
"""

import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

G = 9.81                      # m s^-2, as every module here used it

# Stockdon et al. (2006): setup Eq. 10, incident swash Eq. 11, infragravity swash Eq. 12
STOCKDON_SETUP_K = 0.35
STOCKDON_SWASH_INC = 0.75
STOCKDON_SWASH_IG = 0.06

SETUP_MODELS = ("c", "stockdon")

# Depth (m) the wave height of each wave currency (extract_elevation_contours.wave_columns: the
# '<kind>' of its tag) was measured or converted at, for deshoaling to Stockdon's H0. None = deep
# water, no deshoaling. marconi: hs_best, ADCP wh_4061 currency (21 m off Marconi; adcp_to_navd88.py:
# depth - water_level = 21.02 m); adcp_equiv: historical_forcing.py's ADCP-equivalent Hs; buoy: the
# raw NDBC buoy, taken as deep water.
DESHOAL_DEPTH_BY_CURRENCY = {"marconi": 21.0, "adcp_equiv": 21.0, "buoy": None}

# The foreshore slope: plane-fit window (cells; 10 m at the DEM's 2 m cells, foreshore_slope.py's
# default), the fewest DEM cells a swash band must hold (foreshore_slope.py --iterate stops below
# 20), the passes of the band <-> slope iteration, the slopes taken as a foreshore, and how far (days)
# a DEM's window centre may lie from a frame's day before the frame falls back (a week: the 25-26 Sep
# 2026 storm moved the beach by up to ~1 m in places within one).
SLOPE_WINDOW_CELLS = 5
MIN_BAND_CELLS = 20
SLOPE_PASSES = 3
BETA_RANGE = (0.01, 0.30)
SLOPE_MAX_DAYS = 7.0

# dem_from_contours.py --series-dir names its dated DEMs dem_<window end>_<days>d_dem.asc
SERIES_NAME = re.compile(r"^dem_(\d{4}-\d{2}-\d{2})_(\d+)d_dem\.asc$")


# ---------------------------------------------------------------------
# The wave scale and the station's model (bit-identical to the copies
# this replaced: same operations, same order)
# ---------------------------------------------------------------------

def hs_l0(hs, tp):
    """Hs * L0 (m^2), L0 = g Tp^2 / (2 pi), computed as (Hs * g * Tp^2) / (2 pi)."""
    return hs * G * tp ** 2 / (2 * np.pi)


def wave_scale(hs, tp):
    """sqrt(Hs * L0) (m): the wave scale of the setup (Stockdon et al. 2006), e.g. 10.0 m at Hs 1 m,
    Tp 8 s. Scalars or arrays (NaN where Hs or Tp is NaN)."""
    return np.sqrt(hs_l0(hs, tp))


def setup_c(hs, tp, coef):
    """The station's setup, C * sqrt(Hs * L0) (m): extract_elevation_contours.py --setup-coef."""
    return coef * wave_scale(hs, tp)


# ---------------------------------------------------------------------
# Stockdon et al. (2006) and linear-theory deshoaling
# ---------------------------------------------------------------------

def deep_water_wavelength(t):
    """L0 = g T^2 / (2 pi) (m)."""
    return G * t ** 2 / (2 * np.pi)


def stockdon_setup(h0, t0, beta, k=STOCKDON_SETUP_K):
    """Shoreline setup <eta> = k * beta_f * sqrt(H0 L0), Stockdon et al. (2006) Eq. 10 (k = 0.35)."""
    L0 = deep_water_wavelength(t0)
    return k * beta * np.sqrt(h0 * L0)


def stockdon_swash(h0, t0, beta):
    """Significant swash height S = sqrt(S_inc^2 + S_IG^2), S_inc = 0.75 beta_f sqrt(H0 L0), S_IG =
    0.06 sqrt(H0 L0): Stockdon et al. (2006) Eqs. 11-12."""
    L0 = deep_water_wavelength(t0)
    Sinc = STOCKDON_SWASH_INC * beta * np.sqrt(h0 * L0)
    SIG = STOCKDON_SWASH_IG * np.sqrt(h0 * L0)
    return np.sqrt(Sinc ** 2 + SIG ** 2)


def wavenumber(t, depth):
    """Linear-theory wavenumber k (rad/m) of period t (s) in water `depth` (m): the root of
    omega^2 = g k tanh(k h), from Fenton and McKee's (1990) explicit approximation refined by Newton
    steps (machine precision in four)."""
    t, h = np.asarray(t, float), np.asarray(depth, float)
    w2 = (2 * np.pi / t) ** 2
    k0h = w2 / G * h
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        kh = k0h / np.tanh(k0h ** 0.75) ** (2.0 / 3.0)
        for _ in range(4):
            th = np.tanh(kh)
            f = kh * th - k0h
            df = th + kh * (1.0 - th * th)
            kh = kh - f / df
    return kh / h


def shoaling_coefficient(t, depth):
    """Linear-theory shoaling coefficient Ks = H / H0 = sqrt(Cg0 / Cg) of period t (s) at `depth`
    (m), Cg0 = g T / (4 pi), Cg = n C, n = (1 + 2kh / sinh 2kh) / 2. 1 in deep water (kh > 10)."""
    t = np.asarray(t, float)
    k = wavenumber(t, depth)
    kh = k * np.asarray(depth, float)
    with np.errstate(invalid="ignore", over="ignore"):
        deep = kh > 10.0
        khs = np.where(deep, 1.0, kh)
        n = np.where(deep, 0.5, 0.5 * (1.0 + 2.0 * khs / np.sinh(2.0 * khs)))
        c = (2 * np.pi / t) / k
        cg0 = G * t / (4 * np.pi)
        ks = np.where(deep, 1.0, np.sqrt(cg0 / (n * c)))
    return ks


def deshoal(hs, tp, depth):
    """Deep-water equivalent H0 = Hs / Ks of a height measured at `depth` (m); depth None: Hs as given
    (deep water, or no deshoaling)."""
    if depth is None:
        return hs
    return hs / shoaling_coefficient(tp, depth)


def deshoal_depth(option, wave_tag):
    """The deshoaling depth for --setup-deshoal-depth `option` ('auto', 'none' or metres) and a wave
    file's currency tag (extract_elevation_contours.wave_columns: '<kind>:<height>/<period>').
    -> (depth or None, how). ValueError for an option that is none of these (or not above 0 m)."""
    opt = str(option if option is not None else "auto").strip().lower()
    if opt == "none":
        return None, "H0 = Hs (not deshoaled: --setup-deshoal-depth none)"
    if opt != "auto":
        d = float(opt)                       # ValueError for anything else
        if not d > 0:
            raise ValueError(f"deshoaling depth {option!r}: not a depth")
        return d, f"H0 = Hs / Ks, deshoaled from {d:g} m (--setup-deshoal-depth)"
    kind = str(wave_tag or "").split(":")[0]
    d = DESHOAL_DEPTH_BY_CURRENCY.get(kind)
    if d is None:
        return None, (f"H0 = Hs (the {kind or 'unknown'} wave currency is taken as deep water)")
    return d, f"H0 = Hs / Ks, deshoaled from {d:g} m (the depth of the {kind} wave currency)"


# ---------------------------------------------------------------------
# The foreshore slope of the evolving DEM
# ---------------------------------------------------------------------

def slope_grid(dem, cell, window=SLOPE_WINDOW_CELLS):
    """
    Slope magnitude of each cell along the direction of steepest descent: the gradient of a plane
    fitted by least squares to the finite cells of the window about it, NaN where the cell is NaN or
    the window holds fewer than six cells (or only collinear ones). The same as
    foreshore_slope.surface_slope (to ~1e-12; test_setup_model.py checks it), as array sums over the
    window's offsets instead of one least-squares solve per cell: a year of dated DEMs then takes
    seconds, not minutes.
    """
    dem = np.asarray(dem, float)
    rows, cols = dem.shape
    half = max(1, int(window // 2))
    pad = np.full((rows + 2 * half, cols + 2 * half), np.nan)
    pad[half:half + rows, half:half + cols] = dem
    n = np.zeros(dem.shape)
    sx, sy, sxx, syy, sxy, sz, sxz, syz = (np.zeros(dem.shape) for _ in range(8))
    for dy in range(-half, half + 1):
        for dx in range(-half, half + 1):
            z = pad[half + dy:half + dy + rows, half + dx:half + dx + cols]
            m = np.isfinite(z)
            zz = np.where(m, z, 0.0)
            mf = m.astype(float)
            n += mf
            sx += mf * dx
            sy += mf * dy
            sxx += mf * dx * dx
            syy += mf * dy * dy
            sxy += mf * dx * dy
            sz += zz
            sxz += zz * dx
            syz += zz * dy
    out = np.full(dem.shape, np.nan)
    ok = np.isfinite(dem) & (n >= 6)
    nn = np.where(ok, n, 1.0)
    cxx = sxx - sx * sx / nn
    cyy = syy - sy * sy / nn
    cxy = sxy - sx * sy / nn
    cxz = sxz - sx * sz / nn
    cyz = syz - sy * sz / nn
    det = cxx * cyy - cxy * cxy
    ok &= det > 1e-9 * np.maximum(cxx * cyy, 1e-12)
    d = np.where(ok, det, 1.0)
    a = (cyy * cxz - cxy * cyz) / d
    b = (cxx * cyz - cxy * cxz) / d
    out[ok] = np.hypot(a[ok], b[ok]) / float(cell)
    return out


def read_dem(path):
    """An ESRI ASCII DEM (dem_from_contours.py) -> (grid north to south, cell size), NaN = nodata."""
    from asc_to_geotiff import read_asc
    grid, head = read_asc(path)
    return grid, float(head["cellsize"])


class ForeshoreSlope:
    """The foreshore slope of one DEM over any elevation band: the slope grid is computed once, its
    cells sorted by DEM elevation, and a band's slope is the median slope of the cells inside it."""

    def __init__(self, dem, cell, name="", window=SLOPE_WINDOW_CELLS):
        dem = np.asarray(dem, float)
        s = slope_grid(dem, cell, window)
        ok = np.isfinite(dem) & np.isfinite(s)
        o = np.argsort(dem[ok], kind="mergesort")
        self.z, self.s, self.name = dem[ok][o], s[ok][o], name
        self.cells = int(ok.sum())

    @classmethod
    def from_file(cls, path, window=SLOPE_WINDOW_CELLS):
        grid, cell = read_dem(path)
        return cls(grid, cell, Path(path).name, window)

    def band(self, lo, hi, min_cells=MIN_BAND_CELLS):
        """(median slope of the DEM cells with lo <= z <= hi, cells) or (None, cells) when fewer than
        min_cells."""
        a, b = np.searchsorted(self.z, lo, "left"), np.searchsorted(self.z, hi, "right")
        if b - a < min_cells:
            return None, int(b - a)
        return float(np.median(self.s[a:b])), int(b - a)


class EvolvingSlope:
    """
    The DEM each frame's slope is measured on. `source`: a DEM (<stem>_dem.asc, or its stem), used for
    every frame; or a folder of dated DEMs (dem_from_contours.py --series-dir: dem_<window end>_<days>d_
    dem.asc), where a frame takes the DEM whose window centre lies nearest its UTC day (the DEM windows
    are UTC days, dem_from_contours.py --last-days), if within max_days. Each DEM's slope grid is
    computed the first time a frame needs it.
    """

    def __init__(self, source, max_days=SLOPE_MAX_DAYS, window=SLOPE_WINDOW_CELLS):
        self.max_days, self.window, self._cache = float(max_days), window, {}
        self.series = []                                   # [(centre epoch, path)]
        self.single = None
        p = Path(str(source))
        if p.is_dir():
            for f in sorted(p.iterdir()):
                m = SERIES_NAME.match(f.name)
                if m:
                    end = datetime.strptime(m.group(1), "%Y-%m-%d").replace(tzinfo=timezone.utc).timestamp()
                    days = int(m.group(2))
                    # the window's days run end-days+1 .. end; its centre is the middle of them (noon UTC)
                    centre = end + 43200.0 - (days - 1) / 2.0 * 86400.0
                    self.series.append((centre, f))
            self.series.sort()
            self.what = f"{len(self.series)} dated DEM(s) in {p}"
        else:
            # a DEM file as given, else the <stem>_dem.asc of a stem
            f = p if (p.is_file() or p.name.endswith("_dem.asc")) else Path(str(p) + "_dem.asc")
            self.single = f if f.is_file() else None
            self.what = str(f) + ("" if self.single else " (NOT FOUND)")
        self.centres = np.array([c for c, _ in self.series], float)

    def _load(self, path):
        if path not in self._cache:
            try:
                self._cache[path] = ForeshoreSlope.from_file(path, self.window)
            except (OSError, ValueError, KeyError) as exc:
                print(f"WARNING: cannot read the DEM {path} for the foreshore slope ({exc})")
                self._cache[path] = None
        return self._cache[path]

    def for_epoch(self, epoch):
        """(ForeshoreSlope or None, why None)."""
        if self.single is not None:
            fs = self._load(self.single)
            return fs, ("" if fs else "its DEM could not be read")
        if not self.series:
            return None, "no DEM"
        c = self.centres
        i = int(np.argmin(np.abs(c - float(epoch))))
        if abs(c[i] - float(epoch)) > self.max_days * 86400.0:
            return None, f"no DEM within {self.max_days:g} days"
        fs = self._load(self.series[i][1])
        return fs, ("" if fs else "its DEM could not be read")


def stockdon_frame(still_water, hs, tp, slope, fallback_beta, k=STOCKDON_SETUP_K, depth=None,
                   passes=SLOPE_PASSES, why_none="no DEM"):
    """
    One frame's slope-aware setup. still_water: the frame's water level (m NAVD88); slope: the
    ForeshoreSlope of its DEM, or None (why_none says why); fallback_beta: the slope used where the
    DEM gives none (C / k). The band <-> slope iteration of foreshore_slope.py --iterate: the band is
    still water + setup +/- S/2 at the current slope; a pass that finds too few cells, or a slope
    outside BETA_RANGE, keeps the previous slope. -> (setup m, beta_f, source: 'dem' or the reason
    for the fallback).
    """
    h0 = deshoal(float(hs), float(tp), depth)
    beta = float(fallback_beta)
    src = why_none if slope is None else "too few DEM cells in the swash band"
    if slope is not None:
        for _ in range(passes):
            setup = float(stockdon_setup(h0, tp, beta, k))
            half = float(stockdon_swash(h0, tp, beta)) / 2.0
            b, _n = slope.band(still_water + setup - half, still_water + setup + half)
            if b is None:
                break
            if not BETA_RANGE[0] <= b <= BETA_RANGE[1]:
                if src != "dem":
                    src = f"DEM slope {b:.3f} outside {BETA_RANGE[0]}-{BETA_RANGE[1]}"
                break
            beta, src = b, "dem"
    return float(stockdon_setup(h0, tp, beta, k)), beta, src


def model_tag(k, depth):
    """The setup_model value of a slope-aware row, e.g. 'stockdon2006:k=0.35:H0=Hs/Ks@21m'."""
    return f"stockdon2006:k={k:g}:" + (f"H0=Hs/Ks@{depth:g}m" if depth is not None else "H0=Hs")

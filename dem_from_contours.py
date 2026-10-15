#!/usr/bin/env python3
"""
Intertidal DEM From Waterline Contours
-----------------------------------------
Grids georectified waterline points into a beach elevation model.

THE IDEA: each detected waterline point is a place where the water
surface met the sand, so the BEACH elevation there equals the WATER
elevation, which GNSS-R measured independently. Over a tidal cycle the
waterline sweeps across the intertidal zone and traces out the beach
profile. This grids those traces into a surface.

WHAT THIS DEM IS AND IS NOT:
  * It covers the INTERTIDAL ZONE ONLY -- between the lowest and
    highest water levels observed. There is no information above the
    high-water line or below the low-water line, and none is invented.
  * Cells are filled ONLY where waterlines actually passed. Gaps are
    left as nodata rather than interpolated, because a smooth surface
    across an unmeasured gap looks like data and is not. Use
    --fill-gaps if you need a continuous surface and accept that.
  * Every cell carries a point count and an elevation spread, written
    as companion grids. A cell built from 2 points spanning 0.4 m is
    not the same measurement as one from 200 points spanning 0.03 m,
    and the DEM alone cannot show that difference.

SOURCES OF ERROR, roughly in order of size:
  * Waterline detection: the largest term, and it varies across the
    frame. The far field is worse than the near field.
  * Georectification: 0.33 m horizontal RMS against surveyed control,
    as of the Nov 2025 calibration. Unverified since.
  * Water level: GNSS-R, a measurement rather than a model.
  * Wave setup: the detected edge sits ABOVE the still-water line by an
    amount that grows with wave height. When the contours carry
    beach_elevation_navd88 (water level + C*sqrt(Hs*L0), C from
    extract_elevation_contours.py --setup-coef; the station's since
    6 Oct 2026, C = SETUP_COEF of waterline_timex_cron.sh), that is what
    is gridded, so the setup is included; how well depends on C
    (--fit-setup checks it against repeat crossings). Without it the edge
    is given the still-water elevation and the DEM reads LOW where
    rough-water frames land; --max-hs leaves the roughest frames out,
    using the offshore wave record (fetch_buoy_waves.py).

Usage:
    python3 dem_from_contours.py contour_points_ground.csv dem_out \\
        [--cell 2.0] [--min-points 3] [--camera c1|c2|both]
        [--max-spread 0.5] [--fill-gaps] [--start-date ...] [--end-date ...]

Writes dem_out_{dem,spread,count}.asc (and _source.asc with
--interpolate-edge), dem_out_info.json (dates, frames per camera and
what the filters left out -- what the page needs and the grids do not
hold; with --max-day-offset also each camera-day left out with its
offset, and the day test's first-pass and final offsets) and the page
dem_out_dem.png, drawn by dem_figure.py. To redraw the page from those
files without rebuilding: python3 dem_figure.py dem_out

THE DAY TEST (--max-day-offset) is iterative: the camera-day furthest
over the limit from the DEM of the other days is left out, the others
are measured again without it, and so on until none is over; a day left
out on the way is taken back if it lies within the limit of the DEM of
the kept days and its return puts no kept day over the limit (else it
stays out: a refused take-back). One pass let the days it rejected shape
the DEM it judged
the others by (day_consistency() has the details). A day is left out
only where the test can tell which day is off; days over the limit it
cannot judge (two days that only see each other, the last two days)
keep the one-pass verdict, and a WARNING names those the final DEM leaves
out, or keeps over the limit or untested. A camera whose day offsets run
one way with time gets a WARNING too: real change, or an error drifting
with time (pointing, lens, a detection bias that follows the tides or
waves), and the test then trims an end of it. The rule: 5 or more
testable days, a range over the limit, Spearman's rho of offset against
date above 0.8 and, up to 8 days, one that days in random order reach
at most 2 % of the time (5 days need 1.00, 6 days 0.89). When the one
day the test left out of a camera's days is its first or last, the other
days must still run the same way (rho above 0.8, however few they are);
if on their own they would not warn, the WARNING says the trend rests on
that day.

THE PAGE IS ALWAYS THIS RUN'S. The old dem_out_dem.png is removed first
(unless --no-plot); a window with too few cells, or no points at all,
gets a short page saying so, with the dates and what the filters left
out. Exit status: 0 built (page drawn; a short one if too few cells),
4 no points in the window, or every camera-day rejected by the
--max-day-offset test (short page drawn, no grids written), 3 grids
written but the page could not be drawn (no PNG), 1 any other failure
(no PNG). The email attaches the PNG by name, so a stale one would pass
for current.
"""

import sys
import csv
import json
import argparse
from pathlib import Path

import numpy as np


def carries_setup(text):
    """True when a setup_correction_m cell holds a wave setup: blank is a frame with no
    wave record, 0 a build with C = 0 -- in both the elevation is the water level alone."""
    try:
        return float(text) != 0.0
    except (TypeError, ValueError):
        return False


def load_points(path, camera=None, start_date=None, end_date=None, max_hs=None, exclude=None,
                stats=None):
    """Reads georectified contour points. Rows without ground coordinates are skipped.
    `exclude`: regular expressions; frames whose source_file matches any are left out.
    `stats`: a dict to receive the number of frames each filter left out, and the
    share of points whose elevation includes the wave setup (for the page)."""
    import re
    exclude = [re.compile(x) for x in (exclude or [])]
    excluded = set()
    E, N, Z, cams, dates, frames = [], [], [], [], [], []
    ucol = []                     # pixel column: the order of a line's points (cross-shore rule)
    missing_ground = 0
    rough_frames, unknown_hs_frames = set(), set()
    with_setup = 0
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        if "easting_utm19" not in (reader.fieldnames or []):
            print("ERROR: no 'easting_utm19' column. Run georectify.py on the contour file first.")
            sys.exit(1)
        has_source = "source_file" in (reader.fieldnames or [])
        if max_hs is not None and "offshore_hs_m" not in (reader.fieldnames or []):
            print("ERROR: --max-hs needs an 'offshore_hs_m' column. Run "
                  "extract_elevation_contours.py with --waves first.")
            sys.exit(1)
        if not has_source:
            print("WARNING: no 'source_file' column -- cannot tell which points share a "
                  "frame, so every point is counted as an independent sample.")
        for r in reader:
            if not r.get("easting_utm19") or not r.get("northing_utm19"):
                missing_ground += 1
                continue
            if camera and camera != "both" and r["camera"] != camera:
                continue
            if exclude and any(x.search(r.get("source_file", "")) for x in exclude):
                excluded.add(r.get("source_file", ""))
                continue
            day = r.get("capture_time_utc", "")[:10]
            if start_date and day < start_date:
                continue
            if end_date and day > end_date:
                continue
            if max_hs is not None:
                hs = r.get("offshore_hs_m", "")
                if hs == "":
                    unknown_hs_frames.add(r.get("source_file", ""))
                elif float(hs) > max_hs:
                    rough_frames.add(r.get("source_file", ""))
                    continue
            E.append(float(r["easting_utm19"]))
            N.append(float(r["northing_utm19"]))
            # Beach elevation = water level + setup when the contours were
            # extracted with --setup-coef; otherwise the water level.
            beach = r.get("beach_elevation_navd88")
            with_setup += bool(beach) and carries_setup(r.get("setup_correction_m"))
            Z.append(float(beach or r["tide_elevation_navd88"]))
            cams.append(r["camera"])
            dates.append(day)
            frames.append(r["source_file"] if has_source else f"__point{len(frames)}")
            try:
                ucol.append(float(r.get("pixel_column") or "nan"))
            except ValueError:
                ucol.append(float("nan"))
    if exclude:
        print(f"Excluded by name  : {len(excluded)} frame(s) (--exclude)")
    if max_hs is not None:
        print(f"Wave filter       : {len(rough_frames)} frame(s) left out, offshore Hs > "
              f"{max_hs} m; {len(unknown_hs_frames)} frame(s) with no wave record kept")
    if stats is not None:
        stats.update(excluded_frames=len(excluded), rough_frames=len(rough_frames),
                     setup_share=with_setup / len(Z) if Z else 0.0,
                     pixel_column=np.array(ucol, float))
    return (np.array(E), np.array(N), np.array(Z),
            np.array(cams), np.array(dates), np.array(frames), missing_ground)


def latest_date(path, camera=None):
    """Latest UTC capture date (YYYY-MM-DD) among georectified rows."""
    latest = None
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if not r.get("easting_utm19"):
                continue
            if camera and camera != "both" and r["camera"] != camera:
                continue
            d = r.get("capture_time_utc", "")[:10]
            if d and (latest is None or d > latest):
                latest = d
    return latest


def build_grid(E, N, Z, cell, min_points, max_spread, frames=None):
    """
    Bins points into cells and takes the MEDIAN elevation of each.

    Median rather than mean: a cell can catch an outlier from a single
    bad detection, and one wild value would drag a mean while barely
    moving a median. With repeat tidal crossings most cells hold many
    samples, so the median is well determined.

    Samples are counted per FRAME, not per point. In the far field one
    waterline puts many adjacent pixel columns into the same cell, all
    with the same elevation (one water level per frame). Counted as
    separate points they satisfied min_points with zero spread, so a
    cell resting on a single detection looked perfectly repeatable.
    Each frame's points in a cell are therefore collapsed to one sample
    first; count, min_points and spread then refer to independent
    crossings. `frames` gives the frame of each point; None falls back
    to treating every point as its own frame.
    """
    e0 = np.floor(E.min() / cell) * cell
    n0 = np.floor(N.min() / cell) * cell
    ncols = int(np.ceil((E.max() - e0) / cell)) + 1
    nrows = int(np.ceil((N.max() - n0) / cell)) + 1

    col = ((E - e0) / cell).astype(int)
    row = ((N - n0) / cell).astype(int)
    flat = row * ncols + col

    if frames is not None:
        # One sample per (cell, frame): the median of that frame's
        # points in the cell.
        _, frame_id = np.unique(frames, return_inverse=True)
        key = flat.astype(np.int64) * (int(frame_id.max()) + 1) + frame_id
        korder = np.argsort(key, kind="stable")
        key_s, kz = key[korder], Z[korder]
        kedges = np.flatnonzero(np.diff(key_s)) + 1
        kstarts = np.concatenate([[0], kedges])
        kends = np.concatenate([kedges, [len(key_s)]])
        Z = np.array([np.median(kz[a:b]) for a, b in zip(kstarts, kends)])
        flat = flat[korder][kstarts]

    order = np.argsort(flat)
    flat_s, Z_s = flat[order], Z[order]
    edges = np.flatnonzero(np.diff(flat_s)) + 1
    starts = np.concatenate([[0], edges])
    ends = np.concatenate([edges, [len(flat_s)]])

    dem = np.full(nrows * ncols, np.nan)
    count = np.zeros(nrows * ncols, dtype=int)
    spread = np.full(nrows * ncols, np.nan)

    for s, e in zip(starts, ends):
        idx = flat_s[s]
        vals = Z_s[s:e]
        count[idx] = len(vals)
        if len(vals) >= min_points:
            dem[idx] = np.median(vals)
            # Spread as the 16th-84th percentile range -- a robust
            # stand-in for +/-1 sigma that a couple of outliers cannot
            # inflate the way a standard deviation can.
            spread[idx] = (np.percentile(vals, 84) - np.percentile(vals, 16)) if len(vals) > 2 \
                else float(vals.max() - vals.min())

    if max_spread and max_spread > 0:
        with np.errstate(invalid="ignore"):      # NaN spread: no cell (numpy < 1.18 warns)
            too_noisy = np.isfinite(spread) & (spread > max_spread)
        dem[too_noisy] = np.nan

    return (dem.reshape(nrows, ncols), count.reshape(nrows, ncols),
            spread.reshape(nrows, ncols), e0, n0, ncols, nrows)


# ---------------------------------------------------------------------
# Frame-level analyses: setup-coefficient fit and day consistency.
#
# Both work on (cell, frame) PAIRS -- one sample per frame per cell, as
# build_grid counts them -- on a fixed grid origin, so a subset of the
# data can be compared against the rest cell by cell.
# ---------------------------------------------------------------------

def load_frame_info(path):
    """Per source_file: capture epoch, water level, offshore Hs and Tp."""
    info = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            k = r.get("source_file")
            if k in info or k is None:
                continue
            num = lambda v: float(v) if v not in ("", None) else np.nan
            info[k] = {"epoch": num(r.get("capture_epoch")),
                       "tide": num(r.get("tide_elevation_navd88")),
                       "hs": num(r.get("offshore_hs_m")),
                       "tp": num(r.get("offshore_tp_s")),
                       "setup": num(r.get("setup_correction_m"))}
    return info


def local_day(epoch):
    import pandas as pd
    return pd.Timestamp(int(epoch), unit="s", tz="UTC").tz_convert("America/New_York").strftime("%Y-%m-%d")


def frame_pairs(E, N, frames, cell):
    """Unique (cell, frame) pairs on a fixed origin. Returns (pair_cell, pair_frame_name)."""
    e0 = np.floor(E.min() / cell) * cell
    n0 = np.floor(N.min() / cell) * cell
    ncols = int(np.ceil((E.max() - e0) / cell)) + 1
    flat = ((N - n0) / cell).astype(np.int64) * ncols + ((E - e0) / cell).astype(np.int64)
    names, fid = np.unique(frames, return_inverse=True)
    key = np.unique(flat * len(names) + fid)
    return key // len(names), names[key % len(names)]


def cell_groups(pair_cell):
    order = np.argsort(pair_cell, kind="stable")
    c = pair_cell[order]
    edges = np.flatnonzero(np.diff(c)) + 1
    return order, c, np.r_[0, edges], np.r_[edges, len(c)]


def median_spread_by_cell(groups, vals, min_frames):
    """Median of the 16-84 percentile spread over cells with >= min_frames samples."""
    order, c, starts, ends = groups
    v = vals[order]
    spreads = []
    for a, b in zip(starts, ends):
        x = v[a:b]
        x = x[np.isfinite(x)]
        if len(x) >= min_frames:
            spreads.append(np.percentile(x, 84) - np.percentile(x, 16))
    return float(np.median(spreads)) if spreads else np.nan, len(spreads)


# --fit-setup: the C range whose median spread lies within FIT_SETUP_NEAR_M of the best is
# printed (the near-best band); FIT_SETUP_FLAT_WIDTH or wider is a flat minimum. On one week
# of lines the curve is noisy at the few-mm level (on the filtered 29 Sep - 5 Oct 2026 C = 0
# lines, --max-hs 1.5, it wanders 0.248-0.256 m over C 0-0.04), so neither the best C nor the
# exact band width means much by itself: a C is offered to apply only when C = 0 AND the C in
# use both lie OUTSIDE the near-best band (each spreads more than FIT_SETUP_NEAR_M above the best
# and lies outside the band's range) and the band is narrower than FIT_SETUP_FLAT_WIDTH.
# Otherwise the fit says which of the two it cannot tell apart from its best on this window.
FIT_SETUP_NEAR_M = 0.005
FIT_SETUP_FLAT_WIDTH = 0.02
STATION_CRON = Path(__file__).resolve().parent / "waterline_timex_cron.sh"


def carried_setup_coef(info):
    """The C the contours carry: median of setup_correction_m / sqrt(Hs*L0) over the frames
    with a non-zero setup; None when no frame carries one (a C = 0 build, or no setup column)."""
    cs = []
    for v in info.values():
        if np.isfinite(v["setup"]) and v["setup"] != 0.0 and v["hs"] > 0 and v["tp"] > 0:
            cs.append(v["setup"] / np.sqrt(v["hs"] * 9.81 * v["tp"] ** 2 / (2 * np.pi)))
    return round(float(np.median(cs)), 4) if cs else None


def station_setup_coef(path=STATION_CRON):
    """SETUP_COEF of the station cron beside this script ('' = no setup = 0); None if unreadable."""
    import re
    try:
        for line in Path(path).read_text().splitlines():
            m = re.match(r'\s*SETUP_COEF=["\']?([0-9.]*)["\']?\s*(#.*)?$', line)
            if m:
                return float(m.group(1)) if m.group(1) else 0.0
    except (OSError, ValueError):
        pass
    return None


def setup_fit_verdict(valid, spread_at, c_use, use_from=""):
    """
    The lines --fit-setup prints under its curve. valid: [(C, median spread, cells)] on the
    search grid; spread_at(C): the median spread at any C (for C = 0 and the C in use, which
    need not lie on the grid); c_use: the C in use (None or 0: compare with C = 0 only), and
    use_from where that value came from.
    """
    best = min(valid, key=lambda r: r[1])
    near = [r[0] for r in valid if r[1] <= best[1] + FIT_SETUP_NEAR_M]
    lo, hi = min(near), max(near)
    flat = hi - lo >= FIT_SETUP_FLAT_WIDTH
    out = [f"   C {lo:+.3f} to {hi:+.3f}: median spread within {FIT_SETUP_NEAR_M * 1000:.0f} mm of the "
           f"best (the near-best band)" + (" -- a FLAT minimum: this fit alone does not pin C" if flat else "")]
    # a C is told apart from the best only when it lies outside the near-best band: more than 5 mm
    # above the best AND not between two C that are within it (the curve is noisy at the mm level)
    s0 = spread_at(0.0)
    same0 = s0 <= best[1] + FIT_SETUP_NEAR_M or lo <= 0.0 <= hi
    su = spread_at(c_use) if c_use else None
    same_use = su is not None and (su <= best[1] + FIT_SETUP_NEAR_M or lo <= c_use <= hi)
    if c_use:
        out.append(f"   C in use {c_use}" + (f" ({use_from})" if use_from else "")
                   + f": median spread {su:.3f} m ({su - best[1]:+.3f} m from the best)")
    if best[0] in (valid[0][0], valid[-1][0]):
        out.append("   WARNING: best C is at the edge of the search range -- the data do not constrain it; "
                   "do not apply.")
    elif same0 and same_use:
        out.append(f"   Neither C = 0 nor the C in use ({c_use}) can be told from the best on this window "
                   f"(both inside the near-best band): this fit gives no reason to change C.")
    elif same_use:
        out.append(f"   Not distinguishable from the C in use ({c_use}) on this window (inside the near-best "
                   f"band): keep it.")
    elif same0:
        out.append("   Not distinguishable from C = 0 on this window (inside the near-best band): the data do "
                   "not favour a positive setup correction here. One week does not pin C: refit over weeks of "
                   "different waves before changing the C in use.")
    elif best[0] <= 0:
        out.append("   The data do not favour a positive setup correction.")
    elif flat:
        out.append("   Do not apply on this window alone: refit over weeks of different waves (or keep the "
                   "C in use).")
    else:
        out.append("   Beats C = 0" + (f" and the C in use ({c_use})" if c_use else "")
                   + f" by more than {FIT_SETUP_NEAR_M * 1000:.0f} mm on this window; check it on other weeks "
                   "of different waves before adopting it.")
        out.append(f"   To apply: extract_elevation_contours.py ... --waves ... --setup-coef {best[0]}")
    return out


def fit_setup_coefficient(pair_cell, pair_frame, info, min_frames, coefs=None):
    """
    Finds C in beach elevation = water level + C*sqrt(Hs*L0) that makes
    repeat crossings of each cell agree best (smallest median spread).
    Only frames with both Hs and Tp take part, at every C, so the
    comparison is like for like. coefs: the C values to try (default
    -0.04 to +0.12 in steps of 0.002).
    """
    tide = np.array([info[f]["tide"] for f in pair_frame])
    hs = np.array([info[f]["hs"] for f in pair_frame])
    tp = np.array([info[f]["tp"] for f in pair_frame])
    phi = np.sqrt(hs * 9.81 * tp ** 2 / (2 * np.pi))
    ok = np.isfinite(phi) & np.isfinite(tide)
    groups = cell_groups(pair_cell[ok])
    if coefs is None:
        coefs = np.round(np.arange(-0.04, 0.1201, 0.002), 3)
    result = [(c,) + median_spread_by_cell(groups, tide[ok] + c * phi[ok], min_frames)
              for c in coefs]
    return result, int(ok.sum()), len(set(pair_frame[ok]))


class LeaveOneOut:
    """
    Offsets of groups of (cell, frame) samples -- camera-days or frames, the `keys` -- from the
    DEM of other groups: for key k, the median over k's samples of (sample - the median of the
    reference's samples in that cell), over the cells where the reference has >= min_frames
    samples. The reference is every other key, or a chosen subset (offset(k, ref)), so the
    iterative day test can recompute the offsets after each day it leaves out.

    The samples are sorted once (by cell, then value); a cell median is the middle value or the
    mean of the two middle values, as np.median gives it, so the offsets do not depend on the
    order of the samples and equal those of the one-pass loop this replaced. offsets() gives the
    offsets of many keys at once, without a loop over the keys (the iterative day test recomputes
    every kept day after each one it leaves out: on a synthetic year of camera-days, 730 days and
    167 left out, 668k samples, one key at a time took 23 min on the station's Python 3.8 / numpy
    1.17, this about a minute; the time grows with samples x rounds, and at the station's density,
    1.5M samples for 730 days with 98 left out, it took about 2.5 min of CPU there, 3.6-5.3 min on
    a loaded machine, against 5-8 s for one pass);
    offset() gives one key's, the same to the bit.
    """

    def __init__(self, pair_cell, pair_vals, pair_key, min_frames, min_pairs=30):
        self.names, kid = np.unique(pair_key, return_inverse=True)
        kid = np.asarray(kid).ravel()
        order = np.lexsort((pair_vals, pair_cell))          # by cell, then value (stable)
        self.c, self.v, self.k = pair_cell[order], pair_vals[order], kid[order]
        self.min_frames, self.min_pairs = min_frames, min_pairs
        self.index = {n: i for i, n in enumerate(self.names.tolist())}
        # the cells each key has samples in (sorted, unique)
        ko = np.lexsort((self.c, self.k))
        kk, kc = self.k[ko], self.c[ko]
        new = np.r_[True, (kk[1:] != kk[:-1]) | (kc[1:] != kc[:-1])] if len(kk) else np.zeros(0, bool)
        uk, uc = kk[new], kc[new]
        b = np.searchsorted(uk, np.arange(len(self.names) + 1))
        self.cells = [uc[b[i]:b[i + 1]] for i in range(len(self.names))]
        # for offsets(): each sample's cell (0..), each cell's first sample, and the (cell, key)
        # groups -- every sample of one key in one cell -- with each sample's place in its group
        nk, n = max(len(self.names), 1), len(self.c)
        newc = np.r_[True, self.c[1:] != self.c[:-1]] if n else np.zeros(0, bool)
        self.ci = np.cumsum(newc) - 1
        self.cstart = np.flatnonzero(newc)
        gid, g = np.unique(self.ci.astype(np.int64) * nk + self.k, return_inverse=True)
        self.g = np.asarray(g).ravel()
        self.gcell, self.gkey = gid // nk, gid % nk
        self.gsize = np.bincount(self.g, minlength=len(gid))
        self.gorder = np.argsort(self.g, kind="mergesort")  # by group, then value within it
        self.g_sorted = self.g[self.gorder]
        seq = np.empty(n, np.int64)
        seq[self.gorder] = np.arange(n) - np.repeat(np.r_[0, np.cumsum(self.gsize)[:-1]], self.gsize)
        self.seq = seq

    def offsets(self, ref=None, keys=None):
        """{key: (offset, comparable samples)} for `keys` (a boolean per key; default all), each
        against the DEM of the `ref` keys (default all) without its own samples -- what offset(k,
        ref) gives for each k, to the bit, in one pass over the samples. NaN offset with fewer than
        min_pairs comparable samples."""
        nk = len(self.names)
        ref = np.ones(nk, bool) if ref is None else np.asarray(ref, bool)
        keys = np.ones(nk, bool) if keys is None else np.asarray(keys, bool)
        out = {}
        if not len(self.c):
            return {k: (np.nan, 0) for i, k in enumerate(self.names.tolist()) if keys[i]}
        inref = ref[self.k]
        cr = np.r_[0, np.cumsum(inref)]
        before = cr[self.cstart]                          # reference samples before each cell
        ncell = np.r_[before[1:], cr[-1]] - before        # ... and in it
        rv = self.v[inref]                                # still by cell, then value
        # the cell's reference without the group's own samples (when its key is a reference key)
        excl = np.where(ref[self.gkey], self.gsize, 0)
        left = ncell[self.gcell] - excl
        gg = np.flatnonzero(keys[self.gkey] & (left >= max(self.min_frames, 1)))
        r = left[gg]
        jlo, jhi = (r - 1) // 2, r // 2                   # the middle one(s) of what is left
        ex = excl[gg] > 0
        if ex.any():
            # the j-th of the cell's reference samples is the r-th left once the group's own are
            # out, where j = r + (the group's own samples before it): own sample i (in value order)
            # at place p_i has p_i - i others before it, so it comes before the r-th if that is <= r
            sel = np.zeros(len(self.gsize), bool)
            sel[gg[ex]] = True
            pos = self.gorder[sel[self.g_sorted]]
            p = cr[pos + 1] - 1 - before[self.ci[pos]]    # place among the cell's reference samples
            big = int(ncell.max()) + 1
            q = self.g[pos].astype(np.int64) * big + (p - self.seq[pos])
            gx = gg[ex].astype(np.int64) * big
            at0 = np.searchsorted(q, gx, "left")
            jlo[ex] += np.searchsorted(q, gx + jlo[ex], "right") - at0
            jhi[ex] += np.searchsorted(q, gx + jhi[ex], "right") - at0
        base = before[self.gcell[gg]]
        med = np.full(len(self.gsize), np.nan)
        med[gg] = (rv[base + jlo] + rv[base + jhi]) / 2
        ok = np.zeros(len(self.gsize), bool)
        ok[gg] = True
        s = np.flatnonzero(ok[self.g])                    # the samples of `keys` that are compared
        res = self.v[s] - med[self.g[s]]
        sk = self.k[s]
        o = np.lexsort((res, sk))
        res, sk = res[o], sk[o]
        cnt = np.bincount(sk, minlength=nk)
        st = np.r_[0, np.cumsum(cnt)[:-1]]
        for i, name in enumerate(self.names.tolist()):
            if keys[i]:
                m = int(cnt[i])
                out[name] = ((float((res[st[i] + (m - 1) // 2] + res[st[i] + m // 2]) / 2)
                              if m >= self.min_pairs else np.nan), m)
        return out

    def near(self, i):
        """Positions (in the sorted samples) of every sample in the cells key i has samples in."""
        cells = self.cells[i]
        lo = np.searchsorted(self.c, cells, "left")
        n = np.searchsorted(self.c, cells, "right") - lo
        return np.repeat(lo - np.r_[0, np.cumsum(n)[:-1]], n) + np.arange(int(n.sum()))

    def offset(self, name, ref=None):
        """(median offset of key `name` from the DEM of the reference keys, comparable samples).
        ref: a boolean per key (in self.names order) -- the keys whose samples make the DEM;
        None = all. The key's own samples never do. NaN offset with fewer than min_pairs."""
        i = self.index[name]
        pos = self.near(i)
        c, v, kk = self.c[pos], self.v[pos], self.k[pos]
        mine = kk == i
        use = ~mine if ref is None else (~mine & ref[kk])
        rc, rv = c[use], v[use]                           # still sorted by cell, then value
        res = np.zeros(0)
        if len(rc):
            starts = np.flatnonzero(np.r_[True, rc[1:] != rc[:-1]])
            n = np.diff(np.r_[starts, len(rc)])
            ok = n >= self.min_frames
            uc = rc[starts][ok]
            med = ((rv[starts + (n - 1) // 2] + rv[starts + n // 2]) / 2)[ok]
            if len(uc):
                oc, ov = c[mine], v[mine]
                at = np.minimum(np.searchsorted(uc, oc), len(uc) - 1)
                hit = uc[at] == oc
                res = ov[hit] - med[at[hit]]
        return (float(np.median(res)) if len(res) >= self.min_pairs else np.nan, len(res))

    def judge(self, name, ref, limit):
        """Can the test tell whether key `name` or the DEM it is measured against is off? Where it
        is compared, that DEM (the `ref` keys without `name`) holds in each cell either two or more
        other keys -- a majority can speak -- or one key alone. One key alone judges it only when
        that key is itself checked by the others: its offset from the DEM of the `ref` keys without
        `name` is testable and within `limit` (it 'vouches'); a key that only `name` reaches, or
        that is itself off, cannot say which of the two is wrong. `name` can be judged when at
        least half its comparable samples lie in cells judged so.
        -> {judged, share (of its comparable samples in cells that judge it), samples (comparable),
        alone: [the keys that alone make its DEM in cells and do not vouch, by name]}"""
        i = self.index[name]
        out = {"judged": False, "share": 0.0, "samples": 0, "alone": []}
        pos = self.near(i)
        c, kk = self.c[pos], self.k[pos]
        mine = kk == i
        use = ~mine & ref[kk]
        rc, rk = c[use], kk[use]
        if not len(rc):
            return out
        starts = np.flatnonzero(np.r_[True, rc[1:] != rc[:-1]])
        ok = np.diff(np.r_[starts, len(rc)]) >= self.min_frames
        uc = rc[starts]
        o = np.lexsort((rk, rc))
        rc2, rk2 = rc[o], rk[o]
        one = np.r_[True, (rc2[1:] != rc2[:-1]) | (rk2[1:] != rk2[:-1])]       # one per (cell, key)
        pc, pk = rc2[one], rk2[one]
        lo = np.searchsorted(pc, uc, "left")
        nkeys = np.searchsorted(pc, uc, "right") - lo
        sole = pk[lo]                                     # the key, where a cell holds one
        uc, nkeys, sole = uc[ok], nkeys[ok], sole[ok]
        if not len(uc):
            return out
        oc = c[mine]
        at = np.minimum(np.searchsorted(uc, oc), len(uc) - 1)
        at = at[uc[at] == oc]
        if not len(at):
            return out
        multi = nkeys[at] >= 2
        good = multi.copy()
        if multi.mean() < 0.5:                            # only then does a lone key's word count
            without = np.asarray(ref, bool).copy()
            without[i] = False
            alone = []
            for j in np.unique(sole[at][~multi]).tolist():
                off = self.offset(self.names[j], without)[0]
                if np.isfinite(off) and abs(off) <= limit:
                    good |= ~multi & (sole[at] == j)
                else:
                    alone.append(str(self.names[j]))
            out["alone"] = alone
        out.update(judged=bool(good.mean() >= 0.5), share=float(good.mean()), samples=int(len(at)))
        return out


def day_offsets(pair_cell, pair_vals, pair_key, min_frames, min_pairs=30):
    """
    For each (camera, local day): the median difference between that
    day's samples and the DEM built WITHOUT that day, cell by cell.
    Leaving the day out matters: a day with many frames dominates the
    cells it covers and would otherwise hide its own error.
    One pass, every key against all the others: the --max-frame-offset test
    (keys = frames), and the first pass of the day test (day_consistency).
    -> {key: (offset or NaN, comparable samples)}, keys sorted.
    """
    return LeaveOneOut(pair_cell, pair_vals, pair_key, min_frames, min_pairs).offsets()


# --max-day-offset, the day test. The days the cameras saw are compared with one another: a day
# whose waterlines sit, on median, more than the limit from the DEM of the other days is left
# out. ONE pass (each day against the DEM of all the others) lets the days it is about to leave
# out shape that DEM: on the station's whole-archive DEM of 13 Oct 2026 c1's pre-storm days
# (+0.27 to +1.40 m) pulled c1's cells up and four post-storm c1 days were left out as LOW
# (-0.15 to -0.26 m). The test is therefore iterative (day_consistency): the day furthest over the
# limit goes first and the others are measured again without it.
# Fewer than DAY_TEST_MIN_DAYS testable camera-days and the DEM without one day is the other day
# alone: the test cannot tell which one is off. In the first pass it then leaves none out
# (survey_products.py, DEM_DAY_TEST_MIN, does not ask for the test then); later on, the days it
# cannot judge keep the first pass's verdict and a WARNING names them (day_consistency, step 2).
DAY_TEST_MIN_DAYS = 3
# --anchor-date (the survey-date builds, survey_products.py): each camera's testable day nearest the
# survey date, within this many days of it, is never left out. Without it, on a camera whose days
# drift through the week (real change, or a bias that follows the tides), the test keeps whichever
# end of the drift weighs more: c2 of the 23 Oct 2024 build (+0.33, +0.26, +0.11, +0.09, -0.13,
# -0.26 m over 20-26 Oct) kept 24 and 26 Oct only and left out the survey's own day. With the
# anchor the days kept are those that agree with the survey's.
DAY_ANCHOR_REACH = 1
# A camera whose day offsets (final, below) run one way with time gets a WARNING: real beach change
# through the window, or an error that drifts with time (pointing, lens, a detection bias that
# follows the tides or waves), gives that pattern, and the day test then trims whichever end lies
# further out -- the station's 7-day DEM of 5-11 Oct 2026: c1 -0.238, -0.058, -0.048, -0.034,
# +0.016, +0.132 m, the 5th left out and the 11th kept, by one pass and by this test alike.
# Iterating does not cure that: once one end is left out the DEM moves towards the other days, so
# on a steady trend the same end tends to go again. The warning is the remedy; the synthetic Mar
# 2025 fixture (static beach, c2 drifting +0.07 to -0.09 m over 3-9 Mar with the detector) gets it
# too. The rule (day_trends): >= DAY_TREND_MIN_DAYS testable days, max - min offset above the
# limit, |Spearman rho| of offset against date above DAY_TREND_MIN_RHO (the rule asked for) and,
# up to 8 days, so far from 0 that days in random order reach it, in that direction, at most
# DAY_TREND_ALPHA of the time (exact, from every order of the days). That asks rho 1.00 over 5 days
# (0.90, one pair swapped, comes 4 % of the time by chance) and 0.89 over 6; from 7 days on it is
# |rho| > 0.8 itself (0.82 over 7 days, 0.81 over 8). And a lone day the test left out at either
# end of a camera's days cannot make the trend by itself (one bad day at an end ranks first or
# last and lifts rho): the camera's other days must still run the same way, |rho| >
# DAY_TREND_REST_RHO in that direction, whatever their number -- not the chance part, which the 4
# days left of 5 can never pass (asking for it, before 14 Oct 2026, hid every 5-day trend trimmed at
# one end: the station's 7-day c1 seen on one day fewer gave -0.241 ... +0.111 m, strictly rising,
# 5 Oct left out, and no warning). When the other days would not warn on their own the WARNING says
# the trend rests on the day left out (the end of the trend, or a bad day that lies beyond it).
# Simulated full-pipeline windows (limit 0.15 m, fixed seeds), the share that warns -- 7-day
# windows, two cameras (1000 of each kind): no trend (day biases sd 0.04-0.08 m) 2-6 %; one bad day
# (0.4 m) at an end of one camera's days 3-6 % (7-10 % with no end-day rule, 2-5 % with the chance
# part asked of the other days, as before); a trend of 0.03-0.06 m/day on one camera (day noise
# 0.03-0.08 m) 67-87 % (65-86 % before). A camera seen on 5 of the 7 days (c2 steady; 400-1000 of
# each kind): a trend of 0.04-0.08 m/day (an end trimmed in 61-99 % of them) 31-68 % (15-61 %
# before), one bad end day 4-6 % (0-1 %), no trend 1 %; on 6 days: trend 62-94 % (45-91 %), bad
# end day 6 % (0-4 %), no trend 1 %. DAY_TREND_ALPHA = 0.01 (rho 0.89 over 7 days) flagged only
# 46-70 % of the 7-day trends and missed a 0.35 m trend planted on the station's own data (rho
# +0.86 over 7 days, the two days at its low end trimmed); 0.05 is the rule asked for itself from 5
# days on (|rho| > 0.8: 0.90 over 5 days, 0.83 over 6), which flags 51-91 % of those 5-day trends
# and 70-97 % of the 6-day ones, at 6-7 % with a bad end day and 2-5 % with no trend. Each WARNING
# turns the station status to WARN.
DAY_TREND_MIN_DAYS = 5
DAY_TREND_MIN_RHO = 0.8
DAY_TREND_ALPHA = 0.02
DAY_TREND_REST_RHO = 0.8
DAY_TEST_METHOD = ("iterative: the camera-day furthest over the limit is left out and the others' "
                   "offsets are recomputed without it, until none is over; a day is left out only "
                   "when the test can judge it (where it is compared its DEM holds two or more other "
                   "kept days, or one that the others confirm), else it keeps the first pass's verdict "
                   "(a day kept so is tried again in later rounds); "
                   "each left-out day is then tried again against the DEM of the kept days and taken "
                   "back, one at a time (the closest first), when it is within the limit and puts no "
                   "kept day over it")


def ranks(x):
    """Ranks 1..n of x (ties share their mean rank), as Spearman's rho needs them. No scipy."""
    x = np.asarray(x, float)
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    r = np.empty(len(x))
    if not len(x):
        return r
    starts = np.flatnonzero(np.r_[True, xs[1:] != xs[:-1]])
    ends = np.r_[starts[1:], len(x)]
    for a, b in zip(starts, ends):
        r[order[a:b]] = (a + 1 + b) / 2.0
    return r


def spearman(x, y):
    """Spearman's rank correlation: Pearson's r of the ranks (ties averaged). NaN when either
    side is constant or there are fewer than 2 values."""
    if len(x) < 2:
        return np.nan
    rx, ry = ranks(x), ranks(y)
    rx, ry = rx - rx.mean(), ry - ry.mean()
    den = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / den) if den > 0 else np.nan


_SPEARMAN_PERMS = {}
_SPEARMAN_NULL = {}


def spearman_null(n, y=None):
    """Spearman's rho of n values against time in each of their n! orders (n <= 8: 8! = 40320), as
    spearman() computes it (Pearson's r of the ranks, ties averaged). y: the values, whose ties the
    orders keep (the exact null for that tie pattern); None: n distinct values."""
    r = np.arange(1.0, n + 1) if y is None else ranks(y)
    key = (n, tuple(np.sort(r).tolist()))
    if key not in _SPEARMAN_NULL:
        if n not in _SPEARMAN_PERMS:
            from itertools import permutations
            _SPEARMAN_PERMS[n] = np.array(list(permutations(range(n))), dtype=np.int64)
        rc = r - r.mean()
        t = np.arange(n) - (n - 1) / 2.0
        den = np.sqrt((rc ** 2).sum() * (t ** 2).sum())
        _SPEARMAN_NULL[key] = (rc[_SPEARMAN_PERMS[n]].dot(t) / den if den > 0
                               else np.zeros(len(_SPEARMAN_PERMS[n])))
    return _SPEARMAN_NULL[key]


def spearman_tail(rho, n, y=None):
    """The chance that n values in random order give Spearman's rho >= `rho` against time: exact,
    from all n! orders, for n <= 8; None for larger n. y: the values themselves -- with ties among
    them the orders keep the tie pattern, so the tail is the exact one for those values (5 days
    with one tied pair cannot reach rho 1; their largest rho, 0.97, comes in 2 orders of 120, the
    tied pair either way, where the table for distinct values said 1 in 120); without y, n distinct
    values. The null is symmetric (reversing an order turns rho into -rho), so the same chance holds
    for rho <= -`rho`: the chance is ONE-SIDED (in that direction); either way it is twice this."""
    if n > 8:
        return None
    return float((spearman_null(n, y) >= rho - 1e-9).mean())


def trend_rho_needed(n, alpha=DAY_TREND_ALPHA, min_rho=DAY_TREND_MIN_RHO):
    """What day_trends asks of |rho| over n days: up to 8 days the smallest rho that n days in
    random order reach at most `alpha` of the time and above min_rho (flagged when |rho| >= it;
    None: n days cannot reach it); from 9 days on min_rho (flagged when |rho| > it). For days
    without ties; with ties the chance is computed for their pattern (spearman_tail)."""
    if n > 8:
        return min_rho
    ok = [float(v) for v in np.unique(np.round(spearman_null(n), 9))
          if v > min_rho and spearman_tail(v, n) <= alpha]
    return ok[0] if ok else None


def day_trend_rule(alpha=DAY_TREND_ALPHA, min_rho=DAY_TREND_MIN_RHO, rest_rho=DAY_TREND_REST_RHO):
    """The trend rule in words, with the |rho| it asks for by number of days."""
    need = ", ".join(f"{n}: {trend_rho_needed(n, alpha, min_rho):.2f}" if trend_rho_needed(n, alpha, min_rho)
                     else f"{n}: never" for n in range(DAY_TREND_MIN_DAYS, 9))
    return (f">= {DAY_TREND_MIN_DAYS} testable days, range > the limit, |Spearman rho| > {min_rho:g} and (up "
            f"to 8 days) reached by days in random order at most {100 * alpha:g} % of the time in that "
            f"direction (|rho| needed by days without ties {need}; from 9, > {min_rho:g}); when the one day "
            f"left out of a camera's days is its first or last, the other days must run the same way, "
            f"|rho| > {rest_rho:g} (any number of days), and the WARNING says when the trend rests on that day")


def day_trends(offsets, limit, min_days=DAY_TREND_MIN_DAYS, min_rho=DAY_TREND_MIN_RHO,
               alpha=DAY_TREND_ALPHA, left_out=(), rest_rho=DAY_TREND_REST_RHO):
    """
    Cameras whose camera-day offsets run one way with time. offsets: {"<camera> <YYYY-MM-DD>":
    offset (NaN = not testable)}; left_out: the days the test left out. A camera is listed when it
    has >= min_days testable days, max - min offset > limit, |Spearman rho| of offset against date
    > min_rho, and (up to 8 days) days in random order reach that |rho| in that direction at most
    `alpha` of the time (exact, for the days' own ties: spearman_tail).
    When exactly one of its testable days was left out and that day is its first or last, the
    trend may rest on that day (one bad day at an end ranks first or last and lifts rho): the
    camera's other days must still run the same way, |rho| > rest_rho in that direction -- the rho
    the rule asks for, without the chance part and whatever their number (the 4 days left of 5
    can never reach a 2 % chance, so asking for it hid every 5-day trend trimmed at one end). They
    are given as without, with whether they pass the whole rule on their own (alone); when they
    do not, the trend rests on that day as much as on them, and the WARNING says so.
    -> [{camera, days, spearman_rho, chance (that one-sided probability; None from 9 days),
    range_m, first: (day, offset), last: (day, offset), without: None or {day, days, spearman_rho,
    chance (None from 9 days), range_m, alone}}], cameras sorted.
    """
    def passes(y, sign=None):
        rho = spearman(np.arange(len(y)), y)
        if not (np.isfinite(rho) and abs(rho) > min_rho) or (sign is not None and rho * sign <= 0):
            return None
        chance = spearman_tail(abs(rho), len(y), y)
        return None if chance is not None and chance > alpha else (rho, chance)

    by_cam = {}
    for k, off in offsets.items():
        if off is not None and np.isfinite(off):
            cam, day = k.split(" ", 1)
            by_cam.setdefault(cam, []).append((day, float(off)))
    gone = {}
    for k in left_out:
        cam, day = k.split(" ", 1)
        gone.setdefault(cam, set()).add(day)
    out = []
    for cam in sorted(by_cam):
        days = sorted(by_cam[cam])
        if len(days) < min_days:
            continue
        y = np.array([o for _, o in days])
        span = float(y.max() - y.min())
        hit = passes(y) if span > limit else None
        if hit is None:
            continue
        rho, chance = hit
        without = None
        out_days = [i for i, (d, _) in enumerate(days) if d in gone.get(cam, ())]
        if len(out_days) == 1 and out_days[0] in (0, len(days) - 1):
            i = out_days[0]
            rest = np.delete(y, i)
            rho2 = spearman(np.arange(len(rest)), rest)
            if not (np.isfinite(rho2) and rho2 * np.sign(rho) > rest_rho):
                continue                    # the other days do not run that way: one bad end day
            alone = len(rest) >= min_days and float(rest.max() - rest.min()) > limit and \
                passes(rest, np.sign(rho)) is not None
            without = {"day": f"{cam} {days[i][0]}", "days": len(rest), "spearman_rho": round(rho2, 3),
                       "chance": spearman_tail(abs(rho2), len(rest), rest),
                       "range_m": round(float(rest.max() - rest.min()), 3), "alone": bool(alone)}
        out.append({"camera": cam, "days": len(days), "spearman_rho": round(rho, 3),
                    "chance": chance, "range_m": round(span, 4), "first": days[0], "last": days[-1],
                    "min": days[int(np.argmin(y))], "max": days[int(np.argmax(y))], "without": without})
    return out


def anchor_days(first, anchor, reach=DAY_ANCHOR_REACH):
    """The survey-date anchor of the day test (--anchor-date): (kept, reference). reference:
    every camera's testable days within `reach` days of the anchor date -- the DEM the other days
    are measured against first; kept: for each camera, its reference day(s) nearest the date
    (ties: every day at that distance), never left out. first: {"<camera> <day>": (offset,
    samples)}."""
    from datetime import date
    a = date.fromisoformat(anchor)
    best, ref = {}, []
    for k, (off, _) in first.items():
        if not np.isfinite(off):
            continue
        cam, day = k.split(" ", 1)
        dist = abs((date.fromisoformat(day) - a).days)
        if dist > reach:
            continue
        ref.append(k)
        if cam not in best or dist < best[cam][0]:
            best[cam] = (dist, [k])
        elif dist == best[cam][0]:
            best[cam][1].append(k)
    return sorted(k for _, ks in best.values() for k in ks), sorted(ref)


def day_consistency(pair_cell, pair_vals, pair_key, limit, min_frames, min_pairs=30,
                    min_days=DAY_TEST_MIN_DAYS, anchor=None):
    """
    The --max-day-offset test. pair_key names each sample's camera-day ("<camera> <day>").
    anchor: a survey date (--anchor-date); each camera's testable day nearest it, within
    DAY_ANCHOR_REACH days (anchor_days), is never left out, so the days the test keeps are
    those that agree with the survey's own days.

    1. First pass: every camera-day against the DEM of all the other days (day_offsets). With
       fewer than `min_days` testable days (offset not NaN) every day is kept (too_few).
    2. Drop: while kept, testable days are more than `limit` off the DEM of the other kept days,
       the one furthest off (largest |offset|; a tie goes to the first key in sorted order) that
       the test can JUDGE is left out, and the kept days are measured again without it. A day can
       be judged (LeaveOneOut.judge) while `min_days` or more testable days are kept and when, in
       at least half the cells where it is compared, its DEM holds two or more other kept days,
       or one other day that the rest confirm (its offset from them, without the day, is within
       the limit). Otherwise the test cannot tell which day is off: two days that only see each
       other, a day whose cells only a day as far off reaches, the last two days. Days over the
       limit that cannot be judged keep the FIRST PASS's verdict -- the one-pass test's outcome:
       left out (all at once) if over the limit there, kept if within it (not_judged; the command
       line prints a WARNING) -- and the loop goes on with the other days. A day kept so is tried
       again in later rounds: once the test can judge it, it is left out like any other (and is
       no longer listed as not judged).
    3. Take back: every left-out day is measured against the DEM of the kept days. Of those
       within the limit, the closest (smallest |offset|; ties by key) is taken back if that puts
       no kept day over the limit that was within it (else the next closest is tried); then all
       are measured again against the new set, until a round takes none back. One day at a time
       mirrors step 2 and makes the result independent of the order the days are tried in; the
       kept set only grows here, so this ends.
    4. Final offsets: kept days against the DEM of the other kept days; left-out days against
       the DEM of the kept days. Every kept, testable day is then within the limit but those
       step 2 kept unjudged (over_kept); a day testable in the first pass that is not once the
       days left out are gone (too few comparable samples left: the days it was compared with left
       out, or those kept too thin on their own, e.g. two days of 2 frames a cell once a third is
       left out, min_frames 3) is kept untested (untested; a not_judged day kept too). A day kept
       untested keeps its place whatever its first-pass offset: that offset was measured against
       days since left out (the one-pass test would have rejected it when over the limit, and kept
       it when within, as here).

    A not_judged entry records the step where the day was set aside; where it ended up (left
    out, taken back, kept over or within the limit, untested) is in rejected, taken_back,
    over_kept and untested, and day_test_fate() says it. A day left out unjudged that the final
    DEM can judge (over the limit of the kept days, and judge() against them says so) is marked
    judged_in_end: its verdict no longer rests on the first pass alone.
    With the days split into two groups of about equal weight (say three at 0 m and three at
    +0.6 m on the same cells), the first day left out decides which group stays -- whichever the
    noise puts further out -- where one pass rejected both; nothing warns of that: read a day
    left out at little over the limit, or a camera's days left out together, with care.

    Day offsets never depend on the order of the samples (medians), and every choice above is
    made by |offset| then key, so the result is deterministic.
    -> dict: first, final ({key: (offset or NaN, samples)}, keys sorted), steps ([(key, offset,
    samples)] in the order left out), not_judged ([{day, offset (then), first, left_out (bool:
    the first pass's verdict, applied then), why, judged_in_end (bool)}]), judged_later (days set
    aside, kept, then judged and left out at a later step), retests ([[(key, offset, samples)] per
    round, against the kept days]), refused ([(key, offset, [kept days it would put over])]), taken_back
    (keys in order), rejected (sorted keys), untested, over_kept (sorted keys), too_few (bool),
    testable (first-pass count).
    """
    loo = LeaveOneOut(pair_cell, pair_vals, pair_key, min_frames, min_pairs)
    names = loo.names.tolist()
    idx = loo.index
    first = loo.offsets()
    testable = [k for k in names if np.isfinite(first[k][0])]
    anchored, anchor_ref = anchor_days(first, anchor) if anchor else ([], [])
    anchored = set(anchored)
    out = {"first": first, "final": dict(first), "steps": [], "not_judged": [], "judged_later": [], "retests": [],
           "refused": [], "taken_back": [], "rejected": [], "untested": [], "over_kept": [],
           "too_few": len(testable) < min_days, "testable": len(testable),
           "anchor": anchor, "anchored": sorted(anchored), "anchor_reference": anchor_ref}
    if out["too_few"]:
        return out

    def over(off):
        return bool(np.isfinite(off) and abs(off) > limit)

    def rank(k, off):                       # furthest first; a tie goes to the first key
        return (-abs(off), idx[k])

    kept = np.ones(len(names), bool)
    out["vs_anchor"], out["anchor_left"] = {}, []
    if anchored:
        # The days around the survey date (within DAY_ANCHOR_REACH) are the reference: every other
        # day against their DEM, and over the limit, left out (furthest first) before the iterative
        # test, which then goes on with the rest (the reference days among them, but the one(s)
        # nearest the date). A day that shares too few cells with them is left to the iterative test.
        A = np.zeros(len(names), bool)
        A[[idx[k] for k in anchor_ref]] = True
        vs = loo.offsets(A, ~A)
        out["vs_anchor"] = vs
        for k in sorted((k for k in vs if over(vs[k][0])), key=lambda k: rank(k, vs[k][0])):
            out["steps"].append((k,) + tuple(vs[k]))
            out["anchor_left"].append(k)
            kept[idx[k]] = False
    cur = loo.offsets(kept, kept) if out["anchor_left"] else dict(first)
    settled = set()                         # over the limit, not judged, kept by the first pass
    while True:
        test = [k for k in names if kept[idx[k]] and np.isfinite(cur[k][0])]
        bad = sorted((k for k in test if over(cur[k][0]) and k not in anchored),
                     key=lambda k: rank(k, cur[k][0]))
        if not bad:
            break
        worst, why = None, {}
        if len(test) >= min_days:
            for k in bad:                   # judged lazily, furthest first (a day set aside too)
                j = loo.judge(k, kept, limit)
                if j["judged"]:
                    worst = k
                    break
                why[k] = j
        if worst is not None:
            out["steps"].append((worst,) + tuple(cur[worst]))
            kept[idx[worst]] = False
            if worst in settled:            # set aside unjudged earlier, judged now
                settled.discard(worst)
                out["not_judged"] = [r for r in out["not_judged"] if r["day"] != worst]
                out["judged_later"].append(worst)
        else:
            new = [k for k in bad if k not in settled]
            if not new:                     # every day over the limit already set aside
                break
            for k in sorted(new, key=lambda k: idx[k]):
                j = why.get(k)
                if j is None:
                    reason = (f"only {len(test)} testable camera-day(s) kept: its DEM is the other "
                              f"day(s) alone")
                elif j["alone"]:
                    reason = (f"where it is compared its DEM is, in most cells, {' / '.join(j['alone'])} "
                              f"alone, which the other kept days do not confirm")
                else:
                    reason = "where it is compared its DEM holds no day the others confirm"
                left = over(first[k][0])
                out["not_judged"].append({"day": k, "offset": cur[k][0], "first": first[k][0],
                                          "left_out": left, "why": reason})
                if left:
                    out["steps"].append((k,) + tuple(cur[k]))
                    kept[idx[k]] = False
                else:
                    settled.add(k)
            if not any(over(first[k][0]) for k in new):
                break                       # none left out: the next round would see the same
        cur = loo.offsets(kept, kept)

    dropped = [s[0] for s in out["steps"]]
    away = {}                               # left-out days: offset from the DEM of the kept days
    while dropped:
        m = np.zeros(len(names), bool)
        m[[idx[k] for k in dropped]] = True
        res = loo.offsets(kept, m)
        rnd = [(k,) + tuple(res[k]) for k in dropped]
        out["retests"].append(rnd)
        for k, off, n in rnd:
            away[k] = (off, n)
        cands = sorted((r for r in rnd if np.isfinite(r[1]) and abs(r[1]) <= limit),
                       key=lambda r: (abs(r[1]), idx[r[0]]))
        back, refused = None, []
        for k, off, n in cands:
            trial = kept.copy()
            trial[idx[k]] = True
            new = loo.offsets(trial, trial)
            pushed = [j for j in cur if over(new[j][0]) and not over(cur[j][0])]
            if pushed:
                refused.append((k, off, pushed))
                continue
            back, kept, cur = k, trial, new
            break
        if back is None:
            out["refused"] = refused        # the last round's: those stand
            break
        out["taken_back"].append(back)
        dropped.remove(back)
        away.pop(back)
    out["rejected"] = sorted(dropped)
    if out["steps"] or settled:
        final = dict(cur)
        final.update({k: away[k] for k in dropped})
        out["final"] = {k: final[k] for k in names}
    fin = out["final"]
    out["untested"] = [k for k in names if kept[idx[k]] and np.isfinite(first[k][0])
                       and not np.isfinite(fin[k][0])]
    out["over_kept"] = [k for k in names if kept[idx[k]] and over(fin[k][0])]
    # a day left out unjudged that the final DEM can judge -- over the limit of the kept days, and
    # where it is compared they are two or more days or one the others confirm: judged in the end
    ntest = sum(1 for k in names if kept[idx[k]] and np.isfinite(fin[k][0]))
    for r in out["not_judged"]:
        k = r["day"]
        r["judged_in_end"] = bool(k in out["rejected"] and over(fin[k][0]) and ntest + 1 >= min_days
                                  and loo.judge(k, kept, limit)["judged"])
    return out


def day_test_fate(dt, day):
    """Where a camera-day of a day_consistency result ended up: 'left out', 'taken back', 'kept
    over the limit', 'kept untested' or 'kept within the limit' (or not testable)."""
    if day in dt["rejected"]:
        return "left out"
    if day in dt["taken_back"]:
        return "taken back"
    if day in dt["over_kept"]:
        return "kept over the limit"
    if day in dt["untested"]:
        return "kept untested"
    return "kept within the limit" if np.isfinite(dt["final"][day][0]) else "kept, not testable"


def _off(x):
    return f"{x:+.3f} m" if np.isfinite(x) else "   --   "


def day_test_report(dt, limit, n_left_points=None):
    """
    The lines the day test prints, the first block in the format survey_products.py
    (dem_failure_reason) reads: 'Day consistency (...)', then one line per camera-day starting
    with its key and its first-pass offset, '<-- REJECTED' on the days left out, then 'left out
    N camera-day(s), M point(s)'. When the test left any day out on the way (or could not judge
    one), each line also gives the final offset, and the lines after the block say in what order
    and why. The WARNING for days it could not judge is day_test_warning's, and these lines never
    hold that word (the station status counts the lines that do: one per warning).
    n_left_points: points of the rejected days (None: no 'left out' line).
    """
    first, final = dt["first"], dt["final"]
    moved = bool(dt["steps"] or dt["not_judged"])
    order = {s[0]: i + 1 for i, s in enumerate(dt["steps"])}
    nj = {r["day"]: r for r in dt["not_judged"]}
    back = set(dt["taken_back"])
    rej = set(dt["rejected"])
    refused = {r[0] for r in dt["refused"]}
    untested = set(dt["untested"])
    over_kept = set(dt["over_kept"])
    later = set(dt.get("judged_later", ()))
    lines = [f"Day consistency (median offset from the DEM without that day, limit +/-{limit} m):"]
    for k in sorted(first):
        off, n = first[k]
        fo, fn = final[k]
        if dt["too_few"]:
            if not np.isfinite(off):
                note = "  (too few comparable samples -- kept)"
            elif abs(off) > limit:
                note = f"  (off by more than the limit -- kept: fewer than {DAY_TEST_MIN_DAYS} testable days)"
            else:
                note = ""
        elif k in rej and k in dt.get("anchor_left", ()):
            note = (f"  <-- REJECTED (off the days around the survey date by {_off(dt['vs_anchor'][k][0]).strip()}, "
                    f"left out {ordinal(order[k])})")
        elif k in dt.get("anchored", ()):
            note = ("  (nearest the survey date: kept whatever its offset" +
                    (", over the limit of the kept days" if np.isfinite(fo) and abs(fo) > limit else "") + ")")
        elif k in rej and k in nj and nj[k].get("judged_in_end"):
            note = ("  <-- REJECTED (left out unjudged, as the first pass had it; the kept days judge it "
                    "over the limit in the end)")
        elif k in rej and k in nj:
            note = "  <-- REJECTED (not judged: over the limit in the first pass)"
        elif k in rej and k in refused:
            note = (f"  <-- REJECTED (left out {ordinal(order[k])}; within the limit of the kept days, "
                    f"not taken back: see below)")
        elif k in rej:
            note = f"  <-- REJECTED (left out {ordinal(order[k])})"
        elif k in back and k in nj:
            note = ("  (not judged: over the limit in the first pass, left out; taken back: within the "
                    "limit of the kept days)")
        elif k in back:
            note = f"  (left out {ordinal(order[k])}, taken back: within the limit of the kept days)"
        elif k in nj and k in over_kept:
            note = ("  (not judged: within the limit in the first pass, kept -- over the limit of the "
                    "kept days; see the warning below)")
        elif k in nj and k in untested:
            note = ("  (not judged: within the limit in the first pass, kept -- too few comparable "
                    "samples left, untested; see the warning below)")
        elif k in nj:
            note = ("  (not judged: within the limit in the first pass, kept -- within the limit of "
                    "the kept days in the end)")
        elif not np.isfinite(off):
            note = "  (too few comparable samples -- kept)"
        elif k in untested:
            note = "  (too few comparable samples left -- kept untested; see the warning below)"
        elif abs(fo) > limit:
            note = "  (over the limit -- kept)"
        else:
            note = ""
        fin = f"   final {_off(fo)}" if moved else ""
        lines.append(f"   {k}   {_off(off)}   ({n} samples){fin}{note}")
    if dt["too_few"]:
        lines.append(f"   only {dt['testable']} testable camera-day(s): the test cannot tell which day is off; "
                     f"every day kept")
        return lines
    if n_left_points is not None and rej:
        lines.append(f"   left out {len(rej)} camera-day(s), {n_left_points} point(s)")
    if dt.get("anchor"):
        lines += _wrap(f"anchored on the survey date {dt['anchor']} (--anchor-date): "
                       + (f"{', '.join(dt['anchored'])} (each camera's testable day nearest it) never left "
                          f"out; every day further than {DAY_ANCHOR_REACH} day from it measured first against "
                          f"the DEM of the days within {DAY_ANCHOR_REACH} day ({', '.join(dt['anchor_reference'])}), "
                          f"and left out over the limit, then the test below on the rest"
                          + (": " + ", ".join(_unit(f"{k} ({_off(dt['vs_anchor'][k][0]).strip()})")
                                              for k in sorted(dt["vs_anchor"])
                                              if np.isfinite(dt["vs_anchor"][k][0]))
                             if any(np.isfinite(v[0]) for v in (dt.get("vs_anchor") or {}).values()) else
                             ": no other day comparable with them")
                          if dt["anchored"] else
                          f"no camera has a testable day within {DAY_ANCHOR_REACH} day of it: nothing "
                          f"anchored"))
    if not moved:
        return lines
    lines += _wrap(f"({DAY_TEST_METHOD}. First column: each day against the DEM of every other day, "
                   f"the first pass; final: a kept day against the DEM of the other kept days, a day "
                   f"left out against the DEM of the kept days)")
    if dt["steps"]:
        lines += _wrap("order left out: " + ", ".join(
            _unit(f"{i + 1}. {k} ({_off(o)} at that step{', not judged' if k in nj else ''}"
                  f"{', set aside unjudged at an earlier step' if k in later else ''})")
            for i, (k, o, n) in enumerate(dt["steps"])))
    for r in dt["not_judged"]:
        k = r["day"]
        fate = day_test_fate(dt, k)
        end = {"left out": (f"; in the end over the limit of the kept days ({_off(final[k][0])}), which "
                            f"can judge it" if r.get("judged_in_end") else ""),
               "taken back": f", then taken back: within the limit of the kept days ({_off(final[k][0])})",
               "kept over the limit": f"; in the end still over the limit of the kept days ({_off(final[k][0])})",
               "kept untested": "; in the end untested (too few comparable samples left)",
               "kept within the limit": f"; in the end within the limit of the kept days ({_off(final[k][0])})",
               }.get(fate, "")
        lines += _wrap(f"not judged: {_unit(k)} ({_off(r['offset'])} at that step, first pass "
                       f"{_off(r['first'])}): {r['why']}; "
                       + ("left out, as the first pass had it" if r["left_out"] else
                          "kept, as the first pass had it") + end)
    for k, o, pushed in dt["refused"]:
        lines += _wrap(f"not taken back: {_unit(k)} ({_off(o)} from the kept days, within the limit) "
                       f"would put {', '.join(_unit(j) for j in pushed)} over the limit")
    if dt["taken_back"]:
        lines += _wrap("taken back: " + ", ".join(_unit(f"{k} ({_off(final[k][0])})") for k in dt["taken_back"]))
    elif dt["refused"]:
        lines += _wrap("taken back: none (the days within the limit of the kept days would each put a "
                       "kept day over it, above; every other day left out is still over the limit, or "
                       "not comparable, against the DEM of the kept days)")
    elif dt["steps"]:
        lines += _wrap("taken back: none (every day left out is still over the limit, or not comparable, "
                       "against the DEM of the kept days)")
    return lines


def day_test_unjudged(dt):
    """The days the WARNING names, as the final DEM has them: (left out without being judged --
    still left out; kept over the limit of the kept days without being judged; kept untested).
    A day set aside unjudged that was taken back, or that ended within the limit of the kept
    days, is in the DEM and agrees with it; one left out that the kept days can judge in the end
    (judged_in_end) is off them, judged: neither is named."""
    nj = {r["day"] for r in dt["not_judged"]}
    late = {r["day"] for r in dt["not_judged"] if r.get("judged_in_end")}
    gone = [k for k in sorted(nj) if k in dt["rejected"] and k not in late]
    stay = [k for k in dt["over_kept"] if k in nj]
    return gone, stay, list(dt["untested"])


def day_test_warning(dt, limit):
    """The WARNING (one line, the station status counts them; details on the lines after it) when
    the final DEM leaves out a day the test could not judge, keeps one over the limit that it
    could not judge, or keeps a day untested: '' and [] otherwise.
    -> (the WARNING line, [its lines])"""
    if dt["too_few"]:
        return "", []
    gone, stay, untested = day_test_unjudged(dt)
    parts = []
    if gone:
        parts.append(f"left out without being judged: {', '.join(gone)}")
    if stay:
        parts.append(f"kept over the limit without being judged: {', '.join(stay)}")
    if untested:
        parts.append(f"kept untested: {', '.join(untested)}")
    if not parts:
        return "", []
    head = "WARNING: the day test could not judge every camera-day -- " + "; ".join(parts) + "."
    why = []
    if gone or stay:
        why.append(f"Where a day's DEM is one other day that the rest do not confirm, or when fewer than "
                   f"{DAY_TEST_MIN_DAYS} testable days are left, the test cannot tell which day is off: "
                   f"those days keep the first pass's verdict (the one-pass test's: over the limit, "
                   f"left out; within it, kept)")
    if stay:
        why.append(f"Kept although over the limit of the kept days: {', '.join(_unit(k) for k in stay)}")
    if untested:
        why.append("Kept without a final check, as too few comparable samples are left once the days left "
                   "out are gone; in the first pass each was compared with every other day, the ones since "
                   "left out included: " + ", ".join(_unit(f"{k} ({_off(dt['first'][k][0])})") for k in untested))
    if stay or untested:
        why.append("Read the DEM with care where these days lie (the lines above give their offsets)")
    return head, [head] + _wrap(". ".join(why) + ".")


def _unit(text):
    """Text _wrap never breaks (a camera-day and its number stay on one line)."""
    return text.replace(" ", "\0")


def _wrap(text, width=100):
    """An explanatory line of the log, wrapped: indented 3, continued at 5."""
    import textwrap
    return [ln.replace("\0", " ") for ln in
            textwrap.wrap(text, width=width, initial_indent="   ", subsequent_indent="     ",
                          break_on_hyphens=False, break_long_words=False)]


def ordinal(n):
    return f"{n}{'th' if 10 <= n % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')}"


# What the method alone gives (review r1): on a STATIC synthetic beach with a fixed pointing the detection's own
# drift makes the same WARNING, so the warning alone shows neither a camera fault nor a change of the beach
TREND_METHOD_DRIFT = ("The method alone drifts so: on a static synthetic beach with a fixed pointing (the Mar 2025 "
                      "replica) c2's offsets give rho -0.96 and a range of 0.165 m over 7 days, the same warning; "
                      "only a drift clearly larger than that, or one the pointing check or a survey confirms, "
                      "points at the camera or the beach.")


def trend_chance_text(c):
    """The chance part of a trend WARNING: how often days in random order reach that rho IN THAT DIRECTION (the
    one-sided chance the rule uses), and either way (twice that: the null is symmetric). '' without one."""
    if not c:
        return ""
    return (f", which days in random order reach in that direction 1 time in {int(round(1 / c)):,} (either way "
            f"1 in {int(round(1 / (2 * c))):,})")


def trend_span_text(t, limit):
    """', the lowest ... the highest ..., a range of X m > the limit': the range over ALL the camera's days, with
    the days that make it (the first and last day need not), from the unrounded offsets where the record has them."""
    lo, hi = t.get("min"), t.get("max")
    if lo and hi:
        span = float(hi[1]) - float(lo[1])
        where = f"; over all its days {lo[1]:+.3f} m ({lo[0]}) to {hi[1]:+.3f} m ({hi[0]})"
    else:
        span, where = float(t["range_m"]), ""
    return f"{where}, a range of {span:.3f} m > the {float(limit):g} m limit" if limit is not None else \
        f"{where}, a range of {span:.3f} m"


def day_trend_warning(t, limit, rejected):
    """The WARNING lines for one camera's trend (day_trends), naming the days the test left out."""
    cam = t["camera"]
    gone = [k for k in rejected if k.split(" ", 1)[0] == cam]
    (d0, o0), (d1, o1) = t["first"], t["last"]

    chance = trend_chance_text
    w = t.get("without")
    rest, carried = "", ""
    if w and w.get("alone", True):
        rest = (f"; without {w['day']}, the one day left out, at an end, the other {w['days']} days give "
                f"rho {w['spearman_rho']:+.2f}{chance(w.get('chance'))}")
    elif w:
        short = (f"too few to pass the rule alone (it asks {DAY_TREND_MIN_DAYS})" if w["days"] < DAY_TREND_MIN_DAYS
                 else f"not enough to pass the rule alone (range {w['range_m']:.2f} m)" if w["range_m"] <= limit
                 else "not enough to pass the rule alone")
        rest = (f"; it rests on {w['day']}, the one day left out, at an end: the other {w['days']} days "
                f"give rho {w['spearman_rho']:+.2f}{chance(w.get('chance'))}, {short}")
        carried = (f" {_unit(w['day'])} may be the end of the trend or one bad day that happens to lie "
                   f"beyond it: the other days run the same way, but on their own they would not warn.")
    return [f"WARNING: {cam}'s camera-day offsets run one way with time over {t['days']} days "
            f"(Spearman rho {t['spearman_rho']:+.2f}{chance(t.get('chance'))}; {o0:+.3f} m on {d0} to "
            f"{o1:+.3f} m on {d1}{trend_span_text(t, limit)}{rest})."] + _wrap(
            "Real beach change through the window, or an error that drifts with time (the camera's "
            "pointing or lens, a detection bias that follows the tides or the waves), gives this pattern; "
            "the day test then trims the end(s) furthest out"
            + (f" -- left out here: {', '.join(_unit(k) for k in gone)}" if gone else " -- none left out here")
            + "." + carried + f" Check {cam}'s pointing and the beach change before reading this DEM as "
            "one surface. " + TREND_METHOD_DRIFT)


def day_test_record(dt, limit, trends):
    """What <stem>_info.json keeps of the day test (filters), so a page can name the days left out
    and say why. It describes the final DEM: each day the test could not judge with where it ended
    up, the WARNING as printed. Offsets in m (None = not testable), rounded to 0.1 mm."""
    def r(x):
        return round(float(x), 4) if x is not None and np.isfinite(x) else None
    warning, _ = day_test_warning(dt, limit)
    step = {s[0]: s[1] for s in dt["steps"]}
    order = {s[0]: i + 1 for i, s in enumerate(dt["steps"])}
    nj = {x["day"]: x for x in dt["not_judged"]}
    refused = {k: p for k, _, p in dt["refused"]}

    def rejected_row(k):
        fin = r(dt["final"][k][0])
        if k in nj and nj[k].get("judged_in_end"):
            how = "left out"
            why = (f"over the limit in the first pass ({_off(dt['first'][k][0]).strip()}) and left out "
                   f"unjudged; in the end off the DEM of the kept days by more than the limit, which can "
                   f"judge it")
        elif k in nj:
            how = "not judged"
            why = (f"over the limit in the first pass ({_off(dt['first'][k][0]).strip()} from the DEM of every "
                   f"other day); the test could not judge it, so the first pass's verdict stands"
                   + (f"; within the limit of the kept days now, but taking it back would put "
                      f"{', '.join(refused[k])} over the limit" if k in refused else ""))
        elif k in dt.get("anchor_left", ()):
            how = "off the days around the survey date"
            why = (f"off the DEM of the days around the survey date ({', '.join(dt['anchor_reference'])}) by "
                   f"{_off(dt['vs_anchor'][k][0]).strip()}, more than the limit (--anchor-date)")
        elif k in refused:
            how = "refused take-back"
            why = (f"left out {ordinal(order[k])}; within the limit of the kept days, but taking it back "
                   f"would put {', '.join(refused[k])} over the limit")
        else:
            how = "left out"
            why = (f"left out {ordinal(order[k])}: off the DEM of the other kept days by more than the limit")
        return {"day": k, "how": how, "offset_m": fin, "offset_at_step_m": r(step.get(k)),
                "first_pass_offset_m": r(dt["first"][k][0]), "why": why}
    rows = [rejected_row(k) for k in dt["rejected"]]
    return {
        # each day left out: its offset from the DEM of the kept days; where no kept day is
        # comparable with it, the offset it was left out at (over the limit then, from the DEM of
        # the days kept at that step). day_test.rejected gives both, the first pass's, and why.
        "camera_days_rejected_offsets_m": {x["day"]: x["offset_m"] if x["offset_m"] is not None
                                           else x["offset_at_step_m"] for x in rows},
        "day_test": {
            "method": DAY_TEST_METHOD if not dt["too_few"] else
            f"not applied: fewer than {DAY_TEST_MIN_DAYS} testable camera-days, every day kept",
            "limit_m": limit,
            # --anchor-date: the survey date, and the days kept whatever their offset (anchor_days)
            "anchor_date": dt.get("anchor"),
            "anchored_days": list(dt.get("anchored", [])),
            "anchor_reference_days": list(dt.get("anchor_reference", [])),
            "offsets_from_anchored_days_m": {k: r(v[0]) for k, v in (dt.get("vs_anchor") or {}).items()},
            "first_pass_offsets_m": {k: r(v[0]) for k, v in dt["first"].items()},
            "final_offsets_m": {k: r(v[0]) for k, v in dt["final"].items()},
            "samples": {k: int(v[1]) for k, v in dt["first"].items()},
            "left_out_order": [s[0] for s in dt["steps"]],
            # how: 'left out' (judged off the other kept days), 'not judged' (the first pass's
            # verdict) or 'refused take-back' (within the limit of the kept days, but taking it
            # back would put a kept day over): only 'left out' is 'off the rest by > the limit'
            "rejected": rows,
            "taken_back": list(dt["taken_back"]),
            # left out within the limit of the kept days: taking it back would put these over
            "refused_take_back": [{"day": k, "offset_m": r(o), "would_push": list(p)}
                                  for k, o, p in dt["refused"]],
            # over the limit at a step where the test could not tell which day is off: the first
            # pass's verdict was applied then (first_pass_verdict); verdict is where the day ended
            # up: 'left out', 'taken back' (within the limit of the kept days), 'kept over the
            # limit', 'kept untested' or 'kept within the limit'. The warning names those left out
            # (unless the kept days judge them in the end: judged_in_end), kept over the limit, or
            # kept untested
            "not_judged": [{"day": x["day"], "offset_m": r(x["offset"]), "first_pass_offset_m": r(x["first"]),
                            "first_pass_verdict": "left out" if x["left_out"] else "kept",
                            "verdict": day_test_fate(dt, x["day"]), "final_offset_m": r(dt["final"][x["day"]][0]),
                            "judged_in_end": bool(x.get("judged_in_end")), "why": x["why"]}
                           for x in dt["not_judged"]],
            "kept_untested": list(dt["untested"]),
            "kept_over_limit": list(dt["over_kept"]),
            "warning": warning or None,
            "trend_warnings": [{"camera": t["camera"], "days": t["days"], "spearman_rho": t["spearman_rho"],
                                "chance": None if t.get("chance") is None else round(t["chance"], 6),
                                "range_m": t["range_m"], "first_day": t["first"][0],
                                "first_offset_m": r(t["first"][1]), "last_day": t["last"][0],
                                "last_offset_m": r(t["last"][1]),
                                "min_day": (t.get("min") or (None, None))[0],
                                "min_offset_m": r((t.get("min") or (None, None))[1]),
                                "max_day": (t.get("max") or (None, None))[0],
                                "max_offset_m": r((t.get("max") or (None, None))[1]),
                                # the one day left out, at an end, and the camera's other days:
                                # alone, they pass the rule on their own; else the trend rests on it
                                "without": None if not t.get("without") else {
                                    "day": t["without"]["day"], "days": t["without"]["days"],
                                    "spearman_rho": t["without"]["spearman_rho"],
                                    "chance": None if t["without"].get("chance") is None
                                    else round(t["without"]["chance"], 6),
                                    "range_m": t["without"]["range_m"], "alone": t["without"]["alone"]}}
                               for t in trends],
            "trend_from": "final offsets",
            "trend_rule": day_trend_rule(),
        },
    }


def interpolate_between_contours(E, N, Z, frames, cell, max_spread, dem, max_edge):
    """
    Fill cells between waterlines by linear interpolation.

    The upper beach is crossed by few frames (fewer high tides are caught),
    so its cells rarely reach --min-points although the waterlines there
    are clean, nested contours. Every cell crossed by at least one frame
    (and not too noisy) becomes a node; the nodes are triangulated and each
    empty cell inside a triangle is interpolated from its corners -- the
    usual way an intertidal DEM is built from waterlines. Triangles with
    an edge longer than max_edge (m) are dropped, so nothing is bridged
    across a gap in the data or between separate stretches of beach.
    Cells already measured keep their value. Uses matplotlib.tri (no scipy).
    Returns (dem, source) with source 1 = measured, 2 = interpolated.
    """
    import matplotlib.tri as mtri
    node, _, _, e0, n0, ncols, nrows = build_grid(E, N, Z, cell, 1, max_spread, frames)
    source = np.where(np.isfinite(dem), 1, 0)
    ok = np.isfinite(node)
    if ok.sum() < 3:
        return dem, source
    rr, cc = np.nonzero(ok)
    x = e0 + (cc + 0.5) * cell
    y = n0 + (rr + 0.5) * cell
    tri = mtri.Triangulation(x, y)
    t = tri.triangles
    edge = np.max(np.stack([np.hypot(x[t[:, a]] - x[t[:, b]], y[t[:, a]] - y[t[:, b]])
                            for a, b in ((0, 1), (1, 2), (2, 0))]), axis=0)
    tri.set_mask(edge > max_edge)
    interp = mtri.LinearTriInterpolator(tri, node[ok])
    gr, gc = np.nonzero(~np.isfinite(dem))
    vals = np.ma.filled(interp(e0 + (gc + 0.5) * cell, n0 + (gr + 0.5) * cell), np.nan)
    out = dem.copy()
    good = np.isfinite(vals)
    out[gr[good], gc[good]] = vals[good]
    source[gr[good], gc[good]] = 2
    return out, source


def fine_dem(E, N, Z, frames, cell, along, across, min_frames, max_resid):
    """
    A DEM on a fine grid (e.g. 0.25 m) straight from the waterline points.

    WHY. Binning into 2 m cells and interpolating between cell centres
    gives a surface made of 2 m facets: the profile shows steps and kinks
    that are the grid, not the beach (C. Sherwood, Oct 2026, against the
    23 Jan 2025 lidar). Finer bins do not help -- most 0.25 m cells hold
    no waterline at all.

    HOW. At each node a plane z = a + b*du + c*dv is fitted to the nearby
    points by weighted least squares (local linear regression, as in
    Plant et al.'s argus bathymetry). The neighbourhood is an ellipse in
    beach coordinates: `along` metres alongshore, `across` metres
    cross-shore, tricube weights. Long alongshore because the beach varies
    slowly that way and each waterline is an alongshore curve; short
    cross-shore so the profile shape is kept. The beach orientation is the
    principal axis of the points. A plane, not a weighted mean, so the
    slope does not flatten the surface between contours.

    A node is filled only if its neighbourhood holds >= min_frames
    different frames AND points on both its landward and seaward side
    (no extrapolation beyond the highest or lowest waterline), and the
    fit's weighted RMS residual is <= max_resid. Points are first thinned
    to one per frame per half-cell, so a far-field waterline with many
    points per metre does not outweigh a near-field one.

    Returns (dem, resid_rms, n_frames, e0, n0) on a north-up grid.
    """
    # thin: one point (mean position) per frame per half-cell bin
    _, fid = np.unique(frames, return_inverse=True)
    h = cell / 2
    key = np.stack([fid, np.floor(E / h).astype(np.int64), np.floor(N / h).astype(np.int64)], 1)
    _, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.ravel()
    cnt = np.bincount(inv)
    E = np.bincount(inv, E) / cnt
    N = np.bincount(inv, N) / cnt
    Z = np.bincount(inv, Z) / cnt
    fid = np.bincount(inv, fid) / cnt
    fid = np.rint(fid).astype(np.int64)

    # beach axes: u alongshore (largest spread), v cross-shore
    ec, nc = E.mean(), N.mean()
    w_, vecs = np.linalg.eigh(np.cov(np.stack([E - ec, N - nc])))
    ua, va = vecs[:, 1], vecs[:, 0]
    u = (E - ec) * ua[0] + (N - nc) * ua[1]
    v = (E - ec) * va[0] + (N - nc) * va[1]

    e0 = np.floor(E.min() / cell) * cell
    n0 = np.floor(N.min() / cell) * cell
    ncols = int(np.ceil((E.max() - e0) / cell)) + 1
    nrows = int(np.ceil((N.max() - n0) / cell)) + 1
    gx = e0 + (np.arange(ncols) + 0.5) * cell
    gy = n0 + (np.arange(nrows) + 0.5) * cell
    GX, GY = np.meshgrid(gx, gy)
    gu = (GX - ec) * ua[0] + (GY - nc) * ua[1]
    gv = (GX - ec) * va[0] + (GY - nc) * va[1]

    # buckets of one window size: a node's points are all in its 3x3 buckets
    bu, bv = np.floor(u / along).astype(np.int64), np.floor(v / across).astype(np.int64)
    pts = {}
    order = np.lexsort((bv, bu))
    keys = np.stack([bu[order], bv[order]], 1)
    edges = np.flatnonzero(np.any(np.diff(keys, axis=0) != 0, axis=1)) + 1
    for a, b in zip(np.r_[0, edges], np.r_[edges, len(order)]):
        pts[(int(keys[a, 0]), int(keys[a, 1]))] = order[a:b]

    nbu = np.floor(gu / along).astype(np.int64).ravel()
    nbv = np.floor(gv / across).astype(np.int64).ravel()
    dem = np.full(gu.size, np.nan)
    resid = np.full(gu.size, np.nan)
    nfr = np.zeros(gu.size)
    nkeys = np.stack([nbu, nbv], 1)
    norder = np.lexsort((nbv, nbu))
    nk = nkeys[norder]
    nedges = np.flatnonzero(np.any(np.diff(nk, axis=0) != 0, axis=1)) + 1
    guf, gvf = gu.ravel(), gv.ravel()
    for a, b in zip(np.r_[0, nedges], np.r_[nedges, len(norder)]):
        ku, kv = int(nk[a, 0]), int(nk[a, 1])
        cand = [pts[(ku + i, kv + j)] for i in (-1, 0, 1) for j in (-1, 0, 1)
                if (ku + i, kv + j) in pts]
        if not cand:
            continue
        p = np.concatenate(cand)
        nodes = norder[a:b]
        du = u[p][None, :] - guf[nodes][:, None]
        dv = v[p][None, :] - gvf[nodes][:, None]
        r = np.sqrt((du / along) ** 2 + (dv / across) ** 2)
        w = np.where(r < 1, (1 - r ** 3) ** 3, 0.0)
        inside = w > 0
        # distinct frames in each neighbourhood
        lf, lfid = np.unique(fid[p], return_inverse=True)
        onehot = np.zeros((len(p), len(lf)))
        onehot[np.arange(len(p)), lfid] = 1
        nf = ((inside.astype(float) @ onehot) > 0).sum(1)
        both_sides = (np.where(inside, dv, np.inf).min(1) < 0) & \
                     (np.where(inside, dv, -np.inf).max(1) > 0)
        ok = (nf >= min_frames) & both_sides
        if not ok.any():
            continue
        w, du, dv = w[ok], du[ok], dv[ok]
        # weighted normal equations for z = a + b*du + c*dv, batched 3x3 solves
        X = np.stack([np.ones_like(du), du, dv], -1)                 # nodes x pts x 3
        XtW = X * w[..., None]
        A = np.einsum("npi,npj->nij", XtW, X) + np.eye(3) * 1e-9
        rhs = np.einsum("npi,p->ni", XtW, Z[p])
        try:
            coef = np.linalg.solve(A, rhs[..., None])[..., 0]
        except np.linalg.LinAlgError:
            continue
        fit = np.einsum("npi,ni->np", X, coef)
        rms = np.sqrt((w * (Z[p][None, :] - fit) ** 2).sum(1) / w.sum(1))
        sel = nodes[ok]
        dem[sel] = coef[:, 0]
        resid[sel] = rms
        nfr[sel] = nf[ok]
    if max_resid:
        dem[resid > max_resid] = np.nan
    shape = (nrows, ncols)
    return dem.reshape(shape), resid.reshape(shape), nfr.reshape(shape), e0, n0


def fill_small_gaps(dem, max_iterations=3):
    """
    Fills isolated nodata cells from their immediate neighbours.

    Deliberately limited: it closes pinholes inside measured areas, not
    voids between them. Each pass fills only cells with at least five
    of eight neighbours present, so a gap wider than a few cells stays
    open. Filling a real gap would manufacture terrain.
    """
    out = dem.copy()
    for _ in range(max_iterations):
        holes = np.isnan(out)
        if not holes.any():
            break
        padded = np.pad(out, 1, constant_values=np.nan)
        stack = np.stack([padded[a:a + out.shape[0], b:b + out.shape[1]]
                          for a in range(3) for b in range(3)
                          if not (a == 1 and b == 1)])
        neighbours = np.sum(np.isfinite(stack), axis=0)
        with np.errstate(invalid="ignore"):
            mean_nb = np.nanmean(stack, axis=0)
        fill = holes & (neighbours >= 5)
        if not fill.any():
            break
        out[fill] = mean_nb[fill]
    return out


def write_ascii_grid(path, grid, e0, n0, cell, nodata=-9999.0):
    """ESRI ASCII grid -- readable by QGIS, ArcGIS, GDAL, MATLAB."""
    flipped = np.flipud(grid)  # ASCII grids run north to south
    with open(path, "w") as f:
        f.write(f"ncols {grid.shape[1]}\n")
        f.write(f"nrows {grid.shape[0]}\n")
        f.write(f"xllcorner {e0:.3f}\n")
        f.write(f"yllcorner {n0:.3f}\n")
        f.write(f"cellsize {cell}\n")
        f.write(f"NODATA_value {nodata}\n")
        for r in flipped:
            f.write(" ".join("%.4f" % (v if np.isfinite(v) else nodata) for v in r) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("contour_csv")
    ap.add_argument("output_stem")
    ap.add_argument("--cell", type=float, default=2.0,
                    help="Grid cell size in metres (default 2.0). Below the georectification's "
                         "own 0.33 m accuracy there is nothing to gain; well above it the "
                         "beach profile gets smoothed away.")
    ap.add_argument("--min-points", type=int, default=3,
                    help="Minimum FRAMES (independent waterline crossings) for a cell to be filled "
                         "(default 3). One frame is a "
                         "single detection with no way to tell whether it was a good one.")
    ap.add_argument("--max-spread", type=float, default=0.5,
                    help="Blank cells whose 16-84 percentile elevation range exceeds this, in "
                         "metres (default 0.5). A cell where repeat crossings disagree by more "
                         "than half a metre is not measuring one surface.")
    ap.add_argument("--camera", default="both", choices=["c1", "c2", "both"])
    ap.add_argument("--trend-note", action="store_true",
                    help="With --max-day-offset: print the trend warning (a camera's offsets running one "
                         "way with time) as a NOTE line, not a WARNING. For the station's cron DEMs: "
                         "station_status.py counts WARNING lines, and on the whole-archive DEM (months of "
                         "beach pooled by design) or a 7-day DEM (a detection or setup bias that follows "
                         "the spring-neap tidal range gives the same pattern) the trend is expected, not a "
                         "station fault. The record in <stem>_info.json is the same; the day test's other "
                         "WARNING (a day it could not judge) stays a WARNING. Survey-date builds do not "
                         "pass it.")
    ap.add_argument("--anchor-date", default=None,
                    help="With --max-day-offset: a survey date (YYYY-MM-DD). Each camera's testable day "
                         "nearest it (within 1 day) is never left out, so the days kept are those that "
                         "agree with the survey's own days, not whichever end of a drift weighs more. "
                         "The survey-date builds (survey_products.py) pass it; the station's cron does not.")
    ap.add_argument("--no-track-filter", action="store_true",
                    help="Keep the runs of a line that cross the shore on the ground (by default "
                         "dropped: a line's track running more than 45 deg from the shore over 10 m is "
                         "the line tracer carrying on out over the water; waterline_track.py).")
    ap.add_argument("--max-day-offset", type=float, default=None,
                    help="Leave out camera-days whose samples sit, on median, more than this "
                         "many metres from the DEM built without that day (e.g. 0.15). "
                         "Catches days whose waterlines are consistently displaced -- a "
                         "different failure from the tide-direction check, which only sees "
                         "direction. Iterative: the day furthest over the limit is left out "
                         "and the others measured again without it, until none is over; a "
                         "left-out day within the limit of the kept days is then taken back, if "
                         "its return puts no kept day over the limit (else it stays out: a refused "
                         "take-back). "
                         "Days with fewer than 30 comparable samples are kept, and with fewer "
                         "than 3 testable camera-days every day is kept (the DEM without one "
                         "day is then the other day alone, so the test cannot tell which is "
                         "off); later on, days over the limit it cannot judge that way keep "
                         "the one-pass verdict, with a WARNING. Prints the first-pass and "
                         "final offsets, and a WARNING when a camera's offsets run one way "
                         "with time (real change, or an error drifting with time: the test "
                         "then trims an end): " + day_trend_rule().replace("%", "%%") + ".")
    ap.add_argument("--max-frame-offset", type=float, default=None,
                    help="Leave out single frames whose samples sit, on median, more than this "
                         "many metres from the DEM built without that frame (e.g. 0.2). Runs "
                         "before --max-day-offset, so one bad waterline does not cost its whole "
                         "day. Frames with fewer than 10 comparable cells are kept.")
    ap.add_argument("--fit-setup", action="store_true",
                    help="Report the wave-setup coefficient C (extract_elevation_contours.py "
                         "--setup-coef) that makes repeat crossings agree best, from the "
                         "water level and offshore Hs/Tp of each frame. Diagnostic only: the "
                         "DEM is built from whatever correction the contours already carry. "
                         "Run it on contours WITHOUT a correction applied. It offers a C to "
                         "apply only when C = 0 and the C in use both lie outside its near-best "
                         "band (the C whose spread is within 5 mm of its best) and the minimum "
                         "is not flat.")
    ap.add_argument("--setup-coef-in-use", type=float, default=None,
                    help="With --fit-setup: the C in use, compared with the fit's best. "
                         "Default: the C the contours carry (setup_correction_m / "
                         "sqrt(Hs*L0)), else SETUP_COEF of waterline_timex_cron.sh beside "
                         "this script.")
    ap.add_argument("--max-hs", type=float, default=None,
                    help="Leave out frames whose offshore wave height (offshore_hs_m, from "
                         "extract_elevation_contours.py --waves) exceeds this, in metres. In "
                         "big waves the waterline sits above still water by the wave setup, "
                         "so its GNSS-R elevation is too low for where it lies. Frames with no "
                         "wave record are kept. Off by default.")
    ap.add_argument("--start-date"); ap.add_argument("--end-date")
    ap.add_argument("--last-days", type=int, default=None,
                    help="Build the DEM from only the last N days of data, ending on the "
                         "latest date in the file (or --end-date). Pooling weeks of data "
                         "blurs real beach change into 'spread': on C2+C1 in Sep 2026 the "
                         "whole-archive spread was 0.278 m against 0.231 m for Sep 11-14 "
                         "alone, with more cells filled.")
    ap.add_argument("--series-dir", default=None,
                    help="Also save this DEM as <dir>/dem_<end date>_<N>d_{dem,spread,count}.asc "
                         "(with --last-days), building a dated series that dem_change.py "
                         "differences week to week.")
    ap.add_argument("--interpolate-edge", type=float, default=None,
                    help="Fill cells between waterlines by linear interpolation over a "
                         "triangulation of every crossed cell, using triangles with edges up to "
                         "this many metres (e.g. 6). Measured cells keep their value; "
                         "<stem>_source.asc marks 1 = measured, 2 = interpolated.")
    ap.add_argument("--fine-cell", type=float, default=None,
                    help="Also write <stem>_fine_{dem,resid,frames}.asc on this grid (e.g. 0.25 m), "
                         "fitted directly to the waterline points by local planes -- no 2 m "
                         "facets. See fine_dem().")
    ap.add_argument("--fine-window", default="4,1",
                    help="Fine DEM neighbourhood, metres 'alongshore,cross-shore' (default 4,1)")
    ap.add_argument("--fill-gaps", action="store_true",
                    help="Close isolated single-cell holes from neighbours. Off by default: an "
                         "interpolated cell looks identical to a measured one in the output.")
    ap.add_argument("--exclude", action="append", default=[],
                    help="Leave out frames whose file name matches this regular expression "
                         "(repeatable), e.g. 'Jan.22.*[.]c2[.]' -- for frames shown by their "
                         "overlays to have followed something other than the water's edge.")
    ap.add_argument("--no-plot", action="store_true",
                    help="Skip the page (<stem>_dem.png). The grids and <stem>_info.json are "
                         "still written; 'python3 dem_figure.py <stem>' draws the page later.")
    args = ap.parse_args()

    # The page the daily email attaches. Removed before anything else, so
    # no failure below -- no points, a crash, a page that cannot be drawn
    # -- can leave last run's page beside this run's grids looking current.
    # Every way through main() that is not --no-plot draws a new one.
    page = Path(str(args.output_stem) + "_dem.png")
    if not args.no_plot and page.exists():
        page.unlink()

    def no_points_page(message, stats=None):
        """Nothing to grid: the short page says so, with the window and filters."""
        if args.no_plot:
            return
        stats = stats or {}
        info = {"last_days": args.last_days, "window_start": args.start_date,
                "window_end": args.end_date, "camera": args.camera, "frames": {}, "days": 0,
                "cell": args.cell, "min_points": args.min_points, "max_spread": args.max_spread,
                "filters": {"max_hs": args.max_hs, "rough_frames": stats.get("rough_frames", 0),
                            "excluded_frames": stats.get("excluded_frames", 0)}}
        try:
            import dem_figure
            out = dem_figure.draw_message_page(
                args.output_stem, None, info, dem_figure.STATION, message,
                ["No DEM was built this run, so there is no map; the grids on disk, if any, "
                 "are from an earlier run.",
                 "Likely causes: the cameras, the detector or the water-level record were down "
                 "for the window, or the build's filters left out every frame."],
                cutoff=args.max_spread, min_points=args.min_points, cell=args.cell)
            print(f"wrote {out}  (no points: short page)")
        except Exception as exc:
            print(f"ERROR: the short page was not drawn either: {exc!r}")

    if args.last_days:
        end = args.end_date or latest_date(args.contour_csv, args.camera)
        if end is None:
            print("No georectified points in the file.")
            no_points_page("No waterline points in the file: nothing to map.")
            sys.exit(4)
        from datetime import date, timedelta
        args.end_date = end
        args.start_date = (date.fromisoformat(end) - timedelta(days=args.last_days - 1)).isoformat()
        print(f"Window            : last {args.last_days} day(s), {args.start_date} to {args.end_date}")

    load_stats = {}
    E, N, Z, cams, dates, frames, missing = load_points(
        args.contour_csv, args.camera, args.start_date, args.end_date, args.max_hs, args.exclude,
        stats=load_stats)
    rejected_frames, rejected_days = set(), set()
    day_record = {"camera_days_rejected_offsets_m": {}, "day_test": None}

    if len(E) == 0:
        print("No georectified points matched. Check --camera and the date range.")
        no_points_page("No waterline points in the window: nothing to map.", load_stats)
        sys.exit(4)

    print("=" * 74)
    print("INTERTIDAL DEM")
    print("=" * 74)
    print(f"points            : {len(E)}")
    if missing:
        print(f"  (skipped {missing} row(s) with no ground coordinates)")
    for c in sorted(set(cams)):
        print(f"    {c}: {int((cams == c).sum())}")
    print(f"dates             : {min(dates)} to {max(dates)}  "
          f"({len(set(dates))} day(s))")
    print(f"elevation range   : {Z.min():+.2f} to {Z.max():+.2f} m NAVD88")
    print(f"extent            : {E.max()-E.min():.1f} m E-W by {N.max()-N.min():.1f} m N-S")
    print()

    # Runs of a line across the shore (waterline_track.py): the line tracer carrying on along
    # one pixel row out over the water. The survey-date builds' consistency filter drops them
    # already (their files carry the ground); the live cron's filter runs on the photo rows
    # alone, so its DEMs get them dropped here.
    track_record = {"lines": 0, "points": 0, "by_camera": {}}
    if not args.no_track_filter:
        from waterline_track import cross_shore_drops, axis_bearing
        ucol = load_stats.get("pixel_column")
        if ucol is None or ucol.size != E.size:
            ucol = np.arange(E.size, dtype=float)
        ucol = np.where(np.isfinite(ucol), ucol, np.arange(E.size, dtype=float))
        drop_cs = np.zeros(E.size, bool)
        msgs = []
        for c in sorted(set(cams)):
            m = np.flatnonzero(cams == c)
            fid = np.unique(frames[m], return_inverse=True)[1]
            d, per_line, axis, ratio = cross_shore_drops(fid, ucol[m], E[m], N[m])
            if axis is None:
                msgs.append(f"{c} not judged (no clear shore axis, spread ratio {ratio:.1f})")
                continue
            if d.sum() > 0.5 * m.size:
                msgs.append(f"{c} NOT applied: would drop {100.0 * d.mean():.0f}% of its points "
                            f"(shore axis {axis_bearing(axis):.0f} deg wrong for this view?)")
                continue
            drop_cs[m[d]] = True
            track_record["by_camera"][c] = {"lines": len(per_line), "points": int(d.sum()),
                                            "shore_axis_deg": round(axis_bearing(axis), 1)}
            msgs.append(f"{c} {len(per_line)} line(s), {int(d.sum())} point(s)")
        track_record["lines"] = sum(v["lines"] for v in track_record["by_camera"].values())
        track_record["points"] = int(drop_cs.sum())
        print(f"Cross-shore runs  : {'; '.join(msgs)} dropped (a line's track running more than 45 deg "
              f"from the shore over 10 m: a trace out over the water, not a waterline; waterline_track.py)")
        if drop_cs.any():
            k = ~drop_cs
            E, N, Z, cams, dates, frames = E[k], N[k], Z[k], cams[k], dates[k], frames[k]
        print()

    if args.fit_setup or args.max_day_offset or args.max_frame_offset:
        info = load_frame_info(args.contour_csv)
        pair_cell, pair_frame = frame_pairs(E, N, frames, args.cell)

    if args.max_frame_offset:
        # Same test as the day check below, one frame at a time, and run first:
        # a single bad waterline (e.g. a high-tide line drawn seaward of the
        # low-tide ones) is removed on its own instead of taking its whole day
        # with it.
        z_of = {}
        for f, z in zip(frames, Z):
            z_of.setdefault(f, z)
        pair_vals = np.array([z_of[f] for f in pair_frame])
        offsets = day_offsets(pair_cell, pair_vals, pair_frame, args.min_points, min_pairs=10)
        bad = {k: o for k, (o, n) in offsets.items()
               if np.isfinite(o) and abs(o) > args.max_frame_offset}
        tested = sum(np.isfinite(o) for o, n in offsets.values())
        print(f"Frame consistency (median offset from the DEM without that frame, limit "
              f"+/-{args.max_frame_offset} m): {tested} of {len(offsets)} frames testable, "
              f"{len(bad)} rejected")
        for k in sorted(bad):
            print(f"   {k}   {bad[k]:+.3f} m  <-- REJECTED")
        rejected_frames = set(bad)
        if bad:
            keep = np.array([f not in bad for f in frames])
            E, N, Z, cams, dates, frames = E[keep], N[keep], Z[keep], cams[keep], dates[keep], frames[keep]
            print(f"   left out {len(bad)} frame(s), {int((~keep).sum())} point(s)")
            pair_cell, pair_frame = frame_pairs(E, N, frames, args.cell)
        print()

    if args.fit_setup:
        carried = carried_setup_coef(info)
        if args.setup_coef_in_use is not None:
            c_use, use_from = args.setup_coef_in_use, "--setup-coef-in-use"
        elif carried is not None:
            c_use, use_from = carried, "the contours' setup_correction_m"
        else:
            c_use, use_from = station_setup_coef(), "SETUP_COEF of waterline_timex_cron.sh"
        # the C in use is tried too (0.037 is not on the 0.002 grid), so it can be the best
        grid = np.round(np.arange(-0.04, 0.1201, 0.002), 3)
        if c_use and grid[0] < c_use < grid[-1]:
            grid = np.unique(np.r_[grid, round(c_use, 4)])
        result, n_pairs, n_frames = fit_setup_coefficient(pair_cell, pair_frame, info,
                                                          args.min_points, coefs=grid)
        valid = [r for r in result if np.isfinite(r[1])]
        if not valid:
            print("Setup fit         : no frames with both Hs and Tp -- run "
                  "extract_elevation_contours.py with --waves.")
        else:
            base = next(r for r in valid if r[0] == 0.0)
            best = min(valid, key=lambda r: r[1])
            print(f"Setup fit         : {n_frames} frame(s), {n_pairs} cell samples")
            if carried is not None:
                print(f"   NOTE: these contours carry a setup (C ~{carried}): each line was placed at "
                      f"still water + that setup. Fit on contours built with no setup (C = 0).")
            print(f"   C = 0      median spread {base[1]:.3f} m ({base[2]} cells)")
            print(f"   C = {best[0]:<6} median spread {best[1]:.3f} m  <-- best"
                  f"  (implied slope C/0.35 = {best[0] / 0.35:.3f})")
            for c, sp, n in valid:
                if round(c * 1000) % 20 == 0:
                    print(f"      C {c:+.2f}: {sp:.3f} m")

            def spread_at(c):
                hit = [r[1] for r in valid if abs(r[0] - c) < 1e-9]
                if hit:
                    return hit[0]
                return fit_setup_coefficient(pair_cell, pair_frame, info, args.min_points,
                                             coefs=[c])[0][0][1]
            # how sharp the minimum is, and whether C = 0 or the C in use do as well: on one
            # week of lines the curve is often flat and noisy, and its lowest point alone then
            # says little about C
            for line in setup_fit_verdict(valid, spread_at, c_use or None, use_from):
                print(line)
        print()

    if args.max_day_offset:
        epoch = {f: info[f]["epoch"] for f in info}
        cam_of = dict(zip(frames, cams))
        pair_key = np.array([f"{cam_of[f]} {local_day(epoch[f])}" for f in pair_frame])
        z_of = {}
        for f, z in zip(frames, Z):
            z_of.setdefault(f, z)
        pair_vals = np.array([z_of[f] for f in pair_frame])
        # Iterative (day_consistency): the day furthest over the limit is left out and the
        # others measured again without it; left-out days within the limit of the kept days
        # are taken back. With fewer than 3 testable camera-days "the DEM without that day"
        # is the other day alone: the test cannot tell which day is off, and keeps them all.
        # Days it cannot judge later on keep the first pass's verdict, with a WARNING.
        dt = day_consistency(pair_cell, pair_vals, pair_key, args.max_day_offset, args.min_points,
                             anchor=args.anchor_date)
        rejected = set(dt["rejected"])
        keep = None
        if rejected:
            frame_key = {f: f"{cam_of[f]} {local_day(epoch[f])}" for f in set(frames)}
            keep = np.array([frame_key[f] not in rejected for f in frames])
        for line in day_test_report(dt, args.max_day_offset,
                                    None if keep is None else int((~keep).sum())):
            print(line)
        # days over the limit it could not judge, or a day kept untested: one WARNING line
        for line in day_test_warning(dt, args.max_day_offset)[1]:
            print(line)
        # a camera whose day offsets run one way with time: the test may be trimming real change
        # (or the ends of an error drifting with time), so the DEM is not one surface -- say so
        trends = day_trends({k: v[0] for k, v in dt["final"].items()}, args.max_day_offset,
                            left_out=dt["rejected"])
        for t in trends:
            for line in day_trend_warning(t, args.max_day_offset, dt["rejected"]):
                if args.trend_note and line.startswith("WARNING: "):
                    line = "NOTE: " + line[len("WARNING: "):]   # the cron's DEMs: said, not counted
                print(line)
        day_record = day_test_record(dt, args.max_day_offset, trends)
        rejected_days = set(rejected)
        if keep is not None:
            E, N, Z, cams, dates, frames = E[keep], N[keep], Z[keep], cams[keep], dates[keep], frames[keep]
        print()
        if len(E) == 0:
            print("No points left: the day-consistency test rejected every camera-day "
                  f"({', '.join(sorted(rejected))}). Rerun without --max-day-offset or with a wider window.")
            no_points_page("The day-consistency test rejected every camera-day: nothing to map.", load_stats)
            sys.exit(4)

    dem, count, spread, e0, n0, ncols, nrows = build_grid(
        E, N, Z, args.cell, args.min_points, args.max_spread, frames)

    total_cells = dem.size
    with_any = int((count > 0).sum())
    filled = int(np.isfinite(dem).sum())
    blanked = with_any - filled

    print(f"grid              : {nrows} x {ncols} cells at {args.cell} m")
    print(f"cells with points : {with_any} ({100*with_any/total_cells:.1f}% of grid)")
    print(f"cells filled      : {filled}")
    if blanked:
        print(f"cells blanked     : {blanked}  (fewer than {args.min_points} frames, "
              f"or spread > {args.max_spread} m)")

    occupied = count[count > 0]
    if len(occupied):
        print(f"frames per cell   : median {int(np.median(occupied))}, "
              f"p10 {int(np.percentile(occupied,10))}, p90 {int(np.percentile(occupied,90))}")
    valid_spread = spread[np.isfinite(spread) & np.isfinite(dem)]
    if len(valid_spread):
        print(f"elevation spread  : median {np.median(valid_spread):.3f} m, "
              f"p90 {np.percentile(valid_spread,90):.3f} m")
        print("                    (repeat crossings of the same cell; this is the DEM's")
        print("                     own repeatability, not its accuracy)")

    source = None
    if args.interpolate_edge:
        dem, source = interpolate_between_contours(E, N, Z, frames, args.cell, args.max_spread,
                                                   dem, args.interpolate_edge)
        print(f"interpolated      : +{int((source == 2).sum())} cell(s) between waterlines "
              f"(triangle edges <= {args.interpolate_edge:g} m); {int((source == 1).sum())} measured")

    if args.fill_gaps:
        before = filled
        dem = fill_small_gaps(dem)
        after = int(np.isfinite(dem).sum())
        print(f"gap fill          : +{after-before} cell(s) interpolated from neighbours")

    stem = Path(args.output_stem)
    if args.fine_cell:
        along, across = (float(x) for x in args.fine_window.split(","))
        fdem, fres, fnf, fe0, fn0 = fine_dem(E, N, Z, frames, args.fine_cell, along, across,
                                             args.min_points, args.max_spread / 2)
        ok = np.isfinite(fdem)
        print(f"fine DEM          : {fdem.shape[0]} x {fdem.shape[1]} at {args.fine_cell:g} m, "
              f"{int(ok.sum())} nodes filled ({ok.sum() * args.fine_cell ** 2:.0f} m2); "
              f"window {along:g} m alongshore x {across:g} m cross-shore")
        if ok.any():
            print(f"                    fit residual median {np.nanmedian(fres[ok]):.3f} m, "
                  f"frames per node median {int(np.median(fnf[ok]))}")
        write_ascii_grid(str(stem) + "_fine_dem.asc", fdem, fe0, fn0, args.fine_cell)
        write_ascii_grid(str(stem) + "_fine_resid.asc", np.where(ok, fres, np.nan), fe0, fn0,
                         args.fine_cell)
        write_ascii_grid(str(stem) + "_fine_frames.asc", np.where(ok, fnf, 0.0), fe0, fn0,
                         args.fine_cell, nodata=0.0)
        print(f"wrote {stem}_fine_dem.asc (and _fine_resid, _fine_frames)")
    write_ascii_grid(str(stem) + "_dem.asc", dem, e0, n0, args.cell)
    write_ascii_grid(str(stem) + "_spread.asc", spread, e0, n0, args.cell)
    write_ascii_grid(str(stem) + "_count.asc",
                     count.astype(float), e0, n0, args.cell, nodata=0.0)
    if source is not None:
        write_ascii_grid(str(stem) + "_source.asc", source.astype(float), e0, n0, args.cell,
                         nodata=0.0)
    print()
    print(f"wrote {stem}_dem.asc     (elevation, m NAVD88)")
    print(f"wrote {stem}_spread.asc  (16-84 percentile range, m)")
    print(f"wrote {stem}_count.asc   (samples per cell)")

    sstem = None
    if args.series_dir:
        if not args.last_days:
            print("NOTE: --series-dir needs --last-days (the window defines the date); not saved.")
        else:
            sdir = Path(args.series_dir)
            sdir.mkdir(parents=True, exist_ok=True)
            sstem = sdir / f"dem_{args.end_date}_{args.last_days}d"
            write_ascii_grid(str(sstem) + "_dem.asc", dem, e0, n0, args.cell)
            write_ascii_grid(str(sstem) + "_spread.asc", spread, e0, n0, args.cell)
            write_ascii_grid(str(sstem) + "_count.asc", count.astype(float), e0, n0,
                             args.cell, nodata=0.0)
            print(f"series            : {sstem}_{{dem,spread,count}}.asc (+ _info.json)")

    # What the page needs and the grids do not hold: dates, frames, filters.
    # Written always, so 'python3 dem_figure.py <stem>' can redraw the page
    # later without a rebuild.
    page_info = {
        "first_date": str(min(dates)), "last_date": str(max(dates)),
        "days": int(len(set(dates))),
        "last_days": args.last_days,
        "window_start": args.start_date, "window_end": args.end_date,
        "camera": args.camera,
        "frames": {c: int(len(set(frames[cams == c]))) for c in sorted(set(cams))},
        "points": {c: int((cams == c).sum()) for c in sorted(set(cams))},
        "cell": args.cell, "min_points": args.min_points, "max_spread": args.max_spread,
        # share of the points whose elevation includes the wave setup
        # (beach_elevation_navd88, extract_elevation_contours.py --setup-coef)
        "setup_share": round(load_stats.get("setup_share", 0.0), 4),
        # True: this build wrote <stem>_source.asc (and the page may trust
        # it); False: any _source.asc beside the grids is an older run's.
        "interpolated": source is not None,
        "filters": {
            "max_hs": args.max_hs, "rough_frames": load_stats.get("rough_frames", 0),
            "excluded_frames": load_stats.get("excluded_frames", 0),
            "cross_shore_runs": track_record,
            "max_frame_offset": args.max_frame_offset,
            "frames_rejected": len(rejected_frames),
            "max_day_offset": args.max_day_offset,
            "camera_days_rejected": sorted(rejected_days),
            # with --max-day-offset: the offset of each day left out from the DEM of the kept days
            # (camera_days_rejected_offsets_m), and the test's numbers (day_test: first-pass and
            # final offsets of every camera-day, the order left out, the days taken back and those
            # refused, the days it could not judge, its WARNING, trend warnings); a page can name
            # the days and say why. day_test is None without the test.
            **day_record,
        },
    }
    with open(str(stem) + "_info.json", "w") as f:
        json.dump(page_info, f, indent=1)
    print(f"wrote {stem}_info.json   (dates, frames, filters: for the page)")
    if sstem is not None:
        # the dated copy gets its own, so a page redrawn from it (dem_figure.py
        # <series stem>) credits only the cameras that measured: without it a
        # c1-only week's page said 'seen by c2' and drew a seam
        with open(str(sstem) + "_info.json", "w") as f:
            json.dump(page_info, f, indent=1)

    rc = 0
    if not args.no_plot:
        # A page that cannot be drawn is a failure the cron must see (exit 3:
        # grids written, page not), not a note: the email would otherwise
        # attach whatever page was there before. dem_figure draws a short
        # page itself when there are too few cells to map.
        try:
            import dem_figure
            out = dem_figure.draw_page(stem, info=page_info)
            print(f"wrote {out}")
        except Exception as exc:
            import traceback
            traceback.print_exc()
            print(f"ERROR: page NOT drawn ({exc!r}); the grids and {stem}_info.json are "
                  f"written. No {page.name} this run. Redraw: python3 dem_figure.py {stem}")
            if page.exists():
                page.unlink()                    # a half-written page is not a page
            rc = 3

    print()
    print("REMINDER: intertidal zone only, between the lowest and highest water")
    print("levels observed. No data above or below, and none invented.")
    share = load_stats.get("setup_share", 0.0)
    rough = (f"--max-hs {args.max_hs:g} left out the roughest frames" if args.max_hs is not None
             else "--max-hs would leave out the roughest frames")
    if share > 0:
        print("Elevations include the wave setup (C from extract_elevation_contours.py")
        print(f"--setup-coef) on {100 * share:.0f}% of the points"
              + ("." if share >= 0.995 else " (the rest: no wave record, water level only)."))
        print(f"{rough}, where the setup errs most.")
    else:
        print("The detected edge sits above still water by the wave setup, which grows")
        print("with wave height, but is given the still-water elevation -- so the DEM")
        print("reads LOW where rough-water frames land. Not corrected here (no setup in")
        print(f"these contours); {rough}.")
    if rc:
        sys.exit(rc)


if __name__ == "__main__":
    main()

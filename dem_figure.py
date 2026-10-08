#!/usr/bin/env python3
"""
DEM Page: The Intertidal DEM In The Beach's Own Frame
========================================================
Draws <stem>_dem.png -- the page the daily email attaches -- from the
grids dem_from_contours.py writes next to it (<stem>_dem.asc,
_spread.asc, _count.asc and, with --interpolate-edge, _source.asc).
dem_from_contours.py calls it after every build; run it on its own to
redraw a page from grids already on disk, without rebuilding the DEM:

    python3 dem_figure.py dem_intertidal_7day
    python3 dem_figure.py dem_intertidal --profiles 40,150,260,350
    python3 dem_figure.py dem_intertidal_7day --output /tmp/page.png

WHY A NEW PAGE. The old page drew the grids north-up with equal axes.
The DEM is a strip about 400 m long and 30 m wide running N-NW to
S-SE, so it came out as a thin diagonal sliver with most of the page
empty; the 'terrain' colormap painted sand below 0 m blue and green, so
the lower beach read as water and vegetation; and nothing said which
camera had seen which part, or why a cell was empty -- a gap where no
waterline passed and a cell blanked because its crossings disagreed
looked the same. On the station's 7-day DEM of 29 Sep - 5 Oct 2026 the
southern ~70 m, seen by c1 alone, was mostly blanked, and the old page
could not show that, or why.

WHAT IT SHOWS, AND WHY EACH CHOICE
  * The beach frame. x is alongshore distance from the camera along the
    strip's principal axis; y is cross-shore distance seaward of the
    camera. The maps are drawn as seen from the bluff looking out to
    sea -- sea at the top, the N-NW end (c2's view) on the left and the
    S-SE end (c1's) on the right. That is a rotation of the map, not a
    mirror image, so north is where the arrow says and the alongshore
    axis counts up from right to left. The cross-shore axis is
    stretched by the largest round factor that still fits the strip,
    so it fills the page width; the factor is on the axis
    ('cross-shore x2.5'). Slopes look that much steeper on the maps;
    the profiles (third panel) are the place to read slopes.
  * Elevation in one hue, light = high (dry sand) to dark = low, in
    even lightness steps. A single-hue ramp has one reading -- darker is
    lower -- and no colour that looks like water or vegetation.
    Contours every 0.5 m, traced on a lightly smoothed copy of the grid:
    3 x 3 mean of measured cells, where a cell without a value joins in
    only if at least 5 of its 8 neighbours are measured -- so a contour
    is not cut at every isolated blanked cell, but still breaks at a
    real gap. Pieces shorter than 12 m are left out (a one-cell wiggle
    reads as a dash, and dashes mean photo edges here). Labels sit in a
    gap in their line, wider for a negative value so the minus sign is
    not read as the line's end, and never on a profile line, a photo
    edge or another label. The cells themselves are drawn unsmoothed.
  * Repeatability: the 16-84 percentile spread of the frames in each
    cell, the quantity --max-spread tests, in a second hue on the same
    frame, with the cutoff as the top of its colour scale; cells above
    it are the darkest step and hatched, as they are blanked in the
    elevation map. It is the DEM's own repeatability, not its accuracy
    -- and real change inside the window (a storm cutting the upper
    beach) shows up here too.
  * Empty cells say why they are empty: light grey = crossed by fewer
    than --min-points frames, so no value; hatched = enough frames but
    their spread exceeded --max-spread, so the value was blanked;
    white = no waterline at all. Cells filled by interpolation
    (--interpolate-edge, --fill-gaps) are dotted.
  * The cameras. Position and image edges come from the calibration
    files georectify.py uses (calibration/CACO05_<cam>_20240801_IO.yaml
    and _20251113_EO-CV.yaml); the edges are drawn dashed where they
    cross the strip, the cameras as a triangle at alongshore 0, and a
    bar under the elevation map marks the stretch of beach each camera
    covers. c1's left image edge and c2's right image edge cross the
    beach close together: the seam is the middle of the cells both
    photos hold, and the two sides of it are measured by different
    cameras. Calibration more than 3 km from the DEM is another site's
    and is not drawn.
  * Profiles: 3 or 4 cross-shore profiles at lettered alongshore
    positions (lines on both maps), each the cells within one cell of
    that position, with the foreshore slope as the least-squares line
    through them, given as tan(beta) and 1:N. Four at most: the lines
    cross, so every pair of colours must stay distinct, and a fifth
    categorical colour fails that check against the others. Positions
    are picked automatically -- in the middle 60% of each stretch, at
    least one in each camera's, the profile with the largest elevation
    range along a clean seaward-falling line -- unless --profiles names
    them.
  * The title says what and when: the dates and number of days, the
    frames per camera, the cells filled, blanked and without a value,
    and the frames the build's filters left out. Facts are wrapped
    between, never inside, to the page width; so is the legend.

WHERE THE FACTS COME FROM. Dates, frames and filters are not in the
grids: dem_from_contours.py writes them to <stem>_info.json. Without
that file (grids from an older run) the page is drawn from what the
grids hold and the title says the dates are unknown.

The page is 14 in wide at 150 dpi (2100 px). Text is 9.5 pt or larger,
about 13 px when the email shrinks the page to 1400 px.

Usage:
    python3 dem_figure.py <stem> [--output PNG] [--profiles 40,150,260]
        [--n-profiles 4] [--calibration DIR] [--station NAME]
        [--max-spread 0.5] [--min-points 3]
"""

import sys
import json
import argparse
import warnings
from datetime import date
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STATION = "CACO05 Marconi Beach"
CAMERAS = ("c1", "c2")
IO_FILE = "CACO05_{cam}_20240801_IO.yaml"
EO_FILE = "CACO05_{cam}_20251113_EO-CV.yaml"

# Ink and surfaces: neutral, so the only colour on the page is data.
INK = "#1f1e1c"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e6e5df"
AXIS = "#c3c2b7"
NO_VALUE = "#cfcfcb"
BLANK_FACE = "#f2f1ed"
BLANK_HATCH = "#5f5e58"
# Elevation: one hue (OKLCH h 78, sand), lightness 0.93 (high) -> 0.34 (low).
SAND = ["#f7e5cb", "#e0c399", "#c5a36e", "#a7844d", "#876834", "#674d21", "#483413"]
# Spread: one hue (OKLCH h 295, violet), lightness 0.95 (repeatable) -> 0.38.
VIOLET = ["#efecfb", "#cbc3e7", "#a89cd3", "#8874bd", "#684fa3", "#492e7f"]
VIOLET_OVER = "#2b1850"
# Profiles: categorical slots (blue, orange, aqua, violet). The lines
# cross, so every pair must be told apart, not just neighbours: these
# four pass the all-pairs colour-blind and normal-vision checks; the
# reference palette's 4th slot (yellow) does not, next to orange. Four
# is therefore the most profiles drawn. Each line also carries its
# letter, so colour is never the only key.
PROFILE_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"]
MAX_PROFILES = len(PROFILE_COLOURS)
# Cameras further than this from the DEM are not this site's (another
# station's grids drawn with this folder's calibration): not drawn.
MAX_CAMERA_DISTANCE = 3000.0
# "How to read this page": 10 pt at 1.35 line spacing, and the gap between paragraphs (in).
NOTE_LINE, NOTE_GAP = 0.19, 0.13
PAGE_RC = {"font.size": 10, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
           "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
           "hatch.linewidth": 0.9, "axes.linewidth": 0.8}
# Contour pieces shorter than this (m along the page, cross-shore stretched) are not drawn.
CONTOUR_MIN_M = 12.0
NICE_EXAGGERATION = (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 15, 20)
COMPASS = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
           "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]


# ---------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------

def read_asc(path):
    """ESRI ASCII grid -> (grid with row 0 = south, as dem_from_contours builds it, header)."""
    hdr = {}
    with open(path) as f:
        for _ in range(6):
            k, v = f.readline().split()
            hdr[k.lower()] = float(v)
    g = np.loadtxt(path, skiprows=6, ndmin=2)
    g = np.where(g == hdr.get("nodata_value", -9999.0), np.nan, g)
    return np.flipud(g), hdr


def _same_grid(a, b):
    return all(abs(a[k] - b[k]) < 1e-6 for k in ("ncols", "nrows", "xllcorner", "yllcorner",
                                                 "cellsize"))


def load_grids(stem, interpolated=None):
    """
    The DEM's grids. `interpolated`: True/False when the caller knows
    whether this build interpolated (it then trusts or ignores
    <stem>_source.asc); None (no info) uses the file only if it is on the
    same grid and not older than the DEM, so a _source.asc left from an
    earlier --interpolate-edge run cannot mark this run's cells.
    """
    stem = str(stem)
    dem, h = read_asc(stem + "_dem.asc")
    spread, hs = read_asc(stem + "_spread.asc")
    count, hc = read_asc(stem + "_count.asc")
    for name, hh in (("spread", hs), ("count", hc)):
        if not _same_grid(h, hh):
            raise ValueError(f"{stem}_{name}.asc is not on the same grid as {stem}_dem.asc")
    source = None
    sp = Path(stem + "_source.asc")
    if sp.exists() and interpolated is not False:
        src, hsrc = read_asc(sp)
        fresh = interpolated is True or \
            sp.stat().st_mtime >= Path(stem + "_dem.asc").stat().st_mtime - 120
        if _same_grid(h, hsrc) and fresh:
            source = np.nan_to_num(src, nan=0.0).astype(int)
    return {"dem": dem, "spread": spread, "count": np.nan_to_num(count, nan=0.0).astype(int),
            "source": source, "e0": h["xllcorner"], "n0": h["yllcorner"],
            "cell": h["cellsize"]}


def load_info(stem):
    p = Path(str(stem) + "_info.json")
    if not p.exists():
        return {}
    try:
        with open(p) as f:
            return json.load(f)
    except (OSError, ValueError) as exc:
        print(f"WARNING: {p.name} unreadable ({exc}); drawing without it")
        return {}


def load_cameras(cal_dir):
    """{cam: (io, eo)} for each camera whose calibration is found. {} if none."""
    try:
        from georectify import load_intrinsics, load_extrinsics
    except Exception as exc:                      # pragma: no cover - station has it
        print(f"NOTE: cameras not drawn (georectify.py: {exc})")
        return {}
    out = {}
    for cam in CAMERAS:
        io_p = Path(cal_dir) / IO_FILE.format(cam=cam)
        eo_p = Path(cal_dir) / EO_FILE.format(cam=cam)
        if io_p.exists() and eo_p.exists():
            try:
                out[cam] = (load_intrinsics(io_p), load_extrinsics(eo_p))
            except ValueError as exc:
                print(f"NOTE: {cam} not drawn ({exc})")
    if not out:
        print(f"NOTE: no calibration in {cal_dir} -- cameras and seam not drawn")
    return out


# ---------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------

def beach_frame(E, N, Z, camera_xy):
    """
    Alongshore unit vector `a` (the strip's principal axis, pointing
    from the camera towards most of the strip), seaward unit vector `s`
    (the side where the beach is lower), and the origin (the camera, or
    without one the strip's near end and landward edge).
    `flip_x`: with the sea drawn at the top, a true (unmirrored)
    rotation has the alongshore axis increasing to the left.
    """
    ec, nc = E.mean(), N.mean()
    _, vecs = np.linalg.eigh(np.cov(np.stack([E - ec, N - nc])))
    a = vecs[:, 1].copy()
    s = np.array([a[1], -a[0]])
    ok = np.isfinite(Z)
    if ok.sum() >= 10:
        v = (E[ok] - ec) * s[0] + (N[ok] - nc) * s[1]
        if np.polyfit(v, Z[ok], 1)[0] > 0:
            s = -s
    elif camera_xy is not None and ((ec - camera_xy[0]) * s[0] + (nc - camera_xy[1]) * s[1]) < 0:
        s = -s
    origin = np.array(camera_xy if camera_xy is not None else (ec, nc), float)
    u = (E - origin[0]) * a[0] + (N - origin[1]) * a[1]
    if np.median(u) < 0:
        a = -a
        u = -u
    if camera_xy is None:
        v = (E - origin[0]) * s[0] + (N - origin[1]) * s[1]
        origin = origin + a * u.min() + s * v.min()
    flip_x = (a[0] * s[1] - a[1] * s[0]) < 0
    return {"a": a, "s": s, "origin": origin, "flip_x": bool(flip_x),
            "from_camera": camera_xy is not None}


def to_uv(frame, E, N):
    dE, dN = np.asarray(E) - frame["origin"][0], np.asarray(N) - frame["origin"][1]
    return (dE * frame["a"][0] + dN * frame["a"][1],
            dE * frame["s"][0] + dN * frame["s"][1])


def compass(vec):
    bearing = np.degrees(np.arctan2(vec[0], vec[1])) % 360
    return COMPASS[int(round(bearing / 22.5)) % 16]


def map_limits(u, v, width_in, height_in, pad=4.0):
    """
    Axis limits that fill a width_in x height_in box, and the cross-shore
    stretch: the largest round factor at which the strip still fits.
    The cross-shore extent is taken between the 0.25 and 99.75
    percentiles, so a few stray cells cannot shrink the stretch.
    """
    u_lo, u_hi = u.min() - pad, u.max() + pad
    v_lo, v_hi = np.percentile(v, [0.25, 99.75])
    v_lo, v_hi = v_lo - pad / 2, v_hi + pad / 2
    lu, lv = u_hi - u_lo, v_hi - v_lo
    bound = (height_in / width_in) * lu / lv
    ex = max([x for x in NICE_EXAGGERATION if x <= bound] or [1])
    v_span = lu * (height_in / width_in) / ex
    if v_span < lv:                    # wider than long even at 1:1 -- widen alongshore
        u_mid, u_span = (u_lo + u_hi) / 2, lv * ex * width_in / height_in
        u_lo, u_hi = u_mid - u_span / 2, u_mid + u_span / 2
        v_span = lv
    v_mid = (v_lo + v_hi) / 2
    return (u_lo, u_hi), (v_mid - v_span / 2, v_mid + v_span / 2), ex


def camera_view(cams, frame, Ec, Nc, Zc, core, z_edge):
    """
    For each camera: which cells of `core` (cells with enough frames) it
    sees -- the cell centre at its own elevation projects inside the
    photo -- and its image edges on the ground at z_edge, in beach
    coordinates, kept only where they cross the strip (core grown by 3
    cells).
    """
    from view_reproject import ground_to_pixel
    from georectify import pixel_to_ground
    out = {}
    nrows, ncols = core.shape
    grown = core.copy()
    for _ in range(3):
        g = np.pad(grown, 1)
        grown = g[1:-1, 1:-1] | g[:-2, 1:-1] | g[2:, 1:-1] | g[1:-1, :-2] | g[1:-1, 2:]
    for cam, (io, eo) in cams.items():
        _, _, seen = ground_to_pixel(Ec.ravel(), Nc.ravel(), Zc.ravel(), io, eo)
        nu, nv = io[0], io[1]
        t = np.linspace(0, 1, 600)
        sides = {"left": (np.zeros_like(t), t * (nv - 1)),
                 "right": (np.full_like(t, nu - 1), t * (nv - 1)),
                 "top": (t * (nu - 1), np.zeros_like(t)),
                 "bottom": (t * (nu - 1), np.full_like(t, nv - 1))}
        edges = {}
        for side, (pu, pv) in sides.items():
            E, N = pixel_to_ground(pu, pv, z_edge, io, eo)
            r = np.floor((N - frame["n0"]) / frame["cell"])
            c = np.floor((E - frame["e0"]) / frame["cell"])
            inside = np.isfinite(r) & (r >= 0) & (r < nrows) & (c >= 0) & (c < ncols)
            keep = np.zeros_like(inside)
            keep[inside] = grown[r[inside].astype(int), c[inside].astype(int)]
            if keep.sum() >= 2:
                eu, ev = to_uv(frame, E, N)
                edges[side] = (np.where(keep, eu, np.nan), np.where(keep, ev, np.nan))
        out[cam] = {"seen": seen.reshape(core.shape) & core, "edges": edges}
    return out


def camera_stretches(views, uc):
    """
    The stretch of beach each camera covers, split at the seam: [(cam,
    lo, hi)] in alongshore metres, and the seam position (None for one
    camera). Where two views overlap (the seam is oblique, so they
    overlap a little), the seam is the middle of the cells both see.
    """
    seen = {c: v["seen"] for c, v in views.items() if v["seen"].any()}
    if not seen:
        return [], None
    if len(seen) == 1:
        (c, m), = seen.items()
        return [(c, float(uc[m].min()), float(uc[m].max()))], None
    (c1, m1), (c2, m2) = sorted(seen.items())[:2]
    if np.median(uc[m1]) > np.median(uc[m2]):
        (c1, m1), (c2, m2) = (c2, m2), (c1, m1)      # c1 now the camera on the low-u side
    both = m1 & m2
    if both.any():
        seam = float(np.median(uc[both]))
    else:
        seam = float((uc[m1 & ~m2].max() + uc[m2 & ~m1].min()) / 2)
    return [(c1, float(uc[m1].min()), seam), (c2, seam, float(uc[m2].max()))], seam


# ---------------------------------------------------------------------
# Profiles and contours
# ---------------------------------------------------------------------

def pick_profiles(uc, vc, z, cell, n, stretches=()):
    """
    The most complete clean profile in the middle 60% of each of n
    stretches of the strip, so neighbours stay apart. With camera
    stretches, each camera's stretch gets at least one (so the part only
    c1 sees is always profiled), the rest in proportion to length;
    otherwise n equal stretches.

    'Most complete clean': the largest elevation range explained by a
    seaward-falling line (its slope times its cross-shore length) less
    three times the scatter about it. The plain elevation range would
    prefer a profile with a stray high cell at its seaward end -- the
    kind of wrong waterline a profile should not be chosen to show.
    """
    ok = np.isfinite(z)
    if ok.sum() < 10:
        return []
    lo, hi = np.percentile(uc[ok], [1, 99])
    parts = [(max(a, lo), min(b, hi)) for _, a, b in stretches if min(b, hi) - max(a, lo) > 4 * cell]
    if not parts or len(parts) > n:
        parts = [(lo, hi)]
    total = sum(b - a for a, b in parts)
    share = [max(1, int(round(n * (b - a) / total))) for a, b in parts]
    while sum(share) > n:
        share[int(np.argmax(share))] -= 1
    while sum(share) < n:
        share[int(np.argmax([(b - a) / s for (a, b), s in zip(parts, share)]))] += 1
    segments = []
    for (a, b), k in zip(parts, share):
        e = np.linspace(a, b, k + 1)
        segments += list(zip(e[:-1], e[1:]))
    picks = []
    for a, b in sorted(segments):
        best, best_score = None, -np.inf
        mid = (a + b) / 2
        # the middle 60% of each stretch, so neighbouring picks stay apart
        lo_p, hi_p = mid - 0.3 * (b - a), mid + 0.3 * (b - a)
        for p in np.arange(np.ceil(lo_p), hi_p + 1e-9, cell / 2):
            sel = ok & (np.abs(uc - p) <= cell)
            if sel.sum() < 5 or np.ptp(vc[sel]) < 2 * cell:
                continue
            slope, icpt = np.polyfit(vc[sel], z[sel], 1)
            rms = np.sqrt(np.mean((z[sel] - (icpt + slope * vc[sel])) ** 2))
            score = -slope * np.ptp(vc[sel]) - 3 * rms - 0.001 * abs(p - mid)
            if score > best_score:
                best, best_score = p, score
        if best is not None:
            picks.append(float(np.round(best)))
    return picks


def profile_at(uc, vc, z, p, cell):
    """Cells within one cell of alongshore position p: binned profile and fitted slope."""
    sel = np.isfinite(z) & (np.abs(uc - p) <= cell)
    if sel.sum() < 4 or np.ptp(z[sel]) < 0.2:
        return None
    v, zz = vc[sel], z[sel]
    k = np.round(v / cell).astype(int)
    keys = np.unique(k)
    bv = np.array([v[k == q].mean() for q in keys])
    bz = np.array([np.median(zz[k == q]) for q in keys])
    slope, icpt = np.polyfit(v, zz, 1)
    return {"v": v, "z": zz, "bv": bv, "bz": bz, "tanb": -slope, "icpt": icpt,
            "n": int(sel.sum())}


def smooth_for_contours(dem, need=5):
    """
    3 x 3 mean over measured cells, for tracing contours only. A cell
    with no value of its own gets one only when at least `need` of its 8
    neighbours are measured -- a single blanked or thin cell inside the
    measured beach -- so a contour is not cut at every isolated hole,
    while a real gap (a few cells wide, or the strip's edge) still
    breaks it. The cells themselves are drawn unsmoothed and unfilled.
    """
    nr, nc = dem.shape
    p = np.pad(dem, 1, constant_values=np.nan)
    stack = np.stack([p[a:a + nr, b:b + nc] for a in range(3) for b in range(3)])
    n = np.isfinite(stack).sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        m = np.nanmean(stack, axis=0)
    return np.where(np.isfinite(dem) | (n >= need), m, np.nan)


def contour_segments(x, y, z, levels):
    """
    [(level, [(k, 2) arrays])] traced on a throwaway figure. The page
    draws the lines itself, so it can drop fragments and cut the gaps
    for labels the same way on every matplotlib from 3.3 on.
    """
    import matplotlib.pyplot as plt
    fig = plt.figure()
    try:
        ax = fig.add_subplot(111)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cs = ax.contour(x, y, np.ma.masked_invalid(z), levels=levels)
        return [(float(L), _all_segments(cs, i)) for i, L in enumerate(cs.levels)]
    finally:
        plt.close(fig)


def _cut_for_label(P, centre, half):
    """
    Polyline P (px) split around arc length `centre` +/- `half`: the
    pieces before and after, and the label angle (deg, kept upright).
    None when the polyline is too short on either side.
    """
    d = np.r_[0.0, np.cumsum(np.hypot(*np.diff(P, axis=0).T))]
    if centre - half < 0 or centre + half > d[-1]:
        return None
    a = np.array([np.interp(centre - half, d, P[:, 0]), np.interp(centre - half, d, P[:, 1])])
    b = np.array([np.interp(centre + half, d, P[:, 0]), np.interp(centre + half, d, P[:, 1])])
    before = np.vstack([P[d < centre - half], a])
    after = np.vstack([b, P[d > centre + half]])
    ang = np.degrees(np.arctan2(b[1] - a[1], b[0] - a[0]))
    if ang > 90:
        ang -= 180
    elif ang <= -90:
        ang += 180
    return before, after, ang


def _luminance(rgb):
    lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in rgb[:3]]
    return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]


def _all_segments(cs, i):
    """Contour paths of level i as (k, 2) arrays (matplotlib 3.3 .. 3.11)."""
    try:
        return [s for s in cs.allsegs[i] if len(s) >= 2]
    except (AttributeError, IndexError):     # pragma: no cover
        return []


# ---------------------------------------------------------------------
# Text
# ---------------------------------------------------------------------

def _day(d):
    return date.fromisoformat(d)


def date_range_text(first, last):
    a, b = _day(first), _day(last)
    if a == b:
        return f"{b.day} {b:%b %Y}"
    if (a.year, a.month) == (b.year, b.month):
        return f"{a.day}–{b.day} {b:%b %Y}"
    if a.year == b.year:
        return f"{a.day} {a:%b} – {b.day} {b:%b %Y}"
    return f"{a.day} {a:%b %Y} – {b.day} {b:%b %Y}"


def minus(x, fmt="{:.1f}"):
    """Number with a typographic minus, and never '-0.0'."""
    s = fmt.format(x)
    if float(s) == 0:
        s = fmt.format(0.0)
    return s.replace("-", "−")


def title_lines(info, cells, station, cutoff, min_points, cell):
    """
    The heading. Line 1 says what and when; then two lists of short
    facts -- what went into the grid (frames, cells) and what the build's
    filters left out -- which draw_page joins with ' · ' and wraps to the
    page width, so a long list breaks between facts, not inside one.
    """
    if info.get("first_date") and info.get("last_date"):
        if info.get("last_days") and info.get("window_start") and info.get("window_end"):
            when = date_range_text(info["window_start"], info["window_end"])
        else:
            when = date_range_text(info["first_date"], info["last_date"])
        if info.get("last_days"):
            span = f"last {info['last_days']} days"
        elif info.get("window_start") or info.get("window_end"):
            span = "selected dates"
        else:
            span = "whole archive"
        days = info.get("days")
        frames = info.get("frames") or {}
        nfr = sum(frames.values())
        per_cam = ", ".join(f"{c} {n:,}" for c, n in sorted(frames.items()))
        if info.get("last_days") and days and days < info["last_days"]:
            day_txt = f"on {days} of {info['last_days']} days"
        else:
            day_txt = f"on {days} day{'s' if days != 1 else ''}"
        line1 = f"Intertidal DEM, {station}: {when} ({span})"
        facts = [f"{nfr:,} frames {day_txt} ({per_cam})"]
    else:
        line1 = f"Intertidal DEM, {station}: dates unknown"
        facts = ["dates and frames not recorded (no <stem>_info.json)"]
    facts += [f"{cells['filled']:,} cells of {cell:g} × {cell:g} m with an elevation",
              f"{cells['blanked']:,} blanked: spread > {cutoff:g} m",
              f"{cells['no_value']:,} crossed by fewer than {min_points} frames"]
    if cells.get("interpolated"):
        facts.append(f"{cells['interpolated']:,} interpolated between waterlines")
    left_out = []
    f = info.get("filters") or {}
    if f.get("max_hs") is not None:
        left_out.append(f"{f.get('rough_frames', 0):,} frames with offshore Hs > {f['max_hs']:g} m")
    if f.get("max_frame_offset") is not None:
        left_out.append(f"{f.get('frames_rejected', 0):,} frames off the rest by > "
                        f"{f['max_frame_offset']:g} m")
    if f.get("max_day_offset") is not None:
        n = len(f.get("camera_days_rejected") or [])
        left_out.append(f"{n} camera-day{'' if n == 1 else 's'} off the rest by > "
                        f"{f['max_day_offset']:g} m")
    if f.get("excluded_frames"):
        left_out.append(f"{f['excluded_frames']:,} frames by name (--exclude)")
    if left_out:
        left_out[0] = "Left out before gridding: " + left_out[0]
    return line1, facts, left_out


def wrap_facts(facts, max_px, width_px, sep="  ·  "):
    """Greedy: as many facts per line as fit in max_px; a fact is never split."""
    lines, cur = [], ""
    for f in facts:
        trial = f if not cur else cur + sep + f
        if cur and width_px(trial) > max_px:
            lines.append(cur)
            cur = f
        else:
            cur = trial
    if cur:
        lines.append(cur)
    return lines


def wrap_words(text, max_px, width_px):
    """Greedy word wrap to a measured width (px), not a character count."""
    lines, cur = [], ""
    for w in text.split():
        trial = w if not cur else cur + " " + w
        if cur and width_px(trial) > max_px:
            lines.append(cur)
            cur = w
        else:
            cur = trial
    if cur:
        lines.append(cur)
    return lines


def place_end_labels(ends_px, size_px, bounds_px, dx_px, gap_px=3.0):
    """
    Label centres (px) for letters at the ends of lines. Each goes dx_px
    to the right of its line's end; labels that would overlap are
    gathered into one column right of the rightmost end among them,
    stacked in the order of their ends (top end, top label, so leaders
    do not cross) and centred on them. Kept inside bounds_px (y0, y1).
    """
    ends = np.asarray(ends_px, float)
    w, h = size_px
    y0, y1 = bounds_px
    lab = ends + [dx_px, 0.0]
    group = list(range(len(ends)))

    def root(i):
        while group[i] != i:
            i = group[i]
        return i

    for _ in range(len(ends) + 1):
        merged = False
        for i in range(len(ends)):
            for j in range(i + 1, len(ends)):
                if root(i) != root(j) and abs(lab[i, 0] - lab[j, 0]) < w + gap_px \
                        and abs(lab[i, 1] - lab[j, 1]) < h + gap_px:
                    group[root(j)] = root(i)
                    merged = True
        for g in {root(i) for i in range(len(ends))}:
            members = [i for i in range(len(ends)) if root(i) == g]
            if len(members) < 2:
                continue
            members.sort(key=lambda i: -ends[i, 1])
            step = h + gap_px
            ys = ends[members, 1].mean() + step * ((len(members) - 1) / 2 - np.arange(len(members)))
            shift = max(0.0, (y0 + h / 2) - ys.min()) - max(0.0, ys.max() - (y1 - h / 2))
            ys += shift
            x = ends[members, 0].max() + dx_px
            for i, y in zip(members, ys):
                lab[i] = (x, y)
        if not merged:
            break
    lab[:, 1] = np.clip(lab[:, 1], y0 + h / 2, y1 - h / 2)
    return lab


# ---------------------------------------------------------------------
# The page
# ---------------------------------------------------------------------

def draw_page(stem, output=None, info=None, calibration=None, profiles=None, n_profiles=4,
              station=STATION, max_spread=None, min_points=None, dpi=150):
    """
    Draws the DEM page for <stem> from its grids; returns the PNG path.
    `info`: what dem_from_contours.py knows about the build (dates,
    frames, filters); None reads <stem>_info.json. The page's styling is
    set for this figure only, not left in matplotlib's global settings.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with plt.rc_context(PAGE_RC):
        return _draw_page(stem, output, info, calibration, profiles, n_profiles, station,
                          max_spread, min_points, dpi)


def _draw_page(stem, output, info, calibration, profiles, n_profiles, station, max_spread,
               min_points, dpi):
    import matplotlib.pyplot as plt
    from matplotlib import patheffects
    from matplotlib.backends.backend_agg import RendererAgg
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.font_manager import FontProperties
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    stem = str(stem)
    if info is None:
        info = load_info(stem)
    g = load_grids(stem, interpolated=info.get("interpolated"))
    dem, spread, count, source = g["dem"], g["spread"], g["count"], g["source"]
    cell, e0, n0 = g["cell"], g["e0"], g["n0"]
    nrows, ncols = dem.shape
    cutoff = float(max_spread if max_spread is not None else (info.get("max_spread") or 0.5))
    min_points = int(min_points if min_points is not None else info.get("min_points", 3))

    data = count > 0
    filled = np.isfinite(dem)
    blanked = data & ~filled & np.isfinite(spread)
    no_value = data & ~filled & ~np.isfinite(spread)
    # A value with no spread of its own was not measured in that cell:
    # interpolated (--interpolate-edge, source 2) or gap-filled (--fill-gaps).
    # So was a value over the cutoff: the build blanked it, then a gap
    # fill put a value back.
    interpolated = filled & ~np.isfinite(spread)
    interpolated |= filled & np.isfinite(spread) & (spread > cutoff)
    if source is not None:
        interpolated |= filled & (source == 2)
    core = filled | blanked                 # cells with enough frames to say something
    if core.sum() < 10:
        raise ValueError(f"{stem}: fewer than 10 cells with {min_points} or more frames")
    cells = {"filled": int((filled & ~interpolated).sum()), "blanked": int(blanked.sum()),
             "no_value": int(no_value.sum()), "interpolated": int(interpolated.sum())}

    # cell centres and corners
    rr, cc = np.mgrid[0:nrows, 0:ncols]
    Ec, Nc = e0 + (cc + 0.5) * cell, n0 + (rr + 0.5) * cell
    Ek, Nk = np.meshgrid(e0 + np.arange(ncols + 1) * cell, n0 + np.arange(nrows + 1) * cell)

    # cameras and the beach frame
    cams = load_cameras(calibration or HERE / "calibration")
    framed = info.get("frames") or {}
    if framed:
        cams = {c: v for c, v in cams.items() if framed.get(c, 0) > 0}
    if cams:
        mid = np.array([Ec[core].mean(), Nc[core].mean()])
        far = {c: float(np.hypot(*(v[1][:2] - mid))) for c, v in cams.items()}
        if max(far.values()) > MAX_CAMERA_DISTANCE:
            print(f"NOTE: cameras not drawn -- the calibration's cameras are "
                  f"{max(far.values()) / 1000:.1f} km from this DEM (another site?)")
            cams = {}
    cam_xy = np.mean([v[1][:2] for v in cams.values()], axis=0) if cams else None
    frame = beach_frame(Ec[core], Nc[core], dem[core], cam_xy)
    frame.update(e0=e0, n0=n0, cell=cell)
    uc, vc = to_uv(frame, Ec, Nc)
    uk, vk = to_uv(frame, Ek, Nk)
    zc = np.where(filled, dem, np.nan)
    z_edge = float(np.nanmedian(dem)) if filled.any() else 0.0
    views = camera_view(cams, frame, Ec, Nc, np.where(filled, dem, z_edge), core,
                        z_edge) if cams else {}
    stretches, seam = camera_stretches(views, uc)

    # ---- heading: measured first, so the page grows to fit it -------------
    W = 14.0
    left, map_w, map_h = 0.95, 11.45, 2.75
    l1, facts, left_out = title_lines(info, cells, station, cutoff, min_points, cell)
    measure = RendererAgg(int(W * dpi), 200, dpi)

    def width_px(text, size, weight="normal"):
        return measure.get_text_width_height_descent(
            text, FontProperties(size=size, weight=weight), ismath=False)[0]

    max_px = (W - left - 0.35) * dpi
    fact_lines = wrap_facts(facts, max_px, lambda t: width_px(t, 11))
    out_lines = wrap_facts(left_out, max_px, lambda t: width_px(t, 10)) if left_out else []
    y_text = [0.42]
    for _ in fact_lines:
        y_text.append(y_text[-1] + (0.30 if len(y_text) == 1 else 0.23))
    for _ in out_lines:
        y_text.append(y_text[-1] + (0.25 if len(y_text) == 1 + len(fact_lines) else 0.21))

    # ---- how to read it: composed here, so the page is tall enough for it -----
    ulim, vlim, ex = map_limits(uc[core], vc[core], map_w, map_h)
    xlim = (ulim[1], ulim[0]) if frame["flip_x"] else ulim
    stretch_txt = (f"stretched ×{ex:g}, so slopes look {ex:g} times steeper there than "
                   "they are" if ex > 1 else "drawn to scale")
    cover = ""
    if len(stretches) == 2:
        (ca, a0, a1), (cb_, b0, b1) = stretches
        cover = (f" {ca} alone measured the {a1 - a0:.0f} m {compass(-frame['a'])} of the "
                 f"seam, {cb_} the {b1 - b0:.0f} m {compass(frame['a'])} of it; the seam is "
                 f"where their photos meet.")
    elif len(stretches) == 1:
        cover = f" All of it was measured by {stretches[0][0]}."
    share = float(info.get("setup_share") or 0.0)
    if share >= 0.99:
        level = ("the water level measured at that moment plus the wave setup (a timex "
                 "shows the edge where the waves run up to)")
    elif share > 0:
        level = (f"the water level measured at that moment, plus the wave setup for "
                 f"{100 * share:.0f}% of the points (the rest had no wave record)")
    else:
        level = "the water level measured at that moment"
    notes = [
        f"Each waterline is where the water met the sand, so the beach there is at {level}. "
        "A cell's elevation is the median of the frames that crossed it, one sample per "
        "frame; intertidal zone only.",
        "Spread is the 16–84 percentile range of those frames: the DEM's own repeatability, "
        "not its accuracy. Real change inside the window (a storm cutting the upper beach) "
        "shows up as spread too.",
        f"Maps are drawn as seen from the bluff, sea at the top; cross-shore is "
        f"{stretch_txt}.{cover}",
        f"Profiles: the cells within {cell:g} m of each position; the slope is the "
        "least-squares line through them.",
    ]
    w_c = 6.9
    x_notes = left + w_c + 0.6
    note_px = (W - x_notes - 0.25) * dpi
    note_paras = [wrap_words(t, note_px, lambda x: width_px(x, 10)) for t in notes]
    notes_h = sum(NOTE_LINE * len(w) + NOTE_GAP for w in note_paras)

    # ---- profiles: positions ---------------------------------------------
    if profiles is None:
        profiles = pick_profiles(uc[core], vc[core], zc[core], cell, n_profiles, stretches)
    off_map = [p for p in profiles if not ulim[0] <= p <= ulim[1]]
    if off_map:
        print("NOTE: profile position(s) off the strip, not drawn: "
              + ", ".join(f"{p:g} m" for p in off_map))
    profiles = [p for p in profiles if ulim[0] <= p <= ulim[1]]
    if len(profiles) > MAX_PROFILES:
        print(f"NOTE: {len(profiles)} profiles asked for; the first {MAX_PROFILES} drawn "
              f"(the most whose colours stay distinct)")
    profs = []
    for i, p in enumerate(profiles[:MAX_PROFILES]):
        profs.append((p, profile_at(uc.ravel(), vc.ravel(), zc.ravel(), p, cell),
                      chr(ord("A") + i), PROFILE_COLOURS[i]))

    cam_names = sorted(views) or sorted(cams)
    cam_label = ("cameras " + " + ".join(cam_names)) if len(cam_names) > 1 else \
        ("camera " + "".join(cam_names) if cam_names else "camera")

    # ---- legend entries: in one row if they fit the page, else two -------------
    handles = [Patch(facecolor=NO_VALUE, edgecolor="none",
                     label=f"no value: crossed by fewer than {min_points} frames"),
               Patch(facecolor=BLANK_FACE, edgecolor=BLANK_HATCH, hatch="////", lw=0,
                     label=f"blanked: spread > {cutoff:g} m")]
    if interpolated.any():
        handles.append(Patch(facecolor=SAND[2], edgecolor=INK2, hatch="..", lw=0,
                             label="interpolated, not measured"))
    if views:
        handles.append(Line2D([], [], color=INK2, lw=1.1, ls=(0, (5, 3)),
                              label="edge of a camera's photo"))
    if frame["from_camera"]:
        handles.append(Line2D([], [], color=INK, marker="^", ls="none", ms=9,
                              label=f"{cam_label} (alongshore 0)"))
    if profs:
        handles.append(Line2D([], [], color=INK, lw=0.9,
                              label=("profile A–" + profs[-1][2]) if len(profs) > 1
                              else "profile A"))
    em = 10 * dpi / 72
    leg_w = sum(width_px(h.get_label(), 10) + (2.0 + 0.6 + 1.5) * em for h in handles)
    leg_rows = 1 if leg_w <= max_px else 2
    leg_cols = int(np.ceil(len(handles) / leg_rows))

    # ---- layout, in inches from the top of the page ----------------------
    y_leg = y_text[-1] + 0.36              # legend row(s)
    cb_x, cb_w = left + map_w + 0.25, 0.17
    y_a = y_leg + 0.78 + 0.26 * (leg_rows - 1)   # elevation map
    y_cov = y_a + map_h                   # camera stretches, under it
    y_b = y_cov + 0.95                    # spread map
    y_c = y_b + map_h + 1.05              # profiles
    h_c = 2.15
    H = max(y_c + h_c + 0.62, y_c + notes_h + 0.25)

    def box(x, top, w, h):
        return [x / W, 1 - (top + h) / H, w / W, h / H]

    fig = plt.figure(figsize=(W, H), dpi=dpi, facecolor="white")

    fig.text(left / W, 1 - y_text[0] / H, l1, fontsize=16, fontweight="bold", color=INK)
    for i, t in enumerate(fact_lines):
        fig.text(left / W, 1 - y_text[1 + i] / H, t, fontsize=11, color=INK2)
    for i, t in enumerate(out_lines):
        fig.text(left / W, 1 - y_text[1 + len(fact_lines) + i] / H, t, fontsize=10, color=MUTED)

    # colour scales
    zmin, zmax = np.percentile(dem[filled], [1, 99]) if filled.any() else (0.0, 1.0)
    vmin, vmax = np.floor(zmin * 2) / 2, np.ceil(zmax * 2) / 2
    if vmax - vmin < 0.5:
        vmax = vmin + 0.5
    sand = LinearSegmentedColormap.from_list("sand_low_dark", SAND[::-1])
    znorm = Normalize(vmin, vmax)
    violet = LinearSegmentedColormap.from_list("violet", VIOLET)
    if hasattr(violet, "with_extremes"):                 # matplotlib >= 3.4
        violet = violet.with_extremes(over=VIOLET_OVER)
    else:
        violet.set_over(VIOLET_OVER)
    snorm = Normalize(0, cutoff)

    def quads(mask):
        r, c = np.nonzero(mask)
        return np.stack([np.stack([uk[r, c], vk[r, c]], -1),
                         np.stack([uk[r, c + 1], vk[r, c + 1]], -1),
                         np.stack([uk[r + 1, c + 1], vk[r + 1, c + 1]], -1),
                         np.stack([uk[r + 1, c], vk[r + 1, c]], -1)], 1)

    ylab = (f"cross-shore ×{ex:g}\n(m seaward of camera)" if frame["from_camera"]
            else f"cross-shore ×{ex:g}\n(m seaward)")

    def base_map(ax):
        ax.set_xlim(*xlim)
        ax.set_ylim(*vlim)
        ax.set_facecolor("white")
        ax.grid(True, color=GRID, lw=0.6, zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(length=3, labelsize=9.5)
        ax.set_ylabel(ylab, fontsize=10)
        if no_value.any():
            ax.add_collection(PolyCollection(quads(no_value), facecolors=NO_VALUE,
                                             edgecolors=NO_VALUE, linewidths=0.3, zorder=1))
        for vw in views.values():
            for eu, ev in vw["edges"].values():
                ax.plot(eu, ev, color=INK2, lw=1.1, ls=(0, (5, 3)), zorder=6)
        if frame["from_camera"] and ulim[0] <= 0 <= ulim[1]:
            ax.plot([0], [vlim[0]], marker="^", ms=10, color=INK, clip_on=False, zorder=9)

    def mark_profiles(ax, letters):
        for p, pr, letter, _ in profs:
            if not ulim[0] <= p <= ulim[1]:
                continue
            top = min(np.nanmax(pr["v"]) + 2, vlim[1]) if pr else vlim[1]
            bot = max(np.nanmin(pr["v"]) - 2, vlim[0]) if pr else vlim[0]
            ax.plot([p, p], [bot, top], color=INK, lw=0.9, zorder=7, solid_capstyle="butt",
                    path_effects=[patheffects.withStroke(linewidth=2.4, foreground="white")])
            if letters:
                ax.annotate(letter, (p, vlim[1]), xytext=(0, 4), textcoords="offset points",
                            ha="center", va="bottom", fontsize=11, fontweight="bold",
                            color=INK, annotation_clip=False)

    # ---- elevation map -----------------------------------------------------
    ax_a = fig.add_axes(box(left, y_a, map_w, map_h))
    base_map(ax_a)
    ax_a.tick_params(labelbottom=False)
    mesh = ax_a.pcolormesh(uk, vk, np.ma.masked_invalid(zc), cmap=sand, norm=znorm,
                           shading="flat", zorder=2, edgecolors="face", linewidth=0.15)
    if blanked.any():
        ax_a.add_collection(PolyCollection(quads(blanked), facecolors=BLANK_FACE,
                                           edgecolors=BLANK_HATCH, linewidths=0, hatch="////",
                                           zorder=3))
    if interpolated.any():
        ax_a.add_collection(PolyCollection(quads(interpolated), facecolors="none",
                                           edgecolors=INK2, linewidths=0, hatch="..", zorder=3))
    levels = np.arange(np.ceil(vmin * 2) / 2, vmax + 1e-9, 0.5)
    level_ink = {round(L, 2): (INK if _luminance(sand(znorm(L))) > 0.30 else "white")
                 for L in levels}
    # where a contour label would sit on a profile line or a photo edge
    prof_px = ax_a.transData.transform(np.column_stack([[p for p, _, _, _ in profs],
                                                        np.full(len(profs), vlim[0])]))[:, 0] \
        if profs else np.zeros(0)
    # photo edges as px segments (the samples along an image side land
    # sparsely on the ground, so test against the segments, not the points)
    edge_seg = []
    for vw in views.values():
        for eu, ev in vw["edges"].values():
            q = ax_a.transData.transform(np.column_stack([eu, ev]))
            good = np.isfinite(q).all(axis=1)
            both = good[:-1] & good[1:]
            edge_seg += [np.hstack([q[:-1][both], q[1:][both]])]
    edge_seg = np.vstack(edge_seg) if edge_seg else np.zeros((0, 4))
    ax_px = ax_a.get_window_extent()

    def label_clear(P, half):
        """A label of half-width `half` (px) centred at P (px) clears the frame,
        the profile lines and the photo edges."""
        if not (ax_px.x0 + half + 4 <= P[0] <= ax_px.x1 - half - 4
                and ax_px.y0 + 10 <= P[1] <= ax_px.y1 - 10):
            return False
        if prof_px.size and np.min(np.abs(prof_px - P[0])) < half + 8:
            return False
        if not edge_seg.size:
            return True
        a, b = edge_seg[:, :2], edge_seg[:, 2:]
        ab = b - a
        t = np.clip(((P - a) * ab).sum(1) / np.maximum((ab * ab).sum(1), 1e-9), 0, 1)
        return np.min(np.hypot(*(a + t[:, None] * ab - P).T)) >= half + 6

    # Contours, drawn by hand: fragments shorter than CONTOUR_MIN_M on the
    # page are left out (a 1-cell wiggle reads as a dash, and dashes mean
    # photo edges here); labels go on the longest pieces, at most three
    # per level, never on a profile line, a photo edge or another label.
    zs = smooth_for_contours(zc)
    drawn_levels = []
    if np.isfinite(zs).sum() > 20 and len(levels):
        to_px = ax_a.transData.transform
        to_data = ax_a.transData.inverted().transform
        placed = []                                   # label centres, px
        lab_h = 9.5 * dpi / 72
        for L, segs in contour_segments(uc, vc, zs, levels):
            colour = level_ink.get(round(L, 2), INK)
            text = minus(L)
            # the gap around a label; wider for a negative one, whose minus
            # sign would otherwise read as the end of the line it sits in
            half = width_px(text, 9.5, "bold") / 2 + (15 if L < 0 else 8)
            segs = [sg for sg in segs
                    if np.hypot(np.diff(sg[:, 0]), np.diff(sg[:, 1]) * ex).sum() >= CONTOUR_MIN_M]
            if not segs:
                continue
            drawn_levels.append(L)
            segs.sort(key=lambda sg: -np.hypot(np.diff(sg[:, 0]), np.diff(sg[:, 1]) * ex).sum())
            n_lab = 0
            for k, sg in enumerate(segs):
                P = to_px(sg)
                pieces = [P]
                length = np.hypot(np.diff(sg[:, 0]), np.diff(sg[:, 1]) * ex).sum()
                span = np.ptp(sg[:, 0])
                targets = []
                # long pieces get labels; a level broken into short pieces
                # still gets one, on its longest piece if that can hold it
                if n_lab < 3 and (length >= 40 or (k == 0 and length >= 20)):
                    targets = [0.28, 0.72] if (k == 0 and span > 200) else [0.5]
                for fr in targets:
                    if n_lab >= 3:
                        break
                    d = np.r_[0.0, np.cumsum(np.hypot(*np.diff(P, axis=0).T))]
                    ok = np.array([label_clear(q, half) for q in P])
                    ok &= (d - half >= 0) & (d + half <= d[-1])
                    for q in placed:
                        ok &= (np.abs(P[:, 0] - q[0]) > 2.5 * half) | \
                              (np.abs(P[:, 1] - q[1]) > 1.6 * lab_h)
                    if not ok.any():
                        continue
                    target = sg[:, 0].min() + fr * span
                    i = int(np.argmin(np.where(ok, np.abs(sg[:, 0] - target), np.inf)))
                    # cut the piece of line that holds this point
                    for j, piece in enumerate(pieces):
                        dj = np.r_[0.0, np.cumsum(np.hypot(*np.diff(piece, axis=0).T))]
                        hit = np.flatnonzero(np.hypot(*(piece - P[i]).T) < 0.5)
                        if not hit.size:
                            continue
                        cut = _cut_for_label(piece, dj[hit[0]], half)
                        if cut is None:
                            break
                        before, after, ang = cut
                        pieces[j:j + 1] = [before, after]
                        x, y = to_data(P[i])
                        ax_a.text(x, y, text, rotation=ang, rotation_mode="anchor",
                                  ha="center", va="center", fontsize=9.5, fontweight="bold",
                                  color=colour, zorder=5)
                        placed.append(P[i])
                        n_lab += 1
                        break
                for piece in pieces:
                    if len(piece) >= 2:
                        xy = to_data(piece)
                        ax_a.plot(xy[:, 0], xy[:, 1], color=colour, lw=0.9, zorder=4,
                                  solid_capstyle="round")
    mark_profiles(ax_a, letters=True)
    fig.text(left / W, 1 - (y_a - 0.3) / H, "Elevation (m NAVD88), contours every 0.5 m",
             fontsize=12, fontweight="bold", color=INK, va="bottom")
    fig.text((left + map_w) / W, 1 - (y_a - 0.3) / H, "light = high, dark = low",
             fontsize=10, color=INK2, va="bottom", ha="right")

    cax = fig.add_axes(box(cb_x, y_a, cb_w, map_h))
    ticks = np.arange(vmin, vmax + 1e-9, 0.5)
    cb = fig.colorbar(mesh, cax=cax, ticks=ticks)
    cb.ax.set_yticklabels([minus(t) for t in ticks])
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=9.5, length=2)
    cb.set_label("elevation (m NAVD88)", fontsize=10, color=INK2)
    for L in drawn_levels:                       # the contour levels, on the scale
        f = (L - vmin) / (vmax - vmin)
        cb.ax.plot([0, 1], [f, f], transform=cb.ax.transAxes,
                   color=level_ink.get(round(L, 2), INK), lw=1.0)

    # north arrow, beside the elevation map's title. On the page north is
    # its alongshore and cross-shore components times each axis scale
    # (inches per metre; cross-shore stretched).
    nvec = np.array([0.0, 1.0])
    dx = (nvec @ frame["a"]) * map_w / (ulim[1] - ulim[0]) * (-1 if frame["flip_x"] else 1)
    dy = (nvec @ frame["s"]) * map_h / (vlim[1] - vlim[0])
    nrm = np.hypot(dx, dy)
    ax_n = fig.add_axes(box(cb_x - 0.15, y_a - 0.62, 0.55, 0.55))
    ax_n.set_xlim(-1, 1)
    ax_n.set_ylim(-1, 1)
    ax_n.axis("off")
    ax_n.annotate("", xy=(0.55 * dx / nrm, 0.55 * dy / nrm),
                  xytext=(-0.55 * dx / nrm, -0.55 * dy / nrm),
                  arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.3, mutation_scale=12))
    ax_n.text(0.9 * dx / nrm, 0.9 * dy / nrm, "N", ha="center", va="center", fontsize=10,
              fontweight="bold", color=INK)

    # ---- which camera sees which stretch, under the elevation map ----------
    ax_cov = fig.add_axes(box(left, y_cov + 0.12, map_w, 0.42), sharex=ax_a)
    ax_cov.set_ylim(0, 1)
    ax_cov.axis("off")
    # One wording for every camera: the longest that leaves each bar's
    # ends and the seam tick showing.
    px_per_m = map_w * dpi / (ulim[1] - ulim[0])

    def cov_labels(form):
        out = []
        for cam, _, _ in stretches:
            n = framed.get(cam)
            if "{n}" in form and not n:          # frames not recorded (no info file)
                form = form.replace(" · {n} frames", "")
            out.append(form.format(cam=cam, n=f"{n:,}" if n else ""))
        return out

    forms = ["seen by {cam} · {n} frames", "{cam} · {n} frames", "seen by {cam}", "{cam}"]
    room = [(min(b, ulim[1]) - max(a, ulim[0])) * px_per_m - 0.3 * dpi for _, a, b in stretches]
    labs = next((ls for ls in map(cov_labels, forms)
                 if all(width_px(f" {t} ", 10.5) <= r for t, r in zip(ls, room))),
                cov_labels("{cam}"))
    for (cam, a, b), lab in zip(stretches, labs):
        y = 0.45
        ax_cov.plot([a, b], [y, y], color=INK2, lw=1.2, solid_capstyle="butt")
        for x in (a, b):
            ax_cov.plot([x, x], [y - 0.22, y + 0.22], color=INK2, lw=1.2)
        ax_cov.text((max(a, ulim[0]) + min(b, ulim[1])) / 2, y, f" {lab} ", ha="center",
                    va="center", fontsize=10.5, color=INK,
                    bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none"))
    if seam is not None:
        ax_cov.annotate(f"seam, {seam:.0f} m", (seam, 0.45), xytext=(0, -12),
                        textcoords="offset points", ha="center", va="top", fontsize=9.5,
                        color=INK2)
    if not views:
        ax_cov.text(0.5, 0.45, "cameras not drawn: calibration not found",
                    transform=ax_cov.transAxes, ha="center", va="center", fontsize=9.5,
                    color=MUTED)

    # ---- spread map -----------------------------------------------------------
    ax_b = fig.add_axes(box(left, y_b, map_w, map_h), sharex=ax_a)
    base_map(ax_b)
    smesh = ax_b.pcolormesh(uk, vk, np.ma.masked_invalid(np.where(core, spread, np.nan)),
                            cmap=violet, norm=snorm, shading="flat", zorder=2,
                            edgecolors="face", linewidth=0.15)
    if blanked.any():
        ax_b.add_collection(PolyCollection(quads(blanked), facecolors="none", edgecolors="white",
                                           linewidths=0, hatch="////", zorder=3))
    mark_profiles(ax_b, letters=False)
    ax_b.set_xlabel(("alongshore distance from the camera (m), increasing towards "
                     if frame["from_camera"] else "alongshore distance (m), increasing towards ")
                    + compass(frame["a"]), fontsize=10)
    fig.text(left / W, 1 - (y_b - 0.12) / H,
             "Repeatability: 16–84 percentile spread of the frames in each cell (m)",
             fontsize=12, fontweight="bold", color=INK, va="bottom")
    fig.text((left + map_w) / W, 1 - (y_b - 0.12) / H,
             f"cells above {cutoff:g} m (--max-spread) are blanked in the elevation map",
             fontsize=10, color=INK2, va="bottom", ha="right")
    cbx = fig.add_axes(box(cb_x, y_b, cb_w, map_h))
    scb = fig.colorbar(smesh, cax=cbx, extend="max", ticks=np.linspace(0, cutoff, 6))
    scb.ax.set_yticklabels([f"{t:.2g}" for t in np.linspace(0, cutoff, 6)])
    scb.outline.set_visible(False)
    scb.ax.tick_params(labelsize=9.5, length=2)
    scb.set_label("spread (m)", fontsize=10, color=INK2)
    fig.text((cb_x - 0.05) / W, 1 - (y_b - 0.1) / H, f"> {cutoff:g}: blanked",
             fontsize=9.5, color=INK, va="bottom", ha="left")

    # ---- legend row ---------------------------------------------------------------
    fig.legend(handles=handles, loc="center left",
               bbox_to_anchor=(left / W, 1 - (y_leg + 0.13 * (leg_rows - 1)) / H),
               ncol=leg_cols, frameon=False, fontsize=10, handlelength=2.0,
               columnspacing=1.5, handletextpad=0.6, borderaxespad=0)

    # ---- profiles panel -------------------------------------------------------------
    ax_c = fig.add_axes(box(left, y_c, w_c, h_c))
    ax_c.grid(True, color=GRID, lw=0.6)
    ax_c.set_axisbelow(True)
    for sp in ("top", "right"):
        ax_c.spines[sp].set_visible(False)
    ax_c.tick_params(length=3, labelsize=9.5)
    drawn = [(p, pr, L, col) for p, pr, L, col in profs if pr is not None]
    for p, pr, letter, col in drawn:
        ax_c.plot(pr["v"], pr["z"], "o", ms=3, color=col, alpha=0.35, mec="none", zorder=2)
        slope = (f"tanβ {pr['tanb']:.3f} (1:{1 / pr['tanb']:.0f})" if pr["tanb"] > 0
                 else "no seaward slope")
        ax_c.plot(pr["bv"], pr["bz"], "-", color=col, lw=2, zorder=3,
                  label=f"{letter}  {p:.0f} m:  {slope}")
    if drawn:
        ax_c.legend(loc="best", frameon=False, fontsize=9.5, handlelength=1.6,
                    title="alongshore position: foreshore slope", title_fontsize=9.5)
        xs = np.concatenate([pr["v"] for _, pr, _, _ in drawn])
        zz = np.concatenate([pr["z"] for _, pr, _, _ in drawn])
        ax_c.set_xlim(xs.min() - 3, xs.max() + 8)
        ax_c.set_ylim(zz.min() - 0.15, zz.max() + 0.15)
        # Letters just beyond the seaward ends; where the lines converge
        # they are moved apart, with a leader back to their line.
        ends = np.array([[pr["bv"][-1], pr["bz"][-1]] for _, pr, _, _ in drawn])
        ends_px = ax_c.transData.transform(ends)
        bb = ax_c.get_window_extent()
        lab_px = place_end_labels(ends_px, (11 * dpi / 72, 12 * dpi / 72), (bb.y0, bb.y1),
                                  dx_px=9 * dpi / 72)
        for (_, _, letter, _), end, e, lp in zip(drawn, ends, ends_px, lab_px):
            off = (lp - e) * 72 / dpi
            arrow = (dict(arrowstyle="-", color=INK2, lw=0.6, shrinkA=6, shrinkB=2)
                     if abs(off[1]) > 4 or off[0] > 12 else None)
            ax_c.annotate(letter, tuple(end), xytext=tuple(off), textcoords="offset points",
                          ha="center", va="center", fontsize=10.5, fontweight="bold",
                          color=INK, arrowprops=arrow)
        h_m = (zz.max() - zz.min() + 0.3) / h_c          # metres per inch, vertical
        w_m = (xs.max() - xs.min() + 11) / w_c            # metres per inch, horizontal
        title_c = (f"Cross-shore profiles at {drawn[0][2]}–{drawn[-1][2]}"
                   if len(drawn) > 1 else f"Cross-shore profile at {drawn[0][2]}")
        title_c += f" (vertical ×{w_m / h_m:.0f})" if w_m / h_m >= 1.5 else ""
    else:
        ax_c.text(0.5, 0.5, "no profile with enough cells", transform=ax_c.transAxes,
                  ha="center", va="center", color=MUTED)
        title_c = "Cross-shore profiles"
    ax_c.set_xlabel("cross-shore distance seaward of the camera (m)" if frame["from_camera"]
                    else "cross-shore distance (m)", fontsize=10)
    ax_c.set_ylabel("elevation (m NAVD88)", fontsize=10)
    fig.text(left / W, 1 - (y_c - 0.12) / H, title_c, fontsize=12, fontweight="bold",
             color=INK, va="bottom")

    # ---- how to read it -----------------------------------------------------------------
    fig.text(x_notes / W, 1 - (y_c - 0.12) / H, "How to read this page", fontsize=12,
             fontweight="bold", color=INK, va="bottom")
    yy = y_c + 0.02
    for wrapped in note_paras:
        fig.text(x_notes / W, 1 - yy / H, "\n".join(wrapped), fontsize=10, color=INK2,
                 va="top", linespacing=1.35)
        yy += NOTE_LINE * len(wrapped) + NOTE_GAP

    out = Path(output) if output else Path(stem + "_dem.png")
    fig.savefig(out, dpi=dpi, facecolor="white")
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser(description="Draw the DEM page from existing grids.")
    ap.add_argument("stem", help="DEM stem, e.g. dem_intertidal_7day for dem_intertidal_7day_dem.asc")
    ap.add_argument("--output", default=None, help="PNG to write (default <stem>_dem.png)")
    ap.add_argument("--profiles", default=None,
                    help="Alongshore positions of the profiles, m from the camera, "
                         "comma-separated (default: picked automatically)")
    ap.add_argument("--n-profiles", type=int, default=4,
                    help=f"Profiles to pick (3-{MAX_PROFILES}, default 4)")
    ap.add_argument("--calibration", default=str(HERE / "calibration"),
                    help="Folder with CACO05_<cam>_20240801_IO.yaml and _20251113_EO-CV.yaml")
    ap.add_argument("--station", default=STATION)
    ap.add_argument("--max-spread", type=float, default=None,
                    help="Cutoff the DEM was built with (default: from <stem>_info.json, else 0.5)")
    ap.add_argument("--min-points", type=int, default=None,
                    help="--min-points the DEM was built with (default: from the info file, else 3)")
    args = ap.parse_args()
    if not Path(args.stem + "_dem.asc").exists():
        sys.exit(f"not found: {args.stem}_dem.asc")
    prof = [float(x) for x in args.profiles.split(",")] if args.profiles else None
    out = draw_page(args.stem, output=args.output, calibration=args.calibration, profiles=prof,
                    n_profiles=max(3, min(MAX_PROFILES, args.n_profiles)), station=args.station,
                    max_spread=args.max_spread, min_points=args.min_points)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

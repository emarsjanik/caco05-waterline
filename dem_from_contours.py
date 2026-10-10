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
    amount that grows with wave height, but is assigned the still-water
    elevation -- so the DEM reads LOW where rough-water frames land.
    Not corrected here; --max-hs leaves those frames out instead, using
    the offshore wave record (fetch_buoy_waves.py).

Usage:
    python3 dem_from_contours.py contour_points_ground.csv dem_out \\
        [--cell 2.0] [--min-points 3] [--camera c1|c2|both]
        [--max-spread 0.5] [--fill-gaps] [--start-date ...] [--end-date ...]

Writes dem_out_{dem,spread,count}.asc (and _source.asc with
--interpolate-edge), dem_out_info.json (dates, frames per camera and
what the filters left out -- what the page needs and the grids do not
hold) and the page dem_out_dem.png, drawn by dem_figure.py. To redraw
the page from those files without rebuilding: python3 dem_figure.py dem_out

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
            with_setup += bool(beach)
            Z.append(float(beach or r["tide_elevation_navd88"]))
            cams.append(r["camera"])
            dates.append(day)
            frames.append(r["source_file"] if has_source else f"__point{len(frames)}")
    if exclude:
        print(f"Excluded by name  : {len(excluded)} frame(s) (--exclude)")
    if max_hs is not None:
        print(f"Wave filter       : {len(rough_frames)} frame(s) left out, offshore Hs > "
              f"{max_hs} m; {len(unknown_hs_frames)} frame(s) with no wave record kept")
    if stats is not None:
        stats.update(excluded_frames=len(excluded), rough_frames=len(rough_frames),
                     setup_share=with_setup / len(Z) if Z else 0.0)
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
                       "tp": num(r.get("offshore_tp_s"))}
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


def fit_setup_coefficient(pair_cell, pair_frame, info, min_frames):
    """
    Finds C in beach elevation = water level + C*sqrt(Hs*L0) that makes
    repeat crossings of each cell agree best (smallest median spread).
    Only frames with both Hs and Tp take part, at every C, so the
    comparison is like for like.
    """
    tide = np.array([info[f]["tide"] for f in pair_frame])
    hs = np.array([info[f]["hs"] for f in pair_frame])
    tp = np.array([info[f]["tp"] for f in pair_frame])
    phi = np.sqrt(hs * 9.81 * tp ** 2 / (2 * np.pi))
    ok = np.isfinite(phi) & np.isfinite(tide)
    groups = cell_groups(pair_cell[ok])
    coefs = np.round(np.arange(-0.04, 0.1201, 0.002), 3)
    result = [(c,) + median_spread_by_cell(groups, tide[ok] + c * phi[ok], min_frames)
              for c in coefs]
    return result, int(ok.sum()), len(set(pair_frame[ok]))


def day_offsets(pair_cell, pair_vals, pair_key, min_frames, min_pairs=30):
    """
    For each (camera, local day): the median difference between that
    day's samples and the DEM built WITHOUT that day, cell by cell.
    Leaving the day out matters: a day with many frames dominates the
    cells it covers and would otherwise hide its own error.
    """
    out = {}
    for k in sorted(set(pair_key)):
        mine = pair_key == k
        groups = cell_groups(pair_cell[~mine])
        order, c, starts, ends = groups
        others = pair_vals[~mine][order]
        med = {}
        for a, b in zip(starts, ends):
            if b - a >= min_frames:
                med[int(c[a])] = float(np.median(others[a:b]))
        res = [v - med[int(cc)] for cc, v in zip(pair_cell[mine], pair_vals[mine]) if int(cc) in med]
        out[k] = (float(np.median(res)) if len(res) >= min_pairs else np.nan, len(res))
    return out


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
    ap.add_argument("--max-day-offset", type=float, default=None,
                    help="Leave out camera-days whose samples sit, on median, more than this "
                         "many metres from the DEM built without that day (e.g. 0.15). "
                         "Catches days whose waterlines are consistently displaced -- a "
                         "different failure from the tide-direction check, which only sees "
                         "direction. Days with fewer than 30 comparable samples are kept, and "
                         "with fewer than 3 testable camera-days every day is kept (the DEM "
                         "without one day is then the other day alone, so the test cannot tell "
                         "which is off).")
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
                         "Run it on contours WITHOUT a correction applied.")
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
        result, n_pairs, n_frames = fit_setup_coefficient(pair_cell, pair_frame, info,
                                                          args.min_points)
        valid = [r for r in result if np.isfinite(r[1])]
        if not valid:
            print("Setup fit         : no frames with both Hs and Tp -- run "
                  "extract_elevation_contours.py with --waves.")
        else:
            base = next(r for r in valid if r[0] == 0.0)
            best = min(valid, key=lambda r: r[1])
            print(f"Setup fit         : {n_frames} frame(s), {n_pairs} cell samples")
            print(f"   C = 0      median spread {base[1]:.3f} m ({base[2]} cells)")
            print(f"   C = {best[0]:<6} median spread {best[1]:.3f} m  <-- best"
                  f"  (implied slope C/0.35 = {best[0] / 0.35:.3f})")
            for c, sp, n in valid:
                if round(c * 1000) % 20 == 0:
                    print(f"      C {c:+.2f}: {sp:.3f} m")
            if best[0] in (valid[0][0], valid[-1][0]):
                print("   WARNING: best C is at the edge of the search range -- the data do "
                      "not constrain it; do not apply.")
            elif best[0] <= 0:
                print("   The data do not favour a positive setup correction.")
            else:
                print(f"   To apply: extract_elevation_contours.py ... --waves ... "
                      f"--setup-coef {best[0]}")
        print()

    if args.max_day_offset:
        epoch = {f: info[f]["epoch"] for f in info}
        cam_of = dict(zip(frames, cams))
        pair_key = np.array([f"{cam_of[f]} {local_day(epoch[f])}" for f in pair_frame])
        z_of = {}
        for f, z in zip(frames, Z):
            z_of.setdefault(f, z)
        pair_vals = np.array([z_of[f] for f in pair_frame])
        offsets = day_offsets(pair_cell, pair_vals, pair_key, args.min_points)
        rejected = set()
        print(f"Day consistency (median offset from the DEM without that day, limit "
              f"+/-{args.max_day_offset} m):")
        # With fewer than 3 testable camera-days "the DEM without that day" is the
        # other day alone: two days that differ by more than the limit get equal
        # and opposite offsets and BOTH would be rejected, leaving no points. The
        # test cannot tell which day is off, so every day is kept.
        testable = [k for k, (off, n) in offsets.items() if np.isfinite(off)]
        too_few_days = len(testable) < 3
        for k, (off, n) in offsets.items():
            if not np.isfinite(off):
                note = "  (too few comparable samples -- kept)"
            elif abs(off) > args.max_day_offset and too_few_days:
                note = "  (off by more than the limit -- kept: fewer than 3 testable days)"
            elif abs(off) > args.max_day_offset:
                rejected.add(k)
                note = "  <-- REJECTED"
            else:
                note = ""
            off_s = f"{off:+.3f} m" if np.isfinite(off) else "   --   "
            print(f"   {k}   {off_s}   ({n} samples){note}")
        if too_few_days:
            print(f"   only {len(testable)} testable camera-day(s): the test cannot tell which day is off; "
                  f"every day kept")
        rejected_days = set(rejected)
        if rejected:
            frame_key = {f: f"{cam_of[f]} {local_day(epoch[f])}" for f in set(frames)}
            keep = np.array([frame_key[f] not in rejected for f in frames])
            E, N, Z, cams, dates, frames = E[keep], N[keep], Z[keep], cams[keep], dates[keep], frames[keep]
            print(f"   left out {len(rejected)} camera-day(s), {int((~keep).sum())} point(s)")
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
            "max_frame_offset": args.max_frame_offset,
            "frames_rejected": len(rejected_frames),
            "max_day_offset": args.max_day_offset,
            "camera_days_rejected": sorted(rejected_days),
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
    print("levels observed. No data above or below, and none invented. The detected")
    print("edge sits above still water by the wave setup, which grows with wave height,")
    print("but is given the still-water elevation -- so the DEM reads LOW where rough-")
    print("water frames land. Not corrected here; --max-hs leaves those frames out.")
    if rc:
        sys.exit(rc)


if __name__ == "__main__":
    main()

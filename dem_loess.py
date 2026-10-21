#!/usr/bin/env python3
"""
Scale-Controlled Loess DEM With Its Prediction Error
-------------------------------------------------------
An OPT-IN second surface beside the binned DEM: the same waterline points
interpolated by quadratic loess with explicit smoothing scales, and the
error of that interpolation at every cell (Plant, Holland & Puleo 2002).
dem_from_contours.py --loess writes it; survey_products.py --dem-loess asks
for it. Off by default. The binned DEM (<stem>_dem.asc and its companions),
its page and every figure are the same with it as without; no figure draws
it.

WHY. The binned DEM takes the median of the frames that crossed each 2 m
cell, and its _sigma.asc is the repeatability of those crossings. Neither
says how well the surface is known BETWEEN crossings, nor which scales of
the beach it resolves: a cell rests on whichever lines happened to cross
it. Plant et al. (2002) showed for nearshore surveys that gridding without
control of scale lets features shorter than the sampling (cusps, a berm
edge) alias into the surface, and that a scale-controlled interpolator --
a local weighted least-squares fit whose weights fall to zero at chosen
distances -- both filters those scales out explicitly and gives an error
estimate at every output point, from the scatter of the data about the
fit and from how the data lie around the point.

HOW (loess_grid). At each cell centre of the DEM's own grid:
  * the beach axes are the principal axes of the points: u alongshore,
    v cross-shore (as fine_dem in dem_from_contours.py);
  * every point within the ellipse of half-widths L_along (u) and L_across
    (v) gets the loess tricube weight w = (1 - r^3)^3, r =
    sqrt((du/L_along)^2 + (dv/L_across)^2) < 1 (Cleveland 1979; Cleveland
    & Devlin 1988), and none beyond. L_along and L_across are the smoothing
    scales: long alongshore, where the beach changes slowly and each line
    runs; short cross-shore, so the profile keeps its shape;
  * a quadratic surface z = b0 + b1 du + b2 dv + b3 du^2 + b4 du dv + b5 dv^2
    is fitted by weighted least squares -- QUADRATIC loess (Schlax & Chelton
    1992), the smoother Plant et al. (2002) used: unlike a weighted mean or
    a plane it follows a curved profile (a berm, a terrace) without
    flattening it, so the smoothing removes short features, not the
    profile's shape. The estimate is b0, the fit at the cell centre.
The estimate is a weighted sum of the data, z_hat = sum_j a_j z_j, with a_j
from the fit alone (the weights and where the points lie, not their
elevations). That is what gives the error:
  * s^2, the variance of the data about the fit: the weighted mean square
    residual. It holds the measurement noise and whatever the quadratic at
    these scales cannot follow (features shorter than the scales);
  * nmse, the NORMALISED MEAN-SQUARE ERROR: the share of s^2 that reaches
    the estimate, sum a_j^2 for independent data (Plant et al. 2002). Near
    0 where many data surround the cell, 1 where it rests on one datum,
    above 1 where the fit is carried beyond its data;
  * err = sqrt(nmse x s^2), the expected error of the estimate, m
    (<stem>_loess_err.asc);
  * skill = 1 - (the weighted mean square residual) / (the weighted
    variance of the data in the window): the share of the local scatter of
    elevations the surface explains -- the usual skill of a fit; near 1 on
    the sloping foreshore, low where the window's elevations are mostly
    noise.

TWO CHOICES FOR WATERLINES, both for the same reason as build_grid's rule
of one sample per frame: the points of one waterline are not independent.
Every point of a frame has the same elevation (one water level, one setup
per frame), so a frame's error in level, setup or placement is shared by
all its points, and a far-field line puts dozens of them in one window.
  * nmse is summed per FRAME: nmse = sum over frames of (sum of that
    frame's a_j)^2 -- the variance an error common to each frame passes to
    the estimate. With one point per frame it is Plant's sum a_j^2; with
    many it does not shrink as the points multiply, only as the frames do.
    For an error that varies along a line (a detection that wanders) it is
    an upper bound.
  * s^2 is the weighted residual sum over what it is expected to be for
    a unit error common to each frame: s^2 = sum_j w_j r_j^2 /
    sum_f (W_f - b_f' A^-1 b_f), with W_f the frame's weight (sum of its
    w_j), b_f the sum of its w_j x_j (x_j the model's terms at point j)
    and A = X'WX. The fit spends its six terms on the frames, and the
    residuals are smaller by what it spends: six frames make a quadratic
    go through them exactly, and the plain mean square residual would give
    an error near 0. The n/(n - 6) of an unweighted fit, n the effective
    number of frames ((sum of W_f)^2 / sum of W_f^2), does not hold for
    tricube weights: it puts s ~2.5 x too high at n 6-8 and ~1.4 x at n
    8-10. A cell needs n > 6 (n > 3 for order 1). For an error that varies
    along a line it is an upper bound too.
Points are first thinned to one per frame per half-cell (their mean), as in
fine_dem, so the fit is not slowed by far-field lines and no stretch of one
line outweighs the others. A week of the station's lines (86 frames) takes
~1 s in the sandbox; 3,000 synthetic frames over 400 m of beach ~1 min and
~300 MB (57 s on Python 3.11, 76 s on the station's Python 3.8 / numpy 1.17).

FILLED ONLY WHERE THE DATA SURROUND THE CELL. A cell is left nodata when its
window holds fewer than min_frames frames or too few for the error (n <= 6),
points on one side only, cross-shore or alongshore (the fit would be carried
beyond the highest or lowest waterline, or past the end of the lines: the
DEM covers the intertidal zone only and invents nothing beyond it), nmse
above max_nmse (LOESS_MAX_NMSE: the estimate no better known than the mean
of two frames), or s above max_sigma (dem_from_contours.py passes
--max-spread / 2, as for fine_dem: repeat crossings that disagree by more
are not measuring one surface). The counts of each are in the log and
<stem>_info.json (loess).

WHAT THE ERROR IS AND IS NOT. err is the random error of the interpolation
from the scatter of the lines about the surface and the way they lie
around the cell. Like _sigma.asc, it leaves out every error the frames
share and that does not average down: the water-level datum (the GNSS-R
antenna height, +/-0.061 m), the setup coefficient C (~0.06 m of DEM per
0.01 of C), a georectification bias (times the beach slope), a camera or
detector bias (the c1 - c2 offset of the DEM's log). It is the error of the
surface SMOOTHED at the chosen scales: a feature shorter than the scales is
not in the surface and not in err either (it is in s, the residual scatter).
And it treats frames as independent, which frames of one day sharing an
error are not. dem_loess.py evaluate checks it against a survey.
  On fixed synthetic lines whose only error is a level error common to each
frame (sd 0.04 m, 20 realisations; test_dem_loess.py) err is the size of
the error: the cells' RMS error over their RMS err is 0.94 (median) and
the truth lies within err 70% of the time, within 2 err 92%. With a 0.3 m
detection scatter along each line as well, the truth lies within err in
59-76% of the cells and within 2 err in 97-100%. On the real lines of the
NUC's survey-date builds (copies of 15 Oct 2026, the dem step's settings),
against the lidar with the median offset removed, the scatter is several
times err: Mar 2025 (c2, 3-9 Mar) NMAD 0.137 m against a median err of
0.038 m (3.6 x; 47% within 2 err), Jan 2025 (c2 alone, 18-23 Jan) 0.123 m
against 0.024 m (5.2 x; 33%). The binned DEM's _sigma.asc on the same
builds: 6.2 x and 4.1 x. So err is a floor on the error, not the accuracy:
the scatter about the lidar also holds what the frames share and does not
average down (a day's setup, the pointing, a detection offset), the
beach's change through the window and the lidar's own error.

THE SCALES (LOESS_SCALES, default 20 m alongshore x 5 m cross-shore, the
half-widths of the window). What they pass along each axis (response(): the
amplitude this 2-D window keeps of a sine of that wavelength running along
the axis, densely sampled data): 90% of a wave 1.84 L long, half of one 1.1
L long, at most 8% of any shorter than 0.8 L -- with L = 5 m, 90% of a 9.2 m
wave and half of a 5.5 m one; with L = 20 m, 90% of 37 m and half of 22 m.
(A plane, order 1, needs a wave 4.8 L long to keep 90%: it flattens the
profile's curvature.) Features shorter than those -- beach cusps (tens of
metres apart) alongshore, steps of a few metres cross-shore -- are smoothed
away, and the binned DEM keeps them where its cells were crossed.
Cross-shore the scale is set by the spacing of the lines: on the live lines
of the 2026-09-29 survey-date build, frames of one camera 30 min apart lie a
median 1.8 m apart cross-shore (0.6-3.3 m, 10th-90th percentile; 75 pairs),
so a 5 m half-width holds several lines on each side. Alongshore each line is
continuous and the scale only sets how much is averaged. dem_loess.py
response L prints the numbers for other half-widths.

References:
  Plant, N.G., Holland, K.T. & Puleo, J.A. (2002). Analysis of the scale of
    errors in nearshore bathymetric data. Marine Geology 191, 71-86.
  Schlax, M.G. & Chelton, D.B. (1992). Frequency domain diagnostics for
    linear smoothers. J. Amer. Statist. Assoc. 87, 1070-1081.
  Cleveland, W.S. (1979). Robust locally weighted regression and smoothing
    scatterplots. J. Amer. Statist. Assoc. 74, 829-836.
  Cleveland, W.S. & Devlin, S.J. (1988). Locally weighted regression: an
    approach to regression analysis by local fitting. J. Amer. Statist.
    Assoc. 83, 596-610.

Usage:
    python3 dem_from_contours.py contour_points_ground.csv dem_out --loess [--loess-scales 20,5]
        (writes dem_out_loess.asc, _loess_err.asc, _loess_nmse.asc, _loess_skill.asc)
    python3 dem_loess.py evaluate dem_out survey.tif|survey.asc|points.csv [--json out.json]
        (the binned DEM and the loess against a survey; how often the survey
        lies within err and 2 err of the loess: the NUC check of the error)
    python3 dem_loess.py response 20,5
"""

import sys
import json
import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Half-widths of the loess window, metres (alongshore, cross-shore): dem_from_contours.py --loess-scales.
LOESS_SCALES = (20.0, 5.0)
# A cell whose estimate keeps more than this share of the scatter about the fit is left nodata: 0.5 is the
# error of the mean of two independent frames.
LOESS_MAX_NMSE = 0.5
# Terms of the local model: 1 = plane (z, du, dv), 2 = quadratic (+ du^2, du dv, dv^2).
LOESS_ORDER = 2
# Grids written beside the DEM: <stem><suffix>.asc
LOESS_FILES = {"z": "_loess", "err": "_loess_err", "nmse": "_loess_nmse", "skill": "_loess_skill"}
# Most (nodes x points) solved at once: the arrays take ~120 bytes per pair (~60 MB here).
_BLOCK = 500000


def tricube(r):
    """The loess weight: (1 - r^3)^3 inside r < 1, 0 outside (Cleveland 1979)."""
    r = np.asarray(r, float)
    return np.where(r < 1, (1 - np.minimum(r, 1) ** 3) ** 3, 0.0)


def n_terms(order):
    return 3 if order == 1 else 6


def design(du, dv, order):
    """The local model's columns at offsets (du, dv) in units of the scales: 1, du, dv (+ du^2, du dv, dv^2)."""
    cols = [np.ones_like(du), du, dv]
    if order == 2:
        cols += [du * du, du * dv, dv * dv]
    return np.stack(cols, -1)


def beach_axes(E, N):
    """(centre E, centre N, alongshore unit vector, cross-shore unit vector): the principal axes of the points,
    as fine_dem."""
    ec, nc = float(np.mean(E)), float(np.mean(N))
    _, vecs = np.linalg.eigh(np.cov(np.stack([E - ec, N - nc])))
    return ec, nc, vecs[:, 1], vecs[:, 0]


def thin(E, N, Z, frames, h):
    """One point per frame per h x h bin, at the mean position and elevation of its points (as fine_dem).
    -> (E, N, Z, frame index, frame names)"""
    names, fid = np.unique(frames, return_inverse=True)
    fid = np.asarray(fid).ravel()
    key = np.stack([fid, np.floor(E / h).astype(np.int64), np.floor(N / h).astype(np.int64)], 1)
    _, inv = np.unique(key, axis=0, return_inverse=True)
    inv = np.asarray(inv).ravel()
    cnt = np.bincount(inv)
    first = np.zeros(len(cnt), np.int64)
    first[inv[::-1]] = np.arange(len(inv))[::-1]               # each bin's first point: its frame
    return (np.bincount(inv, E) / cnt, np.bincount(inv, N) / cnt, np.bincount(inv, Z) / cnt,
            fid[first], names)


def frame_sums(m, fid):
    """The columns of m (nodes x points) summed over the points of each frame: nodes x (frames present), in the
    order of the sorted frame indices. A sum over sorted runs, so the memory is nodes x frames, never points x
    frames (a whole-archive window holds thousands of frames)."""
    order = np.argsort(fid, kind="mergesort")
    f = fid[order]
    starts = np.flatnonzero(np.r_[True, f[1:] != f[:-1]])
    return np.add.reduceat(m[:, order], starts, axis=1)


def _solve_block(du, dv, w, z, fid, order):
    """The local fits of one block of nodes (rows) against their candidate points (columns), offsets in units of
    the scales. -> dict of per-node arrays (NaN where the fit is degenerate)."""
    p = n_terms(order)
    n = du.shape[0]
    out = {k: np.full(n, np.nan) for k in ("z", "nmse", "s2", "skill", "neff")}
    # frames in each window, their weight and (below) the share of the estimate they carry
    Wf = frame_sums(w, fid)                                     # nodes x frames
    out["frames"] = (Wf > 0).sum(1).astype(float)
    sw = Wf.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        neff = sw ** 2 / (Wf ** 2).sum(1)
    out["neff"] = neff
    X = design(du, dv, order)                                   # nodes x points x terms
    XW = X * w[..., None]
    A = np.einsum("npi,npj->nij", XW, X)
    # a degenerate design (points on one line, too few frames) has a near-zero eigenvalue
    ev = np.linalg.eigvalsh(A)
    ok = (sw > 0) & (ev[:, 0] > 1e-9 * np.maximum(ev[:, -1], 1e-300)) & (out["frames"] >= p)
    if not ok.any():
        return out
    Xo, XWo, Ao, wo = X[ok], XW[ok], A[ok], w[ok]
    e1 = np.zeros((int(ok.sum()), p, 1))
    e1[:, 0, 0] = 1.0
    g = np.linalg.solve(Ao, e1)[..., 0]                         # row 0 of A^-1 (A is symmetric)
    a = wo * np.einsum("npi,ni->np", Xo, g)                     # z_hat = sum_j a_j z_j
    beta = np.linalg.solve(Ao, np.einsum("npi,p->ni", XWo, z)[..., None])[..., 0]
    res = z[None, :] - np.einsum("npi,ni->np", Xo, beta)
    swo = wo.sum(1)
    wrss = (wo * res ** 2).sum(1)
    msr = wrss / swo
    zbar = (wo * z[None, :]).sum(1) / swo
    var = (wo * (z[None, :] - zbar[:, None]) ** 2).sum(1) / swo
    ne = neff[ok]
    # what the weighted residual sum expects per unit variance of an error common to each frame:
    # sum_f (W_f - b_f' A^-1 b_f), b_f = sum of the frame's w_j x_j
    Bf = np.stack([frame_sums(XWo[..., i], fid) for i in range(p)], 1)      # nodes x terms x frames
    dof = swo - np.einsum("nif,nif->n", Bf, np.linalg.solve(Ao, Bf))
    with np.errstate(invalid="ignore", divide="ignore"):
        s2 = np.where((ne > p) & (dof > 1e-9 * swo), wrss / dof, np.nan)
        skill = np.where(var > 0, 1.0 - msr / var, np.nan)
    out["z"][ok] = beta[:, 0]
    out["nmse"][ok] = (frame_sums(a, fid) ** 2).sum(1)         # per frame: (sum of its a_j)^2
    out["s2"][ok] = s2
    out["skill"][ok] = skill
    return out


def loess_grid(E, N, Z, frames, e0, n0, ncols, nrows, cell, along=LOESS_SCALES[0], across=LOESS_SCALES[1],
               order=LOESS_ORDER, min_frames=3, max_nmse=LOESS_MAX_NMSE, max_sigma=None, extrapolate=False,
               thin_bin=None):
    """
    The quadratic loess surface and its error on the grid of cells e0 + (col + 0.5) x cell, n0 + (row + 0.5) x
    cell (rows south to north, as build_grid), from the points (E, N, Z) of `frames`. See the module docstring.
    along, across: the half-widths of the window, m (alongshore, cross-shore). order: 2 quadratic, 1 plane.
    min_frames: fewest frames in a window. max_nmse (None: no limit), max_sigma (m; None: no limit): above
    them a cell is left nodata. extrapolate: True also fills cells with points on one side only (for tests of
    how the error grows away from the data; never for a product). thin_bin: m (default cell / 2).
    -> dict: z, err, nmse, sigma (s, m), skill, frames (in the window), neff (effective frames): (nrows,
    ncols) arrays, NaN where not filled; 'raw': the same before the limits (for diagnostics); counts: cells
    with data in reach, filled, and left out by each rule; and the settings.
    """
    E, N, Z = (np.asarray(x, float) for x in (E, N, Z))
    frames = np.asarray(frames)
    h = thin_bin or cell / 2.0
    shape = (nrows, ncols)
    size = nrows * ncols
    raw = {k: np.full(size, np.nan) for k in ("z", "nmse", "s2", "skill", "neff")}
    nfr = np.zeros(size)
    reach = np.zeros(size, bool)                  # a point within the window
    sides = np.zeros(size, bool)                  # points on both sides, cross-shore and alongshore
    names, ua = [], np.array([0.0, 1.0])
    if len(Z) >= n_terms(order):
        E, N, Z, fid, names = thin(E, N, Z, frames, h)
    if len(names) and len(Z) >= n_terms(order):
        ec, nc, ua, va = beach_axes(E, N)
        u = (E - ec) * ua[0] + (N - nc) * ua[1]
        v = (E - ec) * va[0] + (N - nc) * va[1]
        gx = e0 + (np.arange(ncols) + 0.5) * cell
        gy = n0 + (np.arange(nrows) + 0.5) * cell
        GX, GY = np.meshgrid(gx, gy)
        gu = ((GX - ec) * ua[0] + (GY - nc) * ua[1]).ravel()
        gv = ((GX - ec) * va[0] + (GY - nc) * va[1]).ravel()

        # buckets one window wide: a node's points all lie in its own and the 8 buckets around it
        bu, bv = np.floor(u / along).astype(np.int64), np.floor(v / across).astype(np.int64)
        order_p = np.lexsort((bv, bu))
        keys = np.stack([bu[order_p], bv[order_p]], 1)
        edges = np.flatnonzero(np.any(np.diff(keys, axis=0) != 0, axis=1)) + 1
        pts = {(int(keys[a, 0]), int(keys[a, 1])): order_p[a:b]
               for a, b in zip(np.r_[0, edges], np.r_[edges, len(order_p)])}
        nbu, nbv = np.floor(gu / along).astype(np.int64), np.floor(gv / across).astype(np.int64)
        norder = np.lexsort((nbv, nbu))
        nk = np.stack([nbu[norder], nbv[norder]], 1)
        nedges = np.flatnonzero(np.any(np.diff(nk, axis=0) != 0, axis=1)) + 1
        for a0, b0 in zip(np.r_[0, nedges], np.r_[nedges, len(norder)]):
            ku, kv = int(nk[a0, 0]), int(nk[a0, 1])
            cand = [pts[(ku + i, kv + j)] for i in (-1, 0, 1) for j in (-1, 0, 1) if (ku + i, kv + j) in pts]
            if not cand:
                continue
            p = np.concatenate(cand)
            nodes_all = norder[a0:b0]
            step = max(1, _BLOCK // max(len(p), 1))
            for c0 in range(0, len(nodes_all), step):
                nodes = nodes_all[c0:c0 + step]
                du = (u[p][None, :] - gu[nodes][:, None]) / along
                dv = (v[p][None, :] - gv[nodes][:, None]) / across
                w = tricube(np.sqrt(du ** 2 + dv ** 2))
                inside = w > 0
                hit = inside.any(1)
                if not hit.any():
                    continue
                reach[nodes[hit]] = True
                sides[nodes] = ((np.where(inside, dv, np.inf).min(1) < 0) & (np.where(inside, dv, -np.inf).max(1) > 0)
                                & (np.where(inside, du, np.inf).min(1) < 0)
                                & (np.where(inside, du, -np.inf).max(1) > 0))
                cols = inside.any(0)             # only the points some node of the block uses
                sel = nodes[hit]
                r = _solve_block(du[hit][:, cols], dv[hit][:, cols], w[hit][:, cols], Z[p][cols], fid[p][cols],
                                 order)
                for k in raw:
                    raw[k][sel] = r[k]
                nfr[sel] = r["frames"]

    # each cell in reach is left out by the first rule it fails, in this order
    pterms = n_terms(order)
    with np.errstate(invalid="ignore"):
        sigma = np.sqrt(raw["s2"])
        err = np.sqrt(raw["nmse"] * raw["s2"])
        few = reach & ((nfr < max(min_frames, 1)) | ~(raw["neff"] > pterms) | ~np.isfinite(raw["z"])
                       | ~np.isfinite(err))
        one_side = reach & ~few & ~sides & (not extrapolate)
        rest = reach & ~few & ~one_side
        loose = rest & (max_nmse is not None) & ~(raw["nmse"] <= (max_nmse or 0))
        noisy = rest & ~loose & (max_sigma is not None) & ~(sigma <= (max_sigma or 0))
    keep = rest & ~loose & ~noisy

    def grid(x):
        return np.where(keep, x, np.nan).reshape(shape)
    return {
        "z": grid(raw["z"]), "err": grid(err), "nmse": grid(raw["nmse"]), "sigma": grid(sigma),
        "skill": grid(raw["skill"]), "frames": np.where(keep, nfr, 0.0).reshape(shape), "neff": grid(raw["neff"]),
        "raw": {"z": raw["z"].reshape(shape), "err": err.reshape(shape), "nmse": raw["nmse"].reshape(shape),
                "sigma": sigma.reshape(shape), "sides": sides.reshape(shape), "reach": reach.reshape(shape)},
        "counts": {"in_reach": int(reach.sum()), "filled": int(keep.sum()), "too_few_frames": int(few.sum()),
                   "one_sided": int(one_side.sum()), "nmse_over": int(loose.sum()), "sigma_over": int(noisy.sum())},
        "settings": {"along_m": float(along), "across_m": float(across), "order": int(order),
                     "min_frames": int(min_frames), "max_nmse": max_nmse, "max_sigma_m": max_sigma,
                     "thin_bin_m": float(h), "points_fitted": int(len(Z)), "frames": int(len(names)),
                     "extrapolate": bool(extrapolate),
                     "axis_alongshore_deg": round(float(np.degrees(np.arctan2(ua[0], ua[1]))) % 180.0, 1)},
    }


def response(wavelength, scale, order=LOESS_ORDER, n=201):
    """The amplitude the window of half-width `scale` (2-D tricube weights on the ellipse, local polynomial of
    `order` in both axes) keeps of a sine of `wavelength` that runs along one of its axes, on densely and evenly
    sampled data (n x n points over the window): its transfer function at the centre of the window, sum_j a_j
    cos(2 pi x_j / wavelength) (Schlax & Chelton 1992). The same along either axis, each with its own half-width.
    1 passes, 0 removes."""
    t = np.linspace(-1.0, 1.0, n)
    du, dv = (q.ravel() for q in np.meshgrid(t, t))
    w = tricube(np.sqrt(du ** 2 + dv ** 2))
    k = w > 0
    du, dv, w = du[k], dv[k], w[k]
    X = design(du, dv, order)
    A = (X * w[:, None]).T @ X
    a = w * (X @ np.linalg.solve(A, np.eye(X.shape[1])[:, 0]))
    return float(np.sum(a * np.cos(2 * np.pi * dv * scale / wavelength)))


def response_table(scale, order=LOESS_ORDER):
    """(the wavelength kept at 90%, the one kept at 50%, m; the largest share of its amplitude any wavelength
    under 0.8 x scale keeps, sign aside) for a half-width `scale`, m."""
    lam = np.linspace(0.1, 6.0, 591) * scale
    h = np.array([response(L, scale, order) for L in lam])

    def first_over(level):
        i = np.flatnonzero(h >= level)
        return float(lam[i[0]]) if len(i) else float("nan")
    return first_over(0.9), first_over(0.5), float(np.abs(h[lam <= 0.8 * scale]).max())


def summary_lines(lo):
    """What the DEM's log prints of a loess_grid result."""
    s, c = lo["settings"], lo["counts"]
    ok = np.isfinite(lo["z"])
    lines = [f"loess DEM         : {'quadratic' if s['order'] == 2 else 'linear'} loess (Plant et al. 2002), window "
             f"{s['along_m']:g} m alongshore x {s['across_m']:g} m cross-shore (half-widths; shore axis "
             f"{s['axis_alongshore_deg']:g} deg), {s['points_fitted']} thinned points of {s['frames']} frames",
             f"                    {c['filled']} cells filled of {c['in_reach']} with points in reach; left out: "
             f"{c['too_few_frames']} too few frames, {c['one_sided']} with points on one side only (no "
             f"extrapolation), {c['nmse_over']} nmse > {s['max_nmse']}, {c['sigma_over']} scatter > "
             f"{s['max_sigma_m'] if s['max_sigma_m'] is not None else 'no limit'} m"]
    if ok.any():
        lines.append(f"                    prediction error median {np.nanmedian(lo['err']):.3f} m, p90 "
                     f"{np.nanpercentile(lo['err'], 90):.3f} m; nmse median {np.nanmedian(lo['nmse']):.3f}; scatter "
                     f"about the fit median {np.nanmedian(lo['sigma']):.3f} m (random error of the smoothed "
                     f"surface only: not the datum, setup or camera bias)")
    return lines


def info_record(lo):
    """What <stem>_info.json keeps (loess): the settings, the counts and the medians, and what the files are."""
    ok = np.isfinite(lo["z"])

    def med(x):
        return round(float(np.nanmedian(x)), 4) if ok.any() else None

    def p90(x):
        return round(float(np.nanpercentile(x, 90)), 4) if ok.any() else None
    order = lo["settings"]["order"]
    return {"method": f"{'quadratic' if order == 2 else 'linear'} loess, tricube weights, explicit scales (Plant, "
                      f"Holland & Puleo 2002, Marine Geology 191); nmse summed per frame, scatter over its expected "
                      f"share per frame after the {n_terms(order)} terms (dem_loess.py)",
            "files": {k: v + ".asc" for k, v in LOESS_FILES.items()},
            "settings": dict(lo["settings"]), "cells": dict(lo["counts"]),
            "err_median_m": med(lo["err"]), "err_p90_m": p90(lo["err"]), "nmse_median": med(lo["nmse"]),
            "sigma_median_m": med(lo["sigma"]), "skill_median": med(lo["skill"]),
            "definition": "err = sqrt(nmse x s^2): the random error of the surface smoothed at these scales, from "
                          "the scatter of the lines about it; not the datum, setup coefficient, georectification "
                          "or camera bias, which every frame shares"}


def write_grids(stem, lo, e0, n0, cell):
    """<stem>_loess.asc, _loess_err.asc, _loess_nmse.asc, _loess_skill.asc on the DEM's grid. -> paths"""
    from dem_from_contours import write_ascii_grid
    out = []
    for k, suf in LOESS_FILES.items():
        p = f"{stem}{suf}.asc"
        write_ascii_grid(p, lo[k], e0, n0, cell)
        out.append(p)
    return out


# ---------------------------------------------------------------------
# The NUC check: the loess and its error against a survey
# ---------------------------------------------------------------------

def _stats(d):
    d = d[np.isfinite(d)]
    if not len(d):
        return {"n": 0}
    m = float(np.median(d))
    return {"n": int(len(d)), "median_m": round(m, 4), "nmad_m": round(float(1.4826 * np.median(np.abs(d - m))), 4),
            "rms_m": round(float(np.sqrt(np.mean(d ** 2))), 4)}


def built_with_loess(stem):
    """The loess record of <stem>_info.json when THAT build wrote the loess grids (dem_from_contours.py --loess),
    else None: grids beside a DEM rebuilt without --loess are an earlier build's."""
    try:
        info = json.loads(Path(f"{stem}_info.json").read_text())
    except (OSError, ValueError):
        return None
    return info.get("loess") if isinstance(info, dict) else None


def evaluate(stem, survey):
    """
    The binned DEM and the loess of a build (<stem>_dem.asc, <stem>_loess*.asc) against a survey: a DSM grid
    (.tif / .asc: its median in each DEM cell, as survey_compare.py) or points (.csv: each point against the
    cell it lies in). Per surface, on the cells (points) where both surfaces have a value and on all of each:
    n, median, NMAD and RMS of surface - survey. For the loess, how often the survey lies within err and
    2 err of it -- 68% and 95% if err were the whole error and normal -- both as it is and after removing the
    median offset (err leaves out the errors every frame shares: datum, setup, camera bias), and by nmse.
    Cells are not independent (neighbours share frames): read the shares as a check of the size of err, not
    as a test at a stated confidence.
    -> (record, lines)
    """
    from asc_to_geotiff import read_asc
    dem, h = read_asc(f"{stem}_dem.asc")
    z, _ = read_asc(f"{stem}{LOESS_FILES['z']}.asc")
    err, _ = read_asc(f"{stem}{LOESS_FILES['err']}.asc")
    nmse, _ = read_asc(f"{stem}{LOESS_FILES['nmse']}.asc")
    cell, x0 = h["cellsize"], h["xllcorner"]
    ytop = h["yllcorner"] + dem.shape[0] * cell
    sp = Path(survey)
    if sp.suffix.lower() in (".tif", ".tiff", ".asc"):
        from compare_dem_survey import read_survey
        from survey_compare import survey_on_grid
        s, sx0, sy_top, scell = read_survey(sp)
        ref, _, how = survey_on_grid(s, sx0, sy_top, scell, x0, ytop, cell, dem.shape[0], dem.shape[1])
        pick = lambda g: g.ravel()                                          # noqa: E731
        ref = ref.ravel()
        unit = "cells"
    else:
        from survey_compare import read_points
        pts = read_points(sp)
        E_, N_, Z_ = (np.asarray(pts[k], float) for k in ("E", "N", "Z"))
        ci = np.floor((E_ - x0) / cell).astype(int)
        ri = np.floor((ytop - N_) / cell).astype(int)
        inside = (ci >= 0) & (ci < dem.shape[1]) & (ri >= 0) & (ri < dem.shape[0])
        ci, ri, ref = ci[inside], ri[inside], Z_[inside]
        pick = lambda g: g[ri, ci]                                          # noqa: E731
        how, unit = "each survey point against the DEM cell it lies in", "points"
    with np.errstate(invalid="ignore"):
        dz, lz, le, ln = pick(dem) - ref, pick(z) - ref, pick(err), pick(nmse)
    both = np.isfinite(dz) & np.isfinite(lz)
    rec = {"survey": str(sp), "how": how, "unit": unit,
           "dem_minus_survey": _stats(dz), "loess_minus_survey": _stats(lz),
           "on_shared": {"dem": _stats(dz[both]), "loess": _stats(lz[both])}}
    lines = [f"Loess check against {sp.name} ({how})",
             f"   binned DEM - survey : {_fmt(rec['dem_minus_survey'], unit)}",
             f"   loess - survey      : {_fmt(rec['loess_minus_survey'], unit)}",
             f"   on the {int(both.sum())} {unit} both have: DEM {_fmt(rec['on_shared']['dem'], unit, n=False)}; "
             f"loess {_fmt(rec['on_shared']['loess'], unit, n=False)}"]
    with np.errstate(invalid="ignore"):
        ok = np.isfinite(lz) & np.isfinite(le) & (le > 0)
    if ok.any():
        off = float(np.median(lz[ok]))
        cover = {}
        for name, d in (("as_is", lz[ok]), ("median_removed", lz[ok] - off)):
            cover[name] = {"within_err": round(float(np.mean(np.abs(d) <= le[ok])), 3),
                           "within_2err": round(float(np.mean(np.abs(d) <= 2 * le[ok])), 3)}
        rec["error_check"] = dict(cover, n=int(ok.sum()), median_offset_m=round(off, 4),
                                  err_median_m=round(float(np.median(le[ok])), 4),
                                  ratio_nmad_to_err=round(float(1.4826 * np.median(np.abs(lz[ok] - off))
                                                                / np.median(le[ok])), 2))
        lines.append(f"   survey within err / 2 err of the loess: {100 * cover['as_is']['within_err']:.0f}% / "
                     f"{100 * cover['as_is']['within_2err']:.0f}% as it is; {100 * cover['median_removed']['within_err']:.0f}"
                     f"% / {100 * cover['median_removed']['within_2err']:.0f}% with the median offset "
                     f"({off:+.3f} m) removed (68% / 95% if err were the whole random error)")
        lines.append(f"   NMAD of loess - survey / median err: {rec['error_check']['ratio_nmad_to_err']:.2f} "
                     f"(1 = err the size of the scatter; above 1: the surface is further off than err says)")
        bins = []
        for lo_, hi_ in ((0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 1e9)):
            with np.errstate(invalid="ignore"):
                m = ok & np.isfinite(ln) & (ln >= lo_) & (ln < hi_)
            if m.sum() >= 5:
                d = lz[m] - off
                bins.append({"nmse_from": lo_, "nmse_to": None if hi_ > 1e8 else hi_, "n": int(m.sum()),
                             "nmad_m": round(float(1.4826 * np.median(np.abs(d - np.median(d)))), 4),
                             "err_median_m": round(float(np.median(le[m])), 4),
                             "within_2err": round(float(np.mean(np.abs(d) <= 2 * le[m])), 3)})
        rec["by_nmse"] = bins
        for b in bins:
            hi_ = f"{b['nmse_to']:.1f}" if b["nmse_to"] is not None else "up"
            lines.append(f"      nmse {b['nmse_from']:.1f}-{hi_}: {b['n']} {unit}, NMAD {b['nmad_m']:.3f} m, "
                         f"median err {b['err_median_m']:.3f} m, within 2 err {100 * b['within_2err']:.0f}%")
    else:
        lines.append("   no loess value with an error where the survey is")
    return rec, lines


def _fmt(s, unit, n=True):
    if not s.get("n"):
        return f"no {unit}"
    return ((f"{s['n']} {unit}, " if n else "") + f"median {s['median_m']:+.3f} m, NMAD {s['nmad_m']:.3f} m, "
            f"RMS {s['rms_m']:.3f} m")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="what")
    ev = sub.add_parser("evaluate", help="a build's binned DEM and loess against a survey (DSM .tif/.asc or "
                                         "points .csv)")
    ev.add_argument("stem", help="the DEM's stem (<stem>_dem.asc, <stem>_loess.asc, ...: dem_from_contours.py --loess)")
    ev.add_argument("survey")
    ev.add_argument("--json", default=None, help="also write the record here")
    rs = sub.add_parser("response", help="what a window of these half-widths passes, by wavelength")
    rs.add_argument("scales", help="half-widths, m: 'alongshore,cross-shore' or one value")
    rs.add_argument("--order", type=int, default=LOESS_ORDER, choices=(1, 2))
    args = ap.parse_args()
    if args.what == "evaluate":
        if not Path(f"{args.stem}{LOESS_FILES['z']}.asc").exists():
            sys.exit(f"not found: {args.stem}{LOESS_FILES['z']}.asc (build with dem_from_contours.py ... --loess)")
        if built_with_loess(args.stem) is None:
            sys.exit(f"{args.stem}_info.json does not record a loess: the {LOESS_FILES['z']}*.asc grids beside "
                     f"{args.stem}_dem.asc are an earlier build's. Rebuild with dem_from_contours.py ... --loess")
        rec, lines = evaluate(args.stem, args.survey)
        rec["loess"] = built_with_loess(args.stem)
        for ln in lines:
            print(ln)
        if args.json:
            Path(args.json).write_text(json.dumps(rec, indent=1) + "\n")
        return 0
    if args.what == "response":
        for L in sorted({float(x) for x in args.scales.split(",")}):
            p90, p50, side = response_table(L, args.order)
            print(f"half-width {L:g} m ({'quadratic' if args.order == 2 else 'plane'}): 90% kept of a {p90:.1f} m "
                  f"wavelength, 50% of a {p50:.1f} m one; under {0.8 * L:.1f} m at most {100 * side:.0f}%")
        return 0
    ap.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())

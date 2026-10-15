#!/usr/bin/env python3
"""
Waterline Track Direction
-------------------------
A waterline follows the shore. Where a detected line runs ACROSS the shore
on the ground, it is not on the water's edge: on the station's own photos
(23 Oct 2024, 18-23 Jan 2025, 3-9 Mar 2025, 26 Sep - 2 Oct 2026) c2's line
tracer, where the edge it follows fades out, carries on as a straight run
along one pixel row out over the water, which on the ground is a track
running tens of metres seaward from the beach (13 of 46 c2 lines on
20-26 Oct 2024, 15 of 108 lines on 26 Sep - 2 Oct 2026). The consistency
filter rarely sees those runs: they lie in the columns where no other line
reaches. They enter the DEM as elevations far out on the water.

THE RULE (cross_shore_drops): per camera, the shore axis is the principal
axis of all its points on the ground. Along each line (its points in
pixel-column order) the chord from each point to the point 10 m further
along the track is measured; where that chord is at least half as long as
the 10 m walked (the track runs straight, not back and forth) and points
more than 45 deg away from the shore axis, the points it spans run across
the shore. The line's BODY is its longest stretch with none of those
points; the cross-shore points are dropped, and beyond the body, from the
first of them to the line's end (the run out over the water and whatever
the tracer drew after it). A beach cusp bends a line by far less than
45 deg over 10 m; on the photos above no c1 line and no line of a c2
stretch that follows the beach came near it.

Used by waterline_consistency.py (when its input carries easting_utm19 /
northing_utm19: the survey-date products) and dem_from_contours.py (always:
the live cron's DEMs, whose consistency filter runs before the lines are on
the ground).
"""

import numpy as np

TRACK_WINDOW_M = 10.0          # length of track over which the direction is read
TRACK_MAX_ANGLE_DEG = 45.0     # more than this from the shore axis = across the shore
TRACK_MIN_STRAIGHT = 0.5       # chord / track length: below it the track folds, no direction
SHORE_AXIS_MIN_RATIO = 2.0     # alongshore spread / cross-shore spread needed to trust the axis


def shore_axis(E, N):
    """Unit vector of the principal axis of the points (E, N), and the ratio of the
    spread along it to the spread across it; (None, ratio) when the points do not
    lie along a clear axis."""
    ok = np.isfinite(E) & np.isfinite(N)
    if ok.sum() < 10:
        return None, 0.0
    e, n = E[ok], N[ok]
    if e.size > 200000:                       # the axis needs no more
        step = e.size // 200000 + 1
        e, n = e[::step], n[::step]
    c = np.vstack([e - e.mean(), n - n.mean()])
    w, v = np.linalg.eigh(np.cov(c))
    ratio = float(np.sqrt(max(w[1], 0.0) / max(w[0], 1e-12)))
    if ratio < SHORE_AXIS_MIN_RATIO:
        return None, ratio
    return v[:, 1], ratio


def line_cross_shore(E, N, axis, window=TRACK_WINDOW_M, max_angle=TRACK_MAX_ANGLE_DEG):
    """For one line's points in pixel-column order: which points to drop (bool),
    and the length of track (m) that ran across the shore."""
    n = E.size
    drop = np.zeros(n, bool)
    ok = np.isfinite(E) & np.isfinite(N)
    idx = np.flatnonzero(ok)
    if idx.size < 3:
        return drop, 0.0
    P = np.column_stack([E[idx], N[idx]])
    step = np.hypot(*np.diff(P, axis=0).T)
    s = np.r_[0.0, np.cumsum(step)]
    j = np.searchsorted(s, s + window)          # first point at least `window` further along
    i = np.flatnonzero(j < idx.size)
    across = np.zeros(idx.size, bool)
    if i.size:
        dv = P[j[i]] - P[i]
        L = np.hypot(dv[:, 0], dv[:, 1])
        walked = s[j[i]] - s[i]
        with np.errstate(invalid="ignore", divide="ignore"):
            cosang = np.abs(dv @ axis) / L
            hit = (L >= TRACK_MIN_STRAIGHT * walked) & (cosang < np.cos(np.radians(max_angle)))
        if hit.any():
            mark = np.zeros(idx.size + 1, np.int64)
            np.add.at(mark, i[hit], 1)
            np.add.at(mark, j[i[hit]] + 1, -1)
            across = np.cumsum(mark)[:-1] > 0
    if not across.any():
        return drop, 0.0
    across_m = float(step[across[1:] & across[:-1]].sum())
    # the body: the longest stretch of track with no cross-shore point
    best, best_len, k = None, -1.0, 0
    while k < idx.size:
        if across[k]:
            k += 1
            continue
        m = k
        while m + 1 < idx.size and not across[m + 1]:
            m += 1
        if s[m] - s[k] > best_len:
            best, best_len = (k, m), s[m] - s[k]
        k = m + 1
    out = across.copy()
    if best is None:
        out[:] = True                             # the whole line runs across the shore
    else:
        a, b = best
        before = np.flatnonzero(across[:a])
        after = np.flatnonzero(across[b + 1:])
        if before.size:
            out[:before[-1] + 1] = True           # from the line's start to the last run before the body
        if after.size:
            out[b + 1 + after[0]:] = True         # from the first run after the body to the line's end
    drop[idx[out]] = True
    return drop, across_m


def cross_shore_drops(frame, col, E, N, axis=None):
    """Points to drop for every line of ONE camera: frame (line index per point),
    col (pixel column, for the order along the line), E, N (ground).
    Returns (drop mask, {line: (points dropped, metres across)}, axis, ratio);
    axis None -> the points lie along no clear shore axis, nothing judged."""
    ratio = None
    if axis is None:
        axis, ratio = shore_axis(E, N)
    drop = np.zeros(frame.size, bool)
    per_line = {}
    if axis is None:
        return drop, per_line, None, ratio
    order = np.lexsort((col, frame))
    f_sorted = frame[order]
    cuts = np.flatnonzero(np.diff(f_sorted)) + 1
    for a, b in zip(np.r_[0, cuts], np.r_[cuts, order.size]):
        sel = order[a:b]
        d, across_m = line_cross_shore(E[sel], N[sel], axis)
        if d.any():
            drop[sel[d]] = True
            per_line[int(f_sorted[a])] = (int(d.sum()), across_m)
    return drop, per_line, axis, ratio


def axis_bearing(axis):
    """Bearing (deg east of north, 0-180) of a shore axis."""
    return float(np.degrees(np.arctan2(axis[0], axis[1])) % 180.0)

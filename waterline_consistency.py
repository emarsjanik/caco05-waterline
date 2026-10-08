#!/usr/bin/env python3
"""
Waterline Consistency Filter
-------------------------------
Removes detected waterlines -- or the parts of them -- that are OUT OF
ORDER with the other waterlines the same camera saw around the same
time, before they reach the elevation maps, georectify.py and the DEM.
Reads contour_points_timex.csv, writes a filtered copy with the same
columns (contour_points_timex_qc.csv), a report of everything dropped,
and optionally a diagnostic picture of each rejection.

WHY: on 29 Sep - 5 Oct 2026 the right third of C1's 7-day map showed
high-tide lines (+1 m, red) lying SEAWARD of low-tide lines (-0.5 m,
blue). No beach does that: those were detection errors (the detector
following a feature on the water instead of the water's edge). The
DEM then saw those lines as +1 m sand where the beach is at -0.5 m,
the cells' spread went over 0.5 m and were blanked -- the southern
~80 m of the DEM, which only C1 sees, came out patchy. C2 had a
-0.75 m line out on the water at the top right of its frame. None of
the existing filters catches this: they judge each line on its own
(signal, straightness, cliffs) or a whole camera-day (tide direction),
and a wrong line can be perfectly smooth and confident.

THE PHYSICAL RULE: at a given moment the beach has one profile, and
the waterline at a higher water level lies further LANDWARD. In one
image column, the detected row is therefore a MONOTONE function of
water elevation. Which way it runs (row rising or falling with the
tide) depends on how that column cuts across the beach, so it is
found from the data, per column. A line that breaks the order by a
lot cannot be right, whatever the image looked like.

HOW, per camera:
  1. Column bins. Each line is reduced to one row per 32-px column bin
     (a least-squares line through its points in the bin, read at the
     bin centre), so lines that cover different parts of a bin are
     compared at the same column. A bin the line does not span is not
     judged for that line.
  2. Window. Each day's lines are judged against all lines from the
     day itself and the 3 days either side, WEIGHTED by closeness in
     time, w = exp(-|dt| / 1 day): the same day weighs ~1, the next
     day 0.37, three days away 0.05. The window supplies water levels a
     single day misses; the weighting makes the lines nearest in time
     the reference, so after real beach change the new beach is
     compared with itself rather than with the old one.
  3. Fit. Row as a monotone function of elevation, by weighted
     pool-adjacent-violators with the weighted MEDIAN as each block's
     value (isotonic L1 regression): a few wrong lines among many
     cannot drag the fit the way a mean would. Then the lines that
     fail (step 5) are set aside, the curve is fitted again without
     them, and every line is judged against that second fit. Without
     the refit, wrong lines at the very top or bottom of the range --
     where a block has few members to outvote them -- pulled the end of
     the curve towards themselves and made the right lines next to
     them look wrong (synthetic test: a correct +1.61 m line on the
     day after an erosion, among wrong +1.64/+1.67 m lines).
  4. Invert. The fitted curve is inverted to the elevation it gives
     each line's row, by linear interpolation between block centres.
     The RESIDUAL is that line's water elevation minus this estimate,
     in METRES -- the unit the DEM and its 0.5 m spread limit are in,
     the same anywhere in the frame. (A pixel threshold would not be:
     10 px is ~0.2 m in C1's near field and several metres in C2's far
     field.) Residual > 0: the line sits where the other lines put
     LOWER water, i.e. seaward of where its water level belongs;
     < 0: landward.
     Beyond the lines at either end a row only BOUNDS the elevation, so
     the residual there is one-sided:
       * landward of every other line: the line is at least as high as
         the highest, and no more can be said -- above the highest
         waterline lies the upper beach or berm, which can be nearly
         flat, so a little more water moves the line a long way (after
         the upper beach erodes, the next high tides land exactly
         there). Only a line whose water is LOWER than the top counts.
       * seaward of every other line: the elevation keeps falling
         seaward, at no less than half the curve's average rate per row
         (--extrapolate-flatten 2). A low-tide terrace is ~0.4x as
         steep as the beach face, but further out each row covers more
         ground, which roughly makes up for it; half leaves margin. So
         a line out on the water, whose level is too HIGH for how far
         out it lies, counts; a lowest-tide line slightly beyond the
         others does not.
  5. Decide. A bin fails when |residual| exceeds the threshold below.
     In a failing bin and the bins either side, each point is then
     judged on its own row against that bin's fit: a wrong segment
     starts part way through a bin, so this drops it up to where it
     leaves the right line and keeps the rest. The whole line is
     dropped only when more than half of its judged bins fail.

THRESHOLDS, and why:
  --max-residual 0.75 m. The errors seen were ~1.5 m out of order (a
      +1 m line among -0.5 m lines), twice this. It has to sit above
      what the beach and the water level do within the ~1-day
      effective window, and 0.5 m (DEM_MAX_SPREAD, tried first) did
      not: in the synthetic test a 0.4 m storm erosion of the upper
      beach, plus the storm's wave setup (which the water level does
      not include unless --setup-coef is used; Sep 23 2026 storm
      frames plotted ~0.5 m high), put correct lines from the day
      after up to 0.67 m out of order with the day before. At 0.75 m
      none of them was touched, nor with a 1.0 m erosion.
  --noise-k 4. Where the line barely moves with the tide (C2's far
      field: ~3 px per metre) a 2 px detection scatter is already
      ~0.7 m of elevation, so a fixed threshold would reject honest
      lines there. The threshold is the larger of --max-residual and
      4 x the robust scatter (1.4826 x MAD) of the residuals of the
      lines kept in that bin and window: honest Gaussian scatter
      exceeds 4 sigma once in ~16,000, while a line on the wrong
      feature is off by metres.
  --window-days 3, --time-scale-days 1. A storm reshapes the upper
      beach within a tide or two. With a 1-day scale the judged day
      carries ~45% of the weight and the days either side ~37%, so
      after a change the lines of the new beach outweigh the old ones
      from the next day on; the days 2-3 away (0.14, 0.05) only fill in
      water levels the nearer days lack (fog, night-time high tides).
  --column-bin 32 px. A bin averages ~32 per-column detections, so the
      bin's row is steadier than any one column, and still a wrong
      segment of ~150 px spans several bins (the C1 errors spanned
      ~500-900 px, 15-28 bins).
  --frame-fraction 0.5. "Most of it": a line that is wrong over more
      than half of what could be judged is not trusted for the rest,
      which mostly could not be judged.
  --min-frames 8 (effective). With fewer lines a block median rests on
      1-3 lines and one wrong line can be the median, so the column is
      not judged and everything in it is KEPT. Counted as the effective
      number (sum w)^2 / sum w^2, so a day with one line among a few
      distant ones (which would otherwise carry most of the weight and
      confirm itself) is not judged either.
  --min-z-range 0.5 m and --min-rank-corr 0.5. The relation must be
      clear in the data over a range of elevations: rank correlation of
      row and elevation at least 0.5 among the lines the fit keeps (the
      wrong lines would otherwise hide the relation in exactly the
      columns they spoil); otherwise the column is not judged.
  Too little data, at any of these steps, means KEEP. The filter only
  removes what the other lines clearly contradict.

  Real change is protected three ways: the time weighting (above), a
  threshold above a storm-sized change, and re-judging: the cron
  rebuilds the contour file from the whole archive every run, so a
  line wrongly dropped on the day of a change is judged again the
  next day, when there are lines from after the change to compare it
  with. Nothing is deleted from the archive.

  Safety: if a camera would lose more than --max-drop-fraction (50%)
  of its points, the filter is NOT applied to that camera (WARNING,
  report rows marked NOT APPLIED): at that point the reference itself
  is suspect -- e.g. a camera knocked out of aim -- and that is for a
  person to look at.

ELEVATION: beach_elevation_navd88 when the contour file has it (water
level + wave setup, --setup-coef), else tide_elevation_navd88 -- the
same choice, row by row, as georectify.py and dem_from_contours.py, so
the filter judges the elevations the DEM will use.

OUTPUTS:
  --output   filtered copy, same columns and row order, rows removed.
  --report   one row per line affected: frame, camera, capture time,
             elevation, points in the line and dropped, bins judged and
             failed, columns dropped, median and largest residual (m),
             threshold (m), action (frame / segments / NOT APPLIED),
             reason.
  --plot     <stem>_<camera>.png: the last --plot-days of lines on a
             photo, kept lines thin, dropped parts thick in their
             elevation colour (dashed = whole line dropped), numbered
             to a table; plus row-vs-elevation for the column bin with
             most rejections, showing the monotone fit they broke.
  The last stdout line starts "CONSISTENCY" and gives the counts per
  camera, for the cron log.

Usage:
    python3 waterline_consistency.py contour_points_timex.csv \\
        --output contour_points_timex_qc.csv \\
        --report waterline_consistency_report.csv \\
        [--plot waterline_consistency --image-dir archive/images_timex --plot-days 7]
"""

import os
import sys
import csv
import argparse
from array import array
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

csv.field_size_limit(10 ** 7)

DAY = 86400.0


# ---------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------

def read_points(path):
    """
    First pass over the contour file: per-row frame index, pixel column
    and row (float32), and per-frame name, camera, capture time, epoch
    and elevation. The file can hold millions of rows, so nothing else
    is kept; the filtered copy is written by a second pass.

    Elevation per row follows dem_from_contours.py exactly:
    beach_elevation_navd88 if present and non-empty, else
    tide_elevation_navd88. A row whose numbers do not parse is kept
    unjudged (NaN column), never dropped.
    """
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        col = {name: i for i, name in enumerate(header)}
        need = ["source_file", "camera", "pixel_column", "pixel_row", "tide_elevation_navd88"]
        missing = [n for n in need if n not in col]
        if missing:
            raise KeyError(f"{path}: missing column(s) {missing}")
        i_src, i_cam = col["source_file"], col["camera"]
        i_u, i_v = col["pixel_column"], col["pixel_row"]
        i_tide = col["tide_elevation_navd88"]
        i_beach = col.get("beach_elevation_navd88")
        i_ep, i_cap = col.get("capture_epoch"), col.get("capture_time_utc")

        frame_of = {}
        names, cams, captures, epochs, elevs, beach_used = [], [], [], [], [], []
        fid, us, vs = array("i"), array("f"), array("f")
        for r in reader:
            key = r[i_src]
            k = frame_of.get(key)
            if k is None:
                k = len(names)
                frame_of[key] = k
                names.append(key)
                cams.append(r[i_cam])
                cap = r[i_cap] if i_cap is not None else ""
                captures.append(cap)
                ep = np.nan
                try:
                    ep = float(r[i_ep]) if i_ep is not None else np.nan
                except ValueError:
                    pass
                if not np.isfinite(ep) and cap:
                    try:
                        ep = datetime.fromisoformat(cap.replace("Z", "+00:00")).timestamp()
                    except ValueError:
                        pass
                epochs.append(ep)
                z = np.nan
                b = r[i_beach] if i_beach is not None else ""
                try:
                    z = float(b or r[i_tide])
                except ValueError:
                    pass
                elevs.append(z)
                beach_used.append(bool(b))
            fid.append(k)
            try:
                us.append(float(r[i_u])); vs.append(float(r[i_v]))
            except ValueError:
                us.append(np.nan); vs.append(np.nan)
    elev_col = "beach_elevation_navd88" if i_beach is not None else "tide_elevation_navd88"
    return dict(header=header, frame=np.frombuffer(fid, dtype=np.int32),
                col=np.frombuffer(us, dtype=np.float32), row=np.frombuffer(vs, dtype=np.float32),
                names=names, cams=np.array(cams), captures=captures,
                epochs=np.array(epochs, float), elevs=np.array(elevs, float),
                beach_used=np.array(beach_used, bool), elev_col=elev_col)


# ---------------------------------------------------------------------
# Monotone fit
# ---------------------------------------------------------------------

def weighted_median(values, weights):
    if len(values) == 1:
        return float(values[0])
    order = np.argsort(values, kind="stable")
    v, w = values[order], weights[order]
    c = np.cumsum(w)
    half = 0.5 * float(c[-1])
    i = int(np.searchsorted(c, half))
    if i < len(v) - 1 and abs(float(c[i]) - half) <= 1e-12 * half:
        return 0.5 * float(v[i] + v[i + 1])         # weight splits exactly between two
    return float(v[min(i, len(v) - 1)])


def isotonic_l1(y, w):
    """
    Weighted L1 isotonic regression (non-decreasing) by
    pool-adjacent-violators, the block value being the weighted median
    of its members (Robertson, Wright & Dykstra 1988: PAV with the
    median solves the L1 problem). y and w are in order of the
    independent variable. Returns block (start, end) index pairs and
    values.
    """
    starts, ends, vals = [], [], []
    for i in range(len(y)):
        starts.append(i); ends.append(i + 1); vals.append(y[i])
        while len(vals) > 1 and vals[-2] > vals[-1]:
            e = ends.pop(); vals.pop(); starts.pop()
            ends[-1] = e
            s = starts[-1]
            vals[-1] = weighted_median(y[s:e], w[s:e])
    return starts, ends, vals


def rank_corr(a, b):
    """Spearman rank correlation (no tie averaging; the inputs are continuous)."""
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    den = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / den) if den > 0 else 0.0


def fit_column(z, y, w):
    """
    Monotone fit of row (y) against elevation (z), weights w.
    Returns (sign, knot_rows, knot_elevations) with knot rows strictly
    increasing after multiplying by sign, or None if there is only one
    block (no relation to invert).
    """
    sgn = 1.0 if rank_corr(z, y) >= 0 else -1.0
    order = np.lexsort((sgn * y, z))       # ties in z: no false violation
    zs, ys, ws = z[order], sgn * y[order], w[order]
    starts, ends, vals = isotonic_l1(ys, ws)
    # Merge equal-valued neighbours so the curve can be inverted.
    groups = []
    for s, e, v in zip(starts, ends, vals):
        if groups and abs(groups[-1][2] - v) <= 1e-9:
            groups[-1][1] = e
        else:
            groups.append([s, e, v])
    if len(groups) < 2:
        return None
    knot_y = np.array([g[2] for g in groups])
    knot_z = np.array([weighted_median(zs[s:e], ws[s:e]) for s, e, _ in groups])
    return sgn, knot_y, knot_z


def residuals(fit, z, y, flatten):
    """
    Water elevation minus the elevation the fitted curve gives each row (m).

    Between the knots: linear interpolation, two-sided.

    Beyond the ends a row only BOUNDS the elevation, so the residual is
    one-sided there:
      * landward of every level (beyond the highest line): the line is at
        least as high as the top knot, and nothing more can be said -- above
        the highest waterline lies the upper beach or berm, which can be
        nearly flat, so a little more water moves the line a long way.
        After the upper beach erodes, the next high tides land exactly
        there. Only a line BELOW the top level is out of order.
      * seaward of every level (beyond the lowest line): the beach goes on
        down seaward, but may flatten (a low-tide terrace is ~0.4x the
        slope of the beach face). The elevation is taken to fall at no
        less than 1/`flatten` of the fit's average rate, so only a line
        whose water level is too HIGH for how far out it lies counts -- a
        line far out on the water, which no flattening explains.
    """
    sgn, ky, kz = fit
    yy = sgn * np.asarray(y, float)
    r = z - np.interp(yy, ky, kz)
    span = ky[-1] - ky[0]
    slope = (kz[-1] - kz[0]) / span if span > 0 else 0.0
    seaward = yy < ky[0]
    z_max = kz[0] - (ky[0] - yy) * slope / flatten
    r = np.where(seaward, np.maximum(0.0, z - z_max), r)
    landward = yy > ky[-1]
    r = np.where(landward, np.minimum(0.0, z - kz[-1]), r)
    return r


def judge_bin(z, y, w, p):
    """
    Fit, set aside the lines that fail, refit without them, and judge every
    line against the refit (one robust iteration). Without the refit, wrong
    lines at the very top or bottom of the elevation range -- where a block
    has few members -- can carry the end of the fit with them and make the
    right lines next to them look wrong.

    The column is judged only if row clearly follows elevation: rank
    correlation at least --min-rank-corr among the lines the refit keeps
    (and at least 0.2 among all, for the direction). Measured on all
    lines, the wrong lines themselves pull it down: on the synthetic test
    the columns where 30% of high-tide lines were wrong fell to 0.47 and
    went unjudged -- exactly the columns that needed it.

    Returns (fit, residuals, threshold) or None when the column cannot be
    judged.
    """
    if abs(rank_corr(z, y)) < 0.2:
        return None
    fit = fit_column(z, y, w)
    if fit is None:
        return None

    def threshold(r, ww):
        med = weighted_median(r, ww)
        sigma = 1.4826 * weighted_median(np.abs(r - med), ww)
        return max(p.max_residual, p.noise_k * sigma)

    r = residuals(fit, z, y, p.extrapolate_flatten)
    thr = threshold(r, w)
    keep = np.abs(r) <= thr
    if not keep.all() and keep.sum() >= p.min_frames \
            and w[keep].sum() ** 2 / (w[keep] ** 2).sum() >= p.min_frames:
        fit2 = fit_column(z[keep], y[keep], w[keep])
        if fit2 is not None:
            fit = fit2
            r = residuals(fit, z, y, p.extrapolate_flatten)
            thr = threshold(r[keep], w[keep])
    if abs(rank_corr(z[keep], y[keep])) < p.min_rank_corr:
        return None
    return fit, r, thr


# ---------------------------------------------------------------------
# Judging one camera
# ---------------------------------------------------------------------

def bin_rows(frame_local, col, row, n_frames, bin_px):
    """
    Each line's row at each column-bin centre: least-squares line through
    the line's points in the bin, read at the centre. Only bins the line
    actually spans (points on both sides of the centre) get a value.
    """
    ok = np.isfinite(col) & np.isfinite(row)
    j = np.floor(col[ok] / bin_px).astype(np.int64)
    n_bins = int(j.max()) + 1 if j.size else 1
    x = col[ok].astype(float) - (j + 0.5) * bin_px
    y = row[ok].astype(float)
    key = frame_local[ok].astype(np.int64) * n_bins + j
    size = n_frames * n_bins
    n = np.bincount(key, minlength=size).astype(float)
    sx = np.bincount(key, x, size); sy = np.bincount(key, y, size)
    sxx = np.bincount(key, x * x, size); sxy = np.bincount(key, x * y, size)
    xmin = np.full(size, np.inf); xmax = np.full(size, -np.inf)
    np.minimum.at(xmin, key, x); np.maximum.at(xmax, key, x)
    with np.errstate(invalid="ignore", divide="ignore"):
        det = n * sxx - sx * sx
        slope = np.where(det > 1e-9, (n * sxy - sx * sy) / det, 0.0)
        centre = (sy - slope * sx) / n
    spans = (n >= 2) & (xmin <= -bin_px / 8.0) & (xmax >= bin_px / 8.0)
    R = np.where(spans, centre, np.nan).reshape(n_frames, n_bins)
    return R, n_bins


def judge_camera(epochs, elevs, frame_local, col, row, n_frames, p):
    """
    Returns per frame x bin: residual (m), threshold (m), the bin rows, and
    the window fits {(day, bin): (fit, frames, weights, threshold)}.
    """
    R, n_bins = bin_rows(frame_local, col, row, n_frames, p.column_bin)
    resid = np.full((n_frames, n_bins), np.nan)
    thresh = np.full((n_frames, n_bins), np.nan)
    days = np.floor(epochs / DAY)
    good_frame = np.isfinite(epochs) & np.isfinite(elevs)
    fits = {}
    for d in np.unique(days[good_frame]):
        judged = np.where(good_frame & (days == d))[0]
        window = np.where(good_frame & (np.abs(days - d) <= p.window_days))[0]
        tc = np.median(epochs[judged])
        w_all = np.exp(-np.abs(epochs[window] - tc) / (p.time_scale_days * DAY))
        bins = np.where(np.isfinite(R[judged]).any(axis=0))[0]
        for b in bins:
            have = np.isfinite(R[window, b])
            idx = window[have]
            if idx.size < p.min_frames:
                continue
            w = w_all[have]
            if w.sum() ** 2 / (w * w).sum() < p.min_frames:
                continue
            z, y = elevs[idx], R[idx, b]
            if z.max() - z.min() < p.min_z_range:
                continue
            judged_bin = judge_bin(z, y, w, p)
            if judged_bin is None:
                continue
            fit, r_all, thr = judged_bin
            mine = np.isin(idx, judged)
            resid[idx[mine], b] = r_all[mine]
            thresh[idx[mine], b] = thr
            fits[(d, b)] = (fit, idx, w, thr)
    return resid, thresh, R, fits


def decide(resid, thresh, p):
    """Failing bins per frame, and the frames to drop whole."""
    judged = np.isfinite(resid)
    fail = judged & (np.abs(resid) > thresh)
    n_judged = judged.sum(axis=1)
    n_fail = fail.sum(axis=1)
    whole = (n_judged >= p.min_judged_bins) & (n_fail > p.frame_fraction * n_judged)
    return fail, whole, n_judged, n_fail


def point_drops(fl, col, row, fail, whole, days, elevs, fits, p):
    """
    Which points to drop. Whole frames: all of them. Otherwise each point
    in a failing bin, or in a bin next to one, is judged on its OWN row
    against that bin's fit and threshold. A wrong segment starts and ends
    part way through a bin, where the bin's single row is a blend of the
    right and the wrong part; judging the points there one by one drops the
    wrong part up to where it leaves the right line, and keeps the rest.
    Points away from any failing bin are never judged one by one, so pixel
    noise on a good line cannot speckle it with drops.
    """
    drop = whole[fl].copy()
    ok = np.isfinite(col)
    n_bins = fail.shape[1]
    j = np.full(col.shape, -1, np.int64)
    j[ok] = np.floor(col[ok] / p.column_bin).astype(np.int64)
    jc = np.clip(j, 0, n_bins - 1)
    near = ok & (j >= 0) & (j < n_bins) & fail[fl, jc]
    near |= ok & (j - 1 >= 0) & (j - 1 < n_bins) & fail[fl, np.clip(j - 1, 0, n_bins - 1)]
    near |= ok & (j + 1 < n_bins) & fail[fl, np.clip(j + 1, 0, n_bins - 1)]
    near &= ~drop
    cand = np.where(near)[0]
    if cand.size == 0:
        return drop
    key_day = days[fl[cand]]
    keys = np.stack([key_day, jc[cand]], axis=1)
    uniq, inv = np.unique(keys, axis=0, return_inverse=True)
    inv = inv.ravel()
    for k, (d, b) in enumerate(uniq):
        hit = fits.get((d, int(b)))
        if hit is None:
            continue                         # bin not judged that day: keep
        fit, _, _, thr = hit
        sel = cand[inv == k]
        r = residuals(fit, elevs[fl[sel]], row[sel].astype(float), p.extrapolate_flatten)
        drop[sel[np.abs(r) > thr]] = True
    return drop


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def run(args):
    data = read_points(args.contour_csv)
    n_rows = data["frame"].size
    keep = np.ones(n_rows, bool)
    cams = data["cams"]
    report = []
    summary = []
    plot_info = {}
    print(f"Contour points     : {n_rows:,} rows, {len(data['names']):,} lines  ({args.contour_csv})")
    print(f"Elevation judged   : {data['elev_col']}"
          + (" (falls back to tide_elevation_navd88 where blank)"
             if data["elev_col"] == "beach_elevation_navd88" else ""))
    print(f"Window             : +/-{args.window_days:g} d, weight exp(-|dt|/{args.time_scale_days:g} d); "
          f"column bins {args.column_bin} px")
    print(f"Thresholds         : |residual| > max({args.max_residual:g} m, {args.noise_k:g} x robust "
          f"scatter); whole line when > {args.frame_fraction:.0%} of judged bins fail")
    print()
    for cam in sorted(set(cams.tolist())):
        if args.camera and cam != args.camera:
            continue
        frames = np.where(cams == cam)[0]
        local = np.full(len(data["names"]), -1, np.int64)
        local[frames] = np.arange(frames.size)
        pts = np.where(local[data["frame"]] >= 0)[0]
        fl = local[data["frame"][pts]]
        col, row = data["col"][pts], data["row"][pts]
        epochs, elevs = data["epochs"][frames], data["elevs"][frames]
        resid, thresh, R, fits = judge_camera(epochs, elevs, fl, col, row, frames.size, args)
        fail, whole, n_judged, n_fail = decide(resid, thresh, args)

        drop = point_drops(fl, col, row, fail, whole, np.floor(epochs / DAY), elevs, fits, args)
        n_cam_pts = pts.size
        n_drop = int(drop.sum())
        applied = n_drop <= args.max_drop_fraction * n_cam_pts
        if applied:
            keep[pts[drop]] = False

        pts_per_frame = np.bincount(fl, minlength=frames.size)
        drop_per_frame = np.bincount(fl[drop], minlength=frames.size)
        affected = np.where(drop_per_frame > 0)[0]
        n_whole = int(whole[affected].sum()) if affected.size else 0
        cmin = np.full(frames.size, np.inf); cmax = np.full(frames.size, -np.inf)
        dc = drop & np.isfinite(col)
        np.minimum.at(cmin, fl[dc], col[dc]); np.maximum.at(cmax, fl[dc], col[dc])
        for f in affected:
            rr = resid[f][fail[f]] if fail[f].any() else resid[f][np.isfinite(resid[f])]
            med = float(np.median(rr)) if rr.size else float("nan")
            worst = float(rr[np.argmax(np.abs(rr))]) if rr.size else float("nan")
            thr = float(np.nanmedian(thresh[f][fail[f]])) if fail[f].any() else float("nan")
            side = "lower" if med > 0 else "higher"
            where = "seaward" if med > 0 else "landward"
            if whole[f]:
                action = "frame"
                reason = (f"{n_fail[f]} of {n_judged[f]} judged column bins out of order "
                          f"(> {args.frame_fraction:.0%}): whole line dropped; median residual "
                          f"{med:+.2f} m, i.e. it sits {where} of where its water level belongs")
            else:
                action = "segments"
                reason = (f"out of order in {n_fail[f]} of {n_judged[f]} judged column bins: the line "
                          f"sits where the other lines put water {abs(med):.2f} m {side} "
                          f"({where} of its level; threshold {thr:.2f} m)")
            if not applied:
                action = "NOT APPLIED"
            report.append(dict(
                source_file=data["names"][frames[f]], camera=cam,
                capture_time_utc=data["captures"][frames[f]],
                elevation_navd88=round(float(elevs[f]), 4),
                elevation_column=("beach_elevation_navd88" if data["beach_used"][frames[f]]
                                  else "tide_elevation_navd88"),
                action=action, points_in_frame=int(pts_per_frame[f]),
                points_dropped=int(drop_per_frame[f]),
                bins_judged=int(n_judged[f]), bins_failed=int(n_fail[f]),
                pixel_column_from=round(float(cmin[f]), 1) if np.isfinite(cmin[f]) else "",
                pixel_column_to=round(float(cmax[f]), 1) if np.isfinite(cmax[f]) else "",
                median_residual_m=round(med, 3), max_residual_m=round(worst, 3),
                threshold_m=round(thr, 3) if np.isfinite(thr) else "",
                reason=reason))
        judged_lines = int((n_judged > 0).sum())
        msg = (f"{cam}: {n_whole} line(s) dropped, {affected.size - n_whole} trimmed, "
               f"{n_drop:,} of {n_cam_pts:,} points ({100.0 * n_drop / max(n_cam_pts, 1):.1f}%); "
               f"{judged_lines} of {frames.size} lines judged")
        if not applied:
            msg = (f"{cam}: NOT APPLIED -- would drop {n_drop:,} of {n_cam_pts:,} points "
                   f"(> {args.max_drop_fraction:.0%}), passed through unfiltered")
            print(f"WARNING: camera {cam}: the filter would drop {100.0 * n_drop / n_cam_pts:.0f}% "
                  f"of its points, more than --max-drop-fraction {args.max_drop_fraction:g}. That")
            print("         is not a few bad detections; the reference itself is suspect (camera")
            print("         moved? detector failing?). Passed through UNFILTERED -- look at the report.")
        print(f"  {msg}")
        summary.append(msg)
        plot_info[cam] = dict(frames=frames, fl=fl, col=col, row=row, drop=drop, whole=whole,
                              resid=resid, thresh=thresh, fail=fail, R=R, fits=fits,
                              epochs=epochs, elevs=elevs, applied=applied)
    return data, keep, report, summary, plot_info


def write_filtered(src, dst, keep):
    """Second pass: copy the rows kept, unchanged, through a temporary file
    that replaces the output only when complete (and is removed if not)."""
    tmp = Path(str(dst) + ".tmp")
    try:
        with open(src, newline="") as fi, open(tmp, "w", newline="") as fo:
            reader = csv.reader(fi)
            writer = csv.writer(fo)
            writer.writerow(next(reader))
            for i, r in enumerate(reader):
                if keep[i]:
                    writer.writerow(r)
        os.replace(tmp, dst)
    finally:
        if tmp.exists():
            tmp.unlink()


REPORT_FIELDS = ["source_file", "camera", "capture_time_utc", "elevation_navd88", "elevation_column",
                 "action", "points_in_frame", "points_dropped", "bins_judged", "bins_failed",
                 "pixel_column_from", "pixel_column_to", "median_residual_m", "max_residual_m",
                 "threshold_m", "reason"]


def write_report(path, report):
    tmp = Path(str(path) + ".tmp")
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=REPORT_FIELDS)
        w.writeheader()
        for r in sorted(report, key=lambda r: (r["camera"], r["capture_time_utc"])):
            w.writerow(r)
    os.replace(tmp, path)


# ---------------------------------------------------------------------
# Diagnostic plot
# ---------------------------------------------------------------------

def _runs(cols, rows, max_gap=40):
    """Split a line into runs at column gaps, so a gap stays a gap."""
    if cols.size == 0:
        return []
    order = np.argsort(cols)
    c, r = cols[order], rows[order]
    cut = np.flatnonzero(np.diff(c) > max_gap) + 1
    return [(c[a:b], r[a:b]) for a, b in zip(np.r_[0, cut], np.r_[cut, c.size]) if b - a >= 2]


def plot_camera(cam, info, data, report, args, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    frames, fl, col, row, drop = info["frames"], info["fl"], info["col"], info["row"], info["drop"]
    epochs, elevs = info["epochs"], info["elevs"]
    last_day = np.floor(np.nanmax(epochs) / DAY)
    in_win = np.floor(epochs / DAY) > last_day - args.plot_days
    win_frames = np.where(in_win)[0]
    first = datetime.fromtimestamp((last_day - args.plot_days + 1) * DAY, tz=timezone.utc).date()
    last = datetime.fromtimestamp(last_day * DAY, tz=timezone.utc).date()

    # Backdrop: the clearest photo of the last day, as the elevation maps choose it.
    image = None
    bg_note = "no photo found"
    if args.image_dir:
        try:
            import cv2
            from daily_elevation_map import pick_background
            cand = {data["names"][frames[f]]: {"capture": data["captures"][frames[f]]}
                    for f in win_frames if np.floor(epochs[f] / DAY) == last_day}
            path, cap, _ = pick_background(cand, args.image_dir, 12, 15, camera=cam)
            if path is not None:
                image = cv2.cvtColor(cv2.imread(str(path)), cv2.COLOR_BGR2RGB)
                bg_note = f"photo {cap[:16].replace('T', ' ')} UTC"
        except Exception as exc:          # a missing photo must not stop the plot
            bg_note = f"no photo ({exc.__class__.__name__})"
    width = image.shape[1] if image is not None else int(np.nanmax(col)) + 1
    height = image.shape[0] if image is not None else int(np.nanmax(row) * 1.1) + 1

    norm = Normalize(vmin=np.nanmin(elevs[win_frames]), vmax=np.nanmax(elevs[win_frames]))
    try:
        cmap = matplotlib.colormaps["turbo"]        # matplotlib >= 3.5
    except AttributeError:
        cmap = plt.get_cmap("turbo")                # 3.3 / 3.4 (the station's floor)

    fig = plt.figure(figsize=(17, 9.6), dpi=110)
    ax = fig.add_axes([0.01, 0.06, 0.62, 0.86])
    if image is not None:
        ax.imshow(image)
    else:
        ax.set_facecolor("0.85")
    rep = {r["source_file"]: r for r in report if r["camera"] == cam}
    listed = []
    for f in win_frames:
        sel = fl == f
        if not sel.any():
            continue
        k = sel & ~drop
        for c, r in _runs(col[k], row[k]):
            ax.plot(c, r, color="white", alpha=0.22, lw=0.5)
    for f in win_frames:
        sel = (fl == f) & drop
        if not sel.any():
            continue
        colour = cmap(norm(elevs[f]))
        whole = bool(info["whole"][f])
        for c, r in _runs(col[sel], row[sel]):
            ax.plot(c, r, color="black", lw=3.6, alpha=0.6, ls="--" if whole else "-")
            ax.plot(c, r, color=colour, lw=2.2, ls="--" if whole else "-")
        listed.append(f)
    # Numbers in two staggered rows above the lines, spread sideways with a
    # leader to each dropped part, so rejections in one place stay readable.
    if listed:
        top = float(np.nanmin(row[drop & np.isin(fl, listed)]))
        anchors = []
        for f in listed:
            sel = (fl == f) & drop
            c, r = col[sel], row[sel]
            i = np.argsort(c)[len(c) // 2]
            anchors.append((float(c[i]), float(r[i])))
        order = list(np.argsort([a[0] for a in anchors]))
        sep = 0.032 * width
        for tier in (0, 1):
            ids = order[tier::2]
            xs = [anchors[n][0] for n in ids]
            for k in range(1, len(xs)):                      # spread rightwards
                xs[k] = max(xs[k], xs[k - 1] + sep)
            if xs and xs[-1] > width - sep / 2:              # then pull back inside
                xs[-1] = width - sep / 2
                for k in range(len(xs) - 2, -1, -1):
                    xs[k] = min(xs[k], xs[k + 1] - sep)
            ly = max(top - (0.05 + 0.035 * tier) * height, 0.03 * height)
            for n, lx in zip(ids, xs):
                lx = max(lx, sep / 2)
                x, y = anchors[n]
                ax.plot([lx, x], [ly, y], color="white", lw=0.5, alpha=0.8)
                ax.text(lx, ly, str(n + 1), fontsize=7.5, fontweight="bold", ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="0.3", lw=0.4, alpha=0.9))
    ax.set_xlim(0, width); ax.set_ylim(height, 0); ax.axis("off")
    status = "" if info["applied"] else "   [NOT APPLIED: would drop too much, see log]"
    ax.set_title(f"{cam.upper()}  {first} to {last}: {len(listed)} line(s) with parts dropped as out of "
                 f"order (thick; dashed = whole line){status}\nthin white = kept lines, {bg_note}",
                 fontsize=10)
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cax = fig.add_axes([0.03, 0.035, 0.25, 0.015])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("water elevation of the line (m NAVD88)", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # Row vs elevation at the column bin with the most rejections in the window.
    ax2 = fig.add_axes([0.68, 0.56, 0.30, 0.36])
    fail = info["fail"][win_frames]
    if fail.any():
        b = int(np.argmax(fail.sum(axis=0)))
        day = max((d for (d, bb) in info["fits"] if bb == b and d <= last_day
                   and d > last_day - args.plot_days),
                  key=lambda d: int(info["fail"][np.floor(epochs / DAY) == d][:, b].sum()), default=None)
        R = info["R"][:, b]
        if day is not None:
            fit, idx, w, _ = info["fits"][(day, b)]
            bad = info["fail"][idx, b] | info["whole"][idx]
            ax2.scatter(elevs[idx][~bad], R[idx][~bad], s=10 + 30 * w[~bad], c="0.45",
                        label="lines in the window (size = weight)")
            ax2.scatter(elevs[idx][bad], R[idx][bad], s=40, marker="x", c="red", label="dropped")
            sgn, ky, kz = fit
            ax2.plot(kz, sgn * ky, "-", color="tab:blue", lw=1.5, label="monotone fit")
            dd = datetime.fromtimestamp(day * DAY, tz=timezone.utc).date()
            ax2.set_title(f"columns {b * args.column_bin}-{(b + 1) * args.column_bin - 1} px: the fit "
                          f"that judged {dd} (window +/-{args.window_days:g} d)", fontsize=9)
            ax2.invert_yaxis()
            ax2.legend(fontsize=7, loc="best")
    else:
        ax2.text(0.5, 0.5, "nothing dropped in this window", ha="center", va="center",
                 transform=ax2.transAxes)
    ax2.set_xlabel("water elevation (m NAVD88)", fontsize=8)
    ax2.set_ylabel("detected row (px), image orientation", fontsize=8)
    ax2.tick_params(labelsize=7)
    ax2.grid(alpha=0.3)

    ax3 = fig.add_axes([0.66, 0.02, 0.33, 0.47]); ax3.axis("off")
    lines = ["#   capture (UTC)     elev    pts   resid  action"]
    for n, f in enumerate(listed, 1):
        r = rep.get(data["names"][frames[f]])
        if r is None:
            continue
        lines.append(f"{n:<3d} {r['capture_time_utc'][5:16].replace('T', ' ')}  {r['elevation_navd88']:+5.2f}"
                     f"  {r['points_dropped']:5d}  {r['median_residual_m']:+5.2f}  {r['action']}")
        if n >= 34:
            lines.append(f"... and {len(listed) - n} more: see the report CSV")
            break
    lines.append("")
    lines.append("resid = water level minus the elevation the other lines give its row (m);")
    lines.append("> 0: line lies seaward of where its water level belongs, < 0 landward.")
    ax3.text(0, 1, "\n".join(lines), va="top", ha="left", family="monospace", fontsize=7.4)
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description="Drop waterlines (or parts) that are out of order with the other lines: "
                    "row must be monotone in water elevation per image column.")
    ap.add_argument("contour_csv")
    ap.add_argument("--output", default=None,
                    help="Filtered copy (default: <input>_qc.csv next to the input).")
    ap.add_argument("--report", default=None,
                    help="CSV of every line affected (default: waterline_consistency_report.csv "
                         "next to the output).")
    ap.add_argument("--plot", default=None,
                    help="Write <PLOT>_<camera>.png diagnostics of the last --plot-days.")
    ap.add_argument("--image-dir", default=None, help="Timex photos for the plot backdrop.")
    ap.add_argument("--plot-days", type=int, default=7, help="Days shown in the plot (default 7).")
    ap.add_argument("--camera", default=None, help="Judge only this camera (others pass through).")
    ap.add_argument("--column-bin", type=int, default=32,
                    help="Width (px) of the column bins lines are compared in (default 32).")
    ap.add_argument("--window-days", type=float, default=3.0,
                    help="Lines from this many days either side form the reference (default 3).")
    ap.add_argument("--time-scale-days", type=float, default=1.0,
                    help="Reference weight exp(-|dt|/this) (default 1 day).")
    ap.add_argument("--max-residual", type=float, default=0.75,
                    help="Elevation residual (m) beyond which a line is out of order (default 0.75; "
                         "raised where the scatter is larger, see --noise-k).")
    ap.add_argument("--extrapolate-flatten", type=float, default=2.0,
                    help="Seaward of every line, elevation falls at no less than 1/this of the "
                         "average rate (default 2).")
    ap.add_argument("--noise-k", type=float, default=4.0,
                    help="Threshold is at least this many robust sigmas of the bin's residuals "
                         "(default 4).")
    ap.add_argument("--min-frames", type=float, default=8.0,
                    help="Effective number of lines needed to judge a bin; fewer = keep (default 8).")
    ap.add_argument("--min-z-range", type=float, default=0.5,
                    help="Elevation range (m) the lines must span to judge a bin (default 0.5).")
    ap.add_argument("--min-rank-corr", type=float, default=0.5,
                    help="Rank correlation of row and elevation needed to judge a bin (default 0.5).")
    ap.add_argument("--min-judged-bins", type=int, default=4,
                    help="Judged bins a line needs before it can be dropped whole (default 4).")
    ap.add_argument("--frame-fraction", type=float, default=0.5,
                    help="Drop the whole line when more than this share of its judged bins fail "
                         "(default 0.5).")
    ap.add_argument("--max-drop-fraction", type=float, default=0.5,
                    help="If a camera would lose more than this share of its points, pass it "
                         "through unfiltered and warn (default 0.5).")
    args = ap.parse_args()

    src = Path(args.contour_csv)
    out = Path(args.output) if args.output else src.with_name(src.stem + "_qc.csv")
    rep_path = Path(args.report) if args.report else out.with_name("waterline_consistency_report.csv")
    if out.resolve() == src.resolve():
        print("ERROR: --output must differ from the input (the input is the unfiltered record).")
        sys.exit(1)

    data, keep, report, summary, plot_info = run(args)
    write_filtered(src, out, keep)
    write_report(rep_path, report)
    print()
    print(f"Filtered copy      : {out}  ({int(keep.sum()):,} of {keep.size:,} rows kept)")
    print(f"Report             : {rep_path}  ({len(report)} line(s) affected)")
    if args.plot:
        for cam, info in plot_info.items():
            png = f"{args.plot}_{cam}.png"
            try:
                plot_camera(cam, info, data, report, args, png)
                print(f"Diagnostic plot    : {png}")
            except Exception as exc:     # the filtered file is the product; a plot failure is not
                print(f"WARNING: diagnostic plot for {cam} failed: {exc!r}")
    print("CONSISTENCY " + "; ".join(summary))


if __name__ == "__main__":
    main()

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
the waterline at a higher water level lies further LANDWARD -- nearer
the cameras on the bluff, so LOWER in the photo. In one image column the
detected row therefore RISES with water elevation (ROW_SIGN; checked
with the calibration for every column of c1 and c2). A line that breaks
the order by a lot cannot be right, whatever the image looked like.

The direction used to be read from the data in each column. That fails
exactly where it matters most: where the detector is wrong for MOST of
a column's lines in the same way, the wrong lines set the direction and
become the reference. A review (Oct 2026) projected the 29 Sep 2026
RTK checkshots into c1: on its right half +1.17 m lands ~33 px below
the floor of the detector's search envelope, so high-water lines there
may have no right row to go to and land seaward, in every high-tide
frame. Synthetic test on the RTK-surveyed beach, every c1 line that
leaves the search envelope redrawn seaward: the data-derived direction
caught 41% of the wrong points (28 of 136 lines); the physical
direction with the bottom-up pass below, 98% (132). In the station's
own 29 Sep - 5 Oct contours it has not gone that far -- rows still rise
with the water in every c1 column, with a quarter to a third of the
high-tide lines seaward of the low-tide ones -- and the filter drops
about what it did before (22 lines and 18 trimmed in c1, against 21
and 17). Physics does not take a vote either way.

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
  2b. Bottom-up, where the order is not clear (rank correlation of row
     and elevation below --sweep-corr 0.9): the lines are judged in
     layers of rising water level (0.25 m), each layer only against
     the lines at LOWER water already kept, with the fit and residual
     of steps 3-4 (bottom_up()). A line can fail here only by sitting
     SEAWARD of where lower water put the waterline -- so where the
     detector got most HIGH-tide lines wrong, the low-tide lines still
     judge them, however many they are -- and only when the lower
     lines BEFORE it in time and those AFTER it both say so: a change
     in the beach (a storm cut) moves the lower lines on one side of a
     line only, so lines from before a cut are never judged by the cut
     beach alone. Where one side has too few lines, the lower lines of
     the line's own tide (+/- 12.42 h) judge it. It must also lie more
     than 4 x 2 px (--noise-k x --row-noise-px) seaward of them in the
     photo, so where the waterline barely moves with the tide (c2's far
     field, a crenulated shoreline) a few pixels are not read as
     metres. The lines this keeps seed step 3.
     A bin whose rank correlation is at or below -0.2 (--reversed-corr)
     is REVERSED: rows run AGAINST the water level, which no beach
     does. It is counted, named in a WARNING with its camera and pixel
     columns, and marked on the diagnostic plot; its lines are judged
     as above (and, if even the kept lines show no clear order, by the
     bottom-up pass alone) -- never kept silently.
  3. Fit. Row as a monotone (rising) function of elevation, by weighted
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
     The most seaward line of a bin is the end of the fit, which cannot
     contradict it: a stray out on the water at the window's LOWEST
     water level became the bottom knot itself, residual 0. So the most
     seaward lines are judged one at a time against all the others,
     from the seaward end in, until one passes (seaward_end()).
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
      feature is off by metres. The bottom-up pass's reference holds
      only lower lines, whose misfit to their own curve says nothing of
      the lines above, so there the scatter is taken out of sample
      (elevation_scatter(): lines of nearly the same row, how far apart
      in elevation) and the line must also lie 4 x --row-noise-px
      seaward in the photo.
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
      row and elevation (in the physical direction) at least 0.5 among
      the lines the fit keeps -- the wrong lines would otherwise hide
      the relation in exactly the columns they spoil; otherwise the
      column is not judged ('unclear').
  --sweep-corr 0.9, --reversed-corr 0.2, --layer-m 0.25. Honest lines
      in a column are ordered: rank correlation 0.93-0.99 in c1 on the
      synthetic RTK-surveyed beach with 1.6 px detection scatter, 0.85-
      0.9 in a few of c2's far-field bins and in the week of a storm
      cut. Below 0.9 the bottom-up pass runs first. The station's real
      lines are far less tidy -- 0.4-0.6 per column over 29 Sep - 5 Oct
      2026 -- so there it runs in ~90% of the bins, and the filter takes
      ~1.7x as long as without it (7 s against 4 s for that week's
      390,000 points). At -0.2 the order is reversed, not noisy. 0.25 m
      layers keep each line within a quarter metre of the lower lines
      that judge it, well inside the 0.75 m threshold.
  --row-noise-px 2. The detector's own row scatter between lines at the
      same water level: c2's bin rows differ by 1.8 px (robust sigma)
      between consecutive slack-water frames of 29 Sep - 5 Oct 2026. The
      bottom-up pass drops a line only when it lies more than --noise-k
      times this seaward of the lower lines in the photo (bottom_up()).
  Too little data, at any of these steps, means KEEP. The filter only
  removes what the other lines clearly contradict.

  Real change is protected four ways: the time weighting (above), a
  threshold above a storm-sized change, the bottom-up pass's two sides
  in time (a line is never dropped by lower lines from across a change
  alone), and re-judging: the cron
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
             reason, and how many of its failing bins were REVERSED.
  --plot     <stem>_<camera>.png: the last --plot-days of lines on a
             photo, kept lines thin, dropped parts thick in their
             elevation colour (dashed = whole line dropped), numbered
             to a table, and a red bar over columns that were REVERSED;
             plus row-vs-elevation for the column bin with most
             rejections, showing the monotone fit they broke.
  The last stdout line starts "CONSISTENCY" and gives the counts per
  camera, for the cron log: lines dropped and trimmed, points dropped,
  and the COLUMN BINS (32 px x day) judged, not judged ('too few lines',
  'unclear': everything in them kept) and REVERSED. Counting lines would
  hide a column nobody could judge: a line counts as judged if any of
  its ~70 bins was. A WARNING line before it names the camera and pixel
  columns of any REVERSED bins.

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
TIDE_CYCLE = 12.42 * 3600.0          # one semidiurnal (M2) tide, s

# Which way the detected row moves as the water rises, in every column of
# every camera here: DOWN the image (row increases). Higher water puts the
# waterline further landward, nearer the cameras on the bluff, and nearer
# ground is lower in a photo taken looking down at it. Checked with the
# calibration (IO 20240801, EO 20251113-CV) on the beach surveyed on 29 Sep
# 2026: row increases with water elevation in 100% of the column steps for
# c1 and c2 alike (>= 34 px/m in c1, >= 4.7 px/m in c2's far field). It was
# once taken from the data in each column, but where most lines are wrong
# the data say the opposite, and the wrong lines became the reference.
ROW_SIGN = 1.0


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
    Monotone fit of row (y) against elevation (z), weights w, in the
    PHYSICAL direction: row rises with elevation (ROW_SIGN). Returns
    (sign, knot_rows, knot_elevations) with knot rows strictly increasing
    after multiplying by sign, or None if there is only one block (no
    relation to invert).
    """
    sgn = ROW_SIGN
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


def n_eff(w):
    """Effective number of lines, (sum w)^2 / sum w^2."""
    return float(w.sum() ** 2 / (w * w).sum()) if w.size else 0.0


def threshold(r, w, p):
    """max(--max-residual, --noise-k x robust scatter (1.4826 x MAD) of r)."""
    med = weighted_median(r, w)
    sigma = 1.4826 * weighted_median(np.abs(r - med), w)
    return max(p.max_residual, p.noise_k * sigma)


def elevation_scatter(z, y):
    """
    Robust scatter (m) of a set of lines' water elevations about their
    trend with row, estimated OUT OF SAMPLE: the lines sorted by row, each
    elevation minus the next one's, 1.4826 x MAD about the median
    difference, divided by sqrt 2. Each difference holds two independent
    lines and the median takes out the trend, so this is how far apart in
    elevation two lines at nearly the same row honestly are -- the noise
    the inversion row -> elevation carries there. It is not the misfit of
    a curve fitted to the same lines, which shrinks with every block the
    fit is free to place. Where the waterline moves 30-170 px per metre
    (c1) it is ~0.05-0.15 m; where it barely moves with the tide (c2's far
    field on a crenulated shoreline, a scarp face) lines metres apart in
    elevation share a row and it is a metre or more.
    """
    if z.size < 3:
        return 0.0
    d = np.diff(z[np.argsort(ROW_SIGN * y, kind="stable")])
    return float(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2.0))


def side_reference(t, z, y, sel, before, p):
    """
    A reference for the bottom-up pass: the kept lower-water lines `sel`
    (indices). before=True / False: all of them before / after the line
    judged, weighted exp(-|dt| / --time-scale-days) by closeness to it --
    taken from the set's own nearest line, which on one side in time differs
    from the judged line's weights by one common factor, and the fit
    (weighted medians), the threshold (weighted median and MAD) and the
    effective number are all unchanged by a common factor, so every line
    between the same two reference lines shares one fit. before=None: the
    lines of one tide cycle either side, all weighted alike.
    Returns (fit, threshold_m), or None when the set has too little to judge
    (--min-frames, --min-z-range): keep. The threshold is the larger of
    threshold() and --noise-k x elevation_scatter() of the set's lines.
    """
    if sel.size < 2:
        return None
    tt = t[sel]
    if before is None:
        w = np.ones(sel.size)                    # one tide cycle: all alike
    else:
        w = np.exp(-np.abs(tt - (tt.max() if before else tt.min())) / (p.time_scale_days * DAY))
    zz, yy = z[sel], y[sel]
    if n_eff(w) < p.min_frames or zz.max() - zz.min() < p.min_z_range:
        return None
    fit = fit_column(zz, yy, w)
    if fit is None:
        return None
    rr = residuals(fit, zz, yy, p.extrapolate_flatten)
    return fit, max(threshold(rr, w, p), p.noise_k * elevation_scatter(zz, yy))


def bottom_up(t, z, y, p):
    """
    Judge one column bin's lines -- ALL of the camera's lines with a row in
    the bin, at times t (epoch s), water elevations z, rows y -- in order of
    water level: each layer (--layer-m, 0.25 m) only against the lines at
    LOWER water already kept. A line fails when it sits SEAWARD of where
    lower water put the waterline, by more than the threshold in elevation
    AND by more than --noise-k x --row-noise-px (4 x 2 px) in the photo --
    and that on BOTH sides in time: against the lower-water lines before it
    (back to --window-days) and against those after it. Where one side has
    too few lines to judge, the lower-water lines of its own tide (+/- one
    M2 cycle, 12.42 h, either side) judge it instead; if they are too few
    too, it is kept. Returns the keep mask and each judged line's residual
    (m) and threshold (m), from the reference closest to passing it (NaN
    where nothing could judge it: kept).

    WHY BOTTOM-UP: the two-sided fit takes the majority of the lines at each
    water level as the reference. Where the detector is wrong for MOST of
    the high-tide lines in a column -- the review's model of c1's right
    half, where the search envelope's floor may sit above the high-tide
    waterline, so every line above ~+0.8 m is put somewhere seaward of it
    -- the majority is the error, and it became the reference (reviewer's
    synthetic test, Oct 2026: 0 of ~90 wrong lines failed in c1's column
    bins 64-72). The lines at lower water are not affected, and physics says
    which way the order must run, so they can judge the ones above them
    however many those are.

    WHY BOTH SIDES IN TIME: lower-water lines from the other side of a
    change in the beach lie somewhere else. After the upper beach is cut,
    the later mid-tide lines lie landward of the earlier high-tide lines,
    so the earlier lines look 'seaward of lower water' although each was
    right on its own day. Judged against one +/-3-day window holding both
    beaches, the day before a storm lost its high-tide lines layer by
    layer, and with each layer the reference that would have kept the next
    -- and once that day is 3 days old its window stops changing, so the
    loss was permanent, and dem_change.py rebuilt the pre-storm week from
    it (review, Oct 2026: 0.8 m cut in 6 h on 1 Oct, four whole c2 lines of
    30 Sep 18:30-22:00 dropped; self-test case 4 on the same week: 19% of
    30 Sep's honest points, 7 lines whole). One change moves the lower
    lines on ONE side of a line only: the lines before a pre-storm line are
    pre-storm too, the lines after a post-storm line post-storm, and a line
    wrong on the photo is wrong against both. When one side cannot judge
    (the start of the record, a sparse far-field column), the other side
    alone could again hold the other beach; the lines of the same tide
    cannot, short of a cut inside those hours.

    WHY THE ROW GUARD: where the honest waterline barely moves with the
    tide -- c2's far field on a crenulated shoreline, a scarp face -- the
    lower lines' curve is nearly flat and inverting it turns 1-4 px into
    1-2 m, while the threshold, set from the lower lines' misfit to their
    own curve, stays at 0.75 m: honest higher lines a few pixels below its
    top failed (review: megacusps of 15 m amplitude, 7.2% of c2's honest
    points dropped against origin/main's 1.5%). In pixels the question has
    an honest answer: is the line further seaward of where lower water put
    it (the reference's row at its elevation; beyond the top only that
    bound) than the detector scatters? The scatter is the DETECTOR's, not
    the lower lines' scatter about their curve: in c1's near field (30-200
    px per metre) a few centimetres of water-level noise spread the lower
    lines' rows by 8-57 px (robust sigma about their neighbours in
    elevation, synthetic c1 columns), and a guard on that let up to a
    quarter of c1's wrong high-tide lines through (synthetic test, Oct
    2026: c1few 97% -> 76% caught). The water level's own noise belongs in
    metres, where elevation_scatter() puts it. --row-noise-px 2: c2's bin
    rows differ by 1.8 px (robust sigma) between consecutive slack-water
    frames of the station's 29 Sep - 5 Oct 2026 contours. A line on the
    wrong feature is tens of pixels out.

    One-sided by construction: a line can fail here only by sitting SEAWARD
    of lower water. A line landward of every lower line is past the end of
    their fit, where only a bound is known, and passes. The lowest layer is
    not judged here (nothing lies below it); judge_bin() judges the most
    seaward lines against all the others (seaward_end()).
    """
    n = z.size
    keep = np.ones(n, bool)
    r = np.full(n, np.nan)
    thr = np.full(n, np.nan)
    if n < 3:
        return keep, r, thr
    order = np.argsort(t, kind="stable")
    t, z, y = t[order], z[order], y[order]
    day = np.floor(t / DAY)
    layer = np.floor((z - z.min()) / p.layer_m).astype(np.int64)
    for L in np.unique(layer)[1:]:
        ref = np.flatnonzero(keep & (layer < L))       # in time order
        if ref.size < 2:
            continue
        tr = t[ref]
        sides_of = {}
        for i in np.flatnonzero(layer == L):
            a = int(np.searchsorted(tr, t[i], "left"))
            b = int(np.searchsorted(tr, t[i], "right"))
            spans = ((int(np.searchsorted(tr, (day[i] - p.window_days) * DAY, "left")), a, True),
                     (b, int(np.searchsorted(tr, (day[i] + p.window_days + 1) * DAY, "left")), False))
            sides = []
            for key in spans:
                if key not in sides_of:
                    sides_of[key] = side_reference(t, z, y, ref[key[0]:key[1]], key[2], p)
                sides.append(sides_of[key])
            if sides[0] is None or sides[1] is None:
                # not both: the lower lines of the same tide, either side
                key = (int(np.searchsorted(tr, t[i] - TIDE_CYCLE, "left")),
                       int(np.searchsorted(tr, t[i] + TIDE_CYCLE, "right")), None)
                if key not in sides_of:
                    sides_of[key] = side_reference(t, z, y, ref[key[0]:key[1]], None, p)
                sides = [sides_of[key]]
            verdicts = []
            for side in sides:
                if side is None:
                    continue
                fit, t_m = side
                ri = float(residuals(fit, z[i:i + 1], y[i:i + 1], p.extrapolate_flatten)[0])
                gap = float(np.interp(z[i], fit[2], fit[1]) - fit[0] * y[i])
                verdicts.append((ri - t_m, ri, t_m, ri > t_m and gap > p.noise_k * p.row_noise_px))
            if not verdicts:
                continue                               # nothing can judge it: keep
            _, r[i], thr[i], _ = min(verdicts)
            keep[i] = not all(v[3] for v in verdicts)
    back = np.empty(n, np.int64)
    back[order] = np.arange(n)
    return keep[back], r[back], thr[back]


def seaward_end(z, y, w, keep, r, thr, p):
    """
    Judge the bin's most seaward lines one at a time against all the OTHER
    kept lines, from the most seaward inwards, stopping at the first that
    passes. Returns (residuals with those judgements in, keep).

    WHY: the most seaward line of a bin is the end of the monotone fit, and a
    line at the end of a fit cannot be contradicted by it. A stray far out
    on the water at the LOWEST water level of the window is consistent with
    the order -- lowest water, most seaward row -- so the fit makes it the
    bottom knot itself, its residual is 0, and the rule for lines seaward of
    every other line (residuals()) never applies to it. Strays happen at low
    water, so the lowest frame is a likely place for one (station, 2 Oct 2026
    15:00, c2: the -0.70 m stray is dropped, but relabelled to -0.89 m, the
    window's lowest, it was kept whole and the 7-day DEM gained 24 m of
    empty cells out on the water). Judged against the others it is beyond
    their most seaward line, where the elevation must keep falling at no less
    than half the fit's rate: a line tens of metres out on the water cannot
    be within centimetres of the lowest water level. The other low-water
    lines stay in the reference -- they are the best evidence against it.
    An honest lowest-tide line is a few pixels beyond the others and passes
    the same rule, which ends the peeling: usually one extra fit a bin.
    """
    keep = keep.copy()
    for i in np.argsort(ROW_SIGN * y, kind="stable"):
        if not keep[i]:
            continue
        rest = keep.copy()
        rest[i] = False
        if rest.sum() < 2 or n_eff(w[rest]) < p.min_frames \
                or z[rest].max() - z[rest].min() < p.min_z_range:
            break
        fit = fit_column(z[rest], y[rest], w[rest])
        if fit is None:
            break
        ri = float(residuals(fit, z[i:i + 1], y[i:i + 1], p.extrapolate_flatten)[0])
        if abs(ri) > abs(r[i]):
            r[i] = ri
        if abs(r[i]) <= thr:
            break
        keep[i] = False
    return r, keep


def judge_bin(z, y, w, p, bottom=None):
    """
    Judge one column bin's lines (rows y, water elevations z, weights w).

    1. Order. Rank correlation of row and elevation, in the physical
       direction. At or below -(--reversed-corr) the bin is REVERSED:
       rows run AGAINST the water level, which no beach does -- most of
       its lines are wrong. Below --sweep-corr the two-sided fit cannot be
       trusted to pick the reference, so step 2 runs first.
    2. Bottom-up (bottom_up(), computed once per bin for all of the
       camera's lines; `bottom` returns its (keep, residual, threshold) for
       these lines): each line against the lines at lower water already
       kept, before and after it in time. Its keep set seeds step 3.
    3. Fit the kept lines, judge every line two-sided against the fit
       (residuals()), set aside the ones that fail, refit without them and
       judge every line against that (one robust iteration). Without the
       refit, wrong lines at the very top or bottom of the elevation range
       -- where a block has few members -- can carry the end of the fit
       with them and make the right lines next to them look wrong.
    3b. The most seaward lines are judged against all the others, one at a
       time from the seaward end (seaward_end()), so a stray at the lowest
       water level cannot be the fit's own bottom knot; if that drops one,
       the fit is made again without it.
    4. The bin counts as judged if the kept lines show row rising with
       elevation (rank correlation >= --min-rank-corr). A REVERSED bin
       that does not get there is still judged by step 2 alone: the lines
       that sit seaward of lower water are dropped, the rest kept -- never
       all kept silently.

    Returns (fit, residuals, threshold, reversed, rank_corr), or None when
    the bin cannot be judged (everything in it is kept). `threshold` is a
    scalar or, for a REVERSED bin judged by step 2 alone, per line.
    """
    rc = rank_corr(z, ROW_SIGN * y)
    rev = rc <= -p.reversed_corr
    seed = np.ones(z.size, bool)
    bu = None
    if rc < p.sweep_corr and bottom is not None:
        bu = bottom()
        seed = bu[0]
    out = None
    if seed.sum() >= 2:
        fit = fit_column(z[seed], y[seed], w[seed])
        if fit is not None:
            r = residuals(fit, z, y, p.extrapolate_flatten)
            thr = threshold(r[seed], w[seed], p)
            keep = np.abs(r) <= thr
            if not keep.all() and keep.sum() >= p.min_frames and n_eff(w[keep]) >= p.min_frames:
                fit2 = fit_column(z[keep], y[keep], w[keep])
                if fit2 is not None:
                    fit = fit2
                    r = residuals(fit, z, y, p.extrapolate_flatten)
                    thr = threshold(r[keep], w[keep], p)
                    keep = np.abs(r) <= thr
            was = keep
            r, keep = seaward_end(z, y, w, keep, r.copy(), thr, p)
            peeled = was & ~keep
            if peeled.any() and keep.sum() >= 2:
                fit3 = fit_column(z[keep], y[keep], w[keep])
                if fit3 is not None:
                    fit = fit3
                    r = np.where(peeled, r, residuals(fit, z, y, p.extrapolate_flatten))
                    keep = np.abs(r) <= thr
            if keep.sum() >= 2 and rank_corr(z[keep], ROW_SIGN * y[keep]) >= p.min_rank_corr:
                out = (fit, r, thr, rev, rc)
    if out is None and rev and bu is not None and np.isfinite(bu[1]).any():
        keep_bu, r_bu, thr_bu = bu
        judged = np.isfinite(r_bu)
        # lines the bottom-up pass could not reach (the lowest layers) pass,
        # and so do lines it kept (on one side in time, or by the row guard)
        r_out = np.where(judged, r_bu, 0.0)
        t_out = np.where(judged & ~keep_bu, thr_bu, np.inf)
        fit = fit_column(z[keep_bu], y[keep_bu], w[keep_bu]) if keep_bu.sum() >= 2 else None
        out = (fit, r_out, t_out, rev, rc)
    return out


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
    Returns per frame x bin: residual (m), threshold (m), the bin rows, the
    window fits {(day, bin): (fit, frames, weights, threshold)}, and the
    status of every (day, bin) the day's lines reach:
    {(day, bin): (status, rank_corr)}, status one of 'judged',
    'reversed' (judged, rows ran against the water level), 'few' (too few
    lines or too little elevation range: kept), 'unclear' (no clear
    relation even after setting the out-of-order lines aside: kept).
    """
    R, n_bins = bin_rows(frame_local, col, row, n_frames, p.column_bin)
    resid = np.full((n_frames, n_bins), np.nan)
    thresh = np.full((n_frames, n_bins), np.nan)
    days = np.floor(epochs / DAY)
    good_frame = np.isfinite(epochs) & np.isfinite(elevs)
    fits = {}
    status = {}
    bottom = {}

    def bottom_for(b):
        """bottom_up() of every line with a row in bin b, computed on first
        need (only bins whose order is unclear in some window need it), as
        arrays over all frames."""
        if b not in bottom:
            ii = np.flatnonzero(good_frame & np.isfinite(R[:, b]))
            k, r, t = bottom_up(epochs[ii], elevs[ii], R[ii, b], p)
            K = np.ones(n_frames, bool); Rr = np.full(n_frames, np.nan); T = np.full(n_frames, np.nan)
            K[ii], Rr[ii], T[ii] = k, r, t
            bottom[b] = (K, Rr, T)
        return bottom[b]

    for d in np.unique(days[good_frame]):
        judged = np.where(good_frame & (days == d))[0]
        window = np.where(good_frame & (np.abs(days - d) <= p.window_days))[0]
        tc = np.median(epochs[judged])
        w_all = np.exp(-np.abs(epochs[window] - tc) / (p.time_scale_days * DAY))
        bins = np.where(np.isfinite(R[judged]).any(axis=0))[0]
        for b in bins:
            have = np.isfinite(R[window, b])
            idx = window[have]
            w = w_all[have]
            z, y = elevs[idx], R[idx, b]
            if idx.size < p.min_frames or n_eff(w) < p.min_frames \
                    or z.max() - z.min() < p.min_z_range:
                status[(d, b)] = ("few", np.nan)
                continue
            judged_bin = judge_bin(z, y, w, p,
                                   bottom=lambda: tuple(a[idx] for a in bottom_for(b)))
            if judged_bin is None:
                status[(d, b)] = ("unclear", rank_corr(z, ROW_SIGN * y))
                continue
            fit, r_all, thr, rev, rc = judged_bin
            status[(d, b)] = ("reversed" if rev else "judged", rc)
            mine = np.isin(idx, judged)
            resid[idx[mine], b] = r_all[mine]
            thresh[idx[mine], b] = thr[mine] if np.ndim(thr) else thr
            if fit is not None:
                if np.ndim(thr):
                    fin = thr[np.isfinite(thr)]
                    t_pts = float(np.median(fin)) if fin.size else p.max_residual
                else:
                    t_pts = thr
                fits[(d, b)] = (fit, idx, w, t_pts)
    return resid, thresh, R, fits, status


def decide(resid, thresh, p):
    """Failing bins per frame, and the frames to drop whole."""
    judged = np.isfinite(resid)
    with np.errstate(invalid="ignore"):            # NaN = not judged (numpy < 1.18 warns)
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


def _col_ranges(bins, bin_px):
    """Consecutive bins as 'a-b' pixel-column ranges."""
    bins = sorted(set(int(b) for b in bins))
    out, start = [], None
    for i, b in enumerate(bins):
        if start is None:
            start = b
        if i == len(bins) - 1 or bins[i + 1] != b + 1:
            out.append((start * bin_px, (b + 1) * bin_px - 1))
            start = None
    return out


def bin_summary(cam, status, p):
    """
    Counts of column bins x days by status, for the CONSISTENCY line, and
    WARNING lines for REVERSED bins: rows running against the water level
    there mean the detector is wrong for most lines, which is for a person
    to look at (search envelope? camera moved?) even though the filter
    drops what it can.
    """
    kinds = [v[0] for v in status.values()]
    n = len(kinds)
    n_j = sum(k in ("judged", "reversed") for k in kinds)
    n_rev = kinds.count("reversed")
    msg = (f"column bins judged {n_j:,} of {n:,} ({p.column_bin} px x day; "
           f"{kinds.count('few'):,} too few lines, {kinds.count('unclear'):,} unclear: kept)")
    if n_rev:
        msg += f", {n_rev:,} REVERSED"
    warn = []
    rev = [(d, b, rc) for (d, b), (k, rc) in status.items() if k == "reversed"]
    if rev:
        days_all = sorted({d for (d, _) in status})
        days_rev = sorted({d for d, _, _ in rev})
        last = datetime.fromtimestamp(days_rev[-1] * DAY, tz=timezone.utc).date()
        first = datetime.fromtimestamp(days_rev[0] * DAY, tz=timezone.utc).date()
        rcs = np.array([rc for _, _, rc in rev])
        ranges = ", ".join(f"{a}-{b}" for a, b in _col_ranges([b for _, b, _ in rev], p.column_bin))
        recent = [b for d, b, _ in rev if d > days_all[-1] - p.plot_days]
        warn.append(f"WARNING: camera {cam}: waterline rows run AGAINST the water level in pixel "
                    f"columns {ranges}")
        warn.append(f"         ({len(rev)} column bin-days on {len(days_rev)} of {len(days_all)} days, "
                    f"{first} to {last}; rank correlation {rcs.min():+.2f} to {rcs.max():+.2f}).")
        if recent:
            warn.append("         In the last {} days: columns {}.".format(
                p.plot_days, ", ".join(f"{a}-{b}" for a, b in _col_ranges(recent, p.column_bin))))
        warn.append("         Higher water must put the line LOWER in the photo (nearer the camera);")
        warn.append("         here the high-tide lines sit seaward of the low-tide ones, so most of")
        warn.append("         them are wrong -- the detector cannot reach the real line (search")
        warn.append("         envelope floor above the high-tide waterline?). Lines seaward of")
        warn.append("         lower-water lines there are dropped; look at the diagnostic plot.")
    return msg, warn


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
        resid, thresh, R, fits, status = judge_camera(epochs, elevs, fl, col, row, frames.size, args)
        fail, whole, n_judged, n_fail = decide(resid, thresh, args)
        bins_msg, warn = bin_summary(cam, status, args)
        for line in warn:
            print(line)

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
        fday = np.floor(epochs / DAY)
        for f in affected:
            n_rev = int(sum(status.get((fday[f], b), ("",))[0] == "reversed"
                            for b in np.flatnonzero(fail[f])))
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
            if n_rev:
                reason += (f"; {n_rev} of the failing bins REVERSED (rows run against the water "
                           f"level there: judged against the lower-water lines only)")
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
                bins_judged=int(n_judged[f]), bins_failed=int(n_fail[f]), bins_reversed=n_rev,
                pixel_column_from=round(float(cmin[f]), 1) if np.isfinite(cmin[f]) else "",
                pixel_column_to=round(float(cmax[f]), 1) if np.isfinite(cmax[f]) else "",
                median_residual_m=round(med, 3), max_residual_m=round(worst, 3),
                threshold_m=round(thr, 3) if np.isfinite(thr) else "",
                reason=reason))
        msg = (f"{cam}: {n_whole} line(s) dropped, {affected.size - n_whole} trimmed, "
               f"{n_drop:,} of {n_cam_pts:,} points ({100.0 * n_drop / max(n_cam_pts, 1):.1f}%); "
               f"{bins_msg}")
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
                              status=status, epochs=epochs, elevs=elevs, applied=applied)
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
                 "bins_reversed", "pixel_column_from", "pixel_column_to", "median_residual_m", "max_residual_m",
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
    # Column bins whose rows ran AGAINST the water level in the window: a
    # red bar along the top of the photo, so where the detector was wrong
    # for most lines is visible even where the dropped lines hide it.
    rev = [(d, b) for (d, b), (k, _) in info.get("status", {}).items()
           if k == "reversed" and last_day - args.plot_days < d <= last_day]
    rev_note = ""
    if rev:
        n_days = len({d for d, _ in rev})
        for a, b in _col_ranges([b for _, b in rev], args.column_bin):
            ax.axvspan(a, b + 1, ymin=0.965, ymax=1.0, color="red", alpha=0.75, lw=0)
            ax.text((a + b) / 2, 0.012 * height, "REVERSED", color="white", fontsize=7.5,
                    fontweight="bold", ha="center", va="top")
        rev_note = (f"\nred bar: columns where rows ran AGAINST the water level on {n_days} "
                    f"day(s) -- most lines there wrong (search envelope?)")
    ax.set_xlim(0, width); ax.set_ylim(height, 0); ax.axis("off")
    status = "" if info["applied"] else "   [NOT APPLIED: would drop too much, see log]"
    ax.set_title(f"{cam.upper()}  {first} to {last}: {len(listed)} line(s) with parts dropped as out of "
                 f"order (thick; dashed = whole line){status}\nthin white = kept lines, {bg_note}"
                 f"{rev_note}", fontsize=10)
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
            kind, rc = info.get("status", {}).get((day, b), ("", np.nan))
            ax2.set_title(f"columns {b * args.column_bin}-{(b + 1) * args.column_bin - 1} px: the fit "
                          f"that judged {dd} (window +/-{args.window_days:g} d)"
                          + (f"\nREVERSED: rank correlation {rc:+.2f}, judged against lower water"
                             if kind == "reversed" else ""), fontsize=9)
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


# ---------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------

def self_test(cal_dir, keep_dir=None):
    """
    python3 waterline_consistency.py --self-test

    Synthetic waterlines through the station's own calibration
    (CACO05_<cam>_20240801_IO.yaml + _20251113_EO-CV.yaml), so it checks
    this filter against the cameras it runs on. Shore-parallel beach with a
    berm, foreshore 1:8 and low-tide terrace 1:30, 30-min frames 11:00-22:30
    UTC, M2 tide unless said otherwise.
      1. PHYSICS: with 5 m alongshore undulation and a 0.4 m erosion of the
         upper beach, the row rises with water elevation in every column of
         c1 and c2 (ROW_SIGN).
      2. SYSTEMATIC ERROR, the review's model of the 29 Sep - 5 Oct 2026
         failure: on c1 EVERY line above +0.8 m is redrawn from column 1400
         on along the -0.4 m contour (seaward of the low-tide lines, inside
         the search band), as a search-envelope floor above the high-tide
         waterline would force. At least 90% of those points must be
         dropped, and c1 columns >= 1400 reported REVERSED.
      3. REAL CHANGE: the 0.4 m erosion of the upper beach on day 4 must
         not be rejected -- the honest lines (all of c2, c1 below +0.8 m)
         lose less than 0.5% of their points.
      4. STORM CUT (review, Oct 2026), on the week the review's generator
         used: the Signature 1000's water level and waves of 16-25 Dec 2024
         (sig1000_waves_ALL.csv) moved to 26 Sep - 5 Oct 2026, so the
         Hs 3.6 m storm lands on 1 Oct and the tide goes from springs to
         neaps across it. The upper beach is LOWERED 0.8 m over 6 h from
         1 Oct 00:00 (the full 0.8 m above +1.0 m, nothing below -0.2 m);
         every line carries its wave setup (0.037 sqrt(Hs L0): ~0.3 m on
         ordinary days, up to 0.74 m in the storm) in its position and in
         beach_elevation_navd88, as the cron writes them. The bottom-up
         pass runs in every bin (the station's lines are untidy enough that
         it runs in ~90% of them). The lines of 30 Sep, the day BEFORE the
         cut, must lose less than 0.5% of their points and none whole
         (HEAD of 8 Oct 2026: 19%, 7 lines whole; origin/main: none), and
         all honest lines less than 0.5%. In the same data a c2 line at the
         LOWEST water level of its window leaves the waterline on the right
         of the photo for a line 40 m out on the water: at least 90% of its
         points 38 m or more out must be dropped (it was the fit's own
         bottom knot and was kept whole).
      5. MEGACUSPS (review, Oct 2026): c2 on a shoreline with 15 m
         crenulations every 140 m, each waterline drawn where the camera
         sees it -- the most seaward edge in each column, which behind a
         horn stops moving with the tide. No defects; honest points dropped
         must stay below 5% (origin/main 3.9%; 7.3% when the bottom-up pass
         judged the higher lines' few-pixel offsets from the lower lines in
         metres).
    Returns 0 if all pass, 1 otherwise.
    """
    import tempfile
    import contextlib
    import io as _io
    from georectify import load_intrinsics, load_extrinsics, pixel_to_ground
    from view_reproject import ground_to_pixel

    W, H = 2448, 2048
    crop = {"c1": (532, 1761), "c2": (61, 1024)}
    io = {c: load_intrinsics(Path(cal_dir) / f"CACO05_{c}_20240801_IO.yaml") for c in ("c1", "c2")}
    eo = {c: load_extrinsics(Path(cal_dir) / f"CACO05_{c}_20251113_EO-CV.yaml") for c in ("c1", "c2")}
    # Shore frame of the surveyed beach: alongshore bearing 345.8 deg (N-NW),
    # seaward 75.8 deg, origin on the wrack line in front of the cameras.
    r0 = np.array([420107.0, 4638322.0])
    t_al = np.array([np.sin(np.radians(345.8)), np.cos(np.radians(345.8))])
    t_cs = np.array([np.sin(np.radians(75.8)), np.cos(np.radians(75.8))])
    t0 = datetime(2026, 9, 28, tzinfo=timezone.utc).timestamp()

    def contour_s(h, a, t, beach):
        """Cross-shore position (m seaward of the wrack line) of the h contour:
        berm +3.8 m to s = 12, foreshore 1:8 to 0 m, terrace 1:30 below.
        beach['cut04'] = t: the upper beach (above ~+0.3 m) cut 0.4 m from t;
        beach['storm'] = (t, seconds, m): the profile LOWERED by up to m,
        ramped in over that time -- the full m above +1.0 m, nothing below
        -0.2 m, linear between -- so the h contour after the cut is where
        the old beach stood at h0, h0 - lowering(h0) = h (the mid-beach
        flattens; the review's generator, Oct 2026); beach['cusp'] = (m, m):
        alongshore undulation amplitude and wavelength (default 5, 150)."""
        if "storm" in beach:
            tc, dur, mag = beach["storm"]
            f = float(np.clip((t - tc) / dur, 0, 1))
            m = mag * f * f * (3 - 2 * f)
            if h > -0.2 and m > 0:
                h = h + m if h > 1.0 - m else (h + 0.2 * m / 1.2) / (1 - m / 1.2)
        if h >= 0:
            s = 12 + (3.75 - h) / 0.125
        else:
            s = 12 + 3.75 / 0.125 - h * 30
        if "cut04" in beach:
            s = s - 0.4 * np.clip((h - 0.3) / 0.7, 0, 1) * (t >= beach["cut04"]) / 0.125
        amp, lam = beach.get("cusp", (5.0, 150.0))
        return s + amp * np.sin(2 * np.pi * a / lam)

    a_grid = np.arange(-250.0, 350.0, 0.5)
    a_fine = np.arange(-450.0, 700.0, 0.2)
    cols = np.arange(0, W, 4)

    def line(cam, h, t, beach, seaward_edge=False, out=None):
        """Rows of the h waterline at `cols`. seaward_edge: where the contour
        folds in the photo (a cusp horn in front of its bay) the camera sees
        the most seaward edge in each column, as the detector does. out(a):
        metres further seaward than the waterline (a line on the water)."""
        a = a_fine if seaward_edge else a_grid
        s = contour_s(h, a, t, beach)
        if out is not None:
            s = s + out(a)
        E = r0[0] + a * t_al[0] + s * t_cs[0]
        N = r0[1] + a * t_al[1] + s * t_cs[1]
        u, v, ok = ground_to_pixel(E, N, np.full(E.size, h), io[cam], eo[cam])
        ok &= (u >= 0) & (u < W) & (v >= 0) & (v < H)
        Eb, Nb = pixel_to_ground(np.where(ok, u, W / 2), np.where(ok, v, H / 2), h, io[cam], eo[cam])
        ok &= np.hypot(Eb - E, Nb - N) < 0.5
        if ok.sum() < 2:
            return np.full(cols.size, np.nan)
        if not seaward_edge:
            o = np.argsort(u[ok])
            return np.interp(cols, u[ok][o], v[ok][o], left=np.nan, right=np.nan)
        best = np.full(cols.size, np.inf)
        np.minimum.at(best, np.clip(np.round(u[ok] / 4).astype(int), 0, cols.size - 1), v[ok])
        has = np.flatnonzero(np.isfinite(best))
        rows = np.interp(cols, cols[has], best[has], left=np.nan, right=np.nan)
        for i, k in zip(has[:-1], has[1:]):
            if k - i > 6:                                  # no bridging over 24 px
                rows[i + 1:k] = np.nan
        return rows

    def tide(ep):
        return 0.1 + 1.3 * np.sin(2 * np.pi * (ep - t0) / (12.42 * 3600))

    tmp = Path(keep_dir) if keep_dir else Path(tempfile.mkdtemp(prefix="wc_selftest_"))
    tmp.mkdir(parents=True, exist_ok=True)

    def run_case(name, frames, extra=()):
        """Write the frames [(name, cam, ep, z_file, cols, rows, flags, tide), ...]
        as a contour file, run the filter (options `extra`), return
        {name: dropped mask}, the run's stdout, summary, plot info and args."""
        src = tmp / f"contour_points_{name}.csv"
        with open(src, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["source_file", "camera", "capture_time_utc", "capture_epoch", "pixel_column",
                        "pixel_row", "tide_elevation_navd88", "beach_elevation_navd88"])
            for fname, cam, ep, z, c, r, _, h in frames:
                cap = datetime.fromtimestamp(ep, tz=timezone.utc).isoformat()
                for cc, rr in zip(c, r):
                    w.writerow([fname, cam, cap, int(ep), int(cc), f"{rr:.1f}", f"{h:.4f}", f"{z:.4f}"])
        args = build_parser().parse_args([str(src)] + list(extra))
        buf = _io.StringIO()
        with contextlib.redirect_stdout(buf):
            data, keep, report, summary, plot_info = run(args)
        frame_of = {n: i for i, n in enumerate(data["names"])}
        dropped = {fname: ~keep[data["frame"] == frame_of[fname]] for fname, *_ in frames}
        return dropped, buf.getvalue(), summary, plot_info, args

    def keep_rows(cam, r, c):
        with np.errstate(invalid="ignore"):
            k = np.isfinite(r) & (r >= crop[cam][0]) & (r <= crop[cam][1])
        return c[k], r[k], k

    failures = []
    # 1. physics
    for cam in ("c1", "c2"):
        for t in (t0, t0 + 4.6 * DAY):
            rows = np.array([line(cam, h, t, {"cut04": t0 + 3.6 * DAY}) for h in np.arange(-1.4, 2.01, 0.2)])
            d = np.diff(rows, axis=0)
            d = d[np.isfinite(d)]
            share = float(np.mean(d * ROW_SIGN > 0)) if d.size else 0.0
            print(f"  physics {cam}: row rises with water level in {share:.2%} of {d.size} column steps")
            if share < 0.999:
                failures.append(f"physics {cam}: only {share:.1%} of column steps rise")

    # 2-3. systematic c1 error, 0.4 m erosion on day 3.6
    rng = np.random.default_rng(5)
    beach = {"cut04": t0 + 3.6 * DAY}
    frames = []
    for d in range(8):
        for k in range(24):
            ep = t0 + d * DAY + 11 * 3600 + k * 1800
            h = tide(ep)
            for cam in ("c1", "c2"):
                r = line(cam, h + rng.normal(0, 0.05), ep, beach) + rng.normal(0, 1.2, cols.size) \
                    + rng.normal(0, 1.5)
                wrong = np.zeros(cols.size, bool)
                if cam == "c1" and h > 0.8:
                    rw = line(cam, -0.4 + rng.normal(0, 0.1), ep, beach) + rng.normal(0, 1.5, cols.size)
                    wrong = (cols >= 1400) & np.isfinite(rw)
                    r = np.where(wrong, rw, r)
                c, rr, k_ = keep_rows(cam, r, cols)
                if c.size >= 20:
                    frames.append((f"{int(ep)}.{cam}.selftest", cam, ep, h, c, rr, wrong[k_], h))
    dropped, out, summary, plot_info, args = run_case("systematic", frames)
    wrong_n = wrong_drop = honest_n = honest_drop = 0
    for fname, cam, ep, z, c, r, wr, h in frames:
        dr = dropped[fname]
        wrong_n += int(wr.sum()); wrong_drop += int((dr & wr).sum())
        if cam == "c2" or h <= 0.8:
            honest_n += dr.size; honest_drop += int(dr.sum())
    catch = wrong_drop / max(wrong_n, 1)
    false = honest_drop / max(honest_n, 1)
    print(f"  systematic c1 error: {wrong_drop:,} of {wrong_n:,} wrong points dropped ({catch:.1%})")
    print(f"  real change (0.4 m erosion): {honest_drop:,} of {honest_n:,} honest points dropped "
          f"({false:.2%})")
    rev_ok = any(l.startswith("WARNING: camera c1: waterline rows run AGAINST") for l in out.splitlines())
    st = plot_info["c1"]["status"]
    rev_cols = sorted({b * args.column_bin for (d, b), (k, _) in st.items() if k == "reversed"})
    print(f"  REVERSED warning for c1: {'yes' if rev_ok else 'NO'}"
          + (f", columns {rev_cols[0]}-{rev_cols[-1] + args.column_bin - 1}" if rev_cols else ""))
    for line_ in summary:
        print(f"  CONSISTENCY {line_}")
    if catch < 0.90:
        failures.append(f"only {catch:.1%} of the systematic wrong points dropped (need 90%)")
    if false >= 0.005:
        failures.append(f"{false:.2%} of honest points dropped (limit 0.5%)")
    if not rev_ok or not rev_cols or rev_cols[0] < 1400 - 3 * args.column_bin:
        failures.append("c1 columns >= 1400 not reported REVERSED (or reported far left of them)")

    # 4. A storm cut on the week the review's generator used: the Signature
    #    1000's water level and waves from 16 Dec 2024 (sig1000_waves_ALL.csv,
    #    +0.09 m to NAVD88), shifted to start 26 Sep 2026, so the Hs 3.6 m
    #    storm lands on 1 Oct and the tide goes from springs to neaps across
    #    it. The upper beach is lowered 0.8 m from 1 Oct 00:00 over 6 h. Every
    #    line carries its wave setup, 0.037 sqrt(Hs L0) -- ~0.3 m on ordinary
    #    days, up to 0.74 m in the storm -- in its position and in
    #    beach_elevation_navd88, with 0.05 m + 0.03 Hs of scatter, and storm
    #    waves widen the detection scatter (x(1 + (Hs - 1.5) / 2)). c2 keeps
    #    60% of its frames, as on the station.
    rng = np.random.default_rng(7)
    t_start = datetime(2026, 9, 26, tzinfo=timezone.utc).timestamp()
    shift = t_start - datetime(2024, 12, 16, tzinfo=timezone.utc).timestamp()
    sig_t, sig_wl, sig_hs, sig_tp = [], [], [], []
    sig_path = Path(__file__).resolve().parent / "sig1000_waves_ALL.csv"
    if sig_path.exists():
        with open(sig_path, newline="") as f:
            for r in csv.DictReader(f):
                sig_t.append(datetime.fromisoformat(r["time"]).replace(tzinfo=timezone.utc).timestamp()
                             + shift)
                sig_wl.append(float(r["water_level"])); sig_hs.append(float(r["wh_4061"]))
                sig_tp.append(float(r["wp_peak"]))
    else:
        print(f"  storm cut: SKIPPED, {sig_path.name} not found next to this script")
    t_cut = datetime(2026, 10, 1, tzinfo=timezone.utc).timestamp()
    beach = {"storm": (t_cut, 6 * 3600.0, 0.8)}
    frames = []
    for d in range(10 if sig_t else 0):
        for k in range(24):
            ep = t_start + d * DAY + 11 * 3600 + k * 1800
            h = float(np.interp(ep, sig_t, sig_wl)) + 0.09
            hs, tp = float(np.interp(ep, sig_t, sig_hs)), float(np.interp(ep, sig_t, sig_tp))
            setup = 0.037 * np.sqrt(hs * 9.81 * tp * tp / (2 * np.pi))
            noise = 1.0 + 0.5 * max(0.0, hs - 1.5)
            for cam in ("c1", "c2"):
                if cam == "c2" and rng.random() > 0.6:
                    continue
                z = h + setup
                r = line(cam, z + rng.normal(0, 0.05 + 0.03 * hs), ep, beach) \
                    + noise * (rng.normal(0, 1.2, cols.size) + rng.normal(0, 1.5))
                c, rr, _ = keep_rows(cam, r, cols)
                if c.size >= 20:
                    frames.append([f"{int(ep)}.{cam}.storm", cam, ep, z, c, rr, None, h])
    if not frames:
        frames = None
    else:
        c2 = [f for f in frames if f[1] == "c2" and 6 <= (f[2] - t_start) / DAY < 8]
        stray = min(c2, key=lambda f: f[3])
        near = [f[3] for f in frames if f[1] == "c2" and abs(f[2] - stray[2]) <= 3.5 * DAY
                and f is not stray]
        # Like the station's stray of 2 Oct 2026 15:00 (c2): on the right of
        # the photo it leaves the waterline for a line on the water, here
        # 40 m seaward of it from 200 m alongshore of the cameras towards the
        # seam (ramping in over 60 m); and it is at the lowest water level of
        # its window.
        def on_water(a):
            return 40.0 * np.clip((260.0 - a) / 60.0, 0, 1)
        dev = line("c2", stray[3], stray[2], beach) - line("c2", stray[3], stray[2], beach, out=on_water)
        dc = np.interp(stray[4], cols, np.nan_to_num(dev, nan=0.0))
        stray[5] = stray[5] - dc
        stray[6] = (dc >= 0.95 * np.nanmax(dev), dc >= 0.25 * np.nanmax(dev))
        stray[3] = min(near) - 0.01
        frames = [tuple(f) for f in frames]
    if frames:
        # The station's lines are far less tidy than these (rank correlation
        # 0.4-0.6 a column, 29 Sep - 5 Oct 2026) and the bottom-up pass runs in
        # ~90% of its bins; here it is made to run in all of them.
        dropped, out, summary, plot_info, args = run_case("storm", frames, ["--sweep-corr", "1.01"])
        day_before = int((t_cut - t_start) // DAY) - 1
        n_b = d_b = n_h = d_h = whole_b = 0
        for fname, cam, ep, z, c, r, flags, h in frames:
            dr = dropped[fname]
            if flags is not None:
                far, ramp = flags
                out_n, out_d = int(far.sum()), int((dr & far).sum())
                ramp_n, ramp_d = int(ramp.sum()), int((dr & ramp).sum())
                continue
            n_h += dr.size; d_h += int(dr.sum())
            if int((ep - t_start) // DAY) == day_before:
                n_b += dr.size; d_b += int(dr.sum()); whole_b += int(dr.all())
        fb, fh = d_b / max(n_b, 1), d_h / max(n_h, 1)
        catch_s = out_d / max(out_n, 1)
        print(f"  storm cut 1 Oct (0.8 m in 6 h, setup +0.15 to +0.74 m): 30 Sep, the day before, "
              f"{d_b:,} of {n_b:,} honest points dropped ({fb:.2%}), {whole_b} line(s) whole; all "
              f"honest lines {d_h:,} of {n_h:,} ({fh:.2%})")
        print(f"  stray at the window's lowest level ({stray[3]:+.2f} m, c2): {out_d:,} of {out_n:,} "
              f"points 38-40 m out on the water dropped ({catch_s:.1%}); {ramp_d:,} of {ramp_n:,} "
              f"10 m or more out")
        for line_ in summary:
            print(f"  CONSISTENCY {line_}")
        if fb >= 0.005 or whole_b:
            failures.append(f"storm: day before the cut lost {fb:.2%} of its honest points, "
                            f"{whole_b} line(s) whole (limit 0.5%, none)")
        if fh >= 0.005:
            failures.append(f"storm: {fh:.2%} of honest points dropped (limit 0.5%)")
        if catch_s < 0.90:
            failures.append(f"storm: only {catch_s:.1%} of the lowest-level stray's points dropped (need 90%)")

    # 5. megacusps, c2
    rng = np.random.default_rng(11)
    beach = {"cusp": (15.0, 140.0)}
    frames = []
    for d in range(8):
        for k in range(24):
            ep = t0 + d * DAY + 11 * 3600 + k * 1800
            h = tide(ep)
            r = line("c2", h + rng.normal(0, 0.05), ep, beach, seaward_edge=True) \
                + rng.normal(0, 1.2, cols.size) + rng.normal(0, 1.5)
            c, rr, _ = keep_rows("c2", r, cols)
            if c.size >= 20:
                frames.append((f"{int(ep)}.c2.cusp", "c2", ep, h, c, rr, None, h))
    dropped, out, summary, plot_info, args = run_case("cusp", frames)
    n_c = sum(v.size for v in dropped.values()); d_c = sum(int(v.sum()) for v in dropped.values())
    fc = d_c / max(n_c, 1)
    print(f"  megacusps (15 m every 140 m), c2: {d_c:,} of {n_c:,} honest points dropped ({fc:.2%})")
    for line_ in summary:
        print(f"  CONSISTENCY {line_}")
    if fc >= 0.05:
        failures.append(f"megacusps: {fc:.2%} of honest points dropped (limit 5%)")

    if not keep_dir:
        for p_ in tmp.iterdir():
            p_.unlink()
        tmp.rmdir()
    print("SELF-TEST " + ("PASSED" if not failures else "FAILED: " + "; ".join(failures)))
    return 0 if not failures else 1


def build_parser():
    ap = argparse.ArgumentParser(
        description="Drop waterlines (or parts) that are out of order with the other lines: "
                    "row must rise with water elevation in every image column.")
    ap.add_argument("contour_csv", nargs="?")
    ap.add_argument("--self-test", action="store_true",
                    help="Run the synthetic self-test through the calibration and exit.")
    ap.add_argument("--calibration", default=str(Path(__file__).resolve().parent / "calibration"),
                    help="Calibration folder for --self-test.")
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
    ap.add_argument("--row-noise-px", type=float, default=2.0,
                    help="The detector's row scatter (px) between lines at the same water level; "
                         "the bottom-up pass drops a line only when it is more than --noise-k x "
                         "this seaward of lower water in the photo (default 2).")
    ap.add_argument("--min-frames", type=float, default=8.0,
                    help="Effective number of lines needed to judge a bin; fewer = keep (default 8).")
    ap.add_argument("--min-z-range", type=float, default=0.5,
                    help="Elevation range (m) the lines must span to judge a bin (default 0.5).")
    ap.add_argument("--min-rank-corr", type=float, default=0.5,
                    help="Rank correlation of row and elevation (physical direction) the kept "
                         "lines need for a bin to be judged two-sided (default 0.5).")
    ap.add_argument("--sweep-corr", type=float, default=0.9,
                    help="Below this rank correlation the lines are first judged bottom-up, "
                         "each against the lines at lower water (default 0.9).")
    ap.add_argument("--reversed-corr", type=float, default=0.2,
                    help="A bin whose rank correlation is at or below minus this is REVERSED: "
                         "warned about, and judged bottom-up (default 0.2).")
    ap.add_argument("--layer-m", type=float, default=0.25,
                    help="Water-level layer (m) of the bottom-up pass (default 0.25).")
    ap.add_argument("--min-judged-bins", type=int, default=4,
                    help="Judged bins a line needs before it can be dropped whole (default 4).")
    ap.add_argument("--frame-fraction", type=float, default=0.5,
                    help="Drop the whole line when more than this share of its judged bins fail "
                         "(default 0.5).")
    ap.add_argument("--max-drop-fraction", type=float, default=0.5,
                    help="If a camera would lose more than this share of its points, pass it "
                         "through unfiltered and warn (default 0.5).")
    return ap


def main():
    ap = build_parser()
    args = ap.parse_args()
    if args.self_test:
        sys.exit(self_test(args.calibration))
    if not args.contour_csv:
        ap.error("contour_csv is required")

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

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
         others does not. A line whose water is HIGHER than the lowest
         line's and yet lies beyond it is out of order already and gets
         no margin: the average rate (residuals()).
     The most seaward line of a bin is the end of the fit, which cannot
     contradict it: a stray out on the water at the window's LOWEST
     water level became the bottom knot itself, residual 0. So the most
     seaward lines are judged one at a time against all the others,
     from the seaward end in, until one passes among them; one that
     passes only beyond their end is held out of the reference while
     the next is judged, so a run of strays (glare or fog over several
     frames) cannot vouch for itself (seaward_end()).
  5. Decide. A bin fails when |residual| exceeds the threshold below
     and the lines on both sides of it in time say so too (5b).
     In a failing bin and the bins either side, each point is then
     judged on its own row against that bin's fit: a wrong segment
     starts part way through a bin, so this drops it up to where it
     leaves the right line and keeps the rest. The whole line is
     dropped only when more than half of its judged bins fail.
  5b. Change, not error (time_test()). A beach that changes puts its
     lines out of order with the beach before it: after a storm cut the
     mid-tide lines lie where the old beach had its high-tide lines, and
     where no post-cut line reaches the levels above them (c1's high
     tides fall below the floor of its search envelope; neaps) the fit
     reads them as the old beach's elevations, a metre landward of their
     level. So each line that fails is judged again against the lines
     BEFORE it and those AFTER it in time, each side its own fit: if
     either agrees with it (within half its threshold), it is kept and
     counted as beach change; it goes only when both sides contradict it,
     or one does and the other cannot say (too few lines there: the
     lines of its own tide judge instead). A line LANDWARD of the beach
     before it is not dropped on the past alone: with nothing after it to
     judge (the last tide of the record) it is kept until there is -- a
     storm cuts a beach in hours, it builds back over days, and the errors
     seen here lie seaward. A line on the wrong feature contradicts the
     lines before and after it alike.

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
      none of them was touched there. It does NOT sit above a storm-
      sized change in general: on the surveyed beach, with c1's high
      tides below the floor of its search envelope after a cut, honest
      post-cut lines are 0.8-1.5 m out of order with the old beach
      (review, Oct 2026) -- as far as the errors. What keeps them is
      step 5b, not the threshold.
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
      bottom-up pass and the seaward end drop a line only when it lies more
      than --noise-k times this seaward of where the other lines put its
      level in the photo (bottom_up(), seaward_end()).
  TWIN_M 0.1 m, SEAWARD_RUN 3 (seaward_end()). Slack low water moves less
      than 0.1 m in an hour, so the frames of one glare or fog lie within
      0.1 m of each other; a run of up to 3 such strays is judged together.
  Too little data, at any of these steps, means KEEP. The filter only
  removes what the other lines clearly contradict.

  (Step 5b only ever KEEPS: where it has too little to judge with, the
  window's verdict stands.)

  Real change is protected three ways: the time weighting (above); the
  two sides in time -- a line is dropped only when the lines before it
  and the lines after it both contradict it (step 5b; the bottom-up
  pass asks the same of its lower-water references); and re-judging:
  the cron rebuilds the contour file from the whole archive every run,
  so a line judged in the last tide of the record, which has no lines
  after it yet, is judged again the next day. Nothing is deleted from
  the archive. Synthetic tests on the surveyed beach (review's
  generator, Oct 2026; 26 scenarios: 0.5-1.0 m cuts at the storm peak,
  after it, at spring tides, in the afternoon, in the last hours of the
  record, in a 60-day record; 0.6-1.0 m accretion): no honest c1 line
  dropped whole (the version before: 63) and at most 0.07% of c1's honest
  points in any one (before: 2.5%); self-test case 7 likewise (before: 9
  lines whole at 0.8 m, 3 at 0.5 m).

  Safety: if a camera would lose more than --max-drop-fraction (50%)
  of its points, the filter is NOT applied to that camera (WARNING,
  report rows marked NOT APPLIED): at that point the reference itself
  is suspect -- e.g. a camera knocked out of aim -- and that is for a
  person to look at.

KNOWN LIMIT -- a camera wrong for MOST lines of a column AND a cut in the
  window (review, Oct 2026; not a regression: 566eef7 does the same). In
  c1's right half after a 0.8 m cut, every line above ~0 m left the search
  envelope and was redrawn on the water (the review's env1_ero08): the
  honest lines there span only -0.6..0 m, the wrong ones 0..+1.8 m on the
  same rows. The lowest wrong lines are only 0.6 m out of order -- under the
  threshold -- and those of the record's last days have nothing after them
  to judge them, so they stay; the window's kept lines then show no clear
  order and the bin is 'unclear: kept' (215 of 760 c1 bin-days; 82% of the
  wrong points dropped, against 96% without the cut). Tried, each on the
  review's scenarios: capping the out-of-sample elevation scatter at what
  2 px of row noise make of the reference (82%); not letting a side that
  passes only by that scatter outvote one that confirms (82%); only lines
  that AGREE with lower water as references for higher layers (86%); the
  bottom-up verdicts in 'unclear' bins (87-88%, also with only failures
  64 px or more seaward counted, or only where the lower lines move 60 px
  or more per metre). Every one dropped more honest c2 lines on a
  crenulated shoreline (cusp15 3,178 -> 3,379..4,391 points; c2-only
  4,881 -> 4,995..7,925), where c2's oblique view folds the high-tide lines
  behind the horns seaward of lower ones -- the same signature in the
  photo. Telling the two apart needs the
  viewing geometry: from the calibration, how far a waterline can move per
  metre of water in each column (>= 34 px/m in c1, >= 4.7 px/m in c2's far
  field, on the surveyed beach). Not done here. What it costs: the post-cut
  7-day DEM south of the seam had 16 of 448 cells off by more than 0.5 m and
  51 blanked (unfiltered: 71 and 410; the station's real 29 Sep - 5 Oct
  contours: no c1 line kept below lower water on the ground).

ELEVATION: beach_elevation_navd88 when the contour file has it (water
level + wave setup, --setup-coef), else tide_elevation_navd88 -- the
same choice, row by row, as georectify.py and dem_from_contours.py, so
the filter judges the elevations the DEM will use.

OUTPUTS:
  --output   filtered copy, same columns and row order, rows removed.
  --report   one row per line affected: frame, camera, capture time,
             elevation, points in the line and dropped, bins judged and
             failed, columns dropped, median residual of the failing bins
             (of the sign most of them have; the reason counts each sign
             when they differ) and the largest (m),
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
AGREE = 0.5                          # time_test(): a reference agrees within this x its threshold
SEAWARD_RUN = 3                      # seaward_end(): a run of up to this many strays is judged
TWIN_M = 0.1                         # seaward_end(): lines this close in level lean on each other

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
            try:                     # both or neither: the arrays must stay aligned
                u_, v_ = float(r[i_u]), float(r[i_v])
            except ValueError:
                u_ = v_ = np.nan
            us.append(u_); vs.append(v_)
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
    n = len(values)
    if n == 1:
        return float(values[0])
    if n <= 32:
        # The same result in plain Python: the monotone fits call this a
        # million times a run on blocks of a few lines, where numpy's
        # per-call overhead was most of the filter's time.
        pairs = sorted(zip(values.tolist(), weights.tolist()), key=lambda q: q[0])
        cum = []
        c = 0.0
        for _, wq in pairs:
            c += wq
            cum.append(c)
        half = 0.5 * c
        i = next((k for k, ck in enumerate(cum) if ck >= half), n - 1)
        if i < n - 1 and abs(cum[i] - half) <= 1e-12 * half:
            return 0.5 * (pairs[i][0] + pairs[i + 1][0])
        return float(pairs[i][0])
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
        slope of the beach face). For a line whose water is LOWER than the
        bottom knot's -- a new lowest tide beyond the others -- the
        elevation is taken to fall at no less than 1/`flatten` of the fit's
        average rate, so only a line whose water level is too HIGH for how
        far out it lies counts: a line far out on the water, which no
        flattening explains.
        A line whose water is HIGHER than the bottom knot's and yet lies
        seaward of it is out of order already, by the height difference,
        and gets no margin for a flatter terrace: the rest of the way out
        is counted at the fit's average rate. Further out each row covers
        more ground, which makes up for the flatter terrace about as much as
        it costs. Review, Oct 2026: a c2 stray 32-44 m out on the water,
        0.42 m above a spring-low line of an hour later and 52 px beyond it,
        was +0.62 m at half the rate -- under the threshold, kept whole --
        and is +0.82 m at the average rate; the bed it was drawn on, 40 m out
        on the 1:67 terrace beyond -1.2 m, is 0.87 m below its water level.
        The local slope of the lowest knots was tried and came out flatter
        still (0.0041-0.0045 m/px against the true 0.0066): the lowest lines
        are the noisiest, 0.1 m of water-level error moving one 7 m on the
        flat terrace.
    """
    sgn, ky, kz = fit
    yy = sgn * np.asarray(y, float)
    z = np.asarray(z, float)
    r = z - np.interp(yy, ky, kz)
    span = ky[-1] - ky[0]
    slope = (kz[-1] - kz[0]) / span if span > 0 else 0.0
    seaward = yy < ky[0]
    rate = np.where(z > kz[0], slope, slope / flatten)
    z_max = kz[0] - (ky[0] - yy) * rate
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
    Judge the bin's most seaward lines one at a time against the OTHER kept
    lines, from the most seaward inwards. A line fails when its residual
    exceeds the threshold and it lies more than --noise-k x --row-noise-px
    (4 x 2 px) seaward of where those lines put its water level. The
    peeling stops at the first line that passes -- unless it passed by
    leaning on a TWIN: the most seaward of the other lines is within 0.1 m
    of its level (TWIN_M) and fewer than 3 of them lie within 8 px of it.
    Then it is held back, out of the reference, while the twin is judged in
    turn; if the twin fails, the lines held back are judged again without
    it (up to SEAWARD_RUN = 3 lines in a run). Returns (residuals with
    those judgements in, keep).

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
    be within centimetres of the lowest water level.

    WHY HOLD BACK: strays come in runs. Glare or fog lasts several frames at
    slack low water (the station's c2 frames of 2 Oct 2026 15:00-16:30 all
    had strays), and the same light comes back at the same hour the next
    days. Stopping at the first line that passed, the stray judged first had
    its twin in the reference -- a few pixels landward at almost the same
    level -- so it passed, and the twin, the line it leaned on, was never
    judged (review, Oct 2026: a -0.60 m stray 25-28 m out alone, 304 of 338
    points dropped; with the next frame's -0.59 m stray, both kept whole;
    the real 2 Oct 15:00 stray relabelled to the window's lowest level
    dropped whole, with a copy 5 min later both kept whole). A run of strays
    from one glare or fog lies at one level -- slack water moves less than
    0.1 m in an hour -- each a few pixels from the next, and as a run they
    are consistent with the order; only the lines landward of the run can
    say that its level does not belong out there. So a line that passes by
    leaning on a line of its own level is held back while that line is
    judged. An honest lowest-tide line either lies among the lines of its
    own slack water (3 or more within 8 px: the peeling ends at once, as
    before) or leans on a line of a different level, and is judged as
    before. Not held back: a line leaning on a line more than 0.1 m from
    its level. Peeling past those too caught more of a scattered week of
    strays (review's 8-stray week: 31% -> 42% of their points) but dropped
    honest spring-low lines where the profile is far flatter per pixel at
    the low end than above it (c2's left edge: -0.71 m and -1.13 m lines
    35-57 px beyond a -0.39 m line, each judged by extrapolation from the
    upper beach; megacusps 3,178 -> 3,310 honest points). So strays on
    different days at levels 0.1 m or more apart, each a few metres beyond
    the next and within the threshold of it, still vouch for one another:
    in that 8-stray week 31% of their points go, as before (the DEM's
    3-frames-per-cell rule kept them out of its cells).
    A stray ABOVE the window's lowest water level but further out than its
    lowest line (review, Oct 2026: 0.3-0.5 m above a spring-low line) is
    out of order with that line; residuals() counts the height difference
    and the rest of the way out at the fit's average rate, and it goes. The
    row guard keeps a few pixels of detector scatter from counting as
    metres where the waterline barely moves with the tide (c2's far field).
    """
    keep = keep.copy()
    guard = p.noise_k * p.row_noise_px
    held = []                         # passed by leaning on a twin, out of the reference
    for i in np.argsort(ROW_SIGN * y, kind="stable"):
        if not keep[i]:
            continue
        rest = keep.copy()
        rest[i] = False
        rest[held] = False
        if rest.sum() < 2 or n_eff(w[rest]) < p.min_frames \
                or z[rest].max() - z[rest].min() < p.min_z_range:
            break
        fit = fit_column(z[rest], y[rest], w[rest])
        if fit is None:
            break
        judge = np.array([i] + held, np.int64)
        ri = residuals(fit, z[judge], y[judge], p.extrapolate_flatten)
        gap = np.interp(z[judge], fit[2], fit[1]) - fit[0] * y[judge]
        out = (ri > thr) & (gap > guard)
        if out[0]:
            # out of order with the lines landward of it -- and so are the
            # lines held back beyond it that the same lines contradict
            r[judge[out]] = np.maximum(r[judge[out]], ri[out])
            keep[judge[out]] = False
            held = [k for k in held if keep[k]]
            continue
        among = int((ROW_SIGN * y[rest] <= ROW_SIGN * y[i] + guard).sum())
        lean = np.flatnonzero(rest)[np.argmin(ROW_SIGN * y[rest])]   # the line it leans on
        if among >= SEAWARD_RUN or abs(z[lean] - z[i]) > TWIN_M or len(held) >= SEAWARD_RUN - 1:
            break                     # passed among the lines landward of it, or on its own merits
        held.append(i)
    return r, keep


def time_test(t, z, y, ref, cand, sign, p):
    """
    Second opinion, in time, on lines the window's fit put out of order:
    `cand` (indices) failed it, SEAWARD (sign +1, residual > threshold) or
    LANDWARD (sign -1). Each is judged again against the reference lines
    (`ref`, a mask: those that passed the fit), split into those BEFORE it
    and those AFTER it in time -- each side its own fit, weighted by
    closeness from its nearest line, with its own threshold
    (side_reference()). A side CONFIRMS the failure (residual beyond its
    threshold, the same way), AGREES with the line (within half its
    threshold), or cannot say (in between; too few lines; or the line lies
    beyond the end of its lines in the direction it failed, where only a
    bound or an extrapolation speaks -- see verdict()). The line is EXCUSED
    -- kept, as beach change -- when a side that can say agrees with it.
    Where no side agrees and one cannot say, the lines of the line's own
    tide (+/- 12.42 h) judge it instead -- never against a side that agrees:
    in the tide of a cut they hold the hours before it too (review, Oct
    2026: a 0.8 m cut on the record's first day, 26 Sep 13:00-19:00; the
    19:30 and 20:00 c1 lines had too few lines before them to judge, the
    lines after them agreed, and their own tide dropped them whole -- for
    good, as nothing ever comes before them). Otherwise the window's verdict
    stands -- except for a LANDWARD failure that nothing after it can judge
    (the end of the record, a column the later lines do not reach): neither
    the past nor its own tide, which holds the hours before a cut in that
    tide, drops it. A storm cuts a beach in hours, and the lines after the
    cut are landward of everything before it; the beach builds back seaward over
    days to weeks, and the errors seen on this station lie seaward (c1's
    envelope floor, c2's lines on the water), so a seaward failure the past
    alone confirms still goes. A line kept this way is judged again in the
    next run, when lines after it exist. (Review's cut in the last hours of
    the record, 0.8 m on 5 Oct 11:00-14:00: the version before dropped 7
    honest c1 lines whole that evening and 4 still the next day; none now.)
    Returns a mask over `cand`: True = confirmed (drop), False = excused.

    WHY: a beach that changes puts its lines out of order with the lines of
    the beach before. After the upper beach is cut, the next days' mid-tide
    lines lie where the old beach had its high-tide lines. Where the window
    holds no post-cut lines higher than them -- c1's high tides fall below
    the floor of its search envelope after a cut; neaps and night-time
    high tides never reach them -- every knot above them comes from before
    the cut, the fit inverts their rows to the OLD beach's elevations, and
    they fail LANDWARD by more than the cut (review, Oct 2026: 0.8 m cut on
    2 Oct, residuals -0.8 to -1.5 m; c1 lost 3 of 16 lines whole that day
    and 8 of 16 the next, 23% of its points; a 0.6 m cut in a 60-day
    record cost 6 of 17 lines that day, still dropped with 35 days of later
    record). The 7-day DEM then kept the pre-storm bed in 44 cells south of
    the seam, 0.4-0.8 m above the real one, and dem_change.py understated
    the erosion. Time weighting cannot outvote lines that do not exist, and
    no threshold that still catches the errors sits above such a cut. But
    the change shows in WHICH lines disagree: those from before the cut
    contradict a post-cut line, those after it agree with it (or all lie
    below it, so it is past the top of their fit). A line on the wrong
    feature contradicts both. The same holds the other way after
    accretion, for lines that look seaward of the old beach. So a cut is a
    change, not an error: a line goes only if the lines on both sides of it
    in time contradict it, or one does and the other cannot say. The
    bottom-up pass (bottom_up()) asks the same of its lower-water
    references; this asks it of the final judgement.

    ORDER: the candidates are judged from the end of the water-level range
    their failure points away from -- landward failures from the highest
    water down, seaward ones from the lowest up -- and each line excused
    joins the reference for the next. After a cut, the highest post-cut
    line is past the top of the lines after it; once excused it is the
    post-cut beach the lower ones are measured against. A line confirmed
    out of order never vouches for another: the references are otherwise
    only lines that passed the fit, since c1's high-tide errors come many
    to a tide and next to each other in time they would vouch for each
    other.

    AGREEMENT, not just 'under the threshold': with two references a line
    gets two chances, and in the station's untidy lines a side 0.4-0.7 m
    out of order passes while the other says 0.9-1.2 m. On the 29 Sep -
    5 Oct 2026 contours, excusing on 'under the threshold' let c1 keep
    lines 0.8 m seaward of all the others at their row (30 Sep 16:00 and
    3 Oct 11:00, +0.53 m) and dropped 12% fewer high-tide points; within
    half the threshold (0.375 m: AGREE) it drops as many as the version
    without this test (28,951 points above +0.5 m, against 28,970), while
    honest post-cut lines agree with the lines after them within 0.1-0.3 m.
    """
    order = np.argsort(t, kind="stable")
    ref = ref.copy()
    state = {}
    guard = p.noise_k * p.row_noise_px       # within this of an end knot: at it

    def reference(key):
        if key not in state["cache"]:
            sel = state["ref_i"][key[0]:key[1]]
            state["cache"][key] = side_reference(t, z, y, sel, key[2], p)
        return state["cache"][key]

    def refresh():
        state["ref_i"] = order[ref[order]]          # reference lines in time order
        state["tr"] = t[state["ref_i"]]
        state["cache"] = {}

    def verdict(side, i, span):
        """True: out of order against this reference by more than its
        threshold (confirms); False: within AGREE x its threshold (agrees:
        the line belongs to the beach this reference describes); None: it
        cannot say -- in between, or the line lies beyond the end of the
        reference's lines (by more than --noise-k x --row-noise-px) in the
        direction it failed, where the residual is only a bound or an
        extrapolation and can confirm but not agree. The exception is a
        line LANDWARD of all of the reference's lines whose water is at
        least as high as their highest: in order, however flat the upper
        beach -- the line after a cut. That holds only where the reference
        lacks higher water because none was seen in that stretch of time
        (`span`: c1's envelope floor, neaps, night). If lines at higher
        water were seen and set aside, the reference's top is where the
        rejections cut it and says nothing (the station's c1: on the right
        of the photo the high-tide lines lie seaward and are set aside, so
        the lines left there end at mid-tide)."""
        if side is None:
            return None
        fit, t_m = side
        ri = sign * float(residuals(fit, z[i:i + 1], y[i:i + 1], p.extrapolate_flatten)[0])
        if ri > t_m:
            return True
        if ri > AGREE * t_m:
            return None             # less out of order, but no agreement either
        yy = fit[0] * y[i]
        if sign > 0 and yy < fit[1][0] - guard:
            return None             # seaward of every line it holds
        if sign < 0 and yy > fit[1][-1] + guard:
            if z[i] < fit[2][-1]:
                return None         # landward of every line, and lower than the highest
            if (span & ~ref & (z > z[i])).any():
                return None         # its top is where the rejections cut it
        return False

    refresh()
    out = np.ones(cand.size, bool)
    # From the end of the water-level range the failures point away from
    # (landward failures from the highest water down, seaward ones from the
    # lowest up), and each line excused joins the reference for the next:
    # after a cut, the highest post-cut line is past the top of the lines
    # after it, and once excused it is the post-cut beach the lower ones
    # are measured against. A line confirmed out of order never vouches.
    for k in np.argsort(sign * z[cand], kind="stable"):
        i = cand[k]
        tr = state["tr"]
        a = int(np.searchsorted(tr, t[i], "left"))
        b = int(np.searchsorted(tr, t[i], "right"))
        votes = [verdict(reference((0, a, True)), i, t < t[i]),
                 verdict(reference((b, tr.size, False)), i, t > t[i])]
        if sign < 0 and votes[1] is None:
            # LANDWARD of the beach before it, and nothing after it can
            # say: the past alone does not drop it -- nor its own tide,
            # which holds the hours before a cut in that tide
            out[k] = False
        elif False in votes:
            # a side that can judge it agrees: it belongs to the beach of
            # that side -- change. Its own tide does not overrule that: in
            # the tide of a cut it holds the beach before the cut too
            out[k] = False
        else:
            if None in votes:
                tide = verdict(reference((int(np.searchsorted(tr, t[i] - TIDE_CYCLE, "left")),
                                          int(np.searchsorted(tr, t[i] + TIDE_CYCLE, "right")), None)),
                               i, np.abs(t - t[i]) <= TIDE_CYCLE)
                if tide is not None:
                    votes = [tide]
            # out of order only if every reference that can judge it says
            # so; if none can, the window's verdict stands
            out[k] = all(v for v in votes if v is not None)
        if not out[k]:
            ref[i] = True
            refresh()
    return out


def judge_bin(t, z, y, w, p, bottom=None):
    """
    Judge one column bin's lines (times t, rows y, water elevations z,
    weights w).

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
    5. Change, not error: every line that failed step 3 or 3b is judged
       again against the lines before it and after it in time
       (time_test()), seaward failures first; a line one side agrees with
       is kept as beach change ('excused', with the direction it failed),
       and joins the reference for the landward ones.

    Returns (fit, residuals, threshold, reversed, rank_corr, excused), or
    None when the bin cannot be judged (everything in it is kept).
    `threshold` is a scalar or, for a REVERSED bin judged by step 2 alone,
    per line; `excused` per line (+1 / -1 / 0).
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
                # 5. out of order with the lines on one side in time only: change
                excused = np.zeros(z.size, np.int8)
                ref = keep.copy()
                for s in (1, -1):
                    cand = np.flatnonzero(~keep & (s * r > thr))
                    if cand.size:
                        ok = cand[~time_test(t, z, y, ref, cand, s, p)]
                        excused[ok] = s
                        ref[ok] = True          # kept as change: part of the beach now
                out = (fit, r, thr, rev, rc, excused)
    if out is None and rev and bu is not None and np.isfinite(bu[1]).any():
        keep_bu, r_bu, thr_bu = bu
        judged = np.isfinite(r_bu)
        # lines the bottom-up pass could not reach (the lowest layers) pass,
        # and so do lines it kept (on one side in time, or by the row guard)
        r_out = np.where(judged, r_bu, 0.0)
        t_out = np.where(judged & ~keep_bu, thr_bu, np.inf)
        fit = fit_column(z[keep_bu], y[keep_bu], w[keep_bu]) if keep_bu.sum() >= 2 else None
        out = (fit, r_out, t_out, rev, rc, np.zeros(z.size, np.int8))
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
    window fits {(day, bin): (fit, frames, weights, threshold)}, the
    status of every (day, bin) the day's lines reach:
    {(day, bin): (status, rank_corr)}, status one of 'judged',
    'reversed' (judged, rows ran against the water level), 'few' (too few
    lines or too little elevation range: kept), 'unclear' (no clear
    relation even after setting the out-of-order lines aside: kept), and
    per frame x bin the direction a line was out of order with the fit but
    NOT with the lines on one side of it in time (+1 seaward, -1 landward,
    0 neither; time_test()): kept as beach change.
    """
    R, n_bins = bin_rows(frame_local, col, row, n_frames, p.column_bin)
    resid = np.full((n_frames, n_bins), np.nan)
    thresh = np.full((n_frames, n_bins), np.nan)
    excused = np.zeros((n_frames, n_bins), np.int8)
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
            judged_bin = judge_bin(epochs[idx], z, y, w, p,
                                   bottom=lambda: tuple(a[idx] for a in bottom_for(b)))
            if judged_bin is None:
                status[(d, b)] = ("unclear", rank_corr(z, ROW_SIGN * y))
                continue
            fit, r_all, thr, rev, rc, exc = judged_bin
            status[(d, b)] = ("reversed" if rev else "judged", rc)
            mine = np.isin(idx, judged)
            resid[idx[mine], b] = r_all[mine]
            thresh[idx[mine], b] = thr[mine] if np.ndim(thr) else thr
            excused[idx[mine], b] = exc[mine]
            if fit is not None:
                if np.ndim(thr):
                    fin = thr[np.isfinite(thr)]
                    t_pts = float(np.median(fin)) if fin.size else p.max_residual
                else:
                    t_pts = thr
                fits[(d, b)] = (fit, idx, w, t_pts)
    return resid, thresh, R, fits, status, excused


def decide(resid, thresh, p, excused=None):
    """Failing bins per frame, and the frames to drop whole. A bin the
    window's fit put out of order but the lines on one side of it in time
    did not (`excused`, time_test()) is not failing: change, not error."""
    judged = np.isfinite(resid)
    with np.errstate(invalid="ignore"):            # NaN = not judged (numpy < 1.18 warns)
        fail = judged & (np.abs(resid) > thresh)
    if excused is not None:
        fail &= excused == 0
    n_judged = judged.sum(axis=1)
    n_fail = fail.sum(axis=1)
    whole = (n_judged >= p.min_judged_bins) & (n_fail > p.frame_fraction * n_judged)
    return fail, whole, n_judged, n_fail


def point_drops(fl, col, row, fail, whole, days, elevs, fits, p, excused=None):
    """
    Which points to drop. Whole frames: all of them. Otherwise each point
    in a failing bin, or in a bin next to one, is judged on its OWN row
    against that bin's fit and threshold. A wrong segment starts and ends
    part way through a bin, where the bin's single row is a blend of the
    right and the wrong part; judging the points there one by one drops the
    wrong part up to where it leaves the right line, and keeps the rest.
    Points away from any failing bin are never judged one by one, so pixel
    noise on a good line cannot speckle it with drops. Nor is a point
    dropped for lying out of order the way the time test excused its line
    in that bin (`excused`): a line kept as beach change in one bin is not
    cut back to the bin next door that failed.
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
        out = np.abs(r) > thr
        if excused is not None:
            e = excused[fl[sel], int(b)]
            out &= ~(((r > 0) & (e > 0)) | ((r < 0) & (e < 0)))
        drop[sel[out]] = True
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


def bin_summary(cam, status, p, change=None):
    """
    Counts of column bins x days by status, for the CONSISTENCY line, and
    WARNING lines for REVERSED bins: rows running against the water level
    there are for a person to look at even though the filter drops what it
    can. Physics says only that the order is broken, not which lines are
    wrong: the high-tide lines may sit seaward (the detector unable to reach
    the real line, c1's search-envelope floor), the low-tide lines landward
    (the wet/dry line on the ebb), or the beach changed within the window.
    `change`: {(day, bin): (lines out of order with the window's fit, of
    those kept as beach change by time_test())}; where in most reversed
    bin-days nothing is out of order beyond the threshold, or most of it is
    kept as change, the first line (the one the cron logs) says 'beach
    change?' (review, Oct 2026: it blamed the high-tide lines and the
    envelope on an honest week after the beach was lowered 0.8 m).
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
        def explained(n_out, n_exc):
            # nothing beyond the threshold (reversed by rank alone), or most
            # of what was beyond it kept as beach change by the test in time
            return n_out == 0 or n_exc > 0.5 * n_out
        mostly = sum(1 for d, b, _ in rev if change and explained(*change.get((d, b), (0, 0))))
        warn.append(f"WARNING: camera {cam}: waterline rows run AGAINST the water level in pixel "
                    f"columns {ranges}"
                    + (" (beach change? in most of them no line is out of order beyond the "
                       "threshold, or it is kept as change)" if mostly > 0.5 * len(rev) else ""))
        warn.append(f"         ({len(rev)} column bin-days on {len(days_rev)} of {len(days_all)} days, "
                    f"{first} to {last}; rank correlation {rcs.min():+.2f} to {rcs.max():+.2f}).")
        if recent:
            warn.append("         In the last {} days: columns {}.".format(
                p.plot_days, ", ".join(f"{a}-{b}" for a, b in _col_ranges(recent, p.column_bin))))
        warn.append("         Higher water must put the line LOWER in the photo (nearer the camera);")
        warn.append("         here it does not. Either the high-tide lines sit seaward (the detector")
        warn.append("         cannot reach the real line: search envelope floor above the high-tide")
        warn.append("         waterline?), or the low-tide lines sit landward (the wet/dry line?), or")
        warn.append(f"         the beach changed in the window ({mostly} of these bin-days: no line out")
        warn.append("         of order beyond the threshold, or most of those kept as change). Lines")
        warn.append("         seaward of lower-water lines there are dropped -- the low-tide lines are")
        warn.append("         taken as right; look at the diagnostic plot.")
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
        resid, thresh, R, fits, status, excused = judge_camera(epochs, elevs, fl, col, row,
                                                               frames.size, args)
        fail, whole, n_judged, n_fail = decide(resid, thresh, args, excused)
        fday = np.floor(epochs / DAY)
        change = {}
        for (d, b), kind in status.items():
            if kind[0] == "reversed":
                on_day = fday == d
                with np.errstate(invalid="ignore"):
                    n_out = int((on_day & (np.abs(resid[:, b]) > thresh[:, b])).sum())
                change[(d, b)] = (n_out, int((on_day & (excused[:, b] != 0)).sum()))
        bins_msg, warn = bin_summary(cam, status, args, change)
        for line in warn:
            print(line)

        drop = point_drops(fl, col, row, fail, whole, np.floor(epochs / DAY), elevs, fits, args,
                           excused)
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
            # failing bins can disagree in sign (seaward in some columns,
            # landward in others): the median of the majority's, and the
            # count of each, rather than a median across both that calls a
            # line seaward when most of it is landward
            n_sea, n_land = int((rr > 0).sum()), int((rr < 0).sum())
            maj = rr[rr > 0] if n_sea >= n_land else rr[rr < 0]
            med = float(np.median(maj)) if maj.size else float("nan")
            worst = float(rr[np.argmax(np.abs(rr))]) if rr.size else float("nan")
            thr = float(np.nanmedian(thresh[f][fail[f]])) if fail[f].any() else float("nan")
            side = "lower" if med > 0 else "higher"
            where = "seaward" if med > 0 else "landward"
            mixed = (f" ({n_sea} bin(s) seaward, {n_land} landward)"
                     if fail[f].any() and n_sea and n_land else "")
            if whole[f]:
                action = "frame"
                reason = (f"{n_fail[f]} of {n_judged[f]} judged column bins out of order "
                          f"(> {args.frame_fraction:.0%}){mixed}: whole line dropped; median "
                          f"residual {med:+.2f} m, i.e. it sits {where} of where its water level "
                          f"belongs")
            else:
                action = "segments"
                reason = (f"out of order in {n_fail[f]} of {n_judged[f]} judged column bins{mixed}: "
                          f"the line sits where the other lines put water {abs(med):.2f} m {side} "
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
        n_change = int(((excused != 0).any(axis=1) & ~whole).sum())
        msg = (f"{cam}: {n_whole} line(s) dropped, {affected.size - n_whole} trimmed, "
               f"{n_drop:,} of {n_cam_pts:,} points ({100.0 * n_drop / max(n_cam_pts, 1):.1f}%); "
               f"{n_change} line(s) out of order with one side in time only, kept as beach change; "
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
                              resid=resid, thresh=thresh, fail=fail, R=R, fits=fits, excused=excused,
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
    # the same words as the elevation maps' colour bar (daily_elevation_map.py)
    elev_words = ("water level + wave setup" if data["elev_col"] == "beach_elevation_navd88"
                  else "water level")
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
    # Numbers in staggered rows above the lines, spread sideways with a
    # leader to each dropped part, so rejections in one place stay readable.
    # As many rows as the numbers need to sit side by side (at least two),
    # stacked upwards from 5% of the photo above the highest dropped line --
    # and, where that would leave the photo, downwards from its top edge
    # instead: both rows were once clamped to the same height there, and on
    # the station's c2 week 24 of 29 numbers overlapped another.
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
        n_tiers = max(2, int(np.ceil(len(order) * sep / (0.95 * width))))
        step = 0.035 * height
        ly_first = top - 0.05 * height                      # the row nearest the lines
        if ly_first - (n_tiers - 1) * step < 0.03 * height:
            ly_first = 0.03 * height + (n_tiers - 1) * step
        for tier in range(n_tiers):
            ids = order[tier::n_tiers]
            xs = [anchors[n][0] for n in ids]
            for k in range(1, len(xs)):                      # spread rightwards
                xs[k] = max(xs[k], xs[k - 1] + sep)
            if xs and xs[-1] > width - sep / 2:              # then pull back inside
                xs[-1] = width - sep / 2
                for k in range(len(xs) - 2, -1, -1):
                    xs[k] = min(xs[k], xs[k + 1] - sep)
            ly = ly_first - tier * step
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
                    f"day(s) -- high-tide lines seaward, low-tide lines landward, or beach change")
    ax.set_xlim(0, width); ax.set_ylim(height, 0); ax.axis("off")
    status = "" if info["applied"] else "   [NOT APPLIED: would drop too much, see log]"
    ax.set_title(f"{cam.upper()}  {first} to {last}: {len(listed)} line(s) with parts dropped as out of "
                 f"order (thick; dashed = whole line){status}\nthin white = kept lines, {bg_note}"
                 f"{rev_note}", fontsize=10)
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cax = fig.add_axes([0.03, 0.035, 0.25, 0.015])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label(f"{elev_words} of the line (m NAVD88)", fontsize=8)
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
    ax2.set_xlabel(f"{elev_words} (m NAVD88)", fontsize=8)
    ax2.set_ylabel("detected row (px), image orientation", fontsize=8)
    ax2.tick_params(labelsize=7)
    ax2.grid(alpha=0.3)

    # The table: what its columns mean first, so it is never the part cut
    # off; then as many rows as the space below the plot holds, and a count
    # of the rest only when there is a rest. (It printed 34 rows and the
    # note after them whatever the space: on a busy week the note fell off
    # the bottom, and '... and 0 more' was printed.)
    ax3 = fig.add_axes([0.66, 0.02, 0.33, 0.47]); ax3.axis("off")
    fs = 7.4
    room = int(0.47 * fig.get_figheight() * 72 / (fs * 1.2))        # text lines that fit
    head = [f"resid = {elev_words} minus the elevation the other lines",
            "give its row (m): > 0 the line lies seaward of where its level",
            "belongs, < 0 landward.",
            "",
            "#   capture (UTC)     elev    pts   resid  action"]
    rows = []
    for n, f in enumerate(listed, 1):
        r = rep.get(data["names"][frames[f]])
        if r is None:
            continue
        rows.append(f"{n:<3d} {r['capture_time_utc'][5:16].replace('T', ' ')}  "
                    f"{r['elevation_navd88']:+5.2f}  {r['points_dropped']:5d}  "
                    f"{r['median_residual_m']:+5.2f}  {r['action']}")
    fit_rows = room - len(head)
    if len(rows) > fit_rows:
        rows = rows[:fit_rows - 1] + [f"... and {len(rows) - (fit_rows - 1)} more: "
                                      "see the report CSV"]
    ax3.text(0, 1, "\n".join(head + rows), va="top", ha="left", family="monospace", fontsize=fs)
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
      6. STRAY ABOVE THE LOWEST LINE (review, Oct 2026): c2 on the
         RTK-style beach (below), the storm week's water levels: a line at
         spring tides 0.3-0.5 m above the lowest line of its window leaves
         the waterline on the right of the photo for a line on the water
         10 m beyond where that lowest level meets the beach (~30 m beyond
         its own). Half the fit's rate beyond the lowest line put it under
         the threshold (88% of that stretch dropped, the ramp into it kept,
         in the version before); at least 95% of the points of its stretch
         on the water must be dropped, and less than 0.5% of the honest
         points.
      7. EROSION ON THE SURVEYED BEACH (review, Oct 2026): c1 alone, the
         RTK-style profile (foreshore 0.13, terrace 0.035 to -1.2 m, 0.015
         beyond, the survey's alongshore variation), the storm week's water
         levels and setup, lowered by 0.8 m and by 0.5 m over 6 h from
         2 Oct 12:00 with the generator's smooth taper; each line loses the
         columns below the floor of c1's search envelope, as the detector
         would. After the cut the high tides fall below the floor, so the
         post-cut mid-tide lines have nothing above them but the old beach.
         No honest line may be dropped whole, and less than 1% of the
         points after the cut (the version before: 9 lines whole at 0.8 m,
         3 at 0.5 m).
      8. CUT ON THE RECORD'S FIRST DAY (review, Oct 2026): case 7's 0.8 m cut
         on 26 Sep 13:00-19:00 instead, so the post-cut lines of that evening
         have too few lines before them to judge and their own tide (+/-
         12.42 h) still holds the beach before the cut. The lines after them
         agree with them; the same limits as case 7. (When the own tide
         overruled an agreeing side, the 19:30 and 20:00 lines were dropped
         whole -- for good, as nothing ever comes before them; the survey-date
         products filter 7-day windows, so every window has such a first day.)
      9. A RUN OF STRAYS (review, Oct 2026): case 4's c2 lines, two
         consecutive frames at slack low water leaving the waterline on the
         right of the photo for lines 40 m and 36 m out on the water, both at
         the window's lowest level. Each must lose at least 90% of its points
         out there (both were kept whole when the first stray, judged against
         the others, leaned on its twin), and the honest lines less than 0.5%.
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
        alongshore undulation amplitude and wavelength (default 5, 150).
        beach['cut'] = (t, seconds, m): as 'storm', with the review's
        generator's smooth taper (smoothstep from -0.2 m to +1.0 m).
        beach['rtk']: the profile the 29 Sep 2026 RTK checkshots describe
        instead -- foreshore 0.13 from the berm crest (+3.75 m, s = 12) to
        0 m, low-tide terrace 0.035 to -1.2 m, 0.015 beyond -- with the
        survey's alongshore variation (the +1.2 m contour ~11 m further
        seaward around 130 m south of the cameras, 3 m undulation every
        140 m)."""
        if "storm" in beach:
            tc, dur, mag = beach["storm"]
            f = float(np.clip((t - tc) / dur, 0, 1))
            m = mag * f * f * (3 - 2 * f)
            if h > -0.2 and m > 0:
                h = h + m if h > 1.0 - m else (h + 0.2 * m / 1.2) / (1 - m / 1.2)
        if "cut" in beach:
            tc, dur, mag = beach["cut"]
            f = float(np.clip((t - tc) / dur, 0, 1))
            m = mag * f * f * (3 - 2 * f)
            if m > 0:
                h0 = np.linspace(-0.2, 4.0, 2101)
                x = np.clip((h0 + 0.2) / 1.2, 0, 1)
                h = float(np.interp(h, h0 - m * x * x * (3 - 2 * x), h0)) if h > -0.2 else h
        if beach.get("rtk"):
            s_f = 12 + 3.75 / 0.13
            if h >= 0:
                s = 12 + (3.75 - h) / 0.13
            elif h >= -1.2:
                s = s_f - h / 0.035
            else:
                s = s_f + 1.2 / 0.035 + (-1.2 - h) / 0.015
            return s + 11.0 * np.exp(-0.5 * ((a + 130.0) / 55.0) ** 2) + 3.0 * np.sin(2 * np.pi * a / 140.0)
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
    base4 = [list(f) for f in frames]           # case 9 starts from the same lines
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

    # 9. A RUN OF STRAYS (review, Oct 2026): case 4's c2 lines, and in two
    #    consecutive frames at slack low water -- glare or fog lasts -- the
    #    line leaves the waterline on the right of the photo for one 40 m and
    #    36 m out on the water, both at the lowest water level of the window.
    #    Each was the other's reference: judged against all the others, the
    #    first leaned on its twin a few pixels landward, passed, and the
    #    peeling stopped before the twin was judged; both were kept whole.
    if base4:
        c2 = sorted((f for f in base4 if f[1] == "c2"), key=lambda f: f[2])
        pairs = [(a, b) for a, b in zip(c2[:-1], c2[1:])
                 if b[2] - a[2] <= 3600 and 6 <= (a[2] - t_start) / DAY < 8]
        if not pairs:
            failures.append("run of strays: no two consecutive c2 frames to make it from")
            pairs = [(None, None)]
    if base4 and pairs[0][0] is not None:
        twin = min(pairs, key=lambda ab: ab[0][3] + ab[1][3])
        low = min(f[3] for f in c2 if abs(f[2] - twin[0][2]) <= 3.5 * DAY
                  and f is not twin[0] and f is not twin[1])
        for k, (f, metres) in enumerate(zip(twin, (40.0, 36.0))):
            def on_water(a, metres=metres):
                return metres * np.clip((260.0 - a) / 60.0, 0, 1)
            dev = line("c2", f[3], f[2], beach) - line("c2", f[3], f[2], beach, out=on_water)
            dc = np.interp(f[4], cols, np.nan_to_num(dev, nan=0.0))
            f[5] = f[5] - dc
            f[6] = dc >= 0.95 * np.nanmax(dev)
            f[3] = low - 0.01 + 0.005 * k
        frames = [tuple(f) for f in c2]
        dropped, out, summary, plot_info, args = run_case("twin", frames)
        n_h = d_h = 0
        caught = []
        for fname, cam, ep, z, c, r, far, h in frames:
            dr = dropped[fname]
            if far is None:
                n_h += dr.size; d_h += int(dr.sum())
            else:
                caught.append((int((dr & far).sum()), int(far.sum()), z))
        print("  two strays in consecutive frames at the window's lowest level (c2): "
              + "; ".join(f"{a:,} of {b:,} points {m:g} m out dropped ({z:+.2f} m)"
                          for (a, b, z), m in zip(caught, (40, 36)))
              + f"; honest lines {d_h:,} of {n_h:,} ({d_h / max(n_h, 1):.2%})")
        for line_ in summary:
            print(f"  CONSISTENCY {line_}")
        for (a, b, z), m in zip(caught, (40, 36)):
            if a < 0.9 * b:
                failures.append(f"run of strays: only {a} of {b} points of the {m} m one dropped (need 90%)")
        if d_h >= 0.005 * n_h:
            failures.append(f"run of strays: {d_h / max(n_h, 1):.2%} of honest points dropped (limit 0.5%)")

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

    # 6. A c2 stray at spring tides, 0.3-0.5 m ABOVE the lowest line of its
    #    window, on the surveyed beach (review, Oct 2026: env05, 26 Sep
    #    21:00, -0.71 m, 0.42 m above a -1.13 m line of an hour later, 32-44 m
    #    out on the water on the right of c2's photo but only ~50 px beyond
    #    that line, kept whole). The storm week's water levels and setup, c2
    #    keeping 60% of its frames.
    if sig_t:
        rng = np.random.default_rng(17)
        beach = {"rtk": True}
        frames = []
        for d in range(8):
            for k in range(24):
                ep = t_start + d * DAY + 11 * 3600 + k * 1800
                if rng.random() > 0.6:
                    continue
                h = float(np.interp(ep, sig_t, sig_wl)) + 0.09
                hs, tp = float(np.interp(ep, sig_t, sig_hs)), float(np.interp(ep, sig_t, sig_tp))
                z = h + 0.037 * np.sqrt(hs * 9.81 * tp * tp / (2 * np.pi))
                noise = 1.0 + 0.5 * max(0.0, hs - 1.5)
                r = line("c2", z + rng.normal(0, 0.05 + 0.03 * hs), ep, beach) \
                    + noise * (rng.normal(0, 1.2, cols.size) + rng.normal(0, 1.5))
                c, rr, _ = keep_rows("c2", r, cols)
                if c.size >= 20:
                    frames.append([f"{int(ep)}.c2.rtk", "c2", ep, z, c, rr, None, h])

        def window_low(f):
            return min(g[3] for g in frames if abs(g[2] - f[2]) <= 3.5 * DAY and g is not f)
        stray = min((f for f in frames if (f[2] - t_start) / DAY < 3),
                    key=lambda f: abs(f[3] - window_low(f) - 0.4))
        above = stray[3] - window_low(stray)

        def on_water(a):
            # 10 m beyond where the window's lowest water level meets the
            # beach (~30 m beyond its own level, ~60 px beyond the lowest
            # line in the photo, as in the review's case), ramping in over
            # 60 m from 200 m alongshore of the cameras towards the seam
            beyond = contour_s(window_low(stray), a, stray[2], beach) - contour_s(stray[3], a, stray[2], beach)
            return (beyond + 10.0) * np.clip((260.0 - a) / 60.0, 0, 1)
        dev = line("c2", stray[3], stray[2], beach) - line("c2", stray[3], stray[2], beach, out=on_water)
        dc = np.interp(stray[4], cols, np.nan_to_num(dev, nan=0.0))
        stray[5] = stray[5] - dc
        far, ramp = dc >= 0.95 * np.nanmax(dev), dc >= 0.25 * np.nanmax(dev)
        frames = [tuple(f) for f in frames]
        dropped, out, summary, plot_info, args = run_case("stray_rtk", frames)
        dr = dropped[stray[0]]
        catch_s = int((dr & far).sum()) / max(int(far.sum()), 1)
        n_h = sum(v.size for k, v in dropped.items() if k != stray[0])
        d_h = sum(int(v.sum()) for k, v in dropped.items() if k != stray[0])
        print(f"  stray {above:.2f} m above its window's lowest line ({stray[3]:+.2f} m, c2, surveyed "
              f"beach), 10 m beyond it: {int((dr & far).sum()):,} of {int(far.sum()):,} points of the "
              f"stretch on the water dropped ({catch_s:.1%}); {int((dr & ramp).sum()):,} of "
              f"{int(ramp.sum()):,} of the ramp in; honest lines {d_h:,} of {n_h:,} "
              f"({d_h / max(n_h, 1):.2%})")
        for line_ in summary:
            print(f"  CONSISTENCY {line_}")
        if catch_s < 0.95:
            failures.append(f"stray above the lowest line: only {catch_s:.1%} of its points on the "
                            f"water dropped (need 95%)")
        if d_h >= 0.005 * n_h:
            failures.append(f"stray case: {d_h / max(n_h, 1):.2%} of honest points dropped (limit 0.5%)")

    # 7. Real erosion on the surveyed beach, mid-record, then neaps (review,
    #    Oct 2026): c1, with its crop and the floor of its search envelope
    #    (a line's columns below the floor are not found), the RTK-style
    #    profile lowered by 0.8 m and by 0.5 m over 6 h from 2 Oct 12:00 --
    #    the day after the storm peak, with the tide going to neaps. After
    #    the cut the high-tide lines fall below the envelope's floor, so no
    #    line after it reaches the levels the old beach had where the new
    #    mid-tide lines now lie. Every line is honest.
    # 8. The same with a 0.8 m cut on the FIRST day of the record (26 Sep
    #    13:00-19:00): nothing before the post-cut lines can vouch for them,
    #    and their own tide still holds the beach before the cut.
    try:                                  # c1's search envelope: no line found below its floor
        from waterline_detector_v5 import CAMERAS
        prof = CAMERAS["CACO05_C1"]
        ex_, hi_ = (np.array([q[0] for q in prof.envelope_points]),
                    np.array([q[2] for q in prof.envelope_points]))
        floor_c1 = crop["c1"][0] + np.interp(cols / (W - 1.0), ex_, hi_) * (crop["c1"][1] - crop["c1"][0]) + 11
    except Exception:                     # no detector here: the photo's crop alone
        floor_c1 = np.full(cols.size, float(crop["c1"][1]))
    cuts = ((0.8, datetime(2026, 10, 2, 12, tzinfo=timezone.utc)),
            (0.5, datetime(2026, 10, 2, 12, tzinfo=timezone.utc)),
            (0.8, datetime(2026, 9, 26, 13, tzinfo=timezone.utc)))
    for mag, when in (cuts if sig_t else ()):
        rng = np.random.default_rng(13)
        t_cut2 = when.timestamp()
        beach = {"rtk": True, "cut": (t_cut2, 6 * 3600.0, mag)}
        frames = []
        for d in range(10):
            for k in range(24):
                ep = t_start + d * DAY + 11 * 3600 + k * 1800
                h = float(np.interp(ep, sig_t, sig_wl)) + 0.09
                hs, tp = float(np.interp(ep, sig_t, sig_hs)), float(np.interp(ep, sig_t, sig_tp))
                z = h + 0.037 * np.sqrt(hs * 9.81 * tp * tp / (2 * np.pi))
                noise = 1.0 + 0.5 * max(0.0, hs - 1.5)
                r = line("c1", z + rng.normal(0, 0.05 + 0.03 * hs), ep, beach) \
                    + noise * (rng.normal(0, 1.2, cols.size) + rng.normal(0, 1.5))
                with np.errstate(invalid="ignore"):
                    r = np.where(r > floor_c1, np.nan, r)
                c, rr, _ = keep_rows("c1", r, cols)
                if c.size >= 20:
                    frames.append((f"{int(ep)}.c1.cut", "c1", ep, z, c, rr, None, h))
        dropped, out, summary, plot_info, args = run_case(f"cut{mag:g}_{when:%m%d}", frames)
        after = [f[0] for f in frames if f[2] >= t_cut2]
        n_a = sum(dropped[k].size for k in after); d_a = sum(int(dropped[k].sum()) for k in after)
        whole = [k for k in dropped if dropped[k].all()]
        n_all = sum(v.size for v in dropped.values()); d_all = sum(int(v.sum()) for v in dropped.values())
        fa = d_a / max(n_a, 1)
        what = (f"erosion {mag:g} m on {when.day} {when:%b %H:%M}"
                + (", the record's first day" if when.day == 26 else ""))
        print(f"  {what}, surveyed beach, c1 with its envelope: after the cut "
              f"{d_a:,} of {n_a:,} honest points dropped ({fa:.2%}), {len(whole)} line(s) whole; "
              f"all {d_all:,} of {n_all:,}")
        for line_ in summary:
            print(f"  CONSISTENCY {line_}")
        if whole:
            failures.append(f"{what}: {len(whole)} honest line(s) dropped whole (limit none)")
        if fa >= 0.01:
            failures.append(f"{what}: {fa:.2%} of the honest points after the cut dropped "
                            f"(limit 1%)")

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
                         "the bottom-up pass and the seaward end drop a line only when it is more "
                         "than --noise-k x this seaward of where the other lines put its level in "
                         "the photo (default 2).")
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

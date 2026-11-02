#!/usr/bin/env python3
"""
Tests For Partial Waterlines (Detector v6.11)
===============================================
Plain python3 (no pytest needed on the station):

    python3 test_partial_lines.py

WHAT.
  * maximal_segments: the Ruzzo & Tompa (1999) example, the recursive
    definition (best segment, then the same either side) on random scores,
    and what a run is: a stretch whose every part from either end is more
    than 80% clear, that no longer such stretch contains; ties.
  * clear_runs: long runs of signal columns kept, short ones and scattered
    columns not; PARTIAL_MAX_GAP_COLUMNS (49 bridged, 50 not) and
    PARTIAL_MIN_RUN_FRACTION (9.9% out, 10% in) at their boundaries; a
    bridged dropout still no-signal; runs cut where they meet noise, not at
    the frame's edge; the haze's chance clear columns neither extend a run
    past the haze's edge nor undo the cut there, nor lose a frame that has
    more signal; --column-limits (inside) keeps and counts runs the same;
    near_edge takes only the first run, from within so many columns of column 0.
  * flat_end: a level, edgeless end (the solver holding its row with no edge
    under it) cut back to where it levels, 12 columns more; a level line on an
    edge, or an edgeless sloping one, not; PARTIAL_FLAT_END_SLOPE and
    PARTIAL_FLAT_END_SHARPNESS at their edges; a line shorter than the window.
  * process_image on a synthetic energy field (the real smoothing, solver,
    signal test, confidence gate, temporal blend and CSV export; only the
    photo and its energy are made up):
      - a frame with signal everywhere: the CSV the same to the byte as the
        v6.10 detector's (git history; skipped without it);
      - the far half hazy (c2 after the Oct 2026 storm): discarded with partial
        lines off, saved as a partial line with them on -- only the clear run
        from the near edge flagged Has_Signal, marked Partial_Line, its rows
        on the true waterline, cut 12 columns from the haze, the hazy half
        not drawn on the overlay;
      - haze at the near end or in the middle: discarded as v6.10 did (the
        near half; the near 35%, for its confidence, never retried), or only
        the run from the near edge kept (a hazy middle); a dropout at the
        near edge does not lose the run; the review's case (a stronger edge
        beside the clear stretch, near haze) not saved;
      - whole -> partial -> whole with the temporal blend on: the whole
        frames are the blend reference, the partial one is not (its state
        untouched), and the third frame's CSV is the v6.10 detector's;
      - a flat end from column 200 cut, the line saved as partial and logged;
        one from column 150 leaving under 30%, discarded saying so; one
        shorter than the window left to the 12-column cut;
      - --column-limits: judged on the columns inside them, so a partial
        line is never saved with fewer than 30% of its columns kept;
      - PARTIAL_MIN_FRACTION at the module's default: 28% kept discarded,
        31.3% saved;
      - only scattered short stretches clear: discarded either way, saying why.
  * extract_elevation_contours.py: a CSV marked Partial_Line is judged by
    --partial-min-coverage (0.25, never above --min-coverage), every other
    frame by --min-coverage as before, however few of its columns carry
    signal; the log's partial count is of the frames used, and says so
    under --include-no-signal.
  * batch_bias_analysis.py skips a partial line (listing it) and
    sweep_profile_scale.py does not score one; the single-pair
    compare_to_ground_truth.py and bias_analysis.py warn on one: they score
    every row, and a partial line's no-signal rows are haze.
Exit status 0 when every check passes.
"""

import ast
import contextlib
import csv
import importlib.util
import io
import subprocess
import sys
import tempfile
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import waterline_detector_v5 as D          # noqa: E402
import extract_elevation_contours as X     # noqa: E402

failures = []
EPOCH = 1792886400                          # 2026-10-25 00:00 UTC
ROWS, COLS = 400, 600
V610_COMMIT = "fd7d021"                     # main with detector v6.10, before partial lines


def check(name, cond, detail=""):
    shown = not cond and not (isinstance(detail, str) and detail == "")
    print(f"  {'ok ' if cond else 'FAIL'} {name}" + (f"  ({detail})" if shown else ""))
    if not cond:
        failures.append(name)


def best_then_sides(x, lo, hi, out):
    """The recursive definition of the maximal segments (Ruzzo & Tompa 1999), by brute force."""
    best, arg = 0.0, None
    for i in range(lo, hi):
        t = 0.0
        for j in range(i, hi):
            t += x[j]
            if t > best:
                best, arg = t, (i, j)
    if arg is None:
        return
    best_then_sides(x, lo, arg[0], out)
    out.append(arg)
    best_then_sides(x, arg[1] + 1, hi, out)


def test_maximal_segments():
    print("maximal_segments")
    check("the paper's example: (4), (3), (1 2 -2 2 -2 1 5)",
          D.maximal_segments([4, -5, 3, -3, 1, 2, -2, 2, -2, 1, 5]) == [(0, 0), (2, 2), (4, 10)])
    rng = np.random.default_rng(11)
    bad = []
    for k in range(300):
        x = list(rng.normal(-0.3, 1.0, int(rng.integers(1, 40))))
        want = []
        best_then_sides(x, 0, len(x), want)
        if D.maximal_segments(x) != want:
            bad.append(k)
    check("300 random sequences: the same as the recursive definition", not bad, bad[:5])
    check("nothing positive: no segment", D.maximal_segments([-1, 0, -2]) == [])
    bad = []
    for k in range(400):
        n = int(rng.integers(1, 22))
        x = [int(v) for v in np.where(rng.random(n) < rng.uniform(0.3, 0.97), 1, -4)]
        ends_clear = [(i, j) for i in range(n) for j in range(i, n)
                      if all(sum(x[i:t + 1]) > 0 for t in range(i, j + 1))
                      and all(sum(x[t:j + 1]) > 0 for t in range(i, j + 1))]
        want = sorted(a for a in ends_clear if not any(b != a and b[0] <= a[0] and a[1] <= b[1] for b in ends_clear))
        if D.maximal_segments(x) != want:
            bad.append(x)
    check("400 random clear/no-signal sequences (+1/-4, i.e. 80%): each segment is a stretch whose every part "
          "from either end is more than 80% clear, and no longer such stretch contains it", not bad, bad[:2])
    check("a tie ends a segment (the paper's R_j >= R_k): a dropout is bridged only when what lies "
          "beyond it outweighs it", D.maximal_segments([1, 1, 1, 1, -2, 2]) == [(0, 3), (5, 5)]
          and D.maximal_segments([1, 1, 1, 1, -2, 3]) == [(0, 5)])


def haze_frame(n, edge, rng=None, p=0.0, every=None):
    """Clear columns 0..edge-1, haze beyond with chance clear columns (rate p, or every `every`th)."""
    s = np.zeros(n, bool)
    s[:edge] = True
    if every:
        s[edge::every] = True
    elif p:
        s[edge:] = rng.random(n - edge) < p
    return s


def test_clear_runs():
    print("clear_runs")
    n = 200
    m, k = D.clear_runs(np.ones(n, bool), min_run=20, max_gap=8, edge_trim=12)
    check("all signal: every column, no trim at the frame's edges", m.all() and k == 1)
    s = np.zeros(n, bool)
    s[60:140] = True
    m, k = D.clear_runs(s, min_run=20, max_gap=8, edge_trim=12)
    check("a run between noise: cut by 12 columns at both ends", k == 1 and m.sum() == 80 - 24
          and m[72] and not m[71] and m[127] and not m[128])
    s = np.zeros(n, bool)
    s[:90] = True
    m, _ = D.clear_runs(s, min_run=20, max_gap=8, edge_trim=12)
    check("a run from the frame's edge: cut only where it meets noise", m[0] and m[77] and not m[78])
    s = np.zeros(n, bool)
    s[20:180] = True
    s[100:106] = False
    m, k = D.clear_runs(s, min_run=20, max_gap=8, edge_trim=12)
    check("a 6-column dropout bridged: one run, the dropout itself still no-signal",
          k == 1 and m[99] and not m[100:106].any() and m[106] and m.sum() == 160 - 6 - 24)
    s = np.zeros(n, bool)
    s[10:60] = True
    s[80:95] = True
    m, k = D.clear_runs(s, min_run=20, max_gap=8, edge_trim=12)
    check("a 15-column run past a 20-column gap: left out (under min_run)",
          k == 1 and m[30] and not m[80:95].any())
    g, trim = D.PARTIAL_MAX_GAP_COLUMNS, D.PARTIAL_EDGE_TRIM_COLUMNS
    s = np.ones(1000, bool)
    s[400:400 + g] = False
    m, k = D.clear_runs(s)
    s9 = np.ones(1000, bool)
    s9[400:401 + g] = False
    m9, k9 = D.clear_runs(s9)
    check(f"PARTIAL_MAX_GAP_COLUMNS ({g}): a {g}-column dropout bridged, a {g + 1}-column one ends the run "
          f"(12 columns cut either side), though the clear columns either side would outweigh it",
          g == 49 and k == 1 and m.sum() == 1000 - g and k9 == 2 and m9.sum() == 1000 - (g + 1) - 2 * trim
          and not m9[400 - trim:401 + g + trim].any() and m9[400 - trim - 1] and m9[401 + g + trim],
          (k, m.sum(), k9, m9.sum()))
    s = np.zeros(n, bool)
    s[::7] = True
    m, k = D.clear_runs(s, min_run=20, max_gap=8, edge_trim=12)
    check("scattered single columns (a comb whose gaps are short enough to bridge): nothing, "
          "under PARTIAL_MIN_RUN_DENSITY", not m.any() and k == 0)
    s = np.ones(n, bool)
    s[::4] = False                                    # 75% dense: under 0.80
    m, k = D.clear_runs(s, min_run=20, max_gap=8, edge_trim=12)
    s2 = np.ones(n, bool)
    s2[5:195:6] = False                               # 83% dense, clear at both frame edges
    m2, k2 = D.clear_runs(s2, min_run=20, max_gap=8, edge_trim=12)
    check("a run 75% clear left out, one 83% clear kept", not m.any() and k2 == 1 and m2.sum() == s2.sum(),
          (m.sum(), k2, m2.sum(), s2.sum()))
    m, k = D.clear_runs(np.zeros(n, bool))
    check("no signal: nothing", not m.any() and k == 0)
    s = np.zeros(1000, bool)
    s[:99] = True                                     # 9.9% of the columns, from the frame's edge
    m, _ = D.clear_runs(s)
    s[99] = True                                      # 10%
    m2, _ = D.clear_runs(s)
    check("PARTIAL_MIN_RUN_FRACTION (10%): a run of 9.9% of the columns left out, one of 10% kept",
          D.PARTIAL_MIN_RUN_FRACTION == 0.10 and not m.any() and m2.sum() == 100 - trim, (m.sum(), m2.sum()))
    s = np.zeros(600, bool)
    s[250:] = True
    s[400:404] = False
    s[500:504] = False                                # two bridged dropouts: one run
    _, k0 = D.clear_runs(s)
    m1, k1 = D.clear_runs(s, inside=np.ones(600, bool))
    lim = np.arange(600) < 300
    m2, k2 = D.clear_runs(s, inside=lim)
    m3, k3 = D.clear_runs(s, inside=np.arange(600) < 255)
    t = np.zeros(1000, bool)
    t[:400] = True
    t[600:900] = True                                 # a near run, haze, a far run
    e1, j1 = D.clear_runs(t, near_edge=8)
    t2 = t.copy()
    t2[:9] = False                                    # the near run starting at column 9
    e2, j2 = D.clear_runs(t2, near_edge=8)
    t3 = t.copy()
    t3[1:6] = False                                   # one clear column, a 5-column dropout, the run at 6
    e3, j3 = D.clear_runs(t3, near_edge=8)
    t4 = t.copy()
    t4[:8] = False                                    # the near run starting at column 8
    e4, j4 = D.clear_runs(t4, near_edge=8)
    t5 = t.copy()
    t5[:D.PARTIAL_MAX_GAP_COLUMNS] = False            # the default allowance: starting at 49
    t6 = t.copy()
    t6[:D.PARTIAL_MAX_GAP_COLUMNS + 1] = False        # starting at 50
    _, j5 = D.clear_runs(t5, near_edge=D.PARTIAL_MAX_GAP_COLUMNS)
    _, j6 = D.clear_runs(t6, near_edge=D.PARTIAL_MAX_GAP_COLUMNS)
    check(f"the detector's allowance (PARTIAL_MAX_GAP_COLUMNS, {D.PARTIAL_MAX_GAP_COLUMNS}): a run starting at "
          f"column {D.PARTIAL_MAX_GAP_COLUMNS} taken, at {D.PARTIAL_MAX_GAP_COLUMNS + 1} not", j5 == 1 and j6 == 0, (j5, j6))
    check("near_edge: only the first run, from within 8 columns of the near edge (a speck before it no matter); "
          "a run starting at column 8 taken, at 9 not",
          j1 == 1 and e1[:400 - trim].all() and not e1[400 - trim:].any() and j2 == 0 and not e2.any()
          and j3 == 1 and e3[6 + trim:400 - trim].all() and not e3[:6 + trim].any()
          and j4 == 1 and e4[8 + trim:400 - trim].all(), (j1, j2, j3, j4))
    check("inside (--column-limits): the same run counted once however the limits fall; only columns inside "
          "kept; a run with none inside not counted",
          k0 == k1 == k2 == 1 and m1.sum() == s[250 + trim:].sum() and not (m2 & ~lim).any()
          and m2.sum() == 300 - 250 - trim and k3 == 0 and not m3.any(), (k0, k1, k2, k3, m2.sum()))

    # The haze's chance clear columns (2.7-12% of a dead far field's columns pass the signal test;
    # bridging them used to carry a run into the haze, or sink it under the density test).
    n, edge, trim = 2448, 980, D.PARTIAL_EDGE_TRIM_COLUMNS
    flat, _ = D.clear_runs(haze_frame(n, edge))
    check("flat haze beyond a 40% clear run: kept up to 12 columns short of the haze",
          flat[:edge - trim].all() and not flat[edge - trim:].any())
    for every in (4, 10):
        m, k = D.clear_runs(haze_frame(n, edge, every=every))
        last = int(np.flatnonzero(m)[-1]) if m.any() else -1
        check(f"every {every}th column clear in the haze: the same run kept (not lost, not extended)",
              k == 1 and m[:edge - trim].all() and edge - trim - 1 <= last <= edge - trim + 1, (k, last))
    for p, max_boundary in ((0.12, 0.01), (0.25, 0.01), (0.5, 0.10)):
        rng = np.random.default_rng(5)
        hazy, boundary, short = 0, 0, 0
        for _ in range(400):
            m, _ = D.clear_runs(haze_frame(n, edge, rng, p))
            hazy += bool(m[edge:].any())
            boundary += int(m[edge - trim:edge].sum()) >= trim // 2
            short += not m[:edge - trim].all()
        check(f"{p:.0%} of the haze's columns clear by chance (400 frames): no hazy column flagged in "
              f"{'any' if p < 0.5 else '99%'}, half the 12 cut columns kept in under {max_boundary:.0%}, "
              f"the clear run whole",
              hazy <= (0 if p < 0.5 else 4) and boundary <= max_boundary * 400 and short == 0,
              (hazy, boundary, short))


def test_flat_end():
    print("flat_end")
    n, w, trim = 1000, D.PARTIAL_FLAT_END_WINDOW, D.PARTIAL_EDGE_TRIM_COLUMNS
    rng = np.random.default_rng(2)
    kept = np.zeros(n, bool)
    kept[:700] = True
    rows = 100 + 0.3 * np.arange(n) + rng.normal(0, 0.4, n)
    sharp = 0.3 + 0.05 * rng.random(n)
    m, cut = D.flat_end(kept, rows, sharp)
    check("a sloping line on an edge: not cut", cut == 0 and (m == kept).all())
    flat_rows = rows.copy()
    flat_rows[500:] = flat_rows[500] + rng.normal(0, 0.3, n - 500)
    flat_sharp = sharp.copy()
    flat_sharp[500:] = 0.004 * rng.random(n - 500)
    m, cut = D.flat_end(kept, flat_rows, flat_sharp)
    end = int(np.flatnonzero(m)[-1]) if m.any() else -1
    check("level and on no edge for its last 200 columns: cut back to where it levels (within a third of the "
          "window), 12 columns more", abs((end + 1) - (500 - trim)) <= w // 3 and cut == 700 - (end + 1), (end, cut))
    _, cut0 = D.flat_end(kept, flat_rows, flat_sharp, edge_trim=0)
    _, cut12 = D.flat_end(kept, flat_rows, flat_sharp)
    check("PARTIAL_EDGE_TRIM_COLUMNS (12) more where the cut line ends", cut12 == cut0 + trim, (cut0, cut12))
    m, cut = D.flat_end(kept, flat_rows, sharp)
    check("level but on an edge (a level stretch of waterline): not cut", cut == 0)
    m, cut = D.flat_end(kept, rows, flat_sharp)
    check("on no edge but sloping: not cut", cut == 0)
    lvl = flat_rows.copy()
    lvl[500:] += 0.04 * np.arange(n - 500)
    m, cut = D.flat_end(kept, lvl, flat_sharp)
    lvl2 = flat_rows.copy()
    lvl2[500:] += 0.06 * np.arange(n - 500)
    m2, cut2 = D.flat_end(kept, lvl2, flat_sharp)
    check("PARTIAL_FLAT_END_SLOPE (0.05 px per column): 0.04 cut, 0.06 not", cut > 0 and cut2 == 0, (cut, cut2))
    s2 = flat_sharp.copy()
    s2[500:] = 0.025
    s3 = flat_sharp.copy()
    s3[500:] = 0.035
    _, c2 = D.flat_end(kept, flat_rows, s2)
    _, c3 = D.flat_end(kept, flat_rows, s3)
    check("PARTIAL_FLAT_END_SHARPNESS (0.03): mean sharpness 0.025 cut, 0.035 not", c2 > 0 and c3 == 0, (c2, c3))
    short = np.zeros(n, bool)
    short[:w - 1] = True
    m, cut = D.flat_end(short, flat_rows * 0, flat_sharp * 0)
    check("a line shorter than the window: not cut", cut == 0 and (m == short).all())


def truth_row(c, shift=0.0):
    return 150.0 + shift + 0.15 * c


def field(clear, shift=0.0):
    """An energy field: a ridge on truth_row where `clear`, flat random haze elsewhere."""
    rng = np.random.default_rng(3)
    e = 0.02 + 0.01 * rng.random((ROWS, COLS))
    r = np.arange(ROWS)[:, None]
    ridge = np.exp(-0.5 * ((r - truth_row(np.arange(COLS), shift)[None, :]) / 3.0) ** 2)
    haze = 0.25 + 0.05 * rng.random((ROWS, COLS))
    return np.where(clear[None, :], e + ridge, haze).astype(np.float32)


def v610_detector():
    """The v6.10 detector (before partial lines) from git history, or None without it."""
    try:
        src = subprocess.run(["git", "-C", str(HERE), "show", f"{V610_COMMIT}:waterline_detector_v5.py"],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                             check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None
    name = "waterline_detector_v610"
    mod = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, loader=None))
    mod.__file__ = str(HERE / "waterline_detector_v5.py")
    sys.modules[name] = mod
    exec(compile(src, f"{V610_COMMIT}:waterline_detector_v5.py", "exec"), mod.__dict__)
    return mod


@contextlib.contextmanager
def detector(mod, tmp, energy, partial_min, blend=False, limits=None):
    """process_image with a made-up photo and energy field, everything else the detector's own."""
    profile = replace(mod.CAMERAS["CACO05_C2"], bias_correction_points=(), bias_correction_points_timex=(),
                      band_contrast_floor=0.0)
    names = ("prepare_image", "build_shoreline_energy", "save_energy_debug", "OUTPUT_FOLDER", "EXPORT_CSV",
             "PARTIAL_MIN_FRACTION", "TEMPORAL_FILTER_ENABLED", "COLUMN_LIMITS")
    saved = {k: getattr(mod, k) for k in names if hasattr(mod, k)}

    def prepare(path):
        d = mod.ImageData()
        d.filename = str(path)
        d.original = np.full((ROWS + 100, COLS, 3), 128, np.uint8)
        d.cropped = np.random.default_rng(1).integers(60, 200, (ROWS, COLS, 3)).astype(np.uint8)
        d.crop_top, d.crop_left, d.profile = 50, 0, profile
        return d

    mod.prepare_image = prepare
    mod.build_shoreline_energy = lambda data: (energy, {"envelope": np.ones((ROWS, COLS), bool)})
    mod.save_energy_debug = lambda *a, **k: None
    mod.OUTPUT_FOLDER, mod.EXPORT_CSV = str(tmp), True
    if "PARTIAL_MIN_FRACTION" in saved and partial_min is not None:   # None: the module's own default
        mod.PARTIAL_MIN_FRACTION = partial_min
    mod.TEMPORAL_FILTER_ENABLED = blend
    mod.COLUMN_LIMITS = dict(limits or {})
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(mod, k, v)


def frame_name(epoch):
    t = datetime.fromtimestamp(epoch, tz=timezone.utc)
    return f"{epoch}.{t:%a.%b.%d_%H_%M_%S}.GMT.{t:%Y}.CACO05.c2.timex.jpg"


def run_frame(tmp, clear, partial_min, epoch=EPOCH, mod=D, blend=False, shift=0.0, limits=None, fresh=True):
    tmp = Path(tmp)
    tmp.mkdir(parents=True, exist_ok=True)
    name = frame_name(epoch)
    if fresh:
        for k in ("previous_shoreline", "previous_confidence", "previous_epoch"):
            getattr(mod.STATE, k).pop("CACO05_C2", None)
    buf = io.StringIO()
    with detector(mod, tmp, field(clear, shift), partial_min, blend, limits), contextlib.redirect_stdout(buf):
        _, status = mod.process_image(tmp / name)
    csv_path = tmp / (Path(name).stem + ".csv")
    rows = list(csv.DictReader(open(csv_path))) if csv_path.exists() else []
    return status, rows, buf.getvalue(), csv_path


def state(mod=D):
    k = "CACO05_C2"
    return (mod.STATE.previous_shoreline.get(k), mod.STATE.previous_confidence.get(k),
            mod.STATE.previous_epoch.get(k))


def red_columns(path):
    import cv2
    img = cv2.imread(str(path))
    if img is None:
        return np.zeros(0, int)
    red = (img[:, :, 2] > 180) & (img[:, :, 1] < 90) & (img[:, :, 0] < 90)
    return np.flatnonzero(red.any(axis=0))


def signal_of(rows):
    return np.array([int(r["Has_Signal"]) for r in rows], bool) if rows else np.zeros(0, bool)


def test_detector():
    print("waterline_detector_v5.py: partial lines")
    old = v610_detector()
    c = np.arange(COLS)
    full = np.ones(COLS, bool)
    far = c < COLS // 2                               # the far (high-column) half hazy: c2 after the Oct storm
    trim = D.PARTIAL_EDGE_TRIM_COLUMNS
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        st_on, rows_on, _, p_on = run_frame(tmp / "on", full, 0.30)
        st_off, _, _, p_off = run_frame(tmp / "off", full, 0.0)
        check("signal everywhere: saved, the same CSV with partial lines on and off, no Partial_Line column",
              st_on == st_off == "saved" and p_on.exists() and p_on.read_bytes() == p_off.read_bytes()
              and rows_on and D.PARTIAL_LINE_COLUMN not in rows_on[0])
        if old is None:
            print(f"  skip the v6.10 comparisons (no git history with {V610_COMMIT} here)")
        else:
            st_old, _, _, p_old = run_frame(tmp / "old", full, 0.30, mod=old)
            check("signal everywhere: the same CSV to the byte as the v6.10 detector",
                  st_old == "saved" and p_old.read_bytes() == p_on.read_bytes())

        st, rows, out, _ = run_frame(tmp / "h_off", far, 0.0)
        check("the far half hazy, partial lines off: discarded as before",
              st == "no_signal" and not rows and "DISCARDED (only" in out, (st, out[-300:]))
        st, rows, out, p_csv = run_frame(tmp / "h_on", far, 0.30)
        sig = signal_of(rows)
        col = np.array([float(r["Column"]) for r in rows]) if rows else np.zeros(0)
        row = np.array([float(r["Row_Precise"]) for r in rows]) if rows else np.zeros(0)
        check("partial lines on: saved as a partial line, logged, all columns written",
              st == "partial" and len(rows) == COLS and "PARTIAL LINE" in out and "Finished." in out,
              (st, out[-400:]))
        check("every row marked Partial_Line = 1",
              bool(rows) and all(r.get(D.PARTIAL_LINE_COLUMN) == "1" for r in rows))
        if len(sig):
            check("only the clear near half flagged, from the near edge, cut 12 columns from the haze",
                  sig[:COLS // 2 - trim].all() and not sig[COLS // 2 - trim:].any(),
                  (np.flatnonzero(sig)[[0, -1]] if sig.any() else None))
            err = row[sig] - (truth_row(col[sig]) + 50)
            check("the kept rows on the true waterline (within 1.5 px)", len(err) and np.abs(err).max() <= 1.5,
                  np.abs(err).max() if len(err) else "none")
            kept = np.flatnonzero(sig)
            red = red_columns(p_csv.with_name(p_csv.stem + "_overlay.png"))
            check("the overlay draws the kept stretch only, not the hazy half",
                  len(kept) and len(red) > 0.9 * len(kept) and red.min() <= 2 and red.max() <= kept[-1] + 2,
                  (len(red), red.max() if len(red) else None, kept[-1] if len(kept) else None))

        # the flat end: the solver's line made level and edgeless from column F on (as on c2's frames)
        real_solver = D.solve_shoreline_vectorized

        def flattened(F):
            def solver(energy, profile):
                res = real_solver(energy, profile)
                res.subpixel_shoreline = res.subpixel_shoreline.copy()
                res.subpixel_shoreline[F:] = res.subpixel_shoreline[F]
                res.peak_sharpness = res.peak_sharpness.copy()
                res.peak_sharpness[F:] = 0.0
                return res
            return solver

        try:
            D.solve_shoreline_vectorized = flattened(200)
            st, rows, out, _ = run_frame(tmp / "flat200", far, None)
            sig = signal_of(rows)
            last = int(np.flatnonzero(sig)[-1]) if sig.any() else -1
            check("a flat end from column 200 (88 of the kept columns): cut where it levels (within a third of "
                  "the window), 12 columns more, the rest saved as a partial line, logged",
                  st == "partial" and abs((last + 1) - (200 - trim)) <= D.PARTIAL_FLAT_END_WINDOW // 3
                  and "a flat end" in out, (st, last))
            D.solve_shoreline_vectorized = flattened(240)
            st, rows, _, _ = run_frame(tmp / "flat240", far, None)
            check("a flat end shorter than the window (48 kept columns): left to the 12-column cut",
                  st == "partial" and signal_of(rows)[:COLS // 2 - trim].all(), st)
            D.solve_shoreline_vectorized = flattened(150)
            st, rows, out, _ = run_frame(tmp / "flat150", far, None)
            check("a flat end from column 150: under 30% left, discarded, saying so",
                  st == "no_signal" and not rows and "after its flat end" in out, (st, out[-300:]))
        finally:
            D.solve_shoreline_vectorized = real_solver

        # haze at the near end or in the middle: the solver's line starts free in it, and a stronger edge
        # beside the clear stretch can take it, which no confidence shows -- discarded, as v6.10 did
        near_half = c >= COLS // 2
        st, rows, out, _ = run_frame(tmp / "near", near_half, None)
        check("the near half hazy: discarded, as v6.10 (no clear run from the near edge)",
              st == "no_signal" and not rows and "clear run from the near edge covers 0.0%" in out, (st, out[-300:]))
        near35 = c >= int(0.35 * COLS)                # 65% clear: whole-line confidence under 0.9
        st, rows, out, _ = run_frame(tmp / "near35", near35, None)
        check("the near 35% hazy (65% clear): discarded for its confidence, as v6.10 -- never retried as a "
              "partial line", st == "low_confidence" and not rows and "PARTIAL LINE" not in out, st)
        middle = (c < 200) | (c >= 500)               # 50%: clear near third and far sixth, haze between
        st, rows, out, _ = run_frame(tmp / "middle", middle, None)
        sig = signal_of(rows)
        check("a hazy middle: the clear run from the near edge kept, the one beyond the haze not",
              st == "partial" and len(sig) and sig[:200 - trim].all() and not sig[200 - trim:].any(), st)
        speck = far.copy()
        speck[1:5] = False                            # a 4-column dropout right at the near edge
        st, rows, _, _ = run_frame(tmp / "speck", speck, None)
        sig = signal_of(rows)
        check("a dropout at the near edge (a single clear column before it): still the run from the near edge",
              st == "partial" and len(sig) and sig[5 + trim:COLS // 2 - trim].all() and not sig[COLS // 2:].any(), st)
        # the review's case: near 30% hazy, a ridge 1.3x the waterline's 80 px off in columns 192-389
        # (v6.10 discarded it; retried as a partial line it kept 205 of 408 columns on that ridge)
        rng = np.random.default_rng(3)
        r = np.arange(ROWS)[:, None]
        ridge = np.exp(-0.5 * ((r - truth_row(c)[None, :]) / 3.0) ** 2)
        other = 1.3 * np.exp(-0.5 * ((r - truth_row(c)[None, :] - 80) / 3.0) ** 2) * ((c >= 192) & (c < 390))[None, :]
        clear = c >= 180
        energy = np.where(clear[None, :], 0.02 + 0.01 * rng.random((ROWS, COLS)) + np.maximum(ridge, other),
                          0.25 + 0.05 * rng.random((ROWS, COLS))).astype(np.float32)
        (tmp / "edge").mkdir()
        buf = io.StringIO()
        for k in ("previous_shoreline", "previous_confidence", "previous_epoch"):
            getattr(D.STATE, k).pop("CACO05_C2", None)
        with detector(D, tmp / "edge", energy, None), contextlib.redirect_stdout(buf):
            _, st = D.process_image(tmp / "edge" / frame_name(EPOCH))
        check("near haze with a stronger edge beside the clear stretch (the review's case): not saved",
              st in ("no_signal", "low_confidence") and not list((tmp / "edge").glob("*.csv")), st)
        if old is not None:
            st85, _, _, p85 = run_frame(tmp / "c85", c >= int(0.15 * COLS), None)
            st85o, _, _, p85o = run_frame(tmp / "c85old", c >= int(0.15 * COLS), None, mod=old)
            stn, _, _, _ = run_frame(tmp / "near35old", near35, None, mod=old)
            check("v6.10 the same: the near 35% hazy discarded for its confidence; 85% clear passes whole, "
                  "the same CSV to the byte", stn == "low_confidence" and st85 == st85o == "saved"
                  and p85.read_bytes() == p85o.read_bytes(), (stn, st85, st85o))

        # whole -> partial -> whole, the temporal blend on (frames 15 min apart; it blends up to 35 min)
        seq = [(full, 0.0), (far, 6.0), (full, 3.0)]
        st1, _, _, _ = run_frame(tmp / "seq", seq[0][0], 0.30, EPOCH, blend=True, shift=seq[0][1])
        after_whole = state()
        check("a whole frame is the next frame's blend reference (line, confidence, time)",
              st1 == "saved" and after_whole[0] is not None and after_whole[1] is not None
              and after_whole[2] == EPOCH, (st1, after_whole[2]))
        st2, _, _, _ = run_frame(tmp / "seq", seq[1][0], 0.30, EPOCH + 900, blend=True, shift=seq[1][1], fresh=False)
        after_partial = state()
        check("a partial line is not: the line, confidence and time it leaves are the whole frame's",
              st2 == "partial" and all(a is b for a, b in zip(after_whole[:2], after_partial[:2]))
              and after_partial[2] == EPOCH, (st2, after_partial[2]))
        st3, _, _, p3 = run_frame(tmp / "seq", seq[2][0], 0.30, EPOCH + 1800, blend=True, shift=seq[2][1], fresh=False)
        if old is not None:
            for i, (clear, shift) in enumerate(seq):
                run_frame(tmp / "seq_old", clear, 0.30, EPOCH + 900 * i, mod=old, blend=True, shift=shift,
                          fresh=(i == 0))
            p3_old = tmp / "seq_old" / p3.name
            check("the third frame: the same CSV to the byte as the v6.10 detector's (which discarded the second)",
                  st3 == "saved" and p3_old.exists() and p3_old.read_bytes() == p3.read_bytes())

        # --column-limits, the near half clear
        for lim, want in (((0.5, 1.0), "no_signal"), ((0.4, 1.0), "no_signal"), ((0.0, 0.6), "partial"),
                          ((0.1, 1.0), "partial")):
            st, rows, out, _ = run_frame(tmp / f"lim{lim[0]}_{lim[1]}", far, 0.30, limits={"CACO05_C2": lim})
            sig = signal_of(rows)
            inside = (c / (COLS - 1) >= lim[0]) & (c / (COLS - 1) <= lim[1])
            good = st == want and (not rows if want == "no_signal" else
                                   sig.mean() >= D.PARTIAL_MIN_FRACTION and not (sig & ~inside).any())
            check(f"--column-limits c2={lim[0]}:{lim[1]} with the near half clear: {want}"
                  + (" (judged on the columns inside the limits)" if want == "no_signal" else
                     ", at least 30% of the columns flagged, all inside the limits"),
                  good and (want == "partial" or "inside --column-limits" in out), (st, sig.mean() if len(sig) else None))

        # PARTIAL_MIN_FRACTION itself, at the module's default
        st28, _, out28, _ = run_frame(tmp / "b28", c < 180, None)      # 28% kept after the cut
        st31, rows31, _, _ = run_frame(tmp / "b31", c < 200, None)     # 31.3%
        sig31 = signal_of(rows31)
        check("PARTIAL_MIN_FRACTION (30%, the default): 28% of the columns kept is discarded, 31.3% saved",
              D.PARTIAL_MIN_FRACTION == 0.30 and st28 == "no_signal" and "covers 28.0%, need 30%" in out28
              and st31 == "partial" and abs(sig31.mean() - 188 / 600) < 1e-9, (st28, st31))

        bits = np.zeros(COLS, bool)
        for a in range(0, COLS, 100):
            bits[a:a + 40] = True                     # 40% of columns, in 40-column stretches
        st, rows, out, _ = run_frame(tmp / "bits", bits, 0.30)
        check("only scattered short stretches clear: discarded, saying the clear run from the near edge falls short",
              st == "no_signal" and not rows and "clear run from the near edge covers" in out, (st, out[-300:]))

    lines = Path(D.__file__).read_text()
    check("v6.11 in the header and in VERSION (the run banner)",
          "Version 6.11" in lines and "PARTIAL LINES" in lines and D.VERSION == "6.11", D.VERSION)


def write_frame(folder, ep, signal, partial=False, shift=0.0):
    t = datetime.fromtimestamp(ep, tz=timezone.utc)
    out = ["Column,Row,Confidence,Row_Precise,Peak_Sharpness,Has_Signal" + (",Partial_Line" if partial else "")]
    for c in range(300):
        row = 500 + shift + c / 10.0 + (c * 7 % 13) / 4.0      # not dead straight (the extractor drops those)
        out.append(f"{c},{int(row)},0.95,{row:.3f},0.5,{int(bool(signal[c]))}" + (",1" if partial else ""))
    (folder / f"{ep}.{t:%a.%b.%d_%H_%M_%S}.GMT.{t:%Y}.CACO05.c2.timex.csv").write_text("\n".join(out) + "\n")


def write_tide(d, start, n, step_m):
    tide = ["time,level"] + [f"{datetime.fromtimestamp(start + k * 600, tz=timezone.utc):%Y-%m-%d %H:%M:%S},"
                             f"{0.2 + step_m * k:.4f}" for k in range(n)]
    (d / "tide.csv").write_text("\n".join(tide) + "\n")


def extract(d, extra=()):
    out = d / "contours.csv"
    if out.exists():
        out.unlink()
    r = subprocess.run([sys.executable, str(HERE / "extract_elevation_contours.py"), str(d / "tide.csv"),
                        "--processed-dir", str(d / "processed"), "--output", str(out)] + list(extra),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
    used = set()
    if out.exists():
        with open(out) as f:
            used = {row["source_file"][:10] for row in csv.DictReader(f)}
    return r.returncode, used, r.stdout


def run_of(a, b):
    s = np.zeros(300, bool)
    s[a:b] = True
    return s


def test_extractor():
    print("extract_elevation_contours.py: the partial-line floor")
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "processed").mkdir()
        frames = {
            "full": (EPOCH, np.ones(300, bool), False),
            "partial40": (EPOCH + 1800, run_of(150, 270), True),      # 40%: a partial line's clear run
            "partial27": (EPOCH + 3600, run_of(200, 280), True),      # 27%: 22.7% (68/300) left after the fragment trim
            "full62": (EPOCH + 5400, run_of(100, 286), False),        # 62%: a whole line, 58% (174/300) left after the trim
            # whole lines with few columns flagged, NOT partial lines: detected with a lower
            # --min-signal-fraction (process_chelsea.py, detect_original_view.py: 0.15), or cut by
            # --column-limits; the share of Has_Signal columns does not make a line partial
            "near27": (EPOCH + 7200, run_of(0, 80), False),          # 23% (69/300) left after the fragment trim
            "scatter59": (EPOCH + 9000, run_of(0, 115) | (np.arange(300) % 3 == 0) & (np.arange(300) >= 115), False),
        }
        for ep, sig, partial in frames.values():
            write_frame(d / "processed", ep, sig, partial)
        write_tide(d, EPOCH - 3600, 30, 0.002)
        ep = {k: str(v[0]) for k, v in frames.items()}

        rc, used, out = extract(d)
        check("the defaults: the whole line and the 40% partial line used, the one left with 22.7% not",
              rc == 0 and ep["full"] in used and ep["partial40"] in used and ep["partial27"] not in used,
              (rc, used, out[-600:]))
        check("the log counts the partial lines used", "(1 partial line(s)" in out, out[-400:])
        check("an unmarked line with 59% of its columns flagged (under 0.50 left after the fragment trim): judged "
              "by --min-coverage 0.50 as before, not used", ep["scatter59"] not in used, used)
        rc, used, _ = extract(d, ["--partial-min-coverage", "0.5"])
        check("--partial-min-coverage 0.5: the 40% partial line not used", ep["partial40"] not in used, used)
        rc, used, _ = extract(d, ["--min-coverage", "0.6"])
        check("a whole line (no Partial_Line) still judged by --min-coverage: 0.6 drops the 62% line, "
              "the partial line kept by its own floor",
              ep["full62"] not in used and ep["partial40"] in used and ep["full"] in used, used)
        rc, used, _ = extract(d, ["--partial-min-coverage", "0"])
        check("--partial-min-coverage 0: no coverage floor for partial lines", ep["partial27"] in used, used)
        rc, used, _ = extract(d, ["--min-coverage", "0.15"])
        check("--min-coverage 0.15 (the Chelsea runs): an unmarked 27% line (23% left) used as before, and a partial "
              "line never held to more than --min-coverage", ep["near27"] in used and ep["partial27"] in used, used)
        rc, used, _ = extract(d, ["--min-coverage", "0"])
        check("--min-coverage 0 turns the floor off for every line, partial ones too",
              all(ep[k] in used for k in ("near27", "partial27", "scatter59")), used)
        rc, used, out = extract(d, ["--include-no-signal"])
        check("--include-no-signal: the log says the partial lines' no-signal columns were kept, not dropped",
              "kept by --include-no-signal" in out and "the hazy rest dropped" not in out, out[-400:])

    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "processed").mkdir()
        # day 1: 8 partial lines whose rows zig-zag while the tide rises steadily (the tide-agreement
        # filter drops them); day 2: one whole line
        for k in range(8):
            write_frame(d / "processed", EPOCH + 43200 + k * 1800, run_of(0, 120), True, shift=10 * (k % 2))
        write_frame(d / "processed", EPOCH + 43200 + 86400, np.ones(300, bool))
        write_tide(d, EPOCH + 43200 - 3600, 170, 0.1)
        rc, used, out = extract(d)
        line = next((ln for ln in out.splitlines() if ln.startswith("Frames used")), "")
        n_used = int(line.split(":")[1].split()[0]) if line else -1
        check("the log's partial count is of the frames used, after the tide-agreement filter",
              rc == 0 and n_used == len(used) and ("partial line" not in line if n_used == 1 else True)
              and len(used) < 9, (line, used))

    check("the detector and the extractor name the same Partial_Line column",
          X.PARTIAL_LINE_COLUMN == D.PARTIAL_LINE_COLUMN)


def test_ground_truth_tools():
    print("ground-truth tools: a partial line's hazy rows are not scored silently")
    import batch_bias_analysis as B
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "processed").mkdir()
        (d / "gt").mkdir()
        whole, part = EPOCH, EPOCH + 1800
        write_frame(d / "processed", whole, np.ones(300, bool))
        write_frame(d / "processed", part, run_of(150, 270), True)
        for ep in (whole, part):
            src = next((d / "processed").glob(f"{ep}.*.csv"))
            (d / "gt" / src.name).write_text(src.read_text().split("\n", 1)[0].split(",Partial_Line")[0] + "\n"
                                             + "\n".join(f"{c},{500 + c // 10}" for c in range(300)) + "\n")
        saved = B.PROCESSED_FOLDER, B.GROUND_TRUTH_FOLDER
        B.PROCESSED_FOLDER, B.GROUND_TRUTH_FOLDER = str(d / "processed"), str(d / "gt")
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                B.main()
        finally:
            B.PROCESSED_FOLDER, B.GROUND_TRUTH_FOLDER = saved
        out = buf.getvalue()
        check("batch_bias_analysis.py: the partial line listed as skipped, the whole one scored",
              "Images compared     : 1" in out and "1 detector output(s) are partial lines" in out
              and f"  - {next((d / 'processed').glob(f'{part}.*.csv')).stem}" in out, out[-500:])
        for tool in ("compare_to_ground_truth.py", "bias_analysis.py"):
            warned = []
            for ep in (whole, part):
                det = next((d / "processed").glob(f"{ep}.*.csv"))
                r = subprocess.run([sys.executable, str(HERE / tool), str(det), str(d / "gt" / det.name)],
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
                warned.append(r.returncode == 0 and "is a partial line (Partial_Line)" in r.stdout)
            check(f"{tool}: a WARNING on the partial line, none on the whole one", warned == [False, True], warned)

    import sweep_profile_scale as S
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        img = d / frame_name(EPOCH)
        img.write_bytes(b"not read: prepare_image is the test's")
        gt = d / "gt.csv"
        gt.write_text("Column,Row\n" + "\n".join(f"{c},{truth_row(c) + 50:.2f}" for c in range(COLS)) + "\n")
        outs = {}
        for name, clear in (("whole", np.ones(COLS, bool)), ("far half hazy", np.arange(COLS) < COLS // 2)):
            buf = io.StringIO()
            argv = sys.argv
            sys.argv = ["sweep_profile_scale.py", str(img), str(gt), "c2", "--param", "vertical_penalty",
                        "--values", "0.35"]
            try:
                with detector(D, d / "unused", field(clear), None), contextlib.redirect_stdout(buf):
                    S.main()
            finally:
                sys.argv = argv
            outs[name] = buf.getvalue()
        check("sweep_profile_scale.py: a partial line not scored, a whole one scored",
              "partial line, not scored" in outs["far half hazy"] and "Best:" not in outs["far half hazy"]
              and "Best:" in outs["whole"], {k: v[-200:] for k, v in outs.items()})


def test_syntax_38():
    print("Python 3.8")
    for f in ("waterline_detector_v5.py", "extract_elevation_contours.py", "test_partial_lines.py",
              "batch_bias_analysis.py", "compare_to_ground_truth.py", "bias_analysis.py",
              "sweep_profile_scale.py"):
        try:
            ast.parse((HERE / f).read_text(), filename=f, feature_version=(3, 8))
            good = True
        except SyntaxError as exc:
            good = False
            print(f"    {exc}")
        check(f"{f} parses as Python 3.8", good)


def main():
    test_maximal_segments()
    test_clear_runs()
    test_flat_end()
    test_detector()
    test_extractor()
    test_ground_truth_tools()
    test_syntax_38()
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        raise SystemExit(1)
    print("all partial-line tests passed")


if __name__ == "__main__":
    main()

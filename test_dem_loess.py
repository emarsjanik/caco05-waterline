#!/usr/bin/env python3
"""
Tests For The Loess DEM And Its Prediction Error
===================================================
Plain python3 (no pytest needed on the station):

    python3 test_dem_loess.py [--keep DIR]

WHAT. dem_loess.py, the opt-in scale-controlled quadratic loess surface of
the DEM's points with its prediction error (Plant, Holland & Puleo 2002),
and its hooks (dem_from_contours.py --loess, survey_products.py
--dem-loess, asc_to_geotiff.py):
  * a quadratic beach is reproduced exactly from noise-free lines (the
    local model is quadratic);
  * a known synthetic beach (curved profile, a gentle alongshore wave)
    seen by waterlines with a water-level error common to each frame and
    a detection scatter along each line is recovered within its error
    bars: unbiased, and the truth within 2 err in >= 90% of the cells of
    every seed;
  * err is the size of the error it predicts: on fixed lines with only a
    level error common to each frame, many times over, the cells' RMS
    error over err is ~1 and the truth within err in ~68% (an err too
    large, as from the n/(n - 6) rule of an unweighted fit, fails);
  * the error grows away from the data: into a gap between the lines
    (nmse rising cell by cell from the last line to the middle of the
    gap, which is left nodata) and beyond the lowest line or the end of
    the lines (left nodata, even with no nmse limit; with extrapolation
    on, for this test only, err rising with distance);
  * nmse is summed per frame: with one point per frame it is Plant's
    sum a_j^2 (checked against a plain weighted least-squares solve), and
    it does not shrink when each line is sampled four times as densely;
  * the window's response: what it keeps of a sine by wavelength;
  * dem_from_contours.py: without --loess every grid, <stem>_info.json
    and the page are what they were; with it the DEM's grids and page
    are the same to the byte and pixel, the loess grids are on the DEM's
    grid, info.json gains 'loess', and bad --loess-scales are refused;
  * dem_loess.py evaluate against a survey grid and survey points (an
    offset planted in the points is the reported median and is removed
    in the 'median removed' shares), and its refusal of loess grids an
    earlier build left beside the DEM;
  * asc_to_geotiff.py: the loess surface carries the NAVD88 vertical keys,
    its error grids do not;
  * survey_products.py: without --dem-loess the dem step's commands and
    files are what they were; with it, --loess, the loess GeoTIFFs, the
    stamp, provenance.json, README and rebuild line saying so; a rebuild
    without it removes the earlier build's loess files.

HOW. Synthetic waterlines with a planted truth, in a temporary folder
(--keep DIR to look at it).

STATION FILES. None are read.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import dem_loess as dl                      # noqa: E402

FAILED = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILED.append(msg)


def run(args, d):
    r = subprocess.run([sys.executable] + [str(a) for a in args], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
    return r.returncode, r.stdout


def read_asc(path):
    a = np.loadtxt(path, skiprows=6)
    return np.where(a == -9999, np.nan, a)


# ---------------------------------------------------------------------
# A synthetic beach and its waterlines
# ---------------------------------------------------------------------

E0, N0 = 420100.0, 4638200.0
# z = Z_TOP + SLOPE x + CURV x^2 + AMP sin(2 pi y / WAVE): x seaward (east), y alongshore (north), m.
# A curved profile (a quadratic, which the local model follows exactly) and an alongshore wave 10 times
# the alongshore half-width (kept at 99.9%), so the truth is all but unsmoothed at the default scales.
Z_TOP, SLOPE, CURV, AMP, WAVE = 1.5, -0.08, 0.0008, 0.1, 200.0


def truth(x, y, amp=AMP):
    return Z_TOP + SLOPE * x + CURV * x * x + amp * np.sin(2 * np.pi * y / WAVE)


def lines(seed, nf=60, lo=-0.3, hi=1.4, sig_e=0.04, sig_x=0.3, gap=None, step=0.5, length=200.0, amp=AMP):
    """nf waterlines at water levels uniform in [lo, hi] m (none inside `gap`), each where the beach is at
    the TRUE level, level + e (e ~ N(0, sig_e): an error of the level, common to the frame), traced with a
    cross-shore scatter sig_x (m) at every point (the detection), every `step` m alongshore. Each point is
    given the frame's level, as the station's rows are. -> E, N, Z, frame names, levels"""
    rng = np.random.RandomState(seed)
    E, N, Z, F, H = [], [], [], [], []
    k = 0
    while k < nf:
        h = rng.uniform(lo, hi)
        if gap and gap[0] < h < gap[1]:
            continue
        e = rng.normal(0, sig_e)
        y = np.arange(0, length, step)
        c = Z_TOP - (h + e - amp * np.sin(2 * np.pi * y / WAVE))
        x = (-SLOPE - np.sqrt(SLOPE * SLOPE - 4 * CURV * c)) / (2 * CURV)
        x = x + rng.normal(0, sig_x, len(x))
        E.append(E0 + x)
        N.append(N0 + y)
        Z.append(np.full(len(x), h))
        F += [f"f{k:03d}"] * len(x)
        H.append(h)
        k += 1
    return np.concatenate(E), np.concatenate(N), np.concatenate(Z), np.array(F), np.array(H)


def grid_of(E, N, cell=2.0):
    """The grid build_grid puts the points on: (e0, n0, ncols, nrows)."""
    e0 = np.floor(E.min() / cell) * cell
    n0 = np.floor(N.min() / cell) * cell
    return e0, n0, int(np.ceil((E.max() - e0) / cell)) + 1, int(np.ceil((N.max() - n0) / cell)) + 1


def cell_xy(e0, n0, ncols, nrows, cell=2.0):
    gx = e0 + (np.arange(ncols) + 0.5) * cell
    gy = n0 + (np.arange(nrows) + 0.5) * cell
    GX, GY = np.meshgrid(gx, gy)
    return GX - E0, GY - N0


def fit(E, N, Z, F, cell=2.0, **kw):
    e0, n0, ncols, nrows = grid_of(E, N, cell)
    lo = dl.loess_grid(E, N, Z, F, e0, n0, ncols, nrows, cell, **kw)
    X, Y = cell_xy(e0, n0, ncols, nrows, cell)
    return lo, X, Y


# ---------------------------------------------------------------------
# The estimator
# ---------------------------------------------------------------------

def test_exact_quadratic():
    print("loess_grid: a quadratic beach from noise-free lines is reproduced exactly")
    E, N, Z, F, _ = lines(1, nf=40, sig_e=0.0, sig_x=0.0, amp=0.0)
    lo, X, Y = fit(E, N, Z, F)
    ok = np.isfinite(lo["z"])
    d = np.abs(lo["z"] - truth(X, Y, amp=0.0))[ok]
    check(ok.sum() > 500 and d.max() < 2e-3 and np.nanmax(lo["err"]) < 2e-3,
          f"{int(ok.sum())} cells, largest |loess - truth| {d.max():.1e} m, largest err {np.nanmax(lo['err']):.1e} m "
          f"(the points are thinned to their mean per half-cell, which moves them off the curve by < 1 mm)")


def test_recovery():
    print("loess_grid: a known synthetic beach recovered within its error bars")
    for seed in range(4):
        E, N, Z, F, _ = lines(seed)
        lo, X, Y = fit(E, N, Z, F)
        ok = np.isfinite(lo["z"])
        d = (lo["z"] - truth(X, Y))[ok]
        e = lo["err"][ok]
        in1, in2 = float(np.mean(np.abs(d) <= e)), float(np.mean(np.abs(d) <= 2 * e))
        rms = float(np.sqrt(np.mean(d ** 2)))
        check(ok.sum() > 1000 and abs(np.median(d)) < 0.015 and rms < 0.03 and in2 >= 0.9 and in1 >= 0.5
              and 0.3 < rms / np.median(e) < 1.5 and np.median(lo["skill"][ok]) > 0.7,
              f"seed {seed}: {int(ok.sum())} cells, loess - truth median {np.median(d):+.4f} m, RMS {rms:.3f} m; "
              f"median err {np.median(e):.3f} m; truth within err {100 * in1:.0f}%, within 2 err {100 * in2:.0f}%; "
              f"median skill {np.median(lo['skill'][ok]):.2f}")
    s = lo["settings"]
    ax = s["axis_alongshore_deg"] % 180.0
    check(s["along_m"] == 20 and s["across_m"] == 5 and s["order"] == 2 and min(ax, 180.0 - ax) < 2,
          f"default scales 20 x 5 m, quadratic; shore axis {s['axis_alongshore_deg']} deg (the lines run north)")


def test_calibration():
    print("loess_grid: err the size of the error it predicts (a level error common to each frame)")
    # the same 30 lines 20 times over, each time with a new level error per frame (sd 0.04 m) and nothing else
    E, N, Z, F, _ = lines(0, nf=30, sig_e=0.0, sig_x=0.0)
    fid = np.unique(F, return_inverse=True)[1]
    rng = np.random.RandomState(1)
    zs, es = [], []
    for _ in range(20):
        lo, X, Y = fit(E, N, Z + rng.normal(0, 0.04, fid.max() + 1)[fid], F)
        zs.append(lo["z"])
        es.append(lo["err"])
    zs, es = np.array(zs), np.array(es)
    ok = np.all(np.isfinite(zs), 0)
    d, e = zs[:, ok] - truth(X, Y)[ok], es[:, ok]
    r = np.sqrt(np.mean(d ** 2, 0)) / np.sqrt(np.mean(e ** 2, 0))
    in1, in2 = float(np.mean(np.abs(d) <= e)), float(np.mean(np.abs(d) <= 2 * e))
    check(ok.sum() > 200 and 0.8 < np.median(r) < 1.25 and 0.55 < in1 < 0.8 and in2 > 0.85,
          f"{int(ok.sum())} cells x 20 realisations: RMS error / RMS err per cell median {np.median(r):.2f} "
          f"(p10 {np.percentile(r, 10):.2f}, p90 {np.percentile(r, 90):.2f}); truth within err {100 * in1:.0f}%, "
          f"within 2 err {100 * in2:.0f}%")


def test_error_grows_into_gap():
    print("loess_grid: the error grows into a gap between the lines")
    # no line between water levels 0.3 and 0.9 m: a gap of ~9 m cross-shore on this profile
    E, N, Z, F, _ = lines(0, gap=(0.3, 0.9))
    lo, X, Y = fit(E, N, Z, F, max_nmse=None)
    tz = truth(X, Y)
    mid = (Y > 40) & (Y < 160)
    bands = [(0.9, 1.0), (0.8, 0.9), (0.7, 0.8)]
    nm = [float(np.nanmedian(lo["raw"]["nmse"][mid & (tz >= a) & (tz < b)])) for a, b in bands]
    er = [float(np.nanmedian(lo["raw"]["err"][mid & (tz >= a) & (tz < b)])) for a, b in bands]
    dense = float(np.nanmedian(lo["nmse"][mid & (tz >= 1.0) & (tz < 1.3)]))
    check(dense < nm[0] < nm[1] < nm[2] and er[0] < er[1] < er[2],
          "from the lines into the gap (beach 1.0-0.7 m): nmse " + " < ".join(f"{v:.2f}" for v in nm)
          + f" (among the lines {dense:.2f}); err " + " < ".join(f"{v:.3f}" for v in er) + " m")
    centre = mid & (tz > 0.5) & (tz < 0.7)
    check(centre.sum() > 50 and not np.isfinite(lo["z"][centre]).any(),
          f"the middle of the gap ({int(centre.sum())} cells, beach 0.5-0.7 m) is left nodata")
    lo2, _, _ = fit(E, N, Z, F)
    check(lo2["counts"]["nmse_over"] > 0 and lo2["counts"]["filled"] < lo["counts"]["filled"]
          and np.nanmax(lo2["nmse"]) <= dl.LOESS_MAX_NMSE,
          f"with the default limit (nmse <= {dl.LOESS_MAX_NMSE}) {lo2['counts']['nmse_over']} more cells at the "
          f"gap's edges are left out")


def test_error_grows_beyond_data():
    print("loess_grid: the error grows beyond the lowest line; no extrapolation in a product")
    # a straight shore (no alongshore wave), dense lines from 0.2 to 1.4 m, and a grid that runs 10 m on seaward
    # and 10 m on past the north end of the lines
    E, N, Z, F, H = lines(2, nf=200, lo=0.2, hi=1.4, amp=0.0)
    cell = 2.0
    e0, n0, ncols, nrows = grid_of(E, N, cell)
    ncols += 5
    nrows += 5
    lo = dl.loess_grid(E, N, Z, F, e0, n0, ncols, nrows, cell, max_nmse=None)
    lox = dl.loess_grid(E, N, Z, F, e0, n0, ncols, nrows, cell, extrapolate=True, max_nmse=None)
    X, Y = cell_xy(e0, n0, ncols, nrows, cell)
    tz = truth(X, Y, amp=0.0)
    xlow = (-SLOPE - np.sqrt(SLOPE * SLOPE - 4 * CURV * (Z_TOP - H.min()))) / (2 * CURV)   # the lowest line, m
    mid = (Y > 40) & (Y < 160)
    bands = [(-4.0, -2.0), (-2.0, 0.0), (0.0, 2.0), (2.0, 4.0)]
    sel = [mid & (X > xlow + a) & (X <= xlow + b) for a, b in bands]
    north = (Y > N.max() - N0 + 2.0) & (X < xlow - 4.0)
    check(lo["counts"]["one_sided"] > 0 and not np.isfinite(lo["z"][sel[3]]).any()
          and not np.isfinite(lo["z"][mid & (X > xlow + 2.0)]).any()
          and north.sum() > 20 and not np.isfinite(lo["z"][north]).any(),
          f"no nmse limit, and even so nodata more than 2 m seaward of the lowest line and past the north end of "
          f"the lines ({lo['counts']['one_sided']} cells with points on one side only left out)")
    er = [float(np.nanmedian(lox["err"][m])) for m in sel]
    nm = [float(np.nanmedian(lox["nmse"][m])) for m in sel]
    ae = [float(np.nanmedian(np.abs(lox["z"] - tz)[m])) for m in sel]
    check(all(np.isfinite(lox["z"][m]).mean() > 0.5 for m in sel) and er[0] < er[1] < er[2] < er[3]
          and nm[0] < nm[1] < nm[2] < nm[3] and nm[3] > 1 and lox["counts"]["one_sided"] == 0,
          "extrapolated (tests only), from 4 m inside to 4 m beyond the lowest line in 2 m steps: err "
          + " < ".join(f"{v:.4f}" for v in er) + " m, nmse " + " < ".join(f"{v:.2f}" for v in nm)
          + " (|loess - truth| " + ", ".join(f"{v:.3f}" for v in ae) + " m)")
    far = mid & (X > xlow + 4.0)
    check(far.sum() > 0 and not np.isfinite(lox["z"][far]).any(),
          "further out, too few lines in the window to say anything: nodata even so")


def test_nmse_per_frame():
    print("loess_grid: nmse summed per frame (Plant's sum a_j^2 with one point per frame)")
    rng = np.random.RandomState(5)
    n = 40
    E = E0 + rng.uniform(-4, 4, n)
    N = N0 + rng.uniform(-15, 15, n)
    Z = 0.5 + 0.1 * (E - E0) + rng.normal(0, 0.05, n)
    F = np.array([f"p{i}" for i in range(n)])
    cell = 2.0
    e0, n0 = E0 - 1.0, N0 - 1.0                   # one cell centred on (E0, N0)
    lo = dl.loess_grid(E, N, Z, F, e0, n0, 1, 1, cell, along=20.0, across=5.0, thin_bin=1e-6, max_nmse=None)
    # the same estimate by a plain weighted least squares at the node
    ec, nc, ua, va = dl.beach_axes(E, N)
    du = ((E - E0) * ua[0] + (N - N0) * ua[1]) / 20.0
    dv = ((E - E0) * va[0] + (N - N0) * va[1]) / 5.0
    w = dl.tricube(np.sqrt(du ** 2 + dv ** 2))
    X = dl.design(du, dv, 2)
    A = (X * w[:, None]).T @ X
    a = w * (X @ np.linalg.solve(A, np.eye(6)[:, 0]))
    check(abs(lo["raw"]["z"][0, 0] - float(a @ Z)) < 1e-9 and abs(lo["raw"]["nmse"][0, 0] - float(np.sum(a ** 2))) < 1e-9,
          f"one point per frame: z_hat = sum a_j z_j ({float(a @ Z):.4f} m) and nmse = sum a_j^2 "
          f"({float(np.sum(a ** 2)):.4f}), as a plain solve gives")
    # the same lines sampled four times as densely along each line: the frames, not the points, count
    E1, N1, Z1, F1, _ = lines(3, step=1.0)
    E4, N4, Z4, F4, _ = lines(3, step=0.25)
    l1, X1, Y1 = fit(E1, N1, Z1, F1, thin_bin=0.25)
    l4, _, _ = fit(E4, N4, Z4, F4, thin_bin=0.25)
    both = np.isfinite(l1["nmse"]) & np.isfinite(l4["nmse"])
    r = float(np.median(l4["nmse"][both] / l1["nmse"][both]))
    check(both.sum() > 500 and 0.8 < r < 1.25,
          f"each line sampled every 0.25 m instead of 1 m (4x the points, no thinning): nmse x {r:.2f} (per point "
          f"it would fall ~4x)")


def test_response():
    print("response: what the window keeps of a sine, by wavelength")
    L = 5.0
    p90, p50, side = dl.response_table(L)
    check(abs(p90 / L - 1.84) < 0.03 and abs(p50 / L - 1.1) < 0.03 and side < 0.1
          and dl.response(20 * L, L) > 0.99 and abs(dl.response(0.5 * L, L)) < 0.1,
          f"half-width {L:g} m: 90% of a {p90:.1f} m wave, 50% of a {p50:.1f} m one, at most {100 * side:.0f}% under "
          f"{0.8 * L:g} m; {dl.response(20 * L, L):.3f} of a {20 * L:g} m one (the docstring's numbers)")
    p90_1, _, _ = dl.response_table(L, order=1)
    check(p90_1 > 2 * p90, f"a plane (order 1) needs a {p90_1:.1f} m wave to keep 90%: it flattens the curvature")


# ---------------------------------------------------------------------
# dem_from_contours.py --loess, dem_loess.py evaluate, asc_to_geotiff.py
# ---------------------------------------------------------------------

HEADER = "source_file,camera,capture_time_utc,capture_epoch,pixel_column,easting_utm19,northing_utm19,tide_elevation_navd88"


def ground_csv(path, seed=4, nf=60):
    """The synthetic lines as georectified contour rows (one camera, frames 30 min apart over 6 days)."""
    E, N, Z, F, _ = lines(seed, nf=nf, step=1.0, length=120.0)
    names = sorted(set(F.tolist()))
    t0 = int(datetime(2026, 3, 2, 12, 0, tzinfo=timezone.utc).timestamp())
    epoch = {f: t0 + 86400 * (i // 10) + 1800 * (i % 10) for i, f in enumerate(names)}
    rows = [HEADER]
    col = {}
    for e, n, z, f in zip(E, N, Z, F):
        ep = epoch[f]
        col[f] = col.get(f, -1) + 1
        t = datetime.fromtimestamp(ep, tz=timezone.utc)
        rows.append(f"{ep}.c2.timex,c2,{t:%Y-%m-%dT%H:%M:%S}+00:00,{ep},{col[f]},{e:.3f},{n:.3f},{z:.4f}")
    Path(path).write_text("\n".join(rows) + "\n")
    return E, N, Z, F


def test_dem_from_contours(d):
    print("dem_from_contours.py --loess: the DEM, its info and page unchanged; the loess grids beside them")
    d = Path(d) / "dfc"
    d.mkdir(parents=True, exist_ok=True)
    g = d / "g.csv"
    ground_csv(g)
    rc0, log0 = run([HERE / "dem_from_contours.py", g, d / "plain"], d)
    rc1, log1 = run([HERE / "dem_from_contours.py", g, d / "lo", "--loess"], d)
    (d / "plain.log").write_text(log0)
    (d / "lo.log").write_text(log1)
    same = rc0 == 0 and rc1 == 0 and all((d / f"plain_{k}.asc").read_bytes() == (d / f"lo_{k}.asc").read_bytes()
                                         for k in ("dem", "spread", "count", "sigma"))
    check(same, f"exit {rc0} / {rc1}; _dem, _spread, _count, _sigma the same to the byte with --loess as without")
    i0 = json.loads((d / "plain_info.json").read_text())
    i1 = json.loads((d / "lo_info.json").read_text())
    lo_rec = i1.pop("loess", None)
    check("loess" not in i0 and i0 == i1 and lo_rec and lo_rec["files"]["err"] == "_loess_err.asc"
          and lo_rec["cells"]["filled"] > 0 and lo_rec["settings"]["max_sigma_m"] == 0.25
          and lo_rec["settings"]["max_nmse"] == dl.LOESS_MAX_NMSE,
          f"info.json: without --loess no 'loess' key; with it the same + 'loess' ({lo_rec and lo_rec['cells']['filled']} "
          f"cells, err median {lo_rec and lo_rec['err_median_m']} m; scatter limit --max-spread / 2 = 0.25 m)")
    check(not any((d / f"plain{s}.asc").exists() for s in dl.LOESS_FILES.values())
          and all((d / f"lo{s}.asc").exists() for s in dl.LOESS_FILES.values())
          and "loess DEM" not in log0 and "loess DEM         :" in log1 and "WARNING" not in log1,
          "the loess grids only with --loess; its lines in the log")
    head = [ln.split()[1] for ln in (d / "lo_dem.asc").read_text().splitlines()[:5]]
    check(all([ln.split()[1] for ln in (d / f"lo{s}.asc").read_text().splitlines()[:5]] == head
              for s in dl.LOESS_FILES.values()),
          f"on the DEM's own grid (ncols, nrows, corner, cell: {' '.join(head)})")
    from matplotlib.image import imread
    p0, p1 = imread(str(d / "plain_dem.png")), imread(str(d / "lo_dem.png"))
    check(p0.shape == p1.shape and np.array_equal(p0, p1),
          f"the DEM page is the same to the pixel ({p0.shape[1]} x {p0.shape[0]}): no figure draws the loess")
    dem, z, err = read_asc(d / "lo_dem.asc"), read_asc(d / "lo_loess.asc"), read_asc(d / "lo_loess_err.asc")
    both = np.isfinite(dem) & np.isfinite(z)
    x = (z - dem)[both]
    check(both.sum() > 100 and abs(np.median(x)) < 0.02 and np.isfinite(err).sum() == np.isfinite(z).sum(),
          f"the loess and the binned DEM agree where both are filled ({int(both.sum())} cells, median "
          f"{np.median(x):+.3f} m); an error on every loess cell")
    rcs, logs = run([HERE / "dem_from_contours.py", g, d / "bad", "--loess", "--loess-scales", "20", "--no-plot"], d)
    rcn, logn = run([HERE / "dem_from_contours.py", g, d / "bad", "--loess-scales", "20,5", "--no-plot"], d)
    rco, logo = run([HERE / "dem_from_contours.py", g, d / "s10", "--loess", "--loess-scales", "30,8", "--no-plot"], d)
    i10 = json.loads((d / "s10_info.json").read_text()) if rco == 0 else {}
    check(rcs == 2 and "--loess-scales" in logs and rcn == 2 and "goes with --loess" in logn
          and not (d / "bad_dem.asc").exists() and rco == 0
          and i10["loess"]["settings"]["along_m"] == 30 and i10["loess"]["settings"]["across_m"] == 8,
          "--loess-scales: one value, or without --loess, refused before any work; 30,8 used as given")
    return d


def test_evaluate(d):
    print("dem_loess.py evaluate: the loess and its error against a survey")
    d = Path(d)
    from dem_from_contours import write_ascii_grid
    dem, h = None, {}
    lines_ = (d / "lo_dem.asc").read_text().splitlines()[:6]
    for ln in lines_:
        k, v = ln.split()[:2]
        h[k.lower()] = float(v)
    dem = read_asc(d / "lo_dem.asc")
    cell, e0, n0 = h["cellsize"], h["xllcorner"], h["yllcorner"]
    nrows, ncols = dem.shape
    X, Y = cell_xy(e0, n0, ncols, nrows, cell)
    # the survey: the true beach on the DEM's grid at 0.5 m (a DSM finer than half a cell: its median per cell)
    fine = 0.5
    k = int(cell / fine)
    xs = e0 + (np.arange(ncols * k) + 0.5) * fine - E0
    ys = n0 + (np.arange(nrows * k) + 0.5) * fine - N0
    FX, FY = np.meshgrid(xs, ys)
    write_ascii_grid(str(d / "survey.asc"), truth(FX, FY), e0, n0, fine)
    rc, log = run([HERE / "dem_loess.py", "evaluate", d / "lo", d / "survey.asc", "--json", d / "ev.json"], d)
    ev = json.loads((d / "ev.json").read_text()) if rc == 0 else {}
    ec = ev.get("error_check", {})
    check(rc == 0 and ev["loess_minus_survey"]["n"] > 100 and abs(ev["loess_minus_survey"]["median_m"]) < 0.02
          and ec["median_removed"]["within_2err"] >= 0.85 and ev.get("loess", {}).get("files")
          and "survey within err / 2 err of the loess" in log,
          f"against the true beach (a 0.5 m grid): loess - survey median {ev.get('loess_minus_survey', {}).get('median_m')} m "
          f"on {ev.get('loess_minus_survey', {}).get('n')} cells; within 2 err {ec.get('median_removed', {}).get('within_2err')}, "
          f"NMAD / err {ec.get('ratio_nmad_to_err')}")
    rng = np.random.RandomState(8)
    pe = E0 + rng.uniform(5, 30, 60)
    pn = N0 + rng.uniform(10, 110, 60)
    pe, pn = e0 + (np.floor((pe - e0) / cell) + 0.5) * cell, n0 + (np.floor((pn - n0) / cell) + 0.5) * cell
    pts = d / "points.csv"
    # at cell centres, 0.3 m above the true beach: the loess 0.3 m below them, within 2 err once that is removed
    pts.write_text("E,N,Z\n" + "\n".join(f"{a:.3f},{b:.3f},{truth(a - E0, b - N0) + 0.3:.4f}" for a, b in zip(pe, pn))
                   + "\n")
    rcp, logp = run([HERE / "dem_loess.py", "evaluate", d / "lo", pts, "--json", d / "evp.json"], d)
    evp = json.loads((d / "evp.json").read_text()) if rcp == 0 else {}
    lp, ecp = evp.get("loess_minus_survey", {}), evp.get("error_check", {})
    check(rcp == 0 and "points, median" in logp and lp.get("n", 0) > 40 and abs(lp["median_m"] + 0.3) < 0.02
          and lp["nmad_m"] < 0.03 and abs(ecp["median_offset_m"] + 0.3) < 0.02
          and ecp["as_is"]["within_2err"] < 0.2 and ecp["median_removed"]["within_2err"] >= 0.8,
          f"survey points at cell centres, 0.3 m above the beach, each against the cell it lies in: loess - points "
          f"median {lp.get('median_m')} m, NMAD {lp.get('nmad_m')} m on {lp.get('n')} points; within 2 err "
          f"{ecp.get('as_is', {}).get('within_2err')} as they are, {ecp.get('median_removed', {}).get('within_2err')} "
          f"with the median offset removed")
    # a DEM rebuilt without --loess beside the old loess grids: refused, not compared
    shutil.copy(str(d / "plain_info.json"), str(d / "lo_info.json"))
    rcs, logs = run([HERE / "dem_loess.py", "evaluate", d / "lo", d / "survey.asc"], d)
    check(rcs != 0 and "earlier build" in logs, "loess grids an earlier build left beside the DEM: refused")
    rcn, logn = run([HERE / "dem_loess.py", "evaluate", d / "plain", d / "survey.asc"], d)
    check(rcn != 0 and "--loess" in logn, "no loess grids: refused, saying how to build them")


def test_geotiff(d):
    print("asc_to_geotiff.py: the loess surface is NAVD88, its error grids are not")
    d = Path(d)
    files = [d / f"lo{s}.asc" for s in dl.LOESS_FILES.values()]
    rc, log = run([HERE / "asc_to_geotiff.py"] + files, d)
    from compare_dem_survey import read_tiff
    tags = {s: read_tiff(d / f"lo{s}.tif")[1] for s in dl.LOESS_FILES.values()} if rc == 0 else {}

    def vertical(t):
        k = list(t.get(34735) or [])
        return 4096 in k[4::4]
    check(rc == 0 and vertical(tags["_loess"]) and not any(vertical(tags[s]) for s in ("_loess_err", "_loess_nmse",
                                                                                       "_loess_skill"))
          and "loess elevation, m NAVD88" in str(tags["_loess"].get(270))
          and "prediction error" in str(tags["_loess_err"].get(270)),
          "_loess.tif carries the vertical keys (EPSG:5703) and says what it is; _loess_err/_nmse/_skill do not")
    z = read_asc(d / "lo_loess.asc")
    t = read_tiff(d / "lo_loess.tif")[0].astype(float)
    t[t == -9999] = np.nan
    check(np.allclose(np.where(np.isnan(t), -1.0, t), np.where(np.isnan(z), -1.0, z), atol=1e-4),
          "the GeoTIFF holds the grid's values")


# ---------------------------------------------------------------------
# survey_products.py --dem-loess
# ---------------------------------------------------------------------

def test_survey_products(d):
    print("survey_products.py --dem-loess: off by default, the dem step as it was")
    import survey_products as sp
    plan = {"date": "2025-03-06", "first": "2025-03-03", "last": "2025-03-09", "out": Path(d) / "sp",
            "cams": {"c2": {"station": "CACO04", "eo_file": "x_20250219_EO.yaml"}}}
    base = dict(dem_cell=2.0, dem_min_points=3, dem_max_spread=0.5, dem_max_hs=1.5, dem_max_day_offset=0.15,
                geotiff_epsg=sp.GEOTIFF_EPSG)
    _, cmd0, tif0, page0 = sp.dem_cmds(plan, SimpleNamespace(**base))
    _, cmdf, tiff, pagef = sp.dem_cmds(plan, SimpleNamespace(dem_loess=False, **base))
    stem, cmd1, tif1, page1 = sp.dem_cmds(plan, SimpleNamespace(dem_loess=True, **base))
    want = [f"{stem}_{k}.asc" for k in ("dem", "spread", "count", "sigma")]
    check([str(c) for c in cmd0] == [str(c) for c in cmdf] and "--loess" not in [str(c) for c in cmd0]
          and [str(c) for c in tif0][2:6] == want and str(tif0[6]) == "--epsg" and tif0 == tiff and page0 == pagef,
          "without --dem-loess (or from options that predate it): no --loess, the four grids converted, as before")
    check([str(c) for c in cmd1] == [str(c) for c in cmd0] + ["--loess"]
          and [str(c) for c in tif1][2:10] == want + [f"{stem}_{k}.asc" for k in sp.DEM_LOESS_GRIDS]
          and page1 == page0,
          "with --dem-loess: dem_from_contours.py --loess, and the loess grids converted too; the page as before")
    st = sp.dem_loess_settings()
    check(sp.dem_grids(SimpleNamespace()) == ("dem", "spread", "count", "sigma")
          and st["half_widths_m"] == {"alongshore": 20.0, "cross_shore": 5.0} and "dem_loess.py evaluate" in st["note"]
          and tuple(f"_{k}" for k in sp.DEM_LOESS_GRIDS) == tuple(dl.LOESS_FILES.values()),
          "the grids it names are dem_loess.py's; its stamp record says how to check the error")
    rc, out = run([HERE / "survey_products.py", "--help"], d)
    check(rc == 0 and "--dem-loess" in out, "--dem-loess is an option of survey_products.py")


def test_survey_products_build(d):
    print("survey_products.py --dem-loess: a dem step built with it, then without it")
    import survey_products as sp
    d = Path(d) / "sp_build"
    (d / "tables").mkdir(parents=True)
    (d / "tables" / "survey_dates.csv").write_text(
        "date,camera,enabled,station,first_day,last_day,eo_file,envelope_survey,utc_hours,era,notes\n"
        "2026-03-02,c2,1,CACO05,2026-03-02,2026-03-07,CACO05_c2_20251113_EO-CV.yaml,,0-24,live,a test row\n")
    (d / "tables" / "surveys.csv").write_text("date,name,path,survey_type,survey_date,label,why,label_override_reason\n")
    cfg = sp.read_table(d / "tables" / "survey_dates.csv", sp.CONFIG_FIELDS)
    srv = sp.read_table(d / "tables" / "surveys.csv", sp.SURVEY_FIELDS)
    args = SimpleNamespace(
        config=str(d / "tables" / "survey_dates.csv"), surveys=str(d / "tables" / "surveys.csv"), force_disabled=False,
        window=None, utc_hours=None, calibration=str(HERE / "calibration"), survey_dirs=[str(d / "no_surveys")],
        output_root=str(d / "out"), setup_coef=sp.SETUP_COEF, dem_cell=2.0, dem_min_points=3, dem_max_spread=0.5,
        dem_max_hs=0, dem_max_day_offset=0.15, geotiff_epsg=sp.GEOTIFF_EPSG, force=False, dem_loess=True,
        photo_roots=[str(d / "no_photos")], live_cron=str(sp.LIVE_CRON), pointing_frames=3,
        pointing_fail=sp.POINT_FAIL_DEG, pointing_log_dir=str(d / "no_log"), keep_moved_days=False, no_download=True)
    for k in ("live_contours", "gnssr_spline", "gauge_csv", "setup_fitted_to", "synthetic", "eo", "calibration_variant",
              "replace", "dry_run", "keep_detector_debug", "adcp", "adcp_navd88", "chatham", "wis", "ndbc"):
        setattr(args, k, None)
    plan = sp.build_plan("2026-03-02", cfg, srv, args)
    (plan["out"] / "waterlines").mkdir(parents=True)
    ground_csv(plan["out"] / "waterlines" / "contour_points_ground.csv")
    stem = sp.dem_stem(plan)
    loess_files = [Path(f"{stem}_{k}.{e}") for k in sp.DEM_LOESS_GRIDS for e in ("asc", "tif")]

    def build():
        st = {"steps": {}, "failures": [], "not_applied": []}
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()):
            res = sp.step_dem(plan, args, st)
            st["steps"] = {"dem": res}
            sp.write_provenance(plan, args, st, True)
        stamp = sp.read_stamp(plan["out"], "dem")
        prov = json.loads((plan["out"] / "provenance.json").read_text())
        return res, stamp, prov, (plan["out"] / "README.txt").read_text(), Path(f"{stem}_dem.asc").read_bytes()

    res1, st1, prov1, rd1, dem1 = build()
    check(res1 == "built" and all(f.exists() for f in loess_files) and st1["settings"]["loess"]["half_widths_m"]
          == {"alongshore": 20.0, "cross_shore": 5.0} and "dem_loess.py" in st1["software"]["scripts"]
          and "--loess" in st1["sig"]["cmd"],
          "built with --dem-loess: the loess .asc and .tif beside the DEM; its stamp records the loess and "
          "dem_loess.py's sha256")
    check(prov1["dem"].get("loess") and "loess_err" in prov1["dem"]["files"] and "_loess GeoTIFF" in prov1["dem"]["geotiff"]
          and "loess (--dem-loess): dem_loess.py half-widths 20 m alongshore x 5 m cross-shore" in rd1
          and "_loess/_loess_err/_loess_nmse/_loess_skill .asc and .tif" in rd1
          and "--dem-loess" in prov1["commands"]["survey_products"],
          "provenance.json and the README say so, and the rebuild line carries --dem-loess")
    args.dem_loess = False
    res2, st2, prov2, rd2, dem2 = build()
    check(res2 == "built" and not any(f.exists() for f in loess_files) and "loess" not in st2["settings"]
          and "dem_loess.py" not in st2["software"]["scripts"] and "--loess" not in st2["sig"]["cmd"]
          and dem1 == dem2,
          "rebuilt without it: the loess files of the earlier build removed, nothing of it in the stamp; the DEM "
          "the same to the byte")
    check("loess" not in prov2["dem"] and set(prov2["dem"]["files"]) == {"dem", "spread", "count"}
          and "loess" not in rd2.replace(str(d), "").replace(str(HERE), "")
          and "--dem-loess" not in prov2["commands"]["survey_products"],
          "provenance.json, the README and the rebuild line as without the option")


def test_py38():
    print("Python 3.8 syntax")
    import ast
    for f in ("dem_loess.py", "test_dem_loess.py", "dem_from_contours.py", "survey_products.py", "asc_to_geotiff.py"):
        src = (HERE / f).read_text(encoding="utf-8")
        try:
            ast.parse(src, feature_version=(3, 8))
        except TypeError:                        # Python 3.7 has no feature_version
            ast.parse(src)
        check(True, f"{f} parses as Python 3.8")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="keep the outputs, in a new folder inside this one")
    args = ap.parse_args()
    if args.keep:
        Path(args.keep).mkdir(parents=True, exist_ok=True)
    d = os.path.realpath(tempfile.mkdtemp(prefix="dem_loess_", dir=args.keep))
    print(f"test folder: {d}")
    try:
        test_exact_quadratic()
        test_recovery()
        test_calibration()
        test_error_grows_into_gap()
        test_error_grows_beyond_data()
        test_nmse_per_frame()
        test_response()
        dd = test_dem_from_contours(d)
        test_evaluate(dd)
        test_geotiff(dd)
        test_survey_products(d)
        test_survey_products_build(d)
        test_py38()
    finally:
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED")
        return 1
    print("all loess DEM tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

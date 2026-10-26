#!/usr/bin/env python3
"""
Tests For The Validation Statistics
=====================================
Plain python3 (no pytest needed on the station):

    python3 test_validation_stats.py [--keep DIR]

WHAT.
  * survey_compare.py's intervals: the block bootstrap's statistics are
    those of each resample spelled out; values that err together are worth
    one value (the effective sample size falls to the number of blocks)
    and their intervals cover the truth where the plain n-value interval
    does not; the residuals' correlogram finds how far a planted error
    runs along the shore; the direction along the shore is perpendicular
    to the survey's slope where compared;
  * the waterline strata: a planted water-level timing error (lines high
    on the flood, low on the ebb) and a planted wave dependence come back
    with intervals that exclude 0, none planted gives intervals over 0,
    and a range dependence shows by range band; a pointing drift while the
    daylight frames move from ebb to flood (no timing error) leaves the
    within-camera-day flood - ebb at 0, where pooled medians read it as
    timing; slack frames, Hs thirds and camera-days on the local day;
  * a whole comparison: the INTERVAL line, the UNCERTAINTY table and the
    strata are in the text and the json, the figure and the CSV still
    drawn and written as before; intervals, or the frames' rates and Hs
    for the strata, that cannot be computed never stop it;
  * the labels' edge cases (survey_compare.audit_label, suggested_label
    and survey_products.honest_label): calibration targets a calibration
    names are PARTLY-CIRCULAR, as its date says, also end to end through
    compare(); a pointing 'fitted to the Jan lidar' (no file named) is
    unresolved, never 'another survey', while a survey named by its stem
    still is a file; a sentence in provenance.json naming the survey is a
    rule, a file a step used a chain fact; and the rules that were right
    stay right (same stem another format, a short stem, the label only gets
    worse, the worst finding wins, CROSS-VALIDATED ranks between, a label
    that is no label is refused, a comparison's provenance entries are
    ignored);
  * dem_from_contours.day_block_sigma (--day-block-se): the independent SE
    again with one crossing per day, k values' worth when each day's
    crossings sit on one side, none on one day; no other grid changes;
  * dem_change.py: the standard errors read from each DEM's _sigma.asc
    (and the formula where there is none), a file of the same definition
    leaving the _sig and _diff grids byte for byte those of the formula on
    noisy grids at the level of detection, and --day-block-lod printed
    beside the map's numbers without changing a grid, also with --rebuild.

HOW. Synthetic grids, points and ground files in a temporary folder
(--keep DIR to look at them); the scripts run as the station runs them.

STATION FILES. None are read: the station's survey folders are pointed at
folders of the test's own, and the lens files a label test needs are
copied next to its EO files from this repository's calibration/.
"""

import argparse
import csv
import re
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import contextlib
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import survey_compare as sc                                   # noqa: E402

FAILED = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILED.append(msg)


def quiet(fn, *a, **k):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = fn(*a, **k)
    return out, buf.getvalue()


def run(args, d):
    r = subprocess.run([sys.executable] + [str(a) for a in args], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
    return r.returncode, r.stdout


def read_asc(path):
    a = np.loadtxt(str(path), skiprows=6, ndmin=2)
    return np.where(a == -9999, np.nan, a)


# ---------------------------------------------------------------------
# The block bootstrap
# ---------------------------------------------------------------------

def test_boot_stats_exact():
    print("survey_compare.boot_stats: the statistics of each resample spelled out")
    rng = np.random.RandomState(4)
    d = rng.normal(0.1, 0.3, 57)
    blk = rng.randint(0, 9, 57)
    draws = sc.block_draws(9, 200, 3)
    fast = sc.boot_stats(d, blk, draws)
    worst = 0.0
    for r in range(200):
        full = np.concatenate([np.tile(d[blk == k], int(c)) for k, c in enumerate(draws[r])])
        want = (np.median(full), full.mean(), np.sqrt(np.mean(full ** 2)))
        worst = max(worst, float(np.max(np.abs(np.array(want) - fast[r]))))
    check(worst < 1e-12, f"median (two middle values averaged on an even count), mean and RMSE equal np.median etc. "
                         f"of the resample written out, on 200 resamples (largest difference {worst:.1e})")
    only0 = sc.boot_stats(d[blk == 0], np.zeros(int((blk == 0).sum()), int), draws[:, :1])
    gone = draws[:, 0] == 0
    check(gone.any() and np.isnan(only0[gone]).all() and np.isfinite(only0[~gone]).all(),
          "a resample that draws none of the values gives NaN, not a number")
    a = sc.block_draws(7, 50, 1)
    check(np.array_equal(a, sc.block_draws(7, 50, 1)) and (a.sum(axis=1) == 7).all(),
          "the draws are the same on every run (fixed seed), 7 blocks drawn per resample")


def test_block_intervals():
    print("survey_compare.block_intervals: the effective sample size and the interval")
    rng = np.random.RandomState(7)
    iid = rng.normal(0, 0.2, 300)
    r = sc.block_intervals(iid, rng.randint(0, 30, 300), "cells", "block(s)")
    check(r["n"] == 300 and 0.6 * 300 < r["n_effective"] <= 300 and r["median_ci"][0] < 0 < r["median_ci"][1],
          f"independent values in 30 blocks: n_eff {r['n_effective']:.0f} of 300 (about n), interval over the truth")
    block = np.repeat(rng.normal(0, 0.2, 20), 15) + rng.normal(0, 0.01, 300)
    lab = np.repeat(np.arange(20), 15)
    r2 = sc.block_intervals(block, lab, "cells", "block(s)")
    naive = 1.96 * block.std(ddof=1) / np.sqrt(300)
    check(r2["n_effective"] < 40 and (r2["mean_ci"][1] - r2["mean_ci"][0]) / 2 > 3 * naive,
          f"20 blocks of 15 equal errors: n_eff {r2['n_effective']:.0f} (about the 20 blocks, not 300); the interval's "
          f"half-width {(r2['mean_ci'][1] - r2['mean_ci'][0]) / 2:.3f} m against {naive:.3f} m for 300 independent "
          f"values")
    r3 = sc.block_intervals(block[:60], lab[:60], "cells", "stretch(es) of 50 m along the shore")
    check("median_ci" not in r3 and r3["blocks"] == 4 and "only 4" in r3["note"],
          f"4 blocks: no interval, said ({r3['note']})")
    r4 = sc.block_intervals(block[:120], lab[:120], "cells", "block(s)")
    check(r4.get("rough") is True and "median_ci" in r4, "8 blocks: an interval, marked rough")
    r5 = sc.block_intervals(np.full(10, 0.2), np.arange(10), "points", "block(s)")
    check(r5["n_effective"] == 10, f"values that do not scatter: n_eff {r5['n_effective']:.0f} of 10, not rounding noise")
    # coverage: the truth (0) inside the 95% interval of the mean
    cover_b = cover_n = 0
    sims = 150
    for s in range(sims):
        g = np.random.RandomState(100 + s)
        v = np.repeat(g.normal(0, 0.2, 12), 10) + g.normal(0, 0.05, 120)
        rr = sc.block_intervals(v, np.repeat(np.arange(12), 10), "cells", "block(s)", n_boot=400, seed=s)
        cover_b += rr["mean_ci"][0] <= 0 <= rr["mean_ci"][1]
        h = 1.96 * v.std(ddof=1) / np.sqrt(len(v))
        cover_n += v.mean() - h <= 0 <= v.mean() + h
    check(cover_b / sims >= 0.8 and cover_n / sims <= 0.6,
          f"12 blocks of shared errors, {sims} simulations: the block interval holds the truth {100 * cover_b / sims:.0f}% "
          f"of the time, the interval of 120 independent values {100 * cover_n / sims:.0f}%")


def test_contrast():
    print("survey_compare.block_contrast: the difference of two groups' medians, resampled by blocks")
    rng = np.random.RandomState(2)
    day = np.repeat(np.arange(12), 10)
    a = np.tile(np.r_[np.ones(5, bool), np.zeros(5, bool)], 12)
    d = np.repeat(rng.normal(0, 0.05, 12), 10) + np.where(a, 0.1, -0.1) + rng.normal(0, 0.02, 120)
    c = sc.block_contrast(d, day, a, ~a)
    check(abs(c["difference_m"] - 0.2) < 0.03 and c["ci"][0] > 0 and c["blocks"] == 12,
          f"+0.2 planted between the groups: {c['difference_m']:+.3f} m [{c['ci'][0]:+.3f}, {c['ci'][1]:+.3f}]")
    c0 = sc.block_contrast(np.repeat(rng.normal(0, 0.05, 12), 10) + rng.normal(0, 0.02, 120), day, a, ~a)
    check(c0["ci"][0] < 0 < c0["ci"][1], f"none planted: the interval holds 0 ({c0['ci'][0]:+.3f}, {c0['ci'][1]:+.3f})")
    check(sc.block_contrast(d, day, a, np.zeros(120, bool)) is None, "an empty group: None")
    big = np.repeat(rng.normal(0, 0.3, 12), 10) + np.where(a, 0.1, -0.1) + rng.normal(0, 0.02, 120)
    cb = sc.block_contrast(big, day, a, ~a)
    check((cb["ci"][1] - cb["ci"][0]) / 2 < 0.1,
          f"days with 0.3 m errors, both groups on every day: the two groups drawn from the same days (half-width "
          f"{(cb['ci'][1] - cb['ci'][0]) / 2:.3f} m, against ~0.15 m with independent draws for the two)")
    w = sc.within_block_contrast(big, day, a, ~a)
    check(abs(w["difference_m"] - 0.2) < 0.03 and w["ci"][0] > 0.15 and w["blocks"] == 12
          and w["blocks_one_group"] == 0 and len(w["per_block"]) == 12,
          f"within blocks: {w['difference_m']:+.3f} m [{w['ci'][0]:+.3f}, {w['ci'][1]:+.3f}] over 12 blocks, the day "
          f"errors gone")
    one = (day < 3) | a                                       # days 0-2 hold only group a
    w2 = sc.within_block_contrast(big, day, a, ~a & ~one)
    check(w2["blocks"] == 9 and w2["blocks_one_group"] == 3, "blocks holding one group only: left out, counted")
    w3 = sc.within_block_contrast(big, day, a & (day < 2), ~a & (day >= 2))
    check(w3["blocks"] == 0 and w3["difference_m"] is None and "ci" not in w3 and w3["note"],
          f"no block holds both groups: no difference, said ({w3['note']})")
    w4 = sc.within_block_contrast(big, day, a & (day < 4), ~a & (day < 4))
    check(w4["blocks"] == 4 and "ci" not in w4 and "only 4" in w4["note"], f"4 blocks: no interval ({w4['note']})")


def test_correlogram():
    print("survey_compare.alongshore_correlogram: how far the errors run along the shore")
    rng = np.random.RandomState(5)
    s = np.arange(0.0, 600.0, 1.0)
    d = np.repeat(rng.normal(0, 0.2, 10), 60) + rng.normal(0, 0.02, 600)      # 60 m stretches
    gram, L = sc.alongshore_correlogram(d, s)
    check(L is not None and 25 < L < 120 and gram[0]["correlation"] > 0.8,
          f"errors shared over 60 m stretches: correlation {gram[0]['correlation']:+.2f} at 0-5 m, falls to 0.1 at "
          f"~{L:.0f} m")
    gram2, L2 = sc.alongshore_correlogram(rng.normal(0, 0.2, 600), s)
    check(L2 == 0.0, "independent errors: correlated over 0 m")
    check(sc.alongshore_correlogram(np.zeros(50), s[:50]) == ([], None), "no scatter: no correlogram")
    big = np.repeat(rng.normal(0, 0.2, 50), 100)
    g3, L3 = sc.alongshore_correlogram(big + rng.normal(0, 0.02, 5000), np.arange(5000.0) * 0.12)
    check(L3 is not None and 6 < L3 < 24, f"5000 values (a fixed subset of 2000 used): ~{L3:.0f} m for 12 m stretches")


def test_dem_intervals():
    print("survey_compare.dem_intervals / interval_headline: the stretch length and the INTERVAL line")
    rng = np.random.RandomState(6)

    def iv(d, N):
        n = len(d)
        return sc.dem_intervals(np.asarray(d, float), np.zeros(n), np.asarray(N, float), np.array(["c1"] * n),
                                np.full(n, 120.0), (0.0, 1.0), "north", "cells")[0]
    N = np.arange(0.0, 2000.0, 2.0)
    r = iv(np.repeat(rng.normal(0, 0.2, 10), 100) + rng.normal(0, 0.02, 1000), N)     # 200 m stretches
    check(r["correlation_length_m"] is not None and r["correlation_length_m"] > 60
          and r["block_m"] == r["correlation_length_m"] and r["rows"][0]["blocks"] == int(np.ptp(N) // r["block_m"]) + 1,
          f"errors shared over 200 m: the stretches resampled are as long as the correlation, "
          f"{r['block_m']:.0f} m (not the 50 m floor), {r['rows'][0]['blocks']} of them")
    sparse = iv(rng.normal(0.1, 0.05, 10), np.arange(10) * 40.0)
    hs = sc.interval_headline(sparse, "points")
    check(not sparse["correlation_measured"] and sparse["rows"][0].get("median_ci") and "too narrow" not in hs
          and "too few pairs" in hs,
          "10 points 40 m apart: the correlation cannot be measured, and the INTERVAL line says so (not 'too narrow')")
    ind = iv(rng.normal(0, 0.2, 300), np.linspace(0.0, 190.0, 300))
    hi = sc.interval_headline(ind, "cells")
    check(ind["correlation_length_m"] == 0.0 and "median_ci" not in ind["rows"][0] and "look independent" in hi
          and "effective n about the 300 cells" in hi,
          "independent errors over 190 m (4 stretches, no interval): the INTERVAL line says they look independent")


def test_alongshore_axis():
    print("survey_compare.alongshore_axis: along the shore = across the survey's slope where compared")
    cell, xll, ytop = 2.0, 420000.0, 4638600.0
    rr, cc = np.mgrid[0:150, 0:30]
    E, N = (xll + (cc + 0.5) * cell).ravel(), (ytop - (rr + 0.5) * cell).ravel()
    ue, un, how = sc.alongshore_axis(E, N, 3.0 - 0.05 * (E - xll))
    check(abs(un) > 0.999 and "1:20" in how and "along 0 deg" in how and "survey" in how,
          f"falling east: along the shore is north ({how})")
    a = np.radians(60.0)                                      # falling toward the bearing 60 deg
    ue, un, how = sc.alongshore_axis(E, N, 3.0 - 0.08 * ((E - xll) * np.sin(a) + (N - ytop) * np.cos(a)))
    check("along 150 deg" in how, f"falling toward 60 deg: along 150 deg ({how})")
    ue, un, how = sc.alongshore_axis(E, N, np.full(E.shape, 1.0))
    check(abs(un) > 0.99 and "long axis" in how and "flat" not in how,
          f"no slope, 300 m N by 60 m E: the long axis, without calling the beach flat ({how})")


# ---------------------------------------------------------------------
# Waterline strata
# ---------------------------------------------------------------------

def strata_frames(days=10, stage=0.0, hs_slope=0.0, seed=1):
    """Two cameras, 12 frames a day each, alternately on the flood (+0.5 m/h) and the ebb (-0.5 m/h);
    each frame's value: its day's error + stage x sign(dz/dt) + hs_slope x (Hs - 1) + noise.
    -> (d, frame, cam, range, rate_of, hs_of, block_of)"""
    rng = np.random.RandomState(seed)
    d, fr, cam, rng_m = [], [], [], []
    rate_of, hs_of, block_of = {}, {}, {}
    for day in range(days):
        hs = 0.5 + 0.1 * day
        for c in ("c1", "c2"):
            day_err = rng.normal(0, 0.03)
            for k in range(12):
                f = f"d{day}{c}f{k}"
                rate = 0.5 if k % 2 else -0.5
                rate_of[f], hs_of[f], block_of[f] = rate, hs, f"{c} 2026-10-{day + 1:02d}"
                d.append(day_err + stage * np.sign(rate) + hs_slope * (hs - 1.0) + rng.normal(0, 0.02))
                fr.append(f)
                cam.append(c)
                rng_m.append(120.0)
    return np.array(d), np.array(fr), np.array(cam), np.array(rng_m), rate_of, hs_of, block_of


def drift_frames(days=10, drift=0.02, seed=4):
    """Two cameras, 12 frames a day each; on day i the first 1 + i frames on the flood (+0.5 m/h), the
    rest on the ebb (-0.5 m/h); each frame's value: drift x i (a pointing drift, the same for both
    stages) + a small camera-day error + noise. No timing error. -> as strata_frames"""
    rng = np.random.RandomState(seed)
    d, fr, cam, rng_m = [], [], [], []
    rate_of, hs_of, block_of = {}, {}, {}
    for day in range(days):
        for c in ("c1", "c2"):
            day_err = drift * day + rng.normal(0, 0.005)
            for k in range(12):
                f = f"d{day}{c}f{k}"
                rate_of[f], hs_of[f], block_of[f] = (0.5 if k <= day else -0.5), 0.5 + 0.1 * day, f"{c} 2026-10-{day + 1:02d}"
                d.append(day_err + rng.normal(0, 0.01))
                fr.append(f)
                cam.append(c)
                rng_m.append(120.0)
    return np.array(d), np.array(fr), np.array(cam), np.array(rng_m), rate_of, hs_of, block_of


def test_strata():
    print("survey_compare.waterline_strata: tide stage, waves and range, with intervals")
    import dem_from_contours as dfc
    check(sc.STAGE_RATE == dfc.DIAG_STAGE_RATE,
          f"flood / ebb split at the DEM log's own rate ({sc.STAGE_RATE} m/h = dem_from_contours.DIAG_STAGE_RATE)")
    d, fr, cam, rg, ro, ho, bo = strata_frames(stage=0.1)
    rec, lines = sc.waterline_strata(d, fr, cam, rg, ro, ho, bo)
    rows = {r["group"]: r for r in rec["rows"]}
    fe = rec["flood_minus_ebb"]
    check(abs(rows["flood"]["median"] - 0.1) < 0.03 and abs(rows["ebb"]["median"] + 0.1) < 0.03
          and abs(fe["difference_m"] - 0.2) < 0.03 and fe["ci"][0] > 0.1 and rows["flood"]["camera_days"] == 20,
          f"lines 0.1 m high on the flood and low on the ebb: flood {rows['flood']['median']:+.3f}, ebb "
          f"{rows['ebb']['median']:+.3f}, flood - ebb {fe['difference_m']:+.3f} [{fe['ci'][0]:+.3f}, {fe['ci'][1]:+.3f}] "
          f"m over 20 camera-days")
    bc = rec["by_camera"]
    check(all(abs(bc[c]["flood_minus_ebb"]["difference_m"] - 0.2) < 0.04 and bc[c]["flood_minus_ebb"]["ci"][0] > 0
              for c in ("c1", "c2")), "the same camera by camera")
    txt = "\n".join(lines)
    check("WATERLINES BY TIDE STAGE, WAVES AND RANGE" in txt and "flood - ebb within camera-days: +0." in txt
          and "slack" not in rows, "printed; no slack frame, no slack row")
    d0, fr0, cam0, rg0, ro0, ho0, bo0 = strata_frames(stage=0.0, seed=2)
    rec0, _ = sc.waterline_strata(d0, fr0, cam0, rg0, ro0, ho0, bo0)
    held, sims = 0, 100
    for k in range(sims):
        dk, frk, _, _, rok, _, bok = strata_frames(stage=0.0, seed=100 + k)
        fl = np.array([rok[f] > 0 for f in frk])
        ck = sc.within_block_contrast(dk, np.array([bok[f] for f in frk]), fl, ~fl, n_boot=400, seed=k)["ci"]
        held += ck[0] <= 0 <= ck[1]
    check(held >= 0.88 * sims, f"none planted, {sims} simulations: the flood - ebb interval holds 0 {held}% of the time")
    dh, frh, camh, rgh, roh, hoh, boh = strata_frames(hs_slope=0.2, seed=3)
    rech, _ = sc.waterline_strata(dh, frh, camh, rgh, roh, hoh, boh)
    hl = rech["high_minus_low_hs"]
    check(rech["hs_tercile_cuts_m"] is not None and hl["difference_m"] > 0.1 and hl["ci"][0] > 0,
          f"0.2 m per m of Hs planted: highest - lowest third {hl['difference_m']:+.3f} m [{hl['ci'][0]:+.3f}, "
          f"{hl['ci'][1]:+.3f}], cuts {rech['hs_tercile_cuts_m']}")
    # range: every frame a line of points from 50 to 350 m, the error growing 0.001 m per m of range
    pts_d, pts_f, pts_c, pts_r = [], [], [], []
    for f, c in zip(fr0, cam0):
        for r_ in np.arange(50.0, 350.0, 10.0):
            pts_d.append(0.001 * (r_ - 200.0)); pts_f.append(f); pts_c.append(c); pts_r.append(r_)
    recr, _ = sc.waterline_strata(np.array(pts_d), np.array(pts_f), np.array(pts_c), np.array(pts_r), ro0, ho0, bo0)
    rr_ = [r for r in recr["rows"] if r["group_type"] == "range"]
    meds = [r["median"] for r in rr_]
    check(len(rr_) >= 5 and all(np.diff(meds) > 0) and rr_[0]["frames"] == len(set(fr0.tolist()))
          and recr["rows"][0]["frames"] == len(set(fr0.tolist())),
          f"a range dependence shows band by band ({', '.join(r['group'] + ' ' + format(r['median'], '+.3f') for r in rr_)}); "
          f"one value per frame and band, one per frame in 'all frames'")
    ro_nan = {k: np.nan for k in ro0}
    recn, linesn = sc.waterline_strata(d0, fr0, cam0, rg0, ro_nan, ho0, bo0)
    check(recn["flood_minus_ebb"] is None and recn["frames_without_rate"] == len(ro0)
          and "without a water-level rate" in "\n".join(linesn),
          "no capture epoch (no rates): no stage rows, said")
    # slack: |dz/dt| within STAGE_RATE
    ro_s = dict(ro0)
    slow = [f for f in sorted(ro_s)[:30]]
    for i, f in enumerate(slow):
        ro_s[f] = 0.05 if i % 2 else -0.05
    recs, _ = sc.waterline_strata(d0, fr0, cam0, rg0, ro_s, ho0, bo0)
    rws = {r["group"]: r for r in recs["rows"]}
    check(rws.get("slack", {}).get("frames") == 30 and rws["flood"]["frames"] + rws["ebb"]["frames"] == len(ro0) - 30,
          "frames at +/-0.05 m/h are slack, not flood or ebb")
    # Hs in thirds of the frames
    fh = np.array([ho0[f] for f in sorted(set(fr0.tolist()))])
    hsr = [r for r in rec0["rows"] if r["group_type"] == "hs"]
    check(np.allclose(rec0["hs_tercile_cuts_m"], np.percentile(fh, [100.0 / 3, 200.0 / 3]))
          and len(hsr) == 3 and all(0.25 * len(fh) <= r["frames"] <= 0.42 * len(fh) for r in hsr),
          f"Hs cut at the frames' thirds ({', '.join(str(r['frames']) for r in hsr)} of {len(fh)} frames)")
    # no timing error, but the days drift (pointing) while the daylight frames move from the ebb to the flood
    dd, frd, camd, rgd, rod, hod, bod = drift_frames()
    recd, linesd = sc.waterline_strata(dd, frd, camd, rgd, rod, hod, bod)
    rwd = {r["group"]: r for r in recd["rows"]}
    fed = recd["flood_minus_ebb"]
    fm = np.array([rod[f] > 0 for f in frd])
    pooled = sc.block_contrast(dd, np.array([bod[f] for f in frd]), fm, ~fm)
    check(rwd["flood"]["median"] - rwd["ebb"]["median"] > 0.03 and pooled["ci"][0] > 0,
          f"a 0.02 m/day drift, the frames moving from ebb to flood: the pooled flood - ebb reads "
          f"{pooled['difference_m']:+.3f} m [{pooled['ci'][0]:+.3f}, {pooled['ci'][1]:+.3f}] (the days' error)")
    check(abs(fed["difference_m"]) < 0.015 and fed["ci"][0] < 0 < fed["ci"][1] and fed["blocks"] == 20
          and all(v["flood_minus_ebb"]["ci"][0] < 0 < v["flood_minus_ebb"]["ci"][1] for v in recd["by_camera"].values()),
          f"... within camera-days it holds 0: {fed['difference_m']:+.3f} m [{fed['ci'][0]:+.3f}, {fed['ci'][1]:+.3f}] "
          f"over 20 camera-days, and camera by camera")
    txd = "\n".join(linesd)
    check("flood - ebb within camera-days" in txd and "pool different camera-days" in txd
          and "Between days by construction" in txd,
          "the text says the stage rows pool days, the contrast is within camera-days, the Hs one between days")
    # the camera-day is the local day: frames either side of local midnight (04:00 UTC in EDT)
    t1 = datetime(2026, 10, 2, 3, 30, tzinfo=timezone.utc).timestamp()
    cont = {"frame": np.array(["fa", "fa", "fb", "fc"]), "epoch": np.array([t1, t1, t1 + 3600, np.nan]),
            "Zt": np.array([1.0, 1.0, 1.5, 1.2]), "hs": np.array([0.8, 0.8, 0.9, 0.7]),
            "cam": np.array(["c1", "c1", "c1", "c2"]), "day": np.array(["2026-10-02"] * 4)}
    rate_of, hs_of, block_of = sc.frame_properties(cont)
    check(block_of == {"fa": "c1 2026-10-01", "fb": "c1 2026-10-02", "fc": "c2 2026-10-02"}
          and abs(rate_of["fa"] - 0.5) < 1e-9 and np.isnan(rate_of["fc"]) and hs_of["fb"] == 0.9,
          f"23:30 and 00:30 local fall on two camera-days, a frame without an epoch on its UTC day ({block_of})")


# ---------------------------------------------------------------------
# A whole comparison
# ---------------------------------------------------------------------

X0, Y0 = 420060.0, 4638000.0


def beach(E, N):
    return 3.0 - 0.05 * (np.asarray(E) - X0) + 0.0 * np.asarray(N)


def write_asc(path, grid_ns, xll, yll, cell):
    with open(path, "w") as f:
        f.write(f"ncols {grid_ns.shape[1]}\nnrows {grid_ns.shape[0]}\nxllcorner {xll}\n"
                f"yllcorner {yll}\ncellsize {cell}\nNODATA_value -9999.0\n")
        for r in grid_ns:
            f.write(" ".join("%.4f" % (v if np.isfinite(v) else -9999.0) for v in r) + "\n")


def test_compare(d):
    print("survey_compare.compare: intervals and strata in the text and json; figure and CSV as before")
    from asc_to_geotiff import write_geotiff
    d = Path(d) / "compare"
    d.mkdir(parents=True)
    cell, ncols, nrows = 2.0, 30, 300                          # 60 m across, 600 m along the shore (north)
    xs = X0 + (np.arange(ncols) + 0.5) * cell
    ys = Y0 + nrows * cell - (np.arange(nrows) + 0.5) * cell
    E, N = np.meshgrid(xs, ys)
    rng = np.random.RandomState(9)
    offs = rng.normal(0.0, 0.15, 10)                           # an error shared by each 60 m stretch
    dem = beach(E, N) + 0.1 + offs[np.minimum(((N - Y0) // 60).astype(int), 9)]
    write_asc(d / "long_dem.asc", dem, X0, Y0, cell)
    sx0, sytop, scl = X0 - 10.3, Y0 + nrows * cell + 10.7, 0.5
    nc, nr = int((ncols * cell + 20) / scl), int((nrows * cell + 20) / scl)
    SE, SN = np.meshgrid(sx0 + (np.arange(nc) + 0.5) * scl, sytop - (np.arange(nr) + 0.5) * scl)
    write_geotiff(d / "long_lidar.tif", beach(SE, SN), sx0, sytop - nr * scl, scl, 6348, "synthetic lidar")
    # waterlines: 10 days, 2 cameras, 8 frames a day, 4 on the flood 0.1 m high, 4 on the ebb 0.1 m low
    cp = d / "contours.csv"
    t0 = int(datetime(2026, 10, 1, 13, 0, tzinfo=timezone.utc).timestamp())
    with open(cp, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["source_file", "camera", "capture_time_utc", "capture_epoch", "tide_elevation_navd88",
                    "offshore_hs_m", "offshore_tp_s", "setup_correction_m", "beach_elevation_navd88",
                    "easting_utm19", "northing_utm19"])
        for day in range(10):
            for k in range(8):
                ep = t0 + day * 86400 + k * 1800
                tide = 1.0 + 0.25 * k if k < 4 else 1.75 - 0.25 * (k - 4)   # rising 0.5 m/h for 2 h, then falling
                err = 0.1 if k < 4 else -0.1
                for cam in ("c1", "c2"):
                    n0 = Y0 + (20 if cam == "c1" else 320)
                    e = X0 + (3.0 - tide) / 0.05
                    for n in np.arange(n0, n0 + 250, 5.0):
                        w.writerow([f"{ep}.{cam}.timex", cam, datetime.fromtimestamp(ep, tz=timezone.utc)
                                    .strftime("%Y-%m-%dT%H:%M:%S+00:00"), ep, f"{tide:.4f}", "0.8", "8.0", "0.0",
                                    f"{tide + err:.4f}", f"{e:.3f}", f"{n:.3f}"])
    out = d / "out"
    h, log = quiet(sc.compare, str(d / "long_dem.asc"), str(d / "long_lidar.tif"), "dsm", "INDEPENDENT",
                   "synthetic test", name="long", output_dir=str(out), survey_date="2026-10-05",
                   photo_dates=("2026-10-01", "2026-10-10"), contours=str(cp))
    iv = h["intervals"]
    r0 = iv["rows"][0]
    txt = (out / "long_comparison.txt").read_text()
    check(r0["n"] > 5000 and r0.get("median_ci") and r0["n_effective"] < r0["n"] / 50 and iv["block_m"] >= 50
          and r0["blocks"] >= 8,
          f"{r0['n']} cells, an error shared by 60 m stretches: block {iv['block_m']:.0f} m, {r0['blocks']} blocks, "
          f"n_eff {r0['n_effective']:.0f}, median {r0['median']:+.3f} {sc.interval_text(r0)}")
    check(r0["median_ci"][0] <= r0["median"] <= r0["median_ci"][1] and "INTERVAL   : 95%: median [" in txt
          and "UNCERTAINTY OF THE STATISTICS" in txt and txt.index("HEADLINE") < txt.index("INTERVAL   :")
          < txt.index("checks     :"), "the INTERVAL line right under the HEADLINE, the UNCERTAINTY table below")
    st = h["waterlines"]["strata"]
    fe = st["flood_minus_ebb"]
    check(abs(fe["difference_m"] - 0.2) < 0.02 and fe["ci"][0] > 0.15 and fe["blocks"] == 20
          and "WATERLINES BY TIDE STAGE" in txt and "flood - ebb within camera-days" in txt,
          f"lines 0.1 m high on the flood, low on the ebb, against the lidar: flood - ebb {fe['difference_m']:+.3f} m "
          f"[{fe['ci'][0]:+.3f}, {fe['ci'][1]:+.3f}] over {fe['blocks']} camera-days")
    js = json.loads((out / "long_comparison.json").read_text())
    groups = {r["group_type"] for r in csv.DictReader((out / "long_comparison.csv").read_text().splitlines())}
    check(js["intervals"]["rows"][0]["n"] == r0["n"] and "strata" in js["waterlines"]
          and groups <= {"overall", "elevation_band", "camera", "distance", "survey_cover", "dem_spread", "dem_frames",
                         "waterline_all", "waterline_camera", "waterline_day", "waterline_band"}
          and (out / "long_comparison.png").exists(),
          "json holds them; the CSV keeps its groups (no intervals in it); the figure is drawn")
    h2, _ = quiet(sc.compare, str(d / "long_dem.asc"), str(d / "long_lidar.tif"), "dsm", "INDEPENDENT", "test",
                  name="long_b100", output_dir=str(out), survey_date="2026-10-05", plot=False, block_m=100.0)
    check(h2["intervals"]["block_m"] >= 100 and h2["median"] == h["median"],
          f"--block-m 100: stretches of {h2['intervals']['block_m']:.0f} m; the headline itself unchanged")
    # report only: intervals that cannot be computed never stop a comparison
    keep = sc.dem_intervals

    def broken(*a, **k):
        raise ValueError("synthetic failure")
    sc.dem_intervals = broken
    try:
        h3, _ = quiet(sc.compare, str(d / "long_dem.asc"), str(d / "long_lidar.tif"), "dsm", "INDEPENDENT", "test",
                      name="long_b100", output_dir=str(out / "broken"), survey_date="2026-10-05", plot=False)
    finally:
        sc.dem_intervals = keep
    t3 = (out / "broken" / "long_b100_comparison.txt").read_text()
    check(h3["intervals"] is None and h3["median"] == h["median"] and "INTERVAL   : not computed" in t3
          and (out / "broken" / "long_b100_comparison.csv").read_bytes() == (out / "long_b100_comparison.csv").read_bytes(),
          "intervals that fail: 'not computed' in the text, None in the json, the headline and the CSV as without them")
    # ... and frame properties that cannot be read (e.g. no pandas for the local day) drop the strata only
    keep_fp = sc.frame_properties

    def no_pandas(*a, **k):
        raise ImportError("No module named 'pandas'")
    sc.frame_properties = no_pandas
    try:
        h4, log4 = quiet(sc.compare, str(d / "long_dem.asc"), str(d / "long_lidar.tif"), "dsm", "INDEPENDENT",
                         "synthetic test", name="long", output_dir=str(out / "nofp"), survey_date="2026-10-05",
                         photo_dates=("2026-10-01", "2026-10-10"), contours=str(cp), plot=False)
    finally:
        sc.frame_properties = keep_fp
    t4 = (out / "nofp" / "long_comparison.txt").read_text()
    check(h4["waterlines"]["strata"] is None and h4["waterlines"]["n_independent"] == h["waterlines"]["n_independent"]
          and h4["median"] == h["median"] and "waterline strata  : not computed" in log4
          and "WATERLINES BY TIDE STAGE" not in t4
          and (out / "nofp" / "long_comparison.csv").read_bytes() == (out / "long_comparison.csv").read_bytes(),
          "frame properties that fail: no strata (said), the waterlines, the headline and the CSV as with them")


# ---------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------

def test_labels(d):
    print("labels: edge cases of survey_compare.audit_label (and survey_products.honest_label)")
    import survey_products as sp
    from types import SimpleNamespace
    d = Path(d) / "labels"
    cal = d / "cal"
    cal.mkdir(parents=True)
    for c in ("c1", "c2"):
        shutil.copy(str(HERE / "calibration" / f"CACO05_{c}_20240801_IO.yaml"), str(cal / f"CACO05_{c}_20240801_IO.yaml"))
    dem = d / "x_dem.asc"
    write_asc(dem, np.full((5, 5), 1.0), X0, Y0, 2.0)
    lid = d / "2025005FA_test_Lidar_DSM.tif"
    lid.write_bytes(b"not read")
    base = "x: 420088.0\ny: 4638319.0\nz: 20.8\nazimuth: 72.0\ntilt: 71.6\nroll: -6.6\n"
    tag = Path(d).parent.name                                 # names no station file can have

    def eo(name, notes=""):
        p = cal / name.replace("TAG", tag)
        p.write_text(base + "".join(f"# {n}\n" for n in notes.split("\n") if n))
        return p

    def audit(label, survey, stype, sdate, eos=None, env=None, dem_=dem):
        # the label compare() would print: audit_label's findings, the worst of them (suggested_label)
        f = sc.audit_label(label, str(survey), stype, sdate, sc.load_cameras(eos, {}) if eos else {}, env, str(dem_))
        return sc.suggested_label(label, f), f

    tg = d / "2025-11-13_Marconi_Extrinsic_Targets_c1_xyz.csv"
    tg.write_text("1,420100.0,4638400.0,1.2\n")
    w, f = audit("INDEPENDENT", tg, "points", None, {"c1": eo("CACO05_c1_TAG_EO.yaml", f"solved from {tg.name}")})
    check(w == "PARTLY-CIRCULAR" and not any("not the calibration's own points" in x[1] for x in f)
          and any(x[0] == "PARTLY-CIRCULAR" and x[2] == "rule" and "names these calibration targets" in x[1] for x in f),
          "calibration targets the calibration names in its notes: PARTLY-CIRCULAR (a rule), as its date would say "
          "-- was CIRCULAR, with a note calling them 'not the calibration's own points'")
    w, f = audit("INDEPENDENT", tg, "points", "2025-11-13", {"c1": eo("CACO05_c1_20251113_EO-CV.yaml")})
    check(w == "PARTLY-CIRCULAR" and not any("not the calibration's own" in x[1] for x in f),
          "... by its date alone: PARTLY-CIRCULAR, as before")
    w, f = audit("INDEPENDENT", tg, "points", None, {"c1": eo("CACO05_c1_20250219_EO.yaml")})
    check(w == "INDEPENDENT" and any("not the calibration's own points" in x[1] for x in f),
          "targets, a calibration of another date that does not name them: a note, no downgrade")
    w, f = audit("INDEPENDENT", lid, "dsm", None,
                 {"c1": eo("CACO05_c1_TAG_a_lidar_EO.yaml", "fit_eo_to_survey.py: azimuth +1.0 fitted to the Jan 2025 lidar")})
    check(w == "CIRCULAR" and any(x[2] == "rule" and "not a file name" in x[1] for x in f)
          and not any("another survey" in x[1] for x in f),
          "a pointing 'fitted to the Jan 2025 lidar' (no file named): CIRCULAR as a rule that can misfire -- was "
          "'fitted to another survey (the)', no downgrade")
    w, f = audit("INDEPENDENT", lid, "dsm", None,
                 {"c1": eo("CACO05_c1_TAG_b_lidar_EO.yaml", f"fit_eo_to_survey.py: fitted to {lid.stem}.asc using 40")})
    check(w == "CIRCULAR" and any(x[2] == "chain" for x in f), "fitted to the same survey in another format: CIRCULAR")
    w, f = audit("INDEPENDENT", lid, "dsm", None,
                 {"c1": eo("CACO05_c1_TAG_c_lidar_EO.yaml", "fit_eo_to_survey.py: fitted to other_survey.tif")})
    check(w == "INDEPENDENT" and any("another survey" in x[1] for x in f), "fitted to another survey file: a note")
    # a hand note naming a survey by its stem (no extension) is still a file name: another survey's stem is a
    # note, this survey's stem the chain fact it was
    w, f = audit("INDEPENDENT", lid, "dsm", None,
                 {"c1": eo("CACO05_c1_TAG_g_lidar_EO.yaml", "fit_eo_to_survey.py: fitted to 2025005FA_Marconi_Mar_YSMP")})
    check(w == "INDEPENDENT" and any("another survey" in x[1] for x in f) and not any("not a file name" in x[1] for x in f),
          "fitted to another survey named by its stem (no extension): a note, not 'unresolved'")
    w, f = audit("INDEPENDENT", lid, "dsm", None,
                 {"c1": eo("CACO05_c1_TAG_h_lidar_EO.yaml", f"fit_eo_to_survey.py: fitted to {lid.stem}")})
    check(w == "CIRCULAR" and any(x[0] == "CIRCULAR" and x[2] == "chain" for x in f),
          "fitted to this survey named by its stem: CIRCULAR, a chain fact")
    check([sc.names_file(x) for x in ("the", "Jan", "2025-01-23", "a.tif", "x/y", "2026-09-29_Marconi_Checkshots",
                                      "B.TIFF")] == [False, False, False, True, True, True, True],
          "names_file: words and dates are not file names; names, paths and underscored stems are")
    sub = d / "prov" / "dem"
    sub.mkdir(parents=True)
    shutil.copy(str(dem), str(sub / "x_dem.asc"))
    pj = sub.parent / "provenance.json"
    for prov, want, kind, what in (
            ({"cameras": {"c1": {"calibration": {"history": f"GCPs of that day, independent of {lid.name}"}}}},
             "CIRCULAR", "rule", "a sentence naming it (calibration.history): CIRCULAR as a rule, overrulable"),
            ({"cameras": {"c1": {"calibration": {"eo_path": f"/x/{lid.stem}_EO.yaml"}}}},
             "CIRCULAR", "chain", "a calibration file named after it: CIRCULAR, a chain fact"),
            ({"cameras": {"c1": {"search_envelope": {"survey": str(lid)}}}},
             "PARTLY-CIRCULAR", "chain", "the envelope's survey: PARTLY-CIRCULAR, a chain fact"),
            ({"cameras": {"c1": {"search_envelope": {"how": f"placed with {lid.name}, as the table says"}}}},
             "PARTLY-CIRCULAR", "rule", "the envelope's 'how' sentence: PARTLY-CIRCULAR as a rule"),
            ({"comparisons": [{"survey_path": str(lid), "why": f"{lid.name} fitted nothing"}]},
             "INDEPENDENT", None, "the comparison's own entries: ignored")):
        pj.write_text(json.dumps(prov))
        w, f = audit("INDEPENDENT", lid, "dsm", None, dem_=sub / "x_dem.asc")
        hit = [x for x in f if x[0]]
        check(w == want and (kind is None or (hit and hit[0][2] == kind)), "provenance: " + what)
    pj.unlink()
    rt = d / "rtk.csv"
    rt.write_text("E,N,Z\n420100,4638400,1.2\n")
    w, f = audit("INDEPENDENT", rt, "points", None, {"c1": eo("CACO05_c1_20250220_EO.yaml", "pointing checked with rtk shots")})
    check(w == "INDEPENDENT", "a survey stem of 6 characters or fewer is not looked for inside other words")
    w, f = audit("INDEPENDENT", d / "2025-02-19_profile.csv", "points", "2025-02-19", {"c1": eo("CACO05_c1_20250219_EO.yaml")})
    check(w == "PARTLY-CIRCULAR" and f[0][2] == "rule", "points of the calibration's day, any name: PARTLY-CIRCULAR rule")
    w, f = audit("INDEPENDENT", lid, "dsm", "2025-02-19", {"c1": eo("CACO05_c1_20250219_EO.yaml")})
    check(w == "INDEPENDENT" and f and f[0][0] is None, "a lidar of the calibration's day: a note only")
    w, f = audit("CIRCULAR", lid, "dsm", None, env=str(lid))
    check(w == "CIRCULAR", "passed CIRCULAR, a PARTLY-CIRCULAR finding: stays CIRCULAR (never better)")
    w, f = audit("INDEPENDENT", lid, "dsm", None, {"c1": eo("CACO05_c1_TAG_d_lidar_EO.yaml",
                                                          f"fit_eo_to_survey.py: fitted to {lid.name}")}, env=str(lid))
    check(w == "CIRCULAR" and {x[0] for x in f if x[0]} == {"CIRCULAR", "PARTLY-CIRCULAR"},
          "an envelope (PARTLY) and a pointing (CIRCULAR) both from it: the worst, CIRCULAR")
    # survey_products.honest_label: the same two rules
    gcp = {"name": "t", "survey_type": "points", "survey_date": "", "label": "INDEPENDENT", "why": ""}
    named_eo = eo("CACO05_c1_TAG_e_EO.yaml", f"solved from {tg.name}")
    plan = {"date": "2025-11-13", "era": "live", "first": "2025-11-10", "last": "2025-11-13", "surveys": [], "out": d,
            "cams": {"c1": {"eo_file": named_eo.name, "eo_path": named_eo, "envelope_path": None, "envelope_survey": None}}}
    v = sp.Verdict(gcp)
    sp.honest_label(gcp, plan, SimpleNamespace(setup_coef=0.0, calibration=str(cal)), tg, fits=[], verdict=v)
    check(v.label == "PARTLY-CIRCULAR" and any("names these calibration targets" in r for r in v.reasons),
          "survey_products.honest_label: targets the calibration names -> PARTLY-CIRCULAR (was CIRCULAR)")
    hand = eo("CACO05_c1_TAG_f_lidar_EO.yaml", "fit_eo_to_survey.py: fitted to the Jan 2025 lidar")
    plan["cams"]["c1"].update(eo_file=hand.name, eo_path=hand)
    jan = {"name": "j", "survey_type": "dsm", "survey_date": "2025-01-23", "label": "INDEPENDENT", "why": ""}
    v = sp.Verdict(jan)
    sp.honest_label(jan, plan, SimpleNamespace(setup_coef=0.0, calibration=str(cal)), lid, fits=[], verdict=v)
    v2 = sp.Verdict(dict(jan, label_override_reason="the fit was to the Mar lidar: see its log"))
    sp.honest_label(jan, plan, SimpleNamespace(setup_coef=0.0, calibration=str(cal)), lid, fits=[], verdict=v2)
    check(v.label == "CIRCULAR" and v2.label == "INDEPENDENT" and v2.overridden,
          "... and a pointing 'fitted to the Jan 2025 lidar': CIRCULAR, a rule surveys.csv can overrule")
    # suggested_label, the rule compare() prints: the worst wins, never better than passed, notes and unknown ignored
    fnd = [(None, "a note", "note"), ("PARTLY-CIRCULAR", "p", "rule"), ("NOT-A-LABEL", "x", "rule")]
    check(sc.suggested_label("INDEPENDENT", fnd) == "PARTLY-CIRCULAR"
          and sc.suggested_label("CROSS-VALIDATED", fnd) == "PARTLY-CIRCULAR"
          and sc.suggested_label("CIRCULAR", fnd) == "CIRCULAR"
          and sc.suggested_label("PARTLY-CIRCULAR", [(None, "n", "note")]) == "PARTLY-CIRCULAR"
          and sc.suggested_label("CROSS-VALIDATED", []) == "CROSS-VALIDATED",
          "suggested_label: the worst finding wins (INDEPENDENT or CROSS-VALIDATED + PARTLY -> PARTLY), never "
          "better than passed (CIRCULAR stays), a note or a word that is no label changes nothing")
    try:
        quiet(sc.compare, str(dem), str(tg), "points", "independent", "lower case", output_dir=str(d / "bad"), plot=False)
        bad = False
    except ValueError:
        bad = True
    check(bad, "a label that is not one of LABELS ('independent') is refused, not ranked as the best")
    # end to end: targets on the DEM, the calibration naming them -> compare() warns PARTLY-CIRCULAR, as a rule
    tg2 = d / "2025-11-13_Marconi_Extrinsic_Targets_c2_xyz.csv"
    tg2.write_text("".join(f"{i},{X0 + 1 + 2 * i:.1f},{Y0 + 1 + 2 * i:.1f},1.0\n" for i in range(5)))
    e2 = eo("CACO05_c1_TAG_i_EO.yaml", f"solved from {tg2.name} and 4 more")
    h, _ = quiet(sc.compare, str(dem), str(tg2), "points", "INDEPENDENT", "end to end", name="tg2",
                 output_dir=str(d / "e2e"), camera_eo={"c1": Path(e2)}, plot=False)
    txt = (d / "e2e" / "tg2_comparison.txt").read_text()
    check(h["suggested_label"] == "PARTLY-CIRCULAR" and h["downgrade_warning"] and h["label"] == "INDEPENDENT"
          and "[PARTLY-CIRCULAR, a rule that can misfire]" in txt and "names these calibration targets" in txt
          and "not the calibration's own points" not in txt,
          "compare(): targets the calibration names -> DOWNGRADE WARNING to PARTLY-CIRCULAR, marked a rule; the "
          "label printed as passed")


# ---------------------------------------------------------------------
# Day-block standard error and dem_change.py
# ---------------------------------------------------------------------

def one_cell(values, days):
    """Points of one 2 m cell (one per frame), frames named by index, on the given days."""
    n = len(values)
    E = np.full(n, 420001.0)
    N = np.full(n, 4638001.0)
    frames = np.array([f"f{i}" for i in range(n)])
    return E, N, np.asarray(values, float), frames, {f"f{i}": days[i] for i in range(n)}


def test_day_block_sigma():
    print("dem_from_contours.day_block_sigma: the SE when a day's crossings share an error")
    import dem_from_contours as dfc

    def se_both(values, days):
        E, N, Z, fr, day_of = one_cell(values, days)
        dem, count, spread, e0, n0, nc, nr = dfc.build_grid(E, N, Z, 2.0, 1, 0, fr)
        return (float(dfc.sigma_grid(dem, spread, count)[0, 0]),
                float(dfc.day_block_sigma(E, N, Z, fr, day_of, 2.0, e0, n0, nc, nr, dem, spread, count)[0, 0]),
                float(spread[0, 0]))
    s, sd, _ = se_both([0.0, 0.1, 0.2, 0.35, 0.4], ["a", "b", "c", "d", "e"])
    check(abs(s - sd) < 1e-12, f"one crossing per day (5 days): the independent SE exactly ({s:.4f} = {sd:.4f})")
    s, sd, sp_ = se_both([0.3, 0.32, 0.34, 0.0, 0.02, 0.04, 0.31, 0.03], list("aaabbbcd"))
    s2, sd2, sp2 = se_both([0.3, 0.32, 0.34, 0.0, 0.02, 0.04, 0.31, 0.03], list("abcdefgh"))
    check(sd > 1.4 * s and abs(sd2 - s2) < 0.15 * s2,
          f"days whose crossings sit on one side: {sd:.4f} m against the independent {s:.4f} m; the same values on 8 "
          f"days: {sd2:.4f} (about the independent {s2:.4f})")
    s, sd, sp_ = se_both([1.0, 1.1, 1.2, 0.0, 0.1, 0.2], list("aaabbb"))
    check(abs(sd - 1.25 * (sp_ / 2.0)) < 1e-9,
          f"2 equal days, each wholly on one side: 1.25 sigma / sqrt(k - 1) = {1.25 * sp_ / 2:.4f} ({sd:.4f})")
    s, sd, _ = se_both([0.0, 0.1, 0.2], ["a", "a", "a"])
    check(np.isfinite(s) and np.isnan(sd), "one day: no day-block SE (NaN)")


def ground(path, days, day_bias=0.0, seed=0, first="2026-10-01"):
    """A ground file: per day 3 hours x 2 cameras, frames of 200 points on a 1:50 plane; each day's
    lines share an error (day_bias x a day-to-day normal value)."""
    rng = np.random.RandomState(seed)
    rows = ["source_file,camera,capture_time_utc,capture_epoch,easting_utm19,northing_utm19,"
            "tide_elevation_navd88,beach_elevation_navd88"]
    d0 = datetime.fromisoformat(first + "T00:00:00+00:00").timestamp()
    for day in range(days):
        bias = day_bias * rng.normal()
        for cam, x0 in (("c1", 420000.0), ("c2", 420040.0)):
            for h in (14, 15, 16):
                ep = int(d0 + day * 86400 + h * 3600)
                t = datetime.fromtimestamp(ep, tz=timezone.utc)
                for i in range(200):
                    e, n = x0 + 2.0 * (i % 20) + 1.0, 4638000.0 + 2.0 * (i // 20) + 1.0
                    z = 1.0 - 0.02 * (e - 420000.0) + bias + rng.normal(0, 0.03)
                    rows.append(f"{ep}.{cam}.timex,{cam},{t:%Y-%m-%dT%H:%M:%S}+00:00,{ep},{e:.2f},{n:.2f},"
                                f"{z:.4f},{z:.4f}")
    Path(path).write_text("\n".join(rows) + "\n")
    return path


def test_dem_day_block_se(d):
    print("dem_from_contours.py --day-block-se: a new file, nothing else changes")
    d = Path(d) / "dayse"
    d.mkdir(parents=True)
    g = ground(d / "g.csv", 7, day_bias=0.05, seed=1)
    rc0, log0 = run([HERE / "dem_from_contours.py", g, d / "plain", "--no-plot", "--last-days", "7",
                     "--series-dir", d / "s0"], d)
    rc1, log1 = run([HERE / "dem_from_contours.py", g, d / "day", "--no-plot", "--last-days", "7",
                     "--series-dir", d / "s1", "--day-block-se"], d)
    same = all(np.array_equal(np.loadtxt(d / f"plain_{k}.asc", skiprows=6), np.loadtxt(d / f"day_{k}.asc", skiprows=6))
               for k in ("dem", "spread", "count", "sigma"))
    i0 = json.loads((d / "plain_info.json").read_text())
    i1 = json.loads((d / "day_info.json").read_text())
    check(rc0 == 0 and rc1 == 0 and same and not (d / "plain_sigma_day.asc").exists()
          and {k: v for k, v in i1.items() if k != "sigma_day"} == i0 and "sigma_day" not in i0,
          "the dem, spread, count and sigma grids and the info file are the same with it; without it no "
          "_sigma_day.asc and no 'sigma_day'")
    sd, s = read_asc(d / "day_sigma_day.asc"), read_asc(d / "day_sigma.asc")
    ok = np.isfinite(sd) & np.isfinite(s)
    check(ok.sum() > 50 and np.median(sd[ok] / s[ok]) > 1.3 and i1["sigma_day"]["cells"] == int(np.isfinite(sd).sum())
          and len(list((d / "s1").glob("*_sigma_day.asc"))) == 1 and "day-block SE" in log1,
          f"days sharing a 0.05 m error: the day-block SE {np.median(sd[ok] / s[ok]):.1f} x the independent one "
          f"(median, {int(ok.sum())} cells); in the series copy too; logged")
    rc2, _ = run([HERE / "dem_from_contours.py", g, d / "day", "--no-plot", "--last-days", "7",
                  "--series-dir", d / "s1"], d)
    check(rc2 == 0 and not (d / "day_sigma_day.asc").exists() and not list((d / "s1").glob("*_sigma_day.asc"))
          and (d / "day_sigma.asc").exists(),
          "the same stem rebuilt without it: the earlier _sigma_day.asc (and its series copy) removed, not left "
          "for dem_change.py to read")


def test_dem_change(d):
    print("dem_change.py: standard errors from _sigma.asc; --day-block-lod, report only")
    d = Path(d) / "change"
    d.mkdir(parents=True)
    for name, seed, first in (("A", 2, "2026-10-01"), ("B", 3, "2026-10-08")):
        g = ground(d / f"g{name}.csv", 7, day_bias=0.08, seed=seed, first=first)
        rc, log = run([HERE / "dem_from_contours.py", g, d / name, "--no-plot", "--day-block-se"], d)
        check(rc == 0, f"DEM {name} built")
    # B rises by 0.08 m everywhere: about the level of detection
    b = (d / "B_dem.asc").read_text().splitlines()
    vals = np.loadtxt(d / "B_dem.asc", skiprows=6, ndmin=2)
    vals = np.where(vals == -9999, vals, vals + 0.08)
    with open(d / "B_dem.asc", "w") as f:
        f.write("\n".join(b[:6]) + "\n")
        for r in vals:
            f.write(" ".join("%.4f" % v for v in r) + "\n")
    import dem_change

    def change(tag, *extra):
        out = d / tag
        rc, log = run([HERE / "dem_change.py", "--a", d / "A", "--b", d / "B", "--output-dir", out, "--no-plot"]
                      + list(extra), d)
        return rc, log, read_asc(out / "change_A_to_B_sig.asc"), read_asc(out / "change_A_to_B_diff.asc")
    rc, log, sig, diff = change("files")
    check(rc == 0 and "standard errors         : A from its _sigma.asc; B from its _sigma.asc" in log,
          "both DEMs' _sigma.asc read (and said)")
    shutil.copy(str(d / "A_sigma.asc"), str(d / "A_sigma.keep"))
    (d / "A_sigma.asc").unlink()
    rc2, log2, sig2, diff2 = change("formula")
    check(rc2 == 0 and "A from its spread and count" in log2 and "no A_sigma.asc" in log2
          and np.array_equal(np.isfinite(sig), np.isfinite(sig2)) and np.array_equal(np.nan_to_num(diff), np.nan_to_num(diff2)),
          "without A's file: the formula from its spread and count (said), the same significant cells -- the map "
          "takes the file's SE only where it differs by more than its 4-decimal rounding")
    hdr = [float(x.split()[1]) for x in (d / "A_sigma.keep").read_text().splitlines()[2:5]]   # xll, yll, cell
    write_asc(d / "A_sigma.asc", read_asc(d / "A_sigma.keep") * 3.0, *hdr)
    rc3, log3, sig3, _ = change("tripled")
    check(rc3 == 0 and "differs from 1.25 x spread/2 / sqrt(n)" in log3 and np.isfinite(sig3).sum() < np.isfinite(sig).sum(),
          f"a _sigma.asc three times the formula is used (fewer significant cells: {int(np.isfinite(sig3).sum())} "
          f"against {int(np.isfinite(sig).sum())}) and the log says it differs")
    shutil.copy(str(d / "A_sigma.keep"), str(d / "A_sigma.asc"))
    rc4, log4, sig4, diff4 = change("dayblock", "--day-block-lod")
    blk = log4[log4.index("DAY-BLOCK LEVEL OF DETECTION"):].split("\nwrote")[0] if "DAY-BLOCK" in log4 else ""
    same_files = all((d / "files" / f"change_A_to_B_{k}.asc").read_bytes() ==
                     (d / "dayblock" / f"change_A_to_B_{k}.asc").read_bytes() for k in ("diff", "sig"))
    check(rc4 == 0 and "cells judged" in blk and "significant erosion" in blk and same_files
          and log4.split("DAY-BLOCK")[0].rstrip() == log.split("\nwrote")[0].rstrip(),
          "--day-block-lod: its block printed after the usual lines, which are the same; the _diff and _sig grids "
          "byte for byte the same")
    lod_map = float(re.search(r"level of detection      : median ([0-9.]+)", log).group(1))
    lod_day = float(re.search(r"level of detection    : median ([0-9.]+)", blk).group(1))
    check(lod_day > 1.2 * lod_map,
          f"days sharing a 0.08 m error: the day-block LoD {lod_day:.3f} m over the map's {lod_map:.3f} m")
    (d / "B_sigma_day.asc").unlink()
    rc5, log5, _, _ = change("noday", "--day-block-lod")
    check(rc5 == 0 and "not judged: no B_sigma_day.asc" in log5 and "--day-block-se" in log5,
          "a DEM without _sigma_day.asc: not judged, and how to make it")
    # day-block SEs far below the map's: the LoD stays the map's (the floor); cells of one day (NaN) not judged
    sA_, sB_ = read_asc(d / "A_sigma.keep"), read_asc(d / "B_sigma.asc")
    sA_[:3] = np.nan
    write_asc(d / "A_sigma_day.asc", 0.01 * sA_, *hdr)
    write_asc(d / "B_sigma_day.asc", 0.01 * sB_,
              *[float(x.split()[1]) for x in (d / "B_sigma.asc").read_text().splitlines()[2:5]])
    rc8, log8, _, diff8 = change("tiny", "--day-block-lod", "--min-lod", "0")
    blk8 = log8[log8.index("DAY-BLOCK LEVEL OF DETECTION"):] if "DAY-BLOCK" in log8 else ""
    judged = int((np.isfinite(diff8) & np.isfinite(sA_) & np.isfinite(sB_)).sum())
    n_cmp = int(np.isfinite(diff8).sum())
    check(rc8 == 0 and "1.00 x the map's" in blk8 and 0 < judged < n_cmp
          and f"cells judged          : {judged} of the {n_cmp} compared ({n_cmp - judged} rest on one day" in blk8,
          f"a _sigma_day.asc at 1% of the map's SE: the LoD stays the map's (1.00 x); {judged} of {n_cmp} cells judged, "
          f"the {n_cmp - judged} with no day-block SE not")
    # --series --rebuild: the reference rebuilt gets --day-block-se itself
    g = ground(d / "gall.csv", 14, day_bias=0.04, seed=4)
    rc6, log6 = run([HERE / "dem_from_contours.py", g, d / "win", "--no-plot", "--last-days", "7", "--series-dir",
                     d / "series", "--day-block-se"], d)
    newb = sorted((d / "series").glob("dem_*_7d_dem.asc"))[0]
    rc7, log7 = run([HERE / "dem_change.py", "--series", d / "series", "--days", "7", "--b",
                     str(newb)[:-len("_dem.asc")], "--rebuild", g, "--no-plot", "--min-ref-cells", "1",
                     "--min-ref-range", "0", "--day-block-lod", "--", "--cell", "2"], d)
    check(rc6 == 0 and rc7 == 0 and "DAY-BLOCK LEVEL OF DETECTION" in log7 and "cells judged" in log7
          and "not judged" not in log7,
          "--series --rebuild --day-block-lod: the rebuilt reference has its day-block SE, the block is printed")
    check(dem_change.load_se(d / "nothere", {"ncols": 1, "nrows": 1, "xllcorner": 0, "yllcorner": 0, "cellsize": 2})
          == (None, "no nothere_sigma.asc"), "load_se: no file, said")
    # the default map is the formula's to the bit: a _sigma.asc of the same definition (4 decimals)
    # flips no borderline cell. B = A + noise of about the LoD, rounded to 4 dp, A's spread, count, SE
    nz = d / "noisy"
    for sub in ("files", "formula"):
        (nz / sub).mkdir(parents=True)
    hdrA = [float(x.split()[1]) for x in (d / "A_dem.asc").read_text().splitlines()[2:5]]
    a_dem, a_sp, a_n, a_se = (read_asc(d / f"A_{k}.asc") for k in ("dem", "spread", "count", "sigma"))
    se_f = 1.25 * (np.nan_to_num(a_sp) / 2.0) / np.sqrt(np.maximum(np.nan_to_num(a_n), 1))
    lod_f = 1.96 * np.sqrt(2.0) * se_f
    lod_g = 1.96 * np.sqrt(2.0) * np.where(np.isfinite(a_se), a_se, se_f)
    same, would_flip = [], 0
    for seed in range(8):
        rng = np.random.RandomState(seed)
        b_dem = np.round(a_dem + np.where(rng.rand(*a_dem.shape) < 0.5, -1, 1) * lod_f
                         * rng.uniform(0.97, 1.03, a_dem.shape), 4)
        with np.errstate(invalid="ignore"):
            use = np.isfinite(a_dem) & (np.nan_to_num(a_n) >= 5)
            dd = np.abs(np.round(b_dem, 4) - a_dem)
            would_flip += int((use & ((dd > lod_f) != (dd > lod_g))).sum())
        for sub in ("files", "formula"):
            for nm, dem in (("A", a_dem), ("N", b_dem)):
                write_asc(nz / sub / f"{nm}_dem.asc", dem, *hdrA)
                for k in ("spread", "count") + (("sigma",) if sub == "files" else ()):
                    shutil.copy(str(d / f"A_{k}.asc"), str(nz / sub / f"{nm}_{k}.asc"))
            quiet(dem_change.compare, str(nz / sub / "A"), str(nz / sub / "N"), 0.0,
                  out_dir=str(nz / sub / f"out{seed}"), plot=False)
        same.append(all((nz / "files" / f"out{seed}" / f"change_A_to_N_{k}.asc").read_bytes()
                        == (nz / "formula" / f"out{seed}" / f"change_A_to_N_{k}.asc").read_bytes()
                        for k in ("sig", "diff")))
    check(all(same) and would_flip > 0,
          f"8 noisy B grids at the LoD (--min-lod 0): _sig and _diff grids byte for byte the same with and "
          f"without the _sigma.asc files ({sum(same)} of 8), where the file's rounded SE alone would flip "
          f"{would_flip} borderline cell(s)")


def test_py38():
    print("Python 3.8 syntax")
    import ast
    for f in ("survey_compare.py", "dem_change.py", "dem_from_contours.py", "survey_products.py",
              "test_validation_stats.py"):
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
    d = os.path.realpath(tempfile.mkdtemp(prefix="validation_stats_", dir=args.keep))
    import survey_products as sp
    station = (sp.CHELSEA, sp.WATERLINE, sp.DEFAULT_SURVEY_DIRS)
    sp.CHELSEA, sp.WATERLINE = Path(d) / "station_chelsea", Path(d) / "station_waterline"
    sp.DEFAULT_SURVEY_DIRS = [str(sp.CHELSEA), str(sp.WATERLINE)]
    for x in sp.DEFAULT_SURVEY_DIRS:
        Path(x).mkdir()
    start = os.getcwd()
    os.chdir(d)
    print(f"test folder: {d} (the current folder during the run)")
    try:
        test_boot_stats_exact()
        test_block_intervals()
        test_contrast()
        test_correlogram()
        test_dem_intervals()
        test_alongshore_axis()
        test_strata()
        test_compare(d)
        test_labels(d)
        test_day_block_sigma()
        test_dem_day_block_se(d)
        test_dem_change(d)
        test_py38()
    finally:
        os.chdir(start)
        sp.CHELSEA, sp.WATERLINE, sp.DEFAULT_SURVEY_DIRS = station
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED")
        return 1
    print("all validation statistics tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

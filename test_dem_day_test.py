#!/usr/bin/env python3
"""
Tests For The DEM Day Test (dem_from_contours.py --max-day-offset)
=====================================================================
Plain python3 (no pytest needed on the station):

    python3 test_dem_day_test.py [--keep DIR] [--year]

WHY. The day test decides which camera-days make every DEM: the cron's
whole-archive and 7-day DEMs, dem_change.py's rebuilt reference week and
the survey-date DEMs. Run once (each day against the DEM of all the
others) it let the days it was about to reject shape that DEM: on the
station's whole-archive DEM of 13 Oct 2026 c1's pre-storm days pulled
c1's cells up and four post-storm c1 days were rejected as LOW, and on
the 7-day DEM c1's steadily rising offsets lost the low end only. The
test is now iterative (leave out the day furthest over the limit,
measure the others again without it, repeat; then take back what lies
within the limit of the kept days) and warns when a camera's offsets run
one way with time. These checks plant the answer and ask for it back:
  * one bad day among good ones: only it is left out;
  * two bad days that mask each other in one pass: both are found;
  * a cluster of off days that swamps good ones in one pass (the
    whole-archive case): the good days are kept;
  * a day left out early because a worse day shaped its DEM: taken back;
    one whose return would put a kept day over: kept out, and recorded;
  * a fringe day whose cells only one other day reaches, that day being
    checked by the rest: left out (as one pass did);
  * days the test cannot judge (two days that only see each other, the
    last two days that disagree): the one-pass verdict, and a WARNING;
    every day disagreeing: every day rejected, exit 4, as before;
  * the WARNING and <stem>_info.json describe the final DEM: a day left
    out unjudged and later taken back, or kept unjudged and within the
    limit in the end, is not named; one the kept days judge in the end
    is not either; one kept unjudged is tried again, and left out once
    the test can judge it; no other log line holds the word WARNING;
  * a steady trend: the WARNING, and what the test trims; a trend trimmed
    at one end of a camera's 5 days: the WARNING, saying it rests on that
    day; days without a trend, or one bad day at an end: the WARNING stays
    rare; the chance of a rho exact for tied offsets;
  * the station's 7-day c1/c2 pattern of 5-11 Oct 2026, replayed;
  * fewer than 3 days: none left out;
  * a survey-date DEM (--anchor-date) whose days around the survey are
    outweighed by the others: they are kept and the others, off them, left
    out; no day near the date: as without;
  * the offsets equal the one-pass code's, whatever the order of the
    samples, ties are broken by name, and random scenarios (stretches of
    beach, and two-camera tidal windows) keep the test's bookkeeping;
  * run time: a quarter-year archive with many days left out in seconds
    (--year: a year, 730 camera-days);
  * the command line: the log block survey_products.py reads, and
    <stem>_info.json (camera_days_rejected as before, plus each left-out
    day's offset and the test's numbers).

HOW. Synthetic (cell, frame) samples: elevation = 0.01 x cell + the day's
bias + noise (seeded). The checks call dem_from_contours.day_consistency
directly, and run dem_from_contours.py on a synthetic ground file in a
temporary folder (--keep DIR to look at it).

STATION FILES. None are read.
"""

import os
import io
import re
import sys
import json
import shutil
import argparse
import tempfile
import contextlib
import subprocess
from pathlib import Path
from datetime import datetime, timezone

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import dem_from_contours as dfc                               # noqa: E402

LIMIT = 0.15
MIN_FRAMES = 3


def quiet(fn, *a, **k):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = fn(*a, **k)
    return out, buf.getvalue()


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)
    print(f"  ok  {msg}")


# ---------------------------------------------------------------------
# synthetic camera-days

def pairs(days, sigma=0.02, seed=0):
    """days: [(key, bias, cells, frames per cell)] -> (pair_cell, pair_vals, pair_key): one sample
    per frame per cell, elevation 0.01 x cell + bias + N(0, sigma)."""
    rng = np.random.RandomState(seed)
    pc, pv, pk = [], [], []
    for key, bias, cells, nf in days:
        cells = np.asarray(cells, np.int64)
        for _ in range(nf):
            pc.append(cells)
            pv.append(0.01 * cells + bias + rng.normal(0.0, sigma, len(cells)))
            pk.append(np.array([key] * len(cells)))
    return np.concatenate(pc), np.concatenate(pv), np.concatenate(pk)


def run(days, limit=LIMIT, **k):
    pc, pv, pk = pairs(days, **k)
    return dfc.day_consistency(pc, pv, pk, limit, MIN_FRAMES), dfc.day_offsets(pc, pv, pk, MIN_FRAMES)


def one_pass_rejects(one, limit=LIMIT):
    """What the one-pass test (before Oct 2026) rejected: every testable day over the limit."""
    return sorted(k for k, (o, n) in one.items() if np.isfinite(o) and abs(o) > limit)


def scenario_one_bad():
    """Eight c2 days over the same 100 cells; 4 Mar reads 0.30 m high."""
    return [(f"c2 2025-03-0{d}", 0.30 if d == 4 else 0.0, range(100), 3) for d in range(1, 9)]


def scenario_masking():
    """Six good days over the mid beach (cells 100-199), three of them also on the low beach (cells
    0-59, one frame each). Two bad days only on the low beach, three frames each: 7 Mar +0.60 m, 8 Mar
    +0.20 m. In one pass 8 Mar's DEM there is 7 Mar + the three good days, whose median sits halfway
    (~+0.3 m): 8 Mar reads ~-0.1 m and is kept, masked by 7 Mar."""
    days = []
    for d in range(1, 7):
        days.append((f"c2 2025-03-0{d}", 0.0, list(range(100, 200)) + (list(range(60)) if d <= 3 else []), 1))
        days.append((f"c2 2025-03-0{d}", 0.0, range(100, 200), 1))
    days.append(("c2 2025-03-07", 0.60, range(60), 3))
    days.append(("c2 2025-03-08", 0.20, range(60), 3))
    return days


PRE_STORM = {"2026-09-11": 0.36, "2026-09-12": 0.44, "2026-09-14": 0.27, "2026-09-15": 0.47,
             "2026-09-18": 1.01, "2026-09-19": 0.88, "2026-09-20": 1.14, "2026-09-21": 1.40}
POST_STORM = ["2026-09-29", "2026-09-30", "2026-10-01", "2026-10-02", "2026-10-03", "2026-10-04",
              "2026-10-05", "2026-10-06", "2026-10-07", "2026-10-09", "2026-10-10", "2026-10-11"]
EARLY_OCT = ["2026-10-04", "2026-10-05", "2026-10-06", "2026-10-07"]


def scenario_swamp():
    """The station's whole-archive c1 case, in miniature. 8 pre-storm days (11-21 Sep) read +0.27 to
    +1.40 m against the post-storm beach (the NUC's first-pass numbers), one frame per cell over cells
    0-349. 12 post-storm days (29 Sep - 11 Oct) agree within +/-0.05 m, two frames per cell over the
    upper beach (cells 0-149); four of them (4-7 Oct) also reach the lower beach (cells 150-349), where
    their DEM is mostly pre-storm days. One pass rejects those four as LOW (-0.28 to -0.33 m), as the
    NUC did (-0.15 to -0.26 m)."""
    rng = np.random.RandomState(7)
    days = [(f"c1 {d}", b, range(350), 1) for d, b in PRE_STORM.items()]
    for d in POST_STORM:
        b = float(rng.uniform(-0.05, 0.05))
        days.append((f"c1 {d}", b, list(range(150)) + (list(range(150, 350)) if d in EARLY_OCT else []), 2))
    return days


def scenario_take_back():
    """3 Mar (-0.05 m: fine) shares its cells (0-49) with 5 Mar (+0.25 m, six frames per cell), which
    then makes its DEM: 3 Mar reads -0.30, worse than 5 Mar's +0.25, and goes first. Without 5 Mar it
    is -0.05 m off the kept days: taken back."""
    days = [(f"c2 2025-03-0{d}", 0.0, list(range(50)) + list(range(50, 150)), 1) for d in (1, 2, 4, 6)]
    days += [(f"c2 2025-03-0{d}", 0.0, range(50, 150), 1) for d in (1, 2, 4, 6)]
    days.append(("c2 2025-03-03", -0.05, range(50), 2))
    days.append(("c2 2025-03-05", 0.25, range(50), 6))
    return days


def scenario_refused():
    """On cells 0-49: four good days (one frame each, 0 m), 2 Mar at +0.12 m (kept: within the limit),
    3 Mar at -0.13 m with five frames per cell, and 5 Mar at +0.30 m with six, which first makes 3 Mar
    read -0.42 m: 3 Mar goes first, then 5 Mar. Against the kept days 3 Mar is then -0.13 m, within
    the limit, but its five frames would make the DEM 2 Mar is judged by -0.13 m too, putting 2 Mar
    over the limit: it is not taken back."""
    days = [(f"c2 2025-03-0{d}", 0.0, list(range(50)), 1) for d in (1, 4, 6, 7)]
    days += [(f"c2 2025-03-0{d}", 0.0, range(50, 150), 2) for d in (1, 4, 6, 7)]
    days += [("c2 2025-03-02", 0.12, range(50), 1), ("c2 2025-03-03", -0.13, range(50), 5),
             ("c2 2025-03-05", 0.30, range(50), 6)]
    return days


TREND_DATES = ["2026-10-05", "2026-10-06", "2026-10-07", "2026-10-08", "2026-10-09", "2026-10-10", "2026-10-11"]


def scenario_trend():
    """c1 rises 0.06 m a day for seven days (5-11 Oct, -0.18 to +0.18 m: a beach building up, or a
    camera drifting); c2 is steady. Five frames per cell, 0.12 m noise per sample (as on the NUC),
    c1 on cells 0-119 and c2 on 100-299."""
    days = [(f"c1 {d}", -0.18 + 0.06 * i, range(120), 5) for i, d in enumerate(TREND_DATES)]
    days += [(f"c2 {d}", 0.0, range(100, 300), 5) for d in TREND_DATES]
    return days


# The station's 7-day DEM, 5-11 Oct 2026 (NUC cron of 13 Oct 2026, last_run.log): the first-pass
# offsets and comparable samples it printed. replay_nuc() finds day biases that give the offsets back
# on a synthetic beach laid out like the station's: c1 on cells 0-119, c2 on 100-299 (a seam of 20
# cells), frames per cell from the NUC's sample counts, noise 0.12 m per sample (the NUC DEM's median
# 16-84 % spread was 0.248 m, i.e. about 2 x 0.12).
NUC_7DAY = {"c1 2026-10-05": -0.238, "c1 2026-10-06": -0.058, "c1 2026-10-07": -0.048,
            "c1 2026-10-09": -0.034, "c1 2026-10-10": 0.016, "c1 2026-10-11": 0.132,
            "c2 2026-10-05": -0.002, "c2 2026-10-06": -0.032, "c2 2026-10-07": -0.058,
            "c2 2026-10-08": 0.183, "c2 2026-10-09": 0.037, "c2 2026-10-10": -0.066,
            "c2 2026-10-11": -0.011}
NUC_7DAY_SAMPLES = {"c1 2026-10-05": 798, "c1 2026-10-06": 983, "c1 2026-10-07": 748, "c1 2026-10-09": 648,
                    "c1 2026-10-10": 1049, "c1 2026-10-11": 862, "c2 2026-10-05": 1021, "c2 2026-10-06": 2176,
                    "c2 2026-10-07": 3167, "c2 2026-10-08": 2948, "c2 2026-10-09": 2873, "c2 2026-10-10": 3105,
                    "c2 2026-10-11": 3436}
NUC_SIGMA = 0.12


def nuc_days(bias):
    out = []
    for k, b in sorted(bias.items()):
        cells = range(0, 120) if k.startswith("c1") else range(100, 300)
        out.append((k, b, cells, max(1, int(round(NUC_7DAY_SAMPLES[k] / len(cells))))))
    return out


def replay_nuc(rounds=40):
    """Biases whose first-pass offsets come back as NUC_7DAY, by fixed-point iteration (offsets do not
    see a shift common to every day, so the biases are kept centred on 0). The layout cannot give every
    pattern exactly: each camera's offsets come back within ~5 mm. -> (bias, days, largest error)."""
    bias = dict(NUC_7DAY)
    best = (np.inf, bias)
    for _ in range(rounds):
        pc, pv, pk = pairs(nuc_days(bias), sigma=NUC_SIGMA, seed=3)
        first = dfc.day_offsets(pc, pv, pk, MIN_FRAMES)
        err = {k: NUC_7DAY[k] - first[k][0] for k in bias}
        worst = max(abs(e) for e in err.values())
        if worst < best[0]:
            best = (worst, dict(bias))
        bias = {k: bias[k] + err[k] for k in bias}
        m = float(np.mean(list(bias.values())))
        bias = {k: v - m for k, v in bias.items()}
    return best[1], nuc_days(best[1]), best[0]


# ---------------------------------------------------------------------

def test_engine():
    print("offsets: the same as the one-pass code, whatever the order of the samples")
    pc, pv, pk = pairs(scenario_masking(), seed=11)

    def old_day_offsets(pair_cell, pair_vals, pair_key, min_frames, min_pairs=30):
        out = {}                                  # the code before Oct 2026, verbatim
        for k in sorted(set(pair_key)):
            mine = pair_key == k
            order, c, starts, ends = dfc.cell_groups(pair_cell[~mine])
            others = pair_vals[~mine][order]
            med = {}
            for a, b in zip(starts, ends):
                if b - a >= min_frames:
                    med[int(c[a])] = float(np.median(others[a:b]))
            res = [v - med[int(cc)] for cc, v in zip(pair_cell[mine], pair_vals[mine]) if int(cc) in med]
            out[k] = (float(np.median(res)) if len(res) >= min_pairs else np.nan, len(res))
        return out
    old = old_day_offsets(pc, pv, pk, MIN_FRAMES)
    new = dfc.day_offsets(pc, pv, pk, MIN_FRAMES)
    check(list(old) == list(new) and all(old[k] == new[k] for k in old),
          f"day_offsets equals the one-pass code exactly ({len(new)} days)")
    p = np.random.RandomState(5).permutation(len(pc))
    check(dfc.day_offsets(pc[p], pv[p], pk[p], MIN_FRAMES) == new, "... and does not depend on the sample order")
    frames = np.array([f"{k}#{i % 7}" for i, k in enumerate(pk)])
    check(dfc.day_offsets(pc, pv, frames, MIN_FRAMES, min_pairs=10) ==
          old_day_offsets(pc, pv, frames, MIN_FRAMES, min_pairs=10),
          "... and for frames (the --max-frame-offset test), min_pairs 10")
    loo = dfc.LeaveOneOut(pc, pv, pk, MIN_FRAMES)
    ref = np.array([n not in ("c2 2025-03-07",) for n in loo.names])
    sub = pk != "c2 2025-03-07"
    check(loo.offset("c2 2025-03-08", ref) == dfc.day_offsets(pc[sub], pv[sub], pk[sub], MIN_FRAMES)["c2 2025-03-08"],
          "an offset against a subset of days equals the one-pass offset with the other days removed")
    every = np.ones(len(loo.names), bool)
    j1, j2 = loo.judge("c2 2025-03-08", every, LIMIT), loo.judge("c2 2025-03-08", ref, LIMIT)
    check(j1["judged"] and j1["share"] == 1.0 and j2["judged"] and j2["share"] == 1.0 and not j1["alone"],
          "judge: on the low beach 8 Mar's DEM holds 7 Mar and three good days, three without 7 Mar: "
          "two or more other days in every cell, so the test can tell which is off")
    # offsets(): every key at once, the same to the bit as offset() one key at a time, against any
    # set of reference keys (a key in it leaves its own samples out), here and on random data with ties
    rng = np.random.RandomState(8)
    cases = [(pc, pv, pk, MIN_FRAMES, 30)]
    for _ in range(40):
        n = int(rng.randint(50, 600))
        v = rng.normal(0, 0.3, n)
        cases.append((rng.randint(0, 40, n).astype(np.int64), np.round(v, 1) if rng.rand() < 0.5 else v,
                      np.array(["k%d" % i for i in rng.randint(0, 9, n)]), int(rng.randint(1, 5)),
                      int(rng.randint(1, 30))))
    same, tried = True, 0
    for c_, v_, k_, mf, mp in cases:
        lo = dfc.LeaveOneOut(c_, v_, k_, mf, mp)
        for _ in range(4):
            r_, t_ = rng.rand(len(lo.names)) < 0.7, rng.rand(len(lo.names)) < 0.8
            got = lo.offsets(r_, t_)
            for i, name in enumerate(lo.names.tolist()):
                if t_[i]:
                    a, b = got[name], lo.offset(name, r_)
                    tried += 1
                    same &= a[1] == b[1] and (a[0] == b[0] or (np.isnan(a[0]) and np.isnan(b[0])))
                else:
                    same &= name not in got
    check(same, f"offsets() equals offset() to the bit for every key and reference set ({tried} offsets, "
                f"{len(cases)} sample sets, ties included)")


def test_spearman():
    print("Spearman's rho by ranks (no scipy)")
    check(list(dfc.ranks([3.0, 1.0, 2.0, 2.0])) == [4.0, 1.0, 2.5, 2.5], "ranks: ties share their mean rank")
    check(dfc.spearman([1, 2, 3, 4, 5], [10, 20, 30, 40, 50]) == 1.0
          and dfc.spearman([1, 2, 3, 4, 5], [5, 4, 3, 2, 1]) == -1.0, "rho +1 / -1 for monotonic series")
    rho = dfc.spearman([1, 2, 3, 4, 5], [1, 3, 2, 5, 4])
    check(abs(rho - 0.8) < 1e-12, f"rho of 1 3 2 5 4 = 1 - 6*4/(5*24) = 0.8 ({rho})")
    rho = dfc.spearman([1, 2, 3, 4, 5, 6], [1, 2, 2, 3, 5, 4])           # Pearson's r of the averaged ranks
    rx, ry = np.arange(1, 7.0), np.array([1, 2.5, 2.5, 4, 6, 5])
    check(abs(rho - np.corrcoef(rx, ry)[0, 1]) < 1e-12, f"with ties: Pearson's r of the averaged ranks ({rho:.4f})")
    check(np.isnan(dfc.spearman([1, 2, 3], [2, 2, 2])), "a constant side: NaN, never a trend")
    # the chance of a rho by random order, exact for the values' own ties: their tie pattern is kept
    # in every order (5 days with a tied pair cannot reach rho 1; their largest rho comes in 2 orders)
    from itertools import permutations
    y = [-0.2, -0.1, 0.0, 0.0, 0.2]
    rho = dfc.spearman(np.arange(5), y)
    brute = np.mean([dfc.spearman(np.arange(5), [y[i] for i in q]) >= rho - 1e-9 for q in permutations(range(5))])
    check(abs(rho - 0.975) < 1e-3 and dfc.spearman_tail(rho, 5, y) == brute == 2 / 120.0
          and dfc.spearman_tail(rho, 5) == 1 / 120.0,
          f"one tie over 5 days: rho {rho:.3f}, reached by 2 of the 120 orders (the tied pair either way): "
          f"{dfc.spearman_tail(rho, 5, y):.4f} (the no-tie table said 1/120)")
    off = {f"c1 2026-10-0{i}": o for i, o in enumerate(y, 1)}
    two = {f"c1 2026-10-0{i}": o for i, o in enumerate([-0.2, -0.2, -0.1, -0.1, 0.2], 1)}
    t1, t2 = dfc.day_trends(off, LIMIT), dfc.day_trends(two, LIMIT)
    check(len(t1) == 1 and abs(t1[0]["chance"] - 1 / 60.0) < 1e-12 and t2 == []
          and dfc.spearman_tail(dfc.spearman(np.arange(5), list(two.values())), 5, list(two.values())) == 4 / 120.0,
          "5 days with one tie: still a WARNING (1 in 60, within 2 %); with two ties: 4 orders in 120 (3.3 %), none")


def test_one_bad():
    print("one bad day among good ones: only it is left out")
    dt, one = run(scenario_one_bad())
    check(one_pass_rejects(one) == ["c2 2025-03-04"], "the one-pass test finds it too")
    check(dt["rejected"] == ["c2 2025-03-04"] and [s[0] for s in dt["steps"]] == ["c2 2025-03-04"]
          and not dt["taken_back"], "left out alone, not taken back")
    check(abs(dt["final"]["c2 2025-03-04"][0] - 0.30) < 0.02,
          f"its final offset (from the kept days) is the planted +0.30 m ({dt['final']['c2 2025-03-04'][0]:+.3f})")
    check(all(abs(dt["final"][k][0]) < 0.02 for k in dt["final"] if k != "c2 2025-03-04"),
          "the kept days sit within 0.02 m of their DEM")
    check(dfc.day_trends({k: v[0] for k, v in dt["final"].items()}, LIMIT) == [], "no trend warning")


def test_masking():
    print("two bad days that mask each other in one pass: both found")
    dt, one = run(scenario_masking())
    o8 = one["c2 2025-03-08"][0]
    check(one_pass_rejects(one) == ["c2 2025-03-07"] and abs(o8) < LIMIT,
          f"one pass: only 7 Mar; 8 Mar (planted +0.20 m) reads {o8:+.3f} m and is kept")
    check([s[0] for s in dt["steps"]] == ["c2 2025-03-07", "c2 2025-03-08"]
          and dt["rejected"] == ["c2 2025-03-07", "c2 2025-03-08"],
          "iterative: 7 Mar first, then 8 Mar (measured without 7 Mar)")
    check(abs(dt["final"]["c2 2025-03-08"][0] - 0.20) < 0.02 and abs(dt["final"]["c2 2025-03-07"][0] - 0.60) < 0.02,
          f"final offsets are the planted ones: 7 Mar {dt['final']['c2 2025-03-07'][0]:+.3f}, "
          f"8 Mar {dt['final']['c2 2025-03-08'][0]:+.3f} m")
    check(all(abs(dt["final"][k][0]) < LIMIT for k in dt["final"] if k not in dt["rejected"]),
          "every kept day within the limit")


def test_swamp():
    print("a cluster of off days swamps good ones in one pass (the station's whole-archive DEM)")
    days = scenario_swamp()
    dt, one = run(days)
    pre = sorted("c1 " + d for d in PRE_STORM)
    post = sorted("c1 " + d for d in POST_STORM)
    lost = [k for k in one_pass_rejects(one) if k in post]
    check(lost == ["c1 " + d for d in EARLY_OCT] and set(pre) <= set(one_pass_rejects(one)),
          f"one pass: the 8 pre-storm days AND 4 post-storm days rejected, these as LOW "
          f"({', '.join('%s %+.2f' % (k[3:], one[k][0]) for k in lost)})")
    check(dt["rejected"] == pre and not dt["taken_back"],
          "iterative: exactly the 8 pre-storm days left out (largest first), every post-storm day kept")
    check(all(abs(dt["final"][k][0]) < 0.08 for k in post),
          "the post-storm days (planted within +/-0.05 m) sit within 0.08 m of the kept DEM")
    tr = dfc.day_trends({k: v[0] for k, v in dt["final"].items()}, LIMIT)
    rho = dfc.spearman(np.arange(20), [dt["final"][k][0] for k in sorted(dt["final"])])
    check(tr == [] and abs(rho) < dfc.DAY_TREND_MIN_RHO,
          f"a step (before/after a storm) is not a steady trend: rho {rho:+.2f}, no trend WARNING")


def test_take_back():
    print("a day left out because a worse day shaped its DEM is taken back")
    dt, one = run(scenario_take_back())
    check(one_pass_rejects(one) == ["c2 2025-03-03", "c2 2025-03-05"],
          f"one pass: both rejected (3 Mar {one['c2 2025-03-03'][0]:+.3f}, 5 Mar {one['c2 2025-03-05'][0]:+.3f})")
    check([s[0] for s in dt["steps"]] == ["c2 2025-03-03", "c2 2025-03-05"],
          "iterative: 3 Mar goes first (|-0.30| > |+0.25|), then 5 Mar")
    check(dt["taken_back"] == ["c2 2025-03-03"] and dt["rejected"] == ["c2 2025-03-05"],
          f"3 Mar taken back ({dt['final']['c2 2025-03-03'][0]:+.3f} m from the kept days); 5 Mar stays out")
    lines = dfc.day_test_report(dt, LIMIT, 1234)
    row = [ln for ln in lines if ln.startswith("   c2 2025-03-03")][0]
    check("REJECTED" not in row and "taken back" in row, f"its line says so: '{row.strip()}'")
    dt, one = run(scenario_refused())
    check([s[0] for s in dt["steps"]] == ["c2 2025-03-03", "c2 2025-03-05"] and dt["taken_back"] == []
          and dt["rejected"] == ["c2 2025-03-03", "c2 2025-03-05"],
          "a left-out day within the limit of the kept days stays out when taking it back puts a kept day over")
    o3, o2 = dt["final"]["c2 2025-03-03"][0], dt["final"]["c2 2025-03-02"][0]
    check(abs(o3) <= LIMIT and [r[0] for r in dt["refused"]] == ["c2 2025-03-03"]
          and dt["refused"][0][2] == ["c2 2025-03-02"] and abs(o2) <= LIMIT,
          f"3 Mar is {o3:+.3f} m from the kept days, but would put 2 Mar ({o2:+.3f} m) over the limit")
    text = " ".join(x.strip() for x in dfc.day_test_report(dt, LIMIT, 1))
    check(sum(x.strip().startswith("not taken back:") for x in dfc.day_test_report(dt, LIMIT, 1)) == 1
          and "not taken back: c2 2025-03-03" in text
          and "would put c2 2025-03-02 over the limit" in text, "... and the lines say so, once")
    lines = dfc.day_test_report(dt, LIMIT, 1)
    row = [ln for ln in lines if ln.startswith("   c2 2025-03-03")][0]
    check("REJECTED" in row and "within the limit of the kept days, not taken back" in row
          and "taken back: none (the days within the limit of the kept days would each put a kept day over it"
          in text and "every day left out is still over the limit, or not comparable" not in text.split(
              "taken back: none")[1][:60],
          "its line and the 'taken back: none' line do not claim every day left out is over the limit")
    rec = dfc.day_test_record(dt, LIMIT, [])["day_test"]
    check(rec["refused_take_back"] == [{"day": "c2 2025-03-03", "offset_m": round(o3, 4),
                                        "would_push": ["c2 2025-03-02"]}],
          "<stem>_info.json: day_test.refused_take_back names it, its offset and the day it would push")
    import dem_figure
    info = {"filters": {"max_day_offset": LIMIT, "camera_days_rejected": dt["rejected"],
                        **dfc.day_test_record(dt, LIMIT, [])}}
    left = dem_figure.title_lines(info, {"filled": 1, "blanked": 0, "no_value": 0}, "TEST", 0.5, 3, 2.0)[2]
    check(any("2 camera-days left out by the day test (±0.15 m; 1 within it but at odds with a kept day)" in x
              for x in left) and not any("off the rest by" in x for x in left),
          "the page does not call it 'off the rest by > 0.15 m'")
    check(all(abs(dt["final"][k][0]) <= LIMIT for k in dt["final"] if k not in dt["rejected"]),
          "every kept day within the limit")


def test_trend():
    print("a steady trend: WARNING, and what the test trims")
    dt, one = run(scenario_trend(), sigma=NUC_SIGMA, seed=5)
    tr = dfc.day_trends({k: v[0] for k, v in dt["final"].items()}, LIMIT, left_out=dt["rejected"])
    check(len(tr) == 1 and tr[0]["camera"] == "c1" and tr[0]["spearman_rho"] == 1.0 and tr[0]["days"] == 7
          and tr[0]["without"] is None,
          f"c1 flagged: rho {tr[0]['spearman_rho'] if tr else '-'}, range {tr[0]['range_m'] if tr else '-'} m "
          f"over 7 days; c2 not")
    c1 = sorted(k for k in dt["first"] if k.startswith("c1"))
    gone = [s[0] for s in dt["steps"]]
    check(dt["rejected"] and set(dt["rejected"]) <= set(c1) and set(dt["rejected"]) & {c1[0], c1[-1]},
          f"the test trims the trend from its ends: left out {', '.join(g[3:] for g in gone)} (in that order; "
          f"one pass: {', '.join(k[3:] for k in one_pass_rejects(one)) or 'none'}), "
          f"taken back: {', '.join(k[3:] for k in dt['taken_back']) or 'none'}")
    kept = [k for k in c1 if k not in dt["rejected"]]
    check(all(abs(dt["final"][k][0]) <= LIMIT for k in kept),
          f"what is left of c1 ({', '.join(k[3:] for k in kept)}) agrees within the limit: the DEM keeps "
          f"part of the change and loses the rest -- hence the WARNING")
    w = dfc.day_trend_warning(tr[0], LIMIT, dt["rejected"])
    rest = " ".join(x.strip() for x in w[1:])
    check(w[0].startswith("WARNING: c1's camera-day offsets run one way with time over 7 days")
          and all(k in rest for k in dt["rejected"]) and "drifts with time" in rest
          and sum("WARNING" in x for x in w) == 1 and all(x.startswith("   ") for x in w[1:])
          and "The method alone drifts so: on a static synthetic beach" in rest and "rho -0.96" in rest,
          "one WARNING line (the station status counts them), naming the days left out, and what the method alone "
          "gives on a static synthetic beach (the same warning: no proof of a camera fault or of beach change)")
    # the range over ALL the camera's days, from the unrounded offsets, with the days that make it (review r1: the
    # Mar fixture's printed 0.17 m was round(0.1647, 3) rounded again, next to first/last offsets 0.159 apart)
    mar = {f"c2 2025-03-0{i}": o for i, o in zip(range(3, 10), [0.066, 0.072, 0.03, -0.01, -0.03, -0.06, -0.0927])}
    tm = dfc.day_trends(mar, LIMIT)
    wm = dfc.day_trend_warning(tm[0], LIMIT, [])[0] if tm else ""
    check(tm and "+0.066 m on 2025-03-03 to -0.093 m on 2025-03-09; over all its days -0.093 m (2025-03-09) to "
                 "+0.072 m (2025-03-04), a range of 0.165 m > the 0.15 m limit" in wm,
          "the trend's range is over all its days, with its lowest and highest day, from unrounded offsets: " + wm[:60])
    flat = dfc.day_trends({"c1 2026-10-0%d" % i: o for i, o in enumerate([0.0, 0.03, 0.05, 0.08, 0.1], 1)}, LIMIT)
    check(flat == [], "a monotonic run within the limit (range 0.10 m) is no warning")
    four = dfc.day_trends({"c1 2026-10-0%d" % i: o for i, o in enumerate([-0.2, -0.1, 0.1, 0.2], 1)}, LIMIT)
    check(four == [], "fewer than 5 testable days: no warning")


def test_trend_rule():
    print("the trend rule: what |rho| it asks for, and how often it fires with no trend")
    need = {n: dfc.trend_rho_needed(n) for n in range(5, 10)}
    check([round(need[n], 3) for n in range(5, 10)] == [1.0, 0.886, 0.821, 0.81, 0.8],
          f"|rho| needed: 5 days 1.00 (strictly monotonic), 6 {need[6]:.3f}, 7 {need[7]:.3f}, 8 {need[8]:.3f}, "
          f"from 9 on > 0.8")
    from itertools import permutations
    smallest = {}
    for n in (7, 8):
        q = np.array(list(permutations(range(n))), dtype=np.int64)
        r = np.unique(np.round(1.0 - 6.0 * ((q - np.arange(n)) ** 2).sum(axis=1) / (n * (n * n - 1.0)), 9))
        smallest[n] = float(r[r > dfc.DAY_TREND_MIN_RHO][0])
    check(abs(need[7] - smallest[7]) < 1e-9 and abs(need[8] - smallest[8]) < 1e-9,
          "from 7 days on that is |rho| > 0.8 itself (the rule asked for): the smallest rho above 0.8 that 7 or 8 "
          "days can give")
    p = np.array(list(permutations(range(9))), dtype=np.int64)
    rho9 = 1.0 - 6.0 * ((p - np.arange(9)) ** 2).sum(axis=1) / (9 * 80.0)
    vals = np.unique(np.round(rho9, 9))
    exact9 = min(v for v in vals if (rho9 >= v - 1e-9).mean() <= dfc.DAY_TREND_ALPHA)
    check(exact9 < dfc.DAY_TREND_MIN_RHO,
          f"9 days, all 362880 orders: the {100 * dfc.DAY_TREND_ALPHA:g} % value is {exact9:.3f}, below 0.8, so "
          f"from 9 days on the 0.8 floor is the stricter (the value falls as days are added)")
    five_swap = {f"c1 2026-10-{i:02d}": o for i, o in enumerate([-0.1, -0.2, 0.0, 0.05, 0.2], 5)}
    six_cycle = {f"c1 2026-10-{i:02d}": o for i, o in enumerate([-0.1, 0.0, -0.2, 0.05, 0.15, 0.2], 5)}
    rising = {f"c1 2026-10-{i:02d}": o for i, o in enumerate([-0.2, -0.1, 0.0, 0.05, 0.1, 0.2], 5)}
    check(dfc.day_trends(five_swap, LIMIT) == [] and dfc.day_trends(six_cycle, LIMIT) == []
          and len(dfc.day_trends(five_swap, LIMIT, alpha=1.0)) == 1
          and len(dfc.day_trends(six_cycle, LIMIT, alpha=1.0)) == 1
          and abs(dfc.day_trends(rising, LIMIT)[0]["chance"] - 1 / 720.0) < 1e-12
          and "which days in random order reach in that direction 1 time in 720 (either way 1 in 360);"
          in dfc.day_trend_warning(dfc.day_trends(rising, LIMIT)[0], LIMIT, [])[0],
          "5 days rising with one pair swapped (rho 0.90, reached by 4.2 % of random orders) and 6 days with "
          "three rotated (rho 0.83, 2.9 %): no warning (|rho| > 0.8 alone warned); strictly rising: WARNING, "
          "chance 1/720")
    # the review's planted trend on the station's data (+0.06 m/day on c1, 29 Sep - 5 Oct): the final offsets
    # it gave, the two low-end days left out; DAY_TREND_ALPHA 0.01 (rho 0.89 needed over 7 days) missed it
    planted = dict(zip(["c1 2026-09-29", "c1 2026-09-30"] + [f"c1 2026-10-0{d}" for d in range(1, 6)],
                       [-0.245, -0.151, -0.130, -0.067, 0.103, -0.010, -0.023]))
    tr = dfc.day_trends(planted, LIMIT, left_out=["c1 2026-09-29", "c1 2026-09-30"])
    check(len(tr) == 1 and tr[0]["spearman_rho"] == 0.857 and tr[0]["without"] is None
          and dfc.day_trends(planted, LIMIT, alpha=0.01) == [],
          f"a 0.35 m trend trimmed at its low end (two days left out), rho +0.86 over 7 days: WARNING "
          f"(chance {tr[0]['chance'] if tr else 0:.4f}; at 1 % it was missed)")
    # one bad day at an end ranks first or last and lifts rho: alone it makes no trend. The camera's
    # other days must still run the same way, |rho| > 0.8 whatever their number; when they do but would
    # not warn on their own, the WARNING says it rests on that day
    days7 = [f"c1 2026-10-{d:02d}" for d in range(5, 12)]
    lone = dict(zip(days7, [0.02, -0.02, 0.0, 0.03, 0.05, 0.04, 0.40]))
    t_all = dfc.day_trends(lone, LIMIT)
    check(len(t_all) == 1 and t_all[0]["spearman_rho"] == 0.857
          and dfc.day_trends(lone, LIMIT, left_out=["c1 2026-10-11"]) == [],
          f"6 days within 0.07 m (rho +0.77 among them) and 11 Oct +0.40 m, left out: rho over all 7 is "
          f"{t_all[0]['spearman_rho'] if t_all else 0:+.2f}, but without the one day left out at the end the "
          f"6 others do not run that way (|rho| not > 0.8): no warning")
    rests = dict(zip(days7, [0.0, -0.02, 0.02, 0.01, 0.04, 0.03, 0.40]))
    tr = dfc.day_trends(rests, LIMIT, left_out=["c1 2026-10-11"])
    w = dfc.day_trend_warning(tr[0], LIMIT, ["c1 2026-10-11"]) if tr else [""]
    check(len(tr) == 1 and tr[0]["without"]["spearman_rho"] == 0.829 and not tr[0]["without"]["alone"]
          and "it rests on c1 2026-10-11, the one day left out, at an end: the other 6 days give rho +0.83" in w[0]
          and "not enough to pass the rule alone (range 0.06 m)" in w[0]
          and "may be the end of the trend or one bad day" in " ".join(" ".join(w[1:]).split())
          and sum("WARNING" in x for x in w) == 1,
          "6 days within 0.06 m that run the same way (rho +0.83) and 11 Oct +0.40 m, left out: WARNING, "
          "saying it rests on 11 Oct (the 6 alone would not warn: range 0.06 m, and 0.83 needs 0.89 over 6)")
    # no trend: independent day biases, two cameras, 7 days, as the NUC's c2 (sd ~0.08 m)
    rng = np.random.RandomState(12)
    n, fired, end_on, end_off = 4000, 0, 0, 0
    days = [f"2026-10-{d:02d}" for d in range(5, 12)]
    tail7 = dfc.spearman_tail(need[7], 7)
    bound = 1 - (1 - 2 * tail7) ** 2
    end5 = 0
    for _ in range(n):
        off = {f"{c} {d}": float(rng.normal(0, 0.08)) for c in ("c1", "c2") for d in days}
        fired += bool(dfc.day_trends(off, LIMIT))
        one = {f"c1 {d}": float(rng.normal(0, 0.04)) for d in days}
        one["c1 2026-10-11"] = 0.40                                     # a bad last day, left out
        end_on += bool(dfc.day_trends(one, LIMIT, left_out=["c1 2026-10-11"]))
        end_off += bool(dfc.day_trends(one, LIMIT))
        five = {f"c1 {d}": float(rng.normal(0, 0.04)) for d in days[:4]}
        five["c1 2026-10-09"] = 0.40                                    # 5 days, the last one bad
        end5 += bool(dfc.day_trends(five, LIMIT, left_out=["c1 2026-10-09"]))
    check(fired / n <= bound,
          f"{n} 7-day windows of two cameras with no trend (day biases sd 0.08 m): WARNING in "
          f"{100.0 * fired / n:.1f} % (at most {100 * bound:.1f} %: each camera {100 * tail7:.1f} % a direction)")

    def near(k, p):                                     # k of n, within 4 binomial sd of p
        return abs(k / n - p) <= 4 * np.sqrt(p * (1 - p) / n)
    p6, p4 = (float((dfc.spearman_null(m) > 0.8 + 1e-9).mean()) for m in (6, 4))
    check(near(end_on, p6) and end_off / n > 2 * end_on / n,
          f"one camera, 6 good days and a bad last day (+0.40 m) left out: WARNING in {100.0 * end_on / n:.1f} % "
          f"(the 6 others run up with rho > 0.8 by chance: {100 * p6:.1f} %); counting that day in the rho, "
          f"{100.0 * end_off / n:.1f} %")
    check(near(end5, p4),
          f"5 days, the last one bad (+0.40 m) and left out: WARNING in {100.0 * end5 / n:.1f} % (the 4 others "
          f"rise strictly by chance, 1 in 24: {100 * p4:.1f} %) -- the price of warning on a 5-day trend "
          f"trimmed at one end, which the 4 days left can never show alone")


def test_trend_end():
    print("a trend trimmed at one end of a camera's days: the other days must run the same way")
    rise5 = [(f"c1 2026-10-0{d}", b, range(80), 3) for d, b in zip(range(1, 6), (-0.35, -0.02, 0.0, 0.03, 0.06))]
    dt, _ = run(rise5)
    tr = dfc.day_trends({k: v[0] for k, v in dt["final"].items()}, LIMIT, left_out=dt["rejected"])
    w = dfc.day_trend_warning(tr[0], LIMIT, dt["rejected"]) if tr else [""]
    check(dt["rejected"] == ["c1 2026-10-01"] and len(tr) == 1 and tr[0]["days"] == 5
          and tr[0]["spearman_rho"] == 1.0 and tr[0]["without"]["days"] == 4 and not tr[0]["without"]["alone"]
          and "it rests on c1 2026-10-01" in w[0] and "too few to pass the rule alone (it asks 5)" in w[0]
          and sum("WARNING" in x for x in w) == 1,
          f"5 days at -0.35, -0.02, 0, +0.03, +0.06 m, the low end left out: rho +1.00 over the 5, range "
          f"{tr[0]['range_m'] if tr else 0:.2f} m: WARNING, resting on 1 Oct (before: none, the 4 others being "
          f"too few to carry it)")
    flat5 = [(f"c1 2026-10-0{d}", b, range(80), 3) for d, b in
             zip(range(1, 7), (-0.30, -0.05, -0.06, -0.03, -0.04, 0.10))]
    dt, _ = run(flat5, sigma=0.0)
    fin = {k: v[0] for k, v in dt["final"].items()}
    rest = [fin[k] for k in sorted(fin) if k not in dt["rejected"]]
    check(dt["rejected"] == ["c1 2026-10-01"] and dfc.day_trends(fin, LIMIT) != []
          and abs(dfc.spearman(np.arange(5), rest) - 0.8) < 1e-12
          and dfc.day_trends(fin, LIMIT, left_out=dt["rejected"]) == [],
          "6 days, the low end (-0.30 m) left out, the other 5 flat but the last (rho +0.80 among them): the 6 "
          "pass the rule (rho +0.89), but the 5 others do not run that way (|rho| not > 0.8): no warning -- the "
          "trend rests on the one left-out end day alone")
    rise6 = [(f"c1 2026-10-0{d}", b, range(80), 3) for d, b in
             zip(range(1, 7), (-0.32, -0.06, -0.02, -0.04, 0.01, 0.05))]
    dt, _ = run(rise6, sigma=0.0)
    tr = dfc.day_trends({k: v[0] for k, v in dt["final"].items()}, LIMIT, left_out=dt["rejected"])
    check(dt["rejected"] == ["c1 2026-10-01"] and len(tr) == 1 and tr[0]["without"]["spearman_rho"] == 0.9
          and not tr[0]["without"]["alone"],
          f"6 days, the low end left out, the other 5 rising (rho +0.90, one pair swapped; range "
          f"{tr[0]['without']['range_m'] if tr else 0:.2f} m): WARNING, resting on 1 Oct")


def test_nuc_7day():
    print("the station's 7-day DEM of 5-11 Oct 2026, replayed")
    bias, days, err = replay_nuc()
    dt, one = run(days, sigma=NUC_SIGMA, seed=3)
    check(err < 0.008, f"the replay gives the NUC's first-pass offsets back within {err * 1000:.1f} mm "
                       f"(c1 {' '.join('%+.3f' % one[k][0] for k in sorted(one) if k[:2] == 'c1')})")
    check(one_pass_rejects(one) == ["c1 2026-10-05", "c2 2026-10-08"],
          "one pass: 5 Oct c1 and 8 Oct c2 rejected, 11 Oct c1 kept, as on the NUC")
    check(dt["rejected"] == ["c1 2026-10-05", "c2 2026-10-08"] and not dt["taken_back"],
          f"iterative: the same two left out ({', '.join(s[0] for s in dt['steps'])}, in that order): "
          f"11 Oct c1 is {dt['final']['c1 2026-10-11'][0]:+.3f} m from the kept days, within the limit")
    tr = dfc.day_trends({k: v[0] for k, v in dt["final"].items()}, LIMIT, left_out=dt["rejected"])
    check([t["camera"] for t in tr] == ["c1"] and tr[0]["spearman_rho"] == 1.0,
          f"so the WARNING carries it: c1 rho {tr[0]['spearman_rho']:+.2f}, range {tr[0]['range_m']:.2f} m "
          f"(5 Oct left out, 11 Oct kept: the trend trimmed at its low end only); none for c2")
    w = tr[0]["without"]
    check(w is not None and w["day"] == "c1 2026-10-05" and w["days"] == 5 and w["spearman_rho"] == 1.0
          and w["alone"],
          "5 Oct is c1's one day left out, at an end: the trend must stand without it, and does (the other "
          "5 days rise strictly, rho +1.00, and pass the rule on their own)")
    check(all(abs(dt["final"][k][0]) <= LIMIT for k in dt["final"] if k not in dt["rejected"]),
          "every kept day within the limit of the kept days")
    # the same window with c1 seen on one day fewer: 5 testable days, 5 Oct trimmed, 4 left -- too few
    # for the rule on their own, so before Oct 14 no warning at all (the review's case)
    for drop in ("c1 2026-10-07", "c1 2026-10-09"):
        dt5, _ = run([d for d in days if d[0] != drop], sigma=NUC_SIGMA, seed=3)
        tr5 = dfc.day_trends({k: v[0] for k, v in dt5["final"].items()}, LIMIT, left_out=dt5["rejected"])
        w5 = dfc.day_trend_warning(tr5[0], LIMIT, dt5["rejected"]) if tr5 else [""]
        c1 = [round(dt5["final"][k][0], 3) for k in sorted(dt5["final"]) if k[:2] == "c1"]
        check(dt5["rejected"] == ["c1 2026-10-05", "c2 2026-10-08"] and [t["camera"] for t in tr5] == ["c1"]
              and tr5[0]["days"] == 5 and tr5[0]["spearman_rho"] == 1.0 and tr5[0]["without"]["days"] == 4
              and not tr5[0]["without"]["alone"] and "it rests on c1 2026-10-05" in w5[0]
              and "too few to pass the rule alone (it asks 5)" in w5[0],
              f"without {drop[3:]}: c1's 5 final offsets {' '.join('%+.3f' % o for o in c1)} rise strictly, 5 Oct "
              f"left out: WARNING, resting on 5 Oct (the other 4 days rise strictly too, too few alone)")


def test_too_few():
    print("fewer than 3 testable days: the test cannot tell which day is off, none left out")
    days = [("c2 2025-03-06", 0.0, range(100), 3), ("c2 2025-03-07", 0.3, range(100), 3)]
    dt, one = run(days)
    check(dt["too_few"] and dt["rejected"] == [] and dt["steps"] == [], "two days 0.3 m apart: both kept")
    lines = dfc.day_test_report(dt, LIMIT)
    check(any("kept: fewer than 3 testable days" in x for x in lines) and "every day kept" in lines[-1]
          and not any("REJECTED" in x for x in lines) and dfc.day_test_warning(dt, LIMIT) == ("", []),
          "... and the lines say so (as test_survey_products expects), no WARNING (as before)")
    # untestable days (too few comparable samples) are kept and do not count
    days = [("c2 2025-03-0%d" % d, 0.0, range(100), 3) for d in (1, 2, 3, 4)]
    days.append(("c2 2025-03-05", 0.5, range(5), 3))                     # 15 samples < 30
    dt, one = run(days)
    check(np.isnan(dt["first"]["c2 2025-03-05"][0]) and dt["rejected"] == [] and not dt["untested"],
          "a day with fewer than 30 comparable samples: kept, not testable (and not 'kept untested')")


def scenario_fringe(checked=True):
    """The low beach (cells 0-19) is reached by two c2 days only: 20 Sep (good, 3 frames per cell)
    and 22 Sep (0.57 m LOW, 8 frames). Nine good days cover cells 20-99. With checked, 20 Sep also
    reaches cells 20-49, where the nine days check it; without, 20 and 22 Sep only see each other."""
    days = [(f"c2 2026-09-{d}", 0.0, range(20, 100), 1) for d in range(11, 20)]
    days.append(("c2 2026-09-20", 0.0, range(50 if checked else 20), 3))
    days.append(("c2 2026-09-22", -0.57, range(20), 8))
    return days


def test_fringe():
    print("a fringe day whose cells only one other day reaches")
    dt, one = run(scenario_fringe())
    check(one_pass_rejects(one) == ["c2 2026-09-22"],
          f"one pass: 22 Sep rejected ({one['c2 2026-09-22'][0]:+.3f} m), 20 Sep kept "
          f"({one['c2 2026-09-20'][0]:+.3f} m)")
    loo = dfc.LeaveOneOut(*pairs(scenario_fringe()), MIN_FRAMES)
    j = loo.judge("c2 2026-09-22", np.ones(len(loo.names), bool), LIMIT)
    check(j["judged"] and j["share"] == 1.0 and not j["alone"],
          "judged: its DEM is 20 Sep alone, but the nine other days check 20 Sep where they meet it "
          "(within the limit): 20 Sep vouches")
    check(dt["rejected"] == ["c2 2026-09-22"] and not dt["not_judged"] and not dt["untested"]
          and dfc.day_test_warning(dt, LIMIT) == ("", []),
          "iterative: 22 Sep left out, every other day kept, no WARNING (the guard of d28b5ec kept it, "
          "silently)")
    dt, one = run(scenario_fringe(checked=False))
    check(one_pass_rejects(one) == ["c2 2026-09-20", "c2 2026-09-22"],
          "20 and 22 Sep only see each other: one pass rejected both")
    check(dt["rejected"] == ["c2 2026-09-20", "c2 2026-09-22"]
          and [r["day"] for r in dt["not_judged"]] == ["c2 2026-09-20", "c2 2026-09-22"]
          and all(r["left_out"] for r in dt["not_judged"]),
          "iterative: neither can be judged -- the first pass's verdict, both left out")
    head, lines = dfc.day_test_warning(dt, LIMIT)
    check(head.startswith("WARNING: the day test could not judge every camera-day -- left out without being "
                          "judged: c2 2026-09-20, c2 2026-09-22")
          and sum("WARNING" in x for x in lines) == 1 and all(x.startswith("   ") for x in lines[1:]),
          "one WARNING line names them")
    rep = dfc.day_test_report(dt, LIMIT, 1)
    check(all("REJECTED (not judged: over the limit in the first pass)" in x for x in rep
              if x.startswith("   c2 2026-09-2")),
          "their lines say why")


def test_cannot_judge():
    print("days the test cannot judge keep the first pass's verdict, with a WARNING")
    # c1's two days only see each other (c2 does not share their cells)
    days = [("c1 2026-10-0%d" % d, 0.3 * (d - 1), range(100), 3) for d in (1, 2)]
    days += [("c2 2026-10-0%d" % d, 0.0, range(100, 200), 3) for d in (1, 2, 3, 4, 5)]
    dt, one = run(days)
    check(one_pass_rejects(one) == ["c1 2026-10-01", "c1 2026-10-02"],
          "c1's two days only see each other: one pass rejected BOTH (7 testable days in all)")
    check(dt["rejected"] == ["c1 2026-10-01", "c1 2026-10-02"] and len(dt["not_judged"]) == 2
          and "alone, which the other kept days do not confirm" in dt["not_judged"][0]["why"],
          "iterative: neither can be judged (each one's DEM is the other day alone): both left out, as one pass")
    # the same, with 30 % of c1's cells shared with c2: the overlap tells which c1 day is off
    days = [("c1 2026-10-05", 0.0, range(100), 3), ("c1 2026-10-06", 1.0, range(100), 3)]
    days += [(f"c2 2026-10-0{d}", 0.0, range(70, 300), 3) for d in range(3, 10)]
    dt, one = run(days)
    check(one_pass_rejects(one) == ["c1 2026-10-05", "c1 2026-10-06"] and dt["rejected"] == ["c1 2026-10-06"]
          and not dt["not_judged"],
          "with 30 % of their cells shared with c2: 6 Oct (+1.0 m) left out, 5 Oct kept (one pass: both out); "
          "5 Oct, checked by c2 where they meet, vouches")
    # three days, two off either way: one goes, then two are left that disagree
    days = [("c2 2025-03-06", 0.0, range(100), 3), ("c2 2025-03-07", 0.4, range(100), 3),
            ("c2 2025-03-08", -0.4, range(100), 3)]
    dt, one = run(days)
    nj = {r["day"]: r for r in dt["not_judged"]}
    check(dt["rejected"] == one_pass_rejects(one) == ["c2 2025-03-07", "c2 2025-03-08"]
          and len(dt["steps"]) == 2 and sorted(nj) == ["c2 2025-03-06", dt["steps"][1][0]]
          and not nj["c2 2025-03-06"]["left_out"] and "only 2 testable camera-day(s) kept" in nj["c2 2025-03-06"]["why"],
          f"three days 0, +0.4, -0.4 m: {dt['steps'][0][0][3:]} left out, then two that disagree: the first "
          f"pass's verdict, 6 Mar kept and {dt['steps'][1][0][3:]} left out, as one pass")
    head = dfc.day_test_warning(dt, LIMIT)[0]
    rec = {x["day"]: x for x in dfc.day_test_record(dt, LIMIT, [])["day_test"]["not_judged"]}
    check("kept untested: c2 2025-03-06" in head and "kept over the limit" not in head
          and dt["untested"] == ["c2 2025-03-06"] and not dt["over_kept"]
          and rec["c2 2025-03-06"]["verdict"] == "kept untested" and rec["c2 2025-03-06"]["final_offset_m"] is None
          and rec["c2 2025-03-06"]["first_pass_verdict"] == "kept",
          f"6 Mar ends the only day kept, so untested: the WARNING says so, not 'kept over the limit' "
          f"({head[9:110]}...)")
    # every day disagrees: every day rejected (one pass too), so the command line exits 4 as before
    days = [(f"c1 2026-10-0{d}", 0.4 * (d - 1), range(100), 3) for d in range(1, 5)]
    dt, one = run(days)
    check(dt["rejected"] == one_pass_rejects(one) == sorted(k for k, _, _, _ in days)
          and [r["day"] for r in dt["not_judged"]] == ["c1 2026-10-01", "c1 2026-10-02"],
          "4 days 0, +0.4, +0.8, +1.2 m: 4 and 3 Oct left out, then 1 and 2 Oct cannot be judged and were "
          "over the limit in the first pass: every day rejected, as one pass")
    rep = dfc.day_test_report(dt, LIMIT, 1200)
    check(sum("REJECTED" in x for x in rep[1:5]) == 4, "every day's line says REJECTED (survey_products "
                                                         "reads 'rejected every camera-day' from that)")


def scenario_taken_back_unjudged():
    """Two c2 days 0.4 m high -- 1 Oct on cells 50-99 (two frames a cell), 2 Oct on 20-59 (three) --
    and two good ones: 3 Oct on 0-49 (one frame), 4 Oct on 30-79 (three). 2 Oct goes first; 3 Oct then
    has too few comparable samples left, so only 1 and 4 Oct are testable and the test cannot tell which
    is off: both keep the first pass's verdict, over the limit (4 Oct read -0.40 m against a DEM the
    high days made), and are left out. Against the kept day (3 Oct) 4 Oct is within the limit: taken
    back. (The review's seed 204 of fuzz2, in small.)"""
    return [("c2 2026-10-01", 0.4, range(50, 100), 2), ("c2 2026-10-02", 0.4, range(20, 60), 3),
            ("c2 2026-10-03", 0.0, range(0, 50), 1), ("c2 2026-10-04", 0.0, range(30, 80), 3)]


def scenario_kept_unjudged():
    """1 Oct 0.2 m low on cells 10-49 (three frames a cell), 2 Oct 0.2 m high on 20-69 (three), and two
    good days: 3 Oct on 30-49 (one frame; under 30 comparable samples, not testable) and 4 Oct on 10-49
    (two). 2 Oct goes first; then 1 and 4 Oct only see each other: 1 Oct (over in the first pass) is left
    out, 4 Oct (within it) kept -- and in the end 4 Oct is within the limit of the kept days. (The
    review's seed 20 of fuzz2, in small.)"""
    return [("c2 2026-10-01", -0.2, range(10, 50), 3), ("c2 2026-10-02", 0.2, range(20, 70), 3),
            ("c2 2026-10-03", 0.0, range(30, 50), 1), ("c2 2026-10-04", 0.0, range(10, 50), 2)]


def test_unjudged_fate():
    print("days the test could not judge: the WARNING and <stem>_info.json say where they ended up")
    for sigma in (0.02, 0.0):
        pc, pv, pk = pairs(scenario_taken_back_unjudged(), sigma=sigma)
        dt = dfc.day_consistency(pc, pv, pk, LIMIT, 1)
        nj = {r["day"]: r for r in dt["not_judged"]}
        check(dt["rejected"] == ["c2 2026-10-01", "c2 2026-10-02"] and dt["taken_back"] == ["c2 2026-10-04"]
              and sorted(nj) == ["c2 2026-10-01", "c2 2026-10-04"] and all(r["left_out"] for r in nj.values())
              and dt["untested"] == ["c2 2026-10-03"] and not check_invariants(dt),
              f"noise {sigma} m: 1 and 4 Oct left out unjudged, 4 Oct taken back "
              f"({dt['final']['c2 2026-10-04'][0]:+.3f} m from the kept days), 3 Oct kept untested")
    head, lines = dfc.day_test_warning(dt, LIMIT)
    check(warning_names(head) == {"left out without being judged": ["c2 2026-10-01"],
                                  "kept untested": ["c2 2026-10-03"]},
          f"the WARNING names 1 Oct left out and 3 Oct untested, not 4 Oct, which is in the DEM ({head[9:]})")
    rec = {x["day"]: x for x in dfc.day_test_record(dt, LIMIT, [])["day_test"]["not_judged"]}
    r4 = rec["c2 2026-10-04"]
    check(r4["verdict"] == "taken back" and r4["first_pass_verdict"] == "left out"
          and r4["final_offset_m"] is not None and rec["c2 2026-10-01"]["verdict"] == "left out",
          "<stem>_info.json: day_test.not_judged gives 4 Oct verdict 'taken back' (first pass: 'left out') "
          "and its final offset")
    row = [x for x in dfc.day_test_report(dt, LIMIT, 1) if x.startswith("   c2 2026-10-04")][0]
    check("taken back" in row and "REJECTED" not in row and "WARNING" not in row, f"its line: '{row.strip()}'")
    pc, pv, pk = pairs(scenario_kept_unjudged(), sigma=0.02)
    dt = dfc.day_consistency(pc, pv, pk, LIMIT, 1)
    nj = {r["day"]: r for r in dt["not_judged"]}
    o4 = dt["final"]["c2 2026-10-04"][0]
    check(dt["rejected"] == ["c2 2026-10-01", "c2 2026-10-02"] and "c2 2026-10-04" in nj
          and not nj["c2 2026-10-04"]["left_out"] and abs(o4) <= LIMIT and not dt["over_kept"]
          and not dt["untested"] and not check_invariants(dt),
          f"4 Oct, kept unjudged, ends {o4:+.3f} m from the kept days: within the limit")
    head, lines = dfc.day_test_warning(dt, LIMIT)
    rec = {x["day"]: x for x in dfc.day_test_record(dt, LIMIT, [])["day_test"]["not_judged"]}
    rep = dfc.day_test_report(dt, LIMIT, 1)
    check(warning_names(head) == {"left out without being judged": ["c2 2026-10-01"]}
          and rec["c2 2026-10-04"]["verdict"] == "kept within the limit"
          and not any("WARNING" in x for x in rep),
          "the WARNING names 1 Oct only, not 'kept over the limit: 4 Oct'; its record says 'kept within the "
          "limit'; no report line holds the word WARNING")
    # left out unjudged, but the kept days judge it in the end: not named either. 1 Oct (good, two
    # frames a cell) and 4 Oct (0.2 m low, three) share cells 10-49, where only 1 Oct also meets 2 Oct
    # (cells 40-49); 2 and 3 Oct (good) cover 40-79. 1 and 4 Oct cannot be judged (each one's DEM is
    # the other alone) and were over the limit in the first pass: both left out. Against the kept days
    # 4 Oct is 0.2 m low where 2 Oct, which 3 Oct confirms, meets it: judged in the end. 1 Oct has too
    # few comparable samples left: still left out unjudged.
    days = [("c2 2026-10-01", 0.0, range(10, 50), 2), ("c2 2026-10-02", 0.0, range(40, 80), 3),
            ("c2 2026-10-03", 0.0, range(60, 80), 2), ("c2 2026-10-04", -0.2, range(10, 50), 3)]
    for sigma in (0.02, 0.0):
        pc, pv, pk = pairs(days, sigma=sigma)
        dt = dfc.day_consistency(pc, pv, pk, LIMIT, 1)
        nj = {r["day"]: r for r in dt["not_judged"]}
        head = dfc.day_test_warning(dt, LIMIT)[0]
        t = dfc.day_test_record(dt, LIMIT, [])
        how = {x["day"]: x["how"] for x in t["day_test"]["rejected"]}
        row = [x for x in dfc.day_test_report(dt, LIMIT, 1) if x.startswith("   c2 2026-10-04")][0]
        check(dt["rejected"] == ["c2 2026-10-01", "c2 2026-10-04"] and sorted(nj) == dt["rejected"]
              and nj["c2 2026-10-04"]["judged_in_end"] and not nj["c2 2026-10-01"]["judged_in_end"]
              and warning_names(head) == {"left out without being judged": ["c2 2026-10-01"]}
              and how == {"c2 2026-10-01": "not judged", "c2 2026-10-04": "left out"}
              and t["camera_days_rejected_offsets_m"]["c2 2026-10-01"] is not None
              and "the kept days judge it over the limit in the end" in row and not check_invariants(dt),
              f"noise {sigma} m: 4 Oct, left out unjudged, is judged against the kept days in the end "
              f"({dt['final']['c2 2026-10-04'][0]:+.3f} m): the WARNING names 1 Oct only; 1 Oct, with no "
              f"kept day to compare, still gets an offset in camera_days_rejected_offsets_m (the one it was "
              f"left out at, {t['camera_days_rejected_offsets_m']['c2 2026-10-01']:+.3f} m)")


def test_set_aside_retried():
    print("a day kept unjudged is tried again later")
    # 4 Oct (0.4 m low, cells 20-69) goes first. 6 Oct (0.2 m low, 40-69) then reads -0.21 m, but its
    # DEM is 2 Oct alone (50-89) in most of its cells, unconfirmed: set aside, kept (within the limit in
    # the first pass, where 4 Oct made its DEM); 2 Oct, as unjudged, is left out (over the limit in the
    # first pass). Without 2 Oct, 6 Oct's DEM is 1 and 5 Oct (good, cells 0-49): judged and left out --
    # d28b5ec/86985d3 kept it over the limit, 'without being judged', for good.
    days = [("c2 2026-10-01", 0.0, range(10, 50), 2), ("c2 2026-10-02", 0.0, range(50, 90), 3),
            ("c2 2026-10-03", -0.2, range(20, 40), 1), ("c2 2026-10-04", -0.4, range(20, 70), 3),
            ("c2 2026-10-05", 0.0, range(0, 50), 3), ("c2 2026-10-06", -0.2, range(40, 70), 3)]
    for sigma in (0.02, 0.0):
        pc, pv, pk = pairs(days, sigma=sigma)
        dt = dfc.day_consistency(pc, pv, pk, LIMIT, 1)
        order = " ".join(x.strip() for x in dfc.day_test_report(dt, LIMIT, 1)).split("order left out:")[1]
        check(dt["judged_later"] == ["c2 2026-10-06"] and "c2 2026-10-06" in dt["rejected"]
              and [s[0] for s in dt["steps"]] == ["c2 2026-10-04", "c2 2026-10-02", "c2 2026-10-06"]
              and [r["day"] for r in dt["not_judged"]] == ["c2 2026-10-02"] and not dt["over_kept"]
              and "c2 2026-10-06 (-0.2" in order and "set aside unjudged at an earlier step" in order
              and warning_names(dfc.day_test_warning(dt, LIMIT)[0]) == {
                  "left out without being judged": ["c2 2026-10-02"]} and not check_invariants(dt),
              f"noise {sigma} m: 6 Oct, set aside, is judged once 2 Oct is out, and left out "
              f"({dt['final']['c2 2026-10-06'][0]:+.3f} m); no day is kept over the limit")


def scenario_survey_minority():
    """A week (20-26 Oct) whose three days around the survey date (22-24 Oct) read 0.25-0.30 m higher
    than the four others, all on one stretch of beach: in the test without an anchor the four
    outweigh them, and the survey's own days are the ones left out."""
    bias = {"20": 0.0, "21": 0.0, "22": 0.25, "23": 0.30, "24": 0.28, "25": 0.0, "26": 0.0}
    return [(f"c2 2024-10-{d}", b, range(150), 3) for d, b in bias.items()]


def test_anchor():
    print("survey-date anchor (--anchor-date): the days kept are those around the survey")
    pc, pv, pk = pairs(scenario_survey_minority(), sigma=0.05)
    free = dfc.day_consistency(pc, pv, pk, LIMIT, MIN_FRAMES)
    around = {"c2 2024-10-22", "c2 2024-10-23", "c2 2024-10-24"}
    check(around <= set(free["rejected"]),
          f"without the anchor the days around the survey are the ones left out ({free['rejected']})")
    dt = dfc.day_consistency(pc, pv, pk, LIMIT, MIN_FRAMES, anchor="2024-10-23")
    check(dt["anchored"] == ["c2 2024-10-23"], f"the survey's own day is the anchor ({dt['anchored']})")
    check(dt["anchor_reference"] == sorted(around),
          f"the days within 1 day of it are the reference ({dt['anchor_reference']})")
    check(not (around & set(dt["rejected"])), f"the days around the survey are kept (left out {dt['rejected']})")
    check(set(dt["rejected"]) == {"c2 2024-10-20", "c2 2024-10-21", "c2 2024-10-25", "c2 2024-10-26"}
          and set(dt["anchor_left"]) == set(dt["rejected"]),
          f"the others, off the days around the survey, are left out against them ({dt['anchor_left']})")
    lines = dfc.day_test_report(dt, LIMIT, 0)
    check(any("anchored on the survey date 2024-10-23" in ln for ln in lines)
          and any("nearest the survey date: kept" in ln for ln in lines)
          and any("off the days around the survey date" in ln for ln in lines)
          and not any("WARNING" in ln for ln in lines),
          "the log says what was anchored, why the anchor day is kept and why the others went")
    rec = dfc.day_test_record(dt, LIMIT, [])["day_test"]
    check(rec["anchor_date"] == "2024-10-23" and rec["anchored_days"] == ["c2 2024-10-23"]
          and all(r["how"] == "off the days around the survey date" for r in rec["rejected"]),
          "<stem>_info.json keeps the anchor, and why each day went")
    far = dfc.day_consistency(pc, pv, pk, LIMIT, MIN_FRAMES, anchor="2024-11-15")
    check(far["anchored"] == [] and far["rejected"] == free["rejected"],
          "no testable day within 1 day of the date: nothing anchored, the test as without it")


def test_deterministic():
    print("deterministic: sample order and ties")
    pc, pv, pk = pairs(scenario_masking(), seed=2)
    a = dfc.day_consistency(pc, pv, pk, LIMIT, MIN_FRAMES)
    p = np.random.RandomState(9).permutation(len(pc))
    b = dfc.day_consistency(pc[p], pv[p], pk[p], LIMIT, MIN_FRAMES)
    check(a == b, "the same result from shuffled samples (every offset, step and take-back)")
    days = [(f"c2 2025-03-0{d}", 0.0, range(100), 3) for d in range(1, 6)]
    days += [("c2 2025-03-07", 0.5, range(100), 3), ("c2 2025-03-06", 0.5, range(100, 200), 3)]
    days += [(f"c2 2025-03-0{d}", 0.0, range(100, 200), 3) for d in range(1, 6)]
    dt, one = run(days, sigma=0.0)
    check(dt["steps"][0][1] == dt["steps"][1][1] and [s[0] for s in dt["steps"]] == ["c2 2025-03-06", "c2 2025-03-07"],
          "two days exactly as far off: the first by name goes first")


def test_cli(d):
    print("dem_from_contours.py --max-day-offset: the log block and <stem>_info.json")
    import survey_products as sp
    g = Path(d) / "ground.csv"
    rows = ["source_file,camera,capture_time_utc,capture_epoch,easting_utm19,northing_utm19,"
            "tide_elevation_navd88,beach_elevation_navd88"]
    rng = np.random.RandomState(4)
    plan = {}
    for i, day in enumerate(["2026-10-05", "2026-10-06", "2026-10-07", "2026-10-08", "2026-10-09",
                             "2026-10-10", "2026-10-11"]):
        plan[("c1", day)] = -0.30 + 0.10 * i                   # c1: a steady rise, -0.30 to +0.30 m
        plan[("c2", day)] = 0.25 if day == "2026-10-08" else 0.0  # c2: one high day
    for (cam, day), bias in sorted(plan.items()):
        x0 = 420000.0 if cam == "c1" else 420040.0
        for h in (14, 15, 16):
            t = datetime.fromisoformat(f"{day}T{h:02d}:00:00+00:00")
            f = f"{int(t.timestamp())}.{cam}.timex"
            for i in range(200):
                e, n = x0 + 2.0 * (i % 20) + 1.0, 4638000.0 + 2.0 * (i // 20) + 1.0
                z = 1.0 - 0.02 * (e - 420000.0) + bias + rng.normal(0, 0.02)
                rows.append(f"{f},{cam},{t:%Y-%m-%dT%H:%M:%S}+00:00,{int(t.timestamp())},{e:.2f},{n:.2f},"
                            f"{z:.4f},{z:.4f}")
    g.write_text("\n".join(rows) + "\n")
    stem = Path(d) / "dem"
    r = subprocess.run([sys.executable, str(HERE / "dem_from_contours.py"), str(g), str(stem),
                        "--max-day-offset", str(LIMIT), "--no-plot"], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
    out = r.stdout
    (Path(d) / "dem.log").write_text(out)
    check(r.returncode == 0 and Path(f"{stem}_dem.asc").exists(), f"built (exit {r.returncode})")
    block = out[out.index("Day consistency"):].splitlines()
    day_lines = [ln for ln in block[1:15] if re.match(r"^\s+(\S+ \d{4}-\d\d-\d\d)\s", ln)]
    check(len(day_lines) == 14 and all(re.match(r"^\s+(\S+ \d{4}-\d\d-\d\d)\s", ln) for ln in block[1:15])
          and block[15].strip().startswith("left out"),
          "one line per camera-day right after the header, then 'left out' (the block survey_products.py reads)")
    info = json.loads(Path(f"{stem}_info.json").read_text())
    f = info["filters"]
    rej = f["camera_days_rejected"]
    check(isinstance(rej, list) and rej == sorted(rej) and all(isinstance(x, str) for x in rej)
          and "c2 2026-10-08" in rej, f"camera_days_rejected: a sorted list of names, as before ({', '.join(rej)})")
    check(sorted(f["camera_days_rejected_offsets_m"]) == rej
          and abs(f["camera_days_rejected_offsets_m"]["c2 2026-10-08"] - 0.25) < 0.03,
          f"camera_days_rejected_offsets_m names each with its offset (c2 8 Oct "
          f"{f['camera_days_rejected_offsets_m']['c2 2026-10-08']:+.3f} m)")
    t = f["day_test"]
    check(len(t["first_pass_offsets_m"]) == 14 and len(t["final_offsets_m"]) == 14 and t["limit_m"] == LIMIT
          and t["left_out_order"] and set(t["left_out_order"]) >= set(rej)
          and [w["camera"] for w in t["trend_warnings"]] == ["c1"],
          "day_test: first-pass and final offsets of all 14 days, the order left out, the c1 trend")
    rej_lines = [ln for ln in day_lines if "REJECTED" in ln]
    check(len(rej_lines) == len(rej) and all(any(k in ln for ln in rej_lines) for k in rej),
          "the lines marked REJECTED are exactly the days left out")
    check(sum("WARNING" in ln for ln in out.splitlines()) == 1
          and "c1's camera-day offsets run one way with time over 7 days" in out,
          "one line holds the word WARNING (the station status counts them), for c1")
    check([x["day"] for x in t["rejected"]] == rej
          and all(x["how"] in ("left out", "not judged", "refused take-back") and x["offset_at_step_m"] is not None
                  for x in t["rejected"])
          and t["warning"] is None and "|rho| needed by days without ties 5: 1.00" in t["trend_rule"],
          "day_test.rejected: each day left out, how and its offsets; the trend rule in words")
    check("rejected every camera-day" not in sp.dem_failure_reason(Path(d) / "dem.log"),
          "survey_products.dem_failure_reason reads the block (no 'every camera-day rejected')")
    # the reader of the info file, as before
    import dem_figure
    _, _, left = dem_figure.title_lines(info, {"filled": 1, "blanked": 0, "no_value": 0}, "TEST", 0.5, 3, 2.0)
    check(any(f"{len(rej)} camera-days off the rest by > {LIMIT:g} m" in x for x in left),
          "the page (dem_figure.title_lines) still counts them")
    h = subprocess.run([sys.executable, str(HERE / "dem_from_contours.py"), "--help"], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d)).stdout
    h = " ".join(h.split())
    check("|rho| needed by days without ties 5: 1.00, 6: 0.89, 7: 0.82, 8: 0.81" in h
          and "when the one day left out of a camera's days is its first or last, the other days must run the "
              "same way, |rho| > 0.8 (any number of days)" in h,
          "--help states the trend rule as applied (the |rho| it asks by number of days, the lone end day)")
    # a camera on 5 days, rising, its low end left out: the WARNING says the trend rests on that day,
    # and <stem>_info.json records the other days (without: their rho, alone False)
    rows5 = rows[:1]
    for i, (day, bias) in enumerate(zip(["2026-10-05", "2026-10-06", "2026-10-07", "2026-10-09", "2026-10-11"],
                                        (-0.35, -0.02, 0.0, 0.03, 0.06))):
        for h in (14, 15, 16):
            t5 = datetime.fromisoformat(f"{day}T{h:02d}:00:00+00:00")
            f5 = f"{int(t5.timestamp())}.c1.timex"
            for j in range(200):
                e, n = 420000.0 + 2.0 * (j % 20) + 1.0, 4638000.0 + 2.0 * (j // 20) + 1.0
                z = 1.0 - 0.02 * (e - 420000.0) + bias + rng.normal(0, 0.02)
                rows5.append(f"{f5},c1,{t5:%Y-%m-%dT%H:%M:%S}+00:00,{int(t5.timestamp())},{e:.2f},{n:.2f},"
                             f"{z:.4f},{z:.4f}")
    (Path(d) / "ground5.csv").write_text("\n".join(rows5) + "\n")
    r = subprocess.run([sys.executable, str(HERE / "dem_from_contours.py"), str(Path(d) / "ground5.csv"),
                        str(Path(d) / "five"), "--max-day-offset", str(LIMIT), "--no-plot"], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
    (Path(d) / "five.log").write_text(r.stdout)
    warn = [ln for ln in r.stdout.splitlines() if "WARNING" in ln]
    t5 = json.loads((Path(d) / "five_info.json").read_text())["filters"]["day_test"]
    tw = t5["trend_warnings"]
    check(r.returncode == 0 and len(warn) == 1 and "c1's camera-day offsets run one way with time over 5 days" in warn[0]
          and "it rests on c1 2026-10-05, the one day left out, at an end: the other 4 days give rho +1.00" in warn[0]
          and len(tw) == 1 and tw[0]["without"]["day"] == "c1 2026-10-05" and tw[0]["without"]["days"] == 4
          and tw[0]["without"]["alone"] is False and tw[0]["without"]["spearman_rho"] == 1.0
          and t5["rejected"][0]["day"] == "c1 2026-10-05",
          "5 days rising, 5 Oct (-0.35 m) left out: one WARNING line, resting on 5 Oct; info.json "
          "trend_warnings[].without gives the other 4 days (rho +1.00, alone false)")
    # the cron's DEMs pass --trend-note: the same trend said as a NOTE, which station_status.py does not count
    rn = subprocess.run([sys.executable, str(HERE / "dem_from_contours.py"), str(Path(d) / "ground5.csv"),
                         str(Path(d) / "five_note"), "--max-day-offset", str(LIMIT), "--no-plot", "--trend-note"],
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
    notes = [ln for ln in rn.stdout.splitlines() if ln.startswith("NOTE: ")]
    tn = json.loads((Path(d) / "five_note_info.json").read_text())["filters"]["day_test"]
    check(rn.returncode == 0 and not any("WARNING" in ln for ln in rn.stdout.splitlines())
          and len(notes) == 1 and notes[0][len("NOTE: "):] == warn[0][len("WARNING: "):]
          and tn["trend_warnings"] == tw,
          "--trend-note (the cron): the same trend as one NOTE line, no WARNING line; info.json's "
          "trend_warnings identical")
    # no test: the info file says so
    r = subprocess.run([sys.executable, str(HERE / "dem_from_contours.py"), str(g), str(Path(d) / "plain"),
                        "--no-plot"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True,
                       cwd=str(d))
    f = json.loads((Path(d) / "plain_info.json").read_text())["filters"]
    check(r.returncode == 0 and f["day_test"] is None and f["camera_days_rejected"] == []
          and f["camera_days_rejected_offsets_m"] == {} and "Day consistency" not in r.stdout,
          "without --max-day-offset: no test, day_test None")


def check_invariants(dt, limit=LIMIT):
    """What every result of day_consistency must hold. -> '' or what failed."""
    rej, back = set(dt["rejected"]), set(dt["taken_back"])
    steps = [x[0] for x in dt["steps"]]
    nj = {r["day"]: r for r in dt["not_judged"]}
    refused = {r[0] for r in dt["refused"]}
    fin, first = dt["final"], dt["first"]
    if dt["too_few"]:                       # fewer than 3 testable days: every day kept, as before
        return "" if not (steps or rej or nj) and fin == first else "too few: every day kept"
    if not (rej <= set(steps) and not (rej & back) and set(steps) == rej | back and len(steps) == len(set(steps))):
        return "bookkeeping: left out = rejected + taken back, each once"
    for k, r in nj.items():
        if r["left_out"] != (np.isfinite(first[k][0]) and abs(first[k][0]) > limit):
            return f"{k}: a day not judged keeps the first pass's verdict"
        if r["left_out"] and k not in steps:
            return f"{k}: not judged, over in the first pass, not left out"
    for k in fin:
        o = fin[k][0]
        if k not in rej and np.isfinite(o) and abs(o) > limit and k not in nj:
            return f"{k}: kept over the limit, but judged"
        if k in rej and np.isfinite(o) and abs(o) <= limit and k not in refused:
            return f"{k}: rejected within the limit of the kept days, not a refused take-back"
    if sorted(dt["over_kept"]) != sorted(k for k in fin if k not in rej and np.isfinite(fin[k][0])
                                         and abs(fin[k][0]) > limit):
        return "over_kept"
    if not set(dt["over_kept"]) <= set(nj):
        return "a day kept over the limit was judged"
    for k, r in nj.items():
        if r["judged_in_end"] and not (k in rej and np.isfinite(fin[k][0]) and abs(fin[k][0]) > limit):
            return f"{k}: judged in the end, but not left out over the limit of the kept days"
    if sorted(dt["untested"]) != sorted(k for k in fin if k not in rej and np.isfinite(first[k][0])
                                        and not np.isfinite(fin[k][0])):
        return "untested: every kept day testable in the first pass and not at the end"
    # the WARNING describes the final DEM: what it names left out is left out, what it names kept
    # over the limit is kept over it, what it names untested is untested -- and it names all of them
    head, lines = dfc.day_test_warning(dt, limit)
    named = warning_names(head)
    want = {"left out without being judged": sorted(k for k in nj if k in rej and not nj[k]["judged_in_end"]),
            "kept over the limit without being judged": sorted(dt["over_kept"]),
            "kept untested": sorted(dt["untested"])}
    if named != {k: v for k, v in want.items() if v}:
        return f"the WARNING names {named}, the final DEM has {want}"
    if sum("WARNING" in x for x in lines) > 1:
        return "one WARNING line"
    rep = dfc.day_test_report(dt, limit, 1 if rej else None)
    if sum("REJECTED" in x for x in rep) != len(rej) or not rep[0].startswith("Day consistency"):
        return "the report's REJECTED lines"
    if any("WARNING" in x for x in rep):
        return "a report line holds the word WARNING (the station status counts it)"
    trends = dfc.day_trends({k: v[0] for k, v in fin.items()}, limit, left_out=dt["rejected"])
    rec = dfc.day_test_record(dt, limit, trends)
    json.dumps(rec, allow_nan=False)
    t = rec["day_test"]
    for x in t["not_judged"]:
        k = x["day"]
        fate = ("left out" if k in rej else "taken back" if k in back else "kept over the limit"
                if k in dt["over_kept"] else "kept untested" if k in dt["untested"] else "kept within the limit")
        fo = None if not np.isfinite(fin[k][0]) else round(fin[k][0], 4)
        if x["verdict"] != fate or x["final_offset_m"] != fo:
            return f"{k}: the record's verdict {x['verdict']!r}, the final DEM's {fate!r}"
    if (sorted(rec["camera_days_rejected_offsets_m"]) != sorted(rej)
            or any(v is None for v in rec["camera_days_rejected_offsets_m"].values())
            or [x["day"] for x in t["rejected"]] != dt["rejected"]
            or any(x["how"] != ("not judged" if x["day"] in nj and not nj[x["day"]]["judged_in_end"] else
                                "left out" if x["day"] in nj else "refused take-back" if x["day"] in refused
                                else "left out") for x in t["rejected"])):
        return "the record's rejected days (each with an offset and how)"
    if t["warning"] != (head or None):
        return "the record's warning is the printed one"
    return ""


def warning_names(head):
    """{what: [days]} from the day test's WARNING line."""
    if not head:
        return {}
    body = head.split(" -- ", 1)[1].rstrip(".")
    out = {}
    for part in body.split("; "):
        what, days = part.split(": ", 1)
        out[what] = sorted(days.split(", "))
    return out


def random_days(seed):
    """A random two-camera window: up to 8 days a camera, a quarter of them off by N(0, 0.3) m, on
    random stretches of cells (some that only one day reaches), 1-5 frames a cell."""
    rng = np.random.RandomState(seed)
    days = []
    for cam, c0 in (("c1", 0), ("c2", int(rng.randint(50, 200)))):
        for d in range(int(rng.randint(1, 9))):
            bias = float(rng.choice([0.0, 0.0, 0.0, rng.normal(0, 0.3)]))
            lo = c0 + int(rng.randint(0, 60))
            days.append((f"{cam} 2026-10-{d + 1:02d}", bias, range(lo, lo + int(rng.randint(5, 120))),
                         int(rng.randint(1, 6))))
    return days, float(rng.choice([0.02, 0.12]))


def tidal_days(seed):
    """A window as the cameras see it (after the review's fuzz2): two cameras on overlapping stretches
    of beach (c1 cells 0-119, c2 100-299), 3-14 days each (or 1-3 for a camera seen little), each day
    reaching a random stretch, as the tide allows, with 1-5 frames that each catch 80 % of its cells; a
    fifth of the days off by 1.2-4 x the limit, the rest within a quarter of it; noise 0.03-0.12 m, a
    limit of 0.1-0.2 m, values rounded to 1 cm in a quarter of the windows (ties).
    -> (pair_cell, pair_vals, pair_key, limit, min_frames)"""
    rng = np.random.RandomState(seed)
    ndays = int(rng.randint(3, 15))
    limit = float(rng.choice([0.1, 0.15, 0.2]))
    sigma = float(rng.choice([0.03, 0.08, 0.12]))
    days = []
    for cam, lo, hi in (("c1", 0, 120), ("c2", 100, 300)):
        nd = int(rng.randint(1, 4)) if rng.rand() < 0.15 else ndays
        start = int(rng.randint(0, 3))
        for d in range(nd):
            off = rng.rand() < 0.2
            b = rng.choice([-1, 1]) * rng.uniform(1.2 * limit, 4 * limit) if off else rng.normal(0, limit / 4)
            a = int(rng.randint(lo, hi - 30))
            w = int(rng.randint(20, hi - a + 1))
            days.append((f"{cam} 2026-10-{start + d + 1:02d}", b, np.arange(a, min(hi, a + w)),
                         int(rng.randint(1, 6))))
    pc, pv, pk = [], [], []
    for key, b, cells, nf in days:
        for _ in range(nf):
            sub = cells[rng.rand(len(cells)) < 0.8]
            pc.append(sub)
            pv.append(0.01 * sub + b + rng.normal(0, sigma, len(sub)))
            pk.append(np.array([key] * len(sub)))
    pc, pv, pk = np.concatenate(pc), np.concatenate(pv), np.concatenate(pk)
    if seed % 4 == 0:
        pv = np.round(pv, 2)
    return pc, pv, pk, limit, int(rng.choice([1, 3]))


def test_random(n=200):
    print(f"{n} random windows of two kinds: the bookkeeping holds, the WARNING and the record describe the "
          f"final DEM, the result does not depend on the sample order")
    keys = ("left out", "taken back", "refused", "not judged", "not judged, taken back", "judged in the end",
            "untested", "warnings")
    for kind in ("stretches", "tidal"):
        bad, stats = [], {k: 0 for k in keys}
        for seed in range(n):
            if kind == "stretches":
                days, sigma = random_days(seed)
                pc, pv, pk = pairs(days, sigma=sigma, seed=seed)
                limit, mf = LIMIT, MIN_FRAMES
            else:
                pc, pv, pk, limit, mf = tidal_days(seed)
            dt = dfc.day_consistency(pc, pv, pk, limit, mf)
            p = np.random.RandomState(seed + 1000).permutation(len(pc))
            why = check_invariants(dt, limit)
            if not why and dt != dfc.day_consistency(pc[p], pv[p], pk[p], limit, mf):
                why = "order dependent"
            if why:
                bad.append(f"seed {seed}: {why}")
            stats["left out"] += len(dt["steps"])
            stats["taken back"] += len(dt["taken_back"])
            stats["refused"] += len(dt["refused"])
            stats["not judged"] += len(dt["not_judged"])
            stats["not judged, taken back"] += sum(r["day"] in dt["taken_back"] for r in dt["not_judged"])
            stats["judged in the end"] += sum(r["judged_in_end"] for r in dt["not_judged"])
            stats["untested"] += len(dt["untested"])
            stats["warnings"] += bool(dfc.day_test_warning(dt, limit)[0])
        check(not bad, f"{kind}: invariants hold (" + ", ".join(f"{k} {v}" for k, v in stats.items()) + ")"
              + ("" if not bad else ": " + "; ".join(bad[:5])))


def archive(ndays, npre, seed=1):
    """A synthetic whole archive: two cameras, ndays days each; c1's first npre days 0.3-1.4 m high (a
    beach before a storm), one day in ten off by 0.2-0.5 m either way, the rest within ~0.05 m.
    Frames of 22 (c1) or 32 (c2) cells, 25 or 40 a day. -> (pairs, planted biases)"""
    rng = np.random.RandomState(seed)
    nx = 50
    cams = {"c1": (0, 22, 25), "c2": (18, 50, 40)}
    pc, pv, pk, truth = [], [], [], {}
    for d in range(ndays):
        day = (np.datetime64("2026-07-01") + d).astype(str)
        for cam, (a0, a1, nf) in cams.items():
            if cam == "c1" and d < npre:
                bias = 0.3 + 1.1 * rng.rand()
            elif rng.rand() < 0.1:
                bias = rng.choice([-1, 1]) * (0.2 + 0.3 * rng.rand())
            else:
                bias = rng.normal(0, 0.05)
            truth[f"{cam} {day}"] = bias
            lo = rng.randint(0, 20)
            hi = lo + rng.randint(20, 30)
            cols = np.arange(a0, a1)
            for f in range(nf):
                x = rng.randint(lo, hi)
                cells = cols * nx + x
                fn = rng.normal(0, 0.04)
                pc.append(cells)
                pv.append(0.02 * x + bias + fn + rng.normal(0, 0.03, len(cells)))
                pk.append(np.array([f"{cam} {day}"] * len(cells)))
    return (np.concatenate(pc), np.concatenate(pv), np.concatenate(pk)), truth


def test_timing(year=False):
    import time
    ndays, npre, limit_s = (365, 90, 600) if year else (90, 30, 30)
    print(f"run time: {2 * ndays} camera-days, many left out")
    (pc, pv, pk), truth = archive(ndays, npre)
    t0 = time.time()
    dt = dfc.day_consistency(pc, pv, pk, LIMIT, MIN_FRAMES)
    took = time.time() - t0
    bad = {k for k, b in truth.items() if abs(b) > LIMIT}
    check(took < limit_s and len(dt["rejected"]) >= 0.95 * len(bad) and set(dt["rejected"]) <= bad | set(
              k for k in truth if abs(truth[k]) > 0.1),
          f"{2 * ndays} camera-days, {len(pc):,} samples, {len(dt['steps'])} left out one at a time: "
          f"{took:.1f} s (limit {limit_s} s; the cron's flock skips a run that overlaps the next); "
          f"{len(set(dt['rejected']) & bad)} of the {len(bad)} planted off days found")


HEADER = ("source_file,camera,capture_time_utc,capture_epoch,easting_utm19,northing_utm19,"
          "tide_elevation_navd88,beach_elevation_navd88")


def ground_file(path, frames):
    """frames: [(day of Sep 2026, UTC hour, minute, beach row x, dz)]: a c2 frame of 20 points along
    row x (northing 2x + 1 m) at elevation 0.2 x - 0.5 + dz, the station's ground-file columns."""
    rows = [HEADER]
    for day, h, m, x, dz in frames:
        t = datetime(2026, 9, day, h, m, tzinfo=timezone.utc)
        f = f"f{day}_{h}_{m}_{x}"
        z = 0.2 * x - 0.5 + dz
        for col in range(20):
            rows.append(f"{f},c2,{t:%Y-%m-%dT%H:%M:%S}+00:00,{int(t.timestamp())},{420001.0 + 2 * col:.2f},"
                        f"{4638001.0 + 2 * x:.2f},{z:.3f},{z:.3f}")
    Path(path).write_text("\n".join(rows) + "\n")


def dem_rows(stem):
    """The DEM's rows south to north (beach row x = 0 first), NaN where blank."""
    a = np.loadtxt(f"{stem}_dem.asc", skiprows=6)[::-1]
    return np.where(a == -9999, np.nan, a)


def test_cli_cannot_judge(d):
    print("dem_from_contours.py: a fringe day, and days that all disagree")
    import survey_products as sp
    rng = np.random.RandomState(1)
    # the review's case: 9 good days on rows 2-9; the low beach (rows 0-1) reached by 20 Sep (good,
    # 3 frames a row, also on rows 2-4) and 22 Sep (0.57 m LOW, 16 frames a row)
    frames = [(day, 12 + i, 0, x, rng.normal(0, 0.02)) for day in range(11, 20) for i, x in enumerate(range(2, 10))]
    frames += [(20, 12 + r, 10 * x, x, rng.normal(0, 0.02)) for r in range(3) for x in (0, 1)]
    frames += [(20, 16 + r, 10 * x, x, rng.normal(0, 0.02)) for r in range(3) for x in range(2, 5)]
    frames += [(22, 4 + r, 3 * x, x, -0.57 + rng.normal(0, 0.02)) for r in range(16) for x in (0, 1)]
    g, stem = Path(d) / "fringe.csv", Path(d) / "fringe"
    ground_file(g, frames)
    r = subprocess.run([sys.executable, str(HERE / "dem_from_contours.py"), str(g), str(stem),
                        "--max-day-offset", str(LIMIT), "--no-plot"], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
    (Path(d) / "fringe.log").write_text(r.stdout)
    rej = [ln for ln in r.stdout.splitlines() if "REJECTED" in ln]
    z = dem_rows(stem)
    check(r.returncode == 0 and len(rej) == 1 and "c2 2026-09-22" in rej[0] and "WARNING" not in r.stdout,
          "a day 0.57 m low on cells only one other (checked) day reaches: it alone is REJECTED, no WARNING")
    check(abs(np.nanmean(z[0]) + 0.50) < 0.03 and abs(np.nanmean(z[1]) + 0.30) < 0.03,
          f"the low beach comes out right: rows 0 and 1 {np.nanmean(z[0]):+.3f}, {np.nanmean(z[1]):+.3f} m "
          f"(truth -0.50, -0.30; with that day kept, as d28b5ec did, -1.07 and -0.86)")
    # every day disagrees: one camera, 4 days at 0, +0.4, +0.8, +1.2 m on the same rows
    frames = [(day, 12 + r, 5 * x, x, 0.4 * i + rng.normal(0, 0.02))
              for i, day in enumerate((11, 12, 13, 14)) for r in range(3) for x in range(10)]
    g, stem = Path(d) / "apart.csv", Path(d) / "apart"
    ground_file(g, frames)
    r = subprocess.run([sys.executable, str(HERE / "dem_from_contours.py"), str(g), str(stem),
                        "--max-day-offset", str(LIMIT), "--no-plot"], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
    (Path(d) / "apart.log").write_text(r.stdout)
    warn = [ln for ln in r.stdout.splitlines() if "WARNING" in ln]
    check(r.returncode == 4 and "No points left: the day-consistency test rejected every camera-day" in r.stdout
          and not Path(f"{stem}_dem.asc").exists(),
          "4 days that all disagree: every camera-day rejected, exit 4, no grids (as the one-pass test)")
    check(len(warn) == 1 and "left out without being judged: c2 2026-09-11, c2 2026-09-12" in warn[0],
          "one WARNING line: the last two days could not be judged")
    check("rejected every camera-day (4 of 4)" in sp.dem_failure_reason(Path(d) / "apart.log"),
          "survey_products.dem_failure_reason reads it as before ('rejected every camera-day')")
    # three days at 0, +0.4, -0.4 m: 12 or 13 Sep goes (whichever the noise puts further out), then the
    # other two cannot be judged; the off one (over in the first pass) is left out, 11 Sep kept -- the
    # only day left, so untested
    frames = [(day, 12 + r, 5 * x, x, b + rng.normal(0, 0.02))
              for day, b in ((11, 0.0), (12, 0.4), (13, -0.4)) for r in range(3) for x in range(10)]
    g, stem = Path(d) / "three.csv", Path(d) / "three"
    ground_file(g, frames)
    r = subprocess.run([sys.executable, str(HERE / "dem_from_contours.py"), str(g), str(stem),
                        "--max-day-offset", str(LIMIT), "--no-plot"], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
    (Path(d) / "three.log").write_text(r.stdout)
    warn = [ln for ln in r.stdout.splitlines() if "WARNING" in ln]
    t = json.loads(Path(f"{stem}_info.json").read_text())["filters"]["day_test"]
    nj = {x["day"]: x for x in t["not_judged"]}
    row = [ln for ln in r.stdout.splitlines() if ln.startswith("   c2 2026-09-11")][0]
    second = t["left_out_order"][1] if len(t["left_out_order"]) == 2 else None
    check(r.returncode == 0 and len(warn) == 1 and second in ("c2 2026-09-12", "c2 2026-09-13")
          and warn[0].endswith(f"left out without being judged: {second}; kept untested: c2 2026-09-11.")
          and "see the warning below" in row and t["warning"] == warn[0]
          and nj["c2 2026-09-11"]["verdict"] == "kept untested" and nj[second]["verdict"] == "left out"
          and t["kept_untested"] == ["c2 2026-09-11"] and t["kept_over_limit"] == [],
          f"3 days 0, +0.4, -0.4 m: one line holds the word WARNING (11 Sep's own line says 'see the warning "
          f"below'); it and <stem>_info.json name {second} left out and 11 Sep kept untested")


def test_py38():
    print("Python 3.8 syntax")
    import ast
    for f in ("dem_from_contours.py", "test_dem_day_test.py"):
        src = (HERE / f).read_text(encoding="utf-8")
        try:
            ast.parse(src, feature_version=(3, 8))
        except TypeError:                        # Python 3.7 has no feature_version
            ast.parse(src)
        check(True, f"{f} parses as Python 3.8")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="keep the outputs, in a new folder inside this one")
    ap.add_argument("--year", action="store_true",
                    help="time a year of camera-days (730, 668k samples, ~170 left out; about a minute on the "
                         "station's Python 3.8 / numpy 1.17; at the station's density, 1.5M samples, 2.5 min of "
                         "CPU) instead of a quarter")
    args = ap.parse_args()
    if args.keep:
        Path(args.keep).mkdir(parents=True, exist_ok=True)
    d = os.path.realpath(tempfile.mkdtemp(prefix="dem_day_test_", dir=args.keep))
    start = os.getcwd()
    os.chdir(d)
    print(f"test folder: {d} (the current folder during the run)")
    try:
        test_engine()
        test_spearman()
        test_one_bad()
        test_masking()
        test_swamp()
        test_take_back()
        test_trend()
        test_trend_rule()
        test_trend_end()
        test_nuc_7day()
        test_too_few()
        test_fringe()
        test_cannot_judge()
        test_unjudged_fate()
        test_set_aside_retried()
        test_anchor()
        test_deterministic()
        test_random()
        test_timing(args.year)
        test_cli(d)
        test_cli_cannot_judge(d)
        test_py38()
    finally:
        os.chdir(start)
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    print("all DEM day-test tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

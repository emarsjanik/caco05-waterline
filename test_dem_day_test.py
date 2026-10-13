#!/usr/bin/env python3
"""
Tests For The DEM Day Test (dem_from_contours.py --max-day-offset)
=====================================================================
Plain python3 (no pytest needed on the station):

    python3 test_dem_day_test.py [--keep DIR]

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
  * a steady trend: the WARNING, and what the test trims;
  * the station's 7-day c1/c2 pattern of 5-11 Oct 2026, replayed;
  * fewer than 3 days, or a camera whose days only see each other: the
    test cannot tell which day is off and leaves none of them out;
  * the offsets equal the one-pass code's, whatever the order of the
    samples, and ties are broken by name;
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
from datetime import datetime

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
    check(loo.support("c2 2025-03-08") == 4.0 and loo.support("c2 2025-03-08", ref) == 3.0,
          "support: 8 Mar's DEM holds 4 other days on the low beach (7 Mar and three good days), 3 without 7 Mar")


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
    check(text.count("not taken back:") == 1 and "not taken back: c2 2025-03-03" in text
          and "would put c2 2025-03-02 over the limit" in text, "... and the lines say so, once")
    check(all(abs(dt["final"][k][0]) <= LIMIT for k in dt["final"] if k not in dt["rejected"]),
          "every kept day within the limit")


def test_trend():
    print("a steady trend: WARNING, and what the test trims")
    dt, one = run(scenario_trend(), sigma=NUC_SIGMA, seed=5)
    tr = dfc.day_trends({k: v[0] for k, v in dt["final"].items()}, LIMIT)
    check(len(tr) == 1 and tr[0]["camera"] == "c1" and tr[0]["spearman_rho"] == 1.0 and tr[0]["days"] == 7,
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
          and sum("WARNING" in x for x in w) == 1 and all(x.startswith("   ") for x in w[1:]),
          "one WARNING line (the station status counts them), naming the days left out")
    flat = dfc.day_trends({"c1 2026-10-0%d" % i: o for i, o in enumerate([0.0, 0.03, 0.05, 0.08, 0.1], 1)}, LIMIT)
    check(flat == [], "a monotonic run within the limit (range 0.10 m) is no warning")
    four = dfc.day_trends({"c1 2026-10-0%d" % i: o for i, o in enumerate([-0.2, -0.1, 0.1, 0.2], 1)}, LIMIT)
    check(four == [], "fewer than 5 testable days: no warning")


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
    tr = dfc.day_trends({k: v[0] for k, v in dt["final"].items()}, LIMIT)
    check([t["camera"] for t in tr] == ["c1"] and tr[0]["spearman_rho"] == 1.0,
          f"so the WARNING carries it: c1 rho {tr[0]['spearman_rho']:+.2f}, range {tr[0]['range_m']:.2f} m "
          f"(5 Oct left out, 11 Oct kept: the trend trimmed at its low end only); none for c2")
    check(all(abs(dt["final"][k][0]) <= LIMIT for k in dt["final"] if k not in dt["rejected"]),
          "every kept day within the limit of the kept days")


def test_too_few():
    print("the test cannot tell which day is off: none left out")
    days = [("c2 2025-03-06", 0.0, range(100), 3), ("c2 2025-03-07", 0.3, range(100), 3)]
    dt, one = run(days)
    check(dt["too_few"] and dt["rejected"] == [] and dt["steps"] == [], "two days 0.3 m apart: both kept")
    lines = dfc.day_test_report(dt, LIMIT)
    check(any("kept: fewer than 3 testable days" in x for x in lines) and "every day kept" in lines[-1]
          and not any("REJECTED" in x for x in lines), "... and the lines say so (as test_survey_products expects)")
    days = [("c2 2025-03-06", 0.0, range(100), 3), ("c2 2025-03-07", 0.4, range(100), 3),
            ("c2 2025-03-08", -0.4, range(100), 3)]
    dt, one = run(days)
    check(len(dt["rejected"]) == 1 and "only 2 testable camera-day(s) kept" in dt["stop"],
          f"three days, two off either way: one left out ({dt['rejected'][0]}), then the test stops with two")
    # c1's two days only see each other (c2 does not share their cells): neither can be judged
    days = [("c1 2026-10-0%d" % d, 0.3 * (d - 1), range(100), 3) for d in (1, 2)]
    days += [("c2 2026-10-0%d" % d, 0.0, range(100, 200), 3) for d in (1, 2, 3, 4, 5)]
    dt, one = run(days)
    check(one_pass_rejects(one) == ["c1 2026-10-01", "c1 2026-10-02"],
          "c1's two days only see each other: one pass rejected BOTH (7 testable days in all)")
    check(dt["rejected"] == [] and "fewer than 2 other kept camera-days" in dt["stop"],
          "iterative: neither left out (each one's DEM is the other day alone)")
    lines = dfc.day_test_report(dt, LIMIT)
    text = " ".join(x.strip() for x in lines)
    check("none left out: c1 2026-10-01, c1 2026-10-02 over the limit but kept" in text,
          "... and the lines say why")
    # untestable days (too few comparable samples) are kept and do not count
    days = [("c2 2025-03-0%d" % d, 0.0, range(100), 3) for d in (1, 2, 3, 4)]
    days.append(("c2 2025-03-05", 0.5, range(5), 3))                     # 15 samples < 30
    dt, one = run(days)
    check(np.isnan(dt["first"]["c2 2025-03-05"][0]) and dt["rejected"] == [],
          "a day with fewer than 30 comparable samples: kept, not testable")


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
    check(sum(ln.startswith("WARNING") for ln in out.splitlines()) == 1
          and "c1's camera-day offsets run one way with time over 7 days" in out, "one WARNING line, for c1")
    check("rejected every camera-day" not in sp.dem_failure_reason(Path(d) / "dem.log"),
          "survey_products.dem_failure_reason reads the block (no 'every camera-day rejected')")
    # the reader of the info file, as before
    import dem_figure
    _, _, left = dem_figure.title_lines(info, {"filled": 1, "blanked": 0, "no_value": 0}, "TEST", 0.5, 3, 2.0)
    check(any(f"{len(rej)} camera-days off the rest by > {LIMIT:g} m" in x for x in left),
          "the page (dem_figure.title_lines) still counts them")
    # no test: the info file says so
    r = subprocess.run([sys.executable, str(HERE / "dem_from_contours.py"), str(g), str(Path(d) / "plain"),
                        "--no-plot"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True,
                       cwd=str(d))
    f = json.loads((Path(d) / "plain_info.json").read_text())["filters"]
    check(r.returncode == 0 and f["day_test"] is None and f["camera_days_rejected"] == []
          and f["camera_days_rejected_offsets_m"] == {} and "Day consistency" not in r.stdout,
          "without --max-day-offset: no test, day_test None")


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
        test_nuc_7day()
        test_too_few()
        test_deterministic()
        test_cli(d)
        test_py38()
    finally:
        os.chdir(start)
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    print("all DEM day-test tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

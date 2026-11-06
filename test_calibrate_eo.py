#!/usr/bin/env python3
"""
Tests For The Live-Calibration Re-Solve (calibrate_eo.py)
============================================================
Plain python3 (no pytest needed on the station):

    python3 test_calibrate_eo.py [--keep DIR]

WHAT. Every synthetic test uses the REAL target geometry and lens (the
2025-11-13 targets, the 2024-08-01 IO, the live EO as the start), so it
tests this station's case, not an easy one:
  1 [SYNTH] a known EO (0.3 m, 0.2 deg off the live one) recovered from
    picks with 3 px of noise, within 3 sd of its own covariance; over 200
    seeded trials the Jacobian's 68% / 95% intervals (Student's t) cover
    the truth 60-76% / 90-99% of the time for tilt; the whitening (sqrt of
    the weight, over the pixel sigma), the jackknife's (n - 1)/n scale and
    the sds' units on known cases;
  2 [SYNTH] a 40 px blunder on one target: robust weight < 0.3, the robust
    solve nearer the truth on the ground, the target flagged in LOO;
  3 [SYNTH] multi-start: all starts reach one minimum; a start 5 deg off
    with --starts 1 converges to the same minimum, is said to be one start
    only, and does not "reproduce the start calibration";
  4 [SYNTH] --free auto: with the position 1 m off, "all" wins; with the
    position exact, "angles" wins or ties; a tie within 1 mm goes to the
    fewer parameters; the in-sample errors are the solved EO's;
  5 [SYNTH] a camera-position prior at the truth (0.03 m) shrinks the tilt
    sd; 6 starts and every leave-one-out fold use the same prior and sigma;
    SIGMA_H holds x and y, SIGMA_Z holds z;
  6 [SYNTH] the lens test flags a d1 changed by -0.1 (planted at georectify's
    own index of d1), and with the lens unchanged leaves >= 90% of 50
    seeded trials unflagged; leave-one-out with the lens freed uses each
    fold's own lens; a flag one target carries is said to rest on it;
  7 [SYNTH] the sea horizon on a rendered frame (tilt +0.3, roll -0.2 deg)
    recovered within 0.03 / 0.06 deg; frames 1 day before and 4 or 10 days
    after the EO's day refused, 3 days after used, 4 days after used with
    --horizon-max-days 5, 10 days after used with --horizon-any-epoch (with
    the moved-camera warning); another camera's or station's frame refused
    (also a link with a c1 name to a c2 frame); an inverted frame, a bent
    horizon (over horizon_check.FIT_OK_PX) and a grey frame not clear; a
    half-size frame refused; frames of 23 Oct 2024 (by name, by epoch, and
    through links with innocent names) refused before anything is solved
    and never read (cv2.imread is recorded with the audit hook's opens);
  8 the guards: inputs named for 23 Oct 2024 (also through a link with an
    innocent name, a chain whose middle link alone has the name, a folder
    named for the day, and 'link/..' paths that climb from where a link
    leads), the calibration history, survey files (also through a link):
    refused, exit 2, and never opened (an audit hook records every open);
    a dozen spellings of the day refused and other days not; --out-dir
    calibration/ refused; an existing output refused and left byte for
    byte; calibration names refused; usage errors (also a NaN or infinite
    --camera-position and a --seed out of range) exit 6;
  9 [REAL] the live inputs: the angles solve reproduces the live EO within
    0.10 m / 0.10 deg; leave-one-out RMS c1 angles 0.453, c2 angles 0.410,
    c1 all 0.589, c2 all 1.060 m (within 0.01); --free auto picks the
    angles for both; the cross-shore sigma, model choice, shift from the
    live EO and pick origin at the data's ranges, the angles' sds and the
    lens verdicts (each drop-one tested) are docs/VALIDATION.md 4.6c's; the
    refit EO reads back; nothing in calibration/ changed; with 6 of c1's
    targets the report does not speak of 14;
 10 the pinned modules are the version survey_products.py pins
    (gcpfree_eo.pin_problem); Python 3.8 syntax. test_gcpfree_eo.py
    (~7 min) and test_georectify_uncertainty.py are run separately.

STATION FILES. Reads only the 2025-11-13 CACO05 inputs in calibration/
(and the sigma tables' sha256, to show they are not touched). The held
files' names are tested on synthetic canaries in a temporary folder;
nothing is written outside it.
"""

import argparse
import contextlib
import io as _io
import json
import os
import shutil
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import calibrate_eo as CE                                   # noqa: E402
import georectify_uncertainty as GU                         # noqa: E402
from georectify import IO_KEYS, load_extrinsics, load_intrinsics    # noqa: E402

FAILED = []
WATCH = {"on": False, "opened": []}
CAL = HERE / "calibration"


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILED.append(msg)


def _audit(event, args):
    if WATCH["on"] and event == "open" and args and isinstance(args[0], (str, bytes, os.PathLike)):
        WATCH["opened"].append(os.fsdecode(args[0]))


def cli(argv):
    """calibrate_eo.main in this process (the audit hook sees its opens). -> (exit code, stdout + stderr)."""
    buf = _io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        try:
            rc = CE.main([str(a) for a in argv])
        except SystemExit as exc:
            rc = exc.code
    return rc, buf.getvalue()


def live(cam):
    """The live 2025-11-13 inputs of a camera: paths, io, eo, nums, xyz, uv."""
    paths = {"io": CAL / CE.IO_FILE.format(cam=cam), "eo": CAL / CE.EO_FILE.format(cam=cam),
             "targets": CAL / GU.TARGETS.format(cam=cam), "picks": CAL / GU.PICKS.format(cam=cam)}
    io, eo = load_intrinsics(str(paths["io"])), load_extrinsics(str(paths["eo"]))
    nums, xyz, uv = GU.load_targets(paths["targets"], paths["picks"])
    return paths, io, eo, nums, xyz, uv


def sd_line(sd):
    return ", ".join(f"{k} {v:.4f}" for k, v in sd.items())


def picks_from(io, eo, xyz, seed, noise=3.0):
    U, V = GU.project(xyz, io, eo)
    rs = np.random.RandomState(seed)
    return np.c_[U, V] + rs.normal(0.0, noise, (len(xyz), 2))


def perturbed(eo, dpos=(0.2, -0.2, 0.1), dang=(0.2, -0.2, 0.2)):
    e = eo.copy()
    e[:3] += np.asarray(dpos, float)
    e[3:] += np.deg2rad(np.asarray(dang, float))
    return e


def write_eo(path, eo):
    x, y, z, az, tilt, roll = eo[:3].tolist() + np.degrees(eo[3:]).tolist()
    Path(path).write_text(f"x: {x!r}\ny: {y!r}\nz: {z!r}\nazimuth: {az!r}\ntilt: {tilt!r}\nroll: {roll!r}\n"
                          f"#x - camera Easting (a synthetic start)\n")
    return Path(path)


def write_targets(d, stem, nums, xyz, uv):
    xp, up = Path(d) / f"{stem}_xyz.csv", Path(d) / f"{stem}_uv.csv"
    xp.write_text("".join(f"{n},{a!r},{b!r},{c!r}\n" for n, (a, b, c) in zip(nums, xyz.tolist())))
    up.write_text("num,U,V\n" + "".join(f"{n},{a!r},{b!r}\n" for n, (a, b) in zip(nums, uv.tolist())))
    return xp, up


# ---------------------------------------------------------------------
# 1-6: the solver, on the real geometry
# ---------------------------------------------------------------------

def test_recover():
    print("1 [SYNTH] a known EO recovered, its Jacobian intervals' coverage (c1 targets and lens, 6 parameters)")
    _, io, eo, _, xyz, _ = live("c1")
    truth = perturbed(eo)
    free = CE.FREE_SETS["all"]
    f, s = CE.plain_solve(io, eo, xyz, picks_from(io, truth, xyz, 0), free)
    C, _ = CE.cov_jacobian(f)
    z = np.abs(f["p"] - truth) / np.sqrt(np.diag(C))
    check(f["converged"] and np.all(z <= 3.0) and 1.5 < s < 4.5,
          f"recovered within 3 sd of its own covariance (|error|/sd up to {z.max():.2f}; pixel sigma {s:.2f} px "
          f"for 3 px of noise)")
    t68, t95 = CE.t_quantile(0.6827, 2 * len(xyz) - 6), CE.t_quantile(0.95, 2 * len(xyz) - 6)
    zs = []
    for seed in range(200):
        f, _ = CE.plain_solve(io, eo, xyz, picks_from(io, truth, xyz, seed), free)
        C, _ = CE.cov_jacobian(f)
        zs.append(np.abs(f["p"] - truth) / np.sqrt(np.diag(C)))
    zs = np.array(zs)
    cov = {k: (float(np.mean(zs[:, i] <= t68)), float(np.mean(zs[:, i] <= t95)))
           for k, i in (("azimuth", 3), ("tilt", 4), ("roll", 5))}
    c68, c95 = cov["tilt"]
    check(0.60 <= c68 <= 0.76 and 0.90 <= c95 <= 0.99,
          f"200 trials: the 68% / 95% intervals (t {t68:.3f} / {t95:.3f}) cover the tilt {100 * c68:.1f}% / "
          f"{100 * c95:.1f}% of the time (azimuth {100 * cov['azimuth'][0]:.1f} / {100 * cov['azimuth'][1]:.1f}%, "
          f"roll {100 * cov['roll'][0]:.1f} / {100 * cov['roll'][1]:.1f}%)")
    # the solver's internals the report's numbers rest on
    n, uv0 = len(xyz), picks_from(io, truth, xyz, 0)
    fw = CE.solve(io, eo, xyz, uv0, free, sigma_px=1.0, weights=np.full(n, 4.0), iterations=0)
    fs = CE.solve(io, eo, xyz, uv0, free, sigma_px=2.0, iterations=0)
    check(np.allclose(fw["r"], 2.0 * fw["resid_px"]) and np.allclose(fs["r"], 0.5 * fs["resid_px"]),
          "the whitening: a weight of 4 doubles a residual (sqrt(w), so the robust IRLS is soft-L1) and a pixel "
          "sigma of 2 px halves it")
    x = np.random.RandomState(3).normal(0.0, 1.0, 15)
    means = np.array([[np.delete(x, k).mean(), 2.0 * np.delete(x, k).mean()] for k in range(15)])
    jk = CE.cov_jackknife(means)
    v = np.var(x, ddof=1) / 15
    check(np.allclose(jk, [[v, 2 * v], [2 * v, 4 * v]], rtol=1e-12, atol=0),
          f"the jackknife covariance of leave-one-out means is var(x)/n ({jk[0, 0]:.5f} against {v:.5f}: the "
          f"(n - 1)/n scale)")
    C = np.diag(np.r_[np.array([0.01, 0.02, 0.03]) ** 2, np.deg2rad([0.1, 0.2, 0.3]) ** 2])
    sd = CE.sds(C, free)
    want = {"x": 0.01, "y": 0.02, "z": 0.03, "azimuth": 0.1, "tilt": 0.2, "roll": 0.3}
    check(list(sd) == list(want) and all(abs(sd[k] - want[k]) < 1e-12 for k in want),
          f"sds: metres for x y z, degrees for the angles ({sd_line(sd)})")


def test_blunder():
    print("2 [SYNTH] a 40 px blunder on one target (c1, #117 at 68 m)")
    _, io, eo, nums, xyz, _ = live("c1")
    truth = perturbed(eo, (0, 0, 0))
    uv = picks_from(io, truth, xyz, 11)
    k = nums.index(117)
    uv[k] += [40 / np.sqrt(2), 40 / np.sqrt(2)]
    c = CE.candidate("angles", io, eo, xyz, uv, None, 1, 0)
    rb = CE.robust_solve(io, c["eo"], xyz, uv, c["free"], c["sigma_px"])
    nd = CE.nodes(io, truth, "c1")
    d_plain = CE.rms(CE.cross_shore(io, truth, c["eo"], nd))
    d_rob = CE.rms(CE.cross_shore(io, truth, rb["eo"], nd))
    others = np.delete(rb["weights"], k)
    check(rb["weights"][k] < 0.3 and others.min() > rb["weights"][k],
          f"its robust weight {rb['weights'][k]:.2f} (< 0.3; the others {others.min():.2f}-{others.max():.2f})")
    check(d_rob < d_plain, f"the robust solve is nearer the truth on the ground: {d_rob:.3f} m against the plain "
                           f"{d_plain:.3f} m (cross-shore RMS at 60-90 m)")
    check(k in c["loo_flagged"], f"leave-one-out flags it (over {CE.LOO_FLAG:g} x the median): "
                                 f"{[nums[j] for j in c['loo_flagged']]}")


def test_multistart(d):
    print("3 [SYNTH] multi-start: one minimum; a start 5 deg off with --starts 1")
    paths, io, eo, nums, xyz, _ = live("c1")
    truth = perturbed(eo)
    uv = picks_from(io, truth, xyz, 5)
    c = CE.candidate("all", io, eo, xyz, uv, None, 12, 0)
    ms = c["multistart"]
    check(ms["n"] == 12 and ms["one_minimum"] and not ms["lower_elsewhere"]
          and max(m[1] for m in ms["start_moves"]) > 1.0,
          f"12 starts (up to {max(m[0] for m in ms['start_moves']):.2f} m, {max(m[1] for m in ms['start_moves']):.2f} "
          f"deg off) reach one minimum (spread {ms['spread_m']:.1e} m, {ms['spread_deg']:.1e} deg)")
    bad = eo.copy()
    bad[4] += np.deg2rad(5.0)
    bad[3] -= np.deg2rad(5.0)
    eo_bad = write_eo(Path(d) / "CACO05_c1_20251113_EO-badstart.yaml", bad)
    xp, up = write_targets(d, "synth_c1", nums, xyz, uv)
    out = Path(d) / "out_badstart"
    rc, log = cli(["--out-dir", out, "--cameras", "c1", "--io", paths["io"], "--eo", eo_bad, "--targets", xp,
                   "--picks", up, "--starts", 1, "--samples", 50, "--no-io-test", "--free", "all"])
    j = json.loads((out / "c1_eo_report.json").read_text()) if rc == 0 else {}
    ok = rc == 0 and j["candidates"]["all"]["multistart"]["n"] == 1
    same = ok and max(abs(j["eo_refit"][k] - CE.eo_deg(c["eo"])[k]) for k in CE.EO_NAMES) < 1e-3
    check(ok and same and not j["reproduces_start"] and "1 start only" in log and "does NOT reproduce" in log,
          "a start 5 deg off (azimuth and tilt), --starts 1: it converges to the same minimum (within 0.001), says "
          "'1 start only' and that it does NOT reproduce the start calibration (it moved "
          + (f"{j['candidates']['all']['move_from_start_m_deg'][1]:.2f} deg)" if ok else "?)"))
    free = CE.FREE_SETS["all"]
    stuck = CE.solve(io, bad, xyz, uv, free, sigma_px=c["sigma_px"], iterations=1)     # a solve stopped early
    ms = CE.multistart(io, eo, xyz, uv, free, c["sigma_px"], None, stuck, 4, 0)
    check(ms["lower_elsewhere"] and ms["best"] is not stuck and ms["best"]["cost"] < stuck["cost"]
          and max(abs(CE.eo_deg(ms["best"]["eo"])[k] - CE.eo_deg(c["eo"])[k]) for k in CE.EO_NAMES) < 1e-3,
          f"a solve from the start EO that stopped short of the minimum (cost {stuck['cost']:.0f} against "
          f"{ms['best']['cost']:.1f}): the lower minimum another start reached is used (lower_elsewhere: the report's WARNING line)")


def test_auto():
    print("4 [SYNTH] --free auto takes the set with the lower leave-one-out ground RMS")
    _, io, eo, _, xyz, _ = live("c1")
    uv = picks_from(io, perturbed(eo, (0.577, 0.577, 0.577)), xyz, 7)
    cands = {k: CE.candidate(k, io, eo, xyz, uv, None, 1, 0) for k in ("all", "angles")}
    name, why = CE.choose(cands, "auto")
    check(name == "all" and cands["all"]["loo_rms"] < cands["angles"]["loo_rms"],
          f"the position 1 m off: 'all' wins (LOO all {cands['all']['loo_rms']:.3f} m, angles "
          f"{cands['angles']['loo_rms']:.3f} m)")
    worst = 0.0
    for k, c in cands.items():
        dE, dN = GU.ground_errors(io, c["eo"], xyz, uv)
        worst = max(worst, float(np.max(np.abs(c["in_ground"] - np.hypot(dE, dN)))))
    check(worst < 1e-9 and CE.eo_move(cands["all"]["eo"], eo)[0] > 0.5,
          f"the in-sample ground errors are the solved EO's (GU.ground_errors; off by {worst:.1e} m), not the start "
          f"EO's ({CE.eo_move(cands['all']['eo'], eo)[0]:.2f} m away)")
    uv = picks_from(io, perturbed(eo, (0, 0, 0)), xyz, 0)
    cands = {k: CE.candidate(k, io, eo, xyz, uv, None, 1, 0) for k in ("all", "angles")}
    name, why = CE.choose(cands, "auto")
    check(cands["angles"]["loo_rms"] <= cands["all"]["loo_rms"] + 0.01 and
          name == ("angles" if cands["angles"]["loo_rms"] <= cands["all"]["loo_rms"] + CE.TIE_M else "all"),
          f"the position exact: 'angles' wins or ties (LOO angles {cands['angles']['loo_rms']:.3f} m, all "
          f"{cands['all']['loo_rms']:.3f} m; chosen {name})")
    check(CE.choose(cands, "all")[0] == "all" and CE.choose(cands, "angles")[0] == "angles",
          "--free all / angles take that set")

    def picked(**loo):
        return CE.choose({k: {"loo_rms": v, "free": list(CE.FREE_SETS[k])} for k, v in loo.items()}, "auto")[0]
    cases = (((0.4531, 0.4534), "angles"), ((0.4536, 0.4534), "angles"), ((0.4546, 0.4534), "all"),
             ((0.44, 0.41), "all"), ((0.41, 0.44), "angles"))
    got = [picked(angles=a, all=b) for (a, b), _ in cases]
    check(got == [w for _, w in cases],
          "the tie rule: within 1 mm of the lowest (0.4531 / 0.4534, 0.4536 / 0.4534 m) the fewer parameters win; "
          f"1.2 mm or 3 cm apart the lower wins ({got})")


def test_prior():
    print("5 [SYNTH] a camera-position prior at the truth (0.03 m)")
    _, io, eo, _, xyz, _ = live("c1")
    truth = perturbed(eo)
    uv = picks_from(io, truth, xyz, 3)
    pr = (float(truth[0]), float(truth[1]), float(truth[2]), 0.03, 0.03)
    a = CE.candidate("all", io, eo, xyz, uv, None, 1, 0)
    b = CE.candidate("all+prior", io, eo, xyz, uv, pr, 1, 0)
    ta, tb = CE.sds(a["cov"]["jacobian"], a["free"])["tilt"], CE.sds(b["cov"]["jacobian"], b["free"])["tilt"]
    check(tb < 0.75 * ta and np.linalg.norm(b["eo"][:3] - truth[:3]) < np.linalg.norm(a["eo"][:3] - truth[:3]),
          f"tilt sd {ta:.3f} -> {tb:.3f} deg (Jacobian); the position {np.linalg.norm(a['eo'][:3] - truth[:3]):.3f} -> "
          f"{np.linalg.norm(b['eo'][:3] - truth[:3]):.3f} m from the truth; reduced chi-square {b['chi2_red']:.2f}")
    # the prior is used by every start and every leave-one-out fold, each sigma where it belongs
    m = CE.candidate("all+prior", io, eo, xyz, uv, pr, 6, 0)
    ms, mv = m["multistart"], CE.eo_move(m["eo"], b["eo"])
    check(ms["n"] == 6 and ms["one_minimum"] and not ms["lower_elsewhere"] and mv[0] < 1e-3 and mv[1] < 1e-4,
          f"6 starts with the prior reach the 1-start solve's minimum ({mv[0]:.1e} m, {mv[1]:.1e} deg from it): the "
          f"extra starts are solved with the prior too")
    keep = np.arange(len(xyz)) != 0
    fold = CE.solve(io, b["eo"], xyz[keep], uv[keep], b["free"], sigma_px=b["sigma_px"], prior=pr)
    dp = float(np.max(np.abs(fold["p"] - b["loo"]["sols"][0])))
    check(dp < 1e-9, f"a leave-one-out fold of all+prior is CE.solve with the same pixel sigma and prior on the kept "
                     f"targets (max |dp| {dp:.1e})")
    pr2 = (float(truth[0]), float(truth[1]), float(truth[2]), 0.01, 1.0)
    c = CE.candidate("all+prior", io, eo, xyz, uv, pr2, 1, 0)
    sd = CE.sds(c["cov"]["jacobian"], c["free"])
    check(max(sd["x"], sd["y"]) <= 0.01 * (1 + 1e-6) and 3 * max(sd["x"], sd["y"]) < sd["z"] <= 1.0,
          f"--camera-position E N Z 0.01 1.0: SIGMA_H holds x and y, SIGMA_Z holds z (posterior sd x {sd['x']:.4f}, "
          f"y {sd['y']:.4f}, z {sd['z']:.4f} m)")
    return truth, uv, pr


def test_prior_cli(d, truth, uv, pr):
    print("5b [SYNTH] --camera-position on the command line")
    paths, io, eo, nums, xyz, _ = live("c1")
    xp, up = write_targets(d, "prior_c1", nums, xyz, uv)
    out = Path(d) / "out_prior"
    rc, log = cli(["--out-dir", out, "--cameras", "c1", "--targets", xp, "--picks", up, "--camera-position"]
                  + list(pr) + ["--starts", 2, "--samples", 50, "--no-io-test", "--free", "all"])
    j = json.loads((out / "c1_eo_report.json").read_text()) if rc == 0 else {}
    ok = rc == 0 and j["chosen"] == "all+prior" and set(j["candidates"]) == {"all", "angles", "all+prior"}
    yt = (out / "CACO05_c1_20251113_EO-refit.yaml").read_text() if rc == 0 else ""
    check(ok and "camera-position prior" in yt and "all+prior uses the surveyed one" in log,
          "three candidates (all, angles, all+prior); --free all with a position takes all+prior; the refit EO "
          "names the prior")


def test_lens():
    print("6 [SYNTH] the lens test (c1; d1 freed with the angles and with all six, judged by leave-one-out)")
    _, io, eo, nums, xyz, _ = live("c1")
    d1 = IO_KEYS.index("d1")                     # planted by georectify's own index, not calibrate_eo's
    check(CE.D1 == d1 == 6 and CE.D2 == IO_KEYS.index("d2") == 7,
          f"the lens terms freed are d1 and d2 (CIRN IO indices {CE.D1} and {CE.D2})")
    truth = perturbed(eo, (0, 0, 0))
    io_bad = io.copy()
    io_bad[d1] -= 0.1

    def verdict(lens_io, seed, drop_one=False):
        uv = picks_from(lens_io, truth, xyz, seed)
        cands = {k: CE.candidate(k, io, eo, xyz, uv, None, 1, 0) for k in ("all", "angles")}
        name, _ = CE.choose(cands, "auto")
        tests, verdicts = CE.lens_tests(io, eo, xyz, uv, cands, name, drop_one=drop_one, term_sets=((CE.D1,),))
        return verdicts[0], tests

    v, tests = verdict(io_bad, 0, drop_one=True)
    t = next(x for x in tests if x["set"] == v["with_set"])
    check(v["flagged"] and v["with"] < v["without"] and abs(t["io_values"][0] - io_bad[d1]) < 0.05
          and v["terms"] == ["d1"] and v["carried_by"] == [],
          f"d1 changed by -0.1 (seed 0): flagged, the best LOO {v['without']:.3f} m ({v['without_set']}) -> "
          f"{v['with']:.3f} m ({v['with_set']} + d1), d1 fitted {t['io_values'][0]:+.3f} (true "
          f"{io_bad[d1]:+.3f}); each target dropped in turn keeps the flag")
    # leave-one-out with the lens freed: the left-out target is projected with its fold's own lens and EO
    uv = picks_from(io_bad, truth, xyz, 0)
    ang = CE.FREE_SETS["angles"]
    full = CE.solve(io, eo, xyz, uv, ang, sigma_px=3.0, io_free=(CE.D1,))
    lo = CE.leave_one_out(full["io"], full["eo"], xyz, uv, ang, 3.0, io_free=(CE.D1,))
    keep = np.arange(len(xyz)) != 0
    fold = CE.solve(full["io"], full["eo"], xyz[keep], uv[keep], ang, sigma_px=3.0, io_free=(CE.D1,))
    dE, dN = GU.ground_errors(fold["io"], fold["eo"], xyz[:1], uv[:1])
    leak = max(abs(lo["dE"][0] - dE[0]), abs(lo["dN"][0] - dN[0]))
    check(leak < 1e-9 and abs(fold["io"][d1] - full["io"][d1]) > 1e-4,
          f"leave-one-out with d1 freed projects the left-out target with its fold's own d1 "
          f"({fold['io'][d1]:+.4f}, not the full fit's {full['io'][d1]:+.4f}; off by {leak:.1e} m)")
    # the drop-one verdict on a hand-made case: one target carries the flag
    fake_c = {"angles": {"loo_rms": 0.50}, "all": {"loo_rms": 0.60}}
    rows_a = [(j, 0.50, 0.40, 0.2, True) for j in range(5)]
    rows_a[2] = (2, 0.45, 0.43, 1 - 0.43 / 0.45, False)
    rows_b = [(j, 0.60, 0.55, 1 - 0.55 / 0.60, False) for j in range(5)]
    fake_t = [{"set": "angles", "terms": ["d1"], "loo_rms_with": 0.40, "drop_one": rows_a},
              {"set": "all", "terms": ["d1"], "loo_rms_with": 0.55, "drop_one": rows_b}]
    fv = CE.lens_verdict(fake_c, fake_t, ["d1"])
    txt = CE.drop_text(fv, [101, 102, 103, 104, 105])
    check(fv["flagged"] and fv["carried_by"] == [2] and "RESTS ON ONE TARGET: without 103" in txt
          and "no longer flagged" in txt,
          f"a flag that one target carries (hand-made drop-one rows): carried_by {fv['carried_by']}, '{txt[:58]}...'")
    hits = sum(verdict(io_bad, seed)[0]["flagged"] for seed in range(10))
    check(hits >= 8, f"d1 changed by -0.1: flagged in {hits} of 10 seeded trials (>= 8; a scratch run of 30 seeds: "
                     f"29; where it is not, all six absorb the lens into the position)")
    flagged, t0 = 0, time.time()
    for seed in range(50):
        flagged += verdict(io, 100 + seed)[0]["flagged"]
    check(flagged <= 5, f"the lens unchanged: flagged in {flagged} of 50 seeded trials (<= 5 allowed; "
                        f"{time.time() - t0:.0f} s)")


# ---------------------------------------------------------------------
# 7: the sea horizon
# ---------------------------------------------------------------------

def render_frame(path, io, eo, seed=1, uniform=None, offset=None, invert=False, half=False):
    """A synthetic timex frame: sky above the EO's horizon, sea below (offset(u): px added to the horizon row of
    each column; invert: dark above and bright below; half: at half the IO's size)."""
    import cv2
    from estimate_eo_rotation import horizon_rows
    nu, nv = int(io[0]), int(io[1])
    if uniform is not None:
        img = np.full((nv, nu), uniform, np.uint8)
    else:
        u = np.arange(nu, dtype=float)
        rows = horizon_rows(io, eo, u)
        if offset is not None:
            rows = rows + offset(u)
        h = np.where(np.isfinite(rows), rows, -1e9)[None, :]
        y = np.arange(nv, dtype=float)[:, None]
        frac = np.clip(h - y, 0.0, 1.0)               # share of each pixel above the horizon: sky
        img = 70.0 + 120.0 * frac + np.random.RandomState(seed).normal(0.0, 4.0, (nv, nu))
        img = np.clip(img, 0, 255).astype(np.uint8)
    if invert:
        img = 255 - img
    if half:
        img = cv2.resize(img, (nu // 2, nv // 2))
    cv2.imwrite(str(path), img, [cv2.IMWRITE_JPEG_QUALITY, 92])
    return Path(path)


def frame_name(d, day, cam, station="CACO05", hour=16):
    t = datetime(day[0], day[1], day[2], hour, 0, tzinfo=timezone.utc)
    return Path(d) / f"{int(t.timestamp())}.{t:%a.%b.%d_%H_%M_%S}.GMT.{t:%Y}.{station}.{cam}.timex.jpg"


def opened(path):
    """Whether the watched opens (the audit hook's, and cv2.imread's while test_horizon records them) include
    path, by its name or by the file it reaches."""
    real = os.path.realpath(str(path))
    return any(Path(p).name == Path(path).name or os.path.realpath(p) == real for p in WATCH["opened"])


def test_horizon(d):
    print("7 [SYNTH] the sea horizon (diagnostic only), c1's lens and live EO")
    import cv2
    real_imread = cv2.imread

    def imread(path, *a, **k):            # cv2 reads in C++, unseen by the audit hook: record its reads too
        if WATCH["on"]:
            WATCH["opened"].append(os.fsdecode(str(path)))
        return real_imread(path, *a, **k)
    cv2.imread = imread
    try:
        _test_horizon(d)
    finally:
        cv2.imread = real_imread


def _test_horizon(d):
    d = Path(d) / "frames"
    d.mkdir()
    paths, io, eo, *_ = live("c1")
    tr = eo.copy()
    tr[4] += np.deg2rad(0.3)
    tr[5] += np.deg2rad(-0.2)
    good = render_frame(frame_name(d, (2025, 11, 14), "c1"), io, tr)
    late = render_frame(frame_name(d, (2025, 11, 23), "c1"), io, tr, seed=2)
    grey = render_frame(frame_name(d, (2025, 11, 15), "c1"), io, tr, uniform=128)
    other = render_frame(frame_name(d, (2025, 11, 14), "c2"), io, tr)
    r = CE.horizon_diagnostic("c1", [good, late, grey, other], io, paths["eo"], {"start": (io, eo)})
    st = {Path(f["frame"]).name: f for f in r["frames"]}
    m = r["medians"].get("start", {})
    check(st[good.name]["status"] == "used" and m.get("frames") == 1 and abs(m["dtilt_deg"] - 0.3) <= 0.03
          and abs(m["droll_deg"] + 0.2) <= 0.06,
          f"rendered at tilt +0.3, roll -0.2 deg: read {m.get('dtilt_deg', float('nan')):+.3f} / "
          f"{m.get('droll_deg', float('nan')):+.3f} deg (clear: prominence "
          f"{st[good.name].get('edge', {}).get('prominence', 0):.0f}, sky above "
          f"{st[good.name].get('edge', {}).get('sky_above', 0):.2f}; {st[good.name].get('cols_in_frame')} grid columns "
          f"in the frame)")
    check(st[late.name]["status"] == "refused" and "10 day(s) after" in st[late.name]["doubt"] and not r["warnings"],
          f"a frame 10 days after the EO's day: refused ({st[late.name]['doubt'][:60]}...)")
    check(st[grey.name]["status"] == "not clear", f"a uniform grey frame: not clear ({st[grey.name]['doubt']})")
    check(st[other.name]["status"] == "refused" and "c2" in st[other.name]["doubt"],
          f"a frame of c2: refused ({st[other.name]['doubt']})")
    r2 = CE.horizon_diagnostic("c1", [late], io, paths["eo"], {"start": (io, eo)}, any_epoch=True)
    check(r2["frames"][0]["status"] == "used" and r2["warnings"] and "may have moved" in r2["warnings"][0],
          "with --horizon-any-epoch the late frame is used, with the warning that the camera may have moved")
    # the date window, the station, the clear-horizon tests and the frame size
    before, day3, day4 = (frame_name(d, day, "c1") for day in ((2025, 11, 12), (2025, 11, 16), (2025, 11, 17)))
    caco03 = frame_name(d, (2025, 11, 14), "c1", station="CACO03")
    for p in (before, day3, day4, caco03):
        shutil.copy(str(good), str(p))
    inverted = render_frame(frame_name(d, (2025, 11, 14), "c1", hour=17), io, tr, invert=True)
    bent = render_frame(frame_name(d, (2025, 11, 14), "c1", hour=18), io, tr,
                        offset=lambda u: 150.0 * ((u - 550.0) / 400.0) ** 2 - 50.0)
    half = render_frame(frame_name(d, (2025, 11, 14), "c1", hour=19), io, tr, half=True)
    (d / "links").mkdir()
    to_c2 = frame_name(d / "links", (2025, 11, 14), "c1", hour=20)       # a c1 name for the c2 frame
    os.symlink(str(other), str(to_c2))
    frames = [before, day3, day4, caco03, inverted, bent, half, to_c2]
    r3 = CE.horizon_diagnostic("c1", frames, io, paths["eo"], {"start": (io, eo)})
    st = {Path(f["frame"]).name: f for f in r3["frames"]}
    want = ((before, "refused", "1 day(s) before", "a frame of the day before the EO's"),
            (day3, "used", "", "a frame 3 days after it (the window's end)"),
            (day4, "refused", "4 day(s) after", "a frame 4 days after it"),
            (caco03, "refused", "station CACO03", "a frame of another station (CACO03, 14 Nov 2025)"),
            (inverted, "not clear", "dark above", "an inverted frame (dark above, bright below)"),
            (bent, "not clear", "FIT_OK_PX", "a horizon that tilt and roll cannot explain (bent by up to 100 px)"),
            (half, "refused", "1224 x 1024 px, not the IO's 2448 x 2048", "a frame at half the IO's size"),
            (to_c2, "refused", "a frame of c2", "a link with a c1 name to a c2 frame (screened on the file it "
                                                "reaches)"))
    for p, status, says, what in want:
        f = st[p.name]
        check(f["status"] == status and says in (f.get("doubt") or ""),
              f"{what}: {f['status']}" + (f" ({f['doubt'][:70]})" if f.get("doubt") else ""))
    check(not r3["warnings"] and r3["medians"]["start"]["frames"] == 1, "no frame outside the window is used, and no "
                                                                        "moved-camera warning")
    out5 = Path(d).parent / "out_horizon_5days"
    rc, log = cli(["--out-dir", out5, "--cameras", "c1", "--horizon", day4, "--horizon-max-days", 5, "--starts", 1,
                   "--samples", 50, "--no-io-test"])
    j = json.loads((out5 / "c1_eo_report.json").read_text()) if rc == 0 else {}
    check(rc == 0 and j["horizon"]["frames"][0]["status"] == "used" and not j["horizon"]["warnings"],
          "--horizon-max-days 5 on the command line: the frame 4 days after the EO's day is used")
    out = Path(d).parent / "out_horizon"
    rc, log = cli(["--out-dir", out, "--cameras", "c1", "--horizon", good, late, grey, "--horizon-any-epoch",
                   "--starts", 1, "--samples", 50, "--no-io-test"])
    txt = (out / "c1_eo_report.txt").read_text() if rc == 0 else ""
    check(rc == 0 and "THE SEA HORIZON" in txt and "WARNING" in txt and "may have moved" in txt
          and "never a constraint" in txt and "against the start EO" in txt and "against the solved EO" in txt,
          "the command line (--horizon-any-epoch): the report has the horizon against the start and solved EO, "
          "and the moved-camera warning")
    WATCH["opened"], WATCH["on"] = [], True
    try:
        rc, log = cli(["--out-dir", Path(d).parent / "out_watch", "--cameras", "c1", "--horizon", good, "--starts", 1,
                       "--samples", 20, "--no-io-test"])
    finally:
        WATCH["on"] = False
    check(rc == 0 and opened(good), "the watch sees a frame being read (cv2.imread is recorded)")
    held = Path(d) / "1729699200.Wed.Oct.23_16_00_00.GMT.2024.CACO03.c1.timex.jpg"
    bare = Path(d) / "1729699200.c1.timex.jpg"          # the same day, by its epoch only
    for p in (held, bare):
        shutil.copy(str(good), str(p))
    link_held = frame_name(d / "links", (2025, 11, 15), "c1")          # innocent frame names for them
    link_bare = frame_name(d / "links", (2025, 11, 16), "c1")
    os.symlink(str(held), str(link_held))
    os.symlink(str(bare), str(link_bare))
    for p, target, what in ((held, held, "by its name"), (bare, bare, "by its epoch only"),
                            (link_held, held, "through a link with a 14 Nov 2025 name"),
                            (link_bare, bare, "by its epoch only, through a link with a 14 Nov 2025 name")):
        out = Path(d).parent / "out_heldframe"
        WATCH["opened"], WATCH["on"] = [], True
        try:
            rc, log = cli(["--out-dir", out, "--cameras", "c1", "--horizon", p])
        finally:
            WATCH["on"] = False
        check(rc == 2 and "REFUSED" in log and "23 Oct 2024" in log and not opened(target) and "solving" not in log
              and not out.exists(), f"a frame of 23 Oct 2024 {what}: exit 2 before anything is solved, never "
                                    f"opened (cv2 included)")


# ---------------------------------------------------------------------
# 8: guards
# ---------------------------------------------------------------------

def test_guards(d):
    print("8 the guards: held, history and survey inputs never opened; outputs never overwritten")
    paths, io, eo, nums, xyz, uv = live("c1")
    d = Path(d) / "guards"
    (d / "2024-10-23").mkdir(parents=True)
    canary_t = d / "2024-10-23_Marconi_Extrinsic_Targets_c1_xyz.csv"     # synthetic canaries, not the held files
    shutil.copy(str(paths["targets"]), str(canary_t))
    canary_eo = d / "CACO03_c1_20241023_EO.yaml"
    shutil.copy(str(paths["eo"]), str(canary_eo))
    canary_io = d / "CACO03_c1_20241023_IO.yaml"
    shutil.copy(str(paths["io"]), str(canary_io))
    canary_uv = d / "gcp_uv_c1_2024_10_23.csv"
    shutil.copy(str(paths["picks"]), str(canary_uv))
    link = d / "innocent_targets.csv"
    os.symlink(str(canary_t), str(link))
    link2 = d / "innocent_too.csv"
    os.symlink(str(link), str(link2))
    in_dir = d / "2024-10-23" / "targets.csv"
    shutil.copy(str(paths["targets"]), str(in_dir))
    history = d / "calibration_history.csv"
    history.write_text("a synthetic canary\n")
    tif, asc, lidar = d / "survey.tif", d / "survey.asc", d / "CACO05_c1_2025-01-18_to_2025-01-23_lidar_EO.yaml"
    for p in (tif, asc):
        p.write_text("ncols 1\n")
    shutil.copy(str(paths["eo"]), str(lidar))
    # a chain whose middle link alone carries the day's name, ending at an innocently named file
    plain = d / "plain.csv"
    shutil.copy(str(paths["targets"]), str(plain))
    mid = d / "2024-10-23_mid.csv"
    os.symlink(str(plain), str(mid))
    chain = d / "innocent_chain.csv"
    os.symlink(str(mid), str(chain))
    # an innocent name for a raster
    raster = d / "area_2025.tif"
    raster.write_text("ncols 1\n")
    to_tif = d / "surveylink_targets.csv"
    os.symlink(str(raster), str(to_tif))
    # 'link/..': the kernel climbs from where the link leads, not from where it stands
    for sub in ("hist/sub", "2024-10-23/sub", "surv/sub", "x"):
        (d / sub).mkdir(parents=True, exist_ok=True)
    os.symlink(str(history), str(d / "hist" / "innocent.csv"))
    points = d / "surv" / "lidar_points.csv"
    shutil.copy(str(paths["targets"]), str(points))
    os.symlink(str(d / "hist" / "sub"), str(d / "x" / "L2"))
    os.symlink(str(d / "2024-10-23" / "sub"), str(d / "x" / "L1"))
    os.symlink(str(d / "surv" / "sub"), str(d / "x" / "L3"))
    up_hist = d / "x" / "L2" / ".." / "innocent.csv"
    up_held = d / "x" / "L1" / ".." / "targets.csv"
    up_surv = d / "x" / "L3" / ".." / "lidar_points.csv"
    canaries = {str(p) for p in (canary_t, canary_eo, canary_io, canary_uv, link, link2, in_dir, history, tif, asc,
                                 lidar, plain, mid, chain, raster, to_tif, points, up_hist, up_held, up_surv)}
    base = ["--cameras", "c1", "--io", paths["io"], "--eo", paths["eo"], "--targets", paths["targets"],
            "--picks", paths["picks"], "--starts", 1, "--samples", 20, "--no-io-test"]
    cases = (("--targets", canary_t, "targets named for 2024-10-23"), ("--eo", canary_eo, "an EO *20241023_EO*"),
             ("--io", canary_io, "an IO named 20241023"), ("--picks", canary_uv, "picks named 2024_10_23"),
             ("--targets", link, "a link with an innocent name to them"), ("--targets", link2, "a link to that link"),
             ("--targets", in_dir, "a file in a folder named for the day"),
             ("--targets", history, "calibration_history.csv"), ("--targets", tif, "a .tif as targets"),
             ("--targets", asc, "a .asc as targets"), ("--eo", lidar, "a lidar-fitted EO"),
             ("--targets", chain, "a chain of links whose middle link alone is named for the day"),
             ("--targets", to_tif, "a link with an innocent name to a .tif"),
             ("--targets", up_hist, "link/../innocent.csv reaching calibration_history.csv"),
             ("--targets", up_held, "link/../targets.csv reaching a folder named for the day"),
             ("--targets", up_surv, "link/../lidar_points.csv reaching a survey file"))
    out = d / "out"
    WATCH["opened"], WATCH["on"] = [], True
    try:
        results = []
        for flag, p, what in cases:
            argv = list(base)
            argv[argv.index(flag) + 1] = p
            rc, log = cli(["--out-dir", out] + argv)
            results.append((what, rc, "REFUSED" in log and "not opened" in log))
    finally:
        WATCH["on"] = False
    seen = [p for p in WATCH["opened"] if p in canaries or os.path.realpath(p) in {os.path.realpath(c) for c in canaries}]
    for what, rc, said in results:
        check(rc == 2 and said, f"{what}: exit 2, REFUSED by a guard ('not opened')")
    check(not seen and not out.exists(), f"none of the {len(cases)} was opened (audit hook: {len(seen)} opens of them), "
                                         f"nothing written")
    for p in (canary_t, link):
        try:
            CE.refuse_input(p)
            check(False, f"refuse_input({p.name}) refused")
        except GU.Refused:
            check(True, f"refuse_input({p.name}): Refused (a GU.Refused, as georectify_uncertainty's guard)")
    # the day's other spellings (names only: none of these files exists)
    held_names = ("GCPs_Oct_23_2024.csv", "gcp_oct_23_2024.csv", "October_23_2024.csv", "10-23-2024_gcp.csv",
                  "23Oct2024.csv", "23-Oct-2024.csv", "2024-Oct-23.csv", "23_10_2024_targets.csv", "Oct 23 2024.csv",
                  "2024_10_23.csv", "20241023.csv", "1729699200.Wed.Oct.23_16_00_00.GMT.2024.CACO05.c1.timex.jpg")
    missed = []
    for name in held_names:
        try:
            CE.refuse_input(d / "spellings" / name)
            missed.append(name)
        except GU.Refused:
            pass
    check(not missed, f"{len(held_names)} spellings of 23 Oct 2024 (Oct_23_2024, October_23_2024, 10-23-2024, "
                      f"23Oct2024, 2024-Oct-23, the station's frame names, ...): all refused"
                      + (f"; MISSED {missed}" if missed else ""))
    fine = [str(p) for p in paths.values()] + [str(d / n) for n in ("2025-10-23_targets.csv", "Oct_23_2025.csv",
                                                                    "10-23-2025.csv", "gcp_2024_10_24.csv",
                                                                    "1763136000.Fri.Nov.14_16_00_00.GMT.2025.CACO05."
                                                                    "c1.timex.jpg")]
    wrong = []
    for p in fine:
        try:
            CE.refuse_input(p)
        except GU.Refused:
            wrong.append(Path(p).name)
    check(not wrong, "the live inputs and other days (2025-10-23, Oct_23_2025, 24 Oct 2024, a frame of 14 Nov 2025) "
                     "are not refused" + (f"; REFUSED {wrong}" if wrong else ""))
    # outputs
    rc, log = cli(["--out-dir", CAL] + base)
    rc2, log2 = cli(["--out-dir", CAL / "sub"] + base)
    check(rc == 2 and rc2 == 2 and "calibration/" in log and not (CAL / "c1_eo_report.txt").exists()
          and not (CAL / "sub").exists(), "--out-dir calibration/ (or a folder in it): exit 2, nothing written there")
    keep = d / "keep"
    keep.mkdir()
    precious = keep / "c1_eo_report.txt"
    precious.write_bytes(b"precious\r\nbytes\n")
    rc, log = cli(["--out-dir", keep] + base)
    check(rc == 2 and "exists" in log and precious.read_bytes() == b"precious\r\nbytes\n"
          and sorted(p.name for p in keep.iterdir()) == ["c1_eo_report.txt"],
          "an existing output: exit 2, left byte for byte, and no other output written")
    archived = CE.archived_names()
    for name, why in (("X_EO-CV.yaml", "_EO-CV.yaml"), ("CACO05_c1_20251113_EO.yaml", "_EO.yaml"),
                      ("README.md", "a name in calibration/"),
                      (sorted(archived)[0] if archived else "none", "gcpfree_eo.ARCHIVED_SHA256")):
        try:
            CE.refuse_output(d / "free" / name)
            check(False, f"an output named {name} ({why}) refused")
        except CE.Refused:
            check(True, f"an output named {name} ({why}): refused")
    check(len(archived) == 6, f"gcpfree_eo.ARCHIVED_SHA256 read with ast: {len(archived)} names")
    for argv, what in ((["--out-dir", d / "u", "--io", paths["io"]], "--io with both cameras"),
                       (["--out-dir", d / "u", "--camera-position", 1, 2, 3, 0.03, 0.05], "--camera-position, both"),
                       (["--cameras", "c1"], "no --out-dir"),
                       (["--out-dir", d / "u", "--cameras", "c1", "--camera-position", 1, 2, 3, 0, 0.05],
                        "a zero sigma"),
                       (["--out-dir", d / "u", "--cameras", "c1", "--camera-position", "nan", 0, 0, 1, 1],
                        "a NaN position"),
                       (["--out-dir", d / "u", "--cameras", "c1", "--camera-position", 420084.27, 4638322.10, 20.51,
                         "nan", 0.05], "a NaN sigma"),
                       (["--out-dir", d / "u", "--cameras", "c1", "--camera-position", 420084.27, 4638322.10, "inf",
                         0.03, 0.05], "an infinite height"),
                       (["--out-dir", d / "u", "--cameras", "c1", "--seed", -1], "--seed -1"),
                       (["--out-dir", d / "u", "--cameras", "c1", "--seed", 2 ** 32], "--seed 2**32")):
        rc, log = cli(argv)
        check(rc == 6 and not (d / "u").exists() and "solving" not in log, f"usage error ({what}): exit 6")


# ---------------------------------------------------------------------
# 9: the real inputs
# ---------------------------------------------------------------------

# [REAL] The live run's numbers (docs/VALIDATION.md 4.6c; calibrate_eo.py --out-dir D, 2000 draws, 12 starts, seed 0):
# on z = 0 at the data's ranges, the cross-shore sigma (the Jacobian's, the larger), the model choice (all six
# minus the angles), the signed shift from the live EO (seaward +) and 1 px of pick origin, m; the angles' sds, deg
REAL_GROUND = {
    "c1": {"sigma": (0.172, 0.260, 0.368), "model": (0.009, 0.018, 0.032), "shift": (-0.001, -0.001, -0.002),
           "pick": (0.071, 0.108, 0.152)},
    "c2": {"sigma": (0.143, 0.310, 0.540, 0.835, 1.196), "model": (0.123, 0.619, 1.462, 2.646, 4.162),
           "shift": (0.099, 0.216, 0.372, 0.567, 0.801), "pick": (0.068, 0.149, 0.261, 0.406, 0.583)}}
REAL_SD = {"c1": {"jacobian": (0.052, 0.053, 0.184), "jackknife": (0.066, 0.026, 0.088)},
           "c2": {"jacobian": (0.034, 0.029, 0.175), "jackknife": (0.026, 0.027, 0.084)}}
# the lens verdicts: (flagged, best LOO RMS with the term or None, the targets without which the verdict flips)
REAL_LENS = {"c1": {"d1": (True, 0.400, [114, 115, 210, 211, 214, 216, 217]), "d1+d2": (True, 0.347, [214, 216])},
             "c2": {"d1": (False, 0.410, [212]), "d1+d2": (False, None, [317])}}


def near(got, want, rel, absolute):
    return got is not None and abs(got - want) <= rel * abs(want) + absolute


def test_real(d):
    print("9 [REAL] the live 2025-11-13 inputs, both cameras (--samples 300 --starts 4)")
    want = {"c1": {"angles": 0.453, "all": 0.589}, "c2": {"angles": 0.410, "all": 1.060}}
    before = {}
    for cam in ("c1", "c2"):
        paths = live(cam)[0]
        for p in list(paths.values()) + [CAL / f"{cam}_geom_sigma.csv"]:
            before[str(p)] = GU.sha256(p)
    out = Path(d) / "real"
    t0 = time.time()
    rc, log = cli(["--out-dir", out, "--samples", 300, "--starts", 4])
    check(rc == 0, f"the run: exit {rc} ({time.time() - t0:.0f} s)")
    if rc != 0:
        print(log[-2000:])
        return
    for cam in ("c1", "c2"):
        j = json.loads((out / f"{cam}_eo_report.json").read_text())
        a = j["candidates"]["angles"]
        mv = a["move_from_start_m_deg"]
        check(mv[0] <= 0.10 and mv[1] <= 0.10, f"{cam}: the angles solve reproduces the live EO (moves {mv[0]:.3f} m, "
                                               f"{mv[1]:.3f} deg)")
        for k, v in want[cam].items():
            got = j["candidates"][k]["loo_rms_m"]
            check(abs(got - v) <= 0.01, f"{cam}: leave-one-out ground RMS, {k}: {got:.3f} m (the design critic's {v})")
        check(j["chosen"] == "angles" and j["reproduces_start"], f"{cam}: --free auto picks {j['chosen']}")
        ms = [j["candidates"][k]["multistart"] for k in ("all", "angles")]
        check(all(m["n"] == 4 and m["one_minimum"] for m in ms), f"{cam}: every start reaches one minimum")
        e = load_extrinsics(str(out / f"CACO05_{cam}_20251113_EO-refit.yaml"))
        ref = j["eo_refit"]
        ytxt = (out / f"CACO05_{cam}_20251113_EO-refit.yaml").read_text()
        check(abs(e[0] - ref["x"]) < 1e-3 and abs(np.degrees(e[4]) - ref["tilt"]) < 1e-5
              and f"# {CE.NOT_ADOPTED}\n" in ytxt and "#roll - camera roll" in ytxt,
              f"{cam}: the refit EO reads back (georectify.load_extrinsics), with the NOT ADOPTED line and the live "
              f"file's comment lines")
        g = j["ground"]
        check(all(v is not None and v > 0 for v in g["sigma_cross_shore_m"][g["covariance_used"]].values())
              and len(g["ranges_m"]) == len(CE.DATA_RANGES_M[cam]),
              f"{cam}: cross-shore sigma at " + "/".join(str(r) for r in g["ranges_m"]) + " m: "
              + "/".join(f"{g['sigma_cross_shore_m'][g['covariance_used']][f'{r}|0.0']:.2f}" for r in g["ranges_m"])
              + f" m on z = 0 ({g['covariance_used']})")
        real_ground(cam, j)
        real_lens(cam, j)
        txt = (out / f"{cam}_eo_report.txt").read_text()
        note = ("these 14 targets of 2025-11-13 are an out-of-sample check" if cam == "c1"
                else "Solving all six parameters here moves the camera 0.33 m from that position")
        check(note in txt, f"{cam}: the report's note on the live EO's provenance: '...{note}...'")
    names = sorted(p.name for p in out.iterdir())
    check(names == sorted(["CACO05_c1_20251113_EO-refit.yaml", "CACO05_c2_20251113_EO-refit.yaml",
                           "c1_eo_report.json", "c1_eo_report.txt", "c2_eo_report.json", "c2_eo_report.txt"]),
          f"the out-dir holds the 6 outputs only (no IO written): {names}")
    after = {p: GU.sha256(p) for p in before}
    check(after == before, "the inputs and calibration/<cam>_geom_sigma.csv are unchanged (sha256)")
    keep = {p.name: p.read_bytes() for p in out.iterdir()}
    rc, log = cli(["--out-dir", out, "--samples", 300, "--starts", 4])
    check(rc == 2 and {p.name: p.read_bytes() for p in out.iterdir()} == keep,
          "run again into the same folder: exit 2, every output left as it was")
    paths, io, eo, nums, xyz, uv = live("c1")
    xp, up = write_targets(d, "six_c1", nums[:6], xyz[:6], uv[:6])
    out6 = Path(d) / "real_six"
    rc, log = cli(["--out-dir", out6, "--cameras", "c1", "--targets", xp, "--picks", up, "--starts", 1,
                   "--samples", 20, "--no-io-test"])
    txt = (out6 / "c1_eo_report.txt").read_text() if rc == 0 else ""
    check(rc == 0 and "6 targets," in txt and "The live EO's provenance is not settled" in txt
          and "14 targets" not in txt and "out-of-sample" not in txt,
          "6 of c1's targets with the live EO: the report says 6 targets, and not that 14 are an out-of-sample check")


def real_ground(cam, j):
    """The ground numbers of the live run against docs/VALIDATION.md 4.6c (z = 0), and their make-up."""
    g, w = j["ground"], REAL_GROUND[cam]
    rs = [int(r) for r in g["ranges_m"]]
    used = g["covariance_used"]

    def cells(table, z="0.0"):
        return [table.get(f"{r}|{z}") for r in rs]
    rows = (("sigma", cells(g["sigma_cross_shore_m"][used]), w["sigma"], 0.10, 0.01),
            ("model choice", cells(g["model_choice_all_minus_angles_m"]), w["model"], 0.03, 0.02),
            ("shift from the live EO", cells(g["shift_from_start_m"]), w["shift"], 0.03, 0.02),
            ("1 px of pick origin", cells(g["pick_origin"]["cross_shore_m_per_px"]), w["pick"], 0.0, 0.015))
    for what, got, want, rel, ab in rows:
        check(len(got) == len(want) and all(near(a, b, rel, ab) for a, b in zip(got, want)),
              f"{cam}: {what} at {'/'.join(map(str, rs))} m, z = 0: "
              + "/".join("-" if a is None else f"{a:+.3f}" for a in got) + " m (docs "
              + "/".join(f"{b:+.3f}" for b in want) + f"; within {100 * rel:.0f}% + {ab:g} m)")
    med = {k: float(np.nanmedian(np.array(v, float))) for k, v in g["per_node"]["sigma"].items()}
    z1 = cells(g["sigma_cross_shore_m"][used], "1.0")
    z1_txt = "/".join("-" if x is None else f"{x:.3f}" for x in z1)
    check(used == "jacobian" and med[used] == max(med.values())
          and all(b is not None and b > a for a, b in zip(cells(g["sigma_cross_shore_m"][used]), z1))
          and g["columns"] == [768, 1024, 1280, 1536, 1792] and g["z_planes"] == [0.0, 1.0],
          f"{cam}: the sigma is the larger covariance's ({used}, median {med[used]:.3f} against "
          f"{min(med.values()):.3f} m), larger on z = 1 ({z1_txt} m), over the middle columns {g['columns']}")
    a = j["candidates"]["angles"]
    for k in ("jacobian", "jackknife"):
        sd = a[f"sd_{k}"]
        got = [sd["azimuth"], sd["tilt"], sd["roll"]]
        check(all(abs(x - y) <= 0.002 for x, y in zip(got, REAL_SD[cam][k])) and set(sd) == {"azimuth", "tilt", "roll"},
              f"{cam}: the angles' sd ({k}): azimuth / tilt / roll " + "/".join(f"{x:.4f}" for x in got)
              + " deg (docs " + "/".join(f"{x:.3f}" for x in REAL_SD[cam][k]) + ", +/-0.002)")


def real_lens(cam, j):
    """The lens verdicts of the live run: every test drop-one tested, the flags and the targets that carry them."""
    ts = j["lens_tests"]
    by_terms = {k: [t for t in ts if "+".join(t["terms"]) == k] for k in ("d1", "d1+d2")}
    check(len(ts) == 4 and all(len(v) == 2 for v in by_terms.values()) and all(t["drop_one"] for t in ts)
          and all(v["drop_one"] for v in j["lens_verdicts"]),
          f"{cam}: the lens test ran (d1 and d1+d2, each with the angles and all six) and every test and verdict was "
          f"re-run with each target dropped in turn")
    for v in j["lens_verdicts"]:
        k = "+".join(v["terms"])
        flagged, best, carried = REAL_LENS[cam][k]
        check(v["flagged"] == flagged and (best is None or abs(v["with"] - best) <= 0.005)
              and v["carried_by"] == carried,
              f"{cam}: the {k} verdict: {'flagged' if v['flagged'] else 'not flagged'} ({v['without']:.3f} -> "
              f"{v['with']:.3f} m, {v['with_set']} + {k}), without any one of {v['carried_by']} it flips")
    if cam == "c1":
        t = next(t for t in ts if t["set"] == "all" and t["terms"] == ["d1"])
        check(abs(t["loo_rms_with"] - 0.400) <= 0.005 and abs(t["loo_rms_without"] - 0.589) <= 0.01,
              f"c1: all + d1, leave-one-out {t['loo_rms_without']:.3f} -> {t['loo_rms_with']:.3f} m (docs 0.589 -> "
              f"0.400)")


def test_pins():
    print("10 the pinned modules; Python 3.8 syntax")
    import gcpfree_eo as G
    prob = G.pin_problem()
    check(prob == "", "gcpfree_eo.pin_problem(): this checkout's gcpfree_eo.py and its MODULES are the version "
                      "survey_products.py pins" + (f" -- {prob}" if prob else ""))
    import ast
    for f in ("calibrate_eo.py", "test_calibrate_eo.py"):
        src = (HERE / f).read_text(encoding="utf-8")
        try:
            ast.parse(src, feature_version=(3, 8))
        except TypeError:                        # Python 3.7 has no feature_version
            ast.parse(src)
        check(True, f"{f} parses as Python 3.8")
    print("  (run test_gcpfree_eo.py and test_georectify_uncertainty.py separately: both must still pass)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="keep the outputs, in a new folder inside this one")
    args = ap.parse_args()
    if args.keep:
        Path(args.keep).mkdir(parents=True, exist_ok=True)
    sys.addaudithook(_audit)
    d = os.path.realpath(tempfile.mkdtemp(prefix="calibrate_eo_", dir=args.keep))
    print(f"test folder: {d}")
    t0 = time.time()
    try:
        test_recover()
        test_blunder()
        test_multistart(d)
        test_auto()
        test_prior_cli(d, *test_prior())
        test_lens()
        test_horizon(d)
        test_guards(d)
        test_real(d)
        test_pins()
    finally:
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    print(f"({time.time() - t0:.0f} s)")
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED")
        return 1
    print("all calibrate_eo tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

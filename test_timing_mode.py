#!/usr/bin/env python3
"""
Tests For The Water-Level Timing Switch (TIMING_MODE)
=======================================================
Plain python3 (no pytest needed on the station):

    python3 test_timing_mode.py [--keep DIR] [--before REF]

WHAT. The switch of timex_timing.py, off by default:
  * LEGACY IS BIT FOR BIT THE CODE WITHOUT THE SWITCH. The detector, the
    contour step and the DEM step are run on the same synthetic inputs by
    the code without the switch and by this code with no mode given, with
    TIMING_MODE=legacy, with --timing-mode legacy and (the detector) with
    no calibration folder at all (the built-in values): every detection
    CSV, overlay PNG, contour file, DEM grid and DEM info.json must be
    byte-identical. "Without the switch" is read with git (before_ref):
    where the switch came in is the commit that added timex_timing.py
    when it is on this checkout's own line (its branch, a fast-forward, a
    rebase), else the merge that joined that commit's line; this
    checkout's file with that commit's changes taken back out by git
    merge-file (against its first parent: exactly the tree before the
    switch, the conflict resolutions and every branch merged before it
    included), so a fix another commit made to the same file is in both
    runs and is not mistaken for the switch. Where that does not merge
    cleanly, the file as of that first parent; that copy, or one taken
    from a commit that added timex_timing.py but is not the switch commit
    as written (a squash, which carries other branches' changes too), is
    not exactly the code without the switch, so a difference is then a
    SKIP naming the --before command to decide it, never a failure.
    --before REF compares with REF as committed instead (REF a commit with
    every change of this checkout but the switch's; on the switch's own
    branch, deep/cam-timing, --before 42e3e86); before the switch is
    committed, HEAD. Without git that comparison is skipped (said so at
    the end), and the rest still runs: the legacy timex bias file holds
    exactly the built-in values, legacy levels and waves are those of the
    file-name epoch, the headers are as before, the blend still runs.
  * MID does all three things, together: the detector's blend is off and its
    timex bias is the mid file's (its rows = a legacy-code run with
    --no-temporal-filter and the same corrections, plus the Timing_Mode
    column); without a mid file a mid timex run is refused (for the cameras
    it has photos of whose bias is on: detect_original_view.py's one camera
    with its bias off still runs); the contour step
    takes the level, the waves and the stale-gap check at epoch + 300 s for
    a timex and + 0 s for a snap, keeps capture_time/capture_epoch, adds
    level_time_offset_s last, and leaves out every detection of the other
    mode (in both directions).
  * The DEM's info.json says when the levels were taken (only for such
    contours; nothing for legacy ones) and WARNs on a mix.
  * The one constant: runup_from_timestack.py, timestack_wave_period.py,
    relabel_manifests.py and reaim_figure.py read TIMEX_MID_OFFSET_S and
    keep their 300.0 s.
  * The cron: one TIMING_MODE line, legacy, passed to the detector and the
    contour step, its own archive for mid. process_historical.py and
    process_chelsea.py: one mode, the same to both steps;
    detect_original_view.py unchanged (pinned by the GCP-free
    pre-registration): its two steps inherit TIMING_MODE.
  * recalibrate_timing.py: the exposure check (file write times, timestack
    lengths) passes a start-named 600 s collection and stops on an end-named
    one; its rules give back the v6.8 corrections from the
    v6.8 numbers (c1 constant -11.03, c2's curve with its far field
    clamped); the pass rule; the mid file is written (never the legacy one)
    and read back by the detector; its derive step on synthetic detections;
    its compare step end to end on synthetic waterlines laid on the real
    CACO05 calibration (a 7.9 min lag planted in the legacy lines, none in
    the mid ones): PASS.

HOW. Synthetic photos (480 x 400), detection CSVs, a gnssrefl-format
spline and wave files in a temporary folder (--keep DIR to look at it).

STATION FILES. None are read; the repository's ground_truth/ lines of 11
Sep 2026 13:00 shape the synthetic waterlines of the compare test, and its
calibration/ georectifies them.
"""

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import timex_timing                               # noqa: E402

FAILED = []
SKIPPED = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILED.append(msg)


def skip(msg):
    print("  skip " + msg)
    SKIPPED.append(msg)


def run(args, cwd, env=None):
    e = dict(os.environ)
    e.pop("TIMING_MODE", None)                    # the test sets the mode itself
    # A script copied out of the repository (the code without the switch, a
    # detector without its calibration folder) still imports its companions
    # from here; its own folder comes first, so its own copies win.
    e["PYTHONPATH"] = str(HERE) + (os.pathsep + e["PYTHONPATH"] if e.get("PYTHONPATH") else "")
    e.update(env or {})
    r = subprocess.run([sys.executable] + [str(a) for a in args], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(cwd), env=e)
    return r.returncode, r.stdout


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def folder_digest(d):
    """{file name: sha256} of every file in d."""
    return {p.name: sha(p) for p in sorted(Path(d).iterdir()) if p.is_file()}


def stamp(epoch):
    t = datetime.fromtimestamp(epoch, tz=timezone.utc)
    return f"{epoch}.{t:%a.%b.%d_%H_%M_%S}.GMT.{t:%Y}.CACO05"


# ---------------------------------------------------------------------
# git: the code before the switch
# ---------------------------------------------------------------------

def git(*a):
    try:
        r = subprocess.run(["git", "-C", str(HERE)] + list(a), stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, universal_newlines=True)
    except OSError:
        return None
    return r.stdout if r.returncode == 0 else None


# The switch commit's subject as written: a commit that added timex_timing.py under another one
# (a squash) may carry other changes besides the switch's.
SWITCH_SUBJECT = "Water-level timing fix behind one switch, TIMING_MODE (off by default)"


def switch_commit():
    """The commit that added timex_timing.py (the timing switch), or None (not committed yet, no git)."""
    hashes = (git("log", "--diff-filter=A", "--format=%H", "--", "timex_timing.py") or "").split()
    return hashes[-1] if hashes else None


def has(commit, tip):
    return git("merge-base", "--is-ancestor", commit, tip) is not None


def oldest_with(commit, tip):
    """The oldest commit of tip's first-parent line that has `commit` (all newer ones have it too), or None."""
    line = (git("rev-list", "--first-parent", tip) or "").split()           # newest first
    if not line or not has(commit, line[0]):
        return None
    lo, hi = 0, len(line) - 1
    if has(commit, line[hi]):
        return line[hi]
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if has(commit, line[mid]) else (lo, mid)
    return line[lo]


def switch_entry(switch):
    """
    Where the switch came into this checkout: the switch commit itself when
    it is on HEAD's own first-parent line (its branch, a fast-forward, a
    rebase), else the innermost merge that joined its line (each merge on
    the way followed into the parent that brought it). That commit's first
    parent is the tree just before the switch. -> the commit, or None.
    """
    tip, entry = "HEAD", None
    while True:
        c = oldest_with(switch, tip)
        if c is None or c == switch:
            return entry if c else None
        entry = c
        parents = (git("rev-list", "--parents", "-n", "1", c) or "").split()[2:]
        tip = next((q for q in parents if has(switch, q)), None)
        if tip is None:
            return entry


def before_ref(override=None):
    """
    The code before the timing switch: (ref, entry, note). ref is
    --before, else the first parent of entry (switch_entry: the switch
    commit, or the merge that joined it), else HEAD (the switch not
    committed yet); None without git. entry is the commit whose changes
    old_copy takes out (None: ref as committed). note: why entry's changes
    may be more than the switch's (the switch commit not as written: a
    squash), else None.
    """
    if override:
        ok = git("rev-parse", "--verify", "-q", override + "^{commit}") is not None
        return (override if ok else None), None, None
    switch = switch_commit()
    entry = switch_entry(switch) if switch else None
    ref = entry + "^1" if entry else switch + "^" if switch else "HEAD"
    if git("rev-parse", "--verify", "-q", ref + "^{commit}") is None:
        return None, None, None
    entry = entry or switch
    note = None
    if switch and entry == switch and (git("log", "-1", "--format=%s", switch) or "").strip() != SWITCH_SUBJECT:
        note = (f"{switch[:7]}, which added timex_timing.py, is not the switch commit as written (a squash?), "
                f"so its changes may hold other branches' too")
    return ref, entry, note


def without_switch(now, at_switch, before, scratch):
    """`now` with the changes from `before` to `at_switch` taken out (git merge-file: a 3-way
    merge of at_switch -> before into now). -> the text, or None on a conflict or an error."""
    scratch.mkdir(parents=True, exist_ok=True)
    paths = [scratch / k for k in ("now", "at_switch", "before")]
    for p, text in zip(paths, (now, at_switch, before)):
        p.write_text(text)
    try:
        r = subprocess.run(["git", "merge-file", "-p"] + [str(p) for p in paths], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, universal_newlines=True)
    except OSError:
        return None
    return r.stdout if r.returncode == 0 else None


def old_copy(ref, names, dest, entry=None, note=None):
    """
    The files `names` without the timing switch, into dest. With `entry`
    (ref is then its first parent): each is THIS checkout's file with only
    entry's changes taken out (without_switch), so a later commit's change
    to the same file is in the comparison on both sides. Where a later
    change does not merge cleanly with the switch taken out (it touched the
    same lines), and for every file without `entry`: the file as of ref,
    as committed. -> None if one is not readable, else why the copy may
    not be exactly the code without the switch ('' when it is; with
    `note`, never '').
    """
    dest.mkdir(parents=True, exist_ok=True)
    why = []
    for n in names:
        text = git("show", f"{ref}:{n}")
        if text is None:
            return None
        at_entry = git("show", f"{entry}:{n}") if entry else None
        now = (HERE / n).read_text()
        if at_entry is not None and now != at_entry:
            merged = without_switch(now, at_entry, text, dest / ".merge" / n)
            if merged is None:
                why.append(f"{n} was changed after {entry[:7]} on lines the switch changed, so it is compared "
                           f"with {ref} as committed, which lacks that later change")
                print(f"  note: {why[-1]}")
            else:
                text = merged
        (dest / n).write_text(text)
    return "; ".join(why + ([note] if note else []))


def check_before(cond, msg, inexact):
    """check() of a comparison with the code without the switch; when that code may not be exact
    (old_copy's reason) and they differ, a skip: the difference may be another change."""
    if cond or not inexact:
        check(cond, msg)
    else:
        skip(f"{msg}: they DIFFER, but {inexact}; decide with python3 test_timing_mode.py --before REF, REF a "
             f"commit with every change of this checkout but the timing switch's (docs/VALIDATION.md 2.7)")


# ---------------------------------------------------------------------
# The shared definitions
# ---------------------------------------------------------------------

def test_definitions():
    print("timex_timing.py: the mode, the offsets, the files")
    tt = timex_timing
    check(tt.TIMEX_MID_OFFSET_S == 300 and tt.DEFAULT_TIMING_MODE == "legacy" and tt.TIMING_MODES == ("legacy", "mid"),
          "TIMEX_MID_OFFSET_S 300, modes legacy (default) and mid")
    check(tt.timing_mode(None, {}) == "legacy" and tt.timing_mode(None, {"TIMING_MODE": "mid"}) == "mid"
          and tt.timing_mode("legacy", {"TIMING_MODE": "mid"}) == "legacy"
          and tt.timing_mode(None, {"TIMING_MODE": ""}) == "legacy" and tt.timing_mode(" MID ", {}) == "mid",
          "mode: the option, else TIMING_MODE, else legacy")
    try:
        tt.timing_mode(None, {"TIMING_MODE": "midd"})
        bad = False
    except ValueError as e:
        bad = "not one of legacy, mid" in str(e)
    check(bad, "an unknown mode is an error, not a guess")
    t, s = "1790863200.Thu.Oct.01_14_00_00.GMT.2026.CACO05.c1.timex.csv", "x/1790863200.CACO05.c1.snap.csv"
    check(tt.level_time_offset_s(t, "legacy") == 0 == tt.level_time_offset_s(s, "legacy")
          and tt.level_time_offset_s(t, "mid") == 300 and tt.level_time_offset_s(s, "mid") == 0,
          "level time: legacy +0 s; mid +300 s for a timex, +0 s for a snap")
    check(tt.archive_suffix("legacy") == "" and tt.archive_suffix("mid") == "_mid",
          "archive: processed_timex (legacy), processed_timex_mid (mid)")
    check(tt.timex_bias_path("mid").name == "timex_bias_mid.json"
          and tt.timex_bias_path("legacy") == HERE / "calibration" / "timex_bias_legacy.json",
          "calibration/timex_bias_<mode>.json")


def test_legacy_bias_file():
    print("calibration/timex_bias_legacy.json: exactly the detector's built-in v6.8 values")
    import waterline_detector_v5 as det
    pts = timex_timing.read_timex_bias(HERE / "calibration" / "timex_bias_legacy.json", "legacy")
    same = all(pts[c] == det.CAMERAS[f"CACO05_{c.upper()}"].bias_correction_points_timex for c in ("c1", "c2"))
    check(same and pts["c1"] == ((0.0, -11.03), (1.0, -11.03)) and len(pts["c2"]) == 10
          and pts["c2"][-3:] == ((0.748, -17.22), (0.848, -17.22), (0.949, -17.22)),
          "identical (==, every float) to CameraProfile.bias_correction_points_timex: c1 -11.03, c2's 10 points")
    try:
        timex_timing.read_timex_bias(HERE / "calibration" / "timex_bias_legacy.json", "mid")
        refused = False
    except ValueError as e:
        refused = "not of 'mid'" in str(e)
    check(refused, "the legacy file is refused as a mid file")
    # (a station's own calibration/timex_bias_mid.json is never read by this test: every mid run
    # names its file, or runs from a copy without a calibration folder)


def test_shared_constant():
    print("TIMEX_MID_OFFSET_S: one definition, the four users unchanged")
    import runup_from_timestack as rft
    import timestack_wave_period as twp
    check(rft.BURST_MID_OFFSET_S == 300.0 and isinstance(rft.BURST_MID_OFFSET_S, float)
          and twp.BURST_MID_OFFSET_S == 300.0 and isinstance(twp.BURST_MID_OFFSET_S, float),
          "runup_from_timestack / timestack_wave_period BURST_MID_OFFSET_S: 300.0 (a float, as before)")
    for f, opt in (("relabel_manifests.py", "--offset"), ("reaim_figure.py", "--exposure-mid")):
        src = (HERE / f).read_text()
        check(f'"{opt}", type=float, default=float(TIMEX_MID_OFFSET_S)' in src and "default=300.0" not in src,
              f"{f} {opt}: default float(TIMEX_MID_OFFSET_S) = 300.0")
    import re
    defs = sorted(f.name for f in HERE.glob("*.py") if f.name != "test_timing_mode.py"
                  and re.search(r"^\s*TIMEX_MID_OFFSET_S\s*=", f.read_text(encoding="utf-8", errors="replace"), re.M))
    lits = sorted(f.name for f in HERE.glob("*.py")
                  if re.search(r"^\s*BURST_MID_OFFSET_S\s*=\s*300", f.read_text(encoding="utf-8", errors="replace"), re.M))
    check(defs == ["timex_timing.py"] and not lits, f"defined once: {defs}; no 300 s literal left: {lits}")
    import relabel_manifests as rm
    import pandas as pd
    adcp = pd.DataFrame({"epoch": [0, 3600], "h": [1.0, 2.0], "t": [8.0, 9.0]})
    df = pd.DataFrame({"filename": ["0.x.jpg"], "wave_height_m": [1.0], "wave_period_s": [8.0]})
    new, _, _ = rm.relabel(df, adcp, float(timex_timing.TIMEX_MID_OFFSET_S), 2.0)
    check(abs(new["wave_height_m"].iloc[0] - (1.0 + 300 / 3600)) < 1e-12
          and new["label_method"].iloc[0] == "adcp linear at epoch+300s",
          "relabel at the default offset: the height 300 s in, labelled 'epoch+300s' as before")


# ---------------------------------------------------------------------
# Detector
# ---------------------------------------------------------------------

W, H = 480, 400
T0 = int(datetime(2026, 10, 1, 14, 0, tzinfo=timezone.utc).timestamp())


def scene(cam, k, seed):
    """A synthetic timex: water, foam, wet and dry sand, the waterline inside the camera's search
    envelope, moving 3 (c1) / 2 (c2) px a frame."""
    import cv2
    rng = np.random.RandomState(seed)
    x = np.arange(W) / (W - 1)
    if cam == "c1":
        line = H * (0.26 + 0.6 * (0.10 + 0.03 * np.sin(2 * np.pi * x * 1.5))) + 3.0 * k
    else:
        line = H * (0.03 + 0.47 * (0.12 + 0.75 * x)) + 2.0 * k
    d = np.arange(H)[:, None] - line[None, :]
    img = np.zeros((H, W, 3), np.float32)
    img[:] = (150, 190, 215)                      # dry sand (BGR)
    img[d < 12] = (95, 120, 140)                  # wet sand
    img[d < 0] = (235, 235, 235)                  # foam
    img[d < -8] = (150, 110, 60)                  # water
    img += rng.normal(0, 6, img.shape).astype(np.float32)
    path = scene.out / f"{stamp(T0 + 1800 * k)}.{cam}.timex.jpg"
    ok = cv2.imwrite(str(path), np.clip(img, 0, 255).astype(np.uint8), [cv2.IMWRITE_JPEG_QUALITY, 95])
    assert ok
    # the file's time = its name's, as on the station (the detector WARNs on a mismatch over a day)
    os.utime(str(path), (T0 + 1800 * k, T0 + 1800 * k))


def detect(script, d, name, extra=(), env=None, cwd=None):
    """Runs a detector on the synthetic photos into d/name. -> (exit code, log, output folder)."""
    out = d / name
    shutil.rmtree(str(out), ignore_errors=True)
    rc, log = run([script, "--image-suffix", "timex.jpg", "--source-dir", d / "photos",
                   "--input-dir", out / "in", "--output-dir", out / "det", "--debug-dir", out / "debug",
                   "--min-confidence", "0", "--min-signal-fraction", "0"] + list(extra),
                  cwd or d, env)
    (d / f"{name}.log").write_text(log)
    return rc, log, out / "det"


def test_detector(d, ref, mid_file, entry=None, note=None):
    print("waterline_detector_v5.py: legacy is the old detector bit for bit; mid = blend off + mid bias")
    d = Path(d) / "detector"
    (d / "photos").mkdir(parents=True, exist_ok=True)
    scene.out = d / "photos"
    for k in range(4):
        for i, cam in enumerate(("c1", "c2")):
            scene(cam, k, 10 * k + i)
    new = HERE / "waterline_detector_v5.py"
    runs = {}
    for name, extra, env in (("default", (), None), ("env_legacy", (), {"TIMING_MODE": "legacy"}),
                             ("flag_legacy", ("--timing-mode", "legacy"), {"TIMING_MODE": "mid"})):
        rc, log, out = detect(new, d, name, extra, env)
        runs[name] = (rc, log, folder_digest(out) if out.is_dir() else {})
    rc, log, dig = runs["default"]
    csvs = sorted(n for n in dig if n.endswith(".csv"))
    check(rc == 0 and len(csvs) == 8 and sum(n.endswith(".png") for n in dig) == 8,
          f"default run: 8 detection CSVs and 8 overlays ({len(csvs)} CSVs)")
    check("Timing mode  : legacy" in log and "Temporal blend: on" in log and "timex_bias_legacy.json" in log
          and "WARNING" not in log,
          "its banner: legacy, blend on, the legacy timex bias file; no WARNING")
    check(runs["env_legacy"][2] == dig and runs["flag_legacy"][2] == dig and runs["flag_legacy"][0] == 0,
          "TIMING_MODE=legacy and --timing-mode legacy (over TIMING_MODE=mid): byte-identical outputs")
    header = (d / "default" / "det" / csvs[0]).read_text().splitlines()[0]
    check(header == "Column,Row,Confidence,Row_Precise,Peak_Sharpness,Has_Signal",
          "legacy CSV header as before (no Timing_Mode column)")
    # no calibration folder: the built-in values
    nocal = d / "nocal"
    nocal.mkdir(exist_ok=True)
    for f in ("waterline_detector_v5.py", "timex_timing.py"):
        shutil.copy2(str(HERE / f), str(nocal / f))
    rc, log, out = detect(nocal / "waterline_detector_v5.py", d, "builtin")
    check(rc == 0 and folder_digest(out) == dig and "built-in (v6.8 values" in log,
          "no calibration/timex_bias_legacy.json: the built-in values, byte-identical outputs")
    # an unreadable legacy file: a WARNING and the built-in values, never a stopped station
    broken = d / "broken_code"            # not d / <run name>: detect() empties that
    (broken / "calibration").mkdir(parents=True, exist_ok=True)
    for f in ("waterline_detector_v5.py", "timex_timing.py"):
        shutil.copy2(str(HERE / f), str(broken / f))
    (broken / "calibration" / "timex_bias_legacy.json").write_text('{"mode": "legacy", "c1": [[0.0, -11.03]')
    rc, log, out = detect(broken / "waterline_detector_v5.py", d, "broken")
    check(rc == 0 and folder_digest(out) == dig and "WARNING: timex bias file unreadable" in log
          and "unreadable)" in log, "an unreadable legacy file: a WARNING, the built-in values, byte-identical outputs")
    # the blend still runs in legacy mode
    rc, log, out = detect(new, d, "noblend", ("--no-temporal-filter",))
    check(rc == 0 and folder_digest(out) != dig and sum(folder_digest(out)[n] != dig[n] for n in csvs) >= 6,
          "the default still blends: --no-temporal-filter changes the lines after the first of each camera")
    # the code before the switch
    if ref is None:
        skip("detector before/after: git cannot read the commit before the switch (not a checkout?)")
    else:
        old = d / "old"
        inexact = old_copy(ref, ["waterline_detector_v5.py"], old, entry, note)
        if inexact is not None:
            rc, log, out = detect(old / "waterline_detector_v5.py", d, "before")
            got = folder_digest(out) if out.is_dir() else {}
            check_before(rc == 0 and got == dig,
                         f"the detector without the switch: byte-identical CSVs and overlays "
                         f"({sum(got.get(n) == h for n, h in dig.items())}/{len(dig)} files equal)", inexact)
        else:
            skip(f"detector before/after: {ref}:waterline_detector_v5.py not readable")
    # mid: refused without a mid file
    rc, log, out = detect(new, d, "mid_nofile", ("--timing-mode", "mid", "--timex-bias-file", d / "missing.json"))
    check(rc != 0 and "not found" in log, "mid with a missing --timex-bias-file: refused")
    rc, log, out = detect(nocal / "waterline_detector_v5.py", d, "mid_nocal", env={"TIMING_MODE": "mid"})
    check(rc == 2 and "needs its own timex bias corrections" in log and not out.exists(),
          "TIMING_MODE=mid without calibration/timex_bias_mid.json: refused (exit 2), nothing detected")
    rc, log, out = detect(nocal / "waterline_detector_v5.py", d, "mid_nobias", ("--no-bias-correction",),
                          env={"TIMING_MODE": "mid"})
    check(rc == 0 and len(list(out.glob("*.csv"))) == 8, "... but runs with --no-bias-correction (the derive step)")
    # detect_original_view.py's case: one camera's photos, its bias off by --profile-json
    one = d / "photos_c2"
    one.mkdir(exist_ok=True)
    for p in (d / "photos").glob("*.c2.timex.jpg"):
        shutil.copy2(str(p), str(one / p.name))
    (d / "c2_nobias.json").write_text(json.dumps({"c2": {"no_bias": True}}))
    rc, log, out = detect(nocal / "waterline_detector_v5.py", d, "mid_onecam",
                          ("--profile-json", d / "c2_nobias.json", "--source-dir", one), env={"TIMING_MODE": "mid"})
    check(rc == 0 and len(list(out.glob("*.csv"))) == 4,
          "... and with only c2 photos, c2's bias off (--profile-json no_bias, as detect_original_view.py): runs")
    rc, log, out = detect(nocal / "waterline_detector_v5.py", d, "mid_twocam",
                          ("--profile-json", d / "c2_nobias.json"), env={"TIMING_MODE": "mid"})
    check(rc == 2 and "(CACO05_C1)" in log, "... but c1 photos with c2's bias off: refused for c1")
    rc, log, out = detect(new, d, "mid_legacyfile", ("--timing-mode", "mid", "--timex-bias-file",
                                                     HERE / "calibration" / "timex_bias_legacy.json"))
    check(rc == 2 and "not of 'mid'" in log, "mid given the legacy file: refused")
    # mid with the mid file: blend off, the file's corrections, Timing_Mode=mid on every row
    rc, log, out = detect(new, d, "mid", ("--timing-mode", "mid", "--timex-bias-file", mid_file))
    pts = timex_timing.read_timex_bias(mid_file, "mid")
    (d / "same_points.json").write_text(json.dumps({c: [list(p) for p in pts[c]] for c in ("c1", "c2")}))
    rc2, log2, out2 = detect(new, d, "same", ("--no-temporal-filter", "--bias-correction-file", d / "same_points.json"))
    equal, marked = True, True
    for n in csvs:
        a = (out / n).read_text().splitlines()
        b = (out2 / n).read_text().splitlines()
        marked &= a[0].endswith(",Timing_Mode") and all(r.endswith(",mid") for r in a[1:])
        equal &= [r.rsplit(",", 1)[0] for r in a] == b
    check(rc == 0 and rc2 == 0 and marked and equal and "Temporal blend: off" in log
          and f"Timex bias   : {mid_file}" in log,
          "mid: every row Timing_Mode=mid, and the rows those of a run with --no-temporal-filter and the "
          "mid file's corrections")


# ---------------------------------------------------------------------
# Contour step
# ---------------------------------------------------------------------

TIDE_T0 = int(datetime(2026, 10, 1, 11, 0, tzinfo=timezone.utc).timestamp())
SPLINE_END = TIDE_T0 + 14 * 3600


def tide(t):
    return 0.3 + 1.1 * np.sin(2 * np.pi * (np.asarray(t, float) - TIDE_T0) / (12.42 * 3600))


def write_spline(path, t_start, t_end, step=1800):
    """gnssrefl subdaily spline_out.txt layout, Hortho 19.014 (no datum shift). -> (epochs, levels)."""
    ts = np.arange(t_start, t_end + 1, step)
    zs = np.round(tide(ts), 3)
    with open(path, "w") as f:
        f.write("% station usgs gnssrefl v4.1.5\n")
        f.write("% orthometric height minus RH, where Hortho (m) is    19.014 \n")
        f.write("% MJD RH(m) YYYY MM DD HH MM SS quasi-sea-level(m)\n")
        for t, z in zip(ts, zs):
            u = datetime.fromtimestamp(int(t), tz=timezone.utc)
            f.write(f"{t / 86400.0 + 40587.0:15.7f} {19.014 - z:10.3f} {u.year:4d} {u.month:3d} {u.day:3d} "
                    f"{u.hour:3d} {u.minute:3d} {u.second:3d} {z:10.3f}\n")
    return ts, zs


def write_detection(path, row_at, mode=None, seed=0):
    rng = np.random.RandomState(seed)
    lines = ["Column,Row,Confidence,Row_Precise,Peak_Sharpness,Has_Signal" + (",Timing_Mode" if mode else "")]
    for c in range(400):
        r = row_at(c) + rng.normal(0, 0.6)
        lines.append(f"{c},{int(round(r))},0.95,{r:.3f},0.5,1" + (f",{mode}" if mode else ""))
    path.write_text("\n".join(lines) + "\n")


def detection_set(folder, mode=None):
    """c1 and c2 timex frames every 30 min, 11:00-20:30 UTC, whose rows follow the tide; two c1 snaps;
    one c1 timex 58 min after the spline's last reading (stale at mid-exposure only)."""
    folder.mkdir(parents=True, exist_ok=True)
    names = []
    frames = [(TIDE_T0 + 1800 * k, cam, "timex") for k in range(20) for cam in ("c1", "c2")]
    frames += [(TIDE_T0 + 1800 * k, "c1", "snap") for k in (3, 4)]
    frames += [(SPLINE_END + 58 * 60, "c1", "timex")]
    for i, (ep, cam, prod) in enumerate(frames):
        z = float(tide(ep + timex_timing.level_time_offset_s(f".{prod}.", mode or "legacy")))
        base = 500.0 + 40.0 * z + (30 if cam == "c2" else 0)
        name = f"{stamp(ep)}.{cam}.{prod}.csv"
        write_detection(folder / name, lambda c: base + 3.0 * np.sin(c / 7.0), mode, seed=i)
        names.append(name)
    return names


def test_extract(d, ref, entry=None, note=None):
    print("extract_elevation_contours.py: legacy as before; mid takes level, waves and stale check at mid-exposure")
    d = Path(d) / "extract"
    d.mkdir(parents=True, exist_ok=True)
    sp_t, sp_z = write_spline(d / "spline.txt", TIDE_T0 - 2 * 3600, SPLINE_END)
    hours = list(range(TIDE_T0 - 3600, SPLINE_END + 4 * 3600, 3600))
    (d / "waves_marconi.csv").write_text("time_utc,epoch,hs_best,hs_sigma,tp_s\n" + "".join(
        f"x,{e},{1.0 + 0.01 * k:.2f},0.1,8.0\n" for k, e in enumerate(hours)))
    detection_set(d / "legacy_det")
    detection_set(d / "mid_det", mode="mid")
    common = ["--waves", d / "waves_marconi.csv", "--setup-coef", "0.037"]

    def extract(script, det, out, extra=(), env=None):
        rc, log = run([script, d / "spline.txt", "--processed-dir", d / det, "--output", d / out] + common
                      + list(extra), d, env)
        (d / f"{out}.log").write_text(log)
        rows = list(csv.DictReader(open(d / out))) if rc == 0 and (d / out).exists() else []
        return rc, log, rows

    new = HERE / "extract_elevation_contours.py"
    rc, log, rows = extract(new, "legacy_det", "default.csv")
    outs = {"default": sha(d / "default.csv") if rc == 0 else None}
    for name, extra, env in (("env_legacy", (), {"TIMING_MODE": "legacy"}),
                             ("flag_legacy", ("--timing-mode", "legacy"), {"TIMING_MODE": "mid"})):
        rc2, _, _ = extract(new, "legacy_det", f"{name}.csv", extra, env)
        outs[name] = sha(d / f"{name}.csv") if rc2 == 0 else None
    frames = {r["source_file"] for r in rows}
    check(rc == 0 and len(frames) == 43 and outs["default"] == outs["env_legacy"] == outs["flag_legacy"],
          f"legacy: all 43 frames; default, TIMING_MODE=legacy, --timing-mode legacy byte-identical")
    head = list(rows[0].keys()) if rows else []
    check(head == ["source_file", "camera", "capture_time_utc", "capture_epoch", "pixel_column", "pixel_row",
                   "tide_elevation_navd88", "tide_model_spread", "tide_gap_minutes", "tide_is_stale",
                   "peak_sharpness", "water_level_source", "offshore_hs_m", "offshore_tp_s",
                   "setup_correction_m", "beach_elevation_navd88", "wave_source"],
          "legacy header as before (wave_source last, no level_time_offset_s)")
    first = {}
    for r in rows:
        first.setdefault(r["source_file"], r)
    lev_ok = all(abs(float(r["tide_elevation_navd88"]) - np.interp(int(r["capture_epoch"]), sp_t, sp_z)) < 1.5e-4
                 for r in first.values())
    hour_hs = {e: 1.0 + 0.01 * k for k, e in enumerate(hours)}

    def hs_at(t):             # the nearest hourly record, the earlier on a tie (nearest_within)
        best = min(hours, key=lambda e: (abs(e - t), e))
        return hour_hs[best]
    hs_ok = all(abs(float(r["offshore_hs_m"]) - hs_at(int(r["capture_epoch"]))) < 1e-9 for r in first.values())
    check(lev_ok and hs_ok, "legacy: each frame's level and Hs those at its file-name epoch")
    inexact = None if ref is None else old_copy(ref, ["extract_elevation_contours.py"], d / "old", entry, note)
    if ref is None:
        skip("contour step before/after: git cannot read the commit before the switch")
    elif inexact is not None:
        rc_o, log_o, _ = extract(d / "old" / "extract_elevation_contours.py", "legacy_det", "before.csv")
        check_before(rc_o == 0 and sha(d / "before.csv") == outs["default"],
                     "the contour step without the switch: byte-identical contour file", inexact)
    else:
        skip(f"contour step before/after: {ref}:extract_elevation_contours.py not readable")

    # mid
    rc, log, rows = extract(new, "mid_det", "mid.csv", ("--timing-mode", "mid"))
    first = {}
    for r in rows:
        first.setdefault(r["source_file"], r)
    stale = f"{stamp(SPLINE_END + 58 * 60)}.c1.timex"
    check(rc == 0 and len(first) == 42 and stale not in first and "Skipped (tide gap > 60 min): 1" in log,
          "mid: the frame 58 min past the spline's end is stale at mid-exposure (63 min) and skipped; 42 kept")
    off_ok = all(r["level_time_offset_s"] == ("300" if ".timex" in f else "0") for f, r in first.items())
    lev_ok = all(abs(float(r["tide_elevation_navd88"])
                     - np.interp(int(r["capture_epoch"]) + int(r["level_time_offset_s"]), sp_t, sp_z)) < 1.5e-4
                 for r in first.values())
    hs_ok = all(abs(float(r["offshore_hs_m"]) - hs_at(int(r["capture_epoch"]) + int(r["level_time_offset_s"])))
                < 1e-9 for r in first.values())
    moved = sum(1 for f, r in first.items() if ":30:00" in r["capture_time_utc"] and ".timex" in f)
    leg = {r["source_file"]: r for r in csv.DictReader(open(d / "default.csv"))}
    hs_moved = sum(1 for f, r in first.items() if ":30:00" in r["capture_time_utc"] and ".timex" in f
                   and abs(float(r["offshore_hs_m"]) - float(leg[f]["offshore_hs_m"]) - 0.01) < 1e-9)
    times_ok = all(r["capture_epoch"] == leg[f]["capture_epoch"] and r["capture_time_utc"] == leg[f]["capture_time_utc"]
                   for f, r in first.items())
    check(off_ok and lev_ok and times_ok, "mid: level at epoch + 300 s (timex) / + 0 s (snap); capture time and "
          "epoch unchanged; level_time_offset_s says which")
    check(hs_ok and moved > 0 and hs_moved == moved,
          f"mid: the waves too (a :30 timex takes the next hour's record at :35: {hs_moved}/{moved})")
    check(list(rows[0].keys())[-2:] == ["wave_source", "level_time_offset_s"] if rows else False,
          "mid: level_time_offset_s is the last column, after wave_source")
    # the other mode's detections are left out, both ways
    rc, log, _ = extract(new, "legacy_det", "mid_on_legacy.csv", ("--timing-mode", "mid"))
    check(rc == 1 and "WARNING: 43 detection(s) made in timing mode 'legacy' left out" in log,
          "mid on legacy (blended) detections: every one left out with a WARNING, no contours")
    rc, log, _ = extract(new, "mid_det", "legacy_on_mid.csv")
    check(rc == 1 and "WARNING: 43 detection(s) made in timing mode 'mid' left out" in log,
          "legacy on mid detections: every one left out with a WARNING")
    mixed = d / "mixed_det"
    shutil.copytree(str(d / "legacy_det"), str(mixed))
    for n in sorted((d / "mid_det").glob("*c2.timex.csv"))[:5]:
        shutil.copy2(str(n), str(mixed / n.name))
    rc, log, rows = extract(new, "mixed_det", "mixed.csv")
    check(rc == 0 and len({r["source_file"] for r in rows}) == 38
          and "Skipped (detected in another timing mode than legacy): 5 mid" in log,
          "a folder of both: legacy uses its 38, counts the 5 mid ones out")


# ---------------------------------------------------------------------
# DEM
# ---------------------------------------------------------------------

def test_dem(d, ref, entry=None, note=None):
    print("dem_from_contours.py: legacy as before; info.json water_level_timing only for timed contours")
    d = Path(d) / "dem"
    d.mkdir(parents=True, exist_ok=True)
    import test_dem_diagnostics as tdd
    rows = tdd.ground_rows(tdd.tidal_frames(days=5), "marconi:hs_best/tp_s")

    def build(name, offset_of):
        lines = [rows[0] + ("" if offset_of is None else ",level_time_offset_s")]
        for r in rows[1:]:
            lines.append(r if offset_of is None else r + "," + offset_of(r.split(",")[0]))
        (d / f"{name}.csv").write_text("\n".join(lines) + "\n")
        rc, log = run([HERE / "dem_from_contours.py", d / f"{name}.csv", d / name, "--no-plot"], d)
        info = json.loads((d / f"{name}_info.json").read_text()) if rc == 0 else {}
        return rc, log, info

    rc0, log0, info0 = build("plain", None)
    rc1, log1, info1 = build("mid", lambda f: "300")
    rc2, log2, info2 = build("mixed", lambda f: "300" if f.endswith(".c1.timex") else "")
    nf = len({r.split(",")[0] for r in rows[1:]})
    check(rc0 == 0 and "water_level_timing" not in info0 and "Water-level time" not in log0,
          "legacy contours: no water_level_timing in info.json, nothing in the log")
    w = info1.get("water_level_timing") or {}
    check(rc1 == 0 and w.get("frames") == {"+300 s": nf} and w.get("mixed") is False and "WARNING" not in log1
          and "Water-level time  : +300 s" in log1,
          f"mid contours: {nf} frames at +300 s, not mixed, no WARNING")
    w2 = info2.get("water_level_timing") or {}
    check(rc2 == 0 and w2.get("mixed") is True and set(w2.get("frames", {})) == {"+300 s", "legacy"}
          and any("WARNING" in ln and "mixes two water-level timings" in ln for ln in log2.splitlines()),
          "a mix of legacy and mid frames: mixed, with a WARNING")
    same = all(np.array_equal(np.loadtxt(d / f"plain_{k}.asc", skiprows=6), np.loadtxt(d / f"mid_{k}.asc", skiprows=6))
               for k in ("dem", "spread", "count"))
    check(same, "the grids do not depend on the column")
    # the DEM step before the switch, on the same legacy contours: every file byte-identical
    inexact = None if ref is None else old_copy(ref, ["dem_from_contours.py"], d / "old", entry, note)
    if ref is None:
        skip("DEM step before/after: git cannot read the commit before the switch")
    elif inexact is not None:
        outs = {}
        for name, script in (("now", HERE / "dem_from_contours.py"), ("before", d / "old" / "dem_from_contours.py")):
            (d / name).mkdir(exist_ok=True)
            # the old copy imports its companions (extract_elevation_contours, waterline_track) from
            # here (run's PYTHONPATH)
            rc, log = run([script, d / "plain.csv", d / name / "dem", "--no-plot", "--max-day-offset", "0.15"], d)
            (d / f"{name}.log").write_text(log)
            outs[name] = (rc, folder_digest(d / name), log.replace(str(d / name), "<out>"))
        a, b = outs["now"], outs["before"]
        check_before(a[0] == 0 and b[0] == 0 and len(a[1]) == 5 and a[1] == b[1] and a[2] == b[2],
                     f"the DEM step without the switch: byte-identical grids, sigma and info.json, and log "
                     f"({sum(b[1].get(n) == h for n, h in a[1].items())}/{len(a[1])} files equal)", inexact)
    else:
        skip(f"DEM step before/after: {ref}:dem_from_contours.py not readable")


# ---------------------------------------------------------------------
# Cron, historical paths
# ---------------------------------------------------------------------

def test_cron():
    print("waterline_timex_cron.sh: one switch, legacy")
    text = (HERE / "waterline_timex_cron.sh").read_text()
    assigns = [ln for ln in text.splitlines() if ln.startswith("TIMING_MODE=")]
    check(assigns == ["TIMING_MODE=legacy"], f"one TIMING_MODE line, legacy (the default not flipped): {assigns}")
    det = text[text.index("python3 waterline_detector_v5.py"):]
    det = det[:det.index(">>")]
    ext = text[text.index('python3 "$BASE/extract_elevation_contours.py"'):]
    ext = ext[:ext.index(">>")]
    check('--timing-mode "$TIMING_MODE"' in det and '--timing-mode "$TIMING_MODE"' in ext
          and "export TIMING_MODE" in text,
          "the detector and the contour step both take --timing-mode \"$TIMING_MODE\"; exported for the rest")
    check('mid) ARCHIVE_CSV="$BASE/archive/processed_timex_mid"' in text
          and '*)   ARCHIVE_CSV="$BASE/archive/processed_timex"' in text,
          "legacy keeps archive/processed_timex; mid has archive/processed_timex_mid")
    r = subprocess.run(["bash", "-n", str(HERE / "waterline_timex_cron.sh")], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT)
    check(r.returncode == 0, "bash -n")
    for f in ("process_historical.py", "process_chelsea.py"):
        src = (HERE / f).read_text()
        check(src.count('"--timing-mode", timing') == 2 and "timex_timing.timing_mode(args.timing_mode)" in src,
              f"{f}: one mode (--timing-mode, else TIMING_MODE), the same to the detector and the contour step")
    # detect_original_view.py is a module of the pre-registered GCP-free procedure (survey_products.py
    # GCPFREE_PINS, by sha256), so it is not edited: its detector and contour step are subprocesses that
    # inherit the environment, and so read the same TIMING_MODE
    import survey_products as sp
    src = (HERE / "detect_original_view.py").read_text()
    check(sha(HERE / "detect_original_view.py") == sp.GCPFREE_PINS[-1]["modules"]["detect_original_view.py"]
          and "env=" not in src and "TIMING_MODE" not in src,
          "detect_original_view.py unchanged (GCP-free pin); its detector and contour step inherit TIMING_MODE")


def test_station_status(d):
    print("station_status.py: the detections of the cron's mode")
    d = Path(d) / "status"
    for sub in ("archive/processed_timex", "archive/processed_timex_mid"):
        (d / sub).mkdir(parents=True, exist_ok=True)
    old = d / "archive/processed_timex/a.csv"
    old.write_text("x\n")
    os.utime(str(old), (1.0e9, 1.0e9))
    (d / "archive/processed_timex_mid/b.csv").write_text("x\n")
    lines = {}
    for mode in ("", "legacy", "mid", "cron_mid"):
        if mode == "cron_mid":
            (d / "waterline_timex_cron.sh").write_text("BASE=x\nTIMING_MODE=mid\n")
        env = {"TIMING_MODE": mode} if mode in ("legacy", "mid") else None
        rc, out = run([HERE / "station_status.py", "--base", d], d, env)
        lines[mode] = next((ln for ln in out.splitlines() if "detections" in ln), "")
    check("ALERT" in lines["legacy"] and "OK" in lines["mid"] and "processed_timex_mid" in lines["mid"]
          and "ALERT" in lines[""] and "OK" in lines["cron_mid"] and "processed_timex_mid" in lines["cron_mid"],
          "TIMING_MODE=legacy reads the legacy archive (old: ALERT), mid the mid one (new: OK); by hand the "
          "cron script's mode, never the newest archive (no script: legacy, ALERT; TIMING_MODE=mid in it: OK)")


# ---------------------------------------------------------------------
# recalibrate_timing.py
# ---------------------------------------------------------------------

V68_C2 = ((0.050, 3.87), (0.150, -1.19), (0.249, -8.15), (0.349, -8.49), (0.449, -5.26), (0.548, -4.55),
          (0.648, -8.44), (0.748, -17.22), (0.848, -41.28), (0.949, -57.78))
# c2's per-frame std: ~18 px everywhere (SEM 5.1 px over 13 frames), so above |mean| in every
# segment but the far two's 48.26 / 67.00 too: a rule that clamps by the scatter cascades on it
V68_C2_STD = (17.0, 18.0, 18.4, 18.4, 18.5, 18.6, 19.0, 18.6, 48.26, 67.00)


def table(camera, xs, raw, std, n, n_frames):
    return {"camera": camera, "n_frames": n_frames,
            "segments": [{"x": x, "raw": r, "std": s, "n": k} for x, r, s, k in zip(xs, raw, std, n)]}


def test_rules():
    print("recalibrate_timing.points_from_table: the v6.8 rules give back the v6.8 corrections")
    import recalibrate_timing as rt
    import waterline_detector_v5 as det
    c2_pts, why2 = rt.points_from_table(table("c2", [p[0] for p in V68_C2], [p[1] for p in V68_C2],
                                              V68_C2_STD, [13] * 10, 13))
    check(c2_pts == det.CAMERAS["CACO05_C2"].bias_correction_points_timex and "PER-COLUMN" in why2
          and "spread/SEM 12.0" in why2 and "0.848" in why2 and "0.949" in why2 and "CLAMPED to -17.22" in why2,
          f"c2's v6.8 numbers (spread 61.6 px, SEM 5.1): the per-column curve, last two clamped to -17.22 -- {why2}")
    xs = [p[0] for p in V68_C2]
    c1_raw = [-11.03 + o for o in (-5.6, -3.0, -1.0, 0.0, 1.0, 2.0, 3.0, 5.6, -1.0, -1.0)]
    c1_pts, why1 = rt.points_from_table(table("c1", xs, c1_raw, [14.5] * 10, [12] * 10, 12))
    check(c1_pts == det.CAMERAS["CACO05_C1"].bias_correction_points_timex and "CONSTANT -11.03" in why1
          and "spread/SEM 2.7" in why1,
          f"c1's (spread 11.2 px, SEM 4.2): constant -11.03 -- {why1}")
    pts, _ = rt.points_from_table(table("c2", xs, [p[1] for p in V68_C2], V68_C2_STD,
                                        [13] * 8 + [2, 0], 13))
    check(len(pts) == 9 and pts[-1] == (0.848, -17.22), "a segment no frame reaches is left out; one too few reach clamped")
    # mid-like means (+2.6 px): the same curve shifted, the far two clamped to the 0.748 value, nothing inside
    mid_raw = [p[1] + 2.6 for p in V68_C2]
    pts, why = rt.points_from_table(table("c2", xs, mid_raw, V68_C2_STD, [13] * 10, 13))
    check(pts is not None and [p[1] for p in pts[:8]] == [round(v, 2) for v in mid_raw[:8]]
          and pts[8][1] == pts[9][1] == round(mid_raw[7], 2) and "nothing clamped" not in why,
          f"mid-like c2 means with the real ~18 px std: the measured curve, only x > 0.8 clamped -- {why}")
    # the clamp does not depend on the sample: 300 draws of 13 frames at the real std, mid-like means
    rng = np.random.RandomState(3)
    same = 0
    for _ in range(300):
        draw = [m + rng.normal(0, sd / np.sqrt(13)) for m, sd in zip(mid_raw, V68_C2_STD)]
        sd = [abs(v + rng.normal(0, 2.0)) for v in V68_C2_STD]
        p2, _ = rt.points_from_table(table("c2", xs, draw, sd, [13] * 10, 13))
        same += p2 is not None and len(p2) == 10 and [q[1] for q in p2[:8]] == [round(v, 2) for v in draw[:8]] \
            and p2[8][1] == p2[9][1] == round(draw[7], 2)
    check(same == 300, f"300 draws at the real scatter: always only the far two clamped, never a cascade ({same}/300)")
    # what the rules cannot settle: nothing, and the reason
    pts, why = rt.points_from_table(table("c2", xs, [p[1] for p in V68_C2], V68_C2_STD, [13] * 6 + [4] + [13] * 3, 13))
    check(pts is None and "NEEDS A DECISION" in why and "0.648 (4 frames)" in why,
          f"a segment at x <= 0.8 reached by too few frames: no curve, a person decides -- {why}")
    xs20 = [round(0.025 + 0.05 * i, 3) for i in range(20)]
    pts, why = rt.points_from_table(table("c2", xs20, [-1.0 - 3.0 * i for i in range(20)], [18.4] * 20, [13] * 20, 13))
    check(pts is None and "4 segments beyond x = 0.8" in why,
          "another segmentation (20 segments, 4 beyond x = 0.8): no curve, a person decides")

    def summ(slope, fl, eb):
        return {"slope": slope, "contrast_m": None if fl is None else round(fl - eb, 4)}
    leg = {"c1": summ(0.13, 0.05, -0.05), "c2": summ(0.10, 0.04, -0.03)}
    v, o = rt.judge(leg, {"c1": summ(0.01, 0.004, -0.002), "c2": summ(-0.02, 0.01, 0.0)})
    check(o == "PASS" and v["c1"][0] == v["c2"][0] == "PASS", "both shrink in both cameras: PASS")
    v, o = rt.judge(leg, {"c1": summ(0.01, 0.004, -0.002), "c2": summ(-0.12, 0.01, 0.0)})
    check(o == "FAIL" and v["c2"][0] == "FAIL" and "does NOT shrink" in v["c2"][1],
          "a slope that does not shrink (here flipped and larger): FAIL")
    v, o = rt.judge(leg, {"c1": summ(0.01, 0.06, -0.06), "c2": summ(0.0, 0.0, 0.0)})
    check(o == "FAIL" and v["c1"][0] == "FAIL", "a flood/ebb contrast that does not shrink: FAIL")
    v, o = rt.judge(leg, {"c1": summ(0.01, None, None)})
    check(o == "NOT JUDGED" and v["c1"][0] == "not judged" == v["c2"][0], "no flood/ebb frames, or no camera: NOT JUDGED")
    v, o = rt.judge(leg, {"c1": summ(0.01, 0.004, -0.002), "c2": summ(0.0, None, None)})
    check(o == "PASS" and v["c2"][0] == "not judged", "one camera judged and passing, the other not judged: PASS")
    good = {"c1": summ(0.01, 0.004, -0.002), "c2": summ(-0.02, 0.01, 0.0)}
    v, o = rt.judge(leg, good, {"c1": 0.05, "c2": 0.0})
    check(o == "PASS" and "within 0.1 m" in v["c1"][1], "each camera's DEM moved about as predicted: PASS")
    v, o = rt.judge(leg, good, {"c1": 0.05, "c2": 0.15})
    check(o == "FAIL" and v["c2"][0] == "FAIL" and "suspect" in v["c2"][1],
          "slope and contrast shrink, but c2's DEM moves 0.15 m (a wrong mid bias file): FAIL")
    v, o = rt.judge(leg, good, {"c1": 0.05, "c2": None})
    check(o == "PASS" and v["c2"][0] == "not judged", "no DEM of c2 alone in both modes: c2 not judged")
    cfg = rt.cron_settings(HERE / "waterline_timex_cron.sh")
    check(cfg["SETUP_COEF"] == "0.037" and cfg["DEM_CELL"] == "2.0" and cfg["DEM_MAX_DAY_OFFSET"] == "0.15"
          and cfg["GNSSR_SPLINE"].endswith("usgs_spline_out.txt") and cfg["CONSISTENCY_ENABLE"] == "1",
          "the cron's settings read from waterline_timex_cron.sh")


def test_exposure(d):
    print("recalibrate_timing.py --steps exposure: the file times and the timesteps against the 300 s")
    import cv2
    import recalibrate_timing as rt
    v, _ = rt.exposure_verdict(np.array([640.0, 655.0, 700.0] * 10), np.array([600.0] * 5))
    check(v == "CONSISTENT", "written 640-700 s after the name, stacks 600 s: CONSISTENT")
    v, ln = rt.exposure_verdict(np.array([15.0, 20.0, 40.0] * 10), np.array([]))
    check(v == "MISMATCH" and "BEFORE a 600 s exposure" in ln[0],
          "written 15-40 s after the name (named by its end): MISMATCH")
    v, ln = rt.exposure_verdict(np.array([]), np.array([300.0] * 5))
    check(v == "MISMATCH" and "NOT the 600 s" in ln[1], "stacks of 300 s: MISMATCH")
    v, _ = rt.exposure_verdict(np.array([]), np.array([]))
    check(v == "NOT MEASURED", "nothing to read: NOT MEASURED")
    v, _ = rt.exposure_verdict(np.array([5000.0] * 10), np.array([]))
    check(v == "CONSISTENT", "written long after (the files waited): never a MISMATCH")

    def station(name, delay, rows, extra=()):
        base = Path(d) / "exposure" / name
        for sub in ("archive/images_timex", "archive/ras_c2"):
            (base / sub).mkdir(parents=True, exist_ok=True)
        for k in range(12):
            ep = T0 + 1800 * k
            for f in (base / "archive/images_timex" / f"{stamp(ep)}.c1.timex.jpg",
                      base / "archive/images_timex" / f"{stamp(ep)}.c2.timex.jpg"):
                f.write_bytes(b"x")
                os.utime(str(f), (ep + delay + k, ep + delay + k))
            ras = base / "archive/ras_c2" / f"{stamp(ep)}.c2.ras.tiff"
            cv2.imwrite(str(ras), np.zeros((rows, 8, 3), np.uint8))
        rc, log = run([HERE / "recalibrate_timing.py", "--base", base, "--nice", "0",
                       "--work", base / "work"] + list(extra), base)
        (base / "run.log").write_text(log)
        return rc, log

    rc, log = station("start_named", 645, 1200)
    check("EXPOSURE: CONSISTENT" in log and "timestacks last 600 s" in log and "DERIVE" in log,
          "a station whose timex are written ~645 s after their names, stacks 1200 rows: CONSISTENT, derive runs on")
    rc, log = station("end_named", 25, 1200)
    check(rc == 4 and "EXPOSURE: MISMATCH" in log and "DERIVE" not in log,
          "written ~25 s after their names: MISMATCH, exit 4, nothing else run")
    rc, log = station("end_named_alone", 25, 1200, ("--steps", "exposure"))
    check(rc == 4 and "EXPOSURE: MISMATCH" in log and "STOPPED" not in log,
          "the same with --steps exposure alone: exit 4 too")


def test_write_calibration(d):
    print("recalibrate_timing.write_mid_calibration: only the mid file, read back as the detector reads it")
    import recalibrate_timing as rt
    d = Path(d) / "calib"
    d.mkdir(parents=True, exist_ok=True)
    pts = {"c1": ((0.0, -8.2), (1.0, -8.2)), "c2": ((0.05, 4.0), (0.5, -6.0), (0.95, -15.0))}
    try:
        rt.write_mid_calibration(d / "timex_bias_legacy.json", pts, {})
        refused = False
    except ValueError as e:
        refused = "never changed" in str(e)
    check(refused and not (d / "timex_bias_legacy.json").exists(), "the legacy file name is refused")
    b1 = rt.write_mid_calibration(d / "timex_bias_mid.json", pts, {"what": "test"})
    b2 = rt.write_mid_calibration(d / "timex_bias_mid.json", pts, {"what": "test 2"})
    got = timex_timing.read_timex_bias(d / "timex_bias_mid.json", "mid")
    check(b1 is None and b2 is not None and b2.exists() and got == pts
          and json.loads((d / "timex_bias_mid.json").read_text())["about"]["what"] == "test 2",
          "written, read back exactly; a second write keeps the first as .bak")
    return d / "timex_bias_mid.json"


def gt_lines(cam):
    x = np.arange(W) / (W - 1)
    return 120.0 + 30 * x if cam == "c1" else 40.0 + 140.0 * x


def test_derive(d):
    print("recalibrate_timing.py --steps derive: synthetic detections -> calibration/timex_bias_mid.json")
    d = Path(d) / "derive"
    gt, photos, work, cal = d / "gt", d / "photos", d / "work", d / "cal"
    for p in (gt, photos, work / "derive" / "blend_on" / "detections", work / "derive" / "blend_off" / "detections", cal):
        p.mkdir(parents=True, exist_ok=True)
    rng = np.random.RandomState(5)
    for k in range(7):
        for cam in ("c1", "c2"):
            stem = f"{stamp(T0 + 1800 * k)}.{cam}.timex"
            truth = gt_lines(cam)
            (gt / f"{stem}.csv").write_text("Column,Row\n" + "".join(f"{c},{int(round(r))}\n" for c, r in enumerate(truth)))
            (photos / f"{stem}.jpg").write_bytes(b"")              # found, never read (--reuse-detections)
            x = np.arange(W) / (W - 1)
            err_off = (-8.0 + rng.normal(0, 1.0)) if cam == "c1" else (4.0 - 22.0 * x + rng.normal(0, 1.0))
            for key, shift, mode in (("blend_off", 0.0, "mid"), ("blend_on", -2.8 if cam == "c1" else -2.6, "")):
                rows = np.round(truth) + err_off + shift + rng.normal(0, 0.3, W)
                lines = ["Column,Row,Confidence,Row_Precise,Peak_Sharpness,Has_Signal" + (",Timing_Mode" if mode else "")]
                lines += [f"{c},{int(round(r))},0.95,{r:.3f},0.5,1" + (f",{mode}" if mode else "") for c, r in enumerate(rows)]
                (work / "derive" / key / "detections" / f"{stem}.csv").write_text("\n".join(lines) + "\n")
    def repo_calibration():         # the repository's timex bias files (a station may have its mid one)
        return {n: sha(HERE / "calibration" / n) if (HERE / "calibration" / n).exists() else None
                for n in ("timex_bias_legacy.json", "timex_bias_mid.json")}
    repo_before = repo_calibration()
    args = [HERE / "recalibrate_timing.py", "--steps", "derive", "--reuse-detections", "--work", work,
            "--calibration-dir", cal, "--ground-truth-dir", gt, "--image-dirs", photos, "--nice", "0",
            "--base", d]
    rc, log = run(args, d)
    (d / "derive.log").write_text(log)
    f = cal / "timex_bias_mid.json"
    data = json.loads(f.read_text()) if f.exists() else {}
    c1 = data.get("c1") or [[0, 0]]
    c2 = data.get("c2") or []
    check(rc == 0 and data.get("mode") == "mid" and len(c1) == 2 and abs(c1[0][1] + 8.0) < 1.0
          and "CONSTANT" in data["about"]["c1"],
          f"c1: a constant near the planted -8 px ({c1[0][1]})")
    check(len(c2) == 10 and c2[0][1] > c2[-1][1] and "PER-COLUMN" in data["about"]["c2"]
          and abs(c2[0][1] - (4.0 - 22.0 * c2[0][0])) < 1.5,
          f"c2: the planted per-column shape ({c2[0] if c2 else None} .. {c2[-1] if c2 else None})")
    shifts = data.get("about", {}).get("control_blend_on", {})
    check(abs(shifts.get("c1", {}).get("blend_shift_px_mean", 0) - 2.8) < 0.3
          and abs(shifts.get("c2", {}).get("blend_shift_px_mean", 0) - 2.6) < 0.3
          and "blend lag the v6.8 values absorbed" in log,
          "the blend lag planted in the control (2.8 / 2.6 px) is what mid minus control measures")
    check(repo_calibration() == repo_before,
          "the repository's calibration/ is untouched (the test wrote only into its own folder)")
    try:
        readable = f.exists() and timex_timing.read_timex_bias(f, "mid") is not None
    except ValueError:
        readable = False
    check(readable, "the file written reads as the detector reads a mid timex bias file")
    rc2, log2 = run(args[:-2] + ["--min-frames", "8", "--base", d], d)
    check(rc2 == 1 and "only 7 ground-truth photo(s) found" in log2, "too few ground-truth photos: stopped, nothing written")
    return f


def test_compare(d, mid_file):
    print("recalibrate_timing.py --steps compare: the cron's chain both ways, residuals vs dz/dt, PASS")
    d = Path(d) / "compare"
    work = d / "work"
    start, end = "2026-10-01", "2026-10-04"
    cdir = work / "compare" / f"{start}_{end}"
    t0 = int(datetime(2026, 10, 1, 0, 0, tzinfo=timezone.utc).timestamp())
    w = 2 * np.pi / (12.42 * 3600)

    def z(t):
        return 0.4 + 1.1 * np.sin(w * (np.asarray(t, float) - t0))
    # GNSS-R spline every 30 min over the window
    ts = np.arange(t0 - 3600, t0 + 5 * 86400, 1800)
    d.mkdir(parents=True, exist_ok=True)
    with open(d / "spline.txt", "w") as f:
        f.write("% orthometric height minus RH, where Hortho (m) is    19.014 \n")
        for t in ts:
            u = datetime.fromtimestamp(int(t), tz=timezone.utc)
            zz = float(np.round(z(t), 3))
            f.write(f"{t / 86400.0 + 40587.0:15.7f} {19.014 - zz:10.3f} {u.year:4d} {u.month:3d} {u.day:3d} "
                    f"{u.hour:3d} {u.minute:3d} {u.second:3d} {zz:10.3f}\n")
    # waterlines: a real 11 Sep 2026 hand-clicked line of each camera moved by k px per m of the level
    # it shows (the review's reproduction: c1 41.0, c2 24.8 px/m); legacy lines show the level 472 s
    # BEFORE the file-name epoch (7.9 min older than the level legacy gives them, the steady state of
    # the blend), mid lines the level at mid-exposure (what mid gives them)
    K = {"c1": 41.0, "c2": 24.8}
    tmpl = {}
    for cam in ("c1", "c2"):
        g = HERE / "ground_truth" / f"1789131600.Fri.Sep.11_13_00_00.GMT.2026.CACO05.{cam}.timex.csv"
        rows = np.array([float(r["Row"]) for r in csv.DictReader(open(g))])
        tmpl[cam] = rows + 0.5 * np.sin(2 * np.pi * np.arange(len(rows)) / 40.0)
    cols = np.arange(0, 2448, 4)
    rng = np.random.RandomState(11)
    for mode, shown in (("legacy", -472), ("mid", 300)):
        det = cdir / mode / "detections"
        det.mkdir(parents=True, exist_ok=True)
        for day in range(4):
            for k in range(22):                                  # 11:00 - 21:30 UTC
                ep = t0 + day * 86400 + 11 * 3600 + 1800 * k
                lvl = float(z(ep + shown))
                for cam in ("c1", "c2"):
                    r = tmpl[cam][cols] + K[cam] * (lvl - 0.4) + rng.normal(0, 0.3)
                    m = "" if mode == "legacy" else ",mid"
                    lines = ["Column,Row,Confidence,Row_Precise,Peak_Sharpness,Has_Signal" + (",Timing_Mode" if m else "")]
                    lines += [f"{c},{int(round(v))},0.95,{v:.3f},0.5,1{m}" for c, v in zip(cols, r)]
                    (det / f"{stamp(ep)}.{cam}.timex.csv").write_text("\n".join(lines) + "\n")
    rc, log = run([HERE / "recalibrate_timing.py", "--steps", "compare", "--reuse-detections", "--work", work,
                   "--start-date", start, "--end-date", end, "--gnssr", d / "spline.txt",
                   "--waves", d / "no_waves.csv", "--gauge", d / "no_gauge.csv",
                   "--calibration-dir", HERE / "calibration", "--nice", "0", "--base", d], d)
    (d / "compare.log").write_text(log)
    rep = json.loads((cdir / "report.json").read_text()) if (cdir / "report.json").exists() else {}
    res = rep.get("residuals", {})
    lag = {m: {c: (res.get(m, {}).get(c) or {}).get("lag_min") for c in ("c1", "c2")} for m in ("legacy", "mid")}
    con = {m: {c: (res.get(m, {}).get(c) or {}).get("contrast_m") for c in ("c1", "c2")} for m in ("legacy", "mid")}
    check(rc == 0 and rep.get("overall") == "PASS" and "VERDICT: PASS" in log,
          f"PASS (exit {rc}); lag legacy {lag['legacy']} min, mid {lag['mid']} min; flood-ebb legacy "
          f"{con['legacy']} m, mid {con['mid']} m")
    ok = all(v is not None for m in lag.values() for v in m.values())
    check(ok and all(4.0 < lag["legacy"][c] < 10.0 and abs(lag["mid"][c]) < 1.5 for c in ("c1", "c2"))
          and all(con["legacy"][c] > 0.03 for c in ("c1", "c2")),
          "the planted 7.9 min comes back in legacy (flood above ebb); mid near 0")
    g = rep.get("grids") or {}
    check(g.get("cells_both", 0) > 50 and g.get("mid_minus_legacy_median_m") is not None
          and "DEM: mid minus legacy over" in log, f"the two DEMs compared over {g.get('cells_both')} shared cells")
    gc = rep.get("grids_by_camera") or {}
    check(all((gc.get(c) or {}).get("cells_both", 0) > 20 for c in ("c1", "c2")) and "DEM of c2 alone" in log
          and all("its DEM mid minus legacy" in rep["verdicts"][c][1] for c in ("c1", "c2")),
          "each camera's DEM built alone in both modes and its mid minus legacy judged: "
          + ", ".join(f"{c} {(gc.get(c) or {}).get('mid_minus_legacy_median_m')} m" for c in ("c1", "c2")))
    ex = (cdir / "mid" / "contours.csv").read_text().splitlines()[0] if (cdir / "mid" / "contours.csv").exists() else ""
    check(ex.endswith(",level_time_offset_s") and "--timing-mode mid" in (cdir / "compare.log").read_text(),
          "the mid chain ran its contour step with --timing-mode mid")


def test_py38():
    print("Python 3.8 syntax")
    import ast
    for f in ("timex_timing.py", "recalibrate_timing.py", "test_timing_mode.py", "waterline_detector_v5.py",
              "extract_elevation_contours.py", "dem_from_contours.py", "derive_bias_correction.py",
              "process_historical.py", "process_chelsea.py", "station_status.py", "runup_from_timestack.py",
              "timestack_wave_period.py", "relabel_manifests.py", "reaim_figure.py", "survey_products.py"):
        try:
            ast.parse((HERE / f).read_text(encoding="utf-8"), feature_version=(3, 8))
            check(True, f"{f} parses as Python 3.8")
        except TypeError:                        # Python 3.7 has no feature_version
            ast.parse((HERE / f).read_text(encoding="utf-8"))
        except SyntaxError as e:
            check(False, f"{f}: {e}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="keep the outputs, in a new folder inside this one")
    ap.add_argument("--before", default=None,
                    help="the commit to compare legacy with (default: the one before the timing switch)")
    args = ap.parse_args()
    if args.keep:
        Path(args.keep).mkdir(parents=True, exist_ok=True)
    d = os.path.realpath(tempfile.mkdtemp(prefix="timing_mode_", dir=args.keep))
    print(f"test folder: {d}")
    ref, entry, note = before_ref(args.before)
    print("the code without the switch: " + ("not readable (git)" if ref is None else
          f"this checkout with {entry[:7]}'s changes taken out (else {ref})" if entry else ref)
          + (f" (NOTE: {note})" if note else ""))
    try:
        test_definitions()
        test_legacy_bias_file()
        test_shared_constant()
        mid_file = test_write_calibration(d)
        test_detector(d, ref, mid_file, entry, note)
        test_extract(d, ref, entry, note)
        test_dem(d, ref, entry, note)
        test_cron()
        test_station_status(d)
        test_rules()
        test_exposure(d)
        test_derive(d)
        test_compare(d, mid_file)
        test_py38()
    finally:
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED")
        return 1
    print("all timing-mode tests passed" + (f" ({len(SKIPPED)} skipped: {'; '.join(SKIPPED)})" if SKIPPED else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())

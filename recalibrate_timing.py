#!/usr/bin/env python3
"""
Recalibrate For The Water-Level Timing Fix, And Test It (on the NUC)
=====================================================================
Everything the station needs before waterline_timex_cron.sh may be set to
TIMING_MODE=mid (timex_timing.py), and the test that decides whether it
should be. It changes no product of the station: it writes only into its
own work folder, calibration/timex_bias_mid.json (derive) and, when asked,
archive/processed_timex_mid (backfill), which the cron reads only in mode
mid. calibration/timex_bias_legacy.json is never written.

WHY ON THE NUC. The fix is three changes that only work together (known
finding C-B1; timex_timing.py): the level at mid-exposure, the detector's
temporal blend off, and timex bias corrections re-derived with the blend
off -- the v6.8 ones were fitted with it on, on 11 Sep 2026 frames that
were all on the flood, and absorbed ~2.8 px (c1) / ~2.6 px (c2) of its lag.
Re-deriving needs the photos of the hand-clicked ground truth, and the
test needs a window of photos and the GNSS-R record: they are on the NUC.

STEPS (--steps; default exposure,derive,compare):

  exposure  is the collection what TIMEX_MID_OFFSET_S assumes (10 min,
            named by its start)? Nothing on the station's side has said so
            (the arguseyes collection settings are not in this repository),
            so two things the station's own files show are read:
              * when each archived timex was written: its modification
                time (the cron copies with cp -p, which keeps it) minus the
                time in its name. A collection named by its START is
                written only after it ends, so at least 2 x
                TIMEX_MID_OFFSET_S (600 s) after the name; one written
                sooner was not a 600 s exposure starting at its name;
              * how long the c2 timestacks (archive/ras_c2) last: their
                rows at 2 Hz (runup_from_timestack.py), 1200 = 600 s.
            CONSISTENT, MISMATCH or NOT MEASURED (no file kept its time).
            A MISMATCH stops the run before the other steps (--ignore-
            exposure goes on): the constant is then the one number to
            change, in timex_timing.py. The written-after test cannot tell
            a 10-min exposure named by its start from a shorter one, and
            the timestacks are the c2 timestack's own collection, assumed
            to match the timex's; the docs' section 2.7 has the check by
            hand in the arguseyes settings.

  derive    1. the photos of the hand-clicked timex ground truth
               (<ground truth>/*.timex.csv) are looked up in --image-dirs;
            2. they are detected twice with the bias correction off: with
               the blend ON (as when the v6.8 values were fitted: the
               CONTROL) and with it OFF (--timing-mode mid);
            3. derive_bias_correction.py --suffix .timex. --uncorrected
               measures each camera's raw error on both;
            4. the blend-off table becomes the mid corrections by the rules
               the v6.8 values were set by (points_from_table), written to
               calibration/timex_bias_mid.json with its provenance (a
               previous one is kept beside it as .bak);
            5. legacy, control and mid are printed side by side: control
               minus legacy says whether the fit itself is reproducible,
               mid minus control is the blend lag the v6.8 values absorbed
               (the review's reproduction: about +2.8 px c1, +2.6 px c2).

  compare   the timex photos of a window (--start-date, --end-date, UTC
            days of archive/images_timex) are detected in BOTH modes and
            each run through the cron's own chain, with the cron's settings
            (read from waterline_timex_cron.sh): contours (GNSS-R level,
            waves, setup, GNSS-R QC), the consistency filter, georectify,
            and the DEM with its diagnostics. Then, per camera, the DEM's
            residuals against the tide's rate (dem_from_contours.py
            residual_diagnostic: each frame's median offset from the other
            frames' cell medians, against dz/dt) for legacy and mid, and
            the verdict of judge():

              PASS for a camera when, from legacy to mid, BOTH shrink:
                * |slope| of the residual on dz/dt (m per m/h; x 60 = the
                  minutes a line is older than its level), and
                * |flood - ebb|, the median residual of the flood frames
                  minus that of the ebb frames (|dz/dt| > 0.1 m/h).
              Overall PASS when every camera that could be judged passes.

            The timing error predicts legacy about +4 to +8 min (lines
            older than their levels: flood high, ebb low) and mid about 0.
            The minutes are indicative (a planted 6 min comes back as
            5.5-8 min); the sign, the flood/ebb contrast and which of the
            two is smaller are the test. Also printed: mid minus legacy
            over the DEM cells both have, and each DEM's median spread.
            --fit-setup also refits the setup coefficient C on contours
            built with C = 0 in both modes (dem_from_contours.py
            --fit-setup): C was fitted on legacy lines and may have
            absorbed a little of their timing error (~0.004-0.008).

  backfill  re-detects every photo of archive/images_timex (or of the
            --start-date..--end-date days) in mode mid into
            archive/processed_timex_mid, one UTC day at a time, never
            overwriting a file there; days done are listed in the work
            folder, so a stopped run resumes. Run it once, after derive
            and a PASS, before setting TIMING_MODE=mid in the cron.

Exit status: 0 done (compare: PASS), 2 compare FAIL, 3 compare could not
judge any camera, 4 the exposure check found a MISMATCH (nothing else
run), 1 an error.

Usage (on the NUC, in /mnt/I2Rgus_Data/waterline):
    python3 recalibrate_timing.py --start-date 2026-09-29 --end-date 2026-10-05
    python3 recalibrate_timing.py --steps exposure
    python3 recalibrate_timing.py --steps derive
    python3 recalibrate_timing.py --steps compare --start-date ... --end-date ... --fit-setup
    python3 recalibrate_timing.py --steps backfill
docs/VALIDATION.md section 2.7 is the full procedure.
"""

import os
import re
import sys
import csv
import json
import shutil
import socket
import argparse
import subprocess
from pathlib import Path
from datetime import datetime, timezone

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import timex_timing                                   # noqa: E402

BASE = Path("/mnt/I2Rgus_Data/waterline")
DETECTOR = HERE / "waterline_detector_v5.py"
CAMERAS = ("c1", "c2")

# The v6.5/v6.8 rule for the shape of a correction: a per-column curve only
# when derive_bias_correction.py's raw error varies across the frame by more
# than this many standard errors; otherwise a constant (fitting the wiggle of
# one day once put a 35 px tilt into c1's snap correction that was not in the
# data).
SPREAD_SEM_CURVE = 3.0

# The review's reproduction of what the v6.8 timex corrections absorbed from
# the blend (mid minus control, px), printed beside the measured values.
EXPECTED_BLEND_SHIFT_PX = {"c1": 2.8, "c2": 2.6}

# The cron's settings this script reads from waterline_timex_cron.sh, and
# what it assumes when a line is missing.
CRON_DEFAULTS = {
    "GNSSR_SPLINE": "/home/argus_user/GNSS/v4.1/products/refl_code/Files/usgs/usgs_spline_out.txt",
    "SETUP_COEF": "0.037", "DEM_CELL": "2.0", "DEM_MIN_POINTS": "3", "DEM_MAX_SPREAD": "0.5",
    "DEM_MAX_HS": "1.5", "DEM_MAX_DAY_OFFSET": "0.15", "CONSISTENCY_ENABLE": "1",
    "USE_MARCONI_WAVES": "1", "WAVE_BUOY": "44008", "TIDE_GAUGE": "8447435",
}


def say(text=""):
    print(text, flush=True)


# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------

def cron_settings(path):
    """
    The cron's own values of CRON_DEFAULTS' names (NAME=value lines,
    quotes stripped; $BASE expanded), so the window is processed with the
    settings the station uses. A missing file or line keeps the default.
    """
    out = dict(CRON_DEFAULTS)
    base = str(BASE)
    try:
        text = Path(path).read_text()
    except OSError:
        return out
    m = re.search(r"^BASE=(\S+)\s*$", text, re.M)
    if m:
        base = m.group(1).strip("\"'")
    for name in CRON_DEFAULTS:
        m = re.search(rf"^{name}=(.*?)\s*(?:#.*)?$", text, re.M)
        if m:
            out[name] = m.group(1).strip().strip("\"'").replace("$BASE", base).replace("${BASE}", base)
    return out


def epoch_of(name):
    m = re.match(r"^(\d{9,11})\.", Path(name).name)
    return int(m.group(1)) if m else None


def camera_of(name):
    n = Path(name).name.lower()
    return "c1" if ".c1." in n else "c2" if ".c2." in n else None


def utc_day(epoch):
    return datetime.fromtimestamp(epoch, tz=timezone.utc).strftime("%Y-%m-%d")


def place(src, dst):
    """Hard link, else symbolic link, else copy: the detector only reads it."""
    try:
        os.link(str(src), str(dst))
        return
    except OSError:
        pass
    try:
        os.symlink(str(Path(src).resolve()), str(dst))
    except OSError:
        shutil.copy2(str(src), str(dst))


# ---------------------------------------------------------------------
# Running the station's scripts
# ---------------------------------------------------------------------

def run(cmd, log, progress=None):
    """Runs cmd with its output appended to log. progress = number of photos: every
    tenth 'Processing:' line of the detector is printed, so a long run is never silent.
    -> exit code."""
    cmd = [str(c) for c in cmd]
    with open(log, "a") as f:
        f.write("\n$ " + " ".join(cmd) + "\n")
        f.flush()
        if not progress:
            return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT).returncode
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True, bufsize=1)
        k, every = 0, max(1, progress // 10)
        for line in p.stdout:
            f.write(line)
            if line.startswith("Processing:"):
                k += 1
                if k % every == 0 or k == progress:
                    say(f"      photo {k} of {progress}")
        return p.wait()


def detector_cmd(mode, src, work, timex_bias_file=None, no_bias=False):
    """The detector as the cron runs it, in timing mode `mode`, on the photos in src;
    its folders under work. no_bias: bias correction off (the derive step)."""
    cmd = [sys.executable, DETECTOR, "--image-suffix", "timex.jpg", "--timing-mode", mode,
           "--source-dir", src, "--input-dir", work / "in", "--output-dir", work / "detections",
           "--debug-dir", work / "debug"]
    if no_bias:
        cmd.append("--no-bias-correction")
    elif timex_bias_file is not None:
        cmd += ["--timex-bias-file", timex_bias_file]
    return cmd


def detect(mode, photos, work, log, timex_bias_file=None, no_bias=False, keep_debug=False):
    """Detects `photos` in timing mode `mode` into work/detections. -> list of detection CSVs."""
    src = work / "src"
    for d in (src, work / "in", work / "detections", work / "debug"):
        shutil.rmtree(str(d), ignore_errors=True)
        d.mkdir(parents=True)
    for p in photos:
        place(p, src / Path(p).name)
    rc = run(detector_cmd(mode, src, work, timex_bias_file, no_bias), log, progress=len(photos))
    if not keep_debug:
        shutil.rmtree(str(work / "debug"), ignore_errors=True)
    shutil.rmtree(str(work / "in"), ignore_errors=True)            # the detector's own copies
    if rc != 0:
        raise RuntimeError(f"the detector failed (exit {rc}); see {log}")
    return sorted((work / "detections").glob("*.csv"))


# ---------------------------------------------------------------------
# exposure
# ---------------------------------------------------------------------

# How many of the newest archived files the exposure check reads.
EXPOSURE_SAMPLE = 200
# The timestacks' sampling rate (runup_from_timestack.SAMPLE_HZ).
STACK_HZ = 2.0
# A timestack whose length differs from the assumed exposure by more than
# this share is a MISMATCH.
STACK_TOLERANCE = 0.05
# Slack (s) on "written only after the exposure ended": file times are
# whole seconds and the collection's clock and the file system's may
# differ by a little.
WRITE_SLACK_S = 30.0


def newest(folder, pattern, n=EXPOSURE_SAMPLE):
    """The n files of folder matching pattern with the latest name epochs."""
    files = [p for p in Path(folder).glob(pattern) if epoch_of(p.name) is not None]
    return sorted(files, key=lambda p: epoch_of(p.name))[-n:]


def write_delays(paths):
    """Seconds from each file's name time to its modification time: when it was written (the
    archive copies keep it). A delay below 0 or over a day is left out: that copy did not keep
    its time."""
    out = []
    for p in paths:
        try:
            d = Path(p).stat().st_mtime - epoch_of(Path(p).name)
        except (OSError, TypeError):
            continue
        if 0.0 <= d <= 86400.0:
            out.append(d)
    return np.array(out, float)


def stack_seconds(paths):
    """How long each timestack lasts: its rows (time) at STACK_HZ."""
    import cv2
    out = []
    for p in paths:
        img = cv2.imread(str(p), cv2.IMREAD_UNCHANGED)
        if img is not None and img.ndim >= 2:
            out.append(img.shape[0] / STACK_HZ)
    return np.array(out, float)


def exposure_verdict(delays, stacks, mid=timex_timing.TIMEX_MID_OFFSET_S):
    """
    Whether the files agree with an exposure of 2 x mid seconds named by
    its start. -> (CONSISTENT / MISMATCH / NOT MEASURED, [lines]).
      * timex write delays (write_delays): MISMATCH when the 10th
        percentile is below 2 x mid - WRITE_SLACK_S, i.e. most files were
        written before such an exposure could have ended (the name is then
        not the start, or the exposure is shorter);
      * timestack lengths (stack_seconds): MISMATCH when the median differs
        from 2 x mid by more than STACK_TOLERANCE.
    Either measured and neither a MISMATCH: CONSISTENT. A long delay says
    nothing (the files may have waited to be written), so it is never a
    MISMATCH.
    """
    t_exp = 2.0 * mid
    lines, bad, measured = [], False, False
    if len(delays):
        measured = True
        p10, med, p90 = (float(np.percentile(delays, q)) for q in (10, 50, 90))
        early = p10 < t_exp - WRITE_SLACK_S
        bad |= early
        lines.append(f"timex written {med:.0f} s after the time in its name (median of {len(delays)}; "
                     f"10-90%: {p10:.0f}-{p90:.0f} s): "
                     + (f"BEFORE a {t_exp:.0f} s exposure starting at the name could end -- the name is not the "
                        f"start, or the exposure is shorter" if early else
                        f"after a {t_exp:.0f} s exposure starting at the name ends, as it should be"))
    else:
        lines.append("timex write time: no archived timex kept its modification time (or none found)")
    if len(stacks):
        measured = True
        med = float(np.median(stacks))
        off = abs(med - t_exp) > STACK_TOLERANCE * t_exp
        bad |= off
        lines.append(f"timestacks last {med:.0f} s (median of {len(stacks)}, rows at {STACK_HZ:g} Hz): "
                     + (f"NOT the {t_exp:.0f} s the constant assumes" if off else f"the {t_exp:.0f} s the constant assumes"))
    else:
        lines.append("timestack length: no timestack found")
    verdict = "MISMATCH" if bad else "CONSISTENT" if measured else "NOT MEASURED"
    return verdict, lines


def step_exposure(args):
    say("=" * 78)
    say(f"EXPOSURE: is a timex {2 * timex_timing.TIMEX_MID_OFFSET_S} s long and named by its start "
        f"(TIMEX_MID_OFFSET_S = {timex_timing.TIMEX_MID_OFFSET_S} s)?")
    say("=" * 78)
    delays = write_delays(newest(args.images, "*.timex.jpg"))
    stacks = stack_seconds(newest(args.stacks, "*.ras.tiff", 20))
    verdict, lines = exposure_verdict(delays, stacks)
    for ln in lines:
        say("   " + ln)
    say(f"EXPOSURE: {verdict}" + {
        "CONSISTENT": " with the middle of the exposure at the file-name epoch + "
                      f"{timex_timing.TIMEX_MID_OFFSET_S} s.",
        "MISMATCH": " -- do not go on with this constant: find the exposure in the arguseyes collection "
                    "settings and change TIMEX_MID_OFFSET_S in timex_timing.py (docs/VALIDATION.md section 2.7).",
        "NOT MEASURED": " -- check the exposure by hand in the arguseyes collection settings "
                        "(docs/VALIDATION.md section 2.7)."}[verdict])
    return verdict


# ---------------------------------------------------------------------
# derive
# ---------------------------------------------------------------------

def ground_truth_photos(gt_dir, image_dirs):
    """{camera: [(ground-truth csv, photo or None)]} for every *.timex.csv in gt_dir."""
    out = {c: [] for c in CAMERAS}
    for g in sorted(Path(gt_dir).glob("*.timex.csv")):
        cam = camera_of(g.name)
        if cam not in out:
            continue
        photo = None
        for d in image_dirs:
            p = Path(d) / (g.name[:-len(".csv")] + ".jpg")
            if p.exists():
                photo = p
                break
        out[cam].append((g, photo))
    return out


def points_from_table(table):
    """
    derive_bias_correction.py's table (--json) -> (points, decision), by the
    rules the v6.8 timex corrections were set by:
      * SHAPE: the raw error's spread across the frame (max - min of the
        segment means) against its standard error (the median over the
        segments of std / sqrt(frames)). Spread/SEM above SPREAD_SEM_CURVE
        (3) is a real shape: one point per segment. Otherwise the shape is
        noise and the correction is the mean, constant across the frame
        (c1 in v6.8: 11.2 px against 4.2 px, 2.7, constant -11.03).
      * FAR-FIELD CLAMP (curves only): from the far end inward, segments
        whose error bar is larger than their value (std > |mean|), or that
        too few frames reach (fewer than max(3, half the frames), derive's
        own FEW FRAMES flag), take the value of the nearest segment that is
        neither (c2 in v6.8: -41.28 +/- 48.26 and -57.78 +/- 67.00 became
        -17.22). The near end is kept as measured, as it was.
    Segments no frame reaches are left out (the detector interpolates, and
    holds the end values beyond the last point). Values are rounded as
    derive_bias_correction.py prints them (x to 0.001, px to 0.01).
    """
    segs = [s for s in table["segments"] if s.get("n", 0) > 0 and s.get("raw") is not None]
    if len(segs) < 2:
        raise ValueError(f"{table.get('camera')}: fewer than 2 segments with data")
    x = np.array([s["x"] for s in segs], float)
    raw = np.array([s["raw"] for s in segs], float)
    std = np.array([s["std"] if s.get("std") is not None else 0.0 for s in segs], float)
    n = np.array([s["n"] for s in segs], float)
    sem = float(np.median(std / np.sqrt(n)))
    spread = float(raw.max() - raw.min())
    ratio = spread / sem if sem > 0 else float("inf")
    if ratio <= SPREAD_SEM_CURVE:
        c = round(float(np.mean(raw)), 2)
        return ((0.0, c), (1.0, c)), (f"CONSTANT {c:+.2f} px: spread across the frame {spread:.1f} px against a "
                                      f"{sem:.1f} px standard error (spread/SEM {ratio:.1f}, not above "
                                      f"{SPREAD_SEM_CURVE:g}: the shape is noise); the mean of {len(segs)} segments")
    few = max(3.0, 0.5 * float(table.get("n_frames", 0)))
    untrusted = (std > np.abs(raw)) | (n < few)
    k = len(segs) - 1
    while k > 0 and untrusted[k]:
        k -= 1
    vals = raw.copy()
    vals[k + 1:] = raw[k]
    pts = tuple((round(float(a), 3), round(float(b), 2)) for a, b in zip(x, vals))
    clamped = [f"{x[i]:.3f} ({raw[i]:+.2f} +/- {std[i]:.2f}, {int(n[i])} frames)" for i in range(k + 1, len(segs))]
    return pts, (f"PER-COLUMN, {len(pts)} points: spread {spread:.1f} px against a {sem:.1f} px standard error "
                 f"(spread/SEM {ratio:.1f}, above {SPREAD_SEM_CURVE:g})"
                 + (f"; far-field segment(s) {', '.join(clamped)} CLAMPED to {raw[k]:+.2f} (x {x[k]:.3f}): "
                    f"error bar larger than the value, or too few frames" if clamped else "; nothing clamped"))


def correction_at(points, xs):
    """A correction's value at x fractions (np.interp, as apply_bias_correction)."""
    return np.interp(np.asarray(xs, float), [p[0] for p in points], [p[1] for p in points])


def write_mid_calibration(path, points, about):
    """Writes the mid timex bias file (and keeps a previous one as .bak). Refuses any
    other name: the legacy file holds the v6.8 values and is never written here."""
    path = Path(path)
    if path.name != timex_timing.timex_bias_path("mid").name:
        raise ValueError(f"refusing to write {path}: only {timex_timing.timex_bias_path('mid').name} is "
                         f"written here (the legacy file is never changed)")
    path.parent.mkdir(parents=True, exist_ok=True)
    backup = None
    if path.exists():
        backup = path.with_name(path.name + "." + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + ".bak")
        shutil.copy2(str(path), str(backup))
    data = {"mode": "mid", "c1": [list(p) for p in points["c1"]], "c2": [list(p) for p in points["c2"]],
            "about": about}
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(data, indent=1) + "\n")
    os.replace(str(tmp), str(path))
    timex_timing.read_timex_bias(path, "mid")          # what the detector will read
    return backup


def derive_table(gt_dir, det_dir, cam, out_json, log):
    """derive_bias_correction.py --suffix .timex. --uncorrected on det_dir -> its table."""
    rc = run([sys.executable, HERE / "derive_bias_correction.py", "--ground-truth-dir", gt_dir,
              "--processed-dir", det_dir, "--suffix", ".timex.", "--uncorrected", "--json", out_json, cam], log)
    if rc != 0 or not Path(out_json).exists():
        raise RuntimeError(f"derive_bias_correction.py failed for {cam} on {det_dir} (exit {rc}); see {log}")
    return json.loads(Path(out_json).read_text())


def step_derive(args, work, cal_dir):
    say("=" * 78)
    say("DERIVE: the mid timex bias corrections, from the ground truth detected with the blend off")
    say("=" * 78)
    gt_dir = Path(args.ground_truth_dir)
    gt = ground_truth_photos(gt_dir, args.image_dirs)
    wdir = work / "derive"
    wdir.mkdir(parents=True, exist_ok=True)
    log = wdir / "derive.log"
    photos = []
    for cam in CAMERAS:
        found = [p for _, p in gt[cam] if p is not None]
        missing = [g.name for g, p in gt[cam] if p is None]
        say(f"{cam}: {len(gt[cam])} timex ground-truth line(s) in {gt_dir}, photo found for {len(found)}"
            + (f"; not found: {', '.join(missing[:4])}" + (" ..." if len(missing) > 4 else "") if missing else ""))
        if len(found) < args.min_frames:
            raise RuntimeError(f"{cam}: only {len(found)} ground-truth photo(s) found (--min-frames "
                               f"{args.min_frames}) in {', '.join(map(str, args.image_dirs))}: give the folder "
                               f"that holds them with --image-dirs, or click more (collect_all_ground_truth.py)")
        photos += found
    runs = {"blend_on": ("legacy", "CONTROL: blend ON, bias off (as when the v6.8 values were fitted)"),
            "blend_off": ("mid", "MID: blend OFF, bias off")}
    tables = {}
    for key, (mode, what) in runs.items():
        rdir = wdir / key
        if args.reuse_detections and list((rdir / "detections").glob("*.csv")):
            say(f"{what}: reusing {rdir / 'detections'}")
        else:
            say(f"{what}: detecting {len(photos)} photo(s) (log {log})")
            detect(mode, photos, rdir, log, no_bias=True, keep_debug=args.keep_debug)
        for cam in CAMERAS:
            tables[(key, cam)] = derive_table(gt_dir, rdir / "detections", cam, wdir / f"derive_{cam}_{key}.json", log)

    import waterline_detector_v5 as detector
    points, decisions, summary = {}, {}, {}
    say()
    for cam in CAMERAS:
        legacy = detector.CAMERAS[f"CACO05_{cam.upper()}"].bias_correction_points_timex
        t_on, t_off = tables[("blend_on", cam)], tables[("blend_off", cam)]
        pts, why = points_from_table(t_off)
        ctl_pts, ctl_why = points_from_table(t_on)
        points[cam], decisions[cam] = pts, why
        xs = [s["x"] for s in t_off["segments"]]
        leg_at, ctl_at, mid_at = correction_at(legacy, xs), correction_at(ctl_pts, xs), correction_at(pts, xs)

        def raw_of(t):        # a segment no frame reaches has no raw error (derive writes 0 there)
            return np.array([s["raw"] if s["raw"] is not None and s["n"] > 0 else np.nan
                             for s in t["segments"]], float)
        raw_on, raw_off = raw_of(t_on), raw_of(t_off)
        both = np.isfinite(raw_on) & np.isfinite(raw_off)
        shift = float(np.mean((raw_off - raw_on)[both])) if both.any() else float("nan")
        repro = float(np.mean(ctl_at - leg_at))
        say(f"{cam}: {t_off['n_frames']} frame(s) blend off, {t_on['n_frames']} blend on")
        say(f"   {'x':>6} | {'legacy':>8} | {'control':>8} | {'mid raw':>8} | {'mid':>8} | {'std':>6} | {'n':>3}")
        for i, s in enumerate(t_off["segments"]):
            say(f"   {s['x']:6.3f} | {leg_at[i]:+8.2f} | {raw_on[i]:+8.2f} | {raw_off[i]:+8.2f} | {mid_at[i]:+8.2f} | "
                f"{(s['std'] if s['std'] is not None else float('nan')):6.2f} | {s['n']:3d}")
        say(f"   mid: {why}")
        say(f"   blend lag the v6.8 values absorbed (mid minus control, raw error): {shift:+.2f} px mean "
            f"(the review's reproduction: about {EXPECTED_BLEND_SHIFT_PX[cam]:+.1f} px)")
        say(f"   control by the same rules minus the v6.8 values: {repro:+.2f} px mean"
            + ("  <-- NOTE: the fit with the blend on does not reproduce the v6.8 values to 3 px: the detector, "
               "the ground truth or its photos changed since; mid is still a like-for-like fit"
               if abs(repro) > 3.0 else " (the fit is reproducible)"))
        summary[cam] = {"control_points": [list(p) for p in ctl_pts], "control_decision": ctl_why,
                        "blend_shift_px_mean": round(shift, 2) if np.isfinite(shift) else None, "control_minus_v68_px_mean": round(repro, 2),
                        "frames": t_off["frames"]}
        say()
    if args.dry_run:
        say("--dry-run: calibration NOT written")
        return points
    target = timex_timing.timex_bias_path("mid", cal_dir)
    about = {
        "what": "Timex bias corrections for TIMING_MODE=mid: (x fraction of the width, correction px) points, "
                "interpolated across columns and SUBTRACTED from the detected row, fitted on detections made "
                "with the temporal blend OFF.",
        "used_by": "waterline_detector_v5.py --timing-mode mid (TIMING_MODE=mid), whose lines "
                   "extract_elevation_contours.py --timing-mode mid gives the water level at mid-exposure "
                   f"(the file-name epoch + {timex_timing.TIMEX_MID_OFFSET_S} s).",
        "derived": f"{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')} by recalibrate_timing.py on "
                   f"{socket.gethostname()}, detector v{detector.VERSION}",
        "method": "the hand-clicked timex ground truth detected with --timing-mode mid --no-bias-correction, "
                  "then derive_bias_correction.py --suffix .timex. --uncorrected; points_from_table: a constant "
                  f"unless spread/SEM > {SPREAD_SEM_CURVE:g}, far-field segments clamped where the error bar is "
                  "larger than the value or few frames reach them",
        "ground_truth_dir": str(gt_dir),
        "c1": decisions["c1"], "c2": decisions["c2"],
        "tables": {cam: tables[("blend_off", cam)] for cam in CAMERAS},
        "control_blend_on": summary,
    }
    backup = write_mid_calibration(target, points, about)
    say(f"WROTE {target}" + (f" (the previous one kept as {backup.name})" if backup else ""))
    say("   (calibration/timex_bias_legacy.json is unchanged.) Commit it with the station's other calibrations.")
    return points


# ---------------------------------------------------------------------
# compare
# ---------------------------------------------------------------------

def window_photos(images_dir, start, end):
    """The timex photos of images_dir taken on the UTC days start..end, both cameras."""
    out = []
    for p in sorted(Path(images_dir).glob("*.timex.jpg")):
        e = epoch_of(p.name)
        if e is not None and camera_of(p.name) and start <= utc_day(e) <= end:
            out.append(p)
    return out


def chain_commands(mode, mdir, cfg, waves, gauge, cal_dir, setup_coef):
    """The cron's steps after detection, on mdir/detections, in timing mode `mode`:
    [(name, command, output)]. setup_coef None: no setup (for --fit-setup)."""
    stem = "" if setup_coef is not None else "c0_"
    contours, qc, ground = mdir / f"{stem}contours.csv", mdir / f"{stem}contours_qc.csv", mdir / f"{stem}ground.csv"
    extract = [sys.executable, HERE / "extract_elevation_contours.py", cfg["GNSSR_SPLINE"],
               "--processed-dir", mdir / "detections", "--timing-mode", mode, "--output", contours]
    if waves:
        extract += ["--waves", waves]
        if setup_coef is not None:
            extract += ["--setup-coef", setup_coef]
    if gauge:
        extract += ["--gnssr-qc-reference", gauge]
    steps = [("contours", extract, contours)]
    if cfg["CONSISTENCY_ENABLE"] == "1":
        steps.append(("consistency", [sys.executable, HERE / "waterline_consistency.py", contours, "--output", qc,
                                      "--report", mdir / f"{stem}consistency_report.csv"], qc))
    geo = [sys.executable, HERE / "georectify.py", "<contours>", ground]
    for cam in CAMERAS:
        geo += [f"--io-{cam}", cal_dir / f"CACO05_{cam}_20240801_IO.yaml",
                f"--eo-{cam}", cal_dir / f"CACO05_{cam}_20251113_EO-CV.yaml"]
    steps.append(("georectify", geo, ground))
    dem = [sys.executable, HERE / "dem_from_contours.py", ground, mdir / f"{stem}dem", "--no-plot",
           "--cell", cfg["DEM_CELL"], "--min-points", cfg["DEM_MIN_POINTS"], "--max-spread", cfg["DEM_MAX_SPREAD"]]
    if waves and cfg["DEM_MAX_HS"]:
        dem += ["--max-hs", cfg["DEM_MAX_HS"]]
    if cfg["DEM_MAX_DAY_OFFSET"]:
        dem += ["--max-day-offset", cfg["DEM_MAX_DAY_OFFSET"], "--trend-note"]
    if setup_coef is None:
        dem += ["--fit-setup"]
    steps.append(("dem", dem, mdir / f"{stem}dem_info.json"))
    return steps


def run_chain(mode, mdir, cfg, waves, gauge, cal_dir, setup_coef, log):
    """Runs chain_commands as the cron does (the consistency filter's failure falls back
    to the unfiltered contours). -> the DEM's info.json path."""
    use = None
    for name, cmd, out in chain_commands(mode, mdir, cfg, waves, gauge, cal_dir, setup_coef):
        cmd = [use if c == "<contours>" else c for c in cmd]
        rc = run(cmd, log)
        if name == "contours":
            if rc != 0 or not Path(out).exists():
                raise RuntimeError(f"{mode}: contour extraction failed (exit {rc}); see {log}")
            use = out
        elif name == "consistency":
            if rc == 0 and Path(out).exists() and Path(out).stat().st_size > 0:
                use = out
            else:
                say(f"   {mode}: consistency filter failed (exit {rc}); the unfiltered contours are used, as the cron does")
        elif rc != 0 or not Path(out).exists():
            raise RuntimeError(f"{mode}: {name} failed (exit {rc}); see {log}")
    return Path(out)


def residual_summary(info):
    """{camera: {...}} of a DEM's diagnostics.residuals.by_camera (cameras c1, c2 and both)."""
    out = {}
    for cam, b in ((info.get("diagnostics") or {}).get("residuals") or {}).get("by_camera", {}).items():
        v = b.get("vs_dzdt") or {}
        st = b.get("by_stage") or {}
        fl = (st.get("flood") or {}).get("median_residual_m")
        eb = (st.get("ebb") or {}).get("median_residual_m")
        out[cam] = {"frames": b.get("frames"), "slope": v.get("slope"), "r": v.get("r"), "lag_min": b.get("lag_min"),
                    "flood_m": fl, "ebb_m": eb, "flood_frames": (st.get("flood") or {}).get("frames"),
                    "ebb_frames": (st.get("ebb") or {}).get("frames"),
                    "contrast_m": None if fl is None or eb is None else round(fl - eb, 4)}
    return out


def judge(legacy, mid):
    """
    The pass rule, per camera (c1, c2): PASS when from legacy to mid BOTH
    |slope of the residual on dz/dt| and |flood median - ebb median| shrink;
    FAIL when either does not; 'not judged' when either mode has no fit or no
    flood or ebb frames. -> ({camera: (verdict, why)}, overall), overall
    PASS (every judged camera passes, at least one judged), FAIL, or
    NOT JUDGED.
    """
    out = {}
    for cam in CAMERAS:
        a, b = legacy.get(cam), mid.get(cam)
        if not a or not b or None in (a["slope"], b["slope"], a["contrast_m"], b["contrast_m"]):
            out[cam] = ("not judged", "too few frames with a rate, or no flood or ebb frames, in one of the modes")
            continue
        s_ok = abs(b["slope"]) < abs(a["slope"])
        c_ok = abs(b["contrast_m"]) < abs(a["contrast_m"])
        why = (f"|slope| {abs(a['slope']):.3f} -> {abs(b['slope']):.3f} m per m/h "
               f"({'shrinks' if s_ok else 'does NOT shrink'}), |flood - ebb| {abs(a['contrast_m']):.3f} -> "
               f"{abs(b['contrast_m']):.3f} m ({'shrinks' if c_ok else 'does NOT shrink'})")
        out[cam] = ("PASS" if s_ok and c_ok else "FAIL", why)
    judged = [v for v, _ in out.values() if v != "not judged"]
    overall = "NOT JUDGED" if not judged else "PASS" if all(v == "PASS" for v in judged) else "FAIL"
    return out, overall


def read_grid(path):
    """{(easting index, northing index): value} of an ESRI ASCII grid's finite cells, and the cell size."""
    from asc_to_geotiff import read_asc
    data, head = read_asc(path)
    cell = head["cellsize"]
    rows, cols = data.shape
    out = {}
    for i, j in zip(*np.nonzero(np.isfinite(data))):
        e = head["xllcorner"] + (j + 0.5) * cell
        n = head["yllcorner"] + (rows - i - 0.5) * cell
        out[(int(round(e / cell)), int(round(n / cell)))] = float(data[i, j])
    return out, cell


def grid_comparison(dem_a, dem_b, spread_a, spread_b):
    """Mid minus legacy over the DEM cells both have, and each DEM's median spread."""
    a, _ = read_grid(dem_a)
    b, _ = read_grid(dem_b)
    both = sorted(set(a) & set(b))
    d = np.array([b[k] - a[k] for k in both])
    sa, _ = read_grid(spread_a)
    sb, _ = read_grid(spread_b)
    return {"cells_legacy": len(a), "cells_mid": len(b), "cells_both": len(both),
            "mid_minus_legacy_median_m": round(float(np.median(d)), 4) if len(d) else None,
            "mid_minus_legacy_p16_p84_m": [round(float(np.percentile(d, 16)), 4),
                                           round(float(np.percentile(d, 84)), 4)] if len(d) else None,
            "spread_median_legacy_m": round(float(np.median(list(sa.values()))), 4) if sa else None,
            "spread_median_mid_m": round(float(np.median(list(sb.values()))), 4) if sb else None}


def fmt(v, f="+.3f"):
    return "--" if v is None else format(v, f)


def report_lines(start, end, res, verdicts, overall, grids, setup):
    lines = [f"Water-level timing test, {start} .. {end}: the DEM's residuals against dz/dt "
             f"(dem_from_contours.py diagnostics)",
             "   slope in m per m/h (x 60 = minutes a line is older than its level); flood/ebb: median residual of",
             "   frames with dz/dt above +0.1 / below -0.1 m/h; flood-ebb: their difference",
             f"   {'camera':<6} {'mode':<7} {'frames':>6} {'slope':>7} {'lag min':>8} {'r':>6} {'flood m':>8} "
             f"{'ebb m':>8} {'flood-ebb':>9}"]
    for cam in CAMERAS + ("both",):
        for mode in ("legacy", "mid"):
            s = res[mode].get(cam)
            if not s:
                continue
            lines.append(f"   {cam:<6} {mode:<7} {s['frames'] or 0:6d} {fmt(s['slope']):>7} {fmt(s['lag_min'], '+.1f'):>8} "
                         f"{fmt(s['r'], '+.2f'):>6} {fmt(s['flood_m']):>8} {fmt(s['ebb_m']):>8} {fmt(s['contrast_m']):>9}")
    lines.append("   expected if the timing error is real: legacy lag about +4 to +8 min and flood above ebb; mid near 0")
    for cam in CAMERAS:
        v, why = verdicts[cam]
        lines.append(f"   {cam}: {v} -- {why}")
    if grids:
        lines.append(f"   DEM: mid minus legacy over {grids['cells_both']} shared cells: median "
                     f"{fmt(grids['mid_minus_legacy_median_m'])} m (16-84%: "
                     + (", ".join(fmt(x) for x in grids["mid_minus_legacy_p16_p84_m"])
                        if grids["mid_minus_legacy_p16_p84_m"] else "--")
                     + f"); median spread legacy {fmt(grids['spread_median_legacy_m'], '.3f')} m, "
                       f"mid {fmt(grids['spread_median_mid_m'], '.3f')} m")
    for mode in ("legacy", "mid"):
        if setup.get(mode):
            lines.append(f"   setup refit ({mode}, contours with C = 0):")
            lines += ["      " + ln for ln in setup[mode]]
    lines.append(f"VERDICT: {overall}" + {
        "PASS": " -- the timing fix removes the tide-stage dependence the legacy DEM shows. Next: --steps backfill, "
                "then TIMING_MODE=mid in waterline_timex_cron.sh (docs/VALIDATION.md section 2.7).",
        "FAIL": " -- keep TIMING_MODE=legacy; see docs/VALIDATION.md section 2.7 for what to look at.",
        "NOT JUDGED": " -- no camera had enough flood and ebb frames with a rate: choose a longer window."}[overall])
    return lines


def setup_fit_lines(log_text):
    """The 'Setup fit' block of a dem_from_contours.py --fit-setup log."""
    out, on = [], False
    for ln in log_text.splitlines():
        if ln.startswith("Setup fit"):
            on = True
        elif on and not ln.strip():
            break
        if on:
            out.append(ln.rstrip())
    return out


def step_compare(args, work, cal_dir, cfg):
    say("=" * 78)
    say(f"COMPARE: legacy against mid over {args.start_date} .. {args.end_date}, the cron's chain")
    say("=" * 78)
    if not (args.start_date and args.end_date):
        raise RuntimeError("compare needs --start-date and --end-date (UTC days of the window)")
    mid_file = timex_timing.timex_bias_path("mid", cal_dir)
    cdir = work / "compare" / f"{args.start_date}_{args.end_date}"
    cdir.mkdir(parents=True, exist_ok=True)
    log = cdir / "compare.log"
    waves = args.waves if args.waves and Path(args.waves).exists() else None
    gauge = args.gauge if args.gauge and Path(args.gauge).exists() else None
    if not Path(cfg["GNSSR_SPLINE"]).exists():
        raise RuntimeError(f"GNSS-R spline not found: {cfg['GNSSR_SPLINE']} (--gnssr)")
    say(f"water level {cfg['GNSSR_SPLINE']}; waves {waves or 'none (no setup, no Hs filter)'}; "
        f"GNSS-R QC {gauge or 'none'}; setup C {cfg['SETUP_COEF'] or 'none'}; log {log}")
    photos = None
    infos, res, setup = {}, {}, {}
    for mode in ("legacy", "mid"):
        mdir = cdir / mode
        if args.reuse_detections and list((mdir / "detections").glob("*.csv")):
            say(f"{mode}: reusing the detections in {mdir / 'detections'}")
        else:
            if photos is None:
                photos = window_photos(args.images, args.start_date, args.end_date)
                if not photos:
                    raise RuntimeError(f"no timex photos of {args.start_date} .. {args.end_date} in {args.images}")
                say(f"{len(photos)} photo(s) in the window ({args.images})")
            if mode == "mid" and not mid_file.exists():
                raise RuntimeError(f"{mid_file} not found: run --steps derive first")
            say(f"{mode}: detecting {len(photos)} photo(s)")
            detect(mode, photos, mdir, log, timex_bias_file=mid_file if mode == "mid" else None,
                   keep_debug=args.keep_debug)
        say(f"{mode}: contours, consistency, georectify, DEM")
        info_path = run_chain(mode, mdir, cfg, waves, gauge, cal_dir, cfg["SETUP_COEF"] or None, log)
        infos[mode] = json.loads(info_path.read_text())
        res[mode] = residual_summary(infos[mode])
        if args.fit_setup:
            say(f"{mode}: setup refit on contours with C = 0")
            before = log.read_text() if log.exists() else ""
            run_chain(mode, mdir, cfg, waves, gauge, cal_dir, None, log)
            setup[mode] = setup_fit_lines(log.read_text()[len(before):])
    verdicts, overall = judge(res["legacy"], res["mid"])
    try:
        grids = grid_comparison(cdir / "legacy" / "dem_dem.asc", cdir / "mid" / "dem_dem.asc",
                                cdir / "legacy" / "dem_spread.asc", cdir / "mid" / "dem_spread.asc")
    except (OSError, ValueError, KeyError) as e:
        say(f"   (DEM comparison skipped: {e})")
        grids = None
    lines = report_lines(args.start_date, args.end_date, res, verdicts, overall, grids, setup)
    say()
    for ln in lines:
        say(ln)
    (cdir / "report.txt").write_text("\n".join(lines) + "\n")
    (cdir / "report.json").write_text(json.dumps(
        {"window": [args.start_date, args.end_date], "residuals": res, "verdicts": verdicts, "overall": overall,
         "grids": grids, "setup_fit": setup, "settings": cfg, "waves": waves, "gauge": gauge,
         "timex_bias_mid": str(mid_file)}, indent=1) + "\n")
    say(f"report: {cdir / 'report.txt'}")
    return overall


# ---------------------------------------------------------------------
# backfill
# ---------------------------------------------------------------------

def step_backfill(args, work, cal_dir):
    say("=" * 78)
    say("BACKFILL: the archive's timex photos detected in mode mid, into archive/processed_timex_mid")
    say("=" * 78)
    mid_file = timex_timing.timex_bias_path("mid", cal_dir)
    if not mid_file.exists():
        raise RuntimeError(f"{mid_file} not found: run --steps derive first")
    target = Path(args.archive) / ("processed_timex" + timex_timing.archive_suffix("mid"))
    target.mkdir(parents=True, exist_ok=True)
    bdir = work / "backfill"
    bdir.mkdir(parents=True, exist_ok=True)
    done_file = bdir / "done_days.txt"
    done = set(done_file.read_text().split()) if done_file.exists() and not args.redo else set()
    days = {}
    for p in sorted(Path(args.images).glob("*.timex.jpg")):
        e = epoch_of(p.name)
        if e is None or not camera_of(p.name):
            continue
        d = utc_day(e)
        if (args.start_date and d < args.start_date) or (args.end_date and d > args.end_date):
            continue
        days.setdefault(d, []).append(p)
    todo = [d for d in sorted(days) if d not in done]
    say(f"{len(days)} day(s) of photos in {args.images}, {len(days) - len(todo)} done before; "
        f"{sum(len(days[d]) for d in todo)} photo(s) to detect -> {target}")
    log = bdir / "backfill.log"
    added = 0
    for k, d in enumerate(todo, 1):
        say(f"   {d} ({k}/{len(todo)}): {len(days[d])} photo(s)")
        dets = detect("mid", days[d], bdir / "run", log, timex_bias_file=mid_file, keep_debug=False)
        for f in dets:
            dst = target / f.name
            if not dst.exists():                        # never overwrite (the cron's cp -pn)
                shutil.copy2(str(f), str(dst))
                added += 1
        with open(done_file, "a") as fh:
            fh.write(d + "\n")
    say(f"backfill: {added} detection(s) added to {target} "
        f"({len(list(target.glob('*.csv')))} there now); days done listed in {done_file}")


# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--steps", default="exposure,derive,compare",
                    help="comma-separated: exposure, derive, compare, backfill (default exposure,derive,compare)")
    ap.add_argument("--start-date", default=None, help="compare/backfill: first UTC day, YYYY-MM-DD")
    ap.add_argument("--end-date", default=None, help="compare/backfill: last UTC day, YYYY-MM-DD")
    ap.add_argument("--base", default=str(BASE), help=f"the station folder (default {BASE})")
    ap.add_argument("--work", default=None, help="work folder (default <base>/timing_recalibration)")
    ap.add_argument("--calibration-dir", default=str(timex_timing.CALIBRATION_DIR),
                    help="where timex_bias_mid.json is written and read, and the IO/EO files are "
                         "(default: this repository's calibration/)")
    ap.add_argument("--ground-truth-dir", default=None,
                    help="hand-clicked ground truth (default <base>/ground_truth, else this repository's)")
    ap.add_argument("--image-dirs", nargs="+", default=None,
                    help="folders searched for the ground truth's photos (default <base>/archive/images_timex, "
                         "<base>/input, <base>/input_timex, /mnt/I2Rgus_Data/ImageProducts)")
    ap.add_argument("--images", default=None,
                    help="compare/backfill: the timex photo archive (default <base>/archive/images_timex)")
    ap.add_argument("--archive", default=None, help="backfill: the archive folder (default <base>/archive)")
    ap.add_argument("--stacks", default=None,
                    help="exposure: the c2 timestack archive (default <base>/archive/ras_c2)")
    ap.add_argument("--ignore-exposure", action="store_true",
                    help="go on to the other steps after an exposure MISMATCH")
    ap.add_argument("--cron", default=str(HERE / "waterline_timex_cron.sh"),
                    help="the cron whose settings the window is processed with")
    ap.add_argument("--gnssr", default=None, help="GNSS-R spline (default: the cron's GNSSR_SPLINE)")
    ap.add_argument("--waves", default=None,
                    help="wave file (default: <base>/archive/waves_marconi.csv as the cron with "
                         "USE_MARCONI_WAVES=1, else the buoy archive)")
    ap.add_argument("--gauge", default=None, help="tide-gauge archive for GNSS-R QC (default: the cron's)")
    ap.add_argument("--min-frames", type=int, default=6,
                    help="derive: fewest ground-truth photos per camera (default 6; v6.8 had 12 and 13)")
    ap.add_argument("--fit-setup", action="store_true",
                    help="compare: also refit the setup coefficient on contours with C = 0, both modes")
    ap.add_argument("--reuse-detections", action="store_true",
                    help="use detections already in the work folder instead of detecting again")
    ap.add_argument("--dry-run", action="store_true", help="derive: print the corrections, write nothing")
    ap.add_argument("--redo", action="store_true", help="backfill: detect days already done again")
    ap.add_argument("--keep-debug", action="store_true", help="keep the detector's debug images")
    ap.add_argument("--nice", type=int, default=10, help="run at this lower priority (default 10; 0 = as is)")
    args = ap.parse_args()

    base = Path(args.base)
    steps = [s.strip() for s in args.steps.split(",") if s.strip()]
    bad = [s for s in steps if s not in ("exposure", "derive", "compare", "backfill")]
    if bad or not steps:
        ap.error(f"--steps: unknown {', '.join(bad) or '(none)'}; give exposure, derive, compare and/or backfill")
    if args.nice and hasattr(os, "nice"):
        os.nice(args.nice)
    work = Path(args.work) if args.work else base / "timing_recalibration"
    work.mkdir(parents=True, exist_ok=True)
    cal_dir = Path(args.calibration_dir)
    if args.ground_truth_dir is None:
        args.ground_truth_dir = str(base / "ground_truth") if (base / "ground_truth").is_dir() \
            else str(HERE / "ground_truth")
    if args.image_dirs is None:
        args.image_dirs = [str(base / "archive" / "images_timex"), str(base / "input"),
                           str(base / "input_timex"), "/mnt/I2Rgus_Data/ImageProducts"]
    args.images = args.images or str(base / "archive" / "images_timex")
    args.archive = args.archive or str(base / "archive")
    args.stacks = args.stacks or str(base / "archive" / "ras_c2")
    cfg = cron_settings(args.cron)
    if args.gnssr:
        cfg["GNSSR_SPLINE"] = args.gnssr
    if args.waves is None:
        marconi = base / "archive" / "waves_marconi.csv"
        buoy = base / "archive" / f"waves_{cfg['WAVE_BUOY']}.csv"
        args.waves = str(marconi) if cfg["USE_MARCONI_WAVES"] == "1" and marconi.exists() else str(buoy)
    if args.gauge is None:
        args.gauge = str(base / "archive" / f"gauge_{cfg['TIDE_GAUGE']}.csv")
    # As the cron: with the waves at Marconi in use, the GNSS-R QC's wave allowances use them too
    # (gnssr_qc.default_waves_path reads CACO_LOCAL_WAVES, which the cron exports)
    if cfg["USE_MARCONI_WAVES"] == "1" and Path(args.waves).name == "waves_marconi.csv" \
            and Path(args.waves).exists():
        os.environ["CACO_LOCAL_WAVES"] = "1"
    say(f"recalibrate_timing.py: steps {', '.join(steps)}; work folder {work}; calibration {cal_dir}")
    say(f"   mid-exposure = file-name epoch + {timex_timing.TIMEX_MID_OFFSET_S} s (timex_timing.TIMEX_MID_OFFSET_S)")
    rc = 0
    try:
        if "exposure" in steps:
            if step_exposure(args) == "MISMATCH" and not args.ignore_exposure and steps != ["exposure"]:
                say("STOPPED before " + ", ".join(s for s in steps if s != "exposure")
                    + ": the exposure does not match TIMEX_MID_OFFSET_S (--ignore-exposure goes on)")
                return 4
            say()
        if "derive" in steps:
            step_derive(args, work, cal_dir)
        if "compare" in steps:
            overall = step_compare(args, work, cal_dir, cfg)
            rc = {"PASS": 0, "FAIL": 2, "NOT JUDGED": 3}[overall]
        if "backfill" in steps:
            step_backfill(args, work, cal_dir)
    except (RuntimeError, ValueError, OSError) as e:
        say(f"STOPPED: {e}")
        return 1
    return rc


if __name__ == "__main__":
    sys.exit(main())

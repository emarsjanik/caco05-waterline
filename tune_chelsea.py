#!/usr/bin/env python3
"""
Fine-Tune The Detector On The Chelsea (pre-Nov-2025) Imagery
===============================================================
Ground truth in, measured accuracy and a fitted correction out.

The detector's settings -- search envelope, signal threshold, and above
all the per-column bias correction -- were tuned on 2026 imagery. The
Chelsea frames are winter 2025, taken with the 2025-02-19 camera
pointing and resampled into the 2025-11-13 view (process_chelsea.py).
Whether those settings suit them is measured here, not assumed.

    python3 tune_chelsea.py sample   pick frames to trace
    python3 collect_all_ground_truth.py ...   (you trace them; the command is printed)
    python3 tune_chelsea.py score    run detector settings on those frames, score each
    python3 tune_chelsea.py derive   fit a Chelsea bias correction, tested day by day
    python3 tune_chelsea.py diagnose where and why a scored setting goes wrong

SAMPLE. Frames come from <out>/work/src -- the RESAMPLED frames the
detector actually sees -- so the traced line and the detector's line are
in the same pixels. Per camera, frames are spread evenly over the water
levels (low to high tide) and rotated over the days, so the correction
is not fitted to one tide stage or one morning. Traced on the resampled
view, the blurred margins are areas the old camera never saw: trace only
where the photo is sharp.

SCORE. Each setting runs the detector on the traced frames only (a few
minutes), then per camera reports: frames detected, coverage (share of
traced columns the detector also found), and the error of the detected
line -- in pixels (median |error|, bias, RMS, 90th percentile) and in
METRES on the ground (both lines georectified at the frame's water
level), which is what matters for the DEM. Only columns the detector
flagged as carrying signal are scored, as only those reach the DEM.

DERIVE. From the setting run WITHOUT any bias correction ("no_bias"),
the median signed error in 10 column segments per camera is the
correction. Fitted on all days and saved as JSON for the detector
(--bias-correction-file, also passed by process_chelsea.py). Its benefit
is reported leave-one-day-out -- each day corrected with a fit to the
OTHER days -- so the number is what to expect on frames it has not
seen, not a fit to itself.

DIAGNOSE. For one scored setting (default "current"): each traced
frame's signal fraction and whether the detector kept or discarded it
(and why), with its error; then the error by column segment (tenths of
the image width, left to right). Separates a correctable, consistent
offset from scatter, and a threshold problem from a real absence of
signal.

FILTERS. score and diagnose take --column-limits "c2=0:0.7" (score only
those columns, exactly what the detector's --column-limits keeps) and
--utc-hours "13.5-18" (score only those frames), so a cut can be judged
on runs already made.
"""

import os
import csv
import sys
import json
import shutil
import argparse
import subprocess
from pathlib import Path
from collections import defaultdict

import numpy as np

HERE = Path(__file__).resolve().parent
CAL = HERE / "calibration"
WIDTH = 2448
N_SEGMENTS = 10

BASE = ["--image-suffix", "timex.jpg", "--min-signal-fraction", "0.35"]
CONFIGS = {
    "current":  BASE + ["--envelope-pad", "0.1"],                         # as process_chelsea runs
    "no_bias":  BASE + ["--envelope-pad", "0.1", "--no-bias-correction"],
    "pad_0.05": BASE + ["--envelope-pad", "0.05", "--no-bias-correction"],
    "pad_0.2":  BASE + ["--envelope-pad", "0.2", "--no-bias-correction"],
    # lower signal thresholds: more frames kept -- at what cost in error?
    "signal_0.25": ["--image-suffix", "timex.jpg", "--min-signal-fraction", "0.25", "--envelope-pad", "0.1"],
    "signal_0.15": ["--image-suffix", "timex.jpg", "--min-signal-fraction", "0.15", "--envelope-pad", "0.1"],
    # c2 frames with 50-75% signal were thrown out by the 0.9 confidence gate
    "conf_0.75": ["--image-suffix", "timex.jpg", "--min-signal-fraction", "0.15", "--envelope-pad", "0.1",
                  "--min-confidence", "0.75"],
    "conf_0.6":  ["--image-suffix", "timex.jpg", "--min-signal-fraction", "0.15", "--envelope-pad", "0.1",
                  "--min-confidence", "0.6"],
}


def parse_limits(spec):
    """"c2=0:0.7,c1=0:0.9" -> {"c2": (0.0, 0.7), ...}; None -> {}."""
    out = {}
    for item in (spec or "").split(","):
        if item.strip():
            cam, rng = item.split("=")
            lo, hi = (float(v) for v in rng.split(":"))
            out[cam.strip().lower()] = (lo, hi)
    return out


def frame_hour(name):
    """UTC hour (fractional) from the epoch at the start of an Argus file name."""
    from datetime import datetime, timezone
    t = datetime.fromtimestamp(int(name.split(".")[0]), tz=timezone.utc)
    return t.hour + t.minute / 60


def frame_levels(out):
    levels = {}
    with open(out / "work" / "frame_levels.csv", newline="") as f:
        for r in csv.DictReader(f):
            levels[r["frame"]] = (r["camera"], float(r["water_level_navd88"]))
    return levels


def day_of(name):
    return name.split(".")[2] + "." + name.split(".")[3].split("_")[0]


# --------------------------------------------------------------------- sample

def cmd_sample(args):
    out = Path(args.out)
    levels = frame_levels(out)
    src = out / "work" / "src"
    dst = out / "ground_truth_images"
    dst.mkdir(parents=True, exist_ok=True)
    picked = []
    for cam in ("c1", "c2"):
        frames = sorted((n for n, (c, _) in levels.items() if c == cam and (src / n).exists()),
                        key=lambda n: levels[n][1])
        if not frames:
            continue
        bins = np.array_split(np.array(frames), min(args.per_camera, len(frames)))
        used = defaultdict(int)
        for b in bins:
            # the frame whose day is least used so far, nearest the bin's middle
            mid = len(b) // 2
            best = min(range(len(b)), key=lambda i: (used[day_of(b[i])], abs(i - mid)))
            used[day_of(b[best])] += 1
            picked.append(b[best])
    for n in picked:
        target = dst / n
        if not target.exists():
            try:
                os.link(src / n, target)
            except OSError:
                shutil.copy2(src / n, target)
    by = defaultdict(list)
    for n in picked:
        by[levels[n][0]].append(levels[n][1])
    for cam, z in sorted(by.items()):
        print(f"{cam}: {len(z)} frames, water level {min(z):+.2f} to {max(z):+.2f} m NAVD88, "
              f"{len({day_of(n) for n in picked if levels[n][0] == cam})} day(s)")
    print(f"\nFrames in {dst}. Trace them (MobaXterm opens the window over X11):")
    print(f"  python3 {HERE / 'collect_all_ground_truth.py'} --input-dir {dst} "
          f"--ground-truth-dir {out / 'ground_truth'}")
    print("Click along the water's edge, left to right, only where the photo is sharp "
          "(not in the blurred margins); Save & Next. It resumes where you stop.")


# --------------------------------------------------------------------- score

def load_gt(path):
    cols, rows = [], []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            cols.append(int(r["Column"]))
            rows.append(float(r["Row"]))
    return dict(zip(cols, rows))


def load_det(path):
    det = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if r.get("Has_Signal", "1") not in ("1", "1.0", "True", ""):
                continue
            try:
                row = float(r.get("Row_Precise") or r["Row"])
            except ValueError:
                row = float(r["Row"])
            det[int(r["Column"])] = row
    return det


def calib():
    from georectify import load_extrinsics, load_intrinsics
    return {cam: (load_intrinsics(CAL / f"CACO05_{cam}_20240801_IO.yaml"),
                  load_extrinsics(CAL / f"CACO05_{cam}_20251113_EO-CV.yaml")) for cam in ("c1", "c2")}


def ground_error(cols, gt_rows, det_rows, z, io, eo):
    """Horizontal distance (m) between the traced and detected lines, per column."""
    from georectify import pixel_to_ground
    gx, gy = pixel_to_ground(cols, gt_rows, z, io, eo)
    dx, dy = pixel_to_ground(cols, det_rows, z, io, eo)
    return np.hypot(dx - gx, dy - gy)


def errors_for(det_dir, out, levels, cal, limits=None, hours=None):
    """Per traced frame: columns, gt rows, det rows (np), z, camera. Missing detection -> None.
    limits: {"c2": (lo, hi)} keeps only detected columns in that share of the width (as the
    detector's --column-limits does); hours: (h0, h1) UTC keeps only frames in that range."""
    res = {}
    for gt_path in sorted((out / "ground_truth").glob("*.csv")):
        name = gt_path.stem + ".jpg"
        if name not in levels:
            continue
        if hours and not hours[0] <= frame_hour(name) < hours[1]:
            continue
        cam, z = levels[name]
        gt = load_gt(gt_path)
        det_path = det_dir / gt_path.name
        if not det_path.exists():
            res[name] = dict(cam=cam, z=z, n_gt=len(gt), cols=None)
            continue
        det = load_det(det_path)
        if limits and cam in limits:
            lo, hi = limits[cam]
            det = {c: r for c, r in det.items() if lo <= c / (WIDTH - 1) <= hi}
        common = sorted(set(gt) & set(det))
        res[name] = dict(cam=cam, z=z, n_gt=len(gt), cols=np.array(common, float),
                         gt=np.array([gt[c] for c in common]), det=np.array([det[c] for c in common]))
    return res


def summarise(res, cal, label, lines):
    for cam in ("c1", "c2"):
        fr = [v for v in res.values() if v["cam"] == cam]
        if not fr:
            continue
        found = [v for v in fr if v["cols"] is not None and len(v["cols"])]
        cov = sum(len(v["cols"]) for v in found) / max(sum(v["n_gt"] for v in fr), 1)
        if not found:
            lines.append(f"  {label:10s} {cam}: 0/{len(fr)} frames detected")
            continue
        e = np.concatenate([v["det"] - v["gt"] for v in found])
        io, eo = cal[cam]
        m = np.concatenate([ground_error(v["cols"], v["gt"], v["det"], v["z"], io, eo) for v in found])
        m = m[np.isfinite(m)]
        lines.append(f"  {label:10s} {cam}: {len(found):2d}/{len(fr):2d} frames, coverage {100 * cov:3.0f}%  "
                     f"|err| median {np.median(np.abs(e)):5.1f} px, bias {e.mean():+6.1f}, RMS "
                     f"{np.sqrt(np.mean(e ** 2)):5.1f}, P90 {np.percentile(np.abs(e), 90):5.1f} px  |  "
                     f"ground median {np.median(m):4.1f} m, P90 {np.percentile(m, 90):4.1f} m")


def run_detector(cfg_args, out, name, log):
    d = out / "tuning" / name
    for sub in ("in", "det", "debug"):
        shutil.rmtree(d / sub, ignore_errors=True)
        (d / sub).mkdir(parents=True)
    cmd = [sys.executable, str(HERE / "waterline_detector_v5.py")] + cfg_args + [
        "--source-dir", str(out / "ground_truth_images"), "--input-dir", str(d / "in"),
        "--output-dir", str(d / "det"), "--debug-dir", str(d / "debug")]
    with open(log, "a") as fh:
        fh.write(f"\n=== {name}: {' '.join(cmd)}\n")
        fh.flush()
        ok = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode == 0
    return d / "det" if ok else None


def parse_hours(spec):
    return tuple(float(v) for v in spec.split("-")) if spec else None


def filters_note(limits, hours):
    parts = [f"{c} columns {lo:.2f}-{hi:.2f} of the width" for c, (lo, hi) in sorted(limits.items())]
    if hours:
        parts.append(f"frames {hours[0]:g}-{hours[1]:g} h UTC")
    return f" [{'; '.join(parts)}]" if parts else ""


def cmd_score(args):
    out = Path(args.out)
    if not list((out / "ground_truth").glob("*.csv")):
        sys.exit(f"No traced frames in {out / 'ground_truth'} -- run 'sample' and trace first.")
    levels, cal = frame_levels(out), calib()
    configs = {n: CONFIGS[n] for n in (args.configs or CONFIGS)}
    if args.bias_correction_file:
        configs["chelsea_bias"] = BASE + ["--envelope-pad", "0.1", "--bias-correction-file",
                                          args.bias_correction_file]
    log = out / "tuning" / "tuning.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    n_gt = len(list((out / "ground_truth").glob("*.csv")))
    limits, hours = parse_limits(args.column_limits), parse_hours(args.utc_hours)
    lines = [f"Detector settings scored on {n_gt} traced frame(s)" + filters_note(limits, hours)]
    for name, cfg in configs.items():
        print(f"running {name} ...", flush=True)
        det_dir = run_detector(cfg, out, name, log)
        if det_dir is None:
            lines.append(f"  {name}: detector failed, see {log}")
            continue
        summarise(errors_for(det_dir, out, levels, cal, limits, hours), cal, name, lines)
    report = "\n".join(lines)
    (out / "tuning" / "score_report.txt").write_text(report + "\n")
    print("\n" + report)
    print(f"\nFull detector output: {log}")


# --------------------------------------------------------------------- derive

def fit_correction(frames):
    """Median signed error (px) in N_SEGMENTS column segments -> [(x_fraction, correction)]."""
    x = np.concatenate([v["cols"] for v in frames]) / (WIDTH - 1)
    e = np.concatenate([v["det"] - v["gt"] for v in frames])
    edges = np.linspace(0, 1, N_SEGMENTS + 1)
    pts = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (x >= a) & (x <= b)
        if m.sum() >= 20:
            pts.append((round(float((a + b) / 2), 4), round(float(np.median(e[m])), 1)))
    if not pts:
        return []
    return [(0.0, pts[0][1])] + pts + [(1.0, pts[-1][1])]


def apply(points, v):
    if not points:
        return v["det"]
    xs, cs = zip(*points)
    return v["det"] - np.interp(v["cols"] / (WIDTH - 1), xs, cs)


def cmd_derive(args):
    out = Path(args.out)
    det_dir = out / "tuning" / "no_bias" / "det"
    if not det_dir.exists():
        sys.exit("Run 'score' first (it runs the no_bias setting this needs).")
    levels, cal = frame_levels(out), calib()
    res = {k: v for k, v in errors_for(det_dir, out, levels, cal).items()
           if v["cols"] is not None and len(v["cols"])}
    result, lines = {}, ["Chelsea bias correction (from the no_bias run)"]
    for cam in ("c1", "c2"):
        fr = {k: v for k, v in res.items() if v["cam"] == cam}
        if len(fr) < 4:
            lines.append(f"  {cam}: only {len(fr)} detected traced frame(s) -- not fitted")
            continue
        pts = fit_correction(list(fr.values()))
        result[cam] = pts
        # leave-one-day-out: correct each day with a fit to the other days
        days = sorted({day_of(k) for k in fr})
        before, after, mb, ma = [], [], [], []
        io, eo = cal[cam]
        for d in days:
            train = [v for k, v in fr.items() if day_of(k) != d]
            test = [v for k, v in fr.items() if day_of(k) == d]
            p = fit_correction(train) if len(train) >= 3 else []
            for v in test:
                corrected = apply(p, v)
                before.append(v["det"] - v["gt"])
                after.append(corrected - v["gt"])
                mb.append(ground_error(v["cols"], v["gt"], v["det"], v["z"], io, eo))
                ma.append(ground_error(v["cols"], v["gt"], corrected, v["z"], io, eo))
        b, a = np.concatenate(before), np.concatenate(after)
        gb, ga = np.concatenate(mb), np.concatenate(ma)
        lines += [f"  {cam}: {len(fr)} frames over {len(days)} day(s); correction (x fraction, px): "
                  + ", ".join(f"{x:.2f}:{c:+.0f}" for x, c in pts[1:-1]),
                  f"      leave-one-day-out  median |err| {np.median(np.abs(b)):.1f} -> "
                  f"{np.median(np.abs(a)):.1f} px, bias {b.mean():+.1f} -> {a.mean():+.1f} px, "
                  f"ground median {np.nanmedian(gb):.1f} -> {np.nanmedian(ga):.1f} m"]
        if len(days) < 3:
            lines.append("      (fewer than 3 days: the leave-one-day-out estimate is weak)")
    path = out / "chelsea_bias_correction.json"
    path.write_text(json.dumps(result, indent=1) + "\n")
    lines += ["", f"Saved {path}.",
              "Check it:   python3 tune_chelsea.py score --configs current no_bias "
              f"--bias-correction-file {path}",
              "Use it:     python3 process_chelsea.py ... --bias-correction-file " + str(path)]
    report = "\n".join(lines)
    (out / "tuning" / "derive_report.txt").write_text(report + "\n")
    print(report)


# --------------------------------------------------------------------- diagnose

def parse_log(log, name):
    """frame -> {"signal": % of columns with signal, "status": kept / why discarded},
    from the LAST run of setting `name` in tuning.log."""
    if not log.exists():
        return {}
    text = "\n" + log.read_text(errors="replace")
    start = text.rfind(f"\n=== {name}: ")
    if start < 0:
        return {}
    end = text.find("\n=== ", start + 5)
    info, cur = {}, None
    for line in text[start:end if end > 0 else None].splitlines():
        line = line.strip()
        if line.startswith("Processing:"):
            cur = line.split(":", 1)[1].strip()
            info[cur] = {"signal": None, "status": "kept"}
        elif cur and line.startswith("Signal columns"):
            try:
                info[cur]["signal"] = float(line.split(":", 1)[1].split("%")[0])
            except ValueError:
                pass
        elif cur and line.startswith("DISCARDED"):
            why = line[len("DISCARDED"):].strip()
            info[cur]["status"] = "DISCARDED " + why[1:why.find(")")] if why.startswith("(") else why
    return info


def cmd_diagnose(args):
    out = Path(args.out)
    det_dir = out / "tuning" / args.config / "det"
    if not det_dir.exists():
        sys.exit(f"Run 'score' first (no {args.config} run in {out / 'tuning'}).")
    levels, cal = frame_levels(out), calib()
    limits, hours = parse_limits(args.column_limits), parse_hours(args.utc_hours)
    res = errors_for(det_dir, out, levels, cal, limits, hours)
    log = parse_log(out / "tuning" / "tuning.log", args.config)
    lines = [f"Diagnosis of the '{args.config}' run ({len(res)} traced frames)"
             + filters_note(limits, hours)]
    summarise(res, cal, args.config, lines)
    for cam in ("c1", "c2"):
        fr = {k: v for k, v in res.items() if v["cam"] == cam}
        if not fr:
            continue
        io, eo = cal[cam]
        lines += ["", f"{cam} -- per frame (signal = share of columns the detector finds signal in)",
                  f"  {'frame (UTC)':16s} {'z m':>6s} {'signal':>7s}  {'scored':>9s}  {'|err| px':>8s} "
                  f"{'bias px':>7s} {'ground m':>8s}  status"]
        for k in sorted(fr, key=lambda k: fr[k]["z"]):
            v = fr[k]
            lg = log.get(k, {})
            sig = f"{lg['signal']:5.1f}%" if lg.get("signal") is not None else "     ?"
            status = lg.get("status", "?" if v["cols"] is None else "kept")
            if v["cols"] is not None and len(v["cols"]):
                e = v["det"] - v["gt"]
                g = ground_error(v["cols"], v["gt"], v["det"], v["z"], io, eo)
                err = (f"{len(v['cols']):4d}/{v['n_gt']:<4d}  {np.median(np.abs(e)):8.1f} "
                       f"{np.median(e):+7.1f} {np.nanmedian(g):8.1f}")
            else:
                err = f"{'0':>4s}/{v['n_gt']:<4d}  {'-':>8s} {'-':>7s} {'-':>8s}"
            parts = k.split(".")
            label = ".".join(parts[2:4]) if len(parts) > 4 else k[:16]
            lines.append(f"  {label[:16]:16s} {v['z']:+6.2f} {sig:>7s}  {err}  {status}")
        found = [v for v in fr.values() if v["cols"] is not None and len(v["cols"])]
        lines += ["", f"{cam} -- by column segment (0.0 = left edge, 1.0 = right edge)",
                  f"  {'segment':9s} {'traced':>6s} {'scored':>6s}  {'|err| px':>8s} {'bias px':>7s} "
                  f"{'spread px':>9s}  {'ground m':>8s} {'P90 m':>6s}"]
        if not found:
            lines.append("  (no detected frames)")
            continue
        x = np.concatenate([v["cols"] for v in found]) / (WIDTH - 1)
        e = np.concatenate([v["det"] - v["gt"] for v in found])
        g = np.concatenate([ground_error(v["cols"], v["gt"], v["det"], v["z"], io, eo) for v in found])
        xt = np.concatenate([np.array(list(load_gt(out / "ground_truth" / (Path(k).stem + ".csv"))), float)
                             for k in fr]) / (WIDTH - 1)
        edges = np.linspace(0, 1, N_SEGMENTS + 1)
        for a, b in zip(edges[:-1], edges[1:]):
            m = (x >= a) & (x < b) if b < 1 else (x >= a)
            nt = int(((xt >= a) & (xt < b)).sum() if b < 1 else (xt >= a).sum())
            if m.sum() < 5:
                lines.append(f"  {a:.1f}-{b:.1f}   {nt:6d} {int(m.sum()):6d}")
                continue
            q25, q75 = np.percentile(e[m], [25, 75])
            gm = g[m][np.isfinite(g[m])]
            lines.append(f"  {a:.1f}-{b:.1f}   {nt:6d} {int(m.sum()):6d}  {np.median(np.abs(e[m])):8.1f} "
                         f"{np.median(e[m]):+7.1f} {q75 - q25:9.1f}  "
                         f"{np.median(gm) if len(gm) else np.nan:8.1f} "
                         f"{np.percentile(gm, 90) if len(gm) else np.nan:6.1f}")
    lines += ["", "bias = detected minus traced row (negative: detected line is HIGHER in the image); "
              "spread = interquartile range of the signed error (large spread = scatter, which a "
              "bias correction cannot fix)."]
    report = "\n".join(lines)
    (out / "tuning" / f"diagnose_{args.config}.txt").write_text(report + "\n")
    print(report)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--out", default="/mnt/I2Rgus_Data/Chelsea_calibration",
                    help="process_chelsea.py output folder")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("sample")
    s.add_argument("--per-camera", type=int, default=20)
    sc = sub.add_parser("score")
    sc.add_argument("--configs", nargs="*", choices=list(CONFIGS), default=None)
    sc.add_argument("--bias-correction-file", default=None,
                    help="also score the detector with this correction (e.g. from 'derive')")
    sub.add_parser("derive")
    dg = sub.add_parser("diagnose")
    dg.add_argument("--config", default="current", help="a setting already run by 'score'")
    for p in (sc, dg):
        p.add_argument("--column-limits", default=None,
                       help='score only these columns, e.g. "c2=0:0.7" (as the detector option '
                            'of that name would keep); no detector rerun needed')
        p.add_argument("--utc-hours", default=None,
                       help='score only frames in this UTC range, e.g. "13.5-18"')
    args = ap.parse_args()
    {"sample": cmd_sample, "score": cmd_score, "derive": cmd_derive,
     "diagnose": cmd_diagnose}[args.cmd](args)


if __name__ == "__main__":
    main()

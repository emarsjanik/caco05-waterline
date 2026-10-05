#!/usr/bin/env python3
"""
Wave Period From Camera Timestacks
====================================
Measures the incident wave period at Marconi from the 2 Hz, 10-minute
pixel timestacks (ras.tiff) the c2 camera already records across the surf
zone. Every breaking wave brightens the pixels it passes, so the spectrum
of pixel intensity peaks at the wave period (Stockdon & Holman 2000;
Lippmann & Holman 1991). A single image cannot do this: Chris Sherwood's
period model scored R2 -0.19 on held-back Marconi frames.

HOW. Each timestack pixel is put on the map at that moment's water level
(camera calibration; water level from marconi_water_level.py, GPS first)
and only the SURF ZONE is used: from --margin (15 m) seaward of the
waterline out to --max-seaward. The swash at the waterline moves at the
slower uprush/backwash rhythm, so the first version, which took the most
variable pixels anywhere on the line, read 13-16 s while the buoys said
6-8 s; on a synthetic stack with 40 s swash and 8 s breakers it read 20 s
where this reads 8.0 s. Of the surf-zone pixels the more variable half is
detrended, and their spectra (Welch, 128 s Hann windows, 50% overlap)
normalised and averaged. In the incident band (3-25 s) the peak gives Tp
(parabolic refinement) and the first moment Tm01. The share of variance
below 0.04 Hz (periods over 25 s) is reported as infragravity fraction.
A line with fewer than --min-pixels surf-zone pixels is logged as
"no surf pixels" (the line does not reach past the swash at that tide).

OUTPUT archive/wave_period_<cam>.csv, one row per timestack and line:
  epoch (burst middle), time_utc, filename, line, tp_s, tm01_s,
  ig_fraction, peak_ratio, status (ok / weak peak / dark / no surf pixels),
  n_pixels, surf_from_m, surf_to_m (metres past the waterline), water_level_m

Usage:
    python3 timestack_wave_period.py                       new stacks in archive/ras_c2
    python3 timestack_wave_period.py --reprocess           every archived stack again
    python3 timestack_wave_period.py --compare             camera periods vs buoys 44008/44013
    python3 timestack_wave_period.py --self-test
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import cv2

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
SAMPLE_HZ = 2.0
BURST_MID_OFFSET_S = 300.0
FIELDS = ["epoch", "time_utc", "filename", "line", "tp_s", "tm01_s", "ig_fraction", "peak_ratio", "status",
          "n_pixels", "surf_from_m", "surf_to_m", "water_level_m"]
# beach geometry from the Jan 2025 lidar: the 0 m NAVD88 contour (sea_patch.py)
SHORE_ORIGIN = (420150.0, 4638400.0)
SEAWARD_DEG = 81.0
FORESHORE_SLOPE = 0.08


def welch(x, fs, nseg=256):
    """Mean periodogram of 50%-overlapping Hann segments; x is (time, n) -> (f, psd (n_f, n))."""
    w = np.hanning(nseg)
    step = nseg // 2
    segs = [x[i:i + nseg] for i in range(0, len(x) - nseg + 1, step)]
    if not segs:
        raise ValueError("series shorter than one segment")
    P = 0
    for s in segs:
        s = s - s.mean(axis=0)
        P = P + np.abs(np.fft.rfft(s * w[:, None], axis=0)) ** 2
    P = P / len(segs) / (fs * (w ** 2).sum())
    return np.fft.rfftfreq(nseg, 1 / fs), P


def period_from_stack(gray, cols=None, fs=SAMPLE_HZ, band=(1 / 25, 1 / 3), ig_cut=0.04, min_mean=25.0):
    """
    gray: (time, pixels) along one line -> dict of tp_s, tm01_s, ig_fraction, peak_ratio, status.
    cols: the pixels to use -- the surf zone seaward of the swash (main() picks them on
    the map). Without it, the seaward half by variance (synthetic tests only).
    """
    gray = gray.astype(np.float32)
    if gray.mean() < min_mean:
        return {"status": "dark"}
    std = gray.std(axis=0)
    if cols is None:
        half = gray.shape[1] // 2
        sea = slice(0, half) if std[:half].mean() >= std[half:].mean() else slice(half, None)
        cols = np.arange(gray.shape[1])[sea]
    cols = np.asarray(cols)
    s = std[cols]
    cols = cols[s >= np.median(s)]
    x = gray[:, cols]
    x = x - x.mean(axis=0)
    f, P = welch(x, fs)
    P = P / np.maximum(P.sum(axis=0, keepdims=True), 1e-12)        # each pixel counts equally
    S = P.mean(axis=1)
    inc = (f >= band[0]) & (f <= band[1])
    if inc.sum() < 3:
        return {"status": "weak peak"}
    k = np.flatnonzero(inc)[np.argmax(S[inc])]
    fp = f[k]
    if 0 < k < len(f) - 1:                                       # parabolic refinement
        a, b, c = np.log(S[k - 1:k + 2] + 1e-20)
        den = a - 2 * b + c
        if den < 0:                                              # stay between the neighbours
            fp = f[k] + float(np.clip(0.5 * (a - c) / den, -0.5, 0.5)) * (f[1] - f[0])
    tm01 = (S[inc].sum() / (f[inc] * S[inc]).sum())
    ig = S[(f > 0) & (f < ig_cut)].sum() / max(S[f > 0].sum(), 1e-12)
    ratio = S[k] / max(np.median(S[inc]), 1e-12)
    return {"tp_s": round(1 / fp, 2), "tm01_s": round(float(tm01), 2), "ig_fraction": round(float(ig), 3),
            "peak_ratio": round(float(ratio), 2), "status": "ok" if ratio >= 3 else "weak peak"}


def self_test():
    """Synthetic surf: bright breaking crests every 9.0 s moving shoreward, plus noise and a slow swash."""
    rng = np.random.default_rng(0)
    t = np.arange(1200) / SAMPLE_HZ
    x = np.arange(300)
    worst = 0.0
    for T in (6.0, 9.0, 13.0):
        phase = 2 * np.pi * (t[:, None] / T + x[None, :] / 40.0)
        crest = np.clip(np.cos(phase), 0, None) ** 4 * 120 * (x[None, :] < 200)
        swash = 15 * np.sin(2 * np.pi * t / 90.0)[:, None]
        img = 80 + crest + swash + rng.normal(0, 8, (len(t), len(x)))
        r = period_from_stack(np.clip(img, 0, 255))
        worst = max(worst, abs(r["tp_s"] - T))
        print(f"  true period {T:4.1f} s -> Tp {r['tp_s']:.2f} s, Tm01 {r['tm01_s']:.2f} s, "
              f"IG {r['ig_fraction']:.2f}, peak/median {r['peak_ratio']:.1f} ({r['status']})")
    ok = worst < 0.3
    print(("PASS" if ok else "FAIL") + f": worst error {worst:.2f} s")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--ras-dir", default=str(HERE / "archive" / "ras_c2"))
    ap.add_argument("--camera", default="c2")
    ap.add_argument("--pix", default=None, help="default /home/argus_user/arguseyes/build/<cam>_timestack.pix")
    ap.add_argument("--lines", nargs="+", type=int, default=[1, 2], help="timestack lines (default 1 2)")
    ap.add_argument("--output", default=None, help="default archive/wave_period_<cam>.csv")
    ap.add_argument("--margin", type=float, default=15.0,
                    help="m seaward of the waterline where the surf zone starts (skips the swash; default 15)")
    ap.add_argument("--max-seaward", type=float, default=300.0, help="m past the waterline to stop (default 300)")
    ap.add_argument("--min-pixels", type=int, default=10, help="surf-zone pixels needed on a line (default 10)")
    ap.add_argument("--reprocess", action="store_true", help="recompute every archived stack")
    ap.add_argument("--compare", action="store_true", help="camera periods vs the buoys' periods")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    import pandas as pd
    if args.compare:
        return compare(args)
    from runup_from_timestack import split_lines, epoch_from_name, DEFAULT_PIX_DIR
    import georectify
    from marconi_water_level import WaterLevel
    pix = np.loadtxt(args.pix or f"{DEFAULT_PIX_DIR}/{args.camera}_timestack.pix")[:, :2]
    lines = split_lines(pix)
    io = georectify.load_intrinsics(HERE / "calibration" / f"CACO05_{args.camera}_20240801_IO.yaml")
    eo = georectify.load_extrinsics(HERE / "calibration" / f"CACO05_{args.camera}_20251113_EO-CV.yaml")
    wl = WaterLevel()
    sea = np.radians(SEAWARD_DEG)
    out = Path(args.output or HERE / "archive" / f"wave_period_{args.camera}.csv")
    done = (pd.read_csv(out) if out.exists() and not args.reprocess
            else pd.DataFrame(columns=FIELDS))
    seen = set(done["filename"])
    rows = []
    for p in sorted(Path(args.ras_dir).glob(f"*.{args.camera}.ras.tiff")):
        if p.name in seen:
            continue
        e = epoch_from_name(p.name)
        ras = cv2.imread(str(p), cv2.IMREAD_UNCHANGED)
        if e is None or ras is None or ras.shape[1] != len(pix):
            continue
        gray = cv2.cvtColor(ras, cv2.COLOR_BGR2GRAY) if ras.ndim == 3 else ras
        mid = e + BURST_MID_OFFSET_S
        z = float(wl.at([mid])[0][0])
        if not np.isfinite(z):
            z = 0.0
        waterline = -z / FORESHORE_SLOPE           # 0 m contour moves seaward as the tide falls
        for ln in args.lines:
            if not 1 <= ln <= len(lines):
                continue
            a, b = lines[ln - 1]
            # each pixel on the map at this moment's water level; keep the surf zone:
            # seaward of the swash (waterline + margin), not out to the horizon
            E, N = georectify.pixel_to_ground(pix[a:b, 0], pix[a:b, 1], z, io, eo)
            sw = (np.asarray(E) - SHORE_ORIGIN[0]) * np.sin(sea) + (np.asarray(N) - SHORE_ORIGIN[1]) * np.cos(sea)
            past = sw - waterline
            cols = np.flatnonzero(np.isfinite(past) & (past > args.margin) & (past < args.max_seaward))
            base = {"epoch": int(mid), "time_utc": datetime.fromtimestamp(mid, tz=timezone.utc).isoformat(),
                    "filename": p.name, "line": ln, "n_pixels": len(cols), "water_level_m": round(z, 2),
                    "surf_from_m": round(float(past[cols].min()), 1) if len(cols) else np.nan,
                    "surf_to_m": round(float(past[cols].max()), 1) if len(cols) else np.nan}
            if len(cols) < args.min_pixels:
                rows.append({**base, "status": "no surf pixels"})
                continue
            try:
                r = period_from_stack(gray[:, a:b], cols=cols)
            except ValueError:
                continue
            rows.append({**base, **r})
    if rows:
        pd.concat([done, pd.DataFrame(rows, columns=FIELDS)], ignore_index=True).to_csv(out, index=False)
    ok = [r for r in rows if r.get("status") == "ok"]
    nos = sum(r.get("status") == "no surf pixels" for r in rows)
    print(f"wave period {args.camera}: {len(rows)} new line-stack(s), {len(ok)} with a clear peak"
          + (f", {nos} with no surf-zone pixels on the line" if nos else "") + f" -> {out}"
          + (f"; latest mean period Tm01 {ok[-1]['tm01_s']:.1f} s" if ok else ""))
    return 0


def compare(args):
    """Camera Tp / Tm01 against the buoys' dominant / average period, hour by hour."""
    import pandas as pd
    out = Path(args.output or HERE / "archive" / f"wave_period_{args.camera}.csv")
    if not out.exists():
        sys.exit(f"no {out} yet")
    c = pd.read_csv(out)
    c = c[c["status"] == "ok"].copy()
    if "n_pixels" not in c or c["n_pixels"].isna().all():
        print("NOTE: these rows predate the surf-zone selection -- rerun with --reprocess first")
    c["hour"] = (c["epoch"] / 3600).round() * 3600
    c = c.groupby("hour")[["tp_s", "tm01_s"]].median()
    print(f"camera periods: {len(c)} hours with a clear peak")
    print(f"  {'buoy':<7} {'hours':>5} {'camera Tp':>10} {'buoy DPD':>9} {'bias':>7} {'RMS':>6} {'r':>6}"
          f"   {'camera Tm01':>11} {'buoy APD':>9} {'r':>6}")
    for st in ("44008", "44013"):
        w = HERE / "archive" / f"waves_{st}.csv"
        if not w.exists():
            continue
        b = pd.read_csv(w)
        b["hour"] = (b["epoch"] / 3600).round() * 3600
        b = b.groupby("hour")[["dpd_s", "apd_s"]].median()
        j = c.join(b, how="inner").dropna()
        if len(j) < 5:
            print(f"  {st:<7} {len(j):>5}  (too few hours in common)")
            continue
        d = j["tp_s"] - j["dpd_s"]
        dm = j["tm01_s"] - j["apd_s"]
        print(f"  {st:<7} {len(j):>5} {j['tp_s'].median():>10.1f} {j['dpd_s'].median():>9.1f} {d.mean():>+7.1f} "
              f"{np.sqrt(np.mean(d ** 2)):>6.1f} {j['tp_s'].corr(j['dpd_s']):>6.2f}"
              f"   {j['tm01_s'].median():>11.1f} {j['apd_s'].median():>9.1f} {j['tm01_s'].corr(j['apd_s']):>6.2f}"
              f"  (Tm01 bias {dm.mean():+.1f} s, RMS {np.sqrt(np.mean(dm ** 2)):.1f} s)")
    print("  Tm01 (mean period, from the whole spectrum) is the camera's reported period: with")
    print("  few surf-zone pixels and 10-min stacks the single spectral peak (Tp) is unstable.")
    print("  The buoys are offshore; the surf zone favours the longer swell, so a modest positive")
    print("  bias is expected.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

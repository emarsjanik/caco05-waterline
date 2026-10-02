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

HOW. On each timestack line, the most variable pixels (the surf zone:
top half by temporal variance, seaward half of the line) are detrended and
their spectra (Welch, 128 s Hann windows, 50% overlap) normalised and
averaged. In the incident band (3-25 s) the peak gives Tp (parabolic
refinement) and the first moment Tm01. The share of variance below 0.04 Hz
(periods over 25 s) is reported as infragravity fraction -- high on a
dissipative surf zone in storms, a useful descriptor in its own right.

OUTPUT archive/wave_period_<cam>.csv, one row per timestack and line:
  epoch (burst middle), time_utc, filename, line, tp_s, tm01_s,
  ig_fraction, peak_ratio, status (ok / weak peak / dark)

Usage:
    python3 timestack_wave_period.py                       new stacks in archive/ras_c2
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
FIELDS = ["epoch", "time_utc", "filename", "line", "tp_s", "tm01_s", "ig_fraction", "peak_ratio", "status"]


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


def period_from_stack(gray, fs=SAMPLE_HZ, band=(1 / 25, 1 / 3), ig_cut=0.04, min_mean=25.0):
    """gray: (time, pixels) along one line -> dict of tp_s, tm01_s, ig_fraction, peak_ratio, status."""
    gray = gray.astype(np.float32)
    if gray.mean() < min_mean:
        return {"status": "dark"}
    std = gray.std(axis=0)
    half = gray.shape[1] // 2
    sea = slice(0, half) if std[:half].mean() >= std[half:].mean() else slice(half, None)
    cols = np.arange(gray.shape[1])[sea]
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
        if den < 0:
            fp = f[k] + 0.5 * (a - c) / den * (f[1] - f[0])
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
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    import pandas as pd
    from runup_from_timestack import split_lines, epoch_from_name, DEFAULT_PIX_DIR
    pix = np.loadtxt(args.pix or f"{DEFAULT_PIX_DIR}/{args.camera}_timestack.pix")[:, :2]
    lines = split_lines(pix)
    out = Path(args.output or HERE / "archive" / f"wave_period_{args.camera}.csv")
    done = pd.read_csv(out) if out.exists() else pd.DataFrame(columns=FIELDS)
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
        for ln in args.lines:
            if not 1 <= ln <= len(lines):
                continue
            a, b = lines[ln - 1]
            try:
                r = period_from_stack(gray[:, a:b])
            except ValueError:
                continue
            mid = e + BURST_MID_OFFSET_S
            rows.append({"epoch": int(mid), "time_utc": datetime.fromtimestamp(mid, tz=timezone.utc).isoformat(),
                         "filename": p.name, "line": ln, **r})
    if rows:
        pd.concat([done, pd.DataFrame(rows, columns=FIELDS)], ignore_index=True).to_csv(out, index=False)
    ok = [r for r in rows if r.get("status") == "ok"]
    print(f"wave period {args.camera}: {len(rows)} new line-stack(s), {len(ok)} with a clear peak -> {out}"
          + (f"; latest Tp {ok[-1]['tp_s']:.1f} s" if ok else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())

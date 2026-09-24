#!/usr/bin/env python3
"""
Does Averaging Nearby Predictions Reduce The Error?
======================================================
Wave height changes slowly -- 0.134 m RMS over a full hour in this
ADCP record -- while the model's per-frame error is around 0.29 m. If
those per-frame errors are largely independent, averaging predictions
from frames close in time should cancel much of the noise while
losing very little real signal.

That is a testable claim, not an assumption, and this tests it on
predictions already made. No GPU, no retraining.

TWO THINGS THIS CHECKS BEYOND THE HEADLINE NUMBER:

  * Whether the errors are actually independent. If consecutive
    predictions err the same way -- because the frames look alike, or
    the whole day is hazy -- averaging cannot help, and the measured
    autocorrelation of the error says which it is up front.

  * What the averaging costs. Combining frames spreads the estimate
    over a window during which the sea state genuinely changes, so
    beyond some width the truth being averaged is no longer the truth
    at any one moment. The script reports the error against the
    centre-frame observation, so that cost is included rather than
    hidden.

Averaging is only legitimate within a run of frames from the same day:
a gap of hours means the sea state has moved on. Windows that would
span a gap larger than --max-gap are not formed.

Usage:
    python3 average_predictions.py --validation owg_c2_H_bright2.validation.csv \\
        --manifest manifest_c2_bright_clean2.csv
"""

import re
import sys
import argparse

import numpy as np
import pandas as pd


def epoch_from_id(i):
    m = re.match(r"^(\d{9,11})\.", str(i))
    return int(m.group(1)) if m else np.nan


def rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a, float) - np.asarray(b, float)) ** 2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--validation", required=True,
                    help="<model>.validation.csv, with id, observed, predicted")
    ap.add_argument("--manifest", default=None,
                    help="Optional; only needed if the validation file lacks timestamps.")
    ap.add_argument("--max-gap", type=float, default=45.0,
                    help="Minutes. Frames further apart than this are not averaged together, "
                         "because the sea state has moved on (default 45).")
    ap.add_argument("--windows", type=int, nargs="+", default=[1, 2, 3, 5, 7])
    ap.add_argument("--output", default=None)
    args = ap.parse_args()

    v = pd.read_csv(args.validation)
    for c in ("id", "observed", "predicted"):
        if c not in v.columns:
            sys.exit(f"{args.validation}: missing column '{c}'")
    v["epoch"] = v["id"].map(epoch_from_id)
    if v["epoch"].isna().any():
        sys.exit("could not read a timestamp from every id")
    v = v.sort_values("epoch").reset_index(drop=True)

    print("=" * 66)
    print("TEMPORAL AVERAGING OF PREDICTIONS")
    print("=" * 66)
    print(f"predictions       : {len(v)}")
    print(f"span              : "
          f"{pd.Timestamp(v.epoch.min(), unit='s')} to {pd.Timestamp(v.epoch.max(), unit='s')}")

    v["err"] = v["predicted"] - v["observed"]
    base = rmse(v["predicted"], v["observed"])
    print(f"single-frame RMSE : {base:.3f}")

    # ---- are consecutive errors independent? ---------------------------
    gap = v["epoch"].diff().fillna(np.inf) / 60.0
    adjacent = gap <= args.max_gap
    if adjacent.sum() >= 20:
        e1 = v["err"].to_numpy()[1:][adjacent.to_numpy()[1:]]
        e0 = v["err"].to_numpy()[:-1][adjacent.to_numpy()[1:]]
        r = float(np.corrcoef(e0, e1)[0, 1])
        print(f"error autocorr    : {r:+.2f} between consecutive frames (n={len(e1)})")
        if r > 0.6:
            print("  Errors are strongly correlated between neighbouring frames, so most of")
            print("  the error is shared rather than independent. Averaging can remove only")
            print("  the independent part, so expect little gain.")
        elif r < 0.3:
            print("  Errors are close to independent, which is the case where averaging helps.")
        else:
            print("  Errors are partly correlated: averaging should help, but less than the")
            print("  square-root-of-n ideal.")
    else:
        print("error autocorr    : too few adjacent frames to measure")

    # ---- how fast does the truth itself move? ---------------------------
    if adjacent.sum() >= 20:
        d = v["observed"].diff()[adjacent]
        print(f"observed change   : {float(np.sqrt((d**2).mean())):.3f} per step between "
              f"adjacent frames")
        print("  Averaging is worth it only while this stays well below the model error.")
    print()

    # ---- rolling averages, centred, never spanning a gap ---------------
    ep = v["epoch"].to_numpy()
    pred = v["predicted"].to_numpy(float)
    obs = v["observed"].to_numpy(float)
    rows = []
    print(f"  {'window':>7}{'n':>7}{'RMSE':>9}{'change':>9}   (centre-frame truth)")
    for w in args.windows:
        half = w // 2
        keep_p, keep_o, used = [], [], 0
        for i in range(len(v)):
            lo, hi = i - half, i + half
            if lo < 0 or hi >= len(v):
                continue
            # Every step inside the window must be close in time, or the
            # window straddles a night or a rejected run and averages
            # across a different sea state.
            steps = np.diff(ep[lo:hi + 1]) / 60.0
            if len(steps) and steps.max() > args.max_gap:
                continue
            keep_p.append(pred[lo:hi + 1].mean())
            keep_o.append(obs[i])          # truth AT THE CENTRE, not averaged
            used += 1
        if used < 20:
            print(f"  {w:>7}{used:>7}      too few complete windows")
            continue
        r = rmse(keep_p, keep_o)
        rows.append((w, used, r))
        print(f"  {w:>7}{used:>7}{r:>9.3f}{100*(r-base)/base:>8.1f}%")

    if rows:
        best = min(rows, key=lambda t: t[2])
        print()
        if best[2] < base - 0.005:
            mins = (best[0] - 1) * 30
            print(f"Best: a {best[0]}-frame average, RMSE {best[2]:.3f} "
                  f"({100*(best[2]-base)/base:+.1f}%).")
            print(f"  That spans about {mins} minutes, so the estimate is no longer")
            print(f"  instantaneous -- a fair trade for a slowly varying quantity, but it")
            print(f"  should be stated alongside the accuracy rather than left implicit.")
        else:
            print("No window improved on single frames by a useful margin.")
            print("  That points to the errors being driven by conditions that persist across")
            print("  neighbouring frames -- haze, light, sea state -- rather than by")
            print("  independent per-frame noise. More averaging will not fix that.")

    if args.output and rows:
        pd.DataFrame(rows, columns=["window", "n", "rmse"]).to_csv(args.output, index=False)
        print(f"\nwrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

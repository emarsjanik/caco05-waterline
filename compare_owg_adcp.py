#!/usr/bin/env python3
"""
Optical Wave Gauge Against The ADCP
======================================
Takes the predictions already made by predict_owg.py and produces the
comparison: scatter, time series, error by wave height, and the effect
of averaging neighbouring frames.

Needs no GPU and no TensorFlow -- it works from the prediction file.

THREE RULES IT ENFORCES, each learned the hard way in this project:

  * Frames the model trained on are excluded from every score. They
    are kept in the time series, drawn differently, because the record
    is more legible unbroken -- but they contribute nothing to the
    accuracy figures.

  * Frames rejected by quality control are excluded by default.
    Including them moves RMSE from 0.288 m to 1.974 m, because a dark
    frame normalised to unit variance becomes noise the model reads as
    heavy breaking. The gauge is defined as the model PLUS its filter;
    scoring it without the filter measures something nobody would
    deploy. --include-rejected shows that contrast deliberately.

  * Averaged results are compared against single-frame results ON THE
    SAME ROWS. A centred window only forms inside an unbroken run of
    frames, which are the clear daytime stretches -- the easy ones. A
    naive comparison against the overall single-frame figure credits
    averaging with that selection and overstates it roughly fourfold.

Usage:
    python3 compare_owg_adcp.py --predictions predictions_all.csv \\
        --clean-manifest manifest_c2_bright_clean2.csv \\
        --out-prefix owg_vs_adcp
"""

import re
import sys
import argparse

import numpy as np
import pandas as pd


def epoch_from_id(i):
    m = re.match(r"^(\d{9,11})\.", str(i))
    return int(m.group(1)) if m else np.nan


def rmse(e):
    return float(np.sqrt(np.mean(np.asarray(e, float) ** 2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--clean-manifest", default=None,
                    help="Manifest of frames that passed quality control. Without it, every "
                         "frame in --predictions is treated as having passed.")
    ap.add_argument("--out-prefix", default="owg_vs_adcp")
    ap.add_argument("--include-rejected", action="store_true")
    ap.add_argument("--window", type=int, default=3,
                    help="Frames to average for the smoothed series (default 3, about an hour)")
    ap.add_argument("--max-gap", type=float, default=45.0)
    ap.add_argument("--unit", default="m")
    args = ap.parse_args()

    p = pd.read_csv(args.predictions)
    for c in ("id", "predicted", "observed"):
        if c not in p.columns:
            sys.exit(f"{args.predictions}: missing column '{c}'")
    p["epoch"] = p["id"].map(epoch_from_id)
    p = p.dropna(subset=["epoch"]).sort_values("epoch").reset_index(drop=True)
    p["time"] = pd.to_datetime(p["epoch"], unit="s")
    if "trained_on" not in p.columns:
        p["trained_on"] = False

    if args.clean_manifest:
        q = pd.read_csv(args.clean_manifest)
        ok = set(q["filename"].str.replace(r"\.jpg$", "", regex=True))
        p["passed_qc"] = p["id"].isin(ok)
    else:
        p["passed_qc"] = True

    print("=" * 70)
    print("OPTICAL WAVE GAUGE vs ADCP")
    print("=" * 70)
    print(f"predictions       : {len(p)}")
    print(f"period            : {p.time.min():%Y-%m-%d} to {p.time.max():%Y-%m-%d}")
    print(f"passed QC         : {int(p.passed_qc.sum())}")
    print(f"used in training  : {int(p.trained_on.sum())}  (excluded from all scores)")

    work = p if args.include_rejected else p[p.passed_qc]
    scored = work[~work.trained_on].copy()
    if len(scored) < 10:
        sys.exit("Too few unseen frames to score.")
    scored["err"] = scored["predicted"] - scored["observed"]

    u = args.unit
    print()
    print(f"UNSEEN FRAMES, n = {len(scored)}")
    print(f"  RMSE            : {rmse(scored.err):.3f} {u}")
    print(f"  bias            : {scored.err.mean():+.3f} {u}")
    print(f"  R2              : {1 - scored.err.var()/scored.observed.var():.2f}")
    print(f"  predicting mean : {scored.observed.std():.3f} {u}")

    if args.clean_manifest and not args.include_rejected:
        rej = p[~p.passed_qc & ~p.trained_on]
        if len(rej) > 10:
            e = rej.predicted - rej.observed
            print(f"  (frames rejected by QC, for contrast: n={len(rej)}, "
                  f"RMSE {rmse(e):.3f} {u}, bias {e.mean():+.3f} {u})")

    print()
    print("  By wave height:")
    edges = [0, 0.5, 1.0, 1.5, 2.0, 10]
    cat = pd.cut(scored["observed"], edges)
    print(f"    {'bin':<14}{'n':>6}{'bias':>9}{'rmse':>9}")
    for b in cat.cat.categories:
        g = scored[cat == b]
        if len(g):
            print(f"    {b.left:.1f}-{b.right:.1f} {u:<8}{len(g):>6}"
                  f"{g.err.mean():>9.3f}{rmse(g.err):>9.3f}")

    # ---- averaging, compared like with like ----------------------------
    w = args.window
    half = w // 2
    ep = work["epoch"].to_numpy()
    pred = work["predicted"].to_numpy(float)
    seen = work["trained_on"].to_numpy()
    idx, avg = [], []
    for i in range(half, len(work) - half):
        if (np.diff(ep[i - half:i + half + 1]) / 60.0).max() > args.max_gap:
            continue
        idx.append(i)
        avg.append(pred[i - half:i + half + 1].mean())
    if idx:
        sub = work.iloc[idx].copy()
        sub["averaged"] = avg
        sub = sub[~sub.trained_on]
        if len(sub) > 10:
            single = rmse(sub.predicted - sub.observed)
            smooth = rmse(sub.averaged - sub.observed)
            print()
            print(f"{w}-frame average (~{(w-1)*30} min), n = {len(sub)}")
            print(f"  single-frame on these same rows : {single:.3f} {u}")
            print(f"  averaged                        : {smooth:.3f} {u}  "
                  f"({100*(smooth-single)/single:+.1f}%)")
            print("  Compared on identical rows: a centred window only forms inside an")
            print("  unbroken run of frames, and those are the clear daytime stretches.")
            sub[["id", "time", "observed", "predicted", "averaged"]].to_csv(
                args.out_prefix + "_averaged.csv", index=False)

    scored.to_csv(args.out_prefix + "_scored.csv", index=False)
    print()
    print(f"wrote {args.out_prefix}_scored.csv")

    # ---- figures --------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable -- no figures written")
        return 0

    fig = plt.figure(figsize=(13, 8.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1], hspace=0.32, wspace=0.25)

    ax = fig.add_subplot(gs[0, 0])
    lo = min(scored.observed.min(), scored.predicted.min())
    hi = max(scored.observed.max(), scored.predicted.max())
    ax.plot([lo, hi], [lo, hi], "k-", lw=1)
    ax.scatter(scored.observed, scored.predicted, s=14, alpha=0.45,
               edgecolor="none", color="#1f77b4")
    ax.set_xlabel(f"ADCP $H_s$ ({u})")
    ax.set_ylabel(f"Optical gauge $H_s$ ({u})")
    ax.set_title(f"RMSE {rmse(scored.err):.3f} {u},  "
                 f"R$^2$ {1 - scored.err.var()/scored.observed.var():.2f},  n = {len(scored)}")
    ax.set_aspect("equal", "box")
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[0, 1])
    ax.axhline(0, color="k", lw=1)
    ax.scatter(scored.observed, scored.err, s=14, alpha=0.45,
               edgecolor="none", color="#d62728")
    ax.set_xlabel(f"ADCP $H_s$ ({u})")
    ax.set_ylabel(f"error ({u})")
    ax.set_title("Error against wave height")
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, :])
    # The ADCP is continuous; the optical record is not, and showing
    # them the same way would imply coverage the camera never had.
    ax.plot(p.time, p.observed, "-", color="0.35", lw=0.9, label="ADCP", zorder=1)
    tr = work[work.trained_on]
    ax.scatter(tr.time, tr.predicted, s=9, color="#aec7e8", edgecolor="none",
               label="optical (training frames)", zorder=2)
    ax.scatter(scored.time, scored.predicted, s=12, color="#1f77b4", edgecolor="none",
               label="optical (unseen)", zorder=3)
    ax.set_ylabel(f"$H_s$ ({u})")
    ax.set_xlabel("date (UTC)")
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.3)
    ax.set_title("Gaps in the optical record are night and frames rejected by quality control")
    fig.autofmt_xdate()

    out = args.out_prefix + ".png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

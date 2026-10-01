#!/usr/bin/env python3
"""
Optical Wave Gauge Validation Scatter (Side By Side)
=======================================================
Observed vs predicted wave height for one or more models, each panel
titled like the USGS OWG figures: RMSE, bias, R2, Willmott's d.

WHICH R2. Two different numbers go by "R2":
  * r2  = squared correlation. Ignores bias and scale errors entirely:
          predictions that are all 0.5 m too high, or half the true
          spread, can still score 1.0.
  * NSE = 1 - SS_res / SS_tot (coefficient of determination, Nash-Sutcliffe).
          Penalises bias and scale. train_marconi_owg.py's "R2" is this one.
Both are printed; the panel title uses --r2 (default r2, the usual
scatter-plot convention). Compare like with like.

Willmott's d = 1 - sum (P-O)^2 / sum (|P-Obar| + |O-Obar|)^2  (0..1).

FAIR COMPARISON. Two models are only comparable on the same frames --
the same days, conditions and wave-height mix. With --common, only ids
present in every file are used (files need an 'id' column; ids are
matched on their leading epoch, so different file-name suffixes still
match). Without it each panel uses its own frames and the n differs.

Input CSVs: id, observed, predicted (train_marconi_owg.py writes these
as <model>.validation.csv).

Usage:
    python3 plot_owg_validation.py owg_c2_H_current_C.validation.csv \\
        chris_owg.csv --labels "CACO05 (ours)" "USGS (Chris)" --common \\
        --output owg_compare.png
"""

import re
import sys
import argparse

import numpy as np
import pandas as pd


def metrics(o, p):
    e = p - o
    ob = o.mean()
    return {"n": len(o),
            "rmse": float(np.sqrt((e ** 2).mean())),
            "bias": float(e.mean()),
            "r2": float(np.corrcoef(o, p)[0, 1] ** 2),
            "nse": float(1 - (e ** 2).sum() / ((o - ob) ** 2).sum()),
            "d": float(1 - (e ** 2).sum() / ((np.abs(p - ob) + np.abs(o - ob)) ** 2).sum()),
            "mae": float(np.abs(e).mean())}


def key(i):
    m = re.match(r"^(\d{9,11})", str(i))
    return m.group(1) if m else str(i)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("csv", nargs="+", help="id,observed,predicted files")
    ap.add_argument("--labels", nargs="+", default=None)
    ap.add_argument("--common", action="store_true",
                    help="Use only frames present in every file")
    ap.add_argument("--r2", choices=["r2", "nse"], default="r2",
                    help="Which R2 the panel title shows (default r2 = squared correlation)")
    ap.add_argument("--output", default="owg_validation.png")
    args = ap.parse_args()

    labels = args.labels or [p.rsplit("/", 1)[-1].replace(".validation.csv", "") for p in args.csv]
    dfs = []
    for p in args.csv:
        d = pd.read_csv(p)
        d.columns = [c.strip().lower() for c in d.columns]
        for c in ("observed", "predicted"):
            if c not in d.columns:
                sys.exit(f"{p}: needs columns observed and predicted (got {list(d.columns)})")
        if "id" in d.columns:
            d["key"] = d["id"].map(key)
        dfs.append(d.dropna(subset=["observed", "predicted"]))
    if args.common:
        if not all("key" in d.columns for d in dfs):
            sys.exit("--common needs an 'id' column in every file")
        common = set.intersection(*(set(d["key"]) for d in dfs))
        dfs = [d[d["key"].isin(common)].drop_duplicates("key") for d in dfs]
        print(f"common frames     : {len(common)}")
        if not common:
            sys.exit("no frames in common -- different periods or id formats")

    print(f"{'model':28s} {'n':>5s} {'RMSE':>7s} {'bias':>7s} {'r2':>6s} {'NSE':>6s} "
          f"{'d':>6s} {'MAE':>6s}")
    res = []
    for lab, d in zip(labels, dfs):
        m = metrics(d["observed"].to_numpy(float), d["predicted"].to_numpy(float))
        res.append(m)
        print(f"{lab:28s} {m['n']:5d} {m['rmse']:7.3f} {m['bias']:+7.3f} {m['r2']:6.3f} "
              f"{m['nse']:6.3f} {m['d']:6.3f} {m['mae']:6.3f}")
    print("   r2 = squared correlation (ignores bias); NSE = 1 - SSres/SStot (penalises it)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(dfs)
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 5), dpi=110, squeeze=False)
    lo = min(min(d["observed"].min(), d["predicted"].min()) for d in dfs)
    hi = max(max(d["observed"].max(), d["predicted"].max()) for d in dfs)
    pad = 0.05 * (hi - lo)
    for k, (ax, lab, d, m) in enumerate(zip(axes[0], labels, dfs, res)):
        ax.scatter(d["observed"], d["predicted"], s=12, alpha=0.55, color="#2a78d6", lw=0)
        ax.plot([lo, hi], [lo, hi], color="#1f1f1e", lw=1.2)
        ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
        ax.set_aspect("equal")
        ax.set_xlabel("Observed H (m)"); ax.set_ylabel("Predicted H (m)")
        r2 = m[args.r2]
        ax.set_title(f"RMSE={m['rmse']:.3f} m, bias={m['bias']:+.3f}, "
                     f"R2={r2:.3f}, d={m['d']:.3f}", fontsize=9.5)
        ax.text(0.03, 0.96, f"{chr(97 + k)}  {lab}  (n={m['n']})", transform=ax.transAxes,
                va="top", fontsize=9, fontweight="bold")
        ax.grid(True, color="#e4e3dc"); ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(args.output)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

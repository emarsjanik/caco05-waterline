#!/usr/bin/env python3
"""
Apply A Trained Optical Wave Gauge And Compare Against The ADCP
=================================================================
Runs a trained model over a set of images, writes the predictions, and
where measured wave data is available, scores them and plots the
comparison.

PREPROCESSING IS READ FROM THE REPORT, NOT RETYPED. Training records
the crop, input size, colour mode and target scaling in
<model>.report.json. Inference must reproduce all of them exactly: a
model trained on a 0.08-0.60 crop of a bright composite, fed a
differently cropped image, returns confident nonsense with no error
and no warning. Passing the report rather than command-line arguments
removes the opportunity to get it wrong.

WHAT COUNTS AS AN HONEST COMPARISON. Predictions on images the model
trained on are not evidence of anything: the model has seen them. The
comparison that matters is against held-out data, which is what
<model>.validation.csv already contains. This script will happily run
over training images -- for producing a time series, say -- but it
warns when it detects them, because a scatter plot that mixes trained
and unseen images looks far better than the model actually is, and
that is how implausible accuracy figures get published.

Usage:
    # score and plot the held-out validation set (the honest comparison)
    python3 predict_owg.py --model owg_models/owg_c2_H_bright2 \\
        --manifest manifest_c2_bright_clean2.csv \\
        --image-dir ~/owg_marconi/images_bright \\
        --only-ids owg_models/owg_c2_H_bright2.validation.csv \\
        --output predictions_val.csv --plot compare_val.png

    # run over every image, for a continuous time series
    python3 predict_owg.py --model ... --manifest ... --image-dir ... \\
        --output predictions_all.csv --plot compare_all.png
"""

import os
import re
import sys
import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import cv2


IMAGENET_MEAN = np.array([0.485, 0.456, 0.406], np.float32)
IMAGENET_STD = np.array([0.229, 0.224, 0.225], np.float32)


def load_image(path, width, height, crop, colour):
    """Identical to the training loader; any divergence invalidates the model."""
    img = cv2.imread(str(path), cv2.IMREAD_COLOR)
    if img is None:
        return None
    if crop:
        h, w = img.shape[:2]
        t, b, l, r = crop
        img = img[int(round(t * h)):int(round(b * h)), int(round(l * w)):int(round(r * w))]
        if img.size == 0:
            return None
    if colour:
        out = cv2.cvtColor(img, cv2.COLOR_BGR2RGB).astype(np.float32)
    else:
        b_, g_, r_ = img[:, :, 0], img[:, :, 1], img[:, :, 2]
        out = (0.21 * r_ + 0.72 * g_ + 0.07 * b_).astype(np.float32)
    return cv2.resize(out, (width, height), interpolation=cv2.INTER_AREA)


def normalise(x, imagenet):
    if imagenet:
        return (x / 255.0 - IMAGENET_MEAN) / IMAGENET_STD
    s = float(x.std())
    return (x - float(x.mean())) / (s if s > 1e-6 else 1.0)


def metrics(obs, pred):
    obs, pred = np.asarray(obs, float), np.asarray(pred, float)
    err = pred - obs
    rmse = float(np.sqrt(np.mean(err ** 2)))
    spread = float(obs.std())
    return {
        "n": int(len(obs)),
        "rmse": rmse,
        "bias": float(err.mean()),
        "r2": float(1 - np.var(err) / np.var(obs)) if np.var(obs) > 0 else float("nan"),
        # The error a model would make by always predicting the mean.
        # Quoting RMSE without it hides whether the model is doing
        # anything at all.
        "rmse_predicting_mean": spread,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True,
                    help="Model stem; reads <stem>.report.json and <stem>.weights.h5")
    ap.add_argument("--manifest", required=True,
                    help="Pairing manifest with filename and, if available, wave_height_m")
    ap.add_argument("--image-dir", required=True)
    ap.add_argument("--only-ids", default=None,
                    help="CSV with an 'id' column; restricts prediction to those images. "
                         "Point this at <stem>.validation.csv for the held-out comparison.")
    ap.add_argument("--output", required=True)
    ap.add_argument("--plot", default=None)
    ap.add_argument("--train-ids", default=None,
                    help="The training CSV. If given, predictions on images the model trained "
                         "on are flagged and excluded from the scores.")
    ap.add_argument("--batch", type=int, default=16)
    args = ap.parse_args()

    stem = os.path.expanduser(args.model)
    rep_path = stem + ".report.json"
    w_path = stem + ".weights.h5"
    for p in (rep_path, w_path):
        if not Path(p).exists():
            sys.exit(f"not found: {p}")
    rep = json.load(open(rep_path))

    width = int(rep.get("img_size", 128))
    height = int(rep.get("img_height") or width)
    crop = tuple(rep["crop"]) if rep.get("crop") else None
    colour = bool(rep.get("colour", False))
    y_mean = float(rep.get("y_mean", 0.0))
    y_std = float(rep.get("y_std", 1.0))
    target = rep.get("target", "H")
    arch = rep.get("arch", "mobilenet")
    aux_cols = rep.get("aux_cols") or []

    print("=" * 70)
    print("OPTICAL WAVE GAUGE -- PREDICTION")
    print("=" * 70)
    print(f"model             : {Path(stem).name}  ({arch}, target {target})")
    print(f"preprocessing     : {width}x{height}, "
          f"{'RGB' if colour else 'grey'}, crop {crop}")
    if rep.get("val_rmse"):
        print(f"reported at train : RMSE {rep['val_rmse']:.3f}, R2 {rep.get('val_r2', float('nan')):.2f}")
    if aux_cols:
        sys.exit("This model uses auxiliary inputs; prediction for those is not implemented.")

    m = pd.read_csv(args.manifest)
    if "filename" not in m.columns:
        sys.exit("manifest needs a 'filename' column")
    m["id"] = m["filename"].str.replace(r"\.jpg$", "", regex=True)

    if args.only_ids:
        keep = set(pd.read_csv(args.only_ids)["id"].astype(str))
        m = m[m["id"].isin(keep)]
        print(f"restricted to     : {len(m)} image(s) from {Path(args.only_ids).name}")

    trained_on = set()
    if args.train_ids and Path(args.train_ids).exists():
        trained_on = set(pd.read_csv(args.train_ids)["id"].astype(str))
    m["trained_on"] = m["id"].isin(trained_on)
    if m["trained_on"].any():
        print(f"  NOTE: {int(m['trained_on'].sum())} image(s) were in the training set and are "
              f"excluded from the scores below.")

    image_dir = os.path.expanduser(args.image_dir)
    print(f"images            : {len(m)}")

    import tensorflow as tf
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from train_marconi_owg import build_model
    model, _ = build_model(height, width, arch, float(rep.get("dropout", 0.5)), 0,
                           bool(rep.get("pretrained", False)))
    model.load_weights(w_path)

    preds, ok_idx = [], []
    buf, buf_idx = [], []

    def flush():
        if not buf:
            return
        X = np.stack(buf)
        if X.ndim == 3:
            X = X[..., None]
        p = model.predict(X, verbose=0).ravel() * y_std + y_mean
        preds.extend(p.tolist())
        ok_idx.extend(buf_idx)
        buf.clear()
        buf_idx.clear()

    for i, (idx, row) in enumerate(m.iterrows(), 1):
        img = load_image(os.path.join(image_dir, row["filename"]), width, height, crop, colour)
        if img is None:
            continue
        buf.append(normalise(img, colour))
        buf_idx.append(idx)
        if len(buf) >= args.batch:
            flush()
        if i % 500 == 0:
            print(f"  {i}/{len(m)}")
    flush()

    out = m.loc[ok_idx].copy()
    out["predicted"] = preds
    col = "wave_height_m" if target == "H" else "wave_period_s"
    has_obs = col in out.columns and out[col].notna().any()
    if has_obs:
        out["observed"] = out[col]
        out["error"] = out["predicted"] - out["observed"]
    out.to_csv(args.output, index=False)
    print(f"wrote             : {args.output}  ({len(out)} prediction(s))")

    if not has_obs:
        print("No measured wave data in the manifest, so no comparison is possible.")
        return 0

    scored = out[~out["trained_on"]]
    if len(scored) == 0:
        print("Every image was in the training set; no honest comparison available.")
        return 0

    unit = "m" if target == "H" else "s"
    s = metrics(scored["observed"], scored["predicted"])
    print()
    print(f"Comparison against the ADCP, {s['n']} image(s) not used in training")
    print(f"  RMSE            : {s['rmse']:.3f} {unit}")
    print(f"  bias            : {s['bias']:+.3f} {unit}")
    print(f"  R2              : {s['r2']:.2f}")
    print(f"  predicting mean : {s['rmse_predicting_mean']:.3f} {unit}  "
          f"(what RMSE would be with no model)")

    print()
    print("  By wave height:")
    edges = [0, 0.5, 1.0, 1.5, 2.0, 10] if target == "H" else [0, 6, 8, 10, 12, 30]
    cat = pd.cut(scored["observed"], edges)
    print(f"    {'bin':<14}{'n':>6}{'bias':>9}{'rmse':>9}")
    for b in cat.cat.categories:
        g = scored[cat == b]
        if len(g) == 0:
            continue
        e = g["predicted"] - g["observed"]
        print(f"    {b.left:.1f}-{b.right:.1f} {unit:<8}{len(g):>6}{e.mean():>9.3f}"
              f"{np.sqrt((e**2).mean()):>9.3f}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(13, 5.2))

        ax = axs[0]
        lo = min(scored["observed"].min(), scored["predicted"].min())
        hi = max(scored["observed"].max(), scored["predicted"].max())
        ax.plot([lo, hi], [lo, hi], "k-", lw=1, zorder=1)
        ax.scatter(scored["observed"], scored["predicted"], s=16, alpha=0.5, zorder=2)
        ax.set_xlabel(f"ADCP {target} ({unit})")
        ax.set_ylabel(f"Optical wave gauge {target} ({unit})")
        ax.set_title(f"RMSE = {s['rmse']:.3f} {unit},  R$^2$ = {s['r2']:.2f},  n = {s['n']}")
        ax.set_aspect("equal", "box")
        ax.grid(alpha=0.3)

        ax = axs[1]
        if "time_utc" in scored.columns:
            t = pd.to_datetime(scored["time_utc"])
            o = scored.sort_values("time_utc")
            t = pd.to_datetime(o["time_utc"])
            ax.plot(t, o["observed"], "k.-", ms=3, lw=0.6, label="ADCP")
            ax.plot(t, o["predicted"], "r.", ms=4, alpha=0.7, label="optical")
            ax.set_xlabel("time (UTC)")
            ax.set_ylabel(f"{target} ({unit})")
            ax.legend()
            ax.grid(alpha=0.3)
            fig.autofmt_xdate()
            # Gaps are expected: nights and rejected frames are absent,
            # so the optical series is intermittent by construction.
            ax.set_title("held-out days only; gaps are night and rejected frames")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print(f"\nwrote             : {args.plot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

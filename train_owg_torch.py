#!/usr/bin/env python3
"""
Optical Wave Gauge, Version 2: Pretrained ConvNeXt, Honest Holdback
=====================================================================
Trains the wave-height network the way Chris Sherwood's
OWG_train_holdback_v3.ipynb does (github.com/csherwood-usgs/owg), with the
changes our own tests called for, and exports a model owg_live.py runs.

WHAT IS DIFFERENT FROM train_marconi_owg.py (Run C)
  * Pretrained backbone. Run C trained a 56M-parameter Inception-ResNetV2
    from scratch on ~1100 images. This starts from ImageNet (timm
    convnext_tiny), or from Chris's Duck NC surf-zone weights (--init).
  * Colour, ImageNet normalisation (what the pretrained stem expects)
    instead of grey with per-image standardisation.
  * No zoom augmentation: zoom changes the apparent size of the waves,
    which is the quantity being measured. Shift is WIDER (+/-15%) and
    rotation +/-5 deg: the camera was re-aimed after training (12% of the
    width, 15% of the height), so the model must not care where exactly
    the surf zone sits in the frame. Trained on the ORIGINAL photos, so
    no redraw into today's view (and none of its pointing errors).
  * Rare wave heights are oversampled by bin (training split only). Run C
    reads waves above 1.5 m 0.4-0.6 m low; that is where the error is.
  * Three-way split by TIME: a contiguous holdback block (default
    2025-02-15 .. 2025-03-01, Chris's choice, so the numbers compare)
    never touched until the end; validation for early stopping is whole
    days from the rest; training is everything else.

OUTPUTS (all with the --output stem)
  .pt                 weights (PyTorch)
  .onnx               the model, de-standardised: it returns metres
  .report.json        input size, crop, input format, scores -- owg_live.py
                      reads it ("input": "rgb_imagenet")
  .check.json/.check_inputs.npy  owg_live.py --check proves the station
                      reproduces these numbers
  .validation.csv, .holdback.csv  id, observed, predicted -- compare models
                      on the same frames with plot_owg_validation.py --common

Run in Colab (GPU):
    !pip install -q timm onnx
    !python train_owg_torch.py --labels all_labels.csv --image-dir images \\
        --output "/content/drive/MyDrive/Optical Waves/owg_c2_H_v2"
The labels CSV is id,H (as the existing manifests: id without .jpg). The
train and val manifests can simply be concatenated: this script re-splits
by time itself, and drops repeated ids.
"""

import os
import sys
import json
import time
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import cv2

MEAN = np.array([0.485, 0.456, 0.406], np.float32)
STD = np.array([0.229, 0.224, 0.225], np.float32)


def load_rgb(path, width, height, crop=None):
    """Same steps as owg_live.py's rgb_imagenet input, minus normalisation."""
    img = cv2.imread(str(path), cv2.IMREAD_COLOR)
    if img is None:
        return None
    if crop:
        h, w = img.shape[:2]
        t, b, l, r = crop
        img = img[int(t * h):int(b * h), int(l * w):int(r * w)]
    img = cv2.resize(img, (width, height), interpolation=cv2.INTER_AREA)
    return cv2.cvtColor(img, cv2.COLOR_BGR2RGB)


def to_tensor_np(rgb):
    x = (rgb.astype(np.float32) / 255.0 - MEAN) / STD
    return np.ascontiguousarray(x.transpose(2, 0, 1))


def augment(rgb, rng, shift=0.15, rot=5.0, shear=0.03, light=0.15):
    h, w = rgb.shape[:2]
    M = cv2.getRotationMatrix2D((w / 2, h / 2), rng.uniform(-rot, rot), 1.0)   # no zoom
    M[0, 1] += rng.uniform(-shear, shear)
    M[0, 2] += rng.uniform(-shift, shift) * w
    M[1, 2] += rng.uniform(-shift, shift) * h
    out = cv2.warpAffine(rgb, M, (w, h), borderMode=cv2.BORDER_REFLECT)
    gain = 1.0 + rng.uniform(-light, light)               # exposure and contrast, mildly
    bias = rng.uniform(-light, light) * 40
    return np.clip(out.astype(np.float32) * gain + bias, 0, 255).astype(np.uint8)


def split_by_time(df, hold_start, hold_end, val_frac, seed):
    t = pd.to_datetime(df["epoch"], unit="s", utc=True)
    hold = (t >= pd.Timestamp(hold_start, tz="UTC")) & (t < pd.Timestamp(hold_end, tz="UTC") + pd.Timedelta(days=1))
    rest = df[~hold]
    days = np.array(sorted(pd.to_datetime(rest["epoch"], unit="s", utc=True).dt.strftime("%Y-%m-%d").unique()))
    rng = np.random.default_rng(seed)
    val_days = set(rng.choice(days, max(1, int(round(val_frac * len(days)))), replace=False))
    is_val = pd.to_datetime(rest["epoch"], unit="s", utc=True).dt.strftime("%Y-%m-%d").isin(val_days)
    return rest[~is_val], rest[is_val], df[hold]


def balance(df, target, n_bins, seed):
    """Oversample each target bin to the same size (training split only)."""
    bins = pd.cut(df[target], n_bins)
    n = max(len(g) for _, g in df.groupby(bins, observed=True))
    parts = [g.sample(n, replace=True, random_state=seed) for _, g in df.groupby(bins, observed=True)]
    return pd.concat(parts, ignore_index=True)


def metrics(obs, pred):
    e = pred - obs
    d_den = np.sum((np.abs(pred - obs.mean()) + np.abs(obs - obs.mean())) ** 2)
    return {"n": int(len(obs)), "rmse": float(np.sqrt(np.mean(e ** 2))), "bias": float(e.mean()),
            "r2": float(1 - np.sum(e ** 2) / np.sum((obs - obs.mean()) ** 2)),
            "d": float(1 - np.sum(e ** 2) / d_den)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--labels", required=True, action="append", help="id,H csv (repeatable)")
    ap.add_argument("--image-dir", required=True)
    ap.add_argument("--ext", default=".jpg")
    ap.add_argument("--target", default="H")
    ap.add_argument("--output", required=True)
    ap.add_argument("--arch", default="convnext_tiny", help="any timm model (default convnext_tiny)")
    ap.add_argument("--init", default=None, help="state dict to start from, e.g. Chris's Duck weights "
                                                 "owg_H_convnext_tiny.pt (default: ImageNet)")
    ap.add_argument("--img-size", type=int, default=320, help="input width (default 320)")
    ap.add_argument("--img-height", type=int, default=256, help="input height (default 256)")
    ap.add_argument("--crop", default=None, help="'top,bottom,left,right' fractions (default none)")
    ap.add_argument("--hold-start", default="2025-02-15")
    ap.add_argument("--hold-end", default="2025-03-01")
    ap.add_argument("--val-frac", type=float, default=0.15, help="fraction of the remaining DAYS")
    ap.add_argument("--bins", type=int, default=10)
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--patience", type=int, default=10)
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--lr", type=float, default=1e-4)
    ap.add_argument("--seed", type=int, default=2018)
    ap.add_argument("--device", default=None)
    ap.add_argument("--max-images", type=int, default=None, help="testing only")
    args = ap.parse_args()

    import torch
    import torch.nn as nn
    import timm

    torch.manual_seed(args.seed)
    dev = args.device or ("cuda" if torch.cuda.is_available() else "cpu")
    crop = tuple(float(v) for v in args.crop.split(",")) if args.crop else None
    W, H = args.img_size, args.img_height

    df = pd.concat([pd.read_csv(p) for p in args.labels], ignore_index=True)
    df.columns = [c.strip() for c in df.columns]
    df = df.drop_duplicates("id").dropna(subset=[args.target]).reset_index(drop=True)
    df["epoch"] = df["id"].astype(str).str.split(".").str[0].astype(int)
    df["path"] = [os.path.join(args.image_dir, i if str(i).lower().endswith(args.ext) else f"{i}{args.ext}")
                  for i in df["id"]]
    have = np.array([os.path.exists(p) for p in df["path"]])
    if not have.any():
        sys.exit(f"none of the {len(df)} images found in {args.image_dir} (first expected: {df['path'].iloc[0]})")
    if not have.all():
        print(f"WARNING: {int((~have).sum())} of {len(df)} labelled images not found -- left out")
    df = df[have].sort_values("epoch").reset_index(drop=True)
    if args.max_images:
        df = df.iloc[np.linspace(0, len(df) - 1, args.max_images).astype(int)].reset_index(drop=True)
    tr, va, ho = split_by_time(df, args.hold_start, args.hold_end, args.val_frac, args.seed)
    print(f"images            : {len(df)} ({pd.to_datetime(df.epoch.min(), unit='s'):%Y-%m-%d} to "
          f"{pd.to_datetime(df.epoch.max(), unit='s'):%Y-%m-%d})")
    print(f"split by time     : train {len(tr)}, validation {len(va)} (whole days), holdback {len(ho)} "
          f"({args.hold_start} .. {args.hold_end}, never used until the end)")
    if not len(tr) or not len(va):
        sys.exit("need training and validation images outside the holdback")

    print("loading images ...", flush=True)
    cache = {p: load_rgb(p, W, H, crop) for p in df["path"]}
    y_mean, y_std = float(tr[args.target].mean()), float(tr[args.target].std())
    tr_bal = balance(tr, args.target, args.bins, args.seed)
    print(f"training rows     : {len(tr_bal)} after balancing {args.bins} wave-height bins "
          f"(target mean {y_mean:.2f}, sd {y_std:.2f} m)")

    class Head(nn.Module):
        """Backbone plus de-standardisation, so the exported model returns metres."""
        def __init__(self):
            super().__init__()
            self.net = timm.create_model(args.arch, pretrained=args.init is None, num_classes=1)
            # --init scratch: random start (testing without internet only)
            self.register_buffer("mean", torch.tensor(y_mean))
            self.register_buffer("std", torch.tensor(y_std))

        def forward(self, x):
            return self.net(x).squeeze(1) * self.std + self.mean

    model = Head()
    if args.init and args.init != "scratch":
        state = torch.load(args.init, map_location="cpu")
        missing, unexpected = model.net.load_state_dict(state, strict=False)
        print(f"initialised from  : {args.init} ({len(missing)} missing, {len(unexpected)} unexpected keys)")
    model.to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=3)
    rng = np.random.default_rng(args.seed)

    def batches(frame, train):
        idx = rng.permutation(len(frame)) if train else np.arange(len(frame))
        for s in range(0, len(idx), args.batch):
            rows = frame.iloc[idx[s:s + args.batch]]
            xs = [to_tensor_np(augment(cache[p], rng) if train else cache[p]) for p in rows["path"]]
            yield torch.from_numpy(np.stack(xs)).to(dev), torch.tensor(rows[args.target].to_numpy(np.float32)).to(dev)

    def predict(frame):
        model.eval(); out = []
        with torch.no_grad():
            for x, _ in batches(frame, False):
                out.append(model(x).cpu().numpy())
        return np.concatenate(out) if out else np.array([])

    best, best_state, bad = np.inf, None, 0
    for ep in range(args.epochs):
        t0 = time.time(); model.train(); tot = 0.0
        for x, y in batches(tr_bal, True):
            opt.zero_grad()
            loss = (((model(x) - y) / y_std) ** 2).mean()     # z-scored MSE
            loss.backward(); opt.step(); tot += loss.item() * len(y)
        pv = predict(va)
        v_rmse = float(np.sqrt(np.mean((pv - va[args.target].to_numpy()) ** 2)))
        sched.step(v_rmse)
        flag = ""
        if v_rmse < best:
            best, bad, flag = v_rmse, 0, "  *"
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        else:
            bad += 1
        print(f"epoch {ep:3d}  train {tot / len(tr_bal):.4f}  validation RMSE {v_rmse:.3f} m"
              f"  ({time.time() - t0:.0f} s){flag}", flush=True)
        if bad >= args.patience:
            print("early stop"); break
    model.load_state_dict(best_state)
    torch.save(model.net.state_dict(), args.output + ".pt")

    out = {}
    for name, frame in (("validation", va), ("holdback", ho)):
        if not len(frame):
            continue
        p = predict(frame)
        pd.DataFrame({"id": frame["id"], "observed": frame[args.target], "predicted": p.round(4)}).to_csv(
            f"{args.output}.{name}.csv", index=False)
        out[name] = metrics(frame[args.target].to_numpy(float), p)

    # export: CPU, fixed input size, returns metres
    model.eval().cpu()
    dummy = torch.zeros(1, 3, H, W)
    torch.onnx.export(model, dummy, args.output + ".onnx", input_names=["image"], output_names=["hs_m"],
                      opset_version=17, dynamo=False)
    sample = (ho if len(ho) else va).iloc[:4]
    X = np.stack([to_tensor_np(cache[p]) for p in sample["path"]]).astype(np.float32)
    with torch.no_grad():
        expected = model(torch.from_numpy(X)).numpy()
    np.save(args.output + ".check_inputs.npy", X)
    Path(args.output + ".check.json").write_text(json.dumps(
        {"inputs": Path(args.output).name + ".check_inputs.npy", "expected": expected.tolist(),
         "framework": "pytorch"}))
    report = {"target": args.target, "arch": args.arch, "framework": "pytorch", "input": "rgb_imagenet",
              "layout": "NCHW", "img_size": W, "img_height": H, "crop": list(crop) if crop else None,
              "init": args.init or "imagenet", "holdback": [args.hold_start, args.hold_end],
              "n_train_unique": int(tr["id"].nunique()), "n_val": int(len(va)),
              "val_rmse": out.get("validation", {}).get("rmse"),
              **{f"{k}_{m}": v for k, d in out.items() for m, v in d.items()}}
    Path(args.output + ".report.json").write_text(json.dumps(report, indent=2))

    print("\n" + "=" * 70)
    for k, m in out.items():
        print(f"{k:<11}: RMSE {m['rmse']:.3f} m  bias {m['bias']:+.3f}  R2 {m['r2']:.2f}  d {m['d']:.2f}  "
              f"(n={m['n']})")
    print("Run C for comparison, same holdback: retrain it with these dates, or compare the\n"
          ".holdback.csv files with plot_owg_validation.py --common.")
    print(f"Copy to the station: {Path(args.output).name}.onnx, .report.json, .check.json, .check_inputs.npy")
    return 0


if __name__ == "__main__":
    sys.exit(main())

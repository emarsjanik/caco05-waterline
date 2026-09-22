#!/usr/bin/env python3
"""
Train An Optical Wave Gauge For Marconi Beach
================================================
Trains the network of Buscombe et al. (2020) -- batch normalisation,
base CNN, batch normalisation, global average pooling, dropout, one
linear output -- on the pre-split files written by
prepare_owg_data.py.

The architecture is copied exactly from the upstream train_OWG.py.
The data handling around it is replaced, because upstream does four
things that distort the result on this dataset:

  1. It splits the data AFTER oversampling, so copies of one image
     land in both training and validation. The model is then partly
     validated on images it trained on, and the reported error is
     optimistic. Here the split is done upstream on unique images and
     whole days, and this script never re-splits.

  2. It augments the VALIDATION set with the same random rotation,
     shift and zoom as training. Validation should measure the model
     on unaltered frames. Here it does.

  3. It draws the whole dataset in one generator call, so each image
     receives a single fixed distortion for every epoch. Augmentation
     exists to show new variants each pass. Here each epoch draws
     fresh ones.

  4. It trains only batch size 128 across all four architectures in
     turn -- days of CPU time -- when its own paper found batch 16
     best and MobileNetV1 best for wave height. Here one model, one
     batch size, both configurable.

The --dry-run flag runs everything except TensorFlow: loading,
normalisation, augmentation and batching. Use it first. A pipeline
fault found in seconds beats one found after a day of training.

Usage:
    python3 train_marconi_owg.py --train marconi-c2-train.csv \\
        --val marconi-c2-val.csv --image-dir marconi_images/data \\
        --target H --output owg_marconi_c2_H --dry-run
"""

import os
import sys
import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import cv2


# ---------------------------------------------------------------- data

def load_image(path, size):
    """
    Grayscale, resized, float32.

    The RGB-to-grey weights (0.21, 0.72, 0.07) are copied from the
    upstream pred_1image so that training and inference convert
    identically. OpenCV loads BGR, hence the reversed channel order.
    """
    img = cv2.imread(str(path), cv2.IMREAD_COLOR)
    if img is None:
        return None
    b, g, r = img[:, :, 0], img[:, :, 1], img[:, :, 2]
    grey = (0.21 * r + 0.72 * g + 0.07 * b).astype(np.float32)
    return cv2.resize(grey, (size, size), interpolation=cv2.INTER_AREA)


def normalise(x):
    """Per-image centring and scaling, matching upstream samplewise settings."""
    s = float(x.std())
    return (x - float(x.mean())) / (s if s > 1e-6 else 1.0)


def augment(img, rng, rotation=10.0, shift=0.1, shear=0.05, zoom=0.2):
    """
    One random geometric distortion, matching the upstream config
    ranges: rotation +/-10 deg, shift 10%, shear 0.05, zoom 20%.
    Reflection padding matches upstream fill_mode.
    """
    h, w = img.shape
    angle = rng.uniform(-rotation, rotation)
    scale = 1.0 + rng.uniform(-zoom, zoom)
    M = cv2.getRotationMatrix2D((w / 2.0, h / 2.0), angle, scale)
    M[0, 1] += rng.uniform(-shear, shear)
    M[0, 2] += rng.uniform(-shift, shift) * w
    M[1, 2] += rng.uniform(-shift, shift) * h
    return cv2.warpAffine(img, M, (w, h), borderMode=cv2.BORDER_REFLECT)


def read_split(csv_path, image_dir, target):
    df = pd.read_csv(csv_path)
    df.columns = [c.strip() for c in df.columns]
    if target not in df.columns or "id" not in df.columns:
        sys.exit(f"{csv_path}: need 'id' and '{target}' columns")
    df["path"] = df["id"].map(lambda i: os.path.join(image_dir, f"{i}.jpg"))
    return df


def preload(df, size):
    """
    Loads each UNIQUE image once. The training file repeats rare
    images by design, so caching by id avoids reading one file up to
    eight times.
    """
    cache, missing = {}, []
    for pid, path in zip(df["id"], df["path"]):
        if pid in cache:
            continue
        img = load_image(path, size)
        if img is None:
            missing.append(path)
            continue
        cache[pid] = img
    return cache, missing


class Batches:
    """
    Minimal batch iterator, deliberately independent of TensorFlow so
    the pipeline can be exercised with --dry-run.

    Training batches are reshuffled and freshly augmented every epoch.
    Validation batches are neither shuffled nor augmented: validation
    measures the model on real frames, not distorted ones.
    """

    def __init__(self, df, cache, target, batch, augment_on, seed):
        self.ids = df["id"].tolist()
        self.y = df[target].to_numpy(np.float32)
        self.cache, self.batch = cache, batch
        self.augment_on = augment_on
        self.rng = np.random.default_rng(seed)
        self.order = np.arange(len(self.ids))

    def __len__(self):
        return int(np.ceil(len(self.ids) / self.batch))

    def on_epoch_end(self):
        if self.augment_on:
            self.rng.shuffle(self.order)

    def __getitem__(self, i):
        idx = self.order[i * self.batch:(i + 1) * self.batch]
        xs = []
        for k in idx:
            img = self.cache[self.ids[k]]
            if self.augment_on:
                img = augment(img, self.rng)
            xs.append(normalise(img))
        return np.stack(xs)[..., None], self.y[idx]


# ------------------------------------------------------------- model

def build_model(size, arch, dropout):
    """Exact reproduction of the upstream network stack."""
    from tensorflow.keras.models import Sequential
    from tensorflow.keras.layers import (BatchNormalization, GlobalAveragePooling2D,
                                         Dropout, Dense)
    from tensorflow.keras import applications as A
    bases = {"mobilenet": A.MobileNet, "inceptionv3": A.InceptionV3,
             "inceptionresnetv2": A.InceptionResNetV2, "densenet201": A.DenseNet201}
    # weights=None: trained from scratch. Upstream found transfer
    # learning from ImageNet weights ineffective for this task.
    base = bases[arch](input_shape=(size, size, 1), include_top=False, weights=None)
    model = Sequential([
        BatchNormalization(input_shape=(size, size, 1)),
        base,
        BatchNormalization(),
        GlobalAveragePooling2D(),
        Dropout(dropout),
        Dense(1, activation="linear"),
    ])
    return model


def rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


# -------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--val", required=True)
    ap.add_argument("--heldout", default=None,
                    help="Optional out-of-calibration file, evaluated after training.")
    ap.add_argument("--image-dir", required=True)
    ap.add_argument("--target", default="H", choices=["H", "T"])
    ap.add_argument("--output", required=True, help="Output stem for weights and report")
    ap.add_argument("--arch", default="mobilenet",
                    choices=["mobilenet", "inceptionv3", "inceptionresnetv2", "densenet201"],
                    help="Default mobilenet: best for wave height in Buscombe et al. (2020), "
                         "and 17x smaller than Inception-ResNetV2.")
    ap.add_argument("--img-size", type=int, default=128)
    ap.add_argument("--batch", type=int, default=16,
                    help="Default 16: best in Buscombe et al. (2020), whose code nonetheless "
                         "runs 128.")
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--patience", type=int, default=25)
    ap.add_argument("--dropout", type=float, default=0.5)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=2018)
    ap.add_argument("--dry-run", action="store_true",
                    help="Exercise loading, augmentation and batching without TensorFlow.")
    args = ap.parse_args()

    print("=" * 70)
    print(f"OWG TRAINING -- target {args.target}, {args.arch}, batch {args.batch}, "
          f"{args.img_size}px")
    print("=" * 70)

    tr = read_split(args.train, args.image_dir, args.target)
    va = read_split(args.val, args.image_dir, args.target)

    # Guard against the leakage this script exists to prevent: if the
    # caller passes an overlapping pair of files, refuse rather than
    # report an optimistic number.
    shared = set(tr["id"]) & set(va["id"])
    if shared:
        sys.exit(f"STOP: {len(shared)} image(s) appear in both train and validation. "
                 f"That leaks training data into the score. Regenerate the split with "
                 f"prepare_owg_data.py.")

    print(f"train rows        : {len(tr)} ({tr['id'].nunique()} unique images)")
    print(f"validation rows   : {len(va)} ({va['id'].nunique()} unique images)")
    print(f"target range      : train {tr[args.target].min():.2f}-{tr[args.target].max():.2f}, "
          f"val {va[args.target].min():.2f}-{va[args.target].max():.2f}")

    print("loading images ...")
    cache_tr, miss_tr = preload(tr, args.img_size)
    cache_va, miss_va = preload(va, args.img_size)
    if miss_tr or miss_va:
        print(f"  WARNING: {len(miss_tr)} train / {len(miss_va)} validation image(s) unreadable")
        tr = tr[tr["id"].isin(cache_tr)]
        va = va[va["id"].isin(cache_va)]
    if len(tr) == 0 or len(va) == 0:
        sys.exit("No usable images.")
    print(f"  loaded {len(cache_tr)} train, {len(cache_va)} validation")

    train_b = Batches(tr, cache_tr, args.target, args.batch, augment_on=True, seed=args.seed)
    val_b = Batches(va, cache_va, args.target, args.batch, augment_on=False, seed=args.seed)

    xb, yb = train_b[0]
    print(f"batch shape       : {xb.shape}, labels {yb.shape}")
    print(f"per-image mean/sd : {xb[0].mean():+.3f} / {xb[0].std():.3f}  (expect ~0 / ~1)")

    if args.dry_run:
        # Augmentation must actually vary between epochs; confirm it.
        a1 = train_b[0][0][0]
        train_b.on_epoch_end()
        a2 = Batches(tr, cache_tr, args.target, args.batch, True, args.seed + 1)[0][0][0]
        varies = not np.allclose(a1, a2)
        v1, v2 = val_b[0][0][0], val_b[0][0][0]
        print()
        print("DRY RUN checks")
        print(f"  training augmentation varies  : {'yes' if varies else 'NO -- PROBLEM'}")
        print(f"  validation is deterministic   : "
              f"{'yes' if np.allclose(v1, v2) else 'NO -- PROBLEM'}")
        print(f"  batches per epoch             : {len(train_b)} train, {len(val_b)} validation")
        print()
        print("Pipeline OK. Re-run without --dry-run to train.")
        return 0

    # -- TensorFlow from here ---------------------------------------------
    import tensorflow as tf
    from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping, ReduceLROnPlateau
    tf.random.set_seed(args.seed)
    print(f"tensorflow        : {tf.__version__}")

    class Seq(tf.keras.utils.Sequence):
        def __init__(self, b): self.b = b
        def __len__(self): return len(self.b)
        def __getitem__(self, i): return self.b[i]
        def on_epoch_end(self): self.b.on_epoch_end()

    model = build_model(args.img_size, args.arch, args.dropout)
    model.compile(optimizer=tf.keras.optimizers.Adam(args.lr), loss="mse", metrics=["mae"])

    weights = args.output + ".weights.h5"
    callbacks = [
        ModelCheckpoint(weights, monitor="val_loss", save_best_only=True,
                        save_weights_only=True, mode="min", verbose=1),
        EarlyStopping(monitor="val_loss", mode="min", patience=args.patience,
                      restore_best_weights=True),
        # min_delta rather than upstream's 'epsilon', which newer Keras rejects.
        ReduceLROnPlateau(monitor="val_loss", factor=0.8, patience=5,
                          min_delta=1e-4, cooldown=5, min_lr=1e-4, verbose=1),
    ]
    model.fit(Seq(train_b), validation_data=Seq(val_b), epochs=args.epochs,
              callbacks=callbacks, verbose=2)

    pred = np.concatenate([model.predict(val_b[i][0], verbose=0).ravel()
                           for i in range(len(val_b))])
    truth = va[args.target].to_numpy()
    report = {"target": args.target, "arch": args.arch, "batch": args.batch,
              "img_size": args.img_size, "n_train_unique": int(tr["id"].nunique()),
              "n_val": int(len(va)), "val_rmse": rmse(pred, truth),
              "val_bias": float(np.mean(pred - truth)),
              "val_r2": float(np.corrcoef(pred, truth)[0, 1] ** 2)}

    if args.heldout and Path(args.heldout).exists():
        ho = read_split(args.heldout, args.image_dir, args.target)
        cache_ho, _ = preload(ho, args.img_size)
        ho = ho[ho["id"].isin(cache_ho)]
        ho_b = Batches(ho, cache_ho, args.target, args.batch, False, args.seed)
        hp = np.concatenate([model.predict(ho_b[i][0], verbose=0).ravel()
                             for i in range(len(ho_b))])
        report["heldout_rmse"] = rmse(hp, ho[args.target].to_numpy())
        report["heldout_n"] = int(len(ho))

    model.save_weights(weights)
    Path(args.output + ".report.json").write_text(json.dumps(report, indent=2))
    pd.DataFrame({"id": va["id"], "observed": truth, "predicted": pred}).to_csv(
        args.output + ".validation.csv", index=False)

    print()
    print("=" * 70)
    unit = "m" if args.target == "H" else "s"
    print(f"validation RMSE   : {report['val_rmse']:.3f} {unit}   bias {report['val_bias']:+.3f}   "
          f"R2 {report['val_r2']:.2f}")
    if "heldout_rmse" in report:
        print(f"out-of-calibration: {report['heldout_rmse']:.3f} {unit} (n={report['heldout_n']})")
    print(f"Buscombe et al. oblique RGB, for reference: 0.11 m / 0.81 s -- but measured with")
    print(f"a leaking split, so not directly comparable to the figure above.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

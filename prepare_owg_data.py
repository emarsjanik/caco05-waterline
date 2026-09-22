#!/usr/bin/env python3
"""
Prepare Marconi Data For The Optical Wave Gauge
==================================================
Converts manifest_c2.csv (from pair_images_waves.py) into the format
the Buscombe et al. OWG repository expects, fetches the images, and
builds the train/validation split in a way that avoids two problems
the upstream code has on a dataset shaped like this one.

WHAT IT FIXES, and why each matters here specifically:

  1. OVERSAMPLING. Upstream draws len(df)/2 samples, with replacement,
     from EACH of ten equal-width height bins. On this dataset the top
     bin holds 3 images, so each would appear 383 times -- 10% of the
     whole training set from three frames. The network would memorise
     them. Here each bin is capped at --max-repeat copies of any image,
     so rare conditions are up-weighted without being reduced to rote.

  2. LEAKAGE. If the train/validation split is taken after
     oversampling, duplicates of one image land on both sides and the
     model is partly validated on its own training data. Here the
     split is made on UNIQUE images first, and oversampling is applied
     to the training side only. Validation is never resampled.

  3. TEMPORAL LEAKAGE. Consecutive half-hourly timex frames of the
     same sea state are nearly identical. A random split puts
     neighbours on both sides and flatters the score in the same way.
     Here --split-by-day keeps each day entirely on one side, so
     validation measures performance on genuinely unseen conditions.

  4. HELD-OUT EXTREMES. Upstream withholds the top and bottom 5% for
     out-of-calibration testing. On this dataset the top 5% is most of
     the storm data. Withholding is kept as an option but defaults to
     a smaller fraction, and the script reports what it removes.

Output is the upstream CSV format (id, H, T), so train_OWG.py runs
unchanged apart from the filename handling patched separately.

Usage:
    python3 prepare_owg_data.py --manifest manifest_c2.csv \\
        --s3-prefix s3://cmgp-coastcam/cameras/caco-04/products/ \\
        --image-dir ~/OpticalWaveGauging_DNN/marconi_images/data \\
        --output-stem ~/OpticalWaveGauging_DNN/marconi-c2
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--s3-prefix", default=None,
                    help="If given, images missing from --image-dir are fetched from here.")
    ap.add_argument("--image-dir", required=True)
    ap.add_argument("--output-stem", required=True,
                    help="Writes <stem>.csv (all), <stem>-train.csv and <stem>-val.csv")
    ap.add_argument("--val-fraction", type=float, default=0.25)
    ap.add_argument("--split-by-day", action="store_true", default=True,
                    help="Keep each day wholly in train or validation (default on). "
                         "Consecutive frames of one sea state are near-duplicates, so a "
                         "random split leaks information and overstates accuracy.")
    ap.add_argument("--no-split-by-day", dest="split_by_day", action="store_false")
    ap.add_argument("--withhold-upper", type=float, default=0.0,
                    help="Percent of highest wave heights withheld for out-of-calibration "
                         "testing (default 0). Upstream uses 5, which on this dataset "
                         "removes most of the storm data from training.")
    ap.add_argument("--withhold-lower", type=float, default=0.0)
    ap.add_argument("--bins", type=int, default=10)
    ap.add_argument("--max-repeat", type=int, default=8,
                    help="Cap on how many times any single image may be repeated when "
                         "balancing bins (default 8). Upstream is effectively unbounded -- "
                         "383x for the rarest images on this dataset.")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--no-download", action="store_true")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    m = pd.read_csv(args.manifest)
    need = {"filename", "time_utc", "wave_height_m", "wave_period_s"}
    if not need.issubset(m.columns):
        sys.exit(f"manifest lacks columns: {need - set(m.columns)}")

    print("=" * 70)
    print("PREPARE OWG DATA")
    print("=" * 70)
    print(f"manifest pairs    : {len(m)}")

    # ---- fetch images --------------------------------------------------
    img_dir = Path(os.path.expanduser(args.image_dir))
    img_dir.mkdir(parents=True, exist_ok=True)
    missing = [f for f in m["filename"] if not (img_dir / f).exists()]
    if missing and not args.no_download:
        if not args.s3_prefix:
            sys.exit(f"{len(missing)} image(s) absent from {img_dir} and no --s3-prefix given")
        print(f"fetching          : {len(missing)} image(s) from S3 ...")
        failed = 0
        for i, f in enumerate(missing, 1):
            r = subprocess.run(["aws", "s3", "cp", args.s3_prefix.rstrip("/") + "/" + f,
                                str(img_dir / f), "--only-show-errors"],
                               capture_output=True)
            if r.returncode != 0:
                failed += 1
            if i % 250 == 0:
                print(f"    {i}/{len(missing)}")
        if failed:
            print(f"  WARNING: {failed} download(s) failed -- those pairs are dropped")
    present = m["filename"].map(lambda f: (img_dir / f).exists())
    m = m[present].copy()
    print(f"images on disk    : {len(m)}")
    if len(m) == 0:
        sys.exit("No images available.")

    # ---- upstream id format --------------------------------------------
    # Upstream builds path = image_dir/id + ".jpg", so id is the filename
    # without its extension.
    m["id"] = m["filename"].str.replace(r"\.jpg$", "", regex=True)
    m["H"] = m["wave_height_m"]
    m["T"] = m["wave_period_s"]
    m["day"] = m["time_utc"].str[:10]

    # ---- out-of-calibration hold-out -----------------------------------
    held = pd.DataFrame()
    if args.withhold_upper > 0 or args.withhold_lower > 0:
        hi = np.percentile(m["H"], 100 - args.withhold_upper) if args.withhold_upper else np.inf
        lo = np.percentile(m["H"], args.withhold_lower) if args.withhold_lower else -np.inf
        mask = (m["H"] > hi) | (m["H"] < lo)
        held = m[mask].copy()
        m = m[~mask].copy()
        print(f"withheld          : {len(held)} image(s) outside "
              f"{lo:.2f}-{hi:.2f} m for out-of-calibration testing")
        if args.withhold_upper and hi < 2.5:
            print(f"  WARNING: the upper cut at {hi:.2f} m removes storm conditions from "
                  f"training entirely. The model will not have seen them.")

    # ---- split on UNIQUE images, before any oversampling ---------------
    if args.split_by_day:
        days = np.array(sorted(m["day"].unique()))
        rng.shuffle(days)
        n_val_days = max(1, int(round(len(days) * args.val_fraction)))
        val_days = set(days[:n_val_days])
        is_val = m["day"].isin(val_days)
        print(f"split             : by day -- {len(days) - n_val_days} train days, "
              f"{n_val_days} validation days")
    else:
        is_val = pd.Series(rng.random(len(m)) < args.val_fraction, index=m.index)
        print("split             : random by image (WARNING: neighbouring frames leak)")
    train, val = m[~is_val].copy(), m[is_val].copy()
    print(f"                    {len(train)} train, {len(val)} validation images")

    # ---- capped rebalancing, TRAINING side only ------------------------
    edges = np.linspace(train["H"].min(), train["H"].max(), args.bins + 1)
    train["bin"] = np.clip(np.digitize(train["H"], edges) - 1, 0, args.bins - 1)
    counts = train["bin"].value_counts().reindex(range(args.bins), fill_value=0)
    target = int(counts[counts > 0].median())
    parts = []
    print()
    print(f"rebalancing to ~{target} per bin, at most {args.max_repeat}x any image:")
    print(f"  {'bin (m)':<14}{'unique':>8}{'drawn':>8}{'max rep':>9}")
    for b in range(args.bins):
        grp = train[train["bin"] == b]
        if len(grp) == 0:
            continue
        want = min(target, len(grp) * args.max_repeat)
        # Deterministic tiling, not sampling with replacement. Random
        # draws with replacement bound the TOTAL count but not any one
        # image's count -- an earlier version promised "at most 8x" and
        # delivered 15x through chance clustering. Tiling every image
        # floor(want/n) times and topping up without replacement makes
        # the ceiling exact: no image exceeds ceil(want/n) <= max_repeat.
        full, extra = divmod(want, len(grp))
        tiled = [grp] * full
        if extra:
            tiled.append(grp.sample(n=extra, replace=False,
                                    random_state=int(rng.integers(1 << 31))))
        take = pd.concat(tiled) if tiled else grp.iloc[0:0]
        parts.append(take)
        rep = take["id"].value_counts().max()
        print(f"  {edges[b]:4.2f}-{edges[b+1]:4.2f} m  {len(grp):>8}{want:>8}{rep:>8}x")
    train_bal = pd.concat(parts).sample(frac=1.0, random_state=args.seed)

    # ---- write ----------------------------------------------------------
    stem = Path(os.path.expanduser(args.output_stem))
    cols = ["id", "H", "T"]
    m[cols].to_csv(str(stem) + ".csv", index=False)
    train_bal[cols].to_csv(str(stem) + "-train.csv", index=False)
    val[cols].to_csv(str(stem) + "-val.csv", index=False)
    if len(held):
        held[cols].to_csv(str(stem) + "-heldout.csv", index=False)

    print()
    print(f"wrote {stem}.csv            ({len(m)} unique images)")
    print(f"wrote {stem}-train.csv      ({len(train_bal)} rows after rebalancing)")
    print(f"wrote {stem}-val.csv        ({len(val)} rows, NEVER resampled)")
    if len(held):
        print(f"wrote {stem}-heldout.csv    ({len(held)} out-of-calibration rows)")

    # Verified, not asserted: both checks inspect the written data.
    overlap = set(train_bal["id"]) & set(val["id"])
    print()
    print(f"leakage check     : {len(overlap)} image(s) in both train and validation"
          f"{'  -- OK' if not overlap else '  -- PROBLEM'}")
    if args.split_by_day:
        shared = set(train["day"]) & set(val["day"])
        print(f"                    {len(shared)} day(s) in both train and validation"
              f"{'  -- OK' if not shared else '  -- PROBLEM'}")
    worst = train_bal["id"].value_counts().max()
    print(f"                    most-repeated training image: {worst}x "
          f"(cap {args.max_repeat}){'  -- OK' if worst <= args.max_repeat else '  -- PROBLEM'}")


if __name__ == "__main__":
    main()

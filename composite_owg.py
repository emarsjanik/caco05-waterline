#!/usr/bin/env python3
"""
Multi-Product Composite Images For The Wave Gauge
===================================================
Each Argus collection writes several images of the same 10 minutes, and
each carries different wave information:

  snap    one instant: crests, their spacing and refraction
  timex   the 10-minute mean: where foam persists (the surf zone's extent)
  bright  the brightest value each pixel reached: every breaker that passed
  dark    the darkest value: water between breakers
  var     the variance, where the station writes it: how much each pixel changes

The current models see ONE product (bright). This stacks several into the
three channels a pretrained network takes, so the same trainer
(train_owg_torch.py) and the same holdback can test whether the
combination measures waves better. Brightness is NOT normalised per
image: how bright the foam is is part of the signal.

Default recipe (--recipe): R = snap, G = bright - dark (the range each
pixel swept: breaking intensity), B = timex. With var available,
"snap,var,timex" is the other natural choice. Any product name or a
difference "a-b" can be a channel.

STEPS
  1. list    what products S3 holds for one collection time (names vary by station)
  2. fetch   download the products for every labelled frame of a manifest
  3. build   write <id>.composite.jpg images + a labels CSV (id,H) for train_owg_torch.py
  4. preview one frame: each product and the composite side by side

Usage (on the station, where aws has the S3 credentials):
    python3 composite_owg.py list
    python3 composite_owg.py fetch --products bright dark snap timex
    python3 composite_owg.py build
    python3 composite_owg.py preview
Frames from before 2025-01 are named CACO03; if fetch reports them missing,
run it again with --s3 s3://cmgp-coastcam/cameras/caco-03/products/.

Then copy composite/ (images/, bright/ and the two label files) to Google
Drive and train twice with identical settings -- the composite, and the
bright product alone on exactly the same frames -- so the difference is
the inputs and nothing else:
    !python train_owg_torch.py --labels composite/labels_composite.csv --image-dir composite/images \\
        --output owg_c2_H_composite --img-size 384 --img-height 320 --lr 3e-5 --patience 15
    !python train_owg_torch.py --labels composite/labels_bright_same_frames.csv --image-dir composite/bright \\
        --output owg_c2_H_bright_ref --img-size 384 --img-height 320 --lr 3e-5 --patience 15
"""

import sys
import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import cv2

HERE = Path(__file__).resolve().parent
S3 = "s3://cmgp-coastcam/cameras/caco-04/products/"
OUT = HERE / "composite"


def product_name(filename, product):
    """'<epoch>....c2.bright.jpg' -> the same collection's '<product>' file name."""
    stem = filename[:-4] if filename.endswith(".jpg") else filename
    return stem.rsplit(".", 1)[0] + f".{product}.jpg"


def load_manifest(path, every):
    m = pd.read_csv(path)
    col = "filename" if "filename" in m else "id"
    hcol = "wave_height_m" if "wave_height_m" in m else "H"
    m = m.rename(columns={col: "filename", hcol: "H"})
    m["filename"] = [f if str(f).endswith(".jpg") else f"{f}.jpg" for f in m["filename"]]
    m = m.dropna(subset=["H"]).drop_duplicates("filename").reset_index(drop=True)
    return m.iloc[::every].reset_index(drop=True)


def cmd_list(args):
    """Every object S3 holds for one collection time (prefix = the epoch)."""
    m = load_manifest(args.manifest, 1)
    epoch = args.epoch or m["filename"].iloc[len(m) // 2].split(".")[0]
    r = subprocess.run(["aws", "s3", "ls", args.s3 + str(epoch)], capture_output=True, text=True)
    print(r.stdout or r.stderr or f"nothing under {args.s3}{epoch}")
    print("product names are the word before .jpg; pass the ones you want to fetch --products")
    return 0


def cmd_fetch(args):
    m = load_manifest(args.manifest, args.every)
    todo = []
    for prod in args.products:
        d = args.out / prod
        d.mkdir(parents=True, exist_ok=True)
        todo += [(product_name(f, prod), d) for f in m["filename"] if not (d / product_name(f, prod)).exists()]
    print(f"{len(m)} labelled frames x {len(args.products)} products: {len(todo)} to download "
          f"(~{len(todo) * 1.3 / 1024:.1f} GB)")
    if args.dry_run:
        for name, d in todo[:5]:
            print(f"  aws s3 cp {args.s3}{name} {d}/")
        return 0
    failed = 0
    for i, (name, d) in enumerate(todo, 1):
        r = subprocess.run(["aws", "s3", "cp", "--only-show-errors", args.s3 + name, str(d / name)],
                           capture_output=True, text=True)
        if r.returncode:
            failed += 1
            if failed <= 5:
                print(f"  not fetched: {name} ({(r.stderr or '').strip()[:120]})")
        if i % 200 == 0:
            print(f"  {i}/{len(todo)}", flush=True)
    print(f"done: {len(todo) - failed} fetched, {failed} missing on S3")
    return 0


def channel(spec, imgs):
    """'snap' -> that product in grey; 'bright-dark' -> the difference, clipped to 0..255."""
    if "-" in spec:
        a, b = spec.split("-")
        return np.clip(imgs[a].astype(np.int16) - imgs[b].astype(np.int16), 0, 255).astype(np.uint8)
    return imgs[spec]


def make_composite(filename, recipe, src, width, height):
    """-> (H x W x 3 uint8 in RGB order, None) or (None, missing product)."""
    need = sorted({p for c in recipe for p in c.split("-")})
    imgs = {}
    for p in need:
        g = cv2.imread(str(src / p / product_name(filename, p)), cv2.IMREAD_GRAYSCALE)
        if g is None:
            return None, p
        imgs[p] = cv2.resize(g, (width, height), interpolation=cv2.INTER_AREA)
    return np.dstack([channel(c, imgs) for c in recipe]), None


def cmd_build(args):
    recipe = args.recipe.split(",")
    if len(recipe) != 3:
        sys.exit("--recipe needs three channels (the network takes RGB), e.g. snap,bright-dark,timex")
    m = load_manifest(args.manifest, args.every)
    img_dir = args.out / "images"
    img_dir.mkdir(parents=True, exist_ok=True)
    rows, missing = [], {}
    for f, h in zip(m["filename"], m["H"]):
        cid = f[:-4].rsplit(".", 1)[0] + ".composite"
        dst = img_dir / f"{cid}.jpg"
        if not dst.exists():
            rgb, miss = make_composite(f, recipe, args.out, args.width, args.height)
            if rgb is None:
                missing[miss] = missing.get(miss, 0) + 1
                continue
            cv2.imwrite(str(dst), cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, 95])
        rows.append((cid, h))
    lab = args.out / "labels_composite.csv"
    pd.DataFrame(rows, columns=["id", "H"]).to_csv(lab, index=False)
    # the same frames with the bright product alone: the like-for-like comparison
    base = args.out / "labels_bright_same_frames.csv"
    pd.DataFrame([(c.replace(".composite", ".bright"), h) for c, h in rows], columns=["id", "H"]).to_csv(base, index=False)
    (args.out / "recipe.txt").write_text(f"R,G,B = {recipe}  size {args.width}x{args.height}\n")
    print(f"composites        : {len(rows)} in {img_dir} (R,G,B = {', '.join(recipe)})")
    if missing:
        print("skipped, product missing: " + ", ".join(f"{k} {v}" for k, v in missing.items()))
    print(f"labels            : {lab}\n                    {base} (bright, same frames, for the comparison)")
    return 0


def cmd_preview(args):
    recipe = args.recipe.split(",")
    m = load_manifest(args.manifest, 1)
    first = recipe[0].split("-")[0]
    m = m[[(args.out / first / product_name(f, first)).exists() for f in m["filename"]]]
    m = m.iloc[(m["H"] - args.hs).abs().argsort()]          # a fetched frame near --hs metres
    for f, h in zip(m["filename"], m["H"]):
        rgb, miss = make_composite(f, recipe, args.out, 612, 512)
        if rgb is not None:
            break
    else:
        sys.exit("no frame has all the products yet -- run fetch first")
    tiles = []
    for p in sorted({q for c in recipe for q in c.split("-")}):
        g = cv2.resize(cv2.imread(str(args.out / p / product_name(f, p)), cv2.IMREAD_GRAYSCALE), (612, 512))
        tiles.append((p, cv2.cvtColor(g, cv2.COLOR_GRAY2BGR)))
    for i, c in enumerate(recipe):
        tiles.append((f"channel {'RGB'[i]} = {c}", cv2.cvtColor(rgb[:, :, i], cv2.COLOR_GRAY2BGR)))
    tiles.append(("composite", cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR)))
    while len(tiles) % 4:
        tiles.append(("", np.zeros((512, 612, 3), np.uint8)))
    for t, im in tiles:
        cv2.rectangle(im, (0, 0), (612, 30), (0, 0, 0), -1)
        cv2.putText(im, t, (8, 22), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255, 255, 255), 2, cv2.LINE_AA)
    rows = [np.hstack([im for _, im in tiles[i:i + 4]]) for i in range(0, len(tiles), 4)]
    out = args.out / "composite_preview.jpg"
    cv2.imwrite(str(out), np.vstack(rows), [cv2.IMWRITE_JPEG_QUALITY, 88])
    print(f"preview           : {out}  ({f.split('.')[0]}, Hs {h:.2f} m)")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("list", "fetch", "build", "preview"):
        p = sub.add_parser(name)
        p.add_argument("--manifest", default=str(HERE / "manifest_c2_bright_clean2.csv"),
                       help="labelled frames (default: the QC-passed c2 bright frames of the ADCP winter)")
        p.add_argument("--out", type=Path, default=OUT, help="default composite/")
        p.add_argument("--s3", default=S3)
        p.add_argument("--every", type=int, default=1, help="use every N-th frame (quick tests)")
        p.add_argument("--recipe", default="snap,bright-dark,timex",
                       help="three channels R,G,B: product names or differences a-b")
        if name == "list":
            p.add_argument("--epoch", default=None, help="collection time to list (default: mid-manifest)")
        if name == "fetch":
            p.add_argument("--products", nargs="+", default=["bright", "dark", "snap", "timex"])
            p.add_argument("--dry-run", action="store_true")
        if name == "build":
            p.add_argument("--width", type=int, default=768)
            p.add_argument("--height", type=int, default=640)
        if name == "preview":
            p.add_argument("--hs", type=float, default=1.5, help="show a frame near this wave height")
    args = ap.parse_args()
    return {"list": cmd_list, "fetch": cmd_fetch, "build": cmd_build, "preview": cmd_preview}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())

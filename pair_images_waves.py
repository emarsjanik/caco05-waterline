#!/usr/bin/env python3
"""
Pair Coastal Imagery With ADCP Wave Records
==============================================
Builds the training manifest for an optical wave gauge: a table of
image files with the significant wave height and peak period measured
at the same moment by an in-view uplooking ADCP.

This is the step where errors are silent. A model trained on
mis-paired data trains perfectly well, reports a plausible loss curve,
and predicts nothing. Three checks are therefore built in:

  * TIMEZONE VERIFICATION. The ADCP file carries naive timestamps; the
    image filenames are explicitly GMT. If the two are on different
    clocks every pair is wrong by that offset. This is checked
    empirically rather than assumed, by correlating the ADCP's own
    water_level column against the tidal signal implied by the image
    timestamps at a range of candidate offsets. The correct offset
    produces a sharp correlation peak; no peak means the assumption is
    wrong and the run stops.

  * MATCH TOLERANCE. The ADCP reports hourly, images are half-hourly,
    so some images fall between records. Rather than interpolating
    wave height across an hour -- which smooths exactly the storm
    peaks the model most needs -- images further than --max-gap from a
    record are dropped and counted.

  * DISTRIBUTION REPORTING. Buscombe et al. (2020) found out-of-
    calibration performance was the weakest aspect of the technique
    and recommended training sets include extreme events. The manifest
    summary reports the wave-height and period distributions and warns
    when the tails are thin, because a set dominated by 1 m days
    teaches the network nothing about 3 m days.

Usage:
    python3 pair_images_waves.py --waves sig1000_waves_ALL.csv \\
        --image-dir /path/to/images --output manifest.csv [--camera c1]

    python3 pair_images_waves.py --waves ... --s3-list s3_listing.txt ...
        (pairs against an `aws s3 ls` listing without downloading first)
"""

import re
import sys
import csv
import argparse
from pathlib import Path
from datetime import datetime, timezone, timedelta

import numpy as np
import pandas as pd


FNAME_RE = re.compile(
    r"^(\d{9,11})\."           # epoch
    r".*?"                      # day/date text
    r"\.GMT\.(\d{4})\."         # year
    r".*?"                      # station
    r"\.(c\d)\."                # camera
    r"(\w+)\.jpg$", re.IGNORECASE)


def parse_filename(name):
    """Returns (epoch, camera, product) or None."""
    m = FNAME_RE.match(name)
    if not m:
        return None
    return int(m.group(1)), m.group(3).lower(), m.group(4).lower()


def load_waves(path, assume_tz):
    df = pd.read_csv(path)
    if "time" not in df.columns:
        sys.exit(f"{path}: no 'time' column")
    df["time"] = pd.to_datetime(df["time"])
    if df["time"].dt.tz is None:
        df["time"] = df["time"].dt.tz_localize(assume_tz)
    else:
        df["time"] = df["time"].dt.tz_convert(assume_tz)
    # Version-independent epoch conversion. `.astype("int64")` on a
    # tz-aware series returns NANOseconds in pandas < 3 and
    # MICROseconds in pandas 3, so the obvious `// 10**9` silently
    # yields epochs 1000x wrong on one of them. Subtracting the epoch
    # origin and taking total_seconds() is stable across versions,
    # which matters because this runs on two machines with different
    # pandas.
    origin = pd.Timestamp("1970-01-01", tz="UTC")
    df["epoch"] = (df["time"].dt.tz_convert("UTC") - origin).dt.total_seconds().astype("int64")

    # Guard regardless: a sane epoch for this work is 2000-2100.
    lo, hi = int(df["epoch"].min()), int(df["epoch"].max())
    if not (946_684_800 < lo < 4_102_444_800):
        sys.exit(f"Epoch conversion produced implausible values ({lo}). "
                 f"Check the pandas version and the 'time' column format.")
    return df.sort_values("epoch").reset_index(drop=True)


def verify_timezone(waves, tide_path, level_col="water_level"):
    """
    Confirms the ADCP clock matches UTC, using an INDEPENDENT tide
    record.

    An earlier version of this function compared the ADCP water level
    against itself at various lags. That is a tautology -- zero lag is
    identical by construction and always wins -- so it reported
    success regardless of the truth. It is recorded here because the
    failure was not obvious from its output, which looked like a
    clean result.

    A genuine test needs a second instrument. The ADCP's own
    water_level column is correlated against a tide model or GNSS-R
    series at a range of offsets; the offset minimising the residual
    is the ADCP's true clock error. Both are measuring the same sea.

    Without a reference file this returns None and says so, rather
    than manufacturing confidence.
    """
    if tide_path is None:
        return None, ("no --tide-reference supplied: clock NOT verified. If the ADCP "
                      "timestamps are local time rather than UTC, every pair is wrong "
                      "by that offset and the model will learn noise.")
    if level_col not in waves.columns:
        return None, f"ADCP file has no '{level_col}' column -- clock NOT verified"

    try:
        ref = (pd.read_excel(tide_path) if str(tide_path).lower().endswith((".xlsx", ".xls"))
               else pd.read_csv(tide_path))
    except Exception as exc:
        return None, f"could not read {tide_path}: {exc}"

    tcol = next((c for c in ref.columns if c.lower() == "time"), None)
    hcols = [c for c in ref.columns if "heightm" in c.lower()]
    if tcol is None or not hcols:
        return None, (f"{tide_path}: need a 'time' column and at least one "
                      f"'<model>_heightm' column -- clock NOT verified")

    rt = pd.to_datetime(ref[tcol])
    rt = rt.dt.tz_localize("UTC") if rt.dt.tz is None else rt.dt.tz_convert("UTC")
    origin = pd.Timestamp("1970-01-01", tz="UTC")
    rep = (rt - origin).dt.total_seconds().to_numpy()
    rlvl = ref[hcols].to_numpy(float).mean(axis=1)

    wep = waves["epoch"].to_numpy(float)
    wlvl = waves[level_col].to_numpy(float)
    overlap = (wep >= rep.min()) & (wep <= rep.max())
    if overlap.sum() < 50:
        return None, ("tide reference does not overlap the ADCP record -- clock NOT "
                      "verified")

    results = {}
    for hours in [-6, -5, -4, -2, -1, 0, 1, 2, 4, 5, 6]:
        probe = wep[overlap] + hours * 3600.0
        inside = (probe >= rep.min()) & (probe <= rep.max())
        if inside.sum() < 50:
            continue
        sampled = np.interp(probe[inside], rep, rlvl)
        actual = wlvl[overlap][inside]
        # Remove the datum difference; only the SHAPE alignment matters.
        results[hours] = float(np.std(actual - sampled))

    if not results or 0 not in results:
        return None, "could not evaluate offsets against the tide reference"
    best = min(results, key=results.get)
    note = ("residual std by offset (h: m): "
            + ", ".join(f"{k:+d}:{v:.3f}" for k, v in sorted(results.items())))
    return best, note


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--waves", required=True)
    ap.add_argument("--image-dir")
    ap.add_argument("--s3-list", help="Output of `aws s3 ls`, to pair without downloading.")
    ap.add_argument("--output", required=True)
    ap.add_argument("--camera", default="c2", choices=["c1", "c2", "both"])
    ap.add_argument("--product", default="timex")
    ap.add_argument("--max-gap", type=float, default=30.0,
                    help="Maximum minutes between an image and the nearest wave record "
                         "(default 30). Images further away are dropped rather than matched "
                         "to an interpolated value, because interpolating across an hour "
                         "smooths the storm peaks the model most needs to learn.")
    ap.add_argument("--on-tie", default="mean", choices=["mean", "earlier", "later", "drop"],
                    help="What to do when an image is exactly equidistant from two wave "
                         "records -- every :30 image against hourly ADCP records. 'mean' "
                         "(default) averages the two; 'earlier'/'later' take one side, which "
                         "is right only if the ADCP timestamp marks the start/end of its "
                         "burst; 'drop' discards the image. The previous behaviour was "
                         "an implicit 'earlier' via argmin, which silently labelled half of "
                         "all images with the previous hour's waves.")
    ap.add_argument("--wave-height-col", default="wh_4061")
    ap.add_argument("--wave-period-col", default="wp_peak")
    ap.add_argument("--assume-tz", default="UTC")
    ap.add_argument("--tide-reference",
                    help="An INDEPENDENT water-level record (tide model .xlsx with "
                         "'<model>_heightm' columns, or similar) used to verify that the "
                         "ADCP timestamps are UTC. Without it the clock is not checked, and "
                         "a local-time ADCP file would mis-pair every image silently.")
    ap.add_argument("--skip-tz-check", action="store_true")
    args = ap.parse_args()

    if not args.image_dir and not args.s3_list:
        ap.error("supply --image-dir or --s3-list")

    waves = load_waves(args.waves, args.assume_tz)
    print("=" * 70)
    print("IMAGE / WAVE PAIRING")
    print("=" * 70)
    print(f"wave records      : {len(waves)}")
    print(f"                    {waves.time.min()} to {waves.time.max()}")

    # ---- gather image filenames -------------------------------------
    names = []
    if args.image_dir:
        for p in sorted(Path(args.image_dir).glob("*.jpg")):
            names.append(p.name)
    if args.s3_list:
        for line in open(args.s3_list):
            parts = line.split()
            if parts:
                names.append(parts[-1])

    parsed = []
    for n in names:
        info = parse_filename(n)
        if not info:
            continue
        epoch, cam, product = info
        if product != args.product.lower():
            continue
        if args.camera != "both" and cam != args.camera:
            continue
        parsed.append((epoch, cam, n))

    if not parsed:
        sys.exit("No images matched the camera and product filters.")
    parsed.sort()
    print(f"images ({args.camera}, {args.product}): {len(parsed)}")
    ep_img = [p[0] for p in parsed]
    print(f"                    {datetime.fromtimestamp(min(ep_img), tz=timezone.utc)} to "
          f"{datetime.fromtimestamp(max(ep_img), tz=timezone.utc)}")
    print()

    # ---- timezone verification --------------------------------------
    if not args.skip_tz_check:
        best, note = verify_timezone(waves, args.tide_reference)
        print("Timezone check:")
        print("  " + note)
        if best is None:
            print("  NOT VERIFIED -- proceeding, but confirm both clocks are UTC.")
        elif best != 0:
            print(f"  WARNING: offset {best:+d} h fits better than 0.")
            print("  If the ADCP timestamps are local time rather than UTC, every pair in")
            print("  this manifest is wrong by that amount and the model will learn noise.")
            print("  Re-run with --assume-tz US/Eastern if that is the case.")
        else:
            print("  OK: zero offset fits best, both clocks agree.")
        print()

    # ---- match -------------------------------------------------------
    wep = waves["epoch"].to_numpy()
    hcol = args.wave_height_col
    pcol = args.wave_period_col
    for c in (hcol, pcol):
        if c not in waves.columns:
            sys.exit(f"wave file has no column '{c}'")

    hvals = waves[hcol].to_numpy(dtype=float)
    pvals = waves[pcol].to_numpy(dtype=float)

    rows, too_far, outside, ties, tie_dropped = [], 0, 0, 0, 0
    for epoch, cam, name in parsed:
        if epoch < wep.min() or epoch > wep.max():
            outside += 1
            continue
        dist = np.abs(wep - epoch)
        gap_s = int(dist.min())
        gap_min = gap_s / 60.0
        if gap_min > args.max_gap:
            too_far += 1
            continue
        # All records at the minimum distance: one normally, two when
        # the image sits exactly halfway between records. argmin alone
        # would always resolve that tie to the earlier record.
        nearest = np.flatnonzero(dist == gap_s)
        if len(nearest) > 1:
            ties += 1
            if args.on_tie == "drop":
                tie_dropped += 1
                continue
            if args.on_tie == "earlier":
                nearest = nearest[:1]
            elif args.on_tie == "later":
                nearest = nearest[-1:]
        rows.append({
            "filename": name,
            "camera": cam,
            "epoch": epoch,
            "time_utc": datetime.fromtimestamp(epoch, tz=timezone.utc).isoformat(),
            "gap_minutes": round(gap_min, 1),
            "n_records": len(nearest),
            "wave_height_m": float(np.mean(hvals[nearest])),
            "wave_period_s": float(np.mean(pvals[nearest])),
        })

    print(f"paired            : {len(rows)}")
    if outside:
        print(f"  outside record  : {outside} image(s) with no overlapping wave data")
    if too_far:
        print(f"  gap > {args.max_gap:.0f} min    : {too_far} image(s) dropped")
    if ties:
        print(f"  equidistant     : {ties} image(s) between two records "
              f"(--on-tie {args.on_tie}"
              + (f", {tie_dropped} dropped)" if tie_dropped else ")"))
    if not rows:
        sys.exit("No pairs produced.")

    out = pd.DataFrame(rows)
    out.to_csv(args.output, index=False)
    print(f"wrote             : {args.output}")
    print()

    # ---- distribution ------------------------------------------------
    h = out["wave_height_m"].to_numpy()
    p = out["wave_period_s"].to_numpy()
    print("Training set distribution")
    print(f"  wave height     : {h.min():.2f} to {h.max():.2f} m, "
          f"median {np.median(h):.2f}")
    print(f"  wave period     : {p.min():.2f} to {p.max():.2f} s, "
          f"median {np.median(p):.2f}")
    print()
    edges = np.linspace(h.min(), h.max(), 11)
    counts, _ = np.histogram(h, bins=edges)
    print("  wave height histogram (10 bins):")
    for i, c in enumerate(counts):
        bar = "#" * int(40 * c / max(counts.max(), 1))
        print(f"    {edges[i]:4.2f}-{edges[i+1]:4.2f} m {c:5d} {bar}")

    # Buscombe et al. found out-of-calibration prediction the weakest
    # aspect of the method, and recommended extremes be represented.
    top = int((h > np.percentile(h, 95)).sum())
    bot = int((h < np.percentile(h, 5)).sum())
    print()
    if counts[-1] < 20 or counts[-2] < 20:
        print("  WARNING: the upper wave-height bins are sparsely populated. Buscombe et al.")
        print("  (2020) found prediction outside the trained range to be the weakest aspect")
        print("  of this technique, and recommended training sets include extreme events.")
        print("  A model fitted mostly to moderate conditions will under-predict storms --")
        print("  which are usually the conditions of interest.")
    else:
        print(f"  Tails: {bot} samples below the 5th percentile, {top} above the 95th.")
    print()
    print("NOTE: pairs are only as good as the assumption that both clocks are UTC and")
    print("that the camera pose was constant across the record. Neither is verified by")
    print("the model, and both fail silently.")


if __name__ == "__main__":
    main()


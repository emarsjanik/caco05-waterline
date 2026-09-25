#!/usr/bin/env python3
"""
Beach Change Between Two Intertidal DEMs
===========================================
Differences two DEMs built by dem_from_contours.py (newer minus older)
and separates real change from measurement noise.

WHY. A DEM pooled over weeks blurs real erosion and accretion into its
"spread". Building DEMs over rolling windows (dem_from_contours.py
--last-days 7 --series-dir ...) and differencing them turns that spread
into measured change.

WHAT COUNTS AS CHANGE. Each DEM cell is the median of n repeat
crossings with a 16-84 percentile spread s, i.e. a per-crossing sigma
of about s/2. The standard error of a median is ~1.25 sigma/sqrt(n), so
a difference is called significant only where

    |B - A| > max(1.96 * sqrt(seA^2 + seB^2), --min-lod)

(95% level of detection). --min-lod (default 0.05 m) is a floor for
errors a single DEM cannot see: errors shared by every crossing in a
week (e.g. a week-long water-level or setup bias) do not show in s.
Cells are compared only where BOTH DEMs rest on at least --min-count
(default 5) tide crossings. At the edge of coverage a cell may hold only
3 or 4 crossings; their spread then badly understates the noise, the
level of detection comes out too small, and isolated "changes" of up to
1 m appeared along the edge of the first real Sep 2026 map.

For the same reason a change covering the whole overlap by a similar
amount is flagged: sand does not usually move uniformly; a water-level
or calibration shift does.

Usage:
    python3 dem_change.py --a archive/dems/dem_2026-09-14_7d --b archive/dems/dem_2026-09-21_7d
    python3 dem_change.py --series archive/dems --days 7     (newest vs one week earlier)

Outputs, next to B (or in --output-dir): change_<A>_to_<B>_diff.asc
(all overlapping cells), _sig.asc (significant change only), .png.
"""

import re
import sys
import argparse
from pathlib import Path
from datetime import date, timedelta

import numpy as np

SERIES_RE = re.compile(r"^dem_(\d{4}-\d{2}-\d{2})_(\d+)d_dem\.asc$")


def read_grid(path):
    hdr = {}
    with open(path) as f:
        for _ in range(6):
            k, v = f.readline().split()
            hdr[k.lower()] = float(v)
    g = np.loadtxt(path, skiprows=6, ndmin=2)
    nodata = hdr.get("nodata_value", -9999.0)
    g = np.where(g == nodata, np.nan, g)
    return np.flipud(g), hdr   # row 0 = south, like dem_from_contours


def write_grid(path, grid, x0, y0, cell, nodata=-9999.0):
    with open(path, "w") as f:
        f.write(f"ncols {grid.shape[1]}\nnrows {grid.shape[0]}\n")
        f.write(f"xllcorner {x0:.3f}\nyllcorner {y0:.3f}\ncellsize {cell}\nNODATA_value {nodata}\n")
        for r in np.flipud(grid):
            f.write(" ".join("%.4f" % (v if np.isfinite(v) else nodata) for v in r) + "\n")


def load_dem(stem):
    dem, h = read_grid(f"{stem}_dem.asc")
    spread, _ = read_grid(f"{stem}_spread.asc")
    count, _ = read_grid(f"{stem}_count.asc")
    return dem, spread, np.nan_to_num(count), h


def overlap(ha, hb):
    """Common extent of two grids on the same cell size, as index slices into each."""
    cell = ha["cellsize"]
    if abs(cell - hb["cellsize"]) > 1e-9:
        sys.exit("The two DEMs have different cell sizes.")
    x0 = max(ha["xllcorner"], hb["xllcorner"])
    y0 = max(ha["yllcorner"], hb["yllcorner"])
    x1 = min(ha["xllcorner"] + ha["ncols"] * cell, hb["xllcorner"] + hb["ncols"] * cell)
    y1 = min(ha["yllcorner"] + ha["nrows"] * cell, hb["yllcorner"] + hb["nrows"] * cell)
    nc, nr = int(round((x1 - x0) / cell)), int(round((y1 - y0) / cell))
    if nc <= 0 or nr <= 0:
        sys.exit("The two DEMs do not overlap.")

    def sl(h):
        c0 = int(round((x0 - h["xllcorner"]) / cell))
        r0 = int(round((y0 - h["yllcorner"]) / cell))
        return (slice(r0, r0 + nr), slice(c0, c0 + nc))
    return sl(ha), sl(hb), x0, y0, cell


def compare(stem_a, stem_b, min_lod, out_dir=None, plot=True, min_count=5):
    A, sA, nA, ha = load_dem(stem_a)
    B, sB, nB, hb = load_dem(stem_b)
    ia, ib, x0, y0, cell = overlap(ha, hb)
    A, sA, nA = A[ia], sA[ia], nA[ia]
    B, sB, nB = B[ib], sB[ib], nB[ib]

    in_both = np.isfinite(A) & np.isfinite(B) & (nA > 0) & (nB > 0)
    both = in_both & (nA >= min_count) & (nB >= min_count)
    diff = np.where(both, B - A, np.nan)
    seA = 1.25 * (np.nan_to_num(sA) / 2.0) / np.sqrt(np.maximum(nA, 1))
    seB = 1.25 * (np.nan_to_num(sB) / 2.0) / np.sqrt(np.maximum(nB, 1))
    lod = np.maximum(1.96 * np.sqrt(seA ** 2 + seB ** 2), min_lod)
    sig = both & (np.abs(diff) > lod)
    sig_diff = np.where(sig, diff, np.nan)

    name_a, name_b = Path(stem_a).name, Path(stem_b).name
    out_dir = Path(out_dir) if out_dir else Path(stem_b).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"change_{name_a.replace('dem_', '')}_to_{name_b.replace('dem_', '')}"
    write_grid(f"{out}_diff.asc", diff, x0, y0, cell)
    write_grid(f"{out}_sig.asc", sig_diff, x0, y0, cell)

    n_both = int(both.sum())
    print("=" * 72)
    print(f"BEACH CHANGE  {name_a}  ->  {name_b}")
    print("=" * 72)
    print(f"cells in both DEMs      : {n_both}  ({n_both * cell * cell:.0f} m2) with >= {min_count} "
          f"crossings in each; {int(in_both.sum()) - n_both} thinner cell(s) not compared")
    if n_both == 0:
        print("No overlapping cells -- nothing to compare.")
        return None
    d = diff[both]
    print(f"change, all overlap     : median {np.median(d):+.3f} m, "
          f"mean {np.mean(d):+.3f} m, p10 {np.percentile(d, 10):+.3f}, p90 {np.percentile(d, 90):+.3f}")
    print(f"level of detection      : median {np.median(lod[both]):.3f} m (95%, floor {min_lod} m)")
    ero, acc = sig & (diff < 0), sig & (diff > 0)
    area = cell * cell
    print(f"significant erosion     : {int(ero.sum())} cells, {ero.sum() * area:.0f} m2, "
          f"{np.nansum(diff[ero]) * area:+.1f} m3")
    print(f"significant accretion   : {int(acc.sum())} cells, {acc.sum() * area:.0f} m2, "
          f"{np.nansum(diff[acc]) * area:+.1f} m3")
    print(f"net (significant cells) : {np.nansum(diff[sig]) * area:+.1f} m3 over "
          f"{100 * sig.sum() / n_both:.0f}% of the overlap")
    same_sign = max((d > 0).mean(), (d < 0).mean())
    if abs(np.median(d)) > max(min_lod, 0.05) and same_sign > 0.8:
        print()
        print(f"WARNING: {100 * same_sign:.0f}% of overlapping cells moved the same way, by a median "
              f"{np.median(d):+.2f} m. Uniform change is more typical of a shift in water level, "
              "wave setup or camera calibration between the two weeks than of sand moving. "
              "Check before interpreting as erosion/accretion.")
    print(f"wrote {out}_diff.asc, {out}_sig.asc")

    if plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            ext = [x0, x0 + diff.shape[1] * cell, y0, y0 + diff.shape[0] * cell]
            lim = max(0.1, float(np.nanpercentile(np.abs(d), 98)))
            fig, axes = plt.subplots(1, 2, figsize=(14, 7), dpi=110)
            for ax, grid, title in ((axes[0], diff, "all overlapping cells"),
                                    (axes[1], sig_diff, "significant change only (95% LoD)")):
                im = ax.imshow(grid, origin="lower", extent=ext, cmap="RdBu", vmin=-lim, vmax=lim,
                               aspect="equal")
                ax.set_title(f"{name_b} minus {name_a}\n{title}", fontsize=10)
                ax.set_xlabel("easting (m, UTM 19N)"); ax.set_ylabel("northing (m, UTM 19N)")
                ax.ticklabel_format(useOffset=False, style="plain")
                plt.colorbar(im, ax=ax, shrink=0.8, label="elevation change (m); blue = accretion")
            plt.tight_layout()
            plt.savefig(f"{out}.png", bbox_inches="tight")
            print(f"wrote {out}.png")
        except Exception as exc:
            print(f"(plot skipped: {exc})")
    return out


def pick_series(series_dir, days):
    """Newest dated DEM, and the newest one ending at least `days` earlier."""
    found = []
    for p in Path(series_dir).glob("dem_*_*d_dem.asc"):
        m = SERIES_RE.match(p.name)
        if m and int(m.group(2)) == days:
            found.append((date.fromisoformat(m.group(1)), str(p)[: -len("_dem.asc")]))
    if not found:
        return None, None
    found.sort()
    newest_date, newest = found[-1]
    earlier = [s for d, s in found if d <= newest_date - timedelta(days=days)]
    return (earlier[-1] if earlier else None), newest


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--a", help="Older DEM stem (without _dem.asc)")
    ap.add_argument("--b", help="Newer DEM stem")
    ap.add_argument("--series", help="Folder of dem_<date>_<N>d_* files: compare the newest "
                                     "with the newest one at least --days earlier.")
    ap.add_argument("--days", type=int, default=7)
    ap.add_argument("--min-lod", type=float, default=0.05,
                    help="Floor on the level of detection, m (default 0.05).")
    ap.add_argument("--min-count", type=int, default=5,
                    help="Compare only cells with at least this many tide crossings in BOTH "
                         "DEMs (default 5). Thinner cells have unreliable spreads.")
    ap.add_argument("--output-dir")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    if args.series:
        a, b = pick_series(args.series, args.days)
        if b is None:
            print(f"No {args.days}-day DEMs in {args.series} yet."); return
        if a is None:
            print(f"Only DEMs less than {args.days} days apart so far (newest {Path(b).name}); "
                  "no change map yet."); return
    elif args.a and args.b:
        a, b = args.a, args.b
    else:
        ap.error("give --a and --b, or --series")
    compare(a, b, args.min_lod, args.output_dir, not args.no_plot, args.min_count)


if __name__ == "__main__":
    main()

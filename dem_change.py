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
or calibration shift does. But only on --min-uniform-cells (50) or more
compared cells, and every change statement says how many cells it
rests on.

THE REFERENCE, AND WHY IT IS REBUILT. On 8 Oct 2026 the cron printed
"100% of overlapping cells moved the same way, by a median -0.75 m"
against archive/dems/dem_2026-09-27_7d. That DEM was degenerate -- 144
cells, all +0.58..+0.64 m, from storm days -- and shared 7 cells with
the new one; and it had been built by older processing (before the wave
setup, the consistency filter and the datum fix). Rebuilt like for like
from the current contour file (17-22 Sep against 1-6 Oct), the change
was real and not uniform: median -0.12 m over 354 cells. So:
  * --rebuild GROUND.csv (the cron's): the reference window is built
    again from the CURRENT contour_points_ground.csv, with the same
    dem_from_contours.py settings as the new window DEM (everything
    after '--' on the command line), into a temporary stem. The window
    ends --days before the new one's end, so the two do not overlap.
    The dated copies in the series folder stay as the record of what
    each run built; they are not the reference.
  * A reference with fewer than --min-ref-cells (200) cells of
    --min-count crossings, or an elevation range (5-95 percentile) under
    --min-ref-range (0.5 m), is DEGENERATE: a few cells at one level
    cannot show how a beach changed (the station's 7-day DEM of 29 Sep
    - 5 Oct 2026 has 973 such cells over -0.82 to +1.00 m). It is
    skipped, saying so, and the window ending one day earlier is tried,
    back to --lookback (7) days.
    Without --rebuild the archived DEMs are used, nearest adequate one
    first. A degenerate NEW DEM gets no change statement at all.

Usage:
    python3 dem_change.py --a archive/dems/dem_2026-09-14_7d --b archive/dems/dem_2026-09-21_7d
    python3 dem_change.py --series archive/dems --days 7     (newest vs one week earlier)
    python3 dem_change.py --series archive/dems --days 7 --rebuild contour_points_ground.csv \
        -- --cell 2 --min-points 3 --max-spread 0.5 --max-hs 2.5   (the cron's: like for like)

Outputs, next to B (or in --output-dir): change_<A>_to_<B>_diff.asc
(all overlapping cells), _sig.asc (significant change only), .png.
"""

import re
import csv
import sys
import shutil
import argparse
import subprocess
import tempfile
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


def compare(stem_a, stem_b, min_lod, out_dir=None, plot=True, min_count=5,
            min_uniform_cells=50, name_a=None, ref_note=None):
    """
    B minus A over the cells both rest on >= min_count crossings. `name_a`
    names the reference in the outputs (default: its stem); `ref_note`
    says where it came from (printed and on the figure).
    """
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
    with np.errstate(invalid="ignore"):          # NaN outside the overlap
        sig = both & (np.abs(diff) > lod)
    sig_diff = np.where(sig, diff, np.nan)

    name_a, name_b = name_a or Path(stem_a).name, Path(stem_b).name
    out_dir = Path(out_dir) if out_dir else Path(stem_b).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"change_{name_a.replace('dem_', '')}_to_{name_b.replace('dem_', '')}"
    write_grid(f"{out}_diff.asc", diff, x0, y0, cell)
    write_grid(f"{out}_sig.asc", sig_diff, x0, y0, cell)

    n_both = int(both.sum())
    print("=" * 72)
    print(f"BEACH CHANGE  {name_a}  ->  {name_b}")
    print("=" * 72)
    if ref_note:
        print(f"reference               : {ref_note}")
    print(f"cells in both DEMs      : {n_both}  ({n_both * cell * cell:.0f} m2) with >= {min_count} "
          f"crossings in each; {int(in_both.sum()) - n_both} thinner cell(s) not compared")
    if n_both == 0:
        print("No overlapping cells -- nothing to compare.")
        return None
    d = diff[both]
    print(f"{'change, ' + str(n_both) + ' cells':<24}: median {np.median(d):+.3f} m, "
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
    if n_both < min_uniform_cells:
        # 7 cells of a degenerate reference once read as '100% moved the
        # same way': a handful of cells cannot say whether change was uniform.
        print()
        print(f"NOTE: only {n_both} cell(s) compared (fewer than {min_uniform_cells}): too few to "
              "say whether the change was uniform or where the beach moved; the numbers above "
              f"rest on those {n_both} cell(s) alone.")
    elif abs(np.median(d)) > max(min_lod, 0.05) and same_sign > 0.8:
        print()
        print(f"WARNING: {100 * same_sign:.0f}% of the {n_both} compared cells moved the same way, by "
              f"a median {np.median(d):+.2f} m. Uniform change is more typical of a shift in water "
              "level, wave setup or camera calibration between the two weeks than of sand moving. "
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
                ax.set_title(f"{name_b} minus {name_a}: {n_both} cells compared\n{title}",
                             fontsize=10)
                ax.set_xlabel("easting (m, UTM 19N)"); ax.set_ylabel("northing (m, UTM 19N)")
                ax.ticklabel_format(useOffset=False, style="plain")
                plt.colorbar(im, ax=ax, shrink=0.8, label="elevation change (m); blue = accretion")
            plt.tight_layout()
            if ref_note:
                import textwrap
                fig.text(0.02, 0.0, "\n".join(textwrap.wrap(f"Reference: {ref_note}", 150)),
                         fontsize=9, va="top", ha="left")
            plt.savefig(f"{out}.png", bbox_inches="tight")
            print(f"wrote {out}.png")
        except Exception as exc:
            print(f"(plot skipped: {exc})")
    return out


def series_dems(series_dir, days):
    """[(end date, stem)] of the dated `days`-day DEMs in the folder, oldest first."""
    found = []
    for p in Path(series_dir).glob("dem_*_*d_dem.asc"):
        m = SERIES_RE.match(p.name)
        if m and int(m.group(2)) == days:
            found.append((date.fromisoformat(m.group(1)), str(p)[: -len("_dem.asc")]))
    return sorted(found)


def pick_series(series_dir, days):
    """Newest dated DEM, and the newest one ending at least `days` earlier."""
    found = series_dems(series_dir, days)
    if not found:
        return None, None
    newest_date, newest = found[-1]
    earlier = [s for d, s in found if d <= newest_date - timedelta(days=days)]
    return (earlier[-1] if earlier else None), newest


def dem_summary(stem, min_count):
    """Cells with an elevation and >= min_count crossings, and their 5-95% elevation range."""
    dem, _, count, _ = load_dem(stem)
    ok = np.isfinite(dem) & (count >= min_count)
    n = int(ok.sum())
    if n == 0:
        return n, 0.0, float("nan"), float("nan")
    p5, p95 = np.percentile(dem[ok], [5, 95])
    return n, float(p95 - p5), float(p5), float(p95)


def degenerate(stem, min_count, min_cells, min_range):
    """None if the DEM can serve in a change map, else why not (one line)."""
    n, rng, p5, p95 = dem_summary(stem, min_count)
    if n < min_cells:
        return (f"{n} cell(s) with >= {min_count} crossings (fewer than {min_cells})"
                + (f", {p5:+.2f} to {p95:+.2f} m" if n else ""))
    if rng < min_range:
        return (f"elevations {p5:+.2f} to {p95:+.2f} m (5-95%) over {n} cells: a range of "
                f"{rng:.2f} m, under {min_range} m -- one level of the beach, not a beach")
    return None


def window_subset(ground, start, end, out_csv):
    """
    Copy the rows of `ground` captured from `start` to `end` (UTC dates,
    inclusive) to `out_csv`, header kept. One pass over the whole-archive
    file; each reference candidate is then built from this small file, so
    trying several costs one read of the archive, not several.
    """
    n = 0
    with open(ground, newline="") as f, open(out_csv, "w", newline="") as g:
        r = csv.reader(f)
        w = csv.writer(g)
        header = next(r)
        w.writerow(header)
        if "capture_time_utc" not in header:
            return 0
        k = header.index("capture_time_utc")
        for row in r:
            if len(row) > k and start <= row[k][:10] <= end:
                w.writerow(row)
                n += 1
    return n


def rebuild_reference(ground, new_end, days, build_args, tmp, lookback, min_count, min_cells,
                      min_range):
    """
    Build the reference window again from the current ground file, with
    the new DEM's dem_from_contours.py settings (`build_args`): the window
    ending `days` before `new_end`, or, if that one is degenerate, the
    nearest earlier one that is not (one day at a time, back `lookback`
    days). Returns (stem, start, end) or None, printing what was tried.
    """
    here = Path(__file__).resolve().parent
    last_end = new_end - timedelta(days=days)
    first_start = last_end - timedelta(days=lookback + days - 1)
    sub = Path(tmp) / "reference_rows.csv"
    n = window_subset(ground, first_start.isoformat(), last_end.isoformat(), sub)
    print(f"reference rebuild       : {n:,} row(s) of {Path(ground).name} from "
          f"{first_start} to {last_end}")
    if n == 0:
        print(f"  no waterline points in the {lookback + days} days before {last_end + timedelta(days=1)}")
        return None
    for back in range(lookback + 1):
        end = last_end - timedelta(days=back)
        start = end - timedelta(days=days - 1)
        stem = Path(tmp) / f"dem_{end.isoformat()}_{days}d"
        cmd = [sys.executable, str(here / "dem_from_contours.py"), str(sub), str(stem),
               "--start-date", start.isoformat(), "--end-date", end.isoformat(),
               "--no-plot"] + list(build_args)
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              universal_newlines=True)
        label = f"  {start} to {end}"
        if proc.returncode != 0 or not Path(f"{stem}_dem.asc").exists():
            tail = [ln for ln in proc.stdout.strip().splitlines() if ln.strip()][-1:] or ["?"]
            print(f"{label}: not built (exit {proc.returncode}: {tail[0].strip()})")
            continue
        why = degenerate(stem, min_count, min_cells, min_range)
        if why:
            print(f"{label}: SKIPPED, degenerate -- {why}")
            continue
        n_cells, rng, p5, p95 = dem_summary(stem, min_count)
        print(f"{label}: {n_cells} cells with >= {min_count} crossings, {p5:+.2f} to {p95:+.2f} m"
              + ("" if back == 0 else f"  (nearest adequate window, {back} day(s) earlier)"))
        return str(stem), start, end
    print(f"  no adequate reference window ending {last_end - timedelta(days=lookback)} to {last_end}")
    return None


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
    ap.add_argument("--min-uniform-cells", type=int, default=50,
                    help="Say whether the change was uniform only on at least this many "
                         "compared cells (default 50).")
    ap.add_argument("--min-ref-cells", type=int, default=200,
                    help="A DEM with fewer cells of --min-count crossings is degenerate: no "
                         "reference, no change statement (default 200).")
    ap.add_argument("--min-ref-range", type=float, default=0.5,
                    help="...nor one whose elevations span (5-95%%) less than this, m (default 0.5).")
    ap.add_argument("--rebuild", metavar="GROUND_CSV",
                    help="With --series: build the reference window again from this (current) "
                         "contour_points_ground.csv with the dem_from_contours.py settings given "
                         "after '--', instead of using the archived DEM.")
    ap.add_argument("--lookback", type=int, default=7,
                    help="With --series: if the reference window is degenerate, try windows "
                         "ending up to this many days earlier (default 7).")
    ap.add_argument("build_args", nargs=argparse.REMAINDER,
                    help="After '--': the new DEM's dem_from_contours.py settings, for --rebuild "
                         "(e.g. -- --cell 2 --min-points 3 --max-spread 0.5 --max-hs 2.5).")
    args = ap.parse_args()
    build_args = [x for x in args.build_args if x != "--"]
    adequacy = (args.min_count, args.min_ref_cells, args.min_ref_range)

    if args.a and args.b:
        compare(args.a, args.b, args.min_lod, args.output_dir, not args.no_plot, args.min_count,
                args.min_uniform_cells)
        return
    if not args.series:
        ap.error("give --a and --b, or --series")

    found = series_dems(args.series, args.days)
    if not found:
        print(f"No {args.days}-day DEMs in {args.series} yet."); return
    b_end, b = found[-1]
    why = degenerate(b, *adequacy)
    if why:
        print(f"Newest {args.days}-day DEM {Path(b).name} is degenerate -- {why}. No change "
              "statement this run.")
        return

    if args.rebuild:
        tmp = tempfile.mkdtemp(prefix="dem_change_ref_")
        try:
            ref = rebuild_reference(args.rebuild, b_end, args.days, build_args, tmp,
                                    args.lookback, *adequacy)
            if ref is None:
                print("No change map this run: no adequate reference window.")
                return
            stem, start, end = ref
            archived = Path(args.series) / f"dem_{end.isoformat()}_{args.days}d_dem.asc"
            note = (f"{start} to {end}, rebuilt from the current {Path(args.rebuild).name} with "
                    f"the new DEM's settings ({' '.join(build_args) or 'defaults'})"
                    + (f"; archived {archived.name[:-len('_dem.asc')]} kept as a record, not used"
                       if archived.exists() else ""))
            compare(stem, b, args.min_lod, args.output_dir, not args.no_plot, args.min_count,
                    args.min_uniform_cells, name_a=f"dem_{end.isoformat()}_{args.days}d",
                    ref_note=note)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        return

    earlier = [(d, s) for d, s in found if d <= b_end - timedelta(days=args.days)]
    if not earlier:
        print(f"Only DEMs less than {args.days} days apart so far (newest {Path(b).name}); "
              "no change map yet."); return
    for k, (d, a) in enumerate(reversed(earlier)):
        if (b_end - timedelta(days=args.days)) - d > timedelta(days=args.lookback):
            break
        why = degenerate(a, *adequacy)
        if why:
            print(f"reference {Path(a).name}: SKIPPED, degenerate -- {why}")
            continue
        note = (f"archived {Path(a).name} (built by the processing of its day, which may "
                "differ from today's; --rebuild compares like for like)")
        if k:
            note += f"; nearest adequate, {k} newer one(s) skipped"
        compare(a, b, args.min_lod, args.output_dir, not args.no_plot, args.min_count,
                args.min_uniform_cells, ref_note=note)
        return
    print(f"No change map this run: no adequate {args.days}-day DEM ending {args.days} to "
          f"{args.days + args.lookback} days before {Path(b).name}.")


if __name__ == "__main__":
    main()

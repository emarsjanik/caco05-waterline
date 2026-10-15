#!/usr/bin/env python3
"""
Plain DEM Map
-------------
The DEM as a plain map: north up, true scale, UTM 19N metres (the grids'
frame, NAD83(2011) / UTM 19N, EPSG:6348), elevation on the left with
0.5 m contours, the repeatability (16-84 percentile spread of the frames
in each cell) on the right. Only cells with a value are coloured; nothing
else is drawn. This is the map the station's DEMs had before the
diagnostic page (dem_figure.py), which Dr Sherwood found easier to read
(Oct 2026); survey_products.py writes both.

    python3 dem_map.py <stem> [--title TEXT] [--output PNG]

reads <stem>_dem.asc and <stem>_spread.asc (dem_from_contours.py) and
writes <stem>_map.png.
"""

import argparse
import sys
from pathlib import Path

import numpy as np


def read_asc(path):
    """ESRI ASCII grid -> (grid with NaN for no data, north row first; x0, y_top, cell)."""
    head = {}
    with open(path) as f:
        for _ in range(6):
            k, v = f.readline().split()
            head[k.lower()] = float(v)
        grid = np.loadtxt(f, ndmin=2)
    grid[grid == head["nodata_value"]] = np.nan
    cell = head["cellsize"]
    return grid, head["xllcorner"], head["yllcorner"] + head["nrows"] * cell, cell


def crop(grid, x0, y_top, cell, pad=3):
    """The grid cut to its cells with a value (plus `pad` cells), and its extent."""
    ok = np.isfinite(grid)
    if not ok.any():
        return grid, (x0, x0 + grid.shape[1] * cell, y_top - grid.shape[0] * cell, y_top)
    rows, cols = np.flatnonzero(ok.any(axis=1)), np.flatnonzero(ok.any(axis=0))
    r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad + 1, grid.shape[0])
    c0, c1 = max(cols[0] - pad, 0), min(cols[-1] + pad + 1, grid.shape[1])
    return grid[r0:r1, c0:c1], (x0 + c0 * cell, x0 + c1 * cell, y_top - r1 * cell, y_top - r0 * cell)


def draw(stem, title=None, output=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    stem = Path(stem)
    dem, x0, yt, cell = read_asc(f"{stem}_dem.asc")
    spread_path = Path(f"{stem}_spread.asc")
    spread = read_asc(spread_path)[0] if spread_path.exists() else None
    box = crop(dem, x0, yt, cell)[1]
    sel = lambda g: crop(g, x0, yt, cell)[0] if g is not None else None   # noqa: E731
    d, s = sel(dem), sel(spread)
    n = int(np.isfinite(dem).sum())
    if not title:                         # the window from the DEM's own record, when there is one
        try:
            import json
            info = json.loads(Path(f"{stem}_info.json").read_text())
            a, b = info.get("first_date") or info.get("window_start"), info.get("last_date") or info.get("window_end")
            if a and b:
                title = f"Marconi Beach: {a} to {b}" + (f" (last {info['last_days']} days)" if info.get("last_days") else "")
        except (OSError, ValueError):
            pass
    fig, axes = plt.subplots(1, 2 if s is not None else 1, figsize=(12, 8.5), squeeze=False)
    ax = axes[0][0]
    im = ax.imshow(d, extent=box, origin="upper", cmap="terrain", interpolation="nearest")
    if np.isfinite(d).sum() > 3 and np.nanmax(d) - np.nanmin(d) >= 0.5:
        X = np.linspace(box[0] + cell / 2, box[1] - cell / 2, d.shape[1])
        Y = np.linspace(box[3] - cell / 2, box[2] + cell / 2, d.shape[0])
        lv = np.arange(np.floor(np.nanmin(d) * 2) / 2, np.nanmax(d) + 0.01, 0.5)
        cs = ax.contour(X, Y, np.ma.masked_invalid(d), levels=lv, colors="k", linewidths=0.6)
        ax.clabel(cs, fmt="%+.1f", fontsize=7)
    fig.colorbar(im, ax=ax, shrink=0.7, label="elevation (m NAVD88)")
    ax.set_title(f"{title or stem.name}\nIntertidal DEM, {n} cells at {cell:g} m; contours every 0.5 m",
                 fontsize=10)
    if s is not None:
        ax2 = axes[0][1]
        im2 = ax2.imshow(s, extent=box, origin="upper", cmap="magma", vmin=0, vmax=0.5,
                         interpolation="nearest")
        fig.colorbar(im2, ax=ax2, shrink=0.7, label="spread (m)")
        ax2.set_title("Repeatability (precision, not accuracy)\n16-84 percentile range of the frames in each cell",
                      fontsize=10)
    for a in axes[0]:
        a.set_aspect("equal")
        a.set_xlabel("easting (m, NAD83(2011) / UTM 19N, EPSG:6348)")
        a.set_ylabel("northing (m)")
        a.ticklabel_format(useOffset=False, style="plain")
        a.tick_params(labelsize=8)
        a.annotate("N", xy=(0.93, 0.97), xytext=(0.93, 0.88), xycoords="axes fraction",
                   arrowprops={"arrowstyle": "->"}, ha="center", fontsize=9)
    fig.text(0.5, 0.02, "Elevation = water level + wave setup where the waterlines crossed, the median of 3 or "
             "more frames per cell; only cells with a value are coloured. Grids: <date>_dem.tif, _spread.tif, "
             "_count.tif.", ha="center", fontsize=8, color="0.3")
    out = Path(output) if output else Path(f"{stem}_map.png")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stem", help="the DEM's stem (<stem>_dem.asc, <stem>_spread.asc)")
    ap.add_argument("--title", default=None)
    ap.add_argument("--output", default=None, help="default <stem>_map.png")
    a = ap.parse_args()
    if not Path(f"{a.stem}_dem.asc").exists():
        sys.exit(f"no {a.stem}_dem.asc")
    print(f"wrote {draw(a.stem, a.title, a.output)}")


if __name__ == "__main__":
    main()

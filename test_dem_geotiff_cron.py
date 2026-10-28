#!/usr/bin/env python3
"""
Tests For The Cron's 7-Day DEM GeoTIFF
========================================
Plain python3 (no pytest needed on the station):

    python3 test_dem_geotiff_cron.py

WHAT. waterline_timex_cron.sh writes the 7-day DEM as a GeoTIFF next to its
.asc (dem_geotiff) for the GPS daily email to attach:
  * the cron calls it on the 7-day DEM whenever that DEM's grids were
    written (dem_from_contours.py exit 0 or 3), and nowhere else;
  * the GeoTIFF holds the .asc's values (nodata -9999 where the grid has
    none), its 2 m cells and corner, and says NAD83(2011) / UTM 19N
    (EPSG:6348) with NAVD88 heights (EPSG:5703);
  * a new run replaces the old file; no scratch folder is left behind;
  * a conversion that fails, or no .asc, removes the old GeoTIFF with a
    WARNING, so the email says it is missing rather than sending an
    older surface as the day's.

HOW. The function is taken from the cron's own text and called by bash in
a temporary folder (asc_to_geotiff.py copied in as $BASE), with log()
writing to a file.

STATION FILES. None are read.
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from compare_dem_survey import read_tiff   # noqa: E402

CRON = (HERE / "waterline_timex_cron.sh").read_text()
failures = []


def check(name, cond, detail=""):
    print(f"  {'ok ' if cond else 'FAIL'} {name}" + (f"  ({detail})" if detail != "" and not cond else ""))
    if not cond:
        failures.append(name)


def function_text(name):
    m = re.search(r"^%s\(\) \{\n.*?^\}\n" % re.escape(name), CRON, re.S | re.M)
    return m.group(0) if m else ""


def run(base, asc):
    script = (f'BASE="{base}"\nLOG="{base}/cron.log"\n'
              'log() { echo "LOG $*" >> "$LOG"; }\n' + function_text("dem_geotiff") + f'dem_geotiff "{asc}"\n')
    # python3 in the function is the one running this test (the station's own on the station)
    bindir = Path(base) / "bin"
    bindir.mkdir(exist_ok=True)
    if not (bindir / "python3").exists():
        (bindir / "python3").write_text(f'#!/bin/sh\nexec "{sys.executable}" "$@"\n')
        (bindir / "python3").chmod(0o755)
    env = dict(os.environ, PATH=f"{bindir}:{os.environ.get('PATH', '')}")
    r = subprocess.run(["bash", "-c", script], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True, env=env)
    log = (Path(base) / "cron.log").read_text() if (Path(base) / "cron.log").exists() else ""
    return r.returncode, r.stdout, log


def write_asc(path, grid, xll=420000.0, yll=4638000.0, cell=2.0):
    rows = ["ncols %d" % grid.shape[1], "nrows %d" % grid.shape[0], f"xllcorner {xll}", f"yllcorner {yll}",
            f"cellsize {cell}", "NODATA_value -9999"]
    rows += [" ".join("-9999" if not np.isfinite(v) else f"{v:.3f}" for v in r) for r in grid]
    Path(path).write_text("\n".join(rows) + "\n")


def geokey(tags, key):
    k = tags.get(34735, ())
    for i in range(4, len(k) - 3, 4):
        if k[i] == key:
            return k[i + 3]
    return None


def main():
    print("waterline_timex_cron.sh: the 7-day DEM as a GeoTIFF")
    fn = function_text("dem_geotiff")
    check("the cron defines dem_geotiff", bool(fn))
    calls = [m.start() for m in re.finditer(r"^\s*dem_geotiff ", CRON, re.M)]
    win = CRON.find('if [ $win_rc -eq 0 ] || [ $win_rc -eq 3 ]; then')
    end = CRON.find("\n                    else", win)
    check("called once, on the 7-day DEM, inside the block run when its grids were written (exit 0 or 3)",
          len(calls) == 1 and win >= 0 and win < calls[0] < end
          and 'dem_geotiff "${DEM_STEM}_${DEM_WINDOW_DAYS}day_dem.asc"' in CRON[calls[0]:calls[0] + 80],
          (calls, win, end))

    rng = np.random.default_rng(1)
    grid = np.round(rng.uniform(-1.0, 3.0, (40, 60)), 3)
    grid[:5, :] = np.nan
    grid[10:12, 20:30] = np.nan
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        shutil.copy(HERE / "asc_to_geotiff.py", base / "asc_to_geotiff.py")
        asc = base / "dem_intertidal_7day_dem.asc"
        tif = base / "dem_intertidal_7day_dem.tif"
        write_asc(asc, grid)
        rc, out, log = run(base, asc)
        ok = tif.exists()
        check("written next to the .asc, logged", rc == 0 and ok
              and "LOG window DEM GeoTIFF written: dem_intertidal_7day_dem.tif" in log, (rc, out, log))
        if ok:
            data, tags = read_tiff(tif)
            same = data.shape == grid.shape and np.allclose(np.where(np.isfinite(grid), grid, -9999.0),
                                                            np.where(np.isfinite(data), data, -9999.0), atol=1e-3)
            check("the .asc's values, nodata where it has none", same, (data.shape, grid.shape))
            check("nodata tag -9999", tags.get(42113, "").strip() == "-9999", tags.get(42113))
            scale, tie = tags.get(33550, ()), tags.get(33922, ())
            check("2 m cells, the upper-left corner of the grid",
                  tuple(scale[:2]) == (2.0, 2.0) and tie[3:5] == (420000.0, 4638000.0 + 2.0 * 40), (scale, tie))
            check("NAD83(2011) / UTM 19N with NAVD88 heights (EPSG 6348 + 5703)",
                  geokey(tags, 3072) == 6348 and geokey(tags, 4096) == 5703,
                  (geokey(tags, 3072), geokey(tags, 4096)))
        check("no scratch folder left", not list(base.glob(".geotiff.*")), list(base.glob(".geotiff.*")))

        write_asc(asc, grid + 1.0)
        run(base, asc)
        data, _ = read_tiff(tif) if tif.exists() else (None, None)
        check("a new run replaces it", data is not None and abs(np.nanmax(data) - (np.nanmax(grid) + 1.0)) < 1e-3)

        asc.write_text("not a grid\n")
        (base / "cron.log").unlink()
        rc, out, log = run(base, asc)
        check("a conversion that fails: the old GeoTIFF removed, a WARNING, no scratch folder",
              rc == 0 and not tif.exists() and "LOG WARNING: window DEM GeoTIFF NOT written" in log
              and not list(base.glob(".geotiff.*")), (rc, log))

        write_asc(asc, grid)
        run(base, asc)
        asc.unlink()
        (base / "cron.log").unlink()
        rc, out, log = run(base, asc)
        check("no .asc: the old GeoTIFF removed, a WARNING", rc == 0 and not tif.exists()
              and "LOG WARNING: no dem_intertidal_7day_dem.asc" in log, (rc, log))

    if failures:
        print(f"{len(failures)} check(s) FAILED")
        raise SystemExit(1)
    print("all DEM GeoTIFF cron tests passed")


if __name__ == "__main__":
    main()

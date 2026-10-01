#!/usr/bin/env python3
"""
ASCII Grids To GeoTIFF
=========================
Converts the ESRI ASCII grids written by dem_from_contours.py (<stem>_dem.asc,
_spread.asc, _count.asc, _source.asc) into GeoTIFFs that open directly in
ArcGIS, QGIS, GDAL and rasterio, georeferenced in UTM zone 19N.

Pure numpy: no GDAL, rasterio or scipy needed (the station's scipy does not
import). The files are baseline float32 GeoTIFFs: one band, uncompressed,
GeoKeys for the projected coordinate system (EPSG code), the pixel size and
upper-left corner, and the GDAL nodata tag.

Values: _dem and _spread in metres (elevation NAVD88 for _dem), _count in
frames per cell, _source 1 = measured, 2 = interpolated. Nodata = -9999.

COORDINATE SYSTEM. Default EPSG:32619 (WGS 84 / UTM zone 19N), as the
project's existing .prj file. If the camera survey was in NAD83(2011) use
--epsg 6348 (NAD83(2011) / UTM zone 19N); the two differ by about a metre.

Usage:
    python3 asc_to_geotiff.py /mnt/I2Rgus_Data/Chelsea_calibration/dem
    (every *.asc in the folder -> .tif next to it; or give files)
"""

import sys
import struct
import argparse
from pathlib import Path

import numpy as np

NODATA = -9999.0


def read_asc(path):
    with open(path) as f:
        head = {}
        for _ in range(6):
            k, v = f.readline().split()[:2]
            head[k.lower()] = float(v)
        data = np.loadtxt(f, dtype=np.float64, ndmin=2)
    nd = head.get("nodata_value", NODATA)
    data[data == nd] = np.nan
    return data, head


def write_geotiff(path, data, xll, yll, cell, epsg, description=""):
    """data: rows north to south (as in the ASCII grid). NaN -> nodata."""
    arr = np.where(np.isfinite(data), data, NODATA).astype("<f4")
    h, w = arr.shape
    top = yll + h * cell
    entries = []                                   # (tag, type, count, bytes)
    SHORT, LONG, DOUBLE, ASCII = 3, 4, 12, 2

    def add(tag, typ, values):
        fmt = {SHORT: "H", LONG: "I", DOUBLE: "d"}.get(typ)
        if typ == ASCII:
            raw = values.encode("ascii") + b"\0"
            entries.append((tag, typ, len(raw), raw))
        else:
            values = list(values)
            entries.append((tag, typ, len(values), struct.pack("<" + fmt * len(values), *values)))

    img = arr.tobytes()
    add(256, LONG, [w])
    add(257, LONG, [h])
    add(258, SHORT, [32])
    add(259, SHORT, [1])                           # no compression
    add(262, SHORT, [1])                           # min-is-black
    if description:
        add(270, ASCII, description)
    add(273, LONG, [0])                            # strip offset, patched below
    add(277, SHORT, [1])
    add(278, LONG, [h])
    add(279, LONG, [len(img)])
    add(284, SHORT, [1])
    add(339, SHORT, [3])                           # IEEE float
    add(33550, DOUBLE, [cell, cell, 0.0])          # ModelPixelScale
    add(33922, DOUBLE, [0.0, 0.0, 0.0, xll, top, 0.0])   # ModelTiepoint: UL corner
    add(34735, SHORT, [1, 1, 0, 4,
                       1024, 0, 1, 1,              # GTModelType = projected
                       1025, 0, 1, 1,              # GTRasterType = PixelIsArea
                       3072, 0, 1, int(epsg),      # ProjectedCSType
                       3076, 0, 1, 9001])          # ProjLinearUnits = metre
    add(42113, ASCII, f"{NODATA:g}")               # GDAL_NODATA
    entries.sort(key=lambda e: e[0])

    ifd_off = 8
    ifd_size = 2 + 12 * len(entries) + 4
    extra_off = ifd_off + ifd_size
    blobs, cursor = [], extra_off
    fields = []
    for tag, typ, count, raw in entries:
        if len(raw) <= 4:
            fields.append((tag, typ, count, raw.ljust(4, b"\0")))
        else:
            if cursor % 2:
                blobs.append(b"\0"); cursor += 1
            fields.append((tag, typ, count, struct.pack("<I", cursor)))
            blobs.append(raw); cursor += len(raw)
    if cursor % 4:
        pad = 4 - cursor % 4
        blobs.append(b"\0" * pad); cursor += pad
    img_off = cursor
    fields = [(t, ty, c, struct.pack("<I", img_off) if t == 273 else v) for t, ty, c, v in fields]

    with open(path, "wb") as f:
        f.write(b"II" + struct.pack("<HI", 42, ifd_off))
        f.write(struct.pack("<H", len(fields)))
        for t, ty, c, v in fields:
            f.write(struct.pack("<HHI", t, ty, c) + v)
        f.write(struct.pack("<I", 0))              # no further IFDs
        for b in blobs:
            f.write(b)
        f.write(img)


DESCRIPTIONS = {
    "_dem": "intertidal elevation, m NAVD88",
    "_spread": "16-84 percentile elevation range of repeat crossings, m",
    "_count": "frames (waterline crossings) per cell",
    "_source": "1 = measured, 2 = interpolated between waterlines",
    "_diff": "elevation change, m",
    "_sig": "significant elevation change, m",
    "_slope": "slope",
    "_resid": "RMS residual of the local plane fit, m",
    "_frames": "waterline frames in the fit neighbourhood",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("paths", nargs="+", help="folders (every *.asc in them) or .asc files")
    ap.add_argument("--epsg", type=int, default=32619,
                    help="projected CRS (default 32619 = WGS 84 / UTM 19N; 6348 = NAD83(2011) / UTM 19N)")
    args = ap.parse_args()

    files = []
    for p in map(Path, args.paths):
        files += sorted(p.glob("*.asc")) if p.is_dir() else [p]
    if not files:
        sys.exit("No .asc files found.")
    for asc in files:
        data, h = read_asc(asc)
        desc = next((d for k, d in DESCRIPTIONS.items() if asc.stem.endswith(k)), asc.stem)
        tif = asc.with_suffix(".tif")
        write_geotiff(tif, data, h["xllcorner"], h["yllcorner"], h["cellsize"], args.epsg,
                      f"{asc.stem}: {desc}; EPSG:{args.epsg}")
        print(f"wrote {tif}  ({data.shape[1]} x {data.shape[0]} cells at {h['cellsize']:g} m, "
              f"{int(np.isfinite(data).sum())} with data)")


if __name__ == "__main__":
    main()

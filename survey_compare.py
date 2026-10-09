#!/usr/bin/env python3
"""
DEM Against A Survey, With An Honest Label
=============================================
Compares an intertidal DEM (dem_from_contours.py) with a survey made
near the same time -- a lidar DSM, RTK check shots, or the ground
control points of a calibration -- and prints, next to every number,
how independent that comparison is.

WHY. A DEM that agrees with a survey proves something only if nothing
that built the DEM was tuned to that survey. The Jan 18-23 2025 DEM
once agreed with the 23 Jan 2025 lidar to a few centimetres -- with
camera pointings fitted to that same lidar (fit_eo_to_survey.py). That
is a check of the fit, not of the method. So every comparison carries a
label, chosen by the caller, printed exactly as given at the top of the
report and of the figure, with the reason:

  INDEPENDENT      nothing in the chain was fitted to or placed with
                   this survey
  CROSS-VALIDATED  fitted on other data from the same source, held out
                   here
  PARTLY-CIRCULAR  some step used this survey (e.g. the detection search
                   envelope was placed with it, or the calibration was
                   solved from these very points)
  CIRCULAR         pointing or offsets fitted to this survey

A comparison is never better than its weakest step. The script cannot
know the whole chain, but it checks what it can see and prints a
DOWNGRADE WARNING when the caller's label looks too good: a camera EO
file named or annotated as a survey fit (fit_eo_to_survey.py writes
*_lidar_EO.yaml with 'fitted to <survey>' in its notes), a calibration
dated the day of a point survey (solved from those GCPs?), a GCP target
file as the survey, a search envelope (--envelope-source, or a
provenance.json next to the DEM) placed with this survey. It also
prints the time between the photos and the survey: the beach moves, and
a difference over a week of storms is not all error.

HOW, DSM SURVEYS (lidar GeoTIFF or .asc, read by compare_dem_survey.py).
The 0.25 m lidar is aggregated to the DEM's own grid: each DEM cell gets
the MEDIAN of the survey inside it (8 x 8 sub-cells aligned with the DEM
grid, resampled bilinearly from the lidar's own grid, which is offset
from the DEM's: taking the lidar cells as they fall would bias every
cell up or down the slope by up to half a lidar cell), and only if the
survey covers at least --min-cover (default 50%) of the cell -- a cell
mostly on the lidar's water mask is not compared on the sliver that is
there. A part-covered cell at that edge is still the median of its dry
part, i.e. biased up the slope by up to slope x cell / 4 (2.5 cm on a
1:20 beach): such cells (50-99% covered) are counted, per elevation
band, and the statistics are also given for fully covered cells only.
The lidar's vertical datum is read from its GeoKeys (VerticalCSType);
the 2025 YSMP files state none, so NAVD88 is assumed and the report
says the geoid model is unknown. Then DEM - survey per cell, written as .asc and
GeoTIFF (metres of difference, not a NAVD88 height: no vertical datum
keys). A survey coarser than half a DEM cell is sampled bilinearly at
the cell centres instead. Lidar sees only the beach that was DRY at
flight time, so the report counts the DEM cells with survey data, and
says where the others are (by elevation band and by camera).

HOW, POINT SURVEYS. Three formats, recognised from the file:
  * Emlid Flow CSV (Name, Code, Easting, Northing, Elevation,
    Description, ..., CS name): every point is used; the 'CS name' column
    must say NAVD88 (a warning otherwise); statistics are also given per
    Description (Transect, High water wrack, ...) so wrack lines, which
    mark the swash limit and not a survey of the sand surface, are seen
    separately.
  * CIRN target files (no header: num,E,N,Z), as
    calibration/2025-11-13_Marconi_Extrinsic_Targets_c1_xyz.csv.
  * Any CSV with E/N/Z, x/y/z or Easting/Northing/Elevation headers.
The DEM value at each point is the bilinear interpolation of the finite
cells among the four surrounding cell centres, only where the cell that
contains the point has a value (so nothing is extrapolated beyond the
DEM's own cells). Points outside the DEM are counted and listed by name,
never silently dropped.

STATISTICS. Sign: DEM - survey, so positive = the DEM is too HIGH. n,
median, mean, NMAD (1.4826 x median absolute deviation: a robust
standard deviation), RMSE, p5/p95, and the share within +/-0.10, 0.20,
0.50 m: overall, by 0.5 m band of SURVEY elevation (the reference, not
the thing tested), by camera, by distance from the camera, and with the
DEM's own _spread/_count grids, by the DEM's repeatability. The camera
of a cell or point: with --camera-eo, the photo that sees it nearer its
centre (projected through the calibration); else with --contours, the
camera whose waterlines crossed that DEM cell most; else the live
station's seam at northing 4638415 (c1 to the south, c2 to the north;
right for the 2025-11-13 view only, and the report says so).

WATERLINES (--contours). The georectified waterlines are also compared
with the survey directly:
  * on a DSM, the survey is sampled bilinearly at every waterline point
    (the logic of compare_dem_survey.py);
  * on RTK transects (Emlid points described 'Transect'), the logic of
    compare_rtk.py: consecutive transect points (a jump of more than
    --transect-split m starts a new transect) make a profile; a
    waterline point counts if it lies within --transect-tolerance m of
    the line and between its first and last surveyed point (never
    extrapolated), where the RTK elevation is interpolated ALONG the
    profile; each frame with >= --transect-min-points such points gives
    ONE value, the median of its waterline - RTK. Matching each line to
    the nearest shot instead would compare points up to a metre apart on
    a 1:7 face (~0.14 m of slope) and depend on which few lines pass
    near a shot. Other points of such a file (wrack lines mark the swash
    limit, not the sand) are not used for the waterlines;
  * on isolated points (GCPs, a CSV without transects), the nearest
    point within --tolerance m, one value per frame and point, NOT
    slope-corrected (the report gives the median distance).
The waterline's elevation is beach_elevation_navd88 (still water + wave
setup) where present, else tide_elevation_navd88 (still water only,
which puts a swash-marked line LOW by about the setup); the report says
which, per row count. Per day and camera: frames, values, median, NMAD.

FEW VALUES. NMAD, p5 and p95 are not given for fewer than 3 values (the
NMAD of one value is 0, which would read as perfect agreement), and any
group with fewer than 10 is marked as not an estimate.

OUTPUT (in --output-dir):
  NAME_dem_minus_survey.asc / .tif   (dsm only) the difference grid
  NAME_points.csv     (points only) every point, the DEM value, the
                      difference, and why a point was not compared
  NAME_comparison.txt the label, warnings, headline and tables
  NAME_comparison.csv one row per statistic group
  NAME_comparison.json the headline (for survey_products.py --summary)
  NAME_comparison.png one page: map, histogram, median by elevation
                      band, headline and label

Runs in 2-5 s here (2 s for the Jan 2025 lidar, 4 s with 390k waterline
points); expect ~5-15 s on the station NUC, + ~3 s per 100k waterline
points with --contours.

Usage:
    python3 survey_compare.py --dem dem/2025-01-23_dem.asc \\
        --survey 2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif --survey-type dsm \\
        --name 2025-01-23_jan_lidar --output-dir compare \\
        --label INDEPENDENT --why "GCP calibration of 23 Jan; envelope placed with the Mar lidar" \\
        --survey-date 2025-01-23 --photo-dates 2025-01-18 2025-01-23 \\
        --contours waterlines/contour_points_ground_filtered.csv \\
        --camera-eo c1=calibration/CACO03_c1_20250123_EO.yaml c2=calibration/CACO03_c2_20250123_EO.yaml

    python3 survey_compare.py --dem dem_intertidal_7day_dem.asc \\
        --survey 2026-09-29_Marconi_Checkshots.csv --survey-type points \\
        --name 2026-09-27_rtk --label INDEPENDENT --why "RTK check shots, not used anywhere" \\
        --survey-date 2026-09-29 --photo-dates 2026-09-20 2026-10-04
"""

import re
import sys
import csv
import json
import time
import argparse
import textwrap
import warnings
from datetime import date
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from asc_to_geotiff import read_asc, write_geotiff        # noqa: E402
from compare_dem_survey import read_survey, sample        # noqa: E402

LABELS = ("INDEPENDENT", "CROSS-VALIDATED", "PARTLY-CIRCULAR", "CIRCULAR")   # best -> worst
# The tower: mean of the c1/c2 positions of the live calibration
# (CACO05_<cam>_20251113_EO-CV.yaml). Earlier setups stood within ~20 m,
# which does not matter for 50-100 m distance bins.
TOWER_EN = (420085.05, 4638322.39)
# Live view (2025-11-13 calibration): c1 sees the beach south of this
# northing, c2 north of it.
SEAM_NORTHING = 4638415.0
DISTANCE_BINS = (0, 100, 150, 200, 250, 300, 400, 600, 1e9)
BAND = 0.5                       # elevation band width, m
NOMINAL_EPSG = 32619             # as asc_to_geotiff.py writes the DEM
MIN_BAND_N = 10                  # fewer values than this in a group: flagged as few
MIN_SPREAD_N = 3                 # fewer than this: no NMAD / p5 / p95 (one value has NMAD 0)
TIME_GAP_WARN_DAYS = 3
# Waterlines on RTK transects: compare_rtk.py's defaults.
TRANSECT_TOLERANCE, TRANSECT_SPLIT, TRANSECT_MIN_POINTS = 2.0, 15.0, 3
FULL_COVER = 0.999               # a DEM cell the survey covers entirely (8 x 8 sub-cells)
VERTICAL_UNKNOWN = ("vertical datum not stated in the survey file; assumed NAVD88 (geoid model "
                    "unknown; GEOID12B and GEOID18 differ by a few cm here)")

# Ink and surfaces as dem_figure.py: neutral, so the only colour is data.
INK = "#1f1e1c"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e6e5df"
AXIS = "#c3c2b7"
BLANK_FACE = "#f2f1ed"
BLANK_HATCH = "#9a9992"
# Diverging blue <-> red with a neutral gray at zero; blue = DEM low, red = DEM high.
DIVERGING = ["#0d366b", "#2a78d6", "#9ec5f4", "#e4e3df", "#f2b6ae", "#e34948", "#8a1c1c"]
# Elevation: one hue (sand), light = high, dark = low (as dem_figure.py).
SAND = ["#f7e5cb", "#e0c399", "#c5a36e", "#a7844d", "#876834", "#674d21", "#483413"]
PAGE_RC = {"font.size": 10, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
           "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
           "hatch.linewidth": 0.6, "axes.linewidth": 0.8}


# ---------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------

def dem_stem(dem_path):
    s = str(Path(dem_path).with_suffix(""))
    return s[:-4] if s.endswith("_dem") else s


def read_grid(path):
    """ESRI ASCII grid -> (rows north to south, xll, y_top, cell, header)."""
    g, h = read_asc(path)
    return g, h["xllcorner"], h["yllcorner"] + g.shape[0] * h["cellsize"], h["cellsize"], h


def same_grid(h1, h2):
    return all(abs(h1[k] - h2[k]) < 1e-6 for k in ("ncols", "nrows", "xllcorner", "yllcorner", "cellsize"))


def survey_crs_text(path):
    """The coordinate system the survey file states, as text ('' if none)."""
    p = Path(path)
    if p.suffix.lower() == ".asc":
        prj = p.with_suffix(".prj")
        if prj.exists():
            m = re.search(r'(?:PROJCS|PROJCRS)\["([^"]+)"', prj.read_text(errors="replace"))
            return m.group(1) if m else prj.read_text(errors="replace")[:80]
        return ""
    # GeoTIFF: the GeoAsciiParams text (e.g. 'NAD83(2011) / UTM zone 19N|NAD83(2011)|')
    raw = p.read_bytes()
    m = re.search(rb"(?:NAD83|WGS ?84|ETRS89)[^|\x00]{0,40}?UTM [Zz]one ?\d+ ?[NS]?", raw) or \
        re.search(rb"UTM [Zz]one ?\d+ ?[NS]?", raw)
    return m.group(0).decode("ascii", "replace").strip() if m else ""


def tiff_geokeys(path):
    """{GeoKey id: value} of a GeoTIFF / BigTIFF's short-valued GeoKeys, from the first IFD only
    (the raster is not read). None if the file cannot be parsed."""
    import struct
    try:
        with open(path, "rb") as f:
            hdr = f.read(16)
            bo = {b"II": "<", b"MM": ">"}.get(hdr[:2])
            if not bo:
                return None
            ver = struct.unpack(bo + "H", hdr[2:4])[0]
            if ver == 42:
                off, head, esize, cf = struct.unpack(bo + "I", hdr[4:8])[0], 2, 12, "I"
            elif ver == 43:
                off, head, esize, cf = struct.unpack(bo + "Q", hdr[8:16])[0], 8, 20, "Q"
            else:
                return None
            f.seek(off)
            n = struct.unpack(bo + ("H" if ver == 42 else "Q"), f.read(head))[0]
            ents = f.read(esize * n)
            cs = struct.calcsize(cf)
            for i in range(n):
                e = ents[esize * i: esize * (i + 1)]
                tag, typ = struct.unpack(bo + "HH", e[:4])
                if tag != 34735 or typ != 3:
                    continue
                cnt = struct.unpack(bo + cf, e[4:4 + cs])[0]
                nb = 2 * cnt
                if nb <= cs:
                    data = e[4 + cs:4 + cs + nb]
                else:
                    f.seek(struct.unpack(bo + cf, e[4 + cs:4 + 2 * cs])[0])
                    data = f.read(nb)
                k = struct.unpack(bo + "H" * cnt, data)
                return {k[j]: k[j + 3] for j in range(4, len(k) - 3, 4) if k[j + 1] == 0}
        return {}
    except (OSError, struct.error, ValueError):
        return None


def survey_vertical_text(path):
    """-> (what the DSM file says of its vertical datum, or '' if nothing; a note when it says
    nothing). GeoTIFF: VerticalCSType (GeoKey 4096; 5703 = NAVD88). .asc: its .prj."""
    p = Path(path)
    if p.suffix.lower() == ".asc":
        prj = p.with_suffix(".prj")
        t = prj.read_text(errors="replace") if prj.exists() else ""
        m = re.search(r'VERTCS\["([^"]+)"', t) or re.search(r"(NAVD ?88[^\"\],]*)", t)
        return (m.group(1), "") if m else ("", VERTICAL_UNKNOWN)
    keys = tiff_geokeys(p)
    if keys is None:
        return "", "vertical datum not checked (GeoKeys unreadable); assumed NAVD88"
    v = keys.get(4096)
    if v is None:
        return "", VERTICAL_UNKNOWN
    if v == 5703:
        return "EPSG:5703 (NAVD88 height)", ""
    return f"EPSG:{v}", (f"the survey states vertical EPSG:{v}, not NAVD88 (5703): the difference includes "
                         f"the datum difference")


def _num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def read_points(path):
    """
    Point survey -> dict(E, N, Z, name, code, desc, format, crs, skipped).
    Emlid Flow CSV, CIRN 'num,E,N,Z' (no header) or a CSV with E/N/Z,
    x/y/z or Easting/Northing/Elevation headers.
    """
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows = [r for r in csv.reader(f) if any(c.strip() for c in r)]
    if not rows:
        sys.exit(f"{path}: empty")
    out = {"E": [], "N": [], "Z": [], "name": [], "code": [], "desc": [], "crs": [],
           "skipped": 0}
    first = rows[0]
    if len(first) >= 4 and all(_num(c) is not None for c in first[:4]):
        out["format"] = "CIRN targets (num,E,N,Z; no header)"
        for r in rows:
            v = [_num(c) for c in r[:4]]
            if len(r) < 4 or any(x is None for x in v[1:4]):
                out["skipped"] += 1
                continue
            out["name"].append(r[0].strip()); out["E"].append(v[1]); out["N"].append(v[2])
            out["Z"].append(v[3]); out["code"].append(""); out["desc"].append("target")
            out["crs"].append("")
    else:
        head = [c.strip() for c in first]
        low = [c.lower() for c in head]

        def col(*names):
            for nm in names:
                if nm in low:
                    return low.index(nm)
            return None
        ie = col("easting", "e", "x", "east", "easting_utm19")
        inn = col("northing", "n", "y", "north", "northing_utm19")
        iz = col("elevation", "z", "elev", "height", "elevation_navd88", "z_navd88", "navd88")
        if ie is None or inn is None or iz is None:
            sys.exit(f"{path}: no Easting/Northing/Elevation (E/N/Z, x/y/z) columns in {head}")
        iname = col("name", "id", "num", "point", "point_id")
        icode = col("code")
        idesc = col("description", "desc", "type", "feature")
        icrs = col("cs name", "crs", "cs_name")
        emlid = "cs name" in low and "elevation" in low
        out["format"] = "Emlid Flow CSV" if emlid else "CSV with headers " + \
            ", ".join(head[i] for i in (ie, inn, iz))
        for k, r in enumerate(rows[1:]):
            g = lambda i: r[i].strip() if (i is not None and i < len(r)) else ""   # noqa: E731
            v = (_num(g(ie)), _num(g(inn)), _num(g(iz)))
            if any(x is None for x in v):
                out["skipped"] += 1
                continue
            out["E"].append(v[0]); out["N"].append(v[1]); out["Z"].append(v[2])
            out["name"].append(g(iname) or str(k + 1)); out["code"].append(g(icode))
            out["desc"].append(g(idesc) or "point"); out["crs"].append(g(icrs))
    for k in ("E", "N", "Z"):
        out[k] = np.array(out[k], float)
    if not len(out["E"]):
        sys.exit(f"{path}: no point with numeric coordinates")
    return out


def read_contours(path, first=None, last=None):
    """Georectified waterline points -> dict of arrays; rows without ground coordinates skipped."""
    E, N, Z, cam, day, frame, used_beach = [], [], [], [], [], [], []
    with open(path, newline="") as f:
        rd = csv.reader(f)
        head = next(rd)
        ix = {c: i for i, c in enumerate(head)}
        if "easting_utm19" not in ix:
            sys.exit(f"{path}: no easting_utm19 column (run georectify.py first)")
        ie, inn = ix["easting_utm19"], ix["northing_utm19"]
        it, ib = ix.get("tide_elevation_navd88"), ix.get("beach_elevation_navd88")
        ic, itime, isrc = ix.get("camera"), ix.get("capture_time_utc"), ix.get("source_file")
        for r in rd:
            if len(r) <= max(ie, inn) or not r[ie]:
                continue
            d = r[itime][:10] if itime is not None else ""
            if (first and d < first) or (last and d > last):
                continue
            zb = _num(r[ib]) if ib is not None and ib < len(r) else None
            zt = _num(r[it]) if it is not None and it < len(r) else None
            z = zb if zb is not None else zt
            if z is None:
                continue
            E.append(float(r[ie])); N.append(float(r[inn])); Z.append(z)
            used_beach.append(zb is not None)
            cam.append(r[ic] if ic is not None else "?"); day.append(d)
            frame.append(r[isrc] if isrc is not None else d)
    return {"E": np.array(E), "N": np.array(N), "Z": np.array(Z), "cam": np.array(cam),
            "day": np.array(day), "frame": np.array(frame), "beach": np.array(used_beach, bool),
            "has_beach_column": ib is not None}


# ---------------------------------------------------------------------
# Cameras
# ---------------------------------------------------------------------

def parse_cam_specs(specs, what):
    out = {}
    for s in specs or []:
        if "=" not in s:
            sys.exit(f"{what}: expected cam=PATH, got '{s}'")
        k, v = s.split("=", 1)
        out[k.strip()] = Path(v.strip())
    return out


def find_io(eo_path, cam):
    """The lens file for an EO: the same setup's IO next to it or in calibration/; else the
    CACO05 IO (the lens has not changed since 2024-08: CACO03 IO = CACO05 IO)."""
    setup = Path(eo_path).name.split("_")[0]
    for d in (Path(eo_path).parent, HERE / "calibration"):
        c = sorted(d.glob(f"{setup}_{cam}_*_IO.yaml"))
        if c:
            return c[-1]
    c = sorted((HERE / "calibration").glob(f"CACO05_{cam}_*_IO.yaml"))
    return c[-1] if c else None


def load_cameras(eo_specs, io_specs):
    from georectify import load_extrinsics, load_intrinsics
    cams = {}
    for cam, eo_p in sorted(eo_specs.items()):
        if not eo_p.exists():
            sys.exit(f"--camera-eo {cam}: {eo_p} not found")
        io_p = io_specs.get(cam) or find_io(eo_p, cam)
        io = None
        if io_p is not None and Path(io_p).exists():
            io = load_intrinsics(io_p)
        cams[cam] = {"eo_path": eo_p, "io_path": io_p, "eo": load_extrinsics(eo_p), "io": io}
    return cams


def history_note(eo_path):
    """The calibration_history.csv note for this EO file ('' if none)."""
    for d in (Path(eo_path).parent, HERE / "calibration"):
        p = d / "calibration_history.csv"
        if p.exists():
            with open(p, newline="") as f:
                for r in csv.DictReader(f):
                    if r.get("file") == Path(eo_path).name:
                        return (r.get("note") or "").strip()
    return ""


def cell_majority_camera(cont, xll, ytop, cell, nrows, ncols):
    """Grid (rows north to south) of the camera whose waterlines crossed each cell most."""
    names = sorted(set(cont["cam"].tolist()))
    col = np.floor((cont["E"] - xll) / cell).astype(int)
    row = np.floor((ytop - cont["N"]) / cell).astype(int)
    ok = (col >= 0) & (col < ncols) & (row >= 0) & (row < nrows)
    gid = row[ok] * ncols + col[ok]
    counts = np.zeros((len(names), nrows * ncols))
    for k, nm in enumerate(names):
        counts[k] = np.bincount(gid[cont["cam"][ok] == nm], minlength=nrows * ncols)
    best = np.argmax(counts, axis=0)
    out = np.array(names, dtype=object)[best]
    out[counts.max(axis=0) == 0] = "no waterline"
    return out.reshape(nrows, ncols)


def assign_camera(E, N, Z, cams=None, cell_cam=None, grid=None):
    """
    -> (camera name per position, how it was decided).
    cams: geometric -- the photo that sees the point nearer its centre column.
    cell_cam + grid (xll, ytop, cell): the camera whose waterlines crossed that DEM cell most.
    Otherwise the live station's seam.
    """
    E, N = np.asarray(E, float), np.asarray(N, float)
    if cams and all(c["io"] is not None for c in cams.values()):
        from view_reproject import ground_to_pixel
        Zs = np.where(np.isfinite(Z), Z, 0.0)
        best = np.full(E.shape, np.inf)
        out = np.full(E.shape, "neither", dtype=object)
        for cam, c in sorted(cams.items()):
            U, _, ok = ground_to_pixel(E, N, Zs, c["io"], c["eo"])
            off = np.abs(U / c["io"][0] - 0.5)
            better = ok & (off < best)
            out[better] = cam
            best[better] = off[better]
        return out, ("calibration footprints (" + ", ".join(
            f"{k}={Path(c['eo_path']).name}" for k, c in sorted(cams.items()))
            + "; where both photos see a spot, the one nearer its centre)")
    if cell_cam is not None:
        xll, ytop, cell = grid
        nrows, ncols = cell_cam.shape
        col = np.floor((E - xll) / cell).astype(int)
        row = np.floor((ytop - N) / cell).astype(int)
        ok = (col >= 0) & (col < ncols) & (row >= 0) & (row < nrows)
        out = np.full(E.shape, "no waterline", dtype=object)
        out[ok] = cell_cam[row[ok], col[ok]]
        return out, "the camera whose waterlines crossed the DEM cell most (--contours)"
    out = np.where(N < SEAM_NORTHING, "c1", "c2").astype(object)
    return out, (f"seam at northing {SEAM_NORTHING:.0f} (c1 south, c2 north): right for the "
                 f"2025-11-13 view only; give --camera-eo for other dates")


# ---------------------------------------------------------------------
# Survey on the DEM grid; DEM at points
# ---------------------------------------------------------------------

def bilinear(grid, x0, y_top, cell, E, N):
    """
    Bilinear value of a grid (rows north to south) at E, N. NaN outside the grid, or where a
    neighbour that carries weight is nodata (never interpolates across a gap). Within half a
    cell of the grid's outer edge the edge value is used.
    """
    nr, nc = grid.shape
    fc = (np.asarray(E, float) - x0) / cell - 0.5
    fr = (y_top - np.asarray(N, float)) / cell - 0.5
    inside = (fc >= -0.5) & (fc <= nc - 0.5) & (fr >= -0.5) & (fr <= nr - 0.5)
    fc = np.clip(fc, 0, nc - 1)
    fr = np.clip(fr, 0, nr - 1)
    c0 = np.clip(np.floor(fc).astype(int), 0, max(nc - 2, 0))
    r0 = np.clip(np.floor(fr).astype(int), 0, max(nr - 2, 0))
    tc, tr = fc - c0, fr - r0
    c1, r1 = np.minimum(c0 + 1, nc - 1), np.minimum(r0 + 1, nr - 1)
    acc = np.zeros(fc.shape)
    bad = ~inside
    for r, c, w in ((r0, c0, (1 - tc) * (1 - tr)), (r0, c1, tc * (1 - tr)),
                    (r1, c0, (1 - tc) * tr), (r1, c1, tc * tr)):
        z = grid[r, c]
        use = w > 0
        bad |= use & ~np.isfinite(z)
        acc += np.where(use & np.isfinite(z), w * np.nan_to_num(z), 0.0)
    return np.where(bad, np.nan, acc)


def block_median(survey, sx0, sy_top, scell, gx0, gy_top, cell, nrows, ncols, chunk=16):
    """
    Median of the survey inside each grid cell, and the share of the cell it covers.
    Rows north to south. The survey is first resampled (bilinearly) onto k x k sub-cells
    ALIGNED with the grid, k = cell / scell (8 for 0.25 m lidar in 2 m cells): the lidar's
    own grid is offset from the DEM's (by 0.676 m E for the 2025 YSMP DSMs), and taking its
    cells as they fall puts the sample's centre up to half a lidar cell off the DEM cell's
    centre -- 4 mm on the synthetic 1:20 beach, ~1 cm on a 1:10 one, the same sign in every
    cell. Aligned sub-cells are symmetric about the centre, so on a plane the median is exact.
    """
    k = max(1, int(round(cell / scell)))
    sub = (np.arange(k) + 0.5) * cell / k
    xs = (gx0 + np.arange(ncols)[:, None] * cell + sub[None, :]).ravel()
    med = np.full((nrows, ncols), np.nan)
    cover = np.zeros((nrows, ncols))
    for r0 in range(0, nrows, chunk):
        r1 = min(r0 + chunk, nrows)
        ys = (gy_top - np.arange(r0, r1)[:, None] * cell - sub[None, :]).ravel()
        X, Y = np.meshgrid(xs, ys)
        v = bilinear(survey, sx0, sy_top, scell, X.ravel(), Y.ravel())
        v = v.reshape(r1 - r0, k, ncols, k).transpose(0, 2, 1, 3).reshape(r1 - r0, ncols, k * k)
        n = np.isfinite(v).sum(axis=2)
        cover[r0:r1] = n / float(k * k)
        if n.any():
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)       # all-NaN cells -> NaN
                med[r0:r1] = np.nanmedian(v, axis=2)
    return med, cover


def survey_on_grid(survey, sx0, sy_top, scell, gx0, gy_top, cell, nrows, ncols):
    """-> (survey value per grid cell, covered share, method text)."""
    if scell <= cell / 2:
        med, cover = block_median(survey, sx0, sy_top, scell, gx0, gy_top, cell, nrows, ncols)
        k = max(1, int(round(cell / scell)))
        return med, cover, (f"median of the {scell:g} m survey inside each {cell:g} m DEM cell "
                            f"({k} x {k} sub-cells aligned with the DEM grid, resampled "
                            f"bilinearly from the survey's own grid)")
    xs = gx0 + (np.arange(ncols) + 0.5) * cell
    ys = gy_top - (np.arange(nrows) + 0.5) * cell
    X, Y = np.meshgrid(xs, ys)
    s = bilinear(survey, sx0, sy_top, scell, X.ravel(), Y.ravel()).reshape(nrows, ncols)
    return s, np.isfinite(s).astype(float), (f"survey ({scell:g} m) coarser than half a DEM "
                                             f"cell: sampled bilinearly at the cell centres")


def dem_at_points(dem, xll, ytop, cell, E, N):
    """
    Bilinear DEM value at each point over the FINITE cells among the four surrounding cell
    centres, only where the cell containing the point has a value. -> (value, status).
    """
    nrows, ncols = dem.shape
    ci = np.floor((E - xll) / cell).astype(int)
    ri = np.floor((ytop - N) / cell).astype(int)
    inside = (ci >= 0) & (ci < ncols) & (ri >= 0) & (ri < nrows)
    has = np.zeros(E.shape, bool)
    has[inside] = np.isfinite(dem[ri[inside], ci[inside]])
    fc = (E - xll) / cell - 0.5
    fr = (ytop - N) / cell - 0.5
    c0, r0 = np.floor(fc).astype(int), np.floor(fr).astype(int)
    tc, tr = fc - c0, fr - r0
    acc = np.zeros(E.shape)
    wsum = np.zeros(E.shape)
    for dc, dr, w in ((0, 0, (1 - tc) * (1 - tr)), (1, 0, tc * (1 - tr)),
                      (0, 1, (1 - tc) * tr), (1, 1, tc * tr)):
        c, r = c0 + dc, r0 + dr
        ok = (c >= 0) & (c < ncols) & (r >= 0) & (r < nrows)
        z = np.full(E.shape, np.nan)
        z[ok] = dem[r[ok], c[ok]]
        ok &= np.isfinite(z)
        acc[ok] += w[ok] * z[ok]
        wsum[ok] += w[ok]
    val = np.full(E.shape, np.nan)
    good = has & (wsum > 0)
    val[good] = acc[good] / wsum[good]
    status = np.where(good, "compared", np.where(inside, "dem_cell_without_value",
                                                  "outside_dem_extent")).astype(object)
    return val, status


def nearest_cell_distance(dem, xll, ytop, cell, E, N):
    """Distance (m) from each point to the nearest DEM cell centre with a value."""
    rr, cc = np.nonzero(np.isfinite(dem))
    if not len(rr):
        return np.full(E.shape, np.nan)
    cx = xll + (cc + 0.5) * cell
    cy = ytop - (rr + 0.5) * cell
    out = np.empty(E.shape)
    for k in range(0, len(E), 256):
        d = np.hypot(E[k:k + 256, None] - cx[None, :], N[k:k + 256, None] - cy[None, :])
        out[k:k + 256] = d.min(axis=1)
    return out


# ---------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------

STAT_KEYS = ("n", "median", "mean", "nmad", "rmse", "p5", "p95", "within_0.10", "within_0.20",
             "within_0.50")


def stats(d):
    """n, median, mean, NMAD, RMSE, p5, p95, shares within 0.1/0.2/0.5 m. NMAD, p5 and p95 are
    NaN below MIN_SPREAD_N values: the NMAD of one value is 0, which would read as perfect."""
    d = np.asarray(d, float)
    d = d[np.isfinite(d)]
    if not len(d):
        return dict(zip(STAT_KEYS, [0] + [float("nan")] * 9))
    med = float(np.median(d))
    a = np.abs(d)
    spread = len(d) >= MIN_SPREAD_N
    nan = float("nan")
    return {"n": int(len(d)), "median": med, "mean": float(d.mean()),
            "nmad": float(1.4826 * np.median(np.abs(d - med))) if spread else nan,
            "rmse": float(np.sqrt(np.mean(d ** 2))),
            "p5": float(np.percentile(d, 5)) if spread else nan,
            "p95": float(np.percentile(d, 95)) if spread else nan,
            # (+0.1 mm: the grids are written to 0.1 mm, so exactly 0.20 m counts as within 0.20)
            "within_0.10": float(np.mean(a <= 0.10 + 1e-4)),
            "within_0.20": float(np.mean(a <= 0.20 + 1e-4)),
            "within_0.50": float(np.mean(a <= 0.50 + 1e-4))}


def band_name(lo):
    return f"{lo:+.1f} to {lo + BAND:+.1f}"


def group_rows(d, groups):
    """groups: list of (group_type, group name, mask) -> list of stat rows (n >= 1 only)."""
    rows = []
    for gt, gname, m in groups:
        s = stats(d[m])
        if s["n"]:
            rows.append(dict({"group_type": gt, "group": gname}, **s))
    return rows


def standard_groups(d, zref, cam, dist, extra=()):
    """Overall, elevation band (of zref), camera, distance -- plus any extra groups."""
    ok = np.isfinite(d)
    g = [("overall", "all", ok)]
    zb = np.floor(zref / BAND) * BAND
    for lo in sorted(set(zb[ok & np.isfinite(zb)].tolist())):
        g.append(("elevation_band", band_name(lo), ok & (zb == lo)))
    for c in sorted(set(cam[ok].tolist())):
        g.append(("camera", c, ok & (cam == c)))
    for lo, hi in zip(DISTANCE_BINS, DISTANCE_BINS[1:]):
        m = ok & (dist >= lo) & (dist < hi)
        if m.any():
            g.append(("distance", f"{lo:.0f}-{hi:.0f} m" if hi < 1e8 else f">{lo:.0f} m", m))
    return group_rows(d, g + list(extra))


def repeatability_groups(d, spread, count):
    g = []
    if spread is not None:
        for lo, hi in ((0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 99)):
            m = np.isfinite(d) & (spread >= lo) & (spread < hi)
            g.append(("dem_spread", f"{lo:.1f}-{hi:.1f} m" if hi < 50 else f">{lo:.1f} m", m))
    if count is not None:
        for lo, hi in ((0, 6), (6, 11), (11, 21), (21, 51), (51, 1e9)):
            m = np.isfinite(d) & (count >= lo) & (count < hi)
            g.append(("dem_frames", f"{lo:.0f}-{hi - 1:.0f}" if hi < 1e8 else f">{lo - 1:.0f}", m))
    return g


# ---------------------------------------------------------------------
# The label
# ---------------------------------------------------------------------

def _rank(label):
    return LABELS.index(label) if label in LABELS else 0


def _date_from_name(name):
    m = re.search(r"(20\d\d)-?(\d\d)-?(\d\d)", name)
    return f"{m.group(1)}-{m.group(2)}-{m.group(3)}" if m else None


def _walk(obj, path=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            for x in _walk(v, f"{path}.{k}" if path else str(k)):
                yield x
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            for x in _walk(v, f"{path}[{i}]"):
                yield x
    elif isinstance(obj, str):
        yield path, obj


def audit_label(label, survey_path, survey_type, survey_date, cams, envelope_source, dem_path):
    """
    What the script can see of the chain. -> list of (suggested label or None, text):
    a suggested label worse than `label` is a DOWNGRADE WARNING; None is a note.
    """
    found = []
    sname, sstem = Path(survey_path).name, Path(survey_path).stem
    gcp_file = survey_type == "points" and re.search(r"target|gcp", sname, re.I)
    gcp_date = _date_from_name(sname) if gcp_file else None
    for cam, c in sorted((cams or {}).items()):
        p = Path(c["eo_path"])
        text = p.read_text(errors="replace")
        notes = " ".join(ln.lstrip("#").strip() for ln in text.splitlines() if ln.strip().startswith("#"))
        fitted = re.findall(r"fitted to (\S+)", notes)
        survey_fit = ("_lidar_EO" in p.name or "_survey" in p.name or "fit_eo_to_survey" in notes)
        if sname in notes or (len(sstem) > 6 and sstem in notes):
            found.append(("CIRCULAR", f"{cam} pointing {p.name} was fitted to this survey "
                                      f"(its notes: '{textwrap.shorten(notes, 170)}')"))
        elif survey_fit and fitted:
            found.append((None, f"{cam} pointing {p.name} was fitted to another survey "
                                f"({', '.join(fitted)}): not circular against this one, but say so"))
        elif survey_fit:
            found.append(("CIRCULAR", f"{cam} pointing {p.name} is a survey fit "
                                      f"(fit_eo_to_survey.py naming) and does not say which survey: "
                                      f"if this one, the comparison is CIRCULAR"))
        eo_date = _date_from_name(p.name)
        hist = history_note(p)
        if survey_type == "points" and eo_date and (eo_date == survey_date or eo_date == gcp_date):
            found.append(("PARTLY-CIRCULAR", f"{cam} calibration {p.name} is dated the day of these "
                          f"points ({eo_date}): if it was solved from them, the pointing is fitted to "
                          f"them (elevations still come from water levels) -> PARTLY-CIRCULAR"
                          + (f". History: '{hist}'" if hist else "")))
        elif survey_type == "dsm" and eo_date and eo_date == survey_date:
            found.append((None, f"{cam} calibration {p.name} is dated the survey day; check it was "
                                f"not fitted to this survey"
                                + (f". History: '{hist}'" if hist else " (no history note)")))
    if gcp_file and not cams:
        found.append(("PARTLY-CIRCULAR", f"{sname} looks like calibration targets (GCPs): if the "
                      f"cameras' calibration was solved from these points, the comparison is "
                      f"PARTLY-CIRCULAR (no --camera-eo to check)"))
    elif gcp_file and not any(s == "PARTLY-CIRCULAR" for s, _ in found):
        found.append((None, f"{sname} looks like calibration targets (GCPs), but no camera EO "
                            f"given is dated {gcp_date or 'their day'}: not the calibration's own points"))
    if envelope_source:
        en = Path(envelope_source)
        if en.name == sname or en.stem == sstem:
            found.append(("PARTLY-CIRCULAR", f"the detection search envelope was placed with this "
                                             f"survey ({en.name})"))
        else:
            found.append((None, f"search envelope placed with {en.name}, not this survey"))
    # provenance.json next to the DEM (survey_products.py): any step that names this survey
    for d in (Path(dem_path).resolve().parent, Path(dem_path).resolve().parent.parent):
        pj = d / "provenance.json"
        if not pj.exists():
            continue
        try:
            prov = json.loads(pj.read_text())
        except (OSError, ValueError):
            found.append((None, f"{pj} unreadable: not checked"))
            break
        for key, val in _walk(prov):
            k = key.lower()
            if "compar" in k or not (sname in val or (len(sstem) > 6 and sstem in val)):
                continue
            if "envelope" in k:
                found.append(("PARTLY-CIRCULAR", f"provenance {pj.name}: {key} = {val}"))
            elif re.search(r"calib|pointing|fit|(^|[._\[])eo([._\]]|$)", k):
                found.append(("CIRCULAR", f"provenance {pj.name}: {key} = {val}"))
            else:
                found.append((None, f"provenance {pj.name} names this survey at {key}"))
        if not envelope_source:
            env = [(k, v) for k, v in _walk(prov) if "envelope" in k.lower()]
            if env:
                found.append((None, "provenance search envelope: " + "; ".join(
                    f"{k} = {v}" for k, v in env[:3])))
        break
    return found


def time_gap(survey_date, photo_dates):
    """-> (days outside the photo window or None, days from its middle or None, text)."""
    if not survey_date:
        return None, None, "survey date not given: the time gap is not known"
    if not photo_dates:
        return None, None, f"survey {survey_date}; photo dates not given: the time gap is not known"
    s = date.fromisoformat(survey_date)
    a, b = date.fromisoformat(photo_dates[0]), date.fromisoformat(photo_dates[1])
    mid = (a.toordinal() + b.toordinal()) / 2.0
    from_mid = s.toordinal() - mid
    outside = 0 if a <= s <= b else ((a - s).days if s < a else (s - b).days)
    span = (b - a).days + 1
    where = ("inside the photo window" if outside == 0 else
             f"{outside} day(s) {'before' if s < a else 'after'} the photo window")
    text = (f"survey {survey_date}; photos {a} to {b} ({span} days): {where}, "
            f"{abs(from_mid):.1f} days {'after' if from_mid > 0 else 'before'} its middle"
            if from_mid else
            f"survey {survey_date}; photos {a} to {b} ({span} days): {where}, at its middle")
    return outside, from_mid, text


# ---------------------------------------------------------------------
# Waterlines against the survey
# ---------------------------------------------------------------------

def waterlines_vs_survey(cont, survey_type, grid=None, pts=None, tolerance=1.0, dist_out=None):
    """waterline - survey at each waterline point (NaN where the survey has nothing). DSM: the
    survey sampled at the point. Points: the nearest survey point within `tolerance` m, one value
    per frame and survey point, NOT slope-corrected (dist_out, a list, gets the distances)."""
    if survey_type == "dsm":
        g, x0, y0, c = grid
        zs = sample(g, x0, y0, c, cont["E"], cont["N"])
        return cont["Z"] - zs
    d = np.full(cont["E"].shape, np.nan)
    match = np.full(cont["E"].shape, -1)
    dist = np.full(cont["E"].shape, np.nan)
    if not len(pts["E"]):
        return d
    lo_e, hi_e = pts["E"].min() - tolerance, pts["E"].max() + tolerance
    lo_n, hi_n = pts["N"].min() - tolerance, pts["N"].max() + tolerance
    near = np.nonzero((cont["E"] >= lo_e) & (cont["E"] <= hi_e) &
                      (cont["N"] >= lo_n) & (cont["N"] <= hi_n))[0]
    for k in range(0, len(near), 4096):
        idx = near[k:k + 4096]
        dd = np.hypot(cont["E"][idx, None] - pts["E"][None, :], cont["N"][idx, None] - pts["N"][None, :])
        j = dd.argmin(axis=1)
        dmin = dd[np.arange(len(idx)), j]
        ok = dmin <= tolerance
        d[idx[ok]] = cont["Z"][idx[ok]] - pts["Z"][j[ok]]
        match[idx[ok]] = j[ok]
        dist[idx[ok]] = dmin[ok]
    # A frame's points all carry its one water level, so every waterline point near the same
    # survey point repeats the same difference: keep one per (frame, survey point), the nearest.
    hit = np.nonzero(match >= 0)[0]
    if len(hit):
        keys = np.char.add(np.char.add(cont["frame"][hit].astype(str), "|"), match[hit].astype(str))
        o = np.lexsort((dist[hit], keys))
        _, first = np.unique(keys[o], return_index=True)
        keep = np.zeros(len(d), bool)
        keep[hit[o[first]]] = True
        d[~keep] = np.nan
    if dist_out is not None:
        dist_out.extend(dist[np.isfinite(d)].tolist())
    return d


def is_transect(desc):
    return str(desc).strip().lower() == "transect"


def waterlines_on_transects(cont, transects, tolerance=TRANSECT_TOLERANCE, min_points=TRANSECT_MIN_POINTS):
    """
    compare_rtk.py's check of the waterlines against RTK transects (transects from
    compare_rtk.load_transects): a waterline point counts if it lies within `tolerance` m of a
    transect line and between its first and last surveyed point; the RTK elevation is
    interpolated along the profile there. ONE value per frame with >= min_points such points:
    the median of waterline - RTK.
    -> (per-frame dict {d, cam, day, frame, points} of arrays, per-point mask of the points used).
    """
    n = len(cont["E"])
    zr = np.full(n, np.nan)
    on = np.zeros(n, bool)
    for t in transects:
        re_, rn = cont["E"] - t["c"][0], cont["N"] - t["c"][1]
        s = re_ * t["d"][0] + rn * t["d"][1]
        lat = np.abs(-re_ * t["d"][1] + rn * t["d"][0])
        m = (lat <= tolerance) & (s >= t["s"][0]) & (s <= t["s"][-1]) & ~on
        zr[m] = np.interp(s[m], t["s"], t["z"])
        on |= m
    fr = {"d": [], "cam": [], "day": [], "frame": [], "points": []}
    used = np.zeros(n, bool)
    idx = np.nonzero(on)[0]
    if len(idx):
        frames = cont["frame"][idx]
        o = np.argsort(frames, kind="stable")
        idx, frames = idx[o], frames[o]
        cuts = np.nonzero(frames[1:] != frames[:-1])[0] + 1
        for grp in np.split(idx, cuts):
            if len(grp) < min_points:
                continue
            fr["d"].append(float(np.median(cont["Z"][grp] - zr[grp])))
            fr["cam"].append(cont["cam"][grp[0]])
            fr["day"].append(cont["day"][grp[0]])
            fr["frame"].append(cont["frame"][grp[0]])
            fr["points"].append(len(grp))
            used[grp] = True
    out = {k: np.array(v, dtype=float if k == "d" else (int if k == "points" else object))
           for k, v in fr.items()}
    return out, used


def waterline_rows(cont, d):
    """Per camera and per day: frames, values, statistics of d (one value per point, or per frame
    for the transect check, where cont holds the per-frame arrays)."""
    ok = np.isfinite(d)
    rows = []
    for cam in sorted(set(cont["cam"][ok].tolist())):
        m = ok & (cont["cam"] == cam)
        rows.append(dict({"group_type": "waterline_camera", "group": cam,
                          "frames": int(len(set(cont["frame"][m].tolist())))}, **stats(d[m])))
    for day in sorted(set(cont["day"][ok].tolist())):
        for cam in sorted(set(cont["cam"][ok].tolist())):
            m = ok & (cont["day"] == day) & (cont["cam"] == cam)
            if m.any():
                rows.append(dict({"group_type": "waterline_day", "group": f"{day} {cam}",
                                  "frames": int(len(set(cont["frame"][m].tolist())))}, **stats(d[m])))
    if ok.any():
        rows.insert(0, dict({"group_type": "waterline_all", "group": "all",
                             "frames": int(len(set(cont["frame"][ok].tolist())))}, **stats(d[ok])))
    return rows


# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

def fv(v, fmt):
    """A number in `fmt`, or '-' when it is missing (NaN: too few values for it)."""
    try:
        return fmt.format(v) if v is not None and np.isfinite(v) else "-"
    except (TypeError, ValueError):
        return "-"


def fmt_row(r, what="cells"):
    if not r["n"]:
        return f"  {r['group']:>16s} {0:>7d}"
    few = f"  (few: n < {MIN_BAND_N}, not an estimate)" if r["n"] < MIN_BAND_N else ""
    return (f"  {r['group']:>16s} {r['n']:>7d} {r['median']:>+7.3f} {r['mean']:>+7.3f} "
            f"{fv(r['nmad'], '{:.3f}'):>6s} {r['rmse']:>6.3f} {fv(r['p5'], '{:+.2f}'):>7s} "
            f"{fv(r['p95'], '{:+.2f}'):>7s} "
            f"{100 * r['within_0.10']:>5.0f}% {100 * r['within_0.20']:>5.0f}% "
            f"{100 * r['within_0.50']:>5.0f}%{few}")


def table(lines, title, rows, what="cells"):
    if not rows:
        return
    lines.append("")
    lines.append(title)
    lines.append(f"  {'':>16s} {what:>7s} {'median':>7s} {'mean':>7s} {'NMAD':>6s} {'RMSE':>6s} "
                 f"{'p5':>7s} {'p95':>7s} {'10cm':>6s} {'20cm':>6s} {'50cm':>6s}")
    for r in rows:
        lines.append(fmt_row(r, what))


def headline_text(s, what):
    return (f"median {s['median']:+.3f} m, NMAD {fv(s['nmad'], '{:.3f} m')}, RMSE {s['rmse']:.3f} m, "
            f"p5/p95 {fv(s['p5'], '{:+.2f}')}/{fv(s['p95'], '{:+.2f}')} m, {100 * s['within_0.20']:.0f}% "
            f"within +/-0.20 m (n = {s['n']} {what})")


def clean_json(o):
    """NaN and infinities -> None: JSON has no NaN, and survey_products.py reads this file."""
    if isinstance(o, dict):
        return {k: clean_json(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean_json(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


# ---------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------

def nice_limit(d):
    a = np.abs(d[np.isfinite(d)])
    if not len(a):
        return 0.5
    p = float(np.percentile(a, 95))
    step = 0.05 if p < 0.5 else (0.1 if p < 1.5 else 0.25)
    return max(step, np.ceil(p / step) * step)


def draw_figure(path, ctx):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, ListedColormap
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator

    div = LinearSegmentedColormap.from_list("dem_minus_survey", DIVERGING, N=256)
    sand = LinearSegmentedColormap.from_list("sand_low_dark", SAND[::-1], N=256)
    lim = ctx["lim"]
    dem, xll, ytop, cell = ctx["dem"], ctx["xll"], ctx["ytop"], ctx["cell"]
    nrows, ncols = dem.shape
    ext = [xll, xll + ncols * cell, ytop - nrows * cell, ytop]
    W, H = 14.0, 9.6
    with plt.rc_context(PAGE_RC):
        fig = plt.figure(figsize=(W, H), dpi=150)
        fig.patch.set_facecolor("white")
        # --- header: title, the label as passed with its reason, warnings
        y = 0.982
        fig.text(0.02, y, ctx["title"], fontsize=14, fontweight="bold", va="top", color=INK)
        y -= 0.034
        lab = f"Label: {ctx['label']}  —  {ctx['why']}"
        for ln in textwrap.wrap(lab, 150):
            fig.text(0.02, y, ln, fontsize=11, va="top", color=INK,
                     fontweight="bold" if ln.startswith("Label:") else "normal")
            y -= 0.024
        for w in ctx["warn_lines"]:
            head_line = w.startswith(("DOWNGRADE", "TIME GAP", "DATUM"))
            for k, ln in enumerate(textwrap.wrap(w, 160)):
                fig.text(0.02 if k == 0 else 0.035, y, ln, fontsize=9.5, va="top",
                         color=INK if head_line else INK2,
                         fontweight="bold" if (head_line and k == 0) else "normal")
                y -= 0.021
        top = y - 0.035
        # --- map (left, with room for its legend below), headline text (top right),
        # histogram + elevation bands (bottom right)
        ax = fig.add_subplot(fig.add_gridspec(1, 1, left=0.085, right=0.35, bottom=0.17,
                                              top=top)[0, 0])
        # the headline text takes the height its lines need (smaller type if many), the two
        # charts the rest, so text never runs into a chart title
        box_text = "\n".join(ctx["box_lines"])
        n_lines = box_text.count("\n") + 1
        box_fs = 9.3 if n_lines <= 18 else max(7.5, 9.3 * 18 / n_lines)
        box_h = n_lines * box_fs * 1.35 / 72.0 / H
        charts_top = top - box_h - 0.07
        gs = fig.add_gridspec(1, 2, left=0.5, right=0.985, bottom=0.07, top=charts_top,
                              wspace=0.34)
        xa, xb, ya, yb = ctx["map_extent"]
        if ctx["survey_type"] == "dsm":
            nosurvey = np.isfinite(dem) & ~np.isfinite(ctx["diff"])
            face = np.where(nosurvey, 1.0, np.nan)
            ax.imshow(face, extent=ext, origin="upper", cmap=ListedColormap([BLANK_FACE]),
                      interpolation="nearest", zorder=1)
            if nosurvey.any():
                k = 4
                up = np.kron(nosurvey.astype(float), np.ones((k, k)))
                xs = xll + (np.arange(ncols * k) + 0.5) * cell / k
                ys = ytop - (np.arange(nrows * k) + 0.5) * cell / k
                ax.contourf(xs, ys, up, levels=[0.5, 1.5], colors="none", hatches=["////"],
                            zorder=2)
                for coll in ax.collections[-1:]:
                    try:
                        coll.set_edgecolor(BLANK_HATCH)
                        coll.set_linewidth(0)
                    except Exception:
                        pass
            im = ax.imshow(ctx["diff"], extent=ext, origin="upper", cmap=div, vmin=-lim, vmax=lim,
                           interpolation="nearest", zorder=3)
            cov, cext = ctx["cover_pad"], ctx["cover_extent"]
            edge = cov is not None and np.nanmax(cov) >= 0.5 and np.nanmin(cov) < 0.5
            if edge:
                cn, cm = cov.shape
                cxs = cext[0] + (np.arange(cm) + 0.5) * cell
                cys = cext[3] - (np.arange(cn) + 0.5) * cell
                ax.contour(cxs, cys, cov, levels=[0.5], colors=INK, linewidths=1.0, zorder=4)
            handles = [Line2D([], [], color=INK, lw=1.0, label="edge of the survey data")] \
                if edge else []
            if nosurvey.any():
                handles.append(Patch(facecolor=BLANK_FACE, edgecolor=BLANK_HATCH, hatch="////",
                                     lw=0.5, label=f"DEM cell, no survey value ({int(nosurvey.sum())})"))
            cb_label = "DEM − survey (m)"
        else:
            ax.imshow(dem, extent=ext, origin="upper", cmap=sand, interpolation="nearest",
                      vmin=ctx["zlim"][0], vmax=ctx["zlim"][1], zorder=1)
            P = ctx["points"]
            m = P["status"] == "compared"
            o = ~m
            ax.scatter(P["E"][o], P["N"][o], s=16, facecolors="none", edgecolors=MUTED, lw=0.8,
                       zorder=3)
            im = ax.scatter(P["E"][m], P["N"][m], c=P["diff"][m], cmap=div, vmin=-lim, vmax=lim,
                            s=34, edgecolors=INK, linewidths=0.5, zorder=4)
            handles = [Line2D([], [], ls="none", marker="o", mfc="#e4e3df", mec=INK, mew=0.5,
                              ms=6, label=f"survey point on the DEM ({int(m.sum())})"),
                       Line2D([], [], ls="none", marker="o", mfc="none", mec=MUTED, ms=5,
                              label=f"survey point off the DEM ({int(o.sum())}"
                              + (f", {ctx['beyond']} of them beyond the map" if ctx.get("beyond")
                                 else "") + ")")]
            cb_label = "DEM − survey (m)"
        cx, cy = ctx["camera_xy"]
        if xa - 50 <= cx <= xb + 50 and ya - 50 <= cy <= yb + 50:
            ax.plot([cx], [cy], marker="^", ms=10, color=INK, zorder=6, ls="none")
            handles.append(Line2D([], [], ls="none", marker="^", ms=8, color=INK,
                                  label="cameras (" + ctx["camera_source"] + ")"))
        if ctx["seam"] is not None:
            ax.axhline(ctx["seam"], color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=5)
            ax.text(xb - 2, ctx["seam"] + 2, "c2 ↑  seam  ↓ c1", ha="right", va="bottom",
                    fontsize=8.5, color=INK2, zorder=7,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.0))
        ax.set_xlim(xa, xb)
        ax.set_ylim(ya, yb)
        ax.set_aspect("equal")
        ax.set_xlabel("Easting (m, UTM 19N)")
        ax.set_ylabel("Northing (m, UTM 19N)")
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.xaxis.set_major_locator(MaxNLocator(3))
        ax.set_title(ctx["map_title"], loc="left", fontsize=10.5, color=INK)
        # colour bars beside the map: the difference on top, the DEM elevation (points) below
        two = ctx["survey_type"] == "points"
        ax.apply_aspect()
        pos = ax.get_position()                  # the map's real box (equal aspect)
        # the legend under the map at a fixed distance, clear of the axis label
        fig.legend(handles=handles, loc="upper left", fontsize=8.5, frameon=False,
                   bbox_to_anchor=(pos.x0 - 0.01, pos.y0 - 0.055), borderaxespad=0.0)
        # colour bars: at least ~0.6 of the page high, whatever the map's shape, so their
        # labels fit beside them
        mid = 0.5 * (pos.y0 + pos.y1)
        y0c = max(min(pos.y0, mid - 0.3), 0.17)
        y1c = min(max(pos.y1, mid + 0.3), top)
        cx0, ch = pos.x1 + 0.022, y1c - y0c
        # arrows on the colour bar where values lie beyond +/-lim (drawn in the end colours)
        dv = ctx["d"][np.isfinite(ctx["d"])]
        hi_, lo_ = bool((dv > lim).any()), bool((dv < -lim).any())
        ext = "both" if (hi_ and lo_) else "max" if hi_ else "min" if lo_ else "neither"
        cb = fig.colorbar(im, cax=fig.add_axes(
            [cx0, y0c + (0.53 if two else 0.25) * ch, 0.01, (0.42 if two else 0.5) * ch]), extend=ext)
        cb.set_label(cb_label + ("; arrows: beyond" if ext != "neither" else ""), color=INK2)
        cb.outline.set_edgecolor(AXIS)
        if two:
            sm = plt.cm.ScalarMappable(cmap=sand, norm=plt.Normalize(*ctx["zlim"]))
            cb2 = fig.colorbar(sm, cax=fig.add_axes([cx0, y0c + 0.05 * ch, 0.01, 0.38 * ch]))
            cb2.set_label("DEM elevation (m NAVD88)", color=INK2)
            cb2.outline.set_edgecolor(AXIS)

        # --- headline numbers
        fig.text(0.5, top, box_text, family="DejaVu Sans Mono", fontsize=box_fs, va="top",
                 ha="left", color=INK, linespacing=1.35)

        # --- histogram, bars coloured as the map
        d = ctx["d"][np.isfinite(ctx["d"])]
        hx = fig.add_subplot(gs[0, 0])
        span = 1.3 * lim
        step = span / 20.0
        edges = np.arange(-span, span + step / 2, step)
        dc = np.clip(d, -span + 1e-9, span - 1e-9)
        cnt, _ = np.histogram(dc, edges)
        mids = 0.5 * (edges[:-1] + edges[1:])
        hx.bar(mids, cnt, width=step * 0.92, color=div(np.clip((mids + lim) / (2 * lim), 0, 1)),
               edgecolor="none", zorder=2)
        hx.axvline(0, color=MUTED, lw=0.9, zorder=1)
        if len(d):
            med = float(np.median(d))
            hx.axvline(med, color=INK, lw=1.4, zorder=3)
            hx.text(med, hx.get_ylim()[1] * 0.985, f" median {med:+.2f} m ", ha="left", va="top",
                    fontsize=8.5, color=INK, zorder=4,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.0))
        hx.set_xlim(-span, span)
        hx.yaxis.set_major_locator(MaxNLocator(5, integer=True))
        hx.set_xlabel("DEM − survey (m); end bars hold everything beyond")
        hx.set_ylabel(ctx["what"])
        hx.set_title("Distribution", loc="left", fontsize=10.5, color=INK)
        hx.grid(True, axis="y", color=GRID)
        hx.set_axisbelow(True)
        for s in ("top", "right"):
            hx.spines[s].set_visible(False)

        # --- median by survey elevation band, NMAD bars
        bx = fig.add_subplot(gs[0, 1])
        bands = [r for r in ctx["rows"] if r["group_type"] == "elevation_band"]
        if bands:
            yc = np.array([float(r["group"].split(" to ")[0]) + BAND / 2 for r in bands])
            med = np.array([r["median"] for r in bands])
            nm = np.nan_to_num(np.array([r["nmad"] for r in bands], float))   # no bar below 3 values
            few = np.array([r["n"] < MIN_BAND_N for r in bands])
            bx.errorbar(med[~few], yc[~few], xerr=nm[~few], fmt="o", color=INK, ms=5,
                        ecolor=INK2, elinewidth=1.2, capsize=0, zorder=3)
            if few.any():
                bx.errorbar(med[few], yc[few], xerr=nm[few], fmt="o", mfc="white", mec=INK2,
                            ms=5, ecolor=AXIS, elinewidth=1.0, capsize=0, zorder=3)
            xm = max(lim, float(np.nanmax(np.abs(med) + nm))) * 1.15
            bx.set_xlim(-xm * 1.45, xm)
            for yy, r in zip(yc, bands):
                bx.text(-xm * 1.42, yy, f"n {r['n']}", va="center", ha="left", fontsize=8,
                        color=INK2)
            bx.set_ylim(yc.min() - BAND, yc.max() + BAND)
        bx.axvline(0, color=MUTED, lw=0.9, zorder=1)
        bx.set_xlabel("median DEM − survey (m), bar = ±NMAD")
        bx.set_ylabel("survey elevation band (m NAVD88)")
        bx.set_title("By elevation (open = fewer than 10)", loc="left", fontsize=10.5, color=INK)
        bx.grid(True, color=GRID)
        bx.set_axisbelow(True)
        for s in ("top", "right"):
            bx.spines[s].set_visible(False)
        if not len(d):
            for a in (hx, bx):
                a.cla()
                a.set_xticks([]); a.set_yticks([])
                for sp in a.spines.values():
                    sp.set_visible(False)
                a.text(0.5, 0.5, f"no {ctx['what']} to compare:\nthe DEM and the survey "
                       f"do not overlap", ha="center", va="center", fontsize=10.5, color=INK2,
                       transform=a.transAxes)
        fig.savefig(path, facecolor="white")
        plt.close(fig)


# ---------------------------------------------------------------------
# The comparison
# ---------------------------------------------------------------------

def compare(dem_path, survey_path, survey_type, label, why, name=None, output_dir=".",
            survey_date=None, photo_dates=None, contours=None, camera_eo=None, camera_io=None,
            spread=None, count=None, envelope_source=None, min_cover=0.5, tolerance=1.0,
            epsg=NOMINAL_EPSG, plot=True, transect_tolerance=TRANSECT_TOLERANCE,
            transect_split=TRANSECT_SPLIT, transect_min_points=TRANSECT_MIN_POINTS):
    """
    Compares a DEM with a survey and writes the outputs. camera_eo / camera_io: {cam: path}.
    Returns the headline (also written as NAME_comparison.json).
    """
    t0 = time.time()
    if label not in LABELS:
        raise ValueError(f"label must be one of {LABELS}")
    if survey_type not in ("dsm", "points"):
        raise ValueError("survey_type must be dsm or points")
    name = name or Path(survey_path).stem
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    stem = dem_stem(dem_path)
    print(f"comparison        : {name}  ({survey_type})")
    print(f"label             : {label} -- {why}")

    dem, xll, ytop, cell, dh = read_grid(dem_path)
    nrows, ncols = dem.shape
    print(f"DEM               : {dem_path}  {ncols} x {nrows} at {cell:g} m, "
          f"{int(np.isfinite(dem).sum())} cells with a value")
    grids = {}
    for key, p in (("spread", spread or stem + "_spread.asc"), ("count", count or stem + "_count.asc")):
        if p and Path(p).exists():
            g, gh = read_asc(p)
            if same_grid(gh, dh):
                grids[key] = g
                print(f"DEM {key:14s}: {p}")
            else:
                print(f"NOTE: {p} is not on the DEM's grid: not used")
    if not photo_dates:
        info_p = Path(stem + "_info.json")
        if info_p.exists():
            try:
                info = json.loads(info_p.read_text())
                if info.get("first_date") and info.get("last_date"):
                    photo_dates = (info["first_date"], info["last_date"])
                    print(f"photo dates       : {photo_dates[0]} to {photo_dates[1]} "
                          f"(from {info_p.name})")
            except (OSError, ValueError):
                pass

    cams = load_cameras(camera_eo or {}, camera_io or {}) if camera_eo else {}
    for cam, c in sorted(cams.items()):
        print(f"camera {cam:11s}: EO {Path(c['eo_path']).name}, IO "
              f"{Path(c['io_path']).name if c['io_path'] else 'not found (footprint not used)'}")
    if cams:
        camera_xy = tuple(np.mean([c["eo"][:2] for c in cams.values()], axis=0))
        camera_source = "from the EO files given"
    else:
        camera_xy, camera_source = TOWER_EN, "live tower position"

    cont = None
    if contours:
        print(f"waterlines        : reading {contours} ...")
        cont = read_contours(contours, *(photo_dates or (None, None)))
        print(f"waterlines        : {len(cont['E'])} georectified points"
              + (f" {photo_dates[0]} to {photo_dates[1]}" if photo_dates else ""))
    cell_cam = cell_majority_camera(cont, xll, ytop, cell, nrows, ncols) \
        if (cont is not None and len(cont["E"]) and not cams) else None

    # -- the label: what can be checked
    findings = audit_label(label, survey_path, survey_type, survey_date, cams, envelope_source,
                           dem_path)
    outside, from_mid, gap_text = time_gap(survey_date, photo_dates)
    worst = label
    for sug, _ in findings:
        if sug and _rank(sug) > _rank(worst):
            worst = sug
    warn_lines = []
    if worst != label:
        warn_lines.append(f"DOWNGRADE WARNING: passed {label}, but the checks below suggest "
                          f"{worst} (a comparison is never better than its weakest step):")
        for sug, txt in findings:
            if sug and _rank(sug) > _rank(label):
                warn_lines.append(f"  - [{sug}] {txt}")
    if outside is not None and outside > TIME_GAP_WARN_DAYS:
        warn_lines.append(f"TIME GAP: {gap_text}. The beach changes; part of the difference may be "
                          f"real change, not error.")
    for w in warn_lines:
        print(w)

    lines = [f"SURVEY COMPARISON: {name}", "=" * (19 + len(name)),
             f"LABEL      : {label}", f"WHY        : {why}"]
    if worst != label:
        lines += warn_lines
    else:
        lines.append("checks     : nothing the script can see contradicts the label"
                     + (" (no camera EO, envelope source or provenance to check)"
                        if not cams and not envelope_source else ""))
    notes = [txt for sug, txt in findings if not (sug and _rank(sug) > _rank(label))]
    for n_ in notes:
        lines.append(f"note       : {n_}")
    lines += [f"time gap   : {gap_text}",
              f"DEM        : {dem_path} ({ncols} x {nrows} cells at {cell:g} m, "
              f"{int(np.isfinite(dem).sum())} with a value)",
              f"survey     : {survey_path} ({survey_type})"]

    head = {"name": name, "label": label, "why": why, "suggested_label": worst,
            "downgrade_warning": worst != label, "warnings": warn_lines,
            "notes": notes, "survey": str(survey_path), "survey_type": survey_type,
            "dem": str(dem_path), "survey_date": survey_date,
            "photo_dates": list(photo_dates) if photo_dates else None,
            "days_outside_photo_window": outside, "days_from_photo_window_middle": from_mid,
            "time_gap": gap_text, "sign": "DEM - survey (positive = DEM too high)"}

    if not np.isfinite(dem).any():
        sys.exit(f"{dem_path}: no cell has a value -- nothing to compare")
    zlo, zhi = np.nanpercentile(dem, 2), np.nanpercentile(dem, 98)
    fin_rr, fin_cc = np.nonzero(np.isfinite(dem))
    bx0, bx1 = xll + fin_cc.min() * cell, xll + (fin_cc.max() + 1) * cell
    by0, by1 = ytop - (fin_rr.max() + 1) * cell, ytop - fin_rr.min() * cell
    files = {}

    if survey_type == "dsm":
        print(f"survey            : reading {survey_path} ...")
        sg, sx0, sy0, sc = read_survey(survey_path)
        crs = survey_crs_text(survey_path)
        vert, vnote = survey_vertical_text(survey_path)
        print(f"survey            : {sg.shape[1]} x {sg.shape[0]} at {sc:g} m"
              + (f", CRS as stated: {crs}" if crs else "") + f"; vertical: {vert or vnote}")
        lines.append(f"survey grid: {sg.shape[1]} x {sg.shape[0]} cells at {sc:g} m; horizontal CRS "
                     f"as stated in the file: {crs or 'not stated'}; compared in the same UTM 19N "
                     f"metres as the DEM, no shift applied (the DEM grid is in the frame of the "
                     f"calibration's GCPs, which were surveyed in NAD83(2011) / UTM 19N like the "
                     f"surveys; an EPSG:{epsg} tag on the GeoTIFFs is"
                     + (" nominal)" if epsg == 32619 else " as given)"))
        lines.append(f"survey vertical datum: {vert}" if vert else f"survey vertical datum: {vnote}")
        if vnote:
            head["vertical_datum_note"] = vnote
        # map extent: the DEM, the cameras if near, padded; the survey is aggregated over it all
        xa, xb, ya, yb = bx0, bx1, by0, by1
        cx, cy = camera_xy
        if np.hypot(cx - (xa + xb) / 2, cy - (ya + yb) / 2) < 600:
            xa, xb, ya, yb = min(xa, cx), max(xb, cx), min(ya, cy), max(yb, cy)
        pad = 12.0
        xa, xb, ya, yb = xa - pad, xb + pad, ya - pad, yb + pad
        pc0 = int(np.floor((xa - xll) / cell)); pc1 = int(np.ceil((xb - xll) / cell))
        pr0 = int(np.floor((ytop - yb) / cell)); pr1 = int(np.ceil((ytop - ya) / cell))
        pc0, pr0 = min(pc0, 0), min(pr0, 0)
        pc1, pr1 = max(pc1, ncols), max(pr1, nrows)
        sv_pad, cov_pad, how = survey_on_grid(sg, sx0, sy0, sc, xll + pc0 * cell, ytop - pr0 * cell,
                                              cell, pr1 - pr0, pc1 - pc0)
        sv = sv_pad[-pr0:-pr0 + nrows, -pc0:-pc0 + ncols]
        cov = cov_pad[-pr0:-pr0 + nrows, -pc0:-pc0 + ncols]
        cover_extent = [xll + pc0 * cell, xll + pc1 * cell, ytop - pr1 * cell, ytop - pr0 * cell]
        print(f"survey on grid    : {how}; a cell needs >= {100 * min_cover:.0f}% covered")
        lines.append(f"survey on the DEM grid: {how}; a cell counts if the survey covers "
                     f">= {100 * min_cover:.0f}% of it")
        has_dem = np.isfinite(dem)
        cmp_ = has_dem & np.isfinite(sv) & (cov >= min_cover)
        diff = np.where(cmp_, dem - sv, np.nan)
        partial = has_dem & (cov > 0) & (cov < min_cover)
        # compared, but the survey covers only part of the cell (the lidar's water edge): the
        # median of the dry part sits up the slope, so these cells read the DEM a little LOW
        cmp_part = cmp_ & (cov < FULL_COVER)
        none_ = has_dem & (cov == 0)
        n_dem = int(has_dem.sum())
        rr, cc = np.nonzero(has_dem)
        E = xll + (cc + 0.5) * cell
        N = ytop - (rr + 0.5) * cell
        cam_all, cam_how = assign_camera(E, N, dem[rr, cc], cams, cell_cam, (xll, ytop, cell))
        cam_grid = np.full(dem.shape, "", dtype=object)
        cam_grid[rr, cc] = cam_all
        dist_grid = np.hypot((xll + (np.arange(ncols) + 0.5) * cell)[None, :] - camera_xy[0],
                             (ytop - (np.arange(nrows) + 0.5) * cell)[:, None] - camera_xy[1])
        # coverage: where the survey has no value
        lines += ["", "COVERAGE (lidar sees only the beach that was dry at flight time)",
                  f"  DEM cells with a value      : {n_dem}",
                  f"  compared (survey >= {100 * min_cover:.0f}%)    : {int(cmp_.sum())} "
                  f"({100 * cmp_.sum() / max(n_dem, 1):.0f}%)",
                  f"    of which partly covered   : {int(cmp_part.sum())} (survey on "
                  f"{100 * min_cover:.0f}-99% of the cell: its value is the median of the covered, "
                  f"higher part, so DEM - survey reads a little low there; see BY SURVEY COVER)",
                  f"  survey on < {100 * min_cover:.0f}% of the cell  : {int(partial.sum())} (not compared)",
                  f"  no survey value at all      : {int(none_.sum())}"]
        if (has_dem & ~cmp_).any():
            zn = dem[has_dem & ~cmp_]
            lines.append(f"  DEM elevation where no survey: median {np.median(zn):+.2f} m, "
                         f"p5-p95 {np.percentile(zn, 5):+.2f} to {np.percentile(zn, 95):+.2f} m")
        if cmp_.any():
            zc = sv[cmp_]
            lines.append(f"  survey elevation where compared: p2 {np.percentile(zc, 2):+.2f} m "
                         f"(about the survey's seaward edge inside the DEM), median "
                         f"{np.median(zc):+.2f} m, max {zc.max():+.2f} m")
        lines.append(f"  {'DEM elevation band':>20s} {'cells':>6s} {'compared':>9s} {'share':>6s} "
                     f"{'partly covered':>15s}")
        zb = np.floor(dem / BAND) * BAND
        cov_rows = []
        for lo in sorted(set(zb[has_dem].tolist())):
            m = has_dem & (zb == lo)
            lines.append(f"  {band_name(lo):>20s} {int(m.sum()):>6d} {int((m & cmp_).sum()):>9d} "
                         f"{100 * (m & cmp_).sum() / m.sum():>5.0f}% {int((m & cmp_part).sum()):>15d}")
            cov_rows.append({"band": band_name(lo), "dem_cells": int(m.sum()),
                             "compared": int((m & cmp_).sum()),
                             "compared_partial_cover": int((m & cmp_part).sum())})
        lines.append(f"  {'camera':>20s} {'cells':>6s} {'compared':>9s} {'share':>6s}")
        for c in sorted(set(cam_all.tolist())):
            m = has_dem & (cam_grid == c)
            lines.append(f"  {c:>20s} {int(m.sum()):>6d} {int((m & cmp_).sum()):>9d} "
                         f"{100 * (m & cmp_).sum() / max(m.sum(), 1):>5.0f}%")
        d = diff[cmp_]
        zref = sv[cmp_]
        cam_v = cam_grid[cmp_]
        dist_v = dist_grid[cmp_]
        part_v = cmp_part[cmp_]
        rows = standard_groups(d, zref, cam_v, dist_v,
                               extra=[("survey_cover", "full", ~part_v),
                                      ("survey_cover", f"{100 * min_cover:.0f}-99%", part_v)])
        rep_g = repeatability_groups(np.where(cmp_, diff, np.nan),
                                     grids.get("spread"), grids.get("count"))
        rows += group_rows(np.where(cmp_, diff, np.nan), rep_g)
        what = "cells"
        full = stats(diff[cmp_ & ~cmp_part])
        head.update({"dem_cells": n_dem, "dem_cells_compared": int(cmp_.sum()),
                     "dem_cells_compared_partial_cover": int(cmp_part.sum()),
                     "full_cover": {k: full[k] for k in STAT_KEYS},
                     "dem_cells_partial_survey": int(partial.sum()),
                     "dem_cells_no_survey": int(none_.sum()), "coverage_by_band": cov_rows,
                     "survey_crs": crs, "survey_on_grid": how})
        # outputs: the difference grid
        if cmp_.any():
            from dem_from_contours import write_ascii_grid
            a = out / f"{name}_dem_minus_survey.asc"
            write_ascii_grid(str(a), np.flipud(diff), xll, ytop - nrows * cell, cell)
            t = out / f"{name}_dem_minus_survey.tif"
            write_geotiff(t, diff, xll, ytop - nrows * cell, cell, epsg,
                          f"{name}: DEM - survey, m (positive = DEM high); a difference, not a "
                          f"NAVD88 height; EPSG:{epsg}")
            files.update({"diff_asc": str(a), "diff_tif": str(t)})
            print(f"wrote {a}\nwrote {t}")
        map_ctx = {"diff": diff, "cover_pad": cov_pad, "cover_extent": cover_extent,
                   "map_extent": (xa, xb, ya, yb),
                   "map_title": "DEM − survey per DEM cell (red = DEM high)"}
        wl_grid = (sg, sx0, sy0, sc)
        pts = None
    else:
        pts = read_points(survey_path)
        crs_set = sorted(set(c for c in pts["crs"] if c))
        not_navd = [c for c in crs_set if "NAVD88" not in c.upper()]
        print(f"survey            : {len(pts['E'])} points, {pts['format']}"
              + (f", {pts['skipped']} rows without coordinates skipped" if pts["skipped"] else ""))
        lines.append(f"survey points: {len(pts['E'])} ({pts['format']})"
                     + (f"; {pts['skipped']} rows without numeric coordinates skipped"
                        if pts["skipped"] else ""))
        if crs_set:
            lines.append("survey CRS (CS name column): " + "; ".join(crs_set))
        if not_navd:
            w = (f"DATUM WARNING: the survey's CS name does not say NAVD88 ({'; '.join(not_navd)}): "
                 f"the vertical difference includes the datum difference")
            warn_lines.append(w); lines.append(w); print(w)
        elif not crs_set:
            lines.append("survey CRS: not stated in the file (assumed UTM 19N, NAVD88)")
        val, status = dem_at_points(dem, xll, ytop, cell, pts["E"], pts["N"])
        diffp = val - pts["Z"]
        near = nearest_cell_distance(dem, xll, ytop, cell, pts["E"], pts["N"])
        camp, cam_how = assign_camera(pts["E"], pts["N"], pts["Z"], cams, cell_cam,
                                      (xll, ytop, cell))
        distp = np.hypot(pts["E"] - camera_xy[0], pts["N"] - camera_xy[1])
        ok = status == "compared"
        lines += ["", "COVERAGE",
                  f"  points                      : {len(ok)}",
                  f"  on the DEM (compared)       : {int(ok.sum())}",
                  f"  in the DEM extent, empty cell: {int((status == 'dem_cell_without_value').sum())}",
                  f"  outside the DEM extent      : {int((status == 'outside_dem_extent').sum())}",
                  "  DEM value: bilinear over the finite cells among the 4 surrounding cell "
                  "centres, only where the point's own cell has a value"]
        descs = sorted(set(pts["desc"]))
        lines.append(f"  {'description':>20s} {'points':>6s} {'compared':>9s}")
        for ds in descs:
            m = np.array([x == ds for x in pts["desc"]])
            lines.append(f"  {ds[:20]:>20s} {int(m.sum()):>6d} {int((m & ok).sum()):>9d}")
        desc_arr = np.array(pts["desc"], dtype=object)
        extra = [("description", ds, ok & (desc_arr == ds)) for ds in descs]
        codes = sorted(set(c for c in pts["code"] if c))
        if len(codes) > 1:
            code_arr = np.array(pts["code"], dtype=object)
            extra += [("code", c, ok & (code_arr == c)) for c in codes]
        d_all = np.where(ok, diffp, np.nan)
        rows = standard_groups(d_all, pts["Z"], camp, distp, extra)
        if grids:
            ci = np.clip(np.floor((pts["E"] - xll) / cell).astype(int), 0, ncols - 1)
            ri = np.clip(np.floor((ytop - pts["N"]) / cell).astype(int), 0, nrows - 1)
            rows += group_rows(d_all, repeatability_groups(
                d_all, grids["spread"][ri, ci] if "spread" in grids else None,
                grids["count"][ri, ci] if "count" in grids else None))
        d = diffp[ok]
        what = "points"
        # every point, compared or not, with why
        pp = out / f"{name}_points.csv"
        with open(pp, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["name", "code", "description", "easting", "northing", "survey_z_navd88",
                        "dem_z_navd88", "dem_minus_survey_m", "status", "nearest_dem_cell_m",
                        "camera", "distance_from_camera_m"])
            for k in range(len(ok)):
                w.writerow([pts["name"][k], pts["code"][k], pts["desc"][k], f"{pts['E'][k]:.3f}",
                            f"{pts['N'][k]:.3f}", f"{pts['Z'][k]:.3f}",
                            f"{val[k]:.3f}" if np.isfinite(val[k]) else "",
                            f"{diffp[k]:+.3f}" if np.isfinite(diffp[k]) else "", status[k],
                            f"{near[k]:.1f}", camp[k], f"{distp[k]:.1f}"])
        files["points_csv"] = str(pp)
        print(f"wrote {pp}")
        off = np.nonzero(~ok)[0]
        lines += ["", f"POINTS NOT COMPARED ({len(off)}; name, description, z, why, nearest DEM "
                      f"cell with a value)"]
        for k in off:
            lines.append(f"  {pts['name'][k]:>14s}  {pts['desc'][k][:18]:18s} z {pts['Z'][k]:+6.2f}  "
                         f"{status[k]:24s} {near[k]:7.1f} m")
        head.update({"points": int(len(ok)), "points_compared": int(ok.sum()),
                     "points_outside_dem": int((status == "outside_dem_extent").sum()),
                     "points_in_empty_dem_cells": int((status == "dem_cell_without_value").sum()),
                     "survey_crs": "; ".join(crs_set), "survey_format": pts["format"]})
        # map extent: the DEM, points within 150 m of it, the cameras if near
        xa, xb, ya, yb = bx0, bx1, by0, by1
        nearp = (pts["E"] > bx0 - 150) & (pts["E"] < bx1 + 150) & \
                (pts["N"] > by0 - 150) & (pts["N"] < by1 + 150)
        if nearp.any():
            xa, xb = min(xa, pts["E"][nearp].min()), max(xb, pts["E"][nearp].max())
            ya, yb = min(ya, pts["N"][nearp].min()), max(yb, pts["N"][nearp].max())
        cx, cy = camera_xy
        if np.hypot(cx - (xa + xb) / 2, cy - (ya + yb) / 2) < 600:
            xa, xb, ya, yb = min(xa, cx), max(xb, cx), min(ya, cy), max(yb, cy)
        pad = 10.0
        map_ctx = {"points": {"E": pts["E"], "N": pts["N"], "diff": diffp, "status": status},
                   "map_extent": (xa - pad, xb + pad, ya - pad, yb + pad), "zlim": (zlo, zhi),
                   "map_title": "Survey points over the DEM (red = DEM high)",
                   "beyond": int((~nearp).sum())}
        wl_grid = None

    lines.append(f"camera of a cell/point: {cam_how}")
    lines.append(f"distance from the cameras: from {camera_xy[0]:.1f} E, {camera_xy[1]:.1f} N "
                 f"({camera_source})")
    s_all = stats(d)
    head.update({k: s_all[k] for k in STAT_KEYS})
    if s_all["n"]:
        hl = headline_text(s_all, what)
        if s_all["n"] < MIN_BAND_N:
            hl += f" -- ONLY {s_all['n']}: not a robust estimate"
        print(f"HEADLINE          : DEM - survey {hl}")
    else:
        hl = f"no {what} to compare -- the DEM and the survey do not overlap"
        print(f"HEADLINE          : {hl}")
    lines[4:4] = ["", f"HEADLINE   : DEM - survey {hl}",
                  "             positive = DEM too HIGH; NMAD = 1.4826 x median |d - median| "
                  "(robust sd)", ""]
    table(lines, "OVERALL", [r for r in rows if r["group_type"] == "overall"], what)
    table(lines, f"BY SURVEY ELEVATION ({BAND} m bands, m NAVD88)",
          [r for r in rows if r["group_type"] == "elevation_band"], what)
    table(lines, "BY CAMERA", [r for r in rows if r["group_type"] == "camera"], what)
    table(lines, "BY DISTANCE FROM THE CAMERAS", [r for r in rows if r["group_type"] == "distance"], what)
    table(lines, "BY DESCRIPTION", [r for r in rows if r["group_type"] == "description"], what)
    table(lines, "BY CODE", [r for r in rows if r["group_type"] == "code"], what)
    table(lines, "BY SURVEY COVER OF THE CELL (partly covered cells read a little low: their survey value "
                 "is the median of the dry, higher part)",
          [r for r in rows if r["group_type"] == "survey_cover"], what)
    table(lines, "BY THE DEM'S OWN SPREAD (16-84 percentile of its frames)",
          [r for r in rows if r["group_type"] == "dem_spread"], what)
    table(lines, "BY THE DEM'S FRAMES PER CELL", [r for r in rows if r["group_type"] == "dem_frames"], what)

    # -- waterlines against the survey
    wl_rows = []
    wl_frames_csv = None
    if cont is not None and len(cont["E"]):
        transects = []
        if survey_type == "points" and sum(is_transect(x) for x in pts["desc"]) >= 3:
            from compare_rtk import load_transects
            try:
                transects = load_transects(survey_path, transect_split)
            except (KeyError, ValueError) as exc:          # not an Emlid file after all
                lines.append(f"  (transects not read: {exc})")
        if transects:
            # compare_rtk.py's check: one value per frame, the RTK interpolated along the profile
            per, okw = waterlines_on_transects(cont, transects, transect_tolerance, transect_min_points)
            dw = per["d"]
            wl_cont = per
            unit = "frames"
            n_other = sum(not is_transect(x) for x in pts["desc"])
            how_w = (f"compare_rtk.py's transect check: {len(transects)} RTK transects (consecutive "
                     f"'Transect' points, split at gaps > {transect_split:g} m); a waterline point "
                     f"counts within {transect_tolerance:g} m of a transect and inside its surveyed "
                     f"stretch (never extrapolated), the RTK interpolated along the profile; ONE value "
                     f"per frame with >= {transect_min_points} points (their median)"
                     + (f"; the {n_other} other points (wrack lines etc.: the swash limit, not the "
                        f"sand) are not used for the waterlines" if n_other else ""))
            method = "frames on RTK transects (compare_rtk.py)"
        else:
            dists = []
            dw = waterlines_vs_survey(cont, survey_type, wl_grid, pts, tolerance, dist_out=dists)
            okw = np.isfinite(dw)
            wl_cont = cont
            unit = "points"
            if survey_type == "dsm":
                how_w = "the survey sampled bilinearly at each waterline point"
                method = "waterline points on the DSM"
            else:
                md = float(np.median(dists)) if dists else float("nan")
                how_w = (f"the nearest survey point within {tolerance:g} m of each waterline point, one "
                         f"value per frame and survey point (a frame's points share one water level); "
                         f"NOT slope-corrected: the pairs lie a median {fv(md, '{:.2f}')} m apart, and "
                         f"1 m along a 1:10 beach face is 0.1 m of elevation")
                method = f"nearest point within {tolerance:g} m, not slope-corrected"
        nb, nt = int((cont["beach"] & okw).sum()), int((~cont["beach"] & okw).sum())
        if nb and not nt:
            setup_mode = "setup-corrected"
            col = "beach_elevation_navd88 (still water + wave setup)"
        elif nt and not nb:
            setup_mode = "NO setup correction"
            col = ("tide_elevation_navd88 (still water only: no setup correction, so lines marked "
                   "by the swash read LOW by about the setup)")
        else:
            setup_mode = "setup on some rows only"
            col = (f"beach_elevation_navd88 for {nb} points, tide_elevation_navd88 (no setup) for "
                   f"{nt}")
        wl_rows = waterline_rows(wl_cont, dw)
        lines += ["", "WATERLINES vs SURVEY (waterline elevation - survey; positive = waterline "
                      "elevation too high)",
                  f"  {contours}: {len(cont['E'])} points, {int(okw.sum())} used"]
        lines += ["  " + ln for ln in textwrap.wrap(how_w, 100)]
        lines.append(f"  elevation used: {col}")
        lines.append(f"  {'day camera':>16s} {'frames':>6s} {unit:>7s} {'median':>7s} "
                     f"{'NMAD':>6s} {'RMSE':>6s} {'20cm':>6s}")
        for r in wl_rows:
            lines.append(f"  {r['group']:>16s} {r['frames']:>6d} {r['n']:>7d} {r['median']:>+7.3f} "
                         f"{fv(r['nmad'], '{:.3f}'):>6s} {r['rmse']:>6.3f} {100 * r['within_0.20']:>5.0f}%"
                         + ("  (few)" if r["n"] < MIN_BAND_N else ""))
        head["waterlines"] = {"file": str(contours), "points": int(len(cont["E"])),
                              "points_on_survey": int(okw.sum()), "elevation_used": col,
                              "setup": setup_mode, "method": method, "how": how_w, "unit": unit}
        if transects:
            wl_frames_csv = out / f"{name}_waterline_frames.csv"
            with open(wl_frames_csv, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["frame", "camera", "day", "points_on_transects", "waterline_minus_rtk_m"])
                for k in range(len(dw)):
                    w.writerow([per["frame"][k], per["cam"][k], per["day"][k], per["points"][k],
                                f"{dw[k]:+.4f}"])
            files["waterline_frames_csv"] = str(wl_frames_csv)
        if wl_rows:
            head["waterlines"].update({k: wl_rows[0][k] for k in STAT_KEYS})
            head["waterlines"]["frames"] = wl_rows[0]["frames"]
            few_w = wl_rows[0]["n"] < MIN_BAND_N
            print(f"waterlines        : waterline - survey median {wl_rows[0]['median']:+.3f} m, "
                  f"NMAD {fv(wl_rows[0]['nmad'], '{:.3f} m')} ({wl_rows[0]['n']} {unit}"
                  + ("" if unit == "frames" else f", {wl_rows[0]['frames']} frames")
                  + f"{'; TOO FEW: not an estimate' if few_w else ''}); "
                  f"{method}; {col.split(' (')[0]}")
        else:
            head["waterlines"].update({"n": 0, "frames": 0})
            print(f"waterlines        : no waterline on the survey ({method})")

    # -- write
    txt = out / f"{name}_comparison.txt"
    txt.write_text("\n".join(lines) + "\n")
    files["txt"] = str(txt)
    cs = out / f"{name}_comparison.csv"
    with open(cs, "w", newline="") as f:
        fields = ["name", "label", "quantity", "group_type", "group", "frames"] + list(STAT_KEYS)
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(dict(r, name=name, label=label, quantity="dem_minus_survey", frames=""))
        for r in wl_rows:
            w.writerow(dict(r, name=name, label=label, quantity="waterline_minus_survey"))
    files["csv"] = str(cs)
    head["by_camera"] = {r["group"]: {k: r[k] for k in ("n", "median", "nmad")}
                         for r in rows if r["group_type"] == "camera"}
    print(f"wrote {txt}\nwrote {cs}")

    if plot:
        lim = nice_limit(d) if len(d) else 0.5
        box = [f"Label: {label}" + (f"   (checks suggest {worst})" if worst != label else ""),
               "",
               f"DEM − survey, {what} compared: {s_all['n']}"
               + (" -- too few for a robust estimate" if 0 < s_all["n"] < MIN_BAND_N else "")]
        if s_all["n"]:
            box += [f"  median {s_all['median']:+.3f} m    mean {s_all['mean']:+.3f} m",
                    f"  NMAD   {fv(s_all['nmad'], '{:.3f} m')}    RMSE {s_all['rmse']:.3f} m",
                    f"  p5 / p95  {fv(s_all['p5'], '{:+.2f}')} / {fv(s_all['p95'], '{:+.2f}')} m",
                    f"  within ±0.10 / 0.20 / 0.50 m:  {100 * s_all['within_0.10']:.0f}% / "
                    f"{100 * s_all['within_0.20']:.0f}% / {100 * s_all['within_0.50']:.0f}%"]
        if survey_type == "dsm":
            rest = head["dem_cells"] - head["dem_cells_compared"]
            box.append(f"  DEM cells with survey: {head['dem_cells_compared']} of {head['dem_cells']}"
                       + (f" ({rest} without: no survey value there)" if rest else ""))
            fc = head.get("full_cover") or {}
            if head.get("dem_cells_compared_partial_cover"):
                box.append(f"  of which partly covered: {head['dem_cells_compared_partial_cover']}; fully "
                           f"covered only: median {fv(fc.get('median'), '{:+.3f}')} m (n {fc.get('n')})")
        else:
            box.append(f"  points on the DEM: {head['points_compared']} of {head['points']} "
                       f"({head['points_outside_dem']} outside it)")
        camr = [r for r in rows if r["group_type"] == "camera"]
        for r in camr:
            box.append(f"  {r['group']:>7s}: median {r['median']:+.3f} m, NMAD {fv(r['nmad'], '{:.3f} m')}, "
                       f"n {r['n']}")
        if wl_rows:
            w0 = wl_rows[0]
            wu = head["waterlines"]["unit"]
            box += ["", f"waterlines − survey: median {w0['median']:+.3f} m, NMAD "
                        f"{fv(w0['nmad'], '{:.3f} m')}",
                    (f"  ({w0['n']} frames, one value each; {head['waterlines']['points_on_survey']} line "
                     f"points on the transects" if wu == "frames" else f"  ({w0['n']} {wu}, {w0['frames']} frames")
                    + ("; TOO FEW: not an estimate" if w0["n"] < MIN_BAND_N else "") + ")",
                    f"  {head['waterlines']['method']}; {head['waterlines']['setup']}"]
        box += ["", "\n".join(textwrap.wrap(gap_text, 78))]
        ctx = dict(map_ctx, title=f"{name}: DEM − survey ({survey_type})", label=label,
                   why=why, warn_lines=warn_lines, lim=lim, dem=dem, xll=xll, ytop=ytop, cell=cell,
                   survey_type=survey_type, camera_xy=camera_xy, camera_source=camera_source,
                   seam=SEAM_NORTHING if (not cams and cell_cam is None) else None,
                   box_lines=box, d=d, rows=rows, what=what)
        png = out / f"{name}_comparison.png"
        draw_figure(png, ctx)
        files["png"] = str(png)
        print(f"wrote {png}")
    head["files"] = files
    js = out / f"{name}_comparison.json"
    files["json"] = str(js)
    js.write_text(json.dumps(clean_json(head), indent=1, allow_nan=False,
                             default=lambda o: None if o is None else str(o)))
    print(f"wrote {js}")
    print(f"done              : {time.time() - t0:.1f} s")
    return head


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--dem", required=True, help="the DEM, STEM_dem.asc (dem_from_contours.py)")
    ap.add_argument("--survey", required=True,
                    help="survey: lidar DSM (.tif/.asc) or points (.csv: Emlid, CIRN num,E,N,Z, E/N/Z)")
    ap.add_argument("--survey-type", choices=["dsm", "points"], default=None,
                    help="default: dsm for .tif/.asc, points otherwise")
    ap.add_argument("--name", default=None, help="output name (default the survey file's stem)")
    ap.add_argument("--output-dir", default=".")
    ap.add_argument("--label", required=True, choices=LABELS,
                    help="how independent the comparison is; printed exactly as given")
    ap.add_argument("--why", required=True, help="the reason for the label, printed next to it")
    ap.add_argument("--survey-date", default=None, help="YYYY-MM-DD the survey was made")
    ap.add_argument("--photo-dates", nargs=2, metavar=("FIRST", "LAST"), default=None,
                    help="first and last photo date (default: from STEM_info.json if there)")
    ap.add_argument("--contours", default=None,
                    help="georectified waterline points (e.g. contour_points_ground_filtered.csv): "
                         "also compared with the survey, per day and camera")
    ap.add_argument("--camera-eo", nargs="+", default=None, metavar="CAM=PATH",
                    help="the EO files the DEM was built with: camera footprints, position, and "
                         "the label checks (e.g. c1=calibration/CACO03_c1_20250123_EO.yaml)")
    ap.add_argument("--camera-io", nargs="+", default=None, metavar="CAM=PATH",
                    help="lens files (default: the EO's setup IO, else CACO05_<cam>_20240801_IO.yaml)")
    ap.add_argument("--spread", default=None, help="STEM_spread.asc (default: next to the DEM)")
    ap.add_argument("--count", default=None, help="STEM_count.asc (default: next to the DEM)")
    ap.add_argument("--envelope-source", default=None,
                    help="the survey the detection search envelope was placed with "
                         "(detect_original_view.py); checked against --survey")
    ap.add_argument("--min-cover", type=float, default=0.5,
                    help="share of a DEM cell the survey must cover to be compared (default 0.5)")
    ap.add_argument("--tolerance", type=float, default=1.0,
                    help="isolated survey points (no transects): max distance (m) from a waterline point "
                         "to a survey point (not slope-corrected)")
    ap.add_argument("--transect-tolerance", type=float, default=TRANSECT_TOLERANCE,
                    help=f"RTK transects: max distance (m) of a waterline point from a transect line "
                         f"(default {TRANSECT_TOLERANCE:g}, as compare_rtk.py)")
    ap.add_argument("--transect-split", type=float, default=TRANSECT_SPLIT,
                    help=f"RTK transects: a gap between consecutive points that starts a new transect, m "
                         f"(default {TRANSECT_SPLIT:g})")
    ap.add_argument("--transect-min-points", type=int, default=TRANSECT_MIN_POINTS,
                    help=f"RTK transects: waterline points a frame needs on them to count "
                         f"(default {TRANSECT_MIN_POINTS})")
    ap.add_argument("--epsg", type=int, default=NOMINAL_EPSG,
                    help="horizontal EPSG written into the difference GeoTIFF (default 32619)")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    st = args.survey_type or ("dsm" if Path(args.survey).suffix.lower() in (".tif", ".tiff", ".asc")
                              else "points")
    for p in (args.dem, args.survey) + ((args.contours,) if args.contours else ()):
        if not Path(p).exists():
            sys.exit(f"not found: {p}")
    print("expected run time : ~2-5 s here; on the station NUC ~5-15 s for a lidar DSM "
          "(mostly reading the GeoTIFF), + ~3 s per 100k waterline points with --contours")
    compare(args.dem, args.survey, st, args.label, args.why, name=args.name,
            output_dir=args.output_dir, survey_date=args.survey_date,
            photo_dates=tuple(args.photo_dates) if args.photo_dates else None,
            contours=args.contours, camera_eo=parse_cam_specs(args.camera_eo, "--camera-eo"),
            camera_io=parse_cam_specs(args.camera_io, "--camera-io"), spread=args.spread,
            count=args.count, envelope_source=args.envelope_source, min_cover=args.min_cover,
            tolerance=args.tolerance, epsg=args.epsg, plot=not args.no_plot,
            transect_tolerance=args.transect_tolerance, transect_split=args.transect_split,
            transect_min_points=args.transect_min_points)
    return 0


if __name__ == "__main__":
    sys.exit(main())

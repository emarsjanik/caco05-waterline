#!/usr/bin/env python3
"""
Wave Runup From Argus Timestacks (ras.tiff)
==============================================
Tracks the landward edge of the swash through each 10-minute runup
timestack and reports where it sat on the beach: the mean swash
position, the 2%-exceedance runup position (R2%) and the swash
excursion, in metres along the timestack line.

WHY. The waterline detector works on timex images, where the bright
band is the MEAN swash position, and assigns it the still-water level.
The difference between the two is wave setup plus part of the swash,
which is the largest uncorrected term in the DEM. A timestack sees
every individual uprush, so it measures that difference instead of
estimating it from Stockdon et al. (2006).

INPUTS.
  * ras.tiff: rows are time (2 Hz, 1200 rows = 10 min), columns are
    the pixels listed in the camera's .pix file, in the same order.
  * .pix file (arguseyes/build/<cam>_timestack.pix): one "u v" image
    coordinate per timestack column. It holds several lines laid end
    to end; a jump of more than 100 px between consecutive points
    starts a new line. Lines are numbered from 1 in file order.
  * IO/EO calibration, to turn the detected pixel into ground
    coordinates.
  * Optionally the GNSS-R spline (still-water level at capture time)
    and the intertidal DEM (to read an elevation at the detected
    positions).

HOW THE EDGE IS FOUND. Each column's bare-sand brightness (its 3rd
percentile over the 10 minutes) is subtracted from it, which removes the
strong seaward-to-landward brightness gradient along the line; in the result every uprush
is a bright tongue of foam with a sharp landward tip, on uniform sand.
Whole-row brightness changes (exposure, passing cloud) are removed
first, using the landward end of the line as reference. A pixel counts
as foam when it is brighter than its column's bare-sand level by more than --k
times the noise measured on the landward 15% of the line. The edge in
each time row is the landward end of the most landward run of at least
--min-run foam pixels, so isolated bright specks do not count. A
5-sample (2.5 s) running median then removes single-row flicker.

An earlier version compared every column against one dry-sand level
taken from the landward end. On real stacks the sand itself brightens
steadily toward the sea, so that marked bare sand as foam and put the
edge well landward of the true tips.

WHICH END IS THE SEA is decided from the data: the end with more
temporal variability is the swash end. It is reported as sea_at so a
wrong call is visible.

WHAT IT CANNOT DO.
  * Ground coordinates are computed on the horizontal plane at the
    still-water level. The swash sits above it, so positions are
    slightly biased along the camera ray. Small for position, but it
    is one more reason elevations here are approximate.
  * Elevations (z_mean, z_r2, setup, R2) come from the DEM, which is
    itself built from waterlines biased by setup. Treat them as
    relative until the DEM is anchored to an independent survey.
  * A line that does not cross the swash at the time (e.g. C2 lines 3
    and 4 except at high water) gives a low valid_fraction; results
    below --min-valid are flagged, not trusted.

Usage:
    python3 runup_from_timestack.py /mnt/I2Rgus_Data/ImageProducts/products/*c2.ras.tiff \\
        --line 1 --gnssr /home/argus_user/GNSS/.../usgs_spline_out.txt \\
        --dem dem_intertidal_dem.asc --output runup_c2_line1.csv --plot-dir runup_plots
"""

import re
import sys
import csv
import argparse
import warnings
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import cv2

import georectify

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_PIX_DIR = "/home/argus_user/arguseyes/build"
SAMPLE_HZ = 2.0
LINE_BREAK_PX = 100.0
IG_CUTOFF_HZ = 0.05          # infragravity / incident split, standard in runup work
BURST_MID_OFFSET_S = 300.0   # filename epoch is the burst start; stats refer to its middle


def camera_from_name(name):
    m = re.search(r"\.(c\d)\.", name)
    return m.group(1) if m else None


def epoch_from_name(name):
    m = re.match(r"^(\d{9,10})\.", name)
    return int(m.group(1)) if m else None


def split_lines(pix):
    """Returns (start, end) column ranges, one per line in the .pix file."""
    jumps = np.flatnonzero(np.hypot(*np.diff(pix, axis=0).T) > LINE_BREAK_PX) + 1
    starts = np.r_[0, jumps]
    ends = np.r_[jumps, len(pix)]
    return list(zip(starts, ends))


def find_edge(gray, k, min_threshold, min_run=15):
    """
    Landward swash edge per time row, as a column index counted from the
    SEA end (column 0 = sea). NaN where no foam tongue is found in that
    row. Also returns the threshold used and the foam-anomaly image (for
    the diagnostic plot).
    """
    g = cv2.GaussianBlur(gray.astype(np.float32), (3, 3), 0)
    n = g.shape[1]
    land = slice(int(0.85 * n), n)

    # Whole-row brightness changes: exposure steps or cloud shadow.
    land_level = np.median(g[:, land], axis=1, keepdims=True)
    g = g - (land_level - np.median(land_level))

    # Bare-sand reference per column: its 3rd-percentile brightness.
    # The median would not do -- in the inner swash a column is under
    # foam more than half the time, so its median IS foam and foam there
    # would never stand out. The 3rd percentile stays on sand for any
    # column exposed at least ~3% of the time (36 of 1200 rows); on a
    # synthetic stack it halved the large errors of a 10th percentile
    # with no extra false detections. In pure noise it sits 1.88 sigma
    # below the median, hence the offset on the threshold.
    anomaly = g - np.percentile(g, 3, axis=0)
    land_part = g[:, land]
    noise = 1.4826 * np.median(np.abs(land_part - np.median(land_part, axis=0)))
    threshold = max(k * noise, min_threshold) + 1.88 * noise
    foam = anomaly > threshold

    # Length of the run of consecutive foam pixels ending at each column.
    run = np.zeros(foam.shape, np.int32)
    run[:, 0] = foam[:, 0]
    for x in range(1, n):
        run[:, x] = np.where(foam[:, x], run[:, x - 1] + 1, 0)
    ok = run >= min_run

    has = ok.any(axis=1)
    last = n - 1 - np.argmax(ok[:, ::-1], axis=1)
    edge = np.where(has, last, np.nan).astype(float)

    # 5-sample running median, ignoring gaps.
    padded = np.pad(edge, 2, mode="edge")
    stack = np.stack([padded[i:i + len(edge)] for i in range(5)])
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        smoothed = np.nanmedian(stack, axis=0)
    smoothed[~has] = np.nan
    return smoothed, threshold, anomaly


def runup_peaks(series, half_width=4):
    """Local maxima (one per uprush) of a runup series, NaNs ignored."""
    s = np.where(np.isfinite(series), series, -np.inf)
    padded = np.pad(s, half_width, constant_values=-np.inf)
    window = np.stack([padded[i:i + len(s)] for i in range(2 * half_width + 1)])
    is_peak = (s == window.max(axis=0)) & np.isfinite(series)
    # Collapse plateaus to their first sample.
    is_peak[1:] &= ~(is_peak[:-1] & (s[1:] == s[:-1]))
    return series[is_peak]


def band_sig(series, lo_hz, hi_hz):
    """Significant excursion (4 x std) of the part of the signal in [lo, hi) Hz."""
    x = series[np.isfinite(series)]
    if len(x) < 64:
        return np.nan
    x = x - x.mean()
    f = np.fft.rfftfreq(len(x), 1.0 / SAMPLE_HZ)
    X = np.fft.rfft(x)
    X[(f < lo_hz) | (f >= hi_hz)] = 0
    return 4.0 * np.std(np.fft.irfft(X, len(x)))


def read_dem(path):
    grid, hdr = None, {}
    with open(path) as f:
        for _ in range(6):
            key, val = f.readline().split()
            hdr[key.lower()] = float(val)
    grid = np.loadtxt(path, skiprows=6)
    grid = np.where(grid == hdr.get("nodata_value", -9999.0), np.nan, grid)
    return np.flipud(grid), hdr


def dem_at(dem, e, n):
    if dem is None or not (np.isfinite(e) and np.isfinite(n)):
        return np.nan
    grid, hdr = dem
    c = int((e - hdr["xllcorner"]) // hdr["cellsize"])
    r = int((n - hdr["yllcorner"]) // hdr["cellsize"])
    if 0 <= r < grid.shape[0] and 0 <= c < grid.shape[1]:
        return float(grid[r, c])
    return np.nan


def save_plot(path, anomaly, edge, threshold):
    """
    Left: the foam anomaly (each column minus its median), contrast-
    stretched so the uprush tongues stand out. Right: the same with the
    detected edge as a thin red line. The left panel is never drawn on,
    so the tips can always be checked against the line.
    """
    scale = 128.0 / max(4.0 * threshold, 1.0)
    grey = np.clip(128 + anomaly * scale, 0, 255).astype(np.uint8)
    left = cv2.cvtColor(grey, cv2.COLOR_GRAY2BGR)
    right = left.copy()
    pts = [(int(round(x)), t) for t, x in enumerate(edge) if np.isfinite(x)]
    for (x0, t0), (x1, t1) in zip(pts[:-1], pts[1:]):
        if t1 - t0 == 1:
            cv2.line(right, (x0, t0), (x1, t1), (0, 0, 255), 1)
    gap = np.zeros((left.shape[0], 8, 3), np.uint8)
    img = np.hstack([left, gap, right])
    w = 1400
    h = int(w * img.shape[0] / img.shape[1])
    cv2.imwrite(str(path), cv2.resize(img, (w, max(h, 1)), interpolation=cv2.INTER_AREA))


def process(ras_path, args, pix_cache, gnssr, dem):
    name = Path(ras_path).name
    cam = args.camera or camera_from_name(name)
    epoch = epoch_from_name(name)
    if cam is None or epoch is None:
        print(f"  {name}: cannot read camera/epoch from filename -- skipped")
        return None

    if cam not in pix_cache:
        pix = np.loadtxt(args.pix or f"{DEFAULT_PIX_DIR}/{cam}_timestack.pix")[:, :2]
        io = georectify.load_intrinsics(args.io or SCRIPT_DIR / f"calibration/CACO05_{cam}_20240801_IO.yaml")
        eo = georectify.load_extrinsics(args.eo or SCRIPT_DIR / f"calibration/CACO05_{cam}_20251113_EO-CV.yaml")
        pix_cache[cam] = (pix, split_lines(pix), io, eo)
    pix, lines, io, eo = pix_cache[cam]

    if not 1 <= args.line <= len(lines):
        sys.exit(f"{cam} has {len(lines)} line(s); --line {args.line} does not exist")
    a, b = lines[args.line - 1]

    ras = cv2.imread(str(ras_path), cv2.IMREAD_UNCHANGED)
    if ras is None:
        print(f"  {name}: unreadable -- skipped")
        return None
    if ras.shape[1] != len(pix):
        print(f"  {name}: {ras.shape[1]} columns but .pix has {len(pix)} -- skipped")
        return None

    seg = ras[:, a:b]
    gray = (cv2.cvtColor(seg, cv2.COLOR_BGR2GRAY) if seg.ndim == 3 else seg).astype(np.float32)
    line_uv = pix[a:b]

    # Sea end = the end with more temporal variability.
    std = gray.std(axis=0)
    n10 = max(len(std) // 7, 1)
    sea_at = "start" if std[:n10].mean() >= std[-n10:].mean() else "end"
    if sea_at == "end":
        gray, seg, line_uv = gray[:, ::-1], seg[:, ::-1], line_uv[::-1]

    edge, threshold, anomaly = find_edge(gray, args.k, args.min_threshold, args.min_run)
    valid = float(np.isfinite(edge).mean())

    mid_epoch = epoch + BURST_MID_OFFSET_S
    wl = np.nan
    if gnssr is not None:
        g_ep, g_lv = gnssr
        if g_ep[0] <= mid_epoch <= g_ep[-1]:
            wl = float(np.interp(mid_epoch, g_ep, g_lv))
    z_plane = wl if np.isfinite(wl) else args.z_plane

    # Ground coordinates of every pixel on the line, and distance along
    # it from the sea end.
    E, N = georectify.pixel_to_ground(line_uv[:, 0], line_uv[:, 1], z_plane, io, eo)
    step = np.hypot(np.diff(E), np.diff(N))
    dist = np.r_[0.0, np.nancumsum(step)]

    idx = np.arange(len(dist))
    pos = np.interp(edge, idx, dist)          # NaN edge stays NaN
    pos[~np.isfinite(edge)] = np.nan

    peaks = runup_peaks(pos)
    enough = len(peaks) >= args.min_peaks and valid >= args.min_valid
    mean_pos = float(np.nanmean(pos)) if np.isfinite(pos).any() else np.nan
    r2_pos = float(np.percentile(peaks, 98)) if enough else np.nan

    def at(p):
        if not np.isfinite(p):
            return np.nan, np.nan
        return float(np.interp(p, dist, E)), float(np.interp(p, dist, N))

    mE, mN = at(mean_pos)
    rE, rN = at(r2_pos)
    z_mean, z_r2 = dem_at(dem, mE, mN), dem_at(dem, rE, rN)

    if args.plot_dir:
        Path(args.plot_dir).mkdir(parents=True, exist_ok=True)
        save_plot(Path(args.plot_dir) / (Path(name).stem + f".line{args.line}.png"),
                  anomaly, edge, threshold)

    return {
        "source_file": name,
        "camera": cam,
        "line": args.line,
        "capture_epoch": epoch,
        "capture_time_utc": datetime.fromtimestamp(mid_epoch, tz=timezone.utc).isoformat(),
        "water_level_navd88": round(wl, 4) if np.isfinite(wl) else "",
        "sea_at": sea_at,
        "threshold": round(threshold, 2),
        "valid_fraction": round(valid, 3),
        "n_peaks": len(peaks),
        "trusted": int(enough),
        "mean_pos_m": round(mean_pos, 2),
        "sd_pos_m": round(float(np.nanstd(pos)), 2),
        "r2_pos_m": round(r2_pos, 2),
        "max_pos_m": round(float(np.nanmax(pos)), 2) if np.isfinite(pos).any() else np.nan,
        "swash_sig_m": round(band_sig(pos, 0.0, SAMPLE_HZ), 2),
        "swash_ig_sig_m": round(band_sig(pos, 1e-9, IG_CUTOFF_HZ), 2),
        "swash_inc_sig_m": round(band_sig(pos, IG_CUTOFF_HZ, SAMPLE_HZ), 2),
        "mean_easting": round(mE, 2), "mean_northing": round(mN, 2),
        "r2_easting": round(rE, 2), "r2_northing": round(rN, 2),
        "z_mean_dem": round(z_mean, 3), "z_r2_dem": round(z_r2, 3),
        "setup_m": round(z_mean - wl, 3) if np.isfinite(z_mean) and np.isfinite(wl) else "",
        "r2_m": round(z_r2 - wl, 3) if np.isfinite(z_r2) and np.isfinite(wl) else "",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("ras", nargs="+", help="One or more *.ras.tiff timestacks.")
    ap.add_argument("--line", type=int, default=1,
                    help="Which line of the .pix file to analyse (1-based, default 1). "
                         "For C2, lines 1 and 2 cross the swash; 3 and 4 only at high water.")
    ap.add_argument("--camera", choices=["c1", "c2"], help="Override camera from filename.")
    ap.add_argument("--pix", help=f"Default {DEFAULT_PIX_DIR}/<cam>_timestack.pix")
    ap.add_argument("--io", help="Default calibration/CACO05_<cam>_20240801_IO.yaml")
    ap.add_argument("--eo", help="Default calibration/CACO05_<cam>_20251113_EO-CV.yaml")
    ap.add_argument("--gnssr", help="gnssrefl spline file for the still-water level.")
    ap.add_argument("--z-plane", type=float, default=0.0,
                    help="Plane elevation (m NAVD88) for georectifying when no GNSS-R "
                         "value is available (default 0).")
    ap.add_argument("--dem", help="Intertidal DEM (.asc) to read elevations from.")
    ap.add_argument("--k", type=float, default=4.0,
                    help="Foam threshold in multiples of the dry-sand noise (default 4).")
    ap.add_argument("--min-threshold", type=float, default=4.0,
                    help="Floor on the foam threshold, in grey levels (default 4).")
    ap.add_argument("--min-run", type=int, default=15,
                    help="Minimum length, in pixels, of a foam tongue (default 15).")
    ap.add_argument("--min-valid", type=float, default=0.5,
                    help="Minimum fraction of time rows with a detected edge (default 0.5).")
    ap.add_argument("--min-peaks", type=int, default=20,
                    help="Minimum number of uprush peaks for an R2%% value (default 20).")
    ap.add_argument("--output", required=True, help="CSV to append results to.")
    ap.add_argument("--plot-dir", help="Write a diagnostic PNG per stack with the edge in red.")
    args = ap.parse_args()

    gnssr = None
    if args.gnssr:
        from extract_elevation_contours import load_gnssr_spline
        g_ep, g_lv, _, _ = load_gnssr_spline(args.gnssr)
        order = np.argsort(g_ep)
        gnssr = (np.asarray(g_ep, float)[order], np.asarray(g_lv, float)[order])
    dem = read_dem(args.dem) if args.dem else None

    out = Path(args.output)
    done = set()
    if out.exists():
        with open(out, newline="") as f:
            done = {(r["source_file"], r["line"]) for r in csv.DictReader(f)}

    pix_cache, rows = {}, []
    for ras_path in sorted(args.ras):
        if (Path(ras_path).name, str(args.line)) in done:
            continue
        row = process(ras_path, args, pix_cache, gnssr, dem)
        if row is None:
            continue
        rows.append(row)
        flag = "" if row["trusted"] else "  (not trusted)"
        print(f"  {row['source_file']}: mean {row['mean_pos_m']} m, R2% {row['r2_pos_m']} m, "
              f"valid {row['valid_fraction']:.0%}, {row['n_peaks']} uprushes{flag}")

    if not rows:
        print("No new timestacks processed.")
        return
    write_header = not out.exists()
    with open(out, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        if write_header:
            w.writeheader()
        w.writerows(rows)
    print(f"Appended {len(rows)} row(s) to {out}")


if __name__ == "__main__":
    main()

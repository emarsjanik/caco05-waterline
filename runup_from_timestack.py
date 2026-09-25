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
  * Optionally the GNSS-R spline (still-water level at capture time),
    used as the plane the positions are georectified onto.

HOW THE EDGE IS FOUND. Each column's bare-sand brightness (its 3rd
percentile within ~100 s blocks, split also wherever the light on the
water changes suddenly, so such a change partway through the stack
does not turn everything after it into "foam")
is subtracted from it, which removes the
strong seaward-to-landward brightness gradient along the line; in the result every uprush
is a bright tongue of foam with a sharp landward tip, on uniform sand.
Whole-row brightness changes (exposure, passing cloud) are removed
first, using the landward end of the line as reference. A pixel counts
as foam when it is brighter than its column's bare-sand level by more than --k
times the noise measured on the landward 15% of the line. The edge in
each time row is the landward end of the most landward run of at least
--min-run foam pixels, so isolated bright specks do not count. A
5-sample (2.5 s) running median then removes single-row flicker.

Backwash is thin, clear water with little foam, so in some rows the
foam test briefly loses it and the edge would jump most of the way to
the sea for half a second. No real backwash moves that fast, so the
edge may retreat by at most --max-retreat-speed (default 3 m/s)
between rows; uprush is not limited.

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
  * NO ELEVATIONS. An earlier version read setup and R2 off the
    intertidal DEM, but that DEM is built by giving each timex
    waterline -- which sits at the mean swash position -- the
    still-water level. At the mean swash position the DEM therefore
    equals still water by construction and "setup" came out near zero
    whatever the waves did. Elevations need a beach profile measured
    independently of the waterlines (RTK, UAS or lidar survey).
  * Swash statistics are HORIZONTAL distances along the line (the
    *_horiz_m columns), not the vertical S of Stockdon et al. (2006).
    Converting needs the same independent profile.
  * A line that does not cross the swash at the time (e.g. C2 lines 3
    and 4 except at high water), or a stack with faint foam, gives a
    low valid_fraction; results below --min-valid are flagged, not
    trusted.

RUNNING IT ROUTINELY. The station's cleanup moves ras.tiff files off
the machine, while GNSS-R lags about 2 days, so waterline_timex_cron.sh
copies the C2 stacks into archive/ras_c2 and runs this script over that
archive with --require-water-level: each stack is processed once, on
the first run after GNSS-R covers it.

Usage:
    python3 runup_from_timestack.py /mnt/I2Rgus_Data/ImageProducts/products/*c2.ras.tiff \\
        --line 1 --gnssr /home/argus_user/GNSS/.../usgs_spline_out.txt \\
        --output runup_c2_line1.csv --plot-dir runup_plots
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


def light_steps(g, threshold, half=60, min_gap=60):
    """
    Rows where the brightness of the seaward half of the line jumps
    suddenly (light on the water changing), found by comparing the
    median level over `half` rows before and after each row. Returns
    the row indices of jumps larger than `threshold`.
    """
    sea = np.median(g[:, : g.shape[1] // 2], axis=1)
    T = len(sea)
    jump = np.zeros(T)
    for t in range(half, T - half):
        jump[t] = abs(np.median(sea[t:t + half]) - np.median(sea[t - half:t]))
    steps = []
    for t in np.argsort(jump)[::-1]:
        if jump[t] <= threshold:
            break
        if all(abs(t - u) >= min_gap for u in steps):
            steps.append(int(t))
    return sorted(steps)


def local_percentile(g, pct, block, steps=(), min_len=60):
    """
    Per-row baseline: the percentile of each column within time segments
    of about `block` rows, split additionally at `steps` (sudden light
    changes). Constant within a segment -- interpolating between
    segments smeared a sudden change over ~100 rows.
    """
    T = g.shape[0]
    nb = max(int(round(T / block)), 1)
    cuts = set(np.linspace(0, T, nb + 1).astype(int)[1:-1])
    cuts = sorted(c for c in cuts if all(abs(c - s) >= min_len for s in steps))
    cuts = sorted(set(cuts) | set(steps))
    bounds = [0] + [c for c in cuts if min_len <= c <= T - min_len] + [T]
    out = np.empty_like(g)
    for a, b in zip(bounds[:-1], bounds[1:]):
        out[a:b] = np.percentile(g[a:b], pct, axis=0)
    return out


def find_edge(gray, k, min_threshold, min_run=15, block=200):
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

    # Bare-sand reference per column: its 3rd-percentile brightness,
    # taken over ~100 s blocks rather than the whole stack -- a change in
    # the light on the water partway through otherwise made the whole
    # brighter half register as foam.
    # The median would not do -- in the inner swash a column is under
    # foam more than half the time, so its median IS foam and foam there
    # would never stand out. The 3rd percentile stays on sand for any
    # column exposed at least ~3% of the time (36 of 1200 rows); on a
    # synthetic stack it halved the large errors of a 10th percentile
    # with no extra false detections. In pure noise it sits 1.88 sigma
    # below the median, hence the offset on the threshold.
    land_part = g[:, land]
    noise = 1.4826 * np.median(np.abs(land_part - np.median(land_part, axis=0)))
    threshold = max(k * noise, min_threshold) + 1.88 * noise
    anomaly = g - local_percentile(g, 3, block, light_steps(g, threshold))
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


MAX_FILL_SAMPLES = 10    # gaps up to 5 s are bridged before spectral analysis


def band_sig(series, lo_hz, hi_hz):
    """
    Significant excursion (4 x std) of the part of the signal in
    [lo, hi) Hz. Gaps of up to MAX_FILL_SAMPLES are bridged linearly;
    with any longer gap the series is no longer an evenly sampled 2 Hz
    record and NaN is returned. (An earlier version dropped the NaNs and
    joined the remainder, which shifts every frequency.)
    """
    x = np.asarray(series, float)
    bad = ~np.isfinite(x)
    if bad.all() or bad.sum() > 0.1 * len(x):
        return np.nan
    if bad.any():
        edges = np.flatnonzero(np.diff(np.r_[0, bad.astype(int), 0]))
        if (edges[1::2] - edges[::2]).max() > MAX_FILL_SAMPLES:
            return np.nan
        idx = np.arange(len(x))
        x = np.interp(idx, idx[~bad], x[~bad])
    if len(x) < 64:
        return np.nan
    x = x - x.mean()
    f = np.fft.rfftfreq(len(x), 1.0 / SAMPLE_HZ)
    X = np.fft.rfft(x)
    X[(f < lo_hz) | (f >= hi_hz)] = 0
    return 4.0 * np.std(np.fft.irfft(X, len(x)))


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


WL_MAX_GAP_S = 3600.0    # same rule as extract_elevation_contours.py --max-gap-minutes


def water_level_at(gnssr, epoch):
    """
    GNSS-R level at `epoch`, linearly interpolated -- or NaN if `epoch`
    is outside the record or either bracketing reading is more than
    WL_MAX_GAP_S away (i.e. the line would be drawn across an outage).
    """
    g_ep, g_lv = gnssr
    if not g_ep[0] <= epoch <= g_ep[-1]:
        return np.nan
    i = int(np.searchsorted(g_ep, epoch))
    if g_ep[i] == epoch:
        return float(g_lv[i])
    if max(epoch - g_ep[i - 1], g_ep[i] - epoch) > WL_MAX_GAP_S:
        return np.nan
    return float(np.interp(epoch, g_ep, g_lv))


def process(ras_path, args, pix_cache, gnssr):
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
    wl = water_level_at(gnssr, mid_epoch) if gnssr is not None else np.nan
    z_plane = wl if np.isfinite(wl) else args.z_plane

    # Ground coordinates of every pixel on the line, and distance along
    # it from the sea end.
    E, N = georectify.pixel_to_ground(line_uv[:, 0], line_uv[:, 1], z_plane, io, eo)
    step = np.hypot(np.diff(E), np.diff(N))
    dist = np.r_[0.0, np.nancumsum(step)]

    idx = np.arange(len(dist))
    pos = np.interp(edge, idx, dist)          # NaN edge stays NaN
    pos[~np.isfinite(edge)] = np.nan

    # Limit how fast the edge may retreat seaward (see docstring).
    step_m = args.max_retreat_speed / SAMPLE_HZ
    for t in range(1, len(pos)):
        if np.isfinite(pos[t - 1]) and np.isfinite(pos[t]) and pos[t] < pos[t - 1] - step_m:
            pos[t] = pos[t - 1] - step_m
    edge = np.where(np.isfinite(pos), np.interp(pos, dist, idx), np.nan)

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
        "swash_sig_horiz_m": round(band_sig(pos, 0.0, SAMPLE_HZ), 2),
        "swash_ig_sig_horiz_m": round(band_sig(pos, 1e-9, IG_CUTOFF_HZ), 2),
        "swash_inc_sig_horiz_m": round(band_sig(pos, IG_CUTOFF_HZ, SAMPLE_HZ), 2),
        "mean_easting": round(mE, 2), "mean_northing": round(mN, 2),
        "r2_easting": round(rE, 2), "r2_northing": round(rN, 2),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("ras", nargs="+",
                    help="*.ras.tiff timestacks, or folders of them (all *.ras.tiff inside "
                         "matching --camera, or the camera in the filename).")
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
    ap.add_argument("--k", type=float, default=4.0,
                    help="Foam threshold in multiples of the dry-sand noise (default 4).")
    ap.add_argument("--min-threshold", type=float, default=4.0,
                    help="Floor on the foam threshold, in grey levels (default 4).")
    ap.add_argument("--min-run", type=int, default=15,
                    help="Minimum length, in pixels, of a foam tongue (default 15).")
    ap.add_argument("--max-retreat-speed", type=float, default=3.0,
                    help="Fastest the swash edge may move seaward, m/s (default 3). "
                         "Stops single-row drop-outs during backwash reading as retreats "
                         "of tens of metres.")
    ap.add_argument("--min-valid", type=float, default=0.8,
                    help="Minimum fraction of time rows with a detected edge for a stack to "
                         "be trusted (default 0.8). Stacks with faint foam came in at 65-67%% "
                         "and gave fragmented edges; clear ones at 92-100%%.")
    ap.add_argument("--min-peaks", type=int, default=20,
                    help="Minimum number of uprush peaks for an R2%% value (default 20).")
    ap.add_argument("--require-water-level", action="store_true",
                    help="Skip (without recording) stacks that --gnssr does not yet cover, so "
                         "a later run processes them once GNSS-R catches up (~2 days).")
    ap.add_argument("--reprocess", action="store_true",
                    help="Redo stacks already in --output, replacing their rows. Use after "
                         "the method or its settings change.")
    ap.add_argument("--output", required=True,
                    help="CSV of results. Stacks already in it are skipped unless --reprocess.")
    ap.add_argument("--plot-dir", help="Write a diagnostic PNG per stack with the edge in red.")
    args = ap.parse_args()

    gnssr = None
    if args.gnssr:
        from extract_elevation_contours import load_gnssr_spline
        g_ep, g_lv, _, _ = load_gnssr_spline(args.gnssr)
        order = np.argsort(g_ep)
        gnssr = (np.asarray(g_ep, float)[order], np.asarray(g_lv, float)[order])

    out = Path(args.output)
    existing = []
    if out.exists():
        with open(out, newline="") as f:
            existing = list(csv.DictReader(f))
    done = {(r["source_file"], r["line"]) for r in existing}

    paths = []
    for item in args.ras:
        item = Path(item)
        if item.is_dir():
            pattern = f"*.{args.camera}.ras.tiff" if args.camera else "*.ras.tiff"
            paths.extend(item.glob(pattern))
        else:
            paths.append(item)

    pix_cache, rows, waiting = {}, [], 0
    for ras_path in sorted(paths, key=lambda q: q.name):
        name = Path(ras_path).name
        if (name, str(args.line)) in done and not args.reprocess:
            continue
        if args.require_water_level and gnssr is not None:
            epoch = epoch_from_name(name)
            mid = (epoch or 0) + BURST_MID_OFFSET_S
            if epoch is None or not gnssr[0][0] <= mid <= gnssr[0][-1]:
                waiting += 1
                continue
        row = process(ras_path, args, pix_cache, gnssr)
        if row is None:
            continue
        rows.append(row)
        flag = "" if row["trusted"] else "  (not trusted)"
        print(f"  {row['source_file']}: mean {row['mean_pos_m']} m, R2% {row['r2_pos_m']} m, "
              f"valid {row['valid_fraction']:.0%}, {row['n_peaks']} uprushes{flag}")

    if waiting:
        print(f"{waiting} stack(s) not yet covered by GNSS-R -- left for a later run.")
    if not rows:
        print("No new timestacks processed.")
        return

    # Rewrite the whole file: kept rows plus new ones, in time order.
    # Rows being redone (--reprocess) are replaced, never duplicated.
    redone = {(r["source_file"], str(r["line"])) for r in rows}
    kept = [r for r in existing if (r["source_file"], r["line"]) not in redone]
    merged = kept + rows
    merged.sort(key=lambda r: (int(r["line"]), int(r["capture_epoch"])))
    fields = list(rows[0].keys())
    tmp = out.with_suffix(out.suffix + ".tmp")
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(merged)
    tmp.replace(out)
    print(f"Wrote {len(rows)} new row(s) to {out} ({len(merged)} total)")


if __name__ == "__main__":
    main()

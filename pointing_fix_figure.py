#!/usr/bin/env python3
"""
Camera Pointing Fix: One-Page Before And After
=================================================
Shows, on one page, what the camera pointing correction changed and why:
the same waterline pixels placed on the ground with the OLD pointing and
with the FITTED pointing, against a lidar survey of the same days.

THE IDEA IN ONE SENTENCE. A waterline is a contour of known elevation (the
water level when the photo was taken); the camera's pointing decides where
on the beach that contour lands. With a wrong pointing every contour is
drawn in the wrong place, so on a sloping beach it is given the wrong
ground elevation: vertical error ~ beach slope x horizontal misplacement
(1:20 beach, 30 m -> 1.5 m).

WHAT IS COMPUTED. Every waterline point of the period (up to --per-frame
per photo) is traced from its pixel to the ground at its water level twice
-- once with the old pointing, once with the fitted one -- with exactly the
geometry fit_eo_to_survey.py used (same lens, same camera position, same
undoing of the resampling into today's view). At each ground position the
lidar elevation is read; a correctly placed waterline point has

    residual = lidar elevation - water level = 0.

Panels:
  A  map: the lidar (grey), the waterlines with the old pointing (orange)
     and the fitted one (blue), one photo's line drawn bold;
  B  distribution of the residuals, old and fitted, with median, NMAD and
     the share of points that land on the survey at all;
  C  residual against distance from the camera (median and middle half per
     band) -- a pointing error grows with range, a datum error does not;
  D  how far a small pointing error moves a waterline, by distance from the
     camera, for this camera (cross-shore metres; vertical = slope x shift).

HONEST LABEL. If the fitted pointing was fitted TO this survey (its EO file
says so), the "after" agreement shows that a pointing error explains the
offset -- it is a consistency check, not an independent accuracy number,
and the page says so. Run it against a survey the pointing was NOT fitted
to (e.g. another period, or RTK) for an independent number.

INPUTS. The EO file written by fit_eo_to_survey.py records the pointing it
started from, the survey and the dates, so by default nothing else is
needed:

    # c2, Jan 18-23 2025, lines detected in frames redrawn into today's view
    python3 pointing_fix_figure.py /mnt/I2Rgus_Data/Chelsea_calibration/contour_points_ground.csv \\
        --survey /mnt/I2Rgus_Data/Chelsea_calibration/2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif \\
        --after-eo calibration/CACO05_c2_2025-01-18_to_2025-01-23_lidar_EO.yaml

    # c2, Mar 2025, lines detected in the ORIGINAL photos (detect_original_view.py)
    python3 pointing_fix_figure.py /mnt/I2Rgus_Data/Chelsea_calibration/original_view_c2_mar/contour_points_ground.csv \\
        --survey /mnt/I2Rgus_Data/Chelsea_calibration/2025005FA_Marconi_Mar_YSMP_Lidar_DSM_25cm.tif \\
        --after-eo calibration/CACO05_c2_2025-03-01_to_2025-03-10_lidar_EO.yaml --original-view

    # independent: the station's own calibration of the old setup (calibration/README.md)
    python3 pointing_fix_figure.py /mnt/I2Rgus_Data/Chelsea_calibration/original_view_c2/contour_points_ground.csv \\
        --survey /mnt/I2Rgus_Data/Chelsea_calibration/2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif \\
        --before-eo calibration/CACO05_c2_2025-01-08_to_2025-01-23_EO.yaml \\
        --after-eo calibration/CACO03_c2_20250123_EO.yaml --original-view \\
        --camera c2 --start-date 2025-01-18 --end-date 2025-01-23
    The second pointing is then named 'CACO03 calibration of 2025-01-23' and the
    page says it was not fitted to the survey; a camera move is reported too.

--before-eo overrides the starting pointing; --resampled-with names the
pointing the frames were redrawn with, if the contour file was made after
chelsea_setups.csv already pointed at the fitted EO (then the lines were
redrawn with the fitted pointing and that is what must be undone).
Writes <output>.png and <output>.csv (the numbers on the page).
"""

import re
import sys
import csv
import argparse
from pathlib import Path
from collections import defaultdict

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from georectify import load_extrinsics, load_intrinsics, pixel_to_ground   # noqa: E402
from view_reproject import dst_to_src_points                              # noqa: E402
from compare_dem_survey import read_survey, sample                        # noqa: E402

CAL = HERE / "calibration"
EO_NEW = "20251113"                          # today's view, as in fit_eo_to_survey.py
BANDS = [(0, 100), (100, 200), (200, 300), (300, 450), (450, 700)]
BEFORE, AFTER = "#eb6834", "#2a78d6"         # orange = old pointing, blue = fitted
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dc"
SHORE_NORMAL_DEG = 80.0                      # Marconi faces ~80 deg: cross-shore direction


def eo_note(path):
    """The '# ...' note fit_eo_to_survey.py writes into an EO file, or ''."""
    for line in Path(path).read_text().splitlines():
        if line.startswith("# fit_eo_to_survey.py:"):
            return line[2:].strip()
    return ""


def parse_fit_note(note):
    """'fit_eo_to_survey.py: OLD.yaml with azimuth +a, tilt +t, roll +r deg fitted to SURVEY using
    N waterline frames D1 to D2; ...' -> dict (missing parts absent)."""
    out = {}
    m = re.search(r"fit_eo_to_survey\.py:\s*(\S+)\s+with", note)
    if m:
        out["before"] = m.group(1)
    m = re.search(r"fitted to (\S+)", note)
    if m:
        out["survey"] = m.group(1)
    m = re.search(r"frames (\d{4}-\d{2}-\d{2}) to (\d{4}-\d{2}-\d{2})", note)
    if m:
        out["dates"] = (m.group(1), m.group(2))
    return out


def angles_deg(eo):
    return np.degrees(np.asarray(eo, float)[3:6])


def stats(r):
    """Residual summary over points that landed on the survey."""
    on = np.isfinite(r)
    x = r[on]
    if not len(x):
        return {"n": 0, "on_pct": 0.0, "median": np.nan, "nmad": np.nan, "rmse": np.nan, "p90_abs": np.nan}
    med = float(np.median(x))
    return {"n": int(len(x)), "on_pct": 100.0 * on.mean(), "median": med,
            "nmad": float(1.4826 * np.median(np.abs(x - med))),
            "rmse": float(np.sqrt(np.mean(x ** 2))), "p90_abs": float(np.percentile(np.abs(x), 90))}


def sensitivity(io, eo, z):
    """Median cross-shore shift (m) of the ground point under each pixel, for small pointing errors,
    by distance band. -> {(name, deg): [shift per band]}"""
    U, V = np.meshgrid(np.arange(25.0, io[0], 50.0), np.arange(5.0, io[1], 10.0))
    U, V = U.ravel(), V.ravel()
    E0, N0 = pixel_to_ground(U, V, z, io, eo)
    rng = np.hypot(E0 - eo[0], N0 - eo[1])
    sn = np.array([np.sin(np.radians(SHORE_NORMAL_DEG)), np.cos(np.radians(SHORE_NORMAL_DEG))])
    table = {}
    for name, k in (("pan", 3), ("tilt", 4), ("roll", 5)):
        for d in (0.1, 1.0):
            e = np.array(eo, float).copy()
            e[k] += np.radians(d)
            E1, N1 = pixel_to_ground(U, V, z, io, e)
            cross = np.abs((E1 - E0) * sn[0] + (N1 - N0) * sn[1])
            row = []
            for a, b in BANDS:
                m = np.isfinite(rng) & np.isfinite(cross) & (rng >= a) & (rng < b)
                row.append(float(np.median(cross[m])) if m.sum() > 5 else np.nan)
            table[(name, d)] = row
    return table


def load_points(path, camera, d1, d2, per_frame):
    """-> U, V, Z, frame index, frame names, epochs (water level from beach_ or tide_elevation_navd88)."""
    frames = defaultdict(list)
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if r.get("camera") != camera:
                continue
            day = r.get("capture_time_utc", "")[:10]
            if (d1 and day < d1) or (d2 and day > d2):
                continue
            frames[r["source_file"]].append(r)
    U, V, Z, F, names = [], [], [], [], []
    for k, (name, rows) in enumerate(sorted(frames.items())):
        names.append(name)
        idx = np.linspace(0, len(rows) - 1, min(per_frame, len(rows))).astype(int)
        for i in idx:
            z = rows[i].get("beach_elevation_navd88") or rows[i]["tide_elevation_navd88"]
            U.append(float(rows[i]["pixel_column"])); V.append(float(rows[i]["pixel_row"]))
            Z.append(float(z)); F.append(k)
    return np.array(U), np.array(V), np.array(Z), np.array(F, int), names


def binned(dist, r):
    """Per distance band: (n, median, p25, p75) of the residuals on the survey."""
    out = []
    for a, b in BANDS:
        m = np.isfinite(r) & (dist >= a) & (dist < b)
        x = r[m]
        out.append((int(m.sum()),) + ((float(np.median(x)), float(np.percentile(x, 25)),
                                       float(np.percentile(x, 75))) if len(x) >= 5 else (np.nan,) * 3))
    return out


def figure(out_png, survey, cam_xy, P0, P1, r0, r1, dist, F, best_frame, s0, s1, b0, b1, sens,
           title, subtitle, label, after_label="fitted pointing"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    grid, x0, y0, cell = survey
    fig = plt.figure(figsize=(15, 10.5), dpi=110)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], height_ratios=[1, 1], hspace=0.38, wspace=0.18,
                          left=0.05, right=0.96, top=0.79, bottom=0.21)

    def style(ax):
        ax.grid(True, color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelsize=9)

    # --- A: map -------------------------------------------------------------------------------
    ax = fig.add_subplot(gs[:, 0])
    on = np.isfinite(r0) | np.isfinite(r1)
    allE = np.concatenate([P0[0][on], P1[0][on]]); allN = np.concatenate([P0[1][on], P1[1][on]])
    if len(allE) < 10:
        allE, allN = np.concatenate([P0[0], P1[0]]), np.concatenate([P0[1], P1[1]])
    ok = np.isfinite(allE) & np.isfinite(allN)
    e_lo, e_hi = np.percentile(allE[ok], [1, 99]); n_lo, n_hi = np.percentile(allN[ok], [1, 99])
    pad = 30.0
    e_lo, e_hi, n_lo, n_hi = e_lo - pad, e_hi + pad, n_lo - pad, n_hi + pad
    c_lo = int(max(0, np.floor((e_lo - x0) / cell))); c_hi = int(min(grid.shape[1], np.ceil((e_hi - x0) / cell)))
    r_lo = int(max(0, np.floor((y0 - n_hi) / cell))); r_hi = int(min(grid.shape[0], np.ceil((y0 - n_lo) / cell)))
    if c_hi > c_lo and r_hi > r_lo:
        sub = grid[r_lo:r_hi, c_lo:c_hi]
        step = max(1, int(np.ceil(max(sub.shape) / 1400)))
        sub = sub[::step, ::step]
        fin = sub[np.isfinite(sub)]
        if len(fin):
            vmin, vmax = np.percentile(fin, [2, 98])
            im = ax.imshow(sub, cmap="Greys_r", vmin=vmin, vmax=vmax, interpolation="nearest",
                           extent=(x0 + c_lo * cell, x0 + c_hi * cell, y0 - r_hi * cell, y0 - r_lo * cell))
            cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.01)
            cb.set_label("lidar elevation (m NAVD88)", color=MUTED, fontsize=9)
            cb.ax.tick_params(colors=MUTED, labelsize=8)
    ax.scatter(P0[0], P0[1], s=2, color=BEFORE, alpha=0.25, lw=0, rasterized=True)
    ax.scatter(P1[0], P1[1], s=2, color=AFTER, alpha=0.25, lw=0, rasterized=True)
    if best_frame is not None:
        m = F == best_frame
        ax.plot(P0[0][m], P0[1][m], ".", ms=4, color=BEFORE)
        ax.plot(P1[0][m], P1[1][m], ".", ms=4, color=AFTER)
    ax.plot(*cam_xy, marker="^", ms=10, color=INK, ls="none")
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=BEFORE, marker="o", ls="none", ms=7, label="old pointing"),
               Line2D([], [], color=AFTER, marker="o", ls="none", ms=7, label=after_label),
               Line2D([], [], color=INK, marker="^", ls="none", ms=8, label="camera")]
    ax.set_xlim(min(e_lo, cam_xy[0] - 10), max(e_hi, cam_xy[0] + 10))
    ax.set_ylim(min(n_lo, cam_xy[1] - 10), max(n_hi, cam_xy[1] + 10))
    ax.set_aspect("equal")
    ax.set_xlabel("UTM 19N easting (m)", color=MUTED, fontsize=9)
    ax.set_ylabel("northing (m)", color=MUTED, fontsize=9)
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.legend(handles=handles, loc="upper left", frameon=True, framealpha=0.9, fontsize=9)
    ax.set_title("A  The waterlines on the lidar\n    bold: one photo; faint: all photos",
                 loc="left", fontsize=10.5, color=INK)
    style(ax)

    # --- B: residual distribution -------------------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    bins = np.arange(-3.0, 3.0001, 0.1)
    for r, s, col, name in ((r0, s0, BEFORE, "old pointing"), (r1, s1, AFTER, after_label)):
        x = np.clip(r[np.isfinite(r)], -3.0, 3.0)
        if len(x):
            ax.hist(x, bins=bins, histtype="step", lw=2, color=col,
                    label=f"{name}: median {s['median']:+.2f} m, NMAD {s['nmad']:.2f} m, "
                          f"{s['on_pct']:.0f}% on the survey")
    ax.axvline(0, color=INK, lw=1)
    ax.set_xlabel("lidar elevation - water level at the waterline (m)   [0 = in the right place]",
                  color=MUTED, fontsize=9)
    ax.set_ylabel("waterline points", color=MUTED, fontsize=9)
    ax.set_ylim(top=ax.get_ylim()[1] * 1.45)          # room for the legend above the bars
    ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor=GRID, framealpha=1, fontsize=8.5)
    ax.set_title("B  How far off the lines are, vertically", loc="left", fontsize=10.5, color=INK)
    style(ax)

    # --- C: by distance -----------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    centres = np.array([(a + b) / 2 for a, b in BANDS])
    for rows, col, name, dx in ((b0, BEFORE, "old pointing", -6), (b1, AFTER, after_label, 6)):
        n = np.array([x[0] for x in rows]); med = np.array([x[1] for x in rows])
        lo = np.array([x[2] for x in rows]); hi = np.array([x[3] for x in rows])
        good = np.isfinite(med)
        if good.any():
            ax.errorbar(centres[good] + dx, med[good], yerr=[med[good] - lo[good], hi[good] - med[good]],
                        fmt="o", ms=7, color=col, lw=2, capsize=0, label=name)
            ax.plot(centres[good] + dx, med[good], color=col, lw=2, alpha=0.6)
    ax.axhline(0, color=INK, lw=1)
    ax.set_xticks(centres)
    ax.set_xticklabels([f"{a}-{b} m" for a, b in BANDS])
    ax.set_xlabel("distance from the camera (dot = median, bar = middle half of points)", color=MUTED, fontsize=9)
    ax.set_ylabel("lidar - water level (m)", color=MUTED, fontsize=9)
    ax.legend(loc="best", frameon=False, fontsize=9)
    ax.set_title("C  A pointing error grows with distance; a datum error would not", loc="left",
                 fontsize=10.5, color=INK)
    style(ax)

    # --- D: sensitivity table, under panel C ------------------------------------------------
    rows = [f"{nm} +{d:g} deg" for (nm, d) in sens]
    cells = [[("-" if not np.isfinite(v) else f"{v:.1f}") for v in sens[k]] for k in sens]
    txt = ("D  How far a small pointing error moves a waterline (cross-shore metres), this camera\n"
           + f"{'':13s}" + "".join(f"{f'{a}-{b} m':>10s}" for a, b in BANDS) + "\n"
           + "\n".join(f"{r:13s}" + "".join(f"{c:>10s}" for c in row) for r, row in zip(rows, cells))
           + "\nElevation error = beach slope x shift (on a 1:20 beach, 10 m -> 0.5 m).")
    fig.text(0.535, 0.015, txt, family="monospace", fontsize=8.3, color=INK, va="bottom")

    fig.suptitle(title, x=0.05, y=0.985, ha="left", fontsize=14, color=INK, weight="bold")
    fig.text(0.05, 0.94, subtitle, ha="left", va="top", fontsize=9.5, color=MUTED)
    fig.text(0.05, 0.868, label, ha="left", va="top", fontsize=9.5,
             color=INK, bbox=dict(boxstyle="round,pad=0.35", fc="#f3f2ec", ec=GRID))
    fig.savefig(out_png)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("contours", help="waterline points (process_chelsea.py / detect_original_view.py CSV)")
    ap.add_argument("--survey", required=True, help="lidar DSM (.tif or .asc), NAVD88")
    ap.add_argument("--after-eo", required=True,
                    help="the new pointing: a fit_eo_to_survey.py output, or a station calibration "
                         "(e.g. CACO03_c2_20250123_EO.yaml)")
    ap.add_argument("--after-label", default=None,
                    help="what to call it on the page (default: 'fitted pointing' for a survey fit, "
                         "'<setup> calibration of <date>' for a station calibration)")
    ap.add_argument("--before-eo", default=None,
                    help="the old pointing (default: the one named in --after-eo's note)")
    ap.add_argument("--resampled-with", default=None,
                    help="pointing the frames were redrawn into today's view with (default: --before-eo)")
    ap.add_argument("--original-view", action="store_true",
                    help="pixels are in the ORIGINAL photos (detect_original_view.py): nothing to undo")
    ap.add_argument("--camera", default=None, help="c1 or c2 (default: from --after-eo's name)")
    ap.add_argument("--start-date"); ap.add_argument("--end-date")
    ap.add_argument("--per-frame", type=int, default=200)
    ap.add_argument("--output", default=None, help="stem for .png and .csv (default pointing_fix_<cam>_<dates>)")
    args = ap.parse_args()

    after_path = Path(args.after_eo) if Path(args.after_eo).exists() else CAL / args.after_eo
    if not after_path.exists():
        sys.exit(f"{args.after_eo}: not found")
    note = parse_fit_note(eo_note(after_path))
    m = re.search(r"CACO05_(c\d)_(\d{4}-\d{2}-\d{2})_to_(\d{4}-\d{2}-\d{2})", after_path.name)
    mc = re.search(r"_(c\d)_", after_path.name)
    cam = args.camera or (mc.group(1) if mc else None)
    cal = re.match(r"^(CACO\d+)_c\d_(\d{4})(\d{2})(\d{2})_EO(?:-CV)?\.yaml$", after_path.name)
    after_label = args.after_label or (
        "fitted pointing" if note.get("survey") else
        f"{cal.group(1)} calibration of {cal.group(2)}-{cal.group(3)}-{cal.group(4)}" if cal else after_path.stem)
    if cam not in ("c1", "c2"):
        sys.exit("--camera c1 or c2 (it could not be read from the EO file name)")
    d1 = args.start_date or (m.group(2) if m else note.get("dates", (None, None))[0])
    d2 = args.end_date or (m.group(3) if m else note.get("dates", (None, None))[1])
    before_name = args.before_eo or note.get("before")
    if not before_name:
        sys.exit(f"{after_path.name} does not name the pointing it started from; give --before-eo")
    before_path = Path(before_name) if Path(before_name).exists() else CAL / before_name
    resampled_path = Path(args.resampled_with) if args.resampled_with else before_path
    if args.resampled_with and not resampled_path.exists():
        resampled_path = CAL / args.resampled_with
    for p in (before_path, resampled_path):
        if not p.exists():
            sys.exit(f"{p}: not found")

    io = load_intrinsics(CAL / f"CACO05_{cam}_20240801_IO.yaml")
    eo_before, eo_after = load_extrinsics(before_path), load_extrinsics(after_path)
    U, V, Z, F, names = load_points(args.contours, cam, d1, d2, args.per_frame)
    if not len(U):
        sys.exit(f"no {cam} waterline points between {d1} and {d2} in {args.contours}")

    # back to the original photo's pixels, exactly as fit_eo_to_survey.py does
    if args.original_view:
        u0, v0, ok = U, V, np.ones(len(U), bool)
    else:
        eo_res = load_extrinsics(resampled_path)
        eo_new = load_extrinsics(CAL / f"CACO05_{cam}_{EO_NEW}_EO-CV.yaml")
        u0, v0, ok = dst_to_src_points(U, V, Z, io, eo_res, eo_new)
    u0, v0, Z, F = u0[ok], v0[ok], Z[ok], F[ok]

    survey = read_survey(args.survey)
    E0, N0 = pixel_to_ground(u0, v0, Z, io, eo_before)
    E1, N1 = pixel_to_ground(u0, v0, Z, io, eo_after)
    r0 = sample(*survey, E0, N0) - Z
    r1 = sample(*survey, E1, N1) - Z
    dist = np.hypot(E1 - eo_after[0], N1 - eo_after[1])
    s0, s1 = stats(r0), stats(r1)
    b0, b1 = binned(dist, r0), binned(dist, r1)
    moved = np.hypot(E1 - E0, N1 - N0)
    sens = sensitivity(io, eo_after, float(np.median(Z)))

    counts = np.bincount(F[np.isfinite(r1)], minlength=len(names)) if len(F) else np.array([])
    best = int(np.argmax(counts)) if len(counts) and counts.max() > 0 else None

    a0, a1 = angles_deg(eo_before), angles_deg(eo_after)
    dA = a1 - a0
    survey_name = Path(args.survey).name
    fitted_to_this = note.get("survey") == survey_name or survey_name in eo_note(after_path)
    if fitted_to_this:
        label = (f"CONSISTENCY CHECK, NOT AN INDEPENDENT ACCURACY NUMBER: the {after_label} was fitted to "
                 f"this survey ({survey_name}).\nIt shows the offset is explained by the camera pointing. "
                 "For independent accuracy, compare against a survey the pointing was not fitted to.")
    else:
        label = (f"INDEPENDENT CHECK: the {after_label} was not fitted to {survey_name}"
                 + (f" (it was fitted to {note['survey']})." if note.get("survey") else
                    " (a station calibration from its own ground control)." if cal else "."))
    title = f"Camera {cam} pointing correction, {d1} to {d2}: waterlines against the lidar"
    subtitle = (f"Same {len(Z)} waterline points from {len(names)} photos, placed with the old pointing "
                f"({before_path.name}) and the {after_label} ({after_path.name}).\n"
                f"Old: azimuth {a0[0]:.2f}, tilt {a0[1]:.2f}, roll {a0[2]:.2f} deg.   "
                f"{after_label[0].upper() + after_label[1:]}: azimuth {a1[0]:.2f}, tilt {a1[1]:.2f}, roll {a1[2]:.2f} deg.\n"
                f"Change: azimuth {dA[0]:+.2f}, tilt {dA[1]:+.2f}, roll {dA[2]:+.2f} deg"
                + (f", camera moved {np.hypot(*(eo_after[:2] - eo_before[:2])):.1f} m"
                   if np.hypot(*(eo_after[:2] - eo_before[:2])) > 0.05 else "") + "; "
                f"the lines move median {np.nanmedian(moved):.1f} m (p90 {np.nanpercentile(moved, 90):.1f} m).")

    stem = args.output or f"pointing_fix_{cam}_{d1}_to_{d2}"
    figure(stem + ".png", survey, (eo_after[0], eo_after[1]), (E0, N0), (E1, N1), r0, r1, dist, F, best,
           s0, s1, b0, b1, sens, title, subtitle, label, after_label=after_label)

    with open(stem + ".csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pointing", "band_m", "n_on_survey", "median_m", "p25_m", "p75_m",
                    "nmad_m", "rmse_m", "p90_abs_m", "on_survey_pct"])
        for name, s, rows in (("old", s0, b0), ("fitted", s1, b1)):
            w.writerow([name, "all", s["n"], f"{s['median']:.3f}", "", "", f"{s['nmad']:.3f}",
                        f"{s['rmse']:.3f}", f"{s['p90_abs']:.3f}", f"{s['on_pct']:.1f}"])
            for (a, b), (n, med, lo, hi) in zip(BANDS, rows):
                w.writerow([name, f"{a}-{b}", n, f"{med:.3f}", f"{lo:.3f}", f"{hi:.3f}", "", "", "", ""])

    print(title)
    print(subtitle)
    print(label)
    for name, s in (("old pointing", s0), (after_label, s1)):
        print(f"  {name}: median {s['median']:+.3f} m, NMAD {s['nmad']:.3f} m, RMSE {s['rmse']:.3f} m, "
              f"{s['on_pct']:.0f}% of points on the survey ({s['n']})")
    print(f"wrote {stem}.png and {stem}.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())

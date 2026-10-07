#!/usr/bin/env python3
"""
Camera Re-Aim: One Survey Drawn On A Photo Before And After
==============================================================
Shows whether a camera was turned between two photos, with nothing but
the photos, a lidar survey and the candidate pointings -- no waterline
detection, no fitting to waterlines.

THE IDEA. Elevation contours of a lidar survey (e.g. 0 m and +3 m NAVD88)
are projected into a photo with a pointing (EO). With the right pointing
each line lies where the beach has that elevation: the 0 m line on the
lower beach, a little landward of the water's edge when the water is
below 0 m; +3 m on the upper beach. With a pointing ~20 deg off the lines
miss by tens of metres -- out on the sea, or up the dune -- which anyone
can see. Doing it on a photo from just before and one from just after a
suspected re-aim, with the same two pointings, shows the change directly:
if one pointing fits before and the other after, the camera was turned
in between.

FOR CACO05, JAN 2025. The 23 Jan 2025 lidar, photos of 22 Jan (CACO03)
and 25 Jan (CACO04); yellow = the 2025-02-19 calibration's aim (each
photo's horizon-fitted pointing), cyan = the pointing fitted to the
23 Jan lidar (fit_eo_to_survey.py, ~22 deg more pan):

    python3 reaim_figure.py --camera c2 \\
        --survey /mnt/I2Rgus_Data/Chelsea_calibration/2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif \\
        --survey-label "23 Jan 2025 lidar" \\
        --before <...Jan.22_17_00_00.GMT.2025.CACO03.c2.timex.jpg> \\
        --after  <...Jan.25_17_00_00.GMT.2025.CACO04.c2.timex.jpg> \\
        --calibration-eo CACO05_c2_2025-01-08_to_2025-01-23_EO.yaml CACO05_c2_2025-01-24_to_2025-02-15_EO.yaml \\
        --fitted-eo CACO05_c2_2025-01-18_to_2025-01-23_lidar_EO.yaml

The water level at each photo (ADCP, NAVD88; --water-level) is printed
on the panel so the reader can judge where the water's edge should be.
The survey is of one day: the further a photo is from it, the more the
beach itself may have changed (metres, not the tens of metres of a
~20 deg turn). Writes <output>.png.
"""

import re
import sys
import argparse
from pathlib import Path
from datetime import datetime, timezone

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from georectify import load_extrinsics, load_intrinsics          # noqa: E402
from compare_dem_survey import read_survey                       # noqa: E402
from project_survey import contour_segments, project_pieces      # noqa: E402

CAL = HERE / "calibration"
CHELSEA = Path("/mnt/I2Rgus_Data/Chelsea_calibration")
INK, MUTED = "#1f1f1e", "#6b6a64"
# Lines go on photos of sand and sea: the brightest hues, each with a black edge.
COLOURS = ["#ffd400", "#00e5ff", "#ff4fd8"]
STYLES = ["-", (0, (6, 3)), (0, (1.5, 2.5)), (0, (8, 3, 1.5, 3))]
EPOCH_RE = re.compile(r"^(\d{9,11})\.")
STATION_RE = re.compile(r"\.(CACO\d+)\.c\d\.", re.I)


def eo_path(p):
    return Path(p) if Path(p).exists() else CAL / p


def photo_time(path):
    m = EPOCH_RE.match(Path(path).name)
    return int(m.group(1)) if m else None


def water_level_source(path):
    """(epochs, levels, description) from an adcp_water_level_navd88.csv, or from the
    repo's ADCP file with the fixed MSL-to-NAVD88 offset; None if neither exists."""
    if path and Path(path).exists():
        from process_chelsea import load_water_level
        ep, lv = load_water_level(path)
        return ep, lv, f"ADCP on NAVD88 ({Path(path).name})"
    raw = HERE / "sig1000_waves_ALL.csv"
    if raw.exists():
        import csv
        ep, lv = [], []
        with open(raw, newline="") as f:
            for r in csv.DictReader(f):
                t = datetime.strptime(r["time"], "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
                ep.append(t.timestamp())
                lv.append(float(r["water_level"]) + 0.09)
        o = np.argsort(ep)
        return np.array(ep)[o], np.array(lv)[o], "ADCP + 0.09 m to NAVD88 (sig1000_waves_ALL.csv)"
    return None


def angles(eo):
    return np.degrees(eo[3:6])


def draw_panel(ax, img, pieces_by_pointing, levels, title, subtitle):
    import matplotlib.patheffects as pe
    ax.imshow(img)
    for k, pieces in enumerate(pieces_by_pointing):
        for lv, U, V in pieces:
            ax.plot(U, V, color=COLOURS[k], lw=2.0, ls=STYLES[levels.index(lv) % len(STYLES)],
                    solid_capstyle="round",
                    path_effects=[pe.Stroke(linewidth=4.2, foreground="black"), pe.Normal()])
    ax.set_xlim(0, img.shape[1])
    ax.set_ylim(img.shape[0], 0)
    ax.axis("off")
    ax.set_title(title, loc="left", fontsize=12, color=INK, weight="bold", pad=22)
    ax.text(0, 1.012, subtitle, transform=ax.transAxes, fontsize=9.5, color=MUTED,
            ha="left", va="bottom")


def main():
    import cv2
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from matplotlib.lines import Line2D

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--camera", required=True, choices=["c1", "c2"])
    ap.add_argument("--survey", required=True, help="lidar DSM (GeoTIFF or .asc), NAVD88")
    ap.add_argument("--survey-label", default=None, help='e.g. "23 Jan 2025 lidar"')
    ap.add_argument("--before", required=True, help="original photo before the suspected re-aim")
    ap.add_argument("--after", required=True, help="original photo after it")
    ap.add_argument("--calibration-eo", nargs="+", required=True, metavar="EO",
                    help="the calibration's aim: one EO for both photos, or one for the before "
                         "photo and one for the after photo (each period's horizon-fitted EO)")
    ap.add_argument("--fitted-eo", nargs="+", required=True, metavar="EO",
                    help="the other candidate (e.g. the lidar-fitted EO): one, or before and after")
    ap.add_argument("--labels", nargs=2, default=None, metavar=("CALIBRATION", "FITTED"),
                    help="legend names of the two pointings")
    ap.add_argument("--levels", nargs="+", type=float, default=[0.0, 3.0],
                    help="contour elevations, m NAVD88 (default 0 3); line style per level")
    ap.add_argument("--range", type=float, default=250.0,
                    help="only survey within this distance of the camera, m (default 250)")
    ap.add_argument("--water-level", default=str(CHELSEA / "adcp_water_level_navd88.csv"),
                    help="adcp_water_level_navd88.csv (adcp_to_navd88.py); falls back to the "
                         "repo's ADCP file + 0.09 m")
    ap.add_argument("--exposure-mid", type=float, default=300.0,
                    help="seconds from the file-name time to the middle of the exposure (default 300)")
    ap.add_argument("--title", default=None, help="page title (default: what is drawn)")
    ap.add_argument("--note", default=None, help="one line under the title, e.g. the conclusion")
    ap.add_argument("--output", default=None, help="output stem (default reaim_<cam>_<before>_<after>)")
    args = ap.parse_args()

    def pair(files, what):
        if len(files) not in (1, 2):
            ap.error(f"{what}: give one EO, or two (before, after)")
        eos = [load_extrinsics(eo_path(p)) for p in files]
        return (eos[0], eos[-1]), (Path(files[0]).name, Path(files[-1]).name)

    cal_eos, cal_names = pair(args.calibration_eo, "--calibration-eo")
    fit_eos, fit_names = pair(args.fitted_eo, "--fitted-eo")
    io = load_intrinsics(CAL / f"CACO05_{args.camera}_20240801_IO.yaml")
    levels = sorted(set(args.levels))

    imgs, times, stations = [], [], []
    for p in (args.before, args.after):
        im = cv2.imread(str(p))
        if im is None:
            sys.exit(f"cannot read {p}")
        if (im.shape[1], im.shape[0]) != (int(io[0]), int(io[1])):
            sys.exit(f"{p}: {im.shape[1]} x {im.shape[0]} px, but the lens model is "
                     f"{int(io[0])} x {int(io[1])} -- use the original photo, not a redrawn one")
        imgs.append(cv2.cvtColor(im, cv2.COLOR_BGR2RGB))
        times.append(photo_time(p))
        m = STATION_RE.search(Path(p).name)
        stations.append(m.group(1).upper() if m else None)

    wl_src = water_level_source(args.water_level)
    levels_at = []
    for t in times:
        wl = None
        if wl_src is not None and t is not None:
            from process_chelsea import level_at
            wl = level_at(wl_src[0], wl_src[1], t + args.exposure_mid, 3 * 3600)
        levels_at.append(wl)

    grid, x0, y0, cell = read_survey(args.survey)
    cx, cy = cal_eos[0][0], cal_eos[0][1]
    segs = contour_segments(grid, x0, y0, cell, levels,
                            (cx - args.range, cx + args.range, cy - args.range, cy + args.range))
    print("contours: " + ", ".join(f"{lv:+g} m ({len(s)} pieces)" for lv, s in segs.items()))
    if not any(segs.values()):
        sys.exit("no contours within --range of the camera: wrong survey or levels?")

    survey_label = args.survey_label or Path(args.survey).stem
    cal_az, fit_az = angles(cal_eos[1])[0], angles(fit_eos[0])[0]
    labels = args.labels or [f"2025-02-19 calibration's aim (azimuth {cal_az:.1f}°)",
                             f"fitted to the lidar (azimuth {fit_az:.1f}°)"]

    fig = plt.figure(figsize=(16, 9.4), dpi=150, facecolor="white")
    gs = fig.add_gridspec(1, 2, left=0.03, right=0.97, top=0.83, bottom=0.17, wspace=0.05)
    for i, (when, img, t, st, wl) in enumerate(zip(("Before", "After"), imgs, times, stations,
                                                   levels_at)):
        pieces = [project_pieces(segs, io, cal_eos[i]), project_pieces(segs, io, fit_eos[i])]
        for name, pc in zip((cal_names[i], fit_names[i]), pieces):
            print(f"  {when.lower()}: {name}: {len(pc)} contour piece(s) in the frame")
        stamp = (datetime.fromtimestamp(t, tz=timezone.utc).strftime("%d %b %Y, %H:%M UTC")
                 if t is not None else Path((args.before, args.after)[i]).name)
        title = f"{when}: {stamp}" + (f"  (station ID {st})" if st else "")
        sub = (f"water level {wl:+.2f} m NAVD88" if wl is not None else "water level not available")
        draw_panel(fig.add_subplot(gs[0, i]), img, pieces, levels, title, sub)
        for name, pc in zip(labels, pieces):
            if not pc:
                print(f"  NOTE: {when.lower()}: no {name} contour falls in the frame")

    title = args.title or (f"Camera {args.camera}: the {survey_label} drawn on a photo from "
                           f"before and after, with two candidate pointings")
    fig.text(0.03, 0.975, title, ha="left", va="top", fontsize=15, color=INK, weight="bold")
    y = 0.925
    if args.note:
        fig.text(0.03, y, args.note, ha="left", va="top", fontsize=11, color=INK)
        y -= 0.035
    fig.text(0.03, y,
             f"Lines: elevation contours of the {survey_label}, projected into each photo with each "
             f"pointing. With the right pointing a line lies where the beach has that elevation "
             f"(the water's edge is at the water level given for each photo). A turn moves the "
             f"lines nearest the camera most; far lines on a straight beach move little.",
             ha="left", va="top", fontsize=9.5, color=MUTED, wrap=True)

    stroke = [pe.Stroke(linewidth=4.2, foreground="black"), pe.Normal()]
    handles = [Line2D([], [], color=c, lw=2.4, path_effects=stroke) for c in COLOURS[:2]]
    handles += [Line2D([], [], color="#9a9a94", lw=2.0, ls=STYLES[k % len(STYLES)])
                for k in range(len(levels))]
    names = labels + [f"{lv:+g} m NAVD88" if lv else "0 m NAVD88" for lv in levels]
    fig.legend(handles, names, loc="lower left", bbox_to_anchor=(0.03, 0.075), ncol=2 + len(levels),
               frameon=False, fontsize=10, handlelength=3.2, labelcolor=INK)
    files = (f"Pointings: calibration {cal_names[0]}" + (f" / {cal_names[1]}" if cal_names[1] != cal_names[0] else "")
             + f"; fitted {fit_names[0]}" + (f" / {fit_names[1]}" if fit_names[1] != fit_names[0] else "")
             + f".  Photos: {Path(args.before).name}, {Path(args.after).name}."
             + (f"  Water level: {wl_src[2]}, at mid-exposure." if wl_src else ""))
    fig.text(0.03, 0.03, files, ha="left", va="bottom", fontsize=7.5, color=MUTED, wrap=True)

    stem = args.output or "reaim_{}_{}_{}".format(
        args.camera,
        *[datetime.fromtimestamp(t, tz=timezone.utc).strftime("%Y-%m-%d") if t else x
          for t, x in zip(times, ("before", "after"))])
    out = Path(stem).with_suffix(".png") if not str(stem).endswith(".png") else Path(stem)
    fig.savefig(out, facecolor="white")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Camera Re-Aim: One Survey Drawn On A Photo Before And After
==============================================================
Shows whether a camera was turned between two photos, with nothing but
the photos, a lidar survey and the candidate pointings -- no waterline
detection, no fitting to waterlines.

THE IDEA. Elevation contours of a lidar survey are projected into a photo
with a pointing (EO). With the right pointing each line lies where the
beach has that elevation. The solid line is the contour at the water
level measured when the photo was taken, so it should run along the
water's edge -- on it, or a little seaward, because in a time-exposure
wave setup and run-up lift the visible edge above still water by
decimetres. A dashed +3 m line falls on the upper beach. With a pointing
~20 deg off the lines miss by tens of metres -- out on the sea, or up the
dune. Doing it on a photo from just before and one from just after a
suspected re-aim, with the same two pointings, shows the change directly:
if one pointing fits before and the other after, the camera was turned in
between. A turn moves the lines nearest the camera most; far lines on a
straight beach move little, so judge by the near ones.

HONEST LABEL. A pointing fitted to this survey (fit_eo_to_survey.py; its
EO file says so) fits photos of the days it was fitted on by
construction. The page says so; the independent test is the other photo.

FOR CACO05, JAN 2025. The 23 Jan 2025 lidar, photos of 22 Jan (CACO03)
and 25 Jan (CACO04); yellow = the 2025-02-19 calibration's aim (each
photo's horizon-fitted pointing), cyan = the station's own calibration of
the old setup, made on the day of the lidar (CACO03_c2_20250123_EO.yaml:
camera ~5.8 m away, turned ~15 deg). It was not fitted to the lidar, so
both panels are independent tests. (The lidar-fitted pointing,
CACO05_c2_2025-01-18_to_2025-01-23_lidar_EO.yaml, ~22 deg of pan with the
position kept, agrees with it to ~1 m cross-shore.)

    python3 reaim_figure.py --camera c2 \\
        --survey /mnt/I2Rgus_Data/Chelsea_calibration/2025005FA_Marconi_Jan_YSMP_Lidar_DSM_25cm.tif \\
        --survey-label "23 Jan 2025 lidar" \\
        --before <...Jan.22_17_00_00.GMT.2025.CACO03.c2.timex.jpg> \\
        --after  <...Jan.25_17_00_00.GMT.2025.CACO04.c2.timex.jpg> \\
        --calibration-eo CACO05_c2_2025-01-08_to_2025-01-23_EO.yaml CACO05_c2_2025-01-24_to_2025-02-15_EO.yaml \\
        --fitted-eo CACO03_c2_20250123_EO.yaml

--levels takes elevations in m NAVD88 and the word "water" (the photo's
measured water level; default: water 3). The water level comes from the
ADCP on NAVD88 (--water-level), at the middle of the exposure. The survey
is of one day: the further a photo is from it, the more the beach itself
may have changed (metres, not the tens of metres of a ~20 deg turn).
Writes <output>.png.
"""

import re
import sys
import textwrap
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
COLOURS = ["#ffd400", "#00e5ff"]
STYLES = ["-", (0, (6, 3)), (0, (1.5, 2.5)), (0, (8, 3, 1.5, 3))]
EPOCH_RE = re.compile(r"^(\d{9,11})\.")
STATION_RE = re.compile(r"\.(CACO\d+)\.c\d\.", re.I)
STATION_CAL_RE = re.compile(r"^(CACO\d+)_c\d_(\d{8})_EO(?:-CV)?\.yaml$")
FITTED_RE = re.compile(r"fit_eo_to_survey\.py:.*fitted to (\S+) using \d+ waterline frames "
                       r"(\d{4}-\d{2}-\d{2}) to (\d{4}-\d{2}-\d{2})")


def eo_path(p):
    return Path(p) if Path(p).exists() else CAL / p


def eo_note(p):
    """The comment lines of an EO file: the tools that write EOs say there how they made it."""
    try:
        return " ".join(l[1:].strip() for l in Path(p).read_text().splitlines() if l.startswith("#"))
    except OSError:
        return ""


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


def describe(files, cam, survey_name):
    """Default legend name of a pointing, from its file names and notes, and the survey
    fit it came from (survey, first day, last day) if it was fitted to this survey."""
    notes = [eo_note(eo_path(f)) for f in files]
    cal = f"CACO05_{cam}_20250219_EO.yaml"
    fits = [FITTED_RE.search(n) for n in notes]
    if all(m for m in fits):
        survey = fits[0].group(1)
        same = all(m.group(1) == survey_name for m in fits)
        name = "fitted to this survey" if same else f"fitted to {survey}"
        return name, ((fits[0].group(2), fits[0].group(3)) if same else None)
    if all(Path(f).name == cal or n.startswith(f"horizon_check.py: {cal}")
           for f, n in zip(files, notes)):
        return "2025-02-19 calibration's aim", None
    setups = [STATION_CAL_RE.match(Path(f).name) for f in files]
    if all(setups) and len({m.groups() for m in setups}) == 1:
        st, d = setups[0].groups()
        return f"{st} calibration of {d[:4]}-{d[4:6]}-{d[6:]}", None
    return " / ".join(dict.fromkeys(Path(f).name for f in files)), None


def azimuth_text(pair):
    a, b = (np.degrees(e[3]) for e in pair)
    return f"azimuth {a:.1f}°" if abs(a - b) < 0.05 else f"azimuth {a:.1f}° before, {b:.1f}° after"


def parse_levels(tokens, ap):
    out = []
    for t in tokens:
        if t.lower() == "water":
            out.append("water")
        else:
            try:
                out.append(float(t))
            except ValueError:
                ap.error(f"--levels: {t!r} is neither a number nor 'water'")
    out = list(dict.fromkeys(out))
    if len(out) > len(STYLES):
        ap.error(f"--levels: at most {len(STYLES)} (one line style each)")
    return out


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
                    help="legend names of the two pointings (default: from the EO files)")
    ap.add_argument("--levels", nargs="+", default=["water", "3"],
                    help="contours: m NAVD88, or 'water' for the photo's measured water level "
                         "(default: water 3); one line style each")
    ap.add_argument("--range", type=float, default=250.0,
                    help="only survey within this distance of the camera, m (default 250)")
    ap.add_argument("--water-level", default=str(CHELSEA / "adcp_water_level_navd88.csv"),
                    help="adcp_water_level_navd88.csv (adcp_to_navd88.py); falls back to the "
                         "repo's ADCP file + 0.09 m")
    ap.add_argument("--exposure-mid", type=float, default=300.0,
                    help="seconds from the file-name time to the middle of the exposure (default 300)")
    ap.add_argument("--title", default=None, help="page title (default: what is drawn)")
    ap.add_argument("--note", default=None, help="one line under the title, e.g. the conclusion")
    ap.add_argument("--output", default=None,
                    help="output file (.png added if missing; default reaim_<cam>_<before>_<after>.png)")
    args = ap.parse_args()

    def pair(files, what):
        if len(files) not in (1, 2):
            ap.error(f"{what}: give one EO, or two (before, after)")
        eos = [load_extrinsics(eo_path(p)) for p in files]
        return (eos[0], eos[-1]), (Path(files[0]).name, Path(files[-1]).name)

    tokens = parse_levels(args.levels, ap)
    pointings = [pair(args.calibration_eo, "--calibration-eo"), pair(args.fitted_eo, "--fitted-eo")]
    io = load_intrinsics(CAL / f"CACO05_{args.camera}_20240801_IO.yaml")

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
    water = []
    for t in times:
        wl = None
        if wl_src is not None and t is not None:
            from process_chelsea import level_at
            wl = level_at(wl_src[0], wl_src[1], t + args.exposure_mid, 3 * 3600)
        water.append(wl)
    if "water" in tokens and None in water:
        print("NOTE: no water level for a photo; its 'water' line is drawn at 0 m instead")

    survey_name = Path(args.survey).name
    survey_label = args.survey_label or Path(args.survey).stem
    grid, x0, y0, cell = read_survey(args.survey)
    cx, cy = pointings[0][0][0][0], pointings[0][0][0][1]
    box = (cx - args.range, cx + args.range, cy - args.range, cy + args.range)

    described = [describe(files, args.camera, survey_name)
                 for files in (args.calibration_eo, args.fitted_eo)]
    labels = [f"{lab} ({azimuth_text(eos)})" for lab, eos in
              zip(args.labels or [d[0] for d in described], (p[0] for p in pointings))]

    fig = plt.figure(figsize=(16, 9.6), dpi=150, facecolor="white")
    stroke = [pe.Stroke(linewidth=4.2, foreground="black"), pe.Normal()]

    # header: title, optional note, what the lines are, and the honest label
    lines = [(args.title or f"Camera {args.camera}: the {survey_label} drawn on a photo from "
              f"before and after, with two candidate pointings", 15, INK, "bold")]
    if args.note:
        lines.append((args.note, 11, INK, "normal"))
    lines.append((f"Lines: contours of the {survey_label}, projected into each photo with each "
                  f"pointing; with the right pointing a line lies where the beach has that elevation. "
                  f"The solid line is the contour at the water level measured when the photo was "
                  f"taken: it should run along the water's edge, on it or a little seaward, since in "
                  f"a 10-minute exposure wave run-up lifts the visible edge above still water. A turn "
                  f"moves the lines nearest the camera most.", 9.5, MUTED, "normal"))
    for k, (_, fitted) in enumerate(described):
        if fitted:
            in_fit = [t is not None and fitted[0] <= datetime.fromtimestamp(t, tz=timezone.utc)
                      .strftime("%Y-%m-%d") <= fitted[1] for t in times]
            which = " and ".join(w for w, f in zip(("before", "after"), in_fit) if f)
            lines.append((f"The {'yellow' if k == 0 else 'cyan'} pointing was fitted to this survey "
                          f"from waterlines of {fitted[0]} to {fitted[1]}, so its fit on photos of those "
                          f"days is expected" + (f" (the {which} photo)" if which else "")
                          + "; the independent test is a photo outside them.", 9.5, INK, "normal"))
    y = 0.978
    for text, size, colour, weight in lines:
        wrapped = textwrap.fill(text, width=int(2150 / size))
        fig.text(0.03, y, wrapped, ha="left", va="top", fontsize=size, color=colour, weight=weight)
        y -= (wrapped.count("\n") + 1) * size * 1.45 / 72 / 9.6 + 0.008
    top = y - 0.065

    gs = fig.add_gridspec(1, 2, left=0.03, right=0.97, top=top, bottom=0.17, wspace=0.05)
    for i, (when, img, t, st, wl) in enumerate(zip(("Before", "After"), imgs, times, stations, water)):
        values = [(wl if wl is not None else 0.0) if tk == "water" else tk for tk in tokens]
        keep = [k for k, v in enumerate(values) if v not in values[:k]]
        segs = contour_segments(grid, x0, y0, cell, [values[k] for k in keep], box)
        print(f"{when.lower()}: contours " + ", ".join(f"{lv:+.2f} m ({len(s)} pieces)"
                                                       for lv, s in segs.items()))
        if not any(segs.values()):
            sys.exit("no contours within --range of the camera: wrong survey, camera or levels?")
        ax = fig.add_subplot(gs[0, i])
        ax.imshow(img)
        for k, (eos, names) in enumerate(pointings):
            pieces = project_pieces(segs, io, eos[i])
            print(f"  {names[i]}: {len(pieces)} contour piece(s) in the frame")
            for lv, U, V in pieces:
                style = STYLES[keep[int(np.argmin([abs(values[j] - lv) for j in keep]))]]
                ax.plot(U, V, color=COLOURS[k], lw=2.0, ls=style, solid_capstyle="round",
                        path_effects=stroke)
        ax.set_xlim(0, img.shape[1])
        ax.set_ylim(img.shape[0], 0)
        ax.axis("off")
        stamp = (datetime.fromtimestamp(t, tz=timezone.utc).strftime("%d %b %Y, %H:%M UTC")
                 if t is not None else Path((args.before, args.after)[i]).name)
        ax.set_title(f"{when}: {stamp}" + (f"  (station ID {st})" if st else ""),
                     loc="left", fontsize=12, color=INK, weight="bold", pad=22)
        sub = (f"water level {wl:+.2f} m NAVD88 (still water, ADCP)" if wl is not None
               else "water level not available")
        ax.text(0, 1.012, sub, transform=ax.transAxes, fontsize=9.5, color=MUTED,
                ha="left", va="bottom")

    # legend: the pointings on one row, the line styles on the next
    leg1 = fig.legend([Line2D([], [], color=c, lw=2.4, path_effects=stroke) for c in COLOURS], labels,
                      loc="lower left", bbox_to_anchor=(0.03, 0.105), ncol=2, frameon=False,
                      fontsize=10, handlelength=3.2)
    style_names = ["contour at the photo's water level" if tk == "water"
                   else ("0 m NAVD88" if tk == 0 else f"{tk:+g} m NAVD88") for tk in tokens]
    fig.legend([Line2D([], [], color="#8a8a84", lw=2.0, ls=STYLES[k]) for k in range(len(tokens))],
               style_names, loc="lower left", bbox_to_anchor=(0.03, 0.07), ncol=len(tokens),
               frameon=False, fontsize=10, handlelength=3.2)
    fig.add_artist(leg1)
    files = ("Pointings: yellow " + " / ".join(dict.fromkeys(pointings[0][1]))
             + "; cyan " + " / ".join(dict.fromkeys(pointings[1][1]))
             + f".  Photos: {Path(args.before).name}, {Path(args.after).name}."
             + (f"  Water level: {wl_src[2]}, at mid-exposure." if wl_src else ""))
    fig.text(0.03, 0.03, textwrap.fill(files, width=260), ha="left", va="bottom", fontsize=7.5,
             color=MUTED)

    stem = args.output or "reaim_{}_{}_{}".format(
        args.camera,
        *[datetime.fromtimestamp(t, tz=timezone.utc).strftime("%Y-%m-%d") if t else x
          for t, x in zip(times, ("before", "after"))])
    out = Path(stem) if str(stem).lower().endswith(".png") else Path(f"{stem}.png")
    fig.savefig(out, facecolor="white")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

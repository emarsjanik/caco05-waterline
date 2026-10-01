#!/usr/bin/env python3
"""
Waterlines And DEM Against RTK Transects
===========================================
Checks the camera's beach elevations against an RTK GNSS survey of
cross-shore transects (Emlid Flow CSV export: Name, Easting, Northing,
Elevation, Description; NAD83(2011) / UTM 19N + NAVD88).

WATERLINES, NOT ONLY THE DEM. The survey walks the dry beach; the camera
sees only the band the tide reached in daylight. Their overlap is the
lower end of each transect. The DEM can only be checked where the two
grids overlap, but every waterline whose PHYSICAL position falls on a
surveyed stretch of transect can be checked directly: the RTK profile
gives the true sand elevation at that spot, the waterline was given the
water level. If the waterlines sit too high on the beach (wave setup and
swash) or the water level is on the wrong datum, every such pair shows
it -- even when the DEM and the survey hardly overlap.

For each transect (consecutive 'Transect' points; a jump of more than
--split m starts a new one) the profile is the RTK points joined in
order of distance along the line. A waterline point counts if it lies
within --tolerance m of the line and between its first and last surveyed
point (no extrapolation). Per frame: the median of RTK - waterline
elevation; then the same waves-versus-datum split as
compare_dem_survey.py (offset = a + b*sqrt(Hs L0)).

--dem samples a DEM (.asc) at every RTK point too.
--plot draws each transect: the RTK profile, the waterline points at
their positions along it (coloured by offshore Hs), and the DEM.

Usage:
    python3 compare_rtk.py 2026-09-29_Marconi_Checkshots.csv \\
        --contours /mnt/I2Rgus_Data/waterline/contour_points_ground.csv \\
        --start-date 2026-09-29 --end-date 2026-09-29 \\
        --dem dem_2026-09-29_fine_dem.asc --plot rtk_check.png
"""

import sys
import csv
import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from compare_dem_survey import waves_vs_datum, sample      # noqa: E402
from asc_to_geotiff import read_asc                        # noqa: E402


def load_transects(path, split):
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            if r.get("Description", "").strip().lower() != "transect":
                continue
            rows.append((float(r["Easting"]), float(r["Northing"]), float(r["Elevation"]),
                         r["Name"]))
    lines, cur = [], []
    for p in rows:
        if cur and np.hypot(p[0] - cur[-1][0], p[1] - cur[-1][1]) > split:
            lines.append(cur); cur = []
        cur.append(p)
    if cur:
        lines.append(cur)
    out = []
    for pts in lines:
        if len(pts) < 3:
            continue
        P = np.array([p[:3] for p in pts])
        c = P[:, :2].mean(0)
        _, vecs = np.linalg.eigh(np.cov((P[:, :2] - c).T))
        d = vecs[:, 1]
        # point the axis seaward: towards the lowest end
        s = (P[:, :2] - c) @ d
        if np.corrcoef(s, P[:, 2])[0, 1] > 0:
            d = -d; s = -s
        o = np.argsort(s)
        out.append({"name": f"T{len(out) + 1} ({pts[0][3]}-{pts[-1][3]})", "c": c, "d": d,
                    "s": s[o], "z": P[o, 2], "E": P[o, 0], "N": P[o, 1]})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("rtk", help="RTK points CSV (Emlid Flow export)")
    ap.add_argument("--contours", default=str(HERE / "contour_points_ground.csv"),
                    help="Georectified waterline points (georectify.py output)")
    ap.add_argument("--start-date"); ap.add_argument("--end-date")
    ap.add_argument("--camera", default="both")
    ap.add_argument("--tolerance", type=float, default=2.0,
                    help="Max distance of a waterline point from a transect line, m (default 2)")
    ap.add_argument("--split", type=float, default=15.0,
                    help="Gap between consecutive RTK points that starts a new transect, m")
    ap.add_argument("--min-points", type=int, default=3,
                    help="Waterline points a frame needs on transects to count (default 3)")
    ap.add_argument("--dem", help="DEM (.asc) to sample at the RTK points")
    ap.add_argument("--output", default=None, help="Per-frame CSV (default <rtk stem>_frames.csv)")
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()

    T = load_transects(args.rtk, args.split)
    if not T:
        sys.exit("no transects (Description 'Transect') in the RTK file")
    print(f"RTK transects     : {len(T)}")
    for t in T:
        print(f"   {t['name']:16s} {len(t['s']):2d} pts, {t['s'][-1] - t['s'][0]:5.1f} m long, "
              f"z {t['z'].min():+.2f} to {t['z'].max():+.2f}")

    frames = {}
    with open(args.contours, newline="") as f:
        for r in csv.DictReader(f):
            if not r.get("easting_utm19"):
                continue
            day = r.get("capture_time_utc", "")[:10]
            if (args.start_date and day < args.start_date) or (args.end_date and day > args.end_date):
                continue
            if args.camera != "both" and r.get("camera") != args.camera:
                continue
            fr = frames.setdefault(r["source_file"], {
                "time": r.get("capture_time_utc", ""), "cam": r.get("camera", ""),
                "E": [], "N": [], "Z": [],
                "tide": float(r["tide_elevation_navd88"]),
                "hs": float(r["offshore_hs_m"]) if r.get("offshore_hs_m") else np.nan,
                "tp": float(r["offshore_tp_s"]) if r.get("offshore_tp_s") else np.nan})
            fr["E"].append(float(r["easting_utm19"]))
            fr["N"].append(float(r["northing_utm19"]))
            fr["Z"].append(float(r.get("beach_elevation_navd88") or r["tide_elevation_navd88"]))
    print(f"waterline frames  : {len(frames)} in {args.start_date or 'all'} to {args.end_date or 'all'}")
    if not frames:
        sys.exit("no georectified waterlines in that date range -- has the station processed "
                 "those days yet? (GNSS-R water levels run ~2 days behind)")

    rows, hits = [], []           # hits: (transect index, s, z_waterline, hs) for the plot
    for name, fr in sorted(frames.items(), key=lambda kv: kv[1]["time"]):
        E, N, Z = np.array(fr["E"]), np.array(fr["N"]), np.array(fr["Z"])
        diffs = []
        for k, t in enumerate(T):
            rel = np.stack([E - t["c"][0], N - t["c"][1]], 1)
            s = rel @ t["d"]
            lat = np.abs(rel @ np.array([-t["d"][1], t["d"][0]]))
            m = (lat <= args.tolerance) & (s >= t["s"][0]) & (s <= t["s"][-1])
            if not m.any():
                continue
            zr = np.interp(s[m], t["s"], t["z"])
            diffs.append(zr - Z[m])
            hits += [(k, ss, zz, fr["hs"]) for ss, zz in zip(s[m], Z[m])]
        d = np.concatenate(diffs) if diffs else np.array([])
        if len(d) < args.min_points:
            continue
        rows.append({"frame": name, "time_utc": fr["time"], "camera": fr["cam"],
                     "tide_navd88": fr["tide"], "hs_m": fr["hs"], "tp_s": fr["tp"],
                     "n_points": len(d), "offset_m": float(np.median(d)),
                     "spread_m": float(np.percentile(d, 84) - np.percentile(d, 16))})

    print(f"frames on transects: {len(rows)} (waterline within {args.tolerance:g} m of a "
          f"transect, inside its surveyed stretch)")
    if rows:
        off = np.array([r["offset_m"] for r in rows])
        print(f"OFFSET RTK - waterline elevation: median {np.median(off):+.3f} m over frames, "
              f"p10-p90 {np.percentile(off, 10):+.2f} to {np.percentile(off, 90):+.2f} m")
        print("   positive = waterline elevations are LOW (lines sit higher on the beach "
              "than their water level)")
        for cam in sorted({r["camera"] for r in rows}):
            o = [r["offset_m"] for r in rows if r["camera"] == cam]
            print(f"   {cam}: {len(o)} frames, median {np.median(o):+.3f} m")
        print()
        a, b, phi, w = waves_vs_datum(off, np.array([r["hs_m"] for r in rows]),
                                      np.array([r["tp_s"] for r in rows]),
                                      np.array([r["tide_navd88"] for r in rows]))
        out = args.output or str(Path(args.rtk).with_suffix("")) + "_frames.csv"
        with open(out, "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=list(rows[0]))
            wr.writeheader()
            wr.writerows(rows)
        print(f"\nwrote {out}")
    else:
        print("   none -- the waterlines do not reach the surveyed parts of the transects. "
              "The --plot shows where they fall.")

    dem_prof = {}
    if args.dem:
        g, h = read_asc(args.dem)
        x0, y0, c = h["xllcorner"], h["yllcorner"] + g.shape[0] * h["cellsize"], h["cellsize"]
        allE = np.concatenate([t["E"] for t in T]); allN = np.concatenate([t["N"] for t in T])
        allZ = np.concatenate([t["z"] for t in T])
        zd = sample(g, x0, y0, c, allE, allN)
        ok = np.isfinite(zd)
        print(f"\nDEM at RTK points : {int(ok.sum())} of {len(allZ)} points on the DEM")
        if ok.any():
            dd = allZ[ok] - zd[ok]
            print(f"   RTK - DEM median {np.median(dd):+.3f} m, mean {dd.mean():+.3f}, "
                  f"RMS {np.sqrt((dd ** 2).mean()):.3f} m (n={ok.sum()})")
        for k, t in enumerate(T):
            ss = np.linspace(t["s"][0] - 30, t["s"][-1] + 30, 400)
            Es, Ns = t["c"][0] + ss * t["d"][0], t["c"][1] + ss * t["d"][1]
            dem_prof[k] = (ss, sample(g, x0, y0, c, Es, Ns))

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        n = len(T)
        ncol = 3
        nrow = int(np.ceil(n / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(15, 3.8 * nrow), dpi=100, squeeze=False,
                                 constrained_layout=True)
        hits = np.array(hits) if hits else np.zeros((0, 4))
        hs_all = hits[:, 3] if len(hits) else np.array([0.0])
        vmin, vmax = np.nanmin(hs_all) if len(hits) else 0, np.nanmax(hs_all) if len(hits) else 1
        sc = None
        for k, t in enumerate(T):
            ax = axes[k // ncol][k % ncol]
            ax.plot(t["s"], t["z"], "-o", color="#1f1f1e", ms=3, lw=1.5, label="RTK")
            if k in dem_prof:
                ax.plot(dem_prof[k][0], dem_prof[k][1], color="#eb6834", lw=2, label="camera DEM")
            m = hits[:, 0] == k if len(hits) else np.array([], bool)
            if m.any():
                sc = ax.scatter(hits[m, 1], hits[m, 2], c=hits[m, 3], cmap="viridis", s=8,
                                vmin=vmin, vmax=vmax, label="waterlines (z = water level)")
            ax.set_title(t["name"], loc="left", fontsize=9)
            ax.set_xlabel("distance seaward (m)"); ax.set_ylabel("m NAVD88")
            ax.grid(True, color="#e4e3dc"); ax.set_axisbelow(True)
            if k == 0:
                ax.legend(fontsize=7, frameon=False, loc="upper right")
        for k in range(n, nrow * ncol):
            axes[k // ncol][k % ncol].axis("off")
        if sc is not None:
            fig.colorbar(sc, ax=axes.ravel().tolist(), label="offshore Hs (m)", shrink=0.6)
        fig.savefig(args.plot, bbox_inches="tight")
        print(f"wrote {args.plot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

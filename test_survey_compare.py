#!/usr/bin/env python3
"""
Tests For survey_compare.py
==============================
Plain python3 (no pytest needed on the station):

    python3 test_survey_compare.py [--keep DIR]

WHY. The comparison is what the survey products are judged by, so its
arithmetic has to be right by construction, not by inspection: a
synthetic beach (a plane, so block medians and bilinear values are
exact) with the DEM put exactly +0.20 m above it must come back as
+0.20 m, from a GeoTIFF lidar on a grid offset from the DEM's, from an
Emlid RTK file, from CIRN target files and from a plain x,y,z CSV.
Water-masked lidar cells must leave DEM cells uncompared, points off
the DEM must be listed, and the label checks must catch a pointing
fitted to the survey, an envelope placed with it, a GCP calibration
dated the survey day and a long time gap.

HOW. Every test writes into a temporary folder (--keep DIR to look at
the outputs, figures included) and asserts on the returned headline and
on the files.
"""

import os
import sys
import csv
import json
import shutil
import argparse
import tempfile
import contextlib
import io
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import survey_compare as sc                                   # noqa: E402
from asc_to_geotiff import write_geotiff, read_asc            # noqa: E402
from compare_dem_survey import read_survey, read_tiff         # noqa: E402

OFFSET = 0.20
# The synthetic beach: rises 0.05 m per metre landward (west), 3 m at x = X0.
X0, Y0 = 420060.0, 4638340.0


def beach(E, N):
    return 3.0 - 0.05 * (E - X0) + 0.0 * N


def write_asc(path, grid_ns, xll, yll, cell):
    """grid rows north to south."""
    with open(path, "w") as f:
        f.write(f"ncols {grid_ns.shape[1]}\nnrows {grid_ns.shape[0]}\nxllcorner {xll}\n"
                f"yllcorner {yll}\ncellsize {cell}\nNODATA_value -9999.0\n")
        for r in grid_ns:
            f.write(" ".join("%.4f" % (v if np.isfinite(v) else -9999.0) for v in r) + "\n")


def make_dem(d, offset=OFFSET, ncols=40, nrows=60, cell=2.0):
    """DEM: the beach + offset at cell centres, x 420060-420140, y 4638340-4638460."""
    xs = X0 + (np.arange(ncols) + 0.5) * cell
    ys = Y0 + nrows * cell - (np.arange(nrows) + 0.5) * cell
    E, N = np.meshgrid(xs, ys)
    dem = beach(E, N) + offset
    dem[:3, :] = np.nan                       # a gap, as real DEMs have
    p = Path(d) / "synth_dem.asc"
    write_asc(p, dem, X0, Y0, cell)
    cnt = np.where(np.isfinite(dem), 12.0, np.nan)
    write_asc(Path(d) / "synth_count.asc", cnt, X0, Y0, cell)
    spr = np.where(np.isfinite(dem), 0.15, np.nan)
    write_asc(Path(d) / "synth_spread.asc", spr, X0, Y0, cell)
    return p, dem


def make_lidar(d, water_below=1.0, scell=0.25):
    """A 0.25 m lidar of the same beach on a grid offset 0.676 m from the DEM's (as the real
    2025 YSMP lidar is), water-masked below `water_below` m; a GeoTIFF."""
    x0, ytop = X0 - 30.324, Y0 + 160.676
    nc, nr = int(180 / scell), int(200 / scell)
    xs = x0 + (np.arange(nc) + 0.5) * scell
    ys = ytop - (np.arange(nr) + 0.5) * scell
    E, N = np.meshgrid(xs, ys)
    z = beach(E, N)
    z[z < water_below] = np.nan
    p = Path(d) / "synth_lidar_DSM.tif"
    write_geotiff(p, z, x0, ytop - nr * scell, scell, 6348, "synthetic lidar")
    return p


def quiet(fn, *a, **k):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = fn(*a, **k)
    return out, buf.getvalue()


def check(cond, msg):
    if not cond:
        raise AssertionError(msg)
    print(f"  ok  {msg}")


# ---------------------------------------------------------------------

def test_block_median(d):
    print("block_median")
    s = np.arange(64, dtype=float).reshape(8, 8)          # one 2 m cell of 0.25 m cells
    med, cov = sc.block_median(s, 0.0, 2.0, 0.25, 0.0, 2.0, 2.0, 1, 1)
    check(abs(med[0, 0] - np.median(s)) < 1e-12 and cov[0, 0] == 1.0, "median of 64 cells, full cover")
    s2 = s.copy(); s2[:5, :] = np.nan                     # 24 of 64 left: 37.5%
    med, cov = sc.block_median(s2, 0.0, 2.0, 0.25, 0.0, 2.0, 2.0, 1, 1)
    check(abs(cov[0, 0] - 24 / 64) < 1e-12 and abs(med[0, 0] - np.median(s2[5:])) < 1e-12,
          "cover share and median of the finite cells")
    med, cov = sc.block_median(s, 100.0, 2.0, 0.25, 0.0, 2.0, 2.0, 1, 1)
    check(np.isnan(med[0, 0]) and cov[0, 0] == 0, "no overlap -> NaN, cover 0")


def test_dem_at_points(d):
    print("dem_at_points")
    dem = np.array([[1.0, 2.0], [3.0, np.nan]])          # rows north to south
    v, st = sc.dem_at_points(dem, 0.0, 2.0, 1.0, np.array([0.5, 1.5, 1.5, 5.0, 0.75]),
                             np.array([1.5, 1.5, 0.5, 1.0, 1.25]))
    check(abs(v[0] - 1.0) < 1e-12 and abs(v[1] - 2.0) < 1e-12, "at a cell centre: that cell")
    check(st[2] == "dem_cell_without_value" and np.isnan(v[2]), "point in a NaN cell: not compared")
    check(st[3] == "outside_dem_extent", "point outside: flagged, not dropped")
    check(abs(v[4] - (1.0 * 0.75 * 0.75 + 2.0 * 0.25 * 0.75 + 3.0 * 0.75 * 0.25) / (1 - 0.25 * 0.25))
          < 1e-12, "bilinear renormalised over the finite neighbours")


def test_dsm_offset(d):
    print("DSM survey, DEM = beach + 0.20 m")
    dem_p, dem = make_dem(d)
    lid = make_lidar(d, water_below=1.0)
    # waterline points: on the beach, with setup (beach column) 0.10 m high for c1,
    # tide-only rows for c2 0.30 m low
    cp = Path(d) / "contours.csv"
    rng = np.random.default_rng(1)
    with open(cp, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["source_file", "camera", "capture_time_utc", "tide_elevation_navd88",
                    "beach_elevation_navd88", "easting_utm19", "northing_utm19"])
        for k in range(40):
            cam = "c1" if k % 2 else "c2"
            day = f"2025-01-{18 + k % 6:02d}"
            E = rng.uniform(X0 + 5, X0 + 35, 50)
            N = rng.uniform(Y0 + 10, Y0 + 110, 50)
            for e, n in zip(E, N):
                zt = beach(e, n)
                if cam == "c1":
                    w.writerow([f"f{k}", cam, day + "T15:00:00+00:00", f"{zt - 0.3:.4f}",
                                f"{zt + 0.10:.4f}", f"{e:.3f}", f"{n:.3f}"])
                else:
                    w.writerow([f"f{k}", cam, day + "T15:00:00+00:00", f"{zt - 0.3:.4f}", "",
                                f"{e:.3f}", f"{n:.3f}"])
    out = Path(d) / "out_dsm"
    h, log = quiet(sc.compare, str(dem_p), str(lid), "dsm", "INDEPENDENT", "synthetic test",
                   name="synth_dsm", output_dir=str(out), survey_date="2025-01-23",
                   photo_dates=("2025-01-18", "2025-01-23"), contours=str(cp))
    check(abs(h["median"] - OFFSET) < 0.005, f"median DEM - survey = {h['median']:+.4f} (+0.20)")
    # cells the lidar covers completely: exactly +0.20 (a part-covered cell at the water edge
    # takes the median of its dry part, biased up the slope by up to slope x cell / 4)
    g_d, _ = read_asc(out / "synth_dsm_dem_minus_survey.asc")
    full_cells = np.isfinite(g_d) & ((dem - OFFSET) > 1.0 + 0.05 + 1e-6)
    check(full_cells.sum() > 500 and np.allclose(g_d[full_cells], OFFSET, atol=2e-4),
          f"{int(full_cells.sum())} fully covered cells: exactly +0.200")
    check(h["nmad"] < 0.01 and abs(h["mean"] - OFFSET) < 0.01, "NMAD ~ 0, mean +0.20")
    check(h["within_0.10"] == 0.0 and h["within_0.50"] == 1.0, "shares within +/-0.10 / 0.50")
    # coverage: the lidar is masked below 1.0 m, the DEM goes below that
    zdem = dem - OFFSET
    full = np.isfinite(zdem) & (zdem > 1.0 + 0.05)
    check(h["dem_cells_compared"] >= full.sum() and h["dem_cells_compared"] < np.isfinite(dem).sum(),
          f"water-masked cells not compared ({h['dem_cells_compared']} of {h['dem_cells']})")
    check(h["dem_cells_no_survey"] > 0, "cells without survey counted")
    a = out / "synth_dsm_dem_minus_survey.asc"
    t = out / "synth_dsm_dem_minus_survey.tif"
    check(a.exists() and t.exists(), "difference .asc and .tif written")
    g, ha = read_asc(a)
    gt, x0, ytop, c = read_survey(t)
    check(np.allclose(g, gt, equal_nan=True, atol=1e-4) and abs(x0 - X0) < 1e-6
          and abs(ytop - (Y0 + 120)) < 1e-6, "GeoTIFF = .asc, same georeferencing")
    _, tags = read_tiff(t)
    keys = tags.get(34735) or ()
    check(4096 not in keys[4::4] and 3072 in keys[4::4], "no vertical datum key on a difference grid")
    wl = h["waterlines"]
    rows = list(csv.DictReader(open(out / "synth_dsm_comparison.csv")))
    c1 = [r for r in rows if r["quantity"] == "waterline_minus_survey" and r["group"] == "c1"][0]
    c2 = [r for r in rows if r["quantity"] == "waterline_minus_survey" and r["group"] == "c2"][0]
    check(abs(float(c1["median"]) - 0.10) < 0.005 and abs(float(c2["median"]) + 0.30) < 0.005,
          "waterlines: c1 (setup column) +0.10, c2 (tide only) -0.30")
    check("beach_elevation_navd88 for" in wl["elevation_used"], "mixed setup columns reported")
    ws = wl.get("without_setup") or {}
    check(abs(ws["by_camera"]["c1"]["median"] + 0.30) < 0.005 and abs(wl["by_camera"]["c1"]["median"] - 0.10) < 0.005
          and "overstates" in ws["how"],
          "C = 0 sensitivity: c1 at its still-water level -0.30 (with the setup +0.10), per camera, said to "
          "overstate without a calibration")
    txt = (out / "synth_dsm_comparison.txt").read_text()
    check("SENSITIVITY: WATERLINES vs SURVEY WITHOUT THE WAVE SETUP" in txt, "the sensitivity is in the txt")
    check(txt.splitlines()[2] == "LABEL      : INDEPENDENT" and "WHY        : synthetic test" in txt,
          "label at the top of the txt, exactly as passed")
    check("DOWNGRADE" not in txt, "no downgrade warning when nothing contradicts the label")
    check(any(r["group_type"] == "elevation_band" for r in rows)
          and any(r["group_type"] == "distance" for r in rows)
          and any(r["group_type"] == "camera" for r in rows)
          and any(r["group_type"] == "dem_frames" for r in rows), "groups: band, distance, camera, frames")
    check((out / "synth_dsm_comparison.png").exists() and (out / "synth_dsm_comparison.json").exists(),
          "figure and json written")
    js = json.loads((out / "synth_dsm_comparison.json").read_text())
    check(js["label"] == "INDEPENDENT" and abs(js["median"] - OFFSET) < 0.005, "json headline")
    # cells the lidar covers only partly (its water edge) are counted, and the fully covered ones
    # alone give the exact offset
    check(js["dem_cells_compared_partial_cover"] > 0 and abs(js["full_cover"]["median"] - OFFSET) < 2e-4
          and js["full_cover"]["n"] + js["dem_cells_compared_partial_cover"] == js["dem_cells_compared"],
          f"{js['dem_cells_compared_partial_cover']} partly covered cells counted; full cover only: "
          f"{js['full_cover']['median']:+.4f}")
    check("partly covered" in txt and any(r["group_type"] == "survey_cover" for r in rows),
          "partial cover in the coverage table and as a group")
    # the synthetic lidar states no vertical datum (as the 2025 YSMP files): said so
    check("geoid model unknown" in (js.get("vertical_datum_note") or ""), "no vertical key: NAVD88 assumed, said so")
    lid2 = Path(d) / "lidar_navd88.tif"
    write_geotiff(lid2, np.zeros((4, 4)), X0, Y0, 1.0, 6348, "test", vertical_epsg=5703)
    check(sc.survey_vertical_text(lid2) == ("EPSG:5703 (NAVD88 height)", ""), "a vertical key 5703 is read")
    # the CLI end to end
    rc = os.system(f"{sys.executable} {HERE / 'survey_compare.py'} --dem {dem_p} --survey {lid} "
                   f"--name cli --output-dir {out} --label CROSS-VALIDATED --why 'cli test' "
                   f"--no-plot > {out / 'cli.log'} 2>&1")
    check(rc == 0 and (out / "cli_comparison.txt").exists(), "CLI runs (survey type from the extension)")


def test_points(d):
    print("point surveys, DEM = beach + 0.20 m")
    dem_p, dem = make_dem(d)
    rng = np.random.default_rng(2)
    E = rng.uniform(X0 + 3, X0 + 77, 30)
    N = rng.uniform(Y0 + 3, Y0 + 110, 30)
    off_E, off_N = np.array([X0 - 50, X0 + 500]), np.array([Y0 + 20, Y0 + 20])
    emlid = Path(d) / "rtk_emlid.csv"
    with open(emlid, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Name", "Code", "Code description", "Easting", "Northing", "Elevation",
                    "Description", "CS name"])
        cs = "NAD83(2011) / UTM zone 19N + NAVD88(GEOID18) height"
        for k, (e, n) in enumerate(zip(E, N)):
            w.writerow([str(k + 1), "P", "Point", f"{e:.3f}", f"{n:.3f}", f"{beach(e, n):.4f}",
                        "Transect" if k % 3 else "High water wrack", cs])
        for k, (e, n) in enumerate(zip(off_E, off_N)):
            w.writerow([f"far{k}", "P", "Point", f"{e:.3f}", f"{n:.3f}", "12.0", "Check", cs])
    out = Path(d) / "out_pts"
    h, log = quiet(sc.compare, str(dem_p), str(emlid), "points", "INDEPENDENT", "rtk test",
                   name="synth_rtk", output_dir=str(out))
    check(abs(h["median"] - OFFSET) < 1e-3 and h["nmad"] < 1e-3, f"Emlid: median {h['median']:+.4f}")
    check(h["points"] == 32 and h["points_outside_dem"] == 2, "2 points outside: counted")
    rows = list(csv.DictReader(open(out / "synth_rtk_points.csv")))
    far = [r for r in rows if r["name"].startswith("far")]
    check(len(rows) == 32 and all(r["status"] == "outside_dem_extent" for r in far),
          "every point in the points CSV, the far ones flagged")
    txt = (out / "synth_rtk_comparison.txt").read_text()
    check("far0" in txt and "far1" in txt, "points not compared listed by name in the txt")
    groups = list(csv.DictReader(open(out / "synth_rtk_comparison.csv")))
    descs = {r["group"] for r in groups if r["group_type"] == "description"}
    check({"Transect", "High water wrack"} <= descs, "wrack reported separately")
    check("DATUM WARNING" not in txt, "NAVD88 CS name: no datum warning")

    # CIRN targets (no header) and a generic x,y,z CSV
    cirn = Path(d) / "2025-01-23_Marconi_Extrinsic_Targets_c1_xyz.csv"
    with open(cirn, "w") as f:
        for k, (e, n) in enumerate(zip(E[:10], N[:10])):
            f.write(f"{100 + k},{e:.3f},{n:.3f},{beach(e, n):.4f}\n")
    p = sc.read_points(cirn)
    check(p["format"].startswith("CIRN") and len(p["E"]) == 10 and p["name"][0] == "100",
          "CIRN num,E,N,Z read without a header")
    gen = Path(d) / "generic.csv"
    with open(gen, "w") as f:
        f.write("id,x,y,z\n")
        for k, (e, n) in enumerate(zip(E[:5], N[:5])):
            f.write(f"p{k},{e:.3f},{n:.3f},{beach(e, n):.4f}\n")
    p = sc.read_points(gen)
    check(len(p["E"]) == 5 and p["name"][0] == "p0", "generic x,y,z CSV")
    h2, _ = quiet(sc.compare, str(dem_p), str(gen), "points", "INDEPENDENT", "generic",
                  name="synth_gen", output_dir=str(out), plot=False)
    check(abs(h2["median"] - OFFSET) < 1e-3, "generic CSV: +0.20")

    # a survey not on NAVD88
    bad = Path(d) / "rtk_ellipsoid.csv"
    with open(bad, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Name", "Easting", "Northing", "Elevation", "Description", "CS name"])
        w.writerow(["1", f"{E[0]:.3f}", f"{N[0]:.3f}", "1.0", "Transect", "NAD83(2011) / UTM zone 19N"])
    h3, _ = quiet(sc.compare, str(dem_p), str(bad), "points", "INDEPENDENT", "datum",
                  name="synth_bad", output_dir=str(out), plot=False)
    check(any("DATUM WARNING" in w for w in h3["warnings"]), "CS name without NAVD88: datum warning")

    # the GCPs the calibration was solved from: PARTLY-CIRCULAR
    eo = Path(d) / "CACO03_c1_20250123_EO.yaml"
    shutil.copy(HERE / "calibration" / "CACO03_c1_20250123_EO.yaml", eo)
    h4, _ = quiet(sc.compare, str(dem_p), str(cirn), "points", "INDEPENDENT", "gcp test",
                  name="synth_gcp", output_dir=str(out), plot=False, camera_eo={"c1": eo})
    check(h4["suggested_label"] == "PARTLY-CIRCULAR" and h4["downgrade_warning"],
          "GCPs + calibration dated their day -> PARTLY-CIRCULAR warning")
    h5, _ = quiet(sc.compare, str(dem_p), str(cirn), "points", "PARTLY-CIRCULAR", "gcp test",
                  name="synth_gcp2", output_dir=str(out), plot=False, camera_eo={"c1": eo})
    check(not h5["downgrade_warning"], "labelled PARTLY-CIRCULAR already: no warning")


def test_label_checks(d):
    print("label checks")
    dem_p, dem = make_dem(d)
    lid = make_lidar(d)
    out = Path(d) / "out_label"
    io_src = HERE / "calibration" / "CACO05_c1_20240801_IO.yaml"
    cal = Path(d) / "cal"
    cal.mkdir(exist_ok=True)
    shutil.copy(io_src, cal / "CACO05_c1_20240801_IO.yaml")
    fit = cal / "CACO05_c1_2025-01-18_to_2025-01-23_lidar_EO.yaml"
    fit.write_text("x: 420088.0\ny: 4638319.0\nz: 20.8\nazimuth: 72.0\ntilt: 71.6\nroll: -6.6\n"
                   f"# fit_eo_to_survey.py: CACO03_c1_20250123_EO.yaml with azimuth +1.0 fitted to "
                   f"{lid.name} using 40 waterline frames 2025-01-18 to 2025-01-23; position kept.\n")
    h, _ = quiet(sc.compare, str(dem_p), str(lid), "dsm", "INDEPENDENT", "claimed independent",
                 name="lab1", output_dir=str(out), camera_eo={"c1": fit}, survey_date="2025-01-23")
    check(h["suggested_label"] == "CIRCULAR" and h["downgrade_warning"],
          "EO fitted to this survey -> CIRCULAR downgrade warning")
    txt = (out / "lab1_comparison.txt").read_text()
    check(txt.splitlines()[2] == "LABEL      : INDEPENDENT" and "DOWNGRADE WARNING" in txt,
          "label still printed as passed, with the warning")
    # apply_pointing_correction.py: the fit's change carried to another period (*_corr_EO.yaml)
    corr = cal / "CACO05_c1_20250219_corr_EO.yaml"
    corr.write_text("x: 420088.0\ny: 4638319.0\nz: 20.8\nazimuth: 72.0\ntilt: 71.6\nroll: -6.6\n"
                    f"# apply_pointing_correction.py: CACO04_c1_20250219_EO.yaml + (+1.000, +0.000, +0.000) deg "
                    f"from {fit.name}; 2024-10-01..2025-01-17\n")
    h, _ = quiet(sc.compare, str(dem_p), str(lid), "dsm", "INDEPENDENT", "carried", name="lab1c",
                 output_dir=str(out), camera_eo={"c1": corr}, plot=False)
    f_ = [f for f in h["findings"] if f["suggested"] == "CIRCULAR"]
    check(h["suggested_label"] == "CIRCULAR" and f_ and f_[0]["kind"] == "chain" and "carries" in f_[0]["text"],
          "a pointing carried from a fit to this survey (*_corr_EO.yaml) -> CIRCULAR, a chain fact")
    lost = cal / "CACO05_c2_20250219_corr_EO.yaml"
    lost.write_text(corr.read_text().replace(fit.name, "CACO05_c2_gone_lidar_EO.yaml"))
    shutil.copy(HERE / "calibration" / "CACO05_c2_20240801_IO.yaml", cal / "CACO05_c2_20240801_IO.yaml")
    sh = sc.audit_label("INDEPENDENT", str(lid), "dsm", "2025-01-23", sc.load_cameras({"c2": lost}, {}), None,
                        str(dem_p))
    check(any(x[0] == "CIRCULAR" and x[2] == "rule" and "cannot be" in x[1] for x in sh),
          "its fit not on this computer: CIRCULAR as a rule that can misfire")
    other = cal / "CACO05_c1_2025-03-01_to_2025-03-06_lidar_EO.yaml"
    other.write_text(fit.read_text().replace(lid.name, "2025005FA_Marconi_Mar_YSMP_Lidar_DSM_25cm.tif"))
    h, _ = quiet(sc.compare, str(dem_p), str(lid), "dsm", "INDEPENDENT", "other survey",
                 name="lab2", output_dir=str(out), camera_eo={"c1": other}, plot=False)
    check(not h["downgrade_warning"] and any("another survey" in n for n in h["notes"]),
          "EO fitted to another survey: a note, not a downgrade")
    h, _ = quiet(sc.compare, str(dem_p), str(lid), "dsm", "INDEPENDENT", "env", name="lab3",
                 output_dir=str(out), envelope_source=str(lid), plot=False)
    check(h["suggested_label"] == "PARTLY-CIRCULAR", "envelope placed with this survey -> PARTLY-CIRCULAR")
    # provenance.json next to the DEM naming the survey as the envelope source
    sub = Path(d) / "prov" / "dem"
    sub.mkdir(parents=True, exist_ok=True)
    shutil.copy(dem_p, sub / "synth_dem.asc")
    (sub.parent / "provenance.json").write_text(json.dumps(
        {"detect": {"search_envelope_source": str(lid)}, "compare": {"survey": str(lid)}}))
    h, _ = quiet(sc.compare, str(sub / "synth_dem.asc"), str(lid), "dsm", "INDEPENDENT", "prov",
                 name="lab4", output_dir=str(out), plot=False)
    check(h["suggested_label"] == "PARTLY-CIRCULAR", "provenance envelope = this survey -> PARTLY-CIRCULAR")
    # time gap
    h, _ = quiet(sc.compare, str(dem_p), str(lid), "dsm", "INDEPENDENT", "gap", name="lab5",
                 output_dir=str(out), survey_date="2025-02-10", photo_dates=("2025-01-18", "2025-01-23"),
                 plot=False)
    check(h["days_outside_photo_window"] == 18 and any("TIME GAP" in w for w in h["warnings"]),
          "survey 18 days after the photos: TIME GAP warning")
    o, m, t = sc.time_gap("2026-09-29", ("2026-09-20", "2026-10-04"))
    check(o == 0 and abs(m - 2.0) < 1e-9 and "inside" in t, "inside the window, 2 days after its middle")


def test_transects(d):
    print("waterlines on RTK transects: one value per frame, the RTK interpolated along the profile")
    dem_p, dem = make_dem(d)
    rtk = Path(d) / "rtk_transects.csv"
    cs = "NAD83(2011) / UTM zone 19N + NAVD88(GEOID18) height"
    # three cross-shore transects 30 m apart (alongshore = N), points every 4 m; a wrack line
    with open(rtk, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Name", "Code", "Easting", "Northing", "Elevation", "Description", "CS name"])
        k = 0
        for n0 in (Y0 + 30, Y0 + 60, Y0 + 90):
            for e in np.arange(X0 + 4, X0 + 60, 4.0):
                k += 1
                w.writerow([str(k), "T", f"{e:.3f}", f"{n0:.3f}", f"{beach(e, n0):.4f}", "Transect", cs])
        for n in np.arange(Y0 + 20, Y0 + 100, 10.0):
            k += 1
            w.writerow([str(k), "W", f"{X0 + 10:.3f}", f"{n:.3f}", f"{beach(X0 + 10, n) + 0.5:.4f}",
                        "High water wrack", cs])
    # 12 frames: each an alongshore line at its own elevation, 0.15 m LOW (its elevation below
    # the beach where it lies), points every 0.5 m; plus one frame that misses the transects
    cp = Path(d) / "contours_rtk.csv"
    with open(cp, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["source_file", "camera", "capture_time_utc", "tide_elevation_navd88",
                    "beach_elevation_navd88", "easting_utm19", "northing_utm19"])
        for k in range(13):
            z_line = 1.0 + 0.1 * k
            e = X0 + (3.0 - z_line) / 0.05                  # where the beach is at z_line
            ns = np.arange(Y0 + 5, Y0 + 115, 0.5) if k < 12 else np.arange(Y0 + 200, Y0 + 210, 0.5)
            for n in ns:
                w.writerow([f"f{k}", "c1", "2026-09-29T15:00:00+00:00", f"{z_line - 0.4:.4f}",
                            f"{z_line - 0.15:.4f}", f"{e:.3f}", f"{n:.3f}"])
    out = Path(d) / "out_rtk"
    h, log = quiet(sc.compare, str(dem_p), str(rtk), "points", "CIRCULAR", "transect test", name="synth_tr",
                   output_dir=str(out), contours=str(cp), photo_dates=("2026-09-29", "2026-09-29"))
    wl = h["waterlines"]
    check(wl["method"].startswith("frames on RTK transects") and wl["n"] == 12 and wl["frames"] == 12,
          f"12 frames on the transects, one value each (n {wl['n']}); the frame off them not counted")
    check(abs(wl["median"] + 0.15) < 1e-3 and wl["nmad"] < 1e-3, f"waterline - RTK = -0.15 exactly ({wl['median']:+.4f})")
    check("wrack" in wl["how"] and (out / "synth_tr_waterline_frames.csv").exists(),
          "wrack points not used for the waterlines; per-frame values written")
    # isolated points (no transects): nearest point, not slope-corrected, and said so
    gen = Path(d) / "isolated.csv"
    gen.write_text("E,N,Z\n" + "".join(f"{X0 + 30:.3f},{Y0 + n:.3f},{beach(X0 + 30, 0):.4f}\n"
                                         for n in (20, 50, 80)))
    h2, _ = quiet(sc.compare, str(dem_p), str(gen), "points", "INDEPENDENT", "iso", name="synth_iso",
                  output_dir=str(out), contours=str(cp), plot=False)
    check("NOT slope-corrected" in h2["waterlines"]["how"], "isolated points: matched to the nearest line point, "
          "NOT slope-corrected (said)")
    # few values: no NMAD, valid JSON (no NaN)
    s1 = sc.stats(np.array([0.3]))
    check(s1["n"] == 1 and np.isnan(s1["nmad"]) and np.isnan(s1["p5"]), "one value: NMAD and p5/p95 not given")
    raw = (out / "synth_iso_comparison.json").read_text()
    check("NaN" not in raw and json.loads(raw) is not None, "the json holds no NaN (null instead)")


def test_no_overlap(d):
    print("no overlap")
    dem_p, dem = make_dem(d)
    lid = make_lidar(d)
    far = Path(d) / "far_dem.asc"
    g, h = read_asc(dem_p)
    write_asc(far, g, X0 + 5000, Y0, 2.0)               # 5 km east of the survey
    out = Path(d) / "out_far"
    h1, _ = quiet(sc.compare, str(far), str(lid), "dsm", "INDEPENDENT", "far", name="far",
                  output_dir=str(out))
    check(h1["n"] == 0 and h1["dem_cells_compared"] == 0 and (out / "far_comparison.png").exists(),
          "DSM not overlapping: n = 0, report and figure still written")
    pts = Path(d) / "far_pts.csv"
    pts.write_text("E,N,Z\n420000,4630000,1.0\n")
    h2, _ = quiet(sc.compare, str(dem_p), str(pts), "points", "INDEPENDENT", "far", name="far_pts",
                  output_dir=str(out))
    check(h2["n"] == 0 and h2["points_outside_dem"] == 1, "point survey off the DEM: n = 0, listed")


def test_py38():
    print("Python 3.8 syntax")
    import ast
    for f in ("survey_compare.py", "test_survey_compare.py"):
        src = (HERE / f).read_text()
        try:
            ast.parse(src, feature_version=(3, 8))
        except TypeError:                        # Python 3.7 has no feature_version
            ast.parse(src)
        check(True, f"{f} parses as Python 3.8")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="write the outputs here and keep them")
    args = ap.parse_args()
    d = args.keep or tempfile.mkdtemp(prefix="survey_compare_test_")
    Path(d).mkdir(parents=True, exist_ok=True)
    try:
        test_block_median(d)
        test_dem_at_points(d)
        test_dsm_offset(d)
        test_points(d)
        test_label_checks(d)
        test_transects(d)
        test_no_overlap(d)
        test_py38()
    finally:
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    print("all survey_compare tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

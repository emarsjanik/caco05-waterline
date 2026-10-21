#!/usr/bin/env python3
"""
Tests For evaluate_setup_model.py: Which Setup Model, On Synthetic Surveys
=============================================================================
Plain python3 (no pytest needed on the station):

    python3 test_evaluate_setup_model.py [--keep DIR]

WHAT. Three synthetic survey-date products, each a planar beach of its own
slope (0.13, 0.07, 0.10: the Jan 2025, Mar 2025 and Sep-Oct 2026 faces),
seen through the station's real c1 calibration; two lidar-like DSMs (the
second dry only above -0.5 m: the lidar's water edge) and one RTK file of
cross-shore transects stopping at 0.0 m. Every line is planted where the
beach meets still water + a TRUE setup, then given pixels by projecting it
into the photo:
  * truth C (C = 0.045 on every beach): the C refit finds 0.045 in-sample
    and held out, Stockdon's refit does not hold across the slopes, and
    the recommendation is to keep C;
  * truth Stockdon (k = 0.35, H0 deshoaled from 21 m, beta_f the beach's):
    the k refit finds 0.35 (and the paper's k scores near the noise), the
    C refit does worse held out, and the recommendation is Stockdon with
    that k;
  * truth C = 0.037, but each survey's lines carry a level offset that
    lines up with its slope (+0.10 m on the 0.13 beach, -0.08 m on the
    0.07 one: a pointing or still-water datum error): Stockdon wins the
    absolute cross-validation, loses it with a free offset per survey and
    camera, and the recommendation is to keep C; each held-out coefficient
    is the fit on the other two surveys (a run on those two alone); with
    no DEM in a product every Stockdon frame takes C / k (the not-deshoaled
    one then IS the station's C);
  * one test of the recommendation: a mean gain AND a majority of the
    held-out surveys, a survey on which Stockdon keeps too few frames
    counted against it, every survey scored;
  * two surveys only: no recommendation, the station keeps C;
  * two planted gross failures (lines 1.5 m off) are counted and left out;
    lines whose plane lies below the DSM's water edge or the transects'
    lowest shot are not scored (the survey-floor rule);
  * the ground point at any plane, from the two planes the script projects
    at, equals georectify.pixel_to_ground at that plane;
  * a BLIND date (2024-10-23, survey_products.BLIND_DATES) in the table,
    with a product folder and a survey file: neither is opened (an audit
    hook watches every open of this process); nor is a survey file of
    another date's row NAMED like a survey of the blind date;
  * the --json record and the exit codes (0 a report, 4 nothing scored).

STATION FILES. None are read (the calibration files of this repository give
the camera geometry).
"""

import argparse
import contextlib
import csv
import io
import json
import os
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import setup_model as sm          # noqa: E402

FAILED = []
OPENED = []
IO_FILE = "CACO05_c1_20240801_IO.yaml"
EO_FILE = HERE / "calibration" / "CACO05_c1_20251113_EO-CV.yaml"
TAG = "marconi:hs_best/tp_s"


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILED.append(msg)


def _audit(event, args):
    if event == "open" and args and isinstance(args[0], (str, bytes, os.PathLike)):
        OPENED.append(os.fsdecode(args[0]))


def camera():
    from georectify import load_extrinsics, load_intrinsics
    return load_intrinsics(HERE / "calibration" / IO_FILE), load_extrinsics(EO_FILE)


def write_asc(path, grid_ns, xll, yll, cell):
    """rows north to south, NaN -> nodata"""
    with open(path, "w") as f:
        f.write(f"ncols {grid_ns.shape[1]}\nnrows {grid_ns.shape[0]}\nxllcorner {xll:.3f}\nyllcorner {yll:.3f}\n"
                f"cellsize {cell}\nNODATA_value -9999\n")
        for r in grid_ns:
            f.write(" ".join("%.4f" % (v if np.isfinite(v) else -9999.0) for v in r) + "\n")


def make_product(root, sdir, date, name, beta, truth, kind, seed, gross=0, dsm_edge=None, low=(), offset=0.0):
    """One synthetic survey-date product (provenance, waterlines, DEM) and its survey file.
    truth: ('c', C) or ('stockdon', k). The beach: z = beta * (E0 - E), E0 40 m east of the camera. The first
    `gross` lines are planted 1.5 m LOW (a detection on the wet sand far seaward of the water's edge); `low`: the
    still-water levels of the last frames (low tides, below a survey's floor); `offset`: added to every line
    (a level error of the survey date, not a setup)."""
    from view_reproject import ground_to_pixel
    from georectify import pixel_to_ground
    io_, eo_ = camera()
    ec, nc = float(eo_[0]), float(eo_[1])
    e0 = ec + 40.0
    rng = np.random.RandomState(seed)
    pdir = root / date
    (pdir / "waterlines").mkdir(parents=True)
    (pdir / "dem").mkdir()
    rows = []
    t0 = int(datetime.strptime(date, "%Y-%m-%d").replace(tzinfo=timezone.utc).timestamp()) - 3 * 86400
    nframes = 30
    for i in range(nframes):
        sw = round(float(rng.uniform(-0.6, 1.2)), 4)
        if i >= nframes - len(low):
            sw = low[i - (nframes - len(low))]
        hs, tp = round(float(rng.uniform(0.4, 2.5)), 2), round(float(rng.uniform(6.0, 12.0)), 1)
        if truth[0] == "c":
            s_true = float(sm.setup_c(hs, tp, truth[1]))
        else:
            h0 = float(sm.deshoal(hs, tp, 21.0))
            s_true = float(sm.stockdon_setup(h0, tp, beta, truth[1]))
        z_line = sw + s_true + float(rng.normal(0, 0.02)) + (-1.5 if i < gross else 0.0) + offset
        e_line = e0 - z_line / beta
        N = nc + np.arange(5.0, 100.0, 0.25)
        U, V, ok = ground_to_pixel(np.full(len(N), e_line), N, np.full(len(N), z_line), io_, eo_)
        U, V = U[ok] + rng.normal(0, 0.3, ok.sum()), V[ok] + rng.normal(0, 0.3, ok.sum())
        ep = t0 + i * 7200
        setup = round(float(sm.setup_c(hs, tp, 0.037)), 4)
        E, Nn = pixel_to_ground(U, V, sw + setup, io_, eo_)
        tu = datetime.fromtimestamp(ep, tz=timezone.utc).isoformat()
        for u, v, e, n in zip(U, V, E, Nn):
            rows.append([f"{ep}.c1.timex", "c1", tu, ep, f"{u:.3f}", f"{v:.3f}", sw, hs, tp, setup,
                         round(sw + setup, 4), TAG, f"{e:.3f}", f"{n:.3f}"])
    ground = pdir / "waterlines" / "contour_points_ground_filtered.csv"
    with open(ground, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["source_file", "camera", "capture_time_utc", "capture_epoch", "pixel_column", "pixel_row",
                    "tide_elevation_navd88", "offshore_hs_m", "offshore_tp_s", "setup_correction_m",
                    "beach_elevation_navd88", "wave_source", "easting_utm19", "northing_utm19"])
        w.writerows(rows)
    # the product's DEM (2 m) over the intertidal, the survey (1 m DSM or transects)
    xll, yll = ec - 10.0, nc - 20.0
    ex = xll + (np.arange(60) + 0.5) * 2.0
    dem = np.tile(beta * (e0 - ex), (80, 1))
    dem[(dem < -1.8) | (dem > 3.0)] = np.nan
    write_asc(pdir / "dem" / f"{date}_dem.asc", dem, xll, yll, 2.0)
    if kind == "dsm":
        sx = xll + (np.arange(120) + 0.5) * 1.0
        g = np.tile(beta * (e0 - sx), (160, 1))
        if dsm_edge is not None:
            g[g < dsm_edge] = np.nan
        spath = sdir / f"{name}.asc"
        write_asc(spath, g, xll, yll, 1.0)
    else:
        spath = sdir / f"{name}.csv"
        with open(spath, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Name", "Code", "Easting", "Northing", "Elevation", "Description", "CS name"])
            k = 0
            for nt in (25.0, 45.0, 65.0, 85.0):
                for z in np.arange(3.0, -0.001, -0.2):
                    k += 1
                    w.writerow([f"P{k}", "", f"{e0 - z / beta:.3f}", f"{nc + nt:.3f}", f"{z:.3f}", "Transect",
                                "NAD83(2011) / UTM zone 19N + NAVD88 height"])
    prov = {"date": date, "status": "complete", "blind": None,
            "cameras": {"c1": {"calibration": {"eo_path": str(EO_FILE), "io_file": IO_FILE}}},
            "filter": {"file": str(ground)},
            "dem": {"files": {"dem": str(pdir / "dem" / date) + "_dem"}, "built_from": {"filters": {}}}}
    (pdir / "provenance.json").write_text(json.dumps(prov))
    return {"date": date, "name": name, "path": spath.name, "survey_type": "dsm" if kind == "dsm" else "points",
            "survey_date": date, "label": "INDEPENDENT", "why": "synthetic"}


def fixture(d, truth, gross=0, offsets=(0.0, 0.0, 0.0)):
    root, sdir = Path(d) / "products", Path(d) / "surveys"
    sdir.mkdir(parents=True)
    root.mkdir()
    rows = [make_product(root, sdir, "2031-01-23", "jan_dsm", 0.13, truth, "dsm", 1, gross=gross,
                         offset=offsets[0]),
            make_product(root, sdir, "2031-03-06", "mar_dsm", 0.07, truth, "dsm", 2, dsm_edge=-0.5,
                         low=(-1.1, -1.0, -0.95), offset=offsets[1]),
            make_product(root, sdir, "2031-09-29", "rtk", 0.10, truth, "transects", 3, offset=offsets[2])]
    # a BLIND date: its row, a survey file and a product folder, none of which may be opened
    (sdir / "held_blind_points.csv").write_text("Name,Easting,Northing,Elevation\n")
    (root / "2024-10-23" / "waterlines").mkdir(parents=True)
    (root / "2024-10-23" / "provenance.json").write_text(json.dumps({"blind": {"held_by": "someone"}}))
    (root / "2024-10-23" / "waterlines" / "contour_points_ground_filtered.csv").write_text("x\n")
    rows.append({"date": "2024-10-23", "name": "oct_points", "path": "held_blind_points.csv", "survey_type": "points",
                 "survey_date": "2024-10-23", "label": "PARTLY-CIRCULAR", "why": "held"})
    # a row of a date with a product whose survey file is NAMED like a survey of the blind date (a stand-in,
    # written here): never opened either
    (sdir / "2024-10-23_lookalike_GCP.csv").write_text("Name,Easting,Northing,Elevation\n")
    rows.append({"date": "2031-09-29", "name": "lookalike", "path": "2024-10-23_lookalike_GCP.csv",
                 "survey_type": "points", "survey_date": "2031-09-29", "label": "INDEPENDENT", "why": "synthetic"})
    table = Path(d) / "surveys.csv"
    with open(table, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["date", "name", "path", "survey_type", "survey_date", "label", "why"])
        w.writeheader()
        w.writerows(rows)
    return root, sdir, table


def run_eval(root, sdir, table, extra=()):
    import evaluate_setup_model as ev
    out = io.StringIO()
    js = Path(str(root) + ".json")
    with contextlib.redirect_stdout(out):
        rc = ev.main(["--output-root", str(root), "--surveys", str(table), "--survey-dirs", str(sdir),
                      "--json", str(js)] + list(extra))
    res = json.loads(js.read_text()) if js.exists() else {}
    if js.exists():
        js.unlink()
    return rc, out.getvalue(), res


def test_two_planes():
    print("the ground point at any plane from the two planes evaluate_setup_model.py projects at")
    from georectify import pixel_to_ground
    io_, eo_ = camera()
    rng = np.random.RandomState(5)
    U, V = rng.uniform(0, 2447, 2000), rng.uniform(900, 2047, 2000)
    sw = 0.3
    Ea, Na = pixel_to_ground(U, V, sw, io_, eo_)
    Eb, Nb = pixel_to_ground(U, V, sw + 1.0, io_, eo_)
    worst = 0.0
    for s in (-0.4, 0.27, 0.83, 1.6):
        Ex, Nx = pixel_to_ground(U, V, sw + s, io_, eo_)
        ok = np.isfinite(Ex) & np.isfinite(Ea) & np.isfinite(Eb)
        worst = max(worst, float(np.max(np.abs(Ea[ok] + s * (Eb[ok] - Ea[ok]) - Ex[ok]))),
                    float(np.max(np.abs(Na[ok] + s * (Nb[ok] - Na[ok]) - Nx[ok]))))
    check(worst < 1e-6, f"straight rays: interpolated = pixel_to_ground to {worst:.1e} m (setups -0.4 to 1.6 m)")


def test_truth_c(d):
    print("truth C = 0.045 on beaches of slope 0.13, 0.07, 0.10")
    root, sdir, table = fixture(Path(d) / "truth_c", ("c", 0.045), gross=2)
    del OPENED[:]
    rc, out, res = run_eval(root, sdir, table)
    blind_opened = [p for p in OPENED if "2024-10-23" in p or "held_blind" in p]
    check(rc == 0 and "2024-10-23 oct_points: BLIND" in out and not blind_opened,
          f"the BLIND row is skipped and nothing of it is opened ({len(OPENED)} files opened, "
          f"{len(blind_opened)} of the blind date)")
    check("2031-09-29 lookalike: 2024-10-23_lookalike_GCP.csv is named like a survey of a BLIND date" in out
          and not any("lookalike" in p for p in OPENED),
          "a row of another date whose survey file is named like a survey of the blind date: skipped, never opened")
    sv = {s["name"]: s for s in res.get("surveys", [])}
    check(set(sv) == {"jan_dsm", "mar_dsm", "rtk"} and all(s["used"] for s in sv.values()),
          "three surveys scored: " + ", ".join(f"{k} {v['scored']} of {v['frames']}" for k, v in sorted(sv.items())))
    check(sv.get("jan_dsm", {}).get("gross_failures") == 2, "the two lines planted 1.5 m low: 2 gross failures, "
          "left out of every score")
    check(sv.get("mar_dsm", {}).get("reached", 99) < sv.get("mar_dsm", {}).get("frames", 0)
          and sv.get("rtk", {}).get("reached", 99) < sv.get("rtk", {}).get("frames", 0),
          f"lines below the survey's floor are not scored (mar DSM, dry above -0.5 m: {sv['mar_dsm']['reached']} of "
          f"{sv['mar_dsm']['frames']} reach it; RTK stopping at 0.0 m: {sv['rtk']['reached']} of {sv['rtk']['frames']})")
    c = res["refit"]["c"]
    check(abs(c["all"] - 0.045) <= 0.002 and all(abs(v["coef"] - 0.045) <= 0.003 for v in c["loso"].values()),
          f"C refit {c['all']} in-sample, held out " + ", ".join(f"{k} {v['coef']}" for k, v in sorted(c["loso"].items())))
    check(all(v["rmse"] < 0.05 for v in c["loso"].values()),
          "held-out RMSE at the noise: " + ", ".join(f"{k} {v['rmse']:.3f} m" for k, v in sorted(c["loso"].items())))
    rec = res["recommendation"]
    check(rec["model"] == "c" and res["refit"]["stockdon"]["loso_rmse_mean"] > c["loso_rmse_mean"] + 0.01,
          f"keep C: Stockdon's cross-validated RMSE {res['refit']['stockdon']['loso_rmse_mean']:.3f} m against C's "
          f"{c['loso_rmse_mean']:.3f} m")
    check("Keep the station's model (SETUP_MODEL=\"\", C = 0.037)" in out and "LEAVE-ONE-SURVEY-OUT" in out,
          "the report says so")
    rc2, out2, res2 = run_eval(root, sdir, table, ["--dates", "2031-01-23", "2031-03-06"])
    check(rc2 == 0 and res2["recommendation"]["model"] == "c" and "a recommendation needs 3" in out2,
          "two surveys only: no recommendation, the station keeps C")
    rc3, out3, _ = run_eval(root, sdir, table, ["--dates", "2030-01-01"])
    check(rc3 == 4, "no product for the dates asked: exit 4")


def test_truth_stockdon(d):
    print("truth Stockdon k = 0.35 (H0 deshoaled from 21 m) on beaches of slope 0.13, 0.07, 0.10")
    root, sdir, table = fixture(Path(d) / "truth_s", ("stockdon", 0.35))
    rc, out, res = run_eval(root, sdir, table)
    s, c = res["refit"]["stockdon"], res["refit"]["c"]
    check(rc == 0 and abs(s["all"] - 0.35) <= 0.02 and all(abs(v["coef"] - 0.35) <= 0.03 for v in s["loso"].values()),
          f"k refit {s['all']} in-sample, held out " + ", ".join(f"{k} {v['coef']}" for k, v in sorted(s["loso"].items())))
    paper = res["fixed"]["Stockdon k = 0.35, H0 = Hs/Ks (the paper)"]["per_survey"]
    check(all(v["rmse"] < 0.05 for v in paper.values()) and all(
        set(v.get("beta_sources", {})) == {"dem"} for v in paper.values()),
          "the paper's k scores at the noise, every beta_f from the product's DEM: "
          + ", ".join(f"{k} {v['rmse']:.3f} m" for k, v in sorted(paper.items())))
    h0 = res["refit"]["stockdon_h0hs"]
    check(h0["all"] > s["all"], f"H0 = Hs (not deshoaled) needs a larger k: {h0['all']} against {s['all']}")
    rec = res["recommendation"]
    check(rec["model"] == "stockdon" and abs(rec["k"] - s["all"]) < 1e-9
          and c["loso_rmse_mean"] > s["loso_rmse_mean"] + 0.01,
          f"recommend Stockdon (k {rec['k']}): cross-validated RMSE {s['loso_rmse_mean']:.3f} m against C's "
          f"{c['loso_rmse_mean']:.3f} m")
    check(f'SETUP_MODEL="stockdon" and SETUP_STOCKDON_K="{rec["k"]:g}"' in out, "the report names the cron lines")


def test_offsets(d):
    print("truth C = 0.037, each survey's lines offset in level with its slope (+0.10 m at 0.13, -0.08 m at 0.07)")
    root, sdir, table = fixture(Path(d) / "offsets", ("c", 0.037), offsets=(0.10, -0.08, 0.0))
    rc, out, res = run_eval(root, sdir, table)
    rec = res["recommendation"]
    ab, fr = rec["absolute"], rec["free_offset"]
    check(rc == 0 and ab.get("gain_m", 0) >= 0.01 and len(ab["stockdon_better_on"]) >= 2,
          f"the absolute cross-validation alone would switch: Stockdon {ab.get('stockdon_rmse_mean', 0):.3f} m "
          f"against C's {ab.get('c_rmse_mean', 0):.3f} m, lower on {len(ab['stockdon_better_on'])} of 3")
    cf = res["refit"]["c"]["free_offset"]
    check(fr.get("gain_m", 1) < 0 and all(abs(v["coef"] - 0.037) <= 0.003 for v in cf["loso"].values()),
          f"with a free offset per survey and camera Stockdon loses ({fr.get('stockdon_rmse_mean', 0):.3f} m against "
          f"C's {fr.get('c_rmse_mean', 0):.3f} m) and C held out is the truth: "
          + ", ".join(f"{k} {v['coef']}" for k, v in sorted(cf["loso"].items())))
    check(rec["model"] == "c" and "Keep the station's model" in out and "does not hold with a free offset" in out,
          "keep C, and the report says why")
    # each held-out coefficient is the fit on the other two surveys
    rc2, _out2, res2 = run_eval(root, sdir, table, ["--dates", "2031-03-06", "2031-09-29"])
    pairs = [(fam, res["refit"][fam]["loso"]["jan_dsm"]["coef"], res2["refit"][fam]["all"]) for fam in res["refit"]]
    pairs += [(fam + " free offset", res["refit"][fam]["free_offset"]["loso"]["jan_dsm"]["coef"],
               res2["refit"][fam]["free_offset"]["all"]) for fam in res["refit"]]
    check(rc2 == 0 and all(a == b for _f, a, b in pairs) and
          any(a != res["refit"][f.split()[0]]["all"] for f, a, _b in pairs[:3]),
          "jan held out = the fit on mar and rtk alone: " + ", ".join(f"{f} {a}/{b}" for f, a, b in pairs))
    # no DEM in a product: C / k for every Stockdon frame
    (root / "2031-01-23" / "dem" / "2031-01-23_dem.asc").unlink()
    rc3, out3, res3 = run_eval(root, sdir, table, ["--dates", "2031-01-23"])
    paper = res3["fixed"]["Stockdon k = 0.35, H0 = Hs/Ks (the paper)"]["per_survey"]["jan_dsm"]
    station = res3["fixed"]["C = 0.037 (the station's)"]["per_survey"]["jan_dsm"]
    h0 = res3["refit"]["stockdon_h0hs"]["in_sample"]["jan_dsm"]
    check(rc3 == 0 and "no DEM in the product" in out3 and paper.get("beta_sources") == {"no DEM in the product": 30}
          and abs(h0["rmse"] - station["rmse"]) < 1e-9 and abs(h0["bias"] - station["bias"]) < 1e-9,
          f"no DEM: every frame beta_f = C / k; H0 = Hs then gives the station's C at any k (RMSE {h0['rmse']:.4f} m "
          f"= {station['rmse']:.4f} m)")


def test_verdict():
    print("one test of the recommendation (evaluate_setup_model.verdict): mean gain, majority, frames kept")
    import evaluate_setup_model as ev
    v = lambda c, s, keeps=True: {"n": 20, "c_rmse": c, "stockdon_rmse": s, "stockdon_keeps": keeps}  # noqa: E731
    ok1, _t1, _n1 = ev.verdict({"a": v(0.08, 0.05), "b": v(0.08, 0.05), "c": v(0.02, 0.03)}, 3, 0.01)
    ok2, t2, n2 = ev.verdict({"a": v(0.08, 0.05), "b": v(0.08, 0.05, False), "c": v(0.02, 0.03)}, 3, 0.01)
    ok3, _t3, _n3 = ev.verdict({"a": v(0.20, 0.05), "b": v(0.02, 0.03), "c": v(0.02, 0.03)}, 3, 0.01)
    ok4, _t4, _n4 = ev.verdict({"a": v(0.08, 0.05), "b": None, "c": v(0.02, 0.03)}, 3, 0.01)
    ok5, _t5, _n5 = ev.verdict({"a": v(0.03, 0.029), "b": v(0.03, 0.029), "c": v(0.02, 0.03)}, 3, 0.0)
    check(ok1 and not ok2 and n2["stockdon_better_on"] == ["a"] and "on b (a loss)" in t2 and not ok3 and not ok4
          and not ok5,
          "lower on 2 of 3 by 0.02 m: met; the same with b keeping too few frames: b is a loss; a large gain on 1 "
          "of 3: no majority; a survey not scored: not met; lower on 2 of 3 but higher on average: not met")


def test_py38():
    print("Python 3.8 syntax")
    import ast
    for f in ("evaluate_setup_model.py", "test_evaluate_setup_model.py"):
        src = (HERE / f).read_text(encoding="utf-8")
        try:
            ast.parse(src, feature_version=(3, 8))
        except TypeError:                        # Python 3.7 has no feature_version
            ast.parse(src)
        check(True, f"{f} parses as Python 3.8")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="keep the outputs, in a new folder inside this one")
    args = ap.parse_args()
    if args.keep:
        Path(args.keep).mkdir(parents=True, exist_ok=True)
    d = os.path.realpath(tempfile.mkdtemp(prefix="evaluate_setup_", dir=args.keep))
    print(f"test folder: {d}")
    sys.addaudithook(_audit)
    try:
        test_two_planes()
        test_truth_c(d)
        test_truth_stockdon(d)
        test_offsets(d)
        test_verdict()
        test_py38()
    finally:
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED")
        return 1
    print("all setup-model evaluation tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

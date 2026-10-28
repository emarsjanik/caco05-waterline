#!/usr/bin/env python3
"""
Tests For The Measured Setup (measured_setup.py, evaluate_measured_setup.py, --setup-source)
===============================================================================================
Plain python3 (no pytest needed on the station):

    python3 test_measured_setup.py

WHAT. Checked here, on synthetic data with a planted truth [SYNTH]:
  * the QC of a line, each rule alone; the median over the lines and its
    uncertainty (the sampling term, the between-line term, one line);
  * the sampling error of a 10-minute mean against Gaussian records of
    known spectrum (flat, empty below 0.005 Hz, level below 1/T);
  * the strictly-earlier DEM: never the stack's own day, never a window
    holding it, the latest earlier one, its age limit (exactly 3 days); a
    DEM file whose window is unknown is refused; what a row records of it;
    the MODEL-BUILT rule: a DEM whose _info.json counts measured (or
    unrecorded) frames, has no setup, or has no info file is refused,
    dem_from_contours.py's own _info.json read both ways; --dem's default
    (archive/dems_model when it exists);
  * the chained generations: the one-generation map measured on the scene
    (unit gain), the loop the rule forbids (drift, the erosion erased), the
    deployed chain on the model-only series (a constant offset, the change
    kept), and the _info.json rule alone (no drift, measurement starved);
  * the output file: appended without duplicates, rewritten only by
    --reprocess, another header refused; the extractor's reading of it
    (QC-passed collections only, the 300 s window);
  * end to end on a synthetic c2 scene (timestack_wave_height.py's: the
    real c2 calibration, a planar beach): the setup on the DEM of the day
    before, two lines, the waits (the GPS level waited for wherever it is
    the allowed source), the sharpness and still-water rules, the edge
    contrast logged and not a gate, and the INHERITANCE: a DEM 0.10 m too
    high or too low moves the setup by 0.10 / (1 - beta / tan(delta)), which
    profile_gain gives from the ray's geometry, and sigma_profile_m;
    timestack_wave_height.estimate_stack's rows are the same with the
    series hook on; --audit (the noise table, the DEMs checked again, the
    sample drawn with the edge);
  * extract_elevation_contours.py: the default (--setup-source model, or
    no option) writes the same file and log, byte for byte, as the code
    this branch started from (git) and as its pinned sha256; 'measured'
    and 'auto': which frames get which setup, both cameras, the 300 s
    window, the buoy outage, the slope-aware model's columns, the
    setup_source column last, the refusals before any frame is read; the
    default does not even import measured_setup.py;
  * the readers of a contour file with setup_source: georectify.py,
    waterline_consistency.py, evaluate_setup_model.read_lines,
    dem_from_contours.py (counts in _info.json and the log only when the
    column is there), survey_products.apply_setup (drops it: C's setup);
  * the cron: SETUP_SOURCE="" or "model" builds the same options;
    'measured' / 'auto' add the two options only with a real measured file
    (missing, empty or malformed: a WARNING), an unknown value a WARNING;
    while it is on, the model-only chain (the station's steps again with
    the model's options, into archive/dems_model); nothing runs the
    measurement;
  * evaluate_measured_setup.py: the fit of C with its day-jackknife
    interval (its coverage at 3 and 5 days, the t table); the still-water
    reference of a product and the corrected model (known answers); the
    measured setup on a survey's own surface (known setup); the lines
    scored against a survey with the model's, the corrected model's and
    the measured setup, with the paired intervals (truth C = 0.045, model
    0.037); a date no stack of which can be read skipped, not a crash; the
    window's sharpness from a period's own images; a BLIND row never
    opened;
  * Python 3.8 syntax (the station's python3).
Exit status 0 when every check passes.
"""

import ast
import contextlib
import csv
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import measured_setup as ms                 # noqa: E402
import timestack_wave_height as thw         # noqa: E402

failures = []
# The commit this branch started from (origin/main, 27 Oct 2026): the extractor's default must write what it wrote
BASE_COMMIT = "33221c8"
# sha256 of the default extractor's output (--waves --setup-coef 0.037) on extract_fixture(), made with BASE_COMMIT's
# extractor; the fixture is exact decimals, so the file is the same on any platform
PINNED_DEFAULT_SHA256 = "779d309b1e1d7935954a259bd52840fa7799f66a2a36b3800da2bbfa4ef1d72f"


def check(name, cond, detail=None):
    print(f"  {'ok ' if cond else 'FAIL'} {name}" + (f"  ({detail})" if detail is not None and not cond else ""))
    if not cond:
        failures.append(name)


def quiet(fn, *a, **k):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = fn(*a, **k)
    return out, buf.getvalue()


def run(args, cwd, env=None):
    r = subprocess.run([sys.executable] + [str(a) for a in args], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True, cwd=str(cwd), env=env)
    return r.returncode, r.stdout


MODEL_INFO = {"setup_share": 1.0}            # what dem_from_contours.py's _info.json says of a model-built DEM


def write_info(dem_path, info=MODEL_INFO):
    """The _info.json beside a series DEM (dem_from_contours.py --series-dir writes one): None writes none."""
    p = Path(dem_path)
    if info is not None and p.name.endswith("_dem.asc"):
        p.with_name(p.name[:-len("_dem.asc")] + "_info.json").write_text(json.dumps(info))
    return p


def write_plane_dem(path, beta, off=0.0, zlo=-1.5, zhi=2.6, info=MODEL_INFO, sigma=None):
    """The synthetic scene's planar beach z = -beta x (x seaward of the Jan 2025 0 m contour) as a 2 m DEM
    between zlo and zhi, raised by `off` (a DEM off the beach by that much); its _info.json (model-built by
    default) and, with sigma, a _sigma.asc of that constant."""
    x0, y0, cell, nx, ny = 419900.0, 4638300.0, 2.0, 200, 200
    XX, YY = np.meshgrid(x0 + (np.arange(nx) + 0.5) * cell, (y0 + (np.arange(ny) + 0.5) * cell)[::-1])
    z = -beta * thw.seaward_distance(XX, YY)
    bad = (z < zlo) | (z > zhi)
    z = z + off
    sg = np.where(bad, -9999, sigma if sigma is not None else 0.0)
    z[bad] = -9999
    head = f"ncols {nx}\nnrows {ny}\nxllcorner {x0}\nyllcorner {y0}\ncellsize {cell}\nNODATA_value -9999\n"
    with open(path, "w") as f:
        f.write(head)
        for row in z:
            f.write(" ".join(f"{v:.4f}" for v in row) + "\n")
    if sigma is not None and str(path).endswith("_dem.asc"):
        with open(str(path)[:-len("_dem.asc")] + "_sigma.asc", "w") as f:
            f.write(head)
            for row in sg:
                f.write(" ".join(f"{v:.4f}" for v in row) + "\n")
    return write_info(path, info)


def line_row(**kw):
    r = {"line": "1", "status": "ok", "note": "", "window_sharpness": "55.0", "valid_fraction": "0.95",
         "n_peaks": "60", "off_profile_fraction": "0.0", "z_extrapolated_fraction": "0.02", "setup_m": "0.31",
         "water_level_m": "0.3", "water_level_source": "gps"}
    r.update(kw)
    return r


# ---------------------------------------------------------------------

def test_qc():
    print("QC of a line (rule 3), each rule alone")
    check("a clean line passes", ms.qc_line(line_row()) == "")
    cases = [({"status": "dark"}, "status dark"),
             ({"status": "low contrast", "note": "window sharpness 4.1 (< 10)"}, "status low contrast (window sharpness"),
             ({"window_sharpness": "9.9"}, "window sharpness 9.9 < 10"),
             ({"window_sharpness": ""}, "no window sharpness"),
             ({"valid_fraction": "0.7"}, "edge in 0.7 of the rows"),
             ({"n_peaks": "12"}, "12 uprushes"),
             ({"off_profile_fraction": "0.03"}, "off profile 0.03 > 0.02"),
             ({"z_extrapolated_fraction": "0.11"}, "extrapolated 0.11 > 0.1"),
             ({"z_extrapolated_fraction": ""}, "extrapolated"),
             ({"setup_m": ""}, "no setup"),
             ({"water_level_source": "chatham"}, "still water from chatham (not gps)")]
    for kw, want in cases:
        got = ms.qc_line(line_row(**kw))
        check(f"{kw} fails: '{want}'", want in got, got)
    check("--allow-no-sharpness: a line without the OWG frame passes, a blurred one still fails",
          ms.qc_line(line_row(window_sharpness=""), ms.Qc(require_sharpness=False)) == "" and
          "9" in ms.qc_line(line_row(window_sharpness="9"), ms.Qc(require_sharpness=False)))
    check("the limits are exactly inclusive: sharpness 10, off profile 0.02, extrapolated 0.10 pass",
          ms.qc_line(line_row(window_sharpness="10", off_profile_fraction="0.02", z_extrapolated_fraction="0.1")) == "")
    check("another level source allowed: 'level csv'",
          ms.qc_line(line_row(water_level_source="level csv"), ms.Qc(level_sources=("level csv",))) == "")


def test_combine():
    print("the collection: median over the lines and its uncertainty (rules 4, 5)")
    a = {"line": 1, "setup": 0.30, "why": "", "se": 0.02, "wl": 0.3, "wl_src": "gps"}
    b = {"line": 2, "setup": 0.36, "why": "", "se": 0.03, "wl": 0.3, "wl_src": "gps"}
    c = ms.combine_lines([a, b])
    check("two lines: the setup is their median (= mean) 0.33", abs(c["setup"] - 0.33) < 1e-12 and c["n"] == 2
          and c["used"] == "1+2" and c["status"] == "ok", c)
    check("between-line term sd/sqrt(n) = half the difference, 0.030", abs(c["lines"] - 0.03) < 1e-12, c["lines"])
    check("sampling term: the lines' mean (fully correlated, not reduced), 0.025", abs(c["sampling"] - 0.025) < 1e-12)
    check("setup_sigma = sqrt(0.025^2 + 0.030^2)", abs(c["sigma"] - np.hypot(0.025, 0.03)) < 1e-12)
    bad = dict(b, why="window sharpness 4.1 < 10")
    c = ms.combine_lines([a, bad])
    check("a failed line does not count: one line, its setup, sigma = its sampling term, noted",
          c["setup"] == 0.30 and c["n"] == 1 and c["sigma"] == 0.02 and np.isnan(c["lines"])
          and "one line" in c["note"], c)
    check("--min-lines 2 with one passing: failed qc", ms.combine_lines([a, bad], 2)["status"] == "failed qc")
    c = ms.combine_lines([dict(a, why="x"), bad])
    check("no line passing: failed qc, no setup", c["status"] == "failed qc" and np.isnan(c["setup"])
          and c["note"] == "no line passed")
    c = ms.combine_lines([a, dict(b, setup=0.20, line=2), dict(b, setup=0.90, line=3)])
    check("three lines: the median (0.30), not the mean", c["setup"] == 0.30)
    row = ms.collection_row("1792500000.c2.ras.tiff", 1792500000, [a, bad], {"dem": "dem_x", "dem_window": "w",
                                                                              "dem_age_days": 1, "dem_sha256": "ab"})
    check("the row: line_setups and line_qc for every line, the DEM fields, the start epoch",
          row["line_setups"] == "1:0.300;2:0.360" and row["line_qc"] == "1:ok;2:window sharpness 4.1 < 10"
          and row["dem"] == "dem_x" and row["collection_epoch"] == 1792500000 and row["setup_m"] == 0.3
          and list(row) == ms.FIELDS, row)


def _gen(spec, n_long, rng, fs=2.0):
    f = np.fft.rfftfreq(n_long, 1.0 / fs)
    X = np.sqrt(spec(f)) * (rng.normal(size=len(f)) + 1j * rng.normal(size=len(f)))
    X[0] = 0
    return np.fft.irfft(X, n_long)


def test_sampling_se():
    print("the sampling error of a 10-minute mean, on Gaussian records of known spectrum")
    specs = {"flat to zero frequency": (lambda f: ((f > 0) & (f < 0.05)).astype(float), 0.85, 1.15),
             "empty below 0.005 Hz": (lambda f: ((f >= 0.005) & (f < 0.05)).astype(float), 1.5, 6.0),
             "level below 1/T, incident peak": (lambda f: np.where((f > 0) & (f < 0.05), 1 / np.maximum(f, 0.003), 0)
                                                + 30 * np.exp(-0.5 * ((f - 0.1) / 0.02) ** 2), 0.8, 1.15)}
    rng = np.random.default_rng(1)
    n = 1200
    for name, (sp, lo, hi) in specs.items():
        means, ests = [], []
        for _ in range(40):
            x = _gen(sp, n * 20, rng)
            for k in range(20):
                w = x[k * n:(k + 1) * n]
                means.append(w.mean())
                ests.append(ms.mean_sampling_se(w))
        ratio = float(np.sqrt(np.mean(np.square(ests))) / np.std(means))
        print(f"    {name:<34} predicted / true sd of the 10-min mean: {ratio:.2f}")
        check(f"{name}: {lo:g} <= ratio <= {hi:g}", lo <= ratio <= hi, ratio)
    ratios = []
    for _ in range(100):
        x = _gen(specs["flat to zero frequency"][0], n, rng)
        g = x.copy()
        s0 = int(rng.integers(0, n - 40))
        g[s0:s0 + 40] = np.nan
        ratios.append(ms.mean_sampling_se(g) / ms.mean_sampling_se(x))
    check("a gappy record (20 s bridged): the same value on median (within 5%; one record within ~20%)",
          abs(float(np.median(ratios)) - 1) < 0.05, float(np.median(ratios)))
    check("under half the samples: no value", np.isnan(ms.mean_sampling_se(np.where(np.arange(n) < 700, np.nan, x))))


def _tiny_asc(path, info=MODEL_INFO):
    Path(path).write_text("ncols 2\nnrows 2\nxllcorner 0\nyllcorner 0\ncellsize 2\nNODATA_value -9999\n1 2\n3 4\n")
    return write_info(path, info)


def test_prior_dem():
    print("the strictly-earlier, model-built DEM (rule 2)")
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp) / "dems"
        d.mkdir()
        for n in ("dem_2026-10-17_7d_dem.asc", "dem_2026-10-19_7d_dem.asc", "dem_2026-10-19_14d_dem.asc",
                  "dem_2026-10-20_7d_dem.asc", "dem_2026-10-22_7d_dem.asc", "dem_2026-10-19_7d_spread.asc",
                  "other_dem.asc"):
            _tiny_asc(d / n)
        pd_ = ms.PriorDem(d)
        e = 1792500000                          # 2026-10-20 12:40 UTC
        f, info = pd_.pick(e)
        check("the stack of 20 Oct: the DEM ending 19 Oct (7 days, not the 14-day one; not 20 or 22 Oct)",
              f is not None and f.name == "dem_2026-10-19_7d_dem.asc", info)
        check("it records the window, the age, the sha256 and the DEM's setup share",
              info.get("dem_window") == "2026-10-13..2026-10-19" and info.get("dem_age_days") == 1 and
              info.get("dem_sha256") == hashlib.sha256((d / "dem_2026-10-19_7d_dem.asc").read_bytes()).hexdigest()[:12]
              and info.get("dem_setup_share") == 1.0, info)
        f, info = pd_.pick(e + 86400)
        check("21 Oct: the one ending 20 Oct (which holds 20 Oct, not 21)", f.name == "dem_2026-10-20_7d_dem.asc")
        f, info = pd_.pick(e + 4 * 86400)
        check("24 Oct: the one ending 22 Oct (age 2)", f.name == "dem_2026-10-22_7d_dem.asc" and info["dem_age_days"] == 2)
        f, info = pd_.pick(e + 5 * 86400)
        check("25 Oct: still the one ending 22 Oct, at the age limit (3 days)",
              f is not None and f.name == "dem_2026-10-22_7d_dem.asc" and info["dem_age_days"] == 3, info)
        f, why = pd_.pick(e + 7 * 86400)
        check("27 Oct: none within 3 days: no DEM, and why", f is None and "1-3 day(s) before 2026-10-27" in why, why)
        f, why = pd_.pick(e - 4 * 86400)
        check("16 Oct: none earlier", f is None)
        bad = []
        for t in np.arange(e - 6 * 86400, e + 8 * 86400, 3 * 3600):
            f, info = pd_.pick(t)
            if f is None:
                continue
            day = ms.utc_day(t)
            first, last = [datetime.strptime(x, "%Y-%m-%d").date() for x in info["dem_window"].split("..")]
            if not last < day or first <= day <= last:
                bad.append((day, f.name))
        check("over two weeks of stack times every 3 h: no DEM's window ever holds the stack's day", not bad, bad)
        check("for_epoch returns the DEM estimate_stack reads (its name), none where there is none",
              pd_.for_epoch(e).name == "dem_2026-10-19_7d_dem.asc" and pd_.for_epoch(e + 7 * 86400) is None)
        one = ms.PriorDem(d / "other_dem.asc")
        check("a single DEM file without a series name: refused (its window is not known)", one.refused)
        one = ms.PriorDem(d / "dem_2026-10-19_7d_dem.asc")
        check("a single dated DEM file: accepted, the same rule", not one.refused and one.pick(e)[0] is not None
              and one.pick(e - 86400)[0] is None)
        args = ms.build_parser().parse_args(["--dem", str(d / "other_dem.asc"), "--output", str(Path(tmp) / "o.csv"),
                                             "--ras-dir", str(Path(tmp))])
        try:
            quiet(ms.run, args)
            refused = False
        except SystemExit as exc:
            refused = "strictly-earlier rule cannot be applied" in str(exc)
        check("the run refuses it before any stack is read", refused)

        # rule 2b: only a DEM built with the model's setup, by its _info.json
        cases = [("built with measured setups", {"setup_share": 1.0, "setup_source": {"measured": 12, "model": 3}},
                  "built with 12 measured frame(s)"),
                 ("frames of unknown source", {"setup_share": 1.0, "setup_source": {"model": 5, "unrecorded": 2}},
                  "built with 2 unrecorded frame(s)"),
                 ("no setup on its lines", {"setup_share": 0.0}, "no setup on its lines"),
                 ("no setup_share at all", {"frames": {"c2": 3}}, "no setup on its lines"),
                 ("no _info.json", None, "no _info.json"),
                 ("an unreadable _info.json", "{not json", "unreadable")]
        for what, info, want in cases:
            d2 = Path(tmp) / ("b_" + what.replace(" ", "_"))
            d2.mkdir()
            _tiny_asc(d2 / "dem_2026-10-18_7d_dem.asc")
            _tiny_asc(d2 / "dem_2026-10-19_7d_dem.asc", info if isinstance(info, dict) else None)
            if isinstance(info, str):
                (d2 / "dem_2026-10-19_7d_info.json").write_text(info)
            p2 = ms.PriorDem(d2)
            f, inf = p2.pick(e)
            check(f"a DEM {what}: refused, the model-built one a day older taken (age 2), the refusal listed",
                  f is not None and f.name == "dem_2026-10-18_7d_dem.asc" and inf["dem_age_days"] == 2
                  and want in p2.not_model_built.get("dem_2026-10-19_7d_dem.asc", ""), (f, p2.not_model_built))
        d3 = Path(tmp) / "accepted"
        d3.mkdir()
        _tiny_asc(d3 / "dem_2026-10-19_7d_dem.asc", {"setup_share": 0.97, "setup_source": {"model": 40, "none": 2}})
        f, inf = ms.PriorDem(d3).pick(e)
        check("auto's DEM without a measured frame (model and none only): accepted, its setup share recorded",
              f is not None and inf["dem_setup_share"] == 0.97, inf)
        d4 = Path(tmp) / "all_refused"
        d4.mkdir()
        _tiny_asc(d4 / "dem_2026-10-19_7d_dem.asc", {"setup_share": 1.0, "setup_source": {"measured": 4}})
        p4 = ms.PriorDem(d4)
        f, why = p4.pick(e)
        check("only a DEM built with measured setups: no DEM, and why (named)", f is None and p4.for_epoch(e) is None
              and "no model-built DEM whose window ends 1-3 day(s) before 2026-10-20" in why
              and "dem_2026-10-19_7d_dem.asc: built with 4 measured frame(s)" in why, why)
        base = Path(tmp) / "station"
        (base / "archive" / "dems").mkdir(parents=True)
        first = ms.default_dem(base)
        (base / "archive" / "dems_model").mkdir()
        check("--dem's default: archive/dems_model when it exists (the cron's model-only series), else archive/dems",
              first == base / "archive" / "dems" and ms.default_dem(base) == base / "archive" / "dems_model"
              and ms.build_parser().parse_args([]).dem is None)


def test_csv():
    print("the output file: appended, no duplicates, --reprocess, another header refused; the extractor's reading")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "setup_measured_c2.csv"

        def row(ep, setup, status="ok"):
            r = {k: "" for k in ms.FIELDS}
            r.update(collection_epoch=ep, filename=f"{ep}.c2.ras.tiff", status=status, setup_m=setup)
            return r
        n1 = ms.write_rows(out, [row(1792500000, 0.3), row(1792501800, 0.4)])
        n2 = ms.write_rows(out, [row(1792501800, 0.9), row(1792503600, 0.5), row(1792503600, 0.6)])
        rows = list(csv.DictReader(open(out)))
        check("appended: 2 then 1 new (a collection already there and a duplicate in the batch are not)",
              n1 == 2 and n2 == 1 and [r["setup_m"] for r in rows] == ["0.3", "0.4", "0.5"], [r["setup_m"] for r in rows])
        ms.write_rows(out, [row(1792501800, 0.45)], reprocess=True)
        rows = list(csv.DictReader(open(out)))
        check("--reprocess replaces that collection only, the file in time order",
              [r["setup_m"] for r in rows] == ["0.3", "0.45", "0.5"])
        (Path(tmp) / "other.csv").write_text("a,b\n1,2\n")
        try:
            ms.write_rows(Path(tmp) / "other.csv", [row(1, 0.1)])
            refused = False
        except SystemExit:
            refused = True
        check("a file with other columns is refused", refused and (Path(tmp) / "other.csv").read_text() == "a,b\n1,2\n")
        ms.write_rows(out, [row(1792505400, "", "failed qc"), row(1792507200, 0.7)])
        t = ms.MeasuredSetups(out)
        check("the extractor's table: the QC-passed collections only", t.rows == 5 and t.ok == 4)
        check("a frame within 300 s of a collection's start gets its setup (0 / 250 / 300 s), not at 301 s or "
              "beside a failed one",
              t.at(1792500000) == 0.3 and t.at(1792500250) == 0.3 and t.at(1792501800 - 300) == 0.45
              and t.at(1792500000 + 1800 + 301 + 600) is None and t.at(1792505400) is None,
              (t.at(1792500250), t.at(1792505400)))
        try:
            ms.MeasuredSetups(Path(tmp) / "other.csv")
            ok = False
        except ValueError:
            ok = True
        check("a file that is not a measured-setup file: ValueError", ok)


# ---------------------------------------------------------------------
# The synthetic scene (test_timestack_wave_height.py's)
# ---------------------------------------------------------------------

def scene(tmp, beta=0.10, h0=1.5, t0=10.0, sharp=55.0):
    import test_timestack_wave_height as T
    sc = T.make_scene(tmp, h0, t0, beta)
    sc["epoch"] = T.EPOCH
    owg = Path(tmp) / "owg.csv"
    owg.write_text("filename,epoch,time_utc,brightness,sharpness,saturation,status\n"
                   f"{T.EPOCH}.c2.bright.jpg,{T.EPOCH},x,160,{sharp},0.0,ok\n")
    sc["owg"] = owg
    return sc


def scene_args(sc, dems, out, extra=()):
    return ms.build_parser().parse_args(["--ras-dir", str(sc["ras"]), "--pix", str(sc["pix"]), "--dem", str(dems),
                                         "--level-csv", str(sc["level"]), "--owg-csv", str(sc["owg"]),
                                         "--lines", "1", "--output", str(out), "--wait-days", "0"] + list(extra))


def measure(sc, dems, extra=(), now=None, lines=(1,)):
    args = scene_args(sc, dems, Path(sc["ras"]).parent / "unused.csv", extra)
    ctx, _ = quiet(ms.make_context, args, ms.PriorDem(args.dem, args.max_dem_age_days))
    r, _ = quiet(ms.measure_stack, Path(sc["ras"]) / sc["name"], ctx, list(lines), now)
    return r


def two_line_stack(sc, tmp):
    """The scene's stack and .pix with the line twice (lines 1 and 2: the same swash)."""
    import cv2
    img = cv2.imread(str(Path(sc["ras"]) / sc["name"]), cv2.IMREAD_UNCHANGED)
    ras2 = Path(tmp) / "ras2"
    ras2.mkdir()
    cv2.imwrite(str(ras2 / sc["name"]), np.hstack([img, img]))
    pix = np.loadtxt(sc["pix"], ndmin=2)
    pix2 = Path(tmp) / "two.pix"
    np.savetxt(pix2, np.vstack([pix, pix]), fmt="%d")
    return dict(sc, ras=ras2, pix=pix2)


def test_scene():
    print("end to end on a synthetic c2 scene: the setup on the DEM of the day before")
    import runup_from_timestack as rft
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = scene(tmp)
        truth = float(sc["z_rel"].mean())
        dems = tmp / "dems"
        dems.mkdir()
        write_plane_dem(dems / "dem_2026-10-19_7d_dem.asc", 0.10, sigma=0.05)
        r = measure(sc, dems)
        ok = isinstance(r, dict) and r["status"] == "ok"
        check("status ok, on the DEM ending the day before (age 1), its window, sha256 and setup share recorded",
              ok and r["dem"] == "dem_2026-10-19_7d_dem.asc" and r["dem_window"] == "2026-10-13..2026-10-19"
              and r["dem_age_days"] == 1 and len(r["dem_sha256"]) == 12 and r["dem_setup_share"] == 1.0, r)
        if not ok:
            return
        setup = float(r["setup_m"])
        print(f"    truth (the swash's own mean) {truth:.3f} m, measured {setup:.3f} +/- {r['setup_sigma_m']} m "
              f"(sampling {r['sigma_sampling_m']}); profile term {r['sigma_profile_m']} m (gain {r['profile_gain']})")
        check("the setup within the detector's bias of the truth (+0.00 to +0.12 m: the edge trails a fast tip)",
              0.0 <= setup - truth <= 0.12, setup - truth)
        check("its sampling error is finite and small (0.005-0.08 m)", 0.005 <= float(r["sigma_sampling_m"]) <= 0.08,
              r["sigma_sampling_m"])
        check("still water from the level file, as the QC allows for --level-csv",
              r["water_level_source"] == "level csv" and r["line_qc"] == "1:ok")
        check("the diagnostics: the window's sharpness (55) and the line's edge contrast, logged",
              r["window_sharpness"] == 55.0 and r["line_edge_contrast"].startswith("1:")
              and np.isfinite(float(r["line_edge_contrast"][2:])), (r["window_sharpness"], r["line_edge_contrast"]))
        check("the profile term = the gain x the DEM's sigma grid (0.05 m)",
              abs(float(r["sigma_profile_m"]) - 0.05 * float(r["profile_gain"])) <= 0.0015, r["sigma_profile_m"])

        # the series hook: estimate_stack's row is the same with it, and the series is the one averaged
        args = scene_args(sc, dems, tmp / "x.csv")
        ctx, _ = quiet(ms.make_context, args, ms.PriorDem(dems))
        rft._PIX_FILES.clear()
        with_hook, _ = quiet(thw.estimate_stack, Path(sc["ras"]) / sc["name"], 1, ctx)
        series = ctx.series_out.get((sc["name"], 1))
        del ctx.series_out
        rft._PIX_FILES.clear()
        without, _ = quiet(thw.estimate_stack, Path(sc["ras"]) / sc["name"], 1, ctx)
        check("timestack_wave_height.estimate_stack: the same row with and without the series hook",
              with_hook == without)
        check("the series captured is the one setup_m is the mean of",
              series is not None and abs(float(np.nanmean(series["eta"])) - float(with_hook["setup_m"])) <= 0.0005)

        # rule 1b: the edge-contrast label is off, the contrast only logged -- a sharp window with a real stack's
        # contrast (0.00-0.02 on every day of Oct 2026) passes
        check("rule 1b: the context's contrast limit is off (-inf)", ctx.min_contrast == float("-inf"))
        real = thw.edge_contrast
        thw.edge_contrast = lambda gray, edge, *a, **k: 0.01
        try:
            rft._PIX_FILES.clear()
            lowc = measure(sc, dems)
        finally:
            thw.edge_contrast = real
        check("a sharp window with the real stacks' edge contrast (0.01): passes, the contrast logged",
              isinstance(lowc, dict) and lowc["status"] == "ok" and lowc["line_edge_contrast"] == "1:0.010",
              lowc if not isinstance(lowc, dict) else (lowc["status"], lowc["line_qc"], lowc["line_edge_contrast"]))

        # two lines (the same swash twice): the median, n 2, the between-line term 0
        sc2 = two_line_stack(sc, tmp)
        rft._PIX_FILES.clear()
        r2 = measure(sc2, dems, lines=(1, 2))
        rft._PIX_FILES.clear()
        check("two lines: both pass, lines_used 1+2, the same setup, between-line term 0",
              isinstance(r2, dict) and r2["n_lines"] == 2 and r2["lines_used"] == "1+2"
              and r2["setup_m"] == r["setup_m"] and float(r2["sigma_lines_m"]) == 0.0, r2)

        # the inheritance: a DEM off the beach by d moves the setup by d / (1 - beta / tan(delta))
        shifts = {}
        for off in (0.10, -0.10):
            d = tmp / f"dems{off:+.2f}"
            d.mkdir()
            write_plane_dem(d / "dem_2026-10-19_7d_dem.asc", 0.10, off)
            ro = measure(sc, d)
            shifts[off] = float(ro["setup_m"]) - setup if isinstance(ro, dict) and ro["setup_m"] != "" else np.nan
        print(f"    a DEM 0.10 m above the beach: setup {shifts[0.10]:+.3f} m; 0.10 m below: {shifts[-0.10]:+.3f} m "
              f"(tan(delta) = beta / (1 - 0.10 / shift) = {0.10 / (1 - 0.10 / shifts[0.10]):.2f})")
        check("INHERITED: a DEM 0.10 m too high reads the setup 0.12-0.20 m higher (it is read off the DEM's level)",
              0.12 <= shifts[0.10] <= 0.20, shifts)
        check("... and 0.10 m too low as much lower (beach change since the DEM is read as setup)",
              abs(shifts[-0.10] + shifts[0.10]) <= 0.01, shifts)
        check("profile_gain, from the ray's geometry, is that response: 1 / (1 - beta / tan(delta)) = shift / 0.10 "
              "within 3%", abs(float(r["profile_gain"]) / (shifts[0.10] / 0.10) - 1) <= 0.03,
              (r["profile_gain"], shifts[0.10] / 0.10))

        # rule 2b: a DEM built with measured setups is never stood on
        mdir = tmp / "measured_built"
        mdir.mkdir()
        write_plane_dem(mdir / "dem_2026-10-19_7d_dem.asc", 0.10,
                        info={"setup_share": 1.0, "setup_source": {"measured": 30, "model": 4}})
        ro = measure(sc, mdir)
        check("only a DEM built with measured setups: not used ('status no slope', the reason in the note)",
              isinstance(ro, dict) and ro["status"] == "failed qc" and "status no slope" in ro["line_qc"]
              and "built with 30 measured frame(s)" in ro["note"], ro)

        # the DEM of the stack's own day only: none earlier -> no measurement (old) / wait (young)
        own = tmp / "own"
        own.mkdir()
        write_plane_dem(own / "dem_2026-10-20_7d_dem.asc", 0.10)
        ro = measure(sc, own)
        check("only the DEM of the stack's own day: not used; 'status no slope', the reason in the note",
              isinstance(ro, dict) and ro["status"] == "failed qc" and "status no slope" in ro["line_qc"]
              and "no DEM whose window ends 1-3 day(s) before 2026-10-20" in ro["note"], ro)
        ro = measure(sc, own, ["--wait-days", "3"], now=sc["epoch"] + 3600)
        check("... a young stack waits for the DEM of the day before", ro == "wait", ro)
        ro = measure(sc, dems, ["--wait-days", "3"], now=sc["epoch"] + 3600)
        check("... and with it, a young stack is measured at once", isinstance(ro, dict) and ro["status"] == "ok")

        # the window and still-water rules
        sc["owg"].write_text("filename,epoch,time_utc,brightness,sharpness,saturation,status\n"
                             f"{sc['epoch']}.c2.bright.jpg,{sc['epoch']},x,160,4.1,0.0,blurred\n")
        ro = measure(sc, dems)
        check("a blurred window (sharpness 4.1): failed qc", ro["status"] == "failed qc"
              and "window sharpness 4.1" in ro["line_qc"], ro["line_qc"])
        sc["owg"].unlink()
        ro = measure(sc, dems)
        check("no OWG frame of the collection: failed qc (the window unknown)", ro["status"] == "failed qc"
              and "no window sharpness" in ro["line_qc"], ro["line_qc"])
        ro = measure(sc, dems, ["--allow-no-sharpness"])
        check("... passed with --allow-no-sharpness", ro["status"] == "ok")
        ro = measure(sc, dems, ["--allow-no-sharpness", "--level-sources", "gps"])
        check("a still water other than the one the frames get: failed qc", ro["status"] == "failed qc"
              and "still water from level csv (not gps)" in ro["line_qc"], ro["line_qc"])
        ro = measure(sc, dems, ["--allow-no-sharpness", "--level-sources", "gps", "--wait-days", "3"],
                     now=sc["epoch"] + 3600)
        check("rule 6: a YOUNG stack on a non-GPS level where the GPS is the allowed source waits (a run whose "
              "spline could not be read does not write it off)", ro == "wait", ro)

        # the run over the archive: appended once, --reprocess
        args = scene_args(sc, dems, tmp / "out.csv", ["--allow-no-sharpness"])
        quiet(ms.run, args)
        quiet(ms.run, args)
        rows = list(csv.DictReader(open(tmp / "out.csv")))
        args.reprocess = True
        quiet(ms.run, args)
        rows2 = list(csv.DictReader(open(tmp / "out.csv")))
        head = next(csv.reader(open(tmp / "out.csv")))
        check("the run: one row per collection, a second run adds none, --reprocess keeps one; its columns FIELDS",
              len(rows) == 1 and len(rows2) == 1 and head == ms.FIELDS and rows[0]["status"] == "ok", (len(rows), len(rows2)))

        # the audit (criterion 4): the lines' noise table, the failures, and the sample drawn with the edge on
        adir = tmp / "audit"
        rft._PIX_FILES.clear()
        rc, txt = quiet(ms.main, ["--ras-dir", str(sc["ras"]), "--pix", str(sc["pix"]), "--dem", str(dems),
                                  "--level-csv", str(sc["level"]), "--owg-csv", str(tmp / "none.csv"),
                                  "--allow-no-sharpness", "--lines", "1", "--output", str(tmp / "out.csv"),
                                  "--audit", "5", "--audit-dir", str(adir)])
        idx = list(csv.DictReader(open(adir / "index.csv"))) if (adir / "index.csv").exists() else []
        check("--audit: the noise table and a passed collection drawn with the edge (a PNG, index.csv); the "
              "measurement file untouched",
              rc == 0 and "|line 1 - line 2|" in txt and len(idx) == 1 and (adir / idx[0]["png"]).exists()
              and idx[0]["setup_now_m"] == idx[0]["setup_m"]
              and list(csv.DictReader(open(tmp / "out.csv"))) == rows2, (rc, txt[-400:], idx))
        check("--audit, rule 2b afterwards: the DEM the passed collection stood on is model-built (0 not)",
              "stood on 1 DEM(s)" in txt and "not model-built now: 0 (0 collection(s); must be 0)" in txt, txt[-600:])
        write_info(dems / "dem_2026-10-19_7d_dem.asc", {"setup_share": 1.0, "setup_source": {"measured": 9}})
        rc, txt = quiet(ms.main, ["--dem", str(dems), "--output", str(tmp / "out.csv"), "--audit", "0",
                                  "--audit-dir", str(adir)])
        write_info(dems / "dem_2026-10-19_7d_dem.asc")
        check("... and a DEM that now says it was built with measured setups is named (1 collection)",
              rc == 0 and "not model-built now: 1 (1 collection(s); must be 0)" in txt
              and "dem_2026-10-19_7d_dem.asc: built with 9 measured frame(s)" in txt, txt[-600:])
        noise = ms.line_noise([{"status": "ok", "lines_used": "1+2", "line_setups": "1:0.300;2:0.360",
                                "window_sharpness": "60", "line_edge_contrast": "1:0.010;2:0.020"},
                               {"status": "ok", "lines_used": "1+2", "line_setups": "1:0.300;2:0.200",
                                "window_sharpness": "15", "line_edge_contrast": "1:0.0"},
                               {"status": "ok", "lines_used": "2", "line_setups": "1:0.900;2:0.250",
                                "window_sharpness": "70", "line_edge_contrast": ""}])
        check("the noise table: |line 1 - line 2| of the collections where both passed, by the window's sharpness",
              noise["clear (>= 25)"]["n"] == 1 and abs(noise["clear (>= 25)"]["median"] - 0.06) < 1e-9
              and noise["low (10-25)"]["n"] == 1 and abs(noise["low (10-25)"]["median"] - 0.10) < 1e-9
              and abs(noise["clear (>= 25)"]["contrast_median"] - 0.015) < 1e-9, noise)


CHAIN_DAYS, CHAIN_WINDOW, STORM_DAY, STORM_DROP = 60, 7, 20, -0.30


def chain_series(F, base, use_files, d0=datetime(2026, 11, 1, tzinfo=timezone.utc), drop=STORM_DROP):
    """
    The cron's 7-day window DEMs at the c2 line over CHAIN_DAYS days, each day's frames given the setup measured
    on the DEM measured_setup.py picks (F: a DEM error at the line -> that day's line error), the beach dropping
    by `drop` from STORM_DAY. use_files None: the loop rule 2b forbids (each day on the previous day's DEM
    built the same way). Else the station's folders under `base`: archive/dems (the cron's own series, built
    with measured setups from day 1, its _info.json counting them) and, if use_files == 'model series', also
    archive/dems_model (the model-only series); each day ms.PriorDem(ms.default_dem(base)) picks.
    -> (DEM level at the line by window end, true beach by day, days with a measured setup, names picked).
    """
    T = {k: (0.0 if k < STORM_DAY else drop) for k in range(-CHAIN_WINDOW - 5, CHAIN_DAYS + 1)}
    model = {D: float(np.mean([T[k] for k in range(D - CHAIN_WINDOW + 1, D + 1)])) for D in range(-4, CHAIN_DAYS + 1)}
    level = {D: model[D] for D in range(-4, 1)}                 # before the switch: model-built
    measured_day, picked = {}, []
    dems, dmod = Path(base) / "archive" / "dems", Path(base) / "archive" / "dems_model"
    names = {}
    if use_files is not None:
        dems.mkdir(parents=True, exist_ok=True)
        if use_files == "model series":
            dmod.mkdir(parents=True, exist_ok=True)
        for D in range(-4, 1):
            for folder in ([dems, dmod] if use_files == "model series" else [dems]):
                f = _tiny_asc(folder / f"dem_{(d0 + timedelta(days=D)).date()}_{CHAIN_WINDOW}d_dem.asc")
                names[(folder, f.name)] = model[D]
    point = {}
    for D in range(1, CHAIN_DAYS + 1):
        if use_files is None:
            point[D] = T[D] + F(level[D - 1] - T[D])
            measured_day[D] = True
        else:
            pd_ = ms.PriorDem(ms.default_dem(base))
            f, _info = pd_.pick((d0 + timedelta(days=D, hours=12)).timestamp())
            if f is None:
                point[D], measured_day[D] = T[D], False      # no DEM: the model's setup (auto), its error 0 here
            else:
                picked.append((f.parent.name, f.name))
                point[D], measured_day[D] = T[D] + F(names[(f.parent, f.name)] - T[D]), True
        win = [k for k in range(D - CHAIN_WINDOW + 1, D + 1)]
        level[D] = float(np.mean([point[k] if k >= 1 else T[k] for k in win]))
        if use_files is not None:
            n_meas = sum(1 for k in win if k >= 1 and measured_day[k])
            src = {"measured": n_meas, "model": CHAIN_WINDOW - n_meas} if n_meas else {"model": CHAIN_WINDOW}
            day = (d0 + timedelta(days=D)).date()
            f = _tiny_asc(dems / f"dem_{day}_{CHAIN_WINDOW}d_dem.asc", {"setup_share": 1.0, "setup_source": src})
            names[(dems, f.name)] = level[D]
            if use_files == "model series":
                f = _tiny_asc(dmod / f"dem_{day}_{CHAIN_WINDOW}d_dem.asc")
                names[(dmod, f.name)] = model[D]
    return level, T, measured_day, picked


def test_chain():
    print("chained generations (rule 2b): measured setups on DEMs built with measured setups, and on model-built ones")
    import georectify
    import runup_from_timestack as rft
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = scene(tmp)
        truth = float(sc["z_rel"].mean())
        s_of = {}
        for off in (-0.10, 0.0, 0.10):
            d = tmp / f"d{off:+.2f}"
            d.mkdir()
            write_plane_dem(d / "dem_2026-10-19_7d_dem.asc", 0.10, off)
            rft._PIX_FILES.clear()
            r = measure(sc, d)
            s_of[off] = float(r["setup_m"]) if isinstance(r, dict) and r["setup_m"] != "" else np.nan
        g_s = (s_of[0.10] - s_of[-0.10]) / 0.20
        # the timex line placed exactly at the mean swash position on the TRUE beach (the measured setup's best
        # case), given still water + the measured setup and re-projected on its ray onto that plane
        io_ = georectify.load_intrinsics(HERE / "calibration" / "CACO05_c2_20240801_IO.yaml")
        eo_ = georectify.load_extrinsics(HERE / "calibration" / "CACO05_c2_20251113_EO-CV.yaml")
        pix = np.loadtxt(sc["pix"], ndmin=2)
        zp = np.asarray(sc["z_pix"], float)
        swl = 0.3
        fin = np.flatnonzero(np.isfinite(zp))
        order = fin[np.argsort(zp[fin])]
        k = float(np.interp(swl + truth, zp[order], order.astype(float)))
        i0 = int(np.floor(k))
        u = pix[i0, 0] + (k - i0) * (pix[i0 + 1, 0] - pix[i0, 0])
        v = pix[i0, 1] + (k - i0) * (pix[i0 + 1, 1] - pix[i0, 1])

        def line_error(setup):
            zz = swl + setup
            E, N = georectify.pixel_to_ground(np.array([u]), np.array([v]), zz, io_, eo_)
            return zz - (-0.10 * float(thw.seaward_distance(np.asarray(E, float), np.asarray(N, float))[0]))

        def F(e):
            return line_error(s_of[0.0] + g_s * e)
        gain = (F(0.10) - F(-0.10)) / 0.20
        print(f"    setup read on a DEM off by e: {s_of[0.0]:.3f} + {g_s:.3f} e (truth {truth:.3f}); the line given it: "
              f"error {F(0.0):+.4f} + {gain:.3f} e (true setup: {line_error(truth):+.4f})")
        check("the line given the TRUE setup lands on the beach (within 0.01 m): the geometry is right",
              abs(line_error(truth)) <= 0.01, line_error(truth))
        check("UNIT GAIN: the setup reads e x 1.4-1.7, the re-projected line turns it back into e (0.95-1.05)",
              1.4 <= g_s <= 1.7 and 0.95 <= gain <= 1.05, (g_s, gain))
        off0 = F(0.0)

        calm, Tc, _m, _p = chain_series(F, tmp / "unused", None, drop=0.0)
        loop, T, _m, _p = chain_series(F, tmp / "unused", None)
        print("    the loop rule 2b forbids (each day on the last DEM built from measured setups), no storm: error at "
              "the line " + ", ".join(f"day {D} {calm[D] - Tc[D]:+.3f}" for D in (7, 14, 30, 60))
              + f"; with the erosion on day {STORM_DAY}: change day 15 -> 35 {loop[35] - loop[15]:+.3f} m (true "
              f"{STORM_DROP:+.2f})")
        check("... drifts: more than 0.5 m off after 60 days with no storm, and the 0.30 m erosion is not seen "
              "(change > -0.10 m)", calm[60] - Tc[60] > 0.5 and loop[35] - loop[15] > -0.10,
              (calm[60] - Tc[60], loop[35] - loop[15]))

        dep, T, meas, picked = chain_series(F, tmp / "station_a", "model series")
        print(f"    deployed (archive/dems_model, the default --dem): error {dep[15] - T[15]:+.3f} (day 15), "
              f"{dep[60] - T[60]:+.3f} (day 60); change 15 -> 35 {dep[35] - dep[15]:+.3f} m; measured on "
              f"{sum(meas.values())} of {CHAIN_DAYS} days")
        check("deployed with the model-only series: every pick a model-built DEM of dems_model, every day measured",
              len(picked) == CHAIN_DAYS and all(p[0] == "dems_model" for p in picked), picked[:3])
        check("... the error stays one generation's (within 5 mm of it on day 15 and day 60: no drift)",
              abs(dep[15] - T[15] - off0) <= 0.005 and abs(dep[60] - T[60] - off0) <= 0.005,
              (dep[15] - T[15], dep[60] - T[60], off0))
        check("... and the erosion is kept (change 15 -> 35 within 0.02 m of -0.30)",
              abs(dep[35] - dep[15] - STORM_DROP) <= 0.02, dep[35] - dep[15])

        own, T, meas, picked = chain_series(F, tmp / "station_b", "own series only")
        worst = max(abs(own[D] - T[D]) for D in range(1, CHAIN_DAYS + 1) if not STORM_DAY <= D < STORM_DAY + 10)
        print(f"    no model-only series (archive/dems alone): measured on {sum(meas.values())} of {CHAIN_DAYS} days, "
              f"worst error away from the storm {worst:.3f} m; change 15 -> 35 {own[35] - own[15]:+.3f} m")
        check("without the model-only series the _info.json rule alone still stops the drift: no pick of a DEM built "
              "with measured setups, the error never above one generation's (+0.005 m) away from the storm",
              worst <= abs(off0) + 0.005 and all(p[0] == "dems" for p in picked)
              and abs(own[35] - own[15] - STORM_DROP) <= 0.05, (worst, own[35] - own[15]))
        check("... but it starves the measurement (fewer than half the days): hence the cron's model-only series",
              sum(meas.values()) < CHAIN_DAYS / 2, sum(meas.values()))


# ---------------------------------------------------------------------
# extract_elevation_contours.py
# ---------------------------------------------------------------------

T0 = int(datetime(2026, 10, 9, 12, 0, tzinfo=timezone.utc).timestamp())
WAVES = [(0.42, 5.3), (0.87, 6.1), (1.13, 7.7), (1.38, None), (1.66, 9.4), (1.91, 10.6)]   # frame 3: no Tp (an outage)


def extract_fixture(d):
    """6 collections 2 h apart, a c1 and a c2 timex line each, the tide and the waves (frame 3 without Tp): every
    number an exact decimal, so the extractor's output is the same on any platform."""
    d = Path(d)
    folder = d / "processed"
    folder.mkdir(parents=True, exist_ok=True)
    for i in range(len(WAVES)):
        ep = T0 + i * 7200
        t = datetime.fromtimestamp(ep, tz=timezone.utc)
        for cam, du in (("c1", 0), ("c2", 1)):
            lines = ["Column,Row,Confidence,Row_Precise,Peak_Sharpness,Has_Signal"]
            for c in range(300):
                row = 500 + 10 * i + du + (c * 7 % 13) / 4.0
                lines.append(f"{c},{int(row)},0.95,{row:.3f},0.5,1")
            (folder / f"{ep}.{t:%a.%b.%d_%H_%M_%S}.GMT.{t:%Y}.CACO05.{cam}.timex.csv").write_text("\n".join(lines) + "\n")
    tide = ["time,level"] + [f"{datetime.fromtimestamp(T0 + k * 600, tz=timezone.utc):%Y-%m-%d %H:%M:%S},"
                             f"{0.2 + 0.002 * k:.4f}" for k in range(-6, 80)]
    (d / "tide.csv").write_text("\n".join(tide) + "\n")
    (d / "waves_marconi.csv").write_text("time_utc,epoch,hs_best,hs_sigma,tp_s\n" + "".join(
        f"x,{T0 + i * 7200},{h:.2f},0.1,{'' if p is None else f'{p:.1f}'}\n" for i, (h, p) in enumerate(WAVES)))
    dems = d / "dems"
    dems.mkdir(exist_ok=True)
    with open(dems / "dem_2026-10-08_7d_dem.asc", "w") as f:        # a 0.08 plane for the slope-aware model
        f.write("ncols 60\nnrows 40\nxllcorner 420000\nyllcorner 4638000\ncellsize 2\nNODATA_value -9999\n")
        for _ in range(40):
            f.write(" ".join(f"{2.0 - 0.16 * k:.4f}" for k in range(60)) + "\n")
    # the measured setups: frame 0 at its start, frame 1 250 s off, frame 2 400 s off (outside), frame 3 (no Tp)
    # measured, frame 4 failed QC, frame 5 none
    meas = d / "setup_measured_c2.csv"
    rows = []
    for i, off, setup, st in ((0, 0, 0.25, "ok"), (1, 250, 0.3, "ok"), (2, 400, 0.9, "ok"), (3, -120, 0.41, "ok"),
                              (4, 0, "", "failed qc")):
        r = {k: "" for k in ms.FIELDS}
        r.update(collection_epoch=T0 + i * 7200 + off, filename=f"{T0 + i * 7200 + off}.c2.ras.tiff", status=st,
                 setup_m=setup)
        rows.append(r)
    ms.write_rows(meas, rows)
    return {"meas": meas, "dems": dems}


def read_rows(p):
    with open(p, newline="") as f:
        rd = csv.DictReader(f)
        return list(rd.fieldnames or []), list(rd)


def per_frame(rows):
    out = {}
    for r in rows:
        out.setdefault((int(r["capture_epoch"]), r["camera"]), r)
    return out


def old_extractor(dst):
    """BASE_COMMIT's extractor in dst (None when the history is not here)."""
    try:
        txt = subprocess.run(["git", "-C", str(HERE), "show", f"{BASE_COMMIT}:extract_elevation_contours.py"],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None
    p = Path(dst) / "old" / "extract_elevation_contours.py"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(txt)
    return p


def test_extractor():
    print("extract_elevation_contours.py --setup-source: the default as before, byte for byte; measured and auto")
    new = HERE / "extract_elevation_contours.py"
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        fx = extract_fixture(d)
        base = [d / "tide.csv", "--processed-dir", d / "processed"]
        W = ["--waves", d / "waves_marconi.csv"]
        C = W + ["--setup-coef", "0.037"]
        S = C + ["--setup-model", "stockdon", "--setup-slope-dem", fx["dems"]]
        env = dict(os.environ, PYTHONPATH=str(HERE))
        old = old_extractor(d)
        sets = [("no waves", [], []), ("--waves", W, W), ("--setup-coef", C, C),
                ("--setup-model stockdon", S, S),
                ("--setup-source model", C, C + ["--setup-source", "model"]),
                ("--setup-source model with stockdon", S, S + ["--setup-source", "model"]),
                ("--setup-measured without --setup-source", C, C + ["--setup-measured", fx["meas"]])]
        out = d / "c.csv"
        if old is None:
            print(f"    NOTE: commit {BASE_COMMIT} is not in this repository's history: the pinned sha256 below "
                  f"is the check of the default")
        for name, old_args, new_args in sets:
            if old is None:
                break
            rc0, log0 = run([old] + base + old_args + ["--output", out], d, env)
            b0 = out.read_bytes() if rc0 == 0 else b""
            out.unlink()
            rc1, log1 = run([new] + base + new_args + ["--output", out], d, env)
            b1 = out.read_bytes() if rc1 == 0 else b""
            out.unlink()
            check(f"{name}: the same file and log as {BASE_COMMIT}'s extractor ({len(b0)} bytes)",
                  rc0 == rc1 == 0 and b0 == b1 and len(b0) > 1000 and log0 == log1,
                  (rc0, rc1, len(b0), len(b1), log0 == log1))
        rc, log = run([new] + base + C + ["--output", out], d, env)
        sha = hashlib.sha256(out.read_bytes()).hexdigest() if rc == 0 else ""
        check(f"the default's output has the pinned sha256 ({sha[:16]}...)", sha == PINNED_DEFAULT_SHA256, sha)
        head0, rows0 = read_rows(out)
        rc, txt = run(["-c", "import sys, runpy; sys.argv = ['x'] + sys.argv[1:]; "
                             f"runpy.run_path({str(new)!r}, run_name='__main__'); "
                             "print('MEASURED_LOADED', 'measured_setup' in sys.modules)"] + base + C
                      + ["--output", d / "c3.csv"], d, env)
        check("the default does not import measured_setup.py (nor the timestack modules)",
              "MEASURED_LOADED False" in txt, txt[-300:])

        model = per_frame(rows0)
        results = {}
        for src in ("measured", "auto"):
            o = d / f"{src}.csv"
            rc, log = run([new] + base + C + ["--setup-source", src, "--setup-measured", fx["meas"], "--output", o], d, env)
            head, rows = read_rows(o) if rc == 0 else ([], [])
            results[src] = (rc, log, head, per_frame(rows), rows)
        for src in ("measured", "auto"):
            rc, log, head, fr, rows = results[src]
            check(f"{src}: exit 0, the columns of the default plus setup_source, last",
                  rc == 0 and head == head0 + ["setup_source"], (rc, head[-4:], log[-500:]))
            if rc != 0:
                continue
            for cam in ("c1", "c2"):
                for i, want in ((0, 0.25), (1, 0.3), (3, 0.41)):
                    r = fr[(T0 + i * 7200, cam)]
                    check(f"{src}, {cam} frame {i}: the collection's measured setup {want} ("
                          + ("at its start" if i == 0 else "250 s off" if i == 1 else "a buoy outage: no Tp") + ")",
                          float(r["setup_correction_m"]) == want and r["setup_source"] == "measured" and
                          float(r["beach_elevation_navd88"]) == round(float(r["tide_elevation_navd88"]) + want, 4), r)
            for i in (2, 4, 5):
                r = fr[(T0 + i * 7200, "c1")]
                why = {2: "its collection 400 s off", 4: "its collection failed QC", 5: "no collection"}[i]
                if src == "measured":
                    check(f"measured, frame {i} ({why}): no setup, the water level, setup_source none",
                          r["setup_correction_m"] == "" and r["beach_elevation_navd88"] == r["tide_elevation_navd88"]
                          and r["setup_source"] == "none", r)
                else:
                    check(f"auto, frame {i} ({why}): the model's setup, setup_source model",
                          r["setup_correction_m"] == model[(T0 + i * 7200, "c1")]["setup_correction_m"]
                          and r["setup_correction_m"] != "" and r["setup_source"] == "model", r)
            same = len(rows) == len(rows0) and all(
                r[k] == r0[k] for r, r0 in zip(rows, rows0) for k in head0
                if k not in ("setup_correction_m", "beach_elevation_navd88"))
            check(f"{src}: every row and every other column as the default's", same)
            check(f"{src}: the log counts the sources", f"--setup-source {src}: measured 6, "
                  + ("model 0, none 6" if src == "measured" else "model 6, none 0") in log,
                  [ln for ln in log.splitlines() if "Setup source" in ln])
        rc, log, head, fr, _ = results["auto"]
        check("auto: the model's frame without Tp would have no setup; the measured one fills it (H2)",
              model[(T0 + 3 * 7200, "c1")]["setup_correction_m"] == "" and
              fr[(T0 + 3 * 7200, "c1")]["setup_correction_m"] == "0.41")
        # auto with a measured file that lacks the outage's collection: neither a measured setup nor Tp
        part = d / "setup_measured_part.csv"
        ms.write_rows(part, [r for r in csv.DictReader(open(fx["meas"])) if int(r["collection_epoch"]) != T0 + 3 * 7200 - 120])
        o = d / "auto_part.csv"
        rc, log = run([new] + base + C + ["--setup-source", "auto", "--setup-measured", part, "--output", o], d, env)
        fr2 = per_frame(read_rows(o)[1]) if rc == 0 else {}
        check("auto, a frame with neither a measured setup nor Tp: no setup, the water level, setup_source none",
              rc == 0 and all(fr2[(T0 + 3 * 7200, cam)]["setup_correction_m"] == "" and
                              fr2[(T0 + 3 * 7200, cam)]["beach_elevation_navd88"] ==
                              fr2[(T0 + 3 * 7200, cam)]["tide_elevation_navd88"] and
                              fr2[(T0 + 3 * 7200, cam)]["setup_source"] == "none" for cam in ("c1", "c2"))
              and f"--setup-source auto: measured 4, model 6, none 2" in log,
              (rc, [ln for ln in log.splitlines() if "Setup source" in ln]))
        o = d / "s_auto.csv"
        rc, log = run([new] + base + S + ["--setup-source", "auto", "--setup-measured", fx["meas"], "--output", o], d, env)
        head, rows = read_rows(o) if rc == 0 else ([], [])
        fr = per_frame(rows)
        check("auto with the slope-aware model: its columns kept (no slope on a measured frame), setup_source last",
              rc == 0 and head[-6:] == ["setup_correction_m", "beach_elevation_navd88", "setup_beta_f", "setup_model",
                                        "wave_source", "setup_source"]
              and fr[(T0, "c1")]["setup_beta_f"] == "" and fr[(T0, "c1")]["setup_source"] == "measured"
              and fr[(T0 + 2 * 7200, "c1")]["setup_source"] == "model" and fr[(T0 + 2 * 7200, "c1")]["setup_beta_f"] != "",
              (rc, head[-6:], log[-400:]))
        for extra, what in ((C + ["--setup-source", "auto"], "no --setup-measured"),
                            (W + ["--setup-source", "measured", "--setup-measured", fx["meas"]], "no --setup-coef"),
                            (C + ["--setup-source", "auto", "--setup-measured", d / "tide.csv"], "a file of other columns"),
                            (C + ["--setup-source", "auto", "--setup-measured", d / "missing.csv"], "a missing file")):
            o = d / "refused.csv"
            rc, log = run([new] + base + extra + ["--output", o], d, env)
            check(f"{what}: exit 1 before any frame is read", rc == 1 and "ERROR" in log and not o.exists()
                  and "Found" not in log, log[-300:])
        return results["auto"][4], rows0


def ground_of(rows):
    """A ground file from contour rows: 30 points per frame crossing the same cells."""
    head = list(rows[0].keys()) + ["easting_utm19", "northing_utm19"]
    out = []
    frames = sorted({r["source_file"] for r in rows})
    for f in frames:
        r = next(x for x in rows if x["source_file"] == f)
        for j in range(30):
            g = dict(r, easting_utm19=f"{420000.5 + 2.0 * (j % 10):.3f}", northing_utm19=f"{4638000.5 + 2.0 * (j // 10):.3f}")
            out.append(g)
    return head, out


def test_readers(auto_rows, rows0):
    print("the readers of a contour file with setup_source")
    import waterline_consistency as wc
    import survey_products as sp
    from evaluate_setup_model import read_lines
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        src = d / "auto.csv"
        with open(src, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(auto_rows[0].keys()))
            w.writeheader()
            w.writerows(auto_rows)
        cal = HERE / "calibration"
        rc, log = run([HERE / "georectify.py", src, d / "g.csv", "--io-c1", cal / "CACO05_c1_20240801_IO.yaml",
                       "--eo-c1", cal / "CACO05_c1_20251113_EO-CV.yaml", "--io-c2", cal / "CACO05_c2_20240801_IO.yaml",
                       "--eo-c2", cal / "CACO05_c2_20251113_EO-CV.yaml"], d)
        head, rows = read_rows(d / "g.csv") if rc == 0 else ([], [])
        check("georectify.py (pinned, unchanged) keeps setup_source on every row",
              rc == 0 and "setup_source" in head and all(r["setup_source"] in ("measured", "model") for r in rows),
              log[-300:])
        data_ok = True
        try:
            wc.read_points(src)
        except Exception as exc:                     # noqa: BLE001 -- any failure is the finding
            data_ok = False
            print(f"    {exc}")
        n = len(auto_rows)
        wc.write_filtered(src, d / "f.csv", np.ones(n, bool))
        check("waterline_consistency.py reads it by name and copies it through",
              data_ok and read_rows(d / "f.csv")[0] == list(auto_rows[0].keys()))
        L = read_lines(src)
        check("evaluate_setup_model.read_lines reads it by name", len(L["source_file"]) == n)
        out = d / "apply.csv"
        sp.apply_setup(src, out, "c1", 0.037)
        h, r2 = read_rows(out)
        check("survey_products.apply_setup drops setup_source with the setup it recomputes (C's, not 'measured')",
              "setup_source" not in h and h[-2:] == ["setup_correction_m", "beach_elevation_navd88"]
              and all(x["setup_correction_m"] == "" or abs(float(x["setup_correction_m"]) - float(
                  sp.setup_m(x["offshore_hs_m"], x["offshore_tp_s"], 0.037))) < 1e-4 for x in r2))
        got = {}
        for name, rows_in in (("model", rows0), ("auto", auto_rows)):
            gh, grows = ground_of(rows_in)
            g = d / f"ground_{name}.csv"
            with open(g, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=gh)
                w.writeheader()
                w.writerows(grows)
            stem = d / f"dem_{name}"
            rc, log = run([HERE / "dem_from_contours.py", g, stem, "--no-plot", "--no-diagnostics"], d)
            info = json.loads(Path(str(stem) + "_info.json").read_text()) if rc == 0 else {}
            got[name] = (rc, log, info)
        rc, log, info = got["model"]
        check("dem_from_contours.py on the default rows: no setup_source key, no new log line",
              rc == 0 and "setup_source" not in info and "Setup source" not in log)
        rc, log, info = got["auto"]
        check("... on the auto rows: setup_source counts in _info.json and the log",
              rc == 0 and info.get("setup_source") == {"measured": 6, "model": 6}
              and "Setup source      : measured 6 frame(s), model 6 frame(s)" in log, (info.get("setup_source"), rc))
        ok_m, why_m, share_m = ms.dem_provenance(d / "dem_model_dem.asc")
        ok_a, why_a, _s = ms.dem_provenance(d / "dem_auto_dem.asc")
        check("measured_setup.py reads dem_from_contours.py's own _info.json: the model rows' DEM model-built, the "
              "auto rows' DEM refused (rule 2b)", ok_m and share_m > 0 and not ok_a
              and "built with 6 measured frame(s)" in why_a, (ok_m, why_m, share_m, ok_a, why_a))


def test_cron():
    print("waterline_timex_cron.sh: SETUP_SOURCE, and the model-only DEM series while it is on")
    text = (HERE / "waterline_timex_cron.sh").read_text()
    check('the default: SETUP_SOURCE=""', '\nSETUP_SOURCE=""\n' in text)
    a = text.find("    setup_measured_on=0\n")
    block = text[a:text.index("    esac\n", a) + len("    esac\n")] if a >= 0 else ""
    check("one block wires it, the only use of $SETUP_SOURCE outside it its definition",
          bool(block) and 'case "$SETUP_SOURCE" in' in block
          and text.count("$SETUP_SOURCE") == block.count("$SETUP_SOURCE"))
    good = ",".join(ms.FIELDS) + "\n1792500000,x,1792500000.c2.ras.tiff,ok,,0.3\n"
    out = {}
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "archive").mkdir()
        cases = [("", "C", "good"), ("model", "C", "good"), ("Auto", "C", "good"), ("auto", "C", "good"),
                 ("measured", "C", "good"), ("auto", "", "good"), ("auto", "C", "missing"), ("auto", "C", "empty"),
                 ("auto", "C", "malformed"), ("auto", "S", "good")]
        args = {"C": "--setup-coef 0.037", "": "",
                "S": f"--setup-coef 0.037 --setup-model stockdon --setup-slope-dem {d}/archive/dems --setup-stockdon-k 0.3"}
        for src, sa, have in cases:
            f = Path(d) / "archive" / "setup_measured_c2.csv"
            if f.exists():
                f.unlink()
            if have != "missing":
                f.write_text({"good": good, "empty": "", "malformed": "x\n1,2\n"}[have])
            script = (f'log() {{ echo "LOG $*"; }}\nBASE={d}\nDEM_SERIES={d}/archive/dems\n'
                      f'DEM_SERIES_MODEL={d}/archive/dems_model\nSETUP_SOURCE="{src}"\nsetup_arg="{args[sa]}"\n{block}'
                      'echo "ARG $setup_arg"\necho "MODEL $model_setup_arg"\necho "ON $setup_measured_on"\n')
            r = subprocess.run(["bash", "-c", script], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True)
            out[(src, sa, have)] = r.stdout.strip().replace(d, "/b")
    plain = "ARG --setup-coef 0.037\nMODEL --setup-coef 0.037\nON 0"
    check('"" (default) or "model": the options as before, nothing logged, the model-only series off',
          out[("", "C", "good")] == out[("model", "C", "good")] == plain, out[("", "C", "good")])
    check('an unknown value ("Auto"): the options as before, with a WARNING',
          out[("Auto", "C", "good")].endswith(plain) and 'LOG WARNING: SETUP_SOURCE="Auto" is not' in
          out[("Auto", "C", "good")], out[("Auto", "C", "good")])
    au = out[("auto", "C", "good")]
    check('"auto" / "measured" with measured_setup.py\'s file: --setup-source and the file added, logged, the '
          'model-only series on with the model\'s own options',
          "ARG --setup-coef 0.037 --setup-source auto --setup-measured /b/archive/setup_measured_c2.csv" in au
          and "MODEL --setup-coef 0.037\n" in au and au.endswith("ON 1") and "LOG setup: SETUP_SOURCE=auto" in au
          and "--setup-source measured" in out[("measured", "C", "good")], au)
    check("no setup (SETUP_COEF or the waves missing): nothing added, nothing logged, nothing on",
          [ln.strip() for ln in out[("auto", "", "good")].splitlines()] == ["ARG", "MODEL", "ON 0"],
          out[("auto", "", "good")])
    for have in ("missing", "empty", "malformed"):
        o = out[("auto", "C", have)]
        check(f"the measured file {have}: the model's setup and a WARNING (the contour stage never fails for it)",
              o.endswith(plain) and "LOG WARNING: SETUP_SOURCE=auto but" in o and "missing, empty or not" in o, o)
    so = out[("auto", "S", "good")]
    check("with SETUP_MODEL=stockdon the model-only chain takes its slope from its own series (dems_model)",
          "MODEL --setup-coef 0.037 --setup-model stockdon --setup-slope-dem /b/archive/dems_model "
          "--setup-stockdon-k 0.3" in so and "ARG --setup-coef 0.037 --setup-model stockdon --setup-slope-dem "
          "/b/archive/dems --setup-stockdon-k 0.3 --setup-source auto" in so, so)

    # the model-only chain (build_model_series), with a stand-in python3 that records its calls
    fa = text.find("build_model_series() {\n")
    fn = text[fa:text.index("\n}\n", fa) + 3] if fa >= 0 else ""
    check("build_model_series is defined once and called once, only while the measured setup is on",
          bool(fn) and text.count("build_model_series() {") == 1
          and len(re.findall(r"^\s+build_model_series$", text, re.M)) == 1 and
          'if [ "$setup_measured_on" = "1" ] && [ -n "$DEM_WINDOW_DAYS" ]; then\n                    build_model_series\n'
          in text)
    with tempfile.TemporaryDirectory() as d:
        bindir = Path(d) / "bin"
        bindir.mkdir()
        stub = bindir / "python3"
        stub.write_text('#!/bin/bash\necho "CALL $*" >> "$CALLS"\nprev=""\nfor a in "$@"; do\n'
                        '  [ "$prev" = "--output" ] && echo x > "$a"\n  prev="$a"\ndone\n'
                        'case "$1" in *georectify.py) echo x > "$3" ;; esac\nexit 0\n')
        stub.chmod(0o755)
        calls = Path(d) / "calls.txt"
        script = (f'log() {{ echo "LOG $*"; }}\nBASE={d}\nLOG={d}/log.txt\nMODEL_WORK={d}/work\n'
                  f'DEM_SERIES_MODEL={d}/dems_model\nGNSSR_SPLINE=/g/spline.txt\nARCHIVE_CSV=/a/csv\n'
                  f'waves_arg="--waves /a/w.csv"\nmodel_setup_arg="--setup-coef 0.037"\nqc_arg="--gnssr-qc-reference /a/g.csv"\n'
                  f'TIMING_MODE=legacy\nCONSISTENCY_ENABLE=1\nCAL=/cal\nhs_arg="--max-hs 1.5"\nDEM_WINDOW_DAYS=7\n'
                  f'DEM_DAY_BLOCK_LOD=""\nDEM_CELL=2.0\nDEM_MIN_POINTS=3\nDEM_MAX_SPREAD=0.5\n{fn}build_model_series\n')
        env = dict(os.environ, PATH=f"{bindir}:{os.environ.get('PATH', '')}", CALLS=str(calls))
        r = subprocess.run(["bash", "-c", script], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           universal_newlines=True, env=env)
        got = [ln.replace(d, "/b") for ln in calls.read_text().splitlines()] if calls.exists() else []
    want = ["CALL /b/extract_elevation_contours.py /g/spline.txt --processed-dir /a/csv --waves /a/w.csv "
            "--setup-coef 0.037 --gnssr-qc-reference /a/g.csv --timing-mode legacy --output "
            "/b/work/contour_points_model.csv",
            "CALL /b/waterline_consistency.py /b/work/contour_points_model.csv --output "
            "/b/work/contour_points_model_qc.csv --report /b/work/waterline_consistency_report_model.csv",
            "CALL /b/georectify.py /b/work/contour_points_model_qc.csv /b/work/contour_points_ground_model.csv "
            "--io-c1 /cal/CACO05_c1_20240801_IO.yaml --eo-c1 /cal/CACO05_c1_20251113_EO-CV.yaml "
            "--io-c2 /cal/CACO05_c2_20240801_IO.yaml --eo-c2 /cal/CACO05_c2_20251113_EO-CV.yaml",
            "CALL /b/dem_from_contours.py /b/work/contour_points_ground_model.csv /b/work/dem_window --max-hs 1.5 "
            "--last-days 7 --series-dir /b/dems_model --cell 2.0 --min-points 3 --max-spread 0.5 --no-plot"]
    check("the model-only chain: the station's extraction (the model's options, no --setup-source), filter, "
          "georectification and window DEM again, into dems_model; no page, nothing of the station's touched",
          got == want and "LOG model-only window DEM written" in r.stdout, (got, r.stdout[-300:]))
    check("the cron runs neither measured_setup.py nor timestack_wave_height.py",
          'measured_setup.py"' not in text and "python3 \"$BASE/measured_setup" not in text
          and "timestack_wave_height" not in text)
    check("the extractor's command line is as it was", '--processed-dir "$ARCHIVE_CSV" $waves_arg $setup_arg $qc_arg' in text)


# ---------------------------------------------------------------------
# evaluate_measured_setup.py
# ---------------------------------------------------------------------

def test_fit():
    import evaluate_measured_setup as ev
    print("evaluate_measured_setup.py: the C that would match, its day-jackknife interval, the subsets, the reference")
    rng = np.random.default_rng(4)
    colls = []
    for day in range(6):
        for k in range(8):
            w = float(rng.uniform(4, 16))
            colls.append({"m": 0.045 * w + float(rng.normal(0, 0.03)), "w": w, "hs": w / 10, "band": ev.hs_band(w / 10),
                          "phase": "flood" if k % 2 else "ebb", "day": f"2026-09-{26 + day}", "day_offset": day - 3})
    res = ev.compare_model(colls, 0.037)
    a = res["all"]
    check("C_fit 0.045 within 0.002, its 90% day-jackknife interval holds it", abs(a["c_fit"] - 0.045) < 0.002
          and a["c_ci90"][0] < 0.045 < a["c_ci90"][1] and a["c_ci90"][1] - a["c_ci90"][0] > 0, a)
    check("the model (C 0.037) reads low: measured - model positive, = (C_fit - 0.037) x mean w within 0.01",
          a["bias"] > 0 and abs(a["bias"] - (a["c_fit"] - 0.037) * np.mean([c["w"] for c in colls])) < 0.01, a["bias"])
    check("by Hs band, by phase and by day offset", set(res["by_phase"]) == {"flood", "ebb"} and len(res["by_day_offset"])
          == 6 and res["by_hs"])
    two = [c for c in colls if c["day"] in ("2026-09-26", "2026-09-27")]
    check("fewer than 3 days: no interval", all(np.isnan(v) for v in ev.compare_model(two, 0.037)["all"]["c_ci90"]))
    check("fit_c by hand: sum(m w) / sum(w^2)", abs(ev.fit_c([1.0, 2.0], [10.0, 20.0]) - 0.1) < 1e-12)

    # the day jackknife's coverage with the few days a survey window has (the bootstrap it replaced: 0.61-0.64
    # at 3 days); the days share an error, as a day's collections share their beach and waves
    cover = {}
    rng = np.random.default_rng(11)
    for nd in (3, 5):
        hit = 0
        for _ in range(300):
            m, w, dd = [], [], []
            for k in range(nd):
                nk = int(rng.integers(5, 25))
                wk = rng.uniform(6, 20, nk)
                m += list(0.04 * wk + rng.normal(0, 0.05) + rng.normal(0, 0.03, nk))
                w += list(wk)
                dd += [str(k)] * nk
            lo, hi = ev.c_fit_ci(m, w, dd)
            hit += lo <= 0.04 <= hi
        cover[nd] = hit / 300.0
    print(f"    the 90% day-jackknife interval of C_fit covers the true C: {cover[3]:.2f} at 3 days, {cover[5]:.2f} at 5")
    check("its coverage is near the nominal 90% (0.83-0.96) at 3 and 5 days", all(0.83 <= v <= 0.96 for v in cover.values()),
          cover)
    try:
        from scipy import stats as sst
        tab = all(abs(ev.t95(k) - sst.t.ppf(0.95, k)) < 2e-3 for k in range(1, 31))
        tab &= all(ev.t95(k) >= sst.t.ppf(0.95, k) - 2e-3 for k in (35, 50, 90, 200))
        how = "against scipy"
    except ImportError:
        tab, how = abs(ev.t95(2) - 2.920) < 1e-9 and abs(ev.t95(10) - 1.812) < 1e-9 and ev.t95(500) == 1.658, "by hand"
    check(f"Student's t table (no scipy on the NUC), {how}: exact to 30 df, conservative beyond", tab)

    # rule 6: a historical date's collections measured above the ADCP level, where the live model (C 0.037 at
    # the GNSS-R level) is exact: the uncorrected comparison is tilted, the corrected one is not
    ref = {"kind": "historical", "usable": True, "k": 0.010, "x_typ": 10.0, "datum_m": 0.02, "tp_bias_m": 0.02}
    colls2 = [dict(c, m=0.037 * c["w"] + float(ev.reference_delta([c["w"]], ref)[0]) + float(rng.normal(0, 0.01)))
              for c in colls]
    r2 = ev.compare_model(colls2, 0.037, ref)["all"]
    check("rule 6: against the CORRECTED model the bias is ~0 and C_fit at C's reference ~0.037 (its interval holds "
          "it); uncorrected, C_fit is off", abs(r2["ref"]["bias"]) < 0.005 and abs(r2["ref"]["c_fit"] - 0.037) < 0.001
          and r2["ref"]["c_ci90"][0] < 0.037 < r2["ref"]["c_ci90"][1] and abs(r2["c_fit"] - 0.037) > 0.001, r2)
    check("rule 6: delta = k (w - x_typ) - datum - tp_bias, by hand", abs(float(ev.reference_delta([14.0], ref)[0])
          - (0.010 * 4.0 - 0.02 - 0.02)) < 1e-12 and np.isnan(ev.reference_delta([14.0], None)[0]))
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        got = {"none": ev.reference_of(d)}
        for name, fj in (("gnssr", {"water_level": {"method": "GNSS-R (marconi_water_level.py)"}}),
                         ("unmeasured", {"water_level": {"method": "harmonic"},
                                         "waves": {"still_water_reference": {"quantified": False}}}),
                         ("measured", {"water_level": {"method": "harmonic"},
                                       "waves": {"tp_currency_setup_bias_m": 0.0204,
                                                 "still_water_reference": {"quantified": True, "k": 0.011,
                                                                           "share": 0.30, "typical_x": 9.6}}})):
            (d / name / "forcing").mkdir(parents=True)
            (d / name / "forcing" / "forcing.json").write_text(json.dumps(fj))
            got[name] = ev.reference_of(d / name)
    check("reference_of: no forcing -> unknown, not corrected; GNSS-R -> C's own (zero); historical without the "
          "share -> not corrected, cannot count; with it -> k, x_typ, the datum 0.02 and the Tp currency",
          got["none"]["kind"] == "unknown" and not got["none"]["usable"]
          and got["gnssr"]["kind"] == "gnssr" and got["gnssr"]["usable"] and got["gnssr"]["k"] == 0
          and not got["unmeasured"]["usable"] and "cannot count" in got["unmeasured"]["note"]
          and got["measured"]["usable"] and got["measured"]["k"] == 0.011 and got["measured"]["x_typ"] == 9.6
          and got["measured"]["datum_m"] == 0.02 and got["measured"]["tp_bias_m"] == 0.0204, got)

    class Lv:
        def at(self, t):
            t = np.asarray(t, float)
            return 0.001 * (t - t[0]), None
    check("the tide's phase: rising flood, falling ebb, flat slack", ev.tide_phase(Lv(), 1e9) == "flood")


def test_blind():
    print("evaluate_measured_setup.py: a BLIND row and a look-alike survey file are never opened")
    import evaluate_measured_setup as ev
    import test_evaluate_setup_model as TE
    with tempfile.TemporaryDirectory() as d:
        root, sdir, table = TE.fixture(Path(d) / "fx", ("c", 0.045))
        opened = []

        def audit(event, args):
            if event == "open" and args and isinstance(args[0], (str, bytes, os.PathLike)):
                opened.append(os.fsdecode(args[0]))
        sys.addaudithook(audit)
        empty = Path(d) / "noras"
        empty.mkdir()
        rc, out = quiet(ev.main, ["--output-root", str(root), "--surveys", str(table), "--survey-dirs", str(sdir),
                                  "--ras-dir", str(empty), "--out-dir", str(Path(d) / "eval")])
        bad = [p for p in opened if "held_blind" in p or "lookalike" in p or "2024-10-23" in p]
        check("the BLIND row skipped, nothing of it or of the look-alike opened; no stacks: exit 4",
              rc == 4 and "2024-10-23 oct_points: BLIND" in out and "named like a survey of a BLIND date" in out
              and not bad and "no c2 stacks of" in out, (rc, bad, out[-600:]))


def test_lines_known_setup():
    print("rule 5 with a known setup: lines planted at C = 0.045, scored with the model (0.037) and the measured")
    import evaluate_measured_setup as ev
    import test_evaluate_setup_model as TE
    from evaluate_setup_model import load_product
    from setup_model import setup_c
    with tempfile.TemporaryDirectory() as d:
        root, sdir = Path(d) / "p", Path(d) / "s"
        root.mkdir()
        sdir.mkdir()
        row = TE.make_product(root, sdir, "2031-01-23", "jan_dsm", 0.13, ("c", 0.045), "dsm", 1)
        sv, why = quiet(load_product, root / "2031-01-23", row, sdir / row["path"], HERE / "calibration")[0]
        mrows = []
        for i, f in enumerate(sv.frames):
            e = ev.frame_epoch(f["name"])
            r = {k: "" for k in ms.FIELDS}
            # every other frame measured: the truth's setup (the measured setup IS the line's, here)
            r.update(collection_epoch=e, filename=f"{e}.c2.ras.tiff", status="ok" if i % 2 == 0 else "failed qc",
                     setup_m=round(float(setup_c(f["hs"], f["tp"], 0.045)), 3) if i % 2 == 0 else "")
            mrows.append(r)
        mcsv = Path(d) / "m.csv"
        ms.write_rows(mcsv, mrows)
        res = ev.score_lines(sv, ms.MeasuredSetups(mcsv), 0.037)
        mo, me = res["model"], res["measured"]
        print(f"    {res['compared']} frames: model bias {mo['bias']:+.3f} RMSE {mo['rmse']:.3f} m; measured bias "
              f"{me['bias']:+.3f} RMSE {me['rmse']:.3f} m; every frame: model {res['model_all']['rmse']:.3f}, auto "
              f"{res['auto_all']['rmse']:.3f} m")
        check("the frames with a measured setup are compared (about half)", 8 <= res["compared"] <= 16, res["compared"])
        check("measured: at the noise (|bias| < 0.02, RMSE < 0.04 m)", abs(me["bias"]) < 0.02 and me["rmse"] < 0.04, me)
        check("model 0.037 < truth 0.045: the lines read low (bias < -0.02) and worse",
              mo["bias"] < -0.02 and mo["rmse"] > me["rmse"], mo)
        check("auto over every frame scored: between the model and the measured",
              me["rmse"] <= res["auto_all"]["rmse"] < res["model_all"]["rmse"])
        p = res["paired"]
        print(f"    paired, measured - {p['against'][:40]}...: RMSE {p['rmse'][0]:+.3f} [{p['rmse'][1]:+.3f}, "
              f"{p['rmse'][2]:+.3f}] over {p['days']} days")
        check("the paired day-jackknife interval: measured better than the uncorrected model beyond its interval "
              "(its upper end below 0), and it says the model is not corrected (no forcing record)",
              p["days"] >= 3 and p["rmse"][2] < 0 and p["absbias"][2] < 0 and "NOT corrected" in p["against"]
              and res["model_ref"] is None, p)

        # rule 6 with a known answer: a historical product whose forcing record makes the correction exactly the
        # gap between C 0.037 and the truth 0.045 (k 0.008, x_typ 10 m, tp_bias -0.10 m: delta = 0.008 w)
        (root / "2031-01-23" / "forcing").mkdir()
        (root / "2031-01-23" / "forcing" / "forcing.json").write_text(json.dumps(
            {"water_level": {"method": "harmonic"},
             "waves": {"tp_currency_setup_bias_m": -0.10,
                       "still_water_reference": {"quantified": True, "k": 0.008, "share": 0.216, "typical_x": 10.0}}}))
        ref = ev.reference_of(root / "2031-01-23")
        res = ev.score_lines(sv, ms.MeasuredSetups(mcsv), 0.037, ref=ref)
        mr, p = res["model_ref"], res["paired"]
        print(f"    with the reference corrected: model+ref bias {mr['bias']:+.3f} RMSE {mr['rmse']:.3f} m; paired "
              f"RMSE {p['rmse'][0]:+.3f} [{p['rmse'][1]:+.3f}, {p['rmse'][2]:+.3f}]")
        check("rule 6: the model corrected to C's reference scores at the noise where the uncorrected one reads low",
              ref["usable"] and abs(mr["bias"]) < 0.02 and mr["rmse"] < 0.04 and res["model"]["bias"] < -0.02, mr)
        check("... and the paired comparison is against the corrected model: no difference, its interval across 0",
              p["against"] == "the model corrected to C's reference" and abs(p["rmse"][0]) < 0.01
              and p["rmse"][1] <= 0.0 <= p["rmse"][2] + 1e-9, p)


def test_survey_surface():
    print("evaluate_measured_setup.py end to end: the setup measured on a survey's own surface (known setup)")
    import evaluate_measured_setup as ev
    from view_reproject import ground_to_pixel
    import georectify
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = scene(tmp, beta=0.10)
        truth = float(sc["z_rel"].mean())
        date = "2026-10-20"
        # the survey: the planar beach as a 1 m DSM over -2.5..4 m
        sdir = tmp / "surveys"
        sdir.mkdir()
        x0, y0, nx, ny = 420000.0, 4638300.0, 300, 300
        XX, YY = np.meshgrid(x0 + np.arange(nx) + 0.5, (y0 + np.arange(ny) + 0.5)[::-1])
        z = -0.10 * thw.seaward_distance(XX, YY)
        z[(z < -2.5) | (z > 4.0)] = -9999
        with open(sdir / "synthetic_dsm.asc", "w") as f:
            f.write(f"ncols {nx}\nnrows {ny}\nxllcorner {x0}\nyllcorner {y0}\ncellsize 1\nNODATA_value -9999\n")
            for r in z:
                f.write(" ".join(f"{v:.4f}" for v in r) + "\n")
        # the product: provenance, forcing (the scene's level), waterlines planted 0.15 m BELOW the truth's setup
        # for the model and at the truth for the measured (Hs, Tp chosen so the model's setup is truth - 0.15)
        pdir = tmp / "products" / date
        (pdir / "waterlines").mkdir(parents=True)
        (pdir / "forcing").mkdir()
        shutil.copy(sc["level"], pdir / "forcing" / "water_level.csv")
        tp = 9.0
        w_want = (truth - 0.15) / 0.037
        hs = round(w_want ** 2 / (9.81 * tp ** 2 / (2 * np.pi)), 2)
        io_ = georectify.load_intrinsics(HERE / "calibration" / "CACO05_c2_20240801_IO.yaml")
        eo_ = georectify.load_extrinsics(HERE / "calibration" / "CACO05_c2_20251113_EO-CV.yaml")
        sw = 0.3
        zl = sw + truth
        sa, al = np.radians(81.0), np.radians(351.0)
        xs = -zl / 0.10
        a = np.arange(60.0, 240.0, 0.5)
        E = 420150.0 + xs * np.sin(sa) + a * np.sin(al)
        N = 4638400.0 + xs * np.cos(sa) + a * np.cos(al)
        U, V, ok = ground_to_pixel(E, N, np.full(len(E), zl), io_, eo_)
        tu = datetime.fromtimestamp(sc["epoch"], tz=timezone.utc).isoformat()
        wfile = pdir / "waterlines" / "contour_points_ground_filtered.csv"
        with open(wfile, "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(["source_file", "camera", "capture_time_utc", "capture_epoch", "pixel_column", "pixel_row",
                         "tide_elevation_navd88", "offshore_hs_m", "offshore_tp_s", "wave_source"])
            for u, v in zip(U[ok], V[ok]):
                wr.writerow([f"{sc['epoch']}.c2.timex", "c2", tu, sc["epoch"], f"{u:.3f}", f"{v:.3f}", sw, hs, tp,
                             "marconi:hs_best/tp_s"])
        prov = {"date": date, "status": "complete", "blind": None, "filter": {"file": str(wfile)},
                "cameras": {"c2": {"window": [date, date], "calibration": {
                    "eo_path": str(HERE / "calibration" / "CACO05_c2_20251113_EO-CV.yaml"),
                    "io_file": "CACO05_c2_20240801_IO.yaml"}}}}
        (pdir / "provenance.json").write_text(json.dumps(prov))
        table = tmp / "surveys.csv"
        table.write_text("date,name,path,survey_type,survey_date,label,why\n"
                         f"{date},synthetic_dsm,synthetic_dsm.asc,dsm,{date},INDEPENDENT,synthetic\n")
        js = tmp / "r.json"
        import runup_from_timestack as rft
        rft._PIX_FILES.clear()
        rc, out = quiet(ev.main, ["--output-root", str(tmp / "products"), "--surveys", str(table), "--survey-dirs",
                                  str(sdir), "--ras-dir", str(sc["ras"]), "--pix", str(sc["pix"]), "--lines", "1",
                                  "--owg-csv", str(sc["owg"]), "--out-dir", str(tmp / "eval"), "--json", str(js)])
        rft._PIX_FILES.clear()
        res = json.loads(js.read_text())[0] if js.exists() else {}
        mrows = list(csv.DictReader(open(res["measured"]["csv"]))) if res else []
        check("one collection measured on the survey's surface, QC passed, the still water the product's",
              rc == 0 and len(mrows) == 1 and mrows[0]["status"] == "ok"
              and mrows[0]["water_level_source"] == "forcing water level"
              and mrows[0]["dem"] == "synthetic_dsm.asc" and mrows[0]["dem_window"].startswith("survey DSM"),
              (rc, mrows, out[-800:]))
        if not mrows or mrows[0]["status"] != "ok":
            return
        m = float(mrows[0]["setup_m"])
        print(f"    truth {truth:.3f} m, measured on the survey {m:.3f} m; model {res['collections'][0]['model']:.3f} m")
        check("the measured setup within the detector's bias of the truth (+0.00 to +0.12 m)", 0.0 <= m - truth <= 0.12,
              m - truth)
        dems = tmp / "dems"
        dems.mkdir()
        write_plane_dem(dems / "dem_2026-10-19_7d_dem.asc", 0.10, 0.0, -2.5, 4.0)
        rft._PIX_FILES.clear()
        live = measure(sc, dems)
        check("the same surface as a station DEM of the day before gives the same setup (within 5 mm)",
              isinstance(live, dict) and abs(float(live["setup_m"]) - m) <= 0.005, live if not isinstance(live, dict)
              else live["setup_m"])
        c = res["collections"][0]
        cm = res["model_compare"]["all"]
        check("joined to the product's frame: C_fit = m / w by hand, the model C w, no interval on one day",
              abs(cm["c_fit"] - m / c["w"]) < 1e-9 and abs(c["model"] - 0.037 * c["w"]) < 1e-9
              and (cm["c_ci90"][0] is None or np.isnan(cm["c_ci90"][0])), (cm, c))
        check("no forcing.json in this product: the reference is unknown, the model not corrected, and it says so",
              res["reference"]["kind"] == "unknown" and not res["model_compare"]["reference_corrected"]
              and "reference (rule 6): no readable forcing" in out, res["reference"])
        check("the model planted 0.15 m below the truth reads it: measured - model 0.15-0.27 m",
              0.14 <= cm["bias"] <= 0.28, cm["bias"])
        sc_ = res.get("lines") or {}
        print(f"    the line on the survey: model median {sc_.get('model', {}).get('median')}, measured "
              f"{sc_.get('measured', {}).get('median')}")
        check("the line planted at the truth's setup: on the survey it reads nearer with the measured setup than "
              "with the model's", sc_.get("compared") == 1 and abs(sc_["measured"]["median"]) <
              abs(sc_["model"]["median"]) and sc_["model"]["median"] < -0.05, sc_)
        check("the report prints both comparisons", "MEASURED against the MODEL" in out and
              "THE LINES ON THE SURVEY" in out)

        # a date none of whose stacks gives a row (a line the .pix lacks): skipped with the reason, no traceback
        rft._PIX_FILES.clear()
        try:
            rc7, out7 = quiet(ev.main, ["--output-root", str(tmp / "products"), "--surveys", str(table),
                                        "--survey-dirs", str(sdir), "--ras-dir", str(sc["ras"]), "--pix", str(sc["pix"]),
                                        "--lines", "7", "--owg-csv", str(sc["owg"]), "--out-dir", str(tmp / "eval7")])
            crashed = None
        except Exception as exc:                         # noqa: BLE001 -- a crash is the finding
            rc7, out7, crashed = None, "", repr(exc)
        rft._PIX_FILES.clear()
        check("no stack readable along --lines 7: the date skipped with the reason (check --pix / --lines), exit 4, "
              "no traceback", crashed is None and rc7 == 4 and "none of the 1 stack(s) could be read" in out7,
              (crashed, rc7, out7[-400:]))

        # rule 3 for a period without the OWG record: the window's sharpness scored from its own images
        import cv2
        imgs = tmp / "images"
        imgs.mkdir()
        rng = np.random.default_rng(3)
        texture = cv2.GaussianBlur(rng.integers(0, 255, (512, 612)).astype(np.uint8), (0, 0), 1.0)
        cv2.imwrite(str(imgs / f"{sc['epoch']}.c2.snap.jpg"), texture)
        cv2.imwrite(str(imgs / f"{sc['epoch']}.c2.timex.jpg"), cv2.GaussianBlur(texture, (0, 0), 6.0))
        rft._PIX_FILES.clear()
        rcs, outs = quiet(ev.main, ["--output-root", str(tmp / "products"), "--surveys", str(table), "--survey-dirs",
                                    str(sdir), "--ras-dir", str(sc["ras"]), "--pix", str(sc["pix"]), "--lines", "1",
                                    "--owg-csv", str(tmp / "no_owg.csv"), "--sharpness-images", str(imgs),
                                    "--out-dir", str(tmp / "evals"), "--json", str(tmp / "s.json")])
        rft._PIX_FILES.clear()
        shp = tmp / "evals" / f"{date}_synthetic_dsm_sharpness.csv"
        srows = list(csv.DictReader(open(shp))) if shp.exists() else []
        res_s = json.loads((tmp / "s.json").read_text())[0] if (tmp / "s.json").exists() else {}
        mrow = list(csv.DictReader(open(res_s["measured"]["csv"])))[0] if res_s else {}
        check("--sharpness-images: the collection's snap scored (preferred to its timex) and used as the window's "
              "sharpness; the row passes on it, with no --allow-no-sharpness",
              rcs == 0 and len(srows) == 1 and srows[0]["product"] == "snap"
              and mrow.get("status") == "ok" and float(mrow["window_sharpness"]) == float(srows[0]["sharpness"])
              and "scored here from the period's own c2 images" in outs and "snap " in outs, (rcs, srows, mrow, outs[-500:]))


def test_py38():
    print("Python 3.8 syntax")
    for f in ("measured_setup.py", "evaluate_measured_setup.py", "test_measured_setup.py",
              "extract_elevation_contours.py", "dem_from_contours.py", "timestack_wave_height.py", "survey_products.py"):
        try:
            ast.parse((HERE / f).read_text(encoding="utf-8"), filename=f, feature_version=(3, 8))
            good = True
        except SyntaxError as exc:
            good = False
            print(f"    {exc}")
        check(f"{f} parses as Python 3.8", good)


def main():
    test_qc()
    test_combine()
    test_sampling_se()
    test_prior_dem()
    test_csv()
    test_scene()
    test_chain()
    auto_rows, rows0 = test_extractor()
    test_readers(auto_rows, rows0)
    test_cron()
    test_fit()
    test_blind()                # before any other run of the evaluation: it arms the blind guard
    test_lines_known_setup()
    test_survey_surface()
    test_py38()
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        raise SystemExit(1)
    print("all measured setup tests passed")


if __name__ == "__main__":
    main()

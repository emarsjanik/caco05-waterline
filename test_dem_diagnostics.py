#!/usr/bin/env python3
"""
Tests For The DEM's Reporting Additions
==========================================
Plain python3 (no pytest needed on the station):

    python3 test_dem_diagnostics.py [--keep DIR]

WHAT. The additions that report on a DEM without changing any value in it:
  * the wave currency (extract_elevation_contours.py wave_source, and the
    DEM's log and <stem>_info.json wave_currency): which wave columns of
    which kind of file tagged the frames, and a WARNING when that is not
    the currency the setup coefficient was fitted in (the raw buoy archive
    the cron falls back to when waves_marconi.csv is missing);
  * the DEM's grids are the same with them as without: the elevations,
    spread and count grids do not depend on the new columns.

HOW. Synthetic detector CSVs, tide model and wave files, and a synthetic
ground file, in a temporary folder (--keep DIR to look at it); the
scripts run as the cron runs them.

STATION FILES. None are read.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

FAILED = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILED.append(msg)


def run(args, d):
    r = subprocess.run([sys.executable] + [str(a) for a in args], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
    return r.returncode, r.stdout


# ---------------------------------------------------------------------
# Wave currency
# ---------------------------------------------------------------------

def test_wave_columns():
    print("extract_elevation_contours.wave_columns: which columns, which currency")
    import extract_elevation_contours as eec
    check(eec.wave_columns(["time_utc", "epoch", "hs_best", "hs_sigma", "tp_s"])
          == ("hs_best", "tp_s", "marconi:hs_best/tp_s") == ("hs_best", "tp_s", eec.WAVE_CURRENCY_FITTED),
          "waves_marconi.csv (hs_best, tp_s): the currency C was fitted in")
    check(eec.wave_columns(["time_utc", "epoch", "wvht_m", "dpd_s", "apd_s", "mwd_deg"])
          == ("wvht_m", "dpd_s", "buoy:wvht_m/dpd_s"),
          "fetch_buoy_waves.py's archive (wvht_m, dpd_s): the raw buoy")
    check(eec.wave_columns(["time_utc", "epoch", "wvht_m", "dpd_s", "mwd_deg", "source", "hs_sigma_m",
                            "tp_sigma_s"]) == ("wvht_m", "dpd_s", "adcp_equiv:wvht_m/dpd_s"),
          "historical_forcing.py's waves.csv (wvht_m, dpd_s, hs_sigma_m): ADCP-equivalent")
    ok, _ = eec.wave_currency_note(eec.WAVE_CURRENCY_FITTED)
    nb, tb = eec.wave_currency_note("buoy:wvht_m/dpd_s")
    na, ta = eec.wave_currency_note("adcp_equiv:wvht_m/dpd_s")
    check(ok and not nb and not na and "NOT the currency" in tb and "ADCP-equivalent Tp" in ta,
          "wave_currency_note: only marconi:hs_best/tp_s is the fitted currency")


def detector_csvs(folder, t0, n=8, step_s=1800):
    """n c1 frames, 30 min apart, of 400 columns each (rows with a little noise, Has_Signal=1)."""
    folder.mkdir(parents=True, exist_ok=True)
    rng = np.random.RandomState(0)
    epochs = []
    for i in range(n):
        ep = t0 + i * step_s
        epochs.append(ep)
        t = datetime.fromtimestamp(ep, tz=timezone.utc)
        name = f"{ep}.{t:%a.%b.%d_%H_%M_%S}.GMT.{t:%Y}.CACO05.c1.timex.csv"
        lines = ["Column,Row,Confidence,Row_Precise,Peak_Sharpness,Has_Signal"]
        for c in range(400):
            row = 500.0 + 10.0 * i + 3.0 * np.sin(c / 7.0) + rng.normal(0, 0.8)
            lines.append(f"{c},{int(row)},0.95,{row:.3f},0.5,1")
        (folder / name).write_text("\n".join(lines) + "\n")
    return epochs


def test_extract_wave_source(d):
    print("extract_elevation_contours.py --waves: wave_source on every row, and the log")
    d = Path(d) / "extract"
    t0 = int(datetime(2026, 10, 12, 14, 0, tzinfo=timezone.utc).timestamp())
    epochs = detector_csvs(d / "processed", t0)
    tide = ["time,level"] + [f"{datetime.fromtimestamp(t0 + k * 600, tz=timezone.utc):%Y-%m-%d %H:%M:%S},"
                             f"{0.2 + 0.0002 * k * 600 / 60:.4f}" for k in range(-6, 60)]
    (d / "tide.csv").write_text("\n".join(tide) + "\n")
    hours = range(t0 - 7200, t0 + 6 * 3600, 3600)
    (d / "waves_marconi.csv").write_text(
        "time_utc,epoch,hs_best,hs_sigma,tp_s\n"
        + "".join(f"x,{e},1.00,0.1,8.0\n" for e in hours))
    (d / "waves_44008.csv").write_text(
        "time_utc,epoch,wvht_m,dpd_s,apd_s,mwd_deg\n"
        + "".join(f"x,{e},1.60,9.0,6.0,120\n" for e in hours))
    out = {}
    for name in ("waves_marconi.csv", "waves_44008.csv"):
        dst = d / f"contours_{name}"
        rc, log = run([HERE / "extract_elevation_contours.py", d / "tide.csv", "--processed-dir",
                       d / "processed", "--waves", d / name, "--setup-coef", "0.037", "--output", dst], d)
        (d / f"{name}.log").write_text(log)
        import csv
        rows = list(csv.DictReader(open(dst))) if dst.exists() else []
        out[name] = (rc, log, rows)
    rc, log, rows = out["waves_marconi.csv"]
    check(rc == 0 and rows and all(r["wave_source"] == "marconi:hs_best/tp_s" for r in rows)
          and list(rows[0])[-1] == "wave_source" and list(rows[0])[-3:-1] == ["setup_correction_m",
                                                                             "beach_elevation_navd88"],
          f"waves_marconi.csv: {len(rows)} rows, wave_source marconi:hs_best/tp_s, the last column "
          f"(after the setup columns, so nothing positional moves)")
    wave_warn = lambda text: [ln for ln in text.splitlines() if "WARNING" in ln and "wave" in ln.lower()]
    check("Wave file                : " in log and "(columns hs_best, tp_s)" in log and not wave_warn(log),
          "its log names the file and columns, no wave WARNING")
    rc, log, rows = out["waves_44008.csv"]
    warn = wave_warn(log)
    check(rc == 0 and rows and all(r["wave_source"] == "buoy:wvht_m/dpd_s" for r in rows)
          and len(warn) == 1 and "raw buoy" in warn[0] and "C = 0.037" in warn[0],
          "the raw 44008 archive: wave_source buoy:wvht_m/dpd_s and one WARNING line naming C")
    # the setup itself is unchanged by the tag: C*sqrt(Hs*g*Tp^2/2pi) from the row's own Hs, Tp
    r = out["waves_marconi.csv"][2][0]
    want = 0.037 * np.sqrt(1.0 * 9.81 * 8.0 ** 2 / (2 * np.pi))
    check(abs(float(r["setup_correction_m"]) - want) < 1e-4,
          f"setup unchanged: {float(r['setup_correction_m']):.4f} m (C sqrt(Hs L0) = {want:.4f})")
    return epochs


HEADER = ("source_file,camera,capture_time_utc,capture_epoch,pixel_column,easting_utm19,northing_utm19,"
          "tide_elevation_navd88,offshore_hs_m,offshore_tp_s,setup_correction_m,beach_elevation_navd88")


def ground_rows(frames, wave_source=None):
    """frames: [(camera, epoch, tide, hs, dz)]: each frame a line of 20 points along the beach row
    whose elevation is its level (tide + 0.05 setup + dz); the beach rises 0.2 m per 2 m row."""
    head = HEADER + (",wave_source" if wave_source is not None else "")
    rows = [head]
    for cam, ep, tide, hs, dz in frames:
        t = datetime.fromtimestamp(ep, tz=timezone.utc)
        f = f"{ep}.{cam}.timex"
        z = tide + 0.05 + dz
        x = int(round((z + 0.5) / 0.2))                    # the row whose beach is at z
        x0 = 420001.0 if cam == "c1" else 420021.0         # c1 columns 0-19, c2 10-29: 10 shared
        for col in range(20):
            rows.append(f"{f},{cam},{t:%Y-%m-%dT%H:%M:%S}+00:00,{ep},{col},{x0 + 2 * col:.2f},"
                        f"{4638001.0 + 2 * x:.2f},{tide:.4f},{hs:.2f},8.0,0.05,{z:.4f}"
                        + (f",{wave_source}" if wave_source is not None else ""))
    return rows


def tidal_frames(days=5, seed=3, c2_offset=0.0, lag_m_per_mph=0.0, hs_slope=0.0):
    """Both cameras every 30 min, 13-21 UTC, on `days` days from 1 Oct 2026, under a 12.42 h tide of
    1.2 m amplitude. Each frame's line sits at its level + noise + lag x dz/dt (m/h) + hs_slope x
    (Hs - 1) (+ c2_offset on c2): the planted water-level timing error and wave dependence."""
    rng = np.random.RandomState(seed)
    out = []
    w = 2 * np.pi / (12.42 * 3600)
    for day in range(days):
        hs = 0.6 + 0.2 * day
        for k in range(17):
            ep = int(datetime(2026, 10, 1 + day, 13, 0, tzinfo=timezone.utc).timestamp()) + 1800 * k
            tide = 0.6 + 1.2 * np.sin(w * ep)
            rate = 1.2 * w * np.cos(w * ep) * 3600                 # m/h
            for cam in ("c1", "c2"):
                dz = rng.normal(0, 0.02) + lag_m_per_mph * rate + hs_slope * (hs - 1.0) \
                     + (c2_offset if cam == "c2" else 0.0)
                out.append((cam, ep, tide, hs, dz))
    return out


def test_dem_wave_currency(d):
    print("dem_from_contours.py: the wave currency in the log and <stem>_info.json")
    d = Path(d) / "dem_waves"
    d.mkdir(parents=True, exist_ok=True)
    frames = tidal_frames()
    res = {}
    for tag in (None, "marconi:hs_best/tp_s", "buoy:wvht_m/dpd_s"):
        name = (tag or "none").split(":")[0]
        g = d / f"ground_{name}.csv"
        g.write_text("\n".join(ground_rows(frames, tag)) + "\n")
        rc, log = run([HERE / "dem_from_contours.py", g, d / name, "--no-plot"], d)
        (d / f"{name}.log").write_text(log)
        info = json.loads((d / f"{name}_info.json").read_text()) if rc == 0 else {}
        res[name] = (rc, log, info)
    nframes = len({(c, e) for c, e, *_ in frames})
    rc, log, info = res["marconi"]
    wc = info.get("wave_currency", {})
    check(rc == 0 and wc.get("frames") == {"marconi:hs_best/tp_s": nframes} and wc.get("all_fitted_currency") is True
          and wc.get("fitted_currency") == "marconi:hs_best/tp_s" and "WARNING" not in log
          and "Wave currency     : marconi:hs_best/tp_s" in log,
          f"marconi rows: info.json wave_currency counts {nframes} frames in the fitted currency; no WARNING")
    rc, log, info = res["buoy"]
    warn = [ln for ln in log.splitlines() if "WARNING" in ln]
    check(rc == 0 and info["wave_currency"]["frames"] == {"buoy:wvht_m/dpd_s": nframes}
          and info["wave_currency"]["all_fitted_currency"] is False
          and len(warn) == 1 and "raw buoy" in warn[0],
          "raw-buoy rows with a setup: one WARNING line, and info.json says not the fitted currency")
    rc, log, info = res["none"]
    check(rc == 0 and info["wave_currency"]["frames"] == {"unrecorded": nframes} and "WARNING" not in log,
          "contours from before wave_source: 'unrecorded', no WARNING")
    grids = [np.loadtxt(d / f"{n}_{k}.asc", skiprows=6) for n in ("none", "marconi", "buoy")
             for k in ("dem", "spread", "count")]
    check(all(np.array_equal(grids[k], grids[3 + k]) and np.array_equal(grids[k], grids[6 + k]) for k in range(3)),
          "the dem, spread and count grids are identical with and without the column")


def test_py38():
    print("Python 3.8 syntax")
    import ast
    for f in ("dem_from_contours.py", "extract_elevation_contours.py", "test_dem_diagnostics.py"):
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
    d = os.path.realpath(tempfile.mkdtemp(prefix="dem_diagnostics_", dir=args.keep))
    print(f"test folder: {d}")
    try:
        test_wave_columns()
        test_extract_wave_source(d)
        test_dem_wave_currency(d)
        test_py38()
    finally:
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED")
        return 1
    print("all DEM diagnostics tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

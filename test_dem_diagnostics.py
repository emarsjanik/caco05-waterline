#!/usr/bin/env python3
"""
Tests For The DEM's Reporting Additions
==========================================
Plain python3 (no pytest needed on the station):

    python3 test_dem_diagnostics.py [--keep DIR]

WHAT. The additions that report on a DEM without changing any value in it:
  * the wave currency (extract_elevation_contours.py wave_source, and the
    DEM's log and <stem>_info.json wave_currency): which wave columns of
    which kind of file tagged the frames, and a WARNING when that is the
    raw buoy archive the cron falls back to when waves_marconi.csv is
    missing (in the DEM's log only when its elevations carry a setup, a
    NOTE otherwise); a NOTE, not the raw buoy, for a file whose kind its
    columns do not tell, such as process_chelsea.py's ADCP record;
  * the residuals against the tide's stage and the waves (the DEM's log
    and <stem>_info.json diagnostics.residuals): a planted water-level
    timing error (each line given an elevation 0.1 m per m/h of dz/dt too
    high, the line 6 min older than its level) comes back as that slope,
    flood frames high and ebb frames low; a planted wave dependence comes
    back against Hs; none planted, none found;
  * <stem>_sigma.asc: 1.25 x (spread/2) / sqrt(frames) on every cell
    with a value, NaN elsewhere (and on a cell of one frame, which has no
    scatter to measure), in the series copy too, convertible to a GeoTIFF
    whose description says what it is;
  * the c1 - c2 offset in the cells both cameras cross (the DEM's log and
    diagnostics.camera_offset): a planted 0.10 m comes back; none, none;
    one camera, None;
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
    # process_chelsea.py's adcp_waves.csv (the ADCP's own wh_4061 and peak period) has wvht_m / dpd_s
    # but neither apd_s nor hs_sigma_m: it was called the raw buoy, with a WARNING in its logs
    check(eec.wave_columns(["time_utc", "epoch", "wvht_m", "dpd_s", "mwd_deg"])
          == ("wvht_m", "dpd_s", "unknown:wvht_m/dpd_s"),
          "process_chelsea.py's adcp_waves.csv (wvht_m, dpd_s, no apd_s): unknown, not the raw buoy")
    ok, _ = eec.wave_currency_note(eec.WAVE_CURRENCY_FITTED)
    nb, tb = eec.wave_currency_note("buoy:wvht_m/dpd_s")
    na, ta = eec.wave_currency_note("adcp_equiv:wvht_m/dpd_s")
    nu, tu = eec.wave_currency_note("unknown:wvht_m/dpd_s")
    check(ok and not nb and not na and not nu and "NOT the currency" in tb and "ADCP-equivalent Tp" in ta
          and "do not name" in tu and "process_chelsea.py" in tu,
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
    # process_chelsea.py's adcp_waves.csv, as it writes it (the ADCP's own record)
    (d / "adcp_waves.csv").write_text(
        "time_utc,epoch,wvht_m,dpd_s,mwd_deg\n"
        + "".join(f"x,{e},1.10,8.5,120\n" for e in hours))
    out = {}
    import csv
    for name, setup in (("waves_marconi.csv", True), ("waves_44008.csv", True), ("adcp_waves.csv", True),
                        ("waves_44008.csv", False)):
        dst = d / f"contours_{name}_{setup}"
        rc, log = run([HERE / "extract_elevation_contours.py", d / "tide.csv", "--processed-dir",
                       d / "processed", "--waves", d / name] + (["--setup-coef", "0.037"] if setup else [])
                      + ["--output", dst], d)
        (d / f"{name}_{setup}.log").write_text(log)
        rows = list(csv.DictReader(open(dst))) if dst.exists() else []
        out[name if setup else name + " (no setup)"] = (rc, log, rows)
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
    rc, log, rows = out["adcp_waves.csv"]
    check(rc == 0 and rows and all(r["wave_source"] == "unknown:wvht_m/dpd_s" for r in rows)
          and not wave_warn(log) and "NOTE: the wave file's kind is not known" in log,
          "process_chelsea.py's ADCP record: wave_source unknown:wvht_m/dpd_s, a NOTE, no WARNING (it was "
          "tagged as the raw buoy, with the WARNING that waves_marconi.csv was missing)")
    # without a setup the fallback is still the alarm (the cron says nothing else when waves_marconi.csv
    # is missing); the DEM, whose elevations then carry no setup, says it as a NOTE
    rc, log, rows = out["waves_44008.csv (no setup)"]
    warn = wave_warn(log)
    check(rc == 0 and rows and all(r["wave_source"] == "buoy:wvht_m/dpd_s" for r in rows)
          and "setup_correction_m" not in rows[0] and len(warn) == 1 and "raw buoy" in warn[0]
          and "C = " not in warn[0] and "waves_marconi.csv" in warn[0],
          "the raw 44008 archive with no setup: still one WARNING (the fallback), naming no C")
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
    where the beach is at the frame's true level (tide + 0.05 m setup; the beach rises 0.05 m per
    2 m row), given the elevation true level + dz -- dz is the error in the elevation the line is
    given, as a timing error or a camera bias would make it."""
    head = HEADER + (",wave_source" if wave_source is not None else "")
    rows = [head]
    for cam, ep, tide, hs, dz in frames:
        t = datetime.fromtimestamp(ep, tz=timezone.utc)
        f = f"{ep}.{cam}.timex"
        true = tide + 0.05
        x = int(round((true + 1.0) / 0.05))                # the row whose beach is at the true level
        x0 = 420001.0 if cam == "c1" else 420021.0         # c1 columns 0-19, c2 10-29: 10 shared
        for col in range(20):
            rows.append(f"{f},{cam},{t:%Y-%m-%dT%H:%M:%S}+00:00,{ep},{col},{x0 + 2 * col:.2f},"
                        f"{4638001.0 + 2 * x:.2f},{tide + dz:.4f},{hs:.2f},8.0,0.05,{true + dz:.4f}"
                        + (f",{wave_source}" if wave_source is not None else ""))
    return rows


def tidal_frames(days=5, seed=3, c2_offset=0.0, lag_m_per_mph=0.0, hs_slope=0.0, first_hour=11, n=25):
    """Both cameras every 30 min, `n` frames from `first_hour` UTC (default 11-23 UTC), on `days` days from 1 Oct 2026, under a 12.42 h tide of
    1.2 m amplitude. Each frame's line sits at its level + noise + lag x dz/dt (m/h) + hs_slope x
    (Hs - 1) (+ c2_offset on c2): the planted water-level timing error and wave dependence."""
    rng = np.random.RandomState(seed)
    out = []
    w = 2 * np.pi / (12.42 * 3600)
    for day in range(days):
        hs = 0.6 + 0.2 * day
        for k in range(n):
            ep = int(datetime(2026, 10, 1 + day, first_hour, 0, tzinfo=timezone.utc).timestamp()) + 1800 * k
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
    for tag in (None, "marconi:hs_best/tp_s", "buoy:wvht_m/dpd_s", "unknown:wvht_m/dpd_s"):
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
    rc, log, info = res["unknown"]
    check(rc == 0 and info["wave_currency"]["frames"] == {"unknown:wvht_m/dpd_s": nframes}
          and info["wave_currency"]["all_fitted_currency"] is False and "WARNING" not in log
          and any(ln.startswith("NOTE: unknown:wvht_m/dpd_s -- a wave file its columns do not name")
                  for ln in log.splitlines()),
          "rows from a wave file of unknown kind (process_chelsea.py's ADCP record): a NOTE, no WARNING")
    rc, log, info = res["none"]
    check(rc == 0 and info["wave_currency"]["frames"] == {"unrecorded": nframes} and "WARNING" not in log,
          "contours from before wave_source: 'unrecorded', no WARNING")
    grids = [np.loadtxt(d / f"{n}_{k}.asc", skiprows=6) for n in ("none", "marconi", "buoy")
             for k in ("dem", "spread", "count")]
    check(all(np.array_equal(grids[k], grids[3 + k]) and np.array_equal(grids[k], grids[6 + k]) for k in range(3)),
          "the dem, spread and count grids are identical with and without the column")


def build(d, name, frames, extra=()):
    g = Path(d) / f"{name}.csv"
    g.write_text("\n".join(ground_rows(frames, "marconi:hs_best/tp_s")) + "\n")
    rc, log = run([HERE / "dem_from_contours.py", g, Path(d) / name, "--no-plot"] + list(extra), d)
    (Path(d) / f"{name}.log").write_text(log)
    info = json.loads((Path(d) / f"{name}_info.json").read_text()) if rc == 0 else {}
    return rc, log, info


def test_water_level_rates():
    print("dem_from_contours.water_level_rates: dz/dt from the frames either side")
    import dem_from_contours as dfc
    t0 = 1790000000
    info = {f"f{i}": {"epoch": float(t0 + 1800 * i), "tide": 0.1 * i, "hs": 1.0, "tp": 8.0, "setup": 0.0}
            for i in range(4)}
    info["c2_twin"] = dict(info["f1"])                                 # the other camera, same capture
    info["alone"] = {"epoch": float(t0 + 1800 * 3 + 4 * 3600), "tide": 0.0, "hs": 1.0, "tp": 8.0, "setup": 0.0}
    r = dfc.water_level_rates(info)
    check(abs(r["f1"] - 0.2) < 1e-9 and abs(r["f0"] - 0.2) < 1e-9 and abs(r["f3"] - 0.2) < 1e-9
          and r["c2_twin"] == r["f1"] and np.isnan(r["alone"]),
          "0.1 m per 30 min = 0.2 m/h: centred inside, one-sided at the ends, the same for both cameras of "
          "a capture, none for a frame 4 h from any other")


def test_residuals(d):
    print("dem_from_contours.py: residuals vs dz/dt and Hs (the water-level timing test)")
    d = Path(d) / "residuals"
    d.mkdir(parents=True, exist_ok=True)
    # a planted timing error: every line given an elevation 0.1 m per m/h of dz/dt too high (the line
    # 6 min older than its level)
    rc, log, info = build(d, "lag", tidal_frames(days=10, lag_m_per_mph=0.1))
    rd = info.get("diagnostics", {}).get("residuals", {}).get("by_camera", {})
    b = rd.get("both", {})
    v = b.get("vs_dzdt") or {}
    st = b.get("by_stage", {})
    check(rc == 0 and v and abs(v["slope"] - 0.1) < 0.03 and v["r"] > 0.6 and abs(b["lag_min"] - 6) < 2
          and st["flood"]["median_residual_m"] > 0.03 and st["ebb"]["median_residual_m"] < -0.03,
          f"planted +0.1 m per m/h: slope {v.get('slope')}, r {v.get('r')}, lag {b.get('lag_min')} min; flood "
          f"{st.get('flood', {}).get('median_residual_m')} m, ebb {st.get('ebb', {}).get('median_residual_m')} m")
    vd = b.get("day_offset_vs_mean_dzdt")
    check(vd is not None and vd["n"] == 20 and len(b["days"]) == 20
          and all(set(x) == {"offset_m", "mean_dzdt_m_per_h"} for x in b["days"].values()),
          "camera-days: each with its offset and the mean dz/dt of its frames, and the regression over the 20")
    # days whose frames fall mostly on the flood or mostly on the ebb (6 h of daylight a day, the
    # tide 50 min later each day): with the timing error the flood days read high
    rcs, logs, infos = build(d, "lag_short", tidal_frames(days=14, lag_m_per_mph=0.1, first_hour=13, n=13))
    vs = infos["diagnostics"]["residuals"]["by_camera"]["both"]["day_offset_vs_mean_dzdt"]
    check(rcs == 0 and vs and vs["slope"] > 0.03 and vs["r"] > 0.5,
          f"6-h days, mostly flood or mostly ebb: day offset vs mean dz/dt slope {vs['slope']:+.3f} m per m/h, "
          f"r {vs['r']:+.2f} (flood days read high)")
    check("Residuals vs tide stage and waves" in log and "lag +" in log and "WARNING" not in log,
          "printed in the DEM log (no WARNING: reporting only)")
    rc0, log0, info0 = build(d, "none", tidal_frames(days=10))
    v0 = info0["diagnostics"]["residuals"]["by_camera"]["both"]["vs_dzdt"]
    check(rc0 == 0 and abs(v0["slope"]) < 0.02 and abs(v0["r"]) < 0.3,
          f"no timing error: slope {v0['slope']:+.3f}, r {v0['r']:+.2f}")
    rch, logh, infoh = build(d, "hs", tidal_frames(days=10, hs_slope=0.1))
    vh = infoh["diagnostics"]["residuals"]["by_camera"]["both"]["vs_hs"]
    check(rch == 0 and abs(vh["slope"] - 0.1) < 0.03 and vh["r"] > 0.6,
          f"planted +0.1 m per m of Hs: vs_hs slope {vh['slope']:+.3f}, r {vh['r']:+.2f}")
    # the DEM itself does not depend on the diagnostics
    rcn, logn, infon = build(d, "lag_nodiag", tidal_frames(days=10, lag_m_per_mph=0.1), ["--no-diagnostics"])
    same = all(np.array_equal(np.loadtxt(d / f"lag_{k}.asc", skiprows=6),
                              np.loadtxt(d / f"lag_nodiag_{k}.asc", skiprows=6)) for k in ("dem", "spread", "count"))
    check(rcn == 0 and same and infon["diagnostics"] == {} and "Residuals vs" not in logn,
          "--no-diagnostics: identical grids, diagnostics {} and nothing printed")
    rcd, logd, infod = build(d, "lag_daytest", tidal_frames(days=10, lag_m_per_mph=0.1), ["--max-day-offset", "0.15"])
    check(rcd == 0 and "before it left any day out" in infod["diagnostics"]["residuals"]["on"]
          and logd.index("Residuals vs") < logd.index("Day consistency"),
          "with the day test: computed on the frames it judged, printed before its block")


def read_asc(path):
    a = np.loadtxt(path, skiprows=6)
    return np.where(a == -9999, np.nan, a)


def test_sigma(d):
    print("dem_from_contours.py: <stem>_sigma.asc, the standard error of each cell")
    d = Path(d) / "sigma"
    d.mkdir(parents=True, exist_ok=True)
    g = d / "g.csv"
    g.write_text("\n".join(ground_rows(tidal_frames(days=8), "marconi:hs_best/tp_s")) + "\n")
    rc, log = run([HERE / "dem_from_contours.py", g, d / "w", "--no-plot", "--last-days", "7",
                   "--series-dir", d / "series"], d)
    (d / "w.log").write_text(log)
    dem, spread, count, sig = (read_asc(d / f"w_{k}.asc") for k in ("dem", "spread", "count", "sigma"))
    count = np.nan_to_num(np.loadtxt(d / "w_count.asc", skiprows=6))
    want = np.where(np.isfinite(dem), 1.25 * (spread / 2) / np.sqrt(np.maximum(count, 1)), np.nan)
    check(rc == 0 and sig.shape == dem.shape and np.isfinite(sig).sum() == np.isfinite(dem).sum() > 0
          and np.allclose(sig[np.isfinite(dem)], want[np.isfinite(dem)], atol=1e-4),
          f"SE = 1.25 (spread/2)/sqrt(n) on every cell with a value ({int(np.isfinite(sig).sum())}; to the "
          f"grids' 4 decimals), NaN elsewhere")
    info = json.loads((d / "w_info.json").read_text())
    series = sorted((d / "series").glob("*_sigma.asc"))
    check(info["sigma"]["file"] == "_sigma.asc" and abs(info["sigma"]["median_m"] - float(np.nanmedian(sig))) < 1e-4
          and "not accuracy" in info["sigma"]["definition"],
          f"info.json sigma: the definition, median {info['sigma']['median_m']} m")
    check(len(series) == 1 and np.array_equal(*(np.where(np.isnan(x), -1.0, x) for x in (read_asc(series[0]), sig)))
          and "wrote " in log and "_sigma.asc" in log and "WARNING" not in log,
          "the dated series copy too; logged")
    rc2, _ = run([HERE / "asc_to_geotiff.py", d / "w_sigma.asc"], d)
    check(rc2 == 0 and (d / "w_sigma.tif").exists(), "asc_to_geotiff.py converts it (as survey_products.py does)")
    # its GeoTIFF says what it is, as the other grids' do (the description was the bare file name)
    check(rc2 == 0 and b"w_sigma: standard error of the cell median" in (d / "w_sigma.tif").read_bytes(),
          "the _sigma GeoTIFF's description names it (standard error of the cell median, m)")
    # one crossing has no scatter: with --min-points 1 a cell of one frame has a value and a spread of
    # 0 (build_grid), and its SE read 0 m -- a perfect cell
    rc1, log1 = run([HERE / "dem_from_contours.py", g, d / "one", "--no-plot", "--min-points", "1",
                     "--start-date", "2026-10-01", "--end-date", "2026-10-01"], d)
    dem1, sig1 = read_asc(d / "one_dem.asc"), read_asc(d / "one_sigma.asc")
    cnt1 = np.nan_to_num(np.loadtxt(d / "one_count.asc", skiprows=6))
    single = np.isfinite(dem1) & (cnt1 == 1)
    multi = np.isfinite(dem1) & (cnt1 >= 2)
    check(rc1 == 0 and single.sum() > 0 and multi.sum() > 0 and np.isnan(sig1[single]).all()
          and np.isfinite(sig1[multi]).all(),
          f"--min-points 1: the {int(single.sum())} cell(s) of one frame get no SE (NaN, not 0 m); the "
          f"{int(multi.sum())} of 2 or more frames get one")
    import dem_change
    check(not any("sigma" in p.name for p in dem_change.Path(d / "series").glob("dem_*_*d_dem.asc")),
          "dem_change.py's series glob (dem_*_*d_dem.asc) does not pick it up")


def test_camera_offset(d):
    print("dem_from_contours.py: the c1 - c2 offset in the cells both cameras cross")
    d = Path(d) / "cameras"
    d.mkdir(parents=True, exist_ok=True)
    rc, log, info = build(d, "off", tidal_frames(days=8, c2_offset=0.10))
    co = info.get("diagnostics", {}).get("camera_offset") or {}
    check(rc == 0 and co.get("cells", 0) >= 10 and abs(co["median_m"] + 0.10) < 0.02 and co["nmad_m"] < 0.05
          and co["p16_m"] <= co["median_m"] <= co["p84_m"],
          f"c2 planted 0.10 m high: c1 - c2 median {co.get('median_m')} m over {co.get('cells')} shared cells, "
          f"NMAD {co.get('nmad_m')} m")
    check("Camera offset     : c1 - c2 over" in log and "WARNING" not in log, "printed in the DEM log")
    rc0, log0, info0 = build(d, "same", tidal_frames(days=8))
    co0 = info0["diagnostics"]["camera_offset"]
    check(rc0 == 0 and abs(co0["median_m"]) < 0.02, f"no offset planted: median {co0['median_m']:+.3f} m")
    rc1, log1, info1 = build(d, "c1only", [f for f in tidal_frames(days=8) if f[0] == "c1"])
    check(rc1 == 0 and info1["diagnostics"]["camera_offset"] is None and "Camera offset" not in log1,
          "one camera: None, nothing printed")
    same = all(np.array_equal(np.loadtxt(d / f"off_{k}.asc", skiprows=6), np.loadtxt(d / f"same_{k}.asc", skiprows=6))
               for k in ("count",))
    check(same, "the count grid is the same whatever the offset (the cells pool both cameras, as before)")


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
        test_water_level_rates()
        test_residuals(d)
        test_sigma(d)
        test_camera_offset(d)
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

#!/usr/bin/env python3
"""
Self-Test For historical_forcing.py
=====================================
Checks that the water-level and wave forcing for the survey products does
what its report says, on data whose right answer is known.

WHY. The Chatham -> Marconi transfer sets every elevation of the 2024-25
DEMs that the ADCP does not cover, and nothing downstream can tell a wrong
transfer from a real beach change. So the fit must provably recover a
known transfer, and gaps or dates outside the records must give clear
messages instead of silent rows.

HOW. Everything is built from the REAL Signature 1000 ADCP record
(sig1000_waves_ALL.csv, in the repository):
  * a synthetic "Chatham" gauge, 6-min, with a KNOWN transfer:
        Marconi(t) = a * Chatham(t - lag) + b,  a 1.35, lag -54 min, b -0.05,
    plus 1 cm white noise and a 36 h hole;
  * a synthetic WIS hindcast with a known height and period conversion;
  * the ADCP with one week cut out, so the transfer has to fill that week
    and can be scored against the real ADCP levels that were cut.
Then: the linear fit must find a, lag and b; the harmonic transfer must
find the same thing (every constituent ratio a, Marconi ahead by 54 min,
residual scaled by a); the end-to-end run must fill the week to within
the noise, leave the hole empty and say so, and write files that
extract_elevation_contours.py and process_chelsea.py read; dates outside
the records, reversed dates, bad dates and missing files must raise
ForcingError with a message that says what to do; the NOAA download
fallback must be used only when needed (tested with a stand-in fetch).

No pytest needed. Runs in ~30 s here (about a minute on the station NUC).

Usage:
    python3 test_historical_forcing.py
    python3 test_historical_forcing.py --keep /tmp/hf_test     (keep the outputs)
    python3 test_historical_forcing.py --adcp /path/to/sig1000_waves_ALL.csv
"""

import io
import sys
import json
import shutil
import argparse
import tempfile
import subprocess
import contextlib
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import historical_forcing as hf   # noqa: E402

A_TRUE, LAG_TRUE_MIN, B_TRUE, NOISE_M = 1.35, -54, -0.05, 0.01
ADCP_OFFSET = 0.1012                        # adcp_water_level_navd88.csv = water_level + this
CUT = ("2025-02-18", "2025-02-25")          # ADCP week removed: the transfer must fill it
HOLE = ("2025-02-21 00:00", "2025-02-22 12:00")   # Chatham hole inside that week
HS_K, HS_C, TP_R = 0.833, 0.004, 1.06       # synthetic WIS: ADCP Hs = K*WIS + C, ADCP Tp = R*WIS Tp

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""), flush=True)


def ep(text):
    return hf.to_epoch([pd.Timestamp(text)])[0]


def quiet(fn, *a, **k):
    """Run fn with its progress lines captured; returns (result, printed text)."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        r = fn(*a, **k)
    return r, buf.getvalue()


def expect_error(name, fn, must_contain):
    try:
        quiet(fn)
    except hf.ForcingError as e:
        check(name, must_contain.lower() in str(e).lower(), str(e)[:150])
        return
    except Exception as e:     # any other exception is a failure: the user gets a traceback
        check(name, False, f"{type(e).__name__}: {e}")
        return
    check(name, False, "no error raised")


# --------------------------------------------------------------------------
# Fixture
# --------------------------------------------------------------------------

def build_fixture(adcp_path, d):
    """Writes the fixture files into d; returns a dict of paths and the truth series."""
    rng = np.random.RandomState(20250123)
    a = pd.read_csv(adcp_path)
    a_ep = hf.to_epoch(a["time"])
    m = a["water_level"].to_numpy(float) + ADCP_OFFSET        # Marconi truth, NAVD88

    # Chatham(s) = (Marconi(s + lag) - b) / a + noise, on the gauge's 6-min clock
    lag_s = LAG_TRUE_MIN * 60.0
    s0 = np.ceil((a_ep[0] - lag_s) / 360.0) * 360.0
    s = np.arange(s0, a_ep[-1] - lag_s + 1.0, 360.0)
    chat = (np.interp(s + lag_s, a_ep, m) - B_TRUE) / A_TRUE + rng.normal(0, NOISE_M, len(s))
    keep = ~((s >= ep(HOLE[0])) & (s < ep(HOLE[1])))
    g = pd.DataFrame({"time": pd.to_datetime(s[keep], unit="s", utc=True).strftime("%Y-%m-%d %H:%M:%S+00:00"),
                      "level": np.round(chat[keep], 4)})
    chatham = d / "chatham_synthetic.csv"
    g.to_csv(chatham, index=False)

    # The ADCP with one week cut out, in both of its files
    cut = (a_ep >= ep(CUT[0])) & (a_ep < ep(CUT[1]))
    adcp_cut = d / "sig1000_cut.csv"
    a[~cut].to_csv(adcp_cut, index=False)
    navd = d / "adcp_navd88_cut.csv"
    pd.DataFrame({"time": pd.to_datetime(a_ep[~cut], unit="s", utc=True).strftime("%Y-%m-%dT%H:%M:%SZ"),
                  "water_level_navd88": np.round(m[~cut], 4),
                  "water_level_adcp": np.round(a["water_level"].to_numpy(float)[~cut], 4)}).to_csv(navd, index=False)

    # WIS hindcast with a known conversion (and a little noise), over the whole ADCP span
    hs = a["wh_4061"].to_numpy(float)
    tp = a["wp_peak"].to_numpy(float)
    wis = d / "wis_synthetic.csv"
    pd.DataFrame({"datetime": pd.to_datetime(a_ep, unit="s", utc=True).strftime("%Y-%m-%d %H:%M:%S"),
                  "lat": 41.916672, "lon": -69.75,
                  "waveTp": np.round(tp / TP_R, 3),
                  "waveHs": np.round((hs - HS_C) / HS_K * (1 + rng.normal(0, 0.02, len(hs))), 3),
                  "waveTm": np.round(tp / TP_R * 0.8, 3),
                  "waveMeanDirection": a["wvdir"].to_numpy(float)}).to_csv(wis, index=False)
    return {"adcp_full": Path(adcp_path), "adcp": adcp_cut, "adcp_navd88": navd, "chatham": chatham,
            "wis": wis, "truth": (a_ep, m, hs, tp), "gauge": (s[keep], chat[keep])}


# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------

def test_linear_recovery(fx):
    a_ep, m = fx["truth"][0], fx["truth"][1]
    lt, _ = quiet(lambda: hf.LinearTransfer().fit(a_ep, m, fx["gauge"]))
    check("linear: a recovered", abs(lt.a - A_TRUE) < 0.01, f"{lt.a:.4f} vs {A_TRUE}")
    check("linear: lag recovered", lt.lag_min == LAG_TRUE_MIN, f"{lt.lag_min:+.0f} vs {LAG_TRUE_MIN:+d} min")
    check("linear: b recovered", abs(lt.b - B_TRUE) < 0.01, f"{lt.b:+.4f} vs {B_TRUE:+.2f}")
    check("linear: fit RMS at the noise level", lt.rms_fit < 2 * A_TRUE * NOISE_M,
          f"{lt.rms_fit:.4f} m (noise x a = {A_TRUE * NOISE_M:.4f})")


def test_harmonic_recovery(fx):
    """A pure scale + delay is a special case of the harmonic transfer: every constituent's
    ratio must be a, every phase a 54-min lead, and the residual factor a at the same lag."""
    a_ep, m = fx["truth"][0], fx["truth"][1]
    ht, _ = quiet(lambda: hf.HarmonicTransfer().fit(a_ep, m, fx["gauge"]))
    tab = {r["constituent"]: r for r in ht.constituent_table()}
    for c, tol_ratio, tol_min in (("M2", 0.03, 3.0), ("N2", 0.04, 6.0), ("S2", 0.04, 6.0),
                                  ("K1", 0.05, 10.0), ("O1", 0.05, 10.0)):
        r = tab.get(c)
        if r is None:
            check(f"harmonic: {c} fitted", False, "missing")
            continue
        check(f"harmonic: {c} ratio and lead", abs(r["ratio"] / A_TRUE - 1) < tol_ratio
              and abs(r["marconi_leads_min"] + LAG_TRUE_MIN) < tol_min,
              f"ratio {r['ratio']:.3f}, leads {r['marconi_leads_min']:+.1f} min")
    check("harmonic: residual factor ~ a", abs(ht.beta / A_TRUE - 1) < 0.08, f"beta {ht.beta:.3f}")
    check("harmonic: residual lag ~ lag", abs(ht.lag_r_min - LAG_TRUE_MIN) <= 12, f"{ht.lag_r_min:+.0f} min")
    # Prediction from the gauge alone over a week it was not fitted on
    tr = a_ep < ep("2025-02-10")
    ht2, _ = quiet(lambda: hf.HarmonicTransfer().fit(a_ep[tr], m[tr], fx["gauge"]))
    te = (a_ep >= ep("2025-02-24")) & (a_ep < ep("2025-03-03"))
    p, _ = quiet(lambda: ht2.predict(a_ep[te], fx["gauge"]))
    rms = float(np.sqrt(np.nanmean((p - m[te]) ** 2)))
    check("harmonic: out-of-window prediction near the noise", rms < 0.04, f"RMS {rms:.4f} m")


def test_end_to_end(fx, out):
    """The cut week is filled by the transfer and the WIS conversion; the hole stays empty."""
    res, log = quiet(lambda: hf.build_forcing("2025-02-17", "2025-02-25", out, adcp=fx["adcp"],
                                              adcp_navd88=fx["adcp_navd88"], chatham=fx["chatham"],
                                              wis=fx["wis"], ndbc=None, download=False, plot=True))
    wl = pd.read_csv(res["water_level_csv"])
    wv = pd.read_csv(res["waves_csv"])
    js = json.loads((Path(out) / "forcing.json").read_text())
    check("e2e: water_level.csv columns", list(wl.columns) == ["time", "water_level_navd88", "source", "sigma_m"],
          str(list(wl.columns)))
    check("e2e: both sources used", set(wl["source"]) == {"adcp", "chatham_transfer"}, str(set(wl["source"])))
    w_ep = hf.to_epoch(wl["time"])
    tr = (wl["source"] == "chatham_transfer").to_numpy()
    check("e2e: transfer only inside the cut week",
          bool(tr.any()) and w_ep[tr].min() > ep(CUT[0]) - 3600 and w_ep[tr].max() < ep(CUT[1]) + 3600,
          f"{hf.day_str(w_ep[tr].min())} .. {hf.day_str(w_ep[tr].max())}")
    # score the transferred rows against the real ADCP levels that were cut (on the hour)
    a_ep, m = fx["truth"][0], fx["truth"][1]
    on_hour = tr & (w_ep % 3600 == 0)
    truth = np.interp(w_ep[on_hour], a_ep, m)
    err = wl["water_level_navd88"].to_numpy(float)[on_hour] - truth
    rms = float(np.sqrt(np.mean(err ** 2)))
    check("e2e: transferred levels match the cut ADCP", rms < 0.04,
          f"RMS {rms:.4f} m over {on_hour.sum()} h, method {js['water_level']['method']}")
    # Marconi at t needs Chatham at t - lag (54 min later), so the hole in Marconi time is the
    # gauge's hole moved 54 min earlier; 30 min margins allow the harmonic's own residual lag.
    h0 = ep(HOLE[0]) + LAG_TRUE_MIN * 60.0 + 1800
    h1 = ep(HOLE[1]) + LAG_TRUE_MIN * 60.0 - 1800
    hole = (w_ep > h0) & (w_ep < h1)
    check("e2e: no rows inside the Chatham hole", not hole.any(), f"{hole.sum()} rows")
    spans = [[ep(x.strip()) for x in sp.replace(" UTC", "").split("..")]
             for sp in js["water_level"]["uncovered_spans"]]
    check("e2e: the hole is reported as not covered", any(lo <= h0 and hi >= h1 for lo, hi in spans),
          "; ".join(js["water_level"]["uncovered_spans"])[:120])
    check("e2e: warning printed for the hole", "WARNING" in log and "no water level for" in log)
    sig = js["water_level"]["transfer_sigma_m"]
    check("e2e: sigma_m of transferred rows = extrapolation RMS",
          np.allclose(wl["sigma_m"][tr], round(sig, 4)) and sig == max(
              v for v in (js["water_level"]["extrapolation_rms_m"], js["water_level"]["extrapolation_rms_daytime_m"])
              if v is not None), f"{sig}")
    for k in ("method", "params", "cv_rms_m", "extrapolation_rms_m", "sources_used"):
        check(f"e2e: forcing.json water_level.{k}", k in js["water_level"])
    check("e2e: both transfers scored", set(js["water_level"]["methods"]) == {"linear", "harmonic"})
    # waves
    check("e2e: waves.csv columns", {"time_utc", "epoch", "wvht_m", "dpd_s", "source", "hs_sigma_m"} <= set(wv.columns),
          str(list(wv.columns)))
    cut = ((wv["epoch"] >= ep(CUT[0])) & (wv["epoch"] < ep(CUT[1]))).to_numpy()
    check("e2e: cut week waves from WIS, with an error", cut.any() and wv["source"][cut].str.contains("wis").all()
          and (wv["hs_sigma_m"][cut] > 0).all(), str(set(wv["source"][cut])))
    check("e2e: ADCP waves elsewhere, error 0", (wv["source"][~cut] == "adcp").all()
          and (wv["hs_sigma_m"][~cut] == 0).all())
    hs_true = np.interp(wv["epoch"][cut].to_numpy(float), a_ep, fx["truth"][2])
    rms_hs = float(np.sqrt(np.mean((wv["wvht_m"][cut].to_numpy(float) - hs_true) ** 2)))
    check("e2e: converted WIS Hs matches the cut ADCP Hs", rms_hs < 0.08, f"RMS {rms_hs:.3f} m")
    tp_true = np.interp(wv["epoch"][cut].to_numpy(float), a_ep, fx["truth"][3])
    rel = float(np.median(wv["dpd_s"][cut].to_numpy(float) / tp_true))
    check("e2e: WIS Tp mapped to the ADCP's", abs(rel - 1) < 0.03, f"median ratio {rel:.3f}")
    check("e2e: Hs currency stated", "ADCP" in js["waves"]["hs_currency"])
    check("e2e: the cut week is described as partly covered", "% covered by the ADCP"
          in js["water_level"].get("period_vs_adcp", ""), js["water_level"].get("period_vs_adcp"))
    check("e2e: report and figure written", (Path(out) / "forcing_report.txt").exists()
          and (Path(out) / "forcing.png").stat().st_size > 20000)
    rep = (Path(out) / "forcing_report.txt").read_text()
    check("e2e: report names the datum caveat and the C currency",
          "ASSUMES Marconi's mean level equals Chatham's" in rep and "currency" in rep)
    check("e2e: report and forcing.json state C's still-water reference (GNSS-R, surf zone)",
          "STILL-WATER REFERENCE OF C" in rep and "c_fitted_with" in js["waves"]["still_water_reference"])

    # downstream readers
    import extract_elevation_contours as eec
    (t, lv, _), _ = quiet(lambda: eec.load_tide_model(res["water_level_csv"], time_col="time",
                                                      level_col="water_level_navd88"))
    check("readers: extract_elevation_contours reads water_level.csv", len(lv) == len(wl)
          and np.allclose(np.sort(t), np.sort(w_ep)))
    (t, lv, _), _ = quiet(lambda: eec.load_tide_model(res["water_level_csv"]))
    check("readers: ... also with auto-detected columns", np.allclose(lv, wl["water_level_navd88"].to_numpy(float)[
        np.argsort(w_ep, kind="stable")]))
    eh, h, et, tp = eec.load_wave_archive(res["waves_csv"])
    check("readers: extract_elevation_contours reads waves.csv", len(h) == len(wv) and len(tp) == len(wv))
    try:
        import process_chelsea
        e2, l2 = process_chelsea.load_water_level(res["water_level_csv"])
        check("readers: process_chelsea.load_water_level", len(l2) == len(wl))
        v = process_chelsea.level_at(e2, l2, ep("2025-02-20 15:30"), 3600)
        check("readers: process_chelsea.level_at inside the transfer", v is not None, str(v))
    except ImportError as e:
        print(f"SKIP  process_chelsea reader ({e})")
    return js


def test_adcp_only(fx, out):
    res, _ = quiet(lambda: hf.build_forcing("2025-01-16", "2025-01-24", out, adcp=fx["adcp"],
                                            adcp_navd88=fx["adcp_navd88"], chatham=fx["chatham"],
                                            wis=fx["wis"], download=False, plot=False))
    wl = pd.read_csv(res["water_level_csv"])
    wv = pd.read_csv(res["waves_csv"])
    check("adcp-only: every level from the ADCP", (wl["source"] == "adcp").all(), str(set(wl["source"])))
    a_ep, m = fx["truth"][0], fx["truth"][1]
    check("adcp-only: levels equal the NAVD88 file", np.allclose(
        wl["water_level_navd88"], np.interp(hf.to_epoch(wl["time"]), a_ep, m), atol=1e-4))
    check("adcp-only: every wave hour from the ADCP", (wv["source"] == "adcp").all())
    check("adcp-only: period described as inside",
          "inside" in res["water_level"].get("period_vs_adcp", ""), res["water_level"].get("period_vs_adcp"))


def test_partial_and_errors(fx, out):
    # past the end of both records: rows up to the end, the rest reported
    res, log = quiet(lambda: hf.build_forcing("2025-03-09", "2025-03-12", Path(out) / "partial", adcp=fx["adcp"],
                                              adcp_navd88=fx["adcp_navd88"], chatham=fx["chatham"], wis=fx["wis"],
                                              download=False, plot=False, min_wave_cover=0))
    unc = " ".join(res["water_level"]["uncovered_spans"])
    check("partial: hours past the records reported", "2025-03-12" in unc and "WARNING" in log, unc[:100])
    check("partial: --no-download said so", "--no-download" in log)
    kw = dict(adcp=fx["adcp"], adcp_navd88=fx["adcp_navd88"], chatham=fx["chatham"], wis=fx["wis"],
              download=False, plot=False)
    # with the setup on, waves covering under half the daytime hours FAIL (after writing the files):
    # an hour of detection would otherwise end with most frames dropped for want of a setup
    expect_error("error: setup on, too few wave hours",
                 lambda: hf.build_forcing("2025-03-09", "2025-03-12", Path(out) / "few_waves", **kw),
                 "waves cover only")
    check("error: too few wave hours still writes the files",
          (Path(out) / "few_waves" / "waves.csv").exists() and (Path(out) / "few_waves" / "forcing.json").exists())
    res0, _ = quiet(lambda: hf.build_forcing("2025-03-09", "2025-03-12", Path(out) / "few_waves_c0",
                                             **dict(kw, setup_coef=0.0)))
    check("no setup (C = 0): few wave hours are not an error", res0 is not None)
    ep = np.arange(hf.parse_day("2025-03-01", "x"), hf.parse_day("2025-03-03", "x"), 3600.0)
    n_cov, n_day = hf.wave_cover(pd.DataFrame({"ep": ep[:30]}), ep[0], ep[0] + 2 * 86400.0, 3600.0)
    check("wave_cover: daytime hours of day 1 covered, day 2 not", (n_cov, n_day) == (5, 10), (n_cov, n_day))
    expect_error("error: period outside every record",
                 lambda: hf.build_forcing("2023-06-01", "2023-06-03", Path(out) / "e1", **kw),
                 "outside every water-level record")
    expect_error("error: end before start",
                 lambda: hf.build_forcing("2025-02-10", "2025-02-01", Path(out) / "e2", **kw), "before")
    expect_error("error: not a date",
                 lambda: hf.build_forcing("2025-13-45", "2025-02-01", Path(out) / "e3", **kw), "not a date")
    expect_error("error: missing file",
                 lambda: hf.build_forcing("2025-02-01", "2025-02-02", Path(out) / "e4",
                                          **dict(kw, wis=str(Path(out) / "nope.csv"))), "file not found")
    expect_error("error: wrong columns",
                 lambda: hf.build_forcing("2025-02-01", "2025-02-02", Path(out) / "e5",
                                          **dict(kw, chatham=fx["wis"])), "columns")
    # the command line: exit code 2 and one ERROR line, no traceback
    p = subprocess.run([sys.executable, str(HERE / "historical_forcing.py"), "--start", "2023-06-01",
                        "--end", "2023-06-02", "--output-dir", str(Path(out) / "cli"), "--adcp", str(fx["adcp"]),
                        "--adcp-navd88", str(fx["adcp_navd88"]), "--chatham", str(fx["chatham"]),
                        "--no-download", "--no-plot"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=300)
    check("cli: exit 2 with an ERROR line", p.returncode == 2 and p.stderr.startswith("ERROR:")
          and "Traceback" not in p.stderr, f"rc {p.returncode}: {p.stderr.strip()[:120]}")


def test_still_water_share(out):
    """C's still-water reference: a synthetic GNSS-R record holding a KNOWN share (0.5) of the setup
    must give that share back; without the files the share is said to be not measured; the
    report text names the reference for the station's C and calls another C of unknown origin."""
    from datetime import datetime, timezone
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.RandomState(3)
    share, coef = 0.5, hf.SETUP_COEF
    t0 = datetime(2026, 8, 1, tzinfo=timezone.utc).timestamp()
    w = 2 * np.pi / (12.4206 * 3600)

    def gauge(t):
        return 0.55 * np.cos(w * t) + 0.08 * np.cos(2 * w * t + 0.3) + 0.02 * np.sin(t / 86400 / 3)
    g_ep = np.arange(t0, t0 + 40 * 86400, 360.0)
    pd.DataFrame({"epoch": g_ep.astype(int), "level_navd88": np.round(gauge(g_ep), 4)}).to_csv(
        out / "gauge.csv", index=False)
    h_ep = np.arange(t0, t0 + 40 * 86400, 3600.0)
    hs = 0.6 + 0.9 * (1 + np.sin(h_ep / 86400 / 2.3)) + rng.uniform(0, 0.3, len(h_ep))
    tp = 7 + 2.5 * (1 + np.sin(h_ep / 86400 / 3.1 + 1))
    pd.DataFrame({"epoch": h_ep.astype(int), "hs_best": np.round(hs, 2), "tp_s": np.round(tp, 1)}).to_csv(
        out / "waves.csv", index=False)
    s_ep = np.arange(t0, t0 + 40 * 86400, 900.0)
    x = np.sqrt(np.interp(s_ep, h_ep, hs) * hf.G * np.interp(s_ep, h_ep, tp) ** 2 / (2 * np.pi))
    lv = 1.24 * gauge(s_ep + 48 * 60) - 0.10 + share * coef * x + rng.normal(0, 0.02, len(s_ep))
    with open(out / "spline.txt", "w") as f:
        f.write("% synthetic gnssrefl spline\n")
        for e, v in zip(s_ep, lv):
            d = pd.Timestamp(e, unit="s")
            f.write(f"60000.0 10.0 {d.year} {d.month} {d.day} {d.hour} {d.minute} {d.second} {v:.4f}\n")
    r, _ = quiet(lambda: hf.gnssr_setup_share(coef, out / "spline.txt", out / "gauge.csv", out / "waves.csv"))
    check("still water: GNSS-R share of the setup recovered", r.get("quantified") and abs(r["share"] - share) < 0.05,
          f"share {r.get('share')}, typical {r.get('typical_m')} m")
    r2 = hf.gnssr_setup_share(coef, out / "missing.txt", out / "gauge.csv", out / "waves.csv")
    check("still water: not measured without the files, and said why",
          not r2["quantified"] and "missing.txt" in r2["note"], r2["note"][:100])
    lines, _ = hf.currency_note({}, coef, r)
    text = "\n".join(lines)
    check("still water: the report names GNSS-R, the surf zone and the measured share",
          "GNSS-R" in text and "surf zone" in text and "share x (its setup - the mean setup" in text
          and "0.02 m HIGH in average waves" in text and "datum chain" in text and "lack" not in text
          and f"{r['share']:.2f}" in text)
    lines, _ = hf.currency_note({}, 0.05, None)
    check("still water: another C is of unknown origin", "not known here" in "\n".join(lines))


def test_download_fallback(fx, out):
    """No Chatham file anywhere: the record is fetched (stand-in fetch) only because the period needs it."""
    import compare_gnssr_to_gauge as cg
    calls = []
    s, lv = fx["gauge"]

    def fake_fetch(station, start, end):
        calls.append((station, start, end))
        return pd.DataFrame({"time": pd.to_datetime(s, unit="s", utc=True), "level": lv})

    def failing_fetch(station, start, end):
        calls.append((station, start, end))
        sys.exit("NOAA API returned no data")      # what fetch_gauge does on a refusal

    real_fetch, real_dirs = cg.fetch_gauge, hf.SEARCH_DIRS
    empty = Path(out) / "empty"
    empty.mkdir(parents=True, exist_ok=True)
    try:
        hf.SEARCH_DIRS = [empty]
        cg.fetch_gauge = fake_fetch
        res, log = quiet(lambda: hf.build_forcing("2025-01-16", "2025-01-18", Path(out) / "dl0", adcp=fx["adcp"],
                                                  adcp_navd88=fx["adcp_navd88"], plot=False))
        check("download: not used when the ADCP covers the period", not calls and
              not res["water_level"].get("chatham_downloaded"))
        # the water level is the subject here: the cut week also has no waves (min_wave_cover=0)
        res, log = quiet(lambda: hf.build_forcing("2025-02-19", "2025-02-20", Path(out) / "dl1", adcp=fx["adcp"],
                                                  adcp_navd88=fx["adcp_navd88"], plot=False, min_wave_cover=0))
        saved = Path(out) / "dl1" / hf.DOWNLOADED_CHATHAM
        check("download: fetched when the period needs Chatham", len(calls) == 1 and saved.exists()
              and res["water_level"]["chatham_downloaded"]
              and res["water_level"]["inputs"]["chatham"] == str(saved), log.split("Chatham download")[-1][:100])
        check("download: covers the fit window and the target +/- 30 d",
              calls and calls[0][1] <= pd.Timestamp("2024-12-10", tz="UTC")
              and calls[0][2] >= pd.Timestamp("2025-03-10", tz="UTC"), str(calls[0][1:]) if calls else "")
        check("download: transfer filled the cut week",
              res["water_level"]["sources_used"].get("chatham_transfer", 0) > 0, str(res["water_level"]["sources_used"]))
        cg.fetch_gauge = failing_fetch
        expect_error("download: refusal -> clear error outside the ADCP",
                     lambda: hf.build_forcing("2024-10-20", "2024-10-21", Path(out) / "dl2", adcp=fx["adcp"],
                                              adcp_navd88=fx["adcp_navd88"], plot=False),
                     "outside every water-level record")
    finally:
        cg.fetch_gauge, hf.SEARCH_DIRS = real_fetch, real_dirs


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--adcp", default=str(HERE / "sig1000_waves_ALL.csv"), help="the real ADCP record")
    ap.add_argument("--keep", default=None, help="write the fixture and outputs here and keep them")
    args = ap.parse_args()
    if not Path(args.adcp).exists():
        print(f"ERROR: ADCP record {args.adcp} not found (give --adcp)")
        return 2
    work = Path(args.keep) if args.keep else Path(tempfile.mkdtemp(prefix="hf_test_"))
    work.mkdir(parents=True, exist_ok=True)
    t = pd.Timestamp.now()
    print(f"work folder       : {work}")
    print("fixture           : synthetic Chatham (a 1.35, lag -54 min, b -0.05, 1 cm noise, 36 h hole), "
          "synthetic WIS, ADCP with 2025-02-18..24 cut")
    (work / "fixture").mkdir(exist_ok=True)
    fx = build_fixture(args.adcp, work / "fixture")
    test_linear_recovery(fx)
    test_harmonic_recovery(fx)
    test_end_to_end(fx, work / "e2e")
    test_adcp_only(fx, work / "adcp_only")
    test_partial_and_errors(fx, work / "errors")
    test_download_fallback(fx, work / "download")
    test_still_water_share(work / "still_water")
    n_fail = sum(1 for _, ok in RESULTS if not ok)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed in "
          f"{(pd.Timestamp.now() - t).total_seconds():.0f} s" + (f"; {n_fail} FAILED" if n_fail else ""))
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())

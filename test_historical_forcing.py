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

STATION FILES. Nothing here depends on what the computer holds: the run
works in a temporary folder that is also its current folder, and every
build is cut off from the station's files for the whole run (the search
folders point at an empty folder of the test's own, and the GNSS-R spline,
gauge and waves archive of the setup-share check at files that do not
exist), so it runs on the fixture alone and gives the same result on the
station NUC as in a clean checkout. The search itself (an input not given
is found by its default name in the first search folder holding one) is
tested on folders of the test's own, and a last check confirms that no
build named a file outside the test's folder. The command line is run
twice: once through main() with the search cut off as above, and once as
the script itself on an error raised before any input is looked for.

No pytest needed. Runs in ~40 s here (about two minutes on the station NUC).

Usage:
    python3 test_historical_forcing.py
    python3 test_historical_forcing.py --keep /tmp/hf_test     (keep the outputs, in a new folder inside it)
    python3 test_historical_forcing.py --adcp /path/to/sig1000_waves_ALL.csv
"""

import io
import os
import sys
import json
import shutil
import argparse
import tempfile
import subprocess
import traceback
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
# Set by main() for the whole run: the empty folder the station search is pointed at, and the
# build_forcing keywords that point the setup-share check's GNSS-R spline, gauge and waves archive
# at files that do not exist there (their defaults are the station's live files).
NO_STATION = None
ISOLATED = {}


def build(start, end, out, **kw):
    """hf.build_forcing on the fixture alone: never the station's GNSS-R, gauge or waves archive."""
    return hf.build_forcing(start, end, out, **dict(ISOLATED, **kw))


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
    # the printed formula has the sign of the model: lag -54 min = Chatham(t + 54 min), Marconi leads
    d = lt.describe()
    check("linear: describe() prints the shift the model uses",
          f"Chatham(t {-LAG_TRUE_MIN:+d} min)" in d and f"Marconi leads by {-LAG_TRUE_MIN:d} min" in d, d)


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
    res, log = quiet(lambda: build("2025-02-17", "2025-02-25", out, adcp=fx["adcp"],
                                              adcp_navd88=fx["adcp_navd88"], chatham=fx["chatham"],
                                              wis=fx["wis"], ndbc=None, download=False, plot=True))
    wl = pd.read_csv(res["water_level_csv"])
    wv = pd.read_csv(res["waves_csv"])
    js = json.loads((Path(out) / "forcing.json").read_text())
    check("e2e: water_level.csv columns", list(wl.columns) == ["time", "water_level_navd88", "source", "sigma_m"],
          str(list(wl.columns)))
    check("e2e: both sources used", set(wl["source"]) == {"adcp", "chatham_transfer"}, str(sorted(set(wl["source"]))))
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
    # exactly the fixture's WIS: any other source here (a station NDBC file found by the search) is a leak
    check("e2e: cut week waves from WIS, with an error", cut.any() and set(wv["source"][cut]) == {"hs=wis;tp=wis"}
          and (wv["hs_sigma_m"][cut] > 0).all(), str(sorted(set(wv["source"][cut]))))
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
    # the share needs the station's GNSS-R record, which this run never reads (test_still_water_share
    # measures one on a synthetic record): 'not measured', naming the test's own missing files
    swr = js["waves"]["still_water_reference"]
    check("e2e: the GNSS-R share is not measured from the station's files (none given here)",
          swr.get("quantified") is False and "Not measured here: needs " in rep
          and swr.get("inputs") and all(Path(x).parent == NO_STATION for x in swr["inputs"]))

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


def test_c0_same_currency(fx, out):
    """A C = 0 build converts the waves exactly as the station's C does (review r1: with no setup to judge by, the
    conversions fell to the first listed, Hs linear and Tp raw, another wave currency than the C build's)."""
    kw = dict(adcp=fx["adcp"], adcp_navd88=fx["adcp_navd88"], chatham=fx["chatham"], wis=fx["wis"], ndbc=None,
              download=False, plot=False)
    res_c, _ = quiet(lambda: build("2025-02-17", "2025-02-25", Path(out) / "c", **kw))
    res_0, _ = quiet(lambda: build("2025-02-17", "2025-02-25", Path(out) / "c0", **dict(kw, setup_coef=0.0)))
    wc, w0 = pd.read_csv(res_c["waves_csv"]), pd.read_csv(res_0["waves_csv"])
    jc = json.loads((Path(out) / "c" / "forcing.json").read_text())["waves"]
    j0 = json.loads((Path(out) / "c0" / "forcing.json").read_text())["waves"]
    conv = (wc["source"] != "adcp").to_numpy()
    cols = [c for c in wc.columns if c not in ("time", "source")]
    same = (len(wc) == len(w0) and list(wc["source"]) == list(w0["source"]) and list(wc.columns) == list(w0.columns)
            and all(np.allclose(pd.to_numeric(wc[c], errors="coerce").to_numpy(float),
                                pd.to_numeric(w0[c], errors="coerce").to_numpy(float), equal_nan=True) for c in cols))
    check("C = 0: the converted waves are those of the station's C (the same values every hour)",
          bool(conv.any()) and same and j0.get("conversions_judged_with_c") == hf.SETUP_COEF
          and j0.get("setup_coef") == 0.0 and jc.get("conversions_judged_with_c") == hf.SETUP_COEF,
          f"{int(conv.sum())} converted hours; columns {cols}; judged with C = {j0.get('conversions_judged_with_c')}")


def test_adcp_only(fx, out):
    res, _ = quiet(lambda: build("2025-01-16", "2025-01-24", out, adcp=fx["adcp"],
                                            adcp_navd88=fx["adcp_navd88"], chatham=fx["chatham"],
                                            wis=fx["wis"], download=False, plot=False))
    wl = pd.read_csv(res["water_level_csv"])
    wv = pd.read_csv(res["waves_csv"])
    check("adcp-only: every level from the ADCP", (wl["source"] == "adcp").all(), str(sorted(set(wl["source"]))))
    a_ep, m = fx["truth"][0], fx["truth"][1]
    check("adcp-only: levels equal the NAVD88 file", np.allclose(
        wl["water_level_navd88"], np.interp(hf.to_epoch(wl["time"]), a_ep, m), atol=1e-4))
    check("adcp-only: every wave hour from the ADCP", (wv["source"] == "adcp").all())
    check("adcp-only: period described as inside",
          "inside" in res["water_level"].get("period_vs_adcp", ""), res["water_level"].get("period_vs_adcp"))


def test_partial_and_errors(fx, out):
    # past the end of both records: rows up to the end, the rest reported
    res, log = quiet(lambda: build("2025-03-09", "2025-03-12", Path(out) / "partial", adcp=fx["adcp"],
                                              adcp_navd88=fx["adcp_navd88"], chatham=fx["chatham"], wis=fx["wis"],
                                              download=False, plot=False, min_wave_cover=0))
    unc = " ".join(res["water_level"]["uncovered_spans"])
    check("partial: hours past the records reported", "2025-03-12" in unc and "no water level for" in log,
          unc[:100])
    wunc = " ".join(res["waves"]["uncovered_spans"])
    check("partial: hours past the wave records reported too", "2025-03-12" in wunc
          and "no waves for" in log, wunc[:100])
    check("partial: --no-download said so", "--no-download" in log)
    kw = dict(adcp=fx["adcp"], adcp_navd88=fx["adcp_navd88"], chatham=fx["chatham"], wis=fx["wis"],
              download=False, plot=False)
    # with the setup on, waves covering under half the daytime hours FAIL (after writing the files):
    # an hour of detection would otherwise end with most frames dropped for want of a setup
    expect_error("error: setup on, too few wave hours",
                 lambda: build("2025-03-09", "2025-03-12", Path(out) / "few_waves", **kw),
                 "waves cover only")
    fw = Path(out) / "few_waves"
    fj = json.loads((fw / "forcing.json").read_text()) if (fw / "forcing.json").exists() else {}
    check("error: too few wave hours still writes the files (their waves really are short)",
          (fw / "waves.csv").exists() and fj.get("waves", {}).get("uncovered_spans"),
          "; ".join(fj.get("waves", {}).get("uncovered_spans", []))[:100])
    try:
        res0, _ = quiet(lambda: build("2025-03-09", "2025-03-12", Path(out) / "few_waves_c0",
                                      **dict(kw, setup_coef=0.0)))
        ok0, why0 = bool(res0["waves"]["uncovered_spans"]), "; ".join(res0["waves"]["uncovered_spans"])[:100]
    except hf.ForcingError as e:
        ok0, why0 = False, f"ForcingError: {str(e)[:150]}"
    check("no setup (C = 0): few wave hours are not an error", ok0, why0)
    ep = np.arange(hf.parse_day("2025-03-01", "x"), hf.parse_day("2025-03-03", "x"), 3600.0)
    n_cov, n_day = hf.wave_cover(pd.DataFrame({"ep": ep[:30]}), ep[0], ep[0] + 2 * 86400.0, 3600.0)
    check("wave_cover: daytime hours of day 1 covered, day 2 not", (n_cov, n_day) == (5, 10), (n_cov, n_day))
    expect_error("error: period outside every record",
                 lambda: build("2023-06-01", "2023-06-03", Path(out) / "e1", **kw),
                 "outside every water-level record")
    expect_error("error: end before start",
                 lambda: build("2025-02-10", "2025-02-01", Path(out) / "e2", **kw), "before")
    expect_error("error: not a date",
                 lambda: build("2025-13-45", "2025-02-01", Path(out) / "e3", **kw), "not a date")
    expect_error("error: missing file",
                 lambda: build("2025-02-01", "2025-02-02", Path(out) / "e4",
                                          **dict(kw, wis=str(Path(out) / "nope.csv"))), "file not found")
    expect_error("error: wrong columns",
                 lambda: build("2025-02-01", "2025-02-02", Path(out) / "e5",
                                          **dict(kw, chatham=fx["wis"])), "columns")
    # the command line: exit code 2 and one ERROR line, no traceback. A child process does not see this
    # process's search folders, so it runs historical_forcing.main() with the station search pointed at
    # the same empty folder (and the setup-share files at missing ones), as every build here
    child = ("import sys; from pathlib import Path; sys.path.insert(0, sys.argv[1]); "
             "import historical_forcing as hf; hf.SEARCH_DIRS = [Path(sys.argv[2])]; "
             "sys.argv = ['historical_forcing.py'] + sys.argv[3:]; sys.exit(hf.main())")
    p = subprocess.run([sys.executable, "-c", child, str(HERE), str(NO_STATION), "--start", "2023-06-01",
                        "--end", "2023-06-02", "--output-dir", str(Path(out) / "cli"), "--adcp", str(fx["adcp"]),
                        "--adcp-navd88", str(fx["adcp_navd88"]), "--chatham", str(fx["chatham"]),
                        "--gnssr-spline", ISOLATED["gnssr_spline"], "--gauge-archive", ISOLATED["gauge_archive"],
                        "--waves-archive", ISOLATED["waves_archive"], "--no-download", "--no-plot"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=300)
    searched = [ln for ln in p.stdout.splitlines() if ln.strip().startswith(("wis ", "ndbc "))]
    check("cli: exit 2 with an ERROR line", p.returncode == 2 and p.stderr.startswith("ERROR:")
          and "Traceback" not in p.stderr and len(searched) == 2 and all("not found" in ln for ln in searched),
          f"rc {p.returncode}: {p.stderr.strip()[:120]}")
    # ... and the script itself, through its own entry point (`sys.exit(main())`): survey_products' forcing step
    # stops on its exit code alone after 'too few wave hours', whose files exist. An error build_forcing raises
    # before it makes the output folder or looks for any input, so no station file can be reached: the folder
    # not made is part of the check (should the order change, this check would no longer be cut off)
    real_out = Path(out) / "cli_script"
    p = subprocess.run([sys.executable, str(HERE / "historical_forcing.py"), "--start", "2025-02-10", "--end",
                        "2025-02-01", "--output-dir", str(real_out), "--no-download", "--no-plot"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=300)
    check("cli: historical_forcing.py itself exits 2 with an ERROR line, before any input is looked for",
          p.returncode == 2 and p.stderr.startswith("ERROR:") and "before" in p.stderr
          and "Traceback" not in p.stderr and not real_out.exists(),
          f"rc {p.returncode}: {p.stderr.strip()[:120]}; output folder made: {real_out.exists()}")


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
    # what a build given none of the three files reads: the station's (marconi_water_level.py's
    # spline and gauge archive, archive/waves_marconi.csv). Path logic only: nothing is read here
    import marconi_water_level as mwl
    check("still water: by default the station's GNSS-R spline, gauge archive and waves archive",
          hf.setup_share_inputs() == [str(mwl.GNSSR_SPLINE), str(mwl.GAUGE_CSV),
                                      str(hf.HERE / "archive" / "waves_marconi.csv")]
          and hf.setup_share_inputs("a", "b", "c") == ["a", "b", "c"])


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
        res, log = quiet(lambda: build("2025-01-16", "2025-01-18", Path(out) / "dl0", adcp=fx["adcp"],
                                                  adcp_navd88=fx["adcp_navd88"], plot=False))
        check("download: not used when the ADCP covers the period", not calls and
              not res["water_level"].get("chatham_downloaded"))
        # the water level is the subject here: the cut week also has no waves (min_wave_cover=0)
        res, log = quiet(lambda: build("2025-02-19", "2025-02-20", Path(out) / "dl1", adcp=fx["adcp"],
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
                     lambda: build("2024-10-20", "2024-10-21", Path(out) / "dl2", adcp=fx["adcp"],
                                              adcp_navd88=fx["adcp_navd88"], plot=False),
                     "outside every water-level record")
    finally:
        cg.fetch_gauge, hf.SEARCH_DIRS = real_fetch, real_dirs


def test_station_search(fx, out):
    """The station behaviour every other test switches off: an input not given is looked for by its
    default names in the search folders, the first folder holding one wins (whichever of its names),
    a file given wins over the search, and an input found nowhere is not used."""
    out = Path(out)
    empty, s1, s2 = out / "station_empty", out / "station_a", out / "station_b"
    for x in (empty, s1, s2):
        x.mkdir(parents=True, exist_ok=True)
    wis_name, chat = hf.DEFAULT_NAMES["wis"][0], hf.DEFAULT_NAMES["chatham"]
    shutil.copy(str(fx["wis"]), str(s2 / wis_name))
    shutil.copy(str(fx["chatham"]), str(s2 / chat[0]))          # the first name, in the second folder
    shutil.copy(str(fx["chatham"]), str(s1 / chat[1]))          # the second name, in the first folder
    saved = hf.SEARCH_DIRS
    try:
        hf.SEARCH_DIRS = [empty, s1, s2]
        found = {k: hf.find_input(k, None) for k in ("wis", "chatham", "ndbc")}
        check("search: an input not given is found in the first search folder holding one of its names",
              found == {"wis": s2 / wis_name, "chatham": s1 / chat[1], "ndbc": None},
              ", ".join(f"{k} {v.parent.name + '/' + v.name if v else None}" for k, v in found.items()))
        check("search: a file given is used, not the one the search would find",
              hf.find_input("wis", str(fx["wis"])) == Path(fx["wis"]))
        res, _ = quiet(lambda: build("2025-02-19", "2025-02-20", out / "found", adcp=fx["adcp"],
                                     adcp_navd88=fx["adcp_navd88"], download=False, plot=False))
        wv = pd.read_csv(res["waves_csv"])
        check("search: build_forcing reads the files the search found (WIS, Chatham), no NDBC",
              res["waves"]["inputs"]["wis"] == str(s2 / wis_name) and res["waves"]["inputs"]["ndbc"] is None
              and res["water_level"]["inputs"]["chatham"] == str(s1 / chat[1])
              and res["water_level"]["sources_used"].get("chatham_transfer", 0) > 0
              and set(wv["source"]) == {"hs=wis;tp=wis"},
              f"waves {sorted(set(wv['source']))}, water level {res['water_level']['sources_used']}")
    finally:
        hf.SEARCH_DIRS = saved


def test_isolation(work):
    """Every forcing.json this run wrote names only the fixture's files and the test's own folders."""
    work = Path(work).resolve()
    outside, n = [], 0
    for fj in sorted(work.rglob("forcing.json")):
        js = json.loads(fj.read_text())
        n += 1
        wl, wv = js.get("water_level") or {}, js.get("waves") or {}
        paths = (list((wl.get("inputs") or {}).values()) + list((wv.get("inputs") or {}).values())
                 + list((wv.get("still_water_reference") or {}).get("inputs") or []))
        outside += [f"{fj.parent.name}: {x}" for x in paths
                    if x and not str(Path(x).resolve()).startswith(str(work) + os.sep)]
    check(f"isolation: the {n} forcing.json written name only the fixture's and the test's own files",
          n >= 8 and not outside, "; ".join(outside)[:200])


def run_group(name, fn, *a):
    """One group of checks; an exception is a FAIL line (with its traceback), never the end of the run."""
    try:
        return fn(*a)
    except Exception as e:                      # noqa: BLE001 -- reported, and counted as a failure
        traceback.print_exc()
        check(f"{name}: ran to the end", False, f"{type(e).__name__}: {str(e)[:200]}")
        return None


def test_extrapolation_reach():
    """How far the transfer's error was measured, against how far out the Oct 2024 date lies
    (found building it end to end: the report said 'measured 0-9 weeks from its fit window', which
    covered October's 6-7 weeks; the tests predict 0-23 and 0-29 days beyond their fit windows)."""
    reach = hf.extrapolation_reach_days()
    check("reach: the two extrapolation tests predict 0-23 and 0-29 days beyond their fit windows",
          [tuple(r) for r in reach] == [(0, 23), (0, 29)], str(reach))
    a_ep = np.arange(ep("2024-12-09 18:00"), ep("2025-03-10 15:00") + 1, 3600.0)   # the real ADCP span
    t0, t1 = ep("2024-10-20"), ep("2024-10-27")                                  # the Oct 2024 window
    b = hf.beyond_record_days(t0, t1, a_ep)
    check("beyond: 20-26 Oct 2024 lies 44-51 days before the ADCP record",
          b is not None and round(b[0]) == 44 and round(b[1]) == 51, str(b))
    rel = hf.relation(t0, t1, a_ep, 3600.0)
    check("relation names both ends", rel.startswith("44-51 days BEFORE the ADCP record"), rel)
    check("beyond: a window inside the record is not beyond it",
          hf.beyond_record_days(ep("2025-01-16"), ep("2025-01-24"), a_ep) is None)
    note = "\n".join(hf.reach_note_lines({"extrapolation_reach_days": [list(r) for r in reach],
                                          "beyond_tested": "this period lies 44-51 days from the ADCP record"}))
    check("report: the reach as measured, not '0-9 weeks', and FURTHER THAN TESTED when it is",
          "0-9 weeks" not in note and "0-23 days before and 0-29 days after" in note
          and "FURTHER THAN TESTED" in note, note[-160:])
    note = "\n".join(hf.reach_note_lines({}))
    check("report: no FURTHER THAN TESTED line for a period inside the reach", "FURTHER THAN TESTED" not in note)


# The two tests' records as forcing.json keeps them (the real Oct 2024 build's numbers)
REAL_TESTS = [{"n": 534, "rms_m": 0.095, "bias_m": 0.0375, "p95_abs_m": 0.1876, "n_daytime": 111,
               "rms_daytime_m": 0.0849, "test": "fit 2025-01-01..2025-03-11 -> predict 2024-12-09..2025-01-01"},
              {"n": 688, "rms_m": 0.0765, "bias_m": -0.0045, "p95_abs_m": 0.1405, "n_daytime": 142,
               "rms_daytime_m": 0.0772, "test": "fit 2024-12-09..2025-02-10 -> predict 2025-02-10..2025-03-11"}]


def test_transfer_reach_by_side():
    """A period on one side of the fit window is judged by the test predicting on THAT side (found in
    review: October, before the record, was given the pooled RMS 0.085 m, which mixes in the forward
    test's smaller error, with no bias, and the farthest reach of either test, 29 days, as the bound)."""
    wl = {"method": "harmonic", "fit_window": ["2024-12-09 18:00", "2025-03-10 15:00"],
          "methods": {"harmonic": {"extrapolation_tests": REAL_TESTS}}}
    r = hf.transfer_reach(wl, "2024-10-20 00:00", "2024-10-27 00:00")
    check("reach: Oct 2024 lies before the fit window, 44-51 days out",
          r["side"] == "before" and round(r["beyond_days"][0]) == 44 and round(r["beyond_days"][1]) == 51,
          str(r["beyond_days"]))
    check("reach: judged by the backward test only (0-23 days before its fit window)",
          len(r["tests"]) == 1 and r["tests"][0]["reach_days"] == [0, 23] and r["reach_days"] == 23,
          str([t["test"] for t in r["tests"]]))
    check("reach: sigma = the backward test's larger RMS (0.095 m), not the pooled 0.085",
          r["sigma_m"] == 0.095 and "backward" in r["basis"] and "+0.037" in r["basis"], r["basis"][:120])
    check("reach: Oct 2024 is FURTHER THAN TESTED", r["beyond_tested"] and "FURTHER THAN TESTED" in r["basis"])
    txt = hf.beyond_tested_text(r)
    check("reach: its sentence names the backward test's 23 days, RMS and bias, not 0-29",
          "0-23 days before" in txt and "0.095" in txt and "+0.037" in txt and "29" not in txt, txt[:160])
    # epochs work as well as forcing.json's strings
    r2 = hf.transfer_reach(wl, ep("2024-10-20"), ep("2024-10-27"))
    check("reach: epochs or strings, same answer", r2["sigma_m"] == r["sigma_m"] and r2["side"] == "before")
    # a period after the record: the forward test (0-29 days after), inside its reach
    r = hf.transfer_reach(wl, "2025-03-20 00:00", "2025-03-27 00:00")
    check("reach: 20-26 Mar 2025 after the record, forward test, inside its 29 days",
          r["side"] == "after" and r["reach_days"] == 29 and not r["beyond_tested"] and r["sigma_m"] == 0.0772,
          f"{r['side']} {r['reach_days']} {r['sigma_m']}")
    # a period inside the fit window's span (hours the ADCP lacks): the pooled RMS stands
    r = hf.transfer_reach(wl, "2025-01-16 00:00", "2025-01-24 00:00")
    check("reach: a period inside the span has no side and no sigma of its own",
          r["side"] is None and r["sigma_m"] is None and not r["beyond_tested"])
    # no backward test ran (too few ADCP hours): no sigma of its own, said so
    wl2 = dict(wl, methods={"harmonic": {"extrapolation_tests": [{"test": REAL_TESTS[0]["test"],
                                                                  "skipped": "not enough ADCP hours"},
                                                                 REAL_TESTS[1]]}})
    r = hf.transfer_reach(wl2, "2024-10-20 00:00", "2024-10-27 00:00")
    check("reach: no backward test ran -> no side sigma, basis says why",
          r["side"] == "before" and r["sigma_m"] is None and "no extrapolation test" in r["basis"], r["basis"])
    check("reach: no fit window -> nothing", hf.transfer_reach({}, "2024-10-20", "2024-10-27")["side"] is None)


def test_week_mean_error():
    """The error a week-long window shares: each left-out week's mean residual over 13.5-18 UTC."""
    t = np.arange(ep("2025-01-06"), ep("2025-02-03"), 3600.0)      # four ISO weeks, hourly
    wk = (t - t[0]) // (7 * 86400)
    r = np.choose(wk.astype(int), [0.10, -0.05, 0.02, 0.0]) + 0.3 * np.sin(t / 3600.0)   # tide-like wiggle
    day = hf.daytime(t)
    want = [float(np.mean(r[(wk == k) & day])) for k in range(4)]
    w = hf.week_mean_error(t, r)
    check("week means: 4 weeks, min/max/RMS of the daytime week means",
          w["weeks"] == 4 and abs(w["min_m"] - min(want)) < 1e-4 and abs(w["max_m"] - max(want)) < 1e-4
          and abs(w["rms_m"] - float(np.sqrt(np.mean(np.square(want))))) < 1e-4, str(w))
    r2 = r.copy()
    r2[wk == 3] = np.nan
    check("week means: a week with no residuals is not counted", hf.week_mean_error(t, r2)["weeks"] == 3)
    check("week means: none at all -> None", hf.week_mean_error(t, np.full(len(t), np.nan)) is None)


def test_before_record_e2e(fx, out):
    """A period wholly BEFORE the ADCP record, end to end: every transferred row carries the backward
    test's error, flagged as not measured when the period lies further out than that test reaches;
    the report and forcing.json say so. The ADCP here starts 2025-01-01 and the backward test is
    shortened (fit 01-10.., predict 01-01..01-10: 0-9 days), so 10-16 Dec lies 16-22 days out."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    a = pd.read_csv(fx["adcp"])
    keep = hf.to_epoch(a["time"]) >= ep("2025-01-01")
    adcp = out / "sig1000_from_jan.csv"
    a[keep].to_csv(adcp, index=False)
    n = pd.read_csv(fx["adcp_navd88"])
    navd = out / "adcp_navd88_from_jan.csv"
    n[hf.to_epoch(n["time"]) >= ep("2025-01-01")].to_csv(navd, index=False)
    real_tests, real_dirs = hf.EXTRAPOLATION_TESTS, hf.SEARCH_DIRS
    empty = out / "empty"
    empty.mkdir(exist_ok=True)
    try:
        hf.EXTRAPOLATION_TESTS = [("2025-01-10", "2025-03-11", "2025-01-01", "2025-01-10"),
                                  ("2025-01-01", "2025-02-10", "2025-02-10", "2025-03-11")]
        hf.SEARCH_DIRS = [empty]                 # nothing found by default names: only the files given
        res, log = quiet(lambda: build("2024-12-10", "2024-12-16", out / "f", adcp=adcp,
                                        adcp_navd88=navd, chatham=fx["chatham"], wis=fx["wis"],
                                        download=False, plot=False))
    finally:
        hf.EXTRAPOLATION_TESTS, hf.SEARCH_DIRS = real_tests, real_dirs
    wl = pd.read_csv(res["water_level_csv"])
    js = json.loads((out / "f" / "forcing.json").read_text())["water_level"]
    m = js["methods"][js["method"]]
    back = [t for t in m["extrapolation_tests"] if t["test"].startswith("fit 2025-01-10")][0]
    want = max(back["rms_m"], back["rms_daytime_m"])
    tr = (wl["source"] == "chatham_transfer").to_numpy()
    check("before: every row is the transfer", tr.all(), str(set(wl["source"])))
    check("before: sigma_m = the backward test's larger RMS, not the pooled one",
          np.allclose(wl["sigma_m"][tr], round(want, 4)) and js["transfer_sigma_m"] == want,
          f"{js['transfer_sigma_m']} vs backward {want}, pooled {m['extrapolation_pooled']['rms_m']}")
    check("before: flagged as not measured (16-22 days out, the backward test reaches 9)",
          js["transfer_sigma_measured"] is False and "backward" in (js.get("beyond_tested") or "")
          and "0-9 days before" in js["beyond_tested"], str(js.get("beyond_tested"))[:160])
    check("before: the week-mean error is recorded", (js.get("week_mean_error") or {}).get("weeks", 0) >= 4,
          str(js.get("week_mean_error")))
    rep = (out / "f" / "forcing_report.txt").read_text()
    check("before: the report says NOT MEASURED, the basis, the week-mean error and FURTHER THAN TESTED",
          "NOT MEASURED this far out" in rep and "from the backward extrapolation test" in rep
          and "WEEK-MEAN ERROR" in rep and "FURTHER THAN TESTED" in rep)
    # the levels themselves are right (the synthetic transfer is exact up to 1 cm noise)
    a_ep, mt = fx["truth"][0], fx["truth"][1]
    w_ep = hf.to_epoch(wl["time"])
    on = (w_ep % 3600 == 0) & (w_ep >= a_ep[0])
    err = wl["water_level_navd88"].to_numpy(float)[on] - np.interp(w_ep[on], a_ep, mt)
    check("before: the transferred levels match the (hidden) ADCP", float(np.sqrt(np.mean(err ** 2))) < 0.04,
          f"RMS {float(np.sqrt(np.mean(err ** 2))):.4f} m")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--adcp", default=str(HERE / "sig1000_waves_ALL.csv"), help="the real ADCP record")
    ap.add_argument("--keep", default=None, help="write the fixture and outputs here and keep them")
    args = ap.parse_args()
    if not Path(args.adcp).exists():
        print(f"ERROR: ADCP record {args.adcp} not found (give --adcp)")
        return 2
    args.adcp = str(Path(args.adcp).resolve())      # the run's current folder is its work folder (below)
    if args.keep:                       # a new folder each run: never the files an earlier run left
        Path(args.keep).mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="hf_test_", dir=args.keep)).resolve()
    t = pd.Timestamp.now()
    start_dir = os.getcwd()
    os.chdir(str(work))                 # a file name never finds a file of the folder the test was started from
    print(f"work folder       : {work} (the current folder during the run)")
    print("fixture           : synthetic Chatham (a 1.35, lag -54 min, b -0.05, 1 cm noise, 36 h hole), "
          "synthetic WIS, ADCP with 2025-02-18..24 cut")
    print("station files     : none used (search folders and setup-share files pointed at the empty "
          "no_station folder; test_station_search tests the search on folders of its own)")
    global NO_STATION
    NO_STATION = work / "no_station"
    NO_STATION.mkdir()
    ISOLATED.update(gnssr_spline=str(NO_STATION / "usgs_spline_out.txt"),
                    gauge_archive=str(NO_STATION / "gauge_8447435.csv"),
                    waves_archive=str(NO_STATION / "waves_marconi.csv"))
    station_dirs, hf.SEARCH_DIRS = hf.SEARCH_DIRS, [NO_STATION]
    try:
        (work / "fixture").mkdir()
        # the fixture too: if it cannot be built, a FAIL line, the groups that need it are named as not run,
        # the others still run, and the summary line is printed
        fx = run_group("fixture", build_fixture, args.adcp, work / "fixture")
        need_fx = [("linear", test_linear_recovery, (fx,)),
                   ("harmonic", test_harmonic_recovery, (fx,)),
                   ("e2e", test_end_to_end, (fx, work / "e2e")),
                   ("adcp-only", test_adcp_only, (fx, work / "adcp_only")),
                   ("C = 0 currency", test_c0_same_currency, (fx, work / "c0_currency")),
                   ("partial and errors", test_partial_and_errors, (fx, work / "errors")),
                   ("download", test_download_fallback, (fx, work / "download")),
                   ("before record", test_before_record_e2e, (fx, work / "before_record"))]
        if fx is None:
            check("the groups that need the fixture", False,
                  f"not run (no fixture): {', '.join(n for n, _, _ in need_fx)}, search")
        else:
            for name, fn, a in need_fx:
                run_group(name, fn, *a)
        run_group("still water", test_still_water_share, work / "still_water")
        if fx is not None:
            run_group("search", test_station_search, fx, work / "search")
        run_group("reach", test_extrapolation_reach)
        run_group("reach by side", test_transfer_reach_by_side)
        run_group("week mean", test_week_mean_error)
        run_group("isolation", test_isolation, work)
    finally:
        hf.SEARCH_DIRS = station_dirs
        os.chdir(start_dir)
    n_fail = sum(1 for _, ok in RESULTS if not ok)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed in "
          f"{(pd.Timestamp.now() - t).total_seconds():.0f} s" + (f"; {n_fail} FAILED" if n_fail else ""))
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())

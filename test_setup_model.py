#!/usr/bin/env python3
"""
Tests For setup_model.py: One Setup Formula, And The Opt-In Slope-Aware Model
================================================================================
Plain python3 (no pytest needed on the station):

    python3 test_setup_model.py [--keep DIR]

WHAT.
  * ONE COPY, BIT FOR BIT. Every place the setup formula was written out
    (extract_elevation_contours.py, survey_products.setup_m / implied_c /
    compared_frames_scale, historical_forcing.setup_m and its GNSS-R setup
    share, dem_from_contours.carried_setup_coef / fit_setup_coefficient,
    compare_dem_survey.waves_vs_datum, validate_waves.py, gnssr_qc.py,
    runup_inversion.py, foreshore_slope.py) now calls setup_model.py. Each
    is checked against the expression as it was written before (copied
    here literally), over a grid of Hs 0-6 m by 0.01 m and Tp 1-25 s by
    0.1 s (the rows' own rounding), random unrounded values, zeros, NaN and
    negative heights: the same bits, not merely close. And the extractor's
    rows, end to end, carry the setup text the old expression gives.
  * The linear-theory deshoaling: the dispersion relation solved to
    machine precision, Ks at 21 m against the numbers setup_model.py's
    docstring quotes, deep water Ks = 1, the tide's effect, the options.
  * The foreshore slope: slope_grid equals foreshore_slope.surface_slope
    (to 1e-10, same NaN cells); a plane gives its slope.
  * The slope-aware frame (stockdon_frame): a plane DEM gives its slope and
    Stockdon's setup; a two-slope beach gives the upper slope at high water
    and the lower at low water, and at still water on the break the upper
    one (the band is still water + setup +/- S/2); no DEM, too few band
    cells and a slope outside BETA_RANGE fall back to C / k (C's setup, up
    to the deshoaling).
  * EvolvingSlope: the dated DEM nearest a frame's day, none beyond 7 days,
    one DEM by stem or file, a missing one.
  * The extractor: --setup-model c (the default) writes what it wrote
    (same columns, the old expression's setup); --setup-model stockdon
    writes setup_beta_f and setup_model before wave_source, the slope of
    the DEM, Stockdon's setup on the deshoaled H0, and falls back with a
    WARNING when no DEM serves; bad options exit 1 before any frame.
  * dem_from_contours.py records the setup model in <stem>_info.json only
    when the contours carry it (the station's C rows: no new key).
  * survey_products.apply_setup drops the slope-aware columns (its rows get
    C's setup); the cron's SETUP_MODEL switch builds the extractor's
    options, "" and "c" leave them as they were, an unknown value warns.

HOW. Synthetic detector CSVs, tide, waves and DEMs in a temporary folder
(--keep DIR to look at it); the scripts run as the cron runs them.

STATION FILES. None are read.
"""

import argparse
import csv
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
from contextlib import redirect_stdout
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import setup_model as sm          # noqa: E402

FAILED = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILED.append(msg)


def run(args, d):
    r = subprocess.run([sys.executable] + [str(a) for a in args], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
    return r.returncode, r.stdout


def same_bits(a, b):
    """The same float64 bits, NaN included (np.array_equal has no equal_nan on numpy 1.17)."""
    a, b = np.asarray(a, np.float64), np.asarray(b, np.float64)
    return a.shape == b.shape and a.tobytes() == b.tobytes()


def wave_grid():
    """Hs 0-6 m by 0.01 m and Tp 1-25 s by 0.1 s as the rows round them (Python floats), plus random
    unrounded values, and the odd ones (0, NaN, a negative height)."""
    hs = [round(0.01 * i, 2) for i in range(601)]
    tp = [round(1.0 + 0.1 * i, 1) for i in range(241)]
    rng = np.random.RandomState(7)
    extra = list(zip(rng.uniform(0, 8, 3000).tolist(), rng.uniform(0.5, 30, 3000).tolist()))
    extra += [(0.0, 0.0), (1.0, 0.0), (float("nan"), 8.0), (1.0, float("nan")), (-0.5, 8.0)]
    return hs, tp, extra


# ---------------------------------------------------------------------
# One copy, bit for bit
# ---------------------------------------------------------------------

def test_bit_identical():
    print("setup_model.py against every copy it replaced, as written before, bit for bit")
    import survey_products as sp
    import historical_forcing as hf
    import dem_from_contours as dfc
    import compare_dem_survey as cds
    import gnssr_qc
    import runup_inversion
    import foreshore_slope
    hs_g, tp_g, extra = wave_grid()
    pairs = [(h, t) for h in hs_g for t in tp_g] + extra
    coefs = (0.037, 0.0, 0.034, 0.0465)
    bad = {"extract": 0, "setup_m": 0, "setup_m_text": 0, "implied_c": 0, "frames_scale": 0, "beach": 0}
    with np.errstate(invalid="ignore"):
        for h, t in pairs:
            for c in coefs:
                # extract_elevation_contours.py main(): the row's setup and beach elevation
                old = c * np.sqrt(float(h) * 9.81 * float(t) ** 2 / (2 * np.pi))
                new = sm.setup_c(float(h), float(t), c)
                if not same_bits(old, new):
                    bad["extract"] += 1
                if np.isfinite(old) and (round(float(old), 4) != round(float(new), 4)
                                         or round(0.25 + float(old), 4) != round(0.25 + float(new), 4)):
                    bad["beach"] += 1
                # survey_products.setup_m (apply_setup, the live forcing's setups), on text and on floats
                if not same_bits(old, sp.setup_m(h, t, c)):
                    bad["setup_m"] += 1
                if not same_bits(old, sp.setup_m(repr(h), repr(t), c)):
                    bad["setup_m_text"] += 1
            # survey_products.implied_c
            x_old = np.sqrt(h * 9.81 * t ** 2 / (2 * np.pi))
            want = (0.3 / x_old if x_old > 0.05 else None) if np.isfinite(x_old) else None
            got = sp.implied_c(0.3, h, t)
            if not ((want is None and got is None) or (want is not None and got is not None
                                                       and same_bits(want, got))):
                bad["implied_c"] += 1
            # survey_products.compared_frames_scale
            if h == h and t == t and not same_bits(float(np.sqrt(h * 9.81 * t ** 2 / (2 * np.pi))),
                                                   float(sm.wave_scale(h, t))):
                bad["frames_scale"] += 1
    n = len(pairs)
    check(bad["extract"] == 0 and bad["beach"] == 0,
          f"extract_elevation_contours.py: setup_c on {n} Hs/Tp x {len(coefs)} C = the old expression "
          f"to the bit; the rounded setup and beach the same ({bad['extract']} / {bad['beach']} differ)")
    check(bad["setup_m"] == 0 and bad["setup_m_text"] == 0,
          f"survey_products.setup_m (floats and the rows' text) = the old expression ({bad['setup_m']} / "
          f"{bad['setup_m_text']} differ)")
    check(bad["implied_c"] == 0 and bad["frames_scale"] == 0,
          f"survey_products.implied_c and compared_frames_scale's sqrt(Hs L0) = the old ones "
          f"({bad['implied_c']} / {bad['frames_scale']} differ)")

    # arrays: historical_forcing, dem_from_contours, compare_dem_survey, validate_waves, gnssr_qc
    H = np.array([p[0] for p in pairs], float)
    T = np.array([p[1] for p in pairs], float)
    with np.errstate(invalid="ignore"):
        old = 0.037 * np.sqrt(H * 9.81 * T ** 2 / (2 * np.pi))
    check(same_bits(old, hf.setup_m(H, T, 0.037)) and same_bits(old, hf.setup_m(list(H), list(T)))
          and not np.isfinite(hf.setup_m(-0.5, 8.0)),
          "historical_forcing.setup_m (arrays, lists, its default C; NaN for a negative height) = the old one")
    with np.errstate(invalid="ignore"):
        old_x = np.sqrt(np.clip(H, 0, None) * 9.81 * T ** 2 / (2 * np.pi))
        new_x = sm.wave_scale(np.clip(H, 0, None), T)
    check(same_bits(old_x, new_x) and float(hf.G) == 9.81,
          "historical_forcing's GNSS-R setup share: its sqrt(Hs L0) = the old one (and hf.G is still 9.81)")
    with np.errstate(invalid="ignore"):
        phi_old = np.sqrt(H * 9.81 * T ** 2 / (2 * np.pi))
        check(same_bits(phi_old, sm.wave_scale(H, T)), "dem_from_contours.fit_setup_coefficient's phi = the old one")

    # dem_from_contours.carried_setup_coef against the function as it was
    rng = np.random.RandomState(3)
    info = {}
    for i in range(400):
        h, t = round(float(rng.uniform(0, 3)), 2), round(float(rng.uniform(4, 14)), 1)
        info[f"f{i}"] = {"setup": round(float(0.037 * np.sqrt(h * 9.81 * t ** 2 / (2 * np.pi))), 4) if i % 7 else
                         float("nan"), "hs": h, "tp": t}

    def carried_old(info):
        cs = []
        for v in info.values():
            if np.isfinite(v["setup"]) and v["setup"] != 0.0 and v["hs"] > 0 and v["tp"] > 0:
                cs.append(v["setup"] / np.sqrt(v["hs"] * 9.81 * v["tp"] ** 2 / (2 * np.pi)))
        return round(float(np.median(cs)), 4) if cs else None
    check(dfc.carried_setup_coef(info) == carried_old(info) == 0.037,
          f"dem_from_contours.carried_setup_coef = the old one ({dfc.carried_setup_coef(info)})")

    # dem_from_contours.fit_setup_coefficient against the function as it was
    names = np.array(sorted(info))
    pair_cell = rng.randint(0, 40, 4000)
    pair_frame = names[rng.randint(0, len(names), 4000)]
    finfo = {k: {"tide": float(rng.uniform(-1, 1)), "hs": v["hs"], "tp": v["tp"]} for k, v in info.items()}

    def fit_old(pair_cell, pair_frame, info, min_frames, coefs):
        tide = np.array([info[f]["tide"] for f in pair_frame])
        hs = np.array([info[f]["hs"] for f in pair_frame])
        tp = np.array([info[f]["tp"] for f in pair_frame])
        phi = np.sqrt(hs * 9.81 * tp ** 2 / (2 * np.pi))
        ok = np.isfinite(phi) & np.isfinite(tide)
        groups = dfc.cell_groups(pair_cell[ok])
        result = [(c,) + dfc.median_spread_by_cell(groups, tide[ok] + c * phi[ok], min_frames) for c in coefs]
        return result, int(ok.sum()), len(set(pair_frame[ok]))
    grid = np.round(np.arange(-0.04, 0.1201, 0.002), 3)
    check(dfc.fit_setup_coefficient(pair_cell, pair_frame, finfo, 3, coefs=grid)
          == fit_old(pair_cell, pair_frame, finfo, 3, grid),
          "dem_from_contours.fit_setup_coefficient = the old one (every C of its grid)")

    # compare_dem_survey.waves_vs_datum's phi
    k = 300
    off = rng.normal(0, 0.1, k)
    with redirect_stdout(io.StringIO()), np.errstate(invalid="ignore"):
        a, b, phi, has = cds.waves_vs_datum(off, H[-k:], T[-k:], rng.uniform(-1, 1, k))
    with np.errstate(invalid="ignore"):
        check(same_bits(phi, np.sqrt(H[-k:] * 9.81 * T[-k:] ** 2 / (2 * np.pi))),
              "compare_dem_survey.waves_vs_datum's sqrt(Hs L0) = the old one")

    # validate_waves.wave_bias's x, on pandas Series as there
    import pandas as pd
    f = pd.DataFrame({"hs": H[:5000], "tp": np.where(np.arange(5000) % 11 == 0, np.nan, T[:5000])})
    with np.errstate(invalid="ignore"):
        old_v = np.sqrt(f.hs * 9.81 * f.tp.fillna(f.tp.median()) ** 2 / (2 * np.pi))
        new_v = sm.wave_scale(f.hs, f.tp.fillna(f.tp.median()))
    check(isinstance(new_v, pd.Series) and same_bits(old_v.to_numpy(), new_v.to_numpy()),
          "validate_waves.py's x (pandas Series) = the old one")

    # gnssr_qc.wave_allowances against the function as it was
    w_ep = np.arange(len(H), dtype=float) * 600.0
    ep = w_ep + rng.uniform(-400, 400, len(w_ep))

    def allow_old(epochs, waves):
        epochs = np.asarray(epochs, dtype=float)
        w_ep, hs, tp = waves
        j = np.clip(np.searchsorted(w_ep, epochs), 0, len(w_ep) - 1)
        jm = np.clip(j - 1, 0, len(w_ep) - 1)
        k = np.where(np.abs(w_ep[jm] - epochs) < np.abs(w_ep[j] - epochs), jm, j)
        ok = np.abs(w_ep[k] - epochs) <= gnssr_qc.WAVE_MAX_GAP_S
        hl = np.maximum(hs[k] * 9.81 * tp[k] ** 2 / (2 * np.pi), 0.0)
        setup = gnssr_qc.SETUP_COEF * gnssr_qc.SETUP_BETA * np.sqrt(hl)
        r2 = 1.1 * (setup + np.sqrt(hl * (0.563 * gnssr_qc.SETUP_BETA ** 2 + 0.004)) / 2)
        return np.where(ok, setup, 0.0), np.where(ok, r2, 0.0)
    with np.errstate(invalid="ignore"):
        o1, o2 = allow_old(ep, (w_ep, H, T))
        n1, n2 = gnssr_qc.wave_allowances(ep, (w_ep, H, T))
    check(same_bits(o1, n1) and same_bits(o2, n2), "gnssr_qc.wave_allowances (setup and R2%) = the old one")

    # runup_inversion.h0_from_runup against the function as it was
    def h0_old(dz, beta, T0):
        L0 = 9.81 * np.asarray(T0, float) ** 2 / (2 * np.pi)
        beta = np.asarray(beta, float)
        k = 0.35 * beta + 0.5 * np.sqrt(0.563 * beta ** 2 + 0.004)
        dz = np.asarray(dz, float)
        out = np.full(dz.shape, np.nan)
        good = (dz > 0) & np.isfinite(beta) & (k > 0) & np.isfinite(L0)
        out[good] = (dz[good] / (1.1 * k[good])) ** 2 / L0[good]
        return out
    dz, beta = rng.uniform(-0.2, 2.0, len(T)), rng.uniform(0.01, 0.2, len(T))
    with np.errstate(divide="ignore", invalid="ignore"):          # the grid's Tp 0 rows
        check(same_bits(h0_old(dz, beta, T), runup_inversion.h0_from_runup(dz, beta, T)),
              "runup_inversion.h0_from_runup = the old one")

    # foreshore_slope's Stockdon setup and swash against the functions as they were
    def swash_old(H0, T0, beta):
        L0 = 9.81 * T0 ** 2 / (2 * np.pi)
        Sinc = 0.75 * beta * np.sqrt(H0 * L0)
        SIG = 0.06 * np.sqrt(H0 * L0)
        return np.sqrt(Sinc ** 2 + SIG ** 2)

    def setup_old(H0, T0, beta):
        L0 = 9.81 * T0 ** 2 / (2 * np.pi)
        return 0.35 * beta * np.sqrt(H0 * L0)
    Hp, Tp, Bp = np.abs(H[:20000]), T[:20000], rng.uniform(0.01, 0.3, 20000)
    with np.errstate(invalid="ignore"):
        ok_arr = (same_bits(swash_old(Hp, Tp, Bp), foreshore_slope.stockdon_swash(Hp, Tp, Bp))
                  and same_bits(setup_old(Hp, Tp, Bp), foreshore_slope.stockdon_setup(Hp, Tp, Bp)))
    ok_sc = all(same_bits(setup_old(0.9, 6.8, b), foreshore_slope.stockdon_setup(0.9, 6.8, b))
                and same_bits(swash_old(0.9, 6.8, b), foreshore_slope.stockdon_swash(0.9, 6.8, b))
                for b in (0.05, 0.106, 0.13))
    check(ok_arr and ok_sc, "foreshore_slope.stockdon_setup / stockdon_swash (arrays and scalars) = the old ones")
    check(same_bits(9.81 * 6.8 ** 2 / (2 * np.pi), sm.deep_water_wavelength(6.8)),
          "foreshore_slope's L0 (its Iribarren number) = the old one")


# ---------------------------------------------------------------------
# Deshoaling
# ---------------------------------------------------------------------

def test_deshoaling():
    print("linear-theory deshoaling (setup_model.wavenumber, shoaling_coefficient, deshoal_depth)")
    T = np.arange(3.0, 25.01, 0.1)
    for h in (2.0, 21.0, 60.0):
        k = sm.wavenumber(T, h)
        w2 = (2 * np.pi / T) ** 2
        res = np.max(np.abs(sm.G * k * np.tanh(k * h) - w2) / w2)
        check(res < 1e-12, f"dispersion relation solved at {h:g} m (max relative residual {res:.1e})")
    ks = sm.shoaling_coefficient([6.0, 8.0, 9.0, 12.0, 14.0], 21.0)
    want = (0.970, 0.920, 0.913, 0.939, 0.974)
    check(np.all(np.abs(ks - np.array(want)) < 0.0015),
          f"Ks at 21 m, Tp 6/8/9/12/14 s: {', '.join(f'{v:.3f}' for v in ks)} (docstring {want})")
    TT = np.arange(4.0, 16.01, 0.1)
    kk = sm.shoaling_coefficient(TT, 21.0)
    check(8.8 < TT[np.argmin(kk)] < 9.6 and abs(kk.min() - 0.913) < 0.001,
          f"its minimum {kk.min():.4f} at Tp {TT[np.argmin(kk)]:.1f} s")
    r = 1.0 / kk[(TT > 4.95) & (TT < 14.05)]
    check(abs(r.min() - 1.007) < 0.002 and abs(r.max() - 1.095) < 0.002
          and abs(np.sqrt(r.min()) - 1.003) < 0.001 and abs(np.sqrt(r.max()) - 1.047) < 0.001,
          f"Tp 5-14 s: H0/Hs {r.min():.3f}-{r.max():.3f}, sqrt {np.sqrt(r.min()):.4f}-{np.sqrt(r.max()):.4f}")
    r7 = 1.0 / kk[(TT > 6.95) & (TT < 12.05)]
    check(1.06 <= r7.min() and r7.max() <= 1.10, f"Tp 7-12 s: H0/Hs {r7.min():.3f}-{r7.max():.3f} (1.06-1.10)")
    sel = (TT > 5.95) & (TT < 12.05)
    dtide = np.maximum(np.abs(sm.shoaling_coefficient(TT, 19.8) - kk), np.abs(sm.shoaling_coefficient(TT, 22.2) - kk))
    check(dtide[sel].max() < 0.006 and dtide.max() < 0.011,
          f"the tide (+/-1.2 m): Ks moves {dtide[sel].max():.4f} at most over Tp 6-12 s, {dtide.max():.4f} to 16 s")
    check(abs(float(sm.shoaling_coefficient(3.0, 500.0)) - 1.0) < 1e-12,
          "deep water: Ks = 1")
    check(sm.deshoal(1.0, 8.0, None) == 1.0 and abs(float(sm.deshoal(1.0, 8.0, 21.0)) - 1 / ks[1]) < 1e-12,
          "deshoal: depth None leaves Hs; 21 m divides by Ks")
    opts = [sm.deshoal_depth("auto", "marconi:hs_best/tp_s")[0], sm.deshoal_depth("auto", "adcp_equiv:wvht_m/dpd_s")[0],
            sm.deshoal_depth("auto", "buoy:wvht_m/dpd_s")[0], sm.deshoal_depth("none", "marconi:hs_best/tp_s")[0],
            sm.deshoal_depth("18", "buoy:wvht_m/dpd_s")[0], sm.deshoal_depth(None, "")[0]]
    refused = 0
    for o in ("x", "-1", "0", "nan"):
        try:
            sm.deshoal_depth(o, "")
        except ValueError:
            refused += 1
    check(opts == [21.0, 21.0, None, None, 18.0, None] and refused == 4,
          f"deshoal_depth: auto 21 / 21 / deep water by currency, none, metres, unknown tag -> {opts}; "
          f"'x', '-1', '0', 'nan' refused ({refused}/4)")
    check(sm.model_tag(0.35, 21.0) == "stockdon2006:k=0.35:H0=Hs/Ks@21m"
          and sm.model_tag(0.4, None) == "stockdon2006:k=0.4:H0=Hs", "model_tag")


# ---------------------------------------------------------------------
# The foreshore slope
# ---------------------------------------------------------------------

def plane_dem(rows, cols, cell, slope, z_top=3.0, gaps=None):
    """A beach falling seaward along the columns (east) at `slope`, rows south to north."""
    x = (np.arange(cols) + 0.5) * cell
    dem = np.tile(z_top - slope * x, (rows, 1))
    if gaps is not None:
        dem[gaps] = np.nan
    return dem


def test_slope_grid():
    print("setup_model.slope_grid against foreshore_slope.surface_slope")
    import foreshore_slope
    rng = np.random.RandomState(11)
    dem = plane_dem(30, 40, 2.0, 0.1) + rng.normal(0, 0.05, (30, 40)) + 0.002 * np.arange(30)[:, None] ** 2
    dem[rng.uniform(size=dem.shape) < 0.3] = np.nan
    dem[:6, :9] = np.nan
    for window in (3, 5, 7):
        a = sm.slope_grid(dem, 2.0, window)
        b = foreshore_slope.surface_slope(dem, 2.0, window)
        both = np.isfinite(a) & np.isfinite(b)
        diff = float(np.max(np.abs(a[both] - b[both]))) if both.any() else 0.0
        check(np.array_equal(np.isfinite(a), np.isfinite(b)) and diff < 1e-10,
              f"window {window}: same cells ({int(both.sum())}), max difference {diff:.1e}")
    s = sm.slope_grid(plane_dem(20, 20, 2.0, 0.08), 2.0)
    check(np.allclose(s, 0.08, atol=1e-12), "a plane of slope 0.08 gives 0.08 everywhere (edges too)")
    fs = sm.ForeshoreSlope(plane_dem(20, 30, 2.0, 0.08), 2.0)
    b, n = fs.band(0.0, 1.0)
    b2, n2 = fs.band(10.0, 11.0)
    check(abs(b - 0.08) < 1e-12 and n > 20 and b2 is None and n2 == 0,
          f"ForeshoreSlope.band: 0.08 over {n} cells; a band above the DEM: None, 0 cells")


def test_stockdon_frame():
    print("setup_model.stockdon_frame: the slope of the swash band, the fallback")
    C, k = 0.037, sm.STOCKDON_SETUP_K
    fb = C / k
    fs = sm.ForeshoreSlope(plane_dem(40, 60, 2.0, 0.08), 2.0)
    setup, beta, src = sm.stockdon_frame(0.3, 1.2, 9.0, fs, fb, k, 21.0)
    h0 = float(sm.deshoal(1.2, 9.0, 21.0))
    check(src == "dem" and abs(beta - 0.08) < 1e-12
          and abs(setup - 0.35 * 0.08 * np.sqrt(h0 * sm.deep_water_wavelength(9.0))) < 1e-12,
          f"plane 0.08: beta_f {beta:.4f} from the DEM, setup {setup:.4f} m = 0.35 beta_f sqrt(H0 L0), H0 deshoaled")
    # two slopes: steep above +0.5 m, flat below
    x = (np.arange(80) + 0.5) * 2.0
    z = np.where(3.0 - 0.13 * x > 0.5, 3.0 - 0.13 * x, 0.5 - 0.06 * (x - (2.5 / 0.13)))
    fs2 = sm.ForeshoreSlope(np.tile(z, (40, 1)), 2.0)
    s_hi = sm.stockdon_frame(1.4, 0.8, 8.0, fs2, fb, k, None)
    s_lo = sm.stockdon_frame(-0.9, 0.8, 8.0, fs2, fb, k, None)
    check(s_hi[2] == "dem" and s_lo[2] == "dem" and abs(s_hi[1] - 0.13) < 0.005 and abs(s_lo[1] - 0.06) < 0.005
          and s_hi[0] > s_lo[0],
          f"a two-slope beach: high water beta_f {s_hi[1]:.3f} (upper face 0.13), low water {s_lo[1]:.3f} "
          f"(terrace 0.06): setups {s_hi[0]:.3f} / {s_lo[0]:.3f} m")
    # still water AT the break, Hs 1.5 m, Tp 10 s: still water + setup +/- S/2 lies on the upper face; centred on
    # still water, or +/- S wide, the band would reach down onto the terrace (0.069 / 0.087)
    s_br = sm.stockdon_frame(0.5, 1.5, 10.0, fs2, fb, k, None)
    check(s_br[2] == "dem" and abs(s_br[1] - 0.13) < 0.005,
          f"the swash band is still water + setup +/- S/2: still water at the break, beta_f {s_br[1]:.3f} "
          f"(the upper face's 0.13)")
    s_none = sm.stockdon_frame(0.3, 1.2, 9.0, None, fb, k, None, why_none="no DEM within 7 days")
    check(s_none[2] == "no DEM within 7 days" and s_none[1] == fb
          and abs(s_none[0] - float(sm.setup_c(1.2, 9.0, C))) < 1e-12,
          f"no DEM: beta_f = C/k = {fb:.4f}, the setup C's ({s_none[0]:.4f} m), the reason kept")
    s_dry = sm.stockdon_frame(8.0, 1.2, 9.0, fs, fb, k, 21.0)
    check(s_dry[2] == "too few DEM cells in the swash band" and s_dry[1] == fb
          and abs(s_dry[0] - C * np.sqrt(h0 * sm.deep_water_wavelength(9.0))) < 1e-12,
          "a band above the DEM: too few cells, C/k, C * sqrt(H0 L0)")
    cliff = sm.ForeshoreSlope(plane_dem(40, 60, 2.0, 0.5, z_top=20.0), 2.0)
    s_cl = sm.stockdon_frame(0.3, 1.2, 9.0, cliff, fb, k, None)
    check(s_cl[2].startswith("DEM slope 0.500 outside") and s_cl[1] == fb, f"a 0.5 cliff: '{s_cl[2]}', C/k kept")


def test_evolving_slope(d):
    print("setup_model.EvolvingSlope: the dated DEM nearest the frame's day")
    import dem_from_contours as dfc
    sd = Path(d) / "dems"
    sd.mkdir(parents=True)
    for end, slope in (("2026-10-05", 0.08), ("2026-10-12", 0.12)):
        dfc.write_ascii_grid(str(sd / f"dem_{end}_7d_dem.asc"), plane_dem(30, 40, 2.0, slope), 420000.0, 4638000.0, 2.0)
    (sd / "dem_2026-10-12_7d_spread.asc").write_text("not a DEM\n")
    ev = sm.EvolvingSlope(sd)
    ep = lambda s: datetime.strptime(s, "%Y-%m-%d %H").replace(tzinfo=timezone.utc).timestamp()  # noqa: E731
    a, _ = ev.for_epoch(ep("2026-10-01 15"))
    b, _ = ev.for_epoch(ep("2026-10-10 15"))
    c, why = ev.for_epoch(ep("2026-11-30 15"))
    check(len(ev.series) == 2 and a is not None and abs(a.band(0, 1)[0] - 0.08) < 1e-9
          and b is not None and abs(b.band(0, 1)[0] - 0.12) < 1e-9 and c is None and "within 7 days" in why,
          "2 dated DEMs (not the _spread grid): 1 Oct -> the 5 Oct window's (0.08), 10 Oct -> 12 Oct's (0.12), "
          f"30 Nov -> none ('{why}')")
    check(abs(ev.series[0][0] - ep("2026-10-02 12")) < 1e-6, "a 7-day window ending 5 Oct is centred 2 Oct noon UTC")
    one = sm.EvolvingSlope(sd / "dem_2026-10-05_7d")
    two = sm.EvolvingSlope(sd / "dem_2026-10-05_7d_dem.asc")
    gone = sm.EvolvingSlope(sd / "nothing")
    check(one.for_epoch(0)[0] is not None and two.for_epoch(0)[0] is not None
          and gone.for_epoch(0) == (None, "no DEM") and "NOT FOUND" in gone.what,
          "one DEM by stem or by file serves every frame; a missing one serves none and says NOT FOUND")


# ---------------------------------------------------------------------
# The extractor, end to end
# ---------------------------------------------------------------------

WAVES = [(0.42, 5.3), (0.87, 6.1), (1.13, 7.7), (1.38, 8.2), (1.66, 9.4), (1.91, 10.6), (2.24, 11.9), (2.57, 12.8)]


def extract_fixture(d):
    """8 c1 frames 30 min apart (rows rising with the tide), the tide, waves_marconi.csv with other Hs/Tp per frame."""
    d = Path(d)
    folder = d / "processed"
    folder.mkdir(parents=True, exist_ok=True)
    t0 = int(datetime(2026, 10, 9, 14, 0, tzinfo=timezone.utc).timestamp())
    rng = np.random.RandomState(0)
    epochs = []
    for i in range(len(WAVES)):
        ep = t0 + i * 1800
        epochs.append(ep)
        t = datetime.fromtimestamp(ep, tz=timezone.utc)
        lines = ["Column,Row,Confidence,Row_Precise,Peak_Sharpness,Has_Signal"]
        for c in range(400):
            row = 500.0 + 10.0 * i + 3.0 * np.sin(c / 7.0) + rng.normal(0, 0.8)
            lines.append(f"{c},{int(row)},0.95,{row:.3f},0.5,1")
        (folder / f"{ep}.{t:%a.%b.%d_%H_%M_%S}.GMT.{t:%Y}.CACO05.c1.timex.csv").write_text("\n".join(lines) + "\n")
    tide = ["time,level"] + [f"{datetime.fromtimestamp(t0 + k * 600, tz=timezone.utc):%Y-%m-%d %H:%M:%S},"
                             f"{0.2 + 0.0002 * k * 600 / 60:.4f}" for k in range(-6, 60)]
    (d / "tide.csv").write_text("\n".join(tide) + "\n")
    (d / "waves_marconi.csv").write_text("time_utc,epoch,hs_best,hs_sigma,tp_s\n"
                                         + "".join(f"x,{e},{h:.2f},0.1,{p:.1f}\n" for e, (h, p) in zip(epochs, WAVES)))
    return epochs


def read_rows(p):
    with open(p, newline="") as f:
        rd = csv.DictReader(f)
        return list(rd.fieldnames or []), list(rd)


def test_extract(d):
    print("extract_elevation_contours.py --setup-model: c (the default) as before, stockdon opt-in")
    d = Path(d) / "extract"
    extract_fixture(d)
    base = [HERE / "extract_elevation_contours.py", d / "tide.csv", "--processed-dir", d / "processed",
            "--waves", d / "waves_marconi.csv", "--setup-coef", "0.037"]
    rc, log = run(base + ["--output", d / "c.csv"], d)
    head, rows = read_rows(d / "c.csv") if rc == 0 else ([], [])
    want_head = ["offshore_hs_m", "offshore_tp_s", "setup_correction_m", "beach_elevation_navd88", "wave_source"]
    with np.errstate(invalid="ignore"):
        bad = [r for r in rows if r["setup_correction_m"] != str(round(float(
            0.037 * np.sqrt(float(r["offshore_hs_m"]) * 9.81 * float(r["offshore_tp_s"]) ** 2 / (2 * np.pi))), 4))]
    check(rc == 0 and len(rows) > 2000 and head[-5:] == want_head and "setup_beta_f" not in head
          and "setup_model" not in head and not bad and len({r["offshore_hs_m"] for r in rows}) == len(WAVES),
          f"default: {len(rows)} rows over {len(WAVES)} wave states, the same columns as before, every setup the "
          f"old expression's text ({len(bad)} differ)")
    rc2, log2 = run(base + ["--setup-model", "c", "--output", d / "c2.csv"], d)
    check(rc2 == 0 and (d / "c.csv").read_bytes() == (d / "c2.csv").read_bytes(),
          "--setup-model c writes the same file, byte for byte")

    import dem_from_contours as dfc
    sd = d / "dems"
    sd.mkdir()
    dfc.write_ascii_grid(str(sd / "dem_2026-10-12_7d_dem.asc"), plane_dem(40, 60, 2.0, 0.08), 420000.0, 4638000.0, 2.0)
    rc, log = run(base + ["--setup-model", "stockdon", "--setup-slope-dem", sd, "--output", d / "s.csv"], d)
    head, rows = read_rows(d / "s.csv") if rc == 0 else ([], [])
    worst = 0.0
    for r in rows:
        h, t = float(r["offshore_hs_m"]), float(r["offshore_tp_s"])
        want = 0.35 * 0.08 * np.sqrt(float(sm.deshoal(h, t, 21.0)) * sm.deep_water_wavelength(t))
        worst = max(worst, abs(float(r["setup_correction_m"]) - want),
                    abs(float(r["beach_elevation_navd88"]) - float(r["tide_elevation_navd88"]) - want))
    check(rc == 0 and head[-5:] == ["setup_correction_m", "beach_elevation_navd88", "setup_beta_f", "setup_model",
                                    "wave_source"]
          and all(r["setup_beta_f"] == "0.08" and r["setup_model"] == "stockdon2006:k=0.35:H0=Hs/Ks@21m" for r in rows)
          and worst <= 1.0e-4,
          f"stockdon: setup_beta_f 0.08 (the DEM's) and setup_model on every row, before wave_source; setup and "
          f"beach = 0.35 x 0.08 x sqrt(H0 L0), H0 deshoaled from 21 m (worst {worst:.1e} m, the 0.1 mm rounding)")
    check("slope-aware (Stockdon et al. 2006, k = 0.35; H0 = Hs / Ks, deshoaled from 21 m" in log
          and "Foreshore slope          : from 1 dated DEM(s)" in log and "took the fallback" not in log,
          "its log: the model, the deshoaling, the slope's source, no fallback")

    empty = d / "no_dems"
    empty.mkdir()
    rc, log = run(base + ["--setup-model", "stockdon", "--setup-slope-dem", empty, "--setup-deshoal-depth", "none",
                          "--output", d / "f.csv"], d)
    head, rows = read_rows(d / "f.csv") if rc == 0 else ([], [])
    _, crows = read_rows(d / "c.csv")
    same = rc == 0 and len(rows) == len(crows) and all(
        abs(float(a["setup_correction_m"]) - float(b["setup_correction_m"])) <= 1.0e-4 for a, b in zip(rows, crows))
    check(same and all(r["setup_beta_f"] == str(round(0.037 / 0.35, 4))
                       and r["setup_model"] == "stockdon2006:k=0.35:H0=Hs" for r in rows)
          and "WARNING: 8 frame(s) took the fallback slope beta_f = C/k = 0.106" in log and "no DEM 8" in log,
          "no DEM: every frame falls back to C/k = 0.106 with a WARNING naming why; with H0 = Hs its setup is C's "
          "(to the 0.1 mm rounding)")
    for extra, what in ((["--setup-model", "stockdon"], "no --setup-slope-dem"),
                        (["--setup-model", "stockdon", "--setup-slope-dem", sd, "--setup-deshoal-depth", "deep"],
                         "a bad --setup-deshoal-depth"),
                        (["--setup-model", "stockdon", "--setup-slope-dem", sd, "--setup-stockdon-k", "0"], "k = 0"),
                        (["--setup-model", "stockdon", "--setup-slope-dem", sd, "--setup-stockdon-k", "inf"],
                         "k = inf")):
        out = d / "refused.csv"
        rc, log = run(base + extra + ["--output", out], d)
        check(rc == 1 and "ERROR" in log and not out.exists() and "Found" not in log,
              f"{what}: exit 1 before any frame is read")


def ground_rows(model):
    """6 frames x 30 points crossing the same cells (a ground file as georectify.py writes it)."""
    head = ["source_file", "camera", "capture_time_utc", "capture_epoch", "pixel_column", "pixel_row",
            "tide_elevation_navd88", "offshore_hs_m", "offshore_tp_s", "setup_correction_m",
            "beach_elevation_navd88"] + (["setup_beta_f", "setup_model"] if model else []) + \
        ["wave_source", "easting_utm19", "northing_utm19"]
    t0 = int(datetime(2026, 10, 9, 14, 0, tzinfo=timezone.utc).timestamp())
    out = [",".join(head)]
    for i in range(6):
        ep = t0 + i * 3600
        tu = datetime.fromtimestamp(ep, tz=timezone.utc).isoformat()
        for j in range(30):
            vals = [f"f{i}", "c1", tu, str(ep), str(j), "500", "0.2", "1.0", "8.0", "0.37", "0.57"]
            vals += (["0.1", "stockdon2006:k=0.35:H0=Hs/Ks@21m"] if model else [])
            vals += ["marconi:hs_best/tp_s", f"{420000.5 + 2.0 * (j % 10):.3f}", f"{4638000.5 + 2.0 * (j // 10):.3f}"]
            out.append(",".join(vals))
    return "\n".join(out) + "\n"


def test_dem_record(d):
    print("dem_from_contours.py: the setup model in <stem>_info.json only when the contours carry it")
    d = Path(d) / "dem"
    d.mkdir(parents=True)
    got = {}
    for model in (False, True):
        g = d / f"ground_{int(model)}.csv"
        g.write_text(ground_rows(model))
        stem = d / f"dem_{int(model)}"
        rc, log = run([HERE / "dem_from_contours.py", g, stem, "--no-plot", "--no-diagnostics"], d)
        info = json.loads(Path(str(stem) + "_info.json").read_text()) if rc == 0 else {}
        got[model] = (rc, log, info)
    rc, log, info = got[False]
    check(rc == 0 and "setup_model" not in info and "Setup model" not in log,
          "the station's C rows: no setup_model key, no new log line")
    rc, log, info = got[True]
    check(rc == 0 and info.get("setup_model") == {"stockdon2006:k=0.35:H0=Hs/Ks@21m": 6}
          and "Setup model       : stockdon2006:k=0.35:H0=Hs/Ks@21m 6 frame(s)" in log,
          f"slope-aware rows: setup_model {info.get('setup_model')} in _info.json and the log")


def test_apply_setup(d):
    print("survey_products.apply_setup: the slope-aware columns go, C's setup comes")
    import survey_products as sp
    d = Path(d) / "apply"
    d.mkdir(parents=True)
    src = d / "rows.csv"
    src.write_text(ground_rows(True))
    sp.apply_setup(src, d / "out.csv", "c1", 0.037)
    head, rows = read_rows(d / "out.csv")
    want = round(float(sp.setup_m("1.0", "8.0", 0.037)), 4)
    check("setup_beta_f" not in head and "setup_model" not in head and head[-2:] == ["setup_correction_m",
                                                                                     "beach_elevation_navd88"]
          and all(float(r["setup_correction_m"]) == want for r in rows),
          f"setup_beta_f and setup_model dropped, every row C's setup {want} m")
    src0 = d / "rows0.csv"
    src0.write_text(ground_rows(False))
    sp.apply_setup(src0, d / "out0.csv", "c1", 0.037)
    head0, _ = read_rows(d / "out0.csv")
    check(head0 == [h for h in read_rows(src0)[0] if h not in ("setup_correction_m", "beach_elevation_navd88")]
          + ["setup_correction_m", "beach_elevation_navd88"], "the station's rows: the same columns as before")


def test_cron():
    print("waterline_timex_cron.sh: SETUP_MODEL builds the extractor's setup options")
    text = (HERE / "waterline_timex_cron.sh").read_text()
    check('\nSETUP_MODEL=""\n' in text and '\nSETUP_STOCKDON_K=""\n' in text and '\nSETUP_COEF="0.037"\n' in text,
          'defaults: SETUP_MODEL="" and SETUP_STOCKDON_K="" (the station stays on C = 0.037)')
    a = text.index('    setup_arg=""\n')
    b = text.index("    fi\n", a) + len("    fi\n")
    block = text[a:b]
    out = {}
    for model, k, coef in (("", "", "0.037"), ("c", "", "0.037"), ("stockdon", "", "0.037"),
                           ("stockdon", "0.41", "0.037"), ("Stockdon", "", "0.037"), ("stockdon", "", "")):
        script = ('log() { echo "LOG $*"; }\nwaves_arg="--waves /w.csv"\nDEM_SERIES=/b/archive/dems\n'
                  f'SETUP_COEF="{coef}"\nSETUP_MODEL="{model}"\nSETUP_STOCKDON_K="{k}"\n' + block
                  + 'echo "ARG $setup_arg"\n')
        r = subprocess.run(["bash", "-c", script], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           universal_newlines=True)
        out[(model, k, coef)] = r.stdout
    check(out[("", "", "0.037")].strip() == "ARG --setup-coef 0.037" == out[("c", "", "0.037")].strip(),
          'SETUP_MODEL "" or "c": --setup-coef 0.037 alone, nothing logged')
    s = out[("stockdon", "", "0.037")]
    check("ARG --setup-coef 0.037 --setup-model stockdon --setup-slope-dem /b/archive/dems\n" in s
          and "LOG setup: slope-aware (SETUP_MODEL=stockdon, k 0.35)" in s,
          "stockdon: the model and the DEM series, logged with k 0.35")
    check("--setup-stockdon-k 0.41" in out[("stockdon", "0.41", "0.037")], "SETUP_STOCKDON_K passes k")
    check("ARG --setup-coef 0.037\n" in out[("Stockdon", "", "0.037")] and "LOG WARNING: SETUP_MODEL" in
          out[("Stockdon", "", "0.037")], "an unknown value: a WARNING, C's setup")
    check(out[("stockdon", "", "")].strip() == "ARG", 'SETUP_COEF="" (no setup): no setup option, no warning')


def test_py38():
    print("Python 3.8 syntax")
    import ast
    for f in ("setup_model.py", "test_setup_model.py", "evaluate_setup_model.py", "test_evaluate_setup_model.py"):
        p = HERE / f
        if not p.exists():
            check(False, f"{f} missing")
            continue
        try:
            ast.parse(p.read_text(encoding="utf-8"), feature_version=(3, 8))
        except TypeError:                        # Python 3.7 has no feature_version
            ast.parse(p.read_text(encoding="utf-8"))
        check(True, f"{f} parses as Python 3.8")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="keep the outputs, in a new folder inside this one")
    args = ap.parse_args()
    if args.keep:
        Path(args.keep).mkdir(parents=True, exist_ok=True)
    d = os.path.realpath(tempfile.mkdtemp(prefix="setup_model_", dir=args.keep))
    print(f"test folder: {d}")
    try:
        test_bit_identical()
        test_deshoaling()
        test_slope_grid()
        test_stockdon_frame()
        test_evolving_slope(d)
        test_extract(d)
        test_dem_record(d)
        test_apply_setup(d)
        test_cron()
        test_py38()
    finally:
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED")
        return 1
    print("all setup model tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Tests For The Camera-Only Wave Height From The Runup Timestacks
=================================================================
Plain python3 (no pytest needed on the station):

    python3 test_timestack_wave_height.py

WHAT. timestack_wave_height.py turns the c2 swash edge into elevation and
inverts Stockdon et al. (2006) for H0. Checked here, on synthetic data with
a planted truth [SYNTH]:
  * the inversions: H0 -> R2 / S_IG / setup -> H0 round trip, and the
    published Eq. 19 against runup_inversion.h0_from_runup; the
    dissipative forms below xi0 = 0.3;
  * the elevation conversion: the ray / DEM intersection on the real c2
    calibration (a planar beach and its DEM), the extension past the DEM,
    the fractional-pixel elevation of the edge, and why the still-water
    plane would not do (its stretch);
  * the swash statistics: the 0.05 Hz split, the detrend (a tide in the
    record), setup and R2, the effective degrees of freedom of S_IG, an
    uprush past the line's end, edge samples off the profile;
  * the combination (its sigma against the closed form, absolute errors
    at the combined H0, the slope term), the period choice and the error
    model option;
  * status labelling: dark, low contrast (edge contrast; the window's
    sharpness), no period, no slope, and the wait for late inputs;
  * the inputs: the DEM series' 7-day limit, gaps inside a DEM, the
    extension's rising-slope guard, a level record with a gap;
  * the CSV: appended without duplicates, rewritten only by --reprocess,
    a changed header refused before any work, appended in chunks;
  * end to end on a synthetic camera scene (the image the real c2 camera
    would record of a planar beach): R2 and H0 from R2 recovered, with a
    level rising through the burst too; --since / --until;
  * --compare against hand-computed scores (held-out days, the OWG on
    the same stacks, the 50:50 blend, the counts of the decision rule);
  * the refactor of runup_from_timestack.py (swash_edge, peak_mask,
    build_parser, bridge_gaps, LAND_REF_FRACTION) leaves its outputs as
    they were;
  * not wired in: the cron and the plots never call it;
  * Python 3.8 syntax (the station's python3), and the --self-test.
Exit status 0 when every check passes.
"""

import ast
import csv
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import timestack_wave_height as tw          # noqa: E402
import runup_from_timestack as rft          # noqa: E402

failures = []
SH = np.array([420150.0, 4638400.0])        # the Jan 2025 0 m contour reference (sea_patch.py)


def check(name, cond, detail=None):
    print(f"  {'ok ' if cond else 'FAIL'} {name}" + (f"  ({detail})" if detail is not None and not cond else ""))
    if not cond:
        failures.append(name)


def calibration():
    import georectify
    io = georectify.load_intrinsics(HERE / "calibration" / "CACO05_c2_20240801_IO.yaml")
    eo = georectify.load_extrinsics(HERE / "calibration" / "CACO05_c2_20251113_EO-CV.yaml")
    return io, eo


def ground_line(io, eo, beta, along=150.0, d_from=-40.0, d_to=30.0):
    """Integer pixels of a cross-shore line on the planar beach z = -beta x (x seaward of the 0 m
    contour), LAND END FIRST as a .pix line might be, and each pixel's true elevation."""
    from view_reproject import ground_to_pixel
    sa, al = np.radians(81.0), np.radians(351.0)
    d = np.arange(d_from, d_to, 0.02)
    p0 = SH + along * np.array([np.sin(al), np.cos(al)])
    E, N, Z = p0[0] + d * np.sin(sa), p0[1] + d * np.cos(sa), -beta * d
    U, V, ok = ground_to_pixel(E, N, Z, io, eo)
    u, v, z = np.round(U[ok]).astype(int), np.round(V[ok]).astype(int), Z[ok]
    key = np.r_[True, (np.diff(u) != 0) | (np.diff(v) != 0)]
    grp = np.cumsum(key) - 1
    return u[key], v[key], np.bincount(grp, z) / np.bincount(grp)


def write_dem(path, beta, zlo=-1.0, zhi=1.6):
    """The planar beach as a 2 m DEM, kept only between zlo and zhi (a DEM covers the intertidal)."""
    x0, y0, cell, nx, ny = 419900.0, 4638300.0, 2.0, 200, 200
    xc = x0 + (np.arange(nx) + 0.5) * cell
    yc = y0 + (np.arange(ny) + 0.5) * cell
    XX, YY = np.meshgrid(xc, yc[::-1])
    dem = -beta * tw.seaward_distance(XX, YY)
    dem[(dem < zlo) | (dem > zhi)] = -9999
    with open(path, "w") as f:
        f.write(f"ncols {nx}\nnrows {ny}\nxllcorner {x0}\nyllcorner {y0}\ncellsize {cell}\nNODATA_value -9999\n")
        for row in dem:
            f.write(" ".join(f"{x:.4f}" for x in row) + "\n")
    return path


EPOCH = 1792500000          # 2026-10-20 12:40 UTC


def make_scene(tmp, h0=1.5, t0=10.0, beta=0.10, swl=0.3, blur=0.0, haze=0.0, dark=False, epoch=EPOCH,
               seed=5, swl_rise=0.0):
    """A synthetic c2 runup timestack of a planar beach: the swash elevation series (Stockdon-shaped,
    backwash at most 2 m/s on the beach) on a still water that rises swl_rise m through the burst
    (linearly, from swl at its start), every pixel below it foam (bright, textured), above it sand.
    -> dict of paths and the true swash series (m above still water)."""
    import cv2
    tmp = Path(tmp)
    rft._PIX_FILES.clear()           # runup_from_timestack caches the first --pix per camera
    io, eo = calibration()
    u, v, zp = ground_line(io, eo, beta)
    np.savetxt(tmp / "c2_timestack.pix", np.c_[u, v], fmt="%d")
    z_rel, _ = tw.synthetic_swash(h0, t0, beta, seed=seed)
    vmax = 2.0 * beta / tw.SAMPLE_HZ
    for k in range(1, len(z_rel)):
        z_rel[k] = max(z_rel[k], z_rel[k - 1] - vmax)
    level = swl + swl_rise * np.arange(len(z_rel)) / len(z_rel)
    z_abs = level + z_rel
    rng = np.random.default_rng(2)
    nt, npx = len(z_abs), len(u)
    img = np.empty((nt, npx), np.float32)
    tex = rng.normal(0, 12, (nt, npx))
    for k in range(nt):
        wet = zp < z_abs[k]
        img[k] = np.where(wet, 205 + tex[k], 105 + 25 * np.arange(npx) / npx + 3 * rng.standard_normal(npx))
    if blur:
        img = cv2.GaussianBlur(img, (0, 0), sigmaX=blur, sigmaY=0.01)
    if haze:
        img = (1 - haze) * img + haze * 190
    if dark:
        img = img * 0.08
    img = np.clip(img, 0, 255).astype(np.uint8)
    ras = tmp / "ras"
    ras.mkdir(exist_ok=True)
    name = f"{epoch}.c2.ras.tiff"
    cv2.imwrite(str(ras / name), cv2.cvtColor(img, cv2.COLOR_GRAY2BGR))
    dem = write_dem(tmp / "dem_test_dem.asc", beta)
    lv = tmp / "level.csv"
    lv.write_text("epoch,level_navd88\n" + "".join(
        f"{e},{swl + swl_rise * min(max((e - epoch) / 600.0, 0.0), 1.0):.6f}\n" for e in range(epoch - 7200, epoch + 7200, 60)))
    per = tmp / "period.csv"
    per.write_text("epoch,time_utc,filename,line,tp_s,tm01_s,ig_fraction,peak_ratio,status\n"
                   f"{epoch + 300},x,{name},1,{t0},{t0 * 0.8},0.2,5,ok\n")
    return {"ras": ras, "pix": tmp / "c2_timestack.pix", "dem": dem, "level": lv, "period": per,
            "name": name, "z_rel": z_rel, "z_pix": zp}


def run_scene(sc, tmp, out, extra=(), now=None, period=True, wl=None, waves=None):
    argv = ["--ras-dir", str(sc["ras"]), "--pix", str(sc["pix"]), "--lines", "1",
            "--dem", str(sc["dem"]), "--level-csv", str(sc["level"]),
            "--period-csv", str(period if isinstance(period, Path) else
                                sc["period"] if period else Path(tmp) / "none.csv"),
            "--waves-csv", str(waves or Path(tmp) / "none.csv"), "--buoy-csv", str(Path(tmp) / "none.csv"),
            "--owg-csv", str(Path(tmp) / "owg.csv"), "--output", str(out)] + list(extra)
    args = tw.build_parser().parse_args(argv)
    import io
    import contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        tw.run(args, now=now, wl=wl)
    return list(csv.DictReader(open(out))) if Path(out).exists() else []


class FakeLevel:
    """A still-water source with marconi_water_level.WaterLevel's .at() and its GPS attribute."""

    def __init__(self, level, source, gps=True):
        self.level, self.source = level, source
        self.s_ep = np.array([0.0]) if gps else None

    def at(self, epochs):
        t = np.atleast_1d(np.asarray(epochs, float))
        return np.full(len(t), self.level), np.array([self.source] * len(t), dtype=object)


def truth_stats(z_rel):
    pk = rft.peak_mask(z_rel)
    zd = tw.detrend(z_rel)
    return {"setup": float(z_rel.mean()), "r2": float(np.percentile(z_rel[pk], 98)),
            "sig": float(rft.band_sig(zd, 1e-9, 0.05)), "sinc": float(rft.band_sig(zd, 0.05, 2.0))}


# ---------------------------------------------------------------------

def test_inversions():
    print("inversions")
    worst = 0.0
    for h0 in (0.3, 0.8, 1.5, 3.0, 5.0):
        for t0 in (5.0, 8.0, 12.0, 16.0):
            for beta in (0.05, 0.10, 0.14):
                r2 = float(tw.stockdon_r2(h0, t0, beta))
                sig = float(tw.stockdon_sig(h0, t0))
                setup = float(tw.stockdon_setup(h0, t0, beta))
                back = (tw.h0_from_r2(r2, beta, t0), tw.h0_from_sig(sig, t0), tw.h0_from_setup(setup, beta, t0))
                worst = max(worst, max(abs(b / h0 - 1) for b in back))
    check("H0 -> R2 / S_IG / setup -> H0 round trip to 1e-9 (60 cases)", worst < 1e-9, f"{worst:.2e}")
    # Eq. 19 as published, written out here independently
    h0, t0, b = 2.0, 10.0, 0.1
    L0 = 9.81 * t0 ** 2 / (2 * np.pi)
    r2 = 1.1 * (0.35 * b * np.sqrt(h0 * L0) + np.sqrt(h0 * L0 * (0.563 * b ** 2 + 0.004)) / 2)
    check("stockdon_r2 is Eq. 19 as published", abs(float(tw.stockdon_r2(h0, t0, b)) - r2) < 1e-12)
    check("S_IG = 0.06 sqrt(H0 L0), no slope", abs(float(tw.stockdon_sig(h0, t0)) - 0.06 * np.sqrt(h0 * L0)) < 1e-12)
    check("setup = 0.35 beta sqrt(H0 L0)", abs(float(tw.stockdon_setup(h0, t0, b)) - 0.35 * b * np.sqrt(h0 * L0)) < 1e-12)
    # 1.1 (<eta> + S/2) from Eqs. 10-12 rounds differently from Eq. 19: not what is inverted
    from setup_model import stockdon_swash
    alt = 1.1 * (tw.stockdon_setup(h0, t0, b) + stockdon_swash(h0, t0, b) / 2)
    check("Eq. 19 differs from 1.1 (<eta> + S/2) by 0-3% (its rounded swash term)",
          0 < r2 / float(alt) - 1 < 0.03, f"{r2 / float(alt) - 1:.4f}")
    for name, v in (("R2 0", tw.h0_from_r2(0.0, 0.1, 8.0)), ("R2 < 0", tw.h0_from_r2(-0.2, 0.1, 8.0)),
                    ("setup < 0", tw.h0_from_setup(-0.1, 0.1, 8.0)), ("S_IG 0", tw.h0_from_sig(0.0, 8.0)),
                    ("no period", tw.h0_from_sig(0.5, np.nan)), ("no slope", tw.h0_from_setup(0.3, np.nan, 8.0))):
        check(f"no H0 from {name}", np.isnan(v), v)
    check("H0 scales as 1/T0^2 (a 25% longer period: 0.64 x H0)",
          abs(tw.h0_from_sig(0.6, 10.0) / tw.h0_from_sig(0.6, 8.0) - 0.64) < 1e-12)
    xi = tw.iribarren(0.1, 1.0, 8.0)
    check("Iribarren xi0 = beta / sqrt(H0/L0)", abs(xi - 0.1 / np.sqrt(1.0 / (9.81 * 64 / (2 * np.pi)))) < 1e-12)
    # the dissipative forms (xi0 < 0.3): R2 0.043, setup 0.016, S_IG 0.046 x sqrt(H0 L0), no slope
    h0, t0, b = 2.0, 10.0, 0.03
    x = np.sqrt(h0 * 9.81 * t0 ** 2 / (2 * np.pi))
    st = {"r2_m": 0.043 * x, "setup_m": 0.016 * x, "s_ig_m": 0.046 * x}
    check("dissipative forms: 1.1 (0.016 + 0.046 / 2) = 0.043 (Eq. 18)", abs(1.1 * (0.016 + 0.046 / 2) - 0.043) < 0.0005)
    check("dissipative: H0 -> R2 / S_IG / setup -> H0 round trip",
          all(abs(v / h0 - 1) < 1e-12 for v in (tw.h0_from_r2(st["r2_m"], b, t0, "dissipative"),
                                               tw.h0_from_sig(st["s_ig_m"], t0, "dissipative"),
                                               tw.h0_from_setup(st["setup_m"], b, t0, "dissipative"))))
    res = tw.estimate_h0(st, b, t0, True)
    check("near the threshold (true xi0 0.27; the intermediate forms read 0.33): both regimes agree with "
          "themselves, the dissipative estimates agree better: kept, H0 back exactly, noted unsure",
          res["regime"] == "dissipative" and abs(res["h0"] / h0 - 1) < 1e-9 and res["xi0"] < 0.3 and
          "regime unsure" in res["note"], (res["regime"], res["h0"], res["note"]))
    h0, t0, b = 3.0, 8.0, 0.02
    x = np.sqrt(h0 * 9.81 * t0 ** 2 / (2 * np.pi))
    st = {"r2_m": 0.043 * x, "setup_m": 0.016 * x, "s_ig_m": 0.046 * x}
    res = tw.estimate_h0(st, b, t0, True)
    check("well inside (xi0 0.12): only the dissipative forms agree with their own xi0: kept, H0 back, noted",
          res["regime"] == "dissipative" and abs(res["h0"] / h0 - 1) < 1e-9 and "dissipative kept" in res["note"]
          and "unsure" not in res["note"], (res["regime"], res["h0"], res["note"]))
    st = {"r2_m": float(tw.stockdon_r2(1.0, 9.0, 0.10)), "setup_m": float(tw.stockdon_setup(1.0, 9.0, 0.10)),
          "s_ig_m": float(tw.stockdon_sig(1.0, 9.0))}
    res = tw.estimate_h0(st, 0.10, 9.0, True)
    check("Marconi's 0.10 and H0 1 m (xi0 1.1): the intermediate forms, no note",
          res["regime"] == "intermediate" and abs(res["h0"] - 1.0) < 1e-9 and res["note"] == "", res)
    res = tw.estimate_h0(st, 0.10, 9.0, False)
    check("a planar profile: S_IG only, R2 and setup not inverted",
          res["used"] == "sig" and np.isnan(res["est"]["r2"]) and np.isnan(res["est"]["setup"]), res["used"])


def test_combination():
    print("combination")
    h, s, used = tw.combine({"r2": 1.2, "sig": 1.2, "setup": 1.2}, 9.0, 0.1)
    check("equal estimates combine to themselves", abs(h - 1.2) < 1e-9 and used == "r2+setup+sig", (h, used))
    h1, s1, _ = tw.combine({"sig": 1.0}, 9.0, 0.1)
    check("a single estimate passes through", abs(h1 - 1.0) < 1e-12)
    check("three estimates are more certain than one", s < s1, (s, s1))
    h, _, used = tw.combine({"r2": np.nan, "sig": 1.0, "setup": -1.0}, 9.0, 0.1)
    check("NaN and non-positive estimates are left out", used == "sig" and abs(h - 1.0) < 1e-12, used)
    h, s, used = tw.combine({}, 9.0, 0.1)
    check("no estimate: NaN", np.isnan(h) and used == "")
    # the combined value lies between, nearer the estimate given the smaller error
    sc = {"r2": (0.05, 0.0), "sig": (0.5, 0.0), "setup": (0.5, 0.0)}
    h, _, _ = tw.combine({"r2": 1.0, "sig": 2.0}, 9.0, 0.1, sc)
    check("weights follow the error model (nearer the better estimate)", 1.0 < h < 1.05, h)
    # errors evaluated at the combined H0: swapping which estimate is high does not move the weights
    sc = {"r2": (0.2, 0.0), "sig": (0.2, 0.0), "setup": (0.0, 0.3)}
    a, _, _ = tw.combine({"r2": 1.0, "sig": 2.0}, 9.0, 0.1, sc)
    b, _, _ = tw.combine({"r2": 2.0, "sig": 1.0}, 9.0, 0.1, sc)
    check("purely relative errors: a geometric mean, whichever reads high", abs(a - b) < 1e-9 and
          abs(a - np.sqrt(2.0)) < 1e-9, (a, b))
    sc = tw.parse_scatter("r2=0.1:0.2,setup=0:0.5:0.3")
    check("--scatter parses relative:absolute[:slope] per estimate (slope left out: the default's)",
          sc["r2"] == (0.1, 0.2, tw.SLOPE_REL) and sc["setup"] == (0.0, 0.5, 0.3) and sc["sig"] == tw.SCATTER["sig"], sc)
    try:
        tw.parse_scatter("foo=1:2")
        bad = False
    except ValueError:
        bad = True
    check("--scatter refuses an unknown estimate", bad)
    check("the default level error is still water and the DEM's level, 0.28 m",
          abs(tw.LEVEL_SIGMA_M - np.hypot(0.12, 0.258)) < 1e-12 and abs(tw.LEVEL_SIGMA_M - 0.2845) < 0.001)
    from dem_from_contours import SIGMA_SETUP_M
    check("the default error model: R2 26% and the level; S_IG 26%, no level; setup its scatter and the level; "
          "slope 15% each", tw.SCATTER == {"r2": (0.38 / 1.44, tw.LEVEL_SIGMA_M, 0.15), "sig": (0.38 / 1.44, 0.0, 0.15),
                                           "setup": (0.0, float(np.hypot(SIGMA_SETUP_M, tw.LEVEL_SIGMA_M)), 0.15)},
          tw.SCATTER)
    # the sigma of ln H0 against the closed form, written out here
    t0, b = 9.0, 0.10
    L0 = 9.81 * t0 ** 2 / (2 * np.pi)
    root = np.sqrt(0.563 * b ** 2 + 0.004)
    m_r2 = 2 * b * (0.35 + 0.563 * b / (2 * root)) / (0.35 * b + root / 2)
    check("R2's slope sensitivity: 2 beta k'/k = 1.52 at beta 0.10 (H0 ~ beta_f^-1.5)",
          abs(tw.slope_sensitivity("r2", b) - m_r2) < 1e-12 and abs(m_r2 - 1.515) < 0.01, tw.slope_sensitivity("r2", b))
    check("slope sensitivity: setup 2, S_IG 2 (its elevations), dissipative R2 and setup 0",
          tw.slope_sensitivity("setup", b) == 2.0 and tw.slope_sensitivity("sig", b) == 2.0 and
          tw.slope_sensitivity("r2", b, "dissipative") == 0.0 and tw.slope_sensitivity("setup", b, "dissipative") == 0.0)
    extra = 0.1
    _, s_sig, _ = tw.combine({"sig": 1.0}, t0, b, sig_extra=extra)
    want_sig = np.hypot(2 * np.hypot(0.38 / 1.44, extra), 2 * 0.15)
    check("h0_sigma_rel, S_IG alone: hypot(2 x hypot(26%, sampling and off-profile), 2 x 15% slope)",
          abs(s_sig - want_sig) < 1e-12, (s_sig, want_sig))
    q = 1.1 * (0.35 * b * np.sqrt(1.2 * L0) + np.sqrt(1.2 * L0 * (0.563 * b ** 2 + 0.004)) / 2)
    s_r2 = np.hypot(2 * np.hypot(0.38 / 1.44 * q, tw.LEVEL_SIGMA_M) / q, m_r2 * 0.15)
    s_sg = np.hypot(2 * np.hypot(0.38 / 1.44, extra), 2 * 0.15)
    h, s2, used = tw.combine({"r2": 1.2, "sig": 1.2}, t0, b, sig_extra=extra)
    check("h0_sigma_rel, R2 and S_IG: 1 / sqrt(1/s_r2^2 + 1/s_sig^2), each from its closed form",
          abs(s2 - 1 / np.sqrt(1 / s_r2 ** 2 + 1 / s_sg ** 2)) < 1e-12 and abs(h - 1.2) < 1e-12, (s2, s_r2, s_sg))
    # absolute errors, evaluated at the COMBINED H0: both statistics go as sqrt(H0) there, so their
    # weights keep the ratio rho = (R2 / setup)^2 at any H0 and the result has a closed form
    sc = {"r2": (0.0, 0.3, 0.0), "setup": (0.0, 0.3, 0.0), "sig": (0.2, 0.0, 0.0)}
    rho = (float(tw.stockdon_r2(1.0, t0, b)) / float(tw.stockdon_setup(1.0, t0, b))) ** 2
    for hi_lo in ((1.0, 2.0), (2.0, 1.0)):
        h, _, _ = tw.combine({"r2": hi_lo[0], "setup": hi_lo[1]}, t0, b, sc)
        want = np.exp((rho * np.log(hi_lo[0]) + np.log(hi_lo[1])) / (rho + 1))
        check(f"absolute errors at the combined H0 (R2 {hi_lo[0]:g}, setup {hi_lo[1]:g} m): the closed form "
              f"exp((rho ln R2 + ln setup) / (rho + 1))", abs(h - want) < 1e-9, (h, want))
    check("sampling error of S_IG: sqrt(2/nu)/2, 9.3% at 58 degrees of freedom",
          abs(tw.sig_sampling_rel(58) - np.sqrt(2 / 58) / 2) < 1e-12 and abs(tw.sig_sampling_rel(58) - 0.0928) < 0.0005)
    check("S_IG's extra error: its sampling (the record's own nu) and the share off the profile, in quadrature",
          abs(tw.sig_extra_rel({"ig_dof": 18.0, "missing": 0.05}, 1200) - np.hypot(np.sqrt(2 / 18) / 2, 0.05)) < 1e-12
          and abs(tw.sig_extra_rel({"ig_dof": np.nan, "missing": 0.0}, 1200) - np.sqrt(2 / 58) / 2) < 1e-12)


def test_spectrum_and_stats():
    print("swash statistics")
    n = 1200
    t = np.arange(n) / 2.0
    ig, inc = 0.20 * np.sin(2 * np.pi * t / 60.0), 0.10 * np.sin(2 * np.pi * t / 8.0 + 0.4)
    tide = 0.10 * t / t[-1]
    z = 0.5 + ig + inc + tide
    swl = 0.2 + tide                                  # still water rises with the tide
    edge = 100 + 50 * (z - 0.2)                       # any monotone mapping to pixels
    above = np.zeros(n, bool)
    st = tw.swash_stats(z, swl, edge, above, above, 1.0, 20, 0.8)
    check("S_IG = 4 sigma of the < 0.05 Hz part (0.2 m sine: 0.566 m)", abs(st["s_ig_m"] / (4 * 0.2 / np.sqrt(2)) - 1) < 0.01,
          st["s_ig_m"])
    check("S_inc = 4 sigma of the > 0.05 Hz part (0.1 m sine: 0.283 m)", abs(st["s_inc_m"] / (4 * 0.1 / np.sqrt(2)) - 1) < 0.01,
          st["s_inc_m"])
    check("total S^2 = S_inc^2 + S_IG^2", abs(st["s_m"] ** 2 - st["s_inc_m"] ** 2 - st["s_ig_m"] ** 2) < 1e-6)
    check("setup = mean above still water (0.3 m)", abs(st["setup_m"] - 0.3) < 0.005, st["setup_m"])
    st2 = tw.swash_stats(z, np.full(n, np.nan), edge, above, above, 1.0, 20, 0.8)
    check("no water level: no R2, no setup, the spectrum still (detrended)",
          np.isnan(st2["r2_m"]) and np.isnan(st2["setup_m"]) and abs(st2["s_ig_m"] - st["s_ig_m"]) < 0.01,
          (st2["r2_m"], st2["setup_m"], st2["s_ig_m"]))
    drift = z + 0.4 * t / t[-1]                       # no still water to subtract, and a 0.5 m drift
    st3 = tw.swash_stats(drift, np.full(n, np.nan), edge, above, above, 1.0, 20, 0.8)
    raw = rft.band_sig(drift, 1e-9, 0.05)             # band_sig alone removes only the mean
    check("the detrend: a 0.5 m drift in 10 min leaves S_IG within 2%, without it +10% or more",
          abs(st3["s_ig_m"] / st["s_ig_m"] - 1) < 0.02 and raw > 1.10 * st["s_ig_m"], (st3["s_ig_m"], raw))
    check("ig_dof: 58 degrees of freedom in 10 min at 2 Hz", tw.ig_dof(1200) == 58, tw.ig_dof(1200))
    # R2 = the 98th percentile of the uprush maxima: 50 triangular uprushes of known heights
    heights = np.linspace(0.2, 1.2, 50)
    rng = np.random.default_rng(1)
    rng.shuffle(heights)
    tri = np.concatenate([h * (1 - np.abs(np.linspace(-1, 1, 25))) for h in heights])
    e = 10 + 100 * tri
    st = tw.swash_stats(tri, np.zeros(len(tri)), e, np.zeros(len(tri), bool), np.zeros(len(tri), bool), 1.0, 20, 0.8)
    check("R2 = 98th percentile of the uprush maxima", abs(st["r2_m"] - np.percentile(heights, 98)) < 1e-9,
          (st["r2_m"], np.percentile(heights, 98)))
    check("one peak per uprush", st["n_peaks"] == 50, st["n_peaks"])
    # uprushes above the profile: R2 known to be above it -> none
    zz = tri.copy()
    above = np.zeros(len(zz), bool)
    hi = zz > 0.9
    zz[hi] = np.nan
    above[hi] = True
    st = tw.swash_stats(zz, np.zeros(len(zz)), e, above, np.zeros(len(zz), bool), 1.0, 20, 0.8)
    check("the 98th percentile above the profile: no R2", np.isnan(st["r2_m"]) and "R2 above the profile" in st["notes"],
          st["notes"])
    check(f"... and no setup ({hi.mean():.0%} of the edge off the profile, over 2%)",
          hi.mean() > 0.02 and np.isnan(st["setup_m"]), hi.mean())
    h100 = np.linspace(0.2, 1.2, 100)
    rng.shuffle(h100)
    zz = np.concatenate([h * (1 - np.abs(np.linspace(-1, 1, 25))) for h in h100])
    e100 = 10 + 100 * zz
    above = np.zeros(len(zz), bool)
    top = zz > np.sort(h100)[-2] + 1e-9                # only the highest of 100 uprushes runs off it
    zz[top] = np.nan
    above[top] = True
    st = tw.swash_stats(zz, np.zeros(len(zz)), e100, above, np.zeros(len(zz), bool), 1.0, 20, 0.8)
    check("one of 100 uprushes above the profile, the 98th percentile below it: R2 still read",
          abs(st["r2_m"] - np.percentile(h100, 98)) < 1e-9, st["r2_m"])
    st = tw.swash_stats(tri, np.zeros(len(tri)), e, np.zeros(len(tri), bool), np.zeros(len(tri), bool), 1.0, 60, 0.8)
    check("fewer uprushes than min-peaks: no R2", np.isnan(st["r2_m"]))


def test_ig_dof():
    print("S_IG sampling error: effective degrees of freedom")
    rng = np.random.default_rng(7)
    n = 1200
    f = np.fft.rfftfreq(n, 0.5)
    shapes = {"flat": ((f >= 0.005) & (f < 0.05)).astype(float),
              "peaked at 0.01 Hz": np.exp(-0.5 * ((f - 0.01) / 0.004) ** 2) * (f > 0) * (f < 0.05)}
    for name, S in shapes.items():
        sig, pred, nus = [], [], []
        for _ in range(400):
            X = np.sqrt(S) * (rng.normal(size=len(f)) + 1j * rng.normal(size=len(f)))
            zd = tw.detrend(np.fft.irfft(X, n))
            sig.append(rft.band_sig(zd, 1e-9, 0.05))
            nus.append(tw.ig_dof_eff(zd))
            pred.append(tw.sig_sampling_rel(nus[-1]))
        mc = float(np.std(sig) / np.mean(sig))
        print(f"    {name}: S_IG sd/mean {mc:.3f} (Monte Carlo, 400 records), predicted median {np.median(pred):.3f}, "
              f"nu median {np.median(nus):.1f}; the flat-spectrum count's {tw.sig_sampling_rel(tw.ig_dof(n)):.3f}")
        check(f"{name}: the predicted sampling error within 20% of the Monte Carlo's", abs(np.median(pred) / mc - 1) < 0.2,
              (np.median(pred), mc))
        if name.startswith("peaked"):
            check("peaked: well above the flat-spectrum 9% (fewer degrees of freedom)", np.median(pred) > 0.14, np.median(pred))
    t = np.arange(n) / 2.0
    check("a single infragravity sine: about 2 degrees of freedom", tw.ig_dof_eff(np.sin(2 * np.pi * t / 100.0)) < 4,
          tw.ig_dof_eff(np.sin(2 * np.pi * t / 100.0)))
    check("never more than the flat-spectrum count", tw.ig_dof_eff(rng.normal(size=n)) <= tw.ig_dof(n))
    x = np.sin(2 * np.pi * t / 60.0)
    x[100:150] = np.nan
    check("a gap band_sig would not bridge: NaN", np.isnan(tw.ig_dof_eff(x)))


def test_line_end_and_off_profile():
    print("the line's landward end and the edge off the profile")
    n_px = 900
    p = np.arange(n_px, dtype=float)
    z_line = -0.5 + 2.5 * p / (n_px - 1)                     # a profile over the WHOLE line
    swl = 0.8
    z_rel, tg = tw.synthetic_swash(1.5, 10.0, 0.10)
    z_abs = swl + z_rel
    edge = np.interp(z_abs, z_line, p)                        # the camera's edge stops at the line's end
    own = float(np.percentile(z_rel[rft.peak_mask(z_rel)], 98))
    zz, above, below = tw.elevation_at(edge, z_line)
    st = tw.swash_stats(zz, np.full(len(zz), swl), edge, above, below, 1.0, 20, 0.8)
    check("without the line-end rule an uprush past the line reads as the line's end: R2 capped low",
          st["r2_m"] < own - 0.1 and not above.any(), (st["r2_m"], own))
    cut = tw.line_end(n_px)
    check("line_end: find_edge's land reference (landward 15%) less 2 px", cut == int(0.85 * n_px) - 2, cut)
    zz, above, below = tw.elevation_at(edge, z_line, cut)
    past = edge >= cut
    st = tw.swash_stats(zz, np.full(len(zz), swl), edge, above, below, 1.0, 20, 0.8, past)
    check("with it: censored like an uprush above the profile, no R2, 'R2 above the line'",
          np.isnan(st["r2_m"]) and "R2 above the line" in st["notes"] and above[past].all() and np.isnan(zz[past]).all(),
          st["notes"])
    check("... and the censored samples count as off the profile", abs(st["missing"] - past.mean()) < 1e-12 and
          any("off the profile" in x for x in st["notes"]), (st["missing"], past.mean()))
    # S_IG with some of the edge off the profile: noted (it reads low), widened in the combination
    z_rel, _ = tw.synthetic_swash(3.0, 12.0, 0.09)
    zz = z_rel.copy()
    top = zz > np.quantile(zz, 0.95)
    zz[top] = np.nan
    e = np.arange(len(zz), dtype=float)
    above = top.copy()
    st0 = tw.swash_stats(z_rel, np.zeros(len(zz)), e, np.zeros(len(zz), bool), np.zeros(len(zz), bool), 1.0, 20, 0.8)
    st = tw.swash_stats(zz, np.zeros(len(zz)), e, above, np.zeros(len(zz), bool), 1.0, 20, 0.8)
    print(f"    top 5% off the profile: S_IG {st['s_ig_m'] / st0['s_ig_m'] - 1:+.1%}")
    check("5% of the edge off the profile: S_IG still given, 2-6% low, and noted",
          -0.06 < st["s_ig_m"] / st0["s_ig_m"] - 1 < -0.02 and any("S_IG with 5.0% of the edge off the profile" in x
                                                                     for x in st["notes"]), st["notes"])


def test_elevation_conversion():
    print("elevation conversion")
    z_line = np.array([np.nan, np.nan, 0.0, 0.1, 0.25, 0.5, np.nan])
    z, above, below = tw.elevation_at([2.0, 2.5, 3.2, 5.0, 5.5, 0.4, np.nan], z_line)
    check("fractional pixel: linear between pixels", np.allclose(z[:4], [0.0, 0.05, 0.13, 0.5]), z)
    check("off the landward end: NaN, flagged above", np.isnan(z[4]) and above[4] and not below[4])
    check("off the seaward end: NaN, flagged below", np.isnan(z[5]) and below[5] and not above[5])
    check("no edge: NaN, not flagged", np.isnan(z[6]) and not above[6] and not below[6])

    io, eo = calibration()
    beta = 0.10
    u, v, zp = ground_line(io, eo, beta)
    u, v, zp = u[::-1], v[::-1], zp[::-1]                    # sea end first, as swash_edge gives it
    pl = tw.line_profile(u, v, io, eo, beta=beta)
    zplane = -beta * tw.seaward_distance(pl["E"], pl["N"])
    check("planar: each ray meets the plane (z = -beta x at its own point)",
          pl["kind"] == "planar" and np.nanmax(np.abs(pl["z"] - zplane)) < 1e-6)
    dz = np.abs(pl["z"] - zp)
    check("planar: the scene's pixel elevations within 4 cm, median 1 cm (pixels are whole: a ray "
          "at 230 m moves ~1.4 m along itself per pixel)", np.nanmax(dz) < 0.04 and np.nanmedian(dz) < 0.01,
          (np.nanmax(dz), np.nanmedian(dz)))
    with tempfile.TemporaryDirectory() as tmp:
        dem = tw.Dem(write_dem(Path(tmp) / "d_dem.asc", beta))
        pr = tw.line_profile(u, v, io, eo, dem=dem, max_extrapolate=10.0)
    cov = np.isfinite(pr["z"]) & ~pr["extrap"]
    ext = pr["extrap"]
    check("DEM: the ray / DEM crossing is the ray / plane crossing within 1 mm", pr["kind"] == "dem" and
          np.nanmax(np.abs(pr["z"][cov] - pl["z"][cov])) < 0.001, np.nanmax(np.abs(pr["z"][cov] - pl["z"][cov])))
    check("DEM: covers its own range, -1.0 to +1.6 m (strict bilinear: within a cell of its edges)",
          -1.0 <= np.nanmin(pr["z"][cov]) < -0.7 and 1.3 < np.nanmax(pr["z"][cov]) <= 1.6,
          (np.nanmin(pr["z"][cov]), np.nanmax(pr["z"][cov])))
    check("extension past the DEM follows the plane within 2 cm", ext.any() and
          np.nanmax(np.abs(pr["z"][ext] - pl["z"][ext])) < 0.02,
          np.nanmax(np.abs(pr["z"][ext] - pl["z"][ext])) if ext.any() else None)
    with np.errstate(invalid="ignore"):
        land = ext & (pr["d"] > np.nanmax(pr["d"][cov]))
    span = np.nanmax(pr["d"][land]) - np.nanmax(pr["d"][cov])
    check("extension stops 10 m along the line past the DEM (landward)", 9.0 < span <= 10.0 + 1e-6, span)
    run_ = np.flatnonzero(np.isfinite(pr["z"]))
    check("no holes from the DEM's ends to the extensions' (the staircase line steps back at its turns)",
          np.isfinite(pr["z"][run_[0]:run_[-1] + 1]).all())
    check("beyond the extension: no elevation", np.isnan(pr["z"]).any() and np.isnan(pr["z"][-1]) and np.isnan(pr["z"][0]))
    # the line runs along the beach normal, so on the beach its length is the change in -z / beta
    i0 = np.flatnonzero(cov)[0]
    dd = np.abs((pr["d"] - pr["d"][i0]) - (zp - zp[i0]) / beta)
    check("d is distance along the line on the beach (within 0.4 m over 30 m)", np.nanmax(dd[cov | ext]) < 0.4 and
          np.nanmax(pr["d"][cov | ext]) - np.nanmin(pr["d"][cov | ext]) > 25, np.nanmax(dd[cov | ext]))
    s = pr["slope"][cov]
    check("DEM slope under the line is the beach's 0.10 (steepest descent)", abs(np.nanmedian(s) - beta) < 0.002,
          np.nanmedian(s))
    st = {"eta": 0.6, "sigma": 0.3}
    b, src = tw.foreshore_beta(pr, dem, st["eta"], st["sigma"])
    check("foreshore slope over <eta> +/- 2 sigma from the cells under the line", abs(b - beta) < 0.002 and
          src.startswith("dem under the line"), (b, src))
    b, src = tw.foreshore_beta({"kind": "none"}, None, 0.6, 0.3)
    check("no DEM, no --beta: no slope", np.isnan(b) and src == "none")
    b, src = tw.foreshore_beta({"kind": "planar"}, None, 0.6, 0.3, fixed=0.09)
    check("planar: the --beta slope", b == 0.09 and src == "fixed (--beta)")
    # the still-water plane stretches the line: why the edge is put on the beach
    import georectify
    E0, N0 = georectify.pixel_to_ground(u, v, 0.3, io, eo)
    dist = np.r_[0.0, np.nancumsum(np.hypot(np.diff(E0), np.diff(N0)))]
    edge = np.interp(np.linspace(0.0, 1.2, 200), pr["z"][cov], np.flatnonzero(cov))
    st = tw.plane_stretch(dist, pr["d"], edge)
    check("c2 at 230 m: 2-3 metres of path on the plane per beach metre over the swash", 2.0 < st < 3.0, st)
    i0, i1 = int(np.floor(np.percentile(edge, 5))), int(np.ceil(np.percentile(edge, 95)))
    chord = np.hypot(E0[i1] - E0[i0], N0[i1] - N0[i0]) / (pr["d"][i1] - pr["d"][i0])
    check("... of which the plane geometry (its straight-line span) is 1.1-1.8 x; the rest the staircase's zigzag",
          1.1 < chord < 1.8 and chord < st, (chord, st))
    check("plane_stretch: 1 where it cannot be measured", tw.plane_stretch(dist, np.full(len(dist), np.nan), edge) == 1.0)


def test_inputs():
    print("inputs: the DEM series, gaps in a DEM, the extension, a level record")
    io, eo = calibration()
    beta = 0.10
    u, v, zp = ground_line(io, eo, beta)
    u, v = u[::-1], v[::-1]
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        series = tmp / "dems"
        series.mkdir()
        write_dem(series / "dem_2026-10-17_7d_dem.asc", beta)          # window 11-17 Oct: centre 14 Oct 12:00
        src = tw.DemSource(str(series))
        near = tw._day_epoch("2026-10-20") + 12 * 3600.0                 # 6 days from the centre
        far = tw._day_epoch("2026-10-22") + 12 * 3600.0                  # 8 days
        check("DEM series: the DEM within 7 days of its window centre is used",
              src.for_epoch(near) is not None and src.path_for(near).name == "dem_2026-10-17_7d_dem.asc")
        check("DEM series: none more than 7 days away (setup_model.EvolvingSlope's rule, not a copy)",
              src.for_epoch(far) is None and src.path_for(far) is None)
        check("the DEM's band slope is the ForeshoreSlope EvolvingSlope already computed",
              src.for_epoch(near).foreshore() is src.ev.for_epoch(near)[0])
        # gaps inside the DEM along the line: bridged up to 6 m, not beyond
        for lo_z, hi_z, bridged in ((0.30, 0.50, True), (0.20, 1.10, False)):
            d = write_dem(tmp / "gap_dem.asc", beta)
            g, head = __import__("asc_to_geotiff").read_asc(d)
            g = np.where((g > lo_z) & (g < hi_z), -9999, np.nan_to_num(g, nan=-9999))
            with open(d, "w") as f:
                f.write("ncols 200\nnrows 200\nxllcorner 419900.0\nyllcorner 4638300.0\ncellsize 2.0\n"
                        "NODATA_value -9999\n")
                for row in g:
                    f.write(" ".join(f"{x:.4f}" for x in row) + "\n")
            pr = tw.line_profile(u, v, io, eo, dem=tw.Dem(d))
            fin = np.flatnonzero(np.isfinite(pr["z"]))
            inner = np.arange(fin[0], fin[-1] + 1)                  # the DEM's part and its extensions
            gap = inner[(pr["z"][inner] > lo_z + 0.05) & (pr["z"][inner] < hi_z - 0.05)]
            width = (hi_z - lo_z) / beta
            if bridged:
                check(f"a {width:.0f} m gap in the DEM along the line: bridged (straight), flagged extrapolated",
                      np.isfinite(pr["z"][inner]).all() and pr["extrap"][gap].all() and len(gap) > 5, len(gap))
            else:
                check(f"a {width:.0f} m gap: not bridged (no elevation inside it)",
                      (~np.isfinite(pr["z"][inner])).sum() > 20)
        # the extension carries on only a profile that rises landward
        d = tmp / "berm_dem.asc"
        x0, y0, cell, nx, ny = 419900.0, 4638300.0, 2.0, 200, 200
        XX, YY = np.meshgrid(x0 + (np.arange(nx) + 0.5) * cell, (y0 + (np.arange(ny) + 0.5) * cell)[::-1])
        xs = tw.seaward_distance(XX, YY)
        z = np.where(xs > -12.0, -beta * xs, 1.2 - 0.04 * (-12.0 - xs))      # a berm crest at +1.2 m, then down
        z[(z < -1.0) | (xs < -30.0)] = -9999
        with open(d, "w") as f:
            f.write(f"ncols {nx}\nnrows {ny}\nxllcorner {x0}\nyllcorner {y0}\ncellsize {cell}\nNODATA_value -9999\n")
            for row in z:
                f.write(" ".join(f"{x:.4f}" for x in row) + "\n")
        pr = tw.line_profile(u, v, io, eo, dem=tw.Dem(d))
        fin = np.flatnonzero(np.isfinite(pr["z"]))
        check("a profile falling landward (behind a berm) is not extended landward", not pr["extrap"][fin[-1]:].any()
              and not np.isfinite(pr["z"][fin[-1] + 1:]).any())
        # a water-level record with a gap
        lvf = tmp / "lv.csv"
        ep = [1000.0 + 600 * i for i in range(12)] + [1000.0 + 600 * 11 + 7200 + 600 * i for i in range(6)]
        lvf.write_text("epoch,level_navd88\n" + "".join(f"{e},{0.001 * e}\n" for e in ep))
        ls = tw.LevelSeries(lvf)
        lv, srcs = ls.at([1300.0, 1000.0 + 600 * 11 + 3600, 1000.0 + 600 * 11 + 7200 + 300, 0.0])
        check("level record: linear between readings, none across a gap over 1 h or outside it",
              abs(lv[0] - 1.3) < 1e-9 and np.isnan(lv[1]) and abs(lv[2] - 0.001 * (1000.0 + 600 * 11 + 7500)) < 1e-9
              and np.isnan(lv[3]) and srcs[0] == "level csv" and srcs[1] == "", (lv, srcs))
        owg = tmp / "owg.csv"
        owg.write_text("filename,epoch,time_utc,brightness,sharpness,saturation,status\n"
                       f"f,{EPOCH + 7200},x,160,4.1,0.0,blurred\nf,{EPOCH + 300},x,160,55,0.0,ok\n")
        ws = tw.WindowSharpness(owg)
        check("window sharpness: the bright frame of the same collection only (within 10 min)",
              ws.at(EPOCH) == 55.0 and np.isnan(ws.at(EPOCH + 5400)), (ws.at(EPOCH), ws.at(EPOCH + 5400)))


def test_classify():
    print("status labelling")
    c = tw.classify
    check("ok", c(120, 0.2, 60, 0.95, 50, True, 9.0) == "ok")
    check("dark (mean grey < 25)", c(20, 0.2, 60, 0.95, 50, True, 9.0) == "dark")
    check("low contrast from the edge contrast", c(120, 0.03, 60, 0.95, 50, True, 9.0) == "low contrast")
    check("low contrast from the window (owg sharpness < 10)", c(120, 0.2, 4.0, 0.95, 50, True, 9.0) == "low contrast")
    check("a sharp window and no contrast measure is not low contrast", c(120, np.nan, 60, 0.95, 50, True, 9.0) == "ok")
    check("no swash edge (valid fraction)", c(120, 0.2, 60, 0.6, 50, True, 9.0) == "no swash edge")
    check("no swash edge (too few uprushes)", c(120, 0.2, 60, 0.95, 10, True, 9.0) == "no swash edge")
    check("no slope (no profile)", c(120, 0.2, 60, 0.95, 50, False, 9.0) == "no slope")
    check("no period", c(120, 0.2, 60, 0.95, 50, True, np.nan) == "no period")
    check("the statuses are the documented ones",
          set(tw.STATUSES) == {"ok", "low contrast", "no swash edge", "no slope", "no period", "dark"})
    edge = 49.0 + 5 * np.sin(np.arange(300) / 3.0)          # uprushes every ~19 rows
    g = np.full((300, 100), 120.0)
    for k in range(300):
        g[k, :int(round(edge[k])) + 1] = 200.0
    check("edge contrast of a sharp tip: (200-120)/(200+120) = 0.25", abs(tw.edge_contrast(g, edge) - 0.25) < 1e-9,
          tw.edge_contrast(g, edge))
    import cv2
    gb = cv2.GaussianBlur(g.astype(np.float32), (0, 0), sigmaX=12, sigmaY=0.01) * 0.5 + 95
    check("blur and haze drive it below 0.05", tw.edge_contrast(gb, edge) < 0.05, tw.edge_contrast(gb, edge))
    check("too few uprushes for a contrast: NaN", np.isnan(tw.edge_contrast(g, np.full(300, np.nan))))


def test_periods():
    print("period")
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        cam = tmp / "wave_period_c2.csv"
        cam.write_text("epoch,time_utc,filename,line,tp_s,tm01_s,ig_fraction,peak_ratio,status\n"
                       "1000300,x,A.c2.ras.tiff,1,9.0,7,0.1,5,ok\n"
                       "1000300,x,A.c2.ras.tiff,2,11.0,7,0.1,5,ok\n"
                       "1000300,x,A.c2.ras.tiff,3,30.0,7,0.1,1,weak peak\n"
                       "2000300,x,B.c2.ras.tiff,1,12.0,7,0.1,1,weak peak\n")
        wm = tmp / "waves_marconi.csv"
        wm.write_text("time_utc,epoch,hs_best,tp_s\nx,2000000,1.0,8.5\nx,3000000,1.0,\n")
        by = tmp / "waves_44008.csv"
        by.write_text("time_utc,epoch,wvht_m,dpd_s,apd_s,mwd_deg\nx,3000000,1.0,7.7,6,90\n")
        p = tw.Periods(cam, wm, by, "auto")
        check("camera Tp when ok: the median of the stack's ok lines", p.for_stack("A.c2.ras.tiff", 1000300) == (10.0, "camera Tp", False),
              p.for_stack("A.c2.ras.tiff", 1000300))
        check("camera weak peak: the buoy (waves_marconi Tp)", p.for_stack("B.c2.ras.tiff", 2000300)[:2] ==
              (8.5, "buoy Tp (waves_marconi)"), p.for_stack("B.c2.ras.tiff", 2000300))
        r = p.for_stack("C.c2.ras.tiff", 3000300)
        check("no camera row yet: pending, and the buoy 44008 DPD meanwhile", r == (7.7, "buoy 44008 DPD", True), r)
        r = p.for_stack("B.c2.ras.tiff", 9000000)
        check("nothing within 1.5 h: no period", np.isnan(r[0]) and r[1] == "", r)
        r = tw.Periods(cam, wm, by, "buoy").for_stack("A.c2.ras.tiff", 1000300)
        check("--period-source buoy ignores the camera", np.isnan(r[0]), r)
        r = tw.Periods(cam, wm, by, "camera").for_stack("B.c2.ras.tiff", 2000300)
        check("--period-source camera never takes the buoy", np.isnan(r[0]), r)
        r = tw.Periods(tmp / "missing.csv", None, None, "auto").for_stack("A.c2.ras.tiff", 1000300)
        check("no camera period file at all: not pending", r[2] is False and np.isnan(r[0]), r)
        cam.write_text("epoch,time_utc,filename,line,tp_s,tm01_s,ig_fraction,peak_ratio,status\n"
                       "1000300,x,A.c2.ras.tiff,1,9.0,7,0.1,5,ok\n1000300,x,A.c2.ras.tiff,2,10.0,7,0.1,5,ok\n"
                       "1000300,x,A.c2.ras.tiff,3,14.0,7,0.1,5,ok\n")
        r = tw.Periods(cam, wm, by, "auto").for_stack("A.c2.ras.tiff", 1000300)
        check("three ok lines 9 / 10 / 14 s: their median 10 s (not the mean, 11)", r[0] == 10.0, r)


def test_csv():
    print("CSV")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "hs.csv"

        def row(name, line, epoch, hs):
            return {"epoch": epoch, "time_utc": "t", "filename": name, "line": line, "status": "ok", "hs_m": hs}
        n = tw.write_rows(out, [row("b", 1, 20, 1.0), row("a", 1, 10, 0.5), row("a", 2, 10, 0.6)])
        first = out.read_bytes()
        rows = list(csv.DictReader(open(out)))
        check("new file: header and rows in time order", n == 3 and [r["filename"] for r in rows] == ["a", "a", "b"]
              and list(rows[0].keys()) == tw.FIELDS)
        n = tw.write_rows(out, [row("a", 1, 10, 9.9), row("c", 1, 30, 2.0), row("c", 1, 30, 2.0)])
        rows = list(csv.DictReader(open(out)))
        check("append: rows already in the file are not added again, duplicates once", n == 1 and len(rows) == 4,
              (n, len(rows)))
        check("append: the existing rows are untouched, byte for byte", out.read_bytes().startswith(first))
        check("append: an existing row keeps its value", [r["hs_m"] for r in rows if r["filename"] == "a" and r["line"] == "1"] == ["0.5"])
        n = tw.write_rows(out, [row("a", 1, 10, 9.9)], reprocess=True)
        rows = list(csv.DictReader(open(out)))
        check("--reprocess replaces the redone row, keeps the others, no duplicate",
              n == 1 and len(rows) == 4 and [r["hs_m"] for r in rows if r["filename"] == "a" and r["line"] == "1"] == ["9.9"])
        check("--reprocess keeps time order", [int(r["epoch"]) for r in rows] == sorted(int(r["epoch"]) for r in rows))
        bad = Path(tmp) / "old.csv"
        bad.write_text("epoch,filename,line,hs_m\n1,a,1,0.5\n")
        try:
            tw.write_rows(bad, [row("z", 1, 5, 1.0)])
            refused = False
        except SystemExit:
            refused = True
        check("a file with other columns is refused (no shifted rows)", refused and bad.read_text().count("\n") == 2)
        check("nothing to write: no file", tw.write_rows(Path(tmp) / "none.csv", []) == 0 and not (Path(tmp) / "none.csv").exists())
        check("read_done: the (filename, line) keys", tw.read_done(out) == {("a", "1"), ("a", "2"), ("b", "1"), ("c", "1")})
        # run(): the header is checked before any stack is processed, and rows are appended as it goes
        ras = Path(tmp) / "ras"
        ras.mkdir()
        for i in range(5):
            (ras / f"{EPOCH + 600 * i}.c2.ras.tiff").write_bytes(b"")
        (Path(tmp) / "lv.csv").write_text("epoch,level_navd88\n1,0\n2,0\n")
        calls = []
        real = tw.estimate_stack

        seen = []

        def fake(p, ln, ctx, now=None):
            calls.append(p.name)
            o = Path(tmp) / "chunk.csv"
            seen.append(len(list(csv.DictReader(open(o)))) if o.exists() else 0)
            if len(calls) == 4:
                raise RuntimeError("a crash part way")
            e = int(p.name.split(".")[0])
            return {"epoch": e + 300, "time_utc": "t", "filename": p.name, "line": ln, "status": "ok", "hs_m": 1.0}
        argv = ["--ras-dir", str(ras), "--lines", "1", "--dem", str(Path(tmp) / "nodem"), "--level-csv",
                str(Path(tmp) / "lv.csv"), "--period-csv", str(Path(tmp) / "none.csv"), "--waves-csv",
                str(Path(tmp) / "none.csv"), "--buoy-csv", str(Path(tmp) / "none.csv"), "--owg-csv",
                str(Path(tmp) / "none.csv")]
        import io
        import contextlib
        tw.estimate_stack = fake
        chunk0 = tw.CHUNK_ROWS
        try:
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    tw.run(tw.build_parser().parse_args(argv + ["--output", str(bad)]))
                refused = False
            except SystemExit:
                refused = True
            check("run(): a file with other columns is refused BEFORE any stack is processed", refused and not calls,
                  calls)
            tw.CHUNK_ROWS = 2
            crashed = False
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    tw.run(tw.build_parser().parse_args(argv + ["--output", str(Path(tmp) / "chunk.csv")]))
            except RuntimeError:
                crashed = True
            rows = list(csv.DictReader(open(Path(tmp) / "chunk.csv")))
            check("run(): appended in chunks as it goes (2 rows on disk while the 3rd and 4th stacks run), and a "
                  "crash after 3 stacks keeps those 3", crashed and len(rows) == 3 and seen[2:4] == [2, 2], (len(rows), seen))
        finally:
            tw.estimate_stack = real
            tw.CHUNK_ROWS = chunk0


def test_end_to_end():
    print("end to end (synthetic c2 scene, real calibration)")
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = make_scene(tmp, 1.5, 10.0, 0.10)
        tr = truth_stats(sc["z_rel"])
        out = tmp / "hs.csv"
        rows = run_scene(sc, tmp, out, ["--wait-days", "0"])
        r = rows[0] if rows else {}
        check("one row, status ok, on the DEM", len(rows) == 1 and r.get("status") == "ok" and
              r.get("profile", "").startswith("dem:"), r.get("status"))
        if r.get("status") != "ok":
            return
        f = {k: float(r[k]) for k in ("r2_m", "setup_m", "s_ig_m", "h0_r2_m", "h0_sig_m", "h0_setup_m", "h0_m",
                                      "ks_21m", "hs_m", "beta_f", "edge_contrast", "period_s")}
        print(f"    truth: R2 {tr['r2']:.3f} setup {tr['setup']:.3f} S_IG {tr['sig']:.3f} m; read: R2 {f['r2_m']:.3f} "
              f"setup {f['setup_m']:.3f} S_IG {f['s_ig_m']:.3f}; H0 r2 {f['h0_r2_m']:.2f} sig {f['h0_sig_m']:.2f} "
              f"setup {f['h0_setup_m']:.2f} combined {f['h0_m']:.2f} (planted 1.50)")
        check("R2 within 3% of the swash's own", abs(f["r2_m"] / tr["r2"] - 1) < 0.03)
        check("H0 from R2 within 5% of the planted 1.5 m", abs(f["h0_r2_m"] / 1.5 - 1) < 0.05)
        check("setup within 0.12 m (the detector trails a fast tip: +0.06 m here)", abs(f["setup_m"] - tr["setup"]) < 0.12)
        check("S_IG within 25% of the swash's own", abs(f["s_ig_m"] / tr["sig"] - 1) < 0.25)
        check("combined H0 within 30% of the planted 1.5 m", abs(f["h0_m"] / 1.5 - 1) < 0.30)
        check("hs_m = H0 x Ks(T0, 21 m)", abs(f["hs_m"] - f["h0_m"] * f["ks_21m"]) < 0.002 and
              abs(f["ks_21m"] - float(tw.shoaling_coefficient(10.0, 21.0))) < 1e-4)
        check("the foreshore slope read off the DEM: 0.10", abs(f["beta_f"] - 0.10) < 0.003)
        check("clear scene: edge contrast well above 0.05", f["edge_contrast"] > 0.15)
        check("period from the camera file", r["period_source"] == "camera Tp" and f["period_s"] == 10.0)

        n_before = len(open(out).read().splitlines())
        rows = run_scene(sc, tmp, out, ["--wait-days", "0"])
        check("a second run adds nothing (no duplicate)", len(open(out).read().splitlines()) == n_before)
        rows = run_scene(sc, tmp, out, ["--wait-days", "0", "--reprocess"])
        check("--reprocess: still one row per stack and line", len(rows) == 1)

        # the window check: the bright frame of the same collection is blurred
        (tmp / "owg.csv").write_text("filename,epoch,time_utc,brightness,sharpness,saturation,status\n"
                                     f"{EPOCH}.c2.bright.jpg,{EPOCH},x,160,4.1,0.0,blurred (fog, rain or wet lens)\n")
        rows = run_scene(sc, tmp, out, ["--wait-days", "0", "--reprocess"])
        check("blurred window (owg sharpness 4.1): low contrast, numbers kept but labelled",
              rows[0]["status"] == "low contrast" and "window sharpness 4.1" in rows[0]["note"] and rows[0]["hs_m"] != "",
              rows[0]["note"])
        (tmp / "owg.csv").unlink()

        # no period anywhere: 'no period' once old enough, left for later while young
        rows = run_scene(sc, tmp, tmp / "np.csv", ["--wait-days", "0"], period=False)
        check("no period: status 'no period', no height", rows and rows[0]["status"] == "no period" and rows[0]["hs_m"] == "")
        rows = run_scene(sc, tmp, tmp / "np2.csv", ["--wait-days", "3"], now=EPOCH + 3600, period=False)
        check("no period, stack an hour old: not written yet (a later run will have it)", rows == [])
        # no DEM and no --beta: no slope; with --beta: the planar fallback, S_IG only
        rows = run_scene(sc, tmp, tmp / "ns.csv", ["--wait-days", "0", "--dem", str(tmp / "nodem")])
        check("no DEM, no --beta: 'no slope'", rows and rows[0]["status"] == "no slope")
        rows = run_scene(sc, tmp, tmp / "pl.csv", ["--wait-days", "0", "--dem", str(tmp / "nodem"), "--beta", "0.1"])
        check("--beta without a DEM: planar profile, only S_IG combined, no H0 from R2 or setup written",
              rows and rows[0]["status"] == "ok" and rows[0]["profile"].startswith("planar") and
              rows[0]["combined_from"] == "sig" and rows[0]["h0_r2_m"] == "" and rows[0]["h0_setup_m"] == ""
              and rows[0]["hs_sig_m"] == rows[0]["hs_m"], rows[0] if rows else "")
        # the camera period file has not reached the stack yet (pending): wait while young, then the buoy
        empty = tmp / "period_empty.csv"
        empty.write_text("epoch,time_utc,filename,line,tp_s,tm01_s,ig_fraction,peak_ratio,status\n")
        waves = tmp / "waves.csv"
        waves.write_text("time_utc,epoch,hs_best,tp_s\n" + "".join(f"x,{EPOCH - 3600 + 3600 * i},1.0,9.5\n"
                                                                     for i in range(4)))
        rows = run_scene(sc, tmp, tmp / "pe.csv", ["--wait-days", "3"], now=EPOCH + 3600, period=empty, waves=waves)
        check("camera period pending, stack an hour old: not written yet (the buoy's is there, the camera's to come)",
              rows == [])
        rows = run_scene(sc, tmp, tmp / "pe.csv", ["--wait-days", "3"], now=EPOCH + 4 * 86400, period=empty, waves=waves)
        check("... four days old: written with the buoy period", rows and rows[0]["period_source"] ==
              "buoy Tp (waves_marconi)" and rows[0]["period_s"] == "9.5", rows[0] if rows else "")
        # the GPS level arrives 1-2 days late: a young stack on the Chatham transfer waits for it
        rows = run_scene(sc, tmp, tmp / "gps.csv", ["--wait-days", "3"], now=EPOCH + 3600,
                         wl=FakeLevel(0.3, "chatham"))
        check("GPS in use, young stack with only the Chatham level: waits", rows == [])
        rows = run_scene(sc, tmp, tmp / "gps.csv", ["--wait-days", "3"], now=EPOCH + 3600, wl=FakeLevel(0.3, "gps"))
        check("... with the GPS level: written", rows and rows[0]["water_level_source"] == "gps")
        rows = run_scene(sc, tmp, tmp / "gps2.csv", ["--wait-days", "3"], now=EPOCH + 4 * 86400,
                         wl=FakeLevel(0.3, "chatham"))
        check("... four days old on the Chatham level: written with it", rows and rows[0]["water_level_source"] == "chatham")
        rows = run_scene(sc, tmp, tmp / "gps3.csv", ["--wait-days", "3"], now=EPOCH + 3600,
                         wl=FakeLevel(0.3, "chatham", gps=False))
        check("no GPS in use (--no-gps, --level-csv): nothing to wait for", len(rows) == 1)
        # --since / --until select stacks by their UTC day
        day = "2026-10-20"
        for opt, d, n_want in (("--since", "2026-10-21", 0), ("--until", "2026-10-19", 0), ("--since", day, 1),
                               ("--until", day, 1)):
            rows = run_scene(sc, tmp, tmp / f"sel{opt}{d}.csv", ["--wait-days", "0", opt, d])
            check(f"{opt} {d}: {n_want} row(s) (the stack is on {day})", len(rows) == n_want, len(rows))
        # a DEM that covers too little of the swash and no extension: no estimate survives -> 'no swash edge'
        narrow = write_dem(tmp / "narrow_dem.asc", 0.10, zlo=-1.0, zhi=0.35)
        rows = run_scene(sc, tmp, tmp / "nar.csv", ["--wait-days", "0", "--dem", str(narrow), "--max-extrapolate", "0"])
        check("the swash mostly off the profile: no R2, setup or S_IG -> 'no swash edge', no height",
              rows and rows[0]["status"] == "no swash edge" and rows[0]["hs_m"] == "" and
              float(rows[0]["off_profile_fraction"]) > 0.1, rows[0] if rows else "")

    # a still-water level rising 0.3 m through the burst (three times a spring tide's 10 minutes): the
    # level is interpolated across the burst, so R2 is still the swash's own
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = make_scene(tmp, 1.5, 10.0, 0.10, swl_rise=0.3)
        tr = truth_stats(sc["z_rel"])
        rows = run_scene(sc, tmp, tmp / "hs.csv", ["--wait-days", "0"])
        r = rows[0] if rows else {}
        ok = r.get("status") == "ok" and r.get("r2_m", "") != ""
        check("level rising 0.3 m through the burst: R2 within 1% of the swash's own (held at the burst middle it "
              "read +2%), setup within 0.12 m",
              ok and abs(float(r["r2_m"]) / tr["r2"] - 1) < 0.01 and abs(float(r["setup_m"]) - tr["setup"]) < 0.12,
              (r.get("r2_m"), tr["r2"], r.get("setup_m"), tr["setup"]))
        check("... and the level written is the burst middle's (0.45 m)", ok and abs(float(r["water_level_m"]) - 0.45) < 0.01,
              r.get("water_level_m"))

    # high water: the uprush reaches the line's landward 15% (find_edge's land reference) on a DEM that
    # covers the whole line: no R2 ('R2 above the line'), not a value capped by the line
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = make_scene(tmp, 1.5, 10.0, 0.10, swl=1.8)
        dem = write_dem(tmp / "high_dem.asc", 0.10, zlo=-1.0, zhi=4.6)
        rows = run_scene(sc, tmp, tmp / "hs.csv", ["--wait-days", "0", "--dem", str(dem)])
        r = rows[0] if rows else {}
        check("high water, the uprush into the line's landward 15%: 'R2 above the line', no R2, off-profile share "
              "written", r.get("status") == "ok" and r.get("r2_m") == "" and "R2 above the line" in r.get("note", "")
              and float(r.get("off_profile_fraction") or 0) > 0.01, r)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = make_scene(tmp, 1.5, 10.0, 0.10, blur=8.0, haze=0.5)
        rows = run_scene(sc, tmp, tmp / "hs.csv", ["--wait-days", "0"])
        check("blurred and hazy stack: 'low contrast' from its own edge contrast",
              rows and rows[0]["status"] == "low contrast" and float(rows[0]["edge_contrast"]) < 0.05,
              rows[0]["edge_contrast"] if rows else "")
        rows = run_scene(sc, tmp, tmp / "young.csv", ["--wait-days", "3", "--dem", str(tmp / "nodem")], now=EPOCH + 3600)
        check("a young blurred stack without its DEM waits too (its numbers are computed once the DEM is there)",
              rows == [])
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = make_scene(tmp, 1.5, 10.0, 0.10, dark=True)
        rows = run_scene(sc, tmp, tmp / "hs.csv", ["--wait-days", "0"])
        check("night stack: 'dark', nothing computed", rows and rows[0]["status"] == "dark" and rows[0]["hs_m"] == "")


CASES = ((0.8, 8.0, 0.12), (1.5, 10.0, 0.10), (3.0, 12.0, 0.09), (2.0, 9.0, 0.14))


def test_scene_cases():
    """The four synthetic scenes the module docstring quotes (THE DETECTOR, WHAT IT CANNOT DO)."""
    print("synthetic scenes: what the detector does to each statistic")
    print(f"    {'H0 T0 beta':<15} {'R2 read/own':>13} {'setup read-own':>15} {'S_IG read/own':>14} "
          f"{'H0: r2':>7} {'sig':>6} {'setup':>6} {'comb.':>6} {'plane/beach':>11}")
    for h0, t0, beta in CASES:
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sc = make_scene(tmp, h0, t0, beta)
            tr = truth_stats(sc["z_rel"])
            rows = run_scene(sc, tmp, tmp / "hs.csv", ["--wait-days", "0"])
            r = rows[0]
            ap = rft.build_parser()
            a = ap.parse_args([str(sc["ras"]), "--camera", "c2", "--pix", str(sc["pix"]), "--output", "-",
                               "--z-plane", "0.3"])
            sd_plane = float(np.nanstd(rft.swash_edge(sc["ras"] / sc["name"], a, {}, None)["pos"]))
        g = {k: (float(r[k]) if r[k] != "" else np.nan) for k in ("r2_m", "setup_m", "s_ig_m", "h0_r2_m", "h0_sig_m",
                                                                  "h0_setup_m", "h0_m")}
        ratio = sd_plane / float(np.std(-sc["z_rel"] / beta))
        print(f"    {h0:.1f} {t0:4.1f} {beta:.2f}   {g['r2_m'] / tr['r2'] - 1:>+12.1%} {g['setup_m'] - tr['setup']:>+14.3f} "
              f"{g['s_ig_m'] / tr['sig'] - 1:>+13.1%}  {g['h0_r2_m'] / h0 - 1:>+6.0%} {g['h0_sig_m'] / h0 - 1:>+6.0%} "
              f"{g['h0_setup_m'] / h0 - 1:>+6.0%} {g['h0_m'] / h0 - 1:>+6.0%} {ratio:>11.2f}"
              + ("   (" + r["note"] + ")" if r["note"] else ""))
        ok_r2 = np.isnan(g["r2_m"]) and "R2 above the profile" in r["note"] or abs(g["r2_m"] / tr["r2"] - 1) < 0.03
        check(f"H0 {h0} m: status ok, R2 within 3% (or flagged above the profile)", r["status"] == "ok" and ok_r2)
        check(f"H0 {h0} m: H0 from R2 within 5%, combined within 30%",
              (np.isnan(g["h0_r2_m"]) or abs(g["h0_r2_m"] / h0 - 1) < 0.05) and abs(g["h0_m"] / h0 - 1) < 0.30)
        check(f"H0 {h0} m: the still-water plane path stretches the swash 2-3 x", 2.0 < ratio < 3.0, ratio)
        if h0 >= 3.0:
            check(f"H0 {h0} m: part of the edge read on the extension and part off the profile, both written",
                  float(r["z_extrapolated_fraction"]) > 0 and float(r["off_profile_fraction"]) > 0,
                  (r["z_extrapolated_fraction"], r["off_profile_fraction"]))


def test_blur_series():
    print("edge contrast against blur (50% haze), and what the blur does to the heights (same swash)")
    vals, clear = [], None
    for blur, haze in ((0, 0.0), (2, 0.5), (4, 0.5), (6, 0.5), (10, 0.5)):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sc = make_scene(tmp, 1.5, 10.0, 0.10, blur=blur, haze=haze)
            r = run_scene(sc, tmp, tmp / "hs.csv", ["--wait-days", "0"])[0]
        vals.append(float(r["edge_contrast"]))
        hs, hsig = float(r["hs_m"]), float(r["hs_sig_m"])
        clear = clear or (hs, hsig)
        print(f"    blur {blur:>2} px, haze {haze:.0%}: edge contrast {r['edge_contrast']}, status {r['status']}, "
              f"hs_m {r['hs_m']} ({hs / clear[0] - 1:+.0%} on the clear stack's), hs_sig_m {r['hs_sig_m']} "
              f"({hsig / clear[1] - 1:+.0%})")
    check("contrast falls as the blur grows", all(a > b for a, b in zip(vals, vals[1:])), vals)
    check("clear 0.15 or more; 10 px blur under 0.05", vals[0] > 0.15 and vals[-1] < 0.05, vals)


def test_compare():
    print("--compare (against hand-computed scores)")
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        out = tmp / "hs.csv"
        day0 = 1792454400                                  # 2026-10-20 00:00 UTC
        # (day, hour, status, hs_sig, hs_r2, reference, owg, owg status): day 0 is before the split
        # (hours two apart: each has a decoy reference the hour after)
        plan = [(0, 12, "ok", 1.0, 1.1, 1.5, 1.4, "ok"), (0, 14, "ok", 1.2, 1.0, 1.0, 1.1, "ok"),
                (2, 12, "ok", 2.0, 2.6, 2.2, 2.5, "ok"), (2, 14, "ok", 2.6, 2.2, 2.4, 2.3, "ok"),
                (3, 12, "ok", 1.1, 1.3, 1.2, 1.0, "ok"), (3, 14, "ok", 0.9, 1.0, 0.8, 0.0, "blurred"),
                (4, 12, "ok", 2.3, 1.9, 2.1, 1.8, "ok"), (4, 16, "ok", 2.9, 2.2, 2.5, 2.7, "ok"),
                (5, 12, "ok", 1.4, 1.5, 1.6, 1.5, "ok"), (5, 14, "low contrast", 0.5, 0.6, 1.4, 1.3, "ok")]
        ks = 0.92
        rows, waves, owg = [], [], []
        for d, hr, st, hsig, hr2, ref, ow, ost in plan:
            hour = day0 + 86400 * d + 3600 * hr
            start = hour + 2400                                # hh:40: the NEXT hour is the nearer one
            rows.append({"epoch": start + 300, "time_utc": "t", "filename": f"{start}.c2.ras.tiff", "line": 1,
                         "status": st, "profile": "dem:x_dem.asc", "water_level_source": "gps" if hr % 4 else "chatham",
                         "period_source": "camera Tp", "h0_sig_m": hsig / ks, "h0_r2_m": hr2 / ks, "h0_setup_m": 1.0,
                         "h0_m": hsig / ks, "ks_21m": ks, "hs_m": hsig, "hs_sig_m": hsig})
            waves.append(f"x,{hour},1,{ref},1,{ref + 0.05},9\n")
            waves.append(f"x,{hour + 3600},1,{ref + 5.0},1,{ref + 5.0},9\n")   # a decoy the hour after
            owg.append(f"f,{start + 60},x,150,50,0,{ost},{ow}\n")
        tw.write_rows(out, rows)
        (tmp / "waves.csv").write_text("time_utc,epoch,hs_buoy,hs_buoy_marconi,hs_camera,hs_best,tp_s\n" + "".join(waves))
        (tmp / "owg.csv").write_text("filename,epoch,time_utc,brightness,sharpness,saturation,status,hs_m\n" + "".join(owg))
        args = tw.build_parser().parse_args(["--compare", "--output", str(out), "--waves-csv", str(tmp / "waves.csv"),
                                             "--owg-csv", str(tmp / "owg.csv"), "--adcp", str(tmp / "no_adcp.csv"),
                                             "--since", "2026-10-22"])
        import io
        import contextlib
        res, buf = {}, io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = tw.compare(args, res)
        txt = buf.getvalue()
        held = [x for x in plan if x[0] >= 2]
        clear = [x for x in held if x[2] == "ok"]
        est = np.array([x[3] for x in clear])
        ref = np.array([x[5] for x in clear])
        bu = res.get("hs_buoy_marconi", {})
        row = next((x for x in bu.get("alone", []) if x["col"] == "hs_sig" and x["subset"] == "all clear"), {})
        check("--compare --since: only the held-out stacks, against hs_buoy_marconi of the hour the stack STARTS in",
              rc == 0 and row.get("n") == len(clear) and abs(row["bias"] - (est - ref).mean()) < 1e-9 and
              abs(row["rmse"] - np.sqrt(np.mean((est - ref) ** 2))) < 1e-9, row)
        storm = [x for x in clear if x[5] >= 2.0]
        es, rs = np.array([x[3] for x in storm]), np.array([x[5] for x in storm])
        row = next((x for x in bu.get("alone", []) if x["col"] == "hs_sig" and x["subset"] == "clear, >= 2 m"), {})
        check("... and in storms (reference >= 2 m)", row.get("n") == len(storm) and
              abs(row["rmse"] - np.sqrt(np.mean((es - rs) ** 2))) < 1e-9, row)
        # with the OWG: only its ok frames, on the stacks all three have
        both = [x for x in clear if x[7] == "ok"]
        e1 = np.array([x[3] - x[5] for x in both])
        e2 = np.array([x[6] - x[5] for x in both])
        blend = np.array([0.5 * (x[3] + x[6]) - x[5] for x in both])
        for subset, sel in (("all clear", both), ("clear, >= 2 m", [x for x in both if x[5] >= 2.0])):
            f1 = np.array([x[3] - x[5] for x in sel])
            f2 = np.array([x[6] - x[5] for x in sel])
            fb = np.array([0.5 * (x[3] + x[6]) - x[5] for x in sel])
            w = next((x for x in bu.get("with_owg", []) if x["col"] == "hs_sig" and x["subset"] == subset), {})
            check(f"with the OWG ({subset}): the same {len(sel)} stacks for both, its ok frames only; RMSEs, error "
                  f"correlation and the 50:50 blend as computed by hand",
                  w.get("n") == len(sel) and abs(w["rmse"] - np.sqrt(np.mean(f1 ** 2))) < 1e-9 and
                  abs(w["owg_rmse"] - np.sqrt(np.mean(f2 ** 2))) < 1e-9 and
                  abs(w["owg_bias"] - f2.mean()) < 1e-9 and abs(w["blend_rmse"] - np.sqrt(np.mean(fb ** 2))) < 1e-9 and
                  abs(w["err_corr"] - np.corrcoef(f1, f2)[0, 1]) < 1e-9, (w, len(sel)))
        check("(the blend really is 50:50: here it beats both)", np.sqrt(np.mean(blend ** 2)) <
              min(np.sqrt(np.mean(e1 ** 2)), np.sqrt(np.mean(e2 ** 2))))
        c = bu.get("counts", {})
        check("rule 1 counts: held-out clear stacks, their days, storm stacks and separate storms (24 h apart)",
              c == {"stacks": len(clear), "days": 4, "storm_stacks": 4, "storms": 2}, c)
        r2row = next((x for x in bu.get("alone", []) if x["col"] == "hs_r2" and x["subset"] == "all clear"), {})
        check("R2, setup and the combined value are listed as not independent of the reference; S_IG is the candidate",
              "from R2 [not indep.]" == r2row.get("estimate") and "from S_IG (camera only)" in txt and
              "hs_m (combined) [not indep.]" in txt and "not independent of hs_buoy_marconi" in txt, r2row)
        sub = {x["subset"] for x in bu.get("alone", []) if x["col"] == "hs_sig"}
        check("split by water-level source, and blurred apart", {"clear, level gps", "clear, level chatham", "blurred"} <= sub,
              sub)
        check("hs_best: a consistency check only, scored too", "hs_best" in res and "CONTAINS the OWG" in txt)
        best = next((x for x in res["hs_best"]["alone"] if x["col"] == "hs_sig" and x["subset"] == "all clear"), {})
        check("hs_best of the start hour too (not the nearer next hour)",
              abs(best.get("bias", 9) - ((est - ref).mean() - 0.05)) < 1e-9, best)
        # a planar profile's R2 and setup are not scored
        rows[4]["profile"] = "planar:beta 0.1"
        out2 = tmp / "hs2.csv"
        tw.write_rows(out2, rows)
        args.output = str(out2)
        res2 = {}
        with contextlib.redirect_stdout(io.StringIO()):
            tw.compare(args, res2)
        r2b = next((x for x in res2["hs_buoy_marconi"]["alone"] if x["col"] == "hs_r2" and x["subset"] == "all clear"), {})
        check("a planar profile's R2 is left out of 'from R2'", r2b.get("n") == len(clear) - 1, r2b)


def test_runup_refactor():
    print("runup_from_timestack.py refactor")
    rng = np.random.default_rng(0)
    x = np.cumsum(rng.normal(0, 1, 500))
    x[rng.choice(500, 30, replace=False)] = np.nan
    a, b = rft.runup_peaks(x), x[rft.peak_mask(x)]
    check("runup_peaks == series[peak_mask]", len(a) == len(b) and bool(np.all((a == b) | (np.isnan(a) & np.isnan(b)))))
    y = np.sin(np.arange(1200) / 9.0)
    y[[5, 6, 300, 301, 302]] = np.nan
    f = rft.bridge_gaps(y)
    yy = y.copy()
    yy[200:215] = np.nan
    check("bridge_gaps: short gaps bridged linearly (band_sig's), None for a gap over 5 s; band_sig NaN then",
          f is not None and np.isfinite(f).all() and abs(f[6] - (y[4] + 2 * (y[7] - y[4]) / 3)) < 1e-12 and
          rft.bridge_gaps(yy) is None and np.isnan(rft.band_sig(yy, 0, 2)) and rft.LAND_REF_FRACTION == 0.85)
    ap = rft.build_parser()
    a = ap.parse_args(["x.ras.tiff", "--output", "o.csv"])
    check("build_parser keeps the defaults (k 4, min-run 15, retreat 3 m/s, valid 0.8, peaks 20)",
          (a.k, a.min_run, a.max_retreat_speed, a.min_valid, a.min_peaks, a.line) == (4.0, 15, 3.0, 0.8, 20, 1))
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = make_scene(tmp, 1.5, 10.0, 0.10)
        args = ap.parse_args([str(sc["ras"]), "--camera", "c2", "--pix", str(sc["pix"]), "--output", str(tmp / "r.csv")])
        s = rft.swash_edge(sc["ras"] / sc["name"], args, {}, None)
        row = rft.process(sc["ras"] / sc["name"], args, {}, None)
        check("process() reports what swash_edge found", row["mean_pos_m"] == round(float(np.nanmean(s["pos"])), 2) and
              row["valid_fraction"] == round(s["valid"], 3) and row["sea_at"] == s["sea_at"])
        sd_plane = float(np.nanstd(s["pos"]))
    sd_beach = float(np.std(-sc["z_rel"] / 0.10))
    check("the runup script's plane positions spread the swash 2-3 x wider than the beach does",
          2.0 < sd_plane / sd_beach < 3.0, (sd_plane, sd_beach))


def test_not_wired_in():
    print("opt-in")
    for f in ("waterline_timex_cron.sh", "owg.sh", "owg_report.py", "owg_live.py", "marconi_waves.py",
              "station_status.py", "runup_from_timestack.py"):
        txt = (HERE / f).read_text()
        bad = "timestack_wave_height" in txt and f != "runup_from_timestack.py"
        if f == "runup_from_timestack.py":
            bad = "import timestack_wave_height" in txt
        check(f"{f} does not run it", not bad)
    check("its own output file only", tw.build_parser().parse_args([]).output is None and
          "hs_timestack_" in Path(tw.__file__).read_text())


def test_syntax_38():
    print("Python 3.8")
    for f in ("timestack_wave_height.py", "test_timestack_wave_height.py", "runup_from_timestack.py"):
        try:
            ast.parse((HERE / f).read_text(), filename=f, feature_version=(3, 8))
            good = True
        except SyntaxError as exc:
            good = False
            print(f"    {exc}")
        check(f"{f} parses as Python 3.8", good)


def test_self_test():
    print("--self-test")
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = tw.self_test()
    check("Stockdon-generated swash: H0 recovered (exact 2%, noisy 5%)", rc == 0 and "PASS" in buf.getvalue(),
          buf.getvalue()[-300:])


def main():
    test_inversions()
    test_combination()
    test_spectrum_and_stats()
    test_ig_dof()
    test_line_end_and_off_profile()
    test_elevation_conversion()
    test_classify()
    test_periods()
    test_csv()
    test_inputs()
    test_end_to_end()
    test_scene_cases()
    test_blur_series()
    test_compare()
    test_runup_refactor()
    test_not_wired_in()
    test_syntax_38()
    test_self_test()
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        raise SystemExit(1)
    print("all timestack wave height tests passed")


if __name__ == "__main__":
    main()

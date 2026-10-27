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
    published Eq. 19 against runup_inversion.h0_from_runup;
  * the elevation conversion: the ray / DEM intersection on the real c2
    calibration (a planar beach and its DEM), the extension past the DEM,
    the fractional-pixel elevation of the edge, and why the still-water
    plane would not do (its stretch);
  * the swash statistics: the 0.05 Hz split, the detrend (a tide in the
    record), setup and R2;
  * the combination, the period choice and the error model option;
  * status labelling: dark, low contrast (edge contrast; the window's
    sharpness), no period, no slope, and the wait for late inputs;
  * the CSV: appended without duplicates, rewritten only by --reprocess,
    a changed header refused;
  * end to end on a synthetic camera scene (the image the real c2 camera
    would record of a planar beach): R2 and H0 from R2 recovered;
  * the refactor of runup_from_timestack.py (swash_edge, peak_mask,
    build_parser) leaves its outputs as they were;
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
               seed=5):
    """A synthetic c2 runup timestack of a planar beach: the swash elevation series (Stockdon-shaped,
    backwash at most 2 m/s on the beach), every pixel below it foam (bright, textured), above it sand.
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
    z_abs = swl + z_rel
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
    lv.write_text("epoch,level_navd88\n" + "".join(f"{e},{swl}\n" for e in range(epoch - 7200, epoch + 7200, 600)))
    per = tmp / "period.csv"
    per.write_text("epoch,time_utc,filename,line,tp_s,tm01_s,ig_fraction,peak_ratio,status\n"
                   f"{epoch + 300},x,{name},1,{t0},{t0 * 0.8},0.2,5,ok\n")
    return {"ras": ras, "pix": tmp / "c2_timestack.pix", "dem": dem, "level": lv, "period": per,
            "name": name, "z_rel": z_rel, "z_pix": zp}


def run_scene(sc, tmp, out, extra=(), now=None, period=True):
    argv = ["--ras-dir", str(sc["ras"]), "--pix", str(sc["pix"]), "--lines", "1",
            "--dem", str(sc["dem"]), "--level-csv", str(sc["level"]),
            "--period-csv", str(sc["period"] if period else Path(tmp) / "none.csv"),
            "--waves-csv", str(Path(tmp) / "none.csv"), "--buoy-csv", str(Path(tmp) / "none.csv"),
            "--owg-csv", str(Path(tmp) / "owg.csv"), "--output", str(out)] + list(extra)
    args = tw.build_parser().parse_args(argv)
    import io
    import contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        tw.run(args, now=now)
    return list(csv.DictReader(open(out))) if Path(out).exists() else []


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
    sc = tw.parse_scatter("r2=0.1:0.2,setup=0:0.5")
    check("--scatter parses relative:absolute per estimate", sc["r2"] == (0.1, 0.2) and sc["setup"] == (0.0, 0.5)
          and sc["sig"] == tw.SCATTER["sig"], sc)
    try:
        tw.parse_scatter("foo=1:2")
        bad = False
    except ValueError:
        bad = True
    check("--scatter refuses an unknown estimate", bad)
    check("the default level error is still water and the DEM's level, 0.28 m",
          abs(tw.LEVEL_SIGMA_M - np.hypot(0.12, 0.258)) < 1e-12 and abs(tw.LEVEL_SIGMA_M - 0.2845) < 0.001)


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
    check("c2 at 230 m: 2-3 plane metres per beach metre over the swash", 2.0 < st < 3.0, st)
    check("plane_stretch: 1 where it cannot be measured", tw.plane_stretch(dist, np.full(len(dist), np.nan), edge) == 1.0)


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
        check("--beta without a DEM: planar profile, only S_IG combined",
              rows and rows[0]["status"] == "ok" and rows[0]["profile"].startswith("planar") and
              rows[0]["combined_from"] == "sig", rows[0] if rows else "")

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = make_scene(tmp, 1.5, 10.0, 0.10, blur=8.0, haze=0.5)
        rows = run_scene(sc, tmp, tmp / "hs.csv", ["--wait-days", "0"])
        check("blurred and hazy stack: 'low contrast' from its own edge contrast",
              rows and rows[0]["status"] == "low contrast" and float(rows[0]["edge_contrast"]) < 0.05,
              rows[0]["edge_contrast"] if rows else "")
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
        check(f"H0 {h0} m: the still-water plane stretches the swash 2-3 x", 2.0 < ratio < 3.0, ratio)


def test_blur_series():
    print("edge contrast against blur (50% haze)")
    vals = []
    for blur, haze in ((0, 0.0), (2, 0.5), (4, 0.5), (6, 0.5), (10, 0.5)):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sc = make_scene(tmp, 1.5, 10.0, 0.10, blur=blur, haze=haze)
            r = run_scene(sc, tmp, tmp / "hs.csv", ["--wait-days", "0"])[0]
        vals.append(float(r["edge_contrast"]))
        print(f"    blur {blur:>2} px, haze {haze:.0%}: edge contrast {r['edge_contrast']}, status {r['status']}, "
              f"hs_m {r['hs_m']}")
    check("contrast falls as the blur grows", all(a > b for a, b in zip(vals, vals[1:])), vals)
    check("clear 0.15 or more; 10 px blur under 0.05", vals[0] > 0.15 and vals[-1] < 0.05, vals)


def test_compare_smoke():
    print("--compare (smoke)")
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        out = tmp / "hs.csv"
        rows = []
        h0 = 1792497600                                   # an hour boundary (2026-10-20 12:00 UTC)
        for i in range(12):
            e = h0 + 3600 * i + 600
            rows.append({"epoch": e + 300, "time_utc": "t", "filename": f"{e}.c2.ras.tiff", "line": 1,
                         "status": "ok" if i % 4 else "low contrast", "h0_r2_m": 1.0 + 0.1 * i,
                         "h0_sig_m": 1.1 + 0.1 * i, "h0_setup_m": 1.3, "h0_m": 1.1 + 0.1 * i, "ks_21m": 0.92,
                         "hs_m": round((1.1 + 0.1 * i) * 0.92, 3)})
        tw.write_rows(out, rows)
        (tmp / "waves.csv").write_text("time_utc,epoch,hs_buoy,hs_buoy_marconi,hs_camera,hs_best,tp_s\n" +
                                       "".join(f"x,{h0 + 3600 * i},1,{0.9 + 0.1 * i},1,{1.0 + 0.1 * i},9\n"
                                               for i in range(12)))
        (tmp / "owg.csv").write_text("filename,epoch,time_utc,brightness,sharpness,saturation,status,hs_m\n" +
                                     "".join(f"f,{h0 + 3600 * i + 600},x,150,50,0,ok,{1.05 + 0.1 * i}\n" for i in range(12)))
        args = tw.build_parser().parse_args(["--compare", "--output", str(out), "--waves-csv", str(tmp / "waves.csv"),
                                             "--owg-csv", str(tmp / "owg.csv"), "--adcp", str(HERE / "sig1000_waves_ALL.csv")])
        import io
        import contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = tw.compare(args)
        txt = buf.getvalue()
        check("--compare runs: no ADCP overlap in Oct 2026, so the converted buoy, hs_best and the OWG",
              rc == 0 and "hs_buoy_marconi" in txt and "hs_best" in txt and "OWG" in txt and "blurred" in txt, txt[-400:])
        check("--compare says hs_best is not independent", "CONTAINS the OWG" in txt)
        n_all = sum(1 for ln in txt.splitlines() if "all clear" in ln and ln.split()[-4] == "9")
        check("--compare scores the 9 clear stacks against each reference (buoy 5 rows, hs_best 4, OWG 4)",
              n_all == 13, n_all)
        check("--compare: the timestack and the OWG together (error correlation, 50:50 blend)",
              "error correlation" in txt and "50:50 blend" in txt)


def test_runup_refactor():
    print("runup_from_timestack.py refactor")
    rng = np.random.default_rng(0)
    x = np.cumsum(rng.normal(0, 1, 500))
    x[rng.choice(500, 30, replace=False)] = np.nan
    a, b = rft.runup_peaks(x), x[rft.peak_mask(x)]
    check("runup_peaks == series[peak_mask]", len(a) == len(b) and bool(np.all((a == b) | (np.isnan(a) & np.isnan(b)))))
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
    test_elevation_conversion()
    test_classify()
    test_periods()
    test_csv()
    test_end_to_end()
    test_scene_cases()
    test_blur_series()
    test_compare_smoke()
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

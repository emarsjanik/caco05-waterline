#!/usr/bin/env python3
"""
How Far Off Can A Georectified Pixel Be? (Ground Control, Left Out One At A Time)
===================================================================================
Measures the horizontal error the live camera calibration (EO 2025-11-13,
IO 2024-08-01) puts into every waterline point, from the targets that
calibration was solved from, and writes it as a table per camera that
dem_from_contours.py --sigma-obs turns into a per-cell uncertainty.

WHY. georectify.py --self-test projects the 2025-11-13 extrinsic targets
from their picked pixels and finds them 0.33 m (mean), 0.42 m (RMS) from
where they were surveyed. That is IN-SAMPLE: the EO was solved from those
same targets (c1 all six parameters, c2 its three angles: known finding
C-B6), so it says how well the solve fits its own points, not how well
it places a point it did not see. And the targets lie 29-150 m from the
cameras, while the waterlines run out to 200-400 m: there a pointing
error of a tenth of a degree moves the line metres, and nothing at the
targets shows it.

(a) LEAVE-ONE-OUT. The EO is solved again without each target in turn,
as the calibration itself was solved (CIRN C_singleExtrinsicSolution /
extrinsicsSolver.m: least squares on the picked-minus-projected PIXEL
residuals, the lens model applied; here Levenberg-Marquardt with numpy
only, as horizon_check.py does, the station's scipy being unusable), and
the target left out is projected from its picked pixel to its surveyed
height with the solve that never saw it. Its distance from the survey is
its leave-one-out error. Which parameters are solved: the most free
parameters whose solve on all the targets reproduces the calibration
(within FIT_POS_TOL_M and FIT_ANG_TOL_DEG): for c1 all six (the refit
moves 0.04 m, 0.04 deg), for c2 the three angles with the position held
(all six would move it 0.33 m: c2's solve held the position). --free
overrides. In-sample and leave-one-out are printed per target, with
their mean, RMS and worst, and kept in the table's header.

(b) MONTE CARLO. The uncertainty of the solved parameters as a
covariance, two ways: the JACOBIAN of the solve (s^2 (J'J)^-1, s^2 the
pixel residual variance; what CIRN's extrinsicsSolver reports as the
extrinsics error), which assumes the picks' errors are random and the
model right, and the JACKKNIFE of the leave-one-out solves ((n-1)/n x
their scatter; Efron & Tibshirani 1993, ch. 11), which measures how much
the solution depends on which targets it had. The table is made with the
one that gives the larger horizontal error over the waterlines' ranges
(out to RANGE_JUDGE_M), so it errs on the side of too much; --cov picks
one. MC_SAMPLES parameter sets are drawn from it (normal, about the
calibration's own values; position and angles that were not solved stay
fixed) and every pixel of the table is projected with each, through
georectify.pixel_to_ground, at each height of Z_LEVELS (the water level
plus setup a waterline is given spans about -1 to +2 m NAVD88). The
spread of where it lands is its horizontal error.

THE TABLE: calibration/<cam>_geom_sigma.csv. One row per pixel column of
GRID_COLUMNS (every 256 px across the frame), ground range of
GRID_RANGES_M (metres from the camera, 20-1000 m: pixel_to_ground keeps
nothing farther) and height of Z_LEVELS, wherever that column reaches
that range on that plane:
  pixel_row, easting/northing   where it is, with the calibration's EO
  sigma_x_m, sigma_y_m, corr_xy standard deviations of the easting and
                                northing (UTM x and y, as in the EO
                                files) and their correlation
  sigma_h_m                     sqrt(sigma_x^2 + sigma_y^2)
  sigma_cross_shore_m           along the station's seaward direction
                                (dem_figure.SHORE_NORMAL_DEG, 75.8 deg)
  footprint_cross_shore_m       the cross-shore size of the pixel on that
                                plane (its ground parallelogram projected
                                on the seaward direction): how finely the
                                camera can place a waterline there at all
                                (Holman & Stanley 2007; the literature's
                                footprint limit is 0.5-1 m)
  footprint_row_m, _column_m    the ground length of one pixel step down
                                and across the image
  mc_valid_share                share of the drawn pointings whose ray
                                still met the plane (sigma left blank
                                below MC_MIN_VALID)
'#' lines first: provenance (inputs by name and sha256, the solve, the
in-sample and leave-one-out numbers per target, the covariance used and
the parameter standard deviations, samples and seed). The numbers are
reproducible: the same inputs and seed give the same table.

WHAT IT IS NOT. It is the error of the calibration's SOLVE, as far as 13-14
targets 29-150 m away can show it: beyond them it is an extrapolation of
the fit, not a measurement, and it does not include a camera that has
turned since (pointing_check.py, horizon_check.py), the picks' origin
convention (VERIFIED L27), or the waterline detector's own pixel error (the
spread of repeat crossings, dem_from_contours.sigma_grid, has that). It is
for the 2025-11-13 calibration only: frames georectified with another EO
(the 2025 historical builds) need their own targets and their own table.

--evaluate (NUC): is a DEM's <stem>_sigma_obs.asc honest? Its cells
against a survey on NAVD88 (a lidar DSM or grid, .tif/.asc, sampled at
the cell centres; or RTK shots, .csv, each against the cell it falls in):
the share within 1 and 2 sigma_obs (about 68 and 95 % if it is right),
the same for the repeatability alone (<stem>_sigma.asc), both raw and
with the median offset removed (a datum error every cell shares), and the
RMS of the residuals in thirds of sigma_obs (it should rise with
sigma_obs if the geometry term means something). The survey must be of a
window georectified with the 2025-11-13 calibration (the 2025 lidars are
of other calibrations), and reach down into the waterlines: the 29 Sep
2026 RTK shots stop at +1.17 m and meet only 2 cells of that week's DEM,
so the test needs a low-tide survey. This tool never opens the 23 Oct
2024 ground control or survey (held for the blind check of calibration B,
gcpfree_eo.py): a file named for that day is refused.

Usage:
    python3 georectify_uncertainty.py                  # both cameras: report, write the tables
    python3 georectify_uncertainty.py --camera c1 --no-write
    python3 georectify_uncertainty.py --cov jacobian --samples 4000
    python3 dem_from_contours.py contour_points_ground.csv dem_test --sigma-obs --no-plot
    python3 georectify_uncertainty.py --evaluate dem_test --survey <survey>.tif|.asc|.csv
"""

import os
import re
import sys
import csv
import argparse
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from georectify import build_P, load_extrinsics, load_intrinsics, pixel_to_ground   # noqa: E402
from view_reproject import distort_uv                                              # noqa: E402

CAL = HERE / "calibration"
IO_DATE = "20240801"
EO_DATE = "20251113"
TARGETS = "2025-11-13_Marconi_Extrinsic_Targets_{cam}_xyz.csv"
PICKS = "gcp_uv_{cam}.csv"
EO_NAMES = ("x", "y", "z", "azimuth", "tilt", "roll")

# A solve "reproduces the calibration" when, run on all the targets from the calibration's own
# values, it moves the camera no more than this and turns it no more than this: c1's six-parameter
# solve moves 0.04 m and 0.04 deg, c2's 0.33 m and 0.16 deg (its position was held); c2's
# three-angle solve turns 0.08 deg.
FIT_POS_TOL_M = 0.10
FIT_ANG_TOL_DEG = 0.10
# Numerical derivatives of the pixel residuals: steps of 1 mm and 1 microradian.
FIT_STEP = (1e-3, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6)

MC_SAMPLES = 2000
MC_SEED = 0
MC_MIN_VALID = 0.95
GRID_COLUMN_STEP = 256
GRID_RANGES_M = tuple(range(20, 300, 10)) + tuple(range(300, 500, 25)) + tuple(range(500, 1001, 50))
Z_LEVELS = (-1.0, 0.0, 1.0, 2.0)
# The covariance that gives the larger median horizontal error at the nodes out to this range,
# where the waterlines are, makes the table (--cov larger).
RANGE_JUDGE_M = 400.0

# Any input named for 23 Oct 2024: Dr Sherwood's ground control of that day is held for the blind
# check of calibration B (gcpfree_eo.py), and this tool never opens it.
HELD_NAME = re.compile(r"2024-?10-?23", re.I)

COLUMNS = ["camera", "pixel_column", "range_m", "z_navd88", "pixel_row", "easting_utm19", "northing_utm19",
           "sigma_x_m", "sigma_y_m", "corr_xy", "sigma_h_m", "sigma_cross_shore_m", "footprint_cross_shore_m",
           "footprint_row_m", "footprint_column_m", "mc_valid_share"]


class Refused(Exception):
    """An input this tool does not open (the reason is printed)."""


def refuse_held(path):
    """Refuse a file named for 23 Oct 2024 (its own name or the name a link points to)."""
    p = Path(path)
    for name in (p.name, Path(os.path.realpath(str(p))).name):
        if HELD_NAME.search(name):
            raise Refused(f"{path}: named for 23 Oct 2024 -- that day's ground control and survey are held "
                          f"for the blind check of calibration B (gcpfree_eo.py); not opened")
    return p


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def git_head():
    """This checkout's commit, or 'unknown' (no git, or not a checkout)."""
    import subprocess
    try:
        r = subprocess.run(["git", "-C", str(HERE), "rev-parse", "--short", "HEAD"], stdout=subprocess.PIPE,
                           stderr=subprocess.DEVNULL, universal_newlines=True, timeout=10)
        return r.stdout.strip() or "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def load_targets(xyz_path, uv_path):
    """The surveyed targets (num, easting, northing, NAVD88 z; no header) and their picked pixels
    (num, U, V), read as georectify.self_test reads them. -> (nums, xyz (n x 3), uv (n x 2))."""
    xyz = {}
    with open(refuse_held(xyz_path), "r") as f:
        for row in csv.reader(f):
            if len(row) < 4:
                continue
            try:
                xyz[int(float(row[0]))] = (float(row[1]), float(row[2]), float(row[3]))
            except ValueError:
                continue
    uv = {}
    with open(refuse_held(uv_path), "r") as f:
        for row in csv.DictReader(f):
            uv[int(float(row["num"]))] = (float(row["U"]), float(row["V"]))
    nums = sorted(set(xyz) & set(uv))
    return nums, np.array([xyz[k] for k in nums], float).reshape(-1, 3), \
        np.array([uv[k] for k in nums], float).reshape(-1, 2)


def project(xyz, io, eo):
    """World points -> distorted pixels: the forward model the calibration was solved with (CIRN
    xyz2DistUV: P = K R [I | -C], then distortUV; georectify.build_P, view_reproject.distort_uv)."""
    P = build_P(io, eo)
    h = np.c_[xyz, np.ones(len(xyz))] @ P.T
    return distort_uv(h[:, 0] / h[:, 2], h[:, 1] / h[:, 2], io)


def fit_eo(io, eo0, xyz, uv, free, iterations=100):
    """
    The EO that best projects the targets onto their picked pixels:
    least squares on the (U, V) pixel residuals over the `free` parameters
    (indices into x y z azimuth tilt roll), starting from eo0, the others
    held -- as CIRN's extrinsicsSolver.m solves a calibration (nlinfit on
    xyz2DistUV). Levenberg-Marquardt with Marquardt's diagonal scaling (the
    pixels move ~40 px per metre of position and ~2600 per radian of angle
    here), numerical central-difference derivatives, numpy only.
    -> dict eo, resid (projected minus picked, all U then all V; 2n),
    J (2n x p, at the solution), converged
    """
    free = np.asarray(free, int)
    step = np.array(FIT_STEP)[free]

    def res(p):
        e = eo0.copy()
        e[free] = p
        U, V = project(xyz, io, e)
        return np.r_[U - uv[:, 0], V - uv[:, 1]]

    def jac(p):
        cols = []
        for k in range(len(p)):
            d = np.zeros(len(p))
            d[k] = step[k]
            cols.append((res(p + d) - res(p - d)) / (2.0 * step[k]))
        return np.column_stack(cols)

    p = eo0[free].astype(float).copy()
    r = res(p)
    cost = float(r @ r)
    lam = 1e-3
    converged = False
    for _ in range(iterations):
        J = jac(p)
        A, g = J.T @ J, J.T @ r
        moved = False
        while lam < 1e12:
            dp = -np.linalg.solve(A + lam * np.diag(np.diag(A)), g)
            r2 = res(p + dp)
            c2 = float(r2 @ r2)
            if np.isfinite(c2) and c2 <= cost:
                rel = (cost - c2) / max(cost, 1e-30)
                p, r, cost = p + dp, r2, c2
                lam = max(lam / 10.0, 1e-12)
                moved = True
                break
            lam *= 10.0
        if not moved or rel < 1e-14 or np.all(np.abs(dp) < step * 1e-3):
            converged = True
            break
    e = eo0.copy()
    e[free] = p
    return {"eo": e, "resid": r, "J": jac(p), "converged": converged}


def ground_errors(io, eo, xyz, uv):
    """(dE, dN) of each target projected from its picked pixel to its surveyed height."""
    E, N = pixel_to_ground(uv[:, 0], uv[:, 1], xyz[:, 2], io, eo)
    return E - xyz[:, 0], N - xyz[:, 1]


def choose_free(io, eo, xyz, uv, mode="auto"):
    """The parameters to solve (see the docstring): -> (free, label, note, refits)."""
    sets = (("all", [0, 1, 2, 3, 4, 5], "6 parameters (x y z azimuth tilt roll)"),
            ("angles", [3, 4, 5], "3 angles (azimuth tilt roll), the position held"))
    refits = {}
    for name, free, _ in sets:
        e = fit_eo(io, eo, xyz, uv, free)["eo"]
        refits[name] = (float(np.linalg.norm(e[:3] - eo[:3])), float(np.max(np.abs(np.rad2deg(e[3:] - eo[3:])))))
    if mode != "auto":
        name, free, label = next(s for s in sets if s[0] == mode)
        dpos, dang = refits[name]
        return free, label, f"chosen with --free {mode}; its refit on all targets moves {dpos:.3f} m, " \
                            f"{dang:.3f} deg", refits
    for name, free, label in sets:
        dpos, dang = refits[name]
        if dpos <= FIT_POS_TOL_M and dang <= FIT_ANG_TOL_DEG:
            other = "; ".join(f"{n} {refits[n][0]:.3f} m, {refits[n][1]:.3f} deg" for n, _, _ in sets if n != name)
            return free, label, (f"the most free parameters whose refit on all targets reproduces the calibration "
                                 f"(moves {dpos:.3f} m, {dang:.3f} deg; within {FIT_POS_TOL_M:g} m / "
                                 f"{FIT_ANG_TOL_DEG:g} deg; {other})"), refits
    name, free, label = sets[1]
    return free, label, (f"NOT REPRODUCED: no solve of these targets returns the calibration (6 parameters "
                         f"move {refits['all'][0]:.3f} m / {refits['all'][1]:.3f} deg, 3 angles "
                         f"{refits['angles'][1]:.3f} deg); the angles are used"), refits


def leave_one_out(io, eo, xyz, uv, free):
    """Each target left out of the solve in turn: -> (dE, dN, pixel error, solved parameters per left-out
    target (n x p))."""
    n = len(xyz)
    dE, dN, dpx, sols = np.zeros(n), np.zeros(n), np.zeros(n), []
    for k in range(n):
        keep = np.arange(n) != k
        e = fit_eo(io, eo, xyz[keep], uv[keep], free)["eo"]
        de, dn = ground_errors(io, e, xyz[k:k + 1], uv[k:k + 1])
        U, V = project(xyz[k:k + 1], io, e)
        dE[k], dN[k] = de[0], dn[0]
        dpx[k] = float(np.hypot(U[0] - uv[k, 0], V[0] - uv[k, 1]))
        sols.append(e[free])
    return dE, dN, dpx, np.array(sols)


def covariances(io, eo, xyz, uv, free, loo_sols):
    """-> {'jacobian': C, 'jackknife': C} over the free parameters (m, rad)."""
    f = fit_eo(io, eo, xyz, uv, free)
    J, r = f["J"], f["resid"]
    dof = max(len(r) - len(free), 1)
    s2 = float(r @ r) / dof
    cj = s2 * np.linalg.inv(J.T @ J)
    n = len(loo_sols)
    d = loo_sols - loo_sols.mean(axis=0)
    ck = (n - 1.0) / n * (d.T @ d)
    return {"jacobian": cj, "jackknife": ck}


def grid_columns(io):
    nu = int(io[0])
    return np.unique(np.r_[np.arange(0, nu, GRID_COLUMN_STEP), nu - 1]).astype(float)


def table_nodes(io, eo, columns=None, ranges=GRID_RANGES_M, zs=Z_LEVELS):
    """Pixel (column, row) of every (column, range, z) node the camera sees, with the calibration's EO:
    per column the rows from the bottom of the frame up, as far as the ground range keeps growing
    (a lens fold or the horizon ends it). -> arrays u, v, rng, z."""
    columns = grid_columns(io) if columns is None else np.asarray(columns, float)
    nv = int(io[1])
    rows = np.arange(nv - 1, -1, -1, dtype=float)          # bottom of the frame first
    U, V, R, Z = [], [], [], []
    for z in zs:
        for u in columns:
            E, N = pixel_to_ground(np.full(rows.shape, u), rows, z, io, eo, max_range_m=0)
            r = np.hypot(E - eo[0], N - eo[1])
            ok = np.isfinite(r)
            if not ok[0]:
                continue
            with np.errstate(invalid="ignore"):      # NaN past the horizon (numpy < 1.18 warns)
                grow = np.r_[True, np.diff(r) > 0] & ok
            stop = np.flatnonzero(~grow)
            last = stop[0] if len(stop) else len(r)
            rr, vv = r[:last], rows[:last]
            for R0 in ranges:
                if rr[0] <= R0 <= rr[-1]:
                    U.append(u)
                    V.append(float(np.interp(R0, rr, vv)))
                    R.append(float(R0))
                    Z.append(float(z))
    return np.array(U), np.array(V), np.array(R), np.array(Z)


def monte_carlo(io, eo, free, cov, u, v, z, samples=MC_SAMPLES, seed=MC_SEED):
    """Spread of the ground position of pixels (u, v) on planes z when the free EO parameters are drawn
    from N(eo, cov). -> (sxx, syy, sxy, valid share) per pixel."""
    free = np.asarray(free, int)
    try:
        L = np.linalg.cholesky(cov + np.eye(len(free)) * 1e-18)
    except np.linalg.LinAlgError:              # singular (a jackknife of few targets): its square root
        w, vec = np.linalg.eigh(cov)
        L = vec * np.sqrt(np.maximum(w, 0.0))
    E0, N0 = pixel_to_ground(u, v, z, io, eo, max_range_m=0)
    rs = np.random.RandomState(seed)
    n = np.zeros(len(u))
    s = np.zeros((5, len(u)))                                # sums of dx, dy, dx2, dy2, dxdy
    for _ in range(samples):
        e = eo.copy()
        e[free] = eo[free] + L @ rs.standard_normal(len(free))
        E, N = pixel_to_ground(u, v, z, io, e, max_range_m=0)
        dx, dy = E - E0, N - N0
        ok = np.isfinite(dx) & np.isfinite(dy)
        dx, dy = np.where(ok, dx, 0.0), np.where(ok, dy, 0.0)
        n += ok
        s += np.stack([dx, dy, dx * dx, dy * dy, dx * dy])
    with np.errstate(invalid="ignore", divide="ignore"):
        mx, my = s[0] / n, s[1] / n
        sxx = (s[2] - n * mx * mx) / (n - 1)
        syy = (s[3] - n * my * my) / (n - 1)
        sxy = (s[4] - n * mx * my) / (n - 1)
    share = n / float(samples)
    bad = share < MC_MIN_VALID
    for a in (sxx, syy, sxy):
        a[bad] = np.nan
    return sxx, syy, sxy, share


def shore_normal():
    """Unit vector (east, north) of the station's seaward direction (dem_figure.SHORE_NORMAL_DEG)."""
    from dem_figure import SHORE_NORMAL_DEG
    a = np.radians(SHORE_NORMAL_DEG)
    return np.array([np.sin(a), np.cos(a)]), SHORE_NORMAL_DEG


def footprints(io, eo, u, v, z, normal):
    """Ground size of each pixel on its plane: (cross-shore extent, row step length, column step length)."""
    def g(du, dv):
        E, N = pixel_to_ground(u + du, v + dv, z, io, eo, max_range_m=0)
        return np.stack([E, N])
    gu = g(0.5, 0.0) - g(-0.5, 0.0)
    gv = g(0.0, 0.5) - g(0.0, -0.5)
    cross = np.abs(normal @ gu) + np.abs(normal @ gv)
    return cross, np.hypot(*gv), np.hypot(*gu)


def summarise(err):
    err = np.asarray(err, float)
    return {"n": int(len(err)), "mean": float(err.mean()), "median": float(np.median(err)),
            "rms": float(np.sqrt(np.mean(err ** 2))), "max": float(err.max())}


def analyse_camera(cam, io_path, eo_path, xyz_path, uv_path, free_mode="auto", cov_mode="larger",
                   samples=MC_SAMPLES, seed=MC_SEED, columns=None, ranges=GRID_RANGES_M, zs=Z_LEVELS,
                   log=print):
    """Everything for one camera: the leave-one-out check and the table. -> (rows, header lines, result)."""
    io, eo = load_intrinsics(refuse_held(io_path)), load_extrinsics(refuse_held(eo_path))
    nums, xyz, uv = load_targets(xyz_path, uv_path)
    if len(nums) < 5:
        raise ValueError(f"{cam}: only {len(nums)} target(s) in both {xyz_path} and {uv_path}")
    rng_t = np.hypot(xyz[:, 0] - eo[0], xyz[:, 1] - eo[1])
    free, solve, why, _ = choose_free(io, eo, xyz, uv, free_mode)
    dE, dN = ground_errors(io, eo, xyz, uv)
    U, V = project(xyz, io, eo)
    ins = np.hypot(dE, dN)
    ins_px = np.hypot(U - uv[:, 0], V - uv[:, 1])
    lE, lN, lpx, sols = leave_one_out(io, eo, xyz, uv, free)
    loo = np.hypot(lE, lN)
    si, sl = summarise(ins), summarise(loo)
    log("=" * 78)
    log(f"{cam}: {Path(eo_path).name} + {Path(io_path).name}, {len(nums)} targets at {rng_t.min():.0f}-"
        f"{rng_t.max():.0f} m from the camera")
    log(f"  solve: {solve} -- {why}")
    log(f"  {'target':>6} {'range m':>8} {'in-sample m':>12} {'LOO m':>8} {'LOO dE':>8} {'LOO dN':>8} "
        f"{'in px':>7} {'LOO px':>7}")
    for k in range(len(nums)):
        log(f"  {nums[k]:>6} {rng_t[k]:8.1f} {ins[k]:12.3f} {loo[k]:8.3f} {lE[k]:+8.3f} {lN[k]:+8.3f} "
            f"{ins_px[k]:7.1f} {lpx[k]:7.1f}")
    worst = int(np.argmax(loo))
    log(f"  in-sample (the calibration's EO): mean {si['mean']:.3f} m, RMS {si['rms']:.3f} m, max {si['max']:.3f} m; "
        f"pixel RMS {np.sqrt(np.mean(ins_px ** 2)):.1f} px")
    log(f"  leave-one-out                   : mean {sl['mean']:.3f} m, RMS {sl['rms']:.3f} m, max {sl['max']:.3f} m "
        f"(target {nums[worst]}, {rng_t[worst]:.0f} m out); RMS {sl['rms'] / si['rms']:.2f} x in-sample")

    covs = covariances(io, eo, xyz, uv, free, sols)
    u, v, rr, zz = table_nodes(io, eo, columns, ranges, zs)
    normal, normal_deg = shore_normal()
    res = {}
    for name, c in covs.items():
        sxx, syy, sxy, share = monte_carlo(io, eo, free, c, u, v, zz, samples, seed)
        sh = np.sqrt(sxx + syy)
        judge = np.isfinite(sh) & (rr <= RANGE_JUDGE_M)
        res[name] = (sxx, syy, sxy, share, float(np.median(sh[judge])) if judge.any() else np.nan)
    if cov_mode == "larger":
        use = max(res, key=lambda k: res[k][4] if np.isfinite(res[k][4]) else -1.0)
        cov_why = (f"the larger of the two over the nodes out to {RANGE_JUDGE_M:g} m (median horizontal sigma: "
                   + ", ".join(f"{k} {res[k][4]:.3f} m" for k in sorted(res)) + ")")
    else:
        use = cov_mode
        cov_why = f"chosen with --cov {cov_mode} (median horizontal sigma out to {RANGE_JUDGE_M:g} m: " \
                  + ", ".join(f"{k} {res[k][4]:.3f} m" for k in sorted(res)) + ")"
    sxx, syy, sxy, share, _ = res[use]
    sd = np.sqrt(np.diag(covs[use]))
    sd_txt = ", ".join(f"{EO_NAMES[i]} {sd[j]:.3f} m" if i < 3 else f"{EO_NAMES[i]} {np.degrees(sd[j]):.3f} deg"
                       for j, i in enumerate(free))
    log(f"  EO covariance: {use} -- {cov_why}")
    log(f"    standard deviations: {sd_txt}")
    E0, N0 = pixel_to_ground(u, v, zz, io, eo, max_range_m=0)
    fc, fr, fcol = footprints(io, eo, u, v, zz, normal)
    with np.errstate(invalid="ignore"):
        sh = np.sqrt(sxx + syy)
        corr = sxy / np.sqrt(sxx * syy)
        scs = np.sqrt(normal[0] ** 2 * sxx + normal[1] ** 2 * syy + 2 * normal[0] * normal[1] * sxy)
    z_show = zs[min(range(len(zs)), key=lambda i: abs(zs[i]))]
    other = [k for k in sorted(res) if k != use]
    log(f"  on the plane z = {z_show:g} m, median over the columns: sigma_h and cross-shore sigma ({use}), sigma_h "
        f"with the {' / '.join(other)} covariance, cross-shore pixel footprint")
    for R0 in (50, 100, 150, 200, 300, 400):
        m = (rr == R0) & (zz == z_show) & np.isfinite(sh)
        if m.any():
            alt = " / ".join(f"{np.nanmedian(np.sqrt(res[k][0][m] + res[k][1][m])):6.2f} m" for k in other)
            log(f"    {R0:4d} m: {np.median(sh[m]):6.2f} m, {np.median(scs[m]):6.2f} m; {alt}; "
                f"{np.median(fc[m]):5.2f} m  ({int(m.sum())} column(s))")
    rows = []
    fmt = lambda x, p=4: "" if not np.isfinite(x) else f"{x:.{p}f}"
    for k in np.lexsort((rr, u, zz)):
        rows.append([cam, f"{u[k]:.0f}", f"{rr[k]:g}", f"{zz[k]:g}", fmt(v[k], 2), fmt(E0[k], 3), fmt(N0[k], 3),
                     fmt(np.sqrt(sxx[k])), fmt(np.sqrt(syy[k])), fmt(corr[k], 3), fmt(sh[k]), fmt(scs[k]),
                     fmt(fc[k]), fmt(fr[k]), fmt(fcol[k]), f"{share[k]:.3f}"])
    files = (("io", io_path), ("eo", eo_path), ("targets_xyz", xyz_path), ("targets_uv", uv_path))
    head = [f"georectify_uncertainty.py: horizontal error of each pixel of camera {cam} (calibration/README.md)",
            f"written_utc: {datetime.now(timezone.utc):%Y-%m-%dT%H:%M:%SZ}",
            f"git_head: {git_head()}",
            f"script_sha256: {sha256(Path(__file__).resolve())}"]
    head += [f"{k}: {Path(p).name} sha256 {sha256(p)}" for k, p in files]
    head += [f"camera_x_utm19: {eo[0]:.4f}", f"camera_y_utm19: {eo[1]:.4f}", f"camera_z_navd88: {eo[2]:.4f}",
             f"solve: {solve} -- {why}",
             f"in_sample: n {si['n']}, ground mean {si['mean']:.3f} m, RMS {si['rms']:.3f} m, max {si['max']:.3f} m; "
             f"pixel RMS {np.sqrt(np.mean(ins_px ** 2)):.2f} px (the calibration's EO, as georectify.py --self-test)",
             f"leave_one_out: ground mean {sl['mean']:.3f} m, RMS {sl['rms']:.3f} m, max {sl['max']:.3f} m (target "
             f"{nums[worst]} at {rng_t[worst]:.0f} m); RMS {sl['rms'] / si['rms']:.2f} x in-sample"]
    head += [f"loo_target: {nums[k]} range {rng_t[k]:.1f} m, in-sample {ins[k]:.3f} m, LOO {loo[k]:.3f} m "
             f"(dE {lE[k]:+.3f}, dN {lN[k]:+.3f})" for k in range(len(nums))]
    head += [f"covariance: {use} -- {cov_why}", f"eo_sd: {sd_txt}",
             f"monte_carlo: {samples} samples, seed {seed}, normal about the calibration's EO; sigma blank where "
             f"fewer than {MC_MIN_VALID:g} of them meet the plane",
             f"shore_normal_deg: {normal_deg:g} (dem_figure.SHORE_NORMAL_DEG)",
             "grid_columns: " + " ".join(f"{c:.0f}" for c in (grid_columns(io) if columns is None else columns)),
             "grid_ranges_m: " + " ".join(f"{r:g}" for r in ranges),
             "grid_z_navd88: " + " ".join(f"{z:g}" for z in zs),
             "units: metres; sigma_x/sigma_y = standard deviations of the easting/northing (UTM 19N x/y) where the "
             "pixel lands on the plane z, under the EO covariance; read with georectify_uncertainty.read_geom_sigma"]
    result = {"camera": cam, "in_sample": si, "loo": sl, "solve": solve, "free": list(free), "covariance": use,
              "cov_median_sigma_h": {k: res[k][4] for k in res}, "eo_sd": sd, "eo": eo, "io": io,
              "nums": nums, "loo_errors": loo, "in_errors": ins}
    return rows, head, result


def write_table(path, rows, head):
    with open(path, "w", newline="") as f:
        for line in head:
            f.write(f"# {line}\n")
        w = csv.writer(f, lineterminator="\n")
        w.writerow(COLUMNS)
        w.writerows(rows)


def read_geom_sigma(path):
    """
    A table written by this script -> dict: camera, x, y (the camera's
    position), columns, ranges, zs (the grid), sxx, syy, sxy (arrays z x
    range x column, m^2, NaN where the camera does not see that node),
    provenance ({key: text} from the '#' lines), sha256, path.
    """
    prov, lines = {}, []
    with open(path, newline="") as f:
        for line in f:
            if line.startswith("#"):
                k, _, val = line[1:].strip().partition(": ")
                if k and k not in prov:
                    prov[k] = val
            else:
                lines.append(line)
    rows = list(csv.DictReader(lines))
    num = lambda key: np.array([float(r[key]) if r[key] != "" else np.nan for r in rows], float)
    cols = np.array([float(x) for x in prov["grid_columns"].split()])
    rngs = np.array([float(x) for x in prov["grid_ranges_m"].split()])
    zs = np.array([float(x) for x in prov["grid_z_navd88"].split()])
    shape = (len(zs), len(rngs), len(cols))
    out = {k: np.full(shape, np.nan) for k in ("sxx", "syy", "sxy")}
    if rows:
        iz = np.searchsorted(zs, num("z_navd88"))
        ir = np.searchsorted(rngs, num("range_m"))
        ic = np.searchsorted(cols, num("pixel_column"))
        sx, sy, c = num("sigma_x_m"), num("sigma_y_m"), num("corr_xy")
        out["sxx"][iz, ir, ic] = sx * sx
        out["syy"][iz, ir, ic] = sy * sy
        out["sxy"][iz, ir, ic] = c * sx * sy
    out.update(camera=rows[0]["camera"] if rows else None, x=float(prov["camera_x_utm19"]),
               y=float(prov["camera_y_utm19"]), columns=cols, ranges=rngs, zs=zs, provenance=prov,
               sha256=sha256(path), path=str(path))
    return out


def geom_cov(table, u, rng, z):
    """
    The horizontal error covariance (sxx, syy, sxy, m^2) of points at
    pixel column u, ground range rng (m from the camera) and height z, by
    trilinear interpolation in the table (a convex mix of covariances is a
    covariance). Coordinates outside the grid are clamped to its edge (z:
    the plane nearest; range and column: the last node), and only the
    corner nodes the camera sees take part (their weights renormalised),
    so a point at the edge of the view takes the nearest nodes there. NaN
    where no corner is seen, or u, rng or z is not finite.
    """
    u, rng, z = (np.asarray(a, float) for a in np.broadcast_arrays(u, rng, z))
    out = [np.zeros(u.shape) for _ in range(3)]
    wsum = np.zeros(u.shape)
    good = np.isfinite(u) & np.isfinite(rng) & np.isfinite(z)
    axes = []
    for grid, x in ((table["zs"], z), (table["ranges"], rng), (table["columns"], u)):
        xc = np.clip(np.where(good, x, grid[0]), grid[0], grid[-1])
        i0 = np.clip(np.searchsorted(grid, xc, side="right") - 1, 0, max(len(grid) - 2, 0))
        i1 = np.minimum(i0 + 1, len(grid) - 1)
        span = np.where(i1 > i0, grid[i1] - grid[i0], 1.0)
        t = np.where(i1 > i0, (xc - grid[i0]) / span, 0.0)
        axes.append((i0, i1, t))
    (z0, z1, tz), (r0, r1, tr), (c0, c1, tc) = axes
    for iz, wz in ((z0, 1 - tz), (z1, tz)):
        for ir, wr in ((r0, 1 - tr), (r1, tr)):
            for ic, wc in ((c0, 1 - tc), (c1, tc)):
                w = wz * wr * wc
                vals = [table[k][iz, ir, ic] for k in ("sxx", "syy", "sxy")]
                seen = np.isfinite(vals[0]) & np.isfinite(vals[1]) & np.isfinite(vals[2]) & (w > 0)
                wsum += np.where(seen, w, 0.0)
                for o, val in zip(out, vals):
                    o += np.where(seen, w * np.where(seen, val, 0.0), 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        res = [np.where(good & (wsum > 0), o / wsum, np.nan) for o in out]
    return res[0], res[1], res[2]


def evaluate(stem, survey_path, log=print):
    """
    --evaluate: a DEM's <stem>_sigma_obs.asc against a survey (see the
    docstring). A grid (.tif or .asc, compare_dem_survey.read_survey): the
    cells with a DEM value, both sigmas and a survey value under their
    centre (bilinear, compare_dem_survey.sample). A point survey (.csv:
    Emlid Flow, CIRN num,E,N,Z or E/N/Z headers, survey_compare.read_points,
    e.g. the 29 Sep 2026 RTK check shots): each shot against the cell it
    falls in. -> record.
    """
    refuse_held(survey_path)
    from compare_dem_survey import read_survey, sample, read_asc
    grids = {}
    for k in ("dem", "sigma", "sigma_obs"):
        p = Path(f"{stem}_{k}.asc")
        if not p.exists():
            raise FileNotFoundError(f"{p} not found" + (" (build the DEM with dem_from_contours.py --sigma-obs)"
                                                        if k == "sigma_obs" else ""))
        grids[k] = read_asc(p)
    dem, h = grids["dem"]
    for k in ("sigma", "sigma_obs"):
        if grids[k][0].shape != dem.shape:
            raise ValueError(f"{stem}_{k}.asc is not on the grid of {stem}_dem.asc")
    cell = h["cellsize"]
    nr, nc = dem.shape
    se, so = grids["sigma"][0], grids["sigma_obs"][0]
    if Path(survey_path).suffix.lower() == ".csv":
        from survey_compare import read_points
        pts = read_points(survey_path)
        col = np.floor((pts["E"] - h["xllcorner"]) / cell).astype(int)
        row = np.floor((h["yllcorner"] + nr * cell - pts["N"]) / cell).astype(int)   # rows north to south
        inside = (col >= 0) & (col < nc) & (row >= 0) & (row < nr)
        row, col = row[inside], col[inside]
        dem, se, so, sv = dem[row, col], se[row, col], so[row, col], pts["Z"][inside]
        unit = "shot"
    else:
        cc, rr = np.meshgrid(np.arange(nc), np.arange(nr))
        E = h["xllcorner"] + (cc + 0.5) * cell
        N = h["yllcorner"] + (nr - rr - 0.5) * cell                 # rows north to south
        sg, x0, y0, scell = read_survey(survey_path)
        sv = sample(sg, x0, y0, scell, E.ravel(), N.ravel()).reshape(dem.shape)
        unit = "cell"
    with np.errstate(invalid="ignore"):          # NaN cells (numpy < 1.18 warns)
        ok = np.isfinite(dem) & np.isfinite(sv) & np.isfinite(se) & np.isfinite(so) & (so > 0)
    if ok.sum() < 10:
        log(f"Only {int(ok.sum())} {unit}(s) with a DEM value, both sigmas and the survey: nothing to judge.")
        return {"cells": int(ok.sum()), "unit": unit}
    r, se, so = (dem - sv)[ok], se[ok], so[ok]
    off = float(np.median(r))
    rec = {"cells": int(ok.sum()), "unit": unit, "median_m": off, "rmse_m": float(np.sqrt(np.mean(r ** 2))),
           "nmad_m": float(1.4826 * np.median(np.abs(r - off)))}
    log(f"{stem}: {rec['cells']} {unit}(s) against {Path(survey_path).name}: DEM - survey median {off:+.3f} m, "
        f"RMSE {rec['rmse_m']:.3f} m, NMAD {rec['nmad_m']:.3f} m")
    for label, x in (("raw", r), ("median offset removed", r - off)):
        for name, s in (("sigma_obs", so), ("sigma (SE)", se)):
            a1, a2 = float(np.mean(np.abs(x) <= s)), float(np.mean(np.abs(x) <= 2 * s))
            zr = float(np.sqrt(np.mean((x / np.maximum(s, 1e-6)) ** 2)))
            rec[f"{label}|{name}"] = {"within_1": a1, "within_2": a2, "rms_z": zr}
            log(f"  {label:<22} {name:<10}: within 1 sigma {100 * a1:5.1f} %, within 2 {100 * a2:5.1f} % "
                f"(about 68 / 95 if right), RMS of residual/sigma {zr:.2f} (about 1)")
    edges = np.percentile(so, [0, 100 / 3.0, 200 / 3.0, 100])
    log("  RMS of (residual - median offset) in thirds of sigma_obs (it should rise with sigma_obs):")
    rec["thirds"] = []
    for k in range(3):
        m = (so >= edges[k]) & ((so <= edges[k + 1]) if k == 2 else (so < edges[k + 1]))
        if m.any():
            t = {"sigma_obs_median_m": float(np.median(so[m])), "rms_m": float(np.sqrt(np.mean((r[m] - off) ** 2))),
                 "cells": int(m.sum())}
            rec["thirds"].append(t)
            log(f"    sigma_obs {edges[k]:.3f}-{edges[k + 1]:.3f} m (median {t['sigma_obs_median_m']:.3f}): "
                f"RMS {t['rms_m']:.3f} m over {t['cells']} {unit}(s)")
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--camera", default="both", choices=["c1", "c2", "both"])
    ap.add_argument("--io", default=None, help=f"default calibration/CACO05_<cam>_{IO_DATE}_IO.yaml (one camera)")
    ap.add_argument("--eo", default=None, help=f"default calibration/CACO05_<cam>_{EO_DATE}_EO-CV.yaml (one camera)")
    ap.add_argument("--targets-xyz", default=None, help=f"default calibration/{TARGETS.format(cam='<cam>')}")
    ap.add_argument("--targets-uv", default=None, help=f"default calibration/{PICKS.format(cam='<cam>')}")
    ap.add_argument("--free", default="auto", choices=["auto", "all", "angles"],
                    help="parameters solved: auto (default: the most free parameters whose solve reproduces the "
                         "calibration), all six, or the three angles with the position held")
    ap.add_argument("--cov", default="larger", choices=["larger", "jackknife", "jacobian"],
                    help="EO covariance for the Monte Carlo (default: the larger of the two)")
    ap.add_argument("--samples", type=int, default=MC_SAMPLES)
    ap.add_argument("--seed", type=int, default=MC_SEED)
    ap.add_argument("--output-dir", default=str(CAL), help="where <cam>_geom_sigma.csv goes (default calibration/)")
    ap.add_argument("--no-write", action="store_true", help="report only, write no table")
    ap.add_argument("--evaluate", metavar="STEM", default=None,
                    help="judge <STEM>_sigma_obs.asc (dem_from_contours.py --sigma-obs) against --survey")
    ap.add_argument("--survey", default=None, help="with --evaluate: lidar DSM or RTK grid (.tif or .asc), NAVD88")
    args = ap.parse_args()
    try:
        if args.evaluate:
            if not args.survey:
                ap.error("--evaluate needs --survey")
            evaluate(args.evaluate, args.survey)
            return 0
        cams = ["c1", "c2"] if args.camera == "both" else [args.camera]
        if len(cams) > 1 and (args.io or args.eo or args.targets_xyz or args.targets_uv):
            ap.error("--io, --eo, --targets-xyz and --targets-uv are for one camera: add --camera")
        for cam in cams:
            io = args.io or CAL / f"CACO05_{cam}_{IO_DATE}_IO.yaml"
            eo = args.eo or CAL / f"CACO05_{cam}_{EO_DATE}_EO-CV.yaml"
            xyz = args.targets_xyz or CAL / TARGETS.format(cam=cam)
            uv = args.targets_uv or CAL / PICKS.format(cam=cam)
            rows, head, _ = analyse_camera(cam, io, eo, xyz, uv, args.free, args.cov, args.samples, args.seed)
            if not args.no_write:
                out = Path(args.output_dir) / f"{cam}_geom_sigma.csv"
                write_table(out, rows, head)
                print(f"  wrote {out} ({len(rows)} nodes)")
    except Refused as exc:
        print(f"REFUSED: {exc}")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

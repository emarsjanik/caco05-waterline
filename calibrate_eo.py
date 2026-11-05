#!/usr/bin/env python3
"""
Re-Solving The Live Camera Pointing From Its Targets (A Diagnostic)
=====================================================================
Solves the live CACO05 2025-11-13 extrinsics (EO) of a camera again from
the surveyed targets that calibration was made from, and says, in plain
numbers, how well that pointing is known, where it is not, and what the
lens model and (optionally) the sea horizon say about it. It writes a
re-solved EO next to its report, marked NOT ADOPTED; it never writes into
calibration/, never overwrites a file and never writes a lens model (IO).

WHY. georectify_uncertainty.py measures the error of the live solve and
writes the sigma tables the DEM uses; it refits without weights and
chooses its free parameters by "reproduces the calibration". The design
review of Nov 2026 (docs/VALIDATION.md 4.6c) found that a re-solve can
reproduce the calibration but not improve it (the same 27 targets), that
the parameter set that predicts a left-out target best is the three
angles for BOTH cameras (c1 0.453 m against 0.589 m with all six; c2
0.410 m against 1.060 m), that holding c2's position rests on a position
whose source and accuracy are not documented, and that the leftover
misfit looks like the lens model. This script runs those checks with ONE
solver, so every check describes the EO it writes:
  * the solve: whitened least squares on the picked-minus-projected pixel
    residuals (georectify_uncertainty.project: CIRN P = K R [I | -C], then
    the lens), divided by a pixel sigma (the plain solve's own residual RMS
    per coordinate, 2n - p degrees of freedom), plus a prior on the camera
    position when one is surveyed (--camera-position); Levenberg-Marquardt
    with Marquardt scaling and central differences, numpy only;
  * the free set: --free auto solves all six parameters and the three
    angles (position held at the start EO's), leaves each target out with
    this solver, and keeps the one that predicts the left-out targets
    better (the lower leave-one-out ground RMS); the other is a
    sensitivity, with its ground difference at the data's ranges;
  * a robust (soft-L1) re-run with each target's weight, as a check (the
    written EO is the plain solve); twelve starts +/-1 m, +/-2 deg to show
    the minimum is the only one;
  * leave-one-out per target (pixels and metres; a target over 3 x the
    median is flagged), the Jacobian and jackknife covariances and their
    correlations, the cross-shore sigma at the data's ranges (c1 60-90 m,
    c2 100-300 m, z = 0 and 1 m) through georectify_uncertainty.monte_carlo,
    its DEM equivalent at slope 0.10, the model-choice term (all six
    against the angles) and the shift from the start EO;
  * the lens test: d1 (and d1 + d2) freed with the EO, judged by
    leave-one-out per target (not an F-test on 2n pixel coordinates, which
    treats U and V of one target as independent), each target dropped in
    turn to see whether one carries it, and where the targets sit in the
    lens against the frame corners and the horizon;
  * the sea horizon (--horizon, frames of the calibration's days only
    unless --horizon-any-epoch): a DIAGNOSTIC, never a constraint. On the
    frame's one-sided horizon segment tilt, roll and the lens are not
    separable, so a constraint would put a lens error into the pointing.

WHAT IT CANNOT DO. It cannot make the EO much more accurate than the
live solve: the targets are the same. It cannot tell a lens error from a
pointing error from a camera moved since 2025-11-13, and it cannot pin
the tilt where the waterlines run beyond the targets (c2 out to 300 m,
the targets to 150 m). What would: a GNSS survey of the mounts with the
antenna-to-lens offset (--camera-position), targets at 200-400 m in c2's
view, an in-situ lens calibration through the housing window, and frames
of the calibration day for the horizon. It uses no survey of any kind
(no lidar, RTK, .tif, .asc or DEM: an EO fitted to a survey makes that
survey's comparisons circular, docs/VALIDATION.md section 3.0), and never
opens the 23 Oct 2024 ground control held for the blind check of
calibration B (gcpfree_eo.py), calibration_history.csv or anything named
for that day. It never touches <cam>_geom_sigma.csv. Adopting a
calibration is a separate change, recorded as a dated entry in
docs/VALIDATION.md (section 3.0's fix clause applies).

Usage:
    python3 calibrate_eo.py --out-dir /tmp/eo_refit                 # both cameras, the live inputs
    python3 calibrate_eo.py --out-dir D --cameras c2 --free all       # a sensitivity
    python3 calibrate_eo.py --out-dir D --cameras c1 --camera-position E N Z 0.03 0.05
    python3 calibrate_eo.py --out-dir D --cameras c1 --horizon <frames of 13-16 Nov 2025>
    python3 calibrate_eo.py --out-dir D --cameras c1 --io P --eo P --targets P --picks P
-> <out-dir>/CACO05_<cam>_20251113_EO-refit.yaml, <cam>_eo_report.txt, <cam>_eo_report.json.
Exit 0 done, 2 an input or output refused (the message names it), 6 a usage error.
"""

from __future__ import annotations

import ast
import argparse
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import georectify_uncertainty as GU                                              # noqa: E402
from georectify import IO_KEYS, load_extrinsics, load_intrinsics, pixel_to_ground  # noqa: E402

CAL = HERE / "calibration"
IO_FILE = "CACO05_{cam}_" + GU.IO_DATE + "_IO.yaml"
EO_FILE = "CACO05_{cam}_" + GU.EO_DATE + "_EO-CV.yaml"
EO_NAMES = GU.EO_NAMES
FREE_SETS = {"all": (0, 1, 2, 3, 4, 5), "angles": (3, 4, 5)}
FREE_TEXT = {"all": "6 parameters (x y z azimuth tilt roll)",
             "angles": "3 angles (azimuth tilt roll), the position held at the start EO's",
             "all+prior": "6 parameters with the surveyed camera position as a prior"}
# The lens terms the IO test frees (indices into the CIRN IO vector, georectify.IO_KEYS)
D1, D2 = IO_KEYS.index("d1"), IO_KEYS.index("d2")
LENS_STEP = 1e-6
# The data's ranges (m from the camera; the station's waterlines: c1 63-89 m, c2 105-307 m p10-p90,
# docs/VALIDATION.md 4.10) and planes where the ground uncertainty is given; the DEM equivalent's slope
DATA_RANGES_M = {"c1": (60, 75, 90), "c2": (100, 150, 200, 250, 300)}
Z_PLANES = (0.0, 1.0)
DEM_SLOPE = 0.10
# Multi-start: perturbations of the start EO (free parameters only) and what counts as the same minimum
STARTS, START_POS_M, START_ANG_DEG = 12, 1.0, 2.0
SAME_MIN_M, SAME_MIN_DEG = 0.01, 0.001
# Leave-one-out: a target whose ground error is over this many times the median is flagged; --free auto counts
# sets within TIE_M of the lowest leave-one-out RMS as tied, and takes the one with the fewer parameters
LOO_FLAG = 3.0
TIE_M = 0.001
# Robust re-run: soft-L1 on each target's pixel residual norm, scale ROBUST_K x its median (the MAD of
# a deviation from a perfect fit, taken about zero)
ROBUST_K = 1.4826
# Lens test: flagged when the best model with the lens term predicts the left-out targets better than the best
# model without it (over the sets tested) by more than this share of the leave-one-out ground RMS. On c1's real
# geometry with 3 px of noise and the lens unchanged, 2 of 50 seeded trials pass it (c2: 8 of 50); with d1 off
# by -0.1, 29 of 30 (c2: 16 of 30) [SYNTH, scratch study of Nov 2026; test_calibrate_eo.py pins c1's]
LENS_GAIN_MIN = 0.10
# Sea horizon: the station's pointing-step grid (survey_products.py), and how far a frame may be from the EO's day
HORIZON_MAX_DAYS = 3
HORIZON_PRECISION = ("precision about 0.03 deg of tilt and 0.05 deg of roll (c1's clear days of 29 Sep - 2 Oct "
                     "2026, docs/VALIDATION.md 4.6); its absolute accuracy is untested")
MC_SAMPLES, SEED = 2000, 0
# What the repository says of where the live EO came from (the design review's open questions), said when the
# start EO is the live file; the clause on the targets only when the targets and picks are the live 2025-11-13
# ones (live_note)
LIVE_NOTES = {
    "c1": ("The live EO's provenance is not settled: calibration/README.md says c1 carries the 2025-11-04 values, "
           "georectify_uncertainty.py that c1 was solved from the 2025-11-13 targets.",
           " If it came from a 4 Nov solve, these {n} targets of 2025-11-13 are an out-of-sample check of it, not "
           "its own fit."),
    "c2": ("The live EO's provenance is not settled: c2's solve held its position, and where that position came "
           "from, and how accurate it is, is not documented.",
           " Solving all six parameters here moves the camera {all_move:.2f} m from that position.")}
LIVE_ASK = " Ask the calibration provider."
EXIT_OK, EXIT_REFUSED, EXIT_USAGE = 0, 2, 6
NOT_ADOPTED = ("NOT ADOPTED: a diagnostic re-solve; adopting a calibration is a separate change recorded as a "
               "dated entry in docs/VALIDATION.md (section 3.0's fix clause applies).")

# Inputs never opened: GU.refuse_held (named for 23 Oct 2024) and, here, the calibration history, the other
# common spellings of that day, survey files (lidar, RTK, rasters, DEMs) and survey-fitted pointings. The
# spellings: year first (2024-10-23, 2024_10_23, 20241023, 2024-Oct-23), month first (10-23-2024, Oct_23_2024,
# October 23 2024, the station's frame names ...Wed.Oct.23_16_00_00.GMT.2024...) and day first (23-10-2024,
# 23Oct2024, 23_October_2024); '_' counts as a separator. A frame's 10-digit epoch is checked separately
# (refuse_frame_day).
_SEP = r"[\W_]*"
HELD_EXTRA = re.compile(
    r"2024[-_. ]?10[-_. ]?23"
    r"|(?<!\d)10[-_. ]?23[-_. ]?2024(?!\d)"
    r"|(?<!\d)23[-_. ]?10[-_. ]?2024(?!\d)"
    r"|oct(ober)?" + _SEP + r"23(?!\d).*2024"
    r"|(?<!\d)23" + _SEP + r"oct(ober)?" + _SEP + r"2024"
    r"|2024" + _SEP + r"oct(ober)?" + _SEP + r"23(?!\d)", re.I)
HISTORY_NAME = "calibration_history.csv"
SURVEY_SUFFIXES = (".tif", ".tiff", ".asc", ".las", ".laz", ".img", ".nc", ".vrt")
SURVEY_NAME = re.compile(r"lidar|rtk|checkshot|(^|[^a-z])(dsm|dtm|dem)([^a-z]|$)|_corr_eo", re.I)


class Refused(GU.Refused):
    """An input or output this tool does not open or write (the reason is printed; exit 2)."""


# ---------------------------------------------------------------------
# Guards
# ---------------------------------------------------------------------

MAX_LINKS = 40


def _names(path):
    """The names a path is known by, without opening it: (every component of it, of each link it goes through and
    of where it ends; the file names alone: its own, each link's target and the real one's). The path is resolved
    as the kernel opens it, one component at a time: a link is followed where it stands and a '..' after it
    climbs from where the link led (os.path.abspath would drop 'link/..' as text and name another file)."""
    raw = str(path)
    parts, files = [], [os.path.basename(raw.rstrip(os.sep)) or raw]
    hops = [0]

    def walk(base, rest):
        for name in rest.split(os.sep):
            if name in ("", "."):
                continue
            parts.append(name)
            if name == "..":
                base = os.path.dirname(base) or os.sep
                continue
            nxt = os.path.join(base, name)
            if os.path.islink(nxt):
                hops[0] += 1
                if hops[0] > MAX_LINKS:
                    raise Refused(f"{raw}: more than {MAX_LINKS} links to follow; not opened")
                target = os.readlink(nxt)
                files.append(os.path.basename(target.rstrip(os.sep)) or target)
                base = walk(os.sep if os.path.isabs(target) else base, target)
            else:
                base = nxt
        return base

    start = os.sep if os.path.isabs(raw) else os.path.realpath(os.getcwd())
    parts += [n for n in Path(start).parts if n != os.sep]
    end = walk(start, raw)
    real = Path(os.path.realpath(raw))
    parts += list(Path(end).parts) + list(real.parts)
    files += [os.path.basename(end), real.name]
    return [n for n in parts if n not in ("", os.sep, ".", "..")], files


def refuse_input(path, what="input"):
    """The path, or Refused (before anything is opened): named for 23 Oct 2024 (itself, a folder it is in, or
    a link it goes through), the calibration history, a survey file or a survey-fitted pointing."""
    p = Path(path)
    GU.refuse_held(p)
    parts, files = _names(p)
    for n in parts:
        if GU.HELD_NAME.search(n) or HELD_EXTRA.search(n):
            raise Refused(f"{path}: named for 23 Oct 2024 ({n}) -- that day's ground control and survey are held for "
                          f"the blind check of calibration B (gcpfree_eo.py); not opened")
        if n.lower() == HISTORY_NAME:
            raise Refused(f"{path}: {HISTORY_NAME} is never opened by this tool (it holds calibration A); not opened")
    for n in files:
        if n.lower().endswith(SURVEY_SUFFIXES) or SURVEY_NAME.search(n):
            raise Refused(f"{path} ({what}): a survey or a survey-fitted pointing ({n}); this tool takes no survey "
                          f"input of any kind (an EO fitted to a survey makes that survey's comparisons circular, "
                          f"docs/VALIDATION.md 3.0); not opened")
    return p


def archived_names(path=None):
    """gcpfree_eo.ARCHIVED_SHA256's file names, read with ast (never imported)."""
    src = Path(path or HERE / "gcpfree_eo.py").read_text(encoding="utf-8")
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "ARCHIVED_SHA256" for t in node.targets):
            return set(ast.literal_eval(node.value))
    return set()


def inside(path, folder):
    a, b = os.path.realpath(str(path)), os.path.realpath(str(folder))
    return a == b or a.startswith(b.rstrip(os.sep) + os.sep)


def refuse_out_dir(out_dir):
    if inside(out_dir, CAL):
        raise Refused(f"--out-dir {out_dir}: is calibration/ (or inside it); this tool never writes there")
    if os.path.exists(str(out_dir)) and not os.path.isdir(str(out_dir)):
        raise Refused(f"--out-dir {out_dir}: exists and is not a folder")
    return Path(out_dir)


def refuse_output(path, archived=None):
    """Refused unless path is a new file outside calibration/ whose name is not a calibration's."""
    p = Path(path)
    archived = archived_names() if archived is None else archived
    if inside(p.parent, CAL):
        raise Refused(f"{p}: inside calibration/; this tool never writes there")
    if p.name.endswith("_EO-CV.yaml") or p.name.endswith("_EO.yaml"):
        raise Refused(f"{p}: named as a calibration (_EO-CV.yaml / _EO.yaml); a re-solve is never written as one")
    if p.name in archived:
        raise Refused(f"{p}: the name of a calibration pinned for the blind check (gcpfree_eo.ARCHIVED_SHA256)")
    if os.path.lexists(str(CAL / p.name)):
        raise Refused(f"{p}: calibration/ has a file of this name; a re-solve never takes a calibration's name")
    if os.path.lexists(str(p)):
        raise Refused(f"{p}: exists; this tool never overwrites a file")
    return p


def eo_day(eo_path):
    """The day in a calibration's file name (CACO05_c1_20251113_EO-CV.yaml -> 2025-11-13), or None."""
    from pointing_daily import eo_day as _d
    return _d(eo_path)


def eo_station(eo_path):
    from pointing_daily import eo_station as _s
    return _s(eo_path)


def output_names(cam, eo_path):
    day = eo_day(eo_path)
    stem = f"{eo_station(eo_path) or 'CACO05'}_{cam}_{day.replace('-', '') if day else 'undated'}_EO-refit.yaml"
    return {"eo": stem, "txt": f"{cam}_eo_report.txt", "json": f"{cam}_eo_report.json"}


# ---------------------------------------------------------------------
# The solver
# ---------------------------------------------------------------------

def solve(io, eo0, xyz, uv, free, sigma_px=1.0, weights=None, prior=None, io_free=(), iterations=200):
    """
    The EO (and the lens terms io_free) that best project the targets onto
    their picks: whitened least squares, residuals (projected minus picked,
    all U then all V) x sqrt(weight) / sigma_px, plus (x - E)/sh, (y - N)/sh,
    (z - Z)/sz for the free position parameters when prior = (E, N, Z, sh, sz).
    Levenberg-Marquardt with Marquardt's diagonal scaling and central
    differences, as georectify_uncertainty.fit_eo (whose result this is
    with no weights and no prior). free: indices into x y z azimuth tilt
    roll (radians). -> dict eo, io, p, r (whitened), resid_px (2n), J
    (whitened, at the solution), cost (r.r), converged, iterations.
    """
    free = np.asarray(free, int)
    io_free = np.asarray(io_free, int)
    nf = len(free)
    step = np.r_[np.array(GU.FIT_STEP)[free], np.full(len(io_free), LENS_STEP)]
    w = np.ones(len(xyz)) if weights is None else np.asarray(weights, float)
    sw = np.sqrt(w) / float(sigma_px)
    prow = [(j, int(i)) for j, i in enumerate(free) if i < 3] if prior is not None else []

    def unpack(p):
        e, o = eo0.copy(), io.copy()
        e[free] = p[:nf]
        o[io_free] = p[nf:]
        return e, o

    def raw(p):
        e, o = unpack(p)
        U, V = GU.project(xyz, o, e)
        return np.r_[U - uv[:, 0], V - uv[:, 1]]

    def res(p):
        r = raw(p) * np.r_[sw, sw]
        if prow:
            r = np.r_[r, [(p[j] - prior[i]) / (prior[3] if i < 2 else prior[4]) for j, i in prow]]
        return r

    def jac(p):
        cols = []
        for k in range(len(p)):
            d = np.zeros(len(p))
            d[k] = step[k]
            cols.append((res(p + d) - res(p - d)) / (2.0 * step[k]))
        return np.column_stack(cols)

    p = np.r_[eo0[free], io[io_free]].astype(float)
    r = res(p)
    cost = float(r @ r)
    lam, converged, it = 1e-3, False, 0
    for it in range(1, iterations + 1):
        J = jac(p)
        A, g = J.T @ J, J.T @ r
        moved, rel, dp = False, 0.0, np.zeros(len(p))
        while lam < 1e12:
            try:
                dp = -np.linalg.solve(A + lam * np.diag(np.diag(A)), g)
            except np.linalg.LinAlgError:
                lam *= 10.0
                continue
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
    e, o = unpack(p)
    return {"eo": e, "io": o, "p": p, "r": r, "resid_px": raw(p), "J": jac(p), "cost": cost,
            "converged": converged, "iterations": it, "free": free, "io_free": io_free}


def pixel_sigma(fit, n, p):
    """The plain solve's residual RMS per coordinate (2n - p degrees of freedom), px."""
    rp = fit["resid_px"]
    return float(np.sqrt(float(rp @ rp) / max(2 * n - p, 1)))


def target_px(resid_px):
    """Per-target pixel residual norm from a 2n residual (all U then all V)."""
    n = len(resid_px) // 2
    return np.hypot(resid_px[:n], resid_px[n:])


def eo_move(a, b):
    """(position change m, largest angle change deg) between two EOs."""
    return float(np.linalg.norm(a[:3] - b[:3])), float(np.max(np.abs(np.rad2deg(a[3:] - b[3:]))))


def cov_jacobian(fit):
    """(J'J)^-1 of the whitened Jacobian (its scale is the pixel sigma's) and the reduced chi-square."""
    J, r = fit["J"], fit["r"]
    A = J.T @ J
    try:
        C = np.linalg.inv(A)
    except np.linalg.LinAlgError:
        C = np.linalg.pinv(A)
    dof = max(len(r) - J.shape[1], 1)
    return C, float(r @ r) / dof


def cov_jackknife(sols):
    n = len(sols)
    d = sols - sols.mean(axis=0)
    return (n - 1.0) / n * (d.T @ d)


def t_quantile(prob, dof):
    """Student's t quantile (two-sided interval half-width for probability prob, e.g. 0.95 -> 2.074 on 22 dof),
    numpy only: Simpson's rule on the density, then bisection. The Jacobian's intervals use it, as the pixel
    sigma is estimated from the same residuals (2n - p degrees of freedom)."""
    import math
    nu = float(dof)
    c = math.exp(math.lgamma((nu + 1) / 2) - math.lgamma(nu / 2)) / math.sqrt(nu * math.pi)

    def cdf_half(x):            # P(0 < T < x)
        t = np.linspace(0.0, x, 4001)
        f = c * (1.0 + t * t / nu) ** (-(nu + 1) / 2)
        h = x / 4000.0
        return h / 3.0 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum())
    lo, hi = 0.0, 60.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if 2.0 * cdf_half(mid) < prob:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def sds(C, free):
    """{name: sd} in m and deg."""
    s = np.sqrt(np.maximum(np.diag(C), 0.0))
    return {EO_NAMES[i]: float(s[j]) if i < 3 else float(np.rad2deg(s[j])) for j, i in enumerate(free)}


def correlations(C, free):
    s = np.sqrt(np.maximum(np.diag(C), 1e-300))
    R = C / np.outer(s, s)
    names = [EO_NAMES[i] for i in free]
    return {f"{names[a]}-{names[b]}": float(R[a, b]) for a in range(len(free)) for b in range(a + 1, len(free))}


# ---------------------------------------------------------------------
# One free set: the plain solve, its starts, leave-one-out and covariances
# ---------------------------------------------------------------------

def plain_solve(io, eo0, xyz, uv, free, prior=None, io_free=()):
    """The plain (unweighted) solve, its pixel sigma, and, with a prior, the solve whitened by that sigma.
    -> (fit, sigma_px)."""
    f0 = solve(io, eo0, xyz, uv, free, io_free=io_free)
    s = pixel_sigma(f0, len(xyz), len(free) + len(io_free))
    if prior is None:
        f = solve(io, f0["eo"], xyz, uv, free, sigma_px=s, io_free=io_free)
    else:
        f = solve(io, f0["eo"], xyz, uv, free, sigma_px=s, prior=prior, io_free=io_free)
    return f, s


def multistart(io, eo_start, xyz, uv, free, sigma_px, prior, base, starts, seed):
    """base (the solve from the start EO) and starts - 1 more from the start EO perturbed by up to
    +/-START_POS_M (free position parameters) and +/-START_ANG_DEG (free angles), seeded.
    -> dict: best fit, n, costs, spread (m, deg) about the best, one_minimum, lower_elsewhere."""
    rs = np.random.RandomState(seed)
    free = np.asarray(free, int)
    fits, starts_used = [base], [eo_start.copy()]
    for _ in range(max(int(starts), 1) - 1):
        e0 = eo_start.copy()
        for i in free:
            e0[i] += rs.uniform(-START_POS_M, START_POS_M) if i < 3 else np.deg2rad(rs.uniform(-START_ANG_DEG,
                                                                                           START_ANG_DEG))
        starts_used.append(e0)
        fits.append(solve(io, e0, xyz, uv, free, sigma_px=sigma_px, prior=prior))
    costs = np.array([f["cost"] for f in fits])
    k = int(np.argmin(costs))
    best = fits[k]
    dpos = [float(np.max(np.abs(f["eo"][:3] - best["eo"][:3]))) for f in fits]
    dang = [float(np.max(np.abs(np.rad2deg(f["eo"][3:] - best["eo"][3:])))) for f in fits]
    same = [a <= SAME_MIN_M and b <= SAME_MIN_DEG for a, b in zip(dpos, dang)]
    lower = (k != 0 and not same[0] and costs[0] - costs[k] > 1e-9 * max(costs[0], 1e-30))
    return {"best": best if lower else base, "n": len(fits), "costs": costs.tolist(),
            "start_moves": [list(eo_move(s, eo_start)) for s in starts_used],
            "spread_m": max(dpos), "spread_deg": max(dang), "one_minimum": all(same), "n_same": int(sum(same)),
            "lower_elsewhere": bool(lower), "lower_by": float(costs[0] - costs[k]) if lower else 0.0}


def leave_one_out(io, eo_from, xyz, uv, free, sigma_px, prior=None, io_free=()):
    """Each target left out of the solve in turn, with the same solver, sigma and prior, started at the
    full solve; the left-out target projected from its pick to its surveyed height.
    -> dict dE, dN, ground (m), px, sols (n x p, the free EO parameters then the lens terms)."""
    n = len(xyz)
    dE, dN, px, sols = np.zeros(n), np.zeros(n), np.zeros(n), []
    for k in range(n):
        keep = np.arange(n) != k
        f = solve(io, eo_from, xyz[keep], uv[keep], free, sigma_px=sigma_px, prior=prior, io_free=io_free)
        de, dn = GU.ground_errors(f["io"], f["eo"], xyz[k:k + 1], uv[k:k + 1])
        U, V = GU.project(xyz[k:k + 1], f["io"], f["eo"])
        dE[k], dN[k] = de[0], dn[0]
        px[k] = float(np.hypot(U[0] - uv[k, 0], V[0] - uv[k, 1]))
        sols.append(f["p"])
    return {"dE": dE, "dN": dN, "ground": np.hypot(dE, dN), "px": px, "sols": np.array(sols)}


def rms(a):
    a = np.asarray(a, float)
    return float(np.sqrt(np.mean(a * a))) if a.size else float("nan")


def candidate(name, io, eo_start, xyz, uv, prior, starts, seed):
    """Everything about one free set: -> dict."""
    free = FREE_SETS["all"] if name == "all+prior" else FREE_SETS[name]
    pr = prior if name == "all+prior" else None
    base, s = plain_solve(io, eo_start, xyz, uv, free, pr)
    ms = multistart(io, eo_start, xyz, uv, free, s, pr, base, starts, seed)
    fit = ms["best"]
    loo = leave_one_out(io, fit["eo"], xyz, uv, free, s, pr)
    cj, chi2 = cov_jacobian(fit)
    ck = cov_jackknife(loo["sols"])
    dof = max(len(fit["r"]) - len(free), 1)
    dE, dN = GU.ground_errors(io, fit["eo"], xyz, uv)
    med = float(np.median(loo["ground"]))
    return {"name": name, "free": list(free), "prior": pr, "fit": fit, "eo": fit["eo"], "sigma_px": s,
            "multistart": ms, "loo": loo, "loo_rms": rms(loo["ground"]), "loo_px_rms": rms(loo["px"]),
            "loo_flagged": [int(k) for k in np.flatnonzero(loo["ground"] > LOO_FLAG * med)],
            "cov": {"jacobian": cj, "jackknife": ck}, "chi2_red": chi2, "dof": dof,
            "t68": t_quantile(0.6827, dof), "t95": t_quantile(0.95, dof),
            "in_ground": np.hypot(dE, dN), "in_px": target_px(fit["resid_px"]),
            "px_rms": rms(fit["resid_px"]), "move": eo_move(fit["eo"], eo_start)}


def choose(cands, mode):
    """The free set used: --free all / angles, or (auto) the lower leave-one-out ground RMS (a tie within
    1 mm: the fewer parameters). -> (name, why)."""
    if mode != "auto":
        name = "all+prior" if mode == "all" and "all+prior" in cands else mode
        return name, f"chosen with --free {mode}"
    best = min(cands[k]["loo_rms"] for k in cands)
    tied = [k for k in cands if cands[k]["loo_rms"] <= best + TIE_M]
    name = min(tied, key=lambda k: (len(cands[k]["free"]), cands[k]["loo_rms"], k))
    txt = ", ".join(f"{k} {cands[k]['loo_rms']:.3f} m" for k in sorted(cands, key=lambda k: cands[k]["loo_rms"]))
    return name, f"the lower leave-one-out ground RMS of the sets solved (--free auto): {txt}"


def robust_solve(io, eo_from, xyz, uv, free, sigma_px, prior=None, iterations=50):
    """Soft-L1 IRLS on each target's pixel residual norm a (weight 1/sqrt(1 + (a/s)^2), s = ROBUST_K x
    median a, re-estimated each pass). -> dict fit, weights, scale_px, passes."""
    w = np.ones(len(xyz))
    f = solve(io, eo_from, xyz, uv, free, sigma_px=sigma_px, prior=prior)
    s, passes = float("nan"), 0
    for passes in range(1, iterations + 1):
        a = target_px(f["resid_px"])
        s = ROBUST_K * float(np.median(a))
        w_new = 1.0 / np.sqrt(1.0 + (a / max(s, 1e-9)) ** 2)
        f = solve(io, f["eo"], xyz, uv, free, sigma_px=sigma_px, weights=w_new, prior=prior)
        done = np.max(np.abs(w_new - w)) < 1e-4
        w = w_new
        if done:
            break
    return {"fit": f, "eo": f["eo"], "weights": w, "scale_px": s, "passes": passes}


# ---------------------------------------------------------------------
# On the ground, at the data's ranges
# ---------------------------------------------------------------------

def middle_columns(io):
    """The sigma tables' columns (every 256 px) in the middle half of the frame."""
    c = GU.grid_columns(io)
    nu = float(io[0])
    return c[(c >= 0.25 * nu) & (c <= 0.75 * nu)]


def nodes(io, eo, cam, ranges=None):
    ranges = DATA_RANGES_M.get(cam, (60, 100, 150, 200, 300)) if ranges is None else ranges
    u, v, rr, zz = GU.table_nodes(io, eo, middle_columns(io), ranges, Z_PLANES)
    return {"u": u, "v": v, "range": rr, "z": zz, "ranges": tuple(ranges)}


def cross_shore(io, eo_a, eo_b, nd, io_b=None):
    """Signed cross-shore shift (m, seaward +) of the node pixels on their planes, eo_b against eo_a."""
    normal, _ = GU.shore_normal()
    Ea, Na = pixel_to_ground(nd["u"], nd["v"], nd["z"], io, eo_a, max_range_m=0)
    Eb, Nb = pixel_to_ground(nd["u"], nd["v"], nd["z"], io if io_b is None else io_b, eo_b, max_range_m=0)
    return normal[0] * (Eb - Ea) + normal[1] * (Nb - Na)


def by_range(nd, values, absolute=True):
    """{(range, z): median over the columns} of a per-node value."""
    out = {}
    vals = np.abs(values) if absolute else np.asarray(values, float)
    for z in Z_PLANES:
        for R in nd["ranges"]:
            m = (nd["range"] == R) & (nd["z"] == z) & np.isfinite(vals)
            out[(R, z)] = float(np.median(vals[m])) if m.any() else float("nan")
    return out


def ground_sigma(io, eo, free, covs, nd, samples, seed):
    """Cross-shore sigma at the nodes for each covariance (GU.monte_carlo). -> {name: per-node sigma}, the
    larger one's name (the larger median over the nodes)."""
    normal, _ = GU.shore_normal()
    out = {}
    for k, C in covs.items():
        sxx, syy, sxy, _ = GU.monte_carlo(io, eo, free, C, nd["u"], nd["v"], nd["z"], samples, seed)
        with np.errstate(invalid="ignore"):
            out[k] = np.sqrt(normal[0] ** 2 * sxx + normal[1] ** 2 * syy + 2 * normal[0] * normal[1] * sxy)
    med = {k: float(np.nanmedian(v)) if np.isfinite(v).any() else -1.0 for k, v in out.items()}
    return out, max(med, key=med.get)


# ---------------------------------------------------------------------
# The lens model (IO) test
# ---------------------------------------------------------------------

def lens_geometry(io, eo, uv):
    """Where the targets sit in the lens: normalised radius of the picks, the frame corners and the horizon on
    the station's column grid. -> dict."""
    def rad(u, v):
        return np.hypot((np.asarray(u, float) - io[2]) / io[4], (np.asarray(v, float) - io[3]) / io[5])
    nu, nv = float(io[0]), float(io[1])
    from estimate_eo_rotation import horizon_rows
    cols = horizon_columns(io)
    hr = horizon_rows(io, eo, cols)
    fin = np.isfinite(hr)
    inn = fin.copy()
    inn[fin] = (hr[fin] >= 0) & (hr[fin] < nv)
    hrad = rad(cols[inn], hr[inn]) if inn.any() else np.array([np.nan])
    return {"picks_r_max": float(rad(uv[:, 0], uv[:, 1]).max()), "picks_r_median": float(np.median(rad(uv[:, 0],
                                                                                                      uv[:, 1]))),
            "corners_r": float(rad([0, nu - 1, 0, nu - 1], [0, 0, nv - 1, nv - 1]).max()),
            "horizon_r": [float(np.nanmin(hrad)), float(np.nanmax(hrad))] if inn.any() else None,
            "horizon_rows": [float(hr[inn].min()), float(hr[inn].max())] if inn.any() else None,
            "horizon_cols_in_frame": [int(inn.sum()), int(len(cols))], "top_pick_row": float(uv[:, 1].min())}


def lens_test(io, eo_start, xyz, uv, cand, io_free, drop_one=True, cache=None):
    """One lens refit: io_free (d1, or d1 + d2) freed with the candidate's free set. The misfit change, the
    fitted terms, the camera move, and leave-one-out ground RMS with and without the term; flagged when the
    term lowers it by more than LENS_GAIN_MIN. drop_one: each target dropped in turn and the test redone on the
    rest (which target, if any, carries a flag). cache: a dict that keeps the refits without the lens term
    (the same for every term set of one candidate). -> dict."""
    free, pr = cand["free"], cand["prior"]
    f, s_l = plain_solve(io, cand["eo"], xyz, uv, free, pr, io_free)
    loo = leave_one_out(f["io"], f["eo"], xyz, uv, free, s_l, pr, io_free)
    without, with_ = cand["loo_rms"], rms(loo["ground"])
    gain = 1.0 - with_ / without if without > 0 else float("nan")
    out = {"set": cand["name"], "terms": [IO_KEYS[i] for i in io_free], "fit": f,
           "io_values": [float(f["io"][i]) for i in io_free], "io_start": [float(io[i]) for i in io_free],
           "px_rms_without": cand["px_rms"], "px_rms_with": rms(f["resid_px"]),
           "move_vs_plain": eo_move(f["eo"], cand["eo"]), "move_vs_start": eo_move(f["eo"], eo_start),
           "loo_rms_without": without, "loo_rms_with": with_, "gain": gain, "flagged": bool(gain > LENS_GAIN_MIN),
           "loo_with": loo["ground"], "improved_targets": int(np.sum(loo["ground"] < cand["loo"]["ground"])),
           "drop_one": None}
    if drop_one:
        n = len(xyz)
        rows = []
        cache = {} if cache is None else cache
        for j in range(n):
            keep = np.arange(n) != j
            key = (cand["name"], j)
            if key not in cache:
                a, sa = plain_solve(io, cand["eo"], xyz[keep], uv[keep], free, pr)
                cache[key] = rms(leave_one_out(io, a["eo"], xyz[keep], uv[keep], free, sa, pr)["ground"])
            la = cache[key]
            b, sb = plain_solve(f["io"], f["eo"], xyz[keep], uv[keep], free, pr, io_free)
            lb = leave_one_out(b["io"], b["eo"], xyz[keep], uv[keep], free, sb, pr, io_free)
            g = 1.0 - rms(lb["ground"]) / la
            rows.append((j, la, rms(lb["ground"]), g, bool(g > LENS_GAIN_MIN)))
        out["drop_one"] = rows
        # the targets without which the verdict would be the other one
        out["carried_by"] = [j for j, _, _, _, fl in rows if fl != out["flagged"]]
    return out


def lens_tests(io, eo_start, xyz, uv, cands, chosen, drop_one=True, term_sets=((D1,), (D1, D2))):
    """The lens tests of the chosen set, the three angles and all six (d1, and d1 + d2), and per term set the
    verdict: the best leave-one-out ground RMS with the lens term against the best without it, over those sets
    (lens_verdict). drop_one: every test, so every verdict says whether one target carries it. -> (tests,
    verdicts)."""
    sets = list(dict.fromkeys([chosen, "angles", "all"]))
    cache = {}
    tests = [lens_test(io, eo_start, xyz, uv, cands[k], terms, drop_one=drop_one, cache=cache)
             for k in sets for terms in term_sets]
    return tests, [lens_verdict(cands, tests, [IO_KEYS[i] for i in terms]) for terms in term_sets]


def lens_verdict(cands, tests, terms):
    """Does a lens term predict the left-out targets better than any EO alone? The best leave-one-out RMS with
    `terms` free (over the sets tested) against the best without them (the same sets): flagged when it is lower
    by more than LENS_GAIN_MIN. With every set's drop-one run, the targets without which the verdict flips."""
    ts = [t for t in tests if t["terms"] == list(terms)]
    without = min(cands[t["set"]]["loo_rms"] for t in ts)
    best = min(ts, key=lambda t: t["loo_rms_with"])
    gain = 1.0 - best["loo_rms_with"] / without
    v = {"terms": list(terms), "sets": [t["set"] for t in ts], "without": without,
         "without_set": min((t["set"] for t in ts), key=lambda k: cands[k]["loo_rms"]), "with": best["loo_rms_with"],
         "with_set": best["set"], "gain": gain, "flagged": bool(gain > LENS_GAIN_MIN), "drop_one": None,
         "carried_by": None}
    if all(t["drop_one"] is not None for t in ts):
        n = len(ts[0]["drop_one"])
        rows = []
        for j in range(n):
            a = min(t["drop_one"][j][1] for t in ts)
            b = min(t["drop_one"][j][2] for t in ts)
            rows.append((j, a, b, 1.0 - b / a, bool(1.0 - b / a > LENS_GAIN_MIN)))
        v["drop_one"] = rows
        v["carried_by"] = [j for j, _, _, _, fl in rows if fl != v["flagged"]]
    return v


# ---------------------------------------------------------------------
# The sea horizon (diagnostic only)
# ---------------------------------------------------------------------

def horizon_columns(io):
    """survey_products' pointing step's column grid."""
    return np.arange(150, int(io[0]) - 150, 100)


FRAME_EPOCH = re.compile(r"^(\d{10})\.")


def frame_meta(path):
    """(epoch or None, camera or None, station or None) from a frame's name
    ('1763136000.Fri.Nov.14_16_00_00.GMT.2025.CACO05.c1.timex.jpg')."""
    name = Path(path).name
    m = FRAME_EPOCH.match(name)
    parts = name.split(".")
    cams = [p for p in parts if re.fullmatch(r"c\d", p)]
    stations = [p for p in parts if re.fullmatch(r"[A-Z]+\d+", p)]
    return (int(m.group(1)) if m else None, cams[0] if len(cams) == 1 else None,
            stations[0] if len(stations) == 1 else None)


def refuse_frame_day(path):
    """Refused (before it is opened) when a frame's epoch falls on 23 Oct 2024 (UTC): in its own name, the name of
    each link it goes through, or the name of the file it reaches."""
    for name in _names(path)[1]:
        epoch = frame_meta(name)[0]
        if epoch is not None and datetime.fromtimestamp(epoch, tz=timezone.utc).strftime("%Y-%m-%d") == "2024-10-23":
            raise Refused(f"{path}: a frame of 23 Oct 2024 ({name}), the day held for the blind check of "
                          f"calibration B; not opened")
    return Path(path)


def real_frame_meta(path):
    """frame_meta of the file a frame path reaches (links followed): the file that is read."""
    return frame_meta(os.path.realpath(str(path)))


def screen_frame(path, cam, eo_path, max_days, any_epoch):
    """Why a frame is not used, before it is opened ('' when it may be): the epoch, camera and station in the
    name of the file it reaches (a link's own name is not trusted) against the start EO's. -> (reason, days
    after the EO's day or None, utc text)."""
    refuse_frame_day(path)
    epoch, fcam, fst = real_frame_meta(path)
    if epoch is None:
        return "no 10-digit epoch at the start of its name (the file it reaches)", None, ""
    t = datetime.fromtimestamp(epoch, tz=timezone.utc)
    utc = t.strftime("%Y-%m-%d %H:%M UTC")
    if fcam is None:
        return "no camera (c1 / c2) in its name", None, utc
    if fcam != cam:
        return f"a frame of {fcam}, not {cam}", None, utc
    st = eo_station(eo_path)
    if fst is not None and st is not None and fst != st:
        return f"a frame of station {fst}, not the EO's {st}", None, utc
    day = eo_day(eo_path)
    if day is None:
        return ("" if any_epoch else "the start EO's name carries no date to check the frame's epoch against "
                                     "(--horizon-any-epoch overrides)"), None, utc
    days = (t.date() - datetime.strptime(day, "%Y-%m-%d").date()).days
    if (days < 0 or days > max_days) and not any_epoch:
        when = "before" if days < 0 else "after"
        return (f"{abs(days)} day(s) {when} the EO's day {day} (within 0-{max_days} days only; "
                f"--horizon-any-epoch overrides)"), days, utc
    return "", days, utc


def horizon_frame(gray, io, eos):
    """One frame: the observed horizon on the station grid (estimate_eo_rotation.observed_horizon, searched
    about the first EO's prediction), survey_products' clear-horizon tests, and the tilt / roll change that
    explains it against each EO of eos {label: (io, eo)} (horizon_check.fit_tilt_roll). -> dict."""
    import survey_products as SP
    from estimate_eo_rotation import horizon_rows, observed_horizon
    from horizon_check import fit_tilt_roll, FIT_OK_PX
    cols = horizon_columns(io)
    first = next(iter(eos))
    pred = horizon_rows(eos[first][0], eos[first][1], cols)
    fin = np.isfinite(pred)
    inn = fin.copy()
    inn[fin] = (pred[fin] >= 0) & (pred[fin] < float(io[1]))
    rec = {"cols_in_frame": int(inn.sum()), "cols": int(len(cols))}
    obs = observed_horizon(gray, cols, pred, window=SP.HORIZON_WINDOW_PX)
    rec["n_obs"] = int(np.isfinite(obs).sum())
    if rec["n_obs"] < 3:
        rec.update(clear=False, doubt="no horizon edge found (fewer than 3 columns with an edge)")
        return rec
    edge = SP.horizon_edge(gray, cols, pred, obs, window=SP.HORIZON_WINDOW_PX)
    if not edge:
        rec.update(clear=False, doubt="no horizon edge found")
        return rec
    row = {"horizon_prominence": edge["prominence"], "horizon_rival": edge["rival"],
           "horizon_sky_above": edge["sky_above"], "horizon_to_limit_px": edge["to_limit_px"]}
    rec["edge"] = edge
    doubts = SP.horizon_doubts(row)
    fits = {}
    for k, (o, e) in eos.items():
        (dt, dr), res = fit_tilt_roll(o, e, cols, obs)
        d = obs - horizon_rows(o, e, cols)
        fits[k] = {"dtilt_deg": float(dt), "droll_deg": float(dr), "fit_px": float(np.median(res)),
                   "offset_px": float(np.median(np.abs(d[np.isfinite(d)]))) if np.isfinite(d).any() else float("nan")}
    rec["fits"] = fits
    if doubts:
        rec.update(clear=False, doubt="not a clear horizon: " + SP.doubts_text(row))
    elif fits[first]["fit_px"] > FIT_OK_PX:
        rec.update(clear=False, doubt=f"the tilt + roll fit leaves {fits[first]['fit_px']:.1f} px (over "
                                      f"horizon_check.FIT_OK_PX {FIT_OK_PX:g})")
    else:
        rec.update(clear=True, doubt="")
    return rec


def horizon_diagnostic(cam, frames, io, eo_path, eos, max_days=HORIZON_MAX_DAYS, any_epoch=False):
    """The frames of one camera: screened by name, read (cv2), judged and fitted (horizon_frame). eos: {label:
    (io, eo)}, the first the start EO. -> dict frames (each with its status), medians of the clear ones per
    label, warnings."""
    out, warn = [], []
    for p in frames:
        why, days, utc = screen_frame(p, cam, eo_path, max_days, any_epoch)
        rec = {"frame": str(p), "utc": utc, "days_after_eo": days}
        if Path(os.path.realpath(str(p))).name != Path(p).name:
            rec["reaches"] = os.path.realpath(str(p))
        if why:
            rec.update(status="refused", doubt=why)
            out.append(rec)
            continue
        import cv2
        g = cv2.imread(str(p), cv2.IMREAD_GRAYSCALE)
        if g is None:
            rec.update(status="unreadable", doubt="cv2 cannot read it")
            out.append(rec)
            continue
        if g.shape[1] != int(io[0]) or g.shape[0] != int(io[1]):
            rec.update(status="refused", doubt=f"{g.shape[1]} x {g.shape[0]} px, not the IO's {int(io[0])} x "
                                               f"{int(io[1])}")
            out.append(rec)
            continue
        rec.update(horizon_frame(g, io, eos))
        rec["status"] = "used" if rec["clear"] else "not clear"
        out.append(rec)
    used = [r for r in out if r["status"] == "used"]
    med = {}
    for k in eos:
        if used:
            med[k] = {"dtilt_deg": float(np.median([r["fits"][k]["dtilt_deg"] for r in used])),
                      "droll_deg": float(np.median([r["fits"][k]["droll_deg"] for r in used])),
                      "fit_px": float(np.median([r["fits"][k]["fit_px"] for r in used])), "frames": len(used)}
    far = [r for r in used if r.get("days_after_eo") is None or not 0 <= r["days_after_eo"] <= max_days]
    if far:
        span = [r["days_after_eo"] for r in far if r.get("days_after_eo") is not None]
        warn.append(f"{len(far)} frame(s) used with --horizon-any-epoch"
                    + (f", {min(span)} to {max(span)} days from the EO's day" if span else "")
                    + ": the camera may have moved or turned since the calibration (c1 panned -1.21 deg on 26-28 "
                      "Sep 2026, docs/VALIDATION.md 4.6; 11 months separate the 2025-11-13 targets and the 7 Oct "
                      "2026 horizon), so an offset here may be a later pointing, not the calibration's")
    return {"frames": out, "medians": med, "warnings": warn, "precision": HORIZON_PRECISION,
            "cols_in_frame": used[0]["cols_in_frame"] if used else None}


# ---------------------------------------------------------------------
# One camera, end to end
# ---------------------------------------------------------------------

def shown(path):
    """A path as written in the provenance: the file it reaches (links followed), relative to the repository when
    inside it."""
    a = os.path.realpath(str(path))
    return os.path.relpath(a, str(HERE)) if inside(a, HERE) else a


def eo_deg(eo):
    return {k: float(eo[i]) if i < 3 else float(np.rad2deg(eo[i])) for i, k in enumerate(EO_NAMES)}


def analyse(cam, paths, free_mode="auto", prior=None, starts=STARTS, seed=SEED, samples=MC_SAMPLES,
            io_test=True, horizon=None, horizon_max_days=HORIZON_MAX_DAYS, horizon_any_epoch=False,
            drop_one=True, ranges=None):
    """Everything for one camera. paths: {io, eo, targets, picks} (already screened by refuse_input).
    prior: (E, N, Z, sigma_h, sigma_z) or None. horizon: frame paths or None. -> report dict (numpy inside;
    jsonable() for the file)."""
    io, eo0 = load_intrinsics(str(paths["io"])), load_extrinsics(str(paths["eo"]))
    nums, xyz, uv = GU.load_targets(paths["targets"], paths["picks"])
    if len(nums) < 5:
        raise ValueError(f"{cam}: only {len(nums)} target(s) in both {paths['targets']} and {paths['picks']}")
    n = len(nums)
    rng_t = np.hypot(xyz[:, 0] - eo0[0], xyz[:, 1] - eo0[1])
    names = ["all", "angles"] + (["all+prior"] if prior is not None else [])
    cands = {k: candidate(k, io, eo0, xyz, uv, prior, starts, seed) for k in names}
    chosen, why = choose(cands, free_mode)
    c = cands[chosen]
    eo_s = c["eo"]
    rep = {"camera": cam, "n_targets": n, "nums": nums, "target_range_m": rng_t, "io": io, "eo_start": eo0,
           "prior": prior, "candidates": cands, "chosen": chosen, "why": why, "free_mode": free_mode,
           "eo": eo_s, "starts": starts, "seed": seed, "samples": samples}
    rep["reproduces"] = bool(c["move"][0] <= GU.FIT_POS_TOL_M and c["move"][1] <= GU.FIT_ANG_TOL_DEG)
    rep["robust"] = robust_solve(io, eo_s, xyz, uv, c["free"], c["sigma_px"], c["prior"])
    rep["robust"]["move_vs_plain"] = eo_move(rep["robust"]["eo"], eo_s)
    # on the ground at the data's ranges
    nd = nodes(io, eo_s, cam, ranges)
    sig, larger = ground_sigma(io, eo_s, c["free"], c["cov"], nd, samples, seed)
    g = {"nodes": nd, "sigma": sig, "larger": larger, "columns": middle_columns(io)}
    g["model_choice"] = cross_shore(io, cands["angles"]["eo"], cands["all"]["eo"], nd)
    g["change_from_start"] = cross_shore(io, eo0, eo_s, nd)
    g["robust_shift"] = cross_shore(io, eo_s, rep["robust"]["eo"], nd)
    g["others"] = {k: cross_shore(io, eo_s, cands[k]["eo"], nd) for k in cands if k != chosen}
    deg_px = float(np.rad2deg(1.0 / io[5]))
    e_po = eo_s.copy()
    e_po[4] += np.deg2rad(deg_px)
    g["pick_origin"] = {"deg_per_px": deg_px, "shift": cross_shore(io, eo_s, e_po, nd)}
    g["footprint"] = GU.footprints(io, eo_s, nd["u"], nd["v"], nd["z"], GU.shore_normal()[0])[0]
    rep["ground"] = g
    # the lens
    rep["lens_geometry"] = lens_geometry(io, eo_s, uv)
    rep["lens"], rep["lens_verdicts"] = [], []
    if io_test:
        rep["lens"], rep["lens_verdicts"] = lens_tests(io, eo0, xyz, uv, cands, chosen, drop_one)
        from estimate_eo_rotation import horizon_rows
        from horizon_check import fit_tilt_roll
        cols = horizon_columns(io)
        base_rows = horizon_rows(io, eo_s, cols)
        for t in rep["lens"]:
            rows = horizon_rows(t["fit"]["io"], t["fit"]["eo"], cols)
            d = rows - base_rows
            ok = np.isfinite(d)
            t["horizon_shift_px"] = float(np.median(d[ok])) if ok.any() else float("nan")
            if ok.sum() >= 3:
                (dt, dr), _ = fit_tilt_roll(io, eo_s, cols, rows)
                t["horizon_shift_deg"] = [float(dt), float(dr)]
            else:
                t["horizon_shift_deg"] = None
    # the sea horizon
    rep["horizon"] = None
    if horizon:
        eos = {"start": (io, eo0), "solved": (io, eo_s)}
        for t in rep["lens"]:
            if t["terms"] == ["d1"]:
                eos[f"d1 refit ({t['set']})"] = (t["fit"]["io"], t["fit"]["eo"])
        rep["horizon"] = horizon_diagnostic(cam, horizon, io, paths["eo"], eos, horizon_max_days, horizon_any_epoch)
    return rep


# ---------------------------------------------------------------------
# The report and the EO file
# ---------------------------------------------------------------------

def fmt_eo_delta(a, b):
    d = eo_deg(a)
    s = eo_deg(b)
    return ", ".join(f"{k} {d[k] - s[k]:+.3f} {'m' if i < 3 else 'deg'}" for i, k in enumerate(EO_NAMES))


def sd_text(sd):
    return ", ".join(f"{k} {v:.3f} {'m' if k in ('x', 'y', 'z') else 'deg'}" for k, v in sd.items())


def range_table(nd, rows):
    """Lines of a table: one row per (label, per-node values, absolute?), columns the (range, z) cells."""
    head = "  " + f"{'':34}" + "".join(f"{f'{R:.0f} m z{z:+.0f}':>11}" for z in Z_PLANES for R in nd["ranges"])
    out = [head]
    for label, vals, absolute in rows:
        br = by_range(nd, vals, absolute)
        out.append("  " + f"{label:34}" + "".join(f"{br[(R, z)]:11.2f}" if np.isfinite(br[(R, z)]) else f"{'-':>11}"
                                                  for z in Z_PLANES for R in nd["ranges"]))
    return out


def drop_text(t, nums):
    """What dropping each target in turn does to a lens verdict (t: a lens_test or lens_verdict record)."""
    gains = [x[3] for x in t["drop_one"]]
    span = f"gains {100 * min(gains):+.0f}% to {100 * max(gains):+.0f}%"
    flip = t["carried_by"]
    if not flip:
        return f"each target dropped in turn keeps the verdict ({span})"
    who = ", ".join(f"{nums[j]} ({100 * t['drop_one'][j][3]:+.1f}%)" for j in flip)
    after = "it is no longer flagged" if t["flagged"] else "it would be flagged"
    if len(flip) == 1:
        return f"THE VERDICT RESTS ON ONE TARGET: without {who} {after} ({span})"
    return f"THE VERDICT IS MARGINAL: without any one of {who} {after} ({span})"


def report_text(rep, prov):
    cam, c = rep["camera"], rep["candidates"][rep["chosen"]]
    cands, nums = rep["candidates"], rep["nums"]
    L = [f"calibrate_eo.py: camera {cam} -- a diagnostic re-solve of its EO from its surveyed targets",
         NOT_ADOPTED, ""]
    L += [f"inputs: {k} {v}" for k, v in prov["inputs"].items()]
    L += [f"git head {prov['git_head']}; written {prov['written_utc']}", ""]
    r = rep["target_range_m"]
    L.append(f"{rep['n_targets']} targets, {r.min():.0f}-{r.max():.0f} m from the camera.")
    if rep["prior"] is not None:
        E, N, Z, sh, sz = rep["prior"]
        L.append(f"A surveyed camera position is a prior: E {E:.3f}, N {N:.3f}, z {Z:.3f} m, sigma {sh:g} m "
                 f"horizontal and {sz:g} m vertical (with GNSS, include the antenna-to-lens offset).")
    L += ["", "THE SOLVE"]
    L.append(f"  Whitened least squares on the picked-minus-projected pixels (georectify_uncertainty.project), "
             f"Levenberg-Marquardt with Marquardt scaling and central differences. The pixel sigma is each set's "
             f"plain (unweighted) solve's residual RMS per coordinate (2n - p degrees of freedom).")
    L.append(f"  {'set':10} {'params':>6} {'sigma px':>8} {'RMS px':>7} {'chi2/dof':>8} {'moves m':>8} "
             f"{'moves deg':>9} {'in-sample m':>11} {'LOO m':>7} {'LOO px':>7}")
    for k, d in cands.items():
        L.append(f"  {k:10} {len(d['free']):6d} {d['sigma_px']:8.2f} {d['px_rms']:7.2f} {d['chi2_red']:8.3f} "
                 f"{d['move'][0]:8.3f} {d['move'][1]:9.3f} {rms(d['in_ground']):11.3f} {d['loo_rms']:7.3f} "
                 f"{d['loo_px_rms']:7.1f}")
    L.append(f"  Used: {rep['chosen']}, {FREE_TEXT[rep['chosen']]}; {rep['why']}.")
    L.append("  The angles solve holds the camera at the start EO's position and assumes it exact; where that "
             "position came from, and how accurate it is, is not documented (the calibration provider would know).")
    if rep["prior"] is not None:
        L.append("  With --camera-position the angles set still holds the start EO's position; all+prior uses the "
                 "surveyed one.")
    others = [k for k in cands if k != rep["chosen"]]
    for k in others:
        L.append(f"  Sensitivity, {k}: moves the camera {eo_move(cands[k]['eo'], rep['eo'])[0]:.3f} m and turns it "
                 f"{eo_move(cands[k]['eo'], rep['eo'])[1]:.3f} deg from the solve used (table below).")
    L.append(f"  Change from the start EO: {fmt_eo_delta(rep['eo'], rep['eo_start'])}.")
    L.append("  " + ("It REPRODUCES the start calibration" if rep["reproduces"] else "It does NOT reproduce the "
                     "start calibration") + f" (moves {c['move'][0]:.3f} m, {c['move'][1]:.3f} deg; the rule: "
             f"<= {GU.FIT_POS_TOL_M:g} m and <= {GU.FIT_ANG_TOL_DEG:g} deg, georectify_uncertainty).")
    for k, d in cands.items():
        ms = d["multistart"]
        if ms["n"] <= 1:
            txt = "1 start only (the start EO): whether the minimum is the only one is not tested"
        elif ms["one_minimum"]:
            txt = (f"all {ms['n']} starts reached one minimum (spread {ms['spread_m']:.4f} m, "
                   f"{ms['spread_deg']:.5f} deg)")
        else:
            txt = (f"{ms['n_same']} of {ms['n']} starts reached the best minimum (within {SAME_MIN_M:g} m / "
                   f"{SAME_MIN_DEG:g} deg); the others ended up to {ms['spread_m']:.3f} m / {ms['spread_deg']:.3f} "
                   f"deg away")
        if ms["lower_elsewhere"]:
            txt += ("; WARNING: A START REACHED A LOWER COST ELSEWHERE THAN THE SOLVE FROM THE START EO (by "
                    f"{ms['lower_by']:.3g}): that lower minimum is used")
        L.append(f"  Multi-start ({k}; +/-{START_POS_M:g} m, +/-{START_ANG_DEG:g} deg, seed {rep['seed']}): {txt}.")
    L += ["", f"PER TARGET ({rep['chosen']}; leave-one-out with the same solver, sigma and prior)"]
    L.append(f"  {'target':>6} {'range m':>8} {'in m':>7} {'in px':>6} {'LOO m':>7} {'LOO dE':>7} {'LOO dN':>7} "
             f"{'LOO px':>7} {'robust w':>8}")
    lo = c["loo"]
    for k in range(rep["n_targets"]):
        flag = "  <-- over 3 x the median" if k in c["loo_flagged"] else ""
        L.append(f"  {nums[k]:>6} {r[k]:8.1f} {c['in_ground'][k]:7.3f} {c['in_px'][k]:6.1f} {lo['ground'][k]:7.3f} "
                 f"{lo['dE'][k]:+7.3f} {lo['dN'][k]:+7.3f} {lo['px'][k]:7.1f} {rep['robust']['weights'][k]:8.2f}"
                 f"{flag}")
    L.append(f"  In-sample RMS {rms(c['in_ground']):.3f} m; leave-one-out RMS {c['loo_rms']:.3f} m, median "
             f"{np.median(lo['ground']):.3f} m, max {lo['ground'].max():.3f} m (target {nums[int(np.argmax(lo['ground']))]}"
             f"). Flagged (over {LOO_FLAG:g} x the median): "
             + (", ".join(str(nums[k]) for k in c["loo_flagged"]) or "none") + ".")
    rb = rep["robust"]
    low = np.argsort(rb["weights"])[:3]
    L.append(f"  Robust check (soft-L1, scale {ROBUST_K} x the median residual norm = {rb['scale_px']:.1f} px): lowest "
             f"weights " + ", ".join(f"{nums[k]} {rb['weights'][k]:.2f}" for k in low)
             + f"; it moves {rb['move_vs_plain'][0]:.3f} m and {rb['move_vs_plain'][1]:.3f} deg from the plain solve, "
               f"which is the one written (a check, not a gain).")
    L += ["", f"HOW WELL THE POINTING IS KNOWN ({rep['chosen']})"]
    for name in ("jacobian", "jackknife"):
        C = c["cov"][name]
        L.append(f"  {name:9}: {sd_text(sds(C, c['free']))}")
    L.append(f"  The Jacobian is whitened: its scale is the pixel sigma's ({c['sigma_px']:.2f} px); reduced "
             f"chi-square {c['chi2_red']:.3f}" + (" (1 by construction without a prior)" if c["prior"] is None else "")
             + f". Its 68% / 95% intervals are +/-{c['t68']:.3f} / +/-{c['t95']:.3f} sd (Student's t on {c['dof']} "
               f"degrees of freedom: the sigma is estimated from the same residuals).")
    for name in ("jacobian", "jackknife"):
        cor = correlations(c["cov"][name], c["free"])
        big = sorted(cor.items(), key=lambda kv: -abs(kv[1]))[:4]
        L.append(f"  correlations ({name}): " + ", ".join(f"{k} {v:+.2f}" for k, v in big)
                 + (f"; z-tilt {cor['z-tilt']:+.2f}" if "z-tilt" in cor else ""))
    for k in others:
        cor = correlations(cands[k]["cov"]["jacobian"], cands[k]["free"])
        L.append(f"  {k}: tilt sd {sds(cands[k]['cov']['jacobian'], cands[k]['free'])['tilt']:.3f} deg (jacobian), "
                 f"{sds(cands[k]['cov']['jackknife'], cands[k]['free'])['tilt']:.3f} deg (jackknife)"
                 + (f"; z-tilt correlation {cor['z-tilt']:+.2f} (position and angles trade off)" if "z-tilt" in cor
                    else ""))
    po = rep["ground"]["pick_origin"]
    L.append(f"  Not included: the pick origin convention. 1 px of it is {po['deg_per_px']:.4f} deg of tilt (1/fy, "
             f"fy {rep['io'][5]:.1f} px), an excluded error (docs/VALIDATION.md 4.10); nor a camera turned since.")
    g = rep["ground"]
    nd = g["nodes"]
    L += ["", f"ON THE GROUND (cross-shore, m; median over the middle columns "
              f"{' '.join(f'{x:.0f}' for x in g['columns'])}; shore normal {GU.shore_normal()[1]:g} deg)"]
    rows = [(f"sigma ({k}{', used' if k == g['larger'] else ''})", g["sigma"][k], True) for k in ("jacobian",
                                                                                                  "jackknife")]
    rows.append((f"DEM equivalent at slope {DEM_SLOPE:g}", DEM_SLOPE * g["sigma"][g["larger"]], True))
    rows.append(("model choice: all minus angles", g["model_choice"], True))
    with np.errstate(invalid="ignore"):
        rows.append(("sigma (+) model choice", np.sqrt(g["sigma"][g["larger"]] ** 2 + g["model_choice"] ** 2), True))
    rows.append(("shift from the start EO", g["change_from_start"], False))
    rows.append(("robust minus plain", g["robust_shift"], False))
    for k, v in g["others"].items():
        rows.append((f"{k} minus {rep['chosen']}", v, False))
    rows.append(("1 px of pick origin (tilt)", po["shift"], True))
    rows.append(("pixel footprint (cross-shore)", g["footprint"], True))
    L += range_table(nd, rows)
    L.append("  The sigma is the solve's own (the larger of the two covariances, drawn "
             f"{rep['samples']} times, seed {rep['seed']}); beyond the targets' "
             f"{r.max():.0f} m it is an extrapolation of the fit, not a measurement. The model-choice row is how far "
             "the two parameter sets disagree there.")
    lg = rep["lens_geometry"]
    L += ["", "THE LENS MODEL (IO; never refitted into a file)"]
    L.append(f"  The picks reach normalised lens radius {lg['picks_r_max']:.2f} (median {lg['picks_r_median']:.2f}); "
             f"the frame's corners are at {lg['corners_r']:.2f}"
             + (f"; the horizon at {lg['horizon_r'][0]:.2f}-{lg['horizon_r'][1]:.2f} (rows "
                f"{lg['horizon_rows'][0]:.0f}-{lg['horizon_rows'][1]:.0f}, inside the frame on "
                f"{lg['horizon_cols_in_frame'][0]} of {lg['horizon_cols_in_frame'][1]} grid columns)"
                if lg["horizon_r"] else "; the horizon is not inside the frame")
             + f"; the highest pick is at row {lg['top_pick_row']:.0f}. The targets can test the lens model but not "
               f"refit it near the top of the frame.")
    for v in rep.get("lens_verdicts", []):
        terms = "+".join(v["terms"])
        L.append(f"  VERDICT, {terms}: the best leave-one-out ground RMS with it is {v['with']:.3f} m ({v['with_set']} + "
                 f"{terms}), without it {v['without']:.3f} m ({v['without_set']}): {100 * v['gain']:+.0f}%, "
                 + ("FLAGGED: a lens term predicts the left-out targets better than the EO alone"
                    if v["flagged"] else "not flagged") + f" (rule: better by over {100 * LENS_GAIN_MIN:.0f}%).")
        if v["drop_one"] is not None:
            L.append("      " + drop_text(v, nums))
    if rep["lens"]:
        L.append("  Per set (leave-one-out with and without the term, the same set):")
    for t in rep["lens"]:
        terms = "+".join(t["terms"])
        vals = ", ".join(f"{k} {v:+.3f} (IO {s:+.3f})" for k, v, s in zip(t["terms"], t["io_values"], t["io_start"]))
        L.append(f"  {t['set']} + {terms}: pixel RMS {t['px_rms_without']:.2f} -> {t['px_rms_with']:.2f}; {vals}; the "
                 f"camera moves {t['move_vs_plain'][0]:.3f} m, {t['move_vs_plain'][1]:.3f} deg; leave-one-out "
                 f"{t['loo_rms_without']:.3f} -> {t['loo_rms_with']:.3f} m ({100 * t['gain']:+.0f}% gain, better on "
                 f"{t['improved_targets']} of {rep['n_targets']} targets): "
                 + ("better by over " if t["flagged"] else "not better by ") + f"{100 * LENS_GAIN_MIN:.0f}%.")
        if t.get("horizon_shift_deg") is not None:
            L.append(f"      it moves the predicted horizon {t['horizon_shift_px']:+.1f} px (median), as "
                     f"{t['horizon_shift_deg'][0]:+.2f} deg of tilt and {t['horizon_shift_deg'][1]:+.2f} of roll")
        if t["drop_one"] is not None:
            L.append("      " + drop_text(t, nums))
    if rep["lens"]:
        L.append("  A lens term that fits better is evidence about the lens model, not a calibration: the fitted terms "
                 "are not a lens anyone measured, and a new IO needs an in-situ lens calibration.")
    h = rep["horizon"]
    if h is not None:
        L += ["", "THE SEA HORIZON (diagnostic only, never a constraint: on the frame's one-sided horizon segment tilt, "
                  "roll and the lens are not separable)"]
        for fr in h["frames"]:
            line = (f"  {Path(fr['frame']).name}" + (f" (a link to {Path(fr['reaches']).name})" if fr.get("reaches")
                                                     else "") + f"  {fr['utc']}  {fr['status']}")
            if fr.get("fits"):
                line += "; " + "; ".join(f"vs {k}: tilt {v['dtilt_deg']:+.3f}, roll {v['droll_deg']:+.3f} deg, fit "
                                         f"{v['fit_px']:.1f} px" for k, v in fr["fits"].items())
            if fr.get("doubt"):
                line += f"  ({fr['doubt']})"
            L.append(line)
        for k, m in h["medians"].items():
            L.append(f"  median of {m['frames']} clear frame(s) against the {k} EO: tilt {m['dtilt_deg']:+.3f} deg, roll "
                     f"{m['droll_deg']:+.3f} deg (fit {m['fit_px']:.1f} px)")
        if not h["medians"]:
            L.append("  no clear frame: nothing to report")
        if h["cols_in_frame"] is not None:
            L.append(f"  the horizon is inside the frame on {h['cols_in_frame']} of {len(horizon_columns(rep['io']))} "
                     f"grid columns")
        L.append(f"  The horizon reads with {h['precision']}.")
        if any(k.startswith("d1 refit") for k in h["medians"]):
            L.append("  Against a d1 refit (the targets only, the lens term freed) the offset left is its 'd1 refit' "
                     "line: evidence about lens versus pointing, not a decision.")
        for w in h["warnings"]:
            L.append(f"  WARNING: {w}")
    L += ["", "WHAT THIS IS NOT"]
    L.append(f"  It never touches calibration/{cam}_geom_sigma.csv, and nothing reads its output. Adopting a "
             "calibration is a separate change, recorded as a dated entry in docs/VALIDATION.md (section 3.0's fix "
             "clause applies).")
    L.append("  It cannot tell a lens error from a pointing error from a camera moved since the calibration, or pin "
             "the tilt beyond the targets. A GNSS survey of the mount (with the antenna-to-lens offset), far targets "
             "(200-400 m), an in-situ lens calibration and frames of the calibration day would.")
    note = live_note(rep, prov)
    if note:
        L.append("  " + note)
    return L


def live_note(rep, prov):
    """LIVE_NOTES for the camera when the start EO is the live file ('' otherwise): c1's clause on the targets
    only when the targets and picks are the live 2025-11-13 files (their count from the run); c2's move from
    the run's all-six solve."""
    cam = rep["camera"]
    if cam not in LIVE_NOTES or Path(str(prov["eo_path"])).name != EO_FILE.format(cam=cam):
        return ""
    head, clause = LIVE_NOTES[cam]
    live = prov.get("live_inputs", {})
    if cam == "c1":
        tail = clause.format(n=rep["n_targets"]) if live.get("targets") and live.get("picks") else ""
    else:
        tail = clause.format(all_move=float(rep["candidates"]["all"]["move"][0]))
    return head + tail + LIVE_ASK


def jsonable(x):
    if isinstance(x, dict):
        return {str(k) if not isinstance(k, tuple) else "|".join(str(v) for v in k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return jsonable(x.tolist())
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    return x


def report_json(rep, prov):
    """Every number of the report, plain JSON (angles in deg where named so, NaN as null)."""
    def cand_json(d):
        ms = dict(d["multistart"])
        ms.pop("best", None)
        return {"free": [EO_NAMES[i] for i in d["free"]], "prior": d["prior"], "sigma_px": d["sigma_px"],
                "pixel_rms": d["px_rms"], "chi2_red": d["chi2_red"], "eo": eo_deg(d["eo"]),
                "move_from_start_m_deg": d["move"], "converged": d["fit"]["converged"],
                "in_sample_ground_m": d["in_ground"], "in_sample_px": d["in_px"],
                "loo_ground_m": d["loo"]["ground"], "loo_dE_m": d["loo"]["dE"], "loo_dN_m": d["loo"]["dN"],
                "loo_px": d["loo"]["px"], "loo_rms_m": d["loo_rms"], "loo_px_rms": d["loo_px_rms"],
                "loo_flagged": [rep["nums"][k] for k in d["loo_flagged"]], "dof": d["dof"],
                "t_68": d["t68"], "t_95": d["t95"],
                "sd_jacobian": sds(d["cov"]["jacobian"], d["free"]), "sd_jackknife": sds(d["cov"]["jackknife"],
                                                                                        d["free"]),
                "corr_jacobian": correlations(d["cov"]["jacobian"], d["free"]),
                "corr_jackknife": correlations(d["cov"]["jackknife"], d["free"]),
                "cov_jacobian_m_rad": d["cov"]["jacobian"], "cov_jackknife_m_rad": d["cov"]["jackknife"],
                "multistart": ms}
    g, nd = rep["ground"], rep["ground"]["nodes"]
    ground = {"columns": g["columns"], "ranges_m": nd["ranges"], "z_planes": Z_PLANES, "covariance_used": g["larger"],
              "shore_normal_deg": GU.shore_normal()[1], "dem_slope": DEM_SLOPE,
              "nodes": {"u": nd["u"], "v": nd["v"], "range_m": nd["range"], "z": nd["z"]},
              "sigma_cross_shore_m": {k: by_range(nd, v) for k, v in g["sigma"].items()},
              "dem_equivalent_m": by_range(nd, DEM_SLOPE * g["sigma"][g["larger"]]),
              "model_choice_all_minus_angles_m": by_range(nd, g["model_choice"]),
              "shift_from_start_m": by_range(nd, g["change_from_start"], False),
              "robust_minus_plain_m": by_range(nd, g["robust_shift"], False),
              "others_minus_used_m": {k: by_range(nd, v, False) for k, v in g["others"].items()},
              "pick_origin": {"deg_tilt_per_px": g["pick_origin"]["deg_per_px"],
                              "cross_shore_m_per_px": by_range(nd, g["pick_origin"]["shift"])},
              "pixel_footprint_cross_shore_m": by_range(nd, g["footprint"]),
              "per_node": {"sigma": g["sigma"], "model_choice": g["model_choice"],
                           "shift_from_start": g["change_from_start"]}}
    lens = []
    for t in rep["lens"]:
        x = {k: v for k, v in t.items() if k not in ("fit",)}
        x["eo"] = eo_deg(t["fit"]["eo"])
        if t["drop_one"] is not None:
            x["drop_one"] = [{"dropped": rep["nums"][j], "loo_without_m": a, "loo_with_m": b, "gain": gg,
                              "flagged": fl} for j, a, b, gg, fl in t["drop_one"]]
            x["carried_by"] = [rep["nums"][j] for j in t["carried_by"]]
        lens.append(x)
    out = {"camera": rep["camera"], "not_adopted": NOT_ADOPTED, "provenance": prov, "targets": rep["nums"],
           "target_range_m": rep["target_range_m"], "eo_start": eo_deg(rep["eo_start"]), "io": dict(zip(IO_KEYS,
                                                                                                         rep["io"])),
           "prior": rep["prior"], "free_mode": rep["free_mode"], "chosen": rep["chosen"], "why": rep["why"],
           "eo_refit": eo_deg(rep["eo"]), "reproduces_start": rep["reproduces"],
           "reproduce_rule": {"pos_m": GU.FIT_POS_TOL_M, "ang_deg": GU.FIT_ANG_TOL_DEG},
           "candidates": {k: cand_json(d) for k, d in rep["candidates"].items()},
           "robust": {"weights": rep["robust"]["weights"], "scale_px": rep["robust"]["scale_px"],
                      "passes": rep["robust"]["passes"], "eo": eo_deg(rep["robust"]["eo"]),
                      "move_vs_plain_m_deg": rep["robust"]["move_vs_plain"]},
           "ground": ground, "lens_geometry": rep["lens_geometry"], "lens_tests": lens, "lens_gain_min": LENS_GAIN_MIN,
           "lens_verdicts": [dict(v, drop_one=None if v["drop_one"] is None else
                                  [{"dropped": rep["nums"][j], "loo_without_m": a, "loo_with_m": b, "gain": gg,
                                    "flagged": fl} for j, a, b, gg, fl in v["drop_one"]],
                                  carried_by=None if v["carried_by"] is None else
                                  [rep["nums"][j] for j in v["carried_by"]]) for v in rep["lens_verdicts"]],
           "horizon": rep["horizon"], "starts": rep["starts"], "seed": rep["seed"], "samples": rep["samples"],
           "geom_sigma_touched": False}
    return jsonable(out)


def live_paths(cam):
    """The live 2025-11-13 inputs of a camera (the defaults)."""
    return {"io": CAL / IO_FILE.format(cam=cam), "eo": CAL / EO_FILE.format(cam=cam),
            "targets": CAL / GU.TARGETS.format(cam=cam), "picks": CAL / GU.PICKS.format(cam=cam)}


def provenance(cam, paths, rep):
    inputs = {k: f"{shown(paths[k])} sha256 {GU.sha256(paths[k])}" for k in ("io", "eo", "targets", "picks")}
    live = {k: os.path.realpath(str(paths[k])) == os.path.realpath(str(v)) for k, v in live_paths(cam).items()}
    return {"inputs": inputs, "git_head": GU.git_head(), "script_sha256": GU.sha256(Path(__file__).resolve()),
            "written_utc": f"{datetime.now(timezone.utc):%Y-%m-%dT%H:%M:%SZ}", "eo_path": shown(paths["eo"]),
            "live_inputs": live}


def eo_file_text(rep, prov, eo_path):
    """The refit EO: '#' provenance, then the CIRN keys as in the live file, then its comment lines."""
    c = rep["candidates"][rep["chosen"]]
    sj, sk = sds(c["cov"]["jacobian"], c["free"]), sds(c["cov"]["jackknife"], c["free"])
    head = [f"calibrate_eo.py: camera {rep['camera']}, a diagnostic re-solve of {Path(str(eo_path)).name} from its "
            f"surveyed targets", NOT_ADOPTED]
    head += [f"{k}: {v}" for k, v in prov["inputs"].items()]
    head += [f"git_head: {prov['git_head']}", f"script_sha256: {prov['script_sha256']}",
             f"written_utc: {prov['written_utc']}",
             f"free: {rep['chosen']} ({FREE_TEXT[rep['chosen']]}) -- {rep['why']}",
             "solver: whitened least squares on the pixel residuals (Levenberg-Marquardt, Marquardt scaling, central "
             "differences), the plain solve (a soft-L1 re-run is a check in the report)"
             + ("; camera-position prior E N Z %.3f %.3f %.3f, sigma %g m / %g m" % tuple(rep["prior"])
                if c["prior"] is not None else ""),
             f"pixel_sigma_px: {c['sigma_px']:.3f} (the plain solve's residual RMS per coordinate, 2n - p dof)",
             f"loo_ground_rms_m: {c['loo_rms']:.3f} ("
             + ", ".join(f"{k} {d['loo_rms']:.3f}" for k, d in rep["candidates"].items()) + ")",
             f"sd_jacobian: {sd_text(sj)}", f"sd_jackknife: {sd_text(sk)}",
             f"reproduces_start: {'yes' if rep['reproduces'] else 'NO'} (moves {c['move'][0]:.3f} m, "
             f"{c['move'][1]:.3f} deg)",
             f"report: {rep['camera']}_eo_report.txt / .json"]
    e = eo_deg(rep["eo"])
    body = [f"x: {e['x']:.4f}", f"y: {e['y']:.4f}", f"z: {e['z']:.4f}", f"azimuth: {e['azimuth']:.6f}",
            f"tilt: {e['tilt']:.6f}", f"roll: {e['roll']:.6f}"]
    comments = []
    with open(str(eo_path), "r") as f:
        for line in f.read().splitlines():
            if line.strip().startswith("#"):
                comments.append(line.strip())
    return "".join(f"# {h}\n" for h in head) + "".join(f"{b}\n" for b in body) + "".join(f"{x}\n" for x in comments)


def write_new(path, text):
    """Write a NEW file (never an existing one: exclusive creation)."""
    with open(str(path), "x") as f:
        f.write(text)


# ---------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------

class Parser(argparse.ArgumentParser):
    def error(self, message):
        self.print_usage(sys.stderr)
        sys.stderr.write(f"{self.prog}: usage error: {message}\n")
        sys.exit(EXIT_USAGE)


def build_parser():
    ap = Parser(description=__doc__.split("\n")[1])
    ap.add_argument("--out-dir", required=True, help="where the refit EO and the reports go (new files only; never "
                                                      "calibration/)")
    ap.add_argument("--cameras", nargs="+", default=["c1", "c2"], choices=["c1", "c2"])
    ap.add_argument("--free", default="auto", choices=["auto", "all", "angles"],
                    help="auto (default): the set with the lower leave-one-out ground RMS; all six; the three angles")
    ap.add_argument("--io", default=None, help=f"default calibration/{IO_FILE.format(cam='<cam>')} (one camera)")
    ap.add_argument("--eo", default=None, help=f"the start EO, default calibration/{EO_FILE.format(cam='<cam>')} "
                                               f"(one camera)")
    ap.add_argument("--targets", default=None, help=f"default calibration/{GU.TARGETS.format(cam='<cam>')} (one camera)")
    ap.add_argument("--picks", default=None, help=f"default calibration/{GU.PICKS.format(cam='<cam>')} (one camera)")
    ap.add_argument("--camera-position", nargs=5, type=float, default=None,
                    metavar=("E", "N", "Z", "SIGMA_H", "SIGMA_Z"),
                    help="a surveyed lens position (m, UTM 19N / NAVD88) and its sigmas: a prior (one camera)")
    ap.add_argument("--horizon", nargs="+", default=None, metavar="FRAME",
                    help="timex frames for the sea-horizon diagnostic (never a constraint)")
    ap.add_argument("--horizon-max-days", type=int, default=HORIZON_MAX_DAYS,
                    help=f"frames up to this many days after the EO's day are used (default {HORIZON_MAX_DAYS})")
    ap.add_argument("--horizon-any-epoch", action="store_true",
                    help="use frames of any date (the report then warns the camera may have moved since)")
    ap.add_argument("--no-io-test", action="store_true", help="skip the lens model test")
    ap.add_argument("--starts", type=int, default=STARTS)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--samples", type=int, default=MC_SAMPLES)
    return ap


def main(argv=None):
    ap = build_parser()
    args = ap.parse_args(argv)
    cams = list(dict.fromkeys(args.cameras))
    if len(cams) > 1 and (args.io or args.eo or args.targets or args.picks):
        ap.error("--io, --eo, --targets and --picks are for one camera: add --cameras c1 (or c2)")
    if len(cams) > 1 and args.camera_position:
        ap.error("--camera-position is for one camera: add --cameras c1 (or c2)")
    if args.camera_position and not all(np.isfinite(args.camera_position)):
        ap.error("--camera-position: E, N, Z, SIGMA_H and SIGMA_Z must all be finite numbers")
    if args.camera_position and (args.camera_position[3] <= 0 or args.camera_position[4] <= 0):
        ap.error("--camera-position: SIGMA_H and SIGMA_Z must be positive")
    if args.starts < 1 or args.samples < 10 or args.horizon_max_days < 0:
        ap.error("--starts must be >= 1, --samples >= 10, --horizon-max-days >= 0")
    if not 0 <= args.seed < 2 ** 32:
        ap.error("--seed must be between 0 and 2**32 - 1")
    try:
        # every input screened by name before anything is opened
        jobs = []
        for cam in cams:
            given = {"io": args.io, "eo": args.eo, "targets": args.targets, "picks": args.picks}
            paths = {k: Path(given[k]) if given[k] else v for k, v in live_paths(cam).items()}
            for k, p in paths.items():
                refuse_input(p, f"--{k}")
            jobs.append((cam, paths))
        frames = {}
        for p in args.horizon or []:
            refuse_input(p, "--horizon")
            refuse_frame_day(p)
            _, fcam, _ = real_frame_meta(p)
            frames.setdefault(fcam if fcam in cams else cams[0], []).append(Path(p))
        out_dir = refuse_out_dir(args.out_dir)
        archived = archived_names()
        outs = {}
        for cam, paths in jobs:
            outs[cam] = {k: refuse_output(out_dir / v, archived) for k, v in output_names(cam, paths["eo"]).items()}
        for cam, paths in jobs:
            for k, p in paths.items():
                if not Path(p).is_file():
                    raise Refused(f"--{k} {p}: not found")
        reps = []
        for cam, paths in jobs:
            print(f"{cam}: solving ...", flush=True)
            rep = analyse(cam, paths, args.free, tuple(args.camera_position) if args.camera_position else None,
                          args.starts, args.seed, args.samples, not args.no_io_test, frames.get(cam),
                          args.horizon_max_days, args.horizon_any_epoch)
            prov = provenance(cam, paths, rep)
            reps.append((cam, paths, rep, prov))
        out_dir.mkdir(parents=True, exist_ok=True)
        for cam, paths, rep, prov in reps:
            o = outs[cam]
            text = report_text(rep, prov)
            write_new(o["eo"], eo_file_text(rep, prov, paths["eo"]))
            write_new(o["txt"], "\n".join(text) + "\n")
            write_new(o["json"], json.dumps(report_json(rep, prov), indent=1, allow_nan=False) + "\n")
            print("\n".join(text))
            print(f"wrote {o['eo']}, {o['txt']}, {o['json']}")
    except GU.Refused as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return EXIT_REFUSED
    except (OSError, ValueError, KeyError) as exc:
        print(f"REFUSED: an input could not be read: {exc}", file=sys.stderr)
        return EXIT_REFUSED
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Is The Measured Setup Right? The Survey-Date Check
====================================================
The independent test of the measured setup (measured_setup.py): on a
survey date with c2 runup timestacks, the swash edge of every stack is put
on the SURVEY's surface instead of a station DEM, so the setup measured is
independent of the station's DEM and of its setup coefficient C; it is
compared, collection by collection, with the model C * sqrt(Hs * L0) of
the same collection's frames, and both setups are scored where it matters,
on the waterlines against the survey. NUC-run: the stacks, the products
and the surveys are there. It writes a report, --json, one CSV of the
measured setups per date in --out-dir (and, with --sharpness-images, the
window sharpness it scored); nothing else.

WHY. Live, the measured setup stands on a model-built DEM of earlier days,
so it inherits that DEM's level at the c2 lines and reads beach change
since as setup (measured_setup.py, rule 2). On the survey's own surface
neither applies: what is left is the measurement itself (the detector, the
still water) and how the timex line relates to the mean swash edge. This
is the BEST CASE for the measured setup, an upper bound on how well it can
do live (it is measured on the very surface it is scored against, with no
DEM error); the live trial (docs/VALIDATION.md 4.12) comes after it.

THE RULE
 1. THE DATES. Each row of surveys/surveys.csv with a survey file and a
    product of survey_products.py (<output root>/<date>/provenance.json)
    whose window has c2 stacks in --ras-dir. A BLIND date (survey_products
    .BLIND_DATES: 2024-10-23, whose GCPs are held for an independent check)
    is skipped before anything of it is opened, and so is a row whose
    survey file is NAMED like a survey of a blind date; the run arms
    survey_products.arm_blind_guard, which refuses any open of such a file.
    A survey row without a path is never looked up by its date. The dates
    with stacks known here: 2026-09-29 (the RTK checkshots, window 26 Sep -
    2 Oct, if archive/ras_c2 holds them: the cron archives every c2 stack);
    2025-01-23 and 2025-03-06 (the 2025 lidars) only once their c2 stacks
    are fetched back from S3 (cleanup.sh moved them), with that period's
    pixel list (--pix) and calibration (the product's own c2 EO, else --eo).
    A date none of whose stacks gives a row (a wrong --pix for their width,
    a --lines the .pix lacks, unreadable stacks) is skipped with the reason.
 2. THE SURFACE. A lidar DSM (compare_dem_survey.read_survey) cropped to the
    c2 lines' footprint (their pixels on the planes SURFACE_Z, +/-
    CROP_MARGIN_M) and kept between SURFACE_Z (-2.5 to +4.0 m NAVD88: the
    swash zone, and below any terrain behind the camera that a ray marched
    from above could otherwise meet). RTK transects (compare_rtk.
    load_transects): only within --rtk-corridor (5 m) of a transect's axis,
    between its first and last shot, each cell the transect's profile at its
    distance along the axis (the beach taken as alongshore-uniform over those
    5 m, nowhere else). The 29 Sep 2026 transects start at ~+1.2 m: only
    high-water stacks whose line runs inside a corridor can pass QC there.
    Each row records the surface (measured_setup.FixedProfile).
 3. THE MEASUREMENT: measured_setup.measure_stack, its rules 1, 3, 4 and 5
    unchanged (the QC, the median over lines 1-2, the uncertainty); the
    still water is the product's own: forcing/water_level.csv (the ADCP or
    the Chatham transfer, 'forcing water level') on a historical date, the
    station's GNSS-R ('gps', marconi_water_level.py) on a live-era date (its
    forcing.json says GNSS-R), so the measured setup is above the level the
    product's lines were given. THE WINDOW'S SHARPNESS: owg_live.py's record
    (--owg-csv) where it covers the date; for a period without it (2025)
    --sharpness-images scores that period's own c2 images of each collection
    (bright, else snap, else timex: score_image_quality.score, the measure
    owg_live.py uses, on the full frame) into --out-dir/<date>_<name>_
    sharpness.csv and uses that, with the distribution per image product in
    the report (the blur limit, 10, was set on bright frames: on snap or
    timex frames it is PROVISIONAL until that distribution is looked at).
    --allow-no-sharpness, with no blur check at all, is the last resort and
    the report says so.
 4. AGAINST THE MODEL. Each collection that passes is joined to the
    product's frames whose file-name epoch lies within 300 s of its start
    (the extractor's own rule); their Hs and Tp give the wave scale w =
    sqrt(Hs L0) and the model's setup C w (setup_model.setup_c, --station-c
    0.037). The C that would match the measured setups is the least-squares
    fit through the origin, C_fit = sum(m w) / sum(w^2), with a 90% DAY
    JACKKNIFE t interval (leave one whole day out, the t of n_days - 1
    degrees of freedom: the collections of a day share their beach, level
    and waves, so treating them as independent would make the interval too
    narrow; fewer than MIN_DAYS_CI 3 days: no interval; the number of days is
    printed beside it). [SYNTH: m = C w + a day effect (sd 0.02-0.05 m) +
    noise (0.03 m), 5-25 collections a day: it covers C 0.86-0.91 of the
    time at 3-10 days, where the day-block percentile bootstrap it replaces
    covered 0.61-0.64 at 3 days and 0.77-0.84 at 5-10.] Reported overall,
    by Hs band (HS_BANDS of the frames' Hs), by the tide's phase (flood /
    ebb / slack: the still water 15 min after the burst's middle less 15
    min before, +/- TIDE_SLACK_M) and by the day's offset from the survey
    (beach change since the survey moves the measurement too). Collections
    the model cannot serve (no Tp: the buoy outages of finding H2) are
    counted.
 5. ON THE SURVEY (the DEM-relevant score). The product's frames as
    evaluate_setup_model.py scores them (load_product: the DEM's frames with
    Hs and Tp, each line re-projected onto the plane still water + its
    setup, the survey-floor rule; set_reference: frames more than 1 m off
    at the station's C are detection failures, left out), each frame with a
    measured setup scored with the model's setup, with the model's setup
    corrected to C's still-water reference (rule 6) where that differs, and
    with its collection's measured one, on the frames where all reach the
    survey. Per frame the median of (line - survey); bias, median, RMSE and
    NMAD, per camera, by Hs band and by phase; and the PAIRED differences
    measured - model (the corrected model where there is one) of RMSE and
    of |bias| on the same frames, each with its 90% day-jackknife t
    interval [SYNTH: frames correlated within a day, 3-7 days: it covered
    0.86-0.88 (RMSE) and 0.83-0.85 (|bias|) of the time]. 'auto' (measured
    where there is one, the model elsewhere) is scored on every frame. The
    2026-09-29 RTK was the survey C was FITTED to (survey_products
    .SETUP_FITS): there the model is circular and the comparison is flagged;
    the 2025 lidars are the independent dates.
 6. THE STILL-WATER REFERENCE (historical dates). C was fitted with lines at
    the GNSS-R level, which holds a share s of the setup (its footprint is
    the surf zone); a 2025 product gives its lines the ADCP or Chatham level,
    which holds none, and its Tp is the ADCP's (VALIDATION 4.1,
    historical_forcing.py). So the model's lines there sit off what the live
    model would give by about +0.02 m - s (C w - the mean setup) plus the
    period currency's offset, while the measured setup, measured above that
    same level, has no such shortfall (its line is the swash's elevation on
    the survey: the level cancels). Unfixed, that tilt alone could decide a
    comparison, most in storms. So on a date whose forcing.json is not
    GNSS-R, the model is also scored corrected: C w + delta, delta = k (w -
    x_typ) - GNSSR_BELOW_CHATHAM_M - tp_bias, from that product's forcing
    record (waves.still_water_reference: k = s C, the slope of GNSS-R -
    Chatham transfer on sqrt(Hs L0) over the 2026 record, x_typ its median
    sqrt(Hs L0); waves.tp_currency_setup_bias_m, the median setup from the
    ADCP period less that from an open-ocean one; GNSS-R ~0.02 m below
    Chatham in the mean, OPUS +/-0.061 m). C_fit is given at that reference
    too (the fit of m - delta). Without a measured share in the record
    (historical_forcing.py measures it only with the 2026 record on the
    computer) the date is reported uncorrected and says it cannot count for
    the decision rule's criteria 1-2.

WHAT IT CANNOT DO
  * Stacks of the survey dates of 2025 are not on the NUC until fetched
    from S3; the 2026 RTK covers only the upper beach (~+1.0 to +1.6 m).
  * The surface is the beach on the survey day: days before and after
    differ by the beach's change (by day offset, 4.). And it is the TRUE
    surface: the live measured setup stands on a model-built DEM whose
    level error it inherits (x ~1.5), which this check cannot see. A win
    here is necessary, not sufficient.
  * The reference correction of rule 6 rests on a share s measured over
    2026 and applied to 2025, a datum offset known to +/-0.061 m, and a
    median period-currency offset (per frame it varies); it corrects the
    model's tilt, not its scatter. On a live-era date both setups stand
    above GNSS-R and no correction applies.
  * The detector's own bias on the synthetic scenes (+0.04 to +0.10 m:
    measured_setup.py) is in the measured setups; still water is the
    product's (GNSS-R holds part of the setup).
  * Frames without Tp are not in the score of 5 (load_product keeps frames
    with waves): the outages are counted, not scored.
  * A window sharpness scored on snap or timex frames is not owg_live.py's
    bright-frame number: its limit is provisional (rule 3).

Usage (on the NUC, after survey_products.py has built the dates):
    python3 evaluate_measured_setup.py --dates 2026-09-29
    python3 evaluate_measured_setup.py --dates 2025-01-23 --ras-dir /path/to/jan2025_c2_stacks \\
        --pix /path/to/c2_timestack.pix.2025 --sharpness-images /path/to/jan2025_c2_images \\
        --json /tmp/msetup_2025-01-23.json
Tests: python3 test_measured_setup.py
"""

import argparse
import csv
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import measured_setup as ms                                  # noqa: E402
from setup_model import setup_c, wave_scale                  # noqa: E402

SURFACE_Z = (-2.5, 4.0)          # m NAVD88: the part of a survey the swash can reach
CROP_MARGIN_M = 50.0             # around the lines' footprint
RTK_CORRIDOR_M = 5.0             # either side of a transect's axis
RTK_CELL_M = 0.5
HS_BANDS = ((0.0, 1.0), (1.0, 2.0), (2.0, 99.0))
TIDE_SLACK_M = 0.02              # |level change over 30 min| below this: slack
TIDE_HALF_S = 900.0
MIN_DAYS_CI = 3                  # days a jackknife interval needs
LEVEL_LABEL = "forcing water level"
GNSSR_BELOW_CHATHAM_M = 0.02     # GNSS-R below Chatham in the mean after its datum fix (historical_forcing.py)
IMAGE_PRODUCTS = ("bright", "snap", "timex")   # --sharpness-images: the first of these a collection has
# one-sided 95% points of Student's t (a two-sided 90% interval), by degrees of freedom; beyond 30 the
# next tabulated df below is used (conservative)
T95 = {1: 6.314, 2: 2.920, 3: 2.353, 4: 2.132, 5: 2.015, 6: 1.943, 7: 1.895, 8: 1.860, 9: 1.833, 10: 1.812,
       11: 1.796, 12: 1.782, 13: 1.771, 14: 1.761, 15: 1.753, 16: 1.746, 17: 1.740, 18: 1.734, 19: 1.729,
       20: 1.725, 21: 1.721, 22: 1.717, 23: 1.714, 24: 1.711, 25: 1.708, 26: 1.706, 27: 1.703, 28: 1.701,
       29: 1.699, 30: 1.697, 40: 1.684, 60: 1.671, 120: 1.658}


# ---------------------------------------------------------------------
# The survey's surface
# ---------------------------------------------------------------------

def line_footprint(uv_lines, io, eo, z=SURFACE_Z):
    """(E0, E1, N0, N1) of the lines' pixels on the planes z (georectify.pixel_to_ground), or None."""
    from georectify import pixel_to_ground
    E, N = [], []
    for uv in uv_lines:
        for zz in z:
            e, n = pixel_to_ground(uv[:, 0], uv[:, 1], zz, io, eo)
            ok = np.isfinite(e) & np.isfinite(n)
            E.append(np.asarray(e)[ok])
            N.append(np.asarray(n)[ok])
    E, N = np.concatenate(E) if E else np.zeros(0), np.concatenate(N) if N else np.zeros(0)
    if not len(E):
        return None
    return float(E.min()), float(E.max()), float(N.min()), float(N.max())


def stack_lines(ras_path, cam, pix, lines):
    """The (u, v) pixels of each line of a stack (runup_from_timestack's .pix for its width), or None."""
    import cv2
    import runup_from_timestack as rft
    ras = cv2.imread(str(ras_path), cv2.IMREAD_UNCHANGED)
    if ras is None:
        return None
    p = rft.pix_for_width(cam, ras.shape[1], pix)
    if p is None:
        return None
    seg = rft.split_lines(p)
    return [p[a:b] for i, (a, b) in enumerate(seg, 1) if i in lines]


def dsm_surface(grid, x0, ytop, cell, name, box=None, z=SURFACE_Z, margin=CROP_MARGIN_M):
    """A DSM (rows north to south) as a timestack_wave_height.Dem: cropped to box (E0, E1, N0, N1) + margin,
    and NaN outside z."""
    import timestack_wave_height as thw
    g = np.asarray(grid, float)
    r0, c0 = 0, 0
    if box is not None:
        e0, e1, n0, n1 = box
        c0 = max(int(np.floor((e0 - margin - x0) / cell)), 0)
        c1 = min(int(np.ceil((e1 + margin - x0) / cell)), g.shape[1])
        r0 = max(int(np.floor((ytop - (n1 + margin)) / cell)), 0)
        r1 = min(int(np.ceil((ytop - (n0 - margin)) / cell)), g.shape[0])
        if c1 <= c0 or r1 <= r0:
            return None
        g = g[r0:r1, c0:c1]
    g = g.copy()
    with np.errstate(invalid="ignore"):
        g[~((g >= z[0]) & (g <= z[1]))] = np.nan
    if not np.isfinite(g).any():
        return None
    xll = x0 + c0 * cell
    yll = ytop - (r0 + g.shape[0]) * cell
    return thw.Dem(grid=g, xll=xll, yll=yll, cell=cell, name=name)


def transect_surface(transects, name, corridor=RTK_CORRIDOR_M, cell=RTK_CELL_M, z=SURFACE_Z):
    """RTK transects (compare_rtk.load_transects) as a Dem: each cell within `corridor` of a transect's axis
    and between its first and last shot takes that transect's profile at its distance along the axis (the
    nearest transect where corridors overlap); NaN elsewhere."""
    import timestack_wave_height as thw
    if not transects:
        return None
    pts = np.vstack([np.c_[t["E"], t["N"]] for t in transects])
    e0, n0 = pts.min(0) - corridor - cell
    e1, n1 = pts.max(0) + corridor + cell
    nc, nr = int(np.ceil((e1 - e0) / cell)), int(np.ceil((n1 - n0) / cell))
    xc = e0 + (np.arange(nc) + 0.5) * cell
    yc = n1 - (np.arange(nr) + 0.5) * cell
    X, Y = np.meshgrid(xc, yc)
    g = np.full(X.shape, np.nan)
    best = np.full(X.shape, np.inf)
    for t in transects:
        rx, ry = X - t["c"][0], Y - t["c"][1]
        s = rx * t["d"][0] + ry * t["d"][1]
        lat = np.abs(-rx * t["d"][1] + ry * t["d"][0])
        inside = (lat <= corridor) & (s >= t["s"][0]) & (s <= t["s"][-1]) & (lat < best)
        g[inside] = np.interp(s[inside], t["s"], t["z"])
        best[inside] = lat[inside]
    with np.errstate(invalid="ignore"):
        g[~((g >= z[0]) & (g <= z[1]))] = np.nan
    return thw.Dem(grid=g, xll=e0, yll=n1 - nr * cell, cell=cell, name=name) if np.isfinite(g).any() else None


def survey_surface(kind, path, box=None, corridor=RTK_CORRIDOR_M):
    """The survey file as a Dem: kind 'dsm' (a DSM) or 'points' (RTK transects). -> (Dem or None, label)."""
    if kind == "dsm":
        from compare_dem_survey import read_survey
        grid, x0, ytop, cell = read_survey(path)
        return (dsm_surface(grid, x0, ytop, cell, Path(path).name, box),
                f"survey DSM {Path(path).name}, {SURFACE_Z[0]:g} to {SURFACE_Z[1]:g} m")
    from compare_rtk import load_transects
    from survey_compare import TRANSECT_SPLIT
    tr = load_transects(path, TRANSECT_SPLIT)
    return (transect_surface(tr, Path(path).name, corridor),
            f"survey transects {Path(path).name}, within {corridor:g} m of {len(tr)} transect axes")


# ---------------------------------------------------------------------
# The window's sharpness from a period's own images (rule 3)
# ---------------------------------------------------------------------

def image_sharpness_table(dirs, camera, lo, hi, out_csv, products=IMAGE_PRODUCTS):
    """
    Per collection of [lo, hi), the sharpness of its first image of `products` in `dirs` (its file-name
    epoch; score_image_quality.score on the full frame, owg_live.py's measure), written as an owg-style
    table (filename, epoch, time_utc, brightness, sharpness, saturation, status, product) that
    timestack_wave_height.WindowSharpness reads. -> {product: [sharpness, ...]} of the collections scored.
    """
    from score_image_quality import score
    found = {}
    for d in dirs:
        for k, prod in enumerate(products):
            for p in Path(d).glob(f"*.{camera}.{prod}.jpg"):
                e = ms.epoch_from_name(p.name)
                if e is None or not lo <= e < hi:
                    continue
                if e not in found or found[e][0] > k:
                    found[e] = (k, prod, p)
    rows, by = [], {}
    for e in sorted(found):
        _k, prod, p = found[e]
        b, sh, sa = score(p)
        rows.append({"filename": p.name, "epoch": int(e),
                     "time_utc": datetime.fromtimestamp(int(e), tz=timezone.utc).isoformat(),
                     "brightness": round(float(b), 2) if np.isfinite(b) else "",
                     "sharpness": round(float(sh), 2) if np.isfinite(sh) else "",
                     "saturation": round(float(sa), 4) if np.isfinite(sa) else "",
                     "status": "ok" if np.isfinite(sh) else "unreadable", "product": prod})
        if np.isfinite(sh):
            by.setdefault(prod, []).append(float(sh))
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["filename", "epoch", "time_utc", "brightness", "sharpness", "saturation",
                                          "status", "product"])
        w.writeheader()
        w.writerows(rows)
    return by


# ---------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------

def fit_c(m, w):
    """The C of m = C w, least squares through the origin (NaN without data)."""
    m, w = np.asarray(m, float), np.asarray(w, float)
    ok = np.isfinite(m) & np.isfinite(w)
    if not ok.any() or not (w[ok] ** 2).sum() > 0:
        return np.nan
    return float((m[ok] * w[ok]).sum() / (w[ok] ** 2).sum())


def t95(df):
    """Student's t, one-sided 95% (two-sided 90%), for df degrees of freedom (T95; NaN below 1)."""
    df = int(df)
    if df < 1:
        return np.nan
    if df in T95:
        return T95[df]
    return T95[max(k for k in T95 if k <= df)]


def day_jackknife_ci(stat, days, arrays, min_days=MIN_DAYS_CI):
    """
    The 90% day-jackknife t interval of stat(*arrays): leave each whole day out in turn, the jackknife
    standard error sqrt((G-1)/G sum (theta_-d - mean)^2) and the t of G-1 degrees of freedom (G days).
    -> (theta, lo, hi, G); lo, hi NaN with fewer than min_days days.
    """
    days = np.asarray(days, object)
    arrays = [np.asarray(a, float) for a in arrays]
    ud = sorted(set(days.tolist()))
    G = len(ud)
    th = stat(*arrays)
    if G < min_days or not np.isfinite(th):
        return th, np.nan, np.nan, G
    loo = np.array([stat(*[a[days != d] for a in arrays]) for d in ud], float)
    if not np.isfinite(loo).all():
        return th, np.nan, np.nan, G
    se = float(np.sqrt((G - 1) / G * ((loo - loo.mean()) ** 2).sum()))
    t = t95(G - 1)
    return th, th - t * se, th + t * se, G


def c_fit_ci(m, w, days, min_days=MIN_DAYS_CI):
    """C_fit and its 90% day-jackknife interval -> (lo, hi); (NaN, NaN) with fewer than min_days days."""
    _th, lo, hi, _g = day_jackknife_ci(fit_c, days, [m, w], min_days)
    return lo, hi


def _rmse(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.sqrt(np.mean(x ** 2))) if len(x) else np.nan


def _absbias(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(abs(x.mean())) if len(x) else np.nan


def paired_diffs(r_ref, r_mea, days, min_days=MIN_DAYS_CI):
    """Measured - comparator on the same frames: RMSE and |bias| differences with their 90% day-jackknife
    intervals -> {rmse: [theta, lo, hi], absbias: [...], days: G}."""
    out = {}
    g = 0
    for key, f in (("rmse", _rmse), ("absbias", _absbias)):
        th, lo, hi, g = day_jackknife_ci(lambda a, b, f=f: f(b) - f(a), days, [r_ref, r_mea], min_days)
        out[key] = [th, lo, hi]
    out["days"] = g
    return out


def diff_stats(d):
    """n, bias (mean), median, RMSE of the differences d (finite ones)."""
    d = np.asarray(d, float)
    d = d[np.isfinite(d)]
    if not len(d):
        return {"n": 0, "bias": None, "median": None, "rmse": None}
    return {"n": int(len(d)), "bias": float(d.mean()), "median": float(np.median(d)),
            "rmse": float(np.sqrt(np.mean(d ** 2)))}


def hs_band(hs):
    for a, b in HS_BANDS:
        if a <= hs < b:
            return f"{a:g}-{b:g} m" if b < 99 else f">= {a:g} m"
    return "unknown"


def tide_phase(level, mid_epoch, half=TIDE_HALF_S, slack=TIDE_SLACK_M):
    """'flood', 'ebb' or 'slack' from the still water `half` s after and before the burst's middle."""
    lv, _ = level.at([mid_epoch - half, mid_epoch + half])
    d = float(lv[1] - lv[0])
    if not np.isfinite(d):
        return "unknown"
    return "flood" if d > slack else ("ebb" if d < -slack else "slack")


# ---------------------------------------------------------------------
# The still-water reference (rule 6)
# ---------------------------------------------------------------------

def reference_of(pdir):
    """
    A product's still-water reference against C's (rule 6), from its forcing/forcing.json: dict kind
    ('gnssr': C's own; 'historical'; 'unknown'), usable (a correction can be made: always for gnssr, where it
    is zero), k, x_typ, share, datum_m, tp_bias_m, note.
    """
    zero = {"k": 0.0, "x_typ": 0.0, "share": None, "datum_m": 0.0, "tp_bias_m": 0.0}
    try:
        f = json.loads((Path(pdir) / "forcing" / "forcing.json").read_text())
    except (OSError, ValueError):
        return dict(zero, kind="unknown", usable=False,
                    note="no readable forcing/forcing.json: the product's still water is not known; not corrected")
    method = str((f.get("water_level") or {}).get("method") or "")
    if method.startswith("GNSS-R"):
        return dict(zero, kind="gnssr", usable=True, note="GNSS-R, C's own reference: no correction")
    wv = f.get("waves") or {}
    swr = wv.get("still_water_reference") or {}
    k, x = ms._num(swr.get("k")), ms._num(swr.get("typical_x"))
    if not (swr.get("quantified") and np.isfinite(k) and np.isfinite(x)):
        return dict(zero, kind="historical", usable=False,
                    note=(f"still water '{method or '?'}': the GNSS-R share s was not measured in this product's "
                          "forcing (historical_forcing.py measures it with the 2026 record on the computer): the "
                          "model is NOT corrected, and this date cannot count for criteria 1-2"))
    tp = ms._num(wv.get("tp_currency_setup_bias_m"))
    note = (f"still water '{method}': model corrected by k (w - x_typ) - {GNSSR_BELOW_CHATHAM_M:g} - tp_bias, "
            f"k {k:+.4f} (s {ms._num(swr.get('share')):.2f}), x_typ {x:.2f} m, tp_bias "
            + (f"{tp:+.3f} m" if np.isfinite(tp) else "not measured (0)"))
    return {"kind": "historical", "usable": True, "k": float(k), "x_typ": float(x), "share": swr.get("share"),
            "datum_m": GNSSR_BELOW_CHATHAM_M, "tp_bias_m": float(tp) if np.isfinite(tp) else 0.0, "note": note}


def reference_delta(w, ref):
    """delta(w) = k (w - x_typ) - datum_m - tp_bias_m (rule 6): add to the product's model setup C w to put its
    line where the live model would. NaN where the reference is not usable."""
    w = np.asarray(w, float)
    if not ref or not ref.get("usable"):
        return np.full(w.shape, np.nan)
    return ref["k"] * (w - ref["x_typ"]) - ref["datum_m"] - ref["tp_bias_m"]


def compare_model(colls, station_c, ref=None):
    """Rule 4 over the collections joined to their frames: [{m, w, hs, tp, phase, day_offset, ...}] ->
    dict of the fit, its interval and the differences, overall and by subset; with a usable historical
    reference (rule 6) also against the corrected model and C_fit at C's reference."""
    corr = bool(ref and ref.get("usable") and ref.get("kind") == "historical")

    def block(sub):
        m = np.array([c["m"] for c in sub], float)
        w = np.array([c["w"] for c in sub], float)
        days = np.array([c["day"] for c in sub], object)
        cf = fit_c(m, w)
        lo, hi = c_fit_ci(m, w, days)
        model = station_c * w
        ratio = m / w
        out = dict(diff_stats(m - model), c_fit=cf, c_ci90=[lo, hi], days=len(set(days.tolist())),
                   c_ratio_median=float(np.median(ratio)) if len(ratio) else None,
                   measured_median=float(np.median(m)) if len(m) else None,
                   model_median=float(np.median(model)) if len(m) else None)
        if corr:
            dl = reference_delta(w, ref)
            ds = diff_stats(m - (model + dl))
            rlo, rhi = c_fit_ci(m - dl, w, days)
            out["ref"] = {"bias": ds["bias"], "rmse": ds["rmse"], "c_fit": fit_c(m - dl, w), "c_ci90": [rlo, rhi]}
        return out
    use = [c for c in colls if np.isfinite(c["w"])]
    out = {"all": block(use) if use else None, "by_hs": {}, "by_phase": {}, "by_day_offset": {},
           "reference_corrected": corr}
    for key, field in (("by_hs", "band"), ("by_phase", "phase"), ("by_day_offset", "day_offset")):
        for v in sorted(set(c[field] for c in use), key=str):
            out[key][str(v)] = block([c for c in use if c[field] == v])
    return out


# ---------------------------------------------------------------------
# One survey date
# ---------------------------------------------------------------------

EPOCH_RE = re.compile(r"^(\d{9,10})\.")


def frame_epoch(name):
    m = EPOCH_RE.match(str(name))
    return int(m.group(1)) if m else None


def product_frames(wfile):
    """The product's frames: {source_file: {cam, epoch, sw, hs, tp}} (evaluate_setup_model.read_lines)."""
    from evaluate_setup_model import read_lines
    L = read_lines(wfile)
    out = {}
    for i, name in enumerate(L["source_file"].tolist()):
        if name in out:
            continue
        ep = L["capture_epoch"][i]
        e2 = frame_epoch(name)
        out[name] = {"cam": str(L["camera"][i]),
                     "epoch": float(ep) if np.isfinite(ep) else (float(e2) if e2 is not None else np.nan),
                     "sw": float(L["tide_elevation_navd88"][i]), "hs": float(L["offshore_hs_m"][i]),
                     "tp": float(L["offshore_tp_s"][i])}
    return out


def waterline_file(pdir, prov):
    w = Path(str((prov.get("filter") or {}).get("file") or ""))
    if w.is_file():
        return w
    for f in ("contour_points_ground_filtered.csv", "contour_points_ground.csv"):
        if (pdir / "waterlines" / f).is_file():
            return pdir / "waterlines" / f
    return None


def window_of(prov, row, cam="c2"):
    """(first, last) UTC day of the product's window for `cam` (else any camera's; else the survey day)."""
    cams = prov.get("cameras") or {}
    for c in [cam] + sorted(cams):
        w = (cams.get(c) or {}).get("window")
        if w and len(w) == 2 and w[0] and w[1]:
            return str(w[0]), str(w[1])
    d = row.get("survey_date") or row["date"]
    return d, d


def live_level(pdir):
    """True when the product's still water is the station's GNSS-R (a live-era date: its forcing.json says so):
    the stacks then take marconi_water_level.py's GPS level, the one the archive rows were given, not the
    per-frame copy in forcing/water_level.csv."""
    try:
        f = json.loads((pdir / "forcing" / "forcing.json").read_text())
    except (OSError, ValueError):
        return False
    return str((f.get("water_level") or {}).get("method") or "").startswith("GNSS-R")


def _quantiles(v):
    v = np.asarray(v, float)
    return {"n": int(len(v)), "p10": float(np.percentile(v, 10)), "median": float(np.median(v)),
            "p90": float(np.percentile(v, 90))} if len(v) else {"n": 0}


def measure_date(row, survey_path, pdir, args, log=print):
    """Rules 2-4 for one survey row -> (result dict, measured rows) or (None, why)."""
    import timestack_wave_height as thw
    prov_p = pdir / "provenance.json"
    if not prov_p.exists():
        return None, f"no survey-date product ({prov_p} not found): build it with survey_products.py first"
    prov = json.loads(prov_p.read_text())
    if prov.get("blind"):
        return None, "BLIND product (its provenance says so): not read"
    first, last = args.window if args.window else window_of(prov, row, args.camera)
    lo, hi = ms._day_epoch(first), ms._day_epoch(last, end=True)
    stacks = [p for p in sorted(Path(args.ras_dir).glob(f"*.{args.camera}.ras.tiff"), key=lambda q: q.name)
              if ms.epoch_from_name(p.name) is not None and lo <= ms.epoch_from_name(p.name) < hi]
    if not stacks:
        return None, f"no {args.camera} stacks of {first} .. {last} in {args.ras_dir}"
    cal = ((prov.get("cameras") or {}).get(args.camera) or {}).get("calibration") or {}
    eo = args.eo or cal.get("eo_path")
    io = args.io or (str(Path(args.calibration) / cal["io_file"]) if cal.get("io_file") else None)
    from georectify import load_extrinsics, load_intrinsics
    io_, eo_ = (load_intrinsics(io or HERE / f"calibration/CACO05_{args.camera}_20240801_IO.yaml"),
                load_extrinsics(eo or HERE / f"calibration/CACO05_{args.camera}_20251113_EO-CV.yaml"))
    uv = next((u for u in (stack_lines(p, args.camera, args.pix, args.lines) for p in stacks[:10]) if u), None)
    box = line_footprint(uv, io_, eo_) if uv else None
    if box is None:
        log(f"  {row['date']}: NOTE: the lines' footprint is not known (no .pix for the stacks' width?): the whole "
            f"survey is read")
    kind = "dsm" if row["survey_type"] == "dsm" else "points"
    surf, label = survey_surface(kind, survey_path, box, args.rtk_corridor)
    if surf is None:
        return None, f"{survey_path}: no surface between {SURFACE_Z[0]:g} and {SURFACE_Z[1]:g} m near the c2 lines"
    profile = ms.FixedProfile(surf, label, ms.file_sha256(survey_path))
    wl_csv = args.level_csv or (pdir / "forcing" / "water_level.csv" if (pdir / "forcing" / "water_level.csv").is_file()
                                and not live_level(pdir) else None)
    if wl_csv:
        wl = thw.LevelSeries(wl_csv, label=LEVEL_LABEL)
        sources = (LEVEL_LABEL,)
    else:
        from marconi_water_level import WaterLevel
        wl = WaterLevel(args.gnssr_spline, args.gauge_csv)
        sources = ("gps",)
    owg_csv, sharp = args.owg_csv, {"source": "owg_live.py record " + str(args.owg_csv)}
    if args.sharpness_images:
        owg_csv = Path(args.out_dir) / f"{row['date']}_{row['name']}_sharpness.csv"
        by = image_sharpness_table(args.sharpness_images, args.camera, lo, hi, owg_csv)
        sharp = {"source": f"scored here from the period's own {args.camera} images "
                           f"({', '.join(args.sharpness_images)}) -> {owg_csv}",
                 "by_product": {k: _quantiles(v) for k, v in sorted(by.items())}}
    if args.allow_no_sharpness:
        sharp["note"] = "--allow-no-sharpness: a collection without a sharpness passes (no blur check for it)"
    margs = argparse.Namespace(camera=args.camera, pix=args.pix, io=io, eo=eo, max_extrapolate=args.max_extrapolate,
                               max_retreat_speed=None, wait_days=0.0, owg_csv=str(owg_csv), level_csv=None,
                               no_gps=False, gnssr_spline=args.gnssr_spline, gauge_csv=args.gauge_csv,
                               min_sharpness=args.min_sharpness, allow_no_sharpness=args.allow_no_sharpness,
                               max_extrapolated=args.max_extrapolated, min_lines=args.min_lines,
                               level_sources=list(sources))
    ctx = ms.make_context(margs, profile, wl)
    log(f"  {row['date']}: {len(stacks)} stack(s) of {first} .. {last}; surface: {label}; still water: "
        + (str(wl_csv) if wl_csv else "GNSS-R (marconi_water_level.py)") + f"; QC: {ctx.qc.describe()}")
    rows = []
    for p in stacks:
        r = ms.measure_stack(p, ctx, args.lines)
        if isinstance(r, dict):
            rows.append(r)
    return {"prov": prov, "window": [first, last], "surface": label, "level": str(wl_csv) if wl_csv else "gps",
            "qc": ctx.qc.describe(), "wl": wl, "stacks": len(stacks), "sharpness": sharp}, rows


def join_frames(rows, frames, level, survey_day, station_c, tol=ms.MATCH_TOL_S, ref=None):
    """The ok collections with their frames' waves (rule 4) -> (collections, counts)."""
    out, n_frames, outage, no_frame = [], 0, 0, 0
    fe = np.array([f["epoch"] for f in frames.values()], float)
    fl = list(frames.values())
    for r in rows:
        if r["status"] != "ok":
            continue
        e = float(r["collection_epoch"])
        near = [fl[i] for i in np.flatnonzero(np.abs(fe - e) <= tol)] if len(fe) else []
        if not near:
            no_frame += 1
            continue
        n_frames += len(near)
        wv = next((f for f in near if np.isfinite(f["hs"]) and np.isfinite(f["tp"]) and f["tp"] > 0), None)
        if wv is None:
            outage += 1
            continue
        mid = e + ms.BURST_MID_OFFSET_S
        day = ms.utc_day(e)
        w = float(wave_scale(wv["hs"], wv["tp"]))
        model = float(setup_c(wv["hs"], wv["tp"], station_c))
        dl = float(reference_delta([w], ref)[0])
        out.append({"epoch": e, "m": float(r["setup_m"]), "sigma": ms._num(r["setup_sigma_m"]),
                    "sigma_profile": ms._num(r.get("sigma_profile_m")),
                    "w": w, "hs": wv["hs"], "tp": wv["tp"], "model": model,
                    "model_ref": model + dl if np.isfinite(dl) else None, "band": hs_band(wv["hs"]),
                    "phase": tide_phase(level, mid), "day": str(day),
                    "day_offset": (day - survey_day).days})
    return out, {"collections_with_frames": len(out) + outage, "collections_without_frames": no_frame,
                 "collections_no_tp": outage, "frames": n_frames}


def score_lines(sv, table, station_c, phase_of=None, ref=None):
    """Rule 5: the frames scored (sv.ref, evaluate_setup_model.Survey) with the model's setup, with the model's
    corrected to C's reference (rule 6: a usable historical `ref`) and with their collection's measured one
    (table: measured_setup.MeasuredSetups); the paired differences measured - model (the corrected one
    where there is one) -> dict."""
    from evaluate_setup_model import frame_stats
    sv.set_reference(station_c)
    fr = sv.frames
    eps = [frame_epoch(f["name"]) for f in fr]
    s_mod = np.array([float(setup_c(f["hs"], f["tp"], station_c)) for f in fr])
    s_mea = np.array([np.nan if e is None or table.at(e) is None else table.at(e) for e in eps], float)
    corr = bool(ref and ref.get("usable") and ref.get("kind") == "historical")
    has = np.isfinite(s_mea)
    r_mod = sv.residuals(s_mod)
    r_mea = sv.residuals(np.where(has, s_mea, np.nan))
    r_auto = sv.residuals(np.where(has, s_mea, s_mod))
    r_ref = None
    if corr:
        w = np.array([float(wave_scale(f["hs"], f["tp"])) for f in fr])
        r_ref = sv.residuals(s_mod + reference_delta(w, ref))
    ref_idx = np.asarray(sv.ref, int)
    both = np.zeros(len(fr), bool)
    both[ref_idx] = True
    both &= has & np.isfinite(r_mod) & np.isfinite(r_mea)
    if r_ref is not None:
        both &= np.isfinite(r_ref)
    cams = np.array([f["cam"] for f in fr], object)
    bands = np.array([hs_band(f["hs"]) for f in fr], object)
    days = np.array([f.get("day", "") for f in fr], object)
    against = r_ref if r_ref is not None else r_mod
    out = {"frames": len(fr), "scored": int(len(ref_idx)), "gross_failures": sv.gross,
           "with_measured": int(has[ref_idx].sum()), "compared": int(both.sum()),
           "model": frame_stats(r_mod[both], cams[both]), "measured": frame_stats(r_mea[both], cams[both]),
           "model_ref": frame_stats(r_ref[both], cams[both]) if r_ref is not None else None,
           "paired": dict(paired_diffs(against[both], r_mea[both], days[both]),
                          against="the model corrected to C's reference" if r_ref is not None else
                          ("the model" + ("" if (ref or {}).get("kind") == "gnssr" else " (NOT corrected: "
                                          + str((ref or {}).get("note") or "no reference") + ")"))),
           "auto_all": frame_stats(r_auto[ref_idx], cams[ref_idx]), "model_all": frame_stats(r_mod[ref_idx], cams[ref_idx]),
           "by_hs": {}, "by_phase": {}}
    for b in sorted(set(bands[both].tolist())):
        m = both & (bands == b)
        out["by_hs"][b] = {"model": frame_stats(r_mod[m]), "measured": frame_stats(r_mea[m]),
                           "model_ref": frame_stats(r_ref[m]) if r_ref is not None else None}
    if phase_of:
        ph = np.array([phase_of.get(table.nearest(e), "unknown") if e is not None else "unknown" for e in eps], object)
        for p in sorted(set(ph[both].tolist())):
            m = both & (ph == p)
            out["by_phase"][p] = {"model": frame_stats(r_mod[m]), "measured": frame_stats(r_mea[m]),
                                  "model_ref": frame_stats(r_ref[m]) if r_ref is not None else None}
    return out


# ---------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------

def f3(v, signed=True):
    return "   -  " if v is None or not np.isfinite(v) else (f"{v:+.3f}" if signed else f"{v:.3f}")


def ci_text(ci, days, nd=3):
    if ci and ci[0] is not None and np.isfinite(ci[0]):
        return f"[{ci[0]:.{nd}f}, {ci[1]:.{nd}f}; {days} d]"
    return f"[{days} d: fewer than {MIN_DAYS_CI}]"


def report_date(res, station_c):
    L = [f"{res['date']} {res['name']}: {res['survey_kind']} {res['survey']}" +
         (f"  [C = {', '.join(f'{c:g}' for c in res['circular_c'])} was FITTED to this survey: the model is "
          f"CIRCULAR here]" if res.get("circular_c") else "")]
    L.append(f"  stacks {res['stacks']} of {res['window'][0]} .. {res['window'][1]}; collections measured "
             f"{res['measured']['total']} ({res['measured']['ok']} passed QC); surface: {res['surface']}; "
             f"still water: {res['level']}")
    sh = res.get("sharpness") or {}
    L.append(f"  window sharpness: {sh.get('source', '?')}" + (f"; {sh['note']}" if sh.get("note") else ""))
    for prod, q in (sh.get("by_product") or {}).items():
        if q.get("n"):
            L.append(f"      {prod:<7} {q['n']:4d} collection(s): p10 {q['p10']:.1f}, median {q['median']:.1f}, "
                     f"p90 {q['p90']:.1f} (blur limit {ms.MIN_SHARPNESS:g}, set on bright frames)")
    ref = res.get("reference") or {}
    L.append(f"  reference (rule 6): {ref.get('note', '?')}")
    for why, n in sorted(res["measured"]["fail_reasons"].items(), key=lambda kv: -kv[1])[:6]:
        L.append(f"      failed: {n:4d} line(s) {why}")
    j = res["join"]
    L.append(f"  joined to the product's frames: {j['collections_with_frames']} collection(s) ({j['frames']} frames); "
             f"{j['collections_without_frames']} without a frame; {j['collections_no_tp']} whose frames have no Tp "
             f"(the model has no setup there: H2)")
    cm = res["model_compare"]
    if cm and cm.get("all"):
        L.append(f"  MEASURED against the MODEL C x sqrt(Hs L0), C = {station_c:g} (rule 4): subset  n  days  "
                 f"measured-model bias  RMSE  C_fit [90% day jackknife; days]  median m/w")

        def line(name, b):
            return (f"      {name:<16} {b['n']:4d} {b['days']:3d}   {f3(b['bias'])}  {f3(b['rmse'], False)}  "
                    f"{b['c_fit']:.4f} {ci_text(b['c_ci90'], b['days'], 4)}  {b['c_ratio_median']:.4f}")
        L.append(line("all", cm["all"]))
        for key, title in (("by_hs", "Hs"), ("by_phase", ""), ("by_day_offset", "day")):
            for k, b in cm[key].items():
                L.append(line(f"{title} {k}".strip() if key != "by_day_offset" else f"day {int(k):+d}", b))
        if cm.get("reference_corrected"):
            L.append("    against the model CORRECTED to C's reference (rule 6): subset  measured-model bias  RMSE  "
                     "C_fit at C's reference [90%; days]")
            for name, b in [("all", cm["all"])] + [(f"Hs {k}", v) for k, v in cm["by_hs"].items()]:
                r = b["ref"]
                L.append(f"      {name:<16} {f3(r['bias'])}  {f3(r['rmse'], False)}  {r['c_fit']:.4f} "
                         f"{ci_text(r['c_ci90'], b['days'], 4)}")
    sc = res.get("lines")
    if sc:
        L.append(f"  THE LINES ON THE SURVEY (rule 5): {sc['scored']} frame(s) scored, {sc['with_measured']} with a "
                 f"measured setup, {sc['compared']} reaching the survey with each setup; per frame median of (line - "
                 f"survey):  n  bias  median  RMSE  NMAD")
        rows = [("model", sc["model"])] + ([("model+ref", sc["model_ref"])] if sc.get("model_ref") else []) + \
            [("measured", sc["measured"])]
        for k, r in rows:
            L.append(f"      {k:<9} {r['n']:4d}  {f3(r['bias'])}  {f3(r['median'])}  {f3(r['rmse'], False)}  "
                     f"{f3(r['nmad'], False)}" + "".join(f"  [{c} {v['median']:+.3f} ({v['n']})]"
                                                        for c, v in sorted((r.get('by_camera') or {}).items())))
        p = sc["paired"]
        L.append(f"      paired, measured - {p['against']}, 90% day jackknife ({p['days']} d): RMSE "
                 f"{f3(p['rmse'][0])} {ci_text(p['rmse'][1:], p['days'])}, |bias| {f3(p['absbias'][0])} "
                 f"{ci_text(p['absbias'][1:], p['days'])}")
        for b, v in sc["by_hs"].items():
            L.append(f"      Hs {b:<9} model RMSE {f3(v['model']['rmse'], False)} ({v['model']['n']})"
                     + (f", model+ref {f3(v['model_ref']['rmse'], False)}" if v.get("model_ref") else "")
                     + f", measured {f3(v['measured']['rmse'], False)}")
        for ph, v in sc["by_phase"].items():
            L.append(f"      {ph:<12} model RMSE {f3(v['model']['rmse'], False)} ({v['model']['n']})"
                     + (f", model+ref {f3(v['model_ref']['rmse'], False)}" if v.get("model_ref") else "")
                     + f", measured {f3(v['measured']['rmse'], False)}")
        a, m = sc["auto_all"], sc["model_all"]
        L.append(f"      every frame scored: model RMSE {f3(m['rmse'], False)} ({m['n']}), auto (measured where "
                 f"there is one) {f3(a['rmse'], False)} ({a['n']})")
    elif res.get("lines_why"):
        L.append(f"  the lines on the survey: not scored ({res['lines_why']})")
    return L


# ---------------------------------------------------------------------

def evaluate_row(row, survey_path, pdir, args, sp, log=print):
    """One survey row through rules 2-6 -> result dict, or None (why logged)."""
    res, rows = measure_date(row, survey_path, pdir, args, log)
    if res is None:
        log(f"  {row['date']} {row['name']}: {rows}: skipped")
        return None
    if not rows:
        log(f"  {row['date']} {row['name']}: none of the {res['stacks']} stack(s) could be read along lines "
            f"{'/'.join(str(x) for x in args.lines)} (check --pix for their width, and --lines): skipped")
        return None
    prov, wl = res.pop("prov"), res.pop("wl")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mcsv = out_dir / f"{row['date']}_{row['name']}_setup_measured_survey.csv"
    if mcsv.exists():
        mcsv.unlink()
    ms.write_rows(mcsv, rows)
    fails = ms.fail_reasons(rows)
    res.update(date=row["date"], name=row["name"], survey=str(survey_path), survey_kind=row["survey_type"],
               measured={"total": len(rows), "ok": sum(r["status"] == "ok" for r in rows), "fail_reasons": fails,
                         "csv": str(mcsv)})
    sha = sp.sha256(survey_path)
    res["circular_c"] = [c for c, f in sorted(sp.SETUP_FITS.items())
                         if f.get("survey_sha256") == sha or Path(str(survey_path)).stem == f.get("survey_stem")]
    ref = reference_of(pdir)
    res["reference"] = ref
    survey_day = datetime.strptime(row.get("survey_date") or row["date"], "%Y-%m-%d").date()
    wfile = waterline_file(pdir, prov)
    frames = product_frames(wfile) if wfile else {}
    colls, counts = join_frames(rows, frames, wl, survey_day, args.station_c, ref=ref)
    res["join"] = counts
    res["model_compare"] = compare_model(colls, args.station_c, ref) if colls else None
    res["collections"] = colls
    table = ms.MeasuredSetups(mcsv)
    phase_of = {c["epoch"]: c["phase"] for c in colls}
    try:
        from evaluate_setup_model import load_product
        sv, why = load_product(pdir, row, survey_path, args.calibration)
    except (OSError, ValueError, KeyError) as exc:
        sv, why = None, str(exc)
    if sv is not None and table.ok:
        res["lines"] = score_lines(sv, table, args.station_c, phase_of, ref)
    else:
        res["lines_why"] = why or "no collection passed QC"
    return res


def main(argv=None):
    import survey_products as sp
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--output-root", default=str(sp.DEFAULT_OUTPUT_ROOT), help="the survey-date products")
    ap.add_argument("--surveys", default=str(sp.DEFAULT_SURVEYS), help="surveys/surveys.csv")
    ap.add_argument("--survey-dirs", nargs="+", default=list(sp.DEFAULT_SURVEY_DIRS),
                    help="where the survey files are looked for by name (survey_products.py's)")
    ap.add_argument("--calibration", default=str(HERE / "calibration"), help="the lens (IO) files")
    ap.add_argument("--dates", nargs="+", default=None, help="only these product dates")
    ap.add_argument("--station-c", type=float, default=sp.SETUP_COEF, help=f"the model's C (default {sp.SETUP_COEF})")
    ap.add_argument("--ras-dir", default=str(HERE / "archive" / "ras_c2"), help="the c2 runup stacks")
    ap.add_argument("--camera", default="c2")
    ap.add_argument("--lines", nargs="+", type=int, default=list(ms.DEFAULT_LINES))
    ap.add_argument("--pix", default=None, help="the stacks' pixel list (default: runup_from_timestack.py's)")
    ap.add_argument("--io", default=None, help="default: the product's c2 lens file")
    ap.add_argument("--eo", default=None, help="default: the product's c2 calibration")
    ap.add_argument("--window", nargs=2, default=None, metavar=("FIRST", "LAST"),
                    help="stack days (default: the product's c2 window)")
    ap.add_argument("--owg-csv", default=str(HERE / "archive" / "owg_c2_H.csv"))
    ap.add_argument("--sharpness-images", nargs="+", default=None, metavar="DIR",
                    help="score the window's sharpness from the period's own c2 images in these folders (bright, "
                         "else snap, else timex) instead of --owg-csv: for a period without the OWG record (2025)")
    ap.add_argument("--min-sharpness", type=float, default=ms.MIN_SHARPNESS)
    ap.add_argument("--allow-no-sharpness", action="store_true",
                    help="pass lines whose collection has no sharpness at all (no blur check: the last resort)")
    ap.add_argument("--max-extrapolated", type=float, default=ms.MAX_EXTRAPOLATED)
    ap.add_argument("--max-extrapolate", type=float, default=10.0, help="m past the surface's edge (10)")
    ap.add_argument("--min-lines", type=int, default=ms.MIN_LINES)
    ap.add_argument("--level-csv", default=None, help="still water (default: the product's forcing, else GNSS-R)")
    from marconi_water_level import GNSSR_SPLINE, GAUGE_CSV
    ap.add_argument("--gnssr-spline", default=GNSSR_SPLINE)
    ap.add_argument("--gauge-csv", default=GAUGE_CSV)
    ap.add_argument("--rtk-corridor", type=float, default=RTK_CORRIDOR_M)
    ap.add_argument("--out-dir", default=str(HERE / "archive" / "measured_setup_eval"),
                    help="the measured setups per date (CSV)")
    ap.add_argument("--json", default=None, help="write every number here")
    args = ap.parse_args(argv)
    sp.arm_blind_guard()       # no file named like a held survey (BLIND_DATES) is opened by this run
    rows = sp.read_table(args.surveys, sp.SURVEY_FIELDS)
    root = Path(args.output_root)
    print(f"surveys table: {args.surveys} ({len(rows)} rows); products: {root}; stacks: {args.ras_dir}")
    results = []
    for row in rows:
        tag = f"{row['date']} {row['name']}"
        if args.dates and row["date"] not in args.dates:
            continue
        if sp.held_survey(row):
            print(f"  {tag}: BLIND (survey_products.BLIND_DATES): skipped, nothing of it opened")
            continue
        if row["survey_type"] not in ("dsm", "points"):
            print(f"  {tag}: survey type {row['survey_type']!r}: not used")
            continue
        if not row.get("path"):
            print(f"  {tag}: no path in {Path(args.surveys).name}: skipped (never looked up by its date)")
            continue
        if sp.held_name(Path(row["path"]).name):
            print(f"  {tag}: {Path(row['path']).name} is named like a survey of a BLIND date: skipped, never opened")
            continue
        hits = sp.resolve_file(row["path"], args.survey_dirs)
        if not hits:
            print(f"  {tag}: {row['path']} not found in {', '.join(args.survey_dirs)}: skipped")
            continue
        res = evaluate_row(row, hits[0], root / row["date"], args, sp)
        if res is not None:
            results.append(res)
    print("=" * 100)
    print("THE MEASURED SETUP ON THE SURVEYS (evaluate_measured_setup.py): measured on the survey's surface, "
          "against the model")
    print("=" * 100)
    if not results:
        print("No survey date could be evaluated (no product, no stacks of its window or none readable, or no "
              "survey file).")
        return 4
    for res in results:
        for line in report_date(res, args.station_c):
            print(line)
        print()
    print("An upper bound: measured on the true surface, the live setup inherits its DEM's error too (rule 2 of "
          "measured_setup.py).")
    print("Decision rule: docs/VALIDATION.md 4.12 (the station stays on the model until it is met).")
    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=1, default=str) + "\n")
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

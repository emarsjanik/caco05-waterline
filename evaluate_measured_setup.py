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
and the surveys are there. It writes a report, --json and one CSV of the
measured setups per date in --out-dir; nothing else.

WHY. Live, the measured setup stands on a DEM of earlier days, so it
inherits that DEM's level (and the model's setup in it) and reads beach
change since as setup (measured_setup.py, rule 2). On the survey's own
surface neither applies: what is left is the measurement itself (the
detector, the still water) and how the timex line relates to the mean
swash edge. This is the clean test; the live trial (docs/VALIDATION.md
4.12) comes after it.

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
    product's lines were given. Without the OWG record (2025), the window
    sharpness check needs --allow-no-sharpness (the report says so).
 4. AGAINST THE MODEL. Each collection that passes is joined to the
    product's frames whose file-name epoch lies within 300 s of its start
    (the extractor's own rule); their Hs and Tp give the wave scale w =
    sqrt(Hs L0) and the model's setup C w (setup_model.setup_c, --station-c
    0.037). The C that would match the measured setups is the least-squares
    fit through the origin, C_fit = sum(m w) / sum(w^2), with a 90% interval
    from a DAY-BLOCK bootstrap (BOOTSTRAP resamples of whole days: the
    collections of a day share their beach, level and waves, so treating
    them as independent would make the interval too narrow; fewer than 3
    days: no interval). Reported overall, by Hs band (HS_BANDS of the frames'
    Hs), by the tide's phase (flood / ebb / slack: the still water 15 min
    after the burst's middle less 15 min before, +/- TIDE_SLACK_M) and by the
    day's offset from the survey (beach change since the survey moves the
    measurement too). Collections the model cannot serve (no Tp: the buoy
    outages of finding H2) are counted.
 5. ON THE SURVEY (the DEM-relevant score). The product's frames as
    evaluate_setup_model.py scores them (load_product: the DEM's frames with
    Hs and Tp, each line re-projected onto the plane still water + its
    setup, the survey-floor rule; set_reference: frames more than 1 m off
    at the station's C are detection failures, left out), each frame with a
    measured setup scored twice: with the model's setup and with its
    collection's measured one, on the frames where both reach the survey.
    Per frame the median of (line - survey); bias, median, RMSE and NMAD,
    per camera, by Hs band and by phase. 'auto' (measured where there is
    one, the model elsewhere) is scored on every frame. The 2026-09-29 RTK
    was the survey C was FITTED to (survey_products.SETUP_FITS): there the
    model is circular and the comparison is flagged; the 2025 lidars are
    the independent dates.

WHAT IT CANNOT DO
  * Stacks of the survey dates of 2025 are not on the NUC until fetched
    from S3; the 2026 RTK covers only the upper beach (~+1.0 to +1.6 m).
  * The surface is the beach on the survey day: days before and after
    differ by the beach's change (by day offset, 4.).
  * The detector's own bias on the synthetic scenes (+0.04 to +0.10 m:
    measured_setup.py) is in the measured setups; still water is the
    product's (GNSS-R holds part of the setup).
  * Frames without Tp are not in the score of 5 (load_product keeps frames
    with waves): the outages are counted, not scored.

Usage (on the NUC, after survey_products.py has built the dates):
    python3 evaluate_measured_setup.py --dates 2026-09-29
    python3 evaluate_measured_setup.py --dates 2025-01-23 --ras-dir /path/to/jan2025_c2_stacks \\
        --pix /path/to/c2_timestack.pix.2025 --allow-no-sharpness --json /tmp/msetup_2025-01-23.json
Tests: python3 test_measured_setup.py
"""

import argparse
import json
import re
import sys
from datetime import datetime
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
BOOTSTRAP = 2000
MIN_DAYS_CI = 3
LEVEL_LABEL = "forcing water level"


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
# Statistics
# ---------------------------------------------------------------------

def fit_c(m, w):
    """The C of m = C w, least squares through the origin (NaN without data)."""
    m, w = np.asarray(m, float), np.asarray(w, float)
    ok = np.isfinite(m) & np.isfinite(w)
    if not ok.any() or not (w[ok] ** 2).sum() > 0:
        return np.nan
    return float((m[ok] * w[ok]).sum() / (w[ok] ** 2).sum())


def day_block_ci(m, w, days, level=0.90, n_boot=BOOTSTRAP, seed=0, min_days=MIN_DAYS_CI):
    """The `level` interval of fit_c from resampling whole days with replacement; (NaN, NaN) with fewer than
    min_days days."""
    m, w, days = np.asarray(m, float), np.asarray(w, float), np.asarray(days, object)
    ud = sorted(set(days.tolist()))
    if len(ud) < min_days:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    idx = {d: np.flatnonzero(days == d) for d in ud}
    out = []
    for _ in range(int(n_boot)):
        pick = rng.integers(0, len(ud), len(ud))
        sel = np.concatenate([idx[ud[k]] for k in pick])
        out.append(fit_c(m[sel], w[sel]))
    out = np.asarray(out, float)
    out = out[np.isfinite(out)]
    a = (1.0 - level) / 2.0
    return float(np.quantile(out, a)), float(np.quantile(out, 1.0 - a))


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


def compare_model(colls, station_c):
    """Rule 4 over the collections joined to their frames: [{m, w, hs, tp, phase, day_offset, ...}] ->
    dict of the fit, its interval and the differences, overall and by subset."""
    def block(sub):
        m = np.array([c["m"] for c in sub], float)
        w = np.array([c["w"] for c in sub], float)
        days = np.array([c["day"] for c in sub], object)
        cf = fit_c(m, w)
        lo, hi = day_block_ci(m, w, days)
        model = station_c * w
        ratio = m / w
        return dict(diff_stats(m - model), c_fit=cf, c_ci90=[lo, hi], days=len(set(days.tolist())),
                    c_ratio_median=float(np.median(ratio)) if len(ratio) else None,
                    measured_median=float(np.median(m)) if len(m) else None,
                    model_median=float(np.median(model)) if len(m) else None)
    use = [c for c in colls if np.isfinite(c["w"])]
    out = {"all": block(use) if use else None, "by_hs": {}, "by_phase": {}, "by_day_offset": {}}
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
    margs = argparse.Namespace(camera=args.camera, pix=args.pix, io=io, eo=eo, max_extrapolate=args.max_extrapolate,
                               max_retreat_speed=None, wait_days=0.0, owg_csv=args.owg_csv, level_csv=None,
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
            "qc": ctx.qc.describe(), "wl": wl, "stacks": len(stacks)}, rows


def join_frames(rows, frames, level, survey_day, station_c, tol=ms.MATCH_TOL_S):
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
        out.append({"epoch": e, "m": float(r["setup_m"]), "sigma": ms._num(r["setup_sigma_m"]),
                    "w": float(wave_scale(wv["hs"], wv["tp"])), "hs": wv["hs"], "tp": wv["tp"],
                    "model": float(setup_c(wv["hs"], wv["tp"], station_c)), "band": hs_band(wv["hs"]),
                    "phase": tide_phase(level, mid), "day": str(day),
                    "day_offset": (day - survey_day).days})
    return out, {"collections_with_frames": len(out) + outage, "collections_without_frames": no_frame,
                 "collections_no_tp": outage, "frames": n_frames}


def score_lines(sv, table, station_c, phase_of=None):
    """Rule 5: the frames scored (sv.ref, evaluate_setup_model.Survey) with the model's setup and with their
    collection's measured one (table: measured_setup.MeasuredSetups) -> dict."""
    from evaluate_setup_model import frame_stats
    sv.set_reference(station_c)
    fr = sv.frames
    eps = [frame_epoch(f["name"]) for f in fr]
    s_mod = np.array([float(setup_c(f["hs"], f["tp"], station_c)) for f in fr])
    s_mea = np.array([np.nan if e is None or table.at(e) is None else table.at(e) for e in eps], float)
    has = np.isfinite(s_mea)
    r_mod = sv.residuals(s_mod)
    r_mea = sv.residuals(np.where(has, s_mea, np.nan))
    r_auto = sv.residuals(np.where(has, s_mea, s_mod))
    ref = np.asarray(sv.ref, int)
    both = np.zeros(len(fr), bool)
    both[ref] = True
    both &= has & np.isfinite(r_mod) & np.isfinite(r_mea)
    cams = np.array([f["cam"] for f in fr], object)
    bands = np.array([hs_band(f["hs"]) for f in fr], object)
    out = {"frames": len(fr), "scored": int(len(ref)), "gross_failures": sv.gross,
           "with_measured": int(has[ref].sum()), "compared": int(both.sum()),
           "model": frame_stats(r_mod[both], cams[both]), "measured": frame_stats(r_mea[both], cams[both]),
           "auto_all": frame_stats(r_auto[ref], cams[ref]), "model_all": frame_stats(r_mod[ref], cams[ref]),
           "by_hs": {}, "by_phase": {}}
    for b in sorted(set(bands[both].tolist())):
        m = both & (bands == b)
        out["by_hs"][b] = {"model": frame_stats(r_mod[m]), "measured": frame_stats(r_mea[m])}
    if phase_of:
        ph = np.array([phase_of.get(table.nearest(e), "unknown") if e is not None else "unknown" for e in eps], object)
        for p in sorted(set(ph[both].tolist())):
            m = both & (ph == p)
            out["by_phase"][p] = {"model": frame_stats(r_mod[m]), "measured": frame_stats(r_mea[m])}
    return out


# ---------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------

def f3(v, signed=True):
    return "   -  " if v is None or not np.isfinite(v) else (f"{v:+.3f}" if signed else f"{v:.3f}")


def report_date(res, station_c):
    L = [f"{res['date']} {res['name']}: {res['survey_kind']} {res['survey']}" +
         (f"  [C = {', '.join(f'{c:g}' for c in res['circular_c'])} was FITTED to this survey: the model is "
          f"CIRCULAR here]" if res.get("circular_c") else "")]
    L.append(f"  stacks {res['stacks']} of {res['window'][0]} .. {res['window'][1]}; collections measured "
             f"{res['measured']['total']} ({res['measured']['ok']} passed QC); surface: {res['surface']}; "
             f"still water: {res['level']}")
    for why, n in sorted(res["measured"]["fail_reasons"].items(), key=lambda kv: -kv[1])[:6]:
        L.append(f"      failed: {n:4d} line(s) {why}")
    j = res["join"]
    L.append(f"  joined to the product's frames: {j['collections_with_frames']} collection(s) ({j['frames']} frames); "
             f"{j['collections_without_frames']} without a frame; {j['collections_no_tp']} whose frames have no Tp "
             f"(the model has no setup there: H2)")
    cm = res["model_compare"]
    if cm and cm.get("all"):
        L.append(f"  MEASURED against the MODEL C x sqrt(Hs L0), C = {station_c:g} (rule 4): subset  n  days  "
                 f"measured-model bias  RMSE  C_fit [90% day-block]  median m/w")

        def line(name, b):
            ci = b["c_ci90"]
            cis = (f"[{ci[0]:.3f}, {ci[1]:.3f}]" if ci and np.isfinite(ci[0]) else "[fewer than "
                   f"{MIN_DAYS_CI} days]")
            return (f"      {name:<16} {b['n']:4d} {b['days']:3d}   {f3(b['bias'])}  {f3(b['rmse'], False)}  "
                    f"{b['c_fit']:.4f} {cis}  {b['c_ratio_median']:.4f}")
        L.append(line("all", cm["all"]))
        for key, title in (("by_hs", "Hs"), ("by_phase", ""), ("by_day_offset", "day")):
            for k, b in cm[key].items():
                L.append(line(f"{title} {k}".strip() if key != "by_day_offset" else f"day {int(k):+d}", b))
    sc = res.get("lines")
    if sc:
        L.append(f"  THE LINES ON THE SURVEY (rule 5): {sc['scored']} frame(s) scored, {sc['with_measured']} with a "
                 f"measured setup, {sc['compared']} reaching the survey with both; per frame median of (line - "
                 f"survey):  n  bias  median  RMSE  NMAD")
        for k in ("model", "measured"):
            r = sc[k]
            L.append(f"      {k:<9} {r['n']:4d}  {f3(r['bias'])}  {f3(r['median'])}  {f3(r['rmse'], False)}  "
                     f"{f3(r['nmad'], False)}" + "".join(f"  [{c} {v['median']:+.3f} ({v['n']})]"
                                                        for c, v in sorted((r.get('by_camera') or {}).items())))
        for b, v in sc["by_hs"].items():
            L.append(f"      Hs {b:<9} model RMSE {f3(v['model']['rmse'], False)} ({v['model']['n']}), measured "
                     f"{f3(v['measured']['rmse'], False)}")
        for p, v in sc["by_phase"].items():
            L.append(f"      {p:<12} model RMSE {f3(v['model']['rmse'], False)} ({v['model']['n']}), measured "
                     f"{f3(v['measured']['rmse'], False)}")
        a, m = sc["auto_all"], sc["model_all"]
        L.append(f"      every frame scored: model RMSE {f3(m['rmse'], False)} ({m['n']}), auto (measured where "
                 f"there is one) {f3(a['rmse'], False)} ({a['n']})")
    elif res.get("lines_why"):
        L.append(f"  the lines on the survey: not scored ({res['lines_why']})")
    return L


# ---------------------------------------------------------------------

def evaluate_row(row, survey_path, pdir, args, sp, log=print):
    """One survey row through rules 2-5 -> result dict, or None (why logged)."""
    res, rows = measure_date(row, survey_path, pdir, args, log)
    if res is None:
        log(f"  {row['date']} {row['name']}: {rows}: skipped")
        return None
    prov, wl = res.pop("prov"), res.pop("wl")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mcsv = out_dir / f"{row['date']}_{row['name']}_setup_measured_survey.csv"
    if mcsv.exists():
        mcsv.unlink()
    ms.write_rows(mcsv, rows)
    fails = {}
    for r in rows:
        for part in str(r["line_qc"]).split(";"):
            why = part.split(":", 1)[1] if ":" in part else part
            if why and why != "ok":
                key = re.sub(r"[-+]?\d+\.\d+|\b\d+\b", "#", why)
                fails[key] = fails.get(key, 0) + 1
    res.update(date=row["date"], name=row["name"], survey=str(survey_path), survey_kind=row["survey_type"],
               measured={"total": len(rows), "ok": sum(r["status"] == "ok" for r in rows), "fail_reasons": fails,
                         "csv": str(mcsv)})
    sha = sp.sha256(survey_path)
    res["circular_c"] = [c for c, f in sorted(sp.SETUP_FITS.items())
                         if f.get("survey_sha256") == sha or Path(str(survey_path)).stem == f.get("survey_stem")]
    survey_day = datetime.strptime(row.get("survey_date") or row["date"], "%Y-%m-%d").date()
    wfile = waterline_file(pdir, prov)
    frames = product_frames(wfile) if wfile else {}
    colls, counts = join_frames(rows, frames, wl, survey_day, args.station_c)
    res["join"] = counts
    res["model_compare"] = compare_model(colls, args.station_c) if colls else None
    res["collections"] = colls
    table = ms.MeasuredSetups(mcsv)
    phase_of = {c["epoch"]: c["phase"] for c in colls}
    try:
        from evaluate_setup_model import load_product
        sv, why = load_product(pdir, row, survey_path, args.calibration)
    except (OSError, ValueError, KeyError) as exc:
        sv, why = None, str(exc)
    if sv is not None and table.ok:
        res["lines"] = score_lines(sv, table, args.station_c, phase_of)
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
    ap.add_argument("--min-sharpness", type=float, default=ms.MIN_SHARPNESS)
    ap.add_argument("--allow-no-sharpness", action="store_true",
                    help="pass lines whose collection has no OWG frame (2025: no OWG record)")
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
        print("No survey date could be evaluated (no product, no stacks of its window, or no survey file).")
        return 4
    for res in results:
        for line in report_date(res, args.station_c):
            print(line)
        print()
    print("Decision rule: docs/VALIDATION.md 4.12 (the station stays on the model until it is met).")
    if args.json:
        Path(args.json).write_text(json.dumps(results, indent=1, default=str) + "\n")
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

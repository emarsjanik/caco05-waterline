#!/usr/bin/env python3
"""
Which Setup Model? C Against Stockdon, On The Surveys
=======================================================
Refits the wave-setup coefficient of the station's model (C, in
C * sqrt(Hs * L0)) and of the opt-in slope-aware model (Stockdon's k, in
k * beta_f * sqrt(H0 * L0), setup_model.py) against the RTK and lidar
surveys, on the survey-date products survey_products.py built on the
station computer, and says which model to run -- with leave-one-survey-out
cross-validation, so that no model is scored on a survey its coefficient
was fitted to.

WHY. C = 0.037 was fitted to the 29 Sep 2026 RTK shots (survey_products.py
SETUP_FITS), and every refit since was made on one week of lines
(docs/VALIDATION.md 4.1). Stockdon et al. (2006, Eq. 10) put the foreshore
slope into the setup; the beach faces of the survey dates differ (the
2025 lidars' profiles, through the fixture DEMs rendered from them: Jan
~0.13, Mar ~0.07; the real Sep-Oct 2026 DEM 0.09-0.125: VALIDATION 2.3,
4.1), so if the lines follow Stockdon a single C must be wrong one way or
the other on some of them.
Which model is right HERE is an empirical question: video shoreline
elevation models carry a site- and detector-specific wave term that has
to be calibrated against surveys, not taken from the literature (Plant et
al., 2007; Aarninkhof et al., 2003; the review behind this repository:
"calibrate gamma for each site and detector against surveys rather than
taking literature values"). And a coefficient fitted to a survey says
nothing about the model on that survey: only a survey held out does.

WHAT IS COMPARED. For each frame (one detected line) of each survey-date
product: its pixels, its still-water level (tide_elevation_navd88), its
offshore Hs and Tp, its camera's calibration (the product's provenance).
For a setup s the line lies on the plane z = still water + s, and its
ground position is where each pixel's ray meets that plane (georectify.
pixel_to_ground). A higher plane moves the line LANDWARD along the rays,
onto higher beach, so a line keeps only ~0.6-0.8 of its setup in the DEM
(VALIDATION 4.1, re-projection): the comparison re-projects every line at
every setup it scores, as a build with that setup would place it. (The
pinhole ray is straight, so the ground point at any plane lies on the line
through its points at two planes, still water and still water + 1 m: two
pixel_to_ground calls per camera, then exact interpolation, which the
test checks against pixel_to_ground at a third plane.) The frame's
residual is the median of (z - survey) over its points:
  * a lidar DSM: the survey sampled at each point (compare_dem_survey.
    sample), the points at least SURVEY_FLOOR_MARGIN above the survey's
    local floor (its lowest elevation within ~10 m: FLOOR_BLOCK_M blocks
    and their neighbours), at least MIN_POINTS_DSM of them;
  * RTK transects: survey_compare.waterlines_on_transects (compare_rtk.py's
    check), the frame's value over its points at least SURVEY_FLOOR_MARGIN
    above the transect's lowest shot ('d_cov').
The floor rule is survey_compare.py's: a line whose plane lies below where
the survey stops can only read low, and a line that would read high there
lands seaward of the survey and drops out, so a fit over every crossing is
biased (the 2026 RTK: C 0.067 over all crossings, 0.043-0.046 above the
floor; VALIDATION 4.1). The DSM's local floor applies the same rule at the
lidar's water edge.

THE FRAMES SCORED. Per survey, the frames of the product's DEM (its wave
filter and day test: survey_compare.dem_frame_selection, from the DEM's
_info.json) that have Hs and Tp and a finite residual at the station's
setup (C = 0.037), less the GROSS failures: a frame more than GROSS_M off
the survey at the station's setup is a detection failure (snow, a wrong
edge), not a setup error; it is counted and left out of every score. Every
model and coefficient is then scored on THOSE frames; one that loses more
than 1 - MIN_KEEP of them (its lines leave the survey) is not admissible
on that survey.

THE MODELS.
  fixed (nothing fitted here):
    C = 0                     no setup (the floor of every comparison)
    C = 0.037                 the station's (CIRCULAR on the survey it was
                              fitted to: SETUP_FITS, by content or name)
    Stockdon k = 0.35         the paper's coefficient, H0 = Hs / Ks (deep
                              water, from 21 m): no fit at all
  refitted (one coefficient, a grid search):
    C                         over C_GRID
    Stockdon k, H0 = Hs / Ks  over K_GRID, beta_f measured per frame on
                              the product's own DEM as the extractor's
                              --setup-model stockdon does (setup_model.
                              stockdon_frame: the swash band, three
                              passes, C / k where the DEM cannot serve)
    Stockdon k, H0 = Hs       the same, not deshoaled (setup_model.py, THE
                              DESHOALING QUESTION: a refitted k absorbs
                              most of the deshoaling; this says how much)
A refit minimises the mean over its surveys of the mean square of their
frames' residuals (each survey weighs the same, however many frames it
has). Its numbers: the coefficient fitted on every survey (in-sample), the
coefficient each survey alone would pick (how stable it is), and LEAVE-
ONE-SURVEY-OUT: for each survey, the coefficient fitted on the OTHER
surveys, scored on it. Per survey: frames, bias (mean residual), median,
RMSE, NMAD, and the median per camera. A coefficient on the edge of its
grid is flagged (the optimum may lie beyond it).

A FREE OFFSET PER SURVEY AND CAMERA. Every refit is also made, and scored
held out, on the residuals less their mean per survey and camera (over at
least MIN_GROUP frames): what is left is how a line's residual changes with
the waves within its survey, as dem_from_contours.py --fit-setup reads the
setup off repeat crossings rather than off a survey. It is needed because
with a few surveys of one dominant slope each, C and Stockdon differ mostly
by a scale on each survey's setup: a level offset of a survey or camera
(the pointing, the still-water datum -- the ADCP's for 2025, GNSS-R's for
2026 -- a camera's bias) that lines up with the slopes is credited to
Stockdon by the absolute score. An offset does not touch the free-offset
score; a wrong wave dependence does.

THE RECOMMENDATION. With 3 or more surveys, Stockdon (k refitted, H0 = Hs
/ Ks: the form fixed beforehand; H0 = Hs is a sensitivity) is recommended
only when, against the C refit, each at its leave-one-survey-out
coefficient and over the frames both keep on the held-out survey:
  * its RMSE, averaged over the held-out surveys, is at least MARGIN_M
    below C's AND it is lower on a majority of the held-out surveys; and
  * with a free offset per survey and camera (each model refitted so),
    its RMSE averaged over the held-out surveys is below C's AND it is
    lower on a majority of them;
a held-out survey on which Stockdon keeps fewer than MIN_KEEP of the
frames scored counts against it. Otherwise the station keeps C. With
fewer surveys no recommendation is made (the station keeps C). The C the
cross-validation favours is printed as information: changing C means
changing it in all three places (the cron, survey_products.py,
historical_forcing.py) with a SETUP_FITS entry (VALIDATION 4.1), not here.

THE PRODUCTS. Each survey row of surveys/surveys.csv whose date has a
product in --output-root (<root>/<date>/provenance.json: its cameras'
calibrations, its waterline file, its DEM and the DEM's _info.json). A
BLIND date (survey_products.BLIND_DATES, e.g. 2024-10-23, whose survey is
held for an independent check) is skipped before anything of it is opened:
its survey row is never resolved, its product never read; so is a row of
any date whose survey file is NAMED like a survey of a blind date
(survey_products.held_name), and survey_products.arm_blind_guard refuses,
in this process, any Python-level open of such a file. A survey row
without a path is skipped too: this script never looks a survey up by its
date, the way survey_products.py looks for GCP files.

OUTPUT. A report on stdout; --json PATH writes every number. Exit 0 with a
report (whatever it recommends), 4 when no survey could be scored.

RUN TIME. ~5 s for the three synthetic products of the test here. On the
station most of it is reading the products' waterline files (~400k rows a
week) and the lidar GeoTIFFs; not yet run there.

Usage (on the NUC, after the survey-date products are built):
    python3 evaluate_setup_model.py
    python3 evaluate_setup_model.py --output-root /mnt/I2Rgus_Data/survey_products --json setup_eval.json
    python3 evaluate_setup_model.py --dates 2025-01-23 2026-09-29

Self-test: python3 test_evaluate_setup_model.py
"""

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from setup_model import STOCKDON_SETUP_K, ForeshoreSlope, deshoal_depth, setup_c, stockdon_frame  # noqa: E402

# The coefficient grids of the refits (C: the station's 0.037 is on it; k: the paper's 0.35 is on it)
C_GRID = np.round(np.arange(0.0, 0.08005, 0.001), 3)
K_GRID = np.round(np.arange(0.05, 0.80005, 0.005), 3)
# The frames scored (THE FRAMES SCORED above)
GROSS_M = 1.0
MIN_KEEP = 0.8
MIN_FRAMES = 5                 # a survey with fewer scored frames takes no part
MIN_GROUP = 3                  # a survey x camera with fewer takes no part in the free-offset scores
# A DSM's local floor: the survey's lowest elevation over FLOOR_BLOCK_M blocks and their eight neighbours
# (~10 m about a point); the margin is survey_compare.SURVEY_FLOOR_MARGIN (read at run time)
FLOOR_BLOCK_M = 5.0
MIN_POINTS_DSM = 20
# Points per frame on a DSM (evenly spaced along the line; the frame's median barely moves): the cost of
# scoring ~400 setups per survey. On transects every point near a transect is kept (the check needs 3
# within 2 m of the line).
MAX_POINTS_PER_FRAME = 400
NEAR_TRANSECT_M = 30.0
# The recommendation (THE RECOMMENDATION above)
MIN_SURVEYS = 3
MARGIN_M = 0.01

FAMILIES = {
    "c": "C * sqrt(Hs L0) (the station's model)",
    "stockdon": "Stockdon k * beta_f * sqrt(H0 L0), H0 = Hs / Ks (deshoaled from the currency's depth)",
    "stockdon_h0hs": "Stockdon k * beta_f * sqrt(H0 L0), H0 = Hs (not deshoaled)",
}


# ---------------------------------------------------------------------
# A survey and its lines
# ---------------------------------------------------------------------

def per_frame_median(values, idx, n, min_n):
    """The median of `values` per frame index (0..n-1) over the finite ones, NaN with fewer than min_n."""
    out = np.full(n, np.nan)
    ok = np.isfinite(values)
    if not ok.any():
        return out
    v, i = values[ok], idx[ok]
    o = np.lexsort((v, i))
    v, i = v[o], i[o]
    cuts = np.flatnonzero(np.diff(i)) + 1
    for grp_v, grp_i in zip(np.split(v, cuts), np.split(i, cuts)):
        if len(grp_v) >= min_n:
            out[grp_i[0]] = float(np.median(grp_v))
    return out


class DsmFloor:
    """The DSM's local floor: its lowest elevation per FLOOR_BLOCK_M block, then over each block and its
    eight neighbours."""

    def __init__(self, grid, x0, y_top, cell, block_m=FLOOR_BLOCK_M):
        b = max(1, int(round(block_m / cell)))
        nr, nc = -(-grid.shape[0] // b), -(-grid.shape[1] // b)
        pad = np.full((nr * b, nc * b), np.nan)
        pad[:grid.shape[0], :grid.shape[1]] = grid
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)          # all-NaN blocks
            mins = np.nanmin(pad.reshape(nr, b, nc, b).transpose(0, 2, 1, 3).reshape(nr, nc, b * b), axis=2)
        p = np.full((nr + 2, nc + 2), np.nan)
        p[1:-1, 1:-1] = mins
        self.floor = np.fmin.reduce([p[1 + dr:1 + dr + nr, 1 + dc:1 + dc + nc]
                                     for dr in (-1, 0, 1) for dc in (-1, 0, 1)])
        self.x0, self.y_top, self.size = x0, y_top, b * cell

    def at(self, E, N):
        c = np.floor((E - self.x0) / self.size)
        r = np.floor((self.y_top - N) / self.size)
        ok = np.isfinite(c) & np.isfinite(r) & (c >= 0) & (r >= 0) & (c < self.floor.shape[1]) & \
            (r < self.floor.shape[0])
        out = np.full(E.shape, np.nan)
        out[ok] = self.floor[r[ok].astype(int), c[ok].astype(int)]
        return out


class Survey:
    """
    One survey and the lines of its date's product, ready to be scored at any per-frame setup.
    frames: dicts {name, cam, sw, hs, tp, tag, day}; the points' frame index `fi`, their ground position
    at still water (E0, N0) and its change per metre of plane (dE, dN).
    kind 'dsm': grid (rows north to south), x0, y_top, cell and its DsmFloor; 'transects': the
    compare_rtk.load_transects profiles.
    """

    def __init__(self, name, date, kind, label, path, frames, fi, E0, N0, dE, dN, slope, survey, notes):
        self.name, self.date, self.kind, self.label, self.path = name, date, kind, label, path
        self.frames, self.fi = frames, fi
        self.E0, self.N0, self.dE, self.dN = E0, N0, dE, dN
        self.slope, self.survey, self.notes = slope, survey, notes
        self.sw_pt = np.array([f["sw"] for f in frames], float)[fi] if len(frames) else np.zeros(0)
        self.fname_pt = np.array([f["name"] for f in frames], object)[fi] if len(frames) else np.zeros(0, object)
        self.cam_pt = np.array([f["cam"] for f in frames], object)[fi] if len(frames) else np.zeros(0, object)
        self.day_pt = np.array([f["day"] for f in frames], object)[fi] if len(frames) else np.zeros(0, object)
        self.circular_c = ()            # the fixed C values fitted to this survey (SETUP_FITS)
        self.ref = None                 # the frames scored (indices), set by set_reference()
        self.gross = 0
        self._cache = {}

    def residuals(self, setups):
        """Per frame: the median of (plane - survey) over its points, the line re-projected at
        still water + its setup; NaN where the frame does not reach the survey (see the docstring)."""
        from survey_compare import SURVEY_FLOOR_MARGIN, TRANSECT_MIN_POINTS, TRANSECT_TOLERANCE
        s = np.asarray(setups, float)[self.fi]
        E, N, Z = self.E0 + s * self.dE, self.N0 + s * self.dN, self.sw_pt + s
        n = len(self.frames)
        if self.kind == "dsm":
            from compare_dem_survey import sample
            grid, x0, y_top, cell, floor = self.survey
            zs = sample(grid, x0, y_top, cell, E, N)
            with np.errstate(invalid="ignore"):
                ok = np.isfinite(zs) & np.isfinite(Z) & (Z >= floor.at(E, N) + SURVEY_FLOOR_MARGIN)
            return per_frame_median(np.where(ok, Z - zs, np.nan), self.fi, n, MIN_POINTS_DSM)
        from survey_compare import waterlines_on_transects
        good = np.isfinite(E) & np.isfinite(N) & np.isfinite(Z)
        cont = {"E": E[good], "N": N[good], "Z": Z[good], "frame": self.fname_pt[good], "cam": self.cam_pt[good],
                "day": self.day_pt[good]}
        out = np.full(n, np.nan)
        if not good.any():
            return out
        per, _ = waterlines_on_transects(cont, self.survey, TRANSECT_TOLERANCE, TRANSECT_MIN_POINTS,
                                         SURVEY_FLOOR_MARGIN)
        index = {f["name"]: i for i, f in enumerate(self.frames)}
        for fr, dc in zip(per["frame"].tolist(), per["d_cov"].tolist()):
            out[index[fr]] = dc
        return out

    def setups(self, family, coef, station_c):
        """Per frame, the setup of a model: family 'c' (coef = C) or 'stockdon' / 'stockdon_h0hs'
        (coef = k; beta_f from the product's DEM, setup_model.stockdon_frame, C / k where it cannot
        serve) -> (setups, beta_f sources {source: frames} or None)."""
        if family == "c":
            return np.array([float(setup_c(f["hs"], f["tp"], coef)) for f in self.frames]), None
        opt = "none" if family == "stockdon_h0hs" else "auto"
        out, src = [], {}
        for f in self.frames:
            depth = deshoal_depth(opt, f["tag"])[0]
            st, _b, why = stockdon_frame(f["sw"], f["hs"], f["tp"], self.slope, station_c / coef, coef, depth,
                                         why_none="no DEM in the product")
            out.append(st)
            src[why] = src.get(why, 0) + 1
        return np.array(out, float), src

    def score(self, family, coef, station_c):
        """The residuals of the frames scored (self.ref) for one model and coefficient (cached)."""
        key = (family, round(float(coef), 6))
        if key not in self._cache:
            s, src = self.setups(family, coef, station_c)
            r = self.residuals(s)
            self._cache[key] = (r[self.ref] if self.ref is not None else r, src)
        return self._cache[key]

    def set_reference(self, station_c):
        """The frames scored: a finite residual at the station's setup, within GROSS_M of the survey."""
        s, _ = self.setups("c", station_c, station_c)
        r = self.residuals(s)
        ok = np.isfinite(r)
        with np.errstate(invalid="ignore"):
            self.gross = int((ok & (np.abs(r) > GROSS_M)).sum())
            self.ref = np.flatnonzero(ok & (np.abs(r) <= GROSS_M))
        self.ref_cams = np.array([self.frames[i]["cam"] for i in self.ref], object)
        self.reached = int(ok.sum())
        self._cache = {}

    def free_offset(self, r):
        """The residuals of the frames scored less their mean per camera (A FREE OFFSET PER SURVEY AND
        CAMERA, the docstring), NaN for a camera with fewer than MIN_GROUP finite ones."""
        r = np.asarray(r, float)
        out = np.full(r.shape, np.nan)
        for c in sorted(set(self.ref_cams.tolist())):
            m = (self.ref_cams == c) & np.isfinite(r)
            if m.sum() >= MIN_GROUP:
                out[m] = r[m] - r[m].mean()
        return out


# ---------------------------------------------------------------------
# Reading a product
# ---------------------------------------------------------------------

def read_lines(path):
    """The waterline rows a survey product's DEM was built from -> dict of arrays (one per row). pandas, for
    the size of a week of lines (~400k rows), only the columns used."""
    import pandas as pd
    text = ("source_file", "camera", "capture_time_utc", "wave_source")
    nums = ("capture_epoch", "pixel_column", "pixel_row", "tide_elevation_navd88", "offshore_hs_m",
            "offshore_tp_s")
    df = pd.read_csv(path, usecols=lambda c: c in text + nums, dtype=str, keep_default_na=False)
    missing = [k for k in text[:3] + nums if k not in df.columns]
    if missing:
        raise ValueError(f"{path}: no column(s) {', '.join(missing)}")
    out = {k: (df[k].to_numpy(object) if k in df.columns else np.full(len(df), "", object)) for k in text}
    for k in nums:
        out[k] = pd.to_numeric(df[k], errors="coerce").to_numpy(float)
    return out


def load_product(pdir, row, survey_path, calibration_dir, max_points=MAX_POINTS_PER_FRAME):
    """One survey row's product -> (Survey, None) or (None, why not)."""
    from georectify import load_extrinsics, load_intrinsics, pixel_to_ground
    prov_p = pdir / "provenance.json"
    if not prov_p.exists():
        return None, f"no survey-date product here ({prov_p} not found)"
    prov = json.loads(prov_p.read_text())
    if prov.get("blind"):
        return None, "BLIND product (its provenance says so): not read"
    notes = []
    if prov.get("status") not in ("complete", None):
        notes.append(f"product status '{prov.get('status')}'")
    wfile = Path(str((prov.get("filter") or {}).get("file") or ""))
    if not wfile.is_file():
        wdir = pdir / "waterlines"
        f1, f0 = wdir / "contour_points_ground_filtered.csv", wdir / "contour_points_ground.csv"
        wfile = f1 if f1.is_file() else f0
    if not wfile.is_file():
        return None, f"no waterline file in {pdir / 'waterlines'}"
    cams = {}
    for cam, v in sorted((prov.get("cameras") or {}).items()):
        c = v.get("calibration") or {}
        eo, io = Path(str(c.get("eo_path") or "")), Path(calibration_dir) / str(c.get("io_file") or "")
        if eo.is_file() and io.is_file():
            cams[cam] = (load_intrinsics(io), load_extrinsics(eo))
        else:
            notes.append(f"{cam}: calibration not found ({eo if not eo.is_file() else io}): its lines left out")
    if not cams:
        return None, "no camera calibration found"
    dem_stem = ((prov.get("dem") or {}).get("files") or {}).get("dem")
    dem_p = Path(str(dem_stem) + ".asc") if dem_stem else None
    slope = None
    if dem_p is not None and dem_p.is_file():
        slope = ForeshoreSlope.from_file(dem_p)
    else:
        notes.append("no DEM in the product: the Stockdon models use C / k for every frame")
    L = read_lines(wfile)
    n_rows = len(L["source_file"])
    keep = np.ones(n_rows, bool)
    info = (prov.get("dem") or {}).get("built_from")
    if info:
        from survey_compare import dem_frame_selection
        cont = {"E": np.zeros(n_rows), "hs": L["offshore_hs_m"], "frame": L["source_file"], "cam": L["camera"],
                "epoch": L["capture_epoch"]}
        mask, why, _note = dem_frame_selection(cont, info)
        if mask is not None:
            keep &= mask
        if why:
            notes.append("left out as the DEM left them out: " + "; ".join(f"{k} {len(v)}" for k, v in why.items()))
    else:
        notes.append("no DEM _info.json in the provenance: every frame of the window is scored")
    with np.errstate(invalid="ignore"):
        keep &= (np.isfinite(L["tide_elevation_navd88"]) & np.isfinite(L["pixel_column"]) & np.isfinite(L["pixel_row"])
                 & (L["offshore_hs_m"] > 0) & (L["offshore_tp_s"] > 0))
    keep &= np.array([c in cams for c in L["camera"].tolist()], bool)
    kind = "dsm" if row["survey_type"] == "dsm" else "transects"
    if kind == "dsm":
        from compare_dem_survey import read_survey
        grid, x0, y_top, cell = read_survey(survey_path)
        survey = (grid, x0, y_top, cell, DsmFloor(grid, x0, y_top, cell))
    else:
        from compare_rtk import load_transects
        from survey_compare import TRANSECT_SPLIT
        try:
            survey = load_transects(survey_path, TRANSECT_SPLIT)
        except (KeyError, ValueError) as exc:
            return None, f"{survey_path}: no RTK transects read ({exc})"
        if not survey:
            return None, f"{survey_path}: no RTK transects (isolated points are not slope-corrected: not scored)"
    # one entry per frame; its points (decimated on a DSM), projected at two planes
    idx = np.flatnonzero(keep)
    o = idx[np.argsort(L["source_file"][idx], kind="stable")]
    names = L["source_file"][o]
    cuts = np.flatnonzero(names[1:] != names[:-1]) + 1
    frames, fi, E0, N0, dE, dN = [], [], [], [], [], []
    for grp in np.split(o, cuts) if len(o) else []:
        if kind == "dsm" and len(grp) > max_points:
            grp = grp[np.unique(np.linspace(0, len(grp) - 1, max_points).round().astype(int))]
        cam = L["camera"][grp[0]]
        io_, eo_ = cams[cam]
        sw = float(L["tide_elevation_navd88"][grp[0]])
        U, V = L["pixel_column"][grp], L["pixel_row"][grp]
        Ea, Na = pixel_to_ground(U, V, sw, io_, eo_)
        Eb, Nb = pixel_to_ground(U, V, sw + 1.0, io_, eo_)
        ok = np.isfinite(Ea) & np.isfinite(Na) & np.isfinite(Eb) & np.isfinite(Nb)
        if kind == "transects" and ok.any():
            near = np.zeros(len(U), bool)
            for t in survey:
                for Ex, Nx in ((Ea, Na), (Eb, Nb)):
                    re_, rn = Ex - t["c"][0], Nx - t["c"][1]
                    s = re_ * t["d"][0] + rn * t["d"][1]
                    lat = np.abs(-re_ * t["d"][1] + rn * t["d"][0])
                    with np.errstate(invalid="ignore"):
                        near |= (lat <= NEAR_TRANSECT_M) & (s >= t["s"][0] - NEAR_TRANSECT_M) & \
                            (s <= t["s"][-1] + NEAR_TRANSECT_M)
            ok &= near
        if ok.sum() < 3:
            continue
        k = len(frames)
        frames.append({"name": str(L["source_file"][grp[0]]), "cam": str(cam), "sw": sw,
                       "hs": float(L["offshore_hs_m"][grp[0]]), "tp": float(L["offshore_tp_s"][grp[0]]),
                       "tag": str(L["wave_source"][grp[0]] or ""), "day": str(L["capture_time_utc"][grp[0]])[:10]})
        fi.append(np.full(int(ok.sum()), k))
        E0.append(Ea[ok]); N0.append(Na[ok]); dE.append(Eb[ok] - Ea[ok]); dN.append(Nb[ok] - Na[ok])
    if not frames:
        return None, "no frame with waves, a calibration and ground positions"
    if not any(f["tag"] for f in frames):
        notes.append("no wave_source in the waterline file: the deshoaled Stockdon model takes every frame as deep "
                     "water (H0 = Hs), the same as the H0 = Hs one")
    cat = lambda xs: np.concatenate(xs) if xs else np.zeros(0)  # noqa: E731
    sv = Survey(row["name"], row["date"], kind, row.get("label") or "", str(survey_path), frames,
                cat(fi).astype(int), cat(E0), cat(N0), cat(dE), cat(dN), slope, survey, notes)
    sv.waterlines = str(wfile)
    sv.dem = str(dem_p) if dem_p is not None else None
    return sv, None


# ---------------------------------------------------------------------
# Fits and cross-validation
# ---------------------------------------------------------------------

def frame_stats(r, cams=None):
    """n, bias (mean), median, RMSE, NMAD of the per-frame residuals (and the median per camera)."""
    r = np.asarray(r, float)
    ok = np.isfinite(r)
    v = r[ok]
    if not len(v):
        return {"n": 0, "bias": None, "median": None, "rmse": None, "nmad": None}
    med = float(np.median(v))
    out = {"n": int(len(v)), "bias": float(v.mean()), "median": med, "rmse": float(np.sqrt(np.mean(v ** 2))),
           "nmad": float(1.4826 * np.median(np.abs(v - med))) if len(v) >= 3 else None}
    if cams is not None:
        cams = np.asarray(cams, object)
        out["by_camera"] = {c: {"n": int((ok & (cams == c)).sum()),
                                "median": float(np.median(r[ok & (cams == c)]))}
                            for c in sorted(set(cams[ok].tolist()))}
    return out


def mean_square(sv, family, coef, station_c, free=False):
    """The survey's mean square residual at a model (free: with a free offset per camera), or inf when
    it is not admissible there (it keeps fewer than MIN_KEEP of the frames scored)."""
    r, _ = sv.score(family, coef, station_c)
    ok = np.isfinite(r)
    if not len(r) or ok.sum() < MIN_KEEP * len(r):
        return np.inf
    if free:
        r = sv.free_offset(r)
        ok = np.isfinite(r)
        if not ok.any():
            return np.inf
    return float(np.mean(r[ok] ** 2))


def fit(surveys, family, grid, station_c, free=False):
    """The coefficient of the grid that minimises the mean over `surveys` of their mean square residuals
    (each survey weighs the same; free: with a free offset per survey and camera). -> (coefficient or
    None, objective)."""
    if not surveys:
        return None, np.inf
    J = np.array([np.mean([mean_square(sv, family, a, station_c, free) for sv in surveys]) for a in grid])
    if not np.isfinite(J).any():
        return None, np.inf
    i = int(np.argmin(J))
    return float(grid[i]), float(J[i])


def survey_result(sv, family, coef, station_c, free=False):
    r, src = sv.score(family, coef, station_c)
    lost = int(len(sv.ref) - np.isfinite(r).sum())
    out = dict(frame_stats(sv.free_offset(r) if free else r, sv.ref_cams), coef=coef,
               frames_scored=int(len(sv.ref)))
    out["lost"] = lost
    if src is not None:
        out["beta_sources"] = src
    return out


def evaluate(surveys, station_c, c_grid=C_GRID, k_grid=K_GRID):
    """Every model on every survey: fixed ones, refits in-sample, per survey alone and leave-one-survey-out."""
    for sv in surveys:
        sv.set_reference(station_c)
    use = [sv for sv in surveys if len(sv.ref) >= MIN_FRAMES]
    res = {"surveys": [], "fixed": {}, "refit": {}, "station_c": station_c, "used": [sv.name for sv in use]}
    for sv in surveys:
        res["surveys"].append({"name": sv.name, "date": sv.date, "kind": sv.kind, "label": sv.label,
                               "survey": sv.path, "frames": len(sv.frames), "reached": sv.reached,
                               "gross_failures": sv.gross, "scored": int(len(sv.ref)), "used": sv in use,
                               "circular_c": list(sv.circular_c), "notes": sv.notes,
                               "waterlines": getattr(sv, "waterlines", None), "dem": getattr(sv, "dem", None)})
    fixed = [("C = 0 (no setup)", "c", 0.0), (f"C = {station_c:g} (the station's)", "c", station_c),
             (f"Stockdon k = {STOCKDON_SETUP_K:g}, H0 = Hs/Ks (the paper)", "stockdon", STOCKDON_SETUP_K)]
    for title, fam, coef in fixed:
        res["fixed"][title] = {"family": fam, "coef": coef,
                               "per_survey": {sv.name: survey_result(sv, fam, coef, station_c) for sv in use}}
    for fam, grid in (("c", c_grid), ("stockdon", k_grid), ("stockdon_h0hs", k_grid)):
        a_all, J = fit(use, fam, grid, station_c)
        rec = {"family": fam, "what": FAMILIES[fam], "grid": [float(grid[0]), float(grid[-1]),
                                                               float(grid[1] - grid[0])],
               "all": a_all, "all_rmse_mean": float(np.sqrt(J)) if np.isfinite(J) else None,
               "alone": {}, "in_sample": {}, "loso": {}}
        for sv in use:
            rec["alone"][sv.name] = fit([sv], fam, grid, station_c)[0]
            if a_all is not None:
                rec["in_sample"][sv.name] = survey_result(sv, fam, a_all, station_c)
            a_h, _ = fit([o for o in use if o is not sv], fam, grid, station_c)
            rec["loso"][sv.name] = (dict(survey_result(sv, fam, a_h, station_c), coef=a_h) if a_h is not None
                                    else {"coef": None, "n": 0, "rmse": None, "bias": None})
        r = [v["rmse"] for v in rec["loso"].values() if v.get("rmse") is not None]
        rec["loso_rmse_mean"] = float(np.mean(r)) if r and len(r) == len(use) else None
        # the same with a free offset per survey and camera (the docstring)
        a_free, _ = fit(use, fam, grid, station_c, free=True)
        fr = rec["free_offset"] = {"all": a_free, "alone": {}, "loso": {}}
        for sv in use:
            fr["alone"][sv.name] = fit([sv], fam, grid, station_c, free=True)[0]
            a_h, _ = fit([o for o in use if o is not sv], fam, grid, station_c, free=True)
            fr["loso"][sv.name] = (survey_result(sv, fam, a_h, station_c, free=True) if a_h is not None
                                   else {"coef": None, "n": 0, "rmse": None, "nmad": None})
        r = [v["rmse"] for v in fr["loso"].values() if v.get("rmse") is not None]
        fr["loso_rmse_mean"] = float(np.mean(r)) if r and len(r) == len(use) else None
        # a coefficient on the edge of its grid: the optimum may lie beyond it
        edge = (float(grid[0]), float(grid[-1]))
        rec["at_grid_edge"] = [w for w, a in
                               [("all", a_all)] + [(f"alone {n}", a) for n, a in rec["alone"].items()]
                               + [(f"loso {n}", v.get("coef")) for n, v in rec["loso"].items()]
                               + [("free_offset all", a_free)]
                               + [(f"free_offset alone {n}", a) for n, a in fr["alone"].items()]
                               + [(f"free_offset loso {n}", v.get("coef")) for n, v in fr["loso"].items()]
                               if a is not None and a in edge]
        res["refit"][fam] = rec
    res["recommendation"] = recommend(res, use, station_c)
    return res


def head_to_head(use, res, station_c, free):
    """Per held-out survey, the C refit against Stockdon's (H0 = Hs / Ks), each at its leave-one-survey-out
    coefficient (free: those of the free-offset refits, scored so), over the frames BOTH keep there ->
    {survey: {n, c_rmse, stockdon_rmse, stockdon_keeps} or None when a coefficient is missing}."""
    cr, sr = res["refit"]["c"], res["refit"]["stockdon"]
    if free:
        cr, sr = cr["free_offset"], sr["free_offset"]
    out = {}
    for sv in use:
        a_c, a_s = cr["loso"][sv.name].get("coef"), sr["loso"][sv.name].get("coef")
        if a_c is None or a_s is None:
            out[sv.name] = None
            continue
        rc, _ = sv.score("c", a_c, station_c)
        rs, _ = sv.score("stockdon", a_s, station_c)
        keeps = bool(len(rs)) and bool(np.isfinite(rs).sum() >= MIN_KEEP * len(rs))
        both = np.isfinite(rc) & np.isfinite(rs)
        rc, rs = np.where(both, rc, np.nan), np.where(both, rs, np.nan)
        if free:
            rc, rs = sv.free_offset(rc), sv.free_offset(rs)
        m = np.isfinite(rc) & np.isfinite(rs)
        out[sv.name] = {"n": int(m.sum()),
                        "c_rmse": float(np.sqrt(np.mean(rc[m] ** 2))) if m.any() else None,
                        "stockdon_rmse": float(np.sqrt(np.mean(rs[m] ** 2))) if m.any() else None,
                        "stockdon_keeps": keeps}
    return out


def verdict(h2h, n_use, margin):
    """One test of THE RECOMMENDATION over head_to_head's surveys -> (passed, text, numbers)."""
    scored = {n: v for n, v in h2h.items() if v is not None and v["n"]}
    wins = [n for n, v in sorted(scored.items()) if v["stockdon_keeps"] and v["stockdon_rmse"] < v["c_rmse"]]
    weak = sorted(n for n, v in scored.items() if not v["stockdon_keeps"])
    nums = {"per_survey": h2h, "stockdon_better_on": wins, "stockdon_keeps_too_few_on": weak}
    if len(scored) < n_use:
        return False, (f"scored on {len(scored)} of the {n_use} held-out surveys (a coefficient or the frames both "
                       f"models keep are missing on the others)"), nums
    rc = float(np.mean([v["c_rmse"] for v in scored.values()]))
    rs = float(np.mean([v["stockdon_rmse"] for v in scored.values()]))
    gain = rc - rs
    nums.update(c_rmse_mean=rc, stockdon_rmse_mean=rs, gain_m=gain)
    ok = (gain >= margin if margin > 0 else gain > 0) and len(wins) * 2 > n_use
    need = f"at least {margin:g} m" if margin > 0 else "below C's"
    text = (f"RMSE {rs:.3f} m against C's {rc:.3f} m ({gain:+.3f} m; it needs {need}), lower on {len(wins)} of "
            f"{n_use} held-out surveys (it needs a majority); frames compared "
            + ", ".join(f"{n} {v['n']}" for n, v in sorted(scored.items()))
            + (f"; Stockdon keeps fewer than {MIN_KEEP:g} of the frames scored on {', '.join(weak)} (a loss)"
               if weak else ""))
    return ok, text, nums


def recommend(res, use, station_c):
    """THE RECOMMENDATION (docstring) -> dict with 'model' ('c' or 'stockdon'), 'k' and the reasons."""
    out = {"model": "c", "k": None, "why": [], "c_cross_validated": res["refit"]["c"]["all"]}
    if len(use) < MIN_SURVEYS:
        out["why"].append(f"{len(use)} survey(s) with at least {MIN_FRAMES} frames scored; a recommendation needs "
                          f"{MIN_SURVEYS} (each held out in turn): the station keeps C")
        return out
    ok_abs, t_abs, n_abs = verdict(head_to_head(use, res, station_c, False), len(use), MARGIN_M)
    ok_free, t_free, n_free = verdict(head_to_head(use, res, station_c, True), len(use), 0.0)
    out.update(absolute=n_abs, free_offset=n_free, gain_m=n_abs.get("gain_m"),
               stockdon_better_on=n_abs["stockdon_better_on"])
    out["why"].append(f"cross-validated, absolute: Stockdon's {t_abs}" + ("" if ok_abs else ": NOT MET"))
    out["why"].append(f"cross-validated, a free offset per survey and camera (the wave dependence): Stockdon's "
                      f"{t_free}" + ("" if ok_free else ": NOT MET"))
    if ok_abs and ok_free:
        out["model"], out["k"] = "stockdon", res["refit"]["stockdon"]["all"]
    elif ok_abs:
        out["why"].append("Stockdon's absolute gain does not hold with a free offset per survey and camera: it is in "
                          "the surveys' levels (where the pointing, the still-water datum or a camera's bias can put "
                          "it), not in how the lines follow the waves within a survey. The station keeps C")
    else:
        out["why"].append("the station keeps C")
    return out


# ---------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------

def f3(v, signed=True):
    return "   -  " if v is None else (f"{v:+.3f}" if signed else f"{v:.3f}")


def by_camera(r):
    """'  c1 +0.012 (20), c2 -0.031 (8)': the median per camera, when there are two."""
    b = r.get("by_camera") or {}
    return ("  [" + ", ".join(f"{c} {v['median']:+.3f} ({v['n']})" for c, v in sorted(b.items())) + "]"
            if len(b) > 1 else "")


def report(res, root):
    L = ["=" * 100, "WHICH SETUP MODEL: C against Stockdon, on the surveys (evaluate_setup_model.py)", "=" * 100,
         f"products       : {root}",
         "residual       : per frame, median of (line - survey), each line re-projected at still water + its "
         "setup; positive = the line reads HIGH",
         f"frames scored  : the DEM's frames with Hs and Tp, finite at the station's C = {res['station_c']:g} and "
         f"within {GROSS_M:g} m there (MIN_KEEP {MIN_KEEP:g})", ""]
    L.append("SURVEYS")
    for s in res["surveys"]:
        circ = f"; C = {', '.join(f'{c:g}' for c in s['circular_c'])} was FITTED to it (CIRCULAR there)" \
            if s["circular_c"] else ""
        L.append(f"  {s['date']} {s['name']:<22} {s['kind']:<9} frames {s['frames']:4d}, reaching the survey "
                 f"{s['reached']:4d}, gross failures {s['gross_failures']:3d}, scored {s['scored']:4d}"
                 + ("" if s["used"] else f"  (fewer than {MIN_FRAMES}: not used)") + circ)
        for n in s["notes"]:
            L.append(f"      note: {n}")
    names = res["used"]
    if not names:
        L.append("\nNo survey could be scored.")
        return L
    L += ["", "FIXED MODELS (nothing fitted here): per survey  n  bias  median  RMSE  NMAD"]
    for title, rec in res["fixed"].items():
        L.append(f"  {title}")
        for n in names:
            r = rec["per_survey"][n]
            src = r.get("beta_sources")
            L.append(f"      {n:<24} {r['n']:4d}  {f3(r['bias'])}  {f3(r['median'])}  {f3(r['rmse'], False)}  "
                     f"{f3(r['nmad'], False)}" + (f"  (lost {r['lost']}: off the survey)" if r.get("lost") else "")
                     + by_camera(r)
                     + ("   beta_f of its frames: " + ", ".join(f"{k} {v}" for k, v in sorted(src.items()))
                        if src else ""))
    L += ["", "REFITS (one coefficient each; 'EDGE': on the edge of its grid, the optimum may lie beyond)"]
    for fam, rec in res["refit"].items():
        edge = set(rec.get("at_grid_edge") or ())
        E = lambda w: " EDGE" if w in edge else ""  # noqa: E731
        fr = rec["free_offset"]
        L.append(f"  {rec['what']}")
        L.append(f"      fitted on every survey (in-sample): {rec['all']}{E('all')}; each survey alone: "
                 + ", ".join(f"{n} {rec['alone'][n]}{E('alone ' + n)}" for n in names))
        L.append("      LEAVE-ONE-SURVEY-OUT (fitted on the others, scored on it):  coef   n  bias  median  RMSE  NMAD")
        for n in names:
            r = rec["loso"][n]
            L.append(f"      {n:<24} {str(r.get('coef')):>6} {r.get('n', 0):4d}  {f3(r.get('bias'))}  "
                     f"{f3(r.get('median'))}  {f3(r.get('rmse'), False)}  {f3(r.get('nmad'), False)}"
                     + (f"  (lost {r['lost']})" if r.get("lost") else "") + by_camera(r) + E("loso " + n))
        L.append(f"      cross-validated RMSE, mean over the held-out surveys: {f3(rec['loso_rmse_mean'], False)} m")
        L.append(f"      WITH A FREE OFFSET PER SURVEY AND CAMERA: in-sample {fr['all']}{E('free_offset all')}; each "
                 f"survey alone: " + ", ".join(f"{n} {fr['alone'][n]}{E('free_offset alone ' + n)}" for n in names))
        L.append("        held out (fitted on the others, scored on it):              coef   n  RMSE  NMAD")
        for n in names:
            r = fr["loso"][n]
            L.append(f"      {n:<24}                                    {str(r.get('coef')):>6} {r.get('n', 0):4d}  "
                     f"{f3(r.get('rmse'), False)}  {f3(r.get('nmad'), False)}" + E("free_offset loso " + n))
        L.append(f"        cross-validated RMSE, mean over the held-out surveys: {f3(fr['loso_rmse_mean'], False)} m")
    rec = res["recommendation"]
    L += ["", "RECOMMENDATION"]
    if rec["model"] == "stockdon":
        L.append(f"  Switch the station to the slope-aware setup: in waterline_timex_cron.sh set "
                 f"SETUP_MODEL=\"stockdon\" and SETUP_STOCKDON_K=\"{rec['k']:g}\" (k fitted on every survey).")
    else:
        L.append(f"  Keep the station's model (SETUP_MODEL=\"\", C = {res['station_c']:g}).")
    for w in rec["why"]:
        L.append(f"  {w}")
    if rec.get("c_cross_validated") is not None:
        L.append(f"  (C refitted on every survey: {rec['c_cross_validated']:g}. A new C is changed in all three "
                 f"places -- the cron, survey_products.py, historical_forcing.py -- with a SETUP_FITS entry: "
                 f"docs/VALIDATION.md 4.1. Not done here.)")
    return L


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main(argv=None):
    import survey_products as sp
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--output-root", default=str(sp.DEFAULT_OUTPUT_ROOT),
                    help=f"the survey-date products (survey_products.py --output-root; default "
                         f"{sp.DEFAULT_OUTPUT_ROOT})")
    ap.add_argument("--surveys", default=str(sp.DEFAULT_SURVEYS), help="surveys/surveys.csv")
    ap.add_argument("--survey-dirs", nargs="+", default=list(sp.DEFAULT_SURVEY_DIRS),
                    help="where the survey files are looked for by name (survey_products.py's)")
    ap.add_argument("--calibration", default=str(HERE / "calibration"), help="the lens (IO) files")
    ap.add_argument("--dates", nargs="+", default=None, help="only these product dates")
    ap.add_argument("--station-c", type=float, default=sp.SETUP_COEF,
                    help=f"the station's C (default survey_products.SETUP_COEF, {sp.SETUP_COEF})")
    ap.add_argument("--json", default=None, help="write every number here")
    args = ap.parse_args(argv)
    sp.arm_blind_guard()       # no file named like a held survey (BLIND_DATES) is opened by this run
    rows = sp.read_table(args.surveys, sp.SURVEY_FIELDS)
    root = Path(args.output_root)
    surveys = []
    print(f"surveys table: {args.surveys} ({len(rows)} rows); products: {root}")
    for row in rows:
        tag = f"{row['date']} {row['name']}"
        if args.dates and row["date"] not in args.dates:
            continue
        # a BLIND date's survey is never resolved and its product never read (survey_products.BLIND_DATES)
        if sp.held_survey(row):
            print(f"  {tag}: BLIND (survey_products.BLIND_DATES): skipped, nothing of it opened")
            continue
        if row["survey_type"] not in ("dsm", "points"):
            print(f"  {tag}: survey type {row['survey_type']!r}: not scored")
            continue
        if not row.get("path"):
            print(f"  {tag}: no path in {Path(args.surveys).name}: skipped (never looked up by its date)")
            continue
        if sp.held_name(Path(row["path"]).name):
            print(f"  {tag}: {Path(row['path']).name} is named like a survey of a BLIND date "
                  f"(survey_products.BLIND_DATES): skipped, never opened")
            continue
        hits = sp.resolve_file(row["path"], args.survey_dirs)
        if not hits:
            print(f"  {tag}: {row['path']} not found in {', '.join(args.survey_dirs)}: skipped")
            continue
        sv, why = load_product(root / row["date"], row, hits[0], args.calibration)
        if sv is None:
            print(f"  {tag}: {why}: skipped")
            continue
        sha = sp.sha256(hits[0])
        sv.circular_c = tuple(c for c, f in sorted(sp.SETUP_FITS.items())
                              if f.get("survey_sha256") == sha or Path(str(hits[0])).stem == f.get("survey_stem"))
        print(f"  {tag}: {len(sv.frames)} frame(s) from {sv.waterlines}")
        surveys.append(sv)
    if not surveys:
        print("No survey could be scored (no product, no survey file, or no frame).")
        return 4
    res = evaluate(surveys, args.station_c)
    for line in report(res, root):
        print(line)
    if args.json:
        Path(args.json).write_text(json.dumps(res, indent=1, default=str) + "\n")
        print(f"\nwrote {args.json}")
    return 0 if res["used"] else 4


if __name__ == "__main__":
    sys.exit(main())

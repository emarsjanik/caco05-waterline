#!/usr/bin/env python3
"""
Measured Wave Setup Per Collection, From The Runup Timestacks (Opt-In)
=======================================================================
The DEM gives every timex waterline the still-water level plus a setup
MODELLED from the offshore waves, C * sqrt(Hs * L0) (setup_model.py). This
module MEASURES the setup of each 10-minute collection instead, from the c2
runup timestack of the same collection: the mean elevation of the swash
edge above still water, as timestack_wave_height.py reads it, one value
per collection. OPT-IN: it writes only archive/setup_measured_c2.csv; the
contour extractor uses it only with --setup-source measured or auto (its
default, model, writes the station's files bit for bit as before), and
whether the station ever switches is decided by docs/VALIDATION.md 4.12.

THE RULE
 1. MEASUREMENT. For each stack and line (1 and 2, the c2 lines that cross
    the swash), timestack_wave_height.estimate_stack as it is: the swash
    edge of runup_from_timestack.py, each pixel's ray followed down to the
    beach profile, less still water (marconi_water_level.py, GPS first);
    setup_m is the mean of that series (Stockdon et al. 2006's <eta>). Three
    of its settings differ from that script's own run, none of them in the
    setup's arithmetic:
      (a) the profile is the DEM of rule 2, never a plane (no --beta);
      (b) its edge-contrast label is off: on the real stacks of Oct 2026
          edge_contrast had a median of 0.00-0.02 on EVERY day, clear days
          included, so its provisional 0.05 labelled every stack 'low
          contrast' (it cannot tell blur from a clear day there). It is kept
          as a diagnostic (line_edge_contrast), not a gate; the window's
          sharpness (rule 3) is the blur check. The same finding says the
          real edge has no foam/sand contrast to check, so whether a passed
          edge is a real swash edge is for the audit (--audit: the lines'
          disagreement on clear days and a random sample drawn with the edge
          on), VALIDATION 4.12 criterion 4;
      (c) the period is a placeholder (PLACEHOLDER_T0_S): the setup is a
          mean, with no period in it; the period gates only that script's
          inversion for H0, which is neither used nor written here. So a
          collection has a measured setup in a buoy outage too, where the
          model has none (finding H2 of the review register: no Tp, no
          setup).
 2. THE DEM: STRICTLY EARLIER AND MODEL-BUILT (the non-circularity rule).
    The swash elevation is read off a DEM, and a DEM of timex waterlines
    carries the setup its lines were given. That makes two loops, and the
    rule breaks both:
      (a) STRICTLY EARLIER. On the DEM of its own week the measured setup
          tends to the model's C sqrt(Hs L0) itself (timestack_wave_height
          .py, CIRCULARITY): it echoes the model. So the profile comes only
          from a DEM built from days strictly before the stack's UTC day: of
          the dated DEMs (dem_<window end>_<days>d_dem.asc, dem_from_contours
          .py --series-dir; the window is its last <days> UTC days), the one
          whose window ENDS latest while still before that day, and at most
          --max-dem-age-days (MAX_DEM_AGE_DAYS) before it. A DEM whose window
          holds the stack's day is never used -- the DEM being built from the
          frame's own window included.
      (b) MODEL-BUILT. A DEM built from lines that were given MEASURED setups
          passes its own error on with UNIT GAIN. A DEM off the beach by e
          reads the setup off by e / (1 - beta / tan(delta)) (the edge is put
          on the wrong surface: the rays meet it early or late; beta the
          beach's slope along the ray, delta the ray's depression), and the
          waterline of the same collection, given that setup and re-projected
          on its own ray, lands off by e again: the ray factor and its
          inverse, whatever their value ([REAL] VALIDATION 4.1 measures the
          factor at ~0.71 per frame on the 2026 RTK). On top come the
          constant offsets the measurement carries and C absorbs (the
          detector's bias; the gap between the timex line and the mean swash
          edge). Chained through the cron's 7-day windows nothing pulls such
          a DEM back to the beach: every generation adds those offsets again,
          random errors add up as a random walk, and beach change at the
          lines is erased [SYNTH: test_measured_setup.py, test_chain, the
          real c2 calibration, a 1:10 beach, H0 1.5 m, T 10 s: one generation
          0.049 + 0.996 e; chained, the DEM at the c2 line +0.08 / +0.16 /
          +0.35 / +0.70 m off after 7 / 14 / 30 / 60 days, and a 0.30 m
          erosion at day 20 read as +0.24 m of change between days 15 and
          35]. So a DEM is used only when its dem_<end>_<days>d_info.json
          (dem_from_contours.py writes it beside every series DEM) says it
          was built with the model's setup: no setup_source count other than
          MODEL_BUILT_SOURCES ('model', 'none'; extract_elevation_contours.py
          --setup-source labels every frame, and the default rows carry no
          such column), and a setup on its lines at all (setup_share > 0: a
          DEM of lines at still water would put every setup low by the whole
          setup x ~1.5). A DEM with no info file, or one that cannot be read,
          is refused too: its setup is not known. With SETUP_SOURCE on, the
          cron keeps a model-only series beside its own (archive/dems_model,
          a second run of the same chain with the model's setup), and --dem's
          default is that folder when it exists, else archive/dems (the
          cron's own series, model-built while SETUP_SOURCE is off). Every
          measured setup is then ONE generation from the model: the offsets
          stay a constant (+0.049 m on that scene) and the change is kept.
    No DEM left: no measurement ('no slope', the reason in the note). Each
    row records the DEM, its window, its age (days from the window's end to
    the stack's day), its sha256 (the first 12 hex) and its setup_share, so
    the rule can be checked afterwards; the run lists the DEMs it refused.
    WHAT THAT BUYS, AND WHAT IT DOES NOT. The measured setup corrects the
    model's DAY-TO-DAY departures, those the offshore Hs and Tp cannot see:
    a buoy outage, the wave direction (sheltered south-west seas),
    groupiness and infragravity energy. It INHERITS the level of the
    model-built DEM it stands on at the c2 lines: that DEM's error there
    (the model's setup error averaged over its window, its own noise, and
    beach change since its window) reads as setup error x 1 / (1 - beta /
    tan(delta)) [SYNTH: a DEM 0.10 m above the beach read the setup 0.13 /
    0.15 / 0.19 m higher on a 1:14 / 1:10 / 1:8 beach, i.e. a factor ~1.5;
    the 1:10 beach is the test's, 1:14 and 1:8 a run made while writing
    it]. Given back to the waterline of the same collection, which lies
    near the mean swash position, that setup puts the line onto the DEM it
    stood on: AT THE TIMESTACK LINES a DEM built with measured setups
    follows the model-built DEM of the day before, offset by a constant,
    and keeps its change, lagged by the DEM's age; elsewhere in the view
    (c1, the rest of c2) the lines get the c2 lines' setup, its error
    included. It cannot correct C's long-term level: the survey-date check
    (evaluate_measured_setup.py) can.
 3. QC, per line, all of:
      status ok        timestack_wave_height's own: not dark, the window not
                       blurred, an edge in >= 80% of the rows with >= 20
                       uprushes, a profile, an estimate;
      window_sharpness >= 10 (station_status.LENS_BLURRED; owg_live.py's
                       sharpness of the c2 bright frame of the same
                       collection: clear days read 49-88, the blurred window
                       of 20-24 Oct 2026 3.4-4.7). A stack with no such frame
                       fails: its window is unknown (--allow-no-sharpness for
                       a period without the OWG record; evaluate_measured_
                       setup.py --sharpness-images scores that period's own
                       images instead);
      valid edge       valid_fraction >= 0.8 and n_peaks >= 20 (the runup
                       script's own test, checked again here);
      off_profile_fraction <= 0.02 (timestack_wave_height.py's own limit
                       for a setup: samples off the profile are not in the
                       mean);
      z_extrapolated_fraction <= MAX_EXTRAPOLATED (0.10, PROVISIONAL): the
                       share of the edge read on the profile's straight
                       extension past the DEM's end (at most 10 m). A 1:10
                       foreshore that steepens or flattens by 0.03 over that
                       10 m is off by at most 0.15 m at its far end, ~0.07 m
                       on average, so 10% of samples there move the mean by
                       ~7 mm; the NUC checks the share on real stacks
                       (VALIDATION 4.12);
      setup_m          finite;
      still water      from an allowed source (--level-sources; default
                       'gps', the GNSS-R level the contour extractor gives the
                       frames: a setup measured above one level and added to
                       another carries their difference 1:1).
 4. THE COLLECTION: the median of setup_m over the lines that pass (for two
    lines, their mean); status ok when at least --min-lines (1) pass.
 5. ITS UNCERTAINTY. setup_sigma_m = sqrt(sampling^2 + lines^2) is the
    collection's OWN RANDOM error as an estimate of the mean setup of its 10
    minutes, between collections standing on the same profile:
      sampling  the standard error of a 10-minute mean, from the series'
                own low-frequency level: G0, the mean one-sided periodogram
                over 0 < f <= 0.005 Hz (the three lowest bins of 10 minutes;
                gaps bridged linearly), taken as the spectral density at
                zero frequency, SE = sqrt(G0 / (2 T)) (mean_sampling_se). On
                Gaussian records of known spectrum [SYNTH: test_measured_
                setup.py and a larger run made while writing it] its rms is
                1.01-1.10 x the true sd for a spectrum flat to zero
                frequency, 2.5-2.9 x for spectra empty below 0.002-0.005 Hz
                (1/f or flat above), and 0.94-0.96 x for one that keeps its
                level below 1/T (flat below 0.003 Hz, an incident peak
                beside): mostly an upper bound, at most ~6% low on those,
                with ~30% scatter of its own per record. What lies below 1/T
                a single 10-minute record cannot see. The lines' sampling
                errors are taken as FULLY correlated (the infragravity swash
                is coherent over the tens of metres between the lines):
                their mean, not reduced by averaging the lines;
      lines     the standard deviation of the lines' setups / sqrt(n) (for
                two lines half their difference): their alongshore
                difference and their independent errors. With one line it is
                unknown and setup_sigma_m is the sampling term alone (noted).
    It is NOT the error of the setup a waterline is given, and it is the
    smaller part of that. The timex of the same collection sees the same 10
    minutes of swash, so for that frame the sampling term is no error at
    all; the larger terms are elsewhere:
      sigma_profile_m  the profile's term: the DEM's own standard error along
                the swash (the median over the line's pixels whose profile
                lies within the swash's 5-95% elevations, of the dem_<...>_
                sigma_obs.asc beside the DEM where there is one, else its
                _sigma.asc, the cells' repeatability: dem_from_contours.py),
                times profile_gain = 1 / (1 - beta / tan(delta)) at the mean
                swash position (beta the DEM's slope along the ray towards
                the camera over one cell either side, delta the ray's
                depression; ~1.5 on the synthetic 1:10 scene, NaN where the
                ray grazes the beach, beta / tan(delta) >= MAX_RAY_RATIO).
                Common to every collection standing on that DEM (a bias
                between them, not noise), so not in setup_sigma_m; and short
                of the DEM's level error, which no grid holds (the model's
                setup in its lines) [REAL: on the station DEM of 27 Sep - 2
                Oct 2026 --sigma-obs gave a median of 0.258 m per cell,
                VALIDATION 4.10; times ~1.5]. A survey's surface has no grid:
                blank;
      the detector's bias, in no column: on the synthetic scenes the edge
                trails a fast tip, setup +0.04 to +0.10 m (timestack_wave_
                height.py and test_measured_setup.py [SYNTH]).
 6. WAITING. A stack younger than --wait-days (3) whose DEM or water level
    is missing, or whose still water is not the GPS where the GPS is an
    allowed source (--level-sources holds 'gps', the default), is left for a
    later run, as timestack_wave_height.py does: the DEM of the day before a
    stack is built by a later cron run, the GPS level arrives 1-2 days late
    (and a run whose GNSS-R spline cannot be read must not write young
    stacks off as 'still water from chatham'). Older, it is written as it
    is (failed qc, the reason in its note and line_qc), and only
    --reprocess measures it again.

HOW THE EXTRACTOR USES IT (extract_elevation_contours.py --setup-source
measured | auto --setup-measured FILE). A frame whose file-name epoch lies
within MATCH_TOL_S (300 s) of a collection's start gets that collection's
setup_m, when its status is ok (timex, snap and stack of one collection
share the file-name epoch, timex_timing.py; 300 s is half the collection):
the frames of both cameras, one value per collection, as the model gives
one per frame. 'measured': a frame without one gets no setup (blank, as a
frame without Tp); 'auto': the model's. Its setup_source column says which.

THE AUDIT (--audit N; VALIDATION 4.12 criteria 3 and 4). Reads the output
file: the distribution of |line 1 - line 2| over the collections where both
passed, by the window's sharpness (clear, >= CLEAR_SHARPNESS 25, owg_live
.py's own flag; low, 10-25), with the edge contrast beside it; the QC's
failure reasons; every DEM a passed collection stood on, checked again
against rule 2b (its _info.json in the --dem folder now: the count not
model-built must be 0); and a random sample (--audit-seed) of N passed
collections measured again and drawn with the edge on
(runup_from_timestack.save_plot: the foam anomaly, and the edge in red) in
--audit-dir, with an index.csv. It writes nothing else.

WHAT IT CANNOT DO
  * The inheritance of rule 2: the level of the model-built DEM at the c2
    lines, and beach change since it, are read as setup.
  * Still water is not setup-free: the GNSS-R footprint is the surf zone,
    whose level already holds part of the setup (timestack_wave_height.py,
    2.). The measured setup is above that same level, which is what a line
    given that level needs; it is not Stockdon's setup above the tide.
  * The timex waterline is where the detector finds the edge of a
    time-averaged image, not the mean elevation of the instantaneous edge.
    How the two relate (the share of the swash a line holds) is what the
    survey-date check measures (evaluate_measured_setup.py); C, fitted to
    lines against a survey, absorbs it, the measurement does not.
  * One place: c2's lines 1 and 2. The c1 view and the rest of c2's beach
    get the same number (alongshore differences are not followed: the model
    gives one number per frame too).
  * Night stacks are dark, fog and a wet window fail the sharpness check:
    those collections fall back to the model (auto) or have no setup
    (measured). Nothing checks that a passed clear-day edge is the real
    swash edge (rule 1b): the audit is a look, not a test.
  * Evidence here is synthetic [SYNTH]: on the synthetic scene the setup is
    read +0.04 to +0.10 m high (the detector), on the real stacks of 17-18
    Oct 2026 timestack_wave_height.py read 0.16-0.55 m against Stockdon's
    ~0.31 m for H0 ~1 m, T 7 s, beta 0.10 (on the nearest DEM, i.e. the
    circular one; NUC, not reproduced here) [REAL].

OUTPUT archive/setup_measured_<cam>.csv, one row per collection, APPENDED
(every CHUNK_ROWS collections): a collection already there is never
rewritten except by --reprocess, and a file of other columns is refused
before any stack is read.
  collection_epoch (the collection's start, the key the frames carry),
  time_utc, filename, status (ok / failed qc), note, setup_m, setup_sigma_m
  (rule 5: the collection's own random error), sigma_sampling_m,
  sigma_lines_m, sigma_profile_m and profile_gain (rule 5: the profile's
  term; the detector's bias is in no column), n_lines (lines that passed),
  lines_used, line_setups ('1:0.312;2:0.295', every line with a value),
  line_qc ('1:ok;2:window sharpness 4.2 < 10'), line_edge_contrast (rule
  1b, a diagnostic), window_sharpness, water_level_m, water_level_source,
  dem, dem_window (first..last day), dem_age_days, dem_sha256,
  dem_setup_share (its info.json's setup_share)

Usage:
    python3 measured_setup.py                       new stacks of archive/ras_c2
    python3 measured_setup.py --since 2026-10-01    from that UTC day
    python3 measured_setup.py --reprocess           every stack again
    python3 measured_setup.py --audit 30            the real-data checks of criterion 4
Tests: python3 test_measured_setup.py
"""

import sys
import csv
import json
import argparse
import hashlib
from pathlib import Path
from datetime import datetime, timezone, timedelta

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import timestack_wave_height as thw                                          # noqa: E402
from runup_from_timestack import SAMPLE_HZ, BURST_MID_OFFSET_S, epoch_from_name  # noqa: E402
from setup_model import SERIES_NAME                                          # noqa: E402
from station_status import LENS_BLURRED                                      # noqa: E402
from timex_timing import TIMEX_MID_OFFSET_S                                  # noqa: E402

FIELDS = ["collection_epoch", "time_utc", "filename", "status", "note", "setup_m", "setup_sigma_m",
          "sigma_sampling_m", "sigma_lines_m", "sigma_profile_m", "profile_gain", "n_lines", "lines_used",
          "line_setups", "line_qc", "line_edge_contrast", "window_sharpness", "water_level_m", "water_level_source",
          "dem", "dem_window", "dem_age_days", "dem_sha256", "dem_setup_share"]
STATUSES = ("ok", "failed qc")

DEFAULT_LINES = thw.DEFAULT_LINES     # 1 and 2: the c2 lines that cross the swash
MAX_DEM_AGE_DAYS = 3                  # the DEM's window ends 1-3 days before the stack's day
MODEL_BUILT_SOURCES = ("model", "none")   # the setup_source counts a model-built DEM may hold (rule 2b)
MIN_SHARPNESS = LENS_BLURRED          # 10: owg_live.py's 'blurred' (station_status.py)
CLEAR_SHARPNESS = 25.0                # owg_live.py --flag-sharpness: sharper is a clear window (the audit)
MAX_OFF_PROFILE = thw.MISSING_ELEV_MAX    # 0.02: timestack_wave_height's own limit for a setup
MAX_EXTRAPOLATED = 0.10               # share of the edge on the profile's extension: PROVISIONAL (QC)
MIN_VALID, MIN_PEAKS = 0.8, 20        # the runup script's own test (runup_from_timestack.py defaults)
LEVEL_SOURCES = ("gps",)              # the still water the contour extractor gives the frames
MIN_LINES = 1
MATCH_TOL_S = float(TIMEX_MID_OFFSET_S)   # 300 s: a frame and a stack of one collection (half of it)
SAMPLING_F_LOW_HZ = 0.005             # the lowest bins taken as the spectral level at zero frequency
MAX_RAY_RATIO = 0.9                   # beta / tan(delta) above this: the ray grazes the beach, no gain
PLACEHOLDER_T0_S = 10.0               # gates timestack_wave_height's inversion only (rule 1c)
WAIT_DAYS = thw.WAIT_DAYS
CHUNK_ROWS = 25                       # collections per append (a crash loses at most these)
AUDIT_SAMPLE = 30

_num, _r, _day_epoch = thw._num, thw._r, thw._day_epoch


# ---------------------------------------------------------------------
# The profile: the strictly-earlier, model-built DEM (or one fixed surface)
# ---------------------------------------------------------------------

def utc_day(epoch):
    return datetime.fromtimestamp(float(epoch), tz=timezone.utc).date()


def file_sha256(path, n=12):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:n]


def default_dem(base=HERE):
    """--dem's default (rule 2b): archive/dems_model, the model-only series the cron keeps while SETUP_SOURCE is
    on, when that folder exists; else archive/dems, the cron's own series (model-built while it is off)."""
    model = Path(base) / "archive" / "dems_model"
    return model if model.is_dir() else Path(base) / "archive" / "dems"


def dem_provenance(dem_path):
    """
    Was a dated DEM built with the model's setup (rule 2b)? Read from its dem_<end>_<days>d_info.json, which
    dem_from_contours.py writes beside every series DEM. -> (ok, why not, setup_share).
    """
    dem_path = Path(dem_path)
    stem = dem_path.name[:-len("_dem.asc")] if dem_path.name.endswith("_dem.asc") else dem_path.stem
    info_p = dem_path.with_name(stem + "_info.json")
    try:
        info = json.loads(info_p.read_text())
    except FileNotFoundError:
        return False, "no _info.json (its setup is not known)", np.nan
    except (OSError, ValueError) as exc:
        return False, f"_info.json unreadable ({exc})", np.nan
    if not isinstance(info, dict):
        return False, "_info.json is not a record", np.nan
    share = _num(info.get("setup_share"))
    src = info.get("setup_source")
    if src is not None:
        if not isinstance(src, dict):
            return False, "its setup_source is not a count of frames", share
        bad = {str(k): v for k, v in src.items() if str(k) not in MODEL_BUILT_SOURCES and _num(v) != 0}
        if bad:
            return False, "built with " + ", ".join(f"{v} {k}" for k, v in sorted(bad.items())) + " frame(s)", share
    if not share > 0:
        return False, f"no setup on its lines (setup_share {info.get('setup_share')})", share
    return True, "", share


def sigma_beside(dem_path):
    """The DEM's own standard-error grid (rule 5): dem_<...>_sigma_obs.asc where there is one, else _sigma.asc;
    None without either."""
    p = Path(dem_path)
    stem = p.name[:-len("_dem.asc")] if p.name.endswith("_dem.asc") else p.stem
    for suf in ("_sigma_obs.asc", "_sigma.asc"):
        q = p.with_name(stem + suf)
        if q.is_file():
            return q
    return None


class PriorDem:
    """
    The DEM a stack's swash is put on (rule 2): of the dated DEMs of `source` (a folder of
    dem_<window end>_<days>d_dem.asc, or one such file), the one whose window ends latest while still
    BEFORE the stack's UTC day, at most max_age_days before it, and built with the model's setup
    (dem_provenance). for_epoch(epoch) is what timestack_wave_height.estimate_stack calls (-> its Dem or
    None); info(epoch) what a row records; sigma_for(epoch) its standard-error grid (rule 5).
    """

    def __init__(self, source, max_age_days=MAX_DEM_AGE_DAYS):
        self.max_age = int(max_age_days)
        self.series = []                                   # [(window end date, days, path)]
        p = Path(str(source))
        files = sorted(p.iterdir()) if p.is_dir() else ([p] if p.is_file() else [])
        for f in files:
            m = SERIES_NAME.match(f.name)
            if m:
                end = datetime.strptime(m.group(1), "%Y-%m-%d").date()
                self.series.append((end, int(m.group(2)), f))
        self.refused = p.is_file() and not self.series
        self.what = (f"{len(self.series)} dated DEM(s) in {p}" if p.is_dir() else
                     (f"{p} (not a dated DEM: its window is not known, refused)" if self.refused else
                      (str(p) if self.series else f"{p} (NOT FOUND)")))
        self._dems, self._sha, self._prov, self._sig = {}, {}, {}, {}
        self.not_model_built = {}                          # name -> why, every DEM refused by rule 2b

    def provenance(self, f):
        if f not in self._prov:
            self._prov[f] = dem_provenance(f)
            if not self._prov[f][0]:
                self.not_model_built[f.name] = self._prov[f][1]
        return self._prov[f]

    def pick(self, epoch):
        """-> (path, {dem, dem_window, dem_age_days, dem_sha256, dem_setup_share}) or (None, why none)."""
        day = utc_day(epoch)
        cands = [(end, days, f) for end, days, f in self.series if end < day and (day - end).days <= self.max_age]
        if not cands:
            return None, f"no DEM whose window ends 1-{self.max_age} day(s) before {day}"
        refused = []
        # latest end; then the shorter window; then by name
        for end, days, f in sorted(cands, key=lambda c: (-c[0].toordinal(), c[1], c[2].name)):
            ok, why, share = self.provenance(f)
            if not ok:
                refused.append(f"{f.name}: {why}")
                continue
            first = end - timedelta(days=days - 1)
            if f not in self._sha:
                self._sha[f] = file_sha256(f)
            return f, {"dem": f.name, "dem_window": f"{first}..{end}", "dem_age_days": (day - end).days,
                       "dem_sha256": self._sha[f], "dem_setup_share": _r(share, 4)}
        return None, (f"no model-built DEM whose window ends 1-{self.max_age} day(s) before {day} ("
                      + "; ".join(refused) + ")")

    def for_epoch(self, epoch):
        f, _info = self.pick(epoch)
        if f is None:
            return None
        if f not in self._dems:
            try:
                self._dems[f] = thw.Dem(f)
            except (OSError, ValueError, KeyError) as exc:
                print(f"WARNING: cannot read the DEM {f} ({exc})")
                self._dems[f] = None
        return self._dems[f]

    def sigma_for(self, epoch):
        f, _info = self.pick(epoch)
        if f is None:
            return None
        if f not in self._sig:
            q = sigma_beside(f)
            try:
                self._sig[f] = thw.Dem(q) if q is not None else None
            except (OSError, ValueError, KeyError):
                self._sig[f] = None
        return self._sig[f]

    def info(self, epoch):
        f, info = self.pick(epoch)
        if f is None:
            return {"dem": "", "dem_window": "", "dem_age_days": "", "dem_sha256": "", "dem_setup_share": "",
                    "why": info}
        return info


class FixedProfile:
    """One surface for every stack (a survey: evaluate_measured_setup.py), and what a row records of it."""

    def __init__(self, dem, label, sha=""):
        self.dem, self.label, self.sha = dem, label, sha
        self.what = f"{dem.name} ({label})"

    def for_epoch(self, epoch):
        return self.dem

    def sigma_for(self, epoch):
        return None

    def info(self, epoch):
        return {"dem": self.dem.name, "dem_window": self.label, "dem_age_days": "", "dem_sha256": self.sha,
                "dem_setup_share": ""}


class NoPeriod:
    """A placeholder period for every stack (rule 1c): estimate_stack's inversion needs one, the setup does
    not; never pending (the camera's period is not waited for)."""

    def for_stack(self, filename, mid_epoch):
        return PLACEHOLDER_T0_S, "placeholder (not used: the setup has no period in it)", False


# ---------------------------------------------------------------------
# Per line: QC, the sampling error and the profile term; per collection: the median
# ---------------------------------------------------------------------

def mean_sampling_se(series, fs=SAMPLE_HZ, f_low=SAMPLING_F_LOW_HZ):
    """
    The standard error of the mean of a record (rule 5): sqrt(G0 / (2 T)), G0 the mean one-sided
    periodogram (2 |X_k|^2 / (fs N)) of the demeaned record over 0 < f <= f_low, T its length. Gaps are
    bridged linearly (the lowest bins barely move). NaN with fewer than 64 samples or under half of them.
    """
    x = np.asarray(series, float)
    ok = np.isfinite(x)
    if len(x) < 64 or ok.sum() < max(64, 0.5 * len(x)):
        return np.nan
    if not ok.all():
        i = np.arange(len(x))
        x = np.interp(i, i[ok], x[ok])
    x = x - x.mean()
    n = len(x)
    f = np.fft.rfftfreq(n, 1.0 / fs)
    p = 2.0 * np.abs(np.fft.rfft(x)) ** 2 / (fs * n)
    sel = (f > 0) & (f <= f_low)
    if not sel.any():
        return np.nan
    return float(np.sqrt(p[sel].mean() / (2.0 * n / fs)))


def _nearest_cell(grid_dem, E, N):
    """The value of the grid cell holding each point (NaN outside or on nodata)."""
    r, c = grid_dem._rc(E, N)
    return grid_dem._cell(grid_dem.z, np.asarray(r, float) + 0.5, np.asarray(c, float) + 0.5)


def profile_term(cap, io, eo, sigma_grid=None, max_ratio=MAX_RAY_RATIO):
    """
    Rule 5's profile term of one line, from what estimate_stack captured (its profile, the line's pixels,
    the DEM and the swash's elevations): -> (gain, sigma). gain = 1 / (1 - beta / tan(delta)) at the mean
    swash position (the line's profile first reaching the swash's mean elevation from the sea end; its
    pixel's ray: tan(delta) = 1 / its horizontal metres per metre of height, georectify.pixel_to_ground at
    0 and 1 m; beta the DEM's rise along that ray towards the camera, one cell either side); NaN where the
    ray grazes the beach (beta / tan(delta) >= max_ratio) or the DEM is not there. sigma: the median of
    sigma_grid over the line's pixels whose profile lies within the swash's 5-95% elevations (NaN without
    a grid).
    """
    import georectify
    prof, dem = cap.get("prof") or {}, cap.get("dem")
    z_abs = np.asarray(cap.get("z_abs"), float)
    zs = z_abs[np.isfinite(z_abs)]
    z = np.asarray(prof.get("z", []), float)
    uv = np.asarray(cap.get("uv"), float)
    if dem is None or len(zs) < 10 or np.isfinite(z).sum() < 2 or len(uv) != len(z):
        return np.nan, np.nan
    zbar = float(zs.mean())
    with np.errstate(invalid="ignore"):
        hit = np.isfinite(z[:-1]) & np.isfinite(z[1:]) & ((z[:-1] - zbar) * (z[1:] - zbar) <= 0) & (z[:-1] != z[1:])
    gain = np.nan
    idx = np.flatnonzero(hit)
    if len(idx):
        i = int(idx[0])
        fr = (zbar - z[i]) / (z[i + 1] - z[i])
        u = np.array([uv[i, 0] + fr * (uv[i + 1, 0] - uv[i, 0])])
        v = np.array([uv[i, 1] + fr * (uv[i + 1, 1] - uv[i, 1])])
        E0, N0 = georectify.pixel_to_ground(u, v, 0.0, io, eo)
        E1, N1 = georectify.pixel_to_ground(u, v, 1.0, io, eo)
        dE, dN = float(np.asarray(E1)[0] - np.asarray(E0)[0]), float(np.asarray(N1)[0] - np.asarray(N0)[0])
        h = float(np.hypot(dE, dN))                    # horizontal m per m of height: 1 / tan(delta)
        if np.isfinite(h) and h > 0:
            E, N = float(np.asarray(E0)[0]) + zbar * dE, float(np.asarray(N0)[0]) + zbar * dN
            ux, uy, st = dE / h, dN / h, float(dem.cell)
            zp = dem.sample(np.array([E + st * ux]), np.array([N + st * uy]))[0]
            zm = dem.sample(np.array([E - st * ux]), np.array([N - st * uy]))[0]
            ratio = (zp - zm) / (2.0 * st) * h            # beta / tan(delta)
            if np.isfinite(ratio) and ratio < max_ratio:
                gain = float(1.0 / (1.0 - ratio))
    sigma = np.nan
    if sigma_grid is not None:
        lo, hi = np.percentile(zs, [5, 95])
        with np.errstate(invalid="ignore"):
            band = np.isfinite(z) & (z >= lo) & (z <= hi)
        if band.any():
            vals = _nearest_cell(sigma_grid, np.asarray(prof["E"], float)[band], np.asarray(prof["N"], float)[band])
            vals = vals[np.isfinite(vals)]
            if len(vals):
                sigma = float(np.median(vals))
    return gain, sigma


class Qc:
    """The per-line limits of rule 3."""

    def __init__(self, min_sharpness=MIN_SHARPNESS, require_sharpness=True, max_off_profile=MAX_OFF_PROFILE,
                 max_extrapolated=MAX_EXTRAPOLATED, min_valid=MIN_VALID, min_peaks=MIN_PEAKS,
                 level_sources=LEVEL_SOURCES):
        self.min_sharpness, self.require_sharpness = float(min_sharpness), bool(require_sharpness)
        self.max_off_profile, self.max_extrapolated = float(max_off_profile), float(max_extrapolated)
        self.min_valid, self.min_peaks = float(min_valid), int(min_peaks)
        self.level_sources = tuple(level_sources)

    def describe(self):
        return (f"status ok, window sharpness >= {self.min_sharpness:g}"
                + ("" if self.require_sharpness else " where known (--allow-no-sharpness)")
                + f", valid >= {self.min_valid:g} with >= {self.min_peaks} uprushes, off profile <= "
                f"{self.max_off_profile:g}, extrapolated <= {self.max_extrapolated:g}, still water from "
                + "/".join(self.level_sources))


def qc_line(row, qc=None):
    """Why one line's row (timestack_wave_height.estimate_stack) fails rule 3, '' when it passes."""
    qc = qc or Qc()
    st = row.get("status", "")
    if st != "ok":
        note = str(row.get("note") or "").split(";")[0].strip()
        return f"status {st}" + (f" ({note})" if note else "")
    sh = _num(row.get("window_sharpness"))
    if not np.isfinite(sh):
        if qc.require_sharpness:
            return "no window sharpness (no OWG frame of the collection)"
    elif sh < qc.min_sharpness:
        return f"window sharpness {sh:g} < {qc.min_sharpness:g}"
    valid, peaks = _num(row.get("valid_fraction")), _num(row.get("n_peaks"))
    if not (valid >= qc.min_valid and peaks >= qc.min_peaks):
        return f"edge in {valid:g} of the rows, {peaks:g} uprushes"
    off = _num(row.get("off_profile_fraction"))
    if not off <= qc.max_off_profile:
        return f"off profile {off:g} > {qc.max_off_profile:g}"
    xf = _num(row.get("z_extrapolated_fraction"))
    if not xf <= qc.max_extrapolated:
        return f"extrapolated {xf:g} > {qc.max_extrapolated:g}"
    if not np.isfinite(_num(row.get("setup_m"))):
        return "no setup"
    src = str(row.get("water_level_source") or "")
    if src not in qc.level_sources:
        return f"still water from {src or 'nothing'} (not {'/'.join(qc.level_sources)})"
    return ""


def line_result(row, series, qc=None):
    """One line: its setup, sampling error, QC verdict, still water and diagnostics (the profile term is
    filled in by measure_stack)."""
    return {"line": int(row["line"]), "setup": _num(row.get("setup_m")), "why": qc_line(row, qc),
            "se": mean_sampling_se(series) if series is not None else np.nan,
            "wl": _num(row.get("water_level_m")), "wl_src": str(row.get("water_level_source") or ""),
            "contrast": _num(row.get("edge_contrast")), "sharp": _num(row.get("window_sharpness")),
            "gain": np.nan, "sig_dem": np.nan}


def combine_lines(results, min_lines=MIN_LINES):
    """
    The collection from its lines (rules 4 and 5): -> dict setup, sigma, sampling, lines (the between-line
    term), profile and gain (the profile term), n, used, status, note. Only the lines that pass QC count.
    """
    good = [r for r in results if not r["why"] and np.isfinite(r["setup"])]
    out = {"setup": np.nan, "sigma": np.nan, "sampling": np.nan, "lines": np.nan, "profile": np.nan,
           "gain": np.nan, "n": len(good), "used": "+".join(str(r["line"]) for r in good), "status": "failed qc",
           "note": ""}
    if len(good) < max(int(min_lines), 1):
        out["note"] = (f"{len(good)} line(s) passed, {max(int(min_lines), 1)} needed" if good else "no line passed")
        return out
    s = np.array([r["setup"] for r in good], float)
    se = np.array([r["se"] for r in good], float)
    out["setup"] = float(np.median(s))
    out["sampling"] = float(np.mean(se)) if np.isfinite(se).all() else np.nan
    g = np.array([r.get("gain", np.nan) for r in good], float)
    p = g * np.array([r.get("sig_dem", np.nan) for r in good], float)
    out["gain"] = float(np.median(g[np.isfinite(g)])) if np.isfinite(g).any() else np.nan
    out["profile"] = float(np.median(p[np.isfinite(p)])) if np.isfinite(p).any() else np.nan
    notes = []
    if len(good) >= 2:
        out["lines"] = float(np.std(s, ddof=1) / np.sqrt(len(s)))
        out["sigma"] = float(np.hypot(out["sampling"], out["lines"]))
    else:
        out["sigma"] = out["sampling"]
        notes.append("one line: setup_sigma_m is the sampling term alone")
    if not np.isfinite(out["sampling"]):
        notes.append("no sampling error (series too gappy)")
    out["status"] = "ok"
    out["note"] = "; ".join(notes)
    return out


def collection_row(name, start_epoch, results, info, min_lines=MIN_LINES):
    """The output row (FIELDS) of one collection."""
    c = combine_lines(results, min_lines)
    good = [r for r in results if not r["why"]]
    lv = next((r for r in good + results if np.isfinite(r["wl"])), None)
    sh = next((r["sharp"] for r in results if np.isfinite(r.get("sharp", np.nan))), np.nan)
    notes = [c["note"]] if c["note"] else []
    if info.get("why"):
        notes.append(info["why"])
    return {"collection_epoch": int(start_epoch),
            "time_utc": datetime.fromtimestamp(int(start_epoch), tz=timezone.utc).isoformat(),
            "filename": name, "status": c["status"], "note": "; ".join(notes),
            "setup_m": _r(c["setup"], 3), "setup_sigma_m": _r(c["sigma"], 3), "sigma_sampling_m": _r(c["sampling"], 3),
            "sigma_lines_m": _r(c["lines"], 3), "sigma_profile_m": _r(c["profile"], 3),
            "profile_gain": _r(c["gain"], 3), "n_lines": c["n"], "lines_used": c["used"],
            "line_setups": ";".join(f"{r['line']}:{r['setup']:.3f}" for r in results if np.isfinite(r["setup"])),
            "line_qc": ";".join(f"{r['line']}:{r['why'] or 'ok'}" for r in results),
            "line_edge_contrast": ";".join(f"{r['line']}:{r['contrast']:.3f}" for r in results
                                           if np.isfinite(r.get("contrast", np.nan))),
            "window_sharpness": _r(sh, 2),
            "water_level_m": _r(lv["wl"], 3) if lv else "", "water_level_source": lv["wl_src"] if lv else "",
            "dem": info.get("dem", ""), "dem_window": info.get("dem_window", ""),
            "dem_age_days": info.get("dem_age_days", ""), "dem_sha256": info.get("dem_sha256", ""),
            "dem_setup_share": info.get("dem_setup_share", "")}


# ---------------------------------------------------------------------
# One stack (a collection)
# ---------------------------------------------------------------------

def make_context(args, profiles, wl=None):
    """timestack_wave_height's context with rule 1's settings: the profile source `profiles` (PriorDem or
    FixedProfile), no plane, the edge-contrast label off, the placeholder period; this module's QC; and
    rule 6's wait for the GPS level wherever the GPS is an allowed source."""
    targs = thw.build_parser().parse_args([])
    targs.camera, targs.pix, targs.io, targs.eo = args.camera, args.pix, args.io, args.eo
    targs.dem, targs.beta = "", None
    targs.max_extrapolate, targs.max_retreat_speed = args.max_extrapolate, args.max_retreat_speed
    targs.wait_days, targs.owg_csv = args.wait_days, args.owg_csv
    targs.level_csv, targs.no_gps = args.level_csv, args.no_gps
    targs.gnssr_spline, targs.gauge_csv = args.gnssr_spline, args.gauge_csv
    targs.period_csv = targs.waves_csv = targs.buoy_csv = ""
    ctx = thw.make_context(targs, wl)
    ctx.dems = profiles
    ctx.periods = NoPeriod()
    ctx.min_contrast = float("-inf")       # rule 1b: off (a Michelson contrast can be negative)
    ctx.series_out = {}
    ctx.qc = Qc(args.min_sharpness, not args.allow_no_sharpness, MAX_OFF_PROFILE, args.max_extrapolated,
                ctx.rargs.min_valid, ctx.rargs.min_peaks, resolve_level_sources(args))
    # rule 6: a young stack waits for the GPS level where the GPS is an allowed source -- whether or not a
    # GPS series could be read this run (timestack_wave_height waits only where one was loaded)
    ctx.want_gps = "gps" in ctx.qc.level_sources
    ctx.min_lines = args.min_lines
    return ctx


def resolve_level_sources(args):
    """--level-sources, else the label of the level in use: 'level csv' for --level-csv, 'chatham' for
    --no-gps, else 'gps'."""
    if getattr(args, "level_sources", None):
        return tuple(args.level_sources)
    if getattr(args, "level_csv", None):
        return ("level csv",)
    if getattr(args, "no_gps", False):
        return ("chatham",)
    return LEVEL_SOURCES


def measure_stack(path, ctx, lines=DEFAULT_LINES, now=None, captures=None):
    """
    One collection's row (collection_row), 'wait' (an input that arrives late is missing and the stack is
    young: rule 6), or None (no line of the stack could be read: a later run tries again). Each line is
    timestack_wave_height.estimate_stack, its swash series captured for the sampling error and the profile
    term (captures: a dict to keep each line's capture in, for the audit).
    """
    path = Path(path)
    start = epoch_from_name(path.name)
    results = []
    for ln in lines:
        ctx.series_out = {}
        r = thw.estimate_stack(path, ln, ctx, now=now)
        if r == "wait":
            return "wait"
        if r is None:
            continue
        cap = ctx.series_out.get((r["filename"], int(ln)))
        res = line_result(r, cap["eta"] if cap is not None else None, ctx.qc)
        if cap is not None and not res["why"] and cap.get("cam") in ctx.pix_cache:
            io, eo = ctx.pix_cache[cap["cam"]]
            sig = ctx.dems.sigma_for(start + BURST_MID_OFFSET_S) if hasattr(ctx.dems, "sigma_for") else None
            res["gain"], res["sig_dem"] = profile_term(cap, io, eo, sig)
        if captures is not None and cap is not None:
            captures[int(ln)] = cap
        results.append(res)
    if not results:
        return None
    return collection_row(path.name, start, results, ctx.dems.info(start + BURST_MID_OFFSET_S), ctx.min_lines)


# ---------------------------------------------------------------------
# The output file: appended, never rewritten except by --reprocess
# ---------------------------------------------------------------------

def read_done(path):
    rows = thw._read_rows(path)
    return set() if not rows else {r["filename"] for r in rows}


def check_header(path, reprocess=False):
    """A non-empty file whose header is not FIELDS is refused (SystemExit) unless reprocess. -> it exists."""
    path = Path(path)
    exists = path.exists() and path.stat().st_size > 0
    if exists:
        with open(path, newline="") as f:
            header = next(csv.reader(f), [])
        if header != FIELDS and not reprocess:
            raise SystemExit(f"{path}: its columns are not this version's; rerun with --reprocess")
    return exists


def write_rows(path, rows, reprocess=False):
    """
    Appends rows (header first if new); a collection (filename) already in the file is not appended
    again. With reprocess the file is rewritten: the rows these redo replaced, the others kept, all in
    collection_epoch order. -> rows written (0 for no rows: then no file is created).
    """
    path = Path(path)
    if not rows:
        return 0
    rows = sorted(rows, key=lambda r: int(r["collection_epoch"]))
    exists = check_header(path, reprocess)
    if reprocess and exists:
        redo = {r["filename"] for r in rows}
        kept = [r for r in (thw._read_rows(path) or []) if r.get("filename") not in redo]
        merged = sorted(kept + rows, key=lambda r: int(float(r["collection_epoch"])))
        tmp = path.with_suffix(path.suffix + ".tmp")
        with open(tmp, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=FIELDS, restval="", extrasaction="ignore")
            w.writeheader()
            w.writerows(merged)
        tmp.replace(path)
        return len(rows)
    done = read_done(path) if exists else set()
    uniq, seen = [], set(done)
    for r in rows:
        if r["filename"] not in seen:
            seen.add(r["filename"])
            uniq.append(r)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, restval="", extrasaction="ignore")
        if not exists:
            w.writeheader()
        w.writerows(uniq)
    return len(uniq)


class MeasuredSetups:
    """
    The collections of a measured-setup file whose status is ok, for extract_elevation_contours.py
    --setup-source: at(epoch) is the setup_m of the collection whose start lies nearest the frame's
    file-name epoch, within tol_s (MATCH_TOL_S), else None. ValueError for a file without the columns.
    """

    def __init__(self, path, tol_s=MATCH_TOL_S):
        with open(path, newline="") as f:
            rd = csv.DictReader(f)
            missing = [k for k in ("collection_epoch", "status", "setup_m") if k not in (rd.fieldnames or [])]
            if missing:
                raise ValueError(f"{path}: not a measured-setup file (no {', '.join(missing)})")
            rows = list(rd)
        self.path, self.tol, self.rows = str(path), float(tol_s), len(rows)
        pairs = {}
        for r in rows:
            e, s = _num(r.get("collection_epoch")), _num(r.get("setup_m"))
            if r.get("status") == "ok" and np.isfinite(e) and np.isfinite(s):
                pairs.setdefault(e, s)
        ep = np.array(sorted(pairs), float)
        self.ep, self.setup = ep, np.array([pairs[e] for e in ep], float)
        self.ok = len(ep)

    def nearest(self, epoch):
        """The start (epoch) of the ok collection nearest `epoch` within tol_s, else None."""
        if not len(self.ep):
            return None
        e = float(epoch)
        i = int(np.searchsorted(self.ep, e))
        best = None
        for j in (i - 1, i):
            if 0 <= j < len(self.ep) and abs(self.ep[j] - e) <= self.tol and \
                    (best is None or abs(self.ep[j] - e) < abs(self.ep[best] - e)):
                best = j
        return None if best is None else float(self.ep[best])

    def at(self, epoch):
        c = self.nearest(epoch)
        return None if c is None else float(self.setup[int(np.searchsorted(self.ep, c))])


# ---------------------------------------------------------------------
# Over the archive
# ---------------------------------------------------------------------

def output_path(args):
    return Path(args.output or HERE / "archive" / f"setup_measured_{args.camera}.csv")


def open_profiles(args):
    """The PriorDem of --dem (default: default_dem()), refusing a DEM file whose window is not known."""
    dem = args.dem or str(default_dem())
    profiles = PriorDem(dem, args.max_dem_age_days)
    if profiles.refused:
        raise SystemExit(f"--dem {dem}: not a dated DEM (dem_<window end>_<days>d_dem.asc), so the "
                         f"strictly-earlier rule cannot be applied; give the series folder (archive/dems_model "
                         f"or archive/dems)")
    return profiles


def report_refused(profiles):
    if profiles.not_model_built:
        print(f"  {len(profiles.not_model_built)} DEM(s) refused as not model-built (rule 2b):")
        for name, why in sorted(profiles.not_model_built.items())[-10:]:
            print(f"      {name}: {why}")


def run(args, wl=None, now=None):
    out = output_path(args)
    check_header(out, args.reprocess)                # refuse a file of other columns before any work
    profiles = open_profiles(args)
    ctx = make_context(args, profiles, wl)
    done = set() if args.reprocess else read_done(out)
    paths = sorted(Path(args.ras_dir).glob(f"*.{args.camera}.ras.tiff"), key=lambda q: q.name)
    lo = _day_epoch(args.since) if args.since else -np.inf
    hi = _day_epoch(args.until, end=True) if args.until else np.inf
    print(f"measured setup {args.camera}: {len(paths)} stack(s) in {args.ras_dir}; profile: the model-built DEM "
          f"whose window ends 1-{profiles.max_age} day(s) before each stack's day, from {profiles.what}")
    print(f"  QC per line: {ctx.qc.describe()}; a collection needs {args.min_lines} line(s) of "
          + "/".join(str(x) for x in args.lines))
    if hasattr(ctx.wl, "describe"):
        print("  " + ctx.wl.describe())
    chunk, waiting, n, counts, latest = [], 0, 0, {}, None
    try:
        for p in paths:
            e = epoch_from_name(p.name)
            if e is None or not lo <= e < hi or p.name in done:
                continue
            r = measure_stack(p, ctx, args.lines, now=now)
            if r is None:
                continue
            if r == "wait":
                waiting += 1
                continue
            chunk.append(r)
            counts[r["status"]] = counts.get(r["status"], 0) + 1
            if r["status"] == "ok":
                latest = r
            if len(chunk) >= CHUNK_ROWS:
                n += write_rows(out, chunk, reprocess=args.reprocess)
                chunk = []
    finally:
        if chunk:
            n += write_rows(out, chunk, reprocess=args.reprocess)
    print(f"measured setup {args.camera}: {n} new collection(s) -> {out}"
          + (" (" + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())) + ")" if counts else "")
          + (f"; {waiting} wait for a DEM, the water level or the GPS level" if waiting else "")
          + (f"; latest {latest['setup_m']} m +/- {latest['setup_sigma_m']} (profile term "
             f"{latest['sigma_profile_m'] or '-'}; {latest['time_utc'][:16]}Z, on {latest['dem']})" if latest else ""))
    report_refused(profiles)
    return 0


# ---------------------------------------------------------------------
# The audit (VALIDATION 4.12 criterion 4)
# ---------------------------------------------------------------------

def _pairs(text):
    out = {}
    for part in str(text or "").split(";"):
        if ":" in part:
            k, v = part.split(":", 1)
            out[k.strip()] = v.strip()
    return out


def _q(v, q):
    return float(np.percentile(v, q)) if len(v) else np.nan


def line_noise(rows):
    """|line 1 - line 2| (the first two lines used) of the ok collections where two lines passed, by the
    window's sharpness: {band: dict n, median, p90 of the difference, edge contrast median}."""
    bands = {}
    for r in rows:
        if r.get("status") != "ok":
            continue
        used = [x for x in str(r.get("lines_used") or "").split("+") if x]
        vals = _pairs(r.get("line_setups"))
        sh = _num(r.get("window_sharpness"))
        band = ("clear (>= %g)" % CLEAR_SHARPNESS if sh >= CLEAR_SHARPNESS else
                "low (%g-%g)" % (MIN_SHARPNESS, CLEAR_SHARPNESS) if np.isfinite(sh) else "unknown")
        b = bands.setdefault(band, {"diff": [], "contrast": []})
        if len(used) >= 2 and used[0] in vals and used[1] in vals:
            b["diff"].append(abs(float(vals[used[0]]) - float(vals[used[1]])))
        b["contrast"] += [_num(v) for v in _pairs(r.get("line_edge_contrast")).values() if np.isfinite(_num(v))]
    out = {}
    for band, b in sorted(bands.items()):
        d, c = np.array(b["diff"], float), np.array(b["contrast"], float)
        out[band] = {"n": int(len(d)), "median": _q(d, 50), "p90": _q(d, 90),
                     "contrast_median": _q(c, 50), "contrast_n": int(len(c))}
    return out


def fail_reasons(rows):
    """The QC failures of every line, with their numbers blanked: {reason: lines}."""
    import re
    out = {}
    for r in rows:
        for why in _pairs(r.get("line_qc")).values():
            if why and why != "ok":
                key = re.sub(r"[-+]?\d+\.\d+|\b\d+\b", "#", why)
                out[key] = out.get(key, 0) + 1
    return out


def audit(args, wl=None):
    """--audit N: the lines' disagreement by window sharpness, the QC's failures, and N passed collections
    drawn with the edge on (module docstring, THE AUDIT). Writes only --audit-dir."""
    import runup_from_timestack as rft
    out = output_path(args)
    rows = thw._read_rows(out) or []
    if not rows:
        print(f"audit: no rows in {out} (run measured_setup.py first)")
        return 4
    ok = [r for r in rows if r.get("status") == "ok"]
    print(f"audit of {out}: {len(rows)} collection(s), {len(ok)} passed QC")
    print("  |line 1 - line 2| where both passed, by window sharpness (owg_live.py's): n, median, p90 (m); "
          "edge contrast median (rule 1b, a diagnostic)")
    for band, b in line_noise(rows).items():
        print(f"      {band:<14} {b['n']:5d}  {b['median']:.3f}  {b['p90']:.3f}   contrast {b['contrast_median']:+.3f} "
              f"({b['contrast_n']} line(s))")
    print("  QC failures (lines):")
    for why, k in sorted(fail_reasons(rows).items(), key=lambda kv: -kv[1])[:12]:
        print(f"      {k:5d}  {why}")
    # rule 2b afterwards: every DEM a passed collection stood on, read again from its _info.json
    folder = Path(args.dem or default_dem())
    folder = folder if folder.is_dir() else folder.parent
    stood = sorted({r["dem"] for r in ok if r.get("dem")})
    bad = {}
    for name in stood:
        p = folder / name
        good, why, _share = dem_provenance(p) if p.is_file() else (False, f"not in {folder} now", np.nan)
        if not good:
            bad[name] = why
    print(f"  rule 2b: the passed collections stood on {len(stood)} DEM(s) of {folder}; not model-built now: "
          f"{len(bad)} ({sum(1 for r in ok if r.get('dem') in bad)} collection(s); must be 0)")
    for name, why in sorted(bad.items())[:10]:
        print(f"      {name}: {why}")
    n = min(int(args.audit), len(ok))
    if n <= 0:
        return 0
    rng = np.random.default_rng(int(args.audit_seed))
    sample = [ok[i] for i in sorted(rng.choice(len(ok), n, replace=False).tolist())]
    adir = Path(args.audit_dir)
    adir.mkdir(parents=True, exist_ok=True)
    profiles = open_profiles(args)
    args.wait_days = 0.0
    ctx = make_context(args, profiles, wl)
    index = []
    for r in sample:
        path = Path(args.ras_dir) / r["filename"]
        if not path.is_file():
            print(f"      {r['filename']}: not in {args.ras_dir}")
            continue
        caps = {}
        again = measure_stack(path, ctx, [int(x) for x in str(r["lines_used"]).split("+") if x], captures=caps)
        for ln, cap in sorted(caps.items()):
            png = adir / f"{Path(r['filename']).stem}.line{ln}.png"
            rft.save_plot(png, cap["anomaly"], cap["edge"], cap["threshold"])
            index.append({"filename": r["filename"], "time_utc": r["time_utc"], "line": ln,
                          "setup_m": r["setup_m"], "setup_now_m": again["setup_m"] if isinstance(again, dict) else "",
                          "window_sharpness": r.get("window_sharpness", ""),
                          "edge_contrast": _pairs(r.get("line_edge_contrast")).get(str(ln), ""), "png": png.name})
    with open(adir / "index.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["filename", "time_utc", "line", "setup_m", "setup_now_m",
                                          "window_sharpness", "edge_contrast", "png"])
        w.writeheader()
        w.writerows(index)
    print(f"  {len(index)} line(s) of {n} passed collection(s) (seed {args.audit_seed}) drawn with the edge in red "
          f"-> {adir} (index.csv): look at every one; the edge must follow the uprush tips")
    return 0


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--ras-dir", default=str(HERE / "archive" / "ras_c2"))
    ap.add_argument("--camera", default="c2")
    ap.add_argument("--lines", nargs="+", type=int, default=list(DEFAULT_LINES),
                    help="timestack lines (default 1 2: the c2 lines that cross the swash)")
    ap.add_argument("--pix", default=None, help="as timestack_wave_height.py")
    ap.add_argument("--io", default=None, help="as timestack_wave_height.py")
    ap.add_argument("--eo", default=None, help="as timestack_wave_height.py")
    ap.add_argument("--output", default=None, help="default archive/setup_measured_<cam>.csv")
    ap.add_argument("--dem", default=None,
                    help="the dated DEMs (dem_from_contours.py --series-dir), or one dated DEM file: the profile "
                         "is the model-built one whose window ends before the stack's day (rule 2). Default: "
                         "archive/dems_model when it exists (the cron's model-only series), else archive/dems")
    ap.add_argument("--max-dem-age-days", type=int, default=MAX_DEM_AGE_DAYS,
                    help=f"latest the DEM's window may end before the stack's day (default {MAX_DEM_AGE_DAYS})")
    ap.add_argument("--max-extrapolate", type=float, default=thw.MAX_EXTRAPOLATE_M,
                    help="m the DEM profile is carried on past its edge (timestack_wave_height.py's)")
    ap.add_argument("--max-retreat-speed", type=float, default=None, help="as timestack_wave_height.py")
    ap.add_argument("--owg-csv", default=str(HERE / "archive" / "owg_c2_H.csv"),
                    help="owg_live.py's c2 frames: the window sharpness per collection")
    ap.add_argument("--min-sharpness", type=float, default=MIN_SHARPNESS,
                    help=f"window sharpness a line needs (default {MIN_SHARPNESS:g})")
    ap.add_argument("--allow-no-sharpness", action="store_true",
                    help="pass a line whose collection has no OWG frame (a period without the OWG record)")
    ap.add_argument("--max-extrapolated", type=float, default=MAX_EXTRAPOLATED,
                    help=f"largest share of the edge on the profile's extension (default {MAX_EXTRAPOLATED:g})")
    ap.add_argument("--level-sources", nargs="+", default=None,
                    help="still-water sources a line may use (default: gps; 'level csv' with --level-csv, "
                         "'chatham' with --no-gps)")
    ap.add_argument("--min-lines", type=int, default=MIN_LINES,
                    help=f"lines that must pass for a collection (default {MIN_LINES})")
    ap.add_argument("--wait-days", type=float, default=WAIT_DAYS,
                    help=f"leave stacks younger than this without a DEM or water level for later (default {WAIT_DAYS:g})")
    ap.add_argument("--level-csv", default=None, help="as timestack_wave_height.py")
    from marconi_water_level import GNSSR_SPLINE, GAUGE_CSV
    ap.add_argument("--gnssr-spline", default=GNSSR_SPLINE)
    ap.add_argument("--gauge-csv", default=GAUGE_CSV)
    ap.add_argument("--no-gps", action="store_true", help="water level from the Chatham transfer only")
    ap.add_argument("--since", default=None, help="only stacks from this UTC day (YYYY-MM-DD)")
    ap.add_argument("--until", default=None, help="only stacks up to this UTC day (YYYY-MM-DD)")
    ap.add_argument("--reprocess", action="store_true", help="redo the stacks and replace their rows")
    ap.add_argument("--audit", type=int, default=None, metavar="N",
                    help=f"no measuring: the audit of the output file, N passed collections drawn (e.g. {AUDIT_SAMPLE})")
    ap.add_argument("--audit-dir", default=str(HERE / "archive" / "measured_setup_audit"))
    ap.add_argument("--audit-seed", type=int, default=0)
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.audit is not None:
        return audit(args)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())

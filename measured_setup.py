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
          contrast' (it cannot tell blur from a clear day there); the
          window's sharpness (rule 3) is the blur check instead;
      (c) the period is a placeholder (PLACEHOLDER_T0_S): the setup is a
          mean, with no period in it; the period gates only that script's
          inversion for H0, which is neither used nor written here. So a
          collection has a measured setup in a buoy outage too, where the
          model has none (finding H2 of the review register: no Tp, no
          setup).
 2. THE DEM IS STRICTLY EARLIER (the non-circularity rule). The swash
    elevation is read off a DEM, and a DEM of timex waterlines carries the
    setup its lines were given: on the DEM of its own week the measured
    setup tends to the model's C sqrt(Hs L0) itself (timestack_wave_height
    .py, CIRCULARITY), i.e. it echoes the model. So the profile comes ONLY
    from a DEM built from days strictly before the stack's UTC day: of the
    dated DEMs of archive/dems (dem_<window end>_<days>d_dem.asc,
    dem_from_contours.py --series-dir; the window is its last <days> UTC
    days), the one whose window ENDS latest while still before that day,
    and at most --max-dem-age-days (MAX_DEM_AGE_DAYS) before it; none, no
    measurement ('no slope'). A DEM whose window holds the stack's day is
    never used -- the DEM being built from the frame's own window included
    -- so no frame's elevation is ever computed from a DEM that holds that
    frame or any frame of its day. A single DEM file is accepted only under
    a series name (its window must be known). Each row records the DEM,
    its window, its age (days from the window's end to the stack's day) and
    its sha256 (the first 12 hex), so the rule can be checked afterwards.
    WHAT THAT BUYS, AND WHAT IT DOES NOT. Live, the measured setup then
    corrects the model's DAY-TO-DAY departures, those the offshore Hs and Tp
    cannot see: a buoy outage, the wave direction (sheltered south-west
    seas), groupiness and infragravity energy. It INHERITS the long-term
    level of the DEM it stands on: that DEM's lines carried the model's
    setup (or, once measured setups are in use, those of still earlier
    DEMs), so a chain of DEMs built with measured setups keeps the level of
    its first, model-built, DEMs; it cannot correct C's long-term level
    (the survey-date check, evaluate_measured_setup.py, can). And it reads
    BEACH CHANGE since that DEM's window as setup: the edge is put on the
    OLD beach, so where the beach has eroded by d since, every ray meets
    the old (higher) surface first and the setup reads high by d / (1 -
    beta / tan(delta)) on a plane of slope beta seen at a ray depression
    delta (on the synthetic scene of test_measured_setup.py, real c2
    calibration: a DEM 0.10 m above the beach read the setup 0.13 / 0.15 /
    0.19 m higher on a 1:14 / 1:10 / 1:8 beach, i.e. tan(delta) ~0.29
    there; 0.10 m below, as much lower [SYNTH: the 1:10 beach is the test's,
    1:14 and 1:8 a run made while writing it]). Given to the waterline of the same
    collection, which lies near the mean swash position, that setup puts
    the line back onto the old surface: AT THE TIMESTACK LINES a DEM built
    with measured setups follows the DEM it stood on, not the beach, so its
    change from week to week is suppressed there (the rest of the view
    keeps its change relative to the lines). That is why the DEM is the most
    recent earlier one, why the default stays the model, and why a live
    trial is scored against a survey after a storm (VALIDATION 4.12).
 3. QC, per line, all of:
      status ok        timestack_wave_height's own: not dark, the window not
                       blurred, an edge in >= 80% of the rows with >= 20
                       uprushes, a profile, an estimate;
      window_sharpness >= 10 (station_status.LENS_BLURRED; owg_live.py's
                       sharpness of the c2 bright frame of the same
                       collection: clear days read 49-88, the blurred window
                       of 20-24 Oct 2026 3.4-4.7). A stack with no such frame
                       fails: its window is unknown (--allow-no-sharpness for
                       a period without the OWG record, e.g. 2025);
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
 5. ITS UNCERTAINTY, setup_sigma_m = sqrt(sampling^2 + lines^2):
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
                with ~30% scatter of its own per record. What lies below 1/T a single
                10-minute record cannot see. The lines' sampling errors are
                taken as FULLY correlated (the infragravity swash is coherent
                over the tens of metres between the lines): their mean, not
                reduced by averaging the lines;
      lines     the standard deviation of the lines' setups / sqrt(n) (for
                two lines half their difference): their alongshore
                difference and their independent errors. With one line it is
                unknown and setup_sigma_m is the sampling term alone (noted).
    NOT in it: the profile's level error (the DEM's, or a survey's), common
    to every collection standing on that profile (a bias, not noise; rule
    2); nor the detector's own bias (on the synthetic scenes the edge trails
    a fast tip: setup +0.04 to +0.10 m, timestack_wave_height.py and
    test_measured_setup.py [SYNTH]).
 6. WAITING. A stack younger than --wait-days (3) whose DEM, water level or
    (where the GPS is the source in use) GPS level is missing is left for a
    later run, as timestack_wave_height.py does: the DEM of the day before a
    stack is built by a later cron run, the GPS level arrives 1-2 days late.
    Older, it is written as it is (failed qc, the reason in its note and
    line_qc), and only --reprocess measures it again.

HOW THE EXTRACTOR USES IT (extract_elevation_contours.py --setup-source
measured | auto --setup-measured FILE). A frame whose file-name epoch lies
within MATCH_TOL_S (300 s) of a collection's start gets that collection's
setup_m, when its status is ok (timex, snap and stack of one collection
share the file-name epoch, timex_timing.py; 300 s is half the collection):
the frames of both cameras, one value per collection, as the model gives
one per frame. 'measured': a frame without one gets no setup (blank, as a
frame without Tp); 'auto': the model's. Its setup_source column says which.

WHAT IT CANNOT DO
  * The inheritance and the beach change of rule 2.
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
    (measured).
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
  time_utc, filename, status (ok / failed qc), note, setup_m, setup_sigma_m,
  sigma_sampling_m, sigma_lines_m, n_lines (lines that passed), lines_used,
  line_setups ('1:0.312;2:0.295', every line with a value), line_qc ('1:ok;
  2:window sharpness 4.2 < 10'), water_level_m, water_level_source, dem,
  dem_window (first..last day), dem_age_days, dem_sha256

Usage:
    python3 measured_setup.py                       new stacks of archive/ras_c2
    python3 measured_setup.py --since 2026-10-01    from that UTC day
    python3 measured_setup.py --reprocess           every stack again
Tests: python3 test_measured_setup.py
"""

import sys
import csv
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
          "sigma_sampling_m", "sigma_lines_m", "n_lines", "lines_used", "line_setups", "line_qc",
          "water_level_m", "water_level_source", "dem", "dem_window", "dem_age_days", "dem_sha256"]
STATUSES = ("ok", "failed qc")

DEFAULT_LINES = thw.DEFAULT_LINES     # 1 and 2: the c2 lines that cross the swash
MAX_DEM_AGE_DAYS = 3                  # the DEM's window ends 1-3 days before the stack's day
MIN_SHARPNESS = LENS_BLURRED          # 10: owg_live.py's 'blurred' (station_status.py)
MAX_OFF_PROFILE = thw.MISSING_ELEV_MAX    # 0.02: timestack_wave_height's own limit for a setup
MAX_EXTRAPOLATED = 0.10               # share of the edge on the profile's extension: PROVISIONAL (QC)
MIN_VALID, MIN_PEAKS = 0.8, 20        # the runup script's own test (runup_from_timestack.py defaults)
LEVEL_SOURCES = ("gps",)              # the still water the contour extractor gives the frames
MIN_LINES = 1
MATCH_TOL_S = float(TIMEX_MID_OFFSET_S)   # 300 s: a frame and a stack of one collection (half of it)
SAMPLING_F_LOW_HZ = 0.005             # the lowest bins taken as the spectral level at zero frequency
PLACEHOLDER_T0_S = 10.0               # gates timestack_wave_height's inversion only (rule 1c)
WAIT_DAYS = thw.WAIT_DAYS
CHUNK_ROWS = 25                       # collections per append (a crash loses at most these)


# ---------------------------------------------------------------------
# The profile: the strictly-earlier DEM (or one fixed surface)
# ---------------------------------------------------------------------

def utc_day(epoch):
    return datetime.fromtimestamp(float(epoch), tz=timezone.utc).date()


def file_sha256(path, n=12):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:n]


class PriorDem:
    """
    The DEM a stack's swash is put on (rule 2): of the dated DEMs of `source` (a folder of
    dem_<window end>_<days>d_dem.asc, or one such file), the one whose window ends latest while still
    BEFORE the stack's UTC day, and at most max_age_days before it. for_epoch(epoch) is what
    timestack_wave_height.estimate_stack calls (-> its Dem or None); info(epoch) what a row records.
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
        self._dems, self._sha = {}, {}

    def pick(self, epoch):
        """-> (path, {dem, dem_window, dem_age_days, dem_sha256}) or (None, why none)."""
        day = utc_day(epoch)
        cands = [(end, days, f) for end, days, f in self.series if end < day and (day - end).days <= self.max_age]
        if not cands:
            return None, f"no DEM whose window ends 1-{self.max_age} day(s) before {day}"
        end, days, f = max(cands, key=lambda c: (c[0], -c[1], c[2].name))   # latest end; then the shorter
        first = end - timedelta(days=days - 1)
        if f not in self._sha:
            self._sha[f] = file_sha256(f)
        return f, {"dem": f.name, "dem_window": f"{first}..{end}", "dem_age_days": (day - end).days,
                   "dem_sha256": self._sha[f]}

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

    def info(self, epoch):
        f, info = self.pick(epoch)
        if f is None:
            return {"dem": "", "dem_window": "", "dem_age_days": "", "dem_sha256": "", "why": info}
        return info


class FixedProfile:
    """One surface for every stack (a survey: evaluate_measured_setup.py), and what a row records of it."""

    def __init__(self, dem, label, sha=""):
        self.dem, self.label, self.sha = dem, label, sha
        self.what = f"{dem.name} ({label})"

    def for_epoch(self, epoch):
        return self.dem

    def info(self, epoch):
        return {"dem": self.dem.name, "dem_window": self.label, "dem_age_days": "", "dem_sha256": self.sha}


class NoPeriod:
    """A placeholder period for every stack (rule 1c): estimate_stack's inversion needs one, the setup does
    not; never pending (the camera's period is not waited for)."""

    def for_stack(self, filename, mid_epoch):
        return PLACEHOLDER_T0_S, "placeholder (not used: the setup has no period in it)", False


# ---------------------------------------------------------------------
# Per line: QC and the sampling error; per collection: the median
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


def _num(x):
    return thw._num(x)


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
    """One line: its setup, sampling error, QC verdict and still water."""
    return {"line": int(row["line"]), "setup": _num(row.get("setup_m")), "why": qc_line(row, qc),
            "se": mean_sampling_se(series) if series is not None else np.nan,
            "wl": _num(row.get("water_level_m")), "wl_src": str(row.get("water_level_source") or "")}


def combine_lines(results, min_lines=MIN_LINES):
    """
    The collection from its lines (rules 4 and 5): -> dict setup, sigma, sampling, lines (the between-line
    term), n, used, status, note. Only the lines that pass QC count.
    """
    good = [r for r in results if not r["why"] and np.isfinite(r["setup"])]
    out = {"setup": np.nan, "sigma": np.nan, "sampling": np.nan, "lines": np.nan, "n": len(good),
           "used": "+".join(str(r["line"]) for r in good), "status": "failed qc", "note": ""}
    if len(good) < max(int(min_lines), 1):
        out["note"] = (f"{len(good)} line(s) passed, {max(int(min_lines), 1)} needed" if good else "no line passed")
        return out
    s = np.array([r["setup"] for r in good], float)
    se = np.array([r["se"] for r in good], float)
    out["setup"] = float(np.median(s))
    out["sampling"] = float(np.mean(se)) if np.isfinite(se).all() else np.nan
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


def _r(x, nd):
    return thw._r(x, nd)


def collection_row(name, start_epoch, results, info, min_lines=MIN_LINES):
    """The output row (FIELDS) of one collection."""
    c = combine_lines(results, min_lines)
    good = [r for r in results if not r["why"]]
    lv = next((r for r in good + results if np.isfinite(r["wl"])), None)
    notes = [c["note"]] if c["note"] else []
    if info.get("why"):
        notes.append(info["why"])
    return {"collection_epoch": int(start_epoch),
            "time_utc": datetime.fromtimestamp(int(start_epoch), tz=timezone.utc).isoformat(),
            "filename": name, "status": c["status"], "note": "; ".join(notes),
            "setup_m": _r(c["setup"], 3), "setup_sigma_m": _r(c["sigma"], 3), "sigma_sampling_m": _r(c["sampling"], 3),
            "sigma_lines_m": _r(c["lines"], 3), "n_lines": c["n"], "lines_used": c["used"],
            "line_setups": ";".join(f"{r['line']}:{r['setup']:.3f}" for r in results if np.isfinite(r["setup"])),
            "line_qc": ";".join(f"{r['line']}:{r['why'] or 'ok'}" for r in results),
            "water_level_m": _r(lv["wl"], 3) if lv else "", "water_level_source": lv["wl_src"] if lv else "",
            "dem": info.get("dem", ""), "dem_window": info.get("dem_window", ""),
            "dem_age_days": info.get("dem_age_days", ""), "dem_sha256": info.get("dem_sha256", "")}


# ---------------------------------------------------------------------
# One stack (a collection)
# ---------------------------------------------------------------------

class Context:
    pass


def make_context(args, profiles, wl=None):
    """timestack_wave_height's context with rule 1's settings: the profile source `profiles` (PriorDem or
    FixedProfile), no plane, the edge-contrast label off, the placeholder period; and this module's QC."""
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


def measure_stack(path, ctx, lines=DEFAULT_LINES, now=None):
    """
    One collection's row (collection_row), 'wait' (an input that arrives late is missing and the stack is
    young: rule 6), or None (no line of the stack could be read: a later run tries again). Each line is
    timestack_wave_height.estimate_stack, its swash series captured for the sampling error.
    """
    path = Path(path)
    results = []
    for ln in lines:
        ctx.series_out = {}
        r = thw.estimate_stack(path, ln, ctx, now=now)
        if r == "wait":
            return "wait"
        if r is None:
            continue
        series = ctx.series_out.get((r["filename"], int(ln)))
        results.append(line_result(r, series, ctx.qc))
    if not results:
        return None
    start = epoch_from_name(path.name)
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
    collection_epoch order. -> rows written.
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

def _day_epoch(s, end=False):
    t = datetime.strptime(s, "%Y-%m-%d").replace(tzinfo=timezone.utc).timestamp()
    return t + 86400.0 if end else t


def run(args, wl=None, now=None):
    out = Path(args.output or HERE / "archive" / f"setup_measured_{args.camera}.csv")
    check_header(out, args.reprocess)                # refuse a file of other columns before any work
    profiles = PriorDem(args.dem, args.max_dem_age_days)
    if profiles.refused:
        raise SystemExit(f"--dem {args.dem}: not a dated DEM (dem_<window end>_<days>d_dem.asc), so the "
                         f"strictly-earlier rule cannot be applied; give the series folder (archive/dems)")
    ctx = make_context(args, profiles, wl)
    done = set() if args.reprocess else read_done(out)
    paths = sorted(Path(args.ras_dir).glob(f"*.{args.camera}.ras.tiff"), key=lambda q: q.name)
    lo = _day_epoch(args.since) if args.since else -np.inf
    hi = _day_epoch(args.until, end=True) if args.until else np.inf
    print(f"measured setup {args.camera}: {len(paths)} stack(s) in {args.ras_dir}; profile: the DEM whose window "
          f"ends 1-{profiles.max_age} day(s) before each stack's day, from {profiles.what}")
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
          + (f"; latest {latest['setup_m']} m +/- {latest['setup_sigma_m']} ({latest['time_utc'][:16]}Z, on "
             f"{latest['dem']})" if latest else ""))
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
    ap.add_argument("--dem", default=str(HERE / "archive" / "dems"),
                    help="the dated DEMs (dem_from_contours.py --series-dir; the cron's archive/dems), or one "
                         "dated DEM file: the profile is the one whose window ends before the stack's day")
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
    return ap


def main(argv=None):
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())

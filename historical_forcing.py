#!/usr/bin/env python3
"""
Historical Water Level And Waves For A Survey Date
====================================================
Writes the still-water level and the waves at Marconi for any period in
Oct 2024 - Mar 2025, in the files the waterline chain reads, with an
error estimate on every row and a report of how each number was made:

    DIR/water_level.csv   time, water_level_navd88, source, sigma_m
    DIR/waves.csv         time_utc, epoch, wvht_m, dpd_s, mwd_deg, source,
                          hs_sigma_m, tp_sigma_s   (extract_elevation_contours.py --waves)
    DIR/forcing_report.txt, DIR/forcing.png, DIR/forcing.json

WHY. The survey dates before the GNSS-R record (Oct 2024, Jan and Mar
2025) need a water level and waves for every photo. The Signature 1000
ADCP measured both off Marconi from 2024-12-09 18:00 to 2025-03-10 15:00
UTC; outside it the only water level is the NOAA Chatham gauge (8447435,
inside Pleasant Bay / Chatham Harbor, 20 km south), whose tide is smaller
and later than Marconi's, and the only waves are a hindcast (WIS ST63064)
and the Boston buoy (NDBC 44013). Every frame's elevation is the water
level plus the wave setup, so an error here moves the whole DEM.

HOW: WATER LEVEL
  * ADCP where it has a reading within --max-gap-minutes on both sides of
    a time (the ADCP's NAVD88 datum is water_level + an offset DERIVED
    FROM CHATHAM's mean over the deployment, adcp_to_navd88.py: Marconi's
    mean level is assumed equal to Chatham's. The report states this and
    cross-checks it against the tide-model MSL-to-NAVD88 constant).
  * Elsewhere Chatham converted to Marconi, fitted on the ADCP overlap
    (--fit-start/--fit-end; default all of it). Two transfers compete:
      linear    Marconi(t) = a * Chatham(t - lag) + b
      harmonic  least-squares harmonic analysis of both records over the
                fit window (the constituents the window resolves by the
                Rayleigh criterion, at most M2 S2 N2 K1 O1 Q1 L2 MU2 M3 MK3
                MN4 M4 MS4 M6 2MS6); per-constituent amplitude ratio and
                phase lag applied to Chatham's OWN harmonic fit over the
                target period (+/- 30 days), plus Chatham's mean over that
                period, plus Chatham's non-tidal residual scaled by a fitted
                factor (beta, at a fitted lag). Unlike the linear transfer,
                the seasonal mean level is carried 1:1, not multiplied by a,
                and the harbour's own overtides are not amplified.
    Both are scored by leaving out one ISO week at a time (CV) and by two
    EXTRAPOLATION tests, because October lies outside the overlap:
      fit 2025-01-01..03-10, predict 2024-12-09..12-31;
      fit 2024-12-09..2025-02-09, predict 2025-02-10..03-10.
    The transfer with the lower extrapolation RMS is used; each statistic
    is also given for daytime hours only (13.5-18 UTC, the hours whose
    photos are used). sigma_m of a transferred reading is the chosen
    method's extrapolation RMS (the larger of all-hours and daytime);
    sigma_m of an ADCP reading is its white-noise level. Neither includes
    the datum (shared by every row; see the report).

HOW: WAVES, AND THE "CURRENCY" OF THE SETUP COEFFICIENT
  C = 0.037 (waterline_timex_cron.sh, fitted 8 Oct 2026 on 29 Sep - 5 Oct
  2026 frames) was fitted on contours tagged from archive/waves_marconi.csv
  (USE_MARCONI_WAVES=1 since 2 Oct 2026; the contours are rebuilt from the
  whole archive every run). There offshore_hs_m = hs_best: the camera wave
  models, TRAINED ON THIS ADCP's wh_4061, blended with buoys CONVERTED TO
  THIS ADCP's wh_4061 by buoy_transfer.py. So C's wave height is "ADCP
  Hs at 21 m off Marconi" and the ADCP's own Hs needs no conversion. Its
  period (offshore_tp_s = tp_s) is the RAW peak period of NDBC 44008
  (gaps from 44013): an open-ocean buoy peak period, not the ADCP's.
  Therefore:
  * ADCP hours use wh_4061 and wp_peak as measured.
  * Other hours use the hindcast/buoy converted to ADCP-equivalent Hs
    (linear, or buoy_transfer.py's direction-aware form) and ADCP-
    equivalent Tp (raw, scaled, or the ADCP median); for each the form
    whose leave-one-week-out predictions give the smaller SETUP error
    wins, since the setup is what reaches the DEM. With several sources,
    Hs is their inverse-error-variance mean and Tp comes from the source
    whose period gives the smallest setup error. hs_sigma_m / tp_sigma_s
    are those CV RMS values (the measured blend RMS where all are present).
  * The Tp currency difference cannot be removed (44008 had no data in
    the ADCP winter): the report gives the ADCP/WIS peak-period ratio
    (WIS standing in for 44008) and the setup bias in metres it implies.
  * C's STILL-WATER REFERENCE is not this period's either. C was fitted
    with each line at the GNSS-R level, and the GNSS-R footprint is the
    surf zone (gnssr_qc.py; gnssir_reflection_audit.py: reflections
    70-210 m out, the waterline 55-90 m), where breaking waves raise the
    mean level: GNSS-R already holds part of the setup and C only the
    rest. The ADCP (21 m depth) and the Chatham harbour gauge see none.
    The MEAN part of that share is already in the datum chain: GNSS-R
    sits ~0.02 m below Chatham after its +0.349 m datum fix (OPUS
    +/-0.061 m), a comparison of mean levels that holds the GNSS-R's mean
    setup share, and the ADCP datum assumes Marconi = Chatham. So every
    frame here is expected at about +0.02 m - share x (its setup - the
    mean setup of that comparison period): ~0.02 m HIGH in average waves,
    LOW only by the share of the setup above that average. When the 2026
    record is on the computer (the GNSS-R spline, archive/gauge_8447435
    .csv, archive/waves_marconi.csv; --gnssr-spline/--gauge-archive/
    --waves-archive) the report measures that share: (GNSS-R - its
    Chatham transfer) regressed on sqrt(Hs*L0), slope / C, a wave-
    dependent effect (the mean part sits in the transfer's offset and the
    datum chain, and is not seen by the regression). Otherwise it says it
    was not measured. A --setup-coef other than
    0.037 is reported as of unknown currency and reference.

MEASURED ON THE REAL RECORDS (Oct 2024 - Mar 2025 files, Oct 2026):
  linear   a 1.403, lag -72 min, b -0.033: week-out CV RMS 0.166 m (daytime
           0.139), extrapolation RMS 0.174 m (daytime 0.143)
  harmonic M2 ratio 1.48, Marconi 70 min ahead, residual x 0.85: CV RMS
           0.076 m (daytime 0.074), extrapolation RMS 0.085 m (daytime 0.081)
  -> harmonic, sigma_m 0.085 m. Waves: ADCP Tp / WIS Tp 1.06 (r 0.67), so
  setup from the ADCP's own period sits +0.02 m (median) above setup from
  an open-ocean period. WIS ends 2024-12-31, so its conversion rests on
  534 h of December overlap only.

INPUTS ON THE STATION. sig1000_waves_ALL.csv is in the repository and
adcp_water_level_navd88.csv in the Chelsea folder. The Chatham, WIS and
NDBC files are looked for in /mnt/I2Rgus_Data/Chelsea_calibration (then
/home/argus_user, the waterline folder, the repository). If the period
needs Chatham and no file is found (or it does not cover the period), the
6-min NAVD88 record is downloaded from NOAA CO-OPS (compare_gnssr_to_gauge
.fetch_gauge; internet needed, --no-download to forbid) and saved in the
output folder, which forcing.json then names. With the setup on (C > 0) a
period whose waves cover less than half of its daytime hours FAILS (exit 2,
after writing its files): every frame without waves would be left out of
the products after an hour of detection (--min-wave-cover; the October
2024 date needs the WIS ST63064 / NDBC 44013 files, there is no download).

RUN TIME: a few seconds here; ~10-30 s on the station NUC (two cores), most
of it the two transfers' cross-validation and the figure.

Usage:
    python3 historical_forcing.py --start 2024-10-18 --end 2024-10-28 \\
        --output-dir /mnt/I2Rgus_Data/survey_products/2024-10-23/forcing
    (defaults look for the inputs on the station; override any with
     --adcp --adcp-navd88 --chatham --wis --ndbc)
    python3 test_historical_forcing.py        (self-test, ~30 s)

Library:
    from historical_forcing import build_forcing
    info = build_forcing("2025-01-16", "2025-01-24", out_dir)   # forcing.json + file paths
"""

import sys
import json
import time
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Station (NUC) defaults. The 2024-25 files are searched for in these
# folders, in order, so the command runs without options on the station.
CHELSEA = Path("/mnt/I2Rgus_Data/Chelsea_calibration")
SEARCH_DIRS = [CHELSEA, Path("/home/argus_user"), Path("/mnt/I2Rgus_Data/waterline"), HERE]
DEFAULT_NAMES = {             # first existing name in the first folder that has one
    "adcp": ["sig1000_waves_ALL.csv"],
    "adcp_navd88": ["adcp_water_level_navd88.csv"],
    "chatham": ["chatham_2024-10_to_2025-03.csv", "chatham_8447435_navd88.csv"],
    "wis": ["WIS_ST63064_2024-10_to_2025-03.csv"],
    "ndbc": ["NDBC_44013_2024-10_to_2025-03.csv"],
}
CHATHAM_STATION = "8447435"
DOWNLOADED_CHATHAM = "chatham_8447435_navd88.csv"   # written to the output folder

SETUP_COEF = 0.037            # waterline_timex_cron.sh SETUP_COEF
G = 9.81
DAYTIME_UTC = (13.5, 18.0)    # hours whose photos the products use
LAG_GRID_MIN = np.arange(-180, 181, 6)
GAUGE_MAX_GAP_S = 1800.0      # Chatham is 6-minutely; a longer hole is a gap
TARGET_PAD_DAYS = 30.0        # Chatham harmonic fit window around a target period
MIN_FIT_DAYS = 14.0
MIN_AMPLITUDE_M = 0.01        # constituents smaller than this at Chatham give no stable ratio
MSL_TO_NAVD88 = 0.09          # tide-model MSL -> NAVD88 (process_historical.py), datum cross-check
EXTRAPOLATION_TESTS = [       # (fit start, fit end, test start, test end), ends exclusive, UTC
    ("2025-01-01", "2025-03-11", "2024-12-09", "2025-01-01"),
    ("2024-12-09", "2025-02-10", "2025-02-10", "2025-03-11"),
]

# Tidal constituents, speeds in degrees per hour. PRIORITY decides which of
# two unresolvable neighbours is kept.
SPEED = {"M2": 28.9841042, "S2": 30.0, "N2": 28.4397295, "K1": 15.0410686, "O1": 13.9430356,
         "M4": 57.9682084, "MS4": 58.9841042, "MN4": 57.4238337, "Q1": 13.3986609,
         "L2": 29.5284789, "MU2": 27.9682084, "M6": 86.9523127, "2MS6": 87.9682084,
         "MK3": 44.0251729, "M3": 43.4761563}
PRIORITY = ["M2", "S2", "N2", "K1", "O1", "M4", "MS4", "MN4", "Q1", "L2", "MU2", "M6", "2MS6",
            "MK3", "M3"]
T_REF = pd.Timestamp("2000-01-01", tz="UTC")


class ForcingError(Exception):
    """A problem the user must fix (missing file, period outside every record)."""


def say(label, value=""):
    print(f"{label:<18}: {value}" if value != "" else label, flush=True)


# --------------------------------------------------------------------------
# Readers. Every series is (epoch seconds, values), sorted, NaNs dropped.
# --------------------------------------------------------------------------

def to_epoch(times):
    """Seconds since 1970 (UTC; naive stamps are UTC), independent of pandas' time resolution."""
    t = pd.to_datetime(pd.Series(times), utc=True)
    return ((t - pd.Timestamp("1970-01-01", tz="UTC")) // pd.Timedelta(seconds=1)).to_numpy(float)


def iso(ep):
    return pd.to_datetime(np.asarray(ep, float), unit="s", utc=True).strftime("%Y-%m-%dT%H:%M:%SZ")


def day_str(ep):
    return pd.to_datetime(float(ep), unit="s", utc=True).strftime("%Y-%m-%d %H:%M")


def find_input(key, given):
    """The path given, or the first station folder holding the default file name, or None."""
    if given:
        p = Path(given)
        if not p.exists():
            raise ForcingError(f"--{key.replace('_', '-')} {given}: file not found")
        return p
    for d in SEARCH_DIRS:
        for name in DEFAULT_NAMES[key]:
            p = d / name
            if p.exists():
                return p
    return None


def download_chatham(lo_ep, hi_ep, out_dir):
    """The Chatham 6-min NAVD88 record for lo..hi from NOAA CO-OPS, saved as
    out_dir/chatham_8447435_navd88.csv (time, level). Returns the path, or None
    with the reason printed (no internet, API refusal)."""
    from datetime import datetime, timezone
    try:
        from compare_gnssr_to_gauge import fetch_gauge
        lo = datetime.fromtimestamp(lo_ep, tz=timezone.utc)
        hi = datetime.fromtimestamp(hi_ep, tz=timezone.utc)
        say("Chatham download", f"{lo:%Y-%m-%d} .. {hi:%Y-%m-%d} from NOAA CO-OPS (30-day requests, "
                                f"~5 s each) ...")
        g = fetch_gauge(CHATHAM_STATION, lo, hi)
    except (Exception, SystemExit) as exc:      # fetch_gauge exits on an API refusal
        say("Chatham download", f"FAILED ({str(exc).strip()[:200] or type(exc).__name__})")
        return None
    path = Path(out_dir) / DOWNLOADED_CHATHAM
    g.assign(time=g["time"].dt.strftime("%Y-%m-%d %H:%M:%S+00:00")).to_csv(path, index=False)
    say("Chatham download", f"{len(g)} readings -> {path} (copy it to {CHELSEA} to reuse it)")
    return path


def read_adcp(path):
    """sig1000_waves_ALL.csv -> DataFrame(ep, wl, hs, tp, dir), hourly, ADCP's own level zero."""
    a = pd.read_csv(path)
    need = {"time", "water_level", "wh_4061", "wp_peak"}
    if not need <= set(a.columns):
        raise ForcingError(f"{path}: expected columns {sorted(need)}, found {list(a.columns)}")
    out = pd.DataFrame({"ep": to_epoch(a["time"]),
                        "wl": pd.to_numeric(a["water_level"], errors="coerce").to_numpy(float),
                        "hs": pd.to_numeric(a["wh_4061"], errors="coerce").to_numpy(float),
                        "tp": pd.to_numeric(a["wp_peak"], errors="coerce").to_numpy(float),
                        "dir": pd.to_numeric(a.get("wvdir", np.nan), errors="coerce")})
    return out.sort_values("ep").drop_duplicates("ep").reset_index(drop=True)


def read_level_csv(path):
    """A NAVD88 water-level CSV: adcp_water_level_navd88.csv (time, water_level_navd88),
    compare_gnssr_to_gauge.fetch_gauge output (time, level) or a fetch_tide_gauge.py
    archive (epoch, level_navd88)."""
    d = pd.read_csv(path)
    if "epoch" in d and "level_navd88" in d:
        ep, lv = d["epoch"].to_numpy(float), d["level_navd88"]
    elif "time" in d and ("level" in d or "water_level_navd88" in d):
        ep = to_epoch(d["time"])
        lv = d["level"] if "level" in d else d["water_level_navd88"]
    else:
        raise ForcingError(f"{path}: no (time, level), (time, water_level_navd88) or "
                           f"(epoch, level_navd88) columns; found {list(d.columns)}")
    lv = pd.to_numeric(lv, errors="coerce").to_numpy(float)
    ok = np.isfinite(ep) & np.isfinite(lv)
    o = np.argsort(ep[ok], kind="stable")
    ep, lv = ep[ok][o], lv[ok][o]
    keep = np.r_[True, np.diff(ep) > 0]
    return ep[keep], lv[keep]


def read_wis(path):
    """WIS hindcast CSV (datetime, waveHs, waveTp, waveMeanDirection) -> DataFrame(ep, hs, tp, dir)."""
    w = pd.read_csv(path)
    need = {"datetime", "waveHs", "waveTp"}
    if not need <= set(w.columns):
        raise ForcingError(f"{path}: expected WIS columns {sorted(need)}, found {list(w.columns)}")
    out = pd.DataFrame({"ep": to_epoch(w["datetime"]),
                        "hs": pd.to_numeric(w["waveHs"], errors="coerce").to_numpy(float),
                        "tp": pd.to_numeric(w["waveTp"], errors="coerce").to_numpy(float),
                        "dir": pd.to_numeric(w.get("waveMeanDirection", np.nan), errors="coerce")})
    return out.dropna(subset=["hs"]).sort_values("ep").drop_duplicates("ep").reset_index(drop=True)


def read_ndbc(path):
    """NDBC CSV (datetime, WVHT, DPD, APD, MWD; 99/999 = missing) -> hourly DataFrame(ep, hs, tp, dir)."""
    import buoy_transfer
    n = pd.read_csv(path)
    need = {"datetime", "WVHT", "DPD"}
    if not need <= set(n.columns):
        raise ForcingError(f"{path}: expected NDBC columns {sorted(need)}, found {list(n.columns)}")
    hs = pd.to_numeric(n["WVHT"], errors="coerce")
    tp = pd.to_numeric(n["DPD"], errors="coerce")
    dr = pd.to_numeric(n.get("MWD", np.nan), errors="coerce")
    df = pd.DataFrame({"hs": hs.where(hs < 99).to_numpy(), "tp": tp.where(tp < 99).to_numpy(),
                       "dir": dr.where(dr <= 360).to_numpy()},
                      index=pd.to_datetime(n["datetime"], utc=True))
    h = buoy_transfer.hourly(df).dropna(subset=["hs"])
    h[["tp", "dir"]] = h[["tp", "dir"]].ffill(limit=3)      # NDBC often blanks them for an hour
    return pd.DataFrame({"ep": to_epoch(h.index), "hs": h["hs"].to_numpy(float),
                         "tp": h["tp"].to_numpy(float), "dir": h["dir"].to_numpy(float)})


def interp_gapped(ep, val, t, max_gap_s):
    """Linear interpolation; NaN outside the record or where the bracketing samples
    are more than max_gap_s apart (a hole in the record is never bridged)."""
    t = np.asarray(t, float)
    out = np.full(len(t), np.nan)
    if len(ep) < 2:
        return out
    i = np.searchsorted(ep, t)
    inside = (i > 0) & (i < len(ep))
    j = np.clip(i, 0, len(ep) - 1)
    exact = ep[j] == t                      # a reading exactly at t, even at the record's first sample
    ii = np.clip(i, 1, len(ep) - 1)
    ok = inside & ((ep[ii] - ep[ii - 1]) <= max_gap_s)
    out[ok] = np.interp(t[ok], ep, val)
    out[exact] = val[j[exact]]
    return out


def covered(ep, t, max_gap_s):
    """True where both readings bracketing t are within max_gap_s of it (process_chelsea.level_at)."""
    t = np.asarray(t, float)
    if len(ep) == 0:
        return np.zeros(len(t), bool)
    i = np.searchsorted(ep, t)
    j = np.clip(i, 0, len(ep) - 1)
    exact = ep[j] == t
    ii = np.clip(i, 1, len(ep) - 1)
    ok = (i > 0) & (i < len(ep)) & ((t - ep[ii - 1]) <= max_gap_s) & ((ep[ii] - t) <= max_gap_s)
    return ok | exact


def daytime(ep):
    h = (np.asarray(ep, float) % 86400) / 3600.0
    return (h >= DAYTIME_UTC[0]) & (h <= DAYTIME_UTC[1])


def relation(t0, t1, a_ep, max_gap_s):
    """Where [t0, t1) lies relative to the ADCP record, in words (from its coverage, so a
    hole inside the record counts as not measured)."""
    frac = float(covered(a_ep, np.arange(t0, t1, 360.0), max_gap_s).mean())
    if frac >= 1.0:
        return "inside the ADCP record (water level measured)"
    if t1 <= a_ep[0]:
        return f"{(a_ep[0] - t1) / 86400.0:.0f} days BEFORE the ADCP record (transfer extrapolated)"
    if t0 > a_ep[-1]:
        return f"{(t0 - a_ep[-1]) / 86400.0:.0f} days AFTER the ADCP record (transfer extrapolated)"
    return (f"{100 * frac:.0f}% covered by the ADCP (measured there, transferred elsewhere)" if frac > 0
            else "inside a hole in the ADCP record (transfer used)")


def count_daytime(df, t0, t1):
    """Rows per source inside [t0, t1) and 13.5-18 UTC: what the photos actually lean on."""
    e = df["ep"].to_numpy(float)
    m = daytime(e) & (e >= t0) & (e < t1)
    return {str(k): int(v) for k, v in df["source"][m].value_counts().items()}


def stats(res, ep=None):
    """RMS, mean, p95 |r| and n of residuals; with ep also the daytime-only RMS."""
    r = np.asarray(res, float)
    ok = np.isfinite(r)
    out = {"n": int(ok.sum())}
    if ok.sum() == 0:
        return dict(out, rms_m=None, bias_m=None, p95_abs_m=None, rms_daytime_m=None, n_daytime=0)
    out.update(rms_m=round(float(np.sqrt(np.mean(r[ok] ** 2))), 4),
               bias_m=round(float(np.mean(r[ok])), 4),
               p95_abs_m=round(float(np.percentile(np.abs(r[ok]), 95)), 4))
    if ep is not None:
        d = ok & daytime(ep)
        out["n_daytime"] = int(d.sum())
        out["rms_daytime_m"] = (round(float(np.sqrt(np.mean(r[d] ** 2))), 4) if d.sum() else None)
    return out


# --------------------------------------------------------------------------
# Chatham -> Marconi transfers. Both take the full Chatham record as
# `gauge` = (epochs, levels) and predict at any epochs it covers.
# --------------------------------------------------------------------------

class LinearTransfer:
    """Marconi(t) = a * Chatham(t - lag) + b, lag on a 6-min grid, least squares."""
    name = "linear"

    def fit(self, t, y, gauge):
        best = None
        for lag in LAG_GRID_MIN:
            x = interp_gapped(gauge[0], gauge[1], t - lag * 60.0, GAUGE_MAX_GAP_S)
            m = np.isfinite(x) & np.isfinite(y)
            if m.sum() < 48:
                continue
            A = np.c_[x[m], np.ones(m.sum())]
            c, *_ = np.linalg.lstsq(A, y[m], rcond=None)
            rms = float(np.sqrt(np.mean((y[m] - A @ c) ** 2)))
            if best is None or rms < best[0]:
                best = (rms, float(lag), float(c[0]), float(c[1]), int(m.sum()))
        if best is None:
            raise ForcingError("linear transfer: fewer than 48 hours where ADCP and Chatham overlap")
        self.rms_fit, self.lag_min, self.a, self.b, self.n = best
        return self

    def predict(self, t, gauge):
        x = interp_gapped(gauge[0], gauge[1], np.asarray(t, float) - self.lag_min * 60.0, GAUGE_MAX_GAP_S)
        return self.a * x + self.b

    def params(self):
        return {"a": round(self.a, 4), "lag_min": self.lag_min, "b_m": round(self.b, 4),
                "fit_rms_m": round(self.rms_fit, 4), "n_hours": self.n,
                "formula": "Marconi(t) = a * Chatham(t - lag) + b"}

    def describe(self):
        # Marconi(t) = a * Chatham(t - lag): a negative lag reads Chatham LATER, i.e. Marconi leads
        shift = -self.lag_min + 0.0
        who = (f"Marconi leads by {shift:.0f} min" if shift > 0 else
               f"Marconi lags by {-shift:.0f} min" if shift < 0 else "no time shift")
        return (f"Marconi = {self.a:.3f} x Chatham(t {shift:+.0f} min) {self.b:+.3f} m ({who}; "
                f"fit RMS {self.rms_fit:.3f} m, {self.n} h)")


def resolvable(span_hours, names=PRIORITY):
    """Constituents a record of span_hours separates (Rayleigh criterion, factor 1)."""
    keep = []
    for n in names:
        if 360.0 / SPEED[n] > span_hours:
            continue
        if all(abs(SPEED[n] - SPEED[k]) * span_hours >= 360.0 for k in keep):
            keep.append(n)
    return keep


def hours_since_ref(ep):
    return (np.asarray(ep, float) - (T_REF - pd.Timestamp("1970-01-01", tz="UTC")).total_seconds()) / 3600.0


def harmonic_fit(ep, y, names):
    """Least-squares mean + constituents. Returns (mean, {name: complex Z}) with
    signal = Re(Z exp(i w t)), t in hours since 2000-01-01."""
    th = hours_since_ref(ep)
    cols = [np.ones(len(th))]
    for n in names:
        w = np.radians(SPEED[n])
        cols += [np.cos(w * th), np.sin(w * th)]
    c, *_ = np.linalg.lstsq(np.column_stack(cols), y, rcond=None)
    return float(c[0]), {n: complex(c[1 + 2 * k], -c[2 + 2 * k]) for k, n in enumerate(names)}


def harmonic_eval(ep, Z):
    th = hours_since_ref(ep)
    out = np.zeros(len(th))
    for n, z in Z.items():
        out += np.real(z * np.exp(1j * np.radians(SPEED[n]) * th))
    return out


class HarmonicTransfer:
    """Per-constituent amplitude ratio and phase lag + Chatham's residual x beta.

    Marconi(t) = mean_C' + dm + T[tide_C'](t) + beta * resid_C'(t - lag_r) + delta
    where ' marks Chatham's own harmonic fit over the target period +/- 30 days."""
    name = "harmonic"

    def fit(self, t, y, gauge):
        x = interp_gapped(gauge[0], gauge[1], t, GAUGE_MAX_GAP_S)
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() < MIN_FIT_DAYS * 24:
            raise ForcingError(f"harmonic transfer: only {m.sum()} overlapping hours "
                               f"(need {MIN_FIT_DAYS:.0f} days)")
        t, y, x = t[m], y[m], x[m]
        span = (t.max() - t.min()) / 3600.0
        names = resolvable(span)
        mA, ZA = harmonic_fit(t, y, names)
        mC, ZC = harmonic_fit(t, x, names)
        self.names = [n for n in names if abs(ZC[n]) >= MIN_AMPLITUDE_M]
        self.T = {n: ZA[n] / ZC[n] for n in self.names}
        self.dm = mA - mC
        self.amp = {n: (abs(ZA[n]), abs(ZC[n])) for n in self.names}
        # Residuals: Marconi's about its own harmonic fit, Chatham's about its fit
        # with the SAME constituents, on Chatham's 6-min record across the window.
        nA = y - mA - harmonic_eval(t, {n: ZA[n] for n in self.names})
        g = (gauge[0] >= t.min() - 86400) & (gauge[0] <= t.max() + 86400)
        gep = gauge[0][g]
        nC = gauge[1][g] - mC - harmonic_eval(gep, {n: ZC[n] for n in self.names})
        best = None
        for lag in LAG_GRID_MIN:
            z = interp_gapped(gep, nC, t - lag * 60.0, GAUGE_MAX_GAP_S)
            k = np.isfinite(z)
            if k.sum() < 48:
                continue
            A = np.c_[z[k], np.ones(k.sum())]
            c, *_ = np.linalg.lstsq(A, nA[k], rcond=None)
            rms = float(np.sqrt(np.mean((nA[k] - A @ c) ** 2)))
            if best is None or rms < best[0]:
                best = (rms, float(lag), float(c[0]), float(c[1]))
        self.rms_fit, self.lag_r_min, self.beta, self.delta = best
        self.n, self.span_days = int(len(t)), span / 24.0
        self.tide_rms_A = float(np.std(harmonic_eval(t, {n: ZA[n] for n in self.names})))
        self.resid_rms_A = float(np.std(nA))
        return self

    def predict(self, t, gauge):
        t = np.asarray(t, float)
        out = np.full(len(t), np.nan)
        if len(t) == 0:
            return out
        lo, hi = t.min() - TARGET_PAD_DAYS * 86400, t.max() + TARGET_PAD_DAYS * 86400
        g = (gauge[0] >= lo) & (gauge[0] <= hi)
        if g.sum() < 24 * 10 * 5:            # < ~5 days of 6-min data: no stable fit
            return out
        gep, glv = gauge[0][g], gauge[1][g]
        names = [n for n in resolvable((gep.max() - gep.min()) / 3600.0) if n in self.T]
        mC, ZC = harmonic_fit(gep, glv, names)
        nC = glv - mC - harmonic_eval(gep, ZC)
        tide = harmonic_eval(t, {n: self.T[n] * ZC[n] for n in names})
        z = interp_gapped(gep, nC, t - self.lag_r_min * 60.0, GAUGE_MAX_GAP_S)
        # where Chatham itself has a hole, there is no prediction (the residual is unknown)
        hole = ~np.isfinite(interp_gapped(gep, glv, t, GAUGE_MAX_GAP_S))
        out = mC + self.dm + self.delta + tide + self.beta * z
        out[hole] = np.nan
        return out

    def constituent_table(self):
        rows = []
        for n in self.names:
            T = self.T[n]
            lead_min = np.degrees(np.angle(T)) / SPEED[n] * 60.0
            rows.append({"constituent": n, "marconi_amp_m": round(self.amp[n][0], 4),
                         "chatham_amp_m": round(self.amp[n][1], 4),
                         "ratio": round(float(abs(T)), 3), "marconi_leads_min": round(float(lead_min), 1)})
        return rows

    def params(self):
        return {"constituents": self.constituent_table(), "beta": round(self.beta, 4),
                "residual_lag_min": self.lag_r_min, "mean_offset_m": round(self.dm + self.delta, 4),
                "fit_rms_m": round(self.rms_fit, 4), "n_hours": self.n,
                "fit_span_days": round(self.span_days, 1),
                "target_window_pad_days": TARGET_PAD_DAYS,
                "formula": "Marconi = mean_C' + offset + sum_k ratio_k*Chatham_k'(phase+lag_k) "
                           "+ beta*Chatham_residual'(t - residual_lag)"}

    def describe(self):
        m2 = self.T.get("M2")
        m2s = (f"M2 ratio {abs(m2):.3f}, Marconi leads by "
               f"{np.degrees(np.angle(m2)) / SPEED['M2'] * 60:.0f} min; " if m2 is not None else "")
        return (f"{len(self.names)} constituents; {m2s}residual x {self.beta:.2f} at "
                f"{self.lag_r_min:+.0f} min; mean offset {self.dm + self.delta:+.3f} m "
                f"(fit RMS {self.rms_fit:.3f} m, {self.n} h)")


METHODS = {"linear": LinearTransfer, "harmonic": HarmonicTransfer}


def iso_week(ep):
    return pd.to_datetime(np.asarray(ep, float), unit="s", utc=True).strftime("%G-W%V").to_numpy()


def cross_validate(cls, t, y, gauge):
    """Leave one ISO week out; returns predictions at t (NaN where not predicted)."""
    pred = np.full(len(t), np.nan)
    wk = iso_week(t)
    for w in np.unique(wk):
        te = wk == w
        tr = ~te
        if (t[tr].max() - t[tr].min()) / 86400 < MIN_FIT_DAYS:
            continue
        model = cls().fit(t[tr], y[tr], gauge)
        pred[te] = model.predict(t[te], gauge)
    return pred


def extrapolation(cls, t, y, gauge):
    """The two fixed extrapolation tests; returns (predictions at t, list of per-test stats)."""
    pred = np.full(len(t), np.nan)
    tests = []
    for fs, fe, ts, te in EXTRAPOLATION_TESTS:
        f = (t >= to_epoch([fs])[0]) & (t < to_epoch([fe])[0])
        s = (t >= to_epoch([ts])[0]) & (t < to_epoch([te])[0])
        label = f"fit {fs}..{fe} -> predict {ts}..{te}"
        if f.sum() < MIN_FIT_DAYS * 24 or s.sum() < 48:
            tests.append({"test": label, "skipped": "not enough ADCP hours"})
            continue
        model = cls().fit(t[f], y[f], gauge)
        p = model.predict(t[s], gauge)
        pred[s] = p
        tests.append(dict(stats(p - y[s], t[s]), test=label))
    return pred, tests


# --------------------------------------------------------------------------
# Water level
# --------------------------------------------------------------------------

def adcp_on_navd88(adcp, navd88_path, gauge):
    """ADCP level on NAVD88 + a description of its datum. Uses the station's
    adcp_water_level_navd88.csv when given; otherwise adcp_to_navd88.py's rule
    (offset = mean Chatham - mean ADCP over the same hours)."""
    a = adcp.dropna(subset=["wl"])
    x = interp_gapped(gauge[0], gauge[1], a["ep"].to_numpy(), 720.0) if gauge else None
    datum = {}
    if x is not None and np.isfinite(x).sum() > 48:
        m = np.isfinite(x)
        datum["offset_from_chatham_mean_m"] = round(float(np.mean(x[m]) - np.mean(a["wl"].to_numpy()[m])), 4)
    datum["offset_from_msl_constant_m"] = round(MSL_TO_NAVD88 - float(a["wl"].mean()), 4)
    if navd88_path is not None:
        ep, lv = read_level_csv(navd88_path)
        own = interp_gapped(a["ep"].to_numpy(), a["wl"].to_numpy(), ep, 1.0)
        k = np.isfinite(own)
        datum["offset_used_m"] = round(float(np.median(lv[k] - own[k])), 4) if k.any() else None
        datum["source"] = f"{navd88_path} (adcp_to_navd88.py)"
        return ep, lv, datum
    if "offset_from_chatham_mean_m" not in datum:
        raise ForcingError("No adcp_water_level_navd88.csv and no Chatham overlap to derive the ADCP "
                           "datum: give --adcp-navd88 or --chatham")
    off = datum["offset_from_chatham_mean_m"]
    datum["offset_used_m"] = off
    datum["source"] = "derived here as adcp_to_navd88.py does: mean Chatham - mean ADCP, same hours"
    return a["ep"].to_numpy(), a["wl"].to_numpy() + off, datum


def white_noise(y):
    """White-noise level of a smooth series from its second differences (robust)."""
    d2 = y[2:] - 2 * y[1:-1] + y[:-2]
    return float(1.4826 * np.median(np.abs(d2 - np.median(d2))) / np.sqrt(6.0))


def adcp_noise(ep, lv):
    """Random error of one ADCP reading: the white-noise level of its departures from its
    own harmonic fit (the tide's curvature would otherwise dominate the differences)."""
    if len(ep) < 24 * MIN_FIT_DAYS:
        return 0.02
    m, Z = harmonic_fit(ep, lv, resolvable((ep.max() - ep.min()) / 3600.0))
    return white_noise(lv - m - harmonic_eval(ep, Z))


def water_level(t0, t1, a_ep, a_lv, gauge, fit_window, max_gap_s, method_choice, log):
    """Assembles the water-level rows for [t0, t1). Returns (DataFrame, info dict, plot data)."""
    info = {"sources_used": {}, "methods": {}}
    plot = {}
    if gauge is not None:
        ov = np.isfinite(interp_gapped(gauge[0], gauge[1], a_ep, GAUGE_MAX_GAP_S))
    else:
        ov = np.zeros(len(a_ep), bool)
    fw = ov.copy()
    if fit_window[0] is not None:
        fw &= a_ep >= fit_window[0]
    if fit_window[1] is not None:
        fw &= a_ep < fit_window[1]
    chosen = None
    if gauge is not None and fw.sum() >= MIN_FIT_DAYS * 24:
        t, y = a_ep[fw], a_lv[fw]
        info["fit_window"] = [day_str(t.min()), day_str(t.max())]
        say("transfer fit", f"{fw.sum()} ADCP hours, {day_str(t.min())} .. {day_str(t.max())} UTC")
        to, yo = a_ep[ov], a_lv[ov]
        for name, cls in METHODS.items():
            t_start = time.time()
            model = cls().fit(t, y, gauge)
            cvp = cross_validate(cls, t, y, gauge)
            exp_p, tests = extrapolation(cls, to, yo, gauge)
            cv = stats(cvp - y, t)
            pooled = stats(exp_p - yo, to)
            info["methods"][name] = {"params": model.params(), "description": model.describe(),
                                     "cv": cv, "extrapolation_tests": tests,
                                     "extrapolation_pooled": pooled}
            plot.setdefault("cv_resid", {})[name] = (t, cvp - y)
            plot.setdefault("exp_resid", {})[name] = (to, exp_p - yo)
            say(f"  {name}", model.describe())
            say("    CV (week out)", f"RMS {cv['rms_m']:.3f} m (daytime {cv['rms_daytime_m']:.3f}), "
                f"bias {cv['bias_m']:+.3f}, p95 |r| {cv['p95_abs_m']:.3f}")
            for tst in tests:
                if "skipped" in tst:
                    say("    extrapolation", f"{tst['test']}: skipped ({tst['skipped']})")
                else:
                    say("    extrapolation", f"{tst['test']}: RMS {tst['rms_m']:.3f} m (daytime "
                        f"{tst['rms_daytime_m']:.3f}), bias {tst['bias_m']:+.3f}")
            say("    took", f"{time.time() - t_start:.1f} s")
            log.append(name)
            info["methods"][name]["_model"] = model
        scores = {n: (m["extrapolation_pooled"]["rms_m"], m["cv"]["rms_m"])
                  for n, m in info["methods"].items() if m["extrapolation_pooled"]["rms_m"] is not None}
        if method_choice in METHODS:
            chosen = method_choice
            why = "forced with --method"
        elif scores:
            chosen = min(scores, key=lambda n: scores[n])
            why = "lower pooled extrapolation RMS"
            cv_best = min(scores, key=lambda n: scores[n][1])
            if cv_best != chosen:
                why += f" (the leave-one-week-out CV prefers {cv_best})"
        else:
            chosen = min(info["methods"], key=lambda n: info["methods"][n]["cv"]["rms_m"])
            why = "lower CV RMS (extrapolation tests could not run)"
        m = info["methods"][chosen]
        ex = m["extrapolation_pooled"]
        cand = [v for v in (ex.get("rms_m"), ex.get("rms_daytime_m")) if v is not None]
        sigma = max(cand) if cand else m["cv"]["rms_m"]
        info.update(method=chosen, why=why, params=m["params"], cv_rms_m=m["cv"]["rms_m"],
                    cv_rms_daytime_m=m["cv"]["rms_daytime_m"], extrapolation_rms_m=ex.get("rms_m"),
                    extrapolation_rms_daytime_m=ex.get("rms_daytime_m"),
                    extrapolation_bias_m=ex.get("bias_m"), transfer_sigma_m=sigma)
        say("transfer used", f"{chosen} ({why}); sigma {sigma:.3f} m per transferred reading")
    elif gauge is not None:
        say("transfer fit", f"only {fw.sum()} ADCP hours overlap Chatham in the fit window: "
                            f"no transfer fitted")

    # rows: ADCP where it covers, Chatham transfer elsewhere (6-min grid)
    # one reading beyond each end of the period, so a photo near midnight on the last
    # day still has a reading on both sides
    a_sel = (a_ep >= t0 - max_gap_s) & (a_ep <= t1 + max_gap_s)
    noise = adcp_noise(a_ep, a_lv)
    info["adcp_sigma_m"] = round(noise, 4)
    rows = [pd.DataFrame({"ep": a_ep[a_sel], "water_level_navd88": a_lv[a_sel], "source": "adcp",
                          "sigma_m": round(noise, 4)})]
    grid = np.arange(t0, t1 + 1.0, 360.0)
    need = grid[~covered(a_ep, grid, max_gap_s)] if len(a_ep) else grid
    if len(need) and chosen is not None:
        model = info["methods"][chosen]["_model"]
        p = model.predict(need, gauge)
        ok = np.isfinite(p)
        rows.append(pd.DataFrame({"ep": need[ok], "water_level_navd88": p[ok],
                                  "source": "chatham_transfer", "sigma_m": round(info["transfer_sigma_m"], 4)}))
    df = pd.concat(rows, ignore_index=True).sort_values("ep").reset_index(drop=True)
    # times not covered by any row: spans of the 6-min grid with no reading within max_gap_s
    inside = grid[grid < t1]
    nocov = inside[~covered(df["ep"].to_numpy(float), inside, max_gap_s)] if len(df) else inside
    info["uncovered_spans"] = spans(nocov, 360.0)
    for s, n in df["source"].value_counts().items():
        info["sources_used"][s] = int(n)
    plot["rows"] = df
    if chosen is not None:
        model = info["methods"][chosen]["_model"]
        plot["transfer_in_window"] = (grid, model.predict(grid, gauge))
    for n in info["methods"]:
        info["methods"][n].pop("_model", None)
    return df, info, plot


def spans(ep, step):
    """Consecutive runs of a regular grid -> list of 'start .. end' strings (UTC)."""
    if len(ep) == 0:
        return []
    br = np.flatnonzero(np.diff(ep) > step * 1.5)
    starts = np.r_[ep[0], ep[br + 1]]
    ends = np.r_[ep[br], ep[-1]]
    return [f"{day_str(s)} .. {day_str(e)} UTC" for s, e in zip(starts, ends)]


# --------------------------------------------------------------------------
# Waves
# --------------------------------------------------------------------------

def setup_m(hs, tp, coef=SETUP_COEF):
    """Stockdon-form setup C * sqrt(Hs * L0), L0 = g Tp^2 / (2 pi) (extract_elevation_contours.py)."""
    hs, tp = np.asarray(hs, float), np.asarray(tp, float)
    with np.errstate(invalid="ignore"):
        return coef * np.sqrt(hs * G * tp ** 2 / (2 * np.pi))


class HsLinear:
    label = "linear"

    def fit(self, src, adcp_hs):
        self.k, self.c = np.polyfit(src["hs"], adcp_hs, 1)
        return self

    def predict(self, src):
        return np.clip(self.k * src["hs"].to_numpy(float) + self.c, 0.05, None)

    def describe(self):
        return f"Hs_adcp = {self.k:.3f} x Hs {self.c:+.3f} m"


class HsDirectional:
    """buoy_transfer.py's form: log Hs_adcp = c0 + direction harmonics + c5 log Tp + c6 log Hs,
    only in the direction sectors the fit saw (>= 10 h); NaN elsewhere."""
    label = "direction-aware"

    def fit(self, src, adcp_hs):
        import buoy_transfer
        ok = (np.isfinite(src["dir"]) & (src["hs"] > 0.05) & (src["tp"] > 1)).to_numpy()
        X = buoy_transfer.features(src["hs"][ok], src["tp"][ok], src["dir"][ok])
        coef, *_ = np.linalg.lstsq(X, np.log(adcp_hs[ok]), rcond=None)
        d = src["dir"].to_numpy(float)[ok]
        ratio = adcp_hs[ok] / src["hs"].to_numpy(float)[ok]
        sectors = {}
        for lo in range(0, 360, 45):
            m = (d >= lo) & (d < lo + 45)
            if m.sum() >= 10:
                sectors[f"{lo}-{lo + 45}"] = round(float(np.median(ratio[m])), 2)
        self.transfer = {"coef": coef.tolist(), "median_ratio_by_direction": sectors}
        return self

    def predict(self, src):
        import buoy_transfer
        return buoy_transfer.apply(self.transfer, src["hs"], src["tp"], src["dir"])

    def describe(self):
        return ("log-linear in direction, log Tp, log Hs (buoy_transfer.py form); sectors seen: "
                + ", ".join(self.transfer["median_ratio_by_direction"]))


class TpScaled:
    label = "scaled"

    def fit(self, src, adcp_tp):
        self.r = float(np.median(adcp_tp / src["tp"].to_numpy(float)))
        return self

    def predict(self, src):
        return self.r * src["tp"].to_numpy(float)

    def describe(self):
        return f"Tp_adcp = {self.r:.3f} x Tp"


class TpRaw:
    label = "raw"

    def fit(self, src, adcp_tp):
        return self

    def predict(self, src):
        return src["tp"].to_numpy(float)

    def describe(self):
        return "Tp as given"


class TpConstant:
    """The ADCP's median period: what to use when a source's period says little about
    Marconi's (44013 in Massachusetts Bay sees short wind seas the ADCP does not)."""
    label = "adcp-median"

    def fit(self, src, adcp_tp):
        self.t = float(np.median(adcp_tp))
        return self

    def predict(self, src):
        return np.full(len(src), self.t)

    def describe(self):
        return f"Tp = {self.t:.2f} s (ADCP median, ignores the source)"


def cv_wave(cls, src, target, ep):
    """Leave one ISO week out; predictions at every row."""
    pred = np.full(len(src), np.nan)
    wk = iso_week(ep)
    for w in np.unique(wk):
        te = wk == w
        if (~te).sum() < 72:
            continue
        model = cls().fit(src[~te].reset_index(drop=True), target[~te])
        pred[te] = model.predict(src[te].reset_index(drop=True))
    return pred


def align(src, adcp, max_gap_s=1800.0):
    """Source rows at the ADCP hours (nearest within max_gap_s), as (src rows, adcp rows)."""
    if src is None or len(src) == 0:
        return None, None
    i = np.clip(np.searchsorted(src["ep"].to_numpy(), adcp["ep"].to_numpy()), 1, len(src) - 1)
    e = src["ep"].to_numpy()
    j = np.where(np.abs(e[i] - adcp["ep"].to_numpy()) <= np.abs(e[i - 1] - adcp["ep"].to_numpy()), i, i - 1)
    ok = np.abs(e[j] - adcp["ep"].to_numpy()) <= max_gap_s
    ok &= np.isfinite(src["hs"].to_numpy()[j]) & np.isfinite(src["tp"].to_numpy()[j])
    ok &= np.isfinite(adcp["hs"].to_numpy()) & np.isfinite(adcp["tp"].to_numpy())
    s = src.iloc[j[ok]].reset_index(drop=True)
    a = adcp[ok].reset_index(drop=True)
    return s, a


def fit_wave_sources(adcp, sources, coef):
    """For every source with an ADCP overlap: the Hs and Tp models, their week-out CV,
    and the setup error they imply. Returns {source: {...}} with fitted models."""
    out = {}
    for sname, src in sources.items():
        s, a = align(src, adcp)
        if s is None or len(s) < 72:
            out[sname] = {"overlap_hours": 0 if s is None else int(len(s)),
                          "note": "too little overlap with the ADCP to fit"}
            continue
        ep = a["ep"].to_numpy()
        res = {"overlap_hours": int(len(s)), "overlap": [day_str(ep.min()), day_str(ep.max())],
               "raw_hs_rms_m": round(float(np.sqrt(np.mean((s["hs"] - a["hs"]) ** 2))), 3),
               "hs_models": {}, "tp_models": {}}
        hs_cv, tp_cv = {}, {}
        ha, ta = a["hs"].to_numpy(), a["tp"].to_numpy()
        su_a = setup_m(ha, ta, coef)
        # Every model is judged by the SETUP error it alone would cause (the other
        # quantity taken from the ADCP): the setup is what reaches the DEM, and a period
        # bias matters there more than a slightly smaller period RMS.
        for cls in (HsLinear, HsDirectional):
            if cls is HsDirectional and np.isfinite(s["dir"]).sum() < 200:
                continue
            p = cv_wave(cls, s, ha, ep)
            if cls is not HsLinear:
                p = np.where(np.isfinite(p), p, hs_cv["linear"])     # unseen sectors: linear
            hs_cv[cls.label] = p
            model = cls().fit(s, ha)
            st = stats(p - ha)
            st["r"] = round(float(np.corrcoef(p[np.isfinite(p)], ha[np.isfinite(p)])[0, 1]), 3)
            st["setup_rms_m"] = stats(setup_m(p, ta, coef) - su_a)["rms_m"]
            res["hs_models"][cls.label] = dict(st, description=model.describe(), _model=model)
        for cls in (TpRaw, TpScaled, TpConstant):
            p = cv_wave(cls, s, ta, ep)
            tp_cv[cls.label] = p
            model = cls().fit(s, ta)
            st = stats(p - ta)
            st["setup_rms_m"] = stats(setup_m(ha, p, coef) - su_a)["rms_m"]
            res["tp_models"][cls.label] = dict(st, description=model.describe(), _model=model)
        # what the source would do to the setup UNconverted (its own Hs and Tp)
        su_raw = setup_m(s["hs"].to_numpy(float), s["tp"].to_numpy(float), coef)
        res["raw_setup"] = {"bias_median_m": round(float(np.nanmedian(su_raw - su_a)), 3),
                            "rms_m": stats(su_raw - su_a)["rms_m"],
                            "median_m": round(float(np.nanmedian(su_raw)), 3)}
        hs_best = min(res["hs_models"], key=lambda k: res["hs_models"][k]["setup_rms_m"])
        tp_best = min(res["tp_models"], key=lambda k: res["tp_models"][k]["setup_rms_m"])
        hp = hs_cv[hs_best]
        su_p = setup_m(hp, tp_cv[tp_best], coef)
        res.update(hs_best=hs_best, tp_best=tp_best, hs_used_cv=stats(hp - ha),
                   setup_cv=stats(su_p - su_a), setup_median_adcp_m=round(float(np.nanmedian(su_a)), 3),
                   tp_ratio_adcp_over_source=round(float(np.median(a["tp"] / s["tp"])), 3),
                   tp_r=round(float(np.corrcoef(a["tp"], s["tp"])[0, 1]), 3),
                   _cv=(ep, hp, tp_cv[tp_best], ha, ta))
        out[sname] = res
    return out


def combine_sources(fits, coef):
    """How the sources are combined, and what the combination scores on the hours they share.

    Hs: inverse-error-variance mean of every converted source that has the hour (as
    marconi_waves.py blends converted buoys: their errors are only weakly related).
    Tp: the source whose period alone causes the smallest setup error, of those with the hour."""
    names = [n for n in fits if "_cv" in fits[n]]
    res = {"sources": names}
    if not names:
        return res
    sig = {n: fits[n]["hs_used_cv"]["rms_m"] for n in names}
    wsum = sum(1.0 / sig[k] ** 2 for k in names)
    res["hs_sigma_by_source_m"] = sig
    res["hs_weights_when_all_present"] = {n: round(1.0 / sig[n] ** 2 / wsum, 3) for n in names}
    res["tp_order"] = sorted(names, key=lambda n: fits[n]["tp_models"][fits[n]["tp_best"]]["setup_rms_m"])
    if len(names) < 2:
        return res
    common = fits[names[0]]["_cv"][0]
    for n in names[1:]:
        common = np.intersect1d(common, fits[n]["_cv"][0])
    if len(common) < 72:
        return res
    H, T = {}, {}
    for n in names:
        ep, hp, tp, ha, ta = fits[n]["_cv"]
        i = np.searchsorted(ep, common)
        H[n], T[n] = hp[i], tp[i]
        ha_c, ta_c = ha[i], ta[i]
    blend = sum(H[n] / sig[n] ** 2 for n in names) / wsum
    tp = T[res["tp_order"][0]]
    su_a = setup_m(ha_c, ta_c, coef)
    res.update(common_hours=int(len(common)), common_period=[day_str(common.min()), day_str(common.max())],
               hs_cv_on_common_m={n: stats(H[n] - ha_c)["rms_m"] for n in names},
               hs_blend_cv=stats(blend - ha_c),
               setup_cv_on_common_m={n: stats(setup_m(H[n], tp, coef) - su_a)["rms_m"] for n in names},
               setup_blend_cv=stats(setup_m(blend, tp, coef) - su_a))
    return res


def waves(t0, t1, adcp, sources, fits, combo, max_gap_s):
    """Hourly wave rows for [t0, t1): the ADCP where it covers; elsewhere the converted
    sources combined as combine_sources() describes. Returns (rows, uncovered spans)."""
    cols = ["ep", "wvht_m", "dpd_s", "mwd_deg", "source", "hs_sigma_m", "tp_sigma_s"]
    hours = np.arange(np.ceil(t0 / 3600) * 3600, t1 + 1.0, 3600.0)
    a_ok = adcp.dropna(subset=["hs", "tp"]) if adcp is not None else None
    rows = []
    if a_ok is not None and len(a_ok):
        sel = (a_ok["ep"] >= t0 - max_gap_s) & (a_ok["ep"] <= t1 + max_gap_s)
        a = a_ok[sel]
        rows.append(pd.DataFrame({"ep": a["ep"], "wvht_m": a["hs"], "dpd_s": a["tp"], "mwd_deg": a["dir"],
                                  "source": "adcp", "hs_sigma_m": 0.0, "tp_sigma_s": 0.0}))
        need = hours[~covered(a_ok["ep"].to_numpy(), hours, max_gap_s)]
    else:
        need = hours
    names = combo.get("sources", [])
    if len(need) and names:
        H, S, T, TS, D = {}, {}, {}, {}, {}
        for n in names:
            src, f = sources[n], fits[n]
            e = src["ep"].to_numpy()
            i = np.clip(np.searchsorted(e, need), 1, max(len(e) - 1, 1))
            j = np.where(np.abs(e[i] - need) <= np.abs(e[i - 1] - need), i, i - 1)
            near = np.abs(e[j] - need) <= max_gap_s
            sub = src.iloc[j].reset_index(drop=True)
            hm = f["hs_models"][f["hs_best"]]
            hs = hm["_model"].predict(sub)
            if f["hs_best"] != "linear":
                hs = np.where(np.isfinite(hs), hs, f["hs_models"]["linear"]["_model"].predict(sub))
            tp = f["tp_models"][f["tp_best"]]["_model"].predict(sub)
            H[n] = np.where(near, hs, np.nan)
            T[n] = np.where(near, tp, np.nan)
            S[n] = f["hs_used_cv"]["rms_m"]
            TS[n] = f["tp_models"][f["tp_best"]]["rms_m"]
            D[n] = np.where(near, sub["dir"].to_numpy(float), np.nan)
        have = np.column_stack([np.isfinite(H[n]) for n in names])
        w = np.column_stack([np.where(np.isfinite(H[n]), 1.0 / S[n] ** 2, 0.0) for n in names])
        wsum = w.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            hs = np.nansum(np.column_stack([np.nan_to_num(H[n]) for n in names]) * w, axis=1) / wsum
            hs_sig = np.sqrt(1.0 / wsum)            # independent errors (they correlate ~0.2)
        allp = have.all(axis=1)
        if "hs_blend_cv" in combo and len(names) > 1:
            hs_sig[allp] = combo["hs_blend_cv"]["rms_m"]   # measured, where it was measured
        tp = np.full(len(need), np.nan)
        tp_sig = np.full(len(need), np.nan)
        tp_src = np.full(len(need), "", object)
        mwd = np.full(len(need), np.nan)
        for n in reversed(combo["tp_order"]):          # best last, so it wins
            k = np.isfinite(T[n])
            tp[k], tp_sig[k], mwd[k] = T[n][k], TS[n], D[n][k]
            tp_src[k] = n if fits[n]["tp_best"] != TpConstant.label else TpConstant.label
        label = np.array(["hs=" + "+".join(n for n, h in zip(names, row) if h) for row in have], object)
        label = label + ";tp=" + tp_src
        ok = (wsum > 0) & np.isfinite(tp)
        if ok.any():
            rows.append(pd.DataFrame({"ep": need[ok], "wvht_m": hs[ok], "dpd_s": tp[ok], "mwd_deg": mwd[ok],
                                      "source": label[ok], "hs_sigma_m": hs_sig[ok], "tp_sigma_s": tp_sig[ok]}))
        need = need[~ok]
    df = (pd.concat(rows, ignore_index=True).sort_values("ep").reset_index(drop=True)
          if rows else pd.DataFrame(columns=cols))
    return df, spans(need[need < t1], 3600.0)


# --------------------------------------------------------------------------
# Report, figure, driver
# --------------------------------------------------------------------------

SOURCE_NAMES = {"adcp": "ADCP", "wis": "WIS 63064", "ndbc44013": "NDBC 44013", "adcp-median": "ADCP median",
                "chatham_transfer": "Chatham -> Marconi transfer"}


def readable_source(label):
    """'hs=wis+ndbc44013;tp=wis' -> 'Hs: WIS 63064 + NDBC 44013 converted; Tp: WIS 63064'."""
    if label == "adcp":
        return "ADCP (measured at Marconi)"
    if not str(label).startswith("hs="):
        return SOURCE_NAMES.get(label, str(label))
    hs, _, tp = str(label)[3:].partition(";tp=")
    hs = " + ".join(SOURCE_NAMES.get(x, x) for x in hs.split("+"))
    return f"Hs: {hs} converted; Tp: {SOURCE_NAMES.get(tp, tp)}"


def strip_private(d):
    if isinstance(d, dict):
        return {k: strip_private(v) for k, v in d.items() if not str(k).startswith("_")}
    if isinstance(d, list):
        return [strip_private(v) for v in d]
    if isinstance(d, (np.floating,)):
        return float(d)
    if isinstance(d, (np.integer,)):
        return int(d)
    return d


def currency_note(fits, coef, swr=None):
    """The Hs/Tp currency statement, the still-water reference of C, and the setup bias the
    period currency implies."""
    if abs(coef - SETUP_COEF) < 1e-9:
        lines = [
            f"Setup coefficient C = {coef} (waterline_timex_cron.sh, fitted 8 Oct 2026 on 29 Sep - 5 Oct",
            "2026 frames) was fitted with offshore_hs_m = waves_marconi.csv hs_best: camera wave models",
            "trained on THIS ADCP's wh_4061 blended with buoys converted to THIS ADCP's wh_4061",
            "(buoy_transfer.py). Its Hs currency is therefore ADCP Hs at 21 m off Marconi: the ADCP",
            "hours below need no conversion, and every other source is converted to it here.",
            "Its Tp is waves_marconi.csv tp_s = the RAW peak period of NDBC 44008 (open Atlantic).",
            "44008 had no data in the ADCP winter, so the ADCP and 44008 periods cannot be compared",
            "directly. The WIS hindcast (open ocean east of Marconi) stands in for it:",
        ]
    else:
        lines = [
            f"Setup coefficient C = {coef} was given by the caller: where and in which wave currency it",
            f"was fitted is not known here. The statements below assume it was fitted like the live",
            f"C = {SETUP_COEF} (ADCP-currency Hs, NDBC 44008 peak periods, lines at the GNSS-R level).",
            "The WIS hindcast stands in for 44008:",
        ]
    lines_ref = still_water_lines(coef, swr)
    bias = None
    w = fits.get("wis")
    if w and "_cv" in w:
        r = w["tp_ratio_adcp_over_source"]
        lines.append(f"  median ADCP wp_peak / WIS waveTp = {r:.3f} over {w['overlap_hours']} h "
                     f"({w['overlap'][0]} .. {w['overlap'][1]}), r {w['tp_r']:.2f}")
        bias = w.get("_setup_bias")
        if bias is not None:
            lines.append(f"  => setup from the ADCP period is {bias:+.3f} m (median) relative to setup "
                         f"from an open-ocean peak period; historical products are expected to sit")
            lines.append("     that much higher than the live currency would put them (not corrected).")
    else:
        lines.append("  (no WIS overlap given: the Tp currency bias could not be estimated)")
    return lines + [""] + lines_ref, bias


def still_water_lines(coef, swr):
    """C's still-water reference (GNSS-R) against this period's (ADCP / Chatham), in words."""
    L = ["STILL-WATER REFERENCE OF C. C was fitted with each line at the GNSS-R water level",
         "(tide_elevation_navd88 of the live rows). The GNSS-R footprint at Marconi is the surf zone",
         "(gnssr_qc.py; gnssir_reflection_audit.py: reflections 70-210 m from the antenna, the",
         "waterline 55-90 m seaward), where breaking waves raise the mean level: GNSS-R already",
         "contains part of the setup, and C carries only the rest. The ADCP (21 m depth) and the",
         "Chatham harbour gauge see no setup. The MEAN part of the GNSS-R's share is already in",
         "the datum chain: GNSS-R sits ~0.02 m below Chatham after its +0.349 m datum fix (OPUS",
         "+/-0.061 m), a comparison of mean levels that holds that mean share, and the ADCP datum",
         "assumes Marconi = Chatham. Expected, for every frame of this period, against the frame C",
         "was fitted in: about +0.02 m - share x (its setup - the mean setup of the GNSS-R/Chatham",
         "comparison period), i.e. ~0.02 m HIGH in average waves and LOW only by the share of the",
         "setup above that average (the wave-dependent part)."]
    if swr and swr.get("quantified"):
        L.append(f"  Measured on the GNSS-R record ({swr['how']}):")
        L.append(f"    GNSS-R - Chatham transfer = {swr['k']:+.4f} (+/- {swr['k_se']:.4f}) x sqrt(Hs*L0) "
                 f"{swr['c']:+.3f} m, n {swr['n']}, r {swr['r']:.2f}")
        L.append(f"    -> GNSS-R sees {swr['share']:.2f} of the setup (k / C): a wave-dependent effect,")
        L.append(f"       k x (sqrt(Hs*L0) - its 2026 median {swr['typical_x']:.1f} m); at that median the GNSS-R")
        L.append(f"       holds ~{swr['typical_m']:.2f} m of setup, the mean part, which sits in the transfer's offset")
        L.append("       and the datum chain (not seen by this regression, not missing from this period's frames).")
    else:
        L.append("  Not measured here: " + ((swr or {}).get("note") or "no GNSS-R record given") + ".")
        L.append("  On the station: historical_forcing.py measures it from the GNSS-R spline,")
        L.append("  archive/gauge_8447435.csv and archive/waves_marconi.csv when they are there.")
    return L


def setup_share_fit(resid, x):
    """Least squares resid = k*x + c. -> dict(k, k_se, c, n, r) (None if too few)."""
    resid, x = np.asarray(resid, float), np.asarray(x, float)
    ok = np.isfinite(resid) & np.isfinite(x)
    if ok.sum() < 30 or np.ptp(x[ok]) <= 0:
        return None
    A = np.column_stack([x[ok], np.ones(ok.sum())])
    (k, c), *_ = np.linalg.lstsq(A, resid[ok], rcond=None)
    e = resid[ok] - A @ np.array([k, c])
    dof = max(int(ok.sum()) - 2, 1)
    cov = np.linalg.inv(A.T @ A) * float(e @ e) / dof
    return {"k": float(k), "k_se": float(np.sqrt(cov[0, 0])), "c": float(c), "n": int(ok.sum()),
            "r": float(np.corrcoef(x[ok], resid[ok])[0, 1])}


def setup_share_inputs(spline=None, gauge_csv=None, waves_csv=None):
    """The three files gnssr_setup_share() reads, as it resolves them: the GNSS-R spline and the
    Chatham archive (marconi_water_level.py's defaults) and archive/waves_marconi.csv. survey_products.py
    puts their size and time in the forcing step's signature, so a file that arrives or changes
    refreshes the measured share. -> [spline, gauge_csv, waves_csv] as strings (None if unknown)."""
    try:
        from marconi_water_level import GNSSR_SPLINE, GAUGE_CSV
    except Exception:
        GNSSR_SPLINE = GAUGE_CSV = None
    return [str(x) if x else None for x in (spline or GNSSR_SPLINE, gauge_csv or GAUGE_CSV,
                                            waves_csv or HERE / "archive" / "waves_marconi.csv")]


def gnssr_setup_share(coef, spline=None, gauge_csv=None, waves_csv=None):
    """
    How much of the wave setup the GNSS-R level already contains: regress (GNSS-R - Chatham
    transfer) on sqrt(Hs*L0) over the GNSS-R record (marconi_water_level.WaterLevel: the QC'd
    spline and its own Chatham transfer; waves archive/waves_marconi.csv hs_best, tp_s). The
    slope k, divided by C, is the share of the setup the GNSS-R sees (wave-dependent part only:
    the mean part is absorbed by the transfer's offset). -> dict (quantified False + note when
    the files are not on this computer).
    """
    try:
        from marconi_water_level import WaterLevel
    except Exception as exc:
        return {"quantified": False, "note": f"marconi_water_level.py not importable ({exc})"}
    spline, gauge_csv, waves_csv = setup_share_inputs(spline, gauge_csv, waves_csv)
    missing = [str(p) for p in (spline, gauge_csv, waves_csv) if not Path(p).exists()]
    if missing:
        return {"quantified": False, "note": "needs " + ", ".join(missing) + " (not on this computer)",
                "inputs": [str(spline), str(gauge_csv), str(waves_csv)]}
    try:
        wl = WaterLevel(spline=str(spline), gauge_csv=str(gauge_csv))
        if wl.s_ep is None or wl.g_ep is None:
            return {"quantified": False, "note": f"GNSS-R or gauge unreadable ({wl.gps_note})"}
        t = wl.s_ep
        resid = wl.s_lv - wl._chatham(t)
        w = pd.read_csv(waves_csv)
        w = w.dropna(subset=["epoch", "hs_best", "tp_s"]).sort_values("epoch")
        we, hs, tp = w["epoch"].to_numpy(float), w["hs_best"].to_numpy(float), w["tp_s"].to_numpy(float)
        i = np.clip(np.searchsorted(we, t), 1, len(we) - 1)
        j = np.where(np.abs(we[i - 1] - t) < np.abs(we[i] - t), i - 1, i)
        near = np.abs(we[j] - t) <= 3600
        x = np.where(near, np.sqrt(np.clip(hs[j], 0, None) * G * tp[j] ** 2 / (2 * np.pi)), np.nan)
        f = setup_share_fit(resid, x)
    except Exception as exc:                       # a cross-check, never a blocker
        return {"quantified": False, "note": f"failed: {exc}"}
    if not f:
        return {"quantified": False, "note": "too few GNSS-R readings with waves"}
    xt = float(np.nanmedian(x))
    f.update(quantified=True, share=f["k"] / coef if coef else None, typical_x=xt, typical_m=f["k"] * xt,
             how=f"{pd.to_datetime(t[0], unit='s'):%Y-%m-%d} .. {pd.to_datetime(t[-1], unit='s'):%Y-%m-%d}, "
                 f"{wl.describe()}", inputs=[str(spline), str(gauge_csv), str(waves_csv)])
    return f


def write_report(path, args_echo, wl_info, wv_info, fits, wl_df, wv_df, cur_lines):
    L = []
    P = L.append
    P("Historical forcing for the survey products (historical_forcing.py)")
    P("=" * 68)
    for k, v in args_echo.items():
        P(f"{k:<18}: {v}")
    P("")
    P("WATER LEVEL")
    P("-----------")
    d = wl_info.get("adcp_datum", {})
    P(f"ADCP datum        : offset used {d.get('offset_used_m')} m ({d.get('source')})")
    P(f"                    = mean Chatham - mean ADCP over the same hours "
      f"({d.get('offset_from_chatham_mean_m')} m here); the MSL-to-NAVD88 constant (+{MSL_TO_NAVD88} m,")
    P(f"                    process_historical.py) would give {d.get('offset_from_msl_constant_m')} m.")
    P("                    CAVEAT: the ADCP's NAVD88 datum ASSUMES Marconi's mean level equals Chatham's")
    P("                    over the deployment. Any real mean-level difference between the open coast")
    P("                    and the harbour shifts every ADCP and transferred level (and every DEM)")
    P("                    by the same amount; it is NOT in sigma_m.")
    P(f"ADCP sigma_m      : {wl_info.get('adcp_sigma_m')} m (white-noise level of the hourly record)")
    if wl_info.get("methods"):
        P("")
        P("Chatham -> Marconi transfer, fitted on " + " .. ".join(wl_info.get("fit_window", [])) + " UTC")
        P(f"{'method':<10} {'CV RMS':>7} {'day':>6} {'p95':>6} {'bias':>7} | {'extrap RMS':>10} {'day':>6} "
          f"{'p95':>6} {'bias':>7}")
        for n, m in wl_info["methods"].items():
            cv, ex = m["cv"], m["extrapolation_pooled"]
            f = lambda v, s="{:.3f}": ("  n/a" if v is None else s.format(v))
            P(f"{n:<10} {f(cv['rms_m']):>7} {f(cv['rms_daytime_m']):>6} {f(cv['p95_abs_m']):>6} "
              f"{f(cv['bias_m'], '{:+.3f}'):>7} | {f(ex['rms_m']):>10} {f(ex.get('rms_daytime_m')):>6} "
              f"{f(ex['p95_abs_m']):>6} {f(ex['bias_m'], '{:+.3f}'):>7}")
        P("  (m; 'day' = 13.5-18 UTC only; CV = leave one ISO week out; extrap = two tests pooled)")
        for n, m in wl_info["methods"].items():
            P(f"  {n}: {m['description']}")
            for t in m["extrapolation_tests"]:
                if "skipped" in t:
                    P(f"     {t['test']}: skipped ({t['skipped']})")
                else:
                    P(f"     {t['test']}: RMS {t['rms_m']:.3f} (day {t['rms_daytime_m']:.3f}), "
                      f"bias {t['bias_m']:+.3f}, n {t['n']}")
        hm = wl_info["methods"].get("harmonic")
        if hm:
            P("  harmonic constituents (amplitudes over the fit window):")
            P(f"     {'':<5} {'Marconi':>8} {'Chatham':>8} {'ratio':>6} {'Marconi leads':>14}")
            for c in hm["params"]["constituents"]:
                P(f"     {c['constituent']:<5} {c['marconi_amp_m']:>8.3f} {c['chatham_amp_m']:>8.3f} "
                  f"{c['ratio']:>6.3f} {c['marconi_leads_min']:>11.0f} min")
        P(f"USED              : {wl_info.get('method')} -- {wl_info.get('why')}")
        if wl_info.get("period_vs_adcp"):
            P(f"This period       : {wl_info['period_vs_adcp']}")
        P(f"sigma_m transfer  : {wl_info.get('transfer_sigma_m')} m (larger of the pooled extrapolation "
          f"RMS, all hours and daytime)")
        P("NOTE              : the 2026 GNSS-R fit was a 1.24, lag -48 min, b -0.10 (marconi_water_level.py")
        P("                    DEFAULT_FIT): the transfer changes with time (inlet and harbour change),")
        P("                    so the error for a date months from the overlap can exceed the")
        P("                    extrapolation RMS, which is measured 0-9 weeks from its fit window.")
    P("")
    if wl_info.get("chatham_downloaded"):
        P("Chatham record    : downloaded from NOAA CO-OPS for this run (inputs.chatham in forcing.json)")
    P(f"Rows written      : {len(wl_df)}  " + ", ".join(f"{k} {v}" for k, v in wl_info["sources_used"].items()))
    P("Daytime rows      : " + (", ".join(f"{k} {v}" for k, v in wl_info.get("daytime_rows_by_source", {}).items())
                              or "none") + "  (13.5-18 UTC inside the period: what the photos use)")
    if wl_info["uncovered_spans"]:
        P("NOT COVERED       : " + "; ".join(wl_info["uncovered_spans"]))
        P("                    (no ADCP reading within --max-gap-minutes and no Chatham reading: photos")
        P("                    at these times get no water level and are skipped downstream)")
    P("")
    P("WAVES")
    P("-----")
    for line in cur_lines:
        P(line)
    P("")
    for n, f in fits.items():
        if "hs_models" not in f:
            P(f"{n:<6}: {f.get('note')} ({f.get('overlap_hours')} h)")
            continue
        P(f"{n} vs ADCP, {f['overlap_hours']} h ({f['overlap'][0]} .. {f['overlap'][1]}); raw Hs RMS "
          f"{f['raw_hs_rms_m']:.3f} m")
        for k, m in f["hs_models"].items():
            P(f"   Hs {k:<16} CV RMS {m['rms_m']:.3f} m, bias {m['bias_m']:+.3f}, r {m['r']:.2f}; setup error "
              f"from Hs alone {m['setup_rms_m']:.3f} m  [{m['description']}]")
        for k, m in f["tp_models"].items():
            P(f"   Tp {k:<16} CV RMS {m['rms_m']:.2f} s, bias {m['bias_m']:+.2f}; setup error from Tp alone "
              f"{m['setup_rms_m']:.3f} m  [{m['description']}]")
        sc = f["setup_cv"]
        rw = f.get("raw_setup")
        if rw:
            P(f"   unconverted (own Hs and Tp): setup bias {rw['bias_median_m']:+.3f} m (median), RMS "
              f"{rw['rms_m']:.3f} m against the ADCP's setup")
        P(f"   -> uses Hs {f['hs_best']}, Tp {f['tp_best']}: setup (C={wv_info['setup_coef']}) CV RMS "
          f"{sc['rms_m']:.3f} m, bias {sc['bias_m']:+.3f} m (ADCP setup median {f['setup_median_adcp_m']:.3f} m)")
    c = wv_info.get("combination", {})
    if c.get("sources"):
        P("")
        P("Combination for hours without the ADCP:")
        P("   Hs = inverse-error-variance mean of the converted sources present (weights when all are "
          "present: " + ", ".join(f"{k} {v}" for k, v in c["hs_weights_when_all_present"].items()) + ")")
        P("   Tp = first available of " + ", ".join(c["tp_order"]) + " (smallest setup error from Tp first)")
        if "hs_blend_cv" in c:
            P(f"   on the {c['common_hours']} h all sources share ({c['common_period'][0]} .. "
              f"{c['common_period'][1]}), week-out CV:")
            for n in c["sources"]:
                P(f"      {n:<10} Hs RMS {c['hs_cv_on_common_m'][n]:.3f} m   setup RMS "
                  f"{c['setup_cv_on_common_m'][n]:.3f} m")
            P(f"      {'combined':<10} Hs RMS {c['hs_blend_cv']['rms_m']:.3f} m   setup RMS "
              f"{c['setup_blend_cv']['rms_m']:.3f} m, setup bias {c['setup_blend_cv']['bias_m']:+.3f} m")
        P("   hs_sigma_m: the combined CV RMS where every source is present, else sqrt(1/sum(1/sigma^2))")
        P("   of the sources present; tp_sigma_s: the Tp CV RMS of the source used.")
    P("")
    P(f"Rows written      : {len(wv_df)} hours  " + ", ".join(f"{k} {v}" for k, v in wv_info["sources_used"].items()))
    for k, v in wv_info.get("source_labels", {}).items():
        if k != v:
            P(f"   {k} = {v}")
    P("Daytime hours     : " + (", ".join(f"{k} {v}" for k, v in wv_info.get("daytime_rows_by_source", {}).items())
                              or "none") + "  (13.5-18 UTC inside the period)")
    if wv_info["uncovered_spans"]:
        P("NOT COVERED       : " + "; ".join(wv_info["uncovered_spans"]))
        P("                    (no wave height+period: those frames get no setup correction downstream)")
    if wv_info.get("setup_in_window"):
        s = wv_info["setup_in_window"]
        P(f"Setup in window   : median {s['median_m']:.3f} m, range {s['min_m']:.3f} .. {s['max_m']:.3f} m "
          f"(C = {wv_info['setup_coef']}; daytime median {s['daytime_median_m']})")
    for note in wv_info.get("caveats", []):
        P("CAVEAT            : " + note)
    Path(path).write_text("\n".join(L) + "\n")


INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
C_ADCP, C_TRANSFER, C_THIRD = "#2a78d6", "#eb6834", "#1baf7a"


def make_figure(path, t0, t1, wl_plot, wl_df, wv_df, wl_info, wv_info, gauge):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    rc = {"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
          "ytick.color": INK2, "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold",
          "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False}
    dt = lambda ep: pd.to_datetime(np.asarray(ep, float), unit="s", utc=True).tz_convert(None)
    with plt.rc_context(rc):
        fig, axes = plt.subplots(5, 1, figsize=(10, 15), dpi=110)
        fig.patch.set_facecolor("white")

        def shade_day(ax):
            d = np.floor(t0 / 86400) * 86400
            while d < t1:
                ax.axvspan(dt([d + DAYTIME_UTC[0] * 3600])[0], dt([d + DAYTIME_UTC[1] * 3600])[0],
                           color="#f0efec", lw=0, zorder=0)
                d += 86400

        # 1. water level used
        ax = axes[0]
        shade_day(ax)
        if gauge is not None:
            g = (gauge[0] >= t0) & (gauge[0] < t1)
            if g.any():
                ge, gl = gauge[0][g], gauge[1][g]
                gl = np.where(np.r_[False, np.diff(ge) > GAUGE_MAX_GAP_S], np.nan, gl)   # no line across holes
                ax.plot(dt(ge), gl, color="#a3a29d", lw=0.8, label="Chatham gauge as measured")
        for src, col in (("chatham_transfer", C_TRANSFER), ("adcp", C_ADCP)):
            s = wl_df[wl_df["source"] == src]
            if len(s):
                e = s["ep"].to_numpy(float)
                y = s["water_level_navd88"].to_numpy(float)
                brk = np.r_[np.diff(e) > 3 * 3600, False]       # do not draw across gaps
                yy = np.where(np.r_[False, brk[:-1]], np.nan, y)
                lab = ("ADCP (measured)" if src == "adcp"
                       else f"Chatham -> Marconi ({wl_info.get('method')}), +/- sigma")
                ax.plot(dt(e), yy, color=col, lw=1.6, label=lab, zorder=3)
                if src == "chatham_transfer":
                    sg = s["sigma_m"].to_numpy(float)
                    ax.fill_between(dt(e), yy - sg, yy + sg, color=col, alpha=0.18, lw=0, zorder=2)
        tw = wl_plot.get("transfer_in_window")
        if tw is not None and (wl_df["source"] == "adcp").any():
            ax.plot(dt(tw[0]), tw[1], color=C_TRANSFER, lw=1.0, ls="--",
                    label=f"Chatham -> Marconi ({wl_info.get('method')}), for comparison (in-sample)", zorder=2)
        ax.set_ylabel("Water level (m NAVD88)")
        ax.set_title("Still-water level used (shaded: 13.5-18 UTC, the hours whose photos are used)", loc="left")
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + 0.3 * (hi - lo))          # room for the legend above the tides
        ax.legend(loc="upper left", fontsize=8, ncol=3)
        ax.grid(axis="y", color=GRID, lw=0.6)

        # 2. transfer error over the overlap
        ax = axes[1]
        cols = {"linear": "#6d6c68", "harmonic": C_TRANSFER}
        er = wl_plot.get("exp_resid", {})
        if er:
            for n, (t, r) in er.items():
                ok = np.isfinite(r)
                st = wl_info["methods"][n]["extrapolation_pooled"]
                ax.plot(dt(t[ok]), np.where(np.r_[False, np.diff(t[ok]) > 7200], np.nan, r[ok]),
                        color=cols.get(n, C_THIRD), lw=0.9 if n != wl_info.get("method") else 1.3,
                        label=f"{n}: RMS {st['rms_m']:.3f} m, daytime {st['rms_daytime_m']:.3f} m, "
                              f"bias {st['bias_m']:+.3f} m" + ("  (used)" if n == wl_info.get("method") else ""))
            ax.axhline(0, color=INK2, lw=0.8)
            for k, (fs, fe, ts, te) in enumerate(EXTRAPOLATION_TESTS):
                ax.axvline(pd.Timestamp(ts), color="#a3a29d", lw=0.8, ls=":")
                ax.text(pd.Timestamp(ts), 0.02, f" test {k + 1}: fitted on {fs[5:]} .. "
                        f"{(pd.Timestamp(fe) - pd.Timedelta(days=1)):%m-%d}", transform=ax.get_xaxis_transform(),
                        fontsize=7.5, color=INK2, va="bottom")
            lo, hi = ax.get_ylim()
            ax.set_ylim(lo - 0.1 * (hi - lo), hi + 0.25 * (hi - lo))
            ax.legend(loc="upper center", fontsize=8)
            ax.set_title("Transfer error where the ADCP can check it: extrapolation predictions minus ADCP",
                         loc="left")
        else:
            ax.text(0.5, 0.5, "no ADCP/Chatham overlap: transfer not fitted", ha="center", va="center",
                    transform=ax.transAxes, color=INK2)
            ax.set_title("Transfer error", loc="left")
        ax.set_ylabel("Predicted - ADCP (m)")
        ax.grid(axis="y", color=GRID, lw=0.6)

        # 3-5. waves in the window
        srcs = list(dict.fromkeys(wv_df["source"])) if len(wv_df) else []
        palette = [C_ADCP, C_TRANSFER, C_THIRD, "#eda100", "#e87ba4"]
        colour = {s: (C_ADCP if s == "adcp" else palette[1 + [x for x in srcs if x != "adcp"].index(s) % 4])
                  for s in srcs}
        e = wv_df["ep"].to_numpy(float) if len(wv_df) else np.array([])
        su = setup_m(wv_df["wvht_m"], wv_df["dpd_s"], wv_info["setup_coef"]) if len(wv_df) else np.array([])
        # setup sigma from Hs and Tp sigma: d(setup) = setup*(dHs/(2 Hs) + dTp/Tp) in quadrature
        if len(wv_df):
            with np.errstate(invalid="ignore", divide="ignore"):
                ssig = su * np.sqrt((wv_df["hs_sigma_m"].to_numpy(float) / (2 * wv_df["wvht_m"].to_numpy(float))) ** 2
                                    + (wv_df["tp_sigma_s"].to_numpy(float) / wv_df["dpd_s"].to_numpy(float)) ** 2)
        def band_note(col):
            """The title's uncertainty clause: only when a band is drawn, in plain words."""
            if not len(wv_df):
                return ""
            s = wv_df[col].to_numpy(float)
            if np.nanmax(np.where(np.isfinite(s), s, 0)) > 0:
                return " (shaded: \u00b11 sigma of the conversion)"
            return " (measured by the ADCP: no conversion uncertainty)"
        for ax, key, ylab, title in (
                (axes[2], "wvht_m", "Hs (m, ADCP-equivalent)", "Wave height used" + band_note("hs_sigma_m")),
                (axes[3], "dpd_s", "Tp (s, ADCP-equivalent)", "Peak period used" + band_note("tp_sigma_s")),
                (axes[4], "setup", "Setup (m)",
                 f"Wave setup C*sqrt(Hs*L0) added to each frame, C = {wv_info['setup_coef']}"
                 + band_note("hs_sigma_m"))):
            shade_day(ax)
            for s in srcs:
                m = (wv_df["source"] == s).to_numpy()
                y = su[m] if key == "setup" else wv_df[key].to_numpy(float)[m]
                sg = (ssig[m] if key == "setup" else
                      wv_df["hs_sigma_m" if key == "wvht_m" else "tp_sigma_s"].to_numpy(float)[m])
                em = e[m]
                y = np.where(np.r_[False, np.diff(em) > 2 * 3600], np.nan, y)
                ax.plot(dt(em), y, color=colour[s], lw=1.4, label=readable_source(s))
                if np.nanmax(sg) > 0:
                    ax.fill_between(dt(em), y - sg, y + sg, color=colour[s], alpha=0.18, lw=0)
            ax.set_ylabel(ylab)
            ax.set_title(title, loc="left")
            ax.grid(axis="y", color=GRID, lw=0.6)
            if len(srcs) > 1:
                ax.legend(loc="upper left", fontsize=8)
            elif srcs:
                ax.text(0.01, 0.95, f"source: {readable_source(srcs[0])}", transform=ax.transAxes, va="top", color=INK2, fontsize=8)
            if not srcs:
                ax.text(0.5, 0.5, "no wave record covers this period", ha="center", va="center",
                        transform=ax.transAxes, color=INK2)
        for ax in (axes[0], axes[2], axes[3], axes[4]):
            ax.set_xlim(dt([t0])[0], dt([t1])[0])
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        for ax in axes:
            ax.set_xlabel("UTC")
        fig.suptitle(f"Water level and waves at Marconi, {dt([t0])[0]:%Y-%m-%d} .. {dt([t1 - 1])[0]:%Y-%m-%d} (UTC)",
                     x=0.01, y=0.995, ha="left", fontsize=12, fontweight="bold", color=INK)
        import textwrap
        dwl = wl_info.get("daytime_rows_by_source", {})
        dwv = wv_info.get("daytime_rows_by_source", {})
        sub = (("This period is " + wl_info["period_vs_adcp"] + ".  " if wl_info.get("period_vs_adcp") else "")
               + "Daytime (13.5-18 UTC) water-level readings: "
               + (", ".join(f"{readable_source(k)} {v}" for k, v in dwl.items()) or "none")
               + ".  Daytime wave hours: "
               + (", ".join(f"{readable_source(k)} {v}" for k, v in dwv.items()) or "none") + ".")
        fig.text(0.01, 0.978, "\n".join(textwrap.wrap(sub, 150)), ha="left", va="top", fontsize=8.5, color=INK2)
        fig.tight_layout(rect=(0, 0, 1, 0.965))
        fig.savefig(path, facecolor="white")
        plt.close(fig)


def parse_day(text, what):
    """'YYYY-MM-DD' -> epoch of 00:00 UTC that day; ForcingError on anything else."""
    try:
        ts = pd.Timestamp(str(text))
    except (ValueError, TypeError):
        raise ForcingError(f"{what} {text!r} is not a date (expected YYYY-MM-DD, UTC)")
    if ts.tzinfo is not None:
        ts = ts.tz_convert("UTC").tz_localize(None)
    return to_epoch([ts.normalize()])[0]


def wave_cover(wv_df, t0, t1, max_gap_s):
    """(daytime hours of [t0, t1) with a wave row within max_gap_s, daytime hours): the hours the
    photos are taken in, which is what the setup (and so the products) needs waves for."""
    hours = np.arange(np.ceil(t0 / 3600.0) * 3600.0, t1, 3600.0)
    dh = hours[daytime(hours)]
    if not len(dh) or not len(wv_df):
        return 0, int(len(dh))
    e = np.sort(wv_df["ep"].to_numpy(float))
    i = np.clip(np.searchsorted(e, dh), 1, max(len(e) - 1, 1))
    d = np.minimum(np.abs(e[np.clip(i, 0, len(e) - 1)] - dh), np.abs(e[i - 1] - dh))
    return int((d <= max_gap_s).sum()), int(len(dh))


def build_forcing(start, end, out_dir, adcp=None, adcp_navd88=None, chatham=None, wis=None, ndbc=None,
                  fit_start=None, fit_end=None, max_gap_minutes=60.0, method="auto",
                  setup_coef=SETUP_COEF, plot=True, download=True, gnssr_spline=None, gauge_archive=None,
                  waves_archive=None, min_wave_cover=0.5):
    """Writes water_level.csv, waves.csv, forcing_report.txt, forcing.png and forcing.json
    for start..end (whole UTC days, end inclusive) into out_dir. Returns the forcing.json
    content plus 'water_level_csv' and 'waves_csv'. Raises ForcingError on unusable input,
    and (after writing everything, so it can be looked at) when the setup is on (C > 0) and
    the waves cover less than min_wave_cover of the daytime hours: every frame without waves
    would then be left out of the products, after an hour of detection.
    download=False never contacts NOAA for a missing Chatham record."""
    t_all = time.time()
    t0 = parse_day(start, "--start")
    t1 = parse_day(end, "--end") + 86400.0
    if t1 <= t0:
        raise ForcingError(f"--end {end} is before --start {start}")
    if not max_gap_minutes or float(max_gap_minutes) <= 0:
        raise ForcingError(f"--max-gap-minutes must be positive (got {max_gap_minutes})")
    if method not in ("auto",) + tuple(METHODS):
        raise ForcingError(f"--method {method}: expected auto, " + " or ".join(METHODS))
    fw = (parse_day(fit_start, "--fit-start") if fit_start else None,
          parse_day(fit_end, "--fit-end") + 86400.0 if fit_end else None)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    max_gap_s = float(max_gap_minutes) * 60.0
    paths = {k: find_input(k, v) for k, v in
             (("adcp", adcp), ("adcp_navd88", adcp_navd88), ("chatham", chatham), ("wis", wis), ("ndbc", ndbc))}
    say("period", f"{day_str(t0)} .. {day_str(t1)} UTC")
    for k, p in paths.items():
        say(f"  {k}", str(p) if p else "not found (not used)")

    gauge = read_level_csv(paths["chatham"]) if paths["chatham"] else None
    adcp_df = read_adcp(paths["adcp"]) if paths["adcp"] else None
    if gauge is not None:
        say("Chatham", f"{len(gauge[0])} readings, {day_str(gauge[0][0])} .. {day_str(gauge[0][-1])} UTC")
    if adcp_df is not None:
        say("ADCP", f"{len(adcp_df)} hours, {day_str(adcp_df.ep.min())} .. {day_str(adcp_df.ep.max())} UTC")

    # ADCP on NAVD88 (the station's converted file, else derived from Chatham's mean)
    if adcp_df is not None:
        a_ep, a_lv, datum = adcp_on_navd88(adcp_df, paths["adcp_navd88"], gauge)
        say("ADCP datum", f"offset {datum['offset_used_m']:+.4f} m ({datum['source']}); Chatham-mean rule "
            f"{datum.get('offset_from_chatham_mean_m')}, MSL constant rule {datum['offset_from_msl_constant_m']}")
    elif paths["adcp_navd88"]:
        a_ep, a_lv = read_level_csv(paths["adcp_navd88"])
        datum = {"source": f"{paths['adcp_navd88']} (no sig1000 file: waves from the ADCP unavailable)"}
    else:
        a_ep, a_lv, datum = np.array([]), np.array([]), {}

    # Does the period need Chatham beyond what the file (if any) covers? Then
    # download it: the fit needs the ADCP overlap, the harmonic transfer the
    # target +/- TARGET_PAD_DAYS.
    grid = np.arange(t0, t1, 360.0)
    needs_transfer = not covered(a_ep, grid, max_gap_s).all() if len(a_ep) else True
    gauge_short = gauge is None or gauge[0][0] > t0 or gauge[0][-1] < t1 - 360.0
    downloaded = False
    if needs_transfer and gauge_short:
        why = ("no Chatham file found" if gauge is None else
               f"the Chatham file covers only {day_str(gauge[0][0])} .. {day_str(gauge[0][-1])}")
        if download:
            pad = (TARGET_PAD_DAYS + 1) * 86400.0
            lo = min([t0 - pad] + ([a_ep[0]] if len(a_ep) else []))
            hi = max([t1 + pad] + ([a_ep[-1]] if len(a_ep) else []))
            say("Chatham", f"needed ({why}); downloading")
            got = download_chatham(lo, hi, out)
            if got is not None:
                new = read_level_csv(got)
                if gauge is not None:            # keep the file's readings where the download has none
                    keep = ~np.isin(gauge[0], new[0])
                    ep = np.r_[gauge[0][keep], new[0]]
                    lv = np.r_[gauge[1][keep], new[1]]
                    o = np.argsort(ep, kind="stable")
                    new = (ep[o], lv[o])
                gauge, paths["chatham"], downloaded = new, got, True
        else:
            say("Chatham", f"needed ({why}); --no-download given, so not fetched")

    # Is the period inside any record at all?
    rec = []
    if len(a_ep):
        rec.append(("ADCP", a_ep[0], a_ep[-1]))
    if gauge is not None:
        rec.append(("Chatham", gauge[0][0], gauge[0][-1]))
    if not rec:
        raise ForcingError("No water-level record: give --adcp (sig1000_waves_ALL.csv) and/or --chatham "
                           "(Chatham 8447435 6-min NAVD88 CSV with columns time,level)")
    if all(t1 <= lo or t0 > hi for _, lo, hi in rec):
        raise ForcingError(
            f"{day_str(t0)} .. {day_str(t1)} UTC is outside every water-level record ("
            + "; ".join(f"{n} {day_str(lo)} .. {day_str(hi)}" for n, lo, hi in rec)
            + "). Give --chatham with a record that covers it (compare_gnssr_to_gauge.py --save-gauge); "
              "for 2026 dates use the station's GNSS-R (marconi_water_level.py, the live contours).")

    log = []
    say("water level", "fitting the Chatham transfers (two methods, week-out CV and two extrapolation "
                       "tests; ~5-15 s on the station NUC)")
    wl_df, wl_info, wl_plot = water_level(t0, t1, a_ep, a_lv, gauge, fw, max_gap_s, method, log)
    wl_info["adcp_datum"] = datum
    wl_info["chatham_downloaded"] = downloaded
    if len(a_ep):
        wl_info["period_vs_adcp"] = relation(t0, t1, a_ep, max_gap_s)
        say("period vs ADCP", wl_info["period_vs_adcp"])
    if needs_transfer and gauge is None:
        say("WARNING", "the period is not all inside the ADCP record and there is no Chatham record")
    if len(wl_df) == 0 or not (wl_df["ep"] < t1).any() or not (wl_df["ep"] >= t0).any():
        raise ForcingError(f"No water level for {day_str(t0)} .. {day_str(t1)} UTC: "
                           + "; ".join(f"{n} covers {day_str(lo)} .. {day_str(hi)}" for n, lo, hi in rec))
    if wl_info["uncovered_spans"]:
        say("WARNING", "no water level for " + "; ".join(wl_info["uncovered_spans"]))
    wl_info["daytime_rows_by_source"] = count_daytime(wl_df, t0, t1)
    wl_csv = out / "water_level.csv"
    pd.DataFrame({"time": iso(wl_df["ep"]), "water_level_navd88": wl_df["water_level_navd88"].round(4),
                  "source": wl_df["source"], "sigma_m": wl_df["sigma_m"].round(4)}).to_csv(wl_csv, index=False)
    say("water_level.csv", f"{len(wl_df)} rows ("
        + ", ".join(f"{k} {v}" for k, v in wl_info["sources_used"].items()) + f") -> {wl_csv}")

    # Waves
    say("waves", "fitting hindcast/buoy -> ADCP-equivalent Hs and Tp (week-out CV)")
    sources = {}
    if paths["wis"]:
        sources["wis"] = read_wis(paths["wis"])
        say("  WIS", f"{len(sources['wis'])} h, {day_str(sources['wis'].ep.min())} .. "
            f"{day_str(sources['wis'].ep.max())} UTC")
    if paths["ndbc"]:
        sources["ndbc44013"] = read_ndbc(paths["ndbc"])
        say("  NDBC", f"{len(sources['ndbc44013'])} h, {day_str(sources['ndbc44013'].ep.min())} .. "
            f"{day_str(sources['ndbc44013'].ep.max())} UTC")
    adcp_w = adcp_df[["ep", "hs", "tp", "dir"]] if adcp_df is not None else None
    if adcp_w is None and sources:
        say("WARNING", "no ADCP wave record (sig1000_waves_ALL.csv): the hindcast/buoy cannot be "
                       "converted to the setup coefficient's Hs currency, so they are not used")
    fits = fit_wave_sources(adcp_w, sources, setup_coef) if adcp_w is not None else {}
    for n, f in fits.items():
        if "setup_cv" in f:
            say(f"  {n}", f"Hs {f['hs_best']} CV RMS {f['hs_models'][f['hs_best']]['rms_m']:.3f} m, Tp "
                f"{f['tp_best']} CV RMS {f['tp_models'][f['tp_best']]['rms_m']:.2f} s -> setup CV RMS "
                f"{f['setup_cv']['rms_m']:.3f} m, bias {f['setup_cv']['bias_m']:+.3f} m "
                f"({f['overlap_hours']} h overlap)")
    if "wis" in fits and "_cv" in fits["wis"]:
        s, a = align(sources["wis"], adcp_w)
        fits["wis"]["_setup_bias"] = float(np.nanmedian(setup_m(a["hs"], a["tp"], setup_coef)
                                                        - setup_m(a["hs"], s["tp"], setup_coef)))
    combo = combine_sources(fits, setup_coef)
    wv_df, wv_uncov = waves(t0, t1, adcp_w, sources, fits, combo, max_gap_s)
    swr = gnssr_setup_share(setup_coef, gnssr_spline, gauge_archive, waves_archive) if setup_coef else None
    if swr and swr.get("quantified"):
        say("GNSS-R setup", f"GNSS-R sees {swr['share']:.2f} of the setup (k {swr['k']:+.4f}, n {swr['n']}): a "
            f"wave-dependent effect, k x (sqrt(Hs*L0) - {swr['typical_x']:.1f} m); its mean part "
            f"(~{swr['typical_m']:.2f} m) is in the datum chain")
    elif swr:
        say("GNSS-R setup", "share of the setup in the GNSS-R level not measured: " + swr["note"])
    cur_lines, bias = currency_note(fits, setup_coef, swr)
    if "hs_blend_cv" in combo:
        say("  combined", f"Hs {'+'.join(combo['sources'])} CV RMS {combo['hs_blend_cv']['rms_m']:.3f} m, "
            f"Tp from {combo['tp_order'][0]} -> setup CV RMS {combo['setup_blend_cv']['rms_m']:.3f} m, bias "
            f"{combo['setup_blend_cv']['bias_m']:+.3f} m ({combo['common_hours']} h shared)")
    wv_info = {"setup_coef": setup_coef, "sources_used": {k: int(v) for k, v in wv_df["source"].value_counts().items()},
               "daytime_rows_by_source": count_daytime(wv_df, t0, t1),
               "source_labels": {k: readable_source(k) for k in dict.fromkeys(wv_df["source"])},
               "uncovered_spans": wv_uncov, "combination": combo,
               "hs_currency": "ADCP wh_4061 (Signature 1000, 21 m off Marconi) = the currency of hs_best, "
                              "on which C was fitted",
               "tp_currency": "ADCP wp_peak; C was fitted with NDBC 44008 raw peak period (tp_s)",
               "tp_currency_setup_bias_m": None if bias is None else round(bias, 4),
               "still_water_reference": dict(
                   swr or {}, c_fitted_with="GNSS-R (surf-zone footprint: contains part of the setup)",
                   this_period="ADCP at 21 m / Chatham harbour transfer (no setup)",
                   expected="about +0.02 m - share x (setup - the mean setup of the GNSS-R/Chatham comparison "
                            "period): the mean part of the share is in the datum chain (GNSS-R ~0.02 m below "
                            "Chatham in the mean, OPUS +/-0.061 m); only the wave-dependent part shows")}
    used = set(wv_df["source"]) if len(wv_df) else set()
    notes = []
    if any("wis" in u for u in used) and fits.get("wis", {}).get("overlap_hours"):
        notes.append(f"the WIS conversion rests on {fits['wis']['overlap_hours']} h of overlap "
                     f"({fits['wis']['overlap'][0]} .. {fits['wis']['overlap'][1]}; the WIS file ends "
                     f"{day_str(sources['wis'].ep.max())}), so seasonal differences (e.g. October swell "
                     f"directions) are not tested")
    if any(u != "adcp" and "wis" not in u for u in used):
        notes.append("some hours have no hindcast: their Hs is NDBC 44013 converted alone and their Tp "
                     "the best of a weak set (see tp_sigma_s)")
    if notes:
        wv_info["caveats"] = notes
    inside = ((wv_df["ep"] >= t0) & (wv_df["ep"] < t1)).to_numpy() if len(wv_df) else np.zeros(0, bool)
    if inside.any():
        su = setup_m(wv_df["wvht_m"][inside], wv_df["dpd_s"][inside], setup_coef)
        dmask = daytime(wv_df["ep"][inside])
        wv_info["setup_in_window"] = {"median_m": round(float(np.nanmedian(su)), 3),
                                      "min_m": round(float(np.nanmin(su)), 3),
                                      "max_m": round(float(np.nanmax(su)), 3),
                                      "daytime_median_m": (round(float(np.nanmedian(su[dmask])), 3)
                                                           if dmask.any() else None)}
        say("setup", f"C = {setup_coef}: median {wv_info['setup_in_window']['median_m']:.3f} m in the period "
            f"(daytime {wv_info['setup_in_window']['daytime_median_m']})")
    if wv_uncov:
        say("WARNING", "no waves for " + "; ".join(wv_uncov))
    wv_csv = out / "waves.csv"
    pd.DataFrame({"time_utc": iso(wv_df["ep"]), "epoch": wv_df["ep"].astype(np.int64),
                  "wvht_m": np.round(wv_df["wvht_m"].astype(float), 3),
                  "dpd_s": np.round(wv_df["dpd_s"].astype(float), 2),
                  "mwd_deg": np.round(wv_df["mwd_deg"].astype(float), 0), "source": wv_df["source"],
                  "hs_sigma_m": np.round(wv_df["hs_sigma_m"].astype(float), 3),
                  "tp_sigma_s": np.round(wv_df["tp_sigma_s"].astype(float), 2)}).to_csv(wv_csv, index=False)
    say("waves.csv", f"{len(wv_df)} hours (" + ", ".join(f"{k} {v}" for k, v in wv_info["sources_used"].items())
        + f") -> {wv_csv}")

    # report, json, figure
    wv_info["models"] = {n: strip_private(f) for n, f in fits.items()}
    wl_json = {k: v for k, v in wl_info.items()}
    wl_json["inputs"] = {k: (str(paths[k]) if paths[k] else None) for k in ("adcp", "adcp_navd88", "chatham")}
    wv_info["inputs"] = {k: (str(paths[k]) if paths[k] else None) for k in ("adcp", "wis", "ndbc")}
    echo = {"start": start, "end": end, "output": str(out), "max_gap_minutes": max_gap_minutes,
            "method": method, "fit window": f"{fit_start or 'overlap start'} .. {fit_end or 'overlap end'}"}
    echo.update({k: (str(p) if p else "not found") for k, p in paths.items()})
    write_report(out / "forcing_report.txt", echo, wl_info, wv_info, fits, wl_df, wv_df, cur_lines)
    result = strip_private({"period": {"start": start, "end": end, "start_utc": day_str(t0), "end_utc": day_str(t1)},
                            "water_level": wl_json, "waves": wv_info,
                            "software": {"script": "historical_forcing.py", "commit": git_commit()}})
    if plot:
        make_figure(out / "forcing.png", t0, t1, wl_plot, wl_df, wv_df, wl_info, wv_info, gauge)
    (out / "forcing.json").write_text(json.dumps(result, indent=1) + "\n")
    say("report", str(out / "forcing_report.txt"))
    say("done", f"{time.time() - t_all:.0f} s")
    result["water_level_csv"] = str(wl_csv)
    result["waves_csv"] = str(wv_csv)
    n_cov, n_day = wave_cover(wv_df, t0, t1, max_gap_s)
    say("wave cover", f"{n_cov} of {n_day} daytime hours ({DAYTIME_UTC[0]:g}-{DAYTIME_UTC[1]:g} UTC) have waves")
    if setup_coef and min_wave_cover and n_day and n_cov < min_wave_cover * n_day:
        missing = [k for k in ("adcp", "wis", "ndbc") if not paths[k]]
        raise ForcingError(
            f"waves cover only {n_cov} of the {n_day} daytime hours of {start} .. {end}. With the setup on "
            f"(C = {setup_coef}) every frame without a wave height and period is left out, so the products "
            f"would lose {100 * (1 - n_cov / float(n_day)):.0f}% of their photos after the detection. "
            + (f"Wave files not found: {', '.join(missing)} (looked for as "
               + "; ".join(f"{k}: {', '.join(DEFAULT_NAMES[k])}" for k in missing)
               + f" in {', '.join(str(d) for d in SEARCH_DIRS)}; or give --wis / --ndbc / --adcp). "
               if missing else "")
            + f"Or build with --setup-coef 0 (no setup), or lower --min-wave-cover (now {min_wave_cover:g}). "
              f"The files written are in {out} to look at.")
    return result


def git_commit():
    import subprocess
    try:
        r = subprocess.run(["git", "-C", str(HERE), "rev-parse", "--short", "HEAD"],
                           stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True, timeout=10)
        dirty = subprocess.run(["git", "-C", str(HERE), "status", "--porcelain", "--", Path(__file__).name],
                               stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True, timeout=10)
        return r.stdout.strip() + ("+modified" if dirty.stdout.strip() else "") if r.returncode == 0 else None
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--start", required=True, help="first day, YYYY-MM-DD (UTC)")
    ap.add_argument("--end", required=True, help="last day, YYYY-MM-DD (UTC, inclusive)")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--adcp", help="sig1000_waves_ALL.csv (default: the repo's copy)")
    ap.add_argument("--adcp-navd88", help="adcp_water_level_navd88.csv (default: "
                                          f"{CHELSEA}/adcp_water_level_navd88.csv; derived if absent)")
    ap.add_argument("--chatham", help="Chatham 8447435 6-min NAVD88 CSV (time, level; default: looked for "
                                      f"in {CHELSEA}, else downloaded from NOAA when the period needs it)")
    ap.add_argument("--wis", help="WIS ST63064 hourly CSV (datetime, waveHs, waveTp, waveMeanDirection)")
    ap.add_argument("--ndbc", help="NDBC 44013 CSV (datetime, WVHT, DPD, APD, MWD); optional extra wave source")
    ap.add_argument("--fit-start", help="first day of the ADCP overlap used to fit the transfer")
    ap.add_argument("--fit-end", help="last day (inclusive) of the ADCP overlap used to fit the transfer")
    ap.add_argument("--max-gap-minutes", type=float, default=60.0,
                    help="an ADCP reading must lie within this of a time on both sides to be used "
                         "(default 60, as process_chelsea.py)")
    ap.add_argument("--method", choices=["auto", "linear", "harmonic"], default="auto",
                    help="Chatham transfer (default auto: the lower extrapolation RMS)")
    ap.add_argument("--setup-coef", type=float, default=SETUP_COEF,
                    help=f"C used for the setup statistics in the report (default {SETUP_COEF})")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--no-download", action="store_true",
                    help="never download the Chatham record from NOAA (default: only when needed and missing)")
    ap.add_argument("--gnssr-spline", default=None,
                    help="GNSS-R spline for the setup-share check (default marconi_water_level.py's)")
    ap.add_argument("--gauge-archive", default=None,
                    help="Chatham archive for that check (default archive/gauge_8447435.csv)")
    ap.add_argument("--waves-archive", default=None,
                    help="waves for that check (default archive/waves_marconi.csv)")
    ap.add_argument("--min-wave-cover", type=float, default=0.5,
                    help="with --setup-coef > 0, fail (exit 2, after writing the files) when waves cover less "
                         "than this share of the daytime hours (default 0.5; 0 = never)")
    a = ap.parse_args()
    try:
        build_forcing(a.start, a.end, a.output_dir, adcp=a.adcp, adcp_navd88=a.adcp_navd88, chatham=a.chatham,
                      wis=a.wis, ndbc=a.ndbc, fit_start=a.fit_start, fit_end=a.fit_end,
                      max_gap_minutes=a.max_gap_minutes, method=a.method, setup_coef=a.setup_coef,
                      plot=not a.no_plot, download=not a.no_download, gnssr_spline=a.gnssr_spline,
                      gauge_archive=a.gauge_archive, waves_archive=a.waves_archive,
                      min_wave_cover=a.min_wave_cover)
    except ForcingError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Storm Watch: Pre- And Post-Storm Beach Maps For CACO-05 (Marconi Beach)
=========================================================================
A daily job (its own crontab line, opt-in: storm_watch.sh --install-cron)
that reads the forecast for Marconi Beach, and when strong weather is
forecast within the lead time (default 4 days) starts an EVENT:

  1. it fixes a PRE-storm DEM from the calm days before the storm, with
     exactly the cron's 7-day window-DEM method and options
     (dem_from_contours.py, from the saved waterline points), over an
     explicit UTC date window that ends before the storm's onset;
  2. it emails a heads-up, so the team can plan a GPS or drone survey;
  3. it follows the storm on the observed waves (buoy and camera), leaves
     the storm days out of every DEM (wave setup and surge push the
     waterline up), and detects its end;
  4. it builds POST-storm DEMs from post-storm days only: a first one from
     the first 3 usable days, then a 7-day one;
  5. it maps the change, pre against post, with dem_change.py (like for
     like: the pre DEM is rebuilt from the current ground file with the
     post DEM's options), and emails the erosion/accretion map and volumes.
A forecast storm that never arrives is closed as a FALSE ALARM. A storm the
forecast missed is opened from the observed waves, and its pre DEM is built
after the fact from the saved waterline points.

WHAT THE CAMERA DEM COVERS: ONLY THE INTERTIDAL BEACH, between the low and
high waterlines the cameras see (about -0.8 to +1.0 m NAVD88 on a 7-day DEM).
The dune, the bluff and the beach above the high-water line are NOT mapped;
erosion there (the bluff loss of Sep 2024, the dune loss of Jan 2025) is not
measured by this job (a GPS or drone survey is).

FORECAST SOURCES (each fetched every run; a failed source counts as
unavailable for that run: it never triggers and never crashes the job):
  * NWS api.weather.gov: active alerts for STORM_ALERT_ZONES (ANZ254, the
    coastal waters off Marconi, and MAZ022, Barnstable), the gridpoint raw
    data (wind, and waves when the beach cell has them) and the text forecast.
  * Open-Meteo marine (offshore Marconi 41.90N 69.80W and the 44008 position)
    and wind. Free for non-commercial use; data CC BY 4.0, so every email that
    used it says "Weather data by Open-Meteo.com". STORM_SOURCES can drop it.
  * NOAA GFS-Wave point bulletin for 44008 (polar.ncep.noaa.gov), a fallback.
  The parsers were written from the providers' documentation: no real
  response could be fetched where they were written. Capture one of each on
  the NUC (--dry-run --save-raw DIR) and check them (test_storm_watch.py
  --real-captures DIR) before trusting them.
  Every request sends the User-Agent "caco05-storm-watch (USGS coastal camera
  research)"; api.weather.gov requests add NWS_CONTACT when that station
  setting is filled in (empty by default).

TRIGGER (within STORM_LEAD_HOURS, default 96 h; heights OFFSHORE, deep water):
  (a) forecast Hs >= STORM_HS_M (3.5 m) from STORM_DIR_SECTOR (0-160 deg) for
      >= STORM_MIN_HOURS (3 h), on either offshore point (Open-Meteo), else
      the GFS-Wave 44008 bulletin, else the NWS grid waves;
  (b) a strong NWS alert (High Surf, Coastal Flood, Storm Surge, Storm,
      Hurricane Force Wind, Hurricane or Tropical Storm WARNING);
  (c) Hs >= STORM_HS_SUPPORT_M (3.0 m) in the sector AND a support alert
      (Gale Warning, High Surf or Coastal Flood Advisory, ...);
  (d) long-period swell: Hs >= 3.0 m with peak period >= 12 s (mean period
      >= 10 s where only a mean period is forecast);
  (e) only when no wave forecast could be read: NE-E wind of gale force.
OBSERVED WAVES (Marconi currency: hs_best of archive/waves_marconi.csv, else
44013 converted, else 44008 x 0.83): a storm is the run of hours above 1.5 m
(the cron's DEM_MAX_HS) ending at the first of 24 calm hours. A forecast
event is confirmed at >= 2.5 m; a storm the forecast missed opens an event at
>= 3.0 m (2 h >= 2.5 m), or 44008 >= 3.5 m for 2 h. Thresholds and their
basis: docs/VALIDATION.md 2.10 and 4.9, storm_watch_settings.conf.example.

WINDOWS (UTC dates, as dem_from_contours.py filters them):
  pre    the 7 dates before the onset date (onset date left out), extended
         back to 14 if fewer than 3 dates are usable (a camera with >= 4
         frames passing --max-hs);
  storm  the onset date to the last date with storm frames: in no DEM;
  post3  from the first post-storm date to its 3rd usable date;
  post7  the 7 dates from the first post-storm date (cut short at the next
         storm; merged into one COMPOUND event when fewer than 3 usable dates
         lie between them).
GNSS-R runs about 2 days behind, so a date can be built only once the ground
file covers it: the final pre DEM about 2 days after onset, post3 about 5 and
post7 about 9-10 days after the storm ends.

FILES (on the station, BASE = /mnt/I2Rgus_Data/waterline):
  archive/storm_watch/state.json       the events (survives reboots and missed days)
  archive/storm_watch/events/<id>/     one folder per event: pre/, post3/, post7/,
                                       change_post3/, change_post7/, emails/, event.json, log.txt
  archive/storm_watch/forecasts/<date>/  each run's raw responses, gzipped (an audit trail)
  archive/storm_watch/work/            the working copy of the ground rows (removed when idle)
  logs/storm_watch.log                 one line per run
  logs/storm_watch_last_run            the last run's stamp (station_status.py reads it)
  storm_watch_settings.conf            station settings (untracked; see the .example)
It never writes contour_points_*.csv, dem_intertidal_*, archive/dems/,
archive/dems_model/ or logs/timex_cron.log, and changes nothing in the cron.

Usage (on the station; storm_watch.sh runs it from cron at 07:50 local):
    python3 storm_watch.py                     the daily run
    python3 storm_watch.py --dry-run           fetch, evaluate, print; write and send nothing
    python3 storm_watch.py --dry-run --save-raw /tmp/sw_raw   ... and keep the raw responses
    python3 storm_watch.py --status            the events and what each waits for
    python3 storm_watch.py --no-email          build, but only write the email bodies
    python3 storm_watch.py --replay 2026-09-15 2026-10-20 --out-dir /tmp/sw_replay
                                               a past range, day by day (a hindcast; never emails)
    python3 storm_watch.py --climatology --station 44008 --years 2010-2025
                                               events per year under the wave rule (NDBC history)
Exit status: 0 ran (whatever it decided), 1 bad arguments or settings,
2 unexpected error (the traceback is printed; the run line says ERROR).
"""

import os
import re
import io
import sys
import csv
import json
import gzip
import time
import shutil
import socket
import argparse
import subprocess
import traceback
from pathlib import Path
from datetime import datetime, date, timedelta, timezone

import numpy as np

HERE = Path(__file__).resolve().parent
STATION = "CACO-05"
USER_AGENT = "caco05-storm-watch (USGS coastal camera research)"
STATE_VERSION = 1
SUBJECT = "[CACO-05 storm watch]"
INTERTIDAL_NOTE = ("WHAT THE CAMERA DEM COVERS: only the intertidal beach, between the low and high "
                   "waterlines the cameras see (about -0.8 to +1.0 m NAVD88 on a 7-day DEM). The dune, the "
                   "bluff and the beach above the high-water line are not mapped; erosion there is not "
                   "measured (a GPS or drone survey is).")
OM_ATTRIBUTION = "Weather data by Open-Meteo.com (CC BY 4.0)."
MAX_BODY_BYTES = 5 * 1024 * 1024
OUR_FILES = ("storm_watch",)

# --------------------------------------------------------------------------------------------
# time helpers (everything is UTC inside; hours are integer hours since 1970)
# --------------------------------------------------------------------------------------------


def utc(dt):
    """An aware UTC datetime (a naive one is taken as UTC)."""
    return dt.replace(tzinfo=timezone.utc) if dt.tzinfo is None else dt.astimezone(timezone.utc)


def parse_iso(text):
    """ISO 8601 with 'Z' or an offset -> aware UTC datetime; None when empty or unreadable.
    (Python 3.8's fromisoformat does not read a trailing 'Z'.)"""
    if not text or not isinstance(text, str):
        return None
    s = text.strip()
    if s.endswith("Z"):
        s = s[:-1] + "+00:00"
    try:
        return utc(datetime.fromisoformat(s))
    except ValueError:
        return None


def iso_z(dt):
    return utc(dt).strftime("%Y-%m-%dT%H:%M:%SZ") if dt is not None else None


def hour_of(dt):
    return int(utc(dt).timestamp() // 3600)


def hour_dt(h):
    return datetime.fromtimestamp(int(h) * 3600, tz=timezone.utc)


def hour_iso(h):
    return None if h is None else iso_z(hour_dt(h))


def iso_hour(s):
    d = parse_iso(s)
    return None if d is None else hour_of(d)


def hour_date(h):
    return hour_dt(h).strftime("%Y-%m-%d")


def dadd(d, k):
    return (date.fromisoformat(d) + timedelta(days=k)).isoformat()


def ddiff(a, b):
    """Days from date b to date a (a - b)."""
    return (date.fromisoformat(a) - date.fromisoformat(b)).days


def drange(a, b):
    out, d = [], a
    while d <= b:
        out.append(d)
        d = dadd(d, 1)
    return out


def date_hour0(d):
    return hour_of(datetime.fromisoformat(d + "T00:00:00+00:00"))


def eastern(dt):
    """'Tue 17 Nov 09:00 EST' for a UTC datetime (pandas tz_convert: no zoneinfo on Python 3.8)."""
    try:
        import pandas as pd
        t = pd.Timestamp(utc(dt)).tz_convert("America/New_York")
        return t.strftime("%a %d %b %H:%M ") + (t.tzname() or "ET")
    except Exception:
        return ""


def when(dt):
    """'17 Nov 14Z (Tue 17 Nov 09:00 EST)'."""
    if dt is None:
        return "?"
    dt = utc(dt)
    e = eastern(dt)
    return dt.strftime("%d %b %HZ") + (f" ({e})" if e else "")


COMPASS = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]


def compass(deg):
    if deg is None or not np.isfinite(deg):
        return "?"
    return f"{deg:.0f} deg ({COMPASS[int((deg % 360) / 22.5 + 0.5) % 16]})"


def fnum(x, fmt="%.2f"):
    return "?" if x is None or (isinstance(x, float) and not np.isfinite(x)) else fmt % x


# --------------------------------------------------------------------------------------------
# settings: the cron's (read only, never copied) and the station's (an untracked file)
# --------------------------------------------------------------------------------------------

CRON_KEYS = ("DEM_ENABLE", "DEM_CELL", "DEM_MIN_POINTS", "DEM_MAX_SPREAD", "DEM_MAX_HS", "DEM_MAX_DAY_OFFSET",
             "DEM_WINDOW_DAYS", "DEM_DAY_BLOCK_LOD", "GROUND", "LOCK", "GNSSR_START_FILE")
CRON_DEFAULTS = {"DEM_ENABLE": "1", "DEM_CELL": "2.0", "DEM_MIN_POINTS": "3", "DEM_MAX_SPREAD": "0.5",
                 "DEM_MAX_HS": "1.5", "DEM_MAX_DAY_OFFSET": "0.15", "DEM_WINDOW_DAYS": "7", "DEM_DAY_BLOCK_LOD": "",
                 "GROUND": "$BASE/contour_points_ground.csv", "LOCK": "$BASE/.timex_cron.lock",
                 "GNSSR_START_FILE": "$BASE/archive/gnssr_record_start.txt"}
CRON_NUMERIC = {"DEM_CELL": float, "DEM_MIN_POINTS": int, "DEM_MAX_SPREAD": float, "DEM_MAX_HS": float,
                "DEM_MAX_DAY_OFFSET": float, "DEM_WINDOW_DAYS": int}
CRON_MAY_BE_EMPTY = ("DEM_MAX_HS", "DEM_MAX_DAY_OFFSET", "DEM_DAY_BLOCK_LOD")
VAR_RE = re.compile(r"^([A-Z][A-Z0-9_]*)=(.*)$")


def var_value(raw):
    """The value of a shell VAR= line: the quoted text, else up to whitespace or '#'."""
    raw = raw.rstrip("\n")
    if raw[:1] in ("'", '"'):
        q = raw[0]
        end = raw.find(q, 1)
        return raw[1:end] if end > 0 else raw[1:]
    m = re.match(r"[^\s#]*", raw)
    return m.group(0) if m else ""


def var_lines(text):
    """{VAR: [values in order]} of the column-0 VAR= lines of a shell script."""
    out = {}
    for line in text.splitlines():
        m = VAR_RE.match(line)
        if m:
            out.setdefault(m.group(1), []).append(var_value(m.group(2)))
    return out


def cron_settings(base, cron_path=None):
    """The cron's DEM settings, read from its column-0 VAR= lines (first occurrence; each must appear once),
    "$BASE" replaced by `base`. -> (dict, notes). Missing file or key: the coded defaults, with a note."""
    base = str(base)
    path = Path(cron_path) if cron_path else Path(base) / "waterline_timex_cron.sh"
    notes, out = [], {}
    try:
        found = var_lines(path.read_text(errors="replace"))
    except OSError as exc:
        found = None
        notes.append(f"cron settings not read ({path.name}: {exc.strerror or exc}): coded defaults "
                     "2.0/3/0.5/1.5/0.15/7")
    missing = []
    for k in CRON_KEYS:
        if found is not None and k in found:
            v = found[k][0]
            if len(found[k]) > 1:
                notes.append(f"cron {k} is set {len(found[k])} times: the first ({v!r}) is used")
        else:
            v = CRON_DEFAULTS[k]
            if found is not None:
                missing.append(k)
        if k in CRON_NUMERIC and not (v == "" and k in CRON_MAY_BE_EMPTY):
            try:
                CRON_NUMERIC[k](v)
            except ValueError:
                notes.append(f"cron {k}={v!r} is not a number: the default {CRON_DEFAULTS[k]!r} is used")
                v = CRON_DEFAULTS[k]
        out[k] = v.replace("${BASE}", base).replace("$BASE", base)
    if missing:
        notes.append(f"cron settings not read ({', '.join(missing)} missing): coded defaults for those "
                     "(2.0/3/0.5/1.5/0.15/7)")
    return out, notes


def dem_args(cron, has_hs):
    """The cron's window-DEM options, in its order (waterline_timex_cron.sh hs_arg and step 6b), minus
    --last-days and --series-dir: every storm DEM is built with exactly these."""
    a = []
    if has_hs and cron.get("DEM_MAX_HS"):
        a += ["--max-hs", cron["DEM_MAX_HS"]]
    if cron.get("DEM_MAX_DAY_OFFSET"):
        a += ["--max-day-offset", cron["DEM_MAX_DAY_OFFSET"], "--trend-note"]
    if cron.get("DEM_DAY_BLOCK_LOD"):
        a += ["--day-block-se"]
    a += ["--cell", cron["DEM_CELL"], "--min-points", cron["DEM_MIN_POINTS"], "--max-spread", cron["DEM_MAX_SPREAD"]]
    return a


# The station settings: (default, kind, basis). Kinds: f float, i int, s text, sec sector "a-b",
# codes comma list, pt "lat,lon", f0 float or "" (= the cron's value).
SETTINGS = {
    "STORM_LEAD_HOURS": ("96", "i", "lead time of the trigger (GEFS-Wave skill holds to ~day 4)"),
    "STORM_HS_M": ("3.5", "f", "offshore Hs, rule (a): ~Marconi p99 / 0.83; 44008 read 3.8-3.9 m on 23 Sep 2026"),
    "STORM_MIN_HOURS": ("3", "i", "hours at or above, rule (a)/(c)/(d): one 3-hourly step (6 h caught 1 of 5 ADCP storms)"),
    "STORM_HS_SUPPORT_M": ("3.0", "f", "offshore Hs with a support alert, rule (c)"),
    "STORM_SWELL_HS_M": ("3.0", "f", "offshore Hs of long-period swell, rule (d)"),
    "STORM_SWELL_TP_S": ("12", "f", "peak period of rule (d)"),
    "STORM_SWELL_TM_S": ("10", "f", "mean period of rule (d) where a source gives only a mean period (Tm ~0.8 Tp)"),
    "STORM_DIR_SECTOR": ("0-160", "sec", "wave direction FROM, deg true; Marconi faces ~80 deg"),
    "STORM_WIND_DIR_SECTOR": ("0-135", "sec", "wind direction FROM, rule (e)"),
    "STORM_WIND_MS": ("17", "f", "sustained wind, rule (e): gale, 34 kt"),
    "STORM_WIND_MIN_MS": ("12", "f", "sustained wind with gusts, rule (e)"),
    "STORM_WIND_GUST_MS": ("20", "f", "gusts, rule (e): the Sep 2026 forecast had 50-60 mph gusts"),
    "STORM_WIND_MIN_HOURS": ("6", "i", "hours of rule (e)"),
    "STORM_ONSET_HS_M": ("1.8", "f", "offshore Hs that marks the forecast onset (~1.5 / 0.83)"),
    "STORM_CONFIRM_HS": ("2.5", "f", "Marconi Hs that confirms a forecast storm (ADCP winter p98)"),
    "STORM_CONFIRM_BUOY_HS": ("3.0", "f", "raw 44008 Hs in the sector for 2 h that confirms"),
    "STORM_OBS_OPEN_HS": ("3.0", "f", "Marconi peak that opens a storm the forecast missed (with 2 h >= confirm)"),
    "STORM_OBS_OPEN_BUOY_HS": ("3.5", "f", "raw 44008 Hs in the sector for 2 h that opens one"),
    "STORM_OBS_LOOKBACK_DAYS": ("14", "i", "how far back a missed storm is still opened"),
    "STORM_BREAK_HS": ("2.0", "f0", "Marconi peak of a wave event that cuts a post window (survey_products EVENT_HS_M)"),
    "STORM_BREAK_HOURS": ("24", "i0", "rough hours of a wave event that cuts a post window"),
    "STORM_CALM_HS": ("", "f0", "Marconi calm level; empty = the cron's DEM_MAX_HS (1.5)"),
    "STORM_END_CALM_HOURS": ("24", "i", "calm hours that end a storm"),
    "STORM_END_MIN_PRESENT": ("20", "i", "hours of the calm run that must have a value"),
    "STORM_PRE_DAYS": ("", "i0", "pre window dates; empty = the cron's DEM_WINDOW_DAYS (7)"),
    "STORM_PRE_MAX_DAYS": ("14", "i", "how far back a thin pre window is extended"),
    "STORM_POST1_USABLE": ("3", "i", "usable dates of the first post DEM"),
    "STORM_POST_DAYS": ("", "i0", "post window dates; empty = the cron's DEM_WINDOW_DAYS (7)"),
    "STORM_POST_MAX_DAYS": ("14", "i", "how far a thin post window is extended"),
    "STORM_MIN_USABLE_DAYS": ("3", "i", "usable dates a DEM needs (the day test needs 3)"),
    "STORM_MIN_FRAMES_PER_DAY": ("4", "i", "frames passing --max-hs that make a camera's date usable"),
    "STORM_MIN_WAVE_HOURS": ("20", "i", "observed wave hours each date of a window needs before it is built"),
    "STORM_OFFSHORE_TO_MARCONI": ("0.83", "f", "offshore to Marconi, onshore directions (VALIDATION 4.5)"),
    "STORM_SHELTERED_FACTOR": ("0.5", "f", "offshore to Marconi outside the sector (buoy_transfer.py)"),
    "STORM_ALERT_ZONES": ("ANZ254,MAZ022", "codes", "NWS zones: coastal waters off Marconi, Barnstable"),
    "STORM_ALERT_STRONG": ("SU.W,CF.W,SS.W,SR.W,HF.W,HU.W,TR.W", "codes", "alerts that trigger on their own"),
    "STORM_ALERT_SUPPORT": ("SU.Y,CF.Y,CF.A,GL.W,SR.A,HF.A,SE.W,HW.W,TR.A,HU.A,SS.A", "codes",
                            "alerts that trigger with Hs >= STORM_HS_SUPPORT_M"),
    "STORM_ALERT_WAVE_GATED": ("", "codes", "strong alerts that count only as support when a wave forecast exists "
                                            "(e.g. CF.W,SS.W: MAZ022 also covers the bay shore)"),
    "STORM_SOURCES": ("om_marine,gfswave,nws,om_wind", "codes", "forecast sources (drop om_* to stop using Open-Meteo)"),
    "STORM_POINT": ("41.89,-69.96", "pt", "NWS /points and Open-Meteo wind"),
    "STORM_MARINE_POINT": ("41.90,-69.80", "pt", "offshore Marconi, Open-Meteo marine"),
    "STORM_BUOY_POINT": ("40.50,-69.25", "pt", "44008, Open-Meteo marine"),
    "STORM_MIN_CELL_LON": ("-69.93", "f", "an Open-Meteo cell west of this is in Cape Cod Bay or behind the shore"),
    "NWS_CONTACT": ("", "s", "added to the User-Agent of api.weather.gov requests only; empty by default"),
    "STORM_EMAIL": ("1", "i", "0: email bodies are written to the event folder, nothing is sent"),
    "STORM_EMAIL_TO": ("", "s", "space-separated; empty = owg.sh's RECIPIENTS"),
    "STORM_HTTP_TIMEOUT_S": ("45", "f", "per request"),
    "STORM_FETCH_DEADLINE_S": ("300", "f", "all requests of a run"),
    "STORM_LOCK_WAIT_S": ("1800", "f", "wait for the cron's lock"),
    "STORM_LOCK_HOLD_S": ("900", "f", "the cron's lock is not taken when a cron run starts within this"),
    "STORM_OWG_LOCK_WAIT_S": ("600", "f", "wait for owg.sh's lock"),
    "STORM_MAX_WAIT_DAYS": ("21", "i", "a date the ground file never covers: give up after this"),
    "STORM_BUOY_STALE_H": ("6", "f", "fetch 44008 when its archive is older"),
    "STORM_FORECAST_KEEP_DAYS": ("730", "i", "days of archived forecasts kept"),
}


def parse_sector(text):
    a, b = (float(x) for x in str(text).split("-"))
    if not (0 <= a <= 360 and 0 <= b <= 360):
        raise ValueError(text)
    return a, b


def in_sector(deg, sector):
    """Direction FROM inside 'a-b' (wrap allowed: 340-160); vectorised; NaN -> False."""
    a, b = sector
    d = np.asarray(deg, float) % 360
    with np.errstate(invalid="ignore"):
        return ((d >= a) & (d <= b)) if a <= b else ((d >= a) | (d <= b))


def load_settings(path, cron=None):
    """Station settings (KEY=VALUE, '#' comments) over the defaults. -> (cfg, notes). Unknown keys and bad
    values give a note and the default. cfg holds parsed values; "" keys resolved from the cron."""
    cron = cron or CRON_DEFAULTS
    raw, notes = {}, []
    p = Path(path) if path else None
    if p is not None and p.exists():
        for k, vals in var_lines(p.read_text(errors="replace")).items():
            if k not in SETTINGS:
                notes.append(f"setting {k} in {p.name} is not known: ignored")
                continue
            raw[k] = vals[-1]
    cfg = {}
    for k, (default, kind, _) in SETTINGS.items():
        v = raw.get(k, default)
        try:
            cfg[k] = parse_setting(v, kind)
        except (ValueError, TypeError):
            notes.append(f"setting {k}={v!r} is not valid ({kind}): the default {default!r} is used")
            cfg[k] = parse_setting(default, kind)
    calm_from_cron = cron.get("DEM_MAX_HS") or "1.5"
    if cfg["STORM_CALM_HS"] is None:
        cfg["STORM_CALM_HS"] = float(calm_from_cron)
    win = int(cron.get("DEM_WINDOW_DAYS") or 7)
    for k in ("STORM_PRE_DAYS", "STORM_POST_DAYS"):
        if cfg[k] is None:
            cfg[k] = win
    cfg["STORM_PRE_MAX_DAYS"] = max(cfg["STORM_PRE_MAX_DAYS"], cfg["STORM_PRE_DAYS"])
    cfg["STORM_POST_MAX_DAYS"] = max(cfg["STORM_POST_MAX_DAYS"], cfg["STORM_POST_DAYS"])
    return cfg, notes


def parse_setting(v, kind):
    v = str(v).strip()
    if kind == "f":
        return float(v)
    if kind == "i":
        return int(v)
    if kind == "f0":
        return None if v == "" else float(v)
    if kind == "i0":
        return None if v == "" else int(v)
    if kind == "sec":
        return parse_sector(v)
    if kind == "codes":
        return [x.strip().upper() if "." in x or x.strip().isupper() else x.strip()
                for x in v.split(",") if x.strip()]
    if kind == "pt":
        lat, lon = (float(x) for x in v.split(","))
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            raise ValueError(v)
        return (lat, lon)
    return v


def owg_recipients(base):
    """The addresses of owg.sh's RECIPIENTS line (never logged; only counted)."""
    try:
        for line in (Path(base) / "owg.sh").read_text(errors="replace").splitlines():
            m = re.match(r'^RECIPIENTS="([^"]*)"', line)
            if m:
                return m.group(1).split()
    except OSError:
        pass
    return []


def recipients(base, cfg):
    """-> (list, where from). Never print the addresses anywhere but --dry-run's terminal."""
    if cfg.get("STORM_EMAIL_TO"):
        return cfg["STORM_EMAIL_TO"].split(), "STORM_EMAIL_TO"
    return owg_recipients(base), "owg.sh RECIPIENTS"

# --------------------------------------------------------------------------------------------
# HTTP (stdlib urllib only, like the rest of the repository)
# --------------------------------------------------------------------------------------------


def headers_for(url, cfg):
    """The request headers: the generic User-Agent everywhere; for api.weather.gov the station's NWS_CONTACT
    (when set) and Accept geo+json. Nothing else identifies anyone."""
    ua = USER_AGENT
    h = {"User-Agent": ua}
    if "://api.weather.gov/" in url:
        contact = (cfg or {}).get("NWS_CONTACT") or ""
        if contact:
            h["User-Agent"] = f"caco05-storm-watch (USGS coastal camera research; {contact})"
        h["Accept"] = "application/geo+json"
    return h


def urllib_opener(url, headers, timeout):
    """-> (status, bytes). An HTTPError is returned as its status and body; network errors raise OSError."""
    import urllib.request
    import urllib.error
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.read(MAX_BODY_BYTES + 1)
    except urllib.error.HTTPError as e:
        try:
            body = e.read(65536)
        except Exception:
            body = b""
        return e.code, body


def http_get(url, headers, timeout, opener=None, tries=2, wait_s=8, sleep=time.sleep, deadline=None):
    """-> (status, body, why). Status 0 = no answer (network error or timeout). Retries once, after wait_s,
    on no answer, 5xx or 429; never on another 4xx. A body over 5 MB counts as unavailable."""
    opener = opener or urllib_opener
    status, body, why = 0, b"", ""
    for attempt in range(1, max(tries, 1) + 1):
        if deadline is not None and time.time() > deadline:
            return 0, b"", "fetch deadline passed"
        try:
            status, body = opener(url, headers, timeout)
            body = body or b""
            why = "" if status == 200 else f"HTTP {status}" + problem_title(body)
        except (OSError, ValueError) as exc:
            status, body, why = 0, b"", f"no answer ({type(exc).__name__}: {str(exc)[:80]})"
        if status == 200 and len(body) > MAX_BODY_BYTES:
            return 0, b"", "response over 5 MB"
        if status == 200 or (400 <= status < 500 and status != 429):
            return status, body, why
        if attempt < tries:
            sleep(wait_s)
    return status, body, why


def problem_title(body):
    """': <title>' of an application/problem+json body, for the log."""
    try:
        js = json.loads(body.decode("utf-8", "replace"))
        t = js.get("title") or js.get("reason") or ""
        return f": {str(t)[:80]}" if t else ""
    except (ValueError, AttributeError):
        return ""


# --------------------------------------------------------------------------------------------
# forecast parsers: each returns series on an hourly UTC index, or the reason it is unusable
# A series is a dict: h (int hours, contiguous), named float arrays, and labels.
# --------------------------------------------------------------------------------------------

DUR_RE = re.compile(r"^P(?:(\d+)D)?(?:T(?:(\d+)H)?(?:(\d+)M)?)?$")
UOM = {"wmoUnit:km_h-1": ("speed", 1 / 3.6), "wmoUnit:m_s-1": ("speed", 1.0), "wmoUnit:kt": ("speed", 0.514444),
       "wmoUnit:m": ("length", 1.0), "wmoUnit:ft": ("length", 0.3048), "wmoUnit:s": ("time", 1.0),
       "nwsUnit:s": ("time", 1.0), "wmoUnit:degree_(angle)": ("angle", 1.0)}
NWS_LAYERS = {"windSpeed": ("ws", "speed"), "windGust": ("gust", "speed"), "windDirection": ("wdir", "angle"),
              "waveHeight": ("hs", "length"), "wavePeriod": ("tp", "time"), "waveDirection": ("dir", "angle")}
BOUNDS = {"hs": (0.0, 20.0), "tp": (1.0, 30.0), "dir": (0.0, 360.0), "wdir": (0.0, 360.0), "ws": (0.0, 90.0),
          "gust": (0.0, 90.0)}
MIN_SOURCE_HOURS = 24


def parse_duration(text):
    """ISO 8601 duration of the NWS validTime (PT6H, P1D, P1DT12H, P7DT21H) -> hours (minutes round up)."""
    m = DUR_RE.match(text or "")
    if not m or not any(m.groups()):
        raise ValueError(f"duration {text!r}")
    d, h, mi = (int(x) if x else 0 for x in m.groups())
    return d * 24 + h + (1 if mi else 0)


def expand_layer(layer, kind):
    """An NWS gridpoint layer {uom, values:[{validTime:'start/dur', value}]} -> ({hour: value}, why).
    Each value fills every hour of [start, start + duration); null values are missing; the unit comes
    from uom (an unknown or wrong-kind unit makes the layer unusable)."""
    uom = (layer or {}).get("uom", "")
    if uom not in UOM or UOM[uom][0] != kind:
        return {}, f"unit {uom or 'missing'} not understood"
    fac = UOM[uom][1]
    out = {}
    for v in layer.get("values") or []:
        try:
            start, dur = str(v["validTime"]).split("/")
            h0 = hour_of(parse_iso(start))
            n = parse_duration(dur)
        except (KeyError, ValueError, TypeError, AttributeError):
            continue
        val = v.get("value")
        if val is None:
            continue
        for k in range(n):
            out[h0 + k] = float(val) * fac
    return out, ""


def make_series(hours_values, columns, hold=0):
    """{column: {hour: value}} -> series dict on a contiguous hourly index; each sample held `hold` more hours."""
    keys = [k for c in columns for k in hours_values.get(c, {})]
    if not keys:
        return None
    h0, h1 = min(keys), max(keys) + hold
    h = np.arange(h0, h1 + 1)
    s = {"h": h}
    for c in columns:
        a = np.full(len(h), np.nan)
        for k, v in sorted(hours_values.get(c, {}).items()):
            if v is None or not np.isfinite(v):
                continue
            a[k - h0] = v
        if hold:
            last = np.nan
            age = hold + 1
            for i in range(len(a)):
                if np.isfinite(a[i]):
                    last, age = a[i], 0
                else:
                    age += 1
                    if age <= hold:
                        a[i] = last
        s[c] = a
    return s


def check_series(s, cols, now, lead_h):
    """Plausibility gate: every value inside its bounds, and >= 24 hours with a value inside the lead window.
    -> '' or why the source is unusable this run."""
    if s is None:
        return "no values"
    h0 = hour_of(now)
    w = (s["h"] >= h0) & (s["h"] <= h0 + lead_h)
    for c in cols:
        a = s.get(c)
        if a is None:
            continue
        lo, hi = BOUNDS[c]
        with np.errstate(invalid="ignore"):
            bad = np.isfinite(a) & ((a < lo) | (a > hi))
        if bad.any():
            return f"implausible {c} ({np.nanmax(np.abs(a[bad])):.1f} outside {lo:g}-{hi:g})"
    key = cols[0]
    if int(np.isfinite(s[key][w]).sum()) < MIN_SOURCE_HOURS:
        return f"only {int(np.isfinite(s[key][w]).sum())} forecast hours with {key} in the next {lead_h} h"
    return ""


def parse_nws_points(js):
    p = (js or {}).get("properties") or {}
    if not p.get("forecastGridData"):
        raise ValueError("no forecastGridData")
    return {k: p.get(k) for k in ("forecastGridData", "forecast", "gridId", "gridX", "gridY", "forecastZone", "county")}


def parse_nws_grid(js, now, lead_h):
    """-> {'updated', 'stale', 'wind': series|None, 'waves': series|None, 'why_wind', 'why_waves', 'notes'}"""
    p = (js or {}).get("properties")
    if not isinstance(p, dict):
        raise ValueError("no properties")
    upd = parse_iso(p.get("updateTime"))
    out = {"updated": iso_z(upd) if upd else None, "notes": [],
           "stale": upd is None or (utc(now) - upd).total_seconds() > 24 * 3600}
    vals = {}
    for name, (col, kind) in NWS_LAYERS.items():
        if name not in p:
            continue
        d, why = expand_layer(p[name], kind)
        if why:
            out["notes"].append(f"{name}: {why}")
            continue
        vals[col] = d
    out["wind"] = make_series(vals, ["ws", "gust", "wdir"]) if vals.get("ws") else None
    out["waves"] = make_series(vals, ["hs", "tp", "dir"]) if vals.get("hs") else None
    out["why_wind"] = "" if out["wind"] is not None else "no windSpeed values"
    out["why_waves"] = "" if out["waves"] is not None else "no waveHeight values at this grid cell"
    for k, cols in (("wind", ["ws", "gust", "wdir"]), ("waves", ["hs", "tp", "dir"])):
        if out[k] is not None:
            why = check_series(out[k], cols, now, lead_h)
            if why:
                out[k], out["why_" + k] = None, why
    return out


VTEC_RE = re.compile(r"/([OTEX])\.([A-Z]{3})\.([A-Z]{4})\.([A-Z]{2})\.([WAYS])\.(\d{4})\.(\d{6}T\d{4}Z)-(\d{6}T\d{4}Z)/")
EVENT_CODES = {"high surf warning": "SU.W", "high surf advisory": "SU.Y", "coastal flood warning": "CF.W",
               "coastal flood watch": "CF.A", "coastal flood advisory": "CF.Y", "gale warning": "GL.W",
               "storm warning": "SR.W", "storm watch": "SR.A", "hurricane force wind warning": "HF.W",
               "hurricane force wind watch": "HF.A", "hurricane warning": "HU.W", "hurricane watch": "HU.A",
               "tropical storm warning": "TR.W", "tropical storm watch": "TR.A", "storm surge warning": "SS.W",
               "storm surge watch": "SS.A", "hazardous seas warning": "SE.W", "high wind warning": "HW.W"}


def vtec_time(s):
    if s.startswith("000000"):
        return None
    return datetime.strptime(s, "%y%m%dT%H%MZ").replace(tzinfo=timezone.utc)


def parse_alerts(js, zones):
    """NWS /alerts/active GeoJSON -> [alert]: status Actual, not a Cancel, not a VTEC CAN/EXP/UPG, in the zones
    (an alert without UGC is kept: the query was by zone). Code from the VTEC (PP.S), else from the event name."""
    feats = (js or {}).get("features")
    if not isinstance(feats, list):
        raise ValueError("no features")
    out = []
    for f in feats:
        p = (f or {}).get("properties") or {}
        if p.get("status") != "Actual" or p.get("messageType") == "Cancel":
            continue
        ugc = ((p.get("geocode") or {}).get("UGC")) or []
        if ugc and zones and not set(ugc) & set(zones):
            continue
        params = p.get("parameters") or {}
        vt = None
        for v in params.get("VTEC") or []:
            m = VTEC_RE.search(v)
            if m:
                vt = m
                break
        if vt and vt.group(2) in ("CAN", "EXP", "UPG"):
            continue
        code = f"{vt.group(4)}.{vt.group(5)}" if vt else EVENT_CODES.get(str(p.get("event", "")).strip().lower())
        start = parse_iso(p.get("onset")) or parse_iso(p.get("effective"))
        end = parse_iso(p.get("ends")) or (vtec_time(vt.group(8)) if vt else None) or parse_iso(p.get("expires"))
        if start is None:
            start = parse_iso(p.get("sent"))
        if start is None or end is None:
            continue
        out.append({"code": code, "event": p.get("event", ""), "zones": ugc, "start": iso_z(start), "end": iso_z(end),
                    "headline": p.get("headline") or "", "nws_headline": " ".join(params.get("NWSheadline") or []),
                    "vtec_action": vt.group(2) if vt else None})
    return out


def om_objects(js):
    if isinstance(js, dict) and js.get("error"):
        raise ValueError(str(js.get("reason", "error"))[:80])
    return js if isinstance(js, list) else [js]


def om_times(hourly):
    return [hour_of(datetime.strptime(t, "%Y-%m-%dT%H:%M").replace(tzinfo=timezone.utc)) for t in hourly["time"]]


def om_fresh(hours, now):
    """The first time on today's UTC date or earlier, the last >= now + 72 h."""
    if not hours:
        return "no times"
    today0 = date_hour0(utc(now).strftime("%Y-%m-%d"))
    if min(hours) > today0 + 23:
        return "starts after today"
    if max(hours) < hour_of(now) + 72:
        return f"ends {hour_iso(max(hours))}, less than 72 h ahead (an old run?)"
    return ""


def parse_om_marine(js, points, now, lead_h, min_cell_lon):
    """Open-Meteo marine (one object per point, in the order asked) -> [series or why] per point. The Marconi
    point (first) is refused when its grid cell lies west of min_cell_lon (Cape Cod Bay / behind the shore)."""
    objs = om_objects(js)
    out = []
    for i, (lat, lon) in enumerate(points):
        if i >= len(objs):
            out.append("no object for this point")
            continue
        o = objs[i] or {}
        try:
            clat, clon = float(o["latitude"]), float(o["longitude"])
        except (KeyError, TypeError, ValueError):
            out.append("no grid cell")
            continue
        if abs(clat - lat) > 0.3 or abs(clon - lon) > 0.3:
            out.append(f"grid cell {clat:.2f},{clon:.2f} is not near {lat:.2f},{lon:.2f}")
            continue
        if i == 0 and clon < min_cell_lon:
            out.append(f"grid cell {clat:.2f}N {-clon:.2f}W lies west of {-min_cell_lon:.2f}W (Cape Cod Bay "
                       "or behind the shore)")
            continue
        units = o.get("hourly_units") or {}
        if units.get("wave_height") != "m":
            out.append(f"wave_height unit {units.get('wave_height')!r}, not m")
            continue
        hr = o.get("hourly") or {}
        try:
            hours = om_times(hr)
        except (KeyError, ValueError, TypeError):
            out.append("unreadable times")
            continue
        why = om_fresh(hours, now)
        if why:
            out.append(why)
            continue
        vals = {}
        for col, key in (("hs", "wave_height"), ("tp", "wave_period"), ("dir", "wave_direction")):
            arr = hr.get(key) or []
            vals[col] = {h: (np.nan if v is None else float(v)) for h, v in zip(hours, arr)}
        s = make_series(vals, ["hs", "tp", "dir"], hold=2)
        why = check_series(s, ["hs", "tp", "dir"], now, lead_h)
        if why:
            out.append(why)
            continue
        s.update(tstat="Tm", cell=f"{clat:.2f}N {-clon:.2f}W")
        out.append(s)
    return out


def parse_om_wind(js, now, lead_h):
    o = om_objects(js)[0] or {}
    units = o.get("hourly_units") or {}
    if units.get("wind_speed_10m") != "m/s":
        raise ValueError(f"wind unit {units.get('wind_speed_10m')!r}, not m/s")
    hr = o.get("hourly") or {}
    hours = om_times(hr)
    why = om_fresh(hours, now)
    if why:
        raise ValueError(why)
    vals = {}
    for col, key in (("ws", "wind_speed_10m"), ("gust", "wind_gusts_10m"), ("wdir", "wind_direction_10m")):
        vals[col] = {h: (np.nan if v is None else float(v)) for h, v in zip(hours, hr.get(key) or [])}
    s = make_series(vals, ["ws", "gust", "wdir"], hold=2)
    why = check_series(s, ["ws", "gust", "wdir"], now, lead_h)
    if why:
        raise ValueError(why)
    return s


BULL_CYCLE = re.compile(r"Cycle\s*:\s*(\d{8})\s+(\d{1,2})\s*UTC")
BULL_ROW = re.compile(r"^\|\s*(\d{1,2})\s+(\d{1,2})\s*\|\s*([\d.]+)")


def parse_gfs_bull(text, now, lead_h):
    """NOAA GFS-Wave point bulletin (gfswave.<station>.bull) -> (series, cycle datetime). Rows '| dd hh | Hst n x |
    [*] Hs Tp dir | ...': Hst the total height; Tp and direction from the largest partition. Month and year from
    the Cycle line; a day smaller than the row before rolls to the next month. Stale over 36 h (ValueError)."""
    m = BULL_CYCLE.search(text or "")
    if not m:
        raise ValueError("no Cycle line")
    cyc = datetime.strptime(m.group(1) + m.group(2).zfill(2), "%Y%m%d%H").replace(tzinfo=timezone.utc)
    if (utc(now) - cyc).total_seconds() > 36 * 3600:
        raise ValueError(f"stale: cycle {cyc:%Y-%m-%d %H}Z")
    year, month, prev_day = cyc.year, cyc.month, None
    vals = {"hs": {}, "tp": {}, "dir": {}}
    for line in text.splitlines():
        r = BULL_ROW.match(line.strip())
        if not r:
            continue
        day, hr, hst = int(r.group(1)), int(r.group(2)), float(r.group(3))
        if prev_day is not None and day < prev_day:
            month += 1
            if month > 12:
                month, year = 1, year + 1
        prev_day = day
        try:
            t = datetime(year, month, day, hr, tzinfo=timezone.utc)
        except ValueError:
            continue
        best = None
        for part in line.strip().split("|")[3:]:
            nums = part.replace("*", " ").split()
            if len(nums) == 3:
                try:
                    hs, tp, d = (float(x) for x in nums)
                except ValueError:
                    continue
                if best is None or hs > best[0]:
                    best = (hs, tp, d)
        h = hour_of(t)
        vals["hs"][h] = hst
        vals["tp"][h] = best[1] if best else np.nan
        vals["dir"][h] = best[2] if best else np.nan
    s = make_series(vals, ["hs", "tp", "dir"], hold=2)
    if s is None:
        raise ValueError("no rows")
    why = check_series(s, ["hs", "tp", "dir"], now, lead_h)
    if why:
        raise ValueError(why)
    s.update(tstat="Tp", cell="44008 model point")
    return s, cyc

# --------------------------------------------------------------------------------------------
# fetching the forecast (one call per source; every failure is "unavailable this run")
# --------------------------------------------------------------------------------------------

SOURCE_FILES = ("nws_points.json", "nws_grid.json", "nws_forecast.json", "nws_alerts.json", "om_marine.json",
                "om_wind.json", "gfswave_44008.bull")
MARINE_PREFIX = ("AN", "AM", "GM", "PZ", "PK", "PH", "PS", "LE", "LH", "LM", "LO", "LS", "SL", "LC")
GFS_URL = "https://polar.ncep.noaa.gov/waves/WEB/gfswave.latest_run/plots/gfswave.44008.bull"


def coord(x):
    return ("%.4f" % x).rstrip("0").rstrip(".")


def om_marine_url(cfg, days=10, cell_sea=True):
    (a, b), (c, d) = cfg["STORM_MARINE_POINT"], cfg["STORM_BUOY_POINT"]
    return ("https://marine-api.open-meteo.com/v1/marine?latitude=%s,%s&longitude=%s,%s"
            "&hourly=wave_height,wave_direction,wave_period,swell_wave_height,swell_wave_direction,"
            "swell_wave_period,wind_wave_height&timezone=GMT&forecast_days=%d%s"
            % (coord(a), coord(c), coord(b), coord(d), days, "&cell_selection=sea" if cell_sea else ""))


def om_wind_url(cfg):
    lat, lon = cfg["STORM_POINT"]
    return ("https://api.open-meteo.com/v1/forecast?latitude=%s&longitude=%s&hourly=wind_speed_10m,"
            "wind_gusts_10m,wind_direction_10m,pressure_msl&wind_speed_unit=ms&timezone=GMT&forecast_days=10"
            % (coord(lat), coord(lon)))


def zone_url(z):
    return f"https://api.weather.gov/zones/{'coastal' if z[:2] in MARINE_PREFIX else 'forecast'}/{z}"


class Fetcher:
    """Fetches (or reads from --forecast-dir) each source; archives every 200 body gzipped (live runs only)."""

    def __init__(self, ctx):
        self.ctx = ctx
        self.raw = {}
        self.calls = []
        self.deadline = time.time() + ctx.cfg["STORM_FETCH_DEADLINE_S"]
        self.stamp = utc(ctx.now).strftime("%H%MZ")

    def get(self, name, url):
        ctx = self.ctx
        if ctx.forecast_dir:
            body = read_saved(ctx.forecast_dir, name)
            res = (200, body, "") if body is not None else (0, b"", f"no {name} in {ctx.forecast_dir}")
        else:
            res = http_get(url, headers_for(url, ctx.cfg), ctx.cfg["STORM_HTTP_TIMEOUT_S"], opener=ctx.opener,
                           sleep=ctx.sleep, deadline=self.deadline)
            self.calls.append((name, url, res[0]))
        self.raw[name] = res
        if res[0] == 200:
            self.archive(name, res[1])
        return res

    def archive(self, name, body):
        ctx = self.ctx
        if ctx.save_raw:
            try:
                Path(ctx.save_raw).mkdir(parents=True, exist_ok=True)
                (Path(ctx.save_raw) / name).write_bytes(body)
            except OSError as exc:
                ctx.note(f"--save-raw {name}: {exc}")
        if ctx.mode != "live" or ctx.forecast_dir:
            return
        d = ctx.forecasts_dir / utc(ctx.now).strftime("%Y-%m-%d")
        try:
            d.mkdir(parents=True, exist_ok=True)
            with gzip.open(str(d / f"{self.stamp}_{name}.gz"), "wb") as f:
                f.write(body)
        except OSError as exc:
            ctx.note(f"forecast archive {name}: {exc}")


def read_saved(folder, name):
    """<name>, <name>.gz or the latest <HHMMZ>_<name>[.gz] in a folder -> bytes or None."""
    d = Path(folder)
    for cand in (d / name, d / (name + ".gz")):
        if cand.exists():
            b = cand.read_bytes()
            return gzip.decompress(b) if cand.suffix == ".gz" else b
    found = sorted(p for p in d.glob(f"*_{name}*") if p.name.endswith((name, name + ".gz")))
    if found:
        b = found[-1].read_bytes()
        return gzip.decompress(b) if found[-1].suffix == ".gz" else b
    return None


def fetch_forecast(ctx, state):
    """Fetches every source in STORM_SOURCES. -> Fetcher (raw bodies by name). Refreshes the cached NWS
    /points answer weekly or after a 404 from the grid; never hardcodes the grid numbers."""
    cfg, f = ctx.cfg, Fetcher(ctx)
    srcs = cfg["STORM_SOURCES"]
    if "nws" in srcs:
        pts = state.get("nws_points") or {}
        fresh = pts.get("utc") and (utc(ctx.now) - parse_iso(pts["utc"])).total_seconds() < 7 * 86400
        refreshed = False

        def refresh():
            lat, lon = cfg["STORM_POINT"]
            st, body, why = f.get("nws_points.json", f"https://api.weather.gov/points/{coord(lat)},{coord(lon)}")
            if st == 200:
                try:
                    p = parse_nws_points(json.loads(body.decode("utf-8", "replace")))
                    p["utc"] = iso_z(ctx.now)
                    if not ctx.forecast_dir:
                        state["nws_points"] = p
                    return p
                except ValueError as exc:
                    ctx.note(f"nws /points unreadable ({exc})")
            return None
        if ctx.forecast_dir and read_saved(ctx.forecast_dir, "nws_points.json") is not None:
            pts, refreshed = refresh() or pts, True
        elif not fresh:
            pts, refreshed = refresh() or {}, True
        if pts.get("forecastGridData"):
            st, _, _ = f.get("nws_grid.json", pts["forecastGridData"])
            if st == 404 and not refreshed and not ctx.forecast_dir:
                state.pop("nws_points", None)
                pts, refreshed = refresh() or {}, True
                if pts.get("forecastGridData"):
                    f.get("nws_grid.json", pts["forecastGridData"])
            if st == 404:
                state.pop("nws_points", None)
            if pts.get("forecast"):
                f.get("nws_forecast.json", pts["forecast"])
        else:
            f.raw["nws_grid.json"] = (0, b"", "no /points answer (grid URL unknown)")
        zones = cfg["STORM_ALERT_ZONES"]
        f.get("nws_alerts.json", "https://api.weather.gov/alerts/active?zone=" + ",".join(zones))
        if refreshed and not ctx.forecast_dir:
            fz = str(pts.get("forecastZone") or "")
            if fz and not any(fz.endswith("/" + z) for z in zones):
                ctx.note(f"NWS /points names forecast zone {fz.rsplit('/', 1)[-1]}, not in STORM_ALERT_ZONES")
            for z in zones:
                st, _, why = f.get(f"nws_zone_{z}.json", zone_url(z))
                if st == 404:
                    ctx.note(f"NWS zone {z} not found (404): renamed or retired? check STORM_ALERT_ZONES")
    if "om_marine" in srcs:
        st, body, why = f.get("om_marine.json", om_marine_url(cfg))
        if st == 400 and not ctx.forecast_dir:
            st, body, why = f.get("om_marine.json", om_marine_url(cfg, 10, False))
            if st == 400:
                f.get("om_marine.json", om_marine_url(cfg, 7, False))
    if "om_wind" in srcs:
        f.get("om_wind.json", om_wind_url(cfg))
    if "gfswave" in srcs:
        f.get("gfswave_44008.bull", GFS_URL)
    return f


def loads(res):
    st, body, why = res
    if st != 200:
        raise ValueError(why or f"HTTP {st}")
    return json.loads(body.decode("utf-8", "replace"))


def build_forecast(raw, now, cfg):
    """Raw bodies by name -> the forecast: the wave series the rules read (the first fresh group of: both
    Open-Meteo points, the GFS-Wave 44008 bulletin, the NWS grid waves), the others (a NOTE only), the wind
    (NWS grid, else Open-Meteo), the NWS alerts, and each source's status."""
    lead = cfg["STORM_LEAD_HOURS"]
    fc = {"now": iso_z(now), "status": {}, "groups": {}, "wind_all": [], "alerts": None, "text": [], "notes": [],
          "om_used": False}
    st = fc["status"]

    def unavailable(name):
        r = raw.get(name)
        return "not fetched" if r is None else (r[2] or f"HTTP {r[0]}")
    if "nws_grid.json" in raw:
        try:
            g = parse_nws_grid(loads(raw["nws_grid.json"]), now, lead)
            if g["stale"]:
                st["nws_grid"] = f"stale (updated {g['updated']})"
            else:
                st["nws_grid"] = "ok" + ("" if g["waves"] is not None else " (wind only)")
                if g["wind"] is not None:
                    g["wind"].update(name="nws_grid wind")
                    fc["wind_all"].append(g["wind"])
                else:
                    st["nws_grid"] = f"no wind ({g['why_wind']})"
                if g["waves"] is not None:
                    g["waves"].update(name="nws_grid waves", tstat="Tp", cell="NWS grid cell")
                    fc["groups"]["nws_grid"] = [g["waves"]]
            fc["notes"] += [f"nws_grid {n}" for n in g["notes"]]
        except (ValueError, KeyError, TypeError) as exc:
            st["nws_grid"] = f"unavailable ({str(exc)[:80]})"
    if "nws_alerts.json" in raw:
        try:
            fc["alerts"] = parse_alerts(loads(raw["nws_alerts.json"]), cfg["STORM_ALERT_ZONES"])
            codes = sorted({a["code"] or a["event"] for a in fc["alerts"]})
            st["nws_alerts"] = "ok" + (f" ({', '.join(codes)})" if codes else " (none)")
        except (ValueError, KeyError, TypeError) as exc:
            st["nws_alerts"] = f"unavailable ({str(exc)[:80]})"
    if "nws_forecast.json" in raw:
        try:
            per = (loads(raw["nws_forecast.json"]).get("properties") or {}).get("periods") or []
            fc["text"] = [f"{p.get('name', '')}: {p.get('detailedForecast', '')}" for p in per[:2]]
        except (ValueError, AttributeError, TypeError):
            pass
    if "om_marine.json" in raw:
        try:
            pts = [cfg["STORM_MARINE_POINT"], cfg["STORM_BUOY_POINT"]]
            res = parse_om_marine(loads(raw["om_marine.json"]), pts, now, lead, cfg["STORM_MIN_CELL_LON"])
            ok = []
            for (lat, lon), r, label in zip(pts, res, ("offshore Marconi", "44008 point")):
                if isinstance(r, dict):
                    r.update(name=f"om_marine {label} ({r['cell']})")
                    ok.append(r)
                else:
                    fc["notes"].append(f"om_marine {label}: {r}")
            if ok:
                fc["groups"]["om_marine"] = ok
                st["om_marine"] = "ok" + ("" if len(ok) == 2 else " (one point)")
            else:
                st["om_marine"] = "unavailable (" + "; ".join(r for r in res if isinstance(r, str))[:120] + ")"
        except (ValueError, KeyError, TypeError) as exc:
            st["om_marine"] = f"unavailable ({str(exc)[:80]})"
    if "om_wind.json" in raw:
        try:
            w = parse_om_wind(loads(raw["om_wind.json"]), now, lead)
            w.update(name="om_wind")
            fc["wind_all"].append(w)
            st["om_wind"] = "ok"
        except (ValueError, KeyError, TypeError, IndexError) as exc:
            st["om_wind"] = f"unavailable ({str(exc)[:80]})"
    if "gfswave_44008.bull" in raw:
        try:
            r = raw["gfswave_44008.bull"]
            if r[0] != 200:
                raise ValueError(r[2] or f"HTTP {r[0]}")
            s, cyc = parse_gfs_bull(r[1].decode("utf-8", "replace"), now, lead)
            s.update(name=f"gfswave 44008 (cycle {cyc:%d %b %H}Z)")
            fc["groups"]["gfswave_44008"] = [s]
            st["gfswave_44008"] = "ok"
        except (ValueError, TypeError) as exc:
            st["gfswave_44008"] = f"unavailable ({str(exc)[:80]})"
    for name in raw:
        key = name.split(".")[0]
        if key.startswith(("nws_points", "nws_zone", "nws_forecast")) or key in st:
            continue
        st[key] = f"unavailable ({unavailable(name)[:80]})"
    order = [g for g in ("om_marine", "gfswave_44008", "nws_grid") if g in fc["groups"]]
    fc["waves_group"] = order[0] if order else None
    fc["waves"] = fc["groups"][order[0]] if order else []
    fc["waves_other"] = [s for g in order[1:] for s in fc["groups"][g]]
    fc["wind"] = fc["wind_all"][0] if fc["wind_all"] else None
    fc["om_used"] = fc["waves_group"] == "om_marine" or (fc["wind"] is not None and fc["wind"]["name"] == "om_wind")
    fc["blind"] = not fc["waves"] and fc["wind"] is None and fc["alerts"] is None
    return fc


def perfect_forecast(obs, now, cfg):
    """Replay: the observed offshore waves over (now, now + lead] as the forecast (raw 44008, else the
    Marconi-equivalent / 0.83), labelled 'hindcast: observed waves'."""
    lead = cfg["STORM_LEAD_HOURS"]
    h0 = hour_of(now)
    fc = {"now": iso_z(now), "status": {}, "groups": {}, "wind_all": [], "alerts": None, "text": [], "notes": [],
          "om_used": False, "wind": None, "waves_other": [], "hindcast": True}
    if obs is not None and obs.n:
        i = np.arange(obs.n)
        hrs = obs.h0 + i
        sel = (hrs > h0) & (hrs <= h0 + lead)
        hs = np.where(np.isfinite(obs.raw_hs), obs.raw_hs, obs.prim / cfg["STORM_OFFSHORE_TO_MARCONI"])[sel]
        if np.isfinite(hs).sum() >= MIN_SOURCE_HOURS:
            s = {"h": hrs[sel], "hs": hs, "tp": obs.tp[sel], "dir": obs.dir[sel], "tstat": "Tp",
                 "cell": "observed", "name": "hindcast: observed waves (44008, else Marconi / 0.83)"}
            fc["groups"]["observed"] = [s]
    fc["waves_group"] = "observed" if fc["groups"] else None
    fc["waves"] = fc["groups"].get("observed", [])
    fc["status"]["hindcast"] = "observed waves as a perfect forecast" if fc["waves"] else "no observed waves ahead"
    fc["blind"] = not fc["waves"]
    return fc


# --------------------------------------------------------------------------------------------
# the trigger: a pure function of the forecast, the time and the settings
# --------------------------------------------------------------------------------------------


def true_runs(mask):
    """[(i0, i1)] inclusive index runs of True."""
    out, i, n = [], 0, len(mask)
    while i < n:
        if mask[i]:
            j = i
            while j + 1 < n and mask[j + 1]:
                j += 1
            out.append((i, j))
            i = j + 1
        else:
            i += 1
    return out


def first_run(mask, n):
    for a, b in true_runs(mask):
        if b - a + 1 >= n:
            return a, b
    return None


def series_at(s, hours, col):
    """Values of series column at the given hours (NaN outside)."""
    out = np.full(len(hours), np.nan)
    if s is None or col not in s:
        return out
    k = hours - s["h"][0]
    ok = (k >= 0) & (k < len(s["h"]))
    out[ok] = s[col][k[ok]]
    return out


def alert_hours(a):
    return iso_hour(a["start"]), iso_hour(a["end"])


RULE_WORDS = {"a": "waves at or above the storm height", "b": "a strong NWS alert",
              "c": "waves at or above the support height with a supporting NWS alert",
              "d": "long-period swell", "e": "gale-force onshore wind (no wave forecast)"}


def evaluate_trigger(fc, now, cfg):
    """-> {triggered, rules, onset_fc, end_fc, peak, first_trigger_hour, hours_over, outlook, notes, by}."""
    lead = cfg["STORM_LEAD_HOURS"]
    h0 = hour_of(now)
    h1 = h0 + lead
    sector, onset_hs = cfg["STORM_DIR_SECTOR"], cfg["STORM_ONSET_HS_M"]
    T = {"triggered": False, "rules": [], "by": [], "onset_fc": None, "end_fc": None, "peak": None,
         "first_trigger_hour": None, "hours_over": 0, "outlook": [], "notes": [], "alerts": []}
    alerts = [a for a in (fc.get("alerts") or []) if alert_hours(a)[0] <= h1 and alert_hours(a)[1] >= h0]
    have_waves = bool(fc.get("waves"))
    gated = set(cfg["STORM_ALERT_WAVE_GATED"])
    strong_codes = set(cfg["STORM_ALERT_STRONG"]) - (gated if have_waves else set())
    support_codes = set(cfg["STORM_ALERT_SUPPORT"]) | (gated if have_waves else set())
    strong = [a for a in alerts if a["code"] in strong_codes]
    support = [a for a in alerts if a["code"] in support_codes]
    T["alerts"] = [dict(a, role="strong" if a in strong else "support" if a in support else "other") for a in alerts]
    onsets, ends, firsts = [], [], []

    def masks(s):
        w = (s["h"] >= h0) & (s["h"] <= h1)
        d = s["dir"].copy()
        nodir = ~np.isfinite(d)
        if nodir.any() and fc.get("wind") is not None:
            d[nodir] = series_at(fc["wind"], s["h"][nodir], "wdir")
        still = ~np.isfinite(d)
        sec = in_sector(d, sector) | still
        if (still & w & (s["hs"] >= onset_hs)).any():
            T["notes"].append(f"{s['name']}: no wave or wind direction at some rough hours: taken as in the sector")
        return w, sec

    def wave_rules(s, record=True):
        w, sec = masks(s)
        hs = np.nan_to_num(s["hs"], nan=-1.0)
        tp = np.nan_to_num(s.get("tp", np.full(len(hs), np.nan)), nan=-1.0)
        thr_tp = cfg["STORM_SWELL_TP_S"] if s.get("tstat", "Tp") == "Tp" else cfg["STORM_SWELL_TM_S"]
        n = cfg["STORM_MIN_HOURS"]
        found = {}
        ra = first_run(w & sec & (hs >= cfg["STORM_HS_M"]), n)
        if ra:
            found["a"] = ra
        if support:
            rc = first_run(w & sec & (hs >= cfg["STORM_HS_SUPPORT_M"]), n)
            if rc:
                found["c"] = rc
        rd = first_run(w & sec & (hs >= cfg["STORM_SWELL_HS_M"]) & (tp >= thr_tp), n)
        if rd:
            found["d"] = rd
        if record:
            T["hours_over"] = max(T["hours_over"], int((w & sec & (hs >= cfg["STORM_HS_M"])).sum()))
            for rule, (a, b) in found.items():
                if rule not in T["rules"]:
                    T["rules"].append(rule)
                T["by"].append(f"{rule}: {s['name']}")
                i = a
                j0 = i
                while j0 - 1 >= 0 and hs[j0 - 1] > onset_hs and sec[j0 - 1]:
                    j0 -= 1
                j1 = i
                while j1 + 1 < len(hs) and hs[j1 + 1] > onset_hs and sec[j1 + 1]:
                    j1 += 1
                onsets.append(int(s["h"][j0]) - 12)
                ends.append(int(s["h"][j1]))
                firsts.append(int(s["h"][a]))
        return found

    for s in fc.get("waves") or []:
        wave_rules(s)
    for s in fc.get("waves") or []:
        w, _ = masks(s)
        if np.isfinite(s["hs"][w]).any():
            k = int(np.nanargmax(np.where(w, s["hs"], np.nan)))
            if T["peak"] is None or s["hs"][k] > T["peak"]["hs"]:
                T["peak"] = {"hs": round(float(s["hs"][k]), 2), "time": hour_iso(int(s["h"][k])),
                             "dir": None if not np.isfinite(s["dir"][k]) else round(float(s["dir"][k])),
                             "tp": None if not np.isfinite(s["tp"][k]) else round(float(s["tp"][k]), 1),
                             "tstat": s.get("tstat", "Tp"), "source": s["name"]}
    if strong:
        T["rules"].append("b")
        for a in strong:
            T["by"].append(f"b: {a['code']} {a['event']} ({', '.join(a['zones'])})")
            ah0, ah1 = alert_hours(a)
            run = None
            for s in fc.get("waves") or []:
                w, sec = masks(s)
                rough = (np.nan_to_num(s["hs"], nan=-1.0) > onset_hs) & sec
                for i0, i1 in true_runs(rough):
                    if s["h"][i0] <= ah1 and s["h"][i1] >= ah0 - 24:
                        cand = (int(s["h"][i0]), int(s["h"][i1]))
                        run = cand if run is None or cand[0] < run[0] else run
                        break
            onsets.append((run[0] if run else ah0) - 12)
            ends.append(max(ah1, run[1]) if run else ah1)
            firsts.append(max(ah0, h0))
    if not have_waves and fc.get("wind") is not None:
        s = fc["wind"]
        w = (s["h"] >= h0) & (s["h"] <= h1)
        ws, gust = (np.nan_to_num(s[c], nan=-1.0) for c in ("ws", "gust"))
        cond = w & in_sector(s["wdir"], cfg["STORM_WIND_DIR_SECTOR"]) & (
            (ws >= cfg["STORM_WIND_MS"]) | ((ws >= cfg["STORM_WIND_MIN_MS"]) & (gust >= cfg["STORM_WIND_GUST_MS"])))
        r = first_run(cond, cfg["STORM_WIND_MIN_HOURS"])
        if r:
            T["rules"].append("e")
            T["by"].append(f"e: {s['name']}")
            j1 = r[1]
            while j1 + 1 < len(cond) and cond[j1 + 1]:
                j1 += 1
            onsets.append(int(s["h"][r[0]]) - 12)
            ends.append(int(s["h"][j1]))
            firsts.append(int(s["h"][r[0]]))
    if T["rules"]:
        T["triggered"] = True
        T["rules"] = sorted(set(T["rules"]))
        on = min(onsets)
        T["onset_fc"] = hour_iso(on)
        T["end_fc"] = hour_iso(max(max(ends), on + 12) if ends else on + 48)
        T["first_trigger_hour"] = hour_iso(min(firsts))
    elif fc.get("waves_other"):
        for s in fc["waves_other"]:
            if wave_rules(s, record=False):
                T["notes"].append(f"NOTE: {s['name']} would trigger, but it is not the primary wave source "
                                  f"({fc.get('waves_group')}); no trigger")
                break
    for s in fc.get("groups", {}).get("om_marine", []):
        out = (s["h"] > h1) & (s["h"] <= h0 + 240)
        for d in sorted({hour_date(int(x)) for x in s["h"][out]}):
            k = out & (np.array([hour_date(int(x)) == d for x in s["h"]]))
            if np.isfinite(s["hs"][k]).any():
                T["outlook"].append((d, round(float(np.nanmax(s["hs"][k])), 1), s["name"]))
        break
    return T

# --------------------------------------------------------------------------------------------
# observed waves, in the Marconi currency (ADCP-equivalent at 21 m: hs_best, DEM_MAX_HS)
# --------------------------------------------------------------------------------------------

OBS_DAYS = 60
SRC_LABELS = ["", "waves_marconi", "44013 converted", "44008 x 0.83", "camera only"]


def read_rows(path, cols, h_min, h_max):
    """Rows of a wave CSV with epoch in [h_min, h_max] hours, as {col: text}; a truncated last line is skipped."""
    out = []
    p = Path(path)
    if not p.exists():
        return out
    try:
        with open(p, newline="", errors="replace") as f:
            r = csv.reader(f)
            head = next(r, [])
            if "epoch" not in head:
                return out
            idx = {c: head.index(c) for c in cols if c in head}
            ke = head.index("epoch")
            for row in r:
                if len(row) != len(head):
                    continue
                try:
                    h = int(float(row[ke]) // 3600)
                except ValueError:
                    continue
                if h_min <= h <= h_max:
                    d = {c: row[i] for c, i in idx.items()}
                    d["h"] = h
                    out.append(d)
    except OSError:
        pass
    return out


def fnum_or_nan(text):
    try:
        v = float(text)
        return v if np.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan


def hourly_buoy(rows):
    """{hour: (hs, tp, dir)} from fetch_buoy_waves.py rows (mean in the hour; direction a circular mean)."""
    acc = {}
    for r in rows:
        hs = fnum_or_nan(r.get("wvht_m"))
        if not np.isfinite(hs):
            continue
        a = acc.setdefault(r["h"], [[], [], []])
        a[0].append(hs)
        a[1].append(fnum_or_nan(r.get("dpd_s")))
        a[2].append(fnum_or_nan(r.get("mwd_deg")))
    out = {}
    for h, (hs, tp, d) in acc.items():
        d = np.array(d, float)
        d = d[np.isfinite(d)]
        mean_d = float(np.degrees(np.arctan2(np.sin(np.radians(d)).mean(), np.cos(np.radians(d)).mean())) % 360) \
            if len(d) else np.nan
        tp = np.array(tp, float)
        out[h] = (float(np.mean(hs)), float(np.nanmean(tp)) if np.isfinite(tp).any() else np.nan, mean_d)
    return out


class Obs:
    """Hourly observed waves on [h0, h0 + n): prim (Marconi-equivalent from the buoy layers: waves_marconi
    blend/buoy rows, else 44013 converted, else 44008 x 0.83), cam (waves_marconi camera-only rows), val (the
    calm test: the larger of prim and cam), raw 44008 height and direction, period and direction."""

    def __init__(self, h0, n):
        self.h0, self.n = h0, n
        self.prim = np.full(n, np.nan)
        self.psrc = np.zeros(n, int)
        self.cam = np.full(n, np.nan)
        self.raw_hs = np.full(n, np.nan)
        self.raw_dir = np.full(n, np.nan)
        self.tp = np.full(n, np.nan)
        self.dir = np.full(n, np.nan)
        self.notes = []
        self.newest = None

    @property
    def val(self):
        with np.errstate(invalid="ignore"):
            return np.where(np.isfinite(self.prim) & np.isfinite(self.cam), np.maximum(self.prim, self.cam),
                            np.where(np.isfinite(self.prim), self.prim, self.cam))


def offshore_to_marconi(hs, d, cfg):
    sec = in_sector(d, cfg["STORM_DIR_SECTOR"]) | ~np.isfinite(d)
    return np.where(sec, cfg["STORM_OFFSHORE_TO_MARCONI"] * hs, cfg["STORM_SHELTERED_FACTOR"] * hs)


def read_observed(base, now, cfg, h_end=None):
    """Reads archive/waves_marconi.csv, waves_44013.csv and waves_44008.csv over the last 60 days up to `now`
    -> Obs. Reads only; the caller holds owg.sh's lock."""
    base = Path(base)
    h1 = hour_of(now) if h_end is None else h_end
    h0 = h1 - OBS_DAYS * 24
    o = Obs(h0, h1 - h0 + 1)
    m = read_rows(base / "archive/waves_marconi.csv", ["hs_best", "best_from", "qc_flag", "dir_deg", "tp_s"], h0, h1)
    b13 = hourly_buoy(read_rows(base / "archive/waves_44013.csv", ["wvht_m", "dpd_s", "mwd_deg"], h0, h1))
    b08 = hourly_buoy(read_rows(base / "archive/waves_44008.csv", ["wvht_m", "dpd_s", "mwd_deg"], h0, h1))
    L1, cam, raw_m = {}, {}, {}
    for r in m:
        q = fnum_or_nan(r.get("qc_flag"))
        hs = fnum_or_nan(r.get("hs_best"))
        if not np.isfinite(hs) or (np.isfinite(q) and q >= 4):
            continue
        bf = (r.get("best_from") or "").strip()
        if bf in ("blend", "buoy"):
            L1[r["h"]] = (hs, fnum_or_nan(r.get("tp_s")), fnum_or_nan(r.get("dir_deg")))
        elif bf == "camera":
            cam[r["h"]] = hs
        elif bf == "buoy (raw)":
            raw_m[r["h"]] = (hs, fnum_or_nan(r.get("tp_s")), fnum_or_nan(r.get("dir_deg")))
    for h, v in raw_m.items():
        b08.setdefault(h, v)
    L2 = {}
    if b13:
        hs13 = np.array([v[0] for v in b13.values()])
        tp13 = np.array([v[1] for v in b13.values()])
        d13 = np.array([v[2] for v in b13.values()])
        conv = np.full(len(hs13), np.nan)
        trp = base / "calibration/buoy_transfer_44013.json"
        if trp.exists():
            try:
                sys.path.insert(0, str(HERE))
                import buoy_transfer
                conv = np.asarray(buoy_transfer.apply(buoy_transfer.load(trp), hs13, tp13, d13), float)
            except Exception as exc:
                o.notes.append(f"buoy_transfer_44013.json not applied ({exc})")
        sec = in_sector(d13, cfg["STORM_DIR_SECTOR"]) | ~np.isfinite(d13)
        fallback = np.where(sec, 0.827 * hs13 + 0.151, cfg["STORM_SHELTERED_FACTOR"] * hs13)
        conv = np.where(np.isfinite(conv), conv, fallback)
        for (h, v), c in zip(b13.items(), conv):
            L2[h] = (float(c), v[1], v[2])
    for h, (hs, tp, d) in b08.items():
        if h0 <= h <= h1:
            k = h - h0
            o.raw_hs[k], o.raw_dir[k] = hs, d
    for k in range(o.n):
        h = h0 + k
        if h in L1:
            v, src = L1[h], 1
        elif h in L2:
            v, src = L2[h], 2
        elif h in b08:
            hs, tp, d = b08[h]
            v, src = (float(offshore_to_marconi(np.array([hs]), np.array([d]), cfg)[0]), tp, d), 3
        else:
            v, src = None, 0
        if v is not None:
            o.prim[k], o.psrc[k] = v[0], src
            o.tp[k], o.dir[k] = v[1], v[2]
        if h in cam:
            o.cam[k] = cam[h]
            if src == 0:
                o.psrc[k] = 4
        if h in b08:
            if not np.isfinite(o.tp[k]):
                o.tp[k] = b08[h][1]
            if not np.isfinite(o.dir[k]):
                o.dir[k] = b08[h][2]
    have = np.flatnonzero(np.isfinite(o.val))
    o.newest = int(h0 + have[-1]) if len(have) else None
    return o


def obs_hour_value(o, h):
    k = h - o.h0
    return o.val[k] if 0 <= k < o.n else np.nan


def segments(o, cfg):
    """Storm segments of the observed series: a run of rough hours (> calm) ending at the first present calm
    hour t whose next 24 h hold no rough hour and >= 20 values (a missing hour never counts as calm). The peak,
    the confirm and open tests use the buoy-based layers only (never camera-only hours)."""
    calm, n_end, n_min = cfg["STORM_CALM_HS"], cfg["STORM_END_CALM_HOURS"], cfg["STORM_END_MIN_PRESENT"]
    v = o.val
    with np.errstate(invalid="ignore"):
        rough = np.isfinite(v) & (v > calm)
        calm_h = np.isfinite(v) & (v <= calm)
    segs, i, n = [], 0, o.n
    sec_raw = in_sector(o.raw_dir, cfg["STORM_DIR_SECTOR"])
    while i < n:
        if not rough[i]:
            i += 1
            continue
        start, last, end, j = i, i, None, i + 1
        while j < n:
            if rough[j]:
                last = j
            elif calm_h[j] and j + n_end <= n:
                if not rough[j:j + n_end].any() and int(np.isfinite(v[j:j + n_end]).sum()) >= n_min:
                    end = j
                    break
            j += 1
        span = slice(start, last + 1)
        prim = o.prim[span]
        has = np.isfinite(prim)
        seg = {"start": o.h0 + start, "last": o.h0 + last, "end": None if end is None else o.h0 + end,
               "rough_hours": int(rough[span].sum()), "camera_only": not has.any(),
               "peak": None, "peak_h": None, "n_confirm": 0,
               "srcs": sorted({SRC_LABELS[int(x)] for x in o.psrc[span] if x})}
        if has.any():
            k = int(np.nanargmax(np.where(has, prim, np.nan)))
            seg.update(peak=round(float(prim[k]), 2), peak_h=o.h0 + start + k,
                       n_confirm=int((np.nan_to_num(prim, nan=-1) >= cfg["STORM_CONFIRM_HS"]).sum()))
        raw = np.where(sec_raw[span], np.nan_to_num(o.raw_hs[span], nan=-1), -1)
        seg["raw_max"] = round(float(raw.max()), 2) if (raw > 0).any() else None
        seg["raw_open"] = first_run(raw >= cfg["STORM_OBS_OPEN_BUOY_HS"], 2) is not None
        seg["raw_confirm"] = first_run(raw >= cfg["STORM_CONFIRM_BUOY_HS"], 2) is not None
        segs.append(seg)
        i = end if end is not None else n
    return segs


def confirms(s, cfg):
    return (s["peak"] is not None and s["peak"] >= cfg["STORM_CONFIRM_HS"]) or s["raw_confirm"]


def opens(s, cfg):
    return ((s["peak"] is not None and s["peak"] >= cfg["STORM_OBS_OPEN_HS"] and s["n_confirm"] >= 2)
            or s["raw_open"])


def breaks(s, cfg):
    """A wave event that cuts a post window: it confirms, or peaks >= STORM_BREAK_HS, or is rough >= STORM_BREAK_HOURS."""
    if confirms(s, cfg):
        return True
    bh, bn = cfg["STORM_BREAK_HS"], cfg["STORM_BREAK_HOURS"]
    return ((bh is not None and s["peak"] is not None and s["peak"] >= bh)
            or (bn is not None and s["rough_hours"] >= bn))


def seg_span(s, now_h):
    return s["start"], (s["end"] if s["end"] is not None else now_h)


def seg_text(s):
    return (f"{hour_iso(s['start'])[:13]}Z-" + (f"{hour_iso(s['end'])[:13]}Z" if s["end"] is not None else "open")
            + (f", peak {s['peak']:.2f} m" if s["peak"] is not None else ", camera-only")
            + (f", 44008 {s['raw_max']:.1f} m" if s["raw_max"] else ""))


# --------------------------------------------------------------------------------------------
# locks (flock(2), as bash's flock): the cron's, owg.sh's, our own
# --------------------------------------------------------------------------------------------


def take_lock(path, wait_s=0.0, sleep=time.sleep, poll=10.0):
    """flock(LOCK_EX|LOCK_NB) on `path`, polled every `poll` s up to wait_s. -> file object or None."""
    import fcntl
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    f = open(path, "a")
    t_end = time.time() + wait_s
    while True:
        try:
            fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            return f
        except OSError:
            if time.time() >= t_end:
                f.close()
                return None
            sleep(min(poll, max(t_end - time.time(), 0.05)))


def release_lock(f):
    if f is not None:
        import fcntl
        try:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)
        finally:
            f.close()


FALLBACK_JOBS = [("12:30", "waterline_timex_cron.sh"), ("19:25", "waterline_timex_cron.sh"),
                 ("19:45", "cleanup.sh"), ("20:55", "waterline_timex_cron.sh"), ("21:00", "cleanup.sh")]


def station_jobs_due(seconds, local_now):
    """Station jobs starting within `seconds` of the NUC's local time (survey_products.STATION_JOBS_LOCAL)."""
    try:
        import contextlib
        sys.path.insert(0, str(HERE))
        with contextlib.redirect_stdout(io.StringIO()):
            import survey_products
        return survey_products.station_job_overlap(seconds, now=local_now)
    except Exception:
        hits = []
        for day in (0, 1):
            b = (local_now + timedelta(days=day)).replace(second=0, microsecond=0)
            for hm, name in FALLBACK_JOBS:
                t = b.replace(hour=int(hm[:2]), minute=int(hm[3:]))
                if local_now <= t <= local_now + timedelta(seconds=seconds):
                    hits.append(f"{hm} {name}")
        return hits


# --------------------------------------------------------------------------------------------
# the ground file: one pass under the cron's lock, a working copy of the dates needed
# --------------------------------------------------------------------------------------------


def frame_passes(hs_text, max_hs):
    if max_hs is None or hs_text == "":
        return True
    try:
        return float(hs_text) <= max_hs
    except ValueError:
        return True


def scan_ground(path, d_min, d_max, out_csv, max_hs):
    """One pass of the ground file: copies the rows dated [d_min, d_max] to out_csv (tmp + replace) and counts,
    per (date, camera), the frames passing --max-hs, all frames and frames without a wave height; the earliest
    capture per date; the latest georectified date of the WHOLE file; the cameras; the last row's fields."""
    info = {"latest": None, "has_hs": False, "counts": {}, "frames": {}, "blank": {}, "earliest": {},
            "cameras": [], "rows": 0, "last_fields_ok": True, "range": [d_min, d_max]}
    cams, frames = set(), {}
    tmp = Path(str(out_csv) + ".tmp") if out_csv else None
    g = open(tmp, "w", newline="") if tmp else None
    try:
        with open(path, newline="", errors="replace") as f:
            r = csv.reader(f)
            head = next(r, None)
            if not head or "capture_time_utc" not in head or "easting_utm19" not in head:
                raise ValueError("no capture_time_utc / easting_utm19 column")
            w = csv.writer(g) if g else None
            if w:
                w.writerow(head)
            kt, ke, kc = head.index("capture_time_utc"), head.index("easting_utm19"), head.index("camera")
            ks = head.index("source_file") if "source_file" in head else None
            kh = head.index("offshore_hs_m") if "offshore_hs_m" in head else None
            kep = head.index("capture_epoch") if "capture_epoch" in head else None
            info["has_hs"] = kh is not None
            mh = max_hs if kh is not None else None
            nh, last_len, latest = len(head), len(head), ""
            for row in r:
                last_len = len(row)
                if len(row) <= max(kt, ke, kc):
                    continue
                d = row[kt][:10]
                ground = row[ke] != ""
                if ground:
                    cams.add(row[kc])
                    if d > latest:
                        latest = d
                if d_min <= d <= d_max:
                    if w:
                        w.writerow(row)
                    info["rows"] += 1
                    if ground:
                        fr = row[ks] if ks is not None else f"{d}{row[kc]}"
                        if fr not in frames:
                            hs = row[kh] if kh is not None and kh < len(row) else ""
                            frames[fr] = (d, row[kc], frame_passes(hs, mh), hs == "")
                        try:
                            ep = int(float(row[kep])) if kep is not None else \
                                int(parse_iso(row[kt]).timestamp())
                            e0 = info["earliest"].get(d)
                            if e0 is None or ep < e0:
                                info["earliest"][d] = ep
                        except (TypeError, ValueError, AttributeError):
                            pass
            info["last_fields_ok"] = last_len == nh
            info["latest"] = latest or None
    finally:
        if g:
            g.close()
    for fr, (d, c, ok, blank) in frames.items():
        for key, add in (("frames", 1), ("counts", 1 if ok else 0), ("blank", 1 if blank else 0)):
            info[key].setdefault(d, {}).setdefault(c, 0)
            info[key][d][c] += add
    info["cameras"] = sorted(cams)
    if tmp:
        os.replace(str(tmp), str(out_csv))
    return info



def ground_suspect(new, old, cron_log):
    """Why a freshly read ground file looks truncated or half-written (or ''): the last row is short, a
    camera seen before is gone, the passing frames on dates both passes covered dropped by > 20 %, or the
    cron's last run did not end, or logged 'georectification failed'."""
    if not new["last_fields_ok"]:
        return "its last line is short (a write cut off?)"
    if old:
        gone = sorted(set(old.get("cameras") or []) - set(new["cameras"]))
        if gone:
            return f"camera(s) {', '.join(gone)} are gone from it"
        both = [d for d in old.get("counts", {}) if d in new["counts"] or (new["range"][0] <= d <= new["range"][1])]
        both = [d for d in both if old.get("latest") and d < old["latest"]]
        n_old = sum(sum(old["counts"].get(d, {}).values()) for d in both)
        n_new = sum(sum(new["counts"].get(d, {}).values()) for d in both)
        if n_old >= 20 and n_new < 0.8 * n_old:
            return f"its frames on dates read before dropped from {n_old} to {n_new} (> 20 %)"
    p = Path(cron_log) if cron_log else None
    if p is not None and p.exists():
        try:
            lines = p.read_text(errors="replace").splitlines()[-5000:]
        except OSError:
            lines = []
        starts = [i for i, l in enumerate(lines) if "=== run start ===" in l]
        if starts:
            tail = lines[starts[-1]:]
            if any("georectification failed" in l for l in tail):
                return "the cron's last run logged 'georectification failed'"
            if not any("=== run end ===" in l for l in tail):
                return "the cron's last run did not log its end (killed mid-run?)"
    return ""


def usable_dates(gi, d0, d1, cfg, latest=None):
    """Covered dates in [d0, d1] on which a camera has >= STORM_MIN_FRAMES_PER_DAY frames passing --max-hs."""
    latest = latest if latest is not None else gi.get("latest")
    out = []
    if not latest or d0 is None or d1 is None:
        return out
    for d in drange(d0, d1):
        if d >= latest:
            break
        c = gi.get("counts", {}).get(d, {})
        if c and max(c.values()) >= cfg["STORM_MIN_FRAMES_PER_DAY"]:
            out.append(d)
    return out


def covered(gi, d):
    return bool(gi and gi.get("latest") and gi["latest"] > d)


def window_stats(src, start, end, max_hs):
    """Frames per camera, frames without a wave height, and the median setup / Hs / Tp and the water-level range
    of the frames passing --max-hs in a window of a ground file (for the change report)."""
    out = {"frames": {}, "blank_hs": 0, "setup_m": None, "hs_m": None, "tp_s": None, "level_range_m": None,
           "dates": {}}
    try:
        f = open(src, newline="", errors="replace")
    except OSError:
        return out
    seen, setups, hss, tps, lv = {}, [], [], [], []
    with f:
        r = csv.reader(f)
        head = next(r, [])
        k = {c: (head.index(c) if c in head else None) for c in
             ("source_file", "camera", "capture_time_utc", "easting_utm19", "offshore_hs_m", "offshore_tp_s",
              "setup_correction_m", "beach_elevation_navd88", "tide_elevation_navd88")}
        if k["capture_time_utc"] is None or k["easting_utm19"] is None:
            return out
        for row in r:
            if len(row) != len(head) or not row[k["easting_utm19"]]:
                continue
            d = row[k["capture_time_utc"]][:10]
            if not (start <= d <= end):
                continue
            fr = row[k["source_file"]] if k["source_file"] is not None else d
            g = lambda c: row[k[c]] if k[c] is not None else ""   # noqa: E731
            if fr not in seen:
                hs = g("offshore_hs_m")
                ok = frame_passes(hs, max_hs if k["offshore_hs_m"] is not None else None)
                seen[fr] = ok
                cam = g("camera")
                if ok:
                    out["frames"][cam] = out["frames"].get(cam, 0) + 1
                    out["dates"].setdefault(d, 0)
                    out["dates"][d] += 1
                    for lst, c in ((setups, "setup_correction_m"), (hss, "offshore_hs_m"), (tps, "offshore_tp_s")):
                        v = fnum_or_nan(g(c))
                        if np.isfinite(v):
                            lst.append(v)
                out["blank_hs"] += hs == ""
            if seen[fr]:
                v = fnum_or_nan(g("beach_elevation_navd88") or g("tide_elevation_navd88"))
                if np.isfinite(v):
                    lv.append(v)
    for key, lst in (("setup_m", setups), ("hs_m", hss), ("tp_s", tps)):
        out[key] = round(float(np.median(lst)), 3) if lst else None
    out["level_range_m"] = [round(float(np.min(lv)), 2), round(float(np.max(lv)), 2)] if lv else None
    out["frames_total"] = len(seen)
    return out

# --------------------------------------------------------------------------------------------
# the run's context, subprocesses
# --------------------------------------------------------------------------------------------


class SettingsError(Exception):
    pass


def subprocess_runner(argv, log_path=None):
    """Runs a repository script (sys.executable). -> (exit code, stdout+stderr). Output appended to log_path."""
    try:
        p = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True,
                           cwd=str(HERE))
        rc, out = p.returncode, p.stdout or ""
    except OSError as exc:
        rc, out = 1, f"could not run {argv[:2]}: {exc}"
    if log_path:
        try:
            Path(log_path).parent.mkdir(parents=True, exist_ok=True)
            with open(log_path, "a") as f:
                f.write(f"\n$ {' '.join(str(a) for a in argv)}\n{out}\n(exit {rc})\n")
        except OSError:
            pass
    return rc, out


class Ctx:
    def __init__(self, base, now, cfg, cron, mode="live", out_root=None, opener=None, runner=None, mailer=None,
                 sleep=time.sleep, local_now=None, forecast_dir=None, save_raw=None, email=True, fetch=True,
                 coverage_cap=None, ground_src=None, log_path=None, stamp_path=None, msmtp="msmtp",
                 send_waits=(30, 60, 0)):
        self.base = Path(base)
        self.now = utc(now)
        self.cfg, self.cron = cfg, cron
        self.mode = mode
        self.out_root = Path(out_root) if out_root else self.base / "archive/storm_watch"
        self.forecasts_dir = self.base / "archive/storm_watch/forecasts"
        self.events_dir = self.out_root / "events"
        self.state_path = self.out_root / "state.json"
        self.work = self.out_root / "work"
        self.ground_src = Path(ground_src) if ground_src else self.work / "ground_window.csv"
        self.log_path = Path(log_path) if log_path else self.base / "logs/storm_watch.log"
        self.stamp_path = Path(stamp_path) if stamp_path else self.base / "logs/storm_watch_last_run"
        self.opener = opener
        self.runner = runner or subprocess_runner
        self.mailer = mailer
        self.sleep = sleep
        self.local_now = local_now or (lambda: datetime.now())
        self.forecast_dir = forecast_dir
        self.save_raw = save_raw
        self.email = email and cfg["STORM_EMAIL"] == 1 and mode == "live"
        self.fetch = fetch
        self.coverage_cap = coverage_cap
        self.msmtp = msmtp
        self.send_waits = send_waits
        self.notes, self.warnings, self.argvs = [], [], []
        self.has_hs = False
        self.printed = []

    def note(self, text):
        if text not in self.notes:
            self.notes.append(text)

    def warn(self, text):
        if text not in self.warnings:
            self.warnings.append(text)

    def run(self, argv, log_path=None):
        self.argvs.append([str(a) for a in argv])
        return self.runner([str(a) for a in argv], log_path)

    @property
    def writes(self):
        return self.mode in ("live", "replay")

    @property
    def now_h(self):
        return hour_of(self.now)

    @property
    def today(self):
        return self.now.strftime("%Y-%m-%d")


def evdir(ctx, E):
    return ctx.events_dir / E["id"]


# --------------------------------------------------------------------------------------------
# DEM builds (exactly the cron's window-DEM method and options), maps, the change map
# --------------------------------------------------------------------------------------------


def dem_change_mod():
    sys.path.insert(0, str(HERE))
    import dem_change
    return dem_change


def build_dem(ctx, E, stem, start, end, no_plot=False, extra=(), title=None, kind=""):
    """dem_from_contours.py over [start, end] (UTC dates) of the working ground rows with the cron's options.
    -> record: window, stem (relative to the event folder), exit code, cells (>= 5 crossings), range, adequate
    (built and not degenerate for dem_change), frames, its WARNING/NOTE lines."""
    dc = dem_change_mod()
    ed = evdir(ctx, E)
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    argv = ([sys.executable, str(HERE / "dem_from_contours.py"), str(ctx.ground_src), str(stem),
             "--start-date", start, "--end-date", end] + dem_args(ctx.cron, ctx.has_hs) + list(extra)
            + (["--no-plot"] if no_plot else []))
    rc, out = ctx.run(argv, ed / "log.txt")
    rec = {"window": [start, end], "stem": os.path.relpath(str(stem), str(ed)), "rc": rc, "built_utc": iso_z(ctx.now),
           "log_lines": [ln.strip() for ln in out.splitlines() if ln.strip().startswith(("WARNING", "NOTE"))][:20]}
    if rc in (0, 3) and Path(str(stem) + "_dem.asc").exists():
        n, rng, p5, p95 = dc.dem_summary(str(stem), 5)
        why = dc.degenerate(str(stem), 5, 200, 0.5)
        rec.update(cells=n, range=[round(p5, 2), round(p95, 2)] if n else None, adequate=why is None, why=why or "")
        try:
            info = json.loads(Path(str(stem) + "_info.json").read_text())
            rec.update(frames=info.get("frames"), days=info.get("days"), first_date=info.get("first_date"),
                       last_date=info.get("last_date"),
                       camera_days_rejected=(info.get("filters") or {}).get("camera_days_rejected"))
        except (OSError, ValueError):
            pass
        if rc == 3:
            rec["note"] = "grids written, its page not drawn (exit 3)"
        if not no_plot:
            draw_products(ctx, E, stem, title or f"Marconi Beach {kind} DEM, {start} to {end}",
                          f"CACO-05 storm watch {E['id']} {kind} {start}..{end}; intertidal beach only")
    else:
        rec.update(adequate=False, cells=0,
                   why={4: "no waterline points in the window, or the day test rejected every camera-day (exit 4)"}
                   .get(rc, f"dem_from_contours.py failed (exit {rc})"))
    return rec


def draw_products(ctx, E, stem, title, note):
    """The north-up map (dem_map.py) and GeoTIFFs of the folder (asc_to_geotiff.py); neither is fatal."""
    ed = evdir(ctx, E)
    rc, _ = ctx.run([sys.executable, str(HERE / "dem_map.py"), str(stem), "--title", title], ed / "log.txt")
    if rc != 0:
        ctx.note(f"{E['id']}: map of {Path(stem).name} not drawn (dem_map.py exit {rc})")
    rc, _ = ctx.run([sys.executable, str(HERE / "asc_to_geotiff.py"), str(Path(stem).parent), "--note", note],
                    ed / "log.txt")
    if rc != 0:
        ctx.note(f"{E['id']}: GeoTIFFs of {Path(stem).parent.name} not written (exit {rc})")


LOD_RE = re.compile(r"^level of detection\s*:\s*median ([\d.]+) m", re.M)


def change_numbers(stem_out, text):
    """The change map's numbers, recomputed from its _diff/_sig grids (the sums dem_change.py prints)."""
    dc = dem_change_mod()
    diff, h = dc.read_grid(f"{stem_out}_diff.asc")
    sig, _ = dc.read_grid(f"{stem_out}_sig.asc")
    area = h["cellsize"] ** 2
    ok = np.isfinite(diff)
    n = int(ok.sum())
    out = {"cells_compared": n, "area_m2": round(n * area, 1)}
    if n:
        d = diff[ok]
        out.update(median_m=round(float(np.median(d)), 3), mean_m=round(float(np.mean(d)), 3),
                   p10_m=round(float(np.percentile(d, 10)), 3), p90_m=round(float(np.percentile(d, 90)), 3))
    s_ok = np.isfinite(sig)
    with np.errstate(invalid="ignore"):
        ero, acc = s_ok & (sig < 0), s_ok & (sig > 0)
    for key, m in (("erosion", ero), ("accretion", acc)):
        out[key] = {"cells": int(m.sum()), "area_m2": round(float(m.sum() * area), 1),
                    "volume_m3": round(float(np.sum(sig[m]) * area), 1)}
    out["net_m3"] = round(float(np.sum(sig[s_ok]) * area), 1)
    out["significant_share"] = round(float(s_ok.sum()) / n, 3) if n else None
    m = LOD_RE.search(text or "")
    out["lod_median_m"] = float(m.group(1)) if m else None
    out["dem_change_notes"] = [ln.strip() for ln in (text or "").splitlines()
                               if ln.strip().startswith(("NOTE:", "WARNING:"))]
    return out


def dem_drift(stem_a, stem_b):
    """Median (B - A) over the cells both DEMs rest on >= 5 crossings: how far today's processing has moved the
    original pre DEM (M12). -> (median m, cells) or (None, 0)."""
    dc = dem_change_mod()
    try:
        A, _, nA, ha = dc.load_dem(str(stem_a))
        B, _, nB, hb = dc.load_dem(str(stem_b))
        ia, ib, *_ = dc.overlap(ha, hb)
    except (SystemExit, OSError, ValueError):
        return None, 0
    A, nA, B, nB = A[ia], nA[ia], B[ib], nB[ib]
    both = np.isfinite(A) & np.isfinite(B) & (nA >= 5) & (nB >= 5)
    if not both.any():
        return None, 0
    return round(float(np.median((B - A)[both])), 3), int(both.sum())


def pointing_flags(base, dates):
    """{camera: {date: status}} from archive/pointing_daily_<cam>.json for the given dates ('not logged' when
    the log has no entry: pointing_daily.py logs at most 2 days a run)."""
    out = {}
    for cam in ("c1", "c2"):
        p = Path(base) / f"archive/pointing_daily_{cam}.json"
        try:
            days = (json.loads(p.read_text()) or {}).get("days") or {}
        except (OSError, ValueError):
            days = None
        res = {}
        for d in dates:
            if days is None:
                res[d] = "no pointing log"
                continue
            e = days.get(d)
            if not e:
                res[d] = "not logged"
            else:
                st = e.get("status", "?")
                res[d] = st + (f" {e.get('d_deg')}" if st in ("moved", "MOVED") else "")
        out[cam] = res
    return out


def moved(flags, cam, dates):
    return [d for d in dates if str(flags.get(cam, {}).get(d, "")).startswith(("moved", "MOVED"))]


def run_change(ctx, E, part, ref_stem, post_stem, cdir):
    dc = dem_change_mod()
    argv = ([sys.executable, str(HERE / "dem_change.py"), "--a", str(ref_stem), "--b", str(post_stem),
             "--output-dir", str(cdir)] + (["--day-block-lod"] if ctx.cron.get("DEM_DAY_BLOCK_LOD") else []))
    rc, text = ctx.run(argv, evdir(ctx, E) / "log.txt")
    stem_out = Path(cdir) / ("change_" + Path(ref_stem).name.replace("dem_", "") + "_to_"
                             + Path(post_stem).name.replace("dem_", ""))
    del dc
    return rc, text, stem_out


def do_change(ctx, E, part, segs):
    """Pre (rebuilt today, like for like) against post `part`: dem_change.py, change_summary.json."""
    dc = dem_change_mod()
    ed = evdir(ctx, E)
    pre, post = E["pre"], E[part]
    a, b = pre["window"]
    cdir = ed / f"change_{part}"
    if cdir.exists():
        shutil.rmtree(str(cdir), ignore_errors=True)
    cdir.mkdir(parents=True, exist_ok=True)
    ref = cdir / f"pre_{a}_{b}"
    res = {"utc": iso_z(ctx.now), "part": part, "status": "no map", "storm_dates": E.get("storm_dates")}
    rr = build_dem(ctx, E, ref, a, b, no_plot=True, kind="pre (rebuilt)")
    res["pre_rebuilt"] = {k: rr.get(k) for k in ("window", "cells", "range", "frames", "adequate", "why", "log_lines")}
    post_stem = ed / post["stem"]
    if not rr["adequate"]:
        res["why"] = f"the pre DEM rebuilt today from the current ground file is not adequate: {rr['why']}"
    elif not post.get("adequate"):
        res["why"] = f"the post DEM is not adequate: {post.get('why')}"
    else:
        rc, text, stem_out = run_change(ctx, E, part, ref, post_stem, cdir)
        (cdir / "dem_change.txt").write_text(text)
        if rc != 0 or not Path(f"{stem_out}_diff.asc").exists():
            tail = [ln for ln in text.strip().splitlines() if ln.strip()][-1:] or ["?"]
            res["why"] = f"dem_change.py exit {rc}: {tail[0].strip()[:160]}"
        else:
            nums = change_numbers(stem_out, text)
            if not nums["cells_compared"]:
                res["why"] = "no overlapping cells with >= 5 crossings in both DEMs"
                res.update(nums)
            else:
                res.update(nums)
                res["status"] = "map"
                res["stem"] = os.path.relpath(str(stem_out), str(ed))
                res["png"] = res["stem"] + ".png" if Path(f"{stem_out}.png").exists() else None
    # what the numbers rest on, whatever the outcome
    res["pre"] = {"window": pre["window"], "cells": rr.get("cells"), "range": rr.get("range"),
                  "frames": rr.get("frames"), "usable": pre.get("usable"), "record_cells": pre.get("cells")}
    res["post"] = {"window": post["window"], "cells": post.get("cells"), "range": post.get("range"),
                   "frames": post.get("frames"), "usable": post.get("usable")}
    res["rebuilt"] = True
    dm, dn = dem_drift(ed / pre["stem"], ref) if pre.get("stem") and Path(str(ed / pre["stem"]) + "_dem.asc").exists() \
        else (None, 0)
    res["drift_rebuilt_minus_original_pre_m"] = dm
    res["drift_cells"] = dn
    mh = float(ctx.cron["DEM_MAX_HS"]) if (ctx.has_hs and ctx.cron.get("DEM_MAX_HS")) else None
    ws_pre = window_stats(ctx.ground_src, a, b, mh)
    ws_post = window_stats(ctx.ground_src, post["window"][0], post["window"][1], mh)
    res["setup"] = {"pre": {k: ws_pre[k] for k in ("setup_m", "hs_m", "tp_s", "level_range_m")},
                    "post": {k: ws_post[k] for k in ("setup_m", "hs_m", "tp_s", "level_range_m")}}
    if ws_pre["setup_m"] is not None and ws_post["setup_m"] is not None:
        dset = ws_post["setup_m"] - ws_pre["setup_m"]
        sens = 0.25 * abs(dset)
        res["setup"].update(post_minus_pre_m=round(dset, 3), c_error_25pct_m=round(sens, 3),
                            flag=bool(res.get("lod_median_m") and sens > res["lod_median_m"]))
    res["blank_hs_frames"] = {"pre": [ws_pre["blank_hs"], ws_pre.get("frames_total", 0)],
                              "post": [ws_post["blank_hs"], ws_post.get("frames_total", 0)]}
    pre_dates, post_dates = drange(a, b), drange(*post["window"])
    flags = pointing_flags(ctx.base, pre_dates + post_dates)
    res["pointing"] = flags
    res["waves_in_windows"] = {"pre": segs_in(segs, a, b), "post": segs_in(segs, *post["window"])}
    res["build_notes"] = {"pre_rebuilt": rr.get("log_lines", []), "post": post.get("log_lines", [])}
    turned = [c for c in ("c1", "c2") if moved(flags, c, post_dates) and not moved(flags, c, pre_dates)]
    if turned and res["status"] == "map":
        res["per_camera"] = per_camera_change(ctx, E, a, b, post["window"], cdir, turned)
    with open(cdir / "change_summary.json", "w") as f:
        json.dump(res, f, indent=1)
    if res["status"] == "map":
        ctx.run([sys.executable, str(HERE / "asc_to_geotiff.py"), str(cdir), "--note",
                 f"CACO-05 storm watch {E['id']} change pre {a}..{b} to {part} "
                 f"{post['window'][0]}..{post['window'][1]}; intertidal beach only"], ed / "log.txt")
    del dc
    return res


def per_camera_change(ctx, E, a, b, post_window, cdir, turned):
    """A camera reads moved in the post window but not in the pre: the change of each camera on its own."""
    out = {"turned": turned}
    for cam in ("c1", "c2"):
        d = Path(cdir) / "per_camera"
        r1 = build_dem(ctx, E, d / f"pre_{cam}_{a}_{b}", a, b, no_plot=True, extra=["--camera", cam])
        r2 = build_dem(ctx, E, d / f"post_{cam}_{post_window[0]}_{post_window[1]}", post_window[0], post_window[1],
                       no_plot=True, extra=["--camera", cam])
        if not (r1["adequate"] and r2["adequate"]):
            out[cam] = {"why": r1["why"] or r2["why"]}
            continue
        rc, text, stem_out = run_change(ctx, E, "per_camera", evdir(ctx, E) / r1["stem"], evdir(ctx, E) / r2["stem"], d)
        out[cam] = change_numbers(stem_out, text) if rc == 0 and Path(f"{stem_out}_diff.asc").exists() \
            else {"why": f"dem_change.py exit {rc}"}
    return out


def segs_in(segs, a, b):
    """Wave events above the calm level with rough hours on the dates [a, b]."""
    h0, h1 = date_hour0(a), date_hour0(b) + 23
    return [seg_text(s) for s in segs or [] if s["start"] <= h1 and s["last"] >= h0]


# --------------------------------------------------------------------------------------------
# the state file
# --------------------------------------------------------------------------------------------


def fresh_state(now):
    return {"version": STATE_VERSION, "station": STATION, "created_utc": iso_z(now), "last_run": {},
            "nws_points": None, "ground_cache": None, "events": [], "seen_segments": [], "install_utc": iso_z(now)}


def load_state(ctx):
    """-> state. Missing: a fresh one (segments that ended before now are 'seen': no email for old storms).
    Unreadable: moved aside, the previous copy (state.json.prev) restored, else rebuilt from events/*/event.json,
    with a WARNING. An unknown version: SettingsError (exit 1)."""
    p = ctx.state_path
    if not p.exists():
        return fresh_state(ctx.now), True
    try:
        st = json.loads(p.read_text())
        if not isinstance(st, dict):
            raise ValueError("not an object")
    except (OSError, ValueError) as exc:
        bad = Path(f"{p}.bad-{ctx.now:%Y%m%dT%H%MZ}")
        if ctx.writes:
            try:
                os.replace(str(p), str(bad))
            except OSError:
                pass
        ctx.warn(f"WARNING: state file unreadable ({str(exc)[:60]}); moved aside to {bad.name}")
        prev = Path(str(p) + ".prev")
        try:
            st = json.loads(prev.read_text())
            ctx.warn("WARNING: restored the previous state (state.json.prev)")
        except (OSError, ValueError):
            st = fresh_state(ctx.now)
            for ej in sorted(ctx.events_dir.glob("*/event.json")):
                try:
                    st["events"].append(json.loads(ej.read_text()))
                except (OSError, ValueError):
                    pass
            ctx.warn(f"WARNING: state rebuilt from {len(st['events'])} events/*/event.json")
        return st, False
    if st.get("version") != STATE_VERSION:
        raise SettingsError(f"{p} has state version {st.get('version')!r}; this storm_watch.py reads version "
                            f"{STATE_VERSION}")
    return st, False


def save_state(ctx, state):
    if not ctx.writes:
        return
    p = ctx.state_path
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(p) + ".tmp")
    with open(tmp, "w") as f:
        json.dump(state, f, indent=1, default=str)
        f.flush()
        os.fsync(f.fileno())
    if p.exists():
        try:
            shutil.copy2(str(p), str(p) + ".prev")
        except OSError:
            pass
    os.replace(str(tmp), str(p))
    for E in state.get("events", []):
        d = evdir(ctx, E)
        if d.exists():
            t = d / "event.json.tmp"
            t.write_text(json.dumps(E, indent=1, default=str))
            os.replace(str(t), str(d / "event.json"))

# --------------------------------------------------------------------------------------------
# events: dates, windows and the state machine
# --------------------------------------------------------------------------------------------


def new_event(ctx, state, kind, onset_d):
    ids = {e["id"] for e in state["events"]}
    eid, k = f"{onset_d}_{kind}", 2
    while eid in ids or (ctx.writes and (ctx.events_dir / eid).exists()):
        eid, k = f"{onset_d}_{kind}_{k}", k + 1
    E = {"id": eid, "kind": kind, "hindcast": ctx.mode == "replay", "status": "watch", "closed_reason": None,
         "opened_utc": iso_z(ctx.now), "updated_utc": iso_z(ctx.now), "forecast": None,
         "observed": {"onset": None, "peak_hs": None, "peak_time": None, "end": None, "confirmed_utc": None,
                      "sources": [], "segments": [], "cut_segment": None, "cut_date": None, "cut_state": None},
         "storm_dates": None, "first_post": None, "first_post_exact": False,
         "pre": {"status": "planned", "planned": None, "window": None}, "pre_provisional": None,
         "post3": None, "post7": None, "change": {}, "pending": [], "emails": [], "notices": {},
         "prev_event": None, "merged_from": [], "notes": [], "history": [], "waiting": {}}
    state["events"].append(E)
    if ctx.writes:
        evdir(ctx, E).mkdir(parents=True, exist_ok=True)
    return E


def hist(ctx, E, text):
    E["history"].append(f"{ctx.today}: {text}")
    E["history"] = E["history"][-120:]
    E["updated_utc"] = iso_z(ctx.now)
    ctx.touched.setdefault(E["id"], []).append(text)


def queue_email(ctx, E, kind, **data):
    E["emails"].append({"kind": kind, "utc": iso_z(ctx.now), "status": "pending", "tries": 0, "data": data,
                        "subject": None, "body_file": None})


def onset_date(E):
    if E["observed"]["onset"]:
        return E["observed"]["onset"][:10]
    if E["forecast"] and E["forecast"].get("onset_fc"):
        return E["forecast"]["onset_fc"][:10]
    return None


def is_open(E):
    return E["status"] != "closed"


def spans_of(E):
    return [(iso_hour(x["start"]), iso_hour(x["end"]) if x.get("end") else None)
            for x in E["observed"]["segments"]]


def overlaps(a0, a1, b0, b1):
    return a0 <= b1 and a1 >= b0


def seg_overlaps_spans(s, spans, now_h):
    a, b = seg_span(s, now_h)
    return any(overlaps(a, b, x0, x1 if x1 is not None else now_h) for x0, x1 in spans)


def owned_by(state, s, now_h, exclude=None):
    """The events (not false alarms or unverified) owning a segment overlapping s, oldest first."""
    out = []
    for E in state["events"]:
        if E is exclude or E["closed_reason"] in ("false_alarm", "unverified"):
            continue
        if seg_overlaps_spans(s, spans_of(E), now_h):
            out.append(E)
    return out


def claimed(state, s, now_h):
    if owned_by(state, s, now_h):
        return True
    for E in state["events"]:
        c = E["observed"].get("cut_segment")
        if E["status"] == "post" and c and E["observed"].get("cut_state") == "pending":
            if seg_overlaps_spans(s, [(iso_hour(c["start"]), iso_hour(c["end"]) if c.get("end") else None)], now_h):
                return True
    seen = [(iso_hour(x["start"]), iso_hour(x["end"]) if x.get("end") else None) for x in state.get("seen_segments", [])]
    return seg_overlaps_spans(s, seen, now_h)


def span_rec(s):
    return {"start": hour_iso(s["start"]), "end": hour_iso(s["end"]) if s["end"] is not None else None,
            "peak": s["peak"]}


def typical_first_frame(gi):
    """Seconds after UTC midnight of the first frame of a day (median over the dates read), else 12:00."""
    e = (gi or {}).get("earliest") or {}
    secs = [ep - date_hour0(d) * 3600 for d, ep in e.items()]
    return float(np.median(secs)) if secs else 12 * 3600.0


def set_first_post(ctx, E, gi):
    """The first post-storm date: the end's date when no ground frame of that date is earlier than the end,
    else the next date (provisional, from the usual first-frame time, until the ground file covers the date);
    never the onset date. The storm dates run from the onset date to the date before."""
    end_h = iso_hour(E["observed"]["end"])
    d = hour_date(end_h)
    if gi and covered(gi, d) and d in (gi.get("earliest") or {}):
        first, exact = (d if gi["earliest"][d] >= end_h * 3600 else dadd(d, 1)), True
    else:
        sec = end_h * 3600 - date_hour0(d) * 3600
        first, exact = (d if sec <= typical_first_frame(gi) else dadd(d, 1)), False
    od = onset_date(E)
    first = max(first, dadd(od, 1))
    E["first_post"], E["first_post_exact"] = first, exact
    E["storm_dates"] = [od, dadd(first, -1)]


def floor_date(state, E, od):
    best = None
    for X in state["events"]:
        if X is E or not X.get("storm_dates") or X["closed_reason"] == "merged":
            continue
        last = X["storm_dates"][1]
        if last < od:
            best = dadd(last, 1) if best is None else max(best, dadd(last, 1))
    return best


def breaking_dates(segs, cfg, E, now_h):
    """[(first date, last date, segment)] of wave events that cut windows, not owned by E."""
    own = spans_of(E)
    return [(hour_date(s["start"]), hour_date(s["last"]), s) for s in segs
            if breaks(s, cfg) and not seg_overlaps_spans(s, own, now_h)]


def plan_pre(ctx, state, E, segs, gi, od, end=None, extend=True):
    """The pre window for onset date od: the PRE_DAYS dates before it (after the last earlier storm), started
    after a wave event inside it when 3 usable dates remain, extended back (up to PRE_MAX_DAYS, never into a
    storm or a wave event) while fewer than MIN_USABLE dates are usable. -> (start, end, flags)."""
    cfg = ctx.cfg
    end = end or dadd(od, -1)
    floor = floor_date(state, E, od)
    start = dadd(od, -cfg["STORM_PRE_DAYS"])
    if floor and floor > start:
        start = floor
    flags = []
    br = breaking_dates(segs, cfg, E, ctx.now_h)
    for d0, d1, s in sorted(br, key=lambda x: x[0]):
        if d1 < start or d0 > end:
            continue
        after = dadd(d1, 1)
        if gi and after <= end and len(usable_dates(gi, after, end, cfg)) >= cfg["STORM_MIN_USABLE_DAYS"]:
            start = max(start, after)
        else:
            flags.append(f"a wave event ({seg_text(s)}) lies inside the pre window: the pre DEM pools the beach "
                         "before and after it")
    if extend and gi:
        lo = dadd(od, -cfg["STORM_PRE_MAX_DAYS"])
        if floor and floor > lo:
            lo = floor
        while len(usable_dates(gi, start, end, cfg)) < cfg["STORM_MIN_USABLE_DAYS"] and start > lo:
            prev = dadd(start, -1)
            if any(d0 <= prev <= d1 for d0, d1, _ in br):
                break
            start = prev
    return start, end, flags


def planned_pre_fc(ctx, state, E, segs):
    """A forecast event's planned pre window, ended before any observed rough spell that started near the
    forecast onset (a storm that came early, or a sea that never confirmed)."""
    od = onset_date(E)
    start, end, _ = plan_pre(ctx, state, E, segs, None, od, extend=False)
    on = iso_hour(E["forecast"]["onset_fc"])
    for s in segs:
        if s["start"] >= date_hour0(start) and s["start"] <= on + 48 + 12:
            end = min(end, dadd(hour_date(s["start"]), -1))
    return [start, end]


def wave_gaps(obs, a, b, cfg):
    """Dates of [a, b] with fewer than STORM_MIN_WAVE_HOURS observed hours (only dates the record reaches)."""
    out = []
    if obs is None:
        return out
    v = np.isfinite(obs.val)
    for d in drange(a, b):
        h0 = date_hour0(d) - obs.h0
        if h0 + 24 <= 0 or h0 >= obs.n:
            continue
        if int(v[max(h0, 0):min(h0 + 24, obs.n)].sum()) < cfg["STORM_MIN_WAVE_HOURS"]:
            out.append(d)
    return out


def waves_ready(ctx, E, key, obs, a, b):
    """True when the wave record covers [a, b] (or has been waited for 3 days: built with a flag)."""
    gaps = wave_gaps(obs, a, b, ctx.cfg)
    if not gaps:
        E["waiting"].pop("waves_" + key, None)
        return True, ""
    since = E["waiting"].setdefault("waves_" + key, ctx.today)
    if ddiff(ctx.today, since) >= 3:
        return True, f"the wave record has gaps on {', '.join(gaps)}: rough frames there could not be left out"
    return False, f"waiting for the wave record on {', '.join(gaps)}"


def part_title(E, kind, a, b):
    sd = E.get("storm_dates")
    if kind.startswith("pre"):
        return (f"Marconi Beach PRE-storm DEM, {a} to {b} (UTC dates)"
                + (f", before the {sd[0]}..{sd[1]} storm" if sd else
                   f", before the storm forecast for {onset_date(E)}"))
    which = "first 3 usable days" if kind == "post3" else "7 days"
    return f"Marconi Beach POST-storm DEM ({which}), {a} to {b}"


def try_pre(ctx, state, E, segs, gi, obs, provisional=False, window=None):
    """Builds the final pre DEM (or the provisional one). -> True when a decision was made this run."""
    cfg = ctx.cfg
    od = onset_date(E)
    ed = evdir(ctx, E)
    if provisional:
        latest = gi["latest"]
        planned = E["pre"]["planned"]
        end = min(dadd(latest, -1), planned[1])
        start = dadd(end, -(cfg["STORM_PRE_DAYS"] - 1))
        if end < start or not covered(gi, end):
            E["pre_provisional"] = {"status": "none", "why": "the ground file does not reach the window yet"}
            return True
        rec = build_dem(ctx, E, ed / "pre_provisional" / f"pre_{start}_{end}", start, end,
                        title=part_title(E, "pre", start, end) + " (PROVISIONAL)", kind="pre (provisional)")
        rec["usable"] = usable_dates(gi, start, end, cfg)
        rec["status"] = "built" if rec["adequate"] else "inadequate"
        E["pre_provisional"] = rec
        hist(ctx, E, f"provisional pre {start}..{end}: {rec['cells']} cells" + ("" if rec["adequate"] else
                                                                                  f" ({rec['why']})"))
        return True
    if window is not None:
        start, end = window
        flags = []
        if gi:
            lo = dadd(od, -cfg["STORM_PRE_MAX_DAYS"])
            while len(usable_dates(gi, start, end, cfg)) < cfg["STORM_MIN_USABLE_DAYS"] and start > lo:
                start = dadd(start, -1)
    else:
        start, end, flags = plan_pre(ctx, state, E, segs, gi, od)
    if not covered(gi, end):
        return False
    u = usable_dates(gi, start, end, cfg)
    pre = E["pre"]
    pre.update(window=[start, end], usable=u, flags=flags)
    if ctx.cron.get("DEM_ENABLE") != "1":
        pre.update(status="disabled", why=f"DEM stage disabled in the cron (DEM_ENABLE={ctx.cron.get('DEM_ENABLE')})")
        hist(ctx, E, "pre: DEM stage disabled")
        return True
    if len(u) < cfg["STORM_MIN_USABLE_DAYS"]:
        pre.update(status="inadequate", why=f"{len(u)} usable date(s) in {start}..{end} (fewer than "
                                            f"{cfg['STORM_MIN_USABLE_DAYS']}; the cameras, GNSS-R or calm "
                                            "frames were missing)", cells=0, built_utc=iso_z(ctx.now))
        hist(ctx, E, f"pre INADEQUATE: {pre['why']}")
        return True
    ok, why = waves_ready(ctx, E, "pre", obs, start, end)
    if not ok:
        pre["waiting_why"] = why
        return False
    if why:
        flags.append(why)
    lo = dadd(od, -cfg["STORM_PRE_MAX_DAYS"])
    fl = floor_date(state, E, od)
    if fl and fl > lo:
        lo = fl
    for attempt in range(4):
        rec = build_dem(ctx, E, ed / "pre" / f"pre_{start}_{end}", start, end,
                        title=part_title(E, "pre", start, end), kind="pre")
        if rec["adequate"] or attempt == 3 or start <= lo:
            break
        hist(ctx, E, f"pre {start}..{end} not adequate ({rec['why']}): one more day back")
        start = dadd(start, -1)
    rec.update(usable=usable_dates(gi, start, end, cfg), flags=flags,
               status="final" if rec["adequate"] else "inadequate")
    if not rec["adequate"]:
        rec["why"] = rec.get("why") or "not adequate"
    E["pre"].update(rec)
    hist(ctx, E, f"pre {start}..{end} {rec['status']}: {rec.get('cells')} cells"
         + ("" if rec["adequate"] else f" ({rec['why']})"))
    return True


def post_lim(ctx, E):
    lim = dadd(E["first_post"], ctx.cfg["STORM_POST_MAX_DAYS"] - 1)
    cd = E["observed"].get("cut_date")
    if cd and E["observed"].get("cut_state") in ("pending", "cut"):
        lim = min(lim, cd)
    return lim


def try_post3(ctx, state, E, gi, obs, segs):
    cfg = ctx.cfg
    fp = E["first_post"]
    lim = post_lim(ctx, E)
    u = usable_dates(gi, fp, lim, cfg)
    need = cfg["STORM_POST1_USABLE"]
    if len(u) < need:
        if covered(gi, lim) and E["observed"].get("cut_state") != "pending":
            E["post3"] = {"status": "none", "why": f"{len(u)} usable date(s) in {fp}..{lim}", "window": None}
            hist(ctx, E, f"post3 none: {E['post3']['why']}")
            return True
        return False
    d3 = u[need - 1]
    ok, why = waves_ready(ctx, E, "post3", obs, fp, d3)
    if not ok:
        return False
    if ctx.cron.get("DEM_ENABLE") != "1":
        E["post3"] = {"status": "disabled", "window": [fp, d3], "usable": u[:need]}
        hist(ctx, E, "post3: DEM stage disabled")
        return True
    rec = build_dem(ctx, E, evdir(ctx, E) / "post3" / f"post3_{fp}_{d3}", fp, d3,
                    title=part_title(E, "post3", fp, d3), kind="post3")
    rec.update(usable=u[:need], status="built" if rec["adequate"] else "inadequate", flags=[why] if why else [])
    E["post3"] = rec
    hist(ctx, E, f"post3 {fp}..{d3}: {rec.get('cells')} cells" + ("" if rec["adequate"] else f" ({rec['why']})"))
    return True


def try_post7(ctx, state, E, gi, obs):
    """-> True when post7 was decided (built, skipped or none) this run."""
    cfg = ctx.cfg
    fp = E["first_post"]
    ob = E["observed"]
    cap = dadd(fp, cfg["STORM_POST_DAYS"] - 1)
    if ob.get("cut_state") == "pending" and ob.get("cut_date") and ob["cut_date"] <= cap:
        return False
    cut = ob.get("cut_state") == "cut" and ob.get("cut_date") and ob["cut_date"] <= cap
    if cut:
        cap = ob["cut_date"]
    p3 = E["post3"] or {}
    if p3.get("status") == "none":
        E["post7"] = {"status": "none", "why": "no first post DEM: " + p3.get("why", ""), "window": None}
        return True
    if not covered(gi, cap):
        return False
    u = usable_dates(gi, fp, cap, cfg)
    reason = ""
    if len(u) < cfg["STORM_MIN_USABLE_DAYS"] and not cut:
        last = dadd(fp, cfg["STORM_POST_MAX_DAYS"] - 1)
        while len(u) < cfg["STORM_MIN_USABLE_DAYS"] and cap < last and covered(gi, dadd(cap, 1)):
            cap = dadd(cap, 1)
            u = usable_dates(gi, fp, cap, cfg)
        if len(u) < cfg["STORM_MIN_USABLE_DAYS"]:
            if cap >= last:
                E["post7"] = {"status": "none", "why": f"{len(u)} usable date(s) in {fp}..{cap}", "window": None}
                return True
            return False
        reason = f"extended to {cap} for {cfg['STORM_MIN_USABLE_DAYS']} usable dates"
    if p3.get("window") == [fp, cap]:
        why = (f"cut by the next storm (wave event from {ob['cut_date'] and dadd(ob['cut_date'], 1)})" if cut else
               f"the 7-day window holds only the usable dates of the first post DEM ({len(u)})")
        E["post7"] = {"status": "same as post3", "why": why, "window": [fp, cap]}
        return True
    ok, wwhy = waves_ready(ctx, E, "post7", obs, fp, cap)
    if not ok:
        return False
    if ctx.cron.get("DEM_ENABLE") != "1":
        E["post7"] = {"status": "disabled", "window": [fp, cap], "usable": u}
        return True
    rec = build_dem(ctx, E, evdir(ctx, E) / "post7" / f"post7_{fp}_{cap}", fp, cap,
                    title=part_title(E, "post7", fp, cap), kind="post7")
    flags = [f for f in (wwhy, reason, ("cut short by the next storm (from %s)" % dadd(ob["cut_date"], 1)) if cut
                         else "") if f]
    rec.update(usable=u, status="built" if rec["adequate"] else "inadequate", flags=flags)
    E["post7"] = rec
    hist(ctx, E, f"post7 {fp}..{cap}: {rec.get('cells')} cells" + ("" if rec["adequate"] else f" ({rec['why']})"))
    return True


def change_for(ctx, E, part, segs):
    pre = E["pre"]
    if pre.get("status") not in ("final", "inadequate", "disabled"):
        return False
    rec = E[part] or {}
    if pre.get("status") != "final" or rec.get("status") != "built":
        why = (f"no pre-storm DEM: {pre.get('why')}" if pre.get("status") != "final"
               else f"no {part} DEM: {rec.get('why') or rec.get('status')}")
        E["change"][part] = {"status": "no map", "why": why}
    else:
        E["change"][part] = do_change(ctx, E, part, segs)
    c = E["change"][part]
    hist(ctx, E, f"change {part}: " + (f"{c.get('cells_compared')} cells, median {c.get('median_m'):+.2f} m, erosion "
                                       f"{c['erosion']['volume_m3']:+.0f} m3" if c.get("status") == "map"
                                       else f"no map ({c.get('why')})"))
    queue_email(ctx, E, f"report_{part}")
    return True


def merge_or_cut(ctx, state, E, segs, gi):
    cfg = ctx.cfg
    ob = E["observed"]
    u = usable_dates(gi, E["first_post"], ob["cut_date"], cfg)
    cs = ob["cut_segment"]
    if len(u) < cfg["STORM_MIN_USABLE_DAYS"]:
        ob["segments"].append(dict(cs))
        merged = []
        span = [(iso_hour(cs["start"]), iso_hour(cs["end"]) if cs.get("end") else None)]
        for X in state["events"]:
            if X is E or X["status"] == "closed" and X["closed_reason"] in ("merged", "complete"):
                continue
            if any(overlaps(a, b if b is not None else ctx.now_h, span[0][0],
                            span[0][1] if span[0][1] is not None else ctx.now_h) for a, b in spans_of(X)):
                X["status"], X["closed_reason"] = "closed", "merged"
                X["notes"].append(f"merged into {E['id']} (one compound event)")
                hist(ctx, X, f"closed: merged into {E['id']}")
                merged.append(X["id"])
        E["merged_from"] += merged
        E["notes"].append(f"COMPOUND: the next storm ({cs['start'][:13]}Z) came before 3 usable post-storm dates "
                          f"({len(u)}: {', '.join(u) or 'none'}): one event, its storm dates extended")
        ob.update(end=None, cut_segment=None, cut_date=None, cut_state=None)
        E.update(status="storm", first_post=None, first_post_exact=False, storm_dates=None, post3=None, post7=None)
        E["waiting"] = {k: v for k, v in E["waiting"].items() if not k.startswith(("waves_post", "post"))}
        hist(ctx, E, f"MERGED with the storm from {cs['start'][:13]}Z ({len(u)} usable post date(s))")
    else:
        ob["cut_state"] = "cut"
        hist(ctx, E, f"post window CUT at {ob['cut_date']} by the storm from {cs['start'][:13]}Z ({len(u)} usable)")


def waits(ctx, state, E, gi):
    """A date the ground file should have covered by now: a note after 7 days, give up after MAX_WAIT."""
    cfg = ctx.cfg
    need = None
    if E["status"] in ("storm", "post") and E["pre"].get("status") not in ("final", "inadequate", "disabled"):
        need = ("pre", E["pre"].get("window", [None, None])[1] if E["pre"].get("window") else dadd(onset_date(E), -1))
    elif E["status"] == "post" and E["post3"] is None:
        need = ("post3", dadd(E["first_post"], cfg["STORM_POST1_USABLE"] - 1))
    elif E["status"] == "post" and E["post7"] is None:
        need = ("post7", dadd(E["first_post"], cfg["STORM_POST_DAYS"] - 1))
    if not need or covered(gi, need[1]):
        return
    late = ddiff(ctx.today, need[1])
    latest = (gi or {}).get("latest") or "?"
    if late > 7 and not E["notices"].get("waiting_ground"):
        E["notices"]["waiting_ground"] = ctx.today
        queue_email(ctx, E, "waiting_ground", part=need[0], date=need[1], latest=latest)
    if late > cfg["STORM_MAX_WAIT_DAYS"]:
        E["status"], E["closed_reason"] = "closed", "gave_up"
        hist(ctx, E, f"closed: gave up waiting for the ground file to reach {need[1]} (it ends {latest})")
        queue_email(ctx, E, "gave_up", part=need[0], date=need[1], latest=latest)


# ---- updates that need no ground ------------------------------------------------------------


def obs_coverage(obs, h0, h1):
    if obs is None or h1 < h0:
        return 0.0
    a, b = max(h0 - obs.h0, 0), min(h1 - obs.h0, obs.n - 1)
    if b < a:
        return 0.0
    return float(np.isfinite(obs.val[a:b + 1]).sum()) / (h1 - h0 + 1)


def adopt_segment(ctx, E, s):
    ob = E["observed"]
    if not ob["onset"]:
        ob["onset"] = hour_iso(s["start"])
    if not seg_overlaps_spans(s, spans_of(E), ctx.now_h):
        ob["segments"].append(span_rec(s))
    follow_segments(ctx, E, [s])


def follow_segments(ctx, E, segs_now):
    """Refresh the owned segments from today's series: spans, peak, end (the last one's)."""
    ob = E["observed"]
    mine = [s for s in segs_now if seg_overlaps_spans(s, spans_of(E), ctx.now_h)]
    if not mine:
        return None
    recs = []
    for x in ob["segments"]:
        hit = [s for s in mine if seg_overlaps_spans(s, [(iso_hour(x["start"]), iso_hour(x["end"]) if x.get("end")
                                                         else None)], ctx.now_h)]
        recs.append(span_rec(hit[0]) if hit else x)
    ob["segments"] = recs
    for s in mine:
        if s["peak"] is not None and (ob["peak_hs"] is None or s["peak"] > ob["peak_hs"]):
            ob["peak_hs"], ob["peak_time"] = s["peak"], hour_iso(s["peak_h"])
        if s["raw_max"] and (ob.get("raw_max") is None or s["raw_max"] > ob["raw_max"]):
            ob["raw_max"] = s["raw_max"]
        ob["sources"] = sorted(set(ob["sources"]) | set(s["srcs"]))
    last = max(mine, key=lambda s: s["start"])
    return last


def update_watch(ctx, state, E, segs, T, t_hit, obs):
    fc = E["forecast"]
    on, en = iso_hour(fc["onset_fc"]), iso_hour(fc["end_fc"])
    for s in sorted(segs, key=lambda x: x["start"]):
        a, b = seg_span(s, ctx.now_h)
        if overlaps(a, b, on - 24, en + 48) and confirms(s, ctx.cfg) and not owned_by(state, s, ctx.now_h, exclude=E):
            adopt_segment(ctx, E, s)
            E["observed"]["confirmed_utc"] = iso_z(ctx.now)
            E["status"] = "storm"
            hist(ctx, E, f"CONFIRMED by the observed waves ({seg_text(s)})")
            update_storm(ctx, state, E, segs)
            return
    cov = obs_coverage(obs, on - 24, min(en + 48, ctx.now_h))
    if ctx.now_h > en + 48 and not t_hit and cov >= 0.5:
        close_false_alarm(ctx, E, segs, obs, "the observed waves never reached the confirm level")
    elif ctx.now_h > en + 7 * 24 and cov < 0.5:
        E["status"], E["closed_reason"] = "closed", "unverified"
        E["pending"].append({"what": "pre_final", "window": E["pre"]["planned"], "since": ctx.today})
        hist(ctx, E, f"closed UNVERIFIED: observed waves cover {100 * cov:.0f} % of the storm's window")
        queue_email(ctx, E, "unverified", coverage=round(cov, 2))
    elif ddiff(ctx.today, E["opened_utc"][:10]) > 10 and not (t_hit and iso_hour(T["onset_fc"]) <= ctx.now_h + 24):
        close_false_alarm(ctx, E, segs, obs, "the forecast onset kept receding for more than 10 days")


def close_false_alarm(ctx, E, segs, obs, why):
    fc = E["forecast"]
    on, en = iso_hour(fc["onset_fc"]), iso_hour(fc["end_fc"])
    peak = None
    if obs is not None:
        a, b = max(on - 24 - obs.h0, 0), min(en + 48 - obs.h0, obs.n - 1)
        if b >= a and np.isfinite(obs.prim[a:b + 1]).any():
            peak = round(float(np.nanmax(obs.prim[a:b + 1])), 2)
    E["status"], E["closed_reason"] = "closed", "false_alarm"
    E["observed"]["peak_hs"] = peak
    E["pending"].append({"what": "pre_final", "window": E["pre"]["planned"], "since": ctx.today})
    hist(ctx, E, f"closed FALSE ALARM: {why} (observed Marconi peak {fnum(peak)} m)")
    queue_email(ctx, E, "false_alarm", why=why, peak=peak)


def update_storm(ctx, state, E, segs):
    last = follow_segments(ctx, E, segs)
    ob = E["observed"]
    if last is not None and last["end"] is not None and ob["end"] is None:
        ob["end"] = hour_iso(last["end"])
        E["status"] = "post"
        set_first_post(ctx, E, state.get("ground_cache"))
        hist(ctx, E, f"storm over: end {ob['end'][:13]}Z, storm dates {E['storm_dates'][0]}..{E['storm_dates'][1]}, "
                     f"first post date {E['first_post']}")
        return
    if ob["end"] is None and ob.get("peak_time"):
        ph = iso_hour(ob["peak_time"])
        if ctx.now_h > ph + 7 * 24 and not E["notices"].get("not_over"):
            E["notices"]["not_over"] = ctx.today
            queue_email(ctx, E, "not_over")
        if ctx.now_h > ph + 14 * 24:
            E["status"], E["closed_reason"] = "closed", "gave_up"
            hist(ctx, E, "closed: no end of the storm 14 days after its peak")
            queue_email(ctx, E, "gave_up", part="end", date=None, latest=None)


def update_post(ctx, state, E, segs):
    follow_segments(ctx, E, segs)
    ob = E["observed"]
    if ob.get("cut_segment"):
        return
    end_h = iso_hour(ob["end"])
    cap = dadd(E["first_post"], ctx.cfg["STORM_POST_DAYS"] - 1)
    for s in sorted(segs, key=lambda x: x["start"]):
        if s["start"] <= end_h or seg_overlaps_spans(s, spans_of(E), ctx.now_h):
            continue
        if E["first_post"] <= hour_date(s["start"]) <= cap and breaks(s, ctx.cfg):
            ob.update(cut_segment=span_rec(s), cut_date=dadd(hour_date(s["start"]), -1), cut_state="pending")
            hist(ctx, E, f"a wave event inside the post window ({seg_text(s)}): post dates end {ob['cut_date']} "
                         "unless merged")
            return


# ---- the forecast's events ---------------------------------------------------------------------


def event_refs(ctx, E):
    """(onset_ref, end_ref) hours of an open event, for matching a new trigger."""
    ob, fc = E["observed"], E["forecast"] or {}
    if E["status"] == "watch":
        return iso_hour(fc["onset_fc"]), iso_hour(fc["end_fc"])
    on = iso_hour(ob["onset"]) if ob["onset"] else iso_hour(fc.get("onset_fc"))
    if E["status"] == "storm":
        return on, iso_hour(ob["end"]) if ob["end"] else ctx.now_h + 24
    return on, iso_hour(ob["end"])


def trigger_match(ctx, state, T):
    if not T["triggered"]:
        return []
    on, en = iso_hour(T["onset_fc"]), iso_hour(T["end_fc"])
    out = []
    for E in state["events"]:
        if not is_open(E) or E["status"] not in ("watch", "storm", "post"):
            continue
        a, b = event_refs(ctx, E)
        if a is None or b is None:
            continue
        if on <= b + 72 and en >= a - 48:
            out.append(E)
    return out


def fc_record(ctx, fc, T):
    return {"onset_fc": T["onset_fc"], "end_fc": T["end_fc"], "rules": T["rules"], "by": T["by"],
            "peak": T["peak"], "alerts": T["alerts"], "first_trigger_hour": T["first_trigger_hour"],
            "hours_over": T["hours_over"], "outlook": T["outlook"], "om_used": fc.get("om_used", False),
            "sources": fc.get("status", {}), "text": fc.get("text", []), "hindcast": bool(fc.get("hindcast")),
            "wind_peak": wind_peak(fc, ctx.now, ctx.cfg)}


def wind_peak(fc, now, cfg):
    s = fc.get("wind")
    if s is None:
        return None
    h0 = hour_of(now)
    w = (s["h"] >= h0) & (s["h"] <= h0 + cfg["STORM_LEAD_HOURS"])
    if not np.isfinite(s["ws"][w]).any():
        return None
    k = int(np.nanargmax(np.where(w, s["ws"], np.nan)))
    g = s["gust"][w]
    return {"ws": round(float(s["ws"][k]), 1), "time": hour_iso(int(s["h"][k])),
            "dir": None if not np.isfinite(s["wdir"][k]) else round(float(s["wdir"][k])),
            "gust": round(float(np.nanmax(g)), 1) if np.isfinite(g).any() else None, "source": s["name"]}


def record_run(E, ctx, T):
    f = E["forecast"]
    f["runs"] = (f.get("runs") or [])[-19:] + [{"utc": iso_z(ctx.now), "triggered": T["triggered"],
                                                 "peak_hs": (T["peak"] or {}).get("hs"),
                                                 "peak_time": (T["peak"] or {}).get("time"), "rules": T["rules"]}]


def handle_trigger(ctx, state, fc, T, matches, segs):
    if not T["triggered"]:
        for E in state["events"]:
            if is_open(E) and E["forecast"]:
                record_run(E, ctx, T)
        return None
    live = [E for E in matches if E["status"] in ("watch", "storm")]
    if live:
        E = max(live, key=lambda e: e["opened_utc"])
        old = E["forecast"] or {}
        new = fc_record(ctx, fc, T)
        why = escalation(old, new)
        new.update(first_trigger_utc=old.get("first_trigger_utc") or iso_z(ctx.now), last_trigger_utc=iso_z(ctx.now),
                   runs=old.get("runs") or [])
        E["forecast"] = new
        record_run(E, ctx, T)
        if E["status"] == "watch":
            E["pre"]["planned"] = planned_pre_fc(ctx, state, E, segs)
        if why:
            sent = [m for m in E["emails"] if m["kind"] == "update" and m["utc"][:10] == ctx.today]
            if not sent:
                queue_email(ctx, E, "update", why=why)
                hist(ctx, E, "forecast escalated: " + "; ".join(why))
        return E
    posts = [E for E in matches if E["status"] == "post"]
    if posts and all(iso_hour(T["onset_fc"]) <= iso_hour(E["observed"]["end"]) for E in posts):
        for E in posts:
            E["notes"].append(f"{ctx.today}: forecast still triggering (onset ~{T['onset_fc'][:13]}Z) after the storm")
        return None
    E = new_event(ctx, state, "fc", T["onset_fc"][:10])
    E["forecast"] = fc_record(ctx, fc, T)
    E["forecast"].update(first_trigger_utc=iso_z(ctx.now), last_trigger_utc=iso_z(ctx.now), runs=[])
    record_run(E, ctx, T)
    if posts:
        E["prev_event"] = max(posts, key=lambda e: e["opened_utc"])["id"]
    E["pre"]["planned"] = planned_pre_fc(ctx, state, E, segs)
    hist(ctx, E, f"OPENED on the forecast: rules {','.join(T['rules'])}, onset ~{T['onset_fc'][:13]}Z, "
                 f"pre planned {E['pre']['planned'][0]}..{E['pre']['planned'][1]}")
    queue_email(ctx, E, "heads_up")
    return E


def escalation(old, new):
    why = []
    oc = {a["code"] for a in old.get("alerts") or [] if a.get("role") == "strong"}
    nc = {a["code"] for a in new.get("alerts") or [] if a.get("role") == "strong"}
    if nc - oc:
        why.append("new strong alert(s): " + ", ".join(sorted(c for c in nc - oc if c)))
    op, np_ = (old.get("peak") or {}).get("hs"), (new.get("peak") or {}).get("hs")
    if op is not None and np_ is not None and np_ >= op + 1.0:
        why.append(f"forecast peak up from {op:.1f} to {np_:.1f} m")
    oo, no = iso_hour(old.get("onset_fc")), iso_hour(new.get("onset_fc"))
    if oo is not None and no is not None and no <= oo - 24:
        why.append(f"onset now ~{new['onset_fc'][:13]}Z, {(oo - no) / 24:.1f} days earlier")
    return why


def open_observed(ctx, state, segs, obs):
    cfg = ctx.cfg
    for s in sorted(segs, key=lambda x: x["start"]):
        if s["start"] < ctx.now_h - cfg["STORM_OBS_LOOKBACK_DAYS"] * 24 or not opens(s, cfg):
            continue
        if claimed(state, s, ctx.now_h):
            continue
        E = new_event(ctx, state, "obs", hour_date(s["start"]))
        E["status"] = "storm"
        adopt_segment(ctx, E, s)
        E["observed"]["confirmed_utc"] = iso_z(ctx.now)
        hist(ctx, E, f"OPENED from the observed waves (not forecast): {seg_text(s)}")
        queue_email(ctx, E, "missed_storm")
        update_storm(ctx, state, E, segs)


def ground_steps(ctx, state, E, gi, obs, segs):
    """Builds owed by an event, in order; each finished step is saved."""
    if E["status"] == "post" and not E["first_post_exact"]:
        set_first_post(ctx, E, gi)
    ob = E["observed"]
    if E["status"] == "post" and ob.get("cut_state") == "pending" and covered(gi, ob["cut_date"]):
        merge_or_cut(ctx, state, E, segs, gi)
        save_state(ctx, state)
    if E["status"] in ("storm", "post") and E["pre"].get("status") not in ("final", "inadequate", "disabled"):
        if try_pre(ctx, state, E, segs, gi, obs):
            save_state(ctx, state)
    if E["status"] == "watch" and E["pre_provisional"] is None and ctx.cron.get("DEM_ENABLE") == "1":
        try_pre(ctx, state, E, segs, gi, obs, provisional=True)
        save_state(ctx, state)
    if E["status"] == "post" and E["pre"].get("status") in ("final", "inadequate", "disabled"):
        if E["post3"] is None and try_post3(ctx, state, E, gi, obs, segs):
            save_state(ctx, state)
        if E["post3"] is not None and "post3" not in E["change"]:
            if E["post3"].get("status") == "none":
                E["change"]["post3"] = {"status": "no map", "why": E["post3"]["why"]}
            else:
                change_for(ctx, E, "post3", segs)
            save_state(ctx, state)
        if E["post3"] is not None and E["post7"] is None and try_post7(ctx, state, E, gi, obs):
            p7 = E["post7"]
            if p7.get("status") == "built":
                change_for(ctx, E, "post7", segs)
            else:
                E["change"]["post7"] = {"status": "no map", "why": p7.get("why") or p7.get("status")}
                hist(ctx, E, f"post7 {p7.get('status')}: {p7.get('why')}")
                queue_email(ctx, E, "report_post7")
            E["status"], E["closed_reason"] = "closed", "complete"
            hist(ctx, E, "closed: complete")
            save_state(ctx, state)
    if is_open(E):
        waits(ctx, state, E, gi)


def pending_steps(ctx, state, E, gi, obs, segs):
    keep = []
    for p in E["pending"]:
        if p["what"] == "pre_final" and p.get("window"):
            if covered(gi, p["window"][1]):
                if ctx.cron.get("DEM_ENABLE") == "1":
                    try_pre(ctx, state, E, segs, gi, obs, window=list(p["window"]))
                continue
            if ddiff(ctx.today, p["window"][1]) > ctx.cfg["STORM_MAX_WAIT_DAYS"]:
                E["notes"].append(f"{ctx.today}: pre {p['window']} never covered by the ground file: dropped")
                continue
        keep.append(p)
    E["pending"] = keep


def needed_from(ctx, state):
    """The earliest date any open event or pending build needs in the working ground rows."""
    cfg = ctx.cfg
    ds = []
    for E in state["events"]:
        if is_open(E):
            if E["pre"].get("status") in ("final", "inadequate", "disabled") and E["pre"].get("window"):
                ds.append(E["pre"]["window"][0])
            else:
                od = onset_date(E)
                if od:
                    ds.append(dadd(od, -cfg["STORM_PRE_MAX_DAYS"]))
        for p in E["pending"]:
            if p.get("window"):
                ds.append(dadd(p["window"][0], -cfg["STORM_PRE_MAX_DAYS"]))
    return min(ds) if ds else None


#!/usr/bin/env python3
"""
Weekly DEM And Waterline Frames, And Their Videos
==================================================
Every week (Monday to Sunday of UTC DATES, the pipeline's DEM day
convention: dem_from_contours.py filters on capture_time_utc[:10]) gets a
brand new DEM of its own and a plan view of its waterline points, drawn as
two FRAMES of two series that all share one frozen look (extent, colour
scale, pixel size), so the weeks can be played as a slideshow of the
intertidal beach changing:

  DEM frame        the week's DEM (dem_from_contours.py with exactly the
                   cron's window-DEM options, over Monday..Sunday)
  waterline frame  the week's waterline points from both cameras, coloured
                   by elevation on the same colour scale

A Sunday evening capture after 00:00 UTC (20:00 EDT, 19:00 EST) belongs to
the NEXT week, because the week is a week of UTC dates.

The pre- and post-storm DEMs of the storm watch (storm_watch.py) are drawn
in the same frame style by storm_watch.py itself, into
archive/storm_watch/events/<id>/frames/ (with a frames.json), and this
script's video builder puts them into the two videos at their dates
(change maps go into the DEM video after that date's DEM frames).

STATUS OF A WEEK
  current      the week holding today's UTC date: rebuilt by every daily run
               (shown in current/); on Monday a brand new one starts
  provisional  the week has ended, but the ground file does not yet cover
               its Sunday (GNSS-R runs about 2 days behind); its frames are
               stamped "PROVISIONAL: data through <date>"
  final        the ground file covers the Sunday (it has points of a later
               date), or WEEKLY_MAX_WAIT_DAYS (7) days have passed since the
               Sunday (then stamped "data through <date>"); never rebuilt
               (only redrawn by `respec`); copied to frames/<series>/
  empty        a final week with no DEM (no points, the day test rejected
               every camera-day, or fewer than 10 cells): its frames say why
               and keep the timeline continuous

FILES (BASE = the waterline folder, /mnt/I2Rgus_Data/waterline on CACO-05;
this script writes only under archive/weekly/ and logs/):
  archive/weekly/frame_spec.json            the frozen look of every frame
  archive/weekly/state.json                 last daily run, weeks, videos
  archive/weekly/weeks/<monday>/            wk_<monday>_<sunday>_{dem,spread,count,sigma}.asc,
                                            _info.json, points.csv, dem_frame.png,
                                            waterline_frame.png, week.json, dem_log.txt
  archive/weekly/current/                   the week in progress: dem_frame.png, waterline_frame.png
  archive/weekly/frames/dem/<sunday>_week.png, frames/waterlines/<sunday>_week.png
                                            the final weeks' frames (the videos' input)
  archive/weekly/video/dem_weekly.mp4, waterlines_weekly.mp4 (+ <series>_frames.txt, the
                                            ffmpeg concat list used)
  logs/weekly_frames.lock                   this script's own lock
The weekly grids never go into archive/dems (dem_change.py --series reads
every dem_*_7d file there).

Usage:
    python3 weekly_frames.py daily --base BASE [--have-lock] [--force]
        the cron's step (its work once per UTC date; --have-lock: the cron
        already holds its lock, so the ground file is read directly)
    python3 weekly_frames.py backfill --base BASE [--from YYYY-MM-DD] [--to YYYY-MM-DD] [--dry-run]
        every ended week from the ground file's first date (resumable)
    python3 weekly_frames.py respec --base BASE
        re-freeze the look from the settings and the ground file and redraw every frame
    python3 weekly_frames.py status --base BASE
        read only: the look, the weeks by status, the videos
Settings: BASE/weekly_frames_settings.conf (see weekly_frames_settings.conf.example).
A new station's BASE/waterline_station.conf is checked and its exported names
applied by every command, as the cron does, so a run by hand matches the cron's.
The video needs ffmpeg (sudo apt-get install -y ffmpeg); without it the
frames are still made. docs/WEEKLY_FRAMES.md has the details.

Exit status: 0 done (notes included), 1 a fault (a WARNING line), bad
arguments or settings. Python 3.8; numpy and matplotlib only.

BLIND: this script reads only the ground file, its own outputs, the
storm watch's frames.json/frames and the cron's settings files
(waterline_timex_cron.sh, a new station's waterline_station.conf); it never opens, globs or stats anything
under calibration/, surveys/ or ground_truth/, nor anything named for a
BLIND date (survey_products.BLIND_DATES, kept equal to BLIND_DATES here),
and skips every week that contains one.
"""

import os
import re
import sys
import csv
import json
import time
import array
import shutil
import struct
import hashlib
import argparse
import subprocess
import textwrap
from pathlib import Path
from datetime import datetime, date, timedelta, timezone

HERE = Path(__file__).resolve().parent
UTC = timezone.utc

# the blind dates: equal to survey_products.BLIND_DATES' keys (test_weekly_frames.py checks it)
BLIND_DATES = ("2024-10-23",)

SETTINGS_FILE = "weekly_frames_settings.conf"
OUT = "archive/weekly"
SERIES = ("dem", "waterlines")
VIDEO_NAMES = {"dem": "dem_weekly.mp4", "waterlines": "waterlines_weekly.mp4"}
FRAME_PNG = {"dem": "dem_frame.png", "waterlines": "waterline_frame.png"}
SPEC_VERSION = 1
DPI = 100                      # fixed: the frame size is set by WEEKLY_SIZE alone
MIN_CELLS = 10                 # dem_figure's MIN_CORE_CELLS: fewer cells = an empty week
EXTENT_PAD_M = 20.0
EXTENT_ROUND_M = 10.0
Z_ROUND_M = 0.25
CONTOUR_STEP_M = 0.5
LOCK_WAIT_S = 600.0            # waiting for the cron's lock (not when run by the cron: --have-lock)
CRON_HOLD_S = 900              # no ground read when a waterline cron run starts within 15 min
CATCHUP_WEEKS = 4              # missed weeks one daily run builds (it runs inside the cron's lock); backfill: all
DEM_TIMEOUT_S = 3600
FFMPEG_TIMEOUT_S = 1800
DEFAULT_STATION = "CACO-05 Marconi Beach"
NOT_COMPARED = "#c9c8c3"       # dem_change.py's colours
NOT_SIGNIFICANT = "#ebeae6"
NO_DATA = "#e4e3df"
BANNER_RED = "#b3261e"
LOOK_KEYS = ("version", "extent", "z_limits", "change_limit", "size", "dpi", "cell", "epsg", "frame_label",
             "station_name")

# KEY: (default, kind, meaning)
SETTINGS = {
    "WEEKLY_ENABLE": ("1", "bool", "1: the daily step makes the weekly frames; 0: it does nothing"),
    "WEEKLY_FPS": ("2", "fps", "frames per second of the videos (each frame shown 1/fps s)"),
    "WEEKLY_MAX_WAIT_DAYS": ("7", "days", "an ended week is final at the latest this many days after its Sunday"),
    "WEEKLY_EXTENT": ("", "extent", "'' = auto (all usable ground points, padded 20 m, rounded to 10 m), "
                                    "else 'e_min,e_max,n_min,n_max' (m)"),
    "WEEKLY_Z_LIMITS": ("", "zlim", "'' = auto (1st-99th percentile, rounded to 0.25 m), else 'z_min,z_max' "
                                    "(m NAVD88)"),
    "WEEKLY_CHANGE_LIMIT": ("0.5", "pos", "the change frames' symmetric colour limit (m)"),
    "WEEKLY_SIZE": ("1920x1080", "size", "frame size in pixels, WIDTHxHEIGHT (even numbers)"),
}
SPEC_SETTINGS = ("WEEKLY_FPS", "WEEKLY_EXTENT", "WEEKLY_Z_LIMITS", "WEEKLY_CHANGE_LIMIT", "WEEKLY_SIZE")


# --------------------------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------------------------

def utcnow():
    return datetime.now(UTC)


def local_now():
    """The NUC's clock (naive local time), for the station's job table."""
    return datetime.now()


def iso_z(t):
    return t.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def d_iso(d):
    return d.isoformat()


def d_of(s):
    return date.fromisoformat(str(s)[:10])


def monday_of(d):
    d = d_of(d) if isinstance(d, str) else d
    return d - timedelta(days=d.weekday())


def sunday_of(monday):
    return d_iso(d_of(monday) + timedelta(days=6))


def fmt_num(x):
    """2.0 -> '2', 0.5 -> '0.5' (ffmpeg arguments, titles)."""
    return ("%.6f" % float(x)).rstrip("0").rstrip(".")


MONTHS = ("January", "February", "March", "April", "May", "June", "July", "August", "September", "October",
          "November", "December")


def blind_forms(d):
    """The name forms of a blind date that count (docs/VALIDATION.md): 2024-10-23, 20241023, 2024-Oct-23,
    10-23-2024, Oct_23_2024, October 23 2024, 23Oct2024, and frame names ...Oct.23_hh_mm_ss.GMT.2024..."""
    y, m, dd = d.split("-")
    mon, month = MONTHS[int(m) - 1][:3], MONTHS[int(m) - 1]
    return [d, f"{y}{m}{dd}", f"{y}-{mon}-{dd}", f"{m}-{dd}-{y}", f"{mon}_{dd}_{y}", f"{month} {dd} {y}",
            f"{dd}{mon}{y}"], re.compile(rf"{mon}\.{dd}_\d\d_\d\d_\d\d\.GMT\.{y}", re.I)


def blind_named(name):
    """True when a file or folder name carries a blind date in any of its forms."""
    low = str(name).lower()
    for d in BLIND_DATES:
        forms, rx = blind_forms(d)
        if any(f.lower() in low for f in forms) or rx.search(str(name)):
            return True
    return False


def blind_in_week(monday):
    """The blind date inside Monday..Sunday, or None."""
    m = d_iso(d_of(monday))
    s = sunday_of(m)
    for d in BLIND_DATES:
        if m <= d <= s:
            return d
    return None


def week_title(monday, station):
    """'<station>: week of Mon 17 Nov to Sun 23 Nov 2026 (UTC dates)'."""
    a = d_of(monday)
    b = a + timedelta(days=6)
    first = f"Mon {a.day} {a.strftime('%b')}" + (f" {a.year}" if a.year != b.year else "")
    return f"{station}: week of {first} to Sun {b.day} {b.strftime('%b')} {b.year} (UTC dates)"


def storm_frame_name(date, event_id, part):
    """The name a storm frame sorts by among the weekly frames: '<date>_storm_<event id>_<part>.png'
    (its first 10 characters are its date, like '<sunday>_week.png')."""
    return f"{str(date)[:10]}_storm_{event_id}_{part}.png"


def read_json(p, default=None):
    try:
        return json.loads(Path(p).read_text())
    except (OSError, ValueError):
        return default


def write_json(p, obj):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, sort_keys=False) + "\n")
    os.replace(str(tmp), str(p))


def copy_atomic(src, dst):
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.name + ".tmp")
    shutil.copyfile(str(src), str(tmp))
    os.replace(str(tmp), str(dst))


def sw():
    """storm_watch.py, imported when first needed (it imports this module lazily too)."""
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    import storm_watch
    return storm_watch


def png_size(path):
    """(width, height) from a PNG's IHDR, or None."""
    try:
        with open(path, "rb") as f:
            head = f.read(24)
    except OSError:
        return None
    if len(head) < 24 or head[:8] != b"\x89PNG\r\n\x1a\n" or head[12:16] != b"IHDR":
        return None
    return struct.unpack(">II", head[16:24])


def png_text(path):
    """The tEXt chunks of a PNG as {key: value}."""
    out = {}
    try:
        b = Path(path).read_bytes()
    except OSError:
        return out
    i = 8
    while i + 8 <= len(b):
        n, typ = struct.unpack(">I4s", b[i:i + 8])
        if typ == b"tEXt":
            k, _, v = b[i + 8:i + 8 + n].partition(b"\x00")
            out[k.decode("latin-1")] = v.decode("latin-1")
        if typ == b"IEND":
            break
        i += 12 + n
    return out


# --------------------------------------------------------------------------------------------
# settings
# --------------------------------------------------------------------------------------------

def parse_setting(v, kind):
    v = str(v).strip()
    if kind == "bool":
        if v not in ("0", "1"):
            raise ValueError("0 or 1")
        return v == "1"
    if kind == "fps":
        x = float(v)
        if not 0.1 <= x <= 60:
            raise ValueError("0.1 to 60")
        return x
    if kind == "days":
        x = int(v)
        if not 0 <= x <= 60:
            raise ValueError("whole days, 0 to 60")
        return x
    if kind == "pos":
        x = float(v)
        if not 0 < x <= 20:
            raise ValueError("above 0, at most 20")
        return x
    if kind == "extent":
        if v == "":
            return None
        p = [float(x) for x in v.split(",")]
        if len(p) != 4 or not (p[0] < p[1] and p[2] < p[3]):
            raise ValueError("e_min,e_max,n_min,n_max with e_min < e_max and n_min < n_max")
        return p
    if kind == "zlim":
        if v == "":
            return None
        p = [float(x) for x in v.split(",")]
        if len(p) != 2 or not p[0] < p[1]:
            raise ValueError("z_min,z_max with z_min < z_max")
        return p
    if kind == "size":
        m = re.match(r"^(\d+)x(\d+)$", v)
        if not m:
            raise ValueError("WIDTHxHEIGHT")
        w, h = int(m.group(1)), int(m.group(2))
        if not (320 <= w <= 7680 and 240 <= h <= 4320) or w % 2 or h % 2:
            raise ValueError("even numbers, 320x240 to 7680x4320 (H.264 needs even sizes)")
        return [w, h]
    raise ValueError(kind)


def load_settings(base):
    """BASE/weekly_frames_settings.conf (KEY=VALUE, '#' comments) over the defaults. -> (cfg, raw, notes): cfg the
    parsed values, raw the text values given in the file (for the frozen-spec check). An unknown key or a value
    that does not parse is a note (the default, or the earlier valid line, is kept)."""
    raw, notes = {}, []
    p = Path(base) / SETTINGS_FILE
    if p.exists():
        try:
            found = sw().var_lines(p.read_text(errors="replace"))
        except OSError as exc:
            found = {}
            notes.append(f"{SETTINGS_FILE} not read ({exc.strerror or exc}): the defaults are used")
        for k, vals in found.items():
            if k not in SETTINGS:
                notes.append(f"setting {k} in {SETTINGS_FILE} is not known: ignored")
                continue
            for v in vals:
                try:
                    parse_setting(v, SETTINGS[k][1])
                except ValueError as exc:
                    kept = f"{raw[k]!r}" if k in raw else f"the default {SETTINGS[k][0]!r}"
                    notes.append(f"setting {k}={v!r} in {SETTINGS_FILE} is not valid ({exc}): {kept} is used")
                    continue
                raw[k] = v
    cfg = {k: parse_setting(raw.get(k, d), kind) for k, (d, kind, _) in SETTINGS.items()}
    return cfg, raw, notes


# --------------------------------------------------------------------------------------------
# the ground file: one pass (under the cron's lock unless the cron itself runs us)
# --------------------------------------------------------------------------------------------

def frame_passes(hs_text, max_hs):
    if max_hs is None or hs_text == "":
        return True
    try:
        return float(hs_text) <= max_hs
    except ValueError:
        return True


def scan_ground(path, d_min, d_max, out_csv=None, max_hs=None, stats=False):
    """One pass of the ground file. Copies the rows dated [d_min, d_max] (UTC dates; never a blind date's) to
    out_csv (tmp + replace); per date in that range: georectified rows, frames and frames passing --max-hs (only
    when the file has offshore_hs_m), and a digest of the date's rows (a re-processed file with the same counts
    but other values still differs); the first and latest georectified dates of the WHOLE file; with stats: the
    extent and the elevations of every usable point (georectified, passing --max-hs, not a blind date) for the
    frozen look."""
    info = {"rows": 0, "first": None, "latest": None, "has_hs": False, "dates": {}, "copied": 0,
            "range": [d_min, d_max], "stats": None, "digest": {}}
    digests = {}
    st = {"e": [float("inf"), float("-inf")], "n": [float("inf"), float("-inf")], "z": array.array("d"),
          "first": None, "last": None} if stats else None
    frames = {}
    tmp = Path(str(out_csv) + ".tmp") if out_csv else None
    if tmp:
        tmp.parent.mkdir(parents=True, exist_ok=True)
    g = open(str(tmp), "w", newline="") if tmp else None
    try:
        with open(str(path), newline="", errors="replace") as f:
            r = csv.reader(f)
            head = next(r, None)
            if not head or "capture_time_utc" not in head or "easting_utm19" not in head:
                raise ValueError("no capture_time_utc / easting_utm19 column")
            w = csv.writer(g) if g else None
            if w:
                w.writerow(head)
            ix = {k: (head.index(k) if k in head else None) for k in
                  ("capture_time_utc", "easting_utm19", "northing_utm19", "camera", "source_file", "offshore_hs_m",
                   "beach_elevation_navd88", "tide_elevation_navd88")}
            kt, ke, kn = ix["capture_time_utc"], ix["easting_utm19"], ix["northing_utm19"]
            kh, kb, kz, kc, ks = (ix["offshore_hs_m"], ix["beach_elevation_navd88"], ix["tide_elevation_navd88"],
                                  ix["camera"], ix["source_file"])
            info["has_hs"] = kh is not None
            mh = max_hs if kh is not None else None
            first, latest = "", ""
            need = max(k for k in (kt, ke, kn) if k is not None)
            for row in r:
                if len(row) <= need:
                    continue
                info["rows"] += 1
                d = row[kt][:10]
                if d in BLIND_DATES:
                    continue
                ground = row[ke] != "" and (kn is None or row[kn] != "")
                if ground:
                    if d > latest:
                        latest = d
                    if not first or d < first:
                        first = d
                hs = row[kh] if kh is not None and kh < len(row) else ""
                ok = frame_passes(hs, mh)
                if d_min <= d <= d_max:
                    if w:
                        w.writerow(row)
                    info["copied"] += 1
                    h = digests.get(d)
                    if h is None:
                        h = digests[d] = hashlib.sha1()
                    h.update(("\x1f".join(row) + "\n").encode("utf-8", "replace"))
                    if ground:
                        fr = row[ks] if ks is not None and ks < len(row) else f"{row[kt]}|{row[kc] if kc is not None else ''}"
                        rec = info["dates"].setdefault(d, [0, 0, 0])
                        rec[0] += 1
                        key = (d, fr)
                        if key not in frames:
                            frames[key] = ok
                            rec[1] += 1
                            rec[2] += 1 if ok else 0
                if st is not None and ground and ok:
                    try:
                        e, n = float(row[ke]), float(row[kn])
                        zt = row[kb] if kb is not None and kb < len(row) and row[kb] != "" else \
                            (row[kz] if kz is not None and kz < len(row) else "")
                        z = float(zt)
                    except (ValueError, TypeError):
                        continue
                    if not all(map(lambda x: x == x and abs(x) != float("inf"), (e, n, z))):
                        continue
                    st["e"] = [min(st["e"][0], e), max(st["e"][1], e)]
                    st["n"] = [min(st["n"][0], n), max(st["n"][1], n)]
                    st["z"].append(z)
                    if st["first"] is None or d < st["first"]:
                        st["first"] = d
                    if st["last"] is None or d > st["last"]:
                        st["last"] = d
            info["first"], info["latest"] = first or None, latest or None
            info["digest"] = {d: h.hexdigest()[:16] for d, h in digests.items()}
    finally:
        if g:
            g.close()
    if tmp:
        os.replace(str(tmp), str(out_csv))
    info["stats"] = st
    return info


def read_points(ground_csv, start, end, max_hs):
    """The waterline points of [start, end] (UTC dates, inclusive) of a ground file (or of a points.csv this script
    wrote): dicts with easting_utm19, northing_utm19, elevation (beach_elevation_navd88, else
    tide_elevation_navd88), capture_time_utc, camera, source_file, in time order. Rows without coordinates are
    skipped; frames whose offshore_hs_m is above max_hs are dropped (blank heights kept), as the DEM drops them."""
    mh = None if max_hs in (None, "") else float(max_hs)
    out = []
    with open(str(ground_csv), newline="", errors="replace") as f:
        r = csv.DictReader(f)
        fields = r.fieldnames or []
        own = "elevation" in fields
        for row in r:
            t = row.get("capture_time_utc") or ""
            d = t[:10]
            if not d or d < start or d > end or d in BLIND_DATES:
                continue
            e, n = row.get("easting_utm19") or "", row.get("northing_utm19") or ""
            if e == "" or n == "":
                continue
            if mh is not None and not frame_passes(row.get("offshore_hs_m") or "", mh):
                continue
            zt = row.get("elevation") if own else (row.get("beach_elevation_navd88") or row.get("tide_elevation_navd88"))
            try:
                rec = {"easting_utm19": float(e), "northing_utm19": float(n), "elevation": float(zt),
                       "capture_time_utc": t, "camera": row.get("camera") or "",
                       "source_file": row.get("source_file") or ""}
            except (TypeError, ValueError):
                continue
            out.append(rec)
    out.sort(key=lambda p: p["capture_time_utc"])
    return out


POINT_FIELDS = ("capture_time_utc", "camera", "source_file", "easting_utm19", "northing_utm19", "elevation")


def write_points(points, path):
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    with open(str(tmp), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(POINT_FIELDS)
        for p in points:
            w.writerow([p.get("capture_time_utc", ""), p.get("camera", ""), p.get("source_file", ""),
                        "%.2f" % p["easting_utm19"], "%.2f" % p["northing_utm19"], "%.4f" % p["elevation"]])
    os.replace(str(tmp), str(path))


# --------------------------------------------------------------------------------------------
# the frozen look (frame_spec.json)
# --------------------------------------------------------------------------------------------

def spec_path(base):
    return Path(base) / OUT / "frame_spec.json"


def spec_hash(spec):
    look = {k: spec.get(k) for k in LOOK_KEYS}
    return hashlib.sha256(json.dumps(look, sort_keys=True).encode()).hexdigest()[:16]


def load_spec(base):
    """The frozen frame spec of BASE (archive/weekly/frame_spec.json), or None when there is none yet."""
    s = read_json(spec_path(base))
    if not isinstance(s, dict) or not s.get("extent") or not s.get("z_limits") or not s.get("size"):
        return None
    return s


STATION_FILE = "waterline_station.conf"


def apply_station_file(base):
    """As the cron does before its Python steps: a new station's BASE/waterline_station.conf (coastcam-station-
    setup's; CACO-05 has none), checked by waterline_station_check.check_overlay, and the names its export line
    exports put into the environment (the file's values win, as when the cron sources it). So a run by hand
    (backfill, respec, status, daily without the cron) freezes the station's own name and frame and builds its
    DEMs with the station's time zone and waves, like the cron's run. -> (names set, error text or '')."""
    p = Path(base) / STATION_FILE
    if not p.is_file():
        return [], ""
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    import waterline_station_check
    try:
        res = waterline_station_check.check_overlay(p.read_text(errors="replace"))
    except OSError as exc:
        return [], f"{STATION_FILE} not read ({exc.strerror or exc})"
    if not res["ok"]:
        e = res["errors"][0]
        return [], (f"{STATION_FILE} refused by waterline_station_check.py check-overlay (as the cron refuses it): "
                    + (f"line {e['line']}: " if e.get("line") else "") + str(e.get("message")))
    names = [k for k in res["exported"] if k in res["keys"]]
    for k in names:
        os.environ[k] = res["keys"][k]
    return names, ""


def station_label():
    import station_env
    return station_env.station_name(DEFAULT_STATION)


def frame_labels():
    import station_env
    return station_env.geotiff_epsg(), station_env.frame_label()


def make_spec(cfg, stats, cron, ground_name, rows, now=None):
    """A new frame spec from the settings and (for what they leave on auto) the usable ground points.
    -> (spec, why): spec None when the ground file has no usable point and the settings do not give the look."""
    import numpy as np
    z = None
    if stats and len(stats["z"]):
        z = np.frombuffer(stats["z"], dtype=float)
    if cfg["WEEKLY_EXTENT"]:
        ext, ext_src = [float(x) for x in cfg["WEEKLY_EXTENT"]], "WEEKLY_EXTENT"
    elif z is not None:
        r = EXTENT_ROUND_M
        ext = [np.floor((stats["e"][0] - EXTENT_PAD_M) / r) * r, np.ceil((stats["e"][1] + EXTENT_PAD_M) / r) * r,
               np.floor((stats["n"][0] - EXTENT_PAD_M) / r) * r, np.ceil((stats["n"][1] + EXTENT_PAD_M) / r) * r]
        ext, ext_src = [float(x) for x in ext], "auto"
    else:
        return None, "no usable ground point yet (and WEEKLY_EXTENT is auto)"
    if cfg["WEEKLY_Z_LIMITS"]:
        zl, z_src = [float(x) for x in cfg["WEEKLY_Z_LIMITS"]], "WEEKLY_Z_LIMITS"
    elif z is not None:
        lo, hi = np.percentile(z, [1, 99])
        zl = [float(np.floor(lo / Z_ROUND_M) * Z_ROUND_M), float(np.ceil(hi / Z_ROUND_M) * Z_ROUND_M)]
        if zl[1] - zl[0] < 2 * Z_ROUND_M:
            zl = [zl[0] - Z_ROUND_M, zl[1] + Z_ROUND_M]
        z_src = "auto"
    else:
        return None, "no usable ground point yet (and WEEKLY_Z_LIMITS is auto)"
    epsg, label = frame_labels()
    spec = {"version": SPEC_VERSION, "extent": ext, "extent_source": ext_src, "z_limits": zl, "z_source": z_src,
            "change_limit": float(cfg["WEEKLY_CHANGE_LIMIT"]), "size": list(cfg["WEEKLY_SIZE"]), "dpi": DPI,
            "fps": float(cfg["WEEKLY_FPS"]), "cell": float(cron.get("DEM_CELL") or 2.0), "epsg": int(epsg),
            "frame_label": label, "station_name": station_label(),
            "created_utc": iso_z(now or utcnow()),
            "made_from": {"ground_file": ground_name, "rows": int(rows),
                          "usable_points": int(len(z)) if z is not None else 0,
                          "first_date": stats["first"] if stats else None,
                          "last_date": stats["last"] if stats else None}}
    spec["hash"] = spec_hash(spec)
    return spec, ""


def spec_differences(spec, cfg, raw, cron=None):
    """What the settings (or the station's name, frame, the cron's cell) now say that differs from the frozen
    spec: one line each, naming `respec`."""
    out = []
    tail = "not applied to the existing frames (the look is frozen); `weekly_frames.py respec --base BASE` " \
           "re-freezes it and redraws every frame"

    def differ(a, b):
        if isinstance(a, (list, tuple)):
            return len(a) != len(b) or any(abs(float(x) - float(y)) > 1e-9 for x, y in zip(a, b))
        return abs(float(a) - float(b)) > 1e-9
    checks = (("WEEKLY_FPS", "fps"), ("WEEKLY_EXTENT", "extent"), ("WEEKLY_Z_LIMITS", "z_limits"),
              ("WEEKLY_CHANGE_LIMIT", "change_limit"), ("WEEKLY_SIZE", "size"))
    for k, sk in checks:
        if k not in raw or cfg[k] is None:
            continue
        if differ(cfg[k], spec.get(sk)):
            out.append(f"{k}={raw[k]} differs from the frozen frame spec's {sk} {spec.get(sk)}: {tail}")
    name = station_label()
    if name != spec.get("station_name"):
        out.append(f"STATION_NAME {name!r} differs from the frozen frame spec's {spec.get('station_name')!r}: {tail}")
    epsg, _ = frame_labels()
    if int(epsg) != int(spec.get("epsg") or 0):
        out.append(f"GEOTIFF_EPSG {epsg} differs from the frozen frame spec's {spec.get('epsg')}: {tail}")
    if cron and cron.get("DEM_CELL") and differ(float(cron["DEM_CELL"]), spec.get("cell") or 0):
        out.append(f"the cron's DEM_CELL {cron['DEM_CELL']} differs from the frozen frame spec's cell "
                   f"{spec.get('cell')}: {tail}")
    return out


def ensure_spec(base, ground_csv, settings=None, log=None, stats=None, cron=None, raw=None):
    """The frozen frame spec of BASE; made (and written to archive/weekly/frame_spec.json) from the settings and the
    usable points of `ground_csv` when there is none. A setting that differs from the frozen value is NOT applied:
    a NOTE names `respec`. -> dict, or None when the ground file has no usable point and the settings do not give
    the extent and the elevation limits."""
    log = log or (lambda s: print(s, flush=True))
    if settings is None:
        cfg, raw, notes = load_settings(base)
    else:
        cfg = settings
        if raw is None:
            raw = {k: str(v) for k, v in settings.items() if k in SETTINGS and v is not None}
    if cron is None:
        cron, _ = sw().cron_settings(base)
    spec = load_spec(base)
    if spec is not None:
        for line in spec_differences(spec, cfg, raw, cron):
            log("NOTE: weekly frames: " + line)
        return spec
    rows = 0
    if stats is None:
        mh = float(cron["DEM_MAX_HS"]) if cron.get("DEM_MAX_HS") else None
        info = scan_ground(ground_csv, "9999-12-31", "9999-12-31", None, mh, stats=True)
        stats, rows = info["stats"], info["rows"]
    else:
        rows = stats.get("rows", 0)
    spec, why = make_spec(cfg, stats, cron, Path(ground_csv).name, rows)
    if spec is None:
        log(f"NOTE: weekly frames: no frame spec yet: {why}")
        return None
    write_json(spec_path(base), spec)
    log(f"weekly frames: frame spec frozen ({spec['hash']}): extent {spec['extent']} ({spec['extent_source']}), "
        f"elevations {spec['z_limits']} m ({spec['z_source']}), {spec['size'][0]}x{spec['size'][1]} px, "
        f"{fmt_num(spec['fps'])} fps")
    return spec


# --------------------------------------------------------------------------------------------
# drawing (matplotlib Agg; every frame the same size, extent and colour scale; atomic writes)
# --------------------------------------------------------------------------------------------

def read_asc(path):
    """An ESRI ASCII grid -> (grid with row 0 = south, NaN for nodata, header dict)."""
    import numpy as np
    hdr = {}
    with open(str(path)) as f:
        for _ in range(6):
            k, v = f.readline().split()
            hdr[k.lower()] = float(v)
    g = np.loadtxt(str(path), skiprows=6, ndmin=2)
    g = np.where(g == hdr.get("nodata_value", -9999.0), np.nan, g)
    return np.flipud(g), hdr


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _layout(spec):
    w, h = spec["size"]
    return {"s": h / 1080.0, "w": w, "h": h, "map": [0.045, 0.06, 0.78, 0.745], "cbar": [0.875, 0.33, 0.016, 0.47],
            "stamp": (0.85, 0.27)}


def _figure(spec, title, subtitle, banner):
    plt = _mpl()
    L = _layout(spec)
    dpi = float(spec.get("dpi") or DPI)
    # + 0.25 px: matplotlib 3.3 (the NUC's) truncates width/dpi*dpi, so 854 px would come out 853 (8.54 * 100 =
    # 853.99...); later versions round, and 854.25 rounds to 854 too
    fig = plt.figure(figsize=((L["w"] + 0.25) / dpi, (L["h"] + 0.25) / dpi), dpi=dpi)
    fig.patch.set_facecolor("white")
    s = L["s"]
    if banner:
        from matplotlib.patches import Rectangle
        fig.add_artist(Rectangle((0, 0.94), 1, 0.06, transform=fig.transFigure, facecolor=BANNER_RED,
                                 edgecolor="none", zorder=0))
        fig.text(0.5, 0.97, banner, ha="center", va="center", color="white", fontsize=22 * s, fontweight="bold")
    fig.text(0.5, 0.905, title, ha="center", va="center", fontsize=24 * s, fontweight="bold")
    if subtitle:
        fig.text(0.5, 0.872, "\n".join(textwrap.wrap(subtitle, 175)), ha="center", va="top", fontsize=13 * s,
                 color="0.25")
    ax = fig.add_axes(L["map"])
    e0, e1, n0, n1 = spec["extent"]
    ax.set_xlim(e0, e1)
    ax.set_ylim(n0, n1)
    ax.set_aspect("equal", adjustable="box", anchor="C")
    ax.set_facecolor(NO_DATA)
    ax.ticklabel_format(useOffset=False, style="plain")
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(8))
    ax.tick_params(labelsize=11 * s)
    short = spec.get("frame_label", "").split(",")[0].split(" / ")[-1] or "UTM"
    ax.set_xlabel(f"easting (m, {short})", fontsize=12 * s)
    ax.set_ylabel(f"northing (m, {short})", fontsize=12 * s)
    _scale_bar(ax, spec, s)
    _north_arrow(ax, s)
    return plt, fig, ax, L


def _scale_bar(ax, spec, s):
    e0, e1, n0, n1 = spec["extent"]
    want = (e1 - e0) / 5.0
    nice = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 5000]
    L = max([x for x in nice if x <= want] or [nice[0]])
    x0, y0 = e0 + 0.04 * (e1 - e0), n0 + 0.04 * (n1 - n0)
    ax.plot([x0, x0 + L], [y0, y0], color="black", lw=4 * s, solid_capstyle="butt", zorder=6)
    ax.plot([x0, x0 + L / 2.0], [y0, y0], color="white", lw=2 * s, solid_capstyle="butt", zorder=7)
    ax.text(x0 + L / 2.0, y0 + 0.015 * (n1 - n0), f"{fmt_num(L)} m", ha="center", va="bottom", fontsize=11 * s,
            zorder=7)


def _north_arrow(ax, s):
    ax.annotate("N", xy=(0.95, 0.95), xytext=(0.95, 0.85), xycoords="axes fraction", ha="center", va="center",
                fontsize=16 * s, fontweight="bold",
                arrowprops=dict(facecolor="black", width=4 * s, headwidth=14 * s, headlength=14 * s), zorder=8)


def _stamp(fig, L, stamp):
    if not stamp:
        return
    s = L["s"]
    lines = []
    for part in str(stamp).split("\n"):
        lines += textwrap.wrap(part, 22) or [""]
    color = "#a35200" if "PROVISIONAL" in stamp or "IN PROGRESS" in stamp else (
        "0.35" if "EMPTY" in stamp else "0.15")
    fig.text(L["stamp"][0], L["stamp"][1], "\n".join(lines), ha="left", va="top", fontsize=13 * s,
             fontweight="bold", color=color,
             bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor=color, linewidth=1.5 * s))


def _colorbar(fig, L, spec, cmap, lo, hi, label, step=None):
    import numpy as np
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    cax = fig.add_axes(L["cbar"])
    sm = ScalarMappable(norm=Normalize(vmin=lo, vmax=hi), cmap=cmap)
    sm.set_array(np.array([lo, hi]))
    span = hi - lo
    step = step or (0.25 if span <= 2.0 else 0.5 if span <= 5.0 else 1.0)
    ticks = np.arange(np.ceil(lo / step - 1e-9) * step, hi + 1e-9, step)
    cb = fig.colorbar(sm, cax=cax, ticks=ticks)
    cb.ax.tick_params(labelsize=11 * L["s"])
    cb.set_label(label, fontsize=12 * L["s"])
    return cb


def _center_note(ax, text, s):
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center", fontsize=26 * s, color="0.45",
            fontweight="bold", zorder=9)


def _save(plt, fig, out_png, spec, kind):
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_png.with_name(out_png.name + ".tmp")
    try:
        fig.savefig(str(tmp), format="png", dpi=float(spec.get("dpi") or DPI), facecolor="white",
                    metadata={"Software": "weekly_frames.py", "weekly_frames_spec": spec.get("hash") or spec_hash(spec),
                              "weekly_frames_kind": kind})
    finally:
        plt.close(fig)
    os.replace(str(tmp), str(out_png))


def draw_dem_frame(stem, out_png, spec, title, subtitle="", stamp="", banner=""):
    """The DEM frame of `stem` (<stem>_dem.asc; missing: an empty map) on the frozen extent and elevation scale:
    terrain colours, 0.5 m contours, fixed colour bar, scale bar, north arrow, title, subtitle, stamp, banner."""
    import numpy as np
    plt, fig, ax, L = _figure(spec, title, subtitle, banner)
    try:
        lo, hi = spec["z_limits"]
        p = Path(str(stem) + "_dem.asc") if stem else None
        if p is not None and p.exists():
            g, h = read_asc(p)
            c = h["cellsize"]
            x0, y0 = h["xllcorner"], h["yllcorner"]
            ext = [x0, x0 + g.shape[1] * c, y0, y0 + g.shape[0] * c]
            ax.imshow(g, origin="lower", extent=ext, cmap="terrain", vmin=lo, vmax=hi, interpolation="nearest",
                      zorder=2)
            levels = np.arange(np.ceil(lo / CONTOUR_STEP_M) * CONTOUR_STEP_M, hi + 1e-9, CONTOUR_STEP_M)
            if np.isfinite(g).sum() >= 4 and g.shape[0] >= 2 and g.shape[1] >= 2 and len(levels):
                xs = x0 + (np.arange(g.shape[1]) + 0.5) * c
                ys = y0 + (np.arange(g.shape[0]) + 0.5) * c
                try:
                    cs = ax.contour(xs, ys, np.ma.masked_invalid(g), levels=levels, colors="k",
                                    linewidths=0.6 * L["s"], zorder=3)
                    ax.clabel(cs, fmt="%.1f", fontsize=8 * L["s"])
                except (ValueError, TypeError):
                    pass
            ax.set_xlim(spec["extent"][0], spec["extent"][1])
            ax.set_ylim(spec["extent"][2], spec["extent"][3])
            if not np.isfinite(g).any():
                _center_note(ax, "no DEM cells", L["s"])
        else:
            _center_note(ax, "no DEM this week", L["s"])
        _colorbar(fig, L, spec, "terrain", lo, hi, "elevation (m NAVD88)")
        _stamp(fig, L, stamp)
    except Exception:
        plt.close(fig)
        raise
    _save(plt, fig, out_png, spec, "dem")


def draw_waterline_frame(points, out_png, spec, title, subtitle="", stamp="", banner=""):
    """The waterline frame: a plan view of `points` (dicts with easting_utm19, northing_utm19, elevation,
    capture_time_utc) in time order, coloured by elevation with the DEM frame's scale and colour map.
    -> the number of points drawn."""
    import numpy as np
    plt, fig, ax, L = _figure(spec, title, subtitle, banner)
    try:
        lo, hi = spec["z_limits"]
        pts = sorted(points or [], key=lambda p: str(p.get("capture_time_utc", "")))
        E = np.array([float(p["easting_utm19"]) for p in pts], dtype=float)
        N = np.array([float(p["northing_utm19"]) for p in pts], dtype=float)
        Z = np.array([float(p["elevation"]) for p in pts], dtype=float)
        ok = np.isfinite(E) & np.isfinite(N) & np.isfinite(Z)
        n = int(ok.sum())
        if n:
            ax.scatter(E[ok], N[ok], c=Z[ok], cmap="terrain", vmin=lo, vmax=hi, s=4.0 * L["s"] ** 2,
                       linewidths=0, zorder=2)
            ax.set_xlim(spec["extent"][0], spec["extent"][1])
            ax.set_ylim(spec["extent"][2], spec["extent"][3])
        else:
            _center_note(ax, "no waterline points this week", L["s"])
        _colorbar(fig, L, spec, "terrain", lo, hi, "elevation (m NAVD88)")
        _stamp(fig, L, stamp)
    except Exception:
        plt.close(fig)
        raise
    _save(plt, fig, out_png, spec, "waterlines")
    return n


def draw_change_frame(diff_stem, out_png, spec, title, subtitle="", banner=""):
    """The change frame of dem_change.py's <diff_stem>_diff.asc (and _sig.asc) on the frozen extent: RdBu with
    the frozen symmetric change limit; cells within the level of detection (in _diff, not in _sig) light grey, cells
    not compared the panel's grey, as dem_change.py draws them."""
    import numpy as np
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    plt, fig, ax, L = _figure(spec, title, subtitle, banner)
    try:
        lim = float(spec.get("change_limit") or 0.5)
        ax.set_facecolor(NOT_COMPARED)
        p = Path(str(diff_stem) + "_diff.asc")
        if p.exists():
            d, h = read_asc(p)
            sp = Path(str(diff_stem) + "_sig.asc")
            sig = read_asc(sp)[0] if sp.exists() else d
            if sig.shape != d.shape:
                sig = d
            c = h["cellsize"]
            x0, y0 = h["xllcorner"], h["yllcorner"]
            ext = [x0, x0 + d.shape[1] * c, y0, y0 + d.shape[0] * c]
            within = np.isfinite(d) & ~np.isfinite(sig)
            ax.imshow(np.where(within, 1.0, np.nan), origin="lower", extent=ext,
                      cmap=ListedColormap([NOT_SIGNIFICANT]), vmin=0, vmax=1, interpolation="nearest", zorder=2)
            ax.imshow(sig, origin="lower", extent=ext, cmap="RdBu", vmin=-lim, vmax=lim, interpolation="nearest",
                      zorder=3)
            ax.set_xlim(spec["extent"][0], spec["extent"][1])
            ax.set_ylim(spec["extent"][2], spec["extent"][3])
        else:
            _center_note(ax, "no change grid", L["s"])
        ax.legend(handles=[Patch(facecolor=NOT_COMPARED, edgecolor="0.4", label="not compared"),
                           Patch(facecolor=NOT_SIGNIFICANT, edgecolor="0.4",
                                 label="within the level of detection")],
                  loc="lower right", fontsize=11 * L["s"], framealpha=0.9)
        step = 0.1 if lim <= 0.5 else 0.25 if lim <= 1.5 else 0.5
        _colorbar(fig, L, spec, "RdBu", -lim, lim, "elevation change (m); blue = accretion", step=step)
    except Exception:
        plt.close(fig)
        raise
    _save(plt, fig, out_png, spec, "change")


# --------------------------------------------------------------------------------------------
# the run context, logging, locks
# --------------------------------------------------------------------------------------------

def subprocess_runner(argv, log_path=None):
    """Runs argv (output captured); appends the command and its output to log_path. -> (exit code, output)."""
    env = dict(os.environ)
    env.setdefault("MPLBACKEND", "Agg")
    try:
        p = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env,
                           timeout=DEM_TIMEOUT_S, universal_newlines=True, errors="replace")
        rc, out = p.returncode, p.stdout or ""
    except subprocess.TimeoutExpired as exc:
        rc, out = 124, f"timed out after {exc.timeout} s\n"
    except OSError as exc:
        rc, out = 127, f"could not run: {exc}\n"
    if log_path:
        with open(str(log_path), "a") as f:
            f.write(f"$ {' '.join(map(str, argv))}\n{out}\n[exit {rc}]\n")
    return rc, out


class Ctx:
    def __init__(self, base, now=None, runner=None, have_lock=False, quiet=False):
        self.base = Path(base).resolve()
        self.out = self.base / OUT
        self.now = now or utcnow()
        self.today = self.now.astimezone(UTC).date()
        self.runner = runner or subprocess_runner
        self.have_lock = have_lock
        self.warnings = 0
        self.lines = []
        self.quiet = quiet
        self.cfg, self.raw, notes = load_settings(self.base)
        self.cron, cnotes = sw().cron_settings(self.base)
        self.setting_notes = notes + cnotes
        self.has_hs = False
        self.spec = None
        self.state = None

    def say(self, line):
        self.lines.append(line)
        if not self.quiet:
            print(line, flush=True)

    def note(self, msg):
        self.say("NOTE: weekly frames: " + msg)

    def warn(self, msg):
        self.warnings += 1
        self.say("WARNING: weekly frames: " + msg)

    @property
    def max_hs(self):
        return float(self.cron["DEM_MAX_HS"]) if (self.has_hs and self.cron.get("DEM_MAX_HS")) else None


def flock_file(path, wait_s=0.0, poll=5.0, create=True):
    """flock(LOCK_EX|LOCK_NB) on path, polled up to wait_s. -> file object, None when busy, or False when the
    file does not exist and create is False (no lock to take)."""
    import fcntl
    p = Path(path)
    if not p.exists() and not create:
        return False
    if create:
        p.parent.mkdir(parents=True, exist_ok=True)
        f = open(str(p), "a")
    else:
        f = open(str(p), "r")
    t_end = time.time() + wait_s
    while True:
        try:
            fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            return f
        except OSError:
            if time.time() >= t_end:
                f.close()
                return None
            time.sleep(min(poll, max(t_end - time.time(), 0.05)))


def unlock(f):
    if f:
        import fcntl
        try:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)
        finally:
            f.close()


def cron_due():
    """Waterline cron runs starting within CRON_HOLD_S of the NUC's clock."""
    try:
        return [h for h in sw().station_jobs_due(CRON_HOLD_S, local_now()) if "waterline_timex_cron" in h]
    except Exception:
        return []


def read_ground(ctx, d_min, d_max, copy=True, stats=False):
    """The pass over the ground file: under the cron's lock (unless --have-lock: the cron runs us and holds it).
    -> (info, why). copy: the rows of [d_min, d_max] into archive/weekly/work/ground_window.csv."""
    gpath = Path(ctx.cron["GROUND"])
    if not gpath.exists():
        return None, f"no ground file ({gpath.name})"
    lockf = None
    if not ctx.have_lock:
        hits = cron_due()
        if hits:
            return None, f"a waterline cron run is due ({', '.join(hits)} local): the ground file is not read now"
        lockf = flock_file(ctx.cron["LOCK"], LOCK_WAIT_S, create=False)
        if lockf is None:
            return None, f"the cron's lock ({Path(ctx.cron['LOCK']).name}) is held: the ground file is not read now"
    try:
        out = ctx.out / "work" / "ground_window.csv" if copy else None
        mh = float(ctx.cron["DEM_MAX_HS"]) if ctx.cron.get("DEM_MAX_HS") else None
        info = scan_ground(gpath, d_min, d_max, out, mh, stats=stats)
    except (OSError, ValueError) as exc:
        return None, f"the ground file could not be read ({str(exc)[:100]})"
    finally:
        unlock(lockf or None)
    info["copy"] = str(ctx.out / "work" / "ground_window.csv") if copy else None
    if stats and info.get("stats") is not None:
        info["stats"]["rows"] = info["rows"]
    ctx.has_hs = bool(info["has_hs"])
    return info, ""


# --------------------------------------------------------------------------------------------
# state
# --------------------------------------------------------------------------------------------

def load_state(ctx):
    st = read_json(ctx.out / "state.json")
    if not isinstance(st, dict):
        st = {"version": 1, "created_utc": iso_z(ctx.now)}
    st.setdefault("weeks", {})
    st.setdefault("videos", {})
    ctx.state = st
    return st


def save_state(ctx):
    if ctx.spec:
        ctx.state["spec_hash"] = ctx.spec.get("hash")
    write_json(ctx.out / "state.json", ctx.state)


# --------------------------------------------------------------------------------------------
# a week
# --------------------------------------------------------------------------------------------

def week_dir(ctx, monday):
    return ctx.out / "weeks" / monday


def week_stem(ctx, monday):
    return week_dir(ctx, monday) / f"wk_{monday}_{sunday_of(monday)}"


def decide_status(monday, today, latest, max_wait):
    """-> (status, reason): current, provisional or final ('covered' or 'waited')."""
    sun = d_of(sunday_of(monday))
    m = d_of(monday)
    if m <= today <= sun:
        return "current", ""
    if latest and latest > d_iso(sun):
        return "final", "covered"
    if (today - sun).days >= max_wait:
        return "final", "waited"
    return "provisional", ""


def week_fingerprint(ctx, info, monday):
    """The week's rows (counts, last date, and a digest of their content: the cron may re-process the ground file,
    e.g. a timing-mode switch, with the same rows but other elevations) and the DEM options."""
    sun = sunday_of(monday)
    days = {d: v for d, v in info["dates"].items() if monday <= d <= sun}
    rows = sum(v[0] for v in days.values())
    frames = sum(v[1] for v in days.values())
    last = max(days) if days else ""
    args = " ".join(sw().dem_args(ctx.cron, ctx.has_hs))
    dig = info.get("digest") or {}
    content = hashlib.sha1("|".join(f"{d}:{dig[d]}" for d in sorted(dig) if monday <= d <= sun).encode()) \
        .hexdigest()[:16]
    return f"{rows}|{frames}|{last}|{args}|{content}"


def dem_argv(ctx, ground, monday):
    """dem_from_contours.py over Monday..Sunday with exactly the cron's window-DEM options (storm_watch.dem_args:
    the cron's VAR= lines), --no-plot, and no --series-dir / --last-days (never archive/dems)."""
    return ([sys.executable, str(HERE / "dem_from_contours.py"), str(ground), str(week_stem(ctx, monday)),
             "--start-date", monday, "--end-date", sunday_of(monday)] + sw().dem_args(ctx.cron, ctx.has_hs)
            + ["--no-plot"])


def count_cells(stem):
    import numpy as np
    p = Path(str(stem) + "_dem.asc")
    if not p.exists():
        return 0
    try:
        return int(np.isfinite(read_asc(p)[0]).sum())
    except (OSError, ValueError):
        return 0


def frame_texts(spec, W):
    """(title, DEM subtitle, waterline subtitle, DEM stamp, waterline stamp) of a week record."""
    title = week_title(W["monday"], spec.get("station_name") or DEFAULT_STATION)
    fr = W.get("frames") or {}
    nfr = sum(fr.values()) if isinstance(fr, dict) else int(fr or 0)
    cams = ", ".join(f"{c} {n}" for c, n in sorted(fr.items())) if isinstance(fr, dict) and fr else ""
    label = spec.get("frame_label", "")
    if W.get("cells"):
        sub_dem = (f"DEM from {W.get('usable_days') or 0} usable day(s), {nfr} camera frames"
                   + (f" ({cams})" if cams else "") + f", {W['cells']} cells of {fmt_num(spec.get('cell') or 2)} m; "
                   f"the cron's window-DEM options over Monday..Sunday; elevations m NAVD88; {label}")
    else:
        sub_dem = f"no DEM: {W.get('why') or 'no data'}; elevations m NAVD88; {label}"
    mh = W.get("max_hs")
    sub_wl = (f"{W.get('lines') or 0} waterlines, {W.get('points') or 0} points from both cameras, in time order"
              + (f"; frames with offshore Hs above {fmt_num(mh)} m left out, as the DEM leaves them out"
                 if mh is not None else "")
              + f"; the DEM's day test is not applied to these points; elevations m NAVD88; {label}")
    st = W.get("status")
    thru = W.get("data_through") or "none yet"
    if st == "current":
        stamp = f"IN PROGRESS: data through {thru}"
    elif st == "provisional":
        stamp = f"PROVISIONAL: data through {thru}"
    elif W.get("final_reason") == "waited":
        stamp = f"data through {thru}"
    else:
        stamp = ""
    sd = stamp
    if not W.get("cells"):
        sd = (stamp + "\n" if stamp else "") + f"EMPTY: {W.get('why') or 'no DEM'}"
    sw_ = stamp
    if not W.get("points"):
        sw_ = (stamp + "\n" if stamp else "") + "EMPTY: no waterline points"
    return title, sub_dem, sub_wl, sd, sw_


def draw_week(ctx, W):
    """Both frames of a week from its kept grids and points.csv."""
    wd = week_dir(ctx, W["monday"])
    pts = read_points(wd / "points.csv", W["monday"], W["sunday"], None) if (wd / "points.csv").exists() else []
    W["points"] = len(pts)
    title, sub_dem, sub_wl, sd, swl = frame_texts(ctx.spec, W)
    stem = week_stem(ctx, W["monday"]) if W.get("cells") else None
    draw_dem_frame(stem, wd / FRAME_PNG["dem"], ctx.spec, title, sub_dem, sd)
    W["points"] = draw_waterline_frame(pts, wd / FRAME_PNG["waterlines"], ctx.spec, title, sub_wl, swl)
    W["spec_hash"] = ctx.spec.get("hash")
    W["drawn_utc"] = iso_z(ctx.now)


def publish_week(ctx, W):
    """current/ for the current week; frames/<series>/<sunday>_week.png for a final (or empty) one."""
    wd = week_dir(ctx, W["monday"])
    if W["status"] == "current":
        for s in SERIES:
            copy_atomic(wd / FRAME_PNG[s], ctx.out / "current" / FRAME_PNG[s])
    elif W["status"] in ("final", "empty"):
        for s in SERIES:
            copy_atomic(wd / FRAME_PNG[s], ctx.out / "frames" / s / f"{W['sunday']}_week.png")


def build_week(ctx, info, monday):
    """(Re)builds one week that is not final: its DEM when its rows changed, its points, both frames. -> record."""
    wd = week_dir(ctx, monday)
    sun = sunday_of(monday)
    prev = read_json(wd / "week.json", {}) or {}
    status, reason = decide_status(monday, ctx.today, info["latest"], ctx.cfg["WEEKLY_MAX_WAIT_DAYS"])
    fp = week_fingerprint(ctx, info, monday)
    stem = week_stem(ctx, monday)
    W = dict(prev) if prev else {}
    W.update(monday=monday, sunday=sun)
    days = {d: v for d, v in info["dates"].items() if monday <= d <= sun}
    latest_in = max(days) if days else None
    thru = min(info["latest"], sun) if info.get("latest") else None
    rebuilt = False
    if prev.get("fingerprint") != fp or prev.get("dem_rc") not in (0, 3, 4) or (
            prev.get("cells") and not Path(str(stem) + "_dem.asc").exists()):
        wd.mkdir(parents=True, exist_ok=True)
        for p in wd.glob(f"wk_{monday}_{sun}_*"):
            p.unlink()
        log = wd / "dem_log.txt"
        if log.exists():
            log.unlink()
        if days:
            rc, _ = ctx.runner(dem_argv(ctx, info["copy"], monday), log)
        else:
            rc = 4
        cells = count_cells(stem) if rc in (0, 3) else 0
        dinfo = read_json(str(stem) + "_info.json", {}) if rc in (0, 3) else {}
        why = ""
        if rc == 4:
            why = ("no waterline points in this week" if not days else
                   "no waterline points left, or the day test rejected every camera-day")
        elif rc not in (0, 3):
            why = f"dem_from_contours.py failed (exit {rc})"
        elif cells < MIN_CELLS:
            why = f"too few cells to map ({cells}, fewer than {MIN_CELLS})"
            cells = 0
        pts = read_points(info["copy"], monday, sun, ctx.max_hs) if days else []
        write_points(pts, wd / "points.csv")
        W.update(fingerprint=fp, dem_rc=rc, cells=cells, why=why, usable_days=dinfo.get("days") or 0,
                 frames=dinfo.get("frames") or {}, built_utc=iso_z(ctx.now),
                 lines=len({(p["camera"], p["source_file"] or p["capture_time_utc"]) for p in pts}),
                 dem_args=sw().dem_args(ctx.cron, ctx.has_hs), max_hs=ctx.max_hs)
        rebuilt = True
        if rc not in (0, 3, 4):
            ctx.warn(f"week {monday}: dem_from_contours.py failed (exit {rc}; see {wd.relative_to(ctx.base)}/"
                     "dem_log.txt): the week stays open and is built again by the next UTC day's daily run "
                     "(or now: `weekly_frames.py daily --base BASE --force`)")
    W.update(data_through=thru, latest_in_week=latest_in)
    if status == "final":
        if W.get("dem_rc") not in (0, 3, 4):
            status = "provisional"              # a failed build is never final: retried by the next day's run
        else:
            status = "final" if W.get("cells") else "empty"
            W["final_reason"] = reason
            W["finalised_utc"] = iso_z(ctx.now)
    W["status"] = status
    draw_week(ctx, W)
    write_json(wd / "week.json", W)
    publish_week(ctx, W)
    ctx.state["weeks"][monday] = status
    save_state(ctx)
    what = (f"DEM {W['cells']} cells, {W.get('usable_days') or 0} usable day(s)" if W.get("cells")
            else f"no DEM ({W.get('why')})")
    ctx.say(f"weekly frames: week {monday}..{sun} {status.upper()}"
            + (f" ({'covered' if reason == 'covered' else 'max wait'})" if status in ("final", "empty") else "")
            + f": {what}; {W.get('points') or 0} waterline points; data through {thru or 'none'}"
            + ("" if rebuilt else " (rows unchanged: redrawn only)"))
    return W


# --------------------------------------------------------------------------------------------
# the videos
# --------------------------------------------------------------------------------------------

WEEK_RE = re.compile(r"^(\d{4}-\d{2}-\d{2})_week\.png$")
PNG_NAME_RE = re.compile(r"^[A-Za-z0-9_.-]+\.png$")
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


def storm_frames(base):
    """The storm watch's frames (archive/storm_watch/events/*/frames/frames.json, read only), as
    [(series, date, sort name, path)]; an entry that does not check out is left out."""
    out, notes = [], []
    root = Path(base) / "archive" / "storm_watch" / "events"
    if not root.is_dir():
        return out, notes
    for ev in sorted(os.listdir(str(root))):
        if blind_named(ev) or ev.startswith("."):
            continue
        fj = root / ev / "frames" / "frames.json"
        if not fj.is_file():
            continue
        entries = read_json(fj)
        if not isinstance(entries, list):
            notes.append(f"{ev}/frames/frames.json is not a list: its frames are left out")
            continue
        for e in entries:
            if not isinstance(e, dict):
                continue
            series, d, png, part = e.get("series"), str(e.get("date") or ""), str(e.get("png") or ""), \
                str(e.get("part") or "")
            if series not in ("dem", "waterlines", "change") or not DATE_RE.match(d) or not PNG_NAME_RE.match(png) \
                    or d in BLIND_DATES or blind_named(png):
                notes.append(f"{ev}/frames/frames.json: an entry left out ({str(e)[:80]})")
                continue
            p = root / ev / "frames" / png
            if not p.is_file():
                notes.append(f"{ev}/frames/{png} is missing: left out")
                continue
            if series == "change" and not part.startswith("change"):
                part = "change_" + part
            out.append(("dem" if series == "change" else series, d, series == "change",
                        storm_frame_name(d, ev, part or Path(png).stem), p))
    return out, notes


def collect_frames(base, series):
    """The frames of a video in order: final weeks (frames/<series>/<sunday>_week.png) and the storm frames, by
    date (their first 10 characters), then by name; change frames after that date's DEM frames. -> (list of
    (sort name, path), notes)."""
    items = []
    d = Path(base) / OUT / "frames" / series
    if d.is_dir():
        for name in os.listdir(str(d)):
            m = WEEK_RE.match(name)
            if m and not blind_named(name):
                items.append((m.group(1), False, name, d / name))
    storm, notes = storm_frames(base)
    for s, dt, is_change, name, p in storm:
        if s == series:
            items.append((dt, is_change, name, p))
    items.sort(key=lambda t: (t[0], t[1], t[2]))
    return [(name, p) for _, _, name, p in items], notes


def concat_text(frames, fps, series):
    dur = 1.0 / float(fps)
    lines = ["ffconcat version 1.0", f"# weekly_frames.py: {series} video, {len(frames)} frame(s), "
             f"{fmt_num(fps)} fps"]
    for name, p in frames:
        q = str(Path(p).resolve()).replace("'", "'\\''")
        lines += [f"# {name}", f"file '{q}'", f"duration {fmt_num(dur)}"]
    return "\n".join(lines) + "\n"


def ffmpeg_argv(ffmpeg, list_path, out_tmp, spec):
    """Each frame at its slot k/fps: setpts snaps the concat list's times (k x 1/fps) to whole frames. Not N
    (the frame count): ffmpeg rebuilds its filters when the frame size changes (a storm frame drawn before a
    respec), N restarts at 0 there, and the frames after it would get earlier frames' times and be dropped."""
    w, h = spec["size"]
    fps = fmt_num(spec.get("fps") or 2)
    vf = (f"setpts=round(T*{fps})/({fps}*TB),scale={w}:{h}:force_original_aspect_ratio=decrease,"
          f"pad={w}:{h}:(ow-iw)/2:(oh-ih)/2:color=white,format=yuv420p")
    return [ffmpeg, "-hide_banner", "-loglevel", "error", "-y", "-f", "concat", "-safe", "0", "-i", str(list_path),
            "-vf", vf, "-r", fps, "-c:v", "libx264", "-pix_fmt", "yuv420p", "-movflags", "+faststart", "-f", "mp4",
            str(out_tmp)]


def build_videos(ctx, force=False):
    """Each series' MP4 from its frames, when their set or content changed since the last build (or force).
    Never a failure of the run (no WARNING, no exit 1): a missing ffmpeg is a NOTE (once per day), a failed ffmpeg
    a NOTE (the daily step tries once per day)."""
    vstate = ctx.state.setdefault("videos", {})
    ffmpeg = shutil.which("ffmpeg")
    noted = set()
    for series in SERIES:
        frames, notes = collect_frames(ctx.base, series)
        for n in notes:
            if n not in noted:
                noted.add(n)
                ctx.note(n)
        out = ctx.out / "video" / VIDEO_NAMES[series]
        if not frames:
            ctx.say(f"weekly frames: video {VIDEO_NAMES[series]}: no final frame yet")
            continue
        sig_items = []
        for name, p in frames:
            try:
                stt = os.stat(str(p))
                sig_items.append([name, str(p), int(stt.st_mtime_ns), int(stt.st_size)])
            except OSError:
                sig_items.append([name, str(p), 0, 0])
        sig = hashlib.sha256(json.dumps([fmt_num(ctx.spec.get("fps") or 2), ctx.spec["size"], sig_items])
                             .encode()).hexdigest()[:16]
        prev = vstate.get(series) or {}
        if not force and prev.get("sig") == sig and out.exists():
            continue
        lst = ctx.out / "video" / f"{series}_frames.txt"
        lst.parent.mkdir(parents=True, exist_ok=True)
        tmp_l = lst.with_name(lst.name + ".tmp")
        tmp_l.write_text(concat_text(frames, ctx.spec.get("fps") or 2, series))
        os.replace(str(tmp_l), str(lst))
        n_storm = sum(1 for name, _ in frames if "_storm_" in name)
        if not ffmpeg:
            if ctx.state.get("ffmpeg_note_date") != d_iso(ctx.today):
                ctx.state["ffmpeg_note_date"] = d_iso(ctx.today)
                ctx.note("ffmpeg not found: the videos are not made (the frames are; install it with "
                         "`sudo apt-get install -y ffmpeg`)")
            vstate[series] = dict(prev, pending_sig=sig, frames=len(frames))
            continue
        tmp = out.with_name("." + out.name + ".tmp")
        rc, text = ctx.runner(ffmpeg_argv(ffmpeg, lst, tmp, ctx.spec), None)
        if rc != 0 or not tmp.exists() or tmp.stat().st_size == 0:
            last = (text.strip().splitlines() or [""])[-1][:160]
            ctx.note(f"video {VIDEO_NAMES[series]} not rebuilt (ffmpeg exit {rc}: {last}); the frames are fine, "
                     "the earlier video is kept, and the next daily run tries again")
            vstate[series] = dict(prev, failed_sig=sig, failed_utc=iso_z(ctx.now), failed_rc=rc, failed_text=last)
            if tmp.exists():
                tmp.unlink()
            continue
        os.replace(str(tmp), str(out))
        vstate[series] = {"sig": sig, "built_utc": iso_z(ctx.now), "frames": len(frames), "storm_frames": n_storm,
                          "fps": ctx.spec.get("fps"), "path": str(out.relative_to(ctx.base))}
        ctx.say(f"weekly frames: video {VIDEO_NAMES[series]} rebuilt: {len(frames)} frame(s) "
                f"({len(frames) - n_storm} week(s), {n_storm} storm frame(s)) at {fmt_num(ctx.spec.get('fps') or 2)} fps")
    ctx.state["last_video_check_utc"] = iso_z(ctx.now)


# --------------------------------------------------------------------------------------------
# the commands
# --------------------------------------------------------------------------------------------

def take_own_lock(ctx):
    return flock_file(ctx.base / "logs" / "weekly_frames.lock", 0.0)


def clean_work(ctx):
    shutil.rmtree(str(ctx.out / "work"), ignore_errors=True)


def get_spec(ctx, info):
    stats = info.get("stats") if info else None
    if load_spec(ctx.base) is None and stats is None:
        return None
    if stats is not None:
        return ensure_spec(ctx.base, ctx.cron["GROUND"], ctx.cfg, log=ctx.say, stats=stats, cron=ctx.cron,
                           raw=ctx.raw)
    return _existing_spec(ctx)


def _existing_spec(ctx):
    spec = load_spec(ctx.base)
    for line in spec_differences(spec, ctx.cfg, ctx.raw, ctx.cron):
        ctx.note(line)
    return spec


def cmd_daily(ctx, force=False):
    if not ctx.cfg["WEEKLY_ENABLE"]:
        ctx.say("weekly frames: off (WEEKLY_ENABLE=0)")
        return 0
    if str(ctx.cron.get("DEM_ENABLE", "1")) != "1":
        ctx.say(f"weekly frames: off (the cron's DEM_ENABLE={ctx.cron.get('DEM_ENABLE')})")
        return 0
    today = d_iso(ctx.today)
    old = read_json(ctx.out / "state.json")
    if not force and isinstance(old, dict) and old.get("last_daily") == today:
        ctx.say(f"weekly frames: done today already ({today} UTC)")
        return 0
    for n in ctx.setting_notes:
        ctx.note(n)
    own = take_own_lock(ctx)
    if not own:
        ctx.say("weekly frames: another weekly_frames.py run is busy: skipped")
        return 0
    try:
        st = load_state(ctx)
        cur = d_iso(monday_of(ctx.today))
        open_weeks = sorted(m for m, s in st["weeks"].items() if s not in ("final", "empty"))
        todo = set(open_weeks) | {cur}
        known = sorted(st["weeks"])
        # missed weeks (the cron was off, a backfill stopped with --to): CATCHUP_WEEKS a run, oldest first, as
        # this runs inside the cron's lock; the rest wait in state.json's "catchup" for the next days' runs
        missed = set(st.get("catchup") or [])
        if known:
            m = monday_of(d_of(known[-1])) + timedelta(days=7)
            while d_iso(m) < cur:
                (todo if blind_in_week(m) else missed).add(d_iso(m))     # a blind week: only said skipped
                m += timedelta(days=7)
        missed = sorted(m for m in missed if DATE_RE.match(m) and m < cur and m not in st["weeks"]
                        and not blind_in_week(m))
        todo |= set(missed[:CATCHUP_WEEKS])
        st["catchup"] = missed                  # built ones drop out next run (they are in "weeks" then)
        rest = missed[CATCHUP_WEEKS:]
        if rest:
            ctx.note(f"{len(missed)} missed week(s): {CATCHUP_WEEKS} built now, {len(rest)} left for the next "
                     f"days' runs ({rest[0]} to {rest[-1]}); `weekly_frames.py backfill --base BASE` builds them "
                     "all at once")
        weeks = []
        for m in sorted(todo):
            b = blind_in_week(m)
            if b:
                ctx.say(f"weekly frames: week {m}..{sunday_of(m)} skipped: it contains the blind date {b}")
                continue
            weeks.append(m)
        if not weeks:
            st["last_daily"] = today
            save_state(ctx)
            return 0
        need_stats = load_spec(ctx.base) is None
        info, why = read_ground(ctx, weeks[0], today, copy=True, stats=need_stats)
        if info is None:
            ctx.note(why + ": nothing done (tried again next run)")
            return 0
        ctx.spec = get_spec(ctx, info)
        if ctx.spec is None:
            return 0
        for m in weeks:
            build_week(ctx, info, m)
        build_videos(ctx)
        st["last_daily"] = today
        st["last_daily_utc"] = iso_z(ctx.now)
        save_state(ctx)
        return 1 if ctx.warnings else 0
    finally:
        clean_work(ctx)
        unlock(own)


def cmd_backfill(ctx, d_from=None, d_to=None, dry_run=False):
    for n in ctx.setting_notes:
        ctx.note(n)
    last_ended = monday_of(ctx.today) - timedelta(days=7)
    to_m = last_ended if d_to is None else min(monday_of(d_of(d_to)), last_ended)
    if d_to is not None and monday_of(d_of(d_to)) > last_ended:
        ctx.note(f"--to {d_to}: only ended weeks are backfilled (up to the week of {d_iso(last_ended)})")
    own = None
    if not dry_run:
        own = take_own_lock(ctx)
        if not own:
            ctx.say("ERROR: another weekly_frames.py run is busy (logs/weekly_frames.lock): try again later")
            return 1
    try:
        d_lo = d_iso(monday_of(d_of(d_from))) if d_from else "0000-00-00"
        info, why = read_ground(ctx, d_lo, sunday_of(d_iso(to_m)), copy=not dry_run,
                                stats=(not dry_run and load_spec(ctx.base) is None))
        if info is None:
            ctx.say(f"ERROR: {why}")
            return 1
        if not info["first"]:
            ctx.say("weekly frames: backfill: the ground file has no georectified point: nothing to do")
            return 0
        first = monday_of(d_of(d_from)) if d_from else monday_of(d_of(info["first"]))
        st = read_json(ctx.out / "state.json") if dry_run else load_state(ctx)
        st = st if isinstance(st, dict) else {"weeks": {}}
        weeks = []
        m = first
        while m <= to_m:
            weeks.append(d_iso(m))
            m += timedelta(days=7)
        ctx.say(f"weekly frames: backfill {len(weeks)} week(s), {weeks[0] if weeks else '-'} to "
                f"{weeks[-1] if weeks else '-'}" + (" (dry run: nothing is written)" if dry_run else ""))
        if not dry_run:
            ctx.spec = get_spec(ctx, info)
            if ctx.spec is None:
                return 1
        for m in weeks:
            b = blind_in_week(m)
            if b:
                ctx.say(f"weekly frames: week {m}..{sunday_of(m)} skipped: it contains the blind date {b}")
                continue
            if (st.get("weeks") or {}).get(m) in ("final", "empty"):
                ctx.say(f"weekly frames: week {m}..{sunday_of(m)} already final: skipped")
                continue
            if dry_run:
                s, r = decide_status(m, ctx.today, info["latest"], ctx.cfg["WEEKLY_MAX_WAIT_DAYS"])
                n = sum(v[0] for d, v in info["dates"].items() if m <= d <= sunday_of(m))
                ctx.say(f"weekly frames: would build week {m}..{sunday_of(m)} ({s}{', ' + r if r else ''}; "
                        f"{n} ground row(s))")
                continue
            build_week(ctx, info, m)
        if not dry_run:
            build_videos(ctx)
            save_state(ctx)
        return 1 if ctx.warnings else 0
    finally:
        if not dry_run:
            clean_work(ctx)
        unlock(own)


def cmd_respec(ctx):
    for n in ctx.setting_notes:
        ctx.note(n)
    own = take_own_lock(ctx)
    if not own:
        ctx.say("ERROR: another weekly_frames.py run is busy (logs/weekly_frames.lock): try again later")
        return 1
    try:
        load_state(ctx)
        old = load_spec(ctx.base)
        auto = not (ctx.cfg["WEEKLY_EXTENT"] and ctx.cfg["WEEKLY_Z_LIMITS"])
        stats, rows = None, 0
        if auto:
            info, why = read_ground(ctx, "9999-12-31", "9999-12-31", copy=False, stats=True)
            if info is None:
                ctx.say(f"ERROR: {why}")
                return 1
            stats, rows = info["stats"], info["rows"]
        spec, why = make_spec(ctx.cfg, stats, ctx.cron, Path(ctx.cron["GROUND"]).name, rows, ctx.now)
        if spec is None:
            ctx.say(f"ERROR: no frame spec: {why}")
            return 1
        if old:
            write_json(ctx.out / "frame_spec.prev.json", old)
        write_json(spec_path(ctx.base), spec)
        ctx.spec = spec
        ctx.say(f"weekly frames: frame spec re-frozen ({(old or {}).get('hash', 'none')} -> {spec['hash']}): extent "
                f"{spec['extent']}, elevations {spec['z_limits']} m, {spec['size'][0]}x{spec['size'][1]} px, "
                f"{fmt_num(spec['fps'])} fps")
        n = 0
        wroot = ctx.out / "weeks"
        for m in sorted(os.listdir(str(wroot))) if wroot.is_dir() else []:
            W = read_json(wroot / m / "week.json")
            if not isinstance(W, dict) or not DATE_RE.match(m) or blind_in_week(m):
                continue
            draw_week(ctx, W)
            write_json(wroot / m / "week.json", W)
            publish_week(ctx, W)
            n += 1
        ctx.say(f"weekly frames: {n} week(s) redrawn")
        stale = []
        for s in SERIES:
            for name, p in collect_frames(ctx.base, s)[0]:
                if "_storm_" in name and png_text(p).get("weekly_frames_spec") != spec["hash"]:
                    stale.append(name)
        if stale:
            ctx.note(f"{len(stale)} storm frame(s) were drawn by storm_watch.py with the earlier look and keep it "
                     "(this script never writes into archive/storm_watch); the videos fit them to the frame size")
        build_videos(ctx, force=True)
        save_state(ctx)
        return 1 if ctx.warnings else 0
    finally:
        unlock(own)


def cmd_status(ctx):
    """Read only."""
    spec = load_spec(ctx.base)
    st = read_json(ctx.out / "state.json") or {}
    if spec is None and not st:
        print(f"weekly frames: nothing yet under {OUT}/ (the first daily run or a backfill makes it)")
        return 0
    if spec:
        print(f"frame spec {spec.get('hash')} (created {spec.get('created_utc')}): extent {spec.get('extent')} "
              f"({spec.get('extent_source')}), elevations {spec.get('z_limits')} m ({spec.get('z_source')}), "
              f"change +-{spec.get('change_limit')} m, {spec['size'][0]}x{spec['size'][1]} px, "
              f"{fmt_num(spec.get('fps') or 2)} fps, cell {spec.get('cell')} m, {spec.get('frame_label')}; "
              f"made from {spec.get('made_from')}")
        cfg, raw, _ = ctx.cfg, ctx.raw, None
        for line in spec_differences(spec, cfg, raw, ctx.cron):
            print("NOTE: " + line)
    else:
        print("frame spec: none yet")
    print(f"last daily run: {st.get('last_daily') or 'never'}")
    weeks = st.get("weeks") or {}
    for s in ("current", "provisional", "final", "empty"):
        ms = sorted(m for m, v in weeks.items() if v == s)
        print(f"{s:<12}: {len(ms)} week(s)" + (f": {', '.join(ms[-8:])}" + (" ..." if len(ms) > 8 else "")
                                              if ms else ""))
    for series in SERIES:
        v = (st.get("videos") or {}).get(series) or {}
        print(f"video {VIDEO_NAMES[series]}: " + (f"built {v.get('built_utc')}, {v.get('frames')} frame(s) "
                                                  f"({v.get('storm_frames', 0)} storm)" if v.get("built_utc")
                                                  else "not built yet"))
    print(f"ffmpeg: {shutil.which('ffmpeg') or 'not found (sudo apt-get install -y ffmpeg)'}")
    return 0


def parse_args(argv):
    ap = argparse.ArgumentParser(description="Weekly DEM and waterline frames and their videos "
                                             "(docs/WEEKLY_FRAMES.md).")
    sub = ap.add_subparsers(dest="cmd")
    for name, help_ in (("daily", "the cron's step: once per UTC date"),
                        ("backfill", "every ended week from the ground file's first date"),
                        ("respec", "re-freeze the look and redraw every frame"),
                        ("status", "read only: the look, the weeks, the videos")):
        p = sub.add_parser(name, help=help_)
        p.add_argument("--base", required=True, help="the waterline folder (BASE)")
        p.add_argument("--now", help=argparse.SUPPRESS)          # tests: the run's UTC time (ISO)
        if name in ("daily", "backfill", "respec"):
            p.add_argument("--have-lock", action="store_true",
                           help="the caller (the cron) holds the cron's lock: read the ground file directly")
        if name == "daily":
            p.add_argument("--force", action="store_true", help="do the work again today")
        if name == "backfill":
            p.add_argument("--from", dest="d_from", help="first date (YYYY-MM-DD; its week)")
            p.add_argument("--to", dest="d_to", help="last date (YYYY-MM-DD; its week, if it has ended)")
            p.add_argument("--dry-run", action="store_true", help="say what would be built; write nothing")
    a = ap.parse_args(argv)
    if not a.cmd:
        ap.print_help()
        sys.exit(1)
    for k in ("d_from", "d_to"):
        v = getattr(a, k, None)
        if v is not None:
            try:
                d_of(v)
            except ValueError:
                ap.error(f"--{k[2:]} {v!r}: not a date (YYYY-MM-DD)")
            if not DATE_RE.match(v):
                ap.error(f"--{k[2:]} {v!r}: not a date (YYYY-MM-DD)")
    return a


def main(argv=None, runner=None, now=None):
    a = parse_args(sys.argv[1:] if argv is None else argv)
    base = Path(a.base)
    if not base.is_dir():
        print(f"ERROR: --base {a.base}: not a folder")
        return 1
    _, bad = apply_station_file(base)
    if bad:
        print(f"ERROR: {bad}: nothing done")
        return 1
    if now is None and a.now:
        now = datetime.fromisoformat(a.now.replace("Z", "+00:00"))
        if now.tzinfo is None:
            now = now.replace(tzinfo=UTC)
    try:
        ctx = Ctx(base, now=now, runner=runner, have_lock=getattr(a, "have_lock", False))
    except SystemExit as exc:                        # station_env refused a value
        print(str(exc))
        return 1
    if a.cmd == "status":
        return cmd_status(ctx)
    try:
        if a.cmd == "daily":
            return cmd_daily(ctx, force=a.force)
        if a.cmd == "backfill":
            return cmd_backfill(ctx, a.d_from, a.d_to, a.dry_run)
        return cmd_respec(ctx)
    except Exception as exc:                         # never a traceback into the cron's log without a WARNING
        import traceback
        traceback.print_exc(file=sys.stdout)
        print(f"WARNING: weekly frames: {a.cmd} stopped ({type(exc).__name__}: {str(exc)[:160]})", flush=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())

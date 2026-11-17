#!/usr/bin/env python3
"""
Station Settings From The Environment
=======================================
The pipeline's built-in values are CACO-05's. A new station's own values
come from $BASE/waterline_station.conf (coastcam-station-setup writes it;
CACO-05 has none): the cron sources it, and its one export line puts these
names in the environment of every Python step it runs.

  name                     read by                                  CACO-05's value (the default)
  GNSSR_ANTENNA_NAVD88_M   extract_elevation_contours.py (and the    19.014 (m NAVD88, OPUS/GEOID18)
                           two scripts that import its constant)
  STATION_TZ               extract_elevation_contours.py,            America/New_York
                           dem_from_contours.py (local days)
  STATION_WAVE_CURRENCY    extract_elevation_contours.py,            marconi:hs_best/tp_s
                           dem_from_contours.py (the waves C was
                           fitted with)
  STATION_NAME             dem_figure.py, dem_map.py (titles)        each page's own CACO-05 title
  STATION_PROVISIONAL      dem_from_contours.py (_info.json),        none
                           dem_figure.py (the page)
  GEOTIFF_EPSG             asc_to_geotiff.py (--epsg), dem_map.py,   6348 (NAD83(2011) / UTM 19N)
                           dem_change.py (axis labels)
  GNSSR_QC_GROSS_M         gnssr_qc.py (gross range "lo,hi")         -2.5,3.5 (m NAVD88)
  WAVE_BUOY                gnssr_qc.py (waves_<buoy>.csv)            44008

Unset, or set to "" -> CACO-05's value, so with none of them set every
output is what it was. The one exception is GNSSR_ANTENNA_NAVD88_M: an
EMPTY antenna height is refused, like any other value that is not a
number, because CACO-05's 19.014 m would put another station's water
levels metres off without a word. A value that does not parse is never
replaced by CACO-05's: the script stops with an ERROR (SystemExit, exit
code 1) naming the variable. waterline_station_check.py check-overlay
refuses such a file before the cron sources it, so on a station these
errors are a last line of defence.

Python 3.8: no zoneinfo there; the time zone is checked with pytz or
dateutil instead (pandas brings one of them).
"""

import os
import re
import sys

# CACO-05's values: the defaults, and what every output was made with before these settings
CACO05_ANTENNA_NAVD88_M = 19.014
CACO05_TZ = "America/New_York"
CACO05_WAVE_CURRENCY = "marconi:hs_best/tp_s"
CACO05_EPSG = 6348
CACO05_QC_GROSS_M = (-2.5, 3.5)
CACO05_WAVE_BUOY = "44008"

_NUMBER = re.compile(r"^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$")
SOURCE = "waterline_station.conf's export line, written by coastcam-station-setup"


def _refuse(name, value, what):
    """Stops the script: `name`'s value does not parse. Never falls back to CACO-05's value."""
    raise SystemExit(f"ERROR: {name}={value!r} {what}. It comes from the environment ({SOURCE}); "
                     f"fix it there, or unset {name} for CACO-05's value.")


def _text(name):
    """The variable's value without surrounding blanks; '' when unset."""
    return os.environ.get(name, "").strip()


def _number(v):
    """float(v) for a plain decimal number ('19.014', '-2.5', '1e1'), else None (no 'nan', 'inf', '1m')."""
    if not _NUMBER.match(v):
        return None
    x = float(v)
    return x if abs(x) != float("inf") else None


def antenna_navd88_m():
    """GNSSR_ANTENNA_NAVD88_M: the GNSS antenna's height, m NAVD88. Unset: CACO-05's 19.014.
    Set but not a number (an empty value included), or outside -100..5000 m: refused."""
    if "GNSSR_ANTENNA_NAVD88_M" not in os.environ:
        return CACO05_ANTENNA_NAVD88_M
    raw = os.environ["GNSSR_ANTENNA_NAVD88_M"]
    x = _number(raw.strip())
    if x is None or not -100 < x < 5000:
        _refuse("GNSSR_ANTENNA_NAVD88_M", raw, "is not the antenna height in m NAVD88 (a number, e.g. 19.014)")
    return x


def _zone_known(name):
    """True/False: the time zone is known here (zoneinfo, else pytz, else dateutil); None: nothing to ask."""
    try:
        import zoneinfo                                  # Python 3.9+
    except ImportError:
        zoneinfo = None
    if zoneinfo is not None:
        try:
            zoneinfo.ZoneInfo(name)
            return True
        except Exception:                                # ZoneInfoNotFoundError; ValueError on odd names
            return False
    try:
        import pytz                                      # Python 3.8: pandas < 2 converts with pytz
    except ImportError:
        pytz = None
    if pytz is not None:
        try:
            pytz.timezone(name)
            return True
        except Exception:                                # UnknownTimeZoneError
            return False
    try:
        from dateutil import tz
    except ImportError:
        return None
    return tz.gettz(name) is not None


def station_tz():
    """STATION_TZ: the station's IANA time zone name (local days), default CACO-05's America/New_York.
    A name that is not a time zone here is refused. The name is returned (pandas converts with it)."""
    name = _text("STATION_TZ")
    if not name:
        return CACO05_TZ
    if name != CACO05_TZ:
        if not re.match(r"^[A-Za-z][A-Za-z0-9_+-]*(/[A-Za-z0-9_+-]+)*$", name) or _zone_known(name) is False:
            _refuse("STATION_TZ", name, "is not a time zone known here (an IANA name, e.g. America/New_York)")
    return name


def wave_currency():
    """STATION_WAVE_CURRENCY: the wave record the station's setup coefficient C was fitted with, as
    '<kind>:<height column>/<period column>' (extract_elevation_contours.wave_columns' tags); default
    CACO-05's marconi:hs_best/tp_s."""
    v = _text("STATION_WAVE_CURRENCY")
    if not v:
        return CACO05_WAVE_CURRENCY
    if not re.match(r"^[A-Za-z0-9_]+:[A-Za-z0-9_]+/[A-Za-z0-9_]+$", v):
        _refuse("STATION_WAVE_CURRENCY", v, "is not '<kind>:<height>/<period>' (e.g. buoy:wvht_m/dpd_s)")
    return v


def station_name(default):
    """STATION_NAME for titles; `default` (the page's own CACO-05 title) when unset or empty."""
    return _text("STATION_NAME") or default


def station_provisional():
    """STATION_PROVISIONAL: what is provisional about a new station's products (free text), or None."""
    return _text("STATION_PROVISIONAL") or None


def geotiff_epsg():
    """GEOTIFF_EPSG: the EPSG code of the grids' horizontal frame, default CACO-05's 6348."""
    v = _text("GEOTIFF_EPSG")
    if not v:
        return CACO05_EPSG
    if not re.match(r"^\d{4,6}$", v):
        _refuse("GEOTIFF_EPSG", v, "is not an EPSG code (digits, e.g. 6348)")
    return int(v)


def qc_gross_m():
    """GNSSR_QC_GROSS_M: gnssr_qc.py's gross range (lo, hi) in m NAVD88 from 'lo,hi'; default CACO-05's
    (-2.5, 3.5)."""
    v = _text("GNSSR_QC_GROSS_M")
    if not v:
        return CACO05_QC_GROSS_M
    parts = [p.strip() for p in v.split(",")]
    lo, hi = (_number(parts[0]), _number(parts[1])) if len(parts) == 2 else (None, None)
    if lo is None or hi is None or not lo < hi:
        _refuse("GNSSR_QC_GROSS_M", v, "is not 'lo,hi' in m NAVD88 with lo below hi (e.g. -2.5,3.5)")
    return lo, hi


def wave_buoy():
    """WAVE_BUOY: the NDBC station whose archive the cron keeps as waves_<buoy>.csv; default 44008."""
    v = _text("WAVE_BUOY")
    if not v:
        return CACO05_WAVE_BUOY
    if not re.match(r"^[A-Za-z0-9]+$", v):
        _refuse("WAVE_BUOY", v, "is not an NDBC station id (letters and digits, e.g. 44008)")
    return v


def frame_name(epsg):
    """The horizontal frame of an EPSG code, as the pages write it ('NAD83(2011) / UTM 19N' for 6348),
    or None for a code not in this short list (NAD83(2011), NAD83 and WGS 84 UTM zones)."""
    epsg = int(epsg)
    if 6330 <= epsg <= 6348:                             # NAD83(2011) / UTM zone 1N .. 19N
        return f"NAD83(2011) / UTM {epsg - 6329}N"
    if epsg in (6328, 6329):                             # NAD83(2011) / UTM zone 59N, 60N
        return f"NAD83(2011) / UTM {epsg - 6269}N"
    if 26901 <= epsg <= 26923:
        return f"NAD83 / UTM {epsg - 26900}N"
    if 32601 <= epsg <= 32660:
        return f"WGS 84 / UTM {epsg - 32600}N"
    if 32701 <= epsg <= 32760:
        return f"WGS 84 / UTM {epsg - 32700}S"
    return None


def frame_label(epsg=None, short=False):
    """The grids' frame for an axis label. Long: 'NAD83(2011) / UTM 19N, EPSG:6348'; short: 'UTM 19N'
    (the zone only). A code frame_name does not know: 'EPSG:<code>' either way. epsg None: GEOTIFF_EPSG."""
    epsg = geotiff_epsg() if epsg is None else int(epsg)
    name = frame_name(epsg)
    if name is None:
        return f"EPSG:{epsg}"
    return name.split(" / ", 1)[1] if short else f"{name}, EPSG:{epsg}"


def summary():
    """{name: value in use} for every setting (a station check prints it)."""
    return {"GNSSR_ANTENNA_NAVD88_M": antenna_navd88_m(), "STATION_TZ": station_tz(),
            "STATION_WAVE_CURRENCY": wave_currency(), "STATION_NAME": _text("STATION_NAME") or None,
            "STATION_PROVISIONAL": station_provisional(), "GEOTIFF_EPSG": geotiff_epsg(),
            "GNSSR_QC_GROSS_M": qc_gross_m(), "WAVE_BUOY": wave_buoy()}


if __name__ == "__main__":
    for k, v in summary().items():
        print(f"{k}={v}")
    sys.exit(0)

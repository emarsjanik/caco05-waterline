#!/usr/bin/env python3
"""
Marks For The Wave-Gauge Plots: Night, Screened Frames, The Model's Error
==========================================================================
Shared by owg_live.py (the 7-day plot of each model) and owg_report.py (the
emailed report), so the two draw the same thing the same way.

NIGHT. Shaded where the sun is below the horizon at Marconi (NOAA's solar
position equations; -0.833 deg, the sunrise/sunset of the almanacs). Worked
out from the clock, not from the frames, so it is right even when the camera
writes nothing at night. The camera's record stops here by design: the
network cannot see the waves in the dark.

SCREENED. A lighter shade over each daylight frame the quality screen
rejected because of the view -- too dark (dusk, a black squall), fog, rain or
a wet lens, sun glare, an unreadable image, the sea patch out of view. Each
frame covers half the frame interval either side of it (the median spacing
of the frames, 10-60 min), and touching frames merge into one band. A gap
in the camera's line with no shade under it is a gap in the images
themselves (camera or transfer down), not the weather.

ERROR BAND. Around the camera's 2-hour average, +/- the model's tested
error per frame: its holdback RMSE (frames it never saw in training), else
its validation RMSE; for several models averaged, the mean of theirs (the
average of two readings errs no more than their mean error). Per frame, so
for the 2-hour average it is on the safe side: a buoy line inside the band
agrees with the camera within what the camera can tell.
"""

import json
from pathlib import Path

import numpy as np

MARCONI_LAT, MARCONI_LON = 41.891, -69.962   # deg; to the km, ample for the sun to the minute
SUNSET_ELEV_DEG = -0.833                     # the sun's centre at sunrise/sunset (refraction + radius)
SCREENED = ("dark", "glare", "blurred (fog, rain or wet lens)", "unreadable", "patch not in view")
NIGHT_COLOUR, SCREENED_COLOUR = "#e9e8e2", "#f6eedb"
NIGHT_LABEL = "night"
SCREENED_LABEL = "view screened out (dusk, fog, wet lens, glare)"


def sun_elevation_deg(epoch, lat=MARCONI_LAT, lon=MARCONI_LON):
    """The sun's elevation (deg, no refraction) at UTC epoch(s) seen from lat, lon (deg, east positive):
    NOAA's solar calculator equations (Meeus), good to about 0.01 deg for this century."""
    ep = np.asarray(epoch, float)
    jc = (ep / 86400.0 + 2440587.5 - 2451545.0) / 36525.0
    l0 = np.radians((280.46646 + jc * (36000.76983 + jc * 0.0003032)) % 360.0)
    m = np.radians(357.52911 + jc * (35999.05029 - 0.0001537 * jc))
    e = 0.016708634 - jc * (0.000042037 + 0.0000001267 * jc)
    centre = (np.sin(m) * (1.914602 - jc * (0.004817 + 0.000014 * jc))
              + np.sin(2 * m) * (0.019993 - 0.000101 * jc) + np.sin(3 * m) * 0.000289)
    omega = np.radians(125.04 - 1934.136 * jc)
    app_long = np.radians(np.degrees(l0) + centre - 0.00569 - 0.00478 * np.sin(omega))
    obliq = np.radians(23.0 + (26.0 + (21.448 - jc * (46.815 + jc * (0.00059 - jc * 0.001813))) / 60.0) / 60.0
                       + 0.00256 * np.cos(omega))
    decl = np.arcsin(np.sin(obliq) * np.sin(app_long))
    y = np.tan(obliq / 2) ** 2
    eq_time = 4 * np.degrees(y * np.sin(2 * l0) - 2 * e * np.sin(m) + 4 * e * y * np.sin(m) * np.cos(2 * l0)
                             - 0.5 * y * y * np.sin(4 * l0) - 1.25 * e * e * np.sin(2 * m))      # minutes
    solar_min = ((ep % 86400.0) / 60.0 + eq_time + 4.0 * lon) % 1440.0
    hour_angle = np.radians(solar_min / 4.0 - 180.0)
    la = np.radians(lat)
    cos_zen = np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(hour_angle)
    return 90.0 - np.degrees(np.arccos(np.clip(cos_zen, -1.0, 1.0)))


def _runs(mask):
    """(first, last) index of each run of True."""
    m = np.r_[False, np.asarray(mask, bool), False]
    edges = np.flatnonzero(np.diff(m.astype(np.int8)))
    return list(zip(edges[0::2], edges[1::2] - 1))


def _merge(spans):
    out = []
    for a, b in sorted(spans):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def night_spans(t0, t1, step_s=120.0, lat=MARCONI_LAT, lon=MARCONI_LON):
    """[(start, end)] epochs where the sun is below SUNSET_ELEV_DEG between t0 and t1 (to step_s)."""
    t = np.arange(float(t0), float(t1) + step_s, step_s)
    below = sun_elevation_deg(t, lat, lon) < SUNSET_ELEV_DEG
    half = step_s / 2.0
    return [(max(t[a] - half, t0), min(t[b] + half, t1)) for a, b in _runs(below)]


def frame_spacing(epochs, default=1800.0):
    """The frames' usual interval (s): the median step between them, kept within 10-60 min."""
    ep = np.unique(np.asarray(epochs, float))
    d = np.diff(ep)
    d = d[d > 0]
    return float(np.clip(np.median(d), 600.0, 3600.0)) if len(d) else default


def screened_spans(epochs, statuses, t0, t1, lat=MARCONI_LAT, lon=MARCONI_LON):
    """[(start, end)] epochs covered by daylight frames screened out for the view (SCREENED), each +/- half
    the frame interval, touching ones merged, clipped to t0..t1. Night frames are left to the night shade."""
    ep = np.asarray(epochs, float)
    st = np.asarray(statuses, object)
    if not len(ep):
        return []
    half = frame_spacing(ep) / 2.0
    keep = np.isin(st, SCREENED) & (ep >= t0 - half) & (ep <= t1 + half)
    keep &= sun_elevation_deg(ep, lat, lon) >= SUNSET_ELEV_DEG
    spans = _merge([(e - half, e + half) for e in ep[keep]])
    return [(max(a, t0), min(b, t1)) for a, b in spans if b > t0 and a < t1]


def _dt(ep):
    return np.datetime64(int(round(float(ep))), "s")


def shade(ax, night=(), screened=(), labels=True):
    """Night and screened spans as background bands, under the grid and the data."""
    for i, (a, b) in enumerate(night):
        ax.axvspan(_dt(a), _dt(b), color=NIGHT_COLOUR, lw=0, zorder=0,
                   label=NIGHT_LABEL if labels and i == 0 else None)
    for i, (a, b) in enumerate(screened):
        ax.axvspan(_dt(a), _dt(b), color=SCREENED_COLOUR, lw=0, zorder=0,
                   label=SCREENED_LABEL if labels and i == 0 else None)


def background_last(handles, labels):
    """Legend entries with the night and screened shades moved after the data."""
    back = (NIGHT_LABEL, SCREENED_LABEL)
    pairs = [hl for hl in zip(handles, labels) if hl[1] not in back] + \
        [hl for b in back for hl in zip(handles, labels) if hl[1] == b]
    return [h for h, _ in pairs], [lab for _, lab in pairs]


def error_band(ax, pieces, rmse, colour, label=True):
    """+/- rmse around each (times, values) piece of an average, floored at 0 m."""
    if rmse is None or not np.isfinite(rmse) or rmse <= 0:
        return
    for i, (ts, vs) in enumerate(pieces):
        vs = np.asarray(vs, float)
        ax.fill_between(ts, np.maximum(vs - rmse, 0.0), vs + rmse, color=colour, alpha=0.14, lw=0, zorder=1,
                        label=f"tested error per frame, ±{rmse:.2f} m" if label and i == 0 else None)


def model_rmse(stems):
    """The tested error per frame of a model, or of several averaged (the mean of theirs): each report's
    holdback_rmse, else its val_rmse. None when any report or value is missing."""
    vals = []
    for s in stems:
        p = Path(str(s) + ".report.json") if s else None
        if p is None or not p.is_file():
            return None
        try:
            rep = json.loads(p.read_text())
        except (OSError, ValueError):
            return None
        v = rep.get("holdback_rmse")
        v = rep.get("val_rmse") if v is None else v
        try:
            v = float(v)
        except (TypeError, ValueError):
            return None
        if not np.isfinite(v) or v <= 0:
            return None
        vals.append(v)
    return float(np.mean(vals)) if vals else None

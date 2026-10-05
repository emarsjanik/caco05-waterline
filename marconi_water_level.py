#!/usr/bin/env python3
"""
Still-Water Level At Marconi: GPS First, Chatham As Backup
============================================================
One water-level source for the camera products, in this order:

  1. GNSS-R (the station's own GPS reflectometry, on the camera tower):
     the gnssrefl spline on NAVD88, quality-checked by gnssr_qc.py. It
     measures the water exactly where the cameras look. It arrives a day
     or two late (precise orbits), so the newest frames usually fall after
     its last reading.
  2. The NOAA Chatham gauge (20 km away, inside a harbour) converted to
     Marconi as a * Chatham(t - lag) + b, with a, lag and b FITTED TO THE
     GPS over the last 30 days every time this runs (gnssr_qc.fit_reference),
     not fixed numbers. Default (1.24, -48 min, -0.10) only if there is no
     GPS overlap. Covers everything up to the gauge's last reading.
  3. Past the gauge's last reading, up to 2 h: a tide fit to the last 25 h
     carried forward (as owg_live.py did), pinned to the last reading.

Each level comes with its source ("gps", "chatham", "tide fit"), so a
product can say what it used.

Library:
    wl = WaterLevel()                   # station default paths
    levels, sources = wl.at(epochs)
    print(wl.describe())

Command line (prints the current calibration and recent levels):
    python3 marconi_water_level.py
"""

import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

GNSSR_SPLINE = "/home/argus_user/GNSS/v4.1/products/refl_code/Files/usgs/usgs_spline_out.txt"
GAUGE_CSV = str(HERE / "archive" / "gauge_8447435.csv")
DEFAULT_FIT = (1.24, -48 * 60.0, -0.10, 0.13)        # a, lag (s), b, sigma: Sep 2026 fit


class WaterLevel:
    def __init__(self, spline=GNSSR_SPLINE, gauge_csv=GAUGE_CSV, ahead_hours=2.0, max_gap_s=3600,
                 use_gps=True):
        self.ahead = ahead_hours
        self.max_gap = max_gap_s
        self.g_ep = self.g_lv = None
        if gauge_csv and Path(gauge_csv).exists():
            g = pd.read_csv(gauge_csv).dropna(subset=["level_navd88"]).sort_values("epoch")
            self.g_ep, self.g_lv = g["epoch"].to_numpy(float), g["level_navd88"].to_numpy(float)
        self.s_ep = self.s_lv = None
        self.gps_note = "GPS not used"
        if use_gps and spline and Path(spline).exists():
            try:
                import io, contextlib
                from extract_elevation_contours import load_gnssr_spline
                from gnssr_qc import qc_filter
                with contextlib.redirect_stdout(io.StringIO()):       # the loader is chatty
                    ep, lv, _, _ = load_gnssr_spline(spline)
                    if self.g_ep is not None:
                        ep, lv, _ = qc_filter(np.asarray(ep, float), np.asarray(lv, float), gauge_csv)
                self.s_ep, self.s_lv = np.asarray(ep, float), np.asarray(lv, float)
                self.gps_note = (f"GPS {pd.to_datetime(self.s_ep[0], unit='s'):%Y-%m-%d} .. "
                                 f"{pd.to_datetime(self.s_ep[-1], unit='s'):%Y-%m-%d %H:%M}Z")
            except Exception as exc:                       # never let the GPS stop the cameras
                self.gps_note = f"GPS unreadable ({exc})"
        self.fit, self.fit_note = DEFAULT_FIT, "default (no GPS overlap)"
        if self.s_ep is not None and self.g_ep is not None:
            try:
                from gnssr_qc import fit_reference
                f = fit_reference(self.s_ep, self.s_lv, self.g_ep, self.g_lv, np.ones(len(self.s_ep), bool))
                if f:
                    self.fit = f
                    self.fit_note = "fitted to the GPS, last 30 days"
            except Exception as exc:
                self.fit_note = f"default (fit failed: {exc})"

    # -- sources
    def _gps(self, t):
        out = np.full(len(t), np.nan)
        if self.s_ep is None or len(self.s_ep) < 2:
            return out
        ep, lv = self.s_ep, self.s_lv
        i = np.clip(np.searchsorted(ep, t), 1, len(ep) - 1)
        ok = (t >= ep[0]) & (t <= ep[-1]) & (ep[i] - ep[i - 1] <= self.max_gap)
        out[ok] = np.interp(t[ok], ep, lv)
        return out

    def _chatham(self, t):
        out = np.full(len(t), np.nan)
        if self.g_ep is None or len(self.g_ep) < 2:
            return out
        a, lag, b, _ = self.fit
        tt = t - lag
        ep, lv = self.g_ep, self.g_lv
        i = np.clip(np.searchsorted(ep, tt), 1, len(ep) - 1)
        ok = (tt >= ep[0]) & (tt <= ep[-1]) & (ep[i] - ep[i - 1] <= 3 * 3600)
        out[ok] = a * np.interp(tt[ok], ep, lv) + b
        return out

    def _tide_ahead(self, t):
        """Past the gauge's last reading: M2+M4 fit to the last 25 h, pinned to the last reading."""
        out = np.full(len(t), np.nan)
        if self.g_ep is None or not self.ahead:
            return out
        a, lag, b, _ = self.fit
        ep, lv = self.g_ep, self.g_lv
        m = ep > ep[-1] - 25 * 3600
        if m.sum() < 150:
            return out
        tt = t - lag
        ok = (tt > ep[-1]) & (tt - ep[-1] <= self.ahead * 3600)
        if not ok.any():
            return out
        w = 2 * np.pi / (12.4206 * 3600)
        basis = lambda x: np.column_stack([np.ones_like(x), np.cos(w * x), np.sin(w * x),
                                           np.cos(2 * w * x), np.sin(2 * w * x)])
        c, *_ = np.linalg.lstsq(basis(ep[m] - ep[-1]), lv[m], rcond=None)
        f = basis(np.r_[0.0, tt[ok] - ep[-1]]) @ c
        out[ok] = a * (f[1:] + lv[-1] - f[0]) + b
        return out

    def at(self, epochs):
        """-> (levels m NAVD88, sources) for each epoch; NaN / '' where nothing covers it."""
        t = np.atleast_1d(np.asarray(epochs, float))
        lv = self._gps(t)
        src = np.where(np.isfinite(lv), "gps", "").astype(object)
        for fn, name in ((self._chatham, "chatham"), (self._tide_ahead, "tide fit")):
            need = ~np.isfinite(lv)
            if not need.any():
                break
            v = fn(t[need])
            idx = np.flatnonzero(need)[np.isfinite(v)]
            lv[idx] = v[np.isfinite(v)]
            src[idx] = name
        return lv, src

    def describe(self):
        a, lag, b, sig = self.fit
        return (f"water level: {self.gps_note}; backup Chatham x {a:.3f} at {lag / 60:+.0f} min "
                f"{b:+.3f} m (sigma {sig:.3f} m, {self.fit_note})")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--spline", default=GNSSR_SPLINE)
    ap.add_argument("--gauge-csv", default=GAUGE_CSV)
    ap.add_argument("--hours", type=float, default=72, help="show the last N hours (default 72)")
    args = ap.parse_args()
    wl = WaterLevel(args.spline, args.gauge_csv)
    print(wl.describe())
    now = pd.Timestamp.now(tz="UTC").timestamp()
    t = np.arange(now - args.hours * 3600, now, 3 * 3600.0)
    lv, src = wl.at(t)
    for e, v, s in zip(t, lv, src):
        print(f"  {pd.to_datetime(e, unit='s'):%Y-%m-%d %H:%M}Z  {v:+.2f} m  {s}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

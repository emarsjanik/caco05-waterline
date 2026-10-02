#!/usr/bin/env python3
"""
Offshore Buoy To Marconi: Direction-Aware Wave Height Transfer
================================================================
Fits how an offshore buoy's wave height turns into the wave height off
Marconi (the ADCP in 21 m of water, Dec 2024 - Mar 2025, the depth the
camera models learned and the USGS total-water-level forecast uses), as a
function of wave DIRECTION, period and height:

    log Hs_marconi = c0 + c1 cos(D) + c2 sin(D) + c3 cos(2D) + c4 sin(2D)
                     + c5 log(Tp) + c6 log(Hs_buoy)

WHY. Marconi faces ~80 deg. Waves from the south and west arrive at about
half the offshore height (sheltering), so a raw buoy says nothing fair
about the beach on those days. Scored by leaving out whole weeks, against
the ADCP:
    NDBC 44013 (Boston)  raw 0.33 m RMSE -> corrected 0.27 m
    WIS hindcast 63064   raw 0.46 m      -> corrected 0.29 m
and averaged 50:50 with the camera (errors correlate only 0.34), the
per-frame error on the camera's validation frames fell 0.288 -> 0.229 m;
for waves above 1.5 m the corrected buoy (0.38 m) beats the camera (0.65 m).

FIT (on the station: it downloads the buoy's NDBC history for the ADCP months)
    python3 buoy_transfer.py fit --station 44008
    python3 buoy_transfer.py fit --station 44008 --input 44008h2024.txt.gz --input 44008h2025.txt.gz
writes calibration/buoy_transfer_44008.json (coefficients, weekly
cross-validated RMSE). owg_report.py then plots the corrected buoy, scores
the camera against it, and blends the two into a best estimate.

APPLY (library): load(path) -> transfer; apply(transfer, hs, tp, dir) -> Hs at Marconi.
"""

import sys
import gzip
import json
import argparse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
HIST = "https://www.ndbc.noaa.gov/data/historical/stdmet/{st}h{yr}.txt.gz"


def features(hs, tp, d):
    th = np.radians(np.asarray(d, float))
    return np.column_stack([np.ones(len(th)), np.cos(th), np.sin(th), np.cos(2 * th), np.sin(2 * th),
                            np.log(np.asarray(tp, float)), np.log(np.asarray(hs, float))])


def apply(transfer, hs, tp, d):
    """Buoy (Hs m, peak period s, mean direction deg from) -> Hs at Marconi (m); NaN where unknown."""
    hs, tp, d = (np.asarray(v, float) for v in (hs, tp, d))
    out = np.full(len(hs), np.nan)
    ok = np.isfinite(hs) & np.isfinite(tp) & np.isfinite(d) & (hs > 0.05) & (tp > 1)
    # only directions the fit saw (>= 10 h in their 45 deg sector): the harmonics
    # are meaningless outside them
    seen = transfer.get("median_ratio_by_direction", {})
    sector = (np.nan_to_num(d) // 45 * 45).astype(int)
    ok &= np.fromiter((f"{s}-{s + 45}" in seen for s in sector), dtype=bool, count=len(sector))
    if ok.any():
        out[ok] = np.exp(features(hs[ok], tp[ok], d[ok]) @ np.array(transfer["coef"]))
        out[ok] = np.clip(out[ok], 0.1 * hs[ok], 2.0 * hs[ok])      # never wild
    return out


def load(path):
    p = Path(path)
    return json.loads(p.read_text()) if p.exists() else None


def parse_ndbc(text):
    """NDBC standard meteorological text (realtime2 or historical) -> hourly Hs, Tp, MWD."""
    rows = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        p = line.split()
        if len(p) < 12:
            continue
        try:
            yr = int(p[0]); yr = yr + 1900 if yr < 100 else yr
            t = pd.Timestamp(year=yr, month=int(p[1]), day=int(p[2]), hour=int(p[3]), minute=int(p[4]), tz="UTC")
        except ValueError:
            continue
        v = []
        for i in (8, 9, 11):
            try:
                x = float(p[i])
            except ValueError:
                x = np.nan
            missing = x > 360 if i == 11 else x >= 99        # NDBC: 99.00 heights/periods, 999 directions
            v.append(np.nan if missing else x)
        rows.append((t, *v))
    # keep rows without a height: NDBC reports direction on rows of its own
    df = pd.DataFrame(rows, columns=["t", "hs", "tp", "dir"]).set_index("t")
    return hourly(df).dropna(subset=["hs"])


def hourly(df):
    """Hourly means; direction as a circular mean."""
    s, c = np.sin(np.radians(df["dir"])), np.cos(np.radians(df["dir"]))
    h = pd.DataFrame({"hs": df["hs"], "tp": df["tp"], "s": s, "c": c}).resample("1h").mean()
    h["dir"] = np.degrees(np.arctan2(h["s"], h["c"])) % 360
    return h[["hs", "tp", "dir"]]


def read_buoy_csv(path):
    """fetch_buoy_waves.py archive (time_utc, epoch, wvht_m, dpd_s, apd_s, mwd_deg)."""
    w = pd.read_csv(path)
    df = pd.DataFrame({"t": pd.to_datetime(w["epoch"], unit="s", utc=True),
                       "hs": pd.to_numeric(w["wvht_m"], errors="coerce"),
                       "tp": pd.to_numeric(w["dpd_s"], errors="coerce"),
                       "dir": pd.to_numeric(w["mwd_deg"], errors="coerce")}).set_index("t")
    return hourly(df).dropna(subset=["hs"])


def fit(args):
    a = pd.read_csv(args.adcp, parse_dates=["time"])
    a["t"] = pd.to_datetime(a["time"], utc=True)
    a = a.set_index("t")["wh_4061"].rename("A")
    parts = []
    for src in args.input or []:
        raw = Path(src).read_bytes()
        text = gzip.decompress(raw).decode() if src.endswith(".gz") else raw.decode()
        parts.append(parse_ndbc(text))
    if args.buoy_csv:
        parts.append(read_buoy_csv(args.buoy_csv))
    if not parts:
        years = sorted({a.index.min().year, a.index.max().year})
        for yr in years:
            url = HIST.format(st=args.station, yr=yr)
            try:
                raw = urllib.request.urlopen(url, timeout=120).read()
            except Exception as e:
                sys.exit(f"download failed: {url}: {e}\n(download it by hand and pass --input FILE)")
            print(f"downloaded {url} ({len(raw) / 1e6:.1f} MB)")
            parts.append(parse_ndbc(gzip.decompress(raw).decode()))
    b = pd.concat(parts).sort_index()
    b = b[~b.index.duplicated()]
    j = pd.concat([b, a], axis=1, sort=True).dropna()
    j = j[(j["hs"] > 0.05) & (j["A"] > 0.05) & (j["tp"] > 1)]
    if len(j) < 200:
        sys.exit(f"only {len(j)} hours where buoy and ADCP overlap -- check the buoy files cover "
                 f"{a.index.min():%Y-%m-%d} .. {a.index.max():%Y-%m-%d}")

    X, y = features(j["hs"], j["tp"], j["dir"]), np.log(j["A"].to_numpy())
    week = j.index.tz_convert(None).to_period("W")
    pred = np.full(len(j), np.nan)
    for g in week.unique():                                   # leave whole weeks out
        te = np.asarray(week == g)
        c, *_ = np.linalg.lstsq(X[~te], y[~te], rcond=None)
        pred[te] = np.exp(X[te] @ c)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = lambda e: float(np.sqrt(np.mean(e ** 2)))
    raw_rmse, cv_rmse = r(j["hs"] - j["A"]), r(pred - j["A"])
    sectors = {}
    for lo in range(0, 360, 45):
        m = (j["dir"] >= lo) & (j["dir"] < lo + 45)
        if m.sum() >= 10:
            sectors[f"{lo}-{lo + 45}"] = round(float((j["A"][m] / j["hs"][m]).median()), 2)
    out = {"station": args.station, "coef": coef.tolist(),
           "terms": ["1", "cos D", "sin D", "cos 2D", "sin 2D", "log Tp", "log Hs"],
           "n_hours": int(len(j)), "period": [str(j.index.min()), str(j.index.max())],
           "raw_rmse": round(raw_rmse, 3), "cv_rmse": round(cv_rmse, 3),
           "cv_bias": round(float(np.mean(pred - j["A"])), 3),
           "cv_r": round(float(np.corrcoef(pred, j["A"])[0, 1]), 3),
           "median_ratio_by_direction": sectors,
           "reference": "ADCP wh_4061, Signature 1000 in 21 m off Marconi (USGS doi:10.5066/P13KF7UW)"}
    dst = Path(args.output or HERE / "calibration" / f"buoy_transfer_{args.station}.json")
    dst.write_text(json.dumps(out, indent=1) + "\n")
    print(f"buoy {args.station} vs ADCP, {len(j)} hours: raw RMSE {raw_rmse:.2f} m -> corrected "
          f"{cv_rmse:.2f} m (whole weeks left out), r {out['cv_r']:.2f}")
    print("median Marconi/buoy ratio by direction (from): " +
          ", ".join(f"{k}: {v}" for k, v in sectors.items()))
    print(f"wrote {dst}")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fit")
    f.add_argument("--station", default="44008")
    f.add_argument("--adcp", default=str(HERE / "sig1000_waves_ALL.csv"))
    f.add_argument("--input", action="append", help="NDBC stdmet file(s), .txt or .txt.gz")
    f.add_argument("--buoy-csv", help="or a fetch_buoy_waves.py archive covering the ADCP months")
    f.add_argument("--output")
    args = ap.parse_args()
    return fit(args)


if __name__ == "__main__":
    sys.exit(main())

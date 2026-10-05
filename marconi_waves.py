#!/usr/bin/env python3
"""
Best Estimate Of The Waves Off Marconi, Hourly
================================================
One wave record for every system that needs one (OWG report, GNSS-IR
wave-setup checks, DEM setup corrections), built from everything we have:

  * the camera: the frame-by-frame mean of the wave models listed in
    --models (ensemble), averaged within each hour;
  * the offshore buoy, converted to Marconi by buoy_transfer.py (wave
    direction, period and height) when calibration/buoy_transfer_<st>.json
    exists, otherwise not used for the blend;
  * blended by inverse error variance where both exist. Their errors are
    only weakly related (r 0.34 on the validation frames), so the blend
    beats either: 0.288 (camera) / 0.278 (buoy) -> 0.229 m per frame.

The camera's error is not one number: it grows with wave height (shrinkage
toward the training mean). On the validation frames, by the camera's own
reading: < 0.6 m 0.16 m RMSE, 0.6-1.0 m 0.23, 1.0-1.4 m 0.29, > 1.4 m 0.65
(--camera-sigma). The two frames of an hour share most of their error
(r 0.80), so averaging them gains little (--frame-corr). In storms the
buoy therefore carries the blend. hs_sigma is the 1-sigma uncertainty of
hs_best from these errors.

Hours: a frame belongs to the hour it starts in (10:00 and 10:30 -> 10:00),
like the buoy's hourly means.

Period and direction come from the buoy (a single camera frame cannot
measure period -- Chris Sherwood's model scored R2 -0.19 at it).

OUTPUT archive/waves_marconi.csv, one row per hour of buoy record:
  time_utc, epoch, hs_buoy, hs_buoy_marconi, hs_camera, n_frames,
  hs_best, hs_sigma, best_from (blend / camera / buoy / buoy (raw)),
  n_models, qc_flag (1 pass, 3 suspect, 4 fail: IOOS QARTOD convention), qc_note
  (which tests: gross range, camera above training range, camera-buoy disagree,
  spike, flat line), tp_s (buoy PEAK period only), tp_camera_s (camera MEAN period
  Tm01 from the timestacks, timestack_wave_period.py -- a different
  statistic, ~0.7-0.85 x Tp, so never mixed into tp_s), dir_deg

Usage:
    python3 marconi_waves.py                 (owg.sh runs it every hour)
"""

import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import buoy_transfer   # noqa: E402


def camera_sigma(hs, table):
    """Per-frame camera error (m) at the camera's reading, from 'HS:SIGMA,...' (linear between)."""
    pts = np.array([[float(x) for x in t.split(":")] for t in table.split(",")])
    return np.interp(np.asarray(hs, float), pts[:, 0], pts[:, 1])


TRAINING_MAX_HS = 3.2     # largest wave height in the camera models' training labels (ADCP, 3.16 m)


def qc_flags(res, sc, sb):
    """IOOS QARTOD-style flags for each hour: 1 pass, 3 suspect, 4 fail (9 missing rows are not written).

    gross range   hs_best < 0 or > 8 m                                  -> 4
    extrapolation camera reading above the training range               -> 3
    agreement     camera and converted buoy differ by > 3 sigma         -> 3
    spike         hs_best differs from the mean of the hours either side by > 0.6 m -> 3
    flat line     camera reading unchanged (to the cm) for 4 hours or more          -> 3"""
    hs, cam, buo = res["hs_best"].to_numpy(float), res["hs_camera"].to_numpy(float), res["hs_buoy_marconi"].to_numpy(float)
    flag = np.ones(len(res), int)
    note = np.full(len(res), "", object)

    def mark(m, f, why):
        m = np.asarray(m, bool)
        flag[m] = np.maximum(flag[m], f)
        note[m] = [f"{n};{why}" if n else why for n in note[m]]
    mark((hs < 0) | (hs > 8), 4, "gross range")
    mark(cam > TRAINING_MAX_HS, 3, "camera above training range")
    with np.errstate(invalid="ignore"):
        mark(np.abs(cam - buo) > 3 * np.sqrt(sc ** 2 + sb ** 2), 3, "camera-buoy disagree")
    t = res["epoch"].to_numpy()
    prev_ok = np.r_[False, np.diff(t) == 3600]
    next_ok = np.r_[np.diff(t) == 3600, False]
    both = prev_ok & next_ok
    nb = np.full(len(res), np.nan)
    nb[both] = (np.r_[np.nan, hs[:-1]][both] + np.r_[hs[1:], np.nan][both]) / 2
    with np.errstate(invalid="ignore"):
        mark(np.abs(hs - nb) > 0.6, 3, "spike")
    run = (pd.Series(cam).round(2).diff() == 0) & pd.Series(prev_ok)
    grp = (~run).cumsum()
    length = run.groupby(grp).transform("sum").to_numpy()
    mark(run.to_numpy() & (length >= 3), 3, "flat line")
    return flag, note


def camera_hourly(archive, specs):
    """Frame-by-frame mean across models, then the mean of the frames in each hour.

    A frame is used only when every model whose archive spans its time has a
    value for it, so an hour never switches between one model alone and the
    mean (their biases differ). A model installed later does not remove the
    earlier frames: before its first frame it is not required."""
    frames = []
    for spec in specs:
        label, name = spec.split("=")[:2]
        p = Path(archive) / name
        if not p.exists():
            continue
        d = pd.read_csv(p)
        d = d[d["status"] == "ok"][["filename", "epoch", "hs_m"]].rename(columns={"hs_m": label})
        frames.append(d)
    if not frames:
        return pd.DataFrame(columns=["hs_camera", "n_frames"])
    e = frames[0]
    for d in frames[1:]:
        e = e.merge(d.drop(columns="epoch"), on="filename", how="outer")
    e["epoch"] = e["filename"].str.split(".").str[0].astype(int)
    labels = [c for c in e.columns if c not in ("filename", "epoch")]
    need = np.ones(len(e), bool)
    for c in labels:
        have = e[c].notna()
        span = (e["epoch"] >= e.loc[have, "epoch"].min()) & (e["epoch"] <= e.loc[have, "epoch"].max())
        need &= ~(span & ~have).to_numpy()           # inside this model's span but missing: drop
    e = e[need].copy()
    e["hs"] = e[labels].mean(axis=1)
    e["n_models"] = e[labels].notna().sum(axis=1)
    e = e[e["hs"].notna()]
    e["hour"] = pd.to_datetime((e["epoch"] // 3600) * 3600, unit="s", utc=True)
    g = e.groupby("hour")
    return pd.DataFrame({"hs_camera": g["hs"].mean(), "n_frames": g.size(), "n_models": g["n_models"].min()})


def load_buoy(path):
    w = pd.read_csv(path)
    t = pd.to_datetime(w["epoch"], unit="s", utc=True)
    b = pd.DataFrame({"hs": pd.to_numeric(w["wvht_m"], errors="coerce").to_numpy(),
                      "tp": pd.to_numeric(w["dpd_s"], errors="coerce").to_numpy(),
                      "dir": pd.to_numeric(w["mwd_deg"], errors="coerce").to_numpy()}, index=t)
    b = buoy_transfer.hourly(b).dropna(subset=["hs"])
    b[["tp", "dir"]] = b[["tp", "dir"]].ffill(limit=3)      # NDBC often blanks them for an hour
    return b


def build(args):
    # every buoy with an archive; those with a fitted transfer are converted to
    # Marconi and averaged, weighted by their cross-validated error
    raw, conv, w_sum = None, 0.0, 0.0
    converted = []
    for spec in args.buoys:
        st, name = spec.split("=")
        p = Path(args.archive) / name
        if not p.exists():
            continue
        b = load_buoy(p)
        if raw is None:
            raw = b.copy()                                   # first buoy listed: the raw reference shown
        tr = buoy_transfer.load(Path(args.calibration) / f"buoy_transfer_{st}.json")
        if not tr:
            continue
        h = pd.Series(buoy_transfer.apply(tr, b["hs"], b["tp"], b["dir"]), index=b.index).dropna()
        wt = 1.0 / tr["cv_rmse"] ** 2
        converted.append((st, tr["cv_rmse"]))
        conv = h.mul(wt).add(conv if isinstance(conv, pd.Series) else 0.0, fill_value=0.0)
        w_sum = pd.Series(wt, index=h.index).add(w_sum if isinstance(w_sum, pd.Series) else 0.0, fill_value=0.0)
        for col in ("tp", "dir"):                            # fill gaps in the reference's period/direction
            raw[col] = raw[col].fillna(b[col])
    if raw is None:
        sys.exit("no buoy archive found")
    b = raw
    b["hs_marconi"] = (conv / w_sum).reindex(b.index) if converted else np.nan
    # combined error of the converted buoys (independent errors assumed)
    sb = float(1.0 / np.sqrt(sum(1.0 / e ** 2 for _, e in converted))) if converted else np.nan
    tr = bool(converted)
    cam = camera_hourly(args.archive, args.models)
    out = b.join(cam, how="outer").sort_index()
    # wave period measured by the camera's timestacks (timestack_wave_period.py)
    tpc = Path(args.archive) / "wave_period_c2.csv"
    out["tp_camera"] = np.nan
    if tpc.exists():
        p = pd.read_csv(tpc)
        p = p[p["status"] == "ok"]
        if len(p):
            p["hour"] = pd.to_datetime((p["epoch"] // 3600) * 3600, unit="s", utc=True)
            # the MEAN period Tm01: it tracks the buoys (r 0.74 with 44008's APD); the
            # spectral peak from few surf-zone pixels does not
            out = out.join(p.groupby("hour")["tm01_s"].median().rename("tp_cam"), how="outer")
            out["tp_camera"] = out.pop("tp_cam")
    # camera error at its own reading, reduced for n frames sharing most of their error
    n = np.clip(out["n_frames"].fillna(1).to_numpy(float), 1, None)
    sc = camera_sigma(out["hs_camera"].fillna(1.0), args.camera_sigma) * np.sqrt((1 + args.frame_corr * (n - 1)) / n)
    wc = sb ** 2 / (sc ** 2 + sb ** 2)
    sig = np.where(out["hs_camera"].notna() & out["hs_marconi"].notna(), np.sqrt(sc ** 2 * sb ** 2 / (sc ** 2 + sb ** 2)),
                   np.where(out["hs_camera"].notna(), sc, sb))
    best = np.where(out["hs_camera"].notna() & out["hs_marconi"].notna(),
                    wc * out["hs_camera"] + (1 - wc) * out["hs_marconi"],
                    np.where(out["hs_camera"].notna(), out["hs_camera"], out["hs_marconi"]))
    src = np.where(out["hs_camera"].notna() & out["hs_marconi"].notna(), "blend",
                   np.where(out["hs_camera"].notna(), "camera",
                            np.where(out["hs_marconi"].notna(), "buoy", "")))
    # no camera and no conversion (night, before buoy_transfer.py fit): the raw
    # buoy, labelled as such, so downstream users never see a gap
    raw = (src == "") & out["hs"].notna().to_numpy()
    best = np.where(raw, out["hs"], best)
    src = np.where(raw, "buoy (raw)", src)
    sig = np.where(raw, np.nan, sig)                       # unconverted buoy: error unknown
    res = pd.DataFrame({
        "time_utc": out.index.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "epoch": ((out.index - pd.Timestamp("1970-01-01", tz="UTC")) // pd.Timedelta("1s")).astype(int),
        "hs_buoy": out["hs"].round(2), "hs_buoy_marconi": out["hs_marconi"].round(2),
        "hs_camera": out["hs_camera"].round(2), "n_frames": out["n_frames"].fillna(0).astype(int),
        "hs_best": np.round(best.astype(float), 2), "hs_sigma": np.round(sig.astype(float), 2),
        "best_from": src, "n_models": out["n_models"].fillna(0).astype(int),
        "tp_s": out["tp"].round(1), "tp_camera_s": out["tp_camera"].round(1),
        "dir_deg": out["dir"].round(0)})
    keep = (res["best_from"] != "").to_numpy()
    res = res[keep].reset_index(drop=True)
    flag, note = qc_flags(res, np.asarray(sc, float)[keep], sb)
    res["qc_flag"], res["qc_note"] = flag, note
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(args.output, index=False)
    last = res.iloc[-1] if len(res) else None
    print(f"Marconi waves      : {len(res)} hours -> {args.output}"
          + (f"  (converted buoys: {', '.join(st for st, _ in converted)})" if tr
             else "  (no buoy transfer fitted yet: run buoy_transfer.py fit --station 44013)")
          + (f"; latest {last['time_utc']} Hs {last['hs_best']:.2f} m ({last['best_from']})"
             if last is not None else ""))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--archive", default=str(HERE / "archive"))
    ap.add_argument("--buoys", nargs="+", default=["44008=waves_44008.csv", "44013=waves_44013.csv"],
                    help="STATION=ARCHIVE in --archive. The first is the raw reference shown; every one "
                         "with calibration/buoy_transfer_<station>.json is converted to Marconi and "
                         "averaged (44008 had no data in the ADCP winter, so 44013 is fitted)")
    ap.add_argument("--calibration", default=str(HERE / "calibration"))
    ap.add_argument("--models", nargs="+", default=["Run C=owg_c2_H.csv", "v2=owg_c2_H_v2.csv"],
                    help="LABEL=CSV camera archives averaged frame by frame (default Run C and v2)")
    ap.add_argument("--camera-sigma", default="0.3:0.16,0.8:0.23,1.2:0.29,1.7:0.65",
                    help="per-frame camera error by camera Hs, 'HS:SIGMA,...' in m, linear between "
                         "(default: validation RMSE by reading)")
    ap.add_argument("--frame-corr", type=float, default=0.8,
                    help="error correlation of frames in the same hour (validation 0.80)")
    ap.add_argument("--output", default=str(HERE / "archive" / "waves_marconi.csv"))
    return build(ap.parse_args())


if __name__ == "__main__":
    sys.exit(main())

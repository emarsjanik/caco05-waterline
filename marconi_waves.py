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

Period and direction come from the buoy (a single camera frame cannot
measure period -- Chris Sherwood's model scored R2 -0.19 at it).

OUTPUT archive/waves_marconi.csv, one row per hour of buoy record:
  time_utc, epoch, hs_buoy, hs_buoy_marconi, hs_camera, n_frames,
  hs_best, best_from (blend / camera / buoy), tp_s, dir_deg

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


def camera_hourly(archive, specs):
    """Frame-by-frame mean across models, then the mean of the frames in each hour."""
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
    e["hs"] = e[[c for c in e.columns if c not in ("filename", "epoch")]].mean(axis=1)
    e["hour"] = pd.to_datetime((e["epoch"] / 3600).round() * 3600, unit="s", utc=True)
    g = e.groupby("hour")["hs"]
    return pd.DataFrame({"hs_camera": g.mean(), "n_frames": g.size()})


def build(args):
    w = pd.read_csv(args.waves_csv)
    t = pd.to_datetime(w["epoch"], unit="s", utc=True)
    b = pd.DataFrame({"hs": pd.to_numeric(w["wvht_m"], errors="coerce").to_numpy(),
                      "tp": pd.to_numeric(w["dpd_s"], errors="coerce").to_numpy(),
                      "dir": pd.to_numeric(w["mwd_deg"], errors="coerce").to_numpy()}, index=t)
    b = buoy_transfer.hourly(b).dropna(subset=["hs"])
    b[["tp", "dir"]] = b[["tp", "dir"]].ffill(limit=3)      # NDBC often blanks them for an hour

    tr = buoy_transfer.load(args.transfer)
    b["hs_marconi"] = buoy_transfer.apply(tr, b["hs"], b["tp"], b["dir"]) if tr else np.nan
    cam = camera_hourly(args.archive, args.models)
    out = b.join(cam, how="outer").sort_index()

    sc = args.camera_rmse / np.sqrt(np.clip(out["n_frames"].fillna(1), 1, 2))   # 2 frames/h, partly correlated
    sb = tr["cv_rmse"] if tr else np.nan
    wc = sb ** 2 / (sc ** 2 + sb ** 2)
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
    res = pd.DataFrame({
        "time_utc": out.index.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "epoch": ((out.index - pd.Timestamp("1970-01-01", tz="UTC")) // pd.Timedelta("1s")).astype(int),
        "hs_buoy": out["hs"].round(2), "hs_buoy_marconi": out["hs_marconi"].round(2),
        "hs_camera": out["hs_camera"].round(2), "n_frames": out["n_frames"].fillna(0).astype(int),
        "hs_best": np.round(best.astype(float), 2), "best_from": src,
        "tp_s": out["tp"].round(1), "dir_deg": out["dir"].round(0)})
    res = res[res["best_from"] != ""]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(args.output, index=False)
    last = res.iloc[-1] if len(res) else None
    print(f"Marconi waves      : {len(res)} hours -> {args.output}"
          + ("" if tr else "  (no buoy transfer fitted yet: run buoy_transfer.py fit)")
          + (f"; latest {last['time_utc']} Hs {last['hs_best']:.2f} m ({last['best_from']})"
             if last is not None else ""))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--archive", default=str(HERE / "archive"))
    ap.add_argument("--waves-csv", default=str(HERE / "archive" / "waves_44008.csv"))
    ap.add_argument("--transfer", default=str(HERE / "calibration" / "buoy_transfer_44008.json"))
    ap.add_argument("--models", nargs="+", default=["Run C=owg_c2_H.csv", "v2=owg_c2_H_v2.csv"],
                    help="LABEL=CSV camera archives averaged frame by frame (default Run C and v2)")
    ap.add_argument("--camera-rmse", type=float, default=0.27,
                    help="per-frame camera error, m (Run C validation, default 0.27)")
    ap.add_argument("--output", default=str(HERE / "archive" / "waves_marconi.csv"))
    return build(ap.parse_args())


if __name__ == "__main__":
    sys.exit(main())

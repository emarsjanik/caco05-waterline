#!/usr/bin/env python3
"""
Live Optical Wave Gauge: Wave Height From Today's Camera
===========================================================
Reads each new c2 "bright" image the station writes, predicts the
nearshore wave height with the trained network, appends it to a CSV that
only grows, and redraws a 7-day plot against the offshore buoy.

MODEL. Trained in Colab (train_marconi_owg.py) on Dec 2024 - Mar 2025
frames redrawn into today's camera view (resample_owg_images.py), with
the ADCP's wave height as truth; exported by export_owg_onnx.py. Run C
(Inception-ResNetV2, crop 0.078,0.50,0.30,1.0): 0.271 m RMSE on 344
held-out frames, bias +0.01 m. Runs with OpenCV alone -- the station has
no TensorFlow. --check first: it proves the exported file gives the
same numbers as TensorFlow did.

SAME PREPARATION AS TRAINING. Crop, grey conversion, resize and
per-image scaling are imported from train_marconi_owg.py itself, and
the crop and input size are read from the model's report, so they
cannot drift apart.

SAME FRAMES AS TRAINING. The network only learned from frames that
passed the quality screen (score_image_quality.py): brightness >= 45
(not night), sharpness >= 10 (not fog or a wet lens), and less than
55% of the sea crop blown out by sun glint. Other frames are logged
with the reason and not given a height: a confident number from a
frame unlike any it trained on is worse than none.

AVERAGING. Per-frame errors are partly random, so the plot also shows a
5-frame (2-hour) centred average: on the held-out frames that cut the
error 11% (average_predictions.py, 0.219 -> 0.195 m on the same frames).

BUOY. 44008 is offshore (Nantucket Shoals); its height is an index of
the incoming sea, not the height at the beach. Expect the two to rise
and fall together, with the camera's nearshore value lower in big seas
(breaking) and different in direction-sheltered conditions.

Usage (on the station):
    python3 owg_live.py --check
    python3 owg_live.py                       # process new images, redraw plot
    python3 owg_live.py --reprocess           # recompute every image still on disk
"""

import re
import sys
import json
import argparse
from pathlib import Path
from datetime import datetime, timezone, timedelta

import numpy as np
import pandas as pd
import cv2

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from train_marconi_owg import normalise, prepare_image   # noqa: E402
from score_image_quality import score                        # noqa: E402

FIELDS = ["filename", "epoch", "time_utc", "brightness", "sharpness", "saturation",
          "status", "hs_m", "water_level_m", "water_level_source"]
EXPECTED_SIZE = (2448, 2048)


def load_net(stem):
    rep = json.loads(Path(stem + ".report.json").read_text())
    net = cv2.dnn.readNetFromONNX(stem + ".onnx")
    width = int(rep["img_size"])
    height = int(rep.get("img_height") or width)
    crop = tuple(rep["crop"]) if rep.get("crop") else None
    return net, rep, width, height, crop


RGB_MEAN = np.array([0.485, 0.456, 0.406], np.float32)
RGB_STD = np.array([0.229, 0.224, 0.225], np.float32)


def is_rgb(rep):
    """train_owg_torch.py models take colour, ImageNet-normalised, channels first."""
    return rep.get("input") == "rgb_imagenet"


def model_input(img_bgr, rep, width, height, crop):
    """Image (BGR, full frame or sea patch) -> the array the model expects."""
    if not is_rgb(rep):
        g = prepare_image(img_bgr, width, height, crop)
        return None if g is None else normalise(g)
    if crop:
        h, w = img_bgr.shape[:2]
        t, b, l, r = crop
        img_bgr = img_bgr[int(t * h):int(b * h), int(l * w):int(r * w)]
    rgb = cv2.cvtColor(cv2.resize(img_bgr, (width, height), interpolation=cv2.INTER_AREA),
                       cv2.COLOR_BGR2RGB)
    return ((rgb.astype(np.float32) / 255.0 - RGB_MEAN) / RGB_STD).transpose(2, 0, 1)


def run(net, x):
    """x: model input from model_input() (H, W grey or 3, H, W colour) -> wave height (m)."""
    x = x[None, :, :, None] if x.ndim == 2 else x[None]
    net.setInput(np.ascontiguousarray(x, dtype=np.float32))
    return float(net.forward().ravel()[0])


def check(stem, net):
    chk = json.loads(Path(stem + ".check.json").read_text())
    X = np.load(Path(stem).parent / chk["inputs"])
    got = np.array([run(net, X[i] if X.ndim == 4 and X.shape[1] == 3 else X[i, :, :, 0])
                    for i in range(len(X))])
    exp = np.array(chk["expected"])
    diff = float(np.abs(got - exp).max())
    print(f"exported model vs {chk.get('framework', 'TensorFlow')}: largest difference {diff * 1000:.3f} mm "
          f"over {len(X)} test inputs")
    for g, e in zip(got, exp):
        print(f"    {g:+.4f}  expected {e:+.4f}")
    if diff > 1e-3:
        print("FAILED: the exported file does not reproduce the trained model.")
        return 1
    print("OK: the station reproduces the trained model.")
    return 0


def find_images(dirs, camera):
    seen = {}
    for d in dirs:
        for p in Path(d).glob(f"*.{camera}.bright.jpg"):
            seen.setdefault(p.name, p)
    return seen


def averaged(ok, window, max_gap_min):
    """Centred `window`-frame averages, never across a gap; -> (times, values) split at gaps."""
    ep = ok["epoch"].to_numpy(); h = ok["hs_m"].to_numpy(float)
    half = window // 2
    avg_t, avg_v = [], []
    for i in range(half, len(ok) - half):
        if np.diff(ep[i - half:i + half + 1]).max(initial=0) / 60 > max_gap_min:
            continue
        avg_t.append(ok["t"].iloc[i]); avg_v.append(h[i - half:i + half + 1].mean())
    if not avg_v:
        return []
    at = pd.DatetimeIndex(avg_t).tz_convert(None).to_numpy(); av = np.array(avg_v)
    gaps = np.where(np.diff(at).astype("timedelta64[m]").astype(float) > max_gap_min)[0] + 1
    return list(zip(np.split(at, gaps), np.split(av, gaps)))


def plot(df, waves_csv, out, days, window, max_gap_min, model_name, rmse=None, others=()):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    ok = df[df["status"] == "ok"].copy()
    end = datetime.now(timezone.utc)
    start = end - timedelta(days=days)
    ok["t"] = pd.to_datetime(ok["epoch"], unit="s", utc=True)
    ok = ok[ok["t"] >= start].sort_values("t")

    BLUE, ORANGE = "#2a78d6", "#eb6834"
    INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dc"
    fig, ax = plt.subplots(figsize=(12, 4.6), dpi=110)
    fig.patch.set_facecolor("white")

    if waves_csv and Path(waves_csv).exists():
        w = pd.read_csv(waves_csv)
        w = w[pd.to_numeric(w["wvht_m"], errors="coerce").notna()]
        w["t"] = pd.to_datetime(w["epoch"], unit="s", utc=True)
        w = w[w["t"] >= start].sort_values("t")
        if len(w):
            # break the line across gaps of more than 3 hours
            t = w["t"].to_numpy(); v = w["wvht_m"].astype(float).to_numpy()
            gaps = np.where(np.diff(t).astype("timedelta64[m]").astype(float) > 180)[0] + 1
            for i, (ts, vs) in enumerate(zip(np.split(t, gaps), np.split(v, gaps))):
                ax.plot(ts, vs, color=ORANGE, lw=2,
                        label="Buoy 44008 (offshore)" if i == 0 else None)

    n_avg = 0
    if len(ok):
        ax.scatter(ok["t"], ok["hs_m"], s=14, color=BLUE, alpha=0.35, lw=0,
                   label="Camera, single frames")
        # centred average of `window` frames, never across a gap (night, rejected frames)
        pieces = averaged(ok, window, max_gap_min)
        n_avg = sum(len(v) for _, v in pieces)
        for i, (ts, vs) in enumerate(pieces):
            ax.plot(ts, vs, color=BLUE, lw=2,
                    label=f"Camera ({model_name}), {window}-frame average" if i == 0 else None)

    # other models run alongside, for comparison (owg_live.py --also)
    for label, path in others:
        if not Path(path).exists():
            continue
        o = pd.read_csv(path)
        o = o[o["status"] == "ok"].copy()
        o["t"] = pd.to_datetime(o["epoch"], unit="s", utc=True)
        o = o[o["t"] >= start].sort_values("t")
        for i, (ts, vs) in enumerate(averaged(o, window, max_gap_min)):
            ax.plot(ts, vs, color="#1baf7a", lw=2,
                    label=f"Camera ({label}), {window}-frame average" if i == 0 else None)

    ax.set_xlim(start, end)
    if not len(ok) or ok["hs_m"].min() >= 0:
        ax.set_ylim(bottom=0)
    ax.set_ylabel("Significant wave height (m)", color=INK)
    n_rej = int(((df["status"] != "ok") & (df["epoch"] >= start.timestamp())).sum())
    ax.set_title(f"Optical wave gauge, CACO05 c2 -- last {days} days (UTC)\n"
                 f"{len(ok)} frames measured, {n_rej} screened out (night, fog, glare); "
                 f"model {model_name}" + (f", validation RMSE {rmse:.2f} m per frame" if rmse else ""),
                 loc="left", fontsize=10, color=INK)
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED)
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print(f"plot              : {out} ({len(ok)} frames, {n_avg} averaged points)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--model", default=str(HERE / "owg_models" / "owg_c2_H_current_C"),
                    help="Model stem: <stem>.onnx and <stem>.report.json (export_owg_onnx.py)")
    ap.add_argument("--image-dir", nargs="+",
                    default=["/mnt/I2Rgus_Data/ImageProducts/products",
                             "/mnt/I2Rgus_Data/ImageProducts"])
    ap.add_argument("--camera", default="c2")
    ap.add_argument("--output", default=str(HERE / "archive" / "owg_c2_H.csv"))
    ap.add_argument("--plot", default=str(HERE / "owg_c2_H_7day.png"))
    ap.add_argument("--waves-csv", default=str(HERE / "archive" / "waves_44008.csv"))
    ap.add_argument("--days", type=int, default=7)
    ap.add_argument("--window", type=int, default=5,
                    help="Frames in the plotted average (default 5 = 2 h; odd)")
    ap.add_argument("--max-gap", type=float, default=45.0,
                    help="Minutes; frames further apart are not averaged together")
    ap.add_argument("--min-brightness", type=float, default=45.0)
    ap.add_argument("--min-sharpness", type=float, default=10.0)
    ap.add_argument("--max-saturation", type=float, default=0.55)
    ap.add_argument("--check", action="store_true",
                    help="Compare the exported model with TensorFlow's answers, then stop")
    ap.add_argument("--reprocess", action="store_true",
                    help="Recompute images already in the CSV (after a model change)")
    ap.add_argument("--also", action="append", default=[],
                    help="LABEL=CSV of another model's archive to draw in the plot, e.g. "
                         "'patch=archive/owg_c2_H_patch.csv' (repeatable)")
    ap.add_argument("--gauge-csv", default=str(HERE / "archive" / "gauge_8447435.csv"),
                    help="Chatham gauge archive (fetch_tide_gauge.py), for the water level a "
                         "patch model projects at")
    ap.add_argument("--gnssr-spline",
                    default="/home/argus_user/GNSS/v4.1/products/refl_code/Files/usgs/usgs_spline_out.txt",
                    help="the station's GNSS-R water level (gnssrefl spline): used first")
    ap.add_argument("--no-gps", action="store_true",
                    help="Chatham gauge only (the behaviour before the GPS was used)")
    ap.add_argument("--gauge-ahead", type=float, default=2.0,
                    help="hours past the last gauge reading to carry a tide fit forward "
                         "(default 2; 0 = wait for the gauge)")
    ap.add_argument("--screen-crop", default="0.078,0.50,0.30,1.0",
                    help="Region of the live frame the quality screen scores when the model "
                         "uses a sea patch (default: the surf zone in today's c2 view)")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    for ext in (".onnx", ".report.json"):
        if not Path(args.model + ext).exists():
            sys.exit(f"missing {args.model}{ext} -- copy it from Colab (export_owg_onnx.py)")
    net, rep, width, height, crop = load_net(args.model)
    if args.check:
        return check(args.model, net)

    # A model trained on a fixed patch of sea (rectify_owg_images.py) comes with
    # <model>.patch.json: live frames are projected onto that patch with today's
    # calibration, at the water level from the Chatham gauge.
    patch = None
    if Path(args.model + ".patch.json").exists():
        from sea_patch import Patch, rectify
        from georectify import load_extrinsics, load_intrinsics
        patch = Patch.from_json(args.model + ".patch.json")
        p_io = load_intrinsics(HERE / "calibration" / f"CACO05_{args.camera}_20240801_IO.yaml")
        p_eo = load_extrinsics(HERE / "calibration" / f"CACO05_{args.camera}_20251113_EO-CV.yaml")
        # water level: the station's own GPS (GNSS-R, on the camera tower) where it
        # covers the frame, else Chatham converted with a fit to the GPS, else a
        # tide fit up to --gauge-ahead hours past the gauge (marconi_water_level.py)
        from marconi_water_level import WaterLevel
        wl = WaterLevel(None if args.no_gps else args.gnssr_spline, args.gauge_csv,
                        ahead_hours=args.gauge_ahead)
        screen = tuple(float(v) for v in args.screen_crop.split(","))
        print(f"sea patch         : {patch.shape[1]} x {patch.shape[0]} px at {patch.res} m")
        print(f"                    {wl.describe()}")

        def water_level(epoch):
            lv, src = wl.at([epoch])
            return (None, "") if not np.isfinite(lv[0]) else (float(lv[0]), str(src[0]))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(out) if out.exists() else pd.DataFrame(columns=FIELDS)
    images = find_images(args.image_dir, args.camera)
    # frames that waited for a water level are retried: the gauge record catches up
    retry = set(df.loc[df["status"] == "no water level", "filename"]) if len(df) else set()
    todo = sorted(n for n in images if args.reprocess or n in retry or n not in set(df["filename"]))
    print(f"images on disk    : {len(images)} {args.camera} bright, {len(todo)} new")

    rows, warned = [], False
    for name in todo:
        m = re.match(r"^(\d{9,11})\.", name)
        if not m:
            continue
        epoch = int(m.group(1))
        path = images[name]
        row = {"filename": name, "epoch": epoch,
               "time_utc": datetime.fromtimestamp(epoch, tz=timezone.utc).isoformat(),
               "hs_m": np.nan}
        b, sh, sa = score(path, crop=screen if patch is not None else crop)
        row.update(brightness=round(b, 2), sharpness=round(sh, 2), saturation=round(sa, 4))
        if not np.isfinite(b):
            row["status"] = "unreadable"
        elif b < args.min_brightness:
            row["status"] = "dark"
        elif sa > args.max_saturation:
            row["status"] = "glare"
        elif sh < args.min_sharpness:
            row["status"] = "blurred (fog, rain or wet lens)"
        else:
            if not warned:
                hw = cv2.imread(str(path), cv2.IMREAD_REDUCED_GRAYSCALE_8).shape
                if (hw[1] * 8, hw[0] * 8) != EXPECTED_SIZE:
                    print(f"WARNING: frames are about {hw[1] * 8}x{hw[0] * 8}, the model "
                          f"trained on {EXPECTED_SIZE[0]}x{EXPECTED_SIZE[1]}")
                warned = True
            if patch is None:
                full = cv2.imread(str(path))
                img = None if full is None else model_input(full, rep, width, height, crop)
            else:
                img = None
                z, z_src = water_level(epoch)
                full = cv2.imread(str(path))
                if z is None:
                    row["status"] = "no water level"
                elif full is not None:
                    rect, cov = rectify(full, p_io, p_eo, z, patch)
                    row["water_level_m"] = round(z, 3)
                    row["water_level_source"] = z_src
                    if cov < 0.98:
                        row["status"] = "patch not in view"
                    else:
                        # through JPEG, as the training images were
                        ok_, enc = cv2.imencode(".jpg", rect, [cv2.IMWRITE_JPEG_QUALITY, 95])
                        img = model_input(cv2.imdecode(enc, cv2.IMREAD_COLOR), rep, width, height, None)
            if row.get("status"):
                pass
            elif img is None:
                row["status"] = "unreadable"
            else:
                row["hs_m"] = round(run(net, img), 3)
                row["status"] = "ok"
        rows.append(row)

    if rows:
        new = pd.DataFrame(rows, columns=FIELDS)
        df = pd.concat([df[~df["filename"].isin(new["filename"])], new], ignore_index=True)
        df = df.sort_values("epoch")
        df.to_csv(out, index=False)
        counts = new["status"].value_counts().to_dict()
        print(f"processed         : {len(new)} -> " +
              ", ".join(f"{v} {k}" for k, v in counts.items()))
        last = new[new["status"] == "ok"].tail(1)
        if len(last):
            print(f"latest            : {last['time_utc'].iloc[0]}  "
                  f"Hs {last['hs_m'].iloc[0]:.2f} m")
        print(f"archive           : {out} ({len(df)} rows)")
    if not args.no_plot and len(df):
        df["epoch"] = df["epoch"].astype(int)
        others = [tuple(a.split("=", 1)) for a in args.also if "=" in a]
        plot(df, args.waves_csv, args.plot, args.days, args.window, args.max_gap,
             Path(args.model).name, rep.get("val_rmse"), others)
    return 0


if __name__ == "__main__":
    sys.exit(main())

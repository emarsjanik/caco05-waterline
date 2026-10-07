#!/usr/bin/env python3
"""
Optical Wave Gauge Daily Report
=================================
Writes the OWG report the station emails every morning: a short text
synopsis (conditions, what the camera measured, agreement with the buoy,
and pass/fail checks on every link of the data chain), plus the plots
that matter for judging the measurement. owg.sh runs it; --email sends
it with msmtp, the same transport as the GPS health email.

PLOTS (written to --out-dir, attached to the email)
  owg_report_7day.png        wave height from both camera models and the
                             buoy, the buoy's period and direction, and
                             the tide -- the context every reading needs
  owg_report_scatter.png     camera vs buoy for each model: the last
                             7 days over the whole record, with n, bias,
                             RMSE, correlation and Willmott's d
  owg_report_residuals.png   camera minus buoy against the tide, the
                             wave direction and the wave period: a
                             model that depends on any of them is
                             reading the beach, not the waves
  owg_report_frames.png      frames per day by outcome: measured, or
                             screened out (dark, glare, fog, no water
                             level) -- the health of the image supply
  owg_report_view.png        the latest frame with what each model
                             looks at: Run C's crop and the sea patch

THE BUOY IS NOT TRUTH. 44008 is ~50 nm offshore; the camera sees the
nearshore. Expect them to rise and fall together (correlation), with
the camera lower in big seas (breaking) and when the waves come from a
direction Marconi is sheltered from. The bias and RMSE against the buoy
are a consistency check, not the camera's accuracy -- that was measured
against the ADCP in validation (Run C 0.27 m per frame).

CHECKS. Each prints [ OK ] or [FAIL]; the exit status is the number of
failures (0 = all healthy), which owg.sh puts in the email subject.

Usage:
    python3 owg_report.py                       write the report and plots
    python3 owg_report.py --email a@x b@y       ... and send them
    python3 owg_report.py --email a@x --dry-run show what would be sent
"""

import sys
import json
import time
import argparse
import subprocess
from pathlib import Path
from datetime import datetime, timezone, timedelta

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dc"
ORANGE, BLUE, GREEN, GREY = "#eb6834", "#2a78d6", "#1baf7a", "#b9b8b0"
PALETTE = (BLUE, GREEN, "#8a5cd1", "#c8417a", "#3b3a36")
# r and Willmott d are only reported when the buoy spans at least this much (m)
MIN_SPREAD = 1.0
COMPASS = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
           "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]


# ---------------------------------------------------------------- data

def load_archive(path):
    if not Path(path).exists():
        return None
    d = pd.read_csv(path)
    if not len(d):
        return d
    d["t"] = pd.to_datetime(d["epoch"], unit="s", utc=True)
    d["hs_m"] = pd.to_numeric(d["hs_m"], errors="coerce")
    return d.sort_values("epoch").reset_index(drop=True)


def load_waves(path):
    if not Path(path).exists():
        return None
    w = pd.read_csv(path)
    for c in ("wvht_m", "dpd_s", "apd_s", "mwd_deg"):
        if c in w:
            w[c] = pd.to_numeric(w[c], errors="coerce")
    w["t"] = pd.to_datetime(w["epoch"], unit="s", utc=True)
    return w.sort_values("epoch").reset_index(drop=True)


def load_gauge(path):
    if not Path(path).exists():
        return None
    g = pd.read_csv(path)
    g["level_navd88"] = pd.to_numeric(g["level_navd88"], errors="coerce")
    g = g.dropna(subset=["level_navd88"]).sort_values("epoch").reset_index(drop=True)
    g["t"] = pd.to_datetime(g["epoch"], unit="s", utc=True)
    return g


def interp_at(src_ep, src_v, ep, max_gap_s):
    """Linear interpolation, NaN where the bracketing samples are more than max_gap_s apart."""
    src_ep = np.asarray(src_ep, float); src_v = np.asarray(src_v, float)
    ok = np.isfinite(src_v)
    src_ep, src_v = src_ep[ok], src_v[ok]
    ep = np.asarray(ep, float)
    if len(src_ep) < 2:
        return np.full(len(ep), np.nan)
    out = np.interp(ep, src_ep, src_v)
    i = np.clip(np.searchsorted(src_ep, ep), 1, len(src_ep) - 1)
    bad = (ep < src_ep[0]) | (ep > src_ep[-1]) | (src_ep[i] - src_ep[i - 1] > max_gap_s)
    out[bad] = np.nan
    return out


def averaged(ok, window=5, max_gap_min=45.0):
    """Centred `window`-frame averages, never across a gap (as owg_live.py plots them)."""
    ep = ok["epoch"].to_numpy(float); h = ok["hs_m"].to_numpy(float)
    half = window // 2
    t, v = [], []
    for i in range(half, len(ok) - half):
        if np.diff(ep[i - half:i + half + 1]).max(initial=0) / 60 > max_gap_min:
            continue
        t.append(ep[i]); v.append(h[i - half:i + half + 1].mean())
    return np.array(t), np.array(v)


def metrics(cam, ref):
    m = np.isfinite(cam) & np.isfinite(ref)
    cam, ref = cam[m], ref[m]
    n = len(cam)
    if n < 3:
        return {"n": n}
    e = cam - ref
    d_den = np.sum((np.abs(cam - ref.mean()) + np.abs(ref - ref.mean())) ** 2)
    return {"n": n, "bias": float(e.mean()), "rmse": float(np.sqrt(np.mean(e ** 2))),
            "corr": float(np.corrcoef(cam, ref)[0, 1]) if cam.std() > 0 and ref.std() > 0 else np.nan,
            "d": float(1 - np.sum(e ** 2) / d_den) if d_den > 0 else np.nan}


def compass(deg):
    return "" if not np.isfinite(deg) else COMPASS[int((deg % 360) / 22.5 + 0.5) % 16]


def ago(seconds):
    if seconds < 3600:
        return f"{seconds / 60:.0f} min"
    if seconds < 2 * 86400:
        return f"{seconds / 3600:.1f} h"
    return f"{seconds / 86400:.1f} days"


# ---------------------------------------------------------------- plots

def style(ax):
    ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)


def split_gaps(t, v, gap_s):
    t = np.asarray(t, float); v = np.asarray(v, float)
    if not len(t):
        return []
    cut = np.where(np.diff(t) > gap_s)[0] + 1
    return list(zip(np.split(t, cut), np.split(v, cut)))


def to_dt(ep):
    return pd.to_datetime(np.asarray(ep, float), unit="s").to_numpy()


def plot_7day(models, waves, gauge, fit, start, end, out, mw=None):
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    fig, ax = plt.subplots(3, 1, figsize=(12, 9), dpi=110, sharex=True,
                           gridspec_kw={"height_ratios": [2.2, 1, 1]})
    fig.patch.set_facecolor("white")
    s0, s1 = start.timestamp(), end.timestamp()

    a = ax[0]
    if waves is not None:
        w = waves[(waves["epoch"] >= s0) & waves["wvht_m"].notna()]
        for i, (t, v) in enumerate(split_gaps(w["epoch"], w["wvht_m"], 3 * 3600)):
            a.plot(to_dt(t), v, color=ORANGE, lw=2, label="Buoy 44008 (offshore)" if i == 0 else None)
    if mw is not None:
        m = mw[mw["epoch"] >= s0]
        c = m[m["hs_buoy_marconi"].notna()]
        for i, (t, v) in enumerate(split_gaps(c["epoch"], c["hs_buoy_marconi"], 3 * 3600)):
            a.plot(to_dt(t), v, color=ORANGE, lw=1.5, ls="--",
                   label="Buoy converted to Marconi" if i == 0 else None)
        bst = m[m["hs_best"].notna()]
        for i, (t, v) in enumerate(split_gaps(bst["epoch"], bst["hs_best"], 3 * 3600)):
            a.plot(to_dt(t), v, color=INK, lw=2.4,
                   label="Best estimate (camera + converted buoy)" if i == 0 else None)
    for (label, d), colour in zip(models, PALETTE):
        if d is None:
            continue
        ok = d[(d["status"] == "ok") & (d["epoch"] >= s0)]
        if not len(ok):
            continue
        a.scatter(to_dt(ok["epoch"]), ok["hs_m"], s=10, color=colour, alpha=0.3, lw=0)
        t, v = averaged(ok)
        for i, (tt, vv) in enumerate(split_gaps(t, v, 45 * 60)):
            a.plot(to_dt(tt), vv, color=colour, lw=2,
                   label=f"Camera {label} (2-h mean; dots: frames)" if i == 0 else None)
    a.set_ylim(bottom=0)
    a.set_ylabel("Significant wave\nheight (m)", color=INK)
    fig.suptitle(f"Optical wave gauge, CACO05 c2 -- last {round((end - start).total_seconds() / 86400)} "
                 f"days (UTC)", x=0.01, ha="left", fontsize=11, color=INK)
    h, lab = a.get_legend_handles_labels()
    fig.legend(h, lab, loc="upper left", bbox_to_anchor=(0.01, 0.965), ncol=3, frameon=False, fontsize=8.5)

    a = ax[1]
    if waves is not None:
        w = waves[waves["epoch"] >= s0]
        a.plot(to_dt(w["epoch"]), w["dpd_s"], ".", color=ORANGE, ms=4, label="dominant period (s)")
        a.set_ylabel("Period (s)", color=INK)
        b = a.twinx()
        b.plot(to_dt(w["epoch"]), w["mwd_deg"], ".", color=MUTED, ms=4, label="mean direction (from)")
        b.set_ylim(0, 360); b.set_yticks([0, 90, 180, 270, 360])
        b.set_yticklabels(["N", "E", "S", "W", "N"])
        b.tick_params(colors=MUTED, labelsize=9)
        for s in ("top",):
            b.spines[s].set_visible(False)
        b.set_ylabel("Direction from", color=MUTED)
        h1, l1 = a.get_legend_handles_labels(); h2, l2 = b.get_legend_handles_labels()
        a.legend(h1 + h2, l1 + l2, loc="upper left", frameon=False, fontsize=8, ncol=2)

    a = ax[2]
    t = np.arange(s0, s1, 600.0)
    lv, src = fit.at(t)                     # fit: a marconi_water_level.WaterLevel
    gps = np.where(src == "gps", lv, np.nan)
    other = np.where(src != "gps", lv, np.nan)
    a.plot(to_dt(t), gps, color=BLUE, lw=1.8, label="Marconi still-water level, GPS (GNSS-R on the camera tower)")
    a.plot(to_dt(t), other, color=BLUE, lw=1.2, ls="--",
           label="newest hours: Chatham gauge, converted with a fit to the GPS")
    a.set_ylabel("Water level\n(m NAVD88)", color=INK)
    a.legend(loc="upper left", frameon=False, fontsize=8)

    for x in ax:
        style(x)
    ax[-1].set_xlim(to_dt([s0])[0], to_dt([s1])[0])
    ax[-1].xaxis.set_major_locator(mdates.DayLocator())
    ax[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(out)
    plt.close(fig)


def plot_scatter(matched, start, out):
    import matplotlib.pyplot as plt
    labels = [l for l, m in matched if m is not None and len(m)]
    if not labels:
        return False
    fig, axs = plt.subplots(1, len(labels), figsize=(5.6 * len(labels), 5.4), dpi=110, squeeze=False)
    fig.patch.set_facecolor("white")
    hi = max(max(np.nanmax(m["buoy"]), np.nanmax(m["hs_m"])) for l, m in matched if m is not None and len(m))
    hi = max(1.0, np.ceil(hi * 2) / 2)
    for a, (label, m) in zip(axs[0], [(l, m) for l, m in matched if m is not None and len(m)]):
        recent = m["epoch"] >= start.timestamp()
        a.scatter(m.loc[~recent, "buoy"], m.loc[~recent, "hs_m"], s=10, color=GREY, lw=0,
                  label="earlier")
        a.scatter(m.loc[recent, "buoy"], m.loc[recent, "hs_m"], s=14, color=BLUE, lw=0, alpha=0.7,
                  label="last 7 days")
        a.plot([0, hi], [0, hi], color=INK, lw=1)
        a.set_xlim(0, hi); a.set_ylim(0, hi); a.set_aspect("equal")
        s7, sa = metrics(m.loc[recent, "hs_m"].to_numpy(float), m.loc[recent, "buoy"].to_numpy(float)), \
            metrics(m["hs_m"].to_numpy(float), m["buoy"].to_numpy(float))

        def fmt(s, name):
            if s["n"] < 3:
                return f"{name}: n={s['n']}"
            return (f"{name}: n={s['n']}  bias {s['bias']:+.2f}  RMSE {s['rmse']:.2f} m\n"
                    f"      r {s['corr']:.2f}  d {s['d']:.2f}")
        a.text(0.03, 0.97, fmt(s7, "7 days") + "\n" + fmt(sa, "all   "), transform=a.transAxes,
               va="top", fontsize=8.5, family="monospace", color=INK)
        a.set_title(f"Camera {label} vs buoy 44008 (single frames)", loc="left", fontsize=10, color=INK)
        a.set_xlabel("Buoy Hs (m), offshore", color=INK); a.set_ylabel("Camera Hs (m)", color=INK)
        a.legend(loc="lower right", frameon=False, fontsize=8)
        style(a)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return True


def plot_residuals(matched, out):
    import matplotlib.pyplot as plt
    have = [(l, m) for l, m in matched if m is not None and len(m)]
    if not have:
        return False
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), dpi=110, sharey=True)
    fig.patch.set_facecolor("white")
    for (label, m), colour in zip(have, PALETTE):
        r = m["hs_m"] - m["buoy"]
        for a, col in zip(axs, ("level", "mwd", "dpd")):
            a.scatter(m[col], r, s=10, color=colour, alpha=0.5, lw=0, label=f"Camera {label}")
    for a, (xl, title) in zip(axs, (("Marconi water level (m NAVD88)", "vs tide"),
                                    ("Buoy mean direction, from (deg)", "vs wave direction"),
                                    ("Buoy dominant period (s)", "vs wave period"))):
        a.axhline(0, color=INK, lw=1)
        a.set_xlabel(xl, color=INK); a.set_title(f"Camera minus buoy {title}", loc="left",
                                                   fontsize=10, color=INK)
        style(a)
    axs[1].set_xlim(0, 360); axs[1].set_xticks([0, 45, 90, 135, 180, 225, 270, 315, 360])
    axs[0].set_ylabel("Camera - buoy Hs (m)", color=INK)
    axs[0].legend(loc="best", frameon=False, fontsize=8)
    fig.suptitle("A trend in any panel means the reading depends on more than the waves "
                 "(Marconi faces ~80 deg, ENE-E)", x=0.01, ha="left", fontsize=9, color=MUTED)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return True


def plot_frames(models, start, end, out):
    import matplotlib.pyplot as plt
    have = [(l, d) for l, d in models if d is not None and len(d)]
    if not have:
        return False
    cats = [("ok", "measured", BLUE), ("dark", "dark", "#3b3a36"), ("glare", "glare", "#f2b134"),
            ("blurred (fog, rain or wet lens)", "fog / rain / wet lens", GREY),
            ("low sharpness (may read low)", "view degraded (measured, not used)", "#a8a59a"),
            ("no water level", "no water level", "#c8417a"), ("patch not in view", "patch not in view", "#8a5cd1"),
            ("unreadable", "unreadable", "#d23b2a")]
    days = pd.date_range(start.date(), end.date(), freq="D", tz="UTC")
    fig, axs = plt.subplots(1, len(have), figsize=(6.2 * len(have), 3.8), dpi=110, squeeze=False,
                            sharey=True)
    fig.patch.set_facecolor("white")
    for a, (label, d) in zip(axs[0], have):
        d = d[d["t"] >= days[0]]
        day = d["t"].dt.floor("D")
        bottom = np.zeros(len(days))
        for status, name, colour in cats:
            n = np.array([int(((day == x) & (d["status"] == status)).sum()) for x in days])
            if n.any():
                a.bar(days.tz_convert(None), n, bottom=bottom, color=colour, width=0.8, label=name)
                bottom += n
        a.set_title(f"Camera {label}: frames per day by outcome", loc="left", fontsize=10, color=INK)
        a.set_xticks(days.tz_convert(None)); a.set_xticklabels([x.strftime("%b %d") for x in days],
                                                                rotation=0, fontsize=8)
        a.legend(loc="upper left", frameon=False, fontsize=7.5)
        style(a)
    axs[0][0].set_ylabel("frames (UTC day)", color=INK)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return True


def plot_view(model_c, model_patch, image_dirs, out):
    """Latest c2 bright frame with Run C's crop and the sea patch drawn on it, and both inputs."""
    import cv2
    import matplotlib.pyplot as plt
    frames = []
    for d in image_dirs:
        frames += list(Path(d).glob("*.c2.bright.jpg"))
    if not frames:
        return None
    frames.sort(key=lambda p: p.name)
    path = frames[-1]
    img = cv2.imread(str(path))
    if img is None:
        return None
    h, w = img.shape[:2]
    panels = []
    show = img.copy()
    rep = Path(model_c + ".report.json")
    if rep.exists():
        crop = json.loads(rep.read_text()).get("crop")
        if crop:
            t, b, l, r = crop
            cv2.rectangle(show, (int(l * w), int(t * h)), (int(r * w), int(b * h)), (214, 120, 42), 8)
            panels.append(("Run C input (crop)", img[int(t * h):int(b * h), int(l * w):int(r * w)]))
    pj = Path(model_patch + ".patch.json")
    if pj.exists():
        try:
            sys.path.insert(0, str(HERE))
            from sea_patch import Patch, rectify
            from georectify import load_extrinsics, load_intrinsics
            from view_reproject import ground_to_pixel
            patch = Patch.from_json(pj)
            io = load_intrinsics(HERE / "calibration" / "CACO05_c2_20240801_IO.yaml")
            eo = load_extrinsics(HERE / "calibration" / "CACO05_c2_20251113_EO-CV.yaml")
            z = 0.0
            arch = load_archive(HERE / "archive" / "owg_c2_H_patch.csv")
            if arch is not None and len(arch) and "water_level_m" in arch:
                row = arch[arch["filename"] == path.name]
                if len(row) and np.isfinite(pd.to_numeric(row["water_level_m"], errors="coerce").iloc[0]):
                    z = float(row["water_level_m"].iloc[0])
            c = patch.corners()
            U, V, ok = ground_to_pixel(np.r_[c[:, 0], c[0, 0]], np.r_[c[:, 1], c[0, 1]], z, io, eo)
            if ok.all():
                pts = np.c_[U, V].astype(np.int32).reshape(-1, 1, 2)
                cv2.polylines(show, [pts], False, (122, 175, 27), 8)
            rect, _ = rectify(img, io, eo, z, patch)
            panels.append((f"Sea patch input ({patch.along[1] - patch.along[0]:.0f} x "
                           f"{patch.cross[1] - patch.cross[0]:.0f} m, at {z:+.2f} m)", rect))
        except Exception as exc:                      # a picture is not worth failing the report
            print(f"view: sea patch not drawn ({exc})", file=sys.stderr)
    fig = plt.figure(figsize=(12, 5.4), dpi=100)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(len(panels) or 1, 2, width_ratios=[1.5, 1])
    a = fig.add_subplot(gs[:, 0])
    a.imshow(cv2.cvtColor(show, cv2.COLOR_BGR2RGB)); a.axis("off")
    when = datetime.fromtimestamp(int(path.name.split(".")[0]), tz=timezone.utc)
    a.set_title(f"Latest c2 frame, {when:%Y-%m-%d %H:%M}Z  (blue: Run C crop, green: sea patch)",
                loc="left", fontsize=9, color=INK)
    for i, (title, p) in enumerate(panels):
        b = fig.add_subplot(gs[i, 1])
        b.imshow(cv2.cvtColor(p, cv2.COLOR_BGR2RGB)); b.axis("off")
        b.set_title(title, loc="left", fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return when


# ---------------------------------------------------------------- report

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--archive", default=str(HERE / "archive"))
    ap.add_argument("--models", nargs="+",
                    default=["Run C=owg_c2_H.csv=owg_c2_H_current_C", "patch=owg_c2_H_patch.csv=owg_c2_H_patch",
                             "v2=owg_c2_H_v2.csv=owg_c2_H_v2",
                             "composite=owg_c2_H_composite.csv=owg_c2_H_composite"],
                    help="LABEL=CSV[=MODEL]: archive CSV in --archive and model stem in owg_models/, "
                         "in plotting order. A model with neither an archive nor a model file is "
                         "left out, so a model not installed yet is not a failure.")
    ap.add_argument("--ensemble", default="Run C,v2",
                    help="labels averaged, frame by frame, into an 'ensemble' reading when at least "
                         "two of them measured the frame (default 'Run C,v2'; '' for none)")
    ap.add_argument("--model-c", default=None, help="model stem for 'Run C' (overrides --models)")
    ap.add_argument("--model-patch", default=None, help="model stem for 'patch' (overrides --models)")
    ap.add_argument("--waves-csv", default=None, help="default <archive>/waves_44008.csv")
    ap.add_argument("--gauge-csv", default=None, help="default <archive>/gauge_8447435.csv")
    ap.add_argument("--marconi-waves", default=None, help="default <archive>/waves_marconi.csv")
    ap.add_argument("--gnssr-spline",
                    default="/home/argus_user/GNSS/v4.1/products/refl_code/Files/usgs/usgs_spline_out.txt",
                    help="the station's GNSS-R water level, used first (Chatham fills the newest hours)")
    ap.add_argument("--image-dir", nargs="+",
                    default=["/mnt/I2Rgus_Data/ImageProducts/products", "/mnt/I2Rgus_Data/ImageProducts"])
    ap.add_argument("--last-run", default=str(HERE / "logs" / "owg_last_run"),
                    help="stamp owg.sh writes after each run")
    ap.add_argument("--days", type=int, default=7)
    ap.add_argument("--out-dir", default=str(HERE / "reports" / "owg"))
    ap.add_argument("--attach", action="append", default=[],
                    help="extra image to attach (repeatable; glob patterns allowed, missing is noted)")
    ap.add_argument("--email", nargs="+", default=None, metavar="ADDRESS")
    ap.add_argument("--subject-prefix", default="Marconi Beach OWG")
    ap.add_argument("--dry-run", action="store_true", help="with --email: print, do not send")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")

    arch = Path(args.archive)
    waves = load_waves(args.waves_csv or arch / "waves_44008.csv")
    gauge = load_gauge(args.gauge_csv or arch / "gauge_8447435.csv")
    from marconi_water_level import WaterLevel
    fit = WaterLevel(args.gnssr_spline, args.gauge_csv or arch / "gauge_8447435.csv")
    models, stems = [], []
    for spec in args.models:
        label, name, *stem = spec.split("=")
        stem = stem[0] if stem else None
        if stem and not Path(stem).is_absolute():
            stem = str(HERE / "owg_models" / stem)
        if label == "Run C" and args.model_c:
            stem = args.model_c
        if label == "patch" and args.model_patch:
            stem = args.model_patch
        d = load_archive(arch / name)
        if d is None and not (stem and Path(stem + ".onnx").exists()):
            continue
        models.append((label, d)); stems.append(stem)
    # ensemble: the mean of several models' readings of the same frame. Models
    # trained differently err differently, so their average errs less.
    ens = [l.strip() for l in args.ensemble.split(",") if l.strip()]
    parts = [d[d["status"] == "ok"][["filename", "epoch", "t", "hs_m"]].rename(columns={"hs_m": lbl})
             for lbl, d in models if lbl in ens and d is not None and len(d)]
    if len(parts) >= 2:
        e = parts[0]
        for p in parts[1:]:
            e = e.merge(p[["filename", p.columns[-1]]], on="filename", how="outer")
        e["epoch"] = e["filename"].str.split(".").str[0].astype(int)
        e["t"] = pd.to_datetime(e["epoch"], unit="s", utc=True)
        cols = [c for c in e.columns if c in ens]
        n = e[cols].notna().sum(axis=1)
        e = e[n >= 2].copy()
        e["hs_m"] = e[cols].mean(axis=1)
        e["status"] = "ok"
        if len(e):
            models.append(("ensemble", e.sort_values("epoch").reset_index(drop=True))); stems.append(None)

    now = datetime.now(timezone.utc)
    end = now
    start = now - timedelta(days=args.days)
    day_ago = now.timestamp() - 86400
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # hourly best estimate and the buoy converted to Marconi (marconi_waves.py)
    mwp = Path(args.marconi_waves or arch / "waves_marconi.csv")
    mw = pd.read_csv(mwp) if mwp.exists() else None
    conv = mw[mw["hs_buoy_marconi"].notna()] if mw is not None else None
    converted = conv is not None and len(conv) >= 2

    # every OK frame matched to the buoy hour, the tide and the wave climate
    matched = []
    for label, d in models:
        if d is None or not len(d):
            matched.append((label, None))
            continue
        ok = d[d["status"] == "ok"].copy()
        if waves is not None:
            for col, name in (("wvht_m", "buoy"), ("dpd_s", "dpd"), ("mwd_deg", "mwd")):
                ok[name] = interp_at(waves["epoch"], waves[col], ok["epoch"], 90 * 60)
        else:
            ok["buoy"] = ok["dpd"] = ok["mwd"] = np.nan
        if converted:   # score against the buoy as it would be at Marconi, when fitted
            ok["buoy"] = interp_at(conv["epoch"], conv["hs_buoy_marconi"], ok["epoch"], 90 * 60)
        ok["level"] = fit.at(ok["epoch"].to_numpy(float))[0]
        matched.append((label, ok[ok["buoy"].notna()]))

    # ---- checks: every link of the chain, newest first
    checks = []

    def check(ok, text):
        checks.append((bool(ok), text))

    if waves is None or not waves["wvht_m"].notna().any():
        check(False, "buoy 44008: no wave archive (fetch_buoy_waves.py)")
    else:
        last = waves[waves["wvht_m"].notna()]["epoch"].iloc[-1]
        check(now.timestamp() - last < 6 * 3600, f"buoy 44008: latest wave height {ago(now.timestamp() - last)} old")
    if gauge is None or not len(gauge):
        check(False, "tide gauge 8447435: no archive (fetch_tide_gauge.py)")
    else:
        last = gauge["epoch"].iloc[-1]
        check(now.timestamp() - last < 6 * 3600, f"tide gauge 8447435: latest reading {ago(now.timestamp() - last)} old")
    lr = Path(args.last_run)
    if lr.exists():
        age = now.timestamp() - lr.stat().st_mtime
        msg = lr.read_text().strip().splitlines()[0] if lr.read_text().strip() else ""
        check(age < 26 * 3600 and "error" not in msg.lower(),
              f"owg.sh: last run {ago(age)} ago" + (f" ({msg})" if msg else ""))
    else:
        check(False, "owg.sh: has never run (no logs/owg_last_run)")
    for (label, d), stem in zip(models, stems):
        if stem is None:
            continue
        need = [stem + ".onnx", stem + ".report.json"] + ([stem + ".patch.json"] if "patch" in label else [])
        missing = [Path(p).name for p in need if not Path(p).exists()]
        if missing:
            check(False, f"camera {label}: model files missing: {', '.join(missing)}")
            continue
        if d is None or not len(d):
            check(False, f"camera {label}: no frames processed yet")
            continue
        recent = d[d["epoch"] >= day_ago]
        n_ok = int((recent["status"] == "ok").sum())
        why = recent.loc[recent["status"] != "ok", "status"].value_counts()
        detail = ", ".join(f"{v} {k.split(' (')[0]}" for k, v in why.items())
        check(len(recent) > 0 and n_ok > 0,
              f"camera {label}: {n_ok} of {len(recent)} frames measured in 24 h"
              + (f" (screened: {detail})" if detail else ""))
    # the station's GPS water level (GNSS-R): a day or two behind is normal (orbits)
    if fit.s_ep is not None and len(fit.s_ep):
        age = now.timestamp() - fit.s_ep[-1]
        check(age < 4 * 86400, f"GPS water level: latest reading {ago(age)} old (1-2 days is normal)")
    else:
        check(False, f"GPS water level: not available ({fit.gps_note}); using the Chatham gauge")
    # camera pointing (pointing_check.py), when a reference bank is set up
    for cam in ("c1", "c2"):
        pl = arch / f"pointing_{cam}.csv"
        if not pl.exists():
            continue
        p = pd.read_csv(pl)
        p = p[p["epoch"] >= now.timestamp() - 86400]
        m = p[p["status"] != "unmatched"]
        if not len(p):
            check(False, f"pointing {cam}: no frame checked in 24 h")
        elif not len(m):
            check(False, f"pointing {cam}: no frame could be matched in 24 h (fog/rain, or the view changed?)")
        else:
            last = m.iloc[-1]
            worst = m[["d_azimuth", "d_tilt", "d_roll"]].abs().max()
            moved = (m["status"].str.lower() == "moved").sum()
            check(moved == 0, f"pointing {cam}: latest change az {last['d_azimuth']:+.2f}, tilt "
                              f"{last['d_tilt']:+.2f}, roll {last['d_roll']:+.2f} deg; largest 24 h "
                              f"{worst.max():.2f} deg ({len(m)} frames"
                              + (f", {moved} beyond the limit" if moved else "") + ")")
    stuck = [lbl for lbl, d in models if d is not None and len(d) and
             (d[d["epoch"] >= day_ago]["status"] == "no water level").sum() > 3]
    for lbl in stuck:
        check(False, f"camera {lbl}: frames waiting for a water level -- is the gauge download working?")

    # ---- plots
    plots = []
    p = out_dir / "owg_report_7day.png"
    plot_7day(models, waves, gauge, fit, start, end, p, mw); plots.append(p)
    p = out_dir / "owg_report_scatter.png"
    if plot_scatter(matched, start, p):
        plots.append(p)
    p = out_dir / "owg_report_residuals.png"
    if plot_residuals(matched, p):
        plots.append(p)
    p = out_dir / "owg_report_frames.png"
    if plot_frames([(l, d) for l, d in models if l != "ensemble"], start, end, p):
        plots.append(p)
    p = out_dir / "owg_report_view.png"
    try:
        stem_of = dict(zip([l for l, _ in models], stems))
        view_time = plot_view(stem_of.get("Run C") or str(HERE / "owg_models" / "owg_c2_H_current_C"),
                              stem_of.get("patch") or str(HERE / "owg_models" / "owg_c2_H_patch"),
                              args.image_dir, p)
    except Exception as exc:
        view_time = None
        print(f"view: not drawn ({exc})", file=sys.stderr)
    if view_time:
        plots.append(p)

    # ---- text
    L = []
    L.append("=" * 68)
    L.append(f"  Optical wave gauge -- {now:%Y-%m-%d %H:%M}Z")
    L.append("=" * 68)
    L.append("")
    L.append("CONDITIONS (buoy 44008, ~50 nm offshore)")
    if waves is not None and waves["wvht_m"].notna().any():
        w = waves[waves["wvht_m"].notna()]
        r = w.iloc[-1]
        # NDBC often leaves period and direction blank ("MM") in the newest
        # row; give the latest that were reported, and their time if older
        bits = [f"Hs {r['wvht_m']:.1f} m"]
        for col, fmt in (("dpd_s", "period {:.0f} s"), ("mwd_deg", "from {:.0f} deg")):
            x = w[w[col].notna()]
            if len(x):
                v = x.iloc[-1]
                txt = fmt.format(v[col]) + (f" {compass(v[col])}" if col == "mwd_deg" else "")
                if v["epoch"] < r["epoch"]:
                    txt += f" (at {v['t']:%H:%M}Z)"
                bits.append(txt)
        L.append(f"  Latest      : {', '.join(bits)}   ({r['t']:%m-%d %H:%M}Z)")
        for hrs, name in ((24, "Last 24 h"), (24 * args.days, f"Last {args.days} days")):
            x = w[w["epoch"] >= now.timestamp() - hrs * 3600]
            if len(x):
                pk = x.loc[x["wvht_m"].idxmax()]
                L.append(f"  {name:<12}: Hs {x['wvht_m'].min():.1f}-{x['wvht_m'].max():.1f} m, "
                         f"mean {x['wvht_m'].mean():.1f} m; peak {pk['wvht_m']:.1f} m at "
                         f"{pk['t']:%m-%d %H:%M}Z")
    else:
        L.append("  no buoy data")
    if mw is not None and "tp_camera_s" in mw and mw["tp_camera_s"].notna().any():
        c = mw[mw["tp_camera_s"].notna()].iloc[-1]
        L.append(f"  Camera Tm01 : {c['tp_camera_s']:.1f} s mean wave period, c2 timestacks, "
                 f"{c['time_utc'][5:16].replace('T', ' ')}Z")
    if gauge is not None and len(gauge):
        t = np.arange(now.timestamp() - 86400, now.timestamp(), 600.0)
        lv, lsrc = fit.at(t)
        if np.isfinite(lv).any():
            gps_share = float(np.mean(lsrc[np.isfinite(lv)] == "gps"))
            L.append(f"  Tide (24 h) : Marconi {np.nanmin(lv):+.2f} to {np.nanmax(lv):+.2f} m NAVD88 "
                     f"({gps_share:.0%} from the station's GPS, the rest Chatham fitted to the GPS)")
            L.append(f"  {fit.describe()}")
    L.append("")
    L.append("CAMERA (c2): significant wave height in the surf zone")
    for label, d in models:
        if d is None or not len(d):
            L.append(f"  {label:<7}: no data")
            continue
        ok = d[d["status"] == "ok"]
        if not len(ok):
            L.append(f"  {label:<7}: no measured frames")
            continue
        last = ok.iloc[-1]
        t, v = averaged(ok)
        avg = f", 5-frame average {v[-1]:.2f} m" if len(v) and t[-1] >= ok["epoch"].iloc[-1] - 3 * 3600 else ""
        L.append(f"  {label:<7}: latest {last['hs_m']:.2f} m at {last['t']:%m-%d %H:%M}Z{avg}")
        x = ok[ok["epoch"] >= day_ago]
        if len(x):
            L.append(f"  {'':<7}  last 24 h {x['hs_m'].min():.2f}-{x['hs_m'].max():.2f} m, "
                     f"mean {x['hs_m'].mean():.2f} m ({len(x)} frames)")
    L.append("")
    if mw is not None and mw["hs_best"].notna().any():
        b = mw[mw["hs_best"].notna()]
        r = b.iloc[-1]
        L.append("BEST ESTIMATE AT MARCONI (camera blended with the converted buoy, hourly)")
        L.append(f"  Latest      : Hs {r['hs_best']:.2f} m at {r['time_utc'][5:16].replace('T', ' ')}Z "
                 f"({r['best_from']}: camera {r['hs_camera'] if pd.notna(r['hs_camera']) else '--'}, "
                 f"converted buoy {r['hs_buoy_marconi'] if pd.notna(r['hs_buoy_marconi']) else '--'} m)")
        x = b[b["epoch"] >= day_ago]
        if len(x):
            L.append(f"  Last 24 h   : {x['hs_best'].min():.2f}-{x['hs_best'].max():.2f} m, mean "
                     f"{x['hs_best'].mean():.2f} m  (archive/waves_marconi.csv, used by the other systems)")
        if not converted:
            L.append("  (buoy not converted yet -- run: python3 buoy_transfer.py fit --station 44013)")
        L.append("")
    L.append("AGREEMENT WITH THE BUOY " + ("CONVERTED TO MARCONI (by wave direction, period and height)"
                                         if converted else "(single frames, raw offshore buoy)"))
    L.append(f"  {'':<20}{'n':>5}{'bias':>9}{'RMSE':>8}{'r':>7}{'d':>7}")
    flat = False
    for label, m in matched:
        if m is None or not len(m):
            L.append(f"  {label:<20}  no frames matched to the buoy")
            continue
        for name, sel in (("24 h", m["epoch"] >= day_ago),
                          (f"{args.days} days", m["epoch"] >= start.timestamp()),
                          ("full record", np.ones(len(m), bool))):
            s = metrics(m.loc[sel, "hs_m"].to_numpy(float), m.loc[sel, "buoy"].to_numpy(float))
            row = f"  {label + ', ' + name:<20}{s['n']:>5}"
            if s["n"] >= 3:
                row += f"{s['bias']:>+8.2f}m{s['rmse']:>7.2f}m"
                # r and d measure whether the two rise and fall together; with
                # the buoy nearly steady there is nothing to follow, and they
                # come out at random -- even negative
                spread = np.ptp(m.loc[sel, "buoy"].to_numpy(float))
                if spread >= MIN_SPREAD and s["n"] >= 10:
                    row += f"{s['corr']:>7.2f}{s['d']:>7.2f}"
                else:
                    row += f"{'--':>7}{'--':>7}"
                    flat = True
            L.append(row)
    if flat:
        L.append(f"  -- : r and d not given: the buoy varied less than {MIN_SPREAD:.1f} m (or under")
        L.append("       10 frames), too little for them to mean anything.")
    if converted:
        L.append("  The converted buoy is itself uncertain (about +/-0.26-0.29 m against the")
        L.append("  ADCP), so camera and buoy differing by less than ~0.3 m is agreement.")
        L.append("  A trend left in owg_report_residuals.png is something neither explains.")
    else:
        L.append("  The buoy is offshore and the camera nearshore, so they are not")
        L.append("  expected to match: watch r and d (do they rise and fall together?).")
        L.append("  A bias that grows with the waves is breaking; one that changes with")
        L.append("  direction is sheltering (see owg_report_residuals.png).")
    L.append("  The camera's own accuracy, against the ADCP: Run C 0.27 m per frame, patch 0.38 m.")
    have = [(l, m) for l, m in matched if m is not None and len(m)]
    for lbl, m in have[1:]:
        both = have[0][1].merge(m, on="filename", suffixes=("_a", "_b"))
        if len(both) >= 3:
            dd = both["hs_m_b"] - both["hs_m_a"]
            L.append(f"  {lbl} minus {have[0][0]} on the same {len(both)} frames: "
                     f"{dd.mean():+.2f} m (sd {dd.std():.2f} m)")
    L.append("")
    n_fail = sum(1 for ok, _ in checks if not ok)
    L.append("CHECKS")
    for ok, text in checks:
        L.append(f"  [{' OK ' if ok else 'FAIL'}] {text}")
    L.append("")
    L.append("PLOTS (attached)")
    descr = {"owg_report_7day.png": "wave height, period, direction and tide, last 7 days",
             "owg_report_scatter.png": "camera vs buoy, each model",
             "owg_report_residuals.png": "camera minus buoy vs tide, direction, period",
             "owg_report_frames.png": "frames per day: measured or why screened out",
             "owg_report_view.png": "what each model sees, latest frame"
                                     + (f" ({view_time:%H:%M}Z)" if view_time else "")}
    for p in plots:
        L.append(f"  {p.name:<26} {descr.get(p.name, '')}")
    extra, missing = [], []
    for pat in args.attach:
        hits = sorted(Path("/").glob(pat.lstrip("/"))) if any(c in pat for c in "*?[") else \
            ([Path(pat)] if Path(pat).exists() else [])
        if hits:
            extra.append(hits[-1])
            L.append(f"  {hits[-1].name:<26} (from {hits[-1].parent})")
        else:
            missing.append(pat)
    if missing:
        L.append("  not attached, file missing: " + ", ".join(missing))
    L.append("")
    L.append(f"Report and plots: {out_dir}")
    L.append("Rerun on the station: /mnt/I2Rgus_Data/waterline/owg.sh   (add --email to send)")
    text = "\n".join(L)

    (out_dir / "owg_report.txt").write_text(text + "\n")
    stamp = out_dir / "daily"
    stamp.mkdir(exist_ok=True)
    (stamp / f"owg_report_{now:%Y-%m-%d}.txt").write_text(text + "\n")
    print(text)

    if args.email:
        status = "all healthy" if n_fail == 0 else f"{n_fail} problem(s) need attention"
        subject = f"{args.subject_prefix}: {status}"
        opening = ("Surf's up -- or it isn't, and the camera can tell you which!\n\n"
                   "Every check passed: buoy, tide gauge, image screening and both\n"
                   "wave models are working as expected." if n_fail == 0 else
                   f"The optical wave gauge reported {n_fail} problem(s) this morning.\n"
                   "Details under CHECKS below.")
        body = opening + "\n\n" + text + "\n\n--\nSent by owg_report.py on " + \
            subprocess.run(["hostname"], capture_output=True, text=True).stdout.strip() + \
            f" at {now:%Y-%m-%d %H:%M}Z."
        rc = send(subject, args.email, body, plots + extra, args.dry_run)
        if rc:
            return 100 + rc
    return n_fail


def send(subject, to, body, attachments, dry_run):
    from email.message import EmailMessage
    msg = EmailMessage()
    msg["Subject"] = subject
    msg["To"] = ", ".join(to)
    msg.set_content(body)
    for p in attachments:
        data = Path(p).read_bytes()
        sub = "jpeg" if str(p).lower().endswith((".jpg", ".jpeg")) else "png"
        msg.add_attachment(data, maintype="image", subtype=sub, filename=Path(p).name)
    size = len(msg.as_bytes()) / 1e6
    if dry_run:
        print(f"\n--- would send to: {' '.join(to)} ---\nSubject: {subject}\n"
              f"Attachments: {[Path(p).name for p in attachments]} ({size:.1f} MB)")
        return 0
    # retry: one network blip should not turn into a silent morning
    for attempt, wait in enumerate((30, 60, 0), 1):
        try:
            r = subprocess.run(["msmtp"] + list(to), input=msg.as_bytes(), capture_output=True)
        except FileNotFoundError:
            print("msmtp not found -- cannot send", file=sys.stderr)
            return 1
        if r.returncode == 0:
            print(f"\nemail sent to {' '.join(to)} ({size:.1f} MB)")
            return 0
        print(f"msmtp attempt {attempt} failed: {r.stderr.decode(errors='replace').strip()}",
              file=sys.stderr)
        if wait:
            time.sleep(wait)
    return 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Tests For The Wave-Gauge Plots And The Drift Check
====================================================
Plain python3 (no pytest needed on the station):

    python3 test_owg_plots.py

WHAT.
  * owg_plot_marks.py: the sun (sunrise and sunset at Boston against the
    almanac, to 3 minutes; the equinox noon at Marconi), the night spans,
    the frames' interval, the screened spans (daylight only, half an
    interval either side, touching ones merged, clipped), the legend order,
    the model's tested error (holdback first, the mean for several models,
    nothing when a report or value is missing) and its band;
  * owg_live.py's plot: the 2-hour average is drawn in a week with no
    low-sharpness frame (it was drawn only when there was one), with its
    band, and night and screened shading in the legend after the data;
  * owg_report.py: the headline model (composite, else ensemble, else Run C,
    else the first with data), the model stems behind an archive
    (model_id), the 7-day plot (the others grey, one legend entry, the
    low-sharpness frames hollow, a headline past the palette still drawn),
    the drift check (a 0.3 m shift flagged; a small one, a noisy one, too
    few days, and a calm week after a rough month for a camera that reads
    big waves low, not), and the whole report run end to end on a
    synthetic archive: the drift lines in its checks, and the ensemble as
    the headline with its members' error band and screened shading.
Exit status 0 when every check passes.
"""

import ast
import contextlib
import io
import json
import sys
import tempfile
from datetime import datetime, timezone, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import matplotlib            # noqa: E402
matplotlib.use("Agg")
import owg_plot_marks as M   # noqa: E402

failures = []
BOSTON = (42.3601, -71.0589)


def check(name, cond, detail=""):
    print(f"  {'ok ' if cond else 'FAIL'} {name}" + (f"  ({detail})" if detail != "" and not cond else ""))
    if not cond:
        failures.append(name)


def ep(text):
    return datetime.strptime(text, "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc).timestamp()


def hhmm(e):
    return datetime.fromtimestamp(e, tz=timezone.utc).strftime("%H:%M")


@contextlib.contextmanager
def figures():
    """Keeps every figure a plot function saves (and closes) for inspection."""
    from matplotlib.figure import Figure
    kept, save = [], Figure.savefig

    def keep(self, *a, **k):
        kept.append(self)
        return save(self, *a, **k)

    Figure.savefig = keep
    try:
        yield kept
    finally:
        Figure.savefig = save


def legend_texts(fig_or_ax):
    leg = fig_or_ax.get_legend() if hasattr(fig_or_ax, "get_legend") else None
    if leg is None and hasattr(fig_or_ax, "legends") and fig_or_ax.legends:
        leg = fig_or_ax.legends[0]
    return [t.get_text() for t in leg.get_texts()] if leg is not None else []


def test_marks():
    print("owg_plot_marks.py")
    # Boston, 20-21 June 2024: sunset 8:24 pm EDT (00:24 UTC), sunrise 5:07 am EDT (09:07 UTC)
    sp = M.night_spans(ep("2024-06-20 12:00"), ep("2024-06-21 12:00"), 60.0, *BOSTON)
    check("one summer night at Boston", len(sp) == 1, sp)
    if sp:
        a, b = sp[0]
        check("summer sunset 00:24 and sunrise 09:07 UTC, to 3 min",
              abs(a - ep("2024-06-21 00:24")) <= 180 and abs(b - ep("2024-06-21 09:07")) <= 180, (hhmm(a), hhmm(b)))
    # 21-22 Dec 2024: sunset 4:14 pm EST (21:14 UTC), sunrise 7:10 am EST (12:10 UTC)
    sp = M.night_spans(ep("2024-12-21 18:00"), ep("2024-12-22 18:00"), 60.0, *BOSTON)
    if sp:
        a, b = sp[0]
        check("winter sunset 21:14 and sunrise 12:10 UTC, to 3 min",
              len(sp) == 1 and abs(a - ep("2024-12-21 21:14")) <= 180 and abs(b - ep("2024-12-22 12:10")) <= 180,
              (hhmm(a), hhmm(b)))
    else:
        check("a winter night at Boston", False)
    t = np.arange(ep("2026-03-20 14:00"), ep("2026-03-20 19:00"), 60.0)
    el = M.sun_elevation_deg(t)
    check("equinox noon at Marconi: 90 - 41.89 = 48.1 deg, near 16:47 UTC",
          abs(el.max() - 48.1) < 0.3 and abs(t[el.argmax()] - ep("2026-03-20 16:47")) <= 300,
          (round(float(el.max()), 2), hhmm(t[el.argmax()])))
    check("spans clipped to the window", all(ep("2026-10-28 00:00") <= a < b <= ep("2026-10-28 12:00")
                                            for a, b in M.night_spans(ep("2026-10-28 00:00"), ep("2026-10-28 12:00"))))

    check("frames every 30 min: interval 1800 s", M.frame_spacing(np.arange(0, 36000, 1800.0)) == 1800.0)
    check("one frame: the default 1800 s", M.frame_spacing([5.0]) == 1800.0)
    check("the interval kept within 10-60 min", M.frame_spacing([0, 60, 120]) == 600.0 and
          M.frame_spacing([0, 20000, 40000]) == 3600.0)

    day = ep("2026-10-27 15:00")                      # mid-afternoon at Marconi: daylight
    frames = np.array([day - 3600, day - 1800, day, day + 1800, day + 3600])
    st = np.array(["ok", "blurred (fog, rain or wet lens)", "glare", "ok", "no water level"], object)
    sp = M.screened_spans(frames, st, day - 7200, day + 7200)
    check("two touching screened daylight frames: one span of an hour",
          sp == [(day - 2700, day + 900)], sp)
    check("a status not about the view (no water level) is not shaded", all(b < day + 3600 for _, b in sp))
    night = ep("2026-10-28 04:00")
    sp = M.screened_spans(np.array([night - 1800, night, night + 1800]),
                          np.array(["dark", "dark", "dark"], object), night - 7200, night + 7200)
    check("dark frames at night: left to the night shade", sp == [], sp)
    sp = M.screened_spans(frames, st, day - 1000, day + 500)
    check("screened spans clipped to the window", sp == [(day - 1000, day + 500)], sp)
    check("no frames: no spans", M.screened_spans([], [], 0, 1) == [])

    h, lab = M.background_last(["n", "b", "s", "c"], [M.NIGHT_LABEL, "buoy", M.SCREENED_LABEL, "camera"])
    check("legend: the data first, then night, then screened",
          lab == ["buoy", "camera", M.NIGHT_LABEL, M.SCREENED_LABEL] and h == ["b", "c", "n", "s"], lab)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "a.report.json").write_text(json.dumps({"holdback_rmse": 0.32, "val_rmse": 0.25}))
        (tmp / "b.report.json").write_text(json.dumps({"holdback_rmse": 0.34}))
        (tmp / "c.report.json").write_text(json.dumps({"val_rmse": 0.27}))
        (tmp / "d.report.json").write_text(json.dumps({"arch": "x"}))
        (tmp / "e.report.json").write_text("not json")
        r = lambda *n: M.model_rmse([str(tmp / x) for x in n])   # noqa: E731
        check("tested error: the holdback RMSE, not the validation one", r("a") == 0.32, r("a"))
        check("two models averaged: the mean of theirs", abs(r("a", "b") - 0.33) < 1e-12, r("a", "b"))
        check("an older report: its validation RMSE", r("c") == 0.27, r("c"))
        check("no RMSE, an unreadable or a missing report: none (no band)",
              r("d") is None and r("e") is None and r("a", "missing") is None and M.model_rmse([]) is None)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ts = np.array(["2026-10-27T12:00", "2026-10-27T13:00"], dtype="datetime64[s]")
    M.error_band(ax, [(ts, np.array([0.1, 1.0]))], None, "b")
    check("no tested error: no band", len(ax.collections) == 0)
    M.error_band(ax, [(ts, np.array([0.1, 1.0]))], 0.3, "b")
    lo = ax.collections[0].get_paths()[0].vertices[:, 1].min() if ax.collections else None
    check("the band floored at 0 m", lo is not None and abs(lo) < 1e-12, lo)
    plt.close(fig)


def live_frames(now, low_sharpness):
    """A week of 30-min frames: dark at night, a glare frame each morning, the rest ok (some low-sharpness)."""
    import owg_live as L
    t0 = (now - timedelta(days=7)).timestamp()
    e = np.arange(int(t0 // 1800) * 1800, now.timestamp(), 1800.0)
    el = M.sun_elevation_deg(e)
    st = np.where(el < 0, "dark", "ok").astype(object)
    day_i = np.flatnonzero(el >= 0)
    st[day_i[::15]] = "glare"
    if low_sharpness:
        st[(st == "ok") & (np.arange(len(e)) % 7 == 0)] = L.LOW_SHARPNESS
    hs = np.where(np.isin(st, ["ok", L.LOW_SHARPNESS]), 1.0 + 0.2 * np.sin(e / 20000.0), np.nan)
    return pd.DataFrame({"epoch": e.astype(int), "status": st, "hs_m": hs})


def test_live_plot():
    print("owg_live.py: the 7-day plot")
    import owg_live as L
    now = datetime.now(timezone.utc)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        w = pd.DataFrame({"epoch": np.arange(now.timestamp() - 8 * 86400, now.timestamp(), 3600.0)})
        w["wvht_m"] = 1.2
        w.to_csv(tmp / "w.csv", index=False)
        for low in (False, True):
            out, buf = tmp / f"p{int(low)}.png", io.StringIO()
            with figures() as kept, contextlib.redirect_stdout(buf):
                L.plot(live_frames(now, low), str(tmp / "w.csv"), str(out), 7, 5, 45.0, "m", "s", rmse=0.3)
            ax = kept[0].axes[0] if kept else None
            texts = legend_texts(ax) if ax is not None else []
            avg = [t for t in texts if "5-frame average" in t]
            tag = "with" if low else "without"
            check(f"a week {tag} low-sharpness frames: the 2-hour average is drawn",
                  out.exists() and bool(avg) and "averaged points" in buf.getvalue()
                  and " 0 averaged points" not in buf.getvalue(), (texts, buf.getvalue()))
            check(f"... {tag}: its band, then night and screened last in the legend",
                  any("tested error per frame, ±0.30 m" in t for t in texts)
                  and texts[-2:] == [M.NIGHT_LABEL, M.SCREENED_LABEL], texts)
        with figures() as kept, contextlib.redirect_stdout(io.StringIO()):
            L.plot(live_frames(now, False), str(tmp / "w.csv"), str(tmp / "n.png"), 7, 5, 45.0, "m", "s")
        texts = legend_texts(kept[0].axes[0]) if kept else []
        check("no tested error known: no band, the rest drawn", not any("tested error" in t for t in texts)
              and any("5-frame average" in t for t in texts), texts)


def report_models(now, n_days=7):
    t0 = (now - timedelta(days=n_days)).timestamp()
    e = np.arange(int(t0 // 1800) * 1800, now.timestamp(), 1800.0)
    st = np.where(M.sun_elevation_deg(e) < 0, "dark", "ok").astype(object)
    st[np.flatnonzero(st == "ok")[::12]] = "blurred (fog, rain or wet lens)"
    out = []
    for k, label in enumerate(("Run C", "patch", "v2", "composite")):
        d = pd.DataFrame({"epoch": e.astype(int), "status": st, "hs_m": np.where(st == "ok", 1.0 + 0.05 * k, np.nan)})
        d["t"] = pd.to_datetime(d["epoch"], unit="s", utc=True)
        out.append((label, d))
    return out


def test_report_parts():
    print("owg_report.py: headline, model stems, 7-day plot")
    import owg_report as R
    now = datetime.now(timezone.utc)
    s0 = (now - timedelta(days=7)).timestamp()
    models = report_models(now)
    empty = models[0][1].iloc[0:0]
    check("headline: the composite models first", R.headline_of(models, s0) == "composite")
    check("... else the ensemble, else Run C",
          R.headline_of(models[:3] + [("ensemble", models[0][1])], s0) == "ensemble"
          and R.headline_of(models[:3], s0) == "Run C")
    check("... else the first model with data; none without",
          R.headline_of([("patch", models[1][1]), ("v2", models[2][1])], s0) == "patch"
          and R.headline_of([("Run C", empty), ("v2", None)], s0) is None)
    d = pd.DataFrame({"model_id": [np.nan, "owg_c2_H_composite:abc+owg_c2_H_composite_s2:def"]})
    check("the stems behind an archive: its model_id, next to the configured stem",
          R.members_of(d, "/m/owg_c2_H_composite") == ["/m/owg_c2_H_composite", "/m/owg_c2_H_composite_s2"])
    check("... no model_id: the configured stem; no stem: none",
          R.members_of(pd.DataFrame({"x": [1]}), "/m/a") == ["/m/a"] and R.members_of(d, None) == [])

    class Fit:
        def at(self, t):
            t = np.asarray(t, float)
            return np.sin(t / 7000.0), np.full(len(t), "gps", object)

    w = pd.DataFrame({"epoch": np.arange(s0 - 86400, now.timestamp(), 3600.0)})
    w["wvht_m"], w["dpd_s"], w["mwd_deg"] = 1.2, 9.0, 90.0
    with tempfile.TemporaryDirectory() as tmp, figures() as kept:
        R.plot_7day(models, w, None, Fit(), now - timedelta(days=7), now, Path(tmp) / "r.png", None, rmse=0.32)
        texts = legend_texts(kept[0]) if kept else []
    scr = models[3][1].copy()
    scr.loc[scr.index[scr["status"] == "ok"][::5], "status"] = "low sharpness (may read low)"
    six = models + [("a", models[0][1]), ("ensemble", models[1][1])]
    with tempfile.TemporaryDirectory() as tmp, figures() as kept2:
        R.plot_7day(models, w, None, Fit(), now - timedelta(days=7), now, Path(tmp) / "r.png", None, 0.32, scr)
        R.plot_7day(six[:4] + six[5:6] + [("b", models[0][1])] + six[4:5], w, None, Fit(),
                    now - timedelta(days=7), now, Path(tmp) / "s.png", None, None)
    t_low = legend_texts(kept2[0]) if kept2 else []
    t_six = legend_texts(kept2[1]) if len(kept2) > 1 else []
    check("7-day plot: the headline's low-sharpness frames drawn hollow, in the legend",
          "view degraded (may read low; not used)" in t_low, t_low)
    check("... a headline past the five palette colours is still drawn",
          any(t.startswith("Camera composite") for t in t_six) or any(t.startswith("Camera ensemble") for t in t_six),
          t_six)
    check("7-day plot: the headline with its band, the others one grey entry, the shades last",
          "Camera composite (2-h mean; dots: frames)" in texts
          and "other camera models (Run C, patch, v2), 2-h means" in texts
          and "tested error per frame, ±0.32 m" in texts and texts[-2:] == [M.NIGHT_LABEL, M.SCREENED_LABEL]
          and not any(t.startswith("Camera Run C") for t in texts), texts)


def drift_frames(now, shift=0.0, noise=0.03, days=35, per_day=10, recent=7, seed=1):
    """Matched frames: camera - buoy = 0.1 m (+ shift over the last `recent` days) + daily noise."""
    rng = np.random.default_rng(seed)
    today = int(now.timestamp() // 86400)
    rows = []
    for k in range(days):
        dday = today - k
        off = 0.1 + (shift if k < recent else 0.0) + rng.normal(0, noise)
        for j in range(per_day):
            e = dday * 86400 + 12 * 3600 + 1800 * j
            rows.append((e, 1.0 + off, 1.0))
    return pd.DataFrame(rows, columns=["epoch", "hs_m", "buoy"])


def shrink_frames(now, step=0.0, k=0.55, seed=4):
    """A camera that reads big waves low and small ones high: camera = 0.5 + k Hs (+ step over the last
    7 days); the buoy rough (1.0-1.6 m) for four weeks, then a calm week (0.5-0.6 m)."""
    rng = np.random.default_rng(seed)
    today = int(now.timestamp() // 86400)
    rows = []
    for kday in range(35):
        h = 0.55 + 0.05 * np.sin(kday) if kday < 7 else 1.3 + 0.3 * np.sin(kday * 0.9)
        for j in range(10):
            hs = h + rng.normal(0, 0.05)
            cam = 0.5 + k * hs + (step if kday < 7 else 0.0) + rng.normal(0, 0.1)
            rows.append(((today - kday) * 86400 + 12 * 3600 + 1800 * j, cam, hs))
    return pd.DataFrame(rows, columns=["epoch", "hs_m", "buoy"])


def test_drift():
    print("owg_report.py: the drift check")
    import owg_report as R
    now = datetime.now(timezone.utc)
    n = now.timestamp()
    r = R.drift(drift_frames(now), n)
    check("a steady bias: not flagged", r is not None and not r["flag"] and abs(r["shift"]) < 0.05, r)
    r = R.drift(drift_frames(now, shift=0.3), n)
    check("the bias moved 0.3 m in the last week: flagged", r is not None and r["flag"]
          and abs(r["shift"] - 0.3) < 0.05 and r["n_recent"] == 7 and r["n_base"] == 28, r)
    r = R.drift(drift_frames(now, shift=-0.3), n)
    check("... down as well as up", r is not None and r["flag"] and r["shift"] < -0.25, r)
    r = R.drift(drift_frames(now, shift=0.1, noise=0.005), n)
    check("a clear but small shift (0.1 m < 0.15 m): not flagged", r is not None and not r["flag"]
          and r["z"] > 3, r)
    r = R.drift(drift_frames(now, shift=0.3, noise=0.4, seed=3), n)
    check("a 0.3 m shift inside day-to-day noise of 0.4 m: not flagged (under 3 standard errors)",
          r is not None and r["z"] < 3 and not r["flag"], r)
    check("too few days: not judged", R.drift(drift_frames(now, days=10), n) is None)
    sparse = drift_frames(now, shift=0.3, per_day=2)
    check("days with fewer than 3 matched frames do not count", R.drift(sparse, n) is None)
    f = drift_frames(now, shift=0.3)
    f.loc[f.index[:5], "buoy"] = np.nan
    check("frames without a buoy reading are left out", R.drift(f, n) is not None and R.drift(f, n)["flag"])
    check("no frames: not judged", R.drift(None, n) is None and R.drift(f.iloc[0:0], n) is None)
    r = R.drift(shrink_frames(now), n)
    plain = R.drift(drift_frames(now), n)
    check("a calm week after a rough month, the camera unchanged: not flagged (its bias read off the "
          "earlier weeks' line in wave height)", r is not None and not r["flag"] and abs(r["shift"]) < 0.1
          and abs(r["slope"] + 0.45) < 0.15 and r["expected"] > r["base"] + 0.15, r)
    check("... the same camera with a 0.3 m step in that calm week: flagged",
          R.drift(shrink_frames(now, step=0.3), n)["flag"], R.drift(shrink_frames(now, step=0.3), n))
    check("earlier weeks of one wave height (spread < 0.3 m): a plain comparison, slope 0",
          plain is not None and plain["slope"] == 0.0)


def test_report_end_to_end():
    print("owg_report.py end to end: a synthetic archive, the drift lines in the checks")
    now = datetime.now(timezone.utc)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        arch, models_dir = tmp / "archive", tmp / "models"
        arch.mkdir()
        models_dir.mkdir()
        e = np.arange(int((now.timestamp() - 40 * 86400) // 1800) * 1800, now.timestamp(), 1800.0)
        st = np.where(M.sun_elevation_deg(e) < 0, "dark", "ok").astype(object)
        day = (e // 86400).astype(int)
        recent = day > int(now.timestamp() // 86400) - 7
        for name, stem, jump in (("owg_c2_H.csv", "owg_c2_H_current_C", 0.0),
                                 ("owg_c2_H_composite.csv", "owg_c2_H_composite", 0.4)):
            hs = np.where(st == "ok", 1.0 + 0.1 * np.sin(day * 1.3) + np.where(recent, jump, 0.0), np.nan)
            pd.DataFrame({"filename": [f"{int(x)}.c2.bright.jpg" for x in e], "epoch": e.astype(int),
                          "time_utc": "", "status": st, "hs_m": hs, "model_id": f"{stem}:abc"}).to_csv(
                arch / name, index=False)
            (models_dir / f"{stem}.onnx").write_bytes(b"x")
            (models_dir / f"{stem}.report.json").write_text(json.dumps({"holdback_rmse": 0.3}))
        h = np.arange(e[0], now.timestamp(), 3600.0)
        pd.DataFrame({"epoch": h, "wvht_m": 1.0, "dpd_s": 9.0, "apd_s": 7.0, "mwd_deg": 90.0}).to_csv(
            arch / "waves_44008.csv", index=False)
        pd.DataFrame({"epoch": h, "level_navd88": np.sin(h / 7000.0)}).to_csv(arch / "gauge_8447435.csv", index=False)
        sys.argv = ["owg_report.py", "--archive", str(arch), "--out-dir", str(tmp / "out"),
                    "--models", f"Run C=owg_c2_H.csv={models_dir / 'owg_c2_H_current_C'}",
                    f"composite=owg_c2_H_composite.csv={models_dir / 'owg_c2_H_composite'}",
                    "--ensemble", "", "--gnssr-spline", str(tmp / "no_spline.txt"), "--image-dir", str(tmp),
                    "--last-run", str(tmp / "no_last_run")]
        import owg_report as R
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                rc = R.main()
            crash = None
        except Exception as exc:              # noqa: BLE001 -- a crash is the finding
            rc, crash = None, repr(exc)
        text = buf.getvalue()
        lines = [x for x in text.splitlines() if "bias against buoy 44008" in x]
        check("the report runs to the end", crash is None and (tmp / "out" / "owg_report_7day.png").exists(),
              crash or text[-600:])
        check("Run C steady: its drift line OK", any("[ OK ] camera Run C: bias" in x for x in lines), lines)
        check("composite jumped 0.4 m in the last week: its drift line FAIL, with what to do",
              any("[FAIL] camera composite: bias" in x and "drifting" in x for x in lines), lines)
        check("the exit status counts the drift failure with the others", isinstance(rc, int) and rc >= 1, rc)

        # no composite: the ensemble (Run C + v2) leads, with its members' error and screened frames
        e2 = e[e >= now.timestamp() - 8 * 86400]
        st2 = np.where(M.sun_elevation_deg(e2) < 0, "dark", "ok").astype(object)
        st2[np.flatnonzero(st2 == "ok")[::9]] = "blurred (fog, rain or wet lens)"
        for name, stem, err in (("owg_c2_H.csv", "owg_c2_H_current_C", 0.30), ("owg_c2_H_v2.csv", "owg_c2_H_v2", 0.34)):
            pd.DataFrame({"filename": [f"{int(x)}.c2.bright.jpg" for x in e2], "epoch": e2.astype(int),
                          "time_utc": "", "status": st2, "hs_m": np.where(st2 == "ok", 1.0, np.nan),
                          "model_id": f"{stem}:abc"}).to_csv(arch / name, index=False)
            (models_dir / f"{stem}.onnx").write_bytes(b"x")
            (models_dir / f"{stem}.report.json").write_text(json.dumps({"holdback_rmse": err}))
        sys.argv = ["owg_report.py", "--archive", str(arch), "--out-dir", str(tmp / "out2"),
                    "--models", f"Run C=owg_c2_H.csv={models_dir / 'owg_c2_H_current_C'}",
                    f"v2=owg_c2_H_v2.csv={models_dir / 'owg_c2_H_v2'}", "--ensemble", "Run C,v2",
                    "--gnssr-spline", str(tmp / "no_spline.txt"), "--image-dir", str(tmp),
                    "--last-run", str(tmp / "no_last_run")]
        try:
            with figures() as kept, contextlib.redirect_stdout(io.StringIO()):
                R.main()
            crash = None
        except Exception as exc:              # noqa: BLE001
            kept, crash = [], repr(exc)
        texts = legend_texts(kept[0]) if kept else []
        check("the ensemble leads when there is no composite: its band is its members' mean error "
              "(0.30, 0.34 -> 0.32 m), the frames its members screened out shaded",
              crash is None and "Camera ensemble (2-h mean; dots: frames)" in texts
              and "tested error per frame, \u00b10.32 m" in texts and M.SCREENED_LABEL in texts,
              crash or texts)


def test_syntax_38():
    print("Python 3.8")
    for f in ("owg_plot_marks.py", "owg_live.py", "owg_report.py", "test_owg_plots.py"):
        try:
            ast.parse((HERE / f).read_text(), filename=f, feature_version=(3, 8))
            good = True
        except SyntaxError as exc:
            good = False
            print(f"    {exc}")
        check(f"{f} parses as Python 3.8", good)


def main():
    test_marks()
    test_live_plot()
    test_report_parts()
    test_drift()
    test_report_end_to_end()
    test_syntax_38()
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        raise SystemExit(1)
    print("all wave-gauge plot tests passed")


if __name__ == "__main__":
    main()

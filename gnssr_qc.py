#!/usr/bin/env python3
"""
Quality Control For GNSS-R Water Levels (QARTOD-style)
=========================================================
Flags GNSS-R water levels that cannot be trusted, so that waterlines and
timestacks at those times get no water level rather than a wrong one.

WHY. GNSS interferometric reflectometry degrades in big waves: the
reflections scatter off crests and foam. In the Sep 2026 storm the
station's GNSS-R read up to ~4 m above the predicted tide on Sep 25-26
(and 2.5-3.0 m NAVD88 on the evening of Sep 26 while the tide at
Chatham was falling), and those levels went straight into the storm DEM
and the runup records.

TESTS (flags follow IOOS QARTOD: 1 good, 2 not evaluated, 3 suspect,
4 fail). A reading's flag is the worst of its tests.

  gross range   Outside GROSS_MIN..GROSS_MAX m NAVD88 -- beyond anything
                tide plus surge produces here.                  -> FAIL
  spike         More than SPIKE_LIMIT from a quadratic through the
                readings within +/-SPIKE_WINDOW_S (at least 2 on each
                side; readings already failed are left out). A curve
                rather than a median, because the tide itself moves up to
                ~0.9 m in two hours. Catches isolated jumps.       -> FAIL
  reference     More than max(REF_LIMIT, REF_SIGMAS * sigma) from the
                level predicted from the NOAA Chatham gauge:
                    GNSS-R ~ a * gauge(t - lag) + b
                with a, lag, b fitted robustly over the last
                REF_FIT_DAYS days (calm Aug-Sep 2026: a 1.29, lag -48 min,
                residual RMS 0.12 m). Catches SUSTAINED failures that a
                spike test lets through.                         -> FAIL
  rate of change  Faster than MAX_RATE m/h between consecutive readings
                -- more than tide plus surge can move.        -> SUSPECT

WAVE SETUP. Chatham (Lydia Cove) is a harbour gauge; the GNSS-R footprint
at Marconi is the surf zone, where breaking waves raise the mean water
level (setup) by an amount the gauge never records. A reading ABOVE the
gauge-predicted level fails the reference test only if it exceeds the
limit PLUS the largest shoreline setup the waves could produce
(Stockdon et al. 2006: 0.35 * beta * sqrt(Hs * L0), L0 = g Tp^2 / 2 pi,
beta = SETUP_BETA), with Hs/Tp from the buoy archive written by
fetch_buoy_waves.py (archive/waves_44008.csv, found next to the gauge
archive). Inside that allowance it is kept as SUSPECT "possible_setup":
in a storm the GNSS-R may be measuring real water the gauge cannot see.
Readings below the gauge are unaffected. No wave archive -> plain limit.

WAVE RUNUP. Above the setup allowance but within the Stockdon 2% runup
R2 = 1.1 * (setup + sqrt(Hs L0 (0.563 beta^2 + 0.004)) / 2), the GNSS-R
is reading the TOTAL water level at the shore (water running up the
beach), not still water. Here such readings still FAIL -- a waterline
labelled with a runup level would sit metres too high -- but with
reason "possible_runup", so they are counted separately and drawn in
their own colour. GPS_code's export_twl.py keeps them as observations
for total water level model validation.

Only FAIL readings are removed by the pipeline; SUSPECT ones are kept
and reported. Where the gauge has no data the reference test is simply
not applied (the other tests still are).

Usage (report and plot):
    python3 gnssr_qc.py /path/to/usgs_spline_out.txt --reference archive/gauge_8447435.csv \\
        --output archive/gnssr_qc.csv --plot gnssr_qc_7day.png

In the pipeline: extract_elevation_contours.py and runup_from_timestack.py
take --gnssr-qc-reference archive/gauge_8447435.csv.
"""

import csv
import argparse
from datetime import datetime, timezone

import numpy as np

GOOD, NOT_EVALUATED, SUSPECT, FAIL = 1, 2, 3, 4

GROSS_MIN, GROSS_MAX = -2.5, 3.5       # m NAVD88
SPIKE_LIMIT = 0.5                      # m
SPIKE_WINDOW_S = 2 * 3600
REF_LIMIT = 0.5                        # m
REF_SIGMAS = 5.0
REF_FIT_DAYS = 30
MAX_RATE = 1.0                         # m per hour
MAX_LAG_MIN = 180
SETUP_COEF = 0.35                      # Stockdon et al. (2006)
SETUP_BETA = 0.10                      # foreshore slope (winter 2025 Marconi DEM)
WAVE_MAX_GAP_S = 3 * 3600


def load_reference(path):
    ep, lv = [], []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if r.get("level_navd88") not in ("", None):
                ep.append(float(r["epoch"]))
                lv.append(float(r["level_navd88"]))
    order = np.argsort(ep)
    return np.asarray(ep)[order], np.asarray(lv)[order]


def load_waves(path):
    """Buoy archive from fetch_buoy_waves.py -> (epochs, hs, tp), or None."""
    try:
        with open(path, newline="") as f:
            rows = [(float(r["epoch"]), float(r["wvht_m"]), float(r["dpd_s"]))
                    for r in csv.DictReader(f) if r.get("wvht_m") and r.get("dpd_s")]
    except (OSError, KeyError, ValueError):
        return None
    if not rows:
        return None
    a = np.array(sorted(rows))
    return a[:, 0], a[:, 1], a[:, 2]


def default_waves_path(reference_path):
    """The buoy archive the cron keeps next to the gauge archive, if there is one."""
    if not reference_path:
        return None
    from pathlib import Path
    p = Path(reference_path).parent / "waves_44008.csv"
    return p if p.exists() else None


def wave_allowances(epochs, waves):
    """Stockdon (2006) shoreline setup and 2% runup (m) at each epoch; 0 where no wave record."""
    epochs = np.asarray(epochs, dtype=float)
    if waves is None:
        return np.zeros(len(epochs)), np.zeros(len(epochs))
    w_ep, hs, tp = waves
    j = np.clip(np.searchsorted(w_ep, epochs), 0, len(w_ep) - 1)
    jm = np.clip(j - 1, 0, len(w_ep) - 1)
    k = np.where(np.abs(w_ep[jm] - epochs) < np.abs(w_ep[j] - epochs), jm, j)
    ok = np.abs(w_ep[k] - epochs) <= WAVE_MAX_GAP_S
    hl = np.maximum(hs[k] * 9.81 * tp[k] ** 2 / (2 * np.pi), 0.0)
    setup = SETUP_COEF * SETUP_BETA * np.sqrt(hl)
    r2 = 1.1 * (setup + np.sqrt(hl * (0.563 * SETUP_BETA ** 2 + 0.004)) / 2)
    return np.where(ok, setup, 0.0), np.where(ok, r2, 0.0)


def setup_allowance(epochs, waves):
    """Largest plausible shoreline wave setup (m) at each epoch; 0 where no wave record."""
    return wave_allowances(epochs, waves)[0]


def reference_at(r_ep, r_lv, epochs):
    """Gauge level at `epochs`; NaN outside the record or across gaps > 12 min."""
    out = np.interp(epochs, r_ep, r_lv, left=np.nan, right=np.nan)
    if len(r_ep) > 1:
        i = np.clip(np.searchsorted(r_ep, epochs), 1, len(r_ep) - 1)
        out[(r_ep[i] - r_ep[i - 1]) > 720] = np.nan
    return out


def fit_reference(ep, lv, r_ep, r_lv, keep):
    """
    Robust fit of GNSS-R = a * gauge(t - lag) + b over the last
    REF_FIT_DAYS of readings not already failed. Returns
    (a, lag_s, b, sigma) or None if there is too little overlap.
    """
    recent = keep & (ep >= ep.max() - REF_FIT_DAYS * 86400)
    best = None
    for lag_min in range(-MAX_LAG_MIN, MAX_LAG_MIN + 1, 6):
        x = reference_at(r_ep, r_lv, ep - lag_min * 60)
        m = recent & np.isfinite(x)
        if m.sum() < 48:
            continue
        use = m.copy()
        for _ in range(4):                          # iteratively drop outliers
            A = np.c_[x[use], np.ones(use.sum())]
            coef, *_ = np.linalg.lstsq(A, lv[use], rcond=None)
            res = lv - (coef[0] * x + coef[1])
            sigma = max(1.4826 * np.nanmedian(np.abs(res[use])), 0.03)
            use = m & (np.abs(res) < max(3 * sigma, 0.15))
        rms = float(np.sqrt(np.mean(res[use] ** 2)))
        if best is None or rms < best[4]:
            best = (float(coef[0]), lag_min * 60.0, float(coef[1]), float(sigma), rms)
    return None if best is None else best[:4]


def run_qc(ep, lv, reference=None, waves=None):
    """
    Returns (flags, reasons, predicted, fit). `reference` is (r_ep, r_lv)
    or None. `predicted` is the reference-based level (NaN if none).
    `waves` is (epochs, hs, tp) or None -- see WAVE SETUP.
    """
    n = len(lv)
    flags = np.full(n, GOOD, dtype=np.int8)
    reasons = [[] for _ in range(n)]

    def mark(i, flag, why):
        flags[i] = max(flags[i], flag)
        reasons[i].append(why)

    for i in np.flatnonzero((lv < GROSS_MIN) | (lv > GROSS_MAX)):
        mark(i, FAIL, "gross_range")

    # Reference test first, so the spike test below does not judge good
    # readings against neighbours that are themselves failures.
    predicted = np.full(n, np.nan)
    fit = None
    if reference is not None and len(reference[0]) > 1:
        fit = fit_reference(ep, lv, reference[0], reference[1], flags < FAIL)
        if fit is not None:
            a, lag_s, b, sigma = fit
            predicted = a * reference_at(reference[0], reference[1], ep - lag_s) + b
            limit = max(REF_LIMIT, REF_SIGMAS * sigma)
            allowance, runup = wave_allowances(ep, waves)
            for i in np.flatnonzero(np.isfinite(predicted) & (np.abs(lv - predicted) > limit)):
                excess = lv[i] - predicted[i]
                if 0 < excess <= limit + allowance[i]:
                    mark(i, SUSPECT, "possible_setup")
                elif 0 < excess <= limit + runup[i]:
                    mark(i, FAIL, "possible_runup")    # total water level, not still water
                else:
                    mark(i, FAIL, "reference")

    ok = flags < FAIL
    for i in range(n):
        if not ok[i]:
            continue
        near = ok & (np.abs(ep - ep[i]) <= SPIKE_WINDOW_S)
        near[i] = False
        before, after = near & (ep < ep[i]), near & (ep > ep[i])
        if before.sum() < 2 or after.sum() < 2:
            continue                                   # not evaluated at edges and gaps
        tt = (ep[near] - ep[i]) / 3600.0
        coef = np.polyfit(tt, lv[near], 2)
        if abs(lv[i] - coef[2]) > SPIKE_LIMIT:         # coef[2] = curve value at t_i
            mark(i, FAIL, "spike")

    for i in range(1, n):
        dt_h = (ep[i] - ep[i - 1]) / 3600.0
        if 0 < dt_h <= 2 and abs(lv[i] - lv[i - 1]) / dt_h > MAX_RATE:
            mark(i, SUSPECT, "rate_of_change")
    return flags, [";".join(r) for r in reasons], predicted, fit


def qc_filter(ep, lv, reference_path, waves_path=None):
    """For the pipeline: drop FAIL readings. Returns (ep, lv, summary text)."""
    try:
        reference = load_reference(reference_path)
    except OSError:
        reference = None
    waves = load_waves(waves_path or default_waves_path(reference_path) or "")
    flags, reasons, _, fit = run_qc(ep, lv, reference, waves)
    keep = flags < FAIL
    counts = {}
    for r, f in zip(reasons, flags):
        if f >= FAIL:
            for why in r.split(";"):
                counts[why] = counts.get(why, 0) + 1
    text = (f"GNSS-R QC: {int((~keep).sum())} of {len(lv)} reading(s) failed and removed"
            + (f" ({', '.join(f'{k} {v}' for k, v in sorted(counts.items()))})" if counts else "")
            + f"; {int((flags == SUSPECT).sum())} suspect kept"
            + f" ({sum('possible_setup' in r for r in reasons)} within wave setup of the gauge)"
            + ("" if fit else "; reference test NOT applied (no gauge overlap)"))
    return ep[keep], lv[keep], text


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("gnssr", help="gnssrefl spline file")
    ap.add_argument("--reference", help="Gauge archive from fetch_tide_gauge.py")
    ap.add_argument("--waves", help="Buoy archive from fetch_buoy_waves.py (default: waves_44008.csv "
                                    "next to --reference, if present)")
    ap.add_argument("--output", help="Write every reading with its flag and reasons (CSV).")
    ap.add_argument("--plot", help="Figure of the last --plot-days with failed readings marked.")
    ap.add_argument("--plot-days", type=int, default=7)
    args = ap.parse_args()

    from extract_elevation_contours import load_gnssr_spline
    ep, lv, _, _ = load_gnssr_spline(args.gnssr)
    reference = load_reference(args.reference) if args.reference else None
    waves_path = args.waves or default_waves_path(args.reference)
    waves = load_waves(waves_path) if waves_path else None
    print(f"Wave setup allowance: {'from ' + str(waves_path) if waves else 'OFF (no wave archive)'}")
    flags, reasons, predicted, fit = run_qc(ep, lv, reference, waves)

    print()
    print(f"GNSS-R QC over {len(lv)} readings "
          f"({datetime.fromtimestamp(ep[0], tz=timezone.utc):%Y-%m-%d} to "
          f"{datetime.fromtimestamp(ep[-1], tz=timezone.utc):%Y-%m-%d})")
    if fit:
        a, lag_s, b, sigma = fit
        print(f"   reference fit (last {REF_FIT_DAYS} d): a {a:.3f}, lag {lag_s / 60:+.0f} min, "
              f"b {b:+.3f} m, sigma {sigma:.3f} m -> limit {max(REF_LIMIT, REF_SIGMAS * sigma):.2f} m")
    else:
        print("   reference test not applied (no gauge record overlapping)")
    for name, f in (("good", GOOD), ("suspect", SUSPECT), ("fail", FAIL)):
        print(f"   {name:8s}: {int((flags == f).sum())}")
    fail_days = {}
    for e, f in zip(ep, flags):
        if f >= FAIL:
            d = datetime.fromtimestamp(e, tz=timezone.utc).strftime("%Y-%m-%d")
            fail_days[d] = fail_days.get(d, 0) + 1
    if fail_days:
        print("   failed readings by day: " + ", ".join(f"{d} {n}" for d, n in sorted(fail_days.items())))

    if args.output:
        with open(args.output, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["time_utc", "epoch", "level_navd88", "predicted_navd88", "residual_m",
                        "qc_flag", "qc_reasons"])
            for e, l, p, fl, r in zip(ep, lv, predicted, flags, reasons):
                w.writerow([datetime.fromtimestamp(e, tz=timezone.utc).isoformat(), int(e),
                            round(float(l), 3), "" if not np.isfinite(p) else round(float(p), 3),
                            "" if not np.isfinite(p) else round(float(l - p), 3), int(fl), r])
        print(f"wrote {args.output}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        sel = ep >= ep.max() - args.plot_days * 86400
        t = np.array([datetime.fromtimestamp(e, tz=timezone.utc) for e in ep[sel]])
        fig, ax = plt.subplots(figsize=(13, 4.5), dpi=110)
        if fit and reference is not None:
            # Expected level on a 6-minute grid, so it follows the gauge
            # through GNSS-R gaps instead of joining them with straight
            # lines; NaN where the gauge itself has a gap.
            a, lag_s, b, _ = fit
            grid = np.arange(ep[sel][0], ep[sel][-1] + 1, 360.0)
            exp_lv = a * reference_at(reference[0], reference[1], grid - lag_s) + b
            ax.plot([datetime.fromtimestamp(e, tz=timezone.utc) for e in grid], exp_lv,
                    color="#e0a060", lw=1.5, label="expected from Chatham gauge")
        good = flags[sel] < FAIL
        ax.plot(t[good], lv[sel][good], ".", ms=4, color="#1f77b4", label="GNSS-R, passed QC")
        if (~good).any():
            bad = ~good & ~np.array(["possible_runup" in r for r in
                                     np.asarray(reasons, dtype=object)[sel]], bool)
            ax.plot(t[bad], lv[sel][bad], "x", ms=6, color="#c0392b", label="GNSS-R, failed QC (not used)")
        setup = np.array(["possible_setup" in r for r in np.asarray(reasons, dtype=object)[sel]], bool)
        sus = (flags[sel] == SUSPECT) & ~setup
        if sus.any():
            ax.plot(t[sus], lv[sel][sus], "o", ms=7, mfc="none", color="#8e44ad", label="suspect (kept)")
        if setup.any():
            ax.plot(t[setup], lv[sel][setup], "D", ms=6, mfc="none", color="#e67e22",
                    label="above gauge by no more than wave setup (kept)")
        runup_m = np.array(["possible_runup" in r for r in np.asarray(reasons, dtype=object)[sel]], bool)
        if runup_m.any():
            ax.plot(t[runup_m], lv[sel][runup_m], "s", ms=6, mfc="none", color="#7d3c98",
                    label="total water level, within wave runup (not used)")
        ax.set_ylabel("water level (m NAVD88)")
        ax.grid(alpha=0.3)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax.legend(loc="upper left", fontsize=8)
        ax.set_title(f"GNSS-R water level quality control, last {args.plot_days} days (UTC)")
        fig.tight_layout()
        fig.savefig(args.plot)
        print(f"wrote {args.plot}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Wave And Water-Level Validation Against The Marconi ADCP
===========================================================
The Signature 1000 ADCP (sig1000_waves_ALL.csv, Dec 2024 - Mar 2025,
~21 m depth off Marconi) measured waves and water level AT the site.
That makes it ground truth for three things the waterline product
depends on:

A. OFFSHORE WAVE SOURCES. Operationally there is no ADCP: wave height
   comes from a buoy (NDBC 44008, fetch_buoy_waves.py) or a model, and is
   used to (1) drop frames in rough seas from the DEM (--max-hs) and
   (2) tag runup records. For each source this reports correlation,
   bias, RMSE, the best time lag, the ratio that scales it to Marconi,
   and -- what the DEM filter actually needs -- how often it agrees with
   the ADCP on "Hs above the threshold". Sources:
     --wis FILE     USACE Wave Information Study hindcast CSV
                    (datetime, waveHs, waveTp, waveMeanDirection)
     --ndbc FILE    NDBC CSV (datetime, WVHT, DPD, APD, MWD; 99 = missing)
     --ndbc-station ID   download NDBC historical stdmet for the ADCP
                    period (on a machine that can reach ndbc.noaa.gov),
                    e.g. 44008, the station the station cron uses.

B. WAVE BIAS IN THE WATERLINE (needs --contours, the georectified
   contour_points_ground.csv with offshore_hs_m from the ADCP). Waves
   push the timex waterline up the beach (setup + swash), so in rough
   seas the line sits where the beach is HIGHER than the water level it
   is labelled with. Each frame is compared with the elevation the OTHER
   days put at the same 2 m cells (leave-one-day-out), and the frame
   offsets are related to Hs and to sqrt(Hs * L0) (the Stockdon setup
   scale). Result: the size of the bias and the Hs above which it
   exceeds 0.10 m -- a measured --max-hs instead of an assumed one.

C. ADCP CLOCK (needs --tide-model, a tide spreadsheet with 'time' and
   '*_heightm' columns). adcp_to_navd88.py found the ADCP best matched
   Chatham ~25 min later than the GNSS-R did. The tide model is
   independent of Chatham: the lag of the ADCP water level against it
   shows whether the ADCP timestamps are shifted.

Usage:
    python3 validate_waves.py --adcp sig1000_waves_ALL.csv \\
        --wis WIS_ST63064.csv --ndbc NDBC_44013.csv --ndbc-station 44008 \\
        --contours /mnt/I2Rgus_Data/Chelsea_calibration/contour_points_ground.csv \\
        --tide-model marconi_tides_sherwood.xlsx \\
        --output /mnt/I2Rgus_Data/Chelsea_calibration/validation
"""

import io
import gzip
import argparse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

HS_BINS = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 10]
FILTER_HS = 1.5          # the --max-hs the station cron uses
BIAS_LIMIT = 0.10        # m; a frame bias worth filtering out
CELL = 2.0               # m, same as the DEM


def load_adcp(path):
    a = pd.read_csv(path)
    a["time"] = pd.to_datetime(a["time"])
    a = a.set_index("time").sort_index()
    return a.rename(columns={"wh_4061": "hs", "wp_peak": "tp", "wvdir": "dir"})


def load_wis(path):
    w = pd.read_csv(path, parse_dates=["datetime"]).set_index("datetime").sort_index()
    return w.rename(columns={"waveHs": "hs", "waveTp": "tp", "waveMeanDirection": "dir"})[["hs", "tp", "dir"]]


def load_ndbc_csv(path):
    n = pd.read_csv(path, parse_dates=["datetime"]).set_index("datetime").sort_index()
    n = n.rename(columns={"WVHT": "hs", "DPD": "tp", "MWD": "dir"})
    for c, bad in (("hs", 99), ("tp", 99), ("dir", 999)):
        if c in n:
            n[c] = pd.to_numeric(n[c], errors="coerce")
            n.loc[n[c] >= bad, c] = np.nan
    return n[["hs", "tp", "dir"]].resample("1h").mean()


def fetch_ndbc(station, years):
    """NDBC historical standard meteorological files, parsed to hourly hs/tp/dir."""
    frames = []
    for y in years:
        url = f"https://www.ndbc.noaa.gov/data/historical/stdmet/{station}h{y}.txt.gz"
        raw = gzip.decompress(urllib.request.urlopen(url, timeout=60).read()).decode()
        lines = [l for l in raw.splitlines() if l.strip()]
        head = lines[0].lstrip("#").split()
        body = [l.split() for l in lines[1:] if not l.startswith("#")]
        d = pd.DataFrame(body, columns=head).apply(pd.to_numeric, errors="coerce")
        yy = "YY" if "YY" in d else "YYYY"
        d.index = pd.to_datetime(dict(year=d[yy], month=d["MM"], day=d["DD"],
                                      hour=d["hh"], minute=d.get("mm", 0)))
        frames.append(d)
    d = pd.concat(frames).sort_index()
    out = pd.DataFrame({"hs": d["WVHT"], "tp": d["DPD"], "dir": d["MWD"]})
    out.loc[out.hs >= 99, "hs"] = np.nan
    out.loc[out.tp >= 99, "tp"] = np.nan
    out.loc[out["dir"] >= 999, "dir"] = np.nan
    return out.resample("1h").mean()


def compare(truth, src, name, lines):
    """Part A statistics of one source against the ADCP."""
    best = None
    for lag in range(-6, 7):
        j = pd.concat([truth, src.shift(lag, freq="1h")], axis=1, join="inner").dropna()
        if len(j) > 48:
            r = np.corrcoef(j.iloc[:, 0], j.iloc[:, 1])[0, 1]
            if best is None or r > best[0]:
                best = (r, lag)
    j = pd.concat([truth.rename("adcp"), src.rename("src")], axis=1, join="inner").dropna()
    if len(j) < 48:
        have = src.dropna()
        span = (f"it has wave data {have.index.min():%Y-%m-%d} to {have.index.max():%Y-%m-%d}"
                if len(have) else "it reported no wave heights at all in the files read")
        lines.append(f"{name}: only {len(j)} hours with a wave height during the ADCP record -- "
                     f"not compared ({span}; the buoy was probably off station or not reporting "
                     f"waves). Try a neighbouring buoy, e.g. --ndbc-station 44018 or 44020.")
        return None
    x, y = j.src, j.adcp
    ratio = float(np.sum(x * y) / np.sum(x * x))
    hit = ((x > FILTER_HS) & (y > FILTER_HS)).sum()
    rough = (y > FILTER_HS).sum()
    false = ((x > FILTER_HS) & (y <= FILTER_HS)).sum()
    calm = (y <= FILTER_HS).sum()
    lines += [
        f"{name}: {len(j)} overlapping hours, {j.index.min():%Y-%m-%d} to {j.index.max():%Y-%m-%d}",
        f"   r {np.corrcoef(x, y)[0, 1]:.2f}   bias (source - ADCP) {np.mean(x - y):+.2f} m   "
        f"RMSE {np.sqrt(np.mean((x - y) ** 2)):.2f} m   best lag {best[1]:+d} h (r {best[0]:.2f})",
        f"   Marconi Hs ~ {ratio:.2f} x source Hs",
        f"   Hs > {FILTER_HS} m filter: catches {hit}/{rough} rough ADCP hours "
        f"({100 * hit / max(rough, 1):.0f}%), false alarms {false}/{calm} calm hours "
        f"({100 * false / max(calm, 1):.0f}%)",
    ]
    # The cut-off ON THIS SOURCE that flags 90% of the hours rough at
    # Marconi -- what --max-hs should be when this source drives it.
    for thr in np.arange(0.5, 3.01, 0.05):
        caught = ((x > thr) & (y > FILTER_HS)).sum() / max(rough, 1)
        if caught < 0.90:
            thr -= 0.05
            fa = ((x > thr) & (y <= FILTER_HS)).sum() / max(calm, 1)
            lines.append(f"   to catch 90% of Marconi Hs > {FILTER_HS} m hours with this source: "
                         f"--max-hs {thr:.2f} (false alarms {100 * fa:.0f}% of calm hours)")
            break
    by = j.groupby(pd.cut(j.adcp, HS_BINS)).apply(lambda g: pd.Series(
        {"n": len(g), "bias": (g.src - g.adcp).mean()}))
    lines.append("   bias by ADCP Hs: " + ", ".join(
        f"{iv.left:g}-{iv.right:g} m {r.bias:+.2f} (n={int(r.n)})" for iv, r in by.iterrows() if r.n))
    return j


def wave_bias(contours, lines):
    """Part B: frame elevation offset against the other days, versus Hs."""
    g = pd.read_csv(contours)
    if "offshore_hs_m" not in g:
        lines.append("B. contours have no offshore_hs_m column -- rerun process_chelsea.py (it now tags frames)")
        return None
    g = g.dropna(subset=["easting_utm19", "northing_utm19", "tide_elevation_navd88"])
    g["day"] = g.capture_time_utc.str[:10]
    g["cell"] = (np.floor(g.easting_utm19 / CELL).astype(int).astype(str) + "_" +
                 np.floor(g.northing_utm19 / CELL).astype(int).astype(str))
    z = "beach_elevation_navd88" if "beach_elevation_navd88" in g else "tide_elevation_navd88"
    per_frame = []
    for day, gd in g.groupby("day"):
        others = g[g.day != day].groupby("cell")[z].median()
        for src, gf in gd.groupby("source_file"):
            ref = gf.cell.map(others)
            ok = ref.notna()
            if ok.sum() < 20:
                continue
            per_frame.append({"frame": src, "day": day, "n": int(ok.sum()),
                              "offset": float(np.median(gf[z][ok] - ref[ok])),
                              "hs": gf.offshore_hs_m.iloc[0], "tp": gf.offshore_tp_s.iloc[0]})
    f = pd.DataFrame(per_frame).dropna(subset=["hs"])
    if len(f) < 8:
        lines.append(f"B. only {len(f)} frames overlap other days with a wave record -- too few")
        return None
    f["x"] = np.sqrt(f.hs * 9.81 * f.tp.fillna(f.tp.median()) ** 2 / (2 * np.pi))
    r_hs = np.corrcoef(f.hs, f.offset)[0, 1]
    r_x = np.corrcoef(f.x, f.offset)[0, 1]
    slope, icpt = np.polyfit(f.hs, f.offset, 1)
    lines += [f"B. {len(f)} frames: elevation offset against the other days' DEM (leave-one-day-out)",
              f"   offset vs Hs: r {r_hs:+.2f}, {slope:+.3f} m per m of Hs (intercept {icpt:+.3f} m)",
              f"   offset vs sqrt(Hs*L0) (setup scale): r {r_x:+.2f}",
              "   (negative = rough-sea lines sit landward, i.e. labelled lower than the beach they cross)"]
    by = f.groupby(pd.cut(f.hs, HS_BINS))["offset"].agg(["count", "median"])
    lines.append("   median offset by Hs: " + ", ".join(
        f"{iv.left:g}-{iv.right:g} m {r['median']:+.2f} (n={int(r['count'])})"
        for iv, r in by.iterrows() if r["count"]))
    # A wave effect is a CHANGE with Hs: compare rougher bins with the
    # calm-sea baseline (Hs < 1 m), whose own offset is frame-to-frame
    # noise (detector scatter, beach change between days), not waves.
    calm = f[f.hs < 1.0]["offset"]
    spread = float(np.subtract(*np.percentile(f.offset, [84, 16])) / 2)
    lines.append(f"   frame-to-frame scatter (half 16-84% range): {spread:.2f} m")
    if len(calm) < 5:
        lines.append("   => too few calm frames (Hs < 1 m) for a baseline; no recommendation")
        return f
    base = float(calm.median())
    worse = [iv.left for iv, r in by.iterrows()
             if iv.left >= 1.0 and r["count"] >= 3 and abs(r["median"] - base) > BIAS_LIMIT]
    rough_n = int((f.hs >= 1.0).sum())
    if worse:
        lines.append(f"   => recommended --max-hs: {min(worse):g} m (first Hs bin whose median offset "
                     f"differs from the calm baseline {base:+.2f} m by more than {BIAS_LIMIT} m)")
    else:
        lines.append(f"   => no wave dependence detected: rougher bins stay within {BIAS_LIMIT} m of "
                     f"the calm baseline ({base:+.2f} m); {rough_n} frame(s) with Hs >= 1 m"
                     + (" -- too few to be sure" if rough_n < 15 else ""))
    return f


def adcp_clock(adcp, tide_path, lines):
    """Part C: lag of the ADCP water level against the tide model."""
    t = pd.read_excel(tide_path) if str(tide_path).endswith("xlsx") else pd.read_csv(tide_path)
    cols = [c for c in t.columns if str(c).endswith("_heightm")]
    t["time"] = pd.to_datetime(t["time"])
    model = t.set_index("time")[cols].mean(axis=1).resample("6min").mean().interpolate(limit=10)
    wl = adcp["water_level"].dropna()
    best = None
    for lag_min in range(-120, 121, 3):
        m = np.interp(wl.index.astype("int64") / 1e9 - lag_min * 60,
                      model.index.astype("int64") / 1e9, model.values, left=np.nan, right=np.nan)
        ok = np.isfinite(m)
        if ok.sum() < 200:
            continue
        r = np.corrcoef(m[ok], wl.values[ok])[0, 1]
        if best is None or r > best[1]:
            best = (lag_min, r, int(ok.sum()))
    if best is None:
        lines += ["C. the tide model does not cover the ADCP period (Dec 2024 - Mar 2025).",
                  "   Generate the same model output for those months and pass it with --tide-model."]
        return
    lines += [f"C. ADCP water level vs tide model: best lag {best[0]:+d} min (r {best[1]:.3f}, {best[2]} h)",
              "   A lag near 0 means the ADCP timestamps are right (the model is for Marconi itself);",
              "   a lag of 15-30 min means they mark the start or end of an averaging burst, and",
              "   water levels at mid-tide are off by up to ~0.2 m -- shift them with --time-shift."]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--adcp", default="sig1000_waves_ALL.csv")
    ap.add_argument("--wis", action="append", default=[])
    ap.add_argument("--ndbc", action="append", default=[])
    ap.add_argument("--ndbc-station", action="append", default=[])
    ap.add_argument("--contours", help="contour_points_ground.csv with offshore_hs_m (part B)")
    ap.add_argument("--tide-model", help="tide spreadsheet with *_heightm columns (part C)")
    ap.add_argument("--output", default="validation")
    args = ap.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    adcp = load_adcp(args.adcp)
    lines = [f"ADCP truth: {adcp.index.min():%Y-%m-%d %H:%M} to {adcp.index.max():%Y-%m-%d %H:%M}, "
             f"Hs mean {adcp.hs.mean():.2f} m, max {adcp.hs.max():.2f} m, "
             f"{100 * (adcp.hs > FILTER_HS).mean():.0f}% of hours above {FILTER_HS} m", "",
             "A. Offshore wave sources vs the ADCP at Marconi"]

    sources = {}
    for p in args.wis:
        sources[f"WIS {Path(p).stem}"] = load_wis(p)
    for p in args.ndbc:
        sources[f"NDBC {Path(p).stem}"] = load_ndbc_csv(p)
    for st in args.ndbc_station:
        try:
            years = sorted({adcp.index.min().year, adcp.index.max().year})
            sources[f"NDBC {st} (historical)"] = fetch_ndbc(st, years)
        except Exception as exc:
            lines.append(f"NDBC {st}: could not download ({exc})")
    joined = {}
    for name, s in sources.items():
        j = compare(adcp.hs, s.hs, name, lines)
        if j is not None:
            joined[name] = j
        lines.append("")

    frames = wave_bias(args.contours, lines) if args.contours else None
    if args.tide_model:
        lines.append("")
        adcp_clock(adcp, args.tide_model, lines)

    report = "\n".join(lines)
    (out / "validation_report.txt").write_text(report + "\n")
    if frames is not None:
        frames.to_csv(out / "frame_wave_bias.csv", index=False)
    print(report)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = 2 + (frames is not None)
    fig, ax = plt.subplots(n, 1, figsize=(12, 3.6 * n))
    ax[0].plot(adcp.index, adcp.hs, "k", lw=1.4, label="ADCP at Marconi (truth)")
    for name, s in sources.items():
        s = s[(s.index >= adcp.index.min()) & (s.index <= adcp.index.max())]
        ax[0].plot(s.index, s.hs, lw=0.9, alpha=0.8, label=name)
    ax[0].axhline(FILTER_HS, color="grey", ls="--", lw=0.8)
    ax[0].set_ylabel("Hs (m)")
    ax[0].legend(fontsize=8)
    ax[0].set_title("Significant wave height")
    for name, j in joined.items():
        ax[1].plot(j.src, j.adcp, ".", ms=2, alpha=0.4, label=name)
    lim = [0, max(4, float(adcp.hs.max()) + 0.2)]
    ax[1].plot(lim, lim, "k:", lw=0.8)
    ax[1].set_xlim(lim)
    ax[1].set_ylim(lim)
    ax[1].set_xlabel("source Hs (m)")
    ax[1].set_ylabel("ADCP Hs (m)")
    if joined:
        ax[1].legend(fontsize=8)
    if frames is not None:
        ax[2].plot(frames.hs, frames.offset, "o", ms=4)
        ax[2].axhline(0, color="k", lw=0.6)
        ax[2].axhspan(-BIAS_LIMIT, BIAS_LIMIT, color="green", alpha=0.08)
        ax[2].set_xlabel("ADCP Hs at the frame (m)")
        ax[2].set_ylabel("frame offset vs other days (m)")
        ax[2].set_title("Wave bias of the detected waterline (one point per frame)")
    fig.tight_layout()
    fig.savefig(out / "validation.png", dpi=120)
    print(f"\nWrote {out / 'validation_report.txt'} and {out / 'validation.png'}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Tests For Partial Waterlines (Detector v6.11)
===============================================
Plain python3 (no pytest needed on the station):

    python3 test_partial_lines.py

WHAT.
  * clear_runs: long runs of signal columns kept, short ones and scattered
    columns not; dropouts up to PARTIAL_MAX_GAP_COLUMNS bridged but still
    no-signal; runs cut where they meet noise, not at the frame's edge.
  * process_image on a synthetic energy field (the real smoothing, solver,
    signal test, confidence gate and CSV export; only the photo and its
    energy are made up):
      - a frame with signal everywhere: saved as before, the CSV the same to
        the byte with partial lines on and off;
      - the far half hazy: discarded with partial lines off, saved as a
        partial line with them on -- only the clear half flagged Has_Signal,
        its rows on the true waterline, the hazy half dropped, not the next
        frame's temporal-blend reference;
      - only scattered short stretches clear: discarded either way, saying why.
  * extract_elevation_contours.py: a partial line is judged by
    --partial-min-coverage (0.25), every other frame by --min-coverage as
    before; 0 turns the partial floor off.
  * the detector and the extractor agree on the 60% that makes a line partial.
Exit status 0 when every check passes.
"""

import ast
import contextlib
import csv
import io
import subprocess
import sys
import tempfile
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import waterline_detector_v5 as D          # noqa: E402
import extract_elevation_contours as X     # noqa: E402

failures = []
EPOCH = 1792886400                          # 2026-10-25 00:00 UTC
ROWS, COLS = 400, 600


def check(name, cond, detail=""):
    print(f"  {'ok ' if cond else 'FAIL'} {name}" + (f"  ({detail})" if detail != "" and not cond else ""))
    if not cond:
        failures.append(name)


def test_clear_runs():
    print("clear_runs")
    n = 200
    m, k = D.clear_runs(np.ones(n, bool), min_run=20, max_gap=8, edge_trim=12)
    check("all signal: every column, no trim at the frame's edges", m.all() and k == 1)
    s = np.zeros(n, bool)
    s[60:140] = True
    m, k = D.clear_runs(s, min_run=20, max_gap=8, edge_trim=12)
    check("a run between noise: cut by 12 columns at both ends", k == 1 and m.sum() == 80 - 24
          and m[72] and not m[71] and m[127] and not m[128])
    s = np.zeros(n, bool)
    s[:90] = True
    m, _ = D.clear_runs(s, min_run=20, max_gap=8, edge_trim=12)
    check("a run from the frame's edge: cut only where it meets noise", m[0] and m[77] and not m[78])
    s = np.zeros(n, bool)
    s[20:180] = True
    s[100:106] = False
    m, k = D.clear_runs(s, min_run=20, max_gap=8, edge_trim=12)
    check("a 6-column dropout bridged: one run, the dropout itself still no-signal",
          k == 1 and m[99] and not m[100:106].any() and m[106] and m.sum() == 160 - 6 - 24)
    s = np.zeros(n, bool)
    s[10:60] = True
    s[80:95] = True
    m, k = D.clear_runs(s, min_run=20, max_gap=8, edge_trim=12)
    check("a gap over 8 columns splits; a 15-column run under min_run is left out",
          k == 1 and m[30] and not m[80:95].any())
    s = np.zeros(n, bool)
    s[::7] = True
    m, k = D.clear_runs(s, min_run=20, max_gap=8, edge_trim=12)
    check("scattered single columns (a comb whose gaps are short enough to bridge): nothing, "
          "under PARTIAL_MIN_RUN_DENSITY", not m.any() and k == 0)
    s = np.ones(n, bool)
    s[::4] = False                                    # 75% dense: under 0.80
    m, k = D.clear_runs(s, min_run=20, max_gap=8, edge_trim=12)
    s2 = np.ones(n, bool)
    s2[3::6] = False                                  # 83% dense, clear at both frame edges
    m2, k2 = D.clear_runs(s2, min_run=20, max_gap=8, edge_trim=12)
    check("a run 75% clear left out, one 83% clear kept", not m.any() and k2 == 1 and m2.sum() == s2.sum())
    m, k = D.clear_runs(np.zeros(n, bool))
    check("no signal: nothing", not m.any() and k == 0)
    m, _ = D.clear_runs(s)
    check("the defaults: PARTIAL_MIN_RUN_FRACTION of the columns", not m.any())


def truth_row(c):
    return 150.0 + 0.15 * c


def field(clear):
    """An energy field: a ridge on truth_row where `clear`, flat random haze elsewhere."""
    rng = np.random.default_rng(3)
    e = 0.02 + 0.01 * rng.random((ROWS, COLS))
    r = np.arange(ROWS)[:, None]
    ridge = np.exp(-0.5 * ((r - truth_row(np.arange(COLS))[None, :]) / 3.0) ** 2)
    haze = 0.25 + 0.05 * rng.random((ROWS, COLS))
    return np.where(clear[None, :], e + ridge, haze).astype(np.float32)


@contextlib.contextmanager
def detector(tmp, energy, partial_min):
    """process_image with a made-up photo and energy field, everything else the detector's own."""
    profile = replace(D.CAMERAS["CACO05_C2"], bias_correction_points=(), bias_correction_points_timex=(),
                      band_contrast_floor=0.0)
    saved = {k: getattr(D, k) for k in ("prepare_image", "build_shoreline_energy", "save_energy_debug",
                                        "OUTPUT_FOLDER", "EXPORT_CSV", "PARTIAL_MIN_FRACTION",
                                        "TEMPORAL_FILTER_ENABLED")}

    def prepare(path):
        d = D.ImageData()
        d.filename = str(path)
        d.original = np.full((ROWS + 100, COLS, 3), 128, np.uint8)
        d.cropped = np.random.default_rng(1).integers(60, 200, (ROWS, COLS, 3)).astype(np.uint8)
        d.crop_top, d.crop_left, d.profile = 50, 0, profile
        return d

    D.prepare_image = prepare
    D.build_shoreline_energy = lambda data: (energy, {"envelope": np.ones((ROWS, COLS), bool)})
    D.save_energy_debug = lambda *a, **k: None
    D.OUTPUT_FOLDER, D.EXPORT_CSV, D.PARTIAL_MIN_FRACTION = str(tmp), True, partial_min
    D.TEMPORAL_FILTER_ENABLED = False
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(D, k, v)


def run_frame(tmp, clear, partial_min, epoch=EPOCH):
    tmp = Path(tmp)
    tmp.mkdir(parents=True, exist_ok=True)
    name = f"{epoch}.Sun.Oct.25_00_00_00.GMT.2026.CACO05.c2.timex.jpg"
    D.STATE.previous_shoreline.pop("CACO05_C2", None)
    buf = io.StringIO()
    with detector(tmp, field(clear), partial_min), contextlib.redirect_stdout(buf):
        _, status = D.process_image(tmp / name)
    csv_path = tmp / (Path(name).stem + ".csv")
    rows = list(csv.DictReader(open(csv_path))) if csv_path.exists() else []
    return status, rows, buf.getvalue(), csv_path


def test_detector():
    print("waterline_detector_v5.py: partial lines")
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        full = np.ones(COLS, bool)
        st_on, rows_on, _, p_on = run_frame(tmp / "on", full, 0.30)
        st_off, rows_off, _, p_off = run_frame(tmp / "off", full, 0.0)
        check("signal everywhere: saved, the same CSV to the byte with partial lines on and off",
              st_on == st_off == "saved" and p_on.exists() and p_on.read_bytes() == p_off.read_bytes())

        hazy = np.zeros(COLS, bool)
        hazy[COLS // 2:] = True                       # the far half hidden
        st, rows, out, _ = run_frame(tmp / "h_off", hazy, 0.0)
        check("the far half hazy, partial lines off: discarded as before",
              st == "no_signal" and not rows and "DISCARDED (only" in out, (st, out[-300:]))
        before = D.STATE.previous_shoreline.get("CACO05_C2")
        st, rows, out, _ = run_frame(tmp / "h_on", hazy, 0.30)
        sig = np.array([int(r["Has_Signal"]) for r in rows], bool) if rows else np.zeros(0, bool)
        col = np.array([float(r["Column"]) for r in rows]) if rows else np.zeros(0)
        row = np.array([float(r["Row_Precise"]) for r in rows]) if rows else np.zeros(0)
        check("partial lines on: saved as a partial line, logged, all columns written",
              st == "partial" and len(rows) == COLS and "PARTIAL LINE" in out and "Finished." in out,
              (st, out[-400:]))
        if len(sig):
            check("only clear columns flagged, the hazy half none, at least 30% kept",
                  not sig[:COLS // 2].any() and sig.mean() >= 0.30, (sig[:COLS // 2].sum(), sig.mean()))
            err = row[sig] - (truth_row(col[sig]) + 50)
            check("the kept rows on the true waterline (within 1.5 px)", np.abs(err).max() <= 1.5,
                  np.abs(err).max() if len(err) else "none")
            first = int(np.flatnonzero(sig)[0]) if sig.any() else -1
            check("cut where the clear run meets the haze", first >= COLS // 2 + D.PARTIAL_EDGE_TRIM_COLUMNS - 2, first)
        check("a partial line is not the next frame's temporal-blend reference",
              D.STATE.previous_shoreline.get("CACO05_C2") is before)

        bits = np.zeros(COLS, bool)
        for a in range(0, COLS, 100):
            bits[a:a + 40] = True                     # 40% of columns, in 40-column stretches
        st, rows, out, _ = run_frame(tmp / "bits", bits, 0.30)
        check("only scattered short stretches clear: discarded, saying the long runs fall short",
              st == "no_signal" and not rows and "long clear runs cover" in out, (st, out[-300:]))

    lines = Path(D.__file__).read_text()
    check("v6.11 named in the detector's header", "Version 6.11" in lines and "PARTIAL LINES" in lines)


def write_frame(folder, ep, signal):
    t = datetime.fromtimestamp(ep, tz=timezone.utc)
    out = ["Column,Row,Confidence,Row_Precise,Peak_Sharpness,Has_Signal"]
    for c in range(300):
        row = 500 + c / 10.0 + (c * 7 % 13) / 4.0      # not dead straight (the extractor drops those)
        out.append(f"{c},{int(row)},0.95,{row:.3f},0.5,{int(bool(signal[c]))}")
    (folder / f"{ep}.{t:%a.%b.%d_%H_%M_%S}.GMT.{t:%Y}.CACO05.c2.timex.csv").write_text("\n".join(out) + "\n")


def extract(d, extra=()):
    out = d / "contours.csv"
    r = subprocess.run([sys.executable, str(HERE / "extract_elevation_contours.py"), str(d / "tide.csv"),
                        "--processed-dir", str(d / "processed"), "--output", str(out)] + list(extra),
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, cwd=str(d))
    used = set()
    if out.exists():
        with open(out) as f:
            used = {row["source_file"][:10] for row in csv.DictReader(f)}
    return r.returncode, used, r.stdout


def test_extractor():
    print("extract_elevation_contours.py: the partial-line floor")
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "processed").mkdir()
        frames = {}
        s = np.ones(300, bool)
        frames["full"] = (EPOCH, s)
        s = np.zeros(300, bool)
        s[150:270] = True                             # 40%: a partial line's clear run
        frames["partial40"] = (EPOCH + 1800, s)
        s = np.zeros(300, bool)
        s[200:280] = True                             # 27%: ~19% left after the extractor's own fragment trim
        frames["partial27"] = (EPOCH + 3600, s)
        s = np.zeros(300, bool)
        s[100:286] = True                             # 62%: a whole line whose truncated coverage is ~0.54
        frames["full62"] = (EPOCH + 5400, s)
        for ep, sig in frames.values():
            write_frame(d / "processed", ep, sig)
        tide = ["time,level"] + [f"{datetime.fromtimestamp(EPOCH + k * 600 - 3600, tz=timezone.utc):%Y-%m-%d %H:%M:%S},"
                                 f"{0.2 + 0.002 * k:.4f}" for k in range(30)]
        (d / "tide.csv").write_text("\n".join(tide) + "\n")
        ep = {k: str(v[0]) for k, v in frames.items()}

        rc, used, out = extract(d)
        check("the defaults: the whole line and the 40% partial line used, the one left with 19% not",
              rc == 0 and ep["full"] in used and ep["partial40"] in used and ep["partial27"] not in used,
              (rc, used, out[-600:]))
        check("the log counts the partial lines used", "partial line(s)" in out, out[-400:])
        rc, used, _ = extract(d, ["--partial-min-coverage", "0.5"])
        check("--partial-min-coverage 0.5: the 40% partial line not used", ep["partial40"] not in used, used)
        rc, used, _ = extract(d, ["--min-coverage", "0.6"])
        check("a whole line (>= 60% signal) still judged by --min-coverage: 0.6 drops the 62% line, "
              "the partial line kept by its own floor",
              ep["full62"] not in used and ep["partial40"] in used and ep["full"] in used, used)
        rc, used, _ = extract(d, ["--partial-min-coverage", "0"])
        check("--partial-min-coverage 0: no coverage floor for partial lines", ep["partial27"] in used, used)

    check("the detector and the extractor agree on the 60% that makes a line partial",
          X.FULL_LINE_SIGNAL_FRACTION == D.MIN_SIGNAL_COLUMN_FRACTION)


def test_syntax_38():
    print("Python 3.8")
    for f in ("waterline_detector_v5.py", "extract_elevation_contours.py", "test_partial_lines.py"):
        try:
            ast.parse((HERE / f).read_text(), filename=f, feature_version=(3, 8))
            good = True
        except SyntaxError as exc:
            good = False
            print(f"    {exc}")
        check(f"{f} parses as Python 3.8", good)


def main():
    test_clear_runs()
    test_detector()
    test_extractor()
    test_syntax_38()
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        raise SystemExit(1)
    print("all partial-line tests passed")


if __name__ == "__main__":
    main()

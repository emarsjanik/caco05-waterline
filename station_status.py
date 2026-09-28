#!/usr/bin/env python3
"""
Station Health Check
=====================
One summary of whether every data stream is arriving and every stage is
working, written after each cron run.

WHY. Failures on this station are silent: the maps and DEMs keep
rendering from older archived data and look fine when the cameras, the
GNSS-R processing, the buoy download or the detector have stopped.
Each check below catches one of those.

CHECKS (OK / WARN / ALERT):
  timex images      newest archived image          WARN > 26 h, ALERT > 72 h
  detections        newest archived waterline CSV  WARN > 48 h, ALERT > 5 d
  GNSS-R            newest water level             WARN > 3 d,  ALERT > 5 d
                    (gnssrefl normally runs ~2 days behind)
  GNSS-R QC         failed share, last 7 days      WARN > 10 %, ALERT > 30 %
  buoy waves        newest record                  WARN > 1 d,  ALERT > 3 d
  tide gauge        newest record                  WARN > 1 d,  ALERT > 3 d
  disk              free space on the data disk    WARN < 100 GB, ALERT < 20 GB
  last cron run     WARNING/ERROR lines logged     WARN on WARNING, ALERT on ERROR

Writes the report to stdout; exit status 0 all OK, 1 any WARN, 2 any ALERT,
so a mail or heartbeat hook can key on it.

Usage:
    python3 station_status.py --base /mnt/I2Rgus_Data/waterline \\
        --gnssr /home/argus_user/GNSS/.../usgs_spline_out.txt
"""

import csv
import shutil
import argparse
from pathlib import Path
from datetime import datetime, timezone

import numpy as np

LEVELS = {"OK": 0, "WARN": 1, "ALERT": 2}


def age_hours(epoch):
    return (datetime.now(timezone.utc).timestamp() - epoch) / 3600.0


def grade(value, warn, alert, higher_is_worse=True):
    if value is None:
        return "ALERT"
    if higher_is_worse:
        return "ALERT" if value > alert else "WARN" if value > warn else "OK"
    return "ALERT" if value < alert else "WARN" if value < warn else "OK"


def newest_mtime(folder, pattern):
    newest = None
    for p in Path(folder).glob(pattern):
        m = p.stat().st_mtime
        if newest is None or m > newest:
            newest = m
    return newest


def newest_epoch_in_csv(path, column="epoch"):
    if not Path(path).exists():
        return None
    last = None
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            try:
                last = max(last or 0, float(r[column]))
            except (KeyError, ValueError):
                continue
    return last


def fmt_age(h):
    return "missing" if h is None else (f"{h:.0f} h ago" if h < 48 else f"{h / 24:.1f} days ago")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--base", default="/mnt/I2Rgus_Data/waterline")
    ap.add_argument("--gnssr", help="gnssrefl spline file")
    ap.add_argument("--buoy-csv", default="archive/waves_44008.csv")
    ap.add_argument("--gauge-csv", default="archive/gauge_8447435.csv")
    ap.add_argument("--qc-csv", default="archive/gnssr_qc.csv")
    ap.add_argument("--log", default="logs/timex_cron.log")
    args = ap.parse_args()
    base = Path(args.base)
    rows = []

    def add(name, status, detail):
        rows.append((name, status, detail))

    t = newest_mtime(base / "archive/images_timex", "*.timex.jpg")
    h = age_hours(t) if t else None
    add("timex images", grade(h, 26, 72), f"newest archived {fmt_age(h)}")

    t = newest_mtime(base / "archive/processed_timex", "*.csv")
    h = age_hours(t) if t else None
    add("detections", grade(h, 48, 120), f"newest waterline {fmt_age(h)}")

    if args.gnssr:
        try:
            from extract_elevation_contours import load_gnssr_spline
            import contextlib, io
            with contextlib.redirect_stdout(io.StringIO()):
                ep, _, _, _ = load_gnssr_spline(args.gnssr)
            h = age_hours(float(ep[-1]))
        except Exception as exc:          # unreadable file is itself the finding
            h = None
            add("GNSS-R", "ALERT", f"cannot read {args.gnssr}: {exc}")
        else:
            add("GNSS-R", grade(h, 72, 120), f"newest water level {fmt_age(h)} (normal: ~2 days)")

    qc = base / args.qc_csv
    if qc.exists():
        flags, now = [], datetime.now(timezone.utc).timestamp()
        with open(qc, newline="") as f:
            for r in csv.DictReader(f):
                if float(r["epoch"]) >= now - 7 * 86400 - 3 * 86400:   # last week of available data
                    flags.append(int(r["qc_flag"]))
        if flags:
            share = 100.0 * np.mean(np.array(flags) >= 4)
            add("GNSS-R QC", grade(share, 10, 30), f"{share:.0f}% of the last week's readings failed QC")
        else:
            add("GNSS-R QC", "WARN", "no recent readings in the QC report")
    else:
        add("GNSS-R QC", "WARN", f"no QC report ({args.qc_csv})")

    for name, rel in (("buoy waves", args.buoy_csv), ("tide gauge", args.gauge_csv)):
        e = newest_epoch_in_csv(base / rel)
        h = age_hours(e) if e else None
        add(name, grade(h, 24, 72), f"newest record {fmt_age(h)}")

    try:
        free_gb = shutil.disk_usage(base).free / 1e9
        add("disk", grade(free_gb, 100, 20, higher_is_worse=False), f"{free_gb:.0f} GB free")
    except OSError as exc:
        add("disk", "ALERT", str(exc))

    log = base / args.log
    if log.exists():
        lines = log.read_text(errors="replace").splitlines()
        starts = [i for i, l in enumerate(lines) if "=== run start ===" in l]
        last_run = lines[starts[-1]:] if starts else lines[-200:]
        n_err = sum(("ERROR" in l) for l in last_run)
        n_warn = sum(("WARNING" in l) for l in last_run)
        status = "ALERT" if n_err else "WARN" if n_warn else "OK"
        add("last cron run", status, f"{n_err} ERROR line(s), {n_warn} WARNING line(s)")

    worst = max((LEVELS[s] for _, s, _ in rows), default=0)
    overall = [k for k, v in LEVELS.items() if v == worst][0]
    print(f"STATION STATUS {datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC: {overall}")
    for name, status, detail in rows:
        print(f"   {status:5s}  {name:14s} {detail}")
    raise SystemExit(worst)


if __name__ == "__main__":
    main()

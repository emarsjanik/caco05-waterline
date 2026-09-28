#!/usr/bin/env python3
"""
Archive A NOAA Tide Gauge Record In NAVD88
=============================================
Downloads a CO-OPS station's 6-minute water level, requested directly in
NAVD88, and merges it into a local CSV that only ever grows.

WHY. gnssr_qc.py checks every GNSS-R water level against what this
gauge predicts for Marconi. GNSS-IR degrades in big waves: on Sep 25-26
2026 it read up to ~4 m above the tide while Chatham showed a normal
tide with modest surge. Keeping the gauge record locally lets the check
run for any period, including ones older than a single download.

GAUGE. 8447435 Chatham, Lydia Cove, MA (see compare_gnssr_to_gauge.py for
how its tide relates to Marconi's). Data are "preliminary" until NOAA
verifies them; a later download of the same times replaces the rows.

OUTPUT columns: time_utc, epoch, level_navd88.

Usage:
    python3 fetch_tide_gauge.py --output archive/gauge_8447435.csv            (last 7 days)
    python3 fetch_tide_gauge.py --output archive/gauge_8447435.csv --days 120 (backfill)
"""

import csv
import argparse
from pathlib import Path
from datetime import datetime, timedelta, timezone

from compare_gnssr_to_gauge import fetch_gauge, to_epoch


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--station", default="8447435")
    ap.add_argument("--days", type=int, default=7, help="How many days back to download (default 7).")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    end = datetime.now(timezone.utc)
    start = end - timedelta(days=args.days)
    new = fetch_gauge(args.station, start, end)
    epochs = to_epoch(new["time"])

    out = Path(args.output)
    merged = {}
    if out.exists():
        with open(out, newline="") as f:
            for r in csv.DictReader(f):
                merged[int(r["epoch"])] = r
    before = len(merged)
    for t, e, lv in zip(new["time"], epochs, new["level"]):
        merged[int(e)] = {"time_utc": t.isoformat(), "epoch": int(e), "level_navd88": round(float(lv), 3)}

    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".tmp")
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["time_utc", "epoch", "level_navd88"])
        w.writeheader()
        for e in sorted(merged):
            w.writerow(merged[e])
    tmp.replace(out)
    first, last = min(merged), max(merged)
    print(f"Gauge {args.station}: {len(new)} row(s) downloaded, archive {before} -> {len(merged)}, "
          f"{datetime.fromtimestamp(first, tz=timezone.utc):%Y-%m-%d} to "
          f"{datetime.fromtimestamp(last, tz=timezone.utc):%Y-%m-%d %H:%M} UTC")


if __name__ == "__main__":
    main()

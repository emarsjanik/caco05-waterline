#!/usr/bin/env python3
"""
Archive Offshore Wave Observations From An NDBC Buoy
=======================================================
Downloads a buoy's "realtime2" standard meteorological file and merges
it into a local CSV that only ever grows.

WHY. NDBC keeps only the last 45 days in the realtime file, but every
waterline in the archive needs the wave conditions at its capture time:
in big waves the swash sits well above the still-water level (wave
setup), so a waterline's GNSS-R elevation is too low for where it lies
on the beach. On Sep 23 2026 frames taken in 3.8-3.9 m waves at 44008
(four times that week's median) plotted about half a metre high. Keeping
the record lets extract_elevation_contours.py tag each frame and
dem_from_contours.py leave rough-water frames out.

BUOY. 44008 (Nantucket Shoals, ~50 nm SE of Nantucket) faces the open
Atlantic like Marconi. It is offshore, so heights are an index of
conditions, not the height breaking on the beach. 44018 (off
Provincetown) returned 404 in Sep 2026; 44090 is in Cape Cod Bay and
does not see the Atlantic swell.

OUTPUT columns: time_utc, epoch, wvht_m, dpd_s, apd_s, mwd_deg. Missing
values ("MM" in the NDBC file) are left blank. Re-running is safe: rows
are keyed on time and a newer download replaces the same time.

Usage:
    python3 fetch_buoy_waves.py --station 44008 --output archive/waves_44008.csv
    python3 fetch_buoy_waves.py --input 44008.txt --output ...   (offline)
"""

import sys
import csv
import argparse
import urllib.request
from pathlib import Path
from datetime import datetime, timezone

URL = "https://www.ndbc.noaa.gov/data/realtime2/{station}.txt"
FIELDS = ["time_utc", "epoch", "wvht_m", "dpd_s", "apd_s", "mwd_deg"]
# Column positions in the NDBC standard meteorological format:
# #YY MM DD hh mm WDIR WSPD GST WVHT DPD APD MWD PRES ...
COLS = {"wvht_m": 8, "dpd_s": 9, "apd_s": 10, "mwd_deg": 11}


def parse_realtime(text):
    rows = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        p = line.split()
        if len(p) <= max(COLS.values()):
            continue
        try:
            t = datetime(int(p[0]), int(p[1]), int(p[2]), int(p[3]), int(p[4]),
                         tzinfo=timezone.utc)
        except ValueError:
            continue
        row = {"time_utc": t.isoformat(), "epoch": int(t.timestamp())}
        for name, i in COLS.items():
            row[name] = "" if p[i] == "MM" else p[i]
        if row["wvht_m"] == "" and row["dpd_s"] == "":
            continue
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--station", default="44008")
    ap.add_argument("--input", help="Parse a saved realtime2 file instead of downloading.")
    ap.add_argument("--output", required=True, help="Archive CSV to create or extend.")
    args = ap.parse_args()

    if args.input:
        text = Path(args.input).read_text()
    else:
        url = URL.format(station=args.station)
        try:
            text = urllib.request.urlopen(url, timeout=60).read().decode()
        except Exception as e:
            sys.exit(f"Download failed for {url}: {e}")
    new = parse_realtime(text)
    if not new:
        sys.exit("No wave rows found in the buoy file.")

    out = Path(args.output)
    merged = {}
    if out.exists():
        with open(out, newline="") as f:
            for r in csv.DictReader(f):
                merged[int(r["epoch"])] = r
    before = len(merged)
    for r in new:
        merged[r["epoch"]] = r

    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".tmp")
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for epoch in sorted(merged):
            w.writerow(merged[epoch])
    tmp.replace(out)

    first = min(merged)
    last = max(merged)
    print(f"Buoy {args.station}: {len(new)} row(s) downloaded, archive "
          f"{before} -> {len(merged)} row(s), "
          f"{datetime.fromtimestamp(first, tz=timezone.utc):%Y-%m-%d} to "
          f"{datetime.fromtimestamp(last, tz=timezone.utc):%Y-%m-%d %H:%M} UTC")


if __name__ == "__main__":
    main()

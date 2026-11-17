#!/usr/bin/env python3
"""
Tests For The Station Status Lens Check
=========================================
Plain python3 (no pytest needed on the station):

    python3 test_station_status.py

WHAT. station_status.py's "c2 lens" line reads owg_live.py's per-frame
sharpness (archive/owg_c2_H.csv) and grades the camera window:
  * clear days (median daylight sharpness well above 10) are OK;
  * the newest day blurred (median below 10) is a WARN: fog, rain or spray;
  * three blurred days running are an ALERT: clean the window. The 20-24 Oct
    2026 storm, with medians 3.4-4.7, is the case it was written for;
  * a blurred day followed by a clear one is OK again (only the run of days
    ending on the newest day counts);
  * night frames ("dark") and unreadable ones never count, and a day with
    fewer than 4 daylight frames is not judged;
  * no file, or no recent day with enough frames: no line at all.
And its "storm watch" line, from storm_watch.py's stamp (logs/storm_watch_last_run):
  * no stamp (the storm watch not installed): no line;
  * a fresh "ok" stamp is OK; 40 h old is a WARN; 5 days old an ALERT (its crontab line gone?);
  * "blind since <date>" (no forecast source answered) is OK for 2 days, a WARN from 3 days on, so a
    station that cannot reach the forecast hosts does not sit at WARN hiding real ones;
  * an "error" stamp is a WARN; the line never holds the words the cron log counts (WARNING, ERROR).
Exit status 0 when every check passes.
"""

import csv
import sys
import tempfile
from pathlib import Path
from datetime import datetime, timezone, timedelta

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from station_status import lens_check, storm_watch_check   # noqa: E402

NOW = datetime(2026, 10, 25, 22, 0, tzinfo=timezone.utc)
FIELDS = ["filename", "epoch", "time_utc", "brightness", "sharpness", "saturation", "status"]
failures = []


def check(name, cond, detail=""):
    print(f"  {'ok ' if cond else 'FAIL'} {name}" + (f"  ({detail})" if detail and not cond else ""))
    if not cond:
        failures.append(name)


def write(path, days):
    """days: {days before NOW's date: (sharpness of each daylight frame, n dark frames)}."""
    rows = []
    for back, (sharps, n_dark) in days.items():
        day = (NOW - timedelta(days=back)).replace(hour=0, minute=0)
        for i, sh in enumerate(sharps):
            t = day + timedelta(hours=11, minutes=30 * i)
            status = "blurred (fog, rain or wet lens)" if sh < 10 else "ok"
            rows.append([f"{int(t.timestamp())}.c2.bright.jpg", int(t.timestamp()), t.isoformat(),
                         160.0, sh, 0.01, status])
        for i in range(n_dark):
            t = day + timedelta(hours=1, minutes=30 * i)
            rows.append([f"{int(t.timestamp())}.c2.bright.jpg", int(t.timestamp()), t.isoformat(),
                         5.0, 0.5, 0.0, "dark"])
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(FIELDS)
        w.writerows(rows)
    return path


def run(tmp, days):
    return lens_check(write(Path(tmp) / "owg.csv", days), now=NOW.timestamp())


def main():
    clear = [48.8, 52.0, 61.3, 70.1, 55.0]
    storm = [3.4, 3.7, 4.1, 4.6, 3.9, 4.2]
    with tempfile.TemporaryDirectory() as tmp:
        print("lens check:")
        r = run(tmp, {0: (clear, 9), 1: (clear, 9), 2: (clear, 9)})
        check("clear days are OK", r and r[0] == "OK", r)
        check("the detail names the median and the day", r and "55.0" in r[1] and "25 Oct" in r[1], r)

        r = run(tmp, {0: (storm, 9), 1: (clear, 9), 2: (clear, 9)})
        check("the newest day blurred is a WARN", r and r[0] == "WARN", r)

        r = run(tmp, {0: (storm, 9), 1: (storm, 9), 2: (storm, 9)})
        check("three blurred days are an ALERT", r and r[0] == "ALERT" and "3 days running" in r[1], r)
        check("the ALERT says to clean the window", r and "clean the c2 window" in r[1], r)

        r = run(tmp, {0: (clear, 9), 1: (storm, 9), 2: (storm, 9), 3: (storm, 9)})
        check("a clear newest day after blurred ones is OK", r and r[0] == "OK", r)

        r = run(tmp, {0: (clear, 30), 1: (clear, 30)})
        check("dark frames do not pull the median down", r and r[0] == "OK", r)

        r = run(tmp, {0: ([3.0, 3.1, 3.2], 9), 1: (clear, 9)})
        check("a day with under 4 daylight frames is not judged", r and r[0] == "OK" and "24 Oct" in r[1], r)

        r = run(tmp, {10: (storm, 9)})
        check("nothing recent: no line", r is None, r)

        check("no file: no line", lens_check(Path(tmp) / "missing.csv", now=NOW.timestamp()) is None)

        print("storm watch line:")
        import os
        stamp = Path(tmp) / "storm_watch_last_run"
        check("no stamp: no line (not installed)", storm_watch_check(stamp, now=NOW.timestamp()) is None)

        def sw(text, age_h):
            stamp.write_text(text + "\n")
            t = NOW.timestamp() - age_h * 3600
            os.utime(stamp, (t, t))
            return storm_watch_check(stamp, now=NOW.timestamp())
        r = sw("ok (2026-10-25 12:50Z): 0 open event(s)", 9)
        check("a fresh ok stamp is OK", r and r[0] == "OK" and "0 open event" in r[1], r)
        r = sw("ok (2026-10-24 02:50Z): 1 open event(s)", 40)
        check("40 h old is a WARN", r and r[0] == "WARN", r)
        r = sw("ok (2026-10-20 12:50Z): 1 open event(s)", 5 * 24)
        check("5 days old is an ALERT naming --install-cron", r and r[0] == "ALERT" and "--install-cron" in r[1], r)
        r = sw("blind since 2026-10-25: no forecast source (om_marine unavailable)", 9)
        check("blind since today is OK (not yet 3 days)", r and r[0] == "OK" and "blind" in r[1], r)
        r = sw("blind since 2026-10-21: no forecast source (om_marine unavailable)", 9)
        check("blind for 4 days is a WARN", r and r[0] == "WARN" and "4 days" in r[1], r)
        r = sw("error in: run (2026-10-25 12:50Z): ERROR KeyError WARNING", 9)
        check("an error stamp is a WARN", r and r[0] == "WARN", r)
        check("the line holds no WARNING / ERROR (the cron log counts them)",
              r and "WARNING" not in r[1] and "ERROR" not in r[1], r)

    if failures:
        print(f"{len(failures)} station status check(s) FAILED")
        raise SystemExit(1)
    print("all station status tests passed")


if __name__ == "__main__":
    main()

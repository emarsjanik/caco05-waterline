# Storm watch on CACO-05: install, check, and turn off

For the person at the CACO-05 NUC (Marconi Beach), logged in as `argus_user`. About 20 minutes,
plus a few minutes for the replay. Nothing here changes the waterline cron, owg.sh or the GPS jobs.

What it does once installed: every morning at 07:50 (the NUC's clock) it reads the forecast for
Marconi Beach. When a strong storm is forecast within 4 days it emails a heads-up (so a GPS or drone
survey can be planned), builds a beach map of the calm days before the storm, waits for the storm to
pass, builds maps of the days after it, and emails the erosion/accretion map. The camera maps only the
intertidal beach (between the low and high waterlines); the dune and bluff are not measured.

Type each command on its own, one block at a time, and compare with what the text says to expect.

## Not on CACO-05?

This page is for CACO-05 only. A NUC set up by coastcam-station-setup gets the storm watch from setup
(its own facts in `storm_watch_station.conf`, alerts mode until its beach maps exist, its crontab line
from `./setup.sh --only 10`): see `docs/STORM_WATCH_ALERTS_MODE.md`. On CACO-05 none of that applies: it
has no `storm_watch_station.conf`, `STORM_MODE` is `full`, and every step below works as written.

---

## 1. Go to the station folder

```
cd /mnt/I2Rgus_Data/waterline
```

## 2. Check that no tracked file has local edits

```
git status --short
```

Expect no line starting with ` M`. (Lines starting with `??` are fine.) If a tracked file shows ` M`,
stop and ask before going on: the pull would refuse.

## 3. Get the code

```
git pull
```

## 4. Run the tests

```
python3 test_storm_watch.py
```

Takes about 1 minute. The last line must be `all storm watch tests passed`.

```
python3 test_station_status.py
```

The last line must be `all station status tests passed`.

```
bash -n storm_watch.sh && echo "storm watch syntax OK"
```

## 5. Check that the NUC reaches the forecast services (writes and sends nothing)

```
./storm_watch.sh --dry-run --save-raw /tmp/sw_raw
```

If a storm is forecast (or the buoys show one), it first prints each email it would send: each starts
with a line like `=== 2026-11-17_fc: heads_up ===`, then the recipients and the whole text. That is only
on your screen: nothing is sent and nothing is written to a log. Then comes one long line starting with `ok`, `degraded` or `blind`, and the output
ends with lines starting `DRY RUN: nothing written, nothing sent.` The line starting with `ok`,
`degraded` or `blind` is the one that matters:

- `ok`: the wave forecast was read.
- `degraded`: no wave forecast, only the NWS alerts and/or wind.
- `blind`: no forecast service answered (the NUC's network may block them).

Now check the saved answers against the parsers:

```
python3 test_storm_watch.py --real-captures /tmp/sw_raw
```

Copy everything it prints into your report to the developers (it holds no password or address).
The first line starts `REAL captures in /tmp/sw_raw`. Each check line after it starts with `ok`
(indented by two spaces); the `alerts kept: ...` and `trigger: ...` lines are information, not checks.
The last line should be `all storm watch tests passed`. If any line starts with `FAIL` instead, the
last lines count and repeat the failures: each names the service whose answer was not understood (the
storm watch still runs without that service). Send them to the developers as they are.

## 6. Replay the September 2026 storm (a check on real data; writes only to /tmp, sends nothing)

Run this outside the waterline cron's times (12:30, 19:25 and 20:55); it refuses on its own if a cron
run is due within 15 minutes.

```
python3 storm_watch.py --replay 2026-09-15 2026-10-20 --out-dir /tmp/sw_replay
```

Takes a few minutes. Expect a table ending with one event for the storm of 23-28 Sep 2026, a pre-storm
window ending about 22 Sep, post-storm windows from about 28-29 Sep, and a `change post7` line giving
the cells, the median change and the erosion and accretion volumes. If the storm came in two bursts, a
second short event closed as `merged` may also be listed: that is expected. Copy the table into your
report.

```
cat /tmp/sw_replay/replay_summary.json
```

## 7. Count how often it would have triggered (needs the internet; about 2 minutes)

```
python3 storm_watch.py --climatology --station 44008 --years 2010-2025
```

Copy the lines `3.0 m: ...`, `3.5 m: ...` and `4.0 m: ...` into your report. The default trigger is
3.5 m; we expect about 6-10 events a year.

## 8. Optional: a contact for the National Weather Service

The NWS asks automated users to identify themselves. The storm watch sends a generic name; you may
add the team's address for the NWS only. First copy the example settings file:

```
cp -n storm_watch_settings.conf.example storm_watch_settings.conf
```

```
nano storm_watch_settings.conf
```

Find the line `#NWS_CONTACT=""`, remove the `#`, type the team's email address between the two quote
marks, then save (Ctrl+O, Enter) and exit (Ctrl+X). This file is not tracked by git, so it never
blocks a `git pull`. Leave everything else as it is.

## 9. Install the daily run

```
./storm_watch.sh --install-cron
```

Expect it to print:

```
Crontab now (previous copy saved in /mnt/I2Rgus_Data/waterline/logs/crontab.before_storm_watch_install):
50 7 * * * /mnt/I2Rgus_Data/waterline/storm_watch.sh # storm_watch.sh (storm watch): daily forecast check
```

That is the whole cron line it adds (do not type it by hand). Check that the other lines are still
there:

```
crontab -l
```

Expect the three `waterline_timex_cron.sh` lines, the two `owg.sh` lines and the GPS_code lines as
before, plus the new storm watch line.

## 10. The next morning (after 07:50)

```
tail -3 logs/storm_watch.log
```

Expect a line dated today starting with `ok` (or `degraded` / `blind`, see step 5).

```
cat logs/storm_watch_last_run
```

```
./storm_watch.sh --status
```

Expect `no events` on a calm week. `--status` only reads: it never changes a file. The station status
(`station_status.txt` after the next cron run) now has a `storm watch` line.

## 11. Turning it off

Stop the daily run (nothing else changes; any event folders stay in `archive/storm_watch/`):

```
./storm_watch.sh --remove-cron
```

Expect it to print `(no storm watch line: the storm watch is off)`.

To keep it running but send no emails (they are still written to each event's `emails/` folder):

```
echo 'STORM_EMAIL=0' >> storm_watch_settings.conf
```

To send emails again, remove that line:

```
sed -i '/^STORM_EMAIL=0$/d' storm_watch_settings.conf
```

---

Where things are: `logs/storm_watch.log` (one line a day), `logs/storm_watch_detail.log` (everything a
run printed), `archive/storm_watch/events/` (one folder per storm: the maps, GeoTIFFs, change map and
the emails). The full description is in `docs/VALIDATION.md` section 2.10 and at the top of
`storm_watch.py`.

**Storm frames for the weekly videos.** Each storm's final products (the pre-storm DEM, the first and the
7-day post-storm DEMs, and both change maps) are also drawn in the look of the weekly frames, with a red
banner such as `STORM <event>: BEFORE (pre-storm, ...)`, into `archive/storm_watch/events/<event>/frames/`
(with a `frames.json` list). The waterline cron's weekly step (`weekly_frames.py daily`, see
`docs/WEEKLY_FRAMES.md`) puts them into `archive/weekly/video/dem_weekly.mp4` and `waterlines_weekly.mp4`
at their dates. The frames need that step's frozen look (`archive/weekly/frame_spec.json`, made by the
cron's first run after the pull): until it exists the storm watch logs a NOTE and draws them on a later
run while the event is open. A frame that cannot be drawn is a NOTE in the event's `log.txt`; the event and
its emails go on. The video needs ffmpeg (`sudo apt-get install -y ffmpeg`).

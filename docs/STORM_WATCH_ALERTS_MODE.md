# Storm watch at a new station, and alerts mode

The storm watch (`storm_watch.py`, run daily by `storm_watch.sh`) was written for CACO-05 (Marconi
Beach). It now runs at any station: every place fact is a setting, and a new station that has no beach
maps yet runs it in **alerts mode**. On CACO-05 nothing here changes anything: it has no
`storm_watch_station.conf` and no `waterline_station.conf`, `STORM_MODE` is `full`, and every run writes
the same files and sends the same emails as before (the defaults are CACO-05's values).

coastcam-station-setup sets all of this up on a new NUC. This page says what it sets, what the storm
watch then does, and what to check or change by hand.

---

## 1. The two modes

| | `STORM_MODE=full` (the default; CACO-05) | `STORM_MODE=alerts` (a new station, before its beach maps) |
|---|---|---|
| Forecast check, heads-up email | yes | yes |
| Follows the storm on the observed waves (buoy) | yes | yes |
| Pre- and post-storm DEMs, change map | yes | **no**: the ground file and the waterline cron's lock are never touched |
| When the storm is over | waits for the post-storm DEMs (about 9-10 days) | a **storm-over** email at once |
| The event closes | when the 7-day post-storm map is done | **72 h after the storm's end** |
| `archive/storm_watch/events.csv` | not written | **one row per event**, when it closes |
| Email wording | DEM and intertidal-beach notes | no DEM or NAVD88 wording |

Why 72 h: a forecast that still shows the same storm (an older model run, the storm's tail) belongs to an
open event for 72 h after its end, so it updates that event instead of opening a new one with a second
heads-up. After 72 h the event closes and its row goes to `events.csv`.

The other closing rules are the same in both modes: a forecast storm that never comes is closed as a
**false alarm**, one whose waves the buoys could not see as **unverified**, one that came as part of another
event's storm as **merged**, a storm that never calms is **given up** 14 days after its peak. In alerts
mode none of them leaves a DEM to build later.

Switching: coastcam-station-setup writes `STORM_MODE=alerts` at install and `STORM_MODE=full` once the
station's beach maps are set up (`./setup.sh --waterline`). An event still open at the switch to full mode
is then followed as a full-mode event (its DEMs are built if the waterline points exist). Switching back
from full to alerts closes open events 72 h after their storm's end, without their DEMs.

### events.csv

`archive/storm_watch/events.csv`, one row per closed event (alerts mode only), never an address:

`event, kind, mode, station, opened_utc, closed_utc, closed_reason, forecast_onset_utc, forecast_end_utc,
forecast_peak_hs_m, forecast_peak_utc, forecast_rules, observed_onset_utc, observed_peak_hs_m,
observed_peak_utc, buoy, buoy_max_hs_m, storm_end_utc, storm_dates, emails`

`kind` is `fc` (opened on the forecast) or `obs` (a storm the forecast missed); `closed_reason` is
`complete`, `false_alarm`, `unverified`, `merged` or `gave_up`; `emails` lists each email's kind and
outcome (`heads_up:sent storm_over:sent`, or `:dry` when emails are off).

---

## 2. The settings files

Two files in the station folder (`/mnt/I2Rgus_Data/waterline`), read in this order; **the last value
wins**:

1. `storm_watch_station.conf`: the station's facts. coastcam-station-setup writes it from its
   `station.conf`; do not edit it (run setup again instead). CACO-05 has none.
2. `storm_watch_settings.conf`: your own edits (start from `storm_watch_settings.conf.example`). They win.

Both are untracked by git. Format: `KEY=VALUE` at column 0, no spaces around `=`, a value with spaces in
double quotes, `#` comments.

### The station facts

A `storm_watch_station.conf` must give every one of these (`--check-settings` refuses the file otherwise):

| Key | What it is | CACO-05's value | From setup's `station.conf` |
|---|---|---|---|
| `STORM_STATION` | the station: email subject `[<station> storm watch]`, User-Agent `<station>-storm-watch (...)` | `CACO-05` | the station name |
| `STORM_SITE` | the beach: titles and email text | `Marconi Beach` | the ZIP lookup's place |
| `STORM_TZ` | local times in the emails | `America/New_York` | `STATION_TZ` |
| `STORM_POINT` | the beach (NWS `/points`, Open-Meteo wind) | `41.89,-69.96` | `STORM_POINT` |
| `STORM_MARINE_POINT` | offshore of the beach (Open-Meteo marine) | `41.90,-69.80` | `STORM_MARINE_POINT` |
| `STORM_BUOY_POINT` | the buoy's position (Open-Meteo marine) | `40.50,-69.25` | the buoy table |
| `STORM_BUOY` | the offshore NDBC wave buoy | `44008` | the first of `STORM_BUOYS` |
| `STORM_BUOY2` | a second buoy converted to the beach; empty = none | `44013` | the second of `STORM_BUOYS`, or empty |
| `STORM_ALERT_ZONES` | NWS zones for the alerts | `ANZ254,MAZ022` | `STORM_ALERT_ZONES` |
| `STORM_FACING` | the direction the beach faces, deg true; empty = no facing test | empty | `STORM_FACING` |
| `STORM_MIN_CELL_LON` | Marconi's Open-Meteo longitude test; empty = off | `-69.93` | empty |

The buoy names, one mapping: setup's `STORM_BUOYS` (primary first) gives `STORM_BUOY` and `STORM_BUOY2`
here; the waterline overlay's `WAVE_BUOY` (in `waterline_station.conf`) is the same buoy as `STORM_BUOY`.

Optional: `STORM_GFSWAVE_STATION` (the GFS-Wave point bulletin; empty = `STORM_BUOY`, `none` = no
bulletin, for a buoy NOAA publishes none for), `STORM_MODE`, `STORM_EMAIL`, `STORM_EMAIL_TO` and every
threshold in the `.example`.

### What the facts change

* **Buoys.** The observed waves are read from `archive/waves_<STORM_BUOY>.csv` (times
  `STORM_OFFSHORE_TO_MARCONI`, 0.83, for waves from the sector) and `archive/waves_<STORM_BUOY2>.csv`
  (converted with `calibration/buoy_transfer_<STORM_BUOY2>.json` when that exists, else also times 0.83;
  CACO-05's own 44013 fit is used only on CACO-05). The storm watch fetches `STORM_BUOY` when its archive is
  more than 6 h old, and on a setup-managed station (no owg.sh there) `STORM_BUOY2` as well. The GFS-Wave
  bulletin, its archived copy (`forecasts/<date>/<HHMMZ>_gfswave_<id>.bull.gz`), the email labels and
  `--climatology` all use the buoy's id.
* **The Open-Meteo sea cell.** Open-Meteo picks a grid cell near `STORM_MARINE_POINT`. With
  `STORM_FACING`, a cell whose bearing from the beach lies more than 80 deg off the facing (a bay behind
  the beach, or land) is refused for that run. Marconi's older test (a cell west of `STORM_MIN_CELL_LON` is
  in Cape Cod Bay) is off when `STORM_MIN_CELL_LON` is empty, as setup writes it.
* **Recipients.** owg.sh's `RECIPIENTS` (CACO-05's list) are used only on CACO-05: never when a
  `storm_watch_station.conf` exists or in alerts mode. There, an empty `STORM_EMAIL_TO` means no email:
  the bodies are still written to `archive/storm_watch/events/<id>/emails/`, and the run line says
  `WARNING: no email recipients`.

---

## 3. Checking the settings: --check-settings

```
cd /mnt/I2Rgus_Data/waterline
./storm_watch.sh --check-settings
```

It reads both files strictly and prints what the storm watch will use (station, site, zone, mode, buoys,
bulletin, alert zones, the recipients' count and where they come from; never an address). Exit 1, with an
`ERROR:` line each, for:

* a line that is not `KEY=VALUE` (for example `STORM_HS_M = 4.0`: such a line is never read);
* a key the storm watch does not know;
* a value that does not parse (a number, a sector `a-b`, a point `lat,lon`, a time zone, a buoy id,
  `full`/`alerts`, an open quote, text after a value);
* a `storm_watch_station.conf` without one of the station facts above.

`WARN:` lines (exit 0): a station other than CACO-05 using one of Marconi's facts (its buoys, points,
zones or longitude test; check them: the nearest buoy is not always the right one on a cape), thresholds
still at Marconi's values, no sea-cell test at all, emails on without recipients.

coastcam-station-setup runs it after it writes `storm_watch_station.conf`, and installs no crontab line
unless it passes. A daily run never stops on these: an unknown key or a bad value is noted in
`logs/storm_watch.log` and its default is used.

---

## 4. Thresholds for a new station

The trigger and confirm heights are Marconi's (docs/VALIDATION.md 2.10 and 4.9). Suggestions from the
buoy's NDBC history, by the same rules (needs the internet, about 2 minutes):

```
python3 storm_watch.py --climatology --years 2010-2025 --emit-settings
```

It counts the storms each threshold would have opened per year, then prints `KEY=VALUE` lines, each
with its reason as a comment: `STORM_HS_M` and `STORM_OBS_OPEN_BUOY_HS` at the buoy's 99th percentile;
`STORM_HS_SUPPORT_M`, `STORM_SWELL_HS_M` and `STORM_CONFIRM_BUOY_HS` at its 98th; `STORM_CONFIRM_HS` and
`STORM_OBS_OPEN_HS` at those times 0.83; with `STORM_FACING`, `STORM_DIR_SECTOR` as the facing +- 80 deg.
These are for review by hand: copy the lines you keep into `storm_watch_settings.conf`, then run
`./storm_watch.sh --check-settings`.

---

## 5. Emails on

At install the emails are written, not sent (`STORM_EMAIL=0`). To send them, in
`storm_watch_settings.conf`:

```
STORM_EMAIL=1
STORM_EMAIL_TO="<address> <address>"
```

(coastcam-station-setup's `tools/waterline_setup.sh storm-email on` does the same.) Then
`./storm_watch.sh --check-settings` shows the number of recipients.

---

## 6. When the beach maps are on (full mode at a new station)

* The storm DEMs use the waterline cron's DEM options, and `waterline_station.conf`'s `DEM_MAX_HS` after
  them, as the cron does (a new station's first maps use 1.0 m). The calm level that ends a storm
  (`STORM_CALM_HS`, empty) follows it.
* `storm_watch.sh` checks `waterline_station.conf` with `waterline_station_check.py check-overlay` and
  sources it before a run, so its exported settings (`STATION_TZ`, `GEOTIFF_EPSG`, `STATION_NAME`,
  `STATION_PROVISIONAL`, ...) reach the DEM builds, as in the cron. A file it refuses stops the run: an
  `ERROR` line in `logs/storm_watch.log`, exit 1. `--dry-run`, `--status`, `--check-settings` and
  `--climatology` build nothing and skip this.

---

## 7. The crontab on a setup-managed NUC

Where `~/stationTools/crontab_new_station.txt` exists, coastcam-station-setup builds the crontab, with
storm_watch.sh's own line and tag. `./storm_watch.sh --install-cron` and `--remove-cron` then refuse and
say what to do instead: set `STORM_WATCH="yes"` (or `"no"`) in setup's `station.conf` and run
`./setup.sh --only 10` in the coastcam-station-setup folder.

---

## 8. Where to look

| | |
|---|---|
| `logs/storm_watch.log` | one line per run (`alerts mode (no DEM)` in alerts mode) |
| `logs/storm_watch_last_run` | the stamp station_status.py reads |
| `archive/storm_watch/state.json` | the events |
| `archive/storm_watch/events/<id>/emails/` | every email's text, sent or not |
| `archive/storm_watch/events.csv` | alerts mode: one row per closed event |
| `./storm_watch.sh --status` | the events and what each waits for |

Tests: `python3 test_storm_watch.py` (T14 another station, T15 alerts mode, T16 storm_watch.sh on a
setup-managed NUC).

# Weekly DEM and waterline frames, and their videos

`weekly_frames.py` gives every week a **brand new DEM** and a **plan view of its waterlines**, drawn as
frames that all share one frozen look. Played in order, the frames show the intertidal beach change week by
week (erosion and accretion). Every Monday a new week starts. The past week is set aside: it is marked
provisional at once and becomes final when its data are complete.

There are two series, each with its own MP4 video:

| Series | Each frame shows | Video |
|---|---|---|
| DEM | the week's DEM: `dem_from_contours.py` with exactly the cron's window-DEM options, over Monday..Sunday | `archive/weekly/video/dem_weekly.mp4` |
| Waterlines | the week's waterline points from both cameras, coloured by elevation on the same scale | `archive/weekly/video/waterlines_weekly.mp4` |

The storm watch's pre-storm and post-storm DEMs are drawn in the same frame style, with a red banner. They go
into the videos at their dates (section 5).

---

## 1. The week is a week of UTC dates

A week runs from **Monday to Sunday of UTC dates**. This is the pipeline's DEM day: `dem_from_contours.py`
picks points by the date of `capture_time_utc`. A week is named by its Monday (`YYYY-MM-DD`), and every title
says "(UTC dates)":

    CACO-05 Marconi Beach: week of Mon 16 Nov to Sun 22 Nov 2026 (UTC dates)

**Sunday evening captures after 00:00 UTC belong to the next week.** That is 20:00 EDT in summer and 19:00
EST in winter. The same split already applies to every 7-day DEM.

## 2. When a week is final

| Status | Meaning | Where its frames are |
|---|---|---|
| `current` | the week that holds today's UTC date; rebuilt by every daily run | `archive/weekly/current/` |
| `provisional` | the week has ended, but the ground file does not yet cover its Sunday (GNSS-R runs about 2 days behind). Stamp: `PROVISIONAL: data through <date>` | its week folder only |
| `final` | the ground file covers the Sunday (it holds points of a later date), or `WEEKLY_MAX_WAIT_DAYS` (7) days have passed since the Sunday. In the second case the stamp is `data through <date>` | `archive/weekly/frames/{dem,waterlines}/<sunday>_week.png` (the videos' input) |
| `empty` | a final week with no DEM: no points, the day test rejected every camera-day, or fewer than 10 cells. Its frames say `EMPTY` and why, so the timeline stays continuous | as `final` |

A final week is never rebuilt. `respec` (section 6) only redraws it.

A provisional week is rebuilt only when its ground rows change: their number, or their content (a digest of
the rows, so a ground file the cron re-processed with the same rows but other elevations rebuilds it too). A
failed DEM build (an exit code other than 0, 3 or 4) is a `WARNING`. That week is never finalised; the next
UTC day's daily run builds it again (the work is done once per UTC date), or at once
`python3 weekly_frames.py daily --base BASE --force`.

On the first daily run of a Monday, last week becomes provisional and the current-week frames start again
from nothing. With the GNSS-R lag, the new week's frames say `EMPTY: no waterline points in this week` until
about Wednesday.

## 3. Files

Everything is under `BASE/archive/weekly/` (BASE = `/mnt/I2Rgus_Data/waterline` on CACO-05).

| Path | What |
|---|---|
| `frame_spec.json` | the frozen look of every frame (section 6) |
| `state.json` | the last daily run, the spec's hash, each week's status, the last video build |
| `weeks/<monday>/` | `wk_<monday>_<sunday>_{dem,spread,count,sigma}.asc` and `_info.json` (the week's DEM); `points.csv` (the week's waterline points: time, camera, frame, E, N, elevation); `dem_frame.png`, `waterline_frame.png`; `week.json` (status, data through, usable days, frames, points, cells, the DEM's exit code, why, built and finalised times); `dem_log.txt` (the build's output) |
| `current/` | `dem_frame.png`, `waterline_frame.png` of the week in progress |
| `frames/dem/`, `frames/waterlines/` | `<sunday>_week.png` of every final week |
| `video/` | `dem_weekly.mp4`, `waterlines_weekly.mp4`, and the ffmpeg concat lists used (`<series>_frames.txt`) |
| `work/` | a working copy of the ground rows, removed at the end of each run |

The script also writes `logs/weekly_frames.lock`, its own lock. It writes nothing else.

The weekly grids **never** go into `archive/dems`: `dem_change.py --series` reads every `dem_*_7d` file there,
and the weekly builds use no `--series-dir`.

## 4. The daily step (the cron)

The waterline cron runs, near its end and only while `DEM_ENABLE=1`:

    python3 "$BASE/weekly_frames.py" daily --base "$BASE" --have-lock

- `--have-lock` tells the script that the cron already holds its lock (`.timex_cron.lock`), so the ground
  file is read directly. Run by hand, without `--have-lock`, the script takes that lock just long enough to
  copy the ground rows. It waits up to 10 minutes for it, and never reads the file when a cron run is due
  within 15 minutes.
- The work is done **once per UTC date**. Later runs that day print `weekly frames: done today already` and
  return at once. `--force` redoes the work.
- One step makes the frozen spec if there is none, builds the current week, moves ended weeks to provisional,
  finalises covered ones and rebuilds a video whose frames changed.
- Missed weeks (the cron was off, or a backfill stopped with `--to`) are built at most 4 per daily run, the
  oldest first, because the step runs inside the cron's lock. The rest wait for the next days' runs, and a
  NOTE says how many are left; `backfill` builds them all at once.
- Each action prints one line starting `weekly frames:`. A `NOTE:` line is information. A `WARNING:` line is a
  real fault (a failed DEM build); the station status counts those. The videos never make a WARNING
  (section 5).
- Exit code: 0 when done, 1 after a WARNING. On exit 1 the cron adds its own line,
  `WARNING: weekly frames failed (see above)`, and goes on to the station status (step 8).
- `WEEKLY_ENABLE=0` turns the step off.

## 5. The videos (ffmpeg)

Install ffmpeg once on the NUC:

    sudo apt-get install -y ffmpeg

Without ffmpeg the frames are still made. The log gets one `NOTE` a day (never a WARNING) and the videos wait.
An ffmpeg that fails (for example one built without libx264, or a full disk) is a `NOTE` too, with its last
error line: the run's exit code stays 0, the earlier MP4 is kept, and the next daily run tries again.

A video is rebuilt when the set of its frames, or one of them, changed since its last build (their list of
names, times and sizes). This happens at most once per run. Each frame is shown for `1/WEEKLY_FPS` seconds
(default 2 frames per second); storm frames are shown for the same time. The video is H.264, `yuv420p`, at
the frame size, written to a temporary file and then moved into place. A video whose frames did not change
is not rebuilt. To force one, delete the MP4 and run `python3 weekly_frames.py daily --base BASE --force`
(`respec` always rebuilds both).

What goes in, in date order (the first 10 characters of each name):

- the final weeks, `frames/<series>/<sunday>_week.png`;
- the storm watch's frames, read in place from
  `archive/storm_watch/events/<event>/frames/frames.json` (the storm watch writes only under
  `archive/storm_watch/`). Each entry of that list is
  `{"part", "series": "dem" | "waterlines" | "change", "date", "window", "png"}`:
  - `date` is the window's end date: the pre window's end for the pre-storm frame, the post window's end for
    a change frame;
  - `png` is a file name in that `frames/` folder;
  - a storm frame sorts as `<date>_storm_<event>_<part>.png` (`weekly_frames.storm_frame_name`);
  - **change frames go into the DEM video, after that date's DEM frames**;
  - an entry that does not check out (another series, a path, a missing file) is left out with a NOTE.

Current and provisional weeks are never in a video.

## 6. The frozen look, and `respec`

The first daily run (or a backfill) writes `frame_spec.json`, and every frame after it uses that look:

- **extent**: every usable ground point (georectified, passing `--max-hs`), padded 20 m and rounded outward
  to 10 m, or `WEEKLY_EXTENT`;
- **elevation colours**: the 1st to 99th percentiles of those points, rounded outward to 0.25 m, or
  `WEEKLY_Z_LIMITS`; colour map `terrain`, with the same limits in both series;
- **change limit**: `WEEKLY_CHANGE_LIMIT`, symmetric;
- **pixel size**: `WEEKLY_SIZE`, at a fixed dpi;
- **fps**, the cron's **cell**, the **EPSG code and frame label** (`station_env`), and the **station name**
  (`STATION_NAME`);
- `made_from`: the ground file, its rows, and the first and last dates.

Each PNG carries the spec's hash in its metadata (`weekly_frames_spec`).

**A setting changed later is NOT applied to the existing frames.** The daily step logs a NOTE naming
`respec`. To apply it:

    python3 weekly_frames.py respec --base /mnt/I2Rgus_Data/waterline

This re-freezes the spec from the current settings and the ground file. The old spec is kept as
`frame_spec.prev.json`. Every week is then redrawn from its kept grids and `points.csv` (no DEM is rebuilt),
and the videos are rebuilt.

Storm frames are drawn by `storm_watch.py` and keep the look they were drawn with. The videos scale and pad
them to the frame size, and `respec` says how many there are. (Each frame's time in the video comes from the
concat list, never from ffmpeg's frame count, which restarts where the frame size changes; so a video of
mixed sizes keeps every frame.) The storm watch draws a product's frames in
the run that makes it final, and only once `frame_spec.json` exists: before that it logs a NOTE and tries
again on its next run while the event is open (the storm watch never writes `frame_spec.json` itself).
Their banners use the storm watch's own settings: `BEFORE (pre-storm, <start> to <end>)`, `AFTER (first
<STORM_POST1_USABLE> usable days, ...)`, `AFTER (<STORM_POST_DAYS> days, ...)`. In alerts mode there are
none.

## 7. Backfill (once, over the archive)

    python3 weekly_frames.py backfill --base /mnt/I2Rgus_Data/waterline --dry-run   # what it would do
    python3 weekly_frames.py backfill --base /mnt/I2Rgus_Data/waterline             # do it
    python3 weekly_frames.py backfill --base ... --from 2026-09-28 --to 2026-11-15   # a range (their weeks)

- It covers every **ended** week from the ground file's first UTC date (or `--from`) to the last ended week
  (or `--to`).
- It builds each week that is not final and skips final ones. So it can be **resumed**: run it again after a
  stop.
- It takes the cron's lock only to copy the ground file, then builds the videos.
- A DEM build takes about as long as the cron's 7-day DEM, so allow a few minutes per week of archive. Run it
  away from the cron's times (12:30, 19:25, 20:55 local), e.g. under `nohup` in the morning.
- `--dry-run` writes nothing.

The daily step then carries on from where the backfill stopped.

## 8. Settings

These live in `BASE/weekly_frames_settings.conf` (optional, untracked). Copy
`weekly_frames_settings.conf.example` to start one. Use `KEY=VALUE` lines. An unknown key, or a value that
does not parse, is a NOTE, and the default is used.

| Key | Default | Meaning |
|---|---|---|
| `WEEKLY_ENABLE` | `1` | 0: the daily step does nothing |
| `WEEKLY_FPS` | `2` | video frames per second (0.1 to 60) |
| `WEEKLY_MAX_WAIT_DAYS` | `7` | an ended week is final at the latest this many days after its Sunday |
| `WEEKLY_EXTENT` | `""` (auto) | `e_min,e_max,n_min,n_max` in metres |
| `WEEKLY_Z_LIMITS` | `""` (auto) | `z_min,z_max` in m NAVD88 |
| `WEEKLY_CHANGE_LIMIT` | `0.5` | the change frames' colour limit (m) |
| `WEEKLY_SIZE` | `1920x1080` | frame pixels, even numbers |

The DEM options are not settings here. They are the cron's own: `DEM_MAX_HS`, `DEM_MAX_DAY_OFFSET`,
`DEM_DAY_BLOCK_LOD`, `DEM_CELL`, `DEM_MIN_POINTS` and `DEM_MAX_SPREAD`, read from `waterline_timex_cron.sh`
through `storm_watch.cron_settings` and `storm_watch.dem_args`, so they cannot drift. A new station's
`waterline_station.conf` `DEM_MAX_HS` applies, as it does for the cron.

**A new station's own values in a run by hand.** The cron sources `waterline_station.conf` and exports the
station's names (`STATION_NAME`, `GEOTIFF_EPSG`, `STATION_TZ`, `STATION_WAVE_CURRENCY` and the others on its
`export` line). Every `weekly_frames.py` command does the same by itself: it checks the file as the cron does
(`waterline_station_check.py check-overlay`) and puts the exported values into its environment. So a
`backfill` or `respec` run from a shell freezes the station's own name and frame and builds its DEMs with
the station's time zone and waves, exactly like the cron's run. A file the check refuses stops the command
with an `ERROR` (exit 1) and nothing is written. CACO-05 has no such file.

## 9. What a frame shows

- **DEM frame.** The week's DEM on the frozen extent and colour scale, with:
  - contours every 0.5 m, a colour bar with fixed ticks, a scale bar and a north arrow;
  - the title, a subtitle (usable days, camera frames, cells, NAVD88, the frame label), and the status stamp.
- **Waterline frame.** The week's points (`easting_utm19`, `northing_utm19`), coloured by
  `beach_elevation_navd88` (else `tide_elevation_navd88`) and drawn in time order:
  - frames with offshore Hs above the cron's `DEM_MAX_HS` are left out, as the DEM leaves them out;
  - the DEM's day test is **not** applied to the drawn points, and the subtitle says so;
  - the subtitle gives the counts of waterlines and points.
- **Change frame** (storm watch only). `dem_change.py`'s difference grid, `RdBu` with the frozen limit:
  - cells within the level of detection are light grey, and cells not compared are the panel's grey;
  - the banner gives both windows (`STORM <event>: CHANGE before -> after (pre <a> to <b> -> post <c> to <d>)`),
    the title which post-storm DEM it is (`first 3 usable days` or `7 days` at the storm watch's defaults),
    and the subtitle the volumes and the level of detection.

The cameras' positions are not drawn. They would need the calibration files, which this script never reads.

## 10. Blind date

The script reads only the ground file, its own outputs, the storm watch's `frames.json` and frames, and the
cron's settings files (`waterline_timex_cron.sh`, a new station's `waterline_station.conf`). It
never opens, globs or stats anything under `calibration/`, `surveys/` or `ground_truth/`, nor anything named
for a blind date (`survey_products.BLIND_DATES`, 2024-10-23, in any of its name forms).

Both `backfill` and `daily` skip any week that contains a blind date, and say so. The ground rows of a blind
date are never copied or drawn. `test_weekly_frames.py` checks all of this with an audit hook on every open.

## 11. Checking it

    python3 weekly_frames.py status --base /mnt/I2Rgus_Data/waterline   # read only
    python3 test_weekly_frames.py                                       # ends "all weekly frames tests passed"

The test takes about 45 s. It:

- runs the real `dem_from_contours.py` for two weeks and checks the exact cron options;
- follows a week through current, provisional and final (by coverage and by the maximum wait), plus an empty
  week and a failed build;
- checks the once-per-day stamp, the frozen look and `respec`, and that every frame has the same pixel size
  and hash;
- checks a resumable backfill that skips the blind week;
- checks the order of the concat lists, ffmpeg's arguments (a fake ffmpeg), and ffmpeg missing;
- makes one real tiny MP4 and counts its frames with ffprobe (when ffmpeg is installed), and one of mixed
  frame sizes (a storm frame from before a `respec`) that must keep every frame;
- checks that a failing ffmpeg is a NOTE (exit 0), that odd sizes such as 854x480 come out exact on
  matplotlib 3.3, that a run by hand applies `waterline_station.conf`, that missed weeks are capped at 4 a
  daily run, that re-processed rows rebuild an open week, and when a failed build is retried;
- checks the locks, and that every write lands under `archive/weekly/` or `logs/`.

# The waterline maps at a new station

This repository was written for CACO-05 (Marconi Beach), and its built-in values are CACO-05's. Another
station runs the same code with its own settings, kept in one file beside the cron:
`/mnt/I2Rgus_Data/waterline/waterline_station.conf`. coastcam-station-setup writes that file on a new
NUC (`./setup.sh`, then `./setup.sh --waterline`). This page covers:

- what the file is and what each key does (sections 1-2);
- what stays provisional at a new station, and why (section 3);
- the two jobs left to do by hand: refitting the wave setup coefficient C (section 4) and tuning the
  detector with hand-clicked ground truth (section 5);
- the checks setup runs, `waterline_station_check.py` (section 6);
- what never to run in the NUC's git clone (section 7).

**CACO-05 has no `waterline_station.conf`.** Without the file the cron runs the same commands and writes
the same files and logs as before. `test_station_overlay.py` checks this against the cron as it was
before the station file existed (docs/VALIDATION.md section 2.11); the one step added since on purpose, the
weekly frames (`weekly_frames.py daily`, step 7c), is the only difference it allows. Never run
coastcam-station-setup on CACO-05.

The storm watch has its own station file and an alerts mode; see `docs/STORM_WATCH_ALERTS_MODE.md`.

---

## 1. The station file

### How the cron reads it

On each run, `waterline_timex_cron.sh` does this, in order:

1. Sets every CACO-05 value (its `NAME=` lines), makes `logs/` and defines `log()`.
2. Clears `WATERLINE_STAGE`, `CAMS`, `OWG_ENABLE`, `DETECTOR_ARGS` and the four `CAL_IO_*`/`CAL_EO_*`
   names from the environment the job was started with. Only the file can set them.
3. If `waterline_station.conf` exists, runs
   `python3 waterline_station_check.py check-overlay --conf waterline_station.conf` on it:
   - The check's summary goes to `logs/timex_cron.log`.
   - If the check refuses the file, the run logs
     `FATAL: /mnt/I2Rgus_Data/waterline/waterline_station.conf refused (see above): nothing run ...`,
     runs nothing else and exits 1.
4. Otherwise it sources the file, so the file's values replace CACO-05's.
5. Sets the defaults of the names the file may leave out, for example
   `: "${CAMS=c1 c2}"`. These default lines are not `NAME=` lines. That matters because five other
   scripts read the cron's `VAR=` lines (survey_products.py, recalibrate_timing.py, storm_watch.py,
   station_status.py, dem_from_contours.py), and they read the same values as before.

The file is shell, so `check-overlay` accepts only these lines:

- comments (`# ...`) and blank lines;
- `KEY="value"`:
  - It starts in column 0, with no spaces around `=`.
  - The value is in double quotes, with nothing after the closing quote.
  - Each key appears once.
  - Inside the quotes, anything is allowed except `$`, a backtick, `\`, `"` and control characters.
    So `;`, `(`, `)` and `=` are fine: the shell expands nothing there.
- one line `export NAME NAME ...`, naming only the exported keys of section 2.

The check also refuses these:

- an unknown or forbidden key;
- a value that does not fit its key;
- Windows line ends;
- one-camera calibration names that do not agree with `CAMS`.

A stage-2 file (one without `WATERLINE_STAGE="archive"`) must also:

- name every station fact. Otherwise the cron would quietly use CACO-05's value: its calibration,
  GNSS-R file, gauge, buoy, antenna height and so on.
- list on the export line every exported key it sets.

### Who else reads it

| Reader | What it takes from the file |
|---|---|
| `waterline_timex_cron.sh` | Everything (above). |
| `storm_watch.sh` | Runs `check-overlay`, then sources the file, so the exports reach the storm DEM builds. A refused file stops the storm watch with an ERROR line in `logs/storm_watch.log`. |
| `storm_watch.py` | Only `DEM_MAX_HS`, applied after the cron's lines. |
| `dem_from_contours.py --fit-setup` | Its `SETUP_COEF`, as the C in use (section 4). |
| Every Python step the cron runs | The exported keys, from the environment (`station_env.py`). |

### Tools that never read it

These are CACO-05's tools and never read the station file, so they give CACO-05 answers:

- `survey_products.py` and `historical_forcing.py`: CACO-05's survey dates and forcing, with their own
  C. Their tables are in `surveys/`, which a new NUC does not have.
- `recalibrate_timing.py`: the `TIMING_MODE` switch. It reads the GNSS-R file, buoy and gauge from the
  cron's own lines, which are CACO-05's.

Do not use them at a new station. `TIMING_MODE` stays `legacy` there, and the file may not set it.

### Two stages

| Stage | When | The file | What a cron run does |
|---|---|---|---|
| 1 | From install day until the calibration, OPUS height and GNSS water levels are in | `WATERLINE_STAGE="archive"`, `OWG_ENABLE="0"`, `ARCHIVE_PRODUCTS=""` and `STATION_NAME` | Copies the new timex pictures into `archive/images_timex` (never overwriting one), then logs `stage 1 (WATERLINE_STAGE=archive): timex kept: N; ...` and `=== run end ===`. No detection, maps or DEM; no Python but the check. |
| 2 | After `./setup.sh --waterline` (setup's "continue" step) | The full file (below), without `WATERLINE_STAGE` | The whole pipeline, with the station's values. |

The products of a stage-2 run, in `/mnt/I2Rgus_Data/waterline`: the waterline points
(`contour_points_timex.csv`, `contour_points_ground.csv`), the 7- and 30-day elevation maps per camera
(`elevation_map_<cam>_{7,30}day.png`), the whole-archive and 7-day DEMs (`dem_intertidal*`, with the 7-day
DEM's map and GeoTIFF and its copy in `archive/dems/`), `station_status.txt` (no runup at a new station),
and the
**weekly frames** (`archive/weekly/`, docs/WEEKLY_FRAMES.md): every week (Monday to Sunday, UTC dates) a new
DEM frame and waterline frame in one frozen look, and the two MP4 videos made from the final weeks
(`archive/weekly/video/`; they need ffmpeg). The storm watch in full mode adds its storm frames to those
videos (docs/STORM_WATCH_ALERTS_MODE.md). The weekly step runs only while `DEM_ENABLE=1`; in stage 1 the
cron stops before it.

Setup's full file (the CACO07 names and numbers are examples):

```
# waterline_station.conf -- CACO07's settings for waterline_timex_cron.sh.
# Written by coastcam-station-setup (./setup.sh --waterline); not tracked by git.
# Do not edit: change station.conf and run ./setup.sh --waterline again.
CAMS="c1 c2"
CAL="/mnt/I2Rgus_Data/waterline/station_calibration"
CAL_IO_C1="CACO07_c1_20261201_IO.yaml"
CAL_EO_C1="CACO07_c1_20261201_EO.yaml"
CAL_IO_C2="CACO07_c2_20261201_IO.yaml"
CAL_EO_C2="CACO07_c2_20261201_EO.yaml"
GNSSR_SPLINE="/home/argus_user/GNSS/v4.1/products/refl_code/Files/caco/caco_spline_out.txt"
TIDE_GAUGE="8447435"
GAUGE_CSV="/mnt/I2Rgus_Data/waterline/archive/gauge_8447435.csv"
WAVE_BUOY="44008"
WAVES_CSV="/mnt/I2Rgus_Data/waterline/archive/waves_44008.csv"
DEM_MAX_HS="1.0"
SETUP_COEF=""
DETECTOR_ARGS="--no-bias-correction --profile-json /mnt/I2Rgus_Data/waterline/station_calibration/detector_profile.json"
OWG_ENABLE="0"
RUNUP_ENABLE="0"
POINTING_LOG_ENABLE="0"
ARCHIVE_PRODUCTS=""
RAS_SOURCE=""
GNSSR_ANTENNA_NAVD88_M="12.345"
STATION_TZ="America/New_York"
GEOTIFF_EPSG="6348"
STATION_NAME="CACO07"
STATION_PROVISIONAL="search band from camera geometry; no wave setup (C = 0)"
export GNSSR_ANTENNA_NAVD88_M STATION_TZ GEOTIFF_EPSG STATION_NAME STATION_PROVISIONAL STATION_WAVE_CURRENCY GNSSR_QC_GROSS_M WAVE_BUOY
```

For a one-camera station: `CAMS="c1"`, `CAL_IO_C2=""` and `CAL_EO_C2=""`. Both empty values must be
written: a missing line would make the cron look for CACO-05's c2 files.

### Editing it by hand

The file is setup's: `./setup.sh --waterline` writes it again in full, and setup's final check is
designed to flag a file that differs from what setup would write. The two jobs in sections 4 and 5 do
change it by hand, so before any later `./setup.sh --waterline`:

- keep a copy of the lines you changed;
- compare them with the file again after setup runs.

After every edit, check the file before the next cron run (12:30, 19:25 or 20:55):

```
cd /mnt/I2Rgus_Data/waterline
python3 waterline_station_check.py check-overlay --conf waterline_station.conf > /dev/null
```

It ends `check-overlay: waterline_station.conf OK (stage 2: 24 setting(s), 8 exported)` and exits 0.
Otherwise it prints `REFUSED (N problem(s)):` with one line per problem and exits 1, for example:

```
  line 15: TIMING_MODE may not be set here: the timing mode is the cron's one switch (...)
  line 16: the value must be in double quotes with nothing after them: SETUP_COEF="value"
  file: missing (CACO-05's value would be used without a word): TIDE_GAUGE ...
```

The file is git-ignored, like `storm_watch_station.conf`, `station_calibration/` and
`.coastcam-station-setup` (setup's marker), so `git pull` never touches them.

---

## 2. Every key

"CACO-05" is the value the cron uses when the file does not set the key (and on CACO-05 itself).
"Setup" is what coastcam-station-setup writes at a new station.

**Stage, cameras and calibration**

| Key | What it does | CACO-05 | Setup |
|---|---|---|---|
| `WATERLINE_STAGE` | `"archive"`: stage 1, timex pictures kept and nothing else run. `""`: the full run. | `""` | `"archive"` in stage 1; not written in stage 2 |
| `CAMS` | The cameras mapped (`elevation_map_<cam>_*.png`): `"c1"` or `"c1 c2"`. | `c1 c2` | from the calibration |
| `CAL` | The folder of the calibration files. An absolute path, no spaces. | `$BASE/calibration` | `$BASE/station_calibration` |
| `CAL_IO_C1`, `CAL_EO_C1` | c1's lens (IO) and position (EO) file names in `CAL`. Used by both georectify calls and the calibration check; the EO also by the pointing log. A file name only, no folder. | `CACO05_c1_20240801_IO.yaml`, `CACO05_c1_20251113_EO-CV.yaml` | setup's normalised copies |
| `CAL_IO_C2`, `CAL_EO_C2` | The same for c2. `""` on a one-camera station: georectify then gets no c2 options, and the calibration check needs only c1. | `CACO05_c2_20240801_IO.yaml`, `CACO05_c2_20251113_EO-CV.yaml` | the c2 copies, or `""` |

**Water levels and waves**

| Key | What it does | CACO-05 | Setup |
|---|---|---|---|
| `GNSSR_SPLINE` | GPS_code's gnssrefl spline: the measured water level every waterline is given. | `.../refl_code/Files/usgs/usgs_spline_out.txt` | the station's own spline |
| `TIDE_GAUGE` | The NOAA CO-OPS gauge (7 digits) that GNSS-R QC checks against. Fetched every run. | `8447435` (Chatham) | the nearest NAVD88 gauge |
| `GAUGE_CSV` | Where that gauge's record is kept. | `archive/gauge_8447435.csv` | `archive/gauge_<gauge>.csv` |
| `WAVE_BUOY` | The NDBC buoy fetched every run. Exported, so `gnssr_qc.py` finds `waves_<buoy>.csv`. | `44008` | from setup's lookup |
| `WAVES_CSV` | Where that buoy's record is kept. Each frame is tagged with its wave height and period. On CACO-05 the cron uses owg.sh's `waves_marconi.csv` instead whenever that file exists; a new station has no such file, so it uses this raw offshore buoy record. | `archive/waves_44008.csv` | `archive/waves_<buoy>.csv` |
| `DEM_MAX_HS` | The DEM leaves out frames whose offshore wave height is above this many metres. `""` keeps every frame. The storm watch's DEMs use it too. | `1.5` | `1.0` |
| `SETUP_COEF` | The wave setup coefficient C: each waterline's level becomes the water level + C x sqrt(Hs x L0). `""` means no setup (C = 0). Section 4. | `0.037` | `""` |

**Switches**

| Key | What it does | CACO-05 | Setup |
|---|---|---|---|
| `DETECTOR_ARGS` | Added to the detector's options, word-split. It may not repeat an option the cron sets itself (`--image-suffix`, `--timing-mode`, `--input-dir`, `--output-dir`, `--debug-dir`, `--source-dir`, `--help`, `-h`, or a short form such as `--out`). It may not hold wildcards or `~`. | `""` | `--no-bias-correction --profile-json $BASE/station_calibration/detector_profile.json` |
| `OWG_ENABLE` | `1` runs owg.sh (CACO-05's optical wave gauge, its models and its email list) in step 5c. A new NUC's clone has no owg.sh. | `1` | `0` |
| `RUNUP_ENABLE` | `1` runs the runup from c2 timestacks along CACO-05's lines. | `1` | `0` |
| `POINTING_LOG_ENABLE` | `1` runs the daily camera pointing log (`pointing_daily.py`). Keep it `0` on a one-camera station: pointing_daily.py would still try c2 with CACO-05's EO. | `1` | `0` |
| `ARCHIVE_PRODUCTS` | The other Argus products kept in `archive/images_<product>` (about 4.5 GB a year per product per camera). | `bright dark snap var` | `""` |
| `RAS_SOURCE` | Where the c2 timestacks are copied from (about 100 GB a year). `""` keeps none. | `/mnt/I2Rgus_Data/ImageProducts/products` | `""` |

**Exported: read by the Python steps from the environment** (`station_env.py`, which checks each value)

| Key | Read by | CACO-05 (when unset) | Setup |
|---|---|---|---|
| `GNSSR_ANTENNA_NAVD88_M` | `extract_elevation_contours.py`, and `compare_gnssr_to_tidemodel.py` and `gnssir_reflection_audit.py` through it: the GNSS antenna's height, m NAVD88, which turns reflector heights into water levels. An empty or non-numeric value stops the script with an ERROR; it never falls back to 19.014, which would put another station's levels metres off. | `19.014` | the OPUS (GEOID18) height |
| `STATION_TZ` | `extract_elevation_contours.py`, `dem_from_contours.py`: which local day a frame belongs to. An IANA name. | `America/New_York` | the station's zone |
| `GEOTIFF_EPSG` | `asc_to_geotiff.py` (the GeoTIFF's tag; `--epsg` still wins), and the axis labels of `dem_map.py` and `dem_change.py`. NAD83(2011) UTM zone N is EPSG 6329+N (zones 1-19). | `6348` (UTM 19N) | the calibration's zone |
| `STATION_NAME` | The titles of the DEM page (`dem_figure.py`) and the map (`dem_map.py`), and the weekly frames' titles (frozen in `archive/weekly/frame_spec.json` when it is made; a later change needs `weekly_frames.py respec`). | the built-in CACO-05 titles | the station's name |
| `STATION_PROVISIONAL` | Free text: what is provisional about the products. It goes into `<dem>_info.json` (`station_provisional`) and appears as a `PROVISIONAL: ...` line under the DEM page's heading. Section 3. | none | what is provisional |
| `STATION_WAVE_CURRENCY` | The wave record C was fitted with, as `<kind>:<height>/<period>`. With a setup on, frames tagged with raw buoy waves that are not this record get a WARNING on every run. Section 4. | `marconi:hs_best/tp_s` | not set while C = 0 |
| `GNSSR_QC_GROSS_M` | `gnssr_qc.py`'s gross range for a GNSS-R level, `"lo,hi"` in m NAVD88. | `-2.5,3.5` | set only when the gauge's MLLW or MHHW falls outside -2.5 to +3.5 m |
| `WAVE_BUOY` | (above) `gnssr_qc.py`'s default waves file. | `44008` | |

**Refused**, with a reason:

- **The cron's own settings.** The storm watch, the survey products and the station status read them
  from the cron itself, so a station value would split the station in two: `TIMING_MODE`, `DEM_ENABLE`,
  `DEM_CELL`, `DEM_MIN_POINTS`, `DEM_MAX_SPREAD`, `DEM_MAX_DAY_OFFSET`, `DEM_WINDOW_DAYS`, `GROUND`,
  `LOCK`, `GNSSR_START_FILE`.
- **The cron's folders.** `BASE`, `LOG`, `SCRATCH`, `SCRATCH_INPUT`, `DEBUG_DIR`, `ARCHIVE_CSV`,
  `ARCHIVE_IMG`, `ARCHIVE_RAS`.
- **Shell and loader names.** `PATH`, `IFS`, `PYTHONPATH`, `PYTHONHOME`, `LD_PRELOAD`,
  `LD_LIBRARY_PATH`, `BASH_ENV`, `ENV`, `SHELLOPTS`.

Any other name not listed above is refused too. Settings such as `CONSISTENCY_ENABLE`, `SETUP_MODEL`,
`SETUP_SOURCE` and `DEBUG_KEEP_DAYS` therefore keep CACO-05's values at every station.

---

## 3. What is provisional, and why

A new station's first maps are made before anything has been checked against its own beach.
`STATION_PROVISIONAL` says so on every DEM page. Setup writes, for example,
`search band from camera geometry; no wave setup (C = 0)`.

| What | Why it is provisional | Its effect | What ends it |
|---|---|---|---|
| **The detector's search band** | It comes from `waterline_station_check.py profile`: from a beach survey or DEM if one was uploaded; otherwise from the camera geometry alone, as a level plane 40-600 m from the camera. It is not from hand-clicked waterlines. | A band from geometry is wide, so the detector has more room to lock onto other edges (a wet-sand line, a dune toe). | Ground truth, then a new envelope (section 5). |
| **No bias correction** | `--no-bias-correction`, plus `"no_bias": true` in the profile JSON. CACO-05's per-column corrections belong to its own cameras and beach. | Any steady offset between the detected and the true waterline goes straight into the elevations. | Ground truth from 2 or more days, then `derive_bias_correction.py` (section 5). |
| **CACO-05's other detector tunings** | The profile JSON replaces only crop, envelope, bias, step and profile_scale. Everything else stays CACO05's: band_contrast_floor, the snap settings, the confidence and partial-line thresholds, the temporal blend. | Tuned on Marconi's 2026 imagery; unknown elsewhere. | Only an edit to the detector. That is outside this repository's station file, and the detector must stay byte-identical so `git pull` keeps working. |
| **No wave setup (C = 0)** | C depends on the beach slope and the wave record; nobody has fitted it here. | In waves, a timex shows the water's edge above still water, so the DEM reads low. At Marconi's C (0.037) the setup would be about 0.26 m at Hs 0.5 m and 0.37 m at Hs 1.0 m (Tp 8 s). This is why `DEM_MAX_HS` is 1.0 m, not CACO-05's 1.5 m. | Refitting C (section 4). |
| **Raw offshore buoy waves** | A new station has no wave transfer to the beach (CACO-05's `waves_marconi.csv`, from owg.sh). | The `DEM_MAX_HS` filter, and later the setup, use offshore heights, which may overstate sheltered seas. | A local wave record (not built here). |
| **The water level** | GNSS-R from the OPUS antenna height (`GNSSR_ANTENNA_NAVD88_M`), quality-checked against the gauge. | An error in the OPUS height moves every level and every DEM by the same amount. | A comparison with a survey, below. |
| **No camera pointing log** | `POINTING_LOG_ENABLE="0"`: there is no reference bank for this camera yet. | A camera that turns by a degree or two moves the lines without any warning. `waterline_station_check.py pointing` checks it once at setup, when a calibration-day picture was uploaded. | Not in this repository's station flow (section 8). |
| **No comparison with a survey** | Nothing has measured the DEM against the beach yet. | The DEM's error is unknown. | An RTK beach survey down to the low-tide line, compared with `survey_compare.py` (below). |

**Comparing with a survey.** `survey_compare.py --dem dem_intertidal_7day_dem.asc --survey <rtk.csv>
--survey-type points --label ... --why "..." --survey-date ... --epsg <GEOTIFF_EPSG>` prints, next to every
number, how independent the comparison is. Pass `--epsg`: its default is CACO-05's 6348.

**When the maps stop being provisional.** That is the station owner's decision. The basis is the
comparison with an RTK survey, after the tuning of section 5 and the refit of section 4. Then:

- shorten `STATION_PROVISIONAL` to what is still provisional, or set it to `""`;
- keep the name on the export line;
- run `check-overlay`.

The next run's DEM page has no PROVISIONAL line.

---

## 4. Refitting the wave setup coefficient C

**Before the refit.** Wait for several weeks of waterlines over different waves. One week's curve is
noisy and can sit near 0. With `SETUP_COEF=""` the waterlines carry no setup, which is what
`--fit-setup` needs.

1. **Run the fit.** The fit compares repeat crossings of each 2 m cell for every C. Write its outputs
   under `/tmp`, never over the station's DEM:

   ```
   cd /mnt/I2Rgus_Data/waterline
   python3 dem_from_contours.py contour_points_ground.csv /tmp/fit_setup --fit-setup \
       --setup-coef-in-use 0 --max-hs 1.0 --no-plot
   ```

   - `--setup-coef-in-use 0` is also what it reads from the station file's `SETUP_COEF=""`.
   - `--max-hs` is the station's `DEM_MAX_HS`.
   - Add `--start-date YYYY-MM-DD --end-date YYYY-MM-DD` (UTC dates) to fit one window.

2. **Read the output.** It prints the curve, the near-best band (each C whose median spread is within
   5 mm of the best) and a verdict:
   - It offers a C to apply only when C = 0 and the C in use both lie outside the band, and the
     minimum is not flat.
   - Adopt a C only when it holds over several windows of different waves.
   - `--fit-setup` does not re-project the lines, so it is in principle biased a little low
     (docs/VALIDATION.md section 4.1).
   - An RTK survey down to the low-tide line, on days with different waves, pins C more closely.
     This repository's survey tools for that (`survey_products.py`, `evaluate_setup_model.py`) are set
     up for CACO-05's survey dates, so at a new station that is a separate analysis.

3. **Set it** in `waterline_station.conf`:

   ```
   SETUP_COEF="0.030"
   STATION_WAVE_CURRENCY="buoy:wvht_m/dpd_s"
   ```

   - `STATION_WAVE_CURRENCY` names the record C was fitted with. At a new station that is the raw
     buoy file the cron keeps (`WAVES_CSV`, written by `fetch_buoy_waves.py`), whose tag is
     `buoy:wvht_m/dpd_s`.
   - Without it, every run warns that the frames are not in the record C was fitted with (CACO-05's
     `marconi:hs_best/tp_s`).
   - The name must be on the export line; setup's line already lists it.
   - Shorten `STATION_PROVISIONAL` (drop "no wave setup (C = 0)").
   - With a setup on, `DEM_MAX_HS` may be raised: CACO-05 uses 1.5 m with its C.
   - Then run `check-overlay`.

**What changes on the next run.** It rebuilds the contours, the ground points and the DEMs from the
whole archive with the new C, so every DEM changes. The week-to-week change map rebuilds both weeks
with the current C (`dem_change.py --rebuild`), so the new C does not show there as beach change.

---

## 5. Tuning the detector with hand-clicked ground truth, then detecting again

**What this replaces.** Setup's search band comes from geometry or a survey, with no bias correction
(section 3). Hand-clicked waterlines replace both:

- **the envelope**, from the clicked rows with a margin. Needs 4 or more frames, over different tides.
- **the bias correction**, from the detector's error against the clicks. Needs 2 or more days of
  clicks.

**Where the results go.** Into the station's own files: `station_calibration/detector_profile.json`
and `DETECTOR_ARGS`. Never into `waterline_detector_v5.py`. It must stay byte-identical, or `git pull`
stops. The station's camera is not one of its built-in CACO05 profiles anyway.

**The tools are CACO-05's.** Their folders are fixed: `/mnt/I2Rgus_Data/waterline/input` for pictures to
click and `/mnt/I2Rgus_Data/waterline/ground_truth` for the clicks.
- `input/` is the detector's default input folder, which the cron never uses (its own is
  `input_timex/`).
- `ground_truth/` is one of the folders a new NUC's clone leaves out (section 7). The station's own
  clicks there are untracked files, which git leaves alone.

Never run the detector by hand without `--input-dir`/`--output-dir`: it clears its input and output
folders.

1. **Choose and copy frames.** Daylight timex from `archive/images_timex`, over several days and the
   whole tidal range, for one camera at a time. Copy them; never move them:

   ```
   cd /mnt/I2Rgus_Data/waterline
   mkdir -p input
   cp -p archive/images_timex/<epoch>.*.c1.timex.jpg input/      # one line per chosen frame
   ```

2. **Click.** This needs a display: `ssh -X`, or MobaXterm from Windows.

   ```
   python3 collect_all_ground_truth.py --filter .c1.
   ```

   Click along the true water's edge, then use the on-screen buttons (Save & Next / Skip / Undo / Quit).
   - Each saved frame becomes `ground_truth/<picture name>.csv`: one row per picture column.
   - Frames already clicked are skipped, so a session can stop and resume.

3. **The envelope.** Run:

   ```
   python3 derive_envelope.py --image-height 2048 c1        # 2048: the picture height (the lens's NV)
   ```

   **Do not paste its "Ready-to-paste" block.** Those fractions are relative to CACO05's own crop.
   The station's profile uses the whole picture (`"crop": [0, 1, 0, 1]`). Instead:
   - Take the `min_row(full)` and `max_row(full)` columns, which already include the margin. Divide
     each by the picture height.
   - Write the six `[x_frac, min/height, max/height]` triples as that camera's `"envelope"` in
     `station_calibration/detector_profile.json`. Leave its `"crop"`, `"step"` and `"profile_scale"`.
   - Keep a copy of the file first. Setup's `profile` command wrote it.
   - With fewer than 4 frames the script warns that the envelope has no spread: click more.

4. **The bias correction** (after 2 or more days of clicks). Run it on the cron's own detections,
   which were made with no correction:

   ```
   python3 derive_bias_correction.py --uncorrected --processed-dir archive/processed_timex --suffix .timex. c1
   ```

   - `--uncorrected` is required here. Without it, CACO-05's built-in correction is added to an error
     that never had it subtracted.
   - Write the table's `x_frac` and `NEW corr` columns (px, subtracted from the detected row) as
     `station_calibration/detector_bias.json`, for example
     `{"c1": [[0.05, 3.1], [0.15, 2.4], ...], "c2": [...]}`.

   Then give the detector the file, and drop `"no_bias": true` from that camera's entry in the
   profile JSON. The detector applies `--profile-json` after `--bias-correction-file`, so
   `"no_bias": true` would clear the new correction again. Then set:

   ```
   DETECTOR_ARGS="--no-bias-correction --bias-correction-file /mnt/I2Rgus_Data/waterline/station_calibration/detector_bias.json --profile-json /mnt/I2Rgus_Data/waterline/station_calibration/detector_profile.json"
   ```

   Run `check-overlay` after the edit.

5. **Try it before it goes live.** Run `waterline_station_check.py trial` (section 6) on 5 kept frames,
   with the new profile. It runs in scratch folders and changes nothing of the station's. If you use
   a bias file, the trial cannot apply it: it runs the detector with `--no-bias-correction
   --profile-json` only.

6. **Detect the archive again.** The archive never overwrites a detection (`cp -pn`), so new settings
   reach only new pictures. To apply them to every kept picture, the old detections are moved aside
   and `archive/images_timex` is detected again. coastcam-station-setup's
   `tools/waterline_setup.sh redetect` does this. By hand, the same steps, all under the cron's lock
   (a cron run meanwhile logs `SKIPPED`; about 25 minutes per 1,800 pictures):

   ```
   cd /mnt/I2Rgus_Data/waterline
   flock .timex_cron.lock sh -c '
     mv archive/processed_timex archive/processed_timex.before_tuning_$(date -u +%Y%m%d) &&
     mkdir -p archive/processed_timex work/redetect/in work/redetect/out work/redetect/debug &&
     python3 waterline_detector_v5.py --image-suffix timex.jpg --timing-mode legacy \
         --source-dir archive/images_timex --input-dir work/redetect/in \
         --output-dir work/redetect/out --debug-dir work/redetect/debug \
         --no-bias-correction --profile-json station_calibration/detector_profile.json
     cp -pn work/redetect/out/*.csv archive/processed_timex/'
   rm -rf work/redetect          # its debug images: about 8 MB per frame
   ```

   - Use exactly the options of the station file's `DETECTOR_ARGS`, including `--bias-correction-file`
     when you set one.
   - The next cron run rebuilds the contours, maps and DEMs from the re-detected archive.
   - Earlier dated DEMs in `archive/dems` stay as they were built.
   - Shorten `STATION_PROVISIONAL` (drop "search band from camera geometry").

---

## 6. waterline_station_check.py

This is coastcam-station-setup's interface to this checkout.

**Output.** Every command that reports a result prints:
- JSON on stdout, which setup reads. NaN and infinities are written as `null`.
- a summary for people on stderr.

**Exit status:** 0 ok, 1 refused or failed, 2 usage (a missing file or a wrong option).

**Python.** `capabilities` and `check-overlay` use the standard library only, because the cron runs the
check on every run. The other commands load numpy and this repository's own modules (unchanged) only
when they run: `profile` also needs matplotlib, `pointing` OpenCV.

Examples are run from `/mnt/I2Rgus_Data/waterline`. Setup keeps its normalised copies of the uploaded
files in `work/setup_candidate/`.

### capabilities

```
python3 waterline_station_check.py capabilities
overlay-1 storm-alerts-1 station-setup-1
```

Setup reads these tags to tell an old clone (one to pull) from one that has the commands below.

### check-overlay

```
python3 waterline_station_check.py check-overlay --conf waterline_station.conf
```

- Section 1 lists the rules.
- JSON fields: `ok`, `stage` (`archive` or `full`), `keys` (the names set), `exported`,
  `errors` (`line`, `key`, `message`) and `notes`.
- An unreadable or non-UTF-8 file is refused.

### pair: which lens, position, targets and picks belong together

```
python3 waterline_station_check.py pair --xyz work/setup_candidate/targets_xyz.csv \
    --uv work/setup_candidate/c1_uv.csv \
    --io work/setup_candidate/lens_A.yaml work/setup_candidate/lens_B.yaml \
    --eo work/setup_candidate/pos_1.yaml work/setup_candidate/pos_2.yaml --image 2448x2048
```

**Inputs.** One target-coordinate file and one picks file per call; setup loops over the pairs.
- `--xyz`: header-less `num,E,N,Z`, in UTM metres and m NAVD88.
- `--uv`: `num,U,V` with that header.

**How a set is chosen.**
- Every lens x position combination goes through georectify's own self-test, unchanged. It passes
  with a mean error under 1.0 m.
- A set is chosen when it passes and the next best *position* is at least 3 times worse.
- Two lenses that pass with that position within 20% of each other are an `io_tie` (exit 0,
  `chosen.io_decided` false), for setup to settle by a name hint or by asking. Cameras with the same
  lens model are never refused.
- Two positions that both pass without the ratio are refused, with an `eo_tie` list.
- Fewer than 4 shared targets: not paired. 4 or 5 give a warning.

**Then it checks:**
- the lens's NU x NV against `--image`;
- the picks inside the frame;
- the targets 10-1000 m from the camera, at -5 m to camera z + 5 m.

**JSON fields:** `combinations`, `chosen` (`"none"` or the set, with the camera's position and angles
and the lens's size and focal lengths), `io_tie`, `eo_tie`, `ratio_next_eo`, `next_eo`, `sanity`,
`reason`.

On CACO-05's own 2025-11-13 c1 files, the right set gives a mean of 0.33 m, and the next best position
is 122 times worse:

```
pair: gcp_uv_c1.csv + 2025-11-13_Marconi_Extrinsic_Targets_c1_xyz.csv: 14 target(s) shared, 4 lens x position combination(s)
  CACO05_c1_20240801_IO.yaml + CACO05_c1_20251113_EO-CV.yaml: n=14 mean 0.33 m max 1.181 m PASS
  CACO05_c2_20240801_IO.yaml + CACO05_c1_20251113_EO-CV.yaml: n=14 mean 6.36 m max 9.487 m FAIL
  ...
  chosen: CACO05_c1_20240801_IO.yaml + CACO05_c1_20251113_EO-CV.yaml, mean 0.33 m; next best position CACO05_c2_20251113_EO-CV.yaml 122.4 times worse
  ok   image_size: lens file 2448 x 2048, pictures 2448 x 2048
  ok   target_z: target Z 0.73 to 3.98 m (allowed -5 to 25.53 m NAVD88: camera z + 5)
pair: OK
```

### profile: the detector's search band

```
python3 waterline_station_check.py profile --camera c1 --io $CAL/CACO07_c1_20261201_IO.yaml \
    --eo $CAL/CACO07_c1_20261201_EO.yaml --image-size 2448x2048 \
    --survey-csv beach_survey.csv --mllw -0.95 --mhhw 1.10 --output work/setup_candidate/detector_profile.json
python3 waterline_station_check.py profile --camera c2 ... --geometry --plane-z 0.05 \
    --output work/setup_candidate/detector_profile.json --merge
```

Here `$CAL` stands for the station's calibration folder.

**Choose one source:**
- `--survey-csv`: E,N,Z points, `num,E,N,Z`, or a header naming them. Triangulated onto a 0.5 m
  grid; there is no surface across gaps over 100 m.
- `--dem`: a beach DEM, `.tif` or `.asc`, in the camera's UTM zone and NAVD88.
- `--geometry`: a level plane at `--plane-z` (the gauge's MSL, default 0) 40-600 m from the camera.
  The result is PROVISIONAL.

**The band.** With a survey or DEM it runs from MLLW - 0.5 m to MHHW + 1.0 m (else -1.5 to +2.5 m),
within 300 m.

**The file.** The band is placed by detect_original_view's own `envelope_from_survey`. The file is
written as detect_original_view writes it, as
`{"c1": {"crop": [0,1,0,1], "envelope": [...], "no_bias": true}}`, plus `step` and `profile_scale`
when the far side is on the left. `--merge` keeps the other camera's entry.

**Refusals.** When the band hardly shows in the picture, or a DEM does not cover the camera's
surroundings, it refuses with exit 1 and a reason.

**JSON fields:**
- `status`: `OK` or `PROVISIONAL`;
- `provisional_reason`;
- `note`: which CACO05 tunings stay;
- `station_provisional`: a text `check-overlay` accepts.

### trial: the whole chain on a few kept pictures, in scratch folders

```
python3 waterline_station_check.py trial --camera c1 --io $CAL/CACO07_c1_20261201_IO.yaml \
    --eo $CAL/CACO07_c1_20261201_EO.yaml --profile-json work/setup_candidate/detector_profile.json \
    --spline /home/argus_user/GNSS/v4.1/products/refl_code/Files/caco/caco_spline_out.txt \
    --frames archive/images_timex/<epoch1>.*.c1.timex.jpg ... (5 at most) \
    --work-dir work/trial --antenna-navd88 12.345 \
    --waves archive/waves_44008.csv --gauge-csv archive/gauge_8447435.csv
```

**What it runs.** The pictures are copied into `<work-dir>/trial_c1/src`. Then:
1. the detector, with `--no-bias-correction --profile-json`;
2. `extract_elevation_contours.py`, with `GNSSR_ANTENNA_NAVD88_M` set to `--antenna-navd88`;
3. `georectify.py`, with this camera's calibration.

Every step writes into `trial_c1/`, and the steps' output goes to `trial_c1/trial.log`.

**The work folder.**
- It may not be inside the station's own folders (`input`, `input_timex`, `processed`,
  `processed_timex`, `archive`, `debug`, `calibration`, `station_calibration`) or ImageProducts.
  Such a folder is refused with exit 2.
- `trial_c1/` is replaced on each run, but only if it carries the trial's marker file.

**The result.**
- A frame passes with signal in at least 50% of the columns, 80% or more of its points placed, and
  all of them 30-1000 m from the camera.
- Exit 0 when one frame passes.
- JSON fields: `frames` (each with its reasons), `passed`, `best`, `steps` (each step's exit status),
  `work_dir`, `log`.

### pointing: has the camera turned since the calibration?

```
python3 waterline_station_check.py pointing --io $CAL/CACO07_c1_20261201_IO.yaml \
    --eo $CAL/CACO07_c1_20261201_EO.yaml \
    --calibration-picture ~/waterline_upload/calibration/c1_calibration_day.jpg \
    --timex archive/images_timex/<epoch1>.*.c1.timex.jpg archive/images_timex/<epoch2>.*.c1.timex.jpg
```

- It compares a picture taken on the calibration day with recent timex of the same camera, using
  pointing_check's matcher, unchanged.
- The verdict is from the median rotation of the matched frames: PASS under 0.3 deg, WARN 0.3-1.0 deg
  (exit 0 for both), FAIL over 1.0 deg (exit 1).
- These also exit 1, as JSON with a reason:
  - `UNMATCHED`: too few land features matched; nothing was checked.
  - `ERROR`: a picture of another size than the lens, or one that cannot be read.
  - `UNAVAILABLE`: no OpenCV.
- The trial's 30-1000 m range check cannot see a turn of 1-2 deg; this can.

---

## 7. The git clone on a new NUC: never fetch what it leaves out

**What the clone is.** coastcam-station-setup clones this repository into `/mnt/I2Rgus_Data/waterline`
as a sparse partial clone (`git clone --filter=blob:none`, plus sparse-checkout patterns). It checks
out every top-level file and folder except:

- `calibration/`, `ground_truth/`, `surveys/`, `owg_models/`;
- `owg.sh`;
- `manifest_*.csv`;
- `sig1000_waves_ALL.csv`.

The contents of those files are never downloaded. That keeps CACO-05's calibration and survey files,
and owg.sh's recipient list, off every other NUC.

**Why some commands are dangerous.** A partial clone downloads a file's contents the first time a git
command needs them, and keeps them in `.git`. So on a NUC, never run:

- `git log -p` or `git log --stat`, `git show <commit>` or `git show <commit>:<path>`, or
  `git diff <commit> <commit>` over the left-out paths;
- `git blame`, `git grep <commit>`, `git archive`, or `git checkout <commit> -- <path>` /
  `git restore --source=<commit> <path>` on them;
- `git sparse-checkout disable`, or a `git sparse-checkout set`/`add` that takes the left-out paths
  back in. This downloads and checks out all of them at once.

**What is safe.** `git status`, `git pull`, `git log --oneline` and a plain `git diff` (your own edits
to checked-out files).

**Repairing it.** If one of the commands above was run, the left-out files are in `.git` (and, after
`sparse-checkout disable`, in the folder too). Setting the patterns again removes them from the folder
but not from `.git`. Only a fresh clone does; tell the station's owner.

**The station's own files.** The station file, `storm_watch_station.conf`, `station_calibration/` and
setup's marker are git-ignored. New ground truth in `ground_truth/` is untracked. `git pull` touches
none of them.

**Tests on a new NUC.** `test_waterline_station_check.py` and `test_station_env.py` read CACO-05's
2025-11-13 calibration files, which a new NUC does not have. Run them on CACO-05 or in a development
checkout.

---

## 8. Known gaps at a new station

- **The DEM page draws no cameras.** `dem_figure.py` still draws the cameras from CACO-05's calibration
  file names.
- **Some outputs still say "UTM 19N".**
  - `survey_compare.py` labels its axes "UTM 19N"; pass it `--epsg`.
  - `foreshore_slope.py` writes a "UTM 19N" `.prj`.
- **The pointing log is off** (`POINTING_LOG_ENABLE="0"`). It needs a reference bank for the camera,
  and on a one-camera station it would still try c2.
- **The detector's tunings other than crop, envelope, bias, step and profile_scale** stay CACO05's
  (section 3).
- **CACO-05-only tools.** These read the cron's own lines, not the station file: `survey_products.py`,
  `historical_forcing.py`, `recalibrate_timing.py`, and `owg.sh` (not in the clone).

The tests behind this page: docs/VALIDATION.md section 2.11.

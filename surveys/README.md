# Survey dates: what survey_products.py builds, and against what

Two tables drive `survey_products.py`. Edit them, not the script, to add a
date, switch one on or off, or give a survey file that has arrived.

```
python3 survey_products.py --date 2025-01-23            # one date
python3 survey_products.py --date 2025-01-23 --dry-run  # what it would do: photos per day, files, labels, run time
python3 survey_products.py --all                        # every enabled date, then the summary
python3 survey_products.py --summary                    # one table + one figure across the dates built
```

Products go to `/mnt/I2Rgus_Data/survey_products/<date>/` (`--output-root`),
about 0.3 GB per two-camera week once the detector's debug images are removed
(they are, unless `--keep-detector-debug`: ~2 GB more).
Each step is skipped when its outputs exist, are newer than its inputs and
were made with the same settings, the same scripts (sha256) and the same
upstream builds; `--force` rebuilds, `--steps compare` (etc.) runs only some
steps. Every skip is printed. A rebuilt step makes everything after it
stale: a stale filtered file, DEM or comparison is never used or reported
as current (the README then says NO CURRENT COMPARISON, and why). The setup
coefficient and the DEM settings (cell, minimum frames, spread, Hs and
day-offset limits, GeoTIFF EPSG) reported are those the outputs on disk were
built with, never a later command line's: `--steps compare --dem-cell 1`
warns that the option is not applied and says to rerun `--steps
dem,maps,compare`.

Status and exit code: `complete`, exit 0; `partial`, exit 1: a step (a
comparison included) failed; `partial`, exit 3: a camera contributed nothing,
a part is missing (a camera's pointing check or its waterline map on the
photos: no photo of the window found, e.g. a wrong `--photo-roots`), or the
outputs on disk are not one build; a disabled date exits 2. With `--all` the worst code
is returned and the summary table has a status column.

## survey_dates.csv: one row per (date, camera)

| column | meaning |
|---|---|
| `date` | the survey date the product is for (its folder name) |
| `camera` | `c1` or `c2` |
| `enabled` | `1` = build; `0` = do not (the reason goes in `notes`). A date with no enabled row is refused, with the reason, unless `--force-disabled` |
| `station` | the station ID in the photo file names (`...GMT.2025.CACO03.c1.timex.jpg`). Only those photos are used: the ID changes when the cameras are set up again, so a window spanning a re-set cannot mix two pointings |
| `first_day`, `last_day` | the photo window (UTC days, inclusive) |
| `eo_file` | the calibration (pointing) for that window, in `calibration/` |
| `envelope_survey` | the survey the detection search envelope is placed with (`detect_original_view.py --envelope-survey`). **Never the survey the DEM is compared with** (that would make the comparison PARTLY-CIRCULAR): Jan 2025 photos use the Mar 2025 lidar and the other way round. Empty for the live era (the live detector has its own envelope, from hand-clicked ground truth) |
| `utc_hours` | daytime hours used, e.g. `13.5-18` (as process_chelsea.py); `0-24` = all (the live era, as the live DEM) |
| `era` | where the water level and waves come from: `adcp` (Signature 1000, Dec 2024 - Mar 2025), `chatham` (no ADCP: NOAA Chatham 8447435 transferred to Marconi + WIS/NDBC waves converted to the ADCP's Hs; same code path, `historical_forcing.py`), `live` (the station's own archive: GNSS-R levels and the cron's waves) |
| `notes` | why the row is as it is; printed when a date is refused |

File names (`eo_file`, `envelope_survey`, survey `path`) are looked for as
given, then by name in `--survey-dirs` (default: the Chelsea_calibration
folder, the waterline folder and its calibration folder, the repository's
`calibration/` and `surveys/`).

### The dates now

| date | status | photos | calibration | water level / waves | survey | label |
|---|---|---|---|---|---|---|
| 2024-10-23 | **skipped for now** (user's request) | c1+c2 CACO03, 20-26 Oct | CACO03_<cam>_20241023_EO (solved from that day's GCPs) | Chatham transfer / WIS | GCPs only, file not on the NUC yet | PARTLY-CIRCULAR |
| 2025-01-23 | enabled | c1+c2 CACO03, 18-23 Jan (24 Jan left out: re-set that day) | CACO03_<cam>_20250123_EO (GCPs, not the lidar) | ADCP | Jan 2025 lidar | INDEPENDENT (envelope from the Mar lidar) |
| 2025-03-06 | enabled | c2 only, CACO04, 3-9 Mar | CACO04_c2_20250219_EO | ADCP (ends 10 Mar 15:00) | Mar 2025 lidar | INDEPENDENT (envelope from the Jan lidar) |
| 2025-03-19 | **skipped** (user's request) | none on the NUC within 7 days | - | no ADCP | - | - |
| 2026-09-29 | enabled | c1+c2 CACO05, 26 Sep - 2 Oct (centred on the survey day; leaves out most of the 25-26 Sep storm), the live archive | CACO05_<cam>_20251113_EO-CV | GNSS-R / live waves | RTK check shots of 29 Sep | **CIRCULAR with C = 0.037** (see below) |

The live date was first called 27 Sep; the user confirmed that its survey is
the 29 Sep RTK, so the product is dated 2026-09-29.

**2025-01-23: an earlier check points the other way.** Before this script,
the same photos (18-23 Jan) with the same CACO03_<cam>_20250123 GCP
calibrations, unfiltered, with no setup and the envelope from the earlier
lidar-fitted pointing, gave waterline - Jan lidar **+0.46 m for c1** (NMAD
0.19) and **+0.15 m for c2** (NMAD 0.23): HIGH, not low. The setup (C = 0.037,
~0.28 m that week) raises every line further (by ~0.8 x the setup once the
landward re-projection is counted, as on the 2026 RTK), so unless the new
envelope changes the detections, expect about +0.7 m (c1) and +0.4 m (c2).
The biases listed under Labels that make the 2025 dates read LOW cannot
explain that; a 0.3 m difference between cameras on the same water level
points at the pointing or the lens model. The README of that date says so,
gives the headline per camera, the C = 0 sensitivity (the same lines without
the setup) and each camera's sea-horizon offset against its calibration.

## surveys.csv: one row per survey to compare with

| column | meaning |
|---|---|
| `date` | the product date it belongs to |
| `name` | short name, used in the output file names (`compare/<date>_<name>_*`) |
| `path` | the survey file (lidar DSM `.tif`/`.asc`, or points `.csv`: Emlid, CIRN `num,E,N,Z`, or E/N/Z headers). Empty = not available yet: the comparison is skipped with a message saying what to provide; for a points survey, files named `<survey_date>*xyz*.csv` or `<survey_date>*GCP*.csv` in the survey folders are then picked up by themselves (several, e.g. one per camera, are merged) |
| `survey_type` | `dsm` or `points` |
| `survey_date` | when it was surveyed (the time gap to the photos is reported) |
| `label` | how independent the comparison is, before the script's own checks (below) |
| `why` | the reason, printed next to the label everywhere |
| `label_override_reason` | optional. Empty = every downgrade applies. A downgrade from a rule that can misfire (a point survey on the calibration's day may be its GCPs; a points survey of the date a setup coefficient was fitted on may be those shots; survey_compare.py's name and date checks) is NOT applied when a reason is written here, e.g. `FIXTURE points sampled from the lidar, not the calibration's GCPs`; the reason and the downgrade it overruled are printed next to the label. A downgrade from what the chain actually did (the envelope, a fitted pointing, a setup coefficient fitted to the same file) cannot be overruled |

## Labels

* **INDEPENDENT**: nothing in the chain was fitted to or placed with this survey.
* **CROSS-VALIDATED**: fitted on other data from the same source, held out here.
* **PARTLY-CIRCULAR**: some step used this survey (e.g. the search envelope, or
  the calibration was solved from these very points).
* **CIRCULAR**: pointing or offsets fitted to this survey.

A comparison is never labelled better than its weakest step. The table's
label is a starting point; `survey_products.py` then checks the chain it
actually ran and can only make the label worse, saying why:

* a search envelope placed with this survey -> PARTLY-CIRCULAR;
* a camera pointing fitted to a survey (`fit_eo_to_survey.py`, `*_lidar_EO.yaml`) -> CIRCULAR;
  so is a pointing CARRYING such a fit's change to another period
  (`apply_pointing_correction.py`, `*_corr_EO.yaml`, e.g.
  CACO05_c1_20250219_corr_EO.yaml from the Jan lidar fit): the fit its notes
  name is read in turn (next to it, in `calibration/`, `--calibration` and the
  survey folders) and its survey matched by name or content; when it cannot
  be followed, CIRCULAR as a rule that can misfire (overrulable);
* a GCP survey on the calibration's own day -> PARTLY-CIRCULAR;
* a point survey made on the calibration's own day -> PARTLY-CIRCULAR, whatever
  the file is called (it may be the GCPs the calibration was solved from);
* **a setup coefficient fitted to this survey -> CIRCULAR.** The survey is
  recognised by its content (sha256), not its file name, so a renamed copy is
  still caught; so is a re-export (more than half its points within 0.05 m of
  the fit's shots, when that file is here); a points survey of the same date
  as the fit's (in surveys.csv or in the Emlid file's own 'Averaging start'
  column) is treated as possibly the same shots. C = 0.037 (`waterline_timex_cron.sh`, the station's
  C, used by default) is the median per-frame C of the waterlines lying on the
  2026-09-29 RTK transects, so the 2026-09-29 comparison with those shots
  cannot test the setup (nor the overall level it sets): it is CIRCULAR for as
  long as C = 0.037 is used. Its spread and its dependence on elevation still
  say something. **An independent 2026 check needs another survey** that
  fitted nothing, e.g. a calm low-tide RTK across the intertidal.
* the setup coefficient is read from the waterlines the DEM was built with
  (the detection stamps, checked against the setup the rows imply), never from
  the command line of a later `--steps compare`. A `--setup-coef` other than a
  known fit is refused unless `--setup-fitted-to` says where it came from: the
  survey file it was fitted to (comparisons with that survey become CIRCULAR)
  or `'none:<how>'`, e.g. `'none:repeat crossings, dem_from_contours.py
  --fit-setup on 2025-01-18..23'`. A C of unknown origin is never INDEPENDENT.

The setup itself is not exact either (the station's C is not changed here;
each README states the bias):
* **under-correction.** C = 0.037 was fitted with the RTK elevation where each
  line lay WITHOUT setup, but a line given the setup is re-projected landward
  onto higher beach, which closes only part of the gap. `compare_rtk.py` on
  the 18 frames of 29 Sep - 5 Oct lying on the transects gives RTK - waterline
  +0.33 m without setup and still +0.13 m with C = 0.037 (median setup applied
  0.26 m; Oct 2026): the 2026 lines still read ~0.13 m low on the RTK, about
  half the applied setup. Every product made with C = 0.037 may read low by
  roughly 0.5 x its setup (each README gives the figure for its date). The
  synthetic test beach shows the same geometry: a line mapped onto the
  still-water plane reads low by about 0.65-0.7 x the setup.
* **still-water reference (2025 dates).** C was fitted with each line at the
  GNSS-R water level. The GNSS-R footprint is the surf zone (`gnssr_qc.py`,
  `gnssir_reflection_audit.py`), where breaking waves raise the mean level:
  GNSS-R already contains part of the setup, and C carries only the rest. The
  2025 dates use the ADCP at 21 m depth (or the Chatham harbour transfer),
  which see no setup: their lines are expected to read LOW by the share of
  the setup the GNSS-R sees, and ~0.02 m HIGH from the mean levels (GNSS-R sits
  ~0.02 m below Chatham; the ADCP datum assumes Marconi = Chatham). On the
  station, `historical_forcing.py` measures that share from the 2026 record
  (GNSS-R spline, `archive/gauge_8447435.csv`, `archive/waves_marconi.csv`)
  and the README gives it in metres. To take the question out, fit C in the
  date's own frame: build with `--setup-coef 0`, run `dem_from_contours.py
  --fit-setup` on `waterlines/contour_points_ground.csv` (repeat crossings, no
  survey), then rebuild with `--setup-coef <C> --setup-fitted-to 'none:repeat
  crossings, <window>'`; the lidar comparison stays INDEPENDENT.
* **wave currency.** C was fitted with ADCP-currency Hs and NDBC 44008 peak
  periods; the 2025 dates use the ADCP's periods (or WIS converted to them):
  `historical_forcing.py` states the expected setup bias in metres.
* frames with no wave record (no setup) are left out of the waterlines when
  C > 0 (listed in `waterlines/no_setup_frames.csv` and counted in the README):
  kept, they would read about one setup low.

Coordinates: the grids are in the frame of the calibration's GCPs and of the
surveys, NAD83(2011) / UTM 19N (EPSG:6348). The GeoTIFFs are tagged
EPSG:32619 (WGS 84 / UTM 19N) as asked: a nominal tag (WGS 84 differs by
~1-1.5 m here; nothing is transformed). `--geotiff-epsg 6348` writes the true
code. The 2025 lidar files state no vertical datum: NAVD88 is assumed (the
geoid model is unknown; GEOID12B and GEOID18 differ by a few cm here).

`survey_compare.py`'s own checks (calibration notes, the envelope, the
`provenance.json` next to the DEM) run before the comparison, and the
headline label is the worst of all of them, with every reason. On RTK
transects the waterlines are compared as `compare_rtk.py` does (one value
per frame, the RTK interpolated along the transect); isolated points
(GCPs) are matched to the nearest line point, not slope-corrected. The
headline is given per camera too, and, when the lines carry a setup, again
WITHOUT it (C = 0: each line at its still-water level, re-projected with the
calibration as a C = 0 build would place it) as a sensitivity.

Pointing: each camera's sea horizon over the window is compared with where
its calibration puts it; the constant tilt/roll is printed for every camera
and, above 0.05 deg, given with its DEM sensitivity per range band (0.1 deg
of tilt is ~0.13-0.17 m of DEM at 250-350 m). Only day-to-day CHANGES of
pointing (`horizon_check.py`, 12 px) leave days out.

## Run time on the NUC

Detection dominates: `detect_original_view.py` runs the detector on the
whole 2448 x 2048 frame, measured at 22-26 s per photo here (not yet on the
NUC10i3; its progress lines print the real rate), i.e. **~50 min for a
two-camera week** of daytime photos. Everything else takes a few minutes. The
run lowers its own priority (`--nice 10`, the default) and warns when it
would overlap the station's jobs (waterline cron 12:30, 19:25, 20:55;
cleanup.sh 19:45, 21:00; owg.sh hourly at :40); a start after ~21:30 local
keeps clear of them. `--dry-run` prints the estimate.

## Deploying on the station

The work is on a branch until it is merged; the NUC runs `main`. After the
merge, on the NUC:

```
cd /mnt/I2Rgus_Data/waterline && git pull
python3 survey_products.py --date 2025-01-23 --dry-run   # inputs found, labels, run time
python3 survey_products.py --date 2025-01-23             # then 2025-03-06, 2026-09-29, --summary
```

## When the October GCP file arrives

Save it as `/mnt/I2Rgus_Data/Chelsea_calibration/2024-10-23_Marconi_Extrinsic_Targets_<cam>_xyz.csv`
(no header: `num,E,N,Z`, UTM 19N / NAVD88, like
`calibration/2025-11-13_Marconi_Extrinsic_Targets_c1_xyz.csv`), or write its
path in the `oct_gcps` row of `surveys.csv`. When the October date is
switched on (`enabled` = 1 in `survey_dates.csv`, or `--force-disabled`), the
comparison runs. It is PARTLY-CIRCULAR: the 2024-10-23 calibration was solved
from those points (the pointing is fitted to them; the elevations still come
from the water levels). Days before 23 Oct may have had another pointing
(the 2024-08-27 calibration differs by ~2.3-2.5 deg in azimuth): the pointing
step checks every day against a photo of 23 Oct and leaves out the days that
differ.

The October date also needs its forcing files, which are NOT downloaded:
* WIS ST63064 hourly waves (`WIS_ST63064_2024-10_to_2025-03.csv`) and/or
  NDBC 44013 (`NDBC_44013_2024-10_to_2025-03.csv`) in
  `/mnt/I2Rgus_Data/Chelsea_calibration` (or `--wis` / `--ndbc`): with the
  setup on, `historical_forcing.py` fails (exit 2) when waves cover less than
  half of the window's daytime hours, before any detection;
* the Chatham 8447435 record (`chatham_2024-10_to_2025-03.csv`; downloaded
  from NOAA when missing, unless `--no-download`) and the ADCP files (the
  transfer is fitted on their overlap).

Run `--date 2024-10-23 --force-disabled --dry-run` first: it lists each
forcing file found or NOT FOUND.

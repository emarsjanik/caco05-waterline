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
dem,maps,compare`. The DEM's day test (a camera-day off the DEM made
without it by more than 0.15 m is left out) is not run on fewer than 3
camera-days: with fewer than 3 the test cannot tell which day is off
(`dem_from_contours.py` keeps them all in that case, and exits 4 rather
than crash if a test ever rejects them all), so the option is not passed;
every camera-day is kept and any difference between them stays in the DEM's
spread. The README says so, and a DEM that still fails says why (from
`logs/dem.log`). The photo window and hours and the input
files are likewise those of the build: a `--window` or `--utc-hours` is applied only by a run of every
step from the forcing through the DEM (otherwise, e.g. `--steps compare
--window ...`, the run keeps the window the outputs were built for and the
README says `NOT APPLIED`; a full run whose first step fails, e.g. a live
`--window` with no archive rows, keeps the rows, stamps and window of the
earlier build and says `NOT APPLIED` too), and an input file (`--live-contours`, the
forcing files, `--photo-roots`) whose step failed or did not run is reported
the same way; the "rebuild everything" line always names the window, hours
and files the outputs on disk came from (their build stamps).

Status and exit code: `complete`, exit 0; `partial`, exit 1: a step (a
comparison included) failed; `partial`, exit 3: a camera contributed nothing,
a part is missing (a camera's pointing check or its waterline map on the
photos: no photo of the window found, e.g. a wrong `--photo-roots`; or a
survey whose file `surveys.csv` names was not found, e.g. a misnamed lidar or a
wrong `--survey-dirs`, or has no current comparison), or the outputs on disk
are not one build; a disabled date exits 2; exit 4: not started, because
another build of the same date is running; exit 5: not started, refused by
the run's own settings (a `--setup-coef` other than the one the date's
waterlines carry, without `--replace`; a C of unknown origin, without
`--setup-fitted-to`; `--dry-run` says when a real run would be refused); exit
6: not started, a usage error (an unknown option or a bad value, e.g.
`--steps`, `--window`, `--utc-hours`, a negative `--setup-coef`), never 2. A
row with no survey file yet (the October GCPs) is not a missing part. With
`--all` the code of the worst date is returned, worst by severity, not by
number: 1 if any date failed, else 5 if any was refused, else 4 if any was
busy, else 3 if any was partial, else 0 (disabled dates do not count); the
summary table has a status column. A missing waterline map
on the photos says why: no photo of the window found (the photo roots), or
`daily_elevation_map.py` failed (its exit code and log).

One build of a date at a time: a build holds `<date>/.build.lock` (created
atomically, with its pid, computer and command) while it runs. A second build
of the same date, e.g. started from another SSH session or by `--all` while a
`--date` run is going, refuses to start and says which build holds it (the
two would share every work file and mark each other failed). A lock left by
a build that is no longer running (killed, power cut) is taken over, and said
so.

A build first marks the product `in progress` (provenance.json, and a
`BUILD IN PROGRESS / INTERRUPTED` banner over the previous README) and
rewrites both when it ends. A build that dies on the way (Ctrl-C, a full
disk) writes `status: failed` with the reason, or, if even that cannot be
written, leaves the mark: the previous build's numbers are never presented
as current over files that have since changed. `--summary` re-checks every
comparison against the files on disk (the build stamps) rather than trusting
provenance.json. The detection checks first that the disk has room for it
(~12 MB a photo of debug images until each camera is done, + ~0.3 GB).

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

A text with a comma in it (`notes` here; `why` and `label_override_reason` in
`surveys.csv`) must be inside double quotes, as the shipped rows are:
`...,live,"off for now, waiting for the October recheck"` (a double quote
inside the text is written twice: `""`). Without the quotes the comma ends
the text: `survey_products.py` warns (`... field(s) more than its ... columns
(a comma in a text that is not in double quotes?)`) and drops what follows,
and the last group of `test_survey_products.py` (the station's own tables)
FAILs on that warning. A row that only ends in a comma (an empty field past
the last column, as a spreadsheet may save it) loses nothing: no warning.
In `surveys.csv`, a row that leaves off its optional
last column gives no warning at all: the text after the comma becomes that
row's `label_override_reason`, which overrules the downgrades that column can
overrule (`--dry-run` then prints `[overruled in surveys.csv (<the text>):
<the downgrade>]` next to the label, a build warns `downgrade NOT applied,
overruled in surveys.csv`, and `test_survey_products.py` lists every
`label_override_reason` the tables hold).

File names (`envelope_survey`, survey `path`) are looked for by name in
`--survey-dirs` (default: the Chelsea_calibration folder, the waterline
folder and its calibration folder, the repository's `calibration/` and
`surveys/`), and `eo_file` in `calibration/` (`--calibration`): never in the
folder the command is run from, so a build finds the same files wherever it
is started, and `--dry-run` and `--summary` label a survey from the same
files the build would use. A path with a folder in it (absolute, or e.g. `./x.csv`) is
taken as given when it exists, else by its name in those folders.

### The dates now

| date | status | photos | calibration | water level / waves | survey | label |
|---|---|---|---|---|---|---|
| 2024-10-23 | **skipped for now** (user's request) | c1+c2 CACO03, 20-26 Oct | CACO03_<cam>_20241023_EO (solved from that day's GCPs) | Chatham transfer / WIS | GCPs only, file not on the NUC yet | PARTLY-CIRCULAR |
| 2025-01-23 | enabled | c1+c2 CACO03, 18-23 Jan (24 Jan left out: re-set that day) | CACO03_<cam>_20250123_EO (GCPs, not the lidar) | ADCP | Jan 2025 lidar | INDEPENDENT (envelope from the Mar lidar) |
| 2025-03-06 | enabled | c2 only, CACO04, 3-9 Mar | CACO04_c2_20250219_EO | ADCP (ends 10 Mar 15:00) | Mar 2025 lidar | INDEPENDENT (envelope from the Jan lidar) |
| 2025-03-19 | **skipped** (user's request) | none on the NUC within 7 days | - | no ADCP | - | - |
| 2026-09-29 | enabled | c1+c2 CACO05, 26 Sep - 2 Oct (centred on the survey day; leaves out most of the 25-26 Sep storm, whose last day is in: the README's STORM caveat), the live archive | CACO05_<cam>_20251113_EO-CV | GNSS-R / live waves | RTK check shots of 29 Sep | **CIRCULAR with C = 0.037** (the default); INDEPENDENT with C = 0 (see below) |

The live date was first called 27 Sep; the user confirmed that its survey is
the 29 Sep RTK, so the product is dated 2026-09-29.

**2025-01-23: the earlier check of the real photos** (itself
**PARTLY-CIRCULAR**: its search envelope was placed with this same Jan lidar,
through the earlier lidar-fitted pointing). Before this script,
`pointing_fix_figure.py` on the REAL photos of 18-23 Jan, with the same
CACO03_<cam>_20250123 GCP calibrations, unfiltered and with no setup, printed
*lidar elevation - water level at the waterline* of **+0.46 m for c1** (NMAD
0.19) and **+0.15 m for c2** (NMAD 0.23). Positive there means the lidar
beach lies ABOVE the level each line was given: the lines read **LOW**, i.e.
waterline - Jan lidar **c1 -0.46 m, c2 -0.15 m** without setup. That is the
sign a line given no setup must have (a timex line is marked where the swash
reaches), the sign of the biases listed under Labels and the sign the C = 0
sensitivity of this product should show: no discrepancy of sign. The setup
(C = 0.037, ~0.28 m that week) raises every line; once re-projected landward
each line keeps a share of its setup that depends on the date, the camera and
the beach (0.71 per frame on the 2026 RTK; ~0.6 and ~0.8 in the Jan and Mar
2025 builds), so unless the new envelope changes the detections, expect
waterline - Jan lidar of about **c1 -0.29 to -0.24 m, c2 +0.02 to +0.07 m**
over shares of 0.6-0.8 (c1 -0.18, c2 +0.13 m if the whole setup counted). The
build gives the expectation from its own paired C = 0 sensitivity, per
camera. It already holds the re-projection and any mismatch of C on this
date (below): they are not to be added again. The 0.3 m difference between
the cameras, measured over the same hours, points at the pointing or the lens
model. The README of that date says so, gives the
headline per camera, the C = 0 sensitivity (the same lines without the setup)
and each camera's sea-horizon offset against its calibration.

The real January (and March) photos are on the station computer only. The
numbers in a build made here, from the synthetic test fixture (photos
rendered from the lidar with the ADCP water level and the setup planted), test
the pipeline, not the real beach: never read them as real-data results. The
real numbers for these dates come from running `survey_products.py` on the
station.

## surveys.csv: one row per survey to compare with

| column | meaning |
|---|---|
| `date` | the product date it belongs to |
| `name` | short name, used in the output file names (`compare/<date>_<name>_*`) |
| `path` | the survey file (lidar DSM `.tif`/`.asc`, or points `.csv`: Emlid, CIRN `num,E,N,Z`, or E/N/Z headers). Empty = not available yet: the comparison is skipped with a message saying what to provide; for a points survey, files named `<survey_date>*xyz*.csv` or `<survey_date>*GCP*.csv` in the survey folders are then picked up by themselves (several, e.g. one per camera, are merged) |
| `survey_type` | `dsm` or `points` |
| `survey_date` | when it was surveyed (the time gap to the photos is reported) |
| `label` | how independent the comparison is, before the script's own checks (below) |
| `why` | the reason, printed next to the label everywhere (with a comma in it: inside double quotes, see above) |
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
  cannot test the setup (nor the overall level it sets): this rule makes it
  CIRCULAR for as long as C = 0.037 is used (the table row says INDEPENDENT,
  which a build with C = 0, or with a C fitted to other data, keeps: nothing
  else in its chain used the shots; build it in its own root, see "a build
  with C = 0" below). Its spread still says something, over
  the narrow band the transects cover (they start at ~+1.2 m: the waterline
  check covers ~+1.0 to +1.6 m only). **An independent check of the whole
  2026 intertidal needs another survey** that fitted nothing, e.g. a calm
  low-tide RTK across it.
* the setup coefficient is read from the waterlines the DEM was built with
  (the detection stamps, checked against the setup the rows imply), never from
  the command line of a later `--steps compare`. A `--setup-coef` other than a
  known fit is refused (exit 5) unless `--setup-fitted-to` says where it
  came from: the survey file it was fitted to (comparisons with that survey
  become CIRCULAR) or `'none:<how>'`, e.g. `'none:repeat crossings,
  dem_from_contours.py --fit-setup on 2025-01-18..23'`. A C of unknown origin
  is never INDEPENDENT.

The setup itself is not exact either (the station's C is not changed here;
each README states the bias):
* **what the 2026 RTK can say about C.** The Oct 2026 recheck (the
  station's rows of 29 Sep - 5 Oct, 18 frames lying on the RTK transects both
  ways, the filter of that time) found waterline - RTK -0.33 m without setup
  and -0.13 m with C = 0.037 (median setup applied 0.26 m). Both medians are
  **biased LOW by where that survey stops**: the transects end at their
  lowest shot (+1.17 to +1.80 m NAVD88); a line below it can only read low
  there (the RTK under it is higher), and a line that would read high at that
  water level lands seaward of the survey and is never compared. Over only
  the lines at least 0.1-0.2 m above the lowest shots, that recheck (a
  scratch analysis, not code in this repository; its numbers are held in
  `SETUP_FITS`) read the lines with C = 0.037 at -0.04 to -0.03 m (7-10
  frames: too few to be an estimate), and a fit with each line re-projected
  at still water + setup gave C ~0.043-0.046 (day-block 90% ~0.02-0.055);
  the within-cell slope of repeat crossings (no survey, the same analysis)
  gave 0.043 (0.026-0.053). 0.037 lies inside both 90% intervals, so the station
  keeps it; the RTK, which stops above the low-tide beach, cannot pin C more
  closely (a survey reaching the low-tide line on days with different waves
  would). Each build measures the same thing on its own frames:
  `survey_compare.py` reports, next to the transect headline, the check over
  only the lines at least 0.10 m above the transects' lowest shots (WHERE
  THE SURVEY STOPS; e.g. -0.069 m on 7 frames in a test build of the 26 Sep -
  2 Oct window from a copy of the station's rows that starts on 29 Sep: on
  the station, with 26-28 Sep in, expect another value). Its README gives
  both and says which is this build's; where they differ, read the build's.
  A line given the setup is re-projected
  landward onto higher beach, so it keeps only part of its setup: ~0.71 x per
  frame on those RTK frames, ~0.6 (c1 0.64, c2 0.59) and ~0.8 in the paired
  C = 0 sensitivities of the Jan and Mar 2025 builds (real lidar beaches and
  calibrations, synthetic photos). The share kept is not one number to carry
  to other dates: it depends on the slope, the camera and its distance; each
  README gives its own build's per camera, and in metres. For 2025-01-23 the
  earlier-check expectation above already holds it.
* **still-water reference (2025 dates).** C was fitted with each line at the
  GNSS-R water level. The GNSS-R footprint is the surf zone (`gnssr_qc.py`,
  `gnssir_reflection_audit.py`), where breaking waves raise the mean level:
  GNSS-R contains a share s of the setup, and C carries only the rest. The
  2025 dates use the ADCP at 21 m depth (or the Chatham harbour transfer),
  which see no setup. The MEAN part of the GNSS-R's share is already in the
  datum chain: GNSS-R sits ~0.02 m below Chatham after its +0.349 m datum fix
  (a comparison of mean levels, which holds that mean share), and the ADCP
  datum is tied to Chatham's mean. So a 2025 frame is expected at about
  **+0.02 m - s x (its setup - the mean setup of the GNSS-R/Chatham comparison
  period)**: ~0.02 m HIGH in average waves, LOW only by the share of the setup
  above that average. On the station, `historical_forcing.py` measures s from
  the 2026 record (GNSS-R spline, `archive/gauge_8447435.csv`,
  `archive/waves_marconi.csv`; a wave-dependent effect) and the README gives
  the window's expected effect in metres. For 2025-01-23 the earlier-check
  expectation already holds this effect (that check was made on lines at the
  same ADCP still water); the caveat explains why the lines may differ from
  the lidar, not from that expectation. To take the question out, fit C in the
  date's own frame, each build in its OWN `--output-root` (in the default
  root a build with another C would replace the C = 0.037 product of that
  date, so it is refused there without `--replace`; for the 2025 dates its
  detection runs again, ~50 min a two-camera week):
  ```
  python3 survey_products.py --date 2025-01-23 --setup-coef 0 \
      --output-root /mnt/I2Rgus_Data/survey_products_c0
  python3 dem_from_contours.py \
      /mnt/I2Rgus_Data/survey_products_c0/2025-01-23/waterlines/contour_points_ground_filtered.csv \
      /tmp/fit_2025-01-23 --fit-setup --max-hs 1.5 --no-plot
  ```
  and read what it prints: the median spread of repeat crossings against C,
  the near-best band (the C whose spread lies within 5 mm of its best; 0.02
  or wider is a FLAT minimum), and how C = 0 and the C in use (the station's
  0.037, from `waterline_timex_cron.sh`) compare. It offers a C to apply
  only when C = 0 and the C in use both lie outside that band and the
  minimum is not flat; otherwise it says which of them it cannot tell from
  its best. One week's curve is noisy: on the C = 0 lines of 29 Sep - 5 Oct
  2026 this recipe printed a band of +0.006 to +0.034 (best 0.026, FLAT;
  0.037 spreads 6 mm more than the best), the unfiltered rows, also with
  `--max-hs 1.5`, +0.026 to +0.040 ("not distinguishable from the C in
  use"), and 29 Sep - 2 Oct alone -0.020 to +0.010 ("not distinguishable
  from C = 0"); the same spread with each line re-projected gave 0.029 (90%
  ~0.00-0.04). Adopt a C only if it
  holds over weeks of different waves, with `--setup-coef <C>
  --setup-fitted-to 'none:repeat crossings, <window>' --output-root
  /mnt/I2Rgus_Data/survey_products_c<C>`; the lidar comparison stays
  INDEPENDENT. Otherwise keep the station's C.
* **wave currency.** C was fitted with ADCP-currency Hs and NDBC 44008 peak
  periods; the 2025 dates use the ADCP's periods (or WIS converted to them):
  `historical_forcing.py` states the expected setup bias in metres.
* frames with no wave record (no setup) are left out of the waterlines when
  C > 0 (listed in `waterlines/no_setup_frames.csv` and counted in the README):
  kept, they would read about one setup low.
* **a build with C = 0** (`--setup-coef 0`) writes `setup_correction_m` 0.0 on
  every row: its comparison says "NO setup correction (C = 0)" and has no C = 0
  sensitivity (both sides would be the same lines). This is the build that
  makes the 2026-09-29 RTK comparison INDEPENDENT. Build it in its own root,
  so the station-C product stays:
  ```
  python3 survey_products.py --date 2026-09-29 --setup-coef 0 \
      --output-root /mnt/I2Rgus_Data/survey_products_c0
  python3 survey_products.py --summary --output-root /mnt/I2Rgus_Data/survey_products_c0
  python3 survey_products.py --summary          # the C = 0.037 builds, default root
  ```
  Run in the default root it is REFUSED (exit 5: it would replace the
  C = 0.037 product of that date, and `--summary` would then show only the
  C = 0 build); the refusal names the root that C belongs in (the default
  root for the station's C, `<default root>_c<C>` for any other).
  `--replace` overwrites on purpose. For a 2025 date the build re-runs the
  detection (~50 min a two-camera week: each line is placed at the
  elevation its setup gives it). Its numbers are biased too: every C = 0
  line on the 2026 transects lies below the survey's lowest shot (WHERE THE
  SURVEY STOPS).
* **GNSS-R datum (live dates).** The live water level rests on one survey of
  the antenna (NGS OPUS, GEOID18: +0.349 m, +/-0.061 m, mostly the geoid
  model). It largely cancels against GEOID18 RTK shots, but every absolute
  NAVD88 height of a live DEM shares it.

Coordinates: the grids are in the frame of the calibration's GCPs, ASSUMED
to be that of the surveys, NAD83(2011) / UTM 19N (EPSG:6348): the GCP files
state no CRS. Were the GCPs in WGS 84 / ITRF instead, the grid would sit
~1-1.5 m off the surveys, up to ~0.1-0.15 m of elevation where that shift runs
across a 1:10 foreshore. The GeoTIFFs are tagged EPSG:32619 (WGS 84 / UTM 19N)
as asked: a nominal tag (WGS 84 differs by ~1-1.5 m here; nothing is
transformed). `--geotiff-epsg 6348` writes the true
code. The 2025 lidar files state no vertical datum: NAVD88 is assumed (the
geoid model is unknown; GEOID12B and GEOID18 differ by a few cm here).

`survey_compare.py`'s own checks (calibration notes, the envelope, the
`provenance.json` next to the DEM) run before the comparison, and the
headline label is the worst of all of them, with every reason. On RTK
transects the waterlines are compared as `compare_rtk.py` does (one value
per frame, the RTK interpolated along the transect); isolated points
(GCPs) are matched to the nearest line point, not slope-corrected, and the
statistic is judged on the DISTINCT survey points (many frames see the same
few points: fewer than 10 points is TOO FEW whatever the number of values).
The headline is given per camera, per 0.5 m band of survey elevation, with
WHERE on the beach the compared line points lie (next to all the lines'
range), and, when the lines carry a setup, again WITHOUT it (C = 0: each line
at its still-water level, re-projected with the calibration as a C = 0 build
would place it) as a sensitivity, read on the SAME frames (points) both ways:
lines at still water move seaward and some leave the survey; those are
counted. The line statistics are those of the frames the DEM gridded (its
own wave and day filters, read from `dem/<date>_info.json`); the frames it
left out (offshore Hs > 1.5 m: the largest setups, where C and the still
water err most) are compared apart and counted, with the statistic of all
frames next to them. `--summary` gives each camera's median too (a pooled
median near zero can hide two cameras off in opposite directions: more than
0.1 m apart is flagged when each camera rests on at least 10 frames, cells
or survey points, else 'too few to compare the cameras'; the difference is
put on the cameras, pointing or lens model, only when both were compared at
the same times, else each camera's median setup and its values with and
without the setup are given), with what each rests on (frames or distinct survey
points, 'too few' below 10), the C = 0 sensitivity (marked the same way), the
setup coefficient of each build and, for a label the table gave, why from the
build (e.g. "C = 0: no setup, nothing fitted to this survey"). Dates enabled
in the configuration with no product yet are listed as such.

Beach change: besides the storms measured on the station (the 25-26 Sep 2026
storm), each README looks at the date's own forcing: a span of the window with
offshore Hs >= 2.0 m gets a HIGH WAVES caveat counting the frames the DEM
gridded before, between, during and after it (by capture time; from the DEM's
own record, `dem/<date>_info.json`: its window, wave filter and rejected
camera-days), the frames it left out and why, and where the survey lies
(2025-01-23: Hs 2.9 m on 20 Jan, between the early frames and the lidar
flight); spans less than 48 h apart are one caveat. The STORM caveat says the
same of the storm's tail before the survey: which days the DEM grids, and for
each day it does not, why (left out by the pointing check, every frame taken
by its wave filter or day test, or no waterline). The note on the per-day
rows below is given once for all of them. The per-day rows of
waterlines - survey show whether the lines moved, but a step or trend there
may be change of the beach OR a method error that depends on the conditions:
on the static synthetic fixture (every day rendered from one lidar) the
per-day medians drift by ~0.09 m over Jan 18-23 and ~0.17 m over Mar 3-9,
with a ~0.07 m step across the 20 Jan event, and two Mar days at the same
median elevation differ by ~0.13 m: the error depends on the elevation, the
waves and setup and the tide phase of the sampled hours, which moves from
day to day. A difference between days below ~0.15-0.2 m is not by itself
evidence of change. Each row gives the median survey elevation its values
were taken at ('at z') and the median setup of its lines: compare days at a
similar elevation and setup, or within one elevation band, before calling a
difference change.

Pointing: each camera's sea horizon over the window is compared with where
its calibration puts it; the constant tilt/roll is printed for every camera
and, above 0.05 deg, given with its DEM sensitivity per range band (0.1 deg
of tilt is ~0.13-0.17 m of DEM at 250-350 m). Only day-to-day CHANGES of
pointing leave days out, decided per day from two measurements
(`survey_products.py` pointing_verdicts; the thresholds and where they were
measured: docs/VALIDATION.md 4.6):

* features matched against the reference frame on 40 or more points
  (pointing_check.py) are good to ~0.1 deg, and a turned match is the scene
  (dirt matches at zero turn): a change of tilt or roll over 0.3 deg (the
  step's tilt/roll tolerance; 0.1 deg of tilt is ~0.13-0.17 m of DEM) or of
  pan over 0.5 deg is a DIFFERENT POINTING, whatever the horizon reads; a
  smaller turn over 0.1 deg is a small change (kept), unless a CLEAR horizon
  reads a change over 0.3 deg that the features confirm within their own
  ~0.1 deg (then DIFFERENT POINTING, whatever the number of points); a turn
  of 0.1 deg or
  less keeps the day against a horizon reading that is not clear (the why
  gives both numbers), but against a CLEAR horizon edge that jumps the day
  is UNCONFIRMED: dirt on the housing window matches at zero turn when fog
  hides the scene. While the scene still matches, a second match consensus
  of 40+ points turned elsewhere shows the dirt ('ambiguous'); when it does
  not, a zero turn no further from none than the match's own noise
  (`turn_sigma` under 4: dirt moves with the camera and matches at EXACTLY
  zero, a scene matched between two days shows the camera's own small
  turns) is no reading at all, except on the reference day;
* a sea-horizon reading counts only on a CLEAR edge: one sharp edge (12x the
  median gradient of its +/-250 px search or more, no rival over 0.65 of it),
  bright sky over darker sea on 2/3 of the columns or more, and more than
  50 px from the end of the search. In a storm the strongest edge is the top
  of a foam or fog band (26-27 Sep 2026: 187-233 px below the horizon, fitted
  to 1-7.5 px), and the fit residual cannot tell. A misread never enters the
  window's median, offset or clear-day count;
* every day is judged against the REFERENCE's tilt/roll as the horizon reads
  it (the clear days' horizon less the turn their matched features
  measured; the reference day's own group when another day confirms it), as
  the features are, not against the window's median, which falls between two
  pointings when a camera moved on half the days; a clear reading that its
  own matched features contradict is no clear day;
* without matched features, a jump of tilt or roll over 0.3 deg from the
  reference is a DIFFERENT POINTING only on a clear edge two or more days
  share; otherwise the day is UNCONFIRMED. A day nothing could read (no
  matched features, and no horizon fitted or one that is not clear) next to
  a day measured to differ, or UNCONFIRMED on more than a horizon reading
  that is not clear (a clear-looking jump, an ambiguous match), is
  UNCONFIRMED too, and so are runs of such days from it: nothing shows the
  camera back. So is such a day in a storm: a run of days, each left out on
  its own readings or read by nothing, that holds misread horizons on two
  days or more, or a day with offshore Hs over 1.5 m in the date's own
  forcing (the DEM's limit) or inside a known storm (a camera can turn then:
  c1, 26-28 Sep 2026). On calm days a single misread day never reaches past
  itself. A day whose clear horizon shows tilt and roll back but that has no
  feature reading ('horizon only') is UNCONFIRMED next to a day that measured
  a pan over 0.5 deg: no horizon sees a pan. Nothing reaches the reference
  day (the others are measured against it); the CSV's `next_to` names the
  day a run starts at.

DIFFERENT POINTING and UNCONFIRMED days are both left out of the detection
(`--keep-moved-days` keeps them), whatever share of the window they are:
when the days left out are most of the window, the step and the README say
so with the count; when the days measured to differ agree with each other
and outnumber those at the reference's pointing, they name that other
pointing (the window holds two, and the reference may be the odd one out:
look at the photos), give the `--window` that builds on the other pointing
instead (in another output root), and the build's status is partial. Until Oct 2026
more than half of the measured days differing kept EVERY day, which put two
pointings and the unconfirmed storm days in one product. There is no option that leaves out or keeps single
days (`--window FIRST LAST` narrows the photo days). The README takes which
days were left out or kept, and why, from the detection's own record (a
later `--steps compare` run without the option does not change it; a record
of an older build that does not say whether it ran with `--keep-moved-days`
is never read as having done so); a build whose detection has not run, or
ran with other pointing results (the pointing step run again on its own),
says so, and calls days kept only once the detection has run with these
results. Window dirt is caught only by these two tests (a second consensus, a
turn within the match's noise): a dirt consensus of thousands of points may
read above the noise test, and dirt under a clear horizon that does not jump
can hide a pan (docs/VALIDATION.md 4.6). The printed table, the CSV
(`pointing/pointing_<cam>.csv`: the edge numbers, the inliers,
`turn_sigma`, the runner-up) and the README give each day's verdict and why. The predicted horizon
(`estimate_eo_rotation.horizon_rows`) holds the full Earth-curvature dip with
standard refraction (k = 0.13): a perfectly calibrated camera reads 0.00 deg
(+/-0.01 for k from 0 to 0.25). Until Oct 2026 it held half the dip, which
read +0.06-0.07 deg on a perfect camera (so the earlier "c2 ~0.4 deg off" is
~0.33 deg). On the synthetic test fixture the horizon is an artefact of the
rendering (not a curved sea): the offsets a fixture build reports are not
camera offsets. Each README says where the calibration puts the horizon: on
the live CACO05 20251113 calibrations it falls at the top edge of the frame
(c1: rows ~140 to 1, then off the frame), where the lens model is
extrapolated far beyond its GCPs, and the real 2026 photos read 34-41 px
(~0.6-0.9 deg) off. Lens model or pointing? The 2026 RTK cannot decide it
with C = 0.037, because that C was fitted to the same RTK shots: any level
error on the transects, of the setup or of the pointing, is absorbed into
C. The README therefore judges it on what fitted nothing to the RTK: the
lines WITHOUT setup (C = 0) against the RTK, minus the setup that a C fitted
to repeat crossings explains (no survey; a constant pointing error moves
every frame's line alike and does not show there; Stockdon's 0.35 x a ~1:10
beach face expects about the same). That C is taken as 0.029-0.043, the two
re-projected estimates of the Oct 2026 recheck (a scratch analysis, not code
in this repository): 0.029 for the spread `dem_from_contours.py
--fit-setup` minimises, a flat minimum (90% ~0.00-0.04; `--fit-setup`
itself, which does not re-project, prints a near-best band of +0.006 to
+0.034 there on the C = 0 lines), and 0.043
(0.026-0.053) for the within-cell slope of elevation on sqrt(Hs L0). What is left
is compared with what the full offset would do at the transects' range from
each camera (~55 m from c1, ~110-120 m from c2, where the full offsets would
put the lines ~0.4 and ~0.6 m low: not the metre-level shifts the table gives
at 200-350 m, where no survey exists). On the 2026 lines up to about half of
c1's offset and a quarter of c2's cannot be ruled out (5-7 frames per camera,
too few for an estimate; swash and detector bias at the waterline are
confounded with it): the full offset is unlikely to be pointing, a part of
it cannot be ruled out. Even that part is OVERSTATED: every one of those
C = 0 lines lies below its transect's lowest shot (the survey stops at +1.17
to +1.80 m), where a line can only read low, so what is left after the setup
is biased low by where the survey stops (WHERE THE SURVEY STOPS, above). The
README's DEM-shift table says what the offset WOULD do if it were all
pointing.

Synthetic inputs: a build from the sandbox test fixture (photos rendered from
the lidar with the water level and setup planted) is run with `--synthetic
<note>` (or a file named `SYNTHETIC_FIXTURE` in a photo root): its README,
comparison pages and summary rows then say "SYNTHETIC FIXTURE: tests the
pipeline, not the beach". The live photo maps (`waterlines/<cam>/src`) are
linked again by each run's maps step from the photo roots of that run, so a
photo that turns up later (a corrected `--photo-roots`, a mount back) is
drawn without a new georectification.

## Run time on the NUC

Detection dominates: `detect_original_view.py` runs the detector on the
whole 2448 x 2048 frame, measured at 17-26 s per photo off the station (a
2.8 GHz Xeon), i.e. **~50-70 min for a two-camera week** of 13.5-18 UTC
photos (120-165 photos: the count depends on the date's photo cadence; the
dry run prints it). The rate on the station's NUC10i3 has NOT been
measured: no date with detection has been built there yet (the 2026-09-29
product runs none), so budget 1-2 h and read the real rate from the first
`detecting : photo k of N, ... (S s per photo)` line of the build log.
Everything else takes a few minutes.
The live date (2026-09-29) runs no detection: it reads the station's
contour file (~1 min per GB) and takes a few minutes in all (79 s here for
67 MB of rows). The
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

One setup coefficient per output root: a run that would rebuild a date with
another `--setup-coef` in the same root is refused (exit 5) and names the
root that C belongs in; `--replace` overwrites on purpose. A root's name says
only which C it is meant for (a date's first build in a root is never
refused): when the name gives the root to the very C that was refused there
(e.g. the default root holding a C = 0 build), the message says which C the
root holds and names a root beside it, `<root>_c<C>`, never the same root.
Give a C = 0 or a refitted C its own `--output-root` (above, "a build with
C = 0"). A full rebuild with another `--window` or `--utc-hours` replaces the
date's product (the run says so): build such a recheck in its own
`--output-root` too.

The code is tested here under the station's oldest libraries too (Python
3.8, numpy 1.17.4, matplotlib 3.3.4, pandas 1.1.5). The waterline maps on the
photos need the `daily_elevation_map.py` that falls back to `plt.get_cmap`
on matplotlib 3.3/3.4 (the version merged with the waterline work).

## When the October GCP file arrives

Save it as `/mnt/I2Rgus_Data/Chelsea_calibration/2024-10-23_Marconi_Extrinsic_Targets_<cam>_xyz.csv`
(no header: `num,E,N,Z`, UTM 19N / NAVD88, like
`calibration/2025-11-13_Marconi_Extrinsic_Targets_c1_xyz.csv`), or write its
path in the `oct_gcps` row of `surveys.csv` (one file per camera is fine: they
are merged, a target in both once, for the comparison and the plan view).
When the October date is switched on (`enabled` = 1 in `survey_dates.csv`, or
`--force-disabled`), the comparison runs; on a date already built,
`--date 2024-10-23 --force-disabled --steps maps,compare` (about a minute)
makes it and redraws the plan view with the points. Only the GCPs inside
the DEM's cells are compared (targets on the berm and the upper beach lie
above the waterlines), so expect few points, each off by up to ~0.1-0.2 m
from the slope inside a 2 m cell. It is PARTLY-CIRCULAR: the 2024-10-23 calibration was solved
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
  transfer is fitted on their overlap; without the ADCP the forcing step
  fails, and the dry run says so). 20-26 Oct 2024 lies 44-51 days before the
  ADCP record. A period before the record is judged by the BACKWARD
  extrapolation test (fit 2025-01-01..03-10, predict 0-23 days before it):
  RMS 0.095 m, bias +0.037 m (the transfer above the ADCP: elevations that
  much high), written as each level's `sigma_m`. October lies further out
  than that test reaches, so its error is NOT measured: the forcing report
  (`FURTHER THAN TESTED`), forcing.json (`transfer_sigma_measured: false`)
  and the README say so. The error is shared, not averaged away: the mean
  error of a whole left-out week of the overlap (13.5-18 UTC) ranged
  -0.09..+0.10 m (RMS 0.05 m), and a week-long DEM carries one such mean in
  every frame. A README written by a newer survey_products.py from a
  forcing.json built before it (`--steps maps,compare` after a `git pull`)
  works this out from forcing.json's own record of the tests; the
  `sigma_m` in water_level.csv and the forcing report change only when the
  forcing step is rebuilt (which rebuilds the detection after it).

Run `--date 2024-10-23 --force-disabled --dry-run` first: it lists each
forcing file found or NOT FOUND.

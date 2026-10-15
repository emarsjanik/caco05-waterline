# Survey dates: what survey_products.py builds, and against what

Two tables drive `survey_products.py`. Edit them, not the script, to add a
date, switch one on or off, or give a survey file that has arrived.

```
python3 survey_products.py --date 2025-01-23            # one date
python3 survey_products.py --date 2025-01-23 --dry-run  # what it would do: photos per day, files, labels, run time
python3 survey_products.py --all                        # every enabled date, then the summary
python3 survey_products.py --summary                    # one table + one figure across the dates built
python3 survey_products.py --date 2024-10-23            # BLIND: a DEM only (below, "23 Oct 2024: BLIND")
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
dem,maps,compare`. The DEM's day test (`--dem-max-day-offset`, 0.15 m;
docs/VALIDATION.md section 2.3a) leaves out a camera-day whose waterlines
sit, on median, more than 0.15 m from the DEM of the other kept days. It is
iterative: the day furthest over the limit goes first and the others are
measured again without it, until none is over. A day is left out only when
the test can judge it: 3 or more testable days kept, and where it is
compared its DEM holds two or more other kept days, or one (a lone partner)
that the rest confirm. A day it cannot judge (two days that only see each
other, the last two days) keeps the first pass's verdict (over the limit
there: left out; within: kept) and is tried again in later rounds; a day
left out is taken back, closest first and one at a time, when it is within
0.15 m of the kept days and its return puts no kept day over. One within
the limit whose return would put a kept day over stays out: the README
calls it a "refused take-back", never "off the rest". A day testable in the
first pass with too few comparable samples left at the end is kept
untested. When the final DEM leaves out or keeps over the limit a day it
could not judge, or keeps one untested, the build log and the README print
the day test's WARNING. A second WARNING flags a camera whose day offsets
run one way with time (5 or more days, a range over 0.15 m, |Spearman rho|
> 0.8, and up to 8 days a rho that random orders reach at most 2 % of the
time: 1.00 over 5 days, 0.89 over 6): real beach change or a drift
(pointing, lens, detection), whose far end the test then trims. A single
day left out at the first or last of a camera's days cannot make that
trend alone: the other days must run the same way (|rho| > 0.8), and when
they would not warn on their own the WARNING says the trend "rests on"
that day. With the days in two groups of equal weight (say three at 0 m
and three at +0.6 m) the first day left out decides which group stays and
nothing warns: read a day left out at little over the limit, or one
camera's days left out together, with care. The README's WHAT WAS USED
names every camera-day left out with how (left out, not judged, refused
take-back), its offsets and why; its KNOWN CAVEATS repeat the WARNINGs.
The test is not run on fewer than 3
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
`--setup-fitted-to`; a calibration variant (`--eo`) in a root not named for
it, a missing `--eo` file, or one whose pre-registered acceptance FAILED; a
blind date over a folder holding an earlier comparison; `--dry-run` says when
a real run would be refused); exit
6: not started, a usage error (an unknown option or a bad value, e.g.
`--steps`, `--window`, `--utc-hours`, a negative `--setup-coef`, `--eo` without
`--calibration-variant` or with `--all`), never 2. A
row with no survey file yet is not a missing part, and the held survey of a
BLIND date is never looked for (a blind product with everything else built is
`complete`). With
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
| 2024-10-23 | enabled, **BLIND: DEM only** (user's request, 13 Oct 2026) | c1+c2 CACO03, 20-26 Oct | A: CACO03_<cam>_20241023_EO (solved from that day's GCPs, the ones held); B: GCP-free (`--eo ... --calibration-variant gcpfree`, its own root) | Chatham transfer / WIS | GCPs held by Dr Chris Sherwood (USGS): never on the NUC, never compared here | none here (`blind (DEM only)`); his check of A is PARTLY-CIRCULAR for the geometry, of B INDEPENDENT |
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
across a 1:10 foreshore. The GeoTIFFs are tagged EPSG:6348 (NAD83(2011) /
UTM 19N, the frame of the 2025 lidar files and the 2026 RTK), the default
since Oct 2026; until then they were tagged EPSG:32619 (WGS 84 / UTM 19N), a
nominal tag a GIS places ~1-1.5 m off the surveys (nothing is transformed
either way; a product built before then keeps its tag until its DEM step
reruns, and `--geotiff-epsg 32619` still writes the old one, said nominal). The 2025 lidar files state no vertical datum: NAVD88 is assumed (the
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

* features matched against the reference frame on 40 or more points spread
  over the frame (pointing_check.py) are good to ~0.1 deg on the views tested
  (CACO05 c2 2026 against another day; a match bunched in a small part of the
  frame reads roll less well: its spread is in `inlier_spread_px`, and the
  README says when a deciding match is bunched), and a turned match is the scene
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

Three choices the rule leaves to the owner, stated rather than decided
(docs/VALIDATION.md 4.6): (1) a storm whose every day nothing could read
(the drops or spray breaking every horizon fit and match) is KEPT: no day of
it is left out on its own readings, and leaving such days out on the wave
record alone is the owner's call (the 1.5 m threshold is the DEM's own,
never measured as a sign that a camera turned); (2) when the reference day
itself is the one that moved, the build is partial with a warning that gives
the `--window` building on the other pointing (in another root), and it is
not refused until the owner chooses; (3) a change of exactly the 0.3 deg
tolerance is DIFFERENT POINTING on the days a clear horizon reads it just
over 0.300, a 'small change' on the others: days are lost at the boundary.

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
python3 survey_products.py --date 2024-10-23 --dry-run   # BLIND, calibration A: forcing files, held check
python3 survey_products.py --date 2024-10-23             # A, default root
python3 survey_products.py --date 2024-10-23 --share     # -> survey_products/share/2024-10-23_A_official_calibration.tgz
# B, LATER (after the owner has pre-registered the outcome rule; below, "23 Oct 2024: BLIND"):
#   gcpfree_eo.py --preregister ...; git commit AND git push PREREGISTERED.txt BEFORE any gcpfree_eo.py run
#   (a trial included); its sha256 and commit to Chris before A's DEM is shared
#   gcpfree_eo.py writes calibration/gcpfree/runs/<run id>/ (a B file only for a camera it delivers):
#   commit and push it in one commit, and send Chris each B file's sha256, the moment the run ends
python3 survey_products.py --date 2024-10-23 \
    --eo c1=calibration/gcpfree/runs/<run id>/CACO03_c1_gcpfree-oct2024_EO.yaml \
         c2=calibration/gcpfree/runs/<run id>/CACO03_c2_gcpfree-oct2024_EO.yaml \
    --calibration-variant gcpfree \
    --output-root /mnt/I2Rgus_Data/survey_products_gcpfree --dry-run   # then without --dry-run
#   (built exactly as A was: the same options as A's build, --a-root <A's root> if A is not in the default
#   root; the first build writes calibration/gcpfree/runs/<run id>/B_BUILD.txt: commit it at once)
python3 survey_products.py --summary --output-root /mnt/I2Rgus_Data/survey_products_gcpfree
python3 survey_products.py --date 2024-10-23 --share --output-root /mnt/I2Rgus_Data/survey_products_gcpfree
#   -> survey_products_gcpfree/share/2024-10-23_B_gcpfree.tgz
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
`--output-root` too. One calibration variant per root, too, and stricter: a
`--calibration-variant NAME` build goes only in a root named `<root>_NAME`
(e.g. `/mnt/I2Rgus_Data/survey_products_gcpfree`), a table build never in a
root holding a variant (exit 5 either way, `--replace` does not override it:
the remedy is the right root, so a root's name always says which calibration
built what is in it).

The code is tested here under the station's oldest libraries too (Python
3.8, numpy 1.17.4, matplotlib 3.3.4, pandas 1.1.5). The waterline maps on the
photos need the `daily_elevation_map.py` that falls back to `plt.get_cmap`
on matplotlib 3.3/3.4 (the version merged with the waterline work).

## 23 Oct 2024: BLIND (DEM only), two calibrations

The GCPs surveyed on 23 Oct 2024 are held by Dr Chris Sherwood (USGS) for an
independent check of our DEM of that date. **Never put them on the NUC or in
this repository.** `survey_products.py` treats the date as BLIND
(`BLIND_DATES`): it builds the forcing, pointing, waterlines, DEM, GeoTIFFs,
maps and README only. The date's survey rows never reach a file look-up
(nothing named like `2024-10-23*xyz*.csv` is globbed, stat-ed or opened, in
any folder, whatever `surveys.csv` says), the compare step does not run, the
plan view draws no survey. As a backstop, a runtime guard (an audit hook)
refuses, in survey_products.py's own process, any Python-level open for
reading of a `.csv`/`.txt`/... named like a survey of that date (by its name
or the name a link points to), and no step it starts is handed such a file
(each argument is checked; the photo folders a step is given were scanned
first in survey_products.py's own process, where a photo-named link to such
a file stops the build). The hook does not see readers that open files in C
(pandas' default CSV reader, OpenCV) or the step scripts' own processes:
those are safe because no look-up ever hands them a held survey's path. The
README and provenance say `blind: the survey is held by Dr Chris Sherwood
(USGS) for an independent check; nothing in this build read it`, and what his
check can test with the calibration used. `--summary` lists the date as
`blind (DEM only)`. A folder holding an earlier comparison of the date is
refused (exit 5): build in a fresh root. The date is unblinded only by removing
its `BLIND_DATES` entry (a commit), after he has reported.

Two DEMs are delivered, each in its own root (one calibration per root):

* **A, the table's calibration** (`CACO03_<cam>_20241023_EO.yaml`), solved
  from those very GCPs: the camera geometry is fitted to his points, so his
  check of A tests the water level, the detection and the setup, not the
  geometry (PARTLY-CIRCULAR for the geometry).
  `python3 survey_products.py --date 2024-10-23`
  (default root `/mnt/I2Rgus_Data/survey_products`).
* **B, a calibration made without them**: the 23 Jan 2025 GCP solution of the
  same mount (`CACO03_<cam>_20250123_EO.yaml`) turned back to the October
  pointing by matching Jan and Oct photos (c2 directly, or through the Dec
  2024 photos; or, if Chris sends CACO03 photos of 25-29 Aug 2024, the Aug
  2024 solution turned by matching Aug and Oct photos), its camera position
  that of the solution it was turned from. Computed on the NUC by the
  PRE-REGISTERED GCP-free procedure, `gcpfree_eo.py` (its docstring is the
  procedure), never against A's numbers: it refuses to read calibration A,
  any GCP or target file and the calibration history (by name and by the
  name a link points to), and the Jan, Aug and lens files unless their
  content is the archived one. Expected (stated before any run): c1 NOT
  DELIVERED from Jan photos (its Jan and Oct views share no fixed feature),
  c2 uncertain; Aug photos are the route most likely to work.
  The README traces each file to the GCPs it rests on: INDEPENDENT of his
  GCPs only when the procedure made it (below).

### B, step by step (later, on the NUC)

1. **The owner pre-registers the outcome rule** for the likely case that
   only one camera passes the procedure's gates (c2): `both-or-none` (no B
   unless both cameras are delivered: Chris gets A alone, with the reports)
   or `per-camera` (a c2-only B DEM; its README names c1 NOT DELIVERED with
   the procedure's report). Nobody else chooses it, and it is chosen before
   any run:
   `python3 gcpfree_eo.py --preregister --outcome-if-one-camera <both-or-none|per-camera> --decided-by "<name>"`
   then `git add calibration/gcpfree/PREREGISTERED.txt && git commit`, **and
   `git push` that commit to GitHub at once, before ANY run of
   `gcpfree_eo.py`, a trial included**: the time GitHub receives the push is
   the rule's public timestamp (a commit's own date is set by the computer
   that made it, so it proves nothing on its own); B's README and MANIFEST
   cite the commit and say whether a remote branch holds it, and Chris can
   check on GitHub when it arrived.
   The file also PINS THE PHOTOS: the name, size and sha256 of every Oct, Jan
   and (c2) Dec photo a run selects from the station's three folders; a run
   that selects other photos (one removed, replaced or added after a look) is
   refused. It is the only pre-registration there is (`calibration/gcpfree/`
   next to `gcpfree_eo.py`; there is no option for another): EVERY run, a
   trial (one camera, other `--roots`) too, refuses (exit 2) without it,
   uncommitted, changed since its commit, not pushed (its commit on no
   remote-tracking branch: `git branch -r --contains <commit>` names none; the
   script's own messages print the `git push`), naming another `gcpfree_eo.py` or
   module sha256, pinning no photos, or once any run was made under another
   version of it; so the rule is committed before this code solves any photo.
   Every run, trials included, is recorded in `calibration/gcpfree/runs/`
   (commit each folder), and B's README lists them all with their start
   times. **A trial is a full result (every verdict and angle): never run one
   on the station.** A run is also refused while a commit outside the
   branch's history (another branch, a stash, the reflog: a branch made from
   before the pre-registration, a reset) wrote a pre-registration or run
   record the branch does not hold. `survey_products.py` accepts a B only from
   `calibration/gcpfree/runs/` of this repository, made under this file as it
   was committed before the run.
   **Send Chris the procedure's sha256 before A's DEM is shared** (ideally
   before it is built: `sha256sum gcpfree_eo.py`, = `GCPFREE_PINS[-1]`), and
   the PREREGISTERED.txt sha256 with its commit hash as soon as it is
   committed and pushed (also before A's DEM is shared): he then knows both
   were fixed before he saw A.
2. `python3 test_gcpfree_eo.py` (synthetic, in a temporary folder, ~25-30 min
   here under load; allow 30-40 on the NUC)
   and `python3 gcpfree_eo.py --dry-run` (the pre-registration, the photos
   per day, the Jan 2025 EO against the visible horizon on 23 Jan; of the
   Oct photos it only lists the names and checks the frame size: nothing is
   computed from them; writes nothing; it says `a delivery run would be
   REFUSED (exit 2): ...` when a run would be). Run the dry run early (it needs
   no pre-registration): a missing photo folder or a camera with no Oct or no
   18-23 Jan 2025 photo refuses every delivery run, and is better known
   before anything else.
3. **The delivery run**: `python3 gcpfree_eo.py` (the station's three photo
   folders, both cameras; ~10-30 min a camera). It writes
   `calibration/gcpfree/runs/<UTC run id>/`: per camera its report
   (`CACO03_<cam>_gcpfree-oct2024_report.json/.txt`), a B file
   (`CACO03_<cam>_gcpfree-oct2024_EO.yaml`) only for a camera DELIVERED
   under the rule, and `run.json`. **Commit that folder at once, in ONE
   commit (never amended), push it, and send Chris each B file's sha256 the
   moment the run ends, before the B DEM is built and before anything is
   compared with the A DEM.** `survey_products.py` builds no B whose file, report, `run.json` or
   `run_report.txt` is uncommitted, changed since, or touched by more than one
   commit, nor while any run folder there (a trial, a run cut off) is
   uncommitted. (The A
   DEM is built first, before this run: nothing of it reaches the procedure,
   whose code is pinned by sha256 and whose outcome rule is committed before
   any run; B's README states when A's DEM was built against the pin, the
   pre-registration and the run.) No run ever writes over another; a second
   delivery run
   is refused (B is fixed by the first), except one with `--aug-dir` (below).
   Bad input is refused (exit 2), never turned into a verdict: no Oct or Jan
   photo for a camera, a missing photo folder.
4. Build B with the cameras it delivered (`--eo c2=...` alone for a c2-only
   B: a camera is left out only when its report in the same run folder says
   NOT DELIVERED), `--calibration-variant gcpfree`, in its own root (commands
   above), EXACTLY AS A WAS BUILT, only the calibration different: the
   committed `surveys/*.csv`, every result-changing option as A's build
   recorded it (A's provenance: give `--a-root` when A is in another root),
   the search envelope A's detection used (path and sha256) and A's forcing
   (inputs, then its water level and waves byte for byte). Anything else is
   refused (exit 5, `NOT the pre-registered B build`). Its first build writes
   `B_BUILD.txt` (its root and command) into the run folder: **commit it at
   once**; the B is then built and shared from that root only, with that
   command. `--share` it.
5. Only then compare A and B (`dem_change.py`): a report, nothing is tuned
   to it.

**If Chris sends Aug 2024 photos** (named as the station names them,
`<epoch>.<Day>.<Mon>.<dd>_<HH>_<MM>_<SS>.GMT.<YYYY>.CACO03.<cam>.timex.jpg`,
full frame): `python3 gcpfree_eo.py --aug-dir <folder>` is the ONE more
delivery run allowed, and only when the first left a camera NOT DELIVERED.
It solves only those cameras; a camera delivered in the first run keeps that
run's verdict and B. Fixed now, before any run: when the first run delivered
no camera, this run's B is the B (`--calibration-variant gcpfree`); otherwise
it is a SECOND, separately labelled product (`--calibration-variant
gcpfree-aug`, root `/mnt/I2Rgus_Data/survey_products_gcpfree-aug`, archive
`2024-10-23_B2_gcpfree_aug.tgz`) and the first run's B stays the delivered
one. The B file and its report name the variant; `survey_products.py`
refuses any other.

### What the code enforces, and what rests on the owner

B is independent of Chris's GCPs only if nothing about it was chosen after
anyone saw A or a result. Part of that the code enforces; part of it no
software can prove, because the computer that runs the procedure is the
owner's: its clock dates every commit and run, and its repository can be
re-cloned or have history deleted. B's README (BLIND section) and the share
MANIFEST say both, in these words.

* **Enforced** (`survey_products.py` and `gcpfree_eo.py`): B is made only by
  the pinned procedure (gcpfree_eo.py, its five modules and constants, by
  sha256: `GCPFREE_PINS`); every run, a trial too, is refused unless this
  repository's `calibration/gcpfree/PREREGISTERED.txt` (the outcome rule and
  the photos, by sha256) is committed, unchanged and pushed (its commit on a
  remote-tracking branch: the computer's record that a push was made; when
  GitHub received it, only GitHub shows), and every run is recorded in
  `calibration/gcpfree/runs/`; the first delivery run fixes B;
  B's file, report and run record are committed once and unchanged; B's DEM
  is built as A was (the committed tables, A's options, A's envelope, A's
  forcing), in one root; the procedure never opens calibration A, a GCP or
  target file or the calibration history, and no build looks for the held
  survey; B's ORDER (the pin, the rule, A's DEM first built and last
  rebuilt, every run, B's commit, with NOTEs) is read from git and the run
  records.
* **Resting on the owner**: PREREGISTERED.txt committed AND pushed to GitHub
  before any run of `gcpfree_eo.py`, a trial included; no copy of the
  procedure run anywhere else (another clone, a copied script) before it; no
  run folder, branch or reflog entry deleted, the repository not re-cloned;
  no B file edited before its commit; the computer's clock right.
* **What answers those, from outside this computer**: the sha256s of the
  procedure and of the pre-registration (with its commit hash), sent to Chris
  before A's DEM is shared; each B file's sha256, sent the moment its run
  ends; and the time GitHub received the pre-registration commit (the
  repository's activity on GitHub), which the owner cannot set. B's README
  cites that commit (`PRE-REGISTRATION COMMIT: ...`) and says whether a
  remote branch this computer has fetched holds it.

### What makes a B file the pre-registered procedure's

`survey_products.py` pins the procedure (`GCPFREE_PINS`: the sha256 of
`gcpfree_eo.py`, of the five modules it computes with, and of its
pre-registered constants as a run records them in effect). A turned
calibration is called INDEPENDENT of his GCPs only when all of these hold;
otherwise its README says `NOT the pre-registered procedure` (why) and labels
his check PARTLY-CIRCULAR, never INDEPENDENT, and it is never built or shared
under the procedure's names:
* its `# Made by gcpfree_eo.py (sha256 H) ...; report <file>.json` line
  names a pinned `gcpfree_eo.py`, and its run's modules and constants are the
  pinned ones (a run whose constants, e.g. its photo folders, were changed in
  the running process is not);
* that report lies beside it (the same run folder) and records the same B
  (angles and position), the same verdict, DELIVERED, the B file's sha256,
  a delivery run (not a trial); `run.json` beside it records the file;
* its run folder is `calibration/gcpfree/runs/<id>` of THIS repository, the
  run was made under this repository's `PREREGISTERED.txt` with the sha256
  and the commit it has now (not changed or committed again since), from the
  station's three photo folders, and it is the first finished delivery run
  there (or the one allowed second run, with Aug photos), no run's record
  (`run.json`) there having been removed or changed since its commit, every
  run's record there committed (a trial, a run cut off), and the photos the
  run used are the ones its pre-registration pins;
* no commit outside this branch's history (another branch, a stash, the
  reflog) wrote a pre-registration or run record this branch does not hold;
* the B file, its report, `run.json` and `run_report.txt` are committed in
  this repository, unchanged, in one commit, each touched by that commit
  only;
* its `# Variant: ... --calibration-variant <name>` line is the report's;
* its `# Solved by gcpfree_eo.py: <reference>.yaml rotated by <angle> deg`
  line is true of its numbers (0.01 deg), its camera position is that of the
  reference its `# Position (x, y, z) kept from` line and report name (1 mm),
  its `pan change <deg> deg from <Jan reference>.yaml` is true of the files
  (0.01 deg), and each reference has the content (sha256) the report recorded
  and, for a GCP solution known here, its archived content.
What these checks cannot see: a B edited BEFORE its first commit with its
report and `run.json` rewritten to agree, by an edit that keeps the `rotated
by` magnitude (a small tilt, say); a copy of the procedure run outside this
repository (another clone, a copy of the script) before the pre-registration,
or a whole re-clone of the repository after a first run; a branch or stash
removed together with its reflog entries; a trial or a run cut off whose
folder was deleted before it was committed. Only the sha256s sent to Chris
beforehand (the procedure's before A, the pre-registration's when committed,
each B file's when its run ends) answer those; hence steps 1 and 3.
A calibration turned by hand (e.g. with `estimate_eo_rotation.py`, its
`rotated by` line and an `# acceptance:` line added) is traced the same way
but is never INDEPENDENT on a blind date: such a turn could have been chosen
after seeing A. The names `gcpfree` and `gcpfree-aug` are the procedure's:
no other calibration is built under them (exit 5), all of a B's cameras come
from one run, and a run's B is built only as the variant it names.

The procedure's verdict is its first comment line, `# CALIBRATION B
(GCP-FREE) ..., camera <cam>: DELIVERED (...)` (= PASS) or
`DELIVERED-WIDER (<why>)` (= PASS-WIDER: built, the reason on every output);
`NOT DELIVERED` (= FAIL) writes no file, and a file saying so is refused
(exit 5); the table's calibration is never built in its place under B's
name. On a blind date a file stating no acceptance is refused too.

B's README gives its two STATED uncertainty terms, neither a validated
bound on real photos (both come from synthetic tests): the rotation, max(2 x
the solve's spread, 4 px) on the beach part of the frame, with the DEM effect
of that rotation taken as a pure tilt, per range band, from the date's own DEM
and calibration (the synthetic sets' rotation errors were almost all roll, so
their millimetre DEM errors say nothing of a tilt: 4 px of tilt moves
waterlines metres at 200-400 m); and the inherited camera position, the
largest DEM effect of five synthetic position-inheritance cases (c1 DEM |dz|
median 0.011 m, p90 0.026 m; c2 0.018 m, 0.086 m; the largest of five cases,
not a bound; the cases were written after A's rows had been seen, and c1's
comes from a -0.55 m height case that equals c1's A minus Jan offset: the
EXPOSURE statement says so), while the mount's GCP solutions of Aug 2024 and
Jan 2025 put the camera 0.36 m (c1) / 1.23 m (c2) apart. It also says the run and its
folder, the outcome rule, the pre-registration (its path, sha256 and
commits), the photo folders and the B file's commit; the procedure's EXPOSURE
statement (what its makers had seen of A, and when); the ORDER of the pin, the
pre-registration (the commit that set the rule in force, and every commit of
the file with its rule when there are several), A's DEM (when it was FIRST
built, kept across rebuilds of its DEM step, and when last rebuilt; from A's
root: the root B's is named after, or `--a-root`), the B run and B's commit,
and every run of the procedure in `runs/` with its start, with a NOTE when
the rule was changed, when A's DEM came before the rule in force was
committed, or before B's commit, when a run started before the rule was
committed, when a trial was made before B's run, when a run never finished,
when no A product was found (the order is then not known), or when B was
built with other code than A; each camera's search envelope with its sha256
(the file A's detection used: Chris can check it against his copy of the Jan
2025 lidar); and that the pointing step's
sea-horizon check of the October photos against B is independent of the
photo matches.

Refused before anything is written (exit 5): a `--calibration-variant` build
in a root not named `<root>_<variant>`, a table build in a root holding a
variant (and the reverse), an `--eo` path that does not exist (it is taken
as given: never a file of the same name elsewhere, never the table's), an
acceptance FAIL, a camera left out that the procedure delivered, B files of
two runs, a B under another variant name than its run's, a B whose run
folder is not committed yet (commit it first), a B not built as A was or in a
second root (step 4), a photo folder holding a photo-named file that is (or
links to) a file named like the held survey (refused in the photo scan,
before any step, a dry run too). A B file changed since B was built is a
rebuild refused without `--replace` (stamps record each calibration's
sha256; a file restored byte for byte, e.g. by git, rebuilds nothing).

**Telling A and B apart once shared.** B's DEM and map files carry the
variant in their names (`2024-10-23_gcpfree_dem.tif`, `_spread`, `_count`,
the page; `maps/2024-10-23_gcpfree_waterlines_plan.png`, `_c1`, `_c2`); A's
are `2024-10-23_dem.tif` and `maps/2024-10-23_waterlines_plan.png`. Every
GeoTIFF's description (ImageDescription) and the plan view's title name the
calibration: the table's or the variant, each camera's EO file and sha256
prefix, and his check's label. `--share` packs each into one archive
named by its calibration, `2024-10-23_A_official_calibration.tgz` and
`2024-10-23_B_gcpfree.tgz` (+ `.sha256`), with README, provenance, GeoTIFFs
(EPSG:6348, NAD83(2011) / UTM 19N + NAVD88), the DEM page, the maps, the
calibration files used (the detector reads every lens as
`CACO05_<cam>_20240801_IO.yaml`; a CACO03 product also carries
`CACO03_<cam>_20240801_IO.yaml`, the file `gcpfree_eo.py` solved with, and its
MANIFEST says they are the same lens, the CACO05 file holding its numbers
rounded), for B the procedure's reports, `run.json`, `B_BUILD.txt` and
`PREREGISTERED.txt` (and the report of a camera it did not deliver), and a
MANIFEST of every file's sha256 (for B also its procedure line with the run
folder and commit, the stated uncertainty, the EXPOSURE statement and the
ORDER); never written over (`--replace`). Refused: a B that is no longer the
procedure's, a B from another root than its `B_BUILD.txt` names (or with that
file uncommitted), and any product (A or B) whose calibration file is no
longer the one it was built with (its sha256 against the provenance and the
detection's stamp): an archive never pairs a DEM with another calibration
file.

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

Run each command with `--dry-run` first: it lists each forcing file found or
NOT FOUND, each calibration and what it was solved from, the held check per
camera, and whether a real run would be refused.

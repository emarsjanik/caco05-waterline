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
were made with the same settings; `--force` rebuilds, `--steps compare`
(etc.) runs only some steps. Every skip is printed.

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
| 2026-09-27 | enabled | c1+c2 CACO05, 24-30 Sep, the live archive | CACO05_<cam>_20251113_EO-CV | GNSS-R / live waves | RTK check shots of 29 Sep | **CIRCULAR with C = 0.037** (see below) |

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
* a GCP survey on the calibration's own day -> PARTLY-CIRCULAR;
* **a setup coefficient fitted to this survey -> CIRCULAR.** C = 0.037
  (`waterline_timex_cron.sh`) is the median per-frame C of the waterlines
  lying on the 2026-09-29 RTK transects, so the 2026-09-27 comparison with
  those shots cannot test the setup (nor the overall level it sets) and is
  labelled CIRCULAR, whatever the table says. Its spread and its dependence
  on elevation still say something.
  Note that the fit is not exact either: it took the RTK elevation where each
  line lay WITHOUT setup, but a line given the setup is re-projected landward
  onto higher beach, which closes only part of the gap. `compare_rtk.py` on
  the 18 frames of 29 Sep - 5 Oct lying on the transects gives RTK - waterline
  +0.33 m without setup and still +0.13 m with C = 0.037 (median setup applied
  0.26 m; Oct 2026): about half the applied setup is still missing, so every
  product made with C = 0.037 may read low by roughly 0.5 x its setup (each
  README gives the figure for its date). The synthetic test beach shows the
  same geometry: a line mapped onto the still-water plane reads low by about
  0.65-0.7 x the setup, not the whole setup. A C refitted with the lines
  re-projected would be the cure; that is the cron's setting, not this script's. For an independent check,
  build the date with a C fitted elsewhere (e.g. from repeat crossings,
  `dem_from_contours.py --fit-setup` on contours without setup, which gave
  0.03-0.04): `--setup-coef 0.035`; the RTK comparison is then labelled as
  the table says. The 2025 lidar comparisons are not affected: C was not
  fitted to them.

`survey_compare.py` checks again (calibration notes, the envelope, the
`provenance.json` next to the DEM) and prints its own downgrade warning
when the label looks too good.

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

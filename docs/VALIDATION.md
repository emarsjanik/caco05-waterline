# Validation: waterline filter, DEM page, week-to-week change, GNSS-R guard and survey-date products

This file backs every claim made by the pull request from `claude/zen-hypatia-voks3l`
(waterline products plus survey-date products). For each claim it says what was
measured, on which data, the number, how to reproduce it from this repository,
and how to check it on the station computer (the NUC) with the real data.

## How to read this file

Every number carries one of these tags. Numbers from different tags are never pooled.

| tag | meaning |
|---|---|
| **[REAL]** | Real station data: the live waterlines of 29 Sep - 5 Oct 2026 (`contour_points_timex.csv` rows copied from the NUC), the 29 Sep 2026 RTK check shots, the Jan and Mar 2025 lidar DSMs, the ADCP, Chatham, WIS and NDBC records, and two real photos of 7 Oct 2026. |
| **[SYNTH]** | Synthetic data with a planted truth: the self-tests, a synthetic 24-day station for the cron, and a fixture of photos rendered from the 2025 lidar. These test the code. They say nothing about the real beach. |
| **[SCRATCH]** | Real data analysed with scripts that are **not** in this repository (the Oct 2026 setup-coefficient refit and its adversarial check). The numbers are real, but only the parts that repo tools also reproduce can be re-run from here. |
| **[NUC]** | Can only be checked on the NUC. The real January and March 2025 photos, and the 26-28 Sep 2026 rows and photos, are there and nowhere else. |

Where a check runs a self-test, it is in this repository and runs in a few minutes.
"Here" means the development sandbox, where the checks were run on three Python setups:
3.11; 3.8 with numpy 1.24, pandas 2.0 and matplotlib 3.7; and 3.8 with the station's
oldest libraries (numpy 1.17.4, pandas 1.1.5, matplotlib 3.3.4, OpenCV 4.5.1, no scipy).

---

## 1. Fifteen-minute check on the NUC after the pull

```
cd /mnt/I2Rgus_Data/waterline
mkdir -p ~/before_pr && cp dem_intertidal_7day_dem.png elevation_map_c?_7day.png \
    waterline_consistency_c?.png ~/before_pr/ 2>/dev/null   # BEFORE pulling: today's pictures, for your own before/after
git status --short            # tracked files must show no local edits, or the pull stops
git pull                      # main, after this PR is merged
bash -n waterline_timex_cron.sh && echo "cron syntax OK"
python3 waterline_consistency.py --self-test     # ends "SELF-TEST PASSED" (90-110 s here)
python3 dem_figure.py --self-test                # ends "SELF-TEST PASSED"
python3 daily_elevation_map.py --self-test       # ends "SELF-TEST PASSED"
python3 test_survey_products.py                  # ends "all survey_products tests passed"
python3 test_survey_compare.py                   # ends "all survey_compare tests passed"
python3 test_historical_forcing.py               # ends "66/66 checks passed"
cat archive/gnssr_record_start.txt               # note "YYYY-MM-DD N"; see section 2.5
```

`view_reproject.py --self-test` needs scipy, which does not import on the NUC. That file
is unchanged by this PR, and the survey products use only its `ground_to_pixel`, which
does not need scipy.

After the next cron run (12:30, 19:25 or 20:55 local):

```
grep -nE "run start|waterline consistency|consistency filter:|DEM written|window DEM|BEACH CHANGE|reference rebuild|No change map|WARNING|ERROR|run end" \
    logs/timex_cron.log | tail -60
cat archive/gnssr_record_start.txt               # N must not have fallen
```

Then compare `dem_intertidal_7day_dem.png`, `elevation_map_c1_7day.png`,
`elevation_map_c2_7day.png` and `waterline_consistency_c?.png` with the copies in `~/before_pr/`.
Section 2 says what should have changed.

---

## 2. The station's daily products (cron)

### 2.1 Consistency filter: what it drops and what it keeps

**Claim.** Lines (or parts of lines) that are out of order with the other lines of the same
camera, by the physics of a beach, are left out of the maps, the ground points and the DEM.
Real beach change, a storm cut included, is kept.

**Rule in one paragraph** (full text: `waterline_consistency.py` docstring). In one image
column, a higher water level puts the waterline further landward, i.e. lower in the photo,
and this direction is fixed by the camera geometry, not voted on by the data. Per camera and
32-px column bin, each line is compared with the lines of its day and the 3 days either side
(weight `exp(-|dt|/1 day)`) through a monotone fit of row against elevation. A bin fails
when the elevation residual exceeds `max(0.75 m, 4 x robust scatter)`. A failing line is
judged again against the lines **before** it and those **after** it in time, and goes only
when both sides contradict it, so a storm cut (which moves the lines on one side only) is
kept and counted as "beach change". The most seaward lines are judged one at a time against
the others, so a stray out on the water cannot become the end of the fit.

| what | data | number | reproduce here | check on the NUC |
|---|---|---|---|---|
| Systematic c1 error (high-tide lines drawn seaward on the right half of c1) | [SYNTH] self-test, real calibration | 16,254 of 16,384 wrong points dropped (99.2%); REVERSED warning on c1 columns 1408-2431 | `python3 waterline_consistency.py --self-test` | same command |
| 0.4 m erosion of the upper beach | [SYNTH] self-test | 8 of 161,970 honest points dropped (0.00%) | same | same |
| Storm cut 0.8 m in 6 h, setup +0.15 to +0.74 m | [SYNTH] self-test | day before: 0 of 22,722 honest points dropped; all honest lines 0 of 228,100 | same | same |
| Cuts of 0.8 m and 0.5 m on the surveyed beach, and one on the record's first day | [SYNTH] self-test | 0 honest points dropped after each cut; 13, 14 and 3 lines kept as beach change | same | same |
| Stray at the window's lowest level, 38-40 m out (c2) | [SYNTH] self-test | 136 of 136 points dropped | same | same |
| Two strays in consecutive frames | [SYNTH] self-test | 136 of 136 and 119 of 119 dropped; honest lines 0 of 85,312 | same | same |
| Stray 30 m out, 3 cm below the window's lowest line | [SYNTH] self-test | 204 of 204 dropped; honest 0 of 73,356 | same | same |
| Megacusps (15 m every 140 m) | [SYNTH] self-test | 3,823 of 117,424 honest points dropped (3.26%; the test's limit is 5%) | same | same |
| Whole cron on a synthetic 24-day station (571 frames, 1.28 M points, 31 wrong c1 lines and 3 c2 strays planted) | [SYNTH] cron harness (scratch, not shipped) | all 31 wrong c1 lines hit: 20,214 of ~21.2k defect points dropped; all 3 c2 strays hit: 1,366 of 1,704; 0 honest points dropped | not in the repo | - |
| Real 29 Sep - 5 Oct 2026 contours (393,668 points, 198 lines), setup C = 0.037 added as the cron adds it, before (main `af8e1bd`) and after (this PR) | [REAL] | **main**: c1 22 lines dropped, 19 trimmed, 44,768 of 269,694 points (16.6%); c2 2 dropped, 24 trimmed, 7,122 of 123,974 (5.7%). **this PR**: c1 20 dropped, 17 trimmed, 41,738 (15.5%), 12 lines kept as beach change; c2 2 dropped, 28 trimmed, 8,163 (6.6%), 6 kept as beach change | `python3 waterline_consistency.py <contours> --output qc.csv --report rep.csv --plot wc --image-dir <photos> --plot-days 7` | the cron log line `waterline consistency: c1: ... ; c2: ...` |
| Station log before this PR | [REAL] station log, 8-9 Oct 2026, reported by the owner | c1 drops ~10.8% of its points, mostly high-tide frames placed out on the water (residual +1.4 to +4 m) | - | compare with the log line after the pull |

What changes in the log: the summary line now also counts `N line(s) out of order with one
side in time only, kept as beach change` and the column bins judged. A column whose rows run
against the water level gets one `consistency filter: ... REVERSED ...` WARNING line.

**Known limits** (all [SYNTH] unless marked):
* A single stray about 80 px out at the window's lowest level is kept; a once-a-day stray
  pattern is caught 24% of the time; a scattered 8-stray week 33%.
* After a cut, a line that sits wrongly LANDWARD in columns the later lines do not reach
  (c1's high tides fall below its search-envelope floor) is kept as "beach change" until a
  later line reaches its level there.
* [REAL] 175 c2 and 67 c1 points that their bin's fit puts beyond the threshold are still
  kept: they sit on lines that fail in no bin. Three short dashes remain on the real c2
  7-day map of 29 Sep - 5 Oct (the PR's after figure).
* [REAL] The 2 Oct 11:30-14:30 c2 line ends this PR drops are judged wrong from the
  contours alone (25-50 px seaward of the 3 Oct lines at the same level); no photo of 2 Oct
  was available to confirm. On the NUC: look at `waterline_consistency_c2.png` for a week
  that includes such lines, against the photos of that day.

### 2.2 Elevation maps on the photos

**Claim.** The shaded bands show the 16-84% spread of the lines that agree (3 or more), and
never a box drawn from one stray line. Lines and bands are coloured by the elevation the DEM
uses (water level + wave setup when the contours carry it).

| what | data | number | reproduce here | check on the NUC |
|---|---|---|---|---|
| A stray and two shore lines; two strays 1 px apart; 0.35 m real change; a 245 px band; columns with only two lines | [SYNTH] self-test | band in 0 of 600 columns from the stray; band top at the shore lines (702 vs 700); change band 703-740; band edges move <= 2 px per column; gaps bridged | `python3 daily_elevation_map.py --self-test` | same |
| Real 7-day maps, 29 Sep - 5 Oct 2026, main vs this PR | [REAL] contours; backdrop = the real 7 Oct 16:00 photo | **main** colours by the still-water level (c1 -0.91 to +1.28 m) while the DEM grids water level + setup; a translucent box from stray lines out on the water in c2. **this PR**: colour = water level + setup (c1 -0.66 to +1.55 m, says so in the title and colour bar); no box; the title says `7 days, 6 with lines` for c2 | the commands in section 5 of this file | `elevation_map_c?_7day.png` after the next run, against `~/before_pr/` |

The backdrop used here is a 7 Oct photo linked under a 5 Oct frame name (no photo of the
window was available in the sandbox). On the NUC the backdrop is the clearest photo of the
window's last day.

### 2.3 The 7-day DEM page (the email's picture)

**Claim.** `dem_intertidal_7day_dem.png` (same file name, so the email keeps attaching it) is
now drawn in the beach's own frame: alongshore distance from the cameras, cross-shore
distance seaward, the c1/c2 seam, the camera, "no value" told apart from "blanked: spread >
0.5 m", a repeatability panel and 3-4 cross-shore profiles with their slopes. The DEM grids
themselves change only through the filter.

| what | data | number | reproduce here | check on the NUC |
|---|---|---|---|---|
| Page orientation (sea away from the cameras, NNW to the left) on five shapes, including a flat strip | [SYNTH] self-test | seaward direction within 0.0-0.2 deg of the truth (limits 1-10 deg); the flat strip, which the old rule turned 180 deg, now 0.0 deg off | `python3 dem_figure.py --self-test` | same |
| Orientation of 16 pages (an independent reviewer's check), including the archived real 27 Sep 2026 DEM that was upside down | [REAL] + [SYNTH] | correct signs on all 16; 27 Sep 2026: corr(page x, NNW) +1.000 -> -1.000 | scratch review script | look at the page: sea at the top, cameras at alongshore 0 on the right |
| Real 7-day DEM, 29 Sep - 5 Oct 2026, main vs this PR, cron settings (2 m cells, >= 3 frames, spread <= 0.5 m, Hs <= 1.5 m, day offset 0.15 m) | [REAL] | main 1,577 cells filled of 2,081 crossed; this PR 1,552 of 2,025 (90 blanked for spread, 383 crossed by < 3 frames); median spread 0.249 m both; on the 1,546 common cells the difference is median 0.000 m, 27 cells over 0.05 m, largest 0.20 m | `python3 dem_from_contours.py <ground.csv> dem_intertidal_7day --max-hs 1.5 --max-day-offset 0.15 --last-days 7 --series-dir series --cell 2.0 --min-points 3 --max-spread 0.5` | the cron's `dem_intertidal_7day_dem.png` and `_info.json` |
| Foreshore slope on the real page, profiles A-D, over -0.5 to +0.5 m NAVD88 | [REAL] | tan(beta) 0.105, 0.089, 0.125, 0.104 (1:8 to 1:11) | same build | the page's profile legend |
| The page's REMINDER | [REAL] | main printed "the DEM reads LOW ... Not corrected here" although the contours carry the setup since 6 Oct; this PR prints "Elevations include the wave setup ... on 100% of the points" | same | the cron log after `DEM written` |
| Two camera-days 0.3 m apart with `--max-day-offset 0.15` | [SYNTH] test | both kept, DEM built; the pre-PR script crashed (ValueError) | `python3 test_survey_products.py` (`test_dem_two_days`) | - |
| Exit codes: 0 built, 3 grids but no page, 4 no points (or every camera-day rejected); the old page is removed first | code | the cron logs a distinct WARNING for 3 and for 4 (`waterline_timex_cron.sh`, section 6); the survey products report the reason from `logs/dem.log` | read the cron's branches | the cron log |

### 2.4 Week-to-week change map

**Claim.** The change map compares like with like, never draws a verdict from a handful of
cells, and tells "stable" from "not compared".

* The reference week is rebuilt from the **current** ground file with the same settings as
  the new 7-day DEM (`dem_change.py --rebuild`), so a change of processing (setup, filter,
  datum) does not show as beach change.
* A reference with fewer than 200 cells of >= 5 crossings, or less than 0.5 m of relief
  (5-95%), is degenerate and skipped, saying so; the window ending a day earlier is tried,
  back 7 days.
* A "uniform change" warning needs >= 50 compared cells; every change statement gives the
  number of cells.
* The new DEM is the series copy this run wrote (`--b`), not whatever is newest in the folder.

| what | data | number | reproduce here | check on the NUC |
|---|---|---|---|---|
| The station's own change map of 8 Oct 2026 (main), against the archived `dem_2026-09-27_7d` | [REAL] station output | "100% of overlapping cells moved the same way, median -0.75 m" on **7 cells**; the reference had 144 cells, all +0.58 to +0.64 m (storm days) | - | `archive/dems/change_2026-09-27_7d_to_2026-10-06_7d.png` |
| Like-for-like rebuild on the station's ground file, 17-22 Sep vs 1-6 Oct 2026 (8 Oct 2026) | [REAL] | 354 cells; median -0.12 m, p10 -0.50, p90 +0.32; level of detection 0.156 m; erosion 137 cells -235 m3, accretion 64 cells +89 m3; largest loss (to ~-1 m) at N 4638400-4638470 near the c1/c2 seam: real, not uniform | `dem_change.py --a <A> --b <B>` | `python3 dem_from_contours.py contour_points_ground.csv /tmp/wA --start-date 2026-09-17 --end-date 2026-09-22 --max-hs 1.5 --max-day-offset 0.15 --cell 2.0 --min-points 3 --max-spread 0.5 --no-plot`, the same for 2026-10-01..2026-10-06 into `/tmp/wB`, then `python3 dem_change.py --a /tmp/wA --b /tmp/wB --output-dir /tmp`; expect similar, not identical, numbers (the filter changed) |
| 29-30 Sep vs 1-5 Oct 2026 from the real contours, each version's own DEMs | [REAL] | main: 102 cells, median -0.046 m; this PR: 106 cells, median -0.044 m, LoD 0.097 m, erosion 30 cells -17.0 m3, accretion 8 cells +4.8 m3 | `dem_change.py --a wkA --b wkB` | - |
| The cron's own call on a ground file that starts 29 Sep | [REAL] | main: "Only DEMs less than 7 days apart so far"; this PR: "reference rebuild: 0 rows from 2026-09-15 to 2026-09-28 ... No change map this run: no adequate reference window" (no false map) | the cron's call, section 6b | the cron log after `window DEM written` |
| Whole cron, synthetic 24-day station | [SYNTH] cron harness (scratch) | 989 cells compared, median -0.021 m; a planted 2-day-old `dem_change_ref_*` folder removed | not in the repo | - |

On the NUC after the pull, expect the change map to name the cells it rests on in its title,
grey for "not compared" and light grey for "compared, within the level of detection". If the
archive before the window is thin, the log says "No change map this run: no adequate
reference window" instead of drawing one.

### 2.5 GNSS-R record guard

**Claim.** `archive/gnssr_record_start.txt` ("YYYY-MM-DD N") now keeps N as the **most** days
with readings any accepted spline has had. A spline that starts later than the record, or
has lost more than 3 days against that best, is refused and the previous products are kept.
Deleting the file accepts a shortening made on purpose.

| what | data | number | reproduce here | check on the NUC |
|---|---|---|---|---|
| 20 cases: header-only spline, missing spline, first run, slow loss (59 -> 56 accepted, -> 53 refused, -> 50 refused, back to 57 accepted), files written by the old version, shortening on purpose | [SYNTH] harness that runs the cron with every python3 call faked (scratch, not shipped) | 20 PASS, 0 FAIL on the final cron | not in the repo | `cat archive/gnssr_record_start.txt` before and after the first run: same start date, N equal or higher |
| QC spike test rewritten from O(n^2) to O(n log n) | [SYNTH] + [REAL] Aug-Sep 2026 spline with the Chatham gauge (commit `08a0cfc`) | identical flags, reasons and fit; a 4-year record ~4 s instead of ~17 s | - | `gnssr_qc_7day.png` looks as before |

The station's record started on 2026-07-09 when this was written (start file "2026-07-09 90").

### 2.6 Code health

| what | number | reproduce |
|---|---|---|
| Python 3.8 syntax of every `.py` changed since `af8e1bd` | vermin: minimum 3.8, no violations; `py_compile` passes on the station-floor venv | `pip install vermin && vermin -t=3.8- --no-tips $(git diff --name-only af8e1bd -- '*.py')` |
| Shell syntax | `bash -n` passes on both `.sh` files | `bash -n waterline_timex_cron.sh` |
| Self-tests on three Python setups (section "How to read this file") | all pass, except `view_reproject.py` on the two 3.8 venvs (no scipy; file unchanged) | section 1 |
| Whole cron on the station-floor venv, synthetic station | rc 0 in 381 s, every processing stage rc 0 (the network fetch stubs exit 1 and station_status exits 2, as expected offline) | scratch harness |

---

## 3. Survey-date products (`survey_products.py`)

One command per survey date builds the waterlines, a filtered copy, a DEM (+ GeoTIFFs), maps
and the comparison with that date's survey, with an independence label and its reason.
Output: `/mnt/I2Rgus_Data/survey_products/<date>/` (`README.txt`, `compare/`, `dem/`,
`maps/`, `waterlines/`, `provenance.json`) and, with `--summary`, `summary.txt/.csv/.png` at
the root. Full documentation: `surveys/README.md`.

**Labels.** INDEPENDENT: nothing in the chain was fitted to or placed with this survey.
CROSS-VALIDATED: fitted on other data from the same source, held out. PARTLY-CIRCULAR: some
step used this survey. CIRCULAR: pointing or offsets fitted to this survey. A comparison is
never labelled better than its weakest step, and the reason is printed next to the label.

Sign convention everywhere: **DEM (or waterline) minus survey; negative = the camera product
is LOW.**

### 3.1 2025-01-23: Jan 2025 lidar (c1 + c2, 18-23 Jan)

**Label: INDEPENDENT.** Calibration `CACO03_<cam>_20250123_EO` solved from GCPs surveyed that
day, not from the lidar; the detection search envelope is placed with the **March** lidar;
water level and waves from the ADCP; C fitted in 2026 to other data.

| what | data | number | reproduce |
|---|---|---|---|
| DEM - Jan lidar | [SYNTH] fixture photos rendered from the lidar with the ADCP water level and the setup planted | median -0.043 m, NMAD 0.062, n 526 cells (c1 -0.035 on 176, c2 -0.047 on 344) | fixture build (scratch fixture, not shipped) |
| Waterlines - Jan lidar | [SYNTH] | -0.036 m on 100 frames (c1 -0.039, c2 -0.033) | same |
| Same lines without setup (C = 0 sensitivity) | [SYNTH] | -0.204 m (c1 -0.215, c2 -0.193) | same |
| Fixture's expected answers (from its README) | [SYNTH] | detector with setup: lines c1 -0.035, c2 -0.028 m; DEM -0.034 m. The build lands within 0.01 m of them | - |
| Earlier check of the REAL photos (pointing_fix_figure.py, same GCP calibrations, NO setup, unfiltered) | [REAL], run on the NUC before this PR; itself PARTLY-CIRCULAR (its envelope came from this lidar) | printed "lidar - water level at the waterline" c1 +0.46 m (NMAD 0.19), c2 +0.15 m (NMAD 0.23), i.e. **waterline - lidar c1 -0.46, c2 -0.15 m: the lines read LOW** | held in `survey_products.py` PRIOR_CHECKS |
| Forcing for the window | [REAL] ADCP | setup median 0.291 m (daytime 0.281), range 0.147-0.708 m at C = 0.037 | `historical_forcing.py --start 2025-01-16 --end 2025-01-24 --output-dir <dir>` |

**[NUC] What to expect from the real photos.** With C = 0.037 each line keeps 0.6-0.8 of its
~0.28 m setup once re-projected, so, unless the new search envelope changes the detections,
waterline - lidar about **c1 -0.29 to -0.24 m, c2 +0.02 to +0.07 m**. The build's README
gives its own expectation from its paired C = 0 sensitivity. A ~0.3 m difference between
the cameras, measured over the same hours, points at the pointing or the lens model, not at
C. The ADCP-era still-water question (section 4.1) may move both by a few centimetres.

```
python3 survey_products.py --date 2025-01-23 --dry-run     # photos per day, files found, label, run time
python3 survey_products.py --date 2025-01-23               # 45 min here (120 photos, 20-25 s each); the NUC prints its own rate as it runs
```
Read: `/mnt/I2Rgus_Data/survey_products/2025-01-23/README.txt` (COMPARISON AND ITS LABEL:
`DEM - survey`, `waterlines - survey`, `by camera`, `without the wave setup (C = 0...)`,
and the expectation line), then `compare/2025-01-23_jan_lidar_comparison.png`.

### 3.2 2025-03-06: Mar 2025 lidar (c2 only, 3-9 Mar)

**Label: INDEPENDENT.** Calibration `CACO04_c2_20250219_EO` (GCPs), envelope from the
**January** lidar, ADCP water level and waves. c1 has no photos on the NUC for this period.

| what | data | number |
|---|---|---|
| DEM - Mar lidar | [SYNTH] fixture | median -0.034 m, NMAD 0.057, n 1,263 cells |
| Waterlines - Mar lidar | [SYNTH] | -0.026 m on 63 frames |
| C = 0 sensitivity | [SYNTH] | -0.243 m |
| Fixture's expected answers | [SYNTH] | detector with setup: lines -0.026 m; DEM -0.034 m: matched |

**[NUC]** No earlier real number exists for March. `python3 survey_products.py --date 2025-03-06`
(26 min here, 70 photos). Read `2025-03-06/README.txt` and
`compare/2025-03-06_mar_lidar_comparison.png`.

On the static synthetic beach the per-day medians drift by ~0.09 m (Jan) and ~0.17 m (Mar)
over the week, and two days at the same elevation differ by ~0.13 m: a day-to-day difference
below ~0.15-0.2 m on the real dates is not by itself evidence of beach change.

### 3.3 2026-09-29: RTK check shots of 29 Sep 2026 (live station, window 26 Sep - 2 Oct)

**Label: CIRCULAR with the station's C = 0.037**, because that C was first fitted to these
very shots (`survey_products.py` SETUP_FITS; recognised by content, name, coordinates or
date). The comparison then cannot test the setup or the overall level, only the spread over
the elevations the survey covers. **A build with C = 0 is INDEPENDENT** (nothing else in its
chain used the shots).

| what | data | number |
|---|---|---|
| Waterlines - RTK at C = 0.037 (CIRCULAR) | [REAL] live rows 29 Sep - 2 Oct (the sandbox copy starts 29 Sep) and the RTK | median -0.127 m, NMAD 0.098, 16 frames (c1 -0.164 on 10, c2 -0.121 on 6) |
| Where the survey stops | [REAL] | the transects end at their lowest shot, +1.17 to +1.80 m; 6 of 16 frames lie below it, where a line can only read low. Over only the lines >= 0.10 m above it: -0.069 m on 7 frames (too few to be an estimate) |
| DEM - RTK | [REAL] | n 2 cells: too few (the DEM hardly reaches the transects) |
| C = 0 build (INDEPENDENT) | [REAL] | waterlines - RTK -0.330 m on 12 frames (c1 -0.389 on 7, c2 -0.272 on 5); **all 12 lie below the lowest shot, so this is biased LOW by construction**; DEM n 3 (too few) |
| Same frames both ways (paired) | [REAL] | -0.330 m without the setup, -0.108 m with it (12 frames) |
| Repo tool on the whole 29 Sep - 5 Oct week (`compare_rtk.py`), this PR's filter, C = 0.037 | [REAL] | RTK - waterline +0.140 m over 23 frames (c1 +0.138 on 16, c2 +0.140 on 7); trend with water level: offset = +0.578 - 0.466 x level (the signature of the survey floor) |
| Same tool on the unfiltered still-water rows (C = 0) | [REAL] | +0.333 m over 18 frames (c1 +0.369 on 12, c2 +0.291 on 6) |

In the sandbox this date is `partial` (exit 3) only because there is no photo of the window
there (no pointing check, no photo maps); with two real 7 Oct photos standing in it builds
`complete`. Every output of the final builds was made at commit `3c03f0e` (README `software`
line).

**[NUC]** The station has the 26-28 Sep rows (the storm's last day included) and the photos,
so the numbers there will differ from the sandbox's; the product should be `complete`.

```
python3 survey_products.py --date 2026-09-29        # ~1-2 min (52-79 s here): no detection, reads the station's contour file
python3 survey_products.py --date 2026-09-29 --setup-coef 0 \
    --output-root /mnt/I2Rgus_Data/survey_products_c0   # the INDEPENDENT C = 0 build, in its own root
python3 survey_products.py --summary --output-root /mnt/I2Rgus_Data/survey_products_c0
```
Read: `2026-09-29/README.txt`, the lines `waterlines - survey`, `WHERE THE SURVEY STOPS` and
`without the wave setup`, and `compare/2026-09-29_rtk_2026-09-29_comparison.png`.

**No unbiased real-data level number exists for the 2026 lines yet.** The RTK walks the
upper beach only; an RTK survey reaching the low-tide line, on days with different waves,
would give one (section 6).

### 3.4 Disabled dates and safety behaviour

| what | data | number | reproduce |
|---|---|---|---|
| 2024-10-23 (GCP file not on the NUC yet) and 2025-03-19 (no photos, no ADCP) refused | config | exit 2: "Not built. To build it anyway: --force-disabled" | `python3 survey_products.py --date 2024-10-23` |
| Another C in a root that holds a build | [SYNTH] tests + a real refresh over a C = 0 root | refused, exit 5, nothing changed; the message names the root that C belongs in (`<root>_c<C>`); `--replace` overwrites on purpose | `test_survey_products.py` |
| A `--setup-coef` of unknown origin | tests | exit 5, nothing written, unless `--setup-fitted-to` says where it was fitted | same |
| Two builds of one date at once | tests | the second exits 4 and names the build holding the lock | same |
| `--all` exit code | tests | the worst by severity: failed (1) > refused (5) > busy (4) > partial (3) > 0 | same |
| Rerun on unchanged inputs | [SYNTH] + [REAL] | every step skipped, numbers identical (every date rerun several times on the final code) | rerun the same command |

---

## 4. OFFSETS AND CORRECTIONS

Every offset or correction in the chain from photo to DEM, with its equation, value, reason,
source, how well it is verified, its limits, and how to re-check it on the NUC.

### 4.1 Wave setup on the waterline elevation (the setup coefficient C)

**Equation used.** Each waterline is given

    z_line = eta_still + C * sqrt(Hs * L0),     L0 = g * Tp^2 / (2 * pi)

where `eta_still` is the still-water level (GNSS-R on the live station, the ADCP for the 2025
dates), `Hs` the offshore significant wave height and `Tp` the peak period of the frame
(`extract_elevation_contours.py --setup-coef`; `survey_products.apply_setup` uses the same
formula and rounding). **Value: C = 0.037** (`waterline_timex_cron.sh` SETUP_COEF; the same
value in `survey_products.py` and `historical_forcing.py`, checked by a test). Not changed by
this PR.

**Reference: Stockdon et al. (2006).** From 10 field experiments, Stockdon et al. parameterise
the wave-driven water levels at the shoreline as

* setup:            `<eta> = 0.35 * beta_f * sqrt(H0 * L0)`
* incident swash:   `S_inc = 0.75 * beta_f * sqrt(H0 * L0)`
* infragravity swash: `S_IG = 0.06 * sqrt(H0 * L0)`
* total swash:      `S = sqrt(S_inc^2 + S_IG^2) = sqrt(H0 * L0 * (0.563 * beta_f^2 + 0.004))`
* 2% exceedance runup: `R2% = 1.1 * (<eta> + S / 2)`

with `H0` the deep-water significant wave height, `L0 = g T0^2 / (2 pi)` the deep-water
wavelength from the peak period, and `beta_f` the foreshore beach slope. (For very
dissipative beaches, Iribarren number below 0.3, they give `R2% = 0.043 sqrt(H0 L0)`
instead; Marconi's foreshore is not in that range: at beta_f 0.1, Hs 1 m, Tp 8 s the
Iribarren number is ~1.0.) The repository implements the same equations in
`foreshore_slope.py` and `gnssr_qc.py`.

**What our C means.** Our form is Stockdon's setup with `C = 0.35 * beta_f`. C = 0.037 implies
**beta_f = 0.037 / 0.35 = 0.106 (about 1:9.5)**. The real 7-day DEM of 29 Sep - 5 Oct 2026
measures foreshore slopes of 0.089-0.125 (1:8 to 1:11) over -0.5 to +0.5 m NAVD88 at four
positions (section 2.3), i.e. Stockdon's setup coefficient 0.031-0.044 there. Worked example,
Hs 1 m, Tp 8 s (L0 = 99.9 m, sqrt(H0 L0) = 10.0 m), beta_f 0.106: setup 0.37 m (= C x 10.0 m,
the cron comment's "+0.37 m"), S = 1.02 m, R2% = 0.97 m. A line marked at R2% would need an
effective C of ~0.10; the fits below are near 0.03-0.05, so the lines sit near the mean setup
level, not near the runup maximum.

**Why it is needed.** A timex (10-minute time-exposure) image does not show the still-water
line. It shows the time-averaged edge of the water, which sits above still water by the wave
setup and part of the swash. Video shoreline-elevation models therefore add wave terms to the
tide (Plant and Holman, 1997; Aarninkhof et al., 2003), and detection methods differ in how
far above still water their line lies (Plant et al., 2007). Without a setup term every line
reads low, by more in rough water. **So our C is an empirical coefficient for this detector,
this camera geometry and this wave record, not pure Stockdon setup.** It can absorb part of
the swash, any detector offset that grows with the waves, and it is reduced by whatever setup
the still-water reference already contains (below). It was fitted with `Hs` in ADCP currency
at 21 m depth (`hs_best` of `archive/waves_marconi.csv`) and NDBC 44008 peak periods, not
Stockdon's deshoaled deep-water `H0`.

**Re-projection: why the DEM moves less than the setup.** The setup changes the elevation
plane a pixel is projected onto. A higher plane moves the georectified point landward along
the camera ray, onto higher beach. So a line gains only part of its setup in the DEM: the
fixture gives the error of a no-setup line as `-setup * (1 - beta / tan(ray depression))`,
about -0.65 x setup there. Measured shares kept: ~0.71 per frame on the 2026 RTK frames
[REAL], 0.6 (c1 0.63, c2 0.57) and 0.8 in the Jan and Mar 2025 fixture builds [SYNTH],
0.64-0.75 for most frames in the refit [SCRATCH]. **DEM sensitivity: 0.057 m per 0.01 in C**
on the real 29 Sep - 5 Oct DEM [SCRATCH]; DEM(C = 0) is 0.22 m below DEM(0.037).

**Where 0.037 came from, and every fit made since.**

| fit | data | C | 90% interval | status |
|---|---|---|---|---|
| First fit, 8 Oct 2026: per frame `C = -(line - RTK)/sqrt(Hs L0)`, each line at its still-water position (no re-projection) | [REAL] 14 frames, 29 Sep - 5 Oct 2026, 2026-09-29 RTK | **0.037** (median) | 0.028-0.058 | the value in use; every crossing used, so the survey-floor bias below applies |
| `dem_from_contours.py --fit-setup` (spread of repeat crossings per 2 m cell, no survey, no re-projection), unfiltered still-water rows, Hs <= 1.5 m | [REAL] 193 frames | best 0.030; near-best band +0.026 to +0.040; spread 0.309 m at C = 0, 0.283 at best, 0.285 at 0.037 | (band, not an interval) | prints "Not distinguishable from the C in use (0.037) ... keep it" (re-run for this file) |
| same, filtered C = 0 lines, Hs <= 1.5 m | [REAL] | best 0.026; band +0.006 to +0.034, FLAT; 0.037 spreads 6 mm more | - | "do not apply on this window alone" |
| same, C = 0 lines of 29 Sep - 2 Oct only, filtered, Hs <= 1.5 m | [REAL] | best -0.004; band -0.020 to +0.010, FLAT | - | "not distinguishable from C = 0": one week can sit near 0 |
| Same spread metric with each line re-projected | [SCRATCH] | 0.029 | frame 0.000-0.040; day-block 0.000-0.086 | a flat minimum (spread 0.255 at 0, 0.249 at 0.028, 0.256 at 0.036) |
| Within-cell slope of elevation on sqrt(Hs L0), re-projected (no survey) | [SCRATCH] | **0.043** | frame 0.026-0.053; day-block 0.017-0.069 | the best-behaved internal estimate |
| RTK, re-projected, all crossings | [SCRATCH] | 0.067 | frame 0.056-0.076; day-block 0.045-0.087 | **biased high by where the survey stops** (below) |
| RTK, re-projected, only lines >= 0.1 / 0.2 m above each transect's lowest shot | [SCRATCH] | **0.046 / 0.043** | day-block ~0.02-0.055 | waterline - RTK at 0.037 there: -0.037 / -0.028 m on 10 / 7 frames (too few for an estimate) |
| RTK, two parameters (`z = tide + a + C phi`), all crossings | [SCRATCH] | 0.049, a = +0.11 m | 0.017-0.085 | also affected by the survey floor |

**The survey-floor bias.** The 29 Sep 2026 transects end at their lowest shot, +1.17 to
+1.80 m NAVD88. A line enters the comparison only if it lands on a surveyed stretch. A line
whose plane lies below the lowest shot can only read low there (the RTK under it is higher),
and a line that would read high at that water level lands seaward of the survey and is never
compared. So every all-crossings RTK number is biased low, and a one-parameter fit pushes C
up to compensate. At C = 0.036, 9 of 23 frames (38% of the points) lie entirely below their
transect's lowest shot; their median is -0.24 m, the other 14 frames -0.09 m [SCRATCH]. A
simulation through the real camera model, lines planted on the RTK profiles with a true C of
0.037 [SCRATCH]: with Gaussian frame errors of 0.15-0.25 m the all-crossings fit returned
0.041-0.045; with 20-30% of frames reading 0.3-0.8 m low (as some real crossings do) it
returned 0.053-0.058, with a median residual at 0.037 of -0.11 to -0.15 m (real: -0.14 m),
while the above-floor fit returned 0.042-0.044. A true constant +0.12 m level offset, by
contrast, is not removed by the above-floor cut (0.055 -> 0.051 in the simulation); the real
data fall from 0.067 to 0.043-0.046, which matches the survey floor, not a level offset. The
same bias affects the old first-fit numbers and the survey-product numbers "0.33 m low without
setup, 0.13 m low with it": they are not a measured level error.

**Decision: C stays 0.037.** It lies inside the 90% intervals of every estimate the survey
floor does not bias (~0.00-0.04, 0.026-0.053, ~0.02-0.055) and inside the unfiltered rows'
`--fit-setup` band (the filtered lines' band stops just below it). Moving to ~0.043-0.046
would raise the DEM by ~0.03-0.05 m. **The RTK cannot pin C more closely.** An RTK survey
reaching down to the low-tide line, on two or more days with different waves, would.

**Still-water reference: ADCP (2025) vs GNSS-R (2026).** C was fitted with each line at the
GNSS-R level. The GNSS-R footprint is the surf zone (reflections 70-210 m from the antenna,
waterline 55-90 m: `gnssir_reflection_audit.py`), where breaking waves raise the mean level,
so GNSS-R contains a share s of the setup and C carries only the rest. The 2025 dates use the
ADCP at 21 m depth (outside the surf zone), which sees none. The mean part of the GNSS-R share
is already in the datum chain (GNSS-R sits ~0.02 m below Chatham after its datum fix, and the
ADCP datum is tied to Chatham's mean), so `historical_forcing.py` expects each 2025 frame at
about **+0.02 m - s x (its setup - the mean setup of the GNSS-R/Chatham comparison period)**:
~0.02 m high in average waves, low only in above-average waves. s is not measurable in the
sandbox (the ADCP and GNSS-R records never overlap). [NUC] `historical_forcing.py` measures s
from the 2026 record (GNSS-R spline, `archive/gauge_8447435.csv`, `archive/waves_marconi.csv`)
and writes it into `forcing/forcing_report.txt` of each 2025 product. (The refit's coarser
bound [SCRATCH], which ignores the mean part, put the 2025 DEMs up to 0.03-0.10 m low: an
upper estimate, not a measurement.)

**Wave currency.** The 2025 dates use the ADCP's own Hs (C's currency) and its peak period;
`historical_forcing.py` measured ADCP Tp / WIS Tp = 1.061 (r 0.67, 534 h) [REAL], so setup
from the ADCP period sits +0.020 m (median) above setup from an open-ocean period. Not
corrected; stated in each forcing report.

**Known limits.** C differs with beach slope (the internal slope estimator gives 0.017 for
frames below -0.2 m, 0.045 at -0.2 to 0.6 m, 0.035 above 0.6 m [SCRATCH]), and the slope
changes with storms. The 2025 beach faces differed from the Sep-Oct 2026 one (tan(beta)
0.089-0.125). The profiles of the fixture DEMs, whose shape comes from the real 2025 lidar
(the fixture is rendered from it), give over -0.5 to +0.5 m:
* 18-23 Jan 2025: 0.125-0.135 (1:7-1:8), where Stockdon's form would scale C to ~0.044-0.047;
* 3-9 Mar 2025: 0.067-0.075 (1:13-1:15), i.e. C ~0.024-0.026.

If the real C follows the slope, 0.037 leaves the January lines ~0.03-0.04 m low and lifts the
March lines ~0.05-0.06 m too high (difference in C x the window's median sqrt(Hs L0), 7.8 m in
January and 6.2 m in March, x the share kept). These are estimates, not measurements. They add
to the ADCP still-water effect above. The real builds on the NUC are the test.

The re-projecting fits exist only as scratch scripts; `dem_from_contours.py --fit-setup` does
not re-project and is in principle biased low by roughly the re-projection factor.

**Re-check on the NUC.**
```
python3 survey_products.py --date 2026-09-29 --setup-coef 0 --output-root /mnt/I2Rgus_Data/survey_products_c0
python3 dem_from_contours.py /mnt/I2Rgus_Data/survey_products_c0/2026-09-29/waterlines/contour_points_ground_filtered.csv \
    /tmp/fit_2026-09-29 --fit-setup --max-hs 1.5 --no-plot
```
Read the printed curve, the near-best band and the verdict. It offers a C to apply only when
C = 0 and the C in use both lie outside the band and the minimum is not flat. Adopt a new C
only if it holds over several weeks of different waves; then change it in all three places
(the cron, `survey_products.py`, `historical_forcing.py`) and give it a SETUP_FITS entry.

### 4.2 GNSS-R water level onto NAVD88: +0.349 m

* **Value.** `GNSSR_ANTENNA_NAVD88_M = 19.014` in `extract_elevation_contours.py`; the
  correction applied is 19.014 m minus the height in the spline header (18.665 m), +0.349 m.
* **Reason.** The gnssrefl station height came from CSRS-PPP, which reports CGVD2013 (Canadian
  datum), not NAVD88. NGS OPUS on the same RINEX (2026-09-28; NAD83(2011), GEOID18) gives
  19.014 m NAVD88. The two solutions agree on the ellipsoid height to 2 mm (-10.025 vs
  -10.023 m), so only the height system differed.
* **Source.** Commit `4932100` (this repository) and `fb10545` (gps_code).
* **Verified.** The uncorrected levels were low by an amount consistent with the Chatham
  gauge (-0.37 m including a harbour effect) and a global tide model (0.26 m)
  (`extract_elevation_contours.py` comment). After the fix GNSS-R sits ~0.02 m below Chatham
  (`historical_forcing.py`). Against the EOT20 tide model (Hart-Davis et al., 2021), moved
  +0.090 m onto NAVD88, a 9 Oct 2026 run reported mean -0.047 m, de-meaned RMS 0.111 m,
  r 0.989 (`compare_gnssr_to_tidemodel.py`; not re-run for this PR).
* **Limits.** +/-0.061 m, mostly the geoid model. Every live elevation shares it; it largely
  cancels against RTK shots on GEOID18. A +/-0.06 m datum error is about +/-0.011 in an RTK
  fit of C [SCRATCH].
* **Re-check.** `head -5 /home/argus_user/GNSS/v4.1/products/refl_code/Files/usgs/usgs_spline_out.txt`
  (the header height); if gnssrefl is ever set to 19.014 m the correction falls to 0 by itself.

### 4.3 ADCP water level onto NAVD88: +0.1012 m

* **Value.** `adcp_water_level_navd88.csv` = ADCP `water_level` + 0.1012 m (`adcp_to_navd88.py`).
* **Reason.** The Signature 1000 ADCP measures pressure depth, not a datum height.
* **Source.** The mean Chatham (8447435, NAVD88) minus the mean ADCP over the same hours.
  This **assumes Marconi's mean level equals Chatham's** over the deployment (Dec 2024 - Mar
  2025). The tide-model MSL-to-NAVD88 constant (+0.09 m) would give 0.1127 m.
* **Verified.** `historical_forcing.py` recomputes it (0.101 m) and reports both [REAL].
* **Limits.** A real mean-level difference between the open coast and the harbour shifts
  every ADCP level, and every 2025 DEM, by the same amount; it is not in `sigma_m`.
* **Re-check.** `forcing/forcing_report.txt` of each 2025 product, line `ADCP datum`.

### 4.4 Chatham to Marconi water-level transfer (only for dates outside the ADCP)

* **Used for.** 2024-10-23 (disabled now) and hours the ADCP lacks. The 2025-01-23 and
  2025-03-06 windows are inside the ADCP record: water level measured.
* **Value** [REAL], fitted on the ADCP overlap: harmonic transfer, 13 constituents, M2 ratio
  1.481, Marconi leads by 70 min, residual x 0.85; leave-one-week-out RMS 0.076 m,
  extrapolation RMS 0.085 m (used as `sigma_m`). The linear alternative
  (Marconi = 1.403 x Chatham(t - 72 min) - 0.033 m) scored 0.166 / 0.174 m and is not used.
* **Live station** (`marconi_water_level.py`, for the OWG products, not the waterlines):
  refitted to GNSS-R over the last 30 days on every run; default (1.24, -48 min, -0.10 m,
  sigma 0.13 m) only without GPS overlap.
* **Limits.** The transfer changes with the inlet and harbour; months from the overlap the
  error can exceed the extrapolation RMS (measured 0-9 weeks out).
* **Re-check.** `forcing_report.txt`, block `Chatham -> Marconi transfer`.

### 4.5 Wave height and period conversions (dates or hours without the ADCP)

* WIS ST63064 Hs -> ADCP Hs: linear form `0.833 x Hs + 0.004 m` (r 0.67, CV RMS 0.456 m) or
  the direction-aware form (r 0.84, CV RMS 0.330 m), which wins; WIS Tp scaled by 1.061.
  NDBC 44013: `0.827 x Hs + 0.151 m` or direction-aware (r 0.87). Combined setup CV RMS
  0.079 m [REAL], 534 h of December overlap for WIS. Not used by the two enabled 2025 dates
  (ADCP waves throughout). Source: `historical_forcing.py`, `forcing_report.txt`.

### 4.6 Camera pointing (calibrations) and the sea-horizon check

* **Calibrations per period** (`calibration/calibration_history.csv`, `calibration/README.md`):
  on 24 Jan 2025 the cameras were re-set: c1 moved 5.6 m and turned -15.3 deg, c2 moved 5.8 m
  and turned -15.5 deg (CACO03 20250123 -> CACO04 20250219). Between 23 Oct 2024 and 23 Jan 2025
  the CACO03 cameras turned ~13 deg. The live CACO05 20251113 calibration: c1 has the
  2025-11-04 values, c2 moved 3.1 m and turned -3.1 deg.
* **Why the source of the pointing matters.** A pointing fitted to the lidar
  (`fit_eo_to_survey.py`) keeps the camera position fixed, so it absorbed the 24 Jan move as
  ~+22.7 deg of pan for c1 (+21.9 deg for c2) [REAL]. With the February 2025 pointing applied
  to January (what had been used before), the Jan 18-23 lines read "lidar - water level"
  +1.78 m (c1, spread 0.39 m) and +3.31 m (c2, spread 0.76 m); with the GCP calibration of
  23 Jan, +0.46 / +0.15 m (spread 0.19 / 0.23 m) [REAL, station runs before this PR, no
  setup]. The products therefore use the station's GCP calibration of each period and never a
  lidar-fitted pointing (a survey-fitted or carried pointing makes the comparison CIRCULAR,
  by rule).
* **Sea-horizon check.** Each build compares each camera's sea horizon with where its
  calibration puts it (`estimate_eo_rotation.horizon_rows`, full Earth-curvature dip with
  refraction k = 0.13; until Oct 2026 it held half the dip and read +0.06-0.07 deg on a perfect
  camera). Only day-to-day changes leave days out; a constant offset is reported with its DEM
  sensitivity.
  * [REAL] 2026 cameras (the 7 Oct 2026 photos): c1 41 px (+0.90 deg tilt / +0.23 deg roll),
    c2 34 px (+0.57 / -0.69 deg) against CACO05 20251113. The calibration puts the horizon at
    the top edge of the frame, partly off it, where the lens model is extrapolated far beyond
    its GCPs, so this may be the lens model rather than the pointing.
  * Is it pointing? Judged only on what fitted nothing to the RTK (C = 0 lines minus what a
    repeat-crossing C of 0.029-0.043 explains), up to about 29-49% of c1's and 12-23% of c2's
    full offset cannot be ruled out [REAL + SCRATCH]; that share is **overstated**, because
    every C = 0 line lies below the survey's lowest shot.
  * [NUC] The 2025 horizon offsets can only be measured on the real 2025 photos. The sandbox
    fixture's horizon is an artefact of the rendering (c1 +0.17, c2 +0.14 deg in January): the
    fixture build reports it like any other offset, under its SYNTHETIC FIXTURE banner, and it
    says nothing about the real cameras.
  * Sensitivity: 0.1 deg of tilt moves a line ~6-9 m at 250-350 m range, ~0.13-0.17 m of DEM.
* **Re-check.** Each product README, `sea horizon vs the calibration` and the KNOWN CAVEATS
  pointing block; `pointing/pointing_<cam>.csv`.

### 4.7 Detector edge position

* The detector places its line where a water pixel has foam around it. In the fixture this is
  ~28 rows seaward of the foam/wet-sand edge, and the fixture was drawn to match; on the
  rendered beach the detector then reads -0.04 to -0.07 m on the low beach and +0.01 to +0.05 m
  on the upper beach [SYNTH]. **Its real-world bias is unknown**: it is part of what C absorbs
  and part of the c1/c2 difference the real January build will show.

### 4.8 Survey datums and coordinates

* RTK 29 Sep 2026: NAD83(2011) / UTM 19N + NAVD88 (GEOID18), stated in the file's `CS name`;
  the file's own ellipsoid-minus-elevation is -27.80 m (e.g. RM2: -14.725 - 13.076) [REAL].
* Lidar 2025: no vertical datum in the files; NAVD88 assumed (GEOID12B and GEOID18 differ by a
  few cm here).
* Camera GCP frame: no CRS in the GCP files; assumed NAD83(2011) / UTM 19N (EPSG:6348). Were
  they WGS 84 / ITRF, the grids would sit ~1-1.5 m off the surveys, up to ~0.1-0.15 m of
  elevation across a 1:10 foreshore. GeoTIFFs are tagged EPSG:32619 as asked (a nominal tag;
  `--geotiff-epsg 6348` writes the true code).

### 4.9 Thresholds that move the numbers (not offsets, but each changes a result)

| setting | value | where | why |
|---|---|---|---|
| Max offshore Hs for the DEM | 1.5 m | cron `DEM_MAX_HS` | the setup errs most in big waves |
| Day-offset test | 0.15 m; not run on fewer than 3 camera-days | cron `DEM_MAX_DAY_OFFSET` | catches a camera-day displaced as a whole |
| DEM cell / frames / spread | 2 m / >= 3 / <= 0.5 m | cron | repeatability |
| Change: LoD floor / crossings | 0.05 m / >= 5 in both | `dem_change.py` | errors a single week cannot see |
| Change: degenerate reference | < 200 cells or < 0.5 m relief | `dem_change.py` | the 8 Oct false "uniform" map |
| Change: uniform warning | >= 50 cells | `dem_change.py` | no verdict on a handful of cells |
| Filter: threshold / window / bins | max(0.75 m, 4 x robust scatter) / +/-3 d, exp(-dt/1 d) / 32 px | `waterline_consistency.py` | section 2.1 |
| Survey floor margin | 0.10 m above each transect's lowest shot | `survey_compare.py` SURVEY_FLOOR_MARGIN | section 4.1 |

---

## 5. Reproducing the before/after on the real contours

The before/after figures of the PR were made from the real rows of 29 Sep - 5 Oct 2026, with
the setup added as the cron adds it, with main (`af8e1bd`, in a git worktree) and with this
branch, the same commands for both:

```
python3 waterline_consistency.py contours.csv --output contour_points_timex_qc.csv \
    --report waterline_consistency_report.csv --plot waterline_consistency --image-dir <photos> --plot-days 7
python3 daily_elevation_map.py contour_points_timex_qc.csv <photos> c1 elevation_map_c1_7day.png --days 7   # and c2
python3 georectify.py contour_points_timex_qc.csv contour_points_ground.csv \
    --io-c1 calibration/CACO05_c1_20240801_IO.yaml --eo-c1 calibration/CACO05_c1_20251113_EO-CV.yaml \
    --io-c2 calibration/CACO05_c2_20240801_IO.yaml --eo-c2 calibration/CACO05_c2_20251113_EO-CV.yaml
python3 dem_from_contours.py contour_points_ground.csv dem_intertidal_7day --max-hs 1.5 --max-day-offset 0.15 \
    --last-days 7 --series-dir series --cell 2.0 --min-points 3 --max-spread 0.5
python3 dem_from_contours.py contour_points_ground.csv wkA --max-hs 1.5 --max-day-offset 0.15 \
    --start-date 2026-09-29 --end-date 2026-09-30 --cell 2.0 --min-points 3 --max-spread 0.5 --no-plot   # wkB: 10-01..10-05
python3 dem_change.py --a wkA --b wkB --output-dir .
```

main's `daily_elevation_map.py` needs matplotlib 3.5 or newer (`matplotlib.colormaps`); this
branch's falls back to `plt.get_cmap` on 3.3/3.4. On the NUC the simplest before/after is the
copy of today's pictures made before the pull (section 1).

---

## 6. What is not verified yet, and what would verify it

* **The real 2025 numbers.** Only the NUC has the real January and March 2025 photos. The
  sandbox numbers for those dates are synthetic and test the pipeline only.
* **An unbiased 2026 level.** The 2026 RTK stops at +1.2 to +1.8 m; every C = 0 line, and some
  C = 0.037 lines, lie below it. A calm, spring-low-tide RTK survey reaching the low-tide line,
  on two or more days with different waves, would give an INDEPENDENT level check of the whole
  intertidal and pin C and any level offset separately.
* **The ADCP-era still-water share s.** Measured by `historical_forcing.py` on the NUC.
* **The 2026 pointing.** A sea-horizon (or GCP) check of the CACO05 20251113 calibration would
  tell lens model from pointing.
* **The detector's rate on the NUC10i3.** Measured here at 20-25 s per photo on a 2.8 GHz
  Xeon; the NUC prints its own rate.
* **Reviews.** The waterline work had three review rounds (the last: all lenses "ship"), the
  survey products three verification and three fix rounds, the combined branch two review
  rounds (interactions; operator and documents) and two fix rounds, and the setup refit an
  adversarial check that found the survey-floor bias. The last fix round (`3c03f0e`) was
  checked by its own tests and the full check list, not by a further independent review.

## References

* Aarninkhof, S.G.J., Turner, I.L., Dronkers, T.D.T., Caljouw, M., Nipius, L., 2003. A
  video-based technique for mapping intertidal beach bathymetry. Coastal Engineering 49(4),
  275-289.
* Hart-Davis, M.G., Piccioni, G., Dettmering, D., Schwatke, C., Passaro, M., Seitz, F., 2021.
  EOT20: a global ocean tide model from multi-mission satellite altimetry. Earth System Science
  Data 13, 3869-3884.
* Plant, N.G., Holman, R.A., 1997. Intertidal beach profile estimation using video images.
  Marine Geology 140(1-2), 1-24.
* Plant, N.G., Aarninkhof, S.G.J., Turner, I.L., Kingston, K.S., 2007. The performance of
  shoreline detection models applied to video imagery. Journal of Coastal Research 23(3),
  658-670.
* Stockdon, H.F., Holman, R.A., Howd, P.A., Sallenger, A.H., 2006. Empirical parameterization
  of setup, swash, and runup. Coastal Engineering 53(7), 573-588.
  https://doi.org/10.1016/j.coastaleng.2005.12.005

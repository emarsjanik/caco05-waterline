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
python3 test_survey_products.py                  # ends "all survey_products tests passed" (see below)
python3 test_survey_compare.py                   # ends "all survey_compare tests passed"
python3 test_historical_forcing.py               # ends "75/75 checks passed"
cat archive/gnssr_record_start.txt               # note "YYYY-MM-DD N"; see section 2.5
```

What a survey-date build reruns after this pull: a step reruns when a script in its signature
changes. This PR changes none of the forcing, pointing, detection, filter, DEM or photo-map
scripts; it changes `survey_compare.py`, so the next build of a date already built skips
the forcing, pointing, detection, filter, DEM and photo maps, and reruns only its plan-view
map and its comparisons (a minute or two per date). Checked on a 2026-09-29 build made with
PR #87's code, then rebuilt after copying in the files this PR changes: every step up to the
DEM skipped, the plan-view map and the comparison rebuilt (12 s), a further run skipped all.

`view_reproject.py --self-test` needs scipy, which does not import on the NUC. That file
is unchanged by this PR, and the survey products use only its `ground_to_pixel`, which
does not need scipy.

The self-tests give the same result on the NUC as in a clean checkout: they run on fixtures
and folders of their own, with the station's search folders and archive files switched off
(`test_historical_forcing.py`'s last check confirms that no build named a station file). The
three `test_*.py` files run in a temporary folder that is also their current folder; the three
`--self-test` runs write only into temporary folders of their own and look for nothing in the
current folder but matplotlib's own `matplotlibrc` (checked with strace). What they read of
the station is said.

`test_survey_products.py` checks `surveys/survey_dates.csv` and `surveys/surveys.csv` as the
owner has them now last of all, so a table problem cannot hide the results of the code's
checks, and lists every problem. It asserts only what `survey_products.py` itself needs of
them or does with them: both tables read with no warning (a comma in a text outside double
quotes cuts the text); every date switched off refused, and planned with `--force-disabled`;
every date switched on planned; the calibration of every camera switched on found in
`calibration/`; and the product's own label rules applied to their rows (a survey the search
envelope was placed with, or a points survey made on its calibration's day unless
`label_override_reason` is written, is no better than PARTLY-CIRCULAR). An edit written as
`surveys/README.md` says (a date added, switched on or off, a survey file given, a second
survey of a date, a row without its optional last column, a text with a comma inside double
quotes) never fails it. A comma outside double quotes does, as the README says: the product
then warns and cuts the text; a row that only ends in a comma loses nothing and is not
warned about. A table the product cannot read at all (e.g. saved by a spreadsheet in Latin-1,
not UTF-8) is one of the problems listed, after the code's checks. Every
`label_override_reason` the tables hold is listed (the text after such a comma, in a row
without its last column, would land there unwarned). With the tables, and whether or not they
pass, it checks the live cron as it is on the computer: `SETUP_COEF` in
`waterline_timex_cron.sh`, read as `survey_products.py` reads it (empty is C = 0), must be the
products' C, or the live waterlines and the products carry different setups (the cron's own
comments say `SETUP_COEF=""` for a refit: the test then FAILs, naming the cron's value, as
`survey_products.py` warns).

What was asked for in the tables (23 Oct 2024 and 19 Mar 2025 switched off, the windows, the
29 Sep 2026 RTK row and its text, envelopes from the other lidar, every calibration in the
repository, the Oct GCPs PARTLY-CIRCULAR, every row with all its columns) is checked on the
tables as the repository shipped them with the test: read with git as of the last commit of
`test_survey_products.py`, not of HEAD, so the owner's edits never reach it, even when they
are committed on the station to pull past them. A commit that changes the shipped tables
therefore changes the test's expectations with them; when `surveys/*.csv` at HEAD differ from
the shipped ones, the test says so (`(surveys/... changed in N later commit(s) ...)`): those
later tables are then checked only for what the product needs (above).

The checks of a pointing fit that is not on the computer (a correction carried from a lidar
fit that cannot be found must be CIRCULAR as a rule) name that fit after the test's own
temporary folder: the products also look for such a fit in the repository's `calibration/`,
where the station keeps its own fits, and no station file can have that name.

Three kinds of check depend on the computer, and are skipped, counted and explained when they
cannot run. Two checks of `test_survey_products.py` (the real RTK file renamed, and re-exported
under another date, must still be CIRCULAR) need the real `2026-09-29_Marconi_Checkshots.csv`
(the same rules are also checked, everywhere, on a test fit of the test's own). The test looks
for it as the survey products do (the path `surveys/surveys.csv` gives: one with a folder in it
as given, a bare name in the station's survey folders: `Chelsea_calibration`, `waterline`,
`waterline/calibration`, the repo's `calibration` and `surveys`; `--survey-dirs <folder>` adds
one) and confirms it by its sha256. If it is not
found, the two checks print `skip` and the last line reads `all survey_products tests passed (2
skipped: renamed-RTK check, re-exported-RTK check: the real RTK file
2026-09-29_Marconi_Checkshots.csv is not in ...)`: put the file in one of those folders, or
pass `--survey-dirs`, to run them. A file of that name whose sha256 differs (e.g. re-saved with
other line endings) is named as such. The commit a build stamp records, and the tables as
shipped, are checked when git can read the code's checkout; otherwise `skipped: commit check:
git cannot read a checkout at ...` and `skipped: shipped-tables checks: git cannot read a
checkout at ...` (a copy that is not a checkout, or a user git refuses: "dubious ownership").

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
| Real 29 Sep - 5 Oct 2026 contours (393,668 points, 198 lines), setup C = 0.037 added as the cron adds it, before (main `af8e1bd`) and after (this PR) | [REAL] | **main**: c1 22 lines dropped, 19 trimmed, 44,768 of 269,694 points (16.6%); c2 2 dropped, 24 trimmed, 7,122 of 123,974 (5.7%). **this PR**: c1 20 dropped, 17 trimmed, 41,738 (15.5%), 12 lines kept as beach change; c2 2 dropped, 28 trimmed, 8,163 (6.6%), 6 kept as beach change | `python3 waterline_consistency.py <contours with the setup columns, section 5> --output qc.csv --report rep.csv --plot wc --image-dir <photos> --plot-days 7` | the cron log line `waterline consistency: c1: ... ; c2: ...` |
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
* [SCRATCH] Some points that their bin-day's fit puts beyond its threshold are still kept.
  The waterline track's count (a scratch script, not in this repository, no longer at hand)
  on the real **still-water** contours (`tide_elevation_navd88`, no setup) found 175 c2 and
  67 c1 such points, on lines that fail in no bin; it was not made on the setup-added
  contours of the table above. A re-count with a wider definition (every kept point of a
  judged bin-day whose residual against that fit exceeds the fit's threshold, the time
  test's excused direction aside; a scratch script built on `waterline_consistency.run`)
  gives on the setup-added contours c2 509 of 115,811 kept points (0.44%) and c1 225 of
  227,956 (0.10%), 209 and 110 of them on lines that fail in some other bin; on the
  still-water contours it gives 530 and 229, so the count depends on the definition.
  [REAL] Three short dashes remain on the real c2 7-day map of 29 Sep - 5 Oct (the PR's
  after figure).
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
| Orientation of 16 pages (an independent reviewer's check), including the archived real 27 Sep 2026 DEM that was upside down | [SCRATCH] (real archived DEMs and synthetic pages, through a review script not in this repository) | correct signs on all 16; 27 Sep 2026: corr(page x, NNW) +1.000 -> -1.000 | scratch review script | look at the page: sea at the top, cameras at alongshore 0 on the right |
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
with readings any accepted spline has had. A spline that starts more than a day later than
the record, or has lost more than 3 days against that best, is refused and the previous
products are kept. Deleting the file accepts a shortening made on purpose.

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
| Forcing for the window | [REAL] ADCP | setup median 0.291 m (daytime 0.281), range 0.147-0.708 m at C = 0.037 | `historical_forcing.py --start 2025-01-18 --end 2025-01-23 --output-dir <dir>` (the build's window; line `Setup in window` of `forcing_report.txt`) |

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
so the numbers there will differ from the sandbox's; the product should be `complete`. After a
pull that changes the pointing rule, the same command reruns the pointing step (its version
number changed) and every step whose inputs that changes; each camera's table is followed by
each day's why.

**[NUC] Which days the first NUC build left out, and its README caveats (fixed Oct 2026).** Its
pointing step left out c1 26-28 Sep and c2 26-27 Sep, c2's on a misread horizon: under the rule
of section 4.6 c2 26-27 Sep are kept and c1 26-28 Sep are left out with their measured reason
(UNCONFIRMED, DIFFERENT POINTING, UNCONFIRMED). Its README also said things the DEM did not do.
The STORM caveat said "the frames of 2026-09-28 ... are gridded with the rest", but the DEM gridded
29 Sep - 2 Oct only (c1 28 Sep was left out by the pointing step, every c2 28 Sep frame by the
DEM's wave filter, Hs > 1.5 m); the three HIGH WAVES caveats (27 Sep 11:30, 27 Sep 20:30,
28 Sep 20:00) counted every waterline ("95 after", "4 during") and said the DEM grids "both
sides together", where it gridded 86 frames, all after; and the note on the per-day rows was
printed four times. Both caveats now read the DEM's own record (`dem/<date>_info.json`: its
window, wave filter and rejected camera-days, against the waterline file it read): the STORM
caveat names the tail days the DEM grids and, for each day it does not, why; the HIGH WAVES
caveat counts the frames it gridded and those it left out, before, between, during and after;
spans less than 48 h apart are one caveat; the note is printed once ([SYNTH]
`test_survey_products.py`, test_storm_and_wave_caveats, on a copy of that window's layout). The
comparison's FRAMES line counted the frames the DEM left out from its first GRIDDED day on only
(its photo dates), so on the NUC it said "left out 5" where the DEM's own record says 9 frames with
Hs > 1.5 m; `survey_compare.py` now looks for them over the DEM's whole window
(`window_start`..`window_end`), as the DEM and the HIGH WAVES caveat count them, and keeps its
statistics on the gridded days' frames ([SYNTH] `test_survey_compare.py`,
test_left_out_frames_window).

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
| Another C in a root that holds a build | [SYNTH] tests + refreshes over copies of the real C = 0 and C = 0.037 2026 builds | refused, exit 5, nothing changed, no lock left; the message names the root that C belongs in (`<root>_c<C>`, the base root for the station's C); a root whose name gives it to the refused C although it holds another (a root without suffix holding a C = 0 build; a `..._c0` root holding C = 0.037) is said to hold that other C, and the run is given a root beside it (`<root>_c<C>`), never the same root; the same in the dry run, the compare-step warning, the CIRCULAR caveat and the end-of-run REPLACED line; `--replace` overwrites on purpose | `test_survey_products.py` |
| A `--setup-coef` of unknown origin | tests | exit 5, nothing written, unless `--setup-fitted-to` says where it was fitted | same |
| A usage error (unknown option; a bad `--steps`, `--window` or `--utc-hours`; a negative `--setup-coef`) | tests | exit 6, nothing written: never 2, which means a disabled date | same |
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
[REAL], 0.6 (c1 0.64, c2 0.59) and 0.8 in the Jan and Mar 2025 fixture builds [SYNTH],
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

If the real C follows the slope, 0.037 leaves the January lines ~0.03-0.05 m low and lifts the
March lines ~0.05-0.07 m too high: the difference in C (January 0.0068-0.0103, March -0.0136
to -0.0108) x the window's daytime median sqrt(Hs L0) (daytime setup median / 0.037: 0.281 m ->
7.6 m in January, 0.228 m -> 6.2 m in March; the builds' READMEs use the same) x the share kept
in that month's build (0.6 in January, 0.8 in March; at 0.8 January would be ~0.04-0.06 m).
These are estimates, not measurements. They add to the ADCP still-water effect above. The
real builds on the NUC are the test.

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
  (Marconi = 1.403 x Chatham(t + 72 min) - 0.033 m, i.e. Marconi leads by 72 min; the code's
  `lag` is -72 min in `Marconi(t) = a x Chatham(t - lag) + b`) scored 0.166 / 0.174 m and is
  not used.
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
* Known limit, not fixed in this PR: the conversions are chosen by the setup error they cause,
  so a `--setup-coef 0` build of hours without the ADCP (no setup, no error to choose by)
  takes the first conversion listed (Hs linear, Tp raw), not those of the C = 0.037 build.
  No date switched on has such hours (ADCP waves throughout) and every `chatham`-era date is
  switched off, so no build is concerned now. The fix (judge them with the station's C) belongs
  with the change that switches one on: it changes `historical_forcing.py`, and so reruns the
  forcing and the detection (~50 min per two-camera week) of every date already built.

### 4.6 Camera pointing (calibrations) and the sea-horizon check

* **Calibrations per period** (`calibration/calibration_history.csv`, `calibration/README.md`):
  on 24 Jan 2025 the cameras were re-set: c1 moved 5.6 m and turned -15.3 deg, c2 moved 5.8 m
  and turned -15.5 deg (CACO03 20250123 -> CACO04 20250219). Between 23 Oct 2024 and 23 Jan 2025
  the CACO03 cameras turned ~13 deg. The live CACO05 20251113 calibration: c1 has the
  2025-11-04 values, c2 moved 3.1 m and turned -3.1 deg.
* **Why the source of the pointing matters.** A pointing fitted to the lidar
  (`fit_eo_to_survey.py`) keeps the camera position fixed, so it absorbed the 24 Jan move as
  ~22 deg of pan (`calibration/README.md`, `project_survey.py`) where the calibrations have a
  ~15 deg turn and a 5.6-5.8 m move [REAL]. With the February 2025 pointing applied to January
  (what had been used before), the Jan 18-23 lines read "lidar - water level" +1.78 m (c1) and
  +3.31 m (c2) [REAL, station runs before this PR; no file in this repository holds these
  two]; with the GCP calibration of 23 Jan, +0.46 / +0.15 m (NMAD 0.19 / 0.23 m) [REAL,
  pointing_fix_figure.py on the NUC before this PR, no setup; held in `survey_products.py`
  PRIOR_CHECKS]. The products therefore use the station's GCP calibration of each period and
  never a lidar-fitted pointing (a survey-fitted or carried pointing makes the comparison
  CIRCULAR, by rule).
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
* **Which days the pointing check leaves out (the rule of Oct 2026).** The NUC's first
  2026-09-29 build left out c2 26-27 Sep as DIFFERENT POINTING because the horizon read a
  2.6-3.3 deg change of tilt, although 73 and 375 matched features measured a turn of at most
  0.07 deg. On the real 17:00 UTC frames of 26-28 Sep (sent from the NUC; for c2 26/27 Sep the
  very frames its feature match used) no horizon is visible: the edge picked is the top of a foam
  or fog band 187-233 px below the predicted horizon, near the end of the +/-250 px search, and the
  tilt + roll fit explains it to 1.0-7.5 px (on CACO05 the horizon is inside the frame on 9 of the
  22 columns only, so any smooth edge fits) [REAL]. The step now decides each day from both
  methods (`survey_products.py` pointing_verdicts, its constants and their comments). The numbers
  below come from a scratch investigation (not code here) on the only real frames there are: the
  six storm frames, 7 Oct 2026 16:00 c1/c2 and CACO03 c1 21 Jan 2025 15:00, exact rotations of
  them, planted rotations and synthetic weather [REAL frames, SYNTH changes, SCRATCH analysis].
  * **Matched features decide, both ways, at 40 inliers** (pointing_check.py's own minimum; the
    trust threshold measured): on real pairs with a known answer the rotation erred by <= 0.063 deg
    (69 results at >= 40 inliers); on 402 planted 0.05-2.5 deg rotations by <= 0.053 deg (a frame
    against itself) and <= 0.101 deg (against another day; 40-59 inliers: <= 0.10 deg on 45
    trials); no matched trial under 0.45 deg was called MOVED (0 of 229), every one of 0.55 deg or
    more was (157 of 157). Below 40 the errors grow (20-39 inliers: up to 0.18 deg; 10-19: 0.58;
    under 10: anything). A turn of more than 0.1 deg is the scene (dirt matches at zero) and
    measures the day's change from the reference directly: a change of tilt or roll over 0.3 deg
    (the step's tilt/roll tolerance, as for the horizon) or of pan over 0.5 deg (pointing_check's)
    is a DIFFERENT POINTING whatever the horizon reads; a smaller turn is a 'small change', kept (a
    review probe: 1500 inliers at +0.15 deg against a clear-looking +2.0 deg jump). A turn of
    0.1 deg or less keeps the day ('matches') whatever a horizon reading that is NOT a clear horizon
    says, and the why gives both readings, as long as the turn is more than its own noise (window
    dirt, below); against a CLEAR horizon edge that jumps a zero turn is not enough: the day is
    UNCONFIRMED. The first version of this rule judged a turned match's tilt and roll by 0.5 deg and
    let a clear horizon decide between: review 2 turned the real c2 27 Sep storm frame by exact ray
    rotations, and +0.30, +0.40 and +0.45 deg of tilt and +0.45 deg of roll, read by the features
    as 0.303 / 0.402 / 0.446 / 0.388 deg on 90-105 inliers against the misread horizon, were 'small
    change', kept (0.45 deg of tilt is ~0.6-0.75 m of DEM at 250-350 m; the deployed rule had left
    them out for the misread horizon's jump); and in a clear window the same 0.300 deg read by the
    features was a 'small change' on one day and DIFFERENT POINTING on the next, as the horizon read
    0.298 or 0.300 deg. A turned match's tilt/roll is now judged on its own. An earlier
    proposal to require 60 inliers was not taken: it was chosen to make c1 27 Sep (48 inliers)
    weak, but that day's -1.21 deg of pan is real: tracked dune landmarks and a more sensitive
    match read -1.20 to -1.25 deg on 26-28 Sep, and the camera was back by 29 Sep.
  * **Window dirt** is the one known way a confident match is wrong (drops at the same pixels in
    both frames match at zero turn): a near-zero winner with a runner-up consensus of 40+ inliers
    turned more than 0.3 deg from it is 'ambiguous' and decides nothing (planted heavy drops: 17 of
    120 trials read 0.00 deg on 388-1573 inliers, each with a runner-up of 64-561 at the true turn;
    without dirt the runner-up had <= 35 inliers on 551 trials and the real pairs). The runner-up
    shows the dirt only while the scene still matches. With the land fogged the dirt is the only
    consensus: in a review scenario (the real 7 Oct 16:00 frames with the same synthetic drops on
    every day, an exact +0.6 deg tilt on a fogged 2 Oct) c1 2 Oct read 68 inliers at 0.00 deg with
    a runner-up of 6 while its clear horizon read the true +0.60 deg, and the first version of this
    rule kept the day. A dirt winner can hold 388-1573 inliers, so no count tells it from the
    scene: a zero turn never overrules a clear horizon jump (that day is now UNCONFIRMED, left out,
    as the old rule left it out). [SYNTH] `test_survey_products.py`, test_feature_runner_up, runs
    pointing_check.measure(second=True) on a drawn scene turned 1.0 deg under drawn dirt: the
    scene wins (runner-up: the dirt at zero), heavier dirt wins ('ambiguous'), or the scene is
    lost ('ok' at zero turn, runner-up 2) and the verdict against a clear jump is UNCONFIRMED.
  * **A zero turn within its own noise shows nothing.** Review 2 found the same dirt against a
    horizon that is NOT clear, or with none fitted, which is what a storm gives: in DIRT2 (the real
    c1 26/27 Sep storm frames, c1 turned -1.2 deg, the same 500 synthetic drops on every frame)
    113 and 165 inliers read 0.00 deg with a runner-up of 3-5 and both days were 'matches', kept;
    in DIRT2_ST10 c2 28 Sep, turned +1.0 deg, read 0.00 deg on 244 inliers (runner-up 26 at the
    true turn, under 40) with its horizon fit failed, kept. The deployed rule left out the first
    two. Dirt moves with the camera, so its consensus sits at EXACTLY zero turn, while a scene
    matched between two days reads the camera's own small turns: `pointing_check.turn_sigma` gives
    the fitted turn's distance from none in units of its standard error (from the inliers' residual
    scatter; a chi of 3 degrees of freedom for noise alone, 99.9% below ~4.0), written to the CSV as
    `turn_sigma`. Measured on the station-floor venv [REAL frames, SYNTH dirt, SCRATCH]: the drops
    alone 0.5-2.7 (21 matches of 113-483 inliers, the scene turned or lost); the real c2 26-28 Sep
    storm frames, the camera at rest, against the real 7 Oct 16:00 frame and a copy of it 4.9-16 (6
    matches of 25-102 inliers; the lowest, 28 Sep, is a 0.06-0.07 deg roll, the turn a match
    constrains least); copies of one frame (no true turn, as the sandbox's stand-ins) 1.2-5.0. A
    second measurement with other drops (500 soft blobs of 3-12 px on the land, the same pixels in
    the reference and the probe; system Python) read the drops alone 0.4-3.3 (18 consensuses of
    16-117 inliers, the scene turned 0.6-1.0 deg or lost; 2.6 at most on the 8 of 40 or more
    inliers; drops over a scene at rest, both at zero, 3.2-5.9, where either answer is right), the
    real c2 storm frames against the real 7 Oct frame 6.5-13.4 and against each other
    (26, 27 and 28 Sep, the camera at rest: 64-292 inliers) 15.6-18.5, and copies of the 7 Oct frames
    1.0-2.4: real scene matches between two days read 6.5 or more, the drops 3.3 or less. With
    those standard errors scaled to the NUC's inlier counts, the NUC's own c2 26, 27 and 28 Sep
    matches (0.07 / 0.02 / 0.01 deg on 73 / 375 / 393 inliers; their 3-decimal d columns turned into
    the rotation and weighed by the normal matrix and residual scatter of the same real frames
    matched here, against 7 Oct and against each other) project to ~17-18 / ~9 / ~8-9.5. A zero-turn
    match under 4.0 (TURN_SIGMA_MIN) is therefore no reading at all: the day is judged as one without
    matched features (a clear horizon at the reference's pointing: 'horizon only', kept; a misread
    horizon that jumps: UNCONFIRMED; no horizon fitted: 'unchecked', and UNCONFIRMED next to a
    left-out day). The reference day's own frames are exempt (two frames of one day match at zero
    by nature), and a row of an older build without `turn_sigma` counts as a reading (the NUC's
    2026-09-29 table below). What it does not catch: a dirt consensus of thousands of inliers may
    read up to ~5 by the images' own systematic errors (the copies did, at 3600-4800 inliers) and is
    then still taken for the scene; under a clear horizon that does not jump, dirt can still hide a
    pan (no horizon sees one); a day nothing could read at all (its horizon unfitted, its match
    blind) is 'unchecked' and kept unless a left-out day is next to it (review 2's DIRT: c1 26-28
    Sep, turned -1.2 deg, the dirt also broke all three horizon fits: kept, as by every earlier
    rule; this round's DIRT2_ST10, below: c2 tilted +1.0 deg, the scene and the drops both under 40
    inliers, the fits failed: kept, as by the deployed rule); a real match under 4.0 is not a reading, so its day is judged without it (a lost day at
    worst, not a wrong one). Only the NUC's next build of 2026-09-29 gives the real values for its
    storm days (`turn_sigma` in `pointing_c2.csv`).
  * **A horizon reading counts only on a clear edge**: prominence (its gradient over the median of
    its search) >= 12, the strongest rival edge 20+ px away <= 0.65 of it, sky over sea on >= 2/3
    of the columns, and > 50 px from the end of the search. The six storm frames read 2.2-6.7 /
    0.57-0.96 / 0.00-0.22 (26-27 Sep) / 15-62 px (26-27 Sep) and each fails two tests or more; the
    real clear frames 61-85 / 0.06-0.24 / >= 0.89 / 174-214 px; exact rotations of the clear frames
    by up to 2 deg of tilt or roll pass (not clear: c2 tilted -2 deg, whose horizon leaves the top of
    the frame, and +3 deg of tilt on c2 and the Jan c1, where it reaches the end of the search;
    matched features measure such turns). On synthetic weather put on the clear frames, the tests
    remove 135 of 206 false jumps; a horizon-like fog front can still pass, which is why each clear
    reading is checked against the matched features (below). The fit residual (FIT_OK_PX 8 px) is not tightened: the NUC's
    storm days fitted to 1.6-5.0 px against 0.2-0.4 px on its clear days, but single storm frames
    fit to 0.97 px and the fixture's rendered horizon to 1.4 px.
  * **Each day is judged against the reference, as the features are.** The anchor is the
    reference's tilt/roll as the horizon reads it: each clear day with matched features gives its
    horizon less the tilt/roll its features measured, and the anchor is the median of a group of
    these that agree within 0.3 deg: the reference day's when another day confirms its reading, else
    the largest (else the reference day's own clear reading; else the
    largest group of clear readings, with nothing tying it to the reference; else, no day clear,
    the fitted readings' median). A day of the first kind outside that group is CONTRADICTED (its
    horizon and its features disagree about its change): it is no clear day, and a day without
    features reading the same jump is UNCONFIRMED (a fog front on 26-28 Sep that passes the edge
    tests, review simulation sB: 27 Sep's 154 matched features measure its real pan and no tilt
    where the front reads +2.9 deg). The first version judged every day against the window's
    median, which falls between two pointings when a camera moved on half the clear days: in review
    scenario E (real frames, an exact +0.6 deg tilt on 30 Sep and 1 Oct with the land fogged) c1
    read 0.91/1.51/1.51/0.90 deg, median 1.21, and both moved days were 'horizon only', 'tilt and
    roll agree with the other days'. Against the reference day they are DIFFERENT POINTING. The
    anchor group is the reference day's own whenever another day confirms its reading, and only
    otherwise the largest: in review 2's simX2 (a planted clear-looking edge at +130 px on c2 26, 27,
    28 and 30 Sep with features matched at zero turn; 29 Sep, the reference, and 1 and 2 Oct real
    clear frames) the four-day front outvoted the reference day, so 1 and 2 Oct were UNCONFIRMED and
    left out and the README gave the front as the calibration offset (131 px, +1.96/+0.02 deg; the
    true one is 34 px, +0.56/-0.67). A lone reference reading that the other days outnumber is
    still the odd one (a fog front on the reference day), and the step and the README then warn
    that the days are judged against a group that leaves the reference day's own clear horizon out.
  * **Misread days are not clear days.** The window's offset against the calibration and the
    README's 'N day(s) with a clear horizon' are those of the clear days at the reference's
    pointing (neither jumping nor contradicted). NUC c1: 4 clear days, +0.93/+0.24 deg, 42.5 px (it
    was 7 days, +0.96/+0.23, the three storm misreads included); c2: 4 days, +0.54/-0.68 deg,
    32.5 px (was 6, +0.56/-0.64). A window with more misread days than clear ones no longer flags
    its clear days. With no clear day at all the readings are compared with each other only and no
    offset against the calibration is reported.
  * **Without matched features** a jump of tilt or roll over 0.3 deg from the reference is a
    DIFFERENT POINTING only on a clear edge that two or more days share, none of them contradicted
    (real re-aims persist; the storm misreads changed by 0.6-3.1 deg from day to day); otherwise the
    day is **UNCONFIRMED**: a change one method reported that nothing could confirm or rule out. A
    day on which nothing could be read (no horizon fitted, features unmatched) next to a left-out
    day is UNCONFIRMED too, as are runs of such days: nothing shows the camera back (review 2's simR,
    the real storm frames as 17:00 with perturbed copies as 16:30/17:30: c1 28 Sep's 3-frame fit
    failed at 14.6 px, and the day was 'unchecked', kept and gridded while c1 was turned -1.2 deg;
    on the NUC its fit held, 3.5 px, and the day was UNCONFIRMED anyway). Never the reference day
    itself: the others are measured against its frame, so a neighbour that differs from it says
    nothing of it. An 'unchecked' day whose horizon fitted within 0.3 deg of the reference on an
    edge that is not clear stays kept.
    UNCONFIRMED days are left out like DIFFERENT POINTING ones (`--keep-moved-days` keeps both);
    only measured DIFFERENT POINTING days count towards 'most days differ from the reference'. When
    the days left out are more than half of the window's, the step and the README say so with the
    count, how many of them rest on a horizon reading alone that is not a clear horizon, and the
    options there are: `--keep-moved-days` keeps every such day, and there is no per-day option
    (`--window FIRST LAST` narrows the photo days; the warnings named `--skip-days`, which is
    detect_original_view.py's, not survey_products.py's) (a review probe: a hazy week, one clear day and six hazy
    readings 0.33-0.55 deg off it, lost six of seven days with no word). Such days are still left
    out: keeping them when they are most of the window was tried and dropped, because the sandbox
    build with the six real storm frames (a 4-day window there) then KEPT c1 26-28 Sep, the very
    days c1 turned -1.2 deg. No saved why says 'left out' (with `--keep-moved-days` the README
    quoted a why that said so for a day it called KEPT). Which days the README calls left out or
    KEPT, and why, comes from the detection's own stamp (its left-out days, and whether it ran with
    `--keep-moved-days`), never from the options of the run that rewrites the README: a pointing-only
    or partial build said both 'KEPT in the detection (most days differ ...)' and 'left out of the
    detection' for the same days (review 2); a build whose detection has not run now says so, and
    so does one whose detection ran with other pointing results (the pointing step run again on its
    own since, e.g. after a change of the rule: its stamp's left-out days are not those the current
    results give), naming what that detection left out.
  * **The NUC's 2026-09-29 table under the rule** (`test_survey_products.py`,
    test_pointing_verdicts, on the NUC's rows with the edge numbers of the real frames): c1 26 and
    28 Sep UNCONFIRMED (no matched features, a misread horizon), 27 Sep DIFFERENT POINTING (the
    features' -1.21 deg of pan on 48 inliers; the horizon's +4.82 deg named a misread), 29 Sep -
    2 Oct matches; c2 every day matches (26/27 Sep kept, the horizon overruled; 28 Sep: its fit
    failed, 393 features decide). Left out: c1 26-28 Sep, as before, with the measured reason; c2
    none (was 26-27 Sep). The NUC's rows predate `turn_sigma`; with the values the real storm frames
    project for c2 26-28 Sep (~17 / ~9 / ~9, the window-dirt bullet above) the verdicts are the same,
    and at a dirt-like 2.0 c2 26-28 Sep would be UNCONFIRMED and left out (both pinned in the test).
  * **The six real frames end to end** (a sandbox build of 2026-09-29 with them and a 7 Oct frame
    standing in for 29 Sep as its reference): c1 26-28 Sep UNCONFIRMED (2-9 inliers against 7 Oct; the
    NUC's 29 Sep reference gave 48 on 27 Sep), c2 26 Sep UNCONFIRMED (25 inliers against 7 Oct; 73
    against the NUC's 29 Sep), c2 27/28 Sep matches (102 / 69 inliers, <= 0.08 deg, against the
    7 Oct 16:00 frame itself; a review build that used a gamma and noise copy of it as the
    reference read 106 / 77: the counts depend on the stand-in). The only real clear 2026 frame in
    the sandbox is 7 Oct 16:00 (c1 and c2): the 198 stand-in frames named 29 Sep - 5 Oct are byte
    copies of it (md5; that folder has since been removed, and this round's rebuilds link the
    7 Oct frame as 29 Sep 14:00 instead), and the lens4 'Oct 5' frames are renders. The same days
    come out of this build under the final rule, c2 27/28 Sep at turn_sigma 14.0 / 6.5 (its window
    holds 4 days, so c1's three left out days are most of it and the step says so).
  * **Review scenarios, end to end** (`survey_products.py --steps pointing` on the Python 3.8 /
    numpy 1.17 / OpenCV 4.5 venv; real 7 Oct 16:00 and storm frames, exact ray rotations, fogged
    land, the same synthetic drops on every day; scratch, outside this repository) [REAL frames,
    SYNTH changes]: dirt with a fogged, tilted 2 Oct (G): c1 2 Oct UNCONFIRMED, left out (was
    'matches'); a 0.6 deg tilt on 30 Sep and 1 Oct with the land fogged (E): DIFFERENT POINTING on
    both cameras (were 'horizon only'), the one-day 1 deg roll UNCONFIRMED; dirt on every day with a
    pan, a tilt and a fogged tilt (F): all three found on c2 (the fogged one confirmed by the
    matched tilt of 1 Oct), c1's fogged one UNCONFIRMED. On c1 in E and F only three days are
    measured (the storm days are not), two of them moved, so the existing 'most days differ from
    the reference' guard keeps them, with its warning; c2 leaves them out. Scenarios A-D (moves of
    0.3-1.0 deg on clear and storm days): every clear day moved by more than 0.3 deg of tilt or roll
    or 0.5 deg of pan is DIFFERENT POINTING, and the 0.3 deg ones are 'small change'; storm days
    whose features matched are measured (c2 27/28 Sep in A); storm days without matched features
    are UNCONFIRMED, except a pure 1 deg pan with the horizon unfitted (A, c1 28 Sep: 'unchecked',
    nothing there can see a pan). The 'most days differ' guard turns on where the scenario moved
    most measured days (A c2, D c1).
  * **Review 2 scenarios, end to end with this rule** (`--steps pointing`, the station-floor venv;
    scratch, outside this repository) [REAL frames, SYNTH changes]. The review's own folders were
    removed from the scratchpad, so its scenarios were re-made from its recipe: the six real storm
    frames as 26-28 Sep 17:00, the real 7 Oct 16:00 frame (gamma, noise and JPEG per frame) as 29 Sep
    - 2 Oct at 16:30 and 17:00, exact ray rotations, and for DIRT2* the same 500 soft drops on the
    land in every frame, the reference included. Each was also run with the deployed rule (ddde694)
    on the same photos. c2 27/28 Sep turned +0.40, +0.45 or +1.0 deg of tilt or +0.45 deg of roll
    read 0.39-1.01 deg on 63-95 inliers: DIFFERENT POINTING on both days, left out (deployed, +0.45
    deg of tilt or roll: 27 Sep DIFFERENT POINTING, 28 Sep a 'small change', kept). +0.30 deg: 27 Sep
    read 0.303, DIFFERENT POINTING, 28 Sep 0.297, a small change (deployed: the same). A 0.4 deg pan:
    a small change on both days (the pan is judged by 0.5 deg). A clear-window 0.3 deg tilt on 1-2 Oct
    (T03): a small change on both days on both cameras (deployed: the same); 0.6 deg (T06):
    DIFFERENT POINTING. DIRT2 (c1 turned -1.2 deg on 26-28 Sep, under the drops): c1 27/28 Sep read
    0.00 deg on 50 / 66 inliers at turn_sigma 1.3 / 2.2, UNCONFIRMED and left out with 26 Sep
    (deployed: 28 Sep 'matches', kept); c2, at rest, 26/27 Sep 'matches' at 5.7 / 6.6 and 28 Sep
    blind at 3.3 with its fit failed, 'unchecked', kept. DIRT2_ST10 (DIRT2 with c2 tilted +1.0 deg on
    27-28 Sep): with these drops neither the scene (32 / 29 inliers, at the true turn) nor the dirt
    reached 40 inliers and both horizon fits failed (52 / 71 px): 'unchecked', KEPT, as by the
    deployed rule. A turned camera whose features stay under 40 inliers and whose horizon cannot be
    fitted is not seen by either rule (the review's own DIRT2_ST10, with other drops, had a dirt
    winner of 244 inliers, blind at 2.1, next to an ambiguous 27 Sep: UNCONFIRMED there). Without
    drops c2 27/28 Sep match at turn_sigma 12.2 / 4.9 and stay 'matches'; the clear stand-in days read
    'horizon only' where the first version read 'matches' (copies of one frame: blind, 0.8-3.5), kept
    either way. Earlier in this round: a re-made simR (the real 17:00 frames with perturbed copies at
    16:30/17:30): c1 28 Sep's 3-frame fit failed (11.8 px) and is UNCONFIRMED next to 27 Sep (was
    'unchecked', kept); c2 27/28 Sep read turn_sigma 14.6 / 5.6, 'matches'. A re-made simX2 (a
    planted clear-looking edge at +130 px on c2 26, 27, 28 and 30 Sep): the window's offset is the
    true 34 px over the 3 clear days (was 131 px); the stand-ins' matches are blind, so the four
    front days, sharing a clear jump with no feature reading, read DIFFERENT POINTING and, as 4 of
    the 7 measured days, the 'most days differ' guard keeps them, with its warning (with real
    day-to-day matches, as in the pinned test, they are contradicted and UNCONFIRMED).
  * [NUC] With c2 26-27 Sep kept, their live rows reach the detection, and the DEM unless its wave
    filter (offshore Hs > 1.5 m) drops them: the 2026-09-29 numbers of section 3.3 may move.
  * [SYNTH] The sandbox fixture's rendered horizon sits near the limits: Jan c1 and Mar c2 pass
    (prominence 14.9-22, rival 0.49-0.62), Jan c2 does not (prominence 6.5-8.0), so the fixture's
    Jan c2 days now read 'unchecked' (were 'horizon only') and its rendered offset is no longer
    reported as one. Nothing is left out on either fixture date, as before.
  * **Sandbox rebuilds with this rule** (the commands of section 3, fresh output roots; outputs
    outside this repository): 2025-01-23 and 2025-03-06 [SYNTH] and 2026-09-29 with and without the
    7 Oct stand-in photos [REAL rows] give the same left-out days (none) and the same headline
    numbers as section 3 (every DEM - survey and waterlines - survey value). What differs: the
    pointing CSVs carry the edge and runner-up columns; Jan c2 'unchecked' (above); a day with no
    second photo says so, not 'features unmatched'; the HIGH WAVES caveats count the frames the DEM
    gridded (Jan: 40 before, 60 after, the 20 during left out by its wave filter, where it said "the
    DEM grids the frames of both sides together: ... 20 during"; Mar: 36 after, not 40, and 4 left
    out); the 2026-09-29 STORM caveat says the DEM grids no frame of 26-28 Sep (the sandbox rows
    start 29 Sep). The pointing step took 29 s (Jan) and 39 s (Mar, on a loaded machine) against
    41 and 26 s before: the second consensus and the edge numbers cost little. Rebuilt again with
    the final rule (turn_sigma, review 2; Jan, Mar and 2026-09-29 in their own output roots, the
    detection up to date, so only the pointing step and what follows it reran; plus the photo
    stand-in and six-storm-frame builds): the same left-out days and the same headline lines. The
    fixture's Mar c2 days, whose rendered photos match at 51-59 inliers and 0.001-0.022 deg (no
    drift between days: turn_sigma 0.7-3.4), now read 'horizon only' instead of 'matches', kept
    either way, and the README adds 'pan not measured' for them. A real camera that did not move
    at all between two days reads the same way: its day loses the feature reading, which only
    matters when its horizon also jumps on an edge that is not clear (then UNCONFIRMED).
* **Re-check.** Each product README, `sea horizon vs the calibration` and the KNOWN CAVEATS
  pointing block; `pointing/pointing_<cam>.csv` (with the edge numbers `horizon_edge`,
  `horizon_prominence`, `horizon_rival`, `horizon_sky_above`, `horizon_to_limit_px`, and the
  match's `turn_sigma`, `runner_up_inliers`, `runner_up_d`); the step prints each day's why under
  its table.

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
branch, the same commands for both.

**The setup columns first.** `contours.csv` must carry `setup_correction_m` and
`beach_elevation_navd88`, as `extract_elevation_contours.py --setup-coef 0.037 --waves ...`
writes them on the station: per row, setup = `round(0.037 x sqrt(Hs x 9.81 x Tp^2 / (2 pi)), 4)`
and beach = `round(tide_elevation_navd88 + setup, 4)`; where Hs or Tp is missing the setup is
blank and beach = tide. The sandbox copy of the rows (`contour_points_timex.csv` rows from
before the setup went live) had no such columns; they were added with a scratch script (not in
this repository) that applies exactly that formula and rounding, so it matches the station's
columns to 0.1 mm (393,668 rows, setup 0.171-0.489 m). On the NUC use the station's `contour_points_timex.csv` when its rows of the week carry
`beach_elevation_navd88`; otherwise add the two columns as above. Then:

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
* **An unbiased 2026 level.** The 2026 RTK stops at +1.17 to +1.80 m; every C = 0 line, and some
  C = 0.037 lines, lie below it. A calm, spring-low-tide RTK survey reaching the low-tide line,
  on two or more days with different waves, would give an INDEPENDENT level check of the whole
  intertidal and pin C and any level offset separately.
* **The ADCP-era still-water share s.** Measured by `historical_forcing.py` on the NUC.
* **The 2026 pointing.** A sea-horizon (or GCP) check of the CACO05 20251113 calibration would
  tell lens model from pointing.
* **The pointing rule on real winter photos.** Its thresholds rest on six storm frames of one
  storm, three real clear frames and their rotations (section 4.6). The real Jan and Mar 2025
  horizons and feature matches are on the NUC only: each build's `pointing/pointing_<cam>.csv`
  gives the edge numbers and inliers per day. The NUC's edge numbers are 3-frame medians, where
  the storm values above are single 17:00 frames. Why c1 turned ~1.2 deg on 26-28 Sep 2026 and
  was back by 29 Sep (wind load, or a slip reset) is not known.
* **Window dirt in a storm.** The zero-turn guard (`turn_sigma` under 4.0 shows nothing) rests on
  synthetic drops on real frames and on real cross-day matches of one camera's storm frames (c2 26-28
  Sep against 7 Oct and against each other; section 4.6); real
  salt or drops on the housing window have not been seen here. The NUC's next 2026-09-29 build
  writes the real `turn_sigma` of its storm days (projected ~17 / ~9 / ~9 for c2 26 / 27 / 28 Sep; 27 and 28 Sep have the least margin). Not caught: a
  dirt consensus of thousands of inliers reading above the threshold, and dirt hiding a pan under a
  clear horizon that does not jump.
* **The detector's rate on the NUC10i3.** Measured here at 20-25 s per photo on a 2.8 GHz
  Xeon; the NUC prints its own rate.
* **Reviews.** The waterline work had three review rounds (the last: all lenses "ship"), the
  survey products three verification and three fix rounds, the combined branch two review
  rounds (interactions; operator and documents) and two fix rounds, and the setup refit an
  adversarial check that found the survey-floor bias. After the PR was opened, a final
  verification (a review of `3c03f0e` and a fact-check of this file against the data) found
  one wrong refusal advice (a root whose name does not say which C it holds), a backwards time
  shift in the printed Chatham formula, and wrong or untraceable numbers and tags in this file;
  they were fixed in the commit after `b54c255`. That last fix round was checked by its own
  tests, the full check list and re-runs on the real data, not by a further independent review.

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

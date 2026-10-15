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
python3 test_historical_forcing.py               # ends "101/101 checks passed"
python3 test_dem_day_test.py                     # ends "all DEM day-test tests passed" (section 2.3a)
cat archive/gnssr_record_start.txt               # note "YYYY-MM-DD N"; see section 2.5
```

What a survey-date build reruns after this pull: a step reruns when a script in its signature
changes, and everything after it with it. This PR changes `historical_forcing.py` (forcing),
`pointing_check.py` and the pointing rule (`POINTING_VERSION` 11: pointing), `dem_from_contours.py`,
`dem_figure.py` and `asc_to_geotiff.py` (DEM) and `survey_compare.py` (plan-view map, compare).
So the next build of a 2025 date reruns every step once, its forcing first and so its detection
too: here, under load, a full rebuild took ~38 min for 2025-01-23 (two cameras, a week) and ~20 min
for 2025-03-06 (c2). The live 2026-09-29 date keeps its forcing (the live forcing does not use
`historical_forcing.py`), reruns its pointing check, and reruns its detection (a georectification
of the station's own rows, minutes) only if the new rule leaves out other days than the detection
on disk did; then its DEM, plan view and comparisons. The GeoTIFFs are then tagged EPSG:6348
(section 4.8), and a further run skips every step. (PR #87's own note, that its pull reran only
the plan-view map and the comparisons, was checked on a 2026-09-29 build then.)

`view_reproject.py --self-test` needs scipy, which does not import on the NUC. That file
is unchanged by this PR, and the survey products use only its `ground_to_pixel`, which
does not need scipy. **Do not run `pointing_check.py --self-test` on the NUC**: it rotates a frame
with `view_reproject.reproject_image`, which needs scipy; without a scipy that imports (none, or
one installed but broken) it now says `SKIPPED (nothing checked)` and exits 0 (it needs `--image
<a c2 frame>` too). The survey products'
pointing step does not use it. Here, with scipy, on the real 7 Oct 2026 16:00 c2 frame: four
planted turns (up to 2.5 deg of pan) recovered to 0.000 deg, PASS.

The self-tests give the same result on the NUC as in a clean checkout: they run on fixtures
and folders of their own, with the station's search folders and archive files switched off
(`test_historical_forcing.py`'s last check confirms that no build named a station file). The
four `test_*.py` files above run in a temporary folder that is also their current folder; the
three `--self-test` runs write only into temporary folders of their own and look for nothing in
the current folder but matplotlib's own `matplotlibrc` (checked with strace). What they read of
the station is said. What `test_survey_products.py` reads of the station, all of it: the tables
`surveys/survey_dates.csv` and `surveys/surveys.csv` and the live cron `waterline_timex_cron.sh`
(below), the real 2026-09-29 RTK file where `surveys/surveys.csv` says (below), and the
repository's own files as COMMITTED (git): its fixtures copy the calibration folder from `git
ls-files calibration`, never a file the station keeps there untracked and never a file named like
the held 23 Oct 2024 survey (review round 1 of the integration: a copy of the whole folder had
opened every untracked file there). A test run stopped with Ctrl-C, SIGTERM or SIGHUP removes its
temporary folder (with the synthetic canary files, named like the held survey, that the blind
tests put there); only a SIGKILL leaves it in `/tmp`.

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

What was asked for in the tables (19 Mar 2025 switched off; 23 Oct 2024 switched on and BLIND,
its survey withheld, calibration A from the table; the windows, the 29 Sep 2026 RTK row and its
text, envelopes from the other lidar, every calibration in the repository, the Oct GCPs
PARTLY-CIRCULAR with no path, held by Chris Sherwood, every row with all its columns) is checked on the
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

**The DEM day test's trend message on the station's DEMs is a NOTE, not a WARNING.** The day test
reports when a camera's camera-day offsets run one way with time (section 2.3a). The cron passes
`--trend-note`, so on the whole-archive and 7-day DEMs this is printed as `NOTE: ...` in
`logs/timex_cron.log`, which `station_status.py` does not count; the same record is in each DEM's
`_info.json` (`filters.day_test.trend_warnings`). Why not a WARNING there: the replay of the real
5-11 Oct 2026 week gives the message for c1, and the review of this PR found it would recur on
roughly half of all 7-day windows, following the spring-neap tidal range (the real c1 5-11 Oct and
the static Mar 2025 replica both track it: a detection or setup bias that grows with the tidal range
gives this pattern; not yet tested against the tide record). On the whole-archive DEM a sustained
beach change (a 0.5 m drift over 13 weeks, simulated) gives it on every run for as long as the change
lasts. As WARNINGs these would hold the station status at WARN on every run and hide a real WARN.
The day test's other message (a day it could not judge) stays a WARNING, and survey-date builds
(`survey_products.py`) keep the trend as a WARNING among their caveats. The day-test text says what
the method alone gives on a static synthetic beach (the same pattern at rho -0.96 and 0.165 m over 7
days), so a camera fault or a beach change needs more than the message.

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

**Two checks before the per-line pass** (Oct 2026, after the first survey-date builds on the
station's real photos). What they drop is not in the per-line pass's reference.

* **The day check** (`judge_days`, `--max-day-residual 1.0`): a camera-day whose lines sit, on
  median, more than 1.0 m of water level from where the OTHER days' lines (within 3 days) put
  their level is left out whole, with a WARNING naming the camera and the day. On 22 Jan 2025
  snow lay on the upper beach and both cameras drew every low-tide line of the day along the
  snow's edge, ~250 px landward of where 18-21 Jan put the same level: the lines read 1.7-1.9 m
  low against the 23 Jan lidar, the DEM's c1 cells -1.80 m. The per-line pass kept them: the
  day's own ten lines were each other's reference (weighted by closeness in time), and on the
  window's last day no later line said otherwise; the DEM's day test could not compare them
  (no cell shared with another day). A line landward of every other day's line and LOWER than
  the highest of them is out of order; the beach is taken to keep rising landward at half the
  rate of the other days' curve (as seaward, step 4), so a day of such lines reads about what
  it is off by, not the 0.2 m the bound alone gives.
* **Cross-shore runs** (`waterline_track.py`; here only when the file carries the ground, as
  the survey-date products' files do; `dem_from_contours.py` applies the same rule to every
  DEM, the live cron's included): a line's track that runs more than 45 deg from the shore
  over 10 m on the ground, and the rest of the line beyond it, is dropped. Where the edge it
  follows fades out, c2's line tracer carries on along one pixel row out over the water.

| what | data | number | reproduce here | check on the NUC |
|---|---|---|---|---|
| Snow edge: every c1 line of one day along the +1.8 m contour, on the window's last day and mid-window | [SYNTH] self-test | 14,568 of 14,568 points dropped both times, WARNING given; honest lines 0 of ~188,800 | `python3 waterline_consistency.py --self-test` | same command |
| Cross-shore tails: every 3 h a c2 line runs flat from column 1300 to the right edge | [SYNTH] self-test, on the ground | 4,327 of 4,327 points more than 20 px off dropped; honest 176 of 48,610 (0.36%) | same | same |
| Every earlier case of this table with both checks on | [SYNTH] self-test | unchanged (the storm cuts: 0 honest points dropped) | same | same |
| 2025-01-23 (18-23 Jan), the survey-date build's lines | [REAL] owner's files, 15 Oct 2026 | day check: c1 22 Jan (10 lines, -1.91 m) and c2 22 Jan (2 lines, -2.34 m) left out, nothing else; cross-shore: c2 4 lines, 3,643 points | `python3 waterline_consistency.py <date>/waterlines/contour_points_ground.csv --output f.csv --report r.csv` | `survey_products.py --date 2025-01-23 --steps filter,dem,maps,compare` |
| 2025-03-06, 2024-10-23, 2026-09-29 (the same) | [REAL] | no day left out on any; cross-shore: Mar c2 2 lines (891 points), Oct c2 15 of 47 lines (7,793 points), Sep 2026 c1 4 lines (1,605), c2 19 of 42 (5,560); no Oct c1 line | same | same |
| What it did to the 2025-01-23 DEM against the Jan lidar (measured AFTER the rules were set, from the photos and the contours alone) | [REAL] | DEM - lidar: median +0.023 -> +0.044 m, NMAD 0.271 -> 0.218, RMSE 1.005 -> 0.207 m, within +/-0.20 m 51% -> 64%; by camera c1 -1.799 -> -0.218 (52 cells), c2 +0.117 -> +0.131 (106); lines c1 -1.826 -> -0.180 (16 frames), c2 +0.099 -> +0.135 (21) | `dem_from_contours.py` + `survey_compare.py` as the build runs them, on the filtered file | the rebuilt README |

The snow failure was found BY the Jan lidar comparison, so the Jan numbers after the fix are
not a blind test of these two rules; the rules were set from the photos (the snow is plain in
c1's photo map and in the detector's overlays of 22 Jan) and the contours, before the lidar
comparison of the filtered lines was run, and they change no line of the 2025-03-06 build but
4 cross-shore runs (DEM - Mar lidar -0.293 -> -0.288 m).

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

### 2.3a The DEM's day test (`--max-day-offset`, every DEM: cron and survey dates)

**Claim.** The day test leaves out a camera-day whose waterlines sit, on median, more than the
limit (0.15 m) from the DEM of the other days, judged against days that are themselves kept,
and says when it could not judge or when a camera's days drift. Before this PR it ran once,
each day against the DEM of all the others, so the days it was about to reject shaped that
DEM: on the station's whole-archive DEM of 13 Oct 2026 c1's pre-storm days (+0.27 to +1.40 m)
pulled c1's cells up and four post-storm c1 days were rejected as LOW [REAL].

How it decides (`dem_from_contours.day_consistency`; the rule as printed is
`DAY_TEST_METHOD`, also in `<stem>_info.json` `filters.day_test.method`):
* **Iterative.** First pass: every camera-day against the DEM of all the other days. Then the
  day furthest over the limit is left out, the others are measured again without it, and so
  on until no kept day is over (ties by name; the result does not depend on sample order).
* **Judging.** A day is left out only when the test can judge it: at least 3 testable days
  are kept, and in at least half the cells where it is compared its DEM holds two or more
  other kept days, or one other day (a **lone partner**) that the rest confirm (that day's
  offset from them, without the day judged, is within the limit). Otherwise (two days that
  only see each other, a day whose cells only a day as far off reaches, the last two days)
  the test cannot tell which day is off.
* **Unjudged days** keep the first pass's verdict (the one-pass test's: over the limit there,
  left out; within it, kept) and are tried again in later rounds; once the test can judge
  one it is treated like any other. When the final DEM still leaves out a day it could not
  judge, keeps one over the limit it could not judge, or keeps a day **untested** (testable in
  the first pass, but too few comparable samples are left once the days left out are gone;
  it keeps its place whatever its first-pass offset), the log prints one `WARNING: the day
  test could not judge every camera-day -- ...` line naming them, with their offsets.
* **Take-back.** Each day left out is then measured against the DEM of the kept days; the
  closest one within the limit is taken back, one at a time, if that puts no kept day over the
  limit. A day within the limit whose return would put a kept day over stays out and is
  recorded as a **refused take-back**: it is NOT "off the rest by more than the limit", and no
  page or README says it is (`dem_figure.day_test_reason`: "within +/-0.15 m of the kept days,
  left out as taking them back would put a kept day over").
* **Even split.** With the days in two groups of about equal weight (say three at 0 m and three
  at +0.6 m on the same cells), the first day left out decides which group stays, whichever
  the noise puts further out (one pass rejected both). Nothing warns of this: read a day left
  out at little over the limit, or a camera's days left out together, with care.
* **More days left out than one pass, on long archives.** It leaves out fewer days than one pass
  on 7-day windows (real-week samples: 2.29 vs 2.65 per window at a day scatter of 0.08 m, 3.00
  vs 3.50 at 0.10 m) and in the planted cases above, but on a long archive whose day scatter is
  near the limit it can leave out MORE (26 weeks: 91 vs 76 over the limit in the first pass at
  0.10 m, 55 vs 48 at 0.08 m), because it trims the ends of the scatter as it goes [SYNTH,
  integration review round 1].
* **Trend WARNING.** A camera whose final day offsets run one way with time gets `WARNING:
  <cam>'s camera-day offsets run one way with time ...`: real beach change through the window,
  or an error that drifts with time (pointing, lens, a detection bias that follows the tides
  or waves); the day test then trims whichever end lies further out, and iterating does not
  cure that. The rule (`day_trends`): at least 5 testable days, max - min offset above the
  limit, |Spearman rho| of offset against date above 0.8 and, up to 8 days, a rho that days in
  random order reach in that direction at most 2 % of the time (exact over every order, for
  the offsets' own ties; `DAY_TREND_ALPHA` 0.02). That asks rho 1.00 over 5 days and 0.89 over
  6; from 7 days on it is |rho| > 0.8 itself (0.82 over 7, 0.81 over 8). The 2 % is one-sided:
  as the rule fires for either direction, its false-alarm chance per camera is up to twice that
  (1.7-3.4 % at the thresholds), and the WARNING gives both chances ('1 time in 720 (either way
  1 in 360)'). Its range is over all the camera's days, with the lowest and highest day, from
  the unrounded offsets. **The method alone can give it**: on a static synthetic beach with a
  fixed pointing (the Mar 2025 replica) c2's offsets give rho -0.96 and 0.165 m over 7 days, the
  same warning, so the WARNING says so, and only a drift clearly larger than that, or one the
  pointing check or a survey confirms, points at the camera or the beach.
* **The end-day rule.** One bad day at an end ranks first or last and lifts rho, so when exactly
  one of a camera's days was left out and it is its first or last, the camera's other days
  must still run the same way, |rho| > 0.8 in that direction (not the chance part, which the 4
  days left of 5 can never reach). When they would not pass the whole rule on their own the
  WARNING says the trend **rests on** that day ("it rests on <day>, the one day left out, at an
  end: ... the end of the trend or one bad day beyond it").

* **Survey-date anchor (`--anchor-date`, the survey-date builds only; Oct 2026).** On a camera
  whose days drift through the window the test keeps whichever end weighs more: c2 of the
  2024-10-23 build (first pass +0.33, +0.26, +0.11, +0.09, -0.13, -0.26 m over 20-26 Oct) kept
  24 and 26 Oct only and left out the survey's own day. `survey_products.py` now passes the
  survey date: each camera's testable days within 1 day of it are the reference -- every other
  day is measured first against their DEM and left out over the limit -- and the day(s)
  nearest the date are never left out; the iterative test then goes on with the rest. On the
  owner's real 2024-10-23 lines (after the Oct 2026 filter): c2 keeps 22, 23 and 24 Oct, leaves
  out 20 (+0.32 m), 21 (+0.29) and 26 Oct (-0.18 m from 22-24 Oct); 22 Oct stays +0.21 m off
  23-24 Oct, unjudged (only two c2 days left to judge it by), with the WARNING; 471 DEM cells
  (349 without the anchor). 2025-01-23, 2025-03-06 (no testable day within a day of the
  survey) and 2026-09-29 (no day over the limit): unchanged. Chosen by the owner, 15 Oct 2026,
  before DEM A went to Dr Sherwood; nothing of his survey is read (`test_anchor`: the days
  around the survey kept when the other days outweigh them; a date with no testable day near
  it: the test as without).

Where it shows: the log's `Day consistency` block (first-pass and final offsets, order left
out, not judged, taken back, refused) and its WARNING lines (`station_status.py` counts them,
one per warning); `<stem>_info.json` `filters.day_test` (`rejected`: each day left out with
`how` -- left out, not judged, refused take-back --, `offset_m`, `offset_at_step_m` and `why`;
`warning`; `trend_warnings` with `without.alone`); the DEM page's heading ("N camera-days off
the rest by > 0.15 m", or "left out by the day test (±0.15 m; K within it but at odds with a
kept day)", and the number of WARNINGs); a survey-date README (WHAT WAS USED, under the DEM:
each camera-day left out, how, its offsets and why; KNOWN CAVEATS: the WARNINGs, as
`provenance.json` caveats too); the comparison's FRAMES line (frames left out per reason).

| what | data | number | reproduce here | check on the NUC |
|---|---|---|---|---|
| One bad day; two that mask each other in one pass; a cluster that swamps good days (the whole-archive case, replayed); a day left out because a worse day shaped its DEM (taken back); a take-back that would push a kept day over (refused, recorded) | [SYNTH] test | only the planted days left out in each; the swamp replay keeps every post-storm day (one pass lost 4) | `python3 test_dem_day_test.py` (118 checks) | - |
| Lone partner, unjudged days, untested days, a day judged later, every day disagreeing | [SYNTH] test | a fringe day 0.57 m LOW whose cells one other day reaches is left out (the low beach -0.497 / -0.283 m, truth -0.50 / -0.30; kept, it read -1.07 / -0.86); two days that only see each other keep the first pass's verdict with the WARNING; four days that all disagree: every day rejected, exit 4, `dem_failure_reason` "rejected every camera-day" | same | the cron log |
| The station's 7-day offsets of 5-11 Oct 2026, replayed | [REAL offsets, SYNTH samples] | within 5.4 mm of the station's; 5 Oct c1 and 8 Oct c2 left out (as one pass), 11 Oct kept, the c1 trend WARNING; without 7 or 9 Oct (5 days, the low end trimmed): the WARNING, resting on 5 Oct | same (`test_nuc_7day`, `test_trend_end`) | `grep -n "Day consistency\|WARNING" logs/timex_cron.log` after a run |
| Trend WARNING rates, simulated full-pipeline windows (limit 0.15 m, fixed seeds, 1000 of each kind) | [SYNTH] | 7-day two-camera windows: no trend 2-6 %, one bad day (0.4 m) at an end 3-6 %, a trend of 0.03-0.06 m/day 67-87 %; a camera on 5 of 7 days: trend of 0.04-0.08 m/day 31-68 %, bad end day 4-6 %, no trend 1 %; on 6 days: trend 62-94 %, bad end day 6 %, no trend 1 % | the scratch simulation (numbers quoted in `dem_from_contours.py` above `DAY_TREND_MIN_DAYS`); the exact rule values in `test_dem_day_test.py` (`test_trend_rule`) | - |
| Survey-date fixture builds 2025-01-23 and 2025-03-06 with the first iterative version (daytest-fix d28b5ec); 2025-03-06 rebuilt at the integrated code (651dffb, floor venv) | [SYNTH] | DEM grids byte-identical to the one-pass builds, comparisons unchanged; at the integrated code 2025-03-06's grids identical again (Mar lidar median -0.034 m, NMAD 0.057 m, n 1263; waterlines -0.026 m, n 119643); Mar gets the trend WARNING (c2 +0.066 to -0.093 m over 3-9 Mar, range 0.165 m with 4 Mar's +0.072, rho -0.96, the fixture's known detector drift) | `survey_products.py --date ... --output-root <fixture root>` | the README's WHAT WAS USED and KNOWN CAVEATS |
| Run time | [SYNTH] | a synthetic year (730 camera-days, 668k samples, 167 left out) 68 s on Python 3.8 / numpy 1.17 (23 min for the first iterative version; 175 s CPU for the pre-PR one-pass code, 1-8 s for this PR's one pass); ~2.5 min of CPU at the station's density (1.5M samples). On the cron's whole-archive DEM the cost grows about as samples x days left out: seconds now; 52 weeks of real-week samples 25-163 s, 104 weeks with a seasonal +/-0.3 m swing 765 s (the old one pass ~370 s) [SYNTH, integration review round 1]. Before the archive is ~2 years old, cap the whole-archive DEM (e.g. its last N months) | `python3 test_dem_day_test.py --year` | the cron's run time |
| What goes to Chris says each day's reason | [SYNTH] test | a refused take-back reads "within +/-0.15 m of the kept days, left out as taking them back would put a kept day over", never "off the rest"; README lines per day (how, offsets, why); WARNINGs among the caveats; a trend 'rests on' its trimmed end day only when `without.alone` is false | `python3 test_survey_products.py` (`test_day_test_record`), `python3 test_survey_compare.py` (`test_fix2`) | a built README |

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
| 2025-03-19 (no photos, no ADCP) refused | config | exit 2: "Not built. To build it anyway: --force-disabled" | `python3 survey_products.py --date 2025-03-19` |
| 2024-10-23 is BLIND (DEM only): its survey is held by Dr Chris Sherwood for an independent check | [SYNTH] tests, files named like the Oct GCPs planted in every survey folder, the calibration folder, the output root and the current folder | never globbed, stat-ed or opened (a tracer on pathlib, os and open(), and an audit hook on every open and every command run); any look-up of a row of the date raises BlindSurvey (not an Exception: no `except Exception` swallows it); the compare step does not run; status `complete`; README and provenance say "blind: the survey is held by Dr Chris Sherwood (USGS) for an independent check; nothing in this build read it"; `--summary` lists it as `blind (DEM only)`; as a backstop, a runtime guard (an audit hook) refuses, in survey_products.py's own process, a Python-level open for reading of any `.csv`/`.txt`/... named like a survey of that date (by its name or the name a link points to; the in-process photo scan, run before any step is handed the photo folders, is stopped by a photo-named link to one), and no step is started with such a file as an argument; NOT covered: readers that open files in C (pandas' default CSV engine, OpenCV) and the steps' own processes, safe only because no look-up ever hands them such a path | `test_survey_products.py` (test_blind_and_calibration) |
| One calibration variant per root (`--eo ... --calibration-variant NAME`) | tests | a variant build outside a root named `<root>_NAME`, a table build in a root holding a variant, a missing `--eo` file (never the table's instead), an `--eo` file whose stated acceptance is FAIL (`NOT DELIVERED`, `# acceptance: FAIL`), or (on a blind date) states none: each refused, exit 5, nothing written; another calibration (by sha256) in the same root: refused without `--replace`; `--eo` without `--calibration-variant` or with `--all`: exit 6 | same |
| Another C in a root that holds a build | [SYNTH] tests + refreshes over copies of the real C = 0 and C = 0.037 2026 builds | refused, exit 5, nothing changed, no lock left; the message names the root that C belongs in (`<root>_c<C>`, the base root for the station's C); a root whose name gives it to the refused C although it holds another (a root without suffix holding a C = 0 build; a `..._c0` root holding C = 0.037) is said to hold that other C, and the run is given a root beside it (`<root>_c<C>`), never the same root; the same in the dry run, the compare-step warning, the CIRCULAR caveat and the end-of-run REPLACED line; `--replace` overwrites on purpose | `test_survey_products.py` |
| A `--setup-coef` of unknown origin | tests | exit 5, nothing written, unless `--setup-fitted-to` says where it was fitted | same |
| A usage error (unknown option; a bad `--steps`, `--window` or `--utc-hours`; a negative `--setup-coef`) | tests | exit 6, nothing written: never 2, which means a disabled date | same |
| Two builds of one date at once | tests | the second exits 4 and names the build holding the lock | same |
| `--all` exit code | tests | the worst by severity: failed (1) > refused (5) > busy (4) > partial (3) > 0 | same |
| Rerun on unchanged inputs | [SYNTH] + [REAL] | every step skipped, numbers identical (every date rerun several times on the final code) | rerun the same command |

### 3.5 2024-10-23: BLIND (DEM only), calibrations A and B

The 23 Oct 2024 GCPs are held by Dr Chris Sherwood (USGS) to check our DEM of that date blind.
Two DEMs, each in its own root: **A** with the table's calibration `CACO03_<cam>_20241023_EO`,
solved from those GCPs (his check of A tests the water level, the detection and the setup, not
the camera geometry: PARTLY-CIRCULAR for the geometry); **B** with a calibration made without
them, the 23 Jan 2025 GCP solution of the same mount turned back to the October pointing by
matching photos (c2 directly or through Dec 2024; or the Aug 2024 solution, if Chris sends Aug
photos), its camera position that solution's. B is computed on the NUC, LATER, by the
pre-registered GCP-free procedure `gcpfree_eo.py` (its docstring is the procedure; its gates and
thresholds were fixed on 13-14 Oct 2026, the horizon band and its two window parameters on 14 Oct,
before any run on real October photos; what its makers had seen of A: its EXPOSURE statement), under an outcome rule
the owner pre-registers and commits (`calibration/gcpfree/PREREGISTERED.txt`: both-or-none or
per-camera; the file also pins the photos a run selects, by name, size and sha256). Every run, a
trial too, is made under that committed file (there is no other), so the rule is fixed before the
code solves any photo, and every run is recorded in `calibration/gcpfree/runs/`. Expected before
any run: c1 NOT DELIVERED from the NUC's photos, c2 uncertain.
`survey_products.py` calls B's check INDEPENDENT only when the pinned procedure made it here and
delivered it (`GCPFREE_PINS`: its 'Made by' sha256, its modules and pre-registered constants; the
report beside it with the same B, verdict and file sha256 from a delivery run; true 'rotated by',
'pan change' and 'Position kept from' lines; its run folder this repository's
`calibration/gcpfree/runs/<id>`, under this repository's `PREREGISTERED.txt` unchanged since the
run, from the station's photo folders and the photos it pins, the first delivery run there with
every run's record intact and committed; nothing outside the branch's history (another branch, a
stash, the reflog) holding another pre-registration or run; the B file, report, run.json and
run_report.txt committed once), and only when the B DEM is the pre-registered B BUILD: the
committed tables, every result-changing option as A's build recorded it, the search envelope A's
detection used (path and sha256), A's forcing, one root (`B_BUILD.txt`); anything else is 'NOT the
pre-registered procedure' (or B build), PARTLY-CIRCULAR, and never built or shared under the
procedure's names. What no check can see: a B edited before its first commit with its report and
run.json rewritten to agree by an edit that keeps the 'rotated by' magnitude (a small tilt); a copy
of the procedure run outside this repository (another clone, a copy of the script) before the
pre-registration, or a re-clone after a first run; a branch or stash removed with its reflog; a
trial or cut-off run whose folder was deleted before its commit. Only the sha256s sent to Chris
beforehand (the procedure's before A, the pre-registration's, each B file's when its run ends)
answer those, with the time GitHub received the pre-registration commit: the owner commits AND
pushes `PREREGISTERED.txt` before any run of `gcpfree_eo.py`, a trial included (a commit's own
date is set by the computer that made it; the push time on GitHub is not). B's README (BLIND)
and the share MANIFEST say this plainly -- WHAT THE CODE ENFORCES, WHAT RESTS ON THE OWNER
FOLLOWING THE PROCEDURE, and the PRE-REGISTRATION COMMIT with whether a remote branch this
computer has fetched holds it (`surveys/README.md`, "What the code enforces, and what rests on
the owner").

| what | data | number | reproduce |
|---|---|---|---|
| The procedure recovers a known turn, and says NOT DELIVERED otherwise | [SYNTH] rendered direction-only worlds (lens and EO files of CACO03), the brief's test turns (c1 -12.5, c2 -13.4 deg pan) | c2 direct: 0.0025 deg from the truth, the Dec chain within 0.0005 deg of it, DELIVERED; c1 with Aug photos: Jan-based B primary, B_aug compared on the beach; c1 when Jan shares nothing with Oct but Aug does: B = the Aug solve (0.002 deg), the Aug position; Oct photos of another world: NOT DELIVERED, no B file | `test_gcpfree_eo.py` |
| The horizon check in a band | [SYNTH] the 'Jan solution' 3 deg of tilt off the world's truth | a perfect B: differential -0.084 / +0.727 deg = the measured expected bias -0.082 / +0.716 (the Jan photos turned by Q): PASS in the band, where a gate around 0 (+/-0.25 / 0.35) FAILED it; B off by 0.8 deg roll or 0.6 deg tilt: FAIL; fog: unavailable. (Found while testing: a turned frame's side has no reference content, whose fill was read as a horizon; fixed before any run) | `test_gcpfree_eo.py` |
| Pre-registration and run history | [SYNTH] temporary git repositories | no pre-registration, uncommitted, changed since its commit, committed but not pushed (its commit on no remote-tracking branch), its rule changed after a (trial) run: refused, exit 2 (a TRIAL too: one camera or other folders, before the pre-registration: refused, nothing solved, no run folder; the `--preregistered` option is gone); `--preregister` pins the 39 photos a run selects (name, size, sha256), and refuses with a station folder missing; a pinned photo removed or of other content (one byte added) after the pre-registration: a run refused before anything is solved, and the dry run says so; one-camera or other-folder runs: TRIAL, recorded in the same runs/, reports only, never a B file; a branch made from before the pre-registration (also with the other branch deleted: the reflog): `--preregister` there refused; the same records cherry-picked onto it (a rebase's copy): nothing foreign; run 1 per-camera: a c2 B and c1's NOT DELIVERED report; a second run without Aug photos and a third run: refused; with Aug photos: c1 solved, c2 kept, a SECOND product (gcpfree-aug); both-or-none with one camera passing: no B file (WITHHELD in its report); an unmounted photo folder, no Oct photo, an empty --aug-dir: refused, never a verdict, and `--dry-run` says beforehand that a delivery run would be refused; a link to A under the Jan EO's name: refused; a photo-named link to a forbidden file: never opened, not even for its size; run.json records the pre-registered constants in effect (their sha256 pinned with the code) and every repository module loaded, with its sha256 | `test_gcpfree_eo.py` |
| Never opened by the procedure | [SYNTH] decoys of calibration A, the calibration history and a GCP file in the test's calibration folder | none opened (an audit hook on every open made through Python, by any module, and a watch on `cv2.imread`, whose C++ reads no audit hook sees; both first shown to see a planted read of the decoy) | `test_gcpfree_eo.py` |
| The pre-registered B BUILD | [SYNTH] a scratch git repository with the tables committed, a fake A product (options, envelope and forcing recorded) | B built as A was: no problem; another `--dem-cell`, `--keep-moved-days` with another `--window`, another `--config`, no A product (`--a-root`), the envelope edited since A's detection: refused (`NOT the pre-registered B build`), the edited envelope also PARTLY-CIRCULAR; an envelope that is A's DEM (its GeoTIFF description: calibration A): PARTLY-CIRCULAR; B's forcing not A's byte for byte: the build stops; B's first build records its root and command (`B_BUILD.txt`), once; a second root: refused; `--share` before that record exists or is committed: refused; an uncommitted TRIAL in runs/: B not (yet) pinned, committed: pinned, its ORDER listing the trial with a NOTE; the first run reset away and a new run committed: NOT the pre-registered procedure (the reflog). A's DEM rebuilt: its first build time kept (`logs/first_built.json`), B's ORDER saying 'first built ... (last rebuilt ...)' and keeping the NOTEs; no A product: a NOTE | `test_survey_products.py` (test_gcpfree_b_build) |
| What each calibration rests on | [SYNTH] B files written by `gcpfree_eo.py`'s own writer, hand-made files in both formats | A: the 23 Oct 2024 GCPs, PARTLY-CIRCULAR for the geometry; the procedure's B traced to Jan 2025: INDEPENDENT; the same B with its tilt edited, another 'Made by' sha256, no report beside it, from a trial run, other modules, no 'rotated by' line, or withheld, and (its report and run.json made to agree with the edit) a 'rotated by' line naming another reference than its report's, a turn its numbers do not bear out (> 0.01 deg), its reference here of other content than the report recorded, or not here, and (the run written by `gcpfree_eo.py`'s own writer in a scratch git repository standing for this one) a copy out of the runs folder, a run under a second pre-registration (another folder or repository, the other rule), a run folder not committed, committed in a nested repository, edited after its commit (report and run.json rewritten: 'changed since its commit'; committed again: '2 commits'), photo folders not the station's (ROOTS patched in the running process: its constants not the pinned ones), a camera position 0.5 m off its reference, a false 'pan change' line, a second delivery run without Aug photos, the pre-registration changed and committed again after the run, the first run's record removed (from the working tree or in a commit) or rewritten as a trial: NOT the pre-registered procedure, never INDEPENDENT, not built or shared as `gcpfree`; NOT seen (said in the docs): a 0.05 deg tilt edit made before the first commit with the report and run.json rewritten; a hand-turned file (estimate_eo_rotation.py notes): traced, PARTLY-CIRCULAR; DELIVERED / DELIVERED-WIDER read as PASS / PASS-WIDER, NOT DELIVERED refused; a stated reference sha256, pan change or kept position the files do not bear out: NOT KNOWN | `test_survey_products.py` (test_gcpfree_b, test_gcpfree_where, test_gcpfree_turn_line, test_blind_and_calibration) |
| A one-camera B, and the procedure's names | [SYNTH] | `--eo c2` alone from a run that did not deliver c1: c1 left out, NOT DELIVERED with its report in the README, KNOWN CAVEATS, provenance and the GeoTIFF description; from a run that delivered c1: refused; B files of two runs, a run-1 B as gcpfree-aug, any other file as gcpfree: refused (exit 5) | `test_survey_products.py` (test_gcpfree_b) |
| A and B told apart once shared | [SYNTH] | B's DEM `2024-10-23_gcpfree_dem.tif` and maps `maps/2024-10-23_gcpfree_waterlines_plan.png`, `_c1`, `_c2` (A's `2024-10-23_dem.tif`, `maps/2024-10-23_waterlines_plan.png`), the plan view titled with the calibration identity; the GeoTIFF description read back: EPSG 6348 and 'calibration variant gcpfree (--eo): c2 CACO03_c2_gcpfree-oct2024_EO.yaml sha256 ... INDEPENDENT; c1 NOT DELIVERED ...'; `--share`: `2024-10-23_A_official_calibration.tgz` and `2024-10-23_B_gcpfree.tgz`, their MANIFESTs naming the calibration and the label, B's holding the pre-registration, the procedure line (run folder, pre-registration, commit), the stated uncertainty, EXPOSURE and ORDER; never written over; refused: a B no longer the procedure's (edited, its pre-registration changed after the run), and any product (A or B) whose calibration file is no longer the one it was built with (sha256 against the provenance and the detection's stamp) | `test_survey_products.py` (test_gcpfree_b) |
| B's stated uncertainty | [SYNTH] the estimator's sets (ultra/H) | rotation max(2 x spread, 4 px): the 4 px floor is above every delivered synthetic B's true error (1.8-3.6 px); those errors were almost all roll (tilt 0.008-0.02 deg), so their millimetre DEM errors say nothing of the term as a tilt, and no synthetic DEM figure is given for it: B's README and MANIFEST state the DEM effect of the rotation term taken as a pure tilt, per range band, from the date's own DEM and calibration (4 px is 0.06-0.09 deg of tilt, which moves waterlines metres at 200-400 m from the camera); inherited position: in the largest of five synthetic position cases, c1 median 0.011 m, p90 0.026 m; c2 0.018 m, 0.086 m (not a bound; the cases were written after A's rows had been seen, and c1's comes from a -0.55 m height case equal to c1's A minus Jan offset: said in EXPOSURE). Neither is a validated bound on real photos. Where the date's DEM gives no tilt table (too few cells with a local slope in a camera's view, even from one-sided differences), the README and MANIFEST say why, not 'no DEM yet' | `gcpfree_eo.py` constants; B's README (`test_survey_products.py` test_gcpfree_b) |
| What Chris is told the checks cannot prove | [SYNTH] test | B's README BLIND section and its MANIFEST: what the code enforces; what rests on the owner (the pre-registration committed AND pushed to GitHub before any run, a trial included; no copy run elsewhere; no history deleted; no B edited before its commit; the clock); the pre-registration commit, said to be on no remote branch when none holds it; A's MANIFEST: none of it. A NOTE when A's DEM was rebuilt, or first built, after B's run began (B's options are A's as last built; said 'first built' when it was) | `test_survey_products.py` (test_gcpfree_b, test_gcpfree_b_build) |
| B's camera position | `calibration/calibration_history.csv` [REAL] | the Aug 2024 and Jan 2025 GCP solutions of the same mount differ by 0.36 m (c1: 0.34 m horizontally, 0.11 m in height) and 1.23 m (c2: 1.08 m, 0.59 m); B keeps the Jan (or Aug) position; B's build uses nothing of the Oct solution | B's README, KNOWN CAVEATS |
| B's build never opens A's calibration file | [SYNTH] tests (audit hook) | none opened | `test_survey_products.py` |
| What Chris is told of B's making | [SYNTH] tests | the procedure's EXPOSURE statement (A's rows seen on 13 Oct 2026 before the design; the synthetic turns are A minus Jan; the position-inheritance cases written five minutes after A's rows were seen, three of them matching c1's A minus Jan offset; A's azimuths read on 14 Oct, the day the horizon band and its two window parameters were set) in every B file, report and run.json, and in B's README BLIND section and MANIFEST; the ORDER, with dates: the procedure's pin (its first commit in `survey_products.py`), the pre-registration (the commit that set the rule in force, every commit of the file with its rule), A's DEM (first built, kept across rebuilds, and last rebuilt; from A's root), the B run and B's commit, every run of the procedure in runs/ with its start; a NOTE when the rule was changed, when A's DEM came before the rule in force was committed or before B's commit, when a run started before the rule, a trial before B's run, a run never finished, no A product was found, or B was built with other code than A; each camera's search envelope with its sha256 | `test_survey_products.py` (test_gcpfree_b), `test_gcpfree_eo.py` (test_runs) |
| A and B built end to end on the synthetic October fixture, at oct-blind's code (before the integration; the rows below it give the integrated code) | [SYNTH] the October fixture (photos of 20-26 Oct, the Chatham transfer, files named like the Oct GCPs planted in the survey, calibration and current folders); clean clones of the branch standing for the station's repository; B made there as on the station (deploy notes B1-B3): `gcpfree_eo.py --dry-run`, `--preregister` (rule per-camera, chosen by the test) and its commit, the delivery run with no option, its run folder committed in one commit, the station's three photo folders bound in a private mount namespace to rendered direction-only photos (plumbing: their turns are the brief's test sizes, not the fixture's pointing); every process under strace -f and an audit hook | gcpfree: the dry run without the folders says a delivery run would be REFUSED; both cameras DELIVERED (exit 0); a second set (c1's October photos share nothing with Jan): c1 NOT DELIVERED, c2 DELIVERED (exit 3); a second pre-registration (both-or-none) made after that run in another folder and committed there: its run is a TRIAL (exit 4), reports only, no B file. A: `status complete`, compare not run, PARTLY-CIRCULAR per camera, `dem/2024-10-23_dem.tif` (808 cells) EPSG:6348 with the calibration identity, plan view titled with it. B, both cameras: complete, INDEPENDENT, `dem/2024-10-23_gcpfree_dem.tif` (609 cells), maps `2024-10-23_gcpfree_waterlines_*`, provenance recording B's commit; B, c2 only (`--a-root`): complete, INDEPENDENT, c1 NOT DELIVERED in the README, provenance, GeoTIFF description and MANIFEST (478 cells); their READMEs carry the procedure line (run folder, pre-registration, commit), the rotation term's DEM effect as a tilt per range band (c2, 0.062 deg: waterlines move 0.5-2.6 m, DEM 0.04-0.12 m at 0-250 m), EXPOSURE and ORDER with dates; the DEM page heading fits. Refused (exit 5, nothing written): the c2-only B before its run folder was committed (saying so), A into B's root, B into A's root, B as `gcpfree2`, a delivered camera left out, B files of two repositories' runs, a B copied out of `runs/`, a B edited after the run with its report and run.json rewritten (build and `--share`: 'not the one this product was built with', 'changed since its commit'). `--share`: `2024-10-23_A_official_calibration.tgz` (15 files), `2024-10-23_B_gcpfree.tgz` (22 files; 19 for c2 only), never written over. Reruns at the final code: every step skipped. No process of any build, dry run, share or gcpfree run opened, stat-ed or read a file named like the held survey (strace: 0, after one lstat of a GeoTIFF-description argument string, removed in 6506a38); calibration A opened only by A's own runs and by one refused B dry run (naming the table's calibration of the camera it left out); B builds only stat it. The 2025-01-23 and 2025-03-06 fixture builds rebuilt from the forcing up at the final code (a clone at another path: the step commands name the script's path): DEM, spread and count grids identical cell by cell to the builds of 10 Oct (ultra/C/final), every comparison number and label identical (Jan lidar -0.043 m, NMAD 0.062, n 526; Mar lidar -0.034 m, NMAD 0.057, n 1263); GeoTIFFs EPSG:6348 (was 32619) with the calibration identity; waterline rows' water_level_source `adcp` (was `tide_model`, the Oct 2024 path fix 0f97a27); re-detection rounding in pixel_row (<= 0.007 px), peak_sharpness (<= 0.04) and easting/northing (1-2 mm) on a few % of rows; beach elevations identical | a scratch fixture, not in the repository |
| A and B rebuilt end to end on the synthetic October fixture at the code of the integration's review round 1 (`68253ea`, floor venv) | [SYNTH] the October fixture's photos only (a survey folder holding nothing but the Jan lidar envelope; no file named like the held survey in any folder given); a clone of the branch for the station's repository with a bare repository for GitHub; B made as on the station (dry run; `--preregister` per-camera, committed; a `--camera c2` run then REFUSED, exit 2, 'its commit ... is on no remote-tracking branch (not pushed)'; pushed; dry run 'pushed: on origin/...'; the delivery run, 174 s, both DELIVERED; its folder committed and pushed), the station's photo folders bound in a private mount namespace to rendered direction-only photos (plumbing, as above); an audit recorder in every Python process | A: exit 0, `status complete` in 880 s (at 651dffb/123224b it was PARTIAL, exit 3: c2's two planted pan days 20-21 Oct outnumbered the reference day, and the README advised leaving out the calibration's own day); no other pointing named (2 moved days against the reference day and four zero-turn days); c2 20-21 Oct DIFFERENT POINTING and 22 Oct UNCONFIRMED (its 263 features at 0.00 deg are no reading next to 21 Oct's pan: a day lost, not a wrong one), all left out; `dem/2024-10-23_dem.tif` 660 cells (54 frames: c1 36, c2 18; 808 cells and c2 30 frames at oct-blind's code, which kept 22 Oct); the README's unchecked days say their features matched within their own noise, the converted-waves setup error, the ADCP datum's size and the 0.6-0.8 x DEM share. B (both cameras): exit 0, complete in 892 s, INDEPENDENT, 460 cells (35 frames: c1 17, c2 18), its reference 'the day nearest the survey (this calibration was turned to the window's own photos, pooled: it has no solve day)' with no 'no photo of the calibration day' caveat, its ORDER NOTEs (A before the rule; B committed after A's DEM), the PRE-REGISTRATION COMMIT 'held by origin/...', provenance `dem.files` naming `2024-10-23_gcpfree_*`; the rerun after `B_BUILD.txt`'s commit skips every step. `--share`: A 38 files (12.0 MB), B 48 files (10.8 MB), each with forcing/, pointing/, the waterline file the DEM used and logs; B's with the Jan 2025 GCP solutions it was turned from at the sha256 its reports recorded and its MANIFEST label 'INDEPENDENT (procedure made with stated exposure to A ...)'. `--summary` over each root: exit 0. The audit recorder: 0 of ~50,000 paths opened or listed by any build, run, share or summary named like the held survey; calibration A opened by A's builds only (none by B's build, the procedure or B's share) | a scratch fixture, not in the repository |
| The real A and B DEMs | [NUC] | A: built on the NUC first; B: later, after the owner's pre-registration of the outcome rule (and, ideally, Aug 2024 photos from Chris); to be filled in from their READMEs | the commands in `surveys/README.md` |

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

* **Used for.** 2024-10-23 (BLIND, DEM only: no comparison here) and hours the ADCP lacks. The 2025-01-23 and
  2025-03-06 windows are inside the ADCP record: water level measured.
* **Value** [REAL], fitted on the ADCP overlap: harmonic transfer, 13 constituents, M2 ratio
  1.481, Marconi leads by 70 min, residual x 0.85; leave-one-week-out RMS 0.076 m,
  extrapolation RMS 0.085 m pooled (`sigma_m` inside the record's span). A period wholly
  before the record gets the backward test's error: RMS 0.095 m (daytime 0.085), bias
  +0.037 m (fit 2025-01-01..03-10, predict 2024-12-09..12-31); after it, the forward test's
  0.077 m. The mean error of a whole left-out week (13.5-18 UTC) ranged -0.089..+0.102 m,
  RMS 0.051 m over 13 weeks: a week-long DEM carries one such mean in every frame
  (`forcing.json` `water_level.week_mean_error`). The linear alternative
  (Marconi = 1.403 x Chatham(t + 72 min) - 0.033 m, i.e. Marconi leads by 72 min; the code's
  `lag` is -72 min in `Marconi(t) = a x Chatham(t - lag) + b`) scored 0.166 / 0.174 m and is
  not used.
* **Live station** (`marconi_water_level.py`, for the OWG products, not the waterlines):
  refitted to GNSS-R over the last 30 days on every run; default (1.24, -48 min, -0.10 m,
  sigma 0.13 m) only without GPS overlap.
* **Limits.** The transfer changes with the inlet and harbour; months from the overlap the
  error can exceed the extrapolation RMS, which the two tests measure on predictions 0-23 days
  before and 0-29 days after their fit windows. 20-26 Oct 2024 lies 44-51 days before the ADCP
  record, beyond the backward test's 23 days: further than tested (the forcing report, its
  `transfer_sigma_measured: false` in forcing.json and the product README say so).
  [SCRATCH] Within the backward test the far half (15-23 days back) scored RMS 0.140 m, bias
  +0.109 m, but that is the calendar stretch (9-24 Dec 2024) that every backward fit predicts
  worst (fit from 23 Jan: 30-45 days back RMS 0.114, bias +0.068; fit from 1 Feb: 45-60 days
  back RMS 0.136, bias +0.102, while its 15-45 days back stay at 0.057-0.063): an event of
  those weeks, not an error measured to grow with distance. How it grows beyond 23 days back
  is not known.
* **Re-check.** `forcing_report.txt`, block `Chatham -> Marconi transfer`.

### 4.5 Wave height and period conversions (dates or hours without the ADCP)

* WIS ST63064 Hs -> ADCP Hs: linear form `0.833 x Hs + 0.004 m` (r 0.67, CV RMS 0.456 m) or
  the direction-aware form (r 0.84, CV RMS 0.330 m), which wins; WIS Tp scaled by 1.061.
  NDBC 44013: `0.827 x Hs + 0.151 m` or direction-aware (r 0.87). Combined setup CV RMS
  0.079 m [REAL], 534 h of December overlap for WIS. Not used by the two enabled 2025 dates
  (ADCP waves throughout). Source: `historical_forcing.py`, `forcing_report.txt`.
* The conversions are chosen by the setup error they cause. Until this PR's review round a
  `--setup-coef 0` build of hours without the ADCP (no setup, no error to choose by) took the
  first conversion listed (Hs linear, Tp raw: WIS Tp bias -0.35 s), another wave currency than
  the C = 0.037 build's. That mattered as soon as 2024-10-23 was switched on (`chatham` era,
  WIS/NDBC waves throughout): its README recommends a C = 0 build with `--fit-setup` to refit C,
  whose lines would then have carried other waves. Fixed: with C = 0, `historical_forcing.py`
  judges the conversions with the station's C (`conversions_judged_with_c` in forcing.json), so a
  C = 0 build converts every hour exactly as the C build does [SYNTH] (`test_historical_forcing.py`,
  "C = 0 currency": the same Hs, Tp and sigmas on all 168 converted hours). It costs no extra
  rebuild: `historical_forcing.py` changes in this PR anyway, so every date's forcing and
  detection rerun once after the pull.

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
    trials, 29 of them against another day); no matched trial under 0.45 deg was called MOVED (0 of 229), every one of 0.55 deg or
    more was (157 of 157). Below 40 the errors grow (20-39 inliers: up to 0.18 deg; 10-19: 0.58;
    under 10: anything). A turn of more than 0.1 deg is the scene (dirt matches at zero) and
    measures the day's change from the reference directly: a change of tilt or roll over 0.3 deg
    (the step's tilt/roll tolerance, as for the horizon) or of pan over 0.5 deg (pointing_check's)
    is a DIFFERENT POINTING whatever the horizon reads; a smaller turn is a 'small change', kept (a
    review probe: 1500 inliers at +0.15 deg against a clear-looking +2.0 deg jump). A match reads
    tilt and roll to ~0.1 deg (MATCH_ERR_DEG; the same planted rotations against another day: <= 0.101
    deg at 40-59 inliers, 29 trials; <= 0.094 at 60-99, 41; <= 0.060 at 100-199, 11, over every run
    of that study -- its first summary, without one run, read 28 / 35 / 10; this last bound rests on
    11 trials only), and a clear horizon scatters about as much from day to day at rest (the
    NUC's clear days of 29 Sep - 2 Oct, their features within 0.02 deg: c1 +0.90 to +0.96 deg of tilt,
    +0.16 to +0.26 of roll): a CLEAR horizon that reads a change over 0.3 deg and agrees with the
    features within 0.1 deg in tilt and in roll makes the day DIFFERENT POINTING too, at any inlier
    count (review 3: 0.28 deg read on 50 inliers, the clear
    horizon at +0.37, was a 'small change', kept; the deployed rule left it out). Review 3's fix cut
    this at 100 inliers and confirmed within 0.3 deg; review 4 found it not monotonic: on 50 inliers a
    0.05 deg tilt was confirmed by a clear +0.34 and not by +0.40, and at 150 inliers 0.29 deg was never
    confirmed by a clear +0.36. A horizon further than 0.1 deg from the features is not the change
    they measured, and the features decide: a jump of 0.34 or 0.40 against 0.05 is a 'small change'
    alike, and the band where a horizon confirms is the features' own error, whatever the inlier
    count. What it does not do: a real change of 0.3-0.4 deg read 0.1-0.2 deg apart by the two (each
    within its own scatter) is judged by the features alone (DIFFERENT POINTING if they read over 0.3).
    The cost is at the tolerance itself: a true 0.30 deg tilt (review 4's T03_L2, round 3's T03) read 0.30
    by the features and 0.300-0.303 by the clear horizon is now DIFFERENT POINTING on the days where
    the horizon reads over 0.300 (round 3's T03: c1 1/2 Oct, c2 2 Oct; review 4's T03_L2: c2 2 Oct),
    a 'small change' on the others: a day lost at a change equal to the tolerance. A turn of
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
    0.298 or 0.300 deg. A turned match's tilt/roll is now judged on its own (with a clear horizon's
    confirmation within 0.1 deg, above). An earlier
    proposal to require 60 inliers was not taken: it was chosen to make c1 27 Sep (48 inliers)
    weak, but that day's -1.21 deg of pan is real: tracked dune landmarks and a more sensitive
    match read -1.20 to -1.25 deg on 26-28 Sep, and the camera was back by 29 Sep.
  * **What the ~0.1 deg rests on** (integration review round 1). The planted rotations against
    another day are of CACO05 c2 2026 only; the matches of a frame against itself add CACO05 c1 2026
    and CACO03 c1 Jan 2025; no CACO04 view (2025-03-06) and no CACO03 c2 view was tried. Every match
    against another day had its inliers spread 180-292 px along the minor axis of their scatter
    (292-618 px along the major one). A match whose inliers are bunched reads roll, and tilt, less
    well: a review model (pointing_check.solve_rotation, the real CACO05 c2 lens, 50 rays, 0.015 deg
    noise, kept at >= 40 inliers) read roll to <= 0.022 deg with the inliers over the whole frame,
    but up to 0.17 deg (p90 0.09) in a 300 px patch and 0.12 deg in a 500 px one [SYNTH, SCRATCH].
    So the step now writes each match's spread (`inlier_spread_px` in pointing_<cam>.csv,
    `pointing_check.inlier_spread`), and a deciding match spread under 180 px says so in its why and
    in the README (the days kept on such a match); the verdict is not changed, as there is no
    real-photo evidence to set a threshold by. The whys read 'a match of 40 or more features spread
    over the frame is good to ~0.1 deg on the views tested'.
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
    horizon that jumps: UNCONFIRMED; no horizon fitted: 'unchecked', and UNCONFIRMED next to a day
    left out on firm evidence, below). The reference day's own frames are exempt (two frames of one
    day match at zero by nature), and a row of an older build without `turn_sigma` counts as a reading (the NUC's
    2026-09-29 table below). What it does not catch: a dirt consensus of thousands of inliers may
    read up to ~5 by the images' own systematic errors (the copies did, at 3600-4800 inliers) and is
    then still taken for the scene; under a clear horizon that does not jump, dirt can still hide a
    pan (no horizon sees one); a day nothing could read at all (its horizon unfitted, its match
    blind) is 'unchecked' and kept unless a day left out on firm evidence is next to it (below;
    review 2's DIRT: c1 26-28 Sep, turned -1.2 deg, the dirt also broke all three horizon fits:
    kept, as by every earlier rule; review 2's DIRT2_ST10, below: c2 tilted +1.0 deg, the scene and
    the drops both under 40 inliers, the fits failed: kept, as by the deployed rule); a real match
    under 4.0 is not a reading, so its day is judged without it (a lost day at worst, not a wrong
    one). Only the NUC's next build of 2026-09-29 gives the real values for its
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
    day nothing could read ('unchecked': no feature reading, and no horizon fitted, or one fitted
    that is not a clear horizon) next to a FIRM day is UNCONFIRMED too, as are runs of such days from
    it: nothing shows the camera back (review 2's simR, the real storm frames as 17:00 with perturbed
    copies as 16:30/17:30: c1 28 Sep's 3-frame fit failed at 14.6 px, and the day was 'unchecked',
    kept and gridded while c1 was turned -1.2 deg; on the NUC its fit held, 3.5 px, and the day was
    UNCONFIRMED anyway). Firm: a day measured to differ; UNCONFIRMED on more than a horizon reading
    that is not clear (a clear-looking jump, an ambiguous match); or UNCONFIRMED on such a misread
    horizon alone in a STORM. A misread horizon alone never shows a day moved, so on a calm day one
    never reaches past itself: in review 3's fog week (one misread fog edge with unmatched features on
    19 Jan between matching days, the fits failing on 20-22 Jan) an earlier form of this rule left out
    19-22 Jan, 4 of 7 days, where the same days next to a matching day stayed kept; 20-22 Jan are
    'unchecked', kept, as by the deployed rule. In a storm a camera can turn (c1 turned -1.2 deg over
    the three storm days of 26-28 Sep 2026) while foam or spray hides both the horizon and the land,
    and the storm's days come misread or unread in any order. A storm is a run of consecutive days,
    each left out on its own readings or read by nothing, that holds a firm day of the first two
    kinds, misread horizons on two days or more, or a storm day: inside a known storm (STORMS) or with
    offshore Hs over 1.5 m (STORM_HS_M, the DEM's own limit) in the forcing's waves, which the forcing
    step writes before the pointing step (`storm_days_of`; without them, the STORMS days only; the
    window's storm days are part of the pointing step's signature, so a new wave record reruns the
    verdicts when it changes them). Review 3's fix asked for two misread days side by side; review 4
    found two shapes it kept, turned (exact ray rotations of the six real storm frames, `--steps
    pointing` on the station-floor venv): ST10_S3 (c1 tilted +1.0 deg on 26-28 Sep: 26 Sep misread,
    the fits of 27/28 Sep failed at 33.0 / 10.7 px with 6 / 2 inliers) and S2728_T10_FOG (+1.0 deg on
    27-28 Sep with the land fogged: c1 misread, unread, misread; c2 misread, unread, unread): both are
    now left out (26 Sep is inside the 25-26 Sep storm; c1 26 and 28 Sep are two misread days in one
    run), and the NUC's own c1 rows with 27 Sep made unread (fit failed, no features) leave out 26-28
    Sep (pinned in the test, with ST10-like rows: 'unchecked' on calm days, UNCONFIRMED with STORMS or
    the waves over 1.5 m on a day of the run). What a storm costs on calm-looking rows: a misread day
    with Hs over 1.5 m (2026-09-29 itself read 1.52 m at its peak) passes UNCONFIRMED to the days
    nothing read next to it, lost days, and the DEM's wave filter leaves many of their frames out
    anyway. A 'horizon only' day (a clear horizon at the reference's pointing, no feature reading)
    shows tilt and roll back but never the pan: next to a day that measured a pan over 0.5 deg (or an
    ambiguous match whose runner-up put one forward), it is UNCONFIRMED too (review 4, P10_30_FOG01:
    a 1.0 deg pan on 30 Sep and 1 Oct, 1 Oct fogged over the land; 30 Sep measured, +1.00 / +1.02 deg
    on 1851 / 4239 inliers, 1 Oct 'horizon only', kept on both cameras; P10_L2_FOG02 the same on 1-2
    Oct; by review 4's estimate a 1 deg pan moves ground points ~3.8 m at 150-400 m, ~0.1-0.4 m of
    DEM at beach slopes of 0.05-0.1). Next to a measured tilt or roll with no pan, a clear horizon at the reference's shows
    the camera back, and the day stays 'horizon only' (pinned). Without a reference day between, a
    window whose features match on one day only and show a pan takes in every 'horizon only' day
    connected to it: nothing measured the pan anywhere else. Never the reference day itself: the
    others are measured against its frame, so a neighbour that differs from it says nothing of it,
    and runs stop there. The row's `next_to` column names the day a run starts at.
    UNCONFIRMED days are left out like DIFFERENT POINTING ones (`--keep-moved-days` keeps both),
    whatever share of the window they are. **'Most days differ from the reference' no longer keeps
    days.** That guard (older than this rule) kept EVERY day when more than half of the measured
    days differed from the reference, the UNCONFIRMED ones and days that differ from each other
    included (review 3: H_T06, real frames with an exact +0.6 deg tilt of c1 on 30 Sep - 2 Oct, read
    +0.592-0.598 deg on 1674-1686 inliers with a clear horizon agreeing: 3 of 4 measured days
    differ, all 7 c1 days kept, 26-28 Sep with them, by the previous form of this rule and by the
    deployed one; and a
    NUC-like window with c1's -1.21 deg pan on 27 Sep and an exact +0.45 deg roll on 30 Sep and
    1 Oct: 3 of 5 differ, all kept, the rolled days gridded, ~0.6 m of DEM at 250-350 m). Now every
    such day is left out, and when the days measured to differ agree with each other (0.3 deg of
    tilt and roll, 0.5 deg of pan) and outnumber the days at the reference's pointing, the step and
    the README name that other pointing: the window holds two, and the product rests on the
    reference's; if the reference day is the odd one out, that is the wrong one, which only the
    photos tell (`--keep-moved-days` keeps both). The warning and the build's missing part give the
    way to build on the other pointing instead: a `--window` over the other pointing's days on one
    side of the reference day (the side with more of them), in another output root, makes the day
    nearest the survey in it the reference (review 4, REF_T06: 29 Sep, the survey day and the
    reference, tilted 0.6 deg; every other day measured -0.60 deg and left out, the product resting
    on the one moved day: `--window 2026-09-30 2026-10-02`). Refusing to grid such a window until the
    owner chooses is not done (an owner's ruling). Since the integration's review round 1
    (POINTING_VERSION 11) two things more: (1) the days at the reference's pointing are the days
    measured there AND the days that read nothing else whose 40 or more features matched the
    reference at zero turn within their own noise (`zero_turn_days`: no reading of a turn, but what
    the reference's pointing gives, or dirt): on the synthetic Oct 2024 fixture c2's two planted
    +1.0 deg pan days (20-21 Oct) outnumbered the reference day alone, while 22 and 24-26 Oct matched
    it at 0.00 deg on 222-263 features, and the build was PARTIAL, its README saying 'the
    reference's is the smaller group'; (2) when the reference is the calibration's own solve day
    (`reference_is_calibration_day`: the date in the calibration's file name, which the step takes
    as the reference whenever it has photos), its pointing is the calibration's: the warning says
    the product rests on the calibration's pointing, gives no rebuild without that day (for A on
    2024-10-23 that would put the 23-Oct calibration on the other pointing's days, a pan error no
    horizon sees), and it is a partial reason only while the detection has not run with these
    results or keeps both pointings (`--keep-moved-days`). A day that read nothing but a blind
    zero-turn match next to a measured pan (22 Oct there) stays UNCONFIRMED and left out: a day
    lost, not a wrong one [SYNTH] (`test_survey_products.py` test_oct_pointing_pattern). When the days left out are more than half of the
    window's, the step and the README say so with the count, how many of them rest on a horizon
    reading alone that is not a clear horizon, how many are left out only for the day next to them,
    and the options there are: `--keep-moved-days` keeps every such day, and there is no per-day option
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
    results give), naming what that detection left out. A stamp of an older build (the deployed
    code's) does not record `--keep-moved-days`; it is never read as having run with it: when its
    left-out days are not those the current results leave out by default, the README says it ran
    with other pointing results or with `--keep-moved-days` and that its record does not say which
    (review 3: such a stamp that left out no day, after `--steps pointing` with this rule flagged
    c2 26 Sep, made the README say 'KEPT in the detection (--keep-moved-days)' while the step said
    'LEFT OUT of the detection'). `--keep-moved-days` leaves out no day, so it is offered as the
    reason only for a stamp that left out none; one that left days out ran with other pointing
    results, and the README says that alone. Days these results keep against the horizon are 'kept'
    only once the detection has run with them: after a pointing-only rerun over a deployed build that
    left c2 26/27 Sep out on the old rule, the README said 'kept, although the sea-horizon check alone
    would have left them out' for the days the detection on disk had left out, and the STORM caveat
    called them '(matches) left out by the pointing check' (review 4); it now says these results keep
    them and the detection on disk left them out, kept once it runs again, and the STORM caveat says
    'left out by an earlier pointing result'. The documented full command reruns the detection and
    its README matched a fresh build's (review 4). The console of a pointing-only run (`--steps
    pointing`) said 'LEFT OUT of the detection' as if done; it now says 'to be LEFT OUT of the
    detection ... when the detection runs again with these results', with what the detection on
    disk left out (or that none has run), and 'LEFT OUT' only when the run includes the detection
    (`test_pointing_verdicts`).
  * **The NUC's 2026-09-29 table under the rule** (`test_survey_products.py`,
    test_pointing_verdicts, on the NUC's rows with the edge numbers of the real frames): c1 26 and
    28 Sep UNCONFIRMED (no matched features, a misread horizon), 27 Sep DIFFERENT POINTING (the
    features' -1.21 deg of pan on 48 inliers; the horizon's +4.82 deg named a misread), 29 Sep -
    2 Oct matches; c2 every day matches (26/27 Sep kept, the horizon overruled; 28 Sep: its fit
    failed, 393 features decide). Left out: c1 26-28 Sep, as before, with the measured reason; c2
    none (was 26-27 Sep). The NUC's rows predate `turn_sigma`; with the values the real storm frames
    project for c2 26-28 Sep (~17 / ~9 / ~9, the window-dirt bullet above) the verdicts are the same.
    If 26 and 27 Sep read under 4, both are UNCONFIRMED and left out, as the deployed rule left them
    out; 28 Sep (its fit failed) follows them only if its own match also reads under 4: at its
    projected ~9 its 393-inlier match is a reading and the day matches (3 / 3 / 8.4: 26-27 Sep left
    out; 2.0 on all three: 26-28 Sep left out, one day more than the deployed rule, a day lost, not a
    wrong one; all three cases pinned in the test). A value under 4 loses a feature reading; it never
    keeps a day the deployed rule left out.
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
    measured (the storm days are not), two of them moved, so the 'most days differ from the
    reference' guard kept them, with its warning; c2 left them out. Since review 3 c1 leaves them
    out too, and names the other pointing where the moved days agree (E's two tilted days: pinned in
    the test). Scenarios A-D (moves of
    0.3-1.0 deg on clear and storm days): every clear day moved by more than 0.3 deg of tilt or roll
    or 0.5 deg of pan is DIFFERENT POINTING, and the 0.3 deg ones are 'small change'; storm days
    whose features matched are measured (c2 27/28 Sep in A); storm days without matched features
    are UNCONFIRMED, except a pure 1 deg pan with the horizon unfitted (A, c1 28 Sep: 'unchecked',
    nothing there can see a pan). The 'most days differ' guard turned on where the scenario moved
    most measured days (A c2, D c1) and kept them; since review 3 they are left out.
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
    (T03): a small change on both days on both cameras (deployed: the same; since review 4, whose clear
    horizon confirms within 0.1 deg at any inlier count, DIFFERENT POINTING where the horizon read just
    over 0.3: c1 1/2 Oct and c2 2 Oct, recomputed from the rows); 0.6 deg (T06):
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
    the 7 measured days, the 'most days differ' guard kept them, with its warning (with real
    day-to-day matches, as in the pinned test, they are contradicted and UNCONFIRMED). Since review
    3 the guard keeps nothing: those 4 days, which agree with each other, are left out as the
    window's other pointing and the build is partial (4 days lost to a front that the stand-ins'
    blind matches could not contradict).
  * **Review 3 scenarios, end to end with this rule** (the review's own photo sets: the six real
    storm frames as 26-28 Sep, the real 7 Oct 16:00 frame with the NUC's day-to-day turns as 29 Sep
    - 2 Oct, exact ray rotations; `--steps pointing`, the station-floor venv; the 20 sets run once,
    their verdicts recomputed from the rows with the final rule, and H_T06, REF_T06, S_T06, S_T10 and
    the base set run again with it) [REAL frames, SYNTH changes]. Against the review's run of the
    previous form of this rule: the base set and the seven L_* sets (0.3-1.0 deg of tilt, roll or pan
    on 1-2 Oct) give the same verdicts and the same days left out. H_T06 (c1 and c2 tilted exactly
    +0.6 deg on 30 Sep - 2 Oct): c1 read +0.59-0.60 deg on 1674-1686 inliers, 3 of its 4 measured
    days, so the 'most days differ' guard had kept all 7 c1 days, the UNCONFIRMED storm days with
    them; now 6 of 7 are left out, the other pointing is named (+0.60 deg of tilt, median) and the
    build is partial; c2 (3 differ, 3 at the reference's) as before. REF_T06 (the reference day 29 Sep
    tilted 0.6 deg, every other day at the calibration's pointing, the case the guard was for):
    both cameras name the other pointing (c1 3 days, c2 5) and leave out every day but the reference
    day, so the product rests on the one day that moved, and the build is partial, saying so; the
    guard had kept every day (6 right, 1 wrong). S_* (c1 and c2 turned 0.3-1.0 deg on all three
    storm days): c1 28 Sep, its fit failed (8.3-50.5 px) and its features unmatched, stays UNCONFIRMED
    next to 26/27 Sep, misread two days running, in 8 of 9; S_T10, where 27 Sep's fit failed too,
    kept 27/28 Sep under review 3's rule (its photo set is gone; review 4's ST10_S3, the same shape,
    is now left out, above). A NUC-like window (turn_full: c1 as on the NUC, and an exact +0.45 deg
    roll on 30 Sep and 1 Oct read on 2103-2286 inliers): 3 of 5 measured days differ, not in one way
    (a -1.21 deg pan, a +0.45 deg roll): no other pointing; the guard had kept all 7, now 5 are left
    out (pinned in the test, with H_T06). End to end (the review runner's full build of that window,
    its synthetic storm waterlines): c1's waterlines hold 29 Sep and 2 Oct only, where the DEM had
    gridded 52 c1 frames of all 7 days, the rolled ones among them. A detection stamp written as the deployed code writes it
    (no `--keep-moved-days` recorded, no day left out) under results that leave out c1 26-28 and
    c2 26 Sep: the README says the detection left out no day where these results leave out those,
    by other pointing results or `--keep-moved-days`, never 'KEPT in the detection
    (--keep-moved-days)'. A turned match under 100 inliers that reads 0.28 deg with a clear horizon
    at +0.37 is DIFFERENT POINTING (it was a 'small change'; since review 4 at any inlier count, the
    two agreeing within 0.1 deg).
  * **Review 4 scenarios** (the review's 40 photo sets: the six real storm frames as 26-28 Sep, the
    real 7 Oct 16:00 frame with the NUC's day-to-day turns as 29 Sep - 2 Oct, exact ray rotations of
    0.3-1.0 deg of tilt, roll or pan on storm and clear days, fogged land; `--steps pointing`, the
    station-floor venv, no forcing built, so the storm days are the STORMS ones, 25-26 Sep) [REAL
    frames, SYNTH changes]. Their measured rows were judged again with this rule, and the eight sets
    it changes or that test it (ST10_S3, S2728_T10_FOG, P10_30_FOG01, P10_L2_FOG02, T03_L2, REF_T06,
    T06_30_FOG01, base) were run again end to end: the same verdicts. Turned days kept on either
    camera, over the 33 sets whose turns the review's summary lists: 27 by review 3's rule, 17 by
    this one; the 10 now left out are ST10_S3 c1 27/28 Sep,
    S2728_T10_FOG c1 27 and c2 27/28 Sep (storm runs), P10_30_FOG01 1 Oct and P10_L2_FOG02 2 Oct on
    both cameras ('horizon only' next to a measured pan), and T03_L2 c2 2 Oct (a 0.30 deg tilt, the
    tolerance itself, confirmed by a clear horizon reading 0.303). The 17 still kept: STAY_T06_27
    (10: the camera moved on 27 Sep and stayed, the reference with it: one pointing, the
    reference's), REF_T06's reference day (2: the product rests on it, the build partial, with the
    `--window` that builds on the other pointing), and 0.30 deg turns at the tolerance (ST03_S3 c2
    27/28 Sep, T03_L2 c1 1/2 Oct and c2 1 Oct). No unturned day outside the storm is newly left out.
    A pointing-only rerun over a deployed build (the review's synthB_dep, which left c2 26/27 Sep out
    on the deployed rule) now writes a README that says what the detection on disk did and keeps
    nothing it left out.
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
    matters when its horizon also jumps on an edge that is not clear (then UNCONFIRMED). Rebuilt
    once more with the rule of review 3 (POINTING_VERSION 9; the same roots and commands; the Jan
    and Mar detection up to date, so their pointing step and what follows it reran, 26 s and 16 s;
    the live-era roots rebuild every step): the same left-out days (none on the fixture dates and
    the live roots; the storm build c1 26-28 Sep and c2 26 Sep) and every headline line identical.
    The Mar README now says why the pan was not measured on its 'horizon only' days (their features
    matched within their own noise), and the stand-in and storm builds' WHAT THIS IS line adds the
    DEM's own count (85 of the 90 waterline frames, on 29 Sep - 2 Oct). Rebuilt again with the rule
    of review 4 (POINTING_VERSION 10; the same roots and commands, the Jan and Mar detection up to
    date): the same left-out days and every headline line identical; the pointing CSVs add
    `next_to`, the pointing JSON the window's storm days (Jan 19-20, Mar 6-7: offshore Hs over 1.5 m;
    the 2026 roots 26 Sep, inside STORMS, and 29 Sep, 1.52 m), which change no verdict there (no
    misread day with unread neighbours, no measured pan next to a 'horizon only' day). The review
    runner's full builds of the storm window (real, synth and synthB photos, its synthetic storm
    waterlines) give the same left-out days and headlines; one why adds the measured pan it names
    ('a pan of -1.22 deg, which no horizon sees').
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
  elevation across a 1:10 foreshore. GeoTIFFs are tagged EPSG:6348 by default since Oct 2026
  (the 2025 lidar files and the 2026 RTK are in that frame); before, EPSG:32619 (WGS 84 / UTM 19N),
  a nominal tag a GIS places ~1-1.5 m off; a product keeps its tag until its DEM step reruns.

### 4.9 Thresholds that move the numbers (not offsets, but each changes a result)

| setting | value | where | why |
|---|---|---|---|
| Max offshore Hs for the DEM | 1.5 m | cron `DEM_MAX_HS` | the setup errs most in big waves |
| Day-offset test | 0.15 m, iterative (a day left out only when judged; take-back; section 2.3a); not run on fewer than 3 camera-days | cron `DEM_MAX_DAY_OFFSET` | catches a camera-day displaced as a whole |
| Day-trend WARNING | >= 5 testable days, range > the limit, \|rho\| > 0.8 and, up to 8 days, a chance <= 2 % (1.00 over 5 days, 0.89 over 6); the lone-end-day rule | `dem_from_contours.py` `DAY_TREND_*` | section 2.3a |
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

* **The 23 Oct 2024 DEMs.** Checked blind by Dr Chris Sherwood against the GCPs he holds: A
  tests the water level, detection and setup (its geometry is fitted to his points), B the
  whole chain, geometry included. Nothing on our side may read his file; A and B are compared
  with each other (DEM difference, `dem_change.py`) only after both are built.

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
  writes the real `turn_sigma` of its storm days (projected ~17 / ~9 / ~9 for c2 26 / 27 / 28 Sep; 27 and
  28 Sep have the least margin; real frames here read 28 Sep closest to the limit, 5.3-6.5). Read
  them in `pointing_c2.csv` before relying on c2 26-28 Sep. Below 4 the error is not a harmful one:
  c2 26/27 Sep would be UNCONFIRMED and left out, as the deployed rule left them out, and 28 Sep,
  whose fit failed, with them only if its own match also reads under 4 (days lost, not wrong ones).
  Not caught: a dirt consensus of thousands of inliers reading above the threshold, and dirt hiding
  a pan under a clear horizon that does not jump.
* **Storms the rows cannot see.** A storm is told from a calm fog spell by the wave record (offshore
  Hs over 1.5 m on a day of the run, or a STORMS entry) or by two misread days in one run of days
  left out or unread. A storm whose every day nothing could read (review 2's DIRT: c1 26-28 Sep,
  turned -1.2 deg, the drops breaking every horizon fit and match) is still kept: no day of it is left
  out on its own readings. Leaving such days out on the wave record alone, or keeping calm misread
  days' neighbours whatever the waves, is the owner's call; the 1.5 m threshold is the DEM's own,
  not measured as a storm signal for camera turns (the one storm seen: 1.95 / 3.71 / 2.11 m at its
  peak on 26 / 27 / 28 Sep 2026).
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
  (That is PR #87's history.) THIS PR: three branches, each built and adversarially verified on
  its own before they were integrated: pointing-verdict-fix (four review rounds of the pointing
  rule, reviews 1-4 in section 4.6), oct-blind (two review rounds of the blind 23 Oct 2024 A/B
  path and the GCP-free procedure) and daytest-fix (the iterative day test; its last commit,
  `846a2be` -- the lone-end-day rule and the exact tie chances -- had NO independent review of its
  own). The integration commits that write the day test's record into what Chris reads, the
  frame counts and the B protocol texts (`6096fe1` .. `651dffb`) were then reviewed in one
  integration review round (runner on a simulated NUC, real-data and fixture builds, blind,
  science, station and merge lenses), whose findings were fixed in the commits after `123224b`
  (the `--summary` crash, the Oct fixture's partial A, the pushed pre-registration and its
  re-pin, the share archive, the error budget, the trend warning, the inlier spread). Those
  fixes were checked by their own tests, the full test list on Python 3.11 and the station floor
  and an Oct fixture rebuild, not yet by a further independent review.

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
